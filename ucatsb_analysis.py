"""Shared analysis logic for UCATS-B flight records: masking, cal-event
detection, calibration and its uncertainty, and the two export writers.

Imported by `ucatsb_gui.py`, which is the only user interface. This module is
deliberately **Qt-free** -- it needs pandas and (for the calibration figure)
matplotlib, nothing more -- so a batch script or notebook can reuse
`calibrate_series`, `calibration_uncertainty` or `export_icartt` without
pulling in a GUI toolkit.

It was once `plot_co2_timeseries.py`, named for a standalone CO2-only figure
CLI that lived at the bottom. That CLI was removed when the GUI superseded it:
it hardcoded CO2 and one set of masking settings, and it read the cal pairing
from `cals.yaml`'s `cals:` block, which describes the tanks plumbed in *now*
and so mislabels the tanks on any earlier flight -- with no way to correct it.
Everything here is gas-agnostic.
"""
from pathlib import Path

import pandas as pd
import matplotlib.dates as mdates

LINE_COLOR = "#2a78d6"
RIGHT_AXIS_COLOR = "#8e44ad"   # purple, distinct from the red/orange masking shades
CAL_SHADE_COLOR = "#898781"
PRESSURE_EXCLUDE_COLOR = "#d03b3b"
WARMUP_EXCLUDE_COLOR = "#ffa64d"   # light orange
PUMPS_EXCLUDE_COLOR = "#8e7cc3"    # violet, distinct from the other bands
POST_CAL_FLUSH_COLOR = "#2fa88a"   # teal, distinct from the warm-up/pressure shades
# Manually flagged points. Near-black and drawn as an x, so a struck-out
# reading reads as struck out rather than as one more colored band: it is the
# only removal the user made by hand, and the only one not explained by a
# setting somewhere in the panel.
FLAGGED_COLOR = "#111111"
CAL0_COLOR = "#eda100"   # golden
CAL1_COLOR = "#0d366b"   # dark blue
GRID_COLOR = "#e1e0d9"
AXIS_COLOR = "#c3c2b7"
TEXT_COLOR = "#0b0b0b"
MUTED_COLOR = "#52514e"

# The detector's spec pressure. The *tolerance* around it is not here: it is
# a per-gas setting the user tunes, and lives in the GUI's
# DEFAULT_GAS_SETTINGS. So are the warm-up length and the cal-mean windows,
# which had fixed values here only to serve the removed CO2 CLI.
D1_P_TARGET_MBARS = 140.0
CAL_MERGE_GAP_S = 2   # bridge cal periods split by a single dropped-flag sample
CALS_YAML_PATH = Path(__file__).parent / "cals.yaml"


def drop_presync_rows(df: pd.DataFrame, jump_threshold_s: float = 60.0) -> pd.DataFrame:
    """Drop leading rows recorded before the logger's clock synced.

    On startup the datalogger's system clock is sometimes stale until GPS/NTP
    sync completes, producing a burst of rows with bogus (too-late)
    timestamps followed by a large backward jump to the true time. Any rows
    before the last such jump have unreliable timestamps and are dropped.
    """
    diffs = df["datetime"].diff()
    big_jumps = diffs[diffs < -pd.Timedelta(seconds=jump_threshold_s)]
    if big_jumps.empty:
        return df
    return df.loc[big_jumps.index[-1]:].reset_index(drop=True)


def find_intervals(datetimes, mask):
    """Return a list of (start, end) timestamps for contiguous True runs in mask."""
    intervals = []
    in_span = False
    span_start = None
    for t, flag in zip(datetimes, mask):
        if flag and not in_span:
            in_span = True
            span_start = t
        elif not flag and in_span:
            in_span = False
            intervals.append((span_start, t))
    if in_span:
        intervals.append((span_start, datetimes.iloc[-1]))
    return intervals


def merge_close_intervals(intervals, gap):
    """Bridge intervals separated by a gap smaller than `gap` into one."""
    merged = []
    for start, end in intervals:
        if merged and (start - merged[-1][1]) <= gap:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


# Ozone readings below this are instrument faults, not measurements (the Feb
# 2025 flight has 7 of them, down to -2292 ppb). The floor is deliberately
# well below zero rather than at zero: a real near-zero ozone measurement
# scatters negative, and 168 readings in -15..0 ppb on that flight are the
# sensor's noise about a small true value -- throwing those away would bias
# the low end upward and hide how noisy the instrument actually is.
O3_VALID_MIN_PPB = -15.0

# Same idea for the water vapour instrument. Precautionary rather than
# demonstrated: the Jul 2026 flight has no w_H2Obest reading below -5 ppm at
# all (its minimum is +14), so this removes nothing today and exists to catch
# the fault mode when it appears.
H2O_VALID_MIN_PPM = -5.0


def below_floor_mask(values, floor):
    """Flag readings below a physical floor as instrument faults.

    Only ever a *display/pairing* filter: like every other mask in this
    module it marks rows, and the raw column is never modified. NaN is not
    flagged -- absent is not the same as invalid, and the two get different
    treatment downstream (a gap versus a removal).
    """
    if floor is None:
        return pd.Series(False, index=values.index)
    return values.notna() & (values < floor)


# --------------------------------------------------------------------------
# Manually flagged rows, as run-length-encoded row ranges
#
# The user can strike out points no automatic rule catches (an ozone spike
# well above the record, a stretch of nonsense after a valve glitch) and those
# removals persist in the flight's config. Ranges rather than a list of row
# numbers: one dragged gesture is one entry whatever its width, so a config
# stays readable at thousands of flagged points.
#
# Ranges are INCLUSIVE at both ends, and every function here returns them
# sorted, non-overlapping and adjacency-merged. That canonical form is what
# makes the on-disk representation stable, which matters because _is_dirty()
# deep-compares it: two different gesture orders reaching the same set of rows
# must produce the same YAML, or the app would offer to save nothing.
#
# Set arithmetic is the whole reason ranges beat storing the drawn rectangles:
# unflagging is a difference, with no ordering or paint/erase semantics to get
# wrong, and the flagged set cannot shift underneath the user when an unrelated
# setting changes.
# --------------------------------------------------------------------------

def merge_ranges(ranges):
    """Canonical form: sorted, non-overlapping, adjacent runs joined.

    Adjacent (`hi + 1 == lo`) merges as well as overlapping, so flagging rows
    10-20 and then 21-30 stores one range and not two -- otherwise the same
    set of rows would have more than one representation and the dirty-state
    comparison would report a change where there is none.
    """
    cleaned = []
    for lo, hi in ranges:
        lo, hi = int(lo), int(hi)
        if hi < lo:
            lo, hi = hi, lo
        if hi < 0:
            continue
        cleaned.append((max(lo, 0), hi))
    if not cleaned:
        return []
    cleaned.sort()
    merged = [cleaned[0]]
    for lo, hi in cleaned[1:]:
        last_lo, last_hi = merged[-1]
        if lo <= last_hi + 1:
            merged[-1] = (last_lo, max(last_hi, hi))
        else:
            merged.append((lo, hi))
    return merged


def add_ranges(ranges, lo, hi):
    """Union: flag rows lo..hi in addition to whatever is already flagged."""
    return merge_ranges(list(ranges) + [(lo, hi)])


def subtract_ranges(ranges, lo, hi):
    """Difference: unflag rows lo..hi. A cut inside an existing range splits
    it in two, which is what makes unflagging exact rather than approximate --
    it never has to discard the whole gesture that produced the range."""
    lo, hi = int(lo), int(hi)
    if hi < lo:
        lo, hi = hi, lo
    out = []
    for r_lo, r_hi in merge_ranges(ranges):
        if r_hi < lo or r_lo > hi:      # untouched
            out.append((r_lo, r_hi))
            continue
        if r_lo < lo:                   # keep the head
            out.append((r_lo, lo - 1))
        if r_hi > hi:                   # keep the tail
            out.append((hi + 1, r_hi))
    return merge_ranges(out)


def ranges_to_mask(ranges, index, offset=0):
    """Boolean Series over `index`, True on the flagged rows.

    `offset` is subtracted from the stored row numbers: ranges are recorded
    against the RAW file's row numbering (so they mean the same thing as the
    companion CSV's rows, and survive a change in how many pre-sync rows get
    trimmed), while the analysis frame starts at raw row `presync_dropped`.
    Rows outside the frame are simply not marked rather than raising -- a
    range can legitimately reach into the trimmed head.
    """
    mask = pd.Series(False, index=index)
    n = len(mask)
    for lo, hi in merge_ranges(ranges):
        lo, hi = lo - offset, hi - offset
        if hi < 0 or lo >= n:
            continue
        mask.iloc[max(lo, 0):min(hi, n - 1) + 1] = True
    return mask


def ranges_row_count(ranges):
    """How many rows the ranges cover, for the readouts."""
    return sum(hi - lo + 1 for lo, hi in merge_ranges(ranges))


def cal_switch_mask(datetimes, cal_mask):
    """Flag the switch-over sample at the end of each cal period: the first
    row whose solenoid flag has gone False while the cell still holds cal gas.

    `j_sol_cals`/`j_sol_aircal` describe the *valve*, but the reading carrying
    a given timestamp is of gas that entered the cell a second or so earlier.
    So the first ambient-flagged row after an injection is not air at all --
    it is the tank, at full concentration (206.6 ppm CO2 against 206.51
    assigned, on the Feb 2025 flight), and its detector pressure is back in
    spec by then, so neither the cal mask nor the pressure mask removes it.
    One row per cal event, 18 on that flight, every one of them a tank-composition
    point sitting in the middle of the ambient record.

    Nothing else catches it. `find_intervals` closes a run *at* this row's
    timestamp, and `post_cal_flush_mask` starts strictly after that, so the row
    falls in the gap between the two masks -- and the flush is 0 by default
    anyway, while this row is wrong regardless of what the flush is set to.

    Rows sharing that timestamp are flagged together: the file logs the same
    second more than once, and half a switch-over is not a useful thing to
    flag. A dropped sample mid-injection (flag False for one row inside a cal
    period) is flagged by the same rule, which is also right -- it is cal gas
    too.

    This must be applied to `cal_mask` (what counts as *not air*), never to
    the cal *intervals*: those set `Cal_p`, which every cal-mean window is
    measured from, so extending them would silently move every cal mean.
    """
    cal_mask = cal_mask.astype(bool)
    previous = cal_mask.shift(1, fill_value=False)
    first_false = (~cal_mask) & previous
    if not first_false.any():
        return first_false
    # Group adjacent rows that share a timestamp, then take whole groups: the
    # switch-over may be logged as two rows with identical times.
    groups = (~datetimes.eq(datetimes.shift(1))).cumsum()
    return groups.isin(groups[first_false].unique()) & ~cal_mask


def post_cal_flush_mask(datetimes, cal_intervals, flush_s, cal_mask=None):
    """Flag ambient rows within `flush_s` seconds after each cal period ends.

    The detector cells still hold cal gas when the solenoid switches back to
    ambient, so the air data immediately after an injection reads toward the
    tank rather than the atmosphere. The width is instrument behaviour, not
    something derivable from the file, so it is a user control.

    Rows inside a cal period are never flagged, even when one injection
    follows another closely enough for the flush window to reach into it --
    those rows are already cal, and double-flagging them would shade the same
    span twice and confuse the "air we threw away" note.

    Returns an all-False Series (not a scalar) when flush_s is 0, so callers
    can use it unconditionally.
    """
    flagged = pd.Series(False, index=datetimes.index)
    if not flush_s or not cal_intervals:
        return flagged
    width = pd.Timedelta(seconds=flush_s)
    for _, end in cal_intervals:
        flagged |= (datetimes > end) & (datetimes <= end + width)
    if cal_mask is not None:
        flagged &= ~cal_mask
    return flagged


def box_stats(datetimes, values, t0, t1, y0=None, y1=None):
    """Summarise the points inside a drawn box (Igor-style marquee stats).

    A plain selection tool: whatever is inside the box is what gets counted.
    It applies no masking of its own -- drawing the box is the user's way of
    saying which data they mean, and the vertical bounds already exclude the
    cal dives when the box is drawn around the ambient band.

    `y0`/`y1` are optional: pass neither to select on time alone (used when
    the chosen trace lives on a different Axes from the one the box was drawn
    in, where the box's y-range means nothing).

    Returns `values_stats` plus the box bounds and `n_clipped` -- points
    inside the time range but outside the box vertically. A 2D marquee
    truncates the distribution, which narrows the reported sigma, so the count
    is what makes an accidental clip visible instead of silent.
    """
    in_x = (datetimes >= t0) & (datetimes <= t1)
    finite = values.notna()
    if y0 is None or y1 is None:
        in_y = pd.Series(True, index=values.index)
        lo = hi = None
    else:
        lo, hi = min(y0, y1), max(y0, y1)
        in_y = values.between(lo, hi)

    keep = in_x & in_y & finite
    stats = values_stats(values, keep)
    stats.update({
        "keep": keep,
        "t0": t0, "t1": t1, "y0": lo, "y1": hi,
        "n_in_span": int((in_x & finite).sum()),
        "n_clipped": int((in_x & finite & ~in_y).sum()),
    })
    return stats


def values_stats(values, keep):
    """n / mean / std / min / max of `values` over the rows `keep`
    selects, skipping missing values. Split out of box_stats so a second
    series can be summarised over an existing selection."""
    kept = values[keep].dropna()
    n = len(kept)
    return {
        "n": n,
        "mean": float(kept.mean()) if n else None,
        # ddof=1 needs two points; a one-point box gets None rather than NaN.
        "std": float(kept.std(ddof=1)) if n > 1 else None,
        "vmin": float(kept.min()) if n else None,
        "vmax": float(kept.max()) if n else None,
    }


def shade_intervals(ax, datetimes, mask, color, alpha):
    """Shade contiguous stretches of x where mask is True."""
    for start, end in find_intervals(datetimes, mask):
        ax.axvspan(start, end, color=color, alpha=alpha, linewidth=0)


def bottle_for_interval(df, start, end):
    """Identify which digital cal state (0 or 1, from j_sol_cals) was active
    during [start, end]. This is the raw solenoid state, not a bottle
    identity -- which physical bottle it corresponds to is determined
    separately, by matching measured concentration (see match_cal_serial),
    since the two need not agree with any naming in cals.yaml."""
    sub = df[(df["datetime"] >= start) & (df["datetime"] <= end)]
    return int(round(sub["j_sol_cals"].fillna(0).mean()))


def _read_cals_yaml(path: Path):
    """Parse cals.yaml into a raw dict, or {} if unavailable.

    cals.yaml uses tabs as "key:\tvalue" separators in places, which the YAML
    spec disallows outright; tabs never carry structural indentation here, so
    blanking them out is safe. Keep this workaround even when the file is
    clean -- a hand-edit or re-copy from the acquisition repo can reintroduce
    them at any time.
    """
    import yaml
    if not path.exists():
        return {}
    try:
        return yaml.safe_load(path.read_text().replace("\t", " ")) or {}
    except yaml.YAMLError as e:
        print(f"Warning: could not parse {path}: {e}")
        return {}


def load_cal_assignment(path: Path):
    """The `cals: {cal0: ..., cal1: ...}` block -- which two roster tanks are
    plumbed in -- as {"cal0": serial, "cal1": serial}, or {} if unavailable.

    cals.yaml describes the CURRENT run, so this is only ever the *default*
    for a flight. An older flight flew different tanks; the GUI stores the
    per-flight choice in its <dataset>_conf.yaml and overrides this.
    """
    data = _read_cals_yaml(path)
    cals = data.get("cals")
    return {k: v for k, v in cals.items() if isinstance(v, str)} if isinstance(cals, dict) else {}


def select_cal_bottles(roster, serials):
    """Narrow a full tank roster down to the ones actually plumbed in, as
    {serial: {"CO2": ..., "N2O": ..., ...}}.

    The single place that decides what bottle *matching* is allowed to see,
    and it must stay narrow: an unrelated roster tank could otherwise sit
    closer to a measured value than the real one and produce a wrong
    identification. Unknown serials are dropped rather than raising, and a
    duplicate (the same tank picked for both states) collapses to one entry
    -- calibrate_series degrades to offset-only there rather than failing.
    """
    return {serial: roster[serial] for serial in dict.fromkeys(serials) if serial in roster}


def load_cal_roster(path: Path):
    """Load EVERY tank in cals.yaml's roster, not just the two named in the
    `cals:` block.

    This exists solely for the mismatch advisory (see cal_mismatch_notes) --
    it must never feed bottle matching, which goes through select_cal_bottles
    for the reason documented there.
    """
    data = _read_cals_yaml(path)
    return {k: v for k, v in data.items() if k != "cals" and isinstance(v, dict)}


def match_cal_serial(value, gas_key, cal_bottles):
    """Return whichever bottle serial's nominal concentration for `gas_key`
    is closest to the measured `value`, or None if no match is possible.
    Matching by measured concentration (rather than trusting any assumed
    digital-state-to-serial mapping) is self-correcting if bottles are
    swapped between flights.
    """
    candidates = [
        (serial, nominal[gas_key])
        for serial, nominal in cal_bottles.items()
        if gas_key in nominal
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda c: abs(c[1] - value))[0]


def most_common_serial(cal_points, digital_state):
    """Return the most frequently matched bottle serial among cal_points for
    the given digital_state (0 or 1), or None if none matched."""
    from collections import Counter
    serials = [s for _, _, state, s in cal_points if state == digital_state and s]
    if not serials:
        return None
    return Counter(serials).most_common(1)[0][0]


def mean_std_label(values):
    """Format a "mean ± std" string for a set of cal-point values."""
    import statistics
    if not values:
        return "n/a"
    mean = statistics.mean(values)
    if len(values) > 1:
        return f"{mean:.2f} ± {statistics.stdev(values):.2f}"
    return f"{mean:.2f}"


def cal_mean_points(df, cal_intervals, value_col, cal0_window, cal1_window,
                     cal_bottles=None, gas_key=None, exclude_mask=None):
    """For each calibration interval, average `value_col` over a window
    relative to the interval's last timestamp (cal_p), tag which digital
    state (0 or 1) was active (from j_sol_cals over the whole interval), and
    identify the bottle serial by matching the measured mean against
    cal_bottles' nominal concentrations for gas_key (if provided).

    cal0_window / cal1_window: (start_offset_s, end_offset_s), e.g. (-15, -1),
    applied as [cal_p + start_offset, cal_p + end_offset]. The digital state
    active during the interval selects which window is used.

    exclude_mask: an optional boolean Series aligned with df (True = row is
    masked out, e.g. warm-up or out-of-spec pressure) applied to the raw data
    before averaging. A window left with no unmasked rows is dropped
    entirely -- a cal point can disappear rather than average over bad data.

    Returns a list of (x_time, mean_value, digital_state, serial) where
    digital_state is 0 or 1, serial is a matched bottle serial or None, and
    x_time is the last timestamp within the (masked) window.
    """
    points = []
    for start, end in cal_intervals:
        digital_state = bottle_for_interval(df, start, end)
        offsets = cal0_window if digital_state == 0 else cal1_window
        window_start = end + pd.Timedelta(seconds=offsets[0])
        window_end = end + pd.Timedelta(seconds=offsets[1])
        window = df[(df["datetime"] >= window_start) & (df["datetime"] <= window_end)]
        if exclude_mask is not None:
            window = window[~exclude_mask.loc[window.index]]
        if window.empty:
            continue
        value_mean = window[value_col].mean()
        if pd.isna(value_mean):
            continue
        x_time = window["datetime"].iloc[-1]
        serial = (
            match_cal_serial(value_mean, gas_key, cal_bottles)
            if cal_bottles and gas_key else None
        )
        points.append((x_time, value_mean, digital_state, serial))
    return points


# --------------------------------------------------------------------------
# Calibration: turning the per-injection cal means into a correction for the
# ambient data.
#
# Drift removal and calibration are ONE step here, not two. Each bottle's
# measured response is interpolated in time, and the two-point line is solved
# at every ambient timestamp:
#
#     slope(t)     = (A_hi - A_lo) / (R_hi(t) - R_lo(t))
#     intercept(t) = A_lo - slope(t) * R_lo(t)
#     calibrated   = slope(t) * measured(t) + intercept(t)
#
# where R_b(t) is bottle b's interpolated response and A_b its assigned mole
# fraction from cals.yaml. Slow drift falls out automatically because the cal
# responses themselves carry it -- there is no separate detrending pass. A
# two-point (slope + intercept) form is required rather than a simple offset:
# the measured span error is several percent, so gain must be corrected too.
# --------------------------------------------------------------------------

CAL_DRIFT_MODELS = ("linear", "smooth", "constant")
CAL_DEFAULT_SMOOTH_EVENTS = 3
CAL_GAP_FACTOR = 3.0        # node gap > this x median spacing => extrapolated
CAL_SLOPE_SANE_RANGE = (0.5, 2.0)


def cal_bottle_series(cal_points, cal_bottles, gas_key, drop_mismatched=True):
    """Group cal_points into one response series per bottle, keyed by digital
    state.

    Each state's bottle identity is the *consensus* serial across its points
    (most_common_serial), not any individual point's match. With
    drop_mismatched, a point whose own matched serial disagrees with that
    consensus is rejected: its averaging window straddled a solenoid
    transition, so it actually measured the *other* tank. Letting such a point
    become a drift-model node corrupts the calibration on both sides of it --
    a real occurrence (a 374.66 ppm point tagged to the ~206 ppm bottle).

    Returns {state: {"serial", "assigned", "assigned_unc", "times", "values",
    "rejected"}}, with times/values sorted by time.
    """
    series = {}
    for state in sorted({s for _, _, s, _ in cal_points}):
        consensus = most_common_serial(cal_points, state)
        nominal = (cal_bottles or {}).get(consensus) or {}
        kept, rejected = [], []
        for t, value, point_state, serial in cal_points:
            if point_state != state:
                continue
            if drop_mismatched and consensus and serial and serial != consensus:
                rejected.append((t, value, serial))
            else:
                kept.append((t, value))
        kept.sort(key=lambda p: p[0])
        series[state] = {
            "serial": consensus,
            "assigned": nominal.get(gas_key),
            "assigned_unc": nominal.get(f"{gas_key}_unc"),
            "times": [t for t, _ in kept],
            "values": [v for _, v in kept],
            "rejected": rejected,
        }
    return series


def drift_nodes(times, values, model="linear", window=CAL_DEFAULT_SMOOTH_EVENTS):
    """Reduce a bottle's per-injection means to the nodes the calibration
    interpolates between.

      "linear"   -- the means themselves; interpolation happens at evaluation
                    time. Assumption-free, but passes each event's noise
                    straight through into the calibrated ambient data.
      "smooth"   -- centred rolling mean over `window` events (min_periods=1,
                    so the first/last events keep full coverage rather than
                    dropping out). Suppresses event-to-event scatter, which
                    can be several times the within-event noise, while
                    preserving slow drift.
      "constant" -- the flight mean, emitted at the first and last node times
                    so all three models flow through one evaluation path.
    """
    times, values = list(times), list(values)
    if not times:
        return [], []
    if model == "constant":
        mean = sum(values) / len(values)
        return [times[0], times[-1]], [mean, mean]
    if model == "smooth":
        smoothed = (
            pd.Series(values)
            .rolling(window=max(1, int(window)), center=True, min_periods=1)
            .mean()
        )
        return times, smoothed.tolist()
    return times, values


def interp_hold(node_times, node_values, target_times):
    """Interpolate a node series onto target_times, piecewise-linear in time,
    holding the first/last node value flat outside the node range.

    Returns a Series on target_times' *index*, not its values: the datetime
    column contains duplicate timestamps (over a thousand in a typical
    flight), and while reindexing *by* a duplicated label set is legal,
    building a Series *from* one is not.
    """
    if len(node_times) == 0:
        return pd.Series(float("nan"), index=target_times.index)
    nodes = pd.Series(list(node_values), index=pd.DatetimeIndex(node_times)).sort_index()
    nodes = nodes[~nodes.index.duplicated(keep="last")]
    targets = pd.DatetimeIndex(target_times)
    union = nodes.index.union(targets.unique())
    filled = (
        nodes.reindex(union)
        .interpolate(method="index", limit_area="inside")
        .ffill()    # flat hold past the last node
        .bfill()    # flat hold before the first node
    )
    return pd.Series(filled.reindex(targets).to_numpy(), index=target_times.index)


def loo_residuals(times, values):
    """Leave-one-out residual per node: how far each cal mean sits from what
    that bottle's *other* nodes predict at its time.

    This is the honest scatter metric, and the reason it exists is subtle:
    residuals about the fitted drift model are useless for the default
    "linear" model, which interpolates through every node and so has
    identically zero residual at each. Leave-one-out instead answers "how
    wrong is the calibration *between* cal events", and is insensitive to
    genuine slow drift (which both neighbours share).
    """
    out = []
    for i, (t, v) in enumerate(zip(times, values)):
        others_t = times[:i] + times[i + 1:]
        others_v = values[:i] + values[i + 1:]
        if not others_t:
            out.append(float("nan"))
            continue
        pred = interp_hold(others_t, others_v, pd.Series([t]))
        out.append(v - pred.iloc[0])
    return out


def _gap_spans(times, factor=CAL_GAP_FACTOR):
    """Spans between consecutive nodes that are anomalously long (> factor x
    the median spacing) -- a telemetry dropout or a stretch where cals stopped.
    The calibration is interpolated across them, but shouldn't be trusted
    there, so they get flagged alongside the true extrapolation regions."""
    if len(times) < 3:
        return []
    deltas = [(b - a).total_seconds() for a, b in zip(times, times[1:])]
    median = sorted(deltas)[len(deltas) // 2]
    if median <= 0:
        return []
    return [
        (times[i], times[i + 1])
        for i, d in enumerate(deltas) if d > factor * median
    ]


def cal_mismatch_notes(bottles, gas_key, roster, rel_threshold=0.01):
    """Warn when a bottle's measured response sits far from its assigned
    value, naming the nearest roster tank as a prompt to check cals.yaml.

    Deliberately advisory only, and deliberately never auto-substituting:
    cals.yaml describes the *current* run, so applying it to an older flight
    can silently use the wrong tank -- but a large offset is equally
    consistent with a real instrument gain error. The nearest-roster-tank hint
    is a starting point for a human, not a conclusion; read it together with
    span_gain, which is what actually distinguishes the two cases.
    """
    notes = []
    for state in sorted(bottles):
        info = bottles[state]
        assigned, values = info.get("assigned"), info.get("values")
        if not assigned or not values:
            continue
        measured = sum(values) / len(values)
        rel = (measured - assigned) / assigned
        if abs(rel) <= rel_threshold:
            continue
        note = (f"{info['serial']}: measured {measured:.3f} vs assigned "
                f"{assigned:.3f} ({rel:+.2%}).")
        nearest = match_cal_serial(measured, gas_key, roster or {})
        if nearest and nearest != info["serial"]:
            note += f" Nearest roster tank: {nearest} ({roster[nearest][gas_key]:.3f})."
        note += " Check that cals.yaml matches the tanks actually flown."
        notes.append(note)
    return notes


def calibrate_series(df, value_col, cal_points, cal_bottles, gas_key,
                     model="linear", smooth_window=CAL_DEFAULT_SMOOTH_EVENTS,
                     roster=None, flush_mask=None, cal_mask=None,
                     exclude_mask=None, pressure=None,
                     pressure_target=D1_P_TARGET_MBARS):
    """Build a time-varying two-point calibration from the per-injection cal
    means and apply it to every row of df[value_col].

    Uses cal_mean_points' output directly as the drift-model nodes -- the
    per-injection averaging is not redone here.

    Returns a dict:
      ok            bool                -- whether a calibration was derived
      reason        str|None            -- if not, a displayable sentence
      mode          "two-point"|"offset"
      bottles       {state: {...}}      -- cal_bottle_series output, plus
                                           "nodes": (times, values) post-model
                                           and "loo" residuals
      low_state/high_state  int|None
      slope, intercept, calibrated      -- Series on df.index, NaN where undefined
      pressure_factor  Series|None      -- the applied pressure correction,
                                           None when it was not asked for
      extrapolated  Series[bool]        -- outside the INTERSECTION of the two
                                           bottles' node spans, or in a long gap
      flushed       Series[bool]        -- post-cal flush rows
      in_cal        Series[bool]        -- rows inside a cal period
      excluded      Series[bool]        -- warm-up / out-of-spec pressure rows
      non_ambient   Series[bool]        -- flushed | in_cal
      blanked       Series[bool]        -- non_ambient | excluded; NaN in
                                           `calibrated` (see below)
      residuals     [(time, state, closure, loo)]
      loo_rms       {state: float}      -- QC scatter, gas units
      span_gain     float|None
      warnings      [str]

    `calibrated` is the calibrated *good ambient* record. Three masks are
    blanked from it as the last step: `cal_mask` (rows inside a cal period)
    and `flush_mask` (post_cal_flush_mask), which are not air; and
    `exclude_mask` (warm-up + out-of-spec detector pressure), which is air the
    instrument was in no state to measure. None of the three affects the
    calibration itself -- the nodes, slope, intercept and residuals are all
    derived before this point.

    `exclude_mask` here is the *same* mask that is handed to
    `cal_mean_points`, but the two uses are independent and must stay that
    way: there it drops raw rows before the cal means are estimated (so it can
    remove a cal point); here it only blanks the output series. Passing it
    here does not and must not change `cal_points`, `loo_rms` or `span_gain`.
    Removing these rows here rather than in each caller keeps the plotted
    trace and the exported CSV from ever disagreeing about which rows are good.

    `pressure` (a Series of the gas's own detector pressure, in mbar) turns on
    the pressure correction: the calibrated value is scaled by
    `pressure_target / P`, normalising each reading to the cell's spec
    pressure. It is applied to the *calibrated* series only, deliberately --
    the cal-bottle responses, the drift nodes, the slope/intercept and every
    residual are the uncorrected measurement, so turning the correction on or
    off cannot move the calibration itself, only its product. (Correcting the
    raw signal first would also be defensible, but the two-point calibration
    would then absorb most of it through the bottle responses, which is a
    different and much less legible operation.)

    Rows with no usable pressure (missing, or <= 0) get no corrected value:
    they go NaN in `calibrated` and join `blanked`, so the trace breaks over
    them rather than mixing corrected and uncorrected values in one series.
    """
    import statistics

    nan_series = pd.Series(float("nan"), index=df.index)
    false_series = pd.Series(False, index=df.index)

    def failed(reason):
        return {
            "ok": False, "reason": reason, "mode": None, "bottles": {},
            "low_state": None, "high_state": None,
            "slope": nan_series, "intercept": nan_series, "calibrated": nan_series,
            "pressure_factor": None,
            "extrapolated": false_series, "flushed": false_series,
            "in_cal": false_series, "excluded": false_series,
            "non_ambient": false_series, "blanked": false_series,
            "residuals": [], "loo_rms": {},
            "span_gain": None, "warnings": [],
        }

    if not cal_points:
        return failed("No cal events available for this gas with the current masking.")

    bottles = cal_bottle_series(cal_points, cal_bottles, gas_key)
    for state, info in bottles.items():
        info["nodes"] = drift_nodes(info["times"], info["values"], model, smooth_window)
        info["loo"] = loo_residuals(info["times"], info["values"])

    usable = [s for s, i in bottles.items() if i["assigned"] is not None and i["times"]]
    if not usable:
        serials = sorted({i["serial"] for i in bottles.values() if i["serial"]})
        return failed(
            f"No assigned {gas_key} value in cals.yaml for the cal tanks "
            f"in use ({', '.join(serials) or 'none identified'})."
        )

    warnings = list(cal_mismatch_notes(bottles, gas_key, roster))
    times = df["datetime"]

    # Two distinct bottles with distinct assigned values are needed for a
    # slope; anything less degrades to an offset-only correction rather than
    # refusing outright, but says so loudly.
    distinct = {bottles[s]["assigned"] for s in usable}
    two_point = len(usable) >= 2 and len(distinct) >= 2
    if two_point:
        low_state = min(usable, key=lambda s: bottles[s]["assigned"])
        high_state = max(usable, key=lambda s: bottles[s]["assigned"])
        a_lo, a_hi = bottles[low_state]["assigned"], bottles[high_state]["assigned"]
        r_lo = interp_hold(*bottles[low_state]["nodes"], times)
        r_hi = interp_hold(*bottles[high_state]["nodes"], times)

        d_r = r_hi - r_lo
        # Relative degeneracy test -- a ppb gas has a different natural scale
        # than a ppm one, so an absolute epsilon would be wrong for one of them.
        scale = pd.concat([r_hi.abs(), r_lo.abs()], axis=1).max(axis=1).clip(lower=1.0)
        slope = (a_hi - a_lo) / d_r.where(d_r.abs() >= 1e-6 * scale)
        intercept = a_lo - slope * r_lo
        span_states = [low_state, high_state]
        mean_lo = statistics.mean(bottles[low_state]["values"])
        mean_hi = statistics.mean(bottles[high_state]["values"])
        span_gain = (mean_hi - mean_lo) / (a_hi - a_lo) if a_hi != a_lo else None
        mode = "two-point"
    else:
        state = usable[0]
        low_state = high_state = state
        r_lo = interp_hold(*bottles[state]["nodes"], times)
        slope = pd.Series(1.0, index=df.index)
        intercept = bottles[state]["assigned"] - r_lo
        span_states = [state]
        span_gain = None
        mode = "offset"
        warnings.append(
            "Only one usable cal bottle -- applying an offset-only correction "
            "(no gain/span term). Results are less reliable away from that "
            "bottle's concentration."
        )

    calibrated = slope * df[value_col] + intercept

    # Pressure correction, applied to the calibrated product and nothing else:
    # everything above -- the bottle responses, the nodes, slope/intercept --
    # is left on the uncorrected measurement, so this cannot move the
    # calibration. A pressure of 0 or NaN yields no factor and therefore no
    # value, which is why the resulting NaNs join `blanked` below.
    pressure_factor = None
    if pressure is not None:
        p = pd.to_numeric(pressure.reindex(df.index), errors="coerce")
        pressure_factor = pressure_target / p.where(p > 0)
        calibrated = calibrated * pressure_factor

    # The trustworthy region is the INTERSECTION of the bottles' node spans,
    # not "first to last cal event": one bottle can lose its points to masking
    # while the other keeps going, and the calibration is only bracketed where
    # both are live.
    span_start = max(bottles[s]["times"][0] for s in span_states)
    span_end = min(bottles[s]["times"][-1] for s in span_states)
    extrapolated = (times < span_start) | (times > span_end)
    for state in span_states:
        for gap_start, gap_end in _gap_spans(bottles[state]["times"]):
            extrapolated |= (times > gap_start) & (times < gap_end)
    extrapolated |= slope.isna()

    # Closure residual: what the calibration makes of its own input. Exact
    # zero by construction under "linear" -- a non-zero value there is a bug
    # signal, not a quality metric (that's what loo is for).
    residuals, loo_rms = [], {}
    for state in span_states:
        info = bottles[state]
        node_slope = interp_hold(times, slope, pd.Series(info["times"]))
        node_icept = interp_hold(times, intercept, pd.Series(info["times"]))
        for i, (t, value) in enumerate(zip(info["times"], info["values"])):
            closure = node_slope.iloc[i] * value + node_icept.iloc[i] - info["assigned"]
            residuals.append((t, state, closure, info["loo"][i]))
        finite = [r for r in info["loo"] if pd.notna(r)]
        if finite:
            loo_rms[state] = (sum(r * r for r in finite) / len(finite)) ** 0.5

    if slope.notna().any() and (slope.min() < CAL_SLOPE_SANE_RANGE[0]
                                 or slope.max() > CAL_SLOPE_SANE_RANGE[1]):
        warnings.append(
            f"Calibration slope ranges {slope.min():.3f}-{slope.max():.3f}, "
            f"outside the usual {CAL_SLOPE_SANE_RANGE[0]}-{CAL_SLOPE_SANE_RANGE[1]}. "
            "Data is kept, not discarded -- check the tank assignment."
        )

    # Applied last, and only to the output series: `calibrated` is the
    # calibrated *good ambient* record, so everything that isn't usable
    # atmosphere is blanked -- the cal gas itself, the flush behind it while
    # the cells still hold it, and the warm-up / out-of-spec-pressure rows the
    # masking controls already exclude from the cal means. Nothing is lost by
    # this: slope/intercept are still on every row, so a consumer who wants
    # the calibrated value of a blanked row (to check closure, say) can
    # recompute it.
    #
    # Note this runs after everything above: cal_points, the drift nodes,
    # residuals and span_gain are all already fixed, so `exclude_mask` cannot
    # reach back and change the calibration by being passed here.
    #
    # The rows stay in the frame as NaN rather than being dropped, so the
    # trace breaks visibly instead of interpolating over the removed stretch.
    def _align(mask):
        if mask is None:
            return false_series
        return mask.reindex(df.index).fillna(False).astype(bool)

    flushed = _align(flush_mask)
    in_cal = _align(cal_mask)
    excluded = _align(exclude_mask)
    non_ambient = flushed | in_cal
    blanked = non_ambient | excluded
    if pressure_factor is not None:
        # A row the correction could not be computed for is already NaN in
        # `calibrated`; joining `blanked` is what makes the trace break over
        # it instead of the row being dropped and drawn across.
        blanked = blanked | pressure_factor.isna()
    if blanked.any():
        calibrated = calibrated.mask(blanked)

    return {
        "ok": True, "reason": None, "mode": mode, "bottles": bottles,
        "low_state": low_state, "high_state": high_state,
        "slope": slope, "intercept": intercept, "calibrated": calibrated,
        "pressure_factor": pressure_factor,
        "extrapolated": extrapolated, "flushed": flushed,
        "in_cal": in_cal, "excluded": excluded,
        "non_ambient": non_ambient, "blanked": blanked,
        "residuals": residuals, "loo_rms": loo_rms,
        "span_gain": span_gain, "warnings": warnings,
    }


def linear_fit(x, y):
    """Ordinary least squares of y on x, as {"n", "slope", "intercept",
    "slope_err", "r"}, or None if there is nothing to fit.

    Deliberately plain OLS, not a fit that uses the error bars: the
    calibration uncertainties are almost entirely *systematic* (they shift a
    whole flight together), so weighting by them would not do what weighting
    is for -- and it would make the slope depend on whether a display toggle
    happens to be on. Slope error is the usual residual-based standard error
    and so describes the scatter about the line, nothing else.
    """
    x, y = pd.Series(x).astype(float), pd.Series(y).astype(float)
    keep = x.notna() & y.notna()
    x, y = x[keep], y[keep]
    n = len(x)
    if n < 3:
        return None
    dx, dy = x - x.mean(), y - y.mean()
    sxx = float((dx * dx).sum())
    syy = float((dy * dy).sum())
    sxy = float((dx * dy).sum())
    if sxx <= 0 or syy <= 0:
        return None
    slope = sxy / sxx
    residual_ss = max(syy - slope * sxy, 0.0)
    return {
        "n": n,
        "slope": slope,
        "intercept": float(y.mean() - slope * x.mean()),
        "slope_err": (residual_ss / (n - 2) / sxx) ** 0.5,
        "r": sxy / (sxx * syy) ** 0.5,
    }


def calibration_uncertainty(result):
    """1-sigma uncertainty on each calibrated value, propagated through the
    calibration this `result` describes.

    Returns (sigma, components): a Series on the same index as
    `result["calibrated"]` (NaN wherever that is), and a dict of the scalar
    inputs, for display.

    Writing the two-point calibration as a weighted blend of the two assigned
    values makes the propagation fall out. With
    f = (c - A_lo) / (A_hi - A_lo)  -- where a point sits between the bottles --

        c = (1-f) A_lo + f A_hi
        var(c) = ((1-f) sA_lo)^2 + (f sA_hi)^2
               + (slope (1-f) sR_lo)^2 + (slope f sR_hi)^2

    f is recovered from the calibrated value itself rather than from the
    interpolated responses, which keeps this a pure function of what
    calibrate_series already returns. Note f is not clamped to [0, 1]:
    ambient air usually sits outside the bracket the two tanks span, and the
    uncertainty genuinely does grow as you extrapolate away from them.

    The two input uncertainties per bottle:

    - **sA**, the assigned value's uncertainty, straight from cals.yaml's
      `<GAS>_unc`. Missing means the roster does not state one, and it is
      treated as 0 -- an underestimate, but inventing a number would be worse.
      This term is almost entirely *systematic*: it shifts a whole flight
      together rather than scattering point to point.
    - **sR**, how well the drift model reproduces that bottle's response,
      as sqrt(loo^2 + closure^2). Both halves are needed and neither alone is
      right. The leave-one-out RMS is the honest scatter under a model that
      interpolates through every node (`linear`, where closure is 0 by
      construction). The closure RMS is the error the model *itself*
      introduces by not passing through the nodes -- 0 under `linear`, but
      the dominant term under `constant`. So this is the answer to "does the
      uncertainty depend on the calibration method": it does, through here.

    What this does NOT include, deliberately: the single-sample noise of the
    ambient measurement (nothing in the calibration constrains it -- it would
    have to come from the raw trace's high-frequency scatter), and any
    inflation over the extrapolated spans, which are flagged separately by
    `result["extrapolated"]` rather than being folded into a number that would
    then look like ordinary uncertainty.
    """
    calibrated = result["calibrated"]
    if not result.get("ok"):
        return pd.Series(float("nan"), index=calibrated.index), {}

    bottles, slope = result["bottles"], result["slope"]

    # The pressure correction is a pure scaling of the calibrated value, so it
    # is divided back out before f is recovered (the blend of the two assigned
    # values only holds on the uncorrected scale) and multiplied back into the
    # answer at the end. The factor itself is treated as exact: the detector
    # pressure carries its own measurement error, but nothing in the
    # calibration constrains it, and inventing a number for it here would be
    # the same mistake as inventing a missing <GAS>_unc.
    factor = result.get("pressure_factor")
    uncorrected = calibrated if factor is None else calibrated / factor

    def rescale(sigma):
        return sigma if factor is None else sigma * factor

    closure_rms = {}
    for _, state, closure, _ in result.get("residuals", []):
        closure_rms.setdefault(state, []).append(closure)
    closure_rms = {state: (sum(c * c for c in v) / len(v)) ** 0.5
                   for state, v in closure_rms.items() if v}

    def response_sigma(state):
        loo = result.get("loo_rms", {}).get(state, 0.0) or 0.0
        closure = closure_rms.get(state, 0.0) or 0.0
        return (loo * loo + closure * closure) ** 0.5

    low, high = result["low_state"], result["high_state"]
    components = {
        "mode": result["mode"],
        "assigned_unc": {s: bottles[s].get("assigned_unc") for s in {low, high}},
        "response_sigma": {s: response_sigma(s) for s in {low, high}},
    }

    if result["mode"] == "offset":
        # c = m + (A - R(t)): no gain term, so the two uncertainties simply
        # add in quadrature and the result is the same for every row.
        s_a = bottles[low].get("assigned_unc") or 0.0
        s_r = response_sigma(low)
        sigma = pd.Series((s_a * s_a + s_r * s_r) ** 0.5, index=calibrated.index)
        return rescale(sigma).mask(calibrated.isna()), components

    a_lo, a_hi = bottles[low]["assigned"], bottles[high]["assigned"]
    s_a_lo = bottles[low].get("assigned_unc") or 0.0
    s_a_hi = bottles[high].get("assigned_unc") or 0.0
    s_r_lo, s_r_hi = response_sigma(low), response_sigma(high)

    f = (uncorrected - a_lo) / (a_hi - a_lo)
    g = 1.0 - f
    var = ((g * s_a_lo) ** 2 + (f * s_a_hi) ** 2
           + (slope * g * s_r_lo) ** 2 + (slope * f * s_r_hi) ** 2)
    return rescale(var.pow(0.5)).mask(calibrated.isna()), components


# --------------------------------------------------------------------------
# Exports
#
# Two products, from the same per-gas "blocks" (see export_companion_csv for
# the structure), because they answer different questions:
#
#   * the companion CSV is for working with the flight -- every row of the raw
#     file, every mask, both the raw and the derived value, aimed at Excel and
#     Igor Pro;
#   * the ICARTT file is for delivering it -- the archive format, good ambient
#     data only, everything else expressed as the format's missing value.
#
# Both live here rather than in the GUI so the CLI and any future batch script
# can produce byte-identical files, the same argument that put the masking and
# calibration here.
# --------------------------------------------------------------------------

# -7777 (above the upper limit of detection) and -8888 (below the lower) are
# fixed by the ICARTT standard. Nothing in this pipeline emits them -- no
# detection limit is established for these instruments -- but the header must
# declare them, so they are named here rather than being literals in the
# writer.
#
# The missing-data value is NOT fixed by the standard; it is declared per
# variable on header line 12, and a reader honours whatever is declared. -99999
# rather than the more common -9999 to match the sister UCATS instrument's
# delivered files (SABRE-UCATS-GC_WB57_20230303_R0.ict), so that anyone
# processing both instruments' files for a campaign meets one convention.
ICARTT_MISSING = -99999
ICARTT_ULOD_FLAG = -7777
ICARTT_LLOD_FLAG = -8888

# Mixing-ratio units are written with the trailing "v" in the delivered UCATS
# files (ppmv/ppbv/pptv). The plot labels keep the shorter form, so the
# translation happens on the way into the ICARTT header only.
ICARTT_UNITS = {"ppm": "ppmv", "ppb": "ppbv", "ppt": "pptv"}

# File Format Index 1001: ASCII, one independent variable (time), no bounded
# or auxiliary variables. What airborne 1 Hz in-situ data is archived as.
ICARTT_FFI = 1001

# Everything the ICARTT header carries that cannot be derived from the data.
# Blank rather than plausibly pre-filled: a wrong-but-reasonable PI name or
# mission would travel into an archived file unnoticed, whereas a blank one
# is visible to whoever reads the file. The two exceptions describe this
# repo's instrument rather than any flight, so they can be stated safely.
DEFAULT_ICARTT_META = {
    "data_id": "RASTA",
    "location_id": "",
    "pi_name": "",
    "pi_affiliation": "",
    "data_source": "RASTA airborne in-situ trace gas measurements",
    "mission": "",
    "pi_contact_info": "",
    "platform": "",
    "location": "",
    "associated_data": "N/A",
    "instrument_info": "",
    "data_info": "",
    "uncertainty": "",
    "ulod_value": "N/A",
    "llod_value": "N/A",
    "dm_contact_info": "",
    "project_info": "",
    "stipulations_on_use": "",
    "other_comments": "",
    "revision": "R0",
    # Two fields written VERBATIM, line breaks and all (see _verbatim_lines).
    # `revision_history` is one `R#: description` per line and **accumulates**
    # across revisions -- the delivered R0 file carries both `RA: Preliminary
    # data` and `R0: Revised data`, so this is a block the user maintains, not
    # a note about the current revision alone.
    "revision_history": "",
    # The special-comments section, in full: free text, and in the delivered
    # files it is where the error-estimate explanation and the "contact the
    # PIs" request live. Blank lines inside it are meaningful and preserved.
    # Nothing is appended to it -- see export_icartt.
    "special_comments": "",
    "var_suffix": "RASTA",
}

# Order matters: the ICARTT standard requires these keywords, one per line, in
# this sequence, at the head of the normal comments. (label, meta key) -- the
# LOD flag lines are constants and are inserted by the writer.
ICARTT_KEYWORDS = (
    ("PI_CONTACT_INFO", "pi_contact_info"),
    ("PLATFORM", "platform"),
    ("LOCATION", "location"),
    ("ASSOCIATED_DATA", "associated_data"),
    ("INSTRUMENT_INFO", "instrument_info"),
    ("DATA_INFO", "data_info"),
    ("UNCERTAINTY", "uncertainty"),
)


def _one_line(text):
    """Collapse to a single line.

    The format is line-oriented and the header declares its own line count,
    so a newline pasted into a metadata box would invalidate the file. Every
    user-typed value goes through here on the way out. Commas are left alone:
    the header lines and the keyword values are free text, and "Dutton,
    Geoff" is precisely the form the standard asks the PI name to take.
    """
    if text is None:
        return ""
    return " ".join(str(text).split())


def _field(text):
    """One line, and no commas either -- for the one place a value sits in a
    comma-delimited position: the `name, unit, standard_name, description`
    variable definition lines. A comma in a description there would read as an
    extra field to anything splitting the line."""
    return _one_line(text).replace(",", ";")


def _verbatim_lines(text):
    """Split free-form multi-line metadata into the lines it will occupy.

    For the two blocks the delivered files maintain by hand -- the special
    comments and the revision history -- where the line breaks (and the blank
    line inside the special comments) are the author's, not noise to be
    collapsed. Only trailing whitespace goes; the section's declared line
    count is taken from the length of this list, so the two cannot disagree.
    """
    if not text or not str(text).strip():
        return []
    return [line.rstrip() for line in str(text).replace("\r\n", "\n").split("\n")]


def icartt_time_base(datetimes, skip_leading=0):
    """(start date, seconds-from-midnight Series, usable-row mask).

    Seconds are counted from midnight UTC on the date of the first usable
    sample and keep counting past 86400 for a flight that crosses midnight,
    which is what the format asks for -- wrapping would make time run
    backwards.

    The mask is the rows ICARTT can actually represent: its independent
    variable must increase strictly, and this record does not. The datetime
    column contains duplicate timestamps (1435 of them in the Feb 2025 file,
    the same duplication interp_hold has to work around), and a duplicate or
    backward step makes the file invalid rather than merely untidy -- so the
    offending rows are identified here and reported to the user by the
    writer, instead of being discovered by an archive's validator later.

    **`skip_leading` is not optional in practice: pass the number of rows
    drop_presync_rows removed.** Those rows carry the stale pre-sync clock
    readings, which run *ahead* of the true time; left in, they set the
    running maximum hours into the future and every genuine row after the
    backward jump then fails the strictly-increasing test. On the Feb 2025
    file that is the difference between 1435 rows rejected and 5013.
    """
    times = pd.to_datetime(pd.Series(list(datetimes))).reset_index(drop=True)
    valid = times.notna()
    if skip_leading:
        valid.iloc[:skip_leading] = False
    if not valid.any():
        return None, pd.Series(dtype=float), pd.Series(dtype=bool)
    midnight = times[valid].iloc[0].normalize()
    # NaN on the skipped rows rather than a number: their clock was wrong, so
    # a seconds-from-midnight value for them would be a plausible-looking lie.
    seconds = (times - midnight).dt.total_seconds().where(valid)
    # cummax over the valid values gives the largest second count seen so far;
    # shifted, that is the value the current row has to beat. Invalid rows
    # drop out and leave the running maximum where it was.
    prev_max = seconds.cummax().shift().fillna(float("-inf"))
    keep = valid & (seconds > prev_max)
    return midnight.date(), seconds, keep


# `<gas>_cal` for a calibrated gas, `<gas>_filtered` for one that only has its
# physical floor applied. The two are not the same claim and must not share a
# column name -- a reader who sorts on "_cal" must not pick up a gas that was
# never calibrated. Keyed off final_kind through one function so the column
# the notes describe is always the column that was written.
FINAL_COLUMN_SUFFIX = {"calibrated": "cal", "filtered": "filtered"}


def _column_for(block):
    return f"{block['gas']}_{FINAL_COLUMN_SUFFIX[block['final_kind']]}"


def companion_notes(gas_blocks, source_path=None, presync_rows=0):
    """The provenance block for the companion CSV: what produced it, what the
    columns mean, and how many rows each mask covers."""
    import datetime as _dt

    notes = [
        "UCATS-B derived data - companion to the raw acquisition CSV.",
        f"generated: {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"source: {source_path or 'unknown'}",
        "One row per row of the source CSV, in the same order, so the two "
        "files can be opened side by side or pasted together without any "
        "alignment step.",
        "Mask columns are 1 = true, 0 = false, blank = not determined.",
    ]
    if presync_rows:
        notes.append(
            f"The first {presync_rows} row(s) were recorded before the "
            f"datalogger's clock synced and are excluded from every analysis; "
            f"they are kept here, blank and flagged in `presync`, only so the "
            f"row numbering matches the source file.")
    for block in gas_blocks:
        gas = block["gas"]
        unit = block.get("unit", "")
        if block.get("final_kind") == "calibrated":
            notes.append(
                f"{gas}: {_column_for(block)} is the calibrated GOOD "
                f"AMBIENT record ({unit}) -- blank wherever the instrument was "
                f"not sampling air, or was sampling it in a state the masking "
                f"settings exclude. The raw {block['value_col']} is untouched "
                f"on every row, and {gas}_cal_slope/{gas}_cal_intercept are "
                f"given everywhere, so a blank row can be recomputed.")
            if block.get("pressure_corrected"):
                # Stated as its own sentence, not folded into the one above:
                # it changes what the number IS, and the slope/intercept
                # recipe just given does not reproduce it without this step.
                notes.append(
                    f"{gas}: {_column_for(block)} is pressure-corrected -- "
                    f"scaled by {D1_P_TARGET_MBARS:.0f}/"
                    f"{block.get('pressure_col', 'P')} to normalise every "
                    f"value to the detector's {D1_P_TARGET_MBARS:.0f} mbar "
                    f"spec pressure. Applied after the calibration, so "
                    f"slope*raw+intercept gives the UNcorrected value.")
        elif block.get("final_kind") == "filtered":
            notes.append(
                f"{gas}: {_column_for(block)} is {block['value_col']} "
                f"({unit}) with readings below the sensor's physical floor "
                f"removed. This gas has no cal bottles, so it is not "
                f"calibrated here.")
        else:
            notes.append(f"{gas}: not calibrated -- {block.get('reason', 'no calibration')}")
        for name, mask in block.get("masks", {}).items():
            if mask is not None and bool(mask.any()):
                notes.append(f"{gas}_{name}: {int(mask.sum())} row(s)")
    return notes


def export_companion_csv(path, datetimes, gas_blocks, source_path=None,
                         include_raw=True, include_masks=True,
                         include_coefficients=True, include_uncertainty=True,
                         presync_rows=0, comment_header=False,
                         time_seconds=None):
    """Write the derived record for every gas, one row per row of the RAW CSV.

    A companion to the acquisition file rather than a replacement: it adds the
    masks and the filtered/calibrated values, repeats the raw columns only on
    request, and is the same length as the file it complements, so the two
    open side by side in Excel or Igor with no alignment step. That row-for-row
    promise is why the pre-sync rows dropped by drop_presync_rows are still
    present here, blank and flagged, rather than silently absent: a quiet
    offset of a few dozen rows between two files that look alignable is a very
    effective way to corrupt an analysis.

    Each entry of `gas_blocks` is a dict with `gas`, `value_col`, `unit`,
    `raw`, optionally `final`/`final_kind`/`sigma`/`slope`/`intercept`, a
    `masks` dict of named boolean Series, and `reason` when a gas that should
    have a calibration has none. **Every Series must already be indexed on the
    raw file's row numbering** -- reconciling the trimmed analysis frame with
    the untrimmed file is the caller's job (the GUI's `_to_raw_rows`), since
    only the caller knows how many rows came off the front.

    `comment_header` defaults False: the stated consumers are Excel and Igor
    Pro, and neither skips a leading `#` block without being told to. With it
    off the provenance goes to a sidecar `<stem>_notes.txt` instead, so it is
    never simply lost. Turn it on for a file meant to be read back with
    pd.read_csv(path, comment="#").
    """
    columns = {"datetime": pd.Series(list(datetimes))}
    if time_seconds is not None:
        # Seconds from midnight UTC, the same time base the ICARTT file uses.
        # Kept alongside the timestamp because Igor and Excel both plot a
        # numeric axis far more readily than a parsed date.
        columns["time_s"] = time_seconds
    if presync_rows:
        flag = pd.Series(0, index=columns["datetime"].index)
        flag.iloc[:presync_rows] = 1
        columns["presync"] = flag

    for block in gas_blocks:
        gas = block["gas"]
        if include_raw:
            columns[block["value_col"]] = block["raw"]
        if block.get("final") is not None:
            columns[_column_for(block)] = block["final"]
        if include_uncertainty and block.get("sigma") is not None:
            columns[f"{gas}_cal_unc"] = block["sigma"]
        if include_coefficients and block.get("slope") is not None:
            columns[f"{gas}_cal_slope"] = block["slope"]
            columns[f"{gas}_cal_intercept"] = block["intercept"]
        if include_masks:
            for name, mask in block.get("masks", {}).items():
                if mask is None:
                    continue
                # Int64, not bool: pandas would write True/False, which Igor
                # loads as a text wave and Excel as text. The nullable dtype
                # leaves the pre-sync rows genuinely empty rather than
                # claiming 0 for rows where the mask was never evaluated.
                columns[f"{gas}_{name}"] = mask.astype("boolean").astype("Int64")

    out = pd.DataFrame(columns)
    notes = companion_notes(gas_blocks, source_path, presync_rows)

    path = Path(path)
    with open(path, "w") as fh:
        if comment_header:
            fh.write("".join(f"# {line}\n" for line in notes))
        out.to_csv(fh, index=False)
    sidecar = None
    if not comment_header:
        sidecar = path.with_name(f"{path.stem}_notes.txt")
        sidecar.write_text("\n".join(notes) + "\n")
    return {"path": path, "notes_path": sidecar, "rows": len(out),
            "columns": list(out.columns)}


def icartt_filename(meta, start_date):
    """`dataID_locationID_YYYYMMDD_R#.ict`, the archive naming convention."""
    meta = {**DEFAULT_ICARTT_META, **(meta or {})}

    def slug(value, fallback):
        # Hyphens kept, underscores not: `_` separates the fields of the file
        # name itself, while a hyphenated data ID is normal and is what the
        # delivered UCATS files use (`SABRE-UCATS-GC_WB57_20230303_R0.ict`).
        # Stripping hyphens here silently mangled that into `SABREUCATSGC`.
        cleaned = "".join(c for c in str(value or "") if c.isalnum() or c == "-")
        return cleaned.strip("-") or fallback

    data_id = slug(meta["data_id"], "RASTA")
    location_id = slug(meta["location_id"], "Aircraft")
    revision = slug(meta["revision"], "R0")
    return f"{data_id}_{location_id}_{start_date:%Y%m%d}_{revision}.ict"


def _icartt_variables(gas_blocks, meta, include_sigma):
    """(name, definition-line, values) per dependent variable.

    Only the *final* series is offered -- calibrated where there is a
    calibration, floor-filtered where there is not. The raw counts are
    deliberately absent: ICARTT is the delivery format, and an uncalibrated
    column beside a calibrated one in an archived file is an invitation to
    plot the wrong one.

    The definition line is `name, unit, standard_name, description` -- four
    fields, matching the delivered UCATS files. A gas with no `standard_name`
    gets the field **empty rather than missing**: inventing an entry from the
    controlled vocabulary would be worse than admitting there isn't one, but
    dropping the field outright would shift the description into field 2 and
    make position mean different things on different lines of the same file.

    An uncertainty variable is named `<species>e_<suffix>` (`CO2e_RASTA`),
    the sister instrument's convention, and reuses its parent's standard name
    exactly as those files do.
    """
    suffix = "".join(c for c in str(meta.get("var_suffix") or "") if c.isalnum())
    variables = []

    def define(name, block, description):
        unit = block.get("unit", "")
        return ", ".join([name, ICARTT_UNITS.get(unit, unit),
                          _field(block.get("standard_name") or ""),
                          _field(description)])

    for block in gas_blocks:
        if block.get("final") is None:
            continue
        short = block.get("short", block["gas"])
        name = f"{short}_{suffix}" if suffix else short
        kind = ("calibrated against onboard reference cylinders"
                if block["final_kind"] == "calibrated"
                else "as measured, below-floor readings removed, uncalibrated")
        variables.append((
            name,
            define(name, block, f"{block.get('long_name', short)}; {kind}"),
            block["final"]))
        if include_sigma and block.get("sigma") is not None:
            error_name = f"{short}e_{suffix}" if suffix else f"{short}e"
            variables.append((
                error_name,
                define(error_name, block,
                       f"ERROR 1-sigma on {name}; propagated from the assigned "
                       f"cylinder values and the calibration scatter"),
                block["sigma"]))
    return variables


def _icartt_uncertainty_line(gas_blocks, meta):
    """The UNCERTAINTY keyword value: the median 1-sigma actually computed for
    each gas, with the user's own text after it.

    Derived rather than typed because it is the one required keyword whose
    answer this program already knows -- and a hand-typed uncertainty in an
    archived file goes stale the moment a drift model or a cal window changes.
    """
    parts = []
    for block in gas_blocks:
        sigma = block.get("sigma")
        if sigma is None or not sigma.notna().any():
            continue
        median = float(sigma.dropna().median())
        parts.append(f"{block.get('short', block['gas'])}: {median:.3g} "
                     f"{block.get('unit', '')} (median 1-sigma)")
    typed = _one_line(meta.get("uncertainty"))
    if typed:
        parts.append(typed)
    return "; ".join(parts) or "N/A"


def export_icartt(path, datetimes, gas_blocks, meta=None, include_sigma=True,
                  drop_empty_rows=True, skip_leading=0):
    """Write an ICARTT (.ict) file, format index 1001.

    The delivered record is the good ambient one: every row is written with
    the missing-value flag wherever `calibrate_series` blanked it, which is
    exactly what the format's -99999 means, so no separate mask columns are
    needed and none are written.

    **Everything in the file comes from `meta` or from the data.** Nothing
    about how the analysis was run is written -- no source file name, no
    masking or drift settings, no counts of what was dropped. Those are the
    experimenters' working record, not a data user's business, and they have
    their own homes (the flight's conf file, the companion CSV's notes, and
    the summary dict below).

    Returns a summary dict -- rows written, rows the time base could not
    represent, rows dropped as empty, and the variables emitted -- because
    the caller has to be able to tell the user that a file which validated
    cleanly is nonetheless shorter than the flight. That reporting is now the
    *only* place those counts appear.
    """
    import datetime as _dt

    meta = {**DEFAULT_ICARTT_META, **(meta or {})}
    # skip_leading: the pre-sync rows. See icartt_time_base -- omitting it
    # does not fail loudly, it just quietly rejects most of the flight.
    start_date, seconds, usable = icartt_time_base(datetimes, skip_leading)
    if start_date is None:
        raise ValueError("No usable timestamps: cannot build an ICARTT time base.")

    variables = _icartt_variables(gas_blocks, meta, include_sigma)
    if not variables:
        raise ValueError("No gas has a usable series to export.")

    frame = pd.DataFrame({name: pd.Series(list(values))
                          for name, _, values in variables})
    frame.insert(0, "Time_Start", seconds)
    keep = usable.copy()
    n_unusable = int((~keep).sum())
    n_empty = 0
    if drop_empty_rows:
        has_data = frame[[name for name, *_ in variables]].notna().any(axis=1)
        n_empty = int((keep & ~has_data).sum())
        keep &= has_data
    frame = frame[keep]
    if frame.empty:
        raise ValueError("Nothing left to write: no row has data on a usable "
                         "timestamp.")

    # Declare the interval only when it really is uniform. A gap-free 1 Hz
    # record gets "1"; anything else gets the format's "non-uniform" 0 rather
    # than a nominal rate the file does not actually keep to.
    diffs = frame["Time_Start"].diff().dropna().round(3).unique()
    interval = f"{diffs[0]:g}" if len(diffs) == 1 else "0"
    # Integer seconds are what a 1 Hz logger produces; only fall back to a
    # decimal format if some timestamp actually needs it.
    time_fmt = "%d" if float(frame["Time_Start"].mod(1).abs().max()) == 0 else "%.1f"

    today = _dt.date.today()
    revision = _one_line(meta["revision"]) or "R0"
    normal = [f"{label}: {_one_line(meta[key]) or 'N/A'}"
              for label, key in ICARTT_KEYWORDS]
    # UNCERTAINTY is the one keyword we can answer from the data.
    normal[-1] = f"UNCERTAINTY: {_icartt_uncertainty_line(gas_blocks, meta)}"
    normal += [
        f"ULOD_FLAG: {ICARTT_ULOD_FLAG}",
        f"ULOD_VALUE: {_one_line(meta['ulod_value']) or 'N/A'}",
        f"LLOD_FLAG: {ICARTT_LLOD_FLAG}",
        f"LLOD_VALUE: {_one_line(meta['llod_value']) or 'N/A'}",
        f"DM_CONTACT_INFO: {_one_line(meta['dm_contact_info']) or 'N/A'}",
        f"PROJECT_INFO: {_one_line(meta['project_info']) or 'N/A'}",
        f"STIPULATIONS_ON_USE: {_one_line(meta['stipulations_on_use']) or 'N/A'}",
        f"OTHER_COMMENTS: {_one_line(meta['other_comments']) or 'N/A'}",
        f"REVISION: {revision}",
    ]
    # The revision history is a block the user maintains, one `R#: text` per
    # line, and it ACCUMULATES: the delivered R0 file lists RA above R0. Only
    # if it is empty does the current revision get a placeholder line of its
    # own, so a first export is still a valid file.
    normal += _verbatim_lines(meta["revision_history"]) or [f"{revision}: N/A"]
    # The column-header line is itself the LAST normal comment line -- part of
    # the count on the line above it, not a separate section.
    normal.append(", ".join(["Time_Start"] + [name for name, *_ in variables]))

    # The author's special comments, verbatim, and NOTHING ELSE. This section
    # used to also carry the source file name, the masking/drift settings each
    # gas was analysed with, and the counts of rows dropped for the format's
    # sake. That is all internal to the experimenters: the settings in
    # particular are a working record of how the analysis was tuned, not
    # something a data user needs or should be reading in a delivered file.
    # None of it is lost -- the omission counts are returned in the summary
    # dict and shown to the user after the export, and the settings live in
    # the flight's own conf file and in the companion CSV's notes.
    special = _verbatim_lines(meta["special_comments"])

    after = [
        _one_line(meta["pi_name"]),
        _one_line(meta["pi_affiliation"]),
        _one_line(meta["data_source"]),
        _one_line(meta["mission"]),
        "1, 1",
        f"{start_date:%Y, %m, %d}, {today:%Y, %m, %d}",
        interval,
        # Three fields -- name, unit, description -- as the delivered UCATS
        # files write it, not two with the description crammed into the unit.
        "Time_Start, seconds, ELAPSED TIME IN SECONDS FROM 00:00:00 GMT ON "
        "THE FLIGHT DATE",
        str(len(variables)),
        ", ".join("1" for _ in variables),
        ", ".join(str(ICARTT_MISSING) for _ in variables),
    ]
    after += [definition for _, definition, _ in variables]
    after += [str(len(special))] + special
    after += [str(len(normal))] + normal

    # Line 1 counts itself, so the total is everything after it plus one.
    lines = [f"{len(after) + 1}, {ICARTT_FFI}"] + after

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        for row in frame.itertuples(index=False):
            values = [time_fmt % row[0]]
            for value in row[1:]:
                values.append(str(ICARTT_MISSING) if pd.isna(value)
                              else f"{value:.4f}")
            fh.write(", ".join(values) + "\n")

    return {"path": Path(path), "rows": len(frame), "header_lines": len(lines),
            "unusable_times": n_unusable, "empty_rows": n_empty,
            "start_date": start_date,
            "variables": [name for name, *_ in variables]}


def _shade_flagged(ax, datetimes, mask, color=CAL_SHADE_COLOR):
    """Hatch-shade extrapolated/untrustworthy spans. Hatching (rather than a
    plain fill) keeps these visually distinct from the timeseries plot's
    cal-period shading, which uses the same color for a different meaning."""
    for start, end in find_intervals(datetimes, mask):
        ax.axvspan(start, end, facecolor="none", edgecolor=color,
                   hatch="///", alpha=0.35, linewidth=0)


def _wrap_to_figure(fig, text, fontsize, usable=0.98):
    """Split `text` into lines that fit across `fig`, measured rather than
    guessed at.

    A fixed character count cannot be right for both a laptop and a large
    display: 110 characters left the header wrapping a sentence a few words
    early on a wide window while still overflowing a narrow one. The width of
    the string is measured with the figure's own renderer instead, so the
    header always uses the width it actually has. `usable` is the fraction of
    the figure the header may occupy -- the block is drawn at x=0.01, so it
    keeps the same margin at the other end.

    Falls back to the old fixed width when there is no renderer to ask (a
    Figure with no canvas), which is a degraded answer rather than a crash.
    """
    import textwrap

    get_renderer = getattr(fig.canvas, "get_renderer", None)
    if get_renderer is None:
        return textwrap.wrap(text, width=110)
    # Measured with a throwaway Text in the same font: matplotlib has no
    # public "how wide would this string be" that does not go through an
    # artist, and one probe per warning line is nothing beside the draw.
    probe = fig.text(0, 0, text, fontsize=fontsize, in_layout=False)
    width = probe.get_window_extent(get_renderer()).width
    probe.remove()
    available = fig.bbox.width * usable
    if width <= available:
        return [text]
    # Proportional font, so characters-per-line is derived from this string's
    # own average character width rather than from a nominal one.
    return textwrap.wrap(text, width=max(20, int(len(text) * available / width)))


def plot_calibration_panels(fig, result, gas_key, ylabel, datetimes, unit=""):
    """Draw the calibration diagnostics onto `fig`.

    Three stacked panels sharing a time axis: bottle response against
    assigned values, the derived slope/intercept, and residuals.
    """
    fig.clear()
    if not result.get("ok"):
        ax = fig.add_subplot(111)
        ax.set_axis_off()
        # in_layout=False: an unclipped Text is part of the Axes' tight bbox,
        # so under constrained_layout it sizes the very Axes it is positioned
        # and (with wrap=True) wrapped against, and the two chase each other a
        # little further on every redraw.
        ax.text(0.5, 0.5, result.get("reason") or "No calibration available.",
                ha="center", va="center", color=MUTED_COLOR, fontsize=11,
                wrap=True, transform=ax.transAxes, in_layout=False)
        return [ax]

    gs = fig.add_gridspec(3, 1, height_ratios=[2, 2, 1.5])
    ax_resp = fig.add_subplot(gs[0])
    ax_coef = fig.add_subplot(gs[1], sharex=ax_resp)
    ax_res = fig.add_subplot(gs[2], sharex=ax_resp)
    extrap = result["extrapolated"]

    state_colors = {result["low_state"]: CAL0_COLOR, result["high_state"]: CAL1_COLOR}

    # --- Panel 1: bottle response as a deviation from assigned ------------
    # Plotted as (measured - assigned) rather than absolute mole fraction:
    # the two tanks sit ~200 ppm apart, so on a shared absolute axis each
    # bottle's drift -- the whole point of this panel, and a couple of ppm at
    # most -- collapses to a flat line. Against zero, drift and event scatter
    # are legible and a wrong-tank offset shows up as a whole series sitting
    # far off the line. Absolute values stay in the legend.
    ax_resp.axhline(0.0, color=MUTED_COLOR, linestyle="--", linewidth=1.0)
    for state, info in sorted(result["bottles"].items()):
        color = state_colors.get(state, MUTED_COLOR)
        if not info["times"] or info["assigned"] is None:
            continue
        assigned = info["assigned"]
        measured = sum(info["values"]) / len(info["values"])
        delta = measured - assigned
        label = (f"{info['serial'] or f'state {state}'}  R={measured:.2f}  "
                 f"A={assigned:.3f}  d={delta:+.2f} ({delta / assigned:+.2%})")
        ax_resp.scatter(info["times"], [v - assigned for v in info["values"]],
                        color=color, s=28, zorder=5, edgecolors="none", label=label)
        node_t, node_v = info.get("nodes", ([], []))
        if len(node_t):
            ax_resp.plot(node_t, [v - assigned for v in node_v], color=color,
                         linewidth=1.2, zorder=4)
        for t, value, serial in info["rejected"]:
            ax_resp.scatter([t], [value - assigned], facecolors="none",
                            edgecolors=color, marker="X", s=55, zorder=6,
                            linewidths=1.2)
    ax_resp.set_ylabel(f"measured - assigned ({unit})" if unit else "measured - assigned",
                       color=TEXT_COLOR, fontsize=9)
    ax_resp.set_title("Cal bottle response, as deviation from assigned value",
                      color=TEXT_COLOR, loc="left", fontsize=10)
    ax_resp.legend(loc="upper left", fontsize=7, framealpha=0.9)

    # --- Panel 2: derived coefficients -----------------------------------
    _shade_flagged(ax_coef, datetimes, extrap)
    ax_coef.plot(datetimes, result["slope"], color=LINE_COLOR, linewidth=1.2)
    ax_coef.axhline(1.0, color=MUTED_COLOR, linestyle=":", linewidth=0.9)
    ax_coef.set_ylabel("slope", color=LINE_COLOR, fontsize=9)
    ax_coef.tick_params(axis="y", colors=LINE_COLOR, labelsize=8)
    ax_icept = ax_coef.twinx()
    # Dashed, and the only dashed trace on the panel: the two lines share one
    # frame but not one scale, so the reader has to keep track of which of
    # them the left axis is talking about. Colour alone did that; a second
    # cue does it in a greyscale print, and the dashes are the y-axis's
    # partner (the right spine and its labels carry the same colour).
    ax_icept.plot(datetimes, result["intercept"], color=RIGHT_AXIS_COLOR,
                  linewidth=1.2, linestyle="--", dashes=(5, 2))
    ax_icept.set_ylabel(f"intercept ({unit})" if unit else "intercept",
                        color=RIGHT_AXIS_COLOR, fontsize=9)
    ax_icept.tick_params(axis="y", colors=RIGHT_AXIS_COLOR, labelsize=8)
    ax_icept.spines["right"].set_color(RIGHT_AXIS_COLOR)
    ax_coef.set_title("Calibration coefficients (hatched = extrapolated)",
                      color=TEXT_COLOR, loc="left", fontsize=10)

    # --- Panel 3: residuals ----------------------------------------------
    ax_res.axhline(0.0, color=MUTED_COLOR, linewidth=0.9)
    for state, info in sorted(result["bottles"].items()):
        color = state_colors.get(state, MUTED_COLOR)
        pts = [(t, closure, loo) for t, s, closure, loo in result["residuals"] if s == state]
        if not pts:
            continue
        times = [p[0] for p in pts]
        rms = result["loo_rms"].get(state)
        ax_res.scatter(times, [p[1] for p in pts], color=color, s=18,
                       zorder=5, edgecolors="none", label=f"{info['serial']} closure")
        ax_res.scatter(times, [p[2] for p in pts], facecolors="none",
                       edgecolors=color, s=26, zorder=4, linewidths=1.0,
                       label=f"{info['serial']} leave-one-out"
                             + (f"  RMS={rms:.3f}" if rms is not None else ""))
        if info.get("assigned_unc"):
            ax_res.axhspan(-info["assigned_unc"], info["assigned_unc"],
                           color=color, alpha=0.15, linewidth=0)
    ax_res.set_ylabel(f"residual ({unit})" if unit else "residual",
                      color=TEXT_COLOR, fontsize=9)
    ax_res.set_title("Residuals: closure (filled) vs leave-one-out (hollow)",
                     color=TEXT_COLOR, loc="left", fontsize=10)
    ax_res.legend(loc="upper left", fontsize=7, framealpha=0.9, ncol=2)

    # --- Header block ------------------------------------------------------
    # Drawn as a figure suptitle rather than inside an Axes: these lines are
    # long enough to collide with the traces and legends at any in-axes anchor.
    # One font size for both the measuring and the drawing below: wrapping to
    # a width measured in a different size is wrapping to the wrong width.
    head_size = 8
    head = [f"{gas_key}  |  mode: {result['mode']}"]
    if result.get("span_gain") is not None:
        head[0] += f"  |  span gain: {result['span_gain']:.4f}"
    head[0] += f"  |  extrapolated: {100 * extrap.mean():.0f}% of record"
    head.append(", ".join(
        f"{info['serial'] or f'state {s}'}: {len(info['times'])} events"
        + (f", {len(info['rejected'])} rejected" if info["rejected"] else "")
        for s, info in sorted(result["bottles"].items())
    ))
    for warning in result["warnings"]:
        head.extend(_wrap_to_figure(fig, warning, head_size))
    # No explicit y: constrained_layout only reserves room for the suptitle
    # when it positions it itself, and this block runs to several lines.
    fig.suptitle("\n".join(head), color=MUTED_COLOR, fontsize=head_size,
                 ha="left", x=0.01, linespacing=1.4)
    if result["warnings"]:
        # Color the whole block as a warning -- matplotlib has no per-line
        # coloring in a suptitle, and a visible flag matters more than
        # keeping the neutral lines neutral.
        fig.texts[-1].set_color(PRESSURE_EXCLUDE_COLOR)

    for ax in (ax_resp, ax_coef, ax_res):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=GRID_COLOR, linewidth=0.6)
        for spine in ax.spines.values():
            spine.set_color(AXIS_COLOR)
        ax.tick_params(colors=MUTED_COLOR, labelsize=8)
    ax_resp.tick_params(labelbottom=False)
    ax_coef.tick_params(labelbottom=False)
    ax_res.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_res.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    for label in ax_res.get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment("right")
    return [ax_resp, ax_coef, ax_res, ax_icept]
