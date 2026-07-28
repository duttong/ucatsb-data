#!/usr/bin/env python3
"""Plot the d1_CO2_ppm (uncalibrated CO2) timeseries from a UCATS-B CSV file.

Usage: python3 plot_co2_timeseries.py <csv_file>
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

LINE_COLOR = "#2a78d6"
RIGHT_AXIS_COLOR = "#8e44ad"   # purple, distinct from the red/orange masking shades
CAL_SHADE_COLOR = "#898781"
PRESSURE_EXCLUDE_COLOR = "#d03b3b"
WARMUP_EXCLUDE_COLOR = "#ffa64d"   # light orange
POST_CAL_FLUSH_COLOR = "#2fa88a"   # teal, distinct from the warm-up/pressure shades
CAL0_COLOR = "#eda100"   # golden
CAL1_COLOR = "#0d366b"   # dark blue
GRID_COLOR = "#e1e0d9"
AXIS_COLOR = "#c3c2b7"
TEXT_COLOR = "#0b0b0b"
MUTED_COLOR = "#52514e"

D1_P_TARGET_MBARS = 140.0
D1_P_TOLERANCE_MBARS = 10.0
WARMUP_MINUTES = 30
CAL_MEAN_WINDOW_S = 15
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


def load_cal_bottles(path: Path):
    """Load the two plumbed cal bottles named in cals.yaml's `cals:` block,
    with their nominal gas concentrations, as {serial: {...}} (or {}).

    cals.yaml is a local copy of ~/code/ucats-b/cals.yaml -- resync it by hand
    if the roster changes. This is the default pairing only; see
    load_cal_assignment.
    """
    return select_cal_bottles(load_cal_roster(path), load_cal_assignment(path).values())


def load_cal_roster(path: Path):
    """Load EVERY tank in cals.yaml's roster, not just the two named in the
    `cals:` block.

    This exists solely for the mismatch advisory (see cal_mismatch_notes) --
    it must never feed bottle matching, which stays restricted to
    load_cal_bottles for the reason documented there.
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
                     exclude_mask=None):
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
    """
    import statistics

    nan_series = pd.Series(float("nan"), index=df.index)
    false_series = pd.Series(False, index=df.index)

    def failed(reason):
        return {
            "ok": False, "reason": reason, "mode": None, "bottles": {},
            "low_state": None, "high_state": None,
            "slope": nan_series, "intercept": nan_series, "calibrated": nan_series,
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
    if blanked.any():
        calibrated = calibrated.mask(blanked)

    return {
        "ok": True, "reason": None, "mode": mode, "bottles": bottles,
        "low_state": low_state, "high_state": high_state,
        "slope": slope, "intercept": intercept, "calibrated": calibrated,
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
        return sigma.mask(calibrated.isna()), components

    a_lo, a_hi = bottles[low]["assigned"], bottles[high]["assigned"]
    s_a_lo = bottles[low].get("assigned_unc") or 0.0
    s_a_hi = bottles[high].get("assigned_unc") or 0.0
    s_r_lo, s_r_hi = response_sigma(low), response_sigma(high)

    f = (calibrated - a_lo) / (a_hi - a_lo)
    g = 1.0 - f
    var = ((g * s_a_lo) ** 2 + (f * s_a_hi) ** 2
           + (slope * g * s_r_lo) ** 2 + (slope * f * s_r_hi) ** 2)
    return var.pow(0.5).mask(calibrated.isna()), components


def export_calibrated_csv(path, df, result, value_col, gas_key,
                          source_path=None, settings=None, analysis=None):
    """Write the calibrated series with its provenance and per-row flags.

    Every row is written, flagged, rather than pre-filtering to ambient: the
    consumer can filter, but silently dropping rows from a calibration
    product would hide what was excluded and why. The leading `#` comment
    block records what produced the numbers -- read it back with
    pd.read_csv(path, comment="#").
    """
    import datetime as _dt

    lines = [
        f"# UCATS-B calibrated {gas_key}",
        f"# generated: {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"# source: {source_path or 'unknown'}",
        f"# column: {value_col}  mode: {result['mode']}",
    ]
    for state, info in sorted(result["bottles"].items()):
        lines.append(
            f"# cal state {state}: serial={info['serial']} "
            f"assigned={info['assigned']} unc={info['assigned_unc']} "
            f"events={len(info['times'])} rejected={len(info['rejected'])}"
        )
    if result.get("span_gain") is not None:
        lines.append(f"# span gain: {result['span_gain']:.6f}")
    for state, rms in sorted(result.get("loo_rms", {}).items()):
        lines.append(f"# leave-one-out RMS, state {state}: {rms:.6f}")
    if settings:
        lines.append("# settings: " + ", ".join(f"{k}={v}" for k, v in sorted(settings.items())))
    for warning in result.get("warnings", []):
        lines.append(f"# WARNING: {warning}")
    lines.append("# is_extrapolated=True means the calibration was held flat "
                 "past the last cal event of a bottle, or spans a long gap.")
    flushed = result.get("flushed")
    in_cal = result.get("in_cal")
    lines.append(f"# {value_col}_cal is the calibrated GOOD AMBIENT record. It is "
                 f"blank wherever the instrument was not sampling air (is_cal_period "
                 f"or is_post_cal_flush), wherever the air it was sampling is "
                 f"excluded by the masking settings (is_masked: warm-up or "
                 f"out-of-spec detector pressure), and also wherever the raw "
                 f"{value_col} is itself missing or the calibration is undefined. The raw "
                 f"{value_col} is left untouched on every row, and cal_slope/"
                 f"cal_intercept are given for every row, so the calibrated "
                 f"value of a blanked row can be recomputed if it is wanted.")
    if in_cal is not None and in_cal.any():
        lines.append(f"# is_cal_period=True: cal gas was flowing "
                     f"({int(in_cal.sum())} rows), so it is not air.")
    if flushed is not None and flushed.any():
        lines.append(f"# is_post_cal_flush=True: within the post-cal flush window "
                     f"({int(flushed.sum())} rows) -- the detector was still "
                     f"clearing cal gas, so the air reading is not yet the atmosphere.")
    masked = result.get("excluded")
    if masked is not None and masked.any():
        lines.append(f"# is_masked=True: excluded by the masking settings "
                     f"({int(masked.sum())} rows) -- instrument warm-up or "
                     f"detector pressure outside tolerance.")

    out = pd.DataFrame({
        "datetime": df["datetime"],
        value_col: df[value_col],
        f"{value_col}_cal": result["calibrated"],
        "cal_slope": result["slope"],
        "cal_intercept": result["intercept"],
        "is_extrapolated": result["extrapolated"],
    })
    if flushed is not None:
        out["is_post_cal_flush"] = flushed
    # is_cal_period comes from the calibration result rather than the analysis
    # so it is present even without one -- the header explains a blank
    # <col>_cal in terms of this column, so it must not be able to go missing.
    if in_cal is not None:
        out["is_cal_period"] = in_cal
    # Same reasoning for is_masked: the header explains a blank <col>_cal in
    # terms of it, so it comes from the calibration result (which did the
    # blanking) and only falls back to the analysis dict.
    excluded = result.get("excluded")
    if excluded is None and analysis is not None:
        excluded = analysis["exclude_mask"]
    if excluded is not None:
        out["is_masked"] = excluded

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        out.to_csv(fh, index=False)
    return path


def _shade_flagged(ax, datetimes, mask, color=CAL_SHADE_COLOR):
    """Hatch-shade extrapolated/untrustworthy spans. Hatching (rather than a
    plain fill) keeps these visually distinct from the timeseries plot's
    cal-period shading, which uses the same colour for a different meaning."""
    for start, end in find_intervals(datetimes, mask):
        ax.axvspan(start, end, facecolor="none", edgecolor=color,
                   hatch="///", alpha=0.35, linewidth=0)


def plot_calibration_panels(fig, result, gas_key, ylabel, datetimes, unit=""):
    """Draw the calibration diagnostics onto `fig`.

    Three stacked panels sharing a time axis: bottle response against
    assigned values, the derived slope/intercept, and residuals.
    """
    fig.clear()
    if not result.get("ok"):
        ax = fig.add_subplot(111)
        ax.set_axis_off()
        ax.text(0.5, 0.5, result.get("reason") or "No calibration available.",
                ha="center", va="center", color=MUTED_COLOR, fontsize=11,
                wrap=True, transform=ax.transAxes)
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
    ax_resp.legend(loc="lower right", fontsize=7, framealpha=0.9)

    # --- Panel 2: derived coefficients -----------------------------------
    _shade_flagged(ax_coef, datetimes, extrap)
    ax_coef.plot(datetimes, result["slope"], color=LINE_COLOR, linewidth=1.2)
    ax_coef.axhline(1.0, color=MUTED_COLOR, linestyle=":", linewidth=0.9)
    ax_coef.set_ylabel("slope", color=LINE_COLOR, fontsize=9)
    ax_coef.tick_params(axis="y", colors=LINE_COLOR, labelsize=8)
    ax_icept = ax_coef.twinx()
    ax_icept.plot(datetimes, result["intercept"], color=RIGHT_AXIS_COLOR, linewidth=1.2)
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
    ax_res.legend(loc="upper right", fontsize=7, framealpha=0.9, ncol=2)

    # --- Header block ------------------------------------------------------
    # Drawn as a figure suptitle rather than inside an Axes: these lines are
    # long enough to collide with the traces and legends at any in-axes anchor.
    import textwrap

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
        head.extend(textwrap.wrap(warning, width=110))
    # No explicit y: constrained_layout only reserves room for the suptitle
    # when it positions it itself, and this block runs to several lines.
    fig.suptitle("\n".join(head), color=MUTED_COLOR, fontsize=8,
                 ha="left", x=0.01, linespacing=1.4)
    if result["warnings"]:
        # Colour the whole block as a warning -- matplotlib has no per-line
        # colouring in a suptitle, and a visible flag matters more than
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


def plot_co2_timeseries(csv_path: Path, out_path: Path):
    df = pd.read_csv(
        csv_path,
        usecols=["datetime", "d1_CO2_ppm", "d1_P_mbars", "j_sol_cals", "j_sol_aircal"],
    )
    df["datetime"] = pd.to_datetime(df["datetime"])
    df = drop_presync_rows(df)

    cal_bottles = load_cal_bottles(CALS_YAML_PATH)

    fig, ax = plt.subplots(figsize=(10, 5), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")

    # Shade calibration/cal-gas periods (j_sol_cals or j_sol_aircal == 1) so
    # they aren't mistaken for ambient air measurements.
    cal = (df["j_sol_cals"].fillna(0).astype(bool)
           | df["j_sol_aircal"].fillna(0).astype(bool))
    shade_intervals(ax, df["datetime"], cal, CAL_SHADE_COLOR, alpha=0.3)

    cal_intervals = merge_close_intervals(
        find_intervals(df["datetime"], cal), pd.Timedelta(seconds=CAL_MERGE_GAP_S)
    )

    # Flag periods where cell pressure (d1_P_mbars) is out of spec.
    bad_pressure = (df["d1_P_mbars"] - D1_P_TARGET_MBARS).abs() > D1_P_TOLERANCE_MBARS
    bad_pressure = bad_pressure.fillna(False)

    # Flag the instrument warm-up period at the start of the record.
    warmup_end = df["datetime"].iloc[0] + pd.Timedelta(minutes=WARMUP_MINUTES)
    warmup = df["datetime"] < warmup_end

    # Cal means are estimated from the raw data with these masks applied --
    # a cal point can be dropped entirely if its window has no valid data.
    exclude_mask = bad_pressure | warmup
    cal_points = cal_mean_points(
        df, cal_intervals, "d1_CO2_ppm",
        cal0_window=(-CAL_MEAN_WINDOW_S, -1), cal1_window=(-CAL_MEAN_WINDOW_S, -1),
        cal_bottles=cal_bottles, gas_key="CO2", exclude_mask=exclude_mask,
    )

    shade_intervals(ax, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
    shade_intervals(ax, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)

    plot_data = df[["datetime", "d1_CO2_ppm"]].dropna()
    line, = ax.plot(plot_data["datetime"], plot_data["d1_CO2_ppm"], color=LINE_COLOR, linewidth=1.2)

    cal0_pts = [(t, v) for t, v, state, serial in cal_points if state == 0]
    cal1_pts = [(t, v) for t, v, state, serial in cal_points if state == 1]
    cal0_label = most_common_serial(cal_points, 0) or "Cal 0"
    cal1_label = most_common_serial(cal_points, 1) or "Cal 1"

    handles = [line]
    labels = ["CO2 (ambient)"]
    if cal0_pts:
        xs, ys = zip(*cal0_pts)
        cal0_scatter = ax.scatter(xs, ys, color=CAL0_COLOR, s=40, zorder=5, edgecolors="none")
        handles.append(cal0_scatter)
        labels.append(f"{cal0_label}: {mean_std_label(ys)}")
    if cal1_pts:
        xs, ys = zip(*cal1_pts)
        cal1_scatter = ax.scatter(xs, ys, color=CAL1_COLOR, s=40, zorder=5, edgecolors="none")
        handles.append(cal1_scatter)
        labels.append(f"{cal1_label}: {mean_std_label(ys)}")

    ax.set_ylabel("CO2 (ppm)", color=TEXT_COLOR)
    ax.set_title("UCATS-B CO2 (uncalibrated) timeseries", color=TEXT_COLOR, loc="left")

    date_str = plot_data["datetime"].iloc[0].strftime("%Y-%m-%d") if not plot_data.empty else ""
    ax.set_xlabel(f"Time (UTC-ish, {date_str})", color=MUTED_COLOR)

    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fig.autofmt_xdate(rotation=45)

    ax.grid(True, color=GRID_COLOR, linewidth=0.8)
    for spine in ax.spines.values():
        spine.set_color(AXIS_COLOR)
    ax.tick_params(colors=MUTED_COLOR)

    ax.legend(handles, labels, loc="lower right", fontsize=9, framealpha=0.9)

    # Note shaded/flagged periods, if any occurred
    notes = []
    if cal.any():
        notes.append("gray = calibration/cal-air (j_sol_cals, j_sol_aircal)")
    if warmup.any():
        notes.append(f"orange = excluded (first {WARMUP_MINUTES} min warm-up)")
    if bad_pressure.any():
        notes.append(
            f"light red = excluded (d1_P_mbars outside {D1_P_TARGET_MBARS:.0f}±{D1_P_TOLERANCE_MBARS:.0f} mbar)"
        )
    if notes:
        ax.text(
            0.01, 0.98, "\n".join(notes),
            transform=ax.transAxes, ha="left", va="top",
            color=MUTED_COLOR, fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 plot_co2_timeseries.py <csv_file>", file=sys.stderr)
        sys.exit(1)
    csv_path = Path(sys.argv[1])
    out_path = csv_path.with_name(csv_path.stem + "_CO2_ppm.png")
    plot_co2_timeseries(csv_path, out_path)
