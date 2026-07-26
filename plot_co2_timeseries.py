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
CAL_SHADE_COLOR = "#898781"
PRESSURE_EXCLUDE_COLOR = "#d03b3b"
WARMUP_EXCLUDE_COLOR = "#ffa64d"   # light orange
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


def load_cal_bottles(path: Path):
    """Load cal bottle serials and their nominal gas concentrations from a
    cals.yaml file (a local copy of ~/code/ucats-b/cals.yaml -- resync it by
    hand if bottles change). Returns {serial: {"CO2": ..., "N2O": ..., ...}},
    or {} if unavailable.
    """
    import yaml
    if not path.exists():
        return {}
    try:
        # cals.yaml uses tabs as "key:\tvalue" separators in places, which
        # the YAML spec disallows outright; tabs never carry structural
        # indentation here, so blanking them out is safe.
        text = path.read_text().replace("\t", " ")
        data = yaml.safe_load(text) or {}
    except yaml.YAMLError as e:
        print(f"Warning: could not parse {path}: {e}")
        return {}
    serials = set(data.get("cals", {}).values())
    return {serial: data[serial] for serial in serials if serial in data}


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
