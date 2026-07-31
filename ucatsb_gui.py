#!/usr/bin/env python3
"""Interactive UCATS-B timeseries viewer with adjustable data-masking controls.

Left panel: a Load Data button, gas selector (CO2/N2O), and masking controls
(warm-up exclusion, detector pressure filter, and per-bottle calibration mean
windows). Right panel: the resulting figure.

Usage: python3 ucatsb_gui.py [csv_file]

The CSV can also be picked (or swapped out) from within the GUI via the
"Load Data" button, so the CLI argument is optional -- run with no argument
to start empty and load a file from the file browser.
"""
import copy
import functools
import html
import re
import sys
from pathlib import Path

import pandas as pd
import yaml
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QFontMetrics
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QGroupBox, QComboBox, QDoubleSpinBox, QSpinBox, QLabel, QGridLayout,
    QButtonGroup, QRadioButton, QPushButton, QFileDialog, QMessageBox,
    QTabWidget, QCheckBox, QAction, QStackedWidget, QMenu, QFrame, QStyle,
    QDialog, QDialogButtonBox, QLineEdit, QPlainTextEdit, QScrollArea,
    QListWidget, QListWidgetItem, QToolButton,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import matplotlib.dates as mdates

from ucatsb_analysis import (
    drop_presync_rows, find_intervals, merge_close_intervals,
    shade_intervals, bottle_for_interval, match_cal_serial, cal_mean_points,
    load_cal_roster, load_cal_assignment, select_cal_bottles,
    most_common_serial, mean_std_label, calibrate_series, interp_hold,
    post_cal_flush_mask,
    cal_switch_mask, below_floor_mask, smooth_pressure, pt_correction_factor,
    O3_VALID_MIN_PPB, H2O_VALID_MIN_PPM,
    box_stats, calibration_uncertainty, linear_fit,
    plot_calibration_panels,
    export_companion_csv, export_icartt, icartt_filename, icartt_time_base,
    DEFAULT_ICARTT_META,
    CALS_YAML_PATH, CAL_DRIFT_MODELS, CAL_DEFAULT_SMOOTH_EVENTS,
    CAL_MERGE_GAP_S, D1_P_TARGET_MBARS, T_GAS_TARGET_K, T_GAS_REFERENCE_C,
    T_GAS_REFERENCE_K, KELVIN_OFFSET,
    merge_ranges, add_ranges, subtract_ranges, ranges_to_mask, ranges_row_count,
    # The shared palette. These used to be declared a second time here, with
    # identical values -- two homes for one decision, agreeing only by
    # discipline, so editing one silently made the calibration panels and the
    # timeseries disagree about what a color means. ucatsb_analysis is the
    # single home; only the two GUI-only colors below are declared here.
    LINE_COLOR, RIGHT_AXIS_COLOR, CAL_SHADE_COLOR, PRESSURE_EXCLUDE_COLOR,
    WARMUP_EXCLUDE_COLOR, PUMPS_EXCLUDE_COLOR, POST_CAL_FLUSH_COLOR,
    CAL0_COLOR, CAL1_COLOR, FLAGGED_COLOR,
    GRID_COLOR, AXIS_COLOR, TEXT_COLOR, MUTED_COLOR,
)

# The calibrated overlay gets its own color rather than reusing LINE_COLOR:
# the raw trace stays blue when the overlay is on, so the two are told apart
# by hue, not by which one happens to be faded. Darker/more saturated than the
# 15%-alpha PRESSURE_EXCLUDE_COLOR band so it doesn't read as shading.
CALIBRATED_COLOR = "#c0392b"
STATS_BOX_COLOR = "#111111"

# The "armed" look for the two selector tools (Stats, Flag). A checked
# QToolButton's stock highlight is a grey barely darker than the toolbar, and
# these two are modes -- while one is on, a drag on the canvas does something
# other than nothing, so it has to be obvious at a glance which one is live.
# Light green fill with a dark green label: the pair keeps the text at ~9:1
# contrast, where the stock highlight leaves it at the toolbar's own.
TOOL_ON_BG = "#c8e6c9"
TOOL_ON_BG_HOVER = "#b2dfb4"
TOOL_ON_FG = "#1b5e20"

# `short`, `standard_name` and `long_name` exist for the exports: the ICARTT
# variable name is `<short>_<suffix>` (so Ozone is delivered as O3, the name
# every archive uses), `standard_name` is the controlled-vocabulary field the
# sister UCATS files carry as the third field of each variable definition
# (`Gas_N2O_InSitu_S_DMF`), and `long_name` is the description after it.
# `short`/`long_name` default sensibly from the gas key; `standard_name` has
# no default and is simply omitted from the line when absent, since inventing
# a controlled-vocabulary entry is worse than leaving it out.
#
# The names come from the NASA ESDS Atmospheric Composition Variable Standard
# Names Convention, which ICARTT V2.0 requires:
#   MeasurementCategory_CoreName_AcquisitionMethod_DescriptiveAttributes
# and for the `Gas` category the attributes are MeasurementSpecificity
# (S = single species) then Reporting (DMF = molar fraction wrt DRY air; AMF
# is the same wrt ambient air; DVMR/AVMR are the volumetric mixing ratios on
# those two bases; None = not stated). Two entries here are not `Gas_..._DMF`
# and neither is an oversight:
#
# - **H2O is not a `Gas` name at all.** Water vapour lives in the convention's
#   `Met` category, whose format is Met_CoreName_AcquisitionMethod_None with
#   no descriptive attributes. `H2OMF` is "mole fraction of water vapor with
#   respect to ambient air"; `H2OMRV` (volumetric mixing ratio) and `H2OMR`
#   (mass mixing ratio to dry air) are the alternatives if the sensor turns
#   out to be stated differently.
# - **Ozone reports `AVMR`, not `DMF`** -- a volumetric mixing ratio against
#   AMBIENT air, water vapour included. The 2B monitor measures the sample as
#   it comes in and nothing dries it, so a dry-air basis would be a claim
#   about the sample stream that this pipeline does not make. Settled with
#   Eric, 2026-07-29; it was briefly `DMF` and then `None`.
GASES = {
    "CO2": {"value_col": "d1_CO2_ppm", "ylabel": "CO2 (ppm)", "title": "UCATS-B CO2 (uncalibrated) timeseries", "detector": "d1", "standard_name": "Gas_CO2_InSitu_S_DMF", "long_name": "Carbon dioxide dry air mole fraction"},
    "N2O": {"value_col": "d1_N2O_ppb", "ylabel": "N2O (ppb)", "title": "UCATS-B N2O (uncalibrated) timeseries", "detector": "d1", "standard_name": "Gas_N2O_InSitu_S_DMF", "long_name": "Nitrous oxide dry air mole fraction"},
    "CH4": {"value_col": "d2_CH4_ppb", "ylabel": "CH4 (ppb)", "title": "UCATS-B CH4 (uncalibrated) timeseries", "detector": "d2", "standard_name": "Gas_CH4_InSitu_S_DMF", "long_name": "Methane dry air mole fraction"},
    # Ozone comes from its own dedicated sensor, not an Aeris detector, and
    # isn't run through the cal-bottle system -- has_masking=False skips the
    # warm-up/pressure-tol/cal-window machinery entirely for this gas, and
    # detector=None disables the Detector Pressure/T_gas aux traces (there's
    # no matching column to route to). These two are kept last so they sort
    # to the bottom of the Gas combo box, after the cal-bottle gases.
    # valid_min is a *physical floor*, not one of the maskable settings:
    # below it the sensor is faulting, not measuring. Declared per gas rather
    # than special-cased by name so the filtering code stays generic -- any
    # gas that gains a floor gets the same raw/filtered treatment.
    "Ozone": {"value_col": "oz_o3best", "ylabel": "O3 (ppb)", "title": "UCATS-B O3 timeseries", "detector": None, "has_masking": False, "valid_min": O3_VALID_MIN_PPB, "short": "O3", "standard_name": "Gas_O3_InSitu_S_AVMR", "long_name": "Ozone volume mixing ratio in ambient air"},
    # Water vapour, from its own instrument (`w_*` columns) -- like Ozone it
    # has no cal bottles, so has_masking=False and it is plotted as recorded.
    "H2O": {"value_col": "w_H2Obest", "ylabel": "H2O (ppm)", "title": "UCATS-B H2O timeseries", "detector": None, "has_masking": False, "valid_min": H2O_VALID_MIN_PPM, "standard_name": "Met_H2OMF_InSitu_None", "long_name": "Water vapour mole fraction in ambient air"},
}


def gas_unit(gas_key: str) -> str:
    """The unit out of a gas's ylabel ("CO2 (ppm)" -> "ppm"). One source for
    the plot labels, the cal-window box titles and both export formats, so a
    unit cannot be right on the figure and wrong in a delivered file."""
    ylabel = GASES[gas_key]["ylabel"]
    return ylabel.split("(")[-1].rstrip(")") if "(" in ylabel else ""

REQUIRED_COLUMNS = [
    "datetime", "d1_P_mbars", "d2_P_mbars", "d1_T_gas", "d2_T_gas",
    "j_sol_cals", "j_sol_aircal", "j_pumps",
] + [g["value_col"] for g in GASES.values()]

# Columns already exposed via a specific named control (gas traces, the
# named aux radio options) -- excluded from the "Other" catch-all combo box
# since picking them there would be redundant. This is a narrower set than
# REQUIRED_COLUMNS: j_sol_cals/j_sol_aircal are required for cal-interval
# detection but aren't plotted anywhere by name, so they stay selectable
# via "Other" (e.g. to sanity-check the raw digital flag against a trace).
# oz_o3 is likewise not a named control -- it's reached through "Other" like
# any other raw column (the Ozone *gas* trace is oz_o3best, which is named).
NAMED_TRACE_COLUMNS = {
    "d1_P_mbars", "d2_P_mbars", "d1_T_gas", "d2_T_gas",
} | {g["value_col"] for g in GASES.values()}

AUX_OPTIONS = ["No Figure", "Detector Pressure", "T_gas", "Other"]

# Content width of the left control panel. The scroll area around it reserves
# the scrollbar's width on top of this, so the controls get the full amount
# either way -- several comments below explain layout choices ("a stretching
# field squeezes the label column to nothing") that only make sense against a
# fixed, and fairly narrow, panel.
#
# 312 rather than the 300 it was for a long time: the widest group box now asks
# for 292 and there is nothing to warn you when one asks for more. A panel too
# narrow for its contents does not scroll or wrap, it silently clips its own
# right-hand edge -- the Save button and the Clear button lose their right
# halves and nothing says why. The dozen px is headroom against that, and costs
# ~1% of the figure's width. Anything widened here should be checked against it.
CONTROLS_WIDTH = 312

# The panel layouts' contents margin, 6 rather than Qt's 9 for the reasons
# above. Named because _elide_field has to subtract it to know how much width a
# label actually gets; every setContentsMargins in the two control panels uses
# it, so a change here must be a change there.
CONTROLS_MARGIN = 6


def aux_trace_info(selection: str, gas: str, other_column: str = None):
    """Return (column, ylabel) for the chosen auxiliary trace, given which
    detector the active gas comes from (CO2/N2O -> d1, CH4 -> d2 -- d2 used
    to carry a redundant CO/N2O channel but now carries CH4/H2O instead), or
    None for "No Figure". "Other" plots whatever column the catch-all combo
    box is set to (its units are unknown, so the column name doubles as the
    ylabel); None if no column is selected yet. Ozone has no detector of its
    own, so Detector Pressure/T_gas have nothing to route to and fall back
    to None (same as "No Figure") rather than raising.
    """
    detector = GASES[gas]["detector"]
    if selection == "Detector Pressure":
        if detector is None:
            return None
        col = f"{detector}_P_mbars"
        return col, f"{col} (mbar)"
    if selection == "T_gas":
        if detector is None:
            return None
        col = f"{detector}_T_gas"
        return col, f"{col} (°C)"
    if selection == "Other":
        if other_column is None:
            return None
        return other_column, other_column
    return None

# TWO app-level files beside the script, split by whether the contents are
# worth sharing:
#
# - `ucatsb_gui_config.yaml` holds the ICARTT header metadata and is TRACKED
#   IN GIT. PI, affiliation, project and stipulations describe the campaign,
#   so every machine analysing a SABRE flight wants the same values and a
#   fresh clone should arrive with them already filled in.
# - `.ucatsb_gui_state.yaml` holds the recent-files list and is not. It is
#   absolute paths from one machine, of no use to anyone else, and hidden
#   because the app maintains it -- there is nothing in it to hand-edit.
#
# The split is what makes tracking the first file practical: the recent list
# is rewritten on **every dataset load**, so while the two shared a file that
# file was permanently modified in the working tree. The shared one is now
# written only by the Export tab's "Save defaults" button.
DEFAULT_CONFIG_PATH = Path(__file__).parent / "ucatsb_gui_config.yaml"
DEFAULT_STATE_PATH = Path(__file__).parent / ".ucatsb_gui_state.yaml"
RECENT_FILES_MAX = 10

# What the Correlations tab can color points by: {key: (label, column,
# colorbar label)}. `column` of None means the time axis, which is not a
# plottable column and needs its own numeric conversion and tick format.
# A key whose column is missing from the loaded CSV is simply not offered --
# the schema differs between flights.
CORR_COLOR_BY = {
    "time": ("Time", None, "Time (UTC-ish)"),
    "ozone": ("Ozone", "oz_o3best", "Ozone (ppb)"),
    "oz_p": ("Pressure (oz_p)", "oz_p", "oz_p (mbar)"),
}
# "turbo" rather than "jet" or "rainbow": it is the same rainbow ordering,
# rebuilt with monotonic luminance, so it does not invent bands of false
# structure where the older maps go light-dark-light. A colorbar is always
# drawn with it -- a continuous color encoding with no scale is unreadable.
CORR_COLORMAP = "turbo"

DEFAULT_GAS_SETTINGS = {
    "warmup_min": 30,
    "end_flight_min": 0,
    "require_pumps": False,
    "pressure_tol_mbar": 0.4,
    "pressure_correct": False,
    "pressure_smooth_s": 0,
    "temperature_correct": False,
    "flag_air_s": 0,
    "cal1_window_s": [-15, -1],
    "cal2_window_s": [-15, -1],
    "drift_model": CAL_DRIFT_MODELS[0],
    "drift_smooth_events": CAL_DEFAULT_SMOOTH_EVENTS,
    # 0 = "use the slope the constant model gives" (= 1/span_gain). Read only
    # under the "fixed slope" model, so the default is inert for every config
    # that predates it.
    "fixed_slope": 0.0,
}


def sanitize_config_variant(text: str) -> str:
    """Reduce free text to the part of a filename it is allowed to be.

    The variant is typed by the user but lands straight in a path, so anything
    that could steer it out of the CSV's directory -- separators, `..`, a
    leading dot -- is removed rather than escaped, and everything else
    collapses to underscores so the result is one word. A trailing `conf` is
    dropped too: the suffix is already added below, and a user copying the
    shape of an existing name would otherwise get `..._test_conf_conf.yaml`.
    """
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("._-")
    cleaned = re.sub(r"(?:^|_)conf$", "", cleaned)
    return cleaned.strip("._-")


def flight_config_path(csv_path: Path, variant: str = "") -> Path:
    """Where a dataset's own settings live: <dataset>_conf.yaml, beside the
    CSV, or <dataset>_<variant>_conf.yaml for a second configuration of the
    same flight. Per-flight rather than global because the right warm-up,
    pressure tolerance, cal windows and -- above all -- cal tanks are
    properties of the flight, and re-deriving them every time a file is
    reopened loses work.

    The dataset stem is not decoration: `_config_candidates` finds a flight's
    configs by it, so a name without it is a config the next open cannot see.
    That is why the name is composed here from a variant rather than typed
    whole into a file dialog.
    """
    stem = Path(csv_path).stem
    variant = sanitize_config_variant(variant)
    middle = f"_{variant}" if variant else ""
    return Path(csv_path).with_name(f"{stem}{middle}_conf.yaml")


def config_variant_name(csv_path: Path, config_path: Path) -> str:
    """The variant `flight_config_path` would need to reproduce `config_path`,
    or "" if the name does not follow the scheme (an older config named freely,
    or one opened from elsewhere) -- in which case offering the default name is
    the right fallback."""
    if config_path is None or csv_path is None:
        return ""
    name = Path(config_path).name
    stem = Path(csv_path).stem
    if not name.startswith(stem) or not name.endswith("_conf.yaml"):
        return ""
    return sanitize_config_variant(name[len(stem):-len("_conf.yaml")])


def _read_yaml(path: Path) -> dict:
    if path is None or not path.exists():
        return {}
    try:
        return yaml.safe_load(path.read_text()) or {}
    except (yaml.YAMLError, OSError) as e:
        print(f"Warning: could not read {path}: {e}")
        return {}


def load_config(path: Path) -> dict:
    """Load per-gas control settings, filling anything missing from
    DEFAULT_GAS_SETTINGS. `path=None` (or a file that does not exist) gives
    the defaults outright, which is what a dataset with no saved config gets.

    There is deliberately no template/seeding path: settings are never
    written without an explicit Save, so there is no continuously-updated
    app-level template to inherit from, and inheriting one flight's tuning
    into another silently is what the per-dataset config files exist to
    prevent.

    Gas blocks only. Gases with has_masking=False (Ozone, H2O) never get an
    entry -- they don't use warm-up/pressure-tol/cal-window settings at all.
    The tank selection shares the same file but is kept out of this dict
    deliberately -- see load_cal_selection.
    """
    config = {
        gas: copy.deepcopy(DEFAULT_GAS_SETTINGS)
        for gas, info in GASES.items() if info.get("has_masking", True)
    }
    for gas, settings in config.items():
        for key, value in DEFAULT_GAS_SETTINGS.items():
            settings.setdefault(key, copy.deepcopy(value))
    loaded = _read_yaml(path)
    for gas in config:
        if isinstance(loaded.get(gas), dict):
            config[gas].update(loaded[gas])
    return config


def load_recent_files(path: Path, legacy_path: Path = None) -> list:
    """The recently-opened datasets from the app-level state file, newest first.

    App-level and never per-flight, the mirror image of `cals:` (which is
    flight-only): "files I have been working on" is a property of the person,
    not of any one dataset, and writing it into a flight's conf would put a
    list of unrelated paths beside that flight's settings.

    `legacy_path` is the shared config, where this list used to live before
    the two were split. It is read only when the state file has nothing, and
    never written back there -- the next save lands in the state file and the
    stale block in the shared config is simply ignored from then on. Without
    this a user upgrading loses a list they never asked to lose.
    """
    loaded = _read_yaml(path).get("recent_files")
    if not isinstance(loaded, list) and legacy_path is not None:
        loaded = _read_yaml(legacy_path).get("recent_files")
    if not isinstance(loaded, list):
        return []
    return list(dict.fromkeys(p for p in loaded if isinstance(p, str)))[:RECENT_FILES_MAX]


def load_icartt_meta(path: Path) -> dict:
    """The default ICARTT header metadata from the app-level config, filled
    out from DEFAULT_ICARTT_META.

    App-level rather than per-dataset, and deliberately so for now: the PI,
    affiliation, project and stipulations are the same on every flight of a
    campaign, so making them per-flight would mean retyping them for each
    file. The fields that genuinely vary per flight (the mission, the
    location) are still editable before each export -- they are just not
    saved with the flight. Only unknown keys are dropped, so a config written
    by a later version that adds a field is not corrupted by this one.
    """
    loaded = _read_yaml(path).get("icartt")
    meta = copy.deepcopy(DEFAULT_ICARTT_META)
    if isinstance(loaded, dict):
        # Stripped on the way in so a hand-edited YAML with a trailing space
        # cannot make the form read dirty the moment it is loaded --
        # _icartt_is_dirty compares against stripped widget text.
        meta.update({k: ("" if v is None else str(v).strip())
                     for k, v in loaded.items() if k in DEFAULT_ICARTT_META})
    return meta


def load_cal_selection(path: Path, default: dict) -> dict:
    """The flight's `cals: {cal0: ..., cal1: ...}` choice from its
    <dataset>_conf.yaml, falling back to `default` (cals.yaml's own block).

    Kept out of load_config's dict, and out of the app-level config file, on
    purpose: cals.yaml describes the tanks plumbed in *now*, so it is the only
    safe default for a flight nobody has assigned tanks to yet. Seeding it
    from "whatever the last flight used" instead would silently apply one
    flight's tanks to another, which is exactly the error the Cal Tanks tab
    exists to prevent.
    """
    selection = dict(default)
    loaded = _read_yaml(path).get("cals")
    if isinstance(loaded, dict):
        selection.update({k: v for k, v in loaded.items()
                          if k in ("cal0", "cal1") and isinstance(v, str)})
    return selection


def load_flagged(path: Path, raw_rows: int = None):
    """The flight's manually flagged row ranges, as {gas: [(lo, hi), ...]}.

    Returns ({}, None) when the file has no `flagged:` block, and
    (ranges, complaint) otherwise -- `complaint` being a displayable sentence
    when a gas's stored `rows:` disagrees with the CSV actually loaded.

    Row numbers are the RAW file's, which is the one thing about this format
    worth being careful with. They are exact and permanent -- unlike a stored
    rectangle, a flagged row cannot quietly become a different row when the
    drift model or the cal tanks change -- but they are only meaningful for the
    file they were drawn on. `rows:` is the tripwire for that: a regenerated
    CSV with a different length would otherwise shift every flag silently, and
    silently wrong is the one outcome this feature cannot have. The flags are
    still applied when it trips, because the user is better placed to judge
    than we are; they are just told.
    """
    loaded = _read_yaml(path).get("flagged")
    if not isinstance(loaded, dict):
        return {}, None
    flagged, mismatches = {}, []
    for gas, block in loaded.items():
        if gas not in GASES or not isinstance(block, dict):
            continue
        raw = block.get("ranges")
        if not isinstance(raw, list):
            continue
        pairs = [(r[0], r[1]) for r in raw
                 if isinstance(r, (list, tuple)) and len(r) == 2
                 and all(isinstance(v, int) for v in r)]
        if not pairs:
            continue
        flagged[gas] = merge_ranges(pairs)
        stored_rows = block.get("rows")
        if raw_rows is not None and isinstance(stored_rows, int) and stored_rows != raw_rows:
            mismatches.append(f"{gas} (saved against {stored_rows} rows)")
    complaint = None
    if mismatches:
        complaint = (f"Flagged rows were saved against a different CSV length "
                     f"than the {raw_rows} rows now loaded: "
                     f"{', '.join(mismatches)}. They have been applied as "
                     f"stored — check they still land on the right points.")
    return flagged, complaint


def flagged_to_yaml(flagged: dict, raw_rows: int):
    """The `flagged:` block as plain data, dropping gases with nothing flagged
    so an untouched flight writes no block at all."""
    return {
        gas: {"rows": raw_rows, "ranges": [[lo, hi] for lo, hi in ranges]}
        for gas, ranges in sorted(flagged.items()) if ranges
    }


class _BlockStyleDumper(yaml.SafeDumper):
    """A SafeDumper that writes multi-line strings as YAML block scalars.

    The default quoted style renders an embedded newline as a blank line
    inside a quoted scalar, which round-trips correctly but is close to
    unreadable — and the ICARTT special-comments block is meant to be
    hand-edited in this file. A Dumper subclass rather than a global
    `add_representer` call so nothing else's YAML output changes.
    """


def _represent_str(dumper, data):
    if "\n" in data:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    return dumper.represent_scalar("tag:yaml.org,2002:str", data)


_BlockStyleDumper.add_representer(str, _represent_str)


def save_config(path: Path, config: dict, cal_selection: dict = None,
                flagged: dict = None,
                recent_files: list = None, icartt_meta: dict = None):
    """Write the per-gas blocks, plus whichever of the non-gas blocks belongs
    in this file: the tank selection and the manually flagged rows for a
    flight's own conf, the ICARTT header metadata for the shared app config,
    the recent-file list for the app state file. **One file's blocks per
    call** -- see load_cal_selection, load_flagged, load_icartt_meta and
    load_recent_files for why each lives where it does.

    This writes a *fresh* document, so an omitted block is a deletion, the
    same trap the per-gas settings have with _controls_to_settings(). The
    flight config now carries *two* non-gas blocks, so that trap is live again
    on this path: a save that passes `cal_selection=` but forgets `flagged=`
    silently discards every flagged point. on_save_clicked passes both.
    """
    doc = dict(config)
    if flagged:
        doc = {"flagged": dict(flagged), **doc}
    if cal_selection:
        doc = {"cals": dict(cal_selection), **doc}
    if icartt_meta:
        doc = {"icartt": dict(icartt_meta), **doc}
    if recent_files:
        doc = {"recent_files": list(recent_files), **doc}
    path.write_text(yaml.dump(doc, Dumper=_BlockStyleDumper, sort_keys=False,
                              default_flow_style=False, width=100))


class _NavToolbar(NavigationToolbar):
    """The stock toolbar, with a Home that can be pointed somewhere else.

    Home normally returns to the first entry of the nav stack, and the only
    way to change that entry through the public API is to clear the stack and
    push a new base -- which throws away every zoom and pan the user has done,
    so Back and Forward stop working. Overriding `home` instead leaves the
    stack completely alone: the history is theirs, and only where Home lands
    is ours to redirect.

    The override is a plain (axes, xlim, ylim) triple rather than a callback,
    and it is checked against the figure's live axes because the panes rebuild
    their Figure from scratch on every draw -- an override captured before a
    redraw refers to a destroyed Axes and must not be applied to the new one.
    """

    def __init__(self, canvas, parent=None):
        super().__init__(canvas, parent)
        self._home_override = None

    def set_home_override(self, ax=None, xlim=None, ylim=None):
        """Point Home at `xlim`/`ylim` on `ax`; `ax=None` restores the
        stock behaviour. Nothing is drawn and no limits are touched -- the
        view on screen is not this method's business."""
        self._home_override = None if ax is None else (ax, xlim, ylim)

    def home(self, *args):
        override = self._home_override
        if override is not None and override[0] in self.canvas.figure.axes:
            ax, xlim, ylim = override
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            # Pushed so Back still returns to wherever they were standing.
            self.push_current()
            self.canvas.draw_idle()
            return
        super().home(*args)


class PlotPane(QWidget):
    """One matplotlib Figure with its own toolbar, as a tab page.

    Both views rebuild their Figure from scratch on every draw (the panel
    count varies), so each needs its own toolbar nav-stack reset -- hence a
    widget rather than a bare canvas.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.figure = Figure(constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = _NavToolbar(self.canvas, self)

        # Appended to the stock toolbar rather than declared through
        # NavigationToolbar.toolitems: toolitems entries have to name a method
        # on the toolbar class, which would mean subclassing it just to reach
        # back into the pane.
        self.toolbar.addSeparator()
        self.stats_action = QAction("Stats", self.toolbar)
        self.stats_action.setCheckable(True)
        self.stats_action.setToolTip(
            "Drag a box over the data to get n / mean / std of the points\n"
            "inside it. Cal periods, masked data and post-cal flush are\n"
            "excluded from the statistics and reported separately."
        )
        self.stats_action.toggled.connect(self._on_stats_toggled)
        self.toolbar.addAction(self.stats_action)

        # A second mode on the SAME selector rather than a second selector.
        # attach_stats_selectors() rebuilds self.selectors wholesale on every
        # draw and disconnects the old event handlers; a parallel list would
        # have to duplicate all of that and would leak canvas connections if
        # it ever got out of step. One selector, one mode flag.
        self.flag_action = QAction("Flag", self.toolbar)
        self.flag_action.setCheckable(True)
        self.flag_action.setToolTip(
            "Drag a box to flag errant points, removing them from the\n"
            "calibrated or filtered record. Right-drag to unflag.\n\n"
            "Flagging matches the box against the RAW (blue) trace, so a\n"
            "flag means the same points after a drift-model or cal-tank\n"
            "change. Unflagging ignores the box height and clears the whole\n"
            "time span — a flagged spike is often off the top of the axes.\n\n"
            "Flags are saved with the flight when you press Save."
        )
        self.flag_action.toggled.connect(self._on_flag_toggled)
        self.toolbar.addAction(self.flag_action)

        # After addAction, not before: widgetForAction only has a button to
        # return once the action is on the toolbar.
        for action in (self.stats_action, self.flag_action):
            self._style_toggle(action)

        # A view toggle, not a setting: it overlays this repo's calibrated
        # series on the figure and changes nothing else, so it belongs beside
        # the other things that decide what the figure shows rather than in the
        # per-gas settings panel, where it was taking a row and reading like
        # something that gets saved. Nothing connects it here -- the owner
        # does, and hides it on the panes that have no raw trace to overlay.
        self.toolbar.addSeparator()
        self.calibrated_action = QAction("Calibrated", self.toolbar)
        self.calibrated_action.setCheckable(True)
        self.calibrated_action.setToolTip(
            "Overlay the calibrated series (red) on the raw trace (blue).\n\n"
            "The raw trace stays exactly as it is -- the two are told apart\n"
            "by colour, not by which is faded. Session-only and off at\n"
            "startup, so the app never opens showing calibrated data\n"
            "without being asked."
        )
        self.toolbar.addAction(self.calibrated_action)

        readout = QHBoxLayout()
        self.stats_combo = QComboBox()
        self.stats_combo.setMinimumWidth(190)
        self.stats_combo.setToolTip("Which plotted trace the box statistics apply to")
        self.stats_combo.currentIndexChanged.connect(self._on_trace_changed)
        self.stats_label = QLabel("")
        # Selectable so the numbers can be dragged out even without the button;
        # wrapped rather than elided, since the line is long enough to lose its
        # tail on a narrow window and the tail is where the caveats are.
        self.stats_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.stats_label.setWordWrap(True)
        self.stats_copy_button = QPushButton("Copy")
        self.stats_copy_button.setMaximumWidth(60)
        self.stats_copy_button.clicked.connect(self._copy_stats)
        readout.addWidget(self.stats_combo)
        readout.addWidget(self.stats_label, 1)
        readout.addWidget(self.stats_copy_button)
        self._set_readout_visible(False)

        layout.addWidget(self.toolbar)
        layout.addLayout(readout)
        # Stretch 1 on the canvas and 0 on everything else: the readout row
        # gets its sizeHint and the plot takes all the rest. Without it the
        # row and the canvas are both Preferred with no stretch, so QVBoxLayout
        # splits the spare space EQUALLY and the figure loses half its height
        # the moment the readout appears. That stayed hidden while the row
        # always contained the trace combo -- a QComboBox is Fixed vertically,
        # which caps the whole QHBoxLayout -- and only showed up once the flag
        # tool started showing the label on its own.
        layout.addWidget(self.canvas, 1)

        # on_box/on_flag_box are set by the owner; the selectors themselves
        # are rebuilt on every draw (see attach_stats_selectors).
        self.on_box = None
        self.on_flag_box = None
        self.selectors = []
        # Which tool the shared selector is currently acting as: None, "stats"
        # or "flag". The two toolbar actions are mutually exclusive because
        # both want the same left-drag.
        self.selector_mode = None
        # (axes index, extents) of the live stats box -- only one exists at a
        # time, in whichever panel it was drawn in. Never set in flag mode: a
        # flag is an action, not a standing selection.
        self._box = None
        self._loading_traces = False

    def _style_toggle(self, action):
        """Give a checkable toolbar action a green ON state.

        Styled per button rather than through one `QToolBar QToolButton` rule
        on the toolbar: a stylesheet makes Qt draw the widget itself instead of
        asking the platform style, and a toolbar-wide rule would take the
        native look off the stock nav buttons (home/pan/zoom/save) as well --
        for a colour change these two never asked for.

        Only `:checked` is given a background, so the off state stays whatever
        the platform draws. The padding is set in the base rule rather than
        inside `:checked`, or the button would resize as it is toggled.
        """
        button = self.toolbar.widgetForAction(action)
        if button is None:
            return
        button.setStyleSheet(
            f"QToolButton {{ padding: 3px 8px; border-radius: 4px; }}"
            f"QToolButton:checked {{ background-color: {TOOL_ON_BG};"
            f" color: {TOOL_ON_FG}; }}"
            f"QToolButton:checked:hover {{ background-color: {TOOL_ON_BG_HOVER}; }}"
        )

    def reset_nav(self):
        """Point the toolbar's Home at the newly-built full-scale view; its
        nav stack otherwise still references the just-destroyed Axes.

        Also drops any Home override: a fresh draw has just established what
        full scale means, and an override from before it described a figure
        that no longer exists.
        """
        self.toolbar.set_home_override(None)
        self.toolbar.update()
        self.toolbar.push_current()

    def set_home_view(self, ax, xlim, ylim):
        """Point Home at `xlim`/`ylim` without touching the current view.

        Needed because hiding the flagged points must not replot -- the user
        is typically zoomed in on the very points being hidden -- yet Home
        should then frame what is left rather than a range set by markers that
        are no longer drawn.

        This sets nothing and draws nothing; it only records where Home goes
        (see _NavToolbar). An earlier version applied the range, pushed it as
        a new nav-stack base and put the old view back, which worked but wiped
        the user's zoom/pan history on every toggle and moved the axis limits
        twice per click for no visible reason.
        """
        self.toolbar.set_home_override(ax, xlim, ylim)

    def clear_home_view(self):
        """Give Home back to the nav stack's own base."""
        self.toolbar.set_home_override(None)

    def _set_readout_visible(self, visible, stats_widgets=True):
        """Show the readout row. The label is shared by both tools -- the flag
        tool reports what it just did there -- but the trace combo and Copy
        button are meaningless to it, so they only appear for Stats."""
        self.stats_combo.setVisible(visible and stats_widgets)
        self.stats_label.setVisible(visible)
        self.stats_copy_button.setVisible(visible and stats_widgets)

    def _copy_stats(self):
        QApplication.clipboard().setText(self.stats_label.text())

    def _release_widgetlock(self):
        """Free the canvas widgetlock that pan/zoom hold.

        `_SelectorWidget.ignore()` drops every event while it is held, so a
        selector switched on under an active pan or zoom would look dead
        rather than merely inactive. Re-invoking the toolbar's own toggle is
        what actually releases it. (The reverse needs no handling: clicking
        pan afterwards just makes the selector inert until pan goes off.)
        """
        mode = str(self.toolbar.mode)
        if "pan" in mode:
            self.toolbar.pan()
        elif "zoom" in mode:
            self.toolbar.zoom()

    def _on_stats_toggled(self, checked):
        if checked:
            self._release_widgetlock()
            if self.flag_action.isChecked():
                self.flag_action.setChecked(False)   # re-enters here via flag
        self._set_selector_mode("stats" if checked else None)
        self._set_readout_visible(checked and bool(self.stats_label.text()))
        self.canvas.draw_idle()

    def _on_flag_toggled(self, checked):
        if checked:
            self._release_widgetlock()
            if self.stats_action.isChecked():
                self.stats_action.setChecked(False)
        self._set_selector_mode("flag" if checked else None)

    def _set_selector_mode(self, mode):
        """Switch the shared selector between tools.

        Rebuilt rather than reconfigured: the two modes differ in `props` and
        in `interactive`, neither of which RectangleSelector exposes reliably
        after construction across matplotlib versions, and the rebuild path
        already exists and costs nothing. The stats box's extents survive
        because attach_stats_selectors carries them across.
        """
        if mode == self.selector_mode:
            return
        # A stats readout left standing under the flag tool would describe a
        # box that is no longer on screen, and vice versa.
        self.stats_label.setText("")
        self._set_readout_visible(False)
        self.selector_mode = mode
        axes = [sel.ax for sel in self.selectors]
        if axes:
            self.attach_stats_selectors(axes)
        self.canvas.draw_idle()

    def set_stats_text(self, text):
        self.stats_label.setText(text)
        self._set_readout_visible(bool(text) and self.selector_mode is not None,
                                  stats_widgets=self.selector_mode == "stats")

    def current_trace_key(self):
        return self.stats_combo.currentData()

    def set_stats_traces(self, traces):
        """Populate the trace selector with [(key, label), ...] for whatever is
        currently plotted, keeping the previous choice when it still exists."""
        previous = self.current_trace_key()
        self._loading_traces = True
        self.stats_combo.clear()
        for key, label in traces:
            self.stats_combo.addItem(label, key)
        if previous is not None:
            index = self.stats_combo.findData(previous)
            if index >= 0:
                self.stats_combo.setCurrentIndex(index)
        self._loading_traces = False

    def _on_trace_changed(self, _index):
        """Recompute against the box already on screen -- switching trace
        shouldn't require redrawing the box."""
        if self._loading_traces or self._box is None or self.on_box is None:
            return
        index, extents = self._box
        self.on_box(self.selectors[index].ax, *extents)

    def attach_stats_selectors(self, axes):
        """(Re)bind one box selector per panel in `axes`.

        Must be called after every draw: the panes rebuild their Figure from
        scratch, so a selector from the previous draw holds a dead Axes and
        silently stops responding. The live box is carried across so a
        selection survives a masking tweak instead of vanishing.
        """
        for sel in self.selectors:
            try:
                sel.set_active(False)
                sel.disconnect_events()
            except Exception:
                pass
        self.selectors = []

        # Identity is stale after a redraw (the old Axes are destroyed), so the
        # box is restored by position in the axes list rather than by object.
        box_index, box_extents = (None, None)
        if self._box is not None:
            box_index, box_extents = self._box[0], self._box[1]
        self._box = None

        live = [ax for ax in axes if ax is not None]
        flagging = self.selector_mode == "flag"
        active = self.selector_mode is not None
        for i, ax in enumerate(live):
            # A RectangleSelector adds its rectangle (and, when interactive,
            # its corner handles) to the Axes at the ORIGIN, and those artists
            # enlarge ax.dataLim to include (0, 0). On a tracer-tracer scatter
            # that is catastrophic and visible: N2O spans 304-341 ppb, so
            # reaching back to zero turns the whole correlation into a blob in
            # the corner. The timeseries survived it only by accident -- its
            # view limits are settled by other means -- but its dataLim was
            # being corrupted too, which any later autoscale would have
            # inherited. Snapshot and restore, so the limits describe the data
            # and nothing else.
            saved_limits = ax.dataLim.frozen()
            sel = RectangleSelector(
                ax, functools.partial(self._on_select, ax), useblit=True,
                # Flag mode is non-interactive and takes the right button too:
                # a flag is an action that fires and clears, not a standing
                # selection with drag handles, and right-drag is the unflag.
                interactive=not flagging,
                button=[1, 3] if flagging else [1],
                minspanx=3, minspany=3, spancoords="pixels",
                props=dict(facecolor="none",
                           edgecolor=FLAGGED_COLOR if flagging else STATS_BOX_COLOR,
                           linewidth=1.4, linestyle="-" if flagging else "--"),
            )
            sel.set_active(active)
            if not flagging and i == box_index and box_extents is not None:
                try:
                    sel.extents = box_extents
                    self._box = (i, box_extents)
                except Exception:
                    pass
            sel.set_visible(active and not flagging
                            and i == box_index and self._box is not None)
            # Restoring dataLim is not enough on its own: adding the artists
            # already triggered an autoscale off the polluted limits, and
            # nothing recomputes the view afterwards. autoscale_view() is the
            # safe way to redo it -- it only touches an axis whose autoscale
            # is still on, so an explicitly set range (a preserved view, or
            # the ozone y-framing) is left exactly as it was.
            ax.dataLim.set(saved_limits)
            ax.autoscale_view()
            self.selectors.append(sel)

    def _on_select(self, ax, eclick, erelease):
        if self.selector_mode == "flag":
            # Button 3 is the unflag. Read off the press rather than the
            # release so a drag that starts with the right button and wanders
            # still means "unflag".
            if self.on_flag_box is not None:
                x0, x1 = sorted((eclick.xdata, erelease.xdata))
                y0, y1 = sorted((eclick.ydata, erelease.ydata))
                if None not in (x0, x1, y0, y1):
                    self.on_flag_box(ax, x0, x1, y0, y1, eclick.button == 3)
            # Clear the rectangle: the flag is applied, and leaving the box up
            # would suggest a selection that no longer means anything.
            for sel in self.selectors:
                sel.set_visible(False)
            self.canvas.draw_idle()
            return
        for sel in self.selectors:
            if sel.ax is ax:
                self._box = (self.selectors.index(sel), sel.extents)
            else:
                # Only one box at a time: dragging in the other panel clears
                # the previous one, so the readout always refers to a box that
                # is actually on screen.
                sel.set_visible(False)
        self.canvas.draw_idle()
        if self.on_box is not None and self._box is not None:
            self.on_box(ax, *self._box[1])


class UcatsbGui(QMainWindow):
    def __init__(self, csv_path: Path = None, config_path: Path = DEFAULT_CONFIG_PATH,
                 state_path: Path = DEFAULT_STATE_PATH):
        super().__init__()
        self.setWindowTitle("UCATS-B Viewer")
        self.resize(1300, 750)

        self.csv_path = None
        self.df = None
        # Set by load_csv; see there for why the Export tab needs them.
        self.raw_df = None
        self.presync_dropped = 0
        self.available_gases = {}
        self.other_columns = []

        # THREE config files, with different jobs, and none of them holds
        # analysis settings for another's business:
        #   default_config_path  shared, tracked  -- ICARTT header metadata
        #   state_path           local, hidden    -- the recent-files list
        #   config_path          per dataset      -- this flight's settings
        # config_path is the one Save writes to (and the name its dialog
        # offers); config_loaded_from is the file actually opened, or None
        # when the dataset started from defaults.
        self.default_config_path = config_path
        self.state_path = state_path
        self.config_path = None
        self.config_loaded_from = None
        self.config = load_config(None)
        # legacy_path picks up a list left in the shared config by a version
        # before the two files were split; see load_recent_files.
        self.recent_files = load_recent_files(state_path, legacy_path=config_path)
        # ICARTT header metadata: app-level, not per-dataset (see
        # load_icartt_meta), with its own Save button on the Export tab and
        # its own saved-state snapshot -- the per-dataset Save writes a
        # different file entirely, so one dirty state cannot serve both.
        self.icartt_meta = load_icartt_meta(config_path)
        self._saved_icartt_meta = copy.deepcopy(self.icartt_meta)
        # None until a dataset is open; _is_dirty() stays False before that.
        self._saved_state = None
        self.current_gas = None
        self.aux_selection = "No Figure"
        self.other_column = None
        self.right_axis_column = None
        # Session-only view state, like the calibrated overlay: on at startup,
        # never saved. `_notes_artist` is the Text the toggle reaches for, and
        # is rebuilt (or set to None) by every redraw -- a stale one would hold
        # a destroyed Axes, the same trap the stats selector has.
        self.show_plot_notes = True
        self._notes_artist = None
        self._loading = False
        self._initializing = True
        # Per-gas caches: the Correlations tab needs two gases' calibrations
        # at once, so these are dicts keyed by gas rather than one slot for
        # whichever gas is on display.
        self._analysis = {}
        self._calibration = {}
        self._uncertainty = {}
        # Rebuilt by every redraw(); initialised here so on_stats_box is safe
        # to reach before the first draw.
        self._stats_traces = {}
        self.drift_model = DEFAULT_GAS_SETTINGS["drift_model"]
        self.drift_smooth_events = DEFAULT_GAS_SETTINGS["drift_smooth_events"]
        self.fixed_slope = DEFAULT_GAS_SETTINGS["fixed_slope"]
        self.show_calibrated = False
        # "export" joins the plot panes so the tab's readiness summary is
        # invalidated by refresh() like everything else -- the masking
        # controls stay visible beside it, so a spin-box nudge has to be able
        # to change what it says about which gases will export.
        self._dirty = {"main": True, "cal": True, "corr": True, "export": True}
        self._preserve = {"main": False, "cal": False, "corr": False, "export": False}
        self.cal_roster = load_cal_roster(CALS_YAML_PATH)
        # cals.yaml's own pairing is the default only; a flight's conf file
        # overrides it, and the Cal Tanks tab edits it.
        self.default_cal_selection = load_cal_assignment(CALS_YAML_PATH)
        self.cal_selection = dict(self.default_cal_selection)
        self.cal_bottles = select_cal_bottles(self.cal_roster,
                                              self.cal_selection.values())
        # Manually flagged rows, {gas: [(lo, hi), ...]} in RAW file row
        # numbers. Per flight and per gas; empty is the no-op, so a flight
        # nobody has flagged behaves exactly as it did before the feature.
        self.flagged = {}
        # Session-only undo stack of previous `flagged` dicts. Deliberately
        # not persisted: it is a property of this editing session, and a
        # config that could undo its own contents would be a strange object.
        self._flag_undo = []
        self.ax = None
        self.ax_aux = None
        self.ax_aux2 = None
        self._had_aux_panel = None
        self._last_aux_key = None
        self._last_right_axis_key = None
        self._cal_ax = None
        self._last_cal_key = None
        # Correlations tab state. Session-only, like the calibrated overlay:
        # which pair of tracers you are looking at is a question you are
        # asking right now, not a property of the flight worth persisting.
        self.corr_x_gas = None
        self.corr_y_gas = None
        # Set by redraw_corr; what the flag tool resolves a dragged box
        # against. None until the tab has drawn something flaggable.
        self._corr_plotted = None
        self._corr_flag_scatter = None
        # Session-only, like the calibrated overlay: which markers you have
        # temporarily taken off the plot is not a property of the flight.
        self.corr_hide_flagged = False
        self.corr_marker_size = 4
        self.corr_error_bars = False
        self.corr_show_cals = False
        # None = single-color points; otherwise a key into CORR_COLOR_BY.
        self.corr_color_by = None
        # Off by default: a straight line through a tracer-tracer plot with
        # real structure in it describes almost none of that structure.
        self.corr_fit = False
        self.corr_tooltip_popup = None
        self._corr_ax = None
        self._last_corr_key = None

        # Before _build_controls: the Load Data button attaches recent_menu,
        # and both menus share this one QAction so "Open…" is defined once.
        self.open_action = QAction("Open…", self)
        self.open_action.setShortcut("Ctrl+O")   # Qt maps this to ⌘O on macOS
        self.open_action.triggered.connect(self.on_load_data_clicked)
        self.recent_menu = QMenu(self)
        self.load_config_action = QAction("Load Configuration…", self)
        self.load_config_action.triggered.connect(self.on_load_config_clicked)
        self.save_action = QAction("Save Configuration…", self)
        self.save_action.setShortcut("Ctrl+S")
        self.save_action.triggered.connect(self.on_save_clicked)
        self._build_menu_bar()

        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)
        # Qt's default 9 px all round is 18 px of height the controls panel
        # would rather have -- enough, at the default window size, to be the
        # difference between the panel fitting and a scrollbar appearing for
        # the sake of a few pixels. The tabs and the group boxes provide all
        # the visual separation this edge needs.
        layout.setContentsMargins(4, 4, 4, 4)

        # The control panel is a stack, not one panel: everything in the
        # main panel is per-gas and the Correlations tab is inherently about
        # two gases at once, so its controls replace them rather than sitting
        # alongside and contradicting them. Tabs that share the per-gas
        # controls (Timeseries, Calibration, Cal Tanks) still share one panel.
        self.controls_stack = QStackedWidget()
        self.controls_stack.addWidget(self._build_controls())
        self.controls_stack.addWidget(self._build_corr_controls())

        # The stack goes inside a QScrollArea so the WINDOW's minimum height
        # stops being "however tall the tallest control panel happens to be".
        # Without it the panel's ~1090 px sizeHint set a hard floor of ~1130 px
        # on the window, which does not fit a laptop screen at default scaling
        # -- and every control added since made it worse with no warning. The
        # panel is compact enough that the scrollbar normally never appears;
        # this is the backstop that makes that a layout detail rather than a
        # constraint on the whole app. Same pattern as the Export tab.
        controls_scroll = QScrollArea()
        controls_scroll.setWidget(self.controls_stack)
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setFrameShape(QFrame.NoFrame)
        # Never scroll sideways: the panel is built to a fixed width, so a
        # horizontal bar would only ever mean something has been laid out
        # wrongly, and stealing height for it would be the wrong response.
        controls_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        # Reserve the vertical scrollbar's width on top of the content width,
        # so the controls get their full CONTROLS_WIDTH whether or not the bar
        # is showing and nothing reflows when it appears.
        sb_extent = self.style().pixelMetric(QStyle.PM_ScrollBarExtent)
        controls_scroll.setFixedWidth(CONTROLS_WIDTH + sb_extent)
        layout.addWidget(controls_scroll, 0)
        layout.addWidget(self._build_tabs(), 1)

        self._initializing = False

        if csv_path is not None:
            self.load_csv(Path(csv_path))

    def load_csv(self, csv_path: Path):
        """Load (or replace) the dataset being viewed, refreshing every
        control that depends on which columns/gases this file has. Raises
        before touching any existing state if the file can't be read or has
        none of the expected gas columns, so a bad pick from the file
        browser leaves whatever was already loaded intact.
        """
        csv_path = Path(csv_path)

        # Detector wiring has changed between flights (e.g. d2 used to carry
        # a redundant CO/N2O channel, now carries CH4/H2O instead), so don't
        # assume every column in REQUIRED_COLUMNS exists in a given file.
        available_cols = set(pd.read_csv(csv_path, nrows=0).columns)
        missing = [c for c in REQUIRED_COLUMNS if c not in available_cols]
        if missing:
            print(f"Note: {csv_path.name} is missing columns {missing}; related features unavailable.")

        available_gases = {
            gas: info for gas, info in GASES.items() if info["value_col"] in available_cols
        }
        if not available_gases:
            raise ValueError(f"{csv_path.name} has none of the expected gas columns: "
                              f"{[g['value_col'] for g in GASES.values()]}")

        # Read every column, not just REQUIRED_COLUMNS, so the "Other" aux
        # trace catch-all has the full remaining roster of CSV fields to
        # offer regardless of what any given flight's schema contains.
        raw = pd.read_csv(csv_path)
        raw["datetime"] = pd.to_datetime(raw["datetime"])
        df = drop_presync_rows(raw)

        # Columns not already exposed via a named control -- offered through
        # the "Other" catch-all combo box. "datetime" is excluded too since
        # it's never a meaningful y-trace.
        other_columns = sorted(c for c in df.columns if c not in NAMED_TRACE_COLUMNS and c != "datetime")

        # Validation above passed -- safe to commit the new dataset now.
        self.csv_path = csv_path
        self.df = df
        # The timestamps and gas columns as read, before drop_presync_rows,
        # plus the number of rows it took off the front. Kept for the Export
        # tab alone: the companion CSV promises to be the same length as the
        # file it complements, and drop_presync_rows resets the index, so the
        # offset is the only thing that can put a derived Series back onto the
        # raw file's row numbers. The gas columns are kept too because the
        # exported raw echo must equal the source column *exactly* -- a
        # pre-sync row has an unreliable clock, not an unreliable reading, so
        # blanking its measurement would make the echo disagree with the file
        # it claims to echo. Narrow slice, not the whole frame: nothing else
        # about those rows is ever wanted.
        self.raw_df = raw[["datetime"] + [info["value_col"]
                                          for info in available_gases.values()]].copy()
        self.presync_dropped = len(raw) - len(df)
        self.available_gases = available_gases
        self.other_columns = other_columns
        self.current_gas = next(iter(available_gases))
        self.aux_selection = "No Figure"
        self.other_column = other_columns[0] if other_columns else None
        self.right_axis_column = None

        self._load_flight_config(csv_path)
        # After the commit above, so a file that failed validation never
        # enters the list. main() loads through here too, which is why a
        # dataset opened from the command line is remembered as well.
        self._remember_recent(csv_path)

        self.setWindowTitle(f"UCATS-B Viewer - {csv_path.name}")
        self._update_file_labels()
        # A flight whose schema predates j_pumps can't offer the filter. The
        # explicit disable survives mask_box's setEnabled(has_masking), since
        # Qt restores a child's own enabled state when its parent comes back.
        self.pumps_check.setEnabled("j_pumps" in df.columns)
        if "j_pumps" not in df.columns:
            self.pumps_check.setChecked(False)

        self._populate_corr_combos()
        self._populate_corr_color_combo()

        was_initializing = self._initializing
        self._initializing = True

        self.gas_combo.blockSignals(True)
        self.gas_combo.clear()
        self.gas_combo.addItems(available_gases.keys())
        self.gas_combo.blockSignals(False)

        self.aux_combo.blockSignals(True)
        self.aux_combo.setCurrentText("No Figure")
        self.aux_combo.blockSignals(False)
        self.aux_selection = "No Figure"
        self.other_combo.setEnabled(False)
        self.right_axis_combo.setEnabled(False)

        self.other_combo.blockSignals(True)
        self.other_combo.clear()
        self.other_combo.addItems(other_columns)
        self.other_combo.blockSignals(False)

        self.right_axis_combo.blockSignals(True)
        self.right_axis_combo.clear()
        self.right_axis_combo.addItem("(none)")
        self.right_axis_combo.addItems(other_columns)
        self.right_axis_combo.blockSignals(False)

        self._select_gas(self.current_gas)
        self._initializing = was_initializing
        if not self._initializing:
            self.refresh()

    def _config_candidates(self, csv_path: Path):
        """Config files that belong to this dataset, default first.

        Any `<dataset stem>*.yaml` beside the CSV: Save names a variant, so
        `..._conf.yaml` and `..._tight_conf.yaml` are both this dataset's
        configs, while another flight's are excluded by the stem.

        Deliberately looser than the names Save now composes -- it also matches
        `..._v2.yaml` and anything else stem-prefixed. Save used to accept a
        free-typed filename, so configs in those shapes exist on disk; matching
        only `*_conf.yaml` would strand them.
        """
        default = flight_config_path(csv_path)
        found = sorted(p for p in csv_path.parent.glob(f"{csv_path.stem}*.yaml")
                       if p.is_file())
        return ([default] if default in found else []) + [p for p in found if p != default]

    def _load_flight_config(self, csv_path: Path):
        """Apply a config for the dataset being opened -- **without writing
        anything**.

        Loading a dataset must leave the disk exactly as it found it: the
        whole point of the explicit Save is that you can open a saved
        analysis, experiment, and quit without disturbing it. So there is no
        write here, and no seeding from an app-level template either -- a
        flight with no config starts from DEFAULT_GAS_SETTINGS.

        One config is applied silently; several bring up the chooser, since
        picking for the user would be picking wrong half the time.
        """
        candidates = self._config_candidates(csv_path)
        chosen = candidates[0] if len(candidates) == 1 else None
        if len(candidates) > 1:
            chosen = self._choose_config_file(csv_path, candidates)

        self._apply_config_file(chosen, fallback_name=flight_config_path(csv_path))

    def _apply_config_file(self, path: Path, fallback_name: Path = None):
        """Load settings + tank selection from `path` (None = shipped
        defaults) and take them as the new saved state."""
        self.config = load_config(path) if path else load_config(None)
        self.cal_selection = (load_cal_selection(path, self.default_cal_selection)
                              if path else dict(self.default_cal_selection))
        raw_rows = len(self.raw_df) if self.raw_df is not None else None
        self.flagged, flag_complaint = (load_flagged(path, raw_rows) if path
                                        else ({}, None))
        self._flag_undo = []
        self._rebuild_cal_bottles()
        # config_path is what Save offers as the default filename. With no
        # config loaded that is the conventional name for this dataset, so
        # the first Save lands where the next open will look.
        self.config_path = path or fallback_name
        self.config_loaded_from = path
        self._snapshot_state()
        self._update_file_labels()
        self._apply_cal_selection_to_controls()
        # After the snapshot: a length mismatch is a warning about what was
        # loaded, not an unsaved change to be offered back.
        self._update_flag_readout()
        if flag_complaint:
            QMessageBox.warning(self, "Flagged points", flag_complaint)

    def _choose_config_file(self, csv_path: Path, candidates):
        """Ask which of several configs to open. Returns a path, or None for
        "start from defaults". Cancel keeps the default (first) entry, since
        the dataset is already committed by the time this runs."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Choose configuration")
        layout = QVBoxLayout(dialog)
        label = QLabel(f"{csv_path.name} has {len(candidates)} saved "
                       f"configurations. Which one should be opened?")
        label.setWordWrap(True)
        layout.addWidget(label)

        group = QButtonGroup(dialog)
        buttons = []
        for i, path in enumerate(candidates):
            rb = QRadioButton(path.name + ("   (default name)" if i == 0 else ""))
            rb.setToolTip(str(path))
            rb.setChecked(i == 0)
            group.addButton(rb)
            layout.addWidget(rb)
            buttons.append((rb, path))
        none_rb = QRadioButton("Start from defaults (open nothing)")
        group.addButton(none_rb)
        layout.addWidget(none_rb)

        box = QDialogButtonBox(QDialogButtonBox.Open | QDialogButtonBox.Cancel)
        box.accepted.connect(dialog.accept)
        box.rejected.connect(dialog.reject)
        layout.addWidget(box)

        if dialog.exec_() != QDialog.Accepted:
            return candidates[0]
        if none_rb.isChecked():
            return None
        return next(path for rb, path in buttons if rb.isChecked())

    def _build_corr_controls(self):
        """Controls for the Correlations tab only (see the stack comment in
        __init__ for why they are not merged into the main panel)."""
        panel = QWidget()
        panel.setFixedWidth(CONTROLS_WIDTH)
        vbox = QVBoxLayout(panel)
        # Same tightening as the per-gas panel, and for the same reason: this
        # is the TALLER of the two stack pages, so it -- not the other one --
        # is what the stack asks the window for.
        vbox.setContentsMargins(6, 6, 6, 6)
        vbox.setSpacing(4)

        self.corr_file_label = QLabel("No file loaded")
        self.corr_file_label.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.corr_file_label)

        axes_box = QGroupBox("Tracers")
        axes_form = QFormLayout(axes_box)
        axes_form.setContentsMargins(6, 6, 6, 6)
        axes_form.setSpacing(4)
        self.corr_x_combo = QComboBox()
        self.corr_x_combo.currentTextChanged.connect(
            functools.partial(self.on_corr_gas_changed, "x"))
        axes_form.addRow("X axis:", self.corr_x_combo)
        self.corr_y_combo = QComboBox()
        self.corr_y_combo.currentTextChanged.connect(
            functools.partial(self.on_corr_gas_changed, "y"))
        axes_form.addRow("Y axis:", self.corr_y_combo)

        swap_button = QPushButton("Swap axes")
        swap_button.clicked.connect(self.on_corr_swap_axes)
        axes_form.addRow(swap_button)
        vbox.addWidget(axes_box)

        style_box = QGroupBox("Points")
        style_form = QFormLayout(style_box)
        style_form.setContentsMargins(6, 6, 6, 6)
        style_form.setSpacing(4)
        self.corr_size_spin = QSpinBox()
        self.corr_size_spin.setRange(1, 20)
        self.corr_size_spin.setValue(self.corr_marker_size)
        self.corr_size_spin.setSuffix(" pt")
        self.corr_size_spin.setToolTip(
            "Marker diameter. A whole flight is tens of thousands of points,\n"
            "so small markers show the structure that large ones fill in."
        )
        self.corr_size_spin.valueChanged.connect(self.on_corr_style_changed)
        style_form.addRow("Marker size:", self.corr_size_spin)

        self.corr_error_check = QCheckBox("Error bars (1σ)")
        self.corr_error_check.setToolTip(
            "1-sigma uncertainty propagated from the calibration itself:\n"
            "the tanks' assigned-value uncertainties from cals.yaml, and how\n"
            "well the drift model reproduces each tank's measured response\n"
            "(leave-one-out + closure). It does NOT include the instrument's\n"
            "own single-sample noise, which the calibration says nothing\n"
            "about. Mostly systematic — it shifts points together rather\n"
            "than scattering them."
        )
        self.corr_error_check.toggled.connect(self.on_corr_style_changed)
        style_form.addRow(self.corr_error_check)

        self.corr_cals_check = QCheckBox("Display cals")
        self.corr_cals_check.setToolTip(
            "Diagnostic overlay only. Draw paired cal-mean points on the\n"
            "correlation figure without adding them to the fit, flagging,\n"
            "tooltip search, or exported air data."
        )
        self.corr_cals_check.toggled.connect(self.on_corr_cals_changed)
        style_form.addRow(self.corr_cals_check)

        self.corr_fit_check = QCheckBox("Linear fit (OLS)")
        self.corr_fit_check.setChecked(self.corr_fit)
        self.corr_fit_check.setToolTip(
            "Least-squares line through the plotted points, with its slope,\n"
            "intercept and r. Off by default: these plots usually have real\n"
            "structure (branches, mixing lines, profiles) that a single\n"
            "straight line describes poorly, and its slope then says more\n"
            "about how the flight was sampled than about the tracers."
        )
        self.corr_fit_check.toggled.connect(self.on_corr_style_changed)
        style_form.addRow(self.corr_fit_check)

        # Checkbox + combo on one spanning row (see the drift-model row: in
        # the two-column form the label gets squeezed to nothing by a
        # stretching field).
        color_row = QHBoxLayout()
        color_row.setSpacing(4)
        self.corr_color_check = QCheckBox("Color by")
        self.corr_color_check.setToolTip(
            "Color each point by a third variable, so the scatter shows\n"
            "where in the flight — or at what pressure — each part of the\n"
            "correlation came from. Points with no value for it are dropped."
        )
        self.corr_color_check.toggled.connect(self.on_corr_color_changed)
        color_row.addWidget(self.corr_color_check)
        self.corr_color_combo = QComboBox()
        self.corr_color_combo.setEnabled(False)
        self.corr_color_combo.currentIndexChanged.connect(self.on_corr_color_changed)
        color_row.addWidget(self.corr_color_combo, 1)
        style_form.addRow(color_row)
        vbox.addWidget(style_box)

        # Flagging from this figure needs one thing the timeseries does not:
        # a point here belongs to TWO gases, and striking out an ozone outlier
        # must not discard the perfectly good N2O measured at the same instant.
        # Hence an explicit target, defaulted to the Y tracer (the one usually
        # read as the dependent variable) and repopulated whenever the axis
        # pickers change.
        self.corr_flag_box = QGroupBox("Flagged Points")
        corr_flag_form = QVBoxLayout(self.corr_flag_box)
        corr_flag_form.setContentsMargins(6, 6, 6, 6)
        corr_flag_form.setSpacing(4)
        target_row = QHBoxLayout()
        target_row.setSpacing(4)
        target_row.addWidget(QLabel("Flag:"))
        self.corr_flag_target = QComboBox()
        self.corr_flag_target.setToolTip(
            "Which tracer a flag drawn here applies to. A point is a pair of\n"
            "measurements, and one of them is usually fine — an ozone spike\n"
            "says nothing about the N2O measured at the same moment."
        )
        self.corr_flag_target.currentIndexChanged.connect(self.on_corr_flag_target_changed)
        target_row.addWidget(self.corr_flag_target, 1)
        corr_flag_form.addLayout(target_row)

        self.corr_flag_label = QLabel("No points flagged")
        self.corr_flag_label.setWordWrap(True)
        self.corr_flag_label.setStyleSheet(f"color: {MUTED_COLOR};")
        corr_flag_form.addWidget(self.corr_flag_label)

        # Hiding is a *view* change, not an edit: it toggles the markers'
        # visibility in place and never redraws, because the user is usually
        # zoomed in on the very points being hidden and a replot would throw
        # that away. Session-only, like the calibrated overlay -- what you are
        # looking at right now is not a property of the flight.
        self.corr_hide_check = QCheckBox("Hide flagged points")
        self.corr_hide_check.setToolTip(
            "Take the struck-out markers off the plot without redrawing, so\n"
            "the current zoom is kept exactly. Home then rescales to the\n"
            "unflagged data. The points stay flagged either way — this only\n"
            "changes what is drawn."
        )
        self.corr_hide_check.toggled.connect(self.on_corr_hide_flagged)
        corr_flag_form.addWidget(self.corr_hide_check)

        corr_flag_buttons = QHBoxLayout()
        self.corr_flag_undo_button = QPushButton("Undo")
        self.corr_flag_undo_button.setToolTip(
            "Step back through this session's flagging, on any tab.")
        self.corr_flag_undo_button.clicked.connect(self.on_flag_undo)
        self.corr_flag_clear_button = QPushButton("Clear")
        self.corr_flag_clear_button.setToolTip(
            "Remove every flag on the tracer named above.")
        self.corr_flag_clear_button.clicked.connect(self.on_corr_flag_clear)
        corr_flag_buttons.addWidget(self.corr_flag_undo_button)
        corr_flag_buttons.addWidget(self.corr_flag_clear_button)
        corr_flag_form.addLayout(corr_flag_buttons)
        vbox.addWidget(self.corr_flag_box)

        self.corr_cal_box = QGroupBox("Calibration Settings")
        corr_cal_form = QFormLayout(self.corr_cal_box)
        corr_cal_form.setContentsMargins(6, 6, 6, 6)
        corr_cal_form.setSpacing(4)
        corr_cal_form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.corr_cal_target = QComboBox()
        self.corr_cal_target.setToolTip(
            "Which correlation-axis tracer these calibration settings edit.\n"
            "They are the same per-gas settings used on the Calibration tab."
        )
        self.corr_cal_target.currentIndexChanged.connect(
            self.on_corr_cal_target_changed)
        corr_cal_form.addRow("Tracer:", self.corr_cal_target)

        self.corr_pressure_correct_check = QCheckBox("Pdelt")
        self.corr_pressure_correct_check.setToolTip(self.pressure_correct_check.toolTip())
        self.corr_pressure_correct_check.toggled.connect(self.on_corr_cal_control_changed)
        self.corr_temperature_correct_check = QCheckBox("Tdelt")
        self.corr_temperature_correct_check.setToolTip(
            self.temperature_correct_check.toolTip())
        self.corr_temperature_correct_check.toggled.connect(
            self.on_corr_cal_control_changed)
        corr_correct_row = QHBoxLayout()
        corr_correct_row.setSpacing(6)
        corr_correct_row.addWidget(self.corr_pressure_correct_check)
        corr_correct_row.addWidget(self.corr_temperature_correct_check)
        corr_correct_row.addStretch(1)
        corr_correct_label = QLabel("Correct:")
        corr_correct_label.setToolTip(self.pressure_correct_check.toolTip())
        corr_cal_form.addRow(corr_correct_label, corr_correct_row)

        self.corr_pressure_smooth_spin = QSpinBox()
        self.corr_pressure_smooth_spin.setRange(0, 300)
        self.corr_pressure_smooth_spin.setSuffix(" s")
        self.corr_pressure_smooth_spin.setMaximumWidth(70)
        self.corr_pressure_smooth_spin.setToolTip(self.pressure_smooth_spin.toolTip())
        self.corr_pressure_smooth_spin.valueChanged.connect(
            self.on_corr_cal_control_changed)
        corr_smooth_p_row = QHBoxLayout()
        corr_smooth_p_row.setSpacing(4)
        corr_smooth_p_row.addWidget(self.corr_pressure_smooth_spin)
        corr_smooth_p_row.addStretch(1)
        corr_smooth_p_label = QLabel("Smooth P:")
        corr_smooth_p_label.setToolTip(self.pressure_smooth_spin.toolTip())
        corr_cal_form.addRow(corr_smooth_p_label, corr_smooth_p_row)

        corr_drift_row = QHBoxLayout()
        corr_drift_row.setSpacing(4)
        corr_drift_label = QLabel("Model:")
        corr_drift_label.setToolTip(self.drift_combo.toolTip())
        corr_drift_row.addWidget(corr_drift_label)
        self.corr_drift_combo = QComboBox()
        self.corr_drift_combo.addItems(CAL_DRIFT_MODELS)
        self.corr_drift_combo.setToolTip(self.drift_combo.toolTip())
        self.corr_drift_combo.currentTextChanged.connect(
            self.on_corr_cal_control_changed)
        corr_drift_row.addWidget(self.corr_drift_combo, 1)

        self.corr_smooth_spin = QSpinBox()
        self.corr_smooth_spin.setRange(2, 21)
        self.corr_smooth_spin.setSuffix(" ev")
        self.corr_smooth_spin.setMaximumWidth(72)
        self.corr_smooth_spin.setToolTip(self.smooth_spin.toolTip())
        self.corr_smooth_spin.valueChanged.connect(self.on_corr_cal_control_changed)
        corr_drift_row.addWidget(self.corr_smooth_spin)

        self.corr_fixed_slope_spin = QDoubleSpinBox()
        self.corr_fixed_slope_spin.setRange(0.0, 10.0)
        self.corr_fixed_slope_spin.setDecimals(4)
        self.corr_fixed_slope_spin.setSingleStep(0.001)
        self.corr_fixed_slope_spin.setSpecialValueText("auto")
        self.corr_fixed_slope_spin.setMaximumWidth(88)
        self.corr_fixed_slope_spin.setToolTip(self.fixed_slope_spin.toolTip())
        self.corr_fixed_slope_spin.valueChanged.connect(
            self.on_corr_cal_control_changed)
        self.corr_fixed_slope_spin.setVisible(False)
        corr_drift_row.addWidget(self.corr_fixed_slope_spin)
        self.corr_fixed_slope_reset_button = QToolButton()
        self.corr_fixed_slope_reset_button.setText("R")
        self.corr_fixed_slope_reset_button.setFixedSize(28, 22)
        self.corr_fixed_slope_reset_button.setStyleSheet(
            "QToolButton { padding: 0; color: #222222; }")
        self.corr_fixed_slope_reset_button.setToolTip(
            "Reset fixed slope to the constant-model value")
        self.corr_fixed_slope_reset_button.clicked.connect(
            self.on_corr_fixed_slope_reset)
        self.corr_fixed_slope_reset_button.setVisible(False)
        corr_drift_row.addWidget(self.corr_fixed_slope_reset_button)
        corr_cal_form.addRow(corr_drift_row)

        vbox.addWidget(self.corr_cal_box)

        self.corr_note = QLabel(
            "Calibrated where available.\n"
            "Masks, cals, and flush apply per gas."
        )
        self.corr_note.setToolTip(
            "Correlation uses calibrated values where available.\n"
            "Masks, cal periods, and flush filters are applied per gas."
        )
        self.corr_note.setWordWrap(True)
        self.corr_note.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.corr_note)

        self.corr_stats_label = QLabel()
        self.corr_stats_label.setWordWrap(True)
        self.corr_stats_label.setTextFormat(Qt.PlainText)
        self.corr_stats_label.setStyleSheet("font-family: Menlo, monospace;")
        vbox.addWidget(self.corr_stats_label)

        vbox.addStretch(1)
        return panel

    def _populate_corr_color_combo(self):
        """Offer the z-axis encodings this file can actually supply. Time
        always works; file columns such as Ozone and `oz_p` are offered only
        when the loaded CSV has them."""
        available = [(key, spec) for key, spec in CORR_COLOR_BY.items()
                     if spec[1] is None or spec[1] in self.df.columns]
        combo = self.corr_color_combo
        combo.blockSignals(True)
        combo.clear()
        for key, (label, _, _) in available:
            combo.addItem(label, key)
        combo.blockSignals(False)
        # A coloring whose column this file lacks cannot survive the reload.
        if self.corr_color_by not in [key for key, _ in available]:
            self.corr_color_by = None
            self.corr_color_check.setChecked(False)

    def _populate_corr_combos(self):
        """Offer every gas in this file, Ozone included.

        Ozone has no cal bottles, so there is nothing to calibrate it against
        and it goes on the axis as the instrument's own product (`oz_o3best`).
        That is worth having anyway: O3 against N2O or CO2 is a standard
        tracer-tracer pairing, and it shares the CSV's timestamps, so the rows
        line up with no resampling. The axis label says which it is.
        """
        gases = list(self.available_gases)
        self.corr_x_gas = gases[0] if gases else None
        # Default to a genuine pair rather than a gas against itself, which
        # would draw a diagonal line and look like a bug.
        self.corr_y_gas = gases[1] if len(gases) > 1 else self.corr_x_gas
        for combo, current in ((self.corr_x_combo, self.corr_x_gas),
                               (self.corr_y_combo, self.corr_y_gas)):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(gases)
            if current:
                combo.setCurrentText(current)
            combo.blockSignals(False)
        self._populate_corr_flag_target()
        self._populate_corr_cal_target()

    def _rebuild_cal_bottles(self):
        """The two tanks matching is allowed to consider, from the current
        selection. Called wherever cal_selection changes -- cal_bottles is
        derived state and must never be edited on its own."""
        self.cal_bottles = select_cal_bottles(self.cal_roster,
                                              self.cal_selection.values())

    def on_load_data_clicked(self):
        if not self._confirm_discard("loading another dataset"):
            return
        start_dir = str(self.csv_path.parent) if self.csv_path else str(Path.cwd())
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load UCATS-B CSV", start_dir, "CSV Files (*.csv);;All Files (*)"
        )
        if not path_str:
            return
        self._try_load(Path(path_str))

    def _try_load(self, path: Path, forget_on_failure=False):
        """Load a dataset, reporting a bad file rather than raising.

        `load_csv` validates before it commits, so a failure here leaves
        whatever was already open untouched -- the user just gets a warning.
        A recent-menu entry that fails is dropped from the list: it is a
        shortcut to a file, and a shortcut that does not work is worse than
        no shortcut.
        """
        try:
            self.load_csv(path)
        except (OSError, ValueError, pd.errors.ParserError) as e:
            QMessageBox.warning(self, "Load Data", f"Could not load {path.name}:\n{e}")
            if forget_on_failure:
                self._forget_recent(path)
            return False
        return True

    def _build_menu_bar(self):
        """A real menu bar, which on macOS is the system one at the top of the
        screen -- Open/Open Recent is where the OS trains people to look for
        this, and it brings ⌘O with it."""
        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction(self.open_action)
        self.file_recent_menu = file_menu.addMenu("Open Recent")
        file_menu.addSeparator()
        file_menu.addAction(self.load_config_action)
        file_menu.addAction(self.save_action)
        self._rebuild_recent_menus()

    def _remember_recent(self, csv_path: Path):
        """Move this dataset to the front of the recent list and persist it."""
        resolved = str(Path(csv_path).resolve())
        self.recent_files = ([resolved]
                             + [p for p in self.recent_files if p != resolved])
        del self.recent_files[RECENT_FILES_MAX:]
        self._save_state()
        self._rebuild_recent_menus()

    def _forget_recent(self, csv_path: Path):
        resolved = str(Path(csv_path).resolve())
        self.recent_files = [p for p in self.recent_files if p != resolved]
        self._save_state()
        self._rebuild_recent_menus()

    # Both of these are triggered FROM an action that the work itself deletes:
    # they end in _rebuild_recent_menus(), and QMenu.clear() destroys the
    # QAction whose `triggered` signal is still unwinding on the stack.
    # Returning into a deleted sender is a use-after-free, so the work is
    # deferred by one event-loop turn and the signal gets to finish first.
    def on_clear_recent_files(self):
        QTimer.singleShot(0, self._clear_recent_now)

    def _clear_recent_now(self):
        self.recent_files = []
        self._save_state()
        self._rebuild_recent_menus()

    def _save_state(self):
        """Write the local state file: the recent-files list, and nothing
        else. Called on every dataset load, which is exactly why it is not
        the shared config -- see DEFAULT_STATE_PATH."""
        try:
            save_config(self.state_path, {}, recent_files=self.recent_files)
        except OSError as e:
            print(f"Warning: could not write {self.state_path}: {e}")

    def _save_shared_config(self):
        """Write the shared, git-tracked config: the ICARTT header metadata,
        and nothing else. Only the Export tab's "Save defaults" button calls
        this, so a tracked file does not churn as flights are opened."""
        try:
            # No gas blocks: this file carries no analysis settings.
            save_config(self.default_config_path, {}, icartt_meta=self.icartt_meta)
        except OSError as e:
            print(f"Warning: could not write {self.default_config_path}: {e}")

    def _recent_labels(self):
        """Menu labels for the recent list: the file name, plus its parent
        directory only where the name alone would be ambiguous. Flights are
        named by date, so two directories holding `ucatsb-20250218All.csv` is
        exactly the case where the bare name is useless."""
        names = [Path(p).name for p in self.recent_files]
        return [f"{name} — {Path(p).parent.name}" if names.count(name) > 1 else name
                for p, name in zip(self.recent_files, names)]

    def _rebuild_recent_menus(self):
        """Repopulate both menus from self.recent_files.

        One builder for the button menu and the File > Open Recent submenu, so
        the two can never drift apart. Rebuilt on every change rather than
        patched, which is the same reasoning as refresh()'s unconditional
        cache invalidation: cheap, and it cannot serve a stale entry.
        """
        current = str(self.csv_path.resolve()) if self.csv_path else None
        labels = self._recent_labels()

        self.recent_menu.clear()
        self.recent_menu.addAction(self.open_action)
        self.recent_menu.addSeparator()
        self.file_recent_menu.clear()

        if not self.recent_files:
            for menu in (self.recent_menu, self.file_recent_menu):
                empty = menu.addAction("No recent files")
                empty.setEnabled(False)
            return

        for path_str, label in zip(self.recent_files, labels):
            exists = Path(path_str).exists()
            for menu in (self.recent_menu, self.file_recent_menu):
                action = menu.addAction(label)
                action.setToolTip(path_str)
                action.setStatusTip(path_str)
                # Kept but greyed rather than dropped: an unmounted volume
                # comes back, and silently forgetting the path would look like
                # the app losing track of the user's work.
                action.setEnabled(exists)
                if not exists:
                    action.setText(f"{label} (missing)")
                if path_str == current:
                    action.setCheckable(True)
                    action.setChecked(True)
                action.triggered.connect(
                    functools.partial(self._open_recent, path_str))

        for menu, label in ((self.recent_menu, "Clear Recent Files"),
                            (self.file_recent_menu, "Clear Menu")):
            menu.addSeparator()
            menu.addAction(label).triggered.connect(self.on_clear_recent_files)

    def _open_recent(self, path_str):
        """Deferred for the reason above on_clear_recent_files."""
        QTimer.singleShot(0, functools.partial(self._open_recent_now, path_str))

    def _open_recent_now(self, path_str):
        path = Path(path_str)
        if path == self.csv_path:
            self._reopen_config()
            return
        if not self._confirm_discard("loading another dataset"):
            return
        self._try_load(path, forget_on_failure=True)

    def _reopen_config(self):
        """Picking the already-open dataset: re-choose its configuration.

        This used to return silently, which left the config chooser reachable
        only by loading some other file first and coming back -- and a flight
        with several saved configs is exactly the case where switching between
        them is the point of the menu entry.

        The CSV is deliberately **not** re-read. Nothing about the data has
        changed, so re-parsing would spend a second arriving at the same frame
        and throw away `raw_df` and the flag row-count it validates against.
        Everything after the read is redone, which is what a config switch is:
        the same sequence `on_load_config_clicked` runs.
        """
        if not self._confirm_discard("reopening this dataset"):
            return
        self._load_flight_config(self.csv_path)
        self._select_gas(self.current_gas)
        self._sync_corr_cal_controls()
        self._update_tank_readout()
        self.refresh()

    def _build_controls(self):
        panel = QWidget()
        panel.setFixedWidth(CONTROLS_WIDTH)
        vbox = QVBoxLayout(panel)
        # Tighter than Qt's defaults (9/6). Nine group boxes' worth of gaps and
        # margins is ~70 px of nothing, which is real estate this panel does
        # not have; 6/4 still reads as separated.
        vbox.setContentsMargins(6, 6, 6, 6)
        vbox.setSpacing(4)

        # The menu IS the button's action (no clicked handler): pressing it
        # opens "Open File… / recent files / Clear", and setMenu also gives
        # the button its native dropdown indicator.
        load_row = QHBoxLayout()
        load_row.setSpacing(6)
        load_button = QPushButton("Load Data")
        load_button.setMenu(self.recent_menu)
        load_row.addWidget(load_button, 1)
        self.save_button = QPushButton(self.SAVE_LABEL)
        self.save_button.clicked.connect(self.on_save_clicked)
        load_row.addWidget(self.save_button)
        vbox.addLayout(load_row)

        # Just the name -- the panel is 300 px wide and a full path would
        # either clip or force the panel wider. The full path is the tooltip.
        # Two fixed lines (Filename:/Config:), written by _update_file_labels.
        # Deliberately not word-wrapped -- see there.
        self.file_label = QLabel("No file loaded")
        self.file_label.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.file_label)

        # Gas and the aux-panel pickers share one box: all four rows answer
        # "what is drawn", and two boxes cost 47 px of title-and-margin chrome
        # each. The aux selection is a combo rather than the four radio buttons
        # it used to be -- radios spend a row per option to show three the user
        # is not choosing, which is 100 px this panel cannot spare. It is a
        # single-choice control either way.
        traces_box = QGroupBox("Traces")
        traces_form = QFormLayout(traces_box)
        traces_form.setContentsMargins(6, 6, 6, 6)
        traces_form.setSpacing(4)
        traces_form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.gas_combo = QComboBox()
        self.gas_combo.addItems(self.available_gases.keys())
        self.gas_combo.currentTextChanged.connect(self.on_gas_changed)
        traces_form.addRow("Gas:", self.gas_combo)

        self.aux_combo = QComboBox()
        self.aux_combo.addItems(AUX_OPTIONS)
        self.aux_combo.setToolTip(
            "What to draw in a second panel above the gas trace.\n"
            '"No Figure" gives the gas the whole plot.'
        )
        self.aux_combo.currentTextChanged.connect(self.on_aux_changed)
        traces_form.addRow("Above:", self.aux_combo)

        self.other_combo = QComboBox()
        self.other_combo.addItems(self.other_columns)
        self.other_combo.setEnabled(False)
        self.other_combo.currentTextChanged.connect(self.on_other_changed)
        traces_form.addRow("Other:", self.other_combo)

        self.right_axis_combo = QComboBox()
        self.right_axis_combo.addItem("(none)")
        self.right_axis_combo.addItems(self.other_columns)
        self.right_axis_combo.setEnabled(False)
        self.right_axis_combo.currentTextChanged.connect(self.on_right_axis_changed)
        traces_form.addRow("Right axis:", self.right_axis_combo)

        # The note block's on/off switch, in this box because it is one more
        # answer to "what is drawn" -- and session-only, like the calibrated
        # overlay and Correlations' "Hide flagged points": it is a view toggle,
        # so it is not in DEFAULT_GAS_SETTINGS, not in _controls_to_settings(),
        # and does not dirty the config. Spanning row: a two-column form row
        # would put a label beside a control that already reads as a sentence.
        self.show_notes_check = QCheckBox("Info notes")
        self.show_notes_check.setChecked(True)
        self.show_notes_check.setToolTip(
            "Show the grey key at the bottom of the Timeseries figure -- what\n"
            "each shaded band means, how many readings were removed and why.\n\n"
            "Turn it off for a clean figure to save or hand on. The trace\n"
            "legend above it is not affected.\n\n"
            "Session only: it is not saved with the flight's settings, and it\n"
            "changes nothing about the data."
        )
        self.show_notes_check.toggled.connect(self.on_show_notes_toggled)
        traces_form.addRow(self.show_notes_check)

        vbox.addWidget(traces_box)

        self.mask_box = QGroupBox("Data Masking")
        mask_form = QFormLayout(self.mask_box)
        mask_form.setContentsMargins(6, 6, 6, 6)
        mask_form.setSpacing(4)
        mask_form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)

        # Warm-up and end-of-flight share a row: they are the same setting at
        # the two ends of the record, both in minutes, and they are already
        # OR-ed into one `trimmed` mask drawn as one orange band. Two form rows
        # spent 62 px saying that twice. The descent is the busiest part of a
        # flight and the least like the rest of it, so trimming the tail is as
        # routine as trimming the warm-up.
        # No " min" suffix on either spin box: the row label carries the unit
        # once for both. Two suffixes plus two inline labels plus the row label
        # do not fit across a 300 px panel -- they pushed the box's minimum
        # width to 347 and the panel silently clipped its own right-hand edge.
        self.warmup_spin = QSpinBox()
        self.warmup_spin.setRange(0, 120)
        self.warmup_spin.setSingleStep(1)
        self.warmup_spin.setToolTip(
            "Exclude this many minutes at the START of the record, while the\n"
            "instrument warms up. Reaches the cal means, not just the plot."
        )
        self.warmup_spin.valueChanged.connect(self.on_control_changed)

        self.end_flight_spin = QSpinBox()
        self.end_flight_spin.setRange(0, 120)
        self.end_flight_spin.setToolTip(
            "Exclude this many minutes at the END of the record, the mirror\n"
            "of the warm-up exclusion at the start -- descent, landing and\n"
            "whatever happens on the ground afterwards. Like the warm-up it\n"
            "reaches the cal means, not just the plot. 0 disables it."
        )
        self.end_flight_spin.valueChanged.connect(self.on_control_changed)

        trim_row = QHBoxLayout()
        trim_row.setSpacing(4)
        for text, spin in (("start", self.warmup_spin), ("end", self.end_flight_spin)):
            tag = QLabel(text)
            tag.setStyleSheet(f"color: {MUTED_COLOR};")
            trim_row.addWidget(tag)
            # Capped: a QSpinBox asks for ~85 px whatever it holds, and two of
            # them plus their tags plus the row label overflow the panel. These
            # hold at most "120", so the cap costs no visible digits -- but it
            # has to be a MAXIMUM, since that is the constraint that wins over
            # the size hint. Widen only alongside CONTROLS_WIDTH.
            spin.setMaximumWidth(55)
            trim_row.addWidget(spin, 1)
        mask_form.addRow("Trim (min):", trim_row)

        self.pressure_tol_spin = QDoubleSpinBox()
        self.pressure_tol_spin.setRange(0.0, 10.0)
        self.pressure_tol_spin.setSingleStep(0.05)
        self.pressure_tol_spin.setDecimals(2)
        # The " mbar" suffix became a separate tag when the correction arrived
        # beside it: suffix + spin + checkbox asks 325 px against a 312 px
        # panel, and a row too wide for the panel neither scrolls nor wraps --
        # it clips its own right-hand edge in silence. Capped for the same
        # reason as the trim pair (a spin box asks ~85 px whatever it holds);
        # 70 is what its widest value, "10.00", needs -- measured, not
        # guessed, and 66 visibly ate the last digit.
        self.pressure_tol_spin.setMaximumWidth(70)
        self.pressure_tol_spin.valueChanged.connect(self.on_control_changed)

        # Both corrections sit on a "Correct:" row of their own (2026-07-31).
        # 140/P used to share the tolerance's row, which was right while it was
        # the only one -- it shares that row's input. With T-40 beside it the
        # pair belongs together more than either belongs to the tolerance, and
        # the arithmetic settles it anyway: tolerance + unit + two checkboxes
        # asks ~330 px against the ~300 the group box has, and a row too wide
        # for the panel clips its own right-hand edge in silence.
        #
        # Short labels keep both correction checkboxes inside the fixed-width
        # panel; the tooltips carry the arithmetic.
        self.pressure_correct_check = QCheckBox("Pdelt")
        self.pressure_correct_check.setToolTip(
            f"Scale the MEASUREMENT by {D1_P_TARGET_MBARS:.0f}/P before "
            f"calibrating it, normalising every\nreading to the detector's "
            f"{D1_P_TARGET_MBARS:.0f} mbar spec pressure. P is that gas's own "
            f"detector\npressure (d1 for CO2/N2O, d2 for CH4), so this only "
            f"exists for the Aeris gases.\n\n"
            f"Applied to the whole flight, cal periods included -- so the "
            f"cal-bottle means\nare corrected too, and this DOES move the "
            f"calibration: slope, intercept,\nspan gain and the residuals all "
            f"change when it is toggled. That is the point.\nA cal injection "
            f"measured at 138 mbar is put on the same footing as ambient\nair "
            f"measured at 140 before it becomes a calibration node.\n\n"
            f"Off by default, so an existing config's numbers do not change "
            f"until it is asked for."
        )
        self.pressure_correct_check.toggled.connect(self.on_control_changed)

        # One line, not two: the second line only restated the target, which
        # the tooltip now carries. A wrapped label costs a full 20 px row.
        pressure_label = QLabel("Pressure:")
        pressure_label.setToolTip(
            f"Exclude data whose detector pressure is further than this from\n"
            f"the {D1_P_TARGET_MBARS:.0f} mbar target."
        )
        # The unit as a muted inline tag, like the trim row's start/end: it
        # cannot go in the row label (widening the label column moves every
        # row and costs more than it saves) and it cannot stay in the spin box
        # (the suffix is what pushed this row over the panel's width).
        pressure_unit = QLabel("mbar")
        pressure_unit.setStyleSheet(f"color: {MUTED_COLOR};")
        pressure_row = QHBoxLayout()
        pressure_row.setSpacing(4)
        pressure_row.addWidget(self.pressure_tol_spin)
        pressure_row.addWidget(pressure_unit)
        pressure_row.addStretch(1)
        mask_form.addRow(pressure_label, pressure_row)

        # Its own row rather than joining the one above: that row already asks
        # ~267 px of the ~300 the group box has, and a spin box asks ~85
        # whatever it holds. A row too wide for the panel clips its own
        # right-hand edge in silence.
        self.pressure_smooth_spin = QSpinBox()
        self.pressure_smooth_spin.setRange(0, 300)
        self.pressure_smooth_spin.setSuffix(" s")
        self.pressure_smooth_spin.setToolTip(
            f"Smooth the detector cell pressure over this many seconds before\n"
            f"the {D1_P_TARGET_MBARS:.0f}/P correction divides by it, so the "
            f"pressure sensor's own\nnoise is not added to the mole fraction "
            f"(~0.07 mbar of scatter on d1 is\n~0.2 ppm of CO2). The window is "
            f"a CENTRED mean, so it introduces no\ntime shift, and it runs "
            f"over the whole flight -- air and cal alike, since\nthe "
            f"correction applies to both.\n\n"
            f"Keep it short (10-30 s gets most of the noise reduction there is "
            f"to get). A\nlong window reads between the two levels either side "
            f"of a solenoid\ntransition, where the cell steps by ~2 mbar.\n\n"
            f"The smoothed trace is drawn in red over the raw one whenever "
            f"'Above:' is\nshowing Detector Pressure, so the window can be "
            f"judged against the real\nexcursions it has to keep.\n\n"
            f"0 disables it, so an existing config's numbers do not change "
            f"until it is\nasked for."
        )
        self.pressure_smooth_spin.setMaximumWidth(70)
        self.pressure_smooth_spin.valueChanged.connect(self.on_control_changed)
        smooth_p_row = QHBoxLayout()
        smooth_p_row.setSpacing(4)
        smooth_p_row.addWidget(self.pressure_smooth_spin)
        smooth_p_row.addStretch(1)
        smooth_p_label = QLabel("Smooth P:")
        smooth_p_label.setToolTip(self.pressure_smooth_spin.toolTip())
        mask_form.addRow(smooth_p_label, smooth_p_row)

        # The temperature companion to 140/P. Deliberately NOT smoothed, on the
        # PI's instruction and because the reading does not need it: d1_T_gas
        # moves 0.34 C across
        # the whole Jul 2026 flight where the pressure carries 0.16 mbar of
        # sample-to-sample noise. There is nothing to filter out, so there is
        # no "Smooth T" beside it.
        self.temperature_correct_check = QCheckBox("Tdelt")
        self.temperature_correct_check.setToolTip(
            f"Scale the MEASUREMENT by 1 + (T_K - {T_GAS_REFERENCE_K:.2f})/"
            f"{T_GAS_TARGET_K:.0f}\nbefore calibrating it. T is that gas's own "
            f"detector cell temperature\n(d1_T_gas for CO2/N2O, d2_T_gas for "
            f"CH4), in degrees C.\n\n"
            f"Number density goes as P/T, so a cell running hotter than "
            f"{T_GAS_REFERENCE_C:.0f} C ({T_GAS_REFERENCE_K:.2f} K)\nholds less gas and the measurement "
            f"is scaled up. A row with no reading\ngets no value, like 140/P.\n\n"
            f"Like 140/P it is applied to the whole flight before the "
            f"calibration, so it\nmoves the cal means and the calibration with "
            f"them. The two multiply when\nboth are on, and reach the "
            f"calibrated trace, the Correlations tab and both\nexports.\n\n"
            f"Off by default, so an existing config's numbers do not change "
            f"until it is\nasked for."
        )
        self.temperature_correct_check.toggled.connect(self.on_control_changed)

        correct_row = QHBoxLayout()
        correct_row.setSpacing(6)
        correct_row.addWidget(self.pressure_correct_check)
        correct_row.addWidget(self.temperature_correct_check)
        correct_row.addStretch(1)
        correct_label = QLabel("Correct:")
        correct_label.setToolTip(
            "Normalise the calibrated mole fraction to the detector cell's\n"
            "spec conditions. Both are post-multipliers on the calibrated\n"
            "product and neither can move the calibration itself."
        )
        mask_form.addRow(correct_label, correct_row)

        # Unlike the two masks above, this one does not touch the cal means
        # (the flush window is ambient by definition, never inside a cal
        # period) -- it only blanks the calibrated product. Defaults to 0 so
        # an existing config's output doesn't change until it is asked for.
        self.flag_air_spin = QSpinBox()
        # Upper limit is 90 s, not 30: on the Jul 2026 flight the cells were
        # still ~1 ppm low at 30 s and took ~45 s to settle within 0.2 ppm.
        self.flag_air_spin.setRange(0, 90)
        self.flag_air_spin.setSuffix(" s")
        self.flag_air_spin.setToolTip(
            "Air data immediately after a cal injection still reads toward the\n"
            "tank while the detector cells flush. This many seconds after each\n"
            "cal period end are dropped from the calibrated series and the\n"
            "export; the raw trace keeps them. 0 disables it."
        )
        self.flag_air_spin.setMaximumWidth(62)
        self.flag_air_spin.valueChanged.connect(self.on_control_changed)

        # Pumps toggle and Flag Air share one spanning row, the toggle to the
        # left of the spin box. Spanning rather than a two-column form row for
        # the reason the drift-model row is: a stretching field squeezes the
        # label column to nothing in a 300 px panel.
        air_row = QHBoxLayout()
        air_row.setSpacing(6)
        self.pumps_check = QCheckBox("Pumps on")
        self.pumps_check.setToolTip(
            "Keep only data recorded with the sample pumps running\n"
            "(j_pumps = 1). Air measured with the pumps off is not\n"
            "ambient air.\n\n"
            "Off by default, and it has to be: a lab test or bench\n"
            "calibration runs with the pumps off from end to end (the\n"
            "2026-07-26 file is 100% pumps-off), and enabling this there\n"
            "would leave nothing at all. Turn it on for a real flight.\n\n"
            "Rows with no j_pumps value count as pumps-off -- an unknown\n"
            "pump state is not evidence the pumps were running."
        )
        self.pumps_check.toggled.connect(self.on_control_changed)
        air_row.addWidget(self.pumps_check)
        air_row.addStretch(1)
        air_row.addWidget(QLabel("Flag Air:"))
        air_row.addWidget(self.flag_air_spin)
        mask_form.addRow(air_row)

        # Warm-up, detector pressure and the cal timing are properties of the
        # instrument on this flight rather than of the species, so the same
        # values usually want to apply to every gas -- but the settings stay
        # per-gas, so that a gas that does need its own can still have it.
        # This button is the bridge. It lives in this box but reaches the cal
        # window boxes below too, which the label and tooltip have to say.
        self.copy_mask_button = QPushButton(self.COPY_SETTINGS_LABEL)
        self.copy_mask_button.setToolTip(
            "Copy warm-up, pressure tolerance and correction, Flag Air and\n"
            "both cal mean windows from the current gas to every other\n"
            "calibrated gas (CO2/N2O/CH4). The drift model and smoothing\n"
            "window are left alone."
        )
        self.copy_mask_button.clicked.connect(self.on_copy_masking_to_all)
        mask_form.addRow(self.copy_mask_button)

        vbox.addWidget(self.mask_box)

        # ONE calibration box, holding what used to be three: a mean-window box
        # per tank plus a box for the drift model. Each cost 41 px of title and
        # margin to hold a single row, and all three answer the same question
        # -- how the two-point calibration is built. The grid rows are the two
        # tanks' windows; the drift model spans beneath them.
        self.cal_box = QGroupBox("Calibration")
        cal_grid = QGridLayout(self.cal_box)
        cal_grid.setContentsMargins(6, 6, 6, 6)
        cal_grid.setHorizontalSpacing(4)
        cal_grid.setVerticalSpacing(4)

        # Visual order is unchanged from the two-box version: the 100% bottle's
        # row is above the 50% one. The tank each row belongs to used to be a
        # box title, which had room for "100% Cal (CC302489) 418.947 ppm"; a row
        # label does not, so it carries the short form and the full
        # identification moved to its tooltip -- see _set_cal_row_label, which
        # sets both together. Start/End get no labels of their own for the same
        # reason; the header row above them names the columns once.
        cal_grid.addWidget(QLabel("Window:"), 0, 0)
        for col, name in ((1, "start"), (2, "end")):
            header = QLabel(name)
            header.setAlignment(Qt.AlignHCenter)
            header.setStyleSheet(f"color: {MUTED_COLOR};")
            cal_grid.addWidget(header, 0, col)
        (self.cal2_label, self.cal2_start_spin,
         self.cal2_end_spin) = self._add_cal_window_row(cal_grid, 1, "Cal 2")
        (self.cal1_label, self.cal1_start_spin,
         self.cal1_end_spin) = self._add_cal_window_row(cal_grid, 2, "Cal 1")

        # One row, one label: the smoothing window only means anything for the
        # "smooth" model, so it reads as part of that choice rather than as an
        # independent setting, and it greys out with the model set to anything
        # else. Its own label would just repeat what the suffix already says.
        # Label carried inside the row and added as a spanning row, not via
        # addRow(str, layout): the control panel is a fixed 300 px, and in the
        # two-column form the label column got squeezed to nothing by the
        # combo's stretch -- the label silently vanished.
        drift_row = QHBoxLayout()
        drift_row.setSpacing(4)
        # "Model:" rather than "Drift model:" since "fixed slope" joined the
        # list: that entry says nothing about the drift nodes, it pins the gain,
        # and inside a box already titled "Calibration" the one word is enough.
        # The 46 px it gives back is not decoration either -- "fixed slope" is
        # 83 px of text plus a dropdown arrow, and the slope spin needs room for
        # four decimals. At "Drift model:" the row asked 305 of the 312 the
        # panel has and the spin silently rendered 1.0367 as "1.036". The stored
        # key is still `drift_model`, so no saved config is disturbed.
        drift_tip = (
            "How the per-injection cal means become a calibration.\n\n"
            "linear / smooth / constant shape each bottle's response in time\n"
            "and solve the two-point line at every sample. 'fixed slope' pins\n"
            "the gain and lets each cal set an intercept anchor."
        )
        drift_label = QLabel("Model:")
        drift_label.setToolTip(drift_tip)
        drift_row.addWidget(drift_label)
        self.drift_combo = QComboBox()
        self.drift_combo.addItems(CAL_DRIFT_MODELS)
        self.drift_combo.setToolTip(drift_tip)
        self.drift_combo.currentTextChanged.connect(self.on_control_changed)
        drift_row.addWidget(self.drift_combo, 1)

        self.smooth_spin = QSpinBox()
        self.smooth_spin.setRange(2, 21)
        # "ev" rather than " events": the row has a label, a combo and this
        # spin box to fit into a 300 px panel, and the full word clipped. The
        # tooltip carries the meaning.
        self.smooth_spin.setSuffix(" ev")
        self.smooth_spin.setToolTip(
            "Width of the centred rolling mean over cal events, used by the\n"
            '"smooth" drift model only.'
        )
        self.smooth_spin.valueChanged.connect(self.on_control_changed)
        self.smooth_spin.setMaximumWidth(72)
        drift_row.addWidget(self.smooth_spin)

        # Shares the smoothing window's slot rather than taking one of its own:
        # each model needs at most one number, and the label + combo + TWO spin
        # boxes do not fit the 300 px panel. Exactly one is visible at a time
        # (_show_drift_extra), so the row never grows.
        self.fixed_slope_spin = QDoubleSpinBox()
        self.fixed_slope_spin.setRange(0.0, 10.0)
        self.fixed_slope_spin.setDecimals(4)
        self.fixed_slope_spin.setSingleStep(0.001)
        # 0 reads "auto" rather than as a slope of zero, which would calibrate
        # every reading to a constant. It is the value that means "whatever the
        # constant model gives", and typing it is how you get back there after
        # nudging the number.
        self.fixed_slope_spin.setSpecialValueText("auto")
        self.fixed_slope_spin.setToolTip(
            'The gain held fixed by the "fixed slope" model, with the\n'
            "intercept anchored to each cal mean and interpolated in time.\n\n"
            "Seeded from the slope the constant model gives (= 1/span gain)\n"
            "when you pick the model, and editable from there. 0 shows as\n"
            '"auto" and means exactly that seed, recomputed per flight.\n\n'
            "Each cal injection is forced to match its assigned tank value;\n"
            "the fixed slope controls how the response is carried between\n"
            "those cal anchors."
        )
        self.fixed_slope_spin.valueChanged.connect(self.on_control_changed)
        # 88, not the 72 the smoothing spin uses: four decimals plus the step
        # arrows need it, and at 72 the box rendered 1.0367 as "1.036" -- the
        # exact silent-truncation failure the panel-width note warns about.
        # Measured against the rendered widget, not guessed.
        self.fixed_slope_spin.setMaximumWidth(88)
        self.fixed_slope_spin.setVisible(False)
        drift_row.addWidget(self.fixed_slope_spin)
        self.fixed_slope_reset_button = QToolButton()
        self.fixed_slope_reset_button.setText("R")
        self.fixed_slope_reset_button.setFixedSize(28, 22)
        self.fixed_slope_reset_button.setStyleSheet(
            "QToolButton { padding: 0; color: #222222; }")
        self.fixed_slope_reset_button.setToolTip(
            "Reset fixed slope to the constant-model value")
        self.fixed_slope_reset_button.clicked.connect(self.on_fixed_slope_reset)
        self.fixed_slope_reset_button.setVisible(False)
        drift_row.addWidget(self.fixed_slope_reset_button)
        cal_grid.addLayout(drift_row, 3, 0, 1, 3)

        # "Show calibrated" is NOT here any more -- it is a checkable action on
        # the Timeseries toolbar (see PlotPane.calibrated_action). It changes
        # no setting and is not saved; it only decides what the figure draws,
        # which is what every other control on that toolbar does. It stays
        # session-only and off by default for the old reason: the timeseries is
        # documented as showing uncalibrated data, and a remembered toggle
        # would let the app start up showing calibrated data with no visible
        # reason why.
        #
        # No export button here either: exporting is the Export tab's job,
        # where it can cover every gas at once. This panel is per-gas, and a
        # per-gas button was quietly the reason the old export could only ever
        # describe one of them.
        vbox.addWidget(self.cal_box)

        # Deliberately NOT in the setEnabled(has_masking) list in _select_gas:
        # Ozone is the gas this feature exists for, and it is precisely the
        # one with no masking settings to enable. Its own group box for the
        # same reason -- it is not a masking *setting*, it is a record of
        # points the user struck out by hand.
        self.flag_box = QGroupBox("Flagged Points")
        flag_form = QVBoxLayout(self.flag_box)
        flag_form.setContentsMargins(6, 6, 6, 6)
        flag_form.setSpacing(4)
        self.flag_label = QLabel("No points flagged")
        self.flag_label.setWordWrap(True)
        self.flag_label.setStyleSheet(f"color: {MUTED_COLOR};")
        flag_form.addWidget(self.flag_label)

        self.flag_all_check = QCheckBox("Apply to all gases")
        self.flag_all_check.setToolTip(
            "Spread each new flag across every gas — for an inlet or pump\n"
            "problem that ruins them all at once. Unticked, a flag belongs\n"
            "to the gas on display. Does not affect Clear, which is always\n"
            "the current gas only."
        )
        flag_form.addWidget(self.flag_all_check)

        flag_row = QHBoxLayout()
        self.flag_undo_button = QPushButton("Undo")
        self.flag_undo_button.setToolTip(
            "Step back through this session's flagging. Not saved with the\n"
            "flight — reopening starts a fresh history."
        )
        self.flag_undo_button.clicked.connect(self.on_flag_undo)
        self.flag_clear_button = QPushButton("Clear")
        self.flag_clear_button.setToolTip("Remove every flag on the current gas.")
        self.flag_clear_button.clicked.connect(self.on_flag_clear)
        flag_row.addWidget(self.flag_undo_button)
        flag_row.addWidget(self.flag_clear_button)
        flag_form.addLayout(flag_row)
        vbox.addWidget(self.flag_box)

        vbox.addStretch(1)
        return panel

    # Carried as a tooltip rather than a standing label -- it's reference
    # material you need once, and the control panel has to fit on a laptop.
    CAL_WINDOW_HELP = (
        "Relative to the last point in a cal period (Cal_p),\n"
        "e.g. -10 s to 2 s = [Cal_p-10s, Cal_p+2s].\n"
        "Positive values reach past Cal_p.\n"
        "Settings are saved per-gas."
    )

    def _add_cal_window_row(self, grid, row, fallback):
        """One tank's [start, end] window as a labelled row of the shared grid.

        Start and End have no labels of their own -- the box title says
        "(start, end)" and each spin box carries its own tooltip. Two words per
        row would push the tank label out of a 300 px panel, and the pair reads
        as one interval anyway.
        """
        label = QLabel(fallback)
        label.setToolTip(self.CAL_WINDOW_HELP)
        grid.addWidget(label, row, 0)

        spins = []
        for which, tip in (("start", "Window START, relative to Cal_p."),
                           ("end", "Window END, relative to Cal_p.")):
            spin = QSpinBox()
            spin.setRange(-60, 60)
            spin.setSuffix(" s")
            # See the trim spins: a maximum, not a hint, so the pair fits the
            # panel. "-60 s" is the widest text these ever hold.
            spin.setMaximumWidth(80)
            spin.setToolTip(f"{tip}\n\n{self.CAL_WINDOW_HELP}")
            spin.valueChanged.connect(self.on_control_changed)
            grid.addWidget(spin, row, len(spins) + 1)
            spins.append(spin)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(2, 1)
        return label, spins[0], spins[1]

    def _set_cal_row_label(self, label_widget, serial, fallback):
        """Name a cal-window row for the tank it belongs to.

        The visible text is the short form -- cals.yaml's `info` ("100%") when
        there is one, else the serial -- because a row label has far less room
        than the group-box title this used to be. The full identification
        ("100% Cal (CC302489) 418.947 ppm") moves to the tooltip, so nothing is
        lost, only demoted. `info` stays optional: not every tank in the roster
        has one, and the next one added may well arrive without it.
        """
        nominal = self.cal_bottles.get(serial, {}) if serial else {}
        if not nominal:
            label_widget.setText(fallback)
            label_widget.setToolTip(self.CAL_WINDOW_HELP)
            return
        info = nominal.get("info")
        full = f"{info} Cal ({serial})" if info else f"Cal ({serial})"
        value = nominal.get(self.current_gas)
        if value is not None:
            full += f" {value:g} {gas_unit(self.current_gas)}"
        label_widget.setText(info if info else serial)
        label_widget.setToolTip(f"{full}\n\n{self.CAL_WINDOW_HELP}")

    def _build_tabs(self):
        """Timeseries, Calibration and Cal Tanks as tabs over the shared
        controls.

        The controls panel deliberately stays outside the tabs: every control
        affects both plot views, so moving them inside would mean duplicating
        the gas selector or making one tab depend on state invisible from the
        other. Cal Tanks is the exception that proves it -- the tank pairing
        is one per flight, not one per gas or per view, so it has nowhere
        sensible to live in a per-gas control panel.
        """
        self.main_pane = PlotPane()
        self.cal_pane = PlotPane()
        self.main_pane.on_box = self.on_stats_box
        self.main_pane.on_flag_box = self.on_flag_box
        self.main_pane.calibrated_action.toggled.connect(self.on_calibrated_toggled)
        # The Calibration tab's three panels each mean something different
        # (response deviation, coefficients, residuals), so a single box-stats
        # readout there would be ambiguous -- Timeseries only for now. Flagging
        # is Timeseries-only for a stronger reason: the cal panels plot derived
        # quantities (deviations, coefficients, residuals), not the rows a flag
        # would have to name.
        self.cal_pane.stats_action.setVisible(False)
        self.cal_pane.flag_action.setVisible(False)
        # The overlay is a Timeseries thing: the other two panes draw derived
        # quantities that are calibrated already or not calibrated at all, so
        # there is no raw trace there to lay a calibrated one over.
        self.cal_pane.calibrated_action.setVisible(False)
        # Keep the historical attribute names bound to the timeseries pane so
        # redraw()'s existing body needs no changes.
        self.figure = self.main_pane.figure
        self.canvas = self.main_pane.canvas
        self.toolbar = self.main_pane.toolbar

        self.tanks_pane = self._build_cal_tanks_pane()
        self.corr_pane = PlotPane()
        # Stats stays hidden here: its readout describes one trace over a time
        # span, which a tracer-tracer scatter is not. Flagging is the opposite
        # case and the reason this tab wanted the tool -- an outlier that is
        # obvious against another tracer can be near-impossible to find in the
        # timeseries. The ambiguity Stats could not resolve (a box names two
        # gases' rows) is settled by the explicit target combo instead.
        self.corr_pane.stats_action.setVisible(False)
        self.corr_pane.calibrated_action.setVisible(False)
        self.corr_pane.on_flag_box = self.on_corr_flag_box
        self.corr_pane.canvas.mpl_connect("button_press_event", self.on_corr_tooltip_press)

        self.export_pane = self._build_export_pane()

        self.tabs = QTabWidget()
        self.tabs.addTab(self.main_pane, "Timeseries")
        self.tabs.addTab(self.cal_pane, "Calibration")
        self.tabs.addTab(self.corr_pane, "Correlations")
        self.tabs.addTab(self.tanks_pane, "Cal Tanks")
        self.tabs.addTab(self.export_pane, "Export")
        # Connected after every addTab call -- the first one fires
        # currentChanged before the other tabs exist.
        self.tabs.currentChanged.connect(self.on_tab_changed)
        return self.tabs

    def _build_cal_tanks_pane(self):
        """Which two roster tanks this flight flew.

        Its own tab rather than another group box in the control panel: the
        choice is per *flight*, while everything in that panel is per gas, and
        getting it wrong invalidates every calibrated number for every gas at
        once -- so it gets room for the roster values that let you check it.
        """
        pane = QWidget()
        vbox = QVBoxLayout(pane)

        intro = QLabel(
            "Which two cal tanks were plumbed in for this flight. Defaults to "
            "the pair named in cals.yaml, which describes the tanks in use "
            "now — an older flight almost certainly flew different ones. The "
            "choice is saved in this dataset's own <dataset>_conf.yaml, so it "
            "travels with the flight."
        )
        intro.setWordWrap(True)
        vbox.addWidget(intro)

        box = QGroupBox("Tanks on this flight")
        form = QFormLayout(box)
        self.tank_combos = {}
        for key in ("cal0", "cal1"):
            combo = QComboBox()
            # Order is cosmetic: match_cal_serial identifies the tank in each
            # cal window by measured concentration, so cal0/cal1 name the SET
            # of two tanks, not which solenoid state each is wired to. Said
            # out loud here because the labels invite the opposite assumption.
            combo.setToolTip(
                "cal0/cal1 mirror the cals.yaml key names. Which of the two "
                "is flowing in a given cal window is identified from the "
                "measured concentration, not from this order."
            )
            combo.currentIndexChanged.connect(
                functools.partial(self.on_cal_tank_changed, key))
            form.addRow(f"{key}:", combo)
            self.tank_combos[key] = combo
        vbox.addWidget(box)

        self.tank_warning = QLabel()
        self.tank_warning.setWordWrap(True)
        self.tank_warning.setStyleSheet(f"color: {PRESSURE_EXCLUDE_COLOR};")
        vbox.addWidget(self.tank_warning)

        values_box = QGroupBox("Assigned values (cals.yaml)")
        values_layout = QVBoxLayout(values_box)
        self.tank_values_label = QLabel()
        self.tank_values_label.setTextFormat(Qt.PlainText)
        self.tank_values_label.setStyleSheet("font-family: Menlo, monospace;")
        values_layout.addWidget(self.tank_values_label)
        vbox.addWidget(values_box)

        reset_row = QHBoxLayout()
        reset_button = QPushButton("Reset to cals.yaml default")
        reset_button.clicked.connect(self.on_reset_cal_tanks)
        reset_row.addWidget(reset_button)
        reset_row.addStretch(1)
        vbox.addLayout(reset_row)

        self.tank_conf_label = QLabel()
        self.tank_conf_label.setWordWrap(True)
        self.tank_conf_label.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.tank_conf_label)

        vbox.addStretch(1)
        self._populate_tank_combos()
        self._apply_cal_selection_to_controls()
        return pane

    def _populate_tank_combos(self):
        """Fill both combos with the whole roster. Every tank ever used is
        offered, not just the two in cals.yaml's `cals:` block -- picking a
        tank that is not the current pairing is the entire point of the tab."""
        for combo in self.tank_combos.values():
            combo.blockSignals(True)
            combo.clear()
            for serial in sorted(self.cal_roster):
                info = self.cal_roster[serial].get("info")
                combo.addItem(f"{serial} ({info})" if info else serial, serial)
            combo.blockSignals(False)

    def _apply_cal_selection_to_controls(self):
        """Sync the Cal Tanks tab to self.cal_selection without re-triggering
        a save/redraw per combo (same guard pattern as the control panel's
        _apply_settings_to_controls)."""
        if not hasattr(self, "tank_combos"):
            return
        loading = self._loading
        self._loading = True
        for key, combo in self.tank_combos.items():
            serial = self.cal_selection.get(key)
            index = combo.findData(serial)
            # A conf file naming a tank the roster has since lost leaves the
            # combo where it is rather than silently retargeting it; the
            # warning label below says the selection is incomplete.
            if index >= 0:
                combo.setCurrentIndex(index)
        self._loading = loading
        self._update_tank_readout()

    def _update_tank_readout(self):
        """Warning + assigned-value table for the current pairing."""
        serials = [self.cal_selection.get(k) for k in ("cal0", "cal1")]
        missing = [s for s in serials if s not in self.cal_roster]
        warnings = []
        if missing:
            warnings.append(
                f"Not in cals.yaml's roster: {', '.join(str(s) for s in missing)}. "
                "Add the tank to cals.yaml, or pick another.")
        elif serials[0] == serials[1]:
            warnings.append(
                "Both states point at the same tank, so there is no span: the "
                "calibration degrades to an offset-only correction.")
        for gas in self.available_gases or GASES:
            if not GASES[gas].get("has_masking", True):
                continue
            without = [s for s in serials if s in self.cal_roster
                       and self.cal_roster[s].get(gas) is None]
            if without:
                warnings.append(
                    f"No assigned {gas} value for {', '.join(without)} — "
                    f"{gas} cannot be calibrated with this pairing.")
        self.tank_warning.setText("\n".join(warnings))
        self.tank_warning.setVisible(bool(warnings))

        rows = [f"{'gas':<6}" + "".join(f"{str(s or '-'):>22}" for s in serials)]
        for gas, info in GASES.items():
            if not info.get("has_masking", True):
                continue
            unit = gas_unit(gas)
            cells = []
            for serial in serials:
                nominal = self.cal_roster.get(serial, {})
                value = nominal.get(gas)
                if value is None:
                    cells.append(f"{'--':>22}")
                else:
                    unc = nominal.get(f"{gas}_unc")
                    text = f"{value:g}" + (f" ± {unc:g}" if unc is not None else "")
                    cells.append(f"{text + ' ' + unit:>22}")
            rows.append(f"{gas:<6}" + "".join(cells))
        self.tank_values_label.setText("\n".join(rows))

        if self.csv_path is None:
            self.tank_conf_label.setText(
                "No dataset loaded — this is cals.yaml's default pairing. "
                "Load a flight to save a choice against it.")
        elif self.config_loaded_from is not None:
            self.tank_conf_label.setText(
                f"From {self.config_loaded_from.name} — Save to keep changes")
        else:
            self.tank_conf_label.setText(
                f"Not saved yet — Save… writes {self.config_path.name}")

    # The ICARTT header metadata form: (key, label, height, tooltip), where
    # height 0 is a single-line box and anything else a text area that tall.
    # Ordered as the file itself is -- the four header lines, then the
    # required normal-comment keywords in the sequence the standard fixes --
    # so the panel reads like the file it produces and a field can be checked
    # against a delivered header without hunting.
    #
    # Most multi-line boxes are a convenience only: their value is flattened
    # to one line on write (`_one_line`), since each keyword occupies exactly
    # one line and a stray newline would invalidate the header's line count.
    # **`special_comments` and `revision_history` are the exceptions** -- they
    # are whole sections rather than keyword values, and their line breaks are
    # written through verbatim (`_verbatim_lines`), which is why they get the
    # tallest boxes.
    ICARTT_FIELDS = (
        ("data_id", "Data ID", 0,
         "First part of the file name — the mission/instrument data product. "
         "The instrument is registered with the archive as SABRE-RASTA "
         "(renamed from SABRE-UCATSB in 2026 to end the confusion with the "
         "UCATS gas chromatograph), which is what this should stay for "
         "SABRE. The ICARTT standard requires the "
         "ID to match whatever the archiving data center has registered, so "
         "a new campaign means asking them, not inventing one. Hyphens are "
         "kept; underscores separate the file name's own fields and are "
         "stripped."),
        ("location_id", "Location ID", 0,
         "Second part of the file name — normally the platform, e.g. WB57."),
        ("pi_name", "PI name(s)", 0,
         "Header line 2. LAST, FIRST — several PIs separated by semicolons, "
         "as in HINTSA, ERIC; MOORE, FRED."),
        ("pi_affiliation", "PI affiliation", 0, "Header line 3."),
        ("data_source", "Data source", 0,
         "Header line 4 — the instrument or measurement this file reports."),
        ("mission", "Mission", 0, "Header line 5 — the campaign name."),
        ("pi_contact_info", "PI contact info", 72,
         "Address, phone, email of the PI."),
        ("platform", "Platform", 0, "The aircraft or other platform."),
        ("location", "Location", 0,
         "Where the data were taken. The delivered UCATS files point at the "
         "nav file here: \"Aircraft location in separate file\"."),
        ("associated_data", "Associated data", 72,
         "Other files needed to use this one (e.g. MMS-1HZ), or N/A."),
        ("instrument_info", "Instrument info", 96,
         "How the measurement is made."),
        ("data_info", "Data info", 96,
         "Anything a user needs to know to use the numbers correctly — the "
         "UCATS files state the units of each species here."),
        ("uncertainty", "Uncertainty (extra notes)", 72,
         "The median 1-sigma per gas is computed and written automatically.\n"
         "Anything typed here is appended to it."),
        ("ulod_value", "ULOD value", 0,
         "Upper limit of detection, or N/A. The flag value itself is fixed by "
         "the format."),
        ("llod_value", "LLOD value", 0,
         "Lower limit of detection, or N/A."),
        ("dm_contact_info", "DM contact info", 72,
         "The data manager, if that is not the PI."),
        ("project_info", "Project info", 96,
         "The project, its dates and its purpose."),
        ("stipulations_on_use", "Stipulations on use", 96,
         "Terms the PI puts on using the data."),
        ("other_comments", "Other comments", 96, "Anything else."),
        ("revision", "Revision", 0,
         "The revision THIS file is. RA for preliminary, R0 for the first "
         "revised release, R1 next, and so on. Goes in the file name too."),
        ("revision_history", "Revision history", 120,
         "One `R#: description` per line, written into the file verbatim.\n"
         "It ACCUMULATES — an R0 file still lists its RA line above the R0 "
         "one.\nLeave empty and the current revision gets a placeholder line."),
        ("special_comments", "Special comments", 300,
         "Free text written into the file verbatim, line breaks and blank "
         "lines included.\nThis is where the delivered UCATS files explain "
         "their error estimates and ask users to contact the PIs.\n"
         "This box is the whole section — nothing is added to it."),
        ("var_suffix", "Variable name suffix", 0,
         "Appended to each species name, e.g. CO2_RASTA; the 1-sigma "
         "variable becomes CO2e_RASTA. Convention is the instrument or PI."),
    )

    def _build_export_pane(self):
        """The Export tab: the two delivery products, and control over what
        goes into them.

        Its own tab rather than a button in the control panel, for the same
        reason as Cal Tanks: everything in that panel is per gas, while both
        of these files cover the whole flight and every gas at once. The old
        per-gas "Export calibrated CSV" button is gone -- it could only ever
        describe the gas that happened to be selected.
        """
        pane = QWidget()
        outer = QVBoxLayout(pane)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        scroll.setWidget(inner)
        outer.addWidget(scroll)

        summary_box = QGroupBox("What will be exported")
        summary_layout = QVBoxLayout(summary_box)
        self.export_summary_label = QLabel("No file loaded.")
        self.export_summary_label.setTextFormat(Qt.PlainText)
        self.export_summary_label.setStyleSheet("font-family: Menlo, monospace;")
        self.export_summary_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        summary_layout.addWidget(self.export_summary_label)
        vbox.addWidget(summary_box)

        # ---- companion CSV
        csv_box = QGroupBox("Derived CSV (companion to the raw file)")
        csv_layout = QVBoxLayout(csv_box)
        csv_intro = QLabel(
            "Every gas, and one row for every row of the raw CSV — including "
            "the pre-sync rows, blank — so the two files line up exactly and "
            "can be opened side by side, or pasted together, in Excel or Igor. "
            "Masks are written as 1/0."
        )
        csv_intro.setWordWrap(True)
        csv_layout.addWidget(csv_intro)

        self.csv_raw_check = QCheckBox("Include the raw value columns")
        self.csv_raw_check.setChecked(True)
        self.csv_raw_check.setToolTip(
            "The raw columns are already in the source file. Keeping them "
            "here makes this file usable on its own; unchecking keeps it\n"
            "strictly complementary.")
        self.csv_masks_check = QCheckBox("Include the mask columns")
        self.csv_masks_check.setChecked(True)
        self.csv_unc_check = QCheckBox("Include 1σ uncertainty columns")
        self.csv_unc_check.setChecked(True)
        self.csv_coeff_check = QCheckBox("Include cal slope / intercept columns")
        self.csv_coeff_check.setChecked(True)
        self.csv_coeff_check.setToolTip(
            "Also includes the cal-mean, pressure, temperature, and P/T factor "
            "columns used to check how each calibrated value was made.")
        self.csv_comment_check = QCheckBox("Put the provenance notes in the CSV as # lines")
        self.csv_comment_check.setChecked(False)
        self.csv_comment_check.setToolTip(
            "Off by default: neither Excel nor Igor skips a leading comment "
            "block without being told to.\nWith it off, no notes are written.")
        for widget in (self.csv_raw_check, self.csv_masks_check,
                       self.csv_unc_check, self.csv_coeff_check,
                       self.csv_comment_check):
            csv_layout.addWidget(widget)

        self.export_csv_button = QPushButton("Export CSV…")
        self.export_csv_button.setEnabled(False)
        self.export_csv_button.clicked.connect(self.on_export_csv_clicked)
        csv_row = QHBoxLayout()
        csv_row.addWidget(self.export_csv_button)
        csv_row.addStretch(1)
        csv_layout.addLayout(csv_row)
        vbox.addWidget(csv_box)

        # ---- ICARTT
        ict_box = QGroupBox("ICARTT (.ict)")
        ict_layout = QVBoxLayout(ict_box)
        ict_intro = QLabel(
            "The archive format (file format index 1001). Good ambient data "
            "only: every row this analysis blanked is written as the format's "
            "missing value, -99999, which is what that flag means — so no mask "
            "columns are written or needed. Time is seconds from midnight UTC."
        )
        ict_intro.setWordWrap(True)
        ict_layout.addWidget(ict_intro)

        self.icartt_sigma_check = QCheckBox("Include a 1σ uncertainty variable per gas")
        self.icartt_sigma_check.setChecked(True)
        self.icartt_drop_check = QCheckBox("Drop rows with no value for any gas")
        self.icartt_drop_check.setChecked(True)
        self.icartt_drop_check.setToolTip(
            "A row that is -99999 in every column carries nothing but a "
            "timestamp. Uncheck to keep the time base unbroken.")
        ict_layout.addWidget(self.icartt_sigma_check)
        ict_layout.addWidget(self.icartt_drop_check)

        meta_box = QGroupBox("Header metadata")
        meta_form = QFormLayout(meta_box)
        # QFormLayout defaults to FieldsStayAtSizeHint on macOS, which pins
        # every box to ~200 px and truncates the values regardless of how much
        # room the tab actually has. These fields hold addresses, PI lists and
        # sentences, so they get the full width instead.
        meta_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        meta_note = QLabel(
            "Saved in ucatsb_gui_config.yaml and reused for every flight — "
            "these are properties of the campaign, not of one dataset, so "
            "they are deliberately not written to a <dataset>_conf.yaml. "
            "That file is shared through the repository, so a fresh checkout "
            "arrives with the campaign already filled in. "
            "Edits here are not stored until Save defaults is pressed; "
            "exporting uses what is in the boxes either way."
        )
        meta_note.setWordWrap(True)
        meta_note.setStyleSheet(f"color: {MUTED_COLOR};")
        meta_form.addRow(meta_note)

        self.icartt_widgets = {}
        for key, label, height, tooltip in self.ICARTT_FIELDS:
            if height:
                widget = QPlainTextEdit()
                # Fixed, not a maximum: with only a maximum the form layout
                # hands back whatever it likes (the special-comments box asked
                # for 260 px and got 192, so the seeded text scrolled). The
                # whole pane is inside a QScrollArea, so a tall box costs
                # nothing but scrolling.
                widget.setFixedHeight(height)
                widget.textChanged.connect(self.on_icartt_meta_changed)
            else:
                widget = QLineEdit()
                widget.textChanged.connect(self.on_icartt_meta_changed)
            widget.setToolTip(tooltip)
            meta_form.addRow(f"{label}:", widget)
            self.icartt_widgets[key] = widget
        ict_layout.addWidget(meta_box)

        ict_row = QHBoxLayout()
        self.icartt_save_button = QPushButton(self.ICARTT_SAVE_LABEL)
        self.icartt_save_button.clicked.connect(self.on_save_icartt_meta)
        self.export_icartt_button = QPushButton("Export ICARTT…")
        self.export_icartt_button.setEnabled(False)
        self.export_icartt_button.clicked.connect(self.on_export_icartt_clicked)
        ict_row.addWidget(self.export_icartt_button)
        ict_row.addWidget(self.icartt_save_button)
        ict_row.addStretch(1)
        ict_layout.addLayout(ict_row)
        vbox.addWidget(ict_box)

        vbox.addStretch(1)
        self._apply_icartt_meta_to_controls(self.icartt_meta)
        return pane

    ICARTT_SAVE_LABEL = "Save defaults"

    def _apply_icartt_meta_to_controls(self, meta):
        """Populate the metadata form without marking it dirty -- the same
        _loading guard the per-gas controls use, for the same reason
        (setText fires textChanged immediately)."""
        was_loading = self._loading
        self._loading = True
        for key, widget in self.icartt_widgets.items():
            value = str(meta.get(key, "") or "")
            if isinstance(widget, QPlainTextEdit):
                widget.setPlainText(value)
            else:
                widget.setText(value)
        self._loading = was_loading
        self._update_icartt_save_state()

    def _icartt_meta_from_controls(self):
        """Every field must be listed in ICARTT_FIELDS: this rebuilds the dict
        wholesale from the widgets, so a key with no widget would be dropped
        on the next edit -- the same trap _controls_to_settings() has."""
        return {
            key: (widget.toPlainText() if isinstance(widget, QPlainTextEdit)
                  else widget.text()).strip()
            for key, widget in self.icartt_widgets.items()
        }

    def _icartt_is_dirty(self):
        """A comparison against what was last saved, not a flag -- same
        reasoning as _is_dirty(), and kept separate from it because the two
        Save buttons write different files. Folding this into _current_state()
        would make the per-dataset Save claim to cover app-level metadata it
        does not write."""
        return self._icartt_meta_from_controls() != self._saved_icartt_meta

    def _update_icartt_save_state(self):
        if not hasattr(self, "icartt_save_button"):
            return
        dirty = self._icartt_is_dirty()
        self.icartt_save_button.setText(self.ICARTT_SAVE_LABEL + (" •" if dirty else ""))
        self.icartt_save_button.setToolTip(
            ("Unsaved changes. " if dirty else "No unsaved changes. ")
            + f"Save writes these fields to {self.default_config_path.name} "
              "as the defaults for every flight. Exporting uses whatever is "
              "in the boxes, saved or not.")

    def on_icartt_meta_changed(self):
        if self._loading or self._initializing:
            return
        self._update_icartt_save_state()

    def on_save_icartt_meta(self):
        """Write the metadata to the shared config as the new defaults. The
        only writer of that file, which is what keeps a tracked file from
        churning every time a dataset is opened."""
        self.icartt_meta = self._icartt_meta_from_controls()
        self._save_shared_config()
        self._saved_icartt_meta = copy.deepcopy(self.icartt_meta)
        self._update_icartt_save_state()
        self._flash_button(self.icartt_save_button,
                           f"Saved to {self.default_config_path.name}",
                           self._update_icartt_save_state)
        return True

    def _apply_settings_to_controls(self, settings: dict):
        """Populate the controls from a per-gas settings dict without
        triggering on_control_changed (and re-saving/redrawing) per field."""
        self._loading = True
        self.warmup_spin.setValue(settings["warmup_min"])
        self.end_flight_spin.setValue(settings.get("end_flight_min", 0))
        self.pumps_check.setChecked(bool(settings.get("require_pumps", False)))
        self.pressure_tol_spin.setValue(settings["pressure_tol_mbar"])
        self.pressure_correct_check.setChecked(
            bool(settings.get("pressure_correct", False)))
        self.pressure_smooth_spin.setValue(settings.get("pressure_smooth_s", 0))
        self.temperature_correct_check.setChecked(
            bool(settings.get("temperature_correct", False)))
        self.flag_air_spin.setValue(settings["flag_air_s"])
        self.cal1_start_spin.setValue(settings["cal1_window_s"][0])
        self.cal1_end_spin.setValue(settings["cal1_window_s"][1])
        self.cal2_start_spin.setValue(settings["cal2_window_s"][0])
        self.cal2_end_spin.setValue(settings["cal2_window_s"][1])
        self.drift_combo.setCurrentText(settings["drift_model"])
        self.smooth_spin.setValue(settings["drift_smooth_events"])
        self.fixed_slope_spin.setValue(settings.get("fixed_slope", 0.0))
        self.drift_model = settings["drift_model"]
        self.drift_smooth_events = settings["drift_smooth_events"]
        self.fixed_slope = settings.get("fixed_slope", 0.0)
        self._show_drift_extra(settings["drift_model"])
        self._loading = False

    def _controls_to_settings(self) -> dict:
        """Every persisted setting must be listed here: on_control_changed
        assigns this dict wholesale, so anything missing is silently dropped
        from the config on the next control change."""
        return {
            "warmup_min": self.warmup_spin.value(),
            "end_flight_min": self.end_flight_spin.value(),
            "require_pumps": self.pumps_check.isChecked(),
            "pressure_tol_mbar": self.pressure_tol_spin.value(),
            "pressure_correct": self.pressure_correct_check.isChecked(),
            "pressure_smooth_s": self.pressure_smooth_spin.value(),
            "temperature_correct": self.temperature_correct_check.isChecked(),
            "flag_air_s": self.flag_air_spin.value(),
            "cal1_window_s": [self.cal1_start_spin.value(), self.cal1_end_spin.value()],
            "cal2_window_s": [self.cal2_start_spin.value(), self.cal2_end_spin.value()],
            "drift_model": self.drift_combo.currentText(),
            "drift_smooth_events": self.smooth_spin.value(),
            "fixed_slope": self.fixed_slope_spin.value(),
        }

    def _show_drift_extra(self, model: str):
        """Show whichever number the chosen drift model needs, and only that.

        The two spin boxes share one slot in the row (see where they are
        built), so this is what keeps exactly one of them on screen. `smooth`
        stays visible-but-disabled for the models that need no number at all,
        rather than leaving the slot empty and letting the combo stretch into
        it on every model change.
        """
        pinned = model == "fixed slope"
        self.fixed_slope_spin.setVisible(pinned)
        self.fixed_slope_reset_button.setVisible(pinned)
        self.fixed_slope_reset_button.setEnabled(pinned)
        self.smooth_spin.setVisible(not pinned)
        self.smooth_spin.setEnabled(model == "smooth")

    def _select_gas(self, gas: str):
        """Sync gas-dependent controls (masking/cal boxes, per-gas
        settings) to `gas`. Shared by on_gas_changed and load_csv so a
        freshly loaded file ends up in the same state as if the user had
        picked this gas from the combo box themselves."""
        self.current_gas = gas
        has_masking = GASES[gas].get("has_masking", True)
        self.mask_box.setEnabled(has_masking)
        self.cal_box.setEnabled(has_masking)
        # The overlay toggle moved to the toolbar, so it no longer greys out
        # with the box it used to sit in -- it has to be told. A gas with no
        # cal bottles has nothing to overlay. Its checked state is left alone,
        # exactly as a disabled checkbox's would be, so switching away and back
        # returns it as it was.
        self.main_pane.calibrated_action.setEnabled(has_masking)
        # flag_box is deliberately absent from that list -- see where it is
        # built. Its readout is per gas, so it does have to follow along.
        # The pressure correction needs THIS gas's detector pressure column,
        # and which detector that is varies by gas (and by flight -- the Jul
        # 2026 file's d2 is a different instrument from the Feb 2025 one), so
        # it is checked per gas rather than once at load like `Pumps on`.
        # Set after mask_box's setEnabled, which Qt would otherwise override.
        # Same condition for the smoothing window: it exists to feed that
        # correction (and to draw the trace it divides by), so a gas with no
        # detector of its own has nothing for it to smooth.
        has_detector = self._pressure_column(gas) is not None
        self.pressure_correct_check.setEnabled(has_masking and has_detector)
        self.pressure_smooth_spin.setEnabled(has_masking and has_detector)
        # Checked separately: a flight's schema can carry one of the two
        # columns and not the other, and the Feb 2025 d2 is a different
        # instrument from the Jul 2026 one.
        self.temperature_correct_check.setEnabled(
            has_masking and self._temperature_column(gas) is not None)
        if has_masking:
            self._apply_settings_to_controls(self.config[gas])
        self._update_flag_readout()

    def on_gas_changed(self, new_gas: str):
        if not new_gas:
            return
        self._select_gas(new_gas)
        if self._initializing:
            return
        self.refresh()

    def on_show_notes_toggled(self, checked):
        """Show or hide the Timeseries figure's grey note block.

        Sets the artist's visibility and redraws the canvas, nothing more --
        deliberately not a refresh(). Nothing about the data, the masks or the
        calibration changes, the user is typically zoomed in on something, and
        recomputing five gases' analyses to hide a caption would be absurd.
        The same reasoning as Correlations' "Hide flagged points".
        """
        self.show_plot_notes = checked
        if self._initializing:
            return
        if self._notes_artist is not None:
            self._notes_artist.set_visible(checked)
            self.canvas.draw_idle()

    def on_control_changed(self):
        if self._loading or self._initializing or self.current_gas is None:
            return
        # Seeding runs BEFORE the settings are read back, so the value it puts
        # in the box is the one that gets stored. It only fires on the switch
        # into "fixed slope" with the box still on "auto", so a saved config's
        # own number is never overwritten and neither is a hand-typed one.
        if (self.drift_combo.currentText() == "fixed slope"
                and self.drift_model != "fixed slope"
                and not self.fixed_slope_spin.value()):
            self._seed_fixed_slope()
        settings = self._controls_to_settings()
        self.config[self.current_gas] = settings
        self.drift_model = settings["drift_model"]
        self.drift_smooth_events = settings["drift_smooth_events"]
        self.fixed_slope = settings["fixed_slope"]
        self._show_drift_extra(self.drift_model)
        if hasattr(self, "corr_cal_target") and self._corr_cal_gas() == self.current_gas:
            self._sync_corr_cal_controls()
        self._mark_dirty()
        self.refresh(preserve_view=True)

    def _seed_fixed_slope(self, gas_key=None, spin=None):
        """Put the constant model's slope in the box when "fixed slope" is
        first chosen, so the user starts from the flight's own gain instead of
        an empty field.

        That slope is exactly `1 / span_gain`: the constant model replaces each
        bottle's nodes with its flight mean, so its gain is
        `(A_hi - A_lo) / (mean_hi - mean_lo)`. Taking it from `span_gain` rather
        than recomputing means it comes from the calibration the user is
        looking at -- and `span_gain` is a property of the cal-point means
        alone, so the currently cached result gives the same answer whatever
        model produced it.

        Leaves the box at "auto" when there is no usable span (one bottle, or
        no calibration at all); calibrate_series falls back the same way.
        """
        gas_key = gas_key or self.current_gas
        spin = spin or self.fixed_slope_spin
        result = self._calibration_for(gas_key) or {}
        span_gain = result.get("span_gain") if result.get("ok") else None
        if not span_gain:
            return False
        spin.blockSignals(True)
        spin.setValue(1.0 / span_gain)
        spin.blockSignals(False)
        return True

    def on_fixed_slope_reset(self):
        if self._loading or self._initializing or self.current_gas is None:
            return
        if self._seed_fixed_slope():
            settings = self._controls_to_settings()
            self.config[self.current_gas] = settings
            self.fixed_slope = settings["fixed_slope"]
            if hasattr(self, "corr_cal_target") and self._corr_cal_gas() == self.current_gas:
                self._sync_corr_cal_controls()
            self._mark_dirty()
            self.refresh(preserve_view=True)

    def _current_state(self):
        """Everything a config file holds, as comparable plain data."""
        return copy.deepcopy({"config": self.config, "cals": self.cal_selection,
                              "flagged": self.flagged})

    def _snapshot_state(self):
        """Take the current settings as "saved" -- called after a load or a
        successful save, and nowhere else."""
        self._saved_state = self._current_state()
        self._update_save_state()

    def _is_dirty(self):
        """Settings differ from the last load/save.

        A comparison rather than a flag set by every handler: a flag has to be
        cleared in as many places as it is set, and gets it wrong when a
        change is undone. Editing a spin box back to where it started leaves
        nothing to save, and this says so.
        """
        return self._saved_state is not None and self._current_state() != self._saved_state

    def _update_file_labels(self):
        """Name the dataset *and* the configuration it is being viewed through,
        on both control pages.

        Which config is loaded was otherwise invisible outside the Cal Tanks
        tab, and since Save names a variant a flight can easily have several —
        "which one am I looking at" is a question the panel should answer
        without a tab change. Both pages are written from here for the same
        reason `_rebuild_recent_menus` builds both menus: two writers of one
        fact drift.

        The config name is deliberately the **filename**, not the variant
        recovered from it: the variant is a name for a thing on disk, and when
        the two could disagree — an older config named freely, one opened from
        another directory — the filename is the one that tells you what to go
        and look at.

        **Always two lines, never wrapped.** One row per fact, so the block has
        a fixed height and the second line means the same thing every time
        rather than sometimes being the tail of the first.

        Which makes the values **elided, not clipped**: a config name runs 305
        px for the default and 349 for a short variant, against ~300 px of
        panel, so without wrapping the end of the name simply vanishes — the
        panel-width failure mode, on the one label whose whole job is to name
        a file. `_elide_field` trims it to fit and the tooltip keeps both full
        paths.
        """
        if not hasattr(self, "file_label"):
            return
        if self.csv_path is None:
            name, conf, tip = "No file loaded", "—", ""
        else:
            name = self.csv_path.name
            if self.config_loaded_from is not None:
                conf = self.config_loaded_from.name
                conf_tip = f"Configuration: {self.config_loaded_from}"
            else:
                # Not the same as "no settings": the dataset is on
                # DEFAULT_GAS_SETTINGS, and config_path is where the first Save
                # would put them.
                conf = "defaults"
                conf_tip = ("No configuration file — built-in defaults. "
                            f"Save… writes {self.config_path.name}.")
            tip = f"{self.csv_path}\n{conf_tip}"
        # Rich text: the field names are bold, the names themselves ordinary
        # weight, so the eye lands on the two values rather than on the labels
        # that never change. escape() because a file name is not our text --
        # an `&` in one would otherwise be read as markup and vanish.
        text = (f"<b>Filename:</b> {html.escape(self._elide_field('Filename:', name))}<br>"
                f"<b>Config:</b> {html.escape(self._elide_field('Config:', conf))}")
        for label in (self.file_label, self.corr_file_label):
            label.setText(text)
            label.setToolTip(tip)

    def _elide_field(self, field, value):
        """Trim `value` to what is left of the panel beside a bold `field:`.

        Elided from the **left**, unusually: the config name begins with the
        dataset stem, which the Filename line directly above already shows, so
        the front is the one part guaranteed to be redundant while the tail
        (the variant, and `_conf.yaml`) is what distinguishes one config from
        another. Right-eliding would cut exactly the part worth reading.

        The budget is `CONTROLS_WIDTH` minus the layout margins rather than the
        label's own width, because the label is measured before it is laid out
        (`_update_file_labels` runs during load) and the panel is fixed-width,
        so the static figure is the true one at every moment.
        """
        bold = QFont(self.file_label.font())
        bold.setBold(True)
        used = QFontMetrics(bold).horizontalAdvance(field + " ")
        room = max(40, CONTROLS_WIDTH - 2 * CONTROLS_MARGIN - used)
        return QFontMetrics(self.file_label.font()).elidedText(
            value, Qt.ElideLeft, room)

    def _update_save_state(self):
        """Mark the Save button when there is something to save."""
        if not hasattr(self, "save_button"):
            return
        dirty = self._is_dirty()
        self.save_button.setText(self.SAVE_LABEL + (" •" if dirty else ""))
        name = self.config_path.name if self.config_path else "a new file"
        self.save_button.setToolTip(
            ("Unsaved changes. " if dirty else "No unsaved changes. ")
            + f"Save writes a config file beside the CSV — {name} unless you "
              "name a variant, which saves a second configuration of the same "
              "dataset."
        )

    def _mark_dirty(self):
        """Settings changed. Nothing is written -- this only updates the
        Save button. Replaces the old auto-save so that opening a saved
        analysis, experimenting and quitting leaves the file untouched."""
        self._update_save_state()

    def _choose_config_name(self):
        """Ask what to call this configuration, and return the full path.

        A variant name, not a filename: the dataset stem and the `_conf`
        suffix are added by `flight_config_path`, so every config this writes
        is one `_config_candidates` will find when the flight is next opened.
        The old QFileDialog let the name be typed whole, which meant a config
        saved as `test.yaml` was silently invisible from then on -- a save
        dialog that cannot express the loader's one requirement.

        The cost is that a config can no longer be written to another
        directory. Nothing is lost: such a file was already unreachable by the
        stem search, and "Load configuration…" still opens a config by name
        from anywhere.
        """
        dialog = QDialog(self)
        dialog.setWindowTitle("Save configuration")
        dialog.setMinimumWidth(420)
        layout = QVBoxLayout(dialog)

        label = QLabel(f"Saved beside {self.csv_path.name}. Leave the name "
                       f"empty for this flight's main configuration, or name a "
                       f"variant to keep several.")
        label.setWordWrap(True)
        layout.addWidget(label)

        form = QFormLayout()
        edit = QLineEdit(config_variant_name(self.csv_path, self.config_path))
        edit.setPlaceholderText("(none — the main configuration)")
        form.addRow("Variant name:", edit)
        layout.addLayout(form)

        # The composed name, live. The whole point of taking a variant instead
        # of a filename is that the user no longer types the parts that matter,
        # so the result has to be visible rather than a surprise on disk.
        preview = QLabel()
        preview.setWordWrap(True)
        preview.setStyleSheet(f"color: {MUTED_COLOR};")
        layout.addWidget(preview)

        def update_preview():
            path = flight_config_path(self.csv_path, edit.text())
            note = "  (replaces the existing file)" if path.exists() else ""
            preview.setText(f"Saves as {path.name}{note}")

        edit.textChanged.connect(update_preview)
        update_preview()

        # Existing configs, so re-saving over one is a click rather than
        # remembering how it was spelled.
        existing = self._config_candidates(self.csv_path)
        if existing:
            layout.addWidget(QLabel("Existing configurations:"))
            listing = QListWidget()
            listing.setMaximumHeight(96)
            for path in existing:
                item = QListWidgetItem(path.name)
                item.setData(Qt.UserRole, config_variant_name(self.csv_path, path))
                listing.addItem(item)
            listing.itemClicked.connect(
                lambda item: edit.setText(item.data(Qt.UserRole) or ""))
            layout.addWidget(listing)

        box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        box.accepted.connect(dialog.accept)
        box.rejected.connect(dialog.reject)
        layout.addWidget(box)

        # Looped, because the overwrite prompt QFileDialog gave for free is now
        # ours to ask: answering "no" has to come back to the name, not cancel
        # the save.
        while True:
            edit.setFocus()
            edit.selectAll()
            if dialog.exec_() != QDialog.Accepted:
                return None
            path = flight_config_path(self.csv_path, edit.text())
            if not path.exists():
                return path
            answer = QMessageBox.question(
                self, "Replace configuration",
                f"{path.name} already exists.\n\nReplace it?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if answer == QMessageBox.Yes:
                return path

    def on_save_clicked(self):
        """Write the current settings to a config file of the user's choosing.

        Save-as every time, deliberately: the request was to keep several
        configurations per dataset, so the name is always offered for editing
        rather than silently overwriting whatever was opened.
        """
        if self.df is None:
            return
        path = self._choose_config_name()
        if path is None:
            return False
        try:
            save_config(path, self.config, cal_selection=self.cal_selection,
                        flagged=flagged_to_yaml(self.flagged, len(self.raw_df)))
        except OSError as e:
            QMessageBox.warning(self, "Save configuration",
                                f"Could not write {path.name}:\n{e}")
            return False
        self.config_path = path
        self.config_loaded_from = path
        self._snapshot_state()
        self._update_file_labels()
        self._update_tank_readout()
        return True

    def on_load_config_clicked(self):
        """Open a config file by name -- for one saved somewhere the
        dataset-stem search would not find it, or to switch configurations
        without reloading the CSV."""
        if self.df is None:
            return
        if not self._confirm_discard("opening another configuration"):
            return
        start = str(self.config_path or self.csv_path.parent)
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load configuration", start, "YAML Files (*.yaml);;All Files (*)")
        if not path_str:
            return
        self._apply_config_file(Path(path_str))
        self._select_gas(self.current_gas)
        self._sync_corr_cal_controls()
        self._update_tank_readout()
        self.refresh()

    def _confirm_discard(self, action_label):
        """Ask before losing unsaved settings. True = go ahead.

        Only ever asked when something actually changed, so the prompt stays
        meaningful; "Don't Save" is the stated goal (quit without disturbing
        the state you started from), and Cancel aborts whatever triggered it.
        """
        if not self._is_dirty():
            return True
        name = self.config_path.name if self.config_path else "this configuration"
        answer = QMessageBox.question(
            self, "Unsaved changes",
            f"Settings have changed since {name} was opened.\n\n"
            f"Save them before {action_label}?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
            QMessageBox.Save)
        if answer == QMessageBox.Cancel:
            return False
        if answer == QMessageBox.Save:
            return bool(self.on_save_clicked())
        return True

    def _confirm_discard_icartt(self):
        """Ask before losing unsaved ICARTT metadata. Asked separately from
        _confirm_discard because the two Save buttons write different files:
        one prompt offering to "save" would have to pick one of them, and
        would silently not write the other."""
        if not self._icartt_is_dirty():
            return True
        answer = QMessageBox.question(
            self, "Unsaved ICARTT metadata",
            f"The ICARTT header metadata has changed since it was last saved "
            f"to {self.default_config_path.name}.\n\nSave it before quitting?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
            QMessageBox.Save)
        if answer == QMessageBox.Cancel:
            return False
        if answer == QMessageBox.Save:
            return bool(self.on_save_icartt_meta())
        return True

    def closeEvent(self, event):
        if self._confirm_discard("quitting") and self._confirm_discard_icartt():
            event.accept()
        else:
            event.ignore()

    SAVE_LABEL = "Save…"

    # What the button copies: the masking values and both cal mean windows.
    # The drift model and its two numbers -- the smoothing window and the fixed
    # slope -- are the per-gas settings left out. They are a judgement about
    # that gas's cal record (how noisy its injections are, what gain its tanks
    # ask for), not a description of the flight, and the fixed slope is in the
    # gas's own units besides: copying CO2's 1.0367 onto CH4 would be
    # meaningless.
    COPIED_SETTING_KEYS = ("warmup_min", "end_flight_min", "require_pumps",
                           "pressure_tol_mbar", "pressure_correct",
                           "pressure_smooth_s", "temperature_correct",
                           "flag_air_s",
                           "cal1_window_s", "cal2_window_s")
    COPY_SETTINGS_LABEL = "Copy settings to all gases"

    def on_copy_masking_to_all(self):
        """Apply this gas's masking + cal-window settings to every other
        calibrated gas.

        Ozone is excluded by construction rather than by name -- it has no
        entry in self.config at all (has_masking=False), so "every gas in the
        config" is already the right set and stays right if a gas is added.
        """
        if self.current_gas is None:
            return
        # Read from the live controls, not from self.config: a spin box being
        # edited by keyboard commits its value on focus-out, which the click
        # on this button is, and the ordering of that against valueChanged is
        # not something to bet the copied numbers on.
        live = self._controls_to_settings()
        self.config[self.current_gas] = live
        targets = [gas for gas in self.config if gas != self.current_gas]
        for gas in targets:
            # deepcopy per target, not one shared source dict: cal1_window_s /
            # cal2_window_s are lists, and handing every gas the same list
            # object makes yaml.safe_dump emit anchors (&id001/*id001) into
            # the conf file, and makes a later in-place edit of one gas's
            # window silently change the others.
            self.config[gas].update(
                {key: copy.deepcopy(live[key]) for key in self.COPIED_SETTING_KEYS})
        self._mark_dirty()
        # No refresh: the current gas's own settings are unchanged, so the
        # plots are still correct. The other gases redraw when selected.
        self._flash_button(self.copy_mask_button,
                           f"Copied to {', '.join(targets)}" if targets else "Nothing to copy",
                           self.COPY_SETTINGS_LABEL)

    def on_corr_gas_changed(self, axis, gas_key):
        if self._loading or self._initializing or not gas_key:
            return
        if axis == "x":
            self.corr_x_gas = gas_key
        else:
            self.corr_y_gas = gas_key
        self._populate_corr_flag_target()
        self._populate_corr_cal_target()
        # A different tracer is a different set of numbers on that axis, so
        # the old limits mean nothing -- rescale, as a gas change does on the
        # timeseries.
        self._refresh_corr(preserve_view=False)

    def _populate_corr_cal_target(self):
        """Rebuild the settings target from the current X/Y tracers.

        These controls edit saved per-gas calibration settings, so there is no
        "both" entry: a drift model or P/T correction belongs to one tracer's
        calibration at a time.
        """
        if not hasattr(self, "corr_cal_target"):
            return
        previous = self.corr_cal_target.currentData()
        entries = []
        if self.corr_y_gas:
            entries.append(("y", f"{self.corr_y_gas} (Y)"))
        if self.corr_x_gas and self.corr_x_gas != self.corr_y_gas:
            entries.append(("x", f"{self.corr_x_gas} (X)"))
        loading = self._loading
        self._loading = True
        self.corr_cal_target.clear()
        for role, label in entries:
            self.corr_cal_target.addItem(label, role)
        index = self.corr_cal_target.findData(previous)
        self.corr_cal_target.setCurrentIndex(index if index >= 0 else 0)
        self._loading = loading
        self._sync_corr_cal_controls()

    def _corr_cal_gas(self):
        role = self.corr_cal_target.currentData() or "y"
        return self.corr_x_gas if role == "x" else self.corr_y_gas

    def on_corr_cal_target_changed(self, _index):
        if self._loading or self._initializing:
            return
        self._sync_corr_cal_controls()

    def _sync_corr_cal_controls(self):
        gas = self._corr_cal_gas()
        has_settings = bool(gas and GASES[gas].get("has_masking", True))
        self.corr_cal_box.setEnabled(True)
        if not has_settings:
            for widget in (
                    self.corr_pressure_correct_check,
                    self.corr_pressure_smooth_spin,
                    self.corr_temperature_correct_check,
                    self.corr_drift_combo,
                    self.corr_smooth_spin,
                    self.corr_fixed_slope_spin,
                    self.corr_fixed_slope_reset_button):
                widget.setEnabled(False)
            return
        settings = self.config.get(gas, DEFAULT_GAS_SETTINGS)
        loading = self._loading
        self._loading = True
        self.corr_pressure_correct_check.setChecked(
            bool(settings.get("pressure_correct", False)))
        self.corr_pressure_smooth_spin.setValue(settings.get("pressure_smooth_s", 0))
        self.corr_temperature_correct_check.setChecked(
            bool(settings.get("temperature_correct", False)))
        self.corr_drift_combo.setCurrentText(settings["drift_model"])
        self.corr_smooth_spin.setValue(settings["drift_smooth_events"])
        self.corr_fixed_slope_spin.setValue(settings.get("fixed_slope", 0.0))
        self._show_corr_drift_extra(settings["drift_model"])
        self._loading = loading
        self._update_corr_cal_enabled(gas)

    def _update_corr_cal_enabled(self, gas):
        has_detector = self._pressure_column(gas) is not None
        self.corr_pressure_correct_check.setEnabled(has_detector)
        self.corr_pressure_smooth_spin.setEnabled(has_detector)
        self.corr_temperature_correct_check.setEnabled(
            self._temperature_column(gas) is not None)
        self.corr_drift_combo.setEnabled(True)
        self._show_corr_drift_extra(self.corr_drift_combo.currentText())

    def _show_corr_drift_extra(self, model):
        pinned = model == "fixed slope"
        self.corr_fixed_slope_spin.setVisible(pinned)
        self.corr_fixed_slope_spin.setEnabled(pinned)
        self.corr_fixed_slope_reset_button.setVisible(pinned)
        self.corr_fixed_slope_reset_button.setEnabled(pinned)
        self.corr_smooth_spin.setVisible(not pinned)
        self.corr_smooth_spin.setEnabled(model == "smooth")

    def on_corr_fixed_slope_reset(self):
        if self._loading or self._initializing:
            return
        gas = self._corr_cal_gas()
        if not gas or not GASES[gas].get("has_masking", True):
            return
        if not self._seed_fixed_slope(gas_key=gas, spin=self.corr_fixed_slope_spin):
            return
        settings = copy.deepcopy(self.config[gas])
        settings["fixed_slope"] = self.corr_fixed_slope_spin.value()
        self.config[gas] = settings
        if gas == self.current_gas:
            self._apply_settings_to_controls(settings)
        self._mark_dirty()
        self.refresh(preserve_view=True)

    def on_corr_cal_control_changed(self):
        if self._loading or self._initializing:
            return
        gas = self._corr_cal_gas()
        if not gas or not GASES[gas].get("has_masking", True):
            return
        if (self.corr_drift_combo.currentText() == "fixed slope"
                and self.config[gas].get("drift_model") != "fixed slope"
                and not self.corr_fixed_slope_spin.value()):
            self._seed_fixed_slope(gas_key=gas, spin=self.corr_fixed_slope_spin)
        settings = copy.deepcopy(self.config[gas])
        settings.update({
            "pressure_correct": self.corr_pressure_correct_check.isChecked(),
            "pressure_smooth_s": self.corr_pressure_smooth_spin.value(),
            "temperature_correct": self.corr_temperature_correct_check.isChecked(),
            "drift_model": self.corr_drift_combo.currentText(),
            "drift_smooth_events": self.corr_smooth_spin.value(),
            "fixed_slope": self.corr_fixed_slope_spin.value(),
        })
        self.config[gas] = settings
        self._show_corr_drift_extra(settings["drift_model"])
        if gas == self.current_gas:
            self._apply_settings_to_controls(settings)
        self._mark_dirty()
        self.refresh(preserve_view=True)

    def _populate_corr_flag_target(self):
        """Rebuild the "Flag applies to" choices for the current pair.

        Keeps the previous *role* (Y / X / both) rather than the previous gas
        name: after swapping the axes, "the Y tracer" is still what the user
        meant, and re-resolving to a gas would silently retarget the tool.
        """
        previous = self.corr_flag_target.currentData()
        x_gas, y_gas = self.corr_x_gas, self.corr_y_gas
        if not x_gas or not y_gas:
            return
        entries = [("y", f"{y_gas}  (Y axis)"), ("x", f"{x_gas}  (X axis)")]
        if x_gas != y_gas:
            entries.append(("both", "both tracers"))
        self._loading, loading = True, self._loading
        self.corr_flag_target.clear()
        for role, label in entries:
            self.corr_flag_target.addItem(label, role)
        index = self.corr_flag_target.findData(previous)
        self.corr_flag_target.setCurrentIndex(index if index >= 0 else 0)
        self._loading = loading
        self._update_corr_flag_readout()

    def _corr_flag_gases(self):
        """Which gases a flag drawn on the scatter applies to."""
        role = self.corr_flag_target.currentData() or "y"
        if role == "x":
            return [self.corr_x_gas]
        if role == "both":
            return list(dict.fromkeys([self.corr_y_gas, self.corr_x_gas]))
        return [self.corr_y_gas]

    def on_corr_flag_target_changed(self, _index):
        if self._loading or self._initializing:
            return
        self._update_corr_flag_readout()

    def _update_corr_flag_readout(self):
        """Counts for the tracer(s) the tool is currently pointed at."""
        gases = [g for g in self._corr_flag_gases() if g]
        parts, total = [], 0
        for gas in gases:
            ranges = self.flagged.get(gas, [])
            total += ranges_row_count(ranges)
            if ranges:
                parts.append(f"{gas}: {ranges_row_count(ranges)} in "
                             f"{len(ranges)} region{'' if len(ranges) == 1 else 's'}")
        self.corr_flag_label.setText(
            "No points flagged" if not parts else "  ·  ".join(parts))
        self.corr_flag_clear_button.setEnabled(bool(total))
        self.corr_flag_undo_button.setEnabled(bool(self._flag_undo))
        # Keyed off what is DRAWN, not off the combo's scope: the markers show
        # flags on either axis, so the toggle stays useful while the combo
        # points at a gas with none of its own.
        self.corr_hide_check.setEnabled(self._corr_flag_scatter is not None)

    def on_corr_hide_flagged(self, checked):
        """Show or hide the struck-out markers, without replotting.

        Deliberately not a refresh(): the whole point is that the current zoom
        survives, and the flags themselves are unchanged -- this is a view
        toggle, so nothing is dirtied and nothing is recomputed. Only the
        artist's visibility and the Home target move.
        """
        if self._initializing:
            return
        self.corr_hide_flagged = checked
        if self._corr_flag_scatter is not None:
            self._corr_flag_scatter.set_visible(not checked)
        self._retarget_corr_home()
        self.corr_pane.canvas.draw_idle()

    def _corr_home_limits(self):
        """(xlim, ylim) Home should return to, given what is currently drawn.

        Built from the recorded plotted data rather than from `ax.dataLim`,
        which still carries the hidden markers -- and which the selectors have
        their own history of polluting. 5% margins, matching what matplotlib's
        own autoscale would have produced. None when there is nothing to frame.
        """
        plotted = self._corr_plotted
        if not plotted:
            return None
        keep = plotted["keep"]
        if self.corr_hide_flagged:
            # Only the paired, unflagged record: `keep` already excludes every
            # flagged row, which is exactly what stays on screen.
            shown = keep
        else:
            shown = keep | plotted["flagged"]
        fx, fy = plotted["x"][shown], plotted["y"][shown]
        if not len(fx.dropna()) or not len(fy.dropna()):
            return None

        def span(values):
            lo, hi = float(values.min()), float(values.max())
            pad = 0.05 * (hi - lo) if hi > lo else (abs(hi) * 0.05 or 0.5)
            return lo - pad, hi + pad

        return span(fx), span(fy)

    def _retarget_corr_home(self):
        """Point Home at what is currently drawn -- or hand it back to the nav
        stack when the markers are showing, since then the stock full-scale
        view is already right."""
        if not self.corr_hide_flagged:
            self.corr_pane.clear_home_view()
            return
        limits = self._corr_home_limits()
        if limits and self._corr_ax is not None:
            self.corr_pane.set_home_view(self._corr_ax, *limits)

    def on_corr_flag_clear(self):
        """Clear the flags on whichever tracer the combo names -- the same
        scope the tool writes to, so the button undoes what it does."""
        gases = [g for g in self._corr_flag_gases() if self.flagged.get(g)]
        if not gases:
            return
        self._flag_undo.append(copy.deepcopy(self.flagged))
        for gas in gases:
            self.flagged.pop(gas, None)
        self._after_flag_change()

    def on_corr_tooltip_press(self, event):
        """Persistent readout for the nearest displayed correlation point."""
        if event.button != 3 or event.inaxes is not self._corr_ax:
            return
        text = self._corr_tooltip_text(event)
        if text:
            self._show_corr_tooltip_popup(text, event)

    def _show_corr_tooltip_popup(self, text, event):
        self._close_corr_tooltip_popup()

        popup = QWidget(self, Qt.Tool | Qt.FramelessWindowHint)
        popup.setAttribute(Qt.WA_DeleteOnClose)
        popup.setStyleSheet(
            "QWidget { background: #ffffff; color: #222222; "
            "border: 1px solid #8a8a8a; }"
            "QPushButton { border: none; padding: 1px 6px; font-weight: bold; }"
            "QPushButton:hover { background: #e8e8e8; }"
        )
        layout = QVBoxLayout(popup)
        layout.setContentsMargins(8, 4, 8, 8)
        layout.setSpacing(4)

        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.addStretch(1)
        close_button = QPushButton("x")
        close_button.setFixedSize(20, 20)
        close_button.setToolTip("Close")
        top.addWidget(close_button)
        layout.addLayout(top)

        label = QLabel(text)
        label.setTextFormat(Qt.PlainText)
        label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        label.setStyleSheet("font-family: Menlo, monospace; border: none;")
        layout.addWidget(label)

        close_button.clicked.connect(self._close_corr_tooltip_popup)
        popup.destroyed.connect(
            lambda _=None, p=popup: setattr(
                self, "corr_tooltip_popup", None
            ) if self.corr_tooltip_popup is p else None)

        self.corr_tooltip_popup = popup
        popup.adjustSize()
        qevent = getattr(event, "guiEvent", None)
        if qevent is not None and hasattr(qevent, "globalPos"):
            pos = qevent.globalPos()
        else:
            pos = self.corr_pane.canvas.mapToGlobal(
                self.corr_pane.canvas.rect().center())
        popup.move(pos.x() + 12, pos.y() + 12)
        popup.show()

    def _close_corr_tooltip_popup(self):
        popup = self.corr_tooltip_popup
        self.corr_tooltip_popup = None
        if popup is not None:
            popup.close()

    def _corr_tooltip_text(self, event):
        plotted = self._corr_plotted
        if not plotted or event.x is None or event.y is None:
            return None
        mask = plotted["keep"].copy()
        if not self.corr_hide_flagged:
            mask |= plotted["flagged"]
        if not mask.any():
            return None

        x_vals = plotted["x"][mask]
        y_vals = plotted["y"][mask]
        points = self._corr_ax.transData.transform(
            list(zip(x_vals.to_numpy(), y_vals.to_numpy())))
        dx = points[:, 0] - event.x
        dy = points[:, 1] - event.y
        dist2 = dx * dx + dy * dy
        nearest_pos = int(dist2.argmin())
        max_px = max(10.0, self.corr_marker_size + 8.0)
        if dist2[nearest_pos] > max_px * max_px:
            return None

        row = x_vals.index[nearest_pos]
        lines = [f"row {row + self.presync_dropped}   {self.df['datetime'].loc[row]}"]
        for axis, gas_key, values, sigma, unit in (
                ("X", plotted["x_gas"], plotted["x"], plotted["x_sigma"], plotted["x_unit"]),
                ("Y", plotted["y_gas"], plotted["y"], plotted["y_sigma"], plotted["y_unit"])):
            lines.extend(self._corr_tooltip_gas_lines(
                axis, gas_key, row, values, sigma, unit))
        return "\n".join(lines)

    def _corr_tooltip_gas_lines(self, axis, gas_key, row, values, sigma, unit):
        value = values.loc[row]
        sigma_value = None if sigma is None else sigma.loc[row]
        head = f"{axis} {gas_key}: {self._fmt_value(value)}"
        if sigma_value is not None and not pd.isna(sigma_value):
            head += f" ± {float(sigma_value):.3f}"
        head += f" {unit}".rstrip()

        info = GASES[gas_key]
        if not info.get("has_masking", True):
            return [head, "  P: n/a", "  T: n/a", "  factor: n/a"]

        analysis = self._analysis_for(gas_key)
        p_raw, p_corr = self._corr_tooltip_pressure(gas_key, row, analysis)
        t_raw, t_corr = self._corr_tooltip_temperature(gas_key, row, analysis)
        factor = analysis["correction_factor"]
        factor_value = 1.0 if factor is None else factor.loc[row]
        return [
            head,
            f"  P: raw {p_raw}; used {p_corr}",
            f"  T: raw {t_raw}; ref {t_corr}",
            f"  factor: {self._fmt_value(factor_value)}",
        ]

    def _corr_tooltip_pressure(self, gas_key, row, analysis):
        col = self._pressure_column(gas_key)
        if col is None:
            return "n/a", "n/a"
        raw = self.df[col].loc[row]
        used = analysis["pressure_series"]
        used_value = None if used is None else used.loc[row]
        raw_text = f"{self._fmt_value(raw)} mbar"
        if analysis["pressure_corrected"]:
            src = f"{self._fmt_value(used_value)} mbar" if used is not None else raw_text
            if self.df["j_sol_aircal"].fillna(0).astype(bool).loc[row]:
                src += f" ({analysis['pressure_aircal_smooth_s']} s aircal mean)"
            elif analysis["pressure_smooth_s"]:
                src += f" ({analysis['pressure_smooth_s']} s air mean)"
            return raw_text, src
        return raw_text, "off"

    def _corr_tooltip_temperature(self, gas_key, row, analysis):
        col = self._temperature_column(gas_key)
        if col is None:
            return "n/a", "n/a"
        raw_c = pd.to_numeric(self.df[col].loc[row], errors="coerce")
        raw_k = raw_c + 273.15 if not pd.isna(raw_c) else float("nan")
        raw_text = f"{self._fmt_value(raw_c)} C ({self._fmt_value(raw_k)} K)"
        if analysis["temperature_corrected"]:
            return raw_text, f"{T_GAS_REFERENCE_C:.0f} C"
        return raw_text, "off"

    @staticmethod
    def _fmt_value(value):
        if value is None or pd.isna(value):
            return "n/a"
        return f"{float(value):.6g}"

    def on_corr_flag_box(self, ax, x0, x1, y0, y1, unflag):
        """Flag (or unflag) the scatter points inside a dragged box.

        The box is matched against the values actually plotted -- calibrated
        for a cal-bottle gas, floor-filtered for Ozone -- because on this
        figure those *are* the axes. That is not the rule the timeseries uses
        (raw there, where two traces overlap and one has to be chosen), and it
        costs nothing: a box is resolved to row numbers once, at the moment of
        the drag, and it is the rows that get stored. Nothing re-resolves
        later, so no flag can drift when the calibration changes.

        Unflagging matches the same box rather than ignoring its height, the
        opposite of the timeseries rule -- and for the same underlying reason.
        There, a flagged point is off-screen because the y-range is framed on
        the filtered data; here it is drawn in place as a struck-out marker,
        so there is always a box the user can draw around it.
        """
        plotted = self._corr_plotted
        if self.df is None or plotted is None or ax is not self._corr_ax:
            return
        fx, fy = plotted["x"], plotted["y"]
        inside = (fx.between(min(x0, x1), max(x0, x1))
                  & fy.between(min(y0, y1), max(y0, y1)))
        # Flagged points are drawn outside `keep` (that is the whole point of
        # keeping them visible), so selection may not be restricted to it --
        # but for FLAGGING, only points that are actually part of the plotted
        # record may be taken.
        targets = [g for g in self._corr_flag_gases() if g]
        if unflag:
            selectable = inside & self._any_flag_mask(targets)
        else:
            selectable = inside & plotted["keep"]
        n = int(selectable.sum())
        if not n:
            self.corr_pane.set_stats_text(
                "Nothing flagged in that box." if unflag else
                "No plotted points in that box.")
            return

        rows = [i + self.presync_dropped
                for i in range(len(selectable)) if selectable.iat[i]]
        spans, start, prev = [], rows[0], rows[0]
        for row in rows[1:]:
            if row != prev + 1:
                spans.append((start, prev))
                start = row
            prev = row
        spans.append((start, prev))

        self._flag_undo.append(copy.deepcopy(self.flagged))
        for gas in targets:
            ranges = self.flagged.get(gas, [])
            for lo, hi in spans:
                ranges = (subtract_ranges(ranges, lo, hi) if unflag
                          else add_ranges(ranges, lo, hi))
            if ranges:
                self.flagged[gas] = ranges
            else:
                self.flagged.pop(gas, None)
        self.corr_pane.set_stats_text(
            f"{'Unflagged' if unflag else 'Flagged'} {n} point"
            f"{'' if n == 1 else 's'} on {', '.join(targets)}")
        self._after_flag_change()

    def _any_flag_mask(self, gases):
        """Rows flagged on any of `gases`."""
        mask = pd.Series(False, index=self.df.index)
        for gas in gases:
            mask |= self._flag_mask(gas)
        return mask

    def on_corr_swap_axes(self):
        if self.corr_x_gas is None or self.corr_x_gas == self.corr_y_gas:
            return
        self.corr_x_gas, self.corr_y_gas = self.corr_y_gas, self.corr_x_gas
        loading = self._loading
        self._loading = True
        self.corr_x_combo.setCurrentText(self.corr_x_gas)
        self.corr_y_combo.setCurrentText(self.corr_y_gas)
        self._loading = loading
        # Rebuilt after the swap, and it keeps the *role* rather than the gas
        # -- swapping the axes should not silently retarget the flag tool from
        # the tracer you were working on to the other one.
        self._populate_corr_flag_target()
        self._populate_corr_cal_target()
        self._refresh_corr(preserve_view=False)

    def on_corr_color_changed(self):
        """z-axis coloring on/off, or a different variable to color by."""
        if self._loading or self._initializing:
            return
        self.corr_color_combo.setEnabled(self.corr_color_check.isChecked())
        self.corr_color_by = (self.corr_color_combo.currentData()
                              if self.corr_color_check.isChecked() else None)
        self._refresh_corr(preserve_view=True)

    def on_corr_style_changed(self):
        if self._loading or self._initializing:
            return
        self.corr_marker_size = self.corr_size_spin.value()
        self.corr_error_bars = self.corr_error_check.isChecked()
        self.corr_show_cals = self.corr_cals_check.isChecked()
        self.corr_fit = self.corr_fit_check.isChecked()
        self._refresh_corr(preserve_view=True)

    def on_corr_cals_changed(self):
        if self._loading or self._initializing:
            return
        self.corr_show_cals = self.corr_cals_check.isChecked()
        self._refresh_corr(preserve_view=False)

    def _refresh_corr(self, preserve_view=True):
        """Redraw the correlation pane only.

        Not refresh(): none of these controls changes any mask, calibration
        or setting, so invalidating the shared caches and dirtying the other
        panes would make a marker-size tweak recompute two gases' worth of
        analysis and redraw the timeseries.
        """
        self._dirty["corr"] = True
        if not preserve_view:
            self._preserve["corr"] = False
        self._draw_current_tab()

    def _flash_button(self, button, message, restore, msec=1600):
        """Confirm an action in the button itself. A modal dialog for a
        one-click settings copy would cost more attention than the action is
        worth, and a status bar would be invisible next to the button.

        `restore` may be a callable, for a button whose resting label is
        computed (the Save buttons carry a dirty marker) -- restoring a
        literal string would drop a change made during the flash.
        """
        button.setText(message)
        QTimer.singleShot(msec, restore if callable(restore)
                          else (lambda: button.setText(restore)))

    def on_cal_tank_changed(self, key, _index):
        """A different tank for cal0/cal1 -- rewrites which bottles matching
        may consider, so every cal point, label and calibrated number changes.
        """
        if self._loading or self._initializing:
            return
        serial = self.tank_combos[key].currentData()
        if serial is None or serial == self.cal_selection.get(key):
            return
        self.cal_selection[key] = serial
        self._rebuild_cal_bottles()
        self._update_tank_readout()
        self._mark_dirty()
        # Full rescale, like a gas change and unlike every other control: the
        # Calibration tab's top panel plots measured MINUS ASSIGNED, so a new
        # tank shifts it by the difference in assigned values (11 ppm between
        # CC470901 and CC302489), and the calibrated overlay moves with the
        # intercept. Preserving the view there leaves the new points off-scale.
        self.refresh(preserve_view=False)

    def on_reset_cal_tanks(self):
        if self.cal_selection == self.default_cal_selection:
            return
        self.cal_selection = dict(self.default_cal_selection)
        self._rebuild_cal_bottles()
        self._apply_cal_selection_to_controls()
        self._mark_dirty()
        self.refresh(preserve_view=False)

    def on_flag_box(self, ax, x0, x1, y0, y1, unflag):
        """Flag (or unflag) the rows inside a dragged box.

        **The box is matched against the RAW value column**, never against
        whichever trace happens to be drawn. A flag has to name the same rows
        for the life of the flight: calibrated values move when the drift
        model or the cal tanks change, so a flag resolved against the red
        overlay would quietly come to mean a different set of points. The cost
        is that a box drawn tightly around the calibrated trace can catch
        nothing -- the overlay sits an intercept away from raw, ~10 ppm on
        CO2 -- which is why an empty box says so rather than doing nothing
        visible.

        Only the main axes flags: the aux panel plots a different quantity
        (detector pressure, T_gas) whose rows are not this gas's to remove.
        """
        if self.df is None or ax is not self.ax:
            return
        t0, t1 = (pd.Timestamp(mdates.num2date(x)).tz_localize(None)
                  for x in (x0, x1))
        in_time = self.df["datetime"].between(t0, t1)

        # **Unflagging ignores the y-bounds and clears the whole time span.**
        # Not an oversight and not mere leniency -- with y-bounds it would be
        # unusable in the case this feature exists for. The default y-range is
        # framed on the *filtered* series precisely so one 3500 ppb ozone
        # spike does not squash the real record, which puts the flagged point
        # off-screen: there would be no box the user could draw around a value
        # they cannot see. Time-only also makes unflagging total -- anything
        # flagged can always be removed -- and sidesteps an asymmetry, since
        # `between` is False for NaN and a flag can otherwise cover rows the
        # matching pass would refuse to select back.
        if unflag:
            selected = in_time
            n = int((selected & self._flag_mask(self.current_gas)).sum())
            if not n:
                self.main_pane.set_stats_text(
                    "Nothing flagged in that time span.")
                return
        else:
            values = self.df[GASES[self.current_gas]["value_col"]]
            selected = in_time & values.between(y0, y1)   # NaN -> False, never flagged
            n = int(selected.sum())
            if not n:
                self.main_pane.set_stats_text(
                    "Nothing in that box — the flag tool matches the raw (blue) "
                    "trace, which the calibrated overlay sits an intercept away "
                    "from.")
                return
        inside = selected

        # Positions in the analysis frame, back to RAW file rows: the whole
        # storage format is in raw numbering so it survives a change in how
        # many pre-sync rows get trimmed. Contiguous by construction here (a
        # time span), but merge_ranges makes the stored form canonical anyway.
        pos = [i + self.presync_dropped for i in range(len(inside)) if inside.iat[i]]
        spans = []
        start = prev = pos[0]
        for row in pos[1:]:
            if row != prev + 1:
                spans.append((start, prev))
                start = row
            prev = row
        spans.append((start, prev))

        targets = (list(self.config) + [g for g in self.available_gases
                                        if not GASES[g].get("has_masking", True)]
                   if self.flag_all_check.isChecked() else [self.current_gas])
        self._flag_undo.append(copy.deepcopy(self.flagged))
        for gas in dict.fromkeys(targets):
            ranges = self.flagged.get(gas, [])
            for lo, hi in spans:
                ranges = (subtract_ranges(ranges, lo, hi) if unflag
                          else add_ranges(ranges, lo, hi))
            if ranges:
                self.flagged[gas] = ranges
            else:
                self.flagged.pop(gas, None)
        self.main_pane.set_stats_text(
            f"{'Unflagged' if unflag else 'Flagged'} {n} point"
            f"{'' if n == 1 else 's'}"
            + ("" if len(targets) == 1 else f" on {len(targets)} gases"))
        self._after_flag_change()

    def _after_flag_change(self):
        """Every path that edits self.flagged ends here.

        refresh() and nothing lighter: the flags are inside `exclude_mask`, so
        they change cal means, the calibration, the uncertainty and both
        exports -- and _analysis_for is cached per gas until refresh() clears
        it. preserve_view because flagging is a fine-grained edit made while
        zoomed in on the very points being removed.
        """
        self._update_flag_readout()
        self._update_corr_flag_readout()
        self._mark_dirty()
        self.refresh(preserve_view=True)

    def _update_flag_readout(self):
        """Label + button states for the Flagging box."""
        ranges = self.flagged.get(self.current_gas, [])
        points = ranges_row_count(ranges)
        self.flag_label.setText(
            "No points flagged" if not points else
            f"{points} point{'' if points == 1 else 's'} flagged in "
            f"{len(ranges)} region{'' if len(ranges) == 1 else 's'}")
        self.flag_clear_button.setEnabled(bool(points))
        self.flag_undo_button.setEnabled(bool(self._flag_undo))

    def on_flag_undo(self):
        if not self._flag_undo:
            return
        self.flagged = self._flag_undo.pop()
        self._after_flag_change()

    def on_flag_clear(self):
        """Clear this gas's flags. Scoped to the current gas even when
        "apply to all gases" is ticked: that box describes how a new flag is
        spread, and reading it here would turn one click into a five-gas
        deletion the button never advertised."""
        if not self.flagged.get(self.current_gas):
            return
        self._flag_undo.append(copy.deepcopy(self.flagged))
        self.flagged.pop(self.current_gas, None)
        self._after_flag_change()

    def on_stats_box(self, ax, x0, x1, y0, y1):
        """Report n/mean/std for the points inside a dragged box, for whichever
        plotted trace the readout's combo box names.

        A plain selection tool: no masking is applied. The box IS the user's
        statement of which data they mean, and its vertical bounds already
        exclude the cal dives when it is drawn around the ambient band.

        The y-bounds only apply when the chosen trace lives on the Axes the box
        was drawn in. Selecting a trace from the other panel (or the calibrated
        overlay, which sits an intercept away from the raw trace) makes the
        box's y-range meaningless for it, so the selection falls back to the
        time span alone and the readout says so.
        """
        trace = self._stats_traces.get(self.main_pane.current_trace_key())
        if trace is None:
            return
        t0, t1 = (pd.Timestamp(mdates.num2date(x)).tz_localize(None)
                  for x in sorted((x0, x1)))

        same_axes = trace["axes"] is ax
        stats = box_stats(trace["x"], trace["y"], t0, t1,
                          *( (y0, y1) if same_axes else (None, None) ))

        span = f"{t0:%H:%M:%S}–{t1:%H:%M:%S}"
        if not stats["n"]:
            self.main_pane.set_stats_text(
                f"{span}   {trace['label']}: no points in box")
            return

        unit = f" {trace['unit']}" if trace["unit"] else ""
        if stats["std"] is None:
            body = f"mean {stats['mean']:.4g}{unit}  (std n/a, n=1)"
        else:
            body = f"mean {stats['mean']:.4g} ± {stats['std']:.4g}{unit}"
        parts = [span, trace["label"], f"n={stats['n']}", body,
                 f"min {stats['vmin']:.4g}  max {stats['vmax']:.4g}"]

        if not same_axes:
            parts.append("(time span only — trace is on another axis)")
        elif stats["n_clipped"]:
            parts.append(f"({stats['n_clipped']} in span, outside box vertically)")
        self.main_pane.set_stats_text("   ".join(parts))

    def on_calibrated_toggled(self, checked):
        self.show_calibrated = checked
        if self._loading or self._initializing:
            return
        self.refresh(preserve_view=True)

    # ---------------------------------------------------------------- export

    def _to_raw_rows(self, series):
        """Put an analysis Series back on the raw CSV row numbers.

        The plots use `self.df`, which starts after any pre-sync rows. The CSV
        export uses `self.raw_df`, which still has those rows. This shift keeps
        row N in the export matched to row N in the raw file.
        """
        shifted = series.reset_index(drop=True)
        offset = self.presync_dropped
        if not offset:
            return shifted
        shifted.index = shifted.index + offset
        return shifted.reindex(range(offset + len(shifted)))

    def _export_cal_diagnostics(self, gas_key, analysis):
        """Build the audit columns for one cal-bottle gas.

        These columns answer practical questions about the calibration: which
        rows made each cal mean, what that mean was, and what pressure and
        temperature values were used.
        """
        df = self.df
        settings = self._settings_for(gas_key)
        index = df.index
        used = pd.Series(False, index=index)
        mean_id = pd.Series(pd.NA, index=index, dtype="Int64")
        mean_value = pd.Series(float("nan"), index=index)
        mean_state = pd.Series(pd.NA, index=index, dtype="Int64")
        mean_serial = pd.Series(pd.NA, index=index, dtype="object")

        event_id = 0
        measured = analysis["corrected"]
        for start, end in analysis.get("cal_intervals", []):
            digital_state = bottle_for_interval(df, start, end)
            offsets = (tuple(settings["cal1_window_s"]) if digital_state == 0
                       else tuple(settings["cal2_window_s"]))
            window_start = end + pd.Timedelta(seconds=offsets[0])
            window_end = end + pd.Timedelta(seconds=offsets[1])
            window = df[(df["datetime"] >= window_start)
                        & (df["datetime"] <= window_end)]
            if window.empty:
                continue
            window = window[~analysis["exclude_mask"].loc[window.index]]
            if window.empty:
                continue
            value = measured.loc[window.index].mean()
            if pd.isna(value):
                continue
            event_id += 1
            serial = match_cal_serial(value, gas_key, self.cal_bottles)
            rows = window.index
            used.loc[rows] = True
            mean_id.loc[rows] = event_id
            mean_value.loc[rows] = value
            mean_state.loc[rows] = digital_state
            mean_serial.loc[rows] = serial

        diagnostics = {
            "cal_mean_id": self._to_raw_rows(mean_id),
            "cal_mean": self._to_raw_rows(mean_value),
            "cal_mean_state": self._to_raw_rows(mean_state),
            "cal_mean_serial": self._to_raw_rows(mean_serial),
            "corrected_input": self._to_raw_rows(analysis["corrected"]),
            "pressure_filter_mbar": self._to_raw_rows(
                analysis["pressure_filter_series"]),
        }

        detector_col = analysis.get("detector_col")
        if detector_col is not None and detector_col in df.columns:
            diagnostics["pressure_raw_mbar"] = self._to_raw_rows(
                pd.to_numeric(df[detector_col], errors="coerce"))
        if analysis.get("pressure_series") is not None:
            diagnostics["pressure_for_correction_mbar"] = self._to_raw_rows(
                analysis["pressure_series"])

        temperature_col = self._temperature_column(gas_key)
        if temperature_col is not None and temperature_col in df.columns:
            temperature_c = pd.to_numeric(df[temperature_col], errors="coerce")
            diagnostics["temperature_raw_C"] = self._to_raw_rows(temperature_c)
            if analysis.get("temperature_corrected"):
                diagnostics["temperature_delta_from_40_C"] = self._to_raw_rows(
                    temperature_c - T_GAS_REFERENCE_C)

        masks = {
            "is_cal_mean_window": self._to_raw_rows(used),
            "is_pressure_filtered": self._to_raw_rows(analysis["bad_pressure"]),
        }
        return masks, diagnostics

    def _export_gas_blocks(self):
        """Collect the prepared export columns for every available gas.

        Export is always whole-flight and all-gas. Each gas is analysed with
        its own saved settings, whether or not that gas is currently selected
        in the GUI.
        """
        blocks = []
        for gas_key in self.available_gases:
            info = GASES[gas_key]
            block = {
                "gas": gas_key,
                "short": info.get("short", gas_key),
                "long_name": info.get("long_name", gas_key),
                "standard_name": info.get("standard_name"),
                "value_col": info["value_col"],
                "unit": gas_unit(gas_key),
                "raw": self.raw_df[info["value_col"]],
                "masks": {},
            }
            if not info.get("has_masking", True):
                # Ozone and water do not use cal bottles here. Keep their
                # filtered values, but do not call them calibrated.
                removed = self._removed_mask(gas_key)
                block["final"] = self._to_raw_rows(
                    self.df[info["value_col"]].mask(removed))
                block["final_kind"] = "filtered"
                block["masks"]["below_floor"] = self._to_raw_rows(
                    self._rejected_mask(gas_key))
                block["masks"]["is_flagged"] = self._to_raw_rows(
                    self._flag_mask(gas_key))
                blocks.append(block)
                continue

            analysis = self._analysis_for(gas_key)
            result = self._calibration_for(gas_key)
            if result and result.get("ok"):
                sigma, _ = self._uncertainty_for(gas_key)
                audit_masks, diagnostics = self._export_cal_diagnostics(
                    gas_key, analysis)
                block["final"] = self._to_raw_rows(result["calibrated"])
                block["final_kind"] = "calibrated"
                # These fields describe the delivered values, so both CSV and
                # ICARTT can report the same correction settings.
                block["pressure_corrected"] = analysis["pressure_corrected"]
                block["pressure_col"] = analysis["pressure_col"]
                # Only report pressure smoothing when the pressure correction
                # was actually applied.
                block["pressure_smooth_s"] = (
                    analysis["pressure_smooth_s"]
                    if block["pressure_corrected"] else 0)
                block["temperature_corrected"] = analysis["temperature_corrected"]
                block["temperature_col"] = analysis["temperature_col"]
                # Keep the P/T factor with the slope and intercept so a user
                # can reproduce the calibrated value from the raw value.
                block["pt_corrected"] = result.get("correction_factor") is not None
                block["pt_factor"] = (
                    None if result.get("correction_factor") is None
                    else self._to_raw_rows(result["correction_factor"]))
                block["sigma"] = self._to_raw_rows(sigma)
                block["slope"] = self._to_raw_rows(result["slope"])
                block["intercept"] = self._to_raw_rows(result["intercept"])
                block["diagnostics"] = diagnostics
                block["masks"] = {
                    "is_cal_period": self._to_raw_rows(result["in_cal"]),
                    "is_post_cal_flush": self._to_raw_rows(result["flushed"]),
                    "is_masked": self._to_raw_rows(result["excluded"]),
                    "is_extrapolated": self._to_raw_rows(result["extrapolated"]),
                    # Manual flags are also part of is_masked, but they are
                    # useful to see on their own.
                    "is_flagged": self._to_raw_rows(analysis["flagged"]),
                }
                block["masks"].update(audit_masks)
            else:
                # Even without a usable calibration, export the masks so the
                # raw column still has the air/non-air decisions beside it.
                block["reason"] = (result or {}).get("reason", "no calibration")
                block["masks"] = {
                    "is_cal_period": self._to_raw_rows(analysis["not_air"]),
                    "is_post_cal_flush": self._to_raw_rows(analysis["post_cal_flush"]),
                    "is_masked": self._to_raw_rows(analysis["exclude_mask"]),
                    "is_flagged": self._to_raw_rows(analysis["flagged"]),
                }
            blocks.append(block)
        return blocks

    def _update_export_summary(self):
        """What the exports would contain right now, gas by gas.

        Recomputed through the same _analysis_for/_calibration_for caches the
        plots use, so the tab cannot claim a gas will export calibrated while
        the Calibration tab is showing that it has no usable bottles.
        """
        if not hasattr(self, "export_summary_label"):
            return
        for widget in (self.export_csv_button, self.export_icartt_button):
            widget.setEnabled(self.df is not None)
        if self.df is None:
            self.export_summary_label.setText("No file loaded.")
            return

        lines = []
        for block in self._export_gas_blocks():
            gas = block["gas"]
            if block.get("final_kind") == "calibrated":
                n = int(block["final"].notna().sum())
                smooth_s = block.get("pressure_smooth_s") or 0
                corrections = []
                if block.get("pressure_corrected"):
                    corrections.append(
                        f"P to {D1_P_TARGET_MBARS:.0f} mbar"
                        + (f" [{smooth_s} s mean]" if smooth_s else ""))
                if block.get("temperature_corrected"):
                    corrections.append(f"T from {T_GAS_REFERENCE_C:.0f} C")
                lines.append(f"{gas:<6} calibrated      {n:>7,} good ambient rows"
                             + (f"   (normalised before calibrating: "
                                f"{', '.join(corrections)})"
                                if corrections else ""))
            elif block.get("final_kind") == "filtered":
                n = int(block["final"].notna().sum())
                removed = int(block["masks"]["below_floor"].sum())
                lines.append(f"{gas:<6} filtered only   {n:>7,} rows"
                             + (f"   ({removed} below floor removed)" if removed else ""))
            else:
                lines.append(f"{gas:<6} NOT calibrated  {block.get('reason', '')}")
        rows = len(self.raw_df)
        lines.append("")
        lines.append(f"{rows:,} rows in {self.csv_path.name}"
                     + (f"  ({self.presync_dropped} pre-sync, blank in the CSV export)"
                        if self.presync_dropped else ""))
        _, _, usable = self._export_time_base()
        unusable = int((~usable).sum())
        if unusable:
            lines.append(f"{unusable:,} row(s) have repeated or backward "
                         f"timestamps and cannot go in the ICARTT file.")
        self.export_summary_label.setText("\n".join(lines))

    def _export_time_base(self):
        """(start date, seconds-from-midnight, usable-row mask) for the raw
        row set. Wrapped so the pre-sync count is passed at every call site --
        forgetting it does not raise, it silently rejects most of the flight
        (see icartt_time_base)."""
        return icartt_time_base(self.raw_df["datetime"], self.presync_dropped)

    def _export_default_path(self, suffix):
        if not self.csv_path:
            return ""
        return str(self.csv_path.with_name(f"{self.csv_path.stem}{suffix}"))

    def on_export_csv_clicked(self):
        """Write the companion CSV: every gas, every row of the raw file."""
        if self.df is None:
            return
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Export derived CSV", self._export_default_path("_derived.csv"),
            "CSV Files (*.csv);;All Files (*)")
        if not path_str:
            return
        _, seconds, _ = self._export_time_base()
        try:
            summary = export_companion_csv(
                Path(path_str), self.raw_df["datetime"], self._export_gas_blocks(),
                source_path=self.csv_path,
                include_raw=self.csv_raw_check.isChecked(),
                include_masks=self.csv_masks_check.isChecked(),
                include_coefficients=self.csv_coeff_check.isChecked(),
                include_uncertainty=self.csv_unc_check.isChecked(),
                presync_rows=self.presync_dropped,
                comment_header=self.csv_comment_check.isChecked(),
                time_seconds=seconds,
            )
        except OSError as e:
            QMessageBox.warning(self, "Export CSV", f"Could not write {path_str}:\n{e}")
            return
        message = (f"Wrote {summary['path'].name}\n\n"
                   f"{summary['rows']:,} rows × {len(summary['columns'])} columns.")
        QMessageBox.information(self, "Export CSV", message)

    def on_export_icartt_clicked(self):
        """Write the ICARTT file, after warning about anything the format
        will silently drop."""
        if self.df is None:
            return
        meta = self._icartt_meta_from_controls()
        missing = [label for label, key in
                   (("PI name", "pi_name"), ("PI affiliation", "pi_affiliation"),
                    ("Mission", "mission"), ("Platform", "platform"),
                    ("Location ID", "location_id"))
                   if not str(meta.get(key, "")).strip()]
        if missing:
            answer = QMessageBox.question(
                self, "Incomplete metadata",
                "These header fields are empty and will be written as N/A "
                "(or defaulted in the file name):\n\n  " + "\n  ".join(missing)
                + "\n\nAn archive will usually reject a file like this. "
                  "Export anyway?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if answer != QMessageBox.Yes:
                return

        start_date, _, _ = self._export_time_base()
        default = ""
        if self.csv_path and start_date is not None:
            default = str(self.csv_path.with_name(icartt_filename(meta, start_date)))
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Export ICARTT", default, "ICARTT Files (*.ict);;All Files (*)")
        if not path_str:
            return
        try:
            summary = export_icartt(
                Path(path_str), self.raw_df["datetime"], self._export_gas_blocks(),
                meta=meta,
                include_sigma=self.icartt_sigma_check.isChecked(),
                drop_empty_rows=self.icartt_drop_check.isChecked(),
                skip_leading=self.presync_dropped,
            )
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Export ICARTT", f"Could not write {path_str}:\n{e}")
            return
        message = (f"Wrote {summary['path'].name}\n\n"
                   f"{summary['rows']:,} data rows, {summary['header_lines']} "
                   f"header lines.\nVariables: "
                   + ", ".join(summary["variables"]))
        # Said out loud rather than left in the header: a file that is shorter
        # than the flight is the one thing about this format most likely to be
        # noticed later and mistaken for lost data.
        if summary["unusable_times"]:
            message += (f"\n\n{summary['unusable_times']:,} row(s) omitted — "
                        f"repeated or backward timestamps, which ICARTT's "
                        f"independent variable cannot represent.")
        if summary["empty_rows"]:
            message += (f"\n\n{summary['empty_rows']:,} row(s) omitted as "
                        f"having no value for any variable.")
        QMessageBox.information(self, "Export ICARTT", message)

    def on_aux_changed(self, selection: str):
        # Fires with "" while the combo is being repopulated on a file load;
        # the caller re-selects "No Figure" itself, so there is nothing to do.
        if not selection:
            return
        self.aux_selection = selection
        self.other_combo.setEnabled(self.aux_selection == "Other")
        # The right-axis trace twins whatever's on the aux panel, so it only
        # makes sense once that panel exists at all.
        self.right_axis_combo.setEnabled(self.aux_selection != "No Figure")
        if self._initializing:
            return
        self.refresh(preserve_view=True)

    def on_other_changed(self, text: str):
        self.other_column = text or None
        if self._initializing or self.aux_selection != "Other":
            return
        self.refresh(preserve_view=True)

    def on_right_axis_changed(self, text: str):
        self.right_axis_column = None if text in ("", "(none)") else text
        if self._initializing:
            return
        self.refresh(preserve_view=True)

    def refresh(self, preserve_view=False):
        """Invalidate the cached analysis, redraw the visible pane, and mark
        the other one dirty so it redraws when it is next shown.

        The single entry point for every state change, so the caches have
        exactly one invalidation site. Deliberately invalidates
        unconditionally rather than comparing a composite key of
        (file, gas, warm-up, tolerance, cal windows, ...): such a key is easy
        to get subtly wrong and then serves a stale plot, whereas
        invalidate-everything cannot. The recompute is milliseconds against
        ~100 ms of rendering, and redrawing only the visible pane avoids
        paying even that twice.
        """
        self._analysis.clear()
        self._calibration.clear()
        self._uncertainty.clear()
        for pane in self._dirty:
            self._dirty[pane] = True
            # A requested full rescale has to survive until it is actually
            # honoured: changing gas (rescale) and then nudging a spinbox
            # (preserve) must still rescale the pane that never redrew.
            if not preserve_view:
                self._preserve[pane] = False
        self._draw_current_tab()

    def _draw_current_tab(self):
        """Draw whichever pane is showing, if it is dirty.

        Keyed on the current *widget*, not the tab index: Cal Tanks draws
        nothing, and an index test would have made it redraw the timeseries
        (or, once a fourth tab appears, whatever else fell through the else).
        """
        if self.df is None:
            return
        widget = self.tabs.currentWidget()
        if widget is self.cal_pane:
            name = "cal"
        elif widget is self.corr_pane:
            name = "corr"
        elif widget is self.main_pane:
            name = "main"
        elif widget is self.export_pane:
            name = "export"
        else:
            return
        if not self._dirty.get(name):
            return
        preserve = self._preserve.get(name, False)
        if name == "cal":
            self.redraw_cal(preserve_view=preserve)
        elif name == "corr":
            self.redraw_corr(preserve_view=preserve)
        elif name == "export":
            # Not a plot, but it does read every gas's analysis, so it goes
            # through the same dirty dispatch rather than recomputing five
            # gases on every spin-box nudge while the tab happens to be open.
            self._update_export_summary()
        else:
            self.redraw(preserve_view=preserve)
        self._dirty[name] = False
        self._preserve[name] = True

    def on_tab_changed(self, index):
        # The control panel follows the tab: correlation controls for the
        # Correlations tab, the per-gas panel for everything else.
        self.controls_stack.setCurrentIndex(
            1 if self.tabs.currentWidget() is self.corr_pane else 0)
        if self.tabs.currentWidget() is self.corr_pane:
            self._sync_corr_cal_controls()
        if self._initializing:
            return
        self._draw_current_tab()

    def _settings_for(self, gas_key):
        """The settings that apply to `gas_key`.

        The live controls for the gas on display, `self.config` for any other
        -- the Correlations tab analyses two gases at once, only one of which
        the control panel is showing. The two agree in practice (every edit
        writes `config[current_gas]`), but reading the widgets for the current
        gas keeps the displayed controls authoritative for what is drawn.
        """
        if gas_key == self.current_gas:
            return self._controls_to_settings()
        return self.config.get(gas_key, DEFAULT_GAS_SETTINGS)

    def _get_analysis(self):
        """Masks, cal intervals and per-injection cal means for the gas on
        display. See _analysis_for."""
        return self._analysis_for(self.current_gas)

    def _rejected_mask(self, gas_key):
        """Rows where this gas's sensor read below its physical floor.

        Separate from `_analysis_for`'s masks and not cached with them: it
        depends on nothing the user can change, only on the gas's declared
        `valid_min`, and it applies to gases (Ozone) that have no analysis
        settings at all. Kept floor-only on purpose -- the note text
        distinguishes "below floor" from "flagged by hand", and conflating
        them here would make that impossible.
        """
        floor = GASES[gas_key].get("valid_min")
        return below_floor_mask(self.df[GASES[gas_key]["value_col"]], floor)

    def _flag_mask(self, gas_key):
        """Rows this gas has been manually flagged on, as a boolean Series.

        The inverse of _to_raw_rows: the ranges are stored in the RAW file's
        row numbering and the analysis frame starts at raw row
        `presync_dropped`, so the offset is subtracted here and nowhere else.
        """
        if self.df is None or gas_key is None:
            return pd.Series(False, index=[] if self.df is None else self.df.index)
        return ranges_to_mask(self.flagged.get(gas_key, []), self.df.index,
                              offset=self.presync_dropped)

    def _removed_mask(self, gas_key):
        """Everything removed from a FLOOR gas's record: below-floor faults
        plus manual flags.

        Only for gases with has_masking=False (Ozone, H2O). A cal-bottle gas
        gets its flags through `exclude_mask` in _analysis_for instead, which
        is what makes them drop cal points as well as blanking output; those
        gases never reach here because they have no `valid_min` and no
        calibrate_series-free display path.
        """
        return self._rejected_mask(gas_key) | self._flag_mask(gas_key)

    def _analysis_for(self, gas_key):
        """Masks, cal intervals and per-injection cal means for one gas, from
        the current file and that gas's settings.

        Cached per gas because more than one view reads the same numbers and
        they must not disagree; see refresh() for how the cache is
        invalidated. Returns None when no file is loaded.
        """
        if self.df is None or gas_key is None:
            return None
        if gas_key in self._analysis:
            return self._analysis[gas_key]

        gas = GASES[gas_key]
        settings = self._settings_for(gas_key)
        df = self.df
        has_masking = gas.get("has_masking", True)
        warmup_minutes = settings["warmup_min"]
        pressure_tol = settings["pressure_tol_mbar"]
        flag_air_s = settings["flag_air_s"]

        # bad_pressure/warmup are computed unconditionally -- even for a gas
        # with no masking of its own (Ozone), the aux panel can still show
        # an Aeris trace (e.g. Detector Pressure) that these masks apply to.
        # Only their use on the *main* plot (shading/notes/cal exclusion) is
        # gated on has_masking by the caller.
        cal = (df["j_sol_cals"].fillna(0).astype(bool)
               | df["j_sol_aircal"].fillna(0).astype(bool))
        # The valve flag leads the measurement by about a sample, so the first
        # "ambient" row after an injection is still tank gas. Kept as its own
        # mask and OR-ed into `not_air` below rather than into `cal` itself:
        # `cal` defines the cal INTERVALS, whose ends are the Cal_p every cal
        # mean window is measured from, and moving those would silently shift
        # every cal mean.
        cal_switch = cal_switch_mask(df["datetime"], cal)
        not_air = cal | cal_switch

        # The column the pressure correction will divide by, or None when it
        # is switched off (or the gas has no detector of its own). Resolved
        # here so every view can name it without re-deriving the rule.
        detector_col = self._pressure_column(gas_key)
        pressure_col = (detector_col
                        if settings.get("pressure_correct", False) else None)

        # The pressure used by the correction and, when pressure correction is
        # on, by the pressure filter. The selectable smoothing window is for air
        # rows. During j_sol_aircal periods the transition is sharper and the
        # correction uses a fixed 5 s smooth instead, so a long air window cannot
        # pull cal-pressure transitions into the cal means.
        smooth_s = settings.get("pressure_smooth_s", 0)
        pressure_smoothed = None
        pressure_air_smoothed = None
        pressure_aircal_smoothed = None
        if detector_col is not None:
            raw_pressure = df[detector_col]
            pressure_air_smoothed = smooth_pressure(
                df["datetime"], raw_pressure, smooth_s)
            pressure_aircal_smoothed = smooth_pressure(
                df["datetime"], raw_pressure, 5)
            aircal = df["j_sol_aircal"].fillna(0).astype(bool)
            pressure_smoothed = pressure_air_smoothed.mask(
                aircal, pressure_aircal_smoothed)
        # What the correction actually divides by: smoothed when there is a
        # window, the raw column when there is not, None when the correction
        # is off. One key, so no call site has to re-derive that choice.
        pressure_series = None
        if pressure_col is not None:
            pressure_series = pressure_smoothed

        # The pressure filter follows the same pressure basis that the product is
        # corrected with. With pressure correction off it remains the raw d1
        # detector pressure historical mask; with correction on it filters the
        # corrected pressure series, including the fixed 5 s j_sol_aircal smooth.
        pressure_filter_series = (pressure_series if pressure_series is not None
                                  else pd.to_numeric(df["d1_P_mbars"], errors="coerce"))
        bad_pressure = ((pressure_filter_series - D1_P_TARGET_MBARS).abs()
                        > pressure_tol)
        bad_pressure = bad_pressure.fillna(False)

        # The temperature correction's column, or None when it is off. Not
        # smoothed: the cell temperature is already smooth (see the checkbox's
        # tooltip), so the raw column is what the T-40 correction is taken from.
        temperature_col = (self._temperature_column(gas_key)
                           if settings.get("temperature_correct", False)
                           else None)

        # THE CORRECTION IS APPLIED TO THE MEASUREMENT, before anything else
        # (2026-07-31). `corrected` is what the cal means are averaged from and
        # what gets calibrated, so the bottle responses, the drift nodes,
        # slope/intercept and span_gain are all on the corrected scale -- the
        # correction now moves the calibration, deliberately. Built once, here,
        # and handed to both cal_mean_points and calibrate_series so the two
        # cannot end up on different scales.
        correction = pt_correction_factor(
            pressure=pressure_series,
            temperature_c=(None if temperature_col is None else df[temperature_col]),
        )
        value_col = gas["value_col"]
        corrected = (df[value_col] if correction is None
                     else df[value_col] * correction)

        warmup_end = df["datetime"].iloc[0] + pd.Timedelta(minutes=warmup_minutes)
        warmup = df["datetime"] < warmup_end

        # Measured back from the last timestamp in the record, not from a
        # clock time: "the last 10 minutes" means the last 10 minutes of data
        # there are, whenever the file happens to stop.
        end_flight_minutes = settings.get("end_flight_min", 0)
        end_flight_start = df["datetime"].iloc[-1] - pd.Timedelta(minutes=end_flight_minutes)
        end_flight = ((df["datetime"] > end_flight_start) if end_flight_minutes
                      else pd.Series(False, index=df.index))
        # One band, one note line: both are "the instrument was not doing what
        # the rest of the flight was doing", and the user asked for them to
        # read as a single exclusion at the two ends.
        trimmed = warmup | end_flight

        # Pumps off: not ambient air at all, so this joins exclude_mask -- it
        # drops cal points as well as blanking the output, exactly like the
        # warm-up and pressure masks. A missing j_pumps reading counts as off:
        # an unknown pump state is not evidence the pumps were running.
        require_pumps = bool(settings.get("require_pumps", False)) and "j_pumps" in df.columns
        pumps_off = (df["j_pumps"].fillna(0) != 1 if require_pumps
                     else pd.Series(False, index=df.index))

        # Manually flagged rows join exclude_mask rather than getting a
        # channel of their own, which is the whole design in one line: that
        # mask is already handed BOTH to cal_mean_points (where it drops raw
        # rows before the cal means are estimated, so flagging a visibly bad
        # injection actually changes the calibration) and to calibrate_series
        # (where it blanks the finished output). Flags therefore behave
        # exactly like the warm-up, pressure and pumps masks, which is what
        # the user asked for and what makes them predictable. Kept separately
        # in the dict below for the markers and the note, which must be able
        # to say *why* a row went.
        flagged = self._flag_mask(gas_key)
        exclude_mask = bad_pressure | trimmed | pumps_off | flagged

        cal_intervals, cal_points, display_cal_points = [], [], []
        post_cal_flush = pd.Series(False, index=df.index)
        if has_masking:
            cal_intervals = merge_close_intervals(
                find_intervals(df["datetime"], cal), pd.Timedelta(seconds=CAL_MERGE_GAP_S)
            )
            # Deliberately not folded into exclude_mask: this flags *ambient*
            # rows after an injection, which is disjoint from the cal windows
            # the means are computed over, and folding it in would suggest it
            # can drop a cal point (it cannot).
            post_cal_flush = post_cal_flush_mask(
                df["datetime"], cal_intervals, flag_air_s, cal_mask=not_air
            )
            # Cal means are estimated with these masks applied -- a cal point
            # can be dropped entirely if its window has no valid data -- and
            # from the P/T-CORRECTED measurement, which is what makes the
            # correction reach the calibration itself.
            cal_points = cal_mean_points(
                df, cal_intervals, value_col,
                tuple(settings["cal1_window_s"]),
                tuple(settings["cal2_window_s"]),
                cal_bottles=self.cal_bottles, gas_key=gas_key,
                exclude_mask=exclude_mask, values=corrected,
            )
            # Timeseries cal markers stay in the instrument's recorded units.
            # The calibration itself uses `cal_points` above, which may be
            # P/T-corrected; these points are only for visual reference.
            display_cal_points = cal_mean_points(
                df, cal_intervals, value_col,
                tuple(settings["cal1_window_s"]),
                tuple(settings["cal2_window_s"]),
                cal_bottles=self.cal_bottles, gas_key=gas_key,
                exclude_mask=exclude_mask,
            )

        self._analysis[gas_key] = {
            "cal": cal, "cal_switch": cal_switch, "not_air": not_air,
            "warmup": warmup, "end_flight": end_flight, "trimmed": trimmed,
            "bad_pressure": bad_pressure,
            "pumps_off": pumps_off, "require_pumps": require_pumps,
            "flagged": flagged,
            "exclude_mask": exclude_mask,
            "cal_intervals": cal_intervals, "cal_points": cal_points,
            "display_cal_points": display_cal_points,
            "post_cal_flush": post_cal_flush,
            "has_masking": has_masking,
            "warmup_minutes": warmup_minutes,
            "end_flight_minutes": end_flight_minutes,
            "pressure_tol": pressure_tol,
            # Not a mask, but it belongs to the same reading and every view
            # that describes the calibrated record has to be able to say so
            # without reaching for the calibration itself (which the
            # timeseries only computes when the overlay is on).
            # True when it is actually being applied, not merely asked for:
            # the column has to exist for this gas on this flight.
            "pressure_corrected": pressure_col is not None,
            "pressure_col": pressure_col,
            "detector_col": detector_col,
            "pressure_smooth_s": smooth_s if detector_col is not None else 0,
            "pressure_aircal_smooth_s": 5 if detector_col is not None else 0,
            "pressure_filter_series": pressure_filter_series,
            "pressure_smoothed": pressure_smoothed,
            "pressure_air_smoothed": pressure_air_smoothed,
            "pressure_aircal_smoothed": pressure_aircal_smoothed,
            "pressure_series": pressure_series,
            "temperature_corrected": temperature_col is not None,
            "temperature_col": temperature_col,
            # The multiplier and the corrected measurement it produced. Both
            # travel in the analysis so the calibration, the exports and the
            # notes all read the identical numbers rather than each rebuilding
            # them from the settings.
            "correction_factor": correction,
            "corrected": corrected,
            "flag_air_s": flag_air_s,
        }
        return self._analysis[gas_key]

    def redraw(self, preserve_view=False):
        analysis = self._get_analysis()
        if analysis is None:
            return
        gas = GASES[self.current_gas]
        value_col = gas["value_col"]
        has_masking = analysis["has_masking"]
        warmup_minutes = analysis["warmup_minutes"]
        pressure_tol = analysis["pressure_tol"]

        df = self.df

        aux_info = aux_trace_info(self.aux_selection, self.current_gas, self.other_column)
        has_aux_panel = aux_info is not None
        aux_key = (self.aux_selection, self.other_column)
        # The right-axis trace twins the aux panel, so it can only appear
        # alongside one.
        has_right_axis = has_aux_panel and self.right_axis_column is not None
        right_axis_key = self.right_axis_column if has_right_axis else None

        # Capture the current view before tearing down the old Axes, so a
        # masking/averaging control change -- or switching/adding/removing
        # the upper trace -- can redraw without rescaling the main plot.
        # The aux panel's own y-range is only worth preserving if it's
        # still showing the same trace (its scale means something different
        # for a different trace, so let that one re-autoscale) -- for
        # "Other" that also means the same catch-all column, not just the
        # same radio button. The right-axis trace is tracked independently
        # since it can change without the left trace changing (or vice versa).
        old_main_view = None
        old_aux_ylim = None
        old_right_ylim = None
        if preserve_view and self.ax is not None:
            old_main_view = (self.ax.get_xlim(), self.ax.get_ylim())
            if (self.ax_aux is not None and has_aux_panel
                    and self._last_aux_key == aux_key):
                old_aux_ylim = self.ax_aux.get_ylim()
            if (self.ax_aux2 is not None and has_right_axis
                    and self._last_right_axis_key == right_axis_key):
                old_right_ylim = self.ax_aux2.get_ylim()

        self.figure.clear()
        if has_aux_panel:
            gs = self.figure.add_gridspec(2, 1, height_ratios=[1, 3])
            ax_aux = self.figure.add_subplot(gs[0])
            ax = self.figure.add_subplot(gs[1], sharex=ax_aux)
        else:
            ax_aux = None
            ax = self.figure.add_subplot(111)
        ax.set_facecolor("#fcfcfb")

        cal = analysis["cal"]
        # The band tracks what is treated as not-air, switch-over sample
        # included, so the shading and the blanked calibrated trace agree.
        not_air = analysis["not_air"]
        warmup = analysis["warmup"]
        # Warm-up and end-of-flight share the band and the note: one orange
        # exclusion at each end of the record.
        trimmed = analysis["trimmed"]
        end_flight_minutes = analysis["end_flight_minutes"]
        bad_pressure = analysis["bad_pressure"]
        cal_points = analysis["cal_points"]
        display_cal_points = analysis["display_cal_points"]
        post_cal_flush = analysis["post_cal_flush"]

        if has_masking:
            shade_intervals(ax, df["datetime"], not_air, CAL_SHADE_COLOR, alpha=0.3)
            shade_intervals(ax, df["datetime"], trimmed, WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax, df["datetime"], analysis["pumps_off"],
                            PUMPS_EXCLUDE_COLOR, alpha=0.16)
            # Shaded whether or not the calibrated trace is showing: the band
            # is how you find out these rows exist before turning it on.
            shade_intervals(ax, df["datetime"], post_cal_flush, POST_CAL_FLUSH_COLOR, alpha=0.22)

        # The calibrated trace is opt-in; the raw trace stays visible
        # underneath it (faded) so the correction being applied is always
        # legible rather than silently swapped in.
        calibration = self._get_calibration() if self.show_calibrated else None
        show_cal = bool(calibration and calibration.get("ok"))

        # Registered as they are plotted rather than scraped back off the Axes
        # afterwards: the artists carry no units and no stable identity, and a
        # trace that is conditionally drawn would be easy to miss.
        self._stats_traces = {}
        unit = gas_unit(self.current_gas)
        self._register_stats_trace(
            "main:raw", f"{self.current_gas} (raw)", ax, df["datetime"],
            df[value_col], unit)

        # A gas with a physical floor (Ozone) gets the same two-trace
        # treatment as a calibrated one -- raw blue underneath, the kept data
        # in red on top -- deliberately reusing LINE_COLOR/CALIBRATED_COLOR
        # rather than inventing a palette: on both figures red means "the
        # series you should be reading" and blue means "everything the
        # instrument recorded".
        # For a floor gas this is the only place removals happen, so manual
        # flags have to join here; a cal-bottle gas gets them through
        # exclude_mask -> calibrate_series instead (see _removed_mask).
        rejected = self._rejected_mask(self.current_gas)
        flagged = analysis["flagged"]
        removed = rejected | flagged
        show_filtered = bool(removed.any()) and not has_masking

        plot_data = df[["datetime", value_col]].dropna()
        line, = ax.plot(plot_data["datetime"], plot_data[value_col], color=LINE_COLOR,
                        linewidth=1.2, alpha=0.55 if (show_cal or show_filtered) else 1.0)

        filtered_line = None
        if show_filtered:
            # Rejected rows stay as NaN rather than being dropped, so the red
            # line breaks over them instead of drawing across the removal --
            # same reasoning as the calibrated trace below.
            filtered = df[value_col].mask(removed)
            self._register_stats_trace(
                "main:filtered", f"{self.current_gas} (filtered)", ax,
                df["datetime"], filtered, unit)
            keep = filtered.notna() | removed
            filtered_line, = ax.plot(df["datetime"][keep], filtered[keep],
                                     color=CALIBRATED_COLOR, linewidth=1.2)

        # Manually flagged points, struck out at their RAW values -- which is
        # also the basis the flag box matches against, so the marker sits
        # exactly where the user dragged. Drawn for every gas, and on top of
        # everything, because "I removed this" has to stay visible after the
        # trace it was removed from has broken over the gap.
        flag_scatter = None
        if flagged.any():
            flag_scatter = ax.scatter(
                df["datetime"][flagged], df[value_col][flagged],
                marker="x", s=28, linewidths=1.1, color=FLAGGED_COLOR, zorder=6)

        cal_line = None
        cal_mean_scatter = None
        if show_cal:
            # The calibrated trace is the *good air* record: cal periods, the
            # flush behind them, and the warm-up/bad-pressure masks are all
            # blanked by calibrate_series. Those rows are kept as NaN instead
            # of being dropped, so the line breaks over them; dropping them
            # would draw a straight segment across each gap and hide the
            # removal.
            calibrated = calibration["calibrated"]
            self._register_stats_trace(
                "main:cal", f"{self.current_gas} (calibrated)", ax,
                df["datetime"], calibrated, unit)
            keep = calibrated.notna() | calibration["blanked"]
            cal_df = pd.DataFrame({"datetime": df["datetime"][keep],
                                   "v": calibrated[keep]})
            cal_line, = ax.plot(cal_df["datetime"], cal_df["v"],
                                color=CALIBRATED_COLOR, linewidth=1.2)
            for start, end in find_intervals(df["datetime"], calibration["extrapolated"]):
                ax.axvspan(start, end, facecolor="none", edgecolor=CAL_SHADE_COLOR,
                           hatch="///", alpha=0.30, linewidth=0)
            calibrated_cal_points = self._calibrated_cal_points(calibration)
            if calibrated_cal_points:
                xs, ys, _ = zip(*calibrated_cal_points)
                cal_mean_scatter = ax.scatter(
                    xs, ys, facecolors="none", edgecolors=CALIBRATED_COLOR,
                    s=52, linewidths=1.2, zorder=5)

        # cal0_pts/cal1_pts are guaranteed empty when has_masking is False, so
        # cal0_label/cal1_label are never read below
        # without having been set here first.
        cal0_pts = [(t, v) for t, v, state, serial in display_cal_points if state == 0]
        cal1_pts = [(t, v) for t, v, state, serial in display_cal_points if state == 1]
        if has_masking:
            cal0_label = most_common_serial(cal_points, 0) or "Cal 1"
            cal1_label = most_common_serial(cal_points, 1) or "Cal 2"
            self._set_cal_row_label(self.cal1_label, cal0_label, "Cal 1")
            self._set_cal_row_label(self.cal2_label, cal1_label, "Cal 2")

        handles = [line]
        # Kept short when both traces are shown -- the y-axis already names
        # the gas, and a four-entry legend clips against the axes edge.
        labels = ["raw (ambient)" if (show_cal or show_filtered)
                  else f"{gas['ylabel']} (ambient)"]
        if filtered_line is not None:
            handles.append(filtered_line)
            labels.append(f"filtered (≥ {GASES[self.current_gas]['valid_min']:g} {unit})")
        if cal_line is not None:
            handles.append(cal_line)
            labels.append(f"calibrated ({calibration['mode']})")
        if cal_mean_scatter is not None:
            handles.append(cal_mean_scatter)
            labels.append("cal means (calibrated)")
        if cal0_pts:
            xs, ys = zip(*cal0_pts)
            handles.append(ax.scatter(xs, ys, color=CAL0_COLOR, s=40, zorder=5, edgecolors="none"))
            labels.append(f"{cal0_label}: {mean_std_label(ys)}")
        if cal1_pts:
            xs, ys = zip(*cal1_pts)
            handles.append(ax.scatter(xs, ys, color=CAL1_COLOR, s=40, zorder=5, edgecolors="none"))
            labels.append(f"{cal1_label}: {mean_std_label(ys)}")
        if flag_scatter is not None:
            handles.append(flag_scatter)
            labels.append(f"flagged ({int(flagged.sum())})")

        ax.set_ylabel(gas["ylabel"], color=TEXT_COLOR)
        # Substituted, not appended: the stock title already says
        # "(uncalibrated)", which would otherwise contradict the suffix.
        title = gas["title"].replace("(uncalibrated)", "(calibrated)") if show_cal \
            else gas["title"]
        ax.set_title(title, color=TEXT_COLOR, loc="left")

        date_str = plot_data["datetime"].iloc[0].strftime("%Y-%m-%d") if not plot_data.empty else ""
        ax.set_xlabel(f"Time (UTC-ish, {date_str})", color=MUTED_COLOR)

        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

        ax.grid(True, color=GRID_COLOR, linewidth=0.8)
        for spine in ax.spines.values():
            spine.set_color(AXIS_COLOR)
        ax.tick_params(colors=MUTED_COLOR)
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment("right")

        notes = []
        if has_masking:
            if cal.any():
                notes.append("gray = calibration/cal-air (j_sol_cals, j_sol_aircal"
                             + (", + switch-over sample)"
                                if analysis["cal_switch"].any() else ")"))
            if trimmed.any():
                ends = []
                if warmup.any():
                    ends.append(f"first {warmup_minutes} min warm-up")
                if analysis["end_flight"].any():
                    ends.append(f"last {end_flight_minutes} min of flight")
                notes.append("orange = excluded (" + ", ".join(ends) + ")")
            if bad_pressure.any():
                pressure_basis = (analysis["pressure_col"] or "d1_P_mbars")
                if analysis["pressure_corrected"]:
                    smooth_note = (
                        f"; {analysis['pressure_smooth_s']} s air mean"
                        if analysis["pressure_smooth_s"] else "; raw air")
                    smooth_note += (
                        f", {analysis['pressure_aircal_smooth_s']} s aircal mean")
                    pressure_basis += smooth_note
                notes.append(
                    f"light red = excluded ({pressure_basis} outside "
                    f"{D1_P_TARGET_MBARS:.0f}±{pressure_tol:.2f} mbar)"
                )
            if analysis["pumps_off"].any():
                notes.append(
                    f"violet = excluded ({int(analysis['pumps_off'].sum())} rows with the "
                    f"pumps off, j_pumps ≠ 1)"
                )
            if post_cal_flush.any():
                notes.append(
                    f"teal = air dropped from calibrated only "
                    f"({analysis['flag_air_s']} s detector flush after each cal)"
                )
            if show_cal:
                notes.append("red = calibrated, blue = raw; calibrated shows "
                             "good air only (cal periods, flush and masked "
                             "spans blanked)")
                if cal_mean_scatter is not None:
                    notes.append("open red circles = cal means on the calibrated scale")
            # Said whether or not the overlay is on: it changes the calibrated
            # record itself, which is what both exports ship, so the figure
            # would otherwise show a raw trace with no hint that the product
            # behind it has been rescaled.
            # One line for both corrections, since they are one operation on
            # the number as far as this figure is concerned -- and because a
            # line each would push the note block up over the data.
            if analysis["pressure_corrected"] or analysis["temperature_corrected"]:
                scaled_to, terms = [], []
                if analysis["pressure_corrected"]:
                    # The smoothing is named here rather than in a line of its
                    # own: it is a property of the P this scales by, and off
                    # the Detector Pressure panel it has no other visible
                    # effect.
                    smooth_s = analysis["pressure_smooth_s"]
                    scaled_to.append(f"{D1_P_TARGET_MBARS:.0f} mbar")
                    terms.append(
                        f"{D1_P_TARGET_MBARS:.0f}/{analysis['pressure_col']}"
                        # Just "[30 s mean]": with both corrections on, this
                        # line is already the widest thing on the figure, and
                        # the tooltip and the export notes carry the full
                        # "centred, whole flight" description.
                        + (f" [{smooth_s} s air mean, "
                           f"{analysis['pressure_aircal_smooth_s']} s aircal]"
                           if smooth_s else
                           f" [raw air, {analysis['pressure_aircal_smooth_s']} s aircal]"))
                if analysis["temperature_corrected"]:
                    scaled_to.append(f"T from {T_GAS_REFERENCE_C:.0f} C")
                    terms.append(
                        f"1+({analysis['temperature_col']}(K)-{T_GAS_REFERENCE_K:.2f})/"
                        f"{T_GAS_TARGET_K:.0f}")
                # "measurement ... before calibrating", not "calibrated values
                # scaled": since 2026-07-31 the correction lands on the input,
                # so it moves the cal means and the calibration with them. The
                # The filled cal-mean dots stay on the raw scale; the open red
                # dots show the same cal means after calibration.
                #
                # TWO lines: with both corrections and a smoothing window this
                # runs to ~150 characters, and a note wider than the Axes runs
                # off the right-hand edge (it is in_layout=False, so nothing
                # makes room for it). The widest other line here is ~85.
                notes.append(
                    f"measurement normalised to {' and '.join(scaled_to)} "
                    f"before calibrating; raw trace unchanged"
                )
                notes.append(
                    f"  ×{' ×'.join(terms)}; cal means shown are the "
                    f"raw-window means"
                )
        # Outside the has_masking gate on purpose: Ozone has no masking
        # settings at all, and this is the one thing removing data from its
        # figure, so it is the only line that would explain the red trace.
        if show_filtered:
            floor = GASES[self.current_gas]["valid_min"]
            n_below = int(rejected.sum())
            notes.append(
                f"red = filtered, blue = raw; {n_below} reading"
                f"{'' if n_below == 1 else 's'} below {floor:g} {unit} "
                f"removed (sensor fault, not a measurement)"
            )
        # Its own line, and outside the has_masking gate: a manual flag is the
        # one removal on this figure that no setting in the panel explains, so
        # it has to say so for every gas, floor or cal-bottle.
        n_flagged = int(flagged.sum())
        if n_flagged:
            notes.append(
                f"black x = {n_flagged} point{'' if n_flagged == 1 else 's'} "
                f"flagged by hand in "
                f"{len(self.flagged.get(self.current_gas, []))} region"
                f"{'' if len(self.flagged.get(self.current_gas, [])) == 1 else 's'}"
                f" — removed from the "
                f"{'filtered' if not has_masking else 'calibrated'} record"
            )
        notes_text = None
        if notes:
            notes_text = ax.text(
                0.01, 0.01, "\n".join(notes),
                transform=ax.transAxes, ha="left", va="bottom",
                color=MUTED_COLOR, fontsize=9,
                # Along the bottom the key sits over the cal dives and the cal
                # mean dots rather than over empty sky, so it needs a backing
                # box to stay readable. Semi-transparent, like the legend's
                # frame, so it dims the data underneath without hiding it.
                bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=2),
                # Kept out of constrained_layout's tight bbox; see the note
                # block in redraw_corr for what happens when it is left in.
                in_layout=False,
            )
            # Set BEFORE _text_height_frac below, deliberately: an invisible
            # Text reports a unit bbox, so hiding the notes also hands the
            # legend back the bottom of the axes. (Toggling at runtime does not
            # move the legend -- that waits for the next redraw, which is the
            # right trade for a toggle that must not rescale the view.)
            notes_text.set_visible(self.show_plot_notes)
        # Rebuilt on every draw, like the stats selectors: figure.clear() above
        # destroyed the previous one.
        self._notes_artist = notes_text

        # Placed after the notes, and confined to the axes *above* them.
        # `loc="best"` searches for a gap in the DATA only -- a Text artist is
        # invisible to that search -- so on a flight whose data leaves the
        # bottom-left clear it would sit straight on top of the notes. It does
        # honour bbox_to_anchor, though (`_find_best_position` anchors its
        # candidate boxes inside it), so handing it the region above the notes
        # keeps the placement automatic while making the collision impossible.
        # (x, y, width, HEIGHT) -- the height is what is left above the notes,
        # not the full axes, or the anchor box overhangs the top and the
        # legend is placed outside the plot.
        # Frame the default view on the FILTERED data. A single -2292 ppb
        # fault otherwise sets the y-range and squashes the whole real ozone
        # record into the top fifth of the axes -- which would make masking
        # the fault pointless on the one figure it most needed to help. The
        # raw trace is still drawn and simply runs off-scale; the note says
        # how many readings that is, and zooming out still reaches them.
        # Skipped when a view is being preserved, which reapplies its own
        # limits below.
        if show_filtered and old_main_view is None:
            lo, hi = filtered.min(), filtered.max()
            if pd.notna(lo) and pd.notna(hi) and hi > lo:
                pad = 0.05 * (hi - lo)
                ax.set_ylim(lo - pad, hi + pad)

        reserved = self._text_height_frac(ax, notes_text)
        anchor = (0.0, reserved, 1.0, 1.0 - reserved)
        # in_layout=False for the same reason the note block above carries
        # it: constrained_layout counts an in-axes legend in the Axes' tight
        # bbox, so a legend wider than the panel demands room the panel cannot
        # give -- which on a narrow window ends in "axes sizes collapsed to
        # zero" rather than in a legend that merely looks cramped.
        ax.legend(handles, labels, loc="best", fontsize=9, framealpha=0.9,
                  bbox_to_anchor=anchor,
                  bbox_transform=ax.transAxes).set_in_layout(False)

        ax_aux2 = None
        if ax_aux is not None:
            aux_col, aux_ylabel = aux_info
            ax_aux.set_facecolor("#fcfcfb")
            shade_intervals(ax_aux, df["datetime"], analysis["trimmed"],
                            WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax_aux, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax_aux, df["datetime"], analysis["pumps_off"],
                            PUMPS_EXCLUDE_COLOR, alpha=0.16)

            aux_unit = aux_ylabel.split("(")[-1].rstrip(")") if "(" in aux_ylabel else ""
            self._register_stats_trace("aux:left", f"{aux_col} (above, left)",
                                       ax_aux, df["datetime"], df[aux_col], aux_unit)

            aux_data = df[["datetime", aux_col]].dropna()
            # Faded only when the smoothed trace is drawn over it, the same
            # rule the calibrated overlay follows: the two are told apart by
            # hue, and a lone raw trace is not dimmed for nothing.
            smoothed_p = (analysis["pressure_smoothed"]
                          if aux_col == analysis["detector_col"] else None)
            aux_line, = ax_aux.plot(
                aux_data["datetime"], aux_data[aux_col], color=LINE_COLOR,
                linewidth=1.0, alpha=0.55 if smoothed_p is not None else 1.0)

            ax_aux.set_ylabel(aux_ylabel, color=TEXT_COLOR, fontsize=9)
            aux_title = self.other_column if self.aux_selection == "Other" else self.aux_selection
            ax_aux.set_title(aux_title, color=TEXT_COLOR, loc="left", fontsize=10)
            ax_aux.grid(True, color=GRID_COLOR, linewidth=0.6)
            for spine in ax_aux.spines.values():
                spine.set_color(AXIS_COLOR)
            ax_aux.tick_params(colors=MUTED_COLOR, labelsize=8, labelbottom=False)

            aux_handles = [aux_line]
            aux_labels = [aux_col]
            # On the same axis, deliberately: the point of drawing it is to
            # see how the mean sits against the reading it replaces, which a
            # second scale would make impossible to judge.
            if smoothed_p is not None:
                self._register_stats_trace(
                    "aux:smoothed", f"{aux_col} smoothed (above, left)",
                    ax_aux, df["datetime"], smoothed_p, aux_unit)
                smooth_data = pd.DataFrame(
                    {"datetime": df["datetime"], "p": smoothed_p}).dropna()
                smooth_line, = ax_aux.plot(
                    smooth_data["datetime"], smooth_data["p"],
                    color=CALIBRATED_COLOR, linewidth=1.0)
                aux_handles.append(smooth_line)
                if analysis["pressure_smooth_s"]:
                    aux_labels.append(
                        f"{analysis['pressure_smooth_s']} s air mean; "
                        f"{analysis['pressure_aircal_smooth_s']} s aircal mean")
                else:
                    aux_labels.append(
                        f"raw air; {analysis['pressure_aircal_smooth_s']} s aircal mean")
            if has_right_axis:
                ax_aux2 = ax_aux.twinx()
                self._register_stats_trace(
                    "aux:right", f"{self.right_axis_column} (above, right)",
                    ax_aux2, df["datetime"], df[self.right_axis_column], "")
                right_data = df[["datetime", self.right_axis_column]].dropna()
                right_line, = ax_aux2.plot(
                    right_data["datetime"], right_data[self.right_axis_column],
                    color=RIGHT_AXIS_COLOR, linewidth=1.0,
                )
                ax_aux2.set_ylabel(self.right_axis_column, color=RIGHT_AXIS_COLOR, fontsize=9)
                ax_aux2.spines["right"].set_color(RIGHT_AXIS_COLOR)
                ax_aux2.tick_params(axis="y", colors=RIGHT_AXIS_COLOR, labelsize=8)
                aux_handles.append(right_line)
                aux_labels.append(self.right_axis_column)

            # Drawn whenever there is more than one trace up there to tell
            # apart, which the right-hand axis is no longer the only way to
            # get: the smoothed pressure is a second line on the same axis.
            if len(aux_handles) > 1:
                ax_aux.legend(aux_handles, aux_labels, loc="upper right",
                              fontsize=8, framealpha=0.9).set_in_layout(False)

        # The new Axes were just built and auto-scaled to the full data
        # range -- reset the toolbar's view stack so Home returns to *this*
        # full-scale view, then (optionally) re-apply the pre-redraw zoom on
        # top of it, without pushing that onto the stack.
        self.main_pane.reset_nav()
        if old_main_view is not None:
            ax.set_xlim(old_main_view[0])
            ax.set_ylim(old_main_view[1])
        if old_aux_ylim is not None and ax_aux is not None:
            ax_aux.set_ylim(old_aux_ylim)
        if old_right_ylim is not None and ax_aux2 is not None:
            ax_aux2.set_ylim(old_right_ylim)

        self.ax = ax
        self.ax_aux = ax_aux
        self.ax_aux2 = ax_aux2
        self._had_aux_panel = has_aux_panel
        self._last_aux_key = aux_key
        self._last_right_axis_key = right_axis_key

        # Must happen on every draw: the Figure was cleared above, so any
        # selector from the previous draw is holding a destroyed Axes and
        # would silently stop responding. The aux panel gets its own selector,
        # so a box can be drawn in either.
        self.main_pane.set_stats_traces(
            [(key, t["label"]) for key, t in self._stats_traces.items()])
        self.main_pane.attach_stats_selectors([ax, ax_aux])

        self.canvas.draw()

    def _calibrated_cal_points(self, calibration):
        """Cal-mean points on the calibrated scale for the timeseries overlay.

        The exported calibrated record remains good ambient air only. These
        points are just a visual check of how the selected calibration maps the
        cal injections.
        """
        points = []
        bottles = calibration.get("bottles") or {}
        for t, state, closure, _ in calibration.get("residuals", []):
            assigned = (bottles.get(state) or {}).get("assigned")
            if assigned is None or pd.isna(closure):
                continue
            points.append((t, assigned + closure, state))
        return points

    def _register_stats_trace(self, key, label, axes, x, y, unit):
        """Record a plotted trace so the box-stats combo can offer it."""
        self._stats_traces[key] = {
            "label": label, "axes": axes, "x": x, "y": y, "unit": unit,
        }

    def _get_calibration(self):
        """The calibration for the gas on display. See _calibration_for."""
        return self._calibration_for(self.current_gas)

    def _calibration_for(self, gas_key):
        """The calibration for one gas, from its analysis and drift settings.

        Cached per gas and computed lazily, so a gas with no cal system, or a
        session that never opens the Calibration or Correlations tab, never
        pays for it. Invalidated alongside the analysis in refresh().
        """
        analysis = self._analysis_for(gas_key)
        if analysis is None:
            return None
        if gas_key in self._calibration:
            return self._calibration[gas_key]

        settings = self._settings_for(gas_key)
        if not analysis["has_masking"]:
            self._calibration[gas_key] = {
                "ok": False,
                "reason": f"{gas_key} is not run through the cal-bottle "
                          f"system, so there is nothing to calibrate against.",
            }
            return self._calibration[gas_key]
        if not analysis["cal_points"]:
            self._calibration[gas_key] = {
                "ok": False,
                "reason": (f"No cal events survive the current masking "
                           f"({analysis['warmup_minutes']} min warm-up, "
                           f"±{analysis['pressure_tol']:.2f} mbar)."),
            }
            return self._calibration[gas_key]

        self._calibration[gas_key] = calibrate_series(
            self.df, GASES[gas_key]["value_col"], analysis["cal_points"],
            self.cal_bottles, gas_key,
            model=settings["drift_model"],
            smooth_window=settings["drift_smooth_events"],
            # Only read under model="fixed slope"; 0 there means "use the
            # constant model's slope", which is what the spin box shows as
            # "auto" and what the seeding puts a number to.
            fixed_slope=settings.get("fixed_slope") or None,
            roster=self.cal_roster, flush_mask=analysis["post_cal_flush"],
            cal_mask=analysis["not_air"],
            # Blanks warm-up/bad-pressure rows from the *output* only. The
            # same mask separately fed cal_mean_points above, which is what
            # affects the calibration; this use cannot.
            exclude_mask=analysis["exclude_mask"],
            # The P/T-corrected measurement, and the factor that produced it.
            # The SAME series fed cal_mean_points above, so the nodes and the
            # record being calibrated are on one scale by construction. The
            # factor is only described and used for blanking here -- it has
            # already been applied.
            values=analysis["corrected"],
            correction_factor=analysis["correction_factor"],
        )
        return self._calibration[gas_key]

    def _pressure_column(self, gas_key):
        """The detector pressure column this gas's correction divides by, or
        None when the gas has no detector (Ozone) or the file has no such
        column.

        Per gas rather than the fixed `d1_P_mbars` the masking uses, because
        the correction is arithmetic on the measurement itself: CH4 comes off
        the second Aeris head, and correcting it by the first head's cell
        pressure would be meaningless.
        """
        detector = GASES[gas_key].get("detector")
        if detector is None or self.df is None:
            return None
        col = f"{detector}_P_mbars"
        return col if col in self.df.columns else None

    def _temperature_column(self, gas_key):
        """The detector cell temperature column this gas's T-40 correction is
        taken from, or None when the gas has no detector or the file has no
        such column. Per gas for the same reason as _pressure_column: CH4 comes
        off the second Aeris head, and its cell runs at a different temperature
        (40.8 C against d1's 42.3 C on the Jul 2026 flight).

        The column is in degrees CELSIUS. The correction uses its delta from
        40 C, so this returns the column name and nothing else.
        """
        detector = GASES[gas_key].get("detector")
        if detector is None or self.df is None:
            return None
        col = f"{detector}_T_gas"
        return col if col in self.df.columns else None

    def _uncertainty_for(self, gas_key):
        """(sigma Series, components dict) for one gas, cached per gas."""
        if gas_key not in self._uncertainty:
            result = self._calibration_for(gas_key)
            if result is None:
                return None, {}
            self._uncertainty[gas_key] = calibration_uncertainty(result)
        return self._uncertainty[gas_key]

    def redraw_cal(self, preserve_view=False):
        """Draw the Calibration tab."""
        result = self._get_calibration()
        if result is None:
            return
        gas = GASES[self.current_gas]
        unit = gas_unit(self.current_gas)

        # Only the response/coefficient/residual y-ranges are worth holding,
        # and only while they still mean the same thing -- a different gas or
        # drift model changes that entirely, so let those re-autoscale.
        # fixed_slope belongs in this key like the smoothing window does: it
        # changes every panel on the tab, so leaving it out would hold the old
        # y-limits against a calibration that had moved under them.
        cal_key = (self.current_gas, self.drift_model, self.drift_smooth_events,
                   self.fixed_slope)
        old_ylims = None
        if (preserve_view and self._cal_ax and self._last_cal_key == cal_key):
            old_ylims = [ax.get_ylim() for ax in self._cal_ax]
            old_xlim = self._cal_ax[0].get_xlim()

        axes = plot_calibration_panels(
            self.cal_pane.figure, result, self.current_gas, gas["ylabel"],
            self.df["datetime"], unit=unit,
        )

        self.cal_pane.reset_nav()
        if old_ylims is not None and len(old_ylims) == len(axes):
            for ax, ylim in zip(axes, old_ylims):
                ax.set_ylim(ylim)
            axes[0].set_xlim(old_xlim)

        self._cal_ax = axes
        self._last_cal_key = cal_key
        self.cal_pane.canvas.draw()

    def _corr_axis(self, gas_key):
        """What goes on one correlation axis, as
        (values, sigma, unit, qualifier, reason).

        Two kinds of axis, and the difference is deliberately visible rather
        than smoothed over:

        - A cal-bottle gas contributes its CALIBRATED series -- a
          tracer-tracer slope from uncalibrated counts would carry that
          detector's gain error straight into the slope -- with a 1σ from the
          calibration, and it is already masked down to good air.
        - Ozone has no cal bottles, so `calibrate_series` has nothing to work
          from. It contributes `oz_o3best`, the ozone instrument's own
          product, with no sigma and no masking of its own. It is not
          left out over that: its partner axis is blanked wherever that gas
          was in cal, flushing or masked, and the pairing is an intersection,
          so those rows drop anyway.

        `reason` is non-None only when a gas that *should* have a calibration
        has no usable one -- that is a failure to report, whereas Ozone having
        none is just what Ozone is.
        """
        gas = GASES[gas_key]
        unit = gas_unit(gas_key)
        if not gas.get("has_masking", True):
            # Filtered, matching the red trace on the timeseries: a -2292 ppb
            # fault is not a point on a tracer-tracer plot, it is an outlier
            # that would set the axis range and drag the fit on its own.
            values = self.df[gas["value_col"]].mask(self._removed_mask(gas_key))
            qual = gas["value_col"]
            if gas.get("valid_min") is not None:
                qual += f" ≥ {gas['valid_min']:g}"
            return values, None, unit, qual, None
        result = self._calibration_for(gas_key)
        if not (result or {}).get("ok"):
            return None, None, unit, None, (result or {}).get("reason", "no calibration")
        sigma = self._uncertainty_for(gas_key)[0] if self.corr_error_bars else None
        return result["calibrated"], sigma, unit, "calibrated", None

    def _corr_axis_flagged(self, gas_key):
        """What a MANUALLY FLAGGED row would have contributed to this axis.

        Flagging blanks a row, so a flagged point drops straight out of the
        scatter's pairing and disappears -- which would leave the outlier you
        just struck out invisible, and unreachable for a right-drag unflag.
        This recovers the value it would have had, so it can be drawn as a
        struck-out marker instead.

        Nothing is re-derived to do it: `calibrate_series` deliberately emits
        `cal_slope`/`cal_intercept` on *every* row, blanked ones included,
        exactly so a blanked row's calibrated value can be recomputed. The value
        they multiply must be the same P/T-corrected measurement the calibration
        used, not the raw column, or flagged markers move when corrections are
        on. A floor gas has no calibration to undo, so its raw column is already
        the answer -- but its below-floor faults stay out, since those are a
        sensor fault rather than something the user chose to remove.

        Returns a Series that is NaN everywhere except this gas's flagged rows.
        """
        flags = self._flag_mask(gas_key)
        gas = GASES[gas_key]
        raw = self.df[gas["value_col"]]
        if not gas.get("has_masking", True):
            return raw.mask(~flags | self._rejected_mask(gas_key))
        result = self._calibration_for(gas_key)
        if not (result or {}).get("ok"):
            return pd.Series(float("nan"), index=self.df.index)
        corrected = self._analysis_for(gas_key)["corrected"]
        return (result["slope"] * corrected + result["intercept"]).mask(~flags)

    def _median_sigmas(self, *gas_units):
        """[(gas, median 1σ, unit, percent)] for those of `gas_units` that
        have a calibration to propagate one from, in the order given.

        `percent` is the 1σ as a percentage of the HIGHER cal tank's assigned
        value, or None when that bottle has none. That denominator rather than
        the gas's own mean: ambient sits nearest the high bottle, so it is the
        one a reader is implicitly comparing against when they ask how good
        the number is (see the 1σ panel on the Calibration tab, which puts the
        whole curve on the same scale).

        A gas with has_masking=False (Ozone, H2O) contributes no entry rather
        than a zero -- no cal bottles means there is no uncertainty to report,
        the same rule _corr_axis applies to the error bars themselves. Shared
        by the figure note and the numbers panel so the two cannot disagree,
        and de-duplicated because both axes may be the same gas.

        Cheap despite looking expensive: both calibrations are already built
        by the time either caller runs, and _uncertainty_for is cached per gas.
        """
        out, seen = [], set()
        for gas, unit in gas_units:
            if gas in seen or not GASES[gas].get("has_masking", True):
                continue
            seen.add(gas)
            sigma, _ = self._uncertainty_for(gas)
            if sigma is None:
                continue
            median = sigma.median()
            if pd.isna(median):
                continue
            result = self._calibration_for(gas) or {}
            high = (result.get("bottles") or {}).get(result.get("high_state"), {})
            assigned = high.get("assigned")
            out.append((gas, median, unit,
                        100.0 * median / assigned if assigned else None))
        return out

    def _corr_cal_mean_points_for_gas(self, gas_key):
        """Cal means on the axis scale used by the correlation figure.

        These are display-only diagnostics. They are not part of the ambient
        data pairing, fit, flagging, tooltip search, or export.
        """
        if not GASES[gas_key].get("has_masking", True):
            return []
        analysis = self._analysis_for(gas_key)
        result = self._calibration_for(gas_key) or {}
        if analysis is None or not result.get("ok"):
            return []
        points = []
        for t, value, state, serial in analysis["cal_points"]:
            slope = interp_hold(
                self.df["datetime"], result["slope"], pd.Series([t])).iloc[0]
            intercept = interp_hold(
                self.df["datetime"], result["intercept"], pd.Series([t])).iloc[0]
            if pd.isna(slope) or pd.isna(intercept) or pd.isna(value):
                continue
            points.append({
                "time": pd.Timestamp(t),
                "state": state,
                "serial": serial,
                "value": float(slope * value + intercept),
            })
        return points

    def _corr_cal_mean_pairs(self, x_gas, y_gas):
        """Paired cal means for the correlation diagnostic overlay."""
        x_pts = self._corr_cal_mean_points_for_gas(x_gas)
        if not x_pts:
            return []
        if x_gas == y_gas:
            return [(p["value"], p["value"]) for p in x_pts]

        y_pts = self._corr_cal_mean_points_for_gas(y_gas)
        if not y_pts:
            return []
        used = set()
        pairs = []
        tolerance = pd.Timedelta(seconds=5)
        for xp in x_pts:
            best_i, best_dt = None, None
            for i, yp in enumerate(y_pts):
                if i in used or yp["state"] != xp["state"]:
                    continue
                dt = abs(yp["time"] - xp["time"])
                if best_dt is None or dt < best_dt:
                    best_i, best_dt = i, dt
            if best_i is not None and best_dt <= tolerance:
                used.add(best_i)
                pairs.append((xp["value"], y_pts[best_i]["value"]))
        return pairs

    def redraw_corr(self, preserve_view=False):
        """Draw the Correlations tab: one tracer against another.

        Calibrated wherever a calibration exists -- a tracer-tracer slope from
        uncalibrated counts would carry each detector's gain error into the
        slope, which is the number the plot exists to produce. Ozone is the
        exception it cannot be: no cal bottles, so no calibration to apply.
        See _corr_axis.
        """
        if self.df is None or not self.corr_x_gas or not self.corr_y_gas:
            return
        self._close_corr_tooltip_popup()
        x_gas, y_gas = self.corr_x_gas, self.corr_y_gas
        fig = self.corr_pane.figure

        corr_key = (x_gas, y_gas)
        old_view = None
        if preserve_view and self._corr_ax is not None and self._last_corr_key == corr_key:
            old_view = (self._corr_ax.get_xlim(), self._corr_ax.get_ylim())

        fig.clear()
        ax = fig.add_subplot(111)
        self._style_axes(ax)

        x_vals, x_sigma, x_unit, x_qual, x_reason = self._corr_axis(x_gas)
        y_vals, y_sigma, y_unit, y_qual, y_reason = self._corr_axis(y_gas)
        failed = [(gas, reason) for gas, reason in ((x_gas, x_reason), (y_gas, y_reason))
                  if reason]
        if failed:
            # in_layout=False for the same reason as the note block below --
            # and doubly so with wrap=True, which measures its line width
            # against the Axes it is sizing.
            ax.text(0.5, 0.5,
                    "\n\n".join(f"{gas}: {reason}" for gas, reason in failed),
                    transform=ax.transAxes, ha="center", va="center",
                    color=MUTED_COLOR, fontsize=10, wrap=True, in_layout=False)
            ax.set_xticks([]); ax.set_yticks([])
            self.corr_stats_label.setText("")
            self._corr_ax = None
            self._corr_plotted = None
            self._corr_flag_scatter = None
            self._last_corr_key = corr_key
            self.corr_pane.reset_nav()
            self.corr_pane.canvas.draw()
            return

        # A calibrated series is already good-air-only (its own masking, cal
        # periods and flush are blanked by calibrate_series); Ozone carries no
        # masking of its own. Either way the pairing is the intersection of
        # what each axis has a value for, which is also what makes a
        # calibrated partner impose its masking on the Ozone axis.
        keep = x_vals.notna() & y_vals.notna()

        # The z variable joins the pairing rule: a point with no value for it
        # has no color, and matplotlib would draw it in the colormap's "bad"
        # color (transparent) -- present in the fit but invisible on the
        # plot. Dropped instead, and counted in the note.
        z_vals, z_label = None, None
        if self.corr_color_by in CORR_COLOR_BY:
            _, z_col, z_label = CORR_COLOR_BY[self.corr_color_by]
            if z_col is None:
                # Dates are not numbers; convert once, and format the colorbar
                # ticks back to clock time below.
                z_vals = pd.Series(mdates.date2num(self.df["datetime"]),
                                   index=self.df.index)
            elif z_col in self.df.columns:
                z_vals = pd.to_numeric(self.df[z_col], errors="coerce")
                if self.corr_color_by == "ozone":
                    z_vals = z_vals.mask(self._removed_mask("Ozone"))
            if z_vals is not None:
                keep &= z_vals.notna()

        x, y = x_vals[keep], y_vals[keep]

        if self.corr_error_bars and len(x) and (x_sigma is not None or y_sigma is not None):
            # Drawn under the markers and thin: at flight scale the bars
            # overlap into a band, and a band that hides its own points would
            # misrepresent the density the plot is mostly about. An axis with
            # no calibration passes None -- there is nothing to propagate, and
            # a zero bar would claim a precision nobody established.
            ax.errorbar(x, y,
                        xerr=None if x_sigma is None else x_sigma[keep],
                        yerr=None if y_sigma is None else y_sigma[keep],
                        fmt="none", ecolor=CAL_SHADE_COLOR, elinewidth=0.6,
                        alpha=0.5, zorder=1)

        # s is an area in points^2; the control is a diameter, which is what
        # "marker size" means to anyone looking at the plot.
        if z_vals is None:
            ax.scatter(x, y, s=self.corr_marker_size ** 2, color=LINE_COLOR,
                       alpha=0.55, edgecolors="none", zorder=2)
        else:
            # Less transparent than the single-color case: at 0.55 the hues
            # blend with the white surface and with each other, which is the
            # one thing this encoding cannot afford.
            points = ax.scatter(x, y, s=self.corr_marker_size ** 2,
                                c=z_vals[keep], cmap=CORR_COLORMAP,
                                alpha=0.85, edgecolors="none", zorder=2)
            bar = fig.colorbar(points, ax=ax, pad=0.02)
            bar.set_label(z_label, color=TEXT_COLOR, fontsize=9)
            bar.ax.tick_params(colors=MUTED_COLOR, labelsize=8)
            bar.outline.set_edgecolor(AXIS_COLOR)
            if CORR_COLOR_BY[self.corr_color_by][1] is None:
                bar.ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

        # Manually flagged points, struck out where they would have plotted.
        # They are NaN in x_vals/y_vals -- flagging is what blanked them -- so
        # they would otherwise simply vanish, leaving the outlier you just
        # removed invisible and out of reach of a right-drag unflag. A point
        # is drawn when EITHER axis is flagged and both have a value to place
        # it at; the pair is what has a position on this figure, even though
        # only one gas may be the one struck out.
        fx = x_vals.combine_first(self._corr_axis_flagged(x_gas))
        fy = y_vals.combine_first(self._corr_axis_flagged(y_gas))
        flagged_here = ((self._flag_mask(x_gas) | self._flag_mask(y_gas))
                        & fx.notna() & fy.notna())
        if z_vals is not None:
            flagged_here &= z_vals.notna()
        self._corr_flag_scatter = None
        if flagged_here.any():
            self._corr_flag_scatter = ax.scatter(
                fx[flagged_here], fy[flagged_here], marker="x",
                s=max(28, self.corr_marker_size ** 2), linewidths=1.1,
                color=FLAGGED_COLOR, zorder=4,
                visible=not self.corr_hide_flagged)

        fit = linear_fit(x, y) if self.corr_fit else None
        if fit:
            xs = [x.min(), x.max()]
            ax.plot(xs, [fit["slope"] * v + fit["intercept"] for v in xs],
                    color=CALIBRATED_COLOR, linewidth=1.4, zorder=3)

        cal_overlay_n = 0
        if self.corr_show_cals:
            cal_pairs = self._corr_cal_mean_pairs(x_gas, y_gas)
            cal_overlay_n = len(cal_pairs)
            if cal_pairs:
                cal_x, cal_y = zip(*cal_pairs)
                ax.scatter(cal_x, cal_y, s=6 ** 2, color=CALIBRATED_COLOR,
                           edgecolors="none", zorder=5)

        # The qualifier rides on each axis label rather than only in the
        # title: with one calibrated axis and one raw one, a single title word
        # would have to lie about one of them.
        ax.set_xlabel(f"{x_gas} ({x_unit}), {x_qual}", color=TEXT_COLOR)
        ax.set_ylabel(f"{y_gas} ({y_unit}), {y_qual}", color=TEXT_COLOR)
        date_str = self.df["datetime"].iloc[0].strftime("%Y-%m-%d")
        quals = {x_qual, y_qual}
        title_qual = quals.pop() if len(quals) == 1 else "see axis labels"
        ax.set_title(f"UCATS-B {y_gas} vs {x_gas} ({title_qual}, {date_str})",
                     color=TEXT_COLOR, loc="left")
        ax.grid(True, color=GRID_COLOR, linewidth=0.8)

        # Extrapolated spans are not excluded -- they are real data, and
        # dropping them silently would be worse than saying how much of the
        # plot rests on a held-flat calibration. Only calibrated axes have the
        # flag at all.
        extrapolated = pd.Series(False, index=self.df.index)
        for gas in (x_gas, y_gas):
            result = self._calibration_for(gas)
            if (result or {}).get("ok"):
                extrapolated |= result["extrapolated"]
        extrapolated = extrapolated[keep]
        notes = [f"n = {len(x)} of {len(self.df)} rows (usable data in both tracers"
                 + (f", and in {CORR_COLOR_BY[self.corr_color_by][0].lower()})"
                    if z_vals is not None else ")")]
        # Said plainly, because the axis label alone is easy to skim past and
        # the consequence -- no gain correction on that axis, so the slope
        # inherits that instrument's scale error -- is the whole reason the
        # other axes are calibrated.
        raw_axes = [gas for gas, qual in ((x_gas, x_qual), (y_gas, y_qual))
                    if qual != "calibrated"]
        for gas in raw_axes:
            note = (f"{gas} has no cal bottles — plotted as recorded "
                    f"({GASES[gas]['value_col']}), not calibrated")
            n_rejected = int(self._rejected_mask(gas).sum())
            if n_rejected:
                note += (f"; {n_rejected} below "
                         f"{GASES[gas]['valid_min']:g} removed")
            n_flagged = int(self._flag_mask(gas).sum())
            if n_flagged:
                note += f"; {n_flagged} flagged by hand"
            notes.append(note)
        # One line for both axes when they agree, since the usual case is the
        # correction being on everywhere. `has_masking` first, like every
        # other note here that assumes a calibration exists.
        for key, target, what in (
                ("pressure_corrected", f"{D1_P_TARGET_MBARS:.0f} mbar", "pressure"),
                ("temperature_corrected", f"T from {T_GAS_REFERENCE_C:.0f} C", "temperature")):
            corrected = [gas for gas in dict.fromkeys((x_gas, y_gas))
                         if GASES[gas].get("has_masking", True)
                         and self._analysis_for(gas)[key]]
            if corrected:
                notes.append(f"{', '.join(corrected)} normalised to {target} "
                             f"before calibrating ({what} correction)")
        # The fit summary rides in this block rather than in a legend: a
        # legend has to sit somewhere, and on a scatter that fills one corner
        # it lands either on the data or on this text.
        if fit:
            notes.append(f"red line: OLS  {y_gas} = {fit['slope']:.4g}(±{fit['slope_err']:.2g})"
                         f"·{x_gas} {fit['intercept']:+.4g}    r = {fit['r']:.4f}")
        if self.corr_show_cals:
            notes.append(f"dark red dots = {cal_overlay_n} paired cal means "
                         f"(diagnostic only)")
        if len(x) and extrapolated.any():
            notes.append(f"{extrapolated.mean():.0%} of points fall where at least one "
                         f"calibration is extrapolated")
        # Worth saying on this figure specifically: post-cal flush points run
        # in a line from the tank's composition to the atmosphere's, which
        # looks exactly like a tracer-tracer correlation and drags the fit.
        no_flush = [gas for gas in (x_gas, y_gas)
                    if GASES[gas].get("has_masking", True)
                    and not self._analysis_for(gas)["flag_air_s"]]
        if no_flush:
            notes.append(f"Flag Air is 0 s for {', '.join(dict.fromkeys(no_flush))} — "
                         f"post-cal flush points are included")
        if self.corr_error_bars:
            note = ("error bars: 1σ from the calibration only "
                    "(assigned values + drift-model reproducibility)")
            if raw_axes:
                note += f"; none on {', '.join(raw_axes)}"
            notes.append(note)
        # in_layout=False, and it has to be: a Text is unclipped by default, so
        # constrained_layout counts it in the Axes' tight bbox and reserves
        # room for the part that hangs off the right-hand edge. But the note is
        # anchored at 0.01 of the AXES width, so shrinking the Axes to make
        # room moves the note left by less than the Axes lost -- there is no
        # width at which the two agree until the Axes is as wide as the note.
        # Each draw walks part of the way there (~20 px, decaying), so the plot
        # crept wider on every redraw that changed nothing: toggling "Hide
        # flagged points" seven or eight times visibly grew the x axis, until
        # it saturated. It showed up here because this note is the widest text
        # on any of the figures -- the Ozone line naming oz_o3best, its floor
        # and its hand-flagged count is ~1000 px on its own -- but the same
        # feedback is latent in every note block, so they all set it.
        ax.text(0.01, 0.99, "\n".join(notes), transform=ax.transAxes,
                ha="left", va="top", color=MUTED_COLOR, fontsize=9,
                in_layout=False)

        self._update_corr_stats(fit, x, y, x_gas, y_gas, x_unit, y_unit)

        # What the flag tool resolves a box against. Recorded here rather than
        # recomputed there for the same reason _register_stats_trace exists:
        # the plotted set is the intersection of two masked series and the
        # z-variable's own pairing rule, and rebuilding that from the outside
        # would be a second implementation free to drift from this one.
        x_tip_sigma = (self._uncertainty_for(x_gas)[0]
                       if GASES[x_gas].get("has_masking", True) else None)
        y_tip_sigma = (self._uncertainty_for(y_gas)[0]
                       if GASES[y_gas].get("has_masking", True) else None)
        self._corr_plotted = {"keep": keep, "x": fx, "y": fy,
                              "flagged": flagged_here,
                              "x_gas": x_gas, "y_gas": y_gas,
                              "x_sigma": x_tip_sigma, "y_sigma": y_tip_sigma,
                              "x_unit": x_unit, "y_unit": y_unit}
        self.corr_pane.attach_stats_selectors([ax])

        # Frame a rescale on what is actually visible. A hidden artist keeps
        # its data limits, so the autoscale above still reaches out to flagged
        # points the user cannot see -- on Ozone that is the difference between
        # a 0-3663 ppb axis and a 0-900 one, for markers that are not drawn.
        #
        # Gated on `old_view is None`, i.e. on a full rescale having been asked
        # for, which is what a tracer change does (on_corr_gas_changed /
        # on_corr_swap_axes pass preserve_view=False). Every other route here
        # preserves the view, so this cannot reframe a plot the user has zoomed
        # in on -- and the Hide toggle itself never redraws at all, so toggling
        # flagged points on and off at a zoomed scale keeps working exactly as
        # before. Applied BEFORE reset_nav so this becomes the view Home
        # returns to, rather than a range sitting on top of a nav base that
        # still spans the hidden markers.
        rescale_to_visible = (old_view is None and self.corr_hide_flagged
                              and flagged_here.any())
        if rescale_to_visible:
            limits = self._corr_home_limits()
            if limits:
                ax.set_xlim(limits[0])
                ax.set_ylim(limits[1])

        self.corr_pane.reset_nav()
        if old_view is not None:
            ax.set_xlim(old_view[0])
            ax.set_ylim(old_view[1])
        self._corr_ax = ax
        self._last_corr_key = corr_key
        # After reset_nav, so Home frames what is actually drawn. Still needed
        # when the rescale above ran -- reset_nav captured the right view, but
        # the override is also what keeps Home correct once the user zooms and
        # then toggles the markers back off without a redraw.
        if self.corr_hide_flagged and flagged_here.any():
            self._retarget_corr_home()
        # Last, because the Hide toggle's enabled state keys off the marker
        # artist this draw just created (or didn't).
        self._update_corr_flag_readout()
        self.corr_pane.canvas.draw()

    def _text_height_frac(self, ax, text, cap=0.5):
        """How tall `text` is as a fraction of `ax`, for reserving space.

        Measured rather than assumed: the note block runs from one line to
        five depending on which masks are active, and the same text is a
        different fraction of a resized window. `get_renderer()` is enough --
        no full draw is needed to measure a Text. Falls back to a generous
        fixed strip if the backend won't give one up, since guessing too big
        only nudges the legend, while guessing too small puts it on the text.
        """
        if text is None:
            return 0.0
        try:
            renderer = self.figure.canvas.get_renderer()
            height = (text.get_window_extent(renderer)
                      .transformed(ax.transAxes.inverted()).height)
        except (AttributeError, RuntimeError, ValueError):
            return 0.2
        return min(height + 0.02, cap)

    @staticmethod
    def _style_axes(ax):
        """The muted spine/tick styling the other panels get inline."""
        for spine in ax.spines.values():
            spine.set_color(AXIS_COLOR)
        ax.tick_params(colors=MUTED_COLOR)

    def _update_corr_stats(self, fit, x, y, x_gas, y_gas, x_unit, y_unit):
        """Numbers panel under the correlation controls. Outside the Figure,
        for the same reason the box-stats readout is: it survives a redraw and
        can be read while the plot is zoomed somewhere else.

        The fit block is included only when the fit is switched on. A single
        straight line through a tracer-tracer plot with real structure in it
        describes almost none of that structure, so its slope and r are more
        likely to mislead than inform -- hence off by default, and hidden
        here rather than left sitting there looking authoritative.
        """
        if not len(x):
            self.corr_stats_label.setText("No overlapping points.")
            return
        head = [f"n      {len(x)}"]
        if self.corr_fit and fit:
            head += [
                f"slope  {fit['slope']:.5g} ± {fit['slope_err']:.3g}",
                f"       {y_unit} per {x_unit}",
                f"icept  {fit['intercept']:.6g} {y_unit}",
                f"r      {fit['r']:.5f}",
            ]
        elif self.corr_fit:
            head.append("(too few points to fit)")
        # Not gated on the error-bar toggle, and stated as ± so it cannot be
        # read as the mean ± std line above it: that one is the spread of the
        # atmosphere across the flight, this one is how well the number is
        # known, and they differ by more than an order of magnitude.
        lines = [f"{gas} ±{value:.3g} {unit}" + (f"  ({pct:.2g}%)" if pct else "")
                 for gas, value, unit, pct in self._median_sigmas((x_gas, x_unit),
                                                                  (y_gas, y_unit))]
        sigma_note = ("\nmedian 1σ  " + "\n           ".join(lines)) if lines else ""
        self.corr_stats_label.setText(
            "\n".join(head) + "\n"
            f"{x_gas:6s} {x.mean():.4g} ± {x.std():.3g} {x_unit}\n"
            f"{y_gas:6s} {y.mean():.4g} ± {y.std():.3g} {y_unit}"
            + sigma_note
        )


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None

    app = QApplication(sys.argv)
    window = UcatsbGui(csv_path)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
