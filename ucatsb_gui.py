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
import sys
from pathlib import Path

import pandas as pd
import yaml
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QGroupBox, QComboBox, QDoubleSpinBox, QSpinBox, QLabel,
    QButtonGroup, QRadioButton, QPushButton, QFileDialog, QMessageBox,
    QTabWidget, QCheckBox, QAction, QStackedWidget, QMenu,
    QDialog, QDialogButtonBox, QLineEdit, QPlainTextEdit, QScrollArea,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import matplotlib.dates as mdates

from ucatsb_analysis import (
    drop_presync_rows, find_intervals, merge_close_intervals,
    shade_intervals, cal_mean_points, load_cal_roster, load_cal_assignment,
    select_cal_bottles,
    most_common_serial, mean_std_label, calibrate_series, post_cal_flush_mask,
    cal_switch_mask, below_floor_mask, O3_VALID_MIN_PPB, H2O_VALID_MIN_PPM,
    box_stats, calibration_uncertainty, linear_fit,
    plot_calibration_panels,
    export_companion_csv, export_icartt, icartt_filename, icartt_time_base,
    DEFAULT_ICARTT_META,
    CALS_YAML_PATH, CAL_DRIFT_MODELS, CAL_DEFAULT_SMOOTH_EVENTS,
    CAL_MERGE_GAP_S, D1_P_TARGET_MBARS,
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
    "Ozone": {"value_col": "oz_o3best", "ylabel": "O3 (ppb)", "title": "UCATS-B O3 timeseries", "detector": None, "has_masking": False, "valid_min": O3_VALID_MIN_PPB, "short": "O3", "standard_name": "Gas_O3_InSitu_S_AVMR", "long_name": "Ozone mole fraction"},
    # Water vapour, from its own instrument (`w_*` columns) -- like Ozone it
    # has no cal bottles, so has_masking=False and it is plotted as recorded.
    "H2O": {"value_col": "w_H2Obest", "ylabel": "H2O (ppm)", "title": "UCATS-B H2O timeseries", "detector": None, "has_masking": False, "valid_min": H2O_VALID_MIN_PPM, "standard_name": "Met_H2OMF_InSitu_None", "long_name": "Water vapour mole fraction"},
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
    "pressure_tol_mbar": 10.0,
    "flag_air_s": 0,
    "cal1_window_s": [-15, -1],
    "cal2_window_s": [-15, -1],
    "drift_model": CAL_DRIFT_MODELS[0],
    "drift_smooth_events": CAL_DEFAULT_SMOOTH_EVENTS,
}


def flight_config_path(csv_path: Path) -> Path:
    """Where a dataset's own settings live: <dataset>_conf.yaml, beside the
    CSV. Per-flight rather than global because the right warm-up, pressure
    tolerance, cal windows and -- above all -- cal tanks are properties of the
    flight, and re-deriving them every time a file is reopened loses work."""
    return Path(csv_path).with_name(f"{Path(csv_path).stem}_conf.yaml")


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
        self.toolbar = NavigationToolbar(self.canvas, self)

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

    def reset_nav(self):
        """Point the toolbar's Home at the newly-built full-scale view; its
        nav stack otherwise still references the just-destroyed Axes."""
        self.toolbar.update()
        self.toolbar.push_current()

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
            # enlarge ax.dataLim to include (0, 0). Restoring dataLim is not
            # enough on its own: adding them has already triggered an autoscale
            # off the polluted limits and nothing recomputes the view, hence
            # the explicit autoscale_view() below -- which is safe because it
            # only touches an axis whose autoscale is still on, leaving a
            # preserved view or an explicitly framed y-range exactly as it was.
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
        self.corr_marker_size = 4
        self.corr_error_bars = False
        # None = single-color points; otherwise a key into CORR_COLOR_BY.
        self.corr_color_by = None
        # Off by default: a straight line through a tracer-tracer plot with
        # real structure in it describes almost none of that structure.
        self.corr_fit = False
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

        # The control panel is a stack, not one panel: everything in the
        # main panel is per-gas and the Correlations tab is inherently about
        # two gases at once, so its controls replace them rather than sitting
        # alongside and contradicting them. Tabs that share the per-gas
        # controls (Timeseries, Calibration, Cal Tanks) still share one panel.
        self.controls_stack = QStackedWidget()
        self.controls_stack.addWidget(self._build_controls())
        self.controls_stack.addWidget(self._build_corr_controls())
        self.controls_stack.setFixedWidth(300)
        layout.addWidget(self.controls_stack, 0)
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
        self.file_label.setText(csv_path.name)
        self.file_label.setToolTip(str(csv_path))
        # A flight whose schema predates j_pumps can't offer the filter. The
        # explicit disable survives mask_box's setEnabled(has_masking), since
        # Qt restores a child's own enabled state when its parent comes back.
        self.pumps_check.setEnabled("j_pumps" in df.columns)
        if "j_pumps" not in df.columns:
            self.pumps_check.setChecked(False)

        self.corr_file_label.setText(csv_path.name)
        self.corr_file_label.setToolTip(str(csv_path))
        self._populate_corr_combos()
        self._populate_corr_color_combo()

        was_initializing = self._initializing
        self._initializing = True

        self.gas_combo.blockSignals(True)
        self.gas_combo.clear()
        self.gas_combo.addItems(available_gases.keys())
        self.gas_combo.blockSignals(False)

        for name, rb in self.aux_radios.items():
            rb.setChecked(name == "No Figure")
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

        Any `<dataset stem>*.yaml` beside the CSV: Save lets the name be
        changed freely, so `..._conf.yaml`, `..._tight_conf.yaml` and
        `..._v2.yaml` are all the same dataset's configs, while another
        flight's are excluded by the stem.
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
        vbox = QVBoxLayout(panel)

        self.corr_file_label = QLabel("No file loaded")
        self.corr_file_label.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.corr_file_label)

        axes_box = QGroupBox("Tracers")
        axes_form = QFormLayout(axes_box)
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

        self.corr_note = QLabel(
            "Calibrated data wherever there is a calibration. A point needs a "
            "value in <i>both</i> tracers, so each gas's own masking, cal "
            "periods and post-cal flush all remove points. Ozone has no cal "
            "bottles, so it is plotted as recorded (<tt>oz_o3best</tt>) and "
            "gets no error bars. Settings for each gas stay on the other tabs."
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
        always works; `oz_p` only if the flight's schema has it."""
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
            return
        if not self._confirm_discard("loading another dataset"):
            return
        self._try_load(path, forget_on_failure=True)

    def _build_controls(self):
        panel = QWidget()
        panel.setFixedWidth(300)
        vbox = QVBoxLayout(panel)

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
        self.file_label = QLabel("No file loaded")
        self.file_label.setStyleSheet(f"color: {MUTED_COLOR};")
        vbox.addWidget(self.file_label)

        gas_box = QGroupBox("Gas")
        gas_layout = QVBoxLayout(gas_box)
        self.gas_combo = QComboBox()
        self.gas_combo.addItems(self.available_gases.keys())
        self.gas_combo.currentTextChanged.connect(self.on_gas_changed)
        gas_layout.addWidget(self.gas_combo)
        vbox.addWidget(gas_box)

        aux_box = QGroupBox("Trace Above")
        aux_layout = QVBoxLayout(aux_box)
        self.aux_group = QButtonGroup(aux_box)
        self.aux_radios = {}
        for name in AUX_OPTIONS:
            rb = QRadioButton(name)
            if name == "No Figure":
                rb.setChecked(True)
            rb.toggled.connect(self.on_aux_changed)
            self.aux_group.addButton(rb)
            aux_layout.addWidget(rb)
            self.aux_radios[name] = rb

        self.other_combo = QComboBox()
        self.other_combo.addItems(self.other_columns)
        self.other_combo.setEnabled(False)
        self.other_combo.currentTextChanged.connect(self.on_other_changed)
        aux_layout.addWidget(self.other_combo)

        aux_layout.addWidget(QLabel("Right axis:"))
        self.right_axis_combo = QComboBox()
        self.right_axis_combo.addItem("(none)")
        self.right_axis_combo.addItems(self.other_columns)
        self.right_axis_combo.setEnabled(False)
        self.right_axis_combo.currentTextChanged.connect(self.on_right_axis_changed)
        aux_layout.addWidget(self.right_axis_combo)

        vbox.addWidget(aux_box)

        self.mask_box = QGroupBox("Data Masking")
        mask_form = QFormLayout(self.mask_box)

        self.warmup_spin = QSpinBox()
        self.warmup_spin.setRange(0, 120)
        self.warmup_spin.setSingleStep(1)
        self.warmup_spin.setSuffix(" min")
        self.warmup_spin.valueChanged.connect(self.on_control_changed)
        mask_form.addRow("Warm-up exclude:", self.warmup_spin)

        # The descent is the busiest part of a flight and the least like the
        # rest of it, so trimming the tail is as routine as trimming the
        # warm-up -- same treatment, same orange band, one shared note line.
        self.end_flight_spin = QSpinBox()
        self.end_flight_spin.setRange(0, 120)
        self.end_flight_spin.setSuffix(" min")
        self.end_flight_spin.setToolTip(
            "Exclude this many minutes at the END of the record, the mirror\n"
            "of the warm-up exclusion at the start -- descent, landing and\n"
            "whatever happens on the ground afterwards. Like the warm-up it\n"
            "reaches the cal means, not just the plot. 0 disables it."
        )
        self.end_flight_spin.valueChanged.connect(self.on_control_changed)
        mask_form.addRow("End-flight exclude:", self.end_flight_spin)

        self.pressure_tol_spin = QDoubleSpinBox()
        self.pressure_tol_spin.setRange(0.0, 10.0)
        self.pressure_tol_spin.setSingleStep(0.05)
        self.pressure_tol_spin.setDecimals(2)
        self.pressure_tol_spin.setSuffix(" mbar")
        self.pressure_tol_spin.setMinimumWidth(130)
        self.pressure_tol_spin.valueChanged.connect(self.on_control_changed)
        mask_form.addRow(f"Pressure tol\n(±{D1_P_TARGET_MBARS:.0f} mbar target):", self.pressure_tol_spin)

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
            "Copy warm-up, pressure tolerance, Flag Air and both cal mean\n"
            "windows from the current gas to every other calibrated gas\n"
            "(CO2/N2O/CH4). The drift model and smoothing window are left\n"
            "alone."
        )
        self.copy_mask_button.clicked.connect(self.on_copy_masking_to_all)
        mask_form.addRow(self.copy_mask_button)

        vbox.addWidget(self.mask_box)

        # Visual order swapped: the 100% bottle box appears above the 50% one.
        self.cal2_box, self.cal2_start_spin, self.cal2_end_spin = self._add_cal_window_box(
            vbox, "Cal 2 Mean Window"
        )
        self.cal1_box, self.cal1_start_spin, self.cal1_end_spin = self._add_cal_window_box(
            vbox, "Cal 1 Mean Window"
        )

        self.cal_box = QGroupBox("Calibration")
        cal_form = QFormLayout(self.cal_box)

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
        drift_row.addWidget(QLabel("Drift model:"))
        self.drift_combo = QComboBox()
        self.drift_combo.addItems(CAL_DRIFT_MODELS)
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
        cal_form.addRow(drift_row)

        # Session-only, deliberately not persisted and off by default: the
        # timeseries is documented as showing uncalibrated data, and a
        # remembered toggle would let the app start up showing calibrated
        # data with no visible reason why.
        self.calibrated_check = QCheckBox("Show calibrated on main plot")
        self.calibrated_check.toggled.connect(self.on_calibrated_toggled)
        cal_form.addRow(self.calibrated_check)

        # No export button here any more: exporting is now the Export tab's
        # job, where it can cover every gas at once. This panel is per-gas,
        # and a per-gas button was quietly the reason the old export could
        # only ever describe one of them.
        vbox.addWidget(self.cal_box)

        # Deliberately NOT in the setEnabled(has_masking) list in _select_gas:
        # Ozone is the gas this feature exists for, and it is precisely the
        # one with no masking settings to enable. Its own group box for the
        # same reason -- it is not a masking *setting*, it is a record of
        # points the user struck out by hand.
        self.flag_box = QGroupBox("Flagged Points")
        flag_form = QVBoxLayout(self.flag_box)
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

    def _add_cal_window_box(self, vbox, title):
        """Start/End on one row -- two form rows per box cost more vertical
        space than the control panel can spare on a laptop screen."""
        box = QGroupBox(title)
        box.setToolTip(self.CAL_WINDOW_HELP)
        row = QHBoxLayout(box)

        spins = []
        for label in ("Start:", "End:"):
            spin = QSpinBox()
            spin.setRange(-60, 60)
            spin.setSuffix(" s")
            spin.valueChanged.connect(self.on_control_changed)
            row.addWidget(QLabel(label))
            row.addWidget(spin, 1)
            spins.append(spin)

        vbox.addWidget(box)
        return box, spins[0], spins[1]

    def _cal_box_title(self, label, fallback):
        """Title a cal-window box as "<info> Cal (<serial>) <mole fraction>"
        (e.g. "50% Cal (CB09960) 206.51 ppm") using cals.yaml's info field
        and its assigned value for the active gas, if the serial was
        matched; otherwise fall back. `info` (the rough-percentage label) is
        optional -- not every tank in the roster has one -- so the mole
        fraction still shows up without it rather than losing the title
        entirely."""
        nominal = self.cal_bottles.get(label, {}) if label else {}
        if not nominal:
            return fallback
        info = nominal.get("info")
        title = f"{info} Cal ({label})" if info else f"Cal ({label})"
        value = nominal.get(self.current_gas)
        if value is not None:
            title += f" {value:g} {gas_unit(self.current_gas)}"
        return title

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
        # The Calibration tab's three panels each mean something different
        # (response deviation, coefficients, residuals), so a single box-stats
        # readout there would be ambiguous -- Timeseries only for now. Flagging
        # is Timeseries-only for a stronger reason: the cal panels plot derived
        # quantities (deviations, coefficients, residuals), not the rows a flag
        # would have to name.
        self.cal_pane.stats_action.setVisible(False)
        self.cal_pane.flag_action.setVisible(False)
        # Keep the historical attribute names bound to the timeseries pane so
        # redraw()'s existing body needs no changes.
        self.figure = self.main_pane.figure
        self.canvas = self.main_pane.canvas
        self.toolbar = self.main_pane.toolbar

        self.tanks_pane = self._build_cal_tanks_pane()
        self.corr_pane = PlotPane()
        # Same reason as the Calibration tab: the box-stats readout describes
        # one trace over a time span, which a tracer-tracer scatter is not --
        # and a box on a scatter names two gases' rows at once, not one gas's.
        self.corr_pane.stats_action.setVisible(False)
        self.corr_pane.flag_action.setVisible(False)

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
                f"Not saved yet — Save… offers {self.config_path.name}")

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
         "First part of the file name — the mission/instrument data product, "
         "e.g. SABRE-UCATSB. The ICARTT standard calls it common practice to "
         "prefix the project acronym, and requires the ID to match whatever "
         "the archiving data center has registered, so check it against the "
         "mission's data-management instructions rather than inventing one. "
         "Hyphens are kept; underscores separate the file name's own fields "
         "and are stripped."),
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
         "Appended to each species name, e.g. CO2_UCATSB; the 1-sigma "
         "variable becomes CO2e_UCATSB. Convention is the instrument or PI."),
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
            "Given on every row, so the calibrated value of a blanked row can "
            "be recomputed by anyone who wants to check it.")
        self.csv_comment_check = QCheckBox("Put the provenance notes in the CSV as # lines")
        self.csv_comment_check.setChecked(False)
        self.csv_comment_check.setToolTip(
            "Off by default: neither Excel nor Igor skips a leading comment "
            "block without being told to.\nWith it off the notes are written "
            "to a <name>_notes.txt beside the CSV instead.")
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
        self.flag_air_spin.setValue(settings["flag_air_s"])
        self.cal1_start_spin.setValue(settings["cal1_window_s"][0])
        self.cal1_end_spin.setValue(settings["cal1_window_s"][1])
        self.cal2_start_spin.setValue(settings["cal2_window_s"][0])
        self.cal2_end_spin.setValue(settings["cal2_window_s"][1])
        self.drift_combo.setCurrentText(settings["drift_model"])
        self.smooth_spin.setValue(settings["drift_smooth_events"])
        self.drift_model = settings["drift_model"]
        self.drift_smooth_events = settings["drift_smooth_events"]
        self.smooth_spin.setEnabled(settings["drift_model"] == "smooth")
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
            "flag_air_s": self.flag_air_spin.value(),
            "cal1_window_s": [self.cal1_start_spin.value(), self.cal1_end_spin.value()],
            "cal2_window_s": [self.cal2_start_spin.value(), self.cal2_end_spin.value()],
            "drift_model": self.drift_combo.currentText(),
            "drift_smooth_events": self.smooth_spin.value(),
        }

    def _select_gas(self, gas: str):
        """Sync gas-dependent controls (masking/cal boxes, per-gas
        settings) to `gas`. Shared by on_gas_changed and load_csv so a
        freshly loaded file ends up in the same state as if the user had
        picked this gas from the combo box themselves."""
        self.current_gas = gas
        has_masking = GASES[gas].get("has_masking", True)
        self.mask_box.setEnabled(has_masking)
        self.cal1_box.setEnabled(has_masking)
        self.cal2_box.setEnabled(has_masking)
        self.cal_box.setEnabled(has_masking)
        # flag_box is deliberately absent from that list -- see where it is
        # built. Its readout is per gas, so it does have to follow along.
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

    def on_control_changed(self):
        if self._loading or self._initializing or self.current_gas is None:
            return
        settings = self._controls_to_settings()
        self.config[self.current_gas] = settings
        self.drift_model = settings["drift_model"]
        self.drift_smooth_events = settings["drift_smooth_events"]
        self.smooth_spin.setEnabled(self.drift_model == "smooth")
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

    def _update_save_state(self):
        """Mark the Save button when there is something to save."""
        if not hasattr(self, "save_button"):
            return
        dirty = self._is_dirty()
        self.save_button.setText(self.SAVE_LABEL + (" •" if dirty else ""))
        name = self.config_path.name if self.config_path else "a new file"
        self.save_button.setToolTip(
            ("Unsaved changes. " if dirty else "No unsaved changes. ")
            + f"Save writes a config file — the dialog offers {name}, and any "
              "other name saves a second configuration of the same dataset."
        )

    def _mark_dirty(self):
        """Settings changed. Nothing is written -- this only updates the
        Save button. Replaces the old auto-save so that opening a saved
        analysis, experimenting and quitting leaves the file untouched."""
        self._update_save_state()

    def on_save_clicked(self):
        """Write the current settings to a config file of the user's choosing.

        Save-as every time, deliberately: the request was to keep several
        configurations per dataset, so the filename is always offered for
        editing rather than silently overwriting whatever was opened.
        """
        if self.df is None:
            return
        default = str(self.config_path or flight_config_path(self.csv_path))
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save configuration", default, "YAML Files (*.yaml);;All Files (*)")
        if not path_str:
            return False
        path = Path(path_str)
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
    # The drift model and its smoothing window are the only per-gas settings
    # left out -- they are a judgement about that gas's cal record (how noisy
    # its injections are), not a description of the flight.
    COPIED_SETTING_KEYS = ("warmup_min", "end_flight_min", "require_pumps",
                           "pressure_tol_mbar", "flag_air_s",
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
        # A different tracer is a different set of numbers on that axis, so
        # the old limits mean nothing -- rescale, as a gas change does on the
        # timeseries.
        self._refresh_corr(preserve_view=False)

    def on_corr_swap_axes(self):
        if self.corr_x_gas is None or self.corr_x_gas == self.corr_y_gas:
            return
        self.corr_x_gas, self.corr_y_gas = self.corr_y_gas, self.corr_x_gas
        loading = self._loading
        self._loading = True
        self.corr_x_combo.setCurrentText(self.corr_x_gas)
        self.corr_y_combo.setCurrentText(self.corr_y_gas)
        self._loading = loading
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
        self.corr_fit = self.corr_fit_check.isChecked()
        self._refresh_corr(preserve_view=True)

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
        """Put a Series computed on the trimmed analysis frame back onto the
        RAW file's row numbering.

        drop_presync_rows both removes leading rows and resets the index, so
        self.df row 0 is raw row `presync_dropped`. Everything the Export tab
        writes goes through here, because the companion CSV's entire promise
        is that its row N is the source file's row N -- an unshifted Series
        would line every gas up a few dozen rows early, silently and
        plausibly.
        """
        shifted = series.reset_index(drop=True)
        offset = self.presync_dropped
        if not offset:
            return shifted
        shifted.index = shifted.index + offset
        return shifted.reindex(range(offset + len(shifted)))

    def _export_gas_blocks(self):
        """One block per gas in this file, for either exporter.

        Built for every available gas rather than the one on display -- both
        products are whole-flight deliverables, and the per-gas control panel
        is the only reason the old export could describe just one. Each gas
        is analysed with its own saved settings via _analysis_for /
        _calibration_for, exactly as the Correlations tab does, so a gas
        nobody has selected this session still exports correctly.
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
                # No cal bottles: the physical floor is the only correction
                # there is, and calling the result "calibrated" would claim
                # something nobody established.
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
                block["final"] = self._to_raw_rows(result["calibrated"])
                block["final_kind"] = "calibrated"
                block["sigma"] = self._to_raw_rows(sigma)
                block["slope"] = self._to_raw_rows(result["slope"])
                block["intercept"] = self._to_raw_rows(result["intercept"])
                block["masks"] = {
                    "is_cal_period": self._to_raw_rows(result["in_cal"]),
                    "is_post_cal_flush": self._to_raw_rows(result["flushed"]),
                    "is_masked": self._to_raw_rows(result["excluded"]),
                    "is_extrapolated": self._to_raw_rows(result["extrapolated"]),
                    # A subset of is_masked (flags ride inside exclude_mask),
                    # written separately because it is the only removal a data
                    # user cannot reconstruct from the settings.
                    "is_flagged": self._to_raw_rows(analysis["flagged"]),
                }
            else:
                # No calibration is not a reason to export nothing for this
                # gas: the masks are still the answer to "which rows are air",
                # and they are what the raw column needs beside it.
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
                lines.append(f"{gas:<6} calibrated      {n:>7,} good ambient rows")
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
        if summary["notes_path"]:
            message += (f"\n\nThe provenance and column notes went to "
                        f"{summary['notes_path'].name} rather than into the CSV, "
                        f"so it opens cleanly in Excel and Igor.")
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

    def on_aux_changed(self, checked: bool):
        if not checked:
            return
        for name, rb in self.aux_radios.items():
            if rb.isChecked():
                self.aux_selection = name
                break
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

        bad_pressure = (df["d1_P_mbars"] - D1_P_TARGET_MBARS).abs() > pressure_tol
        bad_pressure = bad_pressure.fillna(False)

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

        cal_intervals, cal_points = [], []
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
            # Cal means are estimated from the raw data with these masks
            # applied -- a cal point can be dropped entirely if its window
            # has no valid data.
            cal_points = cal_mean_points(
                df, cal_intervals, gas["value_col"],
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
            "post_cal_flush": post_cal_flush,
            "has_masking": has_masking,
            "warmup_minutes": warmup_minutes,
            "end_flight_minutes": end_flight_minutes,
            "pressure_tol": pressure_tol,
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

        # cal0_pts/cal1_pts are guaranteed empty when has_masking is False
        # (cal_points == []), so cal0_label/cal1_label are never read below
        # without having been set here first.
        cal0_pts = [(t, v) for t, v, state, serial in cal_points if state == 0]
        cal1_pts = [(t, v) for t, v, state, serial in cal_points if state == 1]
        if has_masking:
            cal0_label = most_common_serial(cal_points, 0) or "Cal 1"
            cal1_label = most_common_serial(cal_points, 1) or "Cal 2"
            self.cal1_box.setTitle(self._cal_box_title(cal0_label, "Cal 1 Mean Window"))
            self.cal2_box.setTitle(self._cal_box_title(cal1_label, "Cal 2 Mean Window"))

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
                notes.append(
                    f"light red = excluded (d1_P_mbars outside {D1_P_TARGET_MBARS:.0f}±{pressure_tol:.2f} mbar)"
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
            )

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
        ax.legend(handles, labels, loc="best", fontsize=9, framealpha=0.9,
                  bbox_to_anchor=anchor, bbox_transform=ax.transAxes)

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
            aux_line, = ax_aux.plot(aux_data["datetime"], aux_data[aux_col], color=LINE_COLOR, linewidth=1.0)

            ax_aux.set_ylabel(aux_ylabel, color=TEXT_COLOR, fontsize=9)
            aux_title = self.other_column if self.aux_selection == "Other" else self.aux_selection
            ax_aux.set_title(aux_title, color=TEXT_COLOR, loc="left", fontsize=10)
            ax_aux.grid(True, color=GRID_COLOR, linewidth=0.6)
            for spine in ax_aux.spines.values():
                spine.set_color(AXIS_COLOR)
            ax_aux.tick_params(colors=MUTED_COLOR, labelsize=8, labelbottom=False)

            aux_handles = [aux_line]
            aux_labels = [aux_col]
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

            if has_right_axis:
                ax_aux.legend(aux_handles, aux_labels, loc="upper right", fontsize=8, framealpha=0.9)

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
            roster=self.cal_roster, flush_mask=analysis["post_cal_flush"],
            cal_mask=analysis["not_air"],
            # Blanks warm-up/bad-pressure rows from the *output* only. The
            # same mask separately fed cal_mean_points above, which is what
            # affects the calibration; this use cannot.
            exclude_mask=analysis["exclude_mask"],
        )
        return self._calibration[gas_key]

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
        cal_key = (self.current_gas, self.drift_model, self.drift_smooth_events)
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

    def _median_sigmas(self, *gas_units):
        """[(gas, median 1σ, unit)] for those of `gas_units` that have a
        calibration to propagate one from, in the order given.

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
            out.append((gas, median, unit))
        return out

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
            ax.text(0.5, 0.5,
                    "\n\n".join(f"{gas}: {reason}" for gas, reason in failed),
                    transform=ax.transAxes, ha="center", va="center",
                    color=MUTED_COLOR, fontsize=10, wrap=True)
            ax.set_xticks([]); ax.set_yticks([])
            self.corr_stats_label.setText("")
            self._corr_ax = None
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
                z_vals = self.df[z_col]
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

        fit = linear_fit(x, y) if self.corr_fit else None
        if fit:
            xs = [x.min(), x.max()]
            ax.plot(xs, [fit["slope"] * v + fit["intercept"] for v in xs],
                    color=CALIBRATED_COLOR, linewidth=1.4, zorder=3)

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
        # The fit summary rides in this block rather than in a legend: a
        # legend has to sit somewhere, and on a scatter that fills one corner
        # it lands either on the data or on this text.
        if fit:
            notes.append(f"red line: OLS  {y_gas} = {fit['slope']:.4g}(±{fit['slope_err']:.2g})"
                         f"·{x_gas} {fit['intercept']:+.4g}    r = {fit['r']:.4f}")
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
        # The size of that uncertainty, in the gas's own units, on the figure
        # rather than only in the numbers panel: at flight scale the bars
        # overlap into a band whose width cannot be read off the axes, and
        # this is the figure that gets saved and shown to someone else.
        # Independent of the error-bar toggle -- it is a property of the
        # calibration, not of whether the bars are drawn.
        sigmas = self._median_sigmas((x_gas, x_unit), (y_gas, y_unit))
        if sigmas:
            notes.append("median 1σ from the calibration: "
                         + ",  ".join(f"{gas} ±{value:.3g} {unit}"
                                      for gas, value, unit in sigmas))
        ax.text(0.01, 0.99, "\n".join(notes), transform=ax.transAxes,
                ha="left", va="top", color=MUTED_COLOR, fontsize=9)

        self._update_corr_stats(fit, x, y, x_gas, y_gas, x_unit, y_unit)

        self.corr_pane.reset_nav()
        if old_view is not None:
            ax.set_xlim(old_view[0])
            ax.set_ylim(old_view[1])
        self._corr_ax = ax
        self._last_corr_key = corr_key
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
        lines = [f"{gas} ±{value:.3g} {unit}"
                 for gas, value, unit in self._median_sigmas((x_gas, x_unit),
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
