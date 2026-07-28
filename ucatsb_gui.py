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
    QTabWidget, QCheckBox, QAction, QStackedWidget,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import matplotlib.dates as mdates

from plot_co2_timeseries import (
    drop_presync_rows, find_intervals, merge_close_intervals,
    shade_intervals, cal_mean_points, load_cal_roster, load_cal_assignment,
    select_cal_bottles,
    most_common_serial, mean_std_label, calibrate_series, post_cal_flush_mask,
    cal_switch_mask,
    box_stats, calibration_uncertainty, linear_fit,
    plot_calibration_panels, export_calibrated_csv,
    CALS_YAML_PATH, CAL_DRIFT_MODELS, CAL_DEFAULT_SMOOTH_EVENTS,
    POST_CAL_FLUSH_COLOR,
)

LINE_COLOR = "#2a78d6"
# The calibrated overlay gets its own colour rather than reusing LINE_COLOR:
# the raw trace stays blue when the overlay is on, so the two are told apart
# by hue, not by which one happens to be faded. Darker/more saturated than the
# 15%-alpha PRESSURE_EXCLUDE_COLOR band so it doesn't read as shading.
CALIBRATED_COLOR = "#c0392b"
RIGHT_AXIS_COLOR = "#8e44ad"   # purple, distinct from the red/orange masking shades
CAL_SHADE_COLOR = "#898781"
PRESSURE_EXCLUDE_COLOR = "#d03b3b"
WARMUP_EXCLUDE_COLOR = "#ffa64d"   # light orange
STATS_BOX_COLOR = "#111111"
CAL0_COLOR = "#eda100"   # golden
CAL1_COLOR = "#0d366b"   # dark blue
GRID_COLOR = "#e1e0d9"
AXIS_COLOR = "#c3c2b7"
TEXT_COLOR = "#0b0b0b"
MUTED_COLOR = "#52514e"

D1_P_TARGET_MBARS = 140.0
CAL_MERGE_GAP_S = 2   # bridge cal periods split by a single dropped-flag sample

GASES = {
    "CO2": {"value_col": "d1_CO2_ppm", "ylabel": "CO2 (ppm)", "title": "UCATS-B CO2 (uncalibrated) timeseries", "detector": "d1"},
    "N2O": {"value_col": "d1_N2O_ppb", "ylabel": "N2O (ppb)", "title": "UCATS-B N2O (uncalibrated) timeseries", "detector": "d1"},
    "CH4": {"value_col": "d2_CH4_ppb", "ylabel": "CH4 (ppb)", "title": "UCATS-B CH4 (uncalibrated) timeseries", "detector": "d2"},
    # Ozone comes from its own dedicated sensor, not an Aeris detector, and
    # isn't run through the cal-bottle system -- has_masking=False skips the
    # warm-up/pressure-tol/cal-window machinery entirely for this gas, and
    # detector=None disables the Detector Pressure/T_gas aux traces (there's
    # no matching column to route to). Kept last so it sorts to the bottom
    # of the Gas combo box.
    "Ozone": {"value_col": "oz_o3best", "ylabel": "O3 (ppb)", "title": "UCATS-B O3 timeseries", "detector": None, "has_masking": False},
}

REQUIRED_COLUMNS = [
    "datetime", "d1_P_mbars", "d2_P_mbars", "d1_T_gas", "d2_T_gas",
    "j_sol_cals", "j_sol_aircal",
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

DEFAULT_CONFIG_PATH = Path(__file__).parent / "ucatsb_gui_config.yaml"

DEFAULT_GAS_SETTINGS = {
    "warmup_min": 30,
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
    if not path.exists():
        return {}
    try:
        return yaml.safe_load(path.read_text()) or {}
    except (yaml.YAMLError, OSError) as e:
        print(f"Warning: could not read {path}: {e}")
        return {}


def load_config(path: Path, template: dict = None) -> dict:
    """Load per-gas control settings, filling in anything missing from
    `template` (the app-level config, for a flight being opened for the first
    time) or from DEFAULT_GAS_SETTINGS. Gases with has_masking=False (Ozone)
    never get an entry -- they don't use warm-up/pressure-tol/cal-window
    settings at all.

    Gas blocks only. The tank selection shares the same file but is kept out
    of this dict deliberately -- see load_cal_selection.
    """
    config = {
        gas: copy.deepcopy((template or {}).get(gas) or DEFAULT_GAS_SETTINGS)
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


def save_config(path: Path, config: dict, cal_selection: dict = None):
    """Write the per-gas blocks, plus the tank selection when one is given
    (i.e. for a flight's own conf file -- the app-level config stays free of
    tank choices, for the reason in load_cal_selection)."""
    doc = dict(config)
    if cal_selection:
        doc = {"cals": dict(cal_selection), **doc}
    path.write_text(yaml.safe_dump(doc, sort_keys=False))


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
        layout.addWidget(self.canvas)

        # on_box is set by the owner; the selectors themselves are rebuilt on
        # every draw (see attach_stats_selectors).
        self.on_box = None
        self.selectors = []
        # (axes, extents) of the live box -- only one exists at a time, in
        # whichever panel it was drawn in.
        self._box = None
        self._loading_traces = False

    def reset_nav(self):
        """Point the toolbar's Home at the newly-built full-scale view; its
        nav stack otherwise still references the just-destroyed Axes."""
        self.toolbar.update()
        self.toolbar.push_current()

    def _set_readout_visible(self, visible):
        self.stats_combo.setVisible(visible)
        self.stats_label.setVisible(visible)
        self.stats_copy_button.setVisible(visible)

    def _copy_stats(self):
        QApplication.clipboard().setText(self.stats_label.text())

    def _on_stats_toggled(self, checked):
        if checked:
            # Pan/zoom hold the canvas widgetlock, and _SelectorWidget.ignore()
            # drops every event while it is held -- the tool would look dead.
            # (The reverse needs no handling: clicking pan later just makes the
            # selector inert until pan is switched off again.)
            mode = str(self.toolbar.mode)
            if "pan" in mode:
                self.toolbar.pan()
            elif "zoom" in mode:
                self.toolbar.zoom()
        for i, sel in enumerate(self.selectors):
            sel.set_active(checked)
            sel.set_visible(checked and self._box is not None and i == self._box[0])
        self.canvas.draw_idle()
        self._set_readout_visible(checked and bool(self.stats_label.text()))

    def set_stats_text(self, text):
        self.stats_label.setText(text)
        self._set_readout_visible(bool(text) and self.stats_action.isChecked())

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
        active = self.stats_action.isChecked()
        for i, ax in enumerate(live):
            sel = RectangleSelector(
                ax, functools.partial(self._on_select, ax), useblit=True,
                interactive=True, button=[1], minspanx=3, minspany=3,
                spancoords="pixels",
                props=dict(facecolor="none", edgecolor=STATS_BOX_COLOR,
                           linewidth=1.4, linestyle="--"),
            )
            sel.set_active(active)
            if i == box_index and box_extents is not None:
                try:
                    sel.extents = box_extents
                    self._box = (i, box_extents)
                except Exception:
                    pass
            sel.set_visible(active and i == box_index and self._box is not None)
            self.selectors.append(sel)

    def _on_select(self, ax, eclick, erelease):
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
    def __init__(self, csv_path: Path = None, config_path: Path = DEFAULT_CONFIG_PATH):
        super().__init__()
        self.setWindowTitle("UCATS-B Viewer")
        self.resize(1300, 750)

        self.csv_path = None
        self.df = None
        self.available_gases = {}
        self.other_columns = []

        # Two config files, with different jobs: the app-level one is the
        # template a never-before-opened flight starts from, and config_path
        # is whatever is currently authoritative -- the flight's own
        # <dataset>_conf.yaml once a dataset is loaded.
        self.default_config_path = config_path
        self.config_path = config_path
        self.config = load_config(config_path)
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
        self._dirty = {"main": True, "cal": True, "corr": True}
        self._preserve = {"main": False, "cal": False, "corr": False}
        self.cal_roster = load_cal_roster(CALS_YAML_PATH)
        # cals.yaml's own pairing is the default only; a flight's conf file
        # overrides it, and the Cal Tanks tab edits it.
        self.default_cal_selection = load_cal_assignment(CALS_YAML_PATH)
        self.cal_selection = dict(self.default_cal_selection)
        self.cal_bottles = select_cal_bottles(self.cal_roster,
                                              self.cal_selection.values())
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
        self._corr_ax = None
        self._last_corr_key = None

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
        df = pd.read_csv(csv_path)
        df["datetime"] = pd.to_datetime(df["datetime"])
        df = drop_presync_rows(df)

        # Columns not already exposed via a named control -- offered through
        # the "Other" catch-all combo box. "datetime" is excluded too since
        # it's never a meaningful y-trace.
        other_columns = sorted(c for c in df.columns if c not in NAMED_TRACE_COLUMNS and c != "datetime")

        # Validation above passed -- safe to commit the new dataset now.
        self.csv_path = csv_path
        self.df = df
        self.available_gases = available_gases
        self.other_columns = other_columns
        self.current_gas = next(iter(available_gases))
        self.aux_selection = "No Figure"
        self.other_column = other_columns[0] if other_columns else None
        self.right_axis_column = None

        self._adopt_flight_config(csv_path)

        self.setWindowTitle(f"UCATS-B Viewer - {csv_path.name}")
        self.file_label.setText(csv_path.name)
        self.file_label.setToolTip(str(csv_path))
        self.corr_file_label.setText(csv_path.name)
        self.corr_file_label.setToolTip(str(csv_path))
        self._populate_corr_combos()

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

    def _adopt_flight_config(self, csv_path: Path):
        """Switch to this dataset's own <dataset>_conf.yaml, seeding it from
        the app-level config the first time the flight is opened, and write it
        out immediately so the file exists (with every gas in it) from the
        moment the data is loaded rather than only once a control is touched.

        If the dataset's directory can't be written to -- a read-only archive
        or a mounted share -- the app falls back to the app-level config
        rather than failing the load, and says so once. Losing per-flight
        persistence is worth less than losing the ability to open the file.
        """
        path = flight_config_path(csv_path)
        template = load_config(self.default_config_path)
        self.config = load_config(path, template=template)
        self.cal_selection = load_cal_selection(path, self.default_cal_selection)
        self._rebuild_cal_bottles()

        self.config_path = path
        try:
            save_config(path, self.config, self.cal_selection)
        except OSError as e:
            print(f"Warning: could not write {path}: {e}")
            self.config_path = self.default_config_path
        self._apply_cal_selection_to_controls()

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
        vbox.addWidget(style_box)

        self.corr_note = QLabel(
            "Calibrated data only. A point needs good air in <i>both</i> "
            "tracers, so each gas's own masking, cal periods and post-cal "
            "flush all remove points. Settings for each gas stay on the other "
            "tabs."
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

    def _populate_corr_combos(self):
        """Offer the calibratable gases in this file. Ozone is excluded: it
        has no cal bottles, so it has no calibrated series to plot, and this
        figure is calibrated data only."""
        gases = [gas for gas in self.available_gases
                 if GASES[gas].get("has_masking", True)]
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
        start_dir = str(self.csv_path.parent) if self.csv_path else str(Path.cwd())
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load UCATS-B CSV", start_dir, "CSV Files (*.csv);;All Files (*)"
        )
        if not path_str:
            return
        try:
            self.load_csv(Path(path_str))
        except (OSError, ValueError, pd.errors.ParserError) as e:
            QMessageBox.warning(self, "Load Data", f"Could not load {Path(path_str).name}:\n{e}")

    def _build_controls(self):
        panel = QWidget()
        panel.setFixedWidth(300)
        vbox = QVBoxLayout(panel)

        load_button = QPushButton("Load Data")
        load_button.clicked.connect(self.on_load_data_clicked)
        vbox.addWidget(load_button)

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
        mask_form.addRow("Flag Air\n(after cal):", self.flag_air_spin)

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

        self.export_button = QPushButton("Export calibrated CSV…")
        self.export_button.clicked.connect(self.on_export_clicked)
        self.export_button.setEnabled(False)
        cal_form.addRow(self.export_button)

        vbox.addWidget(self.cal_box)

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
            unit = GASES[self.current_gas]["ylabel"].split("(")[-1].rstrip(")")
            title += f" {value:g} {unit}"
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
        # The Calibration tab's three panels each mean something different
        # (response deviation, coefficients, residuals), so a single box-stats
        # readout there would be ambiguous -- Timeseries only for now.
        self.cal_pane.stats_action.setVisible(False)
        # Keep the historical attribute names bound to the timeseries pane so
        # redraw()'s existing body needs no changes.
        self.figure = self.main_pane.figure
        self.canvas = self.main_pane.canvas
        self.toolbar = self.main_pane.toolbar

        self.tanks_pane = self._build_cal_tanks_pane()
        self.corr_pane = PlotPane()
        # Same reason as the Calibration tab: the box-stats readout describes
        # one trace over a time span, which a tracer-tracer scatter is not.
        self.corr_pane.stats_action.setVisible(False)

        self.tabs = QTabWidget()
        self.tabs.addTab(self.main_pane, "Timeseries")
        self.tabs.addTab(self.cal_pane, "Calibration")
        self.tabs.addTab(self.corr_pane, "Correlations")
        self.tabs.addTab(self.tanks_pane, "Cal Tanks")
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
            unit = info["ylabel"].split("(")[-1].rstrip(")")
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
        else:
            self.tank_conf_label.setText(f"Saved in {self.config_path.name}")

    def _apply_settings_to_controls(self, settings: dict):
        """Populate the controls from a per-gas settings dict without
        triggering on_control_changed (and re-saving/redrawing) per field."""
        self._loading = True
        self.warmup_spin.setValue(settings["warmup_min"])
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
        if has_masking:
            self._apply_settings_to_controls(self.config[gas])

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
        self._save_settings()
        self.refresh(preserve_view=True)

    def _save_settings(self):
        """Persist to the flight's conf file, and to the app-level config as
        the template the next never-before-opened flight starts from. Both,
        because either alone loses something: the flight file alone means
        every new flight reverts to shipped defaults, and the app-level file
        alone is what the per-flight requirement exists to replace.

        Tank choices only ever reach the flight file (see load_cal_selection).
        """
        targets = [(self.config_path, self.cal_selection)]
        if self.config_path != self.default_config_path:
            targets.append((self.default_config_path, None))
        for path, cal_selection in targets:
            try:
                save_config(path, self.config, cal_selection)
            except OSError as e:
                print(f"Warning: could not write {path}: {e}")

    # What the button copies: the masking values and both cal mean windows.
    # The drift model and its smoothing window are the only per-gas settings
    # left out -- they are a judgement about that gas's cal record (how noisy
    # its injections are), not a description of the flight.
    COPIED_SETTING_KEYS = ("warmup_min", "pressure_tol_mbar", "flag_air_s",
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
        self._save_settings()
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

    def on_corr_style_changed(self):
        if self._loading or self._initializing:
            return
        self.corr_marker_size = self.corr_size_spin.value()
        self.corr_error_bars = self.corr_error_check.isChecked()
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
        worth, and a status bar would be invisible next to the button."""
        button.setText(message)
        QTimer.singleShot(msec, lambda: button.setText(restore))

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
        self._save_settings()
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
        self._save_settings()
        self.refresh(preserve_view=False)

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

    def on_export_clicked(self):
        result = self._get_calibration()
        if not (result and result.get("ok")):
            return
        default = ""
        if self.csv_path:
            default = str(self.csv_path.with_name(
                f"{self.csv_path.stem}_{self.current_gas}_calibrated.csv"))
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Export calibrated CSV", default, "CSV Files (*.csv);;All Files (*)")
        if not path_str:
            return
        try:
            export_calibrated_csv(
                Path(path_str), self.df, result,
                GASES[self.current_gas]["value_col"], self.current_gas,
                self.csv_path, self._controls_to_settings(),
                analysis=self._get_analysis(),
            )
        except OSError as e:
            QMessageBox.warning(self, "Export", f"Could not write {path_str}:\n{e}")
        else:
            QMessageBox.information(self, "Export", f"Wrote {path_str}")

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
        else:
            return
        if not self._dirty.get(name):
            return
        preserve = self._preserve.get(name, False)
        if name == "cal":
            self.redraw_cal(preserve_view=preserve)
        elif name == "corr":
            self.redraw_corr(preserve_view=preserve)
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
                exclude_mask=bad_pressure | warmup,
            )

        self._analysis[gas_key] = {
            "cal": cal, "cal_switch": cal_switch, "not_air": not_air,
            "warmup": warmup, "bad_pressure": bad_pressure,
            "exclude_mask": bad_pressure | warmup,
            "cal_intervals": cal_intervals, "cal_points": cal_points,
            "post_cal_flush": post_cal_flush,
            "has_masking": has_masking,
            "warmup_minutes": warmup_minutes, "pressure_tol": pressure_tol,
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
        bad_pressure = analysis["bad_pressure"]
        cal_points = analysis["cal_points"]
        post_cal_flush = analysis["post_cal_flush"]

        if has_masking:
            shade_intervals(ax, df["datetime"], not_air, CAL_SHADE_COLOR, alpha=0.3)
            shade_intervals(ax, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)
            # Shaded whether or not the calibrated trace is showing: the band
            # is how you find out these rows exist before turning it on.
            shade_intervals(ax, df["datetime"], post_cal_flush, POST_CAL_FLUSH_COLOR, alpha=0.22)

        # The calibrated trace is opt-in; the raw trace stays visible
        # underneath it (faded) so the correction being applied is always
        # legible rather than silently swapped in.
        calibration = self._get_calibration() if self.show_calibrated else None
        show_cal = bool(calibration and calibration.get("ok"))
        self.export_button.setEnabled(bool(
            (self._get_calibration() or {}).get("ok") if has_masking else False))

        # Registered as they are plotted rather than scraped back off the Axes
        # afterwards: the artists carry no units and no stable identity, and a
        # trace that is conditionally drawn would be easy to miss.
        self._stats_traces = {}
        gas_unit = gas["ylabel"].split("(")[-1].rstrip(")") if "(" in gas["ylabel"] else ""
        self._register_stats_trace(
            "main:raw", f"{self.current_gas} (raw)", ax, df["datetime"],
            df[value_col], gas_unit)

        plot_data = df[["datetime", value_col]].dropna()
        line, = ax.plot(plot_data["datetime"], plot_data[value_col], color=LINE_COLOR,
                        linewidth=1.2, alpha=0.55 if show_cal else 1.0)

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
                df["datetime"], calibrated, gas_unit)
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
        labels = ["raw (ambient)" if show_cal else f"{gas['ylabel']} (ambient)"]
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

        ax.legend(handles, labels, loc="lower right", fontsize=9, framealpha=0.9)

        notes = []
        if has_masking:
            if cal.any():
                notes.append("gray = calibration/cal-air (j_sol_cals, j_sol_aircal"
                             + (", + switch-over sample)"
                                if analysis["cal_switch"].any() else ")"))
            if warmup.any():
                notes.append(f"orange = excluded (first {warmup_minutes} min warm-up)")
            if bad_pressure.any():
                notes.append(
                    f"light red = excluded (d1_P_mbars outside {D1_P_TARGET_MBARS:.0f}±{pressure_tol:.2f} mbar)"
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
        if notes:
            ax.text(
                0.01, 0.98, "\n".join(notes),
                transform=ax.transAxes, ha="left", va="top",
                color=MUTED_COLOR, fontsize=9,
            )

        ax_aux2 = None
        if ax_aux is not None:
            aux_col, aux_ylabel = aux_info
            ax_aux.set_facecolor("#fcfcfb")
            shade_intervals(ax_aux, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax_aux, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)

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
        unit = gas["ylabel"].split("(")[-1].rstrip(")")

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

    def redraw_corr(self, preserve_view=False):
        """Draw the Correlations tab: one calibrated tracer against another.

        Calibrated only, by design -- a tracer-tracer slope from uncalibrated
        counts would carry each detector's gain error into the slope, which is
        the number the plot exists to produce.
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

        x_cal = self._calibration_for(x_gas)
        y_cal = self._calibration_for(y_gas)
        failed = [(gas, res) for gas, res in ((x_gas, x_cal), (y_gas, y_cal))
                  if not (res or {}).get("ok")]
        if failed:
            ax.text(0.5, 0.5,
                    "\n\n".join(f"{gas}: {res.get('reason', 'no calibration')}"
                                for gas, res in failed),
                    transform=ax.transAxes, ha="center", va="center",
                    color=MUTED_COLOR, fontsize=10, wrap=True)
            ax.set_xticks([]); ax.set_yticks([])
            self.corr_stats_label.setText("")
            self._corr_ax = None
            self._last_corr_key = corr_key
            self.corr_pane.reset_nav()
            self.corr_pane.canvas.draw()
            return

        # Both series are already "good air only" (each gas's own masking, cal
        # periods and flush are blanked by calibrate_series), so the pairing is
        # just the intersection of what survived on each axis.
        x, y = x_cal["calibrated"], y_cal["calibrated"]
        keep = x.notna() & y.notna()
        x, y = x[keep], y[keep]
        x_unit = GASES[x_gas]["ylabel"].split("(")[-1].rstrip(")")
        y_unit = GASES[y_gas]["ylabel"].split("(")[-1].rstrip(")")

        if self.corr_error_bars and len(x):
            x_sigma, _ = self._uncertainty_for(x_gas)
            y_sigma, _ = self._uncertainty_for(y_gas)
            # Drawn under the markers and thin: at flight scale the bars
            # overlap into a band, and a band that hides its own points would
            # misrepresent the density the plot is mostly about.
            ax.errorbar(x, y, xerr=x_sigma[keep], yerr=y_sigma[keep],
                        fmt="none", ecolor=CAL_SHADE_COLOR, elinewidth=0.6,
                        alpha=0.5, zorder=1)

        # s is an area in points^2; the control is a diameter, which is what
        # "marker size" means to anyone looking at the plot.
        ax.scatter(x, y, s=self.corr_marker_size ** 2, color=LINE_COLOR,
                   alpha=0.55, edgecolors="none", zorder=2)

        fit = linear_fit(x, y)
        if fit:
            xs = [x.min(), x.max()]
            ax.plot(xs, [fit["slope"] * v + fit["intercept"] for v in xs],
                    color=CALIBRATED_COLOR, linewidth=1.4, zorder=3)

        ax.set_xlabel(f"{x_gas} ({x_unit}), calibrated", color=TEXT_COLOR)
        ax.set_ylabel(f"{y_gas} ({y_unit}), calibrated", color=TEXT_COLOR)
        date_str = self.df["datetime"].iloc[0].strftime("%Y-%m-%d")
        ax.set_title(f"UCATS-B {y_gas} vs {x_gas} (calibrated, {date_str})",
                     color=TEXT_COLOR, loc="left")
        ax.grid(True, color=GRID_COLOR, linewidth=0.8)

        # Extrapolated spans are not excluded -- they are real data, and
        # dropping them silently would be worse than saying how much of the
        # plot rests on a held-flat calibration.
        extrapolated = (x_cal["extrapolated"] | y_cal["extrapolated"])[keep]
        notes = [f"n = {len(x)} of {len(self.df)} rows (good air in both tracers)"]
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
                    if not self._analysis_for(gas)["flag_air_s"]]
        if no_flush:
            notes.append(f"Flag Air is 0 s for {', '.join(dict.fromkeys(no_flush))} — "
                         f"post-cal flush points are included")
        if self.corr_error_bars:
            notes.append("error bars: 1σ from the calibration only "
                         "(assigned values + drift-model reproducibility)")
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

    @staticmethod
    def _style_axes(ax):
        """The muted spine/tick styling the other panels get inline."""
        for spine in ax.spines.values():
            spine.set_color(AXIS_COLOR)
        ax.tick_params(colors=MUTED_COLOR)

    def _update_corr_stats(self, fit, x, y, x_gas, y_gas, x_unit, y_unit):
        """Numbers panel under the correlation controls. Outside the Figure,
        for the same reason the box-stats readout is: it survives a redraw and
        can be read while the plot is zoomed somewhere else."""
        if not fit:
            self.corr_stats_label.setText("Not enough overlapping points to fit.")
            return
        sigma_note = ""
        if self.corr_error_bars:
            x_sigma, _ = self._uncertainty_for(x_gas)
            y_sigma, _ = self._uncertainty_for(y_gas)
            sigma_note = (f"\nmedian 1σ  {x_gas} {x_sigma.median():.3g} {x_unit}"
                          f"\n           {y_gas} {y_sigma.median():.3g} {y_unit}")
        self.corr_stats_label.setText(
            f"n      {fit['n']}\n"
            f"slope  {fit['slope']:.5g} ± {fit['slope_err']:.3g}\n"
            f"       {y_unit} per {x_unit}\n"
            f"icept  {fit['intercept']:.6g} {y_unit}\n"
            f"r      {fit['r']:.5f}\n"
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
