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
import sys
from pathlib import Path

import pandas as pd
import yaml
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QGroupBox, QComboBox, QDoubleSpinBox, QSpinBox, QLabel,
    QButtonGroup, QRadioButton, QPushButton, QFileDialog, QMessageBox,
    QTabWidget, QCheckBox,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.dates as mdates

from plot_co2_timeseries import (
    drop_presync_rows, find_intervals, merge_close_intervals,
    shade_intervals, cal_mean_points, load_cal_bottles, load_cal_roster,
    most_common_serial, mean_std_label, calibrate_series,
    plot_calibration_panels, export_calibrated_csv,
    CALS_YAML_PATH, CAL_DRIFT_MODELS, CAL_DEFAULT_SMOOTH_EVENTS,
)

LINE_COLOR = "#2a78d6"
RIGHT_AXIS_COLOR = "#8e44ad"   # purple, distinct from the red/orange masking shades
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
    "oz_o3", "j_sol_cals", "j_sol_aircal",
] + [g["value_col"] for g in GASES.values()]

# Columns already exposed via a specific named control (gas traces, the
# named aux radio options) -- excluded from the "Other" catch-all combo box
# since picking them there would be redundant. This is a narrower set than
# REQUIRED_COLUMNS: j_sol_cals/j_sol_aircal are required for cal-interval
# detection but aren't plotted anywhere by name, so they stay selectable
# via "Other" (e.g. to sanity-check the raw digital flag against a trace).
NAMED_TRACE_COLUMNS = {
    "d1_P_mbars", "d2_P_mbars", "d1_T_gas", "d2_T_gas", "oz_o3",
} | {g["value_col"] for g in GASES.values()}

AUX_OPTIONS = ["No Figure", "Detector Pressure", "T_gas", "oz_o3", "Other"]


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
    if selection == "oz_o3":
        return "oz_o3", "O3 (ppb)"
    if selection == "Other":
        if other_column is None:
            return None
        return other_column, other_column
    return None

DEFAULT_CONFIG_PATH = Path(__file__).parent / "ucatsb_gui_config.yaml"

DEFAULT_GAS_SETTINGS = {
    "warmup_min": 30,
    "pressure_tol_mbar": 10.0,
    "cal1_window_s": [-15, -1],
    "cal2_window_s": [-15, -1],
    "drift_model": CAL_DRIFT_MODELS[0],
    "drift_smooth_events": CAL_DEFAULT_SMOOTH_EVENTS,
}


def load_config(path: Path) -> dict:
    """Load per-gas control settings, filling in defaults for anything
    missing. Gases with has_masking=False (Ozone) never get an entry --
    they don't use warm-up/pressure-tol/cal-window settings at all.
    """
    config = {
        gas: copy.deepcopy(DEFAULT_GAS_SETTINGS)
        for gas, info in GASES.items() if info.get("has_masking", True)
    }
    if path.exists():
        try:
            loaded = yaml.safe_load(path.read_text()) or {}
        except (yaml.YAMLError, OSError) as e:
            print(f"Warning: could not read {path}: {e}")
            loaded = {}
        for gas in config:
            if isinstance(loaded.get(gas), dict):
                config[gas].update(loaded[gas])
    return config


def save_config(path: Path, config: dict):
    path.write_text(yaml.safe_dump(config, sort_keys=False))


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
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

    def reset_nav(self):
        """Point the toolbar's Home at the newly-built full-scale view; its
        nav stack otherwise still references the just-destroyed Axes."""
        self.toolbar.update()
        self.toolbar.push_current()


class UcatsbGui(QMainWindow):
    def __init__(self, csv_path: Path = None, config_path: Path = DEFAULT_CONFIG_PATH):
        super().__init__()
        self.setWindowTitle("UCATS-B Viewer")
        self.resize(1300, 750)

        self.csv_path = None
        self.df = None
        self.available_gases = {}
        self.other_columns = []

        self.config_path = config_path
        self.config = load_config(config_path)
        self.current_gas = None
        self.aux_selection = "No Figure"
        self.other_column = None
        self.right_axis_column = None
        self._loading = False
        self._initializing = True
        self._analysis = None
        self._calibration = None
        self.drift_model = DEFAULT_GAS_SETTINGS["drift_model"]
        self.drift_smooth_events = DEFAULT_GAS_SETTINGS["drift_smooth_events"]
        self.show_calibrated = False
        self._dirty = {"main": True, "cal": True}
        self._preserve = {"main": False, "cal": False}
        self.cal_bottles = load_cal_bottles(CALS_YAML_PATH)
        self.cal_roster = load_cal_roster(CALS_YAML_PATH)
        self.ax = None
        self.ax_aux = None
        self.ax_aux2 = None
        self._had_aux_panel = None
        self._last_aux_key = None
        self._last_right_axis_key = None
        self._cal_ax = None
        self._last_cal_key = None

        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)

        layout.addWidget(self._build_controls(), 0)
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

        self.setWindowTitle(f"UCATS-B Viewer - {csv_path.name}")

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

        self.drift_combo = QComboBox()
        self.drift_combo.addItems(CAL_DRIFT_MODELS)
        self.drift_combo.currentTextChanged.connect(self.on_control_changed)
        cal_form.addRow("Drift model:", self.drift_combo)

        self.smooth_spin = QSpinBox()
        self.smooth_spin.setRange(2, 21)
        self.smooth_spin.setSuffix(" events")
        self.smooth_spin.valueChanged.connect(self.on_control_changed)
        cal_form.addRow("Smooth over:", self.smooth_spin)

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
        box = QGroupBox(title)
        box.setToolTip(self.CAL_WINDOW_HELP)
        form = QFormLayout(box)

        start_spin = QSpinBox()
        start_spin.setRange(-60, 60)
        start_spin.setSuffix(" s")
        start_spin.valueChanged.connect(self.on_control_changed)
        form.addRow("Start:", start_spin)

        end_spin = QSpinBox()
        end_spin.setRange(-60, 60)
        end_spin.setSuffix(" s")
        end_spin.valueChanged.connect(self.on_control_changed)
        form.addRow("End:", end_spin)

        vbox.addWidget(box)
        return box, start_spin, end_spin

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
        """Timeseries and Calibration views as tabs over the shared controls.

        The controls panel deliberately stays outside the tabs: every control
        affects both views, so moving them inside would mean duplicating the
        gas selector or making one tab depend on state invisible from the
        other.
        """
        self.main_pane = PlotPane()
        self.cal_pane = PlotPane()
        # Keep the historical attribute names bound to the timeseries pane so
        # redraw()'s existing body needs no changes.
        self.figure = self.main_pane.figure
        self.canvas = self.main_pane.canvas
        self.toolbar = self.main_pane.toolbar

        self.tabs = QTabWidget()
        self.tabs.addTab(self.main_pane, "Timeseries")
        self.tabs.addTab(self.cal_pane, "Calibration")
        # Connected after both addTab calls -- the first one fires
        # currentChanged before the second tab exists.
        self.tabs.currentChanged.connect(self.on_tab_changed)
        return self.tabs

    def _apply_settings_to_controls(self, settings: dict):
        """Populate the controls from a per-gas settings dict without
        triggering on_control_changed (and re-saving/redrawing) per field."""
        self._loading = True
        self.warmup_spin.setValue(settings["warmup_min"])
        self.pressure_tol_spin.setValue(settings["pressure_tol_mbar"])
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
        save_config(self.config_path, self.config)
        self.refresh(preserve_view=True)

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
        self._analysis = None
        self._calibration = None
        for pane in self._dirty:
            self._dirty[pane] = True
            # A requested full rescale has to survive until it is actually
            # honoured: changing gas (rescale) and then nudging a spinbox
            # (preserve) must still rescale the pane that never redrew.
            if not preserve_view:
                self._preserve[pane] = False
        self._draw_current_tab()

    def _draw_current_tab(self):
        """Draw whichever pane is showing, if it is dirty."""
        if self.df is None:
            return
        name = "cal" if self.tabs.currentIndex() == 1 else "main"
        if not self._dirty.get(name):
            return
        preserve = self._preserve.get(name, False)
        if name == "cal":
            self.redraw_cal(preserve_view=preserve)
        else:
            self.redraw(preserve_view=preserve)
        self._dirty[name] = False
        self._preserve[name] = True

    def on_tab_changed(self, index):
        if self._initializing:
            return
        self._draw_current_tab()

    def _get_analysis(self):
        """Masks, cal intervals and per-injection cal means for the current
        file, gas and control settings.

        Cached because more than one view reads the same numbers and they
        must not disagree; see refresh() for how the cache is invalidated.
        Returns None when no file is loaded.
        """
        if self.df is None:
            return None
        if self._analysis is not None:
            return self._analysis

        gas = GASES[self.current_gas]
        df = self.df
        has_masking = gas.get("has_masking", True)
        warmup_minutes = self.warmup_spin.value()
        pressure_tol = self.pressure_tol_spin.value()

        # bad_pressure/warmup are computed unconditionally -- even for a gas
        # with no masking of its own (Ozone), the aux panel can still show
        # an Aeris trace (e.g. Detector Pressure) that these masks apply to.
        # Only their use on the *main* plot (shading/notes/cal exclusion) is
        # gated on has_masking by the caller.
        cal = (df["j_sol_cals"].fillna(0).astype(bool)
               | df["j_sol_aircal"].fillna(0).astype(bool))

        bad_pressure = (df["d1_P_mbars"] - D1_P_TARGET_MBARS).abs() > pressure_tol
        bad_pressure = bad_pressure.fillna(False)

        warmup_end = df["datetime"].iloc[0] + pd.Timedelta(minutes=warmup_minutes)
        warmup = df["datetime"] < warmup_end

        cal_intervals, cal_points = [], []
        if has_masking:
            cal_intervals = merge_close_intervals(
                find_intervals(df["datetime"], cal), pd.Timedelta(seconds=CAL_MERGE_GAP_S)
            )
            # Cal means are estimated from the raw data with these masks
            # applied -- a cal point can be dropped entirely if its window
            # has no valid data.
            cal_points = cal_mean_points(
                df, cal_intervals, gas["value_col"],
                (self.cal1_start_spin.value(), self.cal1_end_spin.value()),
                (self.cal2_start_spin.value(), self.cal2_end_spin.value()),
                cal_bottles=self.cal_bottles, gas_key=self.current_gas,
                exclude_mask=bad_pressure | warmup,
            )

        self._analysis = {
            "cal": cal, "warmup": warmup, "bad_pressure": bad_pressure,
            "exclude_mask": bad_pressure | warmup,
            "cal_intervals": cal_intervals, "cal_points": cal_points,
            "has_masking": has_masking,
            "warmup_minutes": warmup_minutes, "pressure_tol": pressure_tol,
        }
        return self._analysis

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
        warmup = analysis["warmup"]
        bad_pressure = analysis["bad_pressure"]
        cal_points = analysis["cal_points"]

        if has_masking:
            shade_intervals(ax, df["datetime"], cal, CAL_SHADE_COLOR, alpha=0.3)
            shade_intervals(ax, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)

        # The calibrated trace is opt-in; the raw trace stays visible
        # underneath it (faded) so the correction being applied is always
        # legible rather than silently swapped in.
        calibration = self._get_calibration() if self.show_calibrated else None
        show_cal = bool(calibration and calibration.get("ok"))
        self.export_button.setEnabled(bool(
            (self._get_calibration() or {}).get("ok") if has_masking else False))

        plot_data = df[["datetime", value_col]].dropna()
        line, = ax.plot(plot_data["datetime"], plot_data[value_col], color=LINE_COLOR,
                        linewidth=1.2, alpha=0.35 if show_cal else 1.0)

        cal_line = None
        if show_cal:
            cal_df = pd.DataFrame({"datetime": df["datetime"],
                                   "v": calibration["calibrated"]}).dropna()
            cal_line, = ax.plot(cal_df["datetime"], cal_df["v"],
                                color=LINE_COLOR, linewidth=1.2)
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
                notes.append("gray = calibration/cal-air (j_sol_cals, j_sol_aircal)")
            if warmup.any():
                notes.append(f"orange = excluded (first {warmup_minutes} min warm-up)")
            if bad_pressure.any():
                notes.append(
                    f"light red = excluded (d1_P_mbars outside {D1_P_TARGET_MBARS:.0f}±{pressure_tol:.2f} mbar)"
                )
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

        self.canvas.draw()

    def _get_calibration(self):
        """The calibration for the current analysis and drift-model settings.

        Cached and computed lazily, so a gas with no cal system, or a session
        that never opens the Calibration tab, never pays for it. Invalidated
        alongside the analysis in refresh().
        """
        analysis = self._get_analysis()
        if analysis is None:
            return None
        if self._calibration is not None:
            return self._calibration

        if not analysis["has_masking"]:
            self._calibration = {
                "ok": False,
                "reason": f"{self.current_gas} is not run through the cal-bottle "
                          f"system, so there is nothing to calibrate against.",
            }
            return self._calibration
        if not analysis["cal_points"]:
            self._calibration = {
                "ok": False,
                "reason": (f"No cal events survive the current masking "
                           f"({analysis['warmup_minutes']} min warm-up, "
                           f"±{analysis['pressure_tol']:.2f} mbar)."),
            }
            return self._calibration

        self._calibration = calibrate_series(
            self.df, GASES[self.current_gas]["value_col"], analysis["cal_points"],
            self.cal_bottles, self.current_gas,
            model=self.drift_model, smooth_window=self.drift_smooth_events,
            roster=self.cal_roster,
        )
        return self._calibration

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


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None

    app = QApplication(sys.argv)
    window = UcatsbGui(csv_path)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
