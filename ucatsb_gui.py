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
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.dates as mdates

from plot_co2_timeseries import (
    drop_presync_rows, find_intervals, merge_close_intervals,
    shade_intervals, cal_mean_points, load_cal_bottles, most_common_serial,
    mean_std_label, CALS_YAML_PATH,
)

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
CAL_MERGE_GAP_S = 2   # bridge cal periods split by a single dropped-flag sample

GASES = {
    "CO2": {"value_col": "d1_CO2_ppm", "ylabel": "CO2 (ppm)", "title": "UCATS-B CO2 (uncalibrated) timeseries", "detector": "d1"},
    "N2O": {"value_col": "d1_N2O_ppb", "ylabel": "N2O (ppb)", "title": "UCATS-B N2O (uncalibrated) timeseries", "detector": "d1"},
    "CH4": {"value_col": "d2_CH4_ppb", "ylabel": "CH4 (ppb)", "title": "UCATS-B CH4 (uncalibrated) timeseries", "detector": "d2"},
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
    ylabel); None if no column is selected yet.
    """
    detector = GASES[gas]["detector"]
    if selection == "Detector Pressure":
        col = f"{detector}_P_mbars"
        return col, f"{col} (mbar)"
    if selection == "T_gas":
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
}


def load_config(path: Path) -> dict:
    """Load per-gas control settings, filling in defaults for anything missing."""
    config = {gas: copy.deepcopy(DEFAULT_GAS_SETTINGS) for gas in GASES}
    if path.exists():
        try:
            loaded = yaml.safe_load(path.read_text()) or {}
        except (yaml.YAMLError, OSError) as e:
            print(f"Warning: could not read {path}: {e}")
            loaded = {}
        for gas in GASES:
            if isinstance(loaded.get(gas), dict):
                config[gas].update(loaded[gas])
    return config


def save_config(path: Path, config: dict):
    path.write_text(yaml.safe_dump(config, sort_keys=False))


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
        self._loading = False
        self._initializing = True
        self.cal_bottles = load_cal_bottles(CALS_YAML_PATH)
        self.ax = None
        self.ax_aux = None
        self._had_aux_panel = None
        self._last_aux_key = None

        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)

        layout.addWidget(self._build_controls(), 0)
        layout.addWidget(self._build_canvas(), 1)

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

        self.other_combo.blockSignals(True)
        self.other_combo.clear()
        self.other_combo.addItems(other_columns)
        self.other_combo.blockSignals(False)

        self._apply_settings_to_controls(self.config[self.current_gas])
        self._initializing = was_initializing
        if not self._initializing:
            self.redraw()

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

        vbox.addWidget(aux_box)

        mask_box = QGroupBox("Data Masking")
        mask_form = QFormLayout(mask_box)

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

        vbox.addWidget(mask_box)

        # Visual order swapped: the 100% bottle box appears above the 50% one.
        self.cal2_box, self.cal2_start_spin, self.cal2_end_spin = self._add_cal_window_box(
            vbox, "Cal 2 Mean Window"
        )
        self.cal1_box, self.cal1_start_spin, self.cal1_end_spin = self._add_cal_window_box(
            vbox, "Cal 1 Mean Window"
        )

        vbox.addWidget(QLabel(
            "Cal windows are relative to the last\n"
            "point in a cal period (Cal_p), e.g.\n"
            "-10 s to 2 s = [Cal_p-10s, Cal_p+2s].\n"
            "Positive values reach past Cal_p.\n"
            "Settings are saved per-gas."
        ))

        vbox.addStretch(1)
        return panel

    def _add_cal_window_box(self, vbox, title):
        box = QGroupBox(title)
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

    def _build_canvas(self):
        container = QWidget()
        layout = QVBoxLayout(container)
        self.figure = Figure(constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        return container

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
        self._loading = False

    def _controls_to_settings(self) -> dict:
        return {
            "warmup_min": self.warmup_spin.value(),
            "pressure_tol_mbar": self.pressure_tol_spin.value(),
            "cal1_window_s": [self.cal1_start_spin.value(), self.cal1_end_spin.value()],
            "cal2_window_s": [self.cal2_start_spin.value(), self.cal2_end_spin.value()],
        }

    def on_gas_changed(self, new_gas: str):
        if not new_gas:
            return
        self.current_gas = new_gas
        self._apply_settings_to_controls(self.config[new_gas])
        if self._initializing:
            return
        self.redraw()

    def on_control_changed(self):
        if self._loading or self._initializing or self.current_gas is None:
            return
        self.config[self.current_gas] = self._controls_to_settings()
        save_config(self.config_path, self.config)
        self.redraw(preserve_view=True)

    def on_aux_changed(self, checked: bool):
        if not checked:
            return
        for name, rb in self.aux_radios.items():
            if rb.isChecked():
                self.aux_selection = name
                break
        self.other_combo.setEnabled(self.aux_selection == "Other")
        if self._initializing:
            return
        self.redraw(preserve_view=True)

    def on_other_changed(self, text: str):
        self.other_column = text or None
        if self._initializing or self.aux_selection != "Other":
            return
        self.redraw(preserve_view=True)

    def redraw(self, preserve_view=False):
        if self.df is None:
            return
        gas = GASES[self.current_gas]
        value_col = gas["value_col"]
        warmup_minutes = self.warmup_spin.value()
        pressure_tol = self.pressure_tol_spin.value()
        cal0_window = (self.cal1_start_spin.value(), self.cal1_end_spin.value())
        cal1_window = (self.cal2_start_spin.value(), self.cal2_end_spin.value())

        df = self.df

        aux_info = aux_trace_info(self.aux_selection, self.current_gas, self.other_column)
        has_aux_panel = aux_info is not None
        aux_key = (self.aux_selection, self.other_column)

        # Capture the current view before tearing down the old Axes, so a
        # masking/averaging control change -- or switching/adding/removing
        # the upper trace -- can redraw without rescaling the main plot.
        # The aux panel's own y-range is only worth preserving if it's
        # still showing the same trace (its scale means something different
        # for a different trace, so let that one re-autoscale) -- for
        # "Other" that also means the same catch-all column, not just the
        # same radio button.
        old_main_view = None
        old_aux_ylim = None
        if preserve_view and self.ax is not None:
            old_main_view = (self.ax.get_xlim(), self.ax.get_ylim())
            if (self.ax_aux is not None and has_aux_panel
                    and self._last_aux_key == aux_key):
                old_aux_ylim = self.ax_aux.get_ylim()

        self.figure.clear()
        if has_aux_panel:
            gs = self.figure.add_gridspec(2, 1, height_ratios=[1, 3])
            ax_aux = self.figure.add_subplot(gs[0])
            ax = self.figure.add_subplot(gs[1], sharex=ax_aux)
        else:
            ax_aux = None
            ax = self.figure.add_subplot(111)
        ax.set_facecolor("#fcfcfb")

        cal = (df["j_sol_cals"].fillna(0).astype(bool)
               | df["j_sol_aircal"].fillna(0).astype(bool))
        shade_intervals(ax, df["datetime"], cal, CAL_SHADE_COLOR, alpha=0.3)

        cal_intervals = merge_close_intervals(
            find_intervals(df["datetime"], cal), pd.Timedelta(seconds=CAL_MERGE_GAP_S)
        )

        bad_pressure = (df["d1_P_mbars"] - D1_P_TARGET_MBARS).abs() > pressure_tol
        bad_pressure = bad_pressure.fillna(False)

        warmup_end = df["datetime"].iloc[0] + pd.Timedelta(minutes=warmup_minutes)
        warmup = df["datetime"] < warmup_end

        shade_intervals(ax, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
        shade_intervals(ax, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)

        # Cal means are estimated from the raw data with these masks applied
        # -- a cal point can be dropped entirely if its window has no valid data.
        exclude_mask = bad_pressure | warmup
        cal_points = cal_mean_points(
            df, cal_intervals, value_col, cal0_window, cal1_window,
            cal_bottles=self.cal_bottles, gas_key=self.current_gas,
            exclude_mask=exclude_mask,
        )

        plot_data = df[["datetime", value_col]].dropna()
        line, = ax.plot(plot_data["datetime"], plot_data[value_col], color=LINE_COLOR, linewidth=1.2)

        cal0_pts = [(t, v) for t, v, state, serial in cal_points if state == 0]
        cal1_pts = [(t, v) for t, v, state, serial in cal_points if state == 1]
        cal0_label = most_common_serial(cal_points, 0) or "Cal 1"
        cal1_label = most_common_serial(cal_points, 1) or "Cal 2"
        self.cal1_box.setTitle(self._cal_box_title(cal0_label, "Cal 1 Mean Window"))
        self.cal2_box.setTitle(self._cal_box_title(cal1_label, "Cal 2 Mean Window"))

        handles = [line]
        labels = [f"{gas['ylabel']} (ambient)"]
        if cal0_pts:
            xs, ys = zip(*cal0_pts)
            handles.append(ax.scatter(xs, ys, color=CAL0_COLOR, s=40, zorder=5, edgecolors="none"))
            labels.append(f"{cal0_label}: {mean_std_label(ys)}")
        if cal1_pts:
            xs, ys = zip(*cal1_pts)
            handles.append(ax.scatter(xs, ys, color=CAL1_COLOR, s=40, zorder=5, edgecolors="none"))
            labels.append(f"{cal1_label}: {mean_std_label(ys)}")

        ax.set_ylabel(gas["ylabel"], color=TEXT_COLOR)
        ax.set_title(gas["title"], color=TEXT_COLOR, loc="left")

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

        if ax_aux is not None:
            aux_col, aux_ylabel = aux_info
            ax_aux.set_facecolor("#fcfcfb")
            shade_intervals(ax_aux, df["datetime"], warmup, WARMUP_EXCLUDE_COLOR, alpha=0.15)
            shade_intervals(ax_aux, df["datetime"], bad_pressure, PRESSURE_EXCLUDE_COLOR, alpha=0.15)

            aux_data = df[["datetime", aux_col]].dropna()
            ax_aux.plot(aux_data["datetime"], aux_data[aux_col], color=LINE_COLOR, linewidth=1.0)

            ax_aux.set_ylabel(aux_ylabel, color=TEXT_COLOR, fontsize=9)
            aux_title = self.other_column if self.aux_selection == "Other" else self.aux_selection
            ax_aux.set_title(aux_title, color=TEXT_COLOR, loc="left", fontsize=10)
            ax_aux.grid(True, color=GRID_COLOR, linewidth=0.6)
            for spine in ax_aux.spines.values():
                spine.set_color(AXIS_COLOR)
            ax_aux.tick_params(colors=MUTED_COLOR, labelsize=8, labelbottom=False)

        # The new Axes were just built and auto-scaled to the full data
        # range -- reset the toolbar's view stack so Home returns to *this*
        # full-scale view, then (optionally) re-apply the pre-redraw zoom on
        # top of it, without pushing that onto the stack.
        self.toolbar.update()
        self.toolbar.push_current()
        if old_main_view is not None:
            ax.set_xlim(old_main_view[0])
            ax.set_ylim(old_main_view[1])
        if old_aux_ylim is not None and ax_aux is not None:
            ax_aux.set_ylim(old_aux_ylim)

        self.ax = ax
        self.ax_aux = ax_aux
        self._had_aux_panel = has_aux_panel
        self._last_aux_key = aux_key

        self.canvas.draw()


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None

    app = QApplication(sys.argv)
    window = UcatsbGui(csv_path)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
