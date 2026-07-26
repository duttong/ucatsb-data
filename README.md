# ucatsb-data

Post-flight data processing and viewing tools for UCATS-B CSV logs (CO2/N2O
from the Aeris `d1` detector, plus supporting pressure/temperature/ozone
traces). Companion to the `ucats-b` data-acquisition repo (`~/code/ucats-b`),
which produces the CSV files this tool reads.

## Requirements

- Python 3.9+
- `pandas>=2.2`
- `matplotlib>=3.9`
- `PyQt5>=5.15`
- `PyYAML>=6.0`

Install with:

```
pip install -r requirements.txt
```

## Usage

### Interactive viewer

```
python3 ucatsb_gui.py /path/to/ucatsb-YYYYMMDDHH.csv
```

Left panel:
- **Gas** — switch the main plot between CO2 (`d1_CO2_ppm`) and N2O (`d1_N2O_ppb`), both uncalibrated.
- **Trace Above** — optionally add a smaller panel above the main plot: Detector Pressure, T_gas, oz_o3, oz_p, or oz_t. Detector Pressure/T_gas pull from `d1_*` when CO2 is selected and `d2_*` when N2O is selected, since N2O has a second channel on the CO/N2O (`d2`) detector. "No Figure" returns to the single full-size plot.
- **Data Masking** — warm-up exclusion (minutes from the start of the record) and detector pressure tolerance (±mbar around 140 mbar). Both are applied to the raw data *before* cal means are computed, not just drawn as bands — a cal point can disappear entirely if its averaging window has no valid data left.
- **Cal Mean Windows** — one box per cal bottle (titled dynamically from `cals.yaml`, e.g. "100% Cal (CC302489)"), each with a start/end offset in seconds relative to the last point in that calibration period (`Cal_p`). Settings are saved per-gas to `ucatsb_gui_config.yaml` and reloaded on next launch.

Zooming/panning (via the matplotlib toolbar) is preserved across masking, averaging, and upper-trace changes — only the Home button resets to full scale.

### Static figure

```
python3 plot_co2_timeseries.py /path/to/ucatsb-YYYYMMDDHH.csv
```

Produces `<csv_stem>_CO2_ppm.png` next to the input file: a fixed CO2 timeseries with calibration shading, mean cal points, and warm-up/pressure exclusion bands, using the same logic the GUI uses (no interactive controls).

## Data assumptions

Expects a UCATS-B CSV with (at least) these columns: `datetime`, `d1_CO2_ppm`,
`d1_N2O_ppb`, `d1_P_mbars`, `d2_P_mbars`, `d1_T_gas`, `d2_T_gas`, `oz_o3`,
`oz_p`, `oz_t`, `j_sol_cals`, `j_sol_aircal`.

On load, `drop_presync_rows` discards any leading rows recorded before the
datalogger's clock synced (a burst of rows with a too-late timestamp,
followed by a large backward jump to the true time — a known artifact of
this logger, not real data).

Calibration bottle identity (which serial number is flowing at a given
`j_sol_cals` state) is determined by matching the measured concentration
against the nominal values in `cals.yaml`, not by trusting any assumed
digital-state-to-serial mapping — this is self-correcting if bottles are
swapped between flights. `cals.yaml` here is a local copy of
`~/code/ucats-b/cals.yaml`; update it by hand if the acquisition repo's
bottles/serials change.

## Files

| File | Purpose |
|---|---|
| `ucatsb_gui.py` | PyQt5 interactive viewer |
| `plot_co2_timeseries.py` | Shared masking/cal-detection logic + standalone static-figure CLI |
| `ucatsb_gui_config.yaml` | Per-gas masking/averaging settings, auto-saved by the GUI |
| `cals.yaml` | Cal bottle serials and nominal concentrations (local copy of `~/code/ucats-b/cals.yaml`) |
