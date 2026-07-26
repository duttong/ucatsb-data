# ucatsb-data

Post-flight data processing and viewing tools for UCATS-B CSV logs (CO2/N2O
from the Aeris `d1` detector, CH4 from `d2`, plus supporting
pressure/temperature/ozone traces). Companion to the `ucats-b`
data-acquisition repo (`~/code/ucats-b`), which produces the CSV files this
tool reads.

Detector wiring is not fixed across flights — e.g. `d2` used to carry a
redundant CO/N2O channel and now carries CH4/H2O instead, after a CH4 Aeris
was installed. The GUI reads each CSV's actual header and only offers gases
whose column is present in that file, rather than assuming a fixed schema.

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
- **Gas** — switch the main plot between CO2 (`d1_CO2_ppm`), N2O (`d1_N2O_ppb`), and CH4 (`d2_CH4_ppb`), all uncalibrated. Only gases whose column exists in the loaded CSV are offered.
- **Trace Above** — optionally add a smaller panel above the main plot: Detector Pressure, T_gas, oz_o3, oz_p, or oz_t. Detector Pressure/T_gas pull from whichever detector the active gas comes from (`d1` for CO2/N2O, `d2` for CH4). "No Figure" returns to the single full-size plot.
- **Data Masking** — warm-up exclusion (minutes from the start of the record) and detector pressure tolerance (±mbar around 140 mbar). Both are applied to the raw data *before* cal means are computed, not just drawn as bands — a cal point can disappear entirely if its averaging window has no valid data left.
- **Cal Mean Windows** — one box per cal bottle (titled dynamically from `cals.yaml`, e.g. "100% Cal (CC302489) 418.947 ppm" — the mole fraction shown is that tank's assigned value for whichever gas is currently active), each with a start/end offset in seconds relative to the last point in that calibration period (`Cal_p`). Settings are saved per-gas to `ucatsb_gui_config.yaml` and reloaded on next launch.

Zooming/panning (via the matplotlib toolbar) is preserved across masking, averaging, and upper-trace changes — only the Home button resets to full scale.

### Static figure

```
python3 plot_co2_timeseries.py /path/to/ucatsb-YYYYMMDDHH.csv
```

Produces `<csv_stem>_CO2_ppm.png` next to the input file: a fixed CO2 timeseries with calibration shading, mean cal points, and warm-up/pressure exclusion bands, using the same logic the GUI uses (no interactive controls).

## Data assumptions

Expects a UCATS-B CSV with `datetime`, `j_sol_cals`, `j_sol_aircal`, and
whichever of `d1_CO2_ppm`, `d1_N2O_ppb`, `d2_CH4_ppb` (plus matching
`d?_P_mbars`/`d?_T_gas`) and `oz_o3`/`oz_p`/`oz_t` are present — a missing
gas column just removes that gas from the selector rather than failing to
load, but a file with none of the three gas columns will raise an error.

On load, `drop_presync_rows` discards any leading rows recorded before the
datalogger's clock synced (a burst of rows with a too-late timestamp,
followed by a large backward jump to the true time — a known artifact of
this logger, not real data).

`cals.yaml` holds the full tank roster (every cal tank ever used, with its
assigned mole fractions and uncertainties per gas), plus a `cals:` block
naming which two of those tanks (`cal0`/`cal1`) are actually plumbed in for
the current run -- update that block by hand when tanks are swapped between
flights. Only the two assigned tanks are used for matching; calibration
bottle identity (which serial number is flowing at a given `j_sol_cals`
state) is determined by matching the measured concentration against those
two tanks' nominal values, not by trusting any assumed digital-state-to-
serial mapping -- this is self-correcting if the two get swapped without
updating `cals:`. `cals.yaml` here is a local copy of
`~/code/ucats-b/cals.yaml`; resync it by hand if the acquisition repo's
roster changes.

## Files

| File | Purpose |
|---|---|
| `ucatsb_gui.py` | PyQt5 interactive viewer |
| `plot_co2_timeseries.py` | Shared masking/cal-detection logic + standalone static-figure CLI |
| `ucatsb_gui_config.yaml` | Per-gas masking/averaging settings, auto-saved by the GUI |
| `cals.yaml` | Full cal tank roster + which two are assigned for the current run (local copy of `~/code/ucats-b/cals.yaml`) |
