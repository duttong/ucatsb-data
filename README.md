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
- **Trace Above** — optionally add a smaller panel above the main plot: Detector Pressure, T_gas, or "Other". Detector Pressure/T_gas pull from whichever detector the active gas comes from (`d1` for CO2/N2O, `d2` for CH4). "Other" opens a combo box listing every remaining column in the loaded CSV (`oz_o3`, `oz_p`, `oz_t`, `j_sol_cals`, …), so anything in the file can be plotted without a dedicated control. A second combo box overlays any of those columns on a right-hand axis. "No Figure" returns to the single full-size plot.
- **Data Masking** — warm-up exclusion (minutes from the start of the record) and detector pressure tolerance (±mbar around 140 mbar). Both are applied to the raw data *before* cal means are computed, not just drawn as bands — a cal point can disappear entirely if its averaging window has no valid data left. **Flag Air** (0–90 s, default 0 = off) additionally drops the air data immediately following each cal injection, while the detector cells are still clearing cal gas; unlike the other two it affects only the calibrated product, never the cal means or the raw trace. See [Post-cal flush](#post-cal-flush-flag-air).
- **Cal Mean Windows** — one box per cal bottle, titled dynamically from `cals.yaml`: e.g. "50% Cal (CB09960) 206.51 ppm", or just "Cal (CC470901) 402.037 ppm" for a tank with no `info` label. The mole fraction shown is that tank's assigned value for whichever gas is currently active. Each box has a start/end offset in seconds relative to the last point in that calibration period (`Cal_p`), e.g. `-10 s` to `2 s` = `[Cal_p-10s, Cal_p+2s]` (positive values are allowed, reaching past `Cal_p`). Settings are saved per-gas to `ucatsb_gui_config.yaml` and reloaded on next launch.
- **Calibration** — drift model and smoothing window for the two-point calibration, a toggle to overlay the calibrated series on the main plot, and a CSV export. See [Calibration](#calibration) below.

The two views (**Timeseries** and **Calibration**) are tabs sharing this one
control panel — every control affects both.

Zooming/panning (via the matplotlib toolbar) is preserved across masking, averaging, and upper-trace changes — only the Home button resets to full scale.

### Box statistics ("Stats")

The **Stats** toggle at the right-hand end of the Timeseries toolbar turns the
cursor into a marquee: drag a box over a stretch of data and a readout appears
under the toolbar with the sample count, mean, standard deviation and range,
plus a **Copy** button.

```
19:34:06–19:52:06   CO2 (raw)   n=978   mean 393.4 ± 0.5685 ppm
     min 391  max 394.7   (152 in span, outside box vertically)
```

It is a plain selection tool — **no masking is applied**. Drawing the box is
how you say which data you mean; because the box has vertical bounds, drawing
it around the ambient band naturally leaves the cal dives outside it. The
count of points that fell in the time span but outside the box vertically is
always reported, so a box that clipped the data rather than framing it says
so. (That matters for σ: a 2D marquee truncates the distribution and narrows
the spread it reports.)

**A combo box on the readout row picks which plotted trace the statistics
apply to**, listing whatever is currently on the figure:

| Entry | When it appears |
|---|---|
| `<gas> (raw)` | always |
| `<gas> (calibrated)` | "Show calibrated on main plot" is on |
| `<column> (above, left)` | a Trace Above panel is shown |
| `<column> (above, right)` | a right-axis trace is selected |

Boxes can be drawn in either panel, and switching the combo re-reports the
existing box without needing to redraw it. When the chosen trace lives on a
different axis from the box — the upper panel's right axis, or a main-plot
trace when the box was drawn above — the vertical bounds mean nothing for it,
so the selection falls back to the time span alone and the readout says
`time span only`.

The box keeps drag handles after release, so it can be nudged or resized with
the numbers updating, and it survives a masking or averaging change rather
than vanishing on redraw. Note that pressing *inside* an existing box moves or
resizes it — to start a fresh selection, drag somewhere clear of it.

Turning Stats on releases pan or zoom if either is active, since they compete
for the same drag. The Calibration tab has no Stats button: its three panels
each mean something different, so a single readout there would be ambiguous.

### Static figure

```
python3 plot_co2_timeseries.py /path/to/ucatsb-YYYYMMDDHH.csv
```

Produces `<csv_stem>_CO2_ppm.png` next to the input file: a fixed CO2 timeseries with calibration shading, mean cal points, and warm-up/pressure exclusion bands, using the same logic the GUI uses (no interactive controls).

## Calibration

The **Calibration** tab turns the in-flight cal-bottle injections into a
time-varying calibration for the active gas, and shows the diagnostics needed
to decide whether to trust it. The controls in the left panel's *Calibration*
box (drift model, smoothing window) feed it; the masking and cal-window
controls above them feed it too, since they determine which injections exist
and what each one averages to.

### How the calibration is built

Two cal tanks of known mole fraction are injected repeatedly through the
flight. Each injection is averaged over its cal-mean window to give one
*measured response* for that tank at that time — these are the same points
drawn as dots on the main timeseries. Interpolating each tank's response
through time gives, at every ambient timestamp:

```
slope(t)     = (A_hi - A_lo) / (R_hi(t) - R_lo(t))
intercept(t) = A_lo - slope(t) * R_lo(t)
calibrated   = slope(t) * measured(t) + intercept(t)
```

where `A` is a tank's assigned value from `cals.yaml` and `R(t)` its
interpolated measured response.

Two consequences are worth being explicit about:

- **Drift removal and calibration are one step, not two.** The cal responses
  themselves drift, so interpolating between them already tracks the
  instrument's drift. There is no separate detrending pass, and adding one
  would double-count.
- **Both a slope and an intercept are needed** — a single offset is not
  enough. The measured span error runs to several percent (span gain 0.95 on
  the Jul 2026 flight, 1.06 on Feb 2025), so gain has to be corrected too.

**Drift model** (per-gas, saved) controls what the interpolation runs through:

| Model | Nodes used | When |
|---|---|---|
| `linear` (default) | the per-injection means themselves | assumption-free; passes each event's noise straight into the calibrated data |
| `smooth` | centred rolling mean over *N* events (**Smooth over**) | event-to-event scatter is several times the within-event noise and you want it suppressed without flattening real drift |
| `constant` | one flight-mean value | a sanity baseline: what a single static calibration would have given |

### Post-cal flush ("Flag Air")

When the solenoid switches back to ambient at the end of an injection, the
detector cells still hold cal gas, so the air data reads toward the tank for
some seconds afterwards. **Flag Air** (in *Data Masking*, 0–90 s, default 0)
drops that stretch.

It behaves differently from the warm-up and pressure masks, deliberately:

- **The raw trace keeps it.** The recovery is left visible in the raw data —
  it is real instrument behaviour, and seeing it is how the right flush length
  gets chosen.
- **The calibrated series drops it** — along with the cal periods themselves,
  since the calibrated series is the *ambient* record — and the calibrated
  line *breaks* over each gap rather than drawing across it. The flagged span
  is shaded teal on the main plot whether or not "Show calibrated" is on, so
  the rows are visible before the overlay is turned on.
- **The cal means are unaffected.** The flush window covers ambient data only,
  which is disjoint from the cal-mean windows, so no cal point can be lost to
  it and the calibration coefficients are identical with it on or off.
- **The export flags rather than deletes.** Flagged rows are written with
  `is_post_cal_flush=True`, a blank `<col>_cal`, and the raw column untouched.

How long the flush actually takes is instrument- and flight-specific — measure
it rather than assuming. On the two reference flights, the median offset from
the settled value after a cal ends:

| s after cal end | Feb 2025 | Jul 2026 |
|---|---|---|
| 0–4 | −24.7 ppm | −152.6 ppm |
| 5–9 | −4.1 | −65.7 |
| 10–14 | −0.3 | −24.9 |
| 15–19 | +0.0 | −11.8 |
| 20–24 | +0.3 | −3.0 |
| 30–34 | +0.2 | −1.1 |
| 45–49 | +0.3 | −0.2 |

Feb 2025 settles inside ~15 s; Jul 2026 needs ~45 s to come within 0.2 ppm and
is still ~1 ppm low at 30 s. Set the value per flight, not once — the range
runs to 90 s because 30 s is not enough for every flight.

### Reading the three panels

**1 — Cal bottle response, as deviation from assigned value.** Each tank's
per-injection means, plotted as `measured − assigned` so zero is "the tank
reads what it should." Deviation rather than absolute mole fraction is
deliberate: the two tanks sit hundreds of ppm apart, so on a shared absolute
axis the drift you actually care about — a couple of ppm — collapses to a flat
line. The line through the dots is the drift model's nodes; the legend carries
the absolute numbers (`R=` mean measured, `A=` assigned, `d=` difference and
percent). Hollow ✕ marks a **rejected** point: its averaging window straddled a
solenoid transition, so it measured the *other* tank and was filtered out by
serial consistency. A whole series sitting well off zero is the "`cals.yaml`
may not describe this flight" signal.

**2 — Calibration coefficients.** `slope(t)` on the left axis with a dotted
reference at 1.0, `intercept(t)` on the right axis. Hatching marks
**extrapolated** spans, where the calibration is held flat rather than
interpolated. The trustworthy region is the *intersection* of the two tanks'
node spans, not "first to last cal event" — one tank can lose points to masking
while the other keeps going, and in that partial-overlap region one side is
interpolating while the other is flat-held. Long gaps between nodes (> 3× the
median spacing) hatch too. Expect a lot of hatching near the ends of a record;
that is honest, not a fault.

**3 — Residuals / QC.** Two markers per injection, answering two different
questions:

- **Filled dot — closure residual.** Push that injection's own mean back
  through the calibration in force at that instant and subtract the assigned
  value. Under `linear` this is **zero by construction**, because the
  calibration is built to pass through exactly those points. It is a
  self-consistency check, not a quality metric: the dots should sit invisibly
  on the zero line, and a non-zero value there is a bug signal.
- **Hollow circle — leave-one-out.** Drop that injection, rebuild the tank's
  response from its *other* nodes, interpolate back to this injection's time,
  and compare. This answers the question that matters — *"if a cal hadn't
  happened right here, how wrong would the calibration have been?"* — which is
  the error carried by ambient data **between** cal events. It is also
  insensitive to genuine slow drift, since both neighbours share it. The
  per-tank **RMS in the legend is the headline QC number** (≈0.75 ppm CO2 on
  the Jul 2026 flight, at both ends of the range).
- **Shaded band — ±the tank's certified uncertainty** from `cals.yaml`. It is
  invisibly thin next to the leave-one-out scatter, which is the intended
  message: tank certification is nowhere near the dominant error term.

Random scatter about zero is noise, and the RMS summarises it. What to look
for instead is *structure* — a run of same-sign points, a trend, a step —
which means drift the model isn't tracking, or something that changed
mid-flight. Note that the first and last node of each tank have no neighbour on
one side, so leave-one-out predicts them by flat-holding; those two points are
pessimistic relative to the interior ones.

### Header block

Above the panels: gas, mode, span gain, and what fraction of the record is
extrapolated; then per-tank event and rejection counts. Any warnings are
appended and turn the whole block red.

The most common warning is a **mismatch advisory** — a tank whose measured
response sits far from its assigned value, naming the nearest tank in the
roster as a prompt to check `cals.yaml`. This is **advisory only and never
auto-substitutes a tank.** `cals.yaml` describes the *current* run, so applying
it to an older flight can genuinely name the wrong tank (it correctly flags
`CC302489` on the Feb 2025 flight), but the same heuristic misfires where the
offset is real gain error rather than a wrong tank (it names `DT0040700` on Jul
2026, where span gain 0.95 is the real story). Read it as "go check", not as a
diagnosis — which is why span gain is printed beside it.

### Using the result

The calibrated series is the calibrated **ambient** record: cal periods and
the post-cal flush behind them are blanked out of it, so what you get is air
and only air. The raw trace still shows everything, so the cal dives stay
visible underneath — the calibrated line simply breaks over them rather than
drawing across. Nothing is thrown away by this: `cal_slope` and
`cal_intercept` are emitted for every row, so the calibrated value of a cal
period can be recomputed if it is ever wanted.

- **Show calibrated on main plot** overlays this calibration on the timeseries
  tab, keeping the raw trace visible underneath at reduced opacity, with
  extrapolated spans hatched. The cal dots stay raw — they are the
  calibration's inputs. The toggle is session-only and defaults off, so the app
  never starts up showing calibrated data without being asked. Note this is
  *this repo's* calibration, not the CSV's own `*c_ppm`/`*c_ppb` columns.
- **Export calibrated CSV…** writes every row (flagged, never pre-filtered)
  with `<col>_cal`, `cal_slope`, `cal_intercept`, `is_extrapolated`,
  `is_post_cal_flush`, `is_cal_period` and `is_masked`, behind a `#` comment
  block recording the source file, tank serials and assigned values, span
  gain, leave-one-out RMS, control settings, any warnings, and what a blank
  `<col>_cal` means. Read it back with `pd.read_csv(path, comment="#")`.

### When there is no calibration

The tab shows a single centred sentence instead of panels, and the export
button is disabled. This happens for Ozone (its own sensor, not run through the
cal-bottle system), when no cal events survive the current masking (the message
quotes the warm-up and tolerance in force, so it is directly actionable), and
when `cals.yaml` has no assigned value for the active gas. A flight with only
one usable tank degrades to an offset-only correction (`mode: offset`, slope
fixed at 1) rather than refusing outright, and says so in the header.

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
