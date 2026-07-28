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
- **Load Data** — opens a menu: **Open File…** for the file browser, then the
  last 10 datasets opened, then **Clear Recent Files**. The same commands are
  in the menu bar as **File > Open…** (⌘O) and **File > Open Recent**. The
  currently-loaded file is check-marked; one whose file has since moved or
  whose volume is unmounted is shown greyed and marked `(missing)` rather than
  forgotten, and is dropped only if opening it actually fails. Two datasets
  with the same name in different directories are disambiguated by their
  parent directory (`ucatsb-20250218All.csv — 20250218`). The list lives in
  `ucatsb_gui_config.yaml`, never in a flight's conf, and opening a recent
  file restores that flight's own settings as any other load does. Launching
  with no argument still starts empty — nothing is opened for you.
- **Gas** — switch the main plot between CO2 (`d1_CO2_ppm`), N2O (`d1_N2O_ppb`), and CH4 (`d2_CH4_ppb`), all uncalibrated. Only gases whose column exists in the loaded CSV are offered.
- **Trace Above** — optionally add a smaller panel above the main plot: Detector Pressure, T_gas, or "Other". Detector Pressure/T_gas pull from whichever detector the active gas comes from (`d1` for CO2/N2O, `d2` for CH4). "Other" opens a combo box listing every remaining column in the loaded CSV (`oz_o3`, `oz_p`, `oz_t`, `j_sol_cals`, …), so anything in the file can be plotted without a dedicated control. A second combo box overlays any of those columns on a right-hand axis. "No Figure" returns to the single full-size plot.
- **Data Masking** — warm-up exclusion (minutes from the start of the record) and detector pressure tolerance (±mbar around 140 mbar). Both are applied to the raw data *before* cal means are computed, not just drawn as bands — a cal point can disappear entirely if its averaging window has no valid data left. **Flag Air** (0–90 s, default 0 = off) additionally drops the air data immediately following each cal injection, while the detector cells are still clearing cal gas; unlike the other two it affects only the calibrated product, never the cal means or the raw trace. See [Post-cal flush](#post-cal-flush-flag-air). **Copy settings to all gases** applies these three values *and both cal mean windows* to every calibrated gas (CO2/N2O/CH4 — Ozone has no masking at all), since they describe the instrument on this flight rather than the species. Only the drift model and its smoothing window are left alone, being a judgement about how noisy that gas's own cal record is. Settings remain per-gas; the button is a shortcut, not a mode.
- **Cal Mean Windows** — one box per cal bottle, titled dynamically from `cals.yaml`: e.g. "50% Cal (CB09960) 206.51 ppm", or just "Cal (CC470901) 402.037 ppm" for a tank with no `info` label. The mole fraction shown is that tank's assigned value for whichever gas is currently active. Each box has a start/end offset in seconds relative to the last point in that calibration period (`Cal_p`), e.g. `-10 s` to `2 s` = `[Cal_p-10s, Cal_p+2s]` (positive values are allowed, reaching past `Cal_p`). Settings are saved per-gas to the flight's own `<dataset>_conf.yaml` (see [Per-flight settings](#per-flight-settings-dataset_confyaml)) and reloaded whenever that dataset is opened again.
- **Calibration** — drift model and smoothing window for the two-point calibration, a toggle to overlay the calibrated series on the main plot, and a CSV export. See [Calibration](#calibration) below.

The **Timeseries**, **Calibration** and **Cal Tanks** tabs share this one
control panel — every control affects all three. The **Correlations** tab
brings its own panel instead (see [Correlations](#correlations)).

Zooming/panning (via the matplotlib toolbar) is preserved across masking, averaging, and upper-trace changes — only the Home button, a gas change, and a cal-tank change reset to full scale.

### Per-flight settings (`<dataset>_conf.yaml`)

Settings belong to the flight, not to the app: the right warm-up, pressure
tolerance, cal windows and — above all — cal tanks are properties of the
dataset. So when a CSV is loaded, the viewer writes a settings file **next to
that CSV**, named after it:

```
/path/to/ucatsb-2025021816.csv
/path/to/ucatsb-2025021816_conf.yaml    <- created on load
```

It holds a block for every gas plus the flight's cal-tank pairing, and is
rewritten on every control change:

```yaml
cals:
  cal0: CC302489
  cal1: CB09960
CO2:
  warmup_min: 12
  pressure_tol_mbar: 0.3
  flag_air_s: 38
  cal1_window_s: [-8, -1]
  cal2_window_s: [-8, -1]
  drift_model: linear
  drift_smooth_events: 3
N2O:
  ...
```

Reopening that CSV — in any session, from any directory — restores exactly
those choices. `ucatsb_gui.py`'s own `ucatsb_gui_config.yaml` keeps working as
before, but its job is now to be the **template** a flight opened for the
first time starts from, so the settings you have converged on carry over to
the next flight instead of reverting to shipped defaults. Two things follow:

- Editing a control writes both files (the flight's, and the template).
- The **tank pairing is never templated** — a flight with no conf file always
  starts from `cals.yaml`, never from the last flight you had open. Tanks get
  swapped between flights, and inheriting the previous flight's pair silently
  is the one error that would corrupt every calibrated number without
  announcing itself.

If the dataset's directory is not writable (a read-only archive or a mounted
share), the viewer says so on stderr and falls back to the app-level config
rather than refusing to open the file.

### Cal Tanks tab

Which two tanks were plumbed in for **this** flight. `cals.yaml` names a
`cal0`/`cal1` pair, but it describes the tanks in use *now* — an older flight
almost certainly flew different ones, and using the wrong assigned values
biases every calibrated number for every gas. The tab lets you pick from the
full roster and saves the choice in the flight's `_conf.yaml`.

Both combo boxes list every tank in `cals.yaml`, labelled with its `info` tag
where it has one (`CB09960 (50%)`). Below them, the assigned value and
uncertainty for each gas are shown side by side for the chosen pair, so the
pick can be checked against what the flight actually measured (the
Calibration tab's header reports `measured R` vs `assigned A` per tank).
Warnings appear for a pairing that cannot work: a serial missing from the
roster, the same tank picked twice (no span — the calibration degrades to
offset-only), or a tank with no assigned value for one of the gases in the
file.

`cal0`/`cal1` mirror the key names in `cals.yaml`; **the order does not
matter.** Which tank is flowing in a given cal window is identified from the
measured concentration, not from this ordering (see [Data
assumptions](#data-assumptions)). What matters is that the *pair* is right.

**Reset to cals.yaml default** puts the pairing back to the `cals:` block.

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
| `smooth` | centred rolling mean over *N* events (the `N ev` spin box beside the model) | event-to-event scatter is several times the within-event noise and you want it suppressed without flattening real drift |
| `constant` | one flight-mean value | a sanity baseline: what a single static calibration would have given |

### The switch-over sample (always flagged)

`j_sol_cals`/`j_sol_aircal` describe the **valve**, but the reading carrying a
given timestamp is of gas that entered the cell a second or so earlier. So the
first ambient-flagged row after each injection is not air at all — it is the
tank, at full concentration:

```
row   datetime             CO2      N2O    j_sol_cals  j_sol_aircal  d1_P_mbars
3192  2025-02-18 16:48:59  206.64  163.79     0.0          1.0         140.64   <- cal
3193  2025-02-18 16:49:00  206.60  163.81     0.0          0.0         140.61   <- "air"
3194  2025-02-18 16:49:01  197.48  153.63     0.0          0.0         139.49
3195  2025-02-18 16:49:02  172.76  125.53     0.0          0.0         122.94   <- pressure mask
```

The cause is a **real ~1 s latency in the serial data** — the reading arrives
after the gas it describes. Measured across all 45 solenoid transitions in the
Feb 2025 file (1 Hz sampling), the response starts moving a median of 1 sample
after the flag changes, and it is symmetric, which is what identifies it as a
fixed transport delay rather than an edge artifact:

| Edge | Samples until the response moves | Seconds |
|---|---|---|
| air → cal (rising) | median 1, range 0–2 | median 1.0, range 0–3 |
| cal → air (falling) | median 1, range 0–2 | median 1.0, range 0–2 |

The rising edge needs no handling: the first cal-flagged rows still read air,
but cal periods run 40–97 s and the mean windows are anchored to `Cal_p` at
the *end*, so even a `-30 s` window never reaches within 2 s of an interval's
start (0 of 41 events on that flight). The falling edge is the one that
escapes into the ambient record.

Nothing else caught it. The cal mask goes by the flag, which is already off;
the pressure mask misses it because the pressure is still in spec at that
sample (the transient hits a second later); and the flush window starts
*after* the interval's end, which is that very row's timestamp. One row per
cal event — 18 on the Feb 2025 flight — sitting in the middle of the ambient
record at the tank's composition. On a timeseries these are single samples
hiding at the trailing edge of each cal dive; on the Correlations tab they
pile up into an obvious tank-composition cluster, which is how they were
found.

These rows are now **always** flagged as not-air, independently of Flag Air,
because they are mislabelled rather than merely early: the gray cal band
extends over them and they are blanked from the calibrated series and the
export (`is_cal_period=True`). This does *not* touch the cal periods
themselves, so `Cal_p`, every cal mean and every calibration coefficient are
unchanged.

### Post-cal flush ("Flag Air")

The sample above is the leading edge of a longer effect. When the solenoid
switches back to ambient at the end of an injection, the
detector cells still hold cal gas, so the air data reads toward the tank for
some seconds afterwards. **Flag Air** (in *Data Masking*, 0–90 s, default 0)
drops that stretch. With Flag Air at 0 the flush is kept, so tank-ward points
remain in the calibrated record for a few seconds after each cal — 21 of them
below 300 ppm on the Feb 2025 CO2 record, versus none at 45 s.

That 21 includes the second sample of a 2 s lag, which the always-on
switch-over mask (one timestamp) does not reach:

```
16:25:44 = 206.6      flag: cal
16:25:45 = 206.6  *   <- switch-over mask
16:25:46 = 206.1  *   <- still tank; needs Flag Air
16:25:47 = 271.4  *
```

So the two mechanisms split the job deliberately: the switch-over mask is the
zero-configuration floor that keeps a full-strength tank reading out of the
air record no matter what, and Flag Air is the tunable part that covers the
rest of the lag *and* the cell-clearing tail behind it.

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

## Correlations

The **Correlations** tab plots one calibrated tracer against another —
`N2O` vs `CO2` and so on. This tab has **its own left panel**: everything in
the usual panel is per-gas, and this figure is about two gases at once, so
the per-gas controls are replaced rather than left there to be misread.
Masking, cal windows and drift model for each gas still live on the other
tabs, and each gas is analysed with *its own* saved settings.

- **X axis / Y axis** — any two calibratable gases in the file (Ozone is not
  offered: no cal bottles, so no calibrated series). **Swap axes** flips them.
- **Marker size** — marker diameter in points. A whole flight is tens of
  thousands of points; small markers show structure that large ones fill in.
- **Error bars (1σ)** — see below.

**It is calibrated data only**, deliberately: a tracer-tracer slope taken from
uncalibrated counts carries each detector's gain error straight into the
slope, which is the number the plot exists to produce. Both series come from
the same `calibrate_series` output the Timeseries tab overlays, so they are
already good-air-only, and a point is plotted when **both** tracers are good
there. The figure says how many rows that left, what fraction sits where a
calibration is extrapolated, and — because it matters here more than anywhere
— when **Flag Air is 0**: post-cal flush points run in a line from the tank's
composition toward the atmosphere's, which looks exactly like a tracer-tracer
correlation and drags the fit. On the Feb 2025 file, raising Flag Air to 45 s
drops that trail (n 11519 → 10088) and moves the slope 0.906 → 0.883.

An ordinary least-squares fit is drawn in red, with slope ± standard error,
intercept and *r* both on the figure and in the panel on the left (which
survives zooming and can be read while the plot is somewhere else).

### Where the error bars come from

Yes — a meaningful uncertainty does follow from the calibration, and that is
exactly what the toggle shows. Writing the two-point calibration as a blend
of the two assigned values, `c = (1-f)·A_lo + f·A_hi` where
`f = (c - A_lo)/(A_hi - A_lo)`:

```
var(c) = ((1-f)·sA_lo)² + (f·sA_hi)² + (slope·(1-f)·sR_lo)² + (slope·f·sR_hi)²
```

with two inputs per tank:

| Term | What it is |
|---|---|
| `sA` | the tank's assigned-value uncertainty, from `cals.yaml`'s `<GAS>_unc`. Missing is treated as 0 — an underestimate, but inventing a number would be worse |
| `sR` | how well the drift model reproduces that tank's measured response: `sqrt(loo² + closure²)` |

`sR` is where **the choice of calibration method shows up in the numbers**.
Leave-one-out RMS is the honest scatter under a model that interpolates
through every node (`linear`, whose closure residual is 0 by construction);
the closure RMS is the error the model itself introduces by *not* passing
through the nodes, which is 0 under `linear` and dominant under `constant`.
On the Feb 2025 flight, switching CO2 from `linear` to `constant` moves the
median 1σ from 2.67 to 3.34 ppm.

Two things it deliberately does **not** include: the instrument's own
single-sample noise (nothing in the calibration constrains it — it would have
to come from the raw trace's high-frequency scatter), and any inflation over
the extrapolated spans, which are reported as a count instead of being buried
in a number that would then look like ordinary uncertainty. Note also that
these uncertainties are mostly **systematic** — they shift a whole flight
together rather than scattering point to point — which is why the fit is plain
OLS rather than weighted by them.

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
| `ucatsb_gui_config.yaml` | Per-gas masking/averaging settings, auto-saved by the GUI; the template a newly-opened flight starts from, plus the recent-files list |
| `<dataset>_conf.yaml` | Written beside each loaded CSV: that flight's own per-gas settings **and** its cal-tank pairing |
| `cals.yaml` | Full cal tank roster + which two are assigned for the current run (local copy of `~/code/ucats-b/cals.yaml`) |
