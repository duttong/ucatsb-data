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
  `.ucatsb_gui_state.yaml`, never in a flight's conf, and opening a recent
  file restores that flight's own settings as any other load does. Launching
  with no argument still starts empty — nothing is opened for you.
- **Save…** — writes the current settings to a config file (see [Per-flight settings](#per-flight-settings-configs)). **Nothing is saved automatically**; a `•` on the button means there are unsaved changes.
- **Gas** — switch the main plot between CO2 (`d1_CO2_ppm`), N2O (`d1_N2O_ppb`), CH4 (`d2_CH4_ppb`), O3 (`oz_o3best`) and H2O (`w_H2Obest`), all uncalibrated. Only gases whose column exists in the loaded CSV are offered. O3 and H2O come from their own instruments rather than an Aeris detector, so they have no cal bottles and no masking controls — each has a validity floor instead (see [Validity floors](#validity-floors)).
- **Trace Above** — optionally add a smaller panel above the main plot: Detector Pressure, T_gas, or "Other". Detector Pressure/T_gas pull from whichever detector the active gas comes from (`d1` for CO2/N2O, `d2` for CH4). "Other" opens a combo box listing every remaining column in the loaded CSV (`oz_o3`, `oz_p`, `oz_t`, `j_sol_cals`, …), so anything in the file can be plotted without a dedicated control. A second combo box overlays any of those columns on a right-hand axis. "No Figure" returns to the single full-size plot.
- **Data Masking** — warm-up exclusion (minutes from the start of the record), **End-flight exclude** (minutes from the *end* of the record — the descent is the busiest, least representative part of a flight, so trimming the tail is as routine as trimming the warm-up; default 0, measured back from the last timestamp there is, and shaded and labelled orange together with the warm-up as one exclusion at each end) and detector pressure tolerance (±mbar around 140 mbar). These are applied to the raw data *before* cal means are computed, not just drawn as bands — a cal point can disappear entirely if its averaging window has no valid data left. **Pumps on** (default off) keeps only data recorded with the sample pumps running (`j_pumps = 1`) — air measured with the pumps off is not ambient air. It has to default to off: a lab test or bench calibration runs pumps-off from end to end (the 2026-07-26 file is 100% pumps-off), and enabling it there would leave nothing at all. A row with no `j_pumps` value counts as pumps-off, since an unknown pump state is not evidence the pumps were running; a file whose schema predates the column greys the toggle out. Like warm-up and pressure it feeds the cal means as well as the plot, and its spans are shaded violet. **Flag Air** (0–90 s, default 0 = off) additionally drops the air data immediately following each cal injection, while the detector cells are still clearing cal gas; unlike the other two it affects only the calibrated product, never the cal means or the raw trace. See [Post-cal flush](#post-cal-flush-flag-air). **Copy settings to all gases** applies these three values *and both cal mean windows* to every calibrated gas (CO2/N2O/CH4 — Ozone has no masking at all), since they describe the instrument on this flight rather than the species. Only the drift model and its smoothing window are left alone, being a judgement about how noisy that gas's own cal record is. Settings remain per-gas; the button is a shortcut, not a mode.
- **Cal Mean Windows** — one box per cal bottle, titled dynamically from `cals.yaml`: e.g. "50% Cal (CB09960) 206.51 ppm", or just "Cal (CC470901) 402.037 ppm" for a tank with no `info` label. The mole fraction shown is that tank's assigned value for whichever gas is currently active. Each box has a start/end offset in seconds relative to the last point in that calibration period (`Cal_p`), e.g. `-10 s` to `2 s` = `[Cal_p-10s, Cal_p+2s]` (positive values are allowed, reaching past `Cal_p`). Settings are saved per-gas to the flight's own `<dataset>_conf.yaml` (see [Per-flight settings](#per-flight-settings-dataset_confyaml)) and reloaded whenever that dataset is opened again.
- **Calibration** — drift model and smoothing window for the two-point calibration, and a toggle to overlay the calibrated series on the main plot. See [Calibration](#calibration) below. Writing files out is the [Export](#export) tab's job, not this panel's — both products cover every gas at once, while everything here is per gas.

The **Timeseries**, **Calibration**, **Cal Tanks** and **Export** tabs share
this one control panel — every control affects all of them. The
**Correlations** tab brings its own panel instead (see
[Correlations](#correlations)).

Zooming/panning (via the matplotlib toolbar) is preserved across masking, averaging, and upper-trace changes — only the Home button, a gas change, and a cal-tank change reset to full scale.

### Per-flight settings (configs)

Settings belong to the flight, not to the app: the right warm-up, pressure
tolerance, cal windows and — above all — cal tanks are properties of the
dataset. They are kept in a YAML file beside the CSV, holding a block per gas
plus that flight's cal-tank pairing and any [flagged points](#flagging-errant-points-flag):

```yaml
cals:
  cal0: CC302489
  cal1: CB09960
flagged:                 # points struck out by hand, as raw-file row ranges
  Ozone:
    rows: 16199          # the CSV length they were drawn against
    ranges:
      - [4271, 4272]
      - [6358, 6358]
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

**Nothing is written unless you press Save.** Opening a saved analysis,
changing things, and quitting leaves the file exactly as you found it. In
particular, loading a dataset creates no file at all.

- **Save…** (button beside Load Data, or `File > Save Configuration…`, ⌘S)
  opens a dialog with `<dataset>_conf.yaml` offered as the name. Changing that
  name saves a **second configuration of the same dataset** — as many as you
  like (`..._tight_conf.yaml`, `..._v2.yaml`).
- A `•` on the Save button marks unsaved changes. It is a comparison, not a
  flag: change a spin box and put it back, and the mark clears.
- **Quitting or loading another dataset with unsaved changes** asks
  Save / Don't Save / Cancel. Don't Save leaves the starting state untouched;
  Cancel aborts the quit or load.

**Opening a dataset that has configs:** with just one, it is applied silently.
With several, a chooser lists them — the conventional `<dataset>_conf.yaml`
first and preselected — plus *Start from defaults (open nothing)*. Any
`<dataset stem>*.yaml` beside the CSV counts as that dataset's config, so
whatever you named it is found; `File > Load Configuration…` opens one by name
from anywhere, or switches configuration without reloading the CSV.

A dataset with no config starts from the **shipped defaults**. There is no
app-level template: neither app-level file holds analysis settings at all —
only the [ICARTT header metadata](#header-metadata) and the recent-files list —
and inheriting one flight's tuning into another silently is exactly what
per-dataset configs exist to prevent.

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

### Flagging errant points ("Flag")

Some readings are visibly wrong and no automatic rule catches them — ozone on
the 2026-07-28 flight spikes to 3500 ppb against a 0–900 ppb record, and no
threshold separates those from real data. The **Flag** toolbar toggle, beside
Stats, is the manual override:

- **Left-drag a box** over the bad points to flag them.
- **Right-drag** to unflag.

Flagged rows are struck out with a black `x` at their raw value and removed
from the calibrated (or, for Ozone and H2O, the filtered) record, from the
Correlations tab and from both export products. The raw blue trace still shows
everything, so nothing is hidden — only excluded.

Two asymmetries are worth knowing, both deliberate:

- **Flagging matches the box against the raw (blue) trace**, so a flag means
  the same points after you change the drift model or swap cal tanks. Drawing
  a tight box around the red calibrated overlay can therefore catch nothing —
  it sits an intercept away, about 10 ppm on CO2 — and the readout will say so.
- **Unflagging ignores the box height** and clears the whole time span you drag
  over. A flagged spike is usually off the top of the axes, because the default
  y-range is framed on the filtered data; requiring you to reach its value
  would make it unremovable.

Flagging works on any gas, and the **Flagged Points** panel on the left shows
the count, with **Undo** (this session only) and **Clear** (current gas only).
Tick **Apply to all gases** to spread each new flag across every species at
once, for an inlet or pump problem that ruins them all together.

Flags on a cal-bottle gas behave like every other mask: a flagged row inside a
cal window is dropped before that injection's mean is computed, so flagging a
visibly bad injection genuinely changes the calibration.

The view is held still while you flag, since you are usually zoomed in on the
points you are removing — press **Home** afterwards to reframe the y-axis on
what is left.

**Flags are saved with the flight, but only when you press Save**, like every
other setting. They live in the flight's config as row ranges, so a few
thousand flagged points cost a handful of lines.

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
  `<gas>_is_post_cal_flush=1`, a blank `<gas>_cal`, and the raw column
  untouched.

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

Writing this calibration out to a file is the [Export](#export) tab's job.

### When there is no calibration

The tab shows a single centred sentence instead of panels, and the Export tab
reports the gas as `NOT calibrated` with the same reason. This happens for
Ozone (its own sensor, not run through the
cal-bottle system), when no cal events survive the current masking (the message
quotes the warm-up and tolerance in force, so it is directly actionable), and
when `cals.yaml` has no assigned value for the active gas. A flight with only
one usable tank degrades to an offset-only correction (`mode: offset`, slope
fixed at 1) rather than refusing outright, and says so in the header.

## Export

The **Export** tab writes two files, both covering **all five gases** (CO2,
N2O, CH4, O3, H2O — whichever the loaded CSV actually has). It is its own tab
rather than a button in the left panel because both products are whole-flight
deliverables, while every control in that panel is per gas.

At the top, **What will be exported** lists each gas as it stands under the
current settings — calibrated, filtered-only, or `NOT calibrated` with the
reason — along with the row counts, so you can see what you are about to write
before you write it.

### Derived CSV — the companion to the raw file

Meant to sit **beside** the acquisition CSV, not replace it. It has **one row
for every row of the raw file, in the same order**, so the two can be opened
side by side, or pasted together, in Excel or Igor Pro with no alignment step
at all. The pre-sync rows the analysis drops (see [Data
assumptions](#data-assumptions)) are still present, blank in the derived
columns and flagged in a `presync` column, precisely so that promise holds — a
quiet offset of a few dozen rows between two files that *look* alignable is an
effective way to corrupt an analysis.

Columns, per gas:

| Column | What it is |
|---|---|
| `<gas>_cal` | the calibrated good-ambient record (blank where it is not good ambient) |
| `<gas>_cal_unc` | 1σ on that value |
| `<gas>_cal_slope`, `<gas>_cal_intercept` | the coefficients, on **every** row, so a blank row can be recomputed |
| `<gas>_is_cal_period`, `<gas>_is_post_cal_flush`, `<gas>_is_masked`, `<gas>_is_extrapolated` | why a row is blank |
| `<gas>_is_flagged` | struck out by hand (see [Flagging errant points](#flagging-errant-points-flag)). A subset of `is_masked` for a calibrated gas, written separately because it is the one removal that cannot be reconstructed from the settings |
| `<gas>_filtered`, `<gas>_below_floor` | for O3 and H2O, which have no cal bottles — the raw value with below-floor faults and flagged points removed |

plus `datetime`, `time_s` (seconds from midnight UTC, which Igor and Excel
both plot far more readily than a parsed date), and the raw value columns
under their original names.

- **Masks are written as 1/0**, not `True`/`False`, which Igor loads as a text
  wave and Excel as text. A blank mask cell means "not evaluated" (a pre-sync
  row), not "false".
- **The raw echo columns are exact copies of the source columns**, pre-sync
  rows included — those rows have an unreliable *clock*, not an unreliable
  *reading*. Uncheck *Include the raw value columns* to keep the file strictly
  complementary.
- **The provenance notes go to a sidecar `<name>_notes.txt` by default**,
  because neither Excel nor Igor skips a leading `#` comment block without
  being told to. Tick *Put the provenance notes in the CSV as # lines* for a
  self-contained file you read back with `pd.read_csv(path, comment="#")`.
- **A gas with no calibration still exports** its raw column and its masks —
  which rows were air is a useful answer whether or not there is a
  calibration — and the notes say why it is uncalibrated.

### ICARTT (.ict)

The archive format (file format index 1001), for delivering the flight. Good
ambient data only: every row the analysis blanked is written as the format's
missing value, `-99999`, which is exactly what that flag means — so no mask
columns are written or needed. Time is seconds from midnight UTC on the data
start date, counting past 86400 for a flight that crosses midnight. A gas with
no calibration is left out of the file entirely, rather than delivered as
uncalibrated counts under a calibrated-looking name.

Two things it will tell you about after writing, and both are worth reading:

- **Rows with repeated or backward timestamps are dropped** (1,435 of them in
  the Feb 2025 file). ICARTT's independent variable must increase strictly, so
  a duplicate makes the file *invalid* rather than merely untidy. Better to
  see the count here than to have an archive's validator find it.
- **Rows with no value for any gas are dropped** by default — a row that is
  `-99999` in every column carries nothing but a timestamp. Untick *Drop rows
  with no value for any gas* to keep the time base unbroken.

The file follows the conventions of the **sister UCATS instrument's delivered
files** rather than generic ICARTT practice, so a campaign's UCATS and UCATS-B
files can be processed the same way: missing data is `-99999`, units carry the
trailing `v` (`ppmv`/`ppbv`), variable lines are
`name, unit, standard_name, description`, and the 1σ variable for a species is
`<species>e_<suffix>` — `CO2e_UCATSB` beside `CO2_UCATSB`.

The `standard_name` field is not free text: it comes from the NASA ESDS
Atmospheric Composition Variable Standard Names Convention that ICARTT V2.0
requires, built as
`MeasurementCategory_CoreName_AcquisitionMethod_DescriptiveAttributes`. For a
trace gas the attributes are the measurement specificity (`S` for a single
species) and how the value is reported (`DMF` = molar fraction with respect to
dry air), giving `Gas_CO2_InSitu_S_DMF` for CO2, N2O and CH4.

Two gases depart from that pattern. **H2O is `Met_H2OMF_InSitu_None`**: water
vapour is not a `Gas` name in the convention at all, but a `Met` one, and
`H2OMF` is a mole fraction with respect to ambient air. **Ozone is
`Gas_O3_InSitu_S_AVMR`**: a volumetric mixing ratio with respect to ambient
air, since nothing dries the 2B monitor's sample stream. Both live in
`GASES[gas]["standard_name"]`.

#### Header metadata

The **Header metadata** form fills in everything the header carries that
cannot be derived from the data, in the order the file itself uses. Two fields
are computed for you rather than typed:

- `UNCERTAINTY` gets the **median 1σ actually computed for each gas**, with
  anything you type appended after it. A hand-typed uncertainty goes stale the
  moment a drift model or cal window changes.
- `ULOD_FLAG`/`LLOD_FLAG` are constants of the format (-7777 / -8888); only
  the LOD *values* are yours to set, and default to `N/A`.

Two fields are written into the file **verbatim**, line breaks and all,
because they are whole sections rather than single keyword values:

- **Special comments** — free text, and in the delivered UCATS files this is
  where the error estimates are explained and users are asked to contact the
  PIs. What you type is the whole section; nothing is appended to it.
- **Revision history** — one `R#: description` per line, and it
  **accumulates**: an R0 file still lists its RA line above the R0 one. Add a
  line here each time you bump `Revision`.

`Data ID`, `Location ID` and `Revision` also build the suggested file name,
`dataID_locationID_YYYYMMDD_R#.ict` (hyphens are kept, so `SABRE-UCATSB`
works; underscores are stripped, since `_` separates the name's own fields). Leaving the PI, mission, platform or
location ID blank is allowed but prompts first, since an archive will normally
reject a file like that.

**The metadata is saved in `ucatsb_gui_config.yaml`, not in the flight's own
config** — the PI, affiliation, project and stipulations are properties of the
campaign, so storing them per flight would mean retyping them for every file.
That file is tracked in the repository, so a fresh clone arrives with the
campaign already filled in rather than needing it retyped per machine.
**Save defaults** is what writes them; a `•` on the button means there are
unsaved changes, and quitting with changes prompts. Exporting always uses
what is in the boxes, saved or not, so a one-off mission name needs no save.

## Correlations

The **Correlations** tab plots one calibrated tracer against another —
`N2O` vs `CO2` and so on. This tab has **its own left panel**: everything in
the usual panel is per-gas, and this figure is about two gases at once, so
the per-gas controls are replaced rather than left there to be misread.
Masking, cal windows and drift model for each gas still live on the other
tabs, and each gas is analysed with *its own* saved settings.

- **X axis / Y axis** — any two gases in the file, Ozone included. **Swap
  axes** flips them.
- **Marker size** — marker diameter in points. A whole flight is tens of
  thousands of points; small markers show structure that large ones fill in.
- **Color by** — a z-axis encoding: **Time** or **Pressure (`oz_p`)**, drawn
  with the `turbo` rainbow and always accompanied by a labelled colorbar (a
  continuous color scale with no key is unreadable). `turbo` rather than
  `jet`: same rainbow ordering, but rebuilt with monotonic luminance, so it
  does not invent bands of false structure. Points with no value for the
  chosen variable are dropped rather than drawn invisibly, and the `n` note
  says so. Only encodings the loaded CSV can supply are offered — a flight
  whose schema lacks `oz_p` gets Time alone.
- **Error bars (1σ)** — see below.
- **Linear fit (OLS)** — off by default. These plots usually have real
  structure (branches, mixing lines, profiles) that a single straight line
  describes poorly, and its slope then says more about how the flight was
  sampled than about the tracers. When on, the fit is drawn in red and its
  slope ± standard error, intercept and *r* appear both on the figure and in
  the panel; when off, the panel still reports `n` and each tracer's
  mean ± standard deviation.

**Each axis is calibrated wherever a calibration exists**, deliberately: a
tracer-tracer slope taken from uncalibrated counts carries that detector's
gain error straight into the slope, which is the number the plot exists to
produce. Those series come from the same `calibrate_series` output the
Timeseries tab overlays, so they are already good-air-only.

**Ozone is the exception it has to be.** It has no cal bottles, so there is
nothing to calibrate it against; it goes on the axis as the ozone
instrument's own product, `oz_o3best`, with its below-floor readings removed
(see [Validity floors](#validity-floors)). The axis label and a line
in the figure's note block both say so, and that axis gets no error bars — there is
no calibration to propagate, and a zero-width bar would claim a precision
nobody established. It needs no time alignment: ozone shares the CSV's
timestamps, so the rows line up directly. It does report on a slower cadence
(~2 s against the Aeris' 1 Hz on the Feb 2025 flight, so only about half the
rows carry a value), which shows up as a lower `n`. Ozone carries no masking
of its own either, but its partner axis is blanked wherever *that* gas was in
cal, flushing or masked, and the pairing is an intersection — so the partner's
masking reaches the ozone axis too.

A point is plotted when **both** tracers have a usable value there. The figure says how many rows that left, what fraction sits where a
calibration is extrapolated, and — because it matters here more than anywhere
— when **Flag Air is 0**: post-cal flush points run in a line from the tank's
composition toward the atmosphere's, which looks exactly like a tracer-tracer
correlation and drags the fit. On the Feb 2025 file, raising Flag Air to 45 s
drops that trail (n 11519 → 10088) and moves the slope 0.906 → 0.883.

An ordinary least-squares fit is available via **Linear fit (OLS)** (off by
default), drawn in red with slope ± standard error, intercept and *r* both on
the figure and in the panel on the left — which survives zooming and can be
read while the plot is somewhere else.

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

### Validity floors

The two instruments with no cal bottles — ozone and water vapour — each
declare a floor below which a reading is a fault rather than a measurement:

| Gas | Column | Floor | Effect on the reference flights |
|---|---|---|---|
| O3 | `oz_o3best` | **−15 ppb** (`O3_VALID_MIN_PPB`) | removes 7 readings on Feb 2025 (−22 to −2292 ppb) |
| H2O | `w_H2Obest` | **−5 ppm** (`H2O_VALID_MIN_PPM`) | removes nothing on Jul 2026 (its minimum is +14 ppm) — precautionary |

The H2O floor is there for the fault mode when it appears; it changes no
current data.

A floor sits well below zero on purpose. A real near-zero ozone measurement
scatters negative, and that flight has 168 readings between −15 and 0 ppb that
are the sensor's noise about a small true value. Clipping at zero would bias
the low end upward and hide how noisy the instrument actually is.

The timeseries of a gas with a floor therefore shows **two traces, in the same colors the Aeris
gases use**: raw in blue underneath (faded), filtered in red on top — the same
"blue is everything recorded, red is the series you should read" convention as
the calibrated overlay. The red line *breaks* over each removed reading rather
than drawing across it, and the note block says how many were removed.

The default y-range is framed on the **filtered** data. One −2292 ppb fault
would otherwise set the scale and squash the entire real ozone record into the
top fifth of the axes; the raw trace is still drawn and simply runs off-scale,
where zooming out reaches it.

These are per-gas properties (`valid_min` in `GASES`), not maskable
settings — it is a physical floor, not a judgement call, so it is not
adjustable from the control panel. The cal-bottle gases declare none, and a
gas without a floor is drawn as a single trace at full opacity.

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
| `ucatsb_analysis.py` | Shared masking, cal-detection, calibration and export logic. A library, never run directly, and deliberately free of any GUI dependency so a script or notebook can reuse it |
| `ucatsb_gui_config.yaml` | App-level and **shared through the repository**: the [ICARTT header metadata](#header-metadata), which describes the campaign rather than one machine. No analysis settings — those live per dataset (see [Per-flight settings](#per-flight-settings-configs)). Written only by **Save defaults** |
| `.ucatsb_gui_state.yaml` | App-level and **local, gitignored, hidden**: the recent-files list, which is absolute paths from whichever machine ran it. Rewritten on every dataset load — which is why it is not in the file above — and recreated when absent |
| `<dataset>_conf.yaml` | Written beside each loaded CSV: that flight's own per-gas settings **and** its cal-tank pairing |
| `cals.yaml` | Full cal tank roster + which two are assigned for the current run (local copy of `~/code/ucats-b/cals.yaml`) |
