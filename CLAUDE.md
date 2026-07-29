# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

Post-flight analysis tools for UCATS-B CSV logs. This is the analysis-side
companion to `~/code/ucats-b`, which runs on the aircraft and produces
`ucatsb-YYYYMMDDHH.csv` files (documented in that repo's `config-plot.yaml`
comment block). This repo does not acquire data — it only reads CSVs already
on disk, from any flight, passed as a CLI argument or picked via the GUI's
"Load Data" button (never hardcode a specific flight's filename or date).

## Requirements

`pandas>=2.2`, `matplotlib>=3.9`, `PyQt5>=5.15`, `PyYAML>=6.0`, Python 3.9+.
See `requirements.txt`.

## Running

- `python3 ucatsb_gui.py [csv_file]` — interactive PyQt5 viewer. The CSV
  argument is optional; a file can also be (re)picked at any time from
  within the GUI via the "Load Data" button, so the app can start with no
  argument at all.
- `python3 plot_co2_timeseries.py <csv_file>` — static PNG, no GUI

No test suite or CI. Both scripts are read-verified by importing the module
under `QT_QPA_PLATFORM=offscreen` and grabbing a screenshot rather than by
unit tests — there isn't a headless-safe way to assert on plot pixels, so
manual visual review (via `Read` on the saved PNG) is how changes here get
checked.

## Architecture

`plot_co2_timeseries.py` is the shared logic module (despite the CO2-specific
name/CLI, its functions are gas-agnostic) and also a standalone CLI producing
one fixed CO2 figure. `ucatsb_gui.py` imports from it rather than
duplicating: `drop_presync_rows`, `find_intervals`, `merge_close_intervals`,
`shade_intervals`, `cal_mean_points`, `load_cal_roster`,
`load_cal_assignment`, `select_cal_bottles`, `most_common_serial`,
`mean_std_label`, `CALS_YAML_PATH`, plus the calibration functions below. Any
change to masking/cal-detection/calibration behavior belongs in that shared
module so both stay in sync.

### Data pipeline (per redraw)

1. `drop_presync_rows` — trims leading rows recorded before the datalogger's
   clock synced. On startup the logger's system clock can be stale for the
   first ~60s until GPS/NTP sync, producing a burst of too-late timestamps
   followed by one large backward jump to the true time; everything before
   that jump is dropped. This is a real recurring artifact of this specific
   logger, not a one-off from a single flight's file.
2. `j_sol_cals` / `j_sol_aircal` flag calibration/cal-air periods.
   `find_intervals` + `merge_close_intervals` turn that boolean mask into
   discrete cal-event windows, bridging events split by a single dropped
   sample (gap ≤ `CAL_MERGE_GAP_S`).
3. Warm-up (first N minutes) and out-of-spec detector pressure
   (`|d1_P_mbars - 140| > tol`) masks are computed and **applied to the raw
   data before cal means are estimated** (`exclude_mask` param on
   `cal_mean_points`) — not just drawn as bands. A cal point can disappear
   entirely if its averaging window has no unmasked data left. This was a
   deliberate correction: the masks used to be purely visual annotations
   that didn't affect the computed means, which was silently wrong. The same
   mask is separately handed to `calibrate_series`, where it blanks those
   rows from the calibrated *output* (see below) — so the calibrated trace
   shows only good air.
4. `cal_mean_points` averages each cal window (offset in seconds relative to
   the interval's last timestamp, independently configurable per bottle) and
   identifies which physical bottle was flowing.
4b. `cal_switch_mask` flags the **switch-over sample**: the first row whose
   solenoid flag has gone False while the cell still holds cal gas. The cause
   is a real ~1 s latency in the serial data — measured across all 45
   transitions in the Feb 2025 file (1 Hz), the response starts moving a
   median of 1 sample after the flag changes, **symmetrically at both edges**
   (rising 0–2 samples, falling 0–2), which is what marks it as a fixed
   transport delay rather than an edge artifact. The rising edge needs no
   handling — cal periods run 40–97 s and the mean windows are anchored to
   `Cal_p` at the end, so even `-30 s` never reaches within 2 s of an
   interval's start (0 of 41 events) — but the falling edge escapes into the
   ambient record. The offending row is the tank at full concentration
   (206.6 ppm against 206.51 assigned, Feb 2025) sitting in the ambient
   record. Three masks all miss it — the cal mask goes by the flag, the
   pressure mask misses it because the transient arrives a sample later, and
   `find_intervals` closes a run *at that row's timestamp* while
   `post_cal_flush_mask` starts strictly after it. It is OR-ed into `not_air`
   (used for shading and for `calibrate_series`' `cal_mask`) and deliberately
   **not** into `cal` itself: `cal` defines the cal intervals, whose ends are
   the `Cal_p` every cal-mean window is measured from, so folding it in there
   would silently shift every cal mean. Verify after touching this:
   `cal_points`, `cal_intervals`, `span_gain`, `loo_rms` and `slope` must be
   identical with the mask and with it stubbed out to all-False.

   It covers **one timestamp, not the whole lag**, and that is a considered
   split rather than an oversight. A 2 s lag leaks a second tank sample
   (16:25:46 = 206.1 ppm on that flight) which this mask does not reach. The
   switch-over mask is the zero-configuration floor that keeps a
   full-strength tank reading out of the air record unconditionally; `Flag
   Air` is the tunable part covering the rest of the lag and the
   cell-clearing tail behind it. Parameterising the lag, or shifting the
   solenoid flags by it, was considered and declined (2026-07-27): shifting
   moves `Cal_p` and therefore every already-saved flight's calibration.
5. `post_cal_flush_mask` ("Flag Air") flags the first N seconds of *ambient*
   data after each cal interval ends — the detector cells are still clearing
   cal gas, so those rows read toward the tank rather than the atmosphere.

**`find_intervals` returns the timestamp of the first row *after* each run as
its `end`**, not the last row inside it. Everything downstream inherits that:
`Cal_p` (so a cal window of `[-15, -1]` ends one sample before the first
ambient row), the flush window's start, and the right-hand edge of the cal
shading. It is why the switch-over sample fell through every mask until
`cal_switch_mask` existed. Changing it would move every cal mean, so leave it
and remember it.

**`calibrate_series`' `calibrated` output is the calibrated GOOD AMBIENT
record.** Three masks are blanked from it as the very last step — `cal_mask`
(rows inside a cal period, plus the switch-over sample: the GUI passes
`not_air`), `flush_mask`, and `exclude_mask` (warm-up +
out-of-spec detector pressure) — so the calibrated trace and the exported
`<col>_cal` column contain good air and nothing else. Nothing is lost by
this: `cal_slope`/`cal_intercept` are still emitted on every row, so the
calibrated value of a blanked row can be recomputed by anyone who wants it
(to check closure, say), and the raw column is untouched. The result exposes
`in_cal`/`flushed`/`excluded` separately, `non_ambient` (= `in_cal | flushed`)
and `blanked` (all three, i.e. exactly the NaN'd rows).

**Blanking is output-only; it must never feed back into the calibration.**
`exclude_mask` has *two independent uses*, and conflating them is the trap
here: passed to `cal_mean_points` it drops raw rows before the cal means are
estimated (so it can remove a cal point), while passed to `calibrate_series`
it only blanks rows of the finished output. `cal_mask` and `flush_mask` have
the output use only — folding either into the `cal_mean_points` call would
corrupt the calibration, obviously so for `cal_mask` (the cal data *is* the
calibration's input), and for the flush because a cal window configured to
reach past `Cal_p` would then start losing points. Verify after touching any
of this, in both directions: `cal_points`, `loo_rms` and `span_gain` must be
identical with the flush at 0 s and at 30 s, **and** identical whether or not
`exclude_mask=` is passed to `calibrate_series`.

Blanked rows stay in the frame as NaN rather than being dropped: `redraw()`
keeps them (`keep = calibrated.notna() | blanked`) precisely so the
calibrated line **breaks** over each gap. Dropping them would make matplotlib
draw a straight segment across the removal and hide it. The raw trace keeps
everything — that contrast is the point of the feature — and the teal flush
band (like the orange/red masking bands) is shaded whether or not "Show
calibrated" is on, since the band is how the user finds out those rows exist
before turning the overlay on.

Every view reads these through `UcatsbGui._analysis_for(gas)` /
`_calibration_for(gas)` / `_uncertainty_for(gas)` — **dicts keyed by gas**,
not one slot for the gas on display, because the Correlations tab needs two
gases' calibrations at once. `_get_analysis()`/`_get_calibration()` remain as
the current-gas shorthands. Their sole invalidation site is `refresh()`, the
entry point every state change calls; it clears all three unconditionally
rather than comparing a composite key of (file, gas, warm-up, tolerance, cal
windows, drift model), since such a key is easy to get subtly wrong and then
serves a stale plot. The recompute is milliseconds against ~100 ms of
rendering.

`_settings_for(gas)` decides *whose* settings a gas is analysed with: the live
control widgets for `current_gas`, `self.config[gas]` for any other. The two
agree in practice (every edit writes `config[current_gas]`), but a
Correlations plot of CO2 against N2O must analyse each with its own saved
warm-up/tolerance/cal windows, not with whatever the panel happens to be
showing. Anything reading the widgets directly instead is a bug for the
non-displayed gas — that is why `_analysis_for` takes settings, not spin-box
values.

### Calibration: drift removal and calibration are ONE step

`calibrate_series` builds a time-varying two-point calibration from the
per-injection means `cal_mean_points` already produces (it does **not**
re-average, and the 4-tuple those come back as is load-bearing — four unpack
sites depend on it). Each bottle's *measured response* is interpolated in
time, and at every ambient timestamp:

    slope(t)     = (A_hi - A_lo) / (R_hi(t) - R_lo(t))
    intercept(t) = A_lo - slope(t) * R_lo(t)
    calibrated   = slope(t) * measured(t) + intercept(t)

Slow drift falls out automatically because the cal responses themselves carry
it — there is no separate detrending pass, and adding one would double-count.
A two-point (slope + intercept) form is required rather than a simple offset:
the measured span error is several percent (span gain 0.96 on the Jul 2026
flight, 1.06 on Feb 2025), so gain must be corrected too.

Non-obvious properties, each of which has bitten a plausible implementation:

- **QC scatter is leave-one-out, not residual-about-the-model.** The default
  `linear` drift model interpolates *through* every node, so its residual at
  each node is identically zero — useless as a quality metric. `loo_residuals`
  predicts each node from that bottle's *other* nodes instead, which also
  cancels genuine slow drift (both neighbours share it). The plain closure
  residual is still computed and plotted, but as a self-consistency check
  that should read ~0 under `linear`; a non-zero value there is a bug signal.
- **The trustworthy region is the INTERSECTION of the two bottles' node
  spans**, not "first to last cal event". A bottle can lose points to masking
  while the other keeps going. In the partial-overlap region one bottle
  interpolates while the other is flat-held, so it is flagged `extrapolated`
  along with the true edges and any anomalously long node gap.
- **Nodes are filtered for serial consistency.** `cal_bottle_series` rejects a
  point whose matched serial disagrees with its state's consensus: its window
  straddled a solenoid transition and it measured the *other* tank. This is
  real — a 374.66 ppm point tagged to the ~206 ppm bottle appears in the Feb
  2025 file with masking off, and would corrupt the calibration either side.
- **`interp_hold` reindexes via a union index** because the datetime column
  contains duplicate timestamps (1435 in one test file). Building a Series
  *from* a duplicated label set raises; reindexing *by* one is fine.
- **numpy is deliberately not used.** The full interpolation is ~9 ms, next to
  ~100 ms of rendering. The pandas path handles the datetime index, duplicate
  labels and the edge-hold policy declaratively.

`cal_mismatch_notes` is **advisory only and must never auto-substitute a
tank.** `cals.yaml` describes the *current* run, so applying it to an older
flight can silently use the wrong tank — it correctly flags `CC302489` on the
Feb 2025 flight. But the same nearest-roster-tank heuristic misfires where the
offset is real gain error (it names `DT0040700` on Jul 2026), so `span_gain`
is reported beside it and the wording only ever asks the reader to check.
`load_cal_roster` (all tanks) exists solely to feed this; bottle *matching*
still uses `load_cal_bottles` (the two plumbed tanks) for the reason in the
section below.

Gases with `has_masking=False` (Ozone) short-circuit to `ok=False` with a
displayable `reason` before any pandas work, as do "no cal events survive the
current masking" and "no assigned value for this gas". A flight with only one
usable bottle degrades to an offset-only correction (`mode="offset"`,
`slope≡1`) rather than refusing, but says so in the header and the legend.

### cals.yaml: full tank roster + a per-run cal0/cal1 assignment

`cals.yaml` (a local copy of `~/code/ucats-b/cals.yaml` — resync by hand if
the acquisition repo's roster changes) has two parts: a top-level block per
serial (the full roster of every cal tank ever used, each with its assigned
mole fraction and uncertainty per gas as flat `<GAS>` / `<GAS>_unc` keys,
e.g. `CO2: 418.947` / `CO2_unc: 0.021`), and a `cals: {cal0: ..., cal1: ...}`
block naming which *two* of those serials are actually plumbed in for the
current run. Bottle *matching* only ever sees two tanks, never the roster:
`select_cal_bottles(roster, serials)` is the single implementation of that
rule, and both entry points go through it — `load_cal_bottles` (cals.yaml's
own `cals:` pair, via `load_cal_assignment`) and the GUI's per-flight
selection. This is deliberate, not an oversight: an unrelated roster tank
could otherwise coincidentally match a measured value more closely and
produce a wrong identification.

`cals.yaml`'s `cals:` block is now only the **default** pairing. Because it
describes the tanks plumbed in *now*, it is wrong for any older flight, so
the GUI's Cal Tanks tab overrides it per flight and stores the choice in
`<dataset>_conf.yaml` (see Config persistence). Editing `cals:` is still what
records a tank swap on the *current* run; it is no longer the only way to
analyse a flight that flew something else.

Not every roster tank has an `info` field (a rough round-number label like
`50%`/`100%` for the original two tanks; the newer ones added don't have an
obvious equivalent) — `_cal_box_title` in `ucatsb_gui.py` treats it as
optional and still shows the serial + mole fraction without it, so don't
reintroduce a hard dependency on `info` being present.

### Cal bottle identity is matched by concentration, not trusted from config

`j_sol_cals` is a raw digital solenoid state (0 or 1) — it is **not** a
reliable bottle identifier. The `cals:` key order in `cals.yaml` does not
necessarily match which serial is actually plumbed to which digital state
on a given flight (confirmed empirically: on the reference flight,
`j_sol_cals==0` measured ~217 ppm CO2, matching `CB09960`'s nominal 206.51
ppm, not `CC302489`'s 418.95 ppm — the reverse of the naive `cal0`→digital-0
reading). `match_cal_serial` instead picks whichever of the two assigned
serials' nominal concentration (for the active gas) is closest to the
measured window mean. This is self-correcting if bottles are swapped
between flights (as long as the flight names the right *pair*, in `cals:` or
in its own conf file) — do not "fix" it to use the `cals.yaml` key order
directly. It also means the Cal Tanks tab's `cal0`/`cal1` combos select a
set, not a wiring: swapping which combo holds which serial changes nothing.

`cals.yaml` has previously contained literal tab characters as
`key:\tvalue` separators, which are invalid YAML syntax. `_read_cals_yaml`
blanks tabs to spaces before parsing rather than erroring — keep that
workaround even though the current file is clean, since a future hand-edit
or re-copy from the acquisition repo could reintroduce them.

### d1 vs d2 detector routing is not fixed across flights

The Aeris `d1` detector measures CO2 and N2O. `d2` is whatever the second
Aeris head happens to be configured for, and that has already changed once:
on the reference (Feb 2025) flight `d2` measured CO plus a redundant N2O
channel; on the July 2026 flight (after a CH4 Aeris was installed) `d2`
instead measures CH4/H2O and has no CO or N2O columns at all. **Don't assume
either layout — check the actual CSV header.** `GASES[gas]["detector"]`
records which detector each gas's main trace and its Detector
Pressure/T_gas aux options come from (`d1` for CO2/N2O, `d2` for CH4);
`aux_trace_info` builds the column name from that instead of a hardcoded
per-gas branch, so adding a new gas is a matter of adding a `GASES` entry
with the right `detector`, not touching `aux_trace_info` itself.

Because the column set genuinely differs between flights,
`UcatsbGui.__init__` reads the CSV header first (`nrows=0`) and only
requests columns that exist, then filters `GASES` down to
`self.available_gases` (only gases whose `value_col` is present) for the Gas
combo box. Downstream code (`redraw`, etc.) still indexes the global `GASES`
dict directly — that's safe because `self.current_gas` is only ever set from
`self.available_gases`' keys, a subset of `GASES`. Don't reintroduce a bare
`usecols=REQUIRED_COLUMNS` read; it will raise if any expected column is
missing from a given file's schema.

Data is plotted **uncalibrated** (`d1_CO2_ppm`, `d1_N2O_ppb`, `d2_CH4_ppb`),
not the `*c_ppm`/`*c_ppb` calibrated columns — this was a deliberate switch;
don't revert to the calibrated columns without being asked. The "Show
calibrated on main plot" toggle does **not** change this: it overlays the
result of *this repo's* `calibrate_series` in red (`CALIBRATED_COLOR`),
keeping the raw trace in its usual blue `LINE_COLOR` underneath at
`alpha=0.55`. The two traces are distinguished by **hue, not by which one is
faded** — recoloring the raw trace when the overlay came on read as the raw
data having changed. It is session-only and defaults off, precisely so
the app never starts up showing calibrated data without the user asking.

### GUI view-preservation (`ucatsb_gui.py` `redraw()`)

`redraw(preserve_view=False)` rebuilds the whole Figure from scratch every
call (`self.figure.clear()` + fresh `add_subplot`/`add_gridspec`) rather than
updating artists in place, because the panel count changes (single axes vs.
main+aux) depending on the "Trace Above" selection. To avoid the zoom
resetting on every masking/averaging tweak or aux-trace change:

- The **main plot's** x/y limits are captured before clearing and reapplied
  after, whenever `preserve_view=True` and a previous draw exists — this
  applies even across aux-panel add/remove/switch, since the user asked for
  the lower plot specifically to hold still through upper-panel changes.
- The **aux panel's** y-limits are only preserved if it's still showing the
  *same* trace as before (tracked via `_last_aux_key`) — a different
  trace has a meaningless old y-range, so it re-autoscales.
- `PlotPane.reset_nav()` (`toolbar.update()` + `push_current()`) resets the
  NavigationToolbar's Home target to the newly-built full-scale view (its
  nav stack otherwise still references the just-destroyed Axes objects).
  Only `on_gas_changed` skips `preserve_view` — switching species changes
  the y-range meaning entirely, so a full rescale there is correct.

### Box stats ("Stats" toolbar toggle)

A checkable `QAction` appended to the stock `NavigationToolbar` (not declared
through `NavigationToolbar.toolitems`, which would require subclassing the
toolbar just to reach back into the pane) drives a `RectangleSelector`; the
box's n/mean/std land in a `QLabel` **outside** the Figure, which is what lets
the readout survive a redraw.

- **The selector must be rebuilt on every draw.** `redraw()` calls
  `figure.clear()`, so a selector from the previous draw holds a destroyed
  Axes and silently stops responding — no error, the tool just goes dead.
  `PlotPane.attach_stats_selector(ax)` is called at the end of `redraw()` for
  exactly this reason, and carries `_box_extents` across so the drawn box
  survives a masking tweak.
- **Pan/zoom hold the canvas widgetlock**, and `_SelectorWidget.ignore()`
  drops every event while it is held, so enabling Stats releases whichever is
  active. The reverse needs no handling: clicking pan afterwards just makes
  the selector inert until pan is switched off.
- **It applies no masking, deliberately.** This is a generic selection tool:
  the box is the user's statement of which data they mean, and its vertical
  bounds already leave the cal dives out when it is drawn around the ambient
  band. An earlier version filtered to ambient+unmasked; that was removed
  because it made a general-purpose tool silently gas-pipeline-specific.
  `n_clipped` (inside the time span, outside the box vertically) is still
  reported — a 2D marquee truncates the distribution and narrows the reported
  sigma, so the count is what makes an accidental clip visible.
- **Traces are registered as they are plotted** (`_register_stats_trace`),
  not scraped back off the Axes: artists carry no units and no stable
  identity, and a conditionally-drawn trace would be easy to miss. The
  registry is rebuilt by every `redraw()` and feeds the combo box; the
  combo's current key is preserved across repopulation when it still exists.
- **`y0`/`y1` are only meaningful on the Axes the box was drawn in.** When
  the chosen trace lives elsewhere (the aux right axis, or a main trace when
  the box is in the upper panel), `box_stats` is called with no y-bounds and
  the readout says so. Without this the y-range of one panel would silently
  filter another panel's data — and even within one Axes it bites: the
  calibrated overlay sits an intercept away from the raw trace (~10 ppm on
  the Jul 2026 CO2 flight), so a box drawn snugly around raw legitimately
  contains almost no calibrated points, which the `outside box vertically`
  count explains rather than hides.
- **One box at a time, tracked by Axes *index*, not identity** — the Axes
  objects are destroyed on redraw, so `PlotPane._box` stores
  `(index into the selector list, extents)`. Selecting in one panel hides the
  other panel's box, so the readout always refers to a box that is on screen.

### Tabs (`PlotPane`, `refresh()` dirty dispatch)

The two views are `PlotPane` instances (Figure + canvas + toolbar) in a
`QTabWidget`. `self.figure`/`self.canvas`/`self.toolbar` stay **bound to the
timeseries pane** rather than being renamed, so `redraw()`'s body needs no
changes; don't rename them while also changing behavior, or a rendering
regression becomes indistinguishable from a refactor slip.

The controls panel stays **outside** the tabs, in a `QStackedWidget` with two
pages: the per-gas panel shared by Timeseries/Calibration/Cal Tanks, and the
Correlations panel. `on_tab_changed` switches the page. Sharing one panel is
still the default and the reason holds (duplicating the gas selector, or
making one tab depend on state invisible from another); Correlations is the
exception because it is inherently about *two* gases, so a per-gas panel
beside it would be actively misleading about what is plotted.

`refresh()` redraws only the visible pane and marks the others dirty, so a
spinbox drag doesn't render a pane nobody is looking at. The `_preserve`
latch matters: a requested full rescale must survive until it is actually
honoured, or changing gas and then nudging a spinbox would leave the
never-drawn pane stuck at a stale scale when first opened. The cal tab's own
"same content" key is `(gas, drift_model, smooth_events)`; the correlation
tab's is `(x_gas, y_gas)`.

Correlation-only controls call **`_refresh_corr`, not `refresh()`**: marker
size, the error-bar toggle and the axis pickers change nothing about any
mask, calibration or setting, so going through `refresh()` would clear both
per-gas caches and dirty every pane — making a marker-size nudge recompute
two gases' analyses and redraw the timeseries.

### Config persistence: explicit save only

Settings live in one block per gas (`warmup_min`, `pressure_tol_mbar`,
`flag_air_s`, `cal1_window_s`, `cal2_window_s`, `drift_model`,
`drift_smooth_events`) plus the flight's `cals:` pairing, in a YAML file
beside the CSV. **Nothing is written without `on_save_clicked`.** The
requirement driving this: open a saved analysis, experiment, quit — and the
file is exactly as it was. So:

- `_load_flight_config` **writes nothing** (the old `_adopt_flight_config`
  created the file at load time; that is gone). It picks a config via
  `_config_candidates` — any `<dataset stem>*.yaml` beside the CSV, default
  name first — applying one silently and asking via `_choose_config_file`
  when there are several, with a "start from defaults" option.
- **There is no template.** `load_config(None)` gives `DEFAULT_GAS_SETTINGS`,
  and `ucatsb_gui_config.yaml` now holds *only* `recent_files`. Seeding one
  flight's settings from another silently is what per-dataset configs exist to
  prevent, and the same argument already applied to the tank pairing.
- `config_path` is the name Save offers; `config_loaded_from` is the file
  actually opened, or None when the dataset started from defaults. They differ
  exactly in that case, which is what the Cal Tanks readout reports.
- Save is **save-as every time**: several configs per dataset was the point,
  so the filename is always offered for editing rather than overwriting what
  was opened.

**Dirty state is a comparison, not a flag.** `_is_dirty()` deep-compares
`_current_state()` (config + cal_selection) against `_saved_state`, snapshotted
by `_snapshot_state()` at load and after a successful save — the only two
places. A flag would need clearing in as many places as it is set and would
still be wrong when a change is *undone*; the comparison reports "nothing to
save" when a spin box goes back to where it started. `_mark_dirty()` therefore
only refreshes the button, and `_confirm_discard()` prompts (Save / Discard /
Cancel) from `closeEvent`, `on_load_data_clicked` and `_open_recent_now`.

**`save_config` writes a *fresh* document, so an omitted argument is a
deletion.** `cals:` goes only to a dataset config; `recent_files:` only to the
app-level file (`_save_app_config`, which passes `{}` for the gas blocks
because that file no longer carries any). An app-level write that forgets
`recent_files=` wipes the list — the same trap `_controls_to_settings()` has
for per-gas keys.

`_initializing`/`_loading` flags exist specifically to suppress redraw/save
during programmatic widget setup (e.g. `setChecked` on a freshly-constructed
radio button fires its signal immediately, before sibling widgets it might
depend on exist yet) — keep that guard pattern when adding new controls.

**Adding a persisted setting requires four edits, not one.**
`on_control_changed` assigns `self.config[gas] = self._controls_to_settings()`
— a *fresh* dict — so a key missing from `_controls_to_settings()` is
silently dropped on the next control change (and from the next save), even
though `load_config`'s `.update()` appeared to preserve it. Touch all of:
`DEFAULT_GAS_SETTINGS`, `_controls_to_settings()`,
`_apply_settings_to_controls()`, and the `setEnabled(has_masking)` list in
`_select_gas()`. (A control added *inside* an existing group box — as
`flag_air_spin` was, inside `mask_box` — inherits that box's `setEnabled`, so
only the first three apply; adding a new group box is what requires the
fourth.)

**"Copy settings to all gases"** (`on_copy_masking_to_all`) writes
`COPIED_SETTING_KEYS` — `warmup_min`, `pressure_tol_mbar`, `flag_air_s`,
`cal1_window_s`, `cal2_window_s` — from the live controls into every *other*
gas block. A new persisted control has to be added to that tuple too (or
consciously left out), or the button will silently leave it behind on the
other gases. Only `drift_model`/`drift_smooth_events` are excluded: those are
a judgement about one gas's cal-record noise rather than a description of the
flight. Ozone drops out by construction rather than by name — it has no entry
in `self.config` at all — so "every gas in the config" stays the right set if
a gas is ever added.

**The copied values are `deepcopy`d per target gas.** `cal1_window_s` /
`cal2_window_s` are lists: handing every gas the same list object makes
`yaml.safe_dump` write anchors (`&id001`/`*id001`) into the conf file, and
makes an in-place edit of one gas's window silently change the others.

New masking settings default to a **no-op value** (`flag_air_s: 0`) rather
than a physically plausible one. `load_config` fills missing keys from
`DEFAULT_GAS_SETTINGS`, so a non-zero default would silently change the
output of every already-saved config on first launch after the upgrade.

**Verification scripts must pass `config_path=` to a scratch file *and* load
a scratch copy of the CSV.** `UcatsbGui` writes the real
`ucatsb_gui_config.yaml` (the recent-files list) whenever a dataset loads, so
point that at a scratch path. Settings themselves are no longer written
without Save, but a test that exercises Save writes a config beside the CSV —
so copy the CSV into the scratchpad rather than pointing the test at
`~/Data/UCATSb/...`. Drive Save by monkeypatching
`QFileDialog.getSaveFileName`, and the config chooser by patching
`UcatsbGui._choose_config_file`; both are what `check_save.py` does.

### Correlations tab and `calibration_uncertainty`

A tracer-tracer scatter, calibrated wherever a calibration exists — and that
is the point rather than a convenience: a slope from uncalibrated counts
carries that detector's gain error into the slope, which is the number the
plot produces. A point survives only where **both** axes have a value, so it
is the intersection of two independently-masked series.

`_corr_axis(gas)` is the single place that decides what an axis contributes,
returning `(values, sigma, unit, qualifier, reason)`. Two kinds:

- a cal-bottle gas → `calibrate_series`' `calibrated`, plus a 1σ, already
  masked to good air;
- **Ozone** → `oz_o3best` with its below-floor readings removed (see below),
  and **no sigma** (nothing to propagate — an axis with no calibration passes
  `None` to `errorbar` rather than a zero bar, which would claim a precision
  nobody established). Its partner axis still imposes its masking through the
  intersection.

`reason` is non-None only when a gas that *should* have a calibration lacks a
usable one; Ozone having none is what Ozone is, not a failure to report. The
qualifier rides on each **axis label** (`calibrated` / `oz_o3best`) because
with one axis of each kind a single title word would have to lie about one of
them; the title falls back to "see axis labels" when they differ. Anything
new that assumes both axes are calibrated — the extrapolated-fraction note,
the Flag Air note, the median-1σ readout — has to check `has_masking` first.

`calibration_uncertainty(result)` (shared module) returns the 1σ on each
calibrated value, propagated by writing the two-point calibration as a blend
of the assigned values, `c = (1-f)A_lo + f A_hi` with
`f = (c - A_lo)/(A_hi - A_lo)`. Non-obvious properties:

- **`f` comes from the calibrated value, not from the interpolated
  responses.** That keeps the function a pure function of what
  `calibrate_series` already returns, with no need to re-expose `R_lo`/`R_hi`.
- **`f` is not clamped to [0, 1].** Ambient air usually sits outside the
  bracket the two tanks span, and the uncertainty genuinely grows as the
  calibration extrapolates away from them.
- **Per-bottle response sigma is `sqrt(loo² + closure²)`, and both halves are
  required.** Under `linear` the closure residual is identically 0 (the model
  interpolates through every node), so closure alone would report zero
  uncertainty; under `constant` the LOO scatter understates the error the
  model introduces by not passing through the nodes. Together they are what
  makes the answer depend on the drift model at all — CO2 on the Feb 2025
  flight moves from 2.67 to 3.34 ppm median 1σ between `linear` and
  `constant`. A change here that makes the two models agree is a bug signal.
- **A missing `<GAS>_unc` counts as 0**, which understates — but inventing a
  number would be worse, and the tank roster is where the fix belongs.
- **Single-sample instrument noise is deliberately absent**: nothing in the
  calibration constrains it. So is any inflation over extrapolated spans,
  which are reported as a percentage of points instead — folding them in
  would disguise "we are guessing here" as ordinary uncertainty.

**The z-axis coloring** (`CORR_COLOR_BY`, `CORR_COLORMAP`) maps each entry to
`(label, column, colorbar label)`, with `column=None` meaning the time axis —
which is not a plottable column and needs `mdates.date2num` going in and a
`DateFormatter` on the colorbar coming out. Three rules:

- **The z variable joins the pairing rule.** A point with no value for it is
  dropped, not drawn: matplotlib would paint it in the colormap's "bad" color
  (transparent), leaving it in the fit but invisible on the plot.
- **Only encodings the loaded CSV can supply are offered** — the schema
  differs between flights, so `oz_p` may not exist. `_populate_corr_color_combo`
  also clears a selection the new file cannot honour.
- **`turbo`, not `jet`/`rainbow`.** Same rainbow ordering, monotonic
  luminance, so it doesn't manufacture bands of false structure. A colorbar is
  always drawn with it; a continuous color encoding with no scale is
  unreadable.

**The OLS fit is opt-in** (`corr_fit`, default False). A straight line through
a tracer-tracer plot with real structure in it describes almost none of that
structure, so leaving it on by default put an authoritative-looking slope on
every plot that mostly measured how the flight was flown. When off, neither
the line, the figure note nor the slope/intercept/r block in the side panel is
produced — `linear_fit` is not even called.

Because these uncertainties are near-entirely *systematic* (they shift a whole
flight together), `linear_fit` is plain OLS and must stay that way: weighting
by them would not do what weighting is for, and would make the reported slope
depend on whether a display toggle is on.

The figure calls out **Flag Air = 0** by name. Post-cal flush points trace a
line from the tank's composition to the atmosphere's, which on this plot looks
exactly like a tracer-tracer correlation and drags the fit (Feb 2025: n 11519
→ 10088 and slope 0.906 → 0.883 when Flag Air goes 0 → 45 s). It is the one
masking setting whose absence is invisible here without saying so.

### Recent files (Load Data menu + File menu)

`self.recent_files` (newest first, capped at `RECENT_FILES_MAX`) is fed by
`_remember_recent`, called from `load_csv` **after the dataset is committed** —
`load_csv` validates first, so a file that fails never enters the list, and
because `main()` goes through it a CLI-opened dataset is remembered too.

`_rebuild_recent_menus` populates the Load Data button's menu and
`File > Open Recent` from one builder, so the two cannot drift. It rebuilds
wholesale rather than patching, for the same reason `refresh()` invalidates
unconditionally. The `QAction` for "Open…" is *shared* by both menus.

**The rebuild deletes the QAction whose signal is running.** `_open_recent`
and `on_clear_recent_files` are triggered from actions that the work destroys
(`QMenu.clear()` deletes them) while `triggered` is still unwinding — a
use-after-free on the sender. Both therefore defer via
`QTimer.singleShot(0, ...)` and let the signal finish first. Any future action
inside these menus that ends in a rebuild needs the same treatment; a
verification script must `processEvents()` after triggering one.

Missing files stay listed but disabled and marked `(missing)` rather than
being pruned on sight — an unmounted volume comes back. They are dropped only
when opening one actually fails (`_try_load(..., forget_on_failure=True)`).

### `valid_min`: a physical floor, not a maskable setting

`GASES["Ozone"]["valid_min"] = O3_VALID_MIN_PPB` (−15 ppb) and
`GASES["H2O"]["valid_min"] = H2O_VALID_MIN_PPM` (−5 ppm) drive
`below_floor_mask`, via `UcatsbGui._rejected_mask(gas)`. Declared per gas
rather than special-cased by name, so any gas that gains a floor gets the same
raw/filtered treatment for free — and *not* a control, because it is a
statement about when the sensor is faulting rather than a judgement the user
should be tuning per flight. The H2O floor is precautionary: the Jul 2026
flight's minimum `w_H2Obest` is +14 ppm, so it currently removes nothing.

The two floor gases are also the two with `has_masking=False` — no cal
bottles, so no warm-up/pressure/cal-window machinery and no calibration.
Anything iterating gases must not assume a cal-bottle gas: check
`has_masking` (see the Correlations tab's `_corr_axis`).

- **The floor is well below zero deliberately.** Real near-zero ozone scatters
  negative; the Feb 2025 flight has 168 readings in −15..0 ppb (noise about a
  small true value) against 7 faults reaching −2292. Clipping at zero would
  bias the low end upward and hide the instrument's true noise. Don't "tidy"
  it to 0.
- **NaN is not flagged.** Absent and invalid get different treatment
  downstream — a gap versus a removal.
- `_rejected_mask` is deliberately **not** part of `_analysis_for`'s cache: it
  depends on nothing the user can change, and it applies to a gas that has no
  analysis settings at all.
- The timeseries draws raw (blue, faded) + filtered (red), reusing
  `LINE_COLOR`/`CALIBRATED_COLOR` so red means the same thing here as on the
  calibrated overlay. The filtered line keeps removed rows as NaN so it
  **breaks** over them, the same trick `calibrate_series` output uses.
- **The default y-range is framed on the filtered series** when anything was
  removed (and only when not preserving a view). A single −2292 ppb fault
  otherwise sets the scale and squashes the real record into the top fifth of
  the axes, which would make masking it pointless on the figure that most
  needed the help. The raw trace still runs off-scale, and the note says how
  many readings that is.

### Cal Tanks tab

Its own tab (`_build_cal_tanks_pane`), not another control-panel group box:
the pairing is one per *flight*, while every control in that panel is per gas.
Two combos over the whole roster (`load_cal_roster`), not just the plumbed
pair — picking a tank outside the current `cals:` block is the entire point.

- `self.cal_bottles` is **derived state**: `_rebuild_cal_bottles` recomputes
  it from `cal_selection` via `select_cal_bottles`, which is the one place
  implementing "matching may only see the plumbed tanks" (see the section on
  `cals.yaml` above — `load_cal_bottles` now goes through it too). Never
  assign `cal_bottles` directly.
- `cal0`/`cal1` name a *set*, not a wiring. `match_cal_serial` identifies the
  tank in each window by measured concentration, so swapping the two combos
  changes nothing — said out loud in the tab's tooltip because the labels
  invite the opposite assumption.
- A tank change refreshes with **`preserve_view=False`**, like a gas change
  and unlike every other control. The Calibration tab's top panel plots
  measured *minus assigned*, so a different tank shifts it by the difference
  in assigned values (11 ppm between `CC470901` and `CC302489`) and the old
  y-range leaves the new points off-scale.
- `_draw_current_tab` dispatches on the current *widget*, not the tab index:
  Cal Tanks draws nothing, and the old `index == 1` test would have made it
  redraw the timeseries pane.
