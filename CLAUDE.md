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

There is **no second entry point**. `ucatsb_analysis.py` is a library and is
never run directly; the CO2-only `plot_co2_timeseries.py <csv>` CLI it used to
carry was removed (2026-07-29) because the GUI supersedes it and it read the
cal pairing from `cals.yaml`'s `cals:` block, mislabelling tanks on any flight
that did not fly the currently-plumbed pair, with no way to override.

No test suite or CI. Changes are read-verified by importing the module under
`QT_QPA_PLATFORM=offscreen` and grabbing a screenshot rather than by unit
tests — there isn't a headless-safe way to assert on plot pixels, so manual
visual review (via `Read` on the saved PNG) is how changes here get checked.

## Architecture

`ucatsb_analysis.py` is the shared logic module — masking, cal detection,
calibration and its uncertainty, and both export writers — and is **Qt-free
on purpose**, so a batch script or notebook can reuse `calibrate_series` or
`export_icartt` without a GUI toolkit. Don't import PyQt5 into it, and don't
move data logic into `ucatsb_gui.py`, which is the GUI and nothing else.
(It was `plot_co2_timeseries.py` until 2026-07-29; the name came from a
standalone CO2 figure CLI that was removed in the same change.)

`ucatsb_gui.py` imports from it rather than duplicating: `drop_presync_rows`,
`find_intervals`, `merge_close_intervals`, `shade_intervals`,
`cal_mean_points`, `load_cal_roster`, `load_cal_assignment`,
`select_cal_bottles`, `most_common_serial`, `mean_std_label`,
`CALS_YAML_PATH`, plus the calibration functions below and the two export
writers (`export_companion_csv`, `export_icartt`, with
`icartt_time_base`/`icartt_filename`/`DEFAULT_ICARTT_META`). Any change to
masking/cal-detection/calibration behavior belongs in that shared module.

**The shared palette lives in `ucatsb_analysis.py` and is imported, not
redeclared.** Thirteen constants (`LINE_COLOR`, `CAL0_COLOR`, `MUTED_COLOR`,
`D1_P_TARGET_MBARS`, …) were once defined in *both* files with identical
values — agreeing by discipline alone, so an edit to one silently made the
calibration panels and the timeseries disagree about what a color meant.
Only `CALIBRATED_COLOR` and `STATS_BOX_COLOR`, which nothing in the analysis
module draws, are declared in the GUI.

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
3. Warm-up (first N minutes), end-of-flight (last M minutes, measured back
   from the record's last timestamp rather than a clock time), out-of-spec
   detector pressure
   (`|d1_P_mbars - 140| > tol`) and — when **Pumps on** is set —
   pumps-off (`j_pumps != 1`) masks are computed and **applied to the raw
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
`<gas>_cal` column contain good air and nothing else. Nothing is lost by
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
still sees only the two plumbed tanks, for the reason in the section below.

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
rule and the only way to build a `cal_bottles` dict. This is deliberate, not
an oversight: an unrelated roster tank could otherwise coincidentally match a
measured value more closely and produce a wrong identification.

A `load_cal_bottles` helper (roster + `cals:` pair, composed) used to sit
beside it as a second entry point. It was deleted with the CO2 CLI that was
its only caller — the GUI reads the pairing per flight, so the `cals:` block
is a *default* that nothing loads directly any more.

`cals.yaml`'s `cals:` block is now only the **default** pairing. Because it
describes the tanks plumbed in *now*, it is wrong for any older flight, so
the GUI's Cal Tanks tab overrides it per flight and stores the choice in
`<dataset>_conf.yaml` (see Config persistence). Editing `cals:` is still what
records a tank swap on the *current* run; it is no longer the only way to
analyse a flight that flew something else.

`info` is a rough round-number label (`50%`/`100%`) used only to name the
cal-window rows. Every tank in the roster carries one as of 2026-07-29, but
it stays **optional**: `_set_cal_row_label` in `ucatsb_gui.py` falls back to
the serial rather than losing the label, so don't reintroduce a hard
dependency on `info` being present — the next tank added to the roster may
well arrive without one. `info` is now the *visible* text (a row label has far
less room than the group-box title this used to be, 2026-07-30); the full
"100% Cal (CC302489) 418.947 ppm" moved to the label's tooltip, which is built
in the same place so the two cannot drift.

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
don't revert to the calibrated columns without being asked. The **Calibrated**
toolbar toggle does **not** change this: it overlays the
result of *this repo's* `calibrate_series` in red (`CALIBRATED_COLOR`),
keeping the raw trace in its usual blue `LINE_COLOR` underneath at
`alpha=0.55`. The two traces are distinguished by **hue, not by which one is
faded** — recoloring the raw trace when the overlay came on read as the raw
data having changed. It is session-only and defaults off, precisely so
the app never starts up showing calibrated data without the user asking.

It lives on the Timeseries **toolbar** (`PlotPane.calibrated_action`), not in
the controls panel where it used to be a checkbox (moved 2026-07-30): it
changes what the figure draws and nothing else, which is what everything else
on that toolbar does, and sitting among the per-gas settings made a
session-only view toggle look like something that gets saved. Two consequences
worth keeping: it is hidden on the Calibration and Correlations panes, which
have no raw trace to overlay; and because it no longer sits inside `cal_box`
it does not grey out with it, so `_select_gas` disables it explicitly for a
gas with `has_masking=False`.

### GUI view-preservation (`ucatsb_gui.py` `redraw()`)

`redraw(preserve_view=False)` rebuilds the whole Figure from scratch every
call (`self.figure.clear()` + fresh `add_subplot`/`add_gridspec`) rather than
updating artists in place, because the panel count changes (single axes vs.
main+aux) depending on the "Above:" trace selection. To avoid the zoom
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

### The controls panel's size is a constraint on the whole window

Both stack pages are `CONTROLS_WIDTH` wide and as tall as their contents, and
the stack is as tall as the **taller** page — currently Correlations, not the
per-gas panel, so a height measurement that only looks at the one you edited
is measuring the wrong thing.

**The stack sits inside a `QScrollArea`** (added 2026-07-30). Without it the
panel's ~1090 px sizeHint set a hard floor of ~1130 px on the *window*, which
does not fit a laptop screen at default scaling — and every control added over
the years pushed it up with nothing to say so. With it the window minimum is
~420 px. Compaction in the same change brought the panel itself to ~715, so
the scrollbar normally never appears; the scroll area is what keeps that a
layout detail rather than something the user hits.

Two failure modes to know, because neither announces itself:

- **Too little height is invisible; too little width is silent corruption.**
  A panel taller than the viewport now scrolls. A panel *wider* than
  `CONTROLS_WIDTH` neither scrolls (the horizontal bar is off by policy) nor
  wraps — it clips its own right-hand edge, and buttons quietly lose their
  right halves. `CONTROLS_WIDTH` is 312 against a widest group box of 292 for
  exactly that headroom. If you add a wide control, measure
  `box.minimumSizeHint().width()`, don't eyeball it.
- **A `QSpinBox` asks for ~85 px whatever it holds.** Rows carrying two of
  them (the trim pair, each cal window) overflow on that alone. They are
  capped with `setMaximumWidth`, which works because Qt's `qSmartMinSize`
  bounds the minimum by the maximum — `setMinimumWidth` would not have done
  it, and the natural size hint ignores that these only ever show three
  characters.

Group boxes are the expensive thing to add: each costs ~41 px of title and
margin before it holds anything, which is why Gas and the aux pickers share
"Traces", and the two cal-window boxes plus the drift model share
"Calibration". Layouts here run 6 px margins / 4 px spacing rather than Qt's
9/6 default for the same reason.

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
  and neither app-level file holds analysis settings at all — only `icartt`
  (see the Export tab) and `recent_files`. Seeding one flight's
  settings from another silently is what per-dataset configs exist to prevent,
  and the same argument already applied to the tank pairing. The ICARTT
  metadata is the deliberate exception, and is not a *setting*: it describes
  the campaign, not the analysis.
- `config_path` is the name Save offers; `config_loaded_from` is the file
  actually opened, or None when the dataset started from defaults. They differ
  exactly in that case, which is what the Cal Tanks readout reports.
- **The file label under Load Data names both**, as two fixed lines —
  `Filename: <dataset>.csv` / `Config: <config>`, the latter `defaults` when
  `config_loaded_from` is None. A flight can have several configs, so which one
  is in force must be answerable without opening the Cal Tanks tab.
  `_update_file_labels` writes **both** control pages (the per-gas panel and
  Correlations) from one place, like `_rebuild_recent_menus`, and is called
  from `load_csv`, `_apply_config_file` and `on_save_clicked` — the three
  points where either name can change. It shows the config *filename* rather
  than the variant recovered from it: the variant is a name for a thing on
  disk, and where the two can disagree (a config predating the scheme, or one
  opened from another directory) the filename is what tells you where to look.
  Rich text, so the field names are bold and the values ordinary weight;
  `html.escape` on both, since a file name is not our markup.

  **Two lines and no word wrap means the values must be elided** — this is the
  panel-width failure mode landing on the one label whose job is to name a
  file. `Config: <default name>` measures 305 px and a short variant 349,
  against ~300 px of panel, so unwrapped they lose their right-hand ends in
  silence. `_elide_field` fits the value to `CONTROLS_WIDTH - 2 *
  CONTROLS_MARGIN` minus the bold field name, and the tooltip keeps both full
  paths. It elides from the **left**, against the usual convention: the config
  name opens with the dataset stem, which the Filename line right above already
  gives, so the front is the redundant part while the tail — the variant, and
  `_conf.yaml` — is what distinguishes one config from another.

  It measures against `CONTROLS_WIDTH`, not the label's own `width()`, because
  it runs during `load_csv` before the panel is laid out (`width()` is still
  the default 100 there). The panel is fixed-width, so the static budget is the
  true one at every moment — but that also means `CONTROLS_MARGIN` is now load
  bearing in two places, the layouts' `setContentsMargins` and here.
- Save is **save-as every time**: several configs per dataset was the point,
  so the name is always offered for editing rather than overwriting what
  was opened.

**Save asks for a *variant name*, not a filename** (`_choose_config_name`,
2026-07-30). `flight_config_path(csv, variant)` composes
`<dataset stem>[_<variant>]_conf.yaml` beside the CSV;
`config_variant_name` is its inverse, recovering the variant to pre-fill the
field (and returning `""` for a name that predates the scheme).

The dataset stem is **structural, not decoration** — `_config_candidates`
finds a flight's configs by globbing it — and the old `QFileDialog` had no way
to say so, so a config saved as `test.yaml` was silently invisible from the
next open onward. Composing the name on this side is what makes "saved" and
"findable" the same thing. Three consequences:

- **`_config_candidates` stays looser than what Save now writes.** It matches
  any `<stem>*.yaml`, not just `*_conf.yaml`, because configs named freely
  under the old dialog exist on disk and tightening it would strand them.
- **The overwrite prompt is ours now.** `QFileDialog` gave it for free; the
  dialog therefore `exec_()`s in a **loop**, so declining to replace returns to
  the name field instead of cancelling the save.
- **Saving to another directory is gone, deliberately.** Such a config was
  already unreachable by the stem search; `on_load_config_clicked` remains the
  way to open one from anywhere.

`sanitize_config_variant` strips path separators and `..` rather than escaping
them — the variant is free text that lands straight in a path — and drops a
trailing `conf` component so copying an existing name's shape doesn't yield
`..._test_conf_conf.yaml`. It is anchored (`(?:^|_)conf$`) so a variant like
`reconf` survives.

**Dirty state is a comparison, not a flag.** `_is_dirty()` deep-compares
`_current_state()` (config + cal_selection) against `_saved_state`, snapshotted
by `_snapshot_state()` at load and after a successful save — the only two
places. A flag would need clearing in as many places as it is set and would
still be wrong when a change is *undone*; the comparison reports "nothing to
save" when a spin box goes back to where it started. `_mark_dirty()` therefore
only refreshes the button, and `_confirm_discard()` prompts (Save / Discard /
Cancel) from `closeEvent`, `on_load_data_clicked` and `_open_recent_now`.

**`save_config` writes a *fresh* document, so an omitted argument is a
deletion.** `cals:` goes only to a dataset config; `recent_files:` and
`icartt:` only to the app-level file (`_save_app_config`, which passes `{}`
for the gas blocks because that file no longer carries any). That file now has
**two** blocks, and `_save_app_config` passes both every time: a write that
forgets `recent_files=` wipes the list, and one that forgets `icartt_meta=`
wipes the metadata — the same trap `_controls_to_settings()` has for per-gas
keys.

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

`require_pumps` (the **Pumps on** toggle) joins `exclude_mask`, so like the
warm-up and pressure masks it drops cal points as well as blanking the output
— pumps-off data is not ambient air at all. Three properties worth keeping:
a missing `j_pumps` value counts as pumps-off (an unknown pump state is not
evidence the pumps were running); a file whose schema lacks the column
disables the control outright, in `load_csv`, and Qt restores that child's own
enabled state when `mask_box` is re-enabled; and it **must** default off,
because a lab/bench run is pumps-off end to end (the 2026-07-26 file is 100%)
and enabling it there excludes every row and every cal point.

Warm-up and end-of-flight are OR-ed into `trimmed`, which is what gets shaded
and what the note describes — deliberately one orange band and one line
naming whichever ends are active, since both say "the instrument was not doing
what the rest of the flight was doing". `analysis` still exposes `warmup` and
`end_flight` separately for anything that needs to tell them apart.

New masking settings default to a **no-op value** (`flag_air_s: 0`,
`require_pumps: False`, `end_flight_min: 0`) rather than a physically plausible one. `load_config` fills missing keys from
`DEFAULT_GAS_SETTINGS`, so a non-zero default would silently change the
output of every already-saved config on first launch after the upgrade.

**Verification scripts must pass `config_path=` *and* `state_path=` to scratch
files *and* load a scratch copy of the CSV.** `UcatsbGui` writes the real
`.ucatsb_gui_state.yaml` (the recent-files list) whenever a dataset loads, so
point that at a scratch path; `config_path` matters less now that only Save
defaults writes it, but a test exercising that button would edit a *tracked*
file. Settings themselves are no longer written
without Save, but a test that exercises Save writes a config beside the CSV —
so copy the CSV into the scratchpad rather than pointing the test at
`~/Data/UCATSb/...`. Drive Save by monkeypatching
`UcatsbGui._choose_config_name` (return a path, or None for cancel), and the
config chooser by patching `UcatsbGui._choose_config_file`.

**Patch `_choose_config_file` *before* constructing `UcatsbGui`, not after.**
It `exec_()`s a modal dialog, which under `QT_QPA_PLATFORM=offscreen` blocks
forever with no output — and it fires from `load_csv` during `__init__`
whenever the scratchpad already holds more than one config for that CSV, which
a previous run of the same script will have left there. The symptom is a script
that simply hangs.

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

**The median 1σ is reported in two places from one source**, `_median_sigmas`:
the figure's note block and the numbers panel. It is deliberately **not** gated
on the "Error bars (1σ)" toggle — it describes the calibration, not the
drawing — and it is written with a `±` because the `mean ± std` line directly
above it in the panel is the *spread of the atmosphere across the flight*, an
entirely different quantity that on this instrument is an order of magnitude
larger (N2O on 2026-07-28: sd 12 ppb, 1σ 0.815 ppb). The number belongs on the
figure as well as the panel because at flight scale the bars overlap into a
band whose width cannot be read off the axes, and the figure is what gets
saved and handed to someone else.

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

**Picking the dataset that is already open re-runs the config choice**
(`_reopen_config`), rather than returning silently as it used to. A flight with
several saved configs is the case where the chooser matters most, and the
silent no-op left it reachable only by loading another file and coming back.
The CSV is **not** re-read — nothing about the data has changed, and re-parsing
would discard `raw_df` (and the row count `load_flagged` validates against) to
rebuild it identically. Everything after the read is redone: the same
`_load_flight_config` → `_select_gas` → `_update_tank_readout` → `refresh()`
sequence `on_load_config_clicked` uses. `_confirm_discard` still guards it, so
a reopen cannot silently drop unsaved settings.

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

### Manually flagged points ("Flag" toolbar toggle)

The manual override for readings no rule catches — the ozone spikes on the Jul
2026 flight reach 3500 ppb against a 0–900 ppb record, and no threshold
separates them from real data. Drag a box on the Timeseries tab to flag,
right-drag to unflag; the rows leave the calibrated (or filtered) series, both
export products and the Correlations tab, while the raw trace keeps showing
everything.

**Stored as run-length-encoded RAW-file row ranges** in a `flagged:` block in
the flight's `<dataset>_conf.yaml`, `{gas: {rows: N, ranges: [[lo, hi], ...]}}`,
inclusive at both ends. Row ranges rather than timestamps or the drawn
rectangles, for three reasons that each killed an alternative:

- **Timestamps do not identify a row** — the Jul 2026 file has 1294 duplicate
  timestamps, so flagging one spike would silently flag its twin.
- **Storing the rectangle makes the flagged set drift.** A box resolved
  against calibrated values catches different rows after a drift-model or
  cal-tank change; row numbers mean the same thing forever.
- **Set arithmetic.** Unflagging is `subtract_ranges`, splitting a range it
  lands inside, with no paint/erase ordering to get wrong. One gesture is one
  entry however wide, so "hundreds or thousands of points" is a handful of
  lines of YAML.

The helpers live in `ucatsb_analysis.py` (`merge_ranges`, `add_ranges`,
`subtract_ranges`, `ranges_to_mask`, `ranges_row_count`). **`merge_ranges`
merges *adjacent* runs as well as overlapping ones**, so a set of rows has
exactly one representation — `_is_dirty()` deep-compares this, and two gesture
orders reaching the same rows must not read as a change.

`rows:` is a **tripwire, not an index**: raw row numbers only mean anything for
the file they were drawn on, so a regenerated CSV of a different length would
shift every flag silently. `load_flagged` still applies them and warns, because
the user is better placed to judge — but silently wrong is the one outcome this
feature cannot have.

**Flags reach the pipeline by two different routes, and both are needed:**

- **Cal-bottle gases** — OR-ed into `exclude_mask` in `_analysis_for`. That one
  line covers both of that mask's documented uses, since it is already handed
  to `cal_mean_points` (dropping raw rows before the cal means, so flagging a
  bad injection genuinely changes the calibration — verified: one flagged
  injection moves CO2's `span_gain` 0.9571 → 0.9645 and drops a cal point)
  *and* to `calibrate_series` (blanking the output).
- **Floor gases (Ozone, H2O)** — never reach `calibrate_series` at all, so
  `_removed_mask(gas)` = `_rejected_mask | _flag_mask` is what the filtered
  trace, `_corr_axis` and the export block use. `_rejected_mask` stays
  floor-only so the note can say *which* removal a row belongs to.

Because flags sit inside `exclude_mask`, **every edit must go through
`refresh()`** — `_analysis_for` is cached per gas and that is its only
invalidation site. `_after_flag_change` is the single funnel.

Non-obvious behaviours, each deliberate:

- **Flagging matches the box against the RAW value column**, never the
  displayed trace, so a flag names the same rows after any recalibration. The
  cost is that a box drawn tightly around the calibrated overlay (an intercept
  away — ~10 ppm on CO2) can catch nothing, which is why an empty box says so
  rather than doing nothing visible.
- **Unflagging ignores the y-bounds entirely** and clears the whole time span.
  Not leniency: the default y-range is framed on the *filtered* series
  precisely so one 3500 ppb spike does not squash the record, which puts the
  flagged point off-screen — there would be no box to draw around a value you
  cannot see. It also makes unflagging total, and sidesteps the NaN asymmetry
  below.
- **NaN rows are never flagged** (`between` is False for NaN), the same
  "absent is not invalid" rule `below_floor_mask` follows. A box over a span
  containing gaps therefore stores several ranges, not one.
- **The Flagging box is deliberately absent from the
  `setEnabled(has_masking)` list** in `_select_gas`. Ozone is the gas this
  exists for and is exactly the one with no masking settings to enable.
- **"Apply to all gases" spreads new flags only; Clear is always current-gas.**
  Reading that checkbox in `on_flag_clear` would turn one click into a
  five-gas deletion the button never advertised.
- **Flagging preserves the view** (`refresh(preserve_view=True)`) — it is a
  fine-grained edit made while zoomed in on the points being removed. The
  consequence is that the y-range does not reframe until Home is pressed,
  which on Ozone is the difference between 0–3663 and −57–929.
- **The Flag tool shares one selector with Stats** via `PlotPane.selector_mode`
  rather than adding a second list; `attach_stats_selectors` rebuilds
  `self.selectors` wholesale and disconnects the old handlers on every draw, so
  a parallel list would duplicate that and leak canvas connections. The two
  toolbar actions are mutually exclusive (both want left-drag), and flag mode
  is built `interactive=False` with `button=[1, 3]` — a flag is an action that
  fires and clears, not a standing selection.
- **`attach_stats_selectors` snapshots and restores `ax.dataLim`, then calls
  `autoscale_view()`.** A RectangleSelector adds its rectangle — and, when
  interactive, three handle Line2Ds — to the Axes *at the origin*, and those
  enlarge `dataLim` to include (0, 0). On a tracer-tracer scatter that is
  ruinous and obvious: N2O spans 304–341 ppb, so reaching back to zero
  squashes the whole correlation into a corner. Restoring `dataLim` alone is
  not enough — adding the artists has already triggered an autoscale off the
  polluted limits and nothing recomputes the view — hence the explicit
  `autoscale_view()`, which is safe because it only touches an axis whose
  autoscale is still on, leaving a preserved view or the ozone y-framing
  alone. The timeseries had the same corruption all along and survived only
  because its view limits are settled by other means; the Correlations tab is
  where it finally showed.
- **Undo is session-only.** A config that could undo its own contents would be
  a strange object.

#### Flagging from the Correlations tab

The reason the tool is there at all: an outlier that is obvious against
another tracer — one ozone point far off the O3/N2O line — can be
near-impossible to find in a timeseries. Stats stays hidden on that pane;
Flag does not.

- **A point belongs to two gases**, so the target is explicit: a
  "Flag applies to" combo offering the Y tracer, the X tracer, or both,
  defaulting to Y. `_populate_corr_flag_target` rebuilds it whenever the axis
  pickers change and **keeps the role, not the gas** — after a swap, "the Y
  tracer" is still what the user meant, and re-resolving to a gas name would
  silently retarget the tool. `on_corr_flag_clear` clears exactly the combo's
  scope, so the button undoes what the tool does.
- **The box is matched against the plotted values** — calibrated, or
  floor-filtered for Ozone — because on this figure those *are* the axes.
  That is the opposite of the timeseries rule (raw there, where two traces
  overlap and one has to be picked) and it costs nothing: a box is resolved to
  row numbers once, at the drag, and rows are what get stored. Nothing
  re-resolves later, so no flag can drift when the calibration changes.
- **Unflagging matches the same box**, again the opposite of the timeseries,
  and for the same underlying reason. There a flagged point is off-screen
  because the y-range is framed on the filtered data; here it is drawn in
  place as a struck-out marker, so there is always a box to draw around it.
- **`_corr_axis_flagged` recovers where a flagged point would have plotted.**
  Flagging blanks the row, so it drops out of the pairing and vanishes —
  leaving the outlier just removed invisible and unreachable. Nothing is
  re-derived to get it back: `calibrate_series` emits `cal_slope`/
  `cal_intercept` on *every* row, blanked ones included, precisely so a
  blanked value can be recomputed. A floor gas needs no undoing, but its
  below-floor faults stay out — those are a sensor fault, not a user choice.
- **Selection for flagging is restricted to `keep`; unflagging is not.**
  Flagged points are deliberately drawn *outside* the plotted record, so
  unflag has to reach them, while flag must only ever take points that are
  actually part of it.
- **"Hide flagged points" never replots.** It sets the marker artist's
  visibility and calls `draw_idle()`, nothing more: the user is typically
  zoomed in on the very points being hidden, and a redraw would throw that
  away. It is a *view* toggle, so it is session-only (like the calibrated
  overlay), does not dirty the config, and does not touch the flags —
  `refresh()` here would be wrong on every count.

  **Home is retargeted instead**, by `_NavToolbar.home` — a subclass whose
  Home consults an `(axes, xlim, ylim)` override before falling back to the
  stock behaviour. `PlotPane.set_home_view` only records that triple: it sets
  no limits and draws nothing.

  **A redraw that was asked to rescale does frame on the visible data**
  (2026-07-30), which is not a contradiction of the above but its complement:
  a hidden artist keeps its data limits, so the autoscale still reaches out to
  markers nobody can see — on Ozone, a 0–3663 ppb axis for a record that tops
  out at 1140 once four flagged spikes are hidden. `redraw_corr` gates this on
  `old_view is None`, i.e. on a *full rescale having been requested*, which is
  what a tracer change is (`on_corr_gas_changed` / `on_corr_swap_axes` pass
  `preserve_view=False`). Every other route preserves the view, and the Hide
  toggle still never redraws at all — so zooming in and then toggling flagged
  points on and off keeps them appearing and disappearing at the zoomed scale,
  which is the whole point of the toggle. The limits are applied **before**
  `reset_nav()`, so the visible frame becomes the nav base rather than sitting
  on top of one that still spans the hidden markers.

  **Overriding `home` rather than rewriting the nav stack is the point.** The
  first version cleared the stack (`toolbar.update()`) and pushed the wanted
  range as a new base, then restored the user's view on top. It worked, but
  it threw away every zoom and pan they had done — Back and Forward stopped
  working — and moved the axis limits twice per click for no visible reason.
  The history is the user's; only where Home lands is ours to redirect. The
  override is checked against `canvas.figure.axes` before use, because the
  panes rebuild their Figure from scratch and a stale override would name a
  destroyed Axes; `reset_nav()` clears it for the same reason.

  `_corr_home_limits` builds that range from `_corr_plotted` rather than from
  `ax.dataLim`, for two reasons: a hidden artist keeps its data limits, so
  `dataLim` still spans the markers; and the selectors have their own history
  of polluting it (see above). 5% margins, matching matplotlib's own autoscale
  — verified against it (N2O data 303.8–341.5 → 301.9–343.4).

  `redraw_corr` re-applies the retarget at the end, after `reset_nav()`, since
  the autoscale that captured is the one that still includes the hidden
  markers.

`save_config` now carries **two** flight-config blocks, so its
omitted-block-is-a-deletion trap is live again on that path: a save passing
`cal_selection=` but forgetting `flagged=` discards every flagged point.
`on_save_clicked` passes both.

### Export tab: two products from one set of "gas blocks"

Its own tab, like Cal Tanks and for the same reason: both files are
whole-flight deliverables covering **every** gas, while every control in the
left panel is per gas. The old per-gas "Export calibrated CSV…" button in the
Calibration group box is gone, and `export_calibrated_csv` with it — a per-gas
button was quietly the reason the old export could only ever describe one gas.

`_export_gas_blocks()` builds one dict per available gas and both writers
consume it, so the two files can never disagree about what a gas's calibrated
record is. Each block carries `raw`, `final` + `final_kind`
(`"calibrated"`, or `"filtered"` for a floor-only gas), optional
`sigma`/`slope`/`intercept`, a `masks` dict, and `reason` when a gas that
should have a calibration lacks one. Each gas is analysed through
`_analysis_for`/`_calibration_for` with **its own** saved settings, exactly as
the Correlations tab does, so a gas nobody selected this session still exports
correctly.

**Every Series in a block is on the RAW file's row numbering**, via
`_to_raw_rows`. `drop_presync_rows` both trims leading rows *and*
`reset_index(drop=True)`, so `self.df` row 0 is raw row `presync_dropped`;
an unshifted Series would line every gas up a few dozen rows early, silently
and plausibly. `self.raw_df` (datetime + the gas value columns only, kept at
load) is what those rows are put back onto.

- **The companion CSV is exactly as long as the raw file**, pre-sync rows
  included — blank in the derived columns and flagged in a `presync` column.
  That row-for-row promise is the whole feature: the two files are meant to be
  opened side by side, or pasted together, in Excel or Igor with no alignment
  step, and a quiet offset between two files that *look* alignable is a very
  effective way to corrupt an analysis.
- **The raw echo columns come from `self.raw_df`, not from `self.df`.** A
  pre-sync row has an unreliable *clock*, not an unreliable *reading*, so
  blanking its measurement would make a column that claims to echo the source
  disagree with it. The derived columns are blank there; the echo is not.
- **Masks are written as `Int64` 1/0, not bool.** pandas writes `True`/`False`,
  which Igor loads as a text wave and Excel as text. The nullable dtype leaves
  the pre-sync rows genuinely empty rather than claiming `0` for rows where
  the mask was never evaluated.
- **`comment_header` defaults False**, unlike the old exporter: neither Excel
  nor Igor skips a leading `#` block without being told to. With it off the
  provenance goes to a sidecar `<stem>_notes.txt`, so it is never simply lost.
- **A gas with no calibration still exports.** The CSV keeps its raw column
  and its masks (the answer to "which rows are air" does not need a
  calibration), and the notes give the reason. The ICARTT file drops the
  variable entirely rather than shipping uncalibrated counts under a
  calibrated-looking name.

#### ICARTT specifics (format index 1001)

The delivered record is the good ambient one: every row `calibrate_series`
blanked is written as `-99999`, which is exactly what that flag means, so no
mask columns are written or needed.

- **`icartt_time_base(datetimes, skip_leading)` — `skip_leading` is not
  optional in practice; pass `presync_dropped`.** Those rows carry the stale
  pre-sync clock, which runs *ahead* of the true time, so left in they set the
  running maximum hours into the future and every genuine row after the
  backward jump then fails the strictly-increasing test. On the Feb 2025 file
  that is 5013 rows rejected instead of 1435. `_export_time_base()` wraps it
  so no call site can forget; `export_icartt` takes the same argument.
- **ICARTT's independent variable must increase strictly, and this record does
  not.** The datetime column has duplicate timestamps (1435 in the Feb 2025
  file — the same duplication `interp_hold` works around), and a duplicate
  makes the file *invalid*, not merely untidy. Those rows are dropped and the
  count is reported in the success dialog and in a special comment, rather
  than being discovered later by an archive's validator.
- **`_one_line` vs `_field`.** Header lines 2–5 and the keyword values are
  free text where commas are legal and expected — `Dutton, Geoff` is the form
  the standard asks the PI name to take — so they only get newlines collapsed.
  Only the `name, unit, description` variable-definition lines sit in a
  comma-delimited position, and only those get commas replaced.
- **Line 1's count includes itself and the column-header line**, and the
  column-header line is the *last normal comment*, inside that section's
  count. Off-by-one here is the most common validator failure.
- **The data interval is declared as `0` (non-uniform) unless the diffs really
  are all equal** — a nominal 1 Hz the file does not keep to would be a claim
  about the data, not a description of it.
- **`UNCERTAINTY` is derived, not typed.** The median 1σ per gas comes from
  `calibration_uncertainty`; the metadata box is appended after it. It is the
  one required keyword this program already knows the answer to, and a
  hand-typed value goes stale the moment a drift model or cal window changes.
- `ULOD_FLAG`/`LLOD_FLAG` (-7777/-8888) are constants of the format and are
  hardcoded; only the LOD *values* are metadata (default `N/A`).

**The writer matches the sister UCATS instrument's delivered files**
(`~/Downloads/SABRE-UCATS-GC_WB57_20230303_R0.ict` was the reference), not
generic ICARTT practice, wherever the two differ — anyone processing a
campaign meets one convention across both instruments:

- **`ICARTT_MISSING = -99999`, not the more common `-9999`.** It is not fixed
  by the standard — it is declared per variable on header line 12 and readers
  honour what is declared — so matching the house style costs nothing.
- **Variable lines carry four fields**, `name, unit, standard_name,
  description`. `standard_name` comes from `GASES[gas]["standard_name"]`
  (`Gas_CO2_InSitu_S_DMF`). A gas without one gets the field **empty, not
  missing**: dropping it would shift the description into field 2 and make
  position mean different things on different lines. Every gas now has one, so
  nothing exercises that path — keep it anyway, since a new gas arrives without
  a name until someone looks it up.

  The names are not free text: they come from the **NASA ESDS Atmospheric
  Composition Variable Standard Names Convention** (ACVSNC), which ICARTT V2.0
  requires, as
  `MeasurementCategory_CoreName_AcquisitionMethod_DescriptiveAttributes`. For
  the `Gas` category the attributes are MeasurementSpecificity (`S` = single
  species, `M` = multiple, `NA`) then Reporting (`DMF` = molar fraction wrt dry
  air, `AMF` = wrt ambient, `DVMR`/`AVMR` the volumetric mixing ratios on those
  two bases, `None` = not stated).
  Two of the five gases are not `Gas_..._DMF`, both deliberately
  (2026-07-29, on the PI's instruction):

  - **H2O is `Met_H2OMF_InSitu_None`** — water vapour is not in the `Gas`
    category at all. It belongs to `Met`, whose format is
    `Met_CoreName_AcquisitionMethod_None` with no descriptive attributes.
    `H2OMF` is "mole fraction of water vapor with respect to ambient air";
    `H2OMRV` (volumetric mixing ratio) and `H2OMR` (mass mixing ratio to dry
    air) are the alternatives if the sensor is ever stated differently. This
    replaces an earlier deliberate blank.
  - **Ozone is `Gas_O3_InSitu_S_AVMR`**, not `..._DMF`: a volumetric mixing
    ratio against *ambient* air, water vapour included, because nothing dries
    the 2B monitor's sample stream. Settled with Eric on 2026-07-29 after the
    entry had been `DMF` and then briefly `None`; unlike the other gases it
    could not be copied from the sister GC file, which has no ozone variable.

  **`long_name` has to agree with the reporting attribute it sits beside.**
  Both of these carry an *ambient* basis, so O3 reads "Ozone volume mixing
  ratio in ambient air" and H2O "Water vapour mole fraction in ambient air" —
  not the bare "mole fraction" they had while their standard names still said
  `DMF`. The two end up on the same variable-definition line, and a
  description contradicting the controlled-vocabulary field beside it is worse
  than either alone.
- **Units gain the trailing `v`** (`ppmv`/`ppbv`) via `ICARTT_UNITS`, for the
  ICARTT header only. Plot labels keep `ppm`/`ppb`.
- **The 1σ variable is `<species>e_<suffix>`** (`CO2e_RASTA`), not
  `_unc`, and reuses its parent's standard name exactly as those files do —
  which the ACVSNC also requires: "when uncertainty is reported as a separate
  variable, the uncertainty variable should share the same standard name as
  the corresponding measurement variable."
- **Header line 9 has three fields** — `Time_Start, seconds, <description>` —
  not two with the description crammed into the unit.
- **Nothing about the analysis is written to the file.** The special-comments
  section carries the author's text and nothing else. It used to also append
  the source file name, every gas's masking/drift settings and the counts of
  rows dropped for the format's sake; that was removed on request
  (2026-07-29) as internal to the experimenters — the settings especially are
  a working record of how the analysis was tuned, not something a data user
  should be reading in a delivered file. `export_icartt` therefore takes
  neither `source_path` nor a settings note. The omission counts still reach
  the user through the returned summary dict and the post-export dialog,
  which is now their only route.
- **`_verbatim_lines` vs `_one_line`.** `special_comments` and
  `revision_history` are whole *sections*, written through with their line
  breaks (and the blank line inside the special comments) intact; the section
  line count is taken from the resulting list so the two cannot disagree.
  Every other field is a keyword *value* occupying exactly one line.
- **The revision history accumulates.** The delivered R0 file lists both
  `RA: Preliminary data` and `R0: Revised data`. It is a block the user
  maintains, not a note about the current revision, so it is free text; only
  when empty does the current revision get a placeholder line.
- **`icartt_filename` keeps hyphens and strips underscores.** `_` separates
  the file name's own fields, while a hyphenated data ID is normal
  (`SABRE-UCATS-GC_WB57_...`); the first version stripped all non-alphanumerics
  and silently mangled it to `SABREUCATSGC`.

#### ICARTT metadata lives in the APP-level config — the SHARED one

There are **two** app-level files beside the script, split by whether their
contents are worth sharing, and the split is what lets the first one be
tracked in git at all:

- **`ucatsb_gui_config.yaml` — ICARTT metadata, tracked.** The campaign's PI,
  affiliation, project and stipulations are the same on every machine, so a
  fresh clone should arrive with them filled in. Written by exactly one
  thing, `_save_shared_config`, called only from the **Save defaults** button.
- **`.ucatsb_gui_state.yaml` — the recent-files list, gitignored.** Absolute
  paths from one machine, useless to anyone else. Hidden because the app
  maintains it and there is nothing in it to hand-edit. Written by
  `_save_state`, called on **every dataset load** — which is precisely why it
  cannot share a file with the tracked one: it would leave that file
  permanently modified in the working tree, and commit one machine's paths.

`load_recent_files(path, legacy_path=)` reads the shared config as a fallback
so an upgrading user does not silently lose the list; the stale
`recent_files:` block there is never written back, and the first
`Save defaults` drops it.

`icartt:` in `ucatsb_gui_config.yaml`, not in a `<dataset>_conf.yaml`
(explicitly requested): PI, affiliation, project and stipulations are
properties of the campaign, so per-flight storage would mean retyping them for
every file. Its own **Save defaults** button on the tab, its own
`_saved_icartt_meta` snapshot, and its own `_confirm_discard_icartt` on close
— kept separate from `_is_dirty()`/`_confirm_discard` because the two Save
buttons write *different files*, and one prompt offering to "save" would have
to pick one and silently not write the other. Dirty state is a comparison, not
a flag, for the same reason as `_is_dirty()`.

`save_config` writes a **fresh document**, so an omitted block is a deletion.
That used to be a live hazard here — while one file held both blocks, every
write had to remember both or delete the other. One block per file retired it;
don't reintroduce a writer that carries two. Adding a metadata field means
adding it to `DEFAULT_ICARTT_META` *and* `ICARTT_FIELDS` — `_icartt_meta_from_controls`
rebuilds the dict wholesale from the widgets, the same trap
`_controls_to_settings()` has.

The Export tab joins the `_dirty`/`_preserve` dispatch in `_draw_current_tab`
even though it draws no figure: it reads every gas's analysis to build its
readiness summary, and going through the dispatch stops a spin-box nudge from
recomputing five gases while the tab merely happens to be open.

### Cal Tanks tab

Its own tab (`_build_cal_tanks_pane`), not another control-panel group box:
the pairing is one per *flight*, while every control in that panel is per gas.
Two combos over the whole roster (`load_cal_roster`), not just the plumbed
pair — picking a tank outside the current `cals:` block is the entire point.

- `self.cal_bottles` is **derived state**: `_rebuild_cal_bottles` recomputes
  it from `cal_selection` via `select_cal_bottles`, which is the one place
  implementing "matching may only see the plumbed tanks" (see the section on
  `cals.yaml` above). Never assign `cal_bottles` directly.
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
