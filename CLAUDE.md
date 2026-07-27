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
`shade_intervals`, `cal_mean_points`, `load_cal_bottles`, `most_common_serial`,
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
   that didn't affect the computed means, which was silently wrong.
4. `cal_mean_points` averages each cal window (offset in seconds relative to
   the interval's last timestamp, independently configurable per bottle) and
   identifies which physical bottle was flowing.

Both views read these through `UcatsbGui._get_analysis()`, a cache whose sole
invalidation site is `refresh()` — the entry point every state change calls.
It invalidates unconditionally rather than comparing a composite key of
(file, gas, warm-up, tolerance, cal windows, drift model): such a key is easy
to get subtly wrong and then serves a stale plot. The recompute is
milliseconds against ~100 ms of rendering.

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
current run. `load_cal_bottles` only returns the two tanks named in `cals:`
(`{serial: data[serial] for serial in set(data["cals"].values())}`) — this
is deliberate, not an oversight: matching should only ever consider the
tanks physically in use this flight, since an unrelated roster tank could
otherwise coincidentally match a measured value more closely and produce a
wrong identification. When a tank is swapped between flights, `cals:` is
what needs editing — the roster entries themselves don't move.

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
between flights (as long as `cals:` is updated to name the new pair) — do
not "fix" it to use the `cals.yaml` key order directly.

`cals.yaml` has previously contained literal tab characters as
`key:\tvalue` separators, which are invalid YAML syntax. `load_cal_bottles`
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
result of *this repo's* `calibrate_series`, keeping the raw trace visible
underneath at `alpha=0.35`. It is session-only and defaults off, precisely so
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

### Tabs (`PlotPane`, `refresh()` dirty dispatch)

The two views are `PlotPane` instances (Figure + canvas + toolbar) in a
`QTabWidget`. `self.figure`/`self.canvas`/`self.toolbar` stay **bound to the
timeseries pane** rather than being renamed, so `redraw()`'s body needs no
changes; don't rename them while also changing behavior, or a rendering
regression becomes indistinguishable from a refactor slip.

The controls panel deliberately stays **outside** the tabs — every control
affects both views, so moving it inside would mean duplicating the gas
selector or making one tab depend on state invisible from the other.

`refresh()` redraws only the visible pane and marks the other dirty, so a
spinbox drag doesn't render a pane nobody is looking at. The `_preserve`
latch matters: a requested full rescale must survive until it is actually
honoured, or changing gas and then nudging a spinbox would leave the
never-drawn pane stuck at a stale scale when first opened. The cal tab's own
"same content" key is `(gas, drift_model, smooth_events)`.

### Config persistence

`ucatsb_gui_config.yaml` (loaded/saved by `load_config`/`save_config`) holds
one settings block per gas (`warmup_min`, `pressure_tol_mbar`,
`cal1_window_s`, `cal2_window_s`, `drift_model`, `drift_smooth_events`),
auto-saved on every control change via `on_control_changed`.
`_initializing`/`_loading` flags exist specifically to suppress redraw/save
during programmatic widget setup (e.g. `setChecked` on a freshly-constructed
radio button fires its signal immediately, before sibling widgets it might
depend on exist yet) — keep that guard pattern when adding new controls.

**Adding a persisted setting requires four edits, not one.**
`on_control_changed` assigns `self.config[gas] = self._controls_to_settings()`
— a *fresh* dict — so a key missing from `_controls_to_settings()` is
silently dropped from the file on the next control change, even though
`load_config`'s `.update()` appeared to preserve it. Touch all of:
`DEFAULT_GAS_SETTINGS`, `_controls_to_settings()`,
`_apply_settings_to_controls()`, and the `setEnabled(has_masking)` list in
`_select_gas()`.

**Verification scripts must pass `config_path=` to a scratch file.**
`UcatsbGui` writes the real `ucatsb_gui_config.yaml` on any programmatic
`setValue`, so an offscreen test that drives controls will otherwise silently
overwrite the user's saved settings — and then the numbers in the next run
won't match, which is confusing to debug.
