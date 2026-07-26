# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

Post-flight analysis tools for UCATS-B CSV logs. This is the analysis-side
companion to `~/code/ucats-b`, which runs on the aircraft and produces
`ucatsb-YYYYMMDDHH.csv` files (documented in that repo's `config-plot.yaml`
comment block). This repo does not acquire data — it only reads CSVs already
on disk, from any flight, passed as a CLI argument (never hardcode a specific
flight's filename or date).

## Requirements

`pandas>=2.2`, `matplotlib>=3.9`, `PyQt5>=5.15`, `PyYAML>=6.0`, Python 3.9+.
See `requirements.txt`.

## Running

- `python3 ucatsb_gui.py <csv_file>` — interactive PyQt5 viewer
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
`mean_std_label`, `CALS_YAML_PATH`. Any change to masking/cal-detection
behavior belongs in that shared module so both stay in sync.

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

### Cal bottle identity is matched by concentration, not trusted from config

`j_sol_cals` is a raw digital solenoid state (0 or 1) — it is **not** a
reliable bottle identifier. `cals.yaml` (a local copy of
`~/code/ucats-b/cals.yaml` — resync by hand if the acquisition repo's
bottles/serials change) maps `cal0`/`cal1` to serials (e.g.
`cal0: CC302489`), but that key order does not necessarily match which
serial is actually plumbed to which digital state on a given flight
(confirmed empirically: on the reference flight, `j_sol_cals==0` measured
~217 ppm CO2, matching `CB09960`'s nominal 206.51 ppm, not `CC302489`'s
418.95 ppm — the reverse of the naive `cal0`→digital-0 reading).
`match_cal_serial` instead picks whichever serial's nominal concentration
(for the active gas) is closest to the measured window mean. This is
self-correcting if bottles are swapped between flights — do not "fix" it to
use the `cals.yaml` key order directly.

`cals.yaml` also contains literal tab characters as `key:\tvalue` separators,
which are invalid YAML syntax. `load_cal_bottles` blanks tabs to spaces
before parsing rather than erroring — don't remove that workaround without
checking the file's current formatting (it'll be reintroduced any time
`cals.yaml` is re-copied from the acquisition repo).

### d1 vs d2 detector routing

The Aeris `d1` detector measures CO2 and N2O; `d2` measures CO and a second,
redundant N2O channel (not CH4 — the CH4 Aeris wasn't installed yet for the
flights this tool was built against, per the CSV column set actually
present). The main plotted trace is always `d1_*` for both gases
(`d1_CO2_ppm`, `d1_N2O_ppb`) — `d2` is only surfaced via the optional "Trace
Above" panel's Detector Pressure/T_gas options, and only for `d2_P_mbars`/
`d2_T_gas` when N2O is the active gas (confirmed with the user; do not swap
this without re-confirming, `d2` has no CO2 channel at all so there is no
symmetric "d2 for CO2" case).

Data is plotted **uncalibrated** (`d1_CO2_ppm`, `d1_N2O_ppb`), not the
`*c_ppm`/`*c_ppb` calibrated columns — this was a deliberate switch; don't
revert to the calibrated columns without being asked.

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
  *same* trace as before (tracked via `_last_aux_selection`) — a different
  trace has a meaningless old y-range, so it re-autoscales.
- `self.toolbar.update()` + `self.toolbar.push_current()` reset the
  NavigationToolbar's Home target to the newly-built full-scale view (its
  nav stack otherwise still references the just-destroyed Axes objects).
  Only `on_gas_changed` skips `preserve_view` — switching species changes
  the y-range meaning entirely, so a full rescale there is correct.

### Config persistence

`ucatsb_gui_config.yaml` (loaded/saved by `load_config`/`save_config`) holds
one settings block per gas (`warmup_min`, `pressure_tol_mbar`,
`cal1_window_s`, `cal2_window_s`), auto-saved on every control change via
`on_control_changed`. `_initializing`/`_loading` flags exist specifically to
suppress redraw/save during programmatic widget setup (e.g. `setChecked` on
a freshly-constructed radio button fires its signal immediately, before
sibling widgets it might depend on exist yet) — keep that guard pattern when
adding new controls.
