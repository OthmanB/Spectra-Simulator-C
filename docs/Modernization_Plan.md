# Spectra-Simulator-C Modernization Plan

This document is a **multi-iteration plan** for improving correctness, maintainability, portability, and documentation for Spectra-Simulator-C.

Scope note: this plan focuses on **non-obsolete** models, defined as the models with configuration examples in `Configurations/examples_cfg/` **excluding** `Configurations/examples_cfg/obselete/`.


## Goals

- Make model selection/configuration **self-validating**, with clear error messages.
- Fix known correctness bugs in non-obsolete models.
- Improve portability (paths, filesystem operations, reproducible RNG).
- Replace ad-hoc `std::cout` logging with a structured **log-level system** (`debug/info/warn/error`).
- Reduce “wire-in-multiple-places” friction by introducing a **single model registry**.
- Produce **model manuals** (1 Markdown file per non-obsolete model) that are configuration-driven and code-driven, with equations.


## Non-Goals (For Now)

- Large physics changes to the models.
- Changing the leakage/apodization handling in `artificial_spectrum.cpp` (deferred; see “Deferred Items”).
- Unifying `nu_spread` semantics across unrelated models (explicitly NOT requested).


## Non-Obsolete Model Inventory (Current)

From `Configurations/examples_cfg/`:

- `generate_cfg_from_synthese_file_Wscaled_Alm` (MS)
  - Example: `Configurations/examples_cfg/main.cfg.Alm`
  - Generator: `models_database.cpp:generate_cfg_from_synthese_file_Wscaled_Alm()`
  - Builder: `artificial_spectrum.cpp:artificial_spectrum_a1Alma3()`

- `generate_cfg_from_synthese_file_Wscaled_aj` (MS)
  - Example: `Configurations/examples_cfg/main.cfg.aj`
  - Generator: `models_database.cpp:generate_cfg_from_synthese_file_Wscaled_aj()`
  - Builder: `artificial_spectrum.cpp:artificial_spectrum_aj()`

- `generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014` (MS)
  - Example: `Configurations/examples_cfg/main.cfg.aj_GRANscaledKallinger`
  - Generator: `models_database.cpp:generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014()`
  - Builder: `artificial_spectrum.cpp:artificial_spectrum_aj(..., noise_modelname="harvey_like")`
  - Noise cfg: `Configurations/noise_Kallinger2014.cfg` (merged into main cfg)

- `asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014` (evolved / mixed modes)
  - Example: `Configurations/examples_cfg/main.cfg.freeDP_curvepmodes.v3_GRANscaled`
  - Generator: `models_database.cpp:asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014()`
  - Builder: `artificial_spectrum.cpp:artificial_spectrum_aj(..., noise_modelname="harvey_like")`
  - Optional override file format: `Configurations/MixedModes_models/star_params.theoretical.EXAMPLE`


## Phase 0: Baseline + Guard Rails

Deliverables:

- A repeatable “smoke run” command per non-obsolete model.
- A minimal test harness (even if only a script) to run these smokes and check that expected output files are created.

Tasks:

1) Add a small doc section or script (e.g. `scripts/smoke_test.sh`) that:
   - builds via CMake
   - runs `specsim` on each file in `Configurations/examples_cfg/`
   - writes outputs to a temporary output directory
   - verifies:
     - `Combinations.txt` created
     - at least 1 spectrum file created
     - `Spectra_info/<id>.in` created

2) Fix/guard critical prerequisites revealed by the smoke tests:
   - Ensure `Configurations/tmp/` exists at runtime (some generators write `Configurations/tmp/*.cfg` and otherwise fail).
   - Ensure `--out_dir` produces `Combinations.txt` under the output directory (avoid missing path separator).
   - Build system should not hard-fail when `gnuplot` is missing (plotting is runtime-only; smoke tests run with `doplots=0`).

3) Ensure smoke runs do not depend on absolute paths embedded in example configs.
   - The smoke harness may patch `extra_params` to point to repo-local `Configurations/infiles/*.in`.

Tests (implemented for Phase 0):

- Unit:
  - `tests/test_smoke_patch_cfg.py` validates `scripts/smoke_test.py:patch_main_cfg()` rewriting invariants (forest type, steps length, template pinning, observation line, inline comments).
- Integration / smoke (black-box binary):
  - `tests/test_specsim_phase0_integration.py` builds `specsim` and runs representative configs in a sandbox to ensure:
    - `Configurations/tmp/` is created automatically
    - `Combinations.txt` is written under `--out_dir`
    - output directory tree exists and basic artifacts are produced
- Property (single-version):
  - `tests/test_specsim_phase0_integration.py` checks frequency-grid invariants (monotonicity, endpoints, constant spacing, expected row count).

Tests (deferred to later phases):

- Cross-version numeric regression of full spectra values is deferred until Phase 5 introduces deterministic seeding (`--seed`).


## Phase 1: Configuration Validation (Correctness + UX)

User decision: OK to implement.

Problem summary:

- `io_star_params.cpp:read_main_cfg()` has an incorrect vector-length check (`&&` where `||` is required).
- Model parameter validation is currently partial; many failures become runtime errors later (or worse: wrong parameter mapping).

Deliverables:

- Each model run fails fast with:
  - “Expected parameters (ordered list)”
  - “Provided labels (from cfg)”
  - missing/extra/duplicate labels
  - wrong vector sizes

Tasks:

1) Fix length check in `io_star_params.cpp:read_main_cfg()`.

2) Add a single validation function, called immediately after `cfg=read_main_cfg(cfg_file)` and after optional noise aggregation, e.g.:
   - `validate_cfg_vectors(cfg)`
   - `validate_cfg_for_model(cfg, model_spec)`

3) Validate that:
   - `cfg.labels`, `cfg.val_min`, `cfg.val_max`, `cfg.step`, `cfg.distrib` are all same length
   - `cfg.labels` contains no duplicates
   - every expected `param_names[]` appears exactly once in `cfg.labels`
   - no unexpected labels remain (unless explicitly allowed by model)

4) Improve error reporting in `ioproc.cpp:order_input_params()`:
   - show the name that failed to match
   - show model name
   - show full `param_names` vs `cfg.labels`

Tests (implemented):

- Unit / integration (black-box):
  - `tests/test_specsim_phase0_integration.py` covers:
    - malformed cfg vectors (mismatched sizes) fail fast with clear message
    - duplicate labels detected
    - missing/extra labels reported with expected/provided lists
    - non-obsolete example cfgs still run after validation is enabled
- Property (cross-version, modulo error fixes):
  - valid configs accepted before must remain accepted
  - invalid configs must fail earlier (but with improved diagnostics)


## Phase 2: Fix Known Non-Obsolete Bug (Width Jitter in `generate_cfg_from_synthese_file_Wscaled_aj`)

User decision: OK to fix with a proper loop.

Problem summary:

- In `models_database.cpp:generate_cfg_from_synthese_file_Wscaled_aj()` the `gamma_spread` logic uses `gamma_star[i]` with `i` not set (undefined behavior).

Deliverables:

- Width jitter applies to all modes deterministically and safely.
- A short note in the model manual documenting how `Gamma_spread` works.

Tasks:

1) Implement a loop across all modes when applying `gamma_spread`.
2) Ensure jitter cannot produce negative widths.
3) Add a small debug/info log line summarizing applied spread (controlled by log level).

Tests (implemented):

- Integration / property (black-box):
  - `tests/test_specsim_phase0_integration.py` includes `test_gamma_spread_applies_per_mode`:
    - `Gamma_spread=0` keeps widths constant (with a constant-width reference infile)
    - `Gamma_spread>0` produces per-mode width diversity and widths remain positive


## Phase 3: Standardize `delta0l_percent` Convention (Kallinger2014 Preferred)

User decision:

- Standardize `delta0l_percent` to the Kallinger2014 convention.
- Do NOT standardize `nu_spread` because it is intentionally different across unrelated models.

Proposed standard:

- `delta0l_percent` is provided as a **positive magnitude in percent** in main config.
- Internally, the generator uses **negative small separations** (consistent with p-mode asymptotics).
  - In practice: `delta0l_percent_internal = -1 * delta0l_percent_input` before computing `d0l = delta0l_percent_internal * Dnu / 100`.

Deliverables:

- All non-obsolete models that use `delta0l_percent` treat the input the same way.
- Example config files in `Configurations/examples_cfg/` remain readable and are updated/documented as needed.

Tasks:

1) Identify all non-obsolete generators that use `delta0l_percent` and unify the sign handling.
2) Update per-model docs to explicitly state the sign convention.
3) Update any example cfg comments that contradict the convention.
4) Add a validation warning if user provides a negative `delta0l_percent` (likely unintended).

Tests (implemented):

- Integration / property (black-box):
  - `tests/test_specsim_phase0_integration.py` includes `test_negative_delta0l_percent_warns`
    - negative input emits a warning and the run completes

Compatibility note:

- This change can alter the generated frequencies for existing configs that previously assumed a different sign.
- Plan: treat this as a deliberate behavior change, document it in `changelog.md`, and verify with smoke tests.


## Phase 4: Template Selection Exclusion Mechanism (Bad Template Skipping)

User decision: OK; add exclusion mechanism with warning.

Problem summary:

- Some templates in `Configurations/templates/` are not parseable by `external/ARMM/bump_DP.cpp:read_templatefile()` which assumes 3 numeric columns in the data table.
- Example known-bad file: `Configurations/templates/KIC12069424.template` has a row with only 2 columns.

Deliverables:

- When selecting templates (especially when `template_files = all`), invalid templates are:
  - skipped
  - logged as `warning`
  - do not abort the full run
- If all candidate templates are invalid, fail with a clear error.

Tasks:

1) Implement `validate_template_file(path)` that checks:
   - required keywords exist (`numax_ref`, `Dnu_ref`, `epsilon_ref`)
   - data table has >=1 row
   - each data row has exactly 3 parseable numeric columns

2) Integrate validation at the template selection step:
   - `iterative_artificial_spectrum.cpp:generate_random()` (and any other template selection location)

3) Add a small helper CLI / developer tool (optional but recommended):
   - `specsim --validate-templates <dir>` or a standalone `scripts/validate_templates.py`

4) Optional follow-up (separate iteration): clean/fix the shipped template files themselves.

Tests (plan):

- Unit:
  - `validate_template_file()` accepts valid templates and rejects malformed ones (wrong columns, missing keywords)
- Integration:
  - evolved model with `template_files=all`:
    - skips invalid templates with warnings
    - succeeds if at least one valid template exists
    - fails clearly if none are valid


## Phase 5: Portability + Robustness Improvements

User decision: Agree.

Deliverables:

- No shelling-out to `rm`, `ls`, `grep` for core operation paths.
- No memory leaks for `getcwd(NULL,0)` usage.
- Output and config handling works on both macOS and Linux with relative paths.
- Reproducibility option via explicit seed.

Tasks:

1) Replace shell calls:
   - `system("rm ...")` -> `std::filesystem::remove` with error handling
   - `list_dir()` implementation that uses `ls|grep` -> `std::filesystem::directory_iterator`

2) Replace `getcwd(NULL,0)` usage with `std::filesystem::current_path()`.

3) Add CLI flag `--seed`:
   - if provided, use deterministic seeding in random sampling and model generators
   - otherwise keep current behavior

4) Remove (or minimize) mixed RNG APIs:
   - stop calling both `srand/rand` and `<random>` in the same model

5) Fix example config portability:
   - replace absolute paths in `Configurations/examples_cfg/*.cfg` with repo-relative paths
   - document how paths are resolved (`--main_dir`, absolute vs relative)

Tests (plan):

- Unit:
  - filesystem helper functions (remove/list) behave correctly on empty/invalid paths
- Property:
  - `--out_dir` works with and without trailing slash
  - `--seed` provides deterministic outputs (byte-for-byte for deterministic artifacts; numeric tolerance for floats if needed)
- Integration:
  - run two identical jobs with same `--seed` and confirm identical outputs


## Phase 6: Structured Logging Upgrade (Info/Debug/Warning/Error)

User requirement: add a proper logging system; stop scattering raw `std::cout` everywhere.

Recommended approach:

- Introduce a thin project wrapper API (`logging.h/.cpp`) with macros/functions:
  - `LOG_DEBUG(...)`, `LOG_INFO(...)`, `LOG_WARN(...)`, `LOG_ERROR(...)`
- Back the wrapper with a dedicated logging library.

Library choice notes:

- Recommended: `spdlog` (widely used, easy levels/sinks/formatting).
  - Due to `cmake_minimum_required(VERSION 3.5)`, prefer a vendored header-only integration under `external/`.
- Alternative: `Boost.Log` (already using Boost, but adds link complexity).

Deliverables:

- CLI option `--log-level` with values like `debug|info|warn|error`.
- Default log level `info`.
- Conversion strategy that does not require a single giant refactor.

Tasks:

1) Add wrapper and minimal initialization (set level, set format, optional file sink).
2) Convert critical paths first:
   - configuration parsing/validation
   - template selection
   - model registry selection
   - fatal errors -> `LOG_ERROR` + controlled exit
3) Convert remaining files gradually.
4) Keep user-facing progress messages at `info`.
5) Use `debug` for verbose internal dumps currently printed unconditionally.

Tests (plan):

- Unit:
  - `--log-level` parsing (invalid level fails with helpful message)
  - log filtering: `debug` messages suppressed at `info`
- Integration:
  - a known error path emits `error`-level log and non-zero exit code


## Phase 7: Model Switching Refactor (Single Registry)

User decision: Yes.

Problem summary:

- Adding a model requires manual wiring in multiple places in `iterative_artificial_spectrum.cpp`:
  - parameter list (`param_names`)
  - `call_model_random()`
  - `call_model_grid()`

Deliverables:

- A single model registry that defines:
  - model name
  - supported forest types (`random`, `grid`)
  - ordered `param_names`
  - whether it needs `extra_params` (.in reference)
  - whether it needs templates
  - whether it needs `noise_Kallinger2014.cfg`
  - generator function pointer
  - builder function pointer
  - any special glue (e.g., appending `Tobs/Cadence` to `input_params`)

Tasks:

1) Create `models_registry.h/.cpp`.
2) Implement lookup by `cfg.model_name`.
3) Move all `param_names` definitions into registry entries.
4) Implement per-model validation based on registry metadata.
5) Replace the big `if(cfg.model_name == ...)` chain with a registry call.
6) Replace `call_model_random/grid` logic with registry dispatch.
7) Add CLI helper:
   - `--list-models`
   - optionally `--describe-model <name>`

Tests (plan):

- Unit:
  - registry lookup finds all non-obsolete models
  - registry param list equals historical `param_names` ordering (for each non-obsolete model)
- Integration:
  - `--list-models` prints expected names
  - each non-obsolete example cfg runs through registry dispatch

Incremental strategy:

- Start by registering only the non-obsolete models.
- Keep legacy chains temporarily behind a compile-time flag until parity is confirmed.


## Phase 8: Documentation Restructure (Per-Model Manuals)

User decision: OK.

Deliverables:

- `docs/models/README.md` (model chooser + MS vs evolved)
- One manual per non-obsolete model:
  - `docs/models/generate_cfg_from_synthese_file_Wscaled_Alm.md`
  - `docs/models/generate_cfg_from_synthese_file_Wscaled_aj.md`
  - `docs/models/generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014.md`
  - `docs/models/asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014.md`

Manual structure (template):

- What the model is for
- Required inputs (main cfg fields, extra files, templates)
- Parameter table in canonical order (matching registry)
- Physics assumptions + key equations (Markdown math)
- How to run (CLI example)
- Outputs produced
- Common pitfalls and debugging tips

Update entry points:

- Add links in `docs/User_Guide.md` to `docs/models/README.md` and to each manual.

Tests (plan):

- Docs lint (optional):
  - link checker for `docs/models/*.md`
  - ensure every non-obsolete model has a manual and is linked from `docs/models/README.md`


## Deferred Items (Explicitly Not Blockers)

These are tracked but not acted on until their original rationale is recovered.

1) Leakage/apodization filter application in `artificial_spectrum.cpp` (previously flagged as “coded but not applied”).
   - Action now: **document as deferred**; do not change behavior.
   - Later investigation: locate any analysis or downstream tool that assumes current behavior.

2) Noise-model mismatch between random vs grid paths.
   - Action now: **document as deferred**; do not change behavior.
   - Later investigation: identify intended differences and whether they are still required.


## Cross-Version Regression Harness (Deferred Until Deterministic Seed)

Objective: ensure invariant behavior between the baseline and improved versions, modulo intentional error fixes.

- Requires Phase 5 (`--seed`) to make stochastic components reproducible.
- Implement a harness that runs two binaries (baseline and candidate) in isolated sandboxes and compares:
  - deterministic artifacts (cfg dumps, `Combinations.txt`, `Spectra_info/*.in`)
  - numeric spectra (`Spectra_ascii/*.data`) either byte-for-byte or within tolerance
  - allowlisted diffs (e.g., improved error messages, changed output paths for corrected bugs)


## Suggested Milestones / PR Breakdown

1) PR A: Phase 1 (cfg validation) + Phase 2 (width jitter bug)
2) PR B: Phase 3 (`delta0l_percent` standard) + update examples/comments
3) PR C: Phase 4 (template exclusion mechanism)
4) PR D: Phase 5 portability (filesystem + seed)
5) PR E: Phase 6 logging (introduce wrapper + start migrations)
6) PR F: Phase 7 model registry (non-obsolete models first)
7) PR G: Phase 8 documentation (per-model manuals)


## Acceptance Criteria (Definition of Done)

- All configs in `Configurations/examples_cfg/` run successfully with smoke tests.
- Non-obsolete models have dedicated manuals and are discoverable from `docs/User_Guide.md`.
- Adding a model requires editing only:
  - the generator implementation, and
  - one registry entry.
- Logging supports at least `debug/info/warn/error` and defaults to `info`.
- Template selection does not abort the run when encountering a malformed template; it warns and skips.
