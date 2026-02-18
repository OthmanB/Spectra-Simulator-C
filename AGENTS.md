# Agent Guide (Spectra-Simulator-C)

This repository is a C++17 CLI tool (`specsim`) that generates synthetic asteroseismic power spectra from configuration files.

Agent-specific rules:

- Copilot instructions exist at `.github/copilot-instructions.md` (follow them when working in this repo).
- No Cursor rules were found (`.cursor/rules/` or `.cursorrules`).

Repo reality note: this project currently uses custom text `.cfg` files (not YAML). Treat “strict config + early validation” guidance as applicable, but do not introduce YAML or reorganize docs paths unless explicitly requested.


## Quick Orientation

- Entry point / CLI: `iterative_artificial_spectrum.cpp` (builds `specsim`)
- Main cfg parsing + I/O: `io_star_params.cpp`, `ioproc.cpp`
- Model generators (write `Configurations/tmp/modes_tmp.cfg` + `Configurations/tmp/noise_tmp.cfg`): `models_database.cpp`, `models_database_grid.cpp`
- Spectrum builders (read tmp cfgs, write outputs): `artificial_spectrum.cpp`
- Mode profiles: `build_lorentzian.cpp`
- Noise models: `noise_models.cpp`
- Evolved stars / mixed modes: `external/ARMM/*`
- Activity / Alm model glue: `external/Alm/Alm_cpp/activity.cpp`
- Modernization roadmap: `docs/Modernization_Plan.md`

Non-obsolete example configs live in `Configurations/examples_cfg/` (exclude `Configurations/examples_cfg/obselete/`).


## Build Commands

Dependencies (typical):

- CMake
- C++17 compiler (macOS Clang or Linux GCC)
- Eigen
- Boost (headers + libs: filesystem, program_options, iostreams; boost_system optional)
- gnuplot (optional at runtime; required only when `doplots=1`)
- OpenMP (optional; controlled by `-DWITH_OPENMP=ON/OFF`)

Build (recommended out-of-tree):

```bash
cmake -S . -B build
cmake --build build -j
```

Linux note: the current `CMakeLists.txt` expects Eigen include path via `EIGEN3_INCLUDE_DIR` env var.
Example: `export EIGEN3_INCLUDE_DIR=/usr/include/eigen3`.


## Run Commands

Show help:

```bash
./build/specsim --help
```

Run an example config (main-sequence aj model) from repo root:

```bash
./build/specsim \
  --main_file main.cfg.aj \
  --main_dir Configurations/examples_cfg \
  --out_dir /tmp/specsim-out \
  --force-create-output-dir 1
```

Notes:

- `--main_file` can be absolute; otherwise it is joined with `--main_dir`.
- `--noise_file` defaults to `noise_Kallinger2014.cfg` and is resolved from `--main_dir` or `Configurations/`.
- Outputs are written under `--out_dir`:
  - `Combinations.txt`
  - `Spectra_ascii/`, `Spectra_info/`, `Spectra_modelfile/`, `Spectra_plot/`


## Test Commands

Current tests are Python `unittest`-based and live under `tests/`.

Run all tests:

```bash
python3 -m unittest discover -s tests
```

Run a single test file:

```bash
python3 -m unittest tests.test_specsim_phase0_integration
```

Run a single test case:

```bash
python3 -m unittest \
  tests.test_specsim_phase0_integration.TestSpecsimPhase0Integration.test_frequency_grid_invariants
```

Smoke test (build + run all non-obsolete example cfgs):

```bash
bash scripts/smoke_test.sh
```

The integration tests build `specsim` into `build-phase0-tests/` (ignored by git via `.gitignore`).


## Lint / Formatting

No repository-wide formatter/linter is enforced (no `.clang-format`, `.clang-tidy`, etc.).

Guidance:

- Avoid drive-by reformatting of large files.
- If you format, do it locally to a small touched region and keep diffs reviewable.
- Prefer adding focused tests over style-only changes.


## Code Style Guidelines (Project-Conformant)

Includes:

- Group includes in this order when touching a file:
  1) standard library (`<...>`)
  2) third-party (`<Eigen/...>`, `<boost/...>`)
  3) local headers (`"..."`)
- Avoid introducing `using namespace std;` in headers.

Types:

- Use explicit Eigen types (`Eigen::VectorXd`, `Eigen::MatrixXd`, `Eigen::VectorXi`) or file-local `using Eigen::...` if already present.
- Prefer `size_t` for indexing STL containers; cast carefully when interacting with Eigen APIs that expect `int`.

Naming:

- Follow existing local conventions of the file you are editing.
- Prefer descriptive names for parameters that map to config labels (keep label spelling stable).
- Keep constants `const` and initialize at declaration.

Error handling:

- For CLI/runtime failures: print a clear message to `std::cerr` with context (cfg path, model name, label) and exit with `EXIT_FAILURE`.
- Validate early (on read/parse) rather than failing deep in model generation.
- Avoid `system("rm ...")` / `system("ls|grep")` in new code; use filesystem APIs instead.

Randomness / determinism:

- Many algorithms are stochastic and currently seed from time.
- Tests should target deterministic invariants (file presence, vector sizes, frequency-grid properties) until a dedicated `--seed` exists.

Performance:

- Keep hot loops simple; avoid unnecessary allocations inside per-mode/per-frequency loops.
- Prefer Eigen vectorized operations when it does not reduce clarity.


## Configuration / Model Development Notes

- “Models” are selected by `cfg.model_name` and currently require wiring in multiple places in `iterative_artificial_spectrum.cpp`.
- When updating configuration semantics, update:
  - example cfgs under `Configurations/examples_cfg/`
  - the relevant generator in `models_database.cpp`
  - docs and tests that rely on these cfgs
- Template files in `Configurations/templates/` are consumed by ARMM; the data table must be 3 numeric columns.

## Additional rules to follow
Good practices and required coding discipline is further detailed in .github/copilot-instructions.md.
Refer to it before any implementation.
