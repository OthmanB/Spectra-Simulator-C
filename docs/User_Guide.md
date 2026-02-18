# Spectra-Simulator User Guide

This guide is meant to complement the Doxygen API docs with a **practical, end-to-end workflow**: install dependencies, build, run a first simulation, and understand configuration files.

---

## 1) What Spectra-Simulator does

`specsim` generates synthetic stellar power spectra for asteroseismology, including mode profiles and noise realizations.

At a high level, each run:

1. Reads a main configuration (`main.cfg`-style file).
2. Optionally reads a separate noise configuration (`noise_Kallinger2014.cfg`-style file).
3. Samples parameters from either:
   - a random sampling mode (`forest_type = random`), or
   - a grid (`forest_type = grid`).
4. Writes spectra and metadata to output subfolders.

---

## 2) Dependencies

The project requires:

- C++17 compiler (GCC/Clang)
- CMake (recommended)
- Boost (`system`, `filesystem`, `iostreams`, `program_options`)
- Eigen3
- gnuplot (required by CMake configuration)
- OpenMP (optional, enabled by default)

> Linux note: the CMake file expects `EIGEN3_INCLUDE_DIR` to be available in the environment on Linux.

---

## 3) Build instructions

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

This creates the executable `build/specsim`.

If you prefer the legacy workflow, you can still copy the binary to the repository root; however, running directly from `build/specsim` is usually cleaner.

---

## 4) Command-line usage

Run help:

```bash
./build/specsim --help
```

Main options:

- `-h, --help`: print help
- `-v, --version`: print version/build info
- `-f, --main_file`: main config filename (default: `main.cfg`)
- `-n, --noise_file`: noise config filename (default: `noise_Kallinger2014.cfg`)
- `-g, --main_dir`: directory containing config files (default: `Configurations/`)
- `-o, --out_dir`: output directory (default: `Data/`)
- `--force-create-output-dir`: set to `1` to create missing output folders

When `--force-create-output-dir=1`, the program creates these subdirectories under `out_dir`:

- `Spectra_ascii`
- `Spectra_info`
- `Spectra_modelfile`
- `Spectra_plot`

---

## 5) Quick start (recommended first run)

Use the provided example configuration:

```bash
./build/specsim \
  --main_file main.cfg.aj_GRANscaledKallinger \
  --main_dir Configurations/examples_cfg/ \
  --noise_file noise_Kallinger2014.cfg \
  --out_dir Data/ \
  --force-create-output-dir 1
```

What this does:

- Uses one of the shipped random-model examples.
- Loads Kallinger noise settings.
- Creates output folders (if needed).
- Produces simulated spectra and associated metadata.

---

## 6) Understanding the main configuration file

A typical main configuration file contains:

1. **Forest specification** (`random` or `grid`) and related controls.
2. **Model name** plus optional model input file (e.g., an `.in` file).
3. **Template selection** (`NONE`, `all`, or specific template names).
4. Parameter block:
   - names,
   - `val_min`,
   - `val_max`,
   - random/grid selector line.
5. Observation setup (`Tobs`, `Cadence`, `Naverage`, `Nrealisation`).
6. Output behavior flags (`erase_old_files`, plots, model files, etc.).

### Important behavior flags

- `erase_old_files = 1`: restarts identifiers and overwrites combination/history outputs.
- `erase_old_files = 0`: appends to existing combinations and continues numbering.

---

## 7) Understanding the noise configuration file

`noise_Kallinger2014.cfg` supports multiple distribution types per parameter:

- `Uniform`
- `Gaussian`
- `Fix`

The file includes both `forest_type = random` and `forest_type = grid` sections. The simulator selects the matching one.

For random mode, the final line acts as a Gaussian uncertainty scaling coefficient (`kerror`) where relevant.

---

## 8) Common model names

The simulator currently handles model names such as:

- `generate_cfg_asymptotic_act_asym_Hgauss`
- `generate_cfg_from_refstar_HWscaled`
- `generate_cfg_from_refstar_HWscaled_GRANscaled`
- `generate_cfg_from_synthese_file_Wscaled_act_asym_a1ovGamma`
- `generate_cfg_from_synthese_file_Wscaled_a1a2a3asymovGamma`
- `generate_cfg_from_synthese_file_Wscaled_Alm`
- `generate_cfg_from_synthese_file_Wscaled_aj`
- `generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled`
- `generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014`
- `asymptotic_mm_v1`, `asymptotic_mm_v2`, `asymptotic_mm_v3`
- `asymptotic_mm_freeDp_numaxspread_curvepmodes_v1`
- `asymptotic_mm_freeDp_numaxspread_curvepmodes_v2`
- `asymptotic_mm_freeDp_numaxspread_curvepmodes_v3`
- `asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014`

If you are unsure which one to use, start from an existing file in `Configurations/examples_cfg/` and adjust parameters incrementally.

---

## 9) Output files and folders

By default (or under your chosen `out_dir`), you get:

- `Combinations.txt`: summary of generated parameter combinations.
- `Spectra_ascii/`: synthetic spectra data files.
- `Spectra_info/`: per-spectrum parameter/config snapshots.
- `Spectra_modelfile/`: optional model files for downstream tools.
- `Spectra_plot/`: plots when plotting is enabled.

---

## 10) Doxygen API documentation

To build the API docs:

```bash
cd docs
doxygen Doxyfile
```

Then open `docs/html/index.html`.

Use this **User Guide** for workflows and configuration, and Doxygen for low-level function/class references.

---

## 11) Troubleshooting

- **“Output directory does not exist”**:
  - Re-run with `--force-create-output-dir 1`, or create directories manually.
- **CMake cannot find Eigen on Linux**:
  - Export/include `EIGEN3_INCLUDE_DIR` before configuring.
- **gnuplot not found**:
  - Install gnuplot; CMake treats it as required.
- **Unexpected overwrite behavior**:
  - Verify `erase_old_files` in your main config.

---

## 12) Suggested learning path

1. Run the quick start command unchanged.
2. Inspect `Data/Combinations.txt` and one file from `Spectra_info/`.
3. Change only one parameter range in the example config.
4. Re-run and compare output distributions.
5. Move to a different model family once the first pipeline is understood.

