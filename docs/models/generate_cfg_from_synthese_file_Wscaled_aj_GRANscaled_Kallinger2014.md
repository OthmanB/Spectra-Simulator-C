# generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014

## What this model is for

This model extends `generate_cfg_from_synthese_file_Wscaled_aj` by adding a granulation background based on the Kallinger+2014 noise model. It rescales a reference `.in` file to a new p-mode pattern and derives the noise profile from `noise_Kallinger2014.cfg`.

## Required inputs

- `main.cfg` line for `model_name` must include a path to a reference `.in` file (the `extra_params` argument).
- `template_files` line must exist but should be `NONE` (templates are not used by this model).
- `noise_Kallinger2014.cfg` is required (passed via `--noise_file` or resolved from `Configurations/`).
- Forest support: `random` only.

Example config: `Configurations/examples_cfg/main.cfg.aj_GRANscaledKallinger`

## Parameter order (registry)

The noise parameters are appended automatically from `noise_Kallinger2014.cfg` and are not listed in the main cfg labels line.

| # | Label | Meaning / notes |
|---|-------|-----------------|
| 1 | Dnu | Large separation (uHz). |
| 2 | epsilon | Asymptotic phase term. |
| 3 | delta0l_percent | Small separation in percent (input magnitude; internal sign is negative). |
| 4 | HNR | Target max height-to-noise ratio relative to the reference star. |
| 5 | a1ovGamma | Target ratio `a1 / Gamma` at `numax`. |
| 6 | Gamma_at_numax | Mode width at `numax`. |
| 7 | a2 | Rotation/asphericity coefficient. |
| 8 | a3 | Latitudinal rotation coefficient. |
| 9 | a4 | Rotation/asphericity coefficient. |
| 10 | a5 | Rotation/asphericity coefficient. |
| 11 | a6 | Rotation/asphericity coefficient. |
| 12 | beta_asym | Line-profile asymmetry parameter. |
| 13 | i | Inclination (degrees). |
| 14 | numax_spread | Percent spread applied to `numax`. |
| 15 | H_spread | Percent height jitter (uniform, applied per mode). |
| 16 | nu_spread | Percent frequency jitter (uniform, applied per mode). |
| 17 | k_Agran | Kallinger noise parameter (from noise cfg). |
| 18 | s_Agran | Kallinger noise parameter (from noise cfg). |
| 19 | k_taugran | Kallinger noise parameter (from noise cfg). |
| 20 | s_taugran | Kallinger noise parameter (from noise cfg). |
| 21 | c0 | Kallinger noise parameter (from noise cfg). |
| 22 | ka | Kallinger noise parameter (from noise cfg). |
| 23 | ks | Kallinger noise parameter (from noise cfg). |
| 24 | k1 | Kallinger noise parameter (from noise cfg). |
| 25 | s1 | Kallinger noise parameter (from noise cfg). |
| 26 | c1 | Kallinger noise parameter (from noise cfg). |
| 27 | k2 | Kallinger noise parameter (from noise cfg). |
| 28 | s2 | Kallinger noise parameter (from noise cfg). |
| 29 | c2 | Kallinger noise parameter (from noise cfg). |
| 30 | N0 | White noise level (from noise cfg). |

Note: `Tobs` and `Cadence` are appended internally for the Kallinger noise conversion and are not part of the main cfg labels.

## Physics assumptions and equations

- Frequencies follow the asymptotic p-mode relation:
  - `nu(n,l) = (n + epsilon + l/2) * Dnu - l(l+1) * d0l`
  - `d0l = delta0l_percent * Dnu / 100` (internal sign is negative)
- `numax` is derived from `Dnu` (Stello 2009 relation) with an optional `numax_spread`.
- Granulation noise follows Kallinger+2014 and is converted to Harvey-like terms internally.

## How to run

```bash
./build/specsim \
  --main_file main.cfg.aj_GRANscaledKallinger \
  --main_dir Configurations/examples_cfg \
  --noise_file noise_Kallinger2014.cfg \
  --out_dir /tmp/specsim-out \
  --force-create-output-dir 1
```

## Outputs

- `Combinations.txt`
- `Spectra_ascii/`, `Spectra_info/`, optional `Spectra_modelfile/`, `Spectra_plot/`

## Common pitfalls and tips

- `extra_params` must point to a valid reference `.in` file.
- Ensure `noise_Kallinger2014.cfg` is available and readable.
- `delta0l_percent` is treated as a magnitude; negative values will be flipped internally.
- Spreads (`numax_spread`, `H_spread`, `nu_spread`) are percentages (0-100).
