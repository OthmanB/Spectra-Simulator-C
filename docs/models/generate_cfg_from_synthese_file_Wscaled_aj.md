# generate_cfg_from_synthese_file_Wscaled_aj

## What this model is for

This model rescales a reference star (`.in`) to a new asymptotic p-mode pattern using `Dnu`, `epsilon`, and `delta0l_percent`. It then applies rotational splitting coefficients (`a2`..`a6`), height/width rescaling, and optional random spreads on heights, frequencies, and widths.

## Required inputs

- `main.cfg` line for `model_name` must include a path to a reference `.in` file (the `extra_params` argument).
- `template_files` line must exist but should be `NONE` (templates are not used by this model).
- No external noise configuration is required.
- Forest support: `random` and `grid`.

Example config: `Configurations/examples_cfg/main.cfg.aj`

## Parameter order (registry)

| # | Label | Meaning / notes |
|---|-------|-----------------|
| 1 | Dnu | Large separation (uHz). |
| 2 | epsilon | Asymptotic phase term. |
| 3 | delta0l_percent | Small separation in percent (sign is respected). |
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
| 14 | H_spread | Percent height jitter (uniform, applied per mode). |
| 15 | nu_spread | Percent frequency jitter (uniform, applied per mode). |
| 16 | Gamma_spread | Percent width jitter (uniform, applied per mode). |
| 17 | do_flat_noise | If >0, flatten the noise profile. |

## Physics assumptions and equations

- Frequencies are rescaled from the reference star using the asymptotic relation:
  - `nu(n,l) = (n + epsilon + l/2) * Dnu - l(l+1) * d0l`
- `d0l = delta0l_percent * Dnu / 100` (sign as provided)
- Heights are rescaled to match `HNR` relative to the reference star.
- Widths are scaled to match `Gamma_at_numax` and `a1 = a1ovGamma * Gamma_at_numax`.

## How to run

```bash
./build/specsim \
  --main_file main.cfg.aj \
  --main_dir Configurations/examples_cfg \
  --out_dir /tmp/specsim-out \
  --force-create-output-dir 1
```

## Outputs

- `Combinations.txt`
- `Spectra_ascii/`, `Spectra_info/`, optional `Spectra_modelfile/`, `Spectra_plot/`

## Common pitfalls and tips

- `extra_params` must point to a valid reference `.in` file.
- Set `template_files` to `NONE` (required by the parser even if unused).
- `delta0l_percent` sign is used as provided (positive and negative values are allowed).
- Spreads (`H_spread`, `nu_spread`, `Gamma_spread`) are percentages (0-100).
