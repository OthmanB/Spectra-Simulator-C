# generate_cfg_from_synthese_file_Wscaled_Alm

## What this model is for

This model uses a reference star file (`.in`) as a template for frequencies and mode profiles, then rescales heights and widths to target a desired height-to-noise ratio and line width at `numax`. It also applies an activity term (Alm) using the `epsilon_nl`, `theta0`, and `delta` parameters.

## Required inputs

- `main.cfg` line for `model_name` must include a path to a reference `.in` file (the `extra_params` argument).
- `template_files` line must exist but should be `NONE` (templates are not used by this model).
- No external noise configuration is required.
- Forest support: `random` and `grid`.

Example config: `Configurations/examples_cfg/main.cfg.Alm`

## Parameter order (registry)

| # | Label | Meaning / notes |
|---|-------|-----------------|
| 1 | HNR | Target max height-to-noise ratio relative to the reference star. |
| 2 | a1ovGamma | Target ratio `a1 / Gamma` at `numax` (sets splitting from the width). |
| 3 | Gamma_at_numax | Mode width at `numax` (used to scale widths). |
| 4 | epsilon_nl | Activity term amplitude for Alm. |
| 5 | theta0 | Activity term parameter (used by Alm). |
| 6 | delta | Activity term parameter (used by Alm). |
| 7 | a3 | Latitudinal rotation coefficient (units consistent with reference frequencies). |
| 8 | beta_asym | Line-profile asymmetry parameter. |
| 9 | i | Inclination (degrees). |

## Physics assumptions and equations

- Frequencies and baseline mode profiles are taken from the reference `.in` file.
- Heights are rescaled so that the maximum height-to-noise ratio matches `HNR`.
- Widths are scaled to satisfy `Gamma_at_numax` and `a1 = a1ovGamma * Gamma_at_numax`.
- An activity perturbation is applied using `epsilon_nl`, `theta0`, and `delta` (Alm model).

## How to run

```bash
./build/specsim \
  --main_file main.cfg.Alm \
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
- Ensure `HNR`, `Gamma_at_numax`, and `i` are non-negative.
