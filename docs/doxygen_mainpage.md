# Solar-Like Spectrum Simulator Documentation

This is the API reference for `specsim`, a C++17 CLI tool that generates synthetic asteroseismic power spectra from configuration files.

If you are looking for practical usage and configuration examples, start with the user documentation:

- Online docs (MkDocs): https://othmanb.github.io/Spectra-Simulator-C/
- Repository User Guide: `docs/User_Guide.md`
- Model manuals: `docs/models/README.md`

## Quick start

Build (out-of-tree recommended):

```bash
cmake -S . -B build
cmake --build build -j
```

Show CLI help:

```bash
./build/specsim --help
```

Run an example configuration:

```bash
./build/specsim \
  --main_file main.cfg.aj \
  --main_dir Configurations/examples_cfg \
  --out_dir /tmp/specsim-out \
  --force-create-output-dir 1
```

## What this API reference covers

- Mode/profile builders (Lorentzians, visibilities, rotation-related utilities)
- Model generators that produce intermediate mode/noise configuration files
- Spectrum builders that assemble model components onto a frequency grid

For end-to-end semantics (configuration labels, model parameters, expected outputs), prefer the User Guide and model manuals linked above.
