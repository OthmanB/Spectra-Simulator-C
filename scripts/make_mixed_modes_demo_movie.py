#!/usr/bin/env python3
"""Generate a short mixed-modes demo movie from specsim.

This script runs the mixed-mode example model for a sequence of parameters that mimic
an evolution from high numax (subgiant-ish) to low numax (RGB-ish), then renders a
2-panel figure per frame:

  Left: zoomed power spectrum around numax (with full-spectrum inset)
  Right: echelle diagram (built from the model spectrum column)

It is meant for outreach/demo purposes (LinkedIn, docs). It does not try to be a
physically-perfect evolutionary track.

Dependencies (recommended in a venv):
  python3 -m venv .venv
  source .venv/bin/activate
  pip install numpy matplotlib

Requirements:
  - A built specsim binary (default: ./build/specsim)
  - ffmpeg installed (optional, only needed to encode mp4)
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


def _require_imports():
    try:
        import numpy as np  # noqa: F401
        import matplotlib.pyplot as plt  # noqa: F401
    except Exception as e:
        raise SystemExit(
            "Missing python deps. Create a venv and install: pip install numpy matplotlib\n"
            f"Original import error: {e}"
        )


def _setup_logging(level: str) -> logging.Logger:
    lvl = getattr(logging, level.upper(), None)
    if not isinstance(lvl, int):
        raise SystemExit(f"Invalid --log-level: {level}")
    logging.basicConfig(
        level=lvl,
        format="%(asctime)s %(levelname)s %(funcName)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("specsim_demo")


@dataclass
class FrameParams:
    numax_uhz: float
    dnu_uhz: float
    dp1_s: float
    nurot_env_uhz: float
    nurot_core_uhz: float
    max_gamma_uhz: float


def stello_dnu_from_numax(numax_uhz: float) -> float:
    # Stello+2009: Dnu ~ beta0 * numax^alpha0
    beta0 = 0.263
    alpha0 = 0.77
    return beta0 * (numax_uhz ** alpha0)


def dp1_from_numax_piecewise(numax_uhz: float) -> float:
    # Simple piecewise mapping for demo purposes.
    # - High numax: DP1 ~ 800s
    # - Rapid drop to ~150s by numax~600
    # - Slow evolution down to ~110s by numax~400
    if numax_uhz >= 1000.0:
        return 800.0
    if numax_uhz >= 600.0:
        # 1000 -> 800, 600 -> 150
        x0, y0 = 1000.0, 800.0
        x1, y1 = 600.0, 150.0
        t = (numax_uhz - x1) / (x0 - x1)
        return y1 + t * (y0 - y1)
    if numax_uhz <= 400.0:
        return 110.0
    # 600 -> 150, 400 -> 110
    x0, y0 = 600.0, 150.0
    x1, y1 = 400.0, 110.0
    t = (numax_uhz - x1) / (x0 - x1)
    return y1 + t * (y0 - y1)


def log_interp(a: float, b: float, t: float) -> float:
    return math.exp(math.log(a) * (1.0 - t) + math.log(b) * t)


def build_frames(nframes: int, numax_start: float, numax_end: float) -> list[FrameParams]:
    frames: list[FrameParams] = []
    for i in range(nframes):
        t = 0.0 if nframes <= 1 else i / (nframes - 1)
        numax = log_interp(numax_start, numax_end, t)
        dnu = stello_dnu_from_numax(numax)
        dp1 = dp1_from_numax_piecewise(numax)
        nurot_env = (1.0 - t) * 1.0 + t * 0.1
        # Keep ratio simple for demo; can be swapped to a ramped ratio if desired.
        nurot_core = 5.0 * nurot_env
        max_gamma = (1.0 - t) * 0.8 + t * 0.1
        frames.append(
            FrameParams(
                numax_uhz=numax,
                dnu_uhz=dnu,
                dp1_s=dp1,
                nurot_env_uhz=nurot_env,
                nurot_core_uhz=nurot_core,
                max_gamma_uhz=max_gamma,
            )
        )
    return frames


def _split_comment(line: str) -> tuple[str, str]:
    if "#" in line:
        pre, post = line.split("#", 1)
        return pre.rstrip(), "#" + post
    return line.rstrip(), ""


def _fmt_row(vals: list[float], comment: str) -> str:
    # Keep it compact and readable; specsim parser is whitespace-based.
    s = " ".join(f"{v:.8g}" for v in vals)
    return (s + (" " if comment else "") + comment).rstrip() + "\n"


def write_cfg_for_frame(
    base_cfg_path: Path,
    out_cfg_path: Path,
    template_name: str,
    frame: FrameParams,
    *,
    q: float,
    beta_p_star: float,
    epsilon: float,
    delta0l_percent: float,
    alpha: float,
    snr: float,
    vl1: float,
    vl2: float,
    vl3: float,
    numax_spread_percent: float,
    nmax_spread: float,
    h0_spread_percent: float,
    hfactor: float,
    wfactor: float,
    tobs_days: float,
    cadence_s: float,
    naverage: int,
    nrealisation: int,
    doplots: int,
    write_inmodel: int,
    limit_data_range: int,
    do_modelfiles: int,
) -> None:
    lines = base_cfg_path.read_text(encoding="utf-8").splitlines(True)

    # Find the key block by scanning for the labels line (first non-comment that contains 'nurot_env').
    # Expected structure:
    #   forest line
    #   model name line
    #   template line
    #   labels line
    #   val_min line
    #   val_max line
    #   step line
    idx_labels = None
    for i, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("nurot_env"):
            idx_labels = i
            break
    if idx_labels is None or idx_labels < 3:
        raise SystemExit(f"Could not locate labels line in cfg: {base_cfg_path}")

    idx_template = idx_labels - 1
    idx_model = idx_labels - 2
    idx_forest = idx_labels - 3
    idx_val_min = idx_labels + 1
    idx_val_max = idx_labels + 2
    idx_step = idx_labels + 3

    # Patch forest to run exactly one sample.
    forest_pre, forest_c = _split_comment(lines[idx_forest])
    forest_tokens = forest_pre.split()
    if not forest_tokens:
        raise SystemExit("Invalid forest line")
    forest_tokens[0] = "random"
    if len(forest_tokens) == 1:
        forest_tokens.append("1")
    else:
        forest_tokens[1] = "1"
    lines[idx_forest] = " ".join(forest_tokens) + (" " if forest_c else "") + forest_c.strip("\n") + "\n"

    # Patch template selection.
    tmpl_pre, tmpl_c = _split_comment(lines[idx_template])
    _ = tmpl_pre  # not used
    lines[idx_template] = f"{template_name}" + (" " if tmpl_c else "") + tmpl_c.strip("\n") + "\n"

    # Parse labels and build val_min/val_max vectors.
    labels_pre, _labels_c = _split_comment(lines[idx_labels])
    labels = labels_pre.split()
    if not labels:
        raise SystemExit("Empty labels line")

    # Build desired values.
    values_by_label: dict[str, float] = {
        "nurot_env": frame.nurot_env_uhz,
        "nurot_core": frame.nurot_core_uhz,
        "a2_l1_core": 0.0,
        "a2_l1_env": 0.0,
        "a2_l2_env": 0.0,
        "a2_l3_env": 0.0,
        "a3_l2_env": 0.0,
        "a3_l3_env": 0.0,
        "a4_l2_env": 0.0,
        "a4_l3_env": 0.0,
        "a5_l3_env": 0.0,
        "a6_l3_env": 0.0,
        "Dnu": frame.dnu_uhz,
        "epsilon": epsilon,
        "delta0l_percent": delta0l_percent,
        "beta_p_star": beta_p_star,
        "nmax_spread": nmax_spread,
        "DP1": frame.dp1_s,
        "alpha": alpha,
        "q": q,
        "SNR": snr,
        "maxGamma": frame.max_gamma_uhz,
        "numax_spread": numax_spread_percent,
        "Vl1": vl1,
        "Vl2": vl2,
        "Vl3": vl3,
        "H0_spread": h0_spread_percent,
        "Hfactor": hfactor,
        "Wfactor": wfactor,
    }

    vals = []
    for lab in labels:
        if lab not in values_by_label:
            raise SystemExit(f"Label '{lab}' not handled by demo script; update values_by_label")
        vals.append(float(values_by_label[lab]))

    # Set val_min == val_max for determinism.
    _vm_pre, vm_c = _split_comment(lines[idx_val_min])
    _vx_pre, vx_c = _split_comment(lines[idx_val_max])
    lines[idx_val_min] = _fmt_row(vals, vm_c)
    lines[idx_val_max] = _fmt_row(vals, vx_c)

    # Set step line to 0 (all constants) to avoid any randomization.
    _step_pre, step_c = _split_comment(lines[idx_step])
    step_vals = [0.0 for _ in labels]
    lines[idx_step] = _fmt_row(step_vals, step_c)

    # Patch the bottom control lines. The example cfg may be missing the final lines expected by the parser.
    # Find the 'Tobs   Cadence' header.
    idx_tobs_hdr = None
    for i, raw in enumerate(lines):
        if raw.strip().startswith("Tobs") and "Cadence" in raw:
            idx_tobs_hdr = i
            break
    if idx_tobs_hdr is None:
        raise SystemExit("Could not find Tobs/Cadence header")
    idx_tobs_vals = idx_tobs_hdr + 1
    # Values: Tobs Cadence Naverage Nspectra? (cfg.Nspectra) and Nrealisation
    # In current code: cfg.Nspectra=tmp[2], cfg.Nrealisation=tmp[3]
    lines[idx_tobs_vals] = f"{tobs_days:.8g} {cadence_s:.8g} {int(naverage)} {int(nrealisation)}\n"

    # Ensure following lines exist: erase_old_files, doplots, write_inmodel, limit_data_range, do_modelfiles.
    tail_needed = 5
    tail_start = idx_tobs_vals + 1
    # Extend file if too short.
    while len(lines) < tail_start + tail_needed:
        lines.append("\n")

    lines[tail_start + 0] = "1\n"  # erase_old_files
    lines[tail_start + 1] = f"{int(doplots)}\n"
    lines[tail_start + 2] = f"{int(write_inmodel)}\n"
    lines[tail_start + 3] = f"{int(limit_data_range)}\n"
    lines[tail_start + 4] = f"{int(do_modelfiles)}\n"

    out_cfg_path.write_text("".join(lines), encoding="utf-8")


def run_specsim(
    specsim: Path,
    cfg_path: Path,
    out_dir: Path,
    seed: int,
    noise_file: str,
    *,
    cwd: Path,
    timeout_s: int,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(specsim),
        "--main_file",
        str(cfg_path),
        "--main_dir",
        ".",
        "--noise_file",
        noise_file,
        "--out_dir",
        str(out_dir),
        "--force-create-output-dir",
        "1",
        "--seed",
        str(seed),
    ]
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=str(cwd),
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired as e:
        raise SystemExit(f"specsim timed out after {timeout_s}s: {e.cmd}")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"specsim failed with exit code {e.returncode}: {e.cmd}")
    # We run one combination and one realisation; file naming is 0000001.0.data
    data_file = out_dir / "Spectra_ascii" / "0000001.0.data"
    if not data_file.exists():
        raise SystemExit(f"Expected spectrum file not found: {data_file}")
    return data_file


def read_spectrum_3col(path: Path):
    import numpy as np

    xs = []
    ys = []
    ms = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
            ms.append(float(parts[2]))
    return np.asarray(xs), np.asarray(ys), np.asarray(ms)


def read_template_table(path: Path):
    """Read a *.template file.

    Returns: (numax_ref, dnu_ref, freq[], height[], width[])
    """
    import numpy as np

    numax_ref = None
    dnu_ref = None
    freq = []
    height = []
    width = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("numax_ref"):
                numax_ref = float(line.split("=", 1)[1].strip())
                continue
            if line.startswith("Dnu_ref"):
                dnu_ref = float(line.split("=", 1)[1].strip())
                continue
            # Data row: 3 numeric columns.
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                f0 = float(parts[0])
                h0 = float(parts[1])
                w0 = float(parts[2])
            except Exception:
                continue
            freq.append(f0)
            height.append(h0)
            width.append(w0)

    if numax_ref is None or dnu_ref is None or len(freq) < 5:
        raise ValueError("invalid template header or too few rows")

    freq = np.asarray(freq, dtype=float)
    height = np.asarray(height, dtype=float)
    width = np.asarray(width, dtype=float)
    # Ensure sorted by frequency.
    order = np.argsort(freq)
    return float(numax_ref), float(dnu_ref), freq[order], height[order], width[order]


def pick_template_auto(repo_root: Path, *, edge_dnu: float = 6.0) -> str:
    """Pick a template with the fastest height falloff away from numax.

    Heuristic: minimize max(height at numax +/- edge_dnu*Dnu_ref) / max(height).
    """
    import numpy as np

    templates_dir = repo_root / "Configurations" / "templates"
    candidates = sorted(templates_dir.glob("*.template"))
    if not candidates:
        raise SystemExit(f"No templates found under {templates_dir}")

    best = None
    for p in candidates:
        try:
            numax_ref, dnu_ref, f, h, _w = read_template_table(p)
        except Exception:
            continue
        hmax = float(np.max(h))
        if not math.isfinite(hmax) or hmax <= 0:
            continue
        f_lo = numax_ref - edge_dnu * dnu_ref
        f_hi = numax_ref + edge_dnu * dnu_ref
        h_lo = float(np.interp(f_lo, f, h, left=h[0], right=h[-1]))
        h_hi = float(np.interp(f_hi, f, h, left=h[0], right=h[-1]))
        score = max(h_lo, h_hi) / hmax
        # Prefer smaller score; break ties by larger hmax.
        key = (score, -hmax)
        if best is None or key < best[0]:
            best = (key, p.name)

    if best is None:
        raise SystemExit("No valid templates found (all malformed?)")
    return best[1]


def gaussian_kernel1d(sigma: float, radius: int | None = None):
    import numpy as np

    if sigma <= 0:
        return np.asarray([1.0])
    if radius is None:
        radius = int(math.ceil(4.0 * sigma))
    x = np.arange(-radius, radius + 1)
    k = np.exp(-0.5 * (x / sigma) ** 2)
    k /= np.sum(k)
    return k


def smooth_1d(y, sigma_pix: float):
    import numpy as np

    k = gaussian_kernel1d(sigma_pix)
    return np.convolve(y, k, mode="same")


def smooth_2d_separable(z, sigma_y: float, sigma_x: float):
    import numpy as np

    out = z.copy()

    if sigma_x > 0:
        kx = gaussian_kernel1d(sigma_x)
        for i in range(out.shape[0]):
            row = out[i, :]
            m = np.isfinite(row)
            if not np.any(m):
                continue
            # Fill NaNs for convolution, but keep mask.
            filled = row.copy()
            filled[~m] = np.nanmedian(row[m])
            out[i, :] = np.convolve(filled, kx, mode="same")

    if sigma_y > 0:
        ky = gaussian_kernel1d(sigma_y)
        for j in range(out.shape[1]):
            col = out[:, j]
            m = np.isfinite(col)
            if not np.any(m):
                continue
            filled = col.copy()
            filled[~m] = np.nanmedian(col[m])
            out[:, j] = np.convolve(filled, ky, mode="same")

    return out


def read_modes_from_info_file(path: Path):
    """Read modes from the per-run Spectra_info/*.in file.

    Expected section:
      # Input mode parameters. degree / freq / H / W / ...
      <numeric rows>
      # Configuration of mode parameters. H0 , tau_0 , ... (noise section)

    Returns Nx2 array [l, nu].
    """
    import numpy as np

    if not path.exists():
        return np.zeros((0, 2))

    l_list = []
    nu_list = []
    in_modes = False
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "Input mode parameters" in line:
                    in_modes = True
                    continue
                if "Configuration of mode parameters" in line and in_modes:
                    break
                continue
            if not in_modes:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                l_list.append(int(round(float(parts[0]))))
                nu_list.append(float(parts[1]))
            except Exception:
                continue

    if not l_list:
        return np.zeros((0, 2))
    return np.column_stack([np.asarray(l_list, dtype=int), np.asarray(nu_list, dtype=float)])


def render_frame(
    data_file: Path,
    modes_l_nu,
    out_png: Path,
    *,
    numax: float,
    dnu: float,
    nurot_env: float,
    nurot_core: float,
    dp1: float,
    max_gamma: float,
    zoom_dnu_halfwidth: float,
):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import transforms
    from matplotlib.lines import Line2D
    from matplotlib.gridspec import GridSpec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    freq, spec, model = read_spectrum_3col(data_file)
    df = float(freq[1] - freq[0]) if len(freq) > 1 else 1.0

    # Left panel range
    fmin = max(0.0, numax - zoom_dnu_halfwidth * dnu)
    fmax = numax + zoom_dnu_halfwidth * dnu
    m = (freq >= fmin) & (freq <= fmax)

    fig = plt.figure(figsize=(12.5, 5.2), dpi=150)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.05, 1.0])
    axL = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1])

    # Avoid per-frame layout jitter from tight_layout().
    # Use fixed margins and fixed colorbar axes.
    fig.subplots_adjust(left=0.065, right=0.92, bottom=0.11, top=0.90, wspace=0.22)

    # Left panel: show model (clean) + noisy (faint)
    # Smooth the noisy PSD for display (very mild): ~0.01 * maxGamma.
    sigma_freq = max(1e-6, 0.01 * max_gamma)
    sigma_pix = sigma_freq / max(df, 1e-12)
    spec_s = smooth_1d(spec, sigma_pix)

    # Plot: noisy spectrum first, then model on top for visibility.
    axL.plot(freq[m], spec_s[m], color="0.25", lw=0.9, alpha=0.85, label="_nolegend_", zorder=2)
    axL.plot(freq[m], model[m], color="#1f77b4", lw=1.0, alpha=0.95, label="_nolegend_", zorder=3)
    axL.set_xlim(fmin, fmax)
    axL.set_yscale("log")
    axL.set_xlabel(r"Frequency $\,\mu$Hz")
    axL.set_ylabel(r"PSD (ppm$^2\,\mu$Hz$^{-1}$)")
    axL.set_title(r"Zoom: $\nu_{max}\pm %.1f\,\Delta\nu$" % zoom_dnu_halfwidth)

    # Overlay mode frequencies if available
    # Mode markers at the bottom of the zoomed spectrum.
    if modes_l_nu.shape[0] > 0:
        colors = {0: "#000000", 1: "#1f77b4", 2: "#d62728", 3: "#8c564b"}
        trans = transforms.blended_transform_factory(axL.transData, axL.transAxes)
        for l in (0, 1, 2, 3):
            mm = (modes_l_nu[:, 0] == l) & (modes_l_nu[:, 1] >= fmin) & (modes_l_nu[:, 1] <= fmax)
            if np.any(mm):
                axL.scatter(
                    modes_l_nu[mm, 1],
                    [0.05] * int(np.count_nonzero(mm)),
                    s=22,
                    color=colors[l],
                    marker="o",
                    edgecolors="none",
                    alpha=0.95,
                    zorder=4,
                    transform=trans,
                    clip_on=False,
                )

        legend_handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[l], markeredgecolor="none", markersize=5, label=f"l={l}")
            for l in (0, 1, 2, 3)
        ]
        axL.legend(legend_handles, [f"l={l}" for l in (0, 1, 2, 3)], loc="upper left", frameon=True, fontsize=8)

    # Full spectrum inset
    # Smaller inset, pushed to the corner to avoid hiding modes.
    ax_in = inset_axes(axL, width="38%", height="30%", loc="upper right", borderpad=0.6)
    # Inset: show noisy PSD (smoothed) in log-log to reveal the low-frequency background.
    # Use a fixed smoothing scale (0.05 microHz) for consistent appearance.
    sigma_inset_freq = 0.05
    sigma_inset_pix = sigma_inset_freq / max(df, 1e-12)
    spec_in = smooth_1d(spec, sigma_inset_pix)
    ax_in.plot(freq, spec_in, color="0.25", lw=0.7)
    ax_in.set_yscale("log")
    ax_in.set_xscale("log")
    ax_in.set_xlim(max(freq[1], 1e-3), min(freq[-1], 5000.0))
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.axvspan(fmin, fmax, color="#d62728", alpha=0.12)

    # Right panel: echelle from the noisy spectrum with mild 2D smoothing.
    # Build stacked orders between fmin..fmax (use a wider window for echelle by default)
    ech_min = max(0.0, numax - 5.0 * dnu)
    ech_max = numax + 5.0 * dnu
    mm_e = (freq >= ech_min) & (freq <= ech_max)
    freq_e = freq[mm_e]
    spec_e = spec[mm_e]

    n0 = int(math.floor(ech_min / dnu))
    n1 = int(math.ceil(ech_max / dnu))
    n_orders = max(1, n1 - n0)
    nx = 480
    xgrid = np.linspace(0.0, dnu, nx)
    z = np.full((n_orders, nx), np.nan)

    for k in range(n_orders):
        a = (n0 + k) * dnu
        b = a + dnu
        mk = (freq_e >= a) & (freq_e < b)
        if not np.any(mk):
            continue
        x = freq_e[mk] - a
        y = spec_e[mk]
        # ensure sorted
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        z[k, :] = np.interp(xgrid, x, y, left=np.nan, right=np.nan)

    # Extend the x-axis a bit beyond [0, Dnu] by wrapping columns. This helps show mixed modes near boundaries.
    pad_frac = 0.25
    pad = max(1, int(pad_frac * nx))
    z_ext = np.concatenate([z[:, -pad:], z, z[:, :pad]], axis=1)
    x_ext = np.concatenate([
        xgrid[-pad:] - dnu,
        xgrid,
        xgrid[:pad] + dnu,
    ])

    # Mild vertical+horizontal smoothing (image-like). Values are in pixels.
    z_sm = smooth_2d_separable(z_ext, sigma_y=1.0, sigma_x=1.2)

    # Contrast: log10 with percentile clipping
    z_log = np.log10(np.clip(z_sm, 1e-30, np.nanmax(z_sm)))
    vmin = np.nanpercentile(z_log, 15)
    vmax = np.nanpercentile(z_log, 99.5)
    # Prefer a blue-green-yellow-red palette (turbo); fall back to Spectral_r if unavailable.
    cmap_name = "turbo"
    try:
        plt.get_cmap(cmap_name)
    except Exception:
        cmap_name = "Spectral_r"

    im = axR.imshow(
        z_log,
        aspect="auto",
        origin="lower",
        extent=[float(x_ext[0]), float(x_ext[-1]), n0, n0 + n_orders],
        cmap=cmap_name,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    axR.set_xlabel(r"$\nu\;\mathrm{mod}\;\Delta\nu\,(\mu\mathrm{Hz})$")
    axR.set_ylabel("Radial order")
    axR.set_title("Echelle Diagram")


    # Small annotation bar
    title = (
        rf"$\nu_{{\max}}\approx {numax:6.1f}\,\mu\mathrm{{Hz}}"
        rf"\quad \Delta\nu\approx {dnu:5.2f}\,\mu\mathrm{{Hz}}"
        rf"\quad \Delta\Pi_1\approx {dp1:5.0f}\,\mathrm{{s}}"
        rf"\quad \nu_{{\mathrm{{env}}}}/\nu_{{\mathrm{{core}}}}={nurot_env:.2f}/{nurot_core:.2f}\,\mu\mathrm{{Hz}}"
        rf"\quad \Gamma(\nu_{{\max}})\approx {max_gamma:.2f}\,\mu\mathrm{{Hz}}$"
    )
    fig.suptitle(title, fontsize=10, y=0.98)

    # Fixed-position colorbar to prevent layout shifts as tick labels change.
    cax = fig.add_axes([0.935, 0.18, 0.015, 0.64])
    fig.colorbar(im, cax=cax, label=r"$\log_{10}$ PSD")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def encode_mp4(frames_dir: Path, fps: int, out_mp4: Path, *, timeout_s: int) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise SystemExit("ffmpeg not found in PATH. Install ffmpeg or use --no-encode")
    cmd = [
        ffmpeg,
        "-y",
        "-framerate",
        str(fps),
        "-i",
        str(frames_dir / "frame_%04d.png"),
        "-pix_fmt",
        "yuv420p",
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        str(out_mp4),
    ]
    try:
        subprocess.run(cmd, check=True, timeout=timeout_s)
    except subprocess.TimeoutExpired as e:
        raise SystemExit(f"ffmpeg timed out after {timeout_s}s: {e.cmd}")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"ffmpeg failed with exit code {e.returncode}: {e.cmd}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--specsim", default="./build/specsim", help="Path to specsim binary")
    ap.add_argument(
        "--base-cfg",
        default="Configurations/examples_cfg/main.cfg.freeDP_curvepmodes.v3_GRANscaled",
        help="Base cfg to patch per frame",
    )
    ap.add_argument(
        "--template",
        default="auto",
        help="Template filename under Configurations/templates (or 'auto' to pick a steep envelope)",
    )
    ap.add_argument("--noise-file", default="noise_Kallinger2014.cfg", help="Noise cfg filename")
    ap.add_argument("--out", default="mixed_modes_demo.mp4", help="Output mp4 path")
    ap.add_argument(
        "--workdir",
        default="_mixed_modes_demo_work",
        help="Working directory for intermediate cfg/runs/frames (set empty to use a temp dir)",
    )
    ap.add_argument("--overwrite-workdir", action="store_true", help="Delete workdir if it already exists")
    ap.add_argument("--resume", action="store_true", help="Allow existing workdir and skip existing frames")
    ap.add_argument("--frames", type=int, default=200, help="Number of frames")
    ap.add_argument("--frame-start", type=int, default=0, help="Start frame index (0-based)")
    ap.add_argument("--frame-count", type=int, default=None, help="Number of frames to render (default: all)")
    ap.add_argument("--fps", type=int, default=10, help="Frames per second")
    ap.add_argument("--seed", type=int, default=123, help="Base RNG seed")
    ap.add_argument(
        "--seed-mode",
        choices=["fixed", "dancing"],
        default="fixed",
        help="fixed: same seed each frame; dancing: seed+frame",
    )
    ap.add_argument("--numax-start", type=float, default=1300.0)
    ap.add_argument("--numax-end", type=float, default=400.0)
    ap.add_argument("--zoom-dnu", type=float, default=4.0, help="Half-width of zoom window in units of Dnu")
    ap.add_argument("--no-encode", action="store_true", help="Only render PNG frames")
    ap.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)")
    ap.add_argument("--specsim-timeout", type=int, default=120, help="Timeout per specsim run (seconds)")
    ap.add_argument("--ffmpeg-timeout", type=int, default=600, help="Timeout for ffmpeg encoding (seconds)")

    # Fixed demo params
    ap.add_argument("--vl1", type=float, default=1.5)
    ap.add_argument("--vl2", type=float, default=0.5)
    ap.add_argument("--vl3", type=float, default=0.08)
    ap.add_argument("--q", type=float, default=0.20)
    ap.add_argument("--beta-p", type=float, default=0.02)
    ap.add_argument("--epsilon", type=float, default=0.25)
    ap.add_argument("--delta0l-percent", type=float, default=1.0)
    ap.add_argument("--alpha", type=float, default=0.3)
    ap.add_argument("--snr", type=float, default=100.0)
    ap.add_argument("--numax-spread", type=float, default=0.0)
    ap.add_argument("--nmax-spread", type=float, default=0.0)
    ap.add_argument("--h0-spread", type=float, default=0.0)
    ap.add_argument("--hfactor", type=float, default=1.0)
    ap.add_argument("--wfactor", type=float, default=0.5)
    ap.add_argument("--tobs", type=float, default=1200.0)
    ap.add_argument("--cadence", type=float, default=120.0)
    ap.add_argument("--naverage", type=int, default=20)

    args = ap.parse_args()

    log = _setup_logging(args.log_level)

    # Import heavy plotting deps after argparse so `--help` works without numpy/matplotlib.
    _require_imports()

    repo_root = Path(__file__).resolve().parents[1]
    specsim = (repo_root / args.specsim).resolve() if not os.path.isabs(args.specsim) else Path(args.specsim)
    if not specsim.exists():
        raise SystemExit(f"specsim not found: {specsim}")

    base_cfg = (repo_root / args.base_cfg).resolve() if not os.path.isabs(args.base_cfg) else Path(args.base_cfg)
    if not base_cfg.exists():
        raise SystemExit(f"base cfg not found: {base_cfg}")

    template_name = args.template
    if template_name.lower() == "auto":
        template_name = pick_template_auto(repo_root)
        log.info("Auto-selected template: %s", template_name)

    out_mp4 = Path(args.out).resolve()

    workdir = Path(args.workdir).resolve() if args.workdir else None
    created_temp_workdir = False
    if workdir is None:
        # Use a non-auto-cleaned temp dir so `--no-encode` / partial frame runs can be
        # encoded manually afterward.
        workdir = Path(tempfile.mkdtemp(prefix="specsim-mixedmodes-demo-"))
        created_temp_workdir = True
        log.info("Using temporary workdir: %s", workdir)
    else:
        if workdir.exists():
            if args.overwrite_workdir:
                shutil.rmtree(workdir)
                workdir.mkdir(parents=True, exist_ok=True)
            elif not args.resume:
                raise SystemExit(f"workdir already exists: {workdir} (use --overwrite-workdir or --resume)")
        else:
            workdir.mkdir(parents=True, exist_ok=True)

    try:
        frames_dir = workdir / "frames"
        runs_dir = workdir / "runs"
        cfg_dir = workdir / "cfg"
        frames_dir.mkdir(parents=True, exist_ok=True)
        runs_dir.mkdir(parents=True, exist_ok=True)
        cfg_dir.mkdir(parents=True, exist_ok=True)

        frames = build_frames(args.frames, args.numax_start, args.numax_end)
        start = max(0, args.frame_start)
        end = len(frames)
        if args.frame_count is not None:
            end = min(end, start + max(0, args.frame_count))

        for i in range(start, end):
            fr = frames[i]
            cfg_path = cfg_dir / f"frame_{i:04d}.cfg"
            write_cfg_for_frame(
                base_cfg,
                cfg_path,
                template_name,
                fr,
                q=args.q,
                beta_p_star=args.beta_p,
                epsilon=args.epsilon,
                delta0l_percent=args.delta0l_percent,
                alpha=args.alpha,
                snr=args.snr,
                vl1=args.vl1,
                vl2=args.vl2,
                vl3=args.vl3,
                numax_spread_percent=args.numax_spread,
                nmax_spread=args.nmax_spread,
                h0_spread_percent=args.h0_spread,
                hfactor=args.hfactor,
                wfactor=args.wfactor,
                tobs_days=args.tobs,
                cadence_s=args.cadence,
                naverage=args.naverage,
                nrealisation=1,
                doplots=0,
                write_inmodel=1,
                limit_data_range=0,
                do_modelfiles=0,
            )

            out_png = frames_dir / f"frame_{i:04d}.png"
            if args.resume and out_png.exists():
                continue

            seed = args.seed if args.seed_mode == "fixed" else (args.seed + i)
            out_dir = runs_dir / f"run_{i:04d}"
            data_file = run_specsim(
                specsim,
                cfg_path,
                out_dir,
                seed=seed,
                noise_file=args.noise_file,
                cwd=repo_root,
                timeout_s=args.specsim_timeout,
            )
            info_file = out_dir / "Spectra_info" / "0000001.in"
            modes = read_modes_from_info_file(info_file)
            render_frame(
                data_file,
                modes,
                out_png,
                numax=fr.numax_uhz,
                dnu=fr.dnu_uhz,
                nurot_env=fr.nurot_env_uhz,
                nurot_core=fr.nurot_core_uhz,
                dp1=fr.dp1_s,
                max_gamma=fr.max_gamma_uhz,
                zoom_dnu_halfwidth=args.zoom_dnu,
            )

            if (i + 1) % 10 == 0 or (i + 1) == len(frames):
                log.info("Rendered %d/%d frames", i + 1, len(frames))

        partial = (start != 0 or end != len(frames))
        if args.no_encode or partial:
            if partial and not args.no_encode:
                log.warning("Partial frame range rendered; skipping encoding.")
            log.info("Frames written to: %s", frames_dir)
            log.info("(Encoding disabled; use ffmpeg manually if desired)")
            return 0

        encode_mp4(frames_dir, fps=args.fps, out_mp4=out_mp4, timeout_s=args.ffmpeg_timeout)
        log.info("Wrote: %s", out_mp4)
        encoded = True
        return 0
    finally:
        # Keep temporary workdir when frames are the final output (no-encode / partial).
        # Only cleanup if we successfully encoded an mp4.
        if created_temp_workdir and locals().get("encoded", False):
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
