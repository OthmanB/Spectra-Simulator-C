#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def is_comment_or_blank(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#")


def find_next_data_line(lines, start_idx):
    i = start_idx
    while i < len(lines) and is_comment_or_blank(lines[i]):
        i += 1
    return i


def strip_inline_comment(s: str) -> str:
    # Mimic io_star_params.cpp:rem_comments(line, "#") behavior for the main.cfg parsing.
    return s.split("#", 1)[0].strip()


def patch_main_cfg(
    src_cfg: Path,
    repo_root: Path,
    fixed_template: str,
    tobs_days: float,
    cadence_sec: float,
    nspectra: int,
    nrealisation: int,
    disable_plots: bool,
    disable_modelfiles: bool,
):
    """Return patched cfg content as list of lines.

    This produces a minimal, fast-running random configuration, while preserving the model_name.
    """

    text = src_cfg.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines(True)  # keep line endings

    # Identify key lines in the main.cfg format as currently parsed by io_star_params.cpp:read_main_cfg
    i = 0
    forest_idx = find_next_data_line(lines, i)
    if forest_idx >= len(lines):
        raise RuntimeError(f"No forest_type line found in {src_cfg}")

    # forest_type line -> force minimal random
    lines[forest_idx] = "random 1\n"

    model_idx = find_next_data_line(lines, forest_idx + 1)
    if model_idx >= len(lines):
        raise RuntimeError(f"No model_name line found in {src_cfg}")

    model_line = strip_inline_comment(lines[model_idx])
    model_tokens = model_line.split()
    model_name = model_tokens[0]
    extra_tokens = model_tokens[1:]

    # Patch extra_params when it refers to known repo files.
    # Prefer: if an *.in filename is present and exists under Configurations/infiles/, use that.
    new_extra = ""
    if extra_tokens:
        m = None
        for t in extra_tokens:
            mm = re.search(r"([0-9]+)\.in$", t)
            if mm:
                m = mm
                break
        if m:
            infile_name = f"{m.group(1)}.in"
            candidate = repo_root / "Configurations" / "infiles" / infile_name
            if candidate.exists():
                new_extra = str(candidate)
            else:
                # If the candidate doesn't exist, keep the original tokens.
                new_extra = " ".join(extra_tokens)
        else:
            new_extra = " ".join(extra_tokens)

    # Rewrite model line
    if new_extra:
        lines[model_idx] = f"{model_name}   {new_extra}\n"
    else:
        lines[model_idx] = f"{model_name}\n"

    # template selection line
    template_idx = find_next_data_line(lines, model_idx + 1)
    if template_idx >= len(lines):
        raise RuntimeError(f"No template line found in {src_cfg}")
    template_line = strip_inline_comment(lines[template_idx])
    if template_line.lower() in ("all", "*"):
        lines[template_idx] = f"{fixed_template}\n"

    # labels line
    labels_idx = find_next_data_line(lines, template_idx + 1)
    if labels_idx >= len(lines):
        raise RuntimeError(f"No labels line found in {src_cfg}")
    labels = strip_inline_comment(lines[labels_idx]).split()
    nlabels = len(labels)
    if nlabels == 0:
        raise RuntimeError(f"Empty labels line in {src_cfg}")

    # val_min, val_max
    valmin_idx = find_next_data_line(lines, labels_idx + 1)
    valmax_idx = find_next_data_line(lines, valmin_idx + 1)
    if valmax_idx >= len(lines):
        raise RuntimeError(f"Missing val_min/val_max lines in {src_cfg}")

    # step line -> force all constants (random mode requires step in {0,1}; we choose 0)
    step_idx = find_next_data_line(lines, valmax_idx + 1)
    if step_idx >= len(lines):
        raise RuntimeError(f"Missing step line in {src_cfg}")
    lines[step_idx] = " ".join(["0"] * nlabels) + "\n"

    # skip the 'Tobs Cadence ...' label line, then patch numeric line
    obs_label_idx = find_next_data_line(lines, step_idx + 1)
    obs_val_idx = find_next_data_line(lines, obs_label_idx + 1)
    if obs_val_idx >= len(lines):
        raise RuntimeError(f"Missing observation values line in {src_cfg}")
    lines[obs_val_idx] = f"{tobs_days} {cadence_sec} {nspectra} {nrealisation}\n"

    # Ensure the config tail exists and is minimal.
    # read_main_cfg() expects 5 additional lines after the observation line:
    #   erase_old_files, doplots, write_inmodel, limit_data_range, do_modelfiles [modelname]
    erase_old_files = "1"
    doplots = "0" if disable_plots else "1"
    write_inmodel = "0"
    limit_data_range = "0"
    do_modelfiles = "0" if disable_modelfiles else "1"

    lines = lines[: obs_val_idx + 1]
    lines.extend(
        [
            erase_old_files + "\n",
            doplots + "\n",
            write_inmodel + "\n",
            limit_data_range + "\n",
            do_modelfiles + "\n",
        ]
    )
    return lines


def validate_outputs(out_dir: Path):
    combi = out_dir / "Combinations.txt"
    if not combi.exists():
        raise RuntimeError(f"Missing {combi}")

    for sub in ("Spectra_ascii", "Spectra_info", "Spectra_modelfile", "Spectra_plot"):
        p = out_dir / sub
        if not p.exists() or not p.is_dir():
            raise RuntimeError(f"Missing output subdirectory {p}")

    # At least one spectrum and one info file
    spectra_files = list((out_dir / "Spectra_ascii").glob("*"))
    info_files = list((out_dir / "Spectra_info").glob("*.in"))
    if not spectra_files:
        raise RuntimeError(f"No spectra files produced in {out_dir / 'Spectra_ascii'}")
    if not info_files:
        raise RuntimeError(f"No .in info files produced in {out_dir / 'Spectra_info'}")


def run_one(
    specsim: Path, cfg_path: Path, out_dir: Path, main_dir: Path, noise_cfg: Path
):
    cmd = [
        str(specsim),
        "--main_file",
        str(cfg_path),
        "--main_dir",
        str(main_dir),
        "--noise_file",
        str(noise_cfg),
        "--out_dir",
        str(out_dir),
        "--force-create-output-dir",
        "1",
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, p.stdout


def main():
    ap = argparse.ArgumentParser(
        description="Run minimal smoke tests for non-obsolete example configs"
    )
    ap.add_argument(
        "--specsim",
        default="build/specsim",
        help="Path to specsim binary (default: build/specsim)",
    )
    ap.add_argument(
        "--examples-dir",
        default="Configurations/examples_cfg",
        help="Examples directory (default: Configurations/examples_cfg)",
    )
    ap.add_argument(
        "--noise-cfg",
        default="Configurations/noise_Kallinger2014.cfg",
        help="Noise configuration file path (default: Configurations/noise_Kallinger2014.cfg)",
    )
    ap.add_argument(
        "--fixed-template",
        default="12508433.template",
        help="Template used when example requests 'all' (default: 12508433.template)",
    )
    ap.add_argument(
        "--tobs", type=float, default=30.0, help="Tobs (days) for smoke cfgs"
    )
    ap.add_argument(
        "--cadence", type=float, default=120.0, help="Cadence (sec) for smoke cfgs"
    )
    ap.add_argument(
        "--nspectra", type=int, default=1, help="Nspectra (averaging) for smoke cfgs"
    )
    ap.add_argument(
        "--nrealisation", type=int, default=1, help="Nrealisation for smoke cfgs"
    )
    ap.add_argument(
        "--out-root",
        default=None,
        help="Root output directory. If omitted, a temp dir is used.",
    )
    ap.add_argument("--keep", action="store_true", help="Keep temporary files")
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    specsim = (
        (repo_root / args.specsim).resolve()
        if not os.path.isabs(args.specsim)
        else Path(args.specsim)
    )
    examples_dir = (
        (repo_root / args.examples_dir).resolve()
        if not os.path.isabs(args.examples_dir)
        else Path(args.examples_dir)
    )
    noise_cfg = (
        (repo_root / args.noise_cfg).resolve()
        if not os.path.isabs(args.noise_cfg)
        else Path(args.noise_cfg)
    )

    if not specsim.exists():
        print(f"ERROR: specsim not found at {specsim}", file=sys.stderr)
        sys.exit(2)
    if not examples_dir.exists():
        print(f"ERROR: examples dir not found: {examples_dir}", file=sys.stderr)
        sys.exit(2)
    if not noise_cfg.exists():
        print(f"ERROR: noise cfg not found: {noise_cfg}", file=sys.stderr)
        sys.exit(2)

    # Discover example cfgs (exclude obsolete directory)
    cfgs = sorted([p for p in examples_dir.glob("main.cfg*") if p.is_file()])
    cfgs = [p for p in cfgs if "obselete" not in str(p)]
    if not cfgs:
        print(f"ERROR: no cfg files found under {examples_dir}", file=sys.stderr)
        sys.exit(2)

    # Output root
    if args.out_root:
        out_root = Path(args.out_root).resolve()
        out_root.mkdir(parents=True, exist_ok=True)
        tmp_ctx = None
    else:
        tmp_ctx = tempfile.TemporaryDirectory(prefix="specsim-smoke-")
        out_root = Path(tmp_ctx.name)

    work_dir = out_root / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    failures = 0
    for cfg in cfgs:
        model_out = out_root / ("out-" + cfg.name.replace("/", "_"))
        model_out.mkdir(parents=True, exist_ok=True)

        patched_lines = patch_main_cfg(
            cfg,
            repo_root=repo_root,
            fixed_template=args.fixed_template,
            tobs_days=args.tobs,
            cadence_sec=args.cadence,
            nspectra=args.nspectra,
            nrealisation=args.nrealisation,
            disable_plots=True,
            disable_modelfiles=True,
        )

        patched_cfg = work_dir / (cfg.name + ".smoke")
        patched_cfg.write_text("".join(patched_lines), encoding="utf-8")

        rc, out = run_one(specsim, patched_cfg, model_out, examples_dir, noise_cfg)
        if rc != 0:
            failures += 1
            print("=" * 80)
            print(f"FAIL: {cfg} (rc={rc})")
            print(f"patched cfg: {patched_cfg}")
            print(f"out dir   : {model_out}")
            print("--- output ---")
            print(out)
            continue

        try:
            validate_outputs(model_out)
        except Exception as e:
            failures += 1
            print("=" * 80)
            print(f"FAIL: {cfg} (output validation)")
            print(f"patched cfg: {patched_cfg}")
            print(f"out dir   : {model_out}")
            print(f"error     : {e}")
            print("--- output ---")
            print(out)
            continue

        print(f"PASS: {cfg.name} -> {model_out}")

    if failures:
        print(f"\nSmoke tests: {failures} failure(s)")
        if args.keep or args.out_root:
            print(f"Outputs kept under: {out_root}")
        sys.exit(1)

    print("\nSmoke tests: all passed")
    if args.keep or args.out_root:
        print(f"Outputs kept under: {out_root}")
    else:
        # temp dir will be cleaned by context manager
        pass


if __name__ == "__main__":
    main()
