#!/usr/bin/env python3

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from smoke_test import patch_main_cfg, find_next_data_line, strip_inline_comment  # noqa: E402


EXAMPLE_FILES = [
    "main.cfg.Alm",
    "main.cfg.aj",
    "main.cfg.aj_GRANscaledKallinger",
    "main.cfg.freeDP_curvepmodes.v3_GRANscaled",
]

ZERO_SPREAD_LABELS = {
    "H_spread",
    "nu_spread",
    "Gamma_spread",
    "numax_spread",
    "nmax_spread",
    "H0_spread",
}


def zero_spread_params(lines):
    idx = 0
    forest_idx = find_next_data_line(lines, idx)
    model_idx = find_next_data_line(lines, forest_idx + 1)
    template_idx = find_next_data_line(lines, model_idx + 1)
    labels_idx = find_next_data_line(lines, template_idx + 1)
    valmin_idx = find_next_data_line(lines, labels_idx + 1)
    valmax_idx = find_next_data_line(lines, valmin_idx + 1)

    labels = strip_inline_comment(lines[labels_idx]).split()
    val_min = strip_inline_comment(lines[valmin_idx]).split()
    val_max = strip_inline_comment(lines[valmax_idx]).split()
    if len(labels) != len(val_min) or len(labels) != len(val_max):
        return lines

    for i, label in enumerate(labels):
        if label in ZERO_SPREAD_LABELS:
            val_min[i] = "0"
            val_max[i] = "0"

    lines[valmin_idx] = " ".join(val_min) + "\n"
    lines[valmax_idx] = " ".join(val_max) + "\n"
    return lines


def find_repo_root(specsim_path: Path) -> Path:
    for parent in [specsim_path.parent, specsim_path.parent.parent, specsim_path.parent.parent.parent]:
        if (parent / "Configurations" / "examples_cfg").exists():
            return parent
    raise RuntimeError(f"Could not locate repo root for {specsim_path}")


def run_specsim(specsim: Path, cwd: Path, cfg_path: Path, out_dir: Path, noise_cfg: Path, seed: int):
    tmp_dir = cwd / "Configurations" / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(specsim),
        "--main_file",
        str(cfg_path),
        "--main_dir",
        str(cwd),
        "--noise_file",
        str(noise_cfg),
        "--out_dir",
        str(out_dir),
        "--force-create-output-dir",
        "1",
        "--seed",
        str(seed),
    ]
    p = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return p.returncode, p.stdout


def read_text_lines(path: Path):
    return [line.rstrip() for line in path.read_text(encoding="utf-8", errors="replace").splitlines()]


def compare_text_files(path_a: Path, path_b: Path):
    lines_a = [line.strip() for line in read_text_lines(path_a) if line.strip()]
    lines_b = [line.strip() for line in read_text_lines(path_b) if line.strip()]
    return lines_a == lines_b


def find_combinations_file(out_dir: Path):
    candidate = out_dir / "Combinations.txt"
    if candidate.exists():
        return candidate
    legacy = Path(str(out_dir) + "Combinations.txt")
    if legacy.exists():
        return legacy
    return None


def extract_mode_rows(model_path: Path):
    lines = read_text_lines(model_path)
    start_idx = None
    end_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("# Eigen solution input parameters"):
            start_idx = idx + 1
        if start_idx is not None and line.startswith("# Noise parameters"):
            end_idx = idx
            break
    if start_idx is None:
        raise RuntimeError(f"Missing mode block header in {model_path}")
    if end_idx is None:
        end_idx = len(lines)
    rows = []
    for line in lines[start_idx:end_idx]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 6:
            raise RuntimeError(f"Unexpected mode row in {model_path}: {line}")
        rows.append([float(x) for x in parts[:6]])
    return rows


def is_close(a: float, b: float, rtol: float, atol: float) -> bool:
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def compare_mode_rows(rows_a, rows_b, rtol: float, atol: float, exact: bool):
    if len(rows_a) != len(rows_b):
        return False, [f"Row count differs: {len(rows_a)} vs {len(rows_b)}"]
    diffs = []
    for i, (ra, rb) in enumerate(zip(rows_a, rows_b)):
        if len(ra) != len(rb):
            diffs.append(f"Row {i} length differs: {len(ra)} vs {len(rb)}")
            continue
        if exact:
            if ra != rb:
                diffs.append(f"Row {i} differs (exact): {ra} vs {rb}")
            continue
        for j, (a, b) in enumerate(zip(ra, rb)):
            if not is_close(a, b, rtol, atol):
                diffs.append(f"Row {i} col {j} differs: {a} vs {b}")
                if len(diffs) >= 10:
                    return False, diffs
    return len(diffs) == 0, diffs


def compare_model_files(dir_a: Path, dir_b: Path, rtol: float, atol: float, exact: bool):
    files_a = sorted(dir_a.glob("*.model"))
    files_b = sorted(dir_b.glob("*.model"))
    if len(files_a) != len(files_b):
        return False, [f"Model file count differs: {len(files_a)} vs {len(files_b)}"]
    diffs = []
    for fa, fb in zip(files_a, files_b):
        rows_a = extract_mode_rows(fa)
        rows_b = extract_mode_rows(fb)
        ok, local_diffs = compare_mode_rows(rows_a, rows_b, rtol, atol, exact)
        if not ok:
            diffs.append(f"Differences in {fa.name} vs {fb.name}:")
            diffs.extend(local_diffs)
            if len(diffs) >= 10:
                break
    return len(diffs) == 0, diffs


def main():
    ap = argparse.ArgumentParser(
        description="Compare model-only outputs across two specsim binaries"
    )
    ap.add_argument("--specsim-a", required=True, help="Path to baseline specsim")
    ap.add_argument("--specsim-b", required=True, help="Path to candidate specsim")
    ap.add_argument(
        "--examples-dir",
        default="Configurations/examples_cfg",
        help="Examples directory (default: Configurations/examples_cfg)",
    )
    ap.add_argument(
        "--noise-cfg",
        default="Configurations/noise_Kallinger2014.cfg",
        help="Noise cfg path (default: Configurations/noise_Kallinger2014.cfg)",
    )
    ap.add_argument(
        "--fixed-template",
        default="12508433.template",
        help="Template used when cfg requests 'all' (default: 12508433.template)",
    )
    ap.add_argument("--tobs", type=float, default=30.0, help="Tobs (days)")
    ap.add_argument("--cadence", type=float, default=120.0, help="Cadence (sec)")
    ap.add_argument("--nspectra", type=int, default=1, help="Nspectra")
    ap.add_argument("--nrealisation", type=int, default=1, help="Nrealisation")
    ap.add_argument("--seed", type=int, default=123, help="Seed for RNG")
    ap.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance")
    ap.add_argument("--atol", type=float, default=1e-8, help="Absolute tolerance")
    ap.add_argument("--exact", action="store_true", help="Exact compare (no tolerance)")
    ap.add_argument(
        "--out-root",
        default=None,
        help="Output root directory (default: temp)",
    )
    ap.add_argument("--keep", action="store_true", help="Keep outputs")
    ap.add_argument(
        "--example",
        action="append",
        default=None,
        help="Limit to specific example cfg (repeatable)",
    )
    args = ap.parse_args()

    specsim_a = Path(args.specsim_a).resolve()
    specsim_b = Path(args.specsim_b).resolve()
    if not specsim_a.exists():
        raise SystemExit(f"specsim-a not found: {specsim_a}")
    if not specsim_b.exists():
        raise SystemExit(f"specsim-b not found: {specsim_b}")

    repo_a = find_repo_root(specsim_a)
    repo_b = find_repo_root(specsim_b)

    examples_dir = (
        (REPO_ROOT / args.examples_dir).resolve()
        if not os.path.isabs(args.examples_dir)
        else Path(args.examples_dir)
    )
    noise_cfg = (
        (REPO_ROOT / args.noise_cfg).resolve()
        if not os.path.isabs(args.noise_cfg)
        else Path(args.noise_cfg)
    )
    if not examples_dir.exists():
        raise SystemExit(f"Examples dir not found: {examples_dir}")
    if not noise_cfg.exists():
        raise SystemExit(f"Noise cfg not found: {noise_cfg}")

    selected = args.example if args.example else EXAMPLE_FILES
    selected_paths = [examples_dir / name for name in selected]
    for p in selected_paths:
        if not p.exists():
            raise SystemExit(f"Example cfg not found: {p}")

    if args.out_root:
        out_root = Path(args.out_root).resolve()
        out_root.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        out_root = Path(tempfile.mkdtemp(prefix="specsim-compare-"))
        cleanup = not args.keep

    failures = 0
    for cfg_path in selected_paths:
        name = cfg_path.name
        example_root = out_root / name
        example_root.mkdir(parents=True, exist_ok=True)

        patched = patch_main_cfg(
            cfg_path,
            repo_root=REPO_ROOT,
            fixed_template=args.fixed_template,
            tobs_days=args.tobs,
            cadence_sec=args.cadence,
            nspectra=args.nspectra,
            nrealisation=args.nrealisation,
            disable_plots=True,
            disable_modelfiles=False,
        )
        patched = zero_spread_params(patched)
        patched_cfg = example_root / "main.cfg"
        patched_cfg.write_text("".join(patched), encoding="utf-8")

        out_a = example_root / "out_a"
        out_b = example_root / "out_b"
        rc_a, log_a = run_specsim(specsim_a, repo_a, patched_cfg, out_a, noise_cfg, args.seed)
        rc_b, log_b = run_specsim(specsim_b, repo_b, patched_cfg, out_b, noise_cfg, args.seed)

        if rc_a != 0 or rc_b != 0:
            failures += 1
            print(f"[{name}] run failure:")
            if rc_a != 0:
                print("  specsim-a failed:")
                print(log_a)
            if rc_b != 0:
                print("  specsim-b failed:")
                print(log_b)
            continue

        combi_a = find_combinations_file(out_a)
        combi_b = find_combinations_file(out_b)
        if combi_a is None or combi_b is None:
            failures += 1
            print(f"[{name}] missing Combinations.txt")
            continue

        if not compare_text_files(combi_a, combi_b):
            failures += 1
            print(f"[{name}] Combinations.txt differs")
            continue

        model_dir_a = out_a / "Spectra_modelfile"
        model_dir_b = out_b / "Spectra_modelfile"
        if not model_dir_a.exists() or not model_dir_b.exists():
            failures += 1
            print(f"[{name}] missing Spectra_modelfile output")
            continue

        ok, diffs = compare_model_files(model_dir_a, model_dir_b, args.rtol, args.atol, args.exact)
        if not ok:
            failures += 1
            print(f"[{name}] model differences detected:")
            for diff in diffs[:10]:
                print(f"  {diff}")
            continue

        print(f"[{name}] OK")

    if cleanup:
        shutil.rmtree(out_root)
    else:
        print(f"Outputs kept under: {out_root}")

    if failures:
        raise SystemExit(f"Comparison failed for {failures} example(s)")


if __name__ == "__main__":
    main()
