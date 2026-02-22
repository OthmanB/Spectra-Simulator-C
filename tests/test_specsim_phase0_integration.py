import os
import hashlib
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


def run(cmd, cwd=None, timeout=None):
    p = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
    )
    return p.returncode, p.stdout


def parse_spectrum_file(path: Path):
    freqs = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            freqs.append(float(parts[0]))
    return freqs


def hash_file(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_output_hashes(out_dir: Path):
    hashes = {}
    for rel in (Path("Combinations.txt"),):
        p = out_dir / rel
        if p.exists():
            hashes[str(rel)] = hash_file(p)

    for sub in ("Spectra_ascii", "Spectra_info"):
        for p in sorted((out_dir / sub).glob("*")):
            if p.is_file():
                rel = str(Path(sub) / p.name)
                hashes[rel] = hash_file(p)
    return hashes


def parse_mode_widths_from_info_file(path: Path):
    """Parse widths (W) from Spectra_info/*.in produced by specsim."""

    widths = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("ID="):
                continue
            parts = s.split()
            if len(parts) < 4:
                continue
            try:
                l = int(float(parts[0]))
            except Exception:
                continue
            if l not in (0, 1, 2, 3):
                continue
            try:
                w = float(parts[3])
            except Exception:
                continue
            widths.append(w)
    return widths


def rewrite_infile_widths_constant(src_in: Path, dst_in: Path, width_value: float):
    """Create a modified *.in file with constant width for all mode lines."""

    lines = src_in.read_text(encoding="utf-8", errors="replace").splitlines(True)
    out = []
    i = 0

    # Copy initial header comments
    while i < len(lines) and lines[i].lstrip().startswith("#"):
        out.append(lines[i])
        i += 1
    if i >= len(lines):
        raise RuntimeError(f"Unexpected format: only comments in {src_in}")

    # Identifier line
    out.append(lines[i])
    i += 1

    # Skip/copy any comment lines following the identifier
    while i < len(lines) and lines[i].lstrip().startswith("#"):
        out.append(lines[i])
        i += 1

    # Optional Tobs/Cadence line: 2 tokens
    if i < len(lines):
        toks = lines[i].split()
        if len(toks) == 2:
            out.append(lines[i])
            i += 1
            while i < len(lines) and lines[i].lstrip().startswith("#"):
                out.append(lines[i])
                i += 1

    # Mode table until the next comment line starts the noise section
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith("#"):
            break
        if not line.strip():
            out.append(line)
            i += 1
            continue
        parts = line.split()
        if len(parts) >= 4:
            parts[3] = f"{width_value}"
            out.append(" ".join(parts) + "\n")
        else:
            out.append(line)
        i += 1

    # Remainder (noise section)
    out.extend(lines[i:])
    dst_in.write_text("".join(out), encoding="utf-8")


class TestSpecsimPhase0Integration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[1]
        cls.build_dir = cls.repo_root / "build-phase0-tests"
        cls.build_dir.mkdir(parents=True, exist_ok=True)

        # Build specsim once for the test suite.
        rc, out = run(["cmake", "-S", str(cls.repo_root), "-B", str(cls.build_dir)])
        if rc != 0:
            raise RuntimeError("CMake configure failed:\n" + out)
        rc, out = run(["cmake", "--build", str(cls.build_dir), "-j"])
        if rc != 0:
            raise RuntimeError("CMake build failed:\n" + out)

        cls.specsim = cls.build_dir / "specsim"
        if not cls.specsim.exists():
            raise RuntimeError(f"specsim binary not found at {cls.specsim}")

        cls.smoke_mod = __import__("scripts.smoke_test", fromlist=["patch_main_cfg"])

    def _make_sandbox(self):
        td = tempfile.TemporaryDirectory(prefix="specsim-phase0-it-")
        root = Path(td.name)

        # minimal runtime tree under the sandbox
        (root / "Configurations").mkdir(parents=True, exist_ok=True)

        # templates are required for evolved-star configs
        (root / "Configurations" / "templates").symlink_to(
            self.repo_root / "Configurations" / "templates"
        )

        # noise cfg may be referenced by CLI
        shutil.copy2(
            self.repo_root / "Configurations" / "noise_Kallinger2014.cfg",
            root / "Configurations" / "noise_Kallinger2014.cfg",
        )

        # external directory is used for some mixed-mode metadata files (writes are non-critical)
        (root / "external").symlink_to(self.repo_root / "external")

        # ensure tmp dir does NOT exist initially (this is what Phase 0 fixes)
        tmp_dir = root / "Configurations" / "tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)

        return td, root

    def _make_sandbox_with_templates(self, templates):
        td = tempfile.TemporaryDirectory(prefix="specsim-phase0-it-templates-")
        root = Path(td.name)

        (root / "Configurations").mkdir(parents=True, exist_ok=True)
        templates_dir = root / "Configurations" / "templates"
        templates_dir.mkdir(parents=True, exist_ok=True)

        for item in templates:
            if isinstance(item, tuple):
                name, content = item
                (templates_dir / name).write_text(content, encoding="utf-8")
            else:
                src = Path(item)
                shutil.copy2(src, templates_dir / src.name)

        shutil.copy2(
            self.repo_root / "Configurations" / "noise_Kallinger2014.cfg",
            root / "Configurations" / "noise_Kallinger2014.cfg",
        )

        (root / "external").symlink_to(self.repo_root / "external")

        tmp_dir = root / "Configurations" / "tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)

        return td, root

    def _run_cfg_in_sandbox(
        self,
        example_cfg: Path,
        fixed_template: str,
        tobs_days: float = 2.0,
        cadence_sec: float = 120.0,
        extra_args=None,
    ):
        patch_main_cfg = getattr(self.smoke_mod, "patch_main_cfg")

        td, root = self._make_sandbox()
        out_dir = root / "out"

        patched_lines = patch_main_cfg(
            example_cfg,
            repo_root=self.repo_root,
            fixed_template=fixed_template,
            tobs_days=tobs_days,
            cadence_sec=cadence_sec,
            nspectra=1,
            nrealisation=1,
            disable_plots=True,
            disable_modelfiles=True,
        )
        patched_cfg = root / "main.cfg"
        patched_cfg.write_text("".join(patched_lines), encoding="utf-8")

        cmd = [
            str(self.specsim),
            "--main_file",
            str(patched_cfg),
            "--noise_file",
            str(root / "Configurations" / "noise_Kallinger2014.cfg"),
            "--out_dir",
            str(out_dir),
            "--force-create-output-dir",
            "1",
        ]
        if extra_args:
            cmd.extend(extra_args)
        rc, out = run(cmd, cwd=root, timeout=600)
        return td, root, out_dir, patched_cfg, rc, out

    def _run_cfg_text_in_sandbox(self, cfg_text: str, extra_args=None):
        td, root = self._make_sandbox()
        out_dir = root / "out"
        cfg_path = root / "main.cfg"
        cfg_path.write_text(cfg_text, encoding="utf-8")

        cmd = [
            str(self.specsim),
            "--main_file",
            str(cfg_path),
            "--noise_file",
            str(root / "Configurations" / "noise_Kallinger2014.cfg"),
            "--out_dir",
            str(out_dir),
            "--force-create-output-dir",
            "1",
        ]
        if extra_args:
            cmd.extend(extra_args)
        rc, out = run(cmd, cwd=root, timeout=600)
        return td, root, out_dir, cfg_path, rc, out

    def test_tmp_dir_created_and_combinations_path(self):
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg, fixed_template="12508433.template"
        )
        try:
            self.assertEqual(rc, 0, msg=out)

            # tmp dir created
            self.assertTrue((root / "Configurations" / "tmp").is_dir())

            # Combinations.txt is inside out_dir
            self.assertTrue((out_dir / "Combinations.txt").exists())

            # Old bug path must not exist
            bug_path = Path(str(out_dir) + "Combinations.txt")
            if bug_path != (out_dir / "Combinations.txt"):
                self.assertFalse(bug_path.exists())

            # Output directories exist
            for sub in (
                "Spectra_ascii",
                "Spectra_info",
                "Spectra_modelfile",
                "Spectra_plot",
            ):
                self.assertTrue((out_dir / sub).is_dir())

            # At least one spectrum file and one info file
            spectra = list((out_dir / "Spectra_ascii").glob("*.data"))
            infos = list((out_dir / "Spectra_info").glob("*.in"))
            self.assertTrue(spectra)
            self.assertTrue(infos)

        finally:
            td.cleanup()

    def test_frequency_grid_invariants(self):
        # Use a simple MS model for deterministic frequency axis properties.
        tobs_days = 2.0
        cadence_sec = 120.0
        Delta = 1e6 / cadence_sec / 2.0
        df = 1e6 / (tobs_days * 86400.0)
        expected_n = int(Delta / df)

        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg,
            fixed_template="12508433.template",
            tobs_days=tobs_days,
            cadence_sec=cadence_sec,
        )
        try:
            self.assertEqual(rc, 0, msg=out)
            spectra = sorted((out_dir / "Spectra_ascii").glob("*.data"))
            self.assertTrue(spectra)

            freqs = parse_spectrum_file(spectra[0])
            self.assertTrue(freqs)

            # Invariants
            self.assertAlmostEqual(freqs[0], 0.0, places=6)
            self.assertAlmostEqual(freqs[-1], Delta, places=4)
            self.assertTrue(all(freqs[i] < freqs[i + 1] for i in range(len(freqs) - 1)))

            # Size is as per implementation (integer truncation)
            self.assertEqual(len(freqs), expected_n)

            # Spacing is (approximately) constant
            steps = [freqs[i + 1] - freqs[i] for i in range(len(freqs) - 1)]
            self.assertGreater(len(steps), 10)
            self.assertLess(max(steps) - min(steps), 1e-3)
        finally:
            td.cleanup()

    def test_evolved_model_smoke(self):
        # Evolved-star pipeline: requires templates and Kallinger noise aggregation.
        example_cfg = (
            self.repo_root
            / "Configurations"
            / "examples_cfg"
            / "main.cfg.freeDP_curvepmodes.v3_GRANscaled"
        )
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg,
            fixed_template="12508433.template",
            tobs_days=2.0,
            cadence_sec=120.0,
        )
        try:
            self.assertEqual(rc, 0, msg=out)
            self.assertTrue((out_dir / "Combinations.txt").exists())
            spectra = list((out_dir / "Spectra_ascii").glob("*.data"))
            infos = list((out_dir / "Spectra_info").glob("*.in"))
            self.assertTrue(spectra)
            self.assertTrue(infos)
        finally:
            td.cleanup()

    def test_gamma_spread_applies_per_mode(self):
        # Regression for generate_cfg_from_synthese_file_Wscaled_aj Gamma_spread:
        # ensure width jitter is applied to all modes (not just one element).

        src_in = self.repo_root / "Configurations" / "infiles" / "8379927.in"
        self.assertTrue(src_in.exists())

        td, root = self._make_sandbox()
        try:
            ref_in = root / "ref_constant_width.in"
            rewrite_infile_widths_constant(src_in, ref_in, width_value=1.0)

            def run_cfg(gamma_spread: float, out_subdir: str):
                out_dir = root / out_subdir
                cfg_text = "".join(
                    [
                        "random 1\n",
                        f"generate_cfg_from_synthese_file_Wscaled_aj {ref_in}\n",
                        "NONE\n",
                        "Dnu epsilon delta0l_percent HNR a1ovGamma Gamma_at_numax a2 a3 a4 a5 a6 beta_asym i H_spread nu_spread Gamma_spread do_flat_noise\n",
                        f"70 0.5 1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 {gamma_spread} 0\n",
                        f"70 0.5 1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 {gamma_spread} 0\n",
                        "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n",
                        "Tobs Cadence Naverage Nrealisation\n",
                        "2 120 1 1\n",
                        "1\n",
                        "0\n",
                        "0\n",
                        "0\n",
                        "0\n",
                    ]
                )
                cfg_path = root / "main.cfg"
                cfg_path.write_text(cfg_text, encoding="utf-8")

                cmd = [
                    str(self.specsim),
                    "--main_file",
                    str(cfg_path),
                    "--noise_file",
                    str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                    "--out_dir",
                    str(out_dir),
                    "--force-create-output-dir",
                    "1",
                ]
                rc, out = run(cmd, cwd=root, timeout=600)
                self.assertEqual(rc, 0, msg=out)

                info_files = list((out_dir / "Spectra_info").glob("*.in"))
                self.assertTrue(info_files)
                widths = parse_mode_widths_from_info_file(info_files[0])
                self.assertTrue(widths)
                self.assertTrue(all(w > 0 for w in widths))
                return widths

            widths0 = run_cfg(0.0, "out-gamma0")
            unique0 = set(round(w, 10) for w in widths0)
            self.assertEqual(len(unique0), 1)

            widths = run_cfg(20.0, "out-gamma20")
            # With constant reference widths and Gamma_spread>0, most widths should differ.
            unique = set(round(w, 10) for w in widths)
            self.assertGreater(len(unique), int(len(widths) * 0.5))
        finally:
            td.cleanup()

    def test_validation_rejects_mismatched_vector_sizes(self):
        # labels has 2 entries but val_min has 1 -> must fail fast in read_main_cfg
        infile = self.repo_root / "Configurations" / "infiles" / "8379927.in"
        self.assertTrue(infile.exists())

        cfg_text = "".join(
            [
                "random 1\n",
                f"generate_cfg_from_synthese_file_Wscaled_aj {infile}\n",
                "NONE\n",
                "Dnu epsilon\n",
                "70\n",
                "70 0.5\n",
                "0 0\n",
                "Tobs Cadence Naverage Nrealisation\n",
                "2 120 1 1\n",
                "1\n",
                "0\n",
                "0\n",
                "0\n",
                "0\n",
            ]
        )
        td, root, out_dir, cfg_path, rc, out = self._run_cfg_text_in_sandbox(cfg_text)
        try:
            self.assertNotEqual(rc, 0)
            self.assertIn("different size", out)
        finally:
            td.cleanup()

    def test_validation_rejects_duplicate_labels(self):
        # Create a valid minimal config, then duplicate a label.
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg, fixed_template="12508433.template"
        )
        try:
            self.assertEqual(rc, 0, msg=out)

            # mutate labels line: replace 'epsilon' with 'Dnu' to create duplicates
            lines = patched_cfg.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines(True)
            idx = None
            for i, line in enumerate(lines):
                if "delta0l_percent" in line and not line.lstrip().startswith("#"):
                    idx = i
                    break
            self.assertIsNotNone(idx)

            lhs, sep, comment = lines[idx].partition("#")
            toks = lhs.split()
            self.assertGreaterEqual(len(toks), 2)
            toks[1] = "Dnu"
            new_lhs = " ".join(toks)
            lines[idx] = new_lhs + (" " + sep + comment if sep else "\n")
            patched_cfg.write_text("".join(lines), encoding="utf-8")

            # rerun
            cmd = [
                str(self.specsim),
                "--main_file",
                str(patched_cfg),
                "--noise_file",
                str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                "--out_dir",
                str(out_dir),
                "--force-create-output-dir",
                "1",
            ]
            rc2, out2 = run(cmd, cwd=root, timeout=600)
            self.assertNotEqual(rc2, 0)
            self.assertIn("Duplicate labels", out2)
            self.assertIn("Dnu", out2)
        finally:
            td.cleanup()

    def test_validation_rejects_missing_or_extra_labels(self):
        # Create a valid minimal config, then introduce a typo.
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg, fixed_template="12508433.template"
        )
        try:
            self.assertEqual(rc, 0, msg=out)

            lines = patched_cfg.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines(True)
            idx = None
            for i, line in enumerate(lines):
                if "delta0l_percent" in line and not line.lstrip().startswith("#"):
                    idx = i
                    break
            self.assertIsNotNone(idx)

            lhs, sep, comment = lines[idx].partition("#")
            toks = lhs.split()
            toks = [
                "delta0l_percent_typo" if t == "delta0l_percent" else t for t in toks
            ]
            new_lhs = " ".join(toks)
            lines[idx] = new_lhs + (" " + sep + comment if sep else "\n")
            patched_cfg.write_text("".join(lines), encoding="utf-8")

            cmd = [
                str(self.specsim),
                "--main_file",
                str(patched_cfg),
                "--noise_file",
                str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                "--out_dir",
                str(out_dir),
                "--force-create-output-dir",
                "1",
            ]
            rc2, out2 = run(cmd, cwd=root, timeout=600)
            self.assertNotEqual(rc2, 0)
            self.assertIn("Missing labels", out2)
            self.assertIn("delta0l_percent", out2)
        finally:
            td.cleanup()

    def test_negative_delta0l_percent_warns(self):
        infile = self.repo_root / "Configurations" / "infiles" / "8379927.in"
        self.assertTrue(infile.exists())

        cfg_text = "".join(
            [
                "random 1\n",
                f"generate_cfg_from_synthese_file_Wscaled_aj {infile}\n",
                "NONE\n",
                "Dnu epsilon delta0l_percent HNR a1ovGamma Gamma_at_numax a2 a3 a4 a5 a6 beta_asym i H_spread nu_spread Gamma_spread do_flat_noise\n",
                "70 0.5 -1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 0 0\n",
                "70 0.5 -1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 0 0\n",
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n",
                "Tobs Cadence Naverage Nrealisation\n",
                "2 120 1 1\n",
                "1\n",
                "0\n",
                "0\n",
                "0\n",
                "0\n",
            ]
        )
        td, root, out_dir, cfg_path, rc, out = self._run_cfg_text_in_sandbox(cfg_text)
        try:
            self.assertEqual(rc, 0, msg=out)
            self.assertIn("delta0l_percent is negative", out)
        finally:
            td.cleanup()

    def test_seed_reproducibility(self):
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td1, root1, out_dir1, patched_cfg1, rc1, out1 = self._run_cfg_in_sandbox(
            example_cfg,
            fixed_template="12508433.template",
            extra_args=["--seed", "123"],
        )
        try:
            self.assertEqual(rc1, 0, msg=out1)
            hashes1 = collect_output_hashes(out_dir1)
        finally:
            td1.cleanup()

        td2, root2, out_dir2, patched_cfg2, rc2, out2 = self._run_cfg_in_sandbox(
            example_cfg,
            fixed_template="12508433.template",
            extra_args=["--seed", "123"],
        )
        try:
            self.assertEqual(rc2, 0, msg=out2)
            hashes2 = collect_output_hashes(out_dir2)
        finally:
            td2.cleanup()

        self.assertEqual(hashes1, hashes2)

    def test_out_dir_trailing_slash(self):
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        patch_main_cfg = getattr(self.smoke_mod, "patch_main_cfg")

        td, root = self._make_sandbox()
        try:
            out_dir = root / "out"
            out_dir_str = str(out_dir) + "/"
            patched_lines = patch_main_cfg(
                example_cfg,
                repo_root=self.repo_root,
                fixed_template="12508433.template",
                tobs_days=2.0,
                cadence_sec=120.0,
                nspectra=1,
                nrealisation=1,
                disable_plots=True,
                disable_modelfiles=True,
            )
            patched_cfg = root / "main.cfg"
            patched_cfg.write_text("".join(patched_lines), encoding="utf-8")

            cmd = [
                str(self.specsim),
                "--main_file",
                str(patched_cfg),
                "--noise_file",
                str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                "--out_dir",
                out_dir_str,
                "--force-create-output-dir",
                "1",
            ]
            rc, out = run(cmd, cwd=root, timeout=600)
            self.assertEqual(rc, 0, msg=out)
            self.assertTrue((out_dir / "Combinations.txt").exists())
            self.assertTrue((out_dir / "Spectra_ascii").is_dir())
        finally:
            td.cleanup()

    def test_invalid_log_level(self):
        example_cfg = self.repo_root / "Configurations" / "examples_cfg" / "main.cfg.aj"
        td, root, out_dir, patched_cfg, rc, out = self._run_cfg_in_sandbox(
            example_cfg,
            fixed_template="12508433.template",
            extra_args=["--log-level", "bogus"],
        )
        try:
            self.assertNotEqual(rc, 0)
            self.assertIn("Invalid --log-level", out)
        finally:
            td.cleanup()

    def test_template_validation_skips_invalid(self):
        valid_template = (
            self.repo_root / "Configurations" / "templates" / "12508433.template"
        )
        self.assertTrue(valid_template.exists())

        invalid_template = (
            "bad.template",
            "".join(
                [
                    "# invalid template for testing\n",
                    "ID_ref= bad\n",
                    "Dnu_ref= 10\n",
                    "epsilon_ref= 1\n",
                    "numax_ref= 100\n",
                    "# Frequency Height Width\n",
                    "100.0 0.1\n",
                ]
            ),
        )

        td, root = self._make_sandbox_with_templates([valid_template, invalid_template])
        try:
            example_cfg = (
                self.repo_root
                / "Configurations"
                / "examples_cfg"
                / "main.cfg.freeDP_curvepmodes.v3_GRANscaled"
            )
            patch_main_cfg = getattr(self.smoke_mod, "patch_main_cfg")
            patched_lines = patch_main_cfg(
                example_cfg,
                repo_root=self.repo_root,
                fixed_template="all",
                tobs_days=2.0,
                cadence_sec=120.0,
                nspectra=1,
                nrealisation=1,
                disable_plots=True,
                disable_modelfiles=True,
            )
            patched_cfg = root / "main.cfg"
            patched_cfg.write_text("".join(patched_lines), encoding="utf-8")

            cmd = [
                str(self.specsim),
                "--main_file",
                str(patched_cfg),
                "--noise_file",
                str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                "--out_dir",
                str(root / "out"),
                "--force-create-output-dir",
                "1",
            ]
            rc, out = run(cmd, cwd=root, timeout=600)
            self.assertEqual(rc, 0, msg=out)
            self.assertIn("Skipping invalid template file", out)
            self.assertIn("bad.template", out)
        finally:
            td.cleanup()

    def test_template_validation_all_invalid_fails(self):
        invalid_template = (
            "bad.template",
            "".join(
                [
                    "# invalid template for testing\n",
                    "ID_ref= bad\n",
                    "Dnu_ref= 10\n",
                    "epsilon_ref= 1\n",
                    "numax_ref= 100\n",
                    "# Frequency Height Width\n",
                    "100.0 0.1\n",
                ]
            ),
        )

        td, root = self._make_sandbox_with_templates([invalid_template])
        try:
            example_cfg = (
                self.repo_root
                / "Configurations"
                / "examples_cfg"
                / "main.cfg.freeDP_curvepmodes.v3_GRANscaled"
            )
            patch_main_cfg = getattr(self.smoke_mod, "patch_main_cfg")
            patched_lines = patch_main_cfg(
                example_cfg,
                repo_root=self.repo_root,
                fixed_template="all",
                tobs_days=2.0,
                cadence_sec=120.0,
                nspectra=1,
                nrealisation=1,
                disable_plots=True,
                disable_modelfiles=True,
            )
            patched_cfg = root / "main.cfg"
            patched_cfg.write_text("".join(patched_lines), encoding="utf-8")

            cmd = [
                str(self.specsim),
                "--main_file",
                str(patched_cfg),
                "--noise_file",
                str(root / "Configurations" / "noise_Kallinger2014.cfg"),
                "--out_dir",
                str(root / "out"),
                "--force-create-output-dir",
                "1",
            ]
            rc, out = run(cmd, cwd=root, timeout=600)
            self.assertNotEqual(rc, 0)
            self.assertIn("No valid template files found", out)
        finally:
            td.cleanup()

    def test_list_models(self):
        rc, out = run([str(self.specsim), "--list-models"], cwd=self.repo_root, timeout=60)
        self.assertEqual(rc, 0, msg=out)
        for name in (
            "generate_cfg_from_synthese_file_Wscaled_Alm",
            "generate_cfg_from_synthese_file_Wscaled_aj",
            "generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014",
            "asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014",
        ):
            self.assertIn(name, out)

    def test_describe_model(self):
        name = "generate_cfg_from_synthese_file_Wscaled_aj"
        rc, out = run(
            [str(self.specsim), "--describe-model", name],
            cwd=self.repo_root,
            timeout=60,
        )
        self.assertEqual(rc, 0, msg=out)
        self.assertIn(f"Model: {name}", out)
        self.assertIn("Parameter order", out)
        self.assertIn("delta0l_percent", out)

    def test_registry_rejects_unknown_model(self):
        cfg_text = "".join(
            [
                "random 1\n",
                "unknown_model\n",
                "NONE\n",
                "Dnu epsilon delta0l_percent HNR a1ovGamma Gamma_at_numax a2 a3 a4 a5 a6 beta_asym i H_spread nu_spread Gamma_spread do_flat_noise\n",
                "70 0.5 1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 0 0\n",
                "70 0.5 1 10 0.6 1 0.1 -0.1 0.15 0.2 0.05 10 60 0 0 0 0\n",
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n",
                "Tobs Cadence Naverage Nrealisation\n",
                "2 120 1 1\n",
                "1\n",
                "0\n",
                "0\n",
                "0\n",
                "0\n",
            ]
        )
        td, root, out_dir, cfg_path, rc, out = self._run_cfg_text_in_sandbox(cfg_text)
        try:
            self.assertNotEqual(rc, 0)
            self.assertIn("registry-only mode", out)
        finally:
            td.cleanup()


if __name__ == "__main__":
    unittest.main()
