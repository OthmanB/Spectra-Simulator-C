import tempfile
import unittest
from pathlib import Path


class TestSmokePatchCfg(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[1]
        cls.smoke_mod = __import__(
            "scripts.smoke_test", fromlist=["patch_main_cfg", "strip_inline_comment"]
        )

    def test_strip_inline_comment(self):
        strip_inline_comment = getattr(self.smoke_mod, "strip_inline_comment")
        self.assertEqual(strip_inline_comment("foo # bar"), "foo")
        self.assertEqual(strip_inline_comment("foo#bar"), "foo")
        self.assertEqual(strip_inline_comment("  foo   #bar"), "foo")
        self.assertEqual(strip_inline_comment("foo"), "foo")
        self.assertEqual(strip_inline_comment("   "), "")

    def test_patch_main_cfg_minimal(self):
        patch_main_cfg = getattr(self.smoke_mod, "patch_main_cfg")

        cfg_text = "".join(
            [
                "# header\n",
                "grid 10 # forest\n",
                "generate_cfg_from_synthese_file_Wscaled_aj   /tmp/does_not_matter/8379927.in  # extra params\n",
                "all\n",
                "Dnu epsilon\n",
                "70 0.5\n",
                "70 0.5\n",
                "0 0\n",
                "Tobs Cadence Naverage Nrealisation\n",
                "100 120 1 1\n",
                "0\n",
                "1\n",
                "1\n",
                "0\n",
                "0 model\n",
            ]
        )

        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "main.cfg"
            src.write_text(cfg_text, encoding="utf-8")

            fixed_template = "12508433.template"
            patched = patch_main_cfg(
                src,
                repo_root=self.repo_root,
                fixed_template=fixed_template,
                tobs_days=2.0,
                cadence_sec=120.0,
                nspectra=1,
                nrealisation=1,
                disable_plots=True,
                disable_modelfiles=True,
            )
            patched_text = "".join(patched)

            # forest line forced to random 1
            self.assertIn("random 1\n", patched_text)

            # template selection pinned
            self.assertIn(f"{fixed_template}\n", patched_text)

            # step vector forced to all zeros and correct length
            self.assertIn("0 0\n", patched_text)

            # observation line patched
            self.assertIn("2.0 120.0 1 1\n", patched_text)

            # The extra_params should be rewritten to the repo-local infiles when an '*.in' token is present
            expected_infile = (
                self.repo_root / "Configurations" / "infiles" / "8379927.in"
            )
            self.assertTrue(expected_infile.exists())
            self.assertIn(str(expected_infile), patched_text)


if __name__ == "__main__":
    unittest.main()
