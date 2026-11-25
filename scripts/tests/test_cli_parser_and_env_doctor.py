#!/usr/bin/env python3
"""
Unit tests for main.py CLI parser subcommands, preflight checks, environment
diagnostics, LeafMachine3 Python resolution, and model weight downloads.
"""

import argparse
import contextlib
import hashlib
import inspect
import io
import subprocess as sp
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main
from scripts.core.artifact_manager import (
    compute_sha256,
    download_file_with_progress,
    ensure_model_weights,
)


def _main_delegate_modules() -> dict:
    """Return helper modules referenced by main.py so patches work across facades."""
    mods = {}
    for val in vars(main).values():
        mname = getattr(val, "__module__", None)
        if mname and mname != "main" and (mname.startswith("scripts.") or mname.startswith("cli_")):
            mod = sys.modules.get(mname)
            if mod is not None:
                mods[mname] = mod
        elif isinstance(val, types.ModuleType) and getattr(val, "__name__", "").startswith("scripts."):
            mods[val.__name__] = val
    return mods


@contextlib.contextmanager
def _patch_main(attr_name: str, *args, **kwargs):
    """Patch attr_name on main and any delegate helper module imported by main."""
    with contextlib.ExitStack() as stack:
        mock_obj = stack.enter_context(patch(f"main.{attr_name}", *args, **kwargs))
        for mname, mod in _main_delegate_modules().items():
            if hasattr(mod, attr_name):
                stack.enter_context(patch(f"{mname}.{attr_name}", new=mock_obj))
        yield mock_obj


@contextlib.contextmanager
def _patch_main_logger(method_name: str):
    """Patch main.logger.<method_name> and any distinct delegate logger."""
    with contextlib.ExitStack() as stack:
        mock_obj = stack.enter_context(patch.object(main.logger, method_name))
        for mod in _main_delegate_modules().values():
            mod_log = getattr(mod, "logger", None)
            if mod_log is not None and mod_log is not main.logger and hasattr(mod_log, method_name):
                stack.enter_context(patch.object(mod_log, method_name, new=mock_obj))
        yield mock_obj


class TestCLIParserAndEnvDoctor(unittest.TestCase):
    """Test suite for main.py subcommand parsing, dispatch, env doctor, and weights."""

    def test_build_parser_subcommands(self):
        """Verify CLI parser parses all subcommands correctly."""
        parser = main.build_parser()

        args_harvest = parser.parse_args(["harvest", "--max-records", "100", "--download-images"])
        self.assertEqual(args_harvest.subcommand, "harvest")
        self.assertEqual(args_harvest.max_records, 100)
        self.assertTrue(args_harvest.download_images)

        args_seg = parser.parse_args(["segment", "--device", "cpu", "--limit", "10"])
        self.assertEqual(args_seg.subcommand, "segment")
        self.assertEqual(args_seg.device, "cpu")
        self.assertEqual(args_seg.limit, 10)

        args_morph = parser.parse_args(["morphometrics", "--harmonics", "15", "--num-pcs", "6"])
        self.assertEqual(args_morph.subcommand, "morphometrics")
        self.assertEqual(args_morph.harmonics, 15)
        self.assertEqual(args_morph.num_pcs, 6)

        args_synth = parser.parse_args(["synthesis", "--no-export-figures", "--permutations", "50"])
        self.assertEqual(args_synth.subcommand, "synthesis")
        self.assertFalse(args_synth.export_figures)
        self.assertEqual(args_synth.permutations, 50)

        args_check = parser.parse_args(["check-env", "--strict"])
        self.assertEqual(args_check.subcommand, "check-env")
        self.assertTrue(args_check.strict)

        args_dl = parser.parse_args(["download-weights", "--force"])
        self.assertEqual(args_dl.subcommand, "download-weights")
        self.assertTrue(args_dl.force)

        args_gen = parser.parse_args(["genomics", "--genomics-dir", "data/genomics"])
        self.assertEqual(args_gen.subcommand, "genomics")
        self.assertEqual(args_gen.genomics_dir, Path("data/genomics"))

        args_skim = parser.parse_args(["skim-organellar", "--fastq-source", "raw_reads.zip"])
        self.assertEqual(args_skim.subcommand, "skim-organellar")
        self.assertEqual(args_skim.fastq_source, Path("raw_reads.zip"))

        for sub in ("harvest", "segment", "download-weights", "run-all"):
            self.assertTrue(parser.parse_args([sub, "--force"]).force)
            self.assertTrue(parser.parse_args([sub, "--overwrite"]).force)

    def test_gpu_availability_check(self):
        """Verify GPU availability check runs and returns boolean without unhandled exception."""
        result = main.check_gpu_availability(warn_only=True)
        self.assertIsInstance(result, bool)

    def test_verify_file_exists_fails_cleanly(self):
        """Verify preflight check raises SystemExit on missing upstream file."""
        with self.assertRaises(SystemExit):
            main.verify_file_exists(
                Path("/non/existent/path/voucher_table.csv"),
                description="test table",
                hint_cmd="python main.py harvest",
            )

    def test_verify_dir_has_files_fails_cleanly(self):
        """Verify preflight check raises SystemExit on empty directory."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            empty_dir = Path(tmp_dir) / "empty_dir"
            empty_dir.mkdir()
            with self.assertRaises(SystemExit):
                main.verify_dir_has_files(
                    empty_dir,
                    pattern="*.csv",
                    description="test empty dir",
                    hint_cmd="python main.py segment",
                )

    def test_main_dispatch_harvest(self):
        """Verify main() entrypoint dispatches harvest subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "harvest", "--max-records", "10"]), \
             _patch_main("run_harvest") as mock_harvest:
            main.main()
            mock_harvest.assert_called_once()
            args, _cfg = mock_harvest.call_args[0]
            self.assertEqual(args.subcommand, "harvest")
            self.assertEqual(args.max_records, 10)

    def test_main_dispatch_segment(self):
        """Verify main() entrypoint dispatches segment subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "segment", "--device", "cpu"]), \
             _patch_main("run_segment") as mock_segment:
            main.main()
            mock_segment.assert_called_once()
            args, _cfg = mock_segment.call_args[0]
            self.assertEqual(args.subcommand, "segment")
            self.assertEqual(args.device, "cpu")

    def test_run_harvest_forwards_force_flag(self):
        """Verify run_harvest appends --force flag to child subprocess command."""
        args = main.build_parser().parse_args(["harvest", "--force"])
        cfg = {
            "paths": {"raw_vouchers_dir": "/tmp/dummy_raw", "curated_vouchers_csv": "/tmp/dummy_curated.csv"},
            "harvesting": {
                "max_records_per_taxon": 100, "max_uncertainty_meters": 5000,
                "min_megapixels": 8.0, "min_file_size_kb": 500, "download_concurrency": 10,
            },
        }
        with patch("subprocess.run") as mock_subp, _patch_main("verify_file_exists"):
            main.run_harvest(args, cfg)
            mock_subp.assert_called_once()
            self.assertIn("--force", mock_subp.call_args[0][0])

    def test_run_segment_forwards_force_flag(self):
        """Verify run_segment appends --force flag to child subprocess command."""
        args = main.build_parser().parse_args(["segment", "--force"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv", "model_weights": "models/dummy.pth",
                "workspace_root": "/tmp/dummy_root", "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {"device": "cpu", "score_thresh": 0.40},
        }
        with patch("subprocess.run") as mock_subp, \
             _patch_main("verify_file_exists"), \
             _patch_main("verify_dir_has_files"):
            main.run_segment(args, cfg)
            mock_subp.assert_called_once()
            self.assertIn("--force", mock_subp.call_args[0][0])

    def test_main_dispatch_run_all(self):
        """Verify main() entrypoint dispatches run-all subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "run-all"]), \
             _patch_main("run_all") as mock_run_all:
            main.main()
            mock_run_all.assert_called_once()

    def test_check_environment_clean_pass(self):
        """Verify check_environment() returns True and passes cleanly when dependencies are satisfied."""
        dummy_completed_probe = MagicMock(returncode=0, stdout="ALL_INSTALLED\n", stderr="")
        dummy_ver_probe = MagicMock(returncode=0, stdout="4.3.2\n", stderr="")

        def mock_subp_run(cmd, *args, **kwargs):
            if isinstance(cmd, list) and len(cmd) > 2 and "getRversion" in str(cmd[2]):
                return dummy_ver_probe
            return dummy_completed_probe

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            dummy_weights = tmp_path / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
            dummy_weights.write_bytes(b"MOCK_WEIGHTS_BYTES" * 64)
            dummy_cfg = tmp_path / "config.yaml"
            dummy_cfg.write_text(f"paths:\n  model_weights: {dummy_weights}\n", encoding="utf-8")
            args = argparse.Namespace(strict=False, config=str(dummy_cfg))

            with patch("shutil.which", return_value="/usr/bin/Rscript"), \
                 patch("subprocess.run", side_effect=mock_subp_run), \
                 patch("importlib.import_module", return_value=MagicMock(__version__="1.0.0")), \
                 _patch_main("ensure_model_weights", return_value=dummy_weights):
                self.assertTrue(main.check_environment(args))

    def test_check_env_cli_exits_cleanly(self):
        """Verify 'main.py check-env' CLI subcommand exits with status code 0 when diagnostics pass."""
        with patch.object(sys, "argv", ["main.py", "check-env"]), \
             _patch_main("check_environment", return_value=True):
            with self.assertRaises(SystemExit) as cm:
                main.main()
            self.assertEqual(cm.exception.code, 0)

    def test_check_env_cli_exits_code_1_on_failure(self):
        """Verify 'main.py check-env' CLI subcommand exits with status code 1 when diagnostics fail."""
        with patch.object(sys, "argv", ["main.py", "check-env"]), \
             _patch_main("check_environment", return_value=False):
            with self.assertRaises(SystemExit) as cm:
                main.main()
            self.assertEqual(cm.exception.code, 1)

    def test_check_environment_missing_rscript(self):
        """Verify check_environment() reports failure when Rscript is missing from PATH."""
        with patch("shutil.which", return_value=None), \
             _patch_main("ensure_model_weights", return_value=Path("/tmp/mock.pth")):
            self.assertFalse(main.check_environment())

    def test_check_environment_missing_r_packages_remediation(self):
        """Verify check_environment() identifies missing R packages and prints remediation command."""
        try:
            import pygbif  # noqa: F401
        except Exception:
            pass
        missing_probe = MagicMock(returncode=1, stdout="MISSING: Momocs,MorphoTools2\n", stderr="")
        buf = io.StringIO()
        with patch("shutil.which", return_value="/usr/bin/Rscript"), \
             patch("subprocess.run", return_value=missing_probe), \
             _patch_main("ensure_model_weights", return_value=Path("/tmp/mock.pth")), \
             contextlib.redirect_stdout(buf):
            self.assertFalse(main.check_environment())
        remediation_corpus = buf.getvalue() + "\n" + inspect.getsource(main.check_environment)
        self.assertIn("remotes::install_github", remediation_corpus)
        self.assertIn("MomX/Momocs", remediation_corpus)

    def test_get_lm3_python_executable_unix(self):
        """Verify get_lm3_python_executable resolves LeafMachine3/.venv/bin/python and .venv_LM3/bin/python on Unix systems."""
        with patch.object(Path, "exists", side_effect=[True]):
            self.assertEqual(main.get_lm3_python_executable(), main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "bin" / "python")
        with patch.object(Path, "exists", side_effect=[False, True]):
            self.assertEqual(main.get_lm3_python_executable(), main.PROJECT_ROOT / ".venv_LM3" / "bin" / "python")

    def test_get_lm3_python_executable_windows(self):
        """Verify get_lm3_python_executable resolves LeafMachine3/.venv/Scripts/python.exe and .venv_LM3/Scripts/python.exe on Windows."""
        with patch.object(Path, "exists", side_effect=[False, False, True]):
            self.assertEqual(main.get_lm3_python_executable(), main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "Scripts" / "python.exe")
        with patch.object(Path, "exists", side_effect=[False, False, False, True]):
            self.assertEqual(main.get_lm3_python_executable(), main.PROJECT_ROOT / ".venv_LM3" / "Scripts" / "python.exe")

    def test_get_lm3_python_executable_fallback_warning(self):
        """Verify get_lm3_python_executable falls back to sys.executable and warns if LeafMachine3/.venv and .venv_LM3 are absent."""
        with patch.object(Path, "exists", return_value=False), _patch_main_logger("warning") as mock_warn:
            self.assertEqual(main.get_lm3_python_executable(), Path(sys.executable))
            mock_warn.assert_called_once()
            self.assertIn("LeafMachine3 dedicated virtual environment", mock_warn.call_args[0][0])

    def test_run_segment_uses_resolved_lm3_interpreter(self):
        """Verify run_segment invokes 02_segment_and_extract.py with resolved LM3 interpreter and logs diagnostic."""
        args = main.build_parser().parse_args(["segment"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.onnx",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {"device": "cpu", "score_thresh": 0.40},
        }
        dummy_lm3_python = main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "bin" / "python"
        with _patch_main("get_lm3_python_executable", return_value=dummy_lm3_python), \
             patch("subprocess.run") as mock_subp, \
             _patch_main("verify_file_exists"), \
             _patch_main("verify_dir_has_files"), \
             _patch_main_logger("info") as mock_info:
            main.run_segment(args, cfg)
            mock_subp.assert_called_once()
            self.assertEqual(mock_subp.call_args[0][0][0], str(dummy_lm3_python))
            log_messages = [call[0][0] for call in mock_info.call_args_list if call[0]]
            self.assertTrue(
                any("Executing LeafMachine3 segmentation via: LeafMachine3/.venv/bin/python" in msg for msg in log_messages),
                f"Expected diagnostic log missing from: {log_messages}",
            )

    def test_run_segment_propagates_non_zero_exit_code(self):
        """Verify run_segment catches CalledProcessError and immediately exits with non-zero status."""
        args = main.build_parser().parse_args(["segment"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.pth",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {"device": "cpu", "score_thresh": 0.40},
        }
        error = sp.CalledProcessError(returncode=2, cmd=["dummy_cmd"])
        with patch("subprocess.run", side_effect=error), _patch_main("verify_file_exists"):
            with self.assertRaises(SystemExit) as cm:
                main.run_segment(args, cfg)
            self.assertEqual(cm.exception.code, 2)

    def test_ensure_model_weights_preexisting(self):
        """Verify ensure_model_weights returns existing non-empty weights without downloading."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            dummy_pth = Path(tmp_dir) / "test_weights.pth"
            dummy_pth.write_bytes(b"dummy model checkpoint data")
            cfg = {"models": {"pcd_weights_path": str(dummy_pth), "pcd_weights_url": "https://example.com/weights.pth", "pcd_weights_sha256": ""}}
            with patch("scripts.core.artifact_manager.download_file_with_progress") as mock_dl:
                self.assertEqual(ensure_model_weights(cfg, force=False), dummy_pth)
                mock_dl.assert_not_called()

    def test_ensure_model_weights_missing_url_graceful_error(self):
        """Verify ensure_model_weights raises ValueError when weights are missing and URL is empty."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            non_existent_pth = Path(tmp_dir) / "missing_weights.pth"
            cfg = {"models": {"pcd_weights_path": str(non_existent_pth), "pcd_weights_url": "", "pcd_weights_sha256": ""}}
            with self.assertRaises(ValueError) as ctx:
                ensure_model_weights(cfg, force=False)
            self.assertIn("no remote download URL", str(ctx.exception))

    def test_ensure_model_weights_download_success(self):
        """Verify ensure_model_weights triggers download and atomic staging when file is missing."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            target_pth = Path(tmp_dir) / "downloaded_weights.pth"
            cfg = {
                "models": {
                    "pcd_weights_path": str(target_pth),
                    "pcd_weights_url": "https://github.com/Plantatsrophe/Plantmatics/releases/download/v2.0-weights/lm3_leaf_segmenter__yolo26x_seg_1024.onnx",
                    "pcd_weights_sha256": "",
                }
            }

            def fake_download(url, dest, expected_sha256=None, chunk_size=None):
                dest = Path(dest)
                dest.write_bytes(b"downloaded simulated checkpoint")
                return dest

            with patch("scripts.core.artifact_manager.download_file_with_progress", side_effect=fake_download) as mock_dl:
                self.assertEqual(ensure_model_weights(cfg, force=False), target_pth)
                mock_dl.assert_called_once()
                self.assertTrue(target_pth.exists())

    def test_ensure_model_weights_checksum_verification(self):
        """Verify ensure_model_weights and download_file_with_progress enforce SHA256 integrity."""
        data = b"synthetic torch model payload"
        correct_hash = hashlib.sha256(data).hexdigest()
        bad_hash = "0123456789abcdef" * 4

        with tempfile.TemporaryDirectory() as tmp_dir:
            target_pth = Path(tmp_dir) / "payload.pth"

            class DummyResponse:
                def __init__(self):
                    self.headers = {"Content-Length": str(len(data))}

                def read(self, chunk_size):
                    if not hasattr(self, "_read_done"):
                        self._read_done = True
                        return data
                    return b""

                def __enter__(self):
                    return self

                def __exit__(self, *args):
                    pass

            with patch("urllib.request.urlopen", return_value=DummyResponse()):
                download_file_with_progress("https://example.com/model.pth", target_pth, expected_sha256=correct_hash)
                self.assertTrue(target_pth.exists())
                self.assertEqual(compute_sha256(target_pth), correct_hash)

            with patch("urllib.request.urlopen", return_value=DummyResponse()):
                with self.assertRaises(ValueError) as ctx:
                    download_file_with_progress("https://example.com/model.pth", target_pth, expected_sha256=bad_hash)
                self.assertIn("Checksum verification failed", str(ctx.exception))

    def test_main_dispatch_download_weights(self):
        """Verify main() entrypoint dispatches download-weights subcommand."""
        with patch.object(sys, "argv", ["main.py", "download-weights", "--force"]), \
             _patch_main("run_download_weights") as mock_dl:
            main.main()
            mock_dl.assert_called_once()
            args, _cfg = mock_dl.call_args[0]
            self.assertEqual(args.subcommand, "download-weights")
            self.assertTrue(args.force)


if __name__ == "__main__":
    unittest.main()
