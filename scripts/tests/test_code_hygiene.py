#!/usr/bin/env python3
"""
Static AST and text/grep hygiene test suite verifying archive isolation,
zero retired segmentation/helper tokens, and strict per-module line ceilings.
"""

from __future__ import annotations

import ast
import re
import sys
import unittest
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests import _hygiene_specs
from scripts.tests._hygiene_specs import (
    ARCHIVE_MARKER,
    EXCLUDED_DIRS,
    EXTRA_500_LINE_BOUNDED_MODULES,
    FORBIDDEN_LEGACY_LM2_TOKENS,
    FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS,
    RAW_STRICT_LINE_CEILINGS,
)

SELF_PATH = Path(__file__).resolve()
HYGIENE_SPECS_PATH = Path(_hygiene_specs.__file__).resolve()


def iter_active_source_files(include_shell: bool = False) -> List[Path]:
    """Return active .py, .R (and optionally .sh) files under scripts/ and root, excluding archive/venv dirs."""
    exts = ("*.py", "*.R", "*.sh") if include_shell else ("*.py", "*.R")
    files: List[Path] = []
    for ext in exts:
        for root_file in sorted(PROJECT_ROOT.glob(ext)):
            resolved = root_file.resolve()
            if resolved.is_file():
                files.append(resolved)
    scripts_dir = PROJECT_ROOT / "scripts"
    for ext in exts:
        for path in sorted(scripts_dir.rglob(ext)):
            resolved = path.resolve()
            if not resolved.is_file():
                continue
            if set(resolved.relative_to(PROJECT_ROOT).parts) & EXCLUDED_DIRS:
                continue
            files.append(resolved)
    return files


class TestCodeHygiene(unittest.TestCase):
    """Static AST and grep hygiene checks for archive isolation, ceilings, and retired token bans."""

    def test_no_active_files_reference_legacy_folders_via_grep(self):
        """Scan all active .py, .R, and .sh files and fail if any references the archive folder."""
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            for lineno, line in enumerate(content.splitlines(), start=1):
                if ARCHIVE_MARKER in line:
                    violations.append(f"{src_path.relative_to(PROJECT_ROOT)}:{lineno}: {line.strip()}")
        self.assertEqual(violations, [], f"Active source files must not reference '{ARCHIVE_MARKER}':\n" + "\n".join(violations))

    def test_no_files_or_folders_named_lm2_or_leafmachine2(self):
        """Verify that zero files or folders in the workspace use 'lm2' or 'leafmachine2' in their naming."""
        skip_roots = {".git", ".venv", ".venv_LM3", "LeafMachine3", "data", "runs", "outputs"}
        bad_paths: List[str] = []
        for entry in sorted(PROJECT_ROOT.iterdir()):
            name_low = entry.name.lower()
            if "lm2" in name_low or "leafmachine2" in name_low:
                bad_paths.append(str(entry.relative_to(PROJECT_ROOT)))
            if entry.name in skip_roots or not entry.is_dir():
                continue
            for sub in entry.rglob("*"):
                sub_low = sub.name.lower()
                if "lm2" in sub_low or "leafmachine2" in sub_low:
                    bad_paths.append(str(sub.relative_to(PROJECT_ROOT)))
        self.assertEqual(bad_paths, [], "No files or folders may use 'lm2' or 'leafmachine2':\n" + "\n".join(bad_paths))

    def test_no_active_files_reference_obsolete_sam2_lm2_or_pruned_helpers(self):
        """Scan active .py, .R, and .sh files (excluding hygiene spec/test files) and assert zero forbidden tokens."""
        violations: List[str] = []
        skip_paths = {SELF_PATH, HYGIENE_SPECS_PATH}
        for src_path in iter_active_source_files(include_shell=True):
            if src_path in skip_paths:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                for token in FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS:
                    if token in ("x_norm", "y_norm", "ingest_genetic_data"):
                        hit = bool(re.search(rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])", line))
                    elif token == "fig06c_":
                        hit = bool(re.search(r"fig06c_(?!cytonuclear_paralog_cytotype)", line))
                    else:
                        hit = token in line
                    if hit:
                        violations.append(f"{rel}:{lineno} [{token}]: {line.strip()}")
        self.assertEqual(violations, [], "Active source files contain forbidden retired tokens:\n" + "\n".join(violations))

    def test_new_public_genomics_modules_respect_500_line_ceiling(self):
        """Verify Step 01 to 07 bounded core, data-prep, analysis, figure, test, and hygiene modules stay within ceilings."""
        strict_ceilings = {PROJECT_ROOT / rel: lim for rel, lim in RAW_STRICT_LINE_CEILINGS}
        bounded_paths = (*strict_ceilings.keys(), *(PROJECT_ROOT / rel for rel in EXTRA_500_LINE_BOUNDED_MODULES))
        for mod_path in bounded_paths:
            self.assertTrue(mod_path.exists(), f"Missing expected bounded module: {mod_path}")
            line_count = len(mod_path.read_text(encoding="utf-8").splitlines())
            self.assertLess(line_count, 500, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds 500-line ceiling ({line_count} lines)")
            if mod_path in strict_ceilings:
                ceiling = strict_ceilings[mod_path]
                self.assertLess(line_count, ceiling, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds {ceiling}-line ceiling ({line_count} lines)")

    def test_step01_and_step06_split_posix_and_mmrr_hygiene(self):
        """Assert Step 01 and Step 06 test monolith removals, pure POSIX regexes in Step 06 R files, and zero abs(lat) IBE forgery."""
        self.assertFalse(
            (PROJECT_ROOT / "scripts" / "tests" / "test_voucher_harvester.py").exists(),
            "Legacy 3,930-line scripts/tests/test_voucher_harvester.py must remain deleted",
        )
        self.assertFalse(
            (PROJECT_ROOT / "scripts" / "tests" / "test_genomic_integration.py").exists(),
            "Legacy 2,426-line scripts/tests/test_genomic_integration.py must remain deleted",
        )
        perl_re = re.compile(r"\\\\[sdwb]|\(\?:")
        step06_r_files = sorted((PROJECT_ROOT / "scripts" / "analysis").glob("*06*.R")) + sorted(
            (PROJECT_ROOT / "scripts" / "figures").glob("*06*.R")
        )
        self.assertTrue(step06_r_files, "Expected Step 06 R files to audit")
        perl_hits: List[str] = []
        for r_path in step06_r_files:
            rel = r_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(r_path.read_text(encoding="utf-8").splitlines(), start=1):
                if perl_re.search(line):
                    perl_hits.append(f"{rel}:{lineno}: {line.strip()}")
        self.assertEqual(perl_hits, [], "Step 06 R files must use pure POSIX ERE instead of Perl escapes:\n" + "\n".join(perl_hits))
        mmrr_text = (PROJECT_ROOT / "scripts" / "analysis" / "step06_spatial_landscape_genomics.R").read_text(encoding="utf-8")
        self.assertNotIn("abs(lat", mmrr_text, "step06_spatial_landscape_genomics.R must never substitute abs(lat into eco_mat")

    def test_data_cytology_directory_stays_flat(self):
        """Verify `data/cytology/` exists and contains zero nested subdirectories."""
        cytology_dir = PROJECT_ROOT / "data" / "cytology"
        self.assertTrue(cytology_dir.is_dir(), f"Missing expected cytology directory: {cytology_dir}")
        subdirs = [p.name for p in cytology_dir.iterdir() if p.is_dir()]
        self.assertEqual(subdirs, [], f"data/cytology/ must remain flat with zero subdirectories: {subdirs}")

    def test_data_taxonomy_directory_stays_flat(self):
        """Verify `data/taxonomy/` contains zero nested subdirectories when present."""
        taxonomy_dir = PROJECT_ROOT / "data" / "taxonomy"
        if taxonomy_dir.exists():
            self.assertTrue(taxonomy_dir.is_dir(), f"Expected directory at {taxonomy_dir}")
            subdirs = [p.name for p in taxonomy_dir.iterdir() if p.is_dir()]
            self.assertEqual(subdirs, [], f"data/taxonomy/ must remain flat with zero subdirectories: {subdirs}")

    def test_no_hardcoded_workstation_legacy_repo_path(self):
        """Scan all active .py, .R, and .sh files and assert zero hardcoded workstation legacy repo paths."""
        legacy_abs_path = "/" + "home/brandon/" + "Packera_dubia_morphometrics"
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                if legacy_abs_path in line:
                    violations.append(f"{rel}:{lineno}: {line.strip()}")
        self.assertEqual(violations, [], "Active source files hardcode the legacy workstation repo path:\n" + "\n".join(violations))

    def test_no_duplicate_python_defs_in_same_scope(self):
        """Assert zero shadowed or duplicate def/async def/class definitions within any module or class body."""
        def _is_overload(node: ast.AST) -> bool:
            return any(
                (isinstance(d, ast.Name) and d.id == "overload") or (isinstance(d, ast.Attribute) and d.attr == "overload")
                for d in getattr(node, "decorator_list", [])
            )

        def _check_body(body: List[ast.stmt], scope_label: str, dupes: List[str]) -> None:
            seen: dict[str, int] = {}
            for stmt in body:
                if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    if _is_overload(stmt):
                        continue
                    if stmt.name in seen:
                        dupes.append(f"{scope_label}:{stmt.lineno} duplicates {stmt.name} (first at L{seen[stmt.name]})")
                    else:
                        seen[stmt.name] = stmt.lineno
                    if isinstance(stmt, ast.ClassDef):
                        _check_body(stmt.body, f"{scope_label}.{stmt.name}", dupes)

        dupes: List[str] = []
        py_files = [PROJECT_ROOT / "main.py"] + sorted(
            p for p in (PROJECT_ROOT / "scripts").rglob("*.py") if ARCHIVE_MARKER not in p.parts
        )
        for py_path in py_files:
            if py_path.exists():
                tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
                _check_body(tree.body, str(py_path.relative_to(PROJECT_ROOT)), dupes)
        self.assertEqual(dupes, [], "Duplicate Python definitions found in same scope:\n" + "\n".join(dupes))

    def test_no_broken_symlinks(self):
        """Assert zero broken symlinks across models/, scripts/, and runs/latest*."""
        broken: List[str] = []
        for root_dir in (PROJECT_ROOT / "models", PROJECT_ROOT / "scripts"):
            if root_dir.exists():
                for entry in root_dir.rglob("*"):
                    if "__pycache__" not in entry.parts and entry.is_symlink() and not entry.exists():
                        broken.append(str(entry.relative_to(PROJECT_ROOT)))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for entry in runs_dir.glob("latest*"):
                if entry.is_symlink() and not entry.exists():
                    broken.append(str(entry.relative_to(PROJECT_ROOT)))
        self.assertEqual(broken, [], "Broken symlinks detected:\n" + "\n".join(broken))

    def test_no_pdf_figure_artifacts(self):
        """Assert zero Rplots.pdf or figure .pdf artifacts at repo root or figure directories."""
        found_pdfs: List[str] = []
        rplots = PROJECT_ROOT / "Rplots.pdf"
        if rplots.exists() or rplots.is_symlink():
            found_pdfs.append("Rplots.pdf")
        out_figs = PROJECT_ROOT / "outputs" / "figures"
        if out_figs.exists():
            found_pdfs.extend(str(p.relative_to(PROJECT_ROOT)) for p in out_figs.rglob("*.pdf"))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for fig_dir in runs_dir.glob("*/outputs/figures"):
                if not fig_dir.is_symlink():
                    found_pdfs.extend(str(p.relative_to(PROJECT_ROOT)) for p in fig_dir.rglob("*.pdf"))
        self.assertEqual(found_pdfs, [], "Unexpected PDF figure artifacts found:\n" + "\n".join(sorted(found_pdfs)))


if __name__ == "__main__":
    unittest.main()
