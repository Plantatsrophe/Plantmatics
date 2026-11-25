#!/usr/bin/env python3
"""Persistent per-stage warm Rscript worker session facade for the Plantmatics test suite."""

from __future__ import annotations

import atexit
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Dict, List, Mapping, Optional, Sequence, Union

try:
    from scripts.tests._r_session_preload import (
        _MANAGED_ENV_KEYS, _SCRIPT_TO_STAGE, _STAGE_CONFIG, _WORKER_R_CODE,
    )
except ImportError:
    from _r_session_preload import (  # type: ignore[no-redef]
        _MANAGED_ENV_KEYS, _SCRIPT_TO_STAGE, _STAGE_CONFIG, _WORKER_R_CODE,
    )

PROJECT_ROOT = Path(__file__).resolve().parents[2]

def rscript_available() -> bool:
    """Return True when Rscript is installed and fast Python-only audit mode is off."""
    for env_key in ("PLANTMATICS_FAST_AUDIT", "PACKERA_FAST_AUDIT"):
        if os.environ.get(env_key, "").strip() == "1":
            return False
    return shutil.which("Rscript") is not None


class WarmRSession:
    """Persistent warm Rscript process preloaded with a stage's R environment."""

    def __init__(self, stage_key: str) -> None:
        if stage_key not in _STAGE_CONFIG:
            raise ValueError(f"Unknown R stage key: {stage_key}")
        self.stage_key = stage_key
        self.preload_rel, self.entry_fn = _STAGE_CONFIG[stage_key]
        self._tmp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
        self._proc: Optional[subprocess.Popen[str]] = None
        self._job_seq = 0

    def _ensure_started(self) -> None:
        if self._proc is not None and self._proc.poll() is None:
            return
        self.close()
        rscript_bin = shutil.which("Rscript")
        if not rscript_bin:
            raise RuntimeError("Rscript binary not found on PATH.")
        self._tmp_dir = tempfile.TemporaryDirectory(prefix=f"packera_r_{self.stage_key}_")
        tmp_path = Path(self._tmp_dir.name)
        worker_file = tmp_path / "warm_worker.R"
        worker_file.write_text(_WORKER_R_CODE, encoding="utf-8")

        boot_env = os.environ.copy()
        boot_env["PACKERA_PRELOAD_R_SCRIPT"] = str(PROJECT_ROOT / self.preload_rel)
        for k in _MANAGED_ENV_KEYS:
            boot_env.pop(k, None)

        self._proc = subprocess.Popen(
            [rscript_bin, "--vanilla", str(worker_file)],
            cwd=str(PROJECT_ROOT),
            env=boot_env,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert self._proc.stdout is not None
        boot_lines: List[str] = []
        while True:
            line = self._proc.stdout.readline()
            if not line:
                rc = self._proc.wait()
                self.close()
                raise RuntimeError(
                    f"Warm R worker ({self.stage_key}) exited during boot (code {rc}):\n"
                    + "".join(boot_lines)
                )
            if "<<<PACKERA_BOOT_READY>>>" in line:
                break
            if line.startswith("<<<PACKERA_BOOT_ERROR>>>:"):
                rc = self._proc.wait()
                self.close()
                raise RuntimeError(
                    f"Warm R worker ({self.stage_key}) failed to preload {self.preload_rel} (code {rc}):\n"
                    + "".join(boot_lines)
                    + line
                )
            boot_lines.append(line)

    def _dispatch(
        self,
        *,
        mode: str,
        cli_args: Sequence[str] = (),
        script_path: Optional[Path] = None,
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, str]] = None,
        skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        self._ensure_started()
        assert self._proc is not None and self._proc.stdin is not None and self._proc.stdout is not None
        assert self._tmp_dir is not None

        self._job_seq += 1
        jid = f"{self._job_seq:04d}"
        jdir = Path(self._tmp_dir.name)

        args_file = jdir / f"job_{jid}.args"
        env_keys_file = jdir / f"job_{jid}.envk"
        env_vals_file = jdir / f"job_{jid}.envv"
        out_file = jdir / f"job_{jid}.out"
        err_file = jdir / f"job_{jid}.err"
        status_file = jdir / f"job_{jid}.status"
        spec_file = jdir / f"job_{jid}.spec"
        done_token = f"<<<PACKERA_RPC_DONE_{jid}>>>"

        args_file.write_text(
            "\n".join(str(a) for a in cli_args) + ("\n" if cli_args else ""),
            encoding="utf-8",
        )

        env_map: Dict[str, str] = {}
        for k in _MANAGED_ENV_KEYS:
            if env is not None and k in env:
                env_map[k] = str(env[k])
            elif env is None and k in os.environ:
                env_map[k] = os.environ[k]
            else:
                env_map[k] = "__UNSET__"
        if skip_figures:
            env_map["PLANTMATICS_SKIP_FIGURES"] = "1"
            env_map["PACKERA_SKIP_FIGURES"] = "1"
        if env is not None:
            for k, v in env.items():
                if k not in os.environ or os.environ.get(k) != v:
                    if "\n" not in k and "\n" not in str(v):
                        env_map[k] = str(v)

        keys_list = list(env_map.keys())
        vals_list = [env_map[k] for k in keys_list]
        env_keys_file.write_text("\n".join(keys_list) + "\n", encoding="utf-8")
        env_vals_file.write_text("\n".join(vals_list) + "\n", encoding="utf-8")
        out_file.write_text("", encoding="utf-8")
        err_file.write_text("", encoding="utf-8")

        req_cwd = str(Path(cwd).resolve()) if cwd is not None else str(PROJECT_ROOT)
        s_path = str(script_path.resolve()) if script_path is not None else ""

        spec_lines = [
            mode, self.entry_fn if mode == "cli" else "", req_cwd, s_path,
            str(args_file), str(env_keys_file), str(env_vals_file),
            str(out_file), str(err_file), str(status_file), done_token,
        ]
        spec_file.write_text("\n".join(spec_lines) + "\n", encoding="utf-8")

        self._proc.stdin.write(str(spec_file) + "\n")
        self._proc.stdin.flush()

        passthrough_lines: List[str] = []
        while True:
            line = self._proc.stdout.readline()
            if not line:
                rc = self._proc.wait() or 1
                out_txt = out_file.read_text(encoding="utf-8", errors="replace") if out_file.exists() else ""
                err_txt = err_file.read_text(encoding="utf-8", errors="replace") if err_file.exists() else ""
                err_txt = err_txt + "".join(passthrough_lines)
                self.close()
                return subprocess.CompletedProcess(["Rscript", self.preload_rel, *cli_args], rc, out_txt, err_txt)
            if done_token in line:
                break
            passthrough_lines.append(line)

        rc_txt = status_file.read_text(encoding="utf-8").strip() if status_file.exists() else "1"
        rc = 0 if rc_txt == "0" else 1
        out_txt = out_file.read_text(encoding="utf-8", errors="replace") if out_file.exists() else ""
        err_txt = err_file.read_text(encoding="utf-8", errors="replace") if err_file.exists() else ""
        if passthrough_lines:
            err_txt = "".join(passthrough_lines) + err_txt

        return subprocess.CompletedProcess(["Rscript", self.preload_rel if mode == "cli" else s_path, *cli_args], rc, out_txt, err_txt)

    def run_cli(
        self, cli_args: Sequence[Union[str, Path]], *, cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, str]] = None, skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        return self._dispatch(mode="cli", cli_args=[str(a) for a in cli_args], cwd=cwd, env=env, skip_figures=skip_figures)

    def run_script(
        self, script_or_code: Union[Path, str], *, cli_args: Sequence[Union[str, Path]] = (),
        cwd: Optional[Union[str, Path]] = None, env: Optional[Mapping[str, str]] = None, skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        self._ensure_started()
        assert self._tmp_dir is not None
        if isinstance(script_or_code, Path):
            s_path = script_or_code
        else:
            s_path = Path(self._tmp_dir.name) / f"inline_{self._job_seq + 1:04d}.R"
            s_path.write_text(script_or_code, encoding="utf-8")
        return self._dispatch(mode="script", cli_args=[str(a) for a in cli_args], script_path=s_path, cwd=cwd, env=env, skip_figures=skip_figures)

    def close(self) -> None:
        if self._proc is not None:
            try:
                if self._proc.poll() is None and self._proc.stdin is not None:
                    try:
                        self._proc.stdin.write("__QUIT__\n")
                        self._proc.stdin.flush()
                        self._proc.stdin.close()
                    except Exception:
                        pass
                    try:
                        self._proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        self._proc.kill()
                        self._proc.wait(timeout=2)
            finally:
                if self._proc.stdout is not None:
                    try:
                        self._proc.stdout.close()
                    except Exception:
                        pass
                self._proc = None
        if self._tmp_dir is not None:
            try:
                self._tmp_dir.cleanup()
            except Exception:
                pass
            self._tmp_dir = None


_SESSIONS: Dict[str, WarmRSession] = {}


def get_warm_r_session(stage_key: str) -> WarmRSession:
    session = _SESSIONS.get(stage_key)
    if session is None:
        session = WarmRSession(stage_key)
        _SESSIONS[stage_key] = session
    return session


def run_warm_r_cmd(
    cmd: Sequence[Union[str, Path]],
    *,
    stage_key: Optional[str] = None,
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Mapping[str, str]] = None,
    skip_figures: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Drop-in warm replacement for subprocess.run(['Rscript', ...])."""
    parts = [str(c) for c in cmd]
    if len(parts) >= 3 and parts[1] == "-e":
        sk = stage_key or "step03"
        return get_warm_r_session(sk).run_script(
            parts[2], cli_args=parts[3:], cwd=cwd, env=env, skip_figures=skip_figures
        )
    if len(parts) >= 2:
        script_name = Path(parts[1]).name
        resolved_stage = _SCRIPT_TO_STAGE.get(script_name)
        if resolved_stage is not None:
            return get_warm_r_session(resolved_stage).run_cli(
                parts[2:], cwd=cwd, env=env, skip_figures=skip_figures
            )
        sk = stage_key or "step04"
        return get_warm_r_session(sk).run_script(
            Path(parts[1]), cli_args=parts[2:], cwd=cwd, env=env, skip_figures=skip_figures
        )
    raise ValueError(f"Unsupported Rscript command sequence: {parts}")


def close_all_warm_r_sessions() -> None:
    for sess in list(_SESSIONS.values()):
        sess.close()
    _SESSIONS.clear()


atexit.register(close_all_warm_r_sessions)

__all__ = [
    "PROJECT_ROOT", "_MANAGED_ENV_KEYS", "_WORKER_R_CODE", "_STAGE_CONFIG", "_SCRIPT_TO_STAGE",
    "rscript_available", "WarmRSession", "_SESSIONS",
    "get_warm_r_session", "run_warm_r_cmd", "close_all_warm_r_sessions",
]
