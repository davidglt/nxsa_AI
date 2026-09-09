#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

VERSION = "v1.7"
DEFAULT_DATASET_DIR = "../dataset"
DEFAULT_SCRIPT = "hdbscan_ccd_v1.6.py"
DEFAULT_LOG_DIR = "logs"
ALL_PROFILES = ["Brecallplus", "Bbalanced", "Brecallpp", "Bprecision"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the HDBSCAN CCD pipeline over all observations in a dataset. "
            "If no profile is given, all profiles are executed for each observation."
        )
    )
    parser.add_argument(
        "--dataset-dir",
        default=DEFAULT_DATASET_DIR,
        help=f"Dataset directory with observation subdirectories (default: {DEFAULT_DATASET_DIR})",
    )
    parser.add_argument(
        "--script",
        default=DEFAULT_SCRIPT,
        help=f"Pipeline script to execute (default: {DEFAULT_SCRIPT})",
    )
    parser.add_argument(
        "--log-dir",
        default=DEFAULT_LOG_DIR,
        help=f"Directory where log files will be written (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--profile",
        choices=ALL_PROFILES,
        help="Run only one profile. By default, all profiles are run.",
    )
    return parser.parse_args()


def run_one(obs: str, profile: str, dataset_dir: Path, pipeline_script: Path, log_dir: Path) -> int:
    log_file = log_dir / f"{obs}_{profile}.log"
    print(f"  - Perfil {profile} > {log_file}")
    cmd = [
        sys.executable,
        str(pipeline_script),
        "--dataset-dir",
        str(dataset_dir),
        "--profile",
        profile,
        obs,
    ]
    with log_file.open("w", encoding="utf-8") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] Fallo en {obs} perfil {profile}. Revisar {log_file}")
    else:
        print(f"  [OK] {obs} perfil {profile}")
    return result.returncode


def main() -> int:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir).expanduser().resolve()
    pipeline_script = Path(args.script).expanduser().resolve()
    log_dir = Path(args.log_dir).expanduser().resolve()

    print(f"run_hdbscan VERSION={VERSION}")
    print(f"dataset_dir={dataset_dir}")
    print(f"pipeline_script={pipeline_script}")
    print(f"log_dir={log_dir}")

    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory does not exist: {dataset_dir}")
    if not dataset_dir.is_dir():
        raise NotADirectoryError(f"Dataset path is not a directory: {dataset_dir}")
    if not pipeline_script.exists():
        raise FileNotFoundError(f"Pipeline script does not exist: {pipeline_script}")

    log_dir.mkdir(parents=True, exist_ok=True)

    profiles = [args.profile] if args.profile else list(ALL_PROFILES)
    if args.profile is None:
        print("[INFO] No se ha indicado perfil. Se ejecutarán todos: " + ", ".join(profiles))
    else:
        print(f"[INFO] Perfil solicitado: {args.profile}")

    obs_dirs = sorted([p for p in dataset_dir.iterdir() if p.is_dir()])
    failures = 0

    for obs_dir in obs_dirs:
        obs = obs_dir.name
        pps_dir = obs_dir / "pps"
        if not pps_dir.exists():
            print(f"[SKIP] {obs} sin pps")
            continue

        print()
        print(f"[OBS] Procesando {obs}")
        for profile in profiles:
            failures += int(run_one(obs, profile, dataset_dir, pipeline_script, log_dir) != 0)

    print()
    print("[INFO] Fin del procesamiento")
    if failures:
        print(f"[WARN] Ejecuciones con error: {failures}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
