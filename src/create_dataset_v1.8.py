#!/usr/bin/env python3
from __future__ import annotations

import argparse
import random
import shutil
from pathlib import Path

from astropy.io import fits

VERSION = "v1.8"
VALID_SUBMODE = "PrimeFullWindow"
EVENT_PATTERN = "P{obs}PNS*PIEVLI0000.FTZ"
DEFAULT_BASE_DIR = "../dataset_all"
DEFAULT_DEST_DIR = "../dataset"
DEFAULT_TEST_DIR = "../dataset_test"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a random XMM-Newton dataset filtered to EPIC-pn PNS "
            "observations with SUBMODE=PrimeFullWindow"
        )
    )
    parser.add_argument(
        "n",
        type=int,
        help="Total number of valid observations to copy (including test observations if enabled)",
    )
    parser.add_argument(
        "--base-dir",
        default=DEFAULT_BASE_DIR,
        help=f"Base directory with observation subdirectories (default: {DEFAULT_BASE_DIR})",
    )
    parser.add_argument(
        "--dest-dir",
        default=DEFAULT_DEST_DIR,
        help=f"Destination dataset directory (default: {DEFAULT_DEST_DIR})",
    )
    parser.add_argument(
        "--test-dir",
        default=DEFAULT_TEST_DIR,
        help=f"Directory with fixed test observations (default: {DEFAULT_TEST_DIR})",
    )
    parser.add_argument(
        "--no-include-test",
        dest="include_test",
        action="store_false",
        help="Do not copy fixed test observations (by default they are included)",
    )
    parser.set_defaults(include_test=True)
    return parser.parse_args()


def read_mode(event_file: Path) -> str:
    with fits.open(event_file, mode="readonly", ignore_missing_end=True) as hdul:
        for hdu in hdul:
            header = hdu.header
            for key in ("SUBMODE", "DATAMODE", "MODE"):
                value = header.get(key)
                if value is not None:
                    return str(value).strip()
    return ""


def find_valid_event(obs_dir: Path) -> tuple[Path | None, str, str]:
    obs = obs_dir.name
    pps_dir = obs_dir / "pps"
    if not pps_dir.exists():
        return None, "sin pps", ""

    candidates = sorted(pps_dir.glob(EVENT_PATTERN.format(obs=obs)))
    if not candidates:
        return None, "sin PNS science event", ""

    event_file = candidates[0]
    try:
        mode = read_mode(event_file)
    except Exception as exc:
        return None, f"error leyendo FITS: {event_file.name} ({exc})", ""

    if not mode:
        return None, f"SUBMODE/DATAMODE/MODE vacío: {event_file.name}", ""

    if mode != VALID_SUBMODE:
        return None, f"modo no {VALID_SUBMODE}: {mode} [{event_file.name}]", mode

    return event_file, "", mode


def clear_destination(dest_dir: Path) -> None:
    if not dest_dir.exists():
        dest_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Destino creado: {dest_dir}")
        return

    print(f"[INFO] Limpiando destino: {dest_dir}")
    for child in list(dest_dir.iterdir()):
        try:
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
        except Exception as exc:
            print(f"[WARN] No se pudo borrar {child}: {exc}")


def copy_observation(obs_dir: Path, dest_dir: Path) -> None:
    target_dir = dest_dir / obs_dir.name
    if target_dir.exists():
        shutil.rmtree(target_dir)
    shutil.copytree(obs_dir, target_dir)


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).expanduser().resolve()
    dest_dir = Path(args.dest_dir).expanduser().resolve()
    test_dir = Path(args.test_dir).expanduser().resolve()

    print(f"create_dataset VERSION={VERSION}")
    print(f"base_dir={base_dir}")
    print(f"dest_dir={dest_dir}")
    print(f"test_dir={test_dir}")
    print(f"include_test={args.include_test}")

    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory does not exist: {base_dir}")
    if not base_dir.is_dir():
        raise NotADirectoryError(f"Base path is not a directory: {base_dir}")

    clear_destination(dest_dir)

    included_obs: set[str] = set()
    selected = 0

    # Copiar primero las observaciones fijas de test, aplicando el mismo filtro científico.
    if args.include_test and test_dir.exists() and test_dir.is_dir():
        fixed_obs_dirs = sorted([p.resolve() for p in test_dir.iterdir() if p.is_dir()])
        for obs_dir in fixed_obs_dirs:
            if selected >= args.n:
                break
            obs = obs_dir.name
            if obs in included_obs:
                continue

            event_file, reason, mode = find_valid_event(obs_dir)
            if event_file is None:
                print(f"[SKIP] {obs} test: {reason}")
                continue

            copy_observation(obs_dir, dest_dir)
            included_obs.add(obs)
            selected += 1
            print(
                f"[OK]   {obs} copiada como test fija ({selected}/{args.n}), "
                f"modo={mode}, evt={event_file.name}"
            )

    # Completar con observaciones aleatorias del repositorio principal, evitando duplicados.
    candidate_obs_dirs = [
        p.resolve()
        for p in base_dir.iterdir()
        if p.is_dir() and p.name not in included_obs
    ]
    random.shuffle(candidate_obs_dirs)

    for obs_dir in candidate_obs_dirs:
        if selected >= args.n:
            break

        obs = obs_dir.name
        if obs in included_obs:
            continue

        event_file, reason, mode = find_valid_event(obs_dir)
        if event_file is None:
            print(f"[SKIP] {obs} {reason}")
            continue

        copy_observation(obs_dir, dest_dir)
        included_obs.add(obs)
        selected += 1
        print(f"[OK]   {obs} copiada ({selected}/{args.n}), modo={mode}, evt={event_file.name}")

    if selected < args.n:
        print(f"[WARN] Solo se copiaron {selected} observaciones válidas de {args.n} solicitadas")
    else:
        print(f"[INFO] Selección completada: {selected} observaciones válidas")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
