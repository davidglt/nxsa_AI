#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

VERSION = "v1.6"
INPUT_PATTERN = "P*EPX*HDBSCANVAL.csv"
OUTPUT_NAME = "HDBSCANVAL_COLLECTION.csv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Colector de productos HDBSCANVAL.csv dentro de un dataset XMM-Newton"
    )
    p.add_argument(
        "dataset_dir",
        help="Directorio raíz del dataset que contiene subdirectorios <OBS>/pps/",
    )
    return p.parse_args()


def find_validation_csvs(dataset_dir: Path) -> list[Path]:
    return sorted(dataset_dir.glob(f"*/pps/{INPUT_PATTERN}"))


def read_csv_rows(csv_path: Path) -> list[dict]:
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def normalize_rows(rows: list[dict], source_csv: Path) -> list[dict]:
    normalized = []
    for i, row in enumerate(rows, start=1):
        r = dict(row)
        r["source_csv"] = str(source_csv)
        r["source_obs_dir"] = source_csv.parent.parent.name
        r["collector_version"] = VERSION
        r["collector_row_index"] = str(i)
        normalized.append(r)
    return normalized


def write_master_csv(out_path: Path, rows: list[dict]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with out_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "source_csv",
                "source_obs_dir",
                "collector_version",
                "collector_row_index",
            ])
        return

    fieldnames = []
    seen = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir).expanduser().resolve()
    if not dataset_dir.exists():
        raise FileNotFoundError(f"No existe el directorio dataset: {dataset_dir}")
    if not dataset_dir.is_dir():
        raise NotADirectoryError(f"La ruta no es un directorio: {dataset_dir}")

    csv_paths = find_validation_csvs(dataset_dir)
    print(f"Colector HDBSCANVAL {VERSION}")
    print(f"dataset_dir={dataset_dir}")
    print(f"CSV de validación encontrados={len(csv_paths)}")

    all_rows = []
    for csv_path in csv_paths:
        rows = read_csv_rows(csv_path)
        if not rows:
            print(f"[WARN] CSV vacío: {csv_path}", file=sys.stderr)
            continue
        norm = normalize_rows(rows, csv_path)
        all_rows.extend(norm)
        print(f"  + {csv_path} -> {len(norm)} fila(s)")

    out_path = dataset_dir / OUTPUT_NAME
    write_master_csv(out_path, all_rows)
    print(f"CSV maestro escrito: {out_path}")
    print(f"Filas totales: {len(all_rows)}")


if __name__ == "__main__":
    main()
