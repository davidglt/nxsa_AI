# nxsa_AI — X-ray Source Detection in XMM-Newton with HDBSCAN

> **Master's Thesis · Máster Universitario en Inteligencia Artificial · UNIR 2026**  
> *Detección de fuentes de rayos X en XMM-Newton con HDBSCAN*  
> **Author:** David González López-Tercero · **Supervisor:** Miryam Liliana Chaves Acero

---

## Overview

This repository contains the full software developed for the Master's Thesis on unsupervised X-ray source detection in XMM-Newton EPIC-pn event files using the **HDBSCAN** (Hierarchical Density-Based Spatial Clustering of Applications with Noise) algorithm.

The goal is to explore whether a density-based machine learning approach can complement the official XMM-Newton **Science Analysis Software (SAS)** pipeline by identifying X-ray source candidates directly from raw event files, without requiring prior knowledge of the number or geometry of the sources.

### Key results (validation campaign — 100 observations, 400 executions)

| Metric | Mean value |
|---|---|
| Precision | 0.22 |
| Recall | 0.25 |
| F1 score | 0.20 |
| Median matching distance | ~3.35 arcsec |

The `Brecallplus` hyperparameter profile achieved the best average F1. A subset of observations reached F1 values close to 0.6, with geometrically consistent matches to SAS detections.

---

## Repository structure

```
nxsa_AI/
├── src/
│   ├── create_dataset_v1.8.py         # Builds the observation dataset index
│   ├── run_hdbscan_v1.7.py            # Batch execution of HDBSCAN across observations
│   ├── hdbscan_ccd_v1.6.py            # Core module: event filtering, feature space, HDBSCAN per CCD, region generation
│   ├── collect_hdbscanval_v1.6.py     # Collects and aggregates per-observation validation results
│   ├── hdbscan_ccd_v2.6.ipynb         # Interactive notebook: single-observation inspection and calibration
│   └── hdbscanval_analysis_v1.5.ipynb # Global analysis notebook: metrics, plots, hyperparameter comparison
├── requirements/
│   └── requirements.txt               # Python dependencies
├── utils/
│   └── download_dataset_files.sh      # Helper script to download test observation files
├── TFM/                               # Master's Thesis document (PDF)
├── LICENSE.md                         # GPLv3
└── README.md
```

---

## Methodology

The workflow processes XMM-Newton EPIC-pn **full-frame** observations through the following stages:

1. **Dataset preparation** — selection and indexing of public observations from the XMM-Newton Science Archive (XSA)
2. **Event filtering** — GTI (Good Time Interval) background filtering + scientific event selection (FLAG, PATTERN, PI)
3. **Feature space construction** — spatiotemporal feature space (detector X/Y, time, PI, PATTERN, LABEL)
4. **Normalization and scaling** — standardisation of heterogeneous variables
5. **HDBSCAN per CCD** — independent clustering on each of the 12 EPIC-pn CCD panels
6. **Cluster filtering** — quality criteria applied to accepted clusters
7. **Region generation** — conversion of clusters to DS9-compatible astronomical regions
8. **Validation against SAS** — optimal linear assignment matching + precision / recall / F1 / angular distance
9. **Global analysis** — aggregated results over 100 observations × 4 hyperparameter profiles

---

## Installation

### Requirements

- Python 3.11
- Linux or macOS recommended (Windows via Git Bash)

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/davidglt/nxsa_AI.git
cd nxsa_AI

# 2. Download test observation files (Linux / Git Bash)
cd utils
./download_dataset_files.sh
cd ..

# 3. Create and activate a virtual environment
python3.11 -m venv .venv
source .venv/bin/activate          # Linux/macOS
# .venv\Scripts\activate.bat       # Windows

# 4. Install dependencies
pip install -r requirements/requirements.txt
```

---

## Usage

### 1 · Build the dataset index

```bash
python src/create_dataset_v1.8.py
```

Creates the CSV index of available observations to be processed.

### 2 · Run the batch validation campaign

```bash
python src/run_hdbscan_v1.7.py
```

Executes `hdbscan_ccd_v1.6.py` across all indexed observations and hyperparameter profiles. Results (metrics CSV + DS9 region files) are written to the output directory.

### 3 · Collect and aggregate results

```bash
python src/collect_hdbscanval_v1.6.py
```

Merges per-observation CSVs into a single consolidated results file.

### 4 · Interactive inspection (notebooks)

```bash
cd src
jupyter notebook
```

- **`hdbscan_ccd_v2.6.ipynb`** — single-observation exploration: event visualisation, GTI filtering, cluster overlay vs. SAS regions
- **`hdbscanval_analysis_v1.5.ipynb`** — global analysis: F1 distribution, precision-recall scatter, per-profile comparison, GTI impact, metric correlation heatmap

---

## Hyperparameter profiles

Four profiles were evaluated in the validation campaign:

| Profile | Description |
|---|---|
| `Bdefault` | Baseline HDBSCAN configuration |
| `Bprecisionplus` | Higher min_cluster_size — favours precision |
| `Brecallplus` | Lower min_cluster_size — favours recall · **best average F1** |
| `Bbalanced` | Intermediate trade-off configuration |

---

## Data

Observation data are publicly available via the [XMM-Newton Science Archive (XSA)](https://www.cosmos.esa.int/web/xmm-newton/xsa). The `utils/download_dataset_files.sh` script automates the download of the test dataset used in the validation campaign.

Input files per observation:
- **FITS event file** — EPIC-pn calibrated events (PPS product)
- **DS9 region file** — SAS source detection regions (reference for validation)

---

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE.md](LICENSE.md) for details.

---

## Acknowledgements

Special thanks to the XMM-Newton Science Operations Centre team at ESAC (Madrid):

- **José Vicente Perea**
- **Pedro Rodríguez**
- **María Santos**
- **Norbert Schartel**

Their scientific and technical expertise in XMM-Newton instrumentation, data analysis, and astrophysical interpretation was essential throughout this work.
