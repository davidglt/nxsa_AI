#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import sys
import time as t
from pathlib import Path
from dataclasses import dataclass
from multiprocessing.pool import ThreadPool as Pool

import hdbscan
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from regions import Regions, CirclePixelRegion, PixCoord
from scipy.optimize import linear_sum_assignment


CCDRAWXMIN, CCDRAWXMAX = 1, 64
CCDRAWYMIN, CCDRAWYMAX = 2, 199
PIMIN, PIMAX = 300, 7000
NUMTHREADS = 8

APPLYBACKGROUNDGTI = True
BKGBINSECONDS = 100
BKGMETHOD = "robust"
BKGFIXEDRATE = 0.4
BKGSIGMA = 3.0
MINGTIBINS = 5

EXCLUDEPIRANGES: list[tuple[float, float]] = []

DETECTIONMODE = "spatiotemporal"  # 'spatial' o 'spatiotemporal'
TIMESCALESECONDS = 5000.0
TIMEWEIGHT = 0.2
HYPERPARAMETERPROFILE = "Brecallplus"
CLUSTERSELECTIONMETHOD = "eom"
METRIC = "euclidean"
MINCLUSTERSIZE = 30
MINSAMPLES = 15

MINPERSISTENCE = 0.03
MINMEANPROBABILITY = 0.08
MINCLUSTEREVENTS = 40
MINRADIUSARCSEC = 6.0
MAXRADIUSARCSEC = 60.0
MATCHTOLERANCEARCSEC = 15.0

@dataclass(frozen=True)
class HyperProfile:
    name: str
    product_code: str
    product_code_all: str
    timeweight: float
    minclustersize: int
    minsamples: int
    minpersistence: float
    minmeanprobability: float
    minclusterevents: int

PROFILES = {
    "Brecallplus": HyperProfile("Brecallplus", "000", "001", 0.2, 30, 15, 0.03, 0.08, 40),
    "Bbalanced":   HyperProfile("Bbalanced",   "010", "011", 0.1, 25, 12, 0.025, 0.07, 30),
    "Brecallpp":   HyperProfile("Brecallpp",   "020", "021", 0.1, 20, 10, 0.02,  0.06, 25),
    "Bprecision":  HyperProfile("Bprecision",  "030", "031", 0.2, 40, 20, 0.03,  0.08, 40),
}
DEFAULT_PROFILE = "Brecallplus"

VERSION = "v1.6"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Pipeline final HDBSCAN por CCD para observaciónes XMM-Newton EPIC-pn"
    )
    p.add_argument("obs", nargs="?", help="Número de observación XMM, p. ej., 0405320501")
    p.add_argument(
        "--dataset-dir",
        default="../dataset",
        help="Directorio raíz que contiene dataset/<OBS>/pps/ (por defecto: ../dataset)",
    )
    p.add_argument(
        "--profile",
        default=DEFAULT_PROFILE,
        choices=sorted(PROFILES.keys()),
        help="Perfil de hiperparámetros HDBSCAN",
    )
    args = p.parse_args()
    if not args.obs:
        p.error("debes indicar el número de observación, p. ej., python hdbscan_ccd.py 0405320501")
    return args


def autodetect_files(obs: str, dataset_dir: str, product_code: str = "000", product_code_all: str = "001"):
    pps = Path(dataset_dir).expanduser().resolve() / obs / "pps"
    if not pps.exists():
        raise FileNotFoundError(f"No existe el directorio PPS: {pps}")

    event_candidates = sorted(pps.glob(f"P{obs}PNS*PIEVLI0000.FTZ"))
    if not event_candidates:
        raise FileNotFoundError(f"No se encontró fichero de eventos PN para OBS={obs} en {pps}")
    if len(event_candidates) > 1:
        print("[INFO] Múltiples ficheros PN detectados; se usa el primero:", file=sys.stderr)
        for c in event_candidates:
            print(f"  - {c.name}", file=sys.stderr)
    event_file = event_candidates[0]

    stem = event_file.name.replace("PIEVLI0000.FTZ", "")
    exposure = stem.split("PNS")[-1]
    region_list = pps / f"P{obs}EPX000REGION0000.ASC"
    hdbscan_region_list = pps / f"P{obs}EPX{product_code}HDBSCANREG.ASC"
    hdbscan_region_list_all = pps / f"P{obs}EPX{product_code_all}HDBSCANREG.ASC"
    hdbscan_results_csv = pps / f"P{obs}EPX{product_code}HDBSCANVAL.csv"

    return {
        "pps": pps,
        "event_file": event_file,
        "exposure": exposure,
        "region_list": region_list,
        "hdbscan_region_list": hdbscan_region_list,
        "hdbscan_region_list_all": hdbscan_region_list_all,
        "hdbscan_results_csv": hdbscan_results_csv,
    }


def filter_background_gti(
    ev_list,
    bin_seconds=100,
    method="robust",
    fixed_rate=0.4,
    sigma_factor=3.0,
    min_gti_bins=5,
    bkg_pi_min=10000,
    bkg_pi_max=12000,
    bkg_pattern=0,
    sci_pi_min=PIMIN,
    sci_pi_max=PIMAX,
    sci_pattern_max=4,
    max_iter=10,
):
    time_all = np.asarray(ev_list["TIME"], dtype=float)
    t0 = float(time_all.min())
    t1 = float(time_all.max())

    bkg_sel = (
        (ev_list["FLAG"] == 0)
        & (ev_list["PATTERN"] == bkg_pattern)
        & (ev_list["PI"] >= bkg_pi_min)
        & (ev_list["PI"] <= bkg_pi_max)
    )
    t_bkg = time_all[np.asarray(bkg_sel, dtype=bool)]

    edges = np.arange(t0, t1 + bin_seconds, bin_seconds)
    counts, _ = np.histogram(t_bkg, bins=edges)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    rates = counts / float(bin_seconds)
    n_bins = len(rates)

    fallback_note = ""
    if method == "fixed":
        quiescent_rate = float(np.median(rates)) if n_bins else 0.0
        robust_sigma = float(1.4826 * np.median(np.abs(rates - quiescent_rate))) if n_bins else 0.0
        threshold = float(fixed_rate)
    elif method == "robust":
        keep = np.ones(n_bins, dtype=bool)
        for _ in range(max_iter):
            med = float(np.median(rates[keep]))
            mad = float(np.median(np.abs(rates[keep] - med)))
            sig = 1.4826 * mad
            if sig <= 0:
                s = float(np.std(rates[keep]))
                sig = s if s > 0 else max(med * 0.1, 1e-6)
            new_keep = rates <= (med + sigma_factor * sig)
            if new_keep.sum() < 3:
                break
            if np.array_equal(new_keep, keep):
                keep = new_keep
                break
            keep = new_keep
        quiescent_rate = float(np.median(rates[keep])) if keep.any() else 0.0
        mad_f = float(np.median(np.abs(rates[keep] - quiescent_rate))) if keep.any() else 0.0
        robust_sigma = 1.4826 * mad_f
        if robust_sigma <= 0:
            s = float(np.std(rates[keep])) if keep.any() else 0.0
            robust_sigma = s if s > 0 else max(quiescent_rate * 0.1, 1e-6)
        threshold = quiescent_rate + sigma_factor * robust_sigma
    else:
        raise ValueError("method debe ser 'robust' o 'fixed'")

    good_bins = rates <= threshold
    if int(good_bins.sum()) < min_gti_bins:
        fallback_note = (
            f"el umbral acepta {int(good_bins.sum())} bins (<{min_gti_bins}) -> salvaguarda con umbral fijo {fixed_rate} ct/s"
        )
        threshold = float(fixed_rate)
        good_bins = rates <= threshold

    bin_idx = np.floor((time_all - t0) / bin_seconds).astype(int)
    bin_idx = np.clip(bin_idx, 0, max(n_bins - 1, 0))
    gti_mask = good_bins[bin_idx] if n_bins else np.zeros_like(time_all, dtype=bool)

    sci_mask = (
        (ev_list["FLAG"] == 0)
        & (ev_list["PI"] >= sci_pi_min)
        & (ev_list["PI"] <= sci_pi_max)
        & (ev_list["PATTERN"] <= sci_pattern_max)
    )
    sci_mask = np.asarray(sci_mask, dtype=bool)
    scientific_before = int(sci_mask.sum())
    final_mask = sci_mask & gti_mask
    good_events = ev_list[final_mask]
    scientific_after = int(final_mask.sum())

    total_span = t1 - t0
    gti_exposure = float(int(good_bins.sum()) * bin_seconds)
    retained_fraction = (scientific_after / scientific_before) if scientific_before else 0.0

    diag = {
        "method": method,
        "bin_seconds": bin_seconds,
        "sigma_factor": sigma_factor,
        "min_gti_bins": min_gti_bins,
        "bkg_band": (bkg_pi_min, bkg_pi_max, bkg_pattern),
        "n_bins": int(n_bins),
        "n_bkg_events": int(t_bkg.size),
        "t0": t0,
        "t1": t1,
        "total_span": float(total_span),
        "threshold": float(threshold),
        "quiescent_rate": float(quiescent_rate),
        "robust_sigma": float(robust_sigma),
        "bin_centers": bin_centers,
        "good_bins": good_bins,
        "n_good_bins": int(good_bins.sum()),
        "n_rejected_bins": int((~good_bins).sum()) if n_bins else 0,
        "gti_exposure": gti_exposure,
        "events_total": int(len(ev_list)),
        "events_after_gti": int(gti_mask.sum()),
        "scientific_before": scientific_before,
        "scientific_after": scientific_after,
        "retained_fraction": float(retained_fraction),
        "fallback_note": fallback_note,
    }
    return good_events, diag


def apply_pi_exclusion(ev, pi_ranges):
    if not pi_ranges:
        return ev
    pi = np.asarray(ev["PI"], dtype=float)
    keep_mask = np.ones(len(pi), dtype=bool)
    for lo, hi in pi_ranges:
        keep_mask &= ~((pi >= lo) & (pi <= hi))
    return ev[keep_mask]


def build_region(members, sky_x, sky_y, min_radius_px, max_radius_px):
    pts = np.unique(np.column_stack([sky_x[members], sky_y[members]]), axis=0)
    cx, cy = np.median(pts[:, 0]), np.median(pts[:, 1])
    saturated = False
    degenerate = False
    if len(pts) < 3:
        r = min_radius_px
        degenerate = True
    else:
        d = np.hypot(pts[:, 0] - cx, pts[:, 1] - cy)
        r = np.percentile(d, 90)
        if r < min_radius_px:
            r = min_radius_px
            degenerate = True
        if r > max_radius_px:
            r = max_radius_px
            saturated = True
    return float(cx), float(cy), float(r), saturated, degenerate


def to_sky_circle(cx, cy, rpx, w, color):
    pix = CirclePixelRegion(center=PixCoord(x=cx, y=cy), radius=rpx)
    sky = pix.to_sky(w)
    sky.visual["color"] = color
    return sky


def validate_regions(sas_regions_sky, hdb_regions_pix, w, tolerance_arcsec):
    n_sas = len(sas_regions_sky)
    n_hdb = len(hdb_regions_pix)
    if n_sas == 0 or n_hdb == 0:
        return dict(TP=0, FP=n_hdb, FN=n_sas, precision=0.0, recall=0.0, F1=0.0, distances=[])

    sas = SkyCoord([r.center for r in sas_regions_sky])
    hcx = np.array([r.center.x for r in hdb_regions_pix])
    hcy = np.array([r.center.y for r in hdb_regions_pix])
    hsky = w.pixel_to_world(hcx, hcy)

    cost = np.zeros((n_sas, n_hdb))
    for i in range(n_sas):
        cost[i, :] = sas[i].separation(hsky).arcsec
    row, col = linear_sum_assignment(cost)
    dists = cost[row, col]
    matched = [float(d) for d in dists if d <= tolerance_arcsec]
    TP = int(sum(d <= tolerance_arcsec for d in dists))
    FP = n_hdb - TP
    FN = n_sas - TP
    P = TP / (TP + FP) if (TP + FP) else 0.0
    R = TP / (TP + FN) if (TP + FN) else 0.0
    F1 = 2 * P * R / (P + R) if (P + R) else 0.0
    return dict(TP=TP, FP=FP, FN=FN, precision=P, recall=R, F1=F1, distances=matched)



def write_product_csv(csv_path, row):
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def main():
    args = parse_args()

    profile = PROFILES[args.profile]
    TIMEWEIGHT = profile.timeweight
    HYPERPARAMETERPROFILE = profile.name
    MINCLUSTERSIZE = profile.minclustersize
    MINSAMPLES = profile.minsamples
    MINPERSISTENCE = profile.minpersistence
    MINMEANPROBABILITY = profile.minmeanprobability
    MINCLUSTEREVENTS = profile.minclusterevents

    obs = args.obs
    inst = "PN"

    files = autodetect_files(obs, args.dataset_dir, profile.product_code, profile.product_code_all)
    event_file = files["event_file"]
    region_list = files["region_list"]
    hdbscan_region_list = files["hdbscan_region_list"]
    hdbscan_region_list_all = files["hdbscan_region_list_all"]
    hdbscan_results_csv = files["hdbscan_results_csv"]

    print(f"Pipeline HDBSCAN {VERSION}")
    print(f"OBS={obs} INST={inst} EXPOSURE={files['exposure']}")
    print(f"event_file={event_file}")
    print(f"region_list={region_list}")

    with fits.open(event_file, mode="readonly", ignore_missing_end=True) as hdul:
        header = hdul[0].header
        ev_list = hdul["EVENTS"].data

    duration = header.get("DURATION", np.nan)
    object_name = header.get("OBJECT", obs)

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [header["REFXCRPX"], header["REFYCRPX"]]
    w.wcs.cdelt = np.array([header["REFXCDLT"], header["REFYCDLT"]])
    w.wcs.crval = [header["REFXCRVL"], header["REFYCRVL"]]
    w.wcs.ctype = [header["REFXCTYP"], header["REFYCTYP"]]
    w.wcs.radesys = header.get("RADECSYS", "ICRS")
    if "EQUINOX" in header:
        w.wcs.equinox = header["EQUINOX"]

    if APPLYBACKGROUNDGTI:
        good_events, gti_diag = filter_background_gti(
            ev_list,
            bin_seconds=BKGBINSECONDS,
            method=BKGMETHOD,
            fixed_rate=BKGFIXEDRATE,
            sigma_factor=BKGSIGMA,
            min_gti_bins=MINGTIBINS,
        )
        print(f"GTI method={gti_diag['method']} bin={gti_diag['bin_seconds']}s bins={gti_diag['n_bins']}")
        print(f"GTI quiescent={gti_diag['quiescent_rate']:.4f} ct/s sigma={gti_diag['robust_sigma']:.4f} threshold={gti_diag['threshold']:.4f}")
        print(f"Bins GTI aceptados/rechazados={gti_diag['n_good_bins']}/{gti_diag['n_rejected_bins']} exposición~={gti_diag['gti_exposure']:.1f} s")
        print(f"Científicos antes/después del GTI={gti_diag['scientific_before']}/{gti_diag['scientific_after']} retenida={gti_diag['retained_fraction']:.3f}")
        if gti_diag["fallback_note"]:
            print(f"GTI fallback: {gti_diag['fallback_note']}")
    else:
        sci = (
            (ev_list["FLAG"] == 0)
            & (ev_list["PI"] >= PIMIN)
            & (ev_list["PI"] <= PIMAX)
            & (ev_list["PATTERN"] <= 4)
        )
        good_events = ev_list[sci]

    cluster_events = apply_pi_exclusion(good_events, EXCLUDEPIRANGES)
    print(f"Eventos para clustering={len(cluster_events)} excluidos_por_PI={len(good_events) - len(cluster_events)}")

    event_ccd = np.asarray(cluster_events["CCDNR"], dtype=int)
    sky_x = np.asarray(cluster_events["X"], dtype=float)
    sky_y = np.asarray(cluster_events["Y"], dtype=float)
    event_time = np.asarray(cluster_events["TIME"], dtype=float)

    arcsec_per_pix = abs(header["REFXCDLT"]) * 3600.0
    psf_radius_pix = 6.0 / arcsec_per_pix
    min_radius_px = MINRADIUSARCSEC / arcsec_per_pix
    max_radius_px = MAXRADIUSARCSEC / arcsec_per_pix

    xs = sky_x / psf_radius_pix
    ys = sky_y / psf_radius_pix
    t_ref = event_time.min() if len(event_time) else 0.0
    ts = ((event_time - t_ref) / TIMESCALESECONDS) * TIMEWEIGHT if len(event_time) else np.array([])

    def features_for_idx(idx):
        if DETECTIONMODE == "spatiotemporal":
            return np.column_stack([xs[idx], ys[idx], ts[idx]])
        return np.column_stack([xs[idx], ys[idx]])

    def cluster_ccd(ccd):
        idx = np.where(event_ccd == ccd)[0]
        if len(idx) < MINCLUSTERSIZE:
            return idx, np.full(len(idx), -1, dtype=int), np.zeros(len(idx)), np.array([])
        clusterer = hdbscan.HDBSCAN(
            cluster_selection_method=CLUSTERSELECTIONMETHOD,
            metric=METRIC,
            min_cluster_size=MINCLUSTERSIZE,
            min_samples=MINSAMPLES,
            core_dist_n_jobs=NUMTHREADS,
        )
        clusterer.fit(features_for_idx(idx))
        return idx, clusterer.labels_, clusterer.probabilities_, clusterer.cluster_persistence_

    t0 = t.time()
    with Pool(NUMTHREADS) as pool:
        per_ccd_results = pool.map(cluster_ccd, range(1, 13))
    elapsed = t.time() - t0
    print(f"Clustering {DETECTIONMODE} en {elapsed:.1f} s")

    global_label = np.full(len(event_ccd), -1, dtype=int)
    global_prob = np.zeros(len(event_ccd))
    raw_clusters_per_ccd = []
    accepted_records = []
    all_records = []
    raw_global_id = 0
    global_id = 0
    n_saturated = 0

    for k, (idx, labels, probs, pers) in enumerate(per_ccd_results):
        ccd = k + 1
        n_here = int(labels.max() + 1) if len(labels) and labels.max() >= 0 else 0
        raw_clusters_per_ccd.append(n_here)
        for L in range(n_here):
            sel = labels == L
            members = idx[sel]
            nev = int(sel.sum())
            meanp = float(probs[sel].mean()) if nev else 0.0
            minp = float(probs[sel].min()) if nev else 0.0
            persistence = float(pers[L]) if L < len(pers) else 0.0
            cx, cy, r, sat, degen = build_region(members, sky_x, sky_y, min_radius_px, max_radius_px)
            rec = dict(
                ccd=ccd,
                local_label=L,
                raw_global_id=raw_global_id,
                members=members.copy(),
                nevents=nev,
                meanprob=meanp,
                minprob=minp,
                persistence=persistence,
                cx=cx,
                cy=cy,
                radiuspx=r,
                saturated=sat,
                degenerate=degen,
            )
            all_records.append(rec)
            raw_global_id += 1
            if (
                persistence >= MINPERSISTENCE
                and meanp >= MINMEANPROBABILITY
                and nev >= MINCLUSTEREVENTS
            ):
                rec_acc = dict(rec, global_id=global_id)
                accepted_records.append(rec_acc)
                global_label[members] = global_id
                global_prob[members] = probs[sel]
                global_id += 1
                if sat:
                    n_saturated += 1

    print(f"Clústeres crudos={len(all_records)} aceptados={len(accepted_records)} saturados_radio_máx={n_saturated}")
    print("CCD eventos ruido crudos aceptados")
    for k, (idx, labels, probs, pers) in enumerate(per_ccd_results):
        ccd = k + 1
        nev = len(idx)
        nnoise = int((labels == -1).sum()) if len(labels) else 0
        nraw = raw_clusters_per_ccd[k]
        nacc = sum(1 for r in accepted_records if r["ccd"] == ccd)
        print(f"{ccd:>3d} {nev:>7d} {nnoise:>7d} {nraw:>6d} {nacc:>9d}")
    print(f"Totales eventos={len(event_ccd)} ruido={int((global_label == -1).sum())} crudos={len(all_records)} aceptados={len(accepted_records)}")

    accepted_sky = [to_sky_circle(r["cx"], r["cy"], r["radiuspx"], w, "red") for r in accepted_records]
    all_sky = [to_sky_circle(r["cx"], r["cy"], r["radiuspx"], w, "orange") for r in all_records]

    if accepted_sky:
        Regions(accepted_sky).write(hdbscan_region_list, format="ds9", overwrite=True)
        print(f"Escrito el catálogo de aceptados {hdbscan_region_list} - {len(accepted_sky)}")
    if all_sky:
        Regions(all_sky).write(hdbscan_region_list_all, format="ds9", overwrite=True)
        print(f"Escrito el catálogo de todos {hdbscan_region_list_all} - {len(all_sky)}")

    if not region_list.exists():
        print(f"[WARN] No existe la lista de regiones SAS: {region_list}", file=sys.stderr)
        return

    sas_regions_radec = Regions.read(region_list, format="ds9")
    hdb_regions_pix = []
    if accepted_sky:
        hdb_regions_radec = Regions.read(hdbscan_region_list, format="ds9")
        for region in hdb_regions_radec:
            hdb_regions_pix.append(region.to_pixel(w))

    res = validate_regions(sas_regions_radec, hdb_regions_pix, w, MATCHTOLERANCEARCSEC)
    print(f"SAS={len(sas_regions_radec)} HDBSCAN={len(hdb_regions_pix)}")
    print(f"TP={res['TP']} FP={res['FP']} FN={res['FN']}")
    print(f"precisión={res['precision']:.3f} recall={res['recall']:.3f} F1={res['F1']:.3f}")
    if res["distances"]:
        dd = np.array(res["distances"])
        median_match = float(np.median(dd))
        max_match = float(dd.max())
        print(f"Distancias TP mediana={median_match:.2f} máx={max_match:.2f} arcsec")
    else:
        median_match = np.nan
        max_match = np.nan

    scientific_before = gti_diag["scientific_before"] if APPLYBACKGROUNDGTI else int(len(good_events))
    scientific_after = gti_diag["scientific_after"] if APPLYBACKGROUNDGTI else int(len(good_events))
    retained_fraction = gti_diag["retained_fraction"] if APPLYBACKGROUNDGTI else 1.0
    gti_threshold = gti_diag["threshold"] if APPLYBACKGROUNDGTI else np.nan
    gti_method = BKGMETHOD if APPLYBACKGROUNDGTI else ""
    gti_bin_seconds = BKGBINSECONDS if APPLYBACKGROUNDGTI else ""

    row = {
        "version": VERSION,
        "obs": obs,
        "object_name": object_name,
        "inst": inst,
        "exposure": files["exposure"],
        "event_file": str(event_file),
        "region_list": str(region_list),
        "productcode": profile.product_code,
        "productcodeall": profile.product_code_all,
        "hdbscan_region_list": str(hdbscan_region_list),
        "hdbscan_region_list_all": str(hdbscan_region_list_all),
        "detection_mode": DETECTIONMODE,
        "time_weight": TIMEWEIGHT,
        "time_scale_seconds": TIMESCALESECONDS,
        "hyperparameter_profile": HYPERPARAMETERPROFILE,
        "cluster_selection_method": CLUSTERSELECTIONMETHOD,
        "metric": METRIC,
        "min_cluster_size": MINCLUSTERSIZE,
        "min_samples": MINSAMPLES,
        "min_persistence": MINPERSISTENCE,
        "min_mean_probability": MINMEANPROBABILITY,
        "min_cluster_events": MINCLUSTEREVENTS,
        "min_radius_arcsec": MINRADIUSARCSEC,
        "max_radius_arcsec": MAXRADIUSARCSEC,
        "match_tolerance_arcsec": MATCHTOLERANCEARCSEC,
        "apply_background_gti": APPLYBACKGROUNDGTI,
        "gti_method": gti_method,
        "gti_bin_seconds": gti_bin_seconds,
        "gti_threshold": gti_threshold,
        "scientific_before": scientific_before,
        "scientific_after": scientific_after,
        "retained_fraction": retained_fraction,
        "events_for_clustering": int(len(cluster_events)),
        "raw_clusters": int(len(all_records)),
        "accepted_clusters": int(len(accepted_records)),
        "saturated_clusters": int(n_saturated),
        "sas_regions": int(len(sas_regions_radec)),
        "hdbscan_regions": int(len(hdb_regions_pix)),
        "TP": int(res["TP"]),
        "FP": int(res["FP"]),
        "FN": int(res["FN"]),
        "precision": float(res["precision"]),
        "recall": float(res["recall"]),
        "F1": float(res["F1"]),
        "median_match_arcsec": median_match,
        "max_match_arcsec": max_match,
        "duration_header": float(duration) if duration is not None else np.nan,
    }
    write_product_csv(hdbscan_results_csv, row)
    print(f"Escrito el resumen CSV {hdbscan_results_csv}")

    print(f"Objeto={object_name} DURACIÓN={duration}")


if __name__ == "__main__":
    main()
