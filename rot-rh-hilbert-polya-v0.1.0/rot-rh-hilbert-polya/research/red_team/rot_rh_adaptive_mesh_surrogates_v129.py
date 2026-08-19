#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v129 — Adaptive-Mesh Surrogate Prime-Cancellation Red-Team
=================================================================

PURPOSE
-------
v123 preregistered PASS found strong prime self-cancellation on the adaptive
crossing mesh. v124/v125 found extreme, partition-robust cross-scale Mangoldt
cancellation. v128 showed the corresponding PNT-error transform has an
extremely small bilinear correlation.

But the spectral cells E_k are themselves defined by crossings of the SAME
total phase. This creates a serious structural control question:

    Would a GENERIC oscillatory arithmetic phase, when allowed to rebuild
    ITS OWN crossing mesh, show the same self-cancellation?

v129 tests that directly.

SURROGATES
----------
Start from the frozen finite-L phase and multiply every arithmetic coefficient
in a dyadic log2(n) band by an independent random sign +/-1.

This preserves:
  * the exact frequency support log n;
  * coefficient magnitudes;
  * exponential regulator;
  * Gamma/AC/PNT-counterterm structure.

It destroys:
  * the true relative signs/alignment between arithmetic scales.

CRUCIALLY, for each surrogate we REBUILD ITS OWN ordered crossing support.
Thus this control preserves the adaptive-mesh mechanism rather than evaluating
randomized phases on the real mesh.

For fresh-depth cells k=1088..1407 the script compares:
  * prime cumulative suppression vs within-sequence permutations;
  * low-frequency suppression;
  * block suppression at 16,32,64;
  * exact dyadic cross-scale cumulative ratio.

If the real Mangoldt phase is not exceptional relative to these remeshed
surrogates, then the observed cancellation may be largely geometric/adaptive.
If the real phase beats them strongly, arithmetic scale alignment survives the
hardest mesh-adaptation control so far.

POST-HOC RED-TEAM. NOT A PREREGISTERED RH TEST.

NO Xi evaluations.
NO zeta evaluations.
NO known zero ordinates.
"""

from __future__ import annotations

import argparse
import copy
import csv
import dataclasses
import importlib
import json
import math
from pathlib import Path

import numpy as np

VERSION = "v129"
BUILD = "2026-08-18-adaptive-mesh-surrogate-prime-cancellation-v129"
EPS = 1e-30


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    if not rows:
        Path(path).write_text("", encoding="utf-8")
        return
    fields, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def cumulative_max(x):
    x = np.asarray(x, float)
    return float(np.max(np.abs(np.cumsum(x)))) if len(x) else 0.0


def lowfreq_fraction(x, bins):
    x = np.asarray(x, float)
    y = x - np.mean(x)
    z = np.fft.rfft(y)
    p = np.abs(z)**2
    if len(p) <= 1 or np.sum(p[1:]) <= EPS:
        return 0.0
    b = min(int(bins), len(p)-1)
    return float(np.sum(p[1:b+1]) / np.sum(p[1:]))


def rms(x):
    x = np.asarray(x, float)
    return float(np.sqrt(np.mean(x*x)))


def percentile_le(samples, value):
    a = np.asarray(samples, float)
    return float(np.mean(a <= value))


def permutation_metrics(x, nperm, rng, lowfreq_bins):
    x = np.asarray(x, float)
    real_c = cumulative_max(x)
    real_l = lowfreq_fraction(x, lowfreq_bins)

    pc = np.empty(int(nperm), float)
    pl = np.empty(int(nperm), float)
    for t in range(int(nperm)):
        y = rng.permutation(x)
        pc[t] = cumulative_max(y)
        pl[t] = lowfreq_fraction(y, lowfreq_bins)

    return {
        "cumulative_max": real_c,
        "lowfreq_fraction": real_l,
        "perm_cumulative_mean": float(np.mean(pc)),
        "cum_over_perm_mean": real_c / max(float(np.mean(pc)), EPS),
        "cumulative_percentile": percentile_le(pc, real_c),
        "lowfreq_percentile": percentile_le(pl, real_l),
    }


def block_sum_rms(x, block):
    x = np.asarray(x, float)
    b = int(block)
    n = len(x)//b
    if n < 1:
        return float("nan")
    z = x[:n*b].reshape(n, b).sum(axis=1)
    return float(np.sqrt(np.mean(z*z)))


def block_ratio(x, block, nperm, rng):
    real = block_sum_rms(x, block)
    vals = np.empty(int(nperm), float)
    for t in range(int(nperm)):
        vals[t] = block_sum_rms(rng.permutation(x), block)
    mean = float(np.mean(vals))
    return {
        "block": int(block),
        "real": real,
        "perm_mean": mean,
        "ratio": real / max(mean, EPS),
        "percentile": percentile_le(vals, real),
    }


def clone_with_phase_weights(phase, new_weights):
    """
    Robustly clone an imported phase object and replace phase_w.
    We deliberately fail loudly if replacement does not alter F.
    """
    new_weights = np.asarray(new_weights, float)

    # dataclass path first.
    if dataclasses.is_dataclass(phase):
        try:
            return dataclasses.replace(phase, phase_w=new_weights.copy())
        except Exception:
            pass

    # deepcopy + setattr.
    q = copy.deepcopy(phase)
    try:
        setattr(q, "phase_w", new_weights.copy())
        return q
    except Exception:
        pass

    # frozen/custom object.
    try:
        object.__setattr__(q, "phase_w", new_weights.copy())
        return q
    except Exception as exc:
        raise RuntimeError(
            "Could not clone phase with replacement phase_w. "
            "Inspect RenormalizedPhase implementation before continuing."
        ) from exc


def prime_endpoint_functions(Evals, phase_w, logs, chunk_terms):
    Evals = np.asarray(Evals, float)
    phase_w = np.asarray(phase_w, float)
    logs = np.asarray(logs, float)

    Fp = np.zeros(len(Evals), float)
    Anti = np.zeros(len(Evals), float)

    n = len(logs)
    for s in range(0, n, int(chunk_terms)):
        e = min(n, s+int(chunk_terms))
        lg = logs[s:e]
        ww = phase_w[s:e]
        ang = Evals[:, None]*lg[None, :]

        Fp += -np.sum(
            np.sin(ang)*ww[None, :], axis=1
        )/math.pi
        Anti += np.sum(
            np.cos(ang)*(ww/lg)[None, :], axis=1
        )/math.pi

    return Fp, Anti


def prime_cell_defects(roots, k_start, k_end, phase_w, logs, chunk_terms):
    ks = list(range(k_start, k_end+1))
    lo = np.asarray([roots[k-1] for k in ks], float)
    hi = np.asarray([roots[k] for k in ks], float)
    h = hi-lo
    endpoints = np.concatenate([[lo[0]], hi])

    Fp, Anti = prime_endpoint_functions(
        endpoints, phase_w, logs, chunk_terms
    )
    q = (
        Anti[1:]-Anti[:-1]
        - 0.5*h*(Fp[:-1]+Fp[1:])
    )
    return q, lo, hi


def dyadic_cross_ratio(q_total, lo, hi, phase_w, logs, chunk_terms):
    mid = 0.5*(lo+hi)
    h = hi-lo
    jband = np.floor(logs/math.log(2.0)+1e-12).astype(int)
    bands = sorted(np.unique(jband).tolist())
    qb = {j: np.zeros(len(lo), float) for j in bands}

    n = len(logs)
    for s in range(0, n, int(chunk_terms)):
        e = min(n, s+int(chunk_terms))
        om = logs[s:e]
        ww = phase_w[s:e]
        jb = jband[s:e]

        delta = 0.5*h[:, None]*om[None, :]
        ph = mid[:, None]*om[None, :]
        kernel = (
            (2.0/math.pi)
            * np.sin(ph)
            * (delta*np.cos(delta)-np.sin(delta))
        )
        contrib = kernel*(ww/om)[None, :]

        for j in np.unique(jb):
            mask = (jb == j)
            qb[int(j)] += np.sum(contrib[:, mask], axis=1)

    recon = np.sum(np.vstack([qb[j] for j in bands]), axis=0)
    rel = float(
        np.linalg.norm(recon-q_total) /
        max(np.linalg.norm(q_total), EPS)
    )
    denom = sum(cumulative_max(qb[j]) for j in bands)
    ratio = cumulative_max(q_total)/max(denom, EPS)
    return ratio, rel


def phase_probe(phase, E):
    x = np.asarray(E, float)
    return np.asarray(phase.F(x), float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--v102-module",
        default="rot_rh_pnt_finite_part_full_pipeline_v102_1",
    )
    ap.add_argument(
        "--base-module",
        default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit",
    )
    ap.add_argument(
        "--v123-prefix",
        default="rot_rh_blind_prime_self_cancellation_v123",
    )
    ap.add_argument("--L", type=int, default=16000)
    ap.add_argument("--tail-factor", type=float, default=18.0)
    ap.add_argument("--counter-q", type=int, default=8192)
    ap.add_argument("--counter-qcheck", type=int, default=16384)
    ap.add_argument("--scan-points", type=int, default=96)
    ap.add_argument("--k-start", type=int, default=1088)
    ap.add_argument("--k-end", type=int, default=1407)
    ap.add_argument("--max-depth", type=int, default=1408)
    ap.add_argument("--surrogates", type=int, default=6)
    ap.add_argument("--permutations", type=int, default=3000)
    ap.add_argument("--lowfreq-bins", type=int, default=8)
    ap.add_argument("--prime-chunk-terms", type=int, default=2048)
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument(
        "--prefix",
        default="rot_rh_adaptive_mesh_surrogates_v129",
    )
    args = ap.parse_args()

    cells_path = Path(
        str(Path(args.v123_prefix).expanduser().resolve()) + "_CELLS.csv"
    )
    if not cells_path.exists():
        raise FileNotFoundError(cells_path)

    raw = read_csv(cells_path)
    real_rows = sorted(
        [
            r for r in raw
            if args.k_start <= int(float(r["k"])) <= args.k_end
        ],
        key=lambda r: int(float(r["k"])),
    )
    if len(real_rows) != args.k_end-args.k_start+1:
        raise RuntimeError("Requested v123 real cell range incomplete.")

    qreal = np.asarray([float(r["q_prime"]) for r in real_rows], float)
    loreal = np.asarray([float(r["lo"]) for r in real_rows], float)
    hireal = np.asarray([float(r["hi"]) for r in real_rows], float)

    print("="*166)
    print("ROT-RH v129 — ADAPTIVE-MESH SURROGATE PRIME-CANCELLATION RED-TEAM")
    print("="*166)
    print("real source                     :", cells_path)
    print("cells                           :", f"{args.k_start}..{args.k_end}")
    print("surrogates                      :", args.surrogates)
    print("permutations/sequence           :", args.permutations)
    print("surrogate type                  : RANDOM_DYADIC_BAND_SIGNS + OWN REMESH")
    print("Xi evaluations                  : NONE")
    print("zeta evaluations                : NONE")
    print("known zero ordinates            : NONE")
    print("="*166)

    v102 = importlib.import_module(args.v102_module)
    base = importlib.import_module(args.base_module)

    print("[1] Build baseline finite-L phase", flush=True)
    phase = v102.RenormalizedPhase.build(
        base, int(args.L), float(args.tail_factor), int(args.counter_q)
    )
    w = np.asarray(phase.phase_w, float)
    logs = np.asarray(phase.logs, float)
    jband = np.floor(logs/math.log(2.0)+1e-12).astype(int)
    bands = sorted(np.unique(jband).tolist())

    # Validate real exact prime reconstruction from stored boundaries.
    endpoints_real = np.concatenate([[loreal[0]], hireal])
    Fpr, Ar = prime_endpoint_functions(
        endpoints_real, w, logs, args.prime_chunk_terms
    )
    hreal = hireal-loreal
    qreal_recon = (
        Ar[1:]-Ar[:-1]
        -0.5*hreal*(Fpr[:-1]+Fpr[1:])
    )
    real_recon_rel = float(
        np.linalg.norm(qreal_recon-qreal) /
        max(np.linalg.norm(qreal), EPS)
    )
    print("  arithmetic terms              :", len(logs))
    print("  dyadic bands                  :", len(bands))
    print("  real qprime reconstruction    :", f"{real_recon_rel:.3e}")

    rng = np.random.default_rng(args.seed)

    print("[2] Real-phase reference metrics")
    real_pm = permutation_metrics(
        qreal, args.permutations, rng, args.lowfreq_bins
    )
    real_b16 = block_ratio(qreal, 16, args.permutations, rng)
    real_b32 = block_ratio(qreal, 32, args.permutations, rng)
    real_b64 = block_ratio(qreal, 64, args.permutations, rng)
    real_cross, real_cross_recon = dyadic_cross_ratio(
        qreal, loreal, hireal, w, logs, args.prime_chunk_terms
    )

    print("  cumMax                         :", f"{real_pm['cumulative_max']:.6e}")
    print("  cum/permMean                   :", f"{real_pm['cum_over_perm_mean']:.6f}")
    print("  cumulative percentile          :", f"{100*real_pm['cumulative_percentile']:.4f}%")
    print("  lowF percentile                :", f"{100*real_pm['lowfreq_percentile']:.4f}%")
    print("  block16/32/64                  :",
          f"{real_b16['ratio']:.4f}, {real_b32['ratio']:.4f}, {real_b64['ratio']:.4f}")
    print("  dyadic cross-scale ratio       :", f"{real_cross:.6e}")
    print("  cross reconstruction relL2     :", f"{real_cross_recon:.3e}")

    if real_recon_rel > 1e-8 or real_cross_recon > 1e-8:
        raise RuntimeError("Real exact-kernel reconstruction failed.")

    print("[3] Adaptive-mesh randomized-band-sign surrogates", flush=True)
    surrogate_rows = []

    probe_E = np.linspace(loreal[0], hireal[-1], 9)
    F_base_probe = phase_probe(phase, probe_E)

    for si in range(args.surrogates):
        signs_by_band = {
            j: float(rng.choice(np.array([-1.0, 1.0])))
            for j in bands
        }

        # Avoid identity surrogate.
        if all(signs_by_band[j] > 0 for j in bands):
            signs_by_band[bands[0]] = -1.0

        signvec = np.asarray([signs_by_band[int(j)] for j in jband], float)
        ws = w*signvec
        surrogate = clone_with_phase_weights(phase, ws)

        # Verify mutation is scientifically active.
        F_sur_probe = phase_probe(surrogate, probe_E)
        phase_change = float(
            np.linalg.norm(F_sur_probe-F_base_probe)
        )
        stored_change = float(
            np.linalg.norm(np.asarray(surrogate.phase_w, float)-w)
        )
        if stored_change <= 1e-12 or phase_change <= 1e-10:
            raise RuntimeError(
                "Surrogate phase_w changed in storage but phase.F did not "
                "respond materially. RenormalizedPhase likely caches the prime "
                "phase elsewhere; inspect implementation before using v129."
            )

        print(
            f"  surrogate {si+1}/{args.surrogates}: "
            f"flippedBands={sum(signs_by_band[j]<0 for j in bands)} "
            f"phaseProbeDelta={phase_change:.3e}",
            flush=True,
        )

        support = v102.ordered_support(
            surrogate,
            depth=int(args.max_depth),
            scan_points=int(args.scan_points),
            qcheck=int(args.counter_qcheck),
        )
        roots = np.asarray(support.roots, float)

        qs, los, his = prime_cell_defects(
            roots,
            args.k_start,
            args.k_end,
            ws,
            logs,
            args.prime_chunk_terms,
        )

        pm = permutation_metrics(
            qs, args.permutations, rng, args.lowfreq_bins
        )
        b16 = block_ratio(qs, 16, args.permutations, rng)
        b32 = block_ratio(qs, 32, args.permutations, rng)
        b64 = block_ratio(qs, 64, args.permutations, rng)
        cross, cross_recon = dyadic_cross_ratio(
            qs, los, his, ws, logs, args.prime_chunk_terms
        )

        row = {
            "surrogate": si+1,
            "flipped_band_count":
                int(sum(signs_by_band[j] < 0 for j in bands)),
            "phase_probe_delta_l2": phase_change,
            "root_min": float(roots[0]),
            "root_max": float(roots[-1]),
            "root_max_residual": float(np.max(support.residuals)),
            "qphase": float(support.quadrature_phase_rel_l2),
            "qder": float(support.quadrature_deriv_rel_l2),
            "prime_rms": rms(qs),
            "prime_cumulative_max": pm["cumulative_max"],
            "cum_over_perm_mean": pm["cum_over_perm_mean"],
            "cumulative_percentile": pm["cumulative_percentile"],
            "lowfreq_percentile": pm["lowfreq_percentile"],
            "block16_ratio": b16["ratio"],
            "block32_ratio": b32["ratio"],
            "block64_ratio": b64["ratio"],
            "dyadic_cross_scale_ratio": cross,
            "cross_reconstruction_rel_l2": cross_recon,
        }
        surrogate_rows.append(row)

        print(
            f"    roots=[{roots[0]:.3f},{roots[-1]:.3f}] "
            f"res={row['root_max_residual']:.2e} "
            f"cum/perm={row['cum_over_perm_mean']:.3f} "
            f"pCum={100*row['cumulative_percentile']:.2f}% "
            f"pLowF={100*row['lowfreq_percentile']:.2f}% "
            f"cross={cross:.3e}"
        )

    print("[4] Real-vs-remeshed-surrogate ranking")
    fields = [
        ("cum_over_perm_mean", real_pm["cum_over_perm_mean"], "lower"),
        ("cumulative_percentile", real_pm["cumulative_percentile"], "lower"),
        ("lowfreq_percentile", real_pm["lowfreq_percentile"], "lower"),
        ("block16_ratio", real_b16["ratio"], "lower"),
        ("block32_ratio", real_b32["ratio"], "lower"),
        ("block64_ratio", real_b64["ratio"], "lower"),
        ("dyadic_cross_scale_ratio", real_cross, "lower"),
    ]

    rank_rows = []
    for key, rval, direction in fields:
        vals = np.asarray([float(r[key]) for r in surrogate_rows], float)
        if direction == "lower":
            frac = float(np.mean(vals <= rval))
        else:
            frac = float(np.mean(vals >= rval))
        row = {
            "metric": key,
            "real_value": rval,
            "surrogate_mean": float(np.mean(vals)),
            "surrogate_median": float(np.median(vals)),
            "surrogate_min": float(np.min(vals)),
            "surrogate_max": float(np.max(vals)),
            "fraction_surrogates_as_good_or_better": frac,
        }
        rank_rows.append(row)
        print(
            f"  {key:<28s}: real={rval:.4e} "
            f"surMed={np.median(vals):.4e} "
            f"asGood={100*frac:.2f}%"
        )

    # Conservative classification.
    qvalid = all(
        r["root_max_residual"] < 1e-8
        and r["cross_reconstruction_rel_l2"] < 1e-8
        for r in surrogate_rows
    )

    real_cross_rank = next(
        r["fraction_surrogates_as_good_or_better"]
        for r in rank_rows if r["metric"] == "dyadic_cross_scale_ratio"
    )
    real_cum_rank = next(
        r["fraction_surrogates_as_good_or_better"]
        for r in rank_rows if r["metric"] == "cum_over_perm_mean"
    )

    if not qvalid:
        mechanism = "SURROGATE_NUMERICS_INVALID__NO_SCIENTIFIC_INTERPRETATION"
    elif real_cross_rank <= 0.10 and real_cum_rank <= 0.10:
        mechanism = (
            "REAL_MANGOLDT_PHASE_BEATS_ADAPTIVE_REMESH_SURROGATES__"
            "ARITHMETIC_ALIGNMENT_REMAINS_SPECIAL"
        )
    elif real_cross_rank >= 0.30 and real_cum_rank >= 0.30:
        mechanism = (
            "ADAPTIVE_MESH_SURROGATES_REPRODUCE_CANCELLATION__"
            "MECHANISM_LARGELY_GEOMETRIC_NOT_ARITHMETIC"
        )
    else:
        mechanism = (
            "ADAPTIVE_MESH_CONTRIBUTES_BUT_REAL_ARITHMETIC_SPECIALNESS_MIXED"
        )

    prefix = Path(args.prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_csv(Path(str(prefix)+"_SURROGATES.csv"), surrogate_rows)
    write_csv(Path(str(prefix)+"_RANKING.csv"), rank_rows)

    packet = {
        "version": VERSION,
        "build": BUILD,
        "verdict": "PASS_V129_ADAPTIVE_MESH_SURROGATE_REDTEAM_COMPLETE"
            if qvalid else "FAIL_V129_SURROGATE_NUMERICAL_SUPPORT_INVALID",
        "mechanism": mechanism,
        "RH_proof": False,
        "Xi_used": False,
        "zeta_used": False,
        "known_zeros_used": False,
        "real_metrics": {
            "cum_over_perm_mean": real_pm["cum_over_perm_mean"],
            "cumulative_percentile": real_pm["cumulative_percentile"],
            "lowfreq_percentile": real_pm["lowfreq_percentile"],
            "block16_ratio": real_b16["ratio"],
            "block32_ratio": real_b32["ratio"],
            "block64_ratio": real_b64["ratio"],
            "dyadic_cross_scale_ratio": real_cross,
        },
        "surrogate_count": args.surrogates,
        "surrogate_rows": surrogate_rows,
        "rank_rows": rank_rows,
        "interpretation": (
            "This is the first control that allows randomized arithmetic "
            "coefficients to construct their own adaptive crossing meshes. "
            "If such surrogates reproduce the real suppression, the staircase "
            "mechanism is substantially geometric. If the real Mangoldt "
            "sequence remains an outlier, the arithmetic cross-scale alignment "
            "survives this adaptive-mesh confound."
        ),
    }

    out = Path(str(prefix)+"_VERDICT.json")
    out.write_text(json.dumps(packet, indent=2), encoding="utf-8")

    print("\n"+"="*166)
    print("FINAL v129 VERDICT")
    print("="*166)
    print("verdict                         :", packet["verdict"])
    print("mechanism                       :", mechanism)
    print("surrogate count                 :", args.surrogates)
    print("real cum/permMean               :", f"{real_pm['cum_over_perm_mean']:.12e}")
    print("real dyadic cross ratio         :", f"{real_cross:.12e}")
    print("surrogates as-good cumulative   :", f"{100*real_cum_rank:.6f}%")
    print("surrogates as-good cross-scale  :", f"{100*real_cross_rank:.6f}%")
    print("RH proof                        : FALSE")
    print("verdict file                    :", out)
    print("="*166)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
