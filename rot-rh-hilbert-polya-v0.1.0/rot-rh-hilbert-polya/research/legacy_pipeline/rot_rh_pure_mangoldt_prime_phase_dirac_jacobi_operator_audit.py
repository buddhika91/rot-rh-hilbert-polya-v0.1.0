#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH / PURE MANGOLDT PRIME-PHASE DIRAC-JACOBI OPERATOR AUDIT
==============================================================

Strict zero-information firewall.

Construction ingredients allowed before freeze
----------------------------------------------
* finite prime sieving and prime powers;
* von Mangoldt weights Lambda(n);
* the archimedean gamma phase (Riemann-Siegel theta);
* a fixed Abel-regularized prime-phase quantization law;
* fixed positive spectral weights for Jacobi realization.

Forbidden before freeze
-----------------------
* zeta zeros;
* Xi, xi, zeta, Xi'/Xi, zeta'/zeta;
* zero-fitted scale, shift, gain, sign, cutoff, depth, or index offset;
* zero-based model selection.

Free archimedean levels
-----------------------
Let

    N_bar(E) = theta(E)/pi + 1,

where

    theta(E) = Im log Gamma(1/4 + i E/2) - (E/2) log pi.

The kth free level E_k^(0) is the fixed solution

    N_bar(E_k^(0)) = k - 1/2.

Prime-phase scattering correction
---------------------------------
For a cutoff X, define epsilon_X = 1/log X and

    S_X(E) =
      -(1/pi) sum_{n <= X}
       [Lambda(n)/(log n * n^(1/2 + epsilon_X))]
       sin(E log n).

The pure arithmetic quantization law is

    N_bar(E) + S_X(E) = k - 1/2.

One root is selected inside the fixed kth free cell, bounded by midpoints of
neighboring archimedean levels.  If multiple roots occur, the root closest to
the free level is selected.  This rule is fixed before any zero scoring.

Primary cutoff
--------------
If --primary-cutoff=0, the primary cutoff is derived without zero information:

    X_primary = ceil((E_depth^(0))^2).

This "quadratic conductor cutoff" depends only on the requested finite operator
depth and the gamma phase.  Stability cutoffs are X/2 and 2X.

Self-adjoint operator realization
---------------------------------
The predicted positive levels are converted to a unique finite Jacobi matrix
using the Stieltjes recurrence with fixed boundary weights

    omega_k proportional to sin^2(k*pi/(d+1)).

The doubled chiral Dirac realization is

    H_D = [[0, J],
           [J, 0]],

whose positive eigenvalues equal the predicted prime-phase levels.

Scoring protocol
----------------
All operators, roots, coefficients, and hashes are written to disk first.
Only after the freeze manifest is complete is mpmath imported to obtain known
zero ordinates for scoring.

This is a numerical falsification audit, not a proof of RH.

Version: 2026-08-04-pure-mangoldt-prime-phase-dirac-jacobi-v1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.linalg import eigh_tridiagonal
    from scipy.optimize import brentq
    from scipy.special import digamma, loggamma
except Exception as exc:
    raise SystemExit(
        "scipy is required. Install with: pip install numpy scipy mpmath matplotlib\n"
        f"Import error: {exc}"
    )

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None


VERSION = "2026-08-04-pure-mangoldt-prime-phase-dirac-jacobi-v1"
EPS = 1.0e-14


# =============================================================================
# Data structures
# =============================================================================

@dataclass
class SpectralMetrics:
    count: int
    rmse: float
    mae: float
    relative_rmse: float
    maximum_absolute_error: float
    spacing_rmse: float
    spacing_relative_rmse: float
    mean_signed_error: float


@dataclass
class QuantizedSpectrum:
    family: str
    cutoff: int
    depth: int
    roots: np.ndarray
    root_residuals: np.ndarray
    root_counts_per_cell: np.ndarray
    free_levels: np.ndarray
    diagonal: np.ndarray
    offdiagonal: np.ndarray
    jacobi_eigenvalues: np.ndarray
    jacobi_self_adjoint_defect: float
    jacobi_minimum_offdiagonal: float
    jacobi_orthogonality_defect: float
    jacobi_eigenvalue_preservation_error: float
    dirac_chiral_defect: float
    operator_hash: str


# =============================================================================
# Generic utilities
# =============================================================================

def json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def sha256_arrays(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        z = np.ascontiguousarray(array)
        digest.update(str(z.dtype).encode("ascii"))
        digest.update(np.asarray(z.shape, np.int64).tobytes())
        digest.update(z.tobytes())
    return digest.hexdigest()


def parse_families(text: str) -> List[str]:
    values = [part.strip() for part in str(text).split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected comma-separated families.")
    return values


def stable_positive_weights(depth: int) -> np.ndarray:
    indices = np.arange(1, depth + 1, dtype=float)
    weights = np.sin(indices * math.pi / (depth + 1.0)) ** 2
    weights /= float(np.sum(weights))
    return weights


# =============================================================================
# Archimedean phase, no zeta or zero access
# =============================================================================

def theta(t: float) -> float:
    z = 0.25 + 0.5j * float(t)
    return float(np.imag(loggamma(z)) - 0.5 * float(t) * math.log(math.pi))


def theta_prime(t: float) -> float:
    z = 0.25 + 0.5j * float(t)
    return float(
        0.5 * np.real(digamma(z)) - 0.5 * math.log(math.pi)
    )


def smooth_count(t: float) -> float:
    return theta(float(t)) / math.pi + 1.0


def free_archimedean_levels(depth: int) -> np.ndarray:
    roots: List[float] = []
    lower_monotone = 7.0
    for k in range(1, depth + 1):
        target = float(k) - 0.5
        upper = max(20.0, 8.0 * float(k) + 20.0)
        while smooth_count(upper) < target:
            upper *= 1.5
            if upper > 1.0e7:
                raise RuntimeError("Failed to bracket archimedean level.")
        root = brentq(
            lambda t: smooth_count(t) - target,
            lower_monotone,
            upper,
            xtol=1e-13,
            rtol=1e-13,
            maxiter=500,
        )
        roots.append(float(root))
    return np.asarray(roots, float)


# =============================================================================
# Prime arithmetic
# =============================================================================

def sieve_primes(limit: int) -> np.ndarray:
    if limit < 2:
        return np.empty(0, dtype=np.int64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(math.isqrt(limit)) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = False
    return np.flatnonzero(sieve).astype(np.int64)


def mangoldt_terms(limit: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    numbers: List[int] = []
    lambdas: List[float] = []
    exponents: List[int] = []
    for p_raw in sieve_primes(limit):
        p = int(p_raw)
        lp = math.log(float(p))
        q = p
        exponent = 1
        while q <= limit:
            numbers.append(q)
            lambdas.append(lp)
            exponents.append(exponent)
            if q > limit // p:
                break
            q *= p
            exponent += 1
    order = np.argsort(numbers)
    return (
        np.asarray(numbers, np.int64)[order],
        np.asarray(lambdas, float)[order],
        np.asarray(exponents, np.int16)[order],
    )


def phase_terms(
    cutoff: int,
    family: str,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    if cutoff < 8:
        raise ValueError("Cutoff must be at least 8.")

    if family == "archimedean_only":
        return (
            np.empty(0, float),
            np.empty(0, float),
            np.empty(0, float),
            {
                "family": family,
                "term_count": 0,
                "epsilon": 0.0,
            },
        )

    if family == "smooth_integer":
        numbers = np.arange(2, cutoff + 1, dtype=np.int64)
        lambdas = np.ones(len(numbers), dtype=float)
        exponents = np.zeros(len(numbers), dtype=np.int16)
    else:
        numbers, lambdas, exponents = mangoldt_terms(cutoff)

    epsilon = 1.0 / math.log(float(cutoff))
    logs = np.log(numbers.astype(float))
    weights = lambdas / (
        logs * numbers.astype(float) ** (0.5 + epsilon)
    )

    if family == "prime_only":
        keep = exponents == 1
        logs = logs[keep]
        weights = weights[keep]

    elif family == "prime_power_only":
        keep = exponents >= 2
        logs = logs[keep]
        weights = weights[keep]

    rng = np.random.default_rng(int(seed))
    phases = np.zeros(len(logs), dtype=float)

    if family == "signal":
        pass
    elif family == "signflip":
        weights = -weights
    elif family == "permuted_weights":
        weights = weights.copy()
        rng.shuffle(weights)
    elif family == "phase_scramble":
        phases = rng.uniform(0.0, 2.0 * math.pi, len(logs))
    elif family in {
        "archimedean_only",
        "smooth_integer",
        "prime_only",
        "prime_power_only",
    }:
        pass
    else:
        raise ValueError(f"Unknown family {family!r}")

    metadata = {
        "family": family,
        "term_count": int(len(logs)),
        "epsilon": float(epsilon),
        "weight_l1": float(np.sum(np.abs(weights))),
        "weight_l2": float(np.linalg.norm(weights)),
        "minimum_log": float(np.min(logs)) if len(logs) else None,
        "maximum_log": float(np.max(logs)) if len(logs) else None,
    }
    return logs, weights, phases, metadata


# =============================================================================
# Prime-phase quantization
# =============================================================================

def prime_phase_correction(
    energies: np.ndarray,
    logs: np.ndarray,
    weights: np.ndarray,
    phases: np.ndarray,
) -> np.ndarray:
    energies = np.atleast_1d(np.asarray(energies, float))
    if len(logs) == 0:
        return np.zeros_like(energies)
    values = np.sin(np.outer(energies, logs) + phases[None, :])
    return -(values @ weights) / math.pi


def quantization_values(
    energies: np.ndarray,
    target: float,
    logs: np.ndarray,
    weights: np.ndarray,
    phases: np.ndarray,
) -> np.ndarray:
    energies = np.atleast_1d(np.asarray(energies, float))
    smooth = np.asarray([smooth_count(value) for value in energies])
    return (
        smooth
        + prime_phase_correction(energies, logs, weights, phases)
        - float(target)
    )


def quantize_levels(
    depth: int,
    cutoff: int,
    family: str,
    seed: int,
    grid_points: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    free_extended = free_archimedean_levels(depth + 1)
    free = free_extended[:depth]

    if family == "archimedean_only":
        return (
            free.copy(),
            np.zeros(depth),
            np.ones(depth, dtype=int),
            free,
            {
                "family": family,
                "term_count": 0,
                "epsilon": 0.0,
            },
        )

    logs, weights, phases, metadata = phase_terms(cutoff, family, seed)

    roots = np.zeros(depth, dtype=float)
    residuals = np.zeros(depth, dtype=float)
    root_counts = np.zeros(depth, dtype=int)

    for index in range(depth):
        k = index + 1
        target = float(k) - 0.5

        if index == 0:
            lower = 7.0
        else:
            lower = 0.5 * (free[index - 1] + free[index])

        if index + 1 < depth:
            upper = 0.5 * (free[index] + free[index + 1])
        else:
            upper = free[index] + 0.5 * (free[index] - free[index - 1])

        grid = np.linspace(lower, upper, int(grid_points))
        values = quantization_values(
            grid,
            target,
            logs,
            weights,
            phases,
        )
        crossings = np.flatnonzero(values[:-1] * values[1:] <= 0.0)
        root_counts[index] = int(len(crossings))

        if len(crossings) == 0:
            best = int(np.argmin(np.abs(values)))
            roots[index] = float(grid[best])
            residuals[index] = float(abs(values[best]))
            continue

        centers = 0.5 * (grid[crossings] + grid[crossings + 1])
        selected = int(
            crossings[np.argmin(np.abs(centers - free[index]))]
        )
        root = brentq(
            lambda energy: float(
                quantization_values(
                    np.asarray([energy]),
                    target,
                    logs,
                    weights,
                    phases,
                )[0]
            ),
            float(grid[selected]),
            float(grid[selected + 1]),
            xtol=1e-13,
            rtol=1e-13,
            maxiter=500,
        )
        roots[index] = float(root)
        residuals[index] = float(
            abs(
                quantization_values(
                    np.asarray([root]),
                    target,
                    logs,
                    weights,
                    phases,
                )[0]
            )
        )

    return roots, residuals, root_counts, free, metadata


# =============================================================================
# Jacobi and Dirac realization
# =============================================================================

def jacobi_from_atoms(
    atoms: np.ndarray,
    weights: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    x = np.asarray(atoms, float)
    w = np.asarray(weights, float)
    w = np.maximum(w, EPS)
    w /= float(np.sum(w))
    depth = len(x)

    q_previous = np.zeros(depth, float)
    q = np.ones(depth, float)
    q /= math.sqrt(float(np.sum(w * q * q)))
    basis: List[np.ndarray] = []

    diagonal = np.zeros(depth, float)
    offdiagonal = np.zeros(max(depth - 1, 0), float)
    beta_previous = 0.0

    for k in range(depth):
        basis.append(q.copy())
        alpha = float(np.sum(w * x * q * q))
        diagonal[k] = alpha
        residual = x * q - alpha * q - beta_previous * q_previous

        for _ in range(2):
            for vector in basis:
                residual -= float(np.sum(w * residual * vector)) * vector

        if k + 1 < depth:
            beta = math.sqrt(max(float(np.sum(w * residual * residual)), 0.0))
            if beta <= EPS or not math.isfinite(beta):
                raise RuntimeError(
                    f"Jacobi recurrence broke at k={k}, beta={beta!r}."
                )
            offdiagonal[k] = beta
            q_previous, q = q, residual / beta
            beta_previous = beta

    gram = np.empty((depth, depth), float)
    for i in range(depth):
        for j in range(depth):
            gram[i, j] = float(np.sum(w * basis[i] * basis[j]))

    eigenvalues = eigh_tridiagonal(
        diagonal,
        offdiagonal,
        eigvals_only=True,
    )
    diagnostics = {
        "orthogonality_defect": float(
            np.linalg.norm(gram - np.eye(depth), "fro")
        ),
        "minimum_offdiagonal": float(np.min(offdiagonal)),
        "eigenvalue_preservation_error": float(
            np.max(np.abs(np.sort(x) - eigenvalues))
        ),
    }
    return diagonal, offdiagonal, diagnostics


def build_quantized_operator(
    family: str,
    cutoff: int,
    depth: int,
    seed: int,
    grid_points: int,
) -> QuantizedSpectrum:
    roots, residuals, root_counts, free, _metadata = quantize_levels(
        depth=depth,
        cutoff=cutoff,
        family=family,
        seed=seed,
        grid_points=grid_points,
    )

    fixed_weights = stable_positive_weights(depth)
    diagonal, offdiagonal, diagnostics = jacobi_from_atoms(
        roots,
        fixed_weights,
    )
    eigenvalues = eigh_tridiagonal(
        diagonal,
        offdiagonal,
        eigvals_only=True,
    )

    jacobi = np.diag(diagonal)
    if depth > 1:
        jacobi += np.diag(offdiagonal, 1) + np.diag(offdiagonal, -1)

    dirac = np.block(
        [
            [np.zeros_like(jacobi), jacobi],
            [jacobi, np.zeros_like(jacobi)],
        ]
    )
    chirality = np.diag(
        np.concatenate([np.ones(depth), -np.ones(depth)])
    )
    chiral_defect = float(
        np.linalg.norm(chirality @ dirac + dirac @ chirality, "fro")
    )

    return QuantizedSpectrum(
        family=family,
        cutoff=int(cutoff),
        depth=int(depth),
        roots=roots,
        root_residuals=residuals,
        root_counts_per_cell=root_counts,
        free_levels=free,
        diagonal=diagonal,
        offdiagonal=offdiagonal,
        jacobi_eigenvalues=eigenvalues,
        jacobi_self_adjoint_defect=float(
            np.linalg.norm(jacobi - jacobi.T, "fro")
        ),
        jacobi_minimum_offdiagonal=float(
            diagnostics["minimum_offdiagonal"]
        ),
        jacobi_orthogonality_defect=float(
            diagnostics["orthogonality_defect"]
        ),
        jacobi_eigenvalue_preservation_error=float(
            diagnostics["eigenvalue_preservation_error"]
        ),
        dirac_chiral_defect=chiral_defect,
        operator_hash=sha256_arrays(
            roots,
            diagonal,
            offdiagonal,
        ),
    )


# =============================================================================
# Scoring after freeze only
# =============================================================================

def load_zero_ordinates(count: int, dps: int) -> np.ndarray:
    import mpmath as mp

    mp.mp.dps = int(dps)
    return np.asarray(
        [
            float(mp.im(mp.zetazero(index)))
            for index in range(1, count + 1)
        ],
        float,
    )


def spectral_metrics(
    target: np.ndarray,
    predicted: np.ndarray,
) -> SpectralMetrics:
    target = np.asarray(target, float)
    predicted = np.asarray(predicted, float)
    error = predicted - target
    rmse = float(np.sqrt(np.mean(error * error)))
    target_scale = max(float(np.sqrt(np.mean(target * target))), EPS)

    if len(target) >= 2:
        target_spacing = np.diff(target)
        predicted_spacing = np.diff(predicted)
        spacing_error = predicted_spacing - target_spacing
        spacing_rmse = float(
            np.sqrt(np.mean(spacing_error * spacing_error))
        )
        spacing_scale = max(
            float(np.sqrt(np.mean(target_spacing * target_spacing))),
            EPS,
        )
    else:
        spacing_rmse = math.nan
        spacing_scale = math.nan

    return SpectralMetrics(
        count=int(len(target)),
        rmse=rmse,
        mae=float(np.mean(np.abs(error))),
        relative_rmse=rmse / target_scale,
        maximum_absolute_error=float(np.max(np.abs(error))),
        spacing_rmse=spacing_rmse,
        spacing_relative_rmse=(
            spacing_rmse / spacing_scale
            if math.isfinite(spacing_rmse)
            else math.nan
        ),
        mean_signed_error=float(np.mean(error)),
    )


# =============================================================================
# Main
# =============================================================================

def run(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    families = parse_families(args.families)
    required = {
        "signal",
        "archimedean_only",
        "signflip",
        "permuted_weights",
        "phase_scramble",
        "smooth_integer",
        "prime_only",
        "prime_power_only",
    }
    missing = sorted(required - set(families))
    if missing:
        raise ValueError(
            "The preregistered audit requires these families: "
            + ", ".join(missing)
        )

    free_for_cutoff = free_archimedean_levels(args.operator_depth)
    if args.primary_cutoff <= 0:
        primary_cutoff = int(
            math.ceil(float(free_for_cutoff[-1]) ** 2)
        )
    else:
        primary_cutoff = int(args.primary_cutoff)

    stability_cutoffs = sorted(
        {
            max(8, int(round(primary_cutoff / 2))),
            primary_cutoff,
            int(round(primary_cutoff * 2)),
        }
    )

    prefix = Path(args.out_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": Path(str(prefix) + "_summary.json"),
        "operators": Path(str(prefix) + "_operators.npz"),
        "roots": Path(str(prefix) + "_quantized_roots.csv"),
        "coefficients": Path(str(prefix) + "_jacobi_coefficients.csv"),
        "controls": Path(str(prefix) + "_controls.csv"),
        "predictions": Path(str(prefix) + "_zero_predictions.csv"),
        "stability": Path(str(prefix) + "_cutoff_stability.csv"),
        "manifest": Path(str(prefix) + "_freeze_manifest.json"),
        "plot": Path(str(prefix) + "_spectrum.png"),
        "report": Path(str(prefix) + "_report.md"),
    }

    print("=" * 164)
    print("ROT-RH / PURE MANGOLDT PRIME-PHASE DIRAC-JACOBI OPERATOR AUDIT")
    print("=" * 164)
    print(f"version                    : {VERSION}")
    print(f"operator depth             : {args.operator_depth}")
    print(f"score zeros                : {args.zero_count}")
    print(f"primary cutoff             : {primary_cutoff}")
    print(f"cutoff derivation          : ceil((free gamma level d)^2)")
    print(f"stability cutoffs          : {stability_cutoffs}")
    print(f"families                   : {families}")
    print("construction inputs        : Gamma phase + primes + Lambda(n)")
    print("Xi/zeta access             : FORBIDDEN")
    print("known zero access          : FORBIDDEN UNTIL AFTER FREEZE")
    print("scale/shift/index fitting  : NONE")
    print("=" * 164)

    operators: Dict[Tuple[str, int], QuantizedSpectrum] = {}

    # Signal at three cutoff scales.
    for cutoff in stability_cutoffs:
        record = build_quantized_operator(
            family="signal",
            cutoff=cutoff,
            depth=args.operator_depth,
            seed=args.seed + cutoff,
            grid_points=args.grid_points,
        )
        operators[("signal", cutoff)] = record
        print(
            f"[construct] family=signal cutoff={cutoff:<7} "
            f"max_resid={np.max(record.root_residuals):.3e} "
            f"min_b={record.jacobi_minimum_offdiagonal:+.3e} "
            f"range=[{record.roots[0]:.6f},{record.roots[-1]:.6f}]"
        )

    # Prespecified controls at primary cutoff.
    for family in families:
        if family == "signal":
            continue
        record = build_quantized_operator(
            family=family,
            cutoff=primary_cutoff,
            depth=args.operator_depth,
            seed=args.seed + 1009 * (families.index(family) + 1),
            grid_points=args.grid_points,
        )
        operators[(family, primary_cutoff)] = record
        print(
            f"[construct] family={family:<20} cutoff={primary_cutoff:<7} "
            f"max_resid={np.max(record.root_residuals):.3e} "
            f"min_b={record.jacobi_minimum_offdiagonal:+.3e} "
            f"range=[{record.roots[0]:.6f},{record.roots[-1]:.6f}]"
        )

    primary = operators[("signal", primary_cutoff)]

    # -------------------------------------------------------------------------
    # Freeze all operators BEFORE importing mpmath or reading zeros.
    # -------------------------------------------------------------------------
    archive_payload: Dict[str, np.ndarray] = {}
    root_rows: List[Dict[str, Any]] = []
    coefficient_rows: List[Dict[str, Any]] = []

    for (family, cutoff), record in sorted(operators.items()):
        token = f"{family}_X{cutoff}_d{record.depth}"
        archive_payload[token + "_roots"] = record.roots
        archive_payload[token + "_free_levels"] = record.free_levels
        archive_payload[token + "_a"] = record.diagonal
        archive_payload[token + "_b"] = record.offdiagonal

        jacobi = np.diag(record.diagonal)
        jacobi += (
            np.diag(record.offdiagonal, 1)
            + np.diag(record.offdiagonal, -1)
        )
        dirac = np.block(
            [
                [np.zeros_like(jacobi), jacobi],
                [jacobi, np.zeros_like(jacobi)],
            ]
        )
        archive_payload[token + "_jacobi"] = jacobi
        archive_payload[token + "_dirac"] = dirac

        for index in range(record.depth):
            root_rows.append(
                {
                    "family": family,
                    "cutoff": cutoff,
                    "index": index + 1,
                    "free_level": float(record.free_levels[index]),
                    "quantized_level": float(record.roots[index]),
                    "prime_phase_shift": float(
                        record.roots[index] - record.free_levels[index]
                    ),
                    "quantization_residual": float(
                        record.root_residuals[index]
                    ),
                    "roots_in_free_cell": int(
                        record.root_counts_per_cell[index]
                    ),
                    "operator_hash": record.operator_hash,
                }
            )
            coefficient_rows.append(
                {
                    "family": family,
                    "cutoff": cutoff,
                    "index": index,
                    "a_diagonal": float(record.diagonal[index]),
                    "b_to_next": (
                        float(record.offdiagonal[index])
                        if index + 1 < record.depth
                        else math.nan
                    ),
                    "operator_hash": record.operator_hash,
                }
            )

    np.savez_compressed(paths["operators"], **archive_payload)

    with paths["roots"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(root_rows[0]))
        writer.writeheader()
        writer.writerows(root_rows)

    with paths["coefficients"].open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(coefficient_rows[0]),
        )
        writer.writeheader()
        writer.writerows(coefficient_rows)

    freeze_manifest = {
        "version": VERSION,
        "primary_family": "signal",
        "primary_cutoff": primary_cutoff,
        "operator_depth": args.operator_depth,
        "operator_hash": primary.operator_hash,
        "operator_file_sha256": hashlib.sha256(
            paths["operators"].read_bytes()
        ).hexdigest(),
        "roots_file_sha256": hashlib.sha256(
            paths["roots"].read_bytes()
        ).hexdigest(),
        "coefficients_file_sha256": hashlib.sha256(
            paths["coefficients"].read_bytes()
        ).hexdigest(),
        "firewall": {
            "known_zeros_used_in_construction": False,
            "xi_used": False,
            "zeta_used": False,
            "rvm_zero_ordinates_used": False,
            "zero_fitted_scale": False,
            "zero_fitted_shift": False,
            "zero_fitted_sign": False,
            "zero_fitted_cutoff": False,
            "zero_fitted_depth": False,
            "zero_fitted_index_offset": False,
        },
        "primary_cutoff_rule": (
            "ceil((operator-depth archimedean gamma level)^2)"
            if args.primary_cutoff <= 0
            else "explicit user-supplied cutoff"
        ),
        "primary_cutoff_inputs": {
            "free_terminal_level": float(free_for_cutoff[-1]),
            "computed_cutoff": primary_cutoff,
        },
        "frozen_at_unix": time.time(),
    }
    paths["manifest"].write_text(
        json.dumps(json_safe(freeze_manifest), indent=2),
        encoding="utf-8",
    )

    print("-" * 164)
    print("OPERATOR FREEZE COMPLETE")
    print(f"primary hash        : {primary.operator_hash}")
    print(
        f"operator archive sha: {freeze_manifest['operator_file_sha256']}"
    )
    print("known zero access now enabled for final scoring only")
    print("-" * 164)

    # -------------------------------------------------------------------------
    # Scoring only after freeze.
    # -------------------------------------------------------------------------
    target = load_zero_ordinates(args.zero_count, args.zero_dps)

    controls: List[Dict[str, Any]] = []
    predictions: List[Dict[str, Any]] = []
    metrics_by_family: Dict[str, SpectralMetrics] = {}

    for family in families:
        record = operators[(family, primary_cutoff)]
        predicted = record.roots[: args.zero_count]
        metrics = spectral_metrics(target, predicted)
        metrics_by_family[family] = metrics
        controls.append(
            {
                "family": family,
                "cutoff": primary_cutoff,
                "operator_depth": args.operator_depth,
                **asdict(metrics),
                "maximum_quantization_residual": float(
                    np.max(record.root_residuals)
                ),
                "minimum_jacobi_offdiagonal": (
                    record.jacobi_minimum_offdiagonal
                ),
                "operator_hash": record.operator_hash,
            }
        )
        for index, (truth, estimate) in enumerate(
            zip(target, predicted), start=1
        ):
            predictions.append(
                {
                    "family": family,
                    "zero_index": index,
                    "target_gamma": float(truth),
                    "predicted_energy": float(estimate),
                    "error": float(estimate - truth),
                    "absolute_error": float(abs(estimate - truth)),
                }
            )

    controls.sort(key=lambda row: (row["rmse"], row["family"]))

    with paths["controls"].open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(controls[0]))
        writer.writeheader()
        writer.writerows(controls)

    with paths["predictions"].open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(predictions[0]),
        )
        writer.writeheader()
        writer.writerows(predictions)

    stability_rows: List[Dict[str, Any]] = []
    for cutoff in stability_cutoffs:
        if cutoff == primary_cutoff:
            continue
        other = operators[("signal", cutoff)]
        difference = (
            other.roots[: args.zero_count]
            - primary.roots[: args.zero_count]
        )
        stability_rows.append(
            {
                "primary_cutoff": primary_cutoff,
                "other_cutoff": cutoff,
                "count": args.zero_count,
                "rmse_drift": float(
                    np.sqrt(np.mean(difference * difference))
                ),
                "relative_l2_drift": float(
                    np.linalg.norm(difference)
                    / max(
                        np.linalg.norm(
                            primary.roots[: args.zero_count]
                        ),
                        EPS,
                    )
                ),
                "maximum_absolute_drift": float(
                    np.max(np.abs(difference))
                ),
            }
        )

    with paths["stability"].open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(stability_rows[0]),
        )
        writer.writeheader()
        writer.writerows(stability_rows)

    signal_metrics = metrics_by_family["signal"]
    archimedean_metrics = metrics_by_family["archimedean_only"]
    non_signal = [
        row
        for row in controls
        if row["family"] not in {"signal", "archimedean_only"}
    ]
    best_disruption = min(non_signal, key=lambda row: row["rmse"])
    maximum_relative_cutoff_drift = max(
        float(row["relative_l2_drift"])
        for row in stability_rows
    )

    structural_gates = {
        "frozen_before_zero_access": True,
        "zero_xi_zeta_free_construction": True,
        "all_primary_cells_have_roots": bool(
            np.all(primary.root_counts_per_cell >= 1)
        ),
        "maximum_quantization_residual": bool(
            np.max(primary.root_residuals)
            <= args.maximum_quantization_residual
        ),
        "jacobi_self_adjoint": bool(
            primary.jacobi_self_adjoint_defect
            <= args.maximum_self_adjoint_defect
        ),
        "positive_jacobi_offdiagonals": bool(
            primary.jacobi_minimum_offdiagonal
            > args.minimum_positive_offdiagonal
        ),
        "jacobi_spectrum_preserved": bool(
            primary.jacobi_eigenvalue_preservation_error
            <= args.maximum_eigenvalue_preservation_error
        ),
        "dirac_chiral_symmetry": bool(
            primary.dirac_chiral_defect
            <= args.maximum_chiral_defect
        ),
        "cutoff_stability": bool(
            maximum_relative_cutoff_drift
            <= args.maximum_relative_cutoff_drift
        ),
    }

    performance_gates = {
        "maximum_direct_relative_rmse": bool(
            signal_metrics.relative_rmse
            <= args.maximum_relative_rmse
        ),
        "maximum_spacing_relative_rmse": bool(
            signal_metrics.spacing_relative_rmse
            <= args.maximum_spacing_relative_rmse
        ),
        "beats_archimedean_free_operator": bool(
            signal_metrics.rmse + args.minimum_control_margin
            <= archimedean_metrics.rmse
        ),
        "beats_best_disruption_control": bool(
            signal_metrics.rmse + args.minimum_control_margin
            <= float(best_disruption["rmse"])
        ),
    }

    gates = {**structural_gates, **performance_gates}
    failed = [name for name, passed in gates.items() if not passed]

    if all(gates.values()):
        verdict = (
            "PASS_PURE_MANGOLDT_PRIME_PHASE_DIRAC_JACOBI_"
            "DIRECT_ZERO_SPECTRUM"
        )
        next_step = (
            "Freeze this exact law. Extend depth and conductor-derived cutoff "
            "without using zero errors, then test whether the predicted levels "
            "and Jacobi prefixes converge. The proof bottleneck becomes deriving "
            "the prime-phase quantization law and its infinite self-adjoint limit."
        )
    elif all(structural_gates.values()):
        verdict = (
            "FAIL_PURE_MANGOLDT_PRIME_PHASE_OPERATOR_ZERO_SPECTRUM"
        )
        next_step = (
            "The zero-free prime-phase operator is structurally valid but does "
            "not reproduce the zeros. Do not tune damping, sign, cutoff, or cell "
            "selection using zero errors; derive a different boundary/scattering law."
        )
    else:
        verdict = (
            "FAIL_PURE_MANGOLDT_PRIME_PHASE_OPERATOR_STRUCTURAL_GATES"
        )
        next_step = (
            "Repair root uniqueness, cutoff stability, or Jacobi realization "
            "before interpreting the zero-spectrum result."
        )

    if plt is not None:
        indices = np.arange(1, args.zero_count + 1)
        plt.figure(figsize=(10, 6))
        plt.plot(indices, target, marker="o", label="Known zeros (score only)")
        plt.plot(
            indices,
            primary.roots[: args.zero_count],
            marker="s",
            label="Pure Mangoldt prime-phase operator",
        )
        plt.plot(
            indices,
            operators[
                ("archimedean_only", primary_cutoff)
            ].roots[: args.zero_count],
            linestyle="--",
            marker=".",
            label="Archimedean free operator",
        )
        plt.xlabel("Index")
        plt.ylabel("Energy")
        plt.title(
            "Pure Mangoldt prime-phase Dirac-Jacobi spectrum "
            f"(X={primary_cutoff}, d={args.operator_depth})"
        )
        plt.legend()
        plt.tight_layout()
        plt.savefig(paths["plot"], dpi=180)
        plt.close()

    summary = {
        "version": VERSION,
        "verdict": verdict,
        "scientific_scope": (
            "Finite zero-information-free prime-phase quantization and "
            "self-adjoint Dirac/Jacobi realization. Not an RH proof."
        ),
        "construction": {
            "free_count": "theta(E)/pi + 1",
            "prime_phase": (
                "-(1/pi) sum Lambda(n)/(log(n)*n^(1/2+1/log X)) "
                "sin(E log n)"
            ),
            "quantization": (
                "free_count(E)+prime_phase(E)=k-1/2 in fixed free cell"
            ),
            "primary_cutoff": primary_cutoff,
            "primary_cutoff_rule": (
                "ceil((depth-th free gamma level)^2)"
                if args.primary_cutoff <= 0
                else "explicit"
            ),
            "jacobi_weights": "sin^2(k*pi/(d+1)), normalized",
            "dirac": "[[0,J],[J,0]]",
            "scale_fit": None,
            "shift_fit": None,
            "index_offset_fit": None,
        },
        "freeze_manifest": freeze_manifest,
        "primary_structure": {
            "operator_hash": primary.operator_hash,
            "maximum_quantization_residual": float(
                np.max(primary.root_residuals)
            ),
            "minimum_roots_per_cell": int(
                np.min(primary.root_counts_per_cell)
            ),
            "maximum_roots_per_cell": int(
                np.max(primary.root_counts_per_cell)
            ),
            "jacobi_self_adjoint_defect": (
                primary.jacobi_self_adjoint_defect
            ),
            "jacobi_minimum_offdiagonal": (
                primary.jacobi_minimum_offdiagonal
            ),
            "jacobi_orthogonality_defect": (
                primary.jacobi_orthogonality_defect
            ),
            "jacobi_eigenvalue_preservation_error": (
                primary.jacobi_eigenvalue_preservation_error
            ),
            "dirac_chiral_defect": primary.dirac_chiral_defect,
            "maximum_relative_cutoff_drift": (
                maximum_relative_cutoff_drift
            ),
        },
        "signal_metrics": asdict(signal_metrics),
        "archimedean_metrics": asdict(archimedean_metrics),
        "best_disruption_control": best_disruption,
        "control_ranking": controls,
        "cutoff_stability": stability_rows,
        "gates": gates,
        "failed_gates": failed,
        "next": next_step,
        "runtime_seconds": time.perf_counter() - started,
        "outputs": {name: str(path) for name, path in paths.items()},
    }
    paths["summary"].write_text(
        json.dumps(json_safe(summary), indent=2),
        encoding="utf-8",
    )

    report = [
        "# Pure Mangoldt prime-phase Dirac-Jacobi operator audit",
        "",
        f"**{verdict}**",
        "",
        "## Firewall",
        "",
        "- Operator roots and matrices were frozen and hashed before zero access.",
        "- No Xi, xi, zeta, zeta derivative, zero list, or RvM zero ordinate entered construction.",
        "- No scale, shift, sign, cutoff, depth, or index offset was fitted from zero errors.",
        "",
        "## Primary result",
        "",
        f"- Cutoff/depth: `{primary_cutoff}/{args.operator_depth}`",
        f"- Operator hash: `{primary.operator_hash}`",
        f"- Signal RMSE: `{signal_metrics.rmse:.12f}`",
        f"- Signal relative RMSE: `{signal_metrics.relative_rmse:.6%}`",
        f"- Signal spacing relative RMSE: `{signal_metrics.spacing_relative_rmse:.6%}`",
        f"- Archimedean RMSE: `{archimedean_metrics.rmse:.12f}`",
        f"- Best disruption control: `{best_disruption['family']}` "
        f"with RMSE `{float(best_disruption['rmse']):.12f}`",
        f"- Maximum quantization residual: `{np.max(primary.root_residuals):.3e}`",
        f"- Minimum Jacobi off-diagonal: `{primary.jacobi_minimum_offdiagonal:.6e}`",
        f"- Relative cutoff drift: `{maximum_relative_cutoff_drift:.6e}`",
        "",
        "## Failed gates",
        "",
        *(
            ["- None"]
            if not failed
            else [f"- `{gate}`" for gate in failed]
        ),
        "",
        "## Next",
        "",
        next_step,
        "",
    ]
    paths["report"].write_text("\n".join(report), encoding="utf-8")

    print("=" * 164)
    print(f"VERDICT          : {verdict}")
    print(
        "PRIMARY OPERATOR : "
        f"X={primary_cutoff} d={args.operator_depth} "
        f"hash={primary.operator_hash[:20]}..."
    )
    print(
        "STRUCTURE        : "
        f"max_resid={np.max(primary.root_residuals):.3e} "
        f"min_b={primary.jacobi_minimum_offdiagonal:+.6e} "
        f"eig_preservation="
        f"{primary.jacobi_eigenvalue_preservation_error:.3e} "
        f"cutoff_drift={maximum_relative_cutoff_drift:.3e}"
    )
    print(
        "ZERO SPECTRUM    : "
        f"RMSE={signal_metrics.rmse:.12f} "
        f"relative={signal_metrics.relative_rmse:.6%} "
        f"spacing_relative="
        f"{signal_metrics.spacing_relative_rmse:.6%}"
    )
    print(
        "ARCHIMEDEAN      : "
        f"RMSE={archimedean_metrics.rmse:.12f} "
        f"gain={archimedean_metrics.rmse - signal_metrics.rmse:+.12f}"
    )
    print(
        "BEST CONTROL     : "
        f"{best_disruption['family']} "
        f"RMSE={float(best_disruption['rmse']):.12f} "
        f"gain={float(best_disruption['rmse']) - signal_metrics.rmse:+.12f}"
    )
    print(
        "FAILED GATES     : "
        + (", ".join(failed) if failed else "none")
    )
    print(f"NEXT             : {next_step}")
    print("OUTPUTS")
    for path in paths.values():
        print(f"  {path}")
    print("=" * 164)

    if args.strict and failed:
        return 2
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build and test a zero-information-free Mangoldt prime-phase "
            "Dirac/Jacobi operator."
        ),
    )
    parser.add_argument("--operator-depth", type=int, default=20)
    parser.add_argument("--zero-count", type=int, default=12)
    parser.add_argument(
        "--primary-cutoff",
        type=int,
        default=0,
        help=(
            "0 uses ceil((depth-th free gamma level)^2), with no zero data."
        ),
    )
    parser.add_argument("--grid-points", type=int, default=2048)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--zero-dps", type=int, default=60)
    parser.add_argument(
        "--families",
        default=(
            "signal,archimedean_only,prime_only,prime_power_only,"
            "signflip,permuted_weights,phase_scramble,smooth_integer"
        ),
    )

    parser.add_argument(
        "--maximum-quantization-residual",
        type=float,
        default=1e-8,
    )
    parser.add_argument(
        "--maximum-self-adjoint-defect",
        type=float,
        default=1e-12,
    )
    parser.add_argument(
        "--minimum-positive-offdiagonal",
        type=float,
        default=1e-12,
    )
    parser.add_argument(
        "--maximum-eigenvalue-preservation-error",
        type=float,
        default=1e-9,
    )
    parser.add_argument(
        "--maximum-chiral-defect",
        type=float,
        default=1e-12,
    )
    parser.add_argument(
        "--maximum-relative-cutoff-drift",
        type=float,
        default=0.005,
    )
    parser.add_argument(
        "--maximum-relative-rmse",
        type=float,
        default=0.005,
    )
    parser.add_argument(
        "--maximum-spacing-relative-rmse",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--minimum-control-margin",
        type=float,
        default=0.0,
    )

    parser.add_argument(
        "--out-prefix",
        default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator",
    )
    parser.add_argument("--strict", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.zero_count > args.operator_depth:
            raise ValueError("--zero-count cannot exceed --operator-depth.")
        return run(args)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
