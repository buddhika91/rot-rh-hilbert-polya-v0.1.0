from __future__ import annotations
import math
import numpy as np


def partial_moments(roots, order: int) -> np.ndarray:
    r = np.asarray(roots, dtype=float)
    if np.any(r <= 0):
        raise ValueError("roots must be positive")
    return np.asarray([np.sum(r ** (-2*m)) for m in range(1, order+1)], dtype=float)


def weyl_tail_moments(Ecut: float, order: int) -> np.ndarray:
    E = float(Ecut)
    if E <= 0:
        raise ValueError("Ecut must be positive")
    vals = []
    for m in range(1, int(order)+1):
        q = 2*m - 1
        vals.append(E**(1-2*m)/(2*math.pi*q) * (math.log(E/(2*math.pi)) + 1/q))
    return np.asarray(vals, dtype=float)


def midpoint_moments(Ecut: float, order: int) -> np.ndarray:
    E = float(Ecut)
    if E <= 0:
        raise ValueError("Ecut must be positive")
    return np.asarray([-0.5 * E**(-2*m) for m in range(1, int(order)+1)], dtype=float)


def finite_det(roots, tgrid) -> np.ndarray:
    r = np.asarray(roots, dtype=float)
    t = np.asarray(tgrid, dtype=float)
    if np.any(r <= 0):
        raise ValueError("roots must be positive")
    inv2 = 1.0/(r*r)
    return np.asarray([np.prod(1.0 - x*x*inv2) for x in t], dtype=float)


def completed_det(raw, tgrid, tail_moments) -> np.ndarray:
    raw = np.asarray(raw, dtype=float)
    tgrid = np.asarray(tgrid, dtype=float)
    tail = np.asarray(tail_moments, dtype=float)
    if raw.shape != tgrid.shape:
        raise ValueError("raw and tgrid must have the same shape")
    out = np.empty_like(raw)
    for i, t in enumerate(tgrid):
        logtail = -sum(float(T) * float(t)**(2*m) / m for m, T in enumerate(tail, start=1))
        if not math.isfinite(logtail) or abs(logtail) > 100:
            raise RuntimeError(f"unstable tail exponent at t={t}: {logtail}")
        out[i] = raw[i] * math.exp(logtail)
    return out
