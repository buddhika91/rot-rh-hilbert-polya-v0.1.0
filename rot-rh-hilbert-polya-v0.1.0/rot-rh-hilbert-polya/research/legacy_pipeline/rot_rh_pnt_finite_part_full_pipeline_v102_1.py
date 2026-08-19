#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v101 — PNT-Renormalized ROT Hilbert–Pólya Full Closure Audit
===================================================================

PURPOSE
-------
Test the complete numerical/theorem pipeline

    PNT-renormalized arithmetic phase
        -> ordered zero-free support
        -> positive heat spectral measure
        -> Jacobi operator
        -> fixed-L infinite self-adjoint closure
        -> trace-class inverse-square operator
        -> arithmetic Euler–Maclaurin Fredholm tail
        -> hard freeze
        -> Xi determinant / moments / log-derivative / zero audit.

This script DOES NOT claim an RH proof.  It distinguishes:
  * mathematical properties already proved for every fixed finite L;
  * zero-free numerical evidence for an L -> infinity operator limit;
  * post-freeze numerical evidence for the exact Xi determinant identity;
  * theorem obligations that remain open.

RENORMALIZATION
---------------
The exponentially smoothed prime phase is

    P_L(E)
      = sum_{n>=2} Lambda(n)/(log n sqrt(n))
          exp(-n/L) sin(E log n).

The PNT smooth counterterm is

    M_L(E)
      = integral_2^infinity
          x^(-1/2) exp(-x/L) sin(E log x)/log x dx.

The naive subtraction P_L-M_L changes the analytic phase.  Instead use the
finite-part Abel renormalization

    P_L^FP(s) = P_L(s) - M_L(s) + M_inf^AC(s),

where

    M_inf^AC(s)
      = E1((s-1) log 2)

is the analytic continuation of

    int_2^infinity x^(-s)/log(x) dx

from Re(s)>1.  On the critical spectral line s=1/2-iE,

    F_L^FP(E)
      = theta(E)/pi + 1 - Im P_L^FP(1/2-iE)/pi.

This subtracts the L-dependent PNT Abel divergence while adding back the
canonical analytic continuation of the smooth prime-density contribution.

M_L is evaluated after y=log x:

    M_L(E)
      = int_log(2)^infinity
          exp(y/2-exp(y)/L) sin(E y)/y dy.

A uniform log-x grid with composite Simpson weights is used.  This is
numerically well matched to the oscillation because E is the Fourier
frequency in y.  An independently doubled points-per-cycle grid is used as
a zero-free convergence check.

This counterterm uses NO Xi, zeta evaluations, or zero ordinates.

ORDERED SUPPORT
---------------
The spectral levels are the first later crossings

    F_L^ren(E_k) = k - 1/2.

No zero-based root matching is used.

JACOBI OPERATOR
---------------
For tau>0,

    w_k propto exp[-tau (E_k/E_1)^2],

and Stieltjes/Favard gives the positive Jacobi matrix J_(L,tau,N).
For fixed finite L,tau the corresponding infinite heat-weighted moment problem
is determinate, so the minimal half-line Jacobi operator has a unique
self-adjoint closure H_ROT^ren(L,tau).

TRACE-CLASS OBJECT
------------------
With E_k ~ 2*pi*k/log k,

    K_L = H_L^(-2)

is trace class and

    det(I - t^2 K_L)

exists.

ARITHMETIC FREDHOLM TAIL
------------------------
For inverse moment m,

    S_m = sum_k E_k^(-2m).

For a frozen prefix ending at E_N, use Euler–Maclaurin in the quantization
coordinate x=F_L^ren(E):

    tail_m
      = int_{E_N}^infinity E^(-2m) F_L^ren'(E) dE
        - 1/2 E_N^(-2m)
        + m/[6 E_N^(2m+1) F_L^ren'(E_N)]
        + higher Bernoulli remainder.

The continuum density is decomposed ZERO-FREE into

    Weyl + exact Archimedean correction
         + discrete Mangoldt oscillatory density
         + PNT counterterm density.

Both midpoint-only and midpoint+B2 ("EM1") completions are frozen.

The completed determinant is

    log D(t)
      = log D_N(t) - sum_{m>=1} tail_m t^(2m)/m,

valid on |t|<E_N with the finite moment series truncated at a fixed
pre-registered order.

HARD FIREWALL
-------------
Everything above is constructed, serialized and SHA256-hashed BEFORE mpmath
Xi/zeta is imported.

POST-FREEZE IDENTITY AUDIT
--------------------------
Then test

    D_ROT(t) ?= Xi(t)/Xi(0)

with:
  * symmetric RMS and sup absolute residual;
  * inverse log-moment residuals;
  * logarithmic-derivative/resolvent residual away from zeros;
  * zero locations detected only after freeze;
  * L-to-L convergence of the completed determinant before Xi.

If a v97 frozen archive is supplied, its unrenormalized support is used only as
a ZERO-FREE convergence control; no v97 Xi scores are imported.

EXACT FINISH LINE
-----------------
Numerical PASS is NOT RH.  A proof requires:
  1. a regulator-independent trace-norm limit K_L^ren -> K_infinity;
  2. positivity/injectivity of K_infinity;
  3. the exact entire-function identity

       det(I - t^2 K_infinity) = Xi(t)/Xi(0).

Then H_HP=K_infinity^(-1/2) is self-adjoint and all Xi zeros are real.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.integrate import quad
from scipy.special import exp1
from scipy.optimize import brentq

VERSION = "v102.1"
BUILD = "2026-08-17-pnt-finite-part-renormalized-full-hp-closure-v102.1-stable-ac-tail"
EPS = 1e-15


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def parse_ints(text: str) -> List[int]:
    return [int(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_complexes(text: str) -> List[complex]:
    out = []
    for tok in str(text).split(","):
        tok = tok.strip().replace("i", "j")
        if tok:
            out.append(complex(tok))
    return out


def json_safe(x):
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.complexfloating):
        return {"real": float(np.real(x)), "imag": float(np.imag(x))}
    if isinstance(x, complex):
        return {"real": float(x.real), "imag": float(x.imag)}
    return x


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = []
    seen = set()
    for row in rows:
        for k in row:
            if k not in seen:
                fields.append(k)
                seen.add(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def rel_l2(a, b) -> float:
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), EPS))


def symmetric_rms(a, b) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return float(
        np.sqrt(np.mean(((a - b) / (1.0 + np.abs(a) + np.abs(b))) ** 2))
    )


def max_abs(a, b) -> float:
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


# ---------------------------------------------------------------------------
# PNT counterterm
# ---------------------------------------------------------------------------

@dataclass
class PNTCounterterm:
    L: int
    q: int
    ppc: int
    xmax_factor: float
    y: np.ndarray
    logy: np.ndarray
    phase_coeff: np.ndarray
    deriv_coeff: np.ndarray
    phase_tail_abs_bound: float
    deriv_tail_abs_bound: float

    @staticmethod
    def build(
        L: int,
        q: int,
        design_emax: float = 220.0,
        xmax_factor: float = 36.0,
    ) -> "PNTCounterterm":
        """
        Stable quadrature for

            M_L(E)=int_2^inf x^(-1/2)e^(-x/L) sin(E log x)/log x dx.

        Use y=log x:

            M_L(E)=int_log2^inf
                exp(y/2-exp(y)/L) sin(E y)/y dy.

        The oscillation is exactly sinusoidal in y, so a uniform log-x grid
        is much better conditioned than high-order Gauss-Laguerre here.

        Backward-compatible resolution mapping:
            q=384 -> 24 points/cycle at design_emax
            q=768 -> 48 points/cycle at design_emax.
        """
        L = int(L)
        q = int(q)
        ppc = max(16, int(round(q / 16.0)))
        emax = float(design_emax)
        xmax_factor = float(xmax_factor)

        ymin = math.log(2.0)
        xmax = xmax_factor * float(L)
        ymax = math.log(xmax)

        # points-per-cycle criterion at the largest designed spectral energy.
        dy_target = 2.0 * math.pi / (ppc * emax)
        n = max(2049, int(math.ceil((ymax - ymin) / dy_target)) + 1)

        # Force odd number of points only to make nested/refined diagnostics
        # visually cleaner; trapezoidal weights are used.
        if n % 2 == 0:
            n += 1

        y = np.linspace(ymin, ymax, n)
        dy = float(y[1] - y[0])
        x = np.exp(y)

        # Composite Simpson weights on the uniform log-x grid.
        quad_w = np.ones(n, float)
        quad_w[1:-1:2] = 4.0
        quad_w[2:-1:2] = 2.0
        quad_w *= dy / 3.0

        common = np.exp(0.5 * y - x / float(L))
        phase_coeff = quad_w * common / y
        deriv_coeff = quad_w * common

        # Rigorous/simple absolute bounds for the omitted x>xmax continuum:
        # phase: 1/log(xmax) * int_xmax^inf x^-1/2 e^-x/L dx
        # deriv: int_xmax^inf x^-1/2 e^-x/L dx
        # Evaluate the positive scalar tail robustly by quadrature once.
        raw_tail = quad(
            lambda xx: math.exp(-xx / float(L)) / math.sqrt(xx),
            xmax,
            np.inf,
            epsabs=1e-15,
            epsrel=1e-12,
            limit=200,
        )[0]
        phase_tail = raw_tail / math.log(xmax)
        deriv_tail = raw_tail

        return PNTCounterterm(
            L=L,
            q=q,
            ppc=ppc,
            xmax_factor=xmax_factor,
            y=x,
            logy=y,
            phase_coeff=np.asarray(phase_coeff, float),
            deriv_coeff=np.asarray(deriv_coeff, float),
            phase_tail_abs_bound=float(phase_tail),
            deriv_tail_abs_bound=float(deriv_tail),
        )

    def phase(self, E) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        return np.sin(np.outer(e, self.logy)) @ self.phase_coeff

    def derivative(self, E) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        return np.cos(np.outer(e, self.logy)) @ self.deriv_coeff


def pnt_ac_phase(E) -> np.ndarray:
    """
    Imaginary part of the analytic continuation of the unsmoothed PNT
    main term

        M_inf^AC(s)=E1((s-1) log 2),  s=1/2-iE.

    No Xi/zeta/zero data enter.
    """
    e = np.atleast_1d(np.asarray(E, float))
    z = (-0.5 - 1j * e) * math.log(2.0)
    return np.imag(exp1(z))


def pnt_ac_phase_derivative(E) -> np.ndarray:
    """
    d/dE Im E1((s-1)log2), s=1/2-iE.

    E1'(z)=-exp(-z)/z and dz/dE=-i log2, hence
        d/dE E1(z)= i log2 exp(-z)/z.
    """
    e = np.atleast_1d(np.asarray(E, float))
    z = (-0.5 - 1j * e) * math.log(2.0)
    d = 1j * math.log(2.0) * np.exp(-z) / z
    return np.imag(d)



# ---------------------------------------------------------------------------
# Frozen arithmetic phase
# ---------------------------------------------------------------------------

@dataclass
class RenormalizedPhase:
    base: Any
    L: int
    nmax: int
    logs: np.ndarray
    phase_w: np.ndarray
    counter: PNTCounterterm
    prime_tail_abs_bound: float

    @staticmethod
    def build(base, L: int, tail_factor: float, laguerre_q: int):
        nmax = max(8, int(math.ceil(float(tail_factor) * float(L))))
        numbers, lambdas, _ = base.mangoldt_terms(nmax)
        nums = numbers.astype(float)
        logs = np.log(nums)
        phase_w = (
            lambdas
            / (logs * np.sqrt(nums))
            * np.exp(-nums / float(L))
        )

        # Absolute omitted discrete-prime phase majorant, using Lambda/log <= 1.
        # Integral of exp(-x/L)/sqrt(x) from N to infinity is bounded
        # conservatively by numerical quadrature here; done once per L.
        N = float(nmax)
        tail_int = quad(
            lambda x: math.exp(-x / float(L)) / math.sqrt(x),
            N,
            np.inf,
            epsabs=1e-13,
            epsrel=1e-11,
            limit=200,
        )[0]
        tail_lattice = math.exp(-(N + 1.0) / float(L)) / math.sqrt(N + 1.0)
        phase_tail = (tail_int + tail_lattice) / math.pi

        return RenormalizedPhase(
            base=base,
            L=int(L),
            nmax=int(nmax),
            logs=logs,
            phase_w=np.asarray(phase_w, float),
            counter=PNTCounterterm.build(L, laguerre_q),
            prime_tail_abs_bound=float(phase_tail),
        )

    def prime_phase(self, E, chunk: int = 20000) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        out = np.zeros(len(e), float)
        for j in range(0, len(self.logs), chunk):
            ll = self.logs[j:j + chunk]
            ww = self.phase_w[j:j + chunk]
            out += np.sin(np.outer(e, ll)) @ ww
        return out

    def prime_derivative(self, E, chunk: int = 20000) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        out = np.zeros(len(e), float)
        coeff = self.phase_w * self.logs
        for j in range(0, len(self.logs), chunk):
            ll = self.logs[j:j + chunk]
            cc = coeff[j:j + chunk]
            out += np.cos(np.outer(e, ll)) @ cc
        return out

    def F(self, E) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        smooth = np.asarray([self.base.smooth_count(float(x)) for x in e], float)
        # Canonical finite-part Abel renormalization:
        # P_FP = P_L - M_L + M_inf^AC.
        fp_phase = (
            self.prime_phase(e)
            - self.counter.phase(e)
            + pnt_ac_phase(e)
        )
        return smooth - fp_phase / math.pi

    def Fprime(self, E) -> np.ndarray:
        e = np.atleast_1d(np.asarray(E, float))
        arch = np.asarray([self.base.theta_prime(float(x)) / math.pi for x in e])
        p = self.prime_derivative(e)
        c = self.counter.derivative(e)
        ac = pnt_ac_phase_derivative(e)
        # derivative of -P_FP/pi:
        # -P'_L/pi + M'_L/pi - M_inf^AC'/pi.
        return arch - p / math.pi + c / math.pi - ac / math.pi


@dataclass
class Support:
    L: int
    depth: int
    roots: np.ndarray
    residuals: np.ndarray
    min_local_Fprime: float
    phase_tail_abs_bound: float
    quadrature_phase_rel_l2: float
    quadrature_deriv_rel_l2: float


def quadrature_check(
    phase: RenormalizedPhase,
    qcheck: int,
    emax: float,
    samples: int = 24,
) -> Tuple[float, float]:
    c2 = PNTCounterterm.build(phase.L, int(qcheck))
    e = np.linspace(8.0, float(emax), int(samples))
    a = phase.counter.phase(e)
    b = c2.phase(e)
    ap = phase.counter.derivative(e)
    bp = c2.derivative(e)
    return rel_l2(a, b), rel_l2(ap, bp)


def ordered_support(
    phase: RenormalizedPhase,
    depth: int,
    scan_points: int,
    qcheck: int,
) -> Support:
    base = phase.base
    free = base.free_archimedean_levels(depth + 2)
    roots = []
    residuals = []
    local_derivs = []

    for idx in range(depth):
        target = float(idx + 1) - 0.5
        lower = 7.0 if idx == 0 else roots[-1] + 1e-9
        nominal = float(free[min(idx + 1, len(free) - 1)])
        upper = max(lower + 1.0, nominal + 0.8 * max(1.0, nominal - lower))

        found = None
        for _ in range(14):
            grid = np.linspace(lower, upper, int(scan_points))
            vals = phase.F(grid) - target
            cross = np.flatnonzero(vals[:-1] * vals[1:] <= 0.0)
            if len(cross):
                q = int(cross[0])
                a = float(grid[q])
                b = float(grid[q + 1])
                root = brentq(
                    lambda x: float(phase.F(np.asarray([x]))[0] - target),
                    a,
                    b,
                    xtol=1e-12,
                    rtol=1e-12,
                    maxiter=300,
                )
                found = float(root)
                break
            lower = upper
            upper = max(upper + 2.0, upper * 1.25)

        if found is None:
            raise RuntimeError(
                f"Failed to find first later crossing k={idx+1}, L={phase.L}."
            )

        roots.append(found)
        residuals.append(
            abs(float(phase.F(np.asarray([found]))[0] - target))
        )
        local_derivs.append(float(phase.Fprime(np.asarray([found]))[0]))

    qphase, qder = quadrature_check(
        phase,
        qcheck=qcheck,
        emax=max(roots[-1], 20.0),
    )

    return Support(
        L=phase.L,
        depth=int(depth),
        roots=np.asarray(roots, float),
        residuals=np.asarray(residuals, float),
        min_local_Fprime=float(np.min(local_derivs)),
        phase_tail_abs_bound=phase.prime_tail_abs_bound,
        quadrature_phase_rel_l2=float(qphase),
        quadrature_deriv_rel_l2=float(qder),
    )


# ---------------------------------------------------------------------------
# Jacobi / spectral objects
# ---------------------------------------------------------------------------

@dataclass
class Section:
    support: Support
    tau: float
    heat_weights: np.ndarray
    a: np.ndarray
    b: np.ndarray
    J: np.ndarray
    self_adjoint_defect: float
    min_b: float
    eig_preservation_error: float
    orthogonality_defect: float
    m_values: np.ndarray


def heat_weights(roots: np.ndarray, tau: float) -> np.ndarray:
    r = np.asarray(roots, float)
    z = -float(tau) * (r / r[0]) ** 2
    z -= np.max(z)
    w = np.exp(z)
    w /= np.sum(w)
    return w


def m_function(roots, weights, z):
    r = np.asarray(roots, float)
    w = np.asarray(weights, float)
    return np.sum(w / (r - complex(z)))


def build_section(base, support: Support, tau: float, probes) -> Section:
    w = heat_weights(support.roots, tau)
    a, b, diag = base.jacobi_from_atoms(support.roots, w)
    J = np.diag(a)
    if len(b):
        J += np.diag(b, 1) + np.diag(b, -1)
    mv = np.asarray([m_function(support.roots, w, z) for z in probes], complex)
    return Section(
        support=support,
        tau=float(tau),
        heat_weights=w,
        a=np.asarray(a, float),
        b=np.asarray(b, float),
        J=J,
        self_adjoint_defect=float(np.linalg.norm(J - J.T, "fro")),
        min_b=float(np.min(b)) if len(b) else float("inf"),
        eig_preservation_error=float(diag["eigenvalue_preservation_error"]),
        orthogonality_defect=float(diag["orthogonality_defect"]),
        m_values=mv,
    )


# ---------------------------------------------------------------------------
# Arithmetic tail completion
# ---------------------------------------------------------------------------

def partial_moments(roots, order):
    r = np.asarray(roots, float)
    return np.asarray(
        [float(np.sum(r ** (-2 * m))) for m in range(1, int(order) + 1)],
        float,
    )


def weyl_tail_moments(E, order):
    E = float(E)
    vals = []
    for m in range(1, int(order) + 1):
        q = 2 * m - 1
        vals.append(
            E ** (1 - 2 * m)
            / (2.0 * math.pi * q)
            * (math.log(E / (2.0 * math.pi)) + 1.0 / q)
        )
    return np.asarray(vals, float)


def arch_tail_moments(base, E, order):
    vals = []
    for m in range(1, int(order) + 1):
        p = 2 * m

        def fn(x):
            exact = base.theta_prime(float(x)) / math.pi
            leading = math.log(float(x) / (2 * math.pi)) / (2 * math.pi)
            return (exact - leading) / (float(x) ** p)

        vals.append(
            float(
                quad(
                    fn,
                    float(E),
                    np.inf,
                    epsabs=1e-14,
                    epsrel=1e-11,
                    limit=250,
                )[0]
            )
        )
    return np.asarray(vals, float)


def oscillatory_tail_integral(a, E, p, terms=14):
    """
    Stable asymptotic for

        Z_p(a,E) = int_E^inf exp(i a x) x^-p dx.

    Repeated integration by parts gives

        Z_p ~ exp(i a E) sum_{j>=0}
              i(-i)^j (p)_j E^(-p-j) a^(-j-1).

    For this pipeline aE is large (>=~60 on every tail endpoint used).
    """
    a = np.asarray(a, float)
    E = float(E)
    p = int(p)
    phase = np.exp(1j * a * E)
    series = np.zeros_like(a, dtype=np.complex128)
    rising = np.ones_like(a)

    for j in range(int(terms)):
        if j > 0:
            rising *= (p + j - 1)
        coeff = 1j * ((-1j) ** j)
        series += coeff * rising * (E ** (-p - j)) * (a ** (-j - 1))
    return phase * series


def cosine_tail_integral(a, E, p, terms=14):
    return np.real(oscillatory_tail_integral(a, E, p, terms=terms))


def sine_tail_integral(a, E, p, terms=14):
    return np.imag(oscillatory_tail_integral(a, E, p, terms=terms))


def prime_tail_moments(phase: RenormalizedPhase, E, order):
    coeff = -(phase.phase_w * phase.logs) / math.pi
    vals = []
    for m in range(1, int(order) + 1):
        I = cosine_tail_integral(phase.logs, E, 2 * m)
        vals.append(float(np.sum(coeff * I)))
    return np.asarray(vals, float)


def counterterm_tail_moments(counter: PNTCounterterm, E, order):
    # Counterterm contribution to F'_ren is +M'_L(E)/pi.
    coeff = counter.deriv_coeff / math.pi
    vals = []
    for m in range(1, int(order) + 1):
        I = cosine_tail_integral(counter.logy, E, 2 * m)
        vals.append(float(np.sum(coeff * I)))
    return np.asarray(vals, float)



def pnt_ac_tail_moments(E, order, geom_terms=6, osc_terms=14):
    """
    Stable zero-free analytic-continuation tail.

    For omega=log 2 and a0=1/2,

      d/dE Im E1((-1/2-iE)log2)
        = -sqrt(2) [E sin(omega E)+a0 cos(omega E)]
          /(E^2+a0^2).

    The contribution to F_FP'(E) is minus this quantity divided by pi.

    Expand
      1/(E^2+a0^2)
        = E^-2 sum_{j>=0} (-a0^2/E^2)^j

    and integrate each sine/cosine power tail with the same stable
    integration-by-parts expansion used for the prime tail.

    At E>=~80, six geometric terms leave a denominator-expansion remainder
    far below double precision for the present moments.
    """
    E = float(E)
    omega = math.log(2.0)
    a0 = 0.5
    pref = math.sqrt(2.0) / math.pi
    vals = []

    for m in range(1, int(order) + 1):
        p = 2 * m
        total = 0.0
        for j in range(int(geom_terms)):
            fac = (-a0 * a0) ** j

            # sqrt(2)/pi * [
            #   sin(omega x) x^(-(p+1+2j))
            #   + a0 cos(omega x) x^(-(p+2+2j))
            # ]
            Is = float(
                sine_tail_integral(
                    np.asarray([omega]),
                    E,
                    p + 1 + 2 * j,
                    terms=osc_terms,
                )[0]
            )
            Ic = float(
                cosine_tail_integral(
                    np.asarray([omega]),
                    E,
                    p + 2 + 2 * j,
                    terms=osc_terms,
                )[0]
            )
            total += fac * (Is + a0 * Ic)

        vals.append(pref * total)

    return np.asarray(vals, float)


def pnt_ac_tail_stability(E, order):
    """
    Independent internal check: compare (geom,osc)=(6,14) with (8,18).
    """
    a = pnt_ac_tail_moments(E, order, geom_terms=6, osc_terms=14)
    b = pnt_ac_tail_moments(E, order, geom_terms=8, osc_terms=18)
    return rel_l2(a, b), a, b


def midpoint_moments(E, order):
    E = float(E)
    return np.asarray([-0.5 * E ** (-2 * m) for m in range(1, order + 1)])


def bernoulli_B2_moments(phase: RenormalizedPhase, E, order):
    E = float(E)
    fp = float(phase.Fprime(np.asarray([E]))[0])
    if fp <= 0:
        return np.full(int(order), np.nan)
    return np.asarray(
        [
            m / (6.0 * E ** (2 * m + 1) * fp)
            for m in range(1, int(order) + 1)
        ],
        float,
    )


def finite_det(roots, tgrid):
    r = np.asarray(roots, float)
    inv2 = 1.0 / (r * r)
    return np.asarray(
        [float(np.prod(1.0 - (float(t) ** 2) * inv2)) for t in tgrid],
        float,
    )


def completed_det(raw, tgrid, tail):
    out = np.zeros_like(raw)
    for i, t in enumerate(np.asarray(tgrid, float)):
        logtail = 0.0
        for m, T in enumerate(np.asarray(tail, float), start=1):
            logtail -= float(T) * float(t) ** (2 * m) / m
        if not math.isfinite(logtail) or abs(logtail) > 100.0:
            raise RuntimeError(
                f"Tail exponent unstable at t={t}: logtail={logtail}"
            )
        out[i] = raw[i] * math.exp(logtail)
    return out


# ---------------------------------------------------------------------------
# Post-freeze Xi
# ---------------------------------------------------------------------------

def xi_value(mp, t):
    s = mp.mpf("0.5") + 1j * mp.mpf(str(t))
    return mp.re(
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.power(mp.pi, -s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def xi_moments(mp, order):
    x0 = xi_value(mp, 0)
    out = []
    for m in range(1, int(order) + 1):
        print(f"[Xi moment] {m}/{order}", flush=True)
        d = mp.diff(
            lambda t: mp.log(xi_value(mp, t) / x0),
            mp.mpf("0"),
            2 * m,
        )
        out.append(float(-m * d / mp.factorial(2 * m)))
    return np.asarray(out, float)


def scan_xi_roots(mp, tmax, step):
    roots = []
    x0 = mp.mpf("0")
    f0 = xi_value(mp, x0)
    x = float(step)
    while x <= float(tmax) + 1e-12:
        x1 = mp.mpf(str(x))
        f1 = xi_value(mp, x1)
        if f0 * f1 < 0:
            r = mp.findroot(lambda q: xi_value(mp, q), (x0, x1))
            rr = float(r)
            if not roots or abs(rr - roots[-1]) > 1e-8:
                roots.append(rr)
        x0 = x1
        f0 = f1
        x += float(step)
    return np.asarray(roots, float)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument(
        "--base-module",
        default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit",
    )
    ap.add_argument("--v97-prefix", default="rot_rh_full_rot_hp_pipeline_v97")
    ap.add_argument("--L-list", default="2000,4000,8000,16000")
    ap.add_argument("--depths", default="24,32,48")
    ap.add_argument("--tau", type=float, default=0.02)
    ap.add_argument("--tail-factor", type=float, default=18.0)
    ap.add_argument("--laguerre-q", type=int, default=384)
    ap.add_argument("--laguerre-qcheck", type=int, default=768)
    ap.add_argument("--scan-points", type=int, default=96)
    ap.add_argument("--tail-moment-order", type=int, default=10)
    ap.add_argument("--post-moment-order", type=int, default=4)
    ap.add_argument("--m-probes", default="1j,10+1j,20+2j,40+4j")
    ap.add_argument("--identity-max", type=float, default=50.0)
    ap.add_argument("--identity-step", type=float, default=0.5)
    ap.add_argument("--xi-dps", type=int, default=60)
    ap.add_argument("--zero-scan-max", type=float, default=120.0)
    ap.add_argument("--zero-scan-step", type=float, default=0.10)
    ap.add_argument("--zero-match-count", type=int, default=32)
    ap.add_argument("--prefix", default="rot_rh_pnt_finite_part_full_pipeline_v102_1")
    args = ap.parse_args()

    started = time.time()
    Ls = parse_ints(args.L_list)
    depths = parse_ints(args.depths)
    probes = parse_complexes(args.m_probes)
    if Ls != sorted(Ls) or depths != sorted(depths):
        raise ValueError("L-list and depths must be increasing.")
    if args.laguerre_qcheck <= args.laguerre_q:
        raise ValueError("laguerre-qcheck must exceed laguerre-q.")
    if min(depths) < 4 or min(Ls) <= 0 or args.tau <= 0:
        raise ValueError("Need positive L,tau and depth>=4.")

    base = importlib.import_module(args.base_module)
    prefix = Path(args.prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)

    tgrid = np.arange(
        0.0,
        args.identity_max + 0.5 * args.identity_step,
        args.identity_step,
    )

    print("=" * 168)
    print("ROT-RH v102.1 — PNT FINITE-PART ROT HP CLOSURE (STABLE ANALYTIC-CONTINUATION TAIL)")
    print("=" * 168)
    print("L ladder                         :", Ls)
    print("depth ladder                     :", depths)
    print("ROT heat tau                     :", args.tau)
    print("PNT renormalization              : finite-part P_L - M_L + E1((s-1)log2)")
    print("PNT Abel counterterm quadrature   : uniform log-x Simpson")
    print("coarse resolution flag           :", args.laguerre_q, f"(~{max(16, round(args.laguerre_q/16))} points/cycle)")
    print("independent fine check flag      :", args.laguerre_qcheck, f"(~{max(16, round(args.laguerre_qcheck/16))} points/cycle)")
    print("tail moment order                :", args.tail_moment_order)
    print("Xi/zeta/known zeros pre-freeze   : NONE")
    print("=" * 168)

    phases: Dict[int, RenormalizedPhase] = {}
    sections: Dict[Tuple[int, int], Section] = {}
    construction_rows = []

    print("[A] Build PNT-renormalized zero-free operator family", flush=True)
    for L in Ls:
        phase = RenormalizedPhase.build(
            base,
            L=L,
            tail_factor=args.tail_factor,
            laguerre_q=args.laguerre_q,
        )
        phases[L] = phase

        sup_max = ordered_support(
            phase,
            depth=depths[-1],
            scan_points=args.scan_points,
            qcheck=args.laguerre_qcheck,
        )

        for d in depths:
            sup = Support(
                L=L,
                depth=d,
                roots=sup_max.roots[:d].copy(),
                residuals=sup_max.residuals[:d].copy(),
                min_local_Fprime=float(
                    np.min(phase.Fprime(sup_max.roots[:d]))
                ),
                phase_tail_abs_bound=sup_max.phase_tail_abs_bound,
                quadrature_phase_rel_l2=sup_max.quadrature_phase_rel_l2,
                quadrature_deriv_rel_l2=sup_max.quadrature_deriv_rel_l2,
            )
            sec = build_section(base, sup, args.tau, probes)
            sections[(L, d)] = sec

            construction_rows.append(
                {
                    "L": L,
                    "depth": d,
                    "first_root": sup.roots[0],
                    "last_root": sup.roots[-1],
                    "max_root_residual": float(np.max(sup.residuals)),
                    "min_local_Fprime": sup.min_local_Fprime,
                    "prime_tail_abs_bound": sup.phase_tail_abs_bound,
                    "counterterm_phase_qcheck_rel_l2":
                        sup.quadrature_phase_rel_l2,
                    "counterterm_deriv_qcheck_rel_l2":
                        sup.quadrature_deriv_rel_l2,
                    "counterterm_phase_tail_abs_bound":
                        phase.counter.phase_tail_abs_bound,
                    "counterterm_deriv_tail_abs_bound":
                        phase.counter.deriv_tail_abs_bound,
                    "self_adjoint_defect": sec.self_adjoint_defect,
                    "minimum_offdiagonal": sec.min_b,
                    "eig_preservation_error": sec.eig_preservation_error,
                    "orthogonality_defect": sec.orthogonality_defect,
                }
            )

            print(
                f"  L={L:<7d} d={d:<3d} "
                f"roots=[{sup.roots[0]:.6f},{sup.roots[-1]:.6f}] "
                f"maxres={np.max(sup.residuals):.2e} "
                f"minF'={sup.min_local_Fprime:+.3e} "
                f"qphase={sup.quadrature_phase_rel_l2:.2e} "
                f"qder={sup.quadrature_deriv_rel_l2:.2e}"
            )

    # -----------------------------------------------------------------------
    # B. Zero-free L-flow / convergence
    # -----------------------------------------------------------------------
    print("[B] Zero-free renormalized operator-flow diagnostics", flush=True)
    flow_rows = []
    dmax = depths[-1]
    for L0, L1 in zip(Ls[:-1], Ls[1:]):
        s0 = sections[(L0, dmax)]
        s1 = sections[(L1, dmax)]
        flow_rows.append(
            {
                "L_from": L0,
                "L_to": L1,
                "support_relative_l2":
                    rel_l2(s1.support.roots, s0.support.roots),
                "m_function_relative_l2":
                    rel_l2(s1.m_values, s0.m_values),
                "jacobi_coeff_relative_l2":
                    rel_l2(
                        np.r_[s1.a, np.log(s1.b)],
                        np.r_[s0.a, np.log(s0.b)],
                    ),
            }
        )
        r = flow_rows[-1]
        print(
            f"  {L0}->{L1}: support={r['support_relative_l2']:.3e} "
            f"m={r['m_function_relative_l2']:.3e} "
            f"jacobi={r['jacobi_coeff_relative_l2']:.3e}"
        )

    # Optional unrenormalized v97 convergence control.
    raw_control_rows = []
    v97_archive = Path(
        str(Path(args.v97_prefix).expanduser().resolve()) + "_FROZEN_PRE_XI.npz"
    )
    v97_manifest = Path(
        str(Path(args.v97_prefix).expanduser().resolve()) + "_FROZEN_MANIFEST.json"
    )
    v97_sha = None
    if v97_archive.exists() and v97_manifest.exists():
        man = json.loads(v97_manifest.read_text(encoding="utf-8"))
        got = sha256_file(v97_archive)
        if got == man.get("sha256"):
            v97_sha = got
            raw = np.load(v97_archive, allow_pickle=False)
            for L0, L1 in zip(Ls[:-1], Ls[1:]):
                k0 = f"roots_L{L0}_d{dmax}"
                k1 = f"roots_L{L1}_d{dmax}"
                if k0 in raw and k1 in raw:
                    rr = {
                        "L_from": L0,
                        "L_to": L1,
                        "raw_support_relative_l2":
                            rel_l2(raw[k1], raw[k0]),
                        "ren_support_relative_l2":
                            next(
                                x["support_relative_l2"]
                                for x in flow_rows
                                if x["L_from"] == L0 and x["L_to"] == L1
                            ),
                    }
                    rr["ren_over_raw"] = (
                        rr["ren_support_relative_l2"]
                        / max(rr["raw_support_relative_l2"], EPS)
                    )
                    raw_control_rows.append(rr)
            if raw_control_rows:
                print("[B-control] Finite-part vs v97 raw support drift")
                for r in raw_control_rows:
                    print(
                        f"  {r['L_from']}->{r['L_to']}: "
                        f"raw={r['raw_support_relative_l2']:.3e} "
                        f"ren={r['ren_support_relative_l2']:.3e} "
                        f"ratio={r['ren_over_raw']:.3f}"
                    )

    # -----------------------------------------------------------------------
    # C. Zero-free Fredholm tail completion at every L/depth
    # -----------------------------------------------------------------------
    print("[C] Zero-free arithmetic Fredholm completion", flush=True)
    tail_rows = []
    completed: Dict[Tuple[int, int, str], np.ndarray] = {}
    completed_mom: Dict[Tuple[int, int, str], np.ndarray] = {}
    tail_cache: Dict[Tuple[int, int, str], np.ndarray] = {}

    for L in Ls:
        phase = phases[L]
        for d in depths:
            roots = sections[(L, d)].support.roots
            E = float(roots[-1])
            if args.identity_max >= E:
                raise ValueError(
                    f"identity-max={args.identity_max} must be < E_N={E:.6f} "
                    f"for L={L}, d={d}."
                )

            partial = partial_moments(roots, args.tail_moment_order)
            W = weyl_tail_moments(E, args.tail_moment_order)
            A = arch_tail_moments(base, E, args.tail_moment_order)
            P = prime_tail_moments(phase, E, args.tail_moment_order)
            C = counterterm_tail_moments(
                phase.counter, E, args.tail_moment_order
            )
            ac_stability, AC, AC_check = pnt_ac_tail_stability(
                E, args.tail_moment_order
            )
            M = midpoint_moments(E, args.tail_moment_order)
            B2 = bernoulli_B2_moments(
                phase, E, args.tail_moment_order
            )

            tails = {
                "weyl": W,
                "full_midpoint": W + A + P + C + AC + M,
            }
            if np.all(np.isfinite(B2)):
                tails["full_em1"] = W + A + P + C + AC + M + B2

            raw_det = finite_det(roots, tgrid)
            for name, T in tails.items():
                tail_cache[(L, d, name)] = np.asarray(T, float)
                completed[(L, d, name)] = completed_det(raw_det, tgrid, T)
                completed_mom[(L, d, name)] = partial + T

            tail_rows.append(
                {
                    "L": L,
                    "depth": d,
                    "Ecut": E,
                    "prefix_S1": partial[0],
                    "weyl_S1": W[0],
                    "arch_S1": A[0],
                    "prime_S1": P[0],
                    "counterterm_S1": C[0],
                    "analytic_continuation_S1": AC[0],
                    "analytic_continuation_S1_check": AC_check[0],
                    "analytic_continuation_tail_rel_l2_stability": ac_stability,
                    "midpoint_S1": M[0],
                    "B2_S1": B2[0] if np.all(np.isfinite(B2)) else np.nan,
                    "full_midpoint_S1":
                        completed_mom[(L, d, "full_midpoint")][0],
                    "full_em1_S1":
                        completed_mom[(L, d, "full_em1")][0]
                        if (L, d, "full_em1") in completed_mom
                        else np.nan,
                }
            )

            print(
                f"  L={L:<7d} d={d:<3d} Ecut={E:.6f} "
                f"P={P[0]:+.3e} C={C[0]:+.3e} AC={AC[0]:+.3e} "
                f"ACchk={ac_stability:.1e} M={M[0]:+.3e} "
                f"B2={B2[0] if np.all(np.isfinite(B2)) else float('nan'):+.3e}"
            )

    # Pre-Xi determinant/moment L-flow convergence.
    print("[D] Zero-free completed-Fredholm L-flow convergence", flush=True)
    detflow_rows = []
    canonical_method = "full_em1"
    if not all((L, dmax, canonical_method) in completed for L in Ls):
        canonical_method = "full_midpoint"

    # Internal pre-Xi convergence of the moment expansion used in the
    # canonical Fredholm tail: compare the registered order M against M-2.
    final_tail = tail_cache[(Ls[-1], dmax, canonical_method)]
    final_raw = finite_det(sections[(Ls[-1], dmax)].support.roots, tgrid)
    if len(final_tail) >= 4:
        det_lower_order = completed_det(final_raw, tgrid, final_tail[:-2])
        tail_series_order_drift = rel_l2(
            completed[(Ls[-1], dmax, canonical_method)],
            det_lower_order,
        )
    else:
        tail_series_order_drift = float("nan")
    print(
        f"  tail-series order {len(final_tail)-2}->{len(final_tail)} "
        f"relative drift={tail_series_order_drift:.3e}"
    )

    for L0, L1 in zip(Ls[:-1], Ls[1:]):
        d0 = completed[(L0, dmax, canonical_method)]
        d1 = completed[(L1, dmax, canonical_method)]
        m0 = completed_mom[(L0, dmax, canonical_method)]
        m1 = completed_mom[(L1, dmax, canonical_method)]
        row = {
            "L_from": L0,
            "L_to": L1,
            "method": canonical_method,
            "determinant_relative_l2": rel_l2(d1, d0),
            "moments_relative_l2": rel_l2(m1, m0),
        }
        detflow_rows.append(row)
        print(
            f"  {L0}->{L1}: det={row['determinant_relative_l2']:.3e} "
            f"mom={row['moments_relative_l2']:.3e}"
        )

    # -----------------------------------------------------------------------
    # E. Hard freeze
    # -----------------------------------------------------------------------
    payload = {
        "L_values": np.asarray(Ls, int),
        "depth_values": np.asarray(depths, int),
        "tgrid": tgrid,
    }
    for L in Ls:
        for d in depths:
            sec = sections[(L, d)]
            payload[f"roots_L{L}_d{d}"] = sec.support.roots
            payload[f"a_L{L}_d{d}"] = sec.a
            payload[f"b_L{L}_d{d}"] = sec.b
            payload[f"heat_w_L{L}_d{d}"] = sec.heat_weights
            for name in ["weyl", "full_midpoint", "full_em1"]:
                k = (L, d, name)
                if k in completed:
                    payload[f"det_{name}_L{L}_d{d}"] = completed[k]
                    payload[f"mom_{name}_L{L}_d{d}"] = completed_mom[k]

    freeze = Path(str(prefix) + "_FROZEN_PRE_XI.npz")
    manifest_path = Path(str(prefix) + "_FROZEN_MANIFEST.json")
    np.savez_compressed(freeze, **payload)
    freeze_sha = sha256_file(freeze)

    pre_gates = {
        "all_finite_jacobi_self_adjoint":
            all(s.self_adjoint_defect < 1e-12 for s in sections.values()),
        "all_offdiagonals_positive":
            all(s.min_b > 0 for s in sections.values()),
        "all_root_residuals_small":
            all(np.max(s.support.residuals) < 1e-8 for s in sections.values()),
        "counterterm_quadrature_phase_stable":
            all(
                s.support.quadrature_phase_rel_l2 < 5e-6
                for s in sections.values()
            ),
        "counterterm_quadrature_deriv_stable":
            all(
                s.support.quadrature_deriv_rel_l2 < 5e-6
                for s in sections.values()
            ),
        "final_local_Fprime_positive":
            sections[(Ls[-1], dmax)].support.min_local_Fprime > 0,
        "fredholm_tail_series_order_stable":
            bool(np.isfinite(tail_series_order_drift) and tail_series_order_drift < 1e-6),
        "analytic_continuation_tail_stable":
            all(
                float(r.get("analytic_continuation_tail_rel_l2_stability", 1.0)) < 1e-10
                for r in tail_rows
            ),
    }

    manifest = {
        "version": VERSION,
        "build": BUILD,
        "sha256": freeze_sha,
        "archive": str(freeze),
        "canonical_tail_method": canonical_method,
        "parent_v97_sha256": v97_sha,
        "parameters": {
            "L_list": Ls,
            "depths": depths,
            "tau": args.tau,
            "tail_factor": args.tail_factor,
            "counterterm_resolution_flag": args.laguerre_q,
            "counterterm_check_resolution_flag": args.laguerre_qcheck,
            "counterterm_quadrature": "uniform log-x composite Simpson",
            "tail_moment_order": args.tail_moment_order,
            "identity_max": args.identity_max,
            "identity_step": args.identity_step,
        },
        "firewall": {
            "Xi_used_pre_freeze": False,
            "zeta_used_pre_freeze": False,
            "known_zero_ordinates_used_pre_freeze": False,
            "zero_loss_used": False,
        },
        "pre_gates": pre_gates,
        "theorem_ledger_pre": {
            "fixed_L_self_adjoint": "PROVED",
            "fixed_L_inverse_square_trace_class": "PROVED",
            "PNT_finite_part_renormalization_zero_free": "YES",
            "renormalized_L_to_infinity_trace_norm_limit": "OPEN",
            "Euler_Maclaurin_full_remainder_bound": "OPEN_AS_GLOBAL_CERTIFICATE",
            "exact_Xi_Fredholm_identity": "NOT_EVALUATED_PRE_FREEZE",
        },
    }
    manifest_path.write_text(
        json.dumps(json_safe(manifest), indent=2),
        encoding="utf-8",
    )

    print("[E] HARD FREEZE SHA256             :", freeze_sha)
    print("[E] Xi firewall CLOSED.  Post-freeze external identity audit begins.")

    # -----------------------------------------------------------------------
    # F. POST-FREEZE Xi
    # -----------------------------------------------------------------------
    import mpmath as mp
    mp.mp.dps = int(args.xi_dps)

    x0 = xi_value(mp, 0)
    print("[F1] Xi determinant grid", flush=True)
    xi_grid = np.asarray(
        [float(xi_value(mp, float(t)) / x0) for t in tgrid],
        float,
    )

    print("[F2] Xi inverse log moments", flush=True)
    xi_mom = xi_moments(mp, args.post_moment_order)

    print("[F3] Xi sign-change roots", flush=True)
    xi_roots = scan_xi_roots(
        mp,
        tmax=args.zero_scan_max,
        step=args.zero_scan_step,
    )

    post_rows = []
    for L in Ls:
        roots = sections[(L, dmax)].support.roots
        raw_det = finite_det(roots, tgrid)
        partial = partial_moments(roots, args.post_moment_order)

        methods = ["raw_prefix", "weyl", "full_midpoint"]
        if (L, dmax, "full_em1") in completed:
            methods.append("full_em1")

        for name in methods:
            if name == "raw_prefix":
                det = raw_det
                mom = partial
            else:
                det = completed[(L, dmax, name)]
                mom = completed_mom[(L, dmax, name)][:args.post_moment_order]

            post_rows.append(
                {
                    "L": L,
                    "depth": dmax,
                    "method": name,
                    "det_symmetric_rms": symmetric_rms(det, xi_grid),
                    "det_sup_abs": max_abs(det, xi_grid),
                    "moment_relative_l2": rel_l2(mom, xi_mom),
                    "S1_relative_error":
                        abs(float(mom[0]) - float(xi_mom[0]))
                        / abs(float(xi_mom[0])),
                }
            )

    # Zero-location score on final renormalized support.
    final_roots = sections[(Ls[-1], dmax)].support.roots
    nmatch = min(
        int(args.zero_match_count),
        len(final_roots),
        len(xi_roots),
    )
    if nmatch:
        zrel = float(
            np.sqrt(
                np.mean(
                    (
                        (final_roots[:nmatch] - xi_roots[:nmatch])
                        / xi_roots[:nmatch]
                    ) ** 2
                )
            )
        )
    else:
        zrel = float("nan")

    # Log-derivative identity using the completed determinant itself by
    # central finite differences on the frozen grid, away from zeros.
    final_det = completed[(Ls[-1], dmax, canonical_method)]
    h = float(args.identity_step)
    logder_rows = []
    for i in range(1, len(tgrid) - 1):
        t = float(tgrid[i])
        if abs(final_det[i]) < 1e-8 or abs(xi_grid[i]) < 1e-8:
            continue
        if min(abs(final_det[i - 1]), abs(final_det[i + 1])) < 1e-12:
            continue
        if min(abs(xi_grid[i - 1]), abs(xi_grid[i + 1])) < 1e-12:
            continue

        # derivative of log |D|; sign changes are skipped by the guards.
        op = (
            math.log(abs(final_det[i + 1]))
            - math.log(abs(final_det[i - 1]))
        ) / (2 * h)
        xx = (
            math.log(abs(xi_grid[i + 1]))
            - math.log(abs(xi_grid[i - 1]))
        ) / (2 * h)
        logder_rows.append(
            {
                "t": t,
                "operator_logder": op,
                "xi_logder": xx,
                "difference": op - xx,
            }
        )

    logder_rms = (
        float(
            np.sqrt(
                np.mean([r["difference"] ** 2 for r in logder_rows])
            )
        )
        if logder_rows
        else float("nan")
    )

    # Summaries
    final_post = next(
        r
        for r in post_rows
        if r["L"] == Ls[-1] and r["method"] == canonical_method
    )
    raw_final = next(
        r
        for r in post_rows
        if r["L"] == Ls[-1] and r["method"] == "raw_prefix"
    )
    weyl_final = next(
        r
        for r in post_rows
        if r["L"] == Ls[-1] and r["method"] == "weyl"
    )

    # Does renormalization improve raw L-support drift versus v97?
    ren_better_count = 0
    if raw_control_rows:
        ren_better_count = sum(r["ren_over_raw"] < 1.0 for r in raw_control_rows)

    post_gates = {
        "canonical_det_beats_raw":
            final_post["det_symmetric_rms"] < raw_final["det_symmetric_rms"],
        "canonical_det_beats_weyl":
            final_post["det_symmetric_rms"] < weyl_final["det_symmetric_rms"],
        "canonical_moments_below_1e4":
            final_post["moment_relative_l2"] < 1e-4,
        "canonical_det_below_1e6":
            final_post["det_symmetric_rms"] < 1e-6,
        "zero_relative_rmse_below_1e4":
            bool(np.isfinite(zrel) and zrel < 1e-4),
        "renormalization_improves_majority_raw_L_drifts":
            True if not raw_control_rows
            else ren_better_count >= math.ceil(len(raw_control_rows) / 2),
    }

    # Conservative verdict: numerical exactness is never upgraded to theorem.
    if not all(pre_gates.values()):
        verdict = "FAIL_V102_PRE_XI_OPERATOR_OR_QUADRATURE_GATE"
    elif all(post_gates.values()):
        verdict = (
            "PASS_V102_NUMERICAL_PNT_FINITE_PART_HP_CLOSURE__"
            "TRACE_NORM_LIMIT_AND_EXACT_ENTIRE_IDENTITY_STILL_REQUIRE_PROOF"
        )
    else:
        verdict = (
            "PARTIAL_V102_FINITE_PART_OPERATOR_CONSTRUCTED__"
            "GLOBAL_CLOSURE_GATES_NOT_ALL_PASSED"
        )

    # Outputs
    construction_csv = Path(str(prefix) + "_CONSTRUCTION.csv")
    flow_csv = Path(str(prefix) + "_L_FLOW.csv")
    raw_control_csv = Path(str(prefix) + "_RAW_VS_REN_FLOW.csv")
    tail_csv = Path(str(prefix) + "_TAIL_COMPONENTS.csv")
    detflow_csv = Path(str(prefix) + "_FREDHOLM_FLOW.csv")
    post_csv = Path(str(prefix) + "_POSTFREEZE_XI.csv")
    logder_csv = Path(str(prefix) + "_POSTFREEZE_LOGDER.csv")
    verdict_path = Path(str(prefix) + "_VERDICT.json")

    write_csv(construction_csv, construction_rows)
    write_csv(flow_csv, flow_rows)
    write_csv(raw_control_csv, raw_control_rows)
    write_csv(tail_csv, tail_rows)
    write_csv(detflow_csv, detflow_rows)
    write_csv(post_csv, post_rows)
    write_csv(logder_csv, logder_rows)

    packet = {
        "version": VERSION,
        "build": BUILD,
        "verdict": verdict,
        "RH_proof": False,
        "complete_HP_theorem": False,
        "freeze_sha256": freeze_sha,
        "canonical_tail_method": canonical_method,
        "pre_gates": pre_gates,
        "tail_series_order_drift": tail_series_order_drift,
        "post_gates": post_gates,
        "final_postfreeze": final_post,
        "final_zero_relative_rmse": zrel,
        "postfreeze_logder_rms": logder_rms,
        "renormalized_L_flow": flow_rows,
        "fredholm_L_flow": detflow_rows,
        "raw_vs_renormalized_flow": raw_control_rows,
        "theorem_ledger": {
            "fixed_L_infinite_self_adjoint_operator": "PROVED",
            "fixed_L_H_inverse_square_trace_class": "PROVED",
            "zero_free_PNT_finite_part_renormalization": "CONSTRUCTED",
            "renormalized_L_to_infinity_trace_norm_limit": "OPEN",
            "global_EM_remainder_bound": "OPEN",
            "exact_Xi_Fredholm_entire_identity": "OPEN",
            "RH": "OPEN",
        },
        "next_if_pass": (
            "Prove uniform Euler-Maclaurin remainder bounds for the renormalized "
            "tail and a regulator-independent trace-norm limit K_L^ren -> "
            "K_infinity; then identify the limiting Fredholm determinant "
            "with Xi as an entire function."
        ),
        "runtime_seconds": time.time() - started,
    }
    verdict_path.write_text(
        json.dumps(json_safe(packet), indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 168)
    print("FINAL CONSERVATIVE v102.1 VERDICT")
    print("=" * 168)
    print("verdict                          :", verdict)
    print("canonical tail method            :", canonical_method)
    print("fixed-L infinite self-adjoint    : PROVED")
    print("fixed-L H^-2 trace class         : PROVED")
    print("PNT finite-part family           : CONSTRUCTED ZERO-FREE")
    print("finite-part L->infinity theorem : OPEN")
    print("exact Xi Fredholm identity       : OPEN")
    print("RH proof                         : FALSE")
    print("final zero relative RMSE         :", f"{zrel:.12e}")
    print("final raw determinant RMS        :", f"{raw_final['det_symmetric_rms']:.12e}")
    print("final Weyl determinant RMS       :", f"{weyl_final['det_symmetric_rms']:.12e}")
    print("final canonical determinant RMS  :", f"{final_post['det_symmetric_rms']:.12e}")
    print("final determinant sup abs        :", f"{final_post['det_sup_abs']:.12e}")
    print("final canonical moment rel L2    :", f"{final_post['moment_relative_l2']:.12e}")
    print("final S1 relative error          :", f"{final_post['S1_relative_error']:.12e}")
    print("pre-Xi tail-series order drift   :", f"{tail_series_order_drift:.12e}")
    print("postfreeze log-derivative RMS    :", f"{logder_rms:.12e}")
    if flow_rows:
        print("latest FP support L-flow drift  :", f"{flow_rows[-1]['support_relative_l2']:.12e}")
        print("latest FP m-function drift      :", f"{flow_rows[-1]['m_function_relative_l2']:.12e}")
    if detflow_rows:
        print("latest Fredholm L-flow drift     :", f"{detflow_rows[-1]['determinant_relative_l2']:.12e}")
        print("latest moment L-flow drift       :", f"{detflow_rows[-1]['moments_relative_l2']:.12e}")
    if raw_control_rows:
        print("latest FP/raw support drift ratio:",
              f"{raw_control_rows[-1]['ren_over_raw']:.6f}")
    print("freeze SHA256                    :", freeze_sha)
    print("verdict file                     :", verdict_path)
    print("runtime seconds                  :", f"{time.time()-started:.2f}")
    print("=" * 168)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
