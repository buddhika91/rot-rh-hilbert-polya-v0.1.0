#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v97 — Full ROT -> Infinite Hilbert-Pólya -> Xi Determinant Closure Harness
================================================================================

PURPOSE
-------
This script implements the complete *theorem ladder* that is presently available
without pretending that the two genuinely open steps are solved.

CONSTRUCTION FIREWALL
---------------------
Before the freeze this file uses only:
  * Riemann-Siegel/Gamma phase;
  * von Mangoldt prime/prime-power arithmetic;
  * a ROT heat/self-observation metric on the spectral measure;
  * Stieltjes/Jacobi operator theory;
  * recursive observation scale L.

It does NOT evaluate Xi, zeta, zeta derivatives, known zero ordinates, or any
zero-based score before the operator family is serialized and hashed.

THEOREM LADDER
--------------
Stage A. Zero-free smooth arithmetic phase

    F_L(E) = theta(E)/pi + 1
             - (1/pi) sum_{n>=2} Lambda(n)/(log n sqrt(n))
                         exp(-n/L) sin(E log n).

The exponential regulator is the one used by the existing fixed-L infinite
operator theorem.  The numerical sum is truncated at n <= tail_factor * L and
a rigorous integral majorant for the omitted absolute tail is reported.

Stage B. Ordered arithmetic support

    E_1(L) < E_2(L) < ...,
    F_L(E_k(L)) = k - 1/2.

The numerical root finder implements the recursive "first later crossing"
definition rather than choosing a root from zero data.

Stage C. ROT heat spectral metric and finite Jacobi sections

    w_k^(tau,L) propto exp[-tau (E_k(L)/E_1(L))^2],

followed by Stieltjes/Lanczos to obtain

    J_(L,tau,N) = tridiagonal(a_n,b_n),   b_n > 0.

Stage D. Infinite fixed-L operator theorem

For every fixed finite L>0 and tau>0, the analytic support asymptotic
E_k(L) ~ 2*pi*k/log(k) plus the heat weight gives a determinate Hamburger
moment problem.  Therefore the minimal Jacobi operator is essentially
self-adjoint and has a unique closure

    H_ROT(L,tau) = closure(J_(L,tau)).

This theorem is *not* inferred from the finite numerics; the numerics merely
audit its finite sections.

Stage E. Compact inverse-square operator

    K_(L,tau) = H_ROT(L,tau)^(-2).

Since E_k(L) ~ 2*pi*k/log(k),

    sum_k E_k(L)^(-2) < infinity,

so K_(L,tau) is trace class.  Hence the ordinary Fredholm determinant exists:

    D_L(t) = det(I - t^2 K_(L,tau))
           = product_k (1 - t^2/E_k(L)^2).

The determinant depends on the spectral support; the heat weights select the
Jacobi/cyclic representation but do not by themselves change the spectrum.
This is recorded explicitly because ROT must become spectrally indispensable
through the recursive L-flow / support law, not merely through reweighting.

Stage F. Recursive-resolution / infinite-resolution audit

For a ladder L_1 < ... < L_m, the script tests only zero-free, gauge-invariant
objects:
  * Weyl/Stieltjes m-functions of the heat measures;
  * inverse even spectral traces;
  * finite canonical determinants;
  * root/support drift;
  * depth drift.

These are numerical diagnostics for the still-open theorem:

    H_ROT(L,tau) -> H_ROT,infinity     (strong/norm resolvent, suitably defined)

or, more directly for the determinant route,

    K_(L,tau) -> K_infinity            in trace norm.

No numerical threshold is promoted to a proof of this limit.

Stage G. HARD FREEZE

All roots, Jacobi coefficients, heat weights, m-function samples, inverse trace
moments, determinant samples, construction parameters and theorem ledger are
serialized and SHA256 hashed.

Stage H. Post-freeze Xi audit

Only after the freeze, evaluate

    Xi(t) = xi(1/2 + i t)

externally and compare the frozen operator determinant with

    Xi(t)/Xi(0).

Also compare post-freeze Xi sign-change roots and inverse logarithmic moments.
No post-freeze scale or parameter is fed back into construction.

Stage I. Exact finish line

The RH proof would require the two open mathematical statements:

  T_LIMIT:
      prove existence/regulator-independence of K_infinity (or H_infinity);

  T_IDENTITY:
      prove as entire functions

          det(I - t^2 K_infinity) = Xi(t)/Xi(0).

If K_infinity is positive compact injective, then

    H_HP = K_infinity^(-1/2)

is self-adjoint and unbounded.  The exact determinant identity then forces every
zero of Xi(t) to be real, which is RH.

This program is a closure harness, not a declaration that T_LIMIT or T_IDENTITY
has already been proved.
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
from scipy.optimize import brentq
from scipy.special import erfc

VERSION = "v97"
BUILD = "2026-08-17-full-rot-hp-xi-closure-v97"
EPS = 1e-14


def parse_ints(text: str) -> List[int]:
    return [int(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_complexes(text: str) -> List[complex]:
    out=[]
    for token in str(text).split(","):
        token=token.strip().replace("i","j")
        if token:
            out.append(complex(token))
    return out


def json_safe(x):
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x,np.ndarray):
        return x.tolist()
    if isinstance(x,(np.bool_,)):
        return bool(x)
    if isinstance(x,np.integer):
        return int(x)
    if isinstance(x,np.floating):
        return float(x)
    if isinstance(x,np.complexfloating):
        return {"real":float(np.real(x)),"imag":float(np.imag(x))}
    if isinstance(x,complex):
        return {"real":float(x.real),"imag":float(x.imag)}
    return x


def write_csv(path: Path, rows: List[Dict[str,Any]]) -> None:
    if not rows:
        path.write_text("",encoding="utf-8")
        return
    fields=[]; seen=set()
    for row in rows:
        for k in row:
            if k not in seen:
                fields.append(k); seen.add(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class SmoothSupport:
    L: int
    depth: int
    nmax: int
    roots: np.ndarray
    residuals: np.ndarray
    root_scan_counts: np.ndarray
    logs: np.ndarray
    weights: np.ndarray
    arithmetic_tail_bound_phase: float


@dataclass
class OperatorSection:
    support: SmoothSupport
    tau: float
    heat_weights: np.ndarray
    diagonal: np.ndarray
    offdiagonal: np.ndarray
    jacobi: np.ndarray
    self_adjoint_defect: float
    minimum_offdiagonal: float
    orthogonality_defect: float
    eigenvalue_preservation_error: float
    carleman_partial_sum: float
    inverse_traces: np.ndarray
    m_values: np.ndarray
    determinant_grid: np.ndarray


def build_exponential_terms(base,L:int,tail_factor:float):
    nmax=max(8,int(math.ceil(float(tail_factor)*float(L))))
    numbers,lambdas,_exponents=base.mangoldt_terms(nmax)
    nums=numbers.astype(float)
    logs=np.log(nums)
    # Lambda(n)/(log n sqrt n) * exp(-n/L)
    weights=lambdas/(logs*np.sqrt(nums))*np.exp(-nums/float(L))

    # Because Lambda(n)/log(n) <= 1 for prime powers, the omitted absolute
    # arithmetic phase is bounded by sum_{n>N} exp(-n/L)/sqrt(n)/pi.
    # Integral majorant plus first omitted lattice term.
    N=float(nmax)
    integral=math.sqrt(math.pi*float(L))*float(erfc(math.sqrt(N/float(L))))
    lattice=math.exp(-(N+1.0)/float(L))/math.sqrt(N+1.0)
    phase_tail=(integral+lattice)/math.pi
    return nmax,logs,weights,float(phase_tail)


def phase_correction(E,logs,weights,chunk=20000):
    e=np.atleast_1d(np.asarray(E,float))
    out=np.zeros(len(e),float)
    for j in range(0,len(logs),chunk):
        ll=logs[j:j+chunk]; ww=weights[j:j+chunk]
        out += np.sin(np.outer(e,ll))@ww
    return -out/math.pi


def F_values(base,E,logs,weights):
    e=np.atleast_1d(np.asarray(E,float))
    smooth=np.asarray([base.smooth_count(float(x)) for x in e],float)
    return smooth+phase_correction(e,logs,weights)


def ordered_support(base,L:int,depth:int,tail_factor:float,scan_points:int) -> SmoothSupport:
    nmax,logs,weights,tail_bound=build_exponential_terms(base,L,tail_factor)
    free=base.free_archimedean_levels(depth+2)
    roots=[]; residuals=[]; scans=[]
    previous=7.0

    for idx in range(depth):
        target=float(idx+1)-0.5
        lower=7.0 if idx==0 else roots[-1]+1e-9
        # Start with an Archimedean-informed upper bound, then expand until a
        # sign change is found.  We take the FIRST crossing after the prior root.
        nominal=float(free[min(idx+1,len(free)-1)])
        upper=max(lower+1.0,nominal+0.75*max(1.0,nominal-lower))
        found=None; total_intervals=0
        for _expand in range(12):
            grid=np.linspace(lower,upper,int(scan_points))
            vals=F_values(base,grid,logs,weights)-target
            cross=np.flatnonzero(vals[:-1]*vals[1:]<=0.0)
            total_intervals += len(grid)-1
            if len(cross):
                q=int(cross[0])
                a=float(grid[q]); b=float(grid[q+1])
                root=brentq(
                    lambda x: float(F_values(base,np.asarray([x]),logs,weights)[0]-target),
                    a,b,xtol=1e-12,rtol=1e-12,maxiter=300,
                )
                found=float(root)
                break
            lower=upper
            upper=max(upper+2.0,upper*1.25)
        if found is None:
            raise RuntimeError(f"Failed to locate ordered crossing k={idx+1} at L={L}.")
        roots.append(found)
        resid=abs(float(F_values(base,np.asarray([found]),logs,weights)[0]-target))
        residuals.append(resid); scans.append(total_intervals)

    return SmoothSupport(
        L=int(L),depth=int(depth),nmax=int(nmax),
        roots=np.asarray(roots,float),residuals=np.asarray(residuals,float),
        root_scan_counts=np.asarray(scans,int),logs=logs,weights=weights,
        arithmetic_tail_bound_phase=float(tail_bound),
    )


def heat_weights(roots:np.ndarray,tau:float)->np.ndarray:
    r=np.asarray(roots,float)
    x=(r/r[0])**2
    logw=-float(tau)*x
    logw-=np.max(logw)
    w=np.exp(logw)
    w/=np.sum(w)
    return w


def m_function(roots,weights,z):
    r=np.asarray(roots,float); w=np.asarray(weights,float)
    return np.sum(w/(r-complex(z)))


def finite_determinant(roots,tgrid):
    r=np.asarray(roots,float)
    out=[]
    inv2=1.0/(r*r)
    for t in np.asarray(tgrid,float):
        out.append(float(np.prod(1.0-(t*t)*inv2)))
    return np.asarray(out,float)


def build_section(base,support:SmoothSupport,tau:float,m_probes,tgrid,moment_order:int)->OperatorSection:
    w=heat_weights(support.roots,tau)
    a,b,diag=base.jacobi_from_atoms(support.roots,w)
    J=np.diag(a)
    if len(b):
        J += np.diag(b,1)+np.diag(b,-1)
    inv=np.asarray([
        float(np.sum(support.roots**(-2*n))) for n in range(1,moment_order+1)
    ])
    mv=np.asarray([m_function(support.roots,w,z) for z in m_probes],complex)
    det=finite_determinant(support.roots,tgrid)
    return OperatorSection(
        support=support,tau=float(tau),heat_weights=w,
        diagonal=np.asarray(a,float),offdiagonal=np.asarray(b,float),jacobi=J,
        self_adjoint_defect=float(np.linalg.norm(J-J.T,"fro")),
        minimum_offdiagonal=float(np.min(b)) if len(b) else float("inf"),
        orthogonality_defect=float(diag["orthogonality_defect"]),
        eigenvalue_preservation_error=float(diag["eigenvalue_preservation_error"]),
        carleman_partial_sum=float(np.sum(1.0/np.maximum(np.asarray(b,float),EPS))) if len(b) else 0.0,
        inverse_traces=inv,m_values=mv,determinant_grid=det,
    )


def rel_l2(a,b):
    a=np.asarray(a); b=np.asarray(b)
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),EPS))


def symmetric_rms(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.sqrt(np.mean(((a-b)/(1.0+np.abs(a)+np.abs(b)))**2)))


def xi_postfreeze(t,mp):
    s=mp.mpf("0.5")+1j*mp.mpf(t)
    val=mp.mpf("0.5")*s*(s-1)*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
    return mp.re(val)


def scan_xi_roots(mp,tmax,step):
    roots=[]
    x0=mp.mpf("0.0"); f0=xi_postfreeze(x0,mp)
    x=float(step)
    while x<=float(tmax)+1e-12:
        x1=mp.mpf(str(x)); f1=xi_postfreeze(x1,mp)
        if f0==0:
            roots.append(float(x0))
        elif f0*f1<0:
            r=mp.findroot(lambda q: xi_postfreeze(q,mp),(x0,x1))
            rr=float(r)
            if not roots or abs(rr-roots[-1])>1e-8:
                roots.append(rr)
        x0=x1; f0=f1; x+=float(step)
    return np.asarray(roots,float)


def xi_log_moments(mp,order):
    x0=xi_postfreeze(mp.mpf("0"),mp)
    vals=[]
    for n in range(1,int(order)+1):
        deriv=mp.diff(lambda t: mp.log(xi_postfreeze(t,mp)/x0),mp.mpf("0"),2*n)
        coeff=deriv/mp.factorial(2*n)
        vals.append(float(-n*coeff))
    return np.asarray(vals,float)


def main():
    ap=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--base-module",default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit")
    ap.add_argument("--L-list",default="2000,4000,8000,16000")
    ap.add_argument("--depths",default="24,32,48")
    ap.add_argument("--tau",type=float,default=0.02)
    ap.add_argument("--tail-factor",type=float,default=18.0)
    ap.add_argument("--scan-points",type=int,default=96)
    ap.add_argument("--moment-order",type=int,default=6)
    ap.add_argument("--m-probes",default="1j,10+1j,20+2j,40+4j")
    ap.add_argument("--identity-max",type=float,default=50.0)
    ap.add_argument("--identity-step",type=float,default=0.25)
    ap.add_argument("--xi-dps",type=int,default=80)
    ap.add_argument("--zero-scan-max",type=float,default=120.0)
    ap.add_argument("--zero-scan-step",type=float,default=0.10)
    ap.add_argument("--zero-match-count",type=int,default=32)
    ap.add_argument("--maximum-tail-phase-bound",type=float,default=1e-6)
    ap.add_argument("--maximum-root-residual",type=float,default=1e-8)
    ap.add_argument("--maximum-selfadjoint-defect",type=float,default=1e-12)
    ap.add_argument("--minimum-offdiagonal",type=float,default=1e-12)
    ap.add_argument("--prefix",default="rot_rh_full_rot_hp_pipeline_v97")
    args=ap.parse_args()

    started=time.time()
    Ls=parse_ints(args.L_list); depths=parse_ints(args.depths)
    if Ls!=sorted(Ls) or depths!=sorted(depths):
        raise ValueError("L-list and depths must be increasing.")
    if min(Ls)<=0 or min(depths)<4 or args.tau<=0:
        raise ValueError("Need L>0, depth>=4, tau>0.")

    base=importlib.import_module(args.base_module)
    probes=parse_complexes(args.m_probes)
    tgrid=np.arange(0.0,args.identity_max+0.5*args.identity_step,args.identity_step)
    prefix=Path(args.prefix).expanduser().resolve(); prefix.parent.mkdir(parents=True,exist_ok=True)

    print("="*164)
    print("ROT-RH v97 — FULL ROT -> INFINITE HP -> Xi DETERMINANT CLOSURE HARNESS")
    print("="*164)
    print("L ladder                      :",Ls)
    print("depth ladder                  :",depths)
    print("ROT heat tau                  :",args.tau)
    print("exp-regulator tail factor     :",args.tail_factor)
    print("construction                  : Gamma + Mangoldt + ROT heat metric")
    print("Xi/zeta/known zeros pre-freeze: NONE")
    print("="*164)

    sections={}; construct_rows=[]
    print("[A-C] Build zero-free supports and ROT heat-Jacobi sections",flush=True)
    for L in Ls:
        # Ordered supports are nested by definition. Build the largest prefix once
        # and slice it for the requested depth ladder; this is faster and also
        # makes the depth-nesting audit exact at the support-construction level.
        sup_max=ordered_support(base,L,depths[-1],args.tail_factor,args.scan_points)
        for d in depths:
            sup=SmoothSupport(
                L=sup_max.L, depth=int(d), nmax=sup_max.nmax,
                roots=sup_max.roots[:d].copy(),
                residuals=sup_max.residuals[:d].copy(),
                root_scan_counts=sup_max.root_scan_counts[:d].copy(),
                logs=sup_max.logs, weights=sup_max.weights,
                arithmetic_tail_bound_phase=sup_max.arithmetic_tail_bound_phase,
            )
            sec=build_section(base,sup,args.tau,probes,tgrid,args.moment_order)
            sections[(L,d)]=sec
            row={
                "L":L,"depth":d,"nmax":sup.nmax,
                "tail_phase_bound":sup.arithmetic_tail_bound_phase,
                "max_root_residual":float(np.max(sup.residuals)),
                "first_root":float(sup.roots[0]),"last_root":float(sup.roots[-1]),
                "self_adjoint_defect":sec.self_adjoint_defect,
                "minimum_offdiagonal":sec.minimum_offdiagonal,
                "orthogonality_defect":sec.orthogonality_defect,
                "eigenvalue_preservation_error":sec.eigenvalue_preservation_error,
                "carleman_partial_sum":sec.carleman_partial_sum,
                "trace_K_partial":float(sec.inverse_traces[0]),
            }
            construct_rows.append(row)
            print(
                f"  L={L:<7d} d={d:<3d} roots=[{sup.roots[0]:.6f},{sup.roots[-1]:.6f}] "
                f"tail<={sup.arithmetic_tail_bound_phase:.2e} maxres={np.max(sup.residuals):.2e} "
                f"minb={sec.minimum_offdiagonal:.2e}"
            )

    print("[D-F] Fixed-L theorem certificates + recursive-resolution diagnostics",flush=True)
    depth_rows=[]; cutoff_rows=[]
    dmax=depths[-1]
    for L in Ls:
        for d0,d1 in zip(depths[:-1],depths[1:]):
            s0=sections[(L,d0)]; s1=sections[(L,d1)]
            n=min(d0,d1)
            depth_rows.append({
                "L":L,"depth_from":d0,"depth_to":d1,
                "support_prefix_relative_l2":rel_l2(s1.support.roots[:n],s0.support.roots[:n]),
                "m_function_relative_l2":rel_l2(s1.m_values,s0.m_values),
                "inverse_trace_relative_l2":rel_l2(s1.inverse_traces,s0.inverse_traces),
                "determinant_relative_l2":rel_l2(s1.determinant_grid,s0.determinant_grid),
            })

    for L0,L1 in zip(Ls[:-1],Ls[1:]):
        s0=sections[(L0,dmax)]; s1=sections[(L1,dmax)]
        cutoff_rows.append({
            "L_from":L0,"L_to":L1,
            "support_relative_l2":rel_l2(s1.support.roots,s0.support.roots),
            "m_function_relative_l2":rel_l2(s1.m_values,s0.m_values),
            "inverse_trace_relative_l2":rel_l2(s1.inverse_traces,s0.inverse_traces),
            "determinant_relative_l2":rel_l2(s1.determinant_grid,s0.determinant_grid),
        })
        print(
            f"  L {L0}->{L1}: support={cutoff_rows[-1]['support_relative_l2']:.3e} "
            f"m={cutoff_rows[-1]['m_function_relative_l2']:.3e} "
            f"traces={cutoff_rows[-1]['inverse_trace_relative_l2']:.3e}"
        )

    final=sections[(Ls[-1],dmax)]

    # Analytic theorem ledger BEFORE Xi access.
    pre_gates={
        "finite_support_strictly_ordered":bool(np.all(np.diff(final.support.roots)>0)),
        "numerical_tail_bound_small":bool(final.support.arithmetic_tail_bound_phase<=args.maximum_tail_phase_bound),
        "root_residuals_small":bool(np.max(final.support.residuals)<=args.maximum_root_residual),
        "finite_jacobi_self_adjoint":bool(final.self_adjoint_defect<=args.maximum_selfadjoint_defect),
        "finite_jacobi_positive_offdiagonal":bool(final.minimum_offdiagonal>=args.minimum_offdiagonal),
    }

    theorem_ledger_pre={
        "T_ARITHMETIC_PHASE":{
            "status":"PROVED_FOR_FIXED_L",
            "statement":"Exponential Mangoldt series is absolutely/uniformly convergent on real E for finite L; F_L is continuous and unbounded above.",
        },
        "T_ORDERED_SUPPORT":{
            "status":"PROVED_FOR_FIXED_L",
            "statement":"Recursive first-later-crossing construction gives an infinite strictly ordered support E_k(L)->infinity.",
        },
        "T_FIXED_L_SELF_ADJOINT":{
            "status":"PROVED_FOR_FIXED_L_TAU",
            "statement":"ROT heat measure has exponential moments and determinate Hamburger moment problem; minimal Jacobi operator is essentially self-adjoint.",
        },
        "T_TRACE_CLASS_INVERSE_SQUARE":{
            "status":"PROVED_FOR_FIXED_L",
            "statement":"E_k(L)~2*pi*k/log(k) implies sum E_k(L)^(-2)<infinity; K_L=H_L^(-2) is trace class and det(I-t^2 K_L) exists.",
        },
        "T_ROT_SPECTRAL_INDISPENSABILITY":{
            "status":"OPEN",
            "statement":"Heat weights change the cyclic/Jacobi representation but not the support. A first-principles ROT recursive law must act on the spectral support/operator geometry and survive controls.",
        },
        "T_L_TO_INFINITY_OPERATOR_LIMIT":{
            "status":"OPEN__NUMERICAL_DIAGNOSTICS_ONLY",
            "statement":"Prove a regulator-independent strong/norm-resolvent H_L limit, or trace-norm K_L limit, as recursive observation resolution L->infinity.",
        },
        "T_ROT_THETA_BRIDGE":{
            "status":"OPEN_FROM_FUNDAMENTAL_ROT_ACTION",
            "statement":"Derive square heat spectrum and modular inversion/self-duality from ROT rather than inserting the theta kernel as an assumption.",
        },
    }

    # HARD FREEZE before mpmath / Xi import.
    archive=Path(str(prefix)+"_FROZEN_PRE_XI.npz")
    manifest_path=Path(str(prefix)+"_FROZEN_MANIFEST.json")
    payload={
        "final_roots":final.support.roots,
        "final_heat_weights":final.heat_weights,
        "final_diagonal":final.diagonal,
        "final_offdiagonal":final.offdiagonal,
        "final_inverse_traces":final.inverse_traces,
        "final_m_values_real":np.real(final.m_values),
        "final_m_values_imag":np.imag(final.m_values),
        "tgrid":tgrid,
        "final_determinant":final.determinant_grid,
        "L_values":np.asarray(Ls,int),
        "depth_values":np.asarray(depths,int),
    }
    for (L,d),sec in sections.items():
        payload[f"roots_L{L}_d{d}"]=sec.support.roots
        payload[f"a_L{L}_d{d}"]=sec.diagonal
        payload[f"b_L{L}_d{d}"]=sec.offdiagonal
        payload[f"w_L{L}_d{d}"]=sec.heat_weights
        payload[f"traces_L{L}_d{d}"]=sec.inverse_traces
    np.savez_compressed(archive,**payload)
    frozen_sha=sha256_file(archive)
    manifest={
        "version":VERSION,"build":BUILD,"archive":str(archive),"sha256":frozen_sha,
        "firewall":{
            "Xi_used_pre_freeze":False,"zeta_used_pre_freeze":False,
            "known_zero_ordinates_used_pre_freeze":False,"zero_loss_used":False,
        },
        "parameters":{
            "L_list":Ls,"depths":depths,"tau":args.tau,"tail_factor":args.tail_factor,
            "scan_points":args.scan_points,"moment_order":args.moment_order,
            "m_probes":[str(z) for z in probes],"identity_max":args.identity_max,
            "identity_step":args.identity_step,
        },
        "pre_gates":pre_gates,"theorem_ledger_pre":theorem_ledger_pre,
    }
    manifest_path.write_text(json.dumps(json_safe(manifest),indent=2),encoding="utf-8")
    print("[G] HARD FREEZE SHA256          :",frozen_sha)
    print("[G] Xi firewall CLOSED. Importing Xi tools only now.")

    # ------------------------------------------------------------------
    # POST-FREEZE Xi / theorem-target diagnostics only.
    # ------------------------------------------------------------------
    import mpmath as mp
    mp.mp.dps=int(args.xi_dps)

    xi0=xi_postfreeze(mp.mpf("0"),mp)
    xi_vals=np.asarray([float(xi_postfreeze(mp.mpf(str(t)),mp)/xi0) for t in tgrid],float)
    det_rms=symmetric_rms(final.determinant_grid,xi_vals)

    xi_roots=scan_xi_roots(mp,args.zero_scan_max,args.zero_scan_step)
    nmatch=min(args.zero_match_count,len(xi_roots),len(final.support.roots))
    if nmatch:
        zr=final.support.roots[:nmatch]; xt=xi_roots[:nmatch]
        zero_rel_rmse=float(np.sqrt(np.mean(((zr-xt)/xt)**2)))
        zero_rmse=float(np.sqrt(np.mean((zr-xt)**2)))
    else:
        zero_rel_rmse=float("nan"); zero_rmse=float("nan")

    xi_mom=xi_log_moments(mp,args.moment_order)
    # The finite support gives partial traces; convergence in depth is separately audited.
    moment_rel=rel_l2(final.inverse_traces,xi_mom)

    # Log-derivative/resolvent identity away from finite operator roots.
    logder_rows=[]
    for t in np.arange(0.5,args.identity_max+1e-12,max(0.5,args.identity_step*2.0)):
        if np.min(np.abs(final.support.roots-t))<0.20:
            continue
        op=-2.0*t*float(np.sum(1.0/(final.support.roots**2-t*t)))
        h=max(mp.mpf("1e-20"),mp.mpf("1e-8")*(1+abs(t)))
        xv=mp.mpf(str(t))
        xi_logder=float(mp.re(mp.diff(lambda q: mp.log(xi_postfreeze(q,mp)/xi0),xv)))
        logder_rows.append({"t":t,"operator_log_derivative":op,"xi_log_derivative":xi_logder,"difference":op-xi_logder})
    logder_rms=float(np.sqrt(np.mean([r["difference"]**2 for r in logder_rows]))) if logder_rows else float("nan")

    theorem_ledger_post=dict(theorem_ledger_pre)
    theorem_ledger_post.update({
        "T_EXACT_XI_DETERMINANT":{
            "status":"OPEN__CURRENT_NUMERICAL_TEST_REPORTED",
            "statement":"Prove det(I-t^2 K_infinity)=Xi(t)/Xi(0) as an entire-function identity; finite-grid agreement alone is never a proof.",
            "finite_symmetric_rms":det_rms,
        },
        "T_RH":{
            "status":"OPEN",
            "statement":"Would follow once K_infinity is positive compact injective and the exact Xi determinant identity is proved, because H=K_infinity^(-1/2) is self-adjoint and all determinant zeros are real.",
        },
    })

    post={
        "xi_root_count_scanned":int(len(xi_roots)),
        "zero_match_count":int(nmatch),
        "zero_rmse":zero_rmse,"zero_relative_rmse":zero_rel_rmse,
        "determinant_symmetric_rms":det_rms,
        "inverse_moment_relative_l2":moment_rel,
        "operator_inverse_traces":final.inverse_traces,
        "xi_inverse_log_moments":xi_mom,
        "log_derivative_rms":logder_rms,
    }

    # Conservative final verdict: numerical data cannot turn open theorems into proofs.
    if not all(pre_gates.values()):
        verdict="FAIL_V97_PRE_XI_FINITE_OR_NUMERICAL_CERTIFICATE"
    elif det_rms<1e-6 and zero_rel_rmse<1e-6:
        verdict="PARTIAL_V97_EXCEPTIONAL_POSTFREEZE_MATCH__LIMIT_AND_EXACT_IDENTITY_STILL_REQUIRE_PROOF"
    else:
        verdict="PARTIAL_V97_FIXED_L_INFINITE_OPERATOR_COMPLETE__GLOBAL_ROT_LIMIT_AND_XI_IDENTITY_OPEN"

    construct_csv=Path(str(prefix)+"_CONSTRUCTION.csv")
    depth_csv=Path(str(prefix)+"_DEPTH_CONVERGENCE.csv")
    cutoff_csv=Path(str(prefix)+"_RESOLUTION_FLOW.csv")
    logder_csv=Path(str(prefix)+"_POSTFREEZE_LOGDER.csv")
    verdict_path=Path(str(prefix)+"_VERDICT.json")
    write_csv(construct_csv,construct_rows); write_csv(depth_csv,depth_rows)
    write_csv(cutoff_csv,cutoff_rows); write_csv(logder_csv,logder_rows)

    packet={
        "version":VERSION,"build":BUILD,"verdict":verdict,"RH_proof":False,
        "frozen_sha256":frozen_sha,"pre_gates":pre_gates,
        "theorem_ledger":theorem_ledger_post,"postfreeze":post,
        "latest_zero_free_resolution_flow":cutoff_rows[-1] if cutoff_rows else None,
        "outputs":{
            "frozen_archive":str(archive),"manifest":str(manifest_path),
            "construction":str(construct_csv),"depth_convergence":str(depth_csv),
            "resolution_flow":str(cutoff_csv),"postfreeze_logder":str(logder_csv),
            "verdict":str(verdict_path),
        },
        "runtime_seconds":time.time()-started,
    }
    verdict_path.write_text(json.dumps(json_safe(packet),indent=2),encoding="utf-8")

    print("\n"+"="*164)
    print("FINAL CONSERVATIVE THEOREM-LADDER VERDICT")
    print("="*164)
    print("verdict                         :",verdict)
    print("fixed-L infinite self-adjoint  : PROVED (analytic theorem for each finite L,tau)")
    print("K_L = H_L^-2 trace class       : PROVED (from E_k ~ 2*pi*k/log k)")
    print("ROT spectral indispensability  : OPEN")
    print("L -> infinity operator limit   : OPEN")
    print("exact Xi determinant identity  : OPEN")
    print("RH proof                        : FALSE")
    print("postfreeze zero relative RMSE  :",f"{zero_rel_rmse:.9e}")
    print("postfreeze determinant sym RMS :",f"{det_rms:.9e}")
    print("postfreeze inverse moment rel  :",f"{moment_rel:.9e}")
    print("postfreeze log-derivative RMS  :",f"{logder_rms:.9e}")
    if cutoff_rows:
        last=cutoff_rows[-1]
        print("latest L-flow support drift     :",f"{last['support_relative_l2']:.9e}")
        print("latest L-flow m-function drift  :",f"{last['m_function_relative_l2']:.9e}")
        print("latest L-flow trace drift       :",f"{last['inverse_trace_relative_l2']:.9e}")
    print("freeze SHA256                   :",frozen_sha)
    print("verdict file                    :",verdict_path)
    print("runtime seconds                 :",f"{time.time()-started:.2f}")
    print("="*164)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
