#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v100 — Full Arithmetic Tail Completion of Frozen Fredholm Determinant
============================================================================

Goal
----
Test whether the SAME zero-free mechanism that repaired S1 in v99 also repairs
higher inverse moments and the full finite-window determinant.

Nothing is tuned from Xi.

Input
-----
The exact frozen v97 operator prefix:
    roots_L{L}_d{d}

Zero-free tail law
------------------
For p = 2m,

    T_m = integral_Ec^inf E^(-p) F'_L(E) dE
          - 1/2 Ec^(-p),

where

    F'_L(E)
      = theta'(E)/pi
        - (1/pi) sum_n w_n log(n) cos(E log n),

and

    w_n = Lambda(n)/(log n sqrt(n)) exp(-n/L).

The continuum piece is decomposed as

    Weyl + exact-Archimedean correction + prime oscillatory correction

and the final -1/2 Ec^(-p) is the same half-level/midpoint correction that
v99 tested for m=1.

For each moment order m this script computes four frozen ablations:

    A: Weyl
    B: Weyl + midpoint
    C: Weyl + exact Archimedean + midpoint
    D: Weyl + exact Archimedean + prime + midpoint   [FULL]

The completed determinant is

    D(t)
      = D_N(t) exp[-sum_{m=1}^M T_m t^(2m)/m].

All four completions are constructed and frozen BEFORE Xi is imported.

Post-freeze tests
-----------------
Compare all four frozen completions against:
  * Xi(t)/Xi(0) over the fixed grid;
  * first few Xi inverse log-moments.

No operator, tail coefficient, moment order, or ablation is selected by Xi.
"""

from __future__ import annotations

import argparse, csv, hashlib, importlib, json, math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.integrate import quad

EPS = 1e-15
VERSION = "v100.1"
BUILD = "2026-08-17-full-arithmetic-tail-fredholm-v100.1-stable-oscillatory-integrals"


def parse_ints(s):
    return [int(x.strip()) for x in str(s).split(",") if x.strip()]


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def json_safe(x):
    if isinstance(x,dict):
        return {str(k):json_safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x,np.ndarray):
        return x.tolist()
    if isinstance(x,(np.bool_,)):
        return bool(x)
    if isinstance(x,np.integer):
        return int(x)
    if isinstance(x,np.floating):
        return float(x)
    return x


def write_csv(path,rows):
    if not rows:
        Path(path).write_text("",encoding="utf-8"); return
    fields=[]; seen=set()
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k); seen.add(k)
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def partial_moments(roots,order):
    r=np.asarray(roots,float)
    return np.asarray([float(np.sum(r**(-2*m))) for m in range(1,order+1)])


def weyl_tail_moments(E,order):
    E=float(E)
    vals=[]
    for m in range(1,order+1):
        q=2*m-1
        vals.append(
            E**(1-2*m)/(2*math.pi*q)
            *(math.log(E/(2*math.pi))+1.0/q)
        )
    return np.asarray(vals,float)


def cosine_tail_integral_even(a,E,p,terms=12):
    """
    Stable large-aE expansion for

        I_p(a,E) = integral_E^inf cos(a x) / x^p dx.

    In this audit a >= log(2) and E >= ~100, hence aE >= ~69
    (and ~97 for the current d=48 endpoint).  The previous upward
    recurrence suffers catastrophic cancellation at high p.

    Repeated integration by parts gives

      integral_E^inf exp(i a x) x^-p dx
        ~ exp(i a E) * sum_{j>=0}
            i(-i)^j (p)_j E^(-p-j) a^(-j-1),

    where (p)_j is the rising factorial.  Taking the real part gives
    the desired cosine tail.  This is an asymptotic expansion in 1/(aE).

    We truncate after a fixed number of terms; no Xi information enters.
    """
    a=np.asarray(a,float)
    E=float(E)
    p=int(p)
    x=a*E
    phase=np.exp(1j*x)

    # j=0
    rising=np.ones_like(a)
    series=np.zeros_like(a,dtype=np.complex128)
    for j in range(int(terms)):
        if j>0:
            rising *= (p+j-1)
        coeff = 1j * ((-1j)**j)
        term = coeff * rising * (E**(-p-j)) * (a**(-j-1))
        series += term
    return np.real(phase*series)


def prime_tail_moments(logs,phase_w,E,order):
    coeff=-(phase_w*logs)/math.pi
    vals=[]
    for m in range(1,order+1):
        I=cosine_tail_integral_even(logs,E,2*m)
        vals.append(float(np.sum(coeff*I)))
    return np.asarray(vals,float)


def arch_tail_moments(base,E,order):
    vals=[]
    for m in range(1,order+1):
        p=2*m
        def fn(x):
            exact=base.theta_prime(float(x))/math.pi
            leading=math.log(float(x)/(2*math.pi))/(2*math.pi)
            return (exact-leading)/(float(x)**p)
        vals.append(float(
            quad(fn,float(E),np.inf,epsabs=1e-14,epsrel=1e-11,limit=250)[0]
        ))
    return np.asarray(vals,float)


def midpoint_moments(E,order):
    E=float(E)
    return np.asarray([-0.5*E**(-2*m) for m in range(1,order+1)],float)


def finite_det(roots,tgrid):
    r=np.asarray(roots,float)
    inv2=1.0/(r*r)
    return np.asarray([float(np.prod(1.0-(t*t)*inv2)) for t in tgrid])


def completed_det(raw,tgrid,tail):
    out=np.zeros_like(raw)
    for i,t in enumerate(tgrid):
        logtail=0.0
        for m,T in enumerate(tail,start=1):
            logtail -= float(T)*(float(t)**(2*m))/m
        if not math.isfinite(logtail) or abs(logtail) > 100.0:
            raise RuntimeError(
                f"Unphysical tail exponent at t={t}: logtail={logtail}. "
                "This indicates numerical instability in a tail coefficient."
            )
        out[i]=raw[i]*math.exp(logtail)
    return out


def rel_l2(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),EPS))


def sym_rms(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.sqrt(np.mean(((a-b)/(1+np.abs(a)+np.abs(b)))**2)))


def xi(mp,t):
    s=mp.mpf("0.5")+1j*mp.mpf(str(t))
    return mp.re(
        mp.mpf("0.5")*s*(s-1)
        *mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
    )


def xi_moments(mp,order):
    x0=xi(mp,0)
    out=[]
    for m in range(1,order+1):
        print(f"[Xi moment] {m}/{order}",flush=True)
        d=mp.diff(lambda t:mp.log(xi(mp,t)/x0),mp.mpf("0"),2*m)
        out.append(float(-m*d/mp.factorial(2*m)))
    return np.asarray(out,float)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--v97-prefix",default="rot_rh_full_rot_hp_pipeline_v97")
    ap.add_argument("--base-module",default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit")
    ap.add_argument("--L",type=int,default=16000)
    ap.add_argument("--depths",default="24,32,48")
    ap.add_argument("--tail-factor",type=float,default=18.0)
    ap.add_argument("--tail-moment-order",type=int,default=12)
    ap.add_argument("--post-moment-order",type=int,default=4)
    ap.add_argument("--identity-max",type=float,default=50.0)
    ap.add_argument("--identity-step",type=float,default=0.5)
    ap.add_argument("--xi-dps",type=int,default=60)
    ap.add_argument("--prefix",default="rot_rh_full_arithmetic_tail_v100")
    args=ap.parse_args()

    vp=Path(args.v97_prefix).expanduser().resolve()
    archive=Path(str(vp)+"_FROZEN_PRE_XI.npz")
    manifest=Path(str(vp)+"_FROZEN_MANIFEST.json")
    man=json.loads(manifest.read_text(encoding="utf-8"))
    sha=sha256_file(archive)
    if sha!=man["sha256"]:
        raise RuntimeError("v97 frozen SHA mismatch.")

    base=importlib.import_module(args.base_module)
    data=np.load(archive,allow_pickle=False)
    depths=parse_ints(args.depths)

    # Reconstruct ONLY the frozen arithmetic density used by v97.
    nmax=max(8,int(math.ceil(args.tail_factor*args.L)))
    numbers,lambdas,_=base.mangoldt_terms(nmax)
    nums=numbers.astype(float)
    logs=np.log(nums)
    phase_w=lambdas/(logs*np.sqrt(nums))*np.exp(-nums/float(args.L))

    tgrid=np.arange(
        0.0,args.identity_max+0.5*args.identity_step,args.identity_step
    )

    print("="*156)
    print("ROT-RH v100.1 — FULL ARITHMETIC TAIL COMPLETION (STABLE OSCILLATORY INTEGRALS)")
    print("="*156)
    print("v97 SHA verified              :",True)
    print("L                             :",args.L)
    print("depths                        :",depths)
    print("tail moment order             :",args.tail_moment_order)
    print("post-freeze Xi moment order   :",args.post_moment_order)
    print("Xi/zeta/known zeros pre-freeze: NONE")
    print("="*156)

    payload={"tgrid":tgrid,"depths":np.asarray(depths,int)}
    pre_rows=[]
    completed={}
    moments={}

    print("[1] Zero-free arithmetic tail construction",flush=True)
    for d in depths:
        roots=np.asarray(data[f"roots_L{args.L}_d{d}"],float)
        E=float(roots[-1])
        if args.identity_max>=E:
            raise ValueError("identity-max must remain below final prefix level.")

        partial=partial_moments(roots,args.tail_moment_order)
        W=weyl_tail_moments(E,args.tail_moment_order)
        A=arch_tail_moments(base,E,args.tail_moment_order)
        P=prime_tail_moments(logs,phase_w,E,args.tail_moment_order)
        M=midpoint_moments(E,args.tail_moment_order)

        tails={
            "weyl":W,
            "weyl_midpoint":W+M,
            "arch_midpoint":W+A+M,
            "full":W+A+P+M,
        }
        raw=finite_det(roots,tgrid)
        completed[d]={}
        moments[d]={}
        for name,T in tails.items():
            completed[d][name]=completed_det(raw,tgrid,T)
            moments[d][name]=partial+T

        payload[f"roots_d{d}"]=roots
        payload[f"partial_moments_d{d}"]=partial
        payload[f"weyl_tail_d{d}"]=W
        payload[f"arch_tail_d{d}"]=A
        payload[f"prime_tail_d{d}"]=P
        payload[f"midpoint_d{d}"]=M
        payload[f"raw_det_d{d}"]=raw
        for name in tails:
            payload[f"{name}_moments_d{d}"]=moments[d][name]
            payload[f"{name}_det_d{d}"]=completed[d][name]

        pre_rows.append({
            "depth":d,
            "Ecut":E,
            "prefix_S1":partial[0],
            "weyl_S1":W[0],
            "arch_S1":A[0],
            "prime_S1":P[0],
            "midpoint_S1":M[0],
            "full_tail_S1":tails["full"][0],
            "full_total_S1":moments[d]["full"][0],
        })
        print(
            f"  d={d:<3d} Ecut={E:.6f} "
            f"W={W[0]:+.9e} A={A[0]:+.3e} "
            f"P={P[0]:+.9e} M={M[0]:+.9e} "
            f"S1={moments[d]['full'][0]:.12e}"
        )

    print("[2] Zero-free depth convergence of FULL completion")
    conv=[]
    for d0,d1 in zip(depths[:-1],depths[1:]):
        r={
            "from":d0,"to":d1,
            "moment_rel_l2":rel_l2(moments[d1]["full"],moments[d0]["full"]),
            "det_rel_l2":rel_l2(completed[d1]["full"],completed[d0]["full"]),
        }
        conv.append(r)
        print(
            f"  {d0}->{d1}: moments={r['moment_rel_l2']:.3e} "
            f"det={r['det_rel_l2']:.3e}"
        )

    # HARD FREEZE
    prefix=Path(args.prefix).expanduser().resolve()
    freeze=Path(str(prefix)+"_FROZEN_PRE_XI.npz")
    fmanifest=Path(str(prefix)+"_FROZEN_MANIFEST.json")
    np.savez_compressed(freeze,**payload)
    fsha=sha256_file(freeze)
    fman={
        "version":VERSION,"build":BUILD,
        "parent_v97_sha256":sha,"sha256":fsha,
        "L":args.L,"depths":depths,
        "tail_moment_order":args.tail_moment_order,
        "identity_max":args.identity_max,
        "identity_step":args.identity_step,
        "formula":"Weyl + exact Archimedean + frozen Mangoldt density + half-level midpoint",
        "Xi_used_pre_freeze":False,
        "known_zeros_used_pre_freeze":False,
    }
    fmanifest.write_text(json.dumps(json_safe(fman),indent=2),encoding="utf-8")
    print("[freeze] SHA256                  :",fsha)
    print("[freeze] Xi firewall CLOSED.")

    # POST-FREEZE
    import mpmath as mp
    mp.mp.dps=args.xi_dps
    x0=xi(mp,0)

    print("[3] Post-freeze Xi determinant grid",flush=True)
    xiv=np.asarray([float(xi(mp,float(t))/x0) for t in tgrid],float)
    print("[4] Post-freeze low-order Xi moments",flush=True)
    xim=xi_moments(mp,args.post_moment_order)

    post_rows=[]
    final_d=depths[-1]
    roots=np.asarray(data[f"roots_L{args.L}_d{final_d}"],float)
    raw=finite_det(roots,tgrid)
    partial=partial_moments(roots,args.post_moment_order)

    for name in ["weyl","weyl_midpoint","arch_midpoint","full"]:
        det=completed[final_d][name]
        mom=moments[final_d][name][:args.post_moment_order]
        row={
            "method":name,
            "det_symmetric_rms":sym_rms(det,xiv),
            "moment_relative_l2":rel_l2(mom,xim),
            "S1":float(mom[0]),
            "S1_relative_error":float(abs(mom[0]-xim[0])/abs(xim[0])),
        }
        post_rows.append(row)

    raw_row={
        "method":"raw_prefix",
        "det_symmetric_rms":sym_rms(raw,xiv),
        "moment_relative_l2":rel_l2(partial,xim),
        "S1":float(partial[0]),
        "S1_relative_error":float(abs(partial[0]-xim[0])/abs(xim[0])),
    }
    post_rows.insert(0,raw_row)

    for r in post_rows:
        print(
            f"  {r['method']:<16} det={r['det_symmetric_rms']:.9e} "
            f"mom={r['moment_relative_l2']:.9e} "
            f"S1rel={r['S1_relative_error']:.9e}"
        )

    lookup={r["method"]:r for r in post_rows}
    full=lookup["full"]
    weyl=lookup["weyl"]
    gates={
        "full_S1_below_1e5":full["S1_relative_error"]<1e-5,
        "full_moments_beat_weyl":
            full["moment_relative_l2"]<weyl["moment_relative_l2"],
        "full_determinant_beats_weyl":
            full["det_symmetric_rms"]<weyl["det_symmetric_rms"],
        "full_depth_moments_stable":
            conv[-1]["moment_rel_l2"]<1e-3,
        "full_depth_determinant_stable":
            conv[-1]["det_rel_l2"]<1e-3,
    }

    if all(gates.values()):
        verdict=(
            "PASS_V100_ARITHMETIC_MIDPOINT_TAIL_REPAIRS_MULTIMOMENT_FREDHOLM__"
            "NEXT_PROVE_EULER_MACLAURIN_ARITHMETIC_TAIL_BOUND"
        )
    elif full["S1_relative_error"]<1e-5:
        verdict=(
            "PARTIAL_V100_S1_MECHANISM_CONFIRMED__"
            "HIGHER_MOMENT_OR_DETERMINANT_CLOSURE_NOT_YET_COMPLETE"
        )
    else:
        verdict=(
            "NO_GO_V100_V99_S1_CANCELLATION_DOES_NOT_GENERALIZE"
        )

    pre_csv=Path(str(prefix)+"_PRE_XI_COMPONENTS.csv")
    conv_csv=Path(str(prefix)+"_DEPTH_CONVERGENCE.csv")
    post_csv=Path(str(prefix)+"_POSTFREEZE.csv")
    verdict_path=Path(str(prefix)+"_VERDICT.json")
    write_csv(pre_csv,pre_rows); write_csv(conv_csv,conv); write_csv(post_csv,post_rows)

    packet={
        "version":VERSION,"build":BUILD,"verdict":verdict,
        "RH_proof":False,"parent_v97_sha256":sha,"freeze_sha256":fsha,
        "gates":gates,"postfreeze":post_rows,"depth_convergence":conv,
        "next_theorem":(
            "derive a rigorous Euler-Maclaurin/explicit-formula bound for the "
            "difference between the actual discrete arithmetic tail and the "
            "frozen continuum+midpoint tail, uniformly on compact t-windows"
        ),
    }
    verdict_path.write_text(json.dumps(json_safe(packet),indent=2),encoding="utf-8")

    print("\n"+"="*156)
    print("FINAL VERDICT")
    print("="*156)
    print("verdict                         :",verdict)
    print("raw determinant RMS             :",f"{lookup['raw_prefix']['det_symmetric_rms']:.9e}")
    print("Weyl determinant RMS            :",f"{weyl['det_symmetric_rms']:.9e}")
    print("FULL determinant RMS            :",f"{full['det_symmetric_rms']:.9e}")
    print("Weyl moment relative L2         :",f"{weyl['moment_relative_l2']:.9e}")
    print("FULL moment relative L2         :",f"{full['moment_relative_l2']:.9e}")
    print("FULL S1 relative error          :",f"{full['S1_relative_error']:.9e}")
    print("gates                           :",gates)
    print("Xi used in tail construction    : FALSE")
    print("known zeros used                : FALSE")
    print("RH proof                        : FALSE")
    print("verdict file                    :",verdict_path)
    print("="*156)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
