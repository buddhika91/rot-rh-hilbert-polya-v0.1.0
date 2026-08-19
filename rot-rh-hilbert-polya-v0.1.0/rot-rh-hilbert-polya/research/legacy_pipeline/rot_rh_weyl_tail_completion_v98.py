#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v98 — Frozen Weyl-Tail Completion of the v97 Fredholm Audit
===================================================================

PURPOSE
-------
v97 compared a finite canonical product built from only N spectral levels with
the full entire function Xi(t)/Xi(0).  That comparison is intentionally
conservative but is NOT a fair approximation to the infinite Fredholm
determinant, because the omitted infinitely many levels contribute a smooth,
non-negligible tail.

v98 DOES NOT CHANGE THE OPERATOR.

It loads the exact v97 frozen archive, verifies its SHA256, and constructs a
zero-free asymptotic tail completion from the already-proved fixed-L counting
law

    E_k ~ 2*pi*k/log(k).

Equivalently, the high-energy Weyl density is

    rho(E) ~ (1/(2*pi)) log(E/(2*pi)).

For a finite prefix ending at E_N, the omitted inverse even moments are
approximated by

    T_m(E_N)
      = integral_{E_N}^infinity E^(-2m) rho(E) dE

      = E_N^(1-2m) / [2*pi*(2m-1)]
        * [ log(E_N/(2*pi)) + 1/(2m-1) ].

This formula uses no Xi, zeta values, or zero ordinates.

The missing determinant tail is then

    log D_tail(t)
       = - sum_{m>=1} T_m(E_N) t^(2m)/m,

valid for |t| < E_N.

Thus

    D_completed,N(t)
       = product_{k<=N}(1-t^2/E_k^2)
         * exp[-sum_m T_m t^(2m)/m].

The tail series is truncated only after a prespecified order.  The script
checks its convergence internally and compares completed results across the
v97 depth ladder BEFORE Xi is evaluated.

Only after the completed objects are frozen and hashed does v98 evaluate Xi
externally.

SCIENTIFIC BOUNDARY
-------------------
A good tail-completed match is much stronger evidence than the raw finite
product, but it is still not the exact determinant theorem.  The remaining
proof obligation is to control the difference between the actual arithmetic
tail and the Weyl tail uniformly and then pass to the true trace-norm
K_infinity limit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

EPS = 1e-15
VERSION = "v98"
BUILD = "2026-08-17-frozen-weyl-tail-completion-v98"


def parse_ints(s: str) -> List[int]:
    return [int(x.strip()) for x in str(s).split(",") if x.strip()]


def json_safe(x):
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k,v in x.items()}
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


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: List[Dict[str,Any]]) -> None:
    if not rows:
        path.write_text("",encoding="utf-8"); return
    fields=[]; seen=set()
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k); seen.add(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def finite_det(roots,tgrid):
    r=np.asarray(roots,float)
    inv2=1.0/(r*r)
    return np.asarray([
        float(np.prod(1.0-(float(t)**2)*inv2))
        for t in np.asarray(tgrid,float)
    ])


def partial_moments(roots,order):
    r=np.asarray(roots,float)
    return np.asarray([
        float(np.sum(r**(-2*m))) for m in range(1,order+1)
    ])


def weyl_tail_moments(Ecut,order):
    E=float(Ecut)
    out=[]
    for m in range(1,int(order)+1):
        q=2*m-1
        val=(
            E**(1-2*m)
            /(2.0*math.pi*q)
            *(math.log(E/(2.0*math.pi))+1.0/q)
        )
        out.append(float(val))
    return np.asarray(out,float)


def completed_det(roots,tgrid,tail_moments):
    raw=finite_det(roots,tgrid)
    out=np.zeros_like(raw)
    for i,t in enumerate(np.asarray(tgrid,float)):
        logtail=0.0
        for m,Tm in enumerate(np.asarray(tail_moments,float),start=1):
            logtail -= Tm*(float(t)**(2*m))/m
        out[i]=raw[i]*math.exp(logtail)
    return raw,out


def sym_rms(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.sqrt(np.mean(((a-b)/(1+np.abs(a)+np.abs(b)))**2)))


def rel_l2(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),EPS))


def xi(mp,t):
    s=mp.mpf("0.5")+1j*mp.mpf(str(t))
    return mp.re(
        mp.mpf("0.5")*s*(s-1)
        *mp.power(mp.pi,-s/2)
        *mp.gamma(s/2)
        *mp.zeta(s)
    )


def xi_moments(mp,order):
    x0=xi(mp,0)
    vals=[]
    for m in range(1,order+1):
        d=mp.diff(lambda t: mp.log(xi(mp,t)/x0),mp.mpf("0"),2*m)
        c=d/mp.factorial(2*m)
        vals.append(float(-m*c))
    return np.asarray(vals,float)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--v97-prefix",default="rot_rh_full_rot_hp_pipeline_v97")
    ap.add_argument("--depths",default="24,32,48")
    ap.add_argument("--L-final",type=int,default=16000)
    ap.add_argument("--tail-moment-order",type=int,default=16)
    ap.add_argument("--tail-order-check",type=int,default=12)
    ap.add_argument("--identity-max",type=float,default=50.0)
    ap.add_argument("--identity-step",type=float,default=0.25)
    ap.add_argument("--xi-dps",type=int,default=100)
    ap.add_argument("--prefix",default="rot_rh_weyl_tail_completion_v98")
    args=ap.parse_args()

    v97=Path(args.v97_prefix).expanduser().resolve()
    archive=Path(str(v97)+"_FROZEN_PRE_XI.npz")
    manifest_path=Path(str(v97)+"_FROZEN_MANIFEST.json")
    if not archive.exists() or not manifest_path.exists():
        raise FileNotFoundError("Could not find v97 frozen archive/manifest.")

    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    actual_sha=sha256_file(archive)
    expected_sha=str(manifest["sha256"])
    if actual_sha!=expected_sha:
        raise RuntimeError(
            f"v97 SHA mismatch: expected {expected_sha}, got {actual_sha}"
        )

    depths=parse_ints(args.depths)
    tgrid=np.arange(
        0.0,args.identity_max+0.5*args.identity_step,args.identity_step
    )

    data=np.load(archive,allow_pickle=False)

    print("="*156)
    print("ROT-RH v98 — FROZEN WEYL-TAIL COMPLETION OF v97")
    print("="*156)
    print("v97 SHA verified              :",True)
    print("v97 SHA                       :",actual_sha)
    print("final L                       :",args.L_final)
    print("depths                        :",depths)
    print("tail moment order             :",args.tail_moment_order)
    print("Xi/zeta/known zeros pre-freeze: NONE")
    print("="*156)

    rows=[]
    payload={}
    completed_by_depth={}
    moments_by_depth={}

    print("[1] Zero-free tail completion across frozen depth ladder")
    for d in depths:
        key=f"roots_L{args.L_final}_d{d}"
        if key not in data:
            raise KeyError(f"Missing {key} from v97 archive.")
        roots=np.asarray(data[key],float)
        if args.identity_max>=roots[-1]:
            raise ValueError(
                f"identity-max={args.identity_max} must be below E_N={roots[-1]}"
            )

        partial=partial_moments(roots,args.tail_moment_order)
        tail=weyl_tail_moments(roots[-1],args.tail_moment_order)
        complete=partial+tail

        raw,comp=completed_det(roots,tgrid,tail)

        # Internal tail-series truncation test.
        q=min(args.tail_order_check,args.tail_moment_order)
        _,comp_q=completed_det(roots,tgrid,tail[:q])
        tail_series_drift=rel_l2(comp,comp_q)

        completed_by_depth[d]=comp
        moments_by_depth[d]=complete

        payload[f"roots_d{d}"]=roots
        payload[f"partial_moments_d{d}"]=partial
        payload[f"weyl_tail_moments_d{d}"]=tail
        payload[f"completed_moments_d{d}"]=complete
        payload[f"raw_det_d{d}"]=raw
        payload[f"completed_det_d{d}"]=comp

        row={
            "depth":d,
            "Ecut":float(roots[-1]),
            "partial_S1":float(partial[0]),
            "weyl_tail_S1":float(tail[0]),
            "completed_S1":float(complete[0]),
            "tail_fraction_of_completed_S1":float(tail[0]/complete[0]),
            "tail_series_order_check_relative_l2":tail_series_drift,
        }
        rows.append(row)
        print(
            f"  d={d:<3d} Ecut={roots[-1]:.6f} "
            f"S1partial={partial[0]:.9e} "
            f"S1tail={tail[0]:.9e} "
            f"S1complete={complete[0]:.9e} "
            f"tailfrac={tail[0]/complete[0]:.3%} "
            f"seriesdrift={tail_series_drift:.2e}"
        )

    print("[2] Zero-free completed-depth convergence")
    conv_rows=[]
    for d0,d1 in zip(depths[:-1],depths[1:]):
        r={
            "depth_from":d0,
            "depth_to":d1,
            "completed_moment_relative_l2":
                rel_l2(moments_by_depth[d1],moments_by_depth[d0]),
            "completed_determinant_relative_l2":
                rel_l2(completed_by_depth[d1],completed_by_depth[d0]),
        }
        conv_rows.append(r)
        print(
            f"  {d0}->{d1}: moments={r['completed_moment_relative_l2']:.3e} "
            f"det={r['completed_determinant_relative_l2']:.3e}"
        )

    # Freeze completed tail objects BEFORE Xi.
    prefix=Path(args.prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True,exist_ok=True)
    freeze_path=Path(str(prefix)+"_FROZEN_PRE_XI.npz")
    freeze_manifest=Path(str(prefix)+"_FROZEN_MANIFEST.json")
    payload["tgrid"]=tgrid
    payload["depths"]=np.asarray(depths,int)
    np.savez_compressed(freeze_path,**payload)
    sha=sha256_file(freeze_path)
    freeze_packet={
        "version":VERSION,
        "build":BUILD,
        "parent_v97_sha256":actual_sha,
        "sha256":sha,
        "method":"zero-free Weyl density tail completion",
        "weyl_density":"rho(E)=(1/(2*pi))*log(E/(2*pi))",
        "tail_formula":"T_m=Ecut^(1-2m)/(2*pi*(2m-1))*(log(Ecut/(2*pi))+1/(2m-1))",
        "Xi_used_pre_freeze":False,
        "zeta_used_pre_freeze":False,
        "known_zeros_used_pre_freeze":False,
        "parameters":{
            "L_final":args.L_final,
            "depths":depths,
            "tail_moment_order":args.tail_moment_order,
            "identity_max":args.identity_max,
            "identity_step":args.identity_step,
        },
    }
    freeze_manifest.write_text(
        json.dumps(json_safe(freeze_packet),indent=2),
        encoding="utf-8"
    )
    print("[freeze] v98 SHA256             :",sha)
    print("[freeze] Xi firewall CLOSED.")

    # POST-FREEZE
    import mpmath as mp
    mp.mp.dps=args.xi_dps
    x0=xi(mp,0)
    xiv=np.asarray([float(xi(mp,t)/x0) for t in tgrid],float)
    xim=xi_moments(mp,args.tail_moment_order)

    post_rows=[]
    print("[3] Post-freeze Xi comparison")
    for d in depths:
        roots=np.asarray(data[f"roots_L{args.L_final}_d{d}"],float)
        partial=partial_moments(roots,args.tail_moment_order)
        raw=finite_det(roots,tgrid)
        comp=completed_by_depth[d]
        cm=moments_by_depth[d]

        raw_det=sym_rms(raw,xiv)
        comp_det=sym_rms(comp,xiv)
        raw_mom=rel_l2(partial,xim)
        comp_mom=rel_l2(cm,xim)

        r={
            "depth":d,
            "raw_det_symmetric_rms":raw_det,
            "completed_det_symmetric_rms":comp_det,
            "det_improvement_factor":raw_det/max(comp_det,EPS),
            "raw_moment_relative_l2":raw_mom,
            "completed_moment_relative_l2":comp_mom,
            "moment_improvement_factor":raw_mom/max(comp_mom,EPS),
            "S1_operator_completed":float(cm[0]),
            "S1_Xi":float(xim[0]),
            "S1_relative_error":float(abs(cm[0]-xim[0])/abs(xim[0])),
        }
        post_rows.append(r)
        print(
            f"  d={d:<3d} det raw={raw_det:.9e} complete={comp_det:.9e} "
            f"x{r['det_improvement_factor']:.1f} | "
            f"mom raw={raw_mom:.9e} complete={comp_mom:.9e} "
            f"x{r['moment_improvement_factor']:.1f}"
        )

    final=post_rows[-1]
    gates={
        "tail_series_numerically_converged":
            rows[-1]["tail_series_order_check_relative_l2"]<1e-6,
        "completed_depth_moments_stabilize":
            conv_rows[-1]["completed_moment_relative_l2"]<1e-3,
        "completed_depth_determinant_stabilizes":
            conv_rows[-1]["completed_determinant_relative_l2"]<1e-2,
        "tail_completion_improves_determinant":
            final["det_improvement_factor"]>10.0,
        "tail_completion_improves_moments":
            final["moment_improvement_factor"]>10.0,
    }

    if all(gates.values()):
        verdict=(
            "PASS_V98_V97_GLOBAL_MISMATCH_WAS_DOMINATED_BY_FINITE_DEPTH_TAIL__"
            "NEXT_PROVE_ARITHMETIC_TAIL_ERROR_AND_L_TO_INFINITY_TRACE_NORM"
        )
    elif (
        final["det_improvement_factor"]>2.0
        and final["moment_improvement_factor"]>2.0
    ):
        verdict=(
            "PARTIAL_V98_WEYL_TAIL_EXPLAINS_MUCH_OF_V97_MISMATCH__"
            "BUT_DEPTH_OR_TAIL_CONTROL_NOT_YET_STRONG"
        )
    else:
        verdict=(
            "NO_GO_V98_V97_GLOBAL_MISMATCH_IS_NOT_EXPLAINED_BY_WEYL_TAIL_TRUNCATION"
        )

    rows_csv=Path(str(prefix)+"_TAIL_COMPLETION.csv")
    conv_csv=Path(str(prefix)+"_DEPTH_CONVERGENCE.csv")
    post_csv=Path(str(prefix)+"_POSTFREEZE_XI.csv")
    verdict_path=Path(str(prefix)+"_VERDICT.json")
    write_csv(rows_csv,rows)
    write_csv(conv_csv,conv_rows)
    write_csv(post_csv,post_rows)

    packet={
        "version":VERSION,
        "build":BUILD,
        "verdict":verdict,
        "RH_proof":False,
        "parent_v97_sha256":actual_sha,
        "v98_sha256":sha,
        "gates":gates,
        "final_postfreeze":final,
        "theorem_boundary":{
            "proved_here":(
                "nothing about RH; this is a numerical/analytic tail-completion "
                "audit using the fixed-L Weyl law"
            ),
            "next_required_theorem":(
                "bound actual arithmetic determinant tail minus Weyl tail "
                "uniformly on compact t-windows and prove trace-norm L->infinity"
            ),
        },
        "outputs":{
            "freeze":str(freeze_path),
            "manifest":str(freeze_manifest),
            "tail_completion":str(rows_csv),
            "depth_convergence":str(conv_csv),
            "postfreeze":str(post_csv),
            "verdict":str(verdict_path),
        },
    }
    verdict_path.write_text(
        json.dumps(json_safe(packet),indent=2),
        encoding="utf-8"
    )

    print("\n"+"="*156)
    print("FINAL CONSERVATIVE VERDICT")
    print("="*156)
    print("verdict                         :",verdict)
    print("final raw determinant RMS       :",f"{final['raw_det_symmetric_rms']:.9e}")
    print("final completed determinant RMS :",f"{final['completed_det_symmetric_rms']:.9e}")
    print("determinant improvement factor  :",f"{final['det_improvement_factor']:.3f}")
    print("final raw moment relative L2    :",f"{final['raw_moment_relative_l2']:.9e}")
    print("final completed moment rel L2   :",f"{final['completed_moment_relative_l2']:.9e}")
    print("moment improvement factor       :",f"{final['moment_improvement_factor']:.3f}")
    print("completed S1 relative error     :",f"{final['S1_relative_error']:.9e}")
    print("gates                           :",gates)
    print("Xi used in tail construction    : FALSE")
    print("known zeros used                : FALSE")
    print("RH proof                        : FALSE")
    print("v98 freeze SHA256               :",sha)
    print("verdict file                    :",verdict_path)
    print("="*156)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
