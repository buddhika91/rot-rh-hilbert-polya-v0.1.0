#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT-RH v99 — Fast arithmetic-tail / midpoint audit for Tr(H^-2)
================================================================

Loads the exact frozen v97 operator prefix.  It does NOT alter the operator.

The v98 Weyl completion used the continuum tail

    int_Ecut^inf E^-2 rho_Weyl(E) dE.

For a discrete half-integer quantization F(E_k)=k-1/2, that continuum integral
starts half a level too early.  The leading Euler/midpoint correction is

    - 1/(2 Ecut^2).

The same frozen arithmetic phase also contributes an oscillatory density

    P'_L(E)
      = -(1/pi) sum_n w_n log(n) cos(E log n),

so its S1 tail contribution is

    int_Ecut^inf P'_L(E) E^-2 dE.

For a=log(n),

    int_Ecut^inf cos(aE)/E^2 dE
      = cos(a Ecut)/Ecut - a [pi/2 - Si(a Ecut)].

All of this is zero-free.

Only after the corrected S1 is frozen do we evaluate the exact Xi S1
coefficient externally.  This is a mechanism audit, not an RH proof.
"""

from __future__ import annotations

import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.special import sici
from scipy.integrate import quad

EPS=1e-15

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def weyl_s1(E):
    E=float(E)
    return (math.log(E/(2*math.pi))+1.0)/(2*math.pi*E)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--v97-prefix",default="rot_rh_full_rot_hp_pipeline_v97")
    ap.add_argument("--base-module",default="rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit")
    ap.add_argument("--L",type=int,default=16000)
    ap.add_argument("--depth",type=int,default=48)
    ap.add_argument("--tail-factor",type=float,default=18.0)
    ap.add_argument("--xi-dps",type=int,default=60)
    ap.add_argument("--prefix",default="rot_rh_tail_midpoint_s1_v99")
    args=ap.parse_args()

    vp=Path(args.v97_prefix).expanduser().resolve()
    archive=Path(str(vp)+"_FROZEN_PRE_XI.npz")
    manifest_path=Path(str(vp)+"_FROZEN_MANIFEST.json")
    man=json.loads(manifest_path.read_text(encoding="utf-8"))
    sha=sha256_file(archive)
    if sha!=man["sha256"]:
        raise RuntimeError("v97 frozen SHA mismatch.")

    import importlib
    base=importlib.import_module(args.base_module)
    data=np.load(archive,allow_pickle=False)
    roots=np.asarray(data[f"roots_L{args.L}_d{args.depth}"],float)
    E=float(roots[-1])
    prefix_s1=float(np.sum(roots**-2))

    # Same exponential arithmetic regulator as v97.
    nmax=max(8,int(math.ceil(args.tail_factor*args.L)))
    numbers,lambdas,_=base.mangoldt_terms(nmax)
    nums=numbers.astype(float)
    logs=np.log(nums)
    phase_w=lambdas/(logs*np.sqrt(nums))*np.exp(-nums/float(args.L))

    # Oscillatory density correction to S1.
    x=logs*E
    Si,_Ci=sici(x)
    I2=np.cos(x)/E - logs*(np.pi/2-Si)
    prime_corr=float(
        np.sum(-(phase_w*logs)/math.pi * I2)
    )

    smooth_weyl=weyl_s1(E)

    # Exact Archimedean-density correction relative to leading Weyl density.
    def density_delta(x):
        exact=base.theta_prime(float(x))/math.pi
        weyl=math.log(float(x)/(2*math.pi))/(2*math.pi)
        return (exact-weyl)/(float(x)**2)

    arch_corr=float(quad(
        density_delta,E,np.inf,epsabs=1e-13,epsrel=1e-11,limit=200
    )[0])

    midpoint=-0.5/(E*E)

    corrected_tail=smooth_weyl+arch_corr+prime_corr+midpoint
    corrected_s1=prefix_s1+corrected_tail

    pre={
        "v97_sha256":sha,
        "L":args.L,
        "depth":args.depth,
        "Ecut":E,
        "prefix_S1":prefix_s1,
        "weyl_tail_S1":smooth_weyl,
        "archimedean_density_correction":arch_corr,
        "prime_density_correction":prime_corr,
        "midpoint_half_level_correction":midpoint,
        "corrected_tail_S1":corrected_tail,
        "corrected_total_S1":corrected_s1,
        "Xi_used":False,
        "known_zeros_used":False,
    }

    prefix=Path(args.prefix).expanduser().resolve()
    freeze=Path(str(prefix)+"_FROZEN_PRE_XI.json")
    freeze.write_text(json.dumps(pre,indent=2),encoding="utf-8")
    fsha=sha256_file(freeze)

    print("="*132)
    print("ROT-RH v99 — FAST ARITHMETIC TAIL + MIDPOINT S1 AUDIT")
    print("="*132)
    print("v97 SHA verified             :",True)
    print("Ecut                         :",f"{E:.9f}")
    print("prefix S1                    :",f"{prefix_s1:.15e}")
    print("Weyl tail                    :",f"{smooth_weyl:.15e}")
    print("Archimedean correction       :",f"{arch_corr:+.15e}")
    print("prime-density correction     :",f"{prime_corr:+.15e}")
    print("midpoint correction          :",f"{midpoint:+.15e}")
    print("corrected tail               :",f"{corrected_tail:.15e}")
    print("corrected total S1           :",f"{corrected_s1:.15e}")
    print("[freeze] SHA256              :",fsha)
    print("[freeze] Xi firewall CLOSED.")

    import mpmath as mp
    mp.mp.dps=args.xi_dps
    def xi(t):
        s=mp.mpf("0.5")+1j*t
        return mp.re(
            mp.mpf("0.5")*s*(s-1)
            *mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
        )
    x0=xi(mp.mpf("0"))
    d2=mp.diff(lambda t: mp.log(xi(t)/x0),mp.mpf("0"),2)
    xi_s1=float(-d2/mp.factorial(2))

    raw_rel=abs(prefix_s1-xi_s1)/abs(xi_s1)
    weyl_total=prefix_s1+smooth_weyl
    weyl_rel=abs(weyl_total-xi_s1)/abs(xi_s1)
    corr_rel=abs(corrected_s1-xi_s1)/abs(xi_s1)

    result={
        **pre,
        "freeze_sha256":fsha,
        "Xi_S1":xi_s1,
        "prefix_relative_error":raw_rel,
        "weyl_completed_relative_error":weyl_rel,
        "arithmetic_midpoint_completed_relative_error":corr_rel,
        "improvement_vs_weyl":weyl_rel/max(corr_rel,EPS),
        "RH_proof":False,
    }
    out=Path(str(prefix)+"_VERDICT.json")
    out.write_text(json.dumps(result,indent=2),encoding="utf-8")

    print("\n"+"="*132)
    print("POST-FREEZE RESULT")
    print("="*132)
    print("Xi S1                        :",f"{xi_s1:.15e}")
    print("raw prefix relative error    :",f"{raw_rel:.12e}")
    print("Weyl-completed rel error     :",f"{weyl_rel:.12e}")
    print("corrected relative error     :",f"{corr_rel:.12e}")
    print("improvement vs Weyl          :",f"{weyl_rel/max(corr_rel,EPS):.3f}x")
    print("RH proof                     : FALSE")
    print("verdict file                 :",out)
    print("="*132)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
