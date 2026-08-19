"""Post-freeze Xi utilities.

Do not import this module from construction code. It is intentionally isolated so the
research pipeline can enforce a target-data firewall.
"""
from __future__ import annotations
import mpmath as mp


def xi(t):
    t = mp.mpf(str(t))
    s = mp.mpf("0.5") + 1j*t
    return mp.re(mp.mpf("0.5")*s*(s-1)*mp.power(mp.pi, -s/2)*mp.gamma(s/2)*mp.zeta(s))


def normalized_xi(t):
    return xi(t)/xi(0)
