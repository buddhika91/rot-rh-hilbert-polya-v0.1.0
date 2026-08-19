from __future__ import annotations
import numpy as np

EPS = 1e-15

def relative_l2(a, b) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    return float(np.linalg.norm(a-b) / max(np.linalg.norm(b), EPS))

def symmetric_rms(a, b) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(((a-b)/(1+np.abs(a)+np.abs(b)))**2)))
