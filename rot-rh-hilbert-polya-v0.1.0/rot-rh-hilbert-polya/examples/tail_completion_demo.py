"""Small, non-RH demo of the determinant-completion API."""
import numpy as np
from rot_rh_hp.tail import finite_det, weyl_tail_moments, completed_det

roots=np.array([14.0, 21.0, 25.0, 30.0])
t=np.linspace(0,10,6)
raw=finite_det(roots,t)
tail=weyl_tail_moments(float(roots[-1]),4)
completed=completed_det(raw,t,tail)
for x,a,b in zip(t,raw,completed):
    print(f"t={x:5.1f} raw={a:+.8e} completed={b:+.8e}")
