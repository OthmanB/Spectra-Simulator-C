#!/usr/bin/env python
import numpy as np
from scipy import special, integrate
import itertools


def Ylm2(_theta, _phi, _l, _m):
    _s = special.sph_harm(_m, _l, _phi, _theta)
    return (_s.real**2 + _s.imag**2) * np.sin(_theta)

lmax = 3
lm_set = [[(l, m) for m in range(-l, l+1)] for l in range(0, lmax + 1)]
lm_set = list(itertools.chain.from_iterable(lm_set))

phi_range = [0, 2.*np.pi]
theta_range = [0, np.pi/4.]

results = [
    (l, m,
     integrate.dblquad(Ylm2,
        phi_range[0], phi_range[1], theta_range[0], theta_range[1], args=(l, m,)))
    for l, m in lm_set
    ]

l=3;
print("Ylm2(", l, "):")
theta=90.;
phi=90./2;
for m in range(-l, l+1):
	r=Ylm2(theta, phi, l, m)
	print(r)
print("------")
print("Integrals:")
for r in results:
    print(r)
