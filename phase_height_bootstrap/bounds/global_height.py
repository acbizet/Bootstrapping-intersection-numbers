"""Analytic global height reduction, evaluated at supplied finite targets.

This is not a local Gaussian optimization.  The proof reduces the entire
three-angle torus to a one-dimensional profile before finding its maximum.
Its inequalities are exact; target sources and scalar bounds are numerically
evaluated rather than interval-certified.
"""
from pathlib import Path
import argparse
import json
import math
import time

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent


def towers(w, h, levels):
    k = np.arange(levels, dtype=float)
    y = math.exp((w+h)/3)
    # Each multiplicity is per adjoint weight.  Charged weights are a pair.
    b = [(np.exp((k+1)*w), k),
         (np.exp(k*w)*y*y, 3*(k+1))]
    f = [(np.exp((k+1)*w+h), k+2),
         (np.exp((k+1)*w)*y, 3*(k+1))]
    return b, f


def tail_bound(w, h, levels, s):
    """Upper bound on omitted loss and on omitted s*loss'(s)."""
    a = math.exp(w)
    y = math.exp((w+h)/3)
    total = 0.
    # (r_N, multiplicity factor, offset, bosonic).
    for r0, mult, off, boson in (
        (math.exp((levels+1)*w), 1., 0., True),
        (math.exp(levels*w)*y*y, 3., 1., True),
        (math.exp((levels+1)*w+h), 1., 2., False),
        (math.exp((levels+1)*w)*y, 3., 1., False)):
        weighted = mult*r0*((levels+off)/(1-a)+a/(1-a)**2)
        den = (1-r0)**2 if boson else 1-4*r0*s
        if den <= 0:
            return math.inf
        total += 4*s*weighted/den
    return total


class Profile:
    def __init__(self, w, h, levels):
        b, f = towers(w, h, levels)
        self.a = np.concatenate([4*r/(1-r)**2 for r, m in b])
        self.mb = np.concatenate([m for r,m in b])
        self.b = np.concatenate([4*r/(1+r)**2 for r,m in f])
        self.mf = np.concatenate([m for r,m in f])
        self.w, self.h, self.levels = w, h, levels
        r, m = b[0]
        self.coarse_C = float(np.sum(m*np.arctanh(r)))

    def loss(self, s):
        return float(np.dot(self.mb, np.log1p(self.a*s))
                     - np.dot(self.mf, np.log1p(-self.b*s)))

    def deficit(self, s):
        return self.loss(s)-math.log(s)

    def g(self, s):
        """s*D'(s), strictly increasing."""
        return float(np.dot(self.mb,self.a*s/(1+self.a*s))
                     + np.dot(self.mf,self.b*s/(1-self.b*s))-1)

    def upper_tail(self, s):
        return tail_bound(self.w,self.h,self.levels,s)

    def stationary_bracket(self):
        root = brentq(self.g, 1e-18, 1., xtol=1e-25, rtol=1e-14)
        # This inflation is many orders above double arithmetic error and
        # the infinite tail.  Explicit sign checks are retained in the data.
        lo, hi = root*(1-1e-9), root*(1+1e-9)
        assert self.g(lo)+self.upper_tail(lo) < 0 < self.g(hi)
        return root, lo, hi


def analyze(target, levels=8192):
    w, h = float(target['positive_omega']), float(target['positive_eta'])
    budget = float(target['B0_minus_ReS'])
    p = Profile(w,h,levels)
    # Global confinement holds without first constraining chemical angles.
    outside_small_gauge_loss = 3*p.coarse_C-math.log(.75)
    assert outside_small_gauge_loss > budget
    root, lo, hi = p.stationary_bracket()
    # L is increasing.  Root lies in [lo,hi], hence this is a lower bound on
    # min D, including the omitted nonnegative oscillator levels.
    d_lower = p.loss(lo)-math.log(hi)
    d_upper = p.deficit(root)+p.upper_tail(root)
    assert outside_small_gauge_loss > d_upper
    max_gap_upper = budget-d_lower
    max_gap_lower = budget-d_upper
    allowed = None
    if max_gap_upper >= 0:
        if p.deficit(root) <= budget:
            # Solve the left root in log(s); at large budget it can be much
            # smaller than any sensible absolute tolerance on s itself.
            left=math.exp(brentq(lambda ls:p.deficit(math.exp(ls))-budget,
                                -745.,math.log(root),xtol=1e-12,rtol=1e-14))
            right=brentq(lambda s:p.deficit(s)-budget,root,.75,
                          xtol=1e-25,rtol=1e-13)
            # A conservative inflated enclosure for prefix sublevel roots.
            left*=1-1e-9
            right*=1+1e-9
            assert p.deficit(left)>budget and p.deficit(right)>budget
            allowed=dict(s_min=left,s_max=right,
                         abs_u_min=2*math.asin(math.sqrt(left)),
                         abs_u_max=2*math.asin(math.sqrt(right)),
                         boundary_deficits=[p.deficit(left),p.deficit(right)])
    return dict(P=float(target['P']),Q=float(target['Q']),n=target['n'],
        levels=levels,source_omega=w,source_eta=h,budget_B0_minus_hres=budget,
        coarse_vector_C=p.coarse_C,
        coarse_loss_at_abs_u_2pi_over3=outside_small_gauge_loss,
        global_positive_chemical_phase_reduction=True,
        minimizing_s=root,minimizing_abs_u=2*math.asin(math.sqrt(root)),
        root_bracket=[lo,hi],
        root_signs=[p.g(lo)+p.upper_tail(lo),p.g(hi)],
        maximum_deficit_lower=d_lower,maximum_deficit_upper=d_upper,
        maximum_height_minus_hres_lower=max_gap_lower,
        maximum_height_minus_hres_upper=max_gap_upper,
        infinite_tail_upper=p.upper_tail(root),
        height_excluded=max_gap_upper<0,
        allowed_gauge_enclosure=allowed,
        global_phase_excluded=None,
        interval_certified=False,
        qualification='Exact analytic global inequalities, numerically evaluated at non-interval-certified targets.')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--levels',type=int,default=8192)
    args=ap.parse_args()
    start=time.perf_counter()
    rows=[]
    for path in sorted((HERE.parent/'targets').glob('target_Q*.json')):
        if '_validation_' in path.name: continue
        row=analyze(json.loads(path.read_text()),args.levels)
        rows.append(row)
        print(json.dumps({k:row[k] for k in ('Q','maximum_deficit_lower',
             'maximum_height_minus_hres_upper','height_excluded','allowed_gauge_enclosure')}),flush=True)
    result=dict(method='global exact oscillator-block inequality followed by monotone scalar root',
                seconds=time.perf_counter()-start,rows=rows)
    (HERE/'global_height_results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
