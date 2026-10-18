"""Vectorized exact oscillator prefixes for U(2) BPS contour tracking.

Independent source coordinates are z=(omega,gamma,u), eta=i*pi+gamma.
All derivatives hold the other two source coordinates fixed.  No Lambda
regulator and no Euler--Maclaurin approximation is used.  N is the number of
oscillator levels in each tower, so convergence must be checked by doubling N.

F includes the Haar factor (1-e^{iu})(1-e^{-iu}), but not the constant Weyl
factor 1/2 or the source-independent coefficient-extraction normalization.
The logarithm is evaluated tower by tower.  Its imaginary part can jump by
2*pi integers; gradient and Hessian are single-valued meromorphic functions.
Track a continuous phase along a path instead of differentiating log branches.
"""
from __future__ import annotations

import numpy as np


def _log_derivatives(s):
    """log(1-exp(s)), its first and second derivatives, without overflow."""
    s = np.asarray(s, dtype=np.complex128)
    logs = np.empty_like(s)
    first = np.empty_like(s)
    second = np.empty_like(s)
    positive = s.real > 0
    a = s[~positive]
    em1 = np.expm1(a)
    e = np.exp(a)
    logs[~positive] = np.log(-em1)
    first[~positive] = e / em1
    second[~positive] = -e / em1**2
    a = s[positive]
    em1 = np.expm1(-a)
    e = np.exp(-a)
    ll = a + np.log(em1)
    # The algebraically factored logarithm can leave the principal sheet.
    logs[positive] = ll.real + 1j*np.arctan2(np.sin(ll.imag), np.cos(ll.imag))
    first[positive] = -1 / em1
    second[positive] = -e / em1**2
    return logs, first, second


class PrefixProduct:
    """Callable returning (F, gradient[3], Hessian[3,3]) as complex128.

    ``evaluate(z, diagnostics=True)`` additionally supplies distances to the
    explicitly retained pole and zero divisors in logarithmic argument space.
    For a gauge pole its actual u-distance is the reported log-argument
    distance (|d s/du|=1).  Charge shifts and omega motion require their own
    source metric.  At gamma=0 the vector numerator and denominator coalesce:
    diagnostics marks these bosonic divisors as removable, not actual poles.
    """
    def __init__(self, levels=1600):
        self.levels = int(levels)
        if self.levels < 1:
            raise ValueError('levels must be positive')
        k = np.arange(self.levels, dtype=float)
        ww, qq, uu, phases, coefficients, labels, indices = [], [], [], [], [], [], []
        # c,d,constant phase,charge,multiplicity per k.
        towers = ((1., 2., 0j, 1., 1., 'vector_numerator'),
                  (1., 0., 0j, 0., -1., 'vector_denominator'),
                  (4/3, 1., -2j*np.pi/3, 1/3, 3., 'chiral_numerator'),
                  (2/3, 1., 2j*np.pi/3, 2/3, -3., 'chiral_denominator'))
        for sign, mult in ((0, 2), (1, 1), (-1, 1)):
            for c,d,phase,charge,coefficient,label in towers:
                ww.append(k+c)
                qq.append(np.full_like(k, charge))
                uu.append(np.full(k.shape, 1j*sign, dtype=complex))
                phases.append(np.full(k.shape, phase, dtype=complex))
                coefficients.append(mult*coefficient*(k+d))
                labels.extend([(label, sign)]*len(k))
                indices.append(k.astype(int))
        self.directions = np.stack((np.concatenate(ww),np.concatenate(qq),np.concatenate(uu)))
        self.phase = np.concatenate(phases)
        self.coefficients = np.concatenate(coefficients)
        self.labels = labels
        self.indices = np.concatenate(indices)
        # Cached pairwise direction products avoid repeated outer products.
        self.pairs = np.stack([self.directions[i]*self.directions[j]
                              for i,j in ((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))])

    def __call__(self, z):
        return self.evaluate(z)

    def evaluate(self, z, diagnostics=False):
        z = np.asarray(z,dtype=complex)
        if z.shape != (3,):
            raise ValueError('z must have shape (3,): omega,gamma,u')
        if z[0].real >= 0:
            raise ValueError('Re(omega) must be negative for the oscillator product')
        # Explicit expressions avoid BLAS startup overhead for tiny contractions.
        s = (self.directions[0]*z[0] + self.directions[1]*z[1]
             + self.directions[2]*z[2] + self.phase)
        logs,K,D = _log_derivatives(s)
        a = self.coefficients
        F = np.sum(a*logs)
        grad = np.sum(self.directions*(a*K)[None,:],axis=1)
        hh = np.sum(self.pairs*(a*D)[None,:],axis=1)
        H = np.array([[hh[0],hh[1],hh[2]],
                      [hh[1],hh[3],hh[4]],
                      [hh[2],hh[4],hh[5]]],dtype=complex)
        # Haar derivatives from the same stable logarithm implementation.
        hlogs,hK,hD = _log_derivatives(np.array([1j*z[2],-1j*z[2]]))
        F += np.sum(hlogs)
        grad[2] += 1j*(hK[0]-hK[1])
        H[2,2] -= np.sum(hD)
        if not diagnostics:
            return F,grad,H
        # Distance to the entire 2*pi*i singular lattice, not just s=0.
        reduced = s.real + 1j*((s.imag+np.pi)%(2*np.pi)-np.pi)
        distances = abs(reduced)
        records = {}
        for kind, mask in (('pole',a<0),('zero',a>0)):
            ids = np.flatnonzero(mask)
            order = ids[np.argsort(distances[ids])[:8]]
            records[kind] = [dict(distance=float(distances[j]),
                                  tower=self.labels[j][0], sign=self.labels[j][1],
                                  k=int(self.indices[j]), power=float(a[j]),
                                  removable_at_exact_index=bool(
                                      z[1] == 0 and self.labels[j][0]=='vector_denominator'))
                             for j in order]
        records['haar_zero_distance'] = float(min(abs(z[2]-2*np.pi*np.rint(z[2].real/(2*np.pi))),
                                                  abs(z[2]+2*np.pi*np.rint(-z[2].real/(2*np.pi)))))
        records['last_level_max_Re_s'] = float(max(s[self.indices==self.levels-1].real))
        records['levels'] = self.levels
        return F,grad,H,records


def microcanonical(engine,z,P,Q):
    """Action and derivatives including eta=i*pi+gamma Legendre term."""
    F,g,H = engine(z)
    z=np.asarray(z,dtype=complex)
    return F-z[0]*P-(1j*np.pi+z[1])*Q, g-np.array([P,Q,0]), H
