"""Arbitrary-precision full-product saddle equations.

mpmath is used for every source, sum, derivative, and Newton solve. The first
oscillator levels are summed exactly, including the resonant level. Only the
smooth, infinite remainder is accelerated by Euler--Maclaurin summation. All
partial source derivatives are formed before imposing gamma=lambda*omega**2.
Use check_resonant_saddles.py for the fitted-lambda grid and reproducible
precision/tail checks. This module is adapted from the preceding numerical
note in this conversation and is included here as a standalone dependency.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from functools import lru_cache
from pathlib import Path

from mpmath import mp

HERE = Path(__file__).resolve().parent
LAMBDAS = [('-.25', '0'), ('-.5', '0'), ('-1', '0'), ('.10', '0'),
           ('.12', '0'), ('-.5', '.1'), ('-.5', '-.1')]


def constants():
    W = -mp.pi * mp.root(4, 3) / 3 * (mp.sqrt(3) + 1j)
    z = mp.exp(2j * mp.pi / 3)
    q0 = 4 * (2 * mp.polylog(2, z) - mp.polylog(2, 1/z)
              - mp.pi**2/6) / W**2
    return W, q0, mp.re(-mp.mpf('1.5') * W)


@lru_cache(None)
def eulerian_rows(max_order):
    """Integer coefficients of Li_{-m}(z)=z*A_m(z)/(1-z)**(m+1)."""
    rows = [(), (1,)]
    for m in range(2, max_order + 1):
        old = rows[-1]
        rows.append(tuple((j + 1) * (old[j] if j < len(old) else 0)
                          + (m - j) * (old[j-1] if j else 0)
                          for j in range(m)))
    return tuple(rows)


def boundary_polylogs(s, max_order, need_log):
    """All integer polylogs needed by one EM boundary; no low-precision cast."""
    z = mp.exp(s)
    den = -mp.expm1(s)
    values = {0: z/den, 1: -mp.log(den)}
    values[2] = mp.polylog(2, z)
    values[3] = mp.polylog(3, z)
    # Extra guard digits in FullProduct protect Horner cancellation at complex z.
    rows = eulerian_rows(max_order)
    denpower = den
    for m in range(1, max_order + 1):
        denpower *= den
        p = mp.mpf(rows[m][-1])
        for a in reversed(rows[m][:-1]):
            p = p*z + a
        values[-m] = z*p/denpower
    return values


def polynomial_derivatives(N, c, d, j):
    """Derivatives at N of p(k)=(k+d)*(k+c)**j, degree j+1."""
    result = []
    for h in range(j + 2):
        v = mp.mpf(0)
        if h <= j:
            v += (N+d) * (math.factorial(j)//math.factorial(j-h)) * (N+c)**(j-h)
        if 1 <= h <= j+1:
            v += h * (math.factorial(j)//math.factorial(j-h+1)) * (N+c)**(j-h+1)
        result.append(v)
    return result


class FullProduct:
    def __init__(self, n, dps=80, extra_prefix=0, extra_terms=0):
        self.n = int(n)
        self.dps = int(dps)
        self.N = max(128, self.n + 2*self.dps) + int(extra_prefix)
        self.terms = (self.dps + 11)//2 + int(extra_terms)
        self.calls = 0
        with mp.workdps(self.dps + 35):
            self.em = [mp.bernoulli(2*j)/mp.factorial(2*j)
                       for j in range(1, self.terms + 1)]

    def tower(self, om, gamma, delta, sign, c, d, phase, charge, need_log):
        """L, L_v, L_w, L_vv, L_wv, L_ww for one unmultiplied tower."""
        n, N = self.n, self.N
        shift = sign*delta + phase + charge*gamma
        vals = [mp.mpc(0) for _ in range(6)]
        # Form the singular argument as (k+c-sign*n)*omega+sign*delta+...
        # so k=n is represented without subtracting nearby computed numbers.
        z = mp.exp((c-sign*n)*om + shift)
        q = mp.exp(om)
        for k in range(N):
            s = (k+c-sign*n)*om + shift
            # Guard digits make the exponential recurrence accurate; expm1 is
            # retained at small arguments, especially the selected resonance.
            den = mp.expm1(s) if abs(s) < mp.mpf('.01') else z-1
            if abs(s) < mp.mpf('.01'):
                z = den + 1
            K = z/den
            D = -z/(den*den)
            w, y = k+d, k+c
            if need_log:
                vals[0] += w*mp.log(-den)
            vals[1] += w*K
            vals[2] += w*y*K
            vals[3] += w*D
            vals[4] += w*y*D
            vals[5] += w*y*y*D
            z *= q
        sN = (N+c-sign*n)*om + shift
        L = boundary_polylogs(sN, 2*self.terms, need_log)
        powers = [om**j for j in range(2*self.terms + 1)]
        for index, j, r in ((0, 0, 0), (1, 0, 1), (2, 1, 1),
                            (3, 0, 2), (4, 1, 2), (5, 2, 2)):
            if index == 0 and not need_log:
                continue
            poly = polynomial_derivatives(N, c, d, j)
            # Integral_N^infinity p(k)*d_s^r log(1-exp(s(k))) dk.
            tail = mp.fsum((-1)**h * ph * L[2-r+h]/om**(h+1)
                           for h, ph in enumerate(poly))
            tail -= poly[0] * L[1-r] / 2
            for ell, bern in enumerate(self.em, 1):
                m = 2*ell-1
                tail += bern * mp.fsum(math.comb(m, h)*ph*powers[m-h]
                                       * L[1-r-m+h]
                                       for h, ph in enumerate(poly) if h <= m)
            vals[index] += tail
        return vals

    def __call__(self, kappa, gamma, delta, need_log=False):
        self.calls += 1
        with mp.workdps(self.dps + 35):
            n = self.n
            om = kappa/n**2
            if mp.re(om) >= 0:
                raise ValueError('Re(omega) must be negative')
            u = 1j*n*om - 1j*delta
            F = mp.log(4*mp.sin(u/2)**2) if need_log else mp.mpc(0)
            P = Q = mp.mpc(0)
            Fu = mp.cot(u/2)
            Fuu = -mp.mpf('.5')/mp.sin(u/2)**2
            Fww = Fwg = Fgg = Fwu = Fgu = mp.mpc(0)
            third = mp.mpf(1)/3
            for sign, multiplicity in ((0, 2), (1, 1), (-1, 1)):
                for c, d, phase, charge, coeff in (
                    (mp.mpf(1), mp.mpf(2), mp.mpf(0), mp.mpf(1), 1),
                    (mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0), -1),
                    (4*third, mp.mpf(1), -2j*mp.pi/3, third, 3),
                    (2*third, mp.mpf(1), 2j*mp.pi/3, 2*third, -3)):
                    L, V, O, VV, OV, OO = self.tower(
                        om, gamma, delta, sign, c, d, phase, charge, need_log)
                    factor = multiplicity*coeff
                    F += factor*L
                    P += factor*O
                    Q += factor*charge*V
                    Fu += factor*1j*sign*V
                    Fuu -= factor*sign*sign*VV
                    Fww += factor*OO
                    Fwg += factor*charge*OV
                    Fgg += factor*charge*charge*VV
                    Fwu += factor*1j*sign*OV
                    Fgu += factor*1j*sign*charge*VV
            return dict(F=F, P=P, Q=Q, Fu=Fu, Fuu=Fuu,
                        Fww=Fww, Fwg=Fwg, Fgg=Fgg, Fwu=Fwu, Fgu=Fgu,
                        omega=om, u=u, prefix=self.N, em_terms=self.terms)


def leading_roots(kappa, lam):
    a = lam*kappa*kappa
    A = -2j*mp.pi/kappa
    disc = mp.sqrt(1 + A*a)
    return [(-1+disc)/A, (-1-disc)/A]


def corrected_root(n, kappa, lam, branch):
    a = lam*kappa*kappa
    A = -2j*mp.pi/kappa
    b = 2 + A*a/n
    disc = mp.sqrt(b*b + 4*A*a*(1-mp.mpf(1)/n))
    rr = [(-b+disc)/(2*A), (-b-disc)/(2*A)]
    x0 = leading_roots(kappa, lam)[branch]
    return min(rr, key=lambda x: abs(x-x0))


def leading_charge(lam):
    W, q0, _ = constants()
    if not mp.im(lam):
        lam = mp.re(lam)
        k, v = -mp.im(W), mp.im(q0)
        c = mp.re(q0) - v/mp.sqrt(3)/(1-8*k*k*v/mp.sqrt(3)*lam)
        t = c-mp.re(q0)
        rho = k*(t*t+v*v)/(2*mp.pi*v)*(mp.sqrt(3)-v/t)
    else:
        def f(c, rho):
            z = lam*W**2*(c-q0)**2 - 2*(c-q0) + 2j*mp.pi*rho/W
            return mp.re(z), mp.im(z)
        c, rho = mp.findroot(f, (mp.mpf('.16'), mp.mpf('4.2')),
                             tol=mp.eps*100)
    if rho <= 0:
        raise ValueError('rho must be positive')
    return c, rho, rho*W, rho*rho/(c-q0)


def gauge_solve(n, kappa, lam, branch, dps=80, start=None,
                extra_prefix=0, extra_terms=0):
    with mp.workdps(dps):
        engine = FullProduct(n, dps, extra_prefix, extra_terms)
        gamma = lam*(kappa/n**2)**2
        X = leading_roots(kappa, lam)[branch] if start is None else +start
        tol = mp.power(10, -(dps-10))
        for iteration in range(15):
            result = engine(kappa, gamma, X/n**3)
            E, Ep = -1j*result['Fu']/n**3, -result['Fuu']/n**6
            if abs(E) < tol:
                break
            X -= E/Ep
        else:
            raise RuntimeError('gauge Newton iteration did not converge')
        result.update(X=X, kappa=kappa, residual=abs(E), iterations=iteration,
                      calls=engine.calls, dps=dps)
        return result


def physical_solve(n, lam, dps=80, start=None, extra_prefix=0, extra_terms=0):
    with mp.workdps(dps):
        lead = leading_charge(lam)
        c0, rho0, kap0, X0 = lead
        kap, X = (kap0, X0) if start is None else (start['kappa'], start['X'])
        engine = FullProduct(n, dps, extra_prefix, extra_terms)
        ps, qs = n**6/rho0**3, n**4/rho0**2
        tol = mp.power(10, -(dps-10))
        for iteration in range(15):
            result = engine(kap, lam*(kap/n**2)**2, X/n**3)
            E = -1j*result['Fu']/n**3
            residual = mp.matrix([mp.re(E), mp.im(E), mp.im(result['P'])/ps,
                                  mp.im(result['Q'])/qs])
            norm = max(abs(v) for v in residual)
            if norm < tol:
                break
            # Holomorphic chain rule for kappa, X, followed by real/imag parts.
            gamma_kap = 2*lam*kap/n**4
            Ek = -mp.j/n**3*(result['Fwu']/n**2 + result['Fgu']*gamma_kap
                          + result['Fuu']*1j/n)
            Ex = -result['Fuu']/n**6
            Pk = result['Fww']/n**2 + result['Fwg']*gamma_kap + result['Fwu']*1j/n
            Px = result['Fwu']*(-1j)/n**3
            Qk = result['Fwg']/n**2 + result['Fgg']*gamma_kap + result['Fgu']*1j/n
            Qx = result['Fgu']*(-1j)/n**3
            J = mp.matrix([[mp.re(Ek), -mp.im(Ek), mp.re(Ex), -mp.im(Ex)],
                           [mp.im(Ek), mp.re(Ek), mp.im(Ex), mp.re(Ex)],
                           [mp.im(Pk)/ps, mp.re(Pk)/ps, mp.im(Px)/ps, mp.re(Px)/ps],
                           [mp.im(Qk)/qs, mp.re(Qk)/qs, mp.im(Qx)/qs, mp.re(Qx)/qs]])
            step = mp.lu_solve(J, residual)
            kap -= mp.mpc(step[0], step[1])
            X -= mp.mpc(step[2], step[3])
        else:
            raise RuntimeError('real-charge Newton iteration did not converge')
        P, Q = mp.re(result['P']), mp.re(result['Q'])
        if P <= 0 or Q <= 0:
            raise RuntimeError('nonpositive charges')
        log_result = engine(kap, lam*(kap/n**2)**2, X/n**3, need_log=True)
        S = mp.re(log_result['F'] - result['omega']*result['P']
                  - (1j*mp.pi + lam*(kap/n**2)**2)*result['Q'])
        result.update(kappa=kap, X=X, residual=norm, iterations=iteration,
                      calls=engine.calls, dps=dps, c_pred=c0, rho_pred=rho0,
                      c=Q/P**(mp.mpf(2)/3), rho=n*n/P**(mp.mpf(1)/3),
                      S_re_normalized=S/P**(mp.mpf(2)/3), leading=lead)
        return result


def fmt(x, digits=None):
    return mp.nstr(x, max(85, mp.dps+5) if digits is None else digits,
                   strip_zeros=False)


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        def cell(k, v):
            # Source labels are the short decimal inputs; extra binary-rounding
            # digits in their printed representation are not measured data.
            if k in ('lambda_re', 'lambda_im'):
                return mp.nstr(mp.mpf(v), 25)
            return fmt(v) if isinstance(v, (mp.mpf, mp.mpc)) else v
        writer.writerows({k: cell(k, v) for k, v in row.items()} for row in rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--dps', type=int, default=80)
    args = parser.parse_args()
    mp.dps = args.dps
    t = time.perf_counter()
    n = 24
    kap, lam = mp.mpc('-3', '-1.7'), mp.mpf('-.25')
    result = gauge_solve(n, kap, lam, 0, args.dps)
    print('gauge', fmt(result['X']), 'residual', fmt(result['residual']),
          'seconds', time.perf_counter()-t, flush=True)
    t = time.perf_counter()
    result = physical_solve(48, lam, args.dps)
    print('physical', fmt(result['c']), 'residual', fmt(result['residual']),
          'seconds', time.perf_counter()-t, flush=True)
