"""Leading large-charge crossing and resonant-family parameter, U(2).

Conventions: P=2J+Q, eta=2*pi*i*alpha, gamma=eta-i*pi,
lambda=gamma/omega**2.  A partial charge derivative is taken before imposing
gamma=lambda*omega**2.  h below is the real counting value of eta, not gamma.
"""
from __future__ import annotations
from mpmath import mp


def constants():
    W=-mp.pi*mp.root(4,3)/3*(mp.sqrt(3)+1j)
    z=mp.exp(2j*mp.pi/3)
    q0=4/W**2*(2*mp.polylog(2,z)-mp.polylog(2,1/z)-mp.pi**2/6)
    return W,q0,mp.re(-mp.mpf(3)/2*W)


def counting_C(h):
    return 4*(mp.zeta(3)-mp.polylog(3,-mp.exp(h))
              +3*mp.polylog(3,mp.exp(2*h/3))
              -3*mp.polylog(3,-mp.exp(h/3)))


def counting_Cp(h):
    return 4*(-mp.polylog(2,-mp.exp(h))
              +2*mp.polylog(2,mp.exp(2*h/3))
              -mp.polylog(2,-mp.exp(h/3)))


def counting_D(h):
    """Coefficient of 1/omega in log Z_BPS at coincident gauge eigenvalues."""
    return 4*(mp.polylog(2,-mp.exp(h))+mp.pi**2/6
              -mp.polylog(2,-mp.exp(h/3))-mp.polylog(2,mp.exp(2*h/3)))


def counting_charge_entropy(h):
    C=counting_C(h)
    c=counting_Cp(h)/(2*C)**(mp.mpf(2)/3)
    s=3*C**(mp.mpf(1)/3)/2**(mp.mpf(2)/3)-h*c
    return c,s


def resonant_c(lam):
    W,q0,_=constants()
    k,v=-mp.im(W),mp.im(q0)
    return mp.re(q0)-v/mp.sqrt(3)/(1-8*k*k*v/mp.sqrt(3)*lam)


def lambda_from_c(c):
    W,q0,_=constants()
    k,v=-mp.im(W),mp.im(q0)
    return (1+v/(mp.sqrt(3)*(c-mp.re(q0))))/(8*k*k*v/mp.sqrt(3))


def resonant_rho(c):
    W,q0,_=constants()
    k,v=-mp.im(W),mp.im(q0)
    t=c-mp.re(q0)
    return k*(t*t+v*v)/(2*mp.pi*v)*(mp.sqrt(3)-v/t)


def asymptotic_match(dps=80):
    with mp.workdps(dps):
        W,q0,s_index=constants()
        h=mp.findroot(lambda h:counting_charge_entropy(h)[1]-s_index,
                      (-7,-5),tol=mp.power(10,-dps+10))
        c,s=counting_charge_entropy(h)
        lam=lambda_from_c(c)
        C=counting_C(h)
        bps_first=-counting_D(h)/(2*C)**(mp.mpf(1)/3)
        index_first=mp.re((16*mp.pi**2/9)/W)
        crossing_first=(index_first-bps_first)/(-h)
        return dict(h=h,c=c,lam=lam,rho=resonant_rho(c),W=W,q0=q0,
                    s_index=s_index,C=C,D=counting_D(h),
                    crossing_first=crossing_first,
                    entropy_equation_residual=abs(s-s_index),
                    resonant_equation_residual=abs(resonant_c(lam)-c))


def index_saddle(P):
    """Exact cubic saddle of the paper's smooth M=1 index free energy.

    This is the + branch; the conjugate branch has the same real entropy.
    Selection is by continuation from W*P^(-1/3), Re omega<0.
    """
    P=mp.mpf(P)
    A=mp.mpf(16)/27
    b=mp.pi*1j
    W,_,_=constants()
    omega=W/(P+A)**(mp.mpf(1)/3)
    for _ in range(80):
        f=(P+A)*omega**3-3*A*b*b*omega-2*A*b**3
        df=3*(P+A)*omega**2-3*A*b*b
        step=f/df
        omega-=step
        if abs(step)<mp.eps*max(1,abs(omega))*10:
            break
    F=-A*(omega+b)**3/omega**2
    S=F-omega*P
    return dict(omega=omega,F=F,S=S,entropy=mp.re(S))
