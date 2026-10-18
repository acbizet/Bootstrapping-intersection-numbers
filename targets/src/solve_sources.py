"""Three-complex-equation Newton solver, with independent omega and gamma.

The physical equations are F_u=0, F_omega=P_target, F_gamma=Q_target.
No lambda path is substituted into charge derivatives. At fixed integer n,
coordinates kappa=n^2*omega, a=n^4*gamma, X=n^3*delta resolve the pole.
The integer n labels the initial resonance; it is held fixed in continuation.
"""
from __future__ import annotations
import time
from mpmath import mp
from resonant_product import FullProduct
from theory import constants


def leading_start(P,Q,n):
    P,Q=mp.mpf(P),mp.mpf(Q)
    W,q0,_=constants()
    rho=n*n/P**(mp.mpf(1)/3)
    c=Q/P**(mp.mpf(2)/3);d=c-q0
    kap=rho*W;X=rho*rho/d
    lam=(2*d-2j*mp.pi*rho/W)/(W*W*d*d)
    return [kap,lam*kap*kap,X]


def equations(engine,z,P,Q,jacobian=True,need_log=False,scales=None):
    kap,a,X=z;n=engine.n
    ps,qs=(P,Q) if scales is None else scales
    r=engine(kap,a/n**4,X/n**3,need_log=need_log)
    E=-1j*r['Fu']/n**3
    residual=mp.matrix([E,(r['P']-P)/ps,(r['Q']-Q)/qs])
    if not jacobian:return residual,r
    Ek=-mp.j/n**3*(r['Fwu']/n**2+mp.j*r['Fuu']/n)
    Ea=-1j*r['Fgu']/n**7
    Ex=-r['Fuu']/n**6
    Pk=(r['Fww']/n**2+1j*r['Fwu']/n)/ps
    Pa=r['Fwg']/n**4/ps
    Px=-1j*r['Fwu']/n**3/ps
    Qk=(r['Fwg']/n**2+1j*r['Fgu']/n)/qs
    Qa=r['Fgg']/n**4/qs
    Qx=-1j*r['Fgu']/n**3/qs
    return residual,mp.matrix([[Ek,Ea,Ex],[Pk,Pa,Px],[Qk,Qa,Qx]]),r


def solve(P,Q,n,dps=50,start=None,extra_prefix=0,extra_terms=0,
          maxsteps=35,progress=None):
    begin=time.perf_counter()
    with mp.workdps(dps):
        P,Q=mp.mpf(P),mp.mpf(Q)
        if P<=0 or Q<=0:raise ValueError('Positive target charges required')
        z=leading_start(P,Q,n) if start is None else [mp.mpc(x) for x in start]
        engine=FullProduct(n,dps,extra_prefix=extra_prefix,extra_terms=extra_terms)
        tolerance=mp.power(10,-dps+12)
        history=[]
        for it in range(maxsteps):
            residual,J,r=equations(engine,z,P,Q)
            norm=max(abs(x) for x in residual)
            if progress:progress(it,norm,z)
            history.append(mp.nstr(norm,15))
            if norm<tolerance:break
            step=mp.lu_solve(J,residual)
            damping=mp.mpf(1)
            for backtrack in range(16):
                trial=[v-damping*s for v,s in zip(z,step)]
                if mp.re(trial[0])<0:
                    try:
                        rr,_=equations(engine,trial,P,Q,jacobian=False)
                        newnorm=max(abs(x) for x in rr)
                    except (ValueError,ZeroDivisionError,OverflowError):
                        newnorm=mp.inf
                    if newnorm<norm*(1-mp.mpf('.0001')*damping) or newnorm<tolerance:
                        z=trial;break
                damping/=2
            else:raise ArithmeticError('Newton line search stalled; residual '+mp.nstr(norm,12))
        else:raise ArithmeticError('Newton exceeded maxsteps; residual '+mp.nstr(norm,12))
        residual,J,r=equations(engine,z,P,Q,need_log=True)
        norm=max(abs(x) for x in residual)
        kap,a,X=z;om=kap/n**2;gamma=a/n**4
        lam=gamma/om**2
        S=r['F']-om*P-(1j*mp.pi+gamma)*Q
        r.update(n=n,coordinates=z,kappa=kap,a=a,X=X,gamma=gamma,
                 lambda_value=lam,target_P=P,target_Q=Q,
                 c=Q/P**(mp.mpf(2)/3),rho=n*n/P**(mp.mpf(1)/3),
                 S=S,residual=norm,iterations=it,calls=engine.calls,
                 dps=dps,seconds=time.perf_counter()-begin,history=history,
                 jacobian_det=mp.det(J),condition=mp.cond(J))
        return r


def serialize(r):
    fields=['target_P','target_Q','P','Q','omega','gamma','u','lambda_value',
            'kappa','a','X','c','rho','S','residual','jacobian_det','condition']
    out={}
    for k in fields:
        v=r[k]
        if isinstance(v,mp.mpc):
            out[k+'_re']=mp.nstr(mp.re(v),r['dps']+5)
            out[k+'_im']=mp.nstr(mp.im(v),r['dps']+5)
        else:out[k]=mp.nstr(v,r['dps']+5)
    for k in ['n','iterations','calls','dps','seconds','prefix','em_terms','history']:out[k]=r[k]
    return out
