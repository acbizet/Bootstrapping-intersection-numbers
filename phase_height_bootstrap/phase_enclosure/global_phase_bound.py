"""Global finite-angle localization and a conservative phase Taylor bound.

Exact inequalities for the infinite product, numerically evaluated at the
saved targets. No interval certificate is claimed. No flow or optimizer is
used to derive the enclosing source-angle box.
"""
from pathlib import Path
from math import comb
import json
from mpmath import mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def tail_moment(x, N, p):
    jmom = [1/(1-x), x/(1-x)**2, x*(1+x)/(1-x)**3,
            x*(1+4*x+x*x)/(1-x)**4,
            x*(1+11*x+11*x*x+x**3)/(1-x)**5]
    return x**N*mp.fsum(comb(p,j)*(N+2)**(p-j)*jmom[j]
                        for j in range(p+1))


def bisect(f, lo, hi, steps=240):
    for _ in range(steps):
        mid=(lo+hi)/2
        if f(mid)>0: hi=mid
        else: lo=mid
    return lo,hi


def angle_box(target, row):
    w=mp.mpf(target['positive_omega']); h=mp.mpf(target['positive_eta'])
    x=mp.exp(w); y=mp.exp((w+h)/3)
    # Outward safety relative to the already conservatively padded scalar
    # enclosures. These paddings are not a substitute for interval arithmetic.
    A=mp.mpf(str(row['maximum_height_minus_hres_upper']))+mp.mpf('1e-7')
    U=mp.mpf(str(row['allowed_gauge_enclosure']['abs_u_max']))*(1+mp.mpf('1e-6'))
    kappa=4*mp.cos(U/2)**2
    aa=(1-x)**2
    def xloss(v):
        # Closed sum: sum_(k>=2) (k-1)x^k(1-cos(k theta)), v=1-cos theta.
        return kappa*x*x*((4+2*x)*v+(4*x*x/aa-2)*v*v)/(aa+2*x*v)**2
    vstar=(1+x/2)*aa
    assert xloss(2)>A and xloss(vstar)>A
    vlo,vhi=bisect(lambda v:xloss(v)-A,mp.mpf(0),vstar)
    X=mp.acos(1-vhi)*(1+mp.mpf('1e-20'))
    # A positive-coefficient chiral-fermion block supplies a first harmonic
    # in theta_y. Its phase origin shifts slightly as theta_x varies.
    W=1/(1-x)**2
    K=3*kappa*x*y/(1+x*y)**8
    argument=1-A/(K*W)
    # The comparison |T|<=W preserves the cosine bound only if the
    # required real part W-A/K is nonnegative. Otherwise retain all phases.
    if argument<=0: Y=mp.pi
    else:
        shift=X+2*mp.atan(x*X/(1-x))
        Y=min(mp.pi,shift+mp.acos(min(mp.mpf(1),argument)))
    return x,y,A,U,X,Y


def phase_bound(target,row,N=5000):
    x,y,A,U,X,Y=angle_box(target,row)
    mu_x=mp.mpf(0);mu_y=mp.mpf(0);third=mp.mpf(0)
    xk=mp.mpf(1)
    for k in range(N):
        for r,m,kx,ky in [(x*xk,k,k+1,0),(y**3*xk,k+2,k,3),
                           (x*y*xk,3*(k+1),k+1,1),(y*y*xk,3*(k+1),k,2)]:
            f=m*r*(1+r)/(1-r)**3
            mu_x+=f*kx;mu_y+=f*ky
            third+=mp.mpf(2)/3*f*(kx*X+ky*Y)**3
        xk*=x
    linear=U*U*(X*mu_x+Y*mu_y)
    rn=x**N;fac=(1+rn)/(1-rn)**3
    moments={p:tail_moment(x,N,p) for p in range(1,5)}
    linear_tail=U*U*fac*(12*X*moments[2]+36*Y*moments[1])
    third_tail=8*fac*(X**3*moments[4]+9*X*X*Y*moments[3]
                     +27*X*Y*Y*moments[2]+27*Y**3*moments[1])
    residual=10*mp.mpf(target['positive_residual'])*(
        (mp.mpf(target['P'])+mp.mpf(target['Q']))*X+3*mp.mpf(target['Q'])*Y)
    total=linear+third+linear_tail+third_tail+residual
    phase=abs(mp.mpf(target['phase_mod_2pi']))
    def fmt(v):return mp.nstr(v,45)
    return dict(Q=int(mp.mpf(target['Q'])),levels=N,dps=mp.dps,
        remaining_height_budget=fmt(A),u_absolute_bound=fmt(U),
        theta_x_absolute_bound=fmt(X),theta_y_absolute_bound=fmt(Y),
        linear_phase_bound=fmt(linear),cubic_phase_bound=fmt(third),
        infinite_tail_bound=fmt(linear_tail+third_tail),
        fock_residual_phase_allowance=fmt(residual),
        absolute_unwrapped_phase_bound=fmt(total),
        saddle_phase_distance=fmt(phase),phase_exclusion_margin=fmt(phase-total),
        phase_excluded=bool(total<phase),interval_certified=False,
        scope='Global angle localization and global Taylor inequality, numerically evaluated inputs and sums.')


def main():
    mp.dps=70
    rows=json.loads((ROOT/'bounds/global_height_results.json').read_text())['rows']
    result=[]
    for row in rows:
        if row['height_excluded']: continue
        Q=int(row['Q']);target=json.loads((ROOT/f'targets/target_Q{Q}.json').read_text())
        r=phase_bound(target,row);result.append(r)
        print(json.dumps(r),flush=True)
    (HERE/'results.json').write_text(json.dumps(dict(rows=result),indent=2)+'\n')


if __name__=='__main__':main()
