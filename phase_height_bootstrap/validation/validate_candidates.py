"""Independent mpmath checks and local KKT refinement; no global certificate.

Reads only the sibling targets/search outputs. Uses a separate FullProduct
instance with different resonance coordinate parameter n, precision and prefix.
"""
from pathlib import Path
import json
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE/'src'))
from mpmath import mp
from resonant_product import FullProduct


def num(x):
    return mp.mpf(str(x))


def nstr(x):
    return mp.nstr(x, 65)


def cstr(x):
    return dict(re=nstr(mp.re(x)), im=nstr(mp.im(x)))


def phase(x):
    return mp.atan2(mp.sin(x), mp.cos(x))


def exact_prefix_action(w,h,P,Q,angles,levels=96):
    """Independent finite product in physical x,y angles; no EM tail.

    Its purpose is an exact-prefix periodicity identity, not a converged
    approximation to the full height. Raw EM in oscillator k is unreliable
    when omega is shifted by a large imaginary period.
    """
    tx,ty,u=angles
    U=mp.exp(1j*u)
    value=mp.log(4*mp.sin(u/2)**2)
    third=mp.mpf(1)/3
    for k in range(levels):
        for radius,angle,power,fermion in (
            (mp.exp((k+1)*w+h),k*tx+3*ty,k+2,True),
            (mp.exp((k+1)*w),(k+1)*tx,-k,False),
            (mp.exp((k+4*third)*w+h/3),(k+1)*tx+ty,3*(k+1),True),
            (mp.exp((k+2*third)*w+2*h/3),k*tx+2*ty,-3*(k+1),False)):
            z=radius*mp.exp(1j*angle)
            sign=1 if fermion else -1
            value+=power*(2*mp.log(1+sign*z)+mp.log(1+sign*z*U)+mp.log(1+sign*z/U))
    return value-P*(w+1j*tx)-Q*(h+1j*(3*ty-tx))


def evaluate(engine, w, h, P, Q, angles):
    tx,ty,u = angles
    om=w+1j*tx
    eta=h+1j*(3*ty-tx)
    z=engine(engine.n**2*om, eta-1j*mp.pi,
             engine.n*om+1j*u, need_log=True)
    S=z['F']-P*om-Q*eta
    g=mp.matrix([z['P']-P,z['Q']-Q,z['Fu']])
    H=mp.matrix([[z['Fww'],z['Fwg'],z['Fwu']],
                 [z['Fwg'],z['Fgg'],z['Fgu']],
                 [z['Fwu'],z['Fgu'],z['Fuu']]])
    D=mp.matrix([[1j,0,0],[-1j,3j,0],[0,0,1]])
    return S,D.T*g,D.T*H*D


def serialize_eval(S,g,H,ref_height,target_phase):
    return dict(S=cstr(S), height_margin=nstr(mp.re(S)-ref_height),
                phase_mod_2pi=nstr(phase(mp.im(S))),
                saddle_phase_mismatch=nstr(phase(mp.im(S)-target_phase)),
                physical_gradient=[cstr(v) for v in g],
                prefix_note='FullProduct exact prefix plus Euler--Maclaurin tail')


def main():
    mp.dps=90
    started=time.perf_counter()
    rows=[]
    for charge in [5883,5940]:
        tfile=ROOT/'targets'/f'target_Q{charge}_validation_75dps.json'
        target=json.loads(tfile.read_text())
        search=json.loads((ROOT/'search'/f'result_Q{charge}.json').read_text())
        w,h=map(num,[target['positive_omega'],target['positive_eta']])
        P,Q=map(num,[target['P'],target['Q']])
        hs=num(target['S']['re']); phi=num(target['phase_mod_2pi'])
        e1=FullProduct(7,55,extra_prefix=256,extra_terms=4)
        e2=FullProduct(11,75,extra_prefix=512,extra_terms=8)
        fixed=[]
        for label in ['positive_peak','winning_phase_candidate']:
            point=search[label]
            a=mp.matrix([num(v) for v in point['physical_angles']])
            S1,g1,H1=evaluate(e1,w,h,P,Q,a)
            S2,g2,H2=evaluate(e2,w,h,P,Q,a)
            fixed.append(dict(label=label,physical_angles=[nstr(v) for v in a],
                **serialize_eval(S2,g2,H2,hs,phi),
                precision_prefix_changes=dict(
                    S_absolute=nstr(abs(S2-S1)),
                    phase_mod_2pi=nstr(abs(phase(mp.im(S2-S1)))),
                    gradient_max=nstr(max(abs(g2[i]-g1[i]) for i in range(3))),
                    hessian_max=nstr(max(abs(H2[i,j]-H1[i,j]) for i in range(3) for j in range(3)))),
                difference_from_float_action=nstr(abs(S2-mp.mpc(num(point['S']['re']),num(point['S']['im']))))))
        # Refine the selected local phase-constrained extremum in the original
        # search whitening. No new search is done and no global conclusion is drawn.
        W=mp.matrix([[num(v) for v in row] for row in search['whitening']])
        center=mp.matrix([num(v) for v in search['center_physical_angles']])
        point=search['winning_phase_candidate']
        p=mp.matrix([num(v) for v in point['p']])
        mu=num(point['lagrange_multiplier'])
        eval0=evaluate(e1,w,h,P,Q,center+W*p)
        branch=phi+2*mp.pi*mp.nint((mp.im(eval0[0])-phi)/(2*mp.pi))
        trace=[]
        for it in range(9):
            angles=center+W*p
            S,g,H=evaluate(e1,w,h,P,Q,angles)
            gp=W.T*g; Hp=W.T*H*W
            gh=mp.matrix([mp.re(v) for v in gp]); gi=mp.matrix([mp.im(v) for v in gp])
            f=mp.matrix([gh[i]-mu*gi[i] for i in range(3)]+[mp.im(S)-branch])
            residual=max(abs(v) for v in f)
            trace.append(dict(iteration=it,residual=nstr(residual)))
            if residual<mp.mpf('1e-44'):
                break
            J=mp.matrix(4)
            for i in range(3):
                for j in range(3):
                    J[i,j]=mp.re(Hp[i,j])-mu*mp.im(Hp[i,j])
                J[i,3]=-gi[i]
                J[3,i]=gi[i]
            step=mp.lu_solve(J,f)
            for i in range(3): p[i]-=step[i]
            mu-=step[3]
        else:
            raise RuntimeError(f'KKT refinement did not converge for Q={charge}')
        angles=center+W*p
        S,g,H=evaluate(e1,w,h,P,Q,angles)
        S2,g2,H2=evaluate(e2,w,h,P,Q,angles)
        gp=W.T*g2; Hp=W.T*H2*W
        gh=mp.matrix([mp.re(v) for v in gp]); gi=mp.matrix([mp.im(v) for v in gp])
        residual=max([abs(gh[i]-mu*gi[i]) for i in range(3)]+[abs(mp.im(S2)-branch)])
        L=mp.matrix([[mp.re(Hp[i,j])-mu*mp.im(Hp[i,j]) for j in range(3)] for i in range(3)])
        normal=gi/mp.sqrt(mp.fdot(gi,gi))
        axis=min(range(3),key=lambda i:abs(normal[i]))
        unit=mp.matrix([int(i==axis) for i in range(3)])
        t1=unit-normal*mp.fdot(unit,normal); t1/=mp.sqrt(mp.fdot(t1,t1))
        t2=mp.matrix([normal[1]*t1[2]-normal[2]*t1[1],
                      normal[2]*t1[0]-normal[0]*t1[2],
                      normal[0]*t1[1]-normal[1]*t1[0]])
        T=mp.matrix([[t1[i],t2[i]] for i in range(3)])
        eig=list(mp.eigsy(T.T*L*T,eigvals_only=True))
        # Exact gauge inversion and chemical conjugation supply distinct checks
        # of all source-coordinate signs and logarithm phase handling.
        inv=mp.matrix([angles[0],angles[1],-angles[2]])
        conj=mp.matrix([-angles[0],-angles[1],angles[2]])
        Si,_,_=evaluate(e2,w,h,P,Q,inv)
        Sc,_,_=evaluate(e2,w,h,P,Q,conj)
        periodic=mp.matrix([angles[0]+2*mp.pi,angles[1],angles[2]])
        # Check periodicity in a separate exact finite product. An EM tail
        # expanded in raw omega cannot be used at omega+2pi*i without first
        # reducing to the correct physical branch.
        E0=exact_prefix_action(w,h,P,Q,angles)
        Ep=exact_prefix_action(w,h,P,Q,periodic)
        row=dict(Q=charge,P=int(P),target_phase=nstr(phi),
            base=dict(dps=55,n=7,prefix=e1.N,em_terms=e1.terms),
            expanded=dict(dps=75,n=11,prefix=e2.N,em_terms=e2.terms),
            fixed_points=fixed,
            refined_local_phase_extremum=dict(physical_angles=[nstr(v) for v in angles],
                whitened_coordinates=[nstr(v) for v in p],
                lagrange_multiplier=nstr(mu),phase_branch=nstr(branch),
                **serialize_eval(S2,g2,H2,hs,phi),
                KKT_residual=nstr(residual),
                constrained_tangent_curvatures=[nstr(v) for v in eig],
                action_precision_change=nstr(abs(S2-S)),trace=trace),
            symmetry_checks=dict(gauge_inversion_height=nstr(abs(mp.re(Si-S2))),
                gauge_inversion_phase=nstr(abs(phase(mp.im(Si-S2)))),
                chemical_conjugation_height=nstr(abs(mp.re(Sc-S2))),
                chemical_conjugation_phase=nstr(abs(phase(mp.im(Sc+S2)))),
                theta_x_2pi_height=nstr(abs(mp.re(Ep-E0))),
                theta_x_2pi_phase=nstr(abs(phase(mp.im(Ep-E0)))),
                periodicity_engine='Independent exact96-level prefix; identity check only, not full-height convergence',
                EM_domain_note='EM evaluations use the small principal omega branch; do not shift raw omega by2pi*i'),
            global_maximum_certified=False,interval_certified=False)
        rows.append(row)
        print(json.dumps(dict(Q=charge,KKT_residual=nstr(residual),
            margin=nstr(mp.re(S2)-hs),precision_change=nstr(abs(S2-S)),
            elapsed=time.perf_counter()-started)),flush=True)
    out=dict(rows=rows,seconds=time.perf_counter()-started,
             scope='Independent mpmath point checks and refinement of known local extrema; no new search and no global upper bound.')
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
