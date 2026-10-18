"""Positive-torus Cauchy bound versus resonant saddle actions.

The inequality is exact algebraically; evaluations use arbitrary-precision
oscillator prefix plus Euler--Maclaurin tail, not interval-certified bounds.
"""
from pathlib import Path
import argparse, csv, json, sys, time
sys.path.insert(0,str(Path(__file__).resolve().parent))
from mpmath import mp
from resonant_product import FullProduct
from solve_sources import solve,serialize
from theory import counting_C,counting_Cp,counting_charge_entropy,asymptotic_match

OUT=Path(__file__).resolve().parent

class Fock:
    def __init__(self,dps=50,extra_prefix=0,extra_terms=0):
        self.engine=FullProduct(0,dps,extra_prefix=extra_prefix,extra_terms=extra_terms)
        self.dps=dps
    def __call__(self,w,h):
        with mp.workdps(self.dps+35):
            F=P=Q=Hww=Hwh=Hhh=mp.mpc(0)
            third=mp.mpf(1)/3
            for c,d,phase,charge,coef in (
                (mp.mpf(1),mp.mpf(2),mp.pi*1j,mp.mpf(1),1),
                (mp.mpf(1),mp.mpf(0),mp.mpf(0),mp.mpf(0),-1),
                (4*third,mp.mpf(1),mp.pi*1j,third,3),
                (2*third,mp.mpf(1),mp.mpf(0),2*third,-3)):
                L,V,O,VV,OV,OO=self.engine.tower(w,h,0,0,c,d,phase,charge,True)
                m=4*coef
                F+=m*L;P+=m*O;Q+=m*charge*V
                Hww+=m*OO;Hwh+=m*charge*OV;Hhh+=m*charge*charge*VV
            return [mp.re(z) for z in (F,P,Q,Hww,Hwh,Hhh)]

def positive_radius(P,Q,dps=50,extra_prefix=0,extra_terms=0):
    eng=Fock(dps,extra_prefix,extra_terms)
    c=Q/P**(mp.mpf(2)/3)
    h=mp.findroot(lambda h:counting_charge_entropy(h)[0]-c,(-7,-4))
    w=-(2*counting_C(h)/P)**(mp.mpf(1)/3)
    tol=mp.power(10,-dps+10)
    for it in range(35):
        F,p,q,A,B,C=eng(w,h)
        res=max(abs(p/P-1),abs(q/Q-1))
        if res<tol:break
        step=mp.lu_solve(mp.matrix([[A,B],[B,C]]),mp.matrix([p-P,q-Q]))
        fac=mp.mpf(1)
        for bt in range(20):
            ww,hh=w-fac*step[0],h-fac*step[1]
            if ww<0 and ww+hh<0:
                v=eng(ww,hh)
                if max(abs(v[1]/P-1),abs(v[2]/Q-1))<res:
                    w,h=ww,hh;break
            fac/=2
        else:raise RuntimeError('positive radius line search')
    else:raise RuntimeError('positive radius iterations')
    return dict(omega=w,eta=h,Fock=F,B=mp.log(4)+F-w*P-h*Q,
                residual=res,iterations=it,prefix=eng.engine.N,terms=eng.engine.terms)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dps',type=int,default=50)
    ap.add_argument('--quick',action='store_true');args=ap.parse_args()
    mp.dps=args.dps
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'data').mkdir(exist_ok=True)
    match=asymptotic_match(args.dps);cx=match['c']
    grids=[mp.mpf('.08'),mp.mpf('.10'),mp.mpf('.12'),cx,mp.mpf('.14'),mp.mpf('.16'),mp.mpf('.20')]
    if args.quick:grids=[mp.mpf('.10'),cx,mp.mpf('.16')]
    cases=[(mp.mpf(100000),12),(mp.mpf(22078774),30)]
    rows=[];details=[]
    for P,n in cases:
        for target_c in grids:
            Q=mp.nint(target_c*P**(mp.mpf(2)/3))
            start=time.perf_counter()
            pos=positive_radius(P,Q,args.dps)
            r=solve(P,Q,n,dps=args.dps)
            gap=pos['B']-mp.re(r['S']);scale=P**(mp.mpf(2)/3)
            row=dict(P=str(P),Q=str(Q),n=n,c=mp.nstr(Q/scale,30),target_c=mp.nstr(target_c,30),
                     B=mp.nstr(pos['B'],40),Re_S_res=mp.nstr(mp.re(r['S']),40),
                     B_minus_ReS=mp.nstr(gap,40),scaled_gap=mp.nstr(gap/scale,30),
                     guaranteed_excluded=(gap<0),positive_omega=mp.nstr(pos['omega'],35),
                     positive_eta=mp.nstr(pos['eta'],35),positive_residual=mp.nstr(pos['residual'],10),
                     saddle_residual=mp.nstr(r['residual'],10),seconds=time.perf_counter()-start)
            rows.append(row);details.append(dict(bound={k:str(v) for k,v in pos.items()},saddle=serialize(r)))
            print(f"P={mp.nstr(P,10)} Q={Q} c={float(Q/scale):.8f} B-S={float(gap):+.8f}",flush=True)
            with (OUT/'data/bounds.json').open('w') as f:json.dump(dict(rows=rows,details=details),f,indent=2)
    # Selected independent precision/tail checks.
    checks=[]
    for idx in [0,len(grids)+len(grids)//2]:
        row=rows[idx];P,Q=mp.mpf(row['P']),mp.mpf(row['Q'])
        with mp.workdps(args.dps+20):
            old=details[idx]['bound'];eng=Fock(args.dps+20,extra_prefix=96,extra_terms=5)
            val=eng(mp.mpf(old['omega']),mp.mpf(old['eta']))
            B=mp.log(4)+val[0]-mp.mpf(old['omega'])*P-mp.mpf(old['eta'])*Q
            checks.append(dict(index=idx,B_change=mp.nstr(abs(B-mp.mpf(row['B'])),20),
                               dps=args.dps+20,prefix=eng.engine.N,em_terms=eng.engine.terms))
    with (OUT/'data/bounds.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    with (OUT/'data/checks.json').open('w') as f:json.dump(checks,f,indent=2)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    for P,n in cases:
        rr=[r for r in rows if mp.mpf(r['P'])==P]
        ax.plot([float(r['c']) for r in rr],[float(r['scaled_gap']) for r in rr],'-o',label=f'P={float(P):.4g}, n={n}')
    ax.axhline(0,color='black',lw=.8)
    ax.axvline(float(cx),color='#a65d00',ls='--',label=r'$c_\infty$')
    ax.set_xlabel(r'$c=Q/P^{2/3}$');ax.set_ylabel(r'$(B-\mathrm{Re}\,\mathcal{S}_{\mathrm{res}})/P^{2/3}$')
    ax.fill_between([.07,.21],[-1,-1],[0,0],color='#d34a4a',alpha=.12,label='Exclusion region: B < Re S')
    ax.set_xlim(.075,.205);ax.set_ylim(-.38,.42);ax.legend(fontsize=9)
    ax.set_title('Original-torus magnitude bound versus resonant saddle')
    fig.tight_layout();fig.savefig(OUT/'bound_gap.png',dpi=180);fig.savefig(OUT/'bound_gap.pdf');plt.close(fig)
    print(json.dumps(dict(rows=len(rows),checks=checks),indent=2))

if __name__=='__main__':main()
