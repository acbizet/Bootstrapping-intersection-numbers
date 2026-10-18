"""Independent mpmath evaluation of the analytic scalar brackets."""
from pathlib import Path
import json
import argparse
from mpmath import mp

HERE=Path(__file__).resolve().parent


def sums(w,h,s,n):
    a=mp.exp(w); y=mp.exp((w+h)/3)
    loss=mp.mpf(0); gg=-mp.mpf(1); tail=mp.mpf(0)
    specs=[(a,1,0,True),(y*y,3,1,True),
           (a*mp.exp(h),1,2,False),(a*y,3,1,False)]
    for r0,mult,off,boson in specs:
        r=r0
        for k in range(n):
            m=mult*(k+off)
            if boson:
                z=4*r*s/(1-r)**2
                loss+=m*mp.log1p(z)
                gg+=m*z/(1+z)
            else:
                z=4*r*s/(1+r)**2
                loss-=m*mp.log1p(-z)
                gg+=m*z/(1-z)
            r*=a
        weighted=mult*r*((n+off)/(1-a)+a/(1-a)**2)
        den=(1-r)**2 if boson else 1-4*r*s
        tail+=4*s*weighted/den
    return loss,gg,tail


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dps',type=int,default=60)
    ap.add_argument('--levels',type=int,default=4096);args=ap.parse_args()
    mp.dps=args.dps
    base=json.loads((HERE/'global_height_results.json').read_text())
    rows=[]
    for row in base['rows']:
        q=int(row['Q']);t=json.loads((HERE.parent/'targets'/f'target_Q{q}.json').read_text())
        w,h=mp.mpf(t['positive_omega']),mp.mpf(t['positive_eta'])
        lo,hi=[mp.mpf(str(x)) for x in row['root_bracket']]
        mid=mp.mpf(str(row['minimizing_s']))
        ll,gl,tl=sums(w,h,lo,args.levels)
        lh,gh,th=sums(w,h,hi,args.levels)
        lm,gm,tm=sums(w,h,mid,args.levels)
        assert gl+tl<0<gh
        lower=ll-mp.log(hi); upper=lm-mp.log(mid)+tm
        budget=mp.mpf(t['B0_minus_ReS'])
        r=dict(Q=q,dps=args.dps,levels=args.levels,
            root_signs=[mp.nstr(gl+tl,45),mp.nstr(gh,45)],
            minimum_deficit_lower=mp.nstr(lower,45),minimum_deficit_upper=mp.nstr(upper,45),
            maximum_height_minus_hres_upper=mp.nstr(budget-lower,45),
            maximum_height_minus_hres_lower=mp.nstr(budget-upper,45),
            infinite_tail_upper=mp.nstr(tm,30),
            agreement_of_lower_bound_with_numpy=mp.nstr(abs(lower-mp.mpf(str(row['maximum_deficit_lower']))),15),
            height_excluded=bool(lower>budget),interval_certified=False)
        rows.append(r); print(json.dumps(r),flush=True)
    (HERE/'scalar_validation.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
