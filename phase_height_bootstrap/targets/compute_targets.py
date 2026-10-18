"""Standalone arbitrary-precision targets for the phase-height bootstrap."""
from pathlib import Path
import argparse
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'src'))
from mpmath import mp
from solve_sources import solve, serialize
from run_bound import positive_radius, Fock


def real(z):
    return mp.nstr(mp.mpf(z), 85)


def complex_obj(z):
    return {'re': real(mp.re(z)), 'im': real(mp.im(z))}


def target_output(r, pos, fock):
    keys = ('Fww', 'Fwg', 'Fwu', 'Fgg', 'Fgu', 'Fuu')
    a, b, c, d, e, f = [r[k] for k in keys]
    H = [[a,b,c], [b,d,e], [c,e,f]]
    phase = mp.im(r['S'])
    wrapped = mp.atan2(mp.sin(phase), mp.cos(phase))
    return dict(P=real(r['target_P']), Q=real(r['target_Q']), n=r['n'],
        c=real(r['c']), coordinate_order=['omega','gamma','u'],
        z=[complex_obj(r[k]) for k in ('omega','gamma','u')],
        **{k:complex_obj(r[k]) for k in ('omega','gamma','u','F','S','lambda_value')},
        H=[[complex_obj(v) for v in row] for row in H],
        residual=real(r['residual']), dps=r['dps'], prefix=r['prefix'],
        em_terms=r['em_terms'], positive_omega=real(pos['omega']),
        positive_eta=real(pos['eta']), positive_gamma=complex_obj(pos['eta']-1j*mp.pi),
        positive_Fock=real(pos['Fock']), B0=real(pos['B']),
        B0_minus_ReS=real(pos['B']-mp.re(r['S'])),
        positive_residual=real(pos['residual']), phase_mod_2pi=real(wrapped),
        exponential_phase=complex_obj(mp.exp(1j*phase)),
        saddle=serialize(r),
        positive_torus=dict(omega=real(pos['omega']), eta=real(pos['eta']),
            gamma=complex_obj(pos['eta']-1j*mp.pi), gauge_Im_u='0',
            B0=real(pos['B']), residual=real(pos['residual']),
            prefix=pos['prefix'], em_terms=pos['terms'],
            Fock=real(fock[0]), P=real(fock[1]), Q=real(fock[2]),
            H=[[real(fock[3]),real(fock[4])],[real(fock[4]),real(fock[5])]]),
        interval_certified=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', type=int, default=60)
    ap.add_argument('--charges', nargs='+', type=int, default=[5880,5882,5883,5884,5940])
    ap.add_argument('--validate', nargs='*', type=int, default=[5883,5940])
    args = ap.parse_args()
    mp.dps = 90
    start_time = time.perf_counter()
    base = {}
    for q in args.charges:
        t = time.perf_counter()
        P,Q = mp.mpf(10_000_000), mp.mpf(q)
        r = solve(P,Q,26,dps=args.dps)
        pos = positive_radius(P,Q,dps=args.dps)
        fock = Fock(args.dps)(pos['omega'],pos['eta'])
        out = target_output(r,pos,fock)
        out['seconds'] = time.perf_counter()-t
        (HERE/f'target_Q{q}.json').write_text(json.dumps(out,indent=2)+'\n')
        base[q] = (r,pos,out)
        print(json.dumps(dict(Q=q,gap=out['B0_minus_ReS'],phase=out['phase_mod_2pi'],seconds=out['seconds'])),flush=True)
    checks = []
    for q in args.validate:
        if q not in base: continue
        old,pos0,out0 = base[q]
        r = solve(mp.mpf(10_000_000),mp.mpf(q),26,dps=75,
            start=old['coordinates'],extra_prefix=128,extra_terms=6)
        pos = positive_radius(mp.mpf(10_000_000),mp.mpf(q),dps=75,
            extra_prefix=128,extra_terms=6)
        fock = Fock(75,extra_prefix=128,extra_terms=6)(pos['omega'],pos['eta'])
        out = target_output(r,pos,fock)
        (HERE/f'target_Q{q}_validation_75dps.json').write_text(json.dumps(out,indent=2)+'\n')
        phase = mp.im(r['S']-old['S'])
        check = dict(Q=q,base_dps=args.dps,validation_dps=75,
            base_prefix=old['prefix'],validation_prefix=r['prefix'],
            base_em_terms=old['em_terms'],validation_em_terms=r['em_terms'],
            source_absolute_changes={k:real(abs(r[k]-old[k])) for k in ('omega','gamma','u')},
            H_absolute_changes={k:real(abs(r[k]-old[k])) for k in ('Fww','Fwg','Fgg','Fwu','Fgu','Fuu')},
            Re_S_absolute_change=real(abs(mp.re(r['S']-old['S']))),
            S_phase_modulo_change=real(abs(mp.atan2(mp.sin(phase),mp.cos(phase)))),
            B0_absolute_change=real(abs(pos['B']-pos0['B'])),
            base_residual=real(old['residual']),validation_residual=real(r['residual']),
            positive_base_residual=real(pos0['residual']),positive_validation_residual=real(pos['residual']),
            interval_certified=False)
        checks.append(check)
        print(json.dumps(check),flush=True)
    (HERE/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    summary = dict(P=10_000_000,n=26,dps=args.dps,seconds=time.perf_counter()-start_time,
        rows=[{k:v for k,v in out.items() if k in ('Q','c','B0','B0_minus_ReS','phase_mod_2pi','residual','positive_residual')}
              for _,_,out in base.values()])
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__': main()
