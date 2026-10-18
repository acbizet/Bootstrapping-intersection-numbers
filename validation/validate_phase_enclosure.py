"""Precision/prefix check of the global phase enclosure, plus algebra checks.

This remains floating-point validation of exact inequalities, not an interval
certificate. The independently audited derivation is in mathematical_audit.md.
"""
from pathlib import Path
import importlib.util
import json
import time
from mpmath import mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
spec=importlib.util.spec_from_file_location('global_phase_bound',ROOT/'phase_enclosure/global_phase_bound.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fmt(x):
    return mp.nstr(x,70)


def main():
    mp.dps=90
    started=time.perf_counter()
    rows=json.loads((ROOT/'bounds/global_height_results.json').read_text())['rows']
    base={int(r['Q']):r for r in json.loads((ROOT/'phase_enclosure/results.json').read_text())['rows']}
    checks=[]
    for row in rows:
        if row['height_excluded']: continue
        Q=int(row['Q'])
        target=json.loads((ROOT/f'targets/target_Q{Q}.json').read_text())
        out=module.phase_bound(target,row,N=6000)
        x,y,A,U,X,Y=module.angle_box(target,row)
        kappa=4*mp.cos(U/2)**2
        W=1/(1-x)**2
        K=3*kappa*x*y/(1+x*y)**8
        argument=1-A/(K*W)
        assert argument>=0, 'Y localization requires W-A/K nonnegative.'
        aa=(1-x)**2
        errors=[]
        for theta in [X/2,X,mp.mpf('.1'),mp.pi/2,mp.pi]:
            v=1-mp.cos(theta)
            rational=kappa*x*x*((4+2*x)*v+(4*x*x/aa-2)*v*v)/(aa+2*x*v)**2
            z=x*mp.exp(1j*theta)
            independent=kappa*(x*x/(1-x)**2-mp.re(z*z/(1-z)**2))
            errors.append(abs(rational-independent))
        differences={key:fmt(abs(mp.mpf(out[key])-mp.mpf(base[Q][key]))) for key in [
            'linear_phase_bound','cubic_phase_bound','absolute_unwrapped_phase_bound','phase_exclusion_margin']}
        checks.append(dict(Q=Q,expanded=out,changes_from_5000level_70dps=differences,
            Y_bound_required_nonnegative_argument=fmt(argument),
            independent_complex_geometric_xloss_max_difference=fmt(max(errors)),
            interval_certified=False))
        print(json.dumps(dict(Q=Q,phase_bound=out['absolute_unwrapped_phase_bound'],
            margin=out['phase_exclusion_margin'],argument=fmt(argument))),flush=True)
    (HERE/'phase_enclosure_validation.json').write_text(json.dumps(dict(checks=checks,
        seconds=time.perf_counter()-started,interval_certified=False),indent=2)+'\n')


if __name__=='__main__': main()
