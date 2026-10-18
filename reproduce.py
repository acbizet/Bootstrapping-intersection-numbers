"""Recompute global bounds and plots; --full also repeats target/search work."""
from pathlib import Path
import argparse
import ast
import json
import subprocess
import sys
from mpmath import mp
from summarize import build

HERE=Path(__file__).resolve().parent


def run(name,*args):
    subprocess.run([sys.executable,str(HERE/name),*args],check=True)


def verify(summary):
    mp.dps=80
    expected={5880:'excluded_by_global_height_bound',5882:'excluded_by_global_height_bound',
              5883:'excluded_by_global_phase_bound',5884:'inconclusive_local_search_only',
              5940:'necessary_height_and_phase_test_passes'}
    for row in summary['rows']:
        assert row['status']==expected[row['Q']]
        t=json.loads((HERE/f'targets/target_Q{row["Q"]}.json').read_text())
        assert mp.mpf(t['residual'])<mp.mpf('1e-50')
        assert mp.mpf(t['positive_residual'])<mp.mpf('1e-50')
    pv=json.loads((HERE/'validation/phase_enclosure_validation.json').read_text())['checks']
    p=next(c for c in pv if c['Q']==5883)
    assert mp.mpf(p['expanded']['phase_exclusion_margin'])>mp.mpf('.54')
    assert mp.mpf(p['changes_from_5000level_70dps']['absolute_unwrapped_phase_bound'])<mp.mpf('1e-15')
    assert mp.mpf(p['Y_bound_required_nonnegative_argument'])>0
    candidates=json.loads((HERE/'validation/results.json').read_text())['rows']
    v=next(r for r in candidates if r['Q']==5940)['refined_local_phase_extremum']
    assert mp.mpf(v['height_margin'])>318
    assert abs(mp.mpf(v['saddle_phase_mismatch']))<mp.mpf('1e-45')
    for r in candidates:
        c=r['refined_local_phase_extremum']
        assert mp.mpf(c['KKT_residual'])<mp.mpf('1e-40')
        assert all(mp.mpf(x)<0 for x in c['constrained_tangent_curvatures'])
    sources=list(HERE.rglob('*.py'))
    for path in sources:ast.parse(path.read_text(),filename=str(path))
    out=dict(saved_consistency_checks='pass',sources_parsed=len(sources),
        excluded_by_global_height=[5880,5882],excluded_by_global_phase=[5883],
        unresolved=[5884],necessary_condition_witness=[5940],
        interval_certified=False,global_nonzero_intersection_established=False,
        main_phase_margin=p['expanded']['phase_exclusion_margin'])
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--full',action='store_true')
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    if args.full:
        run('targets/compute_targets.py')
        run('search/search_phase_height.py')
    run('bounds/global_height.py')
    run('bounds/validate_scalar.py')
    run('phase_enclosure/global_phase_bound.py')
    run('validation/validate_phase_enclosure.py')
    if args.full:run('validation/validate_candidates.py')
    summary=build();verify(summary)
    if not args.no_plots:
        run('plot_bootstrap.py')
        run('bounds/plot_global_height.py')
        run('search/plot_search.py')


if __name__=='__main__':main()
