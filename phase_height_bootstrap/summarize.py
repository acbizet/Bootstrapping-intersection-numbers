"""Assemble the outcome without promoting local maxima to global bounds."""
from pathlib import Path
import json
from mpmath import mp

HERE=Path(__file__).resolve().parent


def build():
    height=json.loads((HERE/'bounds/global_height_results.json').read_text())['rows']
    phase={r['Q']:r for r in json.loads((HERE/'phase_enclosure/results.json').read_text())['rows']}
    search={int(r['Q']):r for r in json.loads((HERE/'search/results.json').read_text())['cases']}
    validation={r['Q']:r for r in json.loads((HERE/'validation/results.json').read_text())['rows']}
    rows=[]
    for h in height:
        q=int(h['Q']);t=json.loads((HERE/f'targets/target_Q{q}.json').read_text())
        p=phase.get(q);s=search[q];v=validation.get(q)
        witness=v['refined_local_phase_extremum'] if v else s['winning_phase_candidate']
        witness_margin=float(witness['height_margin'])
        if h['height_excluded']:status='excluded_by_global_height_bound'
        elif p and p['phase_excluded']:status='excluded_by_global_phase_bound'
        elif witness_margin>0:status='necessary_height_and_phase_test_passes'
        else:status='inconclusive_local_search_only'
        rows.append(dict(P=10_000_000,Q=q,n=26,c=float(t['c']),
            coarse_height_allowance=float(t['B0_minus_ReS']),
            global_height_margin_upper=h['maximum_height_minus_hres_upper'],
            saddle_phase=float(t['phase_mod_2pi']),
            global_phase_absolute_bound=None if p is None else float(p['absolute_unwrapped_phase_bound']),
            global_phase_exclusion_margin=None if p is None else float(p['phase_exclusion_margin']),
            local_same_phase_candidate_margin=witness_margin,status=status,
            global_resonant_intersection_number=None,interval_certified=False))
    result=dict(P=10_000_000,n=26,c_infinity='0.1279839439246513804388040555770683',
        finite_coarse_height_boundary_Q='5879.96141272442083988723049475',rows=rows,
        interpretation='Exact global inequalities evaluated at converged numerical targets; no interval certificate. A positive phase-height witness does not prove attachment.',
        main_result='Q5880 and Q5882 excluded by global height; Q5883 excluded by global phase; Q5884 unresolved; Q5940 has a necessary-condition witness.',
        global_nonzero_intersection_established=False)
    (HERE/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps(build(),indent=2))
