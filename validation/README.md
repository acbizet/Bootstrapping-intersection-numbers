# Independent validation

`validate_candidates.py` independently evaluates and refines the previously
found local phase-constrained candidates at `P=10000000`, `Q=5883,5940`.
It does not perform another search.

The two configurations use `FullProduct` resonance coordinates `n=7` and
`n=11` rather than the physical saddle label `n=26`; this is a change of
coordinates only. They use 55 and 75 requested digits, exact oscillator
prefixes expanded by 256 and 512 levels, and different Euler--Maclaurin
orders. The full original-contour source radii are retained from the
75-digit target files. Independent source derivatives are transformed to
physical angles before the Lagrange equations are solved.

Results:

| Q | Refined local phase-conditioned height minus saddle height | KKT residual |
|---|---:|---:|
| 5883 | -66.0060046852667305902 | 1.79e-47 |
| 5940 | +318.940908957767741345 | 5.06e-58 |

Both constrained tangent curvatures are negative at each refined point.
These are well resolved **local maxima** on the selected phase sheet. They
are not certified global maxima. The `Q=5940` point is an explicit original-
contour point above the saddle height with the required phase; it shows the
necessary phase-height condition is satisfied there. It does not prove that
an upward path from the saddle reaches that point.

Reevaluation with the larger configuration changes these local actions by
less than `5e-85`; these are floating-point convergence checks, not interval
certificates. Agreement this far below the requested digits reflects the
guard digits retained internally by `FullProduct`.

Gauge inversion and chemical conjugation identities are checked separately.
The physical `theta_x -> theta_x+2pi` identity is checked in an independent
exact 96-level product. That finite product is used only for its exact
periodicity identity, not as an approximation to the converged full height.

Important engine restriction: the Euler--Maclaurin remainder is expanded in
the raw oscillator step `omega`, and must be used on a suitable small
principal branch. Evaluating it directly at `omega+2pi*i` gives a spurious
large remainder even though the exact product is periodic. Global angle
searches must use the exact-prefix engine, a controlled Fourier expansion,
or physically reduce the source angles before the EM evaluation. All
validated local candidates have `|Im omega|<2e-5` and are unaffected.

Run from any directory:

```
python /path/to/phase_height_bootstrap/validation/validate_candidates.py
```

Dependencies: `mpmath`. Inputs are the sibling `targets/` and `search/`
outputs; the local `src/` is standalone. `mathematical_audit.md` records the
exclusion theorem, correct periodicity, and an independent proof of the
one-dimensional gauge bootstrap.

## Global phase inequality audit and expanded recomputation

`validate_phase_enclosure.py` independently checks the algebra of the angle
localization and reevaluates `phase_enclosure/global_phase_bound.py` at 90
digits and 6000 oscillator levels, compared with 70 digits and 5000 levels.
The exact proof, factor audit and qualifications are recorded in
`mathematical_audit.md`.

For `Q=5883` the expanded result is

```
global absolute phase bound = 2.17190717653105698595225816956
nearest saddle phase       = 2.71679050775250949813070210493
phase exclusion margin     = 0.54488333122145251217844393537
```

The bound changes by `8.02e-19`; the closed rational vector loss agrees with
an independent complex geometric-sum form within `8.3e-83`. The prerequisite
argument for the chiral angle enclosure is `0.9990047712...>0`. The audit
identified that this positivity is needed; the general implementation now
keeps the full chemical period when it is unavailable.

The expanded bounds for `Q=5884,5940` remain inconclusive, as expected. The
largest numerical bound change among these cases is `4.47e-16`. None of
these calculations is interval-certified. They validate a global analytic
inequality evaluated on numerically supplied, conservatively padded scalar
enclosures, rather than promoting the local optimizer to a global bound.
