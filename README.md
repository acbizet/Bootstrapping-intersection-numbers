# Phase-conditioned contour-height bootstrap

Implemented for the exact free U(2) oscillator product at P=10^7 and resonance
n=26. The new result is a global finite-charge phase exclusion at Q=5883,
in addition to sharper global height exclusions at Q=5880 and Q=5882.
Near c_infinity, Q=5940 has an explicitly validated contour point satisfying
both necessary conditions. No nonzero thimble intersection is established.

## Conventions

    x=exp(omega), y=exp((omega+eta)/3), U=exp(iu), eta=i*pi+gamma,
    P=2J+Q, L=P-Q, R_Q=3Q, lambda=gamma/omega^2.

All targets have physical integer P,Q. The original contour has independent
physical phases theta_x,theta_y,u, each modulo 2*pi. Its real source radii
are optimized using the positive Fock product. The full exponent is

    S=F-P*omega-Q*eta.

F includes Haar and all oscillator factors; the external constant Weyl factor
1/2 is consistently omitted from both heights. No fixed-lambda source path
is imposed, and no winding regulator Lambda is introduced.

The coarse ceiling B0=log4+F_Fock-P*omega_+-Q*eta_+ bounds the whole contour.
It is not its actual maximum. At c_infinity only its leading P^(2/3) term
matches the resonant action. The finite coarse equality is Q=5879.9614...
at P=10^7, whereas c_infinity is approximately .12798394392465.

## What the bootstrap does

1. Vector-boson modulus losses confine every point high enough to meet the
   saddle to |u|<2*pi/3. This is a global inequality, not a sampling claim.
2. In that sector every bosonic adjoint block is largest at zero chemical
   phases. The fermionic four-weight polynomial has nonnegative coefficients
   and obeys the same bound. Thus the global magnitude maximum reduces to a
   one-dimensional positive-chemical gauge profile.
3. With s=sin^2(u/2), the scalar profile has a unique maximum, proved by a
   monotone derivative. Positive omitted losses and explicit geometric tail
   bounds provide conservative height and gauge-angle enclosures.
4. For a surviving target, the remaining height allowance bounds vector and
   chiral chemical-phase losses. Closed geometric sums then enclose every
   competitive point in a small physical theta_x,theta_y,u box.
5. Source-angle Taylor estimates, with exact Fock linear cancellation and
   explicit infinite tails, bound every action phase in that box. Upward
   flow preserves the action phase. A mismatch therefore gives exclusion
   without tracing a complete upward manifold.

The full proofs and formulas are in bounds/README.md and
phase_enclosure/README.md. The logical scope of each numerical result is
preserved in summary.json; no local optimizer is used as a global ceiling.

## Results

| Q | Global maximum minus saddle height | Outcome |
|---:|---:|---|
| 5880 | -12.76402727 | Excluded globally by height |
| 5882 | -0.99794554 | Excluded globally by height |
| 5883 | +4.88434505 | Excluded globally by phase |
| 5884 | +10.76613557 | Global test inconclusive |
| 5940 | +339.35096776 | Necessary height-and-phase test passes |

For Q=5883 every sufficiently high original-contour point satisfies

    |Im S| modulo 2*pi <= 2.171907177,

whereas the saddle phase has distance 2.716790508 from 2*pi*Z. The global
phase gap exceeds .54488 radians. This excludes the tested ordinary resonant
saddle, conditional on the converged numerical target/radius inputs.

For Q=5884 a local constrained optimizer finds its best same-phase candidate
9.74069 below the saddle. The current global phase enclosure is too loose to
turn that into an exclusion. This case is explicitly left unresolved.

For Q=5940, c=.127973420588 is a physical integer-charge target near
c_infinity. An independently refined original-contour point has exactly the
required phase within numerical precision and height 318.94090896 above the
saddle. This point proves that this necessary test cannot exclude the target.
It does not prove an upward path reaches that point or a nonzero net signed
intersection count. There is no established selection of lambda at c_infinity.

## Validation and limits

- Fresh 60-digit saddle solves, with 75-digit enlarged-prefix checks for
  Q5883 and Q5940, reproduce the actions to better than 57 absolute digits.
- Scalar global bounds are independently checked at 60 digits with 4096
  oscillator levels; production uses 8192 levels and tails below 1.9e-38.
- The phase enclosure is repeated at 90 digits and 6000 levels. The Q5883
  bound changes by about 8.1e-19, far below the .54488-radian gap.
- Two independent formula audits check the global loss reductions, angle
  localization, phase Taylor factors and geometric tails.
- Local constrained candidates are refined by independent 55/75-digit
  mpmath calculations, with different internal coordinate parametrizations.

These are analytical global inequalities evaluated at numerically solved
inputs. They are not a complete interval-arithmetic certificate of the exact
critical point and optimized radii. The classification is stable under the
reported precision and cutoff changes.

The global searches use exact oscillator prefixes. The Euler--Maclaurin
evaluator is used only in its small principal-omega sector. A discovered
pitfall is documented: blindly shifting raw omega by 2*pi*i in that tail
evaluator does not give a valid periodicity test. The physical periodicity
checks instead use an independent exact-prefix product. No affected shifted
tail value enters any reported bound or witness.

The new exclusion concerns the ordinary critical-point contribution of the
tested branch. Singular residue cycles, other saddles and cancellations
between different contributions are separate questions. The earlier local
pole-residue-circle intersection is compatible with these results because
that circle is a different cycle from the original coefficient torus.

## Reproduction

Install requirements.txt, then run

    python reproduce.py

The default recomputes the inexpensive global bounds and phase enclosure,
checks the saved data, and redraws the figures. To redo the source solves,
local searches and independent high-precision candidate refinements as well:

    python reproduce.py --full

Every script resolves its paths relative to itself. No previous repository
or network connection is needed after installing dependencies. Individual
folders include their own entry points and detailed limitations.
