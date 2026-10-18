# A finite global phase bootstrap

This argument bounds every original-contour point with height at least the
resonant saddle height. It does not infer a global bound from sampled flows or
from a local optimizer. Exact infinite-product inequalities are evaluated with
converged numerical saddle/radius inputs and explicit positive tail bounds.
There is no interval certificate of the input saddle.

## Remaining height allowance

Write theta_x=arg(x), theta_y=arg(y), with independent physical periods 2*pi,
and reduce the gauge angle to |u|<=pi. The previous block bootstrap confines
every competitive point to |u|<=U<2*pi/3. For fixed u all oscillator blocks
have their largest chemical-phase modulus at theta_x=theta_y=0. Therefore

    B0-Re S = D0(u) + C(theta_x,theta_y,u),   C>=0.

If D_min is a lower bound on min D0 and M=B0-Re S_res, every possible
intersection satisfies C<=A=M-D_min. Scalar prefix/tail estimates from
`../bounds/` provide A and U. The code enlarges their numerical values outward
before using them. These paddings are engineering margins, not interval proofs.

## Global theta_x enclosure without an angular grid

The first bosonic Fourier winding alone implies

    C >= kappa * sum_(k>=2) (k-1)x^k [1-cos(k theta_x)],
    kappa = 4 cos^2(U/2),   x=exp(omega_+)>0.

All omitted Fourier windings and all other blocks have nonnegative excess
loss. With v=1-cos(theta_x) and a=(1-x)^2, the retained sum is exactly

    f(v) = kappa*x^2 * [(4+2x)v+(4x^2/a-2)v^2]/(a+2xv)^2.

Its derivative has the sign of (4+2x)a-4v. Thus f first increases, reaches
one maximum at v_*=(1+x/2)a, then decreases. If f(2)>A, as explicitly checked
for every target here, the entire descending branch is excluded. Solving the
first root f(v)=A yields the global bound |theta_x|<=X, modulo 2*pi. This
rules out remote chemical-phase components rather than merely sampling them.

## Global theta_y enclosure

A fermionic adjoint block is

    D(z)=(1+z)^2(1+2 cos(u) z+z^2).

For |u|<=U<2*pi/3 it has nonnegative coefficients. Normalize its coefficients
at z=r to probabilities p_j. Then

    -log |D(r exp(i phi))/D(r)|
       >= p_0 p_1 (1-cos phi)
       >= kappa*r/(1+r)^8 * (1-cos phi).

The first inequality follows from the exact squared characteristic-function
identity and -log(1-z)>=z. Retain the chiral fermions, with
r=x^(k+1)y, phi=(k+1)theta_x+theta_y and multiplicity 3(k+1). Summing their
geometric first-harmonic bound gives

    C >= K [W-Re(exp(i theta_y) T(theta_x))],
    K=3*kappa*x*y/(1+x*y)^8,
    W=(1-x)^(-2),   T(theta_x)=exp(i theta_x)/(1-x exp(i theta_x))^2.

For |theta_x|<=X,

    |T|<=W,
    |arg T|<=X+2 atan(x X/(1-x)).

When 1-A/(KW)>=0, these inequalities give

    |theta_y| <= Y = X+2 atan(x X/(1-x)) + acos(1-A/(KW)),

with the result capped at pi. If the required real part is negative, the
routine conservatively retains Y=pi. All current targets have a positive
required real part. This step also reduces the full physical theta_y period,
not a restricted omega/eta phase cover.

## Phase in the resulting box

Haar is positive on the real gauge circle and supplies no phase. At positive
chemical sources the oscillator product is real positive. Expand its phase in
the source angles at fixed u. Fock optimization cancels the linear source
terms at u=0, leaving

    Psi = (F_L(u)-L) theta_x + (F_R(u)-R) theta_y
          + a third-order source remainder.

For a tower of multiplicity m, fugacity modulus r, physical source charges
(k_x,k_y), define f=m*r*(1+r)/(1-r)^3. The logarithmic third derivative has
absolute value bounded by r*(1+r)/(1-r)^3 for every intervening angle.
Gauge inversion cancels the first u derivative. Thus

    |Psi| <= U^2 * sum_towers f (k_x X+k_y Y)
             + (4/6) sum_towers f (k_x X+k_y Y)^3
             + explicit tails + numerical Fock-source residual allowance.

The factor U^2 includes the two charged gauge weights and the Taylor factor
1/2. The factor 4/6 includes all four adjoint weights and the third-order
Taylor factor. The physical towers are

| tower | r | m | k_x | k_y |
|---|---|---:|---:|---:|
| vector boson | x^(k+1) | k | k+1 | 0 |
| vector fermion | y^3 x^k | k+2 | k | 3 |
| chiral fermion | x^(k+1)y | 3(k+1) | k+1 | 1 |
| chiral boson | x^k y^2 | 3(k+1) | k | 2 |

Infinite tails are bounded by exact geometric moments through degree four.
The sources and all scalar inequalities remain independent of the resonant
lambda relation. No fixed-lambda path is imposed.

## Decisive target and limitations

For P=10^7,Q=5883,n=26, the global high-height region obeys

    |u|       <= .014169957,
    |theta_x| <= .000077535,
    |theta_y| <= .058698033,
    |Psi|     <= 2.171907177 radians.

The saddle phase has distance 2.716790508 radians from 2*pi*Z. Since the
phase bound is less than pi, no hidden phase wrapping changes the comparison.
The global phase gap is more than .54488 radians. Upward-flow phase
conservation therefore excludes this ordinary resonant saddle, conditional
on the numerically validated target values.

The same conservative bound does not exclude Q5884 or Q5940. A local phase
optimizer finds a negative height margin for Q5884, but this global enclosure
does not certify that optimizer's supremum. Q5940 has a separately validated
phase-matched contour point above the saddle, so the necessary test genuinely
passes there. Passing does not establish an upward attachment or its signed
intersection coefficient.

Run `python global_phase_bound.py` from any working directory. Results are
saved in results.json. Independent longer-prefix checks and formula audits
are recorded under the enclosing validation folder.
