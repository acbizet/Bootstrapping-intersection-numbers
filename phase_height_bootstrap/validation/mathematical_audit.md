# Independent mathematical audit of the phase-height test

## Exclusion statement

Let `S` be the logarithm of a nonzero meromorphic integrand on a regular
neighbourhood of a flow path. Its real part `h=Re S` and the phase `exp(i Im S)`
are independent of logarithm branches. For any constant positive Hermitian
metric, a nonconstant upward gradient trajectory satisfies

\[
\frac{dh}{dt}=\partial S\,M\,\overline{\partial S}>0,
\qquad \frac{d\,\operatorname{Im}S}{dt}=0.
\]

If it emerges from a critical point `sigma`, every other point on the
trajectory has `h>h_sigma` and the same integrand phase. This is true even when
the saddle is approached only as `t -> -infinity`.

Consequently an upper bound `U_phi` for the height of **all** original-contour
points with the saddle phase gives

\[
U_\phi\le h_\sigma
\quad\Longrightarrow\quad
K_\sigma\cap C=\varnothing,
\]

provided the saddle itself is off `C`. Exact equality suffices: ascent is
strict after leaving the saddle. This assertion needs neither a numerical
flow calculation nor a Morse assumption. It does require a globally valid
upper bound. A discovered constrained maximum is a **lower bound** on the
unknown constrained supremum and cannot establish exclusion by itself.

An attained maximum is unnecessary if one has a bound on the supremum. Points
where the integrand vanishes have `h=-infinity`; their undefined phase causes
no difficulty for a finite high-height test. A pole on the original contour
would invalidate the finite-ceiling assumption and must be excluded first.

## Correct periodicity and derivatives

The safe independent physical source angles are `t_x,t_y`:

\[
\omega=w+i t_x,\qquad
\eta=h+i(3t_y-t_x),\qquad
\gamma=\eta-i\pi,\qquad u\in\mathbb R.
\]

Both angles have period `2pi`. The extraction phase is
`-(P-Q)t_x-3Q t_y`; it is periodic when `P-Q` and `3Q` are integers. This
holds for the integer charge targets in this experiment. In contrast,
independent `2pi` shifts of `Im omega` and `Im eta` are not both symmetries:
their shifts `2pi(a,b)` must obey `a+b=0 mod 3`. Conservative independent
periods in those latter coordinates are `6pi`.

For independent-source derivatives from the full product, the contour action
gradient is

\[
\partial_{t_x}S=i(F_\omega-F_\gamma-P+Q),\quad
\partial_{t_y}S=3i(F_\gamma-Q),\quad
\partial_uS=F_u.
\]

The saddle phase must be matched modulo `2pi`, not forced to unity. Direct
principal logarithms can differ by integer multiples of `2pi i`; comparisons
of fixed points use `atan2(sin(delta),cos(delta))`. A local continuous phase
coordinate is adequate only after its branch is explicitly fixed.

## What the old local experiment established

The previous `global_contour_intersections/phase_height` experiment used
`P=22078774,Q=10000`. Its local high-height component has numerical phase
extrema `+/-0.352828393306837...`, separated from the saddle phase
`-0.963265623951152...`. The stationary equations and local tangent
curvatures were checked to high precision. Those are local statements;
neither global periodic coverage nor interval bounds were established.

The present calculation must separately provide: a global localization bound,
coverage of every surviving angle box and periodic image, and validated
upper bounds in those boxes. Finite multi-start searches and dps/cutoff
doubling support numerical evidence only. They do not replace that proof.

## A useful one-factor modulus inequality

For `0<a<1`,

\[
\log\frac{|1-ae^{i\theta}|}{1-a}
=\frac12\log\left(1+\frac{4a\sin^2(\theta/2)}{(1-a)^2}\right),
\]

and for `a>0`,

\[
\log\frac{1+a}{|1+ae^{i\theta}|}
=-\frac12\log\left(1-\frac{4a\sin^2(\theta/2)}{(1+a)^2}\right).
\]

Both are nonnegative. They give exact losses from the positive oscillator
ceiling for bosonic denominator and fermionic numerator factors. Keeping any
selected subset yields an upper bound because all discarded losses remain
nonnegative. The gauge Haar factor is exactly `4 sin^2(u/2)`, so its loss
relative to the crude Haar ceiling `4` is `-log(sin^2(u/2))`. This diverges
where the positive oscillator ceiling would be saturated (`u=0 mod 2pi`).
The resulting strict Haar-oscillator tradeoff is therefore a useful first
bootstrap before imposing the action phase.

These are exact algebraic inequalities. Evaluating their optimized numerical
values with ordinary or arbitrary-precision floating point is not itself an
interval certificate.

## Independent audit of the global gauge bootstrap

For one U(2) adjoint bosonic block of multiplicity `m`, radius `0<r<1`,
chemical phase `phi`, and gauge angle `u`, the loss relative to its positive
oscillator ceiling is

\[
L_b=m\sum_{\ell\ge1}\frac{r^\ell}{\ell}
 [4-(2+2\cos\ell u)\cos\ell\phi].
\]

Because `2+2 cos(ell u)>=0`, this loss is minimized at `phi=0` for fixed `u`.
It follows that

\[
L_b\ge2m\log\frac{|1-r e^{iu}|}{1-r}
\ge2m\operatorname{atanh}(r)(1-\cos u).
\]

The second inequality follows from concavity of `log(1+A s)` on `0<=s<=1`,
with `s=sin^2(u/2)` and `A=4r/(1-r)^2`; it is equality at `u=0,pi`.
Keeping just the vector bosons gives the globally valid loss

\[
L\ge A_x(1-\cos u),\qquad
A_x=2\sum_{k\ge2}(k-1)\operatorname{atanh}(x^k).
\]

All other oscillator losses and the Haar loss remain nonnegative. If a
competitive point must have total loss at most `M` and `M<3 A_x/2`, its gauge
angle is confined to `|u|<2pi/3` modulo `2pi`. No chemical-angle sampling is
involved in this localization.

For each fermionic adjoint block, the product of its four gauge weights is

\[
D(z)=(1+z)^2(1+2\cos u\,z+z^2).
\]

Its five coefficients are

\[
(1,\ 2+2\cos u,\ 2+4\cos u,\ 2+2\cos u,\ 1).
\]

They are all nonnegative in the localized gauge sector. The triangle
inequality therefore gives `|D(r exp(i phi))|<=D(r)`: each fermion block is
also maximized at zero chemical phase, for fixed `u`.

All oscillator blocks can attain their fixed-`u` maxima simultaneously at
`t_x=t_y=0`. Therefore, after excluding the complementary gauge sector with
the coarse loss bound, the global original-torus magnitude maximum is exactly
the maximum of the one-dimensional positive-chemical gauge profile.

This argument is analytic and global. The final numerical maximization of
that one-dimensional profile still needs derivative or interval bounds if
its numerical value is to be used as a rigorous upper bound. A coarse
one-dimensional analytic envelope may already suffice for exclusion near the
boundary without locating the exact maximum.

## Audit of the global phase enclosure

Assume all competitive points have `|u|<=U<2pi/3` and their chemical-angle
height loss above the zero-chemical profile is at most `A`. Set
`kappa=4 cos^2(U/2)`.

The first winding of the vector-boson chemical loss bounds

\[
\kappa\sum_{k\ge2}(k-1)x^k(1-\cos k\theta_x)\le A.
\]

Writing `v=1-cos(theta_x)` and `a=(1-x)^2`, the sum is exactly

\[
\kappa x^2\frac{(4+2x)v+(4x^2/a-2)v^2}{(a+2xv)^2}.
\]

It rises to its sole maximum at `v=(1+x/2)a`, then falls. If its endpoint at
`v=2` still exceeds `A`, the first root alone bounds every surviving physical
angle, modulo `2pi`. The independent validation compares this rational form
against `kappa*(x^2/(1-x)^2-Re[z^2/(1-z)^2])`, `z=x exp(i theta)`.

For a normalized positive-coefficient fermion block `D`, the zero- and
one-power probabilities satisfy

\[
p_0p_1\ge\frac{\kappa r}{(1+r)^8}.
\]

Expanding `|sum p_j exp(i j phi)|^2` and keeping only the `p0 p1` pair gives
`-log|D(r exp(i phi))/D(r)|>=p0 p1(1-cos phi)`. For chiral fermions this
implies

\[
K\left[W-\operatorname{Re}
\frac{e^{i(\theta_y+\theta_x)}}{(1-xe^{i\theta_x})^2}\right]\le A,
\quad W=(1-x)^{-2},\quad K=\frac{3\kappa xy}{(1+xy)^8}.
\]

Using the modulus bound `|e^{i theta_x}/(1-x exp(i theta_x))^2|<=W` gives the
implemented `theta_y` bound **only when** `1-A/(K W)>=0`. Otherwise the safe
fallback is the whole `theta_y` period. This qualification was found in the
independent audit and added to the main implementation; all present targets
satisfy the positive-argument condition by a wide margin.

For the phase step, each adjoint block is chemical-conjugation symmetric at
fixed real `u`. Its phase has no quadratic term about zero chemical angles.
The universal third derivative bound is

\[
|d_\phi^3\log(1\pm r e^{i\phi})|
\le \frac{r(1+r)}{(1-r)^3},\qquad 0<r<1.
\]

Four gauge weights and Taylor's factor `1/6` give exactly the implemented
`2/3` cubic coefficient. The linear phase difference from the optimized Fock
charges is bounded by `u^2` times the same third-derivative sums: the two
opposite gauge weights cancel the term linear in `u`, and their factor two
cancels the second-order Taylor denominator. The neutral weights introduce
no extra factor. The source residual allowance is a conservative factor ten
times `epsilon*((P+Q)X+3QY)`.

Finally, in the omitted oscillator tail each radius is at most `x^k`, each
chemical direction is bounded by `(k+2)X+3Y`, and the total multiplicity is
`8(k+1)`. These inequalities verify that the implemented linear and cubic
tail coefficients are conservative. The four stored geometric moments are
the exact sums of `(k+2)^p x^k` for `k>=N`, `p=1,...,4`.

No uncontrolled phase branch is needed: all surviving chemical angles lie
in the same small physical box about zero. At fixed nonzero real `u`, all
factor radii are below one, so the straight chemical segment encounters no
zero or pole. It defines a continuous phase with value zero at its origin.
The final exclusion compares its absolute bound with the nearest saddle
phase modulo `2pi`.

The derivation is global once the scalar gauge enclosures and remaining-loss
budget are valid. The saved numerical constants and floating-point
evaluations are not interval-certified. This last distinction remains even
after the long-prefix/high-precision validation.
