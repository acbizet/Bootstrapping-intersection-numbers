# Global height bootstrap from oscillator blocks

This calculation replaces a local peak search by an analytic inequality on
the entire original three-angle coefficient torus. It does not integrate flows
or assume that every large point lies in a Gaussian neighbourhood.

The inequalities below are exact for the infinite free U(2) oscillator product.
Their finite-target evaluations use numerical radii and critical values rather
than interval-certified solutions. `scalar_validation.json` independently
evaluates the final scalar brackets with mpmath; this is an arithmetic and
cutoff convergence check, not a formal interval certificate of the saddle.

## Global reduction

Let `s=sin^2(u/2)`, with the real gauge angle reduced to `|u|<=pi`. Let
`B0=log(4)+F_Fock-P*omega-Q*eta`, and let `M=B0-Re(S_res)`.

A bosonic adjoint oscillator of modulus r and chemical phase phi has four
weights, two diagonal and a conjugate charged pair. Its loss relative to the
Fock value is

    sum_{l>=1} r^l/l * [4-(2+2*cos(l*u))*cos(l*phi)].

Every coefficient `2+2*cos(l*u)` is nonnegative. Consequently its modulus
is globally maximal at chemical phase zero for every fixed real u.

A coarse, independent confinement follows from the single-weight inequality

    loss_boson(r,phi) >= atanh(r)*(1-cos(phi)).

The vector bosons alone therefore give

    B0-Re(S) >= 4*C*s-log(s),
    C=sum_{j>=2} (j-1)*atanh(exp(j*omega)).

We retain only a finite positive prefix in C. For every target here, the
right side at s=3/4 is far above M and is increasing for s>=3/4.
Thus every potentially intersecting contour point has `s<3/4`, or
`|u|<2*pi/3`, without any assumption about its chemical angles.

For a fermionic adjoint oscillator put `z=r*exp(i*phi)`. Its product is

    (1+z)^2*(1+2*cos(u)*z+z^2).

The polynomial coefficients are

    [1, 2+2*cos(u), 2+4*cos(u), 2+2*cos(u), 1].

All are nonnegative when `|u|<=2*pi/3`. The triangle inequality therefore
bounds its modulus by the value at phi zero. It follows that, throughout
the only region that can attain the resonant height, **every oscillator
block is simultaneously largest at positive chemical phases**.

The remaining global maximum is one-dimensional. With multiplicity m
per adjoint weight, define

    A_b=4*r_b/(1-r_b)^2,   A_f=4*r_f/(1+r_f)^2,
    D(s)=sum_b m_b*log(1+A_b*s)
         -sum_f m_f*log(1-A_f*s)-log(s).

Then the actual global torus maximum is `B0-min D`, provided its optimizer
lies in the confined region, as checked here. Haar is included in `-log(s)`.

## Unique scalar maximum and finite bounds

Multiplication of the derivative by s gives

    s*D'(s)=-1+sum_b m_b*A_b*s/(1+A_b*s)
                 +sum_f m_f*A_f*s/(1-A_f*s).

This expression is strictly increasing, from -1 to a positive value. The
unique zero gives the unique scalar minimum. No numerical global optimizer
is needed.

All omitted oscillator losses are positive. If the zero of the full derivative
is bracketed by `[s_lo,s_hi]`, a valid lower bound on its minimum is

    min D >= L_prefix(s_lo)-log(s_hi).

The discarded derivative and loss tails have explicit geometric upper bounds.
For a tower with `r_k=r0*a^k` and multiplicity `m*(k+d)`,

    sum_{k>=N} m*(k+d)*r_k
      =m*r_N*((N+d)/(1-a)+a/(1-a)^2).

Multiplication by `4*s/(1-r_N)^2` for bosons, or by
`4*s/(1-4*r_N*s)` for fermions, bounds both the tail loss and its contribution
to `s*D'(s)`. These bounds check the bracket signs and give a two-sided height
enclosure.

The prefix sublevel set `D_prefix(s)<=M` also contains every possible
intersection gauge angle, furnishing the `allowed_gauge_enclosure` output.

## Results and limitations

At P=10^7, n=26, the Q=5880 and Q=5882 targets are globally excluded by height.
Their original coarse ceilings exceeded the critical height, but the global
oscillator/Haar maximum is lower. Q=5883, Q=5884, and Q=5940 are not excluded
by this height test. A phase bootstrap or other information is required there.

An inequality proved on the complete torus does not establish a nonzero
intersection when its height test passes. It only excludes when the maximum
lies strictly below the saddle height.

Run from any working directory:

    python /absolute/path/to/bounds/global_height.py
    python /absolute/path/to/bounds/validate_scalar.py

