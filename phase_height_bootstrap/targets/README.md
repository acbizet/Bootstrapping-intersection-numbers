# Exact-product targets for the phase-height bootstrap

Run from any current directory:

```sh
python /absolute/path/to/targets/compute_targets.py
```

Requires `mpmath`. Source files under `src/` are standalone frozen dependencies;
the script does not import another output repository or read an earlier target.
The default recomputes five targets at 60 decimal digits and repeats Q=5883 and
Q=5940 at 75 digits, with an enlarged exact oscillator prefix and more
Euler--Maclaurin terms. The logarithmic product is evaluated with all source
derivatives independent before any near-index relation is inspected.

Coordinates and conventions are

```
z = (omega, gamma, u)
eta = i*pi + gamma
x = exp(omega), y = exp((omega+eta)/3), U = exp(i*u)
S = F - P*omega - Q*(i*pi+gamma)
```

The full F includes Haar and omits the external constant Weyl factor 1/2, as
does the stored positive-torus ceiling B0. Subtracting log(2) consistently from
both leaves every height gap unchanged. The three complex saddle equations are
F_omega=P, F_gamma=Q, F_u=0. The integer n=26 selects the initial resonance.
No reality condition on lambda=gamma/omega^2 is imposed.

`target_Q*.json` contains full complex sources, action, source Hessian H, wrapped
phase, Fock-optimized positive radii, positive Fock Hessian, ceiling B0, numerical
residuals and truncation parameters. Complex numbers are objects with decimal
string `re` and `im`; reals are decimal strings. The physical coefficient torus
has omega=w0+i*theta_x, eta=h0+i*(3*theta_y-theta_x), and real u.

The optimized ceiling is

```
B0 = log(4) + Fock(w0,h0) - P*w0 - Q*h0.
```

This is an algebraic upper ceiling, not a claim that all individual factor
inequalities saturate at the same contour point. Wrapped action phases are
computed modulo 2*pi; unwrapped logarithms depend on summand branch choices.

`frozen_asymptotic_constants.json` and `frozen_finite_threshold.json` record the
earlier c_infinity and the finite-Q ceiling equality. They are comparison data,
not solver inputs. In particular, c_infinity equality is only leading order;
the finite ceiling equality at P=1e7,n=26 lies near Q=5879.9614.

All precision checks are numerical convergence checks. None is an interval
certificate or a computation of a global contour intersection number.
