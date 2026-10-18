# Local phase-conditioned height search

This standalone finite multistart search compares the resonant saddle height
with heights attained on the original coefficient torus. It does not compute
flows or intersection numbers, and its optimized maxima are not global upper
bounds.

Independent physical angles are theta_x, theta_y and real u:

    omega = w0 + i theta_x
    gamma = h0 + i(3 theta_y - theta_x - pi)
    S = F - P omega - Q(i pi + gamma)

Here w0 and h0 are the source radii in the adjacent target files. The full
oscillator prefix includes Haar and Cartan factors, and omits the common
constant Weyl factor 1/2 just as in the target action. All three sources remain
independent in the analytic gradients and Hessian.

The local positive peak is obtained by solving the gauge derivative and then
checking seven unconstrained optimizer starts. The negative real height Hessian
sets whitening coordinates. Twenty-four sphere directions at radii
2, 5, 10, 20, 40, 70 supply constrained starting points. At least three starts
are tried on each of the unwrapped phase branches phi, phi+2pi, phi-2pi.
Only positive u is searched; gauge inversion is an exact symmetry. Coverage
of the rest of the physical periodic torus is not established by this code.

The phase constraint conserves the integrand phase, not a selected principal
logarithm sheet. Locally the continuous action phase relative to the positive
peak is unambiguous. The phase itself need not be unity. Each retained
candidate is checked against the original resonant phase modulo 2pi, and the
normal-coordinate residual is zero by construction. Analytic KKT residuals
and the two tangent Lagrangian-Hessian eigenvalues are stored.

For P=1e7 and n=26 the computed local margins are:

| Q | local peak minus saddle | same-phase local maximum minus saddle |
|---:|---:|---:|
| 5880 | -12.76402727 | -33.27771366 |
| 5882 | -0.99794554 | -21.50820008 |
| 5883 | 4.88434505 | -66.00600469 |
| 5884 | 10.76613557 | -9.74068558 |
| 5940 | 339.35096775 | 318.94090896 |

Thus the resolved local high-height component has no eligible same-phase point
for 5883 or 5884. This does not by itself exclude other components elsewhere on
the torus. At 5940 a same-phase original-torus point above the resonant saddle
is explicitly found. That fact supplies a necessary condition for an upward
intersection, not a sufficient one.

The odd-Q case 5883 has target phase +2.71679051 rather than approximately
-0.425: both branches of the modulo-2pi condition were included. This explains
its larger local phase cost and is not an optimizer discontinuity.

`results.json` is the aggregate, and `result_Q*.json` retain all phase starts,
sampled points and convergence checks. Doubling the retained oscillator levels
from 3200 to 6400 changes the final complex actions by less than 6.3e-9.
These are floating-point convergence checks, not interval certificates.
Independent arbitrary-precision validation is supplied elsewhere in this
repository.

Recompute from any current directory with numpy, scipy and matplotlib:

    python search_phase_height.py
    python plot_search.py

The target input defaults to ../targets and can be changed with --targets.
A subset can be selected, e.g. --q 5883 5940. Plot regeneration alone reads
the saved aggregate and does not solve or optimize anything.
