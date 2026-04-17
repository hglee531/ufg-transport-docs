# Sn — slab geometry

The $S_N$ solver in UFG-Transport discretises the 1-D slab transport equation in
angle using Gauss–Legendre quadrature and in space using **diamond difference**. It
supports anisotropic scattering to arbitrary Legendre order and three boundary
conditions: reflective, vacuum, and white (isotropic return).

The full derivation — including quadrature weight construction via Golub–Welsch and
the complete sweep algorithm — is in `docs/solver_theory_and_algorithms.md`. This page
summarises the equations, the algorithm, and the knobs the user can turn.

## Governing equation

Within a group $g$, the angular flux $\psi_g(x,\mu)$ along direction $\mu$ satisfies

$$
\mu\,\frac{\partial \psi_g}{\partial x}(x,\mu) + \Sigma_{t,g}(x)\,\psi_g(x,\mu)
\;=\;
\sum_{\ell=0}^{L} \frac{2\ell+1}{2}\,P_\ell(\mu)\,q_{\ell,g}(x),
$$

where the source moment absorbs scattering, fission, and the external source:

$$
q_{\ell,g}(x) \;=\; \sum_{g'}\Sigma_{s,\ell}^{g'\to g}(x)\,\phi_{\ell,g'}(x)
\;+\; \delta_{\ell 0}\,\bigl(\chi_g Q_\text{fis}(x) + Q_{\text{ext},g}(x)\bigr).
$$

The angular flux is integrated to recover scalar moments

$$
\phi_{\ell,g}(x) \;=\; \int_{-1}^{1} P_\ell(\mu)\,\psi_g(x,\mu)\,d\mu
\;\approx\; \sum_{n=1}^{N} w_n\,P_\ell(\mu_n)\,\psi_g(x,\mu_n).
$$

## Angular quadrature

UFG-Transport uses **symmetric Gauss–Legendre** quadrature: the $\mu_n$ are the
roots of $P_N(\mu)$ and the $w_n$ are the matching weights, computed at run time via
the Golub–Welsch algorithm (`src/utils/quadrature.cpp`). $N$ is chosen by the user
via `--sn-order` and is typically 2, 4, 8, or 16.

Because the quadrature is symmetric about $\mu = 0$, the sweep decomposes cleanly
into a **right-going pass** ($\mu_n > 0$) and a **left-going pass** ($\mu_n < 0$).

## Diamond-difference spatial discretization

For cell $i$ with width $\Delta x_i$ and edge fluxes
$\psi_{i-1/2,n},\,\psi_{i+1/2,n}$, the cell-average is

$$
\psi_{i,n} \;=\; \tfrac{1}{2}\bigl(\psi_{i-1/2,n} + \psi_{i+1/2,n}\bigr).
$$

Substituting into the transport equation and re-arranging gives the recurrence for
a right-going ordinate ($\mu_n > 0$):

$$
\psi_{i+1/2,n}
\;=\;
\frac{\bigl(2\mu_n/\Delta x_i - \Sigma_{t,g,i}\bigr)\psi_{i-1/2,n} \;+\; 2\,Q_{i,n}}
     {2\mu_n/\Delta x_i + \Sigma_{t,g,i}},
$$

with the analogous backward recurrence for $\mu_n < 0$. Negative cell-average
fluxes (a known diamond-difference pathology in optically thick cells) are clipped
to zero in post-processing; the sweep itself keeps $\psi_{i-1/2}$ and
$\psi_{i+1/2}$ intact.

## Boundary conditions

| BC | Right-going inlet $\psi(0,\mu_n>0)$ | Left-going inlet $\psi(L,\mu_n<0)$ |
|---|---|---|
| **Vacuum** | 0 | 0 |
| **Reflective** | $\psi(0, -\mu_n)$ | $\psi(L, -\mu_n)$ |
| **White** | $\frac{1}{2}\sum_{\mu_m<0} w_m\lvert\mu_m\rvert \psi(0,\mu_m)$ | symmetric |

## Sweep algorithm

```
repeat until ‖φ_l^(k+1) - φ_l^(k)‖ < tol :
    compute q_l(x) from current φ
    for each μ_n > 0, sweep left → right using the diamond-difference recurrence
    for each μ_n < 0, sweep right → left
    accumulate φ_l(x) = Σ_n w_n P_l(μ_n) ψ_n(x)
```

Convergence is measured on the *scalar* flux $\phi_0$ in a point-wise $L^\infty$
sense. In practice the solver converges in a few dozen source iterations at
UFG resolution, without acceleration; a diffusion synthetic acceleration scheme is a
candidate for future work.

## CLI invocation

```bash
./ufg_1d_app input.i --solver sn --sn-order 4 --legendre 0 \
                    --output-dir ../output/myrun
```

Key flags:

- `--sn-order N` — number of angular ordinates (default 4).
- `--legendre L` — scattering-source expansion order (default 0). Must be ≤ the
  Legendre order used when generating the scattering matrix.
- `--max-iter`, `--tol` — source iteration limits.

## Cross-references

- `docs/solver_theory_and_algorithms.md` — full derivation, including the
  source-moment update and convergence discussion.
- `src/solver/solver_sn.cpp` — implementation.
