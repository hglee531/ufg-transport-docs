# Sn — slab geometry

The slab $S_N$ solver in UFG-Transport discretises the 1-D plane-geometry
transport equation in angle using Gauss–Legendre quadrature and in space
using **diamond difference**. It supports anisotropic scattering to arbitrary
Legendre order and three boundary conditions: reflective, vacuum, and white
(isotropic return).

This is one of the two currently **verified** transport solvers (see
[status](../index.md#status)). The matching $P_N$ solver is documented on
the [slab $P_N$ page](solvers_pn.md).

## Governing equation

Within energy group $g$, the angular flux $\psi_g(x,\mu)$ along direction
cosine $\mu$ satisfies the 1-D plane-geometry transport equation

$$
\mu\,\frac{\partial \psi_g}{\partial x}(x,\mu)
\;+\; \Sigma_{t,g}(x)\,\psi_g(x,\mu)
\;=\;
\sum_{\ell=0}^{L}\,\frac{2\ell+1}{2}\,P_\ell(\mu)\,q_{\ell,g}(x),
$$

with the Legendre-moment source

$$
q_{\ell,g}(x)
\;=\;
\sum_{g'} \Sigma_{s,\ell}^{\,g'\to g}(x)\,\phi_{\ell,g'}(x)
\;+\; \delta_{\ell 0}\,
      \bigl[\chi_g\,Q_\text{fis}(x) + Q_{\text{ext},g}(x)\bigr].
$$

The scalar flux moments are recovered by the angular quadrature,

$$
\phi_{\ell,g}(x) \;=\; \int_{-1}^{1} P_\ell(\mu)\,\psi_g(x,\mu)\,d\mu
\;\approx\; \sum_{n=1}^{N} w_n\,P_\ell(\mu_n)\,\psi_{n,g}(x).
$$

The fission source $Q_\text{fis}$ is fixed externally (UFG-Transport solves
**fixed-source** problems; $k$-effective eigenvalue iteration is on the
roadmap).

## Angular quadrature

UFG-Transport uses **symmetric Gauss–Legendre** quadrature: the $\mu_n$ are
the roots of $P_N(\mu)$ on $[-1, 1]$, and the weights $w_n$ are the
matching Gauss weights. Nodes and weights are computed at run time by the
Golub–Welsch algorithm in `src/utils/quadrature.cpp`. The quadrature order
$N$ is selectable via `--sn-order` (typically 2, 4, 8, 16).

Because Gauss–Legendre is symmetric about $\mu = 0$, the sweep decomposes
cleanly into a **right-going pass** ($\mu_n > 0$) and a **left-going pass**
($\mu_n < 0$).

## Diamond-difference spatial discretisation

For cell $i$ of width $\Delta x_i$ with edge fluxes $\psi_{i-1/2,n}$ and
$\psi_{i+1/2,n}$, the diamond-difference closure sets the cell-average equal
to the arithmetic mean of the two edges,

$$
\psi_{i,n} \;=\; \tfrac{1}{2}\bigl(\psi_{i-1/2,n} + \psi_{i+1/2,n}\bigr).
$$

Substituting into the transport equation and solving for the downwind edge
gives the recurrence for a right-going ordinate ($\mu_n > 0$):

$$
\psi_{i+1/2,n}
\;=\;
\frac{\bigl(2\mu_n/\Delta x_i - \Sigma_{t,g,i}\bigr)\,\psi_{i-1/2,n}
      \;+\; 2\,Q_{i,n}}
     {2\mu_n/\Delta x_i + \Sigma_{t,g,i}},
$$

with the analogous backward recurrence for $\mu_n < 0$. $Q_{i,n}$ is the
directional source evaluated at the cell centre. Diamond difference is known
to produce **negative cell-average fluxes** in optically thick cells; these
are clipped to zero in post-processing, but the sweep itself propagates the
edge fluxes $\psi_{i\pm 1/2}$ intact.

## Boundary conditions

| BC | Right-going inlet, $\psi(0,\mu_n{>}0)$ | Left-going inlet, $\psi(L,\mu_n{<}0)$ |
|---|---|---|
| **Vacuum** | 0 | 0 |
| **Reflective** | $\psi(0,\,-\mu_n)$ | $\psi(L,\,-\mu_n)$ |
| **White** | $\tfrac{1}{\pi}\sum_{\mu_m<0} w_m\,\lvert\mu_m\rvert\,\psi(0,\mu_m)$ | symmetric |

(The white condition returns the total leakage as an isotropic inward flux,
normalised on the half-range.)

## Sweep algorithm

```text
repeat until max_g || φ_0(g)^{k+1} - φ_0(g)^{k} ||_∞ / || φ_0(g)^{k} ||_∞ < tol :
    for each group g (in upscatter order):
        compute moment source q_ℓ(x) from current φ_ℓ
        for each μ_n > 0, sweep left → right (diamond-difference recurrence)
        for each μ_n < 0, sweep right → left
        accumulate  φ_ℓ(x) = Σ_n w_n P_ℓ(μ_n) ψ_n(x)
    if upscatter detected, re-sweep thermal block until inner tolerance
```

Convergence is measured on the **scalar flux** $\phi_0$ in a per-group
relative $L^\infty$ sense. A global criterion masks the thermal groups,
which is why the code uses a per-group metric and a **thermal inner loop**
(`thermal_max_inners`, `thermal_tolerance`) for the up-scatter block — the
fast groups converge in a handful of outer iterations, but the thermal
block needs many sweeps because the source-iteration spectral radius for
up-scatter is near unity.

A diffusion-synthetic acceleration (DSA) scheme is a candidate for future
work if source iteration becomes a bottleneck.

## CLI invocation

```bash
./ufg_1d_app input.i --solver sn --sn-order 4 --legendre 0 \
                     --output-dir ../output/myrun
```

Key flags:

- `--sn-order N` — number of angular ordinates (default 4).
- `--legendre L` — scattering-source expansion order (default 0); must be
  ≤ the Legendre order used when generating the scattering matrix.
- `--max-iters`, `--tolerance` — outer source-iteration limits.
- `--threads N` — OpenMP threads (default 4).

## Files

- `src/solver/solver_sn.cpp`, `include/solver/solver_sn.hpp` — implementation.
- `src/utils/quadrature.cpp` — Golub–Welsch Gauss–Legendre weights.
