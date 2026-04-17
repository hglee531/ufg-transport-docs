# Pn — slab geometry

The $P_N$ solver expands the angular flux in Legendre polynomials rather than
discrete ordinates, yielding a **coupled system of $N+1$ moment equations** per
energy group. UFG-Transport implements $P_1$, $P_3$, $P_5$, and $P_7$ in slab
geometry with cell-centred finite differences and a block Thomas direct solve.

The complete derivation is in `docs/solver_theory_and_algorithms.md`. This page
gives the working equations and the algorithm.

## Moment equations

Expanding
$\psi_g(x,\mu) = \sum_\ell \frac{2\ell+1}{2}\,\phi_{\ell,g}(x)\,P_\ell(\mu)$
and projecting the transport equation onto each $P_\ell$ yields, for $\ell=0,\ldots,N$,

$$
\frac{\ell}{2\ell+1}\,\frac{d\phi_{\ell-1,g}}{dx}
\;+\;
\frac{\ell+1}{2\ell+1}\,\frac{d\phi_{\ell+1,g}}{dx}
\;+\;
\Sigma_{r,\ell,g}(x)\,\phi_{\ell,g}(x)
\;=\;
q_{\ell,g}(x),
$$

with the **removal** cross section

$$
\Sigma_{r,\ell,g}(x) \;=\; \Sigma_{t,g}(x) \;-\; \Sigma_{s,\ell}^{g\to g}(x).
$$

Closure is imposed by setting $\phi_{N+1,g} \equiv 0$ (Mark closure).

## Cell-centred finite difference

All moments live at cell centres. First derivatives use centred differences:

$$
\frac{d\phi_{\ell}}{dx}\biggr|_{x_i}
\;\approx\;
\frac{\phi_{\ell,i+1} - \phi_{\ell,i-1}}{x_{i+1} - x_{i-1}}.
$$

Stacking the $N+1$ moment equations at each cell gives a **block-tridiagonal
linear system**: each diagonal block is an $(N+1)\times(N+1)$ matrix acting on the
cell's moment vector.

## Block Thomas solve

The block-tridiagonal system is solved by a standard block Thomas (Gaussian
elimination with block pivots) implemented in `src/solver/solver_pn.cpp`. Because
scattering couples all moments through $q_{\ell,g}$, each Thomas solve is repeated
as a **source iteration** at the group level, just as for the $S_N$ solver.

## Boundary conditions

UFG-Transport supports **Marshak vacuum**, **reflective**, and **white** boundaries.
Marshak imposes half-space integrals

$$
\int_0^1 \psi_g(0,\mu)\,P_\ell(\mu)\,d\mu \;=\; 0
\qquad \text{for odd } \ell,
$$

which for $P_3$ gives two conditions per boundary (odd moments $\ell=1,3$). These
are applied at the boundary cells by modifying the first and last block rows of the
Thomas system.

## CLI invocation

```bash
./ufg_1d_app input.i --solver pn --pn-order 3 \
                    --output-dir ../output/myrun
```

`--pn-order` accepts 1, 3, 5, or 7. Higher orders improve angular accuracy but
quadruple the cost per Thomas sweep (block size $(N+1)^2$).

## When to pick Pn vs Sn

- $P_N$ handles strongly anisotropic scattering with fewer unknowns than $S_N$ of
  comparable accuracy, and has no ray effect.
- $S_N$ handles heterogeneous streaming and reflecting boundaries more naturally
  and is usually the default in slabs.
- In the UFG-Transport demonstrations the two agree to within a few percent on
  scalar flux for the PWR pin-cell slab problems.

## Cross-references

- `docs/solver_theory_and_algorithms.md` §6 — full $P_N$ derivation with boundary
  closure details.
- `src/solver/solver_pn.cpp` — implementation.
