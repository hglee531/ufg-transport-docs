# Pn — slab geometry

The slab $P_N$ solver expands the angular flux in Legendre polynomials
rather than discrete ordinates, producing a **coupled system of $N+1$
moment equations per energy group**. UFG-Transport implements $P_1$, $P_3$,
$P_5$, and $P_7$ in 1-D slab geometry with cell-centred finite differences
and a block Thomas direct solve.

$P_N$ is the second of the two currently verified transport solvers; the
matching discrete-ordinates solver is [slab $S_N$](solvers_sn.md).

## Moment equations

Expanding the angular flux as

$$
\psi_g(x,\mu) \;=\; \sum_{\ell=0}^{N} \frac{2\ell+1}{2}\,
  \phi_{\ell,g}(x)\,P_\ell(\mu),
$$

and projecting the transport equation onto each $P_\ell$ over $\mu \in [-1,1]$
yields, for $\ell = 0, \ldots, N$,

$$
\frac{\ell}{2\ell+1}\,\frac{d\phi_{\ell-1,g}}{dx}(x)
\;+\;
\frac{\ell+1}{2\ell+1}\,\frac{d\phi_{\ell+1,g}}{dx}(x)
\;+\;
\Sigma_{r,\ell,g}(x)\,\phi_{\ell,g}(x)
\;=\;
q_{\ell,g}(x),
$$

with the **removal cross section**

$$
\Sigma_{r,\ell,g}(x) \;=\;
\Sigma_{t,g}(x) \;-\; \Sigma_{s,\ell}^{\,g\to g}(x),
$$

and the moment source

$$
q_{\ell,g}(x) \;=\;
\sum_{g'\neq g} \Sigma_{s,\ell}^{\,g'\to g}(x)\,\phi_{\ell,g'}(x)
\;+\; \delta_{\ell 0}\,
\bigl[\chi_g\,Q_\text{fis}(x) + Q_{\text{ext},g}(x)\bigr].
$$

**Mark closure** is imposed by setting $\phi_{N+1,g} \equiv 0$, which makes
the coupled moment system finite.

## Cell-centred finite difference

All moments live at cell centres. First derivatives use centred differences
across neighbouring cell centres,

$$
\left.\frac{d\phi_{\ell}}{dx}\right|_{x_i}
\;\approx\;
\frac{\phi_{\ell,i+1} - \phi_{\ell,i-1}}{x_{i+1} - x_{i-1}}.
$$

Stacking the $N+1$ moment equations at each cell gives a
**block-tridiagonal linear system**: each diagonal block is an
$(N+1)\times(N+1)$ matrix acting on the cell's moment vector, coupled to its
neighbours by the derivative blocks.

## Block Thomas solve

The block-tridiagonal system is solved by a standard block Thomas algorithm
(Gaussian elimination with $(N+1)\times(N+1)$ block pivots) in
`src/solver/solver_pn.cpp`. Because scattering couples all moments through
$q_{\ell,g}$, each Thomas solve is repeated as a **source iteration** at the
group level, exactly as for $S_N$. The same thermal-inner-loop treatment
used in $S_N$ applies here: the up-scatter block is re-swept per outer
iteration until it stabilises.

## Boundary conditions

UFG-Transport supports **Marshak vacuum**, **reflective**, and **white**
conditions.

Marshak vacuum imposes half-space integrals

$$
\int_0^1 \psi_g(0,\mu)\,P_\ell(\mu)\,d\mu \;=\; 0
\qquad \text{for odd } \ell,
$$

which for $P_3$ gives two conditions per boundary (odd moments
$\ell = 1, 3$). These are applied at the boundary cells by modifying the
first and last block rows of the Thomas system.

Reflective imposes $\phi_\ell = 0$ for odd $\ell$ at the reflecting edge.
White is implemented by equating the outgoing half-range current with an
isotropic inward flux.

## CLI invocation

```bash
./ufg_1d_app input.i --solver pn --pn-order 3 \
                     --output-dir ../output/myrun
```

Key flags:

- `--pn-order N` — accepts 1, 3, 5, 7 (default 3). Higher orders improve
  angular accuracy but scale the Thomas block cost as $(N+1)^3$.
- `--legendre L` — scattering-source expansion order (default 0).
- `--max-iters`, `--tolerance` — outer source-iteration limits.

## When to pick Pn vs Sn

- $P_N$ has **no ray effect** and handles strongly anisotropic scattering
  with fewer unknowns than $S_N$ of comparable accuracy.
- $S_N$ handles reflective boundaries and heterogeneous streaming more
  naturally and is the usual default in slab geometry.
- In the PWR pin-cell demonstrations the two solvers agree to within a few
  percent on the scalar flux.

## Files

- `src/solver/solver_pn.cpp`, `include/solver/solver_pn.hpp` — implementation.
