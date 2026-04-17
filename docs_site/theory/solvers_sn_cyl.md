# Sn — cylindrical geometry

In cylindrical geometry the $S_N$ streaming operator picks up an **angular
redistribution term** due to the curvature of the coordinate system: as a neutron
flies across a pin, its angle relative to the local $\hat r$ changes even though
its direction in space does not. This is the defining complication of 1-D
cylindrical discrete ordinates.

UFG-Transport implements cylindrical $S_N$ with the classical
**Carlson angular-diamond-difference** treatment and an $\alpha$ recurrence for
the redistribution term. The full derivation (coordinate system, $\alpha$
construction, Carlson absolute-value formulation, sweep order) is in
`docs/Sn_1D_Cylindrical.md`. This page summarises the essentials.

## Governing equation

For azimuthal symmetry and no axial dependence, with direction cosines
$\xi = \cos\phi_\Omega$ and $\eta = \sin\phi_\Omega$ in the transverse plane,
the steady-state transport equation reads

$$
\frac{\xi}{r}\,\frac{\partial}{\partial r}\bigl(r\,\psi\bigr)
\;-\;
\frac{1}{r}\,\frac{\partial}{\partial \phi_\Omega}\bigl(\eta\,\psi\bigr)
\;+\;
\Sigma_t\,\psi
\;=\;
q,
$$

with $q = \sum_\ell \frac{2\ell+1}{2}P_\ell(\mu)q_\ell$ the usual angular-moment
source. The second term on the left is the redistribution: it does not appear in
slab geometry and couples the discrete directions in the sweep.

## Angular discretization: the α recurrence

Discrete angles are chosen as a symmetric Gauss–Legendre set $\{\xi_n, w_n\}$ and
the redistribution derivative $\partial(\eta\psi)/\partial\phi_\Omega$ is
discretised by **diamond difference in angle**:

$$
\frac{1}{r}\,\frac{\partial(\eta\,\psi)}{\partial\phi_\Omega}
\;\approx\;
\frac{\alpha_{n+1/2}\,\psi_{n+1/2} \;-\; \alpha_{n-1/2}\,\psi_{n-1/2}}{r\,w_n}.
$$

The $\alpha_{n\pm1/2}$ are determined by the recurrence

$$
\alpha_{n+1/2} \;=\; \alpha_{n-1/2} \;+\; \xi_n\,w_n,\qquad
\alpha_{1/2} \;=\; 0,
$$

which guarantees that a spatially uniform, isotropic flux gives zero net
redistribution. The angular "half-flux" $\psi_{n\pm1/2}$ is tied to the cell
flux $\psi_n$ by the angular diamond-difference relation
$\psi_n = \tfrac12(\psi_{n-1/2} + \psi_{n+1/2})$.

## Carlson absolute-value formulation

Inside rings close to $r = 0$ the $1/r$ coefficient becomes singular, and the
naive diamond-difference sweep can produce negative angular fluxes. UFG-Transport
adopts Carlson's **absolute-value formulation**, which rewrites the radial
streaming operator so that each ring sees a non-negative contribution from its
neighbours regardless of $\xi_n$'s sign. This is the standard remedy in
cylindrical $S_N$ and is implemented in `src/solver/solver_sn_cyl.cpp`.

## Sweep order and boundary conditions

The sweep proceeds:

1. **Starting direction** ($\xi_n$ most negative) — inward-pointing ray at the
   outer ring, with outer boundary condition imposing $\psi_{N_r+1/2} = $ BC.
2. **Inward radial sweep** to $r = 0$ with the angular-update substituting
   $\psi_{n+1/2}$ from the previous angular index.
3. **Angular reflection** at $\xi_n = 0$: the inward rays rotate to outward rays
   through the symmetric ($r=0$) boundary.
4. **Outward radial sweep** back to the outer ring.

Supported boundary conditions:

- **Inner** ($r=0$) — reflective by symmetry (implicit).
- **Outer** — reflective, vacuum, or white. White is the typical Wigner–Seitz BC.

## CLI invocation

```bash
./ufg_1d_app pwr_pincell_multilevel_sn.i --solver sn --sn-order 8 \
                                        --output-dir ../output/pincell_sn
```

Orders $S_2, S_4, S_8, S_{16}$ are all validated; higher orders are permitted but
diminishing returns set in past $S_{16}$ in the pin-cell demonstrations.

## Implementation notes

The cylindrical $S_N$ solver went through several iterations in the project
history to resolve numerical instabilities near $r=0$ and at ring interfaces.
Relevant design notes live in:

- `docs/sn_cylindrical_formulation_summary.md` — compact summary and the final
  algorithm.
- `docs/report_cylindrical_sn_investigation.md`,
  `docs/report_cylindrical_sn_diagnosis.md` — instability diagnosis and fix.
- `docs/solver_sn_cyl_code_flow.md` — walkthrough of the implementation.

## Cross-references

- `docs/Sn_1D_Cylindrical.md` — full derivation.
- `src/solver/solver_sn_cyl.cpp` — implementation.
