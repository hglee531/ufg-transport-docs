# UFG grid, group averaging, and scattering matrices

This page covers the mechanics of going from **point-wise cross sections** to
the **UFG library** that the transport solvers consume: the group structure,
the weighted integration, the Legendre-expanded scattering matrix, and how to
keep the whole thing tractable for tens of thousands of groups.

Two companion pages cover the specialised kernels that write into the same
library:

- [Resonance scattering kernel](resonance_kernel.md) — the O&S Doppler-broadened
  elastic kernel for heavy nuclides.
- [URR p-tables](urr_ptable.md) — narrow-resonance sampling in the unresolved
  region.

## Group structure

UFG-Transport supports two group-structure modes.

### Uniform lethargy

Given a group count $G$, a top energy $E_\text{max}$, and a bottom energy
$E_\text{min}$, the group boundaries are

$$
u_g \;=\; u_0 + g\,\Delta u, \qquad
\Delta u \;=\; \frac{\ln(E_\text{max}/E_\text{min})}{G}, \qquad
E_g \;=\; E_\text{max}\,e^{-g\,\Delta u}.
$$

Uniform lethargy is natural for the slowing-down region, where the $1/E$
weighting spectrum is flat in lethargy.

### Coarse-partitioned

The user declares coarse energy bins and a per-bin fine-group count:

```text
[GroupStructure]
  energy_bin    = '20.0 1.0e-3 0.5e-6 1.0e-11'
  energy_groups = '200 1500 200 100'
[]
```

Each coarse bin is subdivided uniformly in lethargy. This lets the user put
many groups in the slowing-down region, fewer in the fast and thermal tails,
and align boundaries with thermal-scattering cut-offs, fission thresholds, or
upscatter boundaries.

## Group-averaged cross sections

For a reaction $x$ on nuclide $i$, the group-averaged microscopic cross section
over group $g$ (with boundaries $E_g < E_{g-1}$) is

$$
\sigma_{x,i}^{(g)}
\;=\;
\frac{\displaystyle \int_{E_g}^{E_{g-1}} \sigma_{x,i}(E)\,\phi_w(E)\,dE}
     {\displaystyle \int_{E_g}^{E_{g-1}} \phi_w(E)\,dE},
$$

with a user-chosen weighting spectrum $\phi_w(E)$. UFG-Transport defaults to
$\phi_w(E) \propto 1/E$ (flat in lethargy); flat-in-energy is also available
via `use_1_over_E = false`.

Integration is trapezoidal on the **union** of the ACE point-wise grid and the
UFG boundaries: each point-wise interval is split at the containing group
boundary, so no group ever sees a cross section that interpolates *across* its
edge.

Macroscopic region-level cross sections are assembled by the homogenizer,

$$
\Sigma_x^{(g)} \;=\; \sum_k v_k\,\sum_i N_i^{(k)}\,\sigma_{x,i}^{(g)},
$$

where the $v_k$ are region volume fractions and $N_i^{(k)}$ are per-composition
atom densities.

## The scattering transfer matrix

The quantity the solvers actually need is the **Legendre-expanded transfer
matrix**

$$
\Sigma_{s,\ell}^{\,g \to g'} \;=\; \int_{E_{g}}^{E_{g-1}}\!\!\int_{E_{g'}}^{E_{g'-1}}
  \Sigma_{s,\ell}(E \to E')\,\phi_w(E)\,dE\,dE' \;\Big/\;
  \int_{E_{g}}^{E_{g-1}} \phi_w(E)\,dE,
$$

for $\ell = 0,\ldots,L$.

### Legendre expansion

The differential scattering cross section at fixed incident energy $E$ is
expanded in Legendre polynomials of the scattering cosine $\mu$,

$$
\frac{d\sigma_s}{d\mu}(E,\mu)
\;=\;
\sum_{\ell=0}^{L}\,\frac{2\ell+1}{2}\,\sigma_{s,\ell}(E)\,P_\ell(\mu),
\qquad
\sigma_{s,\ell}(E) \;=\; \int_{-1}^{1}\!\frac{d\sigma_s}{d\mu}(E,\mu)\,P_\ell(\mu)\,d\mu.
$$

The Legendre moments $\sigma_{s,\ell}(E)$ are computed directly from the
OpenMC-tabulated $(\mu,\mathrm{PDF})$ distributions; no assumption of
isotropy-in-CM is baked into the grid averaging.

### Elastic two-body kinematics (stationary target)

For two-body elastic scattering on a stationary target of mass ratio $A$, the
outgoing lab energy is a function of the cosine $\mu_\text{cm}$ of the scattering
angle in the centre-of-mass frame,

$$
E' \;=\; E\,\frac{A^2 + 2A\mu_\text{cm} + 1}{(A+1)^2}.
$$

For each incident energy on the union grid, the XS processor quadratures over
$\mu_\text{cm}$; each sample contributes to the transfer matrix element with
the Legendre factor $P_\ell(\mu_\text{lab}(\mu_\text{cm}))$.

This "cold" kinematic kernel is accurate above ~1 keV. Below that, thermal
motion of the target becomes significant (see
[thermal scattering](thermal_scattering.md)), and for heavy nuclides with
strong resonances, the O&S Doppler-broadened kernel replaces it up to
`resonance_cutoff_ev` (default 1 keV); see
[resonance kernel](resonance_kernel.md).

### Discrete-level inelastic

For an inelastic level with Q-value $Q < 0$ (MT 51 – 90), the reaction has a
threshold at $E_\text{th} = |Q|(A+1)/A$. Above threshold, the centre-of-mass
outgoing energy is

$$
E'_\text{cm} \;=\; \frac{A}{A+1}\bigl(E - E_\text{th}\bigr),
$$

which lands in a narrow transfer band in lab frame that touches only a handful
of outgoing groups per incident energy.

### Continuum inelastic and (n, xn)

MT 91 (continuum inelastic) and (n, 2n), (n, 3n), … use the tabulated
secondary-energy distributions from ACE. These contribute to $\ell = 0$ only
in the current implementation (isotropic CM emission); they feed directly into
$\Sigma_{s,0}^{\,g\to g'}$.

### Fission χ

The fission emission spectrum is built by integrating the published secondary-
energy distribution (Maxwell, Watt, or tabular) over each outgoing group,

$$
\chi_g \;=\; \int_{E_g}^{E_{g-1}} \chi(E')\,dE'.
$$

## Storage: dense vs sparse

For $G = 24{,}000$ and $\ell = 0$, a dense matrix $\Sigma_{s,\ell}^{\,g \to g'}$
is $24{,}000 \times 24{,}000$ doubles — 4.6 GB. For typical UFG configurations
(thousands of groups, $L = 0$, highly diagonal thermal and slowing-down
transfers), UFG-Transport switches to a **sparse banded-row** storage once
$G \geq 5000$:

- each row $g$ stores only the non-zero $g'$ entries above a configurable
  threshold (`scatter_threshold = 1e-6` by default);
- thermal up-scatter bands are retained exactly — the thresholding only trims
  numerical tails;
- the 1-D solvers always see a **dense `RegionXS`**; the sparse matrix is
  densified at the solver boundary, so the sparse/dense choice is invisible to
  the transport loop.

The recommended UFG configuration is $P_0$ scattering, threshold $10^{-6}$,
and $G \approx 24{,}000$.

## References inside the repo

- `src/xs/xs_processor.cpp` — group-averaging and Legendre-moment loops.
- `src/xs/homogenizer.cpp` — multi-nuclide macroscopic assembly.
- `src/xs/scatter_matrix.cpp` — sparse banded storage.
- `include/xs/ufg_grid.hpp` — group-structure construction.
