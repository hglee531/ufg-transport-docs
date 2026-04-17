# UFG grid, group averaging, and scattering matrices

This page covers the mechanics of going from **point-wise cross sections** to the
**UFG library** that the transport solvers consume: the group structure, the weighted
integration, the Legendre-expanded scattering matrix, and how to keep the whole thing
tractable for tens of thousands of groups.

See `REPORT_ANGULAR_MOMENTS_AND_SCATTERING_MATRIX.md` for the full derivation of how
angular distributions are stored in OpenMC HDF5 and turned into Legendre moments.

## Group structure

UFG-Transport supports two group-structure modes.

### Uniform lethargy

Given `num_groups` $G$, top energy $E_\text{max}$, and bottom energy $E_\text{min}$,
the group boundaries are

$$
u_g \;=\; u_0 + g\,\Delta u, \qquad
\Delta u \;=\; \frac{\ln(E_\text{max}/E_\text{min})}{G},
\qquad E_g \;=\; E_\text{max}\,e^{-g\,\Delta u}.
$$

Lethargy spacing is natural for the slowing-down region, where the 1/*E* weighting
spectrum is itself uniform in lethargy.

### Coarse-partitioned

The user declares coarse energy bins and a per-bin fine-group count:

```
[GroupStructure]
  energy_bin    = '20.0 1.0e-3 0.5e-6 1.0e-11'
  energy_groups = '200 1500 200 100'
[]
```

Each coarse bin is subdivided uniformly in lethargy. This gives the flexibility to
put many groups in the slowing-down region, fewer in the fast and thermal tails, and
to align boundaries with thermal-scattering cut-offs.

## Group-averaged cross sections

For a reaction $x$ on nuclide $i$, the group-averaged microscopic cross section is

$$
\sigma_{x,i}^{(g)} \;=\; \frac{\displaystyle\int_{E_g}^{E_{g-1}} \sigma_{x,i}(E)\,\phi_w(E)\,dE}
                              {\displaystyle\int_{E_g}^{E_{g-1}} \phi_w(E)\,dE}
$$

with a user-chosen weighting spectrum $\phi_w(E)$. UFG-Transport defaults to
$\phi_w \propto 1/E$ (flat in lethargy) but also supports a flat-in-energy weighting.

Integration is trapezoidal on the **union** of the ACE point-wise grid and the UFG
boundaries: each pointwise interval is split at the containing group boundary, so no
group ever sees a cross section that interpolates *across* its edge.

Macroscopic cross sections for a homogeneous region $k$ with atom densities
$N_i^{(k)}$ are assembled by the homogenizer:

$$
\Sigma_x^{(k,g)} \;=\; \sum_i N_i^{(k)}\,\sigma_{x,i}^{(g)}.
$$

Material-to-region homogenisation with volume fractions $v_k$ gives

$$
\Sigma_x^{(g)} \;=\; \sum_k v_k\,\Sigma_x^{(k,g)}.
$$

## The scattering transfer matrix

The quantity the solvers actually need is
$\Sigma_{s,\ell}^{(g\to g')}$: the cross section for a neutron in group $g$ to
scatter into group $g'$ with Legendre anisotropy order $\ell$.

### Legendre expansion

For a given incident energy $E$, the differential scattering cross section is
expanded in Legendre polynomials in the scattering cosine $\mu$:

$$
\sigma_s(E,\mu) \;=\; \sum_{\ell=0}^{L}\,
   \frac{2\ell + 1}{2}\,\sigma_{s,\ell}(E)\,P_\ell(\mu),
\qquad
\sigma_{s,\ell}(E) \;=\; \int_{-1}^{1} \sigma_s(E,\mu)\,P_\ell(\mu)\,d\mu.
$$

The Legendre moments $\sigma_{s,\ell}(E)$ are computed directly from the OpenMC
tabulated $(\mu,\,\text{PDF})$ distributions.

### Elastic two-body kinematics

For two-body elastic scattering on a stationary target of mass ratio $A$, the
outgoing lab energy is
$$
E' \;=\; E\,\frac{A^2 + 2A\mu_\text{cm} + 1}{(A+1)^2}.
$$
For each incident energy on the union grid, UFG-Transport quadratures over
$\mu_\text{cm}$; each sample contributes to a specific group-to-group transfer.

### Discrete-level inelastic

For MT 51–90 (discrete inelastic levels), the outgoing CM energy at threshold $Q<0$ is
$$
E'_\text{cm} \;=\; \frac{A}{A+1}\bigl(E - |Q|\bigr),
$$
leading to a narrow transfer band that lands in a small number of groups per incident
energy. MT 91 (continuum inelastic) uses the stored continuous energy–angle
distribution.

### Fission

Fission χ is built by integrating the published secondary-energy distribution
(Maxwell, Watt, or tabular) over each outgoing group:

$$
\chi_g \;=\; \int_{E_g}^{E_{g-1}} \chi(E')\,dE'.
$$

## Storage: dense vs sparse

For $G = 24{,}000$ and $L = 0$, a dense matrix $\Sigma_{s,\ell}^{(g\to g')}$ is
$24{,}000 \times 24{,}000$ doubles — 4.6 GB. For the defaults the project uses
(thousands of groups, $L = 0$, highly diagonal thermal and slowing-down transfers),
UFG-Transport switches to a **sparse banded-row** storage when $G \geq 5000$:

- each row $g$ stores only the non-zero $g'$ entries above a configurable threshold
  (`scatter_threshold = 1e-6` by default);
- thermal up-scatter bands near $E \lesssim E_\text{thermal}$ are retained exactly;
- 1-D solvers receive **dense** `RegionXS` — the sparse matrix is densified at the
  solver boundary, so the sparse/dense choice is invisible to the transport loop.

The recommended UFG configuration is $P_0$ scattering, threshold $10^{-6}$, and
$G \approx 24{,}000$.

## References within the repo

- `src/xs/xs_processor.cpp` — the group-averaging and Legendre-moment loops.
- `src/xs/homogenizer.cpp` — multi-nuclide macroscopic assembly.
- `src/xs/scatter_matrix.cpp` — sparse banded storage.
- `include/xs/ufg_grid.hpp` — group-structure construction.
- `REPORT_ANGULAR_MOMENTS_AND_SCATTERING_MATRIX.md` — full discussion of how OpenMC
  stores angular data and how the Legendre integrals are evaluated.
