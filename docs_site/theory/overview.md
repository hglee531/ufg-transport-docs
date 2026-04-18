# Pipeline overview

UFG-Transport takes **point-wise (ACE) nuclear data** and turns it into a **solved
transport equation** by passing through an *ultra-fine-group* (UFG) representation
of the cross sections — instead of a broad-group multigroup library with heavy
self-shielding corrections.

```
   ACE / OpenMC HDF5       UFG library                Transport solution
  ┌───────────────────┐   ┌────────────────────┐    ┌──────────────────────┐
  │  σ(E) pointwise,  │   │  Σ_x,g, Σ_s,g→g',  │    │  φ(r, g) from        │
  │  angular dists,   │ → │  χ_g, Legendre     │ →  │  slab Sn / slab Pn;  │
  │  thermal S(α,β),  │   │  moments ℓ = 0..L  │    │  optional HFG → UFG  │
  │  URR p-tables     │   │  + resonance kernel│    │  → MG collapse       │
  └───────────────────┘   └────────────────────┘    └──────────────────────┘
```

## Why ultra-fine groups?

Broad-group multigroup libraries (tens to hundreds of groups) rely on heavy
**resonance self-shielding** approximations — Bondarenko tables, equivalence
theory, subgroup methods — to preserve effective reaction rates inside
resonances. At the cost of memory and compute, UFG libraries bypass most of that
machinery:

- Several thousand groups narrow the within-group spectral variation so that a
  simple $1/E$ or flat-in-energy weighting is adequate.
- Resonance peaks in the **resolved resonance region** are resolved directly on
  the group grid, so self-shielding is largely automatic.
- Group-to-group transfer matrices — including thermal up-scatter — become the
  primary data structure; everything else is read-only.

For the **unresolved resonance region (URR)** and for **Doppler-broadened elastic
scattering** at reactor temperatures, the grid alone is not enough: the
group-averaged cross section still depends on the (unknown) local flux through
statistical resonances, and the simple stationary-target kinematic kernel
over-estimates up-scatter in resonances. UFG-Transport adds two targeted models
on top of the UFG machinery:

- **Narrow-resonance p-table sampling** in the URR
  ([URR p-tables](urr_ptable.md)).
- **Ouisloumen–Sanchez (O&S) resonance scattering kernel** for heavy nuclides
  ([Resonance kernel](resonance_kernel.md)).

Both models write into the same UFG library; the transport solver never sees
them directly.

## The three stages

### 1. ACE parsing

The [ACE reader](ace_and_hdf5.md) opens OpenMC-style HDF5 cross-section files
and extracts, for each nuclide and each tabulated temperature:

- the unionised energy grid and point-wise reaction cross sections
  (total, elastic, inelastic levels, fission, capture, …);
- Legendre-coefficient or tabular angular distributions;
- fission $\chi$ distributions (Maxwell, Watt, tabular);
- unresolved-resonance probability tables (ENDF MF 2, LRU = 2);
- thermal $S(\alpha,\beta)$ kernels (incoherent / coherent).

### 2. UFG cross-section processing

The [UFG grid builder](ufg_grid.md) defines group boundaries by uniform lethargy
or by a coarse-partitioned lethargy structure. The XS processor then produces,
for each nuclide at each tabulated temperature, a complete UFG library:

- principal reactions $\sigma_{t,g},\,\sigma_{a,g},\,\sigma_{f,g},\,\nu\sigma_{f,g}$;
- Legendre-expanded scattering moments $\sigma_{s,\ell}^{g\to g'}$ from two-body
  cold-elastic kinematics, discrete-level inelastic, and continuum inelastic;
- thermal scattering from the [three-regime model](thermal_scattering.md): cold
  elastic → free-gas → $S(\alpha,\beta)$;
- optional Doppler-broadened resonance elastic via the
  [O&S kernel](resonance_kernel.md);
- URR-sampled principal cross sections via
  [NR p-tables](urr_ptable.md);
- fission $\chi_g$ by integrating the published secondary-energy distribution.

Cross-temperature interpolation is linear in $\tau = T^{1/3}$; see
[multi-level condensation](multilevel.md#temperature-interpolation).

Mixtures are produced by the homogenizer,

$$
\Sigma_x^{(g)} \;=\; \sum_k v_k \sum_i N_i^{(k)} \,\sigma_{x,i}^{(g)},
$$

where $k$ indexes compositions with volume fractions $v_k$ and $N_i^{(k)}$ are
atom densities.

### 3. Transport solve

UFG-Transport currently supports two **verified** 1-D deterministic solvers:

- **[slab $S_N$](solvers_sn.md)** — diamond-difference spatial discretisation,
  Gauss–Legendre angular quadrature, anisotropic scattering to arbitrary
  Legendre order;
- **[slab $P_N$](solvers_pn.md)** — cell-centred finite difference with block
  Thomas, $P_1$ through $P_7$, Marshak vacuum boundaries.

For problems where even the UFG solve is too expensive to run everywhere, the
[multi-level pipeline](multilevel.md) first solves a **hyperfine-group (HFG)**
problem in 0-D per composition, uses the resulting flux to collapse to a UFG
library, and only then runs 1-D transport. A further UFG → MG collapse produces
a broad-group library (e.g. CASMO-70) for downstream consumers.

An orthogonal research track, [Spectrum Expansion (SpecEx)](../specex/index.md), is
developing a reduced-order replacement for the entire resonance self-shielding
stack.

## Data formats

- **Input**: [HIT](../guide/running.md) (hierarchical input text) — nested
  `[block] … []` sections of `key = value` pairs.
- **Cross sections**: OpenMC `cross_sections.xml` + per-nuclide HDF5 files.
- **Output**: CSV files written under `output/<run-name>/`, plus logged
  quantities (region-averaged flux, Legendre moments, condensed MG XS).

## What this code is not (yet)

- Not a $k$-effective eigenvalue code. Only fixed-source problems are solved;
  power iteration is on the roadmap.
- Not a 2-D / 3-D code. Geometry is infinite, slab, or (for 0-D and for the
  SpecEx sweep) homogeneous mixtures.
- Not a production lattice code. It is a research implementation focused on
  making the UFG-plus-resonance-kernel workflow first-class, and on exploring
  the SpecEx reduced-order idea.

See the [status section on the home page](../index.md#status) for the concrete
list of what is implemented, what is disabled, and what is in flight.
