# Pipeline overview

UFG-Transport is organised around a single idea: go from **point-wise nuclear data**
to a **solved transport equation** by passing through an *ultra-fine-group* (UFG)
representation, rather than a traditional broad-group library.

```
   ACE / OpenMC HDF5       UFG library               Transport solution
  ┌───────────────────┐   ┌──────────────────┐     ┌──────────────────────┐
  │  σ(E) pointwise,  │ → │  Σ_x,g, Σ_s,g→g′,│  →  │  φ(r, g) from Sn/Pn/ │
  │  angular dists,   │   │  χ_g, Legendre   │     │  CPM; optional HFG → │
  │  thermal S(α,β)   │   │  moments L=0..N  │     │  UFG → MG collapse   │
  └───────────────────┘   └──────────────────┘     └──────────────────────┘
```

## Why ultra-fine groups?

Broad-group multigroup libraries (tens to hundreds of groups) rely on heavy
**resonance self-shielding** approximations — Bondarenko tables, equivalence theory,
subgroup methods — to preserve effective reaction rates inside resonances. At the
cost of memory and compute, UFG libraries side-step most of that machinery:

- Several thousand groups reduce the *within-group* spectrum variation enough that a
  simple 1/*E* or flat weighting spectrum is adequate.
- Resonance peaks are resolved directly on the group grid, so self-shielding is
  largely automatic.
- Group-to-group transfer matrices, including thermal up-scattering, become the
  primary data structure — everything else is read-only.

This trades pre-processing cost for modelling cleanliness: the UFG library is
*closer* to the point-wise data, and the deterministic solver that uses it can
remain straightforward.

## The three stages

### 1. ACE parsing

The [ACE reader](ace_and_hdf5.md) opens OpenMC-style HDF5 cross-section files and
extracts, for each nuclide and each ACE temperature:

- union energy grid and pointwise reaction cross sections (total, elastic,
  inelastic levels, fission, capture, …)
- Legendre-coefficient or tabular angular distributions
- fission χ distributions (Maxwell, Watt, tabular)
- thermal S(α,β) kernels (incoherent / coherent)

### 2. UFG group averaging

The [UFG grid builder](ufg_grid.md) defines group boundaries either by uniform
lethargy or by a coarse-partitioned lethargy structure. The XS processor integrates
the point-wise data against a weighting spectrum (default 1/*E*) on each group:

$$
\Sigma_{x,g} \;=\; \frac{\int_{E_g}^{E_{g-1}} \sigma_x(E)\,\phi_w(E)\,dE}
                         {\int_{E_g}^{E_{g-1}} \phi_w(E)\,dE}
$$

Scattering transfer matrices are built with two-body lab-frame kinematics for elastic
scattering and with threshold kinematics for discrete inelastic levels. Thermal
scattering is handled by the [three-regime model](thermal_scattering.md).

Mixtures are produced by the **homogenizer**, which assembles the macroscopic
cross sections of a multi-nuclide region:

$$
\Sigma_x^{(g)} \;=\; \sum_k v_k \sum_i N_i^{(k)} \,\sigma_{x,i}^{(g)}.
$$

### 3. Transport solve

Three solvers share the same region-wise macroscopic library:

- **[S<sub>N</sub> (slab)](solvers_sn.md)** — diamond-difference with Gauss-Legendre
  angular quadrature, anisotropic scattering via Legendre expansion.
- **[P<sub>N</sub> (slab)](solvers_pn.md)** — cell-centred finite difference, block
  Thomas solve, Marshak vacuum boundary conditions.
- **[CPM (cylindrical)](solvers_cpm.md)** — collision probabilities evaluated via
  Bickley–Naylor functions on ring chords.

For problems where even the UFG solve is too expensive everywhere, the
[multi-level pipeline](multilevel.md) solves an even finer HFG library in 0-D per
composition, uses the resulting flux to collapse to a UFG library, and only the UFG
library is carried into the 1-D transport solve. A further UFG → MG collapse produces
a standard broad-group library (e.g. CASMO-70) for downstream use.

## Data formats

- **Input**: [HIT](../guide/running.md) (hierarchical-input-text) files — blocks of
  `key = value` pairs in nested `[block] ... []` sections.
- **Cross sections**: OpenMC `cross_sections.xml` + per-nuclide HDF5 files.
- **Output**: CSV files written per run under `output/<run-name>/`, plus logged
  quantities (region-averaged flux, Legendre moments, condensed MG XS).

## What this code is not

- Not a k-effective eigenvalue code yet — only fixed-source problems are solved.
  Power iteration is on the roadmap.
- Not a 2-D/3-D code. The geometry module supports infinite, slab, and 1-D
  cylindrical ("pin-cell") geometries only.
- Not a production lattice code. It is a research implementation, optimised for
  clarity and extensibility of the UFG-centric workflow.
