# UFG-Transport

**Ultra-fine-group deterministic neutron transport with an ACE → UFG → solver
pipeline.**

UFG-Transport is a research-oriented neutron transport code that processes
point-wise ACE cross-section data directly into ultra-fine multi-group
libraries (up to tens of thousands of groups), then solves a 1-D deterministic
transport problem. It is written in C++17 with a staged CUDA backend.

## What it does

- Reads point-wise cross sections from **OpenMC-style HDF5 ACE libraries**
  (ENDF/B-VIII).
- Builds **ultra-fine group (UFG)** libraries with uniform-lethargy or
  coarse-partitioned group structures — thousands of groups routinely, tens
  of thousands if desired.
- Constructs **Legendre-expanded scattering matrices** from two-body elastic
  kinematics, discrete-level inelastic, continuum inelastic, and fission χ.
- Handles thermal scattering with a **three-regime model**: cold elastic
  above the free-gas cutoff, free-gas in the intermediate range, and
  $S(\alpha,\beta)$ below the thermal cutoff — with numerically stable
  $\ell \geq 1$ moments.
- Applies the **Ouisloumen–Sanchez Doppler-broadened resonance kernel** to
  heavy nuclides in the resolved resonance range.
- Applies **narrow-resonance p-table sampling** in the unresolved resonance
  region.
- Solves 1-D **fixed-source** slab transport with two deterministic solvers:
  **$S_N$** (diamond-difference) and **$P_N$** (cell-centred finite
  difference, $P_1$–$P_7$).
- Runs an end-to-end **multi-level condensation pipeline** — HFG (~24 000)
  → UFG (~2 000) → MG (70) — with flux-weighted group collapsing and
  cube-root-of-$T$ temperature interpolation.
- Carries a parallel research track, **Spectrum Expansion (SpecEx / RSE)**,
  that builds a POD basis over a 10 000-condition UFG spectrum sweep toward
  a reduced-order replacement for the UFG transport solve.

## Who this is for

The code is aimed at students and researchers in reactor physics who want to
experiment with ultra-fine-group treatments, Doppler-broadened resonance
kernels, and reduced-order spectrum methods without the black-box feel of
production lattice codes. The documentation assumes familiarity with the
linear Boltzmann equation and standard multi-group reactor-physics methods.

## Where to start

<div class="grid cards" markdown>

-   :material-book-open-variant: **Read the theory**

    ---

    Start with the [pipeline overview](theory/overview.md), then dive into the
    [slab $S_N$](theory/solvers_sn.md) / [slab $P_N$](theory/solvers_pn.md)
    derivations, the [resonance kernel](theory/resonance_kernel.md), the
    [URR p-table method](theory/urr_ptable.md), the
    [multi-level pipeline](theory/multilevel.md), or the
    The [SpecEx reduced-order method](specex/index.md) has its own section.

-   :material-console: **Build and run**

    ---

    The [install guide](guide/install.md) covers CMake + HDF5 on Linux and
    Windows. The [running guide](guide/running.md) walks through the PWR
    pin-cell multi-level demo.

-   :material-chart-line: **See the demonstrations**

    ---

    [Pin-cell verification](demos/pincell_verification.md) compares
    UFG-Transport against OpenMC; [multi-level flux](demos/multilevel_flux.md)
    visualises the HFG → UFG → MG collapse.

-   :material-github: **Source code**

    ---

    Hosted on GitHub. The main library is `ufg_lib`; driver executables are
    `ufg_transport_app`, `ufg_0d_app`, `ufg_1d_app`, and the SpecEx offline
    tools (`specex_sweep_app`, `specex_pod_app`, `specex_moments_app`).

</div>

## Status

The code is under active research development. The table below reflects the
ground truth — it is intentionally explicit about which pieces are verified,
which are partial, and which are disabled.

### Verified and in routine use

| Capability | Notes |
|---|---|
| ACE / OpenMC HDF5 reader | all distribution types, multi-temperature |
| UFG group-structure builder | uniform lethargy and coarse-partitioned |
| Homogenizer | multi-nuclide, multi-composition, volume fractions |
| Three-regime thermal kernel | cold elastic / free-gas / $S(\alpha,\beta)$, stable $\ell \geq 1$ moments |
| O&S resonance scattering kernel | heavy nuclides in resolved resonance range |
| NR p-table sampling for URR | total / elastic / fission cross sections |
| HFG → UFG → MG multi-level pipeline | flux-weighted collapse, auto-dispatch |
| Cube-root-of-$T$ temperature interpolation | all XS including scattering matrices |
| 0-D fixed-source solver | symmetric Gauss–Seidel with thermal inner loop |
| 1-D slab $S_N$ solver | diamond-difference, Gauss–Legendre, reflective/vacuum/white |
| 1-D slab $P_N$ solver | $P_1$–$P_7$, block Thomas, Marshak vacuum |
| CUDA XS processor (hybrid) | verified against CPU; ~2× end-to-end speedup |

### Partial — offline components only

| Capability | Notes |
|---|---|
| SpecEx Phase 1 (spectrum sweep) | parametric 0-D sweep over 10 000+ conditions |
| SpecEx Phase 2 (POD + WLRA) | block SVD with U-238 $\sigma_a$ weighting |
| SpecEx Phase 3 (XS moments) | $\sigma_{t,kl}$, $\sigma_{s,kl}$ per nuclide/T |

The on-line SpecEx transport solve (Phase 4) is not yet implemented. See the
[SpecEx section](specex/index.md) for the full plan.

### Disabled or not yet implemented

| Item | Status |
|---|---|
| 1-D cylindrical $S_N$ | **disabled** — throws in `main_1d.cpp`; documented in `coldcase/cylindrical_sn_transport_deficit.md` |
| 1-D CPM (cylindrical collision-probabilities) | present in source but **not yet verified** against a Monte Carlo reference to production standards; excluded from the user-facing solver set |
| $k$-effective eigenvalue solver | not implemented; fixed-source only |
| 2-D / 3-D geometry | out of scope for the current code; 2-D $P_N$ prototypes exist but XS integration is incomplete |
| Method of characteristics (MOC) | planned |
| SpecEx online solver (Phase 4) | planned |
| URR p-table temperature interpolation | planned |
| Anisotropic CPM, higher-order CPM | planned research |

### Known limits

- Source iteration has spectral radius ≈ 1 in the thermal up-scatter block.
  The solvers carry a dedicated thermal inner loop
  (`thermal_max_inners`, `thermal_tolerance`); a DSA acceleration scheme is a
  candidate for future work.
- OpenMP default is **4 threads**. Pass `--threads N` to override.
- The GPU path covers principal XS, cold elastic, and discrete-level
  inelastic; thermal, $S(\alpha,\beta)$, the O&S kernel, and continuum
  inelastic are CPU-only in the current hybrid.

## Current development focus

1. Wrap up the SpecEx offline library and start on the on-line solver
   (Phase 4).
2. Add a $k$-effective power iteration on top of the fixed-source driver.
3. Extend the CUDA hybrid to cover thermal / $S(\alpha,\beta)$ and the O&S
   kernel.
4. Re-investigate cylindrical geometry with either a corrected $S_N$ or a
   verified CPM path.

Older milestones are summarised in the [theory overview](theory/overview.md).
