# UFG-Transport

**Ultra-fine-group deterministic neutron transport with an ACE → UFG → solver pipeline.**

UFG-Transport is a research-oriented neutron transport code that processes point-wise ACE
cross-section data directly into ultra-fine multi-group libraries (up to tens of thousands
of groups), then solves the transport equation deterministically. It is written in
C++17 with a staged CUDA backend.

![Pin-cell CPM vs OpenMC](assets/images/pincell_cpm_vs_openmc.png)
/// caption
CPM flux in a PWR pin cell computed by UFG-Transport (2,000-group UFG) compared
with an OpenMC Monte Carlo reference.
///

## What it does

- Reads point-wise cross sections from **OpenMC-style HDF5 ACE libraries** (ENDF/B-VIII).
- Builds **ultra-fine group (UFG)** libraries with either uniform-lethargy or
  coarse-partitioned group structures — thousands of groups routinely, tens of thousands
  if desired.
- Constructs **Legendre-expanded scattering matrices** (elastic two-body kinematics,
  discrete-level inelastic, fission χ).
- Handles thermal scattering with a **three-regime model**: cold elastic above 1 keV,
  free-gas in the intermediate range, and S(α,β) kernels in the thermal region.
- Solves **1-D fixed-source transport** with three deterministic solvers:
  slab **S<sub>N</sub>**, slab **P<sub>N</sub>**, and cylindrical **CPM**.
- Runs a full **multi-level condensation pipeline** — HFG (≈24 000 groups) → UFG
  (≈2 000 groups) → MG (70 groups) — with flux-weighted group collapsing.

## Who this is for

The code is aimed at students and researchers in reactor physics who want to
experiment with ultra-fine-group treatments without the black-box feel of production
lattice codes. The documentation assumes familiarity with the linear Boltzmann equation
and standard multi-group methods.

## Where to start

<div class="grid cards" markdown>

-   :material-book-open-variant: **Read the theory**

    ---

    Start with the [pipeline overview](theory/overview.md), then dive into the
    [solver derivations](theory/solvers_sn.md) or the
    [multi-level condensation scheme](theory/multilevel.md).

-   :material-console: **Build and run**

    ---

    The [install guide](guide/install.md) covers CMake + HDF5 on Linux and Windows.
    The [running guide](guide/running.md) walks through the PWR pin-cell demo.

-   :material-chart-line: **See the demonstrations**

    ---

    [Pin-cell verification](demos/pincell_verification.md) compares UFG-Transport
    against OpenMC; [multi-level flux](demos/multilevel_flux.md) visualises the
    HFG → UFG → MG collapse.

-   :material-github: **Source code**

    ---

    Hosted on GitHub. The main library is `ufg_lib`; three drivers build on top:
    `ufg_transport_app`, `ufg_0d_app`, `ufg_1d_app`.

</div>

## Status

The CPU reference implementation is complete end-to-end: ACE parsing, UFG group
averaging, homogenization, and all three 1-D transport solvers. A k-effective
eigenvalue solver and a CUDA port are the main items still in flight.
