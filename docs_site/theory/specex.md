# Spectrum Expansion (SpecEx / RSE)

**Spectrum Expansion (SpecEx)** — internally called RSE, for *Resonance
Spectrum Expansion* — is an orthogonal research track inside UFG-Transport
that aims to replace the entire resonance self-shielding stack (Bondarenko,
subgroup, HFG condensation, …) with a **reduced-order** model built from
proper orthogonal decomposition (POD) of UFG flux spectra.

The current status is a partial, offline library-generation implementation.
The final on-line transport solve (Phase 4) is not yet coded. This page
describes what SpecEx is, why it is worth exploring, what the code does
today, and what it is expected to do once complete.

The method follows Lee & Kim, *Direct Incorporation of Up-Scattering Sources
in the Resonance Spectrum Expansion (RSE) Method*, *Annals of Nuclear
Energy* (2026).

## Motivation — why a reduced-order method?

The UFG library is already a good trade: it sidesteps most of the
self-shielding machinery at the cost of carrying $G \sim 10^4$ groups through
the transport solve. Two pressure points remain:

1. **Up-scatter coupling at UFG resolution.** The thermal up-scatter block
   couples every pair of groups inside the Maxwellian region. At
   $G \sim 2 \times 10^4$, that block alone has $\sim 10^6$ non-zero entries
   per composition per temperature. Iterating on the full UFG scatter matrix
   under a geometry search (burnup, temperature feedback) is expensive.
2. **Self-shielding still needs flux-weighted collapse.** The HFG → UFG
   condensation step requires a 0-D solve per composition before transport,
   which becomes the dominant cost at $G = 24{,}000$.

The SpecEx observation is that UFG flux spectra, across wide sweeps of
composition, temperature, and background dilution, live on a **low-
dimensional manifold**. A typical LWR analysis needs ~95 POD basis functions
to span that manifold to the $10^{-3}$ collision-rate error level — two
orders of magnitude fewer degrees of freedom than the UFG grid.

If a transport solve can be expressed in the POD basis directly, then:

- the scattering transfer from $10^4 \times 10^4$ groups collapses to a
  $\sim 10^2 \times 10^2$ matrix of **cross-section moments**;
- up-scatter couplings come for free in the same small matrix;
- the library built once offline is reused across all on-line geometry /
  burnup / temperature queries.

## The end-to-end pipeline

```
                        OFFLINE (library generation)
  ┌──────────────────────────────────────────────────────────────────┐
  │   Phase 1            Phase 2            Phase 3                  │
  │   Spectrum           POD / Basis        XS Moment                │
  │   generation   ───►  extraction   ───►  generation               │
  │   (10k+ 0-D          (SVD per block,    (σ_{x,kl} per nuclide/T) │
  │    sweep)             WLRA rank sel.)                            │
  │                                                                  │
  │   spectra.h5         basis.h5            specex_library.h5       │
  └──────────────────────────────────────────────────────────────────┘

                        ONLINE (transport solve)
  ┌──────────────────────────────────────────────────────────────────┐
  │   Phase 4 (future)                                               │
  │   Load library → solve reduced SpecEx FSP for φ_k expansion     │
  │   → reconstruct effective MG / UFG flux → feed solver            │
  └──────────────────────────────────────────────────────────────────┘
```

## Phase 1 — Spectrum generation

A parametric sweep generates 0-D flux spectra across a Cartesian product of
conditions:

| Axis | Typical range | Purpose |
|---|---|---|
| Fuel composition | fresh UO₂, depleted UO₂, Gd-bearing UO₂, MOX | nuclide vector richness |
| Fuel temperature | 600 / 900 / 1200 / 1500 K | Doppler broadening |
| Moderator temperature | 500 / 600 K | thermal spectrum shift |
| Moderator density | 0.66 / 0.71 / 0.76 g/cm³ | void / density feedback |
| Background $\sigma_0$ | $10,\,50,\,100,\,500,\,10^4$ barn | self-shielding strength |

The background is injected as a **flat-in-energy absorber** with
$\sigma_t = \sigma_0$ and $N = 1\,\text{atom/barn-cm}$, the standard technique
for sweeping self-shielding strength in a homogeneous medium. Target library
size is 10 000+ spectra at UFG resolution.

Each sweep point reuses the existing 0-D pipeline: `XsProcessor::process()`
→ `mix()` → `compute_source_spectrum()` → `solve_fixed_source()`. Microscopic
cross sections are cached per `(nuclide, T)` pair; only `mix()` and the
solve run per condition.

The output `spectra.h5` contains the $N \times G$ snapshot matrix plus
metadata: energy boundaries, $\Delta E_g$, the reference U-238 $\sigma_a(g)$
vector used as the WLRA weight, and a condition-metadata table.

## Phase 2 — POD basis extraction

### Snapshot matrix and energy weighting

Assemble the snapshot matrix

$$
\mathbf{A} \;=\;
\begin{bmatrix} \boldsymbol{\varphi}_1 & \boldsymbol{\varphi}_2 & \cdots & \boldsymbol{\varphi}_N \end{bmatrix}^{\!\top},
\qquad \mathbf{A} \in \mathbb{R}^{N \times G}.
$$

Each column is normalised to unit energy-integral so that absolute flux
amplitude does not pollute the SVD. To force the POD basis to be
orthonormal under the continuous $L^2$ inner product, every column $g$ is
multiplied by $\sqrt{\Delta E_g}$ before the SVD, and the resulting basis
vectors are divided back by the same factor:

$$
\tilde{A}_{i,g} \;=\; \sqrt{\Delta E_g}\,\hat{\varphi}_{i}(g),
\qquad
\int b_k(E)\,b_l(E)\,dE \;=\; \delta_{kl}.
$$

### Block structure

The full UFG range is partitioned into ~13 coarse **energy blocks**. An
SVD is run independently inside each block. The reasoning: within a block
the flux shapes are more coherent, so a small number of basis functions is
enough; across blocks the shapes are decoupled by the resonance structure,
so mixing them in one SVD would waste singular values. Empirically a
p-refinement (more basis per block) is more effective than an h-refinement
(more blocks) for the same total rank.

### WLRA — weighted low-rank approximation

The plain rank-truncation SVD minimises $\|\mathbf{A} - \mathbf{A}_m\|_F$ —
the spectral approximation error. What the downstream solver actually cares
about is the **collision-rate error**,

$$
\varepsilon_R \;=\;
\frac{\bigl\| \mathbf{W}(\mathbf{A} - \mathbf{A}_m) \bigr\|_F}
     {\bigl\| \mathbf{W}\mathbf{A}\bigr\|_F},
\qquad
\mathbf{W} \;=\; \operatorname{diag}\bigl(\sigma_a^{\text{U-238}}(g,\,T_\text{ref})\bigr).
$$

The **weighted low-rank approximation (WLRA)** solves

$$
\min \bigl\|\mathbf{W}(\mathbf{A} - \mathbf{A}_m)\bigr\|_F^2
$$

via a projected gradient descent

$$
\mathbf{A}_m^{(i+1)}
\;=\; \mathcal{P}_m\!\left[
  \mathbf{A}_m^{(i)} - t\,\mathbf{W}^{\!\top}\mathbf{W}\bigl(\mathbf{A}_m^{(i)} - \mathbf{A}\bigr)
\right],
$$

where $\mathcal{P}_m$ is rank-$m$ truncation (a plain SVD), $t$ is a learning
rate, and the rank $m$ is incremented until $\varepsilon_R$ drops below a
user-specified threshold (typically 1 %). U-238 $\sigma_a$ is a physically
motivated weight: it emphasises the epithermal resonance region, where most
of the self-shielding error lives.

Typical result: **~95 basis functions total**, distributed across 13 blocks
with 5–12 basis per block. Compare with $G = 24{,}000$ UFG groups — a 250×
reduction in the number of spectral degrees of freedom.

Linear-algebra backend: Eigen 3 (`BDCSVD`) on the CPU by default;
cuSOLVER `gesvdj` on the GPU when CUDA is available.

## Phase 3 — Cross-section moments

For each nuclide and each temperature, the SpecEx library stores small
moment matrices obtained by projecting the point-wise cross sections onto
the basis:

### Total XS moments

$$
\sigma_{t,kl}(T)
\;=\;
\int \sigma_t(E,T)\,b_k(E)\,b_l(E)\,dE
\;\approx\;
\sum_g \sigma_t(g,T)\,b_k(g)\,b_l(g)\,\Delta E_g.
$$

### Scattering moments

A direct double-sum over the $G \times G$ scatter matrix is $O(G^2 K^2)$ and
too expensive at UFG resolution. The computation is factored into two
$O(G K)$ passes:

$$
S_k(E',T) \;=\; \int \sigma_s(E \to E', T)\,b_k(E)\,dE,
\qquad
\sigma_{s,kl}(T) \;=\; \int S_k(E',T)\,b_l(E')\,dE'.
$$

For ~95 basis functions, the final library has on the order of $10^3$ – $10^4$
non-zero moment-matrix entries per nuclide per temperature — compared with
$\sim 10^9$ entries in a full UFG scatter matrix.

The output `specex_library.h5` stores $\sigma_{t,kl},\,\sigma_{s,kl},\,
\sigma_{f,kl},\,\chi_k$ per nuclide per temperature, along with the basis
and energy-block metadata.

## Phase 4 — online transport solve (planned)

The on-line transport solver — the step that would turn SpecEx into a usable
replacement for the multi-group / UFG solve — is the next research item.
The sketch is:

1. Expand the unknown flux in the basis:
   $\varphi(E,\vec r) \approx \sum_k \phi_k(\vec r)\,b_k(E)$.
2. Project the transport equation onto each basis function $b_l$.
3. The streaming term retains its usual spatial discretisation; the
   collision term becomes $\sum_k \Sigma_{t,kl}(\vec r)\,\phi_k(\vec r)$; the
   scattering source becomes $\sum_k \Sigma_{s,kl}(\vec r)\,\phi_k(\vec r)$.
4. Solve a small $K \times K$ per-cell system coupled across space.

The effective multi-group flux for any downstream consumer (MG library,
depletion, OpenMC comparison) is recovered by evaluating
$\varphi(E_g) = \sum_k \phi_k\,b_k(E_g)$ — no separate MG library is needed.

## What SpecEx is expected to do

- Replace HFG → UFG condensation with a single offline basis build.
- Deliver UFG-level accuracy at roughly MG-level solve cost.
- Handle up-scatter coupling natively through the $\sigma_{s,kl}$ moments
  — no special thermal-inner loop needed.
- Serve as a research platform for reduced-order methods in reactor
  physics: basis adaptivity, cross-temperature extrapolation, and
  on-the-fly basis updates during depletion are all natural extensions.

## Current status

| Phase | Role | Status |
|---|---|---|
| 1 | Batch spectrum generation | implemented (`specex_sweep_app`) |
| 2 | POD / basis extraction (SVD + WLRA) | implemented (`specex_pod_app`) |
| 3 | XS moment generation | implemented (`specex_moments_app`) |
| 4 | On-line reduced-order transport solve | **not yet implemented** |

See the [status section on the home page](../index.md#status) for the full
picture of how SpecEx sits next to the rest of the code.

## Files

- `include/specex/specex_sweep.hpp`, `src/specex/…` — Phase 1.
- `include/specex/specex_pod.hpp`, `src/specex/specex_pod.cpp`,
  `src/specex/specex_pod_input.cpp`, `src/main_specex_pod.cpp` — Phase 2.
- `include/specex/specex_moments.hpp` — Phase 3 interface.
- `coldcase/rse_blueprint.md` — full design document, including the sweep
  axes, the WLRA derivation, and the Phase 4 skeleton.

## References

- H. G. Lee and K. S. Kim, *Direct Incorporation of Up-Scattering Sources in
  the Resonance Spectrum Expansion (RSE) Method*, *Annals of Nuclear Energy*,
  2026.
- K. S. Kim and M. L. Williams, *Resonance Spectrum Expansion Method in
  Deterministic Transport Calculations*, *Nuclear Science and Engineering*
  (historical foundation for the method).
