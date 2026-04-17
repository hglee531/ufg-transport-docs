# Multi-level HFG → UFG → MG condensation

For problems where even the UFG solve is too expensive to run *everywhere* in
space, UFG-Transport exposes a **multi-level** workflow that uses a
hyper-fine-group (HFG) library to build a problem-specific weighting spectrum, then
collapses to a UFG library for the transport solve, and optionally collapses
further to a standard multi-group (MG) library for downstream codes.

The canonical resolution hierarchy is

```
HFG   ~ 24,000 groups   (0-D per composition, Σ_t - Σ_s)φ = χ)
UFG   ~  2,000 groups   (1-D transport, Sn / Pn / CPM)
MG    ~     70 groups   (CASMO-70, UO2, etc. — library output only)
```

Full execution flow is documented in `docs/stage1_hfg_ufg_pipeline.md`. This page
gives the rationale, the equations, and the temperature-interpolation scheme.

## Why a multi-level pipeline?

Direct UFG transport on a 1-D pin-cell with $G = 24{,}000$ groups and $S_8$ angular
order produces matrices several GB in size and iteration counts measured in hours.
However, *within a single composition*, the 0-D slowing-down spectrum is already a
very good estimator of the true spatial flux shape for weighting.

The multi-level scheme therefore:

1. Builds an HFG library per **composition** (not per cell).
2. Solves a 0-D fixed-source problem per composition, with a shared fission source.
3. Uses the resulting 0-D flux as the weighting spectrum for HFG → UFG collapse.
4. Runs the 1-D transport problem at UFG resolution, with the composition-specific
   UFG library.

This replaces the naive UFG solve with a cheap 0-D HFG solve plus a much cheaper
UFG transport solve.

## Input format

The multi-level input uses nested sub-blocks inside `[GroupStructure]`:

```text
[GroupStructure]
  [hyperfine]
    type          = hfg
    energy_bin    = '20.0 1.0 1.0e-3 1.0e-5 0.5e-6 1.0e-8 1.0e-11'
    energy_groups = '2250 7000 7000 7000 500 250'
  []
  [ultrafine]
    type          = ufg
    energy_bin    = '20.0 1.0e-3 0.5e-6 1.0e-11'
    energy_groups = '200 1500 200 100'
  []
  [multigroup]
    type    = mg
    library = CASMO-70
  []
[]
```

When `main_1d.cpp` detects `group_structures.is_multilevel()`, it runs the
HFG → UFG pipeline automatically. Single-level inputs use the legacy direct-load
path.

## Stage 1: HFG → UFG

### Step 1. HFG library per composition

For each composition $c$ (fuel, moderator, …), the XS processor builds

$$
\Sigma_{t,c}^{(g)},\quad \Sigma_{a,c}^{(g)},\quad \nu\Sigma_{f,c}^{(g)},\quad
\Sigma_{s,c}^{(g\to g')},\quad \chi_c^{(g)}
$$

at HFG resolution.

### Step 2. Shared χ

A single $\chi_g$ is chosen from the composition with the largest integrated
$\nu\Sigma_f$ (normally the fuel) and broadcast to all compositions. This avoids
mismatched fission spectra between zones and makes the 0-D problems comparable.

### Step 3. 0-D fixed-source per composition

For each composition, solve

$$
\bigl(\Sigma_{t,c}^{(g)} - \Sigma_{s,0,c}^{(g\to g)}\bigr)\phi_c^{(g)}
\;-\; \sum_{g'\ne g} \Sigma_{s,0,c}^{(g'\to g)}\phi_c^{(g')}
\;=\; \chi_\text{shared}^{(g)}.
$$

UFG-Transport uses a **symmetric Gauss–Seidel** sweep (forward + backward per
iteration); at 2000+ groups this converges in ~10 iterations with residual
$\max_g|r_g| < 10^{-10}$.

### Step 4. Flux-weighted collapse

With the HFG flux $\phi_c^{(g)}$ in hand, each UFG group $G$ gets

$$
\Sigma_{x,c}^{(G)} \;=\;
\frac{\displaystyle\sum_{g \in G} \Sigma_{x,c}^{(g)}\,\phi_c^{(g)}\,\Delta u_g}
     {\displaystyle\sum_{g \in G} \phi_c^{(g)}\,\Delta u_g}.
$$

Scattering matrices are collapsed with the **fine-to-coarse map** so that
$\Sigma_{s,c}^{(G\to G')}$ is preserved under the weighted sum. The collapsed
library is ready for a 1-D transport solve.

## Stage 2 (optional): UFG → MG

The same flux-weighted collapse, applied with a broader fine-to-coarse map, takes
the UFG library to a user-specified MG structure. Predefined libraries are
exposed via `include/xs/mg_libraries.hpp`; **CASMO-70** is included.

The MG library is an *output* of the code — useful for feeding into coarse-group
lattice codes or for comparison against published benchmarks.

## Temperature interpolation: cube-root of T

Each composition may be assigned a non-tabulated temperature $T_\text{target}$. The
code finds the bracketing ACE temperatures $T_\text{lo} \leq T_\text{target} \leq T_\text{hi}$
and interpolates the *already-condensed* UFG cross sections in
$\tau = T^{1/3}$:

$$
w \;=\; \frac{T_\text{target}^{1/3} - T_\text{lo}^{1/3}}{T_\text{hi}^{1/3} - T_\text{lo}^{1/3}},
\qquad
\Sigma_x(T_\text{target}) \;=\; (1 - w)\,\Sigma_x(T_\text{lo}) \;+\; w\,\Sigma_x(T_\text{hi}).
$$

The cube-root-of-T variable is a classical NJOY interpolation choice that preserves
resonance-integral behaviour with limited tabulated temperatures. Interpolation is
applied to every XS type including scattering matrices.

## Cross-references

- `docs/stage1_hfg_ufg_pipeline.md` — execution flow, call graph, file map.
- `docs/report_20260305_condenser_fractional_weights.md` — details on partial
  fine-group membership at coarse boundaries.
- `include/xs/condenser.hpp`, `src/xs/condenser.cpp` — implementation.
- `include/xs/temperature_interpolator.hpp` — $\tau = T^{1/3}$ interpolator.
- `include/xs/mg_libraries.hpp` — predefined MG boundaries.
