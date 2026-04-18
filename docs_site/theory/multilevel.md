# Multi-level HFG → UFG → MG condensation

For problems where even the UFG solve is too expensive to run *everywhere* in
space, UFG-Transport exposes a **multi-level** workflow that uses a
hyper-fine-group (HFG) library to build a problem-specific weighting spectrum
per composition, collapses to a UFG library for the transport solve, and
optionally collapses further to a standard multi-group (MG) library for
downstream codes.

The canonical resolution hierarchy is

```text
HFG   ~ 24,000 groups   (0-D per composition, (Σ_t - Σ_s)φ = χ)
UFG   ~  2,000 groups   (1-D transport, Sn / Pn)
MG    ~     70 groups   (CASMO-70, library output only)
```

## Why a multi-level pipeline?

Direct UFG transport on a 1-D pin-cell with $G = 24{,}000$ groups and $S_8$
angular order produces scatter matrices of several GB and iteration counts
measured in hours. However, within a single composition the 0-D
slowing-down spectrum is already a very good estimator of the true spatial
flux shape used for weighting.

The multi-level scheme therefore:

1. Builds an HFG library per **composition** (not per cell).
2. Solves a 0-D fixed-source problem per composition with a shared $\chi$.
3. Uses the resulting 0-D flux as the weighting spectrum for HFG → UFG
   collapse.
4. Runs the 1-D transport problem at UFG resolution with the
   composition-specific UFG library.

This replaces the naive full-UFG transport solve with a cheap 0-D HFG solve
plus a much cheaper UFG transport solve.

An orthogonal track — [Spectrum Expansion](specex.md) — aims to replace the
HFG → UFG step entirely with a reduced-order basis. The two approaches are
complementary: the multi-level pipeline is the practical "works today"
workflow, while SpecEx is the research direction toward an order-of-magnitude
reduction in degrees of freedom.

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
HFG → UFG pipeline automatically. Single-level inputs use the legacy
direct-load path.

## Stage 1: HFG → UFG

### Step 1. HFG library per composition

For each composition $c$ (fuel, moderator, …), the XS processor builds

$$
\Sigma_{t,c}^{(g)},\quad
\Sigma_{a,c}^{(g)},\quad
\nu\Sigma_{f,c}^{(g)},\quad
\Sigma_{s,\ell,c}^{\,g\to g'},\quad
\chi_c^{(g)}
$$

at HFG resolution. The [URR p-table override](urr_ptable.md) and — when
enabled — the [O&S resonance kernel](resonance_kernel.md) both write into
this HFG library.

### Step 2. Shared χ

A single $\chi_g$ is chosen from the composition with the largest integrated
$\nu\Sigma_f$ (normally the fuel) and broadcast to all compositions. This
avoids mismatched fission spectra between zones and makes the 0-D problems
directly comparable as weighting calculations.

### Step 3. 0-D fixed-source solve per composition

For each composition, solve the 0-D slowing-down balance

$$
\bigl(\Sigma_{t,c}^{(g)} - \Sigma_{s,0,c}^{\,g\to g}\bigr)\,\phi_c^{(g)}
\;-\;
\sum_{g' \ne g} \Sigma_{s,0,c}^{\,g'\to g}\,\phi_c^{(g')}
\;=\;
\chi_\text{shared}^{(g)}.
$$

UFG-Transport uses a **symmetric Gauss–Seidel** sweep (forward + backward
per iteration); at tens of thousands of groups this converges in ~10–100
iterations with residual $\max_g |r_g| < 10^{-10}$. Up-scatter is resolved
by a dedicated thermal inner loop that re-sweeps the up-scatter block until
it stabilises, same pattern as the 1-D solvers.

### Step 4. Flux-weighted condensation

With the HFG flux $\phi_c^{(g)}$ in hand, each UFG group $G$ receives

$$
\Sigma_{x,c}^{(G)}
\;=\;
\frac{\displaystyle \sum_{g \in G} \Sigma_{x,c}^{(g)}\,\phi_c^{(g)}\,\Delta u_g}
     {\displaystyle \sum_{g \in G} \phi_c^{(g)}\,\Delta u_g},
$$

with the lethargy widths $\Delta u_g$ ensuring that the weighting is the
natural lethargy-integrated flux. Scattering matrices are collapsed with the
**fine-to-coarse map** so that $\Sigma_{s,\ell,c}^{\,G\to G'}$ is preserved
under the weighted sum. Partial fine-group membership at coarse boundaries
is handled with fractional weights (`docs/report_20260305_condenser_fractional_weights.md`).

The collapsed UFG library is then used directly in the 1-D transport solve.

## Stage 2 (optional): UFG → MG

The same flux-weighted collapse, applied with a broader fine-to-coarse map,
takes the UFG library to a user-specified MG structure. Predefined libraries
are exposed via `include/xs/mg_libraries.hpp`; **CASMO-70** is included out
of the box, and additional libraries can be registered via
`get_library_boundaries(name)`.

The MG library is an **output** of the code — useful for feeding into
coarse-group lattice codes or for comparison against published benchmarks.

## Temperature interpolation: cube root of $T$

Each composition may be assigned a target temperature $T_\text{target}$ that
is not among the ACE-tabulated temperatures. The code then finds the
bracketing ACE temperatures $T_\text{lo} \leq T_\text{target} \leq T_\text{hi}$
and linearly interpolates the **already-condensed** UFG cross sections in
$\tau = T^{1/3}$:

$$
w \;=\;
\frac{T_\text{target}^{1/3} - T_\text{lo}^{1/3}}
     {T_\text{hi}^{1/3} - T_\text{lo}^{1/3}},
\qquad
\Sigma_x(T_\text{target}) \;=\; (1-w)\,\Sigma_x(T_\text{lo}) \;+\; w\,\Sigma_x(T_\text{hi}).
$$

The cube-root-of-$T$ variable is the classical NJOY interpolation choice:
Doppler broadening scales with $\sqrt{T}$, and the empirical fit
$\tau = T^{1/3}$ preserves resonance-integral behaviour over the
reactor-relevant range ($\sim 300$ – $1500$ K) with only two or three
tabulated temperatures.

Interpolation is applied to every XS type, **including all Legendre orders
of the scattering matrix** (sparse and dense formats are handled
transparently). Out-of-range temperatures throw rather than extrapolate —
the user must add a bracketing ACE temperature to the library.

For the URR p-tables, interpolation across temperature is not yet
implemented (see [URR p-tables](urr_ptable.md)); the p-table at the nearest
ACE temperature is used, and the downstream cube-root-of-$T$ interpolation
of the resulting UFG cross sections picks up most of the smooth variation.

## Cross-references

- `docs/stage1_hfg_ufg_pipeline.md` — execution flow, call graph, file map.
- `docs/report_20260305_condenser_fractional_weights.md` — partial fine-group
  membership at coarse boundaries.
- `include/xs/condenser.hpp`, `src/xs/condenser.cpp` — HFG → UFG driver.
- `include/xs/temperature_interpolator.hpp`,
  `src/xs/temperature_interpolator.cpp` — $\tau = T^{1/3}$ interpolator.
- `include/xs/mg_libraries.hpp` — predefined MG boundaries.
