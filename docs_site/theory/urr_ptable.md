# URR probability tables (narrow-resonance sampling)

Between the **resolved resonance region** (RRR, typically up to a few keV for
actinides) and the smooth fast region, nuclear evaluations contain an
**unresolved resonance region (URR)** — an energy range in which individual
resonances overlap, so the published cross section is a statistical average
rather than a point-wise curve. Running a UFG library on the URR-averaged
point-wise data alone under-predicts resonance self-shielding: the dips of
individual (unresolved) resonances are averaged out before the flux-weighting
step can act on them.

ENDF MF 2 / LRU 2 provides **probability tables** (p-tables): at each energy
grid point $E_k$ in the URR, a histogram of $(\sigma_t, \sigma_a, \sigma_e,
\sigma_f)$ values together with a probability (CDF) for each band. UFG-Transport
uses the **narrow-resonance (NR)** approximation to fold these bands into the
UFG library.

## The narrow-resonance approximation

The NR approximation assumes that the neutron flux within one unresolved
resonance is proportional to $1/[\sigma_t(E) \cdot E]$ — i.e. that the
neutron's energy loss per scattering event is much larger than the resonance
width, so a neutron effectively sees a single resonance only once and then
escapes to lower energy. In the p-table language this translates to: the
expected cross section in a narrow URR band is the **CDF-weighted band
average**,

$$
\langle \sigma_x \rangle_g
\;=\;
\sum_{k=1}^{N_\text{band}} w_k\,\sigma_{x,k}(E_g),
\qquad
\sum_k w_k = 1,
$$

where $\sigma_{x,k}$ is the band-$k$ cross section at the energy $E_g$ (the
UFG group midpoint, located by bracketing in the URR energy grid), and $w_k$
is the probability of being in band $k$ (from the CDF column of the p-table).

The `multiply_smooth` option additionally scales each band value by the
smooth point-wise $\sigma_x(E_g)$ from the background file — useful when the
p-tables provide **ratios** rather than absolute cross sections (this is
evaluation-dependent).

### What NR does *not* do

- It does **not** solve a self-consistent slowing-down equation inside the URR
  at the p-table level — that would be the **intermediate-resonance (IR)** or
  **full flux-folding** approach. The NR result feeds the downstream HFG → UFG
  flux-weighted collapse, which does effectively see a self-shielded flux.
- It does **not** currently interpolate p-tables across temperature. For
  target temperatures that don't match an ACE temperature, the p-table is
  read from the nearest tabulated temperature and the downstream
  [cube-root-of-T interpolation](multilevel.md#temperature-interpolation)
  handles the smooth variation.

## What the URR block writes

For each UFG group $g$ whose midpoint falls inside the URR energy range of a
nuclide, the XS processor overrides the smooth-only values of:

- $\sigma_t^{(g)}$,
- $\sigma_e^{(g)}$,
- $\sigma_f^{(g)}$ (where applicable),

with the NR-expected value $\langle \sigma_x \rangle_g$. Capture is derived as
$\sigma_c = \sigma_t - \sigma_e - \sigma_f$ to preserve the balance.

The **scattering transfer matrix** $\Sigma_{s,\ell}^{\,g \to g'}$ is **not**
rebuilt from p-tables; it uses the smooth angular/energy distribution from
the background file. The reasoning is that the URR angular distributions are
statistical and near-isotropic; differentiating them between bands would
require sampling-level detail that ENDF does not publish.

## Pseudocode

```
for each UFG group g in URR [E_min, E_max]:
    (k_lo, k_hi, w_interp) = bracket_urr_energy_points(E_g)
    P = interpolate_ptable(P[k_lo], P[k_hi], w_interp)     # lin-lin or log-log
    for each reaction x in {sigma_t, sigma_e, sigma_f}:
        sigma_bar = sum_k P.probability[k] * P.sigma_x[k]
        if multiply_smooth:
            sigma_bar *= sigma_x_smooth(E_g)
        out.sigma_x[g] = sigma_bar
    out.sigma_c[g] = out.sigma_t[g] - out.sigma_e[g] - out.sigma_f[g]
```

All of this runs inside the normal XS processor loop, immediately after the
smooth group averaging and before the scattering-matrix assembly.

## Configuration

The URR override is on by default whenever p-tables are present in the ACE
library; it can be disabled via `skip_urr = true` in the XS processor config.
The reader handles the p-table dataset directly in `src/ace/ace_reader.cpp`
and exposes it through the `ufg::ace::ProbabilityTable` struct
(`include/ace/ace_table.hpp`).

## Verification and limits

The NR p-table override matches the URR-averaged OpenMC cross sections
point-for-point at each band-sampled energy — this is by construction,
because both are computing $\langle \sigma \rangle$ against the same CDF.

The limits of the NR approximation show up in two places:

1. **Very narrow URR levels at low energy**, where the narrow-resonance
   assumption ($\Delta u_\text{scatter} \gg \Gamma$) starts to break down.
   This is usually already covered by the RRR in modern evaluations.
2. **Heterogeneity self-shielding**, where the NR $1/\Sigma_t E$ flux shape
   inside the resonance does not apply to a pincell geometry. The multi-level
   pipeline's downstream 0-D HFG solve partially compensates, but a full
   IR / subgroup treatment would be more faithful. This is a candidate for a
   follow-up research item.

## Files

- `include/ace/ace_table.hpp` (`ProbabilityTable` struct).
- `src/ace/ace_reader.cpp` — p-table dataset loading.
- `src/xs/xs_processor.cpp` — NR expected-value sampling, called from the
  main XS processor loop when a group falls inside the URR.
- `coldcase/unresolved_resonance.md` — design notes and follow-up items
  (p-table temperature interpolation, IR extension).

## References

- L. B. Levitt, *The Probability Table Method for Treating Unresolved
  Neutron Resonances in Monte Carlo Calculations*, Nuclear Science and
  Engineering **49**, 450 (1972).
- D. E. Cullen, *Application of the Probability Table Method to Multigroup
  Calculations*, UCRL-79761 (1977).
- R. E. MacFarlane, *NJOY: The Nuclear Data Processing System*
  (LA-UR-17-20093), UNRESR module.
- ENDF-6 Manual, Chapter 2 (File 2, LRU = 2).
