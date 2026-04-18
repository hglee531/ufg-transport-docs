# Thermal scattering

Below the slowing-down region, the stationary-target kinematic kernel breaks
down: the target nucleus is no longer at rest, and in moderators it is not even
free — it is bound into a molecule or a crystal lattice. UFG-Transport treats
this regime with a **three-regime model** that covers the full energy range.

## The three regimes

```
    E > E_fg                      ←  cold elastic, stationary target
    E_sab < E ≤ E_fg              ←  free-gas kernel (thermally moving free nucleus)
    E ≤ E_sab                     ←  S(α,β) bound kernel + optional free-gas mix
```

- **Cold elastic** uses the two-body lab kinematics in
  [UFG grid & scattering](ufg_grid.md#elastic-two-body-kinematics-stationary-target).
- **Free-gas** accounts for the Maxwell–Boltzmann velocity distribution of the
  target nucleus at temperature $T$.
- **$S(\alpha,\beta)$** is the ENDF File 7 thermal-scattering law for nuclides
  chemically bound in a specific matrix (H in H₂O, C in graphite, Zr in ZrH, …).

The free-gas upper cutoff $E_\text{fg}$ is controlled by
`free_gas_cutoff_ev` (default **400 eV**, matching OpenMC); the free-gas /
$S(\alpha,\beta)$ transition $E_\text{sab}$ is set per thermal table — typically
4 – 5 eV for common moderators. For heavy nuclides $A > A_\text{max}$ (default
**300**), the cold-elastic kernel is used throughout: heavy nuclei barely move
at reactor temperatures, so the free-gas kernel would only add numerical noise.

The Doppler-broadened resonance kernel ([O&S](resonance_kernel.md)) operates
in the **resonance** energy range ($\lesssim 1$ keV, $A \gtrsim 10$) and can
also substitute for the cold-elastic treatment there when enabled; it is
independent of the free-gas / $S(\alpha,\beta)$ machinery.

## The bound-atom cross section (σ_b)

The free-gas and $S(\alpha,\beta)$ formalisms normalise by the **bound-atom**
cross section

$$
\sigma_b \;=\; \sigma_\text{free}\,\left(\frac{A+1}{A}\right)^2,
$$

where $A$ is the atomic mass ratio. The reduced-mass factor accounts for the
fact that the fundamental scattering amplitude is defined in the
centre-of-mass frame, whereas the tabulated $\sigma_\text{free}$ is a
lab-frame quantity. For hydrogen ($A \approx 1$) the enhancement is ~4×; for
heavy nuclei $\sigma_b \to \sigma_\text{free}$.

## The free-gas kernel

### Differential form

The double-differential free-gas kernel, integrated over scattering angle,
reads

$$
\sigma_s^{\text{fg}}(E \to E')
\;=\; \frac{\sigma_b\,A}{4\,E}\,
  \int_{\alpha_\min}^{\alpha_\max}\!S(\alpha,\beta)\,d\alpha,
$$

with the dimensionless momentum- and energy-transfer variables

$$
\alpha \;=\; \frac{E + E' - 2\mu_0\sqrt{E\,E'}}{A\,k_B T},
\qquad
\beta \;=\; \frac{E' - E}{k_B T},
$$

and the short-collision-time Gaussian scattering law

$$
S(\alpha,\beta) \;=\; \frac{1}{\sqrt{4\pi\alpha}}\,
  \exp\!\left[-\frac{(\alpha + \beta)^2}{4\alpha}\right].
$$

The $\alpha$-integration limits correspond to the $\mu_0 = \pm 1$ kinematic
extremes,

$$
\alpha_\min \;=\; \frac{(\sqrt{E'} - \sqrt{E})^2}{A\,k_B T},
\qquad
\alpha_\max \;=\; \frac{(\sqrt{E'} + \sqrt{E})^2}{A\,k_B T}.
$$

### Legendre moments with stable evaluation

For Legendre order $\ell$, a change of variables $\mu_0 \to t$ reduces the
kernel to an integral of the form $e^{-\beta/2} J_k(\beta)$, where the factor
$e^{-\beta/2}$ can cause catastrophic cancellation at large $|\beta|$ (for
example $\beta \approx -50$ at $E = 2.6$ eV, $T = 600$ K).

UFG-Transport uses the rescaled form

$$
e^{-\beta/2}\,J_k \;=\; \int e^{-u^2}\,dt,
\qquad u \;=\; \tfrac{t}{2} \pm \tfrac{|\beta|}{2t},
$$

with $J_0$ via a stable `erfcx` formula, $J_1$ via differentiation w.r.t.
$a = 1/4$, and higher moments through an erf-based antiderivative recurrence.
This produces numerically stable $\ell \geq 1$ moments throughout the thermal
range.

### Numerical parameters

| Parameter | Default | Role |
|---|---|---|
| `num_quad_points` | 64 | Gauss–Legendre points for the $\alpha$-integral |
| source sub-samples | 16 (≥ `num_quad_points / 4`) | equi-lethargy sub-samples per source group |
| `free_gas_cutoff_ev` | 400 eV | upper energy for the free-gas regime |
| `free_gas_A_max` | 300 | maximum AWR; heavier nuclides use cold elastic |
| `sab_source_subsamples` | 8 | equi-lethargy sub-samples inside $S(\alpha,\beta)$ |

A $k_BT$ fallback

$$
k_B T \;\approx\; 8.617 \times 10^{-5}\,\text{eV/K}\;\times\;T[\text{K}]
$$

is used for nuclides whose HDF5 file is missing the `kTs/<temp>K` dataset
(which historically affected some H-1 and O-16 evaluations).

## S(α,β) for bound scatterers

For bound nuclides the incoherent inelastic contribution is

$$
\frac{d^2\sigma^{\text{sab}}_s}{d\Omega\,dE'}
\;=\;
\frac{\sigma_b}{4\pi\,k_B T}\,\sqrt{\frac{E'}{E}}\,
  e^{-\beta/2}\,S(\alpha,\beta),
$$

with the tabulated $S(\alpha,\beta)$ read from the OpenMC thermal HDF5 file
and the angular integration performed on the stored $(\mu,\mathrm{PDF})$
discrete sets. Coherent elastic is additively included for crystalline binders
(graphite, beryllium).

For $E \leq E_\text{sab}$, the processor can optionally **mix** the
$S(\alpha,\beta)$ contribution with a residual free-gas contribution using the
`sab_fraction` parameter (1.0 = pure $S(\alpha,\beta)$, which is the default
for bound nuclides).

## Why three regimes, and not two?

Merging the free-gas and cold-elastic regimes (i.e. using free-gas everywhere
below some cutoff) produces spurious up-scatter noise for heavy nuclei that
barely thermalise. Using cold elastic too low in energy misses the thermal
up-scatter that drives the Maxwellian spectrum. The three-regime split lets
light moderators see free-gas down to $E_\text{sab}$ and $S(\alpha,\beta)$
below it, while heavy fuel nuclides use cold elastic throughout, at the cost
of one extra threshold.

Unit tests in the project history confirm that the three-regime model
reproduces a smooth Maxwellian flux at 600 K in a UO₂ + H₂O mixture, with the
0-D solver converging to $\max_g |r_g| < 10^{-10}$ in a handful of
symmetric-Gauss–Seidel iterations at UFG resolution.

## Further reading

### Inside the repo

- `include/xs/thermal_kernel.hpp`, `src/xs/thermal_kernel.cpp` — free-gas
  kernel and stable $\ell \geq 1$ moments.
- `src/xs/xs_processor.cpp` — three-regime routing, $S(\alpha,\beta)$ driver.
- `docs/theory_thermal_kernel_prefactor.md` — full derivation of the free-gas
  prefactor.
- `docs/report_free_gas_L1_stable_moments.md` — stable higher-moment formulas.

### External

- R. E. MacFarlane, *NJOY: The Nuclear Data Processing System*
  (LA-UR-17-20093) — chapters on THERMR and GROUPR.
- D. E. Parks, *Thermal Neutron Scattering*, ORNL-TM-1234.
- ENDF-6 Manual, Appendix D: Thermal Scattering Law.
