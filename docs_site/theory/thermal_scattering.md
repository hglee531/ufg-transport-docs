# Thermal scattering

Below the slowing-down region, the classical two-body scattering kernel breaks down:
the target nucleus is no longer stationary, and in moderators it is no longer even
free — it is bound into a molecule or crystal lattice. UFG-Transport treats this
regime with a **three-regime model** that cleanly covers the full energy range.

## The three regimes

```
    E > 1 keV              ←  cold elastic, stationary target (classical)
    E_sab < E ≤ 1 keV      ←  free-gas kernel (thermal agitation of free nucleus)
    E ≤ E_sab              ←  S(α,β) bound kernel + residual free-gas for unbound nuclides
```

- **Cold elastic** uses the two-body lab kinematics described in
  [UFG grid & scattering](ufg_grid.md).
- **Free-gas** accounts for the Maxwellian velocity distribution of the target at
  temperature $T$; mandatory in the range where targets move but are not yet bound.
- **S(α,β)** is the ENDF File 7 thermal-scattering law for nuclides chemically bound
  in a specific matrix (H in H₂O, C in graphite, Zr in ZrH, …).

The transition energy $E_\text{sab}$ is set per thermal table — typically 4–5 eV
for common moderators.

## The bound-atom cross section $\sigma_b$

The S(α,β) and free-gas formalisms normalise by the *bound-atom* cross section

$$
\sigma_b \;=\; \sigma_\text{free}\,\left(\frac{A+1}{A}\right)^2,
$$

where $A$ is the atomic mass ratio. This reduced-mass factor accounts for the fact
that the fundamental scattering amplitude is defined in the centre-of-mass frame,
whereas the experimentally reported $\sigma_\text{free}$ is a lab-frame quantity.
For hydrogen ($A \approx 1$) the enhancement is a factor of 4; for heavy
nuclei $\sigma_b \to \sigma_\text{free}$.

## The free-gas kernel

The double-differential free-gas kernel, angle-integrated, is

$$
\sigma_{s}^\text{fg}(E \to E') \;=\;
\frac{\sigma_b\,A}{4E}\,
F\!\bigl(E, E', T\bigr),
$$

with $F$ built from error functions of the standard thermal variables. The
implementation in `src/xs/xs_processor.cpp` uses:

- a **constant $\sigma_b$** sampled at 1 MeV rather than a re-evaluated
  energy-dependent bound cross section — avoiding artefacts near resonances;
- **8-point source sub-sampling × 4-point Gauss–Legendre dest quadrature** over
  each source/destination group pair;
- a **window floor of 15 $kT$** below the source energy to capture significant
  up-scatter transfers while avoiding ∞-tail cost;
- **mass-ratio gating** ($A \leq 10$): above this, a cold elastic treatment is used
  even inside the nominally-thermal range. Heavy nuclei barely move at reactor
  temperatures, so free-gas would only add numerical noise.

A $k_BT$ fallback

$$
kT \;\approx\; 8.617 \times 10^{-5}\;\text{eV/K}\,\times\,T[\text{K}]
$$

is used for nuclides whose HDF5 file is missing the `kTs/<temp>K` dataset (this
happened historically for H-1 and O-16 in several libraries).

## S(α,β) for bound scatterers

For bound nuclides the incoherent inelastic contribution is

$$
\sigma_{s}^\text{sab}(E\to E',\mu) \;=\;
\frac{\sigma_b}{2kT}\,\sqrt{\frac{E'}{E}}\,
e^{-\beta/2}\,S(\alpha,\beta),
$$

with the usual dimensionless momentum- and energy-transfer variables

$$
\alpha \;=\; \frac{E + E' - 2\mu\sqrt{EE'}}{A\,kT},
\qquad
\beta \;=\; \frac{E' - E}{kT}.
$$

The S(α,β) table is read from the OpenMC thermal HDF5 file; angular integration
uses the tabulated $(\mu,\text{PDF})$ discrete sets. Coherent elastic is additively
included for crystalline binders (graphite, beryllium).

## Why three regimes — and not two?

Merging the free-gas and cold-elastic regimes (i.e. using free-gas everywhere below
some cutoff) causes spurious up-scatter noise for heavy nuclei that barely
thermalise. Conversely, using cold elastic too low in energy misses the thermal
up-scatter that drives the Maxwellian spectrum. The three-regime split lets light
moderators see free-gas down to $E_\text{sab}$ and then S(α,β) below it, while heavy
fuel nuclides use cold elastic throughout — at the cost of one extra threshold.

The unit test results in the project history showed that the three-regime model
(with the 8×4 sub-sampling and mass gating described above) reproduces a smooth
Maxwellian flux at 600 K in a UO₂+H₂O mixture, converging in ~10 symmetric
Gauss–Seidel iterations with residual $\max_g |r_g| < 7\times10^{-11}$.

## Further reading (inside the repo)

- `docs/theory_thermal_kernel_prefactor.md` — full derivation of the free-gas
  prefactor.
- `docs/report_free_gas_L1_stable_moments.md` — higher-order moments and stability.
- `docs/report_analytic_free_gas_kernel.md` — analytic cross-checks.
- `src/xs/thermal_kernel.cpp`, `src/xs/xs_processor.cpp` — implementation.

## External references

- R. E. MacFarlane, *NJOY: The Nuclear Data Processing System* (LA-UR-17-20093) —
  Chapters on THERMR and GROUPR.
- D. E. Parks, *Thermal Neutron Scattering*, ORNL-TM-1234.
- ENDF-6 Manual, Appendix D: Thermal Scattering Law.
