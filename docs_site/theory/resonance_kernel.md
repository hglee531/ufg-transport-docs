# Resonance scattering kernel (Ouisloumen–Sanchez)

For heavy resonant nuclides, the stationary-target elastic kinematic kernel
used in the [UFG scatter matrix](ufg_grid.md#elastic-two-body-kinematics-stationary-target)
is not accurate inside resonances: the target nucleus has a thermal
Maxwell–Boltzmann velocity distribution, and the cross section varies rapidly
across the thermal energy width of the target. A simple cold kinematic kernel
can over- or under-predict up-scatter by a factor of several in the first few
eV above a large resonance, distorting the resonance integral and therefore
$k_\text{eff}$.

UFG-Transport implements the **Ouisloumen–Sanchez (O&S) Doppler-broadened
elastic scattering kernel** (Ouisloumen & Sanchez, *Nucl. Sci. Eng.*, 1991) for
heavy nuclides in the resolved resonance range. It produces the same
Legendre-moment representation as the cold elastic kernel, so the transport
solver is agnostic to whether a row of the scatter matrix came from cold
elastic or from the resonance kernel.

## Why a Doppler-broadened kernel?

The cold kernel assumes that every scattering event happens off a target at
rest. In reality, the target speed $v_T$ follows a Maxwell–Boltzmann
distribution

$$
P(\mathbf v_T) \;=\;
\left(\frac{M}{2\pi k_B T}\right)^{3/2}
\exp\!\left(-\frac{M v_T^2}{2 k_B T}\right),
$$

so the **effective** elastic cross section in the lab frame is the incident
cross section integrated against the relative-speed distribution between
neutron and target. Inside a narrow resonance, where $\sigma(E)$ changes
appreciably over an energy scale comparable to $k_B T \cdot (m / M)$, this
integral cannot be replaced by its value at rest.

For light nuclides ($A \lesssim 10$) the free-gas / $S(\alpha,\beta)$
thermal kernels already include this effect (they are both derived from the
same Maxwell–Boltzmann target-motion picture). For heavy nuclides with strong
resonances — typically the actinides and structural nuclides above
$A \approx 10$ — the dedicated O&S kernel is used instead, over a broader
energy range than the thermal kernels.

## The O&S formulation

Let $E$ and $E'$ be the incident and outgoing neutron energies in the lab
frame, $\mu_L$ the cosine of the lab scattering angle, and $T$ the target
temperature. Introduce the dimensionless variables

$$
x = \sqrt{\frac{A E}{k_B T}},
\qquad
x' = \sqrt{\frac{A E'}{k_B T}},
$$

and the target-speed parameter $t$ (the integration variable). The
Doppler-broadened double-differential elastic cross section can be written
(O&S 1991, Eqs. 1–4) as

$$
\sigma_s(E \to E', \mu_L) \;=\;
\frac{1}{2\pi}\,
\sum_{n=0}^{\infty}\,(2n+1)\,\sigma_{s,n}(E \to E')\,P_n(\mu_L),
$$

with Legendre moments

$$
\sigma_{s,n}(E \to E')
\;=\;
\frac{\sigma_b\,A}{4\,E}\;
\int_{t_\min}^{t_\max} \sigma_e\!\bigl(E_\text{rel}(t)\bigr)\;
    \psi_n(t)\,dt.
$$

Here $\sigma_e$ is the point-wise elastic cross section of the nuclide, the
function $E_\text{rel}(t)$ is the relative-velocity kinematic map for a given
target speed parameter $t$, and $\psi_n(t)$ are the **Legendre transfer
functions**, evaluated by (O&S Eqs. A.1–A.2 for $n = 0, 1$; Eqs. A.3–A.5 for
$n \geq 2$):

- **$\psi_0(t)$** — available analytically, reduces to a difference of
  complementary error functions (`stable_erf_diff` in the code).
- **$\psi_1(t)$** — also analytic (O&S Eq. A.2).
- **$\psi_n(t)$ for $n \geq 2$** — an inner azimuthal integral $Q_n(x,t)$
  followed by an outer integral over a normalised speed variable. Both are
  done numerically by Gauss–Legendre quadrature.

The integration limits $t_\min, t_\max$ are chosen so that $E_\text{rel}(t)$
stays inside the kinematically allowed window
$[(\sqrt{E}-\sqrt{E'})^2,\,(\sqrt{E}+\sqrt{E'})^2]\cdot\tfrac{A}{k_B T}$.

### Numerical stability

Naïve evaluation of $\psi_n$ loses precision at large $x$ (i.e. $E \gg k_B T$).
UFG-Transport uses:

- a **fused `erfcx` evaluation** for $\psi_0,\psi_1$ to avoid overflow when the
  Gaussian factors are tiny;
- **tail truncation** of the $t$-integral at a user-specified number of
  standard deviations;
- quadrature orders that scale with $n$: default 32 GL points for the outer
  $t$-integral and 16 GL points for the azimuthal $\phi$-integral in
  $Q_n(x,t)$, which converged to machine precision in the closure tests.

## Configuration

The O&S kernel is opt-in via the `XsProcessorConfig` struct
(`include/xs/xs_processor.hpp`):

| Field | Default | Meaning |
|---|---|---|
| `use_resonance_kernel` | `false` | master switch |
| `resonance_cutoff_ev` | 1000 eV | upper energy for O&S |
| `resonance_A_min` | 10 | minimum AWR; lighter nuclides use free-gas/$S(\alpha,\beta)$ |
| `resonance_quad_t` | 32 | GL points for the target-speed integral |
| `resonance_quad_phi` | 16 | GL points for the azimuthal integral ($n \geq 2$) |

When enabled, the kernel replaces the cold-elastic contribution for each
nuclide with $A \geq A_\min$ below `resonance_cutoff_ev`. Above the cutoff,
the cold kernel is used; below the thermal cutoff, the free-gas or
$S(\alpha,\beta)$ kernels take over, exactly as in the default three-regime
layout.

## Verification

The implementation was verified against an independent reference
Doppler-broadened elastic cross section on U-238 in the resolved resonance
range. The original investigation is documented in
`coldcase/closed/resonance_kernel_closure.md`, which captures the final
kernel shape and the closed defect list. Briefly:

- the Legendre moments $\sigma_{s,n}(E \to E')$ match reference values to
  better than $10^{-6}$ relative error at the converged quadrature
  (`resonance_quad_t = 32`, `resonance_quad_phi = 16`);
- the resonance-integrated elastic cross section on U-238 at 600 K
  reproduces the NJOY `broadr`-produced Doppler-broadened value to the level
  expected from the UFG grid spacing;
- integration cost is approximately an order of magnitude above cold elastic
  per incident energy, which is acceptable because the kernel only runs
  inside the resonance energy range and only for nuclides with $A \geq 10$.

## Files

- `include/xs/resonance_kernel.hpp` — interface, $\psi_n$ / $Q_n$ declarations.
- `src/xs/resonance_kernel.cpp` — analytic $\psi_0, \psi_1$; numerical
  $\psi_{n\geq 2}$ via Gauss–Legendre quadrature.
- `include/xs/xs_processor.hpp` — configuration fields.
- `coldcase/closed/resonance_kernel_closure.md` — closure report (2026-03-29).

## References

- M. Ouisloumen and R. Sanchez, *A model for neutron scattering off heavy
  nuclei which accounts for the resonance structure*,
  *Nuclear Science and Engineering*, **107**, 189–200 (1991).
- W. Rothenstein, *Neutron scattering kernels in pronounced resonances for
  stochastic Doppler effect calculations*, *Annals of Nuclear Energy*, **23**,
  441–458 (1996) — alternative reference for the same kernel.
- R. E. MacFarlane, *NJOY: The Nuclear Data Processing System*
  (LA-UR-17-20093), BROADR module — baseline against which the O&S kernel is
  cross-checked.
