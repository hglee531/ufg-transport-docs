# CPM — cylindrical geometry

The Collision Probability Method (CPM) discretises the integral form of the
transport equation by tabulating the probability that a neutron born isotropically
and uniformly in one region has its **first collision** in another. In cylindrical
geometry, these probabilities reduce to integrals of **Bickley–Naylor functions**
$\text{Ki}_n$ along chords across the pin.

UFG-Transport implements flat-source $P_0$ CPM for 1-D cylindrical pin-cell
problems. The full derivation is in `docs/solver_theory_and_algorithms.md` §8.

## Governing equation

Under the flat-source and isotropic-scattering assumptions, the scalar flux in
region $i$ satisfies

$$
\phi_i \;V_i\;\Sigma_{t,i}
\;=\;
\sum_{j} P_{j\to i}\;V_j\,Q_j,
$$

where $V_i$ is the region volume (per unit axial length), $Q_j$ is the total
source in region $j$ (scattering + fission + external), and $P_{j\to i}$ is the
first-flight collision probability.

## Collision probabilities via Bickley–Naylor

For an annular-ring discretisation with regions labelled by radial index $i$, the
key building block is the chord transmission. For a chord at impact parameter
$y$ and cumulative optical path $\tau$,

$$
\int_0^{\pi/2}\!e^{-\tau/\cos\theta}\,d\theta \;=\; \text{Ki}_3(\tau),
$$

the third-order Bickley–Naylor function. Region-to-region collision probabilities
are computed by integrating $\text{Ki}_3$ differences over the chord impact
parameters that cross each ring, and then dividing by region volumes.

UFG-Transport evaluates $\text{Ki}_n(x)$ directly by Gauss–Legendre quadrature on
its integral representation (`src/utils/bickley.cpp`); no lookup tables are used.

## Boundary conditions

- **Inner boundary** ($r=0$) — reflective by symmetry (no explicit condition is
  needed: the ring volume integrals handle the axis).
- **Outer boundary** — **white** (isotropic return), appropriate for a
  Wigner–Seitz pin-cell cylindrical approximation of a square pitch lattice.

The white BC is realised by a single extra collision-probability entry that
returns any neutron escaping the outer ring as an isotropic source distributed
over the outer ring by volume.

## Algorithm

```
compute P_{j→i} once per geometry + XS set (N_reg^2 entries)
repeat source iteration until ‖φ^(k+1) - φ^(k)‖ < tol:
    assemble Q_j = Σ_{g'} Σ_s^{g'→g} φ_{g',j} + χ_g Q_fis,j + Q_ext
    φ_i = (Σ_j P_{j→i} V_j Q_j) / (V_i Σ_{t,i})
```

Collision probabilities are built once per group per (ring radii, $\Sigma_{t,g}$)
tuple. At UFG resolution the $P_{j\to i}$ matrices dominate the compute; a dense
table is fine for the handful of rings typical of a pin cell.

## CLI invocation

```bash
./ufg_1d_app input.i --solver cpm --output-dir ../output/myrun
```

CPM has no angular-order flag; it is intrinsically $P_0$ in scatter and in spatial
source. Higher-order CPM (DP$_N$, or $P_1$ correction) is planned but not yet
implemented.

## Strengths and limitations

**Strong points.** Handles the pin cell geometry natively, with exact analytic
treatment of the radial streaming via $\text{Ki}_3$. No ray effect. Excellent for
resonance self-shielding studies when the UFG library has sharp peaks — the
collision probability picks up the resonance absorption in the fuel ring
accurately.

**Limitations.** Flat source per region means fine rings are required inside
regions with strong flux gradients (the fuel periphery, for instance). Isotropic
scatter only — anisotropic CPM would multiply the table by one more dimension per
Legendre order and is a known research topic.

## Cross-references

- `docs/solver_theory_and_algorithms.md` §8 — full CPM derivation with Bickley
  identities.
- `src/solver/solver_cpm.cpp`, `src/utils/bickley.cpp` — implementation.
