# ACE and HDF5 nuclear data

UFG-Transport ingests nuclear data in the **OpenMC HDF5 format**, which is a modern,
self-describing repackaging of the classical ACE (A Compact ENDF) format.
This page summarises what UFG-Transport reads and why the format matters.

For the full, byte-level references (including every reaction MT, angular-distribution
encoding, and URR table layout) see:

- `REPORT_ACE_FORMAT.md` — ASCII ACE structure, NXS/JXS/XSS arrays, scattering laws.
- `REPORT_OPENMC_HDF5_FORMAT.md` — the HDF5 equivalent used by OpenMC ≥ 0.13.

## Why point-wise data?

Nuclear evaluations (ENDF/B-VIII, JEFF-3.3, JENDL-5) publish resonance parameters and
tabulated cross sections from ~10⁻⁵ eV to 20 MeV. Pre-processors (NJOY, OpenMC's
`openmc.data`) Doppler-broaden these evaluations at a fixed temperature and emit
**point-wise** cross sections on a unionised energy grid.

| Method | Typical groups | Self-shielding | Accuracy |
|---|---|---|---|
| Few-group | 2–8 | complex | low–medium |
| Multi-group | 69–2000 | required | medium |
| **Ultra-fine (UFG)** | **10k–100k** | **minimal** | **high** |
| Point-wise Monte Carlo | continuous | none | reference |

UFG-Transport sits in the third row: it uses a group structure fine enough that
resonance peaks are resolved on the grid itself, so explicit self-shielding models
(Bondarenko, subgroup, CENTRM) are largely unnecessary.

## What the reader extracts

For each nuclide at each available temperature the `ufg::ace::AceReader` retrieves:

- **Energy grid** — the unionised energy points $E_k$, in eV.
- **Principal reactions** — total, elastic scatter, absorption, capture, fission,
  plus any discrete-level inelastics (MT 51–91) present in the evaluation.
- **Reaction Q-values and thresholds** — needed for inelastic kinematics.
- **Angular distributions** — per-reaction tabular $\mu,\,\text{PDF},\,\text{CDF}$
  over centre-of-mass angle.
- **Energy distributions** — fission spectra (Maxwell, Watt, continuous tabular)
  and continuous secondary distributions for inelastics.
- **Fission data** — $\nu$ (total or prompt) and the emission spectrum χ.
- **Unresolved resonance probability tables** — Section 33/URR, used for
  p-table sampling in the UFG processor.
- **Thermal scattering** — S(α,β) kernels for bound scatterers
  (`c_H_in_H2O`, `c_C_in_graphite`, …), evaluated as incoherent inelastic scattering
  below a nuclide-dependent cutoff (typically 4–5 eV).

## HDF5 vs ASCII ACE

The reader targets HDF5 rather than raw ASCII ACE because:

| Feature | ASCII ACE | OpenMC HDF5 |
|---|---|---|
| File layout | flat XSS with JXS indexing | named groups/datasets |
| Type identification | implicit JXS block numbers | explicit `type` attribute |
| Multi-temperature | one file per $T$ | all $T$ in one file |
| Energy units | MeV | eV |
| Angular distributions | raw 32-bin or tabulated | all tabulated $\mu,\,\text{PDF},\,\text{CDF}$ |
| Integer indexing | 1-based Fortran | explicit offsets |

HDF5 eliminates most of the manual index arithmetic that plagues ACE readers; reaction
data lives under paths like `/U235/reactions/reaction_018/...` with self-describing
attributes.

## Library layout expected by UFG-Transport

```
xsdata/
  endfb-viii.1-hdf5/
    cross_sections.xml        # index of nuclides → file paths
    neutron/
      H1.h5
      O16.h5
      U235.h5
      ...
    thermal/
      c_H_in_H2O.h5
      c_C_in_graphite.h5
      ...
```

The `cross_sections.xml` index is the OpenMC standard. UFG-Transport parses it to
map nuclide names (`U235`, `H1`) to file paths, and looks up thermal tables by the
`S(α,β) table` name declared in the material input.

## Temperature handling

A single HDF5 file typically contains data at several temperatures
($T \in \{294,\,600,\,900,\,\ldots\}$ K). UFG-Transport does not interpolate the raw
ACE grids directly — that would require re-resolving every resonance. Instead it
processes the UFG library at each bracketing ACE temperature, then **linearly
interpolates the resulting UFG cross sections in $\tau = T^{1/3}$**. See
[Multi-level condensation](multilevel.md) for the cube-root-of-T formulation and
its implementation.

## Further reading

- `REPORT_ACE_FORMAT.md` — complete ACE field-by-field reference (≈20 pages).
- `REPORT_OPENMC_HDF5_FORMAT.md` — HDF5 path reference with every attribute.
- [OpenMC nuclear data documentation](https://docs.openmc.org/en/stable/io_formats/data.html).
- MacFarlane & Muir, *The NJOY Nuclear Data Processing System*, LA-UR-17-20093.
