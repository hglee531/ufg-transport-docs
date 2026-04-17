# Running the code

This page walks through the three drivers and the HIT input format, using the
PWR pin-cell multi-level benchmark as the worked example.

## Input file format (HIT)

Input files are hierarchical `key = value` blocks. A minimal file contains six
top-level blocks:

```text
[CrossSectionLib]    # path to the HDF5 library
[GroupStructure]     # group boundaries (single- or multi-level)
[Materials]          # nuclide mixtures
[Geometry]           # infinite, slab, or cylindrical
[FixedSource]        # source spectrum and spatial profile
[Solve]              # solver selection and tolerances
```

### Path resolution

Paths inside `[CrossSectionLib]` may start with `@`, which means *resolve relative
to the input file*:

```text
[CrossSectionLib]
  directory = '@../xsdata/endfb-viii.1-hdf5'
[]
```

This lets you launch runs from any shell working directory.

### Example: PWR pin-cell, multi-level

`examples/pwr_pincell_multilevel_sn.i`:

```text
# ── Cross-section library ──────────────────────────────────
[CrossSectionLib]
  directory = '@../xsdata/endfb-viii.1-hdf5'
[]

# ── Multi-level group structure ────────────────────────────
[GroupStructure]
  [hyperfine]
    type          = hfg
    energy_bin    = '20.0 1.0 1.0e-3 1.0e-5 0.5e-6 1.0e-8 1.0e-11'
    energy_groups = '500 1500 1500 1500 250 50'
  []
  [ultrafine]
    type          = ufg
    energy_bin    = '20.0 1.0 1.0e-3 1.0e-5 0.5e-6 1.0e-8 1.0e-11'
    energy_groups = '100 250 200 200 150 50'
  []
  [multigroup]
    type    = mg
    library = CASMO-70
  []
[]
```

See [theory/multilevel.md](../theory/multilevel.md) for the full multi-level flow
and the meaning of HFG/UFG/MG.

## The three drivers

### `ufg_transport_app` — single-nuclide verification

Reads one nuclide HDF5 file and prints group-averaged cross sections. Use this to
verify that the ACE reader and XS processor produce sensible numbers for a given
isotope at a given temperature.

```bash
./build/Release/ufg_transport_app nuclide.h5 [temperature_K] [num_groups] [L_max] [threads]
```

Example:

```bash
./build/Release/ufg_transport_app xsdata/endfb-viii.1-hdf5/neutron/U235.h5 293 2000
```

### `ufg_0d_app` — 0-D homogeneous

Parses a `.i` input, builds a homogenised multi-group library, and solves a 0-D
fixed-source problem by symmetric Gauss–Seidel. Output CSVs land in the current
working directory:

```
<stem>_xs.csv        # group-averaged macroscopic XS
<stem>_flux.csv      # converged scalar flux
<stem>_scat_l0.csv   # P0 scattering matrix
```

Example:

```bash
./build/Release/ufg_0d_app examples/uo2_mixed.i --verbose
```

k-eigenvalue mode is not yet implemented; only fixed-source problems run.

### `ufg_1d_app` — 1-D transport

The main executable. Parses a `.i` input, dispatches the chosen solver, writes
per-run CSVs into `--output-dir`.

```bash
./build/Release/ufg_1d_app <input.i> [options]
```

Options:

| Flag | Values | Meaning |
|---|---|---|
| `--solver` | `sn` / `pn` / `cpm` | transport solver |
| `--sn-order` | 2, 4, 8, 16 | angular quadrature order |
| `--pn-order` | 1, 3, 5, 7 | $P_N$ expansion order |
| `--pn-method` | `full` / `spn` | cylindrical $P_N$ or simplified $P_N$ |
| `--legendre` | 0, 1, 2, … | scattering-source expansion order |
| `--stage1` | `0d` / `cpm` / `homogenized` | HFG weighting method in multi-level |
| `--tolerance` | e.g. `1e-6` | source-iteration tolerance |
| `--max-iters` | integer | source-iteration cap |
| `--output-dir` | path | per-run output directory **(always use this)** |
| `--threads` | integer | OpenMP threads (default 4) |
| `--verbose` | flag | extra logging |

!!! warning "Always set `--output-dir`"

    The plotting scripts expect CSVs under `output/<run-name>/`. If you omit
    `--output-dir`, files land in the current working directory and the plot
    scripts will not find them.

### Auto solver selection

If `--solver` is omitted, `ufg_1d_app` picks:

- `cpm` for cylindrical geometry,
- `sn` for slab geometry.

The `[Solve]` block in the input file can also set the default; command-line flags override it.

## Worked example: pin-cell multi-level with Sn

```bash
cd build/Release
./ufg_1d_app ../../examples/pwr_pincell_multilevel_sn.i \
             --solver sn --sn-order 8 \
             --output-dir ../../output/pincell_multilevel_sn
```

On completion, inspect the output directory:

```
output/pincell_multilevel_sn/
  flux_moments.csv
  mg_xs.csv
  hfg_flux_<composition>.csv
  scat_matrix_<composition>.csv
  log.txt
```

Then plot with the standard script:

```bash
python scripts/plot_flux_spectra.py output/pincell_multilevel_sn
```

See [Pin-cell verification](../demos/pincell_verification.md) for the reference
plots and for comparison against OpenMC.

## Flux normalisation (mandatory)

Before comparing solver flux against OpenMC, read
`docs/normalization_protocol.md`. In short:

- **Solver** `phi_l0` is already volume-averaged → convert to per-lethargy by
  dividing by $\Delta u$.
- **OpenMC** `flux_per_u_cellN` is volume-integrated → divide by $V_\text{region}$
  to get volume-averaged.
- Match the two with a single overall constant
  $C = \sum(\text{solver}_\text{fuel} \cdot \Delta u)
       / \sum(\text{openmc}_\text{fuel} \cdot \Delta u)$
  and apply $C$ to OpenMC everywhere.

All comparison scripts (`compare_vs_openmc.py`, `plot_pincell_verification.py`,
`plot_cpm_verification.py`) import from `flux_compare.py` — do not re-implement
normalisation inline.
