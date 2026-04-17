# Install & build

UFG-Transport builds with CMake and produces four executables. The instructions
below cover a clean build on Linux and Windows.

## Requirements

- **CMake** ≥ 3.24
- **C++17** compiler (GCC ≥ 9, Clang ≥ 11, MSVC ≥ 2019)
- **HDF5** with C API (1.14 tested)
- **OpenMP** (optional, enables CPU parallelism)
- Internet access during first configure (GoogleTest is fetched automatically)

You also need an **OpenMC-style HDF5 cross-section library** with a
`cross_sections.xml` index. ENDF/B-VIII.0 or later is recommended. The library
layout expected by the code is

```
xsdata/
  endfb-viii.1-hdf5/
    cross_sections.xml
    neutron/
      H1.h5
      O16.h5
      U235.h5
      ...
    thermal/
      c_H_in_H2O.h5
      ...
```

## Configure

From the repository root:

=== "Linux / macOS"

    ```bash
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    ```

=== "Windows (PowerShell)"

    ```powershell
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release `
          -DHDF5_ROOT="C:/HDF5/1.14.3"
    ```

If HDF5 is not found automatically, pass `-DHDF5_ROOT=<prefix>`. On Windows you
will also need to copy the HDF5 DLLs into the runtime directory (see the troubleshooting
box below).

## Build

=== "Multi-config (Visual Studio)"

    ```powershell
    cmake --build build --config Release
    ```

=== "Single-config (Ninja, Makefiles)"

    ```bash
    cmake --build build -j
    ```

## Test

```bash
ctest --test-dir build -C Release
```

(omit `-C Release` for single-config generators.)

## Executables produced

| Executable | Purpose | Typical input |
|---|---|---|
| `ufg_transport_app` | Single-nuclide verification driver | one nuclide HDF5 |
| `ufg_0d_app` | 0-D homogeneous transport | one `.i` file |
| `ufg_1d_app` | 1-D slab or cylindrical transport | one `.i` file |
| `test_fg_params` | Free-gas parameter comparison utility | developer use |

The executables land under `build/Release/` (Visual Studio) or `build/` (Ninja).

!!! warning "Windows DLL copy"

    On Windows, HDF5 is shipped as DLLs that must be discoverable at run time.
    After the first build, copy the DLLs from the HDF5 `bin/` directory into
    `build/Release/` (or add the HDF5 `bin/` to `PATH`). Expected DLLs include
    `hdf5.dll`, `hdf5_cpp.dll`, `zlib.dll`, and `szip.dll` depending on the
    HDF5 build you chose.

!!! tip "OpenMP threads"

    The code uses **4 threads** by default. Pass `--threads N` to `ufg_1d_app`
    or `ufg_0d_app` to override. The thread count is exposed for the XS processor
    and solver loops; above a certain count returns diminish because the scattering
    matrix construction becomes memory-bound.

## Optional: CUDA

The CMakeLists has a CUDA block that is currently commented out — the GPU kernels
are deferred until the CPU reference implementation is fully validated. Re-enable
the CUDA block only if you plan to work on the GPU port.

## Next

- [Running the code](running.md) walks through the PWR pin-cell multi-level demo.
- [Pin-cell verification](../demos/pincell_verification.md) shows the reference
  comparison plots.
