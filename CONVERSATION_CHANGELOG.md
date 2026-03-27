# Conversation Changelog

This file records the code and build-system changes completed during the current Codex collaboration session.

Date: 2026-03-26
Project: `DFT-FE`

## 0. Reorganization request recorded

### `REORGANIZATION_PLAN.md`
- Added a dedicated planning note for the requested reorganization direction.
- Recorded the requested container split:
  - MFEM types for data that may need MPI transfer
  - `ndarray` for small local arrays
- Recorded the requested implementation order:
  1. data-container and extensibility reorganization
  2. move PAW code into a subfolder
  3. restyle the code according to `STYLE`

## 1. Lifetime and ownership cleanup

### `src/DFTMesh.h`, `src/DFTMesh.cpp`
- Reorganized `DFTMesh` so mesh-related state is not copied accidentally.
- Made `DFTMesh` non-copyable and non-movable.
- Added a safer ownership path so `DFTMesh` can keep a `std::shared_ptr<const Structure>` when the structure lifetime must outlive mesh use.
- Moved heavy implementation out of the header into `DFTMesh.cpp`.

### `src/FEspace.h`, `src/FEspace.cpp`
- Reorganized `DFTGLLHexSpace` so it can safely depend on a shared `DFTMesh`.
- Made `DFTGLLHexSpace` non-copyable and non-movable.
- Changed it to keep stable mesh access instead of relying on fragile copied state.
- Added accessors and helpers for overlap/mass-related true-dof data:
  - mesh access
  - mass diagonal on true dofs
  - inverse square-root diagonal
  - helper application routines
- Moved heavy implementation out of the header into `FEspace.cpp`.

## 2. Build system and compile-speed cleanup

### `src/CMakeLists.txt`
- Created the reusable `dft_core` static library.
- Moved shared source compilation into `dft_core` so tests reuse compiled objects instead of recompiling the same implementation repeatedly.
- Added MFEM/tinyxml2 include paths centrally.
- Enabled a precompiled header for the core library.

### `src/pch.hpp`
- Added a precompiled-header entry point for the heavy MFEM include path.

### `test/CMakeLists.txt`
- Updated tests to link against `dft_core`.
- Kept test executables separate from the core library path.
- Left PCH off the test targets because GCC hit an internal compiler error when PCH was enabled there.

## 3. Overlap/mass path stabilization

### `src/FEspace.cpp`
- Changed mass assembly to use the true-dof system path.
- Kept overlap/mass diagonal extraction alive only while the local MFEM bilinear form was still valid.

### `test/TestOverlap.cpp`
- Reworked the overlap test so it checks:
  - true-dof size consistency
  - positivity of the mass diagonal
  - consistency of `M * (M^{-1/2})^2`
  - helper application routines

## 4. Kinetic operator implementation

### Removed old mixed Hamiltonian draft
- The previous `hamiltonianMatrix` design was not kept as the main implementation path.

### `src/kineticOperator.hpp`, `src/kineticOperator.cpp`
- Added a dedicated `KineticOperator` class.
- Implemented assembly of the true-dof kinetic matrix
  - using the diffusion bilinear form
  - corresponding to the Kohn-Sham kinetic term
- Implemented:
  - matrix storage
  - `Mult(...)`
  - quadratic-form energy evaluation
  - basic diagnostics/sanity checks
- Fixed a lifetime bug from MFEM `FormSystemMatrix(...)` by deep-copying the assembled sparse matrix before the local bilinear form is destroyed.

### `test/TestKinetic.cpp`
- Added a focused kinetic-term test.
- Checks:
  - matrix size matches true dofs
  - the periodic constant vector is near the null space
  - kinetic energy is non-negative

## 5. PAW XML parser

### `third_party/tinyxml2/tinyxml2.h`, `third_party/tinyxml2/tinyxml2.cpp`
- Added vendored `tinyxml2` to the repository for XML parsing.

### `src/PAW.h`, `src/PAW.cpp`
- Built a first PAW data/model layer around the XML setup file.
- Added support for:
  - general XML attribute maps
  - radial functions
  - PAW valence states
  - PAW channels
  - PAW setup storage and registry
- Implemented `LoadPAWSetupXML(...)`.
- The loader now reads:
  - top-level metadata blocks
  - named `radial_grid` definitions
  - named radial datasets
  - valence-state metadata
  - all-electron partial waves
  - pseudo partial waves
  - projector functions
- Explicit radial-grid references like `grid="log1"` are resolved through the named `<radial_grid id="...">` definition instead of assuming one global grid.
- `Dij` was discussed conceptually, but this XML file does not expose an explicit `Dij` dataset in the same way as the radial functions, so no full nonlocal PAW matrix was constructed yet.

### `test/TestPAW.cpp`
- Added a real parser test using `data/C.GGA_PBE-JTH.xml`.
- The test validates the parsed symbol, valence charge, metadata blocks, radial grids, radial datasets, states, and projectors.

## 6. Radial interpolation module

### `src/radialInterpolator.hpp`, `src/radialInterpolator.cpp`
- Added a reusable `RadialInterpolator` class for PAW tabulated radial functions.
- Supports:
  - construction from `RadialFunction`
  - construction from explicit `r` and `values`
  - grid validation
  - linear interpolation
  - end-point clamping

### `test/TestRadialInterpolator.cpp`
- Added a dedicated test using the parsed `zero_potential` radial function from the PAW XML.
- Checks:
  - exact reproduction of tabulated knot values
  - correct midpoint linear interpolation

## 7. Real spherical harmonics module

### `src/sphericalHarmonics.hpp`, `src/sphericalHarmonics.cpp`
- Added a standalone real spherical harmonics module.
- Implemented:
  - factorial helper
  - Cartesian-to-spherical angle conversion
  - associated Legendre polynomials
  - real spherical harmonics from angles
  - real spherical harmonics from Cartesian coordinates
- Current convention is the real harmonic basis with the standard Condon-Shortley phase inherited through the associated Legendre polynomial.

### `test/TestSphericalHarmonics.cpp`
- Added a dedicated regression test for known low-order values on coordinate axes.
- Verifies representative `s`, `p`, and `d` values:
  - `Y00`
  - `Y10`
  - `Y11`
  - `Y1,-1`
  - `Y20`

## 8. Atom-centered PAW basis evaluation

### `src/pawBasisEvaluator.hpp`, `src/pawBasisEvaluator.cpp`
- Added a first `PAWBasisEvaluator` module that combines:
  - atom center
  - radial interpolation
  - real spherical harmonics
- Supports evaluation of
  - `R(r) Y_lm(rhat)`
  - directly from a 3D point or from a centered displacement vector
- Returns zero outside the tabulated radial range.
- Handles the origin safely:
  - finite `s`-like value at the origin
  - zero for higher angular momentum at the origin

### `test/TestPAWBasisEvaluator.cpp`
- Added a dedicated evaluator test.
- Checks:
  - `s` value at the origin
  - `p_z` on the z-axis
  - `p_x` on the x-axis
  - zero outside the radial cutoff

## 9. Test coverage added during this session

The following tests were added or substantially improved:
- `test/TestOverlap.cpp`
- `test/TestKinetic.cpp`
- `test/TestPAW.cpp`
- `test/TestRadialInterpolator.cpp`
- `test/TestSphericalHarmonics.cpp`
- `test/TestPAWBasisEvaluator.cpp`

## 10. Current codebase direction after this session

The project now has a cleaner staged path:

1. stable mesh and FE-space ownership
2. overlap/mass true-dof support
3. kinetic operator
4. PAW XML parsing
5. radial interpolation
6. spherical harmonics
7. atom-centered basis evaluation

This means the next natural development step is to connect parsed PAW states/projectors directly to the 3D basis evaluator and then assemble projector overlaps or nonlocal PAW contributions on the FE basis.

## 11. Periodic KD-tree for grid-point radius search

### `third_party/nanoflann.hpp`, `third_party/nanoflann.COPYING`
- Vendored the header-only `nanoflann` library from GitHub.
- Added the upstream license file alongside the vendored header.

### `src/DFTMesh.h`
- Added lightweight accessors for the attached structure pointer and stored lattice, so periodic geometry can be reused by search utilities.

### `src/FEspace.h`, `src/FEspace.cpp`
- Added `TrueDofCoordinates()` to extract the true-dof grid-point coordinates from the periodic FE space.
- This uses a coordinate `GridFunction` and `GetTrueDofs(...)` to obtain one Cartesian point per scalar true dof.

### `src/periodicKDTree.hpp`, `src/periodicKDTree.cpp`
- Added a periodic KD-tree based radius-search utility for grid points.
- The implementation now follows the Python design direction discussed later in the session:
  - one KD-tree is built only on the base-cell wrapped points
  - periodic images are handled by looping over image shifts and performing equivalent shifted-query searches
  - hits are returned with both the original point index and the periodic image shift
- Added configurable per-axis periodic image depth:
  - for example `{1,1,1}`, `{1,0,0}`, `{2,1,0}`
- Removed the earlier merged/minimum-image result path and kept the image-resolved hit list only, since the PAW use case requires the image shift to survive.

### `test/TestPeriodicKDTree.cpp`
- Added a brute-force validation test for the periodic KD-tree.
- The test checks:
  - image-resolved periodic hits against a brute-force search
  - invariance under lattice translation
  - mixed-periodicity configuration such as `{1,0,0}`

## 12. PAW fixed nonlocal data

### `src/PAW.h`, `src/PAW.cpp`
- Extended `PAWSetup` to store state-resolved radial data:
  - all-electron partial waves by state
  - pseudo partial waves by state
  - projector functions by state
- Added fixed on-site matrices:
  - `kinetic_energy_differences`
  - `static_coulomb_correction`
  - `fixed_nonlocal_correction`
- Added channel-level matrix blocks:
  - `kinetic_correction`
  - `static_coulomb_correction`
  - `fixed_nonlocal_correction`

### Kinetic fixed term
- Parsed `<kinetic_energy_differences>` from the PAW XML as an explicit square matrix in state order.
- Mapped that global state matrix into per-channel sub-blocks.

### Static Coulomb fixed term
- Replaced the earlier local-potential approximation with a form based on the uploaded equations.
- Implemented the spherical `L=0` / `00` contribution only.
- Added helpers for:
  - trapezoidal radial quadrature weights
  - radial moments
  - Coulomb radial inner products
  - radial products and differences
  - `1/r` nuclear radial integral
- Implemented:
  - `Delta^a`
  - `Delta^a_{00,i1,i2}`
  - `g00` construction from the XML `shape_function`

### Shape function handling
- Implemented `shape_function type="sinc"` explicitly.
- Corrected the first version so it now:
  - uses plain `sinc`, not `sinc^2`
  - applies the XML cutoff radius `rc`
  - normalizes after cutoff truncation
- The current implementation uses only the spherical `Y_00 = 1/sqrt(4π)` component.

### `test/TestPAW.cpp`
- Extended the PAW parser test so it checks:
  - fixed-matrix sizes on the full state basis
  - fixed-matrix sizes on the per-channel blocks

## 13. Fixed PAW nonlocal Hamiltonian operator

### `src/pawNonlocalFixedOperator.hpp`, `src/pawNonlocalFixedOperator.cpp`
- Added a first operator for the fixed PAW nonlocal Hamiltonian term.
- The operator:
  - gathers nearby true-dof grid points around an atom with the periodic KD-tree
  - evaluates projector functions on those points
  - forms projector overlaps using the GLL diagonal mass
  - applies the fixed PAW matrix in projector space
  - scatters the resulting contribution back to the true-dof vector

### `src/CMakeLists.txt`
- Added `pawNonlocalFixedOperator.cpp` to the core library build.
- Added `third_party` include usage for the vendored KD-tree header.

### `test/TestPAWNonlocalFixed.cpp`
- Added a smoke test for the fixed PAW nonlocal operator.
- The test checks:
  - matrix/projector size consistency
  - that nearby grid points are found
  - output vector size consistency

## 14. Notes and limitations after the new work

- The periodic KD-tree now returns image-resolved hits only.
- The static Coulomb correction currently uses only the spherical `L=0` part of the uploaded equations.
- The fixed PAW nonlocal operator currently uses a first angular simplification with `m = 0` per state as a scaffold; it is not yet the full projector expansion over all magnetic components.
- Runtime execution of newly built test binaries remained unreliable in this sandbox because some produced artifacts were emitted as non-executable `data`, so several steps were compile-verified but not run here.
