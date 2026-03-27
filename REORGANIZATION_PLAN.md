# Reorganization Plan

Date: 2026-03-27
Project: `DFT-FE`

This file records the requested reorganization direction before implementation.

## Goals

The codebase should be reorganized to make future MPI support and related feature
work more extensible.

## Requested Rules For Data Containers

1. Variables that may need to transfer between MPI processes should be defined
   using MFEM data structures.
2. Large vectors and other MPI-exchanged data should follow the MFEM-managed path.
3. Small local arrays that do not need MPI transfer should use the in-source
   `ndarray` utilities.

## Requested Reorganization Steps

1. Establish a clearer data-ownership split:
   - MFEM types for MPI-facing or potentially distributed data.
   - `ndarray` for small source-local arrays.
2. Move PAW-related files into a dedicated subfolder.
3. Restyle the codebase according to [STYLE](/home/dong/home/DFT-FE/STYLE).

## Notes For Implementation

- The container choice should be driven by whether the data is expected to cross
  MPI process boundaries now or in a near-future extension.
- The PAW subfolder move should preserve build integration and test coverage.
- Style refactoring is intended as a later step after the structural
  reorganization is in place.
