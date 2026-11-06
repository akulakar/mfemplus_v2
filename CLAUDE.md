# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status / build
mfemplus is an add-on library of MFEM integrators (namespace `mfemplus`, C++17), meant to be compiled into
other MFEM projects. There are no tests or lint config. `README.md` is empty.
Dependencies: MFEM (`#include "mfem.hpp"`), MPI/hypre/metis through MFEM.
Build: the root `CMakeLists.txt` finds MFEM/MPI (unless a `mfem` target already exists) and adds `mfemplus/`,
which defines one library target per module (`mfemplus::fracture`). A consuming project does
`add_subdirectory(mfemplus-repo)` and links `mfemplus::fracture`; the include root is the repository root.
Headers and sources live together under `mfemplus/`. Include by path:
`#include "mfemplus/SolidMechanics/Fracture/bilininteg.hpp"` (umbrella: `mfemplus.hpp` at the repo root).
`fracture_utils.hpp` and `solid_mechanics_utils.hpp` are internal (included only by the `.cpp` files).
Headers contain no `using namespace`; qualify with `mfem::`.
The installed MFEM CMake config in this environment fails to import (`adios2::adios2` target missing), so a
CMake build was only checked with a stub `mfem` target; direct `mpicxx` builds with `-I<repo root>` work.

## Layout
- `mfemplus/SolidMechanics`: solid-mechanics modules (Fracture, General, and empty Hyperelasticity /
  LinearElasticity / Viscoelasticity placeholders); see `mfemplus/SolidMechanics/CLAUDE.md` for shared conventions
  (Voigt, isotropic kernels) and `mfemplus/SolidMechanics/Fracture/CLAUDE.md` for the phase-field fracture module.
- `mfemplus/FluidMechanics/CompressibleFlow`: `HLLCFlux` numerical flux for the Euler equations (`hyperbolic.*`),
  plus `reference_files/` (MFEM and CAMR sources, read-only); see `mfemplus/FluidMechanics/CompressibleFlow/CLAUDE.md`.
  It is not yet in CMake (`add_subdirectory(FluidMechanics)` is commented out) and has not been compiled.
- `config`, `documentation` (repo root): empty placeholders.

## Shell commands
When the user runs a command with the `!` prefix (it shows up as `<bash-input>` / `<bash-stdout>`),
don't respond to it. Only comment on the output if the user asks a question about it.