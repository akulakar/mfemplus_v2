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
- `mfemplus/SolidMechanics/Fracture`: variational (phase-field) fracture integrators.
  - `bilininteg.*`: BilinearFormIntegrators (elasticity with degraded stiffness, and damage-equation
    terms of the form `∫ 2W(u) d φ`).
  - `lininteg.*`: matching LinearFormIntegrators (damage RHS, body forces/eigenstrain).
  - `energy.*`: strain energy evaluators and the shared `QuadratureFunction` path.
  - `operators.*` (empty), `fracture_utils.hpp` (internal helpers).
- `mfemplus/SolidMechanics/solid_mechanics_utils.hpp` (internal shared Voigt kernels),
  `mfemplus/SolidMechanics/General` (strain/stress helpers; mostly empty).
- `mfemplus/FluidMechanics`, `mfemplus/SolidMechanics/{Hyperelasticity,LinearElasticity,Viscoelasticity}`,
  `config`, `documentation` (repo root): empty placeholders.

## Conventions that span files
- Integrators take the displacement `GridFunction` + its `FiniteElementSpace` and read element dofs with
  `GetElementVDofs`; the integrator's own element matrix/vector is for the scalar damage space.
- Strains/stresses are Voigt vectors (engineering shear): `Bu` strain, `CBu` stress, ordering (xx,yy,xy) in 2D and
  (xx,yy,zz,yz,xz,xy) in 3D. `planeApprox` 0 = plane strain, 1 = plane stress for 2D.
- Isotropic integrators never form `B` or `C`: they use C = λ I⊗I + 2μ I_sym in closed form
  (`detail::IsotropicLame`, `AddIsotropicStiffness`, `TwiceIsotropicStrainEnergy`, `IsotropicVolumetricStress`).
  Plane stress only changes λ to Eν/(1−ν²). Only the Anisotropic* integrators use the dense Voigt `B` and `C`.
- E and ν are evaluated only at the first integration point of an element (assumed constant per element).
- The `...VolumetricEigenStrain...` pair (`bilininteg.cpp` and `lininteg.cpp`) must stay in sync: both
  assemble 2W = (ε−ε*)ᵀC(ε−ε*) with ε* = α·I, i.e. εᵀCε − 2α(CI)·ε + α²IᵀCI. The quadratic term is
  `0.5 * α * InnerProduct(C·I·α, I)`; `volumetric_stress` already contains one factor of α.
- Integrators assemble 2W, not W (factor absorbed by the caller's coefficients).
- Damage forms depend on u, so the element matrix can be exactly zero at the first assembly; build such
  `BilinearForm`s with `Assemble(0)`/`Finalize(0)` (skip_zeros=0) or MFEM throws
  "Entry for column N is not allocated" on later reassembly.
- The two VolumetricEigenStrain integrators apply a Gershgorin-type check on (stress - volumetric stress)
  (`detail::EffectiveStressBoundIsNegative`); if it fails, the energy is scaled by 1e-20 (effectively dropped).
  The check is on by default and can be disabled with the trailing `use_gershgorin_check` constructor argument
  (matrix, LF and evaluator must be given the same value; with a shared `QuadratureFunction` this is automatic).
  The other integrators have no such check.
- `fracture_utils.hpp` and `solid_mechanics_utils.hpp` (namespace `mfemplus::detail`) hold the shared kernels (closed-form isotropic stiffness/
  energy, Voigt B matrix for the anisotropic path, Gershgorin check); the integrators should use it instead of
  re-implementing them.
- Quadrature default is order 2p + `OrderGrad` for all fracture integrators (bilinear and linear forms must match).
- Damage-equation energy 2W lives in `Fracture/energy.*` (`Isotropic...Evaluator` classes, one per variant: plain,
  Pressure, VolumetricEigenStrain). The matching bilinear and LF integrators each have two constructors: the usual one
  (E, ν, u, ...) evaluates 2W on the fly through an evaluator; the other takes a `QuadratureFunction` filled by
  `Evaluator::Evaluate(W)` (call it after every displacement update) and integrates with the rule of its
  `QuadratureSpace`, so matrix and RHS share one energy and one rule. The matrix and LF therefore always use the same
  2D formulas too.
