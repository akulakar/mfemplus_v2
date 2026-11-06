# CLAUDE.md — SolidMechanics

Guidance shared by all solid-mechanics modules under `mfemplus/SolidMechanics`. Module-specific guidance lives in
each module's own `CLAUDE.md` (e.g. `Fracture/CLAUDE.md`); build and library-wide notes are in the root `CLAUDE.md`.

## Layout
- `Fracture`: variational (phase-field) fracture integrators.
- `General`: strain/stress helpers (mostly empty).
- `solid_mechanics_utils.hpp`: internal shared Voigt kernels (namespace `mfemplus::detail`, included only by `.cpp` files).
- `mfemplus_solids.hpp`: umbrella header for the solid-mechanics modules.
- `Hyperelasticity`, `LinearElasticity`, `Viscoelasticity`: empty placeholders.

## Conventions
- Integrators take the displacement `GridFunction` + its `FiniteElementSpace` and read element dofs with
  `GetElementVDofs`.
- Strains/stresses are Voigt vectors (engineering shear): `Bu` strain, `CBu` stress, ordering (xx,yy,xy) in 2D and
  (xx,yy,zz,yz,xz,xy) in 3D. `planeApprox` 0 = plane strain, 1 = plane stress for 2D.
- Isotropic integrators never form `B` or `C`: they use C = λ I⊗I + 2μ I_sym in closed form
  (`detail::IsotropicLame`, `AddIsotropicStiffness`, `TwiceIsotropicStrainEnergy`, `IsotropicVolumetricStress`).
  Plane stress only changes λ to Eν/(1−ν²). Only the Anisotropic* integrators use the dense Voigt `B` and `C`.
- E and ν are evaluated only at the first integration point of an element (assumed constant per element).
- Use the shared kernels in `solid_mechanics_utils.hpp` (closed-form isotropic stiffness/energy, Voigt B matrix for
  the anisotropic path) instead of re-implementing them.
