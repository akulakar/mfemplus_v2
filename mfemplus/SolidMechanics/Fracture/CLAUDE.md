# CLAUDE.md — Fracture module

Guidance for `mfemplus/SolidMechanics/Fracture` (variational / phase-field fracture integrators, library target
`mfemplus::fracture`). Build notes are in the root `CLAUDE.md`; shared solid-mechanics conventions (Voigt, isotropic kernels) are in
`mfemplus/SolidMechanics/CLAUDE.md`.

## Layout
- `bilininteg.*`: BilinearFormIntegrators (elasticity with degraded stiffness, and damage-equation
  terms of the form `∫ 2W(u) d φ`).
- `lininteg.*`: matching LinearFormIntegrators (damage RHS, body forces/eigenstrain).
- `energy.*`: strain energy evaluators and the shared `QuadratureFunction` path.
- `operators.*` (empty), `fracture_utils.hpp` (internal helpers, namespace `mfemplus::detail`, included only by the `.cpp` files).

## Conventions
- The integrator's own element matrix/vector is for the scalar damage space (the displacement `GridFunction` is read
  via `GetElementVDofs`, see the SolidMechanics guidance).
- Integrators assemble 2W, not W (factor absorbed by the caller's coefficients).
- The `...VolumetricEigenStrain...` pair (`bilininteg.cpp` and `lininteg.cpp`) must stay in sync: both
  assemble 2W = (ε−ε*)ᵀC(ε−ε*) with ε* = α·I, i.e. εᵀCε − 2α(CI)·ε + α²IᵀCI. The quadratic term is
  `0.5 * α * InnerProduct(C·I·α, I)`; `volumetric_stress` already contains one factor of α.
- Damage forms depend on u, so the element matrix can be exactly zero at the first assembly; build such
  `BilinearForm`s with `Assemble(0)`/`Finalize(0)` (skip_zeros=0) or MFEM throws
  "Entry for column N is not allocated" on later reassembly.
- The two VolumetricEigenStrain integrators apply a Gershgorin-type check on (stress - volumetric stress)
  (`detail::EffectiveStressBoundIsNegative`, in `fracture_utils.hpp`); if it fails, the energy is scaled by 1e-20
  (effectively dropped). The check is on by default and can be disabled with the trailing `use_gershgorin_check`
  constructor argument (matrix, LF and evaluator must be given the same value; with a shared `QuadratureFunction`
  this is automatic). The Pressure integrators/evaluator support the same trailing `use_gershgorin_check` argument
  (checked against the pressure vector), but it is OFF by default there. The other integrators have no such check.
- Use the shared kernels in `fracture_utils.hpp` (Gershgorin check, fracture helpers) and `solid_mechanics_utils.hpp`
  instead of re-implementing them.
- Quadrature default is order 2p + `OrderGrad` for all fracture integrators (bilinear and linear forms must match).
- Damage-equation energy 2W lives in `energy.*` (`Isotropic...Evaluator` classes, one per variant: plain,
  Pressure, VolumetricEigenStrain). The matching bilinear and LF integrators each have two constructors: the usual one
  (E, ν, u, ...) evaluates 2W on the fly through an evaluator; the other takes a `QuadratureFunction` filled by
  `Evaluator::Evaluate(W)` (call it after every displacement update) and integrates with the rule of its
  `QuadratureSpace`, so matrix and RHS share one energy and one rule. The matrix and LF therefore always use the same
  2D formulas too.
