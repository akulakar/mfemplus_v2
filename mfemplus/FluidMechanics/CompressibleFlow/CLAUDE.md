# CLAUDE.md — CompressibleFlow module

Guidance for `mfemplus/FluidMechanics/CompressibleFlow` (compressible Euler numerical fluxes for MFEM's hyperbolic
framework). Build notes are in the root `CLAUDE.md`.

## Layout
- `hyperbolic.{hpp,cpp}`: numerical fluxes (namespace `mfemplus`), derived from `mfem::NumericalFlux`. Currently only
  `HLLCFlux`.
- `reference_files/`: read-only third-party references, do not edit and do not build.
  - `hyperbolic.{hpp,cpp}`: MFEM's own hyperbolic module (`FluxFunction`, `NumericalFlux`, `RusanovFlux`, `EulerFlux`,
    `HyperbolicFormIntegrator`).
  - `Godunov_3D.cpp`, `Hydro_cmpflx.H`, `Hydro_riemann.H`: CAMR (AMReX) finite-volume scheme with PLM/PPM
    reconstruction, CTU transverse terms and an acoustic-type Riemann solver (`riemannus`).

## Conventions
- `HLLCFlux` is for the Euler equations, state `(ρ, ρu, E)`, `num_equations = dim + 2`. It takes the `FluxFunction`
  (e.g. `mfem::EulerFlux`) plus γ explicitly, because `EulerFlux` keeps γ private; the two must match.
- Only `Eval()` is implemented. `Grad`, `Average` and `AverageGrad` are not overridden and abort, so use it with
  explicit time integration.
- `nor` is not a unit normal. Fluxes are scaled by `|nor|`, like `FluxFunction::ComputeFluxDotN`. The returned speed
  is unscaled.
- Wave speeds are Einfeldt (Roe-averaged). ρ and p are floored (`small_dens`, `small_pres`, default 1e-12) only in the
  wave speed and star state computation. The upwind fluxes `F(u)n` use the raw states.
- `fluxN1`/`fluxN2` are `mutable` members, so an instance is not thread-safe (the `MFEM_THREAD_SAFE` branches of
  `RusanovFlux` were not carried over).
- The intended use is DG through `HyperbolicFormIntegrator`. The CAMR reconstruction and transverse-correction code
  is finite-volume specific and is not ported.
- Not implemented (ideas from CAMR): passive scalars advected by the sign of S*, pressure-based S_L/S_R estimates,
  general EOS.
- Code in this folder uses 3-space indentation inside `namespace mfemplus` and `using namespace` in the `.cpp`
  (the header still qualifies with `mfem::`).
- Status: `HLLCFlux` has never been compiled or tested. The module is not in CMake yet.
