# CLAUDE.md — CompressibleFlow module

Guidance for `mfemplus/FluidMechanics/CompressibleFlow` (compressible Euler numerical fluxes for MFEM's hyperbolic
framework). Build notes are in the root `CLAUDE.md`.

## Layout
- `flux.{hpp,cpp}`: numerical fluxes (namespace `mfemplus`), derived from `mfem::NumericalFlux`. `HLLCFlux` and `HLLFlux` (HLL,
  same Einfeldt wave speeds and conventions, shared helper `ComputeEulerWaves`).
- `limiter.{hpp,cpp}`: limiters for the DG solution. `PositivityLimiter` (Zhang–Shu) and
  `BarthJespersenLimiter` (slope limiter on face-neighbor means).
- `CMakeLists.txt`: library target `mfemplus::compressible_flow` (`flux.cpp`, `limiter.cpp`).
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
- `PositivityLimiter(fes, gamma, small_dens, small_pres, int_order_offset = 1)` with `int Apply(ParGridFunction &u)`.
  Per element it keeps the mean and scales the solution about it by the smallest θ ∈ [0,1] that gives ρ ≥ floor and
  p ≥ floor at the volume quadrature points and the face quadrature points (order 2p + `int_order_offset`).
  The density step is followed by an exact quadratic solve for the pressure.
  - Element-local, no MPI communication, no cached or `mutable` state.
  - Needs a nodal partition-of-unity basis (checked with `MFEM_VERIFY`).
  - If an element mean is already below a floor, θ is set to 0 (the element becomes its mean) and the element is
    counted in the return value, which is local to the rank. The mean is never floored, to keep conservation.
  - Positivity is guaranteed only at the checked points.
  - Supported geometries: segment, triangle, square, tetrahedron, cube. The face points are mapped from the
    reference face vertices, so no face transformations are needed.
- `BarthJespersenLimiter(pfes, int_order_offset = 1)` with `int Apply(ParGridFunction &u)` (non-const). Any vdim,
  acts on the conserved variables with a single θ per element (min over components and over the volume and face
  quadrature points) from the min/max of the own and face-neighbor element means; scales about the mean. The
  local and shared (MPI) face-neighbor lists are built once in the constructor (static, conforming mesh), and
  `Apply` communicates the neighbor means through `ExchangeFaceNbrData`. Returns the number of limited elements.
  Does not enforce positivity, so apply `PositivityLimiter` afterwards.
- Status: `HLLCFlux`, `HLLFlux`, `PositivityLimiter` and `BarthJespersenLimiter` have never been run. `flux.cpp` and
  `limiter.cpp` pass `mpicxx -fsyntax-only -Wall -Wextra` (they need `-I/opt/homebrew/opt/adios2/include` besides the
  MFEM, hypre and metis includes).
