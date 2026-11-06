//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Classes related to hyperbolic PDEs.
// Currently containts numerical flux functions.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_COMPRESSIBLEFLOW_HYPERBOLIC
#define MFEMPLUS_COMPRESSIBLEFLOW_HYPERBOLIC

#include "mfem.hpp"

namespace mfemplus
{

   /**
    * @brief HLLC (Harten-Lax-van Leer-Contact) approximate Riemann solver for the
    * compressible Euler equations, state u = (ρ, ρu, E).
    *
    * The Riemann fan is approximated by three waves (left, contact, right) with
    * speeds S_L < S* < S_R, so contact and shear waves are resolved, unlike in
    * the HLL or Rusanov fluxes. With the unit normal n̂ and uₙ = u⋅n̂:
    *   - S_L, S_R: Einfeldt (Roe-averaged) estimates
    *       S_L = min(uₙ⁻ - c⁻, ũₙ - c̃),   S_R = max(uₙ⁺ + c⁺, ũₙ + c̃)
    *   - S*: contact speed from the pressure/momentum jump conditions
    *   - F̂ n = F(u⁻)n                         if S_L ≥ 0
    *           F(u⁻)n + S_L(U*⁻ - u⁻)         if S_L < 0 ≤ S*
    *           F(u⁺)n + S_R(U*⁺ - u⁺)         if S* < 0 < S_R
    *           F(u⁺)n                         if S_R ≤ 0
    *
    * @note Requires a FluxFunction for the Euler equations with
    * num_equations = dim + 2 (e.g. mfem::EulerFlux) and the same specific heat
    * ratio γ that the FluxFunction uses.
    * @note Only Eval() is implemented. Grad(), Average() and AverageGrad() are not
    * overridden and abort if called, so use with explicit time integration.
    * @note The density and pressure used for the wave speeds and the star states
    * are floored at small_dens and small_pres, so states with non-positive
    * density or pressure do not produce NaN. The upwind fluxes F(u)n are still
    * evaluated from the unmodified states.
    */
   class HLLCFlux : public mfem::NumericalFlux
   {
   protected:
      const mfem::real_t gamma;      // specific heat ratio, γ
      const int dim;                 // spatial dimension
      const mfem::real_t small_dens; // density floor
      const mfem::real_t small_pres; // pressure floor

      mutable mfem::Vector fluxN1, fluxN2;

   public:
      /**
       * @brief Constructor
       * @param fluxFunction Euler flux function F(u,x)
       * @param specific_heat_ratio specific heat ratio γ (must match the one used
       * by fluxFunction, EulerFlux does not expose it)
       * @param small_dens floor for the density in the wave speed and star state
       * computation
       * @param small_pres floor for the pressure in the wave speed and star state
       * computation
       */
      HLLCFlux(const mfem::FluxFunction &fluxFunction, mfem::real_t specific_heat_ratio,
               mfem::real_t small_dens = 1e-12, mfem::real_t small_pres = 1e-12)
          : NumericalFlux(fluxFunction), gamma(specific_heat_ratio), dim(fluxFunction.dim),
            small_dens(small_dens), small_pres(small_pres)
      {
         MFEM_VERIFY(fluxFunction.num_equations == dim + 2, "HLLCFlux requires the Euler equations, num_equations = dim + 2.");

         fluxN1.SetSize(fluxFunction.num_equations);
         fluxN2.SetSize(fluxFunction.num_equations);
      }

      /**
       * @brief Normal numerical flux F̂(u⁻,u⁺,x) n
       *
       * @param[in] state1 state (u⁻) from the first element (dim + 2)
       * @param[in] state2 state (u⁺) from the second element (dim + 2)
       * @param[in] nor normal vector (not a unit vector) (dim)
       * @param[in] Tr face element transformation
       * @param[out] flux HLLC flux, scaled by |nor| like F(u)n
       * @return max(|S_L|, |S_R|), the maximum wave speed estimate
       */
      mfem::real_t Eval(const mfem::Vector &state1, const mfem::Vector &state2,
                        const mfem::Vector &nor,
                        mfem::FaceElementTransformations &Tr,
                        mfem::Vector &flux) const override;
   };

} // namespace mfemplus

#endif // MFEMPLUS_COMPRESSIBLEFLOW_HYPERBOLIC
