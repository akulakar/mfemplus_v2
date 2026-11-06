//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Classes related to hyperbolic PDEs.
// Currently containts numerical flux functions.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/FluidMechanics/CompressibleFlow/hyperbolic.hpp"

#include <algorithm>
#include <cmath>

using namespace std;
using namespace mfem;

namespace mfemplus
{

   real_t HLLCFlux::Eval(const Vector &state1,
                         const Vector &state2,
                         const Vector &nor,
                         FaceElementTransformations &Tr,
                         Vector &flux) const
   {
      // Scaled normal fluxes F(u)n (scaled by |nor|)
      fluxFunction.ComputeFluxDotN(state1, nor, Tr, fluxN1);
      fluxFunction.ComputeFluxDotN(state2, nor, Tr, fluxN2);

      const real_t nmag = nor.Norml2();
      const int E = dim + 1; // index of the total energy

      // Primitive variables on both sides. Index 1 = left (u⁻), 2 = right (u⁺).
      // Density and pressure are floored so that the wave speeds and star states stay finite.
      const real_t rho1 = std::max(state1(0), small_dens);
      const real_t rho2 = std::max(state2(0), small_dens);
      real_t un1 = 0.0, un2 = 0.0, ke1 = 0.0, ke2 = 0.0;

      for (int d = 0; d < dim; d++)
      {
         const real_t m1 = state1(1 + d), m2 = state2(1 + d);
         un1 += m1 * nor(d);
         un2 += m2 * nor(d);
         ke1 += m1 * m1;
         ke2 += m2 * m2;
      }
      un1 /= rho1 * nmag;
      un2 /= rho2 * nmag;
      ke1 *= 0.5 / rho1; // ½ρ|u|²
      ke2 *= 0.5 / rho2;

      const real_t p1 = std::max((gamma - 1.0) * (state1(E) - ke1), small_pres);
      const real_t p2 = std::max((gamma - 1.0) * (state2(E) - ke2), small_pres);
      const real_t c1 = std::sqrt(gamma * p1 / rho1);
      const real_t c2 = std::sqrt(gamma * p2 / rho2);

      // Roe averages
      const real_t w1 = std::sqrt(rho1), w2 = std::sqrt(rho2);
      const real_t wsum = w1 + w2;
      const real_t H1 = (state1(E) + p1) / rho1;
      const real_t H2 = (state2(E) + p2) / rho2;
      const real_t Hroe = (w1 * H1 + w2 * H2) / wsum;

      real_t unroe = 0.0, u2roe = 0.0;
      for (int d = 0; d < dim; d++)
      {
         const real_t uroe = (state1(1 + d) / w1 + state2(1 + d) / w2) / wsum;
         unroe += uroe * nor(d);
         u2roe += uroe * uroe;
      }
      unroe /= nmag;

      const real_t croe = std::sqrt(std::max((gamma - 1.0) * (Hroe - 0.5 * u2roe), (real_t)0.0));

      // Wave speed estimates (Einfeldt)
      const real_t SL = std::min(un1 - c1, unroe - croe);
      const real_t SR = std::max(un2 + c2, unroe + croe);
      const real_t maxE = std::max(std::fabs(SL), std::fabs(SR));

      // Supersonic cases, pure upwinding
      if (SL >= 0.0)
      {
         flux = fluxN1;
         return maxE;
      }
      if (SR <= 0.0)
      {
         flux = fluxN2;
         return maxE;
      }

      // Contact wave speed
      const real_t dL = rho1 * (SL - un1);
      const real_t dR = rho2 * (SR - un2);
      const real_t Ss = (p2 - p1 + un1 * dL - un2 * dR) / (dL - dR);

      // F̂ n = F(u_K)n + |n| S_K (U*_K - u_K), K = L if S* ≥ 0 else R
      const bool left = (Ss >= 0.0);
      const Vector &state = left ? state1 : state2;
      const Vector &fluxN = left ? fluxN1 : fluxN2;
      const real_t rho = left ? rho1 : rho2;
      const real_t un = left ? un1 : un2;
      const real_t p = left ? p1 : p2;
      const real_t S = left ? SL : SR;

      const real_t fac = rho * (S - un) / (S - Ss);
      flux(0) = fluxN(0) + nmag * S * (fac - state(0));

      for (int d = 0; d < dim; d++)
      {
         const real_t ustar = state(1 + d) / rho + (Ss - un) * nor(d) / nmag;
         flux(1 + d) = fluxN(1 + d) + nmag * S * (fac * ustar - state(1 + d));
      }

      const real_t Estar = fac * (state(E) / rho +
                                  (Ss - un) * (Ss + p / (rho * (S - un))));

      flux(E) = fluxN(E) + nmag * S * (Estar - state(E));

      return maxE;
   }

} // namespace mfemplus
