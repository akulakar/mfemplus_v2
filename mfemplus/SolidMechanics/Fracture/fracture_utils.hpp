//------------------------------------------------------------------------------------------------------------
// Internal helpers shared by the fracture bilinear and linear form integrators (Voigt notation kernels).
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_UTILS
#define MFEMPLUS_FRACTURE_UTILS

#include "mfem.hpp"
#include "mfemplus/SolidMechanics/solid_mechanics_utils.hpp"
#include <algorithm>
#include <cmath>


namespace mfemplus
{
    namespace detail
    {
        /// Degradation function (1 - d)^2 + k_epsilon.
        inline mfem::real_t Degradation(mfem::real_t damage, mfem::real_t k_epsilon)
        {
            return (1.0 - damage) * (1.0 - damage) + k_epsilon;
        }

        /// Gershgorin-type lower bound on the eigenvalues of the effective stress (stress - volumetric_stress), in
        /// Voigt form. Returns true if the bound is negative, i.e. the point is not guaranteed to be in tension
        /// and the damage driving energy is switched off by the eigenstrain integrators.
        /// In 2D the exact smallest principal value is used.
        inline bool EffectiveStressBoundIsNegative(int dim, const mfem::Vector &stress, const mfem::Vector &volumetric_stress)
        {
            if (dim == 2)
            {
                const mfem::real_t mean = (stress(0) - volumetric_stress(0) + stress(1) - volumetric_stress(1)) / 2.0;
                const mfem::real_t radius = std::sqrt(std::pow((stress(0) - stress(1)) / 2.0, 2.0) + std::pow(stress(2), 2.0));
                return (mean - radius) < 0.0;
            }
            const mfem::real_t lambda1 = stress(0) - volumetric_stress(0) - std::abs(stress(5)) - std::abs(stress(4)); // s11 - |s12| - |s13|
            const mfem::real_t lambda2 = stress(1) - volumetric_stress(1) - std::abs(stress(5)) - std::abs(stress(3)); // s22 - |s12| - |s23|
            const mfem::real_t lambda3 = stress(2) - volumetric_stress(2) - std::abs(stress(4)) - std::abs(stress(3)); // s33 - |s13| - |s23|
            return std::min({lambda1, lambda2, lambda3}) < 0.0;
        }

    } // namespace detail
} // namespace mfemplus

#endif
