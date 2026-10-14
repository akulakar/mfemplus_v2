//------------------------------------------------------------------------------------------------------------
// Internal helpers shared by the fracture bilinear and linear form integrators (Voigt notation kernels).
// Not part of the public interface.
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_UTILS
#define MFEMPLUS_FRACTURE_UTILS

#include "mfem.hpp"
#include <algorithm>
#include <cmath>

using namespace mfem;

namespace mfemplus
{
    namespace detail
    {
        /// Number of independent strain components in Voigt form.
        inline int NumStrainComponents(int dim) { return (dim == 2) ? 3 : 6; }

        /// Gathers the (vector) dofs of element `elnum` and the values of `gf` on them.
        /// For a scalar space this is the usual element dofs.
        inline void GatherElementValues(const GridFunction &gf, const FiniteElementSpace &fes,
                                        int elnum, Array<int> &dofs, Vector &values)
        {
            fes.GetElementVDofs(elnum, dofs);
            gf.GetSubVector(dofs, values);
        }

        /// Returns `ir` if set, otherwise a rule of order 2 p + (order of the transformation gradient).
        inline const IntegrationRule &IntRuleOrDefault(const IntegrationRule *ir, const FiniteElement &el,
                                                       ElementTransformation &Tr)
        {
            if (ir != nullptr)
            {
                return *ir;
            }
            return IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderGrad(&el));
        }

        /// Degradation function (1 - d)^2 + k_epsilon.
        inline real_t Degradation(real_t damage, real_t k_epsilon)
        {
            return (1.0 - damage) * (1.0 - damage) + k_epsilon;
        }

        /// Fills the Voigt strain-displacement matrix B (str_comp x dof*dim) from the physical shape function
        /// gradients gshape (dof x dim). B must already be sized and zero-initialized; only its non-zero
        /// entries are written. Displacement dofs are ordered component by component.
        inline void FillStrainDisplacement(const DenseMatrix &gshape, int dim, DenseMatrix &B)
        {
            const int dof = gshape.Height();
            if (dim == 2)
            {
                for (int a = 0; a < dof; a++)
                {
                    B(0, a) = gshape(a, 0);
                    B(1, a + dof) = gshape(a, 1);
                    B(2, a) = gshape(a, 1);
                    B(2, a + dof) = gshape(a, 0);
                }
            }
            else
            {
                for (int a = 0; a < dof; a++)
                {
                    B(0, a) = gshape(a, 0);
                    B(1, a + dof) = gshape(a, 1);
                    B(2, a + 2 * dof) = gshape(a, 2);
                    B(3, a + dof) = gshape(a, 2);
                    B(3, a + 2 * dof) = gshape(a, 1);
                    B(4, a) = gshape(a, 2);
                    B(4, a + 2 * dof) = gshape(a, 0);
                    B(5, a) = gshape(a, 1);
                    B(5, a + dof) = gshape(a, 0);
                }
            }
        }

        /// Fills the isotropic stiffness C in Voigt form. C must be sized and zero-initialized.
        /// In 2D, planeApprox = 0 is plane strain and planeApprox = 1 is plane stress.
        inline void FillIsotropicStiffness(real_t E, real_t NU, int dim, int planeApprox, DenseMatrix &C)
        {
            const real_t G = E / (2.0 * (1.0 + NU));
            if (dim == 2)
            {
                if (planeApprox == 0) // plane strain
                {
                    C(0, 0) = C(1, 1) = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
                    C(0, 1) = C(1, 0) = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
                    C(2, 2) = G;
                }
                else if (planeApprox == 1) // plane stress
                {
                    C(0, 0) = C(1, 1) = E / (1.0 - NU * NU);
                    C(0, 1) = C(1, 0) = E * NU / (1.0 - NU * NU);
                    C(2, 2) = G;
                }
            }
            else
            {
                const real_t lambda_plus_2G = E * (1.0 - NU) / ((1.0 - 2.0 * NU) * (1.0 + NU));
                const real_t lambda = E * NU / ((1.0 - 2.0 * NU) * (1.0 + NU));
                C(0, 0) = C(1, 1) = C(2, 2) = lambda_plus_2G;
                C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = lambda;
                C(3, 3) = C(4, 4) = C(5, 5) = G;
            }
        }

        /// Voigt representation of the identity tensor: (1 1 1 0 0 0) in 3D, (1 1 0) in 2D.
        inline void FillVoigtIdentity(int dim, Vector &identity)
        {
            identity = 0.0;
            identity(0) = 1.0;
            identity(1) = 1.0;
            if (dim == 3)
            {
                identity(2) = 1.0;
            }
        }

        /// Factor applied to a 3D volumetric load in the plane stress approximation.
        inline real_t PlaneStressVolumetricFactor(real_t NU) { return 2.0 - 1.0 / (1.0 - NU); }

        /// Volumetric stress C (alpha I) in Voigt form, for the isotropic eigenstrain alpha I.
        inline void VolumetricStress(const DenseMatrix &C, const Vector &identity, real_t alpha,
                                     Vector &volumetric_stress)
        {
            C.Mult(identity, volumetric_stress);
            volumetric_stress *= alpha;
        }

        /// Twice the elastic energy density, u^T B^T C B u, with strain = B u and stress = C strain.
        /// Strain is returned in Bu and stress in CBu.
        inline real_t TwiceStrainEnergy(const DenseMatrix &C, const DenseMatrix &B,
                                        const Vector &eldofdisp, Vector &Bu, Vector &CBu)
        {
            B.Mult(eldofdisp, Bu);
            C.Mult(Bu, CBu);
            return InnerProduct(CBu, Bu);
        }

        /// Gershgorin-type lower bound on the eigenvalues of the effective stress (stress - volumetric_stress), in
        /// Voigt form. Returns true if the bound is negative, i.e. the point is not guaranteed to be in tension
        /// and the damage driving energy is switched off by the eigenstrain integrators.
        /// In 2D the exact smallest principal value is used.
        inline bool EffectiveStressBoundIsNegative(int dim, const Vector &stress, const Vector &volumetric_stress)
        {
            if (dim == 2)
            {
                const real_t mean = (stress(0) - volumetric_stress(0) + stress(1) - volumetric_stress(1)) / 2.0;
                const real_t radius = std::sqrt(std::pow((stress(0) - stress(1)) / 2.0, 2.0) + std::pow(stress(2), 2.0));
                return (mean - radius) < 0.0;
            }
            const real_t lambda1 = stress(0) - volumetric_stress(0) - std::abs(stress(5)) - std::abs(stress(4)); // s11 - |s12| - |s13|
            const real_t lambda2 = stress(1) - volumetric_stress(1) - std::abs(stress(5)) - std::abs(stress(3)); // s22 - |s12| - |s23|
            const real_t lambda3 = stress(2) - volumetric_stress(2) - std::abs(stress(4)) - std::abs(stress(3)); // s33 - |s13| - |s23|
            return std::min({lambda1, lambda2, lambda3}) < 0.0;
        }

        /// Sets the body pressure (in Voigt form) and twice-the-quadratic pressure energy term used by the
        /// isotropic strain energy with a volumetric pressure. `body_pressure` has all components set (or left at
        /// zero); `quadratic_energy` is only overwritten when a pressure applies.
        inline void SetPressureTerms(int dim, int planeApprox, real_t pressure, real_t E, real_t NU,
                                     Vector &body_pressure, real_t &quadratic_energy)
        {
            if (dim == 2)
            {
                if (planeApprox == 1) // plane stress. The plane strain case is not handled.
                {
                    body_pressure(0) = body_pressure(1) = pressure * PlaneStressVolumetricFactor(NU);
                    quadratic_energy = (pressure * pressure * (1.0 - 2.0 * NU) * (1.0 - 2.0 * NU)) / (E * (1.0 - NU));
                }
            }
            else
            {
                body_pressure(0) = body_pressure(1) = body_pressure(2) = pressure;
                quadratic_energy = (3.0 * pressure * pressure) / ((2.0 * E) / (1.0 - 2.0 * NU));
            }
        }
    } // namespace detail
} // namespace mfemplus

#endif
