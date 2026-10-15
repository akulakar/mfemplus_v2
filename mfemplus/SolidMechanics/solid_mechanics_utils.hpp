//------------------------------------------------------------------------------------------------------------
// Internal helpers shared by the integrators for solid mechanics.
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_SOLID_MECHANICS_UTILS
#define MFEMPLUS_SOLID_MECHANICS_UTILS

#include "mfem.hpp"
#include <algorithm>
#include <cmath>


namespace mfemplus
{
    namespace detail
    {
        /// Number of independent strain components in Voigt form.
        inline int NumStrainComponents(int dim) { return (dim == 2) ? 3 : 6; }

        /// Gathers the (vector) dofs of element `elnum` and the values of `gf` on them.
        /// For a scalar space this is the usual element dofs.
        inline void GatherElementValues(const mfem::GridFunction &gf, const mfem::FiniteElementSpace &fes,
                                        int elnum, mfem::Array<int> &dofs, mfem::Vector &values)
        {
            fes.GetElementVDofs(elnum, dofs);
            gf.GetSubVector(dofs, values);
        }

        /// Returns `ir` if set, otherwise a rule of order 2 p + (order of the transformation gradient).
        inline const mfem::IntegrationRule &IntRuleOrDefault(const mfem::IntegrationRule *ir, const mfem::FiniteElement &el,
                                                       mfem::ElementTransformation &Tr)
        {
            if (ir != nullptr)
            {
                return *ir;
            }
            return mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderGrad(&el));
        }

        /// Fills the Voigt strain-displacement matrix B (str_comp x dof*dim) from the physical shape function
        /// gradients gshape (dof x dim). B must already be sized and zero-initialized; only its non-zero
        /// entries are written. Displacement dofs are ordered component by component.
        inline void FillStrainDisplacement(const mfem::DenseMatrix &gshape, int dim, mfem::DenseMatrix &B)
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

        /// Lame parameters of the isotropic stiffness C = lambda I(x)I + 2 mu I_sym. In 2D, planeApprox = 0 is plane
        /// strain (same lambda as 3D) and planeApprox = 1 is plane stress (lambda = E nu / (1 - nu^2)).
        inline void IsotropicLame(mfem::real_t E, mfem::real_t NU, int dim, int planeApprox, mfem::real_t &lambda, mfem::real_t &mu)
        {
            mu = E / (2.0 * (1.0 + NU));
            lambda = (dim == 2 && planeApprox == 1) ? E * NU / (1.0 - NU * NU)
                                                    : E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
        }

        /// Adds scale * (isotropic element stiffness) to elmat, without forming B or C. For nodes a and b with physical
        /// shape function gradients g_a and g_b the (dim x dim) block is
        ///   K_ab = lambda g_a g_b^T + mu ((g_a . g_b) I + g_b g_a^T).
        /// elmat is (dof*dim) square, displacement dofs ordered component by component (index a + i*dof). Only the
        /// blocks with a <= b are computed; the others are mirrored.
        inline void AddIsotropicStiffness(const mfem::DenseMatrix &gshape, int dim, mfem::real_t scale, mfem::real_t lambda,
                                          mfem::real_t mu, mfem::DenseMatrix &elmat)
        {
            const int dof = gshape.Height();
            const mfem::real_t lam = scale * lambda;
            const mfem::real_t m = scale * mu;
            for (int a = 0; a < dof; a++)
            {
                for (int b = a; b < dof; b++)
                {
                    mfem::real_t dot = 0.0;
                    for (int k = 0; k < dim; k++)
                    {
                        dot += gshape(a, k) * gshape(b, k);
                    }
                    for (int i = 0; i < dim; i++)
                    {
                        for (int j = 0; j < dim; j++)
                        {
                            mfem::real_t v = lam * gshape(a, i) * gshape(b, j) + m * gshape(b, i) * gshape(a, j);
                            if (i == j)
                            {
                                v += m * dot;
                            }
                            elmat(a + i * dof, b + j * dof) += v;
                            if (a != b)
                            {
                                elmat(b + j * dof, a + i * dof) += v;
                            }
                        }
                    }
                }
            }
        }

        /// Voigt strain (engineering shear) at a point from the element displacement dofs `u` (ordered component by
        /// component) and the physical shape function gradients gshape (dof x dim). Equivalent to B u, without B.
        /// Order: (xx, yy, xy) in 2D and (xx, yy, zz, yz, xz, xy) in 3D.
        inline void VoigtStrain(const mfem::DenseMatrix &gshape, int dim, const mfem::Vector &u, mfem::Vector &strain)
        {
            const int dof = gshape.Height();
            mfem::real_t h[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}; // h[i][j] = d u_i / d x_j
            for (int i = 0; i < dim; i++)
            {
                for (int a = 0; a < dof; a++)
                {
                    const mfem::real_t ua = u(a + i * dof);
                    for (int j = 0; j < dim; j++)
                    {
                        h[i][j] += ua * gshape(a, j);
                    }
                }
            }
            strain(0) = h[0][0];
            strain(1) = h[1][1];
            if (dim == 2)
            {
                strain(2) = h[0][1] + h[1][0];
            }
            else
            {
                strain(2) = h[2][2];
                strain(3) = h[1][2] + h[2][1];
                strain(4) = h[0][2] + h[2][0];
                strain(5) = h[0][1] + h[1][0];
            }
        }

        /// Isotropic stress C strain in Voigt form: sigma_ii = lambda tr(eps) + 2 mu eps_ii, sigma_shear = mu gamma.
        inline void IsotropicStress(const mfem::Vector &strain, int dim, mfem::real_t lambda, mfem::real_t mu, mfem::Vector &stress)
        {
            mfem::real_t trace = 0.0;
            for (int i = 0; i < dim; i++)
            {
                trace += strain(i);
            }
            for (int i = 0; i < dim; i++)
            {
                stress(i) = lambda * trace + 2.0 * mu * strain(i);
            }
            for (int i = dim; i < NumStrainComponents(dim); i++)
            {
                stress(i) = mu * strain(i);
            }
        }

        /// Twice the elastic energy density, u^T B^T C B u, for the isotropic stiffness, without B or C.
        /// Strain is returned in `strain` and stress in `stress`.
        inline mfem::real_t TwiceIsotropicStrainEnergy(const mfem::DenseMatrix &gshape, int dim, const mfem::Vector &u, mfem::real_t lambda,
                                                 mfem::real_t mu, mfem::Vector &strain, mfem::Vector &stress)
        {
            VoigtStrain(gshape, dim, u, strain);
            IsotropicStress(strain, dim, lambda, mu, stress);
            return mfem::InnerProduct(stress, strain);
        }

        /// Volumetric stress C (alpha I) in Voigt form for the isotropic stiffness: (dim lambda + 2 mu) alpha on the
        /// normal components, zero shear.
        inline void IsotropicVolumetricStress(int dim, mfem::real_t lambda, mfem::real_t mu, mfem::real_t alpha, mfem::Vector &volumetric_stress)
        {
            volumetric_stress = 0.0;
            for (int i = 0; i < dim; i++)
            {
                volumetric_stress(i) = (dim * lambda + 2.0 * mu) * alpha;
            }
        }

        /// out = B^T s for a stress s that only has normal components (shear entries are ignored):
        /// out(a + i*dof) = s_i g_a,i. Used for the pressure and eigenstrain body forces.
        inline void StrainDisplacementTransposeNormal(const mfem::DenseMatrix &gshape, int dim, const mfem::Vector &s, mfem::Vector &out)
        {
            const int dof = gshape.Height();
            for (int i = 0; i < dim; i++)
            {
                for (int a = 0; a < dof; a++)
                {
                    out(a + i * dof) = s(i) * gshape(a, i);
                }
            }
        }

        /// Voigt representation of the identity tensor: (1 1 1 0 0 0) in 3D, (1 1 0) in 2D.
        inline void FillVoigtIdentity(int dim, mfem::Vector &identity)
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
        inline mfem::real_t PlaneStressVolumetricFactor(mfem::real_t NU) { return 2.0 - 1.0 / (1.0 - NU); }

        /// Twice the elastic energy density, u^T B^T C B u, with strain = B u and stress = C strain.
        /// Strain is returned in Bu and stress in CBu.
        inline mfem::real_t TwiceStrainEnergy(const mfem::DenseMatrix &C, const mfem::DenseMatrix &B,
                                        const mfem::Vector &eldofdisp, mfem::Vector &Bu, mfem::Vector &CBu)
        {
            B.Mult(eldofdisp, Bu);
            C.Mult(Bu, CBu);
            return mfem::InnerProduct(CBu, Bu);
        }

        /// Sets the body pressure (in Voigt form) and twice-the-quadratic pressure energy term used by the
        /// isotropic strain energy with a volumetric pressure. `body_pressure` has all components set (or left at
        /// zero); `quadratic_energy` is only overwritten when a pressure applies.
        inline void SetPressureTerms(int dim, int planeApprox, mfem::real_t pressure, mfem::real_t E, mfem::real_t NU,
                                     mfem::Vector &body_pressure, mfem::real_t &quadratic_energy)
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
