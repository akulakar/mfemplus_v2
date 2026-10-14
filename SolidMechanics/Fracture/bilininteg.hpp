//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_BILININTEG
#define MFEMPLUS_FRACTURE_BILININTEG

#include "mfem.hpp"

using namespace std;
using namespace mfem;
namespace mfemplus
{
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Bilinear form integrators for regularized variational fracture.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------------
    // Integrators for elasticity equation.
    //------------------------------------------------------------------------------------------------------------

    /** Integrator for the damaged elasticity equation with an isotropic material (Young's modulus E, Poisson's ratio nu).
       The stiffness is multiplied by the degradation function (1 - d)^{2} + k_{\epsilon}.
       $$
         a(u,v) =  ((1 - d)^{2} + k_{\epsilon})*(\mathrm{C}_{ijkl} u_{k,l} v_{i,j})
       $$
       The damage d is read from the grid function `damage_gf` and interpolated at each integration point.
       E and nu are evaluated at the first integration point (constant per element). In 2D, `plane_approximation`
       selects plane strain (0) or plane stress (1).
       This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space.
    **/
    class DamagedIsotropicElasticityIntegrator : public BilinearFormIntegrator
    {

    protected:
        Coefficient *young_mod, *poisson_ratio;
        real_t k_epsilon;
        GridFunction *damage_gf;
        FiniteElementSpace *damage_fes;
        Array<int> eldofs; // scalar for damage
        Vector eldofdamage;
        DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        int planeApprox;
        Vector shape;
        DenseMatrix dshape, gshape;
        // Optionally also include a option for the degradation function.

    public:
        DamagedIsotropicElasticityIntegrator(Coefficient &e, Coefficient &nu, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {}

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    /** Integrator for the damaged elasticity equation with a general (anisotropic, linear elastic) stiffness C, given in
       Voigt form as a MatrixCoefficient and evaluated at each integration point. Same form as
       DamagedIsotropicElasticityIntegrator:
       $$
         a(u,v) =  ((1 - d)^{2} + k_{\epsilon})*(\mathrm{C}_{ijkl} u_{k,l} v_{i,j})
       $$
       This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space.
    **/
    class DamagedAnisotropicElasticityIntegrator : public BilinearFormIntegrator
    {

    protected:
        MatrixCoefficient *stiffness;
        real_t k_epsilon;
        GridFunction *damage_gf;
        FiniteElementSpace *damage_fes;
        Array<int> eldofs; // scalar for damage
        Vector eldofdamage;
        DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        DamagedAnisotropicElasticityIntegrator(MatrixCoefficient &CMat, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace) : stiffness(&CMat), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace) {}

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    // Integrators for damage equation.
    //------------------------------------------------------------------------------------------------------------

    /** Integrator for the damage (phase field) equation: isotropic strain energy term.
       $$
         a(d,\phi) =  (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) d \phi = 2W(u) d \phi
       $$
       The displacement u is read from `disp_gf` (frozen at assembly time). E and nu are evaluated at the first
       integration point (constant per element).
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class IsotropicStrainEnergyDamageIntegrator : public BilinearFormIntegrator
    {

    protected:
        Coefficient *young_mod, *poisson_ratio;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B; // Stiffness in Voigt form

        int planeApprox;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyDamageIntegrator(Coefficient &e, Coefficient &nu, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    /** Integrator for the damage (phase field) equation: isotropic strain energy term with a volumetric pressure p.
       The pressure shifts the stress, so twice the energy density is
       $$
         2W = \mathrm{C}_{ijkl} u_{k,l} u_{i,j} - 2 p \, \mathrm{tr}(\varepsilon) + (\text{quadratic term in } p)
       $$
       $$
         a(d,\phi) =  2W(u, p) d \phi
       $$
       The pressure coefficient is evaluated at each integration point. E and nu are evaluated at the first integration
       point. Only the 3D and the 2D plane stress cases carry the pressure terms.
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class IsotropicStrainEnergyPressureDamageIntegrator : public BilinearFormIntegrator
    {

    protected:
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_pressure;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B; // Stiffness in Voigt form
        Vector body_pressure;

        int planeApprox;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyPressureDamageIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_pressure, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), volumetric_pressure(&vol_pressure), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    /** Integrator for the damage (phase field) equation: strain energy term with a general (anisotropic) stiffness C,
       given in Voigt form as a MatrixCoefficient and evaluated at each integration point.
       $$
         a(d,\phi) =  (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) d \phi
       $$
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class AnisotropicStrainEnergyDamageIntegrator : public BilinearFormIntegrator
    {
    protected:
        MatrixCoefficient *stiffness;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B; // Stiffness in Voigt form

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        AnisotropicStrainEnergyDamageIntegrator(MatrixCoefficient &CMat, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp_gridfunc), disp_fes(disp_fespace) {}

        void AssembleElementMatrix(const FiniteElement &el,
                                   ElementTransformation &Tr,
                                   DenseMatrix &elmat) override;
    };

    /** Integrator for the damage (phase field) equation: isotropic strain energy term with a volumetric eigenstrain
       \alpha I, where \alpha is the `vol_strain_coeff` coefficient (evaluated at each integration point).
       $$
         a(d, \phi) =  2W(u, \alpha) d \phi, \quad 2W = (\varepsilon - \alpha I) : \mathrm{C} : (\varepsilon - \alpha I)
       $$
       A Gershgorin-type check on the effective stress (stress - C \alpha I) switches the energy off (scales it by 1e-20)
       where the point is not guaranteed to be in tension, see `detail::EffectiveStressBoundIsNegative`.
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator : public BilinearFormIntegrator
    {
    protected:
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_strain_coeff;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B;                   // Stiffness in Voigt form
        Vector identity, volumetric_stress; // identity and volumetric stress (pressure)

        int planeApprox;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_strain_coeff, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), volumetric_strain_coeff(&vol_strain_coeff), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}

#endif
