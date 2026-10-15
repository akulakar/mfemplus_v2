//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_BILININTEG
#define MFEMPLUS_FRACTURE_BILININTEG

#include "mfem.hpp"
#include "mfemplus/SolidMechanics/Fracture/energy.hpp"

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
    class DamagedIsotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::real_t k_epsilon;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage;
        int planeApprox;
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;
        // Optionally also include a option for the degradation function.

    public:
        DamagedIsotropicElasticityIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
    };

    /** Integrator for the damaged elasticity equation with a general (anisotropic, linear elastic) stiffness C, given in
       Voigt form as a MatrixCoefficient and evaluated at each integration point. Same form as
       DamagedIsotropicElasticityIntegrator:
       $$
         a(u,v) =  ((1 - d)^{2} + k_{\epsilon})*(\mathrm{C}_{ijkl} u_{k,l} v_{i,j})
       $$
       This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space.
    **/
    class DamagedAnisotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::MatrixCoefficient *stiffness;
        mfem::real_t k_epsilon;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage;
        mfem::DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        DamagedAnisotropicElasticityIntegrator(mfem::MatrixCoefficient &CMat, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace) : stiffness(&CMat), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
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
    class IsotropicStrainEnergyDamageIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        IsotropicStrainEnergyDamageIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0) : energy_source(std::make_unique<IsotropicStrainEnergyEvaluator>(e, nu, disp_gridfunc, disp_fespace, plane_approximation)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyDamageIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
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
    class IsotropicStrainEnergyPressureDamageIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        IsotropicStrainEnergyPressureDamageIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_pressure, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0) : energy_source(std::make_unique<IsotropicStrainEnergyPressureEvaluator>(e, nu, vol_pressure, disp_gridfunc, disp_fespace, plane_approximation)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyPressureEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyPressureDamageIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
    };

    /** Integrator for the damage (phase field) equation: strain energy term with a general (anisotropic) stiffness C,
       given in Voigt form as a MatrixCoefficient and evaluated at each integration point.
       $$
         a(d,\phi) =  (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) d \phi
       $$
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class AnisotropicStrainEnergyDamageIntegrator : public mfem::BilinearFormIntegrator
    {
    protected:
        mfem::MatrixCoefficient *stiffness;
        mfem::GridFunction *disp_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::Array<int> eldofs;
        mfem::Vector eldofdisp;
        mfem::Vector CBu, Bu;
        mfem::DenseMatrix C, B; // Stiffness in Voigt form

        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        AnisotropicStrainEnergyDamageIntegrator(mfem::MatrixCoefficient &CMat, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp_gridfunc), disp_fes(disp_fespace) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };

    /** Integrator for the damage (phase field) equation: isotropic strain energy term with a volumetric eigenstrain
       \alpha I, where \alpha is the `vol_strain_coeff` coefficient (evaluated at each integration point).
       $$
         a(d, \phi) =  2W(u, \alpha) d \phi, \quad 2W = (\varepsilon - \alpha I) : \mathrm{C} : (\varepsilon - \alpha I)
       $$
       A Gershgorin-type check on the effective stress (stress - C \alpha I) switches the energy off (scales it by 1e-20; optional, enabled by default)
       where the point is not guaranteed to be in tension, see `detail::EffectiveStressBoundIsNegative`.
       This is a scalar integrator, i.e. defined for scalar FE spaces.
    **/
    class IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        /// `use_gershgorin_check` (default true) enables the Gershgorin-type tension check, see the class description.
        IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_strain_coeff, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0, bool use_gershgorin_check = true) : energy_source(std::make_unique<IsotropicStrainEnergyVolumetricEigenStrainEvaluator>(e, nu, vol_strain_coeff, disp_gridfunc, disp_fespace, plane_approximation, use_gershgorin_check)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyVolumetricEigenStrainEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}

#endif
