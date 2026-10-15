//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_LININTEG
#define MFEMPLUS_FRACTURE_LININTEG

#include "mfem.hpp"
#include "mfemplus/SolidMechanics/Fracture/energy.hpp"

namespace mfemplus
{
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Linear form integrators for regularized variational fracture.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------------
    // Integrators for elasticity equation
    //------------------------------------------------------------------------------------------------------------

    /** DomainLFIntegrator for the damaged elasticity equation: body force from a volumetric pressure p acting as an
     * isotropic stress, scaled by the degradation function (1 - d)^{2} + k_{\epsilon}.
     * $$
     *   l(v) = \int ((1 - d)^{2} + k_{\epsilon}) \, p \, \mathrm{div}(v)
     * $$
     * The damage d is read from `damage_gf` and interpolated at each integration point; p is evaluated at each
     * integration point. This is a 'Vector' integrator (vector valued displacement space).
     * In 2D only plane stress carries the pressure load; plane strain assembles zero.
     **/
    class DamagedPressureBodyForceLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Coefficient *pressure;
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::Vector body_pressure; // volumetric pressure
        mfem::DenseMatrix dshape, gshape;
        mfem::Vector shape;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage, elvec_input;
        mfem::real_t k_epsilon;
        int planeApprox;

    public:
        // Constructs a domain integrator with a given Coefficient

        DamagedPressureBodyForceLFIntegrator(mfem::Coefficient &pressure_coeff, mfem::Coefficient &e, mfem::Coefficient &nu, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace, int plane_approximation)
            : pressure(&pressure_coeff), young_mod(&e), poisson_ratio(&nu), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), k_epsilon(k_eps), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damaged elasticity equation: body force from a volumetric eigenstrain \alpha I
     * (coefficient `volumetric_strain_coeff_in`), i.e. the stress C (\alpha I) acting as an isotropic stress, scaled by
     * the degradation function (1 - d)^{2} + k_{\epsilon}.
     * $$
     *   l(v) = \int ((1 - d)^{2} + k_{\epsilon}) \, (\mathrm{C} \alpha I) : \varepsilon(v)
     * $$
     * The damage d is read from `damage_gf`. This is a 'Vector' integrator (vector valued displacement space).
     **/
    class DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Coefficient *volumetric_strain_coeff;
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::Vector volumetric_stress; // volumetric stress (pressure)
        mfem::DenseMatrix dshape, gshape;
        mfem::Vector shape;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage, elvec_input;
        mfem::real_t k_epsilon;
        int planeApprox;

    public:
        // Constructs a domain integrator with a given Coefficient

        DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &volumetric_strain_coeff_in, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace, int plane_approximation = 0)
            : volumetric_strain_coeff(&volumetric_strain_coeff_in), young_mod(&e), poisson_ratio(&nu), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), k_epsilon(k_eps), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    //------------------------------------------------------------------------------------------------------------
    // Integrators for damage equation
    //------------------------------------------------------------------------------------------------------------

    /** DomainLFIntegrator for the damage (phase field) equation: isotropic strain energy driving force.
     * $$
     *   l(\phi) = \int 2W(u) \phi, \quad 2W = \mathrm{C}_{ijkl} u_{k,l} u_{i,j}
     * $$
     * The displacement u is read from `disp_gf`. E and nu are evaluated at the first integration point (constant per
     * element). This is a scalar integrator, i.e. defined for scalar FE spaces.
     **/
    class IsotropicStrainEnergyDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        IsotropicStrainEnergyDamageLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : energy_source(std::make_unique<IsotropicStrainEnergyEvaluator>(e, nu, disp, disp_fespace, plane_approximation)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyDamageLFIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: isotropic strain energy driving force with a
     * volumetric pressure p, i.e. twice the energy of the stress state (sigma - p I) is integrated against phi.
     * $$
     *   l(\phi) = \int 2W(u, p) \phi, \quad
     *   2W = \mathrm{C}_{ijkl} u_{k,l} u_{i,j} - 2 p \, \mathrm{tr}(\varepsilon) + (\text{quadratic term in } p)
     * $$
     * The pressure is evaluated at each integration point. This is a scalar integrator, i.e. defined for scalar FE spaces.
     **/
    class IsotropicStrainEnergyPressureDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        IsotropicStrainEnergyPressureDamageLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_pressure, mfem::GridFunction &disp, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : energy_source(std::make_unique<IsotropicStrainEnergyPressureEvaluator>(e, nu, vol_pressure, disp, disp_fespace, plane_approximation)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyPressureEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyPressureDamageLFIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: isotropic strain energy driving force with a
     * volumetric eigenstrain \alpha I (coefficient `vol_strain_coeff`, evaluated at each integration point).
     * $$
     *   l(\phi) = \int 2W(u, \alpha) \phi, \quad 2W = (\varepsilon - \alpha I) : \mathrm{C} : (\varepsilon - \alpha I)
     * $$
     * Uses the same energy and Gershgorin-type switch (optional) as IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator, so
     * the matrix and the right hand side of the damage equation are consistent. Scalar integrator.
     **/
    class IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        StrainEnergySource energy_source;
        mfem::Vector energy, shape;

    public:
        /// `use_gershgorin_check` (default true) enables the Gershgorin-type tension check, see the class description.
        IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_strain_coeff, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0, bool use_gershgorin_check = true)
            : energy_source(std::make_unique<IsotropicStrainEnergyVolumetricEigenStrainEvaluator>(e, nu, vol_strain_coeff, disp_gridfunc, disp_fespace, plane_approximation, use_gershgorin_check)) {}

        /// Uses the shared 2W `twice_energy` (filled by `IsotropicStrainEnergyVolumetricEigenStrainEvaluator::Evaluate`) instead of evaluating the energy itself.
        /// The rule of its QuadratureSpace is used for the integration. `twice_energy` must outlive the integrator.
        explicit IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(const mfem::QuadratureFunction &twice_energy) : energy_source(twice_energy) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: strain energy driving force with a general
     * (anisotropic) stiffness C, given in Voigt form as a MatrixCoefficient and evaluated at each integration point.
     * $$
     *   l(\phi) = \int (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) \phi
     * $$
     * The displacement u is read from `disp_gf`. Note `damage_gf` is stored but not used by this integrator.
     * This is a scalar integrator, i.e. defined for scalar FE spaces.
     **/
    class AnisotropicElasticityDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::MatrixCoefficient *stiffness;
        mfem::Vector shape;
        mfem::Vector eldofdisp;
        mfem::Array<int> eldofs;
        mfem::GridFunction *disp_gf;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::DenseMatrix C, B; // stiffness and strain-displacement in Voigt form
        mfem::Vector CBu, Bu;
        mfem::DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        AnisotropicElasticityDamageLFIntegrator(mfem::MatrixCoefficient &CMat, mfem::GridFunction &disp, mfem::GridFunction &damage, mfem::FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp), damage_gf(&damage), disp_fes(disp_fespace) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture linear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}

#endif
