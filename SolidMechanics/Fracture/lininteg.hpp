//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_LININTEG
#define MFEMPLUS_FRACTURE_LININTEG

#include "mfem.hpp"

using namespace std;
using namespace mfem;
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
    class DamagedPressureBodyForceLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Coefficient *pressure;
        Coefficient *young_mod, *poisson_ratio;
        DenseMatrix B;        // stiffness, strain-displacement
        Vector body_pressure; // volumetric pressure
        DenseMatrix dshape, gshape;
        Vector shape;
        GridFunction *damage_gf;
        FiniteElementSpace *damage_fes;
        Array<int> eldofs; // scalar for damage
        Vector eldofdamage, elvec_input;
        real_t k_epsilon;
        int planeApprox;

    public:
        // Constructs a domain integrator with a given Coefficient

        DamagedPressureBodyForceLFIntegrator(Coefficient &pressure_coeff, Coefficient &e, Coefficient &nu, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace, int plane_approximation)
            : pressure(&pressure_coeff), young_mod(&e), poisson_ratio(&nu), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), k_epsilon(k_eps), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damaged elasticity equation: body force from a volumetric eigenstrain \alpha I
     * (coefficient `volumetric_strain_coeff_in`), i.e. the stress C (\alpha I) acting as an isotropic stress, scaled by
     * the degradation function (1 - d)^{2} + k_{\epsilon}.
     * $$
     *   l(v) = \int ((1 - d)^{2} + k_{\epsilon}) \, (\mathrm{C} \alpha I) : \varepsilon(v)
     * $$
     * The damage d is read from `damage_gf`. This is a 'Vector' integrator (vector valued displacement space).
     **/
    class DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Coefficient *volumetric_strain_coeff;
        Coefficient *young_mod, *poisson_ratio;
        DenseMatrix C, B;                   // stiffness, strain-displacement
        Vector identity, volumetric_stress; // identity and volumetric stress (pressure)
        DenseMatrix dshape, gshape;
        Vector shape;
        GridFunction *damage_gf;
        FiniteElementSpace *damage_fes;
        Array<int> eldofs; // scalar for damage
        Vector eldofdamage, elvec_input;
        real_t k_epsilon;
        int planeApprox;

    public:
        // Constructs a domain integrator with a given Coefficient

        DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator(Coefficient &e, Coefficient &nu, Coefficient &volumetric_strain_coeff_in, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace, int plane_approximation = 0)
            : volumetric_strain_coeff(&volumetric_strain_coeff_in), young_mod(&e), poisson_ratio(&nu), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), k_epsilon(k_eps), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
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
    class IsotropicStrainEnergyDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Vector shape;
        Vector eldofdisp;
        Array<int> eldofs;
        Coefficient *young_mod, *poisson_ratio;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B; // stiffness and strain-displacement in Voigt form
        Vector CBu, Bu;
        int planeApprox;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyDamageLFIntegrator(Coefficient &e, Coefficient &nu, GridFunction &disp, FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: isotropic strain energy driving force with a
     * volumetric pressure p, i.e. twice the energy of the stress state (sigma - p I) is integrated against phi.
     * $$
     *   l(\phi) = \int 2W(u, p) \phi, \quad
     *   2W = \mathrm{C}_{ijkl} u_{k,l} u_{i,j} - 2 p \, \mathrm{tr}(\varepsilon) + (\text{quadratic term in } p)
     * $$
     * The pressure is evaluated at each integration point. This is a scalar integrator, i.e. defined for scalar FE spaces.
     **/
    class IsotropicStrainEnergyPressureDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Vector shape;
        Vector eldofdisp;
        Array<int> eldofs;
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_pressure;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B; // stiffness and strain-displacement in Voigt form
        Vector CBu, Bu;
        Vector body_pressure;
        int planeApprox;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyPressureDamageLFIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_pressure, GridFunction &disp, FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : young_mod(&e), poisson_ratio(&nu), volumetric_pressure(&vol_pressure), disp_gf(&disp), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: isotropic strain energy driving force with a
     * volumetric eigenstrain \alpha I (coefficient `vol_strain_coeff`, evaluated at each integration point).
     * $$
     *   l(\phi) = \int 2W(u, \alpha) \phi, \quad 2W = (\varepsilon - \alpha I) : \mathrm{C} : (\varepsilon - \alpha I)
     * $$
     * Uses the same energy and Gershgorin-type switch as IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator, so
     * the matrix and the right hand side of the damage equation are consistent. Scalar integrator.
     **/
    class IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Vector shape;
        Vector eldofdisp;
        Array<int> eldofs;
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_strain_coeff;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B; // stiffness and strain-displacement in Voigt form
        Vector CBu, Bu;
        Vector identity, volumetric_stress; // identity and volumetric stress (pressure)
        int planeApprox;
        DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_strain_coeff, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), volumetric_strain_coeff(&vol_strain_coeff), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect.
        **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** DomainLFIntegrator for the damage (phase field) equation: strain energy driving force with a general
     * (anisotropic) stiffness C, given in Voigt form as a MatrixCoefficient and evaluated at each integration point.
     * $$
     *   l(\phi) = \int (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) \phi
     * $$
     * The displacement u is read from `disp_gf`. Note `damage_gf` is stored but not used by this integrator.
     * This is a scalar integrator, i.e. defined for scalar FE spaces.
     **/
    class AnisotropicElasticityDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        MatrixCoefficient *stiffness;
        Vector shape;
        Vector eldofdisp;
        Array<int> eldofs;
        GridFunction *disp_gf;
        GridFunction *damage_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B; // stiffness and strain-displacement in Voigt form
        Vector CBu, Bu;
        DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        AnisotropicElasticityDamageLFIntegrator(MatrixCoefficient &CMat, GridFunction &disp, GridFunction &damage, FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp), damage_gf(&damage), disp_fes(disp_fespace) {}

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture linear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}

#endif
