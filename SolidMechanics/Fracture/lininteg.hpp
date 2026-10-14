//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form intergrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

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

    // Integrators for variational fracture.
    /** New DomainLFIntegrator for damage force term.
     * For the second (phase field equation)
     **/
    class FractureDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Vector shape;
        Vector eldofdisp, eldofdamage;
        Array<int> eldofs;
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_pressure;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        Vector CBu, Bu;
        Vector body_pressure;
        int planeApprox;
        DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        FractureDamageLFIntegrator(Coefficient &e, Coefficient &nu, GridFunction &disp, FiniteElementSpace *disp_fespace, int plane_approximation = 0, Coefficient *vol_pressure = nullptr)
            : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_pressure(vol_pressure) {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    class IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        Vector shape;
        Vector eldofdisp, eldofdamage;
        Array<int> eldofs;
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_strain_coeff;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        Vector CBu, Bu;
        Vector identity, volumetric_stress; // identity and volumetric stress (pressure)
        int planeApprox;
        DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_strain_coeff, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_strain_coeff(&vol_strain_coeff) {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect.
        **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    class AnisotropicElasticityDamageLFIntegrator : public LinearFormIntegrator
    {
    protected:
        MatrixCoefficient *stiffness;
        Vector shape;
        Vector eldofdisp, eldofdamage;
        Array<int> eldofs;
        GridFunction *disp_gf;
        GridFunction *damage_gf;
        FiniteElementSpace *disp_fes;
        DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        Vector CBu, Bu;
        DenseMatrix dshape, gshape;

    public:
        // Constructs a domain integrator with a given Coefficient
        AnisotropicElasticityDamageLFIntegrator(MatrixCoefficient &CMat, GridFunction &disp, GridFunction &damage, FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp), damage_gf(&damage), disp_fes(disp_fespace) {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;
    };

    /** New DomainLFIntegrator for damaged elasticity force by pressure body force.
     * For the first (damgaged elasticity) equation
     **/
    class PressureBodyForceDamageLFIntegrator : public LinearFormIntegrator
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

        PressureBodyForceDamageLFIntegrator(Coefficient &pressure_coeff, Coefficient &e, Coefficient &nu, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace, int plane_approximation)
            : pressure(&pressure_coeff), young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {};

        void AssembleDevice(const FiniteElementSpace &fes, const Array<int> &markers, Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        virtual void AssembleRHSElementVect(const FiniteElement &el, FaceElementTransformations &Tr, Vector &elvect) override {};

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    /** New DomainLFIntegrator for damaged elasticity force by volumetric eigenstrain which acts as a body force.
     * For the first (damgaged elasticity) equation
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
            : volumetric_strain_coeff(&volumetric_strain_coeff_in), young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override;

        using LinearFormIntegrator::AssembleRHSElementVect;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}