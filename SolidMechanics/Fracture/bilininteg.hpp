//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form intergrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

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

    /** Integrator for variational fracture strain energy term using general stiffness tensor, specifying Young's modulus E and Poisson's ratio nu, and multiplied with the degradation damage function (1 - d^{2}) + k_{\epsilon}.
       $$
         a(u,v) =  ((1 - d)^{2} + k_{\epsilon})*(\mathrm{C}_{ijkl} u_{k,l} v_{i,j})
       $$
       This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space. **/
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
        DamagedIsotropicElasticityIntegrator(Coefficient &e, Coefficient &nu, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {};

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    // For anisotropic (linear elastic) material models.
    class AnisotropicElasticityDamageIntegrator : public BilinearFormIntegrator
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
        AnisotropicElasticityDamageIntegrator(MatrixCoefficient &CMat, real_t k_eps, GridFunction &damage_gridfunc, FiniteElementSpace *damage_fespace) : stiffness(&CMat), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace) {};

        void AssembleElementMatrix(const FiniteElement &el,
                                   ElementTransformation &Tr,
                                   DenseMatrix &elmat) override;
    };

    /** Integrator for variational fracture, elastic energy and damage term, using Young's modulus E and Poisson's ratio nu
       $$
         a(u,v) =  (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) d \phi
       $$
       This is a scalar integrator, i.e. defined for scalar FE spaces. **/
    class IsotropicStrainEnergyDamageIntegrator : public BilinearFormIntegrator
    {

    protected:
        Coefficient *young_mod, *poisson_ratio;
        Coefficient *volumetric_pressure;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        Vector body_pressure;

        int planeApprox, GershgorinCheck;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyDamageIntegrator(Coefficient &e, Coefficient &nu, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0, Coefficient *vol_press = nullptr) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_pressure(vol_press) {};

        void AssembleElementMatrix(const FiniteElement &el,
                                   ElementTransformation &Tr,
                                   DenseMatrix &elmat) override;
    };

    /** Integrator for variational fracture, elastic energy and damage term, using Young's modulus E and Poisson's ratio nu with a volumetric eigenstrain
       $$
         a(d, \phi) =  W(u, \nabla u) d \phi
       $$
       This is a scalar integrator, i.e. defined for scalar FE spaces. **/
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
        DenseMatrix C, B, CB, elmat_input;  // Stiffness in Voigt form
        Vector identity, volumetric_stress; // identity and volumetric stress (pressure)

        int planeApprox, GershgorinCheck;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(Coefficient &e, Coefficient &nu, Coefficient &vol_strain_coeff, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_strain_coeff(&vol_strain_coeff) {};

        void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr, DenseMatrix &elmat) override;
    };

    class AnisotropicStrainEnergyDamageIntegrator : public BilinearFormIntegrator
    {

    protected:
        MatrixCoefficient *stiffness;
        GridFunction *disp_gf;
        FiniteElementSpace *disp_fes;
        Array<int> eldofs;
        Vector eldofdamage;
        Vector eldofdisp;
        Vector CBu, Bu;
        DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form

        int planeApprox, GershgorinCheck;

        Vector shape;
        DenseMatrix dshape, gshape;

    public:
        AnisotropicStrainEnergyDamageIntegrator(MatrixCoefficient &CMat, GridFunction &disp_gridfunc, FiniteElementSpace *disp_fespace) : stiffness(&CMat), disp_gf(&disp_gridfunc), disp_fes(disp_fespace) {};

        void AssembleElementMatrix(const FiniteElement &el,
                                   ElementTransformation &Tr,
                                   DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}