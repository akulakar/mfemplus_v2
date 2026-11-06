//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Nonlinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_NONLININTEG
#define MFEMPLUS_FRACTURE_NONLININTEG

#include "mfem.hpp"
#include "mfemplus/SolidMechanics/Fracture/energy.hpp"

namespace mfemplus
{
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Nonlinear form integrators for regularized variational fracture.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    class RVFT_AT2_Integrator : public mfem::BlockNonlinearFormIntegrator
    {
    private:
        mfem::Coefficient *Ey_coeff, *nu_coeff, *gOverEps_coeff, *gTimesEps_coeff;
        double k_eps;
        mfem::ParFiniteElementSpace *dispfes, *dmgfes;
        mfem::ParGridFunction *disp, *dmg;

        // Scratch storage, sized per element in the methods.
        mfem::DenseMatrix dshape_u, gshape_u, dshape_d;
        mfem::Vector shape_d, grad_d, strain, stress;

    public:
        RVFT_AT2_Integrator(mfem::Coefficient &Ey_, mfem::Coefficient &nu_, mfem::Coefficient &gOverEps_, mfem::Coefficient &gTimesEps_, double k_eps_ = 1e-8)
            : k_eps(k_eps_), dispfes(nullptr), dmgfes(nullptr), disp(nullptr), dmg(nullptr)
        {
            ReadCoefficient(Ey_, nu_, gOverEps_, gTimesEps_);
        }

        void ReadCoefficient(mfem::Coefficient &Ey_, mfem::Coefficient &nu_, mfem::Coefficient &gOverEps_, mfem::Coefficient &gTimesEps_)
        {
            Ey_coeff = &Ey_;
            nu_coeff = &nu_;
            gOverEps_coeff = &gOverEps_;
            gTimesEps_coeff = &gTimesEps_;
        }

        void FESpacesGridFunctions(mfem::ParFiniteElementSpace &dispfes_, mfem::ParFiniteElementSpace &dmgfes_, mfem::ParGridFunction &disp_gf, mfem::ParGridFunction &dmg_gf)
        {
            dispfes = &dispfes_;
            dmgfes = &dmgfes_;
            disp = &disp_gf;
            dmg = &dmg_gf;
        }

        mfem::real_t GetElementEnergy(const mfem::Array<const mfem::FiniteElement *> &el,
                                      mfem::ElementTransformation &Tr,
                                      const mfem::Array<const mfem::Vector *> &elfun) override;

        /// Perform the local action of the NonlinearFormIntegrator
        void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> &el,
                                   mfem::ElementTransformation &Tr,
                                   const mfem::Array<const mfem::Vector *> &elfun,
                                   const mfem::Array<mfem::Vector *> &elvec) override;

        /// Assemble the local gradient matrix
        void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement *> &el,
                                 mfem::ElementTransformation &Tr,
                                 const mfem::Array<const mfem::Vector *> &elfun,
                                 const mfem::Array2D<mfem::DenseMatrix *> &elmats) override;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture nonlinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}

#endif
