//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Nonlinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/SolidMechanics/Fracture/nonlininteg.hpp"
#include "mfemplus/SolidMechanics/Fracture/fracture_utils.hpp"

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Nonlinear form integrators for regularized variational fracture.
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Get RVFT energy for a given state. May be used in modified Newton methods.
    real_t RVFT_AT2_Integrator::GetElementEnergy(const Array<const FiniteElement *> &el, ElementTransformation &Tr, const Array<const Vector *> &elfun)
    {
        if (el.Size() != 2)
            mfem_error("RVFT_AT2_Integrator::GetElementEnergy has incorrect block finite element space size!");

        // Blocks: 0 = displacement (vector), 1 = damage (scalar).
        const int dof_u = el[0]->GetDof();
        const int dof_d = el[1]->GetDof();
        const int dim = el[0]->GetDim();

        // Psi = g(phi) W(u) + gOverEps/2 phi^2 + gTimesEps/2 |grad phi|^2, with W = 0.5 * (2W) from the isotropic kernel.
        const IntegrationRule &ir = detail::IntRuleOrDefault(nullptr, *el[0], Tr);

        dshape_u.SetSize(dof_u, dim);
        gshape_u.SetSize(dof_u, dim);
        dshape_d.SetSize(dof_d, dim);
        shape_d.SetSize(dof_d);
        grad_d.SetSize(dim);

        real_t lambda(0.0), mu(0.0);
        real_t energy = 0.0;
        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            Tr.SetIntPoint(&ip);

            el[0]->CalcDShape(ip, dshape_u);
            Mult(dshape_u, Tr.InverseJacobian(), gshape_u);
            el[1]->CalcPhysShape(Tr, shape_d);
            el[1]->CalcPhysDShape(Tr, dshape_d);

            if (i == 0)
            {
                // E and nu are assumed constant per element; plane strain in 2D.
                detail::IsotropicLame(Ey_coeff->Eval(Tr, ip), nu_coeff->Eval(Tr, ip), dim, 0, lambda, mu);
            }

            const real_t phi = shape_d * *elfun[1];
            dshape_d.MultTranspose(*elfun[1], grad_d);

            const real_t W = 0.5 * detail::TwiceIsotropicStrainEnergy(gshape_u, dim, *elfun[0], lambda, mu, strain, stress);
            const real_t psi = detail::Degradation(phi, k_eps) * W + 0.5 * gOverEps_coeff->Eval(Tr, ip) * phi * phi + 0.5 * gTimesEps_coeff->Eval(Tr, ip) * (grad_d * grad_d);

            energy += ip.weight * Tr.Weight() * psi;
        }
        return energy;
    }

    // Evaluates action of RVFT AT2 weak form on a given state
    void RVFT_AT2_Integrator::AssembleElementVector(const Array<const FiniteElement *> &el,
                                                    ElementTransformation &Tr,
                                                    const Array<const Vector *> &elfun,
                                                    const Array<Vector *> &elvec)
    {
    }

    // Assembled the element Jacobian
    void RVFT_AT2_Integrator::AssembleElementGrad(const Array<const FiniteElement *> &el,
                                                  ElementTransformation &Tr,
                                                  const Array<const Vector *> &elfun,
                                                  const Array2D<DenseMatrix *> &elmats)
    {
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture nonlinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}
