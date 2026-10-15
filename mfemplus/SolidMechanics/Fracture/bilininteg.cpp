//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/SolidMechanics/Fracture/bilininteg.hpp"
#include "mfemplus/SolidMechanics/Fracture/fracture_utils.hpp"

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Bilinear form integrators for regularized variational fracture.
    //
    // The anisotropic integrators use Voigt notation: the element stiffness is \int B^T C B, where B is the
    // strain-displacement matrix (3 rows in 2D, 6 rows in 3D) and C the stiffness in Voigt form.
    // The isotropic integrators use the closed form C = lambda I(x)I + 2 mu I_sym and never form B or C (see
    // detail::AddIsotropicStiffness and detail::TwiceIsotropicStrainEnergy); the Voigt vectors keep the same ordering.
    // Isotropic integrators evaluate E and nu at the first integration point only (constant per element).
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Element matrix of the damaged isotropic elasticity form: sum over quadrature points of
    // w * ((1 - d)^2 + k_eps) * B^T C B, with C built from E and nu (first integration point) and d interpolated from the
    // damage grid function. B^T C B is evaluated with the closed form isotropic stiffness. Returns a (dof * dim) x (dof * dim)
    // matrix.
    void DamagedIsotropicElasticityIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof * dim);
        elmat = 0.0;

        // Damage at the element dofs, interpolated to each quadrature point below.
        detail::GatherElementValues(*damage_gf, *damage_fes, Trans.ElementNo, eldofs, eldofdamage);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans);

        real_t lambda(0.0), mu(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            const real_t degradation = detail::Degradation(InnerProduct(shape, eldofdamage), k_epsilon);
            Mult(dshape, Trans.InverseJacobian(), gshape); // Physical shape function gradients.

            if (i == 0)
            {
                detail::IsotropicLame(young_mod->Eval(Trans, ip), poisson_ratio->Eval(Trans, ip), dim, planeApprox, lambda, mu);
            }
            detail::AddIsotropicStiffness(gshape, dim, w * degradation, lambda, mu, elmat);
        }
    }

    // Same as DamagedIsotropicElasticityIntegrator, but C is the user supplied MatrixCoefficient evaluated at each
    // integration point.
    void DamagedAnisotropicElasticityIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof * dim);
        elmat = 0.0;

        detail::GatherElementValues(*damage_gf, *damage_fes, Trans.ElementNo, eldofs, eldofdamage);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans);

        C.SetSize(str_comp, str_comp);
        B.SetSize(str_comp, dof * dim);
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        C = 0.0;
        B = 0.0;

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            const real_t degradation = detail::Degradation(InnerProduct(shape, eldofdamage), k_epsilon);
            Mult(dshape, Trans.InverseJacobian(), gshape);

            stiffness->Eval(C, Trans, ip); // Stiffness is evaluated at each integration point.
            detail::FillStrainDisplacement(gshape, dim, B);

            Mult(C, B, CB);
            MultAtB(B, CB, elmat_input);
            elmat.Add(w * degradation, elmat_input);
        }
    }

    // Element matrix of the damage equation: sum over quadrature points of w * (u^T B^T C B u) * phi_i phi_j, i.e. twice
    // the strain energy density of the frozen displacement times the mass-like term. Returns a dof x dof matrix.
    void IsotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans), el, Trans, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            AddMult_a_VVt(ip.weight * Trans.Weight() * energy(i), shape, elmat);
        }
    }

    // Same as IsotropicStrainEnergyDamageIntegrator, but the energy density includes the volumetric pressure:
    // 2W = strain energy - 2 (pressure . strain) + 2 (quadratic pressure energy). The pressure terms come from
    // detail::SetPressureTerms and are evaluated at each integration point.
    void IsotropicStrainEnergyPressureDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans), el, Trans, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            AddMult_a_VVt(ip.weight * Trans.Weight() * energy(i), shape, elmat);
        }
    }

    // Element matrix of the damage equation with a volumetric eigenstrain alpha I:
    // 2W = (strain - alpha I)^T C (strain - alpha I) = strain energy - 2 (volumetric stress . strain) + alpha^2 I^T C I,
    // with volumetric stress = C (alpha I). If the Gershgorin check on the effective stress fails at a point, the
    // energy there is scaled by 1e-20 (effectively dropped). Returns a dof x dof matrix.
    void IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans), el, Trans, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            AddMult_a_VVt(ip.weight * Trans.Weight() * energy(i), shape, elmat);
        }
    }

    // Same as IsotropicStrainEnergyDamageIntegrator, but C is the user supplied MatrixCoefficient evaluated at each
    // integration point.
    void AnisotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        detail::GatherElementValues(*disp_gf, *disp_fes, Trans.ElementNo, eldofs, eldofdisp);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Trans), el, Trans);

        C.SetSize(str_comp, str_comp);
        B.SetSize(str_comp, dof * dim);
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        C = 0.0;
        B = 0.0;

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            Mult(dshape, Trans.InverseJacobian(), gshape);

            stiffness->Eval(C, Trans, ip); // Stiffness is evaluated at each integration point.
            detail::FillStrainDisplacement(gshape, dim, B);

            // Twice the strain energy, C_{ijkl} u_{k,l} u_{i,j}.
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);

            AddMult_a_VVt(w * strain_energy, shape, elmat);
        }
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}
