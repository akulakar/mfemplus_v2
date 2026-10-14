//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "bilininteg.hpp"
#include "fracture_utils.hpp"

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Bilinear form integrators for regularized variational fracture.
    //
    // All integrators use Voigt notation. The element stiffness is \int B^T C B, where B is the strain-displacement
    // matrix (3 rows in 2D, 6 rows in 3D) and C the stiffness in Voigt form.
    // Isotropic integrators evaluate E and nu at the first integration point only (constant per element).
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Element matrix of the damaged isotropic elasticity form: sum over quadrature points of
    // w * ((1 - d)^2 + k_eps) * B^T C B, with C built from E and nu (first integration point) and d interpolated from the
    // damage grid function. Returns a (dof * dim) x (dof * dim) matrix.
    void DamagedIsotropicElasticityIntegrator::AssembleElementMatrix(
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

        // Damage at the element dofs, interpolated to each quadrature point below.
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
            Mult(dshape, Trans.InverseJacobian(), gshape); // Physical shape function gradients.

            if (i == 0)
            {
                detail::FillIsotropicStiffness(young_mod->Eval(Trans, ip), poisson_ratio->Eval(Trans, ip), dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            Mult(C, B, CB);
            MultAtB(B, CB, elmat_input);
            elmat.Add(w * degradation, elmat_input);
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

        real_t E(0.0), NU(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            Mult(dshape, Trans.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip);
                detail::FillIsotropicStiffness(E, NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            // Twice the total energy: C_{ijkl} u_{k,l} u_{i,j} - 2 p : strain + 2 (quadratic pressure energy).
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);

            AddMult_a_VVt(w * strain_energy, shape, elmat);
        }
    }

    // Same as IsotropicStrainEnergyDamageIntegrator, but the energy density includes the volumetric pressure:
    // 2W = strain energy - 2 (pressure . strain) + 2 (quadratic pressure energy). The pressure terms come from
    // detail::SetPressureTerms and are evaluated at each integration point.
    void IsotropicStrainEnergyPressureDamageIntegrator::AssembleElementMatrix(
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
        body_pressure.SetSize(str_comp);
        C = 0.0;
        B = 0.0;
        body_pressure = 0.0;

        real_t E(0.0), NU(0.0), quadratic_pressure_energy(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            Mult(dshape, Trans.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip);
                detail::FillIsotropicStiffness(E, NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            detail::SetPressureTerms(dim, planeApprox, volumetric_pressure->Eval(Trans, ip), E, NU, body_pressure, quadratic_pressure_energy);

            // Twice the total energy: C_{ijkl} u_{k,l} u_{i,j} - 2 p : strain + 2 (quadratic pressure energy).
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);
            const real_t pressure_energy = InnerProduct(body_pressure, Bu);
            const real_t total_energy = strain_energy - 2.0 * pressure_energy + 2.0 * quadratic_pressure_energy;

            AddMult_a_VVt(w * total_energy, shape, elmat);
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
        identity.SetSize(str_comp);
        volumetric_stress.SetSize(str_comp);
        C = 0.0;
        B = 0.0;
        volumetric_stress = 0.0;
        detail::FillVoigtIdentity(dim, identity);

        real_t E(0.0), NU(0.0), quadratic_volumetric_strain_energy(0.0), total_energy(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            const real_t w = ip.weight * Trans.Weight();
            Mult(dshape, Trans.InverseJacobian(), gshape);

            const real_t vol_strain_coeff = volumetric_strain_coeff->Eval(Trans, ip);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip);
                detail::FillIsotropicStiffness(E, NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            // Isotropic eigenstrain alpha * I gives the volumetric stress C (alpha I).
            detail::VolumetricStress(C, identity, vol_strain_coeff, volumetric_stress);

            if (dim == 2)
            {
                // 2D is not verified. Plane strain has no quadratic term.
                if (planeApprox == 1) // plane stress
                {
                    volumetric_stress *= detail::PlaneStressVolumetricFactor(NU);
                    quadratic_volumetric_strain_energy = 0.5 * vol_strain_coeff * InnerProduct(volumetric_stress, identity) * ((1.0 - 2.0 * NU) * (1.0 - 2.0 * NU) / (E * (1.0 - NU)));
                }
            }
            else
            {
                // (1/2) alpha^2 I^T C I. Note volumetric_stress already contains one factor of alpha.
                quadratic_volumetric_strain_energy = 0.5 * vol_strain_coeff * InnerProduct(volumetric_stress, identity);
            }

            // Twice the energy density, (strain - alpha I)^T C (strain - alpha I), unless the Gershgorin check on the
            // effective stress fails. In that case the energy is scaled by 1e-20 (the value kept from the previous
            // integration point, 0 at the first one), i.e. effectively dropped.
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);
            if (detail::EffectiveStressBoundIsNegative(dim, CBu, volumetric_stress))
            {
                total_energy *= 1e-20;
            }
            else
            {
                const real_t volumetric_strain_energy = InnerProduct(volumetric_stress, Bu);
                total_energy = strain_energy - 2.0 * volumetric_strain_energy + 2.0 * quadratic_volumetric_strain_energy;
            }

            AddMult_a_VVt(w * total_energy, shape, elmat);
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
