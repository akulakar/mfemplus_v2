//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "lininteg.hpp"
#include "fracture_utils.hpp"

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Linear form integrators for regularized variational fracture.
    //
    // Damage equation right hand sides (scalar damage space) integrate twice the energy density against the shape
    // functions. Body force right hand sides (vector displacement space) integrate B^T (stress) scaled by the degradation.
    // Isotropic integrators evaluate E and nu at the first integration point only (constant per element).
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Element vector of the damage equation: sum over quadrature points of w * 2W * phi_i, with
    // 2W = u^T B^T C B u twice the strain energy density of the frozen displacement.
    void IsotropicStrainEnergyDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        shape.SetSize(dof);
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        elvect.SetSize(dof);
        elvect = 0.0;

        detail::GatherElementValues(*disp_gf, *disp_fes, Tr.ElementNo, eldofs, eldofdisp);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

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
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                detail::FillIsotropicStiffness(young_mod->Eval(Tr, ip), poisson_ratio->Eval(Tr, ip), dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            // Twice the strain energy, C_{ijkl} u_{k,l} u_{i,j}.
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);

            add(elvect, w * strain_energy, shape, elvect);
        }
    }

    // Same as IsotropicStrainEnergyDamageLFIntegrator, but 2W includes the volumetric pressure:
    // 2W = strain energy - 2 (pressure . strain) + 2 (quadratic pressure energy), see detail::SetPressureTerms.
    void IsotropicStrainEnergyPressureDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        shape.SetSize(dof);
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        elvect.SetSize(dof);
        elvect = 0.0;

        detail::GatherElementValues(*disp_gf, *disp_fes, Tr.ElementNo, eldofs, eldofdisp);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

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
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip);
                detail::FillIsotropicStiffness(E, NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            // Volumetric pressure, evaluated at each integration point.
            detail::SetPressureTerms(dim, planeApprox, volumetric_pressure->Eval(Tr, ip), E, NU, body_pressure, quadratic_pressure_energy);

            // Twice the total energy: C_{ijkl} u_{k,l} u_{i,j} - 2 p : strain + 2 (quadratic pressure energy).
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);
            const real_t pressure_energy = InnerProduct(body_pressure, Bu);
            const real_t total_energy = strain_energy - 2.0 * pressure_energy + 2.0 * quadratic_pressure_energy;

            add(elvect, w * total_energy, shape, elvect);
        }
    }

    // Element vector of the damage equation with a volumetric eigenstrain alpha I:
    // sum over quadrature points of w * 2W * phi_i, with 2W = (strain - alpha I)^T C (strain - alpha I), including the
    // same Gershgorin-type switch (energy scaled by 1e-20 where it fails) as the matching bilinear integrator.
    void IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        shape.SetSize(dof);
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        elvect.SetSize(dof);
        elvect = 0.0;

        detail::GatherElementValues(*disp_gf, *disp_fes, Tr.ElementNo, eldofs, eldofdisp);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

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
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            Mult(dshape, Tr.InverseJacobian(), gshape);

            const real_t vol_strain_coeff = volumetric_strain_coeff->Eval(Tr, ip);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip);
                detail::FillIsotropicStiffness(E, NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            // Isotropic eigenstrain alpha * I gives the volumetric stress C (alpha I).
            detail::VolumetricStress(C, identity, vol_strain_coeff, volumetric_stress);

            // 2D is not verified: no quadratic term, and the volumetric stress is not corrected for plane stress here.
            if (dim == 3)
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

            add(elvect, w * total_energy, shape, elvect);
        }
    }

    // Element vector of the damage equation with a general stiffness: sum over quadrature points of
    // w * (u^T B^T C B u) * phi_i, with C evaluated at each integration point.
    void AnisotropicElasticityDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        shape.SetSize(dof);
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        elvect.SetSize(dof);
        elvect = 0.0;

        detail::GatherElementValues(*disp_gf, *disp_fes, Tr.ElementNo, eldofs, eldofdisp);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

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
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            Mult(dshape, Tr.InverseJacobian(), gshape);

            stiffness->Eval(C, Tr, ip); // Stiffness is evaluated at each integration point.
            detail::FillStrainDisplacement(gshape, dim, B);

            // Twice the strain energy, C_{ijkl} u_{k,l} u_{i,j}.
            const real_t strain_energy = detail::TwiceStrainEnergy(C, B, eldofdisp, Bu, CBu);

            add(elvect, w * strain_energy, shape, elvect);
        }
    }

    // Element vector of the pressure body force in the elasticity equation: sum over quadrature points of
    // w * ((1 - d)^2 + k_eps) * B^T p_voigt, where p_voigt holds the pressure on the normal components. Returns a
    // (dof * dim) vector. In 2D, only plane stress assembles a non-zero load.
    void DamagedPressureBodyForceLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elvect.SetSize(dof * dim);
        elvec_input.SetSize(dof * dim);
        elvect = 0.0;

        detail::GatherElementValues(*damage_gf, *damage_fes, Tr.ElementNo, eldofs, eldofdamage);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

        B.SetSize(str_comp, dof * dim);
        body_pressure.SetSize(str_comp);
        B = 0.0;
        body_pressure = 0.0;

        real_t NU(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            const real_t degradation = detail::Degradation(InnerProduct(shape, eldofdamage), k_epsilon);
            const real_t pressure_coeff = pressure->Eval(Tr, ip);
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            if (dim == 2)
            {
                // 2D is not verified: only plane stress is handled; the plane strain load stays zero.
                if (planeApprox == 1)
                {
                    body_pressure(0) = body_pressure(1) = pressure_coeff * detail::PlaneStressVolumetricFactor(NU);
                }
            }
            else
            {
                body_pressure(0) = body_pressure(1) = body_pressure(2) = pressure_coeff;
            }

            B.MultTranspose(body_pressure, elvec_input);
            add(elvect, w * degradation, elvec_input, elvect);
        }
    }

    // Element vector of the eigenstrain body force in the elasticity equation: sum over quadrature points of
    // w * ((1 - d)^2 + k_eps) * B^T (C alpha I), with C built from E and nu (first integration point). Returns a
    // (dof * dim) vector.
    void DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elvect.SetSize(dof * dim);
        elvec_input.SetSize(dof * dim);
        elvect = 0.0;

        detail::GatherElementValues(*damage_gf, *damage_fes, Tr.ElementNo, eldofs, eldofdamage);

        const IntegrationRule &ir = detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr);

        C.SetSize(str_comp, str_comp);
        B.SetSize(str_comp, dof * dim);
        identity.SetSize(str_comp);          // (1 1 1 0 0 0)
        volumetric_stress.SetSize(str_comp); // Stiffness acting on identity, scaled by the strain coefficient
        C = 0.0;
        B = 0.0;
        volumetric_stress = 0.0;
        detail::FillVoigtIdentity(dim, identity);

        real_t NU(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            const real_t w = ip.weight * Tr.Weight();
            const real_t degradation = detail::Degradation(InnerProduct(shape, eldofdamage), k_epsilon);
            const real_t vol_strain_coeff = volumetric_strain_coeff->Eval(Tr, ip);
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                detail::FillIsotropicStiffness(young_mod->Eval(Tr, ip), NU, dim, planeApprox, C);
            }
            detail::FillStrainDisplacement(gshape, dim, B);

            detail::VolumetricStress(C, identity, vol_strain_coeff, volumetric_stress);
            if (dim == 2 && planeApprox == 1)
            {
                volumetric_stress *= detail::PlaneStressVolumetricFactor(NU); // 2D is not verified.
            }

            B.MultTranspose(volumetric_stress, elvec_input);
            add(elvect, w * degradation, elvec_input, elvect);
        }
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture linear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}
