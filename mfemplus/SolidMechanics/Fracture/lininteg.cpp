//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form integrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/SolidMechanics/Fracture/lininteg.hpp"
#include "mfemplus/SolidMechanics/Fracture/fracture_utils.hpp"

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
    // The isotropic integrators use the closed form C = lambda I(x)I + 2 mu I_sym and never form B or C; only the
    // anisotropic one uses the Voigt B and C matrices.
    // Isotropic integrators evaluate E and nu at the first integration point only (constant per element).
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Element vector of the damage equation: sum over quadrature points of w * 2W * phi_i, with
    // 2W = u^T B^T C B u twice the strain energy density of the frozen displacement.
    void IsotropicStrainEnergyDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elvect.SetSize(dof);
        elvect = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr), el, Tr, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            add(elvect, ip.weight * Tr.Weight() * energy(i), shape, elvect);
        }
    }

    // Same as IsotropicStrainEnergyDamageLFIntegrator, but 2W includes the volumetric pressure:
    // 2W = strain energy - 2 (pressure . strain) + 2 (quadratic pressure energy), see detail::SetPressureTerms.
    void IsotropicStrainEnergyPressureDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elvect.SetSize(dof);
        elvect = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr), el, Tr, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            add(elvect, ip.weight * Tr.Weight() * energy(i), shape, elvect);
        }
    }

    // Element vector of the damage equation with a volumetric eigenstrain alpha I:
    // sum over quadrature points of w * 2W * phi_i, with 2W = (strain - alpha I)^T C (strain - alpha I), including the
    // same Gershgorin-type switch (energy scaled by 1e-20 where it fails) as the matching bilinear integrator.
    void IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        const int dof = el.GetDof();

        shape.SetSize(dof);
        elvect.SetSize(dof);
        elvect = 0.0;

        // Twice the energy density at the quadrature points: computed here or read from the shared QuadratureFunction.
        const IntegrationRule &ir = energy_source.GetEnergyAndRule(detail::IntRuleOrDefault(GetIntegrationRule(el, Tr), el, Tr), el, Tr, energy);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);

            add(elvect, ip.weight * Tr.Weight() * energy(i), shape, elvect);
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

        body_pressure.SetSize(str_comp);
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

            detail::StrainDisplacementTransposeNormal(gshape, dim, body_pressure, elvec_input);
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

        volumetric_stress.SetSize(str_comp); // Stiffness acting on identity, scaled by the strain coefficient
        volumetric_stress = 0.0;

        real_t NU(0.0), lambda(0.0), mu(0.0);

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
                detail::IsotropicLame(young_mod->Eval(Tr, ip), NU, dim, planeApprox, lambda, mu);
            }

            detail::IsotropicVolumetricStress(dim, lambda, mu, vol_strain_coeff, volumetric_stress);
            if (dim == 2 && planeApprox == 1)
            {
                volumetric_stress *= detail::PlaneStressVolumetricFactor(NU); // 2D is not verified.
            }

            detail::StrainDisplacementTransposeNormal(gshape, dim, volumetric_stress, elvec_input);
            add(elvect, w * degradation, elvec_input, elvect);
        }
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture linear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}
