//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Strain energy evaluators for the damage (phase field) equation.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/SolidMechanics/Fracture/energy.hpp"
#include "mfemplus/SolidMechanics/Fracture/fracture_utils.hpp"

using namespace mfem;

namespace mfemplus
{
    void StrainEnergyEvaluator::PrepareElement(const FiniteElement &el, ElementTransformation &Tr)
    {
        const int dof = el.GetDof();
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        detail::GatherElementValues(*disp_gf, *disp_fes, Tr.ElementNo, eldofs, eldofdisp);
    }

    // Loops over the elements of the QuadratureSpace of W and stores 2W at each quadrature point. The element
    // transformation comes from the displacement space, whose mesh must be the one of the QuadratureSpace.
    void StrainEnergyEvaluator::Evaluate(QuadratureFunction &W)
    {
        const QuadratureSpaceBase &qs = *W.GetSpace();
        Vector element_energy;
        for (int e = 0; e < qs.GetNE(); e++)
        {
            const IntegrationRule &ir = qs.GetIntRule(e);
            element_energy.SetSize(ir.GetNPoints());
            ElementEnergy(*disp_fes->GetFE(e), *disp_fes->GetElementTransformation(e), ir, element_energy);

            Vector stored;
            W.GetValues(e, stored); // Reference to the values of element e.
            stored = element_energy;
        }
    }

    // 2W = u^T B^T C B u, evaluated with the closed form isotropic stiffness.
    void IsotropicStrainEnergyEvaluator::ElementEnergy(
        const FiniteElement &el, ElementTransformation &Tr, const IntegrationRule &ir, Vector &energy)
    {
        const int dim = el.GetDim();
        PrepareElement(el, Tr);

        real_t lambda(0.0), mu(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                detail::IsotropicLame(young_mod->Eval(Tr, ip), poisson_ratio->Eval(Tr, ip), dim, planeApprox, lambda, mu);
            }
            energy(i) = detail::TwiceIsotropicStrainEnergy(gshape, dim, eldofdisp, lambda, mu, Bu, CBu);
        }
    }

    // 2W = strain energy - 2 (pressure . strain) + 2 (quadratic pressure energy), see detail::SetPressureTerms.
    void IsotropicStrainEnergyPressureEvaluator::ElementEnergy(
        const FiniteElement &el, ElementTransformation &Tr, const IntegrationRule &ir, Vector &energy)
    {
        const int dim = el.GetDim();
        PrepareElement(el, Tr);
        body_pressure.SetSize(detail::NumStrainComponents(dim));
        body_pressure = 0.0;

        real_t E(0.0), NU(0.0), lambda(0.0), mu(0.0), quadratic_pressure_energy(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            Mult(dshape, Tr.InverseJacobian(), gshape);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip);
                detail::IsotropicLame(E, NU, dim, planeApprox, lambda, mu);
            }
            detail::SetPressureTerms(dim, planeApprox, volumetric_pressure->Eval(Tr, ip), E, NU, body_pressure, quadratic_pressure_energy);

            const real_t strain_energy = detail::TwiceIsotropicStrainEnergy(gshape, dim, eldofdisp, lambda, mu, Bu, CBu);
            const real_t pressure_energy = InnerProduct(body_pressure, Bu);
            energy(i) = strain_energy - 2.0 * pressure_energy + 2.0 * quadratic_pressure_energy;
        }
    }

    // 2W = (strain - alpha I)^T C (strain - alpha I) = strain energy - 2 (volumetric stress . strain) + alpha^2 I^T C I,
    // with volumetric stress = C (alpha I). If the Gershgorin check on the effective stress fails at a point, the
    // energy there is scaled by 1e-20 (the value kept from the previous integration point, 0 at the first one).
    // The check is skipped if it is disabled.
    void IsotropicStrainEnergyVolumetricEigenStrainEvaluator::ElementEnergy(
        const FiniteElement &el, ElementTransformation &Tr, const IntegrationRule &ir, Vector &energy)
    {
        const int dim = el.GetDim();
        const int str_comp = detail::NumStrainComponents(dim);
        PrepareElement(el, Tr);
        identity.SetSize(str_comp);
        volumetric_stress.SetSize(str_comp);
        volumetric_stress = 0.0;
        detail::FillVoigtIdentity(dim, identity);

        real_t E(0.0), NU(0.0), lambda(0.0), mu(0.0), quadratic_volumetric_strain_energy(0.0), total_energy(0.0);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            Mult(dshape, Tr.InverseJacobian(), gshape);

            const real_t vol_strain_coeff = volumetric_strain_coeff->Eval(Tr, ip);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip);
                detail::IsotropicLame(E, NU, dim, planeApprox, lambda, mu);
            }

            // Isotropic eigenstrain alpha * I gives the volumetric stress C (alpha I).
            detail::IsotropicVolumetricStress(dim, lambda, mu, vol_strain_coeff, volumetric_stress);

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

            const real_t strain_energy = detail::TwiceIsotropicStrainEnergy(gshape, dim, eldofdisp, lambda, mu, Bu, CBu);
            if (gershgorin_check && detail::EffectiveStressBoundIsNegative(dim, CBu, volumetric_stress))
            {
                total_energy *= 1e-20;
            }
            else
            {
                const real_t volumetric_strain_energy = InnerProduct(volumetric_stress, Bu);
                total_energy = strain_energy - 2.0 * volumetric_strain_energy + 2.0 * quadratic_volumetric_strain_energy;
            }
            energy(i) = total_energy;
        }
    }
}
