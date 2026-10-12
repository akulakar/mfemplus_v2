//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Linear form intergrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "lininteg.hpp"
#include <cmath>
#include <Dense>
#include <SVD>

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Linear form integrators for regularized variational fracture.
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    void FractureDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;

        shape.SetSize(dof); // vector of size dof
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement
        eldofdamage.SetSize(dof);     // scalar valued damage

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Tr);

        if (ir == NULL)
        {
            ir = &IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
        }
        double w, NU, E;

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);

        body_pressure.SetSize(str_comp);
        elvect.SetSize(dof);

        C = 0.0;
        B = 0.0;
        body_pressure = 0.0;
        elvect = 0.0;

        double lambda1, lambda2, lambda3;
        double strain_energy(0.0), pressure_energy(0.0), quadratic_pressure_energy(0.0), total_energy(0.0);
        double pressure_coeff(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
            } // Constant throughout element

            Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                // volumetric pressure
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Tr, ip);

                    switch (planeApprox)
                    {
                    case 1:
                        // Plane stress
                        body_pressure(0) = pressure_coeff * (2.0 - (1.0 / (1.0 - NU)));
                        body_pressure(1) = pressure_coeff * (2.0 - (1.0 / (1.0 - NU)));

                        quadratic_pressure_energy = (pressure_coeff * pressure_coeff * pow(1.0 - 2.0 * NU, 2.0)) / (E * (1.0 - NU));
                        break;
                    }
                }

                if (i == 0)
                {
                    switch (planeApprox)
                    {
                    case 0:
                        // Plane strain
                        C(0, 0) = C(1, 1) = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
                        C(0, 1) = C(1, 0) = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
                        C(2, 2) = E / (2.0 * (1.0 + NU));
                        break;

                    case 1:
                        // Plane stress
                        C(0, 0) = C(1, 1) = (E / (1.0 - pow(NU, 2.0)));
                        C(0, 1) = C(1, 0) = (E * NU / (1.0 - pow(NU, 2.0)));
                        C(2, 2) = E / (2.0 * (1.0 + NU));
                        break;
                    }
                }

                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Tr, ip);
                    body_pressure(0) = pressure_coeff;
                    body_pressure(1) = pressure_coeff;
                    body_pressure(2) = pressure_coeff;

                    quadratic_pressure_energy = (3.0 * pressure_coeff * pressure_coeff) / ((2.0 * E) / (1.0 - 2.0 * NU));
                }

                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }

                // In 3D, we have 6 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = gshape(spf, 2);
                    B(3, spf + 2 * dof) = gshape(spf, 1);
                    B(4, spf) = gshape(spf, 2);
                    B(4, spf + 2 * dof) = gshape(spf, 0);
                    B(5, spf) = gshape(spf, 1);
                    B(5, spf + dof) = gshape(spf, 0);
                }
                break;
            }

            // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
            // This is equivalent to.
            Mult(C, B, CB);          // CB is 6 x (dof * dim)
            CB.Mult(eldofdisp, CBu); // CBu has dimension strain_comps. This is the stress vector.
            // For tensile loading, no need for Gershgorin check.
            B.Mult(eldofdisp, Bu);                             // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = InnerProduct(CBu, Bu);             // This is twice the strain energy
                                                               // if (volumetric_pressure != nullptr)
                                                               // {
            pressure_energy = InnerProduct(body_pressure, Bu); // This is the linear pressure energy.
            // }
            total_energy = strain_energy - 2.0 * pressure_energy + 2.0 * quadratic_pressure_energy; // This is twice the total energy.

            // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            // if (dim == 2)
            // {
            //     // In 2D lambda min is lambda1.
            //     lambda1 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda2 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda3 = lambda1 + 1.0;

            //     // lambda1 = CBu(0) - std::abs(CBu(3));
            //     // lambda2 = CBu(1) - std::abs(CBu(3));
            //     // lambda3 = lambda1 + lambda2; // artificially making it greater than both.
            // }
            // else if (dim == 3)
            // {
            //     lambda1 = CBu(0) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
            //     lambda2 = CBu(1) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
            //     lambda3 = CBu(2) - pressure_coeff - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            // }

            // if (std::min({lambda1, lambda2, lambda3}) < 0)
            // {
            //     total_energy = 0.0;
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            // for now okay, but probably will change it to element average strain energy.
            add(elvect, w * total_energy, shape, elvect);
        }
    };

    void IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;

        shape.SetSize(dof); // vector of size dof
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement
        eldofdamage.SetSize(dof);     // scalar valued damage
        elvect.SetSize(dof);
        elvect = 0.0;

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Tr);

        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }
        double w, NU, E;

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        identity.SetSize(str_comp);
        volumetric_stress.SetSize(str_comp);

        C = 0.0;
        B = 0.0;

        identity = 0.0;
        identity(0) = 1.0;
        identity(1) = 1.0;
        if (dim == 3)
            identity(2) = 1.0;

        volumetric_stress = 0.0;

        // double lambda1(0.0), lambda2(0.0), lambda3(0.0);
        double vol_strain_coeff(0.0);
        double strain_energy(0.0), volumetric_strain_energy(0.0), quadratic_volumetric_strain_energy(0.0), total_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            vol_strain_coeff = volumetric_strain_coeff->Eval(Tr, ip);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
            } // Constant throughout element

            Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                if (i == 0)
                {
                    switch (planeApprox)
                    {
                    case 0:
                        // Plane strain
                        C(0, 0) = C(1, 1) = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
                        C(0, 1) = C(1, 0) = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
                        C(2, 2) = E / (2.0 * (1.0 + NU));
                        break;

                    case 1:
                        // Plane stress
                        C(0, 0) = C(1, 1) = (E / (1.0 - pow(NU, 2.0)));
                        C(0, 1) = C(1, 0) = (E * NU / (1.0 - pow(NU, 2.0)));
                        C(2, 2) = E / (2.0 * (1.0 + NU));
                        break;
                    }
                }

                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }

                // volumetric stress. Let stiffness C act on identity and multiply by strain coeff.
                C.Mult(identity, volumetric_stress);
                volumetric_stress *= vol_strain_coeff;

                break;

            case 3:
                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }

                // In 3D, we have 6 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = gshape(spf, 2);
                    B(3, spf + 2 * dof) = gshape(spf, 1);
                    B(4, spf) = gshape(spf, 2);
                    B(4, spf + 2 * dof) = gshape(spf, 0);
                    B(5, spf) = gshape(spf, 1);
                    B(5, spf + dof) = gshape(spf, 0);
                }

                // volumetric stress. Let stiffness C act on identity and multiply by strain coeff.
                C.Mult(identity, volumetric_stress);
                volumetric_stress *= vol_strain_coeff;
                quadratic_volumetric_strain_energy = (1.0 / 2.0) * vol_strain_coeff * InnerProduct(volumetric_stress, identity);

                break;
            }

            // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
            // This is equivalent to.
            Mult(C, B, CB);          // CB is 6 x (dof * dim)
            CB.Mult(eldofdisp, CBu); // CBu has dimension strain_comps. This is the stress vector.
            // For tensile loading, no need for Gershgorin check.
            B.Mult(eldofdisp, Bu);                                                                                    // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = InnerProduct(CBu, Bu);                                                                    // This is twice the strain energy
            volumetric_strain_energy = InnerProduct(volumetric_stress, Bu);                                           // This is the linear pressure energy
            total_energy = strain_energy - 2.0 * volumetric_strain_energy + 2.0 * quadratic_volumetric_strain_energy; // This is twice the total energy

            // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            // if (dim == 2)
            // {
            //     // In 2D lambda min is lambda1.
            //     lambda1 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda2 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda3 = lambda1 + 1.0;

            //     // lambda1 = CBu(0) - std::abs(CBu(3));
            //     // lambda2 = CBu(1) - std::abs(CBu(3));
            //     // lambda3 = lambda1 + lambda2; // artificially making it greater than both.
            // }
            // else if (dim == 3)
            // {
            //     lambda1 = CBu(0) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
            //     lambda2 = CBu(1) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
            //     lambda3 = CBu(2) - pressure_coeff - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            // }

            // if (std::min({lambda1, lambda2, lambda3}) < 0)
            // {
            //     total_energy = 0.0;
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            // for now okay, but probably will change it to element average strain energy.
            add(elvect, w * total_energy, shape, elvect);
        }
    };

    void AnisotropicElasticityDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;

        shape.SetSize(dof); // vector of size dof
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement
        eldofdamage.SetSize(dof);     // scalar valued damage

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Tr);

        if (ir == NULL)
        {
            ir = &IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
        }
        double w;

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);

        elvect.SetSize(dof);

        C = 0.0;
        B = 0.0;
        elvect = 0.0;

        double lambda1, lambda2, lambda3;
        double strain_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            stiffness->Eval(C, Tr, ip); // The stiffness matrix is evaluated at each integration point.

            Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:
                // In 3D, we have 6 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = gshape(spf, 2);
                    B(3, spf + 2 * dof) = gshape(spf, 1);
                    B(4, spf) = gshape(spf, 2);
                    B(4, spf + 2 * dof) = gshape(spf, 0);
                    B(5, spf) = gshape(spf, 1);
                    B(5, spf + dof) = gshape(spf, 0);
                }
                break;
            }

            // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
            // This is equivalent to.
            Mult(C, B, CB);          // CB is 6 x (dof * dim)
            CB.Mult(eldofdisp, CBu); // CBu has dimension strain_comps. This is the stress vector.
            // For tensile loading, no need for Gershgorin check.
            B.Mult(eldofdisp, Bu);                 // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = InnerProduct(CBu, Bu); // This is twice the strain energy

            // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            // if (dim == 2)
            // {
            //     // In 2D lambda min is lambda1.
            //     lambda1 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda2 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda3 = lambda1 + 1.0;

            //     // lambda1 = CBu(0) - std::abs(CBu(3));
            //     // lambda2 = CBu(1) - std::abs(CBu(3));
            //     // lambda3 = lambda1 + lambda2; // artificially making it greater than both.
            // }
            // else if (dim == 3)
            // {
            //     lambda1 = CBu(0) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
            //     lambda2 = CBu(1) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
            //     lambda3 = CBu(2) - pressure_coeff - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            // }

            // if (std::min({lambda1, lambda2, lambda3}) < 0)
            // {
            //     total_energy = 0.0;
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            // for now okay, but probably will change it to element average strain energy.
            add(elvect, w * strain_energy, shape, elvect);
        }
    };

    void PressureBodyForceDamageLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        {
            int dof = el.GetDof();
            int dim = el.GetDim();
            int str_comp = (dim == 2) ? 3 : 6;
            int elnum = Tr.ElementNo;

            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            shape.SetSize(dof);

            eldofs.SetSize(dof); // scalar for damage
            eldofdamage.SetSize(dof);
            eldofs = 0;
            eldofdamage = 0.0;

            elvect.SetSize(dof * dim);
            B.SetSize(str_comp, dof * dim); // Strain displacement matrix
            body_pressure.SetSize(str_comp);
            elvec_input.SetSize(dof * dim);
            elvect = 0.0;
            B = 0.0;
            body_pressure = 0.0;

            damage_fes->GetElementDofs(elnum, eldofs);

            for (int i = 0; i < eldofdamage.Size(); i++)
            {
                eldofdamage(i) = (*damage_gf)(eldofs[i]);
            }

            // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
            // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

            const IntegrationRule *ir = GetIntegrationRule(el, Tr);
            if (ir == NULL)
            {
                ir = &IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
            }
            double w, NU, E;

            double damage_value(0.0), degradation_constant(0.0);
            double pressure_coeff(0.0);

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);

                el.CalcDShape(ip, dshape);
                Tr.SetIntPoint(&ip);
                el.CalcPhysShape(Tr, shape);

                w = ip.weight * Tr.Weight(); // Quadrature weights

                damage_value = InnerProduct(shape, eldofdamage);
                degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

                pressure_coeff = pressure->Eval(Tr, ip);

                if (i == 0)
                {
                    NU = poisson_ratio->Eval(Tr, ip);
                    E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
                } // Constant throughout element

                Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use Voigt notation to speed up the assembly process.
                // For this, we need the strain displacement matrix B. The element stiffness can be computed as
                // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
                // The B matrix as 3 rows in 2D and 6 rowd in 3D.

                switch (dim)
                {
                case 2:
                    // volumetric pressure. 2 non-zero components
                    switch (planeApprox)
                    {
                    case 1:
                        // plane stress
                        body_pressure(0) = pressure_coeff * (2.0 - (1.0 / (1.0 - NU)));
                        body_pressure(1) = pressure_coeff * (2.0 - (1.0 / (1.0 - NU)));
                        break;
                    }

                    // In 2D, we have 3 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf) = gshape(spf, 1);
                        B(2, spf + dof) = gshape(spf, 0);
                    }
                    break;
                case 3:

                    // volumetric pressure. 3 non zero components
                    body_pressure(0) = pressure_coeff;
                    body_pressure(1) = pressure_coeff;
                    body_pressure(2) = pressure_coeff;

                    // In 3D, we have 6 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf + 2 * dof) = gshape(spf, 2);
                        B(3, spf + dof) = gshape(spf, 2);
                        B(3, spf + 2 * dof) = gshape(spf, 1);
                        B(4, spf) = gshape(spf, 2);
                        B(4, spf + 2 * dof) = gshape(spf, 0);
                        B(5, spf) = gshape(spf, 1);
                        B(5, spf + dof) = gshape(spf, 0);
                    }
                    break;
                }
                B.MultTranspose(body_pressure, elvec_input);
                add(elvect, w * degradation_constant, elvec_input, elvect);
            }
        }
    }

    void DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
    {
        {
            int dof = el.GetDof();
            int dim = el.GetDim();
            int str_comp = (dim == 2) ? 3 : 6;
            int elnum = Tr.ElementNo;

            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            shape.SetSize(dof);
            elvec_input.SetSize(dof * dim);
            elvect = 0.0;

            eldofs.SetSize(dof); // scalar for damage
            eldofdamage.SetSize(dof);
            eldofs = 0;
            eldofdamage = 0.0;

            elvect.SetSize(dof * dim);
            C.SetSize(str_comp, str_comp);       // Stiffness
            B.SetSize(str_comp, dof * dim);      // Strain displacement matrix
            identity.SetSize(str_comp);          // Identity (1 1 1 0 0 0)
            volumetric_stress.SetSize(str_comp); // Stiffness tensor acting on identity multiplied by strain coeff
            C = 0.0;
            B = 0.0;

            identity = 0.0; // Set appropriate components.
            identity(0) = 1.0;
            identity(1) = 1.0;
            if (dim == 3)
                identity(2) = 1.0;

            volumetric_stress = 0.0;

            damage_fes->GetElementDofs(elnum, eldofs);

            for (int i = 0; i < eldofdamage.Size(); i++)
            {
                eldofdamage(i) = (*damage_gf)(eldofs[i]);
            }

            // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
            // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

            const IntegrationRule *ir = GetIntegrationRule(el, Tr);
            if (ir == NULL)
            {
                int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
                ir = &IntRules.Get(el.GetGeomType(), order);
            }
            double w, NU, E;

            double damage_value(0.0), degradation_constant(0.0);
            double vol_strain_coeff(0.0);

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);

                el.CalcDShape(ip, dshape);
                Tr.SetIntPoint(&ip);
                el.CalcPhysShape(Tr, shape);

                w = ip.weight * Tr.Weight(); // Quadrature weights

                damage_value = InnerProduct(shape, eldofdamage);
                degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

                vol_strain_coeff = volumetric_strain_coeff->Eval(Tr, ip);

                if (i == 0)
                {
                    NU = poisson_ratio->Eval(Tr, ip);
                    E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
                }
                // Constant throughout element

                Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use Voigt notation to speed up the assembly process.
                switch (dim)
                {
                case 2:
                    if (i == 0)
                    {
                        switch (planeApprox)
                        {
                        case 0:
                            // Plane strain
                            C(0, 0) = C(1, 1) = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
                            C(0, 1) = C(1, 0) = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
                            C(2, 2) = E / (2.0 * (1.0 + NU));
                            break;

                        case 1:
                            // Plane stress
                            C(0, 0) = C(1, 1) = (E / (1.0 - pow(NU, 2.0)));
                            C(0, 1) = C(1, 0) = (E * NU / (1.0 - pow(NU, 2.0)));
                            C(2, 2) = E / (2.0 * (1.0 + NU));
                            break;
                        }
                    }

                    // In 2D, we have 3 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf) = gshape(spf, 1);
                        B(2, spf + dof) = gshape(spf, 0);
                    }

                    // volumetric stress. Let stiffness C act on identity and multiply by strain coeff.
                    C.Mult(identity, volumetric_stress);
                    volumetric_stress *= vol_strain_coeff;

                    switch (planeApprox)
                    {
                    case 1:
                        // plane stress
                        volumetric_stress *= (2.0 - (1.0 / (1.0 - NU)));
                        break;
                    }

                    break;
                case 3:
                    if (i == 0)
                    {
                        C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                        C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                        C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                    }

                    // In 3D, we have 6 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf + 2 * dof) = gshape(spf, 2);
                        B(3, spf + dof) = gshape(spf, 2);
                        B(3, spf + 2 * dof) = gshape(spf, 1);
                        B(4, spf) = gshape(spf, 2);
                        B(4, spf + 2 * dof) = gshape(spf, 0);
                        B(5, spf) = gshape(spf, 1);
                        B(5, spf + dof) = gshape(spf, 0);
                    }

                    // volumetric stress. Let stiffness C act on identity and multiply by strain coeff.
                    C.Mult(identity, volumetric_stress);
                    volumetric_stress *= vol_strain_coeff;

                    break;
                }
                B.MultTranspose(volumetric_stress, elvec_input);
                add(elvect, w * degradation_constant, elvec_input, elvect);
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture linear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}