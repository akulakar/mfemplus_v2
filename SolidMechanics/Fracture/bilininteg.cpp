//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Bilinear form intergrators for fracture
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "bilininteg.hpp"
#include <cmath>
#include <Dense>
#include <SVD>

using namespace std;
using namespace mfem;

namespace mfemplus
{
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------
    // Bilinear form integrators for regularized variational fracture.
    // -------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------

    // Isotropic material
    void DamagedIsotropicElasticityIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        real_t w, E, NU;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof * dim);
        elmat = 0.0;

        eldofs.SetSize(dof); // scalar for damage
        eldofdamage.SetSize(dof);
        eldofs = 0;
        eldofdamage = 0.0;

        damage_fes->GetElementDofs(elnum, eldofs);

        for (int i = 0; i < eldofdamage.Size(); i++)
        {
            eldofdamage(i) = (*damage_gf)(eldofs[i]);
        }

        // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
        // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

        const IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim); // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        C = 0.0;
        B = 0.0;
        double damage_value, degradation_constant;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            // k_eps_val not expected to change spatially.
            // multiply degradation constant to the final element matrix.

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);
            damage_value = InnerProduct(shape, eldofdamage);                                  // Should this be physical or reference shape values? Should be the same.
            degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

            w = ip.weight * Trans.Weight();                // Quadrature weights
            Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.
            }
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
                        C(0, 0) = C(1, 1) = (E / (1.0 - pow(NU, 2)));
                        C(0, 1) = C(1, 0) = (E * NU / (1.0 - pow(NU, 2)));
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

            Mult(C, B, CB);                                   // CB is 6 x (dof * dim)
            MultAtB(B, CB, elmat_input);                      // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w * degradation_constant, elmat_input); // multiplied by the degradation constant.
        }
    }

    void AnisotropicElasticityDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        real_t w;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof * dim);
        elmat = 0.0;

        eldofs.SetSize(dof); // scalar for damage
        eldofdamage.SetSize(dof);
        eldofs = 0;
        eldofdamage = 0.0;

        damage_fes->GetElementDofs(elnum, eldofs);

        for (int i = 0; i < eldofdamage.Size(); i++)
        {
            eldofdamage(i) = (*damage_gf)(eldofs[i]);
        }

        // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
        // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

        const IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim); // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        C = 0.0;
        B = 0.0;
        double damage_value, degradation_constant;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            // k_eps_val not expected to change spatially.
            // multiply degradation constant to the final element matrix.

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);
            damage_value = InnerProduct(shape, eldofdamage);                                  // Should this be physical or reference shape values? Should be the same.
            degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

            w = ip.weight * Trans.Weight();                // Quadrature weights
            Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            stiffness->Eval(C, Trans, ip); // The stiffness matrix is evaluated at each integration point.

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

            Mult(C, B, CB);                                   // CB is 6 x (dof * dim)
            MultAtB(B, CB, elmat_input);                      // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w * degradation_constant, elmat_input); // multiplied by the degradation constant.
        }
    }

    void IsotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        double w, E, NU;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        eldofs.SetSize(dof * dim); // vector valued for displacement
        eldofdisp.SetSize(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        body_pressure.SetSize(str_comp);
        elmat_input.SetSize(dof, dof);
        C = 0.0;
        B = 0.0;
        body_pressure = 0.0;

        double lambda1(0.0), lambda2(0.0), lambda3(0.0);
        double strain_energy(0.0), pressure_coeff(0.0), pressure_energy(0.0), quadratic_pressure_energy(0.0), total_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            w = ip.weight * Trans.Weight(); // Quadrature weights

            Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at the first integration point.
            }

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                // volumetric pressure. Evaluate at each integration point.
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Trans, ip);

                    // Modify body_pressure according to plane stress or plane strain using the constants.
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
                    pressure_coeff = volumetric_pressure->Eval(Trans, ip);
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

            B.Mult(eldofdisp, Bu);                                                                  // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = InnerProduct(CBu, Bu);                                                  // This is twice the strain energy
            pressure_energy = InnerProduct(body_pressure, Bu);                                      // This is the linear pressure energy
            total_energy = strain_energy - 2.0 * pressure_energy + 2.0 * quadratic_pressure_energy; // This is twice the total energy

            // // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
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
            //     total_energy = total_energy * (1e-20);
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            AddMult_a_VVt(w * total_energy, shape, elmat); // multiplied by twice the strain energy - twice the pressure energy.
        }
    }

    void AnisotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        double w;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        eldofs.SetSize(dof * dim); // vector valued for displacement
        eldofdisp.SetSize(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        elmat_input.SetSize(dof, dof);
        C = 0.0;
        B = 0.0;

        double lambda1(0.0), lambda2(0.0), lambda3(0.0);
        double strain_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            w = ip.weight * Trans.Weight(); // Quadrature weights

            Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            stiffness->Eval(C, Trans, ip); // The stiffness matrix is evaluated at each integration point.

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

            // // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
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
            //     total_energy = total_energy * (1e-20);
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            AddMult_a_VVt(w * strain_energy, shape, elmat); // multiplied by twice the strain energy - twice the pressure energy.
        }
    }

    void IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        double w, E, NU;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        eldofs.SetSize(dof * dim); // vector valued for displacement
        eldofdisp.SetSize(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * el.GetOrder() + Trans.OrderGrad(&el);
            ir = &IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        identity.SetSize(str_comp);
        volumetric_stress.SetSize(str_comp);
        elmat_input.SetSize(dof, dof);
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
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            w = ip.weight * Trans.Weight(); // Quadrature weights

            Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            vol_strain_coeff = volumetric_strain_coeff->Eval(Trans, ip);

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at the first integration point.
            }

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

                // Modify body_pressure according to plane stress or plane strain using the constants.
                switch (planeApprox)
                {
                case 1:
                    // Plane stress
                    volumetric_stress *= (2.0 - (1.0 / (1.0 - NU)));
                    // need to check this.
                    quadratic_volumetric_strain_energy = (1.0 / 2.0) * vol_strain_coeff * InnerProduct(volumetric_stress, identity) * (pow(1.0 - 2.0 * NU, 2.0) / (E * (1.0 - NU)));
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
            volumetric_strain_energy = InnerProduct(volumetric_stress, Bu);                                           // This is the linear volumetric strain energy
            total_energy = strain_energy - 2.0 * volumetric_strain_energy + 2.0 * quadratic_volumetric_strain_energy; // This is twice the total energy

            // // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
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
            //     total_energy = total_energy * (1e-20);
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = InnerProduct(CBu, Bu);
            // }

            AddMult_a_VVt(w * total_energy, shape, elmat); // multiplied by twice the strain energy - twice the pressure energy.
        }
    }

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Regularized variational fracture bilinear form integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
}