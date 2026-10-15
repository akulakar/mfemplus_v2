//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Strain energy evaluators for the damage (phase field) equation.
//
// The damage matrix and right hand side both integrate twice the strain energy density 2W against the damage shape
// functions. An evaluator computes 2W at the quadrature points of an element. It is used in two ways:
//   - internally by the damage integrators (the constructors that take E, nu and the displacement), which evaluate 2W
//     on the fly, and
//   - to fill a QuadratureFunction once per staggered iteration (Evaluate), which the matrix and the right hand side
//     integrators then share (the constructors that take a QuadratureFunction).
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_FRACTURE_ENERGY
#define MFEMPLUS_FRACTURE_ENERGY

#include "mfem.hpp"
#include <memory>

namespace mfemplus
{
    /** Base class of the isotropic strain energy evaluators. Computes twice the strain energy density, 2W, of the
       displacement `disp_gf` at the quadrature points of an element. E and nu are evaluated at the first integration
       point of each element (constant per element). In 2D, `plane_approximation` selects plane strain (0) or plane
       stress (1). Not thread-safe (scratch members).
    **/
    class StrainEnergyEvaluator
    {
    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::GridFunction *disp_gf;
        mfem::FiniteElementSpace *disp_fes;
        int planeApprox;
        mfem::Array<int> eldofs;
        mfem::Vector eldofdisp;
        mfem::Vector CBu, Bu; // stress and strain in Voigt form
        mfem::DenseMatrix dshape, gshape;

        /// Gathers the element displacement and sizes the scratch members.
        void PrepareElement(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr);

    public:
        StrainEnergyEvaluator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation)
            : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation) {}
        virtual ~StrainEnergyEvaluator() = default;

        /// Computes 2W at the points of `ir` for the element with transformation `Tr`. `energy` must have size
        /// ir.GetNPoints().
        virtual void ElementEnergy(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::IntegrationRule &ir, mfem::Vector &energy) = 0;

        /// Fills `W` with 2W at the quadrature points of every element of its QuadratureSpace (vdim = 1). Call it
        /// again whenever the displacement changes, before assembling the integrators that read `W`.
        void Evaluate(mfem::QuadratureFunction &W);
    };

    /** 2W = C_{ijkl} u_{k,l} u_{i,j} for an isotropic material. **/
    class IsotropicStrainEnergyEvaluator : public StrainEnergyEvaluator
    {
    public:
        IsotropicStrainEnergyEvaluator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : StrainEnergyEvaluator(e, nu, disp_gridfunc, disp_fespace, plane_approximation) {}

        void ElementEnergy(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::IntegrationRule &ir, mfem::Vector &energy) override;
    };

    /** 2W = C_{ijkl} u_{k,l} u_{i,j} - 2 p tr(eps) + (quadratic term in p), with the volumetric pressure p evaluated at
       each integration point. Only the 3D and the 2D plane stress cases carry the pressure terms.
    **/
    class IsotropicStrainEnergyPressureEvaluator : public StrainEnergyEvaluator
    {
    protected:
        mfem::Coefficient *volumetric_pressure;
        mfem::Vector body_pressure;

    public:
        IsotropicStrainEnergyPressureEvaluator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_pressure, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0)
            : StrainEnergyEvaluator(e, nu, disp_gridfunc, disp_fespace, plane_approximation), volumetric_pressure(&vol_pressure) {}

        void ElementEnergy(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::IntegrationRule &ir, mfem::Vector &energy) override;
    };

    /** 2W = (eps - alpha I) : C : (eps - alpha I) with the volumetric eigenstrain alpha I (alpha evaluated at each
       integration point), including (optionally, on by default) the Gershgorin-type switch: the energy is scaled by 1e-20
       where the effective stress is not guaranteed to be in tension, see `detail::EffectiveStressBoundIsNegative`.
    **/
    class IsotropicStrainEnergyVolumetricEigenStrainEvaluator : public StrainEnergyEvaluator
    {
    protected:
        mfem::Coefficient *volumetric_strain_coeff;
        mfem::Vector identity, volumetric_stress; // identity and volumetric stress (pressure)
        bool gershgorin_check;

    public:
        /// `use_gershgorin_check` (default true) switches the Gershgorin-type tension check on or off. When off, 2W is
        /// always the full (eps - alpha I) : C : (eps - alpha I).
        IsotropicStrainEnergyVolumetricEigenStrainEvaluator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &vol_strain_coeff, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0, bool use_gershgorin_check = true)
            : StrainEnergyEvaluator(e, nu, disp_gridfunc, disp_fespace, plane_approximation), volumetric_strain_coeff(&vol_strain_coeff), gershgorin_check(use_gershgorin_check) {}

        void ElementEnergy(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::IntegrationRule &ir, mfem::Vector &energy) override;
    };

    /** Where a damage integrator takes 2W from: either an evaluator that computes it on the fly with the integrator's
       own integration rule, or a shared QuadratureFunction (not owned, must outlive the integrator) filled by
       `StrainEnergyEvaluator::Evaluate`. In the second case the QuadratureSpace's rule is used for the integration, so
       the matrix and the right hand side built from the same QuadratureFunction always agree.
    **/
    class StrainEnergySource
    {
    private:
        std::unique_ptr<StrainEnergyEvaluator> evaluator;
        const mfem::QuadratureFunction *shared_energy = nullptr;

    public:
        explicit StrainEnergySource(std::unique_ptr<StrainEnergyEvaluator> eval) : evaluator(std::move(eval)) {}
        explicit StrainEnergySource(const mfem::QuadratureFunction &twice_energy) : shared_energy(&twice_energy) {}

        /// Returns the integration rule to use for the element and fills `energy` with 2W at its points. `ir` is the
        /// rule the integrator would use on its own; it is ignored when a QuadratureFunction is shared.
        const mfem::IntegrationRule &GetEnergyAndRule(const mfem::IntegrationRule &ir, const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &energy)
        {
            if (shared_energy != nullptr)
            {
                const mfem::IntegrationRule &rule = shared_energy->GetIntRule(Tr.ElementNo);
                shared_energy->GetValues(Tr.ElementNo, energy);
                return rule;
            }
            energy.SetSize(ir.GetNPoints());
            evaluator->ElementEnergy(el, Tr, ir, energy);
            return ir;
        }
    };
}

#endif
