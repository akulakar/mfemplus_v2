//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Limiters to avoid spurious oscillations.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEMPLUS_COMPRESSIBLEFLOW_LIMITER
#define MFEMPLUS_COMPRESSIBLEFLOW_LIMITER

#include "mfem.hpp"

#include <memory>
#include <vector>

namespace mfemplus
{

   /**
    * @brief Zhang–Shu positivity-preserving limiter for the compressible Euler
    * equations, state u = (ρ, ρu, E), on a DG (L2) space.
    *
    * Per element, the element mean ū is kept and the smallest θ ∈ [0,1] is found
    * such that ρ ≥ small_dens and p ≥ small_pres at the volume quadrature points
    * and at the quadrature points on the element faces. The solution is then
    * replaced by ū + θ(u - ū), i.e. the high-order modes are scaled by θ:
    *   1. density:  θ₁ = (ρ̄ - ε_ρ)/(ρ̄ - ρ_min) if ρ_min < ε_ρ, else 1
    *   2. pressure: on the density-limited state, the largest t ∈ [0,1] with
    *      p(ū + t d) ≥ ε_p at each point, from a quadratic in t, θ₂ = min t
    *   3. θ = θ₁θ₂
    *
    * The volume and face rules have order 2p + int_order_offset.
    *
    * @note The limiter is element-local, so no communication is done. The
    * element mean must already satisfy the floors. If it does not, no θ can
    * fix it. The element is then set to its mean (θ = 0) and counted in the
    * return value of Apply().
    * @note Requires a nodal basis that is a partition of unity (e.g. L2 spaces
    * with Gauss-Legendre, Gauss-Lobatto or positive bases) so that scaling the
    * dofs about ū is the same as scaling the modes and keeps the mean.
    * @note Positivity is guaranteed at the checked points only.
    */
   class PositivityLimiter
   {
   protected:
      const mfem::FiniteElementSpace &fes; // vector DG space, vdim = dim + 2
      const mfem::real_t gamma;            // specific heat ratio, γ
      const mfem::real_t small_dens;       // density floor
      const mfem::real_t small_pres;       // pressure floor
      const int int_order_offset;          // quadrature order = 2p + int_order_offset
      const int dim;                       // spatial dimension

   public:
      /**
       * @brief Constructor
       * @param fes vector DG space of the state, vdim = dim + 2, must outlive the limiter
       * @param gamma specific heat ratio γ
       * @param small_dens density floor
       * @param small_pres pressure floor
       * @param int_order_offset quadrature order is 2p + int_order_offset
       */
      PositivityLimiter(const mfem::FiniteElementSpace &fes, mfem::real_t gamma,
                        mfem::real_t small_dens, mfem::real_t small_pres,
                        int int_order_offset = 1);

      /**
       * @brief Limit the state in place
       * @param u state (ρ, ρu, E) in the space given to the constructor
       * @return number of local elements whose mean violates a floor (0 = healthy).
       * The count is local to the MPI rank.
       */
      int Apply(mfem::ParGridFunction &u) const;
   };

   /**
    * @brief Barth–Jespersen slope limiter on a DG (L2) space, any vdim.
    *
    * Per element, the element mean ū is kept. For each component, the bounds are the minimum and maximum of
    * ū and the means of the face neighbors. The smallest θ ∈ [0,1] is found such that every component of
    * the solution at the volume and face quadrature points stays within its bounds. All components share
    * this single θ (the minimum over components and points), and the solution is replaced by ū + θ(u - ū):
    *   φ = min(1, (u_max - ū)/(u - ū)) if u > ū,  min(1, (u_min - ū)/(u - ū)) if u < ū,  1 otherwise.
    *
    * The volume and face rules have order 2p + int_order_offset.
    *
    * @note The mesh is assumed static and conforming. The face-neighbor lists are built once in the
    * constructor. For a changed mesh a new limiter must be created.
    * @note Unlike PositivityLimiter, Apply() communicates (neighbor means across MPI ranks).
    * @note Requires a nodal basis that is a partition of unity, as PositivityLimiter does.
    * @note The limiter acts on the conserved variables and does not enforce positivity. Apply
    * PositivityLimiter afterwards.
    * @note Not thread-safe. Apply() uses cached work data.
    */

   class BarthJespersenLimiter
   {
   protected:
      const mfem::ParFiniteElementSpace &fes;               // vector DG space
      const int int_order_offset;                           // quadrature order = 2p + int_order_offset
      const int dim;                                        // spatial dimension
      const int nvar;                                       // number of components (vdim)
      std::vector<std::vector<int>> loc_nbr;                // local face neighbors of each element
      std::vector<std::vector<int>> shr_nbr;                // shared (face-nbr) neighbors of each element
      std::unique_ptr<mfem::FiniteElementCollection> mcoll; // piecewise constant space for the means
      std::unique_ptr<mfem::ParFiniteElementSpace> mfes;    // vdim = nvar
      std::unique_ptr<mfem::ParGridFunction> means;         // element means

   public:
      /**
       * @brief Constructor, builds the face-neighbor lists
       * @param fes DG space of the state, must outlive the limiter, conforming mesh
       * @param int_order_offset quadrature order is 2p + int_order_offset
       */
      BarthJespersenLimiter(const mfem::ParFiniteElementSpace &fes, int int_order_offset = 1);

      /**
       * @brief Limit the state in place
       * @param u state in the space given to the constructor
       * @return number of local elements that were limited (θ < 1), local to the MPI rank
       */
      int Apply(mfem::ParGridFunction &u);
   };

} // namespace mfemplus

#endif // MFEMPLUS_COMPRESSIBLEFLOW_LIMITER
