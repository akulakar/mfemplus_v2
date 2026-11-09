//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// Limiters to avoid spurious oscillations.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "mfemplus/FluidMechanics/CompressibleFlow/limiter.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;
using namespace mfem;

namespace mfemplus
{

   namespace
   {

      // Appends the quadrature points on the faces of the reference element to pts, in element
      // reference coordinates.
      void AddFacePoints(const int geom, const int qorder, vector<IntegrationPoint> &pts)
      {
         const IntegrationRule *verts = Geometries.GetVertices(geom);
         auto vtx = [&](int i) -> const IntegrationPoint &
         { return verts->IntPoint(i); };

         IntegrationPoint ip;

         if (geom == Geometry::SEGMENT)
         {
            // The faces are the two end points.
            pts.push_back(vtx(0));
            pts.push_back(vtx(1));
         }
         else if (geom == Geometry::TRIANGLE || geom == Geometry::SQUARE)
         {
            // The faces are edges, x = v0 + t (v1 - v0).
            const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, qorder);
            int ne;
            if (geom == Geometry::TRIANGLE)
            {
               ne = Geometry::Constants<Geometry::TRIANGLE>::NumEdges;
            }
            else
            {
               ne = Geometry::Constants<Geometry::SQUARE>::NumEdges;
            }
            for (int f = 0; f < ne; f++)
            {
               const int *e;
               if (geom == Geometry::TRIANGLE)
               {
                  e = Geometry::Constants<Geometry::TRIANGLE>::Edges[f];
               }
               else
               {
                  e = Geometry::Constants<Geometry::SQUARE>::Edges[f];
               }
               const IntegrationPoint &v0 = vtx(e[0]), &v1 = vtx(e[1]);
               for (int k = 0; k < ir.GetNPoints(); k++)
               {
                  const real_t t = ir.IntPoint(k).x;
                  ip.Set3(v0.x + t * (v1.x - v0.x), v0.y + t * (v1.y - v0.y), 0.0);
                  pts.push_back(ip);
               }
            }
         }
         else if (geom == Geometry::TETRAHEDRON)
         {
            // Triangular faces, x = v0 + s (v1 - v0) + t (v2 - v0).
            const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, qorder);
            for (int f = 0; f < Geometry::Constants<Geometry::TETRAHEDRON>::NumFaces; f++)
            {
               const int *fv = Geometry::Constants<Geometry::TETRAHEDRON>::FaceVert[f];
               const IntegrationPoint &v0 = vtx(fv[0]), &v1 = vtx(fv[1]), &v2 = vtx(fv[2]);
               for (int k = 0; k < ir.GetNPoints(); k++)
               {
                  const real_t s = ir.IntPoint(k).x, t = ir.IntPoint(k).y;
                  ip.Set3(v0.x + s * (v1.x - v0.x) + t * (v2.x - v0.x),
                          v0.y + s * (v1.y - v0.y) + t * (v2.y - v0.y),
                          v0.z + s * (v1.z - v0.z) + t * (v2.z - v0.z));
                  pts.push_back(ip);
               }
            }
         }
         else if (geom == Geometry::CUBE)
         {
            // Quadrilateral faces with cyclic vertices, bilinear map.
            const IntegrationRule &ir = IntRules.Get(Geometry::SQUARE, qorder);
            for (int f = 0; f < Geometry::Constants<Geometry::CUBE>::NumFaces; f++)
            {
               const int *fv = Geometry::Constants<Geometry::CUBE>::FaceVert[f];
               const IntegrationPoint &v0 = vtx(fv[0]), &v1 = vtx(fv[1]),
                                      &v2 = vtx(fv[2]), &v3 = vtx(fv[3]);
               for (int k = 0; k < ir.GetNPoints(); k++)
               {
                  const real_t s = ir.IntPoint(k).x, t = ir.IntPoint(k).y;
                  const real_t w0 = (1.0 - s) * (1.0 - t), w1 = s * (1.0 - t),
                               w2 = s * t, w3 = (1.0 - s) * t;
                  ip.Set3(w0 * v0.x + w1 * v1.x + w2 * v2.x + w3 * v3.x,
                          w0 * v0.y + w1 * v1.y + w2 * v2.y + w3 * v3.y,
                          w0 * v0.z + w1 * v1.z + w2 * v2.z + w3 * v3.z);
                  pts.push_back(ip);
               }
            }
         }
         else
         {
            MFEM_ABORT("PositivityLimiter: unsupported element geometry.");
         }
      }

      // Largest t in [0,1] such that the pressure of ū + t d is at least eps_p. Requires ρ > 0 on the
      // segment and an admissible ū. With k = eps_p/(γ-1), p ≥ eps_p is equivalent to
      //   f(t) = ρ(t) E(t) - ½|m(t)|² - k ρ(t) = a t² + b t + c ≥ 0.
      real_t PressureTheta(const real_t *ub, const real_t *d, const int dim,
                           const real_t gamma, const real_t eps_p)
      {
         const int iE = dim + 1;
         const real_t k = eps_p / (gamma - 1.0);

         real_t m2 = 0.0, bm = 0.0, d2 = 0.0; // |m̄ + d_m|², m̄·d_m, |d_m|²
         real_t mb2 = 0.0;                    // |m̄|²
         for (int i = 1; i <= dim; i++)
         {
            m2 += (ub[i] + d[i]) * (ub[i] + d[i]);
            bm += ub[i] * d[i];
            d2 += d[i] * d[i];
            mb2 += ub[i] * ub[i];
         }

         // f(1)
         const real_t rho1 = ub[0] + d[0], E1 = ub[iE] + d[iE];
         if (rho1 * E1 - 0.5 * m2 - k * rho1 >= 0.0)
         {
            return 1.0;
         }

         const real_t a = d[0] * d[iE] - 0.5 * d2;
         const real_t b = ub[0] * d[iE] + ub[iE] * d[0] - bm - k * d[0];
         const real_t c = max(ub[0] * ub[iE] - 0.5 * mb2 - k * ub[0], (real_t)0.0);

         real_t t;
         if (abs(a) <= 1e-14 * abs(b))
         {
            // Linear, f(0) ≥ 0 > f(1) implies b < 0
            t = (b < 0.0) ? -c / b : 0.0;
         }
         else
         {
            const real_t disc = b * b - 4.0 * a * c;
            if (disc < 0.0)
            {
               t = 0.0; // only from roundoff
            }
            else
            {
               const real_t q = -0.5 * (b + copysign(sqrt(disc), b));
               const real_t r1 = q / a;
               const real_t r2 = (q != 0.0) ? c / q : r1;
               t = 1.0;
               if (r1 >= 0.0 && r1 < t)
               {
                  t = r1;
               }
               if (r2 >= 0.0 && r2 < t)
               {
                  t = r2;
               }
            }
         }
         return min(max(t, (real_t)0.0), (real_t)1.0);
      }

   } // namespace

   PositivityLimiter::PositivityLimiter(const FiniteElementSpace &fes, real_t gamma,
                                        real_t small_dens, real_t small_pres,
                                        int int_order_offset)
       : fes(fes), gamma(gamma), small_dens(small_dens), small_pres(small_pres),
         int_order_offset(int_order_offset), dim(fes.GetMesh()->Dimension())
   {
      MFEM_VERIFY(fes.GetVDim() == dim + 2,
                  "PositivityLimiter requires the Euler state, vdim = dim + 2.");
   }

   int PositivityLimiter::Apply(ParGridFunction &u) const
   {
      MFEM_VERIFY(u.Size() == fes.GetVSize(),
                  "PositivityLimiter::Apply: u is not in the space given to the constructor.");

      const int nvar = dim + 2;
      const int iE = dim + 1;
      int num_bad_means = 0;

      Array<int> vdofs;
      Vector vals, shape(0), ubar(nvar), d(nvar);
      vector<IntegrationPoint> pts;
      DenseMatrix B, V;

      for (int e = 0; e < fes.GetNE(); e++)
      {
         const FiniteElement &fe = *fes.GetFE(e);
         const int ndof = fe.GetDof();
         const int geom = fe.GetGeomType();
         const int qorder = 2 * fe.GetOrder() + int_order_offset;

         // Points: volume quadrature points first, then the face quadrature points
         const IntegrationRule &ir = IntRules.Get(geom, qorder);
         const int nvol = ir.GetNPoints();
         pts.clear();
         for (int q = 0; q < nvol; q++)
         {
            pts.push_back(ir.IntPoint(q));
         }
         AddFacePoints(geom, qorder, pts);
         const int npts = (int)pts.size();

         B.SetSize(npts, ndof);
         shape.SetSize(ndof);
         for (int i = 0; i < npts; i++)
         {
            fe.CalcShape(pts[i], shape);
            if (i == 0)
            {
               MFEM_VERIFY(abs(shape.Sum() - 1.0) < 1e-8,
                           "PositivityLimiter requires a basis that is a partition of unity.");
            }
            for (int j = 0; j < ndof; j++)
            {
               B(i, j) = shape(j);
            }
         }

         fes.GetElementVDofs(e, vdofs);
         u.GetSubVector(vdofs, vals);
         DenseMatrix U(vals.GetData(), ndof, nvar);

         // State at all points
         V.SetSize(npts, nvar);
         Mult(B, U, V);

         // Element mean from the volume points
         ElementTransformation *T = fes.GetElementTransformation(e);
         ubar = 0.0;
         real_t wsum = 0.0;
         for (int q = 0; q < nvol; q++)
         {
            T->SetIntPoint(&pts[q]);
            const real_t w = pts[q].weight * T->Weight();
            wsum += w;
            for (int v = 0; v < nvar; v++)
            {
               ubar(v) += w * V(q, v);
            }
         }
         ubar /= wsum;

         // The mean must be admissible
         real_t mb2 = 0.0;
         for (int i = 1; i <= dim; i++)
         {
            mb2 += ubar(i) * ubar(i);
         }
         const real_t rho_bar = ubar(0);
         const real_t p_bar = (gamma - 1.0) * (ubar(iE) - 0.5 * mb2 / rho_bar);

         real_t theta;
         if (!(rho_bar >= small_dens) || !(p_bar >= small_pres))
         {
            theta = 0.0;
            num_bad_means++;
         }
         else
         {
            // 1. density
            real_t rho_min = V(0, 0);
            for (int q = 1; q < npts; q++)
            {
               rho_min = min(rho_min, V(q, 0));
            }
            const real_t theta1 = (rho_min < small_dens)
                                      ? (rho_bar - small_dens) / (rho_bar - rho_min)
                                      : 1.0;

            // 2. pressure on the density-limited state
            real_t theta2 = 1.0;
            for (int q = 0; q < npts; q++)
            {
               for (int v = 0; v < nvar; v++)
               {
                  d(v) = theta1 * (V(q, v) - ubar(v));
               }
               theta2 = min(theta2, PressureTheta(ubar.GetData(), d.GetData(), dim, gamma, small_pres));
            }

            theta = theta1 * theta2;
         }

         if (theta < 1.0)
         {
            for (int v = 0; v < nvar; v++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  U(i, v) = ubar(v) + theta * (U(i, v) - ubar(v));
               }
            }
            u.SetSubVector(vdofs, vals);
         }
      }

      return num_bad_means;
   }

   BarthJespersenLimiter::BarthJespersenLimiter(const ParFiniteElementSpace &fes, int int_order_offset)
       : fes(fes), int_order_offset(int_order_offset), dim(fes.GetMesh()->Dimension()),
         nvar(fes.GetVDim())
   {
      ParMesh &pmesh = *const_cast<ParMesh *>(fes.GetParMesh());
      MFEM_VERIFY(!pmesh.Nonconforming(), "BarthJespersenLimiter requires a conforming mesh.");

      const int NE = pmesh.GetNE();
      loc_nbr.assign(NE, {});
      shr_nbr.assign(NE, {});

      // Local interior faces
      for (int f = 0; f < pmesh.GetNumFaces(); f++)
      {
         int e1, e2;
         pmesh.GetFaceElements(f, &e1, &e2);
         if (e1 >= 0 && e2 >= 0)
         {
            loc_nbr[e1].push_back(e2);
            loc_nbr[e2].push_back(e1);
         }
      }

      // Shared faces: element 2 is a face-neighbor element, its index is shifted by NE
      pmesh.ExchangeFaceNbrData();
      for (int sf = 0; sf < pmesh.GetNSharedFaces(); sf++)
      {
         FaceElementTransformations *tr = pmesh.GetSharedFaceTransformations(sf);
         MFEM_VERIFY(tr->Elem1No < NE && tr->Elem2No >= NE,
                     "BarthJespersenLimiter: unexpected shared face element indices.");
         shr_nbr[tr->Elem1No].push_back(tr->Elem2No - NE);
      }

      // Piecewise constant space for the element means
      mcoll = make_unique<L2_FECollection>(0, dim);
      mfes = make_unique<ParFiniteElementSpace>(&pmesh, mcoll.get(), nvar);
      means = make_unique<ParGridFunction>(mfes.get());
   }

   int BarthJespersenLimiter::Apply(ParGridFunction &u)
   {
      MFEM_VERIFY(u.Size() == fes.GetVSize(),
                  "BarthJespersenLimiter::Apply: u is not in the space given to the constructor.");

      const int NE = fes.GetNE();
      int num_limited = 0;

      Array<int> vdofs;
      Vector vals, shape(0), ubar(nvar), umin(nvar), umax(nvar);
      vector<IntegrationPoint> pts;

      // Pass 1: element means. Order 0 L2 with byNODES ordering, so dof of (element e, component v) is e + v NE
      for (int e = 0; e < NE; e++)
      {
         const FiniteElement &fe = *fes.GetFE(e);
         const int ndof = fe.GetDof();
         const int geom = fe.GetGeomType();
         const int qorder = 2 * fe.GetOrder() + int_order_offset;
         const IntegrationRule &ir = IntRules.Get(geom, qorder);
         const int nvol = ir.GetNPoints();

         shape.SetSize(ndof);
         fes.GetElementVDofs(e, vdofs);
         u.GetSubVector(vdofs, vals);
         DenseMatrix U(vals.GetData(), ndof, nvar);

         ElementTransformation *T = fes.GetElementTransformation(e);
         ubar = 0.0;
         real_t wsum = 0.0;
         for (int q = 0; q < nvol; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            fe.CalcShape(ip, shape);
            if (q == 0)
            {
               MFEM_VERIFY(abs(shape.Sum() - 1.0) < 1e-8,
                           "BarthJespersenLimiter requires a basis that is a partition of unity.");
            }
            T->SetIntPoint(&ip);
            const real_t w = ip.weight * T->Weight();
            wsum += w;
            for (int v = 0; v < nvar; v++)
            {
               real_t s = 0.0;
               for (int j = 0; j < ndof; j++)
               {
                  s += shape(j) * U(j, v);
               }
               ubar(v) += w * s;
            }
         }
         for (int v = 0; v < nvar; v++)
         {
            (*means)(e + v * NE) = ubar(v) / wsum;
         }
      }
      means->ExchangeFaceNbrData();
      const Vector &nbr_means = means->FaceNbrData();

      // Pass 2: bounds from the neighbor means and limiting. Neighbor values are means, not u, so
      // modifying u element by element is safe.
      Array<int> nbr_vdofs;
      for (int e = 0; e < NE; e++)
      {
         const FiniteElement &fe = *fes.GetFE(e);
         const int ndof = fe.GetDof();
         const int geom = fe.GetGeomType();
         const int qorder = 2 * fe.GetOrder() + int_order_offset;

         const IntegrationRule &ir = IntRules.Get(geom, qorder);
         pts.clear();
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            pts.push_back(ir.IntPoint(q));
         }
         AddFacePoints(geom, qorder, pts);
         const int npts = (int)pts.size();

         for (int v = 0; v < nvar; v++)
         {
            ubar(v) = (*means)(e + v * NE);
            umin(v) = umax(v) = ubar(v);
            for (int n : loc_nbr[e])
            {
               const real_t m = (*means)(n + v * NE);
               umin(v) = min(umin(v), m);
               umax(v) = max(umax(v), m);
            }
            for (int n : shr_nbr[e])
            {
               mfes->GetFaceNbrElementVDofs(n, nbr_vdofs);
               const real_t m = nbr_means(nbr_vdofs[v]);
               umin(v) = min(umin(v), m);
               umax(v) = max(umax(v), m);
            }
         }

         fes.GetElementVDofs(e, vdofs);
         u.GetSubVector(vdofs, vals);
         DenseMatrix U(vals.GetData(), ndof, nvar);

         shape.SetSize(ndof);
         real_t theta = 1.0;
         for (int i = 0; i < npts && theta > 0.0; i++)
         {
            fe.CalcShape(pts[i], shape);
            for (int v = 0; v < nvar; v++)
            {
               real_t val = 0.0;
               for (int j = 0; j < ndof; j++)
               {
                  val += shape(j) * U(j, v);
               }
               const real_t d = val - ubar(v);
               const real_t tol = 1e-12 * (abs(umax(v)) + abs(umin(v))) + 1e-300;
               if (d > tol)
               {
                  theta = min(theta, (umax(v) - ubar(v)) / d);
               }
               else if (d < -tol)
               {
                  theta = min(theta, (umin(v) - ubar(v)) / d);
               }
            }
         }
         theta = max(theta, (real_t)0.0);

         if (theta < 1.0)
         {
            for (int v = 0; v < nvar; v++)
            {
               for (int j = 0; j < ndof; j++)
               {
                  U(j, v) = ubar(v) + theta * (U(j, v) - ubar(v));
               }
            }
            u.SetSubVector(vdofs, vals);
            num_limited++;
         }
      }

      return num_limited;
   }

} // namespace mfemplus
