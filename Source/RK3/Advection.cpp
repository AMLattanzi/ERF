// Generated by Pankaj Jha on 05/31/2021
#include <RK3.H>

using namespace amrex;

Real
ComputeAdvectedQuantityForMom(const int &i, const int &j, const int &k,
                              const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                              const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                              const enum NextOrPrev &nextOrPrev,
                              const enum AdvectedQuantity &advectedQuantity,
                              const enum AdvectingQuantity &advectingQuantity,
                              const int &spatial_order) {
  Real advectingQty = 0.0;
  Real advectedQty = 1.0;
  if (nextOrPrev == NextOrPrev::next) {
    switch(advectedQuantity) {
    case AdvectedQuantity::u: //x-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::x, spatial_order); // u(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i, j, k)); // Effectively rho_u (i+1/2, j, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::y, spatial_order); // u(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i-1, j+1, k)); // Effectively rho_v (i-1/2, j+1, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::z, spatial_order); // u(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i-1, j, k+1)); // Effectively rho_w (i-1/2, j, k+1)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::v: //y-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::x, spatial_order); // v(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k)); // Effectively rho_u (i+1, j-1/2, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::y, spatial_order); // v(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i, j, k)); // Effectively rho_v (i, j+1/2, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::z, spatial_order); // v(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1)); // Effectively rho_w (i, j-1/2, k+1)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::w: //z-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::x, spatial_order); // w(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)); // Effectively rho_u (i+1, j, k-1/2)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::y, spatial_order); // w(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)); // Effectively rho_v (i, j+1, k-1/2)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::z, spatial_order); // w(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i, j, k)); // Effectively rho_w (i, j, k+1/2)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Advected quantity is unrecognized");
    }
  }
  else { // nextOrPrev == NextOrPrev::prev
    switch(advectedQuantity) {
    case AdvectedQuantity::u: //x-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::x, spatial_order); // u(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i-1, j, k) + rho_u(i, j, k)); // Effectively rho_u (i-1/2, j, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::y, spatial_order); // u(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i-1, j, k)); // Effectively rho_v (i-1/2, j, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::z, spatial_order); // u(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i-1, j, k)); // Effectively rho_w (i-1/2, j, k)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::v: //y-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::x, spatial_order); // v(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i, j, k) + rho_u(i, j-1, k)); // Effectively rho_u (i, j-1/2, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::y, spatial_order); // v(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k)); // Effectively rho_v (i, j-1/2, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::z, spatial_order); // v(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k)); // Effectively rho_w (i, j-1/2, k)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::w: //z-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::x, spatial_order); // w(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i, j, k) + rho_u(i, j, k-1)); // Effectively rho_u (i, j, k-1/2)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::y, spatial_order); // w(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i, j, k-1)); // Effectively rho_v (i, j, k-1/2)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::z, spatial_order); // w(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i, j, k-1)); // Effectively rho_w (i, j, k-1/2)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Advected quantity is unrecognized");
    }
  }

  return advectingQty * advectedQty;
}

Real
AdvectionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const amrex::Geometry &geom,
                            const SolverChoice &solverChoice) {

    const GpuArray<Real, AMREX_SPACEDIM> cellSize = geom.CellSizeArray();
    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real advectionContribution;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real centFluxXXNext, centFluxXXPrev, edgeFluxXYNext, edgeFluxXYPrev, edgeFluxXZNext, edgeFluxXZPrev;
            centFluxXXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            centFluxXXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxXYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxXYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxXZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            edgeFluxXZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            advectionContribution = (centFluxXXNext - centFluxXXPrev) / dx   // Contribution to x-mom eqn from advective flux in x-dir
                                  + (edgeFluxXYNext - edgeFluxXYPrev) / dy   // Contribution to x-mom eqn from advective flux in y-dir
                                  + (edgeFluxXZNext - edgeFluxXZPrev) / dz;  // Contribution to x-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::y:
            Real centFluxYYNext, centFluxYYPrev, edgeFluxYXNext, edgeFluxYXPrev, edgeFluxYZNext, edgeFluxYZPrev;
            edgeFluxYXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxYXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            centFluxYYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            centFluxYYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxYZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            edgeFluxYZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) / dx   // Contribution to y-mom eqn from advective flux in x-dir
                                  + (centFluxYYNext - centFluxYYPrev) / dy   // Contribution to y-mom eqn from advective flux in y-dir
                                  + (edgeFluxYZNext - edgeFluxYZPrev) / dz;  // Contribution to y-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::z:
            Real centFluxZZNext, centFluxZZPrev, edgeFluxZXNext, edgeFluxZXPrev, edgeFluxZYNext, edgeFluxZYPrev;
            edgeFluxZXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxZXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxZYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxZYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            centFluxZZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            centFluxZZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) / dx  // Contribution to z-mom eqn from advective flux in x-dir
                                  + (edgeFluxZYNext - edgeFluxZYPrev) / dy  // Contribution to z-mom eqn from advective flux in y-dir
                                  + (centFluxZZNext - centFluxZZPrev) / dz; // Contribution to z-mom eqn from advective flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return advectionContribution;
}

Real
ComputeAdvectedQuantityForState(const int &i, const int &j, const int &k,
                                const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                                const Array4<Real>& cell_data,
                                const enum NextOrPrev &nextOrPrev,
                                const enum AdvectedQuantity &advectedQuantity,
                                const enum AdvectingQuantity &advectingQuantity,
                                const int &spatial_order) {
  Real advectingQty = 0.0;
  Real advectedQty = 1.0;

  /* For the all three types of state equations (continuity, energy, scalar),
   * the advecting quantities are the same */
  if (nextOrPrev == NextOrPrev::next) {
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      advectingQty = rho_u(i+1, j, k); // rho_u (i+1, j, k)
      break;
    case AdvectingQuantity::rho_v:
      advectingQty = rho_v(i, j+1, k); // rho_v (i, j+1, k)
      break;
    case AdvectingQuantity::rho_w:
      advectingQty = rho_w(i, j, k+1); // rho_w (i, j, k+1)
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
  }
  else { // nextOrPrev == NextOrPrev::prev
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      advectingQty = rho_u(i, j, k); // rho_u (i, j, k)
      break;
    case AdvectingQuantity::rho_v:
      advectingQty = rho_v(i, j, k); // rho_v (i, j, k)
      break;
    case AdvectingQuantity::rho_w:
      advectingQty = rho_w(i, j, k); // rho_w (i, j, k)
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
  }

  // Compute advected quantity for different choice of AdvectingQuantity
  switch(advectedQuantity) {
  case AdvectedQuantity::unity: // continuity equation, reference cell is (i, j, k)
    advectedQty = 1.0; // for AdvectingQuantity = {rho_u, rho_v, rho_w}
    break;
  case AdvectedQuantity::theta:
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      // Get theta (i+1/2,    j, k    ) = theta on face (i+1, j  , k  ) for x-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i-1/2,    j, k    ) = theta on face (i,   j  , k  ) for x-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      break;
    case AdvectingQuantity::rho_v:
      // Get theta (i   , j+1/2, k    ) = theta on face (i  , j+1, k  ) for y-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i   , j-1/2, k    ) = theta on face (i  , j  , k  ) for y-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      break;
    case AdvectingQuantity::rho_w:
      // Get theta (i   , j    , k+1/2) = theta on face (i  , j  , k+1) for z-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i   , j    , k-1/2) = theta on face (i  , j  , k  ) for z-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
    break;
  case AdvectedQuantity::scalar:
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      // Get scalar (i+1/2,    j, k    ) = scalar on face (i+1, j  , k  ) for x-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i-1/2,    j, k    ) = scalar on face (i,   j  , k  ) for x-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      break;
    case AdvectingQuantity::rho_v:
      // Get scalar (i   , j+1/2, k    ) = scalar on face (i  , j+1, k  ) for y-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i   , j-1/2, k    ) = scalar on face (i  , j  , k  ) for y-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      break;
    case AdvectingQuantity::rho_w:
      // Get scalar (i   , j    , k+1/2) = scalar on face (i  , j  , k+1) for z-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i   , j    , k-1/2) = scalar on face (i  , j  , k  ) for z-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
    break;
  default:
    amrex::Abort("Error: Advected quantity is unrecognized");
  }

  // Return the product of advected and advecting quantities
  return advectingQty * advectedQty;
}

Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                              const Array4<Real>& cell_data, const int &qty_index,
                              const amrex::Geometry &geom,
                              const int &spatial_order) {

    const GpuArray<Real, AMREX_SPACEDIM> cellSize = geom.CellSizeArray();
    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real xFaceFluxNext, xFaceFluxPrev, yFaceFluxNext, yFaceFluxPrev, zFaceFluxNext, zFaceFluxPrev;

    switch(qty_index) {
        case Rho_comp: // Continuity
            xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_u, spatial_order);
            xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_u, spatial_order);
            yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_v, spatial_order);
            yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_v, spatial_order);
            zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_w, spatial_order);
            zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_w, spatial_order);
            break;
        case RhoTheta_comp: // Temperature
            xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_u, spatial_order);
            xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_u, spatial_order);
            yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_v, spatial_order);
            yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_v, spatial_order);
            zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_w, spatial_order);
            zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_w, spatial_order);
            break;
        case RhoScalar_comp: // Scalar
            xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_u, spatial_order);
            xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_u, spatial_order);
            yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_v, spatial_order);
            yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_v, spatial_order);
            zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_w, spatial_order);
            zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_w, spatial_order);
            break;
        default:
            amrex::Abort("Error: Conserved quantity index is unrecognized");
    }
    // Assemble advection contribution
    Real advectionContribution = (xFaceFluxNext - xFaceFluxPrev)/ dx  // Advective flux in x-dir
                               + (yFaceFluxNext - yFaceFluxPrev)/ dy  // Advective flux in y-dir
                               + (zFaceFluxNext - zFaceFluxPrev)/ dz; // Advective flux in z-dir

    return advectionContribution;
}