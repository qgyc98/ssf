#include "tetralattice.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"
#include "mathem.h"


tetralattice::tetralattice(void) {
  long i;

  //  number of DOFs per node
  ndof = 6;
  //  number nodes on element
  nne = 4;
  //  number of DOFs on element
  ndofe = nne * ndof;
  //  number of strain/stress components
  tncomp = 3;
  //  number of functions approximated
  napfun = 3;
  //  order of numerical integration of mass matrix
  intordmm = 1;
  //  number of edges on element
  ned = 6;
  //  number of nodes on one edge
  nned = 2;
  //  number of surfaces
  nsurf = 4;
  //  number of nodes on one surface
  nnsurf = 3;
  //  order of numerical integration on element edges (boundaries)
  intordb = 3;
  //  strain/stress state
  ssst = spacestress;

  //  number of blocks (parts of geometric matrix)
  nb = 1;

  //  number of strain/stress components
  ncomp = new long[nb];
  ncomp[0] = 3;

  //  cumulative number of components approximated
  cncomp = new long[nb];
  cncomp[0] = 0;

  //  number of integration points = 12 (fixed)
  tnip = 12;
  nip = new long* [nb];
  for (i = 0; i < nb; i++) {
    nip[i] = new long[nb];
  }
  nip[0][0] = tnip;


  // nodes belonging to edge
  nbe[0][0] = 0; nbe[0][1] = 1;
  nbe[1][0] = 0; nbe[1][1] = 2;
  nbe[2][0] = 1; nbe[2][1] = 2;
  nbe[3][0] = 1; nbe[3][1] = 3;
  nbe[4][0] = 2; nbe[4][1] = 3;
  nbe[5][0] = 0; nbe[5][1] = 3;

  // edge belonging to facet
  ebf[0] = 0;
  ebf[1] = 0;
  ebf[2] = 1;
  ebf[3] = 1;
  ebf[4] = 2;
  ebf[5] = 2;
  ebf[6] = 3;
  ebf[7] = 3;
  ebf[8] = 4;
  ebf[9] = 4;
  ebf[10] = 5;
  ebf[11] = 5;

  // nodes beloging to surface - surfaces indexed according to the oposite node
  nbs[0][0] = 1; nbs[0][1] = 2; nbs[0][2] = 3; // surface oposite node 0
  nbs[1][0] = 0; nbs[1][1] = 3; nbs[1][2] = 2; // surface oposite node 1
  nbs[2][0] = 0; nbs[2][1] = 1; nbs[2][2] = 3; // surface oposite node 2
  nbs[3][0] = 0; nbs[3][1] = 2; nbs[3][2] = 1; // surface oposite node 3

  // edges beloging to surface (with accordance to nbe)
  ebs[0][0] = 2; ebs[0][1] = 4; ebs[0][2] = 3;
  ebs[1][0] = 5; ebs[1][1] = 4; ebs[1][2] = 1;
  ebs[2][0] = 0; ebs[2][1] = 3; ebs[2][2] = 5;
  ebs[3][0] = 1; ebs[3][1] = 2; ebs[3][2] = 0;

  // nodes oposite the edge-node (surfacewise) 
  enbs[0][0] = 4; enbs[0][1] = 3; enbs[0][2] = 2;
  enbs[1][0] = 4; enbs[1][1] = 1; enbs[1][2] = 5;
  enbs[2][0] = 3; enbs[2][1] = 5; enbs[2][2] = 0;
  enbs[3][0] = 2; enbs[3][1] = 0; enbs[3][2] = 1;

  // surface belonging to facet
  sbf[0] = 2;
  sbf[1] = 3;
  sbf[2] = 1;
  sbf[3] = 3;
  sbf[4] = 0;
  sbf[5] = 3;
  sbf[6] = 0;
  sbf[7] = 2;
  sbf[8] = 0;
  sbf[9] = 1;
  sbf[10] = 1;
  sbf[11] = 2;

  /// nodes belonging to facet
  for (i = 0; i < tnip; i++) {
    nbf[i][0] = nbe[ebf[i]][0];
    nbf[i][1] = nbe[ebf[i]][1];
  }
}

tetralattice::~tetralattice(void) {
  long i;

  for (i = 0; i < nb; i++) {
    delete[] nip[i];
    delete[] intordsm[i];
  }
  delete[] nip;
  delete[] intordsm;

  delete[] cncomp;
  delete[] ncomp;

}

/**
   function creates facets in tetrahedral element used by lattice discrete particle model

   @param eid - element ID

   JV, 2.2021
*/
void tetralattice::create_facets(long eid, matrix& mm) {

  // lambda expresion for length
  auto length = [](double x0, double y0, double z0, double x1, double y1, double z1) -> double {
    return sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
  };

  // lambda expression for moments of inertia of tetrahedron to the first node
  auto moi = [](double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double mass, double& ixx, double& iyy, double& izz, double& iyz, double& ixz, double& ixy) {
    double dumx, dumy, dumz, cx, cy, cz, dx, dy, dz;
    dumx = x0 * x0 + x0 * x1 + x1 * x1 + x0 * x2 + x1 * x2 + x2 * x2 + x0 * x3 + x1 * x3 + x2 * x3 + x3 * x3;
    dumy = y0 * y0 + y0 * y1 + y1 * y1 + y0 * y2 + y1 * y2 + y2 * y2 + y0 * y3 + y1 * y3 + y2 * y3 + y3 * y3;
    dumz = z0 * z0 + z0 * z1 + z1 * z1 + z0 * z2 + z1 * z2 + z2 * z2 + z0 * z3 + z1 * z3 + z2 * z3 + z3 * z3;

    cx = 0.25 * (x0 + x1 + x2 + x3);
    cy = 0.25 * (y0 + y1 + y2 + y3);
    cz = 0.25 * (z0 + z1 + z2 + z3);

    dx = cx - x0;
    dy = cy - y0;
    dz = cz - z0;

    ixx = 0.1 * mass * (dumy + dumz) + mass * (dy * dy + dz * dz);
    iyy = 0.1 * mass * (dumx + dumz) + mass * (dx * dx + dz * dz);
    izz = 0.1 * mass * (dumx + dumy) + mass * (dx * dx + dy * dy);

    iyz = 0.05 * mass * (2.0 * y0 * z0 + y1 * z0 + y2 * z0 + y3 * z0 + y0 * z1 + 2.0 * y1 * z1 + y2 * z1 + y3 * z1 +
      y0 * z2 + y1 * z2 + 2.0 * y2 * z2 + y3 * z2 + y0 * z3 + y1 * z3 + y2 * z3 + 2.0 * y3 * z3) +
      mass * dy * dz;

    ixz = 0.05 * mass * (2.0 * x0 * z0 + x1 * z0 + x2 * z0 + x3 * z0 + x0 * z1 + 2.0 * x1 * z1 + x2 * z1 + x3 * z1 +
      x0 * z2 + x1 * z2 + 2.0 * x2 * z2 + x3 * z2 + x0 * z3 + x1 * z3 + x2 * z3 + 2.0 * x3 * z3) +
      mass * dx * dz;
    ixy = 0.05 * mass * (2.0 * x0 * y0 + x1 * y0 + x2 * y0 + x3 * y0 + x0 * y1 + 2.0 * x1 * y1 + x2 * y1 + x3 * y1 +
      x0 * y2 + x1 * y2 + 2.0 * x2 * y2 + x3 * y2 + x0 * y3 + x1 * y3 + x2 * y3 + 2.0 * x3 * y3) +
      mass * dx * dy;
  };


  // according to: Lattice Discrete Particle Model (LDPM) for failure behavior of concrete. I: Theory, Gianluca Cusatis, Daniele Pelessone, Andrea Mencarelli
  long i, j, ipp, n0, n1, ne, ns, dumn;
  double dum1, dum2, dum3;
  double fnx, fny, fnz, flx, fly, flz, fmx, fmy, fmz;
  double vol, density, mass, ixx, iyy, izz, iyz, ixz, ixy;
  bool orderNCES;

  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)); // nodes coordinates
  vector r(ASTCKVEC(nne)), ele(ASTCKVEC(ned)); // node radius, element edge length

  // new nodes
  vector enx(ASTCKVEC(ned)), eny(ASTCKVEC(ned)), enz(ASTCKVEC(ned)); // edge nodes
  vector snx(ASTCKVEC(ned)), sny(ASTCKVEC(ned)), snz(ASTCKVEC(ned)); // surface nodes
  double tx, ty, tz; // element center node

  vector fx(ASTCKVEC(tnip)), fy(ASTCKVEC(tnip)), fz(ASTCKVEC(tnip)); // facet centroids
  vector fa(ASTCKVEC(tnip)); // projected facet area
  ivector nodes(ASTCKIVEC(nne)); // element nodes IDs
  vector dens(ASTCKVEC(nne));  //  nodal densities

  //matrix mm(ASTCKMAT(ndofe, ndofe));

  //  cleaning of mass matrix
  nullm(mm);


  // node IDs
  Mt->give_elemnodes(eid, nodes);
  //  nodal densities
  Mc->give_density(eid, nodes, dens);
  // node coordinates
  Mt->give_node_coord3d(x, y, z, eid);

  //  all node densities are the same
  density = dens[0];

  // particle (node) radii
  for (i = 0; i < nne; i++) {
    r[i] = Mt->node_radius[nodes[i]];
  }

  // edge mid-nodes - E12, E13, E14, E23, E24, E34
  for (i = 0; i < ned; i++) {
    n0 = nbe[i][0];
    n1 = nbe[i][1];

    // edge lengths
    ele[i] = length(x[n0], y[n0], z[n0], x[n1], y[n1], z[n1]);

    dum1 = ele[i] - r[n0] - r[n1]; // gap between particles
    if (dum1 < 0.0) {
      print_err("Wrong 3D tetrahedral lattice element number %ld, negative gap! gap = %e", __FILE__, __LINE__, __func__, eid + 1, dum1);
      abort();
    }
    dum2 = (r[n0] + 0.5 * dum1) / ele[i];

    enx[i] = x[n0] + dum2 * (x[n1] - x[n0]);
    eny[i] = y[n0] + dum2 * (y[n1] - y[n0]);
    enz[i] = z[n0] + dum2 * (z[n1] - z[n0]);
  }

  // element surface nodes - F1, F2, F3, F4 
  for (i = 0; i < nsurf; i++) {
    snx[i] = sny[i] = snz[i] = 0.0;
    for (j = 0; j < nnsurf; j++) {
      n0 = nbs[i][j];
      n1 = enbs[i][j];
      dum3 = length(x[n0], y[n0], z[n0], enx[n1], eny[n1], enz[n1]);

      dum1 = dum3 - r[n0]; // gap between particle and edge
      if (dum1 < 0.0) {
        print_err("Wrong 3D tetrahedral lattice element number %ld, negative gap! gap = %e", __FILE__, __LINE__, __func__, eid + 1, dum1);
        abort();
      }
      dum2 = (r[n0] + 0.5 * dum1) / dum3;

      snx[i] += x[n0] + dum2 * (enx[n1] - x[n0]);
      sny[i] += y[n0] + dum2 * (eny[n1] - y[n0]);
      snz[i] += z[n0] + dum2 * (enz[n1] - z[n0]);
    }
    snx[i] /= nnsurf;
    sny[i] /= nnsurf;
    snz[i] /= nnsurf;
  }


  // element centroid - t
  tx = ty = tz = 0.0;
  for (i = 0; i < nsurf; i++) {
    // surfaces indexed according to the oposite node
    n0 = n1 = i;
    dum3 = length(x[n0], y[n0], z[n0], snx[n1], sny[n1], snz[n1]);

    dum1 = dum3 - r[n0]; // gap between particle and edge
    if (dum1 < 0.0) {
      print_err("Wrong 3D tetrahedral lattice element number %ld, negative gap! gap = %e", __FILE__, __LINE__, __func__, eid + 1, dum1);
      abort();
    }
    dum2 = (r[n0] + 0.5 * dum1) / dum3;

    tx += x[n0] + dum2 * (snx[n1] - x[n0]);
    ty += y[n0] + dum2 * (sny[n1] - y[n0]);
    tz += z[n0] + dum2 * (snz[n1] - z[n0]);
  }
  tx /= nsurf;
  ty /= nsurf;
  tz /= nsurf;


  // facets - centroid + area - 12x
  nullm(mm);
  ipp = Mt->elements[eid].ipp[0][0];
  for (i = 0; i < tnip; i++) {
    ne = ebf[i]; // edge node
    ns = sbf[i]; // surface node

    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node

    // local coordinate system of the projected facet - alligned with the edge, order of axes N-L-M
    fnx = (x[n1] - x[n0]) / ele[ne];
    fny = (y[n1] - y[n0]) / ele[ne];
    fnz = (z[n1] - z[n0]) / ele[ne];

    // facet centroid
    fx[i] = (enx[ne] + snx[ns] + tx) / 3.;
    fy[i] = (eny[ne] + sny[ns] + ty) / 3.;
    fz[i] = (enz[ne] + snz[ns] + tz) / 3.;

    // fl,fm - vectors in the facet plane
    flx = enx[ne] - tx;
    fly = eny[ne] - ty;
    flz = enz[ne] - tz;

    fmx = snx[ns] - tx;
    fmy = sny[ns] - ty;
    fmz = snz[ns] - tz;

    // cross product - fl x fm, vector length = twice the facet area
    dum1 = fly * fmz - flz * fmy;
    dum2 = flz * fmx - flx * fmz;
    dum3 = flx * fmy - fly * fmx;

    // projected facet area into the direction of perpedincular to the edge
    // scalar product of the facet area vector and edge unit vector
    fa[i] = fabs(0.5 * (dum1 * fnx + dum2 * fny + dum3 * fnz));

    // local coordinate system of the projected facet - remaining axis, order of axes N-L-M
    // fl - vector perpedincular to fn by means of already defined vector fm (we could use any random vector instead of fm)
    flx = fmy * fnz - fmz * fny;
    fly = fmz * fnx - fmx * fnz;
    flz = fmx * fny - fmy * fnx;
    dum1 = sqrt(flx * flx + fly * fly + flz * flz);
    flx /= dum1;
    fly /= dum1;
    flz /= dum1;

    // remaining axis m => fm vector perpedincular to fn and fl vectors
    fmx = fny * flz - fnz * fly;
    fmy = fnz * flx - fnx * flz;
    fmz = fnx * fly - fny * flx;

    // volume and moments of inertia
    // tetrahedron nodeOfOriginalElement-centroid(t)-edgeNode(en)-surfaceNode(sn)
    for (j = 0; j < 2; j++) {
      dumn = n0;
      if (j == 1) dumn = n1;

      vol = tetVol(x[dumn], y[dumn], z[dumn], tx, ty, tz, enx[ne], eny[ne], enz[ne], snx[ns], sny[ns], snz[ns]);
      orderNCES = true;
      if (vol < 0.0) {
        // tetrahedron nodeOfOriginalElement-edgeNode(en)-centroid(t)-surfaceNode(sn)
        vol = tetVol(x[dumn], y[dumn], z[dumn], enx[ne], eny[ne], enz[ne], tx, ty, tz, snx[ns], sny[ns], snz[ns]);
        orderNCES = false;
      }
      if (vol <= Mp->zero) {
        print_err("Wrong 3D tetrahedral lattice element number %ld, zero volume! vol = %e", __FILE__, __LINE__, __func__, eid + 1, vol);
        abort();
      }
      // @todo get density
      mass = vol * density;

      if (orderNCES) {
        moi(x[dumn], y[dumn], z[dumn], tx, ty, tz, enx[ne], eny[ne], enz[ne], snx[ns], sny[ns], snz[ns], mass, ixx, iyy, izz, iyz, ixz, ixy);
      } else {
        moi(x[dumn], y[dumn], z[dumn], enx[ne], eny[ne], enz[ne], tx, ty, tz, snx[ns], sny[ns], snz[ns], mass, ixx, iyy, izz, iyz, ixz, ixy);
      }

      // mass
      mm[dumn * 6 + 0][dumn * 6 + 0] += mass;
      mm[dumn * 6 + 1][dumn * 6 + 1] += mass;
      mm[dumn * 6 + 2][dumn * 6 + 2] += mass;

      // inertia
      mm[dumn * 6 + 3][dumn * 6 + 3] += ixx;
      mm[dumn * 6 + 4][dumn * 6 + 4] += iyy;
      mm[dumn * 6 + 5][dumn * 6 + 5] += izz;

      mm[dumn * 6 + 3][dumn * 6 + 4] -= ixy;
      mm[dumn * 6 + 3][dumn * 6 + 5] -= ixz;
      mm[dumn * 6 + 4][dumn * 6 + 5] -= iyz;

      mm[dumn * 6 + 4][dumn * 6 + 3] -= ixy;
      mm[dumn * 6 + 5][dumn * 6 + 3] -= ixz;
      mm[dumn * 6 + 5][dumn * 6 + 4] -= iyz;

      // additonal momentum
      // tetrahedron centroid
      dum1 = 0.25 * (x[dumn] + tx + enx[ne] + snx[ns]);
      dum2 = 0.25 * (y[dumn] + ty + eny[ne] + sny[ns]);
      dum3 = 0.25 * (z[dumn] + tz + enz[ne] + snz[ns]);

      mm[dumn * 6 + 0][dumn * 6 + 4] += mass * (z[dumn] - dum3);
      mm[dumn * 6 + 0][dumn * 6 + 5] -= mass * (y[dumn] - dum2);
      mm[dumn * 6 + 1][dumn * 6 + 3] -= mass * (z[dumn] - dum3);
      mm[dumn * 6 + 1][dumn * 6 + 5] += mass * (x[dumn] - dum1);
      mm[dumn * 6 + 2][dumn * 6 + 3] += mass * (y[dumn] - dum2);
      mm[dumn * 6 + 2][dumn * 6 + 4] -= mass * (x[dumn] - dum1);

      mm[dumn * 6 + 4][dumn * 6 + 0] += mass * (z[dumn] - dum3);
      mm[dumn * 6 + 5][dumn * 6 + 0] -= mass * (y[dumn] - dum2);
      mm[dumn * 6 + 3][dumn * 6 + 1] -= mass * (z[dumn] - dum3);
      mm[dumn * 6 + 5][dumn * 6 + 1] += mass * (x[dumn] - dum1);
      mm[dumn * 6 + 3][dumn * 6 + 2] += mass * (y[dumn] - dum2);
      mm[dumn * 6 + 4][dumn * 6 + 2] -= mass * (x[dumn] - dum1);
    }

    /*---------- DATA SAVING ----------*/

    // save values on IP
    Mm->ip[ipp].other[0] = fx[i];
    Mm->ip[ipp].other[1] = fy[i];
    Mm->ip[ipp].other[2] = fz[i];

    // local coordinates
    Mm->ip[ipp].other[3] = fnx;
    Mm->ip[ipp].other[4] = fny;
    Mm->ip[ipp].other[5] = fnz;

    Mm->ip[ipp].other[6] = flx;
    Mm->ip[ipp].other[7] = fly;
    Mm->ip[ipp].other[8] = flz;

    Mm->ip[ipp].other[9] = fmx;
    Mm->ip[ipp].other[10] = fmy;
    Mm->ip[ipp].other[11] = fmz;

    // facet area
    Mm->ip[ipp].other[12] = fa[i];

    // edge length
    Mm->ip[ipp].other[13] = ele[ebf[i]];

    // volumetric strain
    Mm->ip[ipp].other[14] = 0.0;

    ipp++;
  }


  // diagonalization of the mass matrix
  // @todo
  //diagonalization (mm);

}

/**
   function assembles strain-displacement (geometric) %matrix for all facets in tetrahedral element
   vector of strains has following ordering eps=(e_N, e_L, e_M)
   geometric %matrix size is (tncomp*tnip, 2*ndof)

   @param eid - element ID
   @param gm - geometric %matrix

   JV, 2.2021
*/
void tetralattice::geom_matrix(long eid, matrix& gm) {
  long i, ri, ipp, n0, n1;
  double dum1;
  double* a;
  double fx, fy, fz;
  double fnx, fny, fnz, flx, fly, flz, fmx, fmy, fmz;
  double A04_0, A05_0, A15_0, A04_1, A05_1, A15_1;

  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)); // nodes coordinates

  // node coordinates
  Mt->give_node_coord3d(x, y, z, eid);

  nullm(gm);
  ipp = Mt->elements[eid].ipp[0][0];
  ri = 0;
  for (i = 0; i < tnip; i++) {
    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node

    dum1 = 1. / Mm->ip[ipp].other[13]; // 1. / edge length

    fx = Mm->ip[ipp].other[0];
    fy = Mm->ip[ipp].other[1];
    fz = Mm->ip[ipp].other[2];

    fnx = Mm->ip[ipp].other[3] * dum1;
    fny = Mm->ip[ipp].other[4] * dum1;
    fnz = Mm->ip[ipp].other[5] * dum1;

    flx = Mm->ip[ipp].other[6] * dum1;
    fly = Mm->ip[ipp].other[7] * dum1;
    flz = Mm->ip[ipp].other[8] * dum1;

    fmx = Mm->ip[ipp].other[9] * dum1;
    fmy = Mm->ip[ipp].other[10] * dum1;
    fmz = Mm->ip[ipp].other[11] * dum1;

    A04_0 = (fz - z[n0]);
    A05_0 = (y[n0] - fy);
    A15_0 = (fx - x[n0]);

    A04_1 = (fz - z[n1]);
    A05_1 = (y[n1] - fy);
    A15_1 = (fx - x[n1]);

    // eps_N
    a = gm[ri++];
    a[0] = -fnx;
    a[1] = -fny;
    a[2] = -fnz;
    a[3] = fny * A04_0 + fnz * A05_0;
    a[4] = -fnx * A04_0 + fnz * A15_0;
    a[5] = -fnx * A05_0 - fny * A15_0;

    a[6] = fnx;
    a[7] = fny;
    a[8] = fnz;
    a[9] = -fny * A04_1 - fnz * A05_1;
    a[10] = fnx * A04_1 - fnz * A15_1;
    a[11] = fnx * A05_1 + fny * A15_1;

    // eps_L
    a = gm[ri++];
    a[0] = -flx;
    a[1] = -fly;
    a[2] = -flz;
    a[3] = fly * A04_0 + flz * A05_0;
    a[4] = -flx * A04_0 + flz * A15_0;
    a[5] = -flx * A05_0 - fly * A15_0;

    a[6] = flx;
    a[7] = fly;
    a[8] = flz;
    a[9] = -fly * A04_1 - flz * A05_1;
    a[10] = flx * A04_1 - flz * A15_1;
    a[11] = flx * A05_1 + fly * A15_1;

    // eps_M
    a = gm[ri++];
    a[0] = -fmx;
    a[1] = -fmy;
    a[2] = -fmz;
    a[3] = fmy * A04_0 + fmz * A05_0;
    a[4] = -fmx * A04_0 + fmz * A15_0;
    a[5] = -fmx * A05_0 - fmy * A15_0;

    a[6] = fmx;
    a[7] = fmy;
    a[8] = fmz;
    a[9] = -fmy * A04_1 - fmz * A05_1;
    a[10] = fmx * A04_1 - fmz * A15_1;
    a[11] = fmx * A05_1 + fmy * A15_1;

    ipp++;
  }
}

/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces

   JV 2.2021
*/
void tetralattice::internal_forces(long lcid, long eid, long ri, long ci, vector& ifor) {
  integratedquant iq;
  iq = locstress;

  //  computation of stresses
  compute_nlstress(lcid, eid, ri, ci); // stresses on facets

  //  integration of stresses over the element 
  elem_integration(iq, lcid, eid, ri, ci, ifor);
}

/**
   The function integrates selected stress type quantity over the finite element, i.e.
   it performs \int_{\Omega} \mbf{B}^T \mbf{\sigma} d\Omega. It results in nodal values.

   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values (output)

   @return The function returns nodal values calculated in the %vector nv

   JV 2.2021
*/
void tetralattice::elem_integration(integratedquant iq, long lcid, long eid, long /*ri*/, long /*ci*/, vector& nv) {
  long i, ii, k, l, n0, n1, ipp;
  double vol, ipv0, ipv1, ipv2;
  double* a0, * a1, * a2, * b;
  matrix gm(ASTCKMAT(tncomp * tnip, 2 * ndof)); // geoetric matrix
  vector ipv(ASTCKVEC(tncomp));                 // itegration point values

  nullv(nv);
  ipp = Mt->elements[eid].ipp[0][0];
  geom_matrix(eid, gm);
  ii = 0;
  for (i = 0; i < tnip; i++) {
    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node

    vol = Mm->ip[ipp].other[12] * Mm->ip[ipp].other[13]; // facet area * edge length

    Mm->givequantity(iq, lcid, ipp, cncomp[0], ipv); //  function assembles required quantity at integration point
    ipv0 = ipv[0];
    ipv1 = ipv[1];
    ipv2 = ipv[2];

    a0 = gm[ii++]; // first row of corresponding B matrix
    a1 = gm[ii++]; // second row of corresponding B matrix
    a2 = gm[ii++]; // third row of corresponding B matrix
    // first node
    b = nv.a + n0 * ndof;
    for (k = 0; k < ndof; k++) {
      b[k] += vol * (ipv0 * a0[k] + ipv1 * a1[k] + ipv2 * a2[k]);
    }

    // second node
    b = nv.a + n1 * ndof;
    for (k = ndof, l = 0; l < ndof; k++, l++) {
      b[l] += vol * (ipv0 * a0[k] + ipv1 * a1[k] + ipv2 * a2[k]);
    }

    ipp++;
  }
}

/**
   The function integrates selected stress type quantity over the finite element, i.e.
   it performs \int_{\Omega} \mbf{B}^T \mbf{\sigma} d\Omega. It results in nodal values.

   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values (output)

   @return The function returns nodal values calculated in the %vector nv

   JV 2.2021
*/
void tetralattice::elem_integration2(integratedquant /*iq*/, long /*lcid*/, long eid, long /*ri*/, long /*ci*/, vector& nv, intpoints *ip) {
  long i, ii, k, l, n0, n1, ipp;
  double vol, ipv0, ipv1, ipv2;
  double* a0, * a1, * a2, * b;
  matrix gm(ASTCKMAT(tncomp * tnip, 2 * ndof)); // geoetric matrix
  vector ipv(ASTCKVEC(tncomp));                 // itegration point values

  nullv(nv);
  ipp = Mt->elements[eid].ipp[0][0];
  geom_matrix(eid, gm);
  ii = 0;
  for (i = 0; i < tnip; i++) {
    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node
    
    //  toto ne
    //vol = Mm->ip[ipp].other[12] * Mm->ip[ipp].other[13]; // facet area * edge length
    //  neco takoveho
    vol = ip[i].other[12] * ip[i].other[13]; // facet area * edge length
    
    //  toto ne
    //Mm->givequantity(iq, lcid, ipp, cncomp[0], ipv); //  function assembles required quantity at integration point
    //ipv0 = ipv[0];
    //ipv1 = ipv[1];
    //ipv2 = ipv[2];
    ipv0 = ip[i].stress[0];
    ipv1 = ip[i].stress[1];
    ipv2 = ip[i].stress[2];
    

    a0 = gm[ii++]; // first row of corresponding B matrix
    a1 = gm[ii++]; // second row of corresponding B matrix
    a2 = gm[ii++]; // third row of corresponding B matrix
    // first node
    b = nv.a + n0 * ndof;
    for (k = 0; k < ndof; k++) {
      b[k] += vol * (ipv0 * a0[k] + ipv1 * a1[k] + ipv2 * a2[k]);
    }

    // second node
    b = nv.a + n1 * ndof;
    for (k = ndof, l = 0; l < ndof; k++, l++) {
      b[l] += vol * (ipv0 * a0[k] + ipv1 * a1[k] + ipv2 * a2[k]);
    }

    ipp++;
  }
}

/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices

   JV 2.2021
*/
void tetralattice::compute_nlstress(long /*lcid*/, long eid, long ri, long ci) {
  long i, ipp;

  //  computation of correct stresses
  ipp = Mt->elements[eid].ipp[ri][ci];
  if (Mp->strcomp == 1) {
    for (i = 0; i < tnip; i++) {
      Mm->computenlstresses(ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}

/**
   function for strain calculation

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   JV 2.2021
*/
void tetralattice::strains(long lcid, long eid, long ri, long ci) {

  //double** stra = NULL;
  vector coord, eps;

  if (Mp->strainaver == 0) {
    fprintf(stderr, "\n\n unknown strain averaging is required in function tetralattice::strains (%s, line %d).\n", __FILE__, __LINE__);
  }

  switch (Mm->stra.tape[eid]) {
  case nowhere: {
    break;
  }
  case intpts: {
    //  strains are computed at integration points/facets
    ip_strains(lcid, eid, ri, ci);
    break;
  }
  default: {
    fprintf(stderr, "\n\n unknown strain point is required in function tetralattice::strains (%s, line %d).\n", __FILE__, __LINE__);
  }
  }
  if (Mp->strainaver == 0) {
    // deletion of cretaed vectors
  }
}

/**
   function computes strains in all integration points/facets of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   JV 2.2021
*/
void tetralattice::ip_strains(long lcid, long eid, long ri, long ci) {
  long i, k, l, n0, n1, aca, ipp;
  double dum, epsVol;
  double* a, * b;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)); // coordinates of nodes
  vector r(ASTCKVEC(ndofe));                                   // nodal displacements and rotations
  vector eps(ASTCKVEC(tncomp));                                // facet strain vector
  ivector nodes(ASTCKIVEC(nne));                               // element nodes IDs

  vector aux;
  matrix gm(ASTCKMAT(tncomp * tnip, 2 * ndof)), tmat;

  Mt->give_elemnodes(eid, nodes);
  Mt->give_node_coord3d(x, y, z, eid);
  eldispl(lcid, eid, r.a); // element nodal displacements 24x

  //  transformation of displacement-rotation vector into local coordinate system
  long transf = Mt->locsystems(nodes);
  if (transf > 0) {
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }

  // volumetric strain
  dum = tetVol(x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2], x[3], y[3], z[3]); // initial volume
  epsVol = tetVol(x[0] + r[0], y[0] + r[1], z[0] + r[2],
                  x[1] + r[6], y[1] + r[7], z[1] + r[8],
                  x[2] + r[12], y[2] + r[13], z[2] + r[14],
                  x[3] + r[18], y[3] + r[19], z[3] + r[20]); // new volume

  epsVol = (epsVol - dum) / (dum * 3.0); // hydrostatic strain - in paper labeled as volumetric

  ipp = Mt->elements[eid].ipp[ri][ci];
  geom_matrix(eid, gm);
  for (i = 0; i < tnip; i++) {
    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node

    aca = 0;
    a = gm.a + i * tncomp * gm.n;
    for (k = 0; k < tncomp; k++) {
      dum = 0.0;

      // first node
      b = r.a + n0 * ndof;
      for (l = 0; l < ndof; l++) {
        dum += a[aca] * b[l];
        aca++;
      }

      // second node
      b = r.a + n1 * ndof;
      for (l = 0; l < ndof; l++) {
        dum += a[aca] * b[l];
        aca++;
      }
      eps[k] = dum;
    }

    // update strain vector
    Mm->storestrain(lcid, ipp, eps);

    // update volumetric strain
    Mm->ip[ipp].other[14] = epsVol;

    ipp++;
  }
}

/**
   function computes strains in all integration points/facets of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   JV 2.2021
*/
void tetralattice::ip_strains2(long lcid, long eid, long ri, long ci, intpoints *ip) {
  long i, k, l, n0, n1, aca, ipp;
  double dum, epsVol;
  double* a, * b;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)); // coordinates of nodes
  vector r(ASTCKVEC(ndofe));                                   // nodal displacements and rotations
  vector eps(ASTCKVEC(tncomp));                                // facet strain vector
  ivector nodes(ASTCKIVEC(nne));                               // element nodes IDs

  vector aux;
  matrix gm(ASTCKMAT(tncomp * tnip, 2 * ndof)), tmat;

  Mt->give_elemnodes(eid, nodes);
  Mt->give_node_coord3d(x, y, z, eid);
  eldispl(lcid, eid, r.a); // element nodal displacements 24x

  //  transformation of displacement-rotation vector into local coordinate system
  long transf = Mt->locsystems(nodes);
  if (transf > 0) {
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }

  // volumetric strain
  dum = tetVol(x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2], x[3], y[3], z[3]); // initial volume
  epsVol = tetVol(x[0] + r[0], y[0] + r[1], z[0] + r[2],
                  x[1] + r[6], y[1] + r[7], z[1] + r[8],
                  x[2] + r[12], y[2] + r[13], z[2] + r[14],
                  x[3] + r[18], y[3] + r[19], z[3] + r[20]); // new volume

  epsVol = (epsVol - dum) / (dum * 3.0); // hydrostatic strain - in paper labeled as volumetric

  ipp = Mt->elements[eid].ipp[ri][ci];
  geom_matrix(eid, gm);
  for (i = 0; i < tnip; i++) {
    n0 = nbf[i][0]; // first edge node
    n1 = nbf[i][1]; // second edge node

    aca = 0;
    a = gm.a + i * tncomp * gm.n;
    for (k = 0; k < tncomp; k++) {
      dum = 0.0;

      // first node
      b = r.a + n0 * ndof;
      for (l = 0; l < ndof; l++) {
        dum += a[aca] * b[l];
        aca++;
      }

      // second node
      b = r.a + n1 * ndof;
      for (l = 0; l < ndof; l++) {
        dum += a[aca] * b[l];
        aca++;
      }
      eps[k] = dum;
    }

    // update strain vector
    //  toto ne
    //Mm->storestrain(lcid, ipp, eps);
    //
    //  nejak takto
    for (k=0;k<tncomp;k++){
      ip[i].strain[k]=eps[k];
    }


    // update volumetric strain
    //  toto ne
    //Mm->ip[ipp].other[14] = epsVol;
    //
    //  nejak takto
    ip[i].other[14] = epsVol;
    
    ipp++;
  }
}


/**
   function assembles transformation %matrix

   @param nodes - nodes of element
   @param tmat - transformation %matrix

   JV 2.2021
*/
void tetralattice::transf_matrix(ivector& nodes, matrix& tmat) {
  long i, n, m;

  nullm(tmat);

  n = nodes.n;
  m = tmat.m;
  for (i = 0; i < m; i++) {
    tmat[i][i] = 1.0;
  }

  for (i = 0; i < n; i++) {
    if (Mt->nodes[nodes[i]].transf > 0) {
      tmat[i * 6 + 0][i * 6] = tmat[i * 6 + 3][i * 6 + 3] = Mt->nodes[nodes[i]].e1[0];
      tmat[i * 6 + 1][i * 6] = tmat[i * 6 + 4][i * 6 + 3] = Mt->nodes[nodes[i]].e1[1];
      tmat[i * 6 + 2][i * 6] = tmat[i * 6 + 5][i * 6 + 3] = Mt->nodes[nodes[i]].e1[2];

      tmat[i * 6 + 0][i * 6 + 1] = tmat[i * 6 + 3][i * 6 + 4] = Mt->nodes[nodes[i]].e2[0];
      tmat[i * 6 + 1][i * 6 + 1] = tmat[i * 6 + 4][i * 6 + 4] = Mt->nodes[nodes[i]].e2[1];
      tmat[i * 6 + 2][i * 6 + 1] = tmat[i * 6 + 5][i * 6 + 4] = Mt->nodes[nodes[i]].e2[2];

      tmat[i * 6 + 0][i * 6 + 2] = tmat[i * 6 + 3][i * 6 + 5] = Mt->nodes[nodes[i]].e3[0];
      tmat[i * 6 + 1][i * 6 + 2] = tmat[i * 6 + 4][i * 6 + 5] = Mt->nodes[nodes[i]].e3[1];
      tmat[i * 6 + 2][i * 6 + 2] = tmat[i * 6 + 5][i * 6 + 5] = Mt->nodes[nodes[i]].e3[2];
    }
  }
}

/**
   function calculates the volume of tetrahedron

   @param x0,y0,z0 - coordinates of first node
   @param x1,y1,z1 - coordinates of second node
   @param x2,y2,z2 - coordinates of third node
   @param x3,y3,z3 - coordinates of fourth node

   @return The function returns volume of the tetrahedron defined by 4 nodes

   JV 3.2021
*/
inline double tetralattice::tetVol(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) const
{
  return ((x1 - x0) * (y2 - y0) * (z3 - z0) + (y1 - y0) * (z2 - z0) * (x3 - x0) + (z1 - z0) * (x2 - x0) * (y3 - y0) -
          (z1 - z0) * (y2 - y0) * (x3 - x0) - (x1 - x0) * (z2 - z0) * (y3 - y0) - (y1 - y0) * (x2 - x0) * (z3 - z0)) / 6.0;
}


void tetralattice::largest_eigenfrequency (long lcid, long eid)
{
  long i,k,ri=0,ci=0;
  vector d(ASTCKVEC(ndofe)),nv(ASTCKVEC(ndofe));
  intpoints *ip;
  integratedquant iq;
  iq = locstress;

  ip = new intpoints [tnip];
  
  //  loop over all DOFs
  for (i=0;i<ndofe;i++){
    //  zde je treba vygenerovat spravny vektor uzlovych posunuti d
    

    //  vypocet deformaci
    //  tato funkce ma navic posledni argument
    //  ten tam v potu tvare dodam, ale zaroven to budu
    //  muset udelat na vsech prvcich, aby se to prelozilo
    ip_strains2 (lcid, eid, ri, ci, ip);
    
    
    //  vypocet napeti
    //  nova funkce
    for (k=0;k<tnip;k++)
      Mm->computenlstresses (k,ip[k]);
    
    
    //  vypocet uzlovych sil
    elem_integration2 (iq,lcid,eid,ri,ci,nv,ip);
    
    
    //  zde bude plneni matice tuhosti
    
    
  }
  
  //  zde se sestavi matice hmotnosti
  
  
  //  zde se zavola rutina na vypocet nejvetsi frekvence
  
  
}
