#ifndef TETRADISCR_H
#define TETRADISCR_H

#include "alias.h"
#include "intpoints.h"
struct matrix;
struct vector;
struct ivector;

/**
   class tetradiscr defines tetrahedral elements for lattice discrete models

   8. 2. 2021, JK
*/

class tetralattice 
{
public:
  tetralattice(void);
  ~tetralattice(void);

  ///  number of DOFs per node
  long ndof;
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on one surface
  long nnsurf;
  ///  array of orders of integration of stiffness matrix
  long** intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of numerical interation on surfaces
  long intordb;
  ///  array of numbers of integration points in sets
  long** nip;
  ///  number of blocks
  long nb;
  ///  array of numbers of components of blocks
  long* ncomp;
  ///  cumulative array of numbers of components of blocks
  long* cncomp;
  ///  stress/strain state
  strastrestate ssst;
  /// nodes belonging to edge - matrix [ned][2]
  long nbe[6][2];
  /// edge belonging to facet - vector [tnip] 
  long ebf[12];
  /// nodes belonging to surface of the element - matrix [nsurf][nnsurf] 
  long nbs[4][3];
  /// lines belonging to surface of the element (with accordance to nbe) - matrix [nsurf][nnsurf] 
  long ebs[4][3];
  /// edge-nodes belonging to surface of the element - matrix [nsurf][nnsurf] (corresponding to the oposite node in nbs matrix) 
  long enbs[4][3];
  /// surface belonging to facet - vector [tnip]
  long sbf[12];
  /// nodes belonging to facet - matrix [tnip][2] 
  long nbf[12][2];

  void create_facets(long eid, matrix& mm);
  void geom_matrix(long eid, matrix& B);

  void internal_forces(long lcid, long eid, long ri, long ci, vector& ifor);
  void elem_integration(integratedquant iq, long lcid, long eid, long ri, long ci, vector& nv);
  void elem_integration2(integratedquant iq, long lcid, long eid, long ri, long ci, vector& nv, intpoints *ip);
  void compute_nlstress(long lcid, long eid, long ri, long ci);
  void strains(long lcid, long eid, long ri, long ci);
  void ip_strains(long lcid, long eid, long ri, long ci);
  void ip_strains2(long lcid, long eid, long ri, long ci, intpoints *ip);
  void transf_matrix(ivector& nodes, matrix& tmat);

  void largest_eigenfrequency (long lcid, long eid);

private:
  inline double tetVol(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) const;

};

#endif
