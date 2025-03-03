#ifndef HEXINTERFACE_H
#define HEXINTERFACE_H

#include "alias.h"

struct matrix;
struct vector;
struct ivector;

/**
   Class hexinterface is hexahedral element for simulation of interface between surfaces of 
   linear 3D elements.
   
   basic data
   nne = 8 - number nodes on element
   ndofe = 24 - number of DOFs on element
   
   JK, 8. 7. 2023
*/
class hexinterface
{
 public:
  hexinterface(void);
  ~hexinterface(void);
  
  /// function approximates function defined by the nodal values
  double approx(double xi, double eta, vector &nodval);

  /// function computes global coordinates of the given integration point on element
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);

  /// function computes integration point values from the given nodal values of selected quantity
  void intpointval(long eid, vector &nodval, vector &ipval);

  /// function computes transformation %matrix from the element local coordinate system to the global one
  void coord_transf_matrix(matrix &tran, vector &gx, vector &gy, vector &gz);

  /// function computes transformation %matrix of the element from local coordinate systems at nodes to the global one
  void transf_matrix(ivector &nodes, matrix &tmat);

  /// function transforms node coordinates to the local element coordinate system
  void local_coordinates (vector &gx, vector &gy, vector &gz, vector &lx, vector &ly, double zero);

  /// function computes strain-displacement %matrix in the given point
  void geom_matrix(matrix &gm, vector &gx, vector &gy, vector &gz, double xi, double eta, double &jac);

  /// function computes the stiffness %matrix of the element
  void stiffness_matrix(long eid, matrix &sm);

  /// function computes the stiffness %matrix of the element with transformations to the nodal coordinate systems
  void res_stiffness_matrix(long eid, matrix &sm);  

  /// function computes resulting strains in the main integration points
  void res_mainip_strains(long lcid, long eid);

  /// function computes strains in the main integration points
  void mainip_strains(long lcid, long eid, vector &r);

  /// function computes resulting stresses in the main integration points
  void res_mainip_stresses(long lcid, long eid);
  void nod_strains_ip (long lcid, long eid);
  void nod_strains_comp (long lcid,long eid,double **stra);
  void nod_stresses_ip (long lcid,long eid);

  /// function computes stresses at all integration points of the element
  void compute_nlstress(long lcid, long eid, long ri, long ci);

  /// function computes stresses for nonlocal material models at all integration points of the element
  void compute_nonloc_nlstress(long lcid, long eid, long ri, long ci);

  /// function computes stress resultants 
  void internal_forces(long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y, vector &z);

  /// function computes stress resultants in the case of nonlocal material models
  void nonloc_internal_forces(long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y, vector &z);

  /// function computes stress resultants transformed to the nodal coordinate systems
  void res_internal_forces(long lcid, long eid, vector &ifor);

  /// function computes stress resultants transformed to the nodal coordinate systems in the case of nonlocal material models
  void res_nonloc_internal_forces(long lcid, long eid, vector &ifor);

  /// function integrates required quantity over the element volume, i.e. \int B^T q dV
  void elem_integration(integratedquant iq, long lcid, long eid, long ri, long ci,
                        vector &nv, vector &x, vector &y, vector &z);
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  number of approximated functions on the element
    //long napfun;
  ///  number of edges
    //long ned;
  ///  number of nodes on one edge
    //long nned;
  ///  array containing orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
    //long intordmm;
  ///  order of integration on edges
    //long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
  // array with the order/names of the strain components
  static const mechquant stra_comp_ord[];
};


#endif
