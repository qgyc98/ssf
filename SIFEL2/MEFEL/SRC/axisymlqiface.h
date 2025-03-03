#ifndef AXISYMLQINTERFACE_H
#define AXISYMLQINTERFACE_H

#include "alias.h"

struct matrix;
struct vector;
struct ivector;



/**
   Class axisymlqinterface is axisymmetric quadrilateral element with linear approximation of relative displacements
   for simulation of interface between edges of linear 2D axisymmetric elements.
   
   basic data
   nne = 4 - number nodes on element
   ndofe = 8 - number of DOFs on element
   
   TKo, 2.4.2024
*/
class axisymlqinterface
{
 public:
  axisymlqinterface (void);
  ~axisymlqinterface (void);

  /// function approximates function defined by the nodal values
  double approx (double xi,vector &nodval);

  /// function computes integration point values from the given nodal values of selected quantity
  void intpointval(long eid, vector &nodval, vector &ipval);

  /// assembles natural coordinates of element nodes
  void nodecoord (vector &xi);
  
  // computes global coordinates of the integration point
  void ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  
  // returns natural coordinates of the selected integration point
  void ipncoord (long eid, long ipp, vector &coord);
  
  /// function computes transformation %matrix of the element
  void transf_matrix (ivector &nod,matrix &tmat);

  /// function computes the stiffness %matrix of the element
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);

  /// function computes the stiffness %matrix of the element with transformations to the nodal coordinate systems
  void res_stiffness_matrix(long eid, matrix &sm);

  /// function computes resulting strains in the main integration points
  void res_mainip_strains (long lcid,long eid);

  /// function computes strains in the main integration points
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &r);

  /// function returns integration point numbers closest to the particular element nodes
  void nodipnum (long eid, long ri, long ci, ivector &ipnum);
  
  /// function computes averaged nodal strains from strain values at integration points
  void nod_strains_ip (long lcid,long eid,long ri,long ci);

  /// function computes averaged nodal strains
  void nod_strains_comp (long lcid,long eid);
  
  /// function computes nodal strains without averaging
  void nod_strains (long lcid,long eid);
  
  /// function computes nodal stresses from stress values at integration points
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  
  /// function computes nodal values of other array values at integration points
  void nod_other_ip (long eid);
  
  /// function computes resulting stresses in the main integration points
  void res_mainip_stresses (long lcid,long eid);

  /// function computes strain-displacement %matrix in the given point
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double &jac);

  /// function computes stresses at all integration points of the element
  void compute_nlstress (long lcid,long eid,long ri,long ci);

  /// function computes stresses for nonlocal material models at all integration points of the element
  void compute_nonloc_nlstress(long lcid, long eid, long ri, long ci);

  /// compute eigenstresses due to eigenstrains at all element integration points
  void compute_nlstressincr(long /*lcid*/, long eid, long ri, long ci);
  
  /// compute eigenstresses due to eigenstrains at all element integration points
  void compute_eigstress(long /*lcid*/, long eid, long ri, long ci);
  
  /// function computes stress resultants 
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x, vector &y);

  /// function computes stress resultants in the case of nonlocal material models
  void nonloc_internal_forces(long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y);

  /// function computes stress resultants transformed to the nodal coordinate systems
  void res_internal_forces (long lcid,long eid,vector &ifor);

  /// function computes stress resultants transformed to the nodal coordinate systems in the case of nonlocal material models
  void res_nonloc_internal_forces(long lcid, long eid, vector &ifor);

  // computes forces due to stress increments
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);

  // computes the resulting forces due to stress increments
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  
  // computes nodal forces due to eigenstrains
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y);
  
  // computes the resulting nodal forces due to eigenstrains
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);
  
  /// function integrates required quantity over the element volume, i.e. \int B^T q dV
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);


  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  /// number of edges
  long ned;
  /// number of nodes on edges
  long nned;
  /// number of approximation functions
  long napfun;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};


#endif
