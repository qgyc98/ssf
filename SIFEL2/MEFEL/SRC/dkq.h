#ifndef DKQ_H
#define DKQ_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class dkq defines plate quadrilateral finite element based
   on the Kirchhoff plate theory for thin plates
   
   JK, 20. 10. 2019
*/

class dkq
{
 public:
  dkq (void);
  ~dkq (void);
  
  void auxdata (vector &x,vector &y,vector &l,vector &sx,vector &sy,vector &nx,vector &ny,vector &signs,ivector &nodes);
  double approx (double xi,double eta,vector &nodval);
  void bf_deflection_matrix (matrix &n,double xi,double eta);
  void bf_linear_matrix (matrix &n,double xi,double eta);

  void geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y,vector &sx,vector &sy,vector &nx,vector &ny,vector &l,vector &signs,double &jac);
  void shear_geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y,vector &sx,vector &sy,vector &nx,vector &ny,vector &l,vector &signs);

  void transf_matrix (ivector &nodes,matrix &tmat);

  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);

  void surfload (long eid, double *nodvals, vector &x, vector &y, vector &nf);  
  void node_forces_surf (long eid,double *nodvals,vector &nf);
  void nodeforces (long eid,long *le,double *nv,vector &x,vector &y,vector &nf);
  void node_forces_edge (long eid,long *le,double *nv,vector &nf);

  void ip_curvatures (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_ip_curvatures (long lcid,long eid);
  void moments (long lcid,long eid,long ri,long ci);
  void res_moments (long lcid,long eid);
  void forces (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_forces (long lcid,long eid);

  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);
  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);  
  void res_internal_forces (long lcid,long eid,vector &ifor);
  
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  number of strain/stress components connected with bending
  long bncomp;
  ///  number of strain/stress components connected with shear
  long sncomp;
  ///  number of components for graphic purposes
  long gncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on surface
  long nnsurf;
  ///  array containing orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
  long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;

};

#endif
