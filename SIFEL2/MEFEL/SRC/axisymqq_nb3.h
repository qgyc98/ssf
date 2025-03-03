#ifndef AXISYMQQ_H
#define AXISYMQQ_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class axisymqq defines quadrilateral element for axisymmetric problem
   with quadratic approximation functions
   
   JK
*/

class axisymqq
//  je treba rozmyslet pocty integracnich bodu v nelinearnich
//  ulohach, protoze funce nlstresses neumeji pracovat
//  s nekolika integracnimi body
{
 public:
  axisymqq (void);
  ~axisymqq (void);
  void eleminit (long eid);
  double approx (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);

  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci);
  void nod_strains (long lcid,long eid,long ri,long ci);
  void nod_strains_old (long lcid,long eid,long ri,long ci);
  void elem_strains (double **stra,long lcid,long eid,long ri,long ci);
  void appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps);
  void res_allip_strains (long lcid,long eid);
  void allip_strains (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi,vector &eta);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);

  void res_mainip_stresses (long lcid,long eid);
  void mainip_stresses (long lcid,long eid,long ri,long ci,long ii);
  void nod_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_old (long lcid,long eid,long ri,long ci);
  void elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig);
  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);

  void load_matrix (long eid,matrix &lm);
  void res_temp_forces (long lcid,long eid,vector &nfor);
  void temp_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y);
  void tempstrains (long lcid,long eid,long ipp,double xi,double eta,vector &eps);
  void internal_forces (long lcid,long eid,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);


  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on edge
  long nned;
  ///  array of orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
  long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  array of numbers of components of blocks
  long *ncomp;
  ///  cumulative array of numbers of components of blocks
  long *cncomp;
  ///  stress/strain state
  strastrestate ssst;
};

#endif


