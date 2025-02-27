#ifndef LINHEX_H
#define LINHEX_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class linhex defines hexahedral finite element with
   tri-linear approximation functions
   
   JK
*/

class linhex
{
 public:
  linhex (void);
  ~linhex (void);
  void eleminit (long eid);
  double approx (double xi,double eta,double zeta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
		    double xi,double eta,double zeta,double &jac);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci);
  void nod_strains (long lcid,long eid,long ri,long ci);
  void elem_strains (double **stra,long lcid,long eid,long ri,long ci);
  void appstrain (long lcid,long eid,double xi,double eta,double zeta,long fi,long li,vector &eps);
  void allip_strains (double **stra,long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi,vector &eta,vector &zeta);
  void appval (double xi,double eta,double zeta,long fi,long nc,vector &eps,double **val);

  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses (long lcid,long eid,long ri,long ci);
  void elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,double zeta,long fi,long li,vector &sig);
  void allip_stresses (double **stre,long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);

  void res_temp_forces (long lcid,long eid,vector &nfor);
  void temp_forces (long lcid,long eid,long ri,long ci,vector &nfor);
  void temperaturestrains (long lcid,long eid,long ri,long ci);
  void tempstrains (long lcid,long eid,long ipp,double xi,double eta,double zeta,vector &eps);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void internal_forces2 (long lcid,long eid,long ri,long ci,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void local_values (long lcid,long eid,long ri,long ci);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
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
