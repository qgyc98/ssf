#ifndef PLELEMLQ_H
#define PLELEMLQ_H

#include "alias.h"
#include "galias.h"
struct matrix;
struct vector;
struct ivector;
struct atsel;



/**
   class planeelemlq defines plane quadrilateral element with
   bi-linear approximation functions
   
   basic data
   nne = 4 - number nodes on element
   ndofe = 8 - number of DOFs on element
   ncomp = 3 - number of strain (stress) tensor components
   napfun = 2 - number of approximated functions
   intordmm = 2 - order of numerical integration of mass %matrix (2x2)
   
   nsip - number of sets of integration points
   there are two sets, 
   nip - number of integration points in sets
   intordsm - array containing orders of numerical integration
   
   JK
*/
class planeelemlq
{
 public:
  planeelemlq (void);
  ~planeelemlq (void);

  double approx (double xi,double eta,vector &nodval);
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipcoordblock (long eid,long ri,long ci,double **coord);
  void ipncoord (long eid,long ipp,vector &ncoord);

  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac);
  void bvectors (vector &x,vector &y,double xi,double eta,double &jac,
		 vector &b11,vector &b12,vector &b21,vector &b22);
  void gngeom_matrix (matrix &gm,vector &r,vector &x,vector &y,double xi,double eta,double &jac);
  void gnl_grmatrix (matrix &grm,vector &x,vector &y,double xi,double eta,double &jac);

  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nod,matrix &tmat);

  void gl_stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void gnl_stiffness_matrix (long lcid,long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long lcid,long eid,matrix &sm);
  void bd_matrix (long eid,long ri,long ci,matrix &mbd);
  void dd_matrix (long eid,long ri,long ci,matrix &mdd);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  
  void damping_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_damping_matrix (long eid,matrix &dm);

  void res_ip_strains (long lcid,long eid);
  void gl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void gnl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid);
  void nod_strains (long lcid,long eid);
  void strains (long lcid,long eid,long ri,long ci);

  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre);
  void stresses (long lcid,long eid,long ri,long ci);

  void nod_other_ip (long eid,long ri,long ci);
  
  void load_matrix (long eid,matrix &lm,vector &x,vector &y);
  void res_load_matrix (long eid,matrix &lm);

  void gl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void gnl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);
  void elem_volintegration_quant (long eid, integratedquant iq, long lcid, vector &iv);
  
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
		      vector &t,vector &coef,matrix &km);
  void edgenodeval (long edg,vector &nodval,double *list);

  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  void intpointval (long eid,vector &nodval,vector &ipval);
  
  void define_meaning (long eid);
  
  
  /* termitovo */
  void ntdbr_vector (long eid,vector &ntdbr);
  void ntn_matrix (long eid,matrix &ntn);
  double compute_error (long eid, double &e2, double &u2, double &sizel, vector *rderfull);
  void elchar (long eid, matrix &spsig);
  /* termitovo */
  
  void extract (atsel &at,vector &val);
  void mechq_nodval (long eid,vector &nodval,nontransquant qt);
  void mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt);


  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
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
    //strastrestate ssst;

};

#endif
