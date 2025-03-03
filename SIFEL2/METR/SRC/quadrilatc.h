#ifndef QUADRILATC_H
#define QUADRILATC_H

#include "genfile.h"
#include "alias.h"
#include "aliast.h"

/**
   class quadrilatc defines quadrilateral element for problems of mechanical-transport coupling
   
*/

class quadrilatc
{
 public:
  quadrilatc (void);
  ~quadrilatc (void);
  //void eleminit (long eid);
  void dmatblockcol (long ri,long ci,matrix &d, matrix &dd);
  void dmatblockrow (long ri,long ci,matrix &d, matrix &dd);

  void intpointval (long eid);
  void intpointgrad (long eid);
  void mainip_strains (long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_mainip_strains (long eid);

  void upper_cond_coup_matrix (long eid,long ri,long ci,matrix &vm);
  void lower_cond_coup_matrix (long eid,long ri,long ci,matrix &vm);
  void upper_cap_coup_matrix (long eid,long ri,long ci,matrix &vm);
  void lower_cap_coup_matrix (long eid,long ri,long ci,matrix &vm);

  void res_upper_cond_coup_matrix (long eid,matrix &vm);
  void res_lower_cond_coup_matrix (long eid,matrix &vm);
  void res_upper_cap_coup_matrix (long eid,matrix &vm);
  void res_lower_cap_coup_matrix (long eid,matrix &vm);

  void upper_cond_coup_vector (vector &tvm,vector &nodval,long eid,long ri,long ci);
  void res_upper_cond_coup_vector (vector &f,long eid);
  
  void upper_internal_forces (long ri,long ci,long eid,vector &ifo);
  void res_upper_internal_forces (long eid,vector &ifor);
  void lower_internal_fluxes (long ri,long ci,long eid,vector &ifl);
  void res_lower_internal_fluxes (long eid,vector &elemif);

  //void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  //void res_volume_rhs_vector (vector &f,long eid,long lcid);

  //void mefel_metr (long eid);
  //void trfel_metr (long eid);
  //void trfel_mefel (long eid);

  ///  number of degrees of freedom of mechanical part
  long mndofe;
  ///  number of degrees of freedom of transport part
  long tndofe;
  ///  number of nodes on one element in mechanical problem
  long nnemp;
  ///  number of nodes on one element in transport problem
  long nnetp;
  ///  number of blocks in strain %vector
  long mnb;
  ///  number of transported media
  long ntm;
  ///  number of components of strain %vector
  long tnmcomp;
  ///  array of numbers of components in strain vectors
  long *mncomp;
  ///  stress/strain state
  strastrestate ssst;


  ///  orders of integration of upper coupling %matrices
  long **intordvum;
  ///  orders of integration of lower coupling %matrices
  long **intordvlm;
  ///  numbers of integration points of upper coupling %matrices
  long **nipu;
  ///  numbers of integration points of lower coupling %matrices
  long **nipl;
  ///  number of DOFs for particular medium
  long *dofe;
  ///  array containing ordering of mechanical unknowns
  long *mordering;
  ///  array containing ordering of transport unknowns
  long **tordering;
  /// total number of integration point
  long tnipu,tnipl;

};

#endif
