#ifndef AXIQUADC_H
#define AXIQUADC_H

#include "alias.h"
#include "aliast.h"
#include "genfile.h"

/**
   class of quadrilateral 8-node finite element for thermo-mechanical coupling
   
   JK
*/

class axiquadc
{
 public:
  axiquadc (void);
  ~axiquadc (void);
  
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
