#ifndef AXISYMC_H
#define AXISYMC_H

#include "alias.h"
#include "aliasc.h"
//#include "aliast.h"
//#include "genfile.h"
#include "matrix.h"

/**
   class of quadrilateral 8-node finite element for hydro-mechanical coupling
   
   JK, 2019
*/

class axisymc
{
 public:
  axisymc (void);
  ~axisymc (void);
  
  void codnum (long *cn,long ri);
  
  void stiffness_matrix (long lcid,long eid,matrix &sm);
  void conductivity_matrix (long lcid,long eid,matrix &km);
  void c_pp_matrix (long eid,matrix &cm);
  void c_up_matrix (long eid,matrix &cm);
  void c_pu_matrix (long eid,matrix &cm);

  void zero_order_matrix (long eid,matrix &km);
  void first_order_matrix (long eid,matrix &cm);
  
  void ip_strains (long eid);
  
  ///  the number of DOFs
  ///  there are 8 nodes in mechanics and 4 nodes in transport process
  ///  therefore 8*2+4=20
  long ndofe;
  ///  the number of DOFs in mechanics
  long mndofe;
  ///  the number of DOFs in transport
  long tndofe;
  
  ///  the number of nodes for mechanics
  long mnne;
  ///  the number of nodes for transport
  long tnne;
  
  ///  the number of edges
  long ned;
  ///  the number of nodes on an edge for transport
  long tnned;

  ///  order of integration for linear approximation functions
  long intordlin;
  ///  order of integration for quadratic approximation functions
  long intordquad;
  
  ///  the number of displacement components
  long ncompdispl;
  ///  the number of strain components
  long ncompstr;
  ///  stress/strain state
  strastrestate ssst;
  ///  the number of transported media
  long ntm;
  ///  the number of gradient/flux components
  long ncompgrad;

  ///  the number of integration points
  long nip;
  

  ///  unknown ordering
  long **ordering;
  
};

#endif
