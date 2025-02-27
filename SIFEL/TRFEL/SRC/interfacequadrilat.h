#ifndef CONTACTQUADRILAT_H
#define CONTACTQUADRILAT_H

#include "aliast.h"
#include "genfile.h"

/**
   Class interfacequadrilat is quadrilateral element for the simulation of interfaces.
   The element is considered for connection of edges of 2D elements in TRFEL.
   
   basic data
   nne = 4 - number nodes on element
   ndofe = 4 - number of DOFs on element
   
   JK, 4. 1. 2016
*/
class interfacequadrilat
{
 public:
  interfacequadrilat (void);
  ~interfacequadrilat (void);
  
  void codnum (long *cn,long ri);

  void grad_matrix (matrix &gm,vector &x,vector &y,double &l,double h);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  
  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void res_capacity_matrix (long eid,matrix &cm);
  
  void transq_nodval(long eid, vector &nodval, nonmechquant nmq);
  void transq_init_nodval (long eid,vector &nodval,nonmechquant nmq);

  ///  number of transported matters
  long ntm;
  ///  total number of DOFs on the element
  long ndofe;
  ///  numbers of DOFs for particular problems
  long **dofe;
  ///  number of nodes on one element
  long nne;

  ///  problem dimension
  long ncomp;
  ///  number of integration points
  long **nip;
  ///  unknown ordering
  long **ordering;

};

#endif
