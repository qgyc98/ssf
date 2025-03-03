#ifndef EDGEM_H
#define EDGEM_H

#include <stdio.h>
#include "alias.h"

/**
   class mechanical edges
   
   this class is used in hemivariational inequalities
   
   this class is strongly connected with the class gedge in GEFEL
   class gedge contains node numbers, previous and next edges,
   first and last nodes, direction and normal %vectors
   the class gedge in GEFEL contains problem independent informations
   
   the class edgem contains problem dependent informations
   
   JK, 8.8.2007
*/

class edgem
{
 public:
  edgem (void);
  ~edgem (void);
  void read (FILE *in);
  void init (long edid);

  void nodal_displacements (long lcid);
  void tan_nor_displacements ();
  void compute_jumps ();
  void compute_edge_functional_derivative (long lcid,double *v);


  ///  number of nodes on edge
  long nn;
  
  ///  number of approximated functions
  long napfun;
  
  ///  number of assigned general edge from GEFEL
  long ned;
  
  ///  displacements are discontinuous along the edge, therefore displacements
  ///  from each side of the edge are required in order to compute the jump between them
    
  ///  nodal displacements in the global coordinate system
  double u1,u2,u3,u4,v1,v2,v3,v4;

  ///  tangentional displacements (displacements in the direction defined by the edge)
  double td1,td2,td3,td4;
  
  ///  normal displacements (displacements normal to the direction defined by the edge)
  double nd1,nd2,nd3,nd4;
  
  ///  jumps in the tangential and normal directions
  ///  at each end of the edge
  double jt1,jt2,jn1,jn2;
    
  
  ///  auxiliary array for nodal displacements
  ///  this class serves only for 2D problems
  ///  it means that nodes contain 2 DOFs
  double *r;

  ///  type of material
  ///  in the case of hemivariational inequalities, the material model describes the multifunction b
  mattype *tm;
  ///  number of appropriate material type
  long *idm;

};

#endif
