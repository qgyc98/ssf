#ifndef ENDNODEM_H
#define ENDNODEM_H

#include <stdio.h>
#include "alias.h"
#include "iotools.h"

/**
   class endnodem defines general end node for mechanical problems
   
   this class is used in hemivariational inequalities
   
   this class is strongly connected with the class endnode in GEFEL
   class endnode contains node numbers, first and last nodes,
   the class endnode in GEFEL contains problem independent informations
   
   the class endnodem contains problem dependent informations
   
   JK, 1.3.2009
*/
class endnodem
{
 public:
  endnodem ();
  ~endnodem (void);
  
  void nodal_displacements (long lcid);
  void compute_jumps ();
  void compute_endnode_functional_derivative (long lcid,double *v);
  
  
  ///  number of approximated functions
  long napfun;
  
  ///  number of assigned general endnode from GEFEL
  long nen;
  
  ///  displacements are discontinuous at end node, therefore displacements
  ///  from each side of the end node are required in order to compute the jump between them
    
  ///  nodal displacements in the global coordinate system
  double u1,u2,v1,v2;

  ///  jumps in the displacement field
  double ju,jv;
      
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
