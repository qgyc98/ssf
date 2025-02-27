#ifndef SPRINGEL_H
#define SPRINGEL_H

#include "alias.h"
#include <stdlib.h>
#include "matrix.h"

struct vector;
struct ivector;


/**
   CLASS SPRINGEL
   This class defines special type of onedimensional element which is used
   as spring support in the different directions. The direction of the support
   is given by the element type and can be one of these :

   created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz,18.9.2002
*/
class springel
{
 public:
  springel (void);
  ~springel (void);
  long give_ndofe(long eid);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void strains (long eid,long lcid);
  void stresses (long eid,long lcid);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void ipcoord(long eid, long ipp, vector &ipcoord);

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  order of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  array of numbers of integration points in blocks
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
