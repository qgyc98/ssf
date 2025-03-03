#ifndef LOADN_H
#define LOADN_H

#include <stdio.h>
#include "iotools.h"
struct vector;
struct atsel;



/**
  The class stores data about nodal load.

  Created by JK,
  Modified by Tomas Koudelka,
*/
class loadn
{
 public:
  loadn();
  ~loadn();
  long read(XFILE *in);
  long print(FILE *out);
  long read_prop(FILE *in, long ndof, long lc);
  void assemble (double *rhs);
  void changeparam (atsel &atln,vector &val);
  
  ///  number of loaded node
  long nid;
  /// load case number
  long nlc;
  /// array of load components
  double *f;
};

#endif
