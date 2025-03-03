#ifndef PVALT_H
#define PVALT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"

class pvalt
{
 public:
  pvalt ();
  ~pvalt ();
  long read(XFILE *in);
  long print(FILE *out);
  double getval(double time);
  double getprevval(double t);
  long read_prop(FILE *in, long lc);
  double getipv(void);

  //  type of time function
  generalfunct tfunc;
  // prescribed value
  double v;
  // inital value
  double ipv;

  //  strings with expressions
  char  func[256];
  //  parser binary trees with expressions
  Equation *eq;
  //  array with pointers to expression variables
  variable *var;
  
  //  function defined by table
  tablefunct tabf;

  //number of functions multipling prescribed values
  long nmultpv;
};

#endif
