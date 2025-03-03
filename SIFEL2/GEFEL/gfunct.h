#ifndef GFUNCT_H
#define GFUNCT_H

#include <stdio.h>
#include "iotools.h"
#include "parser.h"
#include "galias.h"
class tablefunct;
class itablefunct;
class gfunct;
struct vector;

/**
  This class represents general type of function which can be described
  by general string expression or by table with function values.
  String expression is processed by the Parser and table is processed by the
  tablefunct. In case of string expression, there is allowed to have several
  expressions for different ranges of the variable.

*/
class gfunct
{
 public:
  gfunct();
  gfunct(generalfunct tf,long nr);
  ~gfunct();
  
  long read(XFILE *in);
  double getval(double t);
  double getinvval(double t);
  double getval(vector &p, const char *namevar[]);
  double getval(double t, double t0);
  long getval_long (double t);
  double getderiv (double t);
  double getinvderiv (double t);

  long read_prop(FILE *in);
  long print(FILE *out);
  void initiate (gfunct &gf);

  void init_tab (long nr);
  /// copies function from the parameter to the actual object
  void copy (gfunct &gf);
  /// merges function from the parameter with the actual object
  void merge (gfunct &gf);
  /// compares function from the parameter to the actual object
  long compare (gfunct &gf);
  
  ///  type of general function
  generalfunct tfunc;
  /// number of expression sets
  long neqs;
  /// number of general functions
  long ngf;
  ///  constant value
  double f;

  ///  strings with expressions
  char  **func;
  ///  binary tree parser for an expressions
  Equation **eq;
  ///  binary tree parser for an expressions
  Equation **deq;
  ///  array with pointers arrays with pointers to expression variables
  variable ***var;
  /// array with limit values for pars_set
  double *limval;
  ///  function defined by table (all values are real numbers)
  tablefunct *tabf;
  ///  function defined by table (x values are real numbers while y values are integers)
  itablefunct *itabf;
  /// function is defined by set of general functions
  gfunct *gfs;
};

#endif
