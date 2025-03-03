#ifndef TABLEFUNCT_H
#define TABLEFUNCT_H

#include <stdio.h>
#include "iotools.h"


enum interpoltype{piecewiselin=1,piecewiseconst=2,lagrange=3,piecewiselin2=4,dirac2=5};///<supported types of interpolation
const enumstr interpoltypestr[5] = {{"linear",1}, {"piecewiseconst",2}, {"lagrange",3}, {"linear2",4}, {"dirac2", 5}};
const kwdset interpoltype_kwdset(sizeof(interpoltypestr)/sizeof(*interpoltypestr), interpoltypestr);

/**
   the class tablefunct deals with function given by a table
   it means, the function is defined in discrete points only
   several interpolation techniques can be used
   it computes also first derivative
   
*/
class tablefunct
{
 public:
  tablefunct();
  tablefunct(long n);
  tablefunct(long n, long it);
  tablefunct(double *px, double *py, long n, interpoltype it=piecewiselin);
  ~tablefunct();

  void read  (XFILE *in);
  void read_data_file(XFILE *in);
  void read_data_file2(XFILE *in);
  void print(FILE *out);
  void print_data_file(FILE *out);
  void print_data_file2(FILE *out);
  void read_prop(FILE *in);
  void readval (XFILE *in);
  void datacheck ();

  double getval(double temp);
  double getinvval(double temp);
  double getval2(double temp, double &k);
  double getval3(double temp, double &k);
  
  double getderiv (double temp);
  double inv_derivative (double temp);
  
  void copy(tablefunct &tf);
  long compare(tablefunct &tf);

  ///  type of interpolation
  interpoltype itype;
  ///  number of function values
  long asize;
  ///  array of independent variable (usually denoted x)
  double *x;
  ///  array of function values (dependent variable, usually denoted y)
  double *y;
  ///  data file for reading
  char *file1;
  /// flag for x and y array deallocation
  bool dealloc;
  
 private:
  double piecewise_const_interpol (double temp);
  double piecewise_linear_interpol (double temp);
  double inverse_piecewise_linear_interpol (double temp);

  double diracinterpol2(double temp);
  double diracinterpol2(double temp, double &k);
  double lininterpol2  (double temp, double &k);
  double lininterpol3  (double temp, double &k);
  double laginterpol  (double temp);
  
  double derivative (double temp);
  
};
#endif
