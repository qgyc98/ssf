#ifndef STRRETALG_H
#define STRRETALG_H

#include "iotools.h"
#include "alias.h"



/**
  Class defines type of stress return algorithm.
  It contains all necessary informations about selected algorithm.
   
  Created by JK,
*/
class strretalg
{
 public:
  strretalg ();
  ~strretalg ();
  void read (XFILE *in);
  void print (FILE *out);
  stressretalgtype give_tsra ();
  long give_ni ();
  double give_err ();
  double give_err_sig ();
  void give_err_stv (double &err_req, double &err_max);
  double give_hmin();
  rktype give_rkt();
  
  ///  type of stress return algorithm
  stressretalgtype tsra;
  ///  computer zero
  double zero;
  ///  maximum number of iterations
  long ni;
  ///  required accuracy
  double err;
  ///  required accuracy of state variables
  double err_stv_req;
  ///  maximum required accuracy of state variables
  double err_stv_max;
  /// minimum step length
  double h_min;
  /// type of Runge-Kutta method
  rktype rkt;
};

#endif
