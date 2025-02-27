#include "hpsolvert.h"
#include "hpnpsolvert.h"
#include "hpnnpsolvert.h"
#include "hptrftiles.h"
#include "hpglobal.h"
#include "globalt.h"
#include "iotools.h"

void solve_trfel_problem_parallel (XFILE *in)
{
  //parallel_trfel_tiles ();
  //parallel_trfel_tiles_parsolver (in);
  
  switch (Mtprob){
  case stationary_problem:
  case nonstationary_problem:{
    parallel_homogenization_lin_nonstat ();
    break;
  }
  case nonlinear_nonstationary_problem:{
    parallel_homogenization_nonlin_nonstat ();
    break;
  }
  case trfel_tiles:{
    parallel_trfel_tiles ();
    break;
  }
  default:{
    par_print_err(Myrank,proc_namet,"unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}
