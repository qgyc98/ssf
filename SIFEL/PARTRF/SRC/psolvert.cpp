#include "pglobalt.h"
#include "psolvert.h"
#include "pspsolvert.h"
#include "pnpsolvert.h"
#include "pnnpsolvert.h"
#include "phsolvert.h"
#include "iotools.h"

void par_solve_trfel_problem ()
{
  switch (Tp->tprob){
  case stationary_problem:{
    par_solve_stationary_problem ();
    break;
  }
  case nonstationary_problem:{
    par_solve_nonstationary_problem ();
    break;
  }
  case nonlinear_nonstationary_problem:{
    //par_solve_nonlinear_nonstationary_problem_dform ();
    //par_solve_nonlinear_nonstationary_problem ();
    par_solve_nonstationary_problem ();//debug!!!???
    break;
  }
  default:{
    par_print_err (Myrank, proc_name,"unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}
