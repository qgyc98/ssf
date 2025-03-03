#include "mpi.h"
#include "psolverm.h"
#include "plssolver.h"
#include "pnssolver.h"
#include "pmtsolver.h"
#include "pglobal.h"
#include "iotools.h"

void par_solve_mefel_problem ()
{
  switch (Mp->tprob){
  case linear_statics:{
    par_solve_linear_statics ();
    break;
  }
  case mat_nonlinear_statics:{
    
    break;
  }
  case mech_timedependent_prob:{
    par_solve_timemech_prob ();
    break;
  }
  default:{
    par_print_err(Myrank,"unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

}
