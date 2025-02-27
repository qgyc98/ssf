#include "hpsolverm.h"
#include "hplssolver.h"
#include "hpmtsolver.h"
#include "hpglobal.h"
#include "global.h"
#include "iotools.h"

void solve_mefel_problem_parallel ()
{
  switch (Mtprobm){
  case linear_statics:{
    //parallel_homogenization_linear_statics_tiles ();
    parallel_homogenization_linear_statics ();
    break;
  }
  case mat_nonlinear_statics:{
    break;
  }
  case mech_timedependent_prob:{
    parallel_homogenization_solve_time_dep_prob ();
    break;
  }
  default:{
    par_print_err(Myrank,proc_namet,"unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}
