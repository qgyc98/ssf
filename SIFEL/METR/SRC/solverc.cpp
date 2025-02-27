#include "solverc.h"
#include "globalc.h"
#include "fcsolver.h"
#include "pcsolver.h"
#include "cpcsolver.h"

void solve_metr_problem ()
{
  
  switch (Cp->tprob){
  case fully_coupled_mech_trans:{
    solve_fcouplprob ();
    break;
  case par_coupl_mech_trans:{
    solve_pcouplprob ();
    break;
  }
  case growing_par_coupl_mech_trans:{
    //solve_gpcouplprob ();
    solve_pcouplprob ();
    break;
  }
  }
  default:{
    print_err("unknown problem type is required in METR",__FILE__,__LINE__,__func__);
  }
  }
}
