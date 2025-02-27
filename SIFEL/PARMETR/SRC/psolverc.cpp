#include "psolverc.h"
#include "pglobalc.h"
#include "globalc.h"
#include "ppcsolver.h"

/**
   parallel solution of coupled problem
*/
void par_solve_metr_problem ()
{
  switch (Cp->tprob){
  case fully_coupled_mech_trans:{
    //  fully coupled hydro-thermo-mechanical problem
    //par_solve_linear_statics ();
    break;
  }
  case par_coupl_mech_trans:{
    //  partially coupled hydro-thermo-mechanical problem
    par_solve_pcouplprob ();
    break;
  }
  default:{
    par_print_err(Myrank,"unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}
