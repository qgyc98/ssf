#include "mpi.h"
#include "pllssolver.h"
#include "pglobal.h"
#include "seqfilesm.h"
#include "genfile.h"
#include <string.h>

void par_solve_layered_linear_statics ()
{
  long i;
  double *lhs,*rhs,*th;
  
  //  actual load case
  i=0;
  
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);

  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  array containing thicknesses
  th = new double [Mt->nn];
  for (i=0;i<Mt->nn;i++){
    th[i]=Mc->give_onethickness (Mt->nodes[i].crst,Mt->nodes[i].idcs);
  }
  
  //  constraint matrix assembling
  Psolm->constr_mat (th,Out);
  
  delete [] th;
  
  
  rhs = Lsrs->give_rhs (0);
  
  Mb->lc[i].assemble (i,rhs+i*Ndofm,NULL,1.0);
  
  
  //  solution of equation system
  Psolm->par_linear_solver (Gtm,Smat,lhs,rhs,Out,Mespr);
  
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){

    compute_req_val (i);

    //  output data
    //Lsrs->output (Out,i);
  }
  
}

