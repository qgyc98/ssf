#define EXTERN


#include "pglobal.h"
#include "seqfilesm.h"
#include "stochdriver.h"
//#include "mefelinit.h"
#include "hpssolver.h"

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <time.h>




int main (int argc, const char *argv[])
{
  time_t bt,et;
  char name[1100];
  stochdriver *stochd;
  stochd = new stochdriver;
  FILE *in;
  
  bt = time (NULL);
  
  //  MPI initialization
  MPI_Init(&argc,(char ***)(&argv));
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  
  
  //  program initiation and data reading
  mefel_init (argc, argv, stochd);
  
  /*
  if (Myrank==0){
    //  sample generation is only on the master processor
    stochd = new stochdriver;
    stochd->read (in);
    St = stochd;
  }
  */
  
  //  solution of the problem
  par_solve_stochastic_problem (stochd);
  
  /*
  et = time (NULL);
  
  fprintf (Out,"\n\n\n Udaje o dobach vypoctu \n");
  fprintf (Out,"\n\n celkova doba vypoctu             %ld",et-bt);
  if (Myrank==0){
    fprintf (stdout,"\n\n\n Udaje o dobach vypoctu \n");
    fprintf (stdout,"\n\n celkova doba vypoctu             %ld",et-bt);
  }
  
  
  fprintf (Out,"\n");
  fclose (in);  fclose (Out);
  */

  fprintf (stdout,"\n");  fprintf (stderr,"\n");
  
  MPI_Finalize ();
  
}
