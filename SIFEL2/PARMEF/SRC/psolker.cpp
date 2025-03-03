#define EXTERN
#include "pglobal.h"
#include "seqfilesm.h"
#include "psolverm.h"
#include "pmefelinit.h"

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

int main (int argc,char *argv[])
{
  time_t bt,et;
  char name[1100];

  
  bt = time (NULL);

  //  MPI initialization
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  

  pmefel_init ();
  
  
  //  solution of the problem
  par_solve_mefel_problem ();
  

  //  testovaci volani
  //Psol->nodesplit (Gtm,Out);



  et = time (NULL);
  
  fprintf (Out,"\n\n\n Udaje o dobach vypoctu \n");
  fprintf (Out,"\n\n celkova doba vypoctu             %ld",et-bt);
  if (Myrank==0){
    fprintf (stdout,"\n\n\n Udaje o dobach vypoctu \n");
    fprintf (stdout,"\n\n celkova doba vypoctu             %ld",et-bt);
  }
  
  
  fprintf (Out,"\n");
  fclose (in);  fclose (Out);
  
  fprintf (stdout,"\n");  fprintf (stderr,"\n");
  
  MPI_Finalize ();
  
}
