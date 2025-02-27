#define EXTERN
#include <mpi.h>
#include "hpglobal.h"
#include "globalg.h"
#include "seqfilesm.h"
#include "hpsolverm.h"
#include "hpmefelinit.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

int main (int argc,char *argv[])
{
  time_t bt,et;
  
  bt = time (NULL);

  //  MPI initialization
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  MPI_Get_processor_name(proc_namet,&nameLengtht);

  //  initialization of the code
  hpmefel_init (argc,(const char **)(argv));
  //hpmefel_init_tiles (argc,(const char **)(argv));
  
  //  solution of the problem
  solve_mefel_problem_parallel ();

  et = time (NULL);
  
  
  //fprintf (Out,"\n\n\n Udaje o dobach vypoctu \n");
  //fprintf (Out,"\n\n celkova doba vypoctu             %ld",et-bt);
  if (Myrank==0){
    fprintf (stdout,"\n\n\n Udaje o dobach vypoctu \n");
    fprintf (stdout,"\n\n celkova doba vypoctu             %ld",et-bt);
  }
  
  //fprintf (Out,"\n");
  //fclose (Out);
  
  fprintf (stdout,"\n");  fprintf (stderr,"\n");
  
  MPI_Finalize ();
  
}
