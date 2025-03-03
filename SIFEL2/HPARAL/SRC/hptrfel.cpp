#define EXTERN
#include "hpglobal.h"
#include "seqfilest.h"
#include "hpsolvert.h"
#include "hptrfelinit.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main (int argc,char *argv[])
{
  time_t bt,et;
  XFILE *in;
  
  bt = time (NULL);

  //  MPI initialization
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  MPI_Get_processor_name(proc_namet,&nameLengtht);

  //  initialization of the code
  hptrfel_init (argc,(const char **)(argv));
  //hptrfel_init_tiles (argc,(const char **)(argv),stochd);
  //hptrfel_init_tiles_parsolver (argc,(const char **)(argv),stochd, in);
 
  //  solution of the problem
  solve_trfel_problem_parallel (in);
  
  et = time (NULL);
  
  //fprintf (Outt,"\n\n\n Udaje o dobach vypoctu \n");
  //fprintf (Outt,"\n\n celkova doba vypoctu             %ld",et-bt);
  if (Myrank==0){
    fprintf (stdout,"\n\n\n Udaje o dobach vypoctu \n");
    fprintf (stdout,"\n\n celkova doba vypoctu             %ld",et-bt);
  }
  
  //fprintf (Outt,"\n");
  //fclose (Outt);
  
  fprintf (stdout,"\n");  fprintf (stderr,"\n");
  
  MPI_Finalize ();
  
}
