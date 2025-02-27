#define EXTERN
#include <mpi.h>
#include "pglobalt.h"
#include "globalg.h"
#include "psolvert.h"
#include "ptrfelinit.h"
#include "seqfilest.h"
#include <stdio.h>
#include <string.h>
#include <time.h>



int main (int argc, const char *argv[])
{
  time_t bt,et;
  long hod,min;
  double sec = clock();
  
  // start time of the job 
  bt = time (NULL);

  memset(proc_name, 0, sizeof(*proc_name)*10000);

  //  MPI initialization
  MPI_Init(&argc,(char***)&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  MPI_Get_processor_name(proc_name, &nameLength);

  //  initialization of the code
  ptrfel_init (argc, argv);
  
  //  solution of the problem
  par_solve_trfel_problem ();

  // end time of the job
  et = time (NULL);
  
  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Outt,"\n\n\n Data about computation time \n");
  fprintf (Outt,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Outt,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);  
  if (Myrank==0){
    fprintf (stdout,"\n -------------------------------------");
    fprintf (stdout,"\n Consumed time by PTRFEL %2ld:%02ld:%05.2f",hod,min,sec);
    fprintf (stdout,"\n -------------------------------------\n");
  }
  
  fprintf (stdout,"\n");  fprintf (stderr,"\n");  fprintf (Outt,"\n");

  MPI_Finalize ();

  // deallocate parallel global variables
  delete Psolt;
  delete Ptp;
  
  // deallocate all sequential global variables
  delete_globt();
  fclose (Outt);

  return 0;
}
