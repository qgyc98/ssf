#define EXTERN
#include <mpi.h>
#include "pglobal.h"
#include "globalg.h"
#include "seqfilesm.h"
#include "psolverm.h"
#include "pmefelinit.h"
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
  
  //  MPI initiation
  MPI_Init(&argc,(char***)&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  MPI_Get_processor_name(proc_name,&nameLength);
    
  fprintf (stdout,"\n MYRANK je %d on %s",Myrank,proc_name);

  //  initiation of the code
  pmefel_init (argc, argv);
  
  //  solution of the problem
  par_solve_mefel_problem ();
  
  // end time of the job
  et = time (NULL);
  
  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Out,"\n\n\n Data about computation time \n");
  fprintf (Out,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Out,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);
  
  if (Myrank==0){
    fprintf (stdout,"\n -------------------------------------");
    fprintf (stdout,"\n Consumed time by PMEFEL %2ld:%02ld:%05.2f",hod,min,sec);
    fprintf (stdout,"\n -------------------------------------\n");
  }

  fprintf (stdout,"\n");  fprintf (stderr,"\n");  fprintf (Out,"\n");
  fclose (Out);
  
  MPI_Finalize ();
  
  // deallocate parallel global variables
  delete Pmp;
  delete Psolm;
  
  // deallocate all sequential global variables
  delete_glob();

  return 0;
}
