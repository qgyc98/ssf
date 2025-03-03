#define EXTERN

#include <mpi.h>
#include "pglobalc.h"
#include "pglobal.h"
#include "pglobalt.h"
#include "globalg.h"

#include "psolverc.h"
#include "pmetrinit.h"

#include "seqfilesm.h"
#include "seqfilest.h"
#include "seqfilesc.h"

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
  MPI_Get_processor_name(proc_name, &nameLength);
  
  // program initiation and data reading
  pmetr_init (argc, argv);
  
  //  solution of the problem
  par_solve_metr_problem ();
  
  // end time of the job
  et = time (NULL);

  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Outc,"\n\n\n Data about computation time \n");
  fprintf (Outc,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Outc,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);
  if (Myrank==0){
    fprintf (stdout,"\n ------------------------------------");
    fprintf (stdout,"\n Consumed time by PMETR %2ld:%02ld:%05.2f",hod,min,sec);
    fprintf (stdout,"\n ------------------------------------\n");
  }

  fprintf (stdout,"\n");  fprintf (stderr,"\n");  fprintf (Outc,"\n");
  fclose (Outc);

  // close MPI
  MPI_Finalize ();
  
  // deallocate parallel global variables
  delete Pcp;
  // deallocate all parallel global variables of MEFEL
  pdelete_glob();
  // deallocate all parallel global variables of TRFEL
  pdelete_globt();
  
  // deallocate all sequential global variables of METR
  delete_globc();
  // deallocate all sequential global variables of MEFEL
  delete_glob();
  // deallocate all sequential global variables of TRFEL
  delete_globt();

  return 0;
}
