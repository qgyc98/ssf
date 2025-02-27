#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "globalg.h"
#include "solverc.h"
#include "metrinit.h"

#include "seqfilesm.h"
#include "seqfilest.h"


int main (int argc, const char *argv[])
{
  time_t bt,et;
  long hod,min;
  double sec = clock();
  
  //setvbuf(stdout, NULL, _IONBF, 0);
  //setvbuf(stderr, NULL, _IONBF, 0);

  // start time of the job
  bt = time (NULL);
  
  //  program initiation and data reading
  metr_init (argc,argv);
  
  //  solution of the problem
  //  the most important function in the code
  solve_metr_problem ();

  // end time of the job
  et = time (NULL);
  
  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Outc,"\n\n\n Data about computation time \n");
  fprintf (Outc,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Outc,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);

  fprintf (stdout,"\n -------------------------------------------");
  fprintf (stdout,"\n Consumed processor time by METR %2ld:%02ld:%05.2f",hod,min,sec);
  fprintf (stdout,"\n -------------------------------------------\n");
  #ifdef INC_OPENMP
   fprintf (stdout, "Time consumed by multithread procedures: %le\n", Omp_wtime);
   fprintf (Out, "Time consumed by multithread procedures: %le\n", Omp_wtime);
  #else
   fprintf (Out, "Time consumed by singlethread procedures: %le\n", Omp_wtime);
  #endif
  fprintf (Outc,"\n -------------------------------------------");
  fprintf (Outc,"\n Consumed processor time by METR %ld:%ld:%5f",hod,min,sec);
  fprintf (Outc,"\n -------------------------------------------\n");
  
  
  fprintf (stdout,"\n\n");  fprintf (stderr,"\n\n");  fprintf (Outc,"\n");
  fclose (Out);
  fclose (Outt);
  fclose (Outc);  
  
  // deallocate all global variables of METR
  delete_globc();
  // deallocate all global variables of MEFEL
  delete_glob();
  // deallocate all global variables of TRFEL
  delete_globt();

  return 0;
}
