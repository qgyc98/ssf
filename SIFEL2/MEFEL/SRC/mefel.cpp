#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#define EXTENDED_GLOBINC
#include "global.h"
#include "globalg.h"

#include "mefelinit.h"
#include "solverm.h"

// for evaluation of hypoplasticity
#include "intpoints.h"
//#include "hypoplunsatexptherm.h"
#include "hypoplunsatexptherm2.h"
//#include "hypoplunsatexptherm2.h"


int main (int argc, const char *argv[])
{
  time_t bt,et;
  long hod,min;
  double sec = clock();

  // start time of the job
  bt = time (NULL);
  
  //  program initiation and data reading
  mefel_init (argc, argv);

  //  solution of the problem
  //  the most important function in the code
  solve_mefel_problem ();
  
  // end time of the job
  et = time (NULL);
  
  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Out,"\n\n\n Data about computation time \n");
  fprintf (Out,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Out,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);
  
  fprintf (stdout,"\n -------------------------------------------");
  fprintf (stdout,"\n Consumed processor time by MEFEL %ld:%ld:%5f", hod, min, sec);
  fprintf (stdout,"\n -------------------------------------------\n");
  
  fprintf (stdout,"\n ---------------------------------------");
  fprintf (stdout,"\n Consumed total time by MEFEL %ld", long(et-bt));
  fprintf (stdout,"\n ---------------------------------------\n");

  /*  
  fprintf (stdout,"\n ----------------------------------------------");
  fprintf (stdout,"\n Peak of consumed memory of matrix class %ld",give_ammax());
  fprintf (stdout,"\n Peak of consumed memory of vector class %ld",give_avmax());
  fprintf (stdout,"\n ----------------------------------------------\n");
  */
  #ifdef INC_OPENMP
   fprintf (stdout, "\n\n Time consumed by multithread procedures: %le\n", Omp_wtime);
   fprintf (Out, "\n Time consumed by multithread procedures: %le\n", Omp_wtime);
  #else
   fprintf (Out, "\n Time consumed by singlethread procedures: %le\n", Omp_wtime);   
  #endif
   
  if (Mm->hypoplustherm)
  {
    fprintf(Out, "\n Total number of hypoplasticity model evalutaion in RKF is %ld\n\n", long(Neval));
    fprintf(Out, "\n Total number of hypoplasticity model evalutaion in RKF is %lg\n\n", Neval);
  }
  fprintf (Out,"\n -------------------------------------------");
  fprintf (Out,"\n Consumed processor time by MEFEL %ld:%ld:%5f",hod,min,sec);
  fprintf (Out,"\n -------------------------------------------\n");

  fprintf (stdout,"\n\n");  fprintf (stderr,"\n\n");  fprintf (Out,"\n");
  fclose (Out);

  // deallocate all global variables
  delete_glob();

  return 0;
}
