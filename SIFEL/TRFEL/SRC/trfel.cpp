#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#include "globalt.h"
#include "globalg.h"
#include "trfelinit.h"
#include "solvert.h"

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
  trfel_init (argc,argv);
  
  
  //  solution of the problem
  //  the most important function in the code
  solve_trfel_problem ();


  // end time of the job
  et = time (NULL);
  
  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;

  fprintf (Outt,"\n\n\n Data about computation time \n");
  fprintf (Outt,"\n Total time of computation      %ld", long(et-bt));
  fprintf (Outt,"\n Processor time of computation  %ld:%ld:%.5f", hod, min, sec);

  fprintf (stdout,"\n ----------------------------------");
  fprintf (stdout,"\n Consumed processor time by TRFEL %2ld:%02ld:%05.2f",hod,min,sec);
  fprintf (stdout,"\n ----------------------------------\n");

  //fprintf (stdout,"\n Acm = %ld \t Acv = %ld \t Aciv = %ld \n",give_acm(),give_acv(),give_aciv());

  fprintf (stdout,"\n\n");  fprintf (stderr,"\n\n");  fprintf (Outt,"\n");
  fclose (Outt);
  //fclose (Outt1);
  //fclose (Outt2);

  // deallocate all global variables
  delete_globt();

  return 0;
}
