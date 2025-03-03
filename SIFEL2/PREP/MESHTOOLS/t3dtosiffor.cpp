#include <stdlib.h>
#include <stdio.h>
#include "siftop.h"

int main (int argc, char *argv[])
{
  siftop *top;
  XFILE *in;
  FILE *out;
  long paral;
  
  if (argc != 4){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s output_file_T3d file_name_SSMF paral(0=seq/1=paral)\n\n", argv[0]);
    return(1);
  }

  in = xfopen(argv[1], "rt");
  if (in==NULL){
    fprintf (stderr,"\n Output file from T3D has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  
  top = new siftop;
  paral = atol(argv[3]);
  top->import_t3d(in, paral);
  xfclose(in);
  
  out = fopen(argv[2], "wt");
  if (out==NULL){
    fprintf (stderr,"\n File with SSMF has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  
  
  top->print (out);
  
  fclose (out);
}
