#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "intools.h"

int main (int argc,char *argv[])
{

  long i, ndom, end;
  FILE *in;
  FILE *out;
  char ofname[2048];
  char *dot;
  char buf[1025];
  char *ptr;

  fprintf (stderr,"\n\n *** PARMEFEL GENERATOR OF LAYERED PROBLEMS ***\n");
  fprintf (stderr," ----------------------------------------------\n");
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name number_of_domains\n\n", argv[0]);
    return(1);
  }
  in = fopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  if (sscanf(argv[2], "%ld", &ndom) != 1)
  {
    fprintf (stderr,"\n Number of domains has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  for (i = 0; i < ndom; i++)
  {
    dot = strrchr(argv[1], '.');
    if (dot == NULL)
    {
      dot = argv[1]+strlen(argv[1]);
      strncpy(ofname, argv[1], size_t(dot-argv[1]));
      sprintf(ofname+(dot-argv[1]), "%ld.pr", i+1);
    }
    else
    {
      strncpy(ofname, argv[1], size_t(dot-argv[1]));
      sprintf(ofname+(dot-argv[1]), "%ld", i+1);
      strcpy(ofname+strlen(ofname), dot);
    }
    out = fopen(ofname, "wt");
    if (out == NULL)
    {
      fprintf(stderr, "\n Output file could not be opened\n\n");
      return(4);
    }
    end = 0;
    fseek(in, 0L, SEEK_SET);
    do
    {
      inputln(in, buf, 1024);
      switch (isgkw(buf, ptr))
      {
        case -1 :
          fprintf(out, "%s\n", buf); // normal line
          break;
        case 0 :
          fprintf(out, "%s\n", buf); // string with 'begin' keyword
          fprintf(out, "%ld\n", i+1);// domain number
          end = 1;
          break;
      }
    } while (! end);

    while ( !feof(in))
    {
      inputln(in, buf, 1024);
      if ( !feof(in))
        fprintf(out, "%s\n", buf);
      else
        break;
    }

    fclose (out);
  }
  fprintf (stdout,"\n");
}
