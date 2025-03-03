#ifndef CONVERTOR_H
#define CONVERTOR_H
#define EXTERN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "prepalias.h"
#include "descrip.h"
#include "globprep.h"
#include "elemhead.h"
#include "kwdset.h"
#include "stacktrace.h"
#include "siftop.h"
#include "gtopology.h"


class connectmeshes{
public:

  // inputline paramerts
  int Argc;
  char **Argv;
  
  XFILE *in,*topf;
  long err;
  descrip d;
  char *logname;
  char *tmp;
  gtopology *Gtop;
  siftop *SifTop;
  long *ltg;
  
  long *ninterfnodes;
  long **interfnodes;
  
  char outputfname[1025]; 
  char topfile[1025]; 
  long nfiles;
  double coordtol;
 
  connectmeshes(int argc,char *argv[]);
  ~connectmeshes();
  long input_siftop(long nfil);
  int read_input_data();
  int read_topology(long nfil);
  void establish_ltg();
  void create_glued_numbering();
  void print();
  void run();

};

#endif

