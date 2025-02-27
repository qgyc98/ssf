#ifndef CONVERTOR_H
#define CONVERTOR_H
#define EXTERN

#include "prepalias.h"
#include "convertoralias.h"
#include "descrip.h"
#include "globprep.h"
#include "elemhead.h"
#include "kwdset.h"
#include "stacktrace.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


 


class convertor{
public:
  

  convertor(int argc, char *argv[]);
  ~convertor();
  int read_input_data();
  int read_topology();
  int run();
  void run_t3d2sifel();
  void run_sifel2gid();
  void run_t3d2gid();
  
  // variables and pointers
  
  // inputline paramerts
  int Argc;
  char **Argv;
  
  XFILE *in,*topf;
  long err;
  descrip d;
  char *logname;
  char *tmp;
  

  long convprocessing;
  char outputfname[1025]; 
  char topfile[1025]; 
  long convType;
  long ndom;


};

#endif
