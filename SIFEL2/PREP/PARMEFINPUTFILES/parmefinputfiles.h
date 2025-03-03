#ifndef PARMEFINPUTFILES_H
#define PARMEFINPUTFILES_H
#define EXTERN
#include "prepalias.h"
#include "descrip.h"
#include "globprep.h"
#include "elemhead.h"
#include "kwdset.h"
#include "stacktrace.h"

class parmefinputfile{
public:
  // variables and pointers
  
  // inputline paramerts
  int Argc;
  const char **Argv;
  
  XFILE *topf;
  long err;
  descrip d;
  char *logname;
  char *tmp;
 
  //functions
  parmefinputfile(int argc, const char *argv[]);
  ~parmefinputfile();
  
  int run();
  // function read input file
  int read_topology();
  // function for pritin parmef input file
  int print_parmefinputfile();
  void print_in(FILE *out);
  void print_parmef_in(FILE *out);
  void parmeftop_print (FILE *out);
  void parmefload_print (FILE *out);
  void dparmefloadcase_print (FILE *out,long lcid);
  void parmefloadcase_print (FILE *out,long lcid);
  void parmefsubloadcase_print (FILE *out,long lcid,long slcid);
  void parmefelemload_print(long lcid,long leid,long deid,FILE *out);
  void parmefselemload_print(long lcid,long slcid,long leid,long deid,FILE *out);
  long give_elemtype(long etype);
};

#endif
