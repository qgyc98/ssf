#ifndef PARTRFINPUTFILES_H
#define PARTRFINPUTFILES_H
#define EXTERN
#include "prepalias.h"
#include "descript.h"
#include "globprept.h"
#include "elemheadt.h"
#include "kwdset.h"
#include "stacktrace.h"

class partrfinputfile{
public:
  // variables and pointers
  
  // inputline paramerts
  int Argc;
  const char **Argv;
  
  XFILE *topf;
  long err;
  descript d;
  char *logname;
  char *tmp;
 
  //functions
  partrfinputfile(int argc, const char *argv[]);
  ~partrfinputfile();
  
  int run();
  // function read input file
  int read_topology();
  // function for pritin partrf input file
  int print_partrfinputfile();
  void print_in(FILE *out);
  void print_partrf_in(FILE *out);
  void partrftop_print (FILE *out);
  void partrfload_print (FILE *out);
  void partrfloadcase_print (FILE *out,long lcid);
  void partrfelemload_print(long lcid,long leid,long deid,FILE *out);
  void partrf_initcondprint (FILE *out);
  long give_elemtype(long etype);
};

#endif
