#ifndef MESHDECOMP_H
#define MESHDECOMP_H
#define EXTERN
#include "prepalias.h"
//#include "input.h"
#include "descrip.h"
#include "globprep.h"
#include "elemhead.h"
#include "kwdset.h"
#include "stacktrace.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "partalias.h"
#include <time.h>
#include "metis.h"
//#include "party_lib.h"
 


class meshdecomp{
public:
  // variables and pointers
  
  // inputline paramerts
  int Argc;
  char **Argv;
  
  XFILE *in,*topf;
  long err;
  descrip d;
  char *logname;
  char *tmp;
  char weightf[1025];   
  char outputfname[1025]; 
  char graphf[1025]; 
  char partf[1025]; 
  char ubecfile[1025]; 
  char tpwgtsfile[1025]; 
  long partType;
  long processing;
  int nparts;
  long parttech;
  long opt;
  long graphWeight;
  long weightOpt;

   
  meshdescription md;
  
  int nvtxs; 
  int nedges;
  int *xadj;
   
  int *adjncy;
  //vertex weight
  int *vwgt;
  //edge weight
  int *adjwgt;
  
  int *part;
  
  // indicator for hybrid parallel homogenization aggregate printing
  int homog_par;
  
  //METIS
  int ncon;
  int *options;
  idxtype numflag, wgtflag, edgecut,volume;
  idxtype metis_nvtxs; 
  idxtype metis_nedges;
  idxtype *metis_xadj;
  float *ubec;
  float *tpwgts;
  
  idxtype *metis_adjncy;
  //vertex weight
  idxtype *metis_vwgt;
  //edge weight
  idxtype *metis_adjwgt;
  
  idxtype *metis_part;

  
  //Chaco
  double *x;
  double *y;
  double *z;
  
  //functions
  meshdecomp(int argc, char *argv[]);
  ~meshdecomp();
  
  int run();
  // function read input file
  int read_input_data();
  int read_topology();
  int run_metis();
  int run_chaco();
  void run_aggregates();
  void read_metis_option();
  void read_chaco_option();
  void part_postprocessing();
  void aggregate_postprocessing();
  void creation_dual_graph();
  void creation_nodal_graph();
  void print_metis_graph();
  void read_partitioning();
  void read_metis_ubec();
  void read_metis_tpwgts();
  void print_chaco_graph();
  void compute_centre_elements();
  void run_jostle();
  void print_jostle_graph();
  int run_party();
  void print_scotch_graph();
  void run_scotch();

};

#endif
