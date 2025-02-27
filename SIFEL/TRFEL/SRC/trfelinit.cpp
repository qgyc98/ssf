#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trfelinit.h"
#include "globalt.h"
#include "globmatt.h"
#include "intpointst.h"
#include "gtopology.h"
#include "transprint.h"
#include "iotools.h"
#include "xfile.h"
#include "stacktrace.h"
#include "libtrace.h"

/**
  The function initializes trfel from the input file.
  It reads data from input file and initiates all necessary variables and arrays.

  @param argc[in] - the number of arguments from command line, i.e. dimension of array argv   
  @param argv[in] - array of strings with arguments from command line   

  @return The function does not return anything but changes stae of global variables.
*/
void trfel_init (int argc, const char *argv[])
{
  long i;
  XFILE *in;
  pkwd_sw kwdsw;

  set_prgname(argv[0]); // sets name of program for stacktrace
  
  fprintf (stdout,"\n\n ******************************************************\n");
  fprintf (stdout,"            ____ ___ _____ _____ _\n");
  fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
  fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
  fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
  fprintf (stdout,"           |____/___|_|   |_____|_____| TRFEL");
  fflush(stdout);
  fprintf (stderr,"\n\n ******************************************************\n");

  // initialize global variables to null values
  initnull_globt();

  if (process_argst(argc, argv, kwdsw))
    abort();

  in = xfopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    abort ();
  }

  Tp = new probdesct;
  Tp->kwdsw = kwdsw; // set keyword processing mode obtained from command line arguments
  
  Gtt = new gtopology;
  Tt = new transtop;
  Tm = new transmat;
  Tc = new transcrsec;
  Tb = new transbclc;
  Outdt = new outdrivert;  

  Outt = fopen ("trfel.log","wt");

  in->warning    = 1;
  in->ignorecase = 1;
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_outd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;
  
  // **************************
  // **  basic problem data  **
  // **************************
  Tp->read (in);
  
  if (Tp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  
  // **************************************
  // **  input of node and element data  **
  // **************************************
  Tt->read (in);
  
  // mesh control TKr
  /*
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    if (Mesprt==1){
      merr=Tt->mesh_check();
      if (merr > 0){
        fprintf (stdout,"Mesh error = %ld\n\n\n",merr);
        exit(0);
      }
    }
  }
  */
  
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->read_gf (in);
    Gtt->alloc_growstr();
  }
  
  //  input of material characteristics
  Tm->read (in);
  
  //  detection of the radiation materials
  radiation_init ();
  
  
  if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem || Gtt->edtype==mater){
    Gtt->auto_subdom_detection (Outt);
    Tt->alloc_enodes ();
    Tt->alloc_edges ();
    Tt->edge_init_edval ();

    Tt->view_factors (Outt);
    Tt->mod_view_factors (Outt);
  }
  
  // *************************************
  //  generation of ordering of unknowns
  //  (generation of code numbers)
  //  number of degrees of freedom
  // *************************************
  Ndoft = Gtt->codenum_generation (Outt);
  if (Mesprt==1)
    fprintf (stdout,"\n the number of degrees of freedom    %ld",Ndoft);

  //  fprintf(stdout, "\n Element renumbering");
  //  Tt->sort_elements();
  //  integration point allocation and initiation
  //Tm->intpointalloc ();
  //Tm->intpointinit ();
  Tm->assemble_dof_nameord(Tp->dofname, Tp->ntm, Tp->var_dofid, Tp->tnkv);

  //  input of cross section characteristics
  Tc->read (in);
  
  //  input of boundary conditions
  Tb->read (in);

  //  prescribed values on elements investigation
  //Tt->elemprescvalue();
  //  searching for elements with nodes where source is defined
  Tb->elemsource ();
  
  //  global vectors allocations
  Lsrst = new lhsrhst;
  Lsrst->alloc ();
  //  initial conditions
  Lsrst->initcond (in);
  
  //  advective terms
  if (Tp->advect==1){
    Tm->advection_velocities (in); //05/06/2017 commented by TKr, because of missing function
    //  circumscribed balls are generated and the radii are stored
    Gtt->construct_circumscribed_balls ();
  }
  
  if (Tp->stochasticcalc){
    Stt = new stochdrivert;
    Stt->read (in);
  }
  else
    Stt = NULL;
  
  if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){

    Tt->enodes_init ();
    Tt->edge_init ();

    //  definition of material types and id on nodes
    //Tt->node_materials ();    
  }
  
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_pd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  // input of output description
  Outdt->read(in);

  if (Tp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  
  fclose (Outt);
  
  char fname[1020];
  if (*Outdt->outfn)
    filename_decomposition(Outdt->outfn,Tp->path,Tp->filename,Tp->suffix);
  else
    filename_decomposition(argv[1],Tp->path,Tp->filename,Tp->suffix);

  sprintf(fname, "%s%s.log", Tp->path, Tp->filename);
  Outt = fopen(fname,"wt");
  
  //  nodes allocation
  Tt->alloc_nodes ();
  
  
  //  some solvers require additional data
  //  if the Schur complement method or FETI method
  //  are used, the variable Ndofm is rewritten
  //  Ndofm is not changed in all other cases
  Tp->ssle->prepare_data (&Ndoft,Gtt,Outt);

  
  //if (Tp->tprob == discont_nonstat_problem){
  //Tt->node_materials ();
  //}

  //  initiation of initial values
  //added by TKr 19/9/2012
  if (Tm->initval==NULL)
    Tm->initval = new double [Tm->tnip];
    
  for (i=0;i<Tm->tnip;i++){
    Tm->initval[i]=0.0;
  }    
  
  xfclose (in);
}

void print_helpt(FILE *out, const char *prgname)
{
  fflush(out);
  fprintf(out, "\nInput file has not been specified\n");
  fprintf(out, " use : %s infile [-kwd=i]\n", prgname);
  fprintf(out, " where :\n *   infile is input file name\n");
  fprintf(out, " *   -kwd={0,1,2,3,4} is optional switch which controls keyword usage in the input file\n");
  fprintf(out, "           0 means no keywords (default value)\n");
  fprintf(out, "           1 means keywords only in the problem description section\n");
  fprintf(out, "           2 means keywords only in the outdriver section\n");
  fprintf(out, "           3 means keywords in the problem description and outdriver sections\n");
  fprintf(out, "           4 means full usage of keywords\n\n");
  abort();
}


long process_argst(int argc, const char *argv[], pkwd_sw &kwdsw)
{
  const char *ptrkwd;
  
  switch (argc)
  {
    case 0:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too few arguments are used\n\n");
      print_helpt(stderr, "ptrfel|trfel");
      return 1;
    case 1:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too few arguments are used\n\n");
      print_helpt(stderr, argv[0]);
      return 1;
    case 2:
      kwdsw = nokwd;
      break;
    case 3:
      // processing of optional keyword switch -kwd={0,1,2,3,4}
      ptrkwd = strstrcis(argv[2], "-kwd=");
      if (ptrkwd)
      {
        kwdsw = pkwd_sw(-1);
        sscanf(ptrkwd+5, "%d", (int *)&kwdsw);
        if ((kwdsw < 0) || (kwdsw > 4)){
          fflush(stdout);
          fprintf(stderr, "\n\nError - unknown -kwd value %d is obtained\n\n", kwdsw);
          print_helpt(stderr, argv[0]);
          return 1;
        }
      }
      else
      {
        fflush(stdout);
        fprintf(stderr, "\n\nError - unknown switch is used\n\n");
        print_helpt(stderr, argv[0]);
        return 1;
      }
      break;
    default:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too many arguments is used\n\n");
      print_helpt(stderr, argv[0]);
      return 1;
  }
  return 0;
}

/**
   function prepares auxiliary data for edge detection used for radiation problems
   function assigns value 1 to the variables auxinf on gelements
   
   JK, 3.6.2011
*/
void radiation_init ()
{
  long i,j,k,nmt;
  
  //  the number of material types used on a single element
  nmt=Tp->ntm*Tp->ntm;
  
  k=0;
  //  loop over the elements
  for (i=0;i<Tt->ne;i++){
    for (j=0;j<nmt;j++){
      if (Tt->elements[i].tm[j]==radiationmater){
	k=1;
	Gtt->edtype=mater;
	break;
      }
    }
    if (k==1){
      break;
    }
  }
  
  //  loop over the elements
  for (i=0;i<Tt->ne;i++){
    for (j=0;j<nmt;j++){
      if (Tt->elements[i].tm[j]==radiationmater){
	Gtt->gelements[i].auxinf=251;
	break;
      }
    }
  }
  
}

