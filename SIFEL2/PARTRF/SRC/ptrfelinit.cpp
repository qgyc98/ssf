#include "ptrfelinit.h"
#include "pglobalt.h"
#include "xfile.h"
#include "stacktrace.h"
#include <mpi.h>



/**
  Function reads data from input file and initiates all
  necessary variables and arrays.

  Parameters:   
  @param argc - number of passed command line parameters
  @param argv - array of pointers to strings with command line parameters passed to program 

  Returns:
  @retval The function does not return anything.

  Created by JK,
  Modified by Tomas Koudelka,
*/
void ptrfel_init (int argc, const char *argv[])
{
  char name[1100];
  char *dot;
  XFILE *in;
  pkwd_sw kwdsw;
  
  // set name of executable for stacktrace output
  set_prgname(argv[0]);

  if (Myrank ==0) {
    fprintf (stdout,"\n\n ******************************************************\n");
    fprintf (stdout,"            ____ ___ _____ _____ _\n");
    fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
    fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
    fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
    fprintf (stdout,"           |____/___|_|   |_____|_____| PARTRFEL");
    fflush(stdout);
    fprintf (stderr,"\n\n ******************************************************\n");
  }

  // initialize sequential global variables to null values
  initnull_globt();
  Ptp = NULL;
  Psolt = NULL;

  /* // stop point for attaching debuggers to particular processes
  if (Myrank==0){
    char bd;
    scanf ("%c",&bd);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  */

  if (process_argst(argc, argv, kwdsw))
    abort();

  //  names of input files
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
  in = xfopen (name,"r");
  if (in==NULL){
    par_print_err(Myrank+1, proc_name, "input file cannot be opened.",__FILE__, __LINE__, __func__);
    abort();
  }

  //  creation of main objects
  Ptp = new pprobdesct;
  Tp = new probdesct;
  Tp->kwdsw = kwdsw; // set keyword processing mode obtained from command line arguments

  Gtt = new gtopology;
  Tt = new transtop;
  Tm = new transmat;
  Tc = new transcrsec;
  Tb = new transbclc;
  Outdt = new outdrivert;
  
  // open log file
  char buf[1001];
  strcpy (buf,argv[1]);
  sprintf (&buf[strlen (argv[1])],"%d.log",Myrank+1);
  Outt = fopen (buf,"w");
  if (Outt==NULL){
    par_print_err (Myrank+1, proc_name,"log file cannot be opened",__FILE__, __LINE__, __func__);
  }

  // set new keyword processing mode for problem description section in the input file
  in->warning = 1;
  in->ignorecase = 1;
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_outd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  //  number of subdomain
  xfscanf (in,"%k%d","domain_num",&Ndom); // read domain number
  skipline(in);  // drop everything until the end of line
  Ndom--;

  fprintf (Outt,"Filename '%s'\n Nproc  %d\n Ndom   %d\n Myrank %d",name, Nproc, Ndom, Myrank);
  fflush(Outt);

  // **************************
  // **  basic problem data  **
  // **************************
  //  basic problem data
  Tp->read (in);

  //  creation of additional object
  Psolt = new psolver (Nproc, Myrank, proc_name, nameLength);
  //  reading data about parallel solver
  Psolt->read (in,Gtt,Ndom,Tp->tstorkm,Mesprt);
  
  // restore original keyword processing mode after reading of problem description section from  input file
  if (Tp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
 
  //  correspondence between processors and subdomains
  Psolt->procdomcorr ();

  
  // **************************************
  // **  input of node and element data  **
  // **************************************

  //  input of node and element data
  Tt->read (in);
  fflush(stdout);
  //  initiation of parallel solver
  Psolt->initiate (Gtt,argc,argv);
  //  reading of the array ltg (local-to-global correspondence)
  Psolt->read_ltg (in,Gtt,Outt);

  fprintf (Outt,"\n Nproc  %d",Nproc);
  fprintf (Outt,"\n Ndom   %d",Ndom);
  fprintf (Outt,"\n Myrank %d",Myrank);

  if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){
    Gtt->auto_subdom_detection (Outt);
    Tt->alloc_enodes ();
    Tt->alloc_edges ();
  }

  // *************************************
  //  generation of ordering of unknowns
  //  (generation of code numbers)
  //  number of degrees of freedom
  // *************************************
  Ndoft=Psolt->ordering (Gtt, Outt, proc_name);

  if (Mesprt==1)
    fprintf (stdout,"\n jouda number of degrees of freedom    %ld",Ndoft);

  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->read_gf (in);
  }

  //  input of material characteristics
  Tm->read (in);
  
  //int. points allocation and initiation
  Tm->intpointalloc ();
  Tm->intpointinit ();

  fprintf(Outt, "\nptrfelinit():\n Tm->ip = %p", (void*)Tm->ip);
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
  //Tt->node_materials ();
  //  toto je pokusne umisteni

  // ******************************************
  // **  definition of output from the code  **
  // ******************************************
  // set new keyword processing mode for data output section in the input file
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_pd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  Outdt->read(in);

  // restore original keyword processing mode after reading of data output section from the input file
  if (Tp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;

  //  names of backup files
  if (Tp->hdbcont.restore_stat())
  {
    strcpy (name,Tp->hdbcont.hdbnamer);
    dot = strchr (name,'.');
    if (dot == NULL)
    {
      print_err("cannot locate '.' in the backup file name '%s'", __FILE__, __LINE__, __func__, name);
      abort();
    }
    sprintf (dot,"%d",Myrank+1);
    sprintf (&name[strlen(name)], "%s", Tp->hdbcont.hdbnamer+int(dot-name));
    strcpy (Tp->hdbcont.hdbnamer, name);
  }
  if (Tp->hdbcont.save_stat())
  {
    strcpy (name,Tp->hdbcont.hdbnames);
    sprintf (&name[strlen(Tp->hdbcont.hdbnames)],"%d",Myrank+1);
    strcpy (Tp->hdbcont.hdbnames, name);
  }

  //  names of output files
  strcpy (name,Outdt->outfn);
  sprintf (&name[strlen (Outdt->outfn)],"%d.out",Myrank+1);
  strcpy (Outdt->outfn, name);
  if (Outdt->gf != grftt_no)
  {
    strcpy (name,Outdt->outgrfn);
    sprintf (&name[strlen (Outdt->outgrfn)],"%d",Myrank+1);
    strcpy (Outdt->outgrfn, name);
  }
  if (Outdt->ndiag > 0)
  {
    strcpy (name,Outdt->outdiagfn);
    sprintf (&name[strlen (Outdt->outdiagfn)],"%d.dat",Myrank+1);
    strcpy (Outdt->outdiagfn, name);
  }  

  //  nodes allocation
  Tt->alloc_nodes ();

  fflush(stdout);
  fflush(stderr);

  Ptp->shift_indices ();
  
  xfclose (in);
}
