#include "pmefelinit.h"
#include "pglobal.h"
#include "outdriverm.h"
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
void pmefel_init (int argc, const char *argv[])
{
  long merr;
  char name[1100], emsg[1000];
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
    fprintf (stdout,"           |____/___|_|   |_____|_____| PARMEFEL");
    fflush(stdout);
    fprintf (stderr,"\n\n ******************************************************\n");
  }

  // initialize global variables to null values
  initnull_glob();
  Pmp = NULL;
  Psolm = NULL;
  
  /* // stop point for attaching debuggers to particular processes
  if (Myrank==0){
    char bd;
    scanf ("%c",&bd);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  */

  if (process_args(argc, argv, kwdsw))
    abort();
  
  // generate name of input file and open it
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
  in = xfopen (name,"r");
  if (in==NULL){
    par_print_err(Myrank+1, proc_name, "input file cannot be opened.",__FILE__, __LINE__, __func__);
    abort();
  }
  
  //  creation of main objects
  Pmp = new pprobdesc;
  Mp  = new probdesc;
  Mp->kwdsw = kwdsw; // set keyword processing mode obtained from command line arguments

  Gtm = new gtopology;
  Mt  = new mechtop;
  Mm  = new mechmat;
  Mc  = new mechcrsec;
  Mb  = new mechbclc;  
  Outdm = new outdriverm;

  // open log file
  char buf[1001];
  strcpy (buf,argv[1]);
  sprintf (&buf[strlen (argv[1])],"%d.log",Myrank+1);
  Out = fopen (buf,"w");
  if (Out==NULL){
    par_print_err (Myrank+1,proc_name,"log file cannot be opened",__FILE__, __LINE__, __func__);
  }

  // set new keyword processing mode for problem description section in the input file
  in->warning = 1;
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_outd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;
  in->ignorecase = 1;

  //  number of subdomain (subdomain id)
  xfscanf (in,"%k%d","domain_num", &Ndom); // read domain number
  skipline(in);  // drop everything until the end of line
  Ndom--;

  // **************************
  // **  basic problem data  **
  // **************************
  //  reading data describing the problem
  Mp->read (in);
  //  creation of object of parallel solver
  Psolm = new psolver (Nproc,Myrank,proc_name,nameLength);  
  //  reading data about parallel solver
  Psolm->read (in,Gtm,Ndom,Mp->tstorsm,Mespr);
  
  // restore original keyword processing mode after reading of problem description section from  input file
  if (Mp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  
  //  correspondence between processors and subdomains
  Psolm->procdomcorr ();
  
  if (Myrank==0){
    fflush(stdout);
    fprintf (stdout,"\n\n kontrola domproc \n");
    fprintf (Out,"\n\n kontrola domproc \n");
    for (long i=0;i<Nproc;i++){
      fprintf (Out,"\n %ld  %ld",i,Psolm->domproc[i]);
      fprintf (stdout,"\n %ld  %ld",i,Psolm->domproc[i]);
    }
    fflush(stdout);
  }
  

  // **************************************
  // **  input of node and element data  **
  // **************************************
  
  //  input of node and element data
  Mt->read (in);
  Pmp->shift_indices ();
  
  // mesh control TKr
  if (Mp->tprob == growing_mech_structure){
    if (Mespr==1){
      merr=Mt->mesh_check();
      if (merr > 0){
        sprintf(emsg, "mesh error = %ld", merr);
	par_print_err(Myrank+1,proc_name,emsg, __FILE__, __LINE__, __func__);
	abort();
      }
    }
  }

  //  initiation of parallel solver
  Psolm->initiate (Gtm,argc,argv);  
  //  reading of the array ltg (local-to-global correspondence)
  Psolm->read_ltg (in,Gtm,Out);  
  
  /*
  // kontrola pro temelin (sdruzene stupne volnosti)
  long *nf = new long[Mt->nn];
  memset(nf, 0, sizeof(*nf)*Mt->nn);
  long i,j, ndofn, dof;
  for (i=0; i<Mt->nn; i++)
  {
    ndofn = Gtm->give_ndofn(i);
    if (Psolm->ltg[i]>=0)
      nf[i] = 2;

    for(j=0; j<ndofn; j++)
    {
      dof = Gtm->give_dof(i, j);
      if (dof > 1)
      {
        nf[i]++;
        break;
      }
    }
  }
  delete [] nf;
  */

  fprintf (Out,"\n Nproc  %d",Nproc);
  fprintf (Out,"\n Ndom   %d",Ndom);
  fprintf (Out,"\n Myrank %d",Myrank);
  fprintf (Out,"\n Name of processor %s",proc_name);
  

  // *************************************
  //  generation of ordering of unknowns
  //  (generation of code numbers)
  //  number of degrees of freedom
  // *************************************
  //  code numbers generation
  
  //Ndofm=Psolm->ptopjk->vse (Psolm->ltg,Gtm,Psolm->domproc,Out,proc_name);

  Ndofm=Psolm->ordering (Gtm,Out,proc_name);

  //  kontrolni tisk
  //fprintf (Out,"\n\n kontrola rohovych uzlu \n");
  //for (long i=0;i<Mt->nn;i++){
  //fprintf (Out,"\n %4ld %4ld",i,Psolm->ltg[i]);
  //}
  //fprintf (Out,"\n\n konec kontroly rohovych uzlu \n\n");

  // provizorni umisteni, po skonceni debuggovani vyhodit
  //Psolm->lpp->ndof=Ndofm;
    
  if (Mespr==1)  
    fprintf (stdout,"\n number of DOFs  %ld",Ndofm);  


  // ***********************************************
  // **  input of material models and parameters  **
  // ***********************************************

  //  reading of material models and parameters
  //  allocation of integation points
  //  definition of the number of material models (nm) on integration points
  //  definition of material types (tm) on integration points
  //  definition of number of material model (idm) on integration points
  Mm->read (in);
  //  definition of variables ncompstr on intergation points
  //  definition of variables ssst on intergation points
  //Mm->init_ip_1 ();


  // **********************************************
  // **  input of cross section characteristics  **
  // **********************************************
  Mc->read (in);
    

  // ************************
  // **  input of loading  **
  // ************************
  Mb->read (in);


  // ********************************************************
  // **  input of additional data necessary for problems   **
  // **  with changing number of nodes, elements and DOFs  **
  // ********************************************************
  if (Mp->tprob == growing_mech_structure){
    Gtm->read_gf (in);
    Mt->alloc_growstr ();
  }


  // *************************************************************
  //  vectors on left and right hand side of system of equations
  // *************************************************************
  Lsrs = new lhsrhs;
  //  allocation of arrays lhs, rhs (tdlhs, stdlhs if necessary)
  Lsrs->alloc ();
  //  function reads initial conditons for dynamics
  Lsrs->initcond ();

  if (Mp->tprob == mech_timedependent_prob)
  {
    Mb->alloc_sumcomp ();  // arrays for resultant force components
    if (Mm->csol)
      Gtm->comp_domain_sizes (); // maximum domain sizes used consolidation model
  }

  // ******************************************
  // **  definition of output from the code  **
  // ******************************************
  // set new keyword processing mode for data output section in the input file
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_pd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  Outdm->read(in);
  
  // restore original keyword processing mode after reading of data output section from the input file
  if (Mp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  

  // ************************************************
  // **  assembling of output files from the code  **
  // ************************************************

  //  names of backup files
  if (Mp->hdbcont.restore_stat())
  {
    strcpy (name,Mp->hdbcont.hdbnamer);
    dot = strchr (name,'.');
    if (dot == NULL)
    {
      print_err("cannot locate '.' in the backup file name '%s'", __FILE__, __LINE__, __func__, name);
      abort();
    }
    sprintf (dot,"%d",Myrank+1);
    sprintf (&name[strlen(name)], "%s", Mp->hdbcont.hdbnamer+int(dot-name));
    strcpy (Mp->hdbcont.hdbnamer, name);
  }
  if (Mp->hdbcont.save_stat())
  {
    strcpy (name,Mp->hdbcont.hdbnames);
    sprintf (&name[strlen(Mp->hdbcont.hdbnames)],"%d",Myrank+1);
    strcpy (Mp->hdbcont.hdbnames, name);
  }

  //  names of output files
  strcpy (name,Outdm->outfn);
  sprintf (&name[strlen (Outdm->outfn)],"%d.out",Myrank+1);
  strcpy (Outdm->outfn, name);
  if (Outdm->gf != grfmt_no)
  {
    strcpy (name,Outdm->outgrfn);
    sprintf (&name[strlen (Outdm->outgrfn)],"%d",Myrank+1);
    strcpy (Outdm->outgrfn, name);
  }
  if (Outdm->ndiag > 0)
  {
    strcpy (name,Outdm->outdiagfn);
    sprintf (&name[strlen (Outdm->outdiagfn)],"%d.dat",Myrank+1);
    strcpy (Outdm->outdiagfn, name);
  }  
  
  //fprintf (stdout,"\n procesor %ld    vystup %s    grafika %s   diagram %s",Myrank,Outdm->outfn,Outdm->outgrfn,Outdm->outdiagfn);
  
  //print_init(-1, "wt");    
  //Outdm->outf = fopen(Outdm->outfn, "w");
  
  //  patek
  //Psolm->nodesplit (Gtm,Outdm->outf);
  //Psolm->nodesplit (Gtm,Out);
  
  //  list of adjacent elements (NUTNO PREDELAT PRO PARALELNI ULOHU)
  Gtm->adjacnodes (Out);
  Gtm->adjacelem (Out);

  //  zacatek JK
  //  toto jsou uzivatelsky definovane body
  //Mm->stra.init(strain);  
  //Mm->stre.init(stress);  
  //  konec JK


  //  zacatek TKo
  //  definition of nonlocal models and models for temperature
  Mm->init_ip_2 ();
  //  konec TKo
  
  //  allocation of arrays defined on nodes
  Mt->alloc_nodes_arrays ();

  //  number of components in arrays used in outdriver
  //  maximum number of components in array strain/stress on nodes and elements
  Mt->give_maxncompstr(Mm->max_ncompstrn, Mm->max_ncompstre);
  //  maximum number of components in array other on nodes and elements
  Mt->give_maxncompo(Mm->max_ncompothern, Mm->max_ncompothere);

  //  zacatek JK
  //if (Mp->straincomp==1){
  //Mm->stra.alloc(Mb->nlc);
  //}
  //if (Mp->stresscomp==1){
  //Mm->stre.alloc(Mb->nlc);
  //}
  //  konec JK

  //  prescribed displacements on elements investigation
  Mt->elemprescdisp ();
  //  prescribed temperatures on elements investigation
  if (Mp->tprob == linear_statics){
    Mt->elempresctemp (0);
  }
  //  reactions on elements investigation
  if (Mp->reactcomp==1){
    Mt->comreac ();
  }

  //  computation of initial values at integration points from initial nodal values
  if (Mb->nico > 0)
    Mb->inicipval();

  // NUTNO PREDELAT PRO PARALELNI ULOHU
  if (Mp->matmodel == nonlocal){    
    //  list of adjacent elements
    Gtm->adjacelem (Out);
    //  list of adjacent integration points was not restored 
    //  from file it is necessary to create it from scratch
    if (restore_adjacip() == 0)
      adjacip ();
    save_adjacip();
    ipvolume ();
  }

  // adaptivity
  if (Mp->adaptivityflag)
    Ada = new adaptivity ();
  else
    Ada = NULL;

  //  nondeterministic computation
  if (Mp->stochasticcalc){
    St = new stochdriver;
    St->read (in);
  }
  else
    St = NULL;

  //  initiation and computation of eigenstrains
  if (Mp->temperature==1){
    Mm->alloceigstrains ();
    Mm->alloctempstrains ();
    Mm->alloceigstresses ();
    
    //Mm->tempr = new double [Mm->tnip];
    //Mm->inittempr = new double [Mm->tnip];
    
    //for (i=0;i<Mm->tnip;i++){
    //Mm->tempr[i]=0.0;
    //Mm->inittempr[i]=0.0;
    //}    
  }

  //Mt->nodedisplacements ();  

  fflush(stdout);
  fflush(stderr);
  
  Pmp->shift_indices ();
  
  xfclose (in);
}
