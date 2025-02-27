#include "mefelinit.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "lhsrhs.h"
#include "outdriverm.h"
#include "stochdriver.h"
#include "globalg.h"
#include "globmat.h"
#include "sequent.h"
#include "intpoints.h"
#include "gtopology.h"
#include "meshtransfer.h"
#include "adaptivity.h"
#include "mechprint.h"
#include "nssolver.h"
#include "element.h"
#include "elemswitch.h"
#include "iotools.h"
#include "xfile.h"
#include "saddpoint.h"
#include "stacktrace.h"
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>




/**
  The function initializes mefel from the input file.
  It reads data from input file and initiates all necessary variables and arrays.

  @param argc[in] - number of passed command line parameters
  @param argv[in] - array of strings with command line arguments passed to program 

  Returns:
  @retval The function does not return anything but changes state of global variables.

  Created by JK,
  Modified by Tomas Koudelka,
*/
void mefel_init (int argc, const char *argv[])
{
  XFILE *in;
  pkwd_sw kwdsw;

  // set name of executable for stacktrace output
  set_prgname(argv[0]);
  
  fprintf (stdout,"\n\n ******************************************************\n");
  fprintf (stdout,"            ____ ___ _____ _____ _\n");
  fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
  fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
  fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
  fprintf (stdout,"           |____/___|_|   |_____|_____| MEFEL");
  fflush(stdout);
  fprintf (stderr,"\n\n ******************************************************\n");
  
  // initialize global variables to null values
  initnull_glob();

  if (process_args(argc, argv, kwdsw))
    abort();

  in = xfopen(argv[1],"r");
  if (in==NULL){
    print_err("input file cannot be opened.",__FILE__,__LINE__,__func__);
    abort();
  }

  Mp  = new probdesc;
  Mp->kwdsw = kwdsw; // set keyword processing mode obtained from command line arguments
  
  Gtm = new gtopology;
  Mt  = new mechtop;
  Mm  = new mechmat;
  Mc  = new mechcrsec;
  Mb  = new mechbclc;
  Outdm = new outdriverm;

  Out  = fopen ("mefel.log","wt");

  in->warning = 1;
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_outd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  
  // **************************
  // **  basic problem data  **
  // **************************
  Mp->read (in);

  if (Mp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;

  // **************************************
  // **  input of node and element data  **
  // **************************************
  Mt->read (in);

  // mesh control TKr
  /*
  if (Mp->tprob == growing_mech_structure){
    if (Mespr==1){
      merr=Mt->mesh_check();
      if (merr > 0){
        sprintf(emsg, "mesh error = %ld", merr);
        print_err(emsg, __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
  */

  // zacatek TKo
  if (Mp->tprob == earth_pressure)
    {
      Gtm->bckcn = new long*[Gtm->nn];
      memset(Gtm->bckcn, 0, sizeof(*Gtm->bckcn)*Gtm->nn);
      Gtm->backup_codnum();
    }
  //  konec TKo
  
  
  if (Mp->tprob==hemivar_inequal || Mp->tprob==lin_floating_subdomain){
    Gtm->auto_subdom_detection (Out);
    Mt->alloc_enodes ();
    Mt->alloc_edges ();
  }
  
  // *************************************
  //  generation of ordering of unknowns
  //  (generation of code numbers)
  //  number of degrees of freedom
  // *************************************
  Ndofm = Gtm->codenum_generation (Out);

  
  if (Mespr==1)
    fprintf (stdout,"\n number of degrees of freedom    %ld",Ndofm);
  //  zacatek TKr
  //Mt->store_code_num_elem(); // enable combination of plane elements, springs and bar elems with beams//tady??!!!
  //  konec TKr
  
  //  zacatek JK
  if (Mp->tprob==layered_linear_statics){
    //  definition of type of layered problem
    if (Mt->elements[0].te<20)  Mp->tlm=1;
    if (Mt->elements[0].te>20)  Mp->tlm=2;
    
    //  generation of ordering of Lagrange multipliers
    Mt->gencodnumlagrmult (Ndofm);
    Gtm->initiate_elemmult ();
    if (Mespr==1)
      fprintf (stdout,"\n total number of degrees of freedom    %ld",Ndofm);
  }


  if (Mp->tprob==lin_floating_subdomain){
    //  auxiliary informations about ordering of unknowns
    Gtm->flsub.ndof_subdom (Gtm->gnodes);
  }  
  //  konec JK  
    

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
  
  if ((Mp->homog==3) || (Mp->homog==5) || (Mp->homog==7) || (Mp->homog==9)){
    //  the first order homogenization based on prescribed macro-stresses
    //
    // The total number of macro-stress/strain components is determined by
    // Mt->max_ncompstr which is initialized in Mt->read function
    
    // Actually, macro-sress component code numbers are supposed to be
    // the same in all load cases => take them from the first load case
    long *mstress_cn = Mb->give_mstress_cn(0);
    ivector cnadd;
    //  makerefv(cnadd, mstress_cn, Mt->max_ncompstr);
    makerefv(cnadd, mstress_cn, Mt->max_ncompstr);

    //  code numbers are copied from nodes to elements
    //  and additional DOFs are added to elements
    // do it only in the case that there are some macro-stress components prescribed (due to Mp->homog==9)
    if (Mb->give_num_mstress_comp(0)){
      for (long i=0; i<Mt->ne; i++)
        Gtm->gelements[i].ndofe += Mt->give_tncomp(i);
      
      Gtm->code_num_mod (cnadd, Ndofm);
    }
    
    //  volume of the domain solved
    Mt->domvol = Mt->give_domain_vol();
    fprintf (stdout,"\n domain volume  %le",Mt->domvol);
  }
  if ((Mp->homog==4) || (Mp->homog==6) || (Mp->homog==8) || (Mp->homog==9)){
    //  the first order homogenization based on prescribed macro-strains
    //
    // The total number of macro-stress/strain components is determined by
    // Mt->max_ncompstr which is initialized in Mt->read function

    //  volume of the domain solved
    if (Mp->homog != 9){
      Mt->domvol = Mt->give_domain_vol();
      fprintf (stdout,"\n domain volume  %le",Mt->domvol);
    }
  }

  // ********************************************************
  // **  input of additional data necessary for problems   **
  // **  with changing number of nodes, elements and DOFs  **
  // ********************************************************
  if (Mp->tprob == growing_mech_structure){
    Gtm->read_gf (in);
    Gtm->alloc_growstr();
    Mt->alloc_growstr();
  }
  
  //  some solvers require additional data
  //  if the Schur complement method or FETI method
  //  are used, the variable Ndofm is rewritten
  //  Ndofm is not changed in all other cases
  Mp->ssle->prepare_data (&Ndofm,Gtm,Out);

  


  // *************************************************************
  //  vectors on left and right hand side of system of equations
  // *************************************************************
  //
  // Ndofm MUST be known here for the correct allocation of Lsrs
  //
  Lsrs = new lhsrhs;
  //  allocation of arrays lhs, rhs (tdlhs, stdlhs if necessary)
  Lsrs->alloc ();
  //  function reads initial conditons for dynamics
  Lsrs->initcond ();

  if ((Mp->tprob == mech_timedependent_prob)||(Mp->tprob == growing_mech_structure))
  {
    Mb->alloc_sumcomp ();  // arrays for resultant force components
    if (Mm->csol)
      Gtm->comp_domain_sizes (); // maximum domain sizes used consolidation model
  }

  // ******************************************
  // **  definition of output from the code  **
  // ******************************************
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_pd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  Outdm->read(in);

  if (Mp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  
  // ************************************************
  // **  assembling of output files from the code  **
  // ************************************************
  
  fclose (Out);
  char fname[1020];
  if (*Outdm->outfn)
    filename_decomposition (Outdm->outfn,Mp->path,Mp->filename,Mp->suffix);
  else
    filename_decomposition (argv[1],Mp->path,Mp->filename,Mp->suffix);
   
  sprintf(fname, "%s%s.log", Mp->path, Mp->filename);
  

  Out = fopen(fname,"wt");
  if (Out == NULL)
    print_err("cannot open log file", __FILE__, __LINE__, __func__);

  //  list of adjacent elements
  //  Gtm->adjacnodes (Out);
  //  Gtm->adjacelem (Out);
  Gtm->adjacnodes (NULL);
  Gtm->adjacelem (NULL);
  
  
  //  zacatek JK
  //  toto jsou uzivatelsky definovane body
  //Mm->stra.init(strain);  
  //Mm->stre.init(stress);  
  //  konec JK
  
  //  zacatek TKo
  // detection of nonlocal models and models for temperature,
  // allocation of stress, strain, eqother, other and nonloc arrays,
  // setting ncompeqother and ncompother
  Mm->init_ip_2 ();
  //  konec TKo
  
  
  //  allocation of arrays defined on nodes
  Mt->alloc_nodes_arrays ();
  
  //  number of components in arrays used in outdriver
  //  maximum number of components in array strain/stress on nodes and elements
  Mt->give_maxncompstr(Mm->max_ncompstrn, Mm->max_ncompstre);
  //  maximum number of components in array other on nodes and elements
  Mt->give_maxncompo(Mm->max_ncompothern, Mm->max_ncompothere);
  //  konec TKo
  
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
  
  
  
  if (Mp->matmodel == nonlocal){    
    //  list of adjacent elements -> it has been already assembled in the above code
    //  Gtm->adjacelem (Out);
    //  list of adjacent integration points was not restored 
    //  from file it is necessary to create it from scratch
    if (restore_adjacip() == 0){
      adjacip2 ();
      /*
      FILE *fadjip = fopen("adjip.log", "wt");
      print_adjacip(fadjip);
      fclose(fadjip);*/
    }
    save_adjacip();
    ipvolume ();
  }
  
  if (Mb->give_num_mstrain_comp(0))
  {
    if (Mp->matmodel != nonlocal) // int. point volumes have been already computed for non-local models
      ipvolume ();
    Mm->allocmacrostresses();
  }
  
  // adaptivity
  //if (Mp->adaptivity)
  //  Ada = new adaptivity ();
  
  //  nondeterministic computation
  if (Mp->stochasticcalc){
    St = new stochdriver;
    St->read (in);
  }
  else
    St = NULL;
  
  //  initiation of eigenstrains
  if (Mp->eigstrains!=0){
    Mm->alloceigstrains ();
    Mm->alloceigstresses ();
  }
  
  if (Mp->pore_press!=0){
    Mm->alloceigstresses ();
  }

  //  initiation of temprstrains
  if (Mp->temperature > 0)
  {
    Mm->alloctempstrains (); 
    Mm->alloceigstresses ();
  }
  
  //Mt->nodedisplacements ();
  
  

  xfclose (in);

  #ifdef INC_OPENMP
  Gtm->initiate_elem_dof_ranges();
  Gtm->initiate_omp_elem_order();
  for (long i = 0; i<Mt->ne; i++){
    // fprintf(Out, "Elem %6ld, mindof %ld, deltadof %ld\n",
    //         Gtm->ompelord[i]+1, Gtm->eldofr[Gtm->ompelord[i]].mindof,
    //         Gtm->eldofr[Gtm->ompelord[i]].deltadof);
    fprintf(Out, "%6ld %6ld\n",  Gtm->ompelord[i]+1, i+1);
  }
  #endif
  Omp_wtime = 0.0;
}



/**
  Function prints help for command line to the opened file out

  @param out[in] - pointer to the opened text file for output of help
  @param prgname[in] - program name

  created 15.10.2007 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz
*/
void print_help(FILE *out, const char *prgname)
{
  fflush(out);
  fprintf(out, "\nInput file has not been specified\n");
  fprintf(out, " use: %s infile [-kwd=i] [-numth=%%d]\n",prgname);
  fprintf(out, " where:\n *   infile is input file name\n");
  fprintf(out, " *   -kwd={0,1,2,3,4} is optional switch which controls keyword usage in the input file\n");
  fprintf(out, "           0 means no keywords (default value)\n");
  fprintf(out, "           1 means keywords only in the problem description section\n");
  fprintf(out, "           2 means keywords only in the outdriver section\n");
  fprintf(out, "           3 means keywords in the problem description and outdriver sections\n");
  fprintf(out, "           4 means full usage of keywords\n\n");
  fprintf(out, " *   -numth={%%d} is optional switch which specifies number of threads used in OpenMP\n");
  abort();
}



/**
  Function processes command line switches.

  @param argc[in] - the total number of command line arguments, i.e. dimension of array argv
  @param argv[in] - pointer to string array with particularcommand line arguments

  @retval 0 - on success
  @retval 1 - error in processing of arguments

  created 15.10.2007 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz
 */
long process_args(int argc, const char *argv[], pkwd_sw &kwdsw)
{
  const char *ptrkwd;
 #if defined(INC_OPENMP) || defined(INC_CHMD)
  // set default number of threads
  Numth = 8;
 #endif
  switch (argc)
  {
    case 0:
      fflush(stderr);
      fflush(stdout);
      print_help(stderr,"pmefel|mefel");
      return 1;
    case 1:
      fflush(stderr);
      fflush(stdout);
      print_help(stderr,argv[0]);
      return 1;
    case 2:
      kwdsw = nokwd;
      break;
    case 3:
    case 4:  
      // processing of optional keyword switch -kwd={0,1,2,3,4}
      // and -numth=%d keyword
      for (long i=2; i<argc; i++){
        ptrkwd = strstrcis(argv[i], "-kwd=");
        if (ptrkwd){
          kwdsw = pkwd_sw(-1);
          sscanf(ptrkwd+5, "%d", (int *)&kwdsw);
          if ((kwdsw < 0) || (kwdsw > 4)){
            fflush(stdout);
            fprintf(stderr, "\n\nError - unknown -kwd value %d is obtained\n\n", kwdsw);
            print_help(stderr,argv[0]);
            return 1;
          }
        continue;
        }
        kwdsw = nokwd;
    #if defined(INC_OPENMP) || defined(INC_CHMD)
        ptrkwd = strstrcis(argv[i], "-numth=");
        if (ptrkwd){
          Numth = -1;
          sscanf(ptrkwd+7, "%d", &Numth);
          if (Numth < 1){
            fflush(stdout);
            fprintf(stderr, "\n\nError - wrong number of threads (%d), it must be > 0\n\n", Numth);
            print_help(stderr,argv[0]);
            return 1;
          }
          continue;
        }
    #else
        ptrkwd = strstrcis(argv[i], "-numth=");
        if (ptrkwd){
          fflush(stdout);
          fprintf(stderr, "\n\nWarning - argument with number of threads is ignored because OpenMP\n"
                  " was not enabled in this executable\n\n");
          continue;
        }
    #endif
      }
      break;
    default:
      fflush(stderr);
      fflush(stdout);
      fprintf(stderr, "\n\nError: too many arguments is used\n");
      print_help(stderr,argv[0]);
      return 1;
  }
 #if defined(INC_OPENMP) || defined(INC_CHMD)
  fprintf(stdout, "\n Number of threads used in MEFEL: %d\n", Numth);
 #endif
  return 0;
}
