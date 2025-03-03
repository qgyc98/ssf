#include "mpi.h"
#include "pmetrinit.h"
#include "metrinit.h"
#include "mefelinit.h"
#include "trfelinit.h"
#include "pglobalc.h"
#include "outdriverm.h"
#include "pglobalt.h"
#include "pglobal.h"
#include "xfile.h"
#include "globalc.h"
#include "globalt.h"
#include "globalc.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "stacktrace.h"
#include "globmatc.h"

void pmetr_init (int argc, const char *argv[])
{
  long i, merr, tmp,n;
  char name[1100];
  XFILE *inm,*intr,*inc;
  pkwd_sw kwdsw;

  // set name of executable for stacktrace output
  set_prgname(argv[0]);
  
  if (Myrank ==0) {
    fprintf (stdout,"\n\n ******************************************************\n");
    fprintf (stdout,"            ____ ___ _____ _____ _\n");
    fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
    fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
    fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
    fprintf (stdout,"           |____/___|_|   |_____|_____| PARMETR");
    fflush(stdout);
    fprintf (stderr,"\n\n ******************************************************\n");
  }

  // initialize parallel global METR variables to null values
  Pcp = NULL;
  // initialize sequential global METR variables to null values
  initnull_globc();
  // initialize sequential global MEFEL variables to null values
  initnull_glob();
  // initialize sequential global TRFEL variables to null values
  initnull_globt();
  fprintf(stderr, "\n\n Myrank %d - Nproc = %d\n", Myrank+1, Nproc);

  // stop point for attaching debuggers to particular processes
  /*
  if (Myrank==0){
    char bd;
    scanf ("%c",&bd);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  */

  if (process_argsc(argc, argv, kwdsw))
    abort();

  //  names of input files
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
  inc = xfopen (name,"r");
  if (inc==NULL){
    par_print_err(Myrank+1, proc_name, "input file cannot be opened.",__FILE__, __LINE__, __func__);
    abort();
  }


  //  creation of main objects
  Pcp = new pprobdescc;
  //  creation of object of parallel solver for mechanical problem
  Psolm = new psolver (Nproc, Myrank, proc_name, nameLength);
  //  creation of object of parallel solver for transport problem
  Psolt = new psolver (Nproc, Myrank, proc_name, nameLength);
  

  //  problem description of all coupled problem
  Cp = new probdescc;
  //  problem description of mechanical part of problem
  Mp = new probdesc;
  Pmp = new pprobdesc;
  //  problem description of thermo-hydro part of problem
  Tp = new probdesct;
  Ptp = new pprobdesct;
  //  topology of mechanical part of problem
  Gtm = new gtopology;
  //  topology of thermo-hydro part of problem
  Gtt = new gtopology;
  //  auxiliary united topology (mechanical + thermo-hydro)
  Gtu = new gtopology;
  //  additional topological data about mechanical part of problem
  Mt = new mechtop;
  //  additional topological data about thermo-hydro part of problem
  Tt = new transtop;
  //  additional topological data about thermo-hydro part of problem
  Ct = new couptop;
  //  materials of mechanical part of problem
  Mm = new mechmat;
  //  materials of thermo-hydro part of problem
  Tm = new transmat;
  //  materials of thermo-hydro part of problem
  Cmu = new coupmatu;
  //  materials of thermo-hydro part of problem
  Cml = new coupmatl;
  //  cross-section in mechanical part of problem
  Mc = new mechcrsec;
  //  cross-section in thermo-hydro part of problem
  Tc = new transcrsec;
  //  boundary conditions and load cases in coupled problem
  Cb = new coupbclc;
  //  boundary conditions and load cases in mechanical part of problem
  Mb = new mechbclc;
  //  boundary conditions and load cases in thermo-hydro part of problem
  Tb = new transbclc;
  //  mechanical outdriver
  Outdm = new outdriverm;
  //  transport outdriver
  Outdt = new outdrivert;
  //  coupled problem outdriver
  Outdc = new outdriverc;
  
  //  docasne umisteni
  if (Bar2d==NULL)    Bar2d = new barel2d;
  if (Barq2d==NULL)   Barq2d = new barelq2d;
  if (Pelq==NULL)     Pelq = new planeelemlq;
  if (Peqq==NULL)     Peqq = new planeelemqq;
  if (Asymlq==NULL)   Asymlq = new axisymlq;
  if (Lhex==NULL)     Lhex = new linhex;
  if (Ltet==NULL)     Ltet = new lintet;
  if (Bar3d==NULL)    Bar3d = new barel3d;
  if (Barq3d==NULL)   Barq3d = new barelq3d;

  //if (Myrank==0){
  //char bd;
  //scanf ("%c",&bd);
  //}
  //MPI_Barrier (MPI_COMM_WORLD);
  


  inc->warning = 1;
  inc->ignorecase = 1;
  if ((Cp->kwdsw == nokwd) || (Cp->kwdsw == kwd_outd))
    inc->kwdmode = ignore_kwd;
  else
    inc->kwdmode = sequent_mode;
  
  //  number of subdomain
  xfscanf (inc,"%d",&Ndom);
  Ndom--;
  
  //  problem description of all coupled problem
  Cp->read (inc);

  if (Cp->tprob!=fully_coupled_material){
    //  input file for mechanical part
    strcpy (name,Cp->minfile);
    sprintf (&name[strlen (Cp->minfile)],"%d.in",Myrank+1);
    inm = xfopen(name,"r");
    if (inm==NULL){
      par_print_err(Myrank+1, proc_name, "input file for mechanical part cannot be opened.",__FILE__, __LINE__, __func__);
      abort ();
    }  
    //  input file for thermo-hydro part
    strcpy (name,Cp->tinfile);
    sprintf (&name[strlen (Cp->tinfile)],"%d.in",Myrank+1);
    intr = xfopen(name,"r");
    if (intr==NULL){
      par_print_err(Myrank+1, proc_name, "input file for thermo-hydro part cannot be opened.",__FILE__, __LINE__, __func__);
      abort ();
    }
  }
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.log",Myrank+1);
  //sprintf (&name[strlen (argv)],"%d.log",Myrank+1);
  Outc = fopen (name,"w");
  if (Outc==NULL){
    par_print_err(Myrank+1, proc_name, "log file for coupled part cannot be opened.",__FILE__, __LINE__, __func__);
    abort ();
  }
  strcpy (name,Cp->minfile);
  sprintf (&name[strlen (Cp->minfile)],"%d.log",Myrank+1);
  Out = fopen (name,"w");
  //Out = fopen ("mefel.log","w");
  if (Out==NULL){
    par_print_err(Myrank+1, proc_name, "log file for mechanical part cannot be opened.",__FILE__, __LINE__, __func__);
    abort ();
  }
  strcpy (name,Cp->tinfile);
  sprintf (&name[strlen (Cp->tinfile)],"%d.log",Myrank+1);
  Outt = fopen (name,"w");
  //Outt = fopen ("trfel.log","w");
  if (Outt==NULL){
    par_print_err(Myrank+1, proc_name, "log file for thermo-hydro part cannot be opened.",__FILE__, __LINE__, __func__);
    abort ();
  }

  inm->warning = 1;
  inm->ignorecase = 1;
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_outd))
    inm->kwdmode = ignore_kwd;
  else
    inm->kwdmode = sequent_mode;

  intr->warning = 1;
  intr->ignorecase = 1;
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_outd))
    intr->kwdmode = ignore_kwd;
  else
    intr->kwdmode = sequent_mode;

  //  number of subdomain can be ommitted
  xfscanf (inm,"%ld", &tmp);
  //  reading data describing the problem
  Mp->read (inm);
  
  //  number of subdomain can be ommitted
  xfscanf (intr,"%ld", &tmp);
  //  basic problem data
  Tp->read (intr);
  
  //  reading data about parallel solver
  Psolm->read (inm,Gtm,Ndom,Mp->tstorsm,Mespr);
  
  //  reading data about parallel solver
  Psolt->read (intr,Gtt,Ndom,Tp->tstorkm,Mesprt);
  
  if ((Cp->tprob==growing_par_coupl_mech_trans) || (Cp->tprob==par_coupl_mech_trans)){
    //  initial time setting
    Cp->time=Cp->timecon.starttime ();
    Tp->time=Tp->timecont.starttime ();
    Mp->time=Mp->timecon.starttime ();
  }

  if (Mp->kwdsw == kwd_full)
    inm->kwdmode = sequent_mode;
  else
    inm->kwdmode = ignore_kwd;
  
  if (Tp->kwdsw == kwd_full)
    intr->kwdmode = sequent_mode;
  else
    intr->kwdmode = ignore_kwd;
  
  if (Cp->kwdsw == kwd_full)
    inc->kwdmode = sequent_mode;
  else
    inc->kwdmode = ignore_kwd;
  

  Psolm->procdomcorr ();
  Psolt->procdomcorr ();
  

  //  input of node and element data
  Mt->read (inm);
  Tt->read (intr);
  Ct->read (inc);

  // mefel mesh control TKr
  if (Mp->tprob == growing_mech_structure){
    if (Mespr==1){
      merr=Mt->mesh_check();
      if (merr > 0){
	fprintf (stdout,"MEFEL Mesh error = %ld\n\n\n",merr);
	exit(0);
      }
    }
  }

    
  // trfel mesh control
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    if (Mesprt==1){
      merr=Tt->mesh_check();
      if (merr > 0){
	fprintf (stdout,"TRFEL Mesh error = %ld\n\n\n",merr);
	exit(0);
      }
    }
  }

  // metr mesh control
  /*
  if ((Cp->tprob==growing_par_coupl_mech_trans) || (Cp->tprob==par_coupl_mech_trans)){
    if (Mesprc==1){
      merr=Ct->mesh_check();
      if (merr > 0){
	fprintf (stdout,"METR Mesh error = %ld\n\n\n",merr);
	exit(0);
      }
    }
  }
  */

  if (Cp->tprob==fully_coupled_mech_trans){
    //  unification of two topologies
    Gtu->comptop (Gtm,Gtt);
  }
  
  //  initiation of the mechanical solver
  Psolm->initiate (Gtm,argc,argv);
  //  initiation of the transport solver
  Psolt->initiate (Gtt,argc,argv);
  

  //  reading of the array ltg (local-to-global correspondence)
  Psolm->read_ltg (inm,Gtm,Out);
  Psolt->read_ltg (intr,Gtt,Outt);
  

  //  generation of code numbers for mechanical problem
  Ndofm=Psolm->ordering (Gtm, Out, proc_name);
  if (Mesprc==1)
    fprintf (stdout,"\n number of degrees of freedom in mechanics    %ld",Ndofm);
  //  generation of code numbers for transport problem
  Ndoft=Psolt->ordering (Gtt, Outt, proc_name);
  if (Mesprc==1)
    fprintf (stdout,"\n number of degrees of freedom in transports    %ld",Ndoft);
  if (Cp->tprob==fully_coupled_mech_trans){
    //  generation of code numbers in both parts
    Ndofc=Gtu->gencodnum ();
    //  distribution of code numbers to local topologies
    Gtu->cndistr (Gtm,Gtt);
    Ndofm=Ndofc;
    Ndoft=Ndofc;
    if (Mesprc==1)
      fprintf (stdout,"\n total number of degrees of freedom    %ld",Ndofc);
  }

  //read trfel function
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->read_gf (intr);
  }

  //  input of material characteristics
  //  input of material characteristics
  if (Cp->tprob!=fully_coupled_material){
    Mm->read (inm);
    Tm->read (intr);
  }
  if (Cp->tprob==fully_coupled_material){
    Cm->read (inc);
  }
  
  //int. points allocation and initiation
  Tm->assemble_dof_nameord(Tp->dofname, Tp->ntm, Tp->var_dofid, Tp->tnkv);

  
  if (Cp->tprob==fully_coupled_mech_trans){
    Cmu->read (inc);
    Cml->read (inc);
  }  
  
  if (Cp->tprob!=fully_coupled_material){
    //  input of cross section characteristics
    Mc->read (inm);
    Tc->read (intr);

    //  input of loading
    Mb->read (inm);
    Tb->read (intr);
  }else{
  }
  

  
  //read mefel function
  if (Mp->tprob == growing_mech_structure){
    Gtm->read_gf (inm);
    Gtm->alloc_growstr();
    Mt->alloc_growstr();
  }

  //  prescribed values on elements investigation
  //Tt->elemprescvalue();
  //  searching for elements with nodes where source is defined
  Tb->elemsource ();

  // allocation of non-mechanical and non-transport arrays
  // definition of their ordering
  n = Mm->search_reqnmq(Mm->nmqo);
  if (n){
    for (i=0; i<n; i++)
      Mm->nmqid[Mm->nmqo[i]-1] = i;
    
    Mm->alloc_nonmechq(n);
  }
  
  n = Tm->search_reqntq(Tm->ntqo);
  if (n){
    for (i=0; i<n; i++)
      Tm->ntqid[Tm->ntqo[i]-1] = i;
    
    Tm->alloc_nontransq(n);
  }

  if (Cp->tprob==fully_coupled_mech_trans){  
    if (Mp->pore_press)   // effective stress concept was defined in the MEFEL
      Mp->pore_press = 3; // indicator of fully coupled approach

    //  left and right hand sides
    Lsrsc = new lhsrhsc;
    Lsrsc->alloc ();
    //Lsrsc->initcond (inc);
    
    Lsrs = new lhsrhs;
    Lsrs->lhsi=Lsrsc->lhsi; //initial conditions
    Lsrs->initcond ();
    
    Lsrst = new lhsrhst;     
    Lsrst->lhsi=Lsrsc->lhsi; //initial conditions
    Lsrst->initcond (intr);    
  }
  
  if ((Cp->tprob==par_coupl_mech_trans) || (Cp->tprob==growing_par_coupl_mech_trans)){
    if (Mp->pore_press)   // effective stress concept was defined in the MEFEL
      Mp->pore_press = 2; // indicator of partially coupled approach

    //  object for left and right hand side of MEFEL problem
    Lsrs = new lhsrhs;
    Lsrs->alloc ();
    Lsrs->initcond ();
    
    //  object for left and right hand side of TRFEL problem
    Lsrst = new lhsrhst;
    Lsrst->alloc ();
    Lsrst->initcond (intr);

  }
  
  Mb->alloc_sumcomp ();  // arrays for resultant force components
  if (Mm->csol)
    Gtm->comp_domain_sizes (); // maximum domain sizes used consolidation model

  // input of output description
  //metr
  if ((Cp->kwdsw == nokwd) || (Cp->kwdsw == kwd_pd))
    inc->kwdmode = ignore_kwd;
  else
    inc->kwdmode = sequent_mode;

  Outdc->read(inc);

  //mefel
  if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_pd))
    inm->kwdmode = ignore_kwd;
  else
    inm->kwdmode = sequent_mode;

  Outdm->read(inm);

  //trfel
  if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_pd))
    intr->kwdmode = ignore_kwd;
  else
    intr->kwdmode = sequent_mode;

  Outdt->read(intr);
  

  // backup file names for MEFEL
  if (Mp->hdbcont.restore_stat()) // restorage is required
  {
    strcpy (name,Mp->hdbcont.hdbnamer);
    sprintf (&name[strlen (Mp->hdbcont.hdbnamer)],"%d.out",Myrank+1);
    strcpy (Mp->hdbcont.hdbnamer, name);
  }
  if (Mp->hdbcont.save_stat()) // saving is required
  {
    strcpy (name,Mp->hdbcont.hdbnames);
    sprintf (&name[strlen (Mp->hdbcont.hdbnames)],"%d.out",Myrank+1);
    strcpy (Mp->hdbcont.hdbnames, name);
  }

  //  names of output MEFEL files
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

  // backup file names for TRFEL
  if (Tp->hdbcont.restore_stat()) // restorage is required
  {
    strcpy (name,Tp->hdbcont.hdbnamer);
    sprintf (&name[strlen (Tp->hdbcont.hdbnamer)],"%d.out",Myrank+1);
    strcpy (Tp->hdbcont.hdbnamer, name);
  }
  if (Tp->hdbcont.save_stat()) // saving is required
  {
    strcpy (name,Tp->hdbcont.hdbnames);
    sprintf (&name[strlen (Tp->hdbcont.hdbnames)],"%d.out",Myrank+1);
    strcpy (Tp->hdbcont.hdbnames, name);
  }

  //  names of output TRFEL files
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


  //  zacatek TKo
  //  definition of nonlocal models and models for temperature
  Mm->init_ip_2 ();
  //  konec TKo

  //  nodes allocation
  Tt->alloc_nodes ();
  //  allocation of arrays defined on nodes
  Mt->alloc_nodes_arrays ();
  //  number of components in arrays used in outdriver
  //  maximum number of components in array strain/stress on nodes and elements
  Mt->give_maxncompstr(Mm->max_ncompstrn, Mm->max_ncompstre);
  //  maximum number of components in array other on nodes and elements
  Mt->give_maxncompo(Mm->max_ncompothern, Mm->max_ncompothere);
  //  konec TKo


  //  prescribed displacements on elements investigation
  Mt->elemprescdisp ();
  //  reactions on elements investigation
  if (Mp->reactcomp==1){
    Mt->comreac ();
  }

  //  computation of initial values at integration points from initial nodal values
  if (Mb->nico > 0)
    Mb->inicipval();
    
  if (Mp->matmodel == nonlocal){ // NUTNO PREDELAT PRO PARALELNI ULOHU
    //  list of adjacent elements
    Gtm->adjacelem (Out);
    //  list of adjacent integration points was not restored 
    //  from file it is necessary to create it from scratch
    if (restore_adjacip() == 0)
      adjacip ();
    save_adjacip();
    ipvolume ();
  }

  //  allocation of eigenstrains
  //  eigenstrains contain strains generated by coupling terms
  if (Mp->eigstrains!=0){
    Mm->alloceigstrains ();
    Mm->alloceigstresses ();
  }

  if (Mp->pore_press!=0){
    Mm->alloceigstresses ();
  }

  //  initiation of temprstrains
  if (Mp->temperature!=0){
    Mm->alloctempstrains ();//tohle jeste rozmyslet??!!
    Mm->alloceigstresses ();
  }

  TM_ip = NULL;
  MT_ip = NULL;
  TM_nod_map = NULL;
  TMipmap = NULL;
  MTipmap = NULL;

  if (Cp->tprob!=fully_coupled_material){
    switch (Cp->dpt){
      case pass_by_closest_ip:
        //data are copied directly between int. points, the SAME int. point must be used on both meshes for the first Tt->ne elements
        //
        // create mapping between TRFEL and MEFEL nodes
        TM_nod_map = new long[Tt->nn];
        create_trfmef_nod_map(TM_nod_map);
        break;
      case pass_by_nodes_comp:
        // data are calculated at nodes with the help of state varibales from the closest int. point,
        // requires the first n elements in both modules to be the same in shape and same ordering of the first m nodes
        if ((Mp->straincomp == 0) || (Mp->strainpos == 1))
          // if not given in the problem description setup, allocate nodal strain arrays
          Mt->alloc_nodal_strains();
      case pass_by_aux_ip:
        // data are calculated in auxiliary int. points,
        // the number of nodes and elements in meshes of particular modules can be independent
        //
        // create and initialize TMipmap and MTipmap arrays
        metr_ip_mapping(Cp->aip_ni, Cp->aip_err);
        // allocate auxiliary integration points in MEFEL
        Mm->alloc_aux_intp(Tm->tnip, TMipmap);
        // allocate auxiliary integration points in TRFEL
        Tm->alloc_aux_intp(Mm->tnip, MTipmap);
        // compute initial values at auxiliary int. points MEFEL
        if (Mb->nico > 0)
          Mm->inic_aipval(Tm->tnip, TMipmap);
        break;
      case pass_by_copy_ip:
        // data are copied directly between int. points, the SAME int. point must be used on both meshes for the first Tt->ne elements
        //
        // create direct mapping between MEFEL and TRFEL integration points
        TM_ip = new long[Tm->tnip];
        MT_ip = new long[Mm->tnip];
        mefel_trfel_ip_mapping_old(TM_ip, MT_ip);
        break;
      default:
        print_err("unknown type of data passing between TRFEL and MEFEL is being rquired (%d)",
                  __FILE__, __LINE__, __func__, int(Cp->dpt));
        abort();
    }
  }
  //Mt->nodedisplacements (); //tohle jeste rozmyslet??!!
  
  xfclose (inc);
  xfclose (inm);
  xfclose (intr);
  
  Pcp->shift_indices ();
  Pmp->shift_indices ();
  Ptp->shift_indices ();
}
