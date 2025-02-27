#include "metrinit.h"

#include "globalc.h"
#include "outdriverm.h"
#include "globmatc.h"
#include "fcsolver.h"
#include "gtopology.h"

#include "seqfilesm.h"
#include "seqfilest.h"

#include "iotools.h"
#include "xfile.h"
#include "globalg.h"
#include "outdriverc.h"
#include "stacktrace.h"
#include "libtrace.h"



void metr_init (int argc, const char *argv[])
{
  long i,n;
  XFILE *in,*inm = NULL,*intr = NULL;
  pkwd_sw kwdsw;

  // ste name of executable for stacktrace output
  set_prgname(argv[0]);
  
  fprintf (stdout,"\n\n ******************************************************\n");
  fprintf (stdout,"            ____ ___ _____ _____ _\n");
  fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
  fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
  fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
  fprintf (stdout,"           |____/___|_|   |_____|_____| METR");
  fflush(stdout);
  fprintf (stderr,"\n\n ******************************************************\n");
    
  // initialize global METR variables to null values
  initnull_globc();
  // initialize global MEFEL variables to null values
  initnull_glob();
  // initialize global TRFEL variables to null values
  initnull_globt();
  
  if (process_argsc(argc, argv, kwdsw))
    abort();
  
  in = xfopen(argv[1],"r");
  if (in==NULL){
    print_err("input file cannot be opened.", __FILE__, __LINE__, __func__);
    abort ();
  }

  //  problem description of all coupled problem
  Cp = new probdescc;
  Cp->kwdsw = kwdsw;
  
  //  problem description of mechanical part of problem
  Mp = new probdesc;
  //  problem description of thermo-hydro part of problem
  Tp = new probdesct;
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
  //  materials of thermo-hydro part of problem
  Cm = new coupmat;
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


  in->warning    = 1;
  in->ignorecase = 1;
  if ((Cp->kwdsw == nokwd) || (Cp->kwdsw == kwd_outd))
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;

  //  problem description of all coupled problem
  Cp->read (in);
  
  if (Cp->tprob != fully_coupled_material){
    
    //  input file for mechanical part
    inm = xfopen(Cp->minfile,"r");
    if (inm==NULL){
      print_err ("input file for mechanical part cannot be opened.", __FILE__, __LINE__, __func__);
      abort ();
    }
    
    //  input file for thermo-hydro part
    intr = xfopen(Cp->tinfile,"r");
    if (intr==NULL){
      print_err ("input file for thermo-hydro part cannot be opened.",__FILE__, __LINE__, __func__);
      abort ();
    }
    
    inm->warning    = 1;
    inm->ignorecase = 1;
    if ((Mp->kwdsw == nokwd) || (Mp->kwdsw == kwd_outd))
      inm->kwdmode = ignore_kwd;
    else
      inm->kwdmode = sequent_mode;
    
    intr->warning    = 1;
    intr->ignorecase = 1;
    if ((Tp->kwdsw == nokwd) || (Tp->kwdsw == kwd_outd))
      intr->kwdmode = ignore_kwd;
    else
      intr->kwdmode = sequent_mode;
    
    //  basic problem data
    Mp->read (inm);
    Tp->read (intr);
    

    Out  = fopen ("mefel.log","wt");
    Outt = fopen ("trfel.log","wt");
    
    if (Mp->kwdsw == kwd_full)
      inm->kwdmode = sequent_mode;
    else
      inm->kwdmode = ignore_kwd;
    
    if (Tp->kwdsw == kwd_full)
      intr->kwdmode = sequent_mode;
    else
      intr->kwdmode = ignore_kwd;
  }
  Outc = fopen ("metr.log","wt");
  
  if (Cp->kwdsw == kwd_full)
    in->kwdmode = sequent_mode;
  else
    in->kwdmode = ignore_kwd;
  
  //  input of node and element data
  if (Cp->tprob != fully_coupled_material){
    Mt->read (inm);
    Tt->read (intr);
  }
  Ct->read (in);
  
  
  // mefel mesh control TKr
  /* if (Mp->tprob == growing_mech_structure){
     if (Mespr==1){
     merr=Mt->mesh_check();
     if (merr > 0){
     fprintf (stdout,"MEFEL Mesh error = %ld\n\n\n",merr);
     exit(0);
     }
     }
     }
  */
    
  // trfel mesh control
  /*if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    if (Mesprt==1){
    merr=Tt->mesh_check();
    if (merr > 0){
    fprintf (stdout,"TRFEL Mesh error = %ld\n\n\n",merr);
    exit(0);
    }
    }
    }
  */

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
    
    //  allocation and initiation of arrays lnso and leso
    //  lnso - list of nodes switched on
    //  leso - list of elements switched on
    //  11/05/2010 TKr
    Gtu->lneso_init ();
  }
  
  if (Cp->tprob != fully_coupled_material){
    
    //  generation of code numbers in mechanical part
    Ndofm = Gtm->codenum_generation (Out);//Ndofm=Gtm->gencodnum ();
    if (Mesprc==1)
      fprintf (stdout,"\n number of degrees of freedom in mechanics    %ld",Ndofm);
    //Mt->store_code_num_elem(); // enable combination of plane elements, springs and bar elems with beams//tady????!!
    
    //  generation of code numbers in transport part
    Ndoft = Gtt->codenum_generation (Outt);//Ndoft=Gtt->gencodnum ();
    if (Mesprc==1)
      fprintf (stdout,"\n number of degrees of freedom in transports    %ld",Ndoft);
    
    if (Cp->tprob==fully_coupled_mech_trans){
      //  generation of code numbers in both parts
      Ndofc = Gtu->codenum_generation (Outc);//Ndofc=Gtu->gencodnum ();
      //  distribution of code numbers to local topologies
      Gtu->cndistr (Gtm,Gtt);
      
      Ndofm=Ndofc;
      Ndoft=Ndofc;
      
      if (Mesprc==1)
	fprintf (stdout,"\n total number of degrees of freedom    %ld",Ndofc);
    }
  }
  
  if (Cp->tprob==fully_coupled_material){
    Ndofc = Gtu->codenum_generation (Outc);
    if (Mesprc==1)
      fprintf (stdout,"\n total number of degrees of freedom    %ld",Ndofc);
  }
  

  /*
  long j, ndofn;
  ivector cn;
  fprintf(Outt, "\n\nTRFEL code numbers:\n");
  for (i=0; i<Gtt->nn; i++)
  {
    ndofn = Gtt->give_ndofn(i);
    reallocv(RSTCKIVEC(ndofn, cn));
    Gtt->give_node_code_numbers(i, cn.a);
    fprintf(Outt, "Node %ld, ndofn=%ld, cn:", i+1, ndofn);
    for (j=0; j<ndofn; j++)
      fprintf(Outt, " %ld", cn[j]);
    fprintf(Outt, "\n");
  }
  fprintf(Out, "\n\nMEFEL code numbers:\n");
  for (i=0; i<Gtm->nn; i++)
  {
    ndofn = Gtm->give_ndofn(i);
    reallocv(RSTCKIVEC(ndofn, cn));
    Gtm->give_node_code_numbers(i, cn.a);
    fprintf(Out, "Node %ld, ndofn=%ld, cn:", i+1, ndofn);
    for (j=0; j<ndofn; j++)
      fprintf(Out, " %ld", cn[j]);
    fprintf(Out, "\n");
  }
  fprintf(Outc, "\n\nMETR code numbers:\n");
  for (i=0; i<Gtu->nn; i++)
  {
    ndofn = Gtu->give_ndofn(i);
    reallocv(RSTCKIVEC(ndofn, cn));
    Gtu->give_node_code_numbers(i, cn.a);
    fprintf(Outc, "Node %ld, ndofn=%ld, cn:", i+1, ndofn);
    for (j=0; j<ndofn; j++)
      fprintf(Outc, " %ld", cn[j]);
    fprintf(Outc, "\n");
  }
  */

  //read trfel function
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->read_gf (intr);
    Gtt->alloc_growstr();
  }

  //  input of material characteristics
  if (Cp->tprob!=fully_coupled_material){
    Mm->read (inm);
    Tm->read (intr);
  }
  if (Cp->tprob==fully_coupled_material){
    Cm->read (in);
  }
  
  //int. points allocation and initiation
  //Tm->intpointalloc ();
  //Tm->intpointinit ();
  Tm->assemble_dof_nameord(Tp->dofname, Tp->ntm, Tp->var_dofid, Tp->tnkv);

  //  allocation of auxiliary arrays on integration points
  //Mm->init_ip_1 ();

  if (Cp->tprob==fully_coupled_mech_trans){
    //Cmu->read (in->file);
    //Cml->read (in->file);
    Cmu->read (in);
    Cml->read (in);
  }

  if (Cp->tprob!=fully_coupled_material){
    //  input of cross section characteristics
    Mc->read (inm);
    Tc->read (intr);
    
    //  input of loading
    Mb->read (inm);
    Tb->read (intr);

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
 
  }else{

  }
  
  //read mefel function
  if (Mp->tprob == growing_mech_structure){
    Gtm->read_gf (inm);
    Gtm->alloc_growstr();
    Mt->alloc_growstr ();
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
    in->kwdmode = ignore_kwd;
  else
    in->kwdmode = sequent_mode;
    
  Outdc->read(in);
   
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

  fclose (Outc);
  fclose (Outt);
  fclose (Out);
  
  char fname[1020];
  if (*Outdc->outfn)
    filename_decomposition(Outdc->outfn,Cp->path,Cp->filename,Cp->suffix);
  else
    filename_decomposition(argv[1],Cp->path,Cp->filename,Cp->suffix);
  sprintf(fname, "%s%s.log", Cp->path, Cp->filename);
  Outc = fopen(fname,"wt");

  //  names of output TRFEL files
  if (*Outdt->outfn)
    filename_decomposition(Outdt->outfn,Tp->path,Tp->filename,Tp->suffix);
  else
    filename_decomposition(Cp->tinfile,Tp->path,Tp->filename,Tp->suffix);
  if ((strlen(Tp->path) == 0) && (strlen(Tp->filename) == 0))
  {
    sprintf(fname, "trfel.log");
    Outt = fopen(fname,"at");
  }
  else
  {
    sprintf(fname, "%s%s.log", Tp->path, Tp->filename);
    Outt = fopen(fname,"wt");
  }
  
  //  names of output MEFEL files
  if (*Outdm->outfn)
    filename_decomposition(Outdm->outfn,Mp->path,Mp->filename,Mp->suffix);
  else
    filename_decomposition(Cp->minfile,Mp->path,Mp->filename,Mp->suffix);
  if ((strlen(Mp->path) == 0) && (strlen(Mp->filename) == 0))
  {
    sprintf(fname, "mefel.log");
    Out = fopen(fname,"at");
  }
  else
  {
    sprintf(fname, "%s%s.log", Mp->path, Mp->filename);
    Out = fopen(fname,"wt");
  }
  fprintf(stdout, "\n Log file name : '%s'", fname);
  
  //Mm->stra.init(strain);
  //Mm->stre.init(stress);

  //  allocation of auxiliary arrays on integration points
  Mm->init_ip_2 ();

  //  nodes allocation
  Tt->alloc_nodes ();
  //  allocation of arrays defined on nodes
  Mt->alloc_nodes_arrays ();
  //  number of components in arrays used in outdriver
  //  maximum number of components in array strain/stress on nodes and elements
  Mt->give_maxncompstr(Mm->max_ncompstrn, Mm->max_ncompstre);
  //  maximum number of components in array other on nodes and elements
  Mt->give_maxncompo(Mm->max_ncompothern, Mm->max_ncompothere);

  //  prescribed displacements on elements investigation
  Mt->elemprescdisp ();
  //  reactions on elements investigation
  if (Mp->reactcomp==1){
    Mt->comreac ();
  }

  //fprintf(Out, "\n\nIntegration point global coordinates:\n\n");
  //matrix ipcoordmat;
  //for(i=0; i<Mt->ne; i++)
  //{    
  //  long tnip = Mt->give_tnip(i);
  //  reallocm(RSTCKMAT(tnip, 3, ipcoordmat));
  //  Mt->give_ipcoord_elem(i, ipcoordmat);
  //  fprintf(Out, "Elem %4ld: %le %le %le\n", i+1, ipcoordmat[0][0], ipcoordmat[0][1], ipcoordmat[0][2]);
  //  for(long j=1; j<tnip; j++)
  //  {      
  //    fprintf(Out, "%11s", " ");
  //    fprintf(Out, "%le %le %le\n", ipcoordmat[j][0], ipcoordmat[j][1], ipcoordmat[j][2]);
  //  }
 // }
 // fflush(Out);

  //  computation of initial values at integration points from initial nodal values
  if (Mb->nico > 0)
    Mb->inicipval();
    
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

  if (Mb->give_num_mstrain_comp(0))
  {
    if (Mp->matmodel != nonlocal) // int. point volumes have been already computed for non-local models
      ipvolume ();
    Mm->allocmacrostresses();
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
        // data are copied directly between nodes, 
        // the first n-th elements must be defined on both meshes with the SAME node numbers, MEFEL elements may use higher number of nodes,
        // but first m nodes on elements must be the SAME
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
        break;
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
        // initialize eigstrid array at auxiliary int. points and compute initial values of eigenstrains or eigenstresses
        if ((Mp->eigstrains > 0) && (Mp->eigstrains != 3)){
          Mm->init_aipeigstr(Tm->tnip, TMipmap, Mp->time);
        }
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

  xfclose (in);
  xfclose (inm);
  xfclose (intr);
}



/**
  The function prints help for command line arguments of metr execuatble to the given output text file.
  @param out[in] - pointer to the opened text file
  @param prgname[in] - string with executable name

  @return The function does not return anything. 

  Created by Tomas Koudelka
*/
void print_helpc(FILE *out, const char *prgname)
{
  fflush(out);
  fprintf(out, "\n Use : %s infile [-kwd=i] [-numth=%%d]\n", prgname);
  fprintf(out, " where :\n *   infile is input file name\n");
  fprintf(out, " *   -kwd={0,1,2,3,4} is optional switch which controls keyword usage in the input file\n");
  fprintf(out, "           0 means no keywords (default value)\n");
  fprintf(out, "           1 means keywords only in the problem description section\n");
  fprintf(out, "           2 means keywords only in the outdriver section\n");
  fprintf(out, "           3 means keywords in the problem description and outdriver sections\n");
  fprintf(out, "           4 means full usage of keywords\n\n");
  fprintf(out, " *   -numth={%%d} is optional switch which specifies number of threads used in OpenMP\n");
}



/**
  The function processes arguments passed by the command line.

  @param argc[in] - the total number of commandline argumets, i.e. dimension of array argv
  @param argv[in] - pointer to the string array with particular command line arguments

  @ratval 0 - on success
  @retval 1 - error in processing of arguments

  Created by Tomas Koudelka
*/
long process_argsc(int argc, const char *argv[], pkwd_sw &kwdsw)
{
 #if defined(INC_OPENMP) || defined(INC_CHMD)
  // set default number of threads
  Numth = 8;
 #endif
  switch (argc)
  {
    case 0:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too few arguments are used\n\n");
      print_helpc(stderr, "pmetr|metr");
      return 1;
    case 1:
      fflush(stdout);
      fprintf(stdout, "\nInput file has not been specified\n");
      print_helpc(stderr,argv[0]);
      return 1;
    case 2:
      kwdsw = nokwd;
      break;
    default:
      if (process_optargsc(argc, argv, kwdsw))
        return 1;
      break;
  }
 #if defined(INC_OPENMP) || defined(INC_CHMD)
  fprintf(stdout, "\n Number of threads used in  METR: %d\n", Numth);
 #endif
  return 0;
}



long process_optargsc(int argc, const char *argv[], pkwd_sw &kwdsw)
{
  const char *ptrkwd;
  if (argc > 2){
    for (long i=2; i<argc; i++){
      ptrkwd = strstrcis(argv[i], "-kwd=");
      if (ptrkwd)
      {
        kwdsw = pkwd_sw(-1);
        sscanf(ptrkwd+5, "%d", (int *)&kwdsw);
        if ((kwdsw < 0) || (kwdsw > 4))
        {
          fflush(stdout);
          fprintf(stderr, "\n\nError - unknown -kwd value %d is obtained\n\n", kwdsw);
          print_helpc(stderr,argv[0]);
          return 1;
        }
        continue;
      }
      else
        kwdsw = nokwd;
    #if defined(INC_OPENMP) || defined(INC_CHMD)
      ptrkwd = strstrcis(argv[i], "-numth=");
      if (ptrkwd)
      {
        Numth = -1;
        sscanf(ptrkwd+7, "%d", &Numth);
        if (Numth < 1)
        {
          fflush(stdout);
          fprintf(stderr, "\n\nError - wrong number of threads (%d), it must be > 0\n\n", Numth);
          print_helpc(stderr,argv[0]);
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
      fflush(stdout);
      fprintf(stderr, "\n\nError - argument %s was not recognized\n", argv[i]);
      print_helpc(stderr,argv[0]);
      return 1;
    }
  }
  else{
    fprintf(stderr, "\n\nError - number of arguments should be more than 2 in this function (argc=%d)\n\n", argc);
    print_helpc(stderr,argv[0]);
    return 1;
  }
  return 0;
}
