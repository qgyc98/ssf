#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#include "globalt.h"

#include "stochdrivert.h"
#include "trfelinit.h"
#include "partrfinputfiles.h"



/**
  Function reads domain topology file in, format of the file is sifel format.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description


  Function returns:
  @retval 0 : on success
*/
long input_siftop(XFILE *in, descript *d)
{
  long ret;

  fprintf(stdout, "\n\n\n Reading of mesh domain topology . . .");

  ret = Top->read(in, d->paral, d->redgn);

  return(ret);
}


partrfinputfile::partrfinputfile(int argc, const char *argv[])
{
  Argc = argc;
  Argv = argv;  
  topf = NULL;
  err = 0;
  logname=NULL;
  tmp=NULL;
}


// descructor
partrfinputfile::~partrfinputfile()
{
  delete Top;
  
}




void partrfinputfile::print_in(FILE *out)
{
  //  basic problem data
  Tp->print (out);

  //  node and element data
  Tt->print (out);

  //  material characteristics
  Tm->print (out);
  
  
  //  cross section characteristics
  Tc->print (out);
  
  //  boundary conditions
  Tb->print (out);

  //  prescribed values on elements investigation
  //Tt->elemprescvalue();
  //  searching for elements with nodes where source is defined
  Tb->elemsource ();

  // initial conditions
  Lsrst->initcondprint (out);

  // output description
  Outdt->print(out);

}


void partrfinputfile::print_partrf_in(FILE *out)
{
  //  domain nuber
  fprintf(out,"%s\n",Argv[4]);

  //  basic problem data
  Tp->print (out);
  
  fprintf(out,"\n\n\n# part describing parallel solver\n");
  fprintf(out,"2 # global numbering of bounding nodes\n");
  fprintf(out,"1 # type of domain decomposition method - Schur complement method\n");
  fprintf(out,"1  3  3# type of solver for reduced system\n");
  //fprintf(out,"10     # type of storage of reduced system matrix = compressed rows\n");
  //fprintf(out,"20 # type of solver for reduced system - conjugate gradient method\n");
  //fprintf(out,"2000 # maximum number of iterations\n");
  //fprintf(out,"1.0e-8 # required error\n");
  //fprintf(out,"0      # no preconditioner will be used \n");
  fprintf(out,"\n\n");

  //  node and element data
  partrftop_print (out);

  //  material characteristics
  Tm->print (out);
  
  
  //  cross section characteristics
  Tc->print (out);
  
  //  boundary conditions
  partrfload_print (out);

  //  prescribed values on elements investigation
  //Tt->elemprescvalue();
  //  searching for elements with nodes where source is defined
  //Tb->elemsource ();

  //  initial conditions
  partrf_initcondprint (out);

  // output description
  Outdt->print(out);

}


/**
   function prints basic data about topology for one domain
   
   @param out - output stream

   19/12/2012, TKr
*/
void partrfinputfile::partrftop_print (FILE *out)
{
  long i,j,el_type;
  
  // ************
  //  node data
  // ************
  fprintf (out,"\n\n%ld #number of nodes",Top->nn);
  
  // the first property of nodes must contain global node number
  for (i=0;i<Top->nn;i++){
    fprintf (out,"\n %ld",i+1);
    Gtt->gnodes[Top->nodes[i].prop[0]-1].print (out);
    Tt->nodes[Top->nodes[i].prop[0]-1].print (out);
  }
  
  // ********************
  //  constrained nodes
  // ********************
  fprintf (out,"\n\n%ld #number of constrained nodes",Top->nn);
  for (i=0;i<Top->nn;i++){
    fprintf (out,"\n %ld",i+1);
    Gtt->gnodes[Top->nodes[i].prop[0]-1].print_constr (Gtt->dofcontr,out);
  }
  
  // ************************************
  //  master node only for homogenization
  // ************************************
  if(Tp->homogt == 1){
    fprintf (out,"\n %ld #master node number\n\n",Tt->mnn);
  }

  // ***************
  //  element data
  // ***************
  fprintf (out,"\n\n %ld #number of elements\n",Top->ne);

  // the first surface property of elements must contain global element number
  for (i=0;i<Top->ne;i++){
    
    el_type = give_elemtype(Top->elements[i].type);
    
    fprintf (out,"%ld %d",i+1, (int)el_type);
    for (j = 0; j < Top->elements[i].nne; j++)
      fprintf (out," %ld",Top->elements[i].nodes[j]+1);
    
    // code numbers on one element
    fprintf (out," %ld ",Gtt->gelements[Top->elements[i].propsurf[0]-1].cne);
    
    if (Gtt->gelements[Top->elements[i].propsurf[0]-1].cne==1 || Gtt->gelements[Top->elements[i].propsurf[0]-1].cne==2){
      for (i=0;i<Gtt->gelements[Top->elements[i].propsurf[0]-1].ndofe;i++){
	fprintf (out," %ld ",Gtt->gelements[Top->elements[i].propsurf[0]-1].cn[i]+1);
      }
    }

    // cross-section priting
    fprintf (out," %d",Tt->elements[Top->elements[i].propsurf[0]-1].crst);
    if (Tt->elements[Top->elements[i].propsurf[0]-1].crst!=0){
      fprintf (out," %ld",Tt->elements[Top->elements[i].propsurf[0]-1].idcs+1);
    }
    
    //  printing of material tpyes and ids
    Tt->elements[Top->elements[i].propsurf[0]-1].printmat (Tp->ntm*Tp->ntm,out);  
  }


  // ******************
  //  paralel code data
  // ******************
  fprintf (out,"\n\n # local to global map\n");
  for (i=0;i<Top->nn;i++){
    fprintf(out,"\n %ld",Top->gnn[i]);
  }
}


/**
   function prints loading one domain
   
   @param out - output stream

   19/12/2012, TKr
*/
void partrfinputfile::partrfload_print (FILE *out)
{
  long i;
  
  fprintf(out,"\n\n## boundary conditions:\n");
  fprintf (out,"\n");
  fprintf (out,"\n %ld #number of loadcases",Tb->nlc);
  for (i=0;i<Tb->nlc;i++){
    partrfloadcase_print (out,i);
  }
}


/**
   function prints load case characteristics for one domain
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   TKr, 3.1.2006
*/
void partrfinputfile::partrfloadcase_print (FILE *out,long lcid)
{
  long i,j,neb;
  
  // ********************
  //  prescribed values
  // ********************
  fprintf(out,"\n## dirichlet b.c.:");
  fprintf (out,"\n%ld\n\n",Tb->lc[lcid].npv);
  for (i=0;i<Tb->lc[lcid].npv;i++){
    Tb->lc[lcid].pv[i].print(out);
  }
  
  // *******************
  //  quantity sources
  // *******************
  //  number of quantity sources //toto nemuze fungovat -> opravit??!!
  fprintf(out,"\n## sources:");
  fprintf (out,"\n%ld\n\n",Tb->lc[lcid].nqs);
  if (Tb->lc[lcid].nqs>0){
    //  models of source
    for (i=0;i<Tb->lc[lcid].nqs;i++){
      fprintf (out,"\n %ld",i+1);
      Tb->lc[lcid].sour[i].print (out);
    }
    
    //  list of nodes with defined source
    fprintf (out,"\n\n %ld",Tb->lc[lcid].nnqs);
    for (i=0;i<Tb->lc[lcid].nnqs;i++){
      fprintf (out,"\n %ld %ld",Tb->lc[lcid].lnqs[i][0]+1,Tb->lc[lcid].lnqs[i][1]+1);
    }
    
    //  list of elements with defined source
    fprintf (out,"\n\n %ld",Tb->lc[lcid].neqs);
    for (i=0;i<Tb->lc[lcid].neqs;i++){
      fprintf (out,"\n %ld %ld",Tb->lc[lcid].leqs[i][0]+1,Tb->lc[lcid].leqs[i][1]+1);
    }
  }
  
  // **********************************************
  //  number of elements with boundary conditions
  // **********************************************
  //computing:
  neb = 0;
  for (i=0;i<Tb->lc[lcid].neb;i++){
    for (j=0;j<Top->ne;j++){
      if(Tb->lc[lcid].elemload[i].eid == (Top->elements[j].propsurf[0]-1))
      neb++;
    }
  }

  // printing:
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",neb);
  for (i=0;i<Tb->lc[lcid].neb;i++){
    for (j=0;j<Top->ne;j++){
      if(Tb->lc[lcid].elemload[i].eid == (Top->elements[j].propsurf[0]-1)){
	//Tb->lc[lcid].elemload[i].print (out,lcid);//tady
	partrfelemload_print(lcid,i,j,out);
      }
    }
  }
 
  // ***********************************
  //  nodal values defined on boundary
  // ***********************************
  fprintf(out,"\n## nodal values:");
  fprintf (out,"\n%ld\n\n",Tb->lc[lcid].nnv);
  for (i=0;i<Tb->lc[lcid].nnv;i++){
    Tb->lc[lcid].nodval[i].print (out);
  }

  // **********************
  //  climatic conditions
  // **********************
  fprintf(out,"\n## climatic conditions:");
  fprintf (out,"\n%ld\n",Tb->lc[lcid].ncc);
  if (Tb->lc[lcid].ncc>0){
    fprintf (out,"%ld\n",Tb->lc[lcid].trcc);
    for (i=0;i<Tb->lc[lcid].ncc;i++){
      fprintf (out,"\n\n %ld",i+1);
      Tb->lc[lcid].climcond[i].print (out);
    }
  }
  
  
  
  //part for multpvalt:
  //  number of time functions
  fprintf(out,"\n## time functions:");
  fprintf (out,"\n%ld\n\n",Tb->lc[lcid].nmultfunc);
  
  for(i=0;i<Tb->lc[lcid].nmultfunc;i++){
    
    fprintf (out,"\n %ld",i+1);
    
    fprintf(out,"\n %d ",(int)Tb->lc[lcid].tfunc);
    switch (Tb->lc[lcid].tfunc){
      
    case tab:{
      Tb->lc[lcid].tabf[i].print(out);
      break;
    }
    default:{
    }
    }
  }

  // ************************
  //  data for homogenization
  // ************************
  fprintf(out,"\n## homogenization:\n");
  if(Tp->homogt == 1){
    fprintf (out,"\n  %lf",Tb->lc[lcid].masterval);
    for (i=0;i<Tb->lc[lcid].ndim;i++)
      fprintf (out,"  %lf",Tb->lc[lcid].mastergrad[i]);
    fprintf (out,"\n");
  }

}


/**
   function prints loading of elements on one domain
   
   @param lcid - loadcase id
   @param leid - load global element eid
   @param deid - domain element eid
   @param out  - output stream

   19/12/2012, TKr
*/
void partrfinputfile::partrfelemload_print(long lcid,long leid,long deid,FILE *out)
{
  long i;
  
  //  domain element id
  fprintf (out,"\n\n %ld ",deid+1);
  
  //  loop over boundary objects (end nodes, edges or surfaces)
  for (i=0;i<Tb->lc[lcid].elemload[leid].nbo;i++){
    if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
      fprintf(out, "%ld %ld %ld %ld\n", 
              //  type of boundary condition
              Tb->lc[lcid].elemload[leid].bcf[i],
              //  id of nodal values
              Tb->lc[lcid].elemload[leid].nvidf[i],
              //  id of transmission coefficients
              Tb->lc[lcid].elemload[leid].trcidf[i],
              //  id of transmission/radiation coefficients
              Tb->lc[lcid].elemload[leid].trridf[i]);
    }
    else{
      fprintf (out,"\n %d", (int)Tb->lc[lcid].elemload[leid].bc[i]);

      if (Tb->lc[lcid].elemload[leid].bc[i]==2){
	//  prescribed flux (Neumann boundary conditon)
	//  id of exterior prescribed values
	fprintf (out," %ld",Tb->lc[lcid].elemload[leid].nvid[i]+1);
      }

      if (Tb->lc[lcid].elemload[leid].bc[i]>10){
	//  prescribed transmission (Newton boundary conditon)
	
	//  id of exterior prescribed values
	fprintf (out," %ld",Tb->lc[lcid].elemload[leid].nvid[i]+1);
	//  id of transmission coefficients
	fprintf (out," %ld",Tb->lc[lcid].elemload[leid].trcid[i]+1);
	//  id of transmission/radiation coefficients
	fprintf (out," %ld",Tb->lc[lcid].elemload[leid].trrid[i]+1);
      }
      if (Tb->lc[lcid].elemload[leid].bc[i]==3 || Tb->lc[lcid].elemload[leid].bc[i]==4){
	// prescribed detail climatic conditions
	
	//  id of detail climatic conditions
	fprintf (out," %ld",Tb->lc[lcid].elemload[leid].nvid[i]+1);
      }
      
    }
  }
}



/**
   function prints initial conditios for one domain
   
   @param out - output stream

   19/12/2012, TKr
*/
void partrfinputfile::partrf_initcondprint (FILE *out)
{
  long i,j,k,ndofn;
  
  // the first property of nodes must contain global node number
  //fprintf (out,"\n\n%ld #initial conditions",Top->nn);
  fprintf (out,"\n\n#initial conditions");
  
  switch (Tp->tprob){
  case nonlinear_stationary_problem:
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
  case discont_nonlin_nonstat_problem:{
    for (i=0;i<Top->nn;i++){
      ndofn=Gtt->give_ndofn (Top->nodes[i].prop[0]-1);
      fprintf (out,"\n %ld  ",i+1);
      for (j=0;j<ndofn;j++){
	k=Gtt->give_dof ((Top->nodes[i].prop[0]-1),j)-1;
	if (k>-1){ fprintf (out," %lf",Lsrst->lhsi[k]);}
	else {fprintf (out," 0.0");}
      }
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function returns element type for TRFEL input file
   
   @param etype - element type in jktk mesh format

   19/12/2012, TKr
*/
long partrfinputfile::give_elemtype(long etype)
{
  long type;
  
  switch(etype){
  case 7:{
    type = 220;
    break;
  }
  default:{
    print_err("unknown type of element is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return(type);
}


int partrfinputfile::print_partrfinputfile()
{
  double ts,te;
  FILE *out;
  char *path,*filename,*suffix;
  const char *dn;

  ts = clock();

  dn = Argv[4];
  
  char fname[1020];
  filename_decomposition(Argv[3],path,filename,suffix);

  sprintf(fname, "%strfel_%s.in", path, filename);
  out = fopen(fname,"wt");

  //  input file name
  fprintf(stdout,"\n Printing of PARTRF domain input file no. %s\n",dn);

  //  printing for sequential input file
  //print_in(out); //commented

  //  printing for parallel input file
  print_partrf_in(out);

  fclose(out);
  te = clock();
  printf("\n Time of Printing of PARTRF domain input file %le s\n\n\n",(te-ts)/(double)CLOCKS_PER_SEC);

  delete [] path; delete [] filename; delete [] suffix;
  return(0);  
}



int partrfinputfile::read_topology()
{
  double ts,te;
  
  // set tprob for reading of mesh
  ts = clock();
  //  topology file name
  topf = xfopen(Argv[3],"r");
  if (d.topf==NULL){
    fprintf (stderr,"\n Input file for domain topology has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  topf->warning = 1;
  topf->kwdmode = ignore_kwd;
  topf->ignorecase = 1;
  err = input_siftop(topf,&d);
  xfclose(topf);
  if (err){
    print_err("\n Reading of mesh topology failed\n", __FILE__, __LINE__, __func__);
    return(3);
  }
  else{
    printf(" OK");
  }
  te = clock();
  printf("\n Time of reading topology %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  return(0);  
}


int partrfinputfile::run()
{
  double ts,te;
  Top = new siftop;
 
  if (Argc < 5){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name [-kwd=0] domain_input_topology_file_name domain_number\n\n",Argv[0]);
    delete Gtt; delete Top;
    return(1);
  }
  
  fprintf(stdout,"\n\n\n Reading of TRFEL main (sequential) input file . . .\n");
  ts = clock();
  //  reading of trfel input file = program initiation and data reading
  trfel_init (Argc-2,Argv);
  te = clock();
  fprintf(stdout,"\n Time of reading TRFEL main (sequential) input file %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);

  //  reading domain topology
  //  this is for paralel version
  d.paral=1;
  d.redgn = 1;
  err=read_topology();

  //  printing of domain PARTRF input file
  err=print_partrfinputfile();

  return(0);
}


int main (int argc,const char *argv[])
{
  if (argc < 5){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name [-kwd=0] domain_input_topology_file_name domain_number\n\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  partrfinputfile *partrffile;
  double tws,twe;
  struct tm *local;
  time_t t;
  
  set_prgname(argv[0]);
  printf("---------------------------------------------------------------\n");
  printf("*** PARTRF INPUT FILE GENERATION FROM SEQUENTIAL INPUT FILE ***\n");
  printf("---------------------------------------------------------------\n");
  t = time(NULL);
  local = localtime(&t);
  printf("Program started  %s\n", asctime(local));
  
  tws = clock();
  partrffile = new partrfinputfile(argc,argv);
  partrffile->run();
  delete partrffile;
  twe = clock();
  printf(" -----------------------------------------------------\n");
  printf ("Whole elapsed time of input file generation %le s\n\n\n",(twe-tws)/(double)CLOCKS_PER_SEC);
}
 
