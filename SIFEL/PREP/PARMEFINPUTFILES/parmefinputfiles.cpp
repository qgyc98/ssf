#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#include "global.h"
#include "stochdriver.h"
#include "mefelinit.h"
#include "parmefinputfiles.h"
#include "node.h"
#include "loadn.h"
#include "element.h"
#include "loadcase.h"
#include "dloadcase.h"

/**
  Function reads domain topology file in, format of the file is sifel format.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description


  Function returns:
  @retval 0 : on success
*/
long input_siftop(XFILE *in, descrip *d)
{
  long ret;

  fprintf(stdout, "\n\n\n Reading of mesh domain topology . . .");

  ret = Top->read(in, d->paral, d->redgn);

  return(ret);
}


parmefinputfile::parmefinputfile(int argc, const char *argv[])
{
  Argc = argc;
  Argv = argv;  
  topf = NULL;
  err = 0;
  logname=NULL;
  tmp=NULL;
}


// descructor
parmefinputfile::~parmefinputfile()
{
  //delete Gtm;
  delete Top;
  
}




void parmefinputfile::print_in(FILE *out)
{
  //  basic problem data
  Mp->print (out);

  //  node and element data
  //Mt->print (out); //dodelat??!!

  //  material characteristics
  Mm->print (out);
  
  //  cross section characteristics
  Mc->print (out);
  
  //  loading
  //Mb->print (out); //dodelat??!!

  // output description
  Outdm->print(out);

}


void parmefinputfile::print_parmef_in(FILE *out)
{
  //  domain nuber
  fprintf(out,"%s\n",Argv[4]);

  //  basic problem data
  Mp->print (out);
  
  fprintf(out,"\n\n\n# part describing parallel solver\n");
  fprintf(out,"2 # global numbering of bounding nodes\n");
  fprintf(out,"1 # type of domain decomposition method - Schur complement method\n");
  fprintf(out,"1  # reduced system of equations is solved sequentially on master processor\n");
  //toto:
  fprintf(out,"2  # type of storage of reduced system matrix = skyline\n");
  fprintf(out,"2  # solver of system of linear equations = LDL method\n");
  //nebo toto:
  //fprintf(out,"10     # type of storage of reduced system matrix = compressed rows\n");
  //fprintf(out,"20     # type of solver for reduced system - conjugate gradient method\n");
  //fprintf(out,"2000   # maximum number of iterations\n");
  //fprintf(out,"1.0e-8 # required error\n");
  //fprintf(out,"0      # no preconditioner will be used \n");
  fprintf(out,"\n\n");

  //  node and element data
  parmeftop_print (out);

  //  material characteristics
  Mm->print (out);
  
  //  cross section characteristics
  Mc->print (out);
  
  //  loading
  parmefload_print (out);

  // output description
  Outdm->print(out);

}


/**
   function prints basic data about topology for one domain
   
   @param out - output stream

   19/12/2012, TKr
*/
void parmefinputfile::parmeftop_print (FILE *out)
{
  long i,j,el_type;
  
  // ************
  //  node data
  // ************
  fprintf (out,"\n\n%ld #number of nodes",Top->nn);
  
  // the first property of nodes must contain global node number
  for (i=0;i<Top->nn;i++){
    fprintf (out,"\n %ld",i+1);
    Gtm->gnodes[Top->nodes[i].prop[0]-1].print (out);
    Mt->nodes[Top->nodes[i].prop[0]-1].print(out);
  }
  
  // ********************
  //  constrained nodes
  // ********************
  fprintf (out,"\n\n%ld #number of constrained nodes",Top->nn);
  for (i=0;i<Top->nn;i++){
    fprintf (out,"\n %ld",i+1);
    Gtm->gnodes[Top->nodes[i].prop[0]-1].print_constr (Gtm->dofcontr,out);
  }
  
  // ***************
  //  element data
  // ***************
  fprintf (out,"\n\n %ld #number of elements\n",Top->ne);

  // the first surface property of elements must contain global element number
  for (i=0;i<Top->ne;i++){
    
    el_type = give_elemtype(Top->elements[i].type);
    
    fprintf (out,"%ld %ld",i+1, el_type);
    for (j = 0; j < Top->elements[i].nne; j++)
      fprintf (out," %ld",Top->elements[i].nodes[j]+1);
    
    // code numbers on one element
    fprintf (out," %ld ",Gtm->gelements[Top->elements[i].propsurf[0]-1].cne);
    
    if (Gtm->gelements[Top->elements[i].propsurf[0]-1].cne==1 || Gtm->gelements[Top->elements[i].propsurf[0]-1].cne==2){
      for (i=0;i<Gtm->gelements[Top->elements[i].propsurf[0]-1].ndofe;i++){
	fprintf (out," %ld ",Gtm->gelements[Top->elements[i].propsurf[0]-1].cn[i]+1);
      }
    }

    // cross-section priting
    fprintf (out," %d",Mt->elements[Top->elements[i].propsurf[0]-1].crst);
    if (Mt->elements[Top->elements[i].propsurf[0]-1].crst!=0){
      fprintf (out," %ld",Mt->elements[Top->elements[i].propsurf[0]-1].idcs+1);
    }
    
    //  printing of material tpyes and ids
    Mt->elements[Top->elements[i].propsurf[0]-1].printmat (out);
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
void parmefinputfile::parmefload_print (FILE *out)
{
  long i;
  
  fprintf(out,"\n\n## loads:\n");
  fprintf (out,"\n");
  fprintf (out,"\n # number of load cases:");

  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<2*Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    Mb->printinic(out);
    break;
  }
  case geom_nonlinear_statics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<2*Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    Mb->printinic(out);
    break;
  }
  case eigen_dynamics:{
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
      dparmefloadcase_print (out,i);
    }
    break;
  }
  case mech_timedependent_prob:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      dparmefloadcase_print (out,i);
    }
    Mb->printinic(out);
    break;
  }
  case growing_mech_structure:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      dparmefloadcase_print (out,i);
    }
    Mb->printinic(out);
    break;
  }
  case earth_pressure:{
    fprintf (out,"\n  %ld",Mb->nlc);
    if (Mp->tepsol < epplast_sol)
    {
      for (i=0;i<Mb->nlc;i++){
        parmefloadcase_print (out,i);
        dparmefloadcase_print (out,i);
      }
      Mb->printinic(out);
    }
    else
    {
      for (i=0;i<2*Mb->nlc;i++){
        parmefloadcase_print (out,i);
      }
      Mb->printinic(out);
    }
    break;
  }
  case layered_linear_statics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }
    
  case lin_floating_subdomain:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }

  case nonlin_floating_subdomain:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<2*Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }

  case nonlinear_dynamics:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
      dparmefloadcase_print (out,i);
    }
    break;
  }
    
  case hemivar_inequal:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }
    
  case load_balancing:{
    fprintf (out,"\n  %ld",Mb->nlc);
    for (i=0;i<Mb->nlc;i++){
      parmefloadcase_print (out,i);
    }
    break;
  }

  default:{
    print_err("unknown problem type", __FILE__, __LINE__, __func__);
  }
  }

  fprintf(out,"\n\n## eigenstrains:\n");
  Mb->print_eigenstrains (out);
}

/**
   function prints time dependent load case characteristics for one domain
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   TKr, 3.1.2006
*/
void parmefinputfile::dparmefloadcase_print (FILE *out,long lcid)
{
  long i,j;
  
  //  type of time-dependent load
  fprintf (out,"\n #type of time-dependent load:");
  fprintf (out,"\n  %d\n",(int)Mb->dlc[lcid].tdl);
  
  switch (Mb->dlc[lcid].tdl){
  case timeindload:{
    //  time-independent load case multiplied by a time function
    //  number of subload cases
    fprintf (out,"\n #time-independent load case:");
    fprintf (out,"\n #number of subloadcases:");
    fprintf (out,"\n\n %ld",Mb->dlc[lcid].nslc);
    
    for (i=0;i<Mb->dlc[lcid].nslc;i++){
      //  subload cases
      parmefsubloadcase_print (out,lcid,i);
      //  time function
      Mb->dlc[lcid].gf[i].print (out);
    }
    break;
  }
  case seismicload:{
    //  seismic load
    fprintf (out,"\n seismic load:");
    //Mb->dlc[lcid].stool.print (out);//doplnit??!!
    break;
  }
  case responsespectrum:{
    //  
    fprintf (out,"\n response spectrum:");
    //Mb->dlc[lcid].stool.print (out);//doplnit??!!
    break;
  }
  case timedepload:{
    //  time-dependent load case
    //  all nodes, elements, prescribed displacements and other
    //  quantities are equiped with their own time functions
    
  
    //  loaded nodes
    fprintf(out,"\n## loaded nodes:");
    fprintf (out,"\n%ld\n\n",Mb->lc[lcid].nln);
    for (i=0;i<Mb->lc[lcid].nln;i++){
      for (j=0;j<Top->nn;j++){
	if(Mb->lc[lcid].lon[i].nid == (Top->nodes[j].prop[0]-1)){
	  Mb->lc[lcid].lon[i].print(out);
	}
      }
    }
    
    /*  
    //  loaded elements
    fscanf (in,"%ld",&nle);
    loe = new dloadel [nle];
    for (i=0;i<nle;i++){
    loe[i].read (in);
    }*/
    

    //  prescribed displacements
    fprintf(out,"\n##   prescribed displacements:");
    fprintf (out,"\n%ld\n\n",Mb->lc[lcid].npd);
    for (i=0;i<Mb->lc[lcid].npd;i++){
      fprintf (out,"\n %ld",i+1);
      fprintf (out,"  %lf",Mb->lc[lcid].pd[i]);
    }
    break;
  }
    
  default:{
    print_err("unknown type of load case is required", __FILE__, __LINE__, __func__);
  }
  }
  
}


/**
   function prints load case characteristics for one domain
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   TKr, 3.1.2006
*/
void parmefinputfile::parmefloadcase_print (FILE *out,long lcid)
{
  long i,j,neb;
  
  //  loaded nodes
  fprintf(out,"\n## loaded nodes:");
  fprintf (out,"\n%ld\n\n",Mb->lc[lcid].nln);
  for (i=0;i<Mb->lc[lcid].nln;i++){
    for (j=0;j<Top->nn;j++){
      if(Mb->lc[lcid].lon[i].nid == (Top->nodes[j].prop[0]-1)){
	Mb->lc[lcid].lon[i].print(out);
      }
    }
  }

  //  loaded elements
  //  number of loaded elements
  //  computing:
  neb = 0;
  for (i=0;i<Mb->lc[lcid].nle;i++){
    for (j=0;j<Top->ne;j++){
      if(Mb->lc[lcid].loe[i].eid == (Top->elements[j].propsurf[0]-1))
	neb++;
    }
  }
  
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",neb);
  for (i=0;i<Mb->lc[lcid].nle;i++){
    for (j=0;j<Top->ne;j++){
      if(Mb->lc[lcid].loe[i].eid == (Top->elements[j].propsurf[0]-1)){
	//Mb->lc[lcid].loe[i].print(out);//tady
	parmefelemload_print(lcid,i,j,out);
      }
    }
  }


  //  prescribed displacements
  fprintf(out,"\n##   prescribed displacements:");
  fprintf (out,"\n%ld\n\n",Mb->lc[lcid].npd);
  for (i=0;i<Mb->lc[lcid].npd;i++){
    fprintf (out,"\n %ld",i+1);
    fprintf (out,"  %lf",Mb->lc[lcid].pd[i]);
  }
  
  //  prescribed temperature changes
  fprintf(out,"\n##   temperature changes:");
  fprintf(out,"\n%ld\n\n",Mb->lc[lcid].tempchang);
  
  if (Mb->lc[lcid].tempchang==1 || Mb->lc[lcid].tempchang==2){
    // the first property of nodes must contain global node number
    for (i=0;i<Top->nn;i++){
      fprintf(out,"\n%lf",Mb->lc[lcid].pt[Top->nodes[i].prop[0]-1]);
    } 
  }
}

/**
   function prints sub-load case characteristics for one domain
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   TKr, 3.1.2006
*/
void parmefinputfile::parmefsubloadcase_print (FILE *out,long lcid, long slcid)
{
  long i,j,neb;
  
  //  loaded nodes
  fprintf(out,"\n## loaded nodes:");
  fprintf (out,"\n%ld\n\n",Mb->dlc[lcid].slc[slcid].nln);
  for (i=0;i<Mb->dlc[lcid].slc[slcid].nln;i++){
    for (j=0;j<Top->nn;j++){
      if(Mb->dlc[lcid].slc[slcid].lon[i].nid == (Top->nodes[j].prop[0]-1)){
	Mb->dlc[lcid].slc[slcid].lon[i].print(out);
      }
    }
  }

  //  loaded elements
  //  number of loaded elements
  //  computing:
  neb = 0;
  for (i=0;i<Mb->dlc[lcid].slc[slcid].nle;i++){
    for (j=0;j<Top->ne;j++){
      if(Mb->dlc[lcid].slc[slcid].loe[i].eid == (Top->elements[j].propsurf[0]-1))
	neb++;
    }
  }
  
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",neb);
  for (i=0;i<Mb->dlc[lcid].slc[slcid].nle;i++){
    for (j=0;j<Top->ne;j++){
      if(Mb->dlc[lcid].slc[slcid].loe[i].eid == (Top->elements[j].propsurf[0]-1)){
	//Mb->slc[lcid].loe[i].print(out);//tady
	parmefselemload_print(lcid,slcid,i,j,out);
      }
    }
  }


  //  prescribed displacements
  fprintf(out,"\n##   prescribed displacements:");
  fprintf (out,"\n%ld\n\n",Mb->dlc[lcid].slc[slcid].npd);
  for (i=0;i<Mb->dlc[lcid].slc[slcid].npd;i++){
    fprintf (out,"\n %ld",i+1);
    fprintf (out,"  %lf",Mb->dlc[lcid].slc[slcid].pd[i]);
  }
  
  //  prescribed temperature changes
  fprintf(out,"\n##   temperature changes:");
  fprintf(out,"\n%ld\n\n",Mb->dlc[lcid].slc[slcid].tempchang);
  
  if (Mb->dlc[lcid].slc[slcid].tempchang==1 || Mb->dlc[lcid].slc[slcid].tempchang==2){
    // the first property of nodes must contain global node number
    for (i=0;i<Top->nn;i++){
      fprintf(out,"\n%lf",Mb->dlc[lcid].slc[slcid].pt[Top->nodes[i].prop[0]-1]);
    } 
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
void parmefinputfile::parmefelemload_print(long lcid,long leid,long deid,FILE *out)
{
  //  domain element id
  fprintf (out,"\n\n %ld ",deid+1);
  
  Mb->lc[lcid].loe[leid].print(out,0);
}

/**
   function prints loading of elements on one domain for sub-loadcase
   
   @param lcid - loadcase id
   @param leid - load global element eid
   @param deid - domain element eid
   @param out  - output stream

   19/12/2012, TKr
*/
void parmefinputfile::parmefselemload_print(long lcid,long slcid,long leid,long deid,FILE *out)
{
  //  domain element id
  fprintf (out,"\n\n %ld ",deid+1);
  
  Mb->dlc[lcid].slc[slcid].loe[leid].print(out,0);
}



/**
   function returns element type for MEFEL input file
   
   @param etype - element type in jktk mesh format

   19/12/2012, TKr
*/
long parmefinputfile::give_elemtype(long etype)
{
  long type;
  
  switch(etype){
  case 7:{
    type = 100;
    break;
  }
  default:{
    print_err("unknown type of element is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return(type);
}


int parmefinputfile::print_parmefinputfile()
{
  double ts,te;
  FILE *out;
  char *path,*filename,*suffix;
  const char *dn;

  ts = clock();

  dn = Argv[4];
  
  char fname[1020];
  filename_decomposition(Argv[3],path,filename,suffix);

  sprintf(fname, "%smefel_%s.in", path, filename);
  out = fopen(fname,"wt");

  //  input file name
  fprintf(stdout,"\n Printing of PARMEF domain input file no. %s\n",dn);

  //  printing for sequential input file
  //print_in(out); //commented

  //  printing for parallel input file
  print_parmef_in(out);

  fclose(out);
  te = clock();
  printf("\n Time of Printing of PARMEF domain input file %le s\n\n\n",(te-ts)/(double)CLOCKS_PER_SEC);

  delete [] path; delete [] filename; delete [] suffix;
  return(0);  
}



int parmefinputfile::read_topology()
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


int parmefinputfile::run()
{
  double ts,te;
  Top = new siftop;
 
  if (Argc < 5){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name [-kwd=0] domain_input_topology_file_name domain_number\n\n",Argv[0]);
    delete Gtm; delete Top;
    return(1);
  }
  
  fprintf(stdout,"\n\n\n Reading of MEFEL main (sequential) input file . . .\n");
  ts = clock();
  //  reading of mefel input file = program initiation and data reading
  mefel_init (Argc-2,Argv);
  te = clock();
  fprintf(stdout,"\n Time of reading MEFEL main (sequential) input file %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);

  //  reading domain topology
  //  this is for paralel version
  d.paral=1;
  d.redgn = 1;
  err=read_topology();

  //  printing of domain PARMEF input file
  err=print_parmefinputfile();

  return(0);
}


int main (int argc,const char *argv[])
{
  if (argc < 5){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name [-kwd=0] domain_input_topology_file_name domain_number\n\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  parmefinputfile *parmeffile;
  double tws,twe;
  struct tm *local;
  time_t t;
  
  set_prgname(argv[0]);
  printf("---------------------------------------------------------------\n");
  printf("*** PARMEF INPUT FILE GENERATION FROM SEQUENTIAL INPUT FILE ***\n");
  printf("---------------------------------------------------------------\n");
  t = time(NULL);
  local = localtime(&t);
  printf("Program started  %s\n", asctime(local));
  
  tws = clock();
  parmeffile = new parmefinputfile(argc,argv);
  parmeffile->run();
  parmeffile->~parmefinputfile();
  twe = clock();
  printf(" -----------------------------------------------------\n");
  printf ("Whole elapsed time of input file generation %le s\n\n\n",(twe-tws)/(double)CLOCKS_PER_SEC);
}
 
