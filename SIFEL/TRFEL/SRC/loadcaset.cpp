#include "loadcaset.h"
#include "globalt.h"
#include "aliast.h"
#include "globmatt.h"
#include "genfile.h"
#include "elemswitcht.h"
#include <math.h>
#include <limits.h>

loadcaset::loadcaset (void)
{
  //  number of element with boundary conditions
  neb = 0;

  //  number of quantity sources
  nqs = 0;
  //  number of nodes with defined source
  nnqs = 0;
  //  number of nodes with defined point source
  nnpqs = 0;
  //  list of nodes with defined source
  lnqs = NULL;
  //  list of nodes with defined point source
  lnpqs = NULL;
  //  number of elements with defined source
  neqs = 0;
  //  list of elements with defined source
  leqs = NULL;
  //  number of elements influenced by nodes with prescribed source
  neins=0;
  //  class contains map between selected nodes and elements
  en = NULL;
  //  source preprocessing
  sip = 0;
  //  array containing nodes in clusters
  lnc = NULL;
  //  number of nodes in clusters
  nc = NULL;
  //number of elements
  ne = 0;
  //  array of quantity sources
  sour = NULL;
  
  //  number of objects describing nodal values
  nnv=0;
  //  nodal values of variable fluxes or external values
  nodval=NULL;
  //  number of objects describing climatic conditions
  ncc=0;
  ncc2=0;
  //  detailed climatic conditions
  climcond=NULL;
  climcond2=NULL;
  
  
  
  //  number of prescribed values
  npv = 0;
  //  array of prescribed values
  pv = NULL;
  
  elemload = NULL;
  climcond = NULL;
  climcond2 = NULL;
  tabf = NULL;
  
  
  masterval = 0.0;
  mastergrad = NULL;
  ndim = 0;
}

loadcaset::~loadcaset (void)
{
  long i;
  
  //  list of nodes with defined source
  for (i=0;i<nnqs;i++){
    delete [] lnqs[i];
  }
  delete [] lnqs;

  //  list of nodes with defined point source
  for (i=0;i<nnpqs;i++){
    delete [] lnpqs[i];
  }
  delete [] lnpqs;
  
  //  list of elements with defined source
  for (i=0;i<neqs;i++){
    delete [] leqs[i];
  }
  delete [] leqs;
  
  //  map between selected nodes and elements
  delete en;
  
  //  array containing nodes in clusters
  if(lnc != NULL){
    for (i=0;i<nqs;i++){
      delete [] lnc[i];
    }
    delete [] lnc;
  } 
  
  //  number of nodes in clusters
  delete [] nc;
  
  //  array of quantity sources
  delete [] sour;
  
  delete [] nodval;

  //  array of prescribed values
  delete [] pv;

  delete [] tabf;
  delete [] elemload;
  delete [] climcond;
  delete [] climcond2;

  delete [] mastergrad;  
}

/**
   function reads load case characteristics
   
   @param in - pointer to input stream
   @param lcid - load case id
   
   JK, 24.7.2001
*/
void loadcaset::read (XFILE *in,long lcid)
{
  long i,j,k,eid,nbo,nnbo;
  XFILE *locin;
  char fname[1020],*path,*name,*suffix;
  
  // ********************
  //  prescribed values
  // ********************
  //  number of prescribed values
  xfscanf (in,"%ld",&npv);
  if (Mesprt==1)
    fprintf (stdout,"\n the number of prescribed values  %ld",npv);
  if (npv<0)
    print_err("negative number of prescribed values",__FILE__,__LINE__,__func__);
  
  pv = new pvalt [npv];
  for (i=0;i<npv;i++){
    pv[i].read(in);
  }
  
  // *******************
  //  quantity sources
  // *******************
  //  number of quantity sources
  xfscanf (in,"%ld",&nqs);
  if (Mesprt==1)
    fprintf (stdout,"\n number of prescribed sources %ld",nqs);
  if (nqs<0)
    print_err("negative number of quantity sources",__FILE__,__LINE__,__func__);
  if (nqs>0){
    //  models of source
    sour = new sourcet [nqs];
    for (i=0;i<nqs;i++){
      xfscanf (in,"%ld",&j);
      if (j<1 || j>nqs)
	print_err("wrong number of quantity source",__FILE__,__LINE__,__func__);
      sour[j-1].read (in);
    }
    
    //  list of nodes with defined source
    xfscanf (in,"%ld",&nnqs);
    if (Mesprt==1)
      fprintf (stdout,"\n the number of nodes with source %ld",nnqs);
    
    if (nnqs<0)
      print_err("negative number of nodes with defined source",__FILE__,__LINE__,__func__);

    lnqs = new long* [nnqs];
    for (i=0;i<nnqs;i++){
      lnqs[i] = new long [2];
      xfscanf (in,"%ld %ld",&lnqs[i][0],&lnqs[i][1]);
      lnqs[i][0]--;
      lnqs[i][1]--;
    }
    
    //  list of elements with defined source
    xfscanf (in,"%ld",&neqs);
    if (Mesprt==1)
      fprintf (stdout,"\n the number of elements with source %ld",neqs);

    if (neqs<0)
      print_err("negative number of elements with defined source",__FILE__,__LINE__,__func__);

    leqs = new long* [neqs];
    for (i=0;i<neqs;i++){
      leqs[i] = new long [2];
      xfscanf (in,"%ld %ld",&leqs[i][0],&leqs[i][1]);
      leqs[i][0]--;
      leqs[i][1]--;
    }
    
    //  list of nodes with defined point source
    xfscanf (in,"%ld",&nnpqs);
    if (Mesprt==1)
      fprintf (stdout,"\n the number of point sources %ld",nnpqs);

    if (nnpqs<0)
      print_err("negative number of nodes with defined point source",__FILE__,__LINE__,__func__);

    lnpqs = new long* [nnpqs];
    for (i=0;i<nnpqs;i++){
      lnpqs[i] = new long [2];
      xfscanf (in,"%ld %ld",&lnpqs[i][0],&lnpqs[i][1]);
      lnpqs[i][0]--;
      lnpqs[i][1]--;
    }
  }
  
  // **********************************************
  //  elements with boundary conditions
  // **********************************************
  //  number of elements with boundary conditions
  xfscanf (in,"%ld",&neb);
  if (neb<0)
    print_err("negative number of elements with boundary conditions",__FILE__,__LINE__,__func__);

  elemload = new loadelt [neb];
  if (Mesprt==1)
    fprintf (stdout,"\n number of elements with boundary edges %ld",neb);

  for (i=0;i<neb;i++){
    //  element id
    xfscanf (in,"%ld",&eid);
    eid--;
    //  check of element id
    if (eid > Tt->ne-1){
      print_err("Element number in boundary conditions is greater than total number of elements", __FILE__, __LINE__, __func__);
      abort();
    }
    if (eid < 0){
      print_err("Element number in boundary conditions is less than 0", __FILE__, __LINE__, __func__);
      abort();
    }
    
    //  number of boundary objects (end nodes, edges, surfaces)
    //  number of nodes on each boundary object
    Tt->give_nbobjects (eid,nbo,nnbo);
    
    elemload[i].read (in,lcid,eid,nbo,nnbo);
    elemload[i].bc_indicators (Tt->elements[eid].transi[lcid]);
  }
  
  // ***********************************
  //  nodal values defined on boundary
  // ***********************************
  //  number of objects storing nodal values
  xfscanf (in,"%ld",&nnv);
  if (nnv<0)
    print_err("negative number of objects with nodal values",__FILE__,__LINE__,__func__);

  if (Mesprt==1)
    fprintf (stdout,"\n number of objects storing nodal values %ld",nnv);

  nodval = new bnodvalt [nnv];
  for (i=0;i<nnv;i++){
    nodval[i].read (in);
  }
  
  // *********************
  //  climatic conditions
  // *********************
  //  number of objects storing climatic conditions
  xfscanf (in,"%ld",&ncc);
  if (ncc<0)
    print_err("negative number of objects with nodal values",__FILE__,__LINE__,__func__);

  if (Mesprt==1)
    fprintf (stdout,"\n number of objects of climatic conditions %ld",ncc);

  climcond = new climatcond [ncc];
  if (ncc>0){
    //  type of reading of climatic conditions
    xfscanf (in,"%ld",&trcc);

    if (trcc==1){
      //  climatic conditions are described in the input file
      for (i=0;i<ncc;i++){
	xfscanf (in,"%ld",&j);
	if (j<1 || j>ncc){
	  print_err("index of climate condition record is out of range",__FILE__,__LINE__,__func__);
	}
	climcond[j-1].read(in);
      }
    }
    if (trcc==2){
      //  climatic conditions are described in additional files which are read

      //  generation of file name
      filename_decomposition(in->fname,path,name,suffix);
      sprintf(fname, "%s%s.txt", path, name);
      //  local input file with climatic conditions
      locin = xfopen(fname,"r");
      delete [] path; delete [] name; delete [] suffix;
      
      locin->warning = in->warning;
      locin->kwdmode = in->kwdmode;
      locin->ignorecase = in->ignorecase;
      
      for (i=0;i<ncc;i++){
	xfscanf (locin,"%ld",&j);
	if (j<1 || j>ncc){
	  print_err("index of climate condition record is out of range",__FILE__,__LINE__,__func__);
	}
	climcond[j-1].read (locin);
      }
      
      xfclose (locin);
    }
  }
  
  //  number of objects storing climatic conditions
  xfscanf (in,"%ld",&ncc2);
  if (ncc2<0)
    print_err("negative number of objects with nodal values",__FILE__,__LINE__,__func__);

  if (Mesprt==1)
    fprintf (stdout,"\n number of objects of climatic conditions type 2 %ld",ncc2);

  climcond2 = new climatcond2 [ncc2];
  if (ncc2>0){
    //  type of reading of climatic conditions
    xfscanf (in,"%ld",&trcc2);

    if (trcc2==1){
      //  climatic conditions are described in the input file
      for (i=0;i<ncc2;i++){
	xfscanf (in,"%ld",&j);
	if (j<1 || j>ncc2){
	  print_err("index of climate condition record is out of range",__FILE__,__LINE__,__func__);
	}
	climcond2[j-1].read (in);
      }
    }
    if (trcc2==2){
      //  climatic conditions are described in additional files which are read

      //  generation of file name
      filename_decomposition(in->fname,path,name,suffix);
      sprintf(fname, "%s%s.txt", path, name);
      //  local input file with climatic conditions
      locin = xfopen(fname,"r");
      delete [] path; delete [] name; delete [] suffix;
      
      locin->warning = in->warning;
      locin->kwdmode = in->kwdmode;
      locin->ignorecase = in->ignorecase;
      
      for (i=0;i<ncc2;i++){
	xfscanf (locin,"%ld",&j);
	if (j<1 || j>ncc2){
	  print_err("index of climate condition record is out of range",__FILE__,__LINE__,__func__);
	}
	climcond2[j-1].read (locin);
      }
      
      xfclose (locin);
    }
  }


  



  
  
  //part for multpvalt:
  //  number of time functions
  xfscanf (in,"%ld",&nmultfunc);
 
  if (Mesprt==1)
    fprintf (stdout,"\n number of objects of time functions %ld",nmultfunc);
  
  tabf = new tablefunct [nmultfunc];
  
  for(i=0;i<nmultfunc;i++){
    
    xfscanf (in,"%ld",&k);
    
    if ((k-1) != i){
      //print error
      fprintf (stderr,"\n\n wrong number of prescribed time function loadcaset::read (file %s, line %d).\n",__FILE__,__LINE__);
      exit(0);
    }
    
    xfscanf(in,"%d",(int*)&tfunc);
    switch (tfunc){
      
    case tab:{
      tabf[i].read(in);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown function is required in function pvalt::read (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  // ************************
  //  data for homogenization
  // ************************
  if(Tp->homogt == 1){
    ndim = Tp->gdim;

    mastergrad = new double [ndim];
    
    xfscanf (in,"%lf",&masterval);
    for (i=0;i<ndim;i++)
      xfscanf (in,"%lf",&mastergrad[i]);
  }
}


/**
   function prints load case characteristics
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   TKr, 3.1.2006
*/
void loadcaset::print (FILE *out,long lcid)
{
  long i;
  
  // ********************
  //  prescribed values
  // ********************
  fprintf(out,"\n## dirichlet b.c.:");
  fprintf (out,"\n%ld\n\n",npv);
  for (i=0;i<npv;i++){
    pv[i].print(out);
  }
  
  // *******************
  //  quantity sources
  // *******************
  //  number of quantity sources
  fprintf(out,"\n## sources:");
  fprintf (out,"\n%ld\n\n",nqs);
  if (nqs>0){
    //  models of source
    for (i=0;i<nqs;i++){
      fprintf (out,"\n %ld",i+1);
      sour[i].print (out);
    }
    
    //  list of nodes with defined source
    fprintf (out,"\n\n %ld",nnqs);
    for (i=0;i<nnqs;i++){
      fprintf (out,"\n %ld %ld",lnqs[i][0]+1,lnqs[i][1]+1);
    }
    
    //  list of elements with defined source
    fprintf (out,"\n\n %ld",neqs);
    for (i=0;i<neqs;i++){
      fprintf (out,"\n %ld %ld",leqs[i][0]+1,leqs[i][1]+1);
    }
    
    //  list of nodes with defined point source
    fprintf (out,"\n\n %ld",nnpqs);
    for (i=0;i<nnpqs;i++){
      fprintf (out,"\n %ld %ld",lnpqs[i][0]+1,lnpqs[i][1]+1);
    }
   
  }
  
  // **********************************************
  //  number of elements with boundary conditions
  // **********************************************
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",neb);
  for (i=0;i<neb;i++){
    elemload[i].print (out,lcid);
  }
  
  // ***********************************
  //  nodal values defined on boundary
  // ***********************************
  fprintf(out,"\n## nodal values:");
  fprintf (out,"\n%ld\n\n",nnv);
  for (i=0;i<nnv;i++){
    nodval[i].print (out);
  }

  // **********************
  //  climatic conditions
  // **********************
  fprintf(out,"\n## climatic conditions:");
  fprintf (out,"\n%ld\n",ncc);
  if (ncc>0){
    fprintf (out,"%ld\n",trcc);
    for (i=0;i<ncc;i++){
      fprintf (out,"\n\n %ld",i+1);
      climcond[i].print (out);
    }
  }
  
  fprintf(out,"\n## climatic conditions2:");
  fprintf (out,"\n%ld\n",ncc2);
  if (ncc>0){
    fprintf (out,"%ld\n",trcc);
    for (i=0;i<ncc;i++){
      fprintf (out,"\n\n %ld",i+1);
      climcond2[i].print (out);
    }
  }
  
  
  
  //part for multpvalt:
  //  number of time functions
  fprintf(out,"\n## time functions:");
  fprintf (out,"\n%ld\n\n",nmultfunc);
  
  for(i=0;i<nmultfunc;i++){
    
    fprintf (out,"\n %ld",i+1);
    
    fprintf(out,"\n %d ",(int)tfunc);
    switch (tfunc){
      
    case tab:{
      tabf[i].print(out);
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
    fprintf (out,"\n  %lf",masterval);
    for (i=0;i<ndim;i++)
      fprintf (out,"  %lf",mastergrad[i]);
    fprintf (out,"\n");
  }

}


/**
   function assembles %vector of right hand side

   function does not compute contribution from
   initial and prescribed values because they are computed
   in a special function assemble_init
   function computes contributions from Neumann and Newton
   boundary conditions, it means from prescribed fluxes
   and transmission on boundaries
   function also computes contributions from sources
   
   in the case of the Newton boundary condition, this function
   computes \int \kappa T_{ext}
   
   @param lcid - load case id
   @param rhs - right hand side %vector
   
   24.7.2001, JK, 
   revised 25.6.2005,
   added hanging nodes extension by TKo, 07.2018
*/
void loadcaset::assemble (long lcid,double *rhs)
{
  long i, ndofe, ndofemn, eid;
  vector r, f, lv, alv, nodval;
  ivector cn;
  matrix km;
  
  // *****************************
  //  contributions from sources
  // *****************************
  source_contrib (lcid,rhs);
  
  
  // ****************************
  //  contributions from fluxes
  // ****************************
  for (i=0;i<neb;i++){
    //  element id (number of element in the list of all elements in the problem)
    eid=elemload[i].eid;
    
    if(eid < 0){
      print_err("element number is less than 0", __FILE__, __LINE__, __func__);
      abort();
    }
    
    if (Gtt->leso[eid]==1){
      //  only elements switched on are processed
      
      //  number of DOFs on element
      ndofe=Tt->give_ndofe (eid);
      reallocv (RSTCKVEC(ndofe,lv));
      fillv (0.0,lv);
      
      if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
	elemload[i].updatevalues(lcid,Tt->elements[eid].transi[lcid]);
      }
      
      //  function computes nodal values on element defined
      //  by prescribed fluxes (Neumann boundary condition)
      //  (function is defined in elemswitcht.cpp)
      elem_neumann_vector (lv,lcid,eid,i);
      
      //  function computes nodal values on element defined
      //  by prescribed transmission (Newton boundary condition)
      //  (function is defined in elemswitcht.cpp)
      elem_newton_vector (lv,lcid,eid,i);
      //fprintf (stdout,"\n newton  %le %le %le %le",lv[0],lv[1],lv[2],lv[3]);
      
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtt->give_ndofe (eid);
      //  code numbers
      reallocv (RSTCKIVEC(ndofemn, cn));
      Tt->give_code_numbers (eid, cn.a);

      if (ndofe == ndofemn){
        //  localization of element values to the global vector for elements with regular nodes
        locglob (rhs, lv.a, cn.a, ndofe);
      }
      else{
        //  this case means hanging nodes
	  
        //  the element contains hanging nodes
        //  the element value vector has to be transformed
        reallocv (RSTCKVEC(ndofemn, alv));
        mtxv (*Tt->elements[i].tmat, lv, alv);
        locglob (rhs, alv.a, cn.a, ndofemn);
      }
    }  
  }
  
  
  /*
  //  contributions from volume integrals
  for (i=0;i<Tt->ne;i++){
  if (Gtt->leso[i]==1){
  //  only elements switched on are processed
  
  //  number of DOFs on element
  ndofe=Tt->give_ndofe (i);
  reallocv (RSTCKVEC(ndofe,lv));
  nullv (lv);
  
  //  function is defined in elemswitcht.cpp
  volume_rhs_vector (lv,lcid,i);
  
  
  //  the number of DOFs on element
  //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
  ndofemn = Gtt->give_ndofe (i);
  reallocv(RSTCKIVEC(ndofemn, cn));
  //  code numbers
  Tt->give_code_numbers (i, cn);
  
  if(ndofe == ndofemn)
    //  localization of element values to the global vector
    locglob (rhs, lv.a, cn.a, ndofemn);
  else{
    //  the element contains hanging nodes
    //  the element value vector has to be transformed
    reallocv (RSTCKVEC(ndofemn, alv));
    mtxv (*Tt->elements[i].tmat, lv, alv);
    //  localization of element values to the global vector
    locglob (rhs, lv.a, cn.a, ndofemn);
    }
  }
  }  
  */
}




/**
   function assembles %vector of flux on boundary (right hand side)
   
   in the case of the Newton boundary condition, this function
   computes \int \kappa (T-T_{ext})
   
   @param lcid - load case id
   @param rhs - %vector of right hand side
   @param n - number of components in the rhs (unmber of unknowns in the problem)
   
   TKr, 20.2.2004
*/
void loadcaset::assemble_flux (long lcid,double *rhs,long n)
{
  long i, ndofe, ndofemn, eid;
  vector r, f, lv, alv, nodval;
  ivector cn;
  matrix km;
  
  //  all components are set to zero
  nullv (rhs,n);
  
  // *****************************
  //  contributions from sources
  // *****************************
  source_contrib (lcid,rhs);
  
  // ****************************
  //  contributions from fluxes
  // ****************************
  for (i=0;i<neb;i++){
    //  element id (number of element in the list of all elements in the problem)
    eid=elemload[i].eid;
    
    if(eid < 0){
      print_err("element number is less than 0", __FILE__, __LINE__, __func__);
      abort();
    }
    
    if (Gtt->leso[eid]==1){
      //  only elements switched on are processed
      
      //  number of DOFs on element
      ndofe=Tt->give_ndofe (eid);
      reallocv (RSTCKVEC(ndofe, lv));
      nullv(lv);
      
      if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
	elemload[i].updatevalues(lcid,Tt->elements[eid].transi[lcid]);
      }
      
      //  function computes nodal values of boundary fluxes on element
      //  caused by prescribed fluxes (Neumann boundary condition)
      //  (function is defined in elemswitcht.cpp)
      elem_neumann_vector (lv,lcid,eid,i);
      
      //  function computes nodal values of boundary fluxes on element
      //  caused by prescribed transmission (Newton boundary condition)
      //  (function is defined in elemswitcht.cpp)
      elem_transmission_flux (lv,lcid,eid,i);
      
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtt->give_ndofe (eid);
      //  code numbers
      reallocv (RSTCKIVEC(ndofemn, cn));
      Tt->give_code_numbers (eid, cn.a);

      if (ndofe == ndofemn){
        //  localization of element values to the global vector for elements with regular nodes
        locglob (rhs, lv.a, cn.a, ndofe);     
      }
      else{
        //  this case means hanging nodes
        //  the element contains hanging nodes
        //  the element value vector has to be transformed
        reallocv (RSTCKVEC(ndofemn, alv));
        mtxv (*Tt->elements[i].tmat, lv, alv);
        locglob (rhs, alv.a, cn.a, ndofemn);
      }
    }
  }
}


double loadcaset::give_fact (double t,long i)
{
  return tabf[i-1].getval(t);
}


/**
   function searches elements with nodes where source of quantity is defined
   
   JK, 25.6.2005
*/
void loadcaset::elemsource ()
{
  long i,j;
  long *aux;
  
  //  auxiliary array
  aux = new long [nnqs];
  for (i=0;i<nnqs;i++){
    aux[i]=lnqs[i][0];
  }
  
  en = new elemnode;
  //  definition of selected nodes
  en->selnode (nnqs,aux);
  
  delete [] aux;
  
  //  construction of map between selected nodes and elements
  en->elemnodes (Gtt);
  
  //  number of elements influenced by nodes with prescribed source
  neins = en->nie;
  
  for (i=0;i<neins;i++){
    j=en->lse[i];
    Tt->elements[j].source=1;
  }
  
  for (i=0;i<neqs;i++){
    j=leqs[i][0];
    Tt->elements[j].source=1;
  }
  

  
  
  /*
  long i,j,k,m,kk,prev,min,e,nn,nne;
  long **alnc;
  ivector nodes;
  
  //  number of all nodes
  nn = Tt->nn;
  //  number of all elements
  ne = Tt->ne;
  
  if (nn == nnqs){
    // *********************************************
    //  source of quantity is defined at all nodes
    // *********************************************
    
    //  check of ordering
    k=0;
    j=lnqs[0][0];
    for (i=1;i<nn;i++){
      if (lnqs[i][0]>j) k++;
      j=lnqs[i][0];
    }
    if (k!=nn-1){
      fprintf (stderr,"\n\n suspicious list of nodes with sources,\n");
      abort ();
    }
    if (k==nn-1){
      for (i=0;i<ne;i++){
	Tt->elements[i].source=1;
      }
    }
    
    
    //  allocation of list of elements with nodes where source is defined
    lenqs = new long* [ne];
    for (i=0;i<ne;i++){
      //  number of nodes on element
      nne = Tt->give_nne (i);
      //  allocation of list of elements with nodes with nodes where source is defined
      lenqs[i] = new long [nne];
      //  allocation of memory
      allocv (nne,nodes);
      //  list of nodes on element
      Tt->give_elemnodes (i,nodes);
      for (j=0;j<nne;j++){
	lenqs[i][j]=nodes[j];
      }
      //  deallocation of memory
      destrv (nodes);
    }
    
  }
  else{
    // *******************************************************
    //  source of quantity is defined only at selected nodes
    // *******************************************************
    
    //  allocation of list of elements with nodes where source is defined
    lenqs = new long* [ne];
    
    for (i=0;i<ne;i++){
      //  number of nodes on element
      nne = Tt->give_nne (i);
      //  allocation of list of elements with nodes with nodes where source is defined
      lenqs[i] = new long [nne];
      //  allocation of memory
      allocv (nne,nodes);
      //  list of nodes on element
      Tt->give_elemnodes (i,nodes);
      
      for (j=0;j<nne;j++){
	lenqs[i][j]=-1;
	for (k=0;k<nnqs;k++){
	  if (nodes[j]==lnqs[k][0]){
	    lenqs[i][j]=k;
	    Tt->elements[i].source=1;
	    break;
	  }
	}
      }
      
      //  deallocation of memory
      destrv (nodes);
    }
  }
  
  if (sip==1){
    nc = new long [nqs];
    for (i=0;i<nqs;i++){
      nc[i]=0;
    }
    
    for (i=0;i<ne;i++){
      //  number of nodes on element
      nne = Tt->give_nne (i);
      for (j=0;j<nne;j++){
	if (lenqs[i][j]>-1){
	  k=lnqs[lenqs[i][j]][1];
	  nc[k]++;
	}
      }
    }
    
    
    //  list of nodes in clusters
    alnc = new long* [nqs];
    for (i=0;i<nqs;i++){
      alnc[i] = new long [nc[i]];
      nc[i]=0;
    }
    for (i=0;i<ne;i++){
      //  number of nodes on element
      nne = Tt->give_nne (i);
      for (j=0;j<nne;j++){
	if (lenqs[i][j]>-1){
	  k=lnqs[lenqs[i][j]][1];
	  m=lnqs[lenqs[i][j]][0];
	  alnc[k][nc[k]]=m;
	  nc[k]++;
	}
      }
    }
    
    
    //  sorting
    for (i=0;i<nqs;i++){
      prev=LONG_MAX;
      for (j=0;j<nc[i];j++){
	min=LONG_MAX;
	for (k=j;k<nc[i];k++){
	  if (alnc[i][k]>-1){
	    if (min>alnc[i][k]){
	      min=alnc[i][k];
	      kk=k;
	    }
	  }
	}
	if (min==prev){
	  alnc[i][kk]=alnc[i][nc[i]-1];
	  nc[i]--;
	  j--;
	}
	else{
	  e=alnc[i][j];
	  alnc[i][j]=min;
	  alnc[i][kk]=e;
	}
	prev=min;
      }
    }
    
    lnc = new long* [nqs];
    for (i=0;i<nqs;i++){
      lnc[i] = new long [nc[i]];
    }
    
    for (i=0;i<nqs;i++){
      for (j=0;j<nc[i];j++){
	lnc[i][j]=alnc[i][j];
      }
    }
    
    for (i=0;i<nqs;i++){
      delete [] alnc[i];
    }
    delete [] alnc;
  }
  
  */
}

/**
   function assembles source nodal values
   
   @param idse - id of selected element - index in array of selected elements
   @param nodval - array of nodal values of source
   
   JK, 26.6.2005
*/
void loadcaset::sourcenodalvalues (long idse,vector &nodval)
{
  long i,j,k,m,nne;
  
  //  number of nodes on element
  nne = nodval.n;
  
  //  contributions from nodes
  for (i=0;i<nne;i++){
    //  index in array of selected nodes - number of selected node
    j=en->elnod[idse][i];
    if (j>-1){
      //  source id
      k=lnqs[j][1];
      //  element id
      m=en->lse[idse];
      nodval[i]=sour[k].giveval (m);
    }
    else{
      nodval[i]=0.0;
    }
  }
  
}

/**
   function assembles contributions from sources
   
   @param lcid - load case id
   @param rhs - pointer to %vector
   
   JK, 20.10.2007
*/
void loadcaset::source_contrib (long lcid,double *rhs)
{
  long i, j, k, nne, ndofe, ndofemn;
  double s;
  ivector cn;
  vector lv, alv, nodval;

  //  source is defined at nodes
  if (nnqs>0){
    for (i=0;i<neins;i++){
      //  loop over elements which are influenced by prescribed source
      
      //  element id
      j=en->lse[i];
      if (Gtt->leso[j]==1){
	//  only elements switched on are performed
	
	//  number of DOFs on element
	ndofe=Tt->give_ndofe (j);
	//  number of nodes on element
	nne=Tt->give_nne (j);
	reallocv (RSTCKVEC(nne,nodval));
	nullv (nodval);
	reallocv (RSTCKVEC(ndofe,lv));
	nullv (lv);
	
	//  nodal values of sources on element
	sourcenodalvalues (i,nodval);
	
	//  computed nodal fluxes
	source_vector (lcid,j,nodval,lv);
	
        //  the number of DOFs on element
        //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
        ndofemn = Gtt->give_ndofe (j);
        //  code numbers
        reallocv (RSTCKIVEC(ndofemn, cn));
        Tt->give_code_numbers (j, cn.a);

        if (ndofe == ndofemn){
          //  localization of element values to the global vector for elements with regular nodes
          locglob (rhs, lv.a, cn.a, ndofe);     
        }
        else{
          //  this case means hanging nodes
          //  the element contains hanging nodes
          //  the element value vector has to be transformed
          reallocv (RSTCKVEC(ndofemn, alv));
          mtxv (*Tt->elements[i].tmat, lv, alv);
          locglob (rhs, alv.a, cn.a, ndofemn);
        }	
      }
    }   
  }
  
  if (neqs>0){
    //  source is defined at elements
    for (i=0;i<neqs;i++){
      //  element id
      j=leqs[i][0];
      
      if (Gtt->leso[j]==1){
	//  only elements switched on are performed
	
	//  number of DOFs on element
	ndofe=Tt->give_ndofe (j);
	//  number of nodes on element
	nne=Tt->give_nne (j);
	reallocv (RSTCKVEC(nne,nodval));
	nullv (nodval);
	reallocv (RSTCKVEC(ndofe,lv));
	nullv (lv);
	
	//  source id
	k = leqs[i][1];
	s = sour[k].giveval (j);
	for (k=0;k<nne;k++){
	  nodval[k]=s;
	}
	
	//  computed nodal fluxes
	source_vector (lcid,j,nodval,lv);
	
        //  the number of DOFs on element
        //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
        ndofemn = Gtt->give_ndofe (j);
        //  code numbers
        reallocv (RSTCKIVEC(ndofemn, cn));
        Tt->give_code_numbers (j, cn.a);

        if (ndofe == ndofemn){
          //  localization of element values to the global vector for elements with regular nodes
          locglob (rhs, lv.a, cn.a, ndofe);     
        }
        else{
          //  this case means hanging nodes
          //  the element contains hanging nodes
          //  the element value vector has to be transformed
          reallocv (RSTCKVEC(ndofemn, alv));
          mtxv (*Tt->elements[i].tmat, lv, alv);
          locglob (rhs, alv.a, cn.a, ndofemn);
        }
      }
    }    
  }
  
  if (nnpqs>0){
    long ndofn,nid,sid,*cnn;
    double s;
    
    cnn = new long [1];
    for (i=0;i<nnpqs;i++){
      //  loop over nodes with prescribed point source
      
      //  node id
      nid=lnpqs[i][0];
      // number of DOFs in the node
      ndofn = Tt->give_ndofn (nid);
      
      if (ndofn!=1){
        print_err("the number of DOFs in nodes in point sources is not 1", __FILE__, __LINE__, __func__);
        abort ();
      }
      
      //  DOF id
      Tt->give_node_code_numbers (nid,cnn);
      
      //  point source id
      sid = lnpqs[i][1];
      
      //  there should be element id as an argument of the function giveval
      //  in the case of source described by a mathematical function, the argument
      //  is not used inside
      s = sour[sid].giveval (0);
      
      if (cnn[0]>0){
        rhs[cnn[0]-1]+=s;
      }

    }
    
    delete [] cnn;
  }
}

