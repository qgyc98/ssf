#include "loadcase.h"
#include "tablefunct.h"
#include "ipmap.h"
#include "iotools.h"
#include "gtopology.h"
#include "matrix.h"
#include "mathem.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "globmat.h"
#include "elemswitch.h"
#include "loadn.h"
#include "loadel.h"
#include "element.h"
#include "node.h"
#include "elemhead.h"
#include "intpoints.h"
#include <math.h>
#include <stdio.h>

/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
loadcase::loadcase (void)
{
  //  number of loaded nodes
  nln = 0;
  //  number of loaded elements
  nle = 0;
  //  number of prescribed displacements
  npd = 0;
  //  number of prescribed temperature changes
  npt = 0;
  tempchang=0;
  
  //  loaded nodes
  lon = NULL;
  //  loaded elements
  loe = NULL;
  //  prescribed displacements
  pd = NULL;
  //  prescribed temperature changes
  pt = NULL;

  // the total number of prescribed macro-values (strains/stresses), it is dimension of below arrays for homogenization
  ncompmv = 0;  
  // array of  component types of prescribed macro-values
  mstrastre = NULL;
  // prescribed macrostress values for homogeniztaion problems (Mp->homog == 3,5,7,9)
  mstress = NULL;
  // the number of prescribed macro-stress components
  nmstrecomp = 0;
  // array of code numbers of prescribed macro-stress components
  mstress_cn = NULL;
  // prescribed macrostrain values for homogeniztaion problems (Mp->homog == 4,6,8,9)
  mstrain = NULL;
  // the number of prescribed macro-strain components
  nmstracomp = 0;
  //  data file for reading
  memset(temp_file, 0, sizeof(*temp_file)*FNAMELEN);
  //  prescribed temperature changes as multiple values according to time
  pts = NULL;
  // time of individual temperature settings (file)
  temp_time = NULL;
}



/**
  Destructor releases allocated memory of the loadcase object.

  Created by JK,
*/
loadcase::~loadcase (void)
{
  long i;

    delete [] lon;
    delete [] loe;
    delete [] pd;
    delete [] pt;

    delete [] mstrastre;
    delete [] mstress;
    delete [] mstress_cn;
    delete [] mstrain;
  
  if(tempchang == 5){
    if (pts!=NULL){
      for(i=0;i<Mt->nn;i++)
        delete pts[i];
    }
  }
  delete [] pts;
}



/**
  Function reads characteristics of static load case from the opened 
  text file.

  @param in - pointer to the opened XFILE
   
  @return The function does not return anything.

  Created by JK, 24.7.2001
*/
void loadcase::read (XFILE *in)
{
  long i,j;
  XFILE *in_data;
  
  //  loaded nodes
  xfscanf (in,"%ld",&nln);
  fprintf (stdout,"\n the number of loaded nodes %ld",nln);
  lon = new loadn [nln];
  for (i=0;i<nln;i++){
    lon[i].read (in);
  }
  
  //  loaded elements
  xfscanf (in,"%ld",&nle);
  fprintf (stdout,"\n the number of loaded elements %ld",nle);
  loe = new loadel [nle];
  for (i=0;i<nle;i++){
    loe[i].read (in);
  }
  
  //  prescribed displacements
  xfscanf (in,"%ld",&npd);
  fprintf (stdout,"\n the number of prescribed displacements %ld",npd);
  pd = new double [npd];
  for (i=0;i<npd;i++){
    xfscanf (in,"%lf",pd+i);
  }
  
  //  prescribed temperature changes
  xfscanf (in,"%ld",&tempchang);
  if (tempchang==0)
    fprintf (stdout,"\n no temperature changes are described");
  switch (tempchang)
    {
    case 0:
      break;
    case 1: // temperature is defined in MEFEL input file
    case 2: // scaled temperature is defined in MEFEL input file
    case 4:
    case 5:
      if (Mp->temperature < 2)
        Mp->temperature=1;
      else
      {
        print_err("type 1 or 2 of temperature load cannot be used together with type 3", __FILE__, __LINE__, __func__);
        abort();
      }
      break;
    case 3: // temperature is calculated by TRFEL
      if (Mp->temperature == 1)
      {
        print_err("type 1 or 2 of temperature load cannot be used together with type 3", __FILE__, __LINE__, __func__);
        abort();
      }
      else
        Mp->temperature=2;
      break;
    default:
      print_err("unknown type of temperature load (=%ld) is required", __FILE__, __LINE__, __func__, tempchang);
      abort();
  }

  if (tempchang==1 || tempchang==2){
    fprintf (stdout,"\n temperature loading from input file");
    npt = Mt->nn;
    pt = new double [npt];
    for (i=0;i<npt;i++){
      xfscanf (in,"%lf",pt+i);
    }
  }

  if (tempchang==3)
    fprintf (stdout,"\n temperatures are computing in TRFEL");

  if (tempchang==4){//reading temperatures from file
    fprintf (stdout,"\n temperature loading - reading temperatures from extra input file");
    xfscanf(in, " %a", temp_file);
    
    in_data = xfopen(temp_file,"r");
    
    npt = Mt->nn;
    pt = new double [npt];
    for (i=0;i<npt;i++){
      xfscanf (in_data,"%lf",pt+i);
    }
    xfclose(in_data);
  }
  
  if (tempchang==5){//reading temperatures from multiple files
    fprintf (stdout,"\n temperature loading - reading temperatures from extra multiple input files");
    xfscanf(in, " %ld",&ntemp_file);
    npt = Mt->nn;
    pts = new double*[npt];
    for (i=0; i<npt; i++)
      pts[i] = new double[ntemp_file];
    temp_time = new double[ntemp_file];
    for (i=0;i<ntemp_file;i++){
      xfscanf(in, " %lf",temp_time+i);
      xfscanf(in, " %a", temp_file);
      in_data = xfopen(temp_file,"r");
      for (j=0;j<npt;j++){
	xfscanf (in_data,"%lf",&pts[j][i]);
      }
      xfclose(in_data);
    }
  }

  if ((Mp->homog == 3) || (Mp->homog == 5) || (Mp->homog == 7))
  {
    ncompmv = nmstrecomp = Mt->max_ncompstr;
    mstress    = new double[nmstrecomp];    
    mstrastre  = new strastre[nmstrecomp];
    mstress_cn = new long[nmstrecomp];
    for (i=0; i<nmstrecomp; i++){
      mstrastre[i] = stress;
      if (Mp->homog == 3)
        xfscanf(in, "%le", mstress+i);
      else
        mstress[i] = 0.0;
      mstress_cn[i] = Ndofm+i+1;
    }
  }
  if ((Mp->homog == 4) || (Mp->homog == 6) || (Mp->homog == 8))
  {
    ncompmv = nmstracomp = Mt->max_ncompstr;
    mstrain   = new double[nmstracomp];
    mstrastre = new strastre[nmstracomp];
    mstress_cn = new long[nmstracomp];
    memset(mstress_cn, 0, sizeof(*mstress_cn)*nmstracomp);
    for (i=0; i<nmstracomp; i++){
      mstrastre[i] = strain;
      if (Mp->homog == 4)
        xfscanf(in, "%le", mstrain+i);
      else
        mstrain[i] = 0.0;
    }
  }
  if (Mp->homog == 9){
    ncompmv = Mt->max_ncompstr;
    
    mstrastre = new strastre[ncompmv];
    memset(mstrastre, 0, sizeof(*mstrastre)*ncompmv);

    mstress = new double[ncompmv];
    memset(mstress, 0, sizeof(*mstress)*ncompmv);

    mstress_cn = new long[ncompmv];
    memset(mstress_cn, 0, sizeof(*mstress_cn)*ncompmv);

    mstrain = new double[ncompmv];
    memset(mstrain, 0, sizeof(*mstrain)*ncompmv);

    for (i=0; i<ncompmv; i++){
      xfscanf(in, "%m", &strastre_kwdset, (int*)(mstrastre+i));
      switch (mstrastre[i]){
        case strain:
          xfscanf(in, "%le", mstrain+i);
          nmstracomp++;
          break;
        case stress:
          xfscanf(in, "%le", mstress+i);
          nmstrecomp++;
          mstress_cn[i] = Ndofm + nmstrecomp;
          break;
        default:
          print_err("unknown type (%d) of macro value is required for the %ld component",
                    __FILE__, __LINE__, __func__, mstrastre[i], i+1);
          abort();
      }
    }
  }
}



/**
  Function prints characteristics of static load case into the opened 
  text file.

  @param out - pointer to the opened FILE
   
  @return The function does not return anything.

  TKr, 07/02/2013 accroding to read (XFILE *in)
*/
void loadcase::print (FILE *out)
{
  long i, j;
  
  //  loaded nodes
  fprintf(out,"\n## loaded nodes:");
  fprintf (out,"\n%ld\n\n",nln);
  for (i=0;i<nln;i++){
    lon[i].print (out);
  }
  
  //  loaded elements
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",nle);
  for (i=0;i<nle;i++){
    loe[i].print (out,0);
  }
  
  //  prescribed displacements
  fprintf(out,"\n##   prescribed displacements:");
  fprintf (out,"\n%ld\n\n",npd);
  for (i=0;i<npd;i++){
    fprintf (out,"\n %le",pd[i]);
  }
  
  //  prescribed temperature changes
  fprintf(out,"\n##   temperature changes:");
  fprintf(out,"\n%ld\n\n",tempchang);
  
  if (tempchang==1 || tempchang==2){
    for (i=0;i<npt;i++){
      fprintf (out,"\n %le",pt[i]);
    }
  }

  if (tempchang==4){//toto jeste opravit TKr 20/05/2014
    fprintf(out, " %s", temp_file);
    fprintf(out, "\n");
  }

  if (tempchang==5){//toto jeste opravit TKr 20/05/2014
  }
  
  if ((Mp->homog == 3) || (Mp->homog == 5) || (Mp->homog == 7))
  {
    fprintf(out,"\n##   macro-stresses:");
    for (i=0; i<ncompmv; i++)
      fprintf(out, "%le ", mstress[i]);
    fprintf(out, "\n");
  }
  if ((Mp->homog == 4) || (Mp->homog == 6) || (Mp->homog == 8))
  {
    fprintf(out,"\n##   macro-strains:");
    for (i=0; i<ncompmv; i++)
      fprintf(out, "%le ", mstrain[i]);
    fprintf(out, "\n");
  }
  if (Mp->homog == 9){
    fprintf(out,"\n##   hybrid macro-stress/macro-strains:");

    for (i=0; i<ncompmv; i++){
      j = 0;
      if (mstrastre[i] == stress){
        fprintf(out, "%d %le\n", (int)(stress), mstress[i]);
        j++;
      }
      if (mstrastre[i] == strain){
        fprintf(out, "%d %le\n", (int)(strain), mstrain[i]);
        j++;
      }
      if (j == 0){
        print_err("no macro-stress nor macro-strain value was prescribed for the %ld-th component",
                  __FILE__, __LINE__, __func__, i+1);
        abort();
      }
      if (j == 2){
        print_err("both macro-stress and macro-strain value cannot be prescribed for the %ld-th component",
                  __FILE__, __LINE__, __func__, i+1);
        abort();
      }
    }
  }
}



/**
  Function assembles %vector of right hand side. 
  It is created by prescribed displacements, forces, moments, macrostrains, etc.   
  In case that the flv is not NULL, %vector of load caused by forces only is stored in flv.

  @param lcid  - load case id
  @param rhs   - pointer to the right hand side
  @param flv   - array of load %vector caused by forces only (default value is NULL)
  @param scale - scale factor, it is used for subload cases
   
  Created by JK, 24.7.2001
  Modifed by Tomas Koudelka, 7.2008
  Rewritten by Tomas Koudelka, 11.2015
*/
void loadcase::assemble (long lcid,double *rhs,double *flv,double scale)
{
  //  array rhs and flv have to be cleaned in the caller (above) function

  // Contributions of all kinds of load except of prescribed displacements are assembled by
  // the following function call.
  // Contributions due to prescribed displacements must be handled separately due to time 
  // dependent problems where the loadcase class is used for time dependent subloadcases.
  // In these problems, PD contributions would be counted repeatedly due to definition of 
  // dloadcase::get_pd which returns sum of all prescribed displacements for all subloadcases.
  assemblewopd (lcid, rhs, flv, scale);

  //  contributions from prescribed displacements
  if (npd>0)
    prdisplcontrib (lcid, rhs);
}



/**
  Function assembles %vector of right hand side WITHOUT contributions due to prescribed displacements. 
  It is created by forces, moments, temperature changes, macrostrains, etc.   
  In case that the flv is not NULL, %vector of load caused by forces only is stored in flv.

  @param lcid  - load case id
  @param rhs   - pointer to the right hand side
  @param flv   - array of load %vector caused by forces only (default value is NULL)
  @param scale - scale factor, it is used for subload cases
   
  Created by JK, 24.7.2001
  Modifed by Tomas Koudelka, 7.2008
  Rewritten by Tomas Koudelka, 11.2015
*/
void loadcase::assemblewopd (long lcid,double *rhs,double *flv,double scale)
{
  long i, j, ncomp, tnip, aux;
  
  //  contributions due to effects of forces
  forcecontrib (lcid, rhs, scale);
  if (flv)
    copyv(rhs, flv, Ndofm);

  //  contributions from temperature changes
  tempercontrib (lcid,rhs,scale);

  // contribution from macro stresses in the homogenization problem
  if (nmstrecomp)
  {
    ncomp = Mt->max_ncompstr;
    for (i=0; i<ncomp; i++)
    {
      if (mstrastre[i] == stress){
        j = mstress_cn[i]-1;
        rhs[j] += scale*mstress[i]*Mt->domvol;
        if (flv)
          flv[j] += scale*mstress[i]*Mt->domvol;
      }
    }
    if (nmstracomp){ // some of macro-strains have been prescribed
      matrix dd(ASTCKMAT(ncomp, ncomp)), gdd(ASTCKMAT(ncomp, ncomp));
      vector ifor(ASTCKVEC(ncomp)), mstra;
      for (i=0; i<Mt->ne; i++){
        ddmatrix(lcid,i,dd); // \int_{\Omega} \mathbf{D} \mathrm{d}V
        addm(gdd, dd, gdd);
      }
      makerefv(mstra, mstrain, ncomp);
      mxv(gdd, mstra, ifor);
      for (i=0; i<ncomp; i++){
        if (mstrastre[i] == stress){
          j = mstress_cn[i]-1;
          rhs[j] -= ifor(i);
          /*          if (flv)
                      flv[j] -= ifor(i);*/
        }
      }
    }
  }
  // contribution from macro strain in the homogenization problem
  if (nmstracomp)
  {
    //  the number of strain components
    ncomp = Mt->max_ncompstr;
    matrix d(ASTCKMAT(ncomp, ncomp));
    vector eps, sig, ifor;
    ivector cn;
    int err;
    //  the number of all integration points
    tnip = Mm->tnip;
    
    for (i=0;i<tnip;i++){
      for (j=0;j<ncomp;j++){
        Mm->ip[i].strain[j]=0.0;
        //  negative strain is stored because the contribution to
        //  the right hand side should be negative, function
        //  internal_forces adds contributions with positive sign
        if (mstrastre[j] == strain)
          Mm->ip[i].strain[j] -= scale*mstrain[j];
      }
      makerefv(eps, Mm->ip[i].strain, ncomp);
      makerefv(sig, Mm->ip[i].stress, ncomp);
      Mm->matstiff(d, i);
      mxv(d, eps, sig);
      /*      for (j=0;j<ncomp;j++){
        if (mstrastre[j] == stress)
          sig(j) = 0.0;
          }*/
    }

    aux = Mp->strcomp;    
    Mp->strcomp = 0;
    long ndofe, ndofemn;    
    for (i=0; i<Mt->ne; i++){
      ndofe=Mt->give_ndofe (i);
      reallocv (RSTCKVEC(ndofe,ifor));
	
      elem_internal_forces (i,lcid,ifor);	
      err = check_math_errel(i);	
      if (err){
        abort();
      }
      
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe (i);
      reallocv (RSTCKIVEC(ndofemn,cn));	
      Mt->give_code_numbers (i,cn.a);
      locglob (rhs, ifor.a, cn.a, ndofe);
    }
    //internal_forces (0,rhs);
    Mp->strcomp = aux;
  }
}



/**
  Function assembles contribution to the right hand side %vector due to effects of forces. 
  It is created by nodal forces, nodal moments, continuous load on elements, etc.   

  @param lcid  - load case id
  @param rhs   - pointer to the right hand side
  @param scale - scale factor, it is used for subload cases
   
  Created by JK, 24.7.2001
  Modifed by Tomas Koudelka, 7.2008
*/
void loadcase::forcecontrib (long /*lcid*/,double *rhs,double scale)
{
  long i,j,k,eid,ndofn,ndofe,ndofemn;
  ivector cn;
  vector anf;
  
  //  array rhs has to be cleaned in the caller (above) function
  
  // ***************************
  //  contributions from nodes
  // ***************************
  for (i=0;i<nln;i++){
    ndofn=Mt->give_ndofn (lon[i].nid);
    if (ndofn<0){
      print_err("loaded node %ld is a hanging node",__FILE__,__LINE__,__func__,i++);
    }
    for (j=0;j<ndofn;j++){
      k=Mt->give_dof (lon[i].nid,j);
      if (k<=0)  continue;
      if (k>0)  rhs[k-1]+=lon[i].f[j];
    }
  }
  
  // ******************************
  //  contributions from elements
  // ******************************
  //  loop over the loaded elements
  for (i=0;i<nle;i++){
    //  element id
    eid=loe[i].eid;
    if (Gtm->leso[eid]==1){
      
      //  the number of DOFs on element
      ndofe=Mt->give_ndofe (eid);
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe (eid);
      
      reallocv (RSTCKIVEC(ndofemn,cn));
      //  code numbers on element
      Mt->give_code_numbers (eid,cn.a);
      
      if (ndofe != ndofemn){
	//  the element contains hanging node
	//  the vector has to be transformed
	reallocv (RSTCKVEC(ndofemn,anf));
	mtxv (*Mt->elements[i].tmat,loe[i].nf,anf);
	locglob (rhs,anf.a,cn.a,ndofemn);
      }else{
	//  the element does not contain hanging node
	locglob (rhs,loe[i].nf.a,cn.a,ndofe);
      }
    }
  }
  cmulv(scale,rhs,Ndofm);
}

  

/**
  Function assembles contribution to the right hand side %vector due to effects of prescribed 
  displacements at nodes. 

  @param lcid  - load case id
  @param rhs   - pointer to the right hand side
   
  Created by JK, 24.7.2001
  Modifed by Tomas Koudelka, 11.2015
*/
void loadcase::prdisplcontrib (long lcid, double *rhs)
{
  long i,ndofe,ndofemn,ne;
  ivector cn;
  vector r,f,af;
  matrix sm;
  
  //  the number of elements
  ne = Mt->ne;
  for (i=0;i<ne;i++){
    if (Mt->elements[i].prescdispl==1){
      //  only elements with prescribed displacements are taken into account
	
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs without hanging nodes
      ndofe=Mt->give_ndofe (i);
      reallocm (RSTCKMAT(ndofe,ndofe,sm));
      //  stiffness matrix of an element
      stiffmat (lcid,i,sm);
	
      reallocv (RSTCKVEC(ndofe,r));
      reallocv (RSTCKVEC(ndofe,f));
      //  prescribed displacements on element
      elprdispl (lcid,i,r.a);
      //  K*r_0 = f_r
      mxv (sm,r,f);
      // rhs -= f_r 
      cmulv (-1.0,f);
	
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe (i);
      reallocv (RSTCKIVEC(ndofemn,cn));
	
      //  code numbers on element
      Mt->give_code_numbers (i,cn.a);
	
      if (ndofe != ndofemn){
        //  the element contains hanging node
        //  the vector f has to be transformed
        reallocv (RSTCKVEC(ndofemn,af));
        mtxv (*Mt->elements[i].tmat,f,af);
        locglob (rhs,af.a,cn.a,ndofemn);
      }else{
        //  the element does not contain hanging node
        locglob (rhs,f.a,cn.a,ndofe);
      }
    }
  }
}



/**
  Function computes contributions from temperature changes
  to the right hand side (%vector of nodal forces).
   
  @param lcid - load case id
  @param rhs - pointer to the right hand side
  @param scale - scale factor, it is used for subload cases
   
  Created by JK, Tomas Koudelka, 21.11.2007
*/
void loadcase::tempercontrib (long lcid,double *rhs,double scale)
{
  long i;

  //  contributions from temperature changes
  if (tempchang==1 || tempchang==4){
    //  contribution to the time-independent problems from temperatures
    //  read from input files
    
    //  type of eigenstrain is temperature strain
    Mm->est=tempstrain;

    //  approximation of nodal temperatures to integration points
    intpointval (pt,temperature,1.0);
    
    //  thermal material models compute temperature strains
    Mm->temprstrains ();

    //  computation of nodal forces caused by temperature
    nodal_eigstrain_forces (lcid,rhs,Mp->time);
  }
  
  if (tempchang==2){
    //  contribution to the time-dependent problems from temperatures
    //  read from input file, temperature is changed by multiplication of
    //  read function
    
    //  type of eigenstrain is temperature strain
    Mm->est=tempstrain;

    //  approximation of nodal temperatures
    intpointval (pt,temperature,scale);
    
    //  thermal material models compute temperature strains
    Mm->temprstrains ();
    
    nodal_eigstrain_forces (lcid,rhs,Mp->time);
  }

  if (tempchang==3){
    //  contribution to the time-dependent coupled problems from
    //  temperatures obtained from heat transfer (from TRFEL)
    
    //  type of eigenstrain is temperature strain
    Mm->est=tempstrain;

    //  thermal material models compute temperature strains
    Mm->temprstrains ();
    
    //  computation of nodal forces caused by temperature
    nodal_eigstrain_forces (lcid,rhs,Mp->time);
  }
  
  //  contributions from temperature changes reading from multiple input files
  if (tempchang==5){
    //  contribution to the time-independent problems from temperatures
    //  read from input files
    
    //  type of eigenstrain is temperature strain
    Mm->est=tempstrain;

    //selection of correct nodal values
    npt = Mt->nn;
    pt = new double[npt];

    for(i=0;i<Mt->nn;i++){
      tablefunct tabf1(temp_time, pts[i], ntemp_file);
      pt[i] = tabf1.getval(Mp->time);
    }

    //  approximation of nodal temperatures to integration points
    intpointval (pt,temperature,1.0);
    
    //  thermal material models compute temperature strains
    Mm->temprstrains ();

    //  computation of nodal forces caused by temperature
    nodal_eigstrain_forces (lcid,rhs,Mp->time);    
  }
  
}



/**
   Function computes reactions.
   
   @param lcid - load case id
   @param scale - scaling factor for contributions of nodal and element load
   @param comp_intfc - flag for computing internal forces contributions to reactions
                       comp_intfc = yes = 1 => contributions will be computed
                       comp_intfc = no = 0 => contributions will NOT be computed
   
   @return The function does not return anything.
   
   Created by JK,
   Modified by Tomas Koudelka, 27.4.2016 (added nodal load contributions to reactions)
*/
void loadcase::compute_reactions (long lcid, double scale, answertype comp_intfc)
{
  long i, j, k, ii;
  long nne, ndofn, ne, ndofe, ndofemn, nid, eid;
  ivector nod;
  vector f, af;
  
  // *************************************************************
  // Compute contributions of element internal forces to reactions
  // *************************************************************

  if (comp_intfc == yes)
  {
    ne = Mt->ne;  //  the number of elements
    for (i=0; i<ne; i++){  //  loop over the number of elements
      if ((Mt->elements[i].react==1) && (Gtm->leso[i]==1)){
        //  only elements connected to a node with reaction are taken into account
      
        //  the number of DOFs on element
        //  this number is generated by elements without hanging nodes
        ndofe=Mt->give_ndofe(i);
        reallocv(RSTCKVEC(ndofe, f));
        //  nodal forces are computed
        elem_internal_forces(i, lcid, f);
      
        //  the number of DOFs on element
        //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
        ndofemn = Gtm->give_ndofe(i);
      
        //  the number of nodes on element
        //  if there are hanging nodes, the variable will be rewritten
        nne=Mt->give_nne(i);
        
        if ((ndofe == ndofemn) || nmstrecomp){
          // there are no hanging nodes or
          // there is homogenization based on prescribed macro-stress component approach
          reallocv(RSTCKIVEC(nne, nod));
          Mt->give_elemnodes(i, nod);
        }
        else{
          //  the element contains hanging node
          //  the vector of internal forces has to be transformed
          reallocv(RSTCKVEC(ndofemn, af));
          //  transformation of the vector f
          mtxv(*Mt->elements[i].tmat, f, af);
          reallocv(RSTCKVEC(ndofemn, f));
          copyv(af, f);
          //  the number of nodes on element
          //  it is equal to the number of master nodes
          nne=Gtm->gelements[i].nmne;
          reallocv(RSTCKIVEC(nne, nod));
          Gtm->gelements[i].give_master_nodes(nod);
        }

        //  at this moment, the arrays nod and f contain the appropriate number of components
        
        ii=0;
        //  loop over the number of nodes on element
        //  in the case of element with hanging nodes, the correct number of nodes is used
        for (j=0; j<nne; j++){
          ndofn=Mt->give_ndofn(nod(j)); //  the number of DOFs in node
          if (Mt->nodes[nod(j)].react == 1){
            //  only nodes with constraint/reaction are taken into account
            
            //  loop over the number of DOFs in the actual node
            for (k=0; k<ndofn; k++){
              Mt->nodes[nod(j)].r[k] += f[ii];
              ii++;
            }//  end of the loop over the number of DOFs in the actual node
          }
          else  ii+=ndofn;
        }//  end of the loop over the number of nodes on element
      }
    }//  end of the loop over the number of elements
  }
  
  
  // *****************************************************************************
  //  contributions from nodal load applied at nodes with prescribed displacements
  // *****************************************************************************
  // prepared for the moving load

  for (i=0; i<nln; i++){ // loop over the number of loaded nodes
    nid   = lon[i].nid; // node id
    ndofn = Mt->give_ndofn(nid);
    if (ndofn<0)
      print_err("loaded node %ld is a hanging node",__FILE__,__LINE__,__func__,nid+1);

    for (j=0; j<ndofn; j++){      
      k = Mt->give_dof(nid, j);
      if (k<=0) Mt->nodes[nid].r[j] -= scale*lon[i].f[j]; // nodal load is being applied at the direction of prescribed displacement
      if (k>0)  continue;
    }
  }

  // **************************************************
  // Compute contributions of element load to reactions
  // **************************************************

  for (i=0; i<nle; i++){  //  loop over the number of loaded elements
    eid=loe[i].eid;  // element id
    if ((Mt->elements[eid].react==1) && (Gtm->leso[eid]==1)){
      //  only switched on elements with reactions are taken into account
      
      //  the number of DOFs on element
      ndofe = Mt->give_ndofe(eid);
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe(eid);

      if ((ndofe == ndofemn) || nmstrecomp){
        // there are no hanging nodes or there is homogenization 

        //  the number of nodes on element
        nne = Mt->give_nne(eid);
        reallocv(RSTCKIVEC(nne, nod));
        Mt->give_elemnodes(eid, nod);
        reallocv(RSTCKVEC(ndofe, f));
        copyv(loe[i].nf, f);
      }
      else{
        //  the element contains hanging node
        //  the vector of internal forces has to be transformed
	
        reallocv(RSTCKVEC(ndofemn, f));
        //  transformation of the vector f
        mtxv(*Mt->elements[eid].tmat, loe[i].nf, f);
	
        //  the number of nodes on element
        //  it is equal to the number of master nodes
        nne = Gtm->gelements[eid].nmne;
        reallocv(RSTCKIVEC(nne, nod));
        Gtm->gelements[eid].give_master_nodes(nod);
      }
      
      //  at this moment, the arrays nod and f contain the appropriate number of components
      
      ii=0;
      //  loop over the number of nodes on element
      //  if hanging nodes are present, the loop takes them into account
      for (j=0; j<nne; j++){
	//  the number of DOFs in node
	ndofn = Mt->give_ndofn(nod(j));
	if (Mt->nodes[nod(j)].react == 1){
	  //  only nodes where reaction is are taken into account
	  for (k=0; k<ndofn; k++){
	    Mt->nodes[nod(j)].r[k] -= f[ii];
	    ii++;
	  }
	}
	else  ii += ndofn;
      }
    }
  }
}



/**
  Function computes contribution to strains due to temperature changes at auxiliary integration points from the given load case.
   
  @param lcid[in]  - load case id
  @param scale[in] - scale factor, it is used for subload cases
  @param n[in]     - the number of required auxiliary points in the mapping array ipm
  @param ipm[in]   - integration point mapping array, 
                     ipm[i].ipp < 0 => auxiliary integration point must be used
                     ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                        no computation of strains is performed, the strains are assumed 
                                        to be computed at the main solution procedure of the problem
   
  Created Tomas Koudelka, 15.06.2018
*/
void loadcase::aip_cumultemperstrains_contrib (long /*lcid*/, double scale, long n, ipmap *ipm)
{
  long i, j, app, eid, nne;
  ivector enod;
  double aip_pt;
  vector natcoord(ASTCKVEC(3));
  vector nodval;
  
  intpoints *tmp_ip;
  long      *tmp_elip;
  double   **tmp_nonmechq;
  double   **tmp_ic;
  double   **tmp_eigstrains;
  double   **tmp_eigstresses;
  double   **tmp_tempstrains;

  
  // swap regular integration and point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip          = Mm->ip;
  tmp_elip        = Mm->elip;
  tmp_nonmechq    = Mm->nonmechq;
  tmp_ic          = Mm->ic;
  tmp_eigstrains  = Mm->eigstrains;
  tmp_eigstresses = Mm->eigstresses;
  tmp_tempstrains = Mm->tempstrains;

  Mm->tnip         = Mm->tnaip;
  Mm->ip           = Mm->aip;
  Mm->elip         = Mm->elaip;
  Mm->nonmechq     = Mm->aip_nonmechq;
  Mm->ic           = Mm->aip_ic;
  Mm->eigstrains   = Mm->aip_eigstrains;
  Mm->eigstresses  = Mm->aip_eigstresses;
  Mm->tempstrains  = Mm->aip_tempstrains;

  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  continue;
    //  contributions from temperature changes
    if (tempchang==1 || tempchang==4)
    {
      //  contribution to the time-independent problems from temperatures
      //  read either from input file directly (1) or from another file (4)
    
      // get element number for the aux. int. point
      eid = Mm->elaip[app];
      // assemble vector of nodal temperature changes at the eid-th element
      nne = Mt->give_nne(eid);
      reallocv(RSTCKIVEC(nne, enod));
      reallocv(RSTCKVEC(nne, nodval));
      Mt->give_elemnodes(eid, enod);
      for (j=0; j<nne; j++)
        nodval[j] = pt[enod[j]];

      // natural coordinates of the aux. int. point
      natcoord(0) = ipm[i].xi;
      natcoord(1) = ipm[i].eta;
      natcoord(2) = ipm[i].zeta;
      // calculate actual temperature at the app-th aux. int. point
      aip_pt = interpolate(eid, nodval, natcoord);

      Mm->storenonmechq(temperature, app, aip_pt);
      //  thermal material models compute temperature strains
      Mm->cumultemprstrainsmat(app);
    }
  
    if (tempchang==2)
    {
      //  contribution to the time-dependent problems from temperatures
      //  read from input file, temperature is changed by multiplication of
      //  read function
    
      // get element number for the aux. int. point
      eid = Mm->elaip[app];
      // assemble vector of nodal temperature changes at the eid-th element
      nne = Mt->give_nne(eid);
      reallocv(RSTCKIVEC(nne, enod));
      Mt->give_elemnodes(eid, enod);
      for (j=0; j<nne; j++)
        nodval[j] = pt[enod[j]];

      // natural coordinates of the aux. int. point
      natcoord(0) = ipm[i].xi;
      natcoord(1) = ipm[i].eta;
      natcoord(2) = ipm[i].zeta;
      // calculate actual temperature at the app-th aux. int. point
      aip_pt = scale*interpolate(eid, nodval, natcoord);

      Mm->storenonmechq(temperature, app, aip_pt);
      //  thermal material models compute temperature strains
      Mm->cumultemprstrainsmat (app);
    }

    if (tempchang==3)
    {
      //  contribution to the time-dependent coupled problems from
      //  temperatures obtained from heat transfer (from TRFEL)    
      //  in this case, it is assumed that actual values of temparature at aux. int points
      //  have already been updated => no temperature approximation is needed

      // thermal material models compute temperature strains
      Mm->cumultemprstrainsmat (app);
    }
  
    //  contributions from temperature changes reading from multiple input files
    if (tempchang==5)
    {
      //  contribution to the time-independent problems from temperatures
      //  read from input files

      //selection of correct nodal values
      npt = Mt->nn;
      pt = new double [npt];

      for(i=0;i<Mt->nn;i++){
        tablefunct tabf1(temp_time, pts[i], ntemp_file);
        pt[i] = tabf1.getval(Mp->time);
      }

      // get element number for the aux. int. point
      eid = Mm->elaip[app];
      // assemble vector of nodal temperature changes at the eid-th element
      nne = Mt->give_nne(eid);
      reallocv(RSTCKIVEC(nne, enod));
      Mt->give_elemnodes(eid, enod);
      for (j=0; j<nne; j++)
        nodval[j] = pt[enod[j]];
      // natural coordinates of the aux. int. point
      natcoord(0) = ipm[i].xi;
      natcoord(1) = ipm[i].eta;
      natcoord(2) = ipm[i].zeta;

      // calculate actual temperature at the app-th aux. int. point
      aip_pt = interpolate(eid, nodval, natcoord);
      // store temperature at aip_nonmechq array 
      // (arrays nonmechq and aip_nonmechq have been swapped at the beginning of this procedure)
      Mm->storenonmechq(temperature, app, aip_pt);
      //  thermal material models compute temperature strains and ADDS them in aip_tempstrains array
      // (arrays tempstrains and aip_tempstrains have been swapped at the beginning of this procedure)
      Mm->cumultemprstrainsmat (app);
    }
  }  

  // swap regular integration and point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Mm->tnip         = Mm->tnrip;
  Mm->ip           = tmp_ip;
  Mm->elip         = tmp_elip;
  Mm->nonmechq     = tmp_nonmechq;
  Mm->ic           = tmp_ic;
  Mm->eigstrains   = tmp_eigstrains;
  Mm->eigstresses  = tmp_eigstresses;
  Mm->tempstrains  = tmp_tempstrains;
}



/**
  The function returns number of prescribed macro-stress components at the given load case.
  It is intended for the homogenization problems (Mp->homog = 3,5,7,9)

  @return The number of prescribed macro-stress components.

  Created by Tomas Koudelka, 02.2020
*/
long loadcase::give_num_mstress_comp()
{
  return nmstrecomp;
}



/**
  The function returns number of prescribed macro-strain components at the given load case.
  It is intended for the homogenization problems (Mp->homog = 4,6,8,9)

  @return The number of prescribed macro-strain components.

  Created by Tomas Koudelka, 02.2020
*/
long loadcase::give_num_mstrain_comp()
{
  return nmstracomp;
}



/**
  The function returns pointer to the array of component types of prescribed macro-value components 
  at the given load case.  It is intended for the homogenization problems (Mp->homog > 2)

  @return The function returns pointer to the answertype array with indicators.

  Created by Tomas Koudelka, 02.2020
*/
strastre* loadcase::give_mstrastre()
{
  return mstrastre;
}



/**
  The function returns pointer to the array of code numbers of macro-stress components 
  at the given load case. It is intended for the homogenization problems (Mp->homog > 2)

  @return The pointer to the array of macro-stress component code (DOF) numbers,
          dimension of the array is the maximum total number of stress/strain components of used elements.
          Particular components of the array contains either positive nonzero value which represents the code (DOF) number
          or zero value for components where no macro-stress value has been prescribed.

  Created by Tomas Koudelka, 02.2020
*/
long* loadcase::give_mstress_cn()
{
  return mstress_cn;
}



/**
  The function returns macro-strain components in the argument mstra.
  
  @param[out] mstra - the resulting vector for prescribed macro-strain components in the given load case,
                      it must be allocated to hold ncompmv components.

  @return The function returns required macro-strain components in the argumnet mstra. 

  Created by Tomas Koudelka, 02.2020
*/
void loadcase::give_mstrains(vector &mstra)
{  
  copyv(mstrain, mstra);
}



/**
  The function returns macro-stress components in the argument mstre.
  
  @param[out] mstre - the resulting vector for prescribed macro-stress components in the given load case,
                      it must be allocated to hold ncompmv components.

  @return The function returns required macro-stress components in the argumnet mstre. 

  Created by Tomas Koudelka, 02.2020
*/
void loadcase::give_mstresses(vector &mstre)
{  
  copyv(mstress, mstre);
}
