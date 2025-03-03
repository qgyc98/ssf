#include <math.h>
#include "dloadcase.h"
#include "loadcase.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "gtopology.h"
#include "matrix.h"
#include "vector.h"
#include "globmat.h"
#include "dloadn.h"
#include "dloadel.h"
#include "dloadpd.h"
#include "element.h"
#include "node.h"
#include "elemswitch.h"
#include "ipmap.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK, Tomas Koudelka,
*/
dloadcase::dloadcase (void)
{
  //  the number of subloadcases
  nslc=0L;
  //  the number of loaded nodes
  nln = 0L;
  //  the number of loaded elements
  nle = 0L;
  //  the number of prescribed displacements
  npd = 0L;
  // the maximum number of prescribed displacements
  max_npd = 0L;
  //  the number of prescribed initial displacements
  npid = 0L;
  //  the number of prescribed initial displacements
  nrpid = 0L;
  //  the number of prescribed initial displacements, tnpid = npid + nrpid, it is dimension of the array pid
  tnpid = 0L;
  
  lon = NULL;  loe = NULL;  pd = NULL;
  pid = NULL;
  //  subload cases
  slc = NULL;
  //  time functions
  gf=NULL;
}



/**
  Destructor releases allocated memory of the dloadcase object.

  Created by JK, Tomas Koudelka
*/
dloadcase::~dloadcase (void)
{
  delete [] lon;  
  delete [] loe;
  delete [] pd;
  delete [] pid;
  delete [] slc;  
  delete [] gf;
}



/**
  Function reads characteristics of time dependent load case from the opened text file.
   
  @param in - pointer to the opened XFILE
   
  Created by JK, Tomas Koudelka, 24.7.2001
*/
void dloadcase::read (XFILE *in)
{
  long i;
  
  //  type of time-dependent load
  //  tdl = 1 - time independent load
  //  tdl = 10 - seismicload
  //  tdl = 11 - responsespectrum
  //  tdl = 20 - time dependent load
  xfscanf (in,"%k%m","type_of_time_load",&dynload_kwdset,(int*)&tdl);
  
  switch (tdl){
  case timeindload:{
    //  time-independent load case multiplied by a time function

    //  number of subload cases
    xfscanf (in,"%ld",&nslc);
    
    slc = new loadcase [nslc];
    gf = new gfunct [nslc];
    
    for (i=0;i<nslc;i++){
      //  subload cases
      slc[i].read (in);
      //  time function
      gf[i].read (in);
    }
    if (check_consistency_macrostr())
      abort();
    break;
  }
  case seismicload:{
    //  seismic load
    stool.read (in);
    break;
  }
  case responsespectrum:{
    //  
    stool.read (in);
    break;
  }
  case timedepload:{
    //  time-dependent load case
    //  all nodes, elements, prescribed displacements and other
    //  quantities are equiped with their own time functions
    
    //  loaded nodes
    xfscanf (in,"%ld",&nln);
    lon = new dloadn [nln];
    for (i=0;i<nln;i++){
      lon[i].read (in);
    }
    
    //  loaded elements
    xfscanf (in,"%ld",&nle);
    loe = new dloadel [nle];
    for (i=0;i<nle;i++)
      loe[i].read (in);

    
    //  prescribed displacements
    xfscanf (in,"%ld",&npd);
    pd = new dloadpd [npd];
    for (i=0;i<npd;i++){
      pd[i].read (in);
    }
    break;
  }
    
  default:{
    print_err("unknown type of load case is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function checks for the consistency of prescribed macro-stress and macro-strain components 
  at particular subloadcases. The same components must be prescribed at all subload cases.

  @retval 0 - no consistency problem
  @retval 1 - in the case of different macro-stress/strain components prescribed at particular 
              subload cases
*/
long dloadcase::check_consistency_macrostr()
{
  long i, j;
  
  if (Mp->homog >= 3){
    for (i=1; i<nslc; i++){
      for (j=0; j<Mt->max_ncompstr; j++){
        if (slc[0].mstrastre[j] == slc[i].mstrastre[j])
          continue;
        else{
          print_err("different types of %ld-th macro-stress components have been prescribed at subload case 1 and %ld",
                    __FILE__, __LINE__, __func__, j+1, i+1);
          return 1;
        }
      }
    }
  }
  return 0;
}



/**
  Function prints characteristics of time dependent load case into the opened text file.
   
  @param out - pointer to the opened FILE
   
  TKr, 07/02/2013 according to read (XFILE *in)
*/
void dloadcase::print (FILE *out)
{
  long i;
  
  //  type of time-dependent load
  fprintf (out,"\n ##type of time-dependent load:");
  fprintf (out,"\n  %d\n",(int)tdl);
  
  switch (tdl){
  case timeindload:{
    //  time-independent load case multiplied by a time function
    //  number of subload cases
    fprintf (out,"\n #time-independent load case:");
    fprintf (out,"\n #number of subloadcases:");
    fprintf (out,"\n\n %ld",nslc);
    
    for (i=0;i<nslc;i++){
      //  subload cases
      slc[i].print (out);
      //  time function
      gf[i].print (out);
    }
    break;
  }
  case seismicload:{
    //  seismic load
    fprintf (out,"\n #seismic load:");
    stool.print (out);
    break;
  }
  case responsespectrum:{
    //  
    fprintf (out,"\n #response spectrum:");
    stool.print (out);
    break;
  }
  case timedepload:{
    //  time-dependent load case
    //  all nodes, elements, prescribed displacements and other
    //  quantities are equiped with their own time functions
    
    //  loaded nodes
    fprintf (out,"\n #time-dependent load case:");
    fprintf (out,"\n #loaded nodes:");
    fprintf (out,"\n %ld",nln);
    for (i=0;i<nln;i++){
      lon[i].print(out);
    }
    
    //  loaded elements
    fprintf (out,"\n #loaded elements:");
    fprintf (out,"\n %ld",nle);
    for (i=0;i<nle;i++)
      loe[i].print(out, 0);
    
    //  prescribed displacements
    fprintf (out,"\n #prescribed displacements:");
    fprintf (out,"\n\n %ld",npd);
    for (i=0;i<npd;i++){
      pd[i].print(out);
    }
    break;
  }
    
  default:{
    print_err("unknown type of load case is required", __FILE__, __LINE__, __func__);
  }
  }
  
}



/**
  Function assembles %vector of right hand side.
  It is created by prescribed displacements, forces, moments, eigenstrains, etc.   
  In case that the flv is not NULL, %vector of load caused by forces only is stored in flv. 
   
  @param lcid - load case id
  @param rhs - right hand side %vector
  @param flv - array of load %vector caused by forces only
  @param n - the number of DOFs in problem
  @param t - actual time

  Created by JK, 26.10.2001
  Modified by TKo 7.2008
*/
void dloadcase::assemble (long lcid,double *rhs,double *flv,long n,double t)
{
  long i,j,k,ndofn,ndofe,slcpd;
  ivector cn;
  vector r,f,af;
  matrix sm;
  double a,*aux,*auxf;
  
  //  arrays rhs and flv have to be cleaned in the caller (above) function
  
  switch (tdl){
  case timeindload:{
    //  load is defined by time independent load case, where all components
    //  are multiplied by scale factor which depends on time

    aux  = new double [n];
    auxf = new double [n];
    slcpd = 0;
    for (i=0;i<nslc;i++){
      nullv (aux,n);
      nullv (auxf,n);
      // indicator of prescribed displacements in particular subloadcases
      if (slc[i].npd)
        slcpd++;
      //  scale factor
      a=gf[i].getval(t);
      
      //  assembling of vector of contributions without prescribed displacements
      slc[i].assemblewopd (lcid, aux, auxf, a);      
      //  updating of rhs array
      //  rhs[j]+=aux[j];
      addv(rhs, aux, rhs, Ndofm);
      // assembling of load vector caused by forces only    
      if (flv)
        //  flv[j]+=auxf[j];
        addv(flv, auxf, flv, Ndofm);
    }

    // ******************************************************************
    //  contributions from prescribed displacements from all subloadcases,
    //  including Mt->nodedispl values and initial displacements on elements
    // ******************************************************************
    if (slcpd>0){
      // contributions can be assembled at once by the call of prdisplcontrib function
      // from arbitrary subloadcase because get_pd function returns sum of all prescribed displacements
      slc[0].prdisplcontrib (lcid, rhs);
    }

    // compute cumulative thermal strains from particular subloadcases

    // clean possible thermal strains at integration points because they might be recalculated
    // (they are added to Mm->tempstrains array by thermal expansion models)
    //Mm->nulltempstrains();

    for (i=0;i<nslc;i++){
      if ((slc[i].tempchang==1) || (slc[i].tempchang==2) || 
          (slc[i].tempchang==4) || (slc[i].tempchang==5))
      {
        // setup of scale factor
        if (slc[i].tempchang==2)
          a=gf[i].getval(t);
        else
          a=1.0;

        //  approximation of nodal temperatures to integration points
        intpointval(slc[i].pt,temperature,a);
    
        //  thermal material models compute temperature strains and accumulate the contributions in mechmat::tempstrains array
        Mm->cumultemprstrains();
      }
      if (slc[i].tempchang==3) // temperatures are defined outside of MEFEL and taken directly from the mechmat::nonmechq
      {
        // just one subloadcase may contain tempchang==3 flag, 
        // it is tested with the help of probdesc::temperature flag in loadcase::read

        //  thermal material models compute temperature strains and rewrites the mechmat::tempstrains array
        //Mm->temprstrains();
        // the subload case with tempchang==3 is unique and therefore we may exit the loop
        break;
      }
    }

    delete [] aux;
    delete [] auxf;
    
    break;
  }
  case seismicload:{
    //  assembles values for seismic load
    
    
    //  auxiliary array
    aux = new double [n];
    
    stool.assemble (rhs,t);

    Mmat->gmxv (rhs,aux);
    
    copyv (aux,rhs,n);
  
    // assembling of load vector caused by forces only    
    if (flv)
      copyv (aux,flv,n);
    
    delete [] aux;
    
    break;
  }
  case responsespectrum:{
    //stool.assemble (rhs,t);
    break;
  }
  case timedepload:{
    //  load is defined by fully time dependent load case
    
    // ***************************
    //  contributions from nodes
    // ***************************
    for (i=0;i<nln;i++){
      //  the number of DOFs on node
      ndofn=Mt->give_ndofn (lon[i].idn);
      if (ndofn<0){
	print_err("loaded node %ld is a hanging node",__FILE__,__LINE__,__func__,i++);
      }
      
      for (j=0;j<ndofn;j++){
	k=Mt->give_dof (lon[i].idn,j);
	if (k<=0)  continue;
	if (k>0)  rhs[k-1]+=lon[i].getval(t,j);
      }
    }//  end of the loop over the number of loaded nodes
    
    
    //  contributions from elements
    for (i=0;i<nle;i++)
    {
      ndofe=Mt->give_ndofe(loe[i].eid);
      reallocv(RSTCKIVEC(ndofe, cn));
      Mt->give_code_numbers (loe[i].eid, cn.a);
      /// computes nodal values from element load
      loe[i].compute_load(t);
      for (j=0;j<ndofe;j++){
        k=cn[j]-1;
        if (k<0)  continue;
        else  rhs[k]+=loe[i].nf[j];
      }
    }
    
    // assembling of load vector caused by forces only    
    if (flv)
      copyv(rhs,flv,n);
    
    // **********************************************
    //  contributions from prescribed displacements
    // **********************************************
    if (npd>0){
      /*
      // commented out due to particle elements which lacks the stiffness matrix, 
      // it results in additional iteration step in the case of linear elastic material model

      //  the number of elements
      long ne = Mt->ne;
      //  loop over the number of elements
      for (i=0;i<ne;i++){
	if (Mt->elements[i].prescdispl==1){
	  //  only elements with prescribed displacements are taken into account
	  
	  //  the number of DOFs on element
	  //  this number is equal to the number of DOFs without hanging nodes
	  ndofe=Mt->give_ndofe (i);
	  reallocm (ndofe,ndofe,sm);
	  //  stiffness matrix of an element
	  stiffmat (lcid,i,sm);
	  
	  reallocv (ndofe,r);
	  reallocv (ndofe,f);

	  //  prescribed displacements on element
	  elprdispl (lcid,i,r.a);
	  //  Kr=f
	  mxv (sm,r,f);
	  cmulv (-1.0,f);
	  
	  //  the number of DOFs on element
	  //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
	  long ndofemn = Gtm->give_ndofe (i);
	  reallocv (ndofemn,cn);
	  
	  //  code numbers on element
	  Mt->give_code_numbers (i,cn.a);
	  
	  if (ndofe != ndofemn){
	    //  the element contains hanging node
	    //  the vector f has to be transformed
	    reallocv (ndofemn,af);
	    mtxv (*Mt->elements[i].tmat,f,af);
	    locglob (rhs,af.a,cn.a,ndofemn);
	  }else{
	    //  the element does not contain hanging node
	    locglob (rhs,f.a,cn.a,ndofe);
	  }
	  
	}
      }
      */
    }
    
    break;
  }
    
  default:{
    print_err("unknown type of load case is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
  Function assembles %vector of right hand side with respect of left hand side
  node forces cause by the attained displacements.
  It is used in epsolver.cpp only.
   
  @param rhs - right hand side %vector
  @param lhs - left hand side
  @param flv   - array of load %vector caused by forces only
   
  Created by Tomas Koudelka, 26.10.2001
  Modified by Tomas Koudelka, 7.2008
*/
void dloadcase::assemble (double *rhs, double *lhs, double *flv)
{
  long i,j,k,ndofn;

  //  array rhs has to be cleaned in the caller (above) function
  
  // ***************************
  //  contributions from nodes
  // ***************************
  for (i=0;i<nln;i++){
    ndofn=Mt->give_ndofn (lon[i].idn);
    for (j=0;j<ndofn;j++){
      k=Mt->give_dof (lon[i].idn,j);
      if (k<=0)  continue;
      if (k>0)  rhs[k-1]+=lon[i].getval(lhs[k-1], j);
    }
  }
  // assembling of load vector caused by forces only    
  if (flv)
    copyv(rhs,flv,Ndofm);
}



/**
  Function evaluates reactions.
   
  @param lcid - load case id

  @return The function does not return anything.

  Created by JK,   
*/
void dloadcase::compute_reactions (long lcid)
{
  long i, j, k, ii, nne, ndofn, ne, ndofe, ndofemn;
  long eid, nid;
  vector f, anf;
  ivector nod;

  if (tdl == timeindload){
    answertype ifc = yes;
    
    for (i=0; i<nslc; i++, ifc=no)
      slc[i].compute_reactions(lcid, gf[i].getval(Mp->time), ifc);
    return;
  }

   
  // *************************************************************
  // Compute contributions of element internal forces to reactions
  // *************************************************************

  ne=Mt->ne;  //  the number of elements
  for (i=0; i<ne; i++){  //  loop over the number of elements
    if ((Mt->elements[i].react == 1) && (Gtm->leso[i] == 1)){
      //  only elements connected to a node with reaction are taken into account
      
      //  the number of DOFs on element
      //  this number is generated by elements without hanging nodes
      ndofe = Mt->give_ndofe(i);
      reallocv(RSTCKVEC(ndofe, f));
      //  nodal forces are computed
      elem_internal_forces(i, lcid, f);
      
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe(i);
      
      //  the number of nodes on element
      //  if there are hanging nodes, the variable will be rewritten
      nne = Mt->give_nne(i);
      
      if ((ndofe == ndofemn) || give_num_mstress_comp()){
        // there are no hanging nodes or there is homogenization with macro-stress
	reallocv(nne, nod);
	Mt->give_elemnodes(i, nod);
      }
      else{
	//  the element contains hanging node
	//  the vector of internal forces has to be transformed
	reallocv(RSTCKVEC(ndofemn,anf));
	//  transformation of the vector f
	mtxv(*Mt->elements[i].tmat, f, anf);
	reallocv(RSTCKVEC(ndofemn, f));
	copyv(anf, f);
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
	//  the number of DOFs in node
	ndofn = Mt->give_ndofn(nod(j));
	if (Mt->nodes[nod(j)].react == 1){
	  //  only nodes with constraint/reaction are taken into account
	  
	  //  loop over the number of DOFs in the actual node
	  for (k=0; k<ndofn; k++){
	    Mt->nodes[nod(j)].r[k]+=f[ii];
	    ii++;
	  }//  end of the loop over the number of DOFs in the actual node
	}
	else  ii+=ndofn;
      }//  end of the loop over the number of nodes on element
      
    }
  }//  end of the loop over the number of elements

  
  // *****************************************************************************
  //  contributions from nodal load applied at nodes with prescribed displacements
  // *****************************************************************************
  // prepared for moving load
  for (i=0; i<nln; i++){
    nid   = lon[i].idn; // node id
    ndofn = Mt->give_ndofn(nid);
    if (ndofn<0)
      print_err("loaded node %ld is a hanging node",__FILE__,__LINE__,__func__,nid+1);

    for (j=0; j<ndofn; j++){      
      k = Mt->give_dof(nid, j);
      if (k<=0) Mt->nodes[nid].r[j] -= lon[i].getval(Mp->time, j); // nodal load is being applied at the direction of prescribed displacement
      if (k>0)  continue;
    }
  }

  
  // **************************************************
  // Compute contributions of element load to reactions
  // **************************************************

  for (i=0; i<nle; i++){    //  loop over the loaded elements
    eid=loe[i].eid;    //  element id
    if ((Mt->elements[eid].react == 1) && (Gtm->leso[eid] == 1)){ 
      //  only switched on elements with reactions are taken into account

      ndofe = Mt->give_ndofe(eid);  //  the number of DOFs on element
      //  the number of DOFs on element
      //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
      ndofemn = Gtm->give_ndofe(eid);
      if ((ndofe == ndofemn) || give_num_mstress_comp()){
	  //  there are no hanging nodes or there is homogenization
	  nne = Mt->give_nne(eid);  //  the number of nodes on element
	  reallocv(RSTCKIVEC(nne, nod));
	  Mt->give_elemnodes(eid, nod);
          reallocv(RSTCKVEC(ndofe, f));
          copyv(loe[i].nf, f);
      }
      else{
        //  the element contains hanging node
        //  the load vector has to be transformed
        reallocv (RSTCKVEC(ndofemn, f));
        //  transformation of the vector f
        mtxv (*Mt->elements[eid].tmat, loe[i].nf, f);
	  
        //  the number of nodes on element
        //  it is equal to the number of master nodes
        nne=Gtm->gelements[eid].nmne;
        reallocv(RSTCKIVEC(nne, nod));
        Gtm->gelements[eid].give_master_nodes(nod);
      }

      //  at this moment, the arrays nod and anf contain the appropriate number of components

      ii=0;
      //  loop over the number of nodes on element
      //  if hanging nodes are present, the loop takes them into account
      for (j=0; j<nne; j++) {
	ndofn=Mt->give_ndofn(nod(j)); //  the number of DOFs in node
	if (Mt->nodes[nod(j)].react == 1){
	  //  only nodes where reaction is are taken into account
	  for (k=0;k<ndofn;k++){
	    Mt->nodes[nod(j)].r[k] -= f[ii];
	    ii++;
	  }
	}
	else  ii+=ndofn;
      }
    }
  }
}



/**
  The function reallocates memory for array of prescribed initial displacements pid.
  The actual dimension of array must be stored at tnpid data member and the array is set to zero after 
  reallocation.
  The array pid is used in course of initial displacement calculation of new activated part 
  of structure at simulations of gradual construction process. 

  @return The function does not return anthing but it changes the content of pid.

  Created by Tomas Koudelka, 08.2016
*/
void dloadcase::realloc_pid_array()
{
  // setup of array pid with prescribed intial displacements of the new active part of structure
  if (pid)
    delete [] pid;
  pid = new double[tnpid];
  memset(pid, 0, sizeof(*pid)*tnpid);
}



/**
  The function computes the maximum number of different prescribed nodal values of displacements at 
  the given time dependnet load case, i.e. it checks either all subloadcases for maximum of npd 
  or it returns dloadcase.npd.

  @return The function does not returns anything but it changes value of dloadcase::max_npd.

  Created by Tomas Koudelka, 08.2016
*/
void dloadcase::compute_max_npd()
{
  long i;

  max_npd = 0L;

  switch (tdl)
  {
    case timeindload:
      for(i=0L; i<nslc; i++)
      {
        if (slc[i].npd > max_npd)
          max_npd = slc[i].npd;
      }
      break;
    case timedepload:
      max_npd = npd;
      break;
    default:
      print_err("unknown type of time dependent load case is required", __FILE__, __LINE__, __func__);
  }
}


/**
  Function returns prescribed displacement for the given time and dof.
   
  @param time - actual time
  @param dof - DOF id

  @return The function returns prescribed displacement.
   
  Created by Tomas Koudelka, 14.11.2007
  Modified by Tomas Koudelka, 07.2016,
  Corrected by TKr, 20/04/2023
*/
double dloadcase::get_pd(double time, long dof)
{
  double ret = 0.0;
  double a;
  long i;
  
  switch (tdl){
  case timeindload:{
    for (i=0; i<nslc; i++)
    {
      //  scale factor
      a=gf[i].getval(time);
      //  assembling of contributions
      if (slc[i].npd && (slc[i].npd >= -dof))
        ret += slc[i].pd[0-dof-1]*a;
    }   
    if (tnpid && (-dof > max_npd)) // initial prescribed displacements at gradual construction problem have been defined
      ret = pid[-(dof+max_npd+1)]; // dof number of prescribed initial displacement is negative and it is less than all dof numbers 
                                   // for prescribed displacements defined in particular subloadcases
    break;
  }
  case timedepload:{
    if (npd)
      //ret = pd[0-dof-1].getval(Mp->time);
      ret = pd[0-dof-1].getval(time); // corrcted for values for the given time
    if (tnpid && (-dof > max_npd)) // initial prescribed displacements at gradual construction problem have been defined
      ret = pid[-(dof+max_npd+1)]; // dof number of prescribed initial displacement is negative and it is less than all dof numbers 
                                   // for prescribed displacements defined in particular subloadcases
    break;
  }
  default:{
    print_err("this type of time dependent load case is not supported", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  return ret;
}



/**
  Function returns cn-th component of macro-stress.

  @param time - actual time
  @param cn - component number

  @return The function returns the value of macro-stress component.
 
  Created by Tomas Koudelka, 16.2.2015
*/
double dloadcase::get_macrostre(double time, long cn)
{
  double ret = 0.0;
  double a;
  long i;

  switch (tdl)
  {
    case timeindload:
      for (i=0; i<nslc; i++)
      {
        //  scale factor
        a=gf[i].getval(time);
        //  assembling of contributions
        if (slc[i].mstrastre[cn] == stress)
          ret += slc[i].mstress[cn]*a;
      }
    default:{
      print_err("this type of load case is not supported for the macro-stress", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  return ret;
}



/**
  Function returns macro-stress %vector.

  @param[in] time - actual time
  @param[out] msig - required macro-stress %vector in Voigt notation.

  @return The function returns the value of macro-stress %vector in the argument msig.
 
  Created by Tomas Koudelka, 09.2023
*/
void dloadcase::get_macrostre(double time, vector &msig)
{
  double a;
  long i, j;

  switch (tdl)
  {
    case timeindload:
      nullv(msig);
      for (i=0; i<nslc; i++)
      {
        //  scale factor
        a=gf[i].getval(time);
        for (j=0; j<msig.n; j++){          
          //  assembling of contributions        
          if (slc[i].mstrastre[j] == stress)
            msig(j) += slc[i].mstress[j]*a;
        }
      }
    default:{
      print_err("this type of load case is not supported for the macro-stress", __FILE__, __LINE__, __func__);
      abort();
    }
  }
}



/**
  Function returns cn-th component of macro-strain.

  @param time - actual time
  @param cn - component number

  @return The function returns the value of macro-strain component.
 
  Created by Tomas Koudelka, 16.2.2015
*/
double dloadcase::get_macrostra(double time, long cn)
{
  double ret = 0.0;
  double a;
  long i;

  switch (tdl)
  {
    case timeindload:
      for (i=0; i<nslc; i++){
        //  scale factor
        a=gf[i].getval(time);
        //  assembling of contributions
        if (slc[i].mstrastre[cn] == strain)
          ret += slc[i].mstrain[cn]*a;
      }
    default:{
      print_err("this type of load case is not supported for the macro-stress", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  return ret;
}



/**
  Function returns cn-th component of macro-strain.

  @param time - actual time
  @param cn - component number

  @return The function returns the value of macro-strain component.
 
  Created by Tomas Koudelka, 16.2.2015
*/
void dloadcase::get_macrostra(double time, vector &meps)
{
  double a;
  long i, j;

  switch (tdl)
  {
    case timeindload:
      nullv(meps);
      for (i=0; i<nslc; i++){
        //  scale factor
        a=gf[i].getval(time);
        for(j=0; j<meps.n; j++){
          //  assembling of contributions
          if (slc[i].mstrastre[j] == strain)
            meps(j) += slc[i].mstrain[j]*a;
        }
      }
    default:{
      print_err("this type of load case is not supported for the macro-stress", __FILE__, __LINE__, __func__);
      abort();
    }
  }
}



/**
  Function computes contributions from temperature changes
  to the nodal forces.

  @param lcid - load case id
  @param rhs - pointer to the right hand side
  @param n - number of components of the array rhs
  @param t - actual time
   
  Created by JK, 21.11.2007
*/
void dloadcase::tempercontrib (long lcid,double *rhs,long n,double t)
{
  long i;
  double a;
  double *aux;

  switch (tdl){
  case timeindload:{
    aux = new double [n];
    for (i=0;i<nslc;i++){
      nullv (aux,n);
      //  scale factor
      a=gf[i].getval(t);
      
      //  assembling of vector of contributions
      slc[i].tempercontrib (lcid,aux,a);
      addv (rhs,aux,n);
    }
    delete [] aux;
    break;
  }
  default:
    print_err("unknown type of load case is required", __FILE__, __LINE__, __func__);
  }

}



/**
  The function computes temperature strains at auxiliary integration points
  for the given load case lcid. Temperature strains are computed cumulatively for all subloadcases
  of the given load case. It is intended for the use in coupled problems esspecially.
  The load case must be type of time dependent load case, dloadcase instance must be defined.

  @param lcid[in] - load case id
  @param n[in]    - the number of required auxiliary points in the mapping array ipm
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of strains is performed, the strains are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it updates content of mechmat::aip_tempstains array.

  Created by Tomas Koudelka, 18.6.2018
*/
void dloadcase::aip_temperstrains(long lcid, long n, ipmap *ipm)
{
  long i;
  double a, t = Mp->time;

  Mm->nullaiptempstrains();

  switch (tdl)
  {
    case timeindload:
    {
      for (i=0;i<nslc;i++)
      {
        if (slc[i].tempchang)
        {
          // setup of scale factor
          if (slc[i].tempchang==2)
            a=gf[i].getval(t);
          else
            a=1.0;
          slc[i].aip_cumultemperstrains_contrib(lcid, a, n, ipm);
        }
      }
    }
    case seismicload:
    case responsespectrum:
    case timedepload:
      break;
    default:
      print_err("unknown type of load case is required",__FILE__,__LINE__,__func__);
      abort();
  }
}



/**
  The function returns number of prescribed macro-stress components at the given load case.
  It is intended for the homogenization problems (Mp->homog = 3,5,7,9)

  @return The number of prescribed macro-stress components.
*/
long dloadcase::give_num_mstress_comp()
{
  long ret=0;
  
  if (slc){
    // all subload cases must have the same number of macro-stress components
    // (it is checked in dloadcase::read function) thus the number can be taken
    // from the arbitrary subload case, e.g. the first one
    ret = slc[0].give_num_mstress_comp();
  }
  return ret;
}
 


/**
  The function returns number of prescribed macro-strain components at the given load case.
  It is intended for the homogenization problems (Mp->homog = 4,6,8,9)

  @return The number of prescribed macro-strain components.
*/
long dloadcase::give_num_mstrain_comp()
{
  long ret=0;
  
  if (slc){
    // all subload cases must have the same number of macro-strain components
    // (it is checked in dloadcase::read function) thus the number can be taken
    // from the arbitrary subload case, e.g. the first one
    ret = slc[0].give_num_mstrain_comp();
  }
  return ret;
}
 


/**
  The function returns pointer to the array of types of prescribed macro-value components 
  at the given load case. It is intended for the homogenization problems (Mp->homog > 2)

  @return The ponter to the array of indivators of macro-stress components.
*/
strastre* dloadcase::give_mstrastre()
{
  strastre *ret = NULL;
  
  if (slc){
    // all subload cases must have the same array of macro-value compoennt types
    // (it is checked in dloadcase::read function) thus the pointer can be taken
    // from the arbitrary subload case, e.g. the first one
    ret = slc[0].give_mstrastre();
  }
  return ret;
} 



/**
  The function returns pointer to the array of code (DOF) numbers of macro-stress components 
  at the given load case. It is intended for the homogenization problems (Mp->homog > 2)

  @return The pointer to the array of macro-stress component code (DOF) numbers,
          dimension of the array is the maximum total number of stress/strain components of used elements.
          Particular components of the array contains either positive nonzero value which represents the code (DOF) number
          or zero value for components where no macro-stress value has been prescribed.
*/
long* dloadcase::give_mstress_cn()
{
  long *ret = NULL;
  
  if (slc){
    // all subload cases must have the same array of macro-value compoennt types
    // (it is checked in dloadcase::read function) thus the pointer can be taken
    // from the arbitrary subload case, e.g. the first one
    ret = slc[0].give_mstress_cn();
  }
  return ret;
} 



/**
  The function returns %vector of actual prescribed macro-strain components 
  at the given load case. It is intended for the homogenization problems (Mp->homog > 2)

  @param[in] time - actual time
  @param[out] mstra - %vector of the resulting macro-strain components at the given time

  @return The function returns value of prescribed macro-strain components at the given time,
*/
void dloadcase::give_mstrains(double time, vector &mstra)
{
  long i, j, ncompmv;
  double a;  

  nullv(mstra);

  if (slc){
    // all subload cases must have the same total number of macro-value components,
    // i.e. this number can be taken from arbitrary subload case
    ncompmv = slc[0].ncompmv;
    for (i=0; i<nslc; i++){
      //  scale factor
      a = gf[i].getval(time);
      for (j=0; j<ncompmv; j++)
        mstra(j) += a*slc[i].mstrain[j];
    }
  }
} 



/**
  The function returns %vector of actual prescribed macro-stress components 
  at the given load case. It is intended for the homogenization problems (Mp->homog > 2)

  @param[in] time - actual time
  @param[out] mstre - %vector of the resulting macro-stress components at the given time

  @return The function returns value of prescribed macro-strain components at the given time,
*/
void dloadcase::give_mstresses(double time, vector &mstre)
{
  long i, j, ncompmv;
  double a;  

  nullv(mstre);

  if (slc){
    // all subload cases must have the same total number of macro-value components,
    // i.e. this number can be taken from arbitrary subload case
    ncompmv = slc[0].ncompmv;
    for (i=0; i<nslc; i++){
      //  scale factor
      a = gf[i].getval(time);
      for (j=0; j<ncompmv; j++)
        mstre(j) += a*slc[i].mstress[j];
    }
  }
} 
