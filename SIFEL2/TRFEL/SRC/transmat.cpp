#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "transmat.h"
#include "globalt.h"
#include "globmatt.h"
#include "intpointst.h"
#include "elementt.h"
#include "nodet.h"
#include "onemedium.h"
#include "twomedia.h"
#include "threemedia.h"
#include "fourmedia.h"
#include "multiphase.h"
#include "gmultiphase.h"
#include "elemswitcht.h"
#include "ipmap.h"

transmat::transmat (void)
{
  //  total number of integration points
  tnip=0;
  // total number of auxiliary integration points
  tnaip=0;
  
  //  integration points
  ip=NULL;
  // auxiliary integration points (for transfer of quantities in coupled problems)
  aip = NULL;
  //  element - integration point map
  elip = NULL;
  //  element - auxiliary integration point map (for transfer of quantities in coupled problems)
  elaip = NULL;

  //  isotropic material
  itrm=NULL;
  //  nonlinear isotropic material
  nlitrm=NULL;
  //  artificial material obtained from homogenization
  hommat=NULL;
  //  isotropic material with optional influence of damage
  damitrm=NULL;
  //  time dependent isotropic diffusion 
  tditrm=NULL;
  
  //  material for contact elements
  ifacemat=NULL;
  
  cernym=NULL;

  sejtkrm=NULL;

  kun=NULL;
  kun2=NULL;
  grunw=NULL;
  dvries=NULL;
  bazped=NULL;
  ped=NULL;
  carb1=NULL;
  sdmat=NULL;
  mill = NULL;
  moisth = NULL;

  salt1 = NULL;
  salt2 = NULL;
  salt3 = NULL;
  salt4 = NULL;

  consol_awf1 = NULL;
  consol_wf1 = NULL;
  consol_wf2 = NULL;
  consol_awf2 = NULL;
  consol_hawf3 = NULL;
  consol_hwf2 = NULL;

  concrete = NULL;
  baroghel = NULL;
  C60baroghel = NULL;
  C30baroghel = NULL;
  o30bazant = NULL;
  C60bazant = NULL;
  C30bazant = NULL;
  tench = NULL;

  soil1 = NULL;
  richar = NULL;
  cemhydr = NULL;
  
  lcmat=NULL;

  mattype = NULL;
  numtype = NULL;

  //  initial temperature at integration points
  initval = NULL;
  //  nontransport quantities
  nontransq = NULL;
  //  array of values of non-transport quantities for auxliary integration points
  aip_nontransq = NULL;
  nntq = 0;
  for (long i=0; i<tnkntq; i++)
    ntqid[i] = -1;
  ntqo = NULL;
}

transmat::~transmat (void)
{
   long i;

  delete [] ip;
  delete [] elip;
  delete [] aip;
  delete [] elaip;

  delete [] itrm;  
  delete [] nlitrm;  
  delete [] hommat;
  delete [] damitrm;
  delete [] tditrm;
  delete [] ifacemat;

  delete [] cernym;

  delete [] sejtkrm;  

  delete [] kun;
  delete [] kun2;
  delete [] grunw;
  delete [] dvries;
  delete [] bazped;
  delete [] ped;
  delete [] carb1;
  delete [] sdmat;
  delete [] mill;
  delete [] moisth;
  
  delete [] salt1;
  delete [] salt2;
  delete [] salt3;
  delete [] salt4;

  delete [] consol_awf1;
  delete [] consol_wf1;
  delete [] consol_wf2;
  delete [] consol_awf2;
  delete [] consol_hawf3;
  delete [] consol_hwf2;

  delete [] concrete;
  delete [] baroghel;
  delete [] C60baroghel;
  delete [] C30baroghel;
  delete [] o30bazant;
  delete [] C60bazant;
  delete [] C30bazant;
  delete [] tench;
  
  delete [] richar;
  
  delete [] soil1;

  delete [] lcmat;

  delete [] mattype;
  delete [] numtype;
  
  delete [] initval;

  if (nontransq != NULL){
    for (i=0;i<nntq;i++){
      delete [] nontransq[i];
    }
    delete [] nontransq;
  }
  if (aip_nontransq != NULL){
    for (i=0;i<nntq;i++){
      delete [] aip_nontransq[i];
    }
    delete [] aip_nontransq;
  }
  delete [] ntqo;
}


/**
   function returns total number of integration points
   function computes number of the first integration
   point on elements
   
   8.7.2001
*/
long transmat::intpnum (void)
{
  long i,j,k,n;
  
  n=0;
  for (i=0;i<Tt->ne;i++){
    Tt->elements[i].ipp = new long* [Tp->ntm];
    for (j=0;j<Tp->ntm;j++){
      Tt->elements[i].ipp[j] = new long [Tp->ntm];
      for (k=0;k<Tp->ntm;k++){
	Tt->elements[i].ipp[j][k]=n;
	n+=Tt->give_nip (i,j,k);
      }
    }
  } 
  return n;
}

/**
   function reads material types and material id
   for each integration point
   
   @param in - input stream
   
   8.7.2001
*/
void transmat::readip (FILE *in)
{
  long i;
  
  ip = new intpointst [tnip];
  
  //  allocation of integration points
  for (i=0;i<tnip;i++)
    ip[i].alloc();
  
  //  reading of input data
  for (i=0;i<tnip;i++){
    ip[i].read (in);
  }
}

/**
   function allocates integration points
*/
void transmat::intpointalloc ()
{
  ip = new intpointst [tnip];
}



/**
  The function allocates array of auxiliary integration points (aip)
  used for the transfer of quantities among different meshes.
  It is used in the coupled problems especially (it is called in metrinit). The array 
  is allocated according to the mapping of integration points which is 
  given in ipm. 
  The function allocates aip only for those ipmap components which
  do not contain direct mapping to the regular integration point array 
  ip, i.e. ipm[i].ipp < 0.

  @param n - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]

  @return The function does not return anything but it changes the class
          data members of aip and elaip arrays and tnaip value.

  Created by Tomas Koudelka, 28.11.2017
*/
void transmat::alloc_aux_intp(long n, ipmap *ipm)
{
  long i;

  // compute the number of aip
  for (i=0, tnaip=0; i<n; i++)
  {
    if (ipm[i].ipp < 0)
    {
      ipm[i].app = tnaip;
      tnaip++;
    }
  }

  // release memory of  
  if (aip)
  {
    delete [] aip;
    delete [] elaip;
    for (i=0; i<nntq; i++)
      delete [] aip_nontransq[i];
    delete [] aip_nontransq;
  }

  if (tnaip)
  {
    aip = new intpointst[tnaip];
    // set material models, number of flux/gradient components
    // allocate flux/grad/other/eqother and non-transport quantity arrays
    aipinit(n, ipm);
  }
  else
  {
    // tnaip=0 =>zero number of auxiliary integration points
    aip = NULL;
    elaip = NULL;
    aip_nontransq = NULL;
  }
}



/**
  The function allocates array of non-transport quantities.
  Number of rows is given by parameter n and number of columns is 
  equal to the total number of integration points.

  @param n - number of non-transport quantities

  @return The function does not return anything.

  Created by Tomas Koudelka, 7.10.2013
*/
void transmat::alloc_nontransq(long n)
{
  long i, j;

  nntq = n;
  nontransq = new double*[nntq];
  
  for (i=0; i<nntq; i++)
  {
    nontransq[i] = new double[tnip];
    // set zero initial values of given quantity
    memset(nontransq[i], 0, sizeof(*nontransq[i])*tnip); 
  }
  
  if (ntqo == NULL)
    {
      ntqo = new nontransquant[nntq];
      j=0;  
      for (i=0; i<tnkntq; i++)
	{
	  if (ntqid[i] > -1)
	    {
	      ntqo[j] = nontransquant(i+1);
	      j++;
	    }
	}
      if (j != nntq)
	{
	  print_err("indices of used non-transport qunatities must be initialized before allocation", __FILE__, __LINE__, __func__);
	  abort();
	}
    }
  return;
}



/**
  Function initializes material models, fluxes/gradient/other/eqother on auxiliary integration points.
  Material models are copied from the corresponding elements.
   
  The function initializes material models only for those auxiliary integration points
  in which ipmap components do not contain direct mapping to the regular integration point array 
  ip, i.e. ipm[i].app < 0.

  @param n - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]

  @return The function does not return anything.

  Created by TKo, 21.11.2017
*/
void transmat::aipinit (long n, ipmap *ipm)
{
  long i, eid, app;
  //  long nm = Tp->ntm*Tp->ntm;
  
  //  element - auxiliary integration point mapping
  elaip = new long[tnaip];
  
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)
      continue;
    eid = ipm[i].eid;
    if (elaip != NULL)
       elaip[app] = eid;
    //    aip[app].tm  = new mattypet[nm];
    //    aip[app].idm = new long[nm];

    // copy the number of flux/gradient components
    aip[app].ncompgrad = Tt->give_ncomp(eid);
    // copy material model definitions
    //    for (j=0; j<nm; j++)
    //    {
    //      aip[app].tm[j]  = Tt->elements[eid].tm[j];
    //      aip[app].idm[j] = Tt->elements[eid].idm[j];
    //    } 
    aip[app].tm  = Tt->elements[eid].tm[0];
    aip[app].idm = Tt->elements[eid].idm[0];
    // allocate arrays for gradients, fluxes and state variables on the given auxiliary integration point
    aip[app].alloc();
  }

  // allocate array of values for non-transport quantities at auxiliary integration points
  if (nntq) 
  {
    // some non-transport quantities are required in the material models
    // their number and ordering is the same as for regular integration point
    aip_nontransq = new double*[nntq];
    for (i=0; i<nntq; i++)
    {
      aip_nontransq[i] = new double[tnaip];
      // set zero initial values of given quantity
      memset(aip_nontransq[i], 0, sizeof(*aip_nontransq[i])*tnaip);    
    }  
  }
}



/**
   function reads advection velocities and stores them to
   
   @param in - input data stream 
   
   JK, 19. 5. 2016
*/
void transmat::advection_velocities (XFILE *in)
{
  long i,j,k,nn,gdim,n;
  double **av;
  
  //  the number of nodes
  nn = Tt->nn;
  //  geometrical dimension
  gdim = Tp->gdim;
  
  //  the number of nontransport quantites is determined
  n = search_reqntq (ntqo);
  if (n){
    for (i=0; i<n; i++){
      //  
      ntqid[ntqo[i]-1] = i;
    }
    //  allocation of the array nontransq
    alloc_nontransq (n);
  }
  
  av = new double* [gdim];
  for (i=0;i<gdim;i++){
    av[i] = new double [nn];
  }
  
  for (i=0;i<nn;i++){
    xfscanf (in,"%ld",&k);
    k--;
    for (j=0;j<gdim;j++){
      xfscanf (in,"%lf",&av[j][k]);
    }
  }
  
  switch (gdim){
  case 1:{
    intpointvalt (av[0],advect_vel_x,1.0);
    break;
  }
  case 2:{
    intpointvalt (av[0],advect_vel_x,1.0);
    intpointvalt (av[1],advect_vel_y,1.0);
    break;
  }
  case 3:{
    intpointvalt (av[0],advect_vel_x,1.0);
    intpointvalt (av[1],advect_vel_y,1.0);
    intpointvalt (av[2],advect_vel_z,1.0);
    break;
  }
  default:{
    print_err("unknown geometrical dimension %ld is required",__FILE__,__LINE__,__func__,gdim);
    abort();
  }
  }
  
  for (i=0;i<gdim;i++){
    delete av[i];
  }
  delete [] av;
}


/**
   function initializes basic data in integration points
*/
void transmat::intpointinit ()
{
  long i,j,k,ii,jj,nb,nip,ipp,ncomp;

  elip = new long[tnip];
 
  for (i=0;i<Tt->ne;i++){
    nb=Tp->ntm;  k=0;
    ncomp=Tt->give_ncomp (i);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
	nip=Tt->give_nip(i,ii,jj);
	ipp=Tt->elements[i].ipp[ii][jj];
	for (j=0;j<nip;j++){
          if (elip != NULL)
            elip[ipp] = i;
          ip[ipp].ncompgrad=ncomp;
	  ip[ipp].tm   = Tt->elements[i].tm[k];
	  ip[ipp].idm  = Tt->elements[i].idm[k];
	  ipp++;
	}
	k++;
      }
    }
  }
  
  for (i=0;i<tnip;i++)
    ip[i].alloc();
}



/**
  The function assembles ordering of nodal uknowns and dof names
  with respect to material model used.

  @param dofname   - array of uknown names for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname
  @param var_dofid - array with nodal dof indices for all implemented types of nodal unknowns (output)
  @param tnkv      - total number of implemented nodal unknown = length of array var_dofid

  @return The function returns dof ordering and dof names in the parameters dofname and var_dofid.

  Created by Tomas Koudelka, 6.12.2013
*/
void transmat::assemble_dof_nameord(namevart *dofname, long ntm, long *var_dofid, long tnkv)
{
  long i, j, id, reord;
  
  reord = 0;
  for(i=0; i<tnip; i++)
  {
    memset(dofname, 0, sizeof(*dofname)*ntm);
    give_dof_names(i, dofname, ntm);

    for(j=0; j<ntm; j++)
    {
      if (dofname[j])
        id = var_dofid[int(dofname[j])-1];
      else // no dof name has been assigned in the material to the given dof, go to the next dof, i.e. medium
        continue;

      if (id < 0) // dof id has not been assigned yet
        var_dofid[int(dofname[j])-1] = j;
      else  // dof id has already been assigned
      {
        if (id != j) // assigned dof id is different from the actual one
          reord = 1; // dof ids must be reordered
      }
    }
  }

  if (reord)
  { // some of dof names has got assigned the same dof id in the material model
    // it can be caused by usage of different material types on an element
    for(i=0, j=0; i<tnkv; i++)
    {
      if (var_dofid[i] >= 0) // if some dof id has been assigned for the i-th type of primary unknown
      {
        var_dofid[i] = j;  // generate new dof id
        j++;
      }
    }
    if (j != ntm)
    {
      print_err("number of nodal unknowns detected in material models(=%ld) is different\n"
                "than the one specified in the probdesct(=%ld)", __FILE__, __LINE__, __func__, j, ntm);
      abort();
    }
    // generate new dof names according to new order od dofs
    for(i=0; i<tnkv; i++)
    {
      if (var_dofid[i] >= 0)
        dofname[var_dofid[i]] = namevart(i+1);
    }
  }
}



/**
  The function  the order of dof names with respect to material type in the given .

  @param ipp - integration point id
  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (fro names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  @return The function returns ordered dof names in the parameter dofname.

  Created by Tomas Koudelka, 6.12.2013
*/
void transmat::give_dof_names(long ipp, namevart *dofname, long ntm)
{
  long i = ip[ipp].idm;

  switch (ip[ipp].tm)
    {
    case consolawf1:
      consol_awf1[i].give_dof_names(dofname, ntm);
      break;

    case consolwf1:
      consol_wf1[i].give_dof_names(dofname, ntm);
      break;      

    case consolwf2:
      consol_wf2[i].give_dof_names(dofname, ntm);
      break;

    case consolawf2:
      consol_awf2[i].give_dof_names(dofname, ntm);
      break;

    case consolhawf3:
      consol_hawf3[i].give_dof_names(dofname, ntm);
      break;

    case consolhwf2:
      consol_hwf2[i].give_dof_names(dofname, ntm);
      break;
      
    case isotransmat:
      itrm[i].give_dof_names(dofname, ntm);
      break;
    case damisotransmat:
      damitrm[i].give_dof_names(dofname, ntm);
      break;
    case nlisotransmat:
      nlitrm[i].give_dof_names(dofname, ntm);
      break;
    case tdisotransmat:
      tditrm[i].give_dof_names(dofname, ntm);
      break;
      
      /*
	case nlisotransmat:
	//nlitrm[i].give_dof_names(dofname, ntm);
	break;
      */
    case homomat:{
      hommat[i].give_dof_names(dofname, ntm);
      break;
    }
    case interfacem:{
      break;
    }

      /*
	case discontisotrmat:
	//ditrm[i].give_dof_names(dofname, ntm);
	break;
	
	case cernyconcrete:
	//cernym[i].give_dof_names(dofname, ntm);
	break;
	
	case bazantpedersen:
	//bazped[i].give_dof_names(dofname, ntm);
	break;
	
	case pedersen:
	//ped[i].give_dof_names(dofname, ntm);
	break;
      */      
    case kunzel:
      kun[i].give_dof_names(dofname, ntm);
      break;
    case moistheat:{
      moisth[i].give_dof_names(dofname, ntm);
      break;
    }
      
    case baroghelB:{
      baroghel[i].give_dof_names(dofname, ntm);
      break; 
    }
      
    case C30bazantB:
      C30bazant[i].give_dof_names(dofname, ntm);
      break;
      
    case C60bazantB:{
      C60bazant[i].give_dof_names(dofname, ntm);
      break;
    }       

      /*      
	      case kunzel2:
	      //kun2[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case grunewald:
	      //grunw[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case simplediscmat:
	      //sdmat[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case devries:
	      //dvries[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case milly:
	      //mill[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case concreteB:
	      //concrete[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case baroghelB:
	      //baroghel[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case C60baroghelB:
	      //C60baroghel[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case C30baroghelB:
	      //C30baroghel[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case o30bazantB:
	      //o30bazant[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case C60bazantB:
	      //C60bazant[i].give_dof_names(dofname, ntm);
	      break;
	      
	      case C30bazantB:
	      //C30bazant[i].give_dof_names(dofname, ntm);
	      break;
	      
      */  
    case glasgow:
      tench[i].give_dof_names(dofname, ntm);
      break;  
    case richardsmat:
      richar[i].give_dof_names(dofname, ntm);
      break;
    case salt1mat:
      salt1[i].give_dof_names(dofname, ntm);
      break;
      
      /*
	case salt2mat:
	//salt2[i].give_dof_names(dofname, ntm);
	break;
	
	case salt3mat:
	//salt3[i].give_dof_names(dofname, ntm);
	break;
	
	case salt4mat:
	//salt4[i].give_dof_names(dofname, ntm);
	break;
	
	case soilmat1:
	//soil1[i].give_dof_names(dofname, ntm);
	break;
	
	case lincoupledmat:
	//lcmat[i].give_dof_names(dofname, ntm);
	break;
      */      
    default:
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
  }
}



/**
   function reads integration points, material characteristics
   
   @param in - input stream
   
   21.7.2001
*/
void transmat::read (XFILE *in)
{
  //  computation of number of all integration points
  tnip = intpnum ();
  if (Mesprt==1){
    fprintf (stdout,"\n number of integration points  %ld",tnip);
  }
  
  //  reading of integration points
  //readip (in);
  //int. points allocation and initiation
  intpointalloc ();
  intpointinit ();
  
  //  reading of material characteristics
  readmatchar (in);
}

/**
   function prints integration points, material characteristics
   
   @param out - output stream
   
   18/12/2012, TKr
*/
void transmat::print (FILE *out)
{
  long i, j;
  //  printing of material characteristics
  fprintf(out,"\n## materials:\n");
  fprintf(out,"\n# number of material types\n");
  fprintf (out,"\n%ld\n",nmt);
  for (i=0;i<nmt;i++)
  {
    fprintf (out,"\n%d %ld", (int)mattype[i], numtype[i]);
    for (j=0;j<numtype[i];j++)
    {
      fprintf (out,"\n%ld ",j+1);
      printmatchar (out, mattype[i], j);
      fprintf(out, "\n");
    }
  }
}

/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   25.9.2001
*/
void transmat::matcond (matrix &d,long ipp,long ri,long ci)
{
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;
    //if (ip[ipp].infl[ri][ci]==1)
    m1.matcond(d,ri,ci,ipp);
    
    break;
  }
  case twomediacoup:{
    med2 m2;
    //if (ip[ipp].infl[ri][ci]==1)
    m2.matcond(d,ri,ci,ipp);
    
    break;
  }
  case threemediacoup:{
    med3 m3;
    //if (ip[ipp].infl[ri][ci]==1)
    m3.matcond(d,ri,ci,ipp);
    
    break;    
  }
  case fourmediacoup:{
    med4 m4;
    //if (ip[ipp].infl[ri][ci]==1)
    m4.matcond(d,ri,ci,ipp);
    
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }
  
  /*
  vector f1(d.m),f2(d.m);
  
  givefluxes (ri,ipp,fl);
  givefluxes (ci,ipp,f2);
  
  if (d.m==1){
    if (f1[0]<f2[0]){
      if (fabs(f1[0]/f2[0])<thres){
	d[0][0]=0.0;
      }
    }else{
      if (fabs(f2[0]/f1[0])<thres){
	d[0][0]=0.0;
      }
    }
    
  }
  if (d.m==2){
    if (f1[0]<f2[0]){
      if (fabs(f1[0]/f2[0])<thres && fabs(f1[1]/f2[1])<thres){
	d[0][0]=0.0;  d[0][1]=0.0;
	d[1][0]=0.0;  d[1][1]=0.0;
      }
    }else{
      if (fabs(f2[0]/f1[0])<thres && fabs(f2[1]/f1[1])<thres){
	d[0][0]=0.0;  d[0][1]=0.0;
	d[1][0]=0.0;  d[1][1]=0.0;
      }
    }
  }
  if (d.m==3){
    if (f1[0]<f2[0]){
      if (fabs(f1[0]/f2[0])<thres && fabs(f1[1]/f2[1])<thres && fabs(f1[2]/f2[2])<thres){
	d[0][0]=0.0;  d[0][1]=0.0;  d[0][2]=0.0;
	d[1][0]=0.0;  d[1][1]=0.0;  d[0][2]=0.0;
	d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=0.0;
      }
    }else{
      if (fabs(f2[0]/f1[0])<thres && fabs(f2[1]/f1[1])<thres && fabs(f2[2]/f1[2])<thres){
	d[0][0]=0.0;  d[0][1]=0.0;  d[0][2]=0.0;
	d[1][0]=0.0;  d[1][1]=0.0;  d[0][2]=0.0;
	d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=0.0;
      }
    }
  }
  */

}

/**
   function computes conductivity %matrix of the material second part - convective terms (only for three media coupling)
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   25.9.2001
*/
void transmat::matcond2 (matrix &d,long ipp,long ri,long ci)
{
  long gdim;
  if (Tp->advect==1){
    gdim = d.n;
    switch(gdim){
    case 1:{
      d[0][0]=givenontransq (advect_vel_x,ipp);
      break;
    }
    case 2:{
      d[0][0]=givenontransq (advect_vel_x,ipp);
      d[0][1]=givenontransq (advect_vel_y,ipp);
      break;
    }
    case 3:{
      d[0][0]=givenontransq (advect_vel_x,ipp);
      d[0][1]=givenontransq (advect_vel_y,ipp);
      d[0][2]=givenontransq (advect_vel_z,ipp);
      break;
    }
    default:{
      print_err("unknown geometrical dimension %ld is required",__FILE__,__LINE__,__func__,gdim);
      abort();
    }
    }
    return;
  }
  
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    m1.matcond2(d,ri,ci,ipp);
    break;
  }
  case twomediacoup:{
    med2 m2;

    m2.matcond2(d,ri,ci,ipp);
    
    break;
  }
  case threemediacoup:{
    med3 m3;

    m3.matcond2(d,ri,ci,ipp);
    
    break;
  }
   case fourmediacoup:{
    med4 m4;
    //if (ip[ipp].infl[ri][ci]==1)
    m4.matcond(d,ri,ci,ipp);

    break;
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes capacity %matrix of the material
   in the required integration point
   
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   25.9.2001

   @retval capcoeff - capacity coefficient
*/
double transmat::capcoeff (long ipp,long ri,long ci)
{
  double c=0.0;
   
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    m1.matcap(c,ri,ci,ipp);
        
    break;
  }
  case twomediacoup:{
    med2 m2;

    m2.matcap(c,ri,ci,ipp);
    
    break;
  }
  case threemediacoup:{
    med3 m3;

    m3.matcap(c,ri,ci,ipp);
   
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    m4.matcap(c,ri,ci,ipp);
   
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }

  return c;
}

/**
   function computes reaction %matrix of the material
   in the required integration point
   
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   JK, 12. 6. 2019

   @retval reactcoeff - reaction coefficient
*/
double transmat::reactcoeff (long ipp,long ri,long ci)
{
  double r=0.0;
  
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;
    
    m1.matreact(r,ri,ci,ipp);
    
    break;
  }
    /*
  case twomediacoup:{
    med2 m2;

    m2.matcap(c,ri,ci,ipp);
    
    break;
  }
  case threemediacoup:{
    med3 m3;

    m3.matcap(c,ri,ci,ipp);
   
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    m4.matcap(c,ri,ci,ipp);
   
    break;    
  }
  */
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }

  return r;
}


/**
   function computes correct fluxes from actual gradients

   @param lcid - load case id
   @param ipp - number of integration point
   
   4.12.2002
   Rewritten by Tomas Koudelka, 6.12.2013
*/
void transmat::computenlfluxes (long lcid,long ipp)
{
  long i,j,k,ncomp;
  ncomp = Tt->give_ncomp (0);
  matrix d(ASTCKMAT(ncomp,ncomp));
  matrix d2(ASTCKMAT(1,ncomp));

  for (i=0;i<ncomp;i++)//null flux
    ip[ipp].fluxes[lcid][i] = 0.0;

  for(i=0; i<Tp->ntm; i++){
    //  loop over the number of media/variables
    
    matcond(d, ipp, lcid, i);
    for (j=0;j<ncomp;j++){
      //  loop over the number of components of the flux vector
      
      //  diffusion contribution
      for (k=0;k<ncomp;k++)
        ip[ipp].fluxes[lcid][j] -= d[j][k]*ip[ipp].grad[i][k];
      
      //  advection contribution
      matcond2(d2,ipp,lcid,i);
      ip[ipp].fluxes[lcid][j] += d2[0][j]*ip[ipp].av[i];
    }    
  }
}

void transmat::computenlcapacities (long lcid,long ipp)
{
  long i;
  double c;

  ip[ipp].cap[lcid] = 0.0;

  for(i=0; i<Tp->ntm; i++){
    //  loop over the number of media/variables
    
    c=capcoeff (ipp, lcid, i);
    ip[ipp].cap[lcid] += c*(ip[ipp].av[i]-ip[ipp].pv[i]);
  }

}


/**
   function computes contributions to fluxes from particular driving forces

   for example: in heat and moisture transfer
   heat flux = k_TT grad T + k_Tw grad w
   moisture flux = k_wT grad T + k_ww grad w
   this function computes the fluxes k_TT grad T, k_Tw grad w,
   k_wT grad T and k_ww grad w
   it is used for adaptive modification of model
   if any of the fluxes is significantly smaller than others, it is
   neglected, the conductivity %matrix is also modified in such case
   
   @param lcid - load case id
   @param ipp - number of integration point
   
   JK, 9.10.2011
*/
void transmat::flux_contributions (long ipp)
{
  long i,j,ncomp;
  double n11,n12,n21,n22;
  //  the number of components (the space dimension, i.e. 1D, 2D or 3D)
  ncomp = Tt->give_ncomp (0);
  //  vectors are set to zero vectors in the constructor
  vector v11(ASTCKVEC(ncomp)),v12(ASTCKVEC(ncomp)),v21(ASTCKVEC(ncomp)),v22(ASTCKVEC(ncomp));
  vector gr(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp));
  
  switch (Tp->tmatt){//  type of transported matter
  case twomediacoup:{
    med2 m2;
    
    //  loop over the number of transported materials
    for (i=0;i<2;i++){
      for (j=0;j<2;j++){
	//  (i,j) block of the conductivity matrix
	m2.matcond(d,i,j,ipp);
	//  gradient of the j-th medium
	givegrad (j,ipp,gr);
	
	if (i==0 && j==0)
	  mxv (d,gr,v11);
	if (i==0 && j==1)
	  mxv (d,gr,v12);
	if (i==1 && j==0)
	  mxv (d,gr,v21);
	if (i==1 && j==1)
	  mxv (d,gr,v22);
      }
    }
    
    //  comparison of contributions
    n11 = normv (v11);
    n12 = normv (v12);
    n21 = normv (v21);
    n22 = normv (v22);

    break;
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
    abort ();
  }
  } 
}



/**
   function computes right-hand side %matrix of the material - contribution into right-hand side from volume integrals
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   @param ncomp - number of components
   
   25.9.2001
*/
void transmat::volume_rhs (matrix &d,long ipp,long ri,long ci,long /*ncomp*/)
{
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;
    
    m1.rhs_volume(d,ri,ci,ipp);
    break;
  }
  case twomediacoup:{
    med2 m2;
    
    m2.rhs_volume(d,ri,ci,ipp);
    
    /* switch (Tm->ip[ipp].tm){//material type
       case consolwf2:{
       long i;
       
       //  material id
       i = Tm->ip[ipp].idm;
       
       Tm->consol_wf2[i].rhs_volume (d,ri,ci,ipp);//tady predelat??!!
       break;
       }
       case consolawf2:{
       long i;
       
       //  material id
       i = Tm->ip[ipp].idm;
       
       Tm->consol_awf2[i].rhs_volume (d,ri,ci,ipp);//tady predelat??!!
       break;
       }
       case consolhwf2:{
       long i;
       
       //  material id
       i = Tm->ip[ipp].idm;
       
       Tm->consol_hwf2[i].rhs_volume (d,ri,ci,ipp);//tady predelat??!!
       break;
       }
       case bazantpedersen:
       case pedersen:
       case homomat:
       case kunzel:
       case grunewald:
       case moistheat:
       case salt1mat:
       break;
       default:{
       print_err("unknown material type is required",__FILE__,__LINE__,__func__);
       }
       }
    */
    break;
  }
  case threemediacoup:{
    med3 m3;

    m3.rhs_volume(d,ri,ci,ipp);
    
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    //m4.rhs_volume(d,ri,ci,ipp);
    
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function computes right-hand side ceofficient of the material 
   - contribution into right-hand side from volume integrals of the second type in the required integration point
   
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   25.9.2001

   @retval volume_rhs2 - volume coefficient
*/
double transmat::volume_rhs2 (long ipp,long ri,long ci)
{
  double c=0.0;
   
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    m1.rhs_volume2(c,ri,ci,ipp);
        
    break;
  }
  case twomediacoup:{
    med2 m2;
    
    m2.rhs_volume2(c,ri,ci,ipp);
    
    break;
  }
  case threemediacoup:{
    med3 m3;
    
    m3.rhs_volume2(c,ri,ci,ipp);
    
    break;    
  }
    /*     case fourmediacoup:{
	   med4 m4;
	   
	   m4.rhs_volume2(c,ri,ci,ipp);
	   
	   break;    
	   }
    */
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }

  return c;
}



/**
   function stores components of gradients into integration points
   
   @param lcid - load case id
   @param ipp - number of integration point
   @param gr - array containing components of gradient
   
   4.12.2002
*/
void transmat::storegrad (long lcid,long ipp,vector &gr)
{
  long i,ncomp;
  ncomp=gr.n;
  
  for (i=0;i<ncomp;i++){
    ip[ipp].grad[lcid][i]=gr[i];
  }
  
}

/**
   function returns components of gradients
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param gr - array containing components of flux

   22.12.2002
*/
void transmat::givegrad (long lcid,long ipp,vector &gr)
{
  long i,ncomp;
  ncomp=gr.n;
  
  for (i=0;i<ncomp;i++){
    gr[i] = ip[ipp].grad[lcid][i];
  }
}

/**
   function stores components of flux to integration points
   
   @param lcid - load case id
   @param ipp - number of integration point
   @param fl - array containing components of flux
   
   4.12.2002
*/
void transmat::storeflux (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    ip[ipp].fluxes[lcid][i]=fl[i];
  }
  
}

/**
   function returns components of flux
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param fl - array containing components of flux

   22.12.2002
*/
void transmat::givefluxes (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    fl[i] = ip[ipp].fluxes[lcid][i];
  }  
}


/**
  Function reads material characteristics from file.
   
  @param in - input file stream
   
  @return The function does not return anything.

  Created by JK, 8.7.2001
  Modified by TKo, TKr
  Modified by TKo 26.6.2014 - added keywords and split into readmatchar and readmattype
*/
void transmat::readmatchar (XFILE *in)
{
  long i;

  xfscanf (in, "%k%ld", "num_mat_types", &nmt);

  numtype = new long [nmt];
  mattype = new mattypet [nmt];  

  if (nmt<1)
    print_err("number of material types is less than 1",__FILE__,__LINE__,__func__);
  
  if (Mesprt==1)  fprintf (stdout,"\n number of different types of materials  %ld",nmt);

  for (i=0;i<nmt;i++)
  {
    xfscanf (in,"%k%m %k%ld", "mattype", &mattypet_kwdset, &mattype[i], "num_inst", &numtype[i]);
    if (numtype[i]<1)
      print_err("number of particular materials is less than 1",__FILE__,__LINE__,__func__);
    
    readmattype(in, mattype[i], numtype[i]);
  }
}



/**
  Function reads material characteristics from file.
   
  @param in - input file stream
  @param mtype - material type
  @param numt - the number of individual sets of material parameters
   
  @return The function does not return anything.

  Created by TKo according to original version of readmatchar, 26.6.2014
*/
void transmat::readmattype(XFILE *in, mattypet mtype, long numt)
{
  long j,k;

  switch (mtype)
  {
    case isotransmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of isotropic material models  %ld",numt);
      itrm = new isotrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of isotropic material is required",__FILE__,__LINE__,__func__);
	}
	itrm[k-1].read (in);
      }
      break;
    }
    case nlisotransmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of nonlinear isotropic material models  %ld",numt);
      nlitrm = new nlisotrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of nonlinear isotropic material is required",__FILE__,__LINE__,__func__);
	}
	nlitrm[k-1].read (in);
      }
      break;
    }
    case damisotransmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of isotropic material models influenced by damage %ld",numt);
      damitrm = new damisotrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of isotropic material influenced by damage is required",__FILE__,__LINE__,__func__);
	}
	damitrm[k-1].read (in);
      }
      break;
    }
    case homomat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models obtained from homogenization %ld",numt);
      hommat = new homogmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material obtained from homogenization is required",__FILE__,__LINE__,__func__);
	}
	//  this material reads only information about homogenized material at mezo level
	hommat[k-1].read (in);
      }
      break;
    }
    case tdisotransmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of time dependent isotropic material models  %ld",numt);
      tditrm = new tdisotrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of time dependent isotropic material is required",__FILE__,__LINE__,__func__);
	}
	tditrm[k-1].read (in);
      }
      break;
    }
    
    case interfacem:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for interface elements  %ld",numt);
      ifacemat = new interfacematt [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for interface elements is required",__FILE__,__LINE__,__func__);
	}
	ifacemat[k-1].read (in);
      }
      break;
    }

    case sejtkr:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numt);
      sejtkrm = new sejtkrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for one-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	sejtkrm[k-1].read (in);
      }
      break;
    }
     
    case consolawf1:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numt);
      consol_awf1 = new con_awf1mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for one-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_awf1[k-1].read (in);
      }
      break;
    }
     
    case consolwf1:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numt);
      consol_wf1 = new con_wf1mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for one-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_wf1[k-1].read (in);
      }
      break;
    }

    case consolwf2:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numt);
      consol_wf2 = new con_wf2mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_wf2[k-1].read (in);
      }
      break;
    }

    case consolawf2:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numt);
      consol_awf2 = new con_awf2mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_awf2[k-1].read (in);
      }
      break;
    }

    case consolhawf3:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for three-phase flow in deforming medium %ld",numt);
      consol_hawf3 = new con_hawf3mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for three-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_hawf3[k-1].read (in);
      }
      break;
    }


    case consolhwf2:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numt);
      consol_hwf2 = new con_hwf2mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_hwf2[k-1].read (in);
      }
      break;
    }

    case discontisotrmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of isotropic material models with discontinuities %ld",numt);
      ditrm = new discisotrmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of isotropic material with discontinuities is required",__FILE__,__LINE__,__func__);
	}
	ditrm[k-1].read (in);
      }
      break;
    }
      
      
    case cernyconcrete:{
      if (Mesprt==1)  fprintf (stdout,"\n number of cerny-concrete material models  %ld",numt);
      cernym = new cernymat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of cerny-concrete material is required",__FILE__,__LINE__,__func__);
	}
	cernym[k-1].read (in);
      }
      break;
    }
      
    case bazantpedersen:{
      if (Mesprt==1)  fprintf (stdout,"\n number of bazant-pedersen material models %ld",numt);
      bazped = new bazpedmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of bazant-pedersen material model is required",__FILE__,__LINE__,__func__);
	}
	bazped[k-1].read (in);
      }
      break;
    }

    case pedersen:{
      if (Mesprt==1)  fprintf (stdout,"\n number of pedersen material models %ld",numt);
      ped = new pedmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of pedersen material model is required",__FILE__,__LINE__,__func__);
	}
	ped[k-1].read (in);
      }
      break;
    }
      
    case kunzel:{
      if (Mesprt==1)  fprintf (stdout,"\n number of kunzel material models %ld",numt);
      kun = new kunmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of kunzel material model is required",__FILE__,__LINE__,__func__);
	}
	kun[k-1].read (in);
      }

      //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;


      break;
    } 
      
    case kunzel2:{
      if (Mesprt==1)  fprintf (stdout,"\n number of kunzel2 material models %ld",numt);
      kun2 = new kunmat2 [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of kunzel2 material model is required",__FILE__,__LINE__,__func__);
	}
	kun2[k-1].read (in);
      }
	   //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;


      break;
    } 


    case grunewald:{
      if (Mesprt==1)  fprintf (stdout,"\n number of grunewald material models %ld",numt);
      grunw = new grunewaldmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of grunewald material model is required",__FILE__,__LINE__,__func__);
	}
	grunw[k-1].read (in);
      }
	  //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;

      break;
    } 

    case simplediscmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of simple discontinous material models for moisture %ld",numt);
      sdmat = new discmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of simple discontinuous material model for moisture is required",__FILE__,__LINE__,__func__);
	}
	sdmat[k-1].read (in);
      }
      //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;
      
      break;
    } 
      
    case devries:{
      if (Mesprt==1)  fprintf (stdout,"\n number of devries material models %ld",numt);
      dvries = new devriesmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of devries material model is required",__FILE__,__LINE__,__func__);
	}
	dvries[k-1].read (in);
      }
      break;
    }

    case milly:{
      if (Mesprt==1)  fprintf (stdout,"\n number of Milly material models %ld",numt);
      mill = new millymat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of Milly material model is required",__FILE__,__LINE__,__func__);
	}
	mill[k-1].read (in);
      }
        //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;
      break;
    }

      
    case concreteB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of general material models for concrete at high temperature from Padova%ld",numt);
      concrete = new concreteBmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of general material model for concrete at high temperature from Padova is required",__FILE__,__LINE__,__func__);
	}
	concrete[k-1].read (in);
      }
      break;
    }

    case baroghelB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of general baroghel material models for concrete  from Padova%ld",numt);
      baroghel = new baroghelmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of general baroghel material model for concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	baroghel[k-1].read (in);
      }
      break;
    }
      

    case C60baroghelB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for high performance C60 concrete  from Padova%ld",numt);
      C60baroghel = new C60barmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of baroghel material model for high performance C60 concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	C60baroghel[k-1].read (in);
      }
      break;
    }

    case C30baroghelB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for ordinary performance C30 concrete  from Padova%ld",numt);
      C30baroghel = new C30barmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of baroghel material model for ordinary performance C30 concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	C30baroghel[k-1].read (in);
      }
      break;
    }


    case o30bazantB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of bazant material models for ordinary-30 concrete  from Padova%ld",numt);
      o30bazant = new o30bazmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of baroghel material model for ordinary-30 concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	o30bazant[k-1].read (in);
      }
      break;
    }

      
    case C60bazantB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of bazant material models for high performance C60 concrete  from Padova%ld",numt);
      C60bazant = new C60bazmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of bazant material model for high performance C60 concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	C60bazant[k-1].read (in);
      }
      break;
    }

      
    case C30bazantB:{
      if (Mesprt==1)  fprintf (stdout,"\n number of bazant material models for high performance C60 concrete  from Padova %ld",numt);
      C30bazant = new C30bazmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of bazant material model for high performance C30 concrete from Padova is required",__FILE__,__LINE__,__func__);
	}
	C30bazant[k-1].read (in);
      }
      break;
    }


    case glasgow:{
      if (Mesprt==1)  fprintf (stdout,"\n number of Glasgow material models for concrete %ld",numt);
      tench = new glasgowmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of Glasgow material model for concrete is required",__FILE__,__LINE__,__func__);
	}
	tench[k-1].read (in);
      }
      break;
    }
    case carb1mat:{
      if (Mesprt==1) fprintf (stdout,"\n number of material models for concrete carbonation %ld", numt);
      carb1 = new carbmat1 [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material model for carbonation is required",__FILE__,__LINE__,__func__);
	}
	carb1[k-1].read (in);
      }
      break;
	}
	case moistheat:{
	  if (Mesprt==1)  fprintf (stdout,"\n number of moistheat material models %ld",numt);
	  moisth = new moistheatmat [numt];
	  for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of moistheat material model is required",__FILE__,__LINE__,__func__);
	}
	moisth[k-1].read (in);
	  }

	  //  this material model needs actual values stored at nodes
	  Tp->nvs=1;
	  //  this material model needs previous values stored at nodes
	  Tp->pnvs=1;


	  break;
	}
      
    case richardsmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of Richards models  %ld",numt);
      richar = new richards [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of Richards material model is required",__FILE__,__LINE__,__func__);
	}
	richar[k-1].read (in);
      }
      break;
    }
      
    case salt1mat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for salt transport %ld",numt);
      salt1 = new saltmat1 [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material model for salt transport is required",__FILE__,__LINE__,__func__);
	}
	salt1[k-1].read (in);
      }
      break;
    }
      
    case salt2mat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for salt transport%ld",numt);
      salt2 = new saltmat2 [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material model for salt transport is required",__FILE__,__LINE__,__func__);
	}
	salt2[k-1].read (in);
      }
      break;
    }

    case salt3mat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for salt transport%ld",numt);
      salt3 = new saltmat3 [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material model for salt transport is required",__FILE__,__LINE__,__func__);
	}
	salt3[k-1].read (in);
      }
      break;
    }
      
    case salt4mat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for salt transport%ld",numt);
      salt4 = new saltmat4 [numt];

      //  this material model needs actual values stored at nodes
      Tp->nvs=1;
      //  this material model needs previous values stored at nodes
      Tp->pnvs=1;

      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material model for salt transport is required",__FILE__,__LINE__,__func__);
	}
	salt4[k-1].read (in);
      }
      break;
    }
      
    case radiationmater:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for radiation %ld",numt);
      radmat = new radiationmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of material models for radiation is required",__FILE__,__LINE__,__func__);
	}
	radmat[k-1].read (in);
      }
      break;
    }
      
    case soilmat1:{
      if (Mesprt==1)  fprintf (stdout,"\n number of soil1mat material models for soils %ld",numt);
      soil1 = new soil1mat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of soil1mat material model for soils is required",__FILE__,__LINE__,__func__);
	}
	soil1[k-1].read (in);
      }
      break;
    }
      
    case cementhydrmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of cemhydmat material models for cement hydration %ld",numt);
      cemhydr = new cemhydmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of cemhydmat material model for cement hydration is required",__FILE__,__LINE__,__func__);
	}
	cemhydr[k-1].read (in);
      }
      break;
    }

    case lincoupledmat:{
      if (Mesprt==1)  fprintf (stdout,"\n number of linear coupled material models %ld",numt);
      lcmat = new lincoupmat [numt];
      for (j=0;j<numt;j++){
	k=numt+1;
	xfscanf (in,"%ld",&k);
	if (k>numt || k<1){
	  print_err("wrong number of linear coupled material model is required",__FILE__,__LINE__,__func__);
	}
	lcmat[k-1].read (in);
      }
      break;
    }

    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
  }
}



/**
   function prints material characteristics
   
   @param out - output stream
   @param mt  - type of printed material
   @param numinst - index of set of parameters for the given material type
   
   TKr 3.1.2006
   Modified by Tomas Koudelka, 4.7.2014
*/
void transmat::printmatchar (FILE *out, mattypet mt, long numinst)
{
  switch (mt)
  {
    case isotransmat:{
      itrm[numinst].print (out);
      break;
    }
    case nlisotransmat:{
      nlitrm[numinst].print (out);
      break;
    }
    case homomat:{
      hommat[numinst].print (out);
      break;
    }
    case damisotransmat:{
      damitrm[numinst].print (out);
      break;
    }
    case discontisotrmat:{
      ditrm[numinst].print (out);
      break;
    }
    case tdisotransmat:{
      tditrm[numinst].print (out);
      break;
    }

    case interfacem:{
      ifacemat[numinst].print (out);
      break;
    }
      
    case cernyconcrete:{
      cernym[numinst].print (out);
      break;
    }
      
    case bazantpedersen:{
      bazped[numinst].print (out);
      break;
    }

    case pedersen:{
      ped[numinst].print (out);
      break;
    }
      
    case kunzel:{
      kun[numinst].print (out);
      break;
    } 

    case kunzel2:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //kun2[numinst].print (out);
      break;
    } 
      
    case grunewald:{
      grunw[numinst].print (out);
      break;
    } 

    case simplediscmat:{
      sdmat[numinst].print (out);
      break;
    } 

    case devries:{
      dvries[numinst].print (out);
      break;
    } 
      
    case milly:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //mill[numinst].print (out);
      break;
    }
      
    case concreteB:{
      concrete[numinst].print (out);
      break;
    }
      
    case baroghelB:{
      baroghel[numinst].print (out);
      break;
    }
      

    case C60baroghelB:{
      C60baroghel[numinst].print (out);
      break;
    }

    case C30baroghelB:{
      C30baroghel[numinst].print (out);
      break;
    }


    case o30bazantB:{
      o30bazant[numinst].print (out);
      break;
    }

      
    case C60bazantB:{
      C60bazant[numinst].print (out);
      break;
    }

      
    case C30bazantB:{
      C30bazant[numinst].print (out);
      break;
    }


    case glasgow:{
      tench[numinst].print (out);
      break;
	}

	case moistheat:{
	  moisth[numinst].print (out);
	  break;
	}
      
    case richardsmat:{
      richar[numinst].print (out);
      break;
    }
      
    case consolawf1:
      consol_awf1[numinst].print(out);
      break;
    case consolwf1:
      consol_wf1[numinst].print(out);
      break;
    case consolwf2:
      consol_wf2[numinst].print(out);
      break;
    case consolawf2:
      consol_awf2[numinst].print(out);
      break;
    case consolhwf2:
      consol_hwf2[numinst].print(out);
      break;
    case salt1mat:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //salt1[numinst].print (out);
      break;
    }

    case salt2mat:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //salt2[numinst].print (out);
      break;
    }

    case salt3mat:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //salt3[numinst].print (out);
      break;
    }

    case salt4mat:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //salt4[numinst].print (out);
      break;
    }

    case soilmat1:{
      print_err("the print method has not yet been implemented for material %s",
                __FILE__,__LINE__,__func__, mattypet_kwdset.get_str(int(mt)));
      //soil1[numinst].print (out);
      break;
    }

    case lincoupledmat:{
      lcmat[numinst].print (out);
      break;
    }
      
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
  }
}


/**
  The function stores all components of array comp to integration point array other, 
  first position in the array other is fi.

  @param ipp - integration point pointer
  @param fi - first index
  @param ncomp - number of required components
  @param comp - %vector containing components
   
  @return The function does not return anything.

  5.11.2002
  Modified by TKo, 19.12.2017
*/
void transmat::storeother (long ipp, long fi, long ncomp,double *comp)
{
  long i, j;

  for (i=fi, j=0; i<fi+ncomp; i++, j++)
    ip[ipp].other[i]=comp[j];
}


/**
  The function stores one component to integration point array other, 
  position in the array other is given by id.

  @param ipp - integration point pointer
  @param id - component index
  @param comp - component value
   
  @return The function does not return anything.

  Created by TKo, 19.12.2017
*/
void transmat::storeother (long ipp, long id, double comp)
{
  ip[ipp].other[id]=comp;
}



/**
  The function stores all components of array comp to integration point array eqother, 
  first position in the array eqother is fi.

  @param ipp - integration point pointer
  @param fi - first index
  @param ncomp - number of required components
  @param comp - %vector containing components
   
  @return The function does not return anything.

  31.1.2007
  Modified by TKo, 19.12.2017
*/
void transmat::storeeqother (long ipp, long fi, long ncompeq, double *compeq)
{
  long i, j;

  for (i=fi, j=0; i<fi+ncompeq; i++, j++)
    ip[ipp].eqother[i]=compeq[j];
}


/**
  The function restores part of array other to the array comp.
  It selects ncomp components from the integration point array other starting from fi-th component, 
  and stores them in comp array.

  @param ipp - integration point pointer
  @param fi - first index
  @param ncomp - number of required components
  @param comp - %vector containing components
   
  @return The function does not return anything.
  5.11.2002
  Modified by TKo, 19.12.2017
*/
void transmat::giveother(long ipp, long fi, long ncomp, double *comp)
{
  long i, j;
  
  for (i=fi, j=0; i<fi+ncomp; i++, j++){
    comp[j]=ip[ipp].other[i];
  }
}


/**
   function restores array eq_other to %vector comp
   (function selects all components of array eqcomp from
   the array eqother, first index in array eqother is fi)
   
   @param ipp - integration point pointer
   @param ncomp - number of components
   @param compeq - %vector containing components
   
   5.11.2002
   Modified by TKo, 19.12.2017
*/
void transmat::giveeqother(long ipp, long fi, long ncompeq, double *compeq)
{
  long i, j;
  
  for (i=fi, j=0; i<fi+ncompeq; i++, j++){
    compeq[j]=ip[ipp].eqother[i];
  }
}




/**
   This function returns the number of components of ipp's other array.
   
   TK, 5.11.2003
*/
long transmat::givencompother ()
{
  long ncompo = 0;
  
  switch (Tp->mednam)//names of transported media
    {
    case water:{
      ncompo=5;
      break;
    }
    case moisture:{
      ncompo=3;
      if(Tp->ntm==2)
	ncompo=6;
      break;
    }
    case heat:
      ncompo=1;
      break;
    case heat_moisture:
      if (Tp->tmatt == twomediacoup)
	ncompo=8;//twomediacoup
      else
	ncompo=9;//threemediacoup
      break;
    case moisture_salt:
      if (Tp->tmatt == twomediacoup)
	ncompo=3;//twomediacoup
      else
	ncompo=8;//threemediacoup
      break;
    case moisture_salt_crystal:
      ncompo=8;//threemediacoup
      break;
    case moisture_salt_crystal_heat:
      ncompo=8;//fourmediacoup
      break;
    default:
      print_err("unknown problem name is required",__FILE__,__LINE__,__func__);
    }
  return ncompo;
}



/**
   This function returns the number of components of ipp's eqother array.
   Material type is given by the ipp parameter, which means index of integration point whose 
   material type is being considred.
   
   @param ipp - integration point id
   @param im - index of material // ???!!! candidate for removal
   
   @return The function returns the number of eqother array components at the given integration point.

   TKr, 12.12.2006
*/
long transmat::givencompeqother (long ipp, long /*im*/)
{
  long ncompeqo=0;  
  mattypet tm = ip[ipp].tm;

  ncompeqo = givencompeqothermat(tm);

  return ncompeqo;
}



/**
   This function returns the number of components of ipp's eqother array.
   Material type is given by the tm parameter.
   
   @param tm - matrial type
   
   @return The function returns the number of eqother array components of the given material type.

   TKo, 25.01.2019
*/
long transmat::givencompeqothermat (mattypet tm)
{
  long ncompeqo=0;  

  switch (tm){
    case isotransmat:
    case nlisotransmat:
    case discontisotrmat:
    case tdisotransmat:
    case interfacem:
    case cernyconcrete:
    case bazantpedersen:
    case pedersen:
    case carb1mat:
    case sejtkr:
      break;
    case consolawf1:
      ncompeqo=2;
      break;
    case consolwf1:
      ncompeqo=2;
      break;
    case lincoupledmat:
    case radiationmater:
    case richardsmat:
      break;
    case damisotransmat:
      ncompeqo=1;
      break;
    case cementhydrmat:
      ncompeqo=6;
      break;
    case homomat:{
      ncompeqo = 15;
      //temporarilly ??!!
      switch (Tp->tmatt){//names of transported media
      case onemedium:{
	ncompeqo = 1;
	break;
      }
      case twomediacoup:{
	if(Tm->kun != NULL)//kunzel mat. type
	  ncompeqo = 5;
	if(Tm->kun2 != NULL)//kunzel2 mat. type
	  ncompeqo = 15;
	if(Tm->moisth != NULL)//moistheat mat. type
	  ncompeqo = 5;
	break;
      }
      default:
	print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
      }
      break;
    }
    case kunzel:{
      ncompeqo = 5;
      break;
    }
    case kunzel2:{
      ncompeqo = 15;
      break;
    }
    case grunewald:{
      ncompeqo = 15;
      break;
    }
    case simplediscmat:{
      ncompeqo = 15;
      break;
    }
    case devries:{
      break;
    }
    case milly:{
      ncompeqo = 10;
      break;
    }
    case moistheat:{
      ncompeqo = 8;
      break;
    }
    case salt1mat:{
      ncompeqo = 1;
      break;
    }
    case salt3mat:{
      ncompeqo = 4;
      break;
    }
    case salt4mat:{
      ncompeqo = 9;
      break;
    }
    case glasgow:
    case salt2mat:
    case concreteB:
    case baroghelB:
    case C60baroghelB:
    case C30baroghelB:
    case o30bazantB:
    case C30bazantB:
    case C60bazantB:
    case soilmat1:
      break;
    case consolwf2:
      ncompeqo = 2;
      break;
    case consolawf2:
      ncompeqo = 2;
      break;
    case consolhawf3:
      ncompeqo = 2;
      break;
    case consolhwf2:
      ncompeqo = 2;
      break;
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
  }
  
  return ncompeqo;
}


/**
   This function returns the components of other array for unknowns.
   
   @param compother - component of other array
   @param ipp       - integration point pointer
   @param r         - vector of unknowns
   
   TKr, 28.6.2004
*/
double transmat::givecompother(long compother,long ipp,double *r)
{
  double compo;
  
  switch (Tp->tmatt){//names of transported media
  case onemedium:{
    med1 m1;

    compo = m1.compute_othervalues(compother,ipp,r);
    break;
  }
  case twomediacoup:{
    med2 m2;
    
    compo = m2.compute_othervalues(compother,ipp,r);
    break;
  }
  case threemediacoup:{
    med3 m3;
    
    compo = m3.compute_othervalues(compother,ipp,r);
    break;
  }
  case fourmediacoup:{
    med4 m4;
    
    compo = m4.compute_othervalues(compother,ipp,r);
    break;
  }
  default:
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  return (compo);
}



/**
   This function prints othervalue name.

   @param out       - output file
   @param compother - component of other array
   @param ipp       - integration point pointer
   
   TKr, 28.6.2004
*/
void transmat::give_othervalue_name(FILE *out,long ipp,long compother)
{
  switch (Tp->tmatt){//names of transported media
  case onemedium:{
    med1 m1;
    
    m1.print_othervaluesnames (out,ipp,compother);    
    break;
  }
  case twomediacoup:{
    med2 m2;
    
    m2.print_othervaluesnames (out,ipp,compother);    
    break;
  }
  case threemediacoup:{
    med3 m3;
    
    m3.print_othervaluesnames (out,ipp,compother);
    break;
  }
  case fourmediacoup:{
    med4 m4;
    
    m4.print_othervaluesnames (out,ipp,compother);
    break;
  }
  default:
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
}


/**
   This function prints eqothervalue name.

   @param out         - output file
   @param compeqother - component of eqother array
   @param ipp         - integration point pointer
   
   TKr, 28.6.2004
*/
void transmat::give_eqothervalue_name(FILE *out,long ipp,long compeqother)
{
  give_othervalue_name(out, ipp, compeqother);
}



/**
   function returns temperature in integration point
   
   @param ipp - integration point id
   
   JK, 24. 10. 2013
*/
double transmat::give_temperature (long ipp)
{
  long idm;
  double t;
  
  idm = ip[ipp].idm;
  
  switch (ip[ipp].tm){
  case isotransmat:{
    t = itrm[idm].give_temperature (ipp);
    break;
  }
  case consolhawf3:{
    t = consol_hawf3[idm].give_temperature (ipp);
    break;
  }
  case consolhwf2:{
    t = consol_hwf2[idm].give_temperature (ipp);
    break;
  }
  case kunzel:{
    t = kun[idm].give_temperature (ipp);
    break;
  }
  case moistheat:{
    t=moisth[idm].give_temperature (ipp);
    break;
  }
  case salt4mat:{
    t = salt4[idm].give_temperature (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required on ip %ld, eid=%ld",__FILE__,__LINE__,__func__, ipp, Tm->elip[ipp]+1);
  }
  }
  
  return t;
}



/**
   function returns temperature in integration point
   
   @param ipp - integration point id
   
   TKr, 17/07/2018
*/
double transmat::give_inittemperature(long ipp)
{
  long idm;
  double t;
  
  idm = ip[ipp].idm;
  
  switch (ip[ipp].tm){
  case isotransmat:{
    t = itrm[idm].give_temperature (ipp);
    break;
  }
  case consolhawf3:{
    t = consol_hawf3[idm].give_temperature (ipp);
    break;
  }
  case consolhwf2:{
    t = consol_hwf2[idm].give_temperature (ipp);
    break;
  }
  case kunzel:{
    t = kun[idm].give_temperature (ipp);
    break;
  }
    case moistheat:{
    t = moisth[idm].give_temperature (ipp);
    break;
  }
  case salt4mat:{
    t = salt4[idm].give_temperature (ipp);
    break;
  }
  default:{
    print_err("unknown material type %d is required on ip %ld, eid=%ld",__FILE__,__LINE__,__func__, ipp, Tm->elip[ipp]+1);
  }
  }
  
  return t;
}


/**
  Function returns required transport quantity stored at integration point
  for MEFEL.
   
  @param qt - type of required quantity
  @param ipp - number of integration point
   
  @return The function returns required %vector of quantitiy in the parameter ipv.
 
  Created by Tomas Krejci according to Tomas Koudelka, 14/10/2013
*/
double transmat::givetransq(nonmechquant qt, long ipp)
{
  double ret=0.0;

  switch (qt){
  case temperature:{
    //  temperature from integration point
    ret = give_temperature(ipp);
    break;
  }
  case rel_hum:{
    //  relative humidity from integration point
    //ret = give_rel_hum(ipp);
    break;
  }
  case initial_temperature:{
    //  initial temperature from integration point
    ret = give_inittemperature(ipp);
    break;
  }
  case water_pressure:{
    //  water pressure from integration point
    ret = give_water_pressure(ipp);
    break;
  }
  case gas_pore_pressure:{
    //  water pressure from integration point
    ret = give_gas_pressure(ipp);
    break;
  }
  case cap_pressure:{
    //  capillary pressure from integration point
    //ret = give_cap_pressure(ipp);
    break;
  }
  case saturation_degree:{
    //  saturation degree from integration point
    ret = give_saturation_degree(ipp);
    break;
  }
  case suction:{
    //  suction pressure from integration point
    ret = give_suction(ipp);
    break;
  }
  case eff_pore_pressure:{
    //  pore pressure from integration point
    ret = give_effective_pore_pressure(ipp);
    break;
  }
  case pore_pressure:{
    //  pore pressure from integration point
    ret = give_pore_pressure(ipp);
    break;
  }
  case vol_moist_cont:{
    // volumetric moisture content
    ret = give_vol_moist_cont(ipp);
    break;
  }
  case volume_change:{
    // volumetric change
    ret = give_volume_change (ipp);
    break;
  }
  default:{
    print_err("unknown type of quantity is required",__FILE__,__LINE__,__func__);
  }
  }

  return ret;
}


/**
  Function returns water pressure especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of water pressure.
  
  Created by Tomas Krejci according to Tomas Koudelka 14/10/2013
*/
double transmat::give_water_pressure(long ipp)
{
  long i;
  double pw;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolawf1:{
    pw = consol_awf1[i].give_water_pressure (ipp);
    break;
  }
  case consolwf1:{
    pw = consol_wf1[i].give_water_pressure (ipp);
    break;
  }
  case consolwf2:{
    pw = consol_wf2[i].give_water_pressure (ipp);
    break;
  }
  case consolawf2:{
    pw = consol_awf2[i].give_water_pressure (ipp);
    break;
  }
  case consolhawf3:{
    pw = consol_hawf3[i].give_water_pressure (ipp);
    break;
  }
  case consolhwf2:{
    pw = consol_hwf2[i].give_water_pressure (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return pw;
}


/**
  Function returns gas (air) pressure especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of water pressure.
  
  Created by Tomas Krejci according to Tomas Koudelka 27/05/2016
*/
double transmat::give_gas_pressure(long ipp)
{
  long i;
  double pg;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolawf1:{
    pg = 0.0;
    break;
  }
  case consolwf1:{
    pg = 0.0;
    break;
  }
  case consolwf2:{
    pg = consol_wf2[i].give_gas_pressure (ipp);
    break;
  }
  case consolawf2:{
    pg = consol_awf2[i].give_gas_pressure (ipp);
    break;
  }
  case consolhawf3:{
    pg = consol_hawf3[i].give_gas_pressure (ipp);
    break;
  }
  case consolhwf2:{
    pg = 0.0;
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return pg;
}

/**
  Function returns pore pressure especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of pore pressure.
  
  Created by Tomas Krejci according to Tomas Koudelka 14/10/2013
*/
double transmat::give_pore_pressure(long ipp)
{
  long i;
  double pp;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolawf1:{
    pp = consol_awf1[i].give_pore_pressure (ipp);
    break;
  }
  case consolwf1:{
    pp = consol_wf1[i].give_pore_pressure (ipp);
    break;
  }
  case consolwf2:{
    pp = consol_wf2[i].give_pore_pressure (ipp);
    break;
  }
  case consolawf2:{
    pp = consol_awf2[i].give_pore_pressure (ipp);
    break;
  }
  case consolhawf3:{
    pp = consol_hawf3[i].give_pore_pressure (ipp);
    break;
  }
  case consolhwf2:{
    pp = consol_hwf2[i].give_pore_pressure (ipp);
    break;
  }
  case salt4mat:{
    pp = salt4[i].give_pore_pressure (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return pp;
}

/**
  Function returns effective pore pressure especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of pore pressure.
  
  Created by Tomas Krejci according to Tomas Koudelka 14/10/2013
*/
double transmat::give_effective_pore_pressure(long ipp)
{
  long i;
  double pp;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolwf1:{
    pp = consol_wf1[i].give_effective_pore_pressure (ipp);
    break;
  }
  case consolawf1:{
    pp = consol_awf1[i].give_effective_pore_pressure (ipp);
    break;
  }
  case consolawf2:{
    pp = consol_awf2[i].give_effective_pore_pressure (ipp);
    break;
  }
  case consolhwf2:{
    pp = consol_hwf2[i].give_effective_pore_pressure (ipp);
    break;
  }
  case consolhawf3:{
    pp = consol_hawf3[i].give_effective_pore_pressure (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return pp;
}

/**
  Function returns suction especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of suction.
  
  Created by Tomas Krejci according to Tomas Koudelka 16/12/2013
*/
double transmat::give_suction(long ipp)
{
  long i;
  double s;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolawf1:{
    s = consol_awf1[i].give_suction (ipp);
    break;
  }
  case consolwf1:{
    s = consol_wf1[i].give_suction (ipp);
    break;
  }
  case consolwf2:{
    s = consol_wf2[i].give_suction (ipp);
    break;
  }
  case consolawf2:{
    s = consol_awf2[i].give_suction (ipp);
    break;
  }
  case consolhawf3:{
    s = consol_hawf3[i].give_suction (ipp);
    break;
  }
  case consolhwf2:{
    s = consol_hwf2[i].give_suction (ipp);
    break;
  }
  default:{
    print_err("unknown material type %d is required in ipp=%ld",__FILE__,__LINE__,__func__,Tm->ip[ipp].tm,ipp);
  }
  }
  
  return s;
}

/**
  Function returns saturation degree especially for models of soils.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of saturation degree
  
  Created by Tomas Krejci according to Tomas Koudelka 16/12/2013
*/
double transmat::give_saturation_degree(long ipp)
{
  long i;
  double s;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case consolawf1:{
    s = consol_awf1[i].give_saturation_degree(ipp);
    break;
  }
  case consolwf1:{
    s = consol_wf1[i].give_saturation_degree(ipp);
    break;
  }
  case consolwf2:{
    s = consol_wf2[i].give_saturation_degree(ipp);
    break;
  }
  case consolawf2:{
    s = consol_awf2[i].give_saturation_degree(ipp);
    break;
  }
  case consolhawf3:{
    s = consol_hawf3[i].give_saturation_degree(ipp);
    break;
  }
  case consolhwf2:{
    s = consol_hwf2[i].give_saturation_degree(ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return s;
}

/**
  Function returns volume change for individual material models.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of volumetric mositure content
  
  Created by Tomas Koudelka 12/02/2019
*/
double transmat::give_volume_change (long ipp)
{
  long i;
  double h;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case moistheat:{
    h = moisth[i].give_volume_change (ipp);
    break;
  }
  default:{
    print_err("unknown material type %d is required",__FILE__,__LINE__,__func__, Tm->ip[ipp].tm);
  }
  }
  
  return h;
}


/**
  Function returns volumetric mositure content for individual material models.
  It is used for transfer required tranpsort values between MEFEL and TRFEL.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of volumetric mositure content
  
  Created by Tomas Koudelka 4/02/2014
*/
double transmat::give_vol_moist_cont(long ipp)
{
  long i;
  double h;
  
  i = ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case isotransmat:
    h = itrm[i].give_vol_moist(ipp);
    break;
  case kunzel:{
    h = kun[i].give_vol_moist(ipp);
    break;
  }
  case moistheat:{
    h = moisth[i].give_vol_moist(ipp);
    break;
  }
  case interfacem:
    h = ifacemat[i].give_vol_moist(ipp);
    break;
  default:{
    print_err("unknown material type %d is required",__FILE__,__LINE__,__func__, Tm->ip[ipp].tm);
  }
  }
  
  return h;
}


/**
  Function returns volumetric moisture content (m3/m3) for individual material models.
  
  @param nid - node id
  @param mattype - type of material

  @return The function returns actual value of volumetric moisture content.
  
  JK, 10. 4. 2014
*/
double transmat::give_nodal_vol_moist_cont (long nid,long mattype)
{
  double w;
  //  long eid,ipp,i,idm;
  //  double *in,*inp, *ineq;
  
  switch (mattype){
  case 154:{
    w = Tt->nodes[nid].eqother[1];
    break;
  }
  case 155:{
    w = Tt->nodes[nid].eqother[0];
    break;
  }
  case 158:{
    w = Tt->nodes[nid].eqother[0];
    break;
  }
  case 180:{
    w = Tt->nodes[nid].eqother[0];
    break;
  }
  case 110:{
    w = 0.0;

    //////////////////////////////////////////
    //zde je nutne predelat??!!
    /*
    //it must be completed here for all possibilities
    //number of element containing node number nid
    //eid = nodel;
    //number of the first integration point of the element
    //ipp=Tt->elements[eid].ipp[0][0];
    //temporarilly??!!
    ipp=Tt->elements[0].ipp[0][0];
    //material number
    i = ip[ipp].idm;
    
    fprintf (stdout,"\n\n i = %ld\n",i);
    
    if(hommat[i].hom_mattype != mattype){
    print_err("Basic material type for homogenization is not corresponding with material type in climatic conditions",__FILE__,__LINE__,__func__);
    abort();
    }
    
    switch (hommat[i].hom_mattype){
    case 180:{
    in = new double[2];
    inp = new double[2];
    ineq = new double[5];
    Tt->nodes[nid].give_values (in,inp,ineq);
    idm = hommat[i].hom_mattype_number;
    
    //fprintf (stdout,"\n\n idm = %ld\n",idm);
    
    if(Tm->moisth == NULL){
    print_err("Basic material for homogenization is not defined in macro problem",__FILE__,__LINE__,__func__);
    abort();
    }
    
    w = moisth[idm].sorpiso.isotherm_value (in[0]); // volumetric moisture content from moisture storage function
    
    fprintf (stdout,"\n\n w = %ld\n",w);
    break;
    }
    default:{
    print_err("unknown type of material for homogenization",__FILE__,__LINE__,__func__);
    }
    }
    */
    //////////////////////////////////////////
    
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mattyp)
  
  return w;
}

/**
  Function returns saturated volumetric moisture content (m3/m3) for individual material models.
  
  @param nid - node id
  @param mattype - type of material
  
  @return The function returns actual value of volumetric moisture content.
  
  JK, 10. 4. 2014
*/
double transmat::give_nodal_sat_vol_moist_cont (long nid,long mattype)
{
  double w;
  //  long eid,ipp,i,idm;
  //  double *in,*inp, *ineq;
  
  switch (mattype){
  case 150:{
    w = 10000;
    break;
  }
  case 154:{
    w = Tt->nodes[nid].eqother[3];
    break;
  }
  case 155:{
    w = Tt->nodes[nid].eqother[2];
    break;
  }
  case 156:{
    w = Tt->nodes[nid].eqother[2];
    break;
  }
  case 158:{
    w = Tt->nodes[nid].eqother[3];
    break;
  }
  case 159:{
	w = Tt->nodes[nid].eqother[3];
	break;
  }
  case 180:{
	w = Tt->nodes[nid].eqother[2];
	break;
  }
  case 203:{
    w = Tt->nodes[nid].eqother[3];
    break;
  }
  case homomat:{
    //////////////////////////////////////////
    //zde je nutne predelat??!!
    /*
    //it must be completed here for all possibilities
    //number of element containing node number nid
    //eid = nodel;
    //number of the first integration point of the element
    //ipp=Tt->elements[eid].ipp[0][0];
    //temporarilly??!!
    ipp=Tt->elements[0].ipp[0][0];
    //material number
    i = ip[ipp].idm;
    
    if(hommat[i].hom_mattype != mattype){
    print_err("Basic material type for homogenization is not corresponding with material type in climatic conditions",__FILE__,__LINE__,__func__);
    abort();
    }
    
    switch (hommat[i].hom_mattype){
    case 180:{
    in = new double[2];
    inp = new double[2];
    ineq = new double[5];
    Tt->nodes[nid].give_values (in,inp,ineq);
    idm = hommat[i].hom_mattype_number;
    
    if(Tm->moisth == NULL){
    print_err("Basic material for homogenization is not defined in macro problem",__FILE__,__LINE__,__func__);
    abort();
    }
    
    w = moisth[idm].sorpiso.isotherm_value (in[0]); // volumetric moisture content from moisture storage function
    break;
    }
    default:{
    print_err("unknown type of material for homogenization",__FILE__,__LINE__,__func__);
    }
    }
    break;
    */
    //////////////////////////////////////////
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mattyp)
  
  return w;
}

/**
  Function returns relative humidity for individual material models.
  
  @param nid - node id
  @param mattype - type of material
  
  @return The function returns actual value of relative humidity.
  
  JK, 10. 4. 2014
*/
double transmat::give_nodal_rel_hum (long nid,long mattype)
{
  double rh;
  
  //  switch with respect to material type
  switch (mattype){
  case 154:{
    //  Milly material model
    rh = Tt->nodes[nid].eqother[2];
    break;
  }
  case 155:{
    //  Kunzel material model
    switch (int(Tt->nodes[nid].eqother[4])){
    case 8:{
      break;
    }
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
  case 156:
  case 157:{
    //  Grunewald and Devries material models
    rh = Tt->nodes[nid].eqother[0];
    break;
  }
  case 159:{
    //  simplediscmat material model
    rh = Tt->nodes[nid].eqother[1];
    break;
  }
   case 180:{
	//  MoistHeat model - modofication of Kunzel's model
   //	tep = Tt->nodes[nid].eqother[3];
	break;
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mattyp)
  
  return rh;
}


// **********************************************************************
// **********************************************************************
//
// Functions for retrieving of non-transport quantities from int. points
//
// **********************************************************************
// **********************************************************************



/**
  The function searches for the total number and type of non-transport quantities 
  that are required by used material models in the problem. Required non-transport quantities
  are searched for at all integration points in all used material models.
  Resulting number of quantities is returned and new array is allocated
  for the required types of quantities if they are required.

  @param - pointer to array of required types of non-transport quantities which is allocated inside 
           function (output parameter)
 
  @return Function returns the number of required non-transport quantities.

  Created by Tomas Krejci, 14/10/2013 according to Tomas Koudelka
*/
long transmat::search_reqntq(nontransquant* &rntq)
{
  long i, j;
  long antq[tnkntq];
  memset(antq, 0, sizeof(*antq)*tnkntq);

  // search and mark required non-transport quantities at all int. points
  for (i=0; i<tnip; i++) 
  {
    give_reqntq(i, antq);
  }
  
  if (Tp->advect==1){
    switch (Tp->gdim){
    case 1:{
      antq[advect_vel_x-1]=1;
      break;
    }
    case 2:{
      antq[advect_vel_x-1]=1;
      antq[advect_vel_y-1]=1;
      break;
    }
    case 3:{
      antq[advect_vel_x-1]=1;
      antq[advect_vel_y-1]=1;
      antq[advect_vel_z-1]=1;
      break;
    }
    default:{
      print_err("unknown geometrical dimension %ld is required",__FILE__,__LINE__,__func__,Tp->gdim);
      abort();
    }
    }

  }
  
  // compute the number of required quantities
  nntq = 0; // declared in transmat
  for (i=0; i<tnkntq; i++)
  {
    if (antq[i] == 1)
      nntq++;
  }

  // no non-mechanical quantities were required
  if (nntq == 0) 
  {
    rntq = NULL;
    return nntq;
  }
  else
  {
    // store types of required non-mechanical quantities
    rntq = new nontransquant[nntq];
    j = 0;
    for (i=0; i<tnkntq; i++)
    {
      if (antq[i] == 1){
        rntq[j] = nontransquant(i+1);
	j++;
      }
    }
  }

  return nntq;
}



/**
  The function marks required types of non-transport quantities by im-th
  material model in the given int. point ipp. Type of required quantities
  are marked by 1 in the array anmq whose length is equaled to the total number
  of known non-transport quantities.

  @param ipp - integration point id
  @param antq - array with flags for used material types
                antq[i] = 1 => qunatity type nontransquant(i+1) is required
                antq[i] = 0 => qunatity type nontransquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.

  Created by Tomas Krejci, 14/10/2013 according to Tomas Koudelka
*/
void transmat::give_reqntq(long ipp, long *antq)
{
  long i;
  
  i = ip[ipp].idm;

  switch(ip[ipp].tm)
  {
  case isotransmat:{
    itrm[i].give_reqntq(antq);
    break;
  }
  case damisotransmat:{
    damitrm[i].give_reqntq(antq);
    break;
  }
  case kunzel:{
	kun[i].give_reqntq(antq);
	break;
  }
  case moistheat:{
	moisth[i].give_reqntq(antq);
	break;
  }
  case homomat:{
    hommat[i].give_reqntq(antq);
    break;
  }

    /* case nlisotransmat:
       case discontisotrmat:
       case cernyconcrete:
       case bazantpedersen:
       case pedersen:
       case carb1mat:
       case sejtkr:
       case consolwf2:
       case consolawf2:
       case consolhawf3:
       case consolhwf2:
       case lincoupledmat:
       case radiationmater:
       case richardsmat:
       case cementhydrmat:
       case kunzel2:
       case grunewald:
       case simplediscmat:
       case devries:
       case milly:
       case salt3mat:
       case glasgow:
       case salt1mat:
       case salt2mat:
       case C60bazantB:
       case soilmat1:
       break;
    */

  case interfacem:
    ifacemat[i].give_reqntq(antq);
    break;
  case concreteB:
    concrete[i].give_reqntq(antq);
    break;
  case baroghelB:
    break;
  case C60baroghelB:
    C60baroghel[i].give_reqntq(antq);
    break;
  case C30baroghelB:
    C30baroghel[i].give_reqntq(antq);
    break;
  case o30bazantB:
    o30bazant[i].give_reqntq(antq);
    break;
  case salt4mat:
    salt4[i].give_reqntq(antq);
    break;
  case consolawf1:
    consol_awf1[i].give_reqntq(antq);
    break;
  case consolwf1:
    consol_wf1[i].give_reqntq(antq);
    break;
  case consolwf2:
    consol_wf2[i].give_reqntq(antq);
    break;
  case consolawf2:
    consol_awf2[i].give_reqntq(antq);
    break;
  case consolhawf3:
    consol_hawf3[i].give_reqntq(antq);
    break;
  case consolhwf2:
    consol_hwf2[i].give_reqntq(antq);
    break;
    
  default:
    print_err("unknown material type is required", __FILE__, __LINE__, __func__);
  }
  
  return;
}


/**
  The function returns required non-transport quantity stored in the nontransq array.
  It is used either in case of prescribed temperature changes in the input file or
  in case of coupled problems for retrieving of values passed between MEFEL and TRFEL.

  @param qt - type of required quantity
  @param ipp - integration point id

  @return The function returns required non-transport quantity.

  Created by Tomas Koudelka, 7.10.2013
*/
double transmat::givenontransq(nontransquant qt, long ipp)
{
  long id = ntqid[qt-1];
  if (id < 0)
  {
    print_err("Required quantity %s is not used in the problem solved", 
              __FILE__, __LINE__, __func__, nontransquantstr[qt-1].alias);
    abort();
  }
  return nontransq[id][ipp];
}



/**
  The function stores required non-transport quantity into the nontransq array.
  It is used either in case of prescribed temperature changes in the input file or
  in case of coupled problems for retrieving of values passed between MEFEL and TRFEL.

  @param qt - type of required quantity
  @param ipp - integration point id
  @param val - stored value

  @return The function stores required non-transport quantity.

  Created by Tomas Koudelka, 7.10.2013
*/
void transmat::storenontransq(nontransquant qt, long ipp, double val)
{
  long id = ntqid[qt-1];
  if (id < 0)
  {
    print_err("Required quantity %s is not used in the problem solved", 
              __FILE__, __LINE__, __func__, nontransquantstr[qt-1].alias);
    abort();
  }
  nontransq[id][ipp] = val;
}



/**
  The function returns status of required non-transport quantity, i.e
  whether the given quantity is defined in the problem or not.

  @param qt - type of required quantity
  
  @retval 0 - the quantity is not defined in the problem solved
  @retval 1 - the quantity is defined in the problem solved

  Created by Tomas Koudelka, 7.10.2013
*/
long transmat::givestatusntq (nontransquant qt)
{
  if (ntqid[qt-1] < 0) // index was not set -> quantity is not defined
    return 0;

  return 1; // quantity is defined in all other cases
}



/**
  The function returns index of required non-transport quantity, i.e
  the firs index of nontransq array.

  @param qt - type of required quantity
  
  @retval -1 - the quantity is not defined in the problem solved
  @retval >-1 - index of required the quantity

  Created by Tomas Koudelka, 7.10.2013
*/
long transmat::givenontransqid (nontransquant qt)
{
  return ntqid[qt-1];
}



/**
   function computes new transmission coefficient for transmission on the boundary
   (third kind of boundary condition)

   @param new_trc - new transmission coefficient
   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nid     - node id
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
void transmat::transmission_transcoeff(double &new_trc,double trc,long ri,long ci,long nid,long bc,long ipp)
{
  new_trc = 0.0;
      
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    new_trc = m1.transmission_transcoeff(trc,ri,ci,nid,bc,ipp,0);//for testing
    break;    
  }
  case twomediacoup:{
    med2 m2;

    new_trc = m2.transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;    
  }
  case threemediacoup:{
    med3 m3;

    new_trc = m3.transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    new_trc = m4.transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }  
}


/**
   function computes new transmission coefficient for transmission on the boundary
   (third kind of boundary condition)

   @param new_trc - new transmission coefficient
   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nid     - node id
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
   @param flag    - coefficient is computing for what 0=matrix,1=loading vector
*/
void transmat::transmission_transcoeff (double &new_trc,double trc,long ri,long ci,long nid,long bc,long ipp,int flag)
{
  new_trc = 0.0;
      
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    new_trc = m1.transmission_transcoeff(trc,ri,ci,nid,bc,ipp,flag);//for testing
    break;    
  }
  case twomediacoup:{
    med2 m2;

    new_trc = m2.transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;    
  }
  case threemediacoup:{
    med3 m3;

    new_trc = m3.transmission_transcoeff(trc,ri,ci,nid,bc,ipp,flag);
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    new_trc = m4.transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }  
}



/**
   function computes new nodal value for transmission on the boundary
   (third kind of boundary condition)

   @param new_nodval - transformed nodal value on boundary
   @param nodval     - nodal value defined on boundary
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if it is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nid        - node id
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
void transmat::transmission_nodval (double &new_nodval,double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp)
{
  new_nodval = 0.0;
      
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    new_nodval = m1.transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;    
  }
  case twomediacoup:{
    med2 m2;

    new_nodval = m2.transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;    
  }
  case threemediacoup:{
    med3 m3;

    new_nodval = m3.transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    new_nodval = m4.transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }  
}



/**
   function computes flux on the boundary for transmission on the boundary
   (third kind of boundary condition)

   @param flux    - transformed nodal value on the boundary
   @param nodval  - prescribed nodal value on boundary
   @param trc2    - second prescribed transmission coefficient on the boundary, 
                    if it is needed (for example heat radiation coef.)
   @param ri      - row index
   @param ci      - column index
   @param nn      - node id
   @param bc      - type of boundary conditions
   @param ipp     - number of first integration point on element
*/
void transmat::transmission_flux(double &flux,double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  flux = 0.0;
      
  switch (Tp->tmatt){//  transported matter
  case onemedium:{
    med1 m1;

    flux = m1.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);//for testing
    break;    
  }
  case twomediacoup:{
    med2 m2;

    flux = m2.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;    
  }
  case threemediacoup:{
    med3 m3;

    flux = m3.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;    
  }
  case fourmediacoup:{
    med4 m4;

    flux = m4.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;    
  }
  default:{
    print_err("unknown number of transported media is required",__FILE__,__LINE__,__func__);
  }
  }  
}


/**
   function corrects values with respect to material model for PUC
   
   08/08/2017 by TKr
*/
void transmat::values_correction_puc (void)
{
  long i,j,ipp,nip;
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      ipp=Tt->elements[i].ipp[0][0];
      nip=Tt->give_tnip (i);
      for (j=0;j<nip;j++){
	values_correction_ipp(ipp);
	ipp++;
      }
    }
  }
}


/**
   function corrects values with respect to material model on integration points
   
   @param ipp - appropriate integration point number
   
   08/08/2017 by TKr
*/
void transmat::values_correction_ipp (long ipp)
{
  long i;
  
  i = ip[ipp].idm;

  switch (ip[ipp].tm){
  case isotransmat:{
    break;
  }
  case moistheat:{
    moisth[i].values_correction_ipp (ipp);
    break;
  }
  default:{
    print_err("unknown matrial type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function corrects values with respect to material model
   
   @param nv - %vector of nodal values
   @param ipp - appropriate integration point number
   
   JK, 14.7.2005
*/
void transmat::values_correction (vector &nv,long ipp)
{
  long i,idm;
  
  i = ip[ipp].idm;

  switch (ip[ipp].tm){
  case isotransmat:{
    break;
  }
  case tdisotransmat:{
    break;
  }
  case nlisotransmat:{
    break;}
  case homomat:{
    ipp=Tt->elements[0].ipp[0][0];
    //material number
    i = ip[ipp].idm;
    
    switch (hommat[i].hom_mattype){
    case 100:{
      break;
    }
    case 180:{
      if(Tm->moisth == NULL){
	print_err("Basic material for homogenization is not defined in macro problem",__FILE__,__LINE__,__func__);
	abort();
      }
      
      idm = hommat[i].hom_mattype_number;
      moisth[idm].values_correction (nv);
      break;
    }
    default:{
      print_err("unknown type of material for homogenization",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  case damisotransmat:
  case discontisotrmat:
  case interfacem:
  case cernyconcrete:
  case pedersen:
  case carb1mat:
  case cementhydrmat:
  case lincoupledmat:
  case richardsmat:
  case radiationmater:{
    break;
  }
  case sejtkr:{
    sejtkrm[i].values_correction (nv, ipp);
    break;
  }
  case consolawf1:{
    consol_awf1[i].values_correction (nv, ipp);
    break;
  }
  case consolwf1:{
    consol_wf1[i].values_correction (nv, ipp);
    break;
  }
  case bazantpedersen:{
    bazped[i].values_correction (nv);
    break;
  }
  case kunzel:{
    kun[i].values_correction (nv);
    break;
  }
  case kunzel2:{
    kun2[i].values_correction (nv, ipp);
    break;
  }
  case grunewald:{
    //grunw[i].values_correction (nv);
    break;
  }
  case simplediscmat:{
    //sdmat[i].values_correction (nv);
    break;
  }
  case devries:{
    //grunw[i].values_correction (nv);
    break;
  }
  case milly:{
    //grunw[i].values_correction (nv);
    break;
  }
  case moistheat:{
    moisth[i].values_correction (nv);
    break;
  }
  case glasgow:{
    break;
  }
  case salt1mat:{
    break;
  }
  case salt2mat:{
    salt2[i].values_correction (nv);
    break;
  }
  case salt3mat:{
    salt3[i].values_correction (nv,ipp);
    break;
  }
  case salt4mat:{
    salt4[i].values_correction (nv,ipp);
    break;
  }
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    mtph.values_correction (nv);
    break;
  }

  case soilmat1:{
    gmultiph gmtph;
    gmtph.values_correction (nv);
    break;
  }
  case consolwf2:{
    consol_wf2[i].values_correction (nv, ipp);
    break;
  }
  case consolawf2:{
    consol_awf2[i].values_correction (nv, ipp);
    break;
  }
  case consolhawf3:{
    consol_hawf3[i].values_correction (nv, ipp);
    break;
  }
  case consolhwf2:{
    consol_hwf2[i].values_correction (nv, ipp);
    break;
  }
  default:{
    print_err("unknown matrial type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function transforms nodal values with respect to material models
   function is used in problems with discontinuities
   
   example: let moisture be assumed as nodal value, moisture along
   material interface contains discontinuity because only relative
   humidity is continuous, this function transforms moisture obtained
   for one material model to new value from second material model
   
   @param nid - node id
   @param mti - type of material model on input
   @param idmi - material id on input
   @param in - input values
   @param mto - type of material model on output
   @param idmo - material id on output
   @param out - recalculated values
   
   7.6.2006, JK
*/
/*
void transmat::values_transformation_old (long nid,long jnid,mattypet mti,long idmi,double *inav,double *inpv,double *ineq,mattypet mto,long idmo,double *out)
{
  double av;
  
  switch (mti){
  case grunewald:{
    av = grunw[idmi].get_rel_hum (inav[0]);
    break;
  }
  case salt4mat:{
    //long i;
    //double *prevval;
     // ndofn = Tt->give_ndofn (nid);
    //prevval = new double [ndofn+2];
    //for (i=0;i<ndofn+2;i++){
    //prevval[i]=Tt->prevval[jnid][i];
    //}

    //av = salt4[idmi].get_rel_hum (inav[0]);
//    salt4[idmi].water_content_relhum (nid,inav,inpv,ineq,out);
   	av = out[10];
    //for (i=0;i<ndofn+2;i++){
    //Tt->prevval[jnid][i]=prevval[i];
    //}
    
    //delete [] prevval;
    break;
  } 
  case devries:{
    av = grunw[idmi].get_rel_hum (inav[0]);
    break;
  }
  case kunzel:{
    av = inav[0];
    break;
  }
  case discontisotrmat:{
    av=ditrm[idmi].compute_rel (inav[0]);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function values_transformation (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  

  switch (mto){
  case grunewald:{
    out[0] = grunw[idmo].get_moisture (av)-inav[0];
    out[1] = 0.0;
    break;
  }
  case salt4mat:{
	  double rhw, wrh;
	  salt4[idmo].get_moisture (nid,av,ineq,wrh);
	  out[11] = av;
	  rhw = wrh;
	  out[12] = rhw;	  ;
	  out[0] = rhw - inav[0];
	  out[1] = 0.0;
	  out[2] = 0.0;
	  out[3] = 0.0;
	  //	fprintf (stdout,"\n Time  %lf Node id %d  vlhkost in %lf relhum av %lf vlhkost out %lf  skok out %lf\n",Tp->time, nid, in[0], av,rhw, out[0]);
	  

    break;
  } 
  case devries:{
    break;
  }

  case kunzel:{
	  out[0] = 0.0;
	  out[1] = 0.0;
    break;
  }
  case discontisotrmat:{
    //ditrm[idmo].correct_val (in,out);
    out[0]=ditrm[idmo].compute_abs (av);
    break;
  }
  case isotransmat:{
    out[0]=-5.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function values_transformation (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  
  //fprintf (Outt,"   in  %lf   av %lf  out %lf\n",in[0],av,out[0]);
}
*/


/**
   function transforms nodal values with respect to material models
   function is used in problems with discontinuities
   
   example: let moisture be assumed as nodal value, moisture along
   material interface contains discontinuity because only relative
   humidity is continuous, this function transforms moisture obtained
   for one material model to new value from second material model
   
   @param mti - type of material model on input
   @param idmi - material id on input
   @param inv - input values
   @param iinv - initial input values
   @param mto - type of material model on output
   @param idmo - material id on output
   @param outv - recalculated values
   @param ioutv - initial values on output
   @param jum - jump on material interface
   
   JK, 8.8.2008
*/
void transmat::values_transformation (mattypet mti,long idmi,double *inv,double *iinv,
				      mattypet mto,long idmo,double *outv,double *ioutv,double *jum)
{
  double av,avi, fii, fio, dfdwi, dfdwo;
  
  switch (mti){
  case milly:
  case moistheat:
  case kunzel:{
    // av = grunw[idmi].get_rel_hum (inv[0]);
    // avi = grunw[idmi].get_rel_hum (iinv[0]);
    break;
  }
  case grunewald:{
    grunw[idmi].get_rel_hum (inv[0], fii, dfdwi);
    grunw[idmo].get_rel_hum (outv[0], fio, dfdwo);
    if (dfdwi > dfdwo)
      {
	av = fio;
      }
    else
      {
	av = fii;
      }
    //av=grunw[idmi].get_rel_hum (inv[0]);
    break;
  }
  case simplediscmat:{
    
    sdmat[idmi].get_rel_hum2 (inv[0], fii, dfdwi);
    sdmat[idmo].get_rel_hum2 (outv[0], fio, dfdwo);
    if (dfdwi > dfdwo)
      {
	av = fio;
      }
    else
      {
	av = fii;
      }
    
    break;
  }
  case discontisotrmat:{
    av =ditrm[idmi].correct_val (inv,iinv);
    break;
  }
  case salt3mat:{
    av = salt3[idmi].get_rel_hum (inv[0]);
    avi = salt3[idmi].get_rel_hum (iinv[0]);
    break;
  }	  
  case salt4mat:{
	/*av = salt4[idmi].get_rel_hum2 (inv[0]);
	avi = salt4[idmi].get_rel_hum2 (iinv[0]);

    double dav, davi;
	salt4[idmi].get_rel_hum2 (inv[0], av, dav);
	salt4[idmi].get_rel_hum2 (iinv[0], avi, davi);
	*/
    break;
  } 

  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  switch (mto){
  case milly:
  case moistheat:
  case kunzel:{
    //jum[0] = grunw[idmo].get_moisture (av)-inv[0];
    //av = grunw[idmo].get_moisture (avi)-iinv[0];
    //jum[0]-=av;
    //jum[0]=0.001;
    jum[0] = 0.0;
    jum[1] = 0.0;
    break;
  }
  case grunewald:{
    if (dfdwi > dfdwo)
      {
	jum[0] = -grunw[idmi].get_moisture (av)+outv[0]+iinv[0]-ioutv[0];
      }
    else
      {
	jum[0] = grunw[idmo].get_moisture (av)-inv[0]-ioutv[0]+iinv[0];
      }
    
    //jum[0] = grunw[idmo].get_moisture (av)-inv[0]-ioutv[0]+iinv[0];
    jum[1] = 0.0;
    break;
  }
  case simplediscmat:{
    if (dfdwi > dfdwo)
      {
	jum[0] = -sdmat[idmi].get_moisture2 (av)+outv[0]+iinv[0]-ioutv[0];
      }
    else
      {
	jum[0] = sdmat[idmo].get_moisture2 (av)-inv[0]-ioutv[0]+iinv[0];
      }
    
    jum[1] = 0.0;
    break;
  }
    
  case discontisotrmat:{
    avi=ditrm[idmo].correct_val (outv,ioutv);
    //jum[0]=avi-av;
    jum[0]=avi-ioutv[0]+iinv[0];
    break;
  }
  case salt3mat:{
    jum[0] = salt3[idmo].get_moisture (av)-inv[0];
    av = salt3[idmo].get_moisture (avi)-iinv[0];
    jum[0]-=av;
    jum[1] = 0.0;
    jum[2] = 0.0;
    break;
  }
  case salt4mat:{
    //jum[0] = salt4[idmo].get_moisture2 (av)-inv[0];
    //av = salt4[idmo].get_moisture2 (avi)-iinv[0];
    jum[0]-=av;
    jum[1] = 0.0;
    jum[2] = 0.0;
    jum[3] = 0.0;
    break;
  } 
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  //fprintf (Outt,"   jump %20.15lf \n",jum[0]);
}


/**
   function updates values at integration points
*/
void transmat::updateipval (void)
{
  long i,j,ipp,nip;
  
  //tady pridano 17.7.2006 TKr
  //updates only used elements(int. points)
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      ipp=Tt->elements[i].ipp[0][0];
      nip=Tt->give_tnip (i);
      for (j=0;j<nip;j++){
	updateipvalmat (ipp,0,0);
	ipp++;
      }
    }
  }
  
  //  loop over the number of elements
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      //  only switched on elements are taken into account
      
      //  id of integration point
      ipp=Tt->elements[i].ipp[0][0];
      //  the number of integration points on element
      nip=Tt->give_tnip (i);
      
      //  loop over the number of integration points on element
      for (j=0;j<nip;j++){
	//fluxcomparing (ipp);
	ipp++;
      }//  end of the loop over the number of integration points on element
    }//  end of the if (Gtt->leso[i]==1) statement
  }//  end of the loop over the number of elements
  
  //for (i=0;i<tnip;i++){
  //fprintf (stdout,"\n int point %6ld   %ld %ld %ld %ld",i,ip[i].infl[0][0],ip[i].infl[0][1],ip[i].infl[1][0],ip[i].infl[1][1]);
  //}
}



/**
  The function updates values at auxiliary integration points.

  @return The function does not return anything but it changes content of auxiliary int. points.

  Created by Tomas Koudelka, 25.5.2018
*/
void transmat::update_aipval (void)
{
  intpointst *tmp_ip;
  long *tmp_elip;
  double **tmp_ntq;
  double *tmp_iv;
  long app, eid;
  
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_ntq = Tm->nontransq;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->nontransq = Tm->aip_nontransq;
  Tm->initval = Tm->aip_initval;

  for (app=0; app<Tm->tnaip; app++)
  {
    eid = Tm->elip[app];
    //  only elements switched on are taken into account
    if (Gtt->leso[eid] == 1)
      updateipvalmat (app,0,0);
  }
  
  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->nontransq = tmp_ntq;
  Tm->initval = tmp_iv;
}



/**
   function updates values at integration points
   
   @param ipp - integration point pointer
   @param im  - index of material type for given ip
   @param ido - index in array eq_other   

   12.12.2006, TKr
*/
void transmat::updateipvalmat (long ipp,long im,long ido)
{
  long i;
  
  i = ip[ipp].idm;
  
  switch (ip[ipp].tm){
  case isotransmat:
  case nlisotransmat:
  case homomat:
  case damisotransmat:
  case discontisotrmat:
  case tdisotransmat:
  case interfacem:
  case cernyconcrete:
  case bazantpedersen:
  case pedersen:
  case carb1mat:
  case sejtkr:
  case lincoupledmat:
  case richardsmat:
  case radiationmater:
    break;
  case consolawf1:
    consol_awf1[i].updateval(ipp);
    break;
  case consolwf1:
    consol_wf1[i].updateval(ipp);
    break;
  case consolwf2:
    consol_wf2[i].updateval(ipp);
    break;
  case consolawf2:
    consol_awf2[i].updateval(ipp);
    break;
  case consolhawf3:
    consol_hawf3[i].updateval(ipp);
    break;
  case consolhwf2:
    consol_hwf2[i].updateval(ipp);
    break;
  case cementhydrmat:
    cemhydr[i].updateval(ipp,im,ido);
    break;
  case kunzel:{
    //kun[i].updateval(ipp,im,ido);
    break;
  }
  case kunzel2:{
    //kun2[i].updateval(ipp,im,ido);
    break;
  }
  case grunewald:{
    //grunw[i].updateval(ipp,im,ido);
    break;
  }
  case simplediscmat:{
    //sdmat[i].updateval(ipp,im,ido);
    break;
  }
  case devries:{
    //grunw[i].updateval(ipp,im,ido);
    break;
  }
  case milly:{
	//grunw[i].updateval(ipp,im,ido);
	break;
  }
  case moistheat:{
	//grunw[i].updateval(ipp,im,ido);
	break;
  }
  case glasgow:
  case salt1mat:
  case salt2mat:
  case salt3mat:
  case salt4mat:
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:
  case soilmat1:
    break;
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function computes and compares contributions to the fluxes
   array infl defined in integration points is assembled
   
   @param ipp - integration point pointer
   
   JK, 8.12.2011
*/
void transmat::fluxcomparing (long ipp)
{
  long i,j,n,ntm;
  double norfl1,norfl2,neglrel,neglabs,zero;
  
  //  relative value for flux neglecting
  neglrel=1.0e-5;
  //  absolute value for flux neglecting
  neglabs=1.0e-9;
  //  computer zero
  zero=Tp->zero;
  
  //  the number of transported media
  ntm = Tp->ntm;
  //  number of gradient and flux components
  n = ip[ipp].ncompgrad;
  
  vector gr(ASTCKVEC(n)),fl1(ASTCKVEC(n)),fl2(ASTCKVEC(n));
  matrix d(ASTCKMAT(n,n));
  
  for (i=0;i<ntm;i++){
    //  conductivity matrix of material
    matcond (d,ipp,i,i);
    //  gradient 
    givegrad (i,ipp,gr);
    //  flux
    mxv (d,gr,fl1);
    //  norm of the flux
    norfl1 = normv (fl1);
    //  influence indicator
    ip[ipp].infl[i][i]=1;
    
    for (j=0;j<ntm;j++){
      if (j==i)
	continue;
      
      //  conductivity matrix of material
      matcond (d,ipp,i,j);
      //  gradient 
      givegrad (j,ipp,gr);
      //  flux
      mxv (d,gr,fl2);
      //  norm of the flux
      norfl2 = normv (fl2);
      
      if (norfl1>zero){
	if (norfl2/norfl1 < neglrel){
	  //  off-diagonal flux is neglected
	  //  influence indicator
	  //fprintf (stdout,"\n\n vypnuto v bode %ld \n",ipp);
	  ip[ipp].infl[i][j]=0;
	}else{
	  //  off-diagonal flux is not neglected
	  //  influence indicator
	  ip[ipp].infl[i][j]=1;
	}
      }else{
	if (norfl2 < neglabs){
	  //  off-diagonal flux is neglected
	  //  influence indicator
	  //fprintf (stdout,"\n\n vypnuto v bode %ld \n",ipp);
	  ip[ipp].infl[i][j]=0;
	}else{
	  //  off-diagonal flux is not neglected
	  //  influence indicator
	  ip[ipp].infl[i][j]=1;
	}
      }
      
      /*
      if (norfl2>zero){
	if (norfl1/norfl2 < neglrel){
	  //  diagonal flux is neglected
	  //  influence indicator
	  //fprintf (stdout,"\n\n vypnuto v bode %ld \n",ipp);
	  ip[ipp].infl[i][i]=0;
	}else{
	  //  diagonal flux is not neglected
	  //  influence indicator
	  ip[ipp].infl[i][i]=1;
	}
      }else{
	if (norfl1 < neglabs){
	  //  diagonal flux is neglected
	  //  influence indicator
	  //fprintf (stdout,"\n\n vypnuto v bode %ld \n",ipp);
	  ip[ipp].infl[i][i]=0;
	}else{
	  //  diagonal flux is not neglected
	  //  influence indicator
	  ip[ipp].infl[i][i]=1;
	}
      }
      */
    }
  }
}


/**
   function updates values at integration points
*/
void transmat::initmaterialmodels (void)
{
  long i,j,ipp,nip;
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      ipp=Tt->elements[i].ipp[0][0];
      nip=Tt->give_tnip (i);
      for (j=0;j<nip;j++){
	initvalues (ipp,0,0);
	ipp++;
      }
    }
  }
}



/**
  Function initializes material models with initial values on auxiliary integration points.
   The function is intended for the transfer of values among meshes in coupled problems especially.

  @param lcid - load case id
  @param n - the number of required auxiliary points in the mapping array ipm
  @param ipm - integration point mapping array, 
               ipm[i].ipp < 0 => auxiliary integration point must be used
               ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                  no computation of strains is performed, the strains are assumed 
                                  to be computed at the main solution procedure of the problem

  @return The function does not return anything but it chages state of auxiliary integration points
          at the array aip.

  Created by Tomas Koudelka, 28.11.2017
*/
void transmat::aip_initmaterialmodels (void)
{
  long app;
  intpointst *tmp_ip;
  long *tmp_elip;
  double **tmp_ntq;
  double *tmp_iv;

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = ip;
  tmp_elip = elip;
  tmp_ntq = nontransq;
  tmp_iv = initval;
  ip = aip;
  elip = elaip;
  nontransq = aip_nontransq;
  initval = aip_initval;

  for (app=0; app<tnaip; app++)
  {
    if (Gtt->leso[elip[app]] == 1)
      initvalues(app, 0, 0);
  }

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  ip = tmp_ip;
  elip = tmp_elip;
  nontransq = tmp_ntq;
  initval = tmp_iv;
}



/**
   function updates values at integration points
   
   @param ipp - integration point pointer
   @param im  - index of material type for given ip
   @param ido - index in array eq_other   

   27.3.2007, TKr
*/
void transmat::initvalues (long ipp,long im,long ido)
{
  long i;
  
  i = ip[ipp].idm;
  
  switch (ip[ipp].tm){
  case isotransmat:
  case nlisotransmat:
  case homomat:
  case damisotransmat:
  case discontisotrmat:
  case tdisotransmat:
  case interfacem:
  case cernyconcrete:
  case bazantpedersen:
  case pedersen:
  case carb1mat:
  case sejtkr:
  case richardsmat:
  case lincoupledmat:
    break;
  case consolawf1:
    consol_awf1[i].initval(ipp);
    break;
  case consolwf1:
    consol_wf1[i].initval(ipp);
    break;
  case consolwf2:
    consol_wf2[i].initval(ipp);
    break;
  case consolawf2:
    consol_awf2[i].initval(ipp);
    break;
  case consolhawf3:
    consol_hawf3[i].initval(ipp);
    break;
  case consolhwf2:
    consol_hwf2[i].initval(ipp);
    break;
  case cementhydrmat:
    cemhydr[i].initvalues(ipp,im,ido);
    break;
  case kunzel:{
    kun[i].initvalues (ipp,ido);
    break;
  }
  case kunzel2:{
    kun2[i].initvalues (ipp,ido);
    break;
  }
  case grunewald:{
    grunw[i].initvalues (ipp,ido);
    break;
  }
  case simplediscmat:{
    sdmat[i].initvalues (ipp,ido);
    break;
  }
  case devries:{
    //grunw[i].updateval(ipp,im,ido);
    break;
  }
  case milly:{
    mill[i].initvalues(ipp,ido);
    break;
  }
  case moistheat:{
	moisth[i].initvalues(ipp,ido);
	break;
  }
  case salt3mat:{
    salt3[i].initvalues (ipp,ido);
    break;
  }
  case salt4mat:{
    //salt4[i].initvalues (ipp,ido);
    break;
  }
  case glasgow:
  case salt1mat:
  case salt2mat:
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:
  case soilmat1:
    break;
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes auxiliary parameters
   
   @param lcid - load case id
   
   JK
*/
void transmat::aux_values (long lcid)
{
  long i,j,ii,idm,nid,nne,ncompeqother,eid,ipp;
  double *in,*out, *inp, *ineq;
  ivector nodes;
  
  //  loop over the number of integration points
  for (i=0;i<tnip;i++){
    
    j = ip[i].idm;
    
    switch (ip[i].tm){
    case isotransmat:
    case nlisotransmat:
    case homomat:
    case damisotransmat:
    case discontisotrmat:
    case tdisotransmat:
    case interfacem:
    case cernyconcrete:
    case bazantpedersen:
    case pedersen:
    case carb1mat:
    case cementhydrmat:
    case richardsmat:
    case radiationmater:{
      break;
    }
    case kunzel:{
      ncompeqother = givencompeqother (i,0);
      in = new double[2];
      inp = new double[2];
      ineq = new double[ncompeqother];
      out = new double[ncompeqother];

      //initialization
      for(ii=0;ii<2;ii++){
	inp[ii] = 0.0;
	in[ii] = 0.0;
      }
      for(ii=0;ii<ncompeqother;ii++){
	ineq[ii] = 0.0;
	out[ii] = 0.0;
      }
      
      kun[j].give_values (i,in,inp,ineq);
      kun[j].aux_values (i,in,inp,ineq,out);
      kun[j].save_values (i,out);

      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case sejtkr:
    case consolawf1:
    case consolwf1:
    case lincoupledmat:
      break;
    case kunzel2:{
      in = new double[2];
      out = new double[6];
      inp = new double[2];
      ineq = new double[6];
      kun2[j].give_values (i,in, inp,ineq);
      kun2[j].aux_values (i,in, inp,ineq,out);
      kun2[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case grunewald:{
      in = new double[2];
      out = new double[6];
      inp = new double[2];
      ineq = new double[6];
      grunw[j].give_values (i,in, inp,ineq);
      grunw[j].aux_values (i,in, inp,ineq,out);
      grunw[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    
    case simplediscmat:{
      in = new double[4];
      inp = new double[4];
      out = new double[8];
      ineq = new double[8];
      sdmat[j].give_values (i,in, inp,ineq);
      sdmat[j].aux_values (in, inp,ineq,out);
      sdmat[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }

    case devries:{
      break;
    }
    case milly:{
      in = new double[2];
      out = new double[10];
      inp = new double[2];
      ineq = new double[10];
      mill[j].give_values (i,in, inp,ineq);
      mill[j].aux_values (i,in, inp,ineq,out);
      mill[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case moistheat:{
      ncompeqother = givencompeqother (i,0);
      in = new double[2];
      inp = new double[2];
      ineq = new double[ncompeqother];
      out = new double[ncompeqother];
      moisth[j].give_values (i,in,inp,ineq);
      moisth[j].aux_values (i,in,inp,ineq,out);
      moisth[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case salt3mat:{
      salt3[j].aux_values (i);
      break;
    }
    case salt4mat:{
      ncompeqother = givencompeqother (i,0);
      in = new double[4];
      inp = new double[4];
      ineq = new double[ncompeqother];
      out = new double[ncompeqother];
      salt4[j].give_values (i,in,inp,ineq);
      salt4[j].aux_values (i,in,inp,ineq,out);
      salt4[j].save_values (i,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case glasgow:
    case salt1mat:
    case salt2mat:
    case concreteB:
    case baroghelB:
    case C60baroghelB:
    case C30baroghelB:
    case o30bazantB:
    case C30bazantB:
    case C60bazantB:
    case soilmat1:
    case consolwf2:
    case consolawf2:
    case consolhawf3:
    case consolhwf2:
      break;
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
    }
    
  }//  end of the statement for (i=0;i<tnip;i++){
  
  
  // ******************************************************************
  //  computation of values necessary for climatic boundary conditions
  // ******************************************************************
  if (Tb->lc[lcid].ncc>0){
    //  loop over the number of elements with prescribed climatic boundary conditions
    for (i=0;i<Tb->lc[lcid].neb;i++){
      
      //  element id
      eid = Tb->lc[lcid].elemload[i].eid;
      
      if (eid >=0){
	//  integration point id
	ipp = Tt->elements[eid].ipp[0][0];
	//  number of nodes on element
	nne = Tt->give_nne (eid);
	//  allocation of array for nodes
	allocv (nne,nodes);
	//  nodes on element
	Tt->give_elemnodes (eid,nodes);
	
	//  id of material
	idm = ip[ipp].idm;
	
	//  loop over the number of nodes on element
	for (j=0;j<nne;j++){
	  //  node id
	  nid = nodes[j];
	  
	  switch (ip[ipp].tm){
	  case isotransmat:
	  case tdisotransmat:
	  case homomat:
	    //here must be completeted
	  case damisotransmat:
	    break;
	  case kunzel:{
	    in = new double[2];
	    inp = new double[2];
	    ineq = new double[5];
	    out = new double[5];
	    Tt->nodes[nid].give_values (in,inp,ineq);
	    kun[idm].aux_values (ipp,in,inp,ineq,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }
	  case kunzel2:{
	    in = new double[2];
	    inp = new double[2];
	    ineq = new double[15];
	    out = new double[15];
	    Tt->nodes[nid].give_values (in, inp,ineq);
	    kun2[idm].aux_values (ipp,in, inp,ineq,out);
	    //	  Tt->nodes[nid].give_values (in);
	    //	  kun2[idm].aux_values (ipp,in,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }
	  case grunewald:{
	    in = new double[2];
	    inp = new double[2];
	    ineq = new double[15];
	    out = new double[15];
	    Tt->nodes[nid].give_values (in, inp,ineq);
	    grunw[idm].aux_values (ipp,in, inp,ineq,out);
	    //	  Tt->nodes[nid].give_values (in);
	    //	  kun2[idm].aux_values (ipp,in,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }  
	  case simplediscmat:{
	    in = new double[2];
	    inp = new double[2];
	    ineq = new double[15];
	    out = new double[15];
	    Tt->nodes[nid].give_values (in, inp,ineq);
	    sdmat[idm].aux_values (in, inp,ineq,out);
	    //	  Tt->nodes[nid].give_values (in);
	    //	  kun2[idm].aux_values (ipp,in,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }  
	  case salt4mat:{
	    ncompeqother = givencompeqother (i,0);
	    in = new double[4];
	    inp = new double[4];
	    ineq = new double[ncompeqother];
	    out = new double[ncompeqother];
	    Tt->nodes[nid].give_values (in,inp,ineq);
	    salt4[idm].aux_values (ipp,in,inp,ineq,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }
	  case milly:{
	    in = new double[2];
	    inp = new double[2];
	    ineq = new double[15];
	    out = new double[15];
	    Tt->nodes[nid].give_values (in, inp,ineq);
	    mill[idm].aux_values (ipp,in, inp,ineq,out);
	    //	  Tt->nodes[nid].give_values (in);
	    //	  kun2[idm].aux_values (ipp,in,out);
	    Tt->nodes[nid].save_values (out);
	    delete [] in;
	    delete [] inp;
	    delete [] ineq;
	    delete [] out;
	    break;
	  }
	  default:{
	    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
	  }
	  }//  end of switch (ip[ipp].tm){
	}//  end of he loop for (j=0;j<nne;j++){
	destrv(nodes);
      }
    }//  end of the loop for (i=0;i<Tb->lc[lcid].neb;i++){
  }//  end of the statement if (Tb->lc[lcid].ncc>0){
}



/**
   function computes auxiliary parameters
   
   @param ipp - integration point id
   
   JK
*/
void transmat::mat_aux_values (long ipp)
{
  long im, ncompeqother;
  double *in,*out, *inp, *ineq;
  
  im = ip[ipp].idm;
    
  switch (ip[ipp].tm){
    case isotransmat:
    case nlisotransmat:
    case homomat:
    case damisotransmat:
    case discontisotrmat:
    case tdisotransmat:
    case interfacem:
    case cernyconcrete:
    case bazantpedersen:
    case pedersen:
    case carb1mat:
    case cementhydrmat:
    case richardsmat:
    case radiationmater:{
      break;
    }
    case kunzel:{
      ncompeqother = givencompeqother (ipp,0);
      in = new double[2];
      inp = new double[2];
      ineq = new double[ncompeqother];
      out = new double[ncompeqother];
      kun[im].give_values (ipp,in,inp,ineq);
      kun[im].aux_values (ipp,in,inp,ineq,out);
      kun[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case sejtkr:
    case consolawf1:
    case consolwf1:
    case lincoupledmat:
      break;
    case kunzel2:{
      in = new double[2];
      out = new double[6];
      inp = new double[2];
      ineq = new double[6];
      kun2[im].give_values (ipp,in, inp,ineq);
      kun2[im].aux_values (ipp,in, inp,ineq,out);
      kun2[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case grunewald:{
      in = new double[2];
      out = new double[6];
      inp = new double[2];
      ineq = new double[6];
      grunw[im].give_values (ipp,in, inp,ineq);
      grunw[im].aux_values (ipp,in, inp,ineq,out);
      grunw[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    
    case simplediscmat:{
      in = new double[4];
      inp = new double[4];
      out = new double[8];
      ineq = new double[8];
      sdmat[im].give_values (ipp,in, inp,ineq);
      sdmat[im].aux_values (in, inp,ineq,out);
      sdmat[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }

    case devries:{
      break;
    }
    case milly:{
      in = new double[2];
      out = new double[10];
      inp = new double[2];
      ineq = new double[10];
      mill[im].give_values (ipp,in, inp,ineq);
      mill[im].aux_values (ipp,in, inp,ineq,out);
      mill[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case salt3mat:{
      salt3[im].aux_values (ipp);
      break;
    }
    case salt4mat:{
      ncompeqother = givencompeqother (ipp,0);
      in = new double[4];
      inp = new double[4];
      ineq = new double[ncompeqother];
      out = new double[ncompeqother];
      salt4[im].give_values (ipp,in,inp,ineq);
      salt4[im].aux_values (ipp,in,inp,ineq,out);
      salt4[im].save_values (ipp,out);
      delete [] in;
      delete [] inp;
      delete [] ineq;
      delete [] out;
      break;
    }
    case glasgow:
    case salt1mat:
    case salt2mat:
    case concreteB:
    case baroghelB:
    case C60baroghelB:
    case C30baroghelB:
    case o30bazantB:
    case C30bazantB:
    case C60bazantB:
    case soilmat1:
    case consolwf2:
    case consolawf2:
    case consolhawf3:
    case consolhwf2:
      break;
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
  }
}



/**
   function detects cycles
   it is intended for freezing cycles
   in the case of new cycle, function returns 1, otherwise, it returns 0
   
   @param r - nodal values
   @param pr - nodal values from the previous time step
   @param ipp - integration point id
   
   11.9.2007, JK
*/
long transmat::cycle_detection (double *r,double *pr,double *ppr,long ipp)
{
  long ret;
  
  switch (ip[ipp].tm){
  case salt4mat:{
    ret=salt4[ip[ipp].idm].cycle_detection (r,pr,ppr);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return ret;
}



//**************************************************
// Functions for backup of regular integration points
//***************************************************



/*
  Function saves data from regular integration points into backup file
  in text format.
   
  @param aux - pointer to auxiliary file
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  JK, 19.9.2004
  TKo, 9.2008
  TKr, 12.2008
*/
void transmat::save_intpointst_txt (FILE *aux, sel &selelems, sel *selother)
{
  long i,j,k,n,nip,ir,ipp;
  sel selno; // selection with no selected components
  int prec = (int)Tp->hdbcont.prec;
  
  selno.n=1;
  fprintf(aux, "\n");
  for (i=0;i<Tt->ne;i++)
  {    
    ipp = Tt->elements[i].ipp[0][0];
    nip = Tt->give_tnip(i);
    selelems.presence_id(i, ir);
    if (ir < 0)
      n = 0;    
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fprintf (aux,"%ld %ld %ld %ld\n",ipp, Tp->ntm, ip[ipp].ncompgrad, n);
      if (ir < 0)
        ip[ipp].save_data_txt(aux, selno);
      else
        ip[ipp].save_data_txt(aux, selother[ir]);
	  
      fprintf(aux, "\n");
      // writing MEFEL quantities
      if (nontransq)
      {
        for (k=0; k<nntq; k++)
          fprintf (aux,"%.*le ",prec, nontransq[k][ipp]);
        fprintf(aux, "\n");
      }
      ipp++;
    }
  }
}


/**
  Function saves data from regular integration points into several backup files
  in text format.
   
  @param ni       - time step id
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  JK, 19.9.2004
  TKo, 9.2008
  TKr, 12.2008
*/
void transmat::save_intpointst_txt (long ni, sel &selelems, sel *selother)
{
  long i,j,k,ir,n,nip,ipp;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Tp->hdbcont.prec;
  
  sprintf(name, "%s.%ld.grad.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fprintf(aux,"%ld %ld %ld\n",i, Tp->ntm, ip[i].ncompgrad);
    for (j=0;j<Tp->ntm;j++)
      for (k=0;k<ip[i].ncompgrad;k++)
        fprintf (aux, "%.*le\n",prec , ip[i].grad[j][k]);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.flux.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fprintf (aux,"%ld %ld %ld\n",i, Tp->ntm, ip[i].ncompgrad);
    for (j=0;j<Tp->ntm;j++)
      for (k=0;k<ip[i].ncompgrad;k++)
        fprintf(aux, "%.*le\n",prec , ip[i].fluxes[j][k]);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.other.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<Tt->ne;i++)
  {    
    ipp = Tt->elements[i].ipp[0][0];
    nip = Tt->give_tnip(i);
    selelems.presence_id(i, ir);    
    if (ir < 0)
      n = 0;
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fprintf (aux,"%ld %ld\n",ipp, n);
      if (n==0) 
      {
        ipp++;
        continue;
      }
      for (k=0;k<ip[ipp].ncompeqother;k++)
      {
        if (selother[ir].presence_id(k))
          fprintf (aux, "%.*le\n",prec , ip[ipp].eqother[k]);
      }
      ipp++;
    }
  }
  fclose(aux);

  // writing MEFEL quantities
  if (nontransq)
  {
    sprintf(name, "%s.%ld.mefquant.bac", Tp->hdbcont.hdbnames, ni);
    aux = fopen(name, "wt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    // write table of indices of non-mechanical quantities
    for (i=0; i<tnkntq; i++)
      fprintf(aux, "%ld ", ntqid[i]);
    fprintf(aux, "\n");

    for (i=0; i<tnip; i++)
    {
      fprintf (aux, "%ld ", i);
      for (j=0; j<nntq; j++)
        fprintf (aux,"%.*le ", prec, nontransq[j][i]);
        
      fprintf(aux, "\n");
    }
    fclose(aux);
  }
}



/**
   Function saves data from regular integration points into backup file
   in binary format.
   
   @param aux - pointer to auxiliary file
   @param selelems - selection of elements whose eqother array will be saved
   @param selother - selection of components of saved eqother array which will be saved on selected elements
                     for each range or list item in selelems an individual selection of eqother components
                     is enabled
   
   TKo, 9.2008
   TKr, 12.2008
*/
void transmat::save_intpointst_bin (FILE *aux, sel &selelems, sel *selother)
{
  long i,j,k,n,nip,ir,ipp;
  sel selno; // selection with no selected components
  
  selno.n=1;
  for (i=0;i<Tt->ne;i++)
  {    
    ipp = Tt->elements[i].ipp[0][0];
    nip = Tt->give_tnip(i);
    selelems.presence_id(i, ir);
      
    if (ir < 0)
      n = 0;
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fwrite(&ipp, sizeof(ipp), 1, aux);
      fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
      fwrite(&ip[ipp].ncompgrad, sizeof(ip[ipp].ncompgrad), 1, aux);
      fwrite(&n, sizeof(n), 1, aux);
      if (ir < 0)
        ip[ipp].save_data_bin(aux, selno);
      else
        ip[ipp].save_data_bin(aux, selother[ir]);
      
      // writing MEFFEL quantities
      if (nontransq)
      {
        for (k=0; k<nntq; k++)
          fwrite (&nontransq[k][ipp], sizeof(nontransq[k][ipp]), 1, aux);
      }
      ipp++;
    }
  }
}



/**
  Function saves data from regular integration points into several backup files
  in binary format.
   
  @param ni       - time step id
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  TKo 9.2008
  TKr, 12.2008
*/
void transmat::save_intpointst_bin (long ni, sel &selelems, sel *selother)
{
  long i, j, k, ir, n, nip, ipp;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  
  sprintf(name, "%s.%ld.grad.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fwrite(&i, sizeof(i), 1, aux);
    fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
    fwrite(&ip[i].ncompgrad, sizeof(ip[i].ncompgrad), 1, aux);
    for (j=0; j<Tp->ntm; j++)
      fwrite(ip[i].grad[j], sizeof(*ip[i].grad[j]), ip[i].ncompgrad, aux);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.flux.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fwrite(&i, sizeof(i), 1, aux);
    fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
    fwrite(&ip[i].ncompgrad, sizeof(ip[i].ncompgrad), 1, aux);
    for (j=0; j<Tp->ntm; j++)
      fwrite(ip[i].fluxes[j], sizeof(*ip[i].fluxes[j]), ip[i].ncompgrad, aux);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.other.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<Tt->ne;i++)
  {    
    ipp = Tt->elements[i].ipp[0][0];
    nip = Tt->give_tnip(i);
    selelems.presence_id(i, ir);    
    if (ir < 0)
      n = 0;
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fwrite(&ipp, sizeof(ipp), 1, aux);
      fwrite(&n, sizeof(n), 1, aux);
      if (n==0) 
      {
        ipp++;
        continue;
      }
      for (k=0;k<ip[ipp].ncompeqother;k++)
      {
        if (selother[ir].presence_id(k))
          fwrite(ip[ipp].eqother+k, sizeof(*ip[ipp].eqother), 1, aux);
      }
      ipp++;
    }
  }
  fclose(aux);

  if (nontransq)  
  {
    // writing MEFEL quantities
    sprintf(name, "%s.%ld.mefquant.bac",Tp->hdbcont.hdbnames, ni);
    aux = fopen(name,"wb");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    // write table of indices of non-mechanical quantities
    fwrite (ntqid, sizeof(*ntqid), tnkntq, aux);

    for (ipp=0;ipp<tnip;ipp++)
    {
      fwrite(&ipp, sizeof(ipp), 1, aux);
      for (j=0; j<nntq; j++)
        fwrite(&nontransq[j][ipp], sizeof(nontransq[j][ipp]), 1, aux);
    }
    fclose(aux);
  }
}



/**
   Function restores data from text backup file into regular integration points.
   
   @param aux - pointer to auxiliary file
   @param selelems - selection of elements whose eqother array will be restored
   @param selother - selection of components of saved eqother array which will be resored on selected elements
                     for each range or list item in selelems an individual selection of eqother components
                     is enabled
   @param selid    - array of indices of positions in eqother array to which will be 
                     restored selected saved eqother components
                     for each range or list item in selother an individual selection of position 
                     of eqother components is enabled
   
  @return The function does not return anything.

  JK, 19.9.2004
  TKo 9.2008
  TKr, 12.2008
*/
void transmat::restore_intpointst_txt (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long i, j, n, ir, ipp;
  long ntm, ncompgrad;
  sel selno; // selection with no selected components
  
  selno.n=1;
  for (i=0;i<tnip;i++){
    fscanf(aux, "%ld %ld %ld %ld", &ipp, &ntm, &ncompgrad, &n);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);
    if (ir < 0)
      ip[ipp].restore_data_txt(aux, n, selno, NULL);
    else
      ip[ipp].restore_data_txt(aux, n, selother[ir], selid[ir]);

    // restoring MEFEL quantities
    if (nontransq)
    {
      for (j=0; j<nntq; j++)
        fscanf(aux,"%le", &nontransq[j][ipp]);
    }
  }
}



/**
   Function restores data from several backup text files into regular integration points.
   
   @param selelems - selection of elements whose eqother array will be restored
   @param selother - selection of components of saved eqother array which will be resored on selected elements
                     for each range or list item in selelems an individual selection of eqother components
                     is enabled
   @param selid    - array of indices of positions in eqother array to which will be 
                     restored selected saved eqother components
                     for each range or list item in selother an individual selection of position 
                     of eqother components is enabled
   
  @return The function does not return anything.

  JK, 19.9.2004
  TKo 9.2008
  TKr, 12.2008
*/
void transmat::restore_intpointst_txt (sel &selelems, sel *selother, long **selid)
{
  long i, j, k, n, ir, ik, is, ipp;
  long ntm, ncompgrad;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  
  // restoring of gradient arrays
  sprintf(name, "%s.grad.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fscanf (aux,"%ld %ld %ld",&ipp, &ntm, &ncompgrad);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0;j<Tp->ntm;j++)
      for (k=0;k<ip[i].ncompgrad;k++)
	fscanf(aux, "%le", ip[i].grad[j]+k);
  }
  // restoring of flux arrays
  sprintf(name, "%s.flux.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fscanf(aux, "%ld %ld %ld", &ipp, &ntm, &ncompgrad);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Tp->ntm; j++)
      for (k=0;k<ip[ipp].ncompgrad;k++)
	fscanf(aux, "%le", ip[ipp].fluxes[j]+k);
  }
  // restoring of other arrays
  sprintf(name, "%s.other.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fscanf(aux,"%ld %ld",&ipp, &n);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);  
    for (k=0;k<n;k++)
    {
      fscanf (aux, "%le", &tmp);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k,ik))
      {
        is = selid[ir][ik]+k-selother[ir].id1[ik]; // storage index in eqother array 
        if (is >= ip[ipp].ncompeqother)
          print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
        else
          ip[ipp].eqother[is] = tmp;
      }
    }
    ipp++;
  }
  fclose(aux);

  // restoring MEFEL quantities
  if (nontransq)
  {
    sprintf(name, "%s.mefquant.bac",Tp->hdbcont.hdbnamer);
    aux = fopen(name,"rt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }

    // read table of indices of non-mechanical quantities
    for (i=0; i<tnkntq; i++)
      fscanf(aux, "%ld ", &ntqid[i]);

    for (i=0;i<tnip;i++)
    {
      fscanf (aux, "%ld", &ipp);
      if ((ipp < 0) || (ipp >= tnip))
      {
        print_err("invalid integration point number", __FILE__, __LINE__, __func__);
        abort();
      }
      for (j=0; j<nntq; j++)
        fscanf (aux,"%le", &nontransq[j][ipp]);
    }
    fclose(aux);
  }
}



/**
   Function restores data from binary backup file into regular integration points.
   
   @param aux - pointer to auxiliary file
   @param selelems - selection of elements whose eqother array will be restored
   @param selother - selection of components of saved eqother array which will be resored on selected elements
                     for each range or list item in selelems an individual selection of eqother components
                     is enabled
   @param selid    - array of indices of positions in eqother array to which will be 
                     restored selected saved eqother components
                     for each range or list item in selother an individual selection of position 
                     of eqother components is enabled

  @return The function does not return anything.

  JK, 19.9.2004
  TKo 9.2008
  TKr, 12.2008
*/
void transmat::restore_intpointst_bin (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long i, j, ipp, ntm, ncompgrad, n;
  long ir;
  sel selno; // selection with no selected components
  
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    fread(&n, sizeof(n), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);
    if (ir < 0)
      ip[ipp].restore_data_bin(aux, n, selno, NULL);
    else
      ip[ipp].restore_data_bin(aux, n, selother[ir], selid[ir]);

    // restoring MEFEL quantities
    if (nontransq)
    {
      for (j=0; j<nntq; j++)
        fread(&nontransq[j][ipp], sizeof(nontransq[j][ipp]), 1, aux);
    }
  }
}



/**
   Function restores data from several binary backup files into regular integration points.
   
   @param selelems - selection of elements whose eqother array will be restored
   @param selother - selection of components of saved eqother array which will be resored on selected elements
                     for each range or list item in selelems an individual selection of eqother components
                     is enabled
   @param selid    - array of indices of positions in eqother array to which will be 
                     restored selected saved eqother components
                     for each range or list item in selother an individual selection of position 
                     of eqother components is enabled

  @return The function does not return anything.

  JK, 19.9.2004
  TKo 9.2008
  TKr, 12.2008
*/
void transmat::restore_intpointst_bin (sel &selelems, sel *selother, long **selid)
{
  long i,j,k,n,ir,ik,is,ipp;
  long ntm, ncompgrad;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of gradient arrays
  sprintf(name, "%s.grad.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0;j<Tp->ntm;j++)
      fread(ip[i].grad[j], sizeof(*ip[i].grad[j]), ip[i].ncompgrad, aux);
  }
  // restoring of flux arrays
  sprintf(name, "%s.flux.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != ip[ipp].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0;j<Tp->ntm;j++)
      fread(ip[i].fluxes[j], sizeof(*ip[i].fluxes[j]), ip[i].ncompgrad, aux);
  }
  // restoring of other arrays
  sprintf(name, "%s.other.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&n, sizeof(n), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);    
    for (k=0;k<n;k++)
    {
      fread(&tmp, sizeof(tmp), 1, aux);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k,ik))
      {
        is = selid[ir][ik]+k-selother[ir].id1[ik]; // storage index in eqother array 
        if (is >= ip[ipp].ncompeqother)
          print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
        else
          ip[ipp].eqother[is] = tmp;
      }
    }   
  }
  fclose(aux);

  // restoring MEFEL quantities
  if (nontransq)
  {
    sprintf(name, "%s.mefquant.bac",Tp->hdbcont.hdbnamer);
    aux = fopen(name,"rb");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
      
    // read table of indices of non-mechanical quantities
    fread (ntqid, sizeof(*ntqid), tnkntq, aux);
      
    for (i=0;i<tnip;i++)
    {
      fread (&ipp, sizeof(ipp), 1, aux);
      if ((ipp < 0) || (ipp >= tnip))
      {
        print_err("invalid integration point number", __FILE__, __LINE__, __func__);
        abort();
      }
      for (k=0; k<nntq; k++)
        fread(&nontransq[k][ipp], sizeof(nontransq[k][ipp]), 1, aux);
    }
    fclose(aux);
  }
}



//*****************************************************
// Functions for backup of auxiliary integration points
//*****************************************************



/*
  Function saves data from auxiliary integration points into backup file
  in text format.
   
  @param aux - pointer to auxiliary file
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 22.12.2017
*/
void transmat::save_auxintpointst_txt (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother)
{
  long i, ir, app, nsc, eid;
  sel selno; // selection with no selected components
  
  selno.n=1;
  fprintf(aux, "\n");
  for (i=0; i<n; i++)
  {    
    app = ipm[i].app;
    if (app < 0)
      continue;
    eid = ipm[i].eid;

    selelems.presence_id(i, ir);
    if (ir < 0)
      nsc = 0;    
    else
      nsc = selother[ir].give_nselcomp(ip[app].ncompeqother);

    fprintf(aux,"%ld %ld %ld %ld\n",app, Tp->ntm, aip[app].ncompgrad, nsc);
    if (ir < 0)
      aip[app].save_data_txt(aux, selno);
    else
      aip[app].save_data_txt(aux, selother[ir]);
	  
    fprintf(aux, "\n");
  }
}



/**
  Function saves data from auxiliary integration points into several backup files
  in text format.
   
  @param ni       - time step id
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled

  @return The function does not return anything.

  Created by Tomas Koudelka, 22.12.2017
*/
void transmat::save_auxintpointst_txt (long ni, long n, ipmap *ipm, sel &selelems, sel *selother)
{
  long i, j, k, ir, eid, nsc, app;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Tp->hdbcont.prec;
  
  sprintf(name, "%s.%ld.agrad.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;

    fprintf(aux, "%ld %ld %ld\n", i, Tp->ntm, aip[app].ncompgrad);
    for (j=0; j<Tp->ntm; j++)
    {
      for (k=0; k<aip[app].ncompgrad; k++)
        fprintf(aux, "%.*le\n", prec, aip[app].grad[j][k]);
    }
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.aflux.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;

    fprintf(aux, "%ld %ld %ld\n", i, Tp->ntm, aip[app].ncompgrad);
    for (j=0;j<Tp->ntm;j++)
    {
      for (k=0; k<aip[app].ncompgrad; k++)
        fprintf(aux, "%.*le\n", prec, aip[app].fluxes[j][k]);
    }
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.aother.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {    
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;

    selelems.presence_id(i, ir);    
    if (ir < 0)
      nsc = 0;
    else
      nsc = selother[ir].give_nselcomp(aip[app].ncompeqother);

    fprintf (aux, "%ld %ld\n", app, nsc);
    if (nsc==0) 
      continue;
    for (k=0; k<aip[app].ncompeqother; k++)
    {
      if (selother[ir].presence_id(k))
        fprintf (aux, "%.*le\n", prec, aip[app].eqother[k]);
    }
  }
  fclose(aux);
}



/**
  Function saves data from auxiliary integration points into backup file
  in binary format.
   
  @param aux - pointer to auxiliary file
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 20.12.2017
*/
void transmat::save_auxintpointst_bin (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother)
{
  long i, ir, eid, nsc, app;
  sel selno; // selection with no selected components
  
  selno.n=1;
  for (i=0; i<n; i++)
  {    
    app = ipm[i].app;
    if (app < 0)
      continue;
    eid = ipm[i].eid;
    selelems.presence_id(i, ir);
      
    if (ir < 0)
      nsc = 0;
    else
      nsc = selother[ir].give_nselcomp(aip[app].ncompeqother);

    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
    fwrite(&aip[app].ncompgrad, sizeof(aip[app].ncompgrad), 1, aux);
    fwrite(&nsc, sizeof(nsc), 1, aux);
    if (ir < 0)
      aip[app].save_data_bin(aux, selno);
    else
      aip[app].save_data_bin(aux, selother[ir]);
	  
  }
}



/**
  Function saves data from auxiliary integration points into several backup files
  in binary format.
   
  @param ni       - time step id
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be saved
  @param selother - selection of components of saved eqother array which will be saved on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
   
  @return The function does not return anything.

  TKo 9.2008
  TKr, 12.2008
*/
void transmat::save_auxintpointst_bin (long ni, long n, ipmap *ipm, sel &selelems, sel *selother)
{
  long i, j, k, ir, eid, nsc, app;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  
  sprintf(name, "%s.%ld.agrad.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point
      continue;

    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
    fwrite(&aip[app].ncompgrad, sizeof(aip[app].ncompgrad), 1, aux);
    for (j=0; j<Tp->ntm; j++)
      fwrite(aip[app].grad[j], sizeof(*aip[app].grad[j]), aip[app].ncompgrad, aux);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.aflux.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point
      continue;

    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Tp->ntm, sizeof(Tp->ntm), 1, aux);
    fwrite(&aip[app].ncompgrad, sizeof(aip[app].ncompgrad), 1, aux);
    for (j=0;j<Tp->ntm;j++)
      fwrite(aip[app].fluxes[j], sizeof(*aip[app].fluxes[j]), aip[app].ncompgrad, aux);
  }
  fclose(aux);
  
  sprintf(name, "%s.%ld.aother.bac", Tp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {    
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;
    selelems.presence_id(i, ir);    
    if (ir < 0)
      nsc = 0;
    else
      nsc = selother[ir].give_nselcomp(aip[app].ncompeqother);
    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&nsc, sizeof(nsc), 1, aux);
    if (nsc==0) 
      continue;
    for (k=0; k<aip[app].ncompeqother; k++)
    {
      if (selother[ir].presence_id(k))
        fwrite(aip[app].eqother+k, sizeof(*aip[app].eqother), 1, aux);
    }
  }
  
  fclose(aux);
}



/**
  Function restores data from text backup file into auxiliary integration points.
   
  @param aux - pointer to auxiliary file
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be restored
  @param selother - selection of components of saved eqother array which will be resored on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
  @param selid    - array of indices of positions in eqother array to which will be 
                    restored selected saved eqother components
                    for each range or list item in selother an individual selection of position 
                    of eqother components is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 22.12.2017
*/
void transmat::restore_auxintpointst_txt (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother, long **selid)
{
  long i, ir, eid, app, tapp;
  long ntm, ncompgrad, ncompother;
  sel selno; // selection with no selected components
  
  selno.n=1;
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;

    fscanf(aux, "%ld %ld %ld %ld", &tapp, &ntm, &ncompgrad, &ncompother);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);
    if (ir < 0)
      aip[app].restore_data_txt(aux, ncompother, selno, NULL);
    else
      aip[app].restore_data_txt(aux, ncompother, selother[ir], selid[ir]);

  }
}



/**
  Function restores data from several text backup files into auxiliary integration points.
   
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be restored
  @param selelems - selection of elements whose eqother array will be restored
  @param selother - selection of components of saved eqother array which will be resored on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
  @param selid    - array of indices of positions in eqother array to which will be 
                    restored selected saved eqother components
                    for each range or list item in selother an individual selection of position 
                    of eqother components is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 22.12.2017
*/
void transmat::restore_auxintpointst_txt (long n, ipmap *ipm, sel &selelems, sel *selother, long **selid)
{
  long i, j, k, ir, ik, is, eid, tapp, app;
  long ntm, ncompgrad, ncompother;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  
  // restoring of gradient arrays
  sprintf(name, "%s.agrad.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;

    fscanf(aux, "%ld %ld %ld", &tapp, &ntm, &ncompgrad);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Tp->ntm; j++)
      for (k=0; k<aip[app].ncompgrad; k++)
	fscanf(aux, "%le", aip[app].grad[j]+k);
  }
  // restoring of flux arrays
  sprintf(name, "%s.aflux.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;

    fscanf(aux, "%ld %ld %ld", &tapp, &ntm, &ncompgrad);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Tp->ntm; j++)
      for (k=0; k<aip[app].ncompgrad; k++)
	fscanf(aux, "%le", aip[app].fluxes[j]+k);
  }
  // restoring of other arrays
  sprintf(name, "%s.aother.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;

    fscanf(aux, "%ld %ld", &tapp, &ncompother);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);  
    for (k=0; k<ncompother; k++)
    {
      fscanf(aux, "%le", &tmp);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k, ik))
      {
        is = selid[ir][ik]+k-selother[ir].id1[ik]; // storage index in eqother array 
        if (is >= aip[app].ncompeqother)
          print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
        else
          aip[app].eqother[is] = tmp;
      }
    }
  }
  fclose(aux);
}



/**
  Function restores data from binary backup file into auxiliary integration points.
   
  @param aux - pointer to auxiliary file
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be restored
  @param selother - selection of components of saved eqother array which will be resored on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
  @param selid    - array of indices of positions in eqother array to which will be 
                    restored selected saved eqother components
                    for each range or list item in selother an individual selection of position 
                    of eqother components is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka 22.12.2017
*/
void transmat::restore_auxintpointst_bin (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother, long **selid)
{
  long i, ntm, ncompgrad, ncompother, eid, app, tapp;
  long ir;
  sel selno; // selection with no selected components
  
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) //direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;

    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    fread(&ncompother, sizeof(ncompother), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);
    if (ir < 0)
      aip[app].restore_data_bin(aux, ncompother, selno, NULL);
    else
      aip[app].restore_data_bin(aux, ncompother, selother[ir], selid[ir]);
  }
}



/**
  Function restores data from several binary backup files into auxliary integration points.
   
  @param n   - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]
  @param selelems - selection of elements whose eqother array will be restored
  @param selother - selection of components of saved eqother array which will be resored on selected elements
                    for each range or list item in selelems an individual selection of eqother components
                    is enabled
  @param selid    - array of indices of positions in eqother array to which will be 
                    restored selected saved eqother components
                    for each range or list item in selother an individual selection of position 
                    of eqother components is enabled
   
  @return The function does not return anything.
   
  Created by Tomas Koudelka, 20.12.2017
*/
void transmat::restore_auxintpointst_bin (long n, ipmap *ipm, sel &selelems, sel *selother, long **selid)
{
  long i, j, k, ir, ik, is, eid, app, tapp;
  long ntm, ncompgrad, ncompother;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of gradient arrays
  sprintf(name, "%s.agrad.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) //direct mapping to the regular integration point
      continue;

    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of transported matter", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0;j<Tp->ntm;j++)
      fread(aip[app].grad[j], sizeof(*aip[app].grad[j]), aip[app].ncompgrad, aux);
  }
  // restoring of flux arrays
  sprintf(name, "%s.aflux.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) //direct mapping to the regular integration point
      continue;

    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&ntm, sizeof(ntm), 1, aux);
    fread(&ncompgrad, sizeof(ncompgrad), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ntm != Tp->ntm)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompgrad != aip[app].ncompgrad)
    {
      print_err("incompatible number of flux/grad components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Tp->ntm; j++)
      fread(aip[app].fluxes[j], sizeof(*aip[app].fluxes[j]), aip[app].ncompgrad, aux);
  }
  // restoring of other arrays
  sprintf(name, "%s.aother.bac", Tp->hdbcont.hdbnamer);
  aux = fopen(name, "rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)  // direct mapping to the regular integration point
      continue;
    eid = ipm[i].eid;

    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&ncompother, sizeof(ncompother), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);    
    for (k=0; k<ncompother; k++)
    {
      fread(&tmp, sizeof(tmp), 1, aux);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k, ik))
      {
        is = selid[ir][ik]+k-selother[ir].id1[ik]; // storage index in eqother array 
        if (is >= aip[app].ncompeqother)
          print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
        else
          aip[app].eqother[is] = tmp;
      }
    }   
  }
  fclose(aux);
}


/**
   function returns the extinction coefficient
   
   @param mt - type of material
   @param matid - id of the material
   
   JK, 26.7.2011
*/
double transmat::give_extinction_coeff (mattypet mt,long /*matid*/)
{
  double ecoeff;
  
  switch (mt){
  case isotransmat:{
    ecoeff = 0.93;
    break;
  }
  case damisotransmat:{
    ecoeff = 0.93;
    break;
  }
  case kunzel:{
    //ecoeff = kun[matid].give_extinction_coeff ();
    break;
  } 
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return ecoeff;
}


void transmat::freezing_thawing ()
{
  long i;
  double w,t,pt,rl,ru,dv,icevol=0.0;

  fprintf (Outt,"%le ",Tp->time);
  for (i=0;i<tnip;i++){
    t = ip[i].av[1];
    pt = ip[i].pv[1];
    
    w = ip[i].eqother[0];
    rl = ip[i].eqother[5];
    ru = ip[i].eqother[6];
    
    //dv = moisth[0].freezing_volume_increment (w,t,pt,rl,ru,icevol);
	dv = 0.0;
    
    if (i<40){
      //fprintf (Outt,"    %15.12le %15.12le",t-273.15,dv);
      //fprintf (Outt,"\n %ld  %le   %15.12le   %15.12le   %15.12le",i,Tp->time,t-273.15,pt,dv);
      //fprintf (Outt,"\n time %le  t %15.12le   pt %15.12le  dv %15.12le",Tp->time,t,pt,dv);
    }
/*    if (dv>1.0e-8)
      fprintf (Outt,"\n %le   %15.12le   %15.12le   %15.12le",Tp->time,t,pt,dv);
*/    
    
    ip[i].eqother[4]=dv;
    ip[i].eqother[5]=rl;
    ip[i].eqother[6]=ru;
    ip[i].eqother[7]=icevol;
    
  }
  //fprintf (Outt,"\n");
  
  actual_previous_change ();

}

