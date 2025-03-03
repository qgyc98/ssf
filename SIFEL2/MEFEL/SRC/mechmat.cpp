#include "mechmat.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "globmat.h"
#include "elemhead.h"
#include "vecttens.h"
#include "intpoints.h"
#include "element.h"
#include "node.h"
#include "elastisomat.h"
#include "elastortomat.h"
#include "elastgmat3d.h"
#include "elastgmat2d.h"
#include "splas1d.h"
#include "j2flow.h"
#include "microM4.h"
#include "microSIM.h"
#include "microfiber.h"
#include "mohrc.h"
#include "mohrcparab.h"
#include "boermat.h"
#include "drprag.h"
#include "doubdp.h"
#include "drprag2.h"
#include "camclay.h"
#include "camclaycoup.h"
#include "bbm.h"
#include "dsm.h"
#include "shefplast.h"
#include "hissplas.h"
#include "glasgmech.h"
#include "glasgowdam.h"
#include "creepdam.h"
#include "timeswmat.h"
#include "effstress.h"
#include "simviscous.h"
#include "lemaitre.h"
#include "ortodam.h"
#include "ortodam2.h"
#include "ortodamrot.h"
#include "anisodam.h"
#include "anisodamrot.h"
#include "scaldam.h"
#include "fixortodam.h"
#include "scaldamcc.h"
#include "aci78.h"
#include "cebfip78.h"
#include "graphmat.h"
#include "varelastisomat.h"
#include "elastisopdmat.h"
#include "geoelast.h"
#include "creep.h"
#include "creep_b3.h"
#include "creep_rspec.h"
#include "creepb.h"
#include "creep_dpl.h"
#include "creep_effym.h"
#include "shrinkmat.h"
#include "consol.h"
#include "winpast.h"
#include "nonlocmicroM4.h"
#include "therisomat.h"
#include "nonlocplast.h"
#include "nonlocdamg.h"
#include "damplast.h"
#include "tensor.h"
#include "visplast.h"
#include "chen.h"
#include "lenjonesmat.h"
#include "contactmat.h"
#include "cebfipcontactmat.h"
#include "damplastifacemat.h"
#include "plastifacemat.h"
#include "elemswitch.h"
#include "relaxeuroc.h"
#include "therisomat.h"
#include "shrinkmat.h"
//#include "hypoplast.h"
//#include "hypoplunsatexptherm.h"
#include "hypoplunsatexptherm2.h"
#include "umatunsat2.h"
#include "layplate.h"
#include "viselast.h"
#include "elasttime.h"
#include "isoviscous.h"
#include "inicd.h"
#include "thervolisomat.h"
#include "cusatismaterial.h"
#include "homogmatm.h"
#include "stacktrace.h"
#include "ipmap.h"
#include "lhsrhs.h"
#include "selection.h"
#include "mathem.h"
#include <math.h>
#include <string.h>




/**
  Constructor initializes data members to zero or default values.

  Created by JK
  Modified by TKo, TKr
*/
mechmat::mechmat ()
{
  //  number of material types
  nmt=0;
  //  total number of integration points
  tnip=0;
  // total number of auxiliary integration points
  tnaip=0;

  //  indicator of plastic regime in domain
  plast=0;

  // maximum number of component of other array at nodes
  max_ncompothern=0;
  // maximum number of component of eqother array on elements
  max_ncompothere=0;
  // maximum number of component of strain/stress arrays at nodes
  max_ncompstrn=0;
  // maximum number of component of strain/stress array on elements
  max_ncompstre=0;
  
  //  intergation points
  ip=NULL;
  // auxiliary integration points (for transfer of quantities in coupled problems)
  aip = NULL;
  //  initial conditions at regular integration points
  ic=NULL; 
  //  initial conditions at auxiliary integration points
  aip_ic=NULL; 
  //  element - integration point map
  elip = NULL;
  //  element - auxiliary integration point map (for transfer of quantities in coupled problems)
  elaip = NULL;


  eliso=NULL;  elgm2d=NULL; elgm3d=NULL;
  spl1d=NULL;  j2f=NULL;   mcoul=NULL;  mcpar=NULL;  
  boerm=NULL;  drprm=NULL;  ddpm=NULL;  drprm2=NULL; cclay=NULL;  cclayc=NULL;   bbm=NULL;  shpl=NULL;  hisspl=NULL;
  glasgmat=NULL;  glasgdam=NULL;  crdam=NULL;  tswmat=NULL; ntswmat=0L; effstr=NULL;
  mpM4=NULL;  mpSIM=NULL;  mpfib=NULL;  nmpM4=NULL;
  svis=NULL;  isovis=NULL;  svipl=NULL;   lmtr=NULL;
  scdam=NULL; scdamcc=NULL; anidam=NULL; anidamrot=NULL;  ortdam=NULL;
  ortdam2=NULL; ortdamrot=NULL;
  aci78mod=NULL;  cebfip78mod=NULL;
  grmat=NULL;   geoel=NULL;  veliso=NULL;  elisopd = NULL;
  crb3=NULL;
  crrs=NULL;
  crdpl=NULL;
  crbaz=NULL;  csol=NULL;   wpast=NULL;
  tidilat=NULL;  tvidilat=NULL;  relaxec=NULL;
  nlplast=NULL;   nldamg=NULL;
  hypopl=NULL;    hypoplustherm=NULL;
  cusmat=NULL;
  
  dampl = NULL;  visplas = NULL;  viselas=NULL; 
  eltimemat = NULL;
  chplast = NULL;
  lenjon = NULL;
  conmat=NULL;
  lplate=NULL;
  cebfipconmat=NULL;

  //  artificial material obtained from homogenization
  hommatm=NULL;

  //  material for shrinkage strain computation
  shmat=NULL;
  // plasticity interface material
  plastifm=NULL;
  // damage plasticity interface material
  damplifm=NULL;

  ipv = NULL;
  
  //  array containing eigenstrains id
  eigstrid = NULL;
  //  eigenstrains at integration points
  eigstrains = NULL;
  //  eigenstresses at integration points
  eigstresses = NULL;
  //  temperature strains at integration points
  tempstrains = NULL;
  //  array containing eigenstrains id at auxiliary int. points
  aip_eigstrid = NULL;
  //  eigenstrains at auxiliary integration points
  aip_eigstrains = NULL;
  //  eigenstresses at auxiliary integration points
  aip_eigstresses = NULL;
  //  temperature strains at auxiliary integration points
  aip_tempstrains = NULL;
  // macro-stresses
  mstress = NULL;
  //  array of values of non-mechanical quantities for regular integration points
  nonmechq = NULL;
  //  array of values of non-mechanical quantities for auxliary integration points
  aip_nonmechq = NULL;
  // the number of nonmechanical quantities
  nnmq = 0;
  // ordering of non-mechanical quantities
  nmqo = NULL;
  for (long i=0; i<tnknmq; i++)
    nmqid[i] = -1;
  
  mtype = NULL;
  numtype = NULL;
}



/**
  Destructor releases allocated memory of the mechmat object.

  Created by JK
  Modified by TKo, TKr
*/
mechmat::~mechmat ()
{
  long i;

  delete [] ip;
  delete [] elip;
  delete [] aip;
  delete [] elaip;
  delete [] eliso;     delete [] elgm2d;     delete [] elgm3d;

  delete [] spl1d;     delete [] j2f;        
  delete [] mcoul;     delete [] mcpar;  

  delete [] boerm;     delete [] drprm;       delete [] ddpm;       delete [] drprm2; delete [] cclay;
  delete [] bbm;
  delete [] cclayc;
  delete [] shpl;
  delete [] hisspl;

  delete [] glasgmat;  delete [] glasgdam;    delete [] crdam;  delete [] tswmat;  delete [] effstr;

  delete [] mpM4;      delete [] mpSIM;       delete [] mpfib;  delete [] nmpM4;

  delete [] svis;  delete [] isovis;
  delete [] lmtr;

  delete [] scdam;     delete [] scdamcc;    delete [] anidam;    delete [] anidamrot;
  delete [] ortdam;    delete [] ortdam2;    delete [] ortdamrot;
  
  delete [] aci78mod;  delete [] cebfip78mod;
  delete [] grmat;     delete [] geoel;       delete [] veliso;  delete [] elisopd;
  delete [] crb3;      delete [] crrs;        delete [] crdpl;
  delete [] crbaz;     delete [] csol;        delete [] wpast;
  delete [] tidilat;   delete [] tvidilat;    delete [] relaxec;
  delete [] nlplast;   delete [] nldamg;
  delete [] dampl;     delete [] visplas;     delete [] viselas; 
  delete [] chplast;
  delete [] lenjon;    delete [] conmat;      delete [] cebfipconmat;
  delete [] lplate;
  delete [] hypopl;
  delete [] hypoplustherm;
  delete [] plastifm;
  delete [] damplifm;
  delete [] cusmat;

  delete [] hommatm;

  //  material for shrinkage strain computation
  delete [] shmat;
  delete [] eltimemat;
  
  if (ic != NULL)
  {
    for (i = 0; i < tnip; i++)
      delete [] ic[i];
  }
  delete [] ic;
  if (aip_ic != NULL)
  {
    for (i = 0; i < tnaip; i++)
      delete [] aip_ic[i];
  }
  delete [] aip_ic;
  
  delete [] ipv;
  delete [] mstress;

  if (nonmechq != NULL){
    for (i=0;i<nnmq;i++){
      delete [] nonmechq[i];
    }
    delete [] nonmechq;
  }

  if (aip_nonmechq != NULL){
    for (i=0;i<nnmq;i++){
      delete [] aip_nonmechq[i];
    }
    delete [] aip_nonmechq;
  }

  delete [] nmqo;

  if (eigstrains != NULL){
    for (i=0;i<tnip;i++){
      delete [] eigstrains[i];
    }
    delete [] eigstrains;
  }

  if (eigstrid != NULL){
    for (i=0;i<tnip;i++){
      delete [] eigstrid[i];
    }
    delete [] eigstrid;
  }

  if (eigstresses != NULL){
    for (i=0;i<tnip;i++){
      delete [] eigstresses[i];
    }
    delete [] eigstresses;
  }
  if (tempstrains != NULL){
    for (i=0;i<tnip;i++){
      delete [] tempstrains[i];
    }
    delete [] tempstrains;
  } 
  
  if (aip_eigstrid != NULL){
    for (i=0;i<tnaip;i++){
      delete [] aip_eigstrid[i];
    }
    delete [] aip_eigstrid;
  }
  if (aip_eigstrains != NULL)
  {
    for (i=0; i<tnaip; i++)
      delete [] aip_eigstrains[i];
    delete [] aip_eigstrains;
  }

  if (aip_eigstresses != NULL){
    for (i=0;i<tnaip;i++){
      delete [] aip_eigstresses[i];
    }
    delete [] aip_eigstresses;
  }

  if (aip_tempstrains != NULL){
    for (i=0;i<tnaip;i++){
      delete [] aip_tempstrains[i];
    }
    delete [] aip_tempstrains;
  } 

  delete [] mtype;
  delete [] numtype;
}



/**
  Function returns number of integration points and determines pointers
  to integration points on elements.
   
  @return The function returns number of integration points.

  Created by JK, 8.7.2001
*/
long mechmat::intpnum (void)
{
  long i,j,k,n,nb;
  
  n=0;
  for (i=0;i<Mt->ne;i++){
    nb=Mt->give_nb (i);
    Mt->elements[i].ipp = new long* [nb];
    for (j=0;j<nb;j++){
      Mt->elements[i].ipp[j] = new long [nb];
    }
    for (j=0;j<nb;j++){
      for (k=0;k<nb;k++){
        Mt->elements[i].ipp[j][k]=n;
        n+=Mt->give_nip (i,j,k);
      }
    }
  }
  return n;
}



/**
  Function allocates integration points.
   
  @return The function does not return anything.

  Created by JK, 28.11.2006
*/
void mechmat::intpointalloc ()
{
  ip = new intpoints [tnip];
}



/**
  Function allocates initial conditions in integration points.
   
  @return The function does not return anything.

  Created by TKr, 08/12/2010
*/
void mechmat::alloc_ic ()
{
  long i,n;

  ic = new double* [tnip];
  for (i=0;i<tnip;i++){
    n = ip[i].ncompstr;
    ic[i] = new double[3+n];
  }
}



/**
  The function allocates array of non-mechanical quantities.
  Number of rows is given by parameter n and number of columns is 
  equal to the total number of integration points.

  @param[in] n - number of non-mechanical quantities

  @return The function does not return anything.

  Created by Tomas Koudelka, 7.10.2013
*/
void mechmat::alloc_nonmechq(long n)
{
  long i, j;

  nnmq = n;
  nonmechq = new double*[nnmq];

  for (i=0; i<nnmq; i++)
  {
    nonmechq[i] = new double[tnip];
    // set zero initial values of given quantity
    memset(nonmechq[i], 0, sizeof(*nonmechq[i])*tnip);    
  }  
  if (nmqo == NULL)
  {
    nmqo = new nonmechquant[nnmq];
    j=0;
    for (i=0; i<tnknmq; i++)
    {
      if (nmqid[i] > -1)
      {
        nmqo[j] = nonmechquant(i+1);
        j++;
      }
    }
    if (j != nnmq)
    {
      print_err("indices of used non-mechanical qunatities must be initialized before allocation", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  return;
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

  @param[in] n   - the number of components in the ipm array
  @param[in] ipm - the array of integration point mapping objects

  @return The function does not return anything but it changes the class
          data members of aip and elaip arrays and tnaip value.

  Created by Tomas Koudelka, 21.11.2017
*/
void mechmat::alloc_aux_intp(long n, ipmap *ipm)
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
    for (i=0; i<nnmq; i++)
      delete [] aip_nonmechq[i];
    delete [] aip_nonmechq;
  }

  if (tnaip)
  {
    aip = new intpoints[tnaip];
    // set material models, number of stress/strain components, stress/strain state indicators,
    // allocate stress/strain/other/eqother and non-mechanical quantity arrays
    aipinit(n, ipm);
  }
  else
  {
    // tnaip=0 =>zero number of auxiliary integration points
    aip = NULL;
    elaip = NULL;
    aip_nonmechq = NULL;
  }
}



/**
  Function reads material types and material id from the 
  opened text file.

  Parameters:
  @param[in] in - pointer to the opened XFILE
   
  Returns:
  @retval The function does not return anything.

  Created by JK, 8.7.2001
*/
void mechmat::readip (XFILE *in)
{
  long i;
  for (i=0;i<tnip;i++){
    ip[i].read (in);
  }
}



/**
  Function defines material models on integration points.
  Material models are read on elements.
   
  @return The function does not return anything.

  Created by JK, 28.11.2006
*/
void mechmat::intpointinit ()
{
  long i, j, ii, jj, k, nb, nip, ipp, ncomp;
  strastrestate ssst;
  
  //  element - integration point mapping
  elip = new long[tnip];
  
  for (i=0;i<Mt->ne;i++){
    //  the number of blocks
    nb=Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      //  number of components in arrays strains/stresses
      ncomp=Mt->give_ncomp (i,ii);
      //  strain/stress state (beam, plane stress, plane strain, etc.
      ssst=Mt->give_ssst (i,ii);

      
      for (jj=0;jj<nb;jj++){
	//  the number of integration points in the ii,jj block
        nip=Mt->give_nip(i,ii,jj);
        //  number of the first integration point in the ii,jj block
        ipp=Mt->elements[i].ipp[ii][jj];
	
        ip[ipp].hmt = 0;
        for (j=0;j<nip;j++){
	  
	  ip[ipp].ncompstr=ncomp;
	  ip[ipp].ssst=ssst;
	  
	  if (elip != NULL)
            elip[ipp] = i;
	  
          ip[ipp].nm  = Mt->elements[i].nm;
          ip[ipp].tm  = new mattype[ip[ipp].nm];
          ip[ipp].idm = new long[ip[ipp].nm];
	  
          if ((ssst == planestrain) || (ssst == planestress))
            ip[ipp].ncompstr=4;
	  
          for (k = 0; k < ip[ipp].nm; k++){
            ip[ipp].tm[k]  = Mt->elements[i].tm[k];
            ip[ipp].idm[k] = Mt->elements[i].idm[k];
            update_hmt(ipp, k);
	  }
	  
          ipp++;
        }
      }
    }
  }
}



/**
  Function initializes material models, stress/strain/other/eqother on auxiliary integration points.
  Material models are copied from the corresponding elements.
   
  The function initializes material models only for those auxiliary integration points
  in which ipmap components do not contain direct mapping to the regular integration point array 
  ip, i.e. ipm[i].ipp < 0.

  @param[in] n   - the number of components in the ipm array
  @param[in] ipm - the array of integration point mapping objects

  @return The function does not return anything but it may change internal data of auxiliary integration points.

  Created by TKo, 21.11.2017
*/
void mechmat::aipinit (long n, ipmap *ipm)
{
  long i, j, ii, eid, nlmid, ncompnl, app;
  strastrestate ssst;
  intpoints *tmp;
  
  //  element - auxiliary integration point mapping
  elaip = new long[tnaip];
  // initial condition values at aip 
  if (ic)
  {
    aip_ic = new double*[tnaip];
    memset(aip_ic, 0, sizeof(*aip_ic)*tnaip);
  }
  
  // swap regular and auxiliary integration point set due to update_hmt and give_ncomp(eq)other functions
  // now, auxiliary integration points will be used in mechmat functions
  tmp = ip;
  ip = aip;
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0)
      continue;
    eid = ipm[i].eid;
    aip[app].hmt = 0;
    if (elaip != NULL)
       elaip[app] = eid;
    aip[app].nm  = Mt->elements[eid].nm;
    aip[app].tm  = new mattype[aip[app].nm];
    aip[app].idm = new long[aip[app].nm];
    // copy material model definitions
    for (j=0; j<aip[app].nm; j++)
    {
      aip[app].tm[j]  = Mt->elements[eid].tm[j];
      aip[app].idm[j] = Mt->elements[eid].idm[j];
      update_hmt(app, j);
    } 
    // copy the number of stress/strain components and stress/strain state indicator
    aip[app].ncompstr = Mt->give_tncomp(eid);
    aip[app].ssst = ssst = Mt->give_ssst(eid, 0);
    if ((ssst == planestrain) || (ssst == planestress))
      aip[app].ncompstr=4;

    // Additional allocation for the nonlocal material models is neccessary
    ii = aip[app].hmt & 2; // the material is non-local
    ncompnl = 0;
    if (ii > 0)
    {
      nlmid = givenonlocid(app);
      ncompnl = give_num_averq(app, nlmid);
    }
    aip[app].alloc(Mb->nlc, app, ncompnl);
  }

  if (eigstrid)
    alloceigstrid(aip_eigstrid, aip, tnaip);

  if (eigstrains)
    alloceigstrains(aip_eigstrains, aip, tnaip);

  if (eigstresses)
    alloceigstresses(aip_eigstresses, aip, tnaip);    

  if (tempstrains)
    alloctempstrains(aip_tempstrains, aip, tnaip);

  // allocate array of values for non-mechanical quantities at auxiliary integration points
  if (nnmq) 
  {
    // some non-mechanical quantities are required in the material models
    // their number and ordering is the same as for regular integration point
    aip_nonmechq = new double*[nnmq];
    for (i=0; i<nnmq; i++)
    {
      aip_nonmechq[i] = new double[tnaip];
      // set zero initial values of given quantity
      memset(aip_nonmechq[i], 0, sizeof(*aip_nonmechq[i])*tnaip);    
    }  
  }
  // swap regular and auxiliary integration point set
  // now, regular integration points will be used in mechmat functions
  ip = tmp;
}



/**
  The function sets the bit flag array hmt according to im-th material model in the model chain
  defined on the given integration point.
  
  @param[in] ipp - integration point id
  @param[in] im  - material model id

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2.2016
*/
void mechmat::update_hmt(long ipp, long im)
{
  switch (ip[ipp].tm[im])
  {
    case therisodilat:
      ip[ipp].hmt ^= 1;
      break;
    case thervolisodilat:
      ip[ipp].hmt ^= 1;
      break;
    case nonlocplastmat:
      ip[ipp].hmt ^= 2;
      break;
    case nonlocdamgmat:
      ip[ipp].hmt ^= 2;
      break;
    case nonlocalmicroM4:
      ip[ipp].hmt ^= 2;
      break;
    case relaxationeuro:
      ip[ipp].hmt ^= 4;
      break;
    case effective_stress:
      ip[ipp].hmt ^= 8;
      break;
    case time_switchmat:
      ip[ipp].hmt = 0;
      if (tswmat[im].ami < 0)
        break;
      for(long i=tswmat[im].ami; i<tswmat[im].nam; i++)   
        update_hmt(ipp, i);
      break;
    default:
      break;
  }
}



/**
  The function initializes  eigenstrain/eigenstress arrays on auxiliary integration point

  @param[in] n - number of auxiliary integration points (dimension of array ipm)
  @param[in] ipm - the array of integration point mapping objects
  @param[in] time - actual time for the evaulation of eigenstrain/eigenstress component functions

  The function does not return anything but it changes content of arrays aip_eigstrid, aip_eigstrains and aip_eigstresses.

  Created by TKo, 10.2021
*/
void mechmat::init_aipeigstr(long n, ipmap *ipm, double time)
{
  long i, ipp, app, eid;
    
  if (aip == NULL) // there are no auxiliary int. points
    return;
  
  for (i = 0; i < n; i++) // for all mapped auxiliary integration points
  {
    if (ipm[i].ipp >= 0)  // direct mapping to regular integration point
      continue;
    app = ipm[i].app;
    eid = ipm[i].eid;
    ipp = Mt->elements[eid].ipp[0][0];    
    if (Gtm->leso[eid]==1){
      // copy function id of all components from the first integration point of the corresponding element
      memcpy(aip_eigstrid[app], eigstrid[ipp], sizeof(*eigstrid[ipp])*aip[app].ncompstr);
    }
  }
  // compute actual values of eigenstrains/eigenstresses
  Mb->aip_eigstrain_computation(n, ipm, time);
}


/**
  Function computes initial values at auxiliary integration points
  from initial nodal values.

  The function computes initial conditions only for those auxiliary integration points
  in which ipmap components do not contain direct mapping to the regular integration point array 
  ip, i.e. ipm[i].ipp < 0.

  @param[in] n - the number of components in the ipm array
  @param[in] ipm - the array of integration point mapping objects
   
  @return The function does not return anything but sets the aip_ic.

  Created by Tomas Koudelka, 28.11.2017
*/
void mechmat::inic_aipval(long n, ipmap *ipm)
{
  long i, j, k, l, nne, nval, eid, aux, app;
  long ncompstr, ncompeqother, nstra, nstre, idic;
  double ipval;
  ivector enod;
  vector anv, natcoord(ASTCKVEC(3));
  matrix  nodval;
  inictype *ictn;   // type of initial condition in the nodes

  if (aip == NULL) // there are no auxiliary int. points
    return;

  for (k = 0; k < n; k++) // for all mapped auxiliary integration points
  {
    if (ipm[k].ipp >= 0)  // direct mapping to regular integration point
      continue;

    app = ipm[k].app;
    eid = ipm[k].eid;
    if (Gtm->leso[eid]==1) 
    {      
      nne = Mt->give_nne(eid);
      reallocv(RSTCKIVEC(nne, enod));
      // find maximum number of initial values prescribed at nodes
      // initial displacement and velocity values are not taken into account
      Mt->give_elemnodes (eid, enod);
      nval = Mb->ico[enod[0]].nval;
      if (Mb->ico[enod[0]].type & inidisp)
        nval -= Gtm->give_ndofn(enod[0]);
      if (Mb->ico[enod[0]].type & iniderdisp)
        nval -= Gtm->give_ndofn(enod[0]);
      for (j = 0; j < nne; j++)
      {
        aux = Mb->ico[enod[j]].nval;
        if (Mb->ico[enod[j]].type & inidisp)
          aux -= Gtm->give_ndofn(enod[j]);
        if (Mb->ico[enod[j]].type & iniderdisp)
          aux -= Gtm->give_ndofn(enod[j]);
        if (nval < aux)
          nval = aux;
      }

      // prepare matrix of given nodal values of initial conditions on the element eid
      reallocm(RSTCKMAT(nne, nval, nodval));
      nullm(nodval);
      ictn = new inictype[nne];
      for (i = 0; i < nne; i++)
      {
        ictn[i] = Mb->ico[enod[i]].type;
        l = 0L;
        if (ictn[i] & inidisp)  // skip initial nodal displacements
          l += Gtm->give_ndofn(enod[i]);
        if (ictn[i] & iniderdisp)  // skip initial time derivatives of nodal displacements
          l += Gtm->give_ndofn(enod[i]);
        for (j = 0; j < nval; j++)
        {
          nodval[i][j] = Mb->ico[enod[i]].val[j+l]; // copy nodal value of initial conditions
        }
      }
      ncompstr = aip[app].ncompstr;
      ncompeqother = aip[app].ncompeqother;
      nstra = nstre = idic = 0;
      for (j = 0; j < nval; j++) // for all initial values
      {
        if ((ictn[0] & inistrain) && (j < ncompstr))
        {
          nstra++;
          continue;
        }
        if ((ictn[0] & inistress) && (j < nstra + ncompstr))
        {
          nstre++;
          continue;
        }
        if ((ictn[0] & iniother) && (j < nstra+nstre+ncompeqother))
          continue;

        // prepare nodal value vector of the j-th initial condition 
        reallocv(RSTCKVEC(nne, anv));  
        for(i = 0; i < nne; i++) 
          anv[i] = nodval[i][j]; // nodal values of j-th inital condition
        
        if ((ictn[0] & inicond) && (j < nval))
        {
          // prepare vector of natural coordinates
          natcoord[0] = ipm[k].xi;
          natcoord[1] = ipm[k].eta;
          natcoord[2] = ipm[k].zeta;

          // approximation of the initial condition of given quantity
          ipval = interpolate(eid, anv, natcoord);

          if (aip_ic[app] == NULL) // allocate initial condition if necessary
          {
            aip_ic[app] = new double[nval-j];
            memset(aip_ic[app], 0, sizeof(*aip_ic[app])*(nval-j)); 
          }
          aip_ic[app][idic] += ipval; // set the initial condition value at the app-th auxiliary integration point
          idic++;
        }
      }
      delete [] ictn;
    }
  }
}



/**
  Function initiates integration point variables.
  Allocation of strain/stress/eqother/other arrays is performed 
  for each int. point. If the temperature models are detected 
  and temprature load was defined in the MEFEL input file then
  nonmechq array is allocated due to storage of temperature and 
  initial temperature.
   
  @return The function does not return anything.

  Created by JK, 21.8.2001
*/
void mechmat::init_ip_2 (void)
{
  long i,n,ii,nlmid,ncompnl;
  


  // setup of non-mechanical quantities
  if ((Mp->temperature == 1) || (Mp->pore_press == 1)) 
  {
    // all non-mechanical quantities are defined in MEFEL
    n = search_reqnmq(nmqo);
    if (n)
    {
      for (i=0; i<n; i++)
        nmqid[nmqo[i]-1] = i;

      alloc_nonmechq(n);
    }
  }
  else
  {
    // allocation of nomechq and initialization is done in METR
  }

  // initialization and allocation for non-local material models
  for (i=0;i<tnip;i++)
  {
    ncompnl = 0;
    ii = ip[i].hmt & 2; // the material is non-local
    if (ii > 0)
    // Additional allocation for the nonlocal material models is neccessary
    {
      nlmid = givenonlocid(i);
      ncompnl = give_num_averq(i, nlmid);
    }

    ip[i].alloc(Mb->nlc,i,ncompnl);
  }
}



/**
  Function reads integration points, material characteristics,
  auxiliary points for strain and stress computation.
   
  @param[in] in - input stream
   
  @return The function does not return anything.

  Created by JK, 21.7.2001
*/
void mechmat::read (XFILE *in)
{

  //  reading of material characteristics
  readmatchar (in);

  //  computation of number of all integration points
  tnip = tnrip = intpnum ();
  if (Mespr==1){
    fprintf (stdout,"\n number of integration points  %ld",tnip);
  }

  //  allocation of integration points
  intpointalloc ();

  //  reading of integration points
  //readip (in);
  intpointinit ();

/*
  //  reading of strain and stress points
  if (Mp->straincomp==1){
    stra.read(in);
  }
  if (Mp->stresscomp==1){
    stre.read(in);
  }
  */
}



/**
  Function reads material characteristics from file.
   
  @param[in] in - input file stream
   
  @return The function does not return anything.

  Created by JK, 8.7.2001
  Modified by TKo, TKr
  Modified by TKo 26.6.2014 - added keywords and split into readmatchar and readmattype
*/
void mechmat::readmatchar (XFILE *in)
{
  long i;

  xfscanf (in, "%k%ld", "num_mat_types", &nmt);

  numtype = new long [nmt];
  mtype = new mattype [nmt];  

  if (nmt<1)
    print_err("wrong number of material types",__FILE__,__LINE__,__func__);
  
  if (Mespr==1)  fprintf (stdout,"\n number of different types of materials  %ld",nmt);

  for (i=0;i<nmt;i++)
  {
    xfscanf (in,"%k%m %k%ld", "mattype", &mattype_kwdset, &mtype[i], "num_inst", &numtype[i]);
    if (numtype[i]<1)
      print_err("wrong number of material characteristics",__FILE__,__LINE__,__func__);
    
    readmattype(in, mtype[i], numtype[i]);
  }
}



/**
  Function reads material characteristics from file.
   
  @param[in] in    - input file stream
  @param[in] mtype - material type
  @param[in] numt  - the number of individual sets of material parameters
   
  @return The function does not return anything.

  Created by TKo according to original version of readmatchar, 26.6.2014
*/
void mechmat::readmattype(XFILE *in, mattype mtype, long numt)
{
  long j,k;

  switch (mtype){
    case elisomat:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic isotropic materials  %ld",numt);
      eliso = new elastisomat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic isotropic material is out of range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        eliso[k-1].read (in);
      }
      break;
    }
    case elortomat:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic orthotropic materials  %ld",numt);
      elorto = new elastortomat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic orthtropic material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        elorto[k-1].read (in);
      }
      break;
    }
    case elgmat3d:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic anisotropic materials  %ld",numt);
      elgm3d = new elastgmat3d [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic anisotropic material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        elgm3d[k-1].read (in);
      }
      break;
    }
    case elgmat2d:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic anisotropic materials for 2D problems%ld",numt);
      elgm2d = new elastgmat2d [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic anisotropic material for 2D problems is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        elgm2d[k-1].read (in);
      }
      break;
    }

    case homomatm:{
      if (Mespr==1)  fprintf (stdout,"\n  number of material models obtained from homogenization %ld",numt);
      hommatm = new homogmatm [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("wrong number of material obtained from homogenization is required",__FILE__,__LINE__,__func__);
          abort();
        }
        hommatm[k-1].read (in);
      }
      break;
    }

    case simplas1d:{
      if (Mespr==1)  fprintf (stdout,"\n number of simple plastic 1D materials   %ld",numt);
      spl1d = new splas1d [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of simple plastic 1D material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        spl1d[k-1].read (in);
      }
      break;
    }
    case jflow:{
      if (Mespr==1)  fprintf (stdout,"\n number of J2 flow materials   %ld",numt);
      j2f = new j2flow [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of J2 flow material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        j2f[k-1].read (in);
      }
      break;
    }
    case microplaneM4:{
      if (Mespr==1)  fprintf (stdout,"\n number of microplane materials   %ld",numt);
      mpM4 = new microM4 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of microplane material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        mpM4[k-1].read (in);
      }
      break;
    }
    case nonlocalmicroM4:{
      if (Mespr==1)  fprintf (stdout,"\n number of nonlocal microplane materials   %ld",numt);
      nmpM4 = new nonlocmicroM4 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of nonlocal microplane material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        nmpM4[k-1].read (in);
      }
      break;
    }
    case microsimp:{
      if (Mespr==1)  fprintf (stdout,"\n number of microplane materials   %ld",numt);
      mpSIM = new microSIM [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of microplane material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        mpSIM[k-1].read (in);
      }
      break;
    }
    case microfibro:{
      if (Mespr==1)  fprintf (stdout,"\n number of microplane materials   %ld",numt);
      mpfib = new microfiber [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of microplane material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        mpfib[k-1].read (in);
      }
      break;
    }
    case mohrcoul:{
      if (Mespr==1)  fprintf (stdout,"\n number of Mohr-Coulomb materials   %ld",numt);
      mcoul = new mohrcoulomb [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Mohr-Coulomb material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        mcoul[k-1].read (in);
      }
      break;
    }
    case mohrcoulparab:{
      if (Mespr==1)  fprintf (stdout,"\n number of parabolic Mohr-Coulomb materials   %ld",numt);
      mcpar = new mohrcoulombpar [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of parabolic Mohr-Coulomb material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        mcpar[k-1].read (in);
      }
      break;
    }
    case boermaterial:{
      if (Mespr==1)  fprintf (stdout,"\n number of Boer materials   %ld",numt);
      boerm = new boermat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Boer material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        boerm[k-1].read (in);
      }
      break;
    }
    case druckerprager:{
      if (Mespr==1)  fprintf (stdout,"\n number of Drucker-Prager materials   %ld",numt);
      drprm = new drprag [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Drucker-Prager material is out of range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        drprm[k-1].read (in);
      }
      break;
    }
    case doubledrprager:{
      if (Mespr==1)  fprintf (stdout,"\n number of double Drucker-Prager materials   %ld",numt);
      ddpm = new doubdp [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of double Drucker-Prager material is out of range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        ddpm[k-1].read (in);
      }
      break;
    }
    case druckerprager2:{
      if (Mespr==1)  fprintf (stdout,"\n number of Drucker-Prager2 materials   %ld",numt);
      drprm2 = new drprag2 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Drucker-Prager material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        drprm2[k-1].read (in);
      }
      break;
    }
    case modcamclaymat:{
      if (Mespr==1)  fprintf (stdout,"\n number of modified cam-clay materials   %ld",numt);
      cclay = new camclay [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of modified cam-clay material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        cclay[k-1].read (in);
      }
      break;
    }
    case modcamclaycoupmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of modified cam-clay materials for coupling %ld",numt);
      cclayc = new camclaycoup [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of modified cam-clay material for coupling is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        cclayc[k-1].read (in);
      }
      break;
    }
    case bbmcoupmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of modified cam-clay materials for coupling %ld",numt);
      bbm = new bbmmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1)
          print_err("wrong number of modified cam-clay material for coupling",__FILE__,__LINE__,__func__);
        bbm[k-1].read (in);
      }
      break;
    }
    case doublestructuremat:{
      if (Mespr==1)  fprintf (stdout,"\n number of double structure materials for clays %ld",numt);
      dsm = new dsmmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1)
          print_err("wrong number of double structure materials for clays",__FILE__,__LINE__,__func__);
        dsm[k-1].read (in);
      }
      break;
    }
    case hypoplastmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of hypoplast materials for coupling %ld",numt);
      hypopl = new hypoplast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of hypoplast material for coupling is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        hypopl[k-1].read (in);
      }
      break;
    }
    case hypoplastusatthermat:{
      if (Mespr==1)  fprintf (stdout,"\n number of hypoplast materials with thermal effects for coupling %ld",numt);
      hypoplustherm = new hypoplast_unsat_exp_therm [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of hypoplast material with thermal effects for coupling is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        hypoplustherm[k-1].read (in);
      }
      break;
    }
    case shefpl:{
      if (Mespr==1)  fprintf (stdout,"\n number of Sheffield plasticity materials   %ld",numt);
      shpl = new shefplast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Sheffiled material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        shpl[k-1].read (in);
      }
      break;
    }
    case hissplasmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of HISS plasticity materials   %ld",numt);
      hisspl = new hissplas [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of HISS plasticity material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        hisspl[k-1].read (in);
      }
      break;
    }
    case chenplast:{
      if (Mespr==1)  fprintf (stdout,"\n number of Chen plasticity materials   %ld",numt);
      chplast = new chen [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Chen material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        chplast[k-1].read (in);
      }
      break;
    }
    case glasgowmechmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of Glasgow plasticity materials   %ld",numt);
      glasgmat = new glasgmech [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Glasgow material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        glasgmat[k-1].read (in);
      }
      break;
    }
    case glasgowdamage:{
      if (Mespr==1)  fprintf (stdout,"\n number of Glasgow damage materials   %ld",numt);
      glasgdam = new glasgowdam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Glasgow damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        glasgdam[k-1].read (in);
      }
      break;
    }
    case creep_damage:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial material creep-damage %ld",numt);
      crdam = new creepdam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial material creep-damge is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
      }
      break;
    }
    case time_switchmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial time-switched material %ld",numt);
      tswmat = new timeswmat [numt];
      ntswmat = numt;
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial time-switched material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        tswmat[k-1].read (in);
      }
      break;
    }
    case effective_stress:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial material for effective stresses %ld",numt);
      effstr = new effstress [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial material for effective stresses is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
	effstr[k-1].read (in);
      }
      break;
    }
    case simvisc:{
      if (Mespr==1)  fprintf (stdout,"\n number of simple viscous materials   %ld",numt);
      svis = new simviscous [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of simple viscous material in function is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        svis[k-1].read (in);
      }
      break;
    }
    case isovisc:{
      if (Mespr==1)  fprintf (stdout,"\n number of viscous materials for isotropic material  %ld",numt);
      isovis = new isoviscous [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of viscous material for isotropic material in function is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        isovis[k-1].read (in);
      }
      break;
    }
    case lemaitr:{
      if (Mespr==1)  fprintf (stdout,"\n number of Lemaitre visco materials   %ld",numt);
      lmtr = new lemaitre [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Lemaitre visco material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        lmtr[k-1].read (in);
      }
      break;
    }
    case scaldamage:{
      if (Mespr==1)  fprintf (stdout,"\n number of simple damage materials   %ld",numt);
      scdam = new scaldam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of simple damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        scdam[k-1].read (in);
      }
      break;
    }
    case fixortodamage:{
      if (Mespr==1)  fprintf (stdout,"\n number of orthotropic damage materials   %ld",numt);
      fixodam = new fixortodam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of orthotropic damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        fixodam[k-1].read (in);
      }
      break;
    }
    case scaldamagecc:{
      if (Mespr==1)  fprintf (stdout,"\n number of simple damage materials with crack closure  %ld",numt);
      scdamcc = new scaldamcc [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of simple damage material with crack closure is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        scdamcc[k-1].read (in);
      }
      break;
    }
    case ortodamage:{
      if (Mespr==1)  fprintf (stdout,"\n number of ortotropic damage materials   %ld",numt);
      ortdam = new ortodam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of ortotropic damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        ortdam[k-1].read (in);
      }
      break;
    }
    case ortodamage2:{
      if (Mespr==1)  fprintf (stdout,"\n number of ortotropic damage materials with general elastic material   %ld",numt);
      ortdam2 = new ortodam2 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of ortotropic damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        ortdam2[k-1].read (in);
      }
      break;
    }
    case ortodamagerot:{
      if (Mespr==1)  fprintf (stdout,"\n number of orthotropic damage materials with crack rotation   %ld",numt);
      ortdamrot = new ortodamrot [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of orthotropic damage material with crack rotation is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        ortdamrot[k-1].read (in);
      }
      break;
    }
    case anisodamage:{
      if (Mespr==1)  fprintf (stdout,"\n number of anisotropic damage materials   %ld",numt);
      anidam = new anisodam [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of anisotropic damage material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        anidam[k-1].read (in);
      }
      break;
    }
    case anisodamagerot:{
      if (Mespr==1)  fprintf (stdout,"\n number of rotational anisotropic damage materials with crack rotation   %ld",numt);
      anidamrot = new anisodamrot [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of anisotropic damage material with crack rotation is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        anidamrot[k-1].read (in);
      }
      break;
    }

    case aci:{
      if (Mespr==1)  fprintf (stdout,"\n number of aci78 materials   %ld",numt);
      aci78mod = new aci78 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of aci78 material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        aci78mod[k-1].read (in);
      }
      break;
    }
    case cebfip:{
      if (Mespr==1)  fprintf (stdout,"\n number of cebfip78 materials   %ld",numt);
      cebfip78mod = new cebfip78 [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of cebfip78 material in function is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        cebfip78mod[k-1].read (in);
      }
      break;
    }
    case graphm:{
      if (Mespr==1)  fprintf (stdout,"\n number of graphic material materials   %ld",numt);
      grmat = new graphmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of graphic material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        grmat[k-1].read (in);
      }
      break;
    }
    case varelisomat:{
      if (Mespr==1)  fprintf (stdout,"\n number of variable elastic isotropic materials  %ld",numt);
      veliso = new varelastisomat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of variable elastic isotropic material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        veliso[k-1].read (in);
      }
      break;
    }
    case elisopdmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic isotropic materials with pressure dependent initial stiffness  %ld",numt);
      elisopd = new elastisopdmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic isotropic material with pressure dependent stiffness is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        elisopd[k-1].read (in);
      }
      break;
    }
    case geoelast:{
      if (Mespr==1)  fprintf (stdout,"\n number of geoelast material materials   %ld",numt);
      geoel = new geoelastmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
	xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of graphic material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        geoel[k-1].read (in);
      }
      break;
    }
    case creepb3:{
      if (Mespr==1)  fprintf (stdout,"\n number of b3 creep material materials   %ld",numt);
      crb3 = new b3mat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of b3 creep material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        crb3[k-1].read (in);
      }
      break;
    }
    case creeprs:{
      if (Mespr==1)  fprintf (stdout,"\n number of b3 creep material materials   %ld",numt);
      crrs = new rspecmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of b3 creep material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        crrs[k-1].read (in);
      }
      break;
    }
    case creepdpl:{
      if (Mespr==1)  fprintf (stdout,"\n number of double-power law creep material materials   %ld",numt);
      crdpl = new dplmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of double-power law creep material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        crdpl[k-1].read (in);
      }
      break;
    }
    case creepbaz:{
      if (Mespr==1)  fprintf (stdout,"\n number of creep material materials   %ld",numt);
      crbaz = new creepb [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of creep material in function is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        crbaz[k-1].read (in);
      }
      break;
    }
    case creepeffym:{
      if (Mespr==1)  fprintf (stdout,"\n number of creep_effym materials   %ld",numt);
      creffym = new creep_effym [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err ("number %ld of creep_effym material in function is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        creffym[k-1].read(in);
      }
      break;
    }
    case consolidation:{
      if (Mespr==1)  fprintf (stdout,"\n number of consolidation material materials   %ld",numt);
      csol = new consol [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of consolidation material model is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        csol[k-1].read (in);
      }
      break;
    }
    case winklerpasternak:{
      if (Mespr==1)  fprintf (stdout,"\n number of Winkler-Pasternak material materials   %ld",numt);
      wpast = new winpast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Winkler-Pasternak material model is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        wpast[k-1].read (in);
      }
      break;
    }
    case therisodilat:{
      if (Mespr==1)  fprintf (stdout,"\n number of isotropic thermal dilatancy materials  %ld",numt);
      tidilat = new therisomat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of isotropic thermal dilatancy material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        tidilat[k-1].read (in);
      }
      break;
    }
    case thervolisodilat:{
      if (Mespr==1)  fprintf (stdout,"\n number of isotropic thermal and volume dilatancy materials  %ld",numt);
      tvidilat = new thervolisomat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of isotropic thermal and volume dilatancy material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        tvidilat[k-1].read (in);
      }
      break;
    }
    case relaxationeuro:{
      if (Mespr==1)  fprintf (stdout,"\n number of models of stress relaxation based on Eurocode2  %ld",numt);
      relaxec = new relaxeuroc [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of models of stress relaxation based on Eurocode2 is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        relaxec[k-1].read (in);
      }
      break;
    }

    case nonlocplastmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of nonlocal plasticity materials   %ld",numt);
      nlplast = new nonlocplast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of nonlocal plasticity material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        nlplast[k-1].read (in);
      }
      // nonlocal material model flag is set
      Mp->matmodel = nonlocal;
      break;
    }
    case nonlocdamgmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of nonlocal damage materials   %ld",numt);
      nldamg = new nonlocdamg [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of nonlocal damage material is out fo range [1;%ld]", __FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        nldamg[k-1].read (in);
      }
      // nonlocal material model flag is set
      Mp->matmodel = nonlocal;
      break;
    }

    case damage_plasticity:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial material damage-plasticity  %ld",numt);
      dampl = new damplast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial material damge-plasticity is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
      }
      break;
    }

    case viscoplasticity:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial material of viscoplasticity  %ld",numt);
      visplas = new visplast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial material of viscoplasticity is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
      }
      break;
    }
      
    case viscoelasticity:{
      if (Mespr==1)  fprintf (stdout,"\n number of artificial material of viscoelasticity   %ld",numt);
      viselas = new viselast [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of artificial material of viscoelasticity is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
      }
      break;
    }

    case lenjonespot:{
      if (Mespr==1)  fprintf (stdout,"\n number of Lennard-Jones interatomic potentials   %ld",numt);
      lenjon = new lenjonesmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of Lennard-Jones interatomic potentials is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
	lenjon[k-1].read (in);
      }
      break;
    }

    case cusatismat:{
      if (Mespr==1)  fprintf (stdout,"\n number of the Cusatis material models   %ld",numt);
      cusmat = new cusatismaterial [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of the Cusatis material models is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
	cusmat[k-1].read (in);
      }
      break;
    }
    
    case contmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of materials for contact problems   %ld",numt);
      conmat = new contactmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of materials for contact problems is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        conmat[k-1].read (in);
      }
      break;
    }

    case cebfipcontmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of CEB-FIP materials for contact problems   %ld",numt);
      cebfipconmat = new cebfipcontactmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of CEB-FIP materials for contact problems is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        cebfipconmat[k-1].read (in);
      }
      break;
    }

    case damplifmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of damage-plasticity materials for interface problems   %ld",numt);
      damplifm = new damplastifacemat[numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of damage-plasticity materials for interface problems is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        damplifm[k-1].read (in);
      }
      break;
    }

    case plastifmat:{
      if (Mespr==1)  fprintf (stdout,"\n number of plasticity materials for interface problems   %ld",numt);
      plastifm = new plastifacemat[numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of plasticity materials for interface problems is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        plastifm[k-1].read (in);
      }
      break;
    }

      
      
    case shrinkagemat:{
      if (Mespr==1)  fprintf (stdout,"\n number of materials for shrinkage description   %ld",numt);
      shmat = new shrinkmat [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of materials for shrinkage description is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        shmat[k-1].read (in);
      }
      break;
    }
      
    case elasttimemat:{
      if (Mespr==1)  fprintf (stdout,"\n number of elastic isotropic material with time dependent elastic modulus %ld",numt);
      eltimemat = new elasttime [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of elastic isotropic material with time dependent elastic modulus is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        eltimemat[k-1].read (in);
      }
      break;
    }

    case layerplate:
      if (Mespr==1)  fprintf (stdout,"\n number of layered plate materials   %ld",numt);
      lplate = new layplate [numt];
      for (j=0;j<numt;j++){
        k=numt+1;
        xfscanf (in,"%ld",&k);
        if (k>numt || k<1){
          print_err("number %ld of layered plate material is out fo range [1;%ld]",__FILE__,__LINE__,__func__, k, numt);
          abort();
        }
        lplate[k-1].read (in);
      }
      break;
    default:
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
}



/**
   The function prints integration points, material characteristics.
   
   @param[in] out - output stream
   
   @return The function does not return anything.

   02/01/2013, TKr
   Modified by TKo, 1.7.2014
*/
void mechmat::print (FILE *out)
{
  long i, j;
  //  printing of material characteristics
  fprintf(out,"\n## materials:\n");
  fprintf(out,"\n# number of material types\n");
  fprintf (out,"\n%ld\n",nmt);
  for (i=0;i<nmt;i++)
  {
    fprintf (out,"\n\n  %d %ld", mtype[i], numtype[i]);
    for (j=0;j<numtype[i];j++)
    {
      fprintf (out,"\n%ld ",j+1);
      printmatchar (out, mtype[i], j);
      fprintf(out, "\n");
    }
  }
}



/**
   The function prints material characteristics of the j-th material type of mt.
   
   @param[in] out - output stream
   @param[in] mt  - type of printed material
   @param[in] numinst - index of set of parameters for the given material type

   @return The function does not return anything.
   
   TKr 02/01/2013
   Modified by TKo, 1.7.2014
*/
void mechmat::printmatchar (FILE *out, mattype mt, long numinst)
{
  switch (mt){
    case elisomat:{
      eliso[numinst].print (out);
      break;
    }
    case elortomat:{
      elorto[numinst].print (out);
      break;
    }
    case elasttimemat:
      eltimemat[numinst].print (out);
      break;
    case simplas1d:{
      spl1d[numinst].print(out);
      break;
    }
    case jflow:{
      j2f[numinst].print(out);
      break;
    }
    case mohrcoul:{
      mcoul[numinst].print(out);
      break;
    }
    case druckerprager:{
      drprm[numinst].print (out);
      break;
    }

    case druckerprager2:{
      drprm2[numinst].print (out);
      break;
    }
    case modcamclaymat:{
      cclay[numinst].print(out);
      break;
    }
    case simvisc:{
      svis[numinst].print(out);
      break;
    }
    case isovisc:{
      isovis[numinst].print(out);
      break;
    }
    case hypoplastusatthermat:{
      hypoplustherm[numinst].print (out);
      break;
    }
    case creeprs:{
      crrs[numinst].print(out);
      break;
    }
    case creepeffym:{
      creffym[numinst].print(out);
      break;
    }
    case scaldamage:{
      scdam[numinst].print (out);
      break;
    }

    case scaldamagecc:{
      scdamcc[numinst].print (out);
      break;
    }

    case ortodamage:{
      ortdam[numinst].print (out);
      break;
    }
    case fixortodamage:{
      fixodam[numinst].print (out);
      break;
    }

    case damplifmat:
      damplifm[numinst].print(out);
      break;

    case plastifmat:
      plastifm[numinst].print(out);
      break;
      
    case nonlocdamgmat:{
      nldamg[numinst].print (out);
      break;
    }
    case winklerpasternak:{
      wpast[numinst].print (out);
      break;
    }
    case therisodilat:{
      tidilat[numinst].print (out);
      break;
    }
    case thervolisodilat:{
      tvidilat[numinst].print (out);
      break;
    }
    case viscoplasticity:{
      break;
    }
    case viscoelasticity:{
      break;
    }
    case time_switchmat:{
      tswmat[numinst].print(out);
      break;
    }
    case relaxationeuro:{
      relaxec[numinst].print (out);
      break;
    }
    case shrinkagemat:{
      shmat[numinst].print (out);
      break;
    }
    case graphm:{
      grmat[numinst].print (out);
      break;
    } 
    case varelisomat:{
      veliso[numinst].print (out);
      break;
    }
    case elisopdmat:{
      elisopd[numinst].print (out);
      break;
    }

    default:
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);

  } 
}



/**
  Function allocates array for function indeces of eigenstrain/eigenstress components.

  @param estrid[in,out] - pointer to array of function indeces for eigenstrain/eigenstress components which will be allocated
  @param iparray[in] - pointer to array of integration points which the estrid array is being associated with
  @param nip - the number of int. points in iparray (iparray dimension)
   
  @return The function returns pointer to the allocated array in the argument estrid.

  Created by TKo 10.2021
*/
void mechmat::alloceigstrid (long** &estrid, intpoints *iparray, long nip)
{
  long i, nc;
  
  if (estrid==NULL){
    estrid = new long*[nip];
    memset(estrid, 0, sizeof(*estrid)*nip);
    for (i=0; i<nip; i++){
      nc = iparray[i].ncompstr;
      estrid[i] = new long[nc];
      memset(estrid[i], 0, sizeof(*estrid[i])*nc);
    }
  }
}



/**
  Function allocates array for eigenstrains of regular int. points.

  @return The function does not return anything but it allocates eigenstrain array.

  Created by JK, 3.3.2004
  Rwritten by TKo, 25.01.2019
*/
void mechmat::alloceigstrains ()
{
  alloceigstrains(eigstrains, ip, tnip);
}



/**
  Function allocates array for eigenstrains.

  @param estra[in,out] - pointer to array of eigenstrains which will be alloacted
  @param iparray[in] - pointer to array of integration points which the eigenstrains are associated with
  @param nip - the number of int. points in iparray (iparray dimension)
   
  @return The function returns pointer to the allocated array of eigenstrains in the argument estra.

  Created by TKo according to JK, 25.01.2019
*/
void mechmat::alloceigstrains (double** &estra, intpoints */*iparray*/, long nip)
{
  long i, nc;
  
  if (estra==NULL){
    estra = new double* [nip];
    for (i=0;i<nip;i++){
      nc=ip[i].ncompstr;
      estra[i] = new double [nc];
      memset(estra[i], 0, sizeof(*estra[i])*nc);
    }
  }
}



/**
  Function set values of eigstrains to zero at regular integration points,
  it is used for nonlinear computations.

  @return The function does not return anything.

  18/9/2012, TKr
*/
void mechmat::nulleigstrains ()
{
  long i, nc;
  
  if (eigstrains != NULL){
    for (i=0;i<tnip;i++){
      nc=ip[i].ncompstr;
      memset(eigstrains[i], 0, sizeof(*eigstrains[i])*nc);
    }
  }
}



/**
  Function allocates array for strains due to temperature at regular int. points.

  @return The function does not return anything but it allocates temperature strain array.

  Created by JK, 3.3.2004
  Rewritten by TKo, 25.01.2019
*/
void mechmat::alloctempstrains ()
{
  alloctempstrains(tempstrains, ip, tnip);
}



/**
  Function allocates array for strains due to temperature.
   
  @param[in,out] tstra - pointer to array of temperature strains which will be alloacted
  @param[in] iparray - pointer to array of integration points which the temparture strains are associated with
  @param[in] nip - the number of int. points in iparray (iparray dimension)
   
  @return The function returns pointer to the allocated array of temperature strains in the argument tstra.

  Created by TKo according to JK, 25.01.2019
*/
void mechmat::alloctempstrains (double** &tstra, intpoints *iparray, long nip)
{
  long i, nc;
  
  if (tstra == NULL){
    tstra = new double* [nip];
    for (i=0;i<nip;i++){
      nc=iparray[i].ncompstr;
      tstra[i] = new double [nc];
      memset(tstra[i], 0, sizeof(*tstra[i])*nc);
    }
  }
}



/**
  Function set values of tempstrains to zero at regular integration points,
  it is used for nonlinear computations.

  @return The function does not return anything.

  13/5/2010, TKr
*/
void mechmat::nulltempstrains ()
{
  long i, nc;
  
  if (tempstrains != NULL){
    for (i=0;i<tnip;i++){
      nc=ip[i].ncompstr;
      memset(tempstrains[i], 0, sizeof(*tempstrains[i])*nc);
    }
  }
}

/**
  Function set values of tempstrains to zero at auxiliary integration points.
  It is intended for coupled problems esspecially.

  @return The function does not return anything.

  Created by Tomas Koudelka, 18.6.2018
*/
void mechmat::nullaiptempstrains ()
{
  long i, nc;
  
  if (aip_tempstrains != NULL){
    for (i=0;i<tnaip;i++){
      nc=aip[i].ncompstr;
      memset(aip_tempstrains[i], 0, sizeof(*aip_tempstrains[i])*nc);
    }
  }
}

void mechmat::nullstrains ()
{
  long i, nc;
  
  for (i=0;i<tnip;i++){
    nc=ip[i].ncompstr;
    memset(ip[i].strain, 0, sizeof(*ip[i].strain)*nc);
  }
}

void mechmat::nullstresses ()
{
  long i, nc;
  
  for (i=0;i<tnip;i++){
    nc=ip[i].ncompstr;
    memset(ip[i].stress, 0, sizeof(*ip[i].stress)*nc);
  }
}


/**
  Function allocates array for eigenstresses at regulat int. points.
   
  @return The function does not return anything but allocates eigenstress array 
          on regular int. points.

  Created by JK, 27.11.2006
  Rewritten by TKo, 25.01.2019
*/
void mechmat::alloceigstresses ()
{
  alloceigstresses(eigstresses, ip, tnip);
}



/**
  Function allocates array for strains due to temperature.
   
  @param[in,out] tstra - pointer to array of temperature strains which will be alloacted
  @param[in] iparray - pointer to array of integration points which the temparture strains are associated with
  @param[in] nip - the number of int. points in iparray (iparray dimension)
   
  @return The function returns pointer to the allocated array of temperature strains in the argument tstra.

  Created by TKo according to JK, 25.01.2019
*/
void mechmat::alloceigstresses (double** &estre, intpoints *iparray, long nip)
{
  long i, nc;
  
  if (estre == NULL){
    estre = new double* [nip];
    for (i=0;i<nip;i++){
      nc=iparray[i].ncompstr;
      estre[i] = new double [nc];
      memset(estre[i], 0, sizeof(*estre[i])*nc);
    }
  }
}



/**
  Function allocates array for macro-stresses. The array is used 
  in the case of macro-strain based approach to homogenization.
   
  @return The function does not return anything, it allocates and cleans the mstress array.

  Created by TKo, 17.2.2015
*/
void mechmat::allocmacrostresses ()
{
  if (mstress == NULL)
  {
    mstress = new double[max_ncompstre];
    memset(mstress, 0, sizeof(*mstress)*max_ncompstre);
  }
}



/**
  Function initializes material models with initial values.

  @param[in] lcid - load case id
  @param[in] rinit - flag for initialization after restorage from hdbackup
   
  @return The function does not return anything.

  Created by JK, 21.6.2004
*/
void mechmat::initmaterialmodels (long lcid, bool rinit)
{
  long i,j,nip,ipp;
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      nip = Mt->give_tnip (i);
      //  number of the first integration point
      ipp=Mt->elements[i].ipp[0][0];
      for(j=0; j<nip; j++){
	initvalues(lcid, ipp+j, 0, 0, rinit);
      }
    }
  }
}



/**
  Function initializes material models with initial values on auxiliary integration points.
  The function is intended for the transfer of values among meshes in coupled problems especially.

  @param[in] lcid - load case id
  @param[in] rinit - flag for initialization after restorage from hdbackup

  @return The function does not return anything but it chages state of auxiliary integration points
          at the array aip.

  Created by Tomas Koudelka, 24.11.2017
*/
void mechmat::aip_initmaterialmodels (long lcid, bool rinit)
{
  long app;
  intpoints *tmp_ip;
  long *tmp_elip;
  double **tmp_nmq;
  double **tmp_ic;
  double **tmp_eigstresses;
  double **tmp_eigstrains;
  double **tmp_tempstrains;
  
  if (aip == NULL) // there are no auxiliary int. points
    return;

  // switch regular integration point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = ip;
  tmp_elip = elip;
  tmp_nmq = nonmechq;
  tmp_ic = ic;
  tmp_eigstresses = eigstresses;
  tmp_eigstrains  = eigstrains;
  tmp_tempstrains = tempstrains;

  ip = aip;
  elip = elaip;
  nonmechq = aip_nonmechq;
  ic = aip_ic;
  eigstresses = aip_eigstresses;
  eigstrains  = aip_eigstrains;
  tempstrains = aip_tempstrains;

  for (app=0; app<tnaip; app++)
  {
    if (Gtm->leso[elip[app]] == 1)
      initvalues(lcid, app, 0, 0, rinit);
  }

  // switch regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  ip = tmp_ip;
  elip = tmp_elip;
  nonmechq = tmp_nmq;
  ic = tmp_ic;
  eigstresses = tmp_eigstresses;
  eigstrains  = tmp_eigstrains;
  tempstrains = tmp_tempstrains;
}



/**
  Function initializes values of internal material variables
  at integration point ipp for given material im.
   
  @param[in] lcid - load case id
  @param[in] ipp  - integration point pointer
  @param[in] im   - index of material type for given ip
  @param[in] ido  - index of internal variables for given material in the ipp eqother array
  @param[in] rinit - flag for initialization after restorage from hdbackup

  @return The function deos not return anything.

  Created by Tomas Koudelka   
*/
void mechmat::initvalues (long lcid, long ipp,long im,long ido, bool rinit)
{
  long ncompo;
  switch (ip[ipp].tm[im]){
    case elisomat:{
      eliso[ip[ipp].idm[im]].initval(ipp, im, ido);
      break;
    }
    case homomatm:
      hommatm[ip[ipp].idm[im]].initval(ipp, im, ido);
      break;
    case elortomat:
    case microplaneM4:
    case nonlocalmicroM4:
    case microsimp:
    case microfibro:
    case elasttimemat:
    case layerplate:
      break;
    case simplas1d:
    case jflow:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
      ncompo = givencompeqother (ipp, im);
      ncompo -= givencompeqother (ipp, im+1);
      initvalues(lcid, ipp, im+1, ido+ncompo, rinit);
      break;
    case modcamclaymat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      cclay[ip[ipp].idm[im]].initval(ipp, ido);
      break;
    case modcamclaycoupmat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      cclayc[ip[ipp].idm[im]].initval(ipp, im, ido);
      break;
    case bbmcoupmat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      bbm[ip[ipp].idm[im]].initval(ipp, im, ido);
      break;
    case doublestructuremat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      dsm[ip[ipp].idm[im]].initval(ipp, im, ido);
      break;
    case hypoplastmat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      hypopl[ip[ipp].idm[im]].initval(ipp, ido);
      break;
    case hypoplastusatthermat:
      if (Mp->strainstate == 0){
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
      hypoplustherm[ip[ipp].idm[im]].initval(ipp, ido, rinit);
      break;
    case shefpl:
    case chenplast:
    case glasgowmechmat:
    case glasgowdamage:
      break;
    case creep_damage:
      crdam[ip[ipp].idm[im]].initvalues(lcid, ipp, im, ido, rinit);
      break;
    case shrinkagemat:
    shmat[ip[ipp].idm[im]].initvalues(lcid, ipp, im, ido, rinit);
      break;
    case time_switchmat:
      tswmat[ip[ipp].idm[im]].initvalues(lcid, ipp, im, ido, rinit);
      break;
    case effective_stress:
      effstr[ip[ipp].idm[im]].initvalues(lcid, ipp, im, ido, rinit);
      break;
    case simvisc:
    case isovisc:
    case lemaitr:
    case scaldamage:
    case scaldamagecc:
    case aci:
    case cebfip:
    case graphm:
    case varelisomat:
    case geoelast:
      break;
    case elisopdmat:
      elisopd[ip[ipp].idm[im]].initval(ipp, ido);
      break;
    case fixortodamage:  
      fixodam[ip[ipp].idm[im]].initval(ipp, ido);
      break;
    case ortodamage:
      ortdam[ip[ipp].idm[im]].initvalues(ipp, ido);
      break;
    case ortodamage2:
      break;
    case ortodamagerot:
      ortdamrot[ip[ipp].idm[im]].initvalues(ipp, ido);
      break;
    case anisodamage:
      anidam[ip[ipp].idm[im]].initvalues(ipp, ido);
      break;
    case anisodamagerot:
      anidamrot[ip[ipp].idm[im]].initvalues(ipp, ido);
      break;
    case contmat:
    case cebfipcontmat:
    case damplifmat:
    case plastifmat:
      break;
    case creeprs:
    case creepb3:
    case creepdpl:{
      creep_initmaterialmodel(ipp,im,ido);
    }
    case creepbaz:
    case consolidation:
    case winklerpasternak:
    case therisodilat:
    case thervolisodilat:
    case relaxationeuro:
    case nonlocplastmat:
    case nonlocdamgmat:
    case damage_plasticity:
    case viscoplasticity:
    case viscoelasticity:
    case creepeffym:
    case rcmatmodnorm:
      break;
    default:
      print_err("unknown material type is required", __FILE__,__LINE__,__func__);
  }
}



/**
  The function returns the number of the normal components in the stress/strain state defined at 
  the given int. point, i.e. number of the normal componnets in int. point arrays strain/stress.

  @param[in] ipp - integration point pointer

  @return The function returns the number of the normal stress/strain components.

  Created by TKo, 02.2024
*/
long mechmat::give_num_norm_comp(long ipp)
{
  strastrestate ssst = ip[ipp].ssst;
  long ret = 0;
  switch(ssst){
    case bar:
      ret = 1;
      break;
      //    case plbeam:
      //    case spacebeam:
    case planestress:
    case planestrain:
      ret = 2;
      break;
    case planecontact:
      ret = 1;
      break;
      //    case platek:
      //    case plates:,
      //    case shell:,
    case axisymm:
    case spacestress:
      ret = 3;
      break;
    default:
      print_err("unknown stress/strain state(%d) is required", __FILE__, __LINE__, __func__, int(ssst));
  }
  return ret;
}






// *****************************************************************************
// *****************************************************************************
//
// Function used for access to stress/strain/eqother/other arrays on int. points
//
// *****************************************************************************
// *****************************************************************************

/**
  The function stores all strains.

  @param[in] lcid - load case id
  @param[in] ipp  - integration point pointer
  @param[in] eps  - %vector containing strain components
  
  @return The function does not return anything.
  
  Created by JK, 17.7.2001
*/
void mechmat::storestrain (long lcid,long ipp,vector &eps)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=eps.n;
  j=0;
  for (i=m;i<m+n;i++){
    ip[ipp].strain[i]=eps[j];
    j++;
  }
}



/**
  Function stores part of strains.

  @param[in] lcid - load case id
  @param[in] ipp  - integration point number
  @param[in] fi   - first index
  @param[in] eps  - %vector containing strain components
   
  @return The function does not return anything.

  Created by JK, 17.7.2001
*/
void mechmat::storestrain (long lcid,long ipp,long fi,vector &eps)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=eps.n;
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    ip[ipp].strain[i]=eps[j];
    j++;
  }
}



/**
   Function stores part of strains.
   ncomp is explicitly defined in order to save part of the %vector eps
   
   @param[in] lcid  - load case id
   @param[in] ipp   - integration point number
   @param[in] fi    - first index
   @param[in] ncomp - number of components
   @param[in] eps   - %vector containing strain components
   
   @return The function does not return anything.
   
   Created by JK, 17.7.2001
*/
void mechmat::storestrain (long lcid,long ipp,long fi,long ncomp,vector &eps)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    ip[ipp].strain[i]=eps[j];
    j++;
  }
}



/**
   Function stores part of strains.
   ncomp is explicitly defined in order to save part of the array eps
   
   @param[in] lcid  - load case id
   @param[in] ipp   - integration point number
   @param[in] fi    - first index
   @param[in] ncomp - number of components
   @param[in] eps   - array containing strain components
   
   @return The function does not return anything.
   
   Created by TKo, 29.11.2013
*/
void mechmat::storestrain (long lcid, long ipp, long fi, long ncomp, double *eps)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    ip[ipp].strain[i]=eps[j];
    j++;
  }
}



/**
  Function restores all strains to %vector eps.

  @param[in]  lcid - load case id
  @param[in]  ipp  - integration point pointer
  @param[out] eps  - %vector containing strain components

  @return The function returns %vector of strains in the parameter eps.
   
  Created by JK, 5.8.2001
*/
void mechmat::givestrain (long lcid,long ipp,vector &eps)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=eps.n;
  j=0;
  for (i=m;i<n+m;i++){
    eps[j]=ip[ipp].strain[i];
    j++;
  }
}



/**
  Function restores part of strains to %vector eps.

  @param[in]  lcid - load case id
  @param[in]  ipp  - integration point number
  @param[in]  fi   - first index
  @param[out] eps  - %vector containing strain components (output)

  @return The function returns part of the strain %vetcor in the parameter eps.
   
  Created by JK, 5.8.2001
*/
void mechmat::givestrain (long lcid,long ipp,long fi,vector &eps)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=eps.n;
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    eps[j]=ip[ipp].strain[i];
    j++;
  }
}



/**
  Function restores part of strains to %vector eps.
  ncomp is explicitly defined in order to restore part of the %vector eps

  Parameters:   
  @param[in]  lcid  - load case id
  @param[in]  ipp   - integration point number
  @param[in]  fi    - first index
  @param[in]  ncomp - number of components
  @param[out] eps   - %vector containing strain components
   
  @return The function returns part of the strain %vector in the parameter eps. 
  
  Created by JK, 5.8.2001
*/
void mechmat::givestrain (long lcid,long ipp,long fi,long ncomp,vector &eps)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    eps[j]=ip[ipp].strain[i];
    j++;
  }
}



/**
  Function cleans array of strains.
  It works for nonlinear problems where only one load case can be used.
  
  @return The function does not return anything.
 
  Created by JK
*/
void mechmat::cleanstrain ()
{
  long i;
  
  for (i=0;i<tnip;i++){
    //  it works for nonlinear problems where only one load case can be used
    ip[i].clean_strains (1);
  }
}



/**
  Function stores all stresses.

  @param[in] lcid - load case id
  @param[in] ipp  - integration point number
  @param[in] sig  - %vector containing stress components

  @return The function does not return anything.
    
  Created by JK, 17.7.2001
*/
void mechmat::storestress (long lcid,long ipp,vector &sig)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=sig.n;
  j=0;
  for (i=m;i<m+n;i++){
    ip[ipp].stress[i]=sig[j];
    j++;
  }
}



/**
  Function stores part of stresses.

  @param[in] lcid - load case id
  @param[in] ipp  - integration point number
  @param[in] fi   - first index
  @param[in] sig  - %vector containing stress components

  @return The function does not return anything.

  Created by JK, 17.7.2001
*/
void mechmat::storestress (long lcid,long ipp,long fi,vector &sig)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=sig.n;
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    ip[ipp].stress[i]=sig[j];
    j++;
  }
}



/**
  Function stores part of stresses.
  ncomp is explicitly defined in order to save part of the %vector sig

  Parameters:   
  @param[in] lcid  - load case id
  @param[in] ipp   - integration point number
  @param[in] fi    - first index
  @param[in] ncomp - number of components
  @param[in] sig   - %vector containing stress components

  @return The function does not return anything.

  Created by JK, 17.7.2001
*/
void mechmat::storestress (long lcid,long ipp,long fi,long ncomp,vector &sig)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    ip[ipp].stress[i]=sig[j];
    j++;
  }
}



/**
  Function stores part of stresses.
  ncomp is explicitly defined in order to save part of the array sig

  Parameters:   
  @param[in] lcid  - load case id
  @param[in] ipp   - integration point number
  @param[in] fi    - first index
  @param[in] ncomp - number of components
  @param[in] sig   - array containing stress components

  @return The function does not return anything.

  Created by TKo, 29.11.2013
*/
void mechmat::storestress (long lcid, long ipp, long fi, long ncomp, double *sig)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    ip[ipp].stress[i]=sig[j];
    j++;
  }
}



/**
  Function restores all stresses to %vector sig.

  @param[in]  lcid - load case id
  @param[in]  ipp  - integration point pointer
  @param[out] sig  - %vector containing stress components

  @return The function returns stress %vector in the parameter sig.
   
  Created by JK, 5.8.2001
*/
void mechmat::givestress (long lcid,long ipp,vector &sig)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=sig.n;
  j=0;
  for (i=m;i<m+n;i++){
    sig[j]=ip[ipp].stress[i];
    j++;
  }
}



/**
  Function restores part of stresses to %vector sig.

  @param[in]  lcid - load case id
  @param[in]  ipp  - integration point number
  @param[in]  fi   - first index
  @param[out] sig  - %vector containing stress components
   
  @return The function returns of the stress %vector in the parameter sig.

  Created by JK, 5.8.2001
*/
void mechmat::givestress (long lcid,long ipp,long fi,vector &sig)
{
  long i,j,m,n;
  m=lcid*ip[ipp].ncompstr;
  n=sig.n;
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    sig[j]=ip[ipp].stress[i];
    j++;
  }
}



/**
  Function restores part of stresses to %vector sig.
  ncomp is explicitly defined in order to restore part of the %vector sig

  Parameters:   
  @param[in]  lcid  - load case id
  @param[in]  ipp   - integration point number
  @param[in]  fi    - first index
  @param[in]  ncomp - number of components
  @param[out] sig   - %vector containing stress components
   
  @return The function returns of the stress %vector in the parameter sig.
   
  Created by JK, 5.8.2001
*/
void mechmat::givestress (long lcid,long ipp,long fi,long ncomp,vector &sig)
{
  long i,j,m;
  m=lcid*ip[ipp].ncompstr;
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    sig[j]=ip[ipp].stress[i];
    j++;
  }
}



/**
  Function stores eigenstrains at integration point.

  @param[in] ipp - number of integration point
  @param[in] eps - %vector of eigenstrains

  @return The function does not return anything.
   
  Created by JK, 3.3.2004
*/
void mechmat::storeeigstrain (long ipp,vector &eps)
{
  long i,nc;
  
  nc=ip[ipp].ncompstr;
  if (est==eigstrain){
    for (i=0;i<nc;i++){
      eigstrains[ipp][i]=eps[i];
    }
  }
  if (est==tempstrain){
    for (i=0;i<nc;i++){
      tempstrains[ipp][i]=eps[i];
    }
  }
}



/**
  Function stores eigenstrains.

  @param[in] ipp   - integration point pointer
  @param[in] fi    - first index
  @param[in] ncomp - number of required components
  @param[in] eps   - %vector containing eigenstrain components
   
  @return The function does not return anything.

  Created by JK, 17.7.2001
*/
void mechmat::storeeigstrain (long ipp,long fi,long ncomp,vector &eps)
{
  long i,j;

  if (est==eigstrain){
    j=0;
    for (i=fi;i<fi+ncomp;i++){
      eigstrains[ipp][i] = eps[j];
      j++;
    }
  }
  if (est==tempstrain){
    j=0;
    for (i=fi;i<fi+ncomp;i++){
      tempstrains[ipp][i] = eps[j];
      j++;
    }
  }
}



/**
  Function returns eigenstrains.

  Parameters:   
  @param[in]  ipp - integration point pointer
  @param[out] eps - %vector containing strain components
   
  @return The function returns %vector of eigenstrains in the parameter eps.

  Created by JK, 17.7.2001
*/
void mechmat::giveeigstrain (long ipp,vector &eps)
{
 long i,n;
  
  n=eps.n;  
  if (est==eigstrain){
    for (i=0;i<n;i++){
      eps[i] = eigstrains[ipp][i];
    }
  }
  if (est==tempstrain){
    for (i=0;i<n;i++){
      eps[i] = tempstrains[ipp][i];
    }
  }
}



/**
  Function restores part of eigenstrains.

  Parameters:   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of selected components
  @param[out] eps   - %vector containing strain components
   
  @return The function returns part of eigenstrain %vector in the parameter eps.

  Created by JK, 17.7.2001
*/
void mechmat::giveeigstrain (long ipp,long fi,long ncomp,vector &eps)
{
  long i,j;
  
  if (est==eigstrain){
    j=0;
    for (i=fi;i<fi+ncomp;i++){
      eps[j] = eigstrains[ipp][i];
      j++;
    }
  }
  if (est==tempstrain){
    j=0;
    for (i=fi;i<fi+ncomp;i++){
      eps[j] = tempstrains[ipp][i];
      j++;
    }
  }
}



/**
  Function stores eigenstresses at integration point.

  @param[in] ipp - number of integration point
  @param[in] sig - %vector of eigenstresses

  @return The function does not return anything.
   
  Created by JK, 27.11.2006
*/
void mechmat::storeeigstress (long ipp,vector &sig)
{
  long i,nc;
  
  nc=sig.n;
  for (i=0;i<nc;i++){
    eigstresses[ipp][i]=sig(i);
  }
  
}



/**
  Function stores eigenstresses.

  Parameters:   
  @param[in] ipp   - integration point pointer
  @param[in] fi    - first index
  @param[in] ncomp - number of required components
  @param[in] sig   - %vector containing eigenstress components
   
  @return The function does not return anything.

  Created by JK, 27.11.2006
*/
void mechmat::storeeigstress (long ipp,long fi,long ncomp,vector &sig)
{
  long i,j;

  j=0;
  for (i=fi;i<fi+ncomp;i++){
    eigstresses[ipp][i] = sig(j);
    j++;
  }
}



/**
  Function returns eigenstresses.
   
  @param[in]  ipp - integration point pointer
  @param[out] sig - %vector containing stress components
   
  @return The function returns %vector of eigenstresses in the parameter sig.

  Created by JK, 27.11.2006
*/
void mechmat::giveeigstress (long ipp,vector &sig)
{
  long i,n;
  
  n=sig.n;  
  for (i=0;i<n;i++){
    sig(i) = eigstresses[ipp][i];
  }
}



/**
  Function restores eigenstress.
   
  @param[in]  ipp - integration point pointer
  @param[in]  fi  - first index
  @param[out] sig - %vector containing stress components

  @return The function returns part of eigenstress %vector in the parameter sig.
   
  Created by JK, 29.4.2008
*/
void mechmat::giveeigstress (long ipp,long fi,vector &sig)
{
  long i,j,ncomp;
  ncomp=sig.n;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    sig(j) = eigstresses[ipp][i];
    j++;
  }
}



/**
  Function restores part of eigenstress %vector.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of selected components
  @param[out] sig   - %vector containing eigenstress components
   
  @return The function returns part of eigenstress %vector in the parameter sig.

  Created by JK, 29.4.2008
*/
void mechmat::giveeigstress (long ipp,long fi,long ncomp,vector &sig)
{
  long i,j;
  j=fi;
  for (i=fi;i<fi+ncomp;i++){
    sig(j) = eigstresses[ipp][i];
    j++;
  }
}



/**
  Function stores components of array eqother.
   
  @param[in] ipp   - integration point pointer
  @param[in] fi    - first index
  @param[in] ncomp - number of required components
  @param[in] comp  - %vector containing components
   
  @return The function does not return anything.

  Created by JK, 5.11.2002
*/
void mechmat::storeeqother (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].eqother[i]=comp[j];
    j++;
  }
}



/**
  Function restores part of array eqother to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of required components
  @param[out] comp  - array containing components

  @return The function returns part of eqother array in the parameter comp.
   
  Created by JK, 5.11.2002
*/
void mechmat::giveeqother (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    comp[j]=ip[ipp].eqother[i];
    j++;
  }
}



/**
  Function restores part of array eqother to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of required components
  @param[out] comp  - %vector for the component storage

  @return The function returns part of eqother array in the parameter comp.
   
  Created by TKo, 09.2020
*/
void mechmat::giveeqother (long ipp,long fi,long ncomp,vector &comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    comp(j)=ip[ipp].eqother[i];
    j++;
  }
}



/**
  Function restores part of array other to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of required components
  @param[out] comp  - array containing components

  @return The function returns part of other array in the parameter comp.
   
  Created by TKo, 6.7.2018
*/
void mechmat::giveother (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    comp[j]=ip[ipp].other[i];
    j++;
  }
}



/**
  Function restores part of array other to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of required components
  @param[out] comp  - %vector for the component storage

  @return The function returns part of other array in the parameter comp.
   
  Created by TKo, 09.2020
*/
void mechmat::giveother (long ipp,long fi,long ncomp,vector &comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    comp(j)=ip[ipp].other[i];
    j++;
  }
}



/**
  Function restores part of array other to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  selcomp   - selection of components
  @param[out] comp  - %vector for the component storage

  @return The function returns part of other array in the argument comp.
   
  Created by TKo, 09.2020
*/
void mechmat::giveother (long ipp, sel &selcomp, vector &comp)
{
  long i,j;
  j = 0;
  if ((selcomp.st == sel_list) || (selcomp.st == sel_single)){
    for (i=0; i<selcomp.n; i++){
      comp(j) = ip[ipp].other[selcomp.id1[i]];
      j++;
    }
  }
  else{
    for (i=0; i<ip[ipp].ncompother; i++){
      if (selcomp.presence_id(i)){
        comp(j) = ip[ipp].other[i];
        j++;
      }
    }
  }
}



/**
  Function stores components of array other.
   
  @param[in] ipp - integration point pointer
  @param[in] fi - first index
  @param[in] ncomp - number of required components
  @param[in] comp - array containing components
   
  @return The function does not return anything.

  Created by TKo, 29.11.2013
*/
void mechmat::storeother (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].other[i]=comp[j];
    j++;
  }
}



/**
  The function returns the number of components of ipp's other array.
  Material type is given by the im parameter, which means index of material
  in the array tm of given integration point ipp.
   
  @param[in] ipp - integration point pointer
  @param[in] im - index of material
   
  @return The function returns total number of components of other array.  

  Created by Tomas Koudelka, 5.11.2003
*/
long mechmat::givencompother (long ipp,long im)
{
  long ncompo = 0;
  long ncompstr = ip[ipp].ncompstr;
  
  switch (ip[ipp].tm[im]){
    case elisomat:
    case elortomat:
    case elgmat3d:
    case elgmat2d:
    case winklerpasternak:
    case graphm:
    case varelisomat:
    case elisopdmat:
    case elasttimemat:
    case geoelast:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case modcamclaymat:
    case modcamclaycoupmat:
    case bbmcoupmat:
    case doublestructuremat:
    case hypoplastmat:
    case hypoplastusatthermat:
    case shefpl:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepbaz:
    case creepeffym:
    case contmat:
    case cebfipcontmat:
    case damplifmat:
    case plastifmat:
    case consolidation:
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
    case glasgowdamage:
    case damage_plasticity:
    case therisodilat:
    case thervolisodilat:
    case relaxationeuro:
    case simvisc:
    case isovisc:
    case lenjonespot:
    case cusatismat:{
      ncompo = givencompeqother (ipp,im);
      break;
    }
    case homomatm:
      break;
    case creeprs:
    case creepb3:
    case creepdpl:{
      ncompo=creep_ncompo(ipp,im);//the same as in eqother
      break;
    }
    case layerplate:
      ncompo = lplate[ip[ipp].idm[im]].compother(ipp);
      break;
    case creep_damage:
      ncompo += givencompother (ipp, im+1); // creep
      ncompo += givencompother (ipp, im+2); // damage
      break;
    case shrinkagemat:
      ncompo += 1 + 2*ncompstr;
      ncompo += givencompother (ipp, im+1);
      break;
    case effective_stress:
      ncompo += 1; // store initial value of effective pore pressure
      ncompo += givencompother (ipp, im+1);
      break;
    case lemaitr:
      ncompo=0;
      break;
    case glasgowmechmat:{	
      ncompo=0;
      break;
    }
    case viscoplasticity:{
      ncompo = givencompeqother (ipp, im+1);  //  viscous material model
      ncompo += givencompeqother (ipp, im+2); //  plasticity material model
      break;
    }
    case viscoelasticity:{
      ncompo=0;
      break;
    }
    case time_switchmat:
      ncompo = tswmat[ip[ipp].idm[im]].givencompother (ipp,im);
      break;
    
    default:
      print_err("unknown material type is required", __FILE__, __LINE__, __func__);
    }

  return ncompo;
}



/**
  The function returns the number of components of ipp's eqother array.
  Material type is given by the im parameter, which means index of material
  in the tm array of the given integration point ipp.
   
  @param[in] ipp - integration point pointer
  @param[in] im - index of material
   
  @return The function returns total number of components of eqother array.  

  Created by Tomas Koudelka, 5.11.2003
*/
long mechmat::givencompeqother (long ipp,long im)
{
  long ncompo = 0;
  long ncompstr = ip[ipp].ncompstr;
  long n;
  
  switch (ip[ipp].tm[im]){
    case elisomat:
      ncompo = 2;//provisionally
      break;
    case elortomat:
      break;
    case elgmat3d:
    case elgmat2d:
      break;
    case elasttimemat:
      ncompo = ncompstr*2;
      break;
      //case geoelmat:
      //ncompo=ncompstr+1;
      //break;
    case layerplate:
      ncompo = lplate[ip[ipp].idm[im]].compeqother(ipp);
      break;
    case winklerpasternak:
      break;
    case graphm:
      break;
    case varelisomat:
      break;
    case elisopdmat:
      ncompo = 2;//provisionally
      break;
    case geoelast:
      ncompo  = 1;
      ncompo += givencompeqother(ipp, im+1); 
      break;
    case simplas1d:
      ncompo=3;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case jflow:
      ncompo=ncompstr+2;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case microplaneM4:
      n=mpM4[ip[ipp].idm[im]].numberOfMicroplanes;
      ncompo=3*n+1+ncompstr;
      break;
    case microsimp:
      n=mpSIM[ip[ipp].idm[im]].numberOfMicroplanes;
      ncompo=n+1+6;
      break;
    case microfibro:
      n=mpfib[ip[ipp].idm[im]].numberOfMicroplanes;
      ncompo=5*n+1+6;
      break;
    case mohrcoul:
      ncompo=ncompstr+1+2;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case mohrcoulparab:
      ncompo=ncompstr+1;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case boermaterial:
      ncompo=ncompstr+1;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case doubledrprager:
      ncompo=ncompstr+4;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case druckerprager:
      ncompo=ncompstr+4;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case druckerprager2:
      ncompo=ncompstr+4;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case chenplast:
      ncompo=ncompstr+3;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case modcamclaymat:
      ncompo=ncompstr+11;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case modcamclaycoupmat:
      ncompo=ncompstr+15+ncompstr*ncompstr;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case bbmcoupmat:
      ncompo=ncompstr+15+ncompstr*ncompstr+ncompstr;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case doublestructuremat:
      ncompo=ncompstr+15+ncompstr*ncompstr+1+ncompstr;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case homomatm:
      break;
    case hypoplastmat:
      ncompo=2*ncompstr+hypopl[ip[ipp].idm[im]].nstatv+1;
      break;
    case hypoplastusatthermat:
      ncompo=hypoplustherm[ip[ipp].idm[im]].givencompeqother(ipp);
      break;
    case shefpl:
      ncompo=ncompstr+2;
      break;
    case glasgowmechmat:
      ncompo=ncompstr*4+5;
      break;
    case glasgowdamage:
      ncompo=2;
      break;
    case creep_damage:
      ncompo += givencompeqother (ipp, im+1); // creep
      ncompo += givencompeqother (ipp, im+2); // damage
      break;
    case shrinkagemat:
      ncompo += 1 + 2*ncompstr;
      ncompo += givencompeqother (ipp, im+1);
      break;
    case time_switchmat:
      ncompo = tswmat[ip[ipp].idm[im]].givencompeqother (ipp,im);
      break;
    case effective_stress:
      ncompo += 1; // store initial value of effective pore pressure
      ncompo += givencompeqother (ipp, im+1);
      break;
    case simvisc:
      ncompo=3*ncompstr;
      break;
    case isovisc:
      ncompo=0;
      break;
    case lemaitr:
      ncompo=3*ncompstr+2;
      break;
    case scaldamage:
      ncompo=3;
      break;
    case fixortodamage:
      ncompo=18;
      break;  
    case scaldamagecc:
      ncompo=2+4;
      break;
    case ortodamage:
      ncompo=19;
      break;
    case ortodamage2:
      ncompo=19;
      break;
    case ortodamagerot:
      ncompo=23;
      break;
    case anisodamage:
//      ncompo=10;
      ncompo=26;
      break;
    case anisodamagerot:
//      ncompo=10;
      ncompo=22;
      break;
    case cebfipcontmat:
    case contmat:
      ncompo = ncompstr;
      break;
    case damplifmat:
      ncompo = 4+ncompstr-1;
      if (damplifm[ip[ipp].idm[im]].coref == yes)
        ncompo += 1;
      break;
    case plastifmat:
      ncompo = 1+ncompstr;
      break;
    case creeprs:
    case creepb3:
    case creepdpl:{
	ncompo=creep_ncompo(ipp,im);
	break;
    }
    case creepbaz:
      //ncompo=ncompstr*9+3+ncompstr;
      n=crbaz[ip[ipp].idm[im]].numberOfCreepb();
      ncompo=ncompstr*(n+3)+3+1;
      break;
    case creepeffym:
      ncompo=1+ncompstr;
      break;
    case consolidation:
      n=csol[ip[ipp].idm[im]].numberOfConsol();
      if(ncompstr==6) ncompo=ncompstr+2*n+1;
      if(ncompstr==5) ncompo=ncompstr+3*n+1;
      if(ncompstr==3) ncompo=ncompstr+3*n+1;
      //if(ncompstr==6) ncompo=ncompstr+2*7+1;
      //if(ncompstr==3) ncompo=ncompstr+3*7+1;
      break;
      
    case nonlocplastmat:
      ncompo = givencompeqother (ipp, im+1);
      break;
    case nonlocdamgmat:
      ncompo = givencompeqother (ipp, im+1);
      break;
    case nonlocalmicroM4:
      n=mpM4[ip[ipp].idm[im]].numberOfMicroplanes;
      ncompo=3*n+1+ncompstr+ncompstr;
      break;
    case damage_plasticity:
      ncompo = givencompeqother (ipp, im+1);
      ncompo += givencompeqother (ipp, im+2);
      break;
    case viscoplasticity:
      ncompo = givencompeqother (ipp, im+1);  //  viscous material model
      ncompo += givencompeqother (ipp, im+2); //  plasticity material model
      break;
    case viscoelasticity:
      ncompo = givencompeqother (ipp, im+1);  //  viscous material model
      break;
    case therisodilat:{
      break;
    }
    case thervolisodilat:{
      break;
    }
    case relaxationeuro:{
      break;
    }
    case lenjonespot:{
      ncompo = 0;
      break;
    }
    case cusatismat:{
      ncompo = 30;
      break;
    }
    default:{
      print_err("unknown material type is required", __FILE__, __LINE__, __func__);
    }
    }
  
  if (((ip[ipp].hmt & 1) && (im < ip[ipp].nm-1)) && (im == 0))
    // thermal dilatancy is present and ncompeqother for different 
    // material type (than thermal dilatancy material) is required
    // and index of this material is 0 (i.e it is main control material)
    {
      switch (ip[ipp].tm[ip[ipp].nm-1])
	{
	case therisodilat:
	  break;
	case thervolisodilat:
	  break;
	default:{
	  print_err("unknown thermal material type is required", __FILE__, __LINE__, __func__);
	}
	}
    }
  return ncompo;
}



/**
  Function stores components of array nonloc.
   
  @param[in] ipp - integration point pointer
  @param[in] fi - first index
  @param[in] ncomp - number of required components
  @param[in] comp - %vector containing components
   
  @return The function does not return anything.

  Created by TKo, 29.11.2013
*/
void mechmat::storenonloc (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].nonloc[i]=comp[j];
    j++;
  }
}



/**
  Function restores part of array nonloc to %vector comp.
   
  @param[in]  ipp   - integration point pointer
  @param[in]  fi    - first index
  @param[in]  ncomp - number of required components
  @param[out] comp  - array containing components

  @return The function returns part of nonloc array in the parameter comp.
   
  Created by TKo, 29.11.2013
*/
void mechmat::givenonloc (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    comp[j]=ip[ipp].nonloc[i];
    j++;
  }
}






// ******************************************************************
// ******************************************************************
//
// Functions for retrieving of mechanical quantities from int. points
//
// ******************************************************************
// ******************************************************************



/**
  Function assembles required quantity stored at integration point.
   
  @param[in]  iq   - id of required quantity
  @param[in]  lcid - load case id
  @param[in]  ipp  - number of integration point
  @param[in]  fi   - first index
  @param[out] ipv  - array of required values from integration point
   
  @return The function returns required %vector of quantitiy in the parameter ipv.
 
  Created by JK, TKo, 29.4.2008
*/
void mechmat::givequantity (integratedquant iq,long lcid,long ipp,long fi,vector &ipv)
{
  switch (iq){
  case locstress:{
    //  stress reading from integration point
    givestress (lcid,ipp,fi,ipv);
    break;
  }
  case nonlocstress:{
    //  stress reading from integration point
    givestress (lcid,ipp,fi,ipv);
    break;
  }
  case stressincr:{
    //  stress increments from integration point
    givestressincr (lcid,ipp,0,0,fi,ipv);
    break;
  }
  case eigstress:{
    //  eigenstress reading from integration point
    giveeigstress (ipp,fi,ipv);
    break;
  }
  case penergydens:{
    long ncompstr = ip[ipp].ncompstr;
    vector eps(ASTCKVEC(ncompstr));
    vector sig(ASTCKVEC(ncompstr));
    givestress (lcid, ipp, fi, sig);
    givestrain (lcid, ipp, fi, eps);
    scprd(sig, eps, ipv(0));
    break;
  }
  default:{
    print_err("unknown type %d of quantity is required on element %ld, ipp=%ld",__FILE__,__LINE__,__func__, int(iq), Mm->elip[ipp]+1, ipp);
  }
  }
}



/**
  Function returns actual value of the Young modulus for given material and
  in the required integration point.
  
  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.
   
  @return Function returns value of Young moduls.

  Created by Tomas Koudelka, 8.8.2005
*/
double mechmat::give_actual_ym(long ipp, long im, long ido)
{
  long ncompo;
  double e = 0.0;

  switch (ip[ipp].tm[im])
  {
    case elisomat:
      e = eliso[ip[ipp].idm[im]].e;
      break;
    case elisopdmat:
      e = elisopd[ip[ipp].idm[im]].give_actual_ym(ipp, ido);
      break;
    case modcamclaymat:
      e = cclay[ip[ipp].idm[im]].give_actual_ym(ipp, im, ido);
      break;
    case modcamclaycoupmat:{
      e = cclayc[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      break;
    }
    case bbmcoupmat:{
      e = bbm[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      break;
    }
    case doublestructuremat:
      {
	e = dsm[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
	break;
      }
    case creepeffym:
      e = creffym[ip[ipp].idm[im]].give_actual_ym();
      break;
    case elasttimemat:
      e = eltimemat[ip[ipp].idm[im]].actual_modulus(ipp);
      break;
    case winklerpasternak:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case shefpl:
    case glasgowmechmat:
    case glasgowdamage:
    case simvisc:
    case lemaitr:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepbaz:
    case consolidation:
      ncompo = givencompeqother(ipp, im); // this material should be last one before elastic
      ncompo -= givencompeqother(ipp, im+1);
      e = give_actual_ym(ipp, im+1, ido+ncompo);
      break;
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e= give_actual_ym(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case viscoplasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e = give_actual_ym(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
      ncompo = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+1);
      e = give_actual_ym(ipp, im+1, ido+ncompo); // these artificial materials sould not have internal variables, i.e. ncompo=0 -> no ido index shift
      break;
    case creep_damage:
      ncompo = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+1); // number of other components of creep material model
      ncompo -= givencompeqother(ipp, im+2); // number of other components of damage material model
      e = give_actual_ym(ipp, im+1, ido+ncompo); // these artificial materials sould not have internal variables, i.e. ncompo=0 -> no ido index shift
      break;
    case shrinkagemat:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the shrinkage material
      e = give_actual_ym (ipp, im+1, ido+ncompo);
      break;
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material
      e = give_actual_ym(ipp, im+1, ido+ncompo);
      break;
    case creepb3:
    case creeprs:
    case creepdpl:
      e = creep_give_actual_ym(ipp, im, ido);
      break;
    case time_switchmat:
      e = tswmat[ip[ipp].idm[im]].give_actual_ym(ipp, im, ido);
      break;
    }
    default:
      print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }

  return e;
}



/**
  Function returns initial value of the Young modulus for given material and
  in the required integration point.
   
  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.

  @return Function returns value of initial Young modulus.

  Created by Tomas Koudelka, 8.2008
*/
double mechmat::give_initial_ym(long ipp, long im, long ido)
{
  long ncompo;
  double e = 0.0;

  switch (ip[ipp].tm[im])
    {
    case elisomat:
      e = eliso[ip[ipp].idm[im]].e;
      break;
    case varelisomat:
      e = veliso[ip[ipp].idm[im]].e;
      break;
    case elisopdmat:
      e = elisopd[ip[ipp].idm[im]].e;
      break;
    case creepeffym:
      e = creffym[ip[ipp].idm[im]].give_initial_ym();
      break;
    case elasttimemat:
      e = eltimemat[ip[ipp].idm[im]].initial_modulus();
      break;
    case winklerpasternak:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
      break;
    case modcamclaymat:{
      e = cclay[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      break;
    }
    case modcamclaycoupmat:{
      e = cclayc[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      break;
    }
    case bbmcoupmat:{
      e = bbm[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      break;
    }
    case doublestructuremat:
      {
	e = dsm[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
	break;
      }
    case shefpl:
    case glasgowmechmat:
    case glasgowdamage:
    case simvisc:
    case lemaitr:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepbaz:
    case consolidation:
      ncompo = givencompeqother(ipp, im); // this material should be last one before elastic
      ncompo -= givencompeqother(ipp, im+1);
      e = give_initial_ym(ipp, im+1, ido+ncompo);
      break;
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e= give_initial_ym(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case viscoplasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e = give_initial_ym(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      e = give_initial_ym(ipp, im+1, ido+ncompo); // these artificial materials should have no internal variables (ncompo=0) -> no ido index shift
      break;
    case creep_damage:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      ncompo -= givencompeqother(ipp, im+2);
      e = give_initial_ym(ipp, im+1, ido+ncompo); // these artificial materials should have no internal variables (ncompo=0) -> no ido index shift
      break;
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material 
      e = give_initial_ym(ipp, im+1, ido+ncompo);
      break;
    case shrinkagemat:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the shrinkage material
      e = give_initial_ym (ipp, im+1, ido+ncompo);
      break;
    case creepdpl:
    case creepb3:
    case creeprs:
      e = creep_compute_inital_ym(ipp, im, ido);
      break;
    case time_switchmat:
      e = tswmat[ip[ipp].idm[im]].give_initial_ym(ipp, im, ido);
      break;
    default:
      print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }
  return e;
}



/**
  Function returns actual value of the Poisson's ratio for given material and
  in the required integration point.

  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.
   
  @return Function returns value of Poisson's ratio.

  Created by Tomas Koudelka, 7.2008
*/
double mechmat::give_actual_nu(long ipp, long im, long ido)
{
  long ncompo;
  double nu = 0.0;

  switch (ip[ipp].tm[im])
  {
    case elisomat:
      nu = eliso[ip[ipp].idm[im]].nu;
      break;
    case elisopdmat:
      nu = elisopd[ip[ipp].idm[im]].nu;
      break;
    case modcamclaymat:
      nu = cclay[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    case modcamclaycoupmat:{
      nu = cclayc[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    }
    case bbmcoupmat:{
      nu = bbm[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    }
    case doublestructuremat:{
      nu = dsm[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    }
    case creepeffym:
      nu = creffym[ip[ipp].idm[im]].give_actual_nu();
      break;
    case elasttimemat:
      nu = eltimemat[ip[ipp].idm[im]].nu;
      break;
    case winklerpasternak:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case shefpl:
    case glasgowmechmat:
    case glasgowdamage:
    case simvisc:
    case lemaitr:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepdpl:
    case creepbaz:
    case creepb3:
    case creeprs:
    case consolidation:
      ncompo  = givencompeqother(ipp, im); // this material should be last one before elastic
      ncompo -= givencompeqother(ipp, im+1); 
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case viscoplasticity:
      ncompo  = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case creep_damage:
      ncompo  = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case shrinkagemat:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the shrinkage stress material
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case effective_stress:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      nu = give_actual_nu(ipp, im+1, ido+ncompo); // these artificial materials should have no internal variables, i.e. ncompo=0 -> no ido index shift
      break;
    case time_switchmat:
      nu = tswmat[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    default:
      print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }
  return nu;
}



/**
  Function returns initial value of the Poisson's ratio for given material and
  in the required integration point.

  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.
   
  @return Function returns value of Poisson's ratio.

  Created by Tomas Koudelka, 7.2008
*/
double mechmat::give_initial_nu(long ipp, long im, long ido)
{
  long ncompo;
  double nu = 0.0;

  switch (ip[ipp].tm[im])
  {
    case elisomat:
      nu = eliso[ip[ipp].idm[im]].nu;
      break;
    case varelisomat:
      nu = veliso[ip[ipp].idm[im]].nu;
      break;
    case elisopdmat:
      nu = elisopd[ip[ipp].idm[im]].nu;
      break;
    case creepeffym:
      nu = creffym[ip[ipp].idm[im]].nu;
      break;
    case winklerpasternak:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case modcamclaymat:
    case modcamclaycoupmat:
    case bbmcoupmat:
    case doublestructuremat:
    case shefpl:
    case glasgowmechmat:
    case glasgowdamage:
    case simvisc:
    case lemaitr:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepdpl:
    case creepbaz:
    case creepb3:
    case creeprs:
    case consolidation:
      ncompo  = givencompeqother(ipp, im); // this material should be last one before elastic
      ncompo -= givencompeqother(ipp, im+1); 
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case viscoplasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case creep_damage:
      ncompo  = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      nu = give_actual_nu(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      break;
    case shrinkagemat:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the shrinkage material
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material
      nu = give_actual_nu(ipp, im+1, ido+ncompo);
      break;
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); 
      nu = give_actual_nu(ipp, im+1, ido+ncompo); // these artificial materials should have no internal variables, i.e. ncompo=0 -> no ido index shift
      break;
    case time_switchmat:
      nu = tswmat[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      break;
    default:
      print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }
  return nu;
}



/**
  Function returns actual value of the tensile strenght for given material and
  in the required integration point.

  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.

  @return Function returns tensile strength.
   
  Created by Tomas Koudelka, 8.8.2005
*/
double mechmat::give_actual_ft(long ipp, long im, long ido)
{
  long ncompo;
  double ft = 0.0;

  switch (ip[ipp].tm[im]){
    case scaldamage:
      ft = scdam[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case ortodamage:
      ft = ortdam[ip[ipp].idm[im]].give_actual_ft(ipp);
      break;
    case ortodamage2:
      ft = ortdam2[ip[ipp].idm[im]].give_actual_ft(ipp);
      break;
    case ortodamagerot:
      ft = ortdamrot[ip[ipp].idm[im]].give_actual_ft(ipp);
      break;
    case anisodamage:
      ft = anidam[ip[ipp].idm[im]].give_actual_ft(ipp);
      break;
    case anisodamagerot:
      ft = anidamrot[ip[ipp].idm[im]].give_actual_ft(ipp);
      break;
    case nonlocdamgmat:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); 
      ft = give_actual_ft(ipp,im+1,ido+ncompo);
      break;
    case damage_plasticity:
      ft = dampl[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case creep_damage:
      ft = crdam[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case shrinkagemat:
      ft = shmat[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); 
      ft = give_actual_ft(ipp,im+1,ido+ncompo);
      break;
    case creepdpl:
    case creepb3:
    case creeprs:
      ft = creep_give_actual_ft(ipp, im, ido);
      break;
    case creepeffym:
      ft = creffym[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case time_switchmat:
      ft = tswmat[ip[ipp].idm[im]].give_actual_ft(ipp, im, ido);
      break;
    case damplifmat:
      ft = damplifm[ip[ipp].idm[im]].give_actual_ft();
      break;
    default:
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  return ft;
}



/**
  Function returns actual value of the compressive strenght for given material and
  in the required integration point.
  
  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.

  @return Function returns tensile strength.
   
  Created by Tomas Koudelka, 8.8.2005
*/
double mechmat::give_actual_fc(long ipp, long im, long ido)
{
  long ncompo;
  double fc = 0.0;

  switch (ip[ipp].tm[im])
  {
    case ortodamage:
      fc = ortdam[ip[ipp].idm[im]].give_actual_fc(ipp);
      break;
    case ortodamage2:
      fc = ortdam2[ip[ipp].idm[im]].give_actual_fc(ipp);
      break;
    case ortodamagerot:
      fc = ortdamrot[ip[ipp].idm[im]].give_actual_fc(ipp);
      break;
    case anisodamage:
      fc = anidam[ip[ipp].idm[im]].give_actual_fc(ipp);
      break;
    case anisodamagerot:
      fc = anidamrot[ip[ipp].idm[im]].give_actual_fc(ipp);
      break;
    case nonlocdamgmat:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); 
      fc = give_actual_fc(ipp,im+1,ido+ncompo);
      break;
    case creep_damage:
      fc = crdam[ip[ipp].idm[im]].give_actual_fc(ipp, im, ido);
      break;
    case shrinkagemat:
      fc = shmat[ip[ipp].idm[im]].give_actual_fc(ipp, im, ido);
      break;
    default:
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  return fc;
}



/**
  The function returns a characteristics of the required quantity.
  The return value is the natural mathematic representation of the quantity - 
  scalar, vector or the second order tensor. Additionally, it returns the maximum 
  number of components used for the qunatity storage in the argument ncomp.

  If the given quantity is the second order tensor, the function returns quantity storage 
  notation in the argument tn and tensor type identifier in the argument tti.

  @param[in] mq - required quantity type.
  @param[in] refmq - referenced mechanical quantity, used only for mq == tensdeviator
  @param[out] ncomp - the maximum number of quantity components.
  @param[out] tn - tensor storarge notation for quantities representded as tensors otherwise 0
  @param[out] tti - tensor type identifier (stress/strain/other) for quantities representded as tensors otherwise -1

  @retval scalq if the given quantity is represented as a scalar.
  @retval vectq if the given quantity is represented as a vector.
  @retval tensq if the given quantity is represented as the second order tensor.

  Created by TKo, 09.2023
*/
quantrep mechmat::give_quant_rep(mechquant mq, mechquant refmq, long &ncomp, tensqnot &tn, strastre &tti)
{
  quantrep ret = undefq;
  tn = undeftn;
  tti = undef_strastre;
  ncomp = 1;
  
  switch (mq){
    case time_q:          // time (constant scalar quantity at whole domain)
    case step_id_q:       // step id (constant scalar quantity at whole domain)
    case load_fact_q:     // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
    case eigval_q:        // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
    case eps_x:           // strain components as scalars
    case eps_y:
    case eps_z:
    case eps_yz:
    case eps_xz:
    case eps_xy:
    case epsp_x:          // components of plastic strains considered as scalars
    case epsp_y:
    case epsp_z:
    case gammap_yz:
    case gammap_xz:
    case gammap_xy:
    case sig_x:            // stress components as scalars
    case sig_y:
    case sig_z:
    case tau_yz:
    case tau_xz:
    case tau_xy:
    case first_inv:       // A11 + A22 + A33
    case second_inv:      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2
    case third_inv:       // det|A|
    case tensor_norm:     // ||A|| = sqrt(a_ij*a_ij)
    case strain_vol:      // eps_v = eps_x + eps_y + eps_z
    case mean_stress:     // sig_m = (sig_x + sig_y + sig_z)/3
    case j2inv:           // negative value of the second invariant of stress deviator, i.e. J2 = 1/2 s_ij s_ij
    case von_mises_stress:// sig_eff = sqrt(3*J2)                
    case cons_param:      // consistency parameter gamma
    case damage_scal:     // scalar damage (omega)
    case damaget_scal:    // scalar damage (omega_t) in tension
    case damagec_scal:    // scalar damage (omega_c) in compression
      ret = scalq;
      break;
    case displ_q:        // displacement vector
    case react_q:        // reactions vector at nodes                
    case load_vect_q:    // load vector at nodes
    case int_force_q:    // internal force vector
    case resid_vect_q:    // residual vector
      ret = vectq;
      ncomp = Mt->give_maxndofn();
      break;
    case nonmech_q:       // nonmechanical quantities at ip
      ret = vectq;
      ncomp = Mm->nnmq;
      break;
    case other_q:        // other array values as a vector
      ret = vectq;
      ncomp = std::max(max_ncompothern, max_ncompothere);
      break;
    case strain_q:       // strain tensor
    case tempr_strain_q: // temperature strain tensor
    case eig_strain_q:   // eigenstrain tensor
    case macrostrain_q:  // macrostrain (constant tensor quantity at whole domain)
    case strain_deviator: // s_ij = sig_ij - delta_ij*sig_m
    case strain_pl:       // plastic strain tensor
      ret = tensq;
      tn = voigtred;
      tti = strain;
      ncomp = std::max(max_ncompstrn, max_ncompstre);
      break;
    case stress_q:       // stress tensor
    case eig_stress_q:   // eigenstress tensor
    case macrostress_q:  // macrostress (constant tensor quantity at whole domain)
    case stress_deviator: // e_ij = eps_ij - delta_ij*eps_v/3
      ret = tensq;
      tn = voigtred;
      tti = stress;
      ncomp = std::max(max_ncompstrn, max_ncompstre);
      break;      
    case tensdeviator:    // D_ij = A_ij - delta_ij*A_kk/3 
      ret = give_quant_rep(mq, refmq, ncomp, tn, tti);
      break;
    case damage_tens:     // damage tensor (Omega)
    case damaget_tens:    // damage tensor (Omega_t) in tension
    case damagec_tens:    // damage tensor (Omega_c) in compression
      ret = tensq;
      tn = fullmtx;
      tti = other;
      ncomp = 3*3;
      break;
    default:
      print_err("unknown type of quantitity (%d) is required\n", __FILE__, __LINE__, __func__, int(mq));
      abort();
  }
  return ret;
}



/**
  The function returns a required scalar qunatity in the given integration point. 
  The value is returned in the argument qv.

  @param[in] ipp - integration point pointer
  @param[in] mqn - quantity name
  @param[in] reftensq - referenced tensor quantity, it is used only if the mqn 
                        is an tensor derived value (invariant, norm,...), otherwise it is ignored
  @param[in] lcid - load case id or eigenvalue id.
  @param[in] selrefcomp - selection of referenced tensor components, it is used only if the mqn 
                          is an tensor derived value (invariant, norm,...), otherwise it is ignored
  @param[in]  ref_strastre - identfier of of stress/strain/other type of tensors given by array, i.e. 
                             reftensq == (other_q || nonmech_q), otherwise it is ignored
  @param[out] qv - required quantity value.

  @return The function returns the quantity value in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechmat::give_quant(long ipp, mechquant mqn, mechquant reftensq, long lcid, sel &selrefcomp,
                         strastre ref_strastre, double &qv)
{
  matrix auxm;
  vector eps, epst, sig;
  qv = 0.0;
  switch (mqn){
    case time_q:          // time (constant scalar quantity at whole domain)
      qv = Mp->time;
      break;
    case step_id_q:       // step id (constant scalar quantity at whole domain)
      qv = Mp->istep;
      break;
    case load_fact_q:     // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
      qv = Mp->lambda;
      break;
    case eigval_q:        // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
      //qv = Lsrs->eigv[lcid];
      break;
    case first_inv:       // A11 + A22 + A33
      reallocm(RSTCKMAT(3, 3, auxm));
      give_quant(ipp, reftensq, nomech_q, lcid, selrefcomp, ref_strastre, auxm);
      qv = first_invar(auxm);;
      break;
    case second_inv:      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2
      reallocm(RSTCKMAT(3, 3, auxm));
      give_quant(ipp, reftensq, nomech_q, lcid, selrefcomp, ref_strastre, auxm);
      qv = second_invar(auxm);;
      break;
    case third_inv:       // det|A|
      reallocm(RSTCKMAT(3, 3, auxm));
      give_quant(ipp, reftensq, nomech_q, lcid, selrefcomp, ref_strastre, auxm);
      qv = third_invar(auxm);;
      break;
    case tensor_norm:     // ||A|| : sqrt(a_ij*a_ij)
      reallocm(RSTCKMAT(3, 3, auxm));
      give_quant(ipp, reftensq, nomech_q, lcid, selrefcomp, ref_strastre, auxm);
      qv = norm(auxm);;
      break;
    case strain_vol:      // eps_v : eps_x + eps_y + eps_z
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = first_invar(eps);
      break;
    case mean_stress:     // sig_m : (sig_x + sig_y + sig_z)/3
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = (1.0/3.0)*first_invar(sig);
      break;
    case j2inv:           // negative value of the second invariant of stress deviator, i.e. J2 : 1/2 s_ij s_ij
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = j2_stress_invar(sig);      
      break;
    case von_mises_stress:// sig_eff : sqrt(3*J2)                
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = 3.0*j2_stress_invar(sig);
      qv = sqrt(qv);
      break;
    case cons_param:      // consistency parameter gamma
      qv = give_consparam(ipp);
      break;
    case damage_scal:     // scalar damage (omega)
      qv = give_dampar(ipp);
      break;
    //case damaget_scal:    // scalar damage (omega_t) in tension
    //case damagec_scal:    // scalar damage (omega_c) in compression      
    case eps_x:           // strain components as scalars
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(0);
      break;
    case eps_y:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(1);
      break;
    case eps_z:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(2);
      break;
    case gamma_yz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(3);
      break;
    case gamma_xz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(4);
      break;
    case gamma_xy:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = eps(5);
      break;
    case eps_yz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = 0.5*eps(3);
      break;
    case eps_xz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = 0.5*eps(4);
      break;
    case eps_xy:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, ip[ipp].strain, ip[ipp].ssst);
      qv = 0.5*eps(5);
      break;
    case sig_x:           // stress components as scalars
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(0);
      break;
    case sig_y:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(1);
      break;
    case sig_z:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(2);
      break;
    case tau_yz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(3);
      break;
    case tau_xz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(4);
      break;
    case tau_xy:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(sig, ip[ipp].stress, ip[ipp].ssst);
      qv = sig(5);
      break;
    case epsp_x:          // components of plastic strains considered as scalars
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(0);
      break;
    case epsp_y:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(1);
      break;
    case epsp_z:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(2);
      break;
    case gammap_yz:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(3);
      break;
    case gammap_xz:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(4);
      break;
    case gammap_xy:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = epst(5);
      break;
    case epsp_yz:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = 0.5*epst(3);
      break;
    case epsp_xz:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = 0.5*epst(4);
      break;
    case epsp_xy:
      reallocv(RSTCKVEC(ip[ipp].ncompstr, eps));
      reallocv(RSTCKVEC(6, epst));
      giveirrstrains(ipp, 0, 0, eps);
      give_full_vector(epst, eps, ip[ipp].ssst);
      qv = 0.5*epst(5);
      break;
    default:
      print_err("unknown scalar quantity type (%d) is required on element %ld, ipp=%ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), elip[ipp]+1, ipp);
      abort();
  }
}



/**
  The function returns required components of a %vector qunatity in the given integration point. 
  The quantity components are returned in the argument qv.

  @param[in] ipp - integration point pointer
  @param[in] mqn - quantity name
  @param[in] lcid - load case id or eigenvalue id.
  @param[in] selcomp - selection of the %vector components
  @param[out] qv - required quantity value.

  @return The function returns the quantity components in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechmat::give_quant(long ipp, mechquant mqn, long lcid, sel &selcomp, vector &qv)
{
  long i, id;

  switch (mqn){
    //case displ_q:        // displacement vector, it will be handled in new materialpoints concept, now displacements are not defined on integration points

    // following options have no meaning for the integration points
    //case react_q:        // reactions vector at nodes                
    //case load_vect_q:    // load vector at nodes
    //case int_force_q:    // internal force vector
    //case resid_vect_q:   // residual vector
    
    case nonmech_q:       // nonmechanical quantities at ip
      id = 0;
      for(i=0; i < nnmq; i++){
        if (selcomp.presence_id(i)){
          qv(id) = nonmechq[i][ipp];
          id++;
        }
      }      
      break;
    case other_q:        // other array values as a vector
      id = 0;
      for(i=0; i < ip[ipp].ncompother; i++){
        if (selcomp.presence_id(i)){
          qv(id) = ip[ipp].other[i];
          id++;
        }
      }      
      break;
    default:
      print_err("unknown vector quantity type (%d) is required on element %ld, ipp=%ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), elip[ipp]+1, ipp);
      abort();
  }
}



/**
  The function returns a required tensor qunatity in the given integration point.
  The quantity components are returned in the argument qv.

  @param[in]  ipp - integration point pointer
  @param[in]  mqn - required quantity name
  @param[in]  reftensq - referenced quantity handled as an ordinary tensor
  @param[in]  lcid - load case id or eigenvalue id.
  @param[in]  selcomp - selection of array components for mqn == (other_q || nonmech_q)
                        or mqn == tensdeviator reftensq == (other_q || nonmech_q), 
                        otherwise it is ignored
  @param[in]  array_strastre - identfier of of stress/strain/other type of tensors given by array, i.e. 
                               mqn == (other_q || nonmech_q) or for referenced tensor reftensq == (other_q || nonmech_q),
                               otherwise it is ignored
  @param[out] qv - required quantity value.

  @return The function returns the quantity components in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechmat::give_quant(long ipp, mechquant mqn, mechquant reftensq, long lcid, sel &selcomp,
                         strastre array_strastre, matrix &qv)
{
  vector eps, epst, eps0, epsp, meps, sig, sig0, msig, other_arr;
  double *epsptr, *sigptr;
  long ncompstr, ncsel, nc;
  double inv;
  
  switch (mqn){
    case strain_q:        // strain tensor
      ncompstr  = ip[ipp].ncompstr;
      epsptr = ip[ipp].strain+lcid*ncompstr;
      makerefv(eps, epsptr, ncompstr);
      vector_tensor(eps, qv, ip[ipp].ssst, strain);
      break;
    case stress_q:        // stress tensor
      ncompstr  = ip[ipp].ncompstr;
      sigptr = ip[ipp].stress+lcid*ncompstr;
      makerefv(sig, sigptr, ncompstr);
      vector_tensor(sig, qv, ip[ipp].ssst, stress);
      break;
    case tempr_strain_q:  // temperature strain tensor
      makerefv(epst, tempstrains[ipp], ip[ipp].ncompstr);
      vector_tensor(epst, qv, ip[ipp].ssst, strain);
      break;
    case eig_strain_q:    // eigenstrain tensor
      makerefv(eps0, eigstrains[ipp], ip[ipp].ncompstr);
      vector_tensor(eps0, qv, ip[ipp].ssst, strain);
      break;
    case eig_stress_q:    // eigenstress tensor
      makerefv(sig0, eigstresses[ipp], ip[ipp].ncompstr);
      vector_tensor(sig, qv, ip[ipp].ssst, stress);
      break;
    case macrostrain_q:   // macrostrain (constant tensor quantity at whole domain)
      reallocv(RSTCKVEC(Mt->max_ncompstr, meps));
      macrostrains(lcid, meps);
      vector_tensor(meps, qv, ip[ipp].ssst, strain);
      break;
    case macrostress_q:   // macrostress (constant tensor quantity at whole domain)
      reallocv(RSTCKVEC(Mt->max_ncompstr, msig));
      macrostresses(lcid, msig);
      vector_tensor(msig, qv, ip[ipp].ssst, stress);
      break;
    case tensdeviator:    // D_ij : A_ij - delta_ij*A_kk/3 
      give_quant(ipp, reftensq, nomech_q, lcid, selcomp, array_strastre, qv);      
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    case strain_deviator: // e_ij : eps_ij - delta_ij*eps_v/3
      makerefv(eps, ip[ipp].strain, ip[ipp].ncompstr);
      vector_tensor(eps, qv, ip[ipp].ssst, strain);
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    case stress_deviator: // s_ij : sig_ij - delta_ij*sig_m
      makerefv(sig, ip[ipp].stress, ip[ipp].ncompstr);
      vector_tensor(sig, qv, ip[ipp].ssst, stress);
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    case strain_pl:       // plastic strain tensor
      reallocv(RSTCKVEC(ip[ipp].ncompstr, epsp));
      giveirrstrains(ipp, 0, 0, epsp);
      vector_tensor(epsp, qv, ip[ipp].ssst, strain);
      break;
    case other_q:
      ncsel = selcomp.give_nselcomp(ip[ipp].ncompother);
      if (ncsel == 9){
        makerefv(other_arr, qv.a, ncsel);
        giveother(ipp, selcomp, other_arr);
      }
      else{
        nc = give_ncompstr(ip[ipp].ssst);
        if (nc == ncsel){
          reallocv(RSTCKVEC(ncsel, other_arr));
          giveother(ipp, selcomp, other_arr);
          vector_tensor(other_arr, qv, ip[ipp].ssst, array_strastre);
        }
        else{
          print_err("number of selected other array components %ld does not correspond\n"
                    "to the number of tensor components in Voigt's notation %ld for the %s\n"
                    "stress/strain state at ipp=%ld, eid=%ld\n", __FILE__, __LINE__, __func__,
                    ncsel, nc, strastrestate_kwdset.get_str(ip[ipp].ssst), ipp, elip[ipp]+1);
          abort();
        }
      }      
      break;
    case nonmech_q:
      ncsel = selcomp.give_nselcomp(nnmq);
      if (ncsel == 9){
        makerefv(other_arr, qv.a, ncsel);
        give_nonmechqcomp(ipp, selcomp, other_arr);
      }
      else{
        nc = give_ncompstr(ip[ipp].ssst);
        if (nc == ncsel){
          reallocv(RSTCKVEC(ncsel, other_arr));
          give_nonmechqcomp(ipp, selcomp, other_arr);
          vector_tensor(other_arr, qv, ip[ipp].ssst, array_strastre);
        }
        else{
          print_err("number of selected nonmechq array components %ld does not correspond\n"
                    "to the number of tensor components in Voigt's notation %ld for the %s\n"
                    "stress/strain state at ipp=%ld, eid=%ld\n", __FILE__, __LINE__, __func__,
                    ncsel, nc, strastrestate_kwdset.get_str(ip[ipp].ssst), ipp, elip[ipp]+1);
          abort();
        }
      }      
      break;
    //case damage_tens:     // damage tensor (Omega)
    //case damaget_tens:    // damage tensor (Omega_t) in tension
    //case damagec_tens:    // damage tensor (Omega_c) in compression
    default:
      print_err("unknown tensor quantity type (%d) is required on element %ld, ipp=%ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), elip[ipp]+1, ipp);
      abort();
  }
}



/**
  The function returns required quantity on the given integration point.
  
  @param[in] ipp - integration point pointer
  @param[in] mq  - required quantity type
  @param[in] reftensq - referenced tensor quantity, it is used only if the mqn 
                        is an tensor derived value (invariant, norm,...)
  @param[in] lcid - load case id
  @param[in] qr  -  original character of the required quantity (scalar/vector/tensor)
  @param[in] selcomp - component selection for a vector quantity
  @param[in]  array_strastre - identfier of of stress/strain/other type of tensors given by array, i.e. 
                               mqn == (other_q || nonmech_q) or for referenced tensor reftensq == (other_q || nonmech_q),
                               otherwise it is ignored
  @param[out] sv - output argument for a scalar quantity
  @param[out] vv - output argument for a vector quantity
  @param[out] tv - output argument for a tensor quantity

  @return The function returns required quantity value(s) in one of the arguments
          sv, vv, tv according to qt argument value.

  Created by Tomas Koudelka, 09.2023
*/
void mechmat::give_quant(long ipp, mechquant mq, mechquant reftensq, long lcid, quantrep qr, sel &selcomp,
                         strastre array_strastre, double &sv, vector &vv, matrix &tv)
{
  switch(qr){
    case scalq:
      give_quant(ipp, mq, reftensq, lcid, selcomp, array_strastre, sv);
      break;
    case vectq:
      give_quant(ipp, mq, lcid, selcomp, vv);
      break;
    case tensq:
      give_quant(ipp, mq, reftensq, lcid, selcomp, array_strastre, tv);
      break;
    default:
      print_err("unknown qunatity character type (%d)\n", __FILE__, __LINE__, __func__, int(qr));
      abort();
  }
}



/**
  Function returns required mechanical quantity stored at integration point
  for TRFEL .
   
  @param[in] qt - type of required quantity
  @param[in] ipp - number of integration point
   
  @return The function returns required %vector of quantitiy in the parameter ipv.
 
  Created by Tomas Koudelka, 8.10.2013
*/
double mechmat::givemechq(nontransquant qt, long ipp)
{
  double ret=0.0;

  switch (qt){
  case saturation_deg:{
    //  saturation degree from mefel integration point
    ret = give_saturation_degree(ipp);
    break;
  }
  case precons_press:{
    //  preconsolidtaion pressure from integration point
    ret = give_preconspress(ipp);
    break;
  }
  case mean_stress_eff:{
    vector sig(ip[ipp].ncompstr);
    vector sigt(ASTCKVEC(6));
    givestress(0,ipp,sig);
    give_full_vector(sigt,sig,ip[ipp].ssst);
    // total mean stress
    ret = first_invar(sigt)/3.0;
    if (givestatusnmq(eff_pore_pressure)){
      // effective stress concept model was employed =>
      // general - given by used transport model in TRFEL or effstress model in MEFEL
      // \sigma_{eff} = \sigma_{tot} - (S_{r,w}*u_w + S_{r,a}*u_a) m^T 
      ret -= givenonmechq(eff_pore_pressure, ipp);
    }
    else{
      if (givestatusnmq(pore_pressure)){
        if (givestatusnmq(saturation_degree)){
          // effective stress concept model was employed => \sigma_{eff} = \sigma_{tot} - S_{rw}*u_w m^T (Schrefler)
          ret -= givenonmechq(saturation_degree, ipp)*givenonmechq(pore_pressure, ipp);
        }
        else{
          // effective stress concept model was employed => \sigma_{eff} = \sigma_{tot} - u_w m^T (Terzaghi)
          ret -= givenonmechq(pore_pressure, ipp);
        }
      }
    }
    break;
  }
  case virgin_porosity:{
    //  virgin porosity from integration point
    ret = give_virgporosity(ipp);
    break;
  }
  case init_porosity:{
    //  initial porosity reading from integration point
    ret = give_iniporosity(ipp);
    break;
  }
  case porosity:{
    //  actual porosity reading from integration point
    ret = give_porosity(ipp);
    break;
  }
  case void_ratio:{
    //  actual void_ratio reading from integration point
    ret = give_void_ratio(ipp);
    break;
  }
  case scal_iso_damage:{
    //  damage parameter omega
    ret = give_dampar (ipp);
    break;
  }
  case proc_zone_length:{
    //  length of process zone
    ret = give_proczonelength (ipp);
    break;
  }
  case crack_width:{
    //  crack width
    ret = give_crackwidth (ipp);
    break;
  }
  case der_saturation_deg:{
    //  derivative of saturation degree with respect to suction from mefel integration point
    ret = give_der_saturation_degree(ipp);
    break;
  }
  case der_saturation_deg_dtemp:{
    // derivative of saturation degree with respect to temperature from mefel integration point
    ret = give_der_saturation_degree_dtemp(ipp);
    break;
  }
  case strain_vol_rate: // rate of the volumetric strain
    ret = give_strain_vol_rate(ipp);
    break;
  case der_saturation_deg_depsv:{
    //  derivative of saturation degree with respect to volumetric strain from mefel integration point
    ret = give_der_saturation_degree_depsv(ipp);
    break;
  }
  case bulk_modulus:
    // bulk modulus
    ret = give_bulk_modulus(ipp);
    break;
  case mmean_stress:{
    vector sig(ip[ipp].ncompstr);
    vector sigt(ASTCKVEC(6));
    givestress(0,ipp,sig);
    give_full_vector(sigt,sig,ip[ipp].ssst);
    // total mean stress
    ret = first_invar(sigt)/3.0;
    break;
  }
  default:{
    print_err("unknown type of quantity is required",__FILE__,__LINE__,__func__);
  }
  }

  return ret;
}



/**
  Function returns saturation degree from hypo-plasticity models.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of preconsolidation pressure.
  
  Created by Tomas Krejci 14/12/2015
*/
double mechmat::give_saturation_degree(long ipp, long im, long ido)
{
  long i, ncompo;
  double s=0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case effective_stress:{
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    s = give_saturation_degree(ipp, im+1, ido+ncompo);
  break;
  }
  case modcamclaymat:{
    s = cclayc[i].give_saturation_degree(ipp,ido);
    break;
  }
  case modcamclaycoupmat:{
    s = cclayc[i].give_saturation_degree(ipp,ido);
    break;
  }
  case bbmcoupmat:{
    s = bbm[i].give_saturation_degree(ipp,ido);
    break;
  }
  case doublestructuremat:{
    s = dsm[i].give_saturation_degree(ipp,ido);
    break;
  }
  case hypoplastusatthermat:{
    s = hypoplustherm[i].give_sr(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return s;
}




/**
  Function returns the derivative of saturation degree with respect to suction.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of preconsolidation pressure.
  
  Created by Tomas Krejci 14/12/2015
*/
double mechmat::give_der_saturation_degree(long ipp, long im, long ido)
{
  long i, ncompo;
  double dsrds=0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case effective_stress:{
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    dsrds = give_der_saturation_degree(ipp, im+1, ido+ncompo);
  break;
  }
    //case modcamclaymat:{
    //dsrds = cclayc[i].give_der_saturation_degree(ipp,ido);
    //break;
    //}
    //case modcamclaycoupmat:{
    //dsrds = cclayc[i].give_der_saturation_degree(ipp,ido);
    //break;
    //}
  case hypoplastusatthermat:{
    dsrds = hypoplustherm[i].give_dsr_ds(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return dsrds;
}


/**
  Function returns the derivative of saturation degree with respect to temperature. 
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of preconsolidation pressure.
  
  Created by Tomas Koudelka 29/06/2017
*/
double mechmat::give_der_saturation_degree_dtemp(long ipp, long im, long ido)
{
  long i;
  double dsrdt=0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case hypoplastusatthermat:{
    dsrdt = hypoplustherm[i].give_dsr_dtemp(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return dsrdt;
}


/**
  Function returns the derivative of saturation degree with respect to volumetric strain.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of preconsolidation pressure.
  
  Created by Tomas Koudelka 29/06/2018
*/
double mechmat::give_der_saturation_degree_depsv(long ipp, long im, long ido)
{
  long i;
  double dsrdt=0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case hypoplastusatthermat:{
    dsrdt = hypoplustherm[i].give_dsr_depsv(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return dsrdt;
}


/**
  Function returns preconsolidation pressure especially for plasticity models.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of preconsolidation pressure.
  
  Created by Tomas Koudelka 8.10.2013
*/
double mechmat::give_preconspress(long ipp, long im, long ido)
{
  long i, ncompo;
  double pc=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case modcamclaymat:{
    pc = cclay[i].give_preconspress(ipp,ido);
    break;
  }
  case modcamclaycoupmat:{
    pc = cclayc[i].give_preconspress(ipp,ido);
    break;
  }
  case bbmcoupmat:{
    pc = bbm[i].give_preconspress(ipp,ido);
    break;
  }
  case doublestructuremat:{
    pc = dsm[i].give_preconspress(ipp,ido);
    break;
  }
  case effective_stress:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    pc = give_preconspress(ipp, im+1, ido+ncompo);
    break;
  case damage_plasticity:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+2);
    pc = give_preconspress(ipp, im+2, ido+ncompo);
    break;
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return pc;
}



/**
  Function returns virgin porosity especially for plasticity models.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of virgin porosity parameter.
  
  Created by Tomas Koudelka 8.10.2013
*/
double mechmat::give_virgporosity(long ipp, long im, long ido)
{
  long i, ncompo;
  double e_lambda1=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case modcamclaymat:{
    e_lambda1 = cclay[i].give_virgporosity(ipp,ido);
    break;
  }
  case modcamclaycoupmat:{
    e_lambda1 = cclayc[i].give_virgporosity(ipp,ido);
    break;
  }
  case bbmcoupmat:{
    e_lambda1 = bbm[i].give_virgporosity(ipp,ido);
    break;
  }
  case doublestructuremat:{
    e_lambda1 = dsm[i].give_virgporosity(ipp,ido);
    break;
  }
  case hypoplastmat:{
    e_lambda1 = hypopl[i].give_virgporosity(ipp,ido);
    break;
  }
  case effective_stress:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    e_lambda1 = give_virgporosity(ipp, im+1, ido+ncompo);
    break;
  case damage_plasticity:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+2);
    e_lambda1 = give_virgporosity(ipp, im+2, ido+ncompo);
    break;
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return e_lambda1;
}



/**
  Function returns actual value of the bulk modulus for given material and
  in the required integration point.
  
  @param[in] ipp - number of integration point
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.
   
  @return Function returns value of bulk moduls.

  Created by Tomas Koudelka, 09.2023
*/
double mechmat::give_bulk_modulus(long ipp, long im, long ido)
{
  long ncompo;
  double e, nu, ret = 0.0;
  switch (ip[ipp].tm[im])
  {
    case elisomat:
      e = eliso[ip[ipp].idm[im]].e;
      nu = eliso[ip[ipp].idm[im]].nu;
      ret = e/(1-2.0*nu)/3.0;
      break;
    case elisopdmat:
      e = elisopd[ip[ipp].idm[im]].give_actual_ym(ipp, ido);
      nu = elisopd[ip[ipp].idm[im]].nu;
      ret = e/(1-2.0*nu)/3.0;
      break;
    case modcamclaymat:
      e = cclay[ip[ipp].idm[im]].give_actual_ym(ipp, im, ido);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case modcamclaycoupmat:{
      e = cclayc[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    }
    case bbmcoupmat:{
      e = bbm[ip[ipp].idm[im]].give_actual_ym(ipp,im,ido);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    }
    case creepeffym:
      e = creffym[ip[ipp].idm[im]].give_actual_ym();
      nu = creffym[ip[ipp].idm[im]].give_actual_nu();
      ret = e/(1-2.0*nu)/3.0;
      break;
    case elasttimemat:
      e = eltimemat[ip[ipp].idm[im]].actual_modulus(ipp);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case winklerpasternak:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case shefpl:
    case glasgowmechmat:
    case glasgowdamage:
    case simvisc:
    case lemaitr:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepbaz:
    case consolidation:
      ncompo = givencompeqother(ipp, im); // this material should be last one before elastic
      ncompo -= givencompeqother(ipp, im+1);
      e = give_bulk_modulus(ipp, im+1, ido+ncompo);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e= give_bulk_modulus(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case viscoplasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2); // ido index shift up to the second nonartificial material
      e = give_bulk_modulus(ipp, im+2, ido+ncompo); // Young modulus will be obtained from the second nonartificial material
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
      ncompo = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+1);
      e = give_bulk_modulus(ipp, im+1, ido+ncompo); // these artificial materials sould not have internal variables, i.e. ncompo=0 -> no ido index shift
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case creep_damage:
      ncompo = givencompeqother(ipp, im); 
      ncompo -= givencompeqother(ipp, im+1); // number of other components of creep material model
      ncompo -= givencompeqother(ipp, im+2); // number of other components of damage material model
      e = give_bulk_modulus(ipp, im+1, ido+ncompo); // these artificial materials sould not have internal variables, i.e. ncompo=0 -> no ido index shift
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case shrinkagemat:
      ncompo = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the shrinkage material
      e = give_bulk_modulus(ipp, im+1, ido+ncompo);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material
      e = give_bulk_modulus(ipp, im+1, ido+ncompo);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case creepb3:
    case creeprs:
    case creepdpl:
      e = creep_give_actual_ym(ipp, im, ido);
      nu = give_actual_nu(ipp);
      ret = e/(1-2.0*nu)/3.0;
      break;
    case time_switchmat:
      e = tswmat[ip[ipp].idm[im]].give_actual_ym(ipp, im, ido);
      nu = tswmat[ip[ipp].idm[im]].give_actual_nu(ipp, im, ido);
      ret = e/(1-2.0*nu)/3.0;
      break;    
    case hypoplastusatthermat:
      ret = hypoplustherm[ip[ipp].idm[im]].give_bulk_modulus(ipp, ido);
      break;
    default:
      print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }

  return ret;
}


/**
  Function returns initial porosity especially for plasticity models.
  It is used for transfer required mechanical values between MEFEL and TRFEL.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of initial porosity.
  
  Created by Tomas Koudelka 8.10.2013
*/
double mechmat::give_iniporosity(long ipp, long im, long ido)
{
  long i, ncompo;
  double n_ini=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case modcamclaymat:{
    n_ini = cclay[i].give_iniporosity(ipp,ido);
    break;
  }
  case modcamclaycoupmat:{
    n_ini = cclayc[i].give_iniporosity(ipp,ido);
    break;
  }
  case bbmcoupmat:{
    n_ini = bbm[i].give_iniporosity(ipp,ido);
    break;
  }
  case doublestructuremat:{
    n_ini = dsm[i].give_iniporosity(ipp,ido);
    break;
  }
  case hypoplastmat:{
    n_ini = hypopl[i].give_iniporosity(ipp,ido);
    break;
  }
  case effective_stress:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    n_ini = give_iniporosity(ipp, im+1, ido+ncompo);
    break;
  case damage_plasticity:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+2);
    n_ini = give_iniporosity(ipp, im+2, ido+ncompo);
    break;
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return n_ini;
}


/**
   Function returns actual porosity especially for soil models.
   It is used for transfer required mechanical values between MEFEL and TRFEL.
   
   @param[in] ipp - integration point number in the mechmat ip array.
   @param[in] im - material index
   @param[in] ido - index of internal variables for given material in the ipp other array
   
   @return The function returns actual value of porosity.
   
   Created by Tomas Krejci 3/10/2016
*/
double mechmat::give_porosity(long ipp, long im, long ido)
{
  long i,ncompo;
  double n=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case effective_stress:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    n = give_porosity(ipp, im+1, ido+ncompo);
    break;    
  case modcamclaymat:{
    n = cclay[i].give_porosity(ipp, ido);
    break;
  }
  case modcamclaycoupmat:{
    n = cclayc[i].give_porosity(ipp, ido);
    break;
  }
  case bbmcoupmat:{
    n = bbm[i].give_porosity(ipp, ido);
    break;
  }
  case doublestructuremat:{
    n = dsm[i].give_porosity(ipp, ido);
    break;
  }
  case hypoplastusatthermat:{
    n = hypoplustherm[i].give_porosity(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  return n;
}

/**
   Function returns actual void ratio for soil models.
   It is used for transfer required mechanical values between MEFEL and TRFEL.
   
   @param[in] ipp - integration point number in the mechmat ip array.
   @param[in] im - material index
   @param[in] ido - index of internal variables for given material in the ipp other array
   
   @return The function returns actual value of void ratio.
   
   Created by Tomas Krejci 3/10/2016
*/
double mechmat::give_void_ratio(long ipp, long im, long ido)
{
  long i;
  double void_r=0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:
  case elisopdmat:{
    break;
  }
  case hypoplastusatthermat:{
    void_r = hypoplustherm[i].give_void_ratio(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return void_r;
}



/**
  The function returns rate of the volumetric strain at the given integration point.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by Tomas Koudelka 05.2018
*/
double mechmat::give_strain_vol_rate(long ipp, long im, long ido)
{
  long i,ncompo;
  double epsvr=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
  case elisomat:{
    epsvr = eliso[i].give_strain_vol(ipp,ido);
    break;
  }
  case effective_stress:
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    epsvr = give_strain_vol_rate(ipp, im+1, ido+ncompo);
    break;    
  case modcamclaymat:{
    epsvr = cclay[i].give_strain_vol_rate(ipp, ido);
    break;
  }
  case modcamclaycoupmat:{
    epsvr = cclayc[i].give_strain_vol_rate(ipp, ido);
    break;
  }
  case bbmcoupmat:{
    epsvr = bbm[i].give_strain_vol_rate(ipp, ido);
    break;
  }
  case doublestructuremat:{
    epsvr = dsm[i].give_strain_vol_rate(ipp, ido);
    break;
  }
  case hypoplastusatthermat:{
    epsvr = hypoplustherm[i].give_strain_vol_rate(ipp,ido);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  return epsvr;
}



/**
  Function returns consistency parameter for models of plasticity.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function returns actual value of consistency parameter.
  
  Created by Tomas Koudelka
*/
double mechmat::give_consparam (long ipp, long im, long ido)
{
  long i,ncompo;
  double cp=0.0;
  
  i = ip[ipp].idm[im];
  
  switch (ip[ipp].tm[im]){
    case simplas1d:{
      cp = spl1d[i].give_consparam (ipp,ido);
      break;
    }
    case jflow:{
      cp = j2f[i].give_consparam (ipp,ido);
      break;
    }
    case mohrcoul:{
      cp = mcoul[i].give_consparam (ipp,ido);
      break;
    }
    case mohrcoulparab:{
      cp = mcpar[i].give_consparam (ipp,ido);
      break;
    }
    case boermaterial:{
      cp = boerm[i].give_consparam (ipp,ido);
      break;
    }
    case druckerprager:{
      cp = drprm[i].give_consparam (ipp,ido);
      break;
    }
    case doubledrprager:{
      // cp = ddpm[i].give_consparam (ipp,ido);
      break;
    }
    case druckerprager2:{
      cp = drprm2[i].give_consparam (ipp,ido);
      break;
    }
    case modcamclaymat:{
      cp = cclay[i].give_consparam (ipp,ido);
      break;
    }
    case modcamclaycoupmat:{
      cp = cclayc[i].give_consparam (ipp,ido);
      break;
    }
    case bbmcoupmat:{
      cp = bbm[i].give_consparam (ipp,ido);
      break;
    }
    case doublestructuremat:{
      cp = dsm[i].give_consparam (ipp,ido);
      break;
    }
    case hissplasmat:{
      cp = hisspl[i].give_consparam (ipp,ido);
      break;
    }
    case damage_plasticity:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+2);
      cp = give_consparam(ipp, im+2, ido+ncompo);
      break;
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      cp = give_consparam(ipp, im+1, ido+ncompo);
      break;
    case damplifmat:
      cp = damplifm[i].give_consparam (ipp,ido);
      break;
    case plastifmat:
      cp = plastifm[i].give_consparam (ipp,ido);
      break;
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
  }
  
  return cp;
}



/**
  Function returns damage parameter omega for given material im.

  @param[in] ipp    - integration point number
  @param[in] im     - index of given material
  @param[in] ido    - index of the first internal variable for given material in the ipp eqother array
   
  @return Function returns damage parameter.

  Created by Tomas Koudelka
*/
double mechmat::give_dampar(long ipp, long im, long ido)
{
  long ncompo;
  double dam = 0.0;
  switch (ip[ipp].tm[im]){
    case elisomat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elisopdmat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elortomat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case scaldamage:{
      dam = scdam[ip[ipp].idm[im]].givedamage(ipp, ido);
      break;
    }
    case nonlocdamgmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      dam = give_dampar(ipp, im+1, ido+ncompo);
      break;
    }
    case shrinkagemat:{
      ncompo = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the shrinkage material
      dam = give_dampar(ipp, im+1, ido+ncompo);
      break;
    }
    case effective_stress:{
      ncompo = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the effective stress material
      dam = give_dampar (ipp, im+1, ido+ncompo);
      break;
    }
    case damplifmat:{
      dam = damplifm[ip[ipp].idm[im]].givedamage(ipp, ido);
      break;
    }

    default:{
      print_err("unknown material model is required",__FILE__,__LINE__,__func__);
    }
  }
  return dam;
}

/**
  Function returns length of process zone for given material im.

  @param[in] ipp    - integration point number
  @param[in] im     - index of given material
  @param[in] ido    - index of the first internal variable for given material in the ipp eqother array
   
  @return Function returns length of process zone.

  Created by Tomas Koudelka
*/
double mechmat::give_proczonelength (long ipp, long im, long ido)
{
  long ncompo;
  double l = 0.0;
  switch (ip[ipp].tm[im])
  {
    case elisomat:
      l = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elisopdmat:
      l = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elortomat:
      l = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case scaldamage:{
      l = scdam[ip[ipp].idm[im]].give_proczonelength(ipp, ido);
      break;
    }
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      l = give_proczonelength (ipp, im+1, ido+ncompo);
      break;
    }
    case nonlocdamgmat:
      l = nldamg[ip[ipp].idm[im]].give_proczonelength (ipp, im, ido);
      break;
    case shrinkagemat:
      ncompo  = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the shrinkage material    
      l = give_proczonelength (ipp, im+1, ido+ncompo);
      break;
    case damage_plasticity:{
      ncompo  = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+2); // number of components due to the damage material 
      l = give_proczonelength (ipp, im+2, ido+ncompo);
      break;
    }
    default:{
      print_err("unknown material model is required",__FILE__,__LINE__,__func__);
    }
  }
  return l;
}


/**
  Function returns crack width

  @param[in] ipp    - integration point number
  @param[in] im     - index of given material
  @param[in] ido    - index of the first internal variable for given material in the ipp eqother array
   
  @return Function returns crack width.

  Created by Tomas Koudelka
*/
double mechmat::give_crackwidth (long ipp, long im, long ido)
{
  long ncompo;
  double dam = 0.0;
  switch (ip[ipp].tm[im]){
    case elisomat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elisopdmat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case elortomat:
      dam = 0.0;   // due to passing data in METR where the elastic material may be assigned to some parts of the structure
      break;
    case scaldamage:{
      dam = scdam[ip[ipp].idm[im]].give_crackwidth (ipp, ido);
      break;
    }
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      dam = give_crackwidth (ipp, im+1, ido+ncompo);
      break;
    }
    case nonlocdamgmat:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      dam = give_crackwidth (ipp, im+1, ido+ncompo);
      break;
    case shrinkagemat:
      ncompo  = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the shrinkage material    
      dam = give_crackwidth (ipp, im+1, ido+ncompo);
      break;
    case damplifmat:
      dam = damplifm[ip[ipp].idm[im]].give_crackwidth (ipp, ido);
      break;
    default:{
      print_err("unknown material model is required",__FILE__,__LINE__,__func__);
    }
  }
  return dam;
}

/**
  Function returns actual value of equivalent strain norm for scalar damage models.
  It is used in nonlocal models and thus it returns value from OTHER array 
  i.e. from the nonequilibrium state !!!

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array
  @param[in] kappa - output parameter

  @return Function does not return anything

  Created by Tomas Koudelka, 10.2009
*/
void mechmat::give_kappa(long ipp,long im,long ido, vector &kappa)
{
  switch (ip[ipp].tm[im])
  {
    case scaldamage:
      kappa[0] = ip[ipp].other[ido];
      break;
    case scaldamagecc:
      kappa[0] = ip[ipp].other[ido];
      break;
    case fixortodamage:
      kappa[0] = ip[ipp].other[ido];
      break;  
    default:
      print_err ("unknown type of damage model is required",__FILE__,__LINE__,__func__);
  }  
  return;
}



/**
  Function returns irreversible strains for given material im.
  
  @param[in]  ipp    - integration point number
  @param[in]  im     - index of given material
  @param[in]  ido    - index of the first internal variable for given material in the ipp eqother array
  @param[out] epsirr - %vector of irreversible strains

  @return Function returns vector of irreversible strains via parameter epsirr.

  Created by Tomas Koudelka
*/
void mechmat::giveirrstrains(long ipp, long im, long ido, vector &epsirr)
{
  long ncompo;
  
  nullv(epsirr);
  switch (ip[ipp].tm[im]){
    case simplas1d:{
      spl1d[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case jflow:{
      j2f[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case mohrcoul:{
      mcoul[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case mohrcoulparab:{
      mcpar[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case boermaterial:{
      boerm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case druckerprager:{
      drprm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case doubledrprager:{
      //ddpm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case druckerprager2:{
      drprm2[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case chenplast:{
      //chplast[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case modcamclaymat:{
      cclay[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case modcamclaycoupmat:{
      cclayc[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case bbmcoupmat:{
      bbm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case doublestructuremat:{
      dsm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case hissplasmat:{
      hisspl[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      giveirrstrains(ipp, im+1, ido+ncompo, epsirr);
      break;
    }
    case creepbaz:{
      crbaz[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    }
    case creepeffym:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      giveirrstrains(ipp, im+1, ido+ncompo, epsirr);
      break;
    }
    case creepb3:
    case creeprs:
    case creepdpl:{
      creep_giveirrstrains(ipp, im, ido, epsirr);
      break;
    }
    case shrinkagemat:
      shmat[ip[ipp].idm[im]].giveirrstrains(ipp, im, ido, epsirr);
      break;
    case damplifmat:
      damplifm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    case plastifmat:
      plastifm[ip[ipp].idm[im]].giveirrstrains(ipp, ido, epsirr);
      break;
    default:
      print_err("unknown material model is required",__FILE__,__LINE__,__func__);
  }
}






// **********************************************************************
// **********************************************************************
//
// Functions for retrieving of non-mechanical quantities from int. points
//
// **********************************************************************
// **********************************************************************



/**
  The function searches for the total number and type of non-mechanical quantities 
  that are required by used material models in the problem. Required non-mechanical quantities
  are searched for at all integration points in all used material models.
  Resulting number of quantities is returned and new array is allocated
  for the required types of quantities if they are required.

  @param[out] rnmq - pointer to array of required types of non-mechanical quantities 
                     which is allocated inside function
 
  @return Function returns the number of required non-mechanical quantities.

  Created by Tomas Koudelka, 11.10.2013
*/
long mechmat::search_reqnmq(nonmechquant* &rnmq)
{
  long i, j, n;
  long anmq[tnknmq];
  memset(anmq, 0, sizeof(*anmq)*tnknmq);

  // search and mark required non-mechanical quantities at all int. points
  for (i=0; i<tnip; i++) 
  {
    for (j=0; j<ip[i].nm; j++)
      give_reqnmq(i, j, anmq);
  }

  // compute the number of required quantities
  n = 0; // decalred in mechmat
  for (i=0; i<tnknmq; i++)
  {
    if (anmq[i] == 1)
      n++;
  }

  // no non-mechanical quantities were required
  if (n == 0) 
  {
    rnmq = NULL;
    return n;
  }
  else
  {
    // store types of required non-mechanical quantities
    rnmq = new nonmechquant[n];
    j = 0;
    for (i=0; i<tnknmq; i++)
    {
      if (anmq[i] == 1){
        rnmq[j] = nonmechquant(i+1);
	j++;
      }
    }
  }

  return n;
}



/**
  The function marks required types of non-mechanical quantities by im-th
  material model in the given int. point ipp. Type of required quantities
  are marked by 1 in the array anmq whose length is equaled to the total number
  of known non-mechanical quantities.

  @param[in] ipp  - integration point id
  @param[in] im   - index of tested material model on the integration point
  @param[in] anmq - array with flags for used material types
                    anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                    anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.

  Created by Tomas Koudelka, 11.10.2013
*/
void mechmat::give_reqnmq(long ipp, long im, long *anmq)
{
  switch(ip[ipp].tm[im])
  {
    case elisomat:
    case elortomat:
    case elgmat3d:
    case elgmat2d:
    case winklerpasternak:
    case graphm:
    case elasttimemat:
    case varelisomat:
    case elisopdmat:
    case geoelast:
    case homomatm:
    case simplas1d:
    case jflow:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case modcamclaymat:
    case shefpl:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case creepbaz:
    case creepeffym:
    case contmat:
    case cebfipcontmat:
    case damplifmat:
    case plastifmat:
    case consolidation:
    case nonlocplastmat:
    case nonlocdamgmat:
    case nonlocalmicroM4:
    case glasgowdamage:
    case relaxationeuro:
    case simvisc:
    case lemaitr:
    case glasgowmechmat:
    case viscoplasticity:
    case viscoelasticity:
    case lenjonespot:
    case cusatismat:
    case creep_damage:
    case creepdpl:
    case damage_plasticity:
      break;

    //
    //  models that depend on non-mechanical quantities
    //
    case creeprs:
      crrs[im].give_reqnmq(anmq);
      break;
    case creepb3:
      crb3[im].give_reqnmq(anmq);
      break;
    case modcamclaycoupmat:
      cclayc[im].give_reqnmq(anmq);
      break;
    case bbmcoupmat:
      bbm[im].give_reqnmq(anmq);
      break;
    case doublestructuremat:
      dsm[im].give_reqnmq(anmq);
      break;
    case hypoplastmat:
      hypopl[im].give_reqnmq(anmq);
      break;
    case hypoplastusatthermat:
      hypoplustherm[im].give_reqnmq(anmq);
      break;
    case therisodilat:
      tidilat[im].give_reqnmq(anmq);
      break;
    case thervolisodilat:
      tvidilat[im].give_reqnmq(anmq);
      break;
    case effective_stress:
      effstr[im].give_reqnmq(anmq);
      break;
    case shrinkagemat:
      shmat[im].give_reqnmq(anmq);
      break;
    case time_switchmat:
      //tswmat[im].give_reqnmq(anmq);
      break;

    default:
      print_err("unknown material type is required", __FILE__, __LINE__, __func__);
  }

  return;
}



/**
  The function returns required non-mechanical quantity stored in the nonmechq array.
  It is used either in case of prescribed temperature changes in the input file or
  in case of coupled problems for retrieving of values passed between MEFEL and TRFEL.

  @param[in] qt  - type of required quantity
  @param[in] ipp - integration point id

  @return The function returns required non-mechanical quantity.

  Created by Tomas Koudelka, 7.10.2013
*/
double mechmat::givenonmechq(nonmechquant qt, long ipp)
{
  long id = nmqid[qt-1];
  if (id < 0)
  {
    print_err("Required quantity %s is not used in the problem solved", 
              __FILE__, __LINE__, __func__, nonmechquantstr[qt-1].alias);
    abort();
  }
  return nonmechq[id][ipp];
}



/** 
  The function returns selected components of array of non-mechanical quantities (nonmechq) 
  on given integration point.

  @param[in] ipp - integration point pointer
  @param[in] selcomp - definition of selected components
  @param[out] nmqcomp - output array of required nonmechanical quantities.

  @return The function returns required components of nonmechnical quantities in the argument nmqcomp.

  Created by Tomas Koudelka, 10.2023
*/
void mechmat::give_nonmechqcomp (long ipp, sel selcomp, vector &nmqcomp)
{
  long i, k;

  k=0;
  for (i=0; i<tnknmq; i++){
    if (selcomp.presence_id(i))
      nmqcomp[k] = nonmechq[i][ipp];
  }
}



/**
  The function stores required non-mechanical quantity into the nonmechq array.
  It is used either in case of prescribed temperature changes in the input file or
  in case of coupled problems for retrieving of values passed between MEFEL and TRFEL.

  @param[in] qt  - type of required quantity
  @param[in] ipp - integration point id
  @param[in] val - stored value

  @return The function stores required non-mechanical quantity.

  Created by Tomas Koudelka, 7.10.2013
*/
void mechmat::storenonmechq(nonmechquant qt, long ipp, double val)
{
  long id = nmqid[qt-1];
  if (id < 0)
  {
    print_err("Required quantity %s is not used in the problem solved", 
              __FILE__, __LINE__, __func__, nonmechquantstr[qt-1].alias);
    abort();
  }
  nonmechq[id][ipp] = val;
}



/**
  The function returns status of required non-mechanical quantity, i.e
  whether the given quantity is defined in the problem or not.

  @param[in] qt - type of required quantity
  
  @retval 0 - the quantity is not defined in the problem solved
  @retval 1 - the quantity is defined in the problem solved

  Created by Tomas Koudelka, 7.10.2013
*/
long mechmat::givestatusnmq (nonmechquant qt)
{
  if (nmqid[qt-1] < 0) // index was not set -> quantity is not defined
    return 0;

  return 1; // quantity is defined in all other cases
}



/**
  The function returns index of required non-mechanical quantity, i.e
  the first index of nonmechq array.

  @param[in] qt - type of required quantity
  
  @retval -1 - the quantity is not defined in the problem solved
  @retval >-1 - index of required the quantity

  Created by Tomas Koudelka, 7.10.2013
*/
long mechmat::givenonmechqid (nonmechquant qt)
{
  return nmqid[qt-1];
}






// **********************************************
// **********************************************
//
// Functions used in course of stress calculation
//
// **********************************************
// **********************************************



/**
  Function computes stiffness %matrix of the material
  in the required integration point.

  @param[out] d   - stiffness %matrix
  @param[in]  ipp - number of integration point
  @param[in]  im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in]  ido - index of internal variables in the ip's ipp other array
   
  @return The function returns stiffness %matrix in the parametr d.

  Created by JK, 17.7.2001
*/
void mechmat::matstiff (matrix &d,long ipp,long im,long ido)
{
  long i;
  long ncompo;
  switch (ip[ipp].tm[im]){
    case elisomat:{
      i=ip[ipp].idm[im];
      eliso[i].matstiff (d,ip[ipp].ssst);
      break;
    }
    case elortomat:{
      i=ip[ipp].idm[im];
      elorto[i].matstiff (d,ipp);
      break;
    }
    case elgmat3d:{
      i=ip[ipp].idm[im];
      elgm3d[i].matstiff (d);
      break;
    }
    case elgmat2d:{
      i=ip[ipp].idm[im];
      elgm2d[i].matstiff (d,ip[ipp].ssst);
      break;
    }
    case simplas1d:{
      i=ip[ipp].idm[im];
      spl1d[i].matstiff (d,ipp,ido);
      break;
    }
    case jflow:{
      i=ip[ipp].idm[im];
      j2f[i].matstiff (d,ipp,ido);
      break;
    }
    case microplaneM4:{
      i=ip[ipp].idm[im];
      mpM4[i].matstiff (d,ipp,ido);
      break;
    }
    case microsimp:{
      i=ip[ipp].idm[im];
      mpSIM[i].matstiff (d,ipp,ido);
      break;
    }
    case microfibro:{
      i=ip[ipp].idm[im];
      mpfib[i].matstiff (d,ipp,ido);
      break;
    }
    case mohrcoul:{
      i=ip[ipp].idm[im];
      mcoul[i].matstiff (d,ipp,ido);
      break;
    }
    case mohrcoulparab:{
      i=ip[ipp].idm[im];
      mcpar[i].matstiff (d,ipp,ido);
      break;
    }
    case boermaterial:{
      i=ip[ipp].idm[im];
      boerm[i].matstiff (d,ipp,ido);
      break;
    }
    case druckerprager:{
      i=ip[ipp].idm[im];
      drprm[i].matstiff (d,ipp,ido);
      break;
    }
    case doubledrprager:{
      i=ip[ipp].idm[im];
      ddpm[i].matstiff (d,ipp,ido);
      break;
    }
    case druckerprager2:{
      i=ip[ipp].idm[im];
      drprm2[i].matstiff (d,ipp,ido);
      break;
    }
    case chenplast:{
      i=ip[ipp].idm[im];
      chplast[i].matstiff (d,ipp,ido);
      break;
    }
    case modcamclaymat:{
      i=ip[ipp].idm[im];
      cclay[i].matstiff (d,ipp,im,ido);
      break;
    }
    case modcamclaycoupmat:{
      i=ip[ipp].idm[im];
      cclayc[i].matstiff (d,ipp,ido);
      break;
    }
    case bbmcoupmat:{
      i=ip[ipp].idm[im];
      bbm[i].matstiff (d,ipp,ido);
      break;
    }
    case doublestructuremat:{
      i=ip[ipp].idm[im];
      dsm[i].matstiff (d,ipp,ido);
      break;
    }
    case homomatm:{
      i=ip[ipp].idm[im];
      hommatm[i].matstiff (d,ip[ipp].ssst);
      break;
    }
    case hypoplastmat:{
      i=ip[ipp].idm[im];
      hypopl[i].matstiff (d,ipp,ido);
      break;
    }
    case hypoplastusatthermat:{
      i=ip[ipp].idm[im];
      hypoplustherm[i].matstiff (d,ipp,ido);
      break;
    }
    case shefpl:{
      i=ip[ipp].idm[im];
      shpl[i].matstiff (d,ipp,ido);
      break;
    }
    case hissplasmat:{
      i=ip[ipp].idm[im];
      hisspl[i].matstiff (d,ipp,ido);
      break;
    }
    case glasgowmechmat:{
      i=ip[ipp].idm[im];
      glasgmat[i].matstiff (d,ipp,ido);
      break;
    }
    case glasgowdamage:{
      i=ip[ipp].idm[im];
      glasgdam[i].matstiff (d,ipp,ido);
      break;
    }
    case creep_damage:{
      i=ip[ipp].idm[im];
      crdam[i].matstiff (d,ipp,im,ido);
      break;
    }
    case time_switchmat:{
      i=ip[ipp].idm[im];
      tswmat[i].matstiff (d,ipp,im,ido);
      break;
    }
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      matstiff (d,ipp,im+1,ido+ncompo);
      break;
    }
      /*
        case lemaitr:{
        i=ip[ipp].idm[im];
        lmtr[i].matstiff (d,ipp,ido);
        break;
        }
      */
    case scaldamage:{
      i=ip[ipp].idm[im];
      scdam[i].matstiff (d,ipp,ido);
      break;
    }
    case fixortodamage:{
      i=ip[ipp].idm[im];
      fixodam[i].matstiff (ipp,ido,d);
      break;
    }
    case scaldamagecc:{
      i=ip[ipp].idm[im];
      scdamcc[i].matstiff (d,ipp,ido);
      break;
    }
    case ortodamage:{
      i=ip[ipp].idm[im];
      ortdam[i].matstiff (ipp,ido,d);
      break;
    }
    case ortodamage2:{
      i=ip[ipp].idm[im];
      ortdam2[i].matstiff (d,ipp,ido);
      break;
    }
    case ortodamagerot:{
      i=ip[ipp].idm[im];
      ortdamrot[i].matstiff (d,ipp,ido);
      break;
    }
    case anisodamage:{
      i=ip[ipp].idm[im];
      anidam[i].matstiff (d,ipp,ido);
      break;
    }
    case anisodamagerot:{
      i=ip[ipp].idm[im];
      anidamrot[i].matstiff (d,ipp,ido);
      break;
    }
    case graphm:{
      i=ip[ipp].idm[im];
      grmat[i].matstiff (d,ipp);
      break;
    }
    case varelisomat:{
      i=ip[ipp].idm[im];
      veliso[i].matstiff (d,ipp);
      break;
    }
    case elisopdmat:{
      i=ip[ipp].idm[im];
      elisopd[i].matstiff (d, ipp, ido);
      break;
    }
    case geoelast:{
      i=ip[ipp].idm[im];
      geoel[i].matstiff (d,ipp,ido);
      break;
    }
    case creepdpl:
    case creeprs:
    case creepb3:{
      creep_matstiff (d,ipp,im,ido);
      break;
    }
    case creepbaz:{
      i=ip[ipp].idm[im];
      crbaz[i].matstiff (d, ipp, ido);
      break;
    }
    case creepeffym:{
      i=ip[ipp].idm[im];
      creffym[i].matstiff (d, ipp, im, ido);
      break;
    }
    case consolidation:{
      i=ip[ipp].idm[im];
      csol[i].matstiff (d,ipp);
      break;
    }
    case winklerpasternak:{
      i=ip[ipp].idm[im];
      wpast[i].matstiff (d,ipp);
      break;
    }

    case nonlocdamgmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      matstiff (d,ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocplastmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      matstiff (d,ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocalmicroM4:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      matstiff (d,ipp,im+1,ido+ncompo);
      break;
    }

    case damage_plasticity:{
      i=ip[ipp].idm[im];
      dampl[i].matstiff(d, ipp, ido);
      break;
    }

    case viscoplasticity:{
      i=ip[ipp].idm[im];
      visplas[i].matstiff (d,ipp,im,ido);
      break;
    }

    case viscoelasticity:{
      i=ip[ipp].idm[im];
      viselas[i].matstiff (d,ipp,im,ido);
      break;
    }

    case contmat:{
      i=ip[ipp].idm[im];
      conmat[i].matstiff (d,ipp);
      break;
    }

    case cebfipcontmat:{
      i=ip[ipp].idm[im];
      cebfipconmat[i].matstiff (d,ipp);
      break;
    }

    case damplifmat:{
      i=ip[ipp].idm[im];
      damplifm[i].matstiff(d, ipp, ido);
      break;
    }

    case plastifmat:{
      i=ip[ipp].idm[im];
      plastifm[i].matstiff(d, ipp, ido);
      break;
    }

    case shrinkagemat:{
      i=ip[ipp].idm[im];
      shmat[i].matstiff (d,ipp,im,ido);
      break;
    }

    case elasttimemat:{
      i=ip[ipp].idm[im];
      eltimemat[i].matstiff (d, ipp);
      break;
    }

    case layerplate:{
      i=ip[ipp].idm[im];
      lplate[i].matstiff (d,ipp,ido);
      break;
    }

    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
  }
}



/**
  Function computes elastic stiffness %matrix of the material
  in the required integration point.

  @param[out] d  - stiffness %matrix
  @param[in] ipp - number of integration point
  @param[in] ido  - index of internal variables in the ip's ipp other array
 
  @return The function returns elastic stiffness %matrix in the parametr d.

  Cretaed by JK, 17.7.2001
  modified by TKo, 7.2008
*/
void mechmat::elmatstiff (matrix &d, long ipp, long ido)
{
  long idem = ip[ipp].gemid();
  long ncompo = givencompother(ipp, 0) - givencompother(ipp, idem);  
  elmatstiff(d, ipp, idem, ido+ncompo, ip[ipp].ssst);
}



/**
  Function computes elastic stiffness %matrix of the material
  in the required integration point.

  @param[out] d   - stiffness %matrix
  @param[in] ipp  - number of integration point
  @param[in] im   - index of the required elastic material model
  @param[in] ido  - index of internal variables in the ip's ipp other array
  @param[in] ssst - indicator of stress state for which the stiffness matrix will be assembled
 
  @return The function returns elastic stiffness %matrix in the parametr d.

  Cretaed by JK, 17.7.2001
  modified by TKo, 7.2008
*/
void mechmat::elmatstiff (matrix &d, long ipp, long im, long ido, strastrestate ssst)
{
  long i;

  switch (ip[ipp].tm[im]){
  case elisomat:{
    i=ip[ipp].idm[im];
    eliso[i].elmatstiff (d, ssst);
    break;
  }
  case elortomat:{
    i=ip[ipp].idm[im];
    elorto[i].matstiff (d, ipp, ssst);
    break;
  }
  case elgmat3d:{
    i=ip[ipp].idm[im];
    elgm3d[i].elmatstiff (d);
    break;
  }
  case elgmat2d:{
    i=ip[ipp].idm[im];
    elgm2d[i].elmatstiff (d, ssst);
    break;
  }
  case varelisomat:{
    i=ip[ipp].idm[im];
    veliso[i].elmatstiff (d,ipp, ssst);
    break;
  }
  case elisopdmat:{
    i=ip[ipp].idm[im];
    long ncompo = givencompeqother(ipp, 0);
    ncompo -= givencompeqother(ipp, im);
    elisopd[i].elmatstiff (d, ipp, ido+ncompo);
    break;
  }
  case elasttimemat:{
    i=ip[ipp].idm[im];
    eltimemat[i].elmatstiff (d, ipp);
    break;
  }/*
  case bbmcoupmat:{
    i=ip[ipp].idm[im];
    bbm[i].elmatstiff (d, ipp, im, ido);
    break;
  }*/
  default:
    print_err("Unknown elastic material type %d is required (im=%ld, ipp=%ld, eid=%ld)",__FILE__,__LINE__,__func__, ip[ipp].tm[im], im, ipp, elip[ipp]+1);
  }
}



/**
  Function computes elastic stiffness %matrix of the material
  in the local coordinate system for the required integration point.

  @param[out] d   - local elastic stiffness %matrix
  @param[in]  ipp - number of integration point
 
  @return The function returns local elastic stiffness %matrix in the argument d.

  Created by Tomas Koudelka, 19.12.2014
*/
void mechmat::loc_elmatstiff (matrix &d,long ipp)
{
  long i,idem;

  idem = ip[ipp].gemid();
  switch (ip[ipp].tm[idem]){
  case elortomat:{
    i=ip[ipp].idm[idem];
    elorto[i].loc_matstiff (d,ip[ipp].ssst);
    break;
  }
  default:
    print_err("Unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function returns transformation %matrix from the local material coordinate system 
  to the global one (x_g = T . x_l).

  @param d[out]  - local elastic stiffness %matrix
  @param ipp[in] - number of integration point
 
  @return The function returns local elastic stiffness %matrix in the argument d.

  Created by Tomas Koudelka, 19.12.2014
*/
void mechmat::loc_transf_mat (matrix &tmat, long ipp)
{
  long i,idem;

  idem = ip[ipp].gemid();
  switch (ip[ipp].tm[idem]){
  case elortomat:{
    i=ip[ipp].idm[idem];
    elorto[i].give_transf_mat (tmat, ipp);
    break;
  }
  default:
    print_err("Unknown elastic material type is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function computes elastic compliance %matrix of the material
  in the required integration point.

  @param[out] c  - compliance %matrix
  @param[in] ipp - number of integration point
  @param[in] ido - index of internal variables in the ip's ipp other array
   
  @return The function returns elastic compliance %matrix in the parametr c.

  Created by JK, 3.2.2003
*/
void mechmat::elmatcompl (matrix &c, long ipp, long ido)
  //  27.10.2001
{
  long i, idem;

  idem = ip[ipp].gemid();
  switch (ip[ipp].tm[idem]){
  case elisomat:{
    i=ip[ipp].idm[idem];
    eliso[i].matcompl (c,ip[ipp].ssst);
    break;
  }
  case elgmat2d:{
    i=ip[ipp].idm[idem];
    elgm2d[i].matcompl (c,ip[ipp].ssst);
    break;
  }
  case varelisomat:{
    i=ip[ipp].idm[idem];
    veliso[i].matcompl (c,ipp);
    break;
  }
  case elisopdmat:{
    i=ip[ipp].idm[idem];
    elisopd[i].matcompl (c, ipp, ido);
    break;
  }
  case elasttimemat:{
    i=ip[ipp].idm[idem];
    eltimemat[i].matcompl (c, ipp);
    break;
  }
  default:
    print_err("unknown elastic material type is required",__FILE__,__LINE__,__func__);

  }
}


/**
  Function computes stiffness %matrix of the material
  in the required integration point.

  @param[out] d - stiffness %matrix
  @param[in]  ipp - number of integration point
  @param[in]  im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in]  ido - index of internal variables in the ip's ipp other array
   
  @return The function returns stiffness %matrix in the parametr d.

  Created by JK, 17.7.2001
*/
void mechmat::matdamp (matrix &d,long ipp,long im,long ido)
{
  long i;

  switch (ip[ipp].tm[im]){
  case viscoelasticity:{
    i=ip[ipp].idm[im];
    viselas[i].matdamp (d,ipp,im,ido);
    break;
  }
  case isovisc:{
    i=ip[ipp].idm[im];
    isovis[i].matdamp (d,ip[ipp].ssst);
    break;
  }
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);

  }
}



/**
   Function computes eigenstresses in required integration point.
   
   @param[in] sig - array of eigenstresses
   @param[in] ipp - number of integration point
   @param[in] id - id of material model
   
   @return The function returns elastic stiffness %matrix in the parametr d.
   
   Created by JK, 13. 6. 2013
*/
void mechmat::eigenstresses (vector &sig,long ipp,long id)
{
  long i,n;
  n=sig.n;
  vector eps(n);
  
  switch (ip[ipp].tm[id]){
  case relaxationeuro:{
    //  id of the relaxation model
    i=ip[ipp].idm[id];
    //  computation of actual eigenstresses
    relaxec[i].stress (sig,eps,ip[ipp].ssst);
    //  storage of actual eigenstrains
    storeeigstrain (ipp,eps);
    //  storage of actual eigenstresses
    storeeigstress (ipp,sig);    
    break;
  }
  default:
    print_err("Unknown model for eigenstresses computation is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function computes correct stresses for computed strains.
   
  @param[in] ipp - integration point pointer
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of the first internal variable in the ipp other arrays for given material model
                   The standard value is zero.
   
  @return The function does not return anything.

  Created by JK, 4.8.2001
  Modified by TKo, TKr
*/
void mechmat::computenlstresses (long ipp, intpoints &aip, long im, long ido)
{
  long ncompo;
  switch (ip[ipp].tm[im]){
    case elisomat:{
      eliso[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case elortomat:{
      elorto[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case elgmat3d:{
      elgm3d[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case elgmat2d:{
      elgm2d[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case simplas1d:{
      spl1d[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case jflow:{
      j2f[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case microplaneM4:{
      mpM4[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case microsimp:{
      mpSIM[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case microfibro:{
      mpfib[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case mohrcoul:{
      mcoul[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case mohrcoulparab:{
      mcpar[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case boermaterial:{
      boerm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case druckerprager:{
      drprm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case doubledrprager:{
      ddpm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case druckerprager2:{
      drprm2[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case chenplast:{
      chplast[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case modcamclaymat:{
      cclay[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case modcamclaycoupmat:{
      cclayc[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case bbmcoupmat:{
      bbm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    } 
    case doublestructuremat:{
      dsm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    } 
    case homomatm:{
      hommatm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case hypoplastmat:{
      hypopl[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case hypoplastusatthermat:{
      hypoplustherm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case hissplasmat:{
      hisspl[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case shefpl:{
      shpl[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case glasgowmechmat:{
      glasgmat[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case glasgowdamage:{
      glasgdam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case creep_damage:{
      crdam[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case shrinkagemat:{
      shmat[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case time_switchmat:{
      tswmat[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case effective_stress:{
      effstr[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
      /*
        case lemaitr:{
        lmtr[ip[ipp].idm[im]].nlstresses (ipp,ido);
        break;
        }
      */
    case scaldamage:{
      scdam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case fixortodamage:{
      fixodam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case scaldamagecc:{
      scdamcc[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamage:{
      ortdam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamage2:{
      ortdam2[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamagerot:{
      ortdamrot[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case anisodamage:{
      anidam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case anisodamagerot:{
      anidamrot[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case graphm:{
      grmat[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case varelisomat:{
      veliso[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case elisopdmat:{
      elisopd[ip[ipp].idm[im]].nlstresses (ipp, ido);
      break;
    }
    case geoelast:{
      geoel[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case creepb3:
    case creeprs:
    case creepdpl:{
      creep_nlstresses(ipp,im,ido);
      break;
    }
    case creepbaz:{
      crbaz[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case creepeffym:{
      creffym[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case consolidation:{
      csol[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case winklerpasternak:{
      wpast[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }

    case nonlocdamgmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      computenlstresses (ipp,aip,im+1,ido+ncompo);
      break;
    }
    case nonlocplastmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      computenlstresses (ipp,aip,im+1,ido+ncompo);
      break;
    }
    case nonlocalmicroM4:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      computenlstresses (ipp,aip,im+1,ido+ncompo);
      break;
    }
    case damage_plasticity:{
      dampl[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case viscoplasticity:{
      visplas[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case viscoelasticity:{
      viselas[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
  
    case elasttimemat:{
      eltimemat[ip[ipp].idm[im]].nlstresses (ipp, ido);
      break;
    }

    case contmat:{
      conmat[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case cebfipcontmat:{
      cebfipconmat[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case damplifmat:{
      damplifm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case plastifmat:{
      plastifm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case layerplate:{
      lplate[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case cusatismat:{
      cusmat[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }

    default:{
      print_err("unknown material type is required", __FILE__,__LINE__,__func__);
      abort();
    }
  }
}



/**
  Function computes correct increments of stresses for computed strains.
   
  @param[in] ipp - integration point pointer
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of the first internal variable in the ipp other arrays for given material model
                   The standard value is zero.
   
  @return The function does not return anything.

  Created by JK, TKo, 29.4.2008
*/
void mechmat::computenlstressesincr (long ipp,long im,long ido)
{
  long ncompo;
  
  switch (ip[ipp].tm[im]){   
  case creep_damage:{
    crdam[ip[ipp].idm[im]].nlstressesincr (ipp, im, ido);
    break;
  }
  case shrinkagemat:{
    shmat[ip[ipp].idm[im]].nlstressesincr (ipp, im, ido);
    break;
  }
  case time_switchmat:{
    tswmat[ip[ipp].idm[im]].nlstressesincr (ipp, im, ido);
    break;
  }
  case creepb3:
  case creeprs:
  case creepdpl:{
    creep_nlstressesincr (ipp,im,ido);
    break;
  }
  case creepeffym:{
    creffym[ip[ipp].idm[im]].nlstressesincr(ipp,im,ido);
    break;
  }
    /*
      case lemaitr:{
      lmtr[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
      }
      case creepbaz:{
      crbaz[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
      }
    */
  case consolidation:{
    csol[ip[ipp].idm[im]].nlstressesincr (ipp);
    break;
  }
  case effective_stress:{
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1); // ido index shift due to the effective stress material
    computenlstressesincr (ipp,im+1,ido+ncompo);
    break;
  }
  case viscoplasticity:{
    visplas[ip[ipp].idm[im]].nlstressesincr (ipp,im,ido);
    break;
  }
  case hypoplastusatthermat:{
    hypoplustherm[ip[ipp].idm[im]].nlstressesincr(ipp,im,ido);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[im]].nlstressesincr(ipp,im,ido);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[im]].nlstressesincr(ipp,im,ido);
    break;
  }
  case elisomat:
  case homomatm:
  case simplas1d:
  case jflow:
  case microplaneM4:
  case microsimp:
  case microfibro:
  case mohrcoul:
  case mohrcoulparab:
  case boermaterial:
  case druckerprager:
  case doubledrprager:
  case druckerprager2:
  case chenplast:
  case modcamclaymat:
  case modcamclaycoupmat:
  case shefpl:
  case glasgowmechmat:
  case glasgowdamage:
  case scaldamage:
  case fixortodamage:
  case scaldamagecc:
  case ortodamage:
  case ortodamage2:
  case ortodamagerot:
  case anisodamage:
  case anisodamagerot:
  case graphm:
  case varelisomat:
  case elisopdmat:
  case geoelast:
  case winklerpasternak:
  case nonlocdamgmat:
  case nonlocplastmat:
  case nonlocalmicroM4:
  case damage_plasticity:
  case cebfipcontmat:
  case contmat:
  case damplifmat:
  case plastifmat:
  case elasttimemat:
  case layerplate:
    if (im == 0)
      ip[ipp].clean_stresses(Mb->nlc);
    break;    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }

}



/**
  Function computes correct stresses for computed strains for non-local material models.

  @param[in] ipp - integration point pointer
  @param[in] im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido - index of the first internal variable in the ipp other arrays for given material model
                   The standard value is zero.

  @return The function does not return anything.

  Created by JK, TKo
*/
void mechmat::compnonloc_nlstresses (long ipp,long im,long ido)
{
  long ncompo;
  mattype mtype = ip[ipp].tm[im];
  switch (mtype){
    case elisomat:{
      eliso[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case simplas1d:{
      if (ip[ipp].hmt & 2) 
        spl1d[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      else
        spl1d[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case jflow:{
      if (ip[ipp].hmt & 2) 
        j2f[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      else
        j2f[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case microplaneM4:{
      if (ip[ipp].hmt & 2) 
        mpM4[ip[ipp].idm[im]].nonloc_nlstresses (ipp,ido);
      else
        mpM4[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case microsimp:{
      //    mpSIM[ip[ipp].idm[0]].nlstresses (ipp);
      break;
    }
    case microfibro:{
      //    mpfib[ip[ipp].idm[0]].nlstresses (ipp);
      break;
    }
    case mohrcoul:{
      if (ip[ipp].hmt & 2) 
        mcoul[ip[ipp].idm[im]].nonloc_nlstresses (ipp,ido);
      else
        mcoul[ip[ipp].idm[im]].nlstresses (ipp,ido);
      break;
    }
    case mohrcoulparab:{
      if (ip[ipp].hmt & 2) 
        mcpar[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      else
        mcpar[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case boermaterial:{
      if (ip[ipp].hmt & 2) 
        boerm[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      else
        boerm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case druckerprager:{
      if (ip[ipp].hmt & 2) 
        drprm[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      else
        drprm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case doubledrprager:{
      //ddpm[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        ddpm[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case druckerprager2:{
      // drprm[ip[ipp].idm[im]].nonloc_nlstresses (ipp,im,ido);
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        drprm2[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case modcamclaymat:{
      // cclay[ip[ipp].idm[im]].nonloc_nlstresses (ipp);
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        cclay[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case modcamclaycoupmat:{
      // cclay[ip[ipp].idm[im]].nonloc_nlstresses (ipp);
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        cclayc[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case hissplasmat:{
      // hisspl[ip[ipp].idm[im]].nonloc_nlstresses (ipp);
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        hisspl[ip[ipp].idm[im]].nlstresses (ipp, ido);
      break;
    }
    case effective_stress:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      compnonloc_nlstresses(ipp, im+1, ido+ncompo);
      break;
    }
      /*
        case lemaitr:{
        // lmtr[ip[ipp].idm[im]].nonloc_nlstresses (ipp);
        if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
        __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
        }
        else
        lmtr[ip[ipp].idm[im]].nlstresses (ipp);
        break;
        }*/
    case scaldamage:{
      scdam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case fixortodamage:{
      fixodam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case scaldamagecc:{
      scdamcc[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamage:{
      ortdam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamage2:{
      ortdam2[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case ortodamagerot:{
      ortdamrot[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case anisodamage:{
      anidam[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case anisodamagerot:{
      anidamrot[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case graphm:{
      grmat[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case varelisomat:{
      veliso[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case elisopdmat:{
      elisopd[ip[ipp].idm[im]].nlstresses (ipp, ido);
      break;
    }
    case creepbaz:{
      //    crbaz[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case consolidation:{
      csol[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }
    case winklerpasternak:{
      wpast[ip[ipp].idm[im]].nlstresses (ipp);
      break;
    }

    case creep_damage:{
      if (ip[ipp].hmt & 2) 
        crdam[ip[ipp].idm[im]].nonloc_nlstresses (ipp, im, ido);
      else
        crdam[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case shrinkagemat:{
      if (ip[ipp].hmt & 2) 
        shmat[ip[ipp].idm[im]].nonloc_nlstresses (ipp, im, ido);
      else
        shmat[ip[ipp].idm[im]].nlstresses (ipp, im, ido);
      break;
    }
    case creepb3:
    case creeprs:
    case creepdpl:{
      creep_nlstresses(ipp,im,ido);
      break;
    }
    case creepeffym:{
      creffym[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }
    case nonlocdamgmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      compnonloc_nlstresses (ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocplastmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      compnonloc_nlstresses (ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocalmicroM4:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      compnonloc_nlstresses (ipp,im+1,ido+ncompo);
      break;
    }
    
    case elasttimemat:{
      eltimemat[ip[ipp].idm[im]].nlstresses (ipp, ido);
      break;
    }
    case contmat:{
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        conmat[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case cebfipcontmat:{
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        cebfipconmat[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case damplifmat:{
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        damplifm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    case plastifmat:{
      if (ip[ipp].hmt & 2){
        print_err("non-local version of the model %s is not yet implemented (ipp=%ld, eid=%ld)",
                  __FILE__, __LINE__, __func__, mattype_kwdset.get_str(mtype), ipp, elip[ipp]+1);
        abort();
      }
      else
        plastifm[ip[ipp].idm[im]].nlstresses (ipp,im,ido);
      break;
    }

    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
  }
}



/**
  The function computes averaged strains.
  !!!??? The function is not used, linhex returns stresses
  !!!??? Candidate for removal.

  @return The function does not return anything

  Created by ???
*/
void mechmat::compute_averstrains ()
{
  long i,ne;
  double volume;
  elemtype te;
  vector averstra (6);
  
  volume=0.0;
  fillv (0.0,averstra);
  
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
      case linearhex:{
	Lhex->aver_strains (0,i,0,0,averstra,volume);
	break;
      }
      case quadrhex:{
	Qhex->aver_strains (0,i,0,0,averstra,volume);
	break;
      }
	
      default:{
	
      }
      }
    }
  }
  
  for (i=0;i<6;i++){
    averstra[i]/=volume;
  }
  
  printf ("\n\n prumerne deformace");
  for (i=0;i<6;i++){
    printf ("\n slozka %ld   %f",i+1,averstra[i]);
  }
  
}



/**
  Function computes strains caused by temperature.
   

  @return The function does not return anything.

  Created by JK,   
*/
void mechmat::temprstrains ()
{
  long i,nm,ii;
  
  for (i=0;i<tnip;i++){
    
    ii = ip[i].hmt & 1; // presence of thermal dilatation material model
    if (ii > 0){
      nm=ip[i].nm-1;
      switch (ip[i].tm[nm]){
      case therisodilat:{
	tidilat[ip[i].idm[nm]].temprstrains (i);
	break;
      }
      case thervolisodilat:{
	tvidilat[ip[i].idm[nm]].cumulstrains (i);
	//tvidilat[ip[i].idm[nm]].strains (i);
	break;
      }
      default:
	print_err("unknown thermal dilatation material is required",__FILE__,__LINE__,__func__);
      }
    }

  }
}



/**
  Function accumulates strains caused by temperature to tempstrains array on all integration points.
   

  @return The function does not return anything.

  Created by TKo, 10.2016   
*/
void mechmat::cumultemprstrains ()
{
  long i, ii;
  
  for (i=0;i<tnip;i++){
    
    ii = ip[i].hmt & 1; // presence of thermal dilatation material model
    if (ii > 0)
      cumultemprstrainsmat(i);
  }
}



/**
  Function accumulates strains caused by temperature to tempstrains array on the given integration point.
   
  @param[in] ipp - integration point pointer

  @return The function does not return anything.

  Created by TKo, 10.2016   
*/
void mechmat::cumultemprstrainsmat (long ipp)
{
  long nm = ip[ipp].nm-1; // index of the thermal dilation material

  switch (ip[ipp].tm[nm])
  {
    case therisodilat:
      tidilat[ip[ipp].idm[nm]].cumultemprstrains(ipp);
      break;      
    case thervolisodilat:
      tvidilat[ip[ipp].idm[nm]].cumulstrains(ipp);
      break;      
    default:
      print_err("unknown thermal dilatation material is required",__FILE__,__LINE__,__func__);
      abort();
  }
}



/**
  Function computes strains without eigenstrains.
   
  @return The function does not return anything.

  Created by JK, 3.3.2004
*/
void mechmat::eigstrmod ()
{
  long i,j,k,nc,ipp,nip;
  
  if (eigstrains != NULL){
    for (i=0;i<Mt->ne;i++){
      if (Gtm->leso[i]==1){
	//  first integration point on element
	ipp=Mt->elements[i].ipp[0][0];
	//  total number of integration points on element
	nip = Mt->give_tnip (i);
	for (j=0;j<nip;j++){
	  //  number of components in strain array
	  nc=ip[ipp].ncompstr;
	  for (k=0;k<nc;k++)
	    ip[ipp].strain[k]-=eigstrains[ipp][k];
	  ipp++;
	}
      }
    }
  }
  if (tempstrains != NULL){
    for (i=0;i<Mt->ne;i++){
      if (Gtm->leso[i]==1){
	//  first integration point on element
	ipp=Mt->elements[i].ipp[0][0];
	//  total number of integration points on element
	nip = Mt->give_tnip (i);
	for (j=0;j<nip;j++){
	  //  number of components in strain array
	  nc=ip[ipp].ncompstr;
	  for (k=0;k<nc;k++)
	    ip[ipp].strain[k]-=tempstrains[ipp][k];
	  ipp++;
	}
      }
    }
  }
}



/**
  Function computes strains without eigenstrains at auxiliary integration points.
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 15.6.2018
*/
void mechmat::aip_eigstrmod ()
{
  long k, eid, nc, app;
  
  if (aip_eigstrains != NULL)
  {
    for (app=0; app < tnaip; app++)
    {
      eid = elaip[app];      
      if (Gtm->leso[eid]==1)
      {
        //  number of components in strain array
        nc = aip[app].ncompstr;
        for (k=0; k<nc; k++)
          aip[app].strain[k] -= aip_eigstrains[app][k];
      }
    }
  }

  if (aip_tempstrains != NULL)
  {
    for (app=0; app<tnaip; app++)
    {
      eid = elaip[app];      
      if (Gtm->leso[eid]==1)
      {
        nc = aip[app].ncompstr;
        for (k=0; k<nc; k++)
          aip[app].strain[k] -= tempstrains[app][k];
      }
    }
  }
}



/**
   Function adds macro-strain components to the strain %vector which
   contains fluctuation strains. Function modifies strains in all integration points.

   @param[in] lcid - load case id
   @param[in] macrostrains - macrostarin %vector
   
   @return The function does not return anything but modifies strain array of all int. points.
   
   Created by JK, 16. 4. 2014
*/
void mechmat::add_macro_strains (long lcid, vector &macrostrains)
{
  long i,j,ncompstr;
  
  for (i=0;i<tnip;i++){
    ncompstr = ip[i].ncompstr;
    if (ncompstr < macrostrains.n){
      print_err("the number of strain components (%ld) in integration point %ld\n"
		"is not equal to the number of macro-strain component (%ld) in homogenization",
		__FILE__, __LINE__, __func__,ncompstr,i,macrostrains.n);
      abort ();
    }
    for (j=0;j<macrostrains.n;j++){
      ip[i].strain[ncompstr*lcid+j] += macrostrains[j];
    }
  }
}



/**
  Function computes the total strains at all integration points.
   
  @return The function does not return anything.

  Created by JK, 3.3.2004
*/
void mechmat::totstrains ()
{
  long i,j,k,nc,ipp,nip;
  
  if (eigstrains != NULL){
    for (i=0;i<Mt->ne;i++){
      if (Gtm->leso[i]==1){
	//  first integration point on element
	ipp=Mt->elements[i].ipp[0][0];
	//  total number of integration points on element
	nip = Mt->give_tnip (i);
	for (j=0;j<nip;j++){
	  //  number of components in strain array
	  nc=ip[ipp].ncompstr;
	  for (k=0;k<nc;k++)
	    ip[ipp].strain[k]+=eigstrains[ipp][k];
	  ipp++;
	}
      }
    }
  }
  if (tempstrains != NULL){
    for (i=0;i<Mt->ne;i++){
      if (Gtm->leso[i]==1){
	//  first integration point on element
	ipp=Mt->elements[i].ipp[0][0];
	//  total number of integration points on element
	nip = Mt->give_tnip (i);
	for (j=0;j<nip;j++){
	  //  number of components in strain array
	  nc=ip[ipp].ncompstr;
	  for (k=0;k<nc;k++)
	    ip[ipp].strain[k]+=tempstrains[ipp][k];
	  ipp++;
	}
      }
    }
  }
}



/**
  Function computes %matrix of thermal dilatancy of the material
  in the required integration point.

  @param[out] d - %matrix of thermal dilatancy
  @param[in] ipp - number of integration point

  @return The function does not return anything.

  Created by JK, 29.11.2002
*/
void mechmat::matdilat (matrix &d,long ipp)
{
  long i, iddm = ip[ipp].nm-1;

  switch (ip[ipp].tm[iddm]){
  case therisodilat:{
    i=ip[ipp].idm[iddm];
    tidilat[i].matdilat (d,ip[ipp].ssst);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }

}



/**
  Function stores volume belonging to the given integration point.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] vol - volume belonging to the integration point  
  
  @return The function does not return anything.

  Created by JK,
*/
void mechmat::storeipvol (long ipp,double vol)
{
  ipv[ipp]=vol;
}



/**
  Function averages local values and stores them in nonlocal values.

  @param[in] ipp - integration point pointer
  @param[in] im - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by JK, TKo,   4.8.2003
*/
void mechmat::nonlocaverage (long ipp,long im,long ido)
{
  long ncompo;

  if (ip[ipp].hmt & 2)
  {
    switch (ip[ipp].tm[im]){
    case nonlocplastmat :
      nlplast[ip[ipp].idm[im]].average(ipp, ido);
      break;
    case nonlocdamgmat :
      nldamg[ip[ipp].idm[im]].average(ipp,im,ido);
      break;
    case nonlocalmicroM4 :
      nmpM4[ip[ipp].idm[im]].average(ipp);
      break;
    case creep_damage:
      ncompo  = givencompeqother(ipp,im);   
      ncompo -= givencompeqother(ipp,im+2); // number of components due to the creep material
      nonlocaverage(ipp, im+2, ido+ncompo);
      break;
    case shrinkagemat:
      ncompo  = givencompeqother(ipp,im);  
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the shrinkage material    
      nonlocaverage(ipp, im+1, ido+ncompo);
      break;
    default:{
      print_err("unknown type of averaging model is required", __FILE__, __LINE__, __func__);
    }
    }
  }
}



/**
  The function returns the number of averaged local values.

  @param[in] ipp - integration point pointer
  @param[in] im  - material index

  @return The function returns the number of averaged local values.


  Created by Tomas Koudelka, 10.2009
*/
long mechmat::give_num_averq(long ipp,long im)
{
  long ret = 0;
  if (ip[ipp].hmt & 2)
  {
    switch (ip[ipp].tm[im]){
    case nonlocalmicroM4 :
      ret = nmpM4[ip[ipp].idm[im]].give_num_averq(ipp);
      break;
    case nonlocplastmat :
      ret = nlplast[ip[ipp].idm[im]].give_num_averq(ipp);
      break;
    case nonlocdamgmat :
      ret = nldamg[ip[ipp].idm[im]].give_num_averq(ipp, im);
      break;
    case creep_damage:
      ret = give_num_averq(ipp, im+2);
      break;
    case shrinkagemat:
      ret = give_num_averq(ipp, im+1);
      break;
    default:{
      print_err("unknown type of averaging model is required", __FILE__, __LINE__, __func__);
    }
    }
  }
  return ret;
}



/**
  The function returns vector of averaged quantities.

  @param[in] ipp - integration point pointer
  @param[in] im  - material index
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 10.2009
*/
void mechmat::give_aver_quantv(long ipp,long im,long ido, vector &qv)
{
  long ncompo;

  if (ip[ipp].hmt & 2)
  {
    switch (ip[ipp].tm[im]){
    case nonlocplastmat :
      nlplast[ip[ipp].idm[im]].give_aver_quantv(ipp, im, ido, qv);
      break;
    case nonlocdamgmat :
      nldamg[ip[ipp].idm[im]].give_aver_quantv(ipp, im, ido, qv);
      break;
    case creep_damage:
      ncompo  = givencompeqother(ipp,im);
      ncompo -= givencompeqother(ipp,im+2); // number of components due to the creep material
      give_aver_quantv(ipp, im+2, ido+ncompo, qv);
      break;
    case shrinkagemat:
      ncompo  = givencompeqother(ipp,im); 
      ncompo -= givencompeqother(ipp,im+1); // number of components due to the shrinkage material
      give_aver_quantv(ipp, im+2, ido+ncompo, qv);
      break;
    default:{
      print_err("unknown type of averaging model is required", __FILE__, __LINE__, __func__);
    }
    }
  }
  return;
}



/**
  Function returns radius of averaged neighbourhood for the given non-local material model.
   
  @param[in] ipp - integration point pointer
  @param[in] im - material index
  
  @return Function returns real number representing radius for nonlocal model.

  Created by Tomas Koudelka
*/
double mechmat::nonlocradius (long ipp,long im)
{
  double r=0.0;
  
  switch (ip[ipp].tm[im]){
  case nonlocplastmat :
  {
    r = nlplast[ip[ipp].idm[im]].r;
    break;
  }
  case nonlocdamgmat :
  {
    r = nldamg[ip[ipp].idm[im]].r;
    break;
  }
  case nonlocalmicroM4 :
  {
    r = nmpM4[ip[ipp].idm[im]].r;
    break;
  }
  default:{
    print_err("wrong material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return r;
}



/**
  Function searches index of nonlocal material type at given integration point.

  @param[in] ipp - integration point id

  @retval >= 0 in case that the nonlocal material type was found
  @retval < 0  in case that the integration point has no nonlocal material type

  Created by Tomas Koudelka, 8.2008
*/
long mechmat::givenonlocid(long ipp)
{
  long i, ret=-1;
  for(i=0; i<ip[ipp].nm; i++)
  {
    switch(ip[ipp].tm[i])
    {
      case nonlocplastmat:
      case nonlocdamgmat:
      case nonlocalmicroM4:
        ret = i;
        break;
      default:
        break;
    }
    if (ret >= 0)
      break;
  }
  return ret;
}



/**
  Function updates values of internal variables (eqother array) at integration points,
  it is used for nonlinear computations.

  @return The function does not return anything.

  Created by JK,  
  Modified by TKr, 17.7.2006 
*/
void mechmat::updateipval (void)
{
  long i,j,ipp,nip;
  
  //updates only used elements(int. points)
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      ipp=Mt->elements[i].ipp[0][0];
      nip=Mt->give_tnip (i);
      for (j=0;j<nip;j++){
	updateipvalmat (ipp,0,0);
	ipp++;
      }
    }
  }
}


/**
  The function updates values at auxiliary integration points.

  @return The function does not return anything but it changes content of auxiliary int. points.

  Created by Tomas Koudelka, 25.5.2018
*/
void mechmat::update_aipval (void)
{
  long app, eid;
  intpoints *tmp_ip;
  long      *tmp_elip;
  double   **tmp_nonmechq;
  double   **tmp_ic;
  double   **tmp_eigstrains;
  double   **tmp_eigstresses;
  double   **tmp_tempstrains;
  
  if (Mm->aip == NULL) // there are no auxiliary int. points
    return;

  // swap regular integration point auxiliary integration point arrays
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

  //updates only used elements(int. points)
  // just call updateipvalmat procedure on all auxiliary integration points
  for (app=0; app<Mm->tnaip; app++)
  {
    eid = Mm->elip[app];
    if (Gtm->leso[eid] == 1)
      updateipvalmat (app,0,0);
  }

  // swap regular integration point and auxiliary integration point arrays
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
  Function updates values of internal variables (eqother array) at integration points,
  it is used for nonlinear computations.
   
  @param[in] ipp - number of integration point
  @param[in] im - index of material type
  @param[in] ido - index in array other
  
  Created by JK,
  Modified by TKo, TKr
*/
void mechmat::updateipvalmat (long ipp,long im,long ido)
{
  long ncompo;
  switch (ip[ipp].tm[im]){
    case elisomat :
      eliso[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    case homomatm:
      hommatm[ip[ipp].idm[im]].updateval(ipp,ido);
      break;
    case elortomat :
      break;
    case winklerpasternak:
      break;
    case simplas1d:{
      spl1d[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    }
    case jflow:{
      j2f[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case mohrcoul:{
      mcoul[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case mohrcoulparab:{
      mcpar[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case boermaterial:{
      boerm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case druckerprager:{
      drprm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case doubledrprager:{
      ddpm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case druckerprager2:{
      drprm2[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case chenplast:{
      chplast[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case modcamclaymat:{
      cclay[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    }
    case modcamclaycoupmat:{
      cclayc[ip[ipp].idm[im]].updateval (ipp, im, ido);
      break;
    }
    case bbmcoupmat:{
      bbm[ip[ipp].idm[im]].updateval (ipp,im, ido);
      break;
    }
    case doublestructuremat:{
      dsm[ip[ipp].idm[im]].updateval (ipp,im, ido);
      break;
    }
    case hypoplastmat:{
      hypopl[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case hypoplastusatthermat:{
      hypoplustherm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case hissplasmat:{
      hisspl[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    }
    case effective_stress:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      updateipvalmat(ipp, im+1, ido+ncompo);
      break;
    case microplaneM4:{
      mpM4[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case microsimp:{
      mpSIM[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case microfibro:{
      mpfib[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case shefpl:{
      shpl[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    }
    case glasgowdamage:{
      glasgdam[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case creep_damage:{
      crdam[ip[ipp].idm[im]].updateval (ipp, im, ido);
      break;
    }
    case shrinkagemat:{
      shmat[ip[ipp].idm[im]].updateval (ipp, im, ido);
      break;
    }
    case time_switchmat:{
      tswmat[ip[ipp].idm[im]].updateval (ipp, im, ido);
      break;
    }
    case graphm:{
      break;
    }
    case varelisomat :
      break;
    case elisopdmat :
      elisopd[ip[ipp].idm[im]].updateval (ipp,ido);
      break;
    case elasttimemat:
      eltimemat[ip[ipp].idm[im]].updateval (ipp, im, ido);
      break;
    case geoelast:{
      geoel[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case scaldamage:{
      scdam[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case fixortodamage:{
      fixodam[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case scaldamagecc:{
      scdamcc[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }    
    case ortodamage:{
      ortdam[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case ortodamage2:{
      ortdam2[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case ortodamagerot:{
      ortdamrot[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case anisodamage:{
      anidam[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    case anisodamagerot:{
      anidamrot[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }

    case nonlocdamgmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      updateipvalmat (ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocplastmat:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      updateipvalmat (ipp,im+1,ido+ncompo);
      break;
    }
    case nonlocalmicroM4:{
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      updateipvalmat (ipp,im+1,ido+ncompo);
      break;
    }
    
    case damage_plasticity:{
      dampl[ip[ipp].idm[im]].updateval (ipp);
      break;
    }
    case simvisc:
    case lemaitr:
    case viscoplasticity:
      visplas[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    case creeprs:
    case creepb3:
    case creepdpl:{
      creep_updateval (ipp,im,ido);
      break;
    }
    case creepbaz:
      crbaz[ip[ipp].idm[im]].updateval ();
      break;
    case creepeffym:
      creffym[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    case lenjonespot:
      break;
    
    case contmat:
      conmat[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;

    case cebfipcontmat:
      cebfipconmat[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    
    case damplifmat:
      damplifm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    
    case plastifmat:
      plastifm[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;

    case layerplate:
      lplate[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    
    case consolidation:{
      csol[ip[ipp].idm[im]].updateval (ipp,im,ido);
      break;
    }
    default:{
      print_err ("unknown material type is required",__FILE__,__LINE__,__func__);
      break;
    }
  }
}



/**
  Function copyies equilibrium values of internal variables (eqother array) to the other array at integration points,
  it is used for nonlinear computations after restart from backup. Actualized values in other array may be used for the 
  computation of tangent stiffness %matrix.

  @return The function does not return anything.

  Created by TKo, 7.6.2018  
*/
void mechmat::updateother (void)
{
  long i,j,ipp,nip;
  
  //updates only used elements(int. points)
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      ipp=Mt->elements[i].ipp[0][0];
      nip=Mt->give_tnip (i);
      for (j=0;j<nip;j++){
	updateothermat (ipp,0,0);
	ipp++;
      }
    }
  }
}



/**
  Function copyies equilibrium values of internal variables (eqother array) to the other array at one integration point,
  it is used for nonlinear computations after restart from backup. Actualized values in other array may be used for the 
  computation of tangent stiffness %matrix.
   
  @param[in] ipp - number of integration point
  @param[in] im - index of material type
  @param[in] ido - index in array other
  
  Created by TKo, 7.6.2018
*/
void mechmat::updateothermat (long ipp,long im,long /*ido*/)
{
  switch (ip[ipp].tm[im]){
  case elisomat :
    break;
  case elortomat :
    break;
  case winklerpasternak:
    break;
  case varelisomat :
    break;
  case elisopdmat :
    break;
  case lenjonespot:
    break;
  case simplas1d:
  case jflow:
  case mohrcoul:
  case mohrcoulparab:
  case boermaterial:
  case druckerprager:
  case doubledrprager:
  case druckerprager2:
  case chenplast:
  case modcamclaymat:
  case modcamclaycoupmat:
  case bbmcoupmat:
  case doublestructuremat:
  case hypoplastmat:
  case hypoplastusatthermat:
  case hissplasmat:
  case effective_stress:
  case microplaneM4:
  case microsimp:
  case microfibro:
  case shefpl:
  case glasgowdamage:
  case creep_damage:
  case shrinkagemat:
  case time_switchmat:
  case graphm:
  case elasttimemat:
  case geoelast:
  case scaldamage:
  case fixortodamage:
  case scaldamagecc:
  case ortodamage:
  case ortodamage2:
  case ortodamagerot:
  case anisodamage:
  case anisodamagerot:
  case nonlocdamgmat:
  case nonlocplastmat:
  case nonlocalmicroM4:    
  case damage_plasticity:
  case simvisc:
  case lemaitr:
  case viscoplasticity:
  case creeprs:
  case creepb3:
  case creepdpl:
  case creepbaz:
  case creepeffym:
  case contmat:
  case cebfipcontmat:
  case damplifmat:
  case plastifmat:
  case consolidation:
    memcpy(ip[ipp].other, ip[ipp].eqother, sizeof(*ip[ipp].other)*ip[ipp].ncompother);
    break;
  default:{
    print_err ("unknown material type is required",__FILE__,__LINE__,__func__);
    break;
  }
  }
}



/**
  The function updates actual material index for time switched materials.
  It must be called in time solvers before each time step and optionally, it performes also 
  initialization of actual material model.

  @param[in] lcid - load case id obtained from solver
  @param[in] init - flag for initialization of actualized material model (0 = no call of initvalues, 1=initval is called)
  @param[in] rinit - flag for initialization after restorage from hdbackup

  @retval 0 - if no material index has been changed
  @retval 1 - if some material indeces have been actualized

  Created by Tomas Koudelka, 16.3.2016
*/
long mechmat::update_actual_mat_id(long lcid, long init, bool rinit)
{
  long i,j,k,ipp,nip,app;
  long oami,ami,ret = 0;
  intpoints *tmp_ip;
  long *tmp_elip;
  double **tmp_nmq, **tmp_ic;
  
  for (k=0; k<ntswmat; k++)
  {
    // save original actual material index
    oami = tswmat[k].ami;
    // it is suppossed that time_switchmat is the first material in the material model chain
    // otherwise the following code must be rewritten
    tswmat[k].update_ami(0);
    ami = tswmat[k].ami;
    
    if ((ami < 0))
    {
      print_err("%ld. time switched material returns no active material at time %le", __FILE__, __LINE__, __func__, k+1, Mp->time);
      abort();
    }

    // init == 1 => initialization of all integration points with the k-th timeswitching material
    for (i=0;i<Mt->ne;i++)
    {
      // updates only used elements(int. points)
      if (Gtm->leso[i]==1)
      {
        ipp=Mt->elements[i].ipp[0][0];
        nip=Mt->give_tnip (i);
        for (j=0; j<nip; j++, ipp++)
        {
          if (ami == oami)
            continue;
          else 
          {
            // it is suppossed that time_switchmat is the first material in the material model chain
            // otherwise the following condition must be rewritten
            if ((ip[ipp].tm[0] == time_switchmat) && (ip[ipp].idm[0] == k))
            {
              ret = 1;
              ip[ipp].hmt = 0;
              update_hmt(ipp, 0);
              ip[ipp].ncompother = givencompother(ipp, ami);
              ip[ipp].ncompeqother = givencompeqother(ipp, ami);
              ip[ipp].clean_other();
              if (tswmat[k].store_eigstr == yes)
                storeeqother(ipp, 0, ip[ipp].ncompstr, ip[ipp].strain); 
              // if the material model has changed (ami != oami) then the 
              // initialization of the material model must be performed if required (init == 1)
              if (init)
                tswmat[ip[ipp].idm[0]].initvalues (lcid, ipp, 0, 0, rinit);
            }
          }
        }
      }
    }

    // swap regular integration point auxiliary integration point arrays
    // from now, all functions will work with AUXILIARY integration points 
    tmp_ip = ip;
    tmp_elip = elip;
    tmp_nmq = nonmechq;
    tmp_ic  = ic;
    ip = aip;
    elip = elaip;
    nonmechq = aip_nonmechq;
    ic = aip_ic;
    // init == 1 => initialization of all integration points with the k-th timeswitching material
    for (app=0; app<tnaip; app++)
    {
      // updates only used elements(int. points)
      if (Gtm->leso[elip[app]]==1)
      {
        if (ami == oami)
          continue;
        else 
        {
          // it is suppossed that time_switchmat is the first material in the material model chain
          // otherwise the following condition must be rewritten
          if ((ip[app].tm[0] == time_switchmat) && (ip[app].idm[0] == k))
          {
            ret = 1;
            ip[app].hmt = 0;
            update_hmt(app, 0);
            ip[app].ncompother = givencompother(app, ami);
            ip[app].ncompeqother = givencompeqother(app, ami);
            ip[app].clean_other();
            if (tswmat[k].store_eigstr == yes)
              storeeqother(app, 0, ip[app].ncompstr, ip[app].strain); 
            // if the material model has changed (ami != oami) then the 
            // initialization of the material model must be performed if required (init == 1)
            if (init)
              tswmat[aip[app].idm[0]].initvalues (lcid, app, 0, 0, rinit);
          }
        }
      }
    }
    // swap regular integration point auxiliary integration point arrays
    // from now, all functions will work with REGULAR integration points 
    ip = tmp_ip;
    elip = tmp_elip;
    nonmechq = tmp_nmq;
    ic = tmp_ic;
  }
  return ret;
}



/**
  Function cleans all arrays defined at integration points.
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void mechmat::clean_ip ()
{
  long i;
  for (i=0;i<tnip;i++){
    ip[i].clean (Mb->nlc);
  }
}



/**
  Function checks time step size requirements which originate in material models.

  @return The function returns minimum ratio of required time step size to the actual one 
          collected from all integration points or dstep_r = 1.0 for no change in time step size.

  Created by Tomas Koudelka, 10.3.2014
*/
double mechmat::dstep_red_mat()
{
  long i;
  double ret = 1.0;
  double dstep_r;
  
  for (i=0; i<tnip; i++)
  {
    dstep_r = dstep_red_ip(i, 0, 0);
    if (dstep_r < ret)
      ret = dstep_r;
  }
  return ret;
}



/**
  Function returns time step size required by material models in the given integration point.

  @param[in] ipp - number of integration point
  @param[in] im - index of material type (default is 0)
  @param[in] ido - index in array other  (default is 0)


  @return The function returns minimum ratio of required time step size to the actual one 
          for the given integration points or dstep_r = 1.0 for no change in time step size.

  Created by Tomas Koudelka, 10.3.2014
*/
double mechmat::dstep_red_ip(long ipp, long im, long ido)
{
  double ret = 1.0;
  long ncompo;
  
  switch (ip[ipp].tm[im])
    {
    case elisomat :
    case elortomat :
    case homomatm:
    case simplas1d:
    case jflow:
    case mohrcoul:
    case mohrcoulparab:
    case boermaterial:
    case druckerprager:
    case doubledrprager:
    case druckerprager2:
    case chenplast:
    case modcamclaymat:
    case hissplasmat:
    case microplaneM4:
    case microsimp:
    case microfibro:
    case shefpl:
    case glasgowdamage:
    case creep_damage:
    case shrinkagemat:
    case time_switchmat:
    case graphm:
    case varelisomat :
    case elisopdmat :
    case geoelast:
    case scaldamage:
    case fixortodamage:
    case scaldamagecc:
    case ortodamage:
    case ortodamage2:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:
    case damage_plasticity:
    case simvisc:
    case lemaitr:
    case viscoplasticity:
    case creeprs:
    case creepb3:
    case creepdpl:
    case creepbaz:
    case creepeffym:
    case lenjonespot:
    case contmat:
    case cebfipcontmat:
    case damplifmat:
    case plastifmat:
    case winklerpasternak:
    case consolidation:
    case elasttimemat:
      ret = 1.0;
      break;

    case effective_stress:
    case nonlocdamgmat:
    case nonlocplastmat:
    case nonlocalmicroM4:
      ncompo  = givencompeqother(ipp, im);
      ncompo -= givencompeqother(ipp, im+1);
      ret = dstep_red_ip(ipp,im+1,ido+ncompo);
      break;

    case modcamclaycoupmat:
      ret = cclayc[ip[ipp].idm[im]].dstep_red(ipp, im, ido);
      break;
    case bbmcoupmat:
      ret = bbm[ip[ipp].idm[im]].dstep_red(ipp, im, ido);
      break;
    case doublestructuremat:
      ret = dsm[ip[ipp].idm[im]].dstep_red(ipp, im, ido);
      break;
    case hypoplastmat:
      ret = hypopl[ip[ipp].idm[im]].dstep_red(ipp,im,ido);
      break;
    case hypoplastusatthermat:
      ret = hypoplustherm[ip[ipp].idm[im]].dstep_red(ipp,im,ido);
      break;
    case layerplate:
      ret = lplate[ip[ipp].idm[im]].dstep_red(ipp,im,ido);
      break;
    default:{
      print_err ("unknown material type is required",__FILE__,__LINE__,__func__);
      break;
    }  
  }

  return ret;
}





// *****************************************************************
// *****************************************************************
//   PART CONTAINING FUNCTIONS DEALING WITH PLASTICITY
// *****************************************************************
// *****************************************************************



/**
  Function evaluates yield function for stresses and internal variables.
  
  @param[in] ipp - integration point pointer
  @param[in] idpm - id of plasticity model
  @param[in] sig - stress components
  @param[in] q - internal variables (hardening parameters)
   
  @return The function returns value of the yield function with respect to 
  the actual values of stresses and internal variables.

  Created by JK, TKo, 28.10.2001
*/
double mechmat::yieldfunction (long ipp, long idpm, vector &sig,vector &q)
{
  double f=0.0;
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    f = spl1d[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case jflow:{
    f = j2f[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case mohrcoulparab:{
    f = mcpar[ip[ipp].idm[idpm]].yieldfunction (sig);
    break;
  }
  case boermaterial:{
    f = boerm[ip[ipp].idm[idpm]].yieldfunction (sig);
    break;
  }
  case druckerprager:{
    f = drprm[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case doubledrprager:{
    f = ddpm[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case chenplast:{
    f = chplast[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case modcamclaymat:{
    f = cclay[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case modcamclaycoupmat:{
    f = cclayc[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case bbmcoupmat:{
    f = bbm[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case doublestructuremat:{
    f = dsm[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  case hissplasmat:{
    f = hisspl[ip[ipp].idm[idpm]].yieldfunction (sig,q);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  return f;
}



/**
  Function evaluates derivates of yield function with respect to stresses.
  Stresses must be assembled in tensor notation and stored in 3x3 %matrix.
  Function overwrites array sig, some material models use the array sig for auxiliary computations.
  Function stores derivatives in %vector dfds.

  @param[in]  ipp  - integration point pointer
  @param[in]  idpm - id of plasticity model
  @param[in]  sig  - stress components
  @param[in]  q    - internal parameters (hardening parameters)
  @param[out] dfds - derivatives of yield function with respect to stresses (output)
   
  @return The function returns resulting %vector of derivatives in the parameter dfds.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dfdsigma (long ipp, long idpm, vector &sig, vector &q, vector &dfds)
{
  vector dfdst(ASTCKVEC(6));

  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigma (sig,dfdst);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigma (sig,dfdst);
    break;
  }
  case mohrcoulparab:{
    mcpar[ip[ipp].idm[idpm]].deryieldfsigma (sig,dfdst);
    break;
  }
  case boermaterial:{
    boerm[ip[ipp].idm[idpm]].deryieldfsigma (sig,dfdst);
    break;
  }
  case druckerprager:{
    drprm[ip[ipp].idm[idpm]].dfdsigma (sig,dfdst,q);
    break;
  }
  case doubledrprager:{
    ddpm[ip[ipp].idm[idpm]].deryieldfsigma (sig,dfdst,q);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigma (sig,q,dfdst);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dfdst);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dfdst);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dfdst);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dfdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  
  //  conversion of full tensor notation to vector notation
  give_red_vector(dfdst,dfds,ip[ipp].ssst);
}



/**
  Function evaluates derivates of plastic potential with respect to stresses.
  Stresses must be assembled in tensor notation and stored in 3x3 %matrix.
  Function overwrites array sig, some material models use the array sig for auxiliary computations.
  Function stores derivatives in %vector dgds.

  @param[in]  ipp  - integration point pointer
  @param[in]  idpm - id of plasticity model
  @param[in]  sig  - stress components
  @param[in]  q    - internal parameters (hardening parameters)
  @param[out] dfds - derivatives of plastic potential with respect to stresses
   
  @return The function returns resulting %vector of derivatives in the parameter dgds.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dgdsigma (long ipp, long idpm, vector &sig, vector &q, vector &dgds)
{
  vector dgdst(ASTCKVEC(6));
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigma (sig,dgdst);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigma (sig,dgdst);
    break;
  }
  case mohrcoulparab:{
    mcpar[ip[ipp].idm[idpm]].derplaspotsigma (sig,dgdst);
    break;
  }
  case boermaterial:{
    boerm[ip[ipp].idm[idpm]].deryieldfsigma (sig,dgdst);
    break;
  }
  case druckerprager:{
    drprm[ip[ipp].idm[idpm]].dgdsigma (sig,dgdst,q);
    break;
  }
  case doubledrprager:{
    ddpm[ip[ipp].idm[idpm]].deryieldfsigma (sig,dgdst,q);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigma (sig,q,dgdst);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dgdst);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dgdst);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dgdst);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfsigma (sig, q, dgdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  
  //  conversion of full tensor notation to vector notation
  give_red_vector(dgdst,dgds,ip[ipp].ssst);
}



/**
  Function evaluates derivates of yield function with respect to internal variables (hardening parameters).
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  The function stores derivatives in %vector dq.
 
  @param[in]  ipp  - integration point pointer
  @param[in]  idpm - id of plasticity model
  @param[in]  sig  - stress components
  @param[in]  q    - internal variables (hardening parameters)
  @param[out] dq   - derivatives of yield function with respect to internal variables (output)
   
  @return The function returns resulting %vector of derivatives in the parameter dq.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dfdqpar (long ipp, long idpm, vector &sig,vector &q,vector &dq)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdqpar (dq);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdqpar (dq);
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    break;
  }
  case doubledrprager:{
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfq (sig, dq);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdq (sig,q,dq);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}



/**
  Function evaluates second derivates of yield function with respect to stresses.   
  Stresses must be assembled in tensor notation and stored in 3x3 %matrix,
  function stores derivatives in %matrix dfdsds.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - stress components
  @param[out] dfdsds - second derivatives of yield function with respect to stresses
   
  @return The function returns resulting %vector of derivatives in the parameter dfdsds.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dfdsigmadsigma (long ipp, long idpm, vector &sig, matrix &dfdsds)
{
  matrix dfdsdst(ASTCKMAT(6,6));
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigmadsigma (dfdsdst);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigmadsigma (sig,dfdsdst);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].dderyieldfsigma (dfdsdst);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].dderyieldfsigma (dfdsdst);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].dderyieldfsigma (dfdsdst);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].dderyieldfsigma (dfdsdst);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigmadsigma (sig,dfdsdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  
  //  conversion from the full tensor notation to matrix notation
  tensor4_matrix(dfdsds, dfdsdst, ip[ipp].ssst);
}



/**
  Function evaluates second derivates of the plastic potential function with respect to stresses.   
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  function stores derivatives in %matrix dfdsds.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - %vector stress components in the full format, i.e. 6 components
  @param[in]  q      - %vector of hardening variable values
  @param[out] dfdsds - second derivatives of yield function with respect to stresses in the reduced format, i.e. ncompstr x ncompstr
   
  @return The function returns resulting %vector of derivatives in the parameter dgdsds.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dgdsigmadsigma (long ipp, long idpm, vector &sig, vector &/*q*/, matrix &dgdsds)
{
  matrix dgdsdst(ASTCKMAT(6,6));
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigmadsigma (dgdsdst);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigmadsigma (sig,dgdsdst);
    break;
  }
  case druckerprager:{
    drprm[ip[ipp].idm[idpm]].dfdsigmadsigma (sig, dgdsdst);
    break; 
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].dderyieldfsigma (dgdsdst);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].dderyieldfsigma (dgdsdst);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].dderyieldfsigma (dgdsdst);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].dderyieldfsigma (dgdsdst);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigmadsigma (sig, dgdsdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  
  //  conversion from the full tensor notation to matrix notation
  tensor4_matrix(dgdsds, dgdsdst, ip[ipp].ssst);
}



/**
  Function evaluates derivates of yield function with respect to stresses and
  internal variables (hardening parameters).

  Stresses must be assembled in tensor notation and stored in 3x3 %matrix,
  function stores derivatives in %matrix dfdsdq.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - stress components
  @param[in]  q      - internal variables (hardening parameters)
  @param[out] dfdsdq - derivatives of yield function with respect to stresses and internal variables (output)
   
  @return The function returns resulting %matrix of derivatives in the parameter dfdsdq.

  Created by JK,TKo,  28.10.2001
*/
void mechmat::dfdsigmadqpar (long ipp, long idpm, vector &sig,vector &q,matrix &dfdsdq)
{
  //  number of hardening parameters
  long nh=dfdsdq.n;
  
  matrix dfdsdqt(ASTCKMAT(6,nh));

  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigmadq (dfdsdqt);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigmadq (dfdsdqt);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfdsdq (dfdsdqt);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfdsdq (dfdsdqt);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfdsdq (dfdsdqt);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfdsdq (dfdsdqt);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigmadq (sig,q,dfdsdqt);
    break;
  }
  default:
    print_err("unknown plasticity model is required", __FILE__,__LINE__,__func__);

  }
  
  //  conversion from full tensorial notation to reduced matrix notation
  red_rows_matrix(dfdsdq, dfdsdqt, ip[ipp].ssst);
}



/**
  Function evaluates second derivates of the plastic potential function with respect to stresses and hardening variables.   
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  function stores derivatives in %matrix dfdsds.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - %vector stress components in the full format, i.e. 6 components
  @param[in]  q      - %vector of hardening variable values
  @param[out] dgdsdq - second derivatives of plastic potential function with respect to stresses in the reduced format, i.e. ncompstr x q.n
   
  @return The function returns resulting %vector of derivatives in the parameter dfdsds.

  Created by TKo,  5.2022
*/
void mechmat::dgdsigmadqpar (long ipp, long idpm, vector &sig, vector &q, matrix &dgdsdq)
{
  matrix dgdsdqt(ASTCKMAT(q.n,6));
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigmadq (dgdsdqt);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigmadq (dgdsdqt);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfdsdq(dgdsdqt);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfdsdq(dgdsdqt);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfdsdq (dgdsdqt);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfdsdq (dgdsdqt);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdsigmadq (sig,q,dgdsdqt);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
    
  //  conversion from tensorial notation to matrix notation
  red_rows_matrix(dgdsdq, dgdsdqt, ip[ipp].ssst);
}



/**
  Function evaluates derivates of the hardening variables with respect to consistency parameter gamma.   
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  The function stores derivatives in %vector dqdg.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - %vector stress components in the full format, i.e. 6 components
  @param[in]  q      - %vector of hardening variable values
  @param[out] dqdg -   %vector of derivatives of hardeining variables with respect to consistency parameter, dimension must be ncomphard = q.n
   
  @return The function returns resulting %vector of derivatives in the parameter dfdsds.

  Created by TKo,  5.2022
*/
void mechmat::dqpardgamma (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, vector &dqdg)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    nullv(dqdg);
    break;
  }
  case jflow:{
    nullv(dqdg);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].der_q_gamma(ipp, ido, sig, q, dqdg);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].der_q_gamma(ipp, ido, sig, q, epsp, dqdg);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].der_q_gamma(ipp, ido, sig, q, epsp, dqdg);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].der_q_gamma(ipp, ido, sig, q, epsp, dqdg);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);
  }    
}



/**
  Function evaluates the second derivates of the hardening variables with respect to consistency parameter gamma and stresses.   
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  The function stores derivatives in %matrix dqdgds.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - %vector stress components in the full format, i.e. 6 components
  @param[in]  q      - %vector of hardening variable values
  @param[out] dqdgds - %matrix of the second derivatives of hardeining variables with respect to consistency parameter and stress components, 
                       dimension must be q.n x ncompstr (ncompstr means reduced stress %vector format)
   
  @return The function returns resulting %matrix of derivatives in the parameter dqdgds.

  Created by TKo,  5.2022
*/
void mechmat::dqpardgammadsigma (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, matrix &dqdgds)
{
  matrix dqdgdst(ASTCKMAT(q.n, 6));
  
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    nullm(dqdgds);
    break;
  }
  case jflow:{
    nullm(dqdgds);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].dqpardgammadsigma(ipp, ido, sig, q, epsp, dqdgdst);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].dqpardgammadsigma(ipp, ido, sig, q, epsp, dqdgdst);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].dqpardgammadsigma(ipp, ido, sig, q, epsp, dqdgdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }    

  //  conversion from the full tensor notation to the reduced matrix notation
  red_cols_matrix(dqdgds, dqdgdst, ip[ipp].ssst);
}



/**
  Function evaluates the second derivates of the hardening parameters with respect to consistency parameter gamma and hardening parameters.   
  Stresses tensor must be stored in the full %vector format, i.e. argument sig has dimension 6,
  The function stores derivatives in %matrix dqdgds.

  @param[in]  ipp    - integration point pointer
  @param[in]  idpm   - id of plasticity model
  @param[in]  sig    - %vector stress components in the full format, i.e. 6 components
  @param[in]  q      - %vector of hardening variable values
  @param[out] dqdgds - %matrix of the second derivatives of hardeining variables with respect to consistency parameter and stress components, 
                       dimension must be q.n x q.n
   
  @return The function returns resulting %matrix of derivatives in the parameter dqdgdq.

  Created by TKo,  5.2022
*/
void mechmat::dqpardgammadqpar (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, matrix &dqdgdq)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    nullm(dqdgdq);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].dqpardgammadqpar(ipp, ido, sig, q, epsp, dqdgdq);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].dqpardgammadqpar(ipp, ido, sig, q, epsp, dqdgdq);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].dqpardgammadqpar(ipp, ido, sig, q, epsp, dqdgdq);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }    
}



/**
  Function evaluates derivates of plastic potential with respect to internal variables (hardening parameters).
  Stresses must be assembled in tensor notation and stored in 3x3 %matrix,
  function stores derivatives in %vector dq.

  @param[in]  ipp  - integration point pointer
  @param[in]  idpm - id of plasticity model
  @param[in]  sig  - stress components
  @param[in]  q    - internal variables (hardening parameters)
  @param[out] dq   - derivatives of plastic potential with respect to internal variables (output)
   
  @return The function returns resulting %vector of derivatives in the parameter dq.

  Created by JK,TKo,  28.10.2001
*/
/*
void mechmat::dgdqpar (long ipp, long idpm, vector &sig, vector &q, vector &dq)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    //spl1d[ip[ipp].idm[idpm]].dfdqpar (dq);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].deryieldfq (dq);
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    break;
  }
  case doubledrprager:{
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfq (sig, dq);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].deryieldfq (sig, q, dq);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].deryieldfdq (sig,q,dq);
    break;
  }
  default:
    print_err ("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}
*/



/**
  Function computes derivatives of hardening function with respect to stresses.

  @param[in]  ipp  - id of integration point
  @param[in]  idpm - id of plasticity model
  @param[in]  sigt - stress components stored in 3x3 %matrix
  @param[in]  q    - internal variables stored in %vector
  @param[out] dhds - derivatives of hardening function with respect to stresses stored in 6 x ncomphard %matrix (output)
   
  @return The function returns resulting %matrix of derivatives in the parameter dhds.

  Created by JK, 7.8.2005
*/
void mechmat::dhdsigma (long ipp, long idpm, vector &/*sig*/, vector &/*q*/, matrix &dhds)
{
  //  number of hardening parameters
  long nh=dhds.m;
  
  matrix dhdst(nh,6);

  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdsigmadq (dhdst);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdsigmadq (dhdst);
    break;
  }
  case chenplast:{
    //chplast[ip[ipp].idm[idpm]].dhdsigma (sig,q,dhdst);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }

  //  conversion from tensorial notation to matrix notation
  red_cols_matrix (dhds,dhdst,ip[ipp].ssst);
}



/**
  Function computes derivatives of hardening function with respect to consistency parameter.

  @param[in]  ipp  - id of integration point
  @param[in]  idpm - id of plasticity model
  @param[out] dhdc - derivatives of hardening function with respect to consistency parameter stored in ncomphard x 1 %vector (output)
   
  @return The function returns resulting %vector of derivatives in the parameter dhdc.

  Created by JK, 7.8.2005
*/
/*
void mechmat::dhdgamma (long ipp,long idpm, vector &epsp, vector &sig,vector &q,vector &dhdc)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dhdgamma (ipp,epsp,sig,dhdc);
    break;
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].dhdgamma (dhdc);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}
*/
void mechmat::hardvect (long ipp,long idpm,vector &hv)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].hardvect (hv);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].hardvect (hv);
    break;
  }
  case druckerprager:{
    drprm[ip[ipp].idm[idpm]].hardvect (hv);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}



/**
  Function computes derivatives of hardening function with respect to hardening parameters.

  @param[in] ipp - id of integration point
  @param[in] idpm - id of plasticity model
  @param[in] sigt - stress components stored in %vector with 6 components
  @param[in] dgamma - increment of consistency parameter
  @param[in] q - internal variables stored in %vector
  @param[in] dhdq - derivatives of hardening function with respect to hardening parameters 
                stored in ncomphard x ncomphard %matrix (output)
   
  @return The function returns resulting %matrix of derivatives in the parameter dhdq.

  Created by JK, 7.8.2005
*/
void mechmat::dhdqpar (long ipp, long idpm, long /*ido*/, vector &/*sigt*/, double /*dgamma*/, vector &/*q*/, matrix &dhdq)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].dfdqpardqpar (dhdq);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].dfdqpardqpar (dhdq);
    break;
  }
  case druckerprager:{
    //drprm[ip[ipp].idm[idpm]].dfdqpardqpar (dhdq);
    break;
  }
  case chenplast:{
    //chplast[ip[ipp].idm[idpm]].dhdqpar (sigt,q,dhdq);
    break;
  }
  case modcamclaymat:{
    //cclay[ip[ipp].idm[idpm]].dhardfdq (ipp, ido, dgamma, q, dhdq);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}



/**
   Function computes hardening/softening function for required parameters.
   
   @param[in]  ipp  - integration point id
   @param[in]  idpm - id of plastic material
   @param[in]  sigt - stress components stored in 3x3 %matrix
   @param[in]  q    - internal parameters
   @param[out] h    - values of hardening function stored in %vector
   
   @return The function returns values of hardening/softening function in the parameter h.

   JK, 3.2.2007
*/
void mechmat::hardsoftfunction (long ipp,long idpm, vector &sig, vector &q, vector &h)
{
  switch (ip[ipp].tm[idpm]){
  case druckerprager:{
    break;  
  }
  case chenplast:{
    chplast[ip[ipp].idm[idpm]].hvalues (sig,q,h);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }  
}



/**
  The function computes matrix of plastic moduli with respect of actual values of 
  stresses and internal variables.

  @param[in]  ipp  - id of integration point
  @param[in]  idpm - id of plasticity model
  @param[in]  sigt - stress components stored in 3x3 %matrix
  @param[in]  epsp - %vetcor of plastic strains
  @param[out] h    -  ncomphard x ncomphard %matrix for plastic moduli (output)
   
  @return The function returns resulting %matrix of moduli in the parameter h.

  Created by JK, 7.8.2005
*/
void mechmat::plasticmoduli (long ipp,long idpm,vector &/*epsp*/,vector &/*sig*/,matrix &h)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    spl1d[ip[ipp].idm[idpm]].plasmod (h);
    break;
  }
  case jflow:{
    //j2f[ip[ipp].idm[idpm]].plasmod (h);
    break;
  }
  case chenplast:{
    //chplast[ip[ipp].idm[idpm]].plasmod (h);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }

}



/*
void mechmat::dfunctdqdq (long ipp, long idpm,matrix &sig,vector &q,matrix &dfdqdq)
{
  switch (ip[ipp].tm[idpm]){
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].deryieldfdqdq (dfdqdq);
    break;
  }
  case chenplast:{
    //chplast[ip[ipp].idm[idpm]].deryieldfdqdq (dfdqdq);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown plasticity model is required in function dfdqdq (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}
*/



/**
  Function assembleses contribution from hardening to the denominator in cutting plane method.

  @param[out] denom - resulting contribution from hardening
  @param[in]  ipp   - integration point id
  @param[in]  idpm  - id of plastic material
  @param[in]  ido   - index of internal variables for given material in the ipp other array
  @param[in]  gamma - consistency parameter
  @param[in]  sig   - %vector of actial stresses
  @param[in]  eps   - %vector of actual strains
  @param[in]  epsp  - %vector of actual plastic strains
  @param[in]  qtr   - %vector of trial hardening parameters

  @return The function returns computed contribution to denominator in the parameter denom.

  Created by JK, TKo
*/
void mechmat::plasmodscalar(double &denom,long ipp, long idpm, long ido, vector &sig, vector &/*eps*/, vector &epsp, vector &qtr,double gamma)
{
  double ret = 0.0;
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    //ret = spl1d[ip[ipp].idm[idpm]].plasmodscalar (qtr);
    ret = spl1d[ip[ipp].idm[idpm]].plasmodscalar (sig,epsp,qtr,gamma);
    break;
  }
  case jflow:{
    ret = j2f[ip[ipp].idm[idpm]].plasmodscalar ();
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    ret = drprm[ip[ipp].idm[idpm]].plasmodscalar(qtr);
    break;
  }
  case doubledrprager:{
    ret = ddpm[ip[ipp].idm[idpm]].plasmodscalar(sig, qtr);
    break;
  }
  case chenplast:{
    //ret = chplast[ip[ipp].idm[idpm]].plasmodscalar (ipp,sig,epsp,qtr);
    if (ret<1.0e-6)
      denom*=2.0;
    //else
    //if (epsp[0]<1.0e-6)
    //denom*=2.0;
    //else
    //denom+=ret;
    return;
    //break;
  }
  case modcamclaymat:{
    ret = cclay[ip[ipp].idm[idpm]].plasmodscalar(ipp, ido, sig, qtr);
    break;
  }
  case modcamclaycoupmat:{
    ret = cclayc[ip[ipp].idm[idpm]].plasmodscalar(ipp, ido, sig, qtr, epsp);
    break;
  }
  case bbmcoupmat:{
    ret = bbm[ip[ipp].idm[idpm]].plasmodscalar(ipp, ido, sig, qtr, epsp);
    break;
  }
  case doublestructuremat:{
    ret = dsm[ip[ipp].idm[idpm]].plasmodscalar(ipp, ido, sig, qtr, epsp);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  denom = ret;
//  return ret;
}



/**
  Function returns number of internal parameters for material models of plasticity.

  @param[in] ipp  - integration point number in the mechmat ip array.
  @param[in] idpm - index of material which is plasticty.

  @return Function returns number of material parameters. 

  Created by JK,
*/
long mechmat::give_num_interparam (long ipp,long idpm)
{
  long np=0;

  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    np = spl1d[ip[ipp].idm[idpm]].give_num_interparam ();
    break;
  }
  case jflow:{
    np = j2f[ip[ipp].idm[idpm]].give_num_interparam ();
    break;
  }
  case mohrcoul:{
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    break;
  }
  case modcamclaymat:{
    break;
  }
  case modcamclaycoupmat:{
    break;
  }
  case bbmcoupmat:
  case doublestructuremat:{
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
  
  return np;
}



/**
  Function gives internal parameters of the required material model.
   
  @param[in]  ipp - number of integration point
  @param[in]  im  - index of the material
  @param[in]  ido - index in the eqother array
  @param[out] q   - array containing internal parameters
    
  @return The function returns %vector of internal parameters in the parameter q. 

  Created by JK,
*/
void mechmat::give_interparam (long ipp,long im,long ido,vector &q)
{
  switch (ip[ipp].tm[im]){
  case simplas1d:{
    spl1d[ip[ipp].idm[im]].give_interparam (ipp,ido,q);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[im]].give_interparam (ipp,ido,q);
    break;
  }
  case mohrcoul:{
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    break;
  }
  case modcamclaymat:{
    break;
  }
  case modcamclaycoupmat:{
    break;
  }
  case bbmcoupmat:
  case doublestructuremat:{
    break;
  }
    /*
  case lemaitr:{
    break;
  }
  */
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);

  }
  
}



/**
  Function computes (updates) new hardening parameters depending on the actual 
  values of the internal parameters (plastic strains) of the given plasticity model.

  @param[in]  ipp - integration point id
  @param[in]  idpm - id of plastic material
  @param[in]  ido - index of internal variables for given material in the ipp other array
  @param[in]  dgamma - increment of consistency parameter
  @param[in]  eps  - %vector of actual strains
  @param[in]  epsp - %vector of actual plastic strains
  @param[in]  sig  - %vector of actual stresses
  @param[out] q    - %vector of resulting hardening parameters

  @return The function returns %vetcor of updated hardening parameters in the parameter q.

  Created by JK,
*/
void mechmat::updateq (long ipp, long idpm, long ido,double dgamma, vector &eps, vector &epsp, vector &sig, vector &q)
{
  switch (ip[ipp].tm[idpm]){
  case simplas1d:{
    //spl1d[ip[ipp].idm[idpm]].updateq(dgamma, q);
    spl1d[ip[ipp].idm[idpm]].updateq (ipp,dgamma, epsp,q);
    break;
  }
  case jflow:{
    j2f[ip[ipp].idm[idpm]].updateq (dgamma, q);
    break;
  }
  case mohrcoul:{
    break;
  }
  case mohrcoulparab:{
    break;
  }
  case boermaterial:{
    break;
  }
  case druckerprager:{
    drprm[ip[ipp].idm[idpm]].updateq(ipp, epsp, q);
    break;
  }
  case doubledrprager:{
    ddpm[ip[ipp].idm[idpm]].updateq(ipp, dgamma, q);
    break;
  }
  case chenplast:{
    //chplast[ip[ipp].idm[idpm]].updateq (ipp,dgamma,epsp,sig,q);
    break;
  }
  case modcamclaymat:{
    cclay[ip[ipp].idm[idpm]].updateq(ipp, ido, eps, epsp, sig, q);
    break;
  }
  case modcamclaycoupmat:{
    cclayc[ip[ipp].idm[idpm]].updateq(ipp, ido, eps, epsp, sig, q);
    break;
  }
  case bbmcoupmat:{
    bbm[ip[ipp].idm[idpm]].updateq(ipp, ido, eps, epsp, sig, q);
    break;
  }
  case doublestructuremat:{
    dsm[ip[ipp].idm[idpm]].updateq(ipp, ido, eps, epsp, sig, q);
    break;
  }
  default:
    print_err("unknown plasticity model is required",__FILE__,__LINE__,__func__);

  }
}



/**
  Function returns stresses on the yield sufrace.
  Cutting plane method is used.

  gamma,epsp and q will be replaced by new values

  @param[in]     ipp   - integration point pointer
  @param[in]     im    - material index
  @param[in]     ido   - index of internal variables for given material in the ipp other array
  @param[in/out] gamma - consistency parameter
  @param[in]     epsn  - total strain components
  @param[in/out] epsp  - plastic strain components
  @param[in/out] q     - hardening parameters
  @param[in]     ni    - maximum number of iterations
  @param[in]     err   - required error
  
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual consistency parameter in the parameter gamma and actual %vector of
          hardening parameters in the parameter q.
 
  Created by JK, 4.8.2001
  Modified by TKo
*/
void mechmat::cutting_plane (long ipp, long im, long ido, double &gamma, vector &epsn, vector &epsp, vector &q, long ni, double err)
{
  long i,j,ncomp=epsn.n;
  double f,denom, denomh, dgamma, nordsig, norsig, rnorsig;
  vector epsa(ASTCKVEC(ncomp)),sig(ASTCKVEC(ncomp)),dfds(ASTCKVEC(ncomp)),dgds(ASTCKVEC(ncomp));
  vector sigt(ASTCKVEC(6)), sigp(ASTCKVEC(ncomp)), dsig(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp));

  //  initialization
  dgamma=0.0;

  //  elastic stiffness matrix
  elmatstiff (d,ipp);
  if (ip[ipp].ssst == planestrain)
  {
     d[0][3] = d[0][1]; d[1][3] = d[1][0];
     d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
/*  if (ip[ipp].ssst == planestress)
  {
    nu = give_actual_nu(ipp);
    epsn[3] = -nu / (1.0 - nu) * (epsn[0]-epsp[0]+epsn[1]-epsp[1]);
  }
*/  
  
  //fprintf (Out,"\n\n int. point %ld \n",ipp);


  //  main iteration loop
  for (i=0;i<ni;i++){

    //fprintf (Out,"\n cutting plane, iterace  %ld",i);

    // save stresses atatined in the previous cutting plane step
    copyv(sig, sigp);
    //  elastic strain
    subv (epsn,epsp,epsa);
    //  trial stress computation
    mxv (d,epsa,sig);
    // add intial stresses (eigenstresses) if defined
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
      for (j=0; j<ncomp; j++)
        sig(j) += eigstresses[ipp][j];
    }
    
    give_full_vector (sigt,sig,ip[ipp].ssst);
    // updating hardening variables from previous step
    updateq(ipp, im, ido, dgamma, epsn, epsp, sig, q);
    // checking yield function
    f = yieldfunction (ipp,im,sigt,q);

    // calcuate relative norm of stress increment from the two successive cutting plane steps,
    // for i == 0, the relative norm should be zero and thus only actual yield function value
    // contributes in the tolerance checking
   
    // calculate stress increment from the two successive cutting plane steps
    if (i > 0)
      subv(sig, sigp, dsig);
    else
      // initial value of the reference stress used in the stress norm
      copyv(sig, sigp);
    norsig  = normv(sigp);
    nordsig = normv(dsig);
    rnorsig = nordsig/norsig;
    if ((f<err) && (rnorsig < 1.0e-6))
      break;
    if (i==ni-1 && f>err){
      print_err("yield surface was not reached on element %ld  (ip=%ld)\n  f=%le rnorsig=%le  ni %ld",__FILE__,__LINE__,__func__, elip[ipp]+1, ipp+1, f, rnorsig, i);
      fprintf (Out,"\n yield surface was not reached");
      fflush(Out);
    }

    dfdsigma (ipp, im, sigt, q, dfds);
    dgdsigma (ipp, im, sigt, q, dgds);


    mxv (d,dgds,epsa);
    scprd (dfds,epsa,denom);

    denomh=0;
    plasmodscalar(denomh,ipp, im, ido, sig, epsn, epsp, q,gamma);
    denom += denomh;

    //  new increment of consistency parameter
    dgamma = f/denom;
    //  new internal variables
    gamma+=dgamma;
    // updating plastic strains
    for (j=0;j<ncomp;j++){
      epsp[j]+=dgamma*dgds[j];
    }
/*    if (ip[ipp].ssst == planestress)
      epsp[3] = epsp[0]+epsp[1]; */
  }
  
  //  elastic strain
  subv (epsn,epsp,epsa);
  //  stress computation
  mxv (d,epsa,sig);
  // add initial stresses (eigenstresses) if defined
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (j=0; j<ncomp; j++)
      sig(j) += eigstresses[ipp][j];
  }
  //  storing resulting stresses
  storestress (0,ipp,sig);
}



/**
  Function returns stresses on the yield sufrace.
  Cutting plane method is used.

  gamma,epsp and q will be replaced by new values

  @param[in]     ipp   - integration point pointer
  @param[in]     im    - material index
  @param[in]     ido   - index of internal variables for given material in the ipp other array
  @param[in/out] gamma - consistency parameter
  @param[in]     epsn  - total strain components
  @param[in/out] epsp  - plastic strain components
  @param[in/out] q     - hardening parameters
  @param[in]     ni    - maximum number of iterations
  @param[in]     err   - required error
  @param[in]     dcomp - flag for evaluation algorithmic consistent stiffness matrix
  @param[out]    cd    - algorithmic consistent stiffness matrix, it must be allocated to dimensions (ncompstr,ncompstr)
  
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual consistency parameter in the parameter gamma and actual %vector of
          hardening parameters in the parameter q and if the dcomp=1 cd argument contains actual 
          algoritmic stiffness %matrix.
  @retval 0 - successful stress return with the prescribed tolerances and number of steps
  @retval 1 - stress return was not successful, given tolerances cannot be attained in the required number of steps
 
  Created by TKo, 5.2022
*/
long mechmat::cutting_plane (long ipp, long im, long ido, double &gamma, vector &epsn, vector &epsp, vector &q,
                             long ni, double err, long dcomp, matrix &cd)
{
  long i,j,ncomp=epsn.n;
  long ncsh = ncomp + q.n; // dimension of extended vector
  double f,denom, denomh, dgamma, nordsig, norsig, rnorsig, rnorf;
  vector epsa(ASTCKVEC(ncomp)),sig(ASTCKVEC(ncomp)),dfds(ASTCKVEC(ncomp)),dgds(ASTCKVEC(ncomp));
  vector sigt(ASTCKVEC(6)), sigp(ASTCKVEC(ncomp)), dsig(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp));
  // auxiliary vectors and matrices used for the calculation of algorithmic consistent stiffness matrix
  matrix ecd, drds, auxm, dauxm;
  vector r, dfdq, dqdg, edfdso, edfdsn, eauxv, dgammadeps;
  ivector id;
  double denomdg;
  strastrestate ssst = ip[ipp].ssst;
  long ret = 1;

  //  initialization
  dgamma=0.0;

  //  elastic stiffness matrix
  elmatstiff (d,ipp);

  // initial value of algorithmic consistent stiffness matrix 
  if (dcomp){
    reallocv(RSTCKVEC(ncsh, r));
    reallocv(RSTCKVEC(q.n, dqdg));
    reallocv(RSTCKVEC(q.n, dfdq));
    reallocv(RSTCKVEC(ncsh, edfdso));
    reallocv(RSTCKVEC(ncsh, edfdsn));
    reallocv(RSTCKVEC(ncomp, dgammadeps));
    
    reallocm(RSTCKMAT(ncsh, ncsh, drds));
    reallocm(RSTCKMAT(ncsh, ncomp, ecd));
    storeblock(ecd, d, 0, 0);
    reallocm(RSTCKMAT(ncsh, ncsh, auxm));
    reallocm(RSTCKMAT(ncsh, ncsh, dauxm));
  }
  
  if (ip[ipp].ssst == planestrain)
  {
     d[0][3] = d[0][1]; d[1][3] = d[1][0];
     d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
/*  if (ip[ipp].ssst == planestress)
  {
    nu = give_actual_nu(ipp);
    epsn[3] = -nu / (1.0 - nu) * (epsn[0]-epsp[0]+epsn[1]-epsp[1]);
  }
*/  
  
  //fprintf (Out,"\n\n int. point %ld \n",ipp);


  //  main iteration loop
  for (i=0;i<ni;i++){

    //fprintf (Out,"\n cutting plane, iterace  %ld",i);

    // save stresses atatined in the previous cutting plane step
    copyv(sig, sigp);
    if ((i == 0) || (dcomp == 0)){
      //  elastic strain
      subv(epsn, epsp, epsa);
      //  trial stress computation
      mxv(d, epsa, sig);
      // add intial stresses (eigenstresses) if defined
      if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
        for (j=0; j<ncomp; j++)
          sig(j) += eigstresses[ipp][j];
      }
    
      give_full_vector(sigt ,sig, ip[ipp].ssst);
      // updating hardening variables from previous step
      updateq(ipp, im, ido, dgamma, epsn, epsp, sig, q);
    }
    // checking yield function
    f = yieldfunction(ipp, im, sigt, q);

    // calcuate relative norm of stress increment from the two successive cutting plane steps,
    // for i == 0, the relative norm should be zero and thus only actual yield function value
    // contributes in the tolerance checking
    //
    // calculate stress increment from the two successive cutting plane steps
    if (i > 0)
      subv(sig, sigp, dsig);
    else
      // initial value of the reference stress used in the stress norm
      copyv(sig, sigp);
    norsig  = normv(sigp);
    nordsig = normv(dsig);
    rnorsig = nordsig/norsig;
    if (i == 0){
      if ((f<err) && (rnorsig < 1.0e-6)){
        ret = 0;
        break;
      }
    }
    else{
      rnorf = sqrt(fabs(f));
      if (norsig > 0.0)
        rnorf /= norsig;
      if ((rnorf < err) && (rnorsig < 1.0e-6)){
        ret = 0;
        break;
      }
    }
    if (i==ni-1 && f>err){
      
      print_err("yield surface was not reached on element %ld  (ip=%ld)\n  f=%le rnorsig=%le  ni %ld",__FILE__,__LINE__,__func__, elip[ipp]+1, ipp+1, f, rnorsig, i);
      fprintf (Out,"\n yield surface was not reached");
      fflush(Out);
    }

    if ((i == 0) || (dcomp == 0)){
      dfdsigma(ipp, im, sigt, q, dfds);
    }
    dgdsigma(ipp, im, sigt, q, dgds);


    mxv(d, dgds, epsa);
    scprd(dfds, epsa, denom);

    denomh=0.0;
    if (dcomp){
      // dfdq and dqdg vectors will be used later in the assembling of the algorithmic stiffness matrix
      // but the computation is equal to the one in plasmodscalar function
      dfdqpar(ipp, ido, sig, q, dfdq);
      dqpardgamma(ipp, ido, im, sig, q, epsp, dqdg);
      scprd(dfdq, dqdg, denomh);
      denomh = -denomh;
    }
    else
      plasmodscalar(denomh, ipp, im, ido, sig, epsn, epsp, q, gamma);

    denom += denomh;
    
    //  new increment of consistency parameter
    dgamma = f/denom;
    //  new internal variables
    gamma+=dgamma;
    // updating plastic strains
    for (j=0;j<ncomp;j++){
      epsp[j]+=dgamma*dgds[j];
    }
    if (dcomp){
      // assemble vector R^T = {D_e*dg/dsig, dq/dgamma} from the values of the actual CPA step
      makerefv(eauxv, r.a, ncomp);
      mxv(d, dgds, eauxv); // shear components of dgds vector should be doubled due to tensor->vector transformation
      makerefv(eauxv, r.a+ncomp, r.n-ncomp);
      extract(eauxv, dqdg, 0, dqdg.n);

      // assemble vector (d\Phi/dSig)^T = {df/dsig, df/dq} from the values of the actual CPA step
      extract(edfdso, dfds, 0, ncomp);
      makerefv(eauxv, edfdso.a+ncomp, edfdso.n-ncomp);
      extract(eauxv, dfdq, 0, dfdq.n);

      //                           / D_e*d(dg/dsig)/dsig,  D_e*d(dg/dsig)/dq \  from values 
      // assemble matrix dR/dSig = \ d(dq/dgamma)/dsig,    d(dq/dgamma)/dq   /  of the actual step
      //
      // compute block D_e*d(dg/dsig)/dsig of dR/dSig
      reallocm(RSTCKMAT(ncomp, ncomp, dauxm));
      dgdsigmadsigma(ipp, im, sigt, q, dauxm); // shear component block of ddgdsig matrix should be doubled due to tensor->vector transformation
      reallocm(RSTCKMAT(ncomp, ncomp, auxm));
      mxm(d, dauxm, auxm);
      storeblock(drds, auxm, 0, 0);
      // compute block D_e*d(dg/dsig)/dq of dR/dSig
      reallocm(RSTCKMAT(ncomp, q.n, dauxm));
      dgdsigmadqpar(ipp, im, sigt, q, dauxm); // shear components of d(dg/dsig)/dq matrix should be doubled due to tensor->vector transformation
      reallocm(RSTCKMAT(ncomp, q.n, auxm));
      mxm(d, dauxm, auxm);
      storeblock(drds, auxm, 0, ncomp);
      // compute block d(dq/dgamma)/dsig of dR/dSig
      reallocm(RSTCKMAT(q.n, ncomp, dauxm));
      dqpardgammadsigma (ipp, ido, im, sigt, q, epsp, dauxm);
      storeblock(drds, dauxm, ncomp, 0);
      // compute block d(dkappa/dgamma)/dkappa of dR/dSig
      reallocm(RSTCKMAT(q.n, q.n, dauxm));
      dqpardgammadqpar(ipp, ido, im, sigt, q, epsp, dauxm);
      storeblock(drds, dauxm, ncomp, ncomp);

      //
      // actualize stresses and hardening variables for new CPA step
      //
      
      //  elastic strain
      subv(epsn, epsp, epsa);
      //  trial stress computation
      mxv(d, epsa, sig);
      // add intial stresses (eigenstresses) if defined
      if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
        for (j=0; j<ncomp; j++)
          sig(j) += eigstresses[ipp][j];
      }
      give_full_vector(sigt, sig, ip[ipp].ssst);
      // updating hardening variables from previous step
      updateq(ipp, im, ido, dgamma, epsn, epsp, sig, q);

      //
      // assemble vector (d\Phi/dSig)^T = {df/dsig, df/dq} from the values of the new CPA step
      //
      // compute actual values of derivatives of yield function
      dfdsigma (ipp, im, sigt, q, dfds);
      dfdqpar(ipp, ido, sig, q, dfdq);
      extract(edfdsn, dfds, 0, ncomp); // shear components of dgds vector should be doubled due to tensor->vector transformation
      makerefv(eauxv, edfdsn.a+ncomp, edfdsn.n-ncomp);
      extract(eauxv, dfdq, 0, q.n);

      // compute d(dgamma)/d(Deps) = [d\Phi/dSig|(k+1) - dgamma.d\Phi/dSig|(k)*dR/dSig|(k)]/(d\Phi/dSig|(k)*R(k))
      reallocv(RSTCKVEC(ncsh, eauxv));
      vxm(edfdso, drds, eauxv);
      cmulv(dgamma, eauxv);
      // auxiliary vector eauxv should have doubled shear components due to tensor->vector transformation
      reallocv(RSTCKIVEC(3, id));
      give_shear_indices(ssst, id);
      for(long k=0; k<id.n; k++)
        eauxv(id(k)) *= 2.0;
      subv(edfdsn, eauxv, eauxv);
      vxm(eauxv, ecd, dgammadeps);
      scprd(edfdso, r, denomdg);
      cmulv(1.0/denomdg, dgammadeps);

      // compute new value of algorithmic stiffness matrix
      cmulm(dgamma, drds);
      // dR/dSig should have doubled shear components due to tensor->vector transformation
      for(long k=0; k<ncsh; k++){
        for(long l=0; l<id.n; l++)
          drds(k,id(l)) *= 2.0;
      }
      // compute dR/dSig*D_c|(k), where D_c|(k) = dSig/dDeps|(k) is the algorithmic stiffness matrix from the k-th step of CPA
      reallocm(RSTCKMAT(ncsh, ncomp, auxm));
      mxm(drds, ecd, auxm);
      // compute D_c|(k) - dR/dSig*D_c|(k)
      subm(ecd, auxm, ecd);
      // compute R|(k)*d(dgamma)/d(Deps)
      vxv(r, dgammadeps, auxm);
      // compute D_c - dR/dSig*D_c|(k) - R|(k)*d(dgamma)/d(Deps)
      subm(ecd, auxm, ecd);
    }
  }

  //  elastic strain
  subv(epsn, epsp, epsa);
  //  stress computation
  mxv(d, epsa, sig);
  // add initial stresses (eigenstresses) if defined
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (j=0; j<ncomp; j++)
      sig(j) += eigstresses[ipp][j];
  }
  //  storing resulting stresses
  storestress(0, ipp, sig);
  if (dcomp){    
    extractm(cd , ecd, 0, ncomp);
  }
    

  if (test_math_errip(Mm->elip[ipp], ipp))
    ret = 1;

  return ret;
}



/**
  Function returns stresses on sufrace of plasticity.
  Cutting plane method is used.
  Special approach for plane stress is implemented.

  gamma,epsp and q will be replaced by new values

  @param[in]     ipp   - integration point pointer
  @param[in]     im    - material index
  @param[in]     ido   - index of internal variables for given material in the ipp other array
  @param[in/out] gamma - consistency parameter
  @param[in]     epsn  - total strain components
  @param[in/out] epsp  - plastic strain components
  @param[in/out] q     - hardening parameters
  @param[in]     ni    - maximum number of iterations
  @param[in]     err   - required error
  
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual consistency parameter in the parameter gamma and actual %vector of
          hardening parameters in the parameter q.
 
  2/2015  J. Fiedler
*/

void mechmat::cutting_plane3 (long ipp,long im,long ido,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err)
{
  long j,n=ip[ipp].ncompstr;
  matrix d(n,n);
  double pom,nu;
  
  if (ip[ipp].ssst == planestress)
  {
    ip[ipp].ssst = axisymm;
    epsn[3] = epsn[2];          // switching total strain components due to change to axisym problem
    
    pom     = epsp[3];          // switching platic strain components due to change to axisym problem
    epsp[3] = epsp[2];
    epsp[2] = pom;
      
    nu = give_actual_nu(ipp);
    epsn[2] = -nu/(1-nu)*(epsn[0]-epsp[0]+epsn[1]-epsp[1])+epsp[2];     // trial transverse strain

    // iteration loop for transverse stress correction
    for (j=0;j<ni;j++)
    {
      cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
      if (abs(ip[ipp].stress[2])<err)      // transverse stress testing
        break;
      if (j==ni-1 && abs(ip[ipp].stress[2])>err)
      {
        print_err("transverse stress correction was not successful",__FILE__,__LINE__,__func__);
        fprintf (Out,"\n transverse stress correction was not successful");
      }
      elmatstiff (d,ipp);
      epsn[2]+= -ip[ipp].stress[2]/d(2,2);      // correction of the transverse strain
    }

    pom               = ip[ipp].stress[3];  // switching stress components according to original plane stress problem
    ip[ipp].stress[3] = ip[ipp].stress[2];
    ip[ipp].stress[2] = pom;

    pom     = epsp[3];      // switching plastic strain components according to original plane stress problem
    epsp[3] = epsp[2];
    epsp[2] = pom;

    ip[ipp].ssst = planestress;
  }
  else
    cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
}


/**
  Function returns stresses on surface of plasticity in multi-surface
  plasticity with hardening is assumed. Cutting plane method is used.   
  Vectors epsp,gamma and q will be overwritten by new values.
   
  @param[in]     ipp   - integration point integer
  @param[in]     epsn  - new total strain components
  @param[in/out] epsp  - plastic strain components
  @param[in/out] gamma - %vector of consistency parameters
  @param[in/out] q     - %vector of hardening parameters
   
  ncomp - number of components of strain tensors
  nyf - number of yield functions
  nhp - number of hardening parameters

  
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual %vector of consistency parameters in the parameter gamma and actual %vector of
          hardening parameters in the parameter q.
   
  Created by JK, 16.8.2001
*/
void mechmat::mult_surf_cutting_plane (long ipp,vector &epsn,vector &epsp,vector &gamma,vector &q)
{
  long ni,i,j,k,nas,ncomp=epsn.n,nhp=q.n,nyf=gamma.n;
  double err;
  vector epsptr(ASTCKVEC(ncomp)),epsa(ASTCKVEC(ncomp)),sig(ASTCKVEC(ncomp)),f(ASTCKVEC(nyf)),qtr(ASTCKVEC(nhp));
  vector dq,dgamma,ff;
  matrix d(ASTCKMAT(ncomp,ncomp)),h(ASTCKMAT(nhp,nhp));
  matrix fsig,fq,cpm,amq,am,acpm;
  ivector stat(ASTCKIVEC(nyf));
  
  //ni = Mp->nicp;
  //err = Mp->errcp;
  ni=0;
  err=10.0;
  
  //  initialization
  copyv (epsp,epsptr);
  copyv (q,qtr);
   
  //  main iteration loop
  for (i=0;i<ni;i++){

    subv (epsn,epsptr,epsa);
    
    //  elastic stiffness matrix
    elmatstiff (d,ipp);
    
    //  trial stress computation
    mxv (d,epsa,sig);
    // add initial stresses (eigenstresses) if defined
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
      for (j=0; j<ncomp; j++)
        sig(j) += eigstresses[ipp][j];
    }
    


    /*
    switch (ip[ipp].tm[0]){
      //case mohrcoul3d:{
      //mc3d[ip[ipp].idm[0]].yieldfunction (sig,f);
      //break;
    }
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function mult_surf_cutting_plane (mechmat.cpp).\n");
    }
    }
    */

    //  yield function control
    nas=0;
    for (j=0;j<nyf;j++){
      if (f[j]<err){      stat[j]=0;  }
      else{  nas++;       stat[j]=1;  }
    }
    if (nas==0)  break;

    
    reallocm (nas,6,fsig);  reallocm (nas,6,am);  reallocm (nas,nas,cpm);
    reallocm (nas,nhp,fq);  reallocm (nhp,nas,amq);  reallocm (nas,nas,acpm);
    reallocv (nas,ff);  reallocv (nas,dgamma);  reallocv (nas,dq);
    
    /*
    switch (ip[ipp].tm[0]){

	case mohrcoul3d:{
        k=0;
	for (j=0;j<6;j++){
        if (stat[j]==1){
	mc3d[ip[ipp].idm].deryieldfsigma (sig,j,fsig[k]);
	ff[k]=f[j];
	k++;
        }
	}
	break;

    }
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function mult_surf_cutting_plane (mechmat.cpp).\n");
    }
    }
    */
    
    mxm (fsig,d,am);
    mxmt (am,fsig,cpm);
    

    mxmt (h,fq,amq);
    mxm (fq,amq,acpm);

    addm (cpm,acpm,cpm);
    
    //  solution of auxiliary system of equations
    gemp (cpm.a,dgamma.a,ff.a,nas,1,Mp->zero,1);

    //  new internal variables
    k=0;
    for (j=0;j<nyf;j++){
      if (stat[j]==1){
        gamma[j]+=dgamma[k];
        k++;
      }
    }
    
    mtxv (fsig,dgamma,sig);
    mxv (amq,dgamma,dq);

    //  new plastic strains
    for (j=0;j<ncomp;j++){
      epsptr[j]+=sig[j];
    }
    //  new hardening parameters
    for (j=0;j<nhp;j++){
      qtr[j]-=dq[j];
    }
  }

  copyv (epsptr,epsp);
  copyv (qtr,q);


  //  correct elastic strain
  subv (epsn,epsp,epsa);

  //  elastic stiffness matrix
  elmatstiff (d,ipp);
  
  //  correct stress components
  mxv (d,epsa,sig);
  // add initial stresses (eigenstresses) if defined
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (j=0; j<ncomp; j++)
      sig(j) += eigstresses[ipp][j];
  }

  storestress (0,ipp,sig);

}



/**
  Function returns stresses on sufrace of plasticity.
  Newton method is used.

  gamma,epsp and q will be replaced by new values

  @param[in]     ipp   - integration point id
  @param[in]     im    - id of material type (for combination of material models, default value is 0)
  @param[in]     ido   - first index in eqother array (for combination of material models, default value is 0)
  @param[in/out] gamma - consistency parameter
  @param[in]     epsn  - total strain components
  @param[in/out] epsp  - plastic strain components
  @param[in/out] q     - hardening parameters
  @param[in]     ni    - maximum number of iterations
  @param[in]     err   - required error
   
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual consistency parameter in the parameter gamma and actual %vector of
          hardening parameters in the parameter q.

  JK, 16.3.2005, revised 7.8.2005, 17. 6. 2015
*/
void mechmat::newton_stress_return (long ipp,long im,long ido,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err)
{
  long i,j,k,n,ncompstr,ncomphard;
  double zero,f,dgamma,norstr,norstrincr,norconspar,norhardpar,norhardparincr;
  strastrestate ssst;
  
  //  computer zero
  zero=Mp->zero;

  //  number of strain components
  ncompstr = epsn.n;
  //  number of hardening components
  ncomphard = q.n;
  //  number of unknowns in the Newton method
  n = ncompstr + 1 + ncomphard;
  //  unknowns in the Newton-Raphson method are organized as follows:
  //  increments of stresses (ncompstr components)
  //  increment of consistency parameter (1 component)
  //  increments of hardening/softening parameters (ncomphard components)
  //
  //  hardening parameters may be omitted
  
  //  strain/stress state in actual integration point
  ssst = ip[ipp].ssst;
  
  //  array containing matrix of the system (Jacobian matrix)
  matrix a(ASTCKMAT(n,n));
  //  array containing solution of the system
  vector x(ASTCKVEC(n));
  //  array containing right hand side
  vector y(ASTCKVEC(n));
  //  stiffness matrix of the material
  matrix d(ASTCKMAT(ncompstr,ncompstr));
  //  compliance matrix of the material
  matrix c(ASTCKMAT(ncompstr,ncompstr));
  //  vector of elastic strains
  vector epse(ASTCKVEC(ncompstr));
  //  vector of backup of plastic strains
  vector epspold(ASTCKVEC(ncompstr));
  //  vector of stresses (in integration point)
  vector sig(ASTCKVEC(ncompstr));
  //  stress tensor stored in a vector with 6 components
  vector sigt(ASTCKVEC(6));
  //  backup of the trial stresses
  vector sigtrial(ASTCKVEC(ncompstr));
  //  backup of attained hardening parameters
  vector qold(ASTCKVEC(ncomphard));
  //  matrix of second derivatives of plastic potential with respect to stresses
  matrix dgdsds(ASTCKMAT(ncompstr,ncompstr));
  //  vector of derivatives of plastic potential with respect to stresses
  vector dgds(ASTCKVEC(ncompstr));
  //  vector of derivatives of yield function with respect to stresses
  vector dfds(ASTCKVEC(ncompstr));
  //  matrix of second derivatives of plastic potential with respect to stresses and hardening parameters
  matrix dgdsdq(ASTCKMAT(ncompstr,ncomphard));
  //  vector of derivatives of yield function with respect to hardening parameters
  vector dfdq(ASTCKVEC(ncomphard));
  //  matrix of derivatives of hardening function with respect to stresses
  matrix dhds (ASTCKMAT(ncomphard,ncompstr));
  //  matrix of derivatives of hardening function with respect to hardening parameters
  matrix dhdq(ASTCKMAT(ncomphard,ncomphard));
  //  hardening/softening parameters
  vector hs(ASTCKVEC(ncomphard));

  //  auxiliary matrices
  matrix auxm(ASTCKMAT(ncompstr,ncompstr)),auxm2(ASTCKMAT(ncompstr,ncomphard));
  //  auxiliary vector
  vector auxv(ASTCKVEC(ncompstr)),auxv2(ASTCKVEC(ncompstr));


  //  elastic stiffness matrix
  elmatstiff (d,ipp);
  //  elastic compliance matrix
  elmatcompl (c,ipp);
  
  //  initialization of increment of consistency parameter
  dgamma=0.0;
  
  //  elastic strain
  subv (epsn,epsp,epse);
  //  trial stress computation
  mxv (d,epse,sig);
  // add initial stresses (eigenstresses) if defined
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (i=0; i<ncompstr; i++)
      sig(i) += eigstresses[ipp][i];
  }
  
  //  initial stresses
  if (ic){
    for (j=0; j < ncompstr; j++)
      sig[j] += ic[ipp][3+j];
  }

  //  conversion of vector to tensor notation
  give_full_vector (sigt,sig,ssst);

  //  copy of plastic strains
  copyv (epsp,epspold);
  //  copy of hardening parameters 
  copyv (q,qold);
  
  //  evaluation of the yield function
  f = yieldfunction (ipp,im,sigt,q);
  
  fprintf (Out,"\n stress return vstup  f %le",f);
  
  if (f>err){
    //  stress is outside from yield surface
    //  stress return algorithm must be used
    
    //  main iteration loop
    for (i=0;i<ni;i++){
      
      //fprintf (Out,"\n stress return %ld",i);
      
      fillm (0.0,a);
      nullv (x);
      nullv (y);
      
      //  Jacobian matrix assembling
      //  submatrix A_11 = C + \gamma d2g/dsds
      
      //  matrix of second derivatives of plastic potential with respect to stresses
      dgdsigmadsigma (ipp, im, sigt, q, dgdsds);
      
      cmulm (dgamma,dgdsds);
      addm (c,dgdsds,dgdsds);
      
      //fprintf (Out,"\n d . dgdsds . dgamma  %le",auxm[0][0]);

      for (j=0;j<ncompstr;j++){
	for (k=0;k<ncompstr;k++){
	  a(j,k)=dgdsds(j,k);
	}
      }
      
      //  submatrix A_12 = dg/ds
      
      //  vector of first derivatives of plastic potential with respect to stresses
      dgdsigma (ipp, im, sigt, q, dgds);
      
      //fprintf (Out,"\n d . dgds   %le",auxv[0]);
      
      for (j=0;j<ncompstr;j++){
	a(j,ncompstr)=dgds(j);
      }
      
      //  submatrix A_21 = df/ds
      
      //  vector of first derivatives of yield function with respect to stresses
      dfdsigma (ipp, im, sigt, q, dfds);
      
      //fprintf (Out,"\n dfds    %le",dfds[0]);

      for (j=0;j<ncompstr;j++){
	a(ncompstr,j)=dfds(j);
      }
      
      //  submatrix A_22 = 0
      //  there is no contribution to the A_22 submatrix
      
      
      //  submatrices generated by hardening terms
      if (ncomphard>0){
	
	//  submatrix A_13 = gamma d2g/dsdq
	
	//  matrix of second derivatives of plastic potential with respect to stresses and hardening parameters
	dgdsigmadqpar (ipp,im,sigt,q,dgdsdq);
	
	for (j=0;j<ncompstr;j++){
	  for (k=0;k<ncomphard;k++){
	    a(j,ncompstr+1+k)=dgdsdq(j,k)*dgamma;
	  }
	}
	
	//fprintf (Out,"\n d . dgdsdq . dgamma    %le",auxm2[0][0]*dgamma);
	
	//  submatrix A_31 = gamma dh/ds

	//  matrix of derivatives of hardening function with respect to stresses
	dhdsigma (ipp,im,sigt,q,dhds);
	
	//fprintf (Out,"\n dhds   %le",dhds[0][0]);
	
	for (j=0;j<ncomphard;j++){
	  for (k=0;k<ncompstr;k++){
	    a(ncompstr+1+j,k)=dhds(j,k)*dgamma;
	  }
	}
	
	//  submatrix A_23 = df/dq
	
	//  vector of derivatives of yield function with respect to hardening parameters
	dfdqpar (ipp,im,sigt,q,dfdq);
	
	//fprintf (Out,"\n dfdq    %le",dfdq[0]);
	
	for (j=0;j<ncomphard;j++){
	  a(ncompstr,ncompstr+1+j)=dfdq(j);
	}
	
	//  submatrix A_32 = h
	
	//  vector of hardening functions
	hardvect (ipp,im,hs);
	
	//fprintf (Out,"\n dhdc   %le",dhdc[0]);

	for (j=0;j<ncomphard;j++){
	  a(ncompstr+1+j,ncompstr)=hs(j);
	}
	

	//  submatrix A_33 = I + gamma dh/dq
	
	//  matrix of derivatives of hardening function with respect to hardening parameters
	dhdqpar (ipp,im,ido,sigt,gamma,q,dhdq);
	
	//fprintf (Out,"\n dhdq   %le",dhdq[0][0]);

	for (j=0;j<ncomphard;j++){
	  for (k=0;k<ncomphard;k++){
	    a(ncompstr+1+j,ncompstr+1+k)=dhdq(j,k)*dgamma;
	  }
	  a(ncompstr+1+j,ncompstr+1+j)+=1.0;
	}
      }//  end of hardening block
      
      
      //  right hand side assembling
      //  subvector y_1
      
      for (j=0;j<ncompstr;j++){
	y[j]=epsp(j)-epspold(j)-dgamma*dgds(j);
      }
      
      //  subvector y_2
      y[ncompstr]=0.0-f;
      
      //  subvector y_3
      for (j=0;j<ncomphard;j++){
	y(ncompstr+1+j)=qold[j]-q[j]-dgamma*hs(j);
      }
      
      //fprintf (Out,"\n prava strana \n");
      //for (j=0;j<n;j++){
      //fprintf (Out,"  %le",y[j]);
      //}
      
      //fprintf (Out,"\n matice");
      //for (j=0;j<n;j++){
      //fprintf (Out,"\n");
      //for (k=0;k<n;k++){
      //fprintf (Out,"  %le",a[j*n+k]);
      //}
      //}

      //  solution of linearized system of equations
      gemp (a.a,x.a,y.a,n,1,1.0e-14,1);
      
      //fprintf (Out,"\n reseni \n");
      //for (j=0;j<n;j++){
      //fprintf (Out,"  %le",x[j]);
      //}
      
      //  norm of stresses
      norstr=0.0;
      //  norm of stress increments
      norstrincr=0.0;
      
      for (j=0;j<ncompstr;j++){
	norstrincr+=x[j]*x[j];
	norstr+=sig[j]*sig[j];
      }
      
      norconspar=x[ncompstr]*x[ncompstr];
      
      //  norm of hardening parameters
      norhardpar=0.0;
      //  norm of increments of hardening parameters
      norhardparincr=0.0;
      
      for (j=0;j<ncomphard;j++){
	norhardpar+=q[j]*q[j];
      }
      for (j=ncompstr+1;j<ncompstr+1+ncomphard;j++){
	norhardparincr+=x[j]*x[j];
      }
      
      
      //  stress update
      for (j=0;j<ncompstr;j++){
	sig[j]+=x[j];
	auxv[j]=x[j];
      }
      
      //  consistency parameter update
      dgamma+=x[ncompstr];
      
      //  hardening parameter update
      for (j=0;j<ncomphard;j++){
	q[j]+=x[ncompstr+1+j];
      }
      
	
      //fprintf (Out,"\n");
      //for (j=0;j<ncompstr;j++){
      //fprintf (Out,"   %le",sig[j]);
      //}
      //fprintf (Out,"       %le",dgamma);
      //for (j=0;j<ncomphard;j++){
      //fprintf (Out,"  %le",q[j]);
      //}
      
      //  plastic strains update
      //mxv (c,auxv,auxv2);
      mxv (c,sig,auxv2);
      
      for (j=0;j<ncompstr;j++){
	epsp[j]=epsn[j]-auxv2[j];
      }


      //  conversion of vector to tensor notation
      give_full_vector (sigt,sig,ssst);

      //updateq (ipp,im,ido,dgamma,epsn,epsp,sigt,q);
	
      //  evaluation of the yield function
      f = yieldfunction (ipp,im,sigt,q);
      
      fprintf (Out,"\n stress return i %4ld  f %le",i,f);

      //  evaluation of hardening/softening function
      //hardsoftfunction (ipp,im,sig,q,hs);


      //for (j=0;j<ncompstr;j++){
      //fprintf (Out,"   %le",epsp[j]);
      //}

      if (f<err)
        break;

      /*
      j=0;  k=0;
      if (norstr>zero){
	j++;
	if (norstrincr/norstr<err)  k++;
      }
      if (norhardpar>zero){
	j++;
	if (norhardparincr/norhardpar<err)  k++;
      }
      if (dgamma>zero){
	j++;
	if (norconspar/dgamma<err)  k++;
      }
      */
      /*
      if (j==k){
      fprintf (stdout,"\n bod %ld  yield function  %le   q %le    epsp %le %le %le   nor %le",
      ipp,f,q[0],epsp[0],epsp[1],epsp[2], epsp[0]*epsp[0]+epsp[1]*epsp[1]+epsp[2]*epsp[2]);
	break;
      }
      */
    }
  }
  
  //  update of consistency parameter
  gamma+=dgamma;
  
  //  storing resulting stresses
  storestress (0,ipp,sig);
  
}



/**
  Function returns stresses on sufrace of plasticity.
  Newton method is used.

  gamma,epsp and q will be replaced by new values

  @param[in]     ipp   - integration point id
  @param[in]     im    - id of material type (for combination of material models, default value is 0)
  @param[in]     ido   - first index in eqother array (for combination of material models, default value is 0)
  @param[in/out] gamma - consistency parameter (input/output)
  @param[in]     epsn  - total strain components
  @param[in/out] epsp  - plastic strain components (input/output)
  @param[in/out] q     - hardening parameters (input/output)
  @param[in]     ni    - maximum number of iterations
  @param[in]     err   - required error
   
  @return The function returns actual plastic strain %vector in the parameter epsp,
          actual consistency parameter in the parameter gamma and actual %vector of
          hardening parameters in the parameter q.

  Created by JK, 16.2.2007
*/
void mechmat::newton_stress_return_2 (long ipp,long im,long /*ido*/,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err)
{
  long i,j,k,n,ncompstr,ncomphard;
  double zero,f,dgamma,norstr,norstrincr,norconspar,norhardpar,norhardparincr;


  //  computer zero
  zero=Mp->zero;

  //  number of strain components
  ncompstr = epsn.n;
  //  number of hardening components
  ncomphard = q.n;
  //  number of unknowns in the Newton method
  n = ncompstr + 1 + ncomphard;
  //  unknowns in the Newton-Raphson method are organized as follows:
  //  increments of stresses (ncompstr components)
  //  increment of consistency parameter (1 component)
  //  increments of hardening/softening parameters (ncomphard components)
  //
  //  hardening parameters may be omitted
  vector qold(ASTCKVEC(ncomphard));

  //  stiffness matrix of the material
  matrix d(ASTCKMAT(ncompstr,ncompstr));
  //  compliance matrix of the material
  matrix c(ASTCKMAT(ncompstr,ncompstr));
  //  vector of elastic strains
  vector epse(ASTCKVEC(ncompstr));
  //  vector of stresses
  vector sig(ASTCKVEC(ncompstr));
  //  vector of stresses
  vector sigt(ASTCKVEC(6));
  //  backup of the trial stresses
  vector sigtrial(ASTCKVEC(ncompstr));
  //  backup of attained hardening parameters
  //vector qold(ncomphard);
  //  matrix of second derivatives of plastic potential with respect to stresses
  matrix dgdsds(ASTCKMAT(ncompstr,ncompstr));
  //  vector of derivatives of plastic potential with respect to stresses
  vector dgds(ASTCKVEC(ncompstr));
  //  vector of derivatives of yield function with respect to stresses
  vector dfds(ASTCKVEC(ncompstr));
  //  matrix of second derivatives of plastic potential with respect to stresses and hardening parameters
  matrix dgdsdq(ASTCKMAT(ncompstr,ncomphard));
  //  vector of derivatives of yield function with respect to hardening parameters
  vector dfdq(ASTCKVEC(ncomphard));
  //  matrix of derivatives of hardening function with respect to stresses
  matrix dhds (ASTCKMAT(ncompstr,ncomphard));
  //  vector of derivatives of hardening function with respect to consistency parameter
  vector dhdc (ASTCKVEC(ncomphard));
  //  matrix of derivatives of hardening function with respect to hardening parameters
  matrix dhdq(ASTCKMAT(ncomphard,ncomphard));
  //  hardening/softening parameters
  vector hs(ASTCKVEC(ncomphard));

  //  auxiliary matrices
  matrix auxm(ASTCKMAT(ncompstr,ncompstr)), auxm2(ASTCKMAT(ncompstr,ncomphard));
  //  auxiliary vector
  vector auxv(ncompstr),auxv2(ncompstr);


  //  allocation of memory for Newton method
  //  array containing matrix of the system (Jacobian matrix)
  vector a(ASTCKVEC(n*n));
  //  array containing solution of the system
  vector x(ASTCKVEC(n));
  //  array containing right hand side
  vector y(ASTCKVEC(n));
  
  //  elastic stiffness matrix
  elmatstiff (d,ipp);
  //  elastic compliance matrix
  elmatcompl (c,ipp);
  
  //  initialization of increment of consistency parameter
  dgamma=0.0;
  
  //  elastic strain
  subv (epsn,epsp,epse);
  //  trial stress computation
  mxv (d,epse,sig);
  
  // add intial stresses (eigenstresses) if defined
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (i=0; i<ncompstr; i++)
      sig(i) += eigstresses[ipp][i];
  }

  //  conversion of vector to tensor notation
  //vector_tensor (sig,sigt,ip[ipp].ssst,stress);
  //  copy of hardening parameters 
  copyv (q,qold);
  
  //  evaluation of the yield function
  f = yieldfunction (ipp,im,sig,q);
  
  //  evaluation of hardening/softening function
  hardsoftfunction (ipp,im,sig,q,hs);

  if (f>err){
    //  stress is outside from yield surface
    //  stress return algorithm must be used
    
    //  backup of the trial stresses (will be used on the right hand side)
    copyv (sig,sigtrial);
    
    //  main iteration loop
    for (i=0;i<ni;i++){
      
      //fprintf (Out,"\n %ld",i);

      nullv(a);
      nullv(x);
      nullv(y);
      
      //  Jacobian matrix assembling
      //  submatrix A_11
      
      //  matrix of second derivatives of plastic potential with respect to stresses
      dgdsigmadsigma (ipp, im, sig, q, dgdsds);
      
      mxm (d,dgdsds,auxm);
      cmulm (dgamma,auxm);
      
      //fprintf (Out,"\n d . dgdsds . dgamma  %le",auxm[0][0]);

      for (j=0;j<ncompstr;j++){
	for (k=0;k<ncompstr;k++){
	  a[j*n+k]=auxm[j][k];
	}
	a[j*n+j]+=1.0;
      }
      
      //  submatrix A_12
      
      //  vector of first derivatives of plastic potential with respect to stresses
      dgdsigma (ipp, im, sig, q, dgds);
      
      mxv (d,dgds,auxv);
      
      //fprintf (Out,"\n d . dgds   %le",auxv[0]);
      
      for (j=0;j<ncompstr;j++){
	a[j*n+ncompstr]=auxv[j];
      }
      
      //  submatrix A_21
      
      //  vector of first derivatives of yield function with respect to stresses
      dfdsigma (ipp, im, sig, q, dfds);
      
      //fprintf (Out,"\n dfds    %le",dfds[0]);

      for (j=0;j<ncompstr;j++){
	a[ncompstr*n+j]=dfds[j];
      }
      
      //  submatrix A_22
      //  there is no contribution to the A_22 submatrix
      
      
      //  submatrices generated by hardening terms
      if (ncomphard>0){
	
	//  submatrix A_13
	
	//  matrix of second derivatives of plastic potential with respect to stresses and hardening parameters
	dgdsigmadqpar (ipp,im,sig,q,dgdsdq);
	
	mxm (d,dgdsdq,auxm2);
	
	for (j=0;j<ncompstr;j++){
	  for (k=0;k<ncomphard;k++){
	    a[j*n+ncompstr+1+k]=auxm2[j][k]*dgamma;
	  }
	}
	
	//fprintf (Out,"\n d . dgdsdq . dgamma    %le",auxm2[0][0]*dgamma);
	
	//  submatrix A_31
	
	//  matrix of derivatives of hardening function with respect to stresses
	dhdsigma (ipp,im,sig,q,dhds);
	
	
	//fprintf (Out,"\n dhds   %le",dhds[0][0]);
	
	for (j=0;j<ncomphard;j++){
	  for (k=0;k<ncompstr;k++){
	    a[(ncompstr+1+j)*n+k]=dhds[k][j]*dgamma;
	  }
	}
	
	//  submatrix A_23
	
	//  vector of derivatives of yield function with respect to hardening parameters
	dfdqpar (ipp,im,sig,q,dfdq);
	
	//fprintf (Out,"\n dfdq    %le",dfdq[0]);
	
	for (j=0;j<ncomphard;j++){
	  a[ncompstr*n+ncompstr+1+j]=dfdq[j];
	}
	
	//  submatrix A_32
	
	//  vector of derivatives of hardening function with respect to consistency parameter
	//dhdgamma (ipp,im,epsp,sig,q,dhdc);
	
	//fprintf (Out,"\n dhdc   %le",dhdc[0]);
	
	for (j=0;j<ncomphard;j++){
	  a[(ncompstr+1+j)*n+ncompstr]=dhdc[j];
	}
	
	//  submatrix A_33
	
	//  matrix of derivatives of hardening function with respect to hardening parameters
	//dhdqpar (ipp,im,ido,sigt,dgamma,q,dhdq);
	
	//fprintf (Out,"\n dhdq   %le",dhdq[0][0]);
	
	for (j=0;j<ncomphard;j++){
	  for (k=0;k<ncomphard;k++){
	    //a[(ncompstr+1+j)*n+ncompstr+1+k]=dhdq[j][k];
	  }
	  a[(ncompstr+1+j)*n+ncompstr+1+j]-=1.0;
	}

      }
      
      //  right hand side assembling
      //  subvector y_1
      
      for (j=0;j<ncompstr;j++){
	y[j]=sigtrial[j]-sig[j]-dgamma*auxv[j];
      }
      
      //  subvector y_2
      y[ncompstr]=0.0-f;
      
      //  subvector y_3
      for (j=0;j<ncomphard;j++){
	//y[ncompstr+1+j]=q[j]-hs[j];
      }
      
      //fprintf (Out,"\n prava strana \n");
      //for (j=0;j<n;j++){
      //fprintf (Out,"  %le",y[j]);
      //}
      
      //fprintf (Out,"\n matice");
      //for (j=0;j<n;j++){
      //fprintf (Out,"\n");
      //for (k=0;k<n;k++){
      //fprintf (Out,"  %le",a[j*n+k]);
      //}
      //}

      //  solution of linearized system of equations
      gemp (a.a,x.a,y.a,n,1,1.0e-10,1);
      
      //fprintf (Out,"\n reseni \n");
      //for (j=0;j<n;j++){
      //fprintf (Out,"  %le",x[j]);
      //}
      
      //  norm of stresses
      norstr=0.0;
      //  norm of stress increments
      norstrincr=0.0;
      
      for (j=0;j<ncompstr;j++){
	norstrincr+=x[j]*x[j];
	norstr+=sig[j]*sig[j];
      }
      
      norconspar=x[ncompstr]*x[ncompstr];
      
      //  norm of hardening parameters
      norhardpar=0.0;
      //  norm of increments of hardening parameters
      norhardparincr=0.0;
      
      for (j=0;j<ncomphard;j++){
	norhardpar+=q[j]*q[j];
      }
      for (j=ncompstr+1;j<ncompstr+1+ncomphard;j++){
	norhardparincr+=x[j]*x[j];
      }
      
      
      //  stress update
      for (j=0;j<ncompstr;j++){
	sig[j]+=x[j];
	auxv[j]=x[j];
      }
      
      //  consistency parameter update
      dgamma+=x[ncompstr];
      
      //  hardening parameter update
      for (j=0;j<ncomphard;j++){
	q[j]+=x[ncompstr+1+j];
      }
      
	
      //fprintf (Out,"\n");
      //for (j=0;j<ncompstr;j++){
      //fprintf (Out,"   %le",sig[j]);
      //}
      //fprintf (Out,"       %le",dgamma);
      //for (j=0;j<ncomphard;j++){
      //fprintf (Out,"  %le",q[j]);
      //}
      
      //  plastic strains update
      mxv (c,auxv,auxv2);
      
      for (j=0;j<ncompstr;j++){
	epsp[j]-=auxv2[j];
      }


      //  conversion of vector to tensor notation
      //vector_tensor (sig,sigt,ip[ipp].ssst,stress);

      //  evaluation of the yield function
      f = yieldfunction (ipp,im,sig,q);
      
      //  evaluation of hardening/softening function
      hardsoftfunction (ipp,im,sig,q,hs);


      //for (j=0;j<ncompstr;j++){
      //fprintf (Out,"   %le",epsp[j]);
      //}



      
      j=0;  k=0;
      if (norstr>zero){
	j++;
	if (norstrincr/norstr<err)  k++;
      }
      if (norhardpar>zero){
	j++;
	if (norhardparincr/norhardpar<err)  k++;
      }
      if (dgamma>zero){
	j++;
	if (norconspar/dgamma<err)  k++;
      }
      
      
      if (j==k){
	/*
	  fprintf (stdout,"\n bod %ld  yield function  %le   q %le    epsp %le %le %le   nor %le",
	  ipp,f,q[0],epsp[0],epsp[1],epsp[2], epsp[0]*epsp[0]+epsp[1]*epsp[1]+epsp[2]*epsp[2]);
	*/
	break;
      }
    }
  }
  
  fprintf (stdout,"\n\n parametr zpevneni  %lf",q[0]);
  fprintf (stdout,"\n parametr konzistence  %lf\n",gamma);
  
  //  update of consistency parameter
  gamma+=dgamma;
  
  //  storing resulting stresses
  storestress (0,ipp,sig);
}



/**
  Function detects evolved plasticity according to actual values of 
  consistency parameter.It is used for adaptivity computations.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by LS
*/
void mechmat::refresh_plast (long im, long ido)
{
  if (plast == 0)
    for (long i=0;i<tnip;i++)
      if ( give_consparam(i,im,ido) ){
	plast = 1;
	return;
      }
}






// *****************************************************************
// *****************************************************************
//   PART CONTAINING FUNCTIONS DEALING WITH DAMAGE
// *****************************************************************
// *****************************************************************


/**
  Function returns damage function parameters for actual values of total strains
  for scalar damage models.

  @param[in]  ipp   - integration point pointer
  @param[in]  eps   - total strain components
  @param[in]  im    - index of material type for given ip
  @param[out] kappa - damage function parameters

  @return The function returns actual value of damage function parameters in the parameter kappa.

  Created by Tomas Koudelka, 27.1.2003
*/
void mechmat::damfuncpar(long ipp, long im, vector &eps, vector &kappa)
{
  switch (ip[ipp].tm[im])
  {
    case scaldamage:
    {
      scdam[ip[ipp].idm[im]].damfuncpar(ipp, eps, kappa);
      break;
    }
    case scaldamagecc:
    {
      kappa(0) = scdamcc[ip[ipp].idm[im]].damfuncpar(ipp, eps);
      break;
    }
    default:
      print_err("unknown scalar damage model is required",__FILE__,__LINE__,__func__);
  }
  return;
}



/**
  Function returns damage function for actual values of kappa
  for scalar damage models.

  @param[in] ipp    - integration point pointer
  @param[in] im     - material index
  @param[in] ido    - index of internal variables for given material in the ipp other array
  @param[in] kappa  - %vector of damage function parameters
  @param[in] eps    - strain %vector
  @param[in] omegao - array of attained values of damage parameters

  @return The function returns actual value of damage parameter.
  
  Created by Tomas Koudelka, 27.1.2003
  Modified by Tomas Koudelka, 4.2.2010
*/
double mechmat::damfunction(long ipp, long im, long /*ido*/,vector &kappa, vector &/*eps*/, vector &omegao)
{
  double ret = 0.0;

  switch (ip[ipp].tm[im])
  {
    case scaldamage:
    {
      ret = scdam[ip[ipp].idm[im]].damfunction(ipp, kappa, omegao);
      break;
    }
    default:
      print_err("unknown scalar damage model is required",__FILE__,__LINE__,__func__);
  }
  return ret;
}



/**
  Function returns limit elastic strain
  for scalar damage models.

  @param[in] ipp - integration point pointer
  @param[in] im  - index of damage material in the ip.tm array

  @return The function returns treshold value of elastic strain.

  Created by Frederic DUFOUR, 23.2.2003
*/
double mechmat::epsefunction(long ipp, long im)
{
  double ret = 0.0;

  switch (ip[ipp].tm[im]){
    case scaldamage:
    {
      ret = scdam[ip[ipp].idm[im]].epsefunction(ipp);
      break;
    }
    case scaldamagecc:
    {
      ret = scdamcc[ip[ipp].idm[im]].epsefunction(ipp);
      break;
    }
    case damplifmat:
    {
      ret = damplifm[ip[ipp].idm[im]].epsefunction(ipp);
      break;
    }
    default:
      print_err("unknown scalar damage model in function",__FILE__,__LINE__,__func__);
  }
  return ret;
}



/**
  Function returns stresses for scalar damage models and value of damage function.
  Original formulation for the nonlocal model used in SIFEL benchmark.

  @param[in]     ipp    - integration point pointer
  @param[in]     im     - material index
  @param[in]     ido    - index of internal variables for given material in the ipp other array
  @param[in]     eps    - total strain components
  @param[in/out] kappa  - damage function parameters from previous stesp (input/output)
  @param[in]     omegao - damage parameters from the previous step
  @param[out]    sigma  - resulting stresses
  @param[in]     d      - elastic stiffness matrix with respect to timedependent Young modulus

  @return The function returns actual value of damage parameter \omega and %vector of actual 
          stresses in the parameter sigma.

  Created by Tomas Koudelka, 27.1.2003
*/
double mechmat::scal_dam_sol (long ipp,long im,long ido,vector &eps,vector &kappa, vector &omegao, vector &sigma, matrix &d)
{
  long nk = kappa.n;
//  long ncomp=ip[ipp].ncompstr;
  double epse, omega;
  vector updkappa(nk);
//  matrix d(ncomp, ncomp);


  if (nk > 1)
  {
    print_err("unsupported number of componetnts of vector kappa", __FILE__, __LINE__,__func__);
    return 0.0;
  }

  if ((Mp->matmodel == local) || ((ip[ipp].hmt & 2) == 0)) // model is not nonlocal
  // local model is used
  {
    damfuncpar(ipp, im, eps, updkappa);
    if (updkappa[0] > kappa[0])
    {
      kappa[0] = updkappa[0];
      epse = epsefunction(ipp,im) ;
      if (kappa[0] > epse)
        omega = damfunction(ipp, im, ido, kappa, eps, omegao);
      else
        omega = 0.0;
      if (omega < omegao[0])
        omega = omegao[0];
    }
    else
      omega = omegao[0];

    mxv(d, eps, sigma);
    cmulv(1.0-omega, sigma);
    return omega;
  }

  if (Mp->matmodel == nonlocal)
  // nonlocal model is used
  {
    if (Mp->nonlocphase == 1)
    // compute local values for averaging -> nonlocphase=1 -> no stresses are required
    {
      damfuncpar(ipp, im, eps, updkappa);
      kappa[0] = updkappa[0]; // the updating of the kappa is required
      nullv(sigma);
      return 0.0;
    }
    if (Mp->nonlocphase == 2)
    // compute stresses from the averaged nonlocal values -> nonlocphase=2
    {
      updkappa[0] = ip[ipp].nonloc[0];
      if (updkappa[0] > kappa[0])
      {
        kappa[0] = updkappa[0];
        epse = epsefunction(ipp, im);
        if (kappa[0] > epse)
          omega = damfunction(ipp, im, ido, kappa, eps, omegao);
        else
          omega = 0.0;
        if (omega < omegao[0])
          omega = omegao[0];
      }
      else
        omega = omegao[0];

      mxv(d, eps, sigma);
      cmulv(1.0-omega, sigma);
      return omega;
    }
  }

  // Following command is only for formal purposes (compiler warnings)
  return 0.0;
}



/**
  Function returns stresses for scalar damage models and value of damage function.
  New version definition of nonlocal model where equivalent strain is calculated 
  from averaged strain tensor. Originally, averageing of equivalent strain directly 
  was defined.

  @param[in]     ipp    - integration point pointer
  @param[in]     im     - material index
  @param[in]     ido    - index of internal variables for given material in the ipp other array
  @param[in]     eps    - total strain components
  @param[in/out] kappa  - damage function parameters from previous stesp (input/output)
  @param[in]     omegao - damage parameters from the previous step
  @param[out]    sigma  - resulting stresses
  @param[in]     d      - elastic stiffness matrix with respect to timedependent Young modulus

  @return The function returns actual value of damage parameter \omega and %vector of actual 
          stresses in the parameter sigma.

  Created by Tomas Koudelka, 16.2.2022
*/
double mechmat::scal_dam_sol2 (long ipp,long im,long ido,vector &eps,vector &kappa, vector &omegao, vector &sigma, matrix &d)
{
  long nk = kappa.n;
//  long ncomp=ip[ipp].ncompstr;
  double epse, omega;
  vector updkappa(nk);
//  matrix d(ncomp, ncomp);


  if (nk > 1)
  {
    print_err("unsupported number of componetnts of vector kappa", __FILE__, __LINE__,__func__);
    return 0.0;
  }

  if ((Mp->matmodel == local) || ((ip[ipp].hmt & 2) == 0)) // model is not nonlocal
  // local model is used
  {
    damfuncpar(ipp, im, eps, updkappa);
    if (updkappa[0] > kappa[0])
    {
      kappa[0] = updkappa[0];
      epse = epsefunction(ipp,im) ;
      if (kappa[0] > epse)
        omega = damfunction(ipp, im, ido, kappa, eps, omegao);
      else
        omega = 0.0;
      if (omega < omegao[0])
        omega = omegao[0];
    }
    else
      omega = omegao[0];

    mxv(d, eps, sigma);
    cmulv(1.0-omega, sigma);
    return omega;
  }

  if (Mp->matmodel == nonlocal)
  {
    if (Mp->nonlocphase == 1){
      // compute local values for averaging in nonlocphase=1 -> no stresses are required
      // in this case, nothing is needed to be computed in the first phase,
      // just strain tensor will be averaged
      return 0.0;
    }
    if (Mp->nonlocphase == 2){
      // nonlocphase = 2 -> compute stresses from the averaged nonlocal values
      // collect averaged strains
      long ncomp = eps.n;
      vector epsn(ASTCKVEC(ncomp));
      for (long i=0; i<ncomp; i++)
        epsn(i) = Mm->ip[ipp].nonloc[i];
      // update equivalent strain
      damfuncpar(ipp, im, epsn, updkappa);
      if (updkappa[0] > kappa[0]){
        // new maximum value of equivalent strain kappa was attained
        kappa[0] = updkappa[0];
        epse = epsefunction(ipp, im);
        if (kappa[0] > epse)
          omega = damfunction(ipp, im, ido, kappa, eps, omegao);
        else
          omega = 0.0;
        if (omega < omegao[0])
          omega = omegao[0];
      }
      else
        omega = omegao[0];

      mxv(d, eps, sigma);
      cmulv(1.0-omega, sigma);
      return omega;
    }
  }

  // Following command is only for formal purposes (compiler warnings)
  return 0.0;
}






// *****************************************************************
// *****************************************************************
//   PART CONTAINING FUNCTIONS DEALING WITH VISCO-PLASTICITY
// *****************************************************************
// *****************************************************************



/**
  Function restores part of stress increments to %vector sig.

  @param[in]  lcid - load case id
  @param[in]  ipp  - integration point number
  @param[in]  im   - index of the material in the tm and idm arrays on integration point
  @param[in]  ido  - index of internal variables in the ip's ipp eqother array
  @param[in]  fi   - first index of the required stress increment component
  @param[out] sig  - %vector containing stress increment components
   
  @return The function returns %vector of stress increments in the parameter sig.

  Created by JK, TKo, 29.4.2008
*/
void mechmat::givestressincr (long lcid,long ipp,long im,long ido,long fi,vector &sig)
{
  long i,ncompo;

  if (im == 0)
    nullv(sig);

  switch (ip[ipp].tm[im]){
  case effective_stress:{
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    givestressincr (lcid, ipp, im+1, ido+ncompo, fi, sig);
    break;
  }
  case shrinkagemat:{
    i=ip[ipp].idm[im];
    shmat[i].givestressincr (lcid, ipp, im, ido, fi, sig);
    break;
  }
  case creep_damage:{
    ncompo  = givencompeqother(ipp, im);
    ncompo -= givencompeqother(ipp, im+1);
    ncompo -= givencompeqother(ipp, im+2);
    givestressincr (lcid, ipp, im+1, ido+ncompo, fi, sig);
    break;
  }
  case time_switchmat:{
    tswmat[ip[ipp].idm[im]].givestressincr (lcid, ipp, im, ido, fi, sig);
    break;
  }
  case viscoplasticity:{
    i=ip[ipp].idm[im];
    visplas[i].givestressincr (ipp,ido,fi,sig);
    break;
  }
  case creepb3:
  case creeprs:
  case creepdpl:{
    creep_givestressincr (ipp,im,ido,fi,sig);
    break;
  }
  case hypoplastmat:{
    i=ip[ipp].idm[im];
    hypopl[i].givestressincr(ipp, im, ido, sig);
    break;
  }
  case hypoplastusatthermat:{
    i=ip[ipp].idm[im];
    hypoplustherm[i].givestressincr(lcid, ipp, im, ido, fi, sig);
    break;
  }
  case bbmcoupmat:{
    i=ip[ipp].idm[im];
    bbm[i].givestressincr(lcid, ipp, im, ido, fi, sig);
  }
  case doublestructuremat:{
    i=ip[ipp].idm[im];
    dsm[i].givestressincr(lcid, ipp, im, ido, fi, sig);
  }
  /*
  case consolidation:{
    csol[ip[ipp].idm[im]].give_dstresses_eqother (ipp,ido,sig);
    break;
  }
  case lemaitr:{
    lmtr[ip[ipp].idm[im]].givestressincr (ipp,ido,sig);
    break;
  }
  case creepbaz:{
    crbaz[ip[ipp].idm[im]].givestressincr (ipp,ido,sig);
    break;
  }
  */
  case elisomat:
  case elortomat:
  case homomatm:
  case simplas1d:
  case simvisc:
  case jflow:
  case microplaneM4:
  case microsimp:
  case microfibro:
  case mohrcoul:
  case mohrcoulparab:
  case boermaterial:
  case druckerprager:
  case doubledrprager:
  case druckerprager2:
  case chenplast:
  case modcamclaymat:
  case modcamclaycoupmat:
  case shefpl:
  case glasgowmechmat:
  case glasgowdamage:
  case scaldamage:
  case fixortodamage:
  case scaldamagecc:
  case ortodamage:
  case ortodamage2:
  case ortodamagerot:
  case anisodamage:
  case anisodamagerot:
  case graphm:
  case varelisomat:
  case elisopdmat:
  case geoelast:
  case creepeffym: 
  case winklerpasternak:
  case nonlocdamgmat:
  case nonlocplastmat:
  case nonlocalmicroM4:
  case damage_plasticity:
  case cebfipcontmat:
  case contmat:
  case damplifmat:
  case plastifmat:
  case elasttimemat:
  case layerplate:
    break;    
  default:{
    print_err("unknown material model is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Function computes D \delta \epsilon
  it is used for right hand side computation in time dependent problems.

  @param[in]  ipp - integration point pointer
  @param[in]  im  - material index
  @param[in]  ido - index of internal variables for given material in the ipp other array
  @param[out] sig - stress increment components
  @param[in]  q   - internal variables
  @param[in]  dt  - time increment

  @return The function returns stress increment %vector in the parameter sig.
   
  Created by JK, 28.10.2001
*/
void mechmat::stiff_deps_vispl (long ipp,long im,long ido,vector &sig,vector &q,double dt)
{
  long i,ncomp,ncomph,ncompo;
  //  number of strain/stress components
  ncomp=sig.n;
  //  number of internal parameters
  ncomph=q.n;
  
  //  number of components of array eqother belonging to viscous material
  ncompo = givencompother(ipp,im+1);
  
  double f,g=0.0,dlambda,cs=0.0;
  vector epsp(ASTCKVEC(ncomp)),dq(ASTCKVEC(ncomph)),dh(ASTCKVEC(ncomph)),dfds(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp)),epspt(ASTCKMAT(3,3)), h(ASTCKMAT(ncomph, ncomph));
  vector sigt(ASTCKVEC(6));

  //  value of yield function
  f = yieldfunction (ipp,im+2,sig,q);
  //  derivative of yield function with respect to stress components
  dfdsigma (ipp,im+2,sig,q,dfds);
  //  derivative of yield function with respect to hardening parameters
  dfdqpar (ipp,im+2,sig,q,dq);
  
  switch (ip[ipp].tm[im+1]){
  case simvisc:{
    //  viscous function value
    g = svis[ip[ipp].idm[im+1]].gfun (f);
    break;
  }
  case lemaitr:{
    //  cumulative deformation
    cs=ip[ipp].other[ido+3*ncomp];
    //  viscous function value
    g = lmtr[ip[ipp].idm[im+1]].gfun (f,cs);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  dlambda = 0.0;
  if (g>0.0){
    //  lambda increment
    dlambda = dt*g;
    //  irreversible strain increment
    for (i=0;i<ncomp;i++){
      epsp[i]=dlambda*dfds[i];
    }

    //  hardening parameter increment
    plasticmoduli (ipp, im+2, epsp, sigt, h);
    mxv (h,dq,dh);
    for (i=0;i<ncomph;i++){
      dh[i]=-dlambda*dh[i];
    }
  }
  
  //  storage of new values
  switch (ip[ipp].tm[im+1]){
  case simvisc:{
    //  irreversible strain increments are added to irreversible strains
    for (i=0;i<ncomp;i++){
      ip[ipp].eqother[ido+ncompo+i]+=epsp[i];
    }
    //  consistency parameter
    ip[ipp].eqother[ido+ncompo+ncomp]+=dlambda;
    //  internal variable increment
    for (i=0;i<ncomph;i++){
      ip[ipp].eqother[ido+ncompo+ncomp+1+i]+=dh[i];
    }
    
    break;
  }
  case lemaitr:{
    //  irreversible strain increment
    for (i=0;i<ncomp;i++){
      ip[ipp].eqother[ido+ncompo+i]+=epsp[i];
    }
    //  cumulative irreversible strain computation
    vector_tensor (epsp,epspt,ip[ipp].ssst,strain);
    //cumulstrain (epspt,cs);
    ip[ipp].eqother[ido+ncomp*2]=cs;

    //  consistency parameter
    ip[ipp].eqother[ido+ncompo+ncomp]+=dlambda;
    
    //  internal variable increment
    for (i=0;i<ncomph;i++){
      ip[ipp].eqother[ido+ncompo+ncomp+1+i]+=dh[i];
    }
    break;
  }
  default:{
    print_err("unknown material type is requiered",__FILE__,__LINE__,__func__);
  }
  }

  //  stress increments
  matstiff (d,ipp);
  mxv (d,epsp,sig);
  
}



/**
  Function computes contribution to the right hand side
  in creep and shrinkage analysis.

  !!!??? unused function candidate for removal
  !!!??? modifies strain in the integration point ?!
   
  @param[in] ipp - integration point pointer

  Created by JK, 9.11.2001
*/
void mechmat::stiff_eps (long ipp)
{
  long i,ncomp=ip[ipp].ncompstr;
  double time,
    complian=0.0,
    strength,
    shrink;
  vector eps(ncomp),sig(ncomp);
  
  time=Mp->time/3600.0/24.0;
  
  switch (ip[ipp].tm[0]){
  case aci:{
    aci78mod[ip[ipp].idm[0]].compliance (time,complian,strength,shrink);
    break;
  }
  case cebfip:{
    cebfip78mod[ip[ipp].idm[0]].compliance (time,complian,strength,shrink);
    break;
  }
  default:
    print_err("unknown material type is requiered",__FILE__,__LINE__,__func__);
  }
  
  
  //  strain components
  for (i=0;i<ncomp;i++){
    ip[ipp].strain[i]=ip[ipp].stress[i]*complian;
  }
}



// *************************************************
// *************************************************
//
// Functions connected with Len-Jons potential model
//
// *************************************************
// *************************************************



/**
  The function returns actual value of of the first derivative of lenjones potential
  with respect to ???

  @param[in] ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of the first derivative.

  Created by JK,
*/
double mechmat::give_first_derivative (long ipp,double r)
{
  double d=0.0;
  
  switch (ip[ipp].tm[0]){
  case lenjonespot:{
    d = lenjon[ip[ipp].idm[0]].first_derivative (r);
    break;
  }
  default:{
    print_err("unknown material model is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return d;
}



/**
  The function returns actual value of the first derivative of lenjones potential
  with respect to ???

  @param[in] ipp - integration point number in the mechmat ip array.

  @return The function returns actual value of the second derivative.

  Created by JK,
*/
double mechmat::give_second_derivative (long ipp,double r)
{
  double d=0.0;
  
  switch (ip[ipp].tm[0]){
  case lenjonespot:{
    d = lenjon[ip[ipp].idm[0]].second_derivative (r);
    break;
  }
  default:{
    print_err("unknown material model is required",__FILE__,__LINE__,__func__);
  }
  }

  return d;
}




// *************************************************
// *************************************************
//
// Funtions for backup of regular int. point content
//
// *************************************************
// *************************************************



/**
  Function saves data from regular integration points into backup file
  in text format.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled

  @return The function does not return anything.
   
  Created by JK, 19.9.2004
  Rewritten by TKo, 9.2008
  Modified by TKo, 7.10.2013
*/
void mechmat::save_intpoints_txt (FILE *aux, sel &selelems, sel *selother)
{
  long i,j,k,n,nip,ir,ipp;
  sel selno; // selection with no selected components
  int prec = (int)Mp->hdbcont.prec;
  
  selno.n=1;
  fprintf(aux, "\n");

  if (nonmechq) // write table of indices of non-mechanical quantities
  {
    for (i=0; i<tnknmq; i++)
      fprintf(aux, "%ld ", nmqid[i]);
    fprintf(aux, "\n");
  }

  for (i=0;i<Mt->ne;i++)
  {    
    ipp = Mt->elements[i].ipp[0][0];
    nip = Mt->give_totnip(i);
    selelems.presence_id(i, ir);
    if (ir < 0)
      n = 0;    
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fprintf (aux,"%ld %ld %ld %ld\n",ipp, Mb->nlc, ip[ipp].ncompstr, n);
      if (ir < 0)
        ip[ipp].save_data_txt(aux,Mb->nlc, selno);
      else
        ip[ipp].save_data_txt(aux,Mb->nlc, selother[ir]);

      // writing TRFEL quantities
      if (nonmechq)
      {
        for (k=0; k<nnmq; k++)
          fprintf (aux,"%.*le ",prec, nonmechq[k][ipp]);
        fprintf(aux, "\n");
      }
      ipp++;
    }
  }
}



/**
  Function saves data from regular integration points into several backup files
  in text format.
   
  @param[in] ni       - time step id
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by JK, 19.9.2004
  Rewritten by TKo, 9.2008
  Modified by TKo, 7.10.2013
*/
void mechmat::save_intpoints_txt (long ni, sel &selelems, sel *selother)
{
  long i,j,k,ir,n,nip,ipp;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Mp->hdbcont.prec;

  sprintf(name, "%s.%ld.strain.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<tnip; i++)
  {
    fprintf (aux, "%ld %ld %ld\n", i, Mb->nlc, ip[i].ncompstr);
    for (j=0; j<Mb->nlc*ip[i].ncompstr; j++)
      fprintf (aux, "%.*le\n", prec, ip[i].strain[j]);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.stress.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<tnip; i++)
  {
    fprintf (aux, "%ld %ld %ld\n", i, Mb->nlc, ip[i].ncompstr);
    for (j=0; j<Mb->nlc*ip[i].ncompstr; j++)
      fprintf (aux, "%.*le\n", prec, ip[i].stress[j]);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.other.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0; i<Mt->ne; i++)
  {    
    ipp = Mt->elements[i].ipp[0][0];
    nip = Mt->give_totnip(i);
    selelems.presence_id(i, ir);    
    if (ir < 0)
      n = 0;
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fprintf (aux, "%ld %ld\n", ipp, n);
      if (n==0) 
      {
        ipp++;
        continue;
      }
      for (k=0; k<ip[ipp].ncompeqother; k++)
      {
        if (selother[ir].presence_id(k))
          fprintf (aux, "%.*le\n", prec, ip[ipp].eqother[k]);
      }
      ipp++;
    }
  }
  fclose(aux);
 
  // writing TRFEL quantities
  if (nonmechq)
  {
    sprintf(name, "%s.%ld.trfquant.bac", Mp->hdbcont.hdbnames, ni);
    aux = fopen(name, "wt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    // write table of indices of non-mechanical quantities
    for (i=0; i<tnknmq; i++)
      fprintf(aux, "%ld ", nmqid[i]);
    fprintf(aux, "\n");

    for (i=0; i<tnip; i++)
    {
      fprintf (aux, "%ld ", i);
      for (j=0; j<nnmq; j++)
        fprintf (aux,"%.*le ", prec, nonmechq[j][i]);
        
      fprintf(aux, "\n");
    }
    fclose(aux);
  }
}



/**
  Function saves data from regular integration points into backup file
  in binary format.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void mechmat::save_intpoints_bin (FILE *aux, sel &selelems, sel *selother)
{
  long i,j,k,n,nip,ir,ipp;
  sel selno; // selection with no selected components
  
  selno.n=1;

  if (nonmechq)  // write table of indices of non-mechanical quantities
    fwrite (nmqid, sizeof(*nmqid), tnknmq, aux);

  for (i=0;i<Mt->ne;i++)
  {    
    ipp = Mt->elements[i].ipp[0][0];
    nip = Mt->give_totnip(i);
    selelems.presence_id(i, ir);

    if (ir < 0)
      n = 0;
    else
      n = selother[ir].give_nselcomp(ip[ipp].ncompeqother);
    for (j=0; j<nip; j++)
    {
      fwrite(&ipp, sizeof(ipp), 1, aux);
      fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
      fwrite(&ip[ipp].ncompstr, sizeof(ip[ipp].ncompstr), 1, aux);
      fwrite(&n, sizeof(n), 1, aux);
      if (ir < 0)
        ip[ipp].save_data_bin(aux,Mb->nlc, selno);
      else
        ip[ipp].save_data_bin(aux,Mb->nlc, selother[ir]);

      // writing TRFEL quantities
      if (nonmechq)
      {
        for (k=0; k<nnmq; k++)
          fwrite (&nonmechq[k][ipp], sizeof(nonmechq[k][ipp]), 1, aux);
      }
      ipp++;
    }
  }
}



/**
  Function saves data from regular integration points into several backup files
  in binary format.
   
  @param[in] ni       - time step id
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka 9.2008
*/
void mechmat::save_intpoints_bin (long ni, sel &selelems, sel *selother)
{
  long i,j,k,ir,n,nip,ipp;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.%ld.strain.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fwrite(&i, sizeof(i), 1, aux);
    fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
    fwrite(&ip[i].ncompstr, sizeof(ip[i].ncompstr), 1, aux);
    fwrite (ip[i].strain, sizeof(*ip[i].strain),Mb->nlc*ip[i].ncompstr, aux);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.stress.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fwrite(&i, sizeof(i), 1, aux);
    fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
    fwrite(&ip[i].ncompstr, sizeof(ip[i].ncompstr), 1, aux);
    fwrite(ip[i].stress, sizeof(*ip[i].stress),Mb->nlc*ip[i].ncompstr, aux);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.other.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<Mt->ne;i++)
  {    
    ipp = Mt->elements[i].ipp[0][0];
    nip = Mt->give_totnip(i);
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
      for (k=0; k<ip[ipp].ncompeqother; k++)
      {
        if (selother[ir].presence_id(k))
	  fwrite (ip[ipp].eqother+k, sizeof(*ip[ipp].eqother), 1, aux);
      }
      ipp++;
    }
  }
  fclose(aux);

  if (nonmechq)  
  {
    // writing TRFEL quantities
    sprintf(name, "%s.%ld.trfquant.bac",Mp->hdbcont.hdbnames, ni);
    aux = fopen(name,"wb");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    // write table of indices of non-mechanical quantities
    fwrite (nmqid, sizeof(*nmqid), tnknmq, aux);

    for (ipp=0;ipp<tnip;ipp++)
    {
      fwrite(&ipp, sizeof(ipp), 1, aux);
      for (j=0; j<nnmq; j++)
        fwrite(&nonmechq[j][ipp], sizeof(nonmechq[j][ipp]), 1, aux);
    }
    fclose(aux);
  }
}



/**
  Function restores data from text backup file into regular integration points.
   
  @param[in] aux - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  Created by JK, 19.9.2004
  Rewritten by TKo 9.2008
  Modified by TKo, 7.10.2013
*/
void mechmat::restore_intpoints_txt (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long i,j,n,ir,ipp;
  long nlc, ncompstr;
  sel selno; // selection with no selected components
  
  selno.n=1;

  if (nonmechq) // read table of indices of non-mechanical quantities
  {
    for (i=0; i<tnknmq; i++)
      fscanf(aux, "%ld", &nmqid[i]);
  }

  for (i=0;i<tnip;i++)
  {
    fscanf (aux,"%ld %ld %ld %ld",&ipp, &nlc, &ncompstr, &n);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != ip[ipp].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);
    if (ir < 0)
      ip[ipp].restore_data_txt(aux, nlc, n, selno, NULL);
    else
      ip[ipp].restore_data_txt(aux, nlc, n, selother[ir], selid[ir]);

    // restoring TRFEL quantities
    if (nonmechq)
    {
      for (j=0; j<nnmq; j++)
        fscanf(aux,"%le", &nonmechq[j][ipp]);
    }
  }
}



/**
  Function restores data from several text backup files into regular integration points.
   
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  @return The function does not return anything.

  Created by JK, 19.9.2004
  Rewritten by TKo 9.2008
*/
void mechmat::restore_intpoints_txt (sel &selelems, sel *selother, long **selid)
{
  long i,j,k,n,ir,ik,is,ipp;
  long nlc, ncompstr;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of strain arrays
  sprintf(name, "%s.strain.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fscanf (aux,"%ld %ld %ld",&ipp, &nlc, &ncompstr);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != ip[ipp].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Mb->nlc*ip[i].ncompstr; j++)
      fscanf (aux, "%le", ip[i].strain+j);
  }
  // restoring of stress arrays
  sprintf(name, "%s.stress.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++){
    fscanf (aux,"%ld %ld %ld",&ipp, &nlc, &ncompstr);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != ip[ipp].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0;j<nlc*ncompstr;j++)
      fscanf (aux, "%le", ip[i].stress+j);
  }
  // restoring of other arrays
  sprintf(name, "%s.other.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (i=0;i<tnip;i++)
  {
    fscanf (aux,"%ld %ld",&ipp, &n);
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

  // restoring TRFEL quantities
  if (nonmechq)
  {
    sprintf(name, "%s.trfquant.bac",Mp->hdbcont.hdbnamer);
    aux = fopen(name,"rt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }

    // read table of indices of non-mechanical quantities
    for (i=0; i<tnknmq; i++)
      fscanf(aux, "%ld ", &nmqid[i]);

    for (i=0;i<tnip;i++)
    {
      fscanf (aux, "%ld", &ipp);
      if ((ipp < 0) || (ipp >= tnip))
      {
        print_err("invalid integration point number", __FILE__, __LINE__, __func__);
        abort();
      }
      for (j=0; j<nnmq; j++)
        fscanf (aux,"%le", &nonmechq[j][ipp]);
    }
    fclose(aux);
  }
}



/**
  Function restores data from binary backup file into regular integration points.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  @return The function does not return anything.

  Created by Tomas Koudelka 9.2008
*/
void mechmat::restore_intpoints_bin (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long i, j, ipp, nlc, ncompstr, n;
  long ir;
  sel selno; // selection with no selected components
  
  if (nonmechq)  // read table of indices of non-mechanical quantities
  {
    fread (nmqid, sizeof(*nmqid), tnknmq, aux);
    // the following  lines are temporary bugfix of nonmechanical quantity storage
    //long *aaa = new long[sizeof(nmqid)];
    //fread (aaa, sizeof(*aaa), sizeof(nmqid)-tnknmq, aux);
    //delete [] aaa;
    // end of bug fix
  }

  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    fread(&n, sizeof(n), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != ip[ipp].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(elip[ipp], ir);
    if (ir < 0)
      ip[ipp].restore_data_bin(aux, Mb->nlc, n, selno, NULL);
    else
      ip[ipp].restore_data_bin(aux, Mb->nlc, n, selother[ir], selid[ir]);

    // restoring TRFEL quantities
    if (nonmechq)
    {
      for (j=0; j<nnmq; j++)
        fread(&nonmechq[j][ipp], sizeof(nonmechq[j][ipp]), 1, aux);
    }
  }
}



/**
  Function restores data from several binary backup files into regular integration points.
   
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled
  
  @return The function does not return anything.
   
  Created by Tomas Koudelka, 9.2008
*/
void mechmat::restore_intpoints_bin (sel &selelems, sel *selother, long **selid)
{
  long i,k,n,ir,ik,is,ipp;
  long nlc, ncompstr;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of strain arrays
  sprintf(name, "%s.strain.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
      {
	print_err("invalid integration point number", __FILE__, __LINE__, __func__);
	abort();
      }
    if (nlc != Mb->nlc)
      {
	print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
	abort();
      }
    if (ncompstr != ip[ipp].ncompstr)
      {
	print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
	abort();
      }
    fread (ip[ipp].strain, sizeof(*ip[ipp].strain),Mb->nlc*ip[ipp].ncompstr, aux);
  }
  // restoring of stress arrays
  sprintf(name, "%s.stress.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
  for (i=0;i<tnip;i++){
    fread(&ipp, sizeof(ipp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
      {
	print_err("invalid integration point number", __FILE__, __LINE__, __func__);
	abort();
      }
    if (nlc != Mb->nlc)
      {
	print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
	abort();
      }
    if (ncompstr != ip[ipp].ncompstr)
      {
	print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
	abort();
      }
    fread (ip[ipp].stress, sizeof(*ip[ipp].stress),Mb->nlc*ip[ipp].ncompstr, aux);
  }
  // restoring of other arrays
  sprintf(name, "%s.other.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
  for (i=0;i<tnip;i++){
    fread (&ipp, sizeof(ipp), 1, aux);
    fread (&n, sizeof(n), 1, aux);
    if ((ipp < 0) || (ipp >= tnip))
      {
	print_err("invalid integration point number", __FILE__, __LINE__, __func__);
	abort();
      }
    selelems.presence_id(elip[ipp], ir);    
    for (k=0;k<n;k++)
      {
	fread (&tmp, sizeof(tmp), 1, aux);
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
  // restoring TRFEL quantities
  if (nonmechq)
  {
    sprintf(name, "%s.trfquant.bac",Mp->hdbcont.hdbnamer);
    aux = fopen(name,"rb");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
      
    // read table of indices of non-mechanical quantities
    fread (nmqid, sizeof(*nmqid), tnknmq, aux);
      
    for (i=0;i<tnip;i++)
    {
      fread (&ipp, sizeof(ipp), 1, aux);
      if ((ipp < 0) || (ipp >= tnip))
      {
        print_err("invalid integration point number", __FILE__, __LINE__, __func__);
        abort();
      }
      for (k=0; k<nnmq; k++)
        fread(&nonmechq[k][ipp], sizeof(nonmechq[k][ipp]), 1, aux);
    }
    fclose(aux);
  }
}



// ***************************************************
// ***************************************************
//
// Funtions for backup of auxiliary int. point content
//
// ***************************************************
// ***************************************************



/**
  Function saves data from auxiliary integration points into backup file
  in text format.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled

  @return The function does not return anything.
   
  Created by TKo, 19.12.2017
*/
void mechmat::save_auxintpoints_txt (FILE *aux, sel &selelems, sel *selother)
{
  long ir, app, eid, nsc;
  sel selno; // selection with no selected components
  
  selno.n=1;
  fprintf(aux, "\n");

  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];

    selelems.presence_id(eid, ir);
    if (ir < 0)
      nsc = 0;    
    else
      nsc = selother[ir].give_nselcomp(aip[app].ncompeqother);

    fprintf (aux, "%ld %ld %ld %ld\n", app, Mb->nlc, aip[app].ncompstr, nsc);
    if (ir < 0)
      aip[app].save_data_txt(aux, Mb->nlc, selno);
    else
      aip[app].save_data_txt(aux, Mb->nlc, selother[ir]);

    // writing TRFEL quantities is not needed, their values can be recalculated after restorage
  }
}



/**
  Function saves data from auxiliary integration points into several backup files
  in text format.
   
  @param[in] ni       - time step id
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by TKo, 20.12.2017
*/
void mechmat::save_auxintpoints_txt (long ni, sel &selelems, sel *selother)
{
  long j, k, ir, eid, nsc, app;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Mp->hdbcont.prec;

  sprintf(name, "%s.%ld.astrain.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }  
  for (app=0; app<tnaip; app++)
  {
    fprintf (aux,"%ld %ld %ld\n", app, Mb->nlc, aip[app].ncompstr);
    for (j=0; j<Mb->nlc*aip[app].ncompstr; j++)
      fprintf (aux, "%.*le\n", prec, aip[app].strain[j]);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.astress.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    fprintf (aux, "%ld %ld %ld\n", app, Mb->nlc, aip[app].ncompstr);
    for (j=0; j<Mb->nlc*aip[app].ncompstr; j++)
      fprintf (aux, "%.*le\n", prec, aip[app].stress[j]);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.aother.bac", Mp->hdbcont.hdbnames, ni);
  aux = fopen(name, "wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];
    selelems.presence_id(eid, ir);    
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
 
  // writing TRFEL quantities (nonmechq) is not needed, their values can be recalculated after restorage
}



/**
  Function saves data from auxiliary integration points into backup file
  in binary format.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 20.12.2017
*/
void mechmat::save_auxintpoints_bin (FILE *aux, sel &selelems, sel *selother)
{
  long ir, eid, nsc, app;
  sel selno; // selection with no selected components
  
  selno.n=1;

  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];
    selelems.presence_id(eid, ir);

    if (ir < 0)
      nsc = 0;
    else
      nsc = selother[ir].give_nselcomp(aip[app].ncompeqother);
    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
    fwrite(&aip[app].ncompstr, sizeof(aip[app].ncompstr), 1, aux);
    fwrite(&nsc, sizeof(nsc), 1, aux);
    if (ir < 0)
      aip[app].save_data_bin(aux, Mb->nlc, selno);
    else
      aip[app].save_data_bin(aux, Mb->nlc, selother[ir]);

    // writing TRFEL quantities is not needed, their values can be recalculated after restorage
  }
}



/**
  Function saves data from auxiliary integration points into several backup files
  in binary format.
   
  @param[in] ni       - time step id
  @param[in] selelems - selection of elements whose eqother array will be saved
  @param[in] selother - selection of components of saved eqother array which will be saved on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
   
  @return The function does not return anything.

  Created by Tomas Koudelka 20.12.2017
*/
void mechmat::save_auxintpoints_bin (long ni, sel &selelems, sel *selother)
{
  long k, ir, eid, nsc, app;
  sel selno; // selection with no selected components
  selno.n=1;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.%ld.astrain.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
    fwrite(&aip[app].ncompstr, sizeof(aip[app].ncompstr), 1, aux);
    fwrite(aip[app].strain, sizeof(*aip[app].strain), Mb->nlc*aip[app].ncompstr, aux);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.astress.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    fwrite(&app, sizeof(app), 1, aux);
    fwrite(&Mb->nlc, sizeof(Mb->nlc), 1, aux);
    fwrite(&aip[app].ncompstr, sizeof(aip[app].ncompstr), 1, aux);
    fwrite(aip[app].stress, sizeof(*aip[app].stress), Mb->nlc*aip[app].ncompstr, aux);
  }
  fclose(aux);

  sprintf(name, "%s.%ld.aother.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    eid = elaip[app];
    selelems.presence_id(eid, ir);    
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
        fwrite (aip[app].eqother+k, sizeof(*aip[app].eqother), 1, aux);
    }
  }
  fclose(aux);

  // writing TRFEL quantities is not needed, their values can be recalculated after restorage
}



/**
  Function restores data from text backup file into auxiliary integration points.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  Created by Tomas Koudelka, 20.12.2017
*/
void mechmat::restore_auxintpoints_txt (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long ir, eid, app, tapp;
  long nlc, ncompstr, ncompother;
  sel selno; // selection with no selected components
  
  selno.n=1;

  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];

    fscanf (aux,"%ld %ld %ld %ld",&tapp, &nlc, &ncompstr, &ncompother);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);
    if (ir < 0)
      aip[app].restore_data_txt(aux, nlc, ncompother, selno, NULL);
    else
      aip[app].restore_data_txt(aux, nlc, ncompother, selother[ir], selid[ir]);

    // restoring TRFEL quantities is not needed, they are recalculated from other stored values
  }
}



/**
  Function restores data from several text backup files into auxiliary integration points.
   
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 20.12.2017
*/
void mechmat::restore_auxintpoints_txt (sel &selelems, sel *selother, long **selid)
{
  long j, k, ir, ik, is, eid, app, tapp;
  long nlc, ncompstr, ncompother;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of strain arrays
  sprintf(name, "%s.astrain.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    fscanf(aux, "%ld %ld %ld",&tapp, &nlc, &ncompstr);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<Mb->nlc*aip[app].ncompstr; j++)
      fscanf(aux, "%le", aip[app].strain+j);
  }
  // restoring of stress arrays
  sprintf(name, "%s.astress.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    fscanf(aux, "%ld %ld %ld", &tapp, &nlc, &ncompstr);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    for (j=0; j<nlc*ncompstr; j++)
      fscanf(aux, "%le", aip[app].stress+j);
  }
  // restoring of other arrays
  sprintf(name, "%s.aother.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {
    eid = elaip[app];

    fscanf (aux, "%ld %ld", &tapp, &ncompother);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);    
    for (k=0; k<ncompother; k++)
    {
      fscanf(aux, "%le", &tmp);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k,ik))
      {
        is = selid[ir][ik]+k-selother[ir].id1[ik]; // storage index in eqother array 
        if (is >= aip[app].ncompeqother)
          print_err("invalid index for eqother restoring is required in app=%ld", __FILE__, __LINE__, __func__, app);
        else
          aip[app].eqother[is] = tmp;
      }
    }
  }
  fclose(aux);
}



/**
  Function restores data from binary backup file into auxiliary integration points.
   
  @param[in] aux      - pointer to auxiliary file
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled for each selection of elements
   
  @return The function does not return anything.

  Created by Tomas Koudelka 20.12.2017
*/
void mechmat::restore_auxintpoints_bin (FILE *aux, sel &selelems, sel *selother, long **selid)
{
  long nlc, ncompstr, ncompother, eid, app, tapp;
  long ir;
  sel selno; // selection with no selected components
  
  if (nonmechq)  // read table of indices of non-mechanical quantities
    fread (nmqid, sizeof(*nmqid), sizeof(nmqid), aux);

  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];

    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    fread(&ncompother, sizeof(ncompother), 1, aux);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);
    if (ir < 0)
      aip[app].restore_data_bin(aux, Mb->nlc, ncompother, selno, NULL);
    else
      aip[app].restore_data_bin(aux, Mb->nlc, ncompother, selother[ir], selid[ir]);
  }
}



/**
  Function restores data from several binary backup files into auxliary integration points.
   
  @param[in] selelems - selection of elements whose eqother array will be restored
  @param[in] selother - selection of components of saved eqother array which will be resored on selected elements
                        for each range or list item in selelems an individual selection of eqother components
                        is enabled
  @param[in] selid    - array of indices of positions in eqother array to which will be 
                        restored selected saved eqother components
                        for each range or list item in selother an individual selection of position 
                        of eqother components is enabled
  
  @return The function does not return anything.
   
  Created by Tomas Koudelka, 20.12.2017
*/
void mechmat::restore_auxintpoints_bin (sel &selelems, sel *selother, long **selid)
{
  long k, ir, ik, is, eid, app, tapp;
  long nlc, ncompstr, ncompother;
  double tmp;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  // restoring of strain arrays
  sprintf(name, "%s.astrain.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {    
    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    if (tapp != app)
    {
      print_err("invalid auxiliary integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    fread (aip[app].strain, sizeof(*aip[app].strain),Mb->nlc*aip[app].ncompstr, aux);
  }
  // restoring of stress arrays
  sprintf(name, "%s.astress.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {    
    fread(&tapp, sizeof(tapp), 1, aux);
    fread(&nlc, sizeof(nlc), 1, aux);
    fread(&ncompstr, sizeof(ncompstr), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    if (nlc != Mb->nlc)
    {
      print_err("invalid number of load cases", __FILE__, __LINE__, __func__);
      abort();
    }
    if (ncompstr != aip[app].ncompstr)
    {
      print_err("incompatible number of stress/strain components", __FILE__, __LINE__, __func__);
      abort();
    }
    fread(aip[app].stress, sizeof(*aip[app].stress),Mb->nlc*aip[app].ncompstr, aux);
  }
  // restoring of other arrays
  sprintf(name, "%s.aother.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  for (app=0; app<tnaip; app++)
  {    
    eid = elaip[app];

    fread (&tapp, sizeof(tapp), 1, aux);
    fread (&ncompother, sizeof(ncompother), 1, aux);
    if (tapp != app)
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    selelems.presence_id(eid, ir);    
    for (k=0;k<ncompother;k++)
    {
      fread (&tmp, sizeof(tmp), 1, aux);
      if (ir < 0)
        continue;
      if (selother[ir].presence_id(k,ik))
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
