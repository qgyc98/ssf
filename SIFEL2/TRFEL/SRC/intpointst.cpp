#include "intpointst.h"
#include "globalt.h"
#include <string.h>

intpointst::intpointst ()
{
  //  number of components in gradients
  ncompgrad=0;
  //  number of components in array other
  ncompother=0;
  //  number of components in array eqother
  ncompeqother=0;
  
  //  type of material
  tm = (mattypet) 0;
  //  material id
  idm=0;

  //  array of influence of particular gradients
  infl=NULL;
  cap = NULL;
  
  //  array of actual values
  av=NULL;
  //  array of previous values (values from previous time step)
  pv=NULL;
  //  array of gradients
  grad=NULL;
  //  array of fluxes
  fluxes=NULL;
  //  auxiliary array other
  other = NULL;
  //  auxiliary array eqother
  eqother = NULL;
}

intpointst::~intpointst (void)
{
  long i;
  delete [] av;  
  delete [] pv;
  for (i=0; i<Tp->ntm; i++)
  {
    if (infl != NULL)
      delete [] infl[i];
    if (grad != NULL)
      delete [] grad[i];
    if (fluxes != NULL)
      delete [] fluxes[i]; 
  }
  delete [] grad;
  delete [] fluxes;
  delete [] infl;
  delete [] cap;

  delete [] other;   
  delete [] eqother;
}

/**
   function reads data connected with integration point
   
   @param in - input file
*/
void intpointst::read (FILE *in)
{
  fscanf (in,"%d %ld",(int*)&tm,&idm);
  idm--;
}


/**
   function allocates memory for actual values, previous values, axiliary values,
   equilibriated auxiliary values, gradients and fluxes
   
*/
void intpointst::alloc ()
{
  long i,ntm;

  //  number of transported media
  ntm = Tp->ntm;
  //  number of components of gradients
  ncompgrad = Tp->gdim;
  //  number of components of array other
  ncompother = Tm->givencompother(); 
  //  number of components of array eqother
  ncompeqother = Tm->givencompeqothermat(tm);

  if (ncompother > 0){
    other = new double [ncompother];
  }
  
  if (ncompeqother > 0){
    eqother = new double [ncompeqother];
  }
  
  //  array of influence of particular gradients
  //infl = new long* [ntm];
  
  //  array for actual values
  av = new double [ntm];
  //  array for previous values (values from the previous time step)
  pv = new double [ntm];
  //  array for gradients
  grad = new double* [ntm];
  //  array for fluxes
  fluxes = new double* [ntm];
  cap = new double [ntm];

  for (i=0;i<ntm;i++){
    //infl[i] = new long [ntm];
    grad[i] = new double [ncompgrad];
    fluxes[i] = new double [ncompgrad];
  }
  
  clean ();
}


/**
   function cleans all arrays defined on integration point
   
   JK
*/
void intpointst::clean ()
{
  long i,j,ntm;
  
  //  number of transported media
  ntm = Tp->ntm;
  
  for (i=0;i<ntm;i++){
    av[i]=0.0;
    pv[i]=0.0;
    cap[i]=0.0;
    for (j=0;j<ncompgrad;j++){
      //infl[i][j]=1;
      grad[i][j]=0.0;
      fluxes[i][j]=0.0;
    }
  }
  
  for (i=0;i<ncompother;i++){
    other[i]=0.0;
  }
  
  for (i=0;i<ncompeqother;i++){
    eqother[i]=0.0;
  }
}



/**
  Function stores %vector of lcid-th quantity gradient to the 
  given integration point.

  @param lcid - load case id
  @param gradv - array with %vector of gradients
   
  Created by Tomas Koudelka, 3.2.2014
*/
void intpointst::storegrad (long lcid, double *gradv)
{
  long i;
  for (i=0; i<ncompgrad; i++)
    grad[lcid][i]=gradv[i];
}



/**
  The function stores copy of the integration point content in ip parameter. 

  @param ip - integration point copy (output)
  @param ntm - number of transported media

  @return The function stores copy in the parameter ip.

  Created by Tomas Koudelka, 10.12.2013
*/
void intpointst::copy(intpointst &ip, long ntm, long realloc)
{
  long i;
  // material type
  tm = ip.tm;

  // material id
  idm = ip.idm;

  // number of component of flux/grad array
  ncompgrad = ip.ncompgrad;

  // number of component of other array
  ncompother = ip.ncompother;

  // number of component of eqother array
  ncompeqother = ip.ncompeqother;
  
  // array of influence of particular gradients
  if (realloc)
  {
    if (infl)
    {
      for(i=0; i<ntm; i++)
        delete [] infl[i];
      delete [] infl;
    }
    if (ip.infl)
    {
      infl = new long*[ntm];
      for(i=0; i<ntm; i++)
        infl[i] = new long[ntm];
    }
    else
      infl = NULL;
  }
  if (ip.infl)
  {
    for(i=0; i<ntm; i++)
      memcpy(infl[i], ip.infl[i], sizeof(*infl[i])*ntm);
  }

  //  array of actual values
  if (realloc)
  {
    if (av)
      delete [] av;
    if (ip.av)
      av = new double[ntm];
    else
      av = NULL;
  }
  if (ip.av)
    memcpy(av, ip.av, sizeof(*av)*ntm);

  // array of previous values
  if (realloc)
  {
    if (pv)
      delete [] pv;
    if (ip.pv)
      pv = new double[ntm];
    else
      pv = NULL;
  }
  if (ip.pv)
    memcpy(pv, ip.pv, sizeof(*pv)*ntm);

  // array of gradients
  if (realloc)
  {
    if (grad)
    {
      for(i=0; i<ntm; i++)
        delete [] grad [i];
      delete [] grad;
    }
    if (ip.grad)
    {    
      grad = new double*[ntm];
      for (i=0; i<ntm; i++)
        grad[i] = new double[ncompgrad];
    }
    else
      grad = NULL;
  }
  if (ip.grad)
  {
    for (i=0; i<ntm; i++)
      memcpy(grad[i], ip.grad[i], sizeof(*grad[i])*ncompgrad);
  }

  // array of fluxes
  if (realloc)
  {
    if (fluxes)
    {
      for(i=0; i<ntm; i++)
        delete [] fluxes [i];
      delete [] fluxes;
    }
    if (ip.fluxes)
    {
      fluxes = new double*[ntm];
      for (i=0; i<ntm; i++)
        fluxes[i] = new double[ncompgrad];
    }
    else
      fluxes = NULL;
  }
  if (ip.fluxes)
  {
    for (i=0; i<ntm; i++)
      memcpy(fluxes[i], ip.fluxes[i], sizeof(*fluxes[i])*ncompgrad);
  }

  // other components
  if (realloc)
  {
    if (other)
      delete [] other;
    if (ip.other)
      other = new double[ncompother];
    else
      other = NULL;
  }
  if (ip.other)
    memcpy(other, ip.other, sizeof(*other)*ncompother);

  // eqother components
  if (realloc)
  {
    if (eqother)
      delete [] eqother;
    if (ip.eqother)
      eqother = new double[ncompeqother];
    else
      eqother = NULL;
  }
  if (ip.eqother)
    memcpy(eqother, ip.eqother, sizeof(*eqother)*ncompeqother);
}



/**
   function saves all arrays into the auxiliary file
   
   @param aux - pointer to auxiliary file
   @param selother - selection of components of eqother array which will be saved
   
   JK, 21.9.2004
   TKr 12.2008
*/
void intpointst::save_data_txt (FILE *aux,sel &selother)
{
  long i,j,ntm,gdim;
  int prec = (int)Tp->hdbcont.prec;
  
  //  number of transported media
  ntm = Tp->ntm;
  //  number of geometrical dimensions
  gdim = Tp->gdim;
  
  for (i=0;i<ntm;i++){
    fprintf (aux,"%e\n",av[i]);
  }
  for (i=0;i<ntm;i++){
    fprintf (aux,"%e\n",pv[i]);
  }
  for (i=0;i<ntm;i++){
    for (j=0;j<gdim;j++){
      fprintf (aux,"%e\n",grad[i][j]);
    }
  }
  for (i=0;i<ntm;i++){
    for (j=0;j<gdim;j++){
      fprintf (aux,"%e\n",fluxes[i][j]);
    }
  }
  
  for (i=0;i<ncompeqother;i++)
    {
      if (selother.presence_id(i))
	fprintf (aux,"%.*le\n",prec,eqother[i]);
    }  
}

/**
   function restores all arrays from the auxiliary file
   
   @param aux - pointer to auxiliary file
   @param nlc - number of load cases
   @param ncompo - number of saved eqother components
   @param selother - selection of components of eqother array to which will be other values restored
   @param selid - array of indeces 
   
   JK, 21.9.2004
   TKr 12.2008
*/
void intpointst::restore_data_txt (FILE *aux,long ncompo, sel &selother, long *selid)
{
  long i,j,ntm,gdim,ik,is;
  double tmp;
  
  //  number of transported media
  ntm = Tp->ntm;
  //  number of geometrical dimesnions
  gdim = Tp->gdim;
  
  for (i=0;i<ntm;i++){
    fscanf (aux,"%le\n",&av[i]);
  }
  for (i=0;i<ntm;i++){
    fscanf (aux,"%le\n",&pv[i]);
  }
  for (i=0;i<ntm;i++){
    for (j=0;j<gdim;j++){
      fscanf (aux,"%le\n",&grad[i][j]);
    }
  }
  for (i=0;i<ntm;i++){
    for (j=0;j<gdim;j++){
      fscanf (aux,"%le\n",&fluxes[i][j]);
    }
  }
  
  for (i=0;i<ncompo;i++)
    {
      fscanf (aux,"%le",&tmp);
      if (selother.presence_id(i,ik))
	{
	  is = selid[ik]+i-selother.id1[ik]; // storage index in eqother array 
	  if (is >= ncompeqother)
	    print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
	  else
	    eqother[is] = tmp;
	}  
    }
}

/**
   function moves array of actual values to array of previous values
   
   JK, 29.5.2007
*/
void intpointst::actual_previous_change ()
{
  long i,ntm;
  
  //  number of transported media
  ntm = Tp->ntm;
  
  for (i=0;i<ntm;i++){
    pv[i]=av[i];
  }
}



/**
   function saves all arrays into the auxiliary file in binary format
   
   @param aux - pointer to auxiliary file
   @param nlc - number of load cases
   @param selother - selection of components of eqother array which will be saved
   
   TKo 9.2008
*/
void intpointst::save_data_bin (FILE *aux,sel &selother)
{
  long i;
  
  //fwrite (strain, sizeof(*strain),nlc*ncompstr, aux);
  //fwrite (stress, sizeof(*stress),nlc*ncompstr, aux);
  for (i=0;i<ncompeqother;i++)
  {
    if (selother.presence_id(i))
      fwrite (eqother+i, sizeof(*eqother), 1, aux);
  }
}


/**
   function restores all arrays from the auxiliary file in binary format
   
   @param aux - pointer to auxiliary file
   @param nlc - number of load cases
   @param selother - selection of components of eqother array to which will be other values restored
   @param selid    - array of indeces of positions in eqother array to which will be 
                     restored selected saved eqother components
   
   TKo 9.2008
*/
void intpointst::restore_data_bin (FILE *aux,long ncompo, sel &selother, long *selid)
{
  long i, ik, is;
  double tmp;
  
  //fread (strain, sizeof(*strain), nlc*ncompstr, aux);
  //fread (stress, sizeof(*stress), nlc*ncompstr, aux);
  for (i=0;i<ncompo;i++)
  {
    fread (&tmp, sizeof(tmp), 1, aux);
    if (selother.presence_id(i,ik))
    {
      is = selid[ik]+i-selother.id1[ik]; // storage index in eqother array 
      if (is >= ncompeqother)
        print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
      else
        eqother[is] = tmp;
    }
  }
}
