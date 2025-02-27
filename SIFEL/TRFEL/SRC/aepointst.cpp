#include "aepointst.h"
#include "globalt.h"
#include "vector.h"


aepointst::aepointst ()
{
  tape=NULL;  ptape=NULL;
}
aepointst::~aepointst ()
{
  delete [] tape;
  delete [] ptape;
}

void aepointst::read (FILE *in)
{
  long i,j,k;
  
  tape = new long [Tt->ne];
  ptape = new long [Tt->ne];
  
  for (i=0;i<Tt->ne;i++){
    ptape[i]=-2;
    fscanf (in,"%ld",tape+i);
    if (tape[i]==userdefinedt){
      fscanf (in,"%ld",ptape+i);
      ptape[i]--;
    }
  }
  
  //  number of point sets
  fscanf (in,"%ld",&nudsets);
  udpc = new double** [nudsets];
  udpa = new long* [nudsets];
  for (i=0;i<nudsets;i++){
    udpa[i] = new long [2];
    //  number of points, number of components
    fscanf (in,"%ld %ld",&udpa[i][0],&udpa[i][1]);
    udpc[i] = new double* [udpa[i][0]];
    for (j=0;j<udpa[i][0];j++){
      udpc[i][j] = new double [udpa[i][1]];
      for (k=0;k<udpa[i][1];k++){
        fscanf (in,"%lf",&udpc[i][j][k]);
      }
    }
  }
  
  //  number of elements with required transformation
  fscanf (in,"%ld",&net);
  ent = new long [net];
  npt = new long [net];
  pt = new long* [net];
  plcs = new long* [net];
  for (i=0;i<net;i++){
    fscanf (in,"%ld %ld",ent+i,npt+i);
    ent[i]--;
    pt[i] = new long [npt[i]];
    plcs[i] = new long [npt[i]];
    for (j=0;j<npt[i];j++){
      fscanf (in,"%ld %ld",&pt[i][j],&plcs[i][j]);
      pt[i][j]--;  plcs[i][j]--;
    }
  }
  
  //  number of local coordinate systems
  if (net>0){
    fscanf (in,"%ld",&nlcs);
    lcs = new double* [nlcs];
    nclcs = new long [nlcs];
    for (i=0;i<nlcs;i++){
      fscanf (in,"%ld",nclcs+i);
      lcs[i] = new double [nclcs[i]];
      for (j=0;j<nclcs[i];j++){
        fscanf (in,"%lf",&lcs[i][j]);
      }
    }
  }
  
}



void aepointst::init (strastret sst)
{
  long i;

  tape = new long [Tt->ne];
  ptape = new long [Tt->ne];

  if (sst==grad)
  {
    for (i=0; i<Tt->ne; i++)
    {
      if (Outdt->eo.selegrad.presence_id(i) || Tp->gradcomp)
        tape[i] = 1;
    }
  }
  
  if (sst==flux)
  {
    for (i=0; i<Tt->ne; i++)
    {
      if (Outdt->eo.seleflux.presence_id(i) || Tp->fluxcomp)
        tape[i] = 1;
    }
  }
  nudsets = 0;
  net = 0;
  nlcs =0;
}



/**
   function returns number of auxiliary points on element
   (number of user defined points, where values will be computed)
   
   @param eid - element id
   
   19.5.2002
*/
long aepointst::give_naep (long eid)
{
  long naep;
  naep = udpa[eid][0];
  return naep;
}
/**
   function returns number of components of quantity computed in auxiliary points
   e.g. number of components of strain tensor
   
   @param eid - element id
   
   19.5.2002
*/
long aepointst::give_ncomp (long eid)
{
  long ncomp;
  ncomp = udpa[eid][1];
  return ncomp;
}

/**
   function returns number of set of auxiliary element points on element
   
   @param eid - element id
   
   19.5.2002
*/
long aepointst::give_sid (long eid)
{
  long sid;
  sid = ptape[eid];
  return sid;
}

/**
   function returns coordinates of auxiliary element points
   
   @param sid - number of set of auxiliary points
   @param pid - number of point in reqiered set
   @param coord - array containing coordinates
   
   19.5.2002
*/
void aepointst::give_aepcoord (long sid,long pid,vector &coord)
{
  long i;
  for (i=0;i<udpa[sid][1];i++){
    coord[i]=udpc[sid][pid][i];
  }
}

/**
   function allocates array ev
   
   @param nlc - number of load cases
   
   19.5.2002
*/
void aepointst::alloc (long nlc)
{
  long i,j,naep,ncomp;
  
  ev = new double** [Tt->ne];
  
  for (i=0;i<Tt->ne;i++){
    if (tape[i]==userdefinedt){
      //  number of auxiliary points on element
      naep=udpa[ptape[i]][0];
      ev[i] = new double* [naep];
      //  number of components
      ncomp=Tt->give_ncomp(i);
      for (j=0;j<naep;j++){
        ev[i][j] = new double [nlc*ncomp];
      }
    }
  }
}

/**
   function stores evaluated values
   
   @param lcid - load case id
   @param eid - element id
   @param pid - auxiliary point id
   @param val - array containing values
   
   JK, 22.2.2002
*/
void aepointst::storevalues (long lcid,long eid,long pid,vector &val)
{
  long i,ncomp;

  if (tape[eid]!=userdefinedt){
    fprintf (stderr,"\n\n wrong application of function aepointst::storevalues (%s, line %d).\n",__FILE__,__LINE__);
  }
  
  ncomp=val.n;
  for (i=0;i<ncomp;i++){
    ev[eid][pid][i]=val[lcid*ncomp+i];
  }
  
}


/**
   function transforms values in auxiliary points to local coordinate systems
   
   @param tt - type of transformation
   
   tt=0 - for stress
   tt=1 - for strain
   
   JK, 22.2.2002
*/
void aepointst::transformvalues (long tt)
{
  tt = tt;
  /*
  long i,j,ii,jj,kk,nc,ncomp;
  vector quant;
  matrix tmat;
  
  for (i=0;i<net;i++){
    //  number of required element
    ii=ent[i];
    //  number of components of required quantity
    ncomp=Tt->give_ncomp(ii);
    allocv (ncomp,quant);

    //  loop over required points on element
    for (j=0;j<npt[i];j++){
      
      //  number of required point
      jj=pt[i][j];
      //  number of required local coordinate system
      kk=plcs[i][j];
      //  number of components of local coordinate system
      nc=nclcs[kk];
      
      //  required data extraction
      restorevalues (ii,jj,quant.a,tt);
      
      //  transformation matrix assembling
      if (nc==4){
        allocm (2,2,tmat);
        tmat[0][0]=lcs[kk][0];  tmat[0][1]=lcs[kk][2];
        tmat[1][0]=lcs[kk][1];  tmat[1][1]=lcs[kk][3];
      }
      if (nc==9){
        allocm (3,3,tmat);
        tmat[0][0]=lcs[kk][0];  tmat[0][1]=lcs[kk][3];  tmat[0][2]=lcs[kk][6];
        tmat[1][0]=lcs[kk][1];  tmat[1][1]=lcs[kk][4];  tmat[1][2]=lcs[kk][7];
        tmat[2][0]=lcs[kk][2];  tmat[2][1]=lcs[kk][5];  tmat[2][2]=lcs[kk][8];
      }
      
      if (tt==1)
        //engtens (quant);
      
      glob_loc_tens_trans (quant,tmat);
      
      if (tt==1)
        //tenseng (quant);
      
      //  required data storing
      storevalues (ii,jj,quant.a,tt);
      
      destrm (tmat);
    }
    destrv (quant);
  }
  */
}

