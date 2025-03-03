#include "springel.h"
#include "global.h"
#include "globmat.h"
#include "element.h"
#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include <math.h>

springel::springel (void)
{
  long i;

  nne=1;  ndofe=-1;  tncomp=1;  napfun=0;
  ssst = bar;

  nb=1;

  ncomp = new long [nb];
  ncomp[0]=1;

  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }

  nip[0][0]=1;
  tnip=1;
  
  intordsm[0][0]=1;
}

springel::~springel (void)
{
  long i;

  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;

  delete [] cncomp;
  delete [] ncomp;
  delete [] nip;
}

void springel::eleminit (long eid)
{
  long ii,jj;
  Mt->elements[eid].nb=nb;
  Mt->elements[eid].intordsm = new long* [nb];
  Mt->elements[eid].nip = new long* [nb];

  for (ii=0;ii<nb;ii++){
    Mt->elements[eid].intordsm[ii] = new long [nb];
    Mt->elements[eid].nip[ii] = new long [nb];
    for (jj=0;jj<nb;jj++){
      Mt->elements[eid].intordsm[ii][jj]=intordsm[ii][jj];
      Mt->elements[eid].nip[ii][jj]=nip[ii][jj];
    }
  }
}

/**
  This function returns correct number of dofs on the spring for given element eid and also setups 
   ndofe attribute to this value.
*/
long springel::give_ndofe (long eid)
{  
  if ((Gtm->give_ndofe(eid) == 0) && (Gtm->give_nne(eid) > 0))
  {  
    ivector enodes(1);
    Mt->give_elemnodes (eid, enodes);
    Gtm->gelements[eid].ndofe = Mt->give_ndofn(enodes[0]);
  }

  return Gtm->give_ndofe(eid);
}

/**
  This function computes stiffness %matrix of given block. The type
  of element determines direction of the spring support and thus the
  stiffness contribution is stored to the appropriate position of the
  matrix sm.

  @param eid - element id.
  @param ri  - block row id
  @param ci  - block column id
  @param sm  - stiffness %matrix where the results are stored.

*/
void springel::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  matrix d(1, 1);
  Mm->matstiff (d,Mt->elements[eid].ipp[ri+0][ci+0]);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      sm[0][0] = d[0][0];
      break;
    case spring_2:
      sm[1][1] = d[0][0];
      break;
    case spring_3:
      sm[2][2] = d[0][0];
      break;
    case spring_4:
      sm[3][3] = d[0][0];
      break;
    case spring_5:
      sm[4][4] = d[0][0];
      break;
    case spring_6:
      sm[5][5] = d[0][0];
      break;
    default:{break;}
  }
}

/**
  This function computes resulting stiffness %matrix

  @param eid - element id.
  @param sm  - stiffness %matrix where the results are stored.

*/
void springel::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

/**
  This function computes mass matrix of given element.

  @param eid - element id.
  @param mm  - mass %matrix where the results are stored
*/
void springel::mass_matrix (long eid,matrix &mm)
{
  fillm (0.0,mm);
}

/**
  This function computes strains. Not yet implemented.
  11.8.2001
*/
void springel::strains (long eid,long lcid)
{
  long n = Mt->give_ndofe(eid);
  vector  r(n);
  vector  sig(n);
  vector  eps(1);
  long ii;
  double temp;

  eldispl (0,eid,r.a);
  ii = Mt->elements[eid].ipp[0][0];
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      eps[0] = r[0];
      break;
    case spring_2:
      eps[0] = r[1];
      break;
    case spring_3:
      eps[0] = r[2];
      break;
    case spring_4:
      eps[0] = r[3];
      break;
    case spring_5:
      eps[0] = r[4];
      break;
    case spring_6:
      eps[0] = r[5];
      break;
    default:{break;}
  }
  Mm->storestrain (lcid,ii,eps);
}

/**
  This function computes stresses. Not yet implemented.

  11.8.2001
*/
void springel::stresses (long eid,long lcid)
{
  long n = Mt->give_ndofe(eid);
  ivector cn(n);
  vector  r(n);
  vector  sig(1);
  vector  eps(1);
  matrix  d(1,1);
  long ii;

  ii = Mt->elements[eid].ipp[0][0];
  Mm->givestrain(lcid,ii,eps);
  Mm->matstiff(d, ii);
  mxv(d, eps, sig);
  Mm->storestress(lcid, ii, sig);
}

/**
  This function computes internal forces of given block.

  @param lcid - load case id
  @param eid - element id
  @param ri - block row id
  @param ci - block column id
  @param ifor - vector of internal forces

  12.8.2001
*/
void springel::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long ii;
  matrix d(1,1);
  vector sig(1);

  ii = Mt->elements[eid].ipp[ri+0][ci+0];
  if (Mp->strcomp==1)
    Mm->computenlstresses (ii);
  Mm->givestress (lcid,ii, sig);
  fillv(0.0, ifor);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      ifor[0] = sig[0];
      break;
    case spring_2:
      ifor[1] = sig[0];
      break;
    case spring_3:
      ifor[2] = sig[0];
      break;
    case spring_4:
      ifor[3] = sig[0];
      break;
    case spring_5:
      ifor[4] = sig[0];
      break;
    case spring_6:
      ifor[5] = sig[0];
      break;
    default:{break;}
  }
/*  fillv(0.0, ifor);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      ifor[0] = d[0][0]*r[0];
      break;
    case spring_2:
      ifor[1] = d[0][0]*r[1];
      break;
    case spring_3:
      ifor[2] = d[0][0]*r[2];
      break;
    case spring_4:
      ifor[3] = d[0][0]*r[3];
      break;
    case spring_5:
      ifor[4] = d[0][0]*r[4];
      break;
    case spring_6:
      ifor[5] = d[0][0]*r[5];
      break;
    default:{break;}
  }*/
}

/**
  This function computes resulting internal forces.

  @param lcid - load case id
  @param eid - element id
  @param ifor - vector of internal forces

  12.8.2001
*/
void springel::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,0,0,ifor);
}

void springel::intpointval (long eid,vector &nodval,vector &ipval)
{
  ipval[0]=nodval[0];
}

