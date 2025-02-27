#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splas1d.h"
#include "global.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"


splas1d::splas1d (void)
{
  fs=0.0;  k=0.0;
}

splas1d::~splas1d (void)
{

}

void splas1d::read (FILE *in)
{
  fscanf (in,"%lf %lf",&fs,&k);
  sra.read (in);
}

double splas1d::yieldfunction (matrix &sig,vector &q)
  //  function evaluates yield function for given stresses
  //  12.8.2001
{
  double f;
  f = fabs(sig[0][0])-(fs+k*q[0]);
  return f;
}

/**
   function evaluates derivatives of yield function
   with respect of stress components
   
   12.8.2001
*/
void splas1d::deryieldfsigma (matrix &sig,matrix &dfds)
{
  if (sig[0][0]<0.0)  dfds[0][0]=-1.0;
  else                dfds[0][0]=1.0;
}

/**
   function evaluates derivatives of yield function
   with respect of internal variable q
   
   21.8.2001
*/
void splas1d::deryieldfq (vector &dq)
{
  dq[0]=k*(-1.0);
}

void splas1d::plasmod (matrix &h)
  //  function assembles matrix of generalized plastic moduli
  //  28.10.2001
{
  h[0][0]=k;
}

double splas1d::plasmodscalar(vector &qtr)
{
  double ret;
  vector dfq(qtr.n), hp(qtr.n);
  matrix h(qtr.n, qtr.n);

  deryieldfq(dfq);
  plasmod (h);
  mxv (h,dfq,hp);
  scprd (hp,dfq,ret);
  return ret;
}

void splas1d::updateq(double dgamma, vector &q)
{
  long j;
  vector dfq(q.n), hp(q.n);
  matrix h(q.n, q.n);

  deryieldfq(dfq);
  plasmod (h);
  mxv (h,dfq,hp);
  for (j=0;j<q.n;j++)
    q[j]-=dgamma*hp[j];
}

void splas1d::nlstresses (long ipp, long im, long ido)
{
  long ni;
  double err;
  vector epsn(1),epsp(1),q(1);
  double gamma;

  //  initial values
  epsn[0] = Mm->ip[ipp].strain[0];
  gamma   = Mm->ip[ipp].eqother[ido+0];
  epsp[0] = Mm->ip[ipp].eqother[ido+1];
  q[0]    = Mm->ip[ipp].eqother[ido+2];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  //  new data storage
  Mm->ip[ipp].other[ido+0]=gamma;
  Mm->ip[ipp].other[ido+1]=epsp[0];
  Mm->ip[ipp].other[ido+2]=q[0];
}

void splas1d::nonloc_nlstresses (long ipp, long im, long ido)
{
  long ni;
  double err;
  vector epsn(1),epsp(1),q(1);
  double gamma;

  //  initial values
  epsn[0] = Mm->ip[ipp].strain[0];
  gamma   = Mm->ip[ipp].eqother[ido+0];
  epsp[0] = Mm->ip[ipp].nonloc[ido+1];
  q[0]    = Mm->ip[ipp].eqother[ido+2];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  //  new data storage
  Mm->ip[ipp].other[ido+0]=gamma;
  Mm->ip[ipp].other[ido+1]=epsp[0];
  Mm->ip[ipp].other[ido+2]=q[0];
}

void splas1d::matstiff (matrix &d,long ipp,long ido)
{
  double denom1,denom2;
  vector sig(1),dq(1),u(1);
  matrix de(1,1),h(1,1),g(1,1),sigt(3,3);

  if (Mp->nlman->stmat==0){
    Mm->elmatstiff (d,ipp);
  }
  else{

    //  elastic stiffness matrix
    Mm->elmatstiff (de,ipp);

    //  stress component
    sig[0]=Mm->ip[ipp].stress[0];
    if (Mm->ip[ipp].eqother[ido+1]<Mp->zero){
      Mm->elmatstiff (d,ipp);
    }
    else{


      deryieldfsigma (sigt,sigt);
      tensor_vector (sig,sigt,Mm->ip[ipp].ssst,stress);

      vxv (sig,sig,h);
      mxm (de,h,g);
      mxm (g,de,h);

      mxv (de,sig,u);
      scprd (u,sig,denom1);

      deryieldfq (dq);
      scprd (dq,dq,denom2);

      cmulm (1.0/(denom1+denom2),h);

      subm (de,h,d);
    }
  }
}

void splas1d::updateval (long ipp, long ido)
{
  Mm->ip[ipp].eqother[ido+0]=Mm->ip[ipp].other[ido+0];
  Mm->ip[ipp].eqother[ido+1]=Mm->ip[ipp].other[ido+1];
  Mm->ip[ipp].eqother[ido+2]=Mm->ip[ipp].other[ido+2];
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  Returns vector of irreversible strains via parameter epsp
*/
void splas1d::giveirrstrains (long ipp, long ido, vector &epsp)
{
  epsp[0] = Mm->ip[ipp].eqother[ido+1];
}



void splas1d::changeparam (atsel &atm,vector &val)
{
  long i;

  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      fs=val[i];
      break;
    }
    case 1:{
      k=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}

double splas1d::give_consparam (long ipp, long ido)
{
  double gamma;
  gamma = Mm->ip[ipp].other[ido+0];
  return gamma;
}

long splas1d::give_num_interparam ()
{
  return 1;
}

void splas1d::give_interparam (long ipp,long ido,vector &q)
{
  long ncompstr=Mm->ip[ipp].ncompstr;
  
  q[0]=Mm->ip[ipp].eqother[ido+ncompstr+1];
}
