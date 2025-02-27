#include "splas1d.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include <math.h>

splas1d::splas1d (void)
{
  //  flow stress
  fs=0.0;
  //  plastic modulus
  k=0.0;
}

splas1d::~splas1d (void)
{

}

/**
   function reads input material parameters
   
   @param in - pointer to input file
   
   JK
*/
void splas1d::read (XFILE *in)
{
  //  yield stress
  //  plastic modulus
  xfscanf (in,"%k%lf %k%lf","fs", &fs, "k", &k);
  //  type of stress return algorithm
  sra.read (in);
  //  type of hardening/softening
  hs.read (in);
}

/**
   The function prints input material parameters to the opened text file.
   
   @param out - pointer to the opened text file for output.
   
   TKo
*/
void splas1d::print(FILE *out)
{
  //  yield stress
  //  plastic modulus
  fprintf(out, "%le %le ", fs, k);
  //  type of stress return algorithm
  sra.print(out);
  fprintf(out, " ");
  //  type of hardening/softening
  hs.print(out);
}

/**
   function evaluates the yield function
   
   engineering components of stress are stored in the %vector
   with 6 components
   the same %vector is used in 2 and 3 D problems
   
   @param sig - stresses
   @param q - internal (hardening/softening) variables
   
   JK, 12.8.2001
*/
double splas1d::yieldfunction (vector &sig,vector &q)
{
  double f;
  f = fabs(sig[0])-(fs+q[0]);
  return f;
}

/**
   function evaluates derivatives of the yield function
   with respect to stress components
   
   engineering components of stress are stored in the %vector
   with 6 components
   the same %vector is used in 2 and 3 D problems
   
   @param sig - stresses
   @param dfds - %vector containing derivatives with respect to stresses
   
   JK, 12.8.2001
*/
void splas1d::dfdsigma (vector &sig,vector &dfds)
{
  if (sig[0]<0.0)  dfds[0] = -1.0;
  else             dfds[0] =  1.0;
}

/**
   function computes second derivatives of the yield function with respect to stresses
   all of them are always equal to zero
   
   JK, 8.8.2005
*/
void splas1d::dfdsigmadsigma(matrix &dfdsds)
{
  fillm (0.0,dfdsds);
}

/**
   function evaluates derivatives of the yield function with respect to internal variable q
   
   @param dq - %vector containing derivatives of yield function with respect to internal variable
   
   JK, 21.8.2001
*/
void splas1d::dfdqpar (vector &dq)
{
  dq[0]=-1.0;
}


/**
   function evaluates derivates of the yield function with respect to stresses and
   internal variables (hardening/softening parameters)
   
   @param dfdsdq - derivatives of yield function with respect to stresses and internal variables
   
   JK, 8.8.2005
*/
void splas1d::dfdsigmadq (matrix &dfdsdq)
{
  fillm (0.0,dfdsdq);
}

/**
   function evaluates derivates of the yield function with respect to stresses and
   internal variables (hardening/softening parameters)
   
   @param dfdsdq - derivatives of yield function with respect to stresses and internal variables
   
   JK, 8.8.2005
*/
void splas1d::dfdqpardqpar (matrix &dfdqdq)
{
  fillm (0.0,dfdqdq);
}



/**
   function computes derivatives of hardening function with respect to consistency parameter
   
   @param ipp - id of integration point
   @param idpm - id of plasticity model
   @param dhdc - %vector containing derivatives of hardening function with respect to consistency parameter
   
   JK, 8.8.2005
*/
/*
void splas1d::dhdgamma (long ipp,vector &epsp, vector &sig,vector &dhdc)
{
  matrix dgds(3,3);
  
  //  derivatives of yield function with respect to stresses
  //dfdsigma (sig,dgds);
  
  //  derivatives of hardening function with respect to consistency parameter
  //hs.dhdgamma (ipp,notdef,epsp,dgds,dhdc);
  
  dhdc[0] *= -1.0;
}
*/

/**
   function assembles the %vector of hardening h which governs
   the evolution of internal parameters
   
   17. 6. 2015
*/
void splas1d::hardvect (vector &hv)
{
  hv[0]=0.0-k;
}




void splas1d::plasmod (matrix &h)
{
  h[0][0]=0.0-k;
}


double splas1d::plasmodscalar(vector &/*sig*/,vector &/*epsp*/,vector &/*qtr*/,double /*gamma*/)
{
  /*
  double ret;
  vector dfq(qtr.n),res(1);
  matrix dfds(1,1),h(qtr.n, qtr.n),hp(1,1);

  deryieldfsigma (sig,dfds);
  deryieldfq(dfq);
  plasmod (h);

  mxm (h,dfds,hp);
  vxm (dfq,hp,res);
  
  ret=-1.0*res[0];
  return ret;
  */
  /*
  //  jine zpevneni
  double ret,a,b,c;
  vector dfq(qtr.n),res(1);
  matrix dfds(1,1),h(qtr.n, qtr.n),hp(1,1);

  deryieldfsigma (sig,dfds);
  deryieldfq(dfq);
  plasmod (h);
  
  a=dfds[0][0];
  b=dfq[0];
  c=h[0][0];
  
  ret = b*2.0*c*a*a*gamma;
  ret=-1.0*ret;
  return ret;
  */


  //  a jeste jine zpevneni
  /*
  double ret,a,b,c,e;
  vector dfq(qtr.n),res(1);
  matrix dfds(1,1),h(qtr.n, qtr.n),hp(1,1);

  //dfdsigma (sig,dfds);
  dfdq(dfq);
  plasmod (h);
  
  a=dfds[0][0];
  b=dfq[0];
  c=h[0][0];
  e=epsp[0];
  
  if (e<1.0e-6)
    e=1.0e-6;
  
  ret = b*0.5*c*a/sqrt(e);
  
  ret=-1.0*ret;
  */
  
  vector dfq(1),hv(1);
  
  hardvect (hv);
  dfdqpar (dfq);

  return (ss (dfq.a,hv.a,hv.n));
}


void splas1d::updateq(long /*ipp*/,double dgamma, vector &/*epsp*/,vector &q)
{
  vector hv(1);
  
  hardvect (hv);
  
  q[0]-=dgamma*hv[0];
}

void splas1d::nlstresses (long ipp, long im, long ido)
{
  long ni;
  double err;
  vector epsn(1),epsp(1),q(1);
  double gamma;

  //  initial values
  epsn[0] = Mm->ip[ipp].strain[0];
  epsp[0] = Mm->ip[ipp].eqother[ido+0];
  gamma   = Mm->ip[ipp].eqother[ido+1];
  q[0]    = Mm->ip[ipp].eqother[ido+2];

  //  stress return algorithm
  switch (sra.tsra){
  case cp:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  case gsra:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  default:{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  }
  
  //  new data storage
  Mm->ip[ipp].other[ido+0]=epsp[0];
  Mm->ip[ipp].other[ido+1]=gamma;
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
  gamma   = Mm->ip[ipp].eqother[ido+1];
  epsp[0] = Mm->ip[ipp].nonloc[ido+0];
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
  Mm->ip[ipp].other[ido+0]=epsp[0];
  Mm->ip[ipp].other[ido+1]=gamma;
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


      //dfdsigma (sigt,sigt);
      tensor_vector (sig,sigt,Mm->ip[ipp].ssst,stress);

      vxv (sig,sig,h);
      mxm (de,h,g);
      mxm (g,de,h);

      mxv (de,sig,u);
      scprd (u,sig,denom1);

      dfdqpar (dq);
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
  epsp[0] = Mm->ip[ipp].eqother[ido+0];
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
      //k=val[i];
      break;
    }
    default:{
    print_err("wrong attribute number",__FILE__,__LINE__,__func__);
    }
    }
  }
}

double splas1d::give_consparam (long ipp, long ido)
{
  double gamma;
  gamma = Mm->ip[ipp].other[ido+1];
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
