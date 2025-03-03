#include "glasgownew.h"
#include "matrix.h"
#include "vector.h"
#include "vecttens.h"
#include "global.h"
#include "intpoints.h"
#include "mechmat.h"

glasgownew::glasgownew ()
{

}

glasgownew::~glasgownew ()
{

}

void glasgownew::read (XFILE *in)
{

}

/**
   function computes material matrix of displacement-temperature block
   it is based on free thermal strains
   
   @param ipp - number of integration point
   @param m - auxiliary matrix
   
   JK, 24.10.2004
*/
void glasgownew::material_matrix_fts (long ipp,matrix &d,strastrestate ssst)
{
  long nmc;
  double t;
  
  //  number of mechanical strain/stress components
  nmc = d.m;

  matrix m(nmc,1),dd(nmc,nmc);
  
  //  auxiliary matrix
  tensor_vector_matrix (ssst,m);
  
  //  elastic stiffness matrix
  Mm->elmatstiff (dd,ipp);

  //  actual temperature in Kelvins
  t = Mm->givenonmechq(temperature, ipp);
  
  //  normalized temperature
  norm_tempr = (t-273.15-20.0)/100.0;
  
  //  coefficient of free thermal strain
  if (0.0<= norm_tempr && norm_tempr<=6.0)
    alpha = 6.0e-5/(7.0-norm_tempr);
  else
    alpha = 0.0;
  
  //  matrix m
  cmulm (alpha,m);
  
  //  material matrix of coupled model
  mxm (dd,m,d);

  return;
}

/**
   thermal damage
   
   JK, 29.10.2004
*/
void glasgownew::material_matrix_td (long ipp,matrix &d,strastrestate ssst)
{
  long nmc;
  double t_old,t_new,dt,chi,htd;
  
  //  number of mechanical strain/stress components
  nmc = d.m;
  
  matrix m(nmc,1),dd(nmc,nmc),toteps(nmc,1),epst(nmc,1),epsel(nmc,1),sig(nmc,1);
  
  //  auxiliary matrix
  tensor_vector_matrix (ssst,m);
  
  //  elastic stiffness matrix
  Mm->elmatstiff (dd,ipp);

  //  actual temperature in Kelvins
  //  t_new = 
  //  previous temperature in Kelvins
  //  t_old = 
  //  the highest reached temperature
  //  kappa =
  //  the highest reached normalized temperature
  //  t_hat = 
  //  total strains
  //Mm->givestrain (0,ipp,vector &eps);
  
  //  normalized temperature
  norm_tempr = (t_new-273.15-20.0)/100.0;
  
  //  coefficient of free thermal strain
  if (0.0<= norm_tempr && norm_tempr<=6.0)
    alpha = 6.0e-5/(7.0-norm_tempr);
  else
    alpha = 0.0;
  
  //  matrix m
  cmulm (alpha,m);
  
  //  increment of temperature
  dt = t_new-t_old;
  
  //  increment of thermal strains
  cmulm (dt,m,epst);
  
  //  elastic strains
  subm(toteps,epst,epsel);
  
  //  effective stress
  mxm(dd,epsel,sig);
  
  if (norm_tempr>t_hat){
    t_hat=norm_tempr;
  }
  
  //  thermal damage parameter
  chi = 0.2*t_hat - 0.01*t_hat*t_hat;
  
  //  damage modulus
  htd = (0.2-0.02*t_hat)/100.0;
  
  
  //  material matrix of coupled model
  mxm (dd,m,d);
  cmulm (htd,sig,sig);
  addm (d,sig,d);

  return;  
}

/**
   function computes free thermal strains
   
   @param ipp - number of integration point
   
   JK, 31.10.2004
*/
void glasgownew::free_thermal_strains (long ipp, vector &epsft)
{
  long nmc, i;
  double t,t_new;
  strastrestate ssst;
  
  //  type of strain/stress state
  ssst = Mm->ip[ipp].ssst;
  
  //  number of mechanical strain/stress components
  nmc = Mm->ip[ipp].ncompstr;

  matrix m(nmc,1);
  
  //  auxiliary matrix
  tensor_vector_matrix (ssst,m);
  
  
  //  actual temperature in Kelvins
  t_new = Mm->givenonmechq(temperature, ipp);
  
  
  //  normalized temperature
  norm_tempr = (t-273.15-20.0)/100.0;
  
  //  coefficient of free thermal strain
  if (0.0<= norm_tempr && norm_tempr<=6.0)
    alpha = 6.0e-5/(7.0-norm_tempr);
  else
    alpha = 0.0;
  
  //  matrix m
  cmulm (alpha,m);

  for (i=0; i<epsft.n; i++)
    epsft[i] = m[i][0];  

  return;
}


/**
  This function computes thermal damage parameter chi which is the result of the
  thermal damage function.

  @param ipp - integration point number
  @param tempr - actual temperature
  @param kappa - %vector of the parameters of thermal damage function
                 it contains the maximum of either the largest value attained by temperature
		 or the reference temperature

  Returns value of thermal damage.

*/
double glasgownew::thermdamfunction (long ipp,double tempr,vector &kappa)
{
  double chi = 0.0;
  double tkappa0 = 20.0;
  
  if (tempr > kappa[0])
    kappa[0] = tempr;
  if (kappa[0] > tkappa0)
    chi = 2.0e-3*(kappa[0] - tkappa0)*(1.0-5.0e-4*(kappa[0]-tkappa0));
  else
    chi = 0.0;
  
  return chi;
}



void glasgownew::nlstresses(long ipp, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  vector eps(ncomp), eps_ft(ncomp), eps_el(ncomp);
  matrix d(ncomp, ncomp);
  vector sig(ncomp);
  strastrestate ssst = Mm->ip[ipp].ssst;
  long i;
  double chi;
  vector kappa(1);

  //  total strains
  for (i=0;i<ncomp;i++)
    eps[i]=Mm->ip[ipp].strain[i];
  // stiffness matrix
//  matstiff(d, ssst);
  // computation of free thermal strains
  free_thermal_strains(ipp, eps_ft);  
  // elastic strains
  subv(eps, eps_ft, eps_el);
  // stress computation
  mxv (d, eps, sig);  
  // restoring the largest attained temperature
  kappa[0] = Mm->ip[ipp].eqother[0+ido];
  // thermal damage computation
  chi = thermdamfunction (ipp, Mm->givenonmechq(temperature, ipp), kappa);
  // stress with thermal damage influence
  cmulv(1.0-chi, sig);

  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].stress[i]=sig[i];  
  
  return;
}

void glasgownew::updateval(long ipp, long ido)
{
  Mm->ip[ipp].eqother[0+ido] = Mm->givenonmechq(temperature, ipp);
  return;
}
