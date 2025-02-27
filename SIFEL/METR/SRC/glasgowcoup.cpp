#include "glasgowcoup.h"
#include "vecttens.h"
#include "global.h"
#include "mechmat.h"

glasgowcoup::glasgowcoup ()
{

}

glasgowcoup::~glasgowcoup ()
{

}

void glasgowcoup::read (XFILE */*in*/)
{

}

/**
   function computes material matrix of displacement-temperature block
   it is based on free thermal strains
   
   @param ipp - number of integration point
   @param m - auxiliary matrix
   
   JK, 24.10.2004
*/
void glasgowcoup::material_matrix_fts (long ipp,matrix &d,strastrestate ssst)
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
  t = 0.0; // provisionally
  
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
}

/**
   thermal damage
   
   JK, 29.10.2004
*/
void glasgowcoup::material_matrix_td (long ipp,matrix &d,strastrestate ssst)
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
  t_new = 0.0; // provisionally
  //  previous temperature in Kelvins
  t_old =  0.0; // provisionally
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
  
}

