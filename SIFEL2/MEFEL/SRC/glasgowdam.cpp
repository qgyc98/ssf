#include "glasgowdam.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>

#define nijac 20
#define limit 1.0e-8



/**
  This constructor inializes attributes to zero values.
  
  TKo
*/
glasgowdam::glasgowdam (void)
{
  st=0.0;
  k=0.0;
  gf0=0.0;
  lc=0.0;
  reftemp = 0.0;
  ft = paramf_type(0);
}



/**
  This destructor is only for the formal purposes.
  
  TKo
*/
glasgowdam::~glasgowdam (void)
{

}



/**
  This function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  TKo
*/
void glasgowdam::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %m", &st, &gf0, &lc, &paramf_type_kwdset, (int *)(&ft));
  xfscanf(in, "%lf %lf", &k, &reftemp);
}



/**
  This function computes parameters for the damage function
  Different types of the parameter computing are given by the
  attribute ft.

  @param  ipp - integration point number
  @param  eps - %vector of the strains
  @param  kappa - %vector where the computed parameters are stored
  
  TKo
*/
void glasgowdam::damfuncpar(long ipp, vector &eps, vector &kappa)
{
  vector poseps;
  vector princeps;
  matrix pvect, epsm;
  vector epst;
  vector dev;
  double i1e, j2e, nu, e, tmp, a, b;
  long ncomp, idem, im;

  switch (ft)
  {
    case vonmises :
      // glasgow modified von Mises norm
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(6, epst));
      reallocv(RSTCKVEC(6, dev));
      give_full_vector(epst, eps, Mm->ip[ipp].ssst);
      i1e = first_invar(epst);
      deviator (epst,dev);
      j2e = j2_strain_invar (dev);
      idem = Mm->ip[ipp].gemid();
      if (Mm->ip[ipp].tm[idem] != elisomat)
      {
        fprintf(stderr, "\n\n Invalid type of elastic material is required");
        fprintf(stderr, "\n  in function glasgowdam::damfuncpar, (file %s, line %d)\n", __FILE__, __LINE__);
      }
      im = Mm->ip[ipp].idm[idem];
      e  = Mm->eliso[im].e;
      nu = Mm->eliso[im].nu;
      a   = (k-1.0)/(2.0*k)/(1.0-2.0*nu);
      b   = 3.0/k/(1+nu)/(1.0+nu);
      tmp = a*a*i1e*i1e + b*j2e;
      kappa[0] = a*i1e+sqrt(tmp);

      /*
      // original von Mises norm
      kappa[0]  = (k-1.0)/(2.0*k)/(1.0-2.0*nu)*i1e;
      tmp = (k-1.0)*(k-1.0)/(1.0-2.0*nu)/(1.0-2.0*nu)*i1e*i1e -12.0*k/(1.0+nu)/(1.0+nu)*j2e;
      kappa[0] += 1.0/(2.0*k)*sqrt(tmp);*/
      break;
    case normazar :
      reallocm(RSTCKMAT(3, 3, epsm));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, princeps));
      reallocv(RSTCKVEC(3,poseps));
      fillv(0.0, princeps);
      vector_tensor(eps, epsm, Mm->ip[ipp].ssst, strain);
      princ_val(epsm, princeps, pvect, nijac, limit, Mp->zero, 3, 1);
      extractposv(princeps, poseps);
      scprd(poseps, poseps, tmp);
      kappa[0] = sqrt(tmp);
      break;
    default :
    {
      fprintf(stderr, "\n\n Unknown damage parameter function type");
      fprintf(stderr, "\n  in function glasgowdam::damfuncpar, (file %s, line %d)\n", __FILE__, __LINE__);
    }
  }
  return;
}



/**
  This function computes damage parameter omega which is the result of the
  damage function.

  @param ipp - integration point number
  @param tempr - actual temperature
  @param chi - thermal damage parameter
  @param kappa - %vector of the parameters of damage function
                 kappa[0] = maximum reached equivalent strain
                 kappa[1] = maximum reached temperature in Celsius
  
  TKo
*/
double glasgowdam::damfunction(long ipp,double tempr,double chi,vector &kappa)
{
  double omega = 0.0, e0, kappa0, tn, gamma, gf,th;
  long idem = Mm->ip[ipp].gemid();
  
  if (Mm->ip[ipp].tm[idem] != elisomat)
  {
    fprintf(stderr, "\n\nError - ip %ld has wrong type of elastic material combined with", ipp);
    fprintf(stderr, "\n scalar damge. The elastic isotropic material is requried");
    fprintf(stderr, "\n Function glasgowdam::damfunction (file %s, line %d)\n", __FILE__, __LINE__);
    return 0.0;
  }

  
  // Young modulus
  e0 = Mm->eliso[Mm->ip[ipp].idm[idem]].e;
  
  //  normalized temperature
  tn = (tempr-20.0)/100.0;
  //  maximum reached normalized temperature
  th = (kappa[1]-20.0)/100.0;
  if (tn>th){
    th=tn;
  }
  //  check of maximum reached normalized temperature
  if ((th < 0.0) || (th > 7.9)){
    fprintf(stderr, "\n\nWarning - largest value of normalised temperature\n");
    fprintf(stderr, "          is out of suggested range 0.0-7.9 (ipp %ld)\n",ipp);
  }
  
  //  threshold of the mechanical damage (depends on temperature)
  //  kappa0 = st/e0*(1.0 - 0.016*th*th)/(1.0 - 0.1*th)/(1.0 - 0.1*th);
  kappa0 = st/e0*(1.0 - 0.016*th*th);
  
  //  fracture energy
  gf = gf0*(1.0+0.39*th-0.07*th*th);
  
  //  softening parameter
  gamma = (1.0-chi)*st*lc/gf;
  
  if (kappa[0] < kappa0)
    // elastic loading
    return 0.0;
  
  //  mechanical damage parameter
  omega = 1.0-(kappa0/kappa[0])*exp(-gamma*(kappa[0]-kappa0));

  return omega;
}



/**
   function computes 
   
   @param t - temeprature
   @param dt - incrment of temperature
   @param eps - thermal strain increment
   
*/
void glasgowdam::compute_thermdilat (double t,double dt,vector &eps)
{
  double auxt,alpha;
  
  auxt = (t-20.0)/100.0;
  if (auxt>=0.0 && auxt<=6.0){
    alpha = 6.0e-5/(7.0-auxt);
  }
  else{  alpha=0.0; }
  
  nullv(eps);
  
  eps(0)=alpha*dt;
  eps(1)=alpha*dt;
  eps(2)=alpha*dt;
  
}



/**
  This function computes thermal damage parameter chi which is the result of the
  thermal damage function.

  @param ipp - integration point number
  @param tempr - actual temperature
  @param kappa - %vector of the parameters of thermal damage function
                 it contains the maximum of either the largest value attained by temperature
		 or the reference temperature

*/
double glasgowdam::thermdamfunction (long /*ipp*/,double tempr,vector &kappa)
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



/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  TKo
*/
void glasgowdam::matstiff (matrix &d,long ipp,long ido)
{
  double omega, chi, dp;
  vector kappa(1);

  switch (Mp->nlman->stmat)
  {
    case initial_stiff :
      Mm->elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff :
      Mm->elmatstiff (d,ipp);
      omega =Mm->ip[ipp].eqother[ido+1];
      chi = 0.0;
      kappa[0] = 20.0;
      chi   = thermdamfunction(ipp, reftemp, kappa);
      dp = (1.0 - chi)*(1.0-omega);
      if (dp < 0.000001)
        dp = 0.000001;
      cmulm (dp,d);
      break;
    default :
      fprintf(stderr, "\n\nError - unknown type of stifness matrix");
      fprintf(stderr, "\n in function glasgowdam::matstiff (%s, line %d)\n", __FILE__, __LINE__);
  }
}



/**
  This function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  TKo, 7.10.2001
*/
void glasgowdam::nlstresses (long ipp, long /*im*/, long ido)
{
  long i,ncomp=Mm->ip[ipp].ncompstr;
  vector epsn(ASTCKVEC(ncomp)), sigma(ASTCKVEC(ncomp)), kappa(ASTCKVEC(2));
  vector updkappa(ASTCKVEC(2)), epsft(ASTCKVEC(ncomp)), eps(ASTCKVEC(ncomp)), tkappa(ASTCKVEC(1));
  matrix d(ASTCKMAT(ncomp, ncomp));
  vector epst(ASTCKVEC(6));
  double omega, chi, dp;

  //  initial values
  for (i=0;i<ncomp;i++){
    epsn[i] = Mm->ip[ipp].strain[i];
  }
  kappa[0] = Mm->ip[ipp].eqother[ido+0];
  updkappa[1] = kappa[1] = reftemp;

  // damage stress solver
  if (Mp->matmodel == local)
  // local model is used
  {
    damfuncpar(ipp, epsn, updkappa);
    if (updkappa[0] > kappa[0])
    {
      kappa[0] = updkappa[0];
      omega = damfunction(ipp, reftemp, 0.0, kappa);
    }
    else
      omega = Mm->ip[ipp].eqother[ido+1];

    copyv(epsn, eps);

    //  free thermal strains
    compute_thermdilat (reftemp, reftemp-20.0, epst);
    //  conversion of tensor notation into vector one
    give_red_vector(epst, epsft, Mm->ip[ipp].ssst);
    subv(epsn, epsft, eps);

    tkappa[0] = 20.0;
    // thermal damage  
    chi = thermdamfunction(ipp, reftemp, tkappa);

//    chi = 0.0;
    dp = (1.0 - omega)*(1.0-chi);
    Mm->elmatstiff (d,ipp);
    mxv(d, eps, sigma);
    cmulv(dp, sigma);
  }
  else{
    print_err("the non-local glasgowdam model has not yet been implemented", __FILE__, __LINE__, __func__);
    abort();
  }

  //  new data storage
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].stress[i]=sigma[i];
  }
  Mm->ip[ipp].other[ido+0]=kappa[0];
  Mm->ip[ipp].other[ido+1]=omega;
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  TKo
*/
void glasgowdam::updateval (long ipp,long im,long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}
