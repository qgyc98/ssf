#include "scaldamcc.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "mathem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nijac 20
#define limit 1.0e-8



/**
  This constructor inializes attributes to zero values.
*/
scaldamcc::scaldamcc(void)
{
  // damage function type
  damfunct = simpexp;
  // simple exponential paparemetrs
  ft = uf = 0.0;
  // mazars exponential function parameters
  k0 = at = bt = ac = bc = beta = k = 0.0;
  // equivalent strain norm type
  eqepsnt = normazar;
}



/**
  This destructor is only for the formal purposes.
*/
scaldamcc::~scaldamcc(void)
{

}



/**
  This function reads material parameters from the opened text file given
  by the argument in.

  @param in[in] - pointer to the opened text input file.

  @return The function does not return anything.

  Created by Tomas Koudelka, 06.2003
  Modified by Tomas Koudelka, 06.2019
*/
void scaldamcc::read(XFILE *in)
{
  xfscanf (in, "%k%m", "damfunct", &damfunc_type_kwdset, (int*)(&damfunct));
  switch(damfunct){
    case simpexp:
      xfscanf(in, "%k%lf %k%lf", "ft", &ft, "uf", &uf);
      xfscanf(in, "%k%lf %k%lf %k%lf", "ac", &ac, "bc", &bc, "beta", &beta);
      break;
    case mazarsexp:
      xfscanf (in,"%k%lf %k%lf %k%lf %k%lf %k%lf %k%lf", 
               "k0", &k0, "at", &at, "bt", &bt, "ac", &ac, "bc", &bc, "beta", &beta); 
      break;
    default:
      print_err("unknown type %d of damage function is required", __FILE__, __LINE__, __func__, int(damfunct));
      abort();
  }
  xfscanf (in,"%k%m", "norm_type", &paramf_type_kwdset, (int *)(&eqepsnt));
  if (eqepsnt == vonmises)
    xfscanf(in, "%k%lf", "k", &k);
}



/**
  This function prints material parameters from the opened text file given
  by the argument out.

  @param out[in] - pointer to the opened text output file.

  @return The function does not return anything.

  Created by Tomas Koudelka, 06.2019
*/
void scaldamcc::print(FILE *out)
{
  fprintf (out, "%d", int(damfunct));
  switch(damfunct){
    case simpexp:
      fprintf(out, " %le %le", ft, uf);
      fprintf(out, " %le %le %le", ac, bc, beta); 
      break;
    case mazarsexp:
      fprintf(out," %le %le %le %le %le %le", k0, at, bt, ac, bc, beta);
      break;
    default:
      print_err("unknown type %d of damage function is required", __FILE__, __LINE__, __func__, int(damfunct));
      abort();
  }
  fprintf(out," %d", int(eqepsnt));
  if (ft == vonmises)
    fprintf(out, " %le", k);
}



/**
  This function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d[out]  - allocated matrix structure for the resulting elastic stiffness %matrix,
  @param ipp[in] - integration point number.

  @return The function returns computed %matrix in the argument d.

  Created by Tomas Koudelka,  08.2008
  Modified by Tomas Koudelka, 06.2019
*/
void scaldamcc::elmatstiff(matrix &d, long ipp)
{
  double e, e0;
  long idem  = Mm->ip[ipp].gemid();             // index of elastic material
  long idoem = Mm->givencompeqother(ipp, 0);    // total number of internal variables
  long tmp   = Mm->givencompeqother(ipp, idem); // number of internal variables for elastic material and thermal material 
  idoem -= tmp;
  
  Mm->elmatstiff(d, ipp);
  e = Mm->give_actual_ym(ipp); // actual Young modulus (it may be computed also in creep model)
  e0 = Mm->give_actual_ym(ipp, idem, idoem); // Young modulus from elastic material
  cmulm(e/e0, d);
}



/**
  This function computes parameter for the damage function, i.e. equivalent strain.
  Type of the evolution function for the equivalent strain computing is given by 
  the attribute eqepsnt.

  @param  ipp[in]   - integration point number,
  @param  eps[in]   - %vector of strains,
  @param  kappa[in] - %vector for the storage of computed equivalent strain.

  @return The function returns actual value of the equivalent strain.

  Created by Tomas Koudelka,  06.2003
  Modified by Tomas Koudelka, 06.2019
*/
double scaldamcc::damfuncpar(long ipp, vector &eps)
{
  vector poseps;
  vector peps;
  matrix pvect;
  matrix epst(ASTCKMAT(3, 3));
  vector epstt;
  double i1e, j2e, nu, e, tmp, kappa;
  long ncomp, idem;

  switch (eqepsnt)
  {
    case vonmises:
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(6, epstt));
      give_full_vector(epstt, eps, Mm->ip[ipp].ssst);
      i1e = first_invar(epstt);
      j2e = j2_strain_invar(epstt);
      idem = Mm->ip[ipp].gemid();
      if (Mm->ip[ipp].tm[idem] != elisomat){
        print_err("invalid type of elastic material %d is required", __FILE__, __LINE__, __func__, int(Mm->ip[ipp].tm[idem]));
        abort();
      }
      e  = Mm->eliso[idem].e;
      nu = Mm->eliso[idem].nu;
      kappa = (k - 1.0)/((2.0*k)*(1.0-2.0*nu))*i1e;
      tmp = sqr(k - 1.0)/sqr(1.0 - 2.0*nu)*sqr(i1e) - 12.0*k/sqr(1.0+nu)*j2e;
      kappa += 1.0/(2.0*k)*sqrt(tmp);
      break;
    case normazar:
      // kappa = \sqrt{\sum_{\alpha=1}^3 \langle\peps_{\alpha}\rangle^2}
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, peps));
      reallocv(RSTCKVEC(3, poseps));
      vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
      princ_val(epst, peps, pvect, nijac, limit, Mp->zero, 3, 1);
      extractposv(peps, poseps);
      scprd(poseps, poseps, tmp);
      kappa = sqrt(tmp);
      break;
    default:
      print_err("unknown damage parameter function type %d is required", __FILE__, __LINE__, __func__, int(eqepsnt));
      abort();
  }

  return kappa;
}



/**
  This function computes damage parameter omega which is the result of the
  damage function.

  @param ipp[in]    - integration point number,
  @param kappa[in]  - equivalent strain,
  @param eps[in]    - %vector of strains,
  @param omegao[in] - attained damage parameters at the last equilibrium state,
  @param dt[in,out] - damage associated with tension,
  @param dc[in,out] - damage associated with compression,
  @param alphat[in,out] - weigth of damage parameter dt in the resulting damage parameter omega,
  @param alphac[in,out] - weigth of damage parameter dc in the resulting damage parameter omega.

  @return The function returns actual value of the damage parameter omega.

  Created by Tomas Koudelka, 06.2019
*/
double scaldamcc::damfunction(long ipp, double kappa, vector &eps,  double omegao, double &dt, double &dc, double &alphat, double &alphac)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  matrix epst, sigt, pvect, d(ASTCKMAT(ncomp, ncomp)), c(ASTCKMAT(6, 6)), cn;
  vector peps;
  vector poseps;
  vector tenseps;
  vector compeps;
  vector sigef(ncomp), psigef, tenssig, compsig;
  vector negeps;
  double omega = 0.0, eps0;
  strastrestate ssst = Mm->ip[ipp].ssst;

  eps0 = epsefunction(ipp);
  if (kappa <= eps0)
    // elastic loading/unloading
    return omegao;

  reallocm(RSTCKMAT(3, 3, epst));
  reallocm(RSTCKMAT(3, 3, sigt));
  reallocm(RSTCKMAT(3, 3, cn));
  reallocv(RSTCKVEC(3, poseps));
  reallocv(RSTCKVEC(3, negeps));
  reallocv(RSTCKVEC(3, tenseps));
  reallocv(RSTCKVEC(3, compeps));
  reallocv(RSTCKVEC(3, peps));
  reallocv(RSTCKVEC(3, psigef));
  reallocv(RSTCKVEC(3, tenssig));
  reallocv(RSTCKVEC(3, compsig));
  reallocm(RSTCKMAT(3, 3, pvect));

  switch(damfunct)
  {
    case simpexp:
      dt = 1.0-(eps0/kappa)*exp(-(kappa-eps0)/(uf-eps0));
      if (bc*(kappa-eps0) < 700)
        dc = 1.0 - k0 * (1.0 - ac) / kappa - ac*exp(-bc*(kappa-eps0));
      else
        dc = 1.0;
      break;
    case mazarsexp:
      if (bt*(kappa-eps0) < 700)
        dt = 1.0 - k0 * (1.0 - at) / kappa - at*exp(-bt*(kappa-eps0));
      else
        dt = 1.0;
      if (bc*(kappa-k0) < 700)
        dc = 1.0 - k0 * (1.0 - ac) / kappa - ac*exp(-bc*(kappa-k0));
      else
        dc = 1.0;
      break;
    default:
      print_err("unknown damage function type %d is required", __FILE__, __LINE__, __func__, int(damfunct));
      abort();
  }
  // compute principal strains peps
  vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
  princ_val(epst, peps, pvect, nijac, limit, Mp->zero, 3, 1);

  // compute effective stresses
  Mm->elmatstiff (d,ipp);
  mxv(d,eps,sigef);
  // compute principal effective stresses
  vector_tensor(sigef, sigt, Mm->ip[ipp].ssst, stress);
  princ_val(sigt, psigef, pvect, nijac, limit, Mp->zero, 3, 1);

  extractposv(psigef,tenssig);
  //extractnegv(psigef,compsig);
  // c is the full compliance matrix
  Mm->ip[ipp].ssst = spacestress;
  Mm->elmatcompl(c, ipp);
  Mm->ip[ipp].ssst = ssst;

  // extract block for normal components of the compliance matrix => cn
  extractm(cn, c, 0, 3);
  // tenseps = tensile strains due to positive principal stresses
  mxv(cn, tenssig, tenseps);
  // compeps = compressive strains due to negative principal stresses
  mxv(cn, compsig, compeps);

  // peps = peps^+ + pesp^-
  extractposv (peps, poseps);  // poseps = positive principal strains due positive stresses
  subv(peps, poseps, negeps);  // negeps = negative principal strains (complementary to poseps)

  scprd(tenseps, poseps, alphat);
  alphat = alphat / sqr(kappa);

  //alphat = 1.0;

  //  scprd(compeps, poseps, alphac);
  //  alphac = alphac / (kappa*kappa);
  alphac = 1.0 - alphat;
  alphat = pow(alphat, beta);
  alphac = pow(alphac, beta);
  

  omega = alphat*dt + alphac*dc;

  if ((omega > 1.0) || (omega < 0.0)) {
    //    print_warning("Error - omega=%g is out of range <0,1>, element %ld, ip=%ld ", __FILE__, __LINE__, __func__, omega, Mm->elip[ipp]+1, ipp);
  }

  if (omega < omegao)
    omega = omegao;

  if(omega > 1.0)
    omega = 1.0 ;

  if(omega < 0.0)
    omega = 0.0;

  return omega;
}



/**
  This function computes material stiffnes %matrix.

  @param d[out]  - allocated matrix structure for material stiffness %matrix
  @param ipp[in] - integration point number
  @param ido[in] - index of internal variables for given material in the ipp other array

  @return The function returns actual stiffness matrix in the argument d.

  Created by Tomas Koudelka,  06.2003
  Modified by Tomas Koudelka, 06.2019
*/
void scaldamcc::matstiff(matrix &d, long ipp, long ido)
{
  double dp;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      Mm->elmatstiff (d,ipp);
      break;
    case tangent_stiff:
    case secant_stiff:
    case incr_tangent_stiff:
    case ijth_tangent_stiff:
      Mm->elmatstiff (d,ipp);
      dp=Mm->ip[ipp].other[ido+1];
      if (dp > 0.999)
        dp = 0.999;
      cmulm (1.0-dp,d);
      break;
    default:
      print_err("unknown type of stifness matrix %d is required", __FILE__, __LINE__, __func__, int(Mp->nlman->stmat));
      abort();
  }
}



/**
  This function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer,
  @param im  - index of material type for given ip,
  @param ido - index of internal variables for given material in the ipp other array.

  Created by Tomas Koudelka,  06.2003
  Modified by Tomas Koudelka, 06.2019
*/
void scaldamcc::nlstresses(long ipp, long im, long ido)
{
  long i,ncomp=Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(ncomp)), sigma(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp, ncomp));
  double nu, epse, omega, omegao;
  double kappa, updkappa;
  double dt, dc, alphat, alphac;

  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }
  //  attained values of local strains
  for (i=0;i<ncomp;i++)
    eps(i) = Mm->ip[ipp].strain[i];

  // attained values of internal variables from the previous equilibrium state
  kappa  = Mm->ip[ipp].eqother[ido+0]; // maximum equivalent strain norm in the history
  omegao = Mm->ip[ipp].eqother[ido+1]; // resulting damage parameter
  dt     = Mm->ip[ipp].eqother[ido+2]; // damage due to tension
  dc     = Mm->ip[ipp].eqother[ido+3]; // damage due to compression
  alphat = Mm->ip[ipp].eqother[ido+4]; // weight coefficient of dt
  alphac = Mm->ip[ipp].eqother[ido+5]; // weight coefficient of dc

  // elastic stiffness matrix
  Mm->elmatstiff(d, ipp);

  if ((Mp->matmodel == local) || ((Mm->ip[ipp].hmt & 2) == 0)){
    //
    // local model is used
    //
    // update equivalent strain
    updkappa = damfuncpar(ipp, eps);
    if (updkappa > kappa){ 
      // new maximum value of equivalent strain kappa was attained
      kappa = updkappa;
      // compute strain of damage treshold
      epse = Mm->epsefunction(ipp, im) ;
      if (kappa > epse) 
        // equivalent strain is out of elastic treshold
        omega = damfunction(ipp, kappa, eps, omegao, dt, dc, alphat, alphac);
      else
        // elastic loading
        omega = 0.0;
      if (omega < omegao)
        omega = omegao;
    }
    else
      // elastic loading/unloading
      omega = omegao;

    // compute stresses (1-\omega) D_e \eps
    mxv(d, eps, sigma);
    cmulv(1.0-omega, sigma);
  }

  if (Mp->matmodel == nonlocal){
    //
    // nonlocal model is used
    //
    if (Mp->nonlocphase == 1){
      // compute local values for averaging in nonlocphase=1 -> no stresses are required
      // in this case, nothing is needed to be computed in the first phase,
      // just strain tensor will be averaged
      /*kappa = damfuncpar(ipp, eps);
      Mm->ip[ipp].other[ido+0]=kappa;   // maximum equivalent strain norm in the history
      nullv(sigma);*/
      return;
    }
    else{
      // nonlocphase = 2 -> compute stresses from the averaged nonlocal values
      // collect averaged strains
      vector epsn(ASTCKVEC(ncomp));
      for (i=0; i<ncomp; i++)
        epsn(i) = Mm->ip[ipp].nonloc[i];
      // update equivalent strain
      updkappa = damfuncpar(ipp, epsn);
      //updkappa = Mm->ip[ipp].nonloc[0];      
      if (updkappa > kappa){
        // new maximum value of equivalent strain kappa was attained
        kappa = updkappa;
        // compute strain of damage treshold
        epse = Mm->epsefunction(ipp, im);
        if (kappa > epse)
          // equivalent strain is out of elastic treshold
          omega = damfunction(ipp, kappa, eps, omegao, dt, dc, alphat, alphac);
        else
          // elastic loading
          omega = 0.0;
        if (omega < omegao)
          omega = omegao;
      }
      else
        // elastic loading/unloading
        omega = omegao;

      // compute stresses from local strains (1-\omega) D_e \eps
      mxv(d, eps, sigma);
      cmulv(1.0-omega, sigma);
    }
  }

  //
  // storage of actual values
  //
  for (i=0;i<ncomp;i++) 
    Mm->ip[ipp].stress[i]=sigma(i); // actual stresses

  Mm->ip[ipp].other[ido+0]=kappa;   // maximum equivalent strain norm in the history
  Mm->ip[ipp].other[ido+1]=omega;   // resulting damage parameter
  Mm->ip[ipp].other[ido+2]=dt;      // damage due to tension
  Mm->ip[ipp].other[ido+3]=dc;      // damage due to compression
  Mm->ip[ipp].other[ido+4]=alphat;  // weight coefficient of dt
  Mm->ip[ipp].other[ido+5]=alphac;  // weight coefficient of dc
}



/**
  Function updates values in the other array reached in the previous equlibrium state to
  ones attained in the actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array,
  @param im  - index of material type for given ip,
  @param ido - index of internal variables for given material in the ipp other array.
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 06.2003
  Modified by Tomas Koudelka, 06.2019
*/
void scaldamcc::updateval(long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}



/**
  This function returns the value of the limit elastic strain.

  @param ipp - integration point pointer

  @return The function returns actual value of limit elastic strain.

  Created by Tomas Koudelka, 10.2004
  Modified by Tomas Koudelka, 06.2019
*/
double scaldamcc::epsefunction (long ipp)
{
  double eps0 = 0.0;
  double e = 0.0;

  switch(damfunct){
    case simpexp:
      e = Mm->give_actual_ym(ipp);
      eps0 = ft/e;
      break;
    case mazarsexp:
      eps0 = k0;
      break;
    default:
      print_err("unknown damage function type %d is required", __FILE__, __LINE__, __func__, int(damfunct));
      abort();
  }
  return eps0;
}
