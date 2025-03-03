#include "scaldam.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>
#include <float.h>

#define nijac 20
#define limit 1.0e-8



/**
  The constructor inializes attributes to zero values.
  
  Created by Tomas Koudelka,
*/
scaldam::scaldam (void)
{
  ft=0.0;  uf = 0.0; k = 0.0; cde = corr_on;
  min_exp_arg = log(DBL_MIN);
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka,
*/
scaldam::~scaldam (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void scaldam::read (XFILE *in)
{
  //  types of norm
  //  norstrain = 1
  //  norenergy = 2
  //  norposstrain = 3
  //  norposenergy = 4
  //  norrankine = 5
  //  norrankinesmooth = 6
  //  normazar = 7
  //  vonmises= 8
  xfscanf (in,"%k%lf %k%lf %k%m %k%m","ft",&ft,"uf",&uf,"norm_type",&paramf_type_kwdset,(int *)(&ftype), 
           "cor_dis_energy", &corr_disip_en_kwdset, (int *)(&cde));
  if (ftype == vonmises)
    xfscanf(in, "%k%lf", "k", &k);
  
  if (cde != corr_off)
    sra.read (in);
/*
  xfscanf(in, "%k%ld", "temp_evf_ft", &cftt);
  if (cftt)
    ft_temp.read(in);
*/
}



/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened output text file
  
  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void scaldam::print (FILE *out)
{
  fprintf (out,"%le %le %d %d ",ft,uf,ftype,cde);
  if (ftype == vonmises)
    fprintf(out, " %e ",k);
  
  if (cde != corr_off)
    sra.print (out);
}



/**
  The function computes parameters for the damage function
  Different types of the parameter computing are given by the
  attribute ft.

  @param  ipp - integration point number
  @param  eps - %vector of the strains
  @param  kappa - %vector where the computed parameters are stored (output)
  
  @return The function returns actual %vetcor of parameters of damage function 
          in the parameter kappa.

  Created by Tomas Koudelka,
*/
void scaldam::damfuncpar(long ipp, vector &eps, vector &kappa)
{
  double epseq;
  vector poseps;
  vector psig;
  vector peps;
  vector pospsig;
  vector tmp;
  vector sig;
  matrix pvect;
  matrix sigt;
  matrix epst;
  vector epsf;
  matrix dev;
  matrix d;
  double i1e, j2e, a, b, e, nu, temp; 
  long ncomp;

  switch (ftype)
  {
    case norstrain:
      // strain norm = |eps| = \sqrt{eps_{ij} eps_{ij}}
      kappa[0] = tensor_strain_norm(eps);
      break;
    case norenergy:
      // energy norm c\sqrt{eps : D_e : eps} where c should be E^{-1/2}
      ncomp=Mm->ip[ipp].ncompstr;
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocv(RSTCKVEC(ncomp, tmp));
      elmatstiff (d,ipp);
      mxv(d, eps, tmp);
      scprd(tmp, eps, epseq);
      e = Mm->give_actual_ym(ipp);
      kappa[0] = sqrt(epseq/e);
      break;
    case norposstrain:
      // positive strain norm = \sqrt{<eps_{ij}><eps_{ij}>}
      reallocv(RSTCKVEC(eps.n,poseps));
      extractposv (eps, poseps);
      kappa[0] = tensor_strain_norm(poseps);
      break;
    case norposenergy:
      // positive energy norm c\sqrt{<eps> : D_e : <eps>} where c should be E_oed^{-1/2}
      ncomp=Mm->ip[ipp].ncompstr;
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocv(RSTCKVEC(ncomp, tmp));
      elmatstiff (d,ipp);
      reallocv(RSTCKVEC(eps.n,poseps));
      extractposv (eps, poseps);
      mxv(d, poseps, tmp);
      scprd(tmp, poseps, epseq);
      kappa[0] = sqrt(epseq/d(0,0));
      break;
    case norrankine:
    // Rankine norm = max(\bar{sig}_1, \bar{sig}_2, \bar{sig}_3)/E
    // where \bar{\sig_{\alpha}} are the effective principal stress components
    // computed from \bar{sig} = D_e : eps
    {
      e = Mm->give_actual_ym(ipp);
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(ncomp, sig));
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocm(RSTCKMAT(3, 3, sigt));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, psig));
      elmatstiff (d,ipp);
      mxv(d, eps, sig);
      vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
      princ_val (sigt,psig,pvect,nijac,limit,Mp->zero,3,1);
      // principal values are sorted -> maximum value of principal stresses is at the last
      // position of vector princsig
      if (psig[2] > 0.0)
        kappa[0] = psig[2] / e;
      else
        kappa[0] = 0.0;
      break;
    }
    case norrankinesmooth:
    // smoothed Rankine norm = \sqrt{<\bar{sig}_1>^2 + <\bar{sig}_2>^2 + <\bar{sig}_3>^2}/E
    // where \bar{\sig_{\alpha}} is the effective principal stress component
    // computed from \bar{sig} = D_e : eps
    {
      e = Mm->give_actual_ym(ipp);
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(ncomp, sig));
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocm(RSTCKMAT(3, 3, sigt));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, psig));
      reallocv(RSTCKVEC(3, pospsig));
      elmatstiff (d,ipp);
      mxv(d, eps, sig);
      vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
      princ_val (sigt,psig,pvect,nijac,limit,Mp->zero,3,1);
      extractposv (psig, pospsig);
      scprd(pospsig, pospsig, epseq);
      kappa[0] = sqrt(epseq) / e;
      break;
    }
    case normazar:
    // Mazar's norm = \sqrt{<eps_1>^2 + <eps_2>^2 + <eps_3>^2}
    // where eps_{alpha} is the principal strain component
    {
      reallocm(RSTCKMAT(3, 3, epst));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, peps));
      reallocv(RSTCKVEC(3,poseps));
      nullv(peps);
      vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
      princ_val (epst,peps,pvect,nijac,limit,Mp->zero,3,1);
      extractposv (peps, poseps);
      scprd(poseps, poseps, epseq);
      kappa[0] = sqrt(epseq);
      break;
    }
    case vonmises:
      // glasgow modified von Mises norm
      ncomp=Mm->ip[ipp].ncompstr;
      //allocm(3, 3, epst);
      reallocv(RSTCKVEC(6, epsf));
      reallocm(RSTCKMAT(3, 3, dev));
      give_full_vector (epsf,eps,Mm->ip[ipp].ssst);
      i1e = first_invar(epsf);
      j2e = j2_strain_invar (epsf);
      e = Mm->give_actual_ym(ipp);
      nu = Mm->give_actual_nu(ipp);
      a   = (k-1.0)/(2.0*k)/(1.0-2.0*nu);
      b   = 3.0/k/(1+nu)/(1.0+nu);
      temp = a*a*i1e*i1e + b*j2e;
      kappa[0] = a*i1e+sqrt(temp);
/*
      // original von Mises norm, it is the same as the above one
      kappa[0]  = (k-1.0)/(2.0*k)/(1.0-2.0*nu)*i1e;
      tmp = (k-1.0)*(k-1.0)/(1.0-2.0*nu)/(1.0-2.0*nu)*i1e*i1e -12.0*k/(1.0+nu)/(1.0+nu)*j2e;
      kappa[0] += 1.0/(2.0*k)*sqrt(tmp);*/
      break;
    default:
      print_err("unknown damage parameter function type", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function computes derivatives of the damage function parameter with respect to 
  strain components. Different types of the damage function parameters are given by the
  attribute ft.

  @param  ipp[in]    - integration point number
  @param  eps[in]    - %vector of the strains
  @param  kappa[out] - %vector where the computed derivatives are stored
  
  @return The function returns actual %vetcor of parameters of damage function 
          in the parameter kappa.

  Created by Tomas Koudelka,
*/
void scaldam::derdamfuncpar(long ipp, vector &eps, vector &dkappa)
{
  long i;
  double epseq;
  vector poseps;
  vector psig;
  vector peps;
  vector pospsig;
  vector tmp;
  vector sig;
  matrix pvect;
  matrix sigt;
  matrix epst;
  vector epsf;
  vector dev;
  matrix d;
  matrix dkappat, tmpm;
  matrix em[3];
  double i1e, j2e, a, b, e, nu, temp;
  long ncomp, nncomp;

  switch (ftype)
  {
    case norstrain:
      // strain norm = |eps| = \sqrt{eps_{ij} eps_{ij}}
      temp = tensor_strain_norm(eps);
      copyv(eps, dkappa);
      cmulv(1.0/temp, dkappa);
      ncomp=Mm->ip[ipp].ncompstr;
      nncomp = Mm->give_num_norm_comp(ipp);
      for(i=nncomp; i<ncomp; i++)
        dkappa[i] *= 0.5;
      break;
    case norenergy:
      // energy norm c\sqrt{eps : D_e : eps} where c should be E^{-1/2}
      ncomp=Mm->ip[ipp].ncompstr;
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      makerefv(sig, dkappa);
      elmatstiff (d,ipp);
      mxv(d, eps, sig);
      scprd(sig, eps, epseq);
      epseq = sqrt(epseq);
      e = Mm->give_actual_ym(ipp);
      cmulv(e*epseq, dkappa);
      break;
    case norposstrain:
      // positive strain norm = \sqrt{<eps_{ij}><eps_{ij}>}
      makerefv(poseps, dkappa);
      extractposv (eps, poseps);      
      epseq = tensor_strain_norm(poseps);
      cmulv(1.0/epseq, dkappa);
      ncomp=Mm->ip[ipp].ncompstr;
      nncomp = Mm->give_num_norm_comp(ipp);
      for(i=nncomp; i<ncomp; i++)
        dkappa[i] *= 0.5;      
      break;
    case norposenergy:
      // positive energy norm c\sqrt{<eps> : D_e : <eps>} where c should be E_oed^{-1/2}
      ncomp=Mm->ip[ipp].ncompstr;
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      makerefv(sig, dkappa);
      elmatstiff (d,ipp);
      reallocv(RSTCKVEC(eps.n,poseps));
      mxv(d, poseps, sig);
      scprd(sig, poseps, epseq);
      epseq = sqrt(epseq);
      temp = d(0,0)*epseq;
      temp = 1.0/temp;
      cmulv(temp, dkappa);
      extractposv (dkappa, dkappa);   
      break;
    case norrankine:
    // Rankine norm = max(\bar{sig}_1, \bar{sig}_2, \bar{sig}_3)/E
    // where \bar{\sig_{\alpha}} are the effective principal stress components
    // computed from \bar{sig} = D_e : eps
    {
      e = Mm->give_actual_ym(ipp);
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(ncomp, sig));
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocm(RSTCKMAT(3, 3, sigt));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, psig));
      elmatstiff (d,ipp);
      mxv(d, eps, sig);
      vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
      princ_val (sigt, psig, pvect, nijac, limit, Mp->zero, 3, 1);
      // principal values are sorted -> maximum value of principal stresses is at the last
      // position of vector psig
      nullv(dkappa);
      if (psig[2] > 0.0){
        e = Mm->give_actual_ym(ipp);
        nu = Mm->give_actual_nu(ipp);
        reallocv(RSTCKVEC(3, tmp));
        for (i=0; i<3; i++){
          reallocm(RSTCKMAT(3, 3, em[i]));
          extractcol(pvect, i, tmp);
          vxv(tmp, tmp, em[i]);
        }
        reallocm(RSTCKMAT(3,3,tmpm));
        addm(em[0], em[1], tmpm);
        temp = (1+nu)*(1-2.0*nu);
        cmulm(nu/temp, tmpm);
        cmulm((1.0-nu)/temp, em[2]);
        addm(em[2], tmpm, tmpm);
        tensor_vector(dkappa, tmpm, Mm->ip[ipp].ssst, stress); // strain option must be reconsidered in the last argument
      }
      break;
    }
    case norrankinesmooth:
    // smoothed Rankine norm = \sqrt{<\bar{sig}_1>^2 + <\bar{sig}_2>^2 + <\bar{sig}_3>^2}/E
    // where \bar{\sig_{\alpha}} is the effective principal stress component
    // computed from \bar{sig} = D_e : eps
    {
      e = Mm->give_actual_ym(ipp);
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(ncomp, sig));
      reallocm(RSTCKMAT(ncomp, ncomp, d));
      reallocm(RSTCKMAT(3, 3, sigt));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, psig));
      elmatstiff (d,ipp);
      mxv(d, eps, sig);
      vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
      princ_val(sigt, psig, pvect, nijac, limit, Mp->zero, 3, 0);
      extractposv (psig, pospsig);
      nullv(dkappa);
      epseq = normv(pospsig);
      if (epseq > 0.0){
        reallocv(RSTCKVEC(3, tmp));
        reallocm(RSTCKMAT(3, 3, dkappat));
        for (i=0; i<3; i++){
          reallocm(RSTCKMAT(3, 3, em[i]));
          extractcol(pvect, i, tmp);
          vxv(tmp, tmp, em[i]);
        }
        e = Mm->give_actual_ym(ipp);
        nu = Mm->give_actual_nu(ipp);
        temp = (1+nu)*(1-2.0*nu);
        reallocm(RSTCKMAT(3,3,tmpm));
        for (i=0; i<3; i++){
          if (pospsig[i] > 0){
            addm(em[(i+1)%3], em[(i+2)%3], tmpm);
            addmultm(em[i], (1.0-nu)/temp, tmpm, nu/temp, dkappat);
          }
        }
        cmulm(1.0/epseq, dkappat);
        tensor_vector(dkappa, dkappat, Mm->ip[ipp].ssst, stress); // strain option must be reconsidered in the last argument
      }
      break;
    }
    case normazar:
    // Mazar's norm = \sqrt{<eps_1>^2 + <eps_2>^2 + <eps_3>^2}
    // where eps_{alpha} is the principal strain component
    {
      reallocm(RSTCKMAT(3, 3, epst));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, peps));
      reallocv(RSTCKVEC(3,poseps));
      nullv(peps);
      vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
      princ_val (epst, peps, pvect, nijac, limit, Mp->zero, 3, 0);
      extractposv (peps, poseps);
      epseq = normv(poseps);
      nullv(dkappa);
      if (epseq > 0.0){
        reallocv(RSTCKVEC(3, tmp));
        reallocm(RSTCKMAT(3, 3, dkappat));
        for (i=0; i<3; i++){
          if (poseps[i] > 0.0){
            reallocm(RSTCKMAT(3, 3, em[i]));
            extractcol(pvect, i, tmp);
            vxv(tmp, tmp, em[i]);
            addm(em[i], dkappat, dkappat);
          }
        }
        cmulm(1.0/epseq, dkappat);
        tensor_vector(dkappa, dkappat, Mm->ip[ipp].ssst, stress); // strain option must be reconsidered in the last argument
      }
      break;
    }
    case vonmises:
      // glasgow modified von Mises norm: \kappa = a.I_{1\eps} + sqrt(a^2.I^2_{1\eps} + b.J_{2\eps})
      // dkappa/deps = a.delta_{ij} + \frac{1}{2.sqrt(a^2.I^2_{1\eps} +
      //               b.J_{2\eps})}[2a^2\delta_{ij} + b.s_e_{ij}:(I_{ijkl} - 1/3\delta_{ij}\delta{kl})]
      ncomp=Mm->ip[ipp].ncompstr;
      reallocv(RSTCKVEC(6, epsf));
      reallocv(RSTCKVEC(6, dev));
      give_full_vector (epsf, eps, Mm->ip[ipp].ssst);
      i1e = first_invar(epsf);
      deviator (epsf, dev);
      j2e = j2_strain_invar (epsf);
      e = Mm->give_actual_ym(ipp);
      nu = Mm->give_actual_nu(ipp);
      a   = (k-1.0)/(2.0*k)/(1.0-2.0*nu);
      b   = 3.0/k/(1+nu)/(1.0+nu);
      temp = a*a*i1e*i1e + b*j2e;
      temp = sqrt(temp);
      // compute I - 1/3(\delta*\delta) = \delta_{ik}\delta_{jl} - 1/3\delta_{ij}\delta_{kl}
      reallocm(RSTCKMAT(6, 6, d));
      d(0,0) = d(1,1) = d(2,2) = 2.0/3.0;
      d(3,3) = d(4,4) = d(5,5) = 1.0;
      d(0,1) = d(0,2) = d(1,0) = d(1,2) = d(2,0) = d(2,1) = -1.0/3.0;
      // compute b.s_e:(I-1/3(\delta*\delta))
      reallocv(RSTCKVEC(6, tmp));
      mxv(d, dev, tmp);
      cmulv(b, tmp);
      // add 2a^2\delta
      tmp(0) += 2.0*a*a;
      tmp(1) += 2.0*a*a;
      tmp(2) += 2.0*a*a;
      // scale by 1/(2.sqrt(a^2I^2_{1\eps} + bJ_{2\eps}))
      cmulv(1.0/(2.0*temp), tmp);
      // add a\delta
      tmp(0) += a;
      tmp(1) += a;
      tmp(2) += a;
      // convert from the full vector length to the reduced vector format
      give_red_vector(tmp, dkappa, Mm->ip[ipp].ssst);
      break;
    default:
      print_err("unknown damage parameter function type", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function computes damage parameter omega which is the result of the
  damage function.

  @param ipp - integration point number
  @param kappa - %vector of the parameters of damage function
  @param omegao - %vector of attained damage parameters
  
  @return Function returns value of the damage parameter omega.
  
  Created by Tomas Koudelka,
*/
double scaldam::damfunction(long ipp, vector &kappa, vector &omegao)
{
  double err,omega = 0.0, dtmp, tmp, h, e, uf_a, aft;
  long i,ni, eid;
  double indam;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  e = Mm->give_actual_ym(ipp);
  // tensile strength
  aft = Mm->give_actual_ft(ipp);
  // initial equivalent strain treshold
  indam = ft/e;
  // actual value of initial crack opening
  uf_a = uf*(aft/ft);
  
  if (kappa[0] < indam)
    // elastic loading
    return 0.0;
  // cut the omega value, speed-up the computation and prevent math domain error of exp() function
  if (1.0-omegao[0] < err){
    return 1.0;
  }
  if ((Mm->ip[ipp].hmt & 2) > 0)
  {
    omega = 1.0;
    if (-(kappa[0]-aft/e)/(uf_a-aft/e) > min_exp_arg){ // prevent exponent overflow of exp() function
      omega -= (aft/e/kappa[0])*exp(-(kappa[0]-aft/e)/(uf_a-aft/e));
    }
    if ((omega > 1.0) || (omega < 0.0))
    {
      print_err("computed omega=%le is out of prescribed range <0,1>\n"
                " elem %ld, ip %ld:\n"
                " kappa=%le, aft=%le, e=%le, uf=%le",
                __FILE__, __LINE__, __func__, omega, Mm->elip[ipp]+1, ipp+1, kappa[0], aft, e, uf_a);
      abort();
    }
    return omega;
  }
  // correction of dissipated energy - mesh independent scalar damage
  eid = Mm->elip[ipp];
  switch (Mt->give_dimension(eid))
  {
    case 1:
      h = Mt->give_length(eid);
      break;
    case 2:
      h = sqrt(Mt->give_area(eid));
      break;
    case 3:
      h = pow(Mt->give_volume(eid), 1.0/3.0);
      break;
    default:
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
  }
  if (cde == corr_off)
  // no correction of dissipated energy only simple scalar damge
  { 
    omega = 1.0;
    if (-(kappa[0]-aft/e)/(uf_a-aft/e) > min_exp_arg){ // prevent exponent overflow of exp() function
      omega -= (aft/e/kappa[0])*exp(-(kappa[0]-aft/e)/(uf_a-aft/e));
    }
    if ((omega > 1.0) || (omega < 0.0))
    {
      print_err("computed omega=%le is out of prescribed range <0,1>\n"
		" elem %ld, ip %ld:\n"
		" kappa=%le, aft=%le, e=%le, uf=%le",
		__FILE__, __LINE__, __func__, omega, Mm->elip[ipp]+1, ipp+1, kappa[0], aft, e, uf_a);
      abort();
    }
    return omega;
  }
  // it is better to start from omega=1 in order to avoid of convergention problems
  omega = 1.0;
  tmp = 0.0;
  if (-h*omega*kappa[0]/uf_a > min_exp_arg){ // prevent exponent overflow of exp() function
    for (i = 0; i < ni; i++)
    {
      dtmp = -e*kappa[0]+aft/uf_a*h*kappa[0]*exp(-h*omega*kappa[0]/uf_a);
      tmp = (1.0-omega)*e*kappa[0]-aft*exp(-h*omega*kappa[0]/uf_a);
      if (fabs(tmp) < err*aft){
        omega += - tmp/dtmp;
        break;
      }
      omega += - tmp/dtmp;
    }
    if ((omega > 1.0) || (omega < 0.0) || (fabs(tmp) > aft*err))
    // convergention problems have been encountered
    // starting from previously attained omega can help
    {
      omega = omegao[0];
      for (i = 0; i < ni; i++){
        dtmp = -e*kappa[0]+aft/uf_a*h*kappa[0]*exp(-h*omega*kappa[0]/uf_a);
        tmp = (1.0-omega)*e*kappa[0]-aft*exp(-h*omega*kappa[0]/uf_a);
        if (fabs(tmp) < err*aft){
          omega += - tmp/dtmp;
          break;
        }
        omega += - tmp/dtmp;
      }
    }
  }
  
  if (fabs(tmp) > aft*err)
  {
    print_err("iteration of omega doesn't converge with required error\n"
              " elem %ld, ip %ld:\n"
              " kappa=%le, aft=%le, e=%le, uf=%le, h=%le",
              __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1, kappa[0], aft, e, uf_a, h);
    abort();
  }
  
  if ((omega > 1.0) || (omega < 0.0))
  {
    print_err("computed omega=%le is out of prescribed range <0,1>\n"
              " elem %ld, ip %ld:\n"
              " kappa=%le, aft=%le, e=%le, uf=%le, h=%le",
              __FILE__, __LINE__, __func__, omega, Mm->elip[ipp]+1, ipp+1, kappa[0], aft, e, uf_a, h);
    abort();
  }
  
  return omega;
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  @return The function returns computed elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka,
*/
void scaldam::elmatstiff (matrix &d,long ipp)
{
  double e, e0;
  long idem = Mm->ip[ipp].gemid();         // index of elastic material
  long idoem=Mm->givencompeqother(ipp,0);  // total number of internal variables
  long tmp=Mm->givencompeqother(ipp,idem); // number of internal variables for elastic material and thermal material 
  idoem -= tmp;
  
  Mm->elmatstiff (d,ipp);
  e = Mm->give_actual_ym(ipp); // actual Young modulus (it can be computed also in creep model)
  e0 = Mm->give_actual_ym(ipp,idem,idoem); // Young modulus from elastic material

  cmulm(e/e0, d);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns computed stiffness %matrix in the parameter d.
  
  Created by Tomas Koudelka,
*/
void scaldam::matstiff (matrix &d,long ipp,long ido)
{
  double dp;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      elmatstiff (d,ipp);
      break;
    case secant_stiff:
      elmatstiff (d,ipp);
      dp=Mm->ip[ipp].eqother[ido+1];
      if (dp > 0.999)
        dp = 0.999;
      cmulm (1.0-dp,d);
      break;
    case tangent_stiff:
    case incr_tangent_stiff:
    case ijth_tangent_stiff:
      elmatstiff (d,ipp);
      dp=Mm->ip[ipp].eqother[ido+1];
      if (dp > 0.999)
        dp = 0.999;
      cmulm (1.0-dp,d);
      break;
    default:
      print_err("unknown type of stifness matrix is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes actual stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 7.10.2001
*/
void scaldam::nlstresses (long ipp, long im, long ido)
{
  long i,ncomp=Mm->ip[ipp].ncompstr;
  double dp, e, nu, u_crack;
  vector epsn(ncomp),sigma(ncomp), kappa(1), omegao(1);
  matrix d(ncomp, ncomp);


  e = Mm->give_actual_ym(ipp);
  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }
  //  initial values of total strains
  for (i=0;i<ncomp;i++){
    epsn[i] = Mm->ip[ipp].strain[i];
  }
  kappa[0]  = Mm->ip[ipp].eqother[ido+0];
  omegao[0] = Mm->ip[ipp].eqother[ido+1];

  elmatstiff(d,ipp);
  // damage stress solver
  dp = Mm->scal_dam_sol (ipp, im, ido, epsn, kappa, omegao, sigma,d);
  // dp = Mm->scal_dam_sol2 (ipp, im, ido, epsn, kappa, omegao, sigma,d);

  //  new stress data storage
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].stress[i]=sigma[i];
  }
  // computation of approximate value of crack opening
  u_crack = dp*kappa[0]*Mm->give_proczonelength(ipp);
  // storage of actual values of internal variables
  Mm->ip[ipp].other[ido+0]=kappa[0];
  Mm->ip[ipp].other[ido+1]=dp;
  Mm->ip[ipp].other[ido+2]=u_crack;
}



/**
  The function updates values in the eqother array containing values from the previous equlibrium state 
  to the values attained in the actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void scaldam::updateval (long ipp,long im,long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }

}



/**
  The function returns the actual value of tensile strength.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns actual value of tensile strength.

  Created by Tomas Koudelka,
*/
double scaldam::give_actual_ft (long /*ipp*/, long /*im*/, long /*ido*/)
{
/* 
  double ret;
  if (cftt)
  {
    ret = ft_temp.getval(Mm->givenonmechq(temperature, ipp));
  }
*/
  return ft;
}



/**
  This function returns the value of the limit elastic strain which is dependent on
  actual value of Young modulus and actual value of tensile strength.

  @param ipp - integration point number in the mechmat ip array.

  @return Function returns actual value of limit elastic strain
  
  Created by Tomas Koudelka,
*/
double scaldam::epsefunction (long ipp)
{
  double e = Mm->give_actual_ym(ipp);
  double ft = Mm->give_actual_ft(ipp);

  return (ft /e) ;
}



/**
  This function returns the value of damage parameter obtained from the actual state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the value of damage parameter omega.

  Created by Tomas Koudelka,
*/
double scaldam::givedamage (long ipp, long ido)
{
  return Mm->ip[ipp].other[ido+1];
}

/**
  This function returns the length of process zone

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the length of process zone

  Created by Tomas Koudelka,
*/
double scaldam::give_proczonelength (long ipp, long /*ido*/)
{
  long eid;
  double h;
  
  // computation of approximate value of crack opening
  eid = Mm->elip[ipp];
  switch (Mt->give_dimension(eid))
  {
    case 1:
      h = Mt->give_length(eid);
      break;
    case 2:
      h = sqrt(Mt->give_area(eid));
      break;
    case 3:
      h = pow(Mt->give_volume(eid), 1.0/3.0);
      break;
    default:
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
  }
  
  return h;
}


/**
  This function returns the crack width from the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the crack width

  Created by Tomas Koudelka,
*/
double scaldam::give_crackwidth (long ipp, long ido)
{
  return Mm->ip[ipp].other[ido+2];
}

/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void scaldam::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      ft=val[i];
      break;
    }
    case 1:{
      uf=val[i];
      break;
    }
    default:{
      print_err("wrong index of changed atribute",__FILE__,__LINE__,__func__);
    }
    }
  }
}

