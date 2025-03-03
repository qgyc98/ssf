#include "ortodam.h"
#include "bmatrix.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>
#include <stdlib.h>

#define nijac 200
#define limit 1.0e-15  // zero limit for principal values of tensors (Jacobi method)
#define qlimit 1.0e-8  // zero limit for switch from quadratuc Bezier to pure quadrat
#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif
#define sqr(x) (x)*(x)



/**
  The constructor inializes attributes to zero values.
  
  Created by Tomas Koudelka, 9.2006
  Modified Tomas Koudelka, 10.2008
*/
ortodam::ortodam (void)
{
  cde = corr_off;
  fat = fatigue_off;
  ac = 0.0;  bc = 0.0;  y0c = 0.0;  at = 0.0;  bt = 0.0;  y0t = 0.0;
  gfc = 0.0; gft = 0.0;
  ft = uft = ult = 0.0;
  fc = ufc = ulc = 0.0;
  pqt = pqc = 0;
  betac = betat = 0.0;
}



/**
  The destructor is only for the formal purposes.
  
  Created by Tomas Koudelka, 9.2006
*/
ortodam::~ortodam (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void ortodam::read (XFILE *in)
{
  xfscanf (in, "%k%m", "dam_evol_f", &dam_evolfunc_kwdset, (int *)&damevf);
  switch(damevf)
  {
    case brittle:
      xfscanf (in,"%k%lf %k%lf", "f_t", &ft, "u_ft", &uft);
      xfscanf (in,"%k%lf %k%lf", "f_c", &fc, "u_fc", &ufc);
      xfscanf (in,"%k%m", "corr_dissip", &corr_disip_en_kwdset, (int *)(&cde));
      if (cde != corr_off)
        sra.read (in);      
      break;
    case quasi_brittle:
      xfscanf (in, "%k%m", "corr_dissip", &corr_disip_en_kwdset, (int *)(&cde));
      switch (cde)
      {
        case corr_off:
          xfscanf (in, "%k%lf %k%lf %k%lf", "A_c", &ac, "B_c", &bc, "Y_0c", &y0c);
          xfscanf (in, "%k%lf %k%lf %k%lf", "A_t", &at, "B_t", &bt, "Y_0t", &y0t);
          break;
        case corr_on:
          xfscanf (in, "%k%lf %k%lf %k%lf", "G_fc", &gfc, "B_c", &bc, "Y_0c", &y0c);
          xfscanf (in, "%k%lf %k%lf %k%lf", "G_ft", &gft, "B_t", &bt, "Y_0t", &y0t);
          break;
        default:
          print_err(" unknown type of correction of disipated energy is required", __FILE__, __LINE__,__func__);
      }
      break;
    case quadbezier:
      xfscanf (in,"%k%lf %k%lf %k%lf", "f_t", &ft, "u_ft", &uft, "u_limt", &ult);
      if (fabs(uft/ult-0.5) < qlimit)
        pqt = 1;
      xfscanf (in,"%k%lf %k%lf %k%lf", "f_c", &fc, "u_fc", &ufc, "u_limc", &ulc);
      if (fabs(ufc/ulc-0.5) < qlimit)
        pqc = 1;
      xfscanf (in,"%k%m", "corr_dissip", &corr_disip_en_kwdset, (int *)(&cde));
      if (cde != corr_off)
        sra.read (in);
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
  xfscanf (in,"%k%m", "dam_eq_strain", &paramf_type_kwdset, (int *)(&dameqstr));
  xfscanf (in,"%k%m", "fatigue", &fatigue_flag_kwdset, (int *)(&fat));
  if (fat == fatigue_on)
    xfscanf (in,"%k%lf %k%lf", "beta_t", &betat, "beta_c", &betac);
}



/**
  The function computes the equivalent strain norm of required type.
  Only strain or energy type of norm are accepted.

  @param nu   - Poissons ratio
  @param i    - principal direction id
  @param peps - %vector of principal strains

  @return The function returns calculated equivalent strain norm.

  Created by Tomas Koudelka 07.2011
*/
double ortodam::dam_eq_strain(double nu, long i, vector &peps)
{
  double kappa = 0.0;
  switch (dameqstr)
  {
    case norstrain:
      // for tension
      if (peps(i) > 0.0)
        kappa = peps(i);
      // for compression
      if (peps(i) < 0.0)
        kappa = fabs(peps(i));
      break;
    case norenergy:
      // for tension
      if (peps(i) > 0.0)
      {
        kappa  = (1-nu)*peps(i)*peps(i);     // for i=0: dsig_I*eps_I/dD_I = ((1-nu)*(<eps_I>)^2 +
        if (peps((i+1)%3) > 0.0)
          kappa += nu*peps(i)*peps((i+1)%3); //                             + nu*<eps_I>*<eps_II>
        if (peps((i+2)%3) > 0.0)
          kappa += nu*peps(i)*peps((i+2)%3); //                             + nu*<eps_I>*<eps_III>) *
        kappa /= (1-nu);                     //                             * (E/((1-2*nu)*(1+nu)))
//        kappa /= (1-2*nu)*(1+nu);            //                          * (E/((1-2*nu)*(1+nu)))
//        kappa = sqrt(kappa);                 // kappa_I = sqrt(dsig_I*eps_I/dD_I / E)
        kappa = sqrt(kappa);                 // kappa_I = sqrt(dsig_I*eps_I/dD_I / (E*(1-nu)/((1-2*nu)*(1+nu))))
      }
      // for compression
      if (peps(i) < 0.0)
      {
        kappa  = (1-nu)*peps(i)*peps(i);     // dsig_I*eps_I/dD_I = ((1-nu)*(<-eps_I>)^2 +
        if (peps((i+1)%3) < 0.0)
          kappa += nu*peps(i)*peps((i+1)%3); //                    + nu*<-eps_I>*<-eps_II>
        if (peps((i+2)%3) < 0.0)
          kappa += nu*peps(i)*peps((i+2)%3); //                    + nu*<-eps_I>*<-eps_III>) *
        kappa /= (1-nu);                     //                    * (E/((1-2*nu)*(1+nu)))
//        kappa /= (1-2*nu)*(1+nu);            //                          * (E/((1-2*nu)*(1+nu)))
//        kappa = sqrt(kappa);                 // kappa_I = sqrt(dsig_I*eps_I/dD_I / E)
        kappa = sqrt(kappa);                 // kappa_I = sqrt(dsig_I*eps_I/dD_I / (E*(1-nu)/((1-2*nu)*(1+nu))))
      }
      break;
    case norrankine:
      // for tension and compression
      kappa  = (1-nu)*peps(i) + nu*(peps((i+1)%3) + peps((i-1)%3));     //        effective sig(i)/E
      kappa /= (1-2*nu)*(1+nu);
      kappa = fabs(kappa);
      break;
    default:
      print_err(" unknown type of equivalent strain norm is required", __FILE__, __LINE__,__func__);
  }
  return kappa;
}



/**
  Function computes value of load function. 

  @param  ipp  - integration point number
  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  damt - %vector actual values of principal damage parameters for tension
  @param  damc - %vector actual values of principal damage parameters for compression
  @param  aat  - actual values of material parameter a for tension 
  @param  aac  - actual values of material parameter a for compression 
  @param  lft  - %vector of resulting values of load function for tension for each principal direction (output)
  @param  lfc  - %vector of resulting values of load function for compression for each principal direction (output)

  @return The function returns load function values separately for tension and compression
          in the parameters lft and lfc for each principal direction.

  Created by Tomas Koudelka, 9.2006
  Modified Tomas Koudelka, 10.2008
  Modified Tomas Koudelka, 07.2011
*/
void ortodam::loadfunc(long ipp, double nu, vector &peps, vector &damt, vector &damc, 
                       double aat, double aac, vector &lft, vector &lfc)
{
  long i;
  double e, af, auf, aul, kappa, tmp1, tmp2, a, b, h;
  vector pyt(3), pyc(3);

  nu = Mm->give_actual_nu(ipp);
  switch(damevf)
  {
    case brittle:
      for (i=0; i<3; i++)
      {
        // for tension
        if (peps(i) >= 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          e = Mm->give_actual_ym(ipp);
          af = Mm->give_actual_ft(ipp);
          if (kappa <= af/e) 
            lft(i) = -1.0;      
          else
          {
            e = Mm->give_actual_ym(ipp);
            // actual tensile strength
            af = Mm->give_actual_ft(ipp);
            // actual value of initial crack opening
            auf = uft*(af/ft);
            tmp1 = -(kappa-af/e)/(auf-af/e);
            lft(i) = 1.0 - af/(e*kappa)*exp(tmp1) - damt(i);
          }
        }
        // for compression
        if (peps(i) < 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          e = Mm->give_actual_ym(ipp);
          af = Mm->give_actual_fc(ipp);
          if (kappa <= af/e) 
            lfc(i) = -1.0;      
          else
          {
            e = Mm->give_actual_ym(ipp);
            // actual compressive strength
            af = Mm->give_actual_fc(ipp);
            // actual value of initial crack opening
            auf = ufc*(af/fc);
            tmp1 = -(kappa-af/e)/(auf-af/e);
            lfc(i) = 1.0 - af/(e*kappa)*exp(tmp1) - damc(i);
          }
        }
      }
      break;
    case quasi_brittle:
      for (i=0; i<3; i++)
      {
        // for tension
        if (peps(i) >= 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          if (kappa <= y0t) 
            lft(i) = -1.0;      
          else
          {
            // compute (yt-y0t)^bt
            tmp1 = kappa - y0t;
            tmp2 = log(tmp1);
            tmp1 = exp(bt*tmp2);
            // compute value of load function for tension in given principal direction i
            lft(i) = (1.0 - damt(i))*(1.0+aat*tmp1) - 1.0;
          }
        }
        // for compression
        if (peps(i) < 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          if (kappa <= y0c) 
            lfc(i) = -1.0;      
          else
          {
            // compute (yc-y0c)^bc
            tmp1 = kappa - y0c;
            tmp2 = log(tmp1);
            tmp1 = exp(bc*tmp2);
            // compute value of load function for compression in given principal direction i
            lfc(i) = (1.0 - damc(i))*(1.0+aac*tmp1) - 1.0;
          }
        }
      }
      break;
    case quadbezier:
      h = 1.0; // here should be size of element but it is not neccesary
      for (i=0; i<3; i++)
      {
        // for tension
        if (peps(i) >= 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          if (kappa <= y0t) 
            lft(i) = -1.0;      
          else
          {
            e = Mm->give_actual_ym(ipp);
            // actual tensile strength
            af = Mm->give_actual_ft(ipp);
            // actual value of initial crack opening
            auf = uft*(af/ft);
            // actual value of limit crack opening
            aul = ult*(af/ft);
            a = 4.0*auf*auf;
            b = (aul-2.0*auf);
            if (pqt)
              lft(i) = e*kappa-damt(i)*e*kappa-af*(1.0-damt(i)*kappa*h/auf+sqr(damt(i)*kappa*h)*(aul-auf)/(auf*aul*aul));
            else
              // following expression should be positive
              lft(i) = e*kappa-damt(i)*e*kappa-af*(1.0-auf*(2.0*aul-3.0*auf)/(b*b)-2.0*(aul-auf)/(b*b)*sqrt(a+b*damt(i)*h*kappa)+damt(i)*h*kappa/b);
          }
        }
        // for compression
        if (peps(i) < 0.0)
        {
          kappa = dam_eq_strain(nu, i, peps);
          if (kappa <= y0c) 
            lfc(i) = -1.0;      
          else
          {
            e = Mm->give_actual_ym(ipp);
            // actual compressive strength
            af = Mm->give_actual_fc(ipp);
            // actual value of initial crack opening
            auf = ufc*(af/fc);
            // actual value of limit crack opening
            aul = ult*(af/ft);
            a = 4.0*auf*auf;
            b = (aul-2.0*auf);
            if (pqc)  // pure quadratic function have to be used
              lfc(i) = e*kappa-damc(i)*e*kappa-af*(1.0-damc(i)*kappa*h/auf+sqr(damc(i)*kappa*h)*(aul-auf)/(auf*aul*aul));
            else
              // following expression should be positive
              lfc(i) = e*kappa-damc(i)*e*kappa-af*(1.0-auf*(2.0*aul-3.0*auf)/(b*b)-2.0*(aul-auf)/(b*b)*sqrt(a+b*damc(i)*h*kappa)+damc(i)*h*kappa/b);
          }
        }
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }  
  return;
}



/**
  The function computes damage parameter which is the result of the
  damage function for brittle type of damage function. Correction of dissipated
  energy is controled by the cde flag.

  @param ipp - integration point number
  @param y -  parameter of damage function=equivalent strain (generalized conjugated dimensionless thermodynamical force)
  @param e - actual value of Young modulus
  @param f - actual value of strength (tensile or compressive)
  @param uf - actual value of initial crack opening
  @param omegao - old value of damage
  
  @return Function returns value of the damage parameter omega.
  
  Created by Tomas Koudelka,
*/
double ortodam::brittle_damage(long ipp, double y, double e, double f, double uf, double omegao)
{
  double err,omega = 0.0, dtmp, tmp, h, indam, kappa;
  long i,ni, eid;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  indam = f/e;
  
  kappa = y;
  if (kappa < indam)
    // elastic loading
    return 0.0;
  if ((Mm->ip[ipp].hmt & 2) > 0)
  {
    omega = 1.0-(f/(e*kappa))*exp(-(kappa-f/e)/(uf-f/e));
    return omega;
  }
  // correction of dispiated energy - mesh independent scalar damage
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
  // no correction of disipated energy only simple scalar damge
  { 
    omega = 1.0-(f/(e*kappa))*exp(-(kappa-f/e)/(uf-f/e));
    return omega;
  }
  // it is better to start from omega=1 in order to avoid of convergention problems
  omega = 1.0;
  for (i = 0; i < ni; i++)
  {
    dtmp = -e*kappa+f/uf*h*kappa*exp(-h*omega*kappa/uf);
    tmp = (1.0-omega)*e*kappa-f*exp(-h*omega*kappa/uf);
    if (fabs(tmp) < ft*err)
    {
      omega += - tmp/dtmp;
      break;
    }
    omega += - tmp/dtmp;
  }
  if ((omega > 1.0) || (omega < 0.0) || (fabs(tmp) > ft*err))
  // convergention problems have been encountered
  // starting from previously attained omega can help
  {
    omega = omegao;
    for (i = 0; i < ni; i++)
    {
      dtmp = -e*kappa+ft/uf*h*kappa*exp(-h*omega*kappa/uf);
      tmp = (1.0-omega)*e*kappa-ft*exp(-h*omega*kappa/uf);
      if (fabs(tmp) < ft*err)
      {
        omega += - tmp/dtmp;
        break;
      }
      omega += - tmp/dtmp;
    }
  }

  if (fabs(tmp) > ft*err)
  {
    print_err("iteration of omega doesn't converge with required error\n"
              " elem %ld, ip %ld:\n"
              " kappa=%le, ft=%le, e=%le, uf=%le, h=%le",
              __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1, kappa, ft, e, uf, h);
    abort();
  }
  
  if ((omega > 1.0) || (omega < 0.0))
  {
    print_err("computed omega=%le is out of prescribed range <0,1>\n"
              " elem %ld, ip %ld:\n"
              " kappa=%le, ft=%le, e=%le, uf=%le, h=%le",
              __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1, kappa, ft, e, uf, h);
    abort();
  }
  
  return omega;
}



/**
  The function computes damage parameter which is the result of the
  damage function for Bezier type of damage function. Correction of dissipated
  energy is controled by the cde flag.

  @param ipp - integration point number
  @param y -  parameter of damage function=equivalent strain (generalized conjugated dimensionless thermodynamical force)
  @param e - actual value of Young modulus
  @param f - actual value of strength (tensile or compressive)
  @param uf - actual value of initial crack opening
  @param ul - actual value of limit crack opening
  @param omegao - old value of damage
  @param pq - indicator of singularity in quadratic Beziere function (pq=1), pure quadratic function have to be used instead of its.
  
  @return Function returns value of the damage parameter omega.
  
  Created by Tomas Koudelka,
*/
double ortodam::qbezier_damage(long ipp, double y, double e, double f, double uf, double ul,double omegao,long pq)
{
  double err,omega = 0.0, dtmp, tmp, h, indam, kappa, a, b;
  long i,ni, eid;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  indam = f/e;
  
  kappa = y;
  if (kappa-indam <= Mp->zero)
    // elastic loading
    return 0.0;  
  a = uf*uf;
  b = (ul-2.0*uf);
  if ((Mm->ip[ipp].hmt & 2) > 0)
  {
    if (kappa >= ul) // limit deformation was reached => full damage, omega = 1.0
      return 1.0;
    omega = 1.0-2.0*uf*(ul-uf)/(b*b)-2.0*(ul-uf)*sqrt(a+b*(kappa-f/e))/(b*b)+(kappa-f/e)/b;
    omega = 1.0-(f/(e*(kappa-f/e)))*omega;
    return omega;
  }
  // correction of dispiated energy - mesh independent scalar damage
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
  // no correction of disipated energy only simple scalar damge
  { 
    if (kappa >= ul) // limit deformation was reached => full damage, omega = 1.0 
      return 1.0;
    if (pq) // pure quadratic function have to be used due to singularity
      omega = (1.0-omega)*e*kappa-f*(1.0-(kappa-f/e)/uf+(ul-uf)/(uf*ul*ul)*(kappa-f/e)*(kappa-f/e));
    else  // quadratic Beziere function
    {
      omega = 1.0-2.0*uf*(ul-uf)/(b*b)-2.0*(ul-uf)*sqrt(a+b*(kappa-f/e))/(b*b)+(kappa-f/e)/b;
      omega = 1.0-(f/(e*(kappa-f/e)))*omega;
    }
    return omega;
  }
  if (omega < 0.5)
    omega = 0.5;
  else
    omega = omegao;
  if (omega < 0.0)
    print_err("previous iteration of omega returned value < 0", __FILE__, __LINE__, __func__);

  if (kappa*omegao*h >= ul) // limit crack opening was reached => full damage, omega = 1.0 
    return 1.0;
  for (i = 0; i < ni; i++)
  {
    if (pq)
    {
      dtmp = 2.0*f*(uf-ul)*sqr(kappa*h)*omega/(uf*ul*ul)+kappa*h*f/uf-e*kappa;
      tmp  = (1.0-omega)*e*kappa-f*(1.0-omega*kappa*h/uf+(ul-uf)*sqr(omega*kappa*h)/(uf*ul*ul));
    }
    else
    {      
      dtmp = f*(ul-uf)*kappa*h/(b*sqrt(a+omega*b*kappa*h))-(e+f*h/b)*kappa;
      tmp = (1.0-omega)*e*kappa-f*(1.0+(2.0*uf*(ul-uf))/(b*b)-(2.0*(ul-uf)*sqrt(a+b*omega*h*kappa))/(b*b)+omega*h*kappa/b);
    }
    omega += - tmp/dtmp;
    if (omega < 0.0)
    {
      print_err("during iteration, omega reaches value < 0", __FILE__, __LINE__, __func__);
      abort();
    }
    if (fabs(tmp) <= ft*err)
      break;
  }
  if (fabs(tmp) > ft*err)
    print_err("iteration of omega doesn't converge with required error", __FILE__, __LINE__, __func__);
  
  if (omega > 1.0)
  { 
    if ((omega-1.0) > err)
      print_err("iteration of returns value > 1", __FILE__, __LINE__, __func__);
    else
      omega = 1.0;
  }
  if (omega < 0.0)
  {
    if (fabs(omega) > err)
      print_err("iteration of omega returns value < 0", __FILE__, __LINE__, __func__);
    else
      omega = 0.0;
  }

  
  return omega;
}



/**
  Function computes values of principal damage parameters Dt and Dc. The damage parameters
  are derived from the damage function.

  @param ipp  - integration point number
  @param peps - %vector of principal strains
  @param aat  - actual value of material parameter at 
  @param aac  - actual value of material parameter ac 
  @param pdamc - %vector of principal damage parameters for tension (output)
  @param pdamt - %vector of principal damage parameters for compression (output)
  
  @return The function returns principal damage parameters in the parameters pdamt and pdamc.

  Created by Tomas Koudelka, 9.2006
*/
void ortodam::princ_dam(long ipp, vector &peps, double aat, double aac, vector &pdamt, vector &pdamc)
{
  long i;
  double tmp, ltmp;
  double e, nu, af, auf, aul;
  double kappa;


  switch(damevf)
  {
    case brittle:
      e = Mm->give_actual_ym(ipp);
      nu = Mm->give_actual_nu(ipp);
      for(i=0; i<3; i++)
      {           
        kappa = dam_eq_strain(nu, i, peps);
        // tension 
        if (peps(i) >= 0.0)
        {
          // actual tensile strength
          af = Mm->give_actual_ft(ipp);
          // actual value of initial crack opening
          auf = uft*(af/ft);
          pdamt(i) = brittle_damage(ipp, kappa, e, af, auf, pdamt(i));
        }
        // compression 
        if (peps(i) < 0.0)
        {
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc*(af/fc);
          pdamc(i) = brittle_damage(ipp, kappa, e, af, auf, pdamc(i));
        }
      }
      break;
    case quasi_brittle:
      nu = Mm->give_actual_nu(ipp);
      for(i=0; i<3; i++)
      {    
        kappa = dam_eq_strain(nu, i, peps);
        // tension 
        if (peps(i) >= 0.0)
        {
          if (kappa > y0t)
          {    
            ltmp = log(kappa - y0t);
            pdamt(i) = aat*exp(ltmp*bt);
            tmp = 1.0/(1.0+pdamt(i));
            pdamt(i) = 1.0 - tmp;
          }
          if (pdamt(i) < 0.0)
          {
            fprintf(stdout, "\nDamage < 0.0 (ipp=%ld)\n", ipp);
            abort();
          }
        }
        // compression 
        if (peps(i) < 0.0)
        {
          if (kappa > y0c)
          {    
            ltmp = log(kappa - y0c);
            pdamc(i) = aac*exp(ltmp*bc);
            tmp = 1.0/(1.0+pdamc(i));
            pdamc(i) = 1.0-tmp;
          }
          if (pdamc(i) < 0.0)
          {
            fprintf(stdout, "\nDamage < 0.0 (ipp=%ld)\n", ipp);
            abort();
          }
        }
      }
      break;
    case quadbezier:
      e = Mm->give_actual_ym(ipp);
      nu = Mm->give_actual_nu(ipp);
      for(i=0; i<3; i++)
      {           
        kappa = dam_eq_strain(nu, i, peps);
        // tension 
        if (peps(i) >= 0.0)
          {
            // actual tensile strength
            af = Mm->give_actual_ft(ipp);
            // actual value of initial crack opening
            auf = uft*(af/ft);
            // actual value of limit crack opening
            aul = ult*(af/ft);
            pdamt(i) = qbezier_damage(ipp, kappa, e, af, auf, aul, pdamt(i),pqt);
          }
        // compression 
        if (peps(i) < 0.0)
        {
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc*(af/fc);
          // actual value of limit crack opening
          aul = ulc*(af/fc);
          pdamc(i) = qbezier_damage(ipp, kappa, e, af, auf, aul, pdamc(i),pqc);
        }
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
}



/**
  Function computes material parameters A, At and Ac. 
  If the correction of dissipated energy is required, they are computed with respect 
  to size of element for given ipp and variable softening modulus technique was used in 
  order to avoid mesh dependence. Otherwise, values read from the input file are used.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  @param aat - actual value of parameter At 
  @param aac - actual value of parameter Ac 

  @return The function returns actual values of material parameters A, A_t and A_c in the 
          parameters aa, aat, aac.
  
  Created by Tomas Koudelka, 9.2006
*/
void ortodam::give_actual_param_a(long ipp, long ido, double &aat, double &aac)
{
  long eid;
  double e, h, tmp;

  aat = Mm->ip[ipp].eqother[ido];
  aac = Mm->ip[ipp].eqother[ido+1];
  if ((aac == 0.0) || (aat == 0.0))
  {
    if (cde == corr_off)
    // no correction of dissipated energy is required
    {
      aat = at;
      aac = ac;
      Mm->ip[ipp].eqother[ido] = aat;
      Mm->ip[ipp].eqother[ido+1] = aac;
      return;
    }
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
        print_err("cannot determine dimension of the element.", __FILE__, __LINE__, __func__);
    }
    e = Mm->give_actual_ym(ipp);
    // variable softening for tension
    tmp = (gft/h - y0t/e)*bt*e*sin(M_PI/bt)/M_PI; 
    if (tmp < 0.0)
      print_err("cannot express variable softening for tensile damage.", __FILE__, __LINE__, __func__);
    aat = pow(tmp, -bt);
    // variable softening for compression
    tmp = (gfc/h - y0c/e)*bc*e*sin(M_PI/bc)/M_PI; 
    if (tmp < 0.0)
      print_err("cannot express variable softening for compressive damage.", __FILE__, __LINE__, __func__);
    aac = pow(tmp, -bc);
    Mm->ip[ipp].eqother[ido] = aat;
    Mm->ip[ipp].eqother[ido+1] = aac;
  }
}



/**
  The function returns the value of tensile strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns tensile strength. 

  10.2008 Created by Tomas Koudelka,
*/
double ortodam::give_actual_ft(long ipp)
{
  return ft;
}



/**
  The function returns the value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns compressive strength. 

  10.2008 Created by Tomas Koudelka,
*/
double ortodam::give_actual_fc(long ipp)
{
  return fc;
}



/**
  Function computes initializes material parameters A, At and Ac. 

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.
  
  Created by Tomas Koudelka, 9.2006
*/
void ortodam::initvalues(long ipp, long ido)
{
  double aat, aac;
  give_actual_param_a(ipp, ido, aac, aat);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka,
*/
void ortodam::matstiff (matrix &d,long ipp,long ido)
{
//  long ncompstr = Mm->ip[ipp].ncompstr;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      elmatstiff(d, ipp);
      break;
    case tangent_stiff:
      elmatstiff(d, ipp);
//      tmatstiff(d, ipp, ido);
      break;
    default:
      print_err("unknown type of stifness matrix is required.", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number

  @return The function returns elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka,
*/
void ortodam::elmatstiff (matrix &d,long ipp)
{
  double e, e0;
  long idem = Mm->ip[ipp].gemid();         // index of elastic material
  long idoem=Mm->givencompeqother(ipp,0);  // total number of internal variables
  long tmp=Mm->givencompeqother(ipp,idem); // number of internal variables for elastic material and thermal material 
  idoem -= tmp;

  Mm->elmatstiff (d,ipp);
  e = Mm->give_actual_ym(ipp);
  e0 = Mm->give_actual_ym(ipp,idem,idoem); // Young modulus from elastic material
  cmulm(e/e0, d);
}




/**
  The function computes material stiffnes %matrix used in principal stress computations.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param e - actual value of Young modulus
  @param nu - actual value of Poissons ratio
  
  @return The function returns elastic stiffness %matrix for principal stress computation 
          in the parameter d.

  Created by Tomas Koudelka,
*/
void ortodam::pelmatstiff (long ipp, matrix &d, double e, double nu)
{
  double b;

  switch (Mm->ip[ipp].ssst)
  {
    case planestress:
      b = e/(1.0-nu*nu);
      d[0][0] = b;     d[0][1] = b*nu;  d[0][2] = b*nu;
      d[1][0] = b*nu;  d[1][1] = b;     d[1][2] = b*nu;
      d[2][0] = b*nu;  d[2][1] = b*nu;  d[2][2] = b;
      break;
    case planestrain:
    case spacestress:
      b = e /((1.0 + nu)*(1.0 - 2.0*nu));
      d[0][0] = d[1][1] = d[2][2] = 1.0 - nu;
      d[0][1] = d[0][2] = d[1][0] = d[1][2] = d[2][0] = d[2][1] = nu;
      cmulm(b, d);
      break;
    default:
      print_err("Unknown type of stress/strain state is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2006
*/
void ortodam::nlstresses (long ipp, long im, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr;

  // actual material parameters a 
  double aat, aac;
  // old values of driving forces
  //double yo;
  vector pyto(3), pyco(3);
  // old values of damage parameters
  vector pdto(3), pdco(3);
  // damage parameters
  vector pdt(3), pdc(3);
  // strains and principal strains
  vector epsn(ncompstr);
  // tensor form of strains and transformation matrix for principal direction of strains
  matrix epst(3,3), t(3,3);
  // stress tensor and stress vector
  matrix sigt(3,3), d(3,3);
  vector sigma(ncompstr);
  // bulk modulus*3, shear modulus, volumetric strain, Young modulus, Poissons ratio
  double e, nu;
  // actual values of load fuctions
  //double lf;
  vector peps(3), psig(3), lft(3), lfc(3),pepsp(3);
  double h;
  long eid;

  // initialize values of material parameter a
  give_actual_param_a(ipp, ido, aat, aac);
  e  = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);

  //  total new strains
  for (i=0;i<ncompstr;i++)
    epsn[i] = Mm->ip[ipp].strain[i];

  // previous principal damage parameters
  for (i=0; i<3; i++)
  {
    // damage parameters
    pdt[i] = pdto[i] = Mm->ip[ipp].eqother[ido+2+i];   // tension
    pdc[i] = pdco[i] = Mm->ip[ipp].eqother[ido+2+3+i]; // compression    
  }
  
  // new driving forces
  vector_tensor(epsn, epst, Mm->ip[ipp].ssst, strain);
  princ_val (epst,peps,t,nijac,limit,Mp->zero,3, 1);
  
  // inelastic strains due to previous damage - fatigue model
  if (fat == fatigue_on)
  {
    for (i=0; i<3; i++)
    {
      pepsp[i]  = 0.5*betat*pdto[i]*pdto[i]/(1.0-pdto[i]);
      pepsp[i] -= 0.5*betac*pdco[i]*pdco[i]/(1.0-pdco[i]);
    }
    subv(peps,pepsp,peps);
  }

  princ_dam(ipp, peps, aat, aac, pdt, pdc);
  for (i=0; i<3; i++)
  {
    if (pdt[i] < pdto[i])
      pdt[i] = pdto[i];
    if (pdc[i] < pdco[i])
      pdc[i] = pdco[i];
  }

  // inelastic strains due to actual damage - fatigue model
  if (fat == fatigue_on)
  {
//    for (i=0; i<3; i++)
//    {
//      pepsp[i]  = 0.5*betat*pdt[i]*pdt[i]/(1.0-pdt[i]) - pepsp[i];
//      pepsp[i] -= 0.5*betac*pdc[i]*pdc[i]/(1.0-pdc[i]);
//    }
//    subv(peps,pepsp,peps);
  }
  // new principal stresses
  pelmatstiff (ipp, d, e, nu);
  mxv(d, peps, psig);
  for (i=0; i<3; i++)
  {
    if (psig(i) > 0.0)
      psig(i) *= 1.0-pdt(i);
    if (psig(i) < 0.0)
      psig(i) *= 1.0-pdc(i);
  }
   
  // transformation of principal stresses to the global coordinate system
  fillm(0.0, sigt);
  for (i = 0; i < 3; i++)
    sigt[i][i] = psig[i];
  lgmatrixtransf(sigt, t);
  tensor_vector(sigma, sigt, Mm->ip[ipp].ssst, stress);
  if (Mm->ip[ipp].ssst == planestress)
    sigma[3] = 0.0;

  //  new data storage
  for (i=0;i<ncompstr;i++)
    Mm->ip[ipp].stress[i]=sigma[i];
  for (i=0; i<3; i++)
  {
    // damage parameters
    Mm->ip[ipp].other[ido+2+i]   = pdt[i]; // tension
    Mm->ip[ipp].other[ido+2+3+i] = pdc[i]; // compression    
  }
  // following is only for debugging purposes
  // storage of tensile damage deviatoric part of tensor
  fillm(0.0, sigt);
  for (i = 0; i < 3; i++)
    sigt[i][i] = pdt[i];
  lgmatrixtransf(sigt, t);
  for (i = 0; i < 3; i++)
    Mm->ip[ipp].other[ido+8+i] = sigt[i][i];
  Mm->ip[ipp].other[ido+11] = sigt[1][2];       
  Mm->ip[ipp].other[ido+12] = sigt[0][2];       
  Mm->ip[ipp].other[ido+13] = sigt[0][1];       
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
  // storage of cracks openings computed from tensile damage tensor
  for (i = 0; i < 3; i++)
  {
    if (peps(i) > 0.0)
      Mm->ip[ipp].other[ido+14+i] = pdt[i]*(peps[i]-psig[i]/e)*h;
    else
      Mm->ip[ipp].other[ido+14+i] = 0.0;
  }

   
}

/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void ortodam::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}
