#include "anisodam.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "bmatrix.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
//#include "elastisomat.h"
#include "global.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>
#include <stdlib.h>

#define nijac 200
#define limit 1.0e-15
#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif





/**
  The constructor inializes attributes to zero values.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
anisodam::anisodam (void)
{
  cde = corr_off;
  ac = 0.0;  bc = 0.0;  y0c = 0.0;  at = 0.0;  bt = 0.0;  y0t = 0.0;
  a = 0.0;  b = 0.0;  y0 = 0.0; gfc = 0.0; gft = 0.0; gf = 0.0;
  ft = uft = 0.0;
  fc = ufc = 0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka, 9.2006
*/
anisodam::~anisodam (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void anisodam::read (XFILE *in)
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
          xfscanf (in, "%k%lf %k%lf %k%lf", "a", &a, "b", &b, "y0", &y0);
          break;
        case corr_on:
          xfscanf (in, "%k%lf %k%lf %k%lf", "G_fc", &gfc, "B_c", &bc, "Y_0c", &y0c);
          xfscanf (in, "%k%lf %k%lf %k%lf", "G_ft", &gft, "B_t", &bt, "Y_0t", &y0t);
          xfscanf (in, "%k%lf %k%lf %k%lf", "G_f", &gf, "b", &b, "y0", &y0);
          break;
        default:
          print_err(" unknown type of correction of disipated energy is required", __FILE__, __LINE__,__func__);
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
}



/** 
  Function computes damage driving force for volumetric damage parameter.

  @param  nu - Poissons ratio
  @param  peps - %vector of principal strains

  @return Function returns real value of damage driving force.

  Created by Tomas Koudelka, - 9.2006
*/
double anisodam::damdrvforce_vol(double nu, vector &peps)
{
  long i;
  double ev = 0.0, ret = 0.0;
  double k0_prime;  
  double g0_prime;  
  // strain tensor due to positive principal strains
  for (i = 0; i < 3; i++)
    ev += peps(i);
  ev /= 3.0;

  // take only positive volumetric strain  
  if (ev <= 0.0)
    return 0.0;
  // positive dimensionless volumetric damage driving force
  k0_prime = 1.0/(1.0-2.0*nu);
  g0_prime = 1.0/(2.0*(1.0+nu));
  ret = 3.0/2.0*(k0_prime-2.0*g0_prime)*ev*ev;  
  return ret;
}



/* 
  Function computes damage driving force for deviatoric damage parameters.

  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  t    - %matrix of principal directions stored in colomns for given
                 principal strains
  @param  pyt  - %vector of principal values of driving forces for tension damage parameters (output)
  @param  pyc  - %vector of principal values of driving forces for compression damage parameters  (output)       

  @return The function returns actual values of driving forces in parameters pyt and pyc.

  Created by Tomas Koudelka, - 9.2006
*/
/*
void anisodam::damdrvforce_dev(double nu, vector &peps, matrix &t, vector &pyt, vector &pyc)
{
  long i;
  matrix et(3,3), ec(3,3), pe(3,3), tyt(3,3), tyc(3,3);  
  double g0_prime;
  // strain tensor due to positive principal strains
  fillm(0.0, et);
  for (i = 0; i < 3; i++)
  {
    if (peps(i) > 0.0)
      et(i,i) = peps(i);
  }
  lgmatrixtransf(et, t);
  // strain tensor due to negative principal strains
  fillm(0.0, ec);
  for (i = 0; i < 3; i++)
  {
    if (peps(i) < 0.0)
      ec(i,i) = peps(i);
  }
  lgmatrixtransf(ec, t);
  // tensor of dimensionless damage driving force in tension 
  mxm(et, et, tyt);
  g0_prime = 1.0/(2.0*(1.0+nu));  
  cmulm(g0_prime, tyt);
  glmatrixtransf(tyt, t); // transformation to the principal directions
  // tensor of dimensionless damage driving force in compression 
  mxm(ec, ec, tyc);
  g0_prime = 1.0/(2.0*(1.0+nu));  
  cmulm(g0_prime, tyc);
  glmatrixtransf(tyc, t);  // transformation to the principal directions
  // transformation to the principal values vector
  for (i = 0; i < 3; i++)
  {   
    pyt(i) = tyt(i,i);
    pyc(i) = tyc(i,i);
  }  
}
*/



/** 
  Function computes damage driving force for deviatoric damage parameters.

  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  pyt  - %vector of principal values of driving forces for tension damage parameters
  @param  pyc  - %vector of principal values of driving forces for compression damage parameters        

  @return The function returns actual values of driving forces in parameters pyt and pyc.

  Recreated by Tomas Koudelka, - 6.2008
*/
void anisodam::damdrvforce_dev(double nu, vector &peps, vector &pyt, vector &pyc)
{
  long i;
  double g0_prime;

  // strain tensor due to positive principal strains
  g0_prime = 1.0/(2.0*(1.0+nu));  
  fillv(0.0, pyt);
  fillv(0.0, pyc);
  for (i = 0; i < 3; i++)
  {
    if (peps(i) > 0.0)
      pyt[i] = g0_prime*peps[i]*peps[i];
    if (peps(i) < 0.0)
      pyc[i] = g0_prime*peps[i]*peps[i];
  }
}



/**
  Function computes value of load function for the volumetric damage.

  @param  ipp  - integration point number
  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  d    - actual value of damage parameter
  @param  aa   - actual value of material parameter a

  @return The function returns value of load function for the volumetric damage.
  
  Created by Tomas Koudelka, 9.2006
*/
double anisodam::loadfuncvol(long ipp, double nu, vector &peps, double d, double aa)
{
  double y, f, kappa;
  double tmp1, tmp2;
  double e, aft, auft;

  y = damdrvforce_vol(nu, peps);
  switch(damevf)
  {
    case brittle:
      kappa = sqrt(y);
      if (kappa <= y0) 
       f = -1.0;      
      else
      {
        e = Mm->give_actual_ym(ipp);
        // actual tensile strength
        aft = Mm->give_actual_ft(ipp);
        // actual value of initial crack opening
        auft = uft/(aft/ft);
        tmp1 = -(kappa-aft/e)/(auft-aft/e);
        f = 1.0 - aft/(e*kappa)*exp(tmp1) - d;
      }
      break;
    case quasi_brittle:
      if (y <= y0) 
       f = -1.0;      
      else
      {
        // compute (y-y0)^b
        tmp1 = y - y0;
        tmp2 = log(tmp1);
        tmp1 = exp(b*tmp2);
        // compute value of load function
        f = (1.0 - d)*(1.0+aa*tmp1) - 1.0;
      }
      break;
    default:
      f = 0.0;
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
  return f;
}



/**
  Function computes value of load function for deviatoric damage

  @param  ipp  - integration point number
  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  damt - %vector actual values of principal damage parameters for tension
  @param  damc - %vector actual values of principal damage parameters for compression
  @param  aat  - actual values of material parameter a for tension 
  @param  aac  - actual values of material parameter a for compression 
  @param  lft  - %vector of resulting values of load function for tension for each principal direction
  @param  lfc  - %vector of resulting values of load function for compression for each principal direction

  @return The function returns load function values separately for tension and compression
          in parameters lft and lfc for each principal direction.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void anisodam::loadfuncdev(long ipp, double nu, vector &peps, vector &damt, vector &damc, 
                           double aat, double aac, vector &lft, vector &lfc)
{
  long i;
  double e, af, auf, kappa, tmp1, tmp2;
  vector pyt(3), pyc(3);

  damdrvforce_dev(nu, peps, pyt, pyc);
  switch(damevf)
  {
    case brittle:
      for (i=0; i<3; i++)
      {
        // for tension
        kappa = sqrt(pyt(i));
        if (kappa <= y0t) 
          lft(i) = -1.0;      
        else
        {
          e = Mm->give_actual_ym(ipp);
          // actual tensile strength
          af = Mm->give_actual_ft(ipp);
          // actual value of initial crack opening
          auf = uft/(af/ft);
          tmp1 = -(kappa-af/e)/(auf-af/e);
          lft(i) = 1.0 - af/(e*kappa)*exp(tmp1) - damt(i);
	}
        // for compression
        kappa = sqrt(pyc(i));
        if (kappa <= y0c) 
          lfc(i) = -1.0;      
        else
        {
          e = Mm->give_actual_ym(ipp);
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc/(af/fc);
          tmp1 = -(kappa-af/e)/(auf-af/e);
          lfc(i) = 1.0 - af/(e*kappa)*exp(tmp1) - damc(i);
	}
      }
      break;
    case quasi_brittle:
      for (i=0; i<3; i++)
      {
        // for tension
        if (pyt(i) <= y0t) 
          lft(i) = -1.0;      
        else
        {
          // compute (yt-y0t)^bt
          tmp1 = pyt(i) - y0t;
          tmp2 = log(tmp1);
          tmp1 = exp(bt*tmp2);
          // compute value of load function for tension in given principal direction i
          lft(i) = (1.0 - damt(i))*(1.0+aat*tmp1) - 1.0;
        }
        // for compression
        if (pyc(i) <= y0c) 
          lfc(i) = -1.0;      
        else
        {
          // compute (yc-y0c)^bc
          tmp1 = pyc(i) - y0c;
          tmp2 = log(tmp1);
          tmp1 = exp(bc*tmp2);
          // compute value of load function for compression in given principal direction i
          lfc(i) = (1.0 - damc(i))*(1.0+aac*tmp1) - 1.0;
        }
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }  
  return;
}



/**
  Function computes increments of damage parameter d. The increment of damage parameter
  is expressed as derivative of damage evolution function with respect to driving force.

  @param ipp - integration point number
  @param y   - driving forces for volumetric damage
  @param dy  - increment of driving forces for volumetric damage
  @param aa  - actual value of material parameter a 
  @param lf  - actual value of load function

  @return The function returns value of the damage parameter increment.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
double anisodam::daminc_vol(long ipp, double y, double dy, double aa, double lf)
{
  double e = Mm->give_actual_ym(ipp);
  double aft, auft, tmp, kappa, ddam = 0.0;

  switch(damevf)
  {
    case brittle:
      if ((lf > 0.0) && (dy > 0.0))
      {  
        kappa = sqrt(y);
        // actual tensile strength
        aft = Mm->give_actual_ft(ipp);
        // actual value of initial crack opening
        auft = uft/(aft/ft);
        tmp = -(kappa-aft/e)/(auft-aft/e);
        ddam = (aft/e)*(kappa-1.0)*exp(tmp)/(kappa*kappa);
      }
      break;
    case quasi_brittle:
      if ((lf > 0.0) && (dy > 0.0))
      {  
        tmp = log(y-y0);
        ddam = aa*b*exp(tmp*(b-1.0));
        tmp = 1.0+a*exp(tmp*b);
        ddam /= e*tmp*tmp;
        ddam *= dy;
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
  return ddam;
}



/**
  Function computes increments of damage parameters Dt and Dc. The increments of damage parameters
  are expressed as derivatives of damage evolution function with respect to driving forces.

  @param ipp - integration point number
  @param pyc - %vector of principal driving forces for tension
  @param pyt - %vector of principal driving forces for tension
  @param dyc - %vector of principal driving forces increments for tension
  @param dyt - %vector of principal driving forces increments for compression
  @param aat - actual value of material parameter at 
  @param aac - actual value of material parameter ac 
  @param lft - %vector of actual load function values for tension
  @param lfc - %vector of actual load function values for compression
  @param pdamc - %vector of principal damage parameters increments for tension
  @param pdamt - %vector of principal damage parameters increments for compression
  
  @return The function returns value of principal damage parameters increments in the parameters dpdamt and dpdamc.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void anisodam::pdaminc_dev(long ipp, vector &pyc, vector &pyt, vector &dyt, vector &dyc, 
                           double aat, double aac, vector &lft, vector &lfc, vector &dpdamt, vector &dpdamc)
{
  long i;
  double e = Mm->give_actual_ym(ipp);
  double kappa, tmp, ltmp;
  double af, auf;

  fillv(0.0, dpdamt);
  fillv(0.0, dpdamc);

  switch(damevf)
  {
    case brittle:
      for(i=0; i<3; i++)
      {    
        // tension 
        kappa = sqrt(pyt(i));
        if ((dyt[i] > 0.0) && (lft[i] > 0.0))
        {
          // actual tensile strength
          af = Mm->give_actual_ft(ipp);
          // actual value of initial crack opening
          auf = uft/(af/ft);
          tmp = -(kappa-af/e)/(auf-af/e);
          dpdamt(i) = (af/e)*(kappa-1.0)*exp(tmp)/(kappa*kappa);
	}
        // compression
        kappa = sqrt(pyc(i));
        if ((dyc[i] > 0.0) && (lfc[i] > 0.0))
        {
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc/(af/fc);
          tmp = -(kappa-af/e)/(auf-af/e);
          dpdamc(i) = (af/e)*(kappa-1.0)*exp(tmp)/(kappa*kappa);
	}
      }
      break;
    case quasi_brittle:
      for(i=0; i<3; i++)
      {    
        // tension 
        if ((dyt[i] > 0.0) && (lft[i] > 0.0))
        {
          ltmp = log(pyt(i) - y0t);
          dpdamt(i) = aat*bt*exp(ltmp*(bt-1.0));
          tmp = 1.0+aat*exp(ltmp*bt);
          dpdamt(i) /= e*(1.0+tmp*tmp);
          dpdamt(i) *= dyt(i);
	}
        // compression 
        if ((dyc[i] > 0.0) && (lfc[i] > 0.0))
        {
          ltmp = log(pyc(i) - y0c);
          dpdamc(i) = aac*bc*exp(ltmp*(bc-1.0));
          tmp = 1.0+aac*exp(ltmp*bc);
          dpdamc(i) /= e*(1.0+tmp*tmp);
          dpdamc(i) *= dyc(i);
	}
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
}



/**
  The function computes damage parameter which is the result of the
  damage function for brittle type of damage function. Correction of dissipated
  energy is controled by the cde flag.

  @param ipp - integration point number
  @param y -  parameter of damage function=equvalent strain (generalized conjugated thermodynamical force)
  @param e - actual value of Young modulus
  @param f - actual value of strength (tensile or compressive)
  @param uf - actual value of initial crack opening
  @param omegao - old value of damage
  
  @return Function returns value of the damage parameter omega.
  
  Created by Tomas Koudelka, 10.2008
*/
double anisodam::brittle_damage(long ipp, double y, double e, double f, double uf, double omegao)
{
  double err,omega = 0.0, dtmp, tmp, h=0.0, indam, kappa;
  long i,ni, eid;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  indam = f/e;
  
  kappa = sqrt(y);
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
  tmp = 0.0;
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
  Function computes value of damage parameter d. The damage parameter
  is derived from the damage function.

  @param ipp - integration point number
  @param y   - driving forces for volumetric damage
  @param aa  - actual value of material parameter a 
  @param dvo - old value of volumteric damage

  @return The function returns value of the volumetric damage parameter.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
double anisodam::dam_vol(long ipp, double y, double aa, double dvo)
{
  double e, tmp, dam = 0.0;
  double aft, auft;

  switch(damevf)
  {
    case brittle:
      e = Mm->give_actual_ym(ipp);
      // actual tensile strength
      aft = Mm->give_actual_ft(ipp);
      // actual value of initial crack opening
      auft = uft/(aft/ft);
      dam = brittle_damage(ipp, y, e, aft, auft, dvo);
      break;
    case quasi_brittle:
      if (y < y0)
        break;
      tmp = log(y-y0);
      dam = aa*exp(tmp*b);
      tmp = 1.0/(1.0+dam);
      dam = 1.0-tmp;
      if (dam < 0.0)
      {
        fprintf(stdout, "\nDamage < 0.0 (ipp=%ld)\n", ipp);
        abort();
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
  return dam;
}



/**
  Function computes values of principal damage parameters Dt and Dc. The damage parameters
  are derived from the damage function. Prameters pdamc and pdamt should contain
  old values of damage parameters attained in the previous steps.

  @param ipp - integration point number
  @param pyc - %vector of principal driving forces for tension
  @param pyt - %vector of principal driving forces for tension
  @param aat - actual value of material parameter at 
  @param aac - actual value of material parameter ac 
  @param pdamc - %vector of principal damage parameters for tension (input/output)
  @param pdamt - %vector of principal damage parameters for compression (input/output)
  
  @return The function returns values of principal damage parameters in the parameters pdamt and pdamc.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void anisodam::pdam_dev(long ipp, vector &pyc, vector &pyt, double aat, double aac, vector &pdamt, vector &pdamc)
{
  long i;
  double tmp, ltmp;
  double e, af, auf;

  fillv(0.0, pdamt);
  fillv(0.0, pdamc);

  switch(damevf)
  {
    case brittle:
      e = Mm->give_actual_ym(ipp);
      for(i=0; i<3; i++)
      {    
        // tension 
        // actual tensile strength
        af = Mm->give_actual_ft(ipp);
        // actual value of initial crack opening
        auf = uft/(af/ft);
        pdamt(i) = brittle_damage(ipp, pyt(i), e, af, auf, pdamt(i));
        // compression 
        // actual compressive strength
        af = Mm->give_actual_fc(ipp);
        // actual value of initial crack opening
        auf = ufc/(af/fc);
        pdamc(i) = brittle_damage(ipp, pyc(i), e, af, auf, pdamc(i));
      }
      break;
    case quasi_brittle:
      for(i=0; i<3; i++)
      {    
        // tension 
        if (pyt(i) > y0t)
        {    
          ltmp = log(pyt(i) - y0t);
          pdamt(i) = aat*exp(ltmp*bt);
          tmp = 1.0/(1.0+pdamt(i));
          pdamt(i) = 1.0 - tmp;
        }
        if (pdamt(i) < 0.0)
        {
          fprintf(stdout, "\nDamage < 0.0 (ipp=%ld)\n", ipp);
          abort();
        }
        // compression 
        if (pyc(i) > y0c)
        {    
          ltmp = log(pyc(i) - y0c);
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
  @param aa -  actual value of parameter A  (output)
  @param aat - actual value of parameter At (output) 
  @param aac - actual value of parameter Ac (output)

  @retval The function returns actual values of material parameters A, A_t and A_c in the
          parameters aa, aat, aac.
  
  Created by Tomas Koudelka, 9.2006
*/
void anisodam::give_actual_param_a(long ipp, long ido, double &aa, double &aat, double &aac)
{
  long eid;
  double e, h, tmp;

  aa  = Mm->ip[ipp].eqother[ido];
  aat = Mm->ip[ipp].eqother[ido+1];
  aac = Mm->ip[ipp].eqother[ido+2];
  if ((aa == 0.0) || (aac == 0.0) || (aat == 0.0))
  {
    if (cde == corr_off)
    // no correction of dissipated energy is required
    {
      aa  = a;
      aat = at;
      aac = ac;
      Mm->ip[ipp].eqother[ido]   = aa;
      Mm->ip[ipp].eqother[ido+1] = aat;
      Mm->ip[ipp].eqother[ido+2] = aac;
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
        h = 0.0;
        print_err("cannot determine dimension of element", __FILE__, __LINE__, __func__);
    }
    e = Mm->give_actual_ym(ipp);
    // variable softening for volumetric damage
    tmp = (gf/h - y0/e)*b*e*sin(M_PI/b)/M_PI; 
    if (tmp < 0.0)
      print_err("cannot express variable softening for volumetric damage", __FILE__, __LINE__, __func__);
    aa = pow(tmp, -b);
    // variable softening for tension
    tmp = (gft/h - y0t/e)*bt*e*sin(M_PI/bt)/M_PI; 
    if (tmp < 0.0)
      print_err("cannot express variable softening for tensile damage", __FILE__, __LINE__, __func__);
    aat = pow(tmp, -bt);
    // variable softening for compression
    tmp = (gfc/h - y0c/e)*bc*e*sin(M_PI/bc)/M_PI; 
    if (tmp < 0.0)
      print_err("cannot express variable softening for compressive damage", __FILE__, __LINE__, __func__);
    aac = pow(tmp, -bc);
    Mm->ip[ipp].eqother[ido]   = aa;
    Mm->ip[ipp].eqother[ido+1] = aat;
    Mm->ip[ipp].eqother[ido+2] = aac;
  }
}



/**
  The function returns the actual value of tensile strength.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns tensile strength. 

  Created by Tomas Koudelka, 10.2008
*/
double anisodam::give_actual_ft(long /*ipp*/)
{
  return ft;
}



/**
  The function returns the actual value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns compressive strength. 

  Created by Tomas Koudelka, 10.2008
*/
double anisodam::give_actual_fc(long /*ipp*/)
{
  return fc;
}



/**
  Function initializes material parameters A, At and Ac. 

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2006
*/
void anisodam::initvalues(long ipp, long ido)
{
  double aa, aat, aac;
  give_actual_param_a(ipp, ido, aa, aac, aat);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.
  
  Created by Tomas Koudelka,  9.2006
*/
void anisodam::matstiff (matrix &d,long ipp,long ido)
{
//  long ncompstr = Mm->ip[ipp].ncompstr;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff:
//      elmatstiff(d, ipp);
      tmatstiff(d, ipp, ido);
      break;
    default:
      print_err("unknown type of stiffness matrix is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 9.2006
*/
void anisodam::elmatstiff (matrix &d, long ipp, long ido)
{
  double e, e0;
  long idem = Mm->ip[ipp].gemid();

  if (Mm->ip[ipp].tm[idem] != elisomat)
    print_err("invalid type of elastic material is required", __FILE__, __LINE__, __func__);

  Mm->elmatstiff (d, ipp, ido);
  e  = Mm->give_actual_ym(ipp);
  e0 = Mm->give_initial_ym(ipp, idem, ido);
  cmulm(e/e0, d);
}




/**
  The function computes tangent material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  @return The function returns tangent stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 9.2006
*/
void anisodam::tmatstiff (matrix &d,long ipp, long ido)
{
  double e, nu, g, k0, dv, tmpt, tmpc;
  long i, j, k, l, m, idem, ncompstr = Mm->ip[ipp].ncompstr;
  vector epsn(ncompstr), peps(3), pdt(3), pdc(3);
  matrix epst(3,3), t(3,3);
  matrix dt(3,3), dc(3,3), omega_t(3,3), omega_c(3,3), dpdi_deps(3,3), dpdj_deps(3,3);
  bmatrix d_c(3,3), d_t(3,3), tmp_dc(3,3), tmp_dt(3,3);
  matrix d_tot(6,6), c_tot(6,6), c_red(3,3);
  double epsv;

  idem = Mm->ip[ipp].gemid();
  if (Mm->ip[ipp].tm[idem] != elisomat)
    print_err("invalid type of elastic material is required", __FILE__, __LINE__, __func__);

  e  = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);

  g = e/(2.0*(1.0+nu));
  k0 = e/(1.0-2.0*nu);
  dv = Mm->ip[ipp].other[ido+3]; // volumetric damage parameter
  for (i=0; i<3; i++)
  {
    pdt[i] = Mm->ip[ipp].eqother[ido+3+1+i];   // deviatoric - tension
    pdc[i] = Mm->ip[ipp].eqother[ido+3+1+3+i]; // deviatoric - compression    
  }
  //  total strains and their principal values and directions
  for (i=0;i<ncompstr;i++)
    epsn[i] = Mm->ip[ipp].strain[i];
  vector_tensor(epsn, epst, Mm->ip[ipp].ssst, strain);
  princ_val (epst,peps,t,nijac,limit,Mp->zero,3,1);
  epsv = peps[0] + peps[1] + peps[2];
  epsv /= 3.0;
  fillm(0.0, d);
  // volumetric part
  for(i=0; i<3; i++)
  {
    d_tot[i][i] = (1/3.0);
    if (epsv > 0.0)
      d_tot[i][i] -= dv*1/3.0;
    d_tot[i][i] *= (k0-2.0*g);
  }
  for(i=0; i<3; i++)
  {
    dt[i][i] = 1.0-pdt[i];
    dc[i][i] = 1.0-pdc[i];
  }  
  lgmatrixtransf(dt, t);
  lgmatrixtransf(dc, t);
  f_tensor(dt, &sqrt, t, omega_t);
  f_tensor(dc, &sqrt, t, omega_c);

  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      allocm(3, 3, d_c[i][j]);
      allocm(3, 3, d_t[i][j]);
      allocm(3, 3, tmp_dc[i][j]);
      allocm(3, 3, tmp_dt[i][j]);
    }
  }
  d_c.gen_indices();
  d_t.gen_indices();
  tmp_dc.gen_indices();
  tmp_dt.gen_indices();

  // computing of deviatoric parts
  for (m=0; m<3; m++) // loop over principal directions
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
        for(k=0; k<3; k++)
        {
          for(l=0; l<3; l++)
	  {
            if (dpdir_da(t, peps, m, i, dpdi_deps))
	    {
              elmatstiff(d, ipp, ido);
              return;
	    }
            dpdir_da(t, peps, m, j, dpdj_deps);
	    {
              elmatstiff(d, ipp, ido);
              return;
	    }
            if (peps[m] >= 0.0)
            {
              tmpt = peps[m]*(dpdj_deps[k][l]*t[i][m]+dpdi_deps[k][l]*t[j][m]);
              tmpt += t[i][m]*t[j][m]*t[k][m]*t[l][m];
	    }
            if (peps[m] < 0.0)
            {
              tmpc = peps[m]*(dpdj_deps[k][l]*t[i][m]+dpdi_deps[k][l]*t[j][m]);
              tmpc += t[i][m]*t[j][m]*t[k][m]*t[l][m];
	    }
            d_c[i][j][k][l] += 2.0*g*tmpc;
            d_t[i][j][k][l] += 2.0*g*tmpt;
	  }
	}
      }
    }
  }
  // computing of (1-d)^(1/2)*de/deps*(1-d)^(1/2)
  fillm(0.0, tmp_dc);
  fillm(0.0, tmp_dt);
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      for(k=0; k<3; k++)
      {
        for(l=0; l<3; l++)
	{
          for(m=0; m<3; m++)
	  {
            tmp_dc[i][j][k][l] += omega_c[i][m] * d_c[m][j][k][l];
            tmp_dt[i][j][k][l] += omega_t[i][m] * d_t[m][j][k][l];
	  }
	}
      }
    }
  }
  fillm(0.0, d_c);
  fillm(0.0, d_t);
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      for(k=0; k<3; k++)
      {
        for(l=0; l<3; l++)
	{
          for(m=0; m<3; m++)
	  {
            d_c[i][j][k][l] += omega_c[m][j] * tmp_dc[i][m][k][l];
            d_t[i][j][k][l] += omega_t[m][j] * tmp_dt[i][m][k][l];
	  }
	}
      }
    }
  } 
  // filling up the stiffness matrix
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
      d_tot[i][j] += d_c[i][i][j][j] + d_t[i][i][j][j];
  }
  for(i=0; i<3; i++)
  {
    d_tot[i][3] += d_c[i][i][1][2] + d_t[i][i][1][2];
    d_tot[i][4] += d_c[i][i][0][2] + d_t[i][i][0][2];
    d_tot[i][5] += d_c[i][i][0][1] + d_t[i][i][0][1];

    d_tot[3][i] += d_c[1][2][i][i] + d_t[1][2][i][i];
    d_tot[4][i] += d_c[0][2][i][i] + d_t[0][2][i][i];
    d_tot[5][i] += d_c[0][1][i][i] + d_t[0][1][i][i];
  }
  d_tot[3][3] += d_c[1][2][1][2] + d_t[1][2][1][2];
  d_tot[3][4] += d_c[1][2][0][2] + d_t[1][2][0][2];
  d_tot[3][5] += d_c[1][2][0][1] + d_t[1][2][0][1];

  d_tot[4][3] += d_c[0][2][1][2] + d_t[0][2][1][2];
  d_tot[4][4] += d_c[0][2][0][2] + d_t[0][2][0][2];
  d_tot[4][5] += d_c[0][2][0][1] + d_t[0][2][0][1];

  d_tot[5][3] += d_c[0][1][1][2] + d_t[0][1][1][2];
  d_tot[5][4] += d_c[0][1][0][2] + d_t[0][1][0][2];
  d_tot[5][5] += d_c[0][1][0][1] + d_t[0][1][0][1];

  switch (Mm->ip[ipp].ssst)
  {
    case planestress:
    {
      invm(d_tot, c_tot, 1e-8);
      for(i=0; i<2; i++)
      {
        for(j=0; j<2; j++)
          c_red[i][j] = c_tot[i][j];
      }
      c_red[0][2] = c_tot[0][5];
      c_red[1][2] = c_tot[1][5];
      c_red[2][2] = c_tot[5][5];

      c_red[2][0] = c_tot[5][0];
      c_red[2][1] = c_tot[5][1];
      invm(c_red, d, 1e-8);
      break;
    }
    case planestrain:
    {
      for(i=0; i<2; i++)
      {
        for(j=0; j<2; j++)
          d[i][j] = d_tot[i][j];
      }
      d[0][2] = d_tot[0][5];
      d[1][2] = d_tot[1][5];
      d[2][2] = d_tot[5][5];

      d[2][0] = d_tot[5][0];
      d[2][1] = d_tot[5][1];
      break;
    }
    case spacestress:
      copym(d_tot, d);
      break;
    default:
      print_err("Unknown type of stress/strain state is required", __FILE__, __LINE__, __func__);  
      abort();
  }
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
*/
void anisodam::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr;

  // actual material parameters a 
  double aa, aat, aac;
  // old values of driving forces
  //double yo;
  vector pyto(3), pyco(3);
  // old values of damage parameters
  double dvo;
  vector pdto(3), pdco(3);
  // increments of driving forces;
  //double dy; 
  vector pdyt(3), pdyc(3);
  // new values of driving forces
  double yn;
  vector pytn(3), pycn(3);
  // damage parameters
  double dv;
  vector pdt(3), pdc(3);
  // strains and principal strains
  vector epsn(ncompstr);
  // tensor form of strains and transformation matrix for principal direction of strains
  matrix epst(3,3), t(3,3);
  // stress tensor and stress vector
  matrix sigt(3,3);
  vector sigma(ncompstr);
  // bulk modulus*3, shear modulus, volumetric strain, Young modulus, Poissons ratio
  double k0, g0, epsv, e, nu;
  // actual values of load fuctions
  //double lf;
  vector peps(3), psig(3), lft(3), lfc(3),pepsp(3);
  double tmp, h;
//  double u_crack;
  long eid;
  // coefficients for inelastic strains evolution
  // double betat = 1.0e-4;
  // double betac = 1.0e-4;
  // double gamma = 1.0e-4;

  // initialize values of material parameter a
  give_actual_param_a(ipp, ido, aa, aat, aac);

  e  = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);

  //  total new strains
  for (i=0;i<ncompstr;i++)
    epsn[i] = Mm->ip[ipp].strain[i];

  // corresponding previous principal damage parameters
  dvo     = Mm->ip[ipp].eqother[ido+3]; // volumetric damage parameter
  for (i=0; i<3; i++)
  {
    // damage parameters
    pdto[i] = Mm->ip[ipp].eqother[ido+3+1+i];   // deviatoric - tension
    pdco[i] = Mm->ip[ipp].eqother[ido+3+1+3+i]; // deviatoric - compression    
  }
  
  // new driving forces
  vector_tensor(epsn, epst, Mm->ip[ipp].ssst, strain);
  princ_val (epst,peps,t,nijac,limit,Mp->zero,3, 1);
  
  /* for (i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
      Mm->ip[ipp].other[ido+3+1+3+3+1+3+3+3*i+j] = t[i][j]; // 
  }  
  // inelastic strains due to previous damage
  for (i=0; i<3; i++)
  {
    pepsp[i]  = betat*pdto[i]*pdto[i]/(1.0-pdto[i]);
    pepsp[i] -= betac*pdco[i]*pdco[i]/(1.0-pdco[i]);
    pepsp[i] -= gamma*dvo*dvo/(1.0-dvo);
  }
  subv(peps,pepsp,peps);*/

  //if (Mp->time > 190800.0){
  //printf("tady!!! \n\n");
  //}

  yn = damdrvforce_vol(nu, peps);
  damdrvforce_dev(nu, peps, pytn, pycn);

  // new damage parameters
  dv = dam_vol(ipp, yn, aa, dvo);
  if (dv < dvo)
    dv = dvo;
  copyv(pdto, pdt);
  copyv(pdco, pdc);
  pdam_dev(ipp, pycn, pytn, aat, aac, pdt, pdc);
  for (i=0; i<3; i++)
  {
    if (pdt[i] < pdto[i])
      pdt[i] = pdto[i];
    if (pdc[i] < pdco[i])
      pdc[i] = pdco[i];
  }

  // new principal stresses
  k0 = e;
  k0 /= (1-2.0*nu);
  g0 = e;
  g0 /= 2.0*(1.0+nu);
  epsv = peps[0]+peps[1]+peps[2];
  epsv /= 3.0;
  psig[0] = k0 - 2.0*g0;
  if (epsv > 0.0)
    epsv -= dv*epsv;
//    psig[0] -= k0*(dv);
  psig[0] *= epsv;
  psig[1] = psig[2] = psig[0];
  for (i=0; i<3; i++)
  {
    tmp = 1.0;
    if (peps[i] > 0.0)
      tmp -= pdt[i];
    if (peps[i] < 0.0)
      tmp -= pdc[i];
    psig[i] += 2.0*g0*tmp*peps[i];
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
  Mm->ip[ipp].other[ido+3] = dv; // volumetric damage parameter
  for (i=0; i<3; i++)
  {
    // damage parameters
    Mm->ip[ipp].other[ido+3+1+i]   = pdt[i]; // deviatoric - tension
    Mm->ip[ipp].other[ido+3+1+3+i] = pdc[i]; // deviatoric - compression    
  }
  // following is only for debugging purposes
  // storage of tensile damage deviatoric part of tensor
  fillm(0.0, sigt);
  for (i = 0; i < 3; i++)
    sigt[i][i] = pdt[i];
  lgmatrixtransf(sigt, t);
  for (i = 0; i < 3; i++)
    Mm->ip[ipp].other[ido+10+i] = sigt[i][i];
  Mm->ip[ipp].other[ido+13] = sigt[1][2];       
  Mm->ip[ipp].other[ido+14] = sigt[0][2];       
  Mm->ip[ipp].other[ido+15] = sigt[0][1];       
  // storage of tensile damage of tensor
  for (i = 0; i < 3; i++)
    Mm->ip[ipp].other[ido+16+i] = sigt[i][i]+dv;       
  Mm->ip[ipp].other[ido+19] = sigt[1][2];       
  Mm->ip[ipp].other[ido+20] = sigt[0][2];       
  Mm->ip[ipp].other[ido+21] = sigt[0][1];       
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
      h = 0.0;
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
  }
  // storage of crack opening in volumetric part
  Mm->ip[ipp].other[ido+22+i] = (sqrt(yn)-(psig[0]+psig[1]+psig[2])/(3.0*e))*h; 
  // storage of cracks openings computed from tensile damage tensor
  for (i = 0; i < 3; i++)
    Mm->ip[ipp].other[ido+23+i] = (sqrt(pytn[i])-psig[i]/e)*h;   
}



/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2006
*/
void anisodam::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}
