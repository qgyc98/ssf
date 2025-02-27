#include "anisodamrot.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>
#include <stdlib.h>

#define nijac 20
#define limit 1.0e-16
#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif




/**
  The constructor inializes attributes to zero values.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
anisodamrot::anisodamrot (void)
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
anisodamrot::~anisodamrot (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
void anisodamrot::read (XFILE *in)
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
  Function computes dimensionless damage driving force for volumetric damage parameter.

  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
    
  @return Function returns value of damage driving force.

  Created by Tomas Koudelka, 9.2006
*/
double anisodamrot::damdrvforce_vol(double nu, vector &peps)
{
  long i;
  double ev = 0.0;
  double k0_prime, g0_prime, ret;  
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



/** 
  Function computes dimensionless damage driving force for deviatoric damage parameters.

  @param  nu   - Poissons ratio
  @param  peps - %vector of principal strains
  @param  pyt  - %vector of principal values of driving forces for tension damage parameters
  @param  pyc  - %vector of principal values of driving forces for compression damage parameters        

  @return The function returns actual values of driving forces in parameters pyt and pyc.

  Created by  Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
void anisodamrot::damdrvforce_dev(double nu, vector &peps, vector &pyt, vector &pyc)
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
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
double anisodamrot::loadfuncvol(long ipp, double nu, vector &peps, double d, double aa)
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
  @param  lft  - %vector of resulting values of load function for tension for each principal direction (output)
  @param  lfc  - %vector of resulting values of load function for compression for each principal direction (output)

  @return The function returns load function values separately for tension and compression
          in parameters lft and lfc for each principal direction.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2008
  Modified by Tomas Koudelka, 10.2009
*/
void anisodamrot::loadfuncdev(long ipp, double nu, vector &peps, vector &damt, vector &damc, 
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
  Function computes value of damage parameter d. The damage parameter
  is derived from the damage function.

  @param ipp - integration point number
  @param y   - driving forces for volumetric damage
  @param aa  - actual value of material parameter a 
  @param dvo - old value of volumteric damage 

  @return The function returns value of the volumetric damage parameter.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
double anisodamrot::dam_vol(long ipp, double y, double aa, double dvo)
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
  are derived from the damage function.

  @param ipp - integration point number
  @param pyc - %vector of principal driving forces for tension
  @param pyt - %vector of principal driving forces for tension
  @param aat - actual value of material parameter at 
  @param aac - actual value of material parameter ac 
  @param pdamc - %vector of principal damage parameters for tension (output)
  @param pdamt - %vector of principal damage parameters for compression (output)

  @return The function returns actual principal values of damage in the parameters pdamt and pdamc.
  
  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
void anisodamrot::pdam_dev(long ipp, vector &pyc, vector &pyt, double aat, double aac, vector &pdamt, vector &pdamc)
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
  
  Created by Tomas Koudelka, 10.2009
*/
double anisodamrot::brittle_damage(long ipp, double y, double e, double f, double uf, double omegao)
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

  @retval The function returns computed material parameters A, Ac and At in the parameters aa, aat and aac.
  
  Created by Tomas Koudelka, 9.2006
*/
void anisodamrot::give_actual_param_a(long ipp, long ido, double &aa, double &aat, double &aac)
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
  
  @return The function returns tensile stregth. 

  Created by Tomas Koudelka, 10.2009
*/
double anisodamrot::give_actual_ft(long /*ipp*/)
{
  return ft;
}



/**
  The function returns the value of compressive strength.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns compressive stregth. 

  Created by Tomas Koudelka, 10.2009
*/
double anisodamrot::give_actual_fc(long /*ipp*/)
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
void anisodamrot::initvalues(long ipp, long ido)
{
  double aa, aat, aac;
  give_actual_param_a(ipp, ido, aa, aac, aat);
}



/**
  This function computes actual material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka,
*/
void anisodamrot::matstiff (matrix &d,long ipp,long ido)
{
//  long ncompstr = Mm->ip[ipp].ncompstr;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff:
      elmatstiff(d, ipp, ido);
//      tmatstiff(d, ipp, ido);
      break;
    default:
      print_err("unknown type of stiffness matrix is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  @return The function returns elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 9.2006
*/
void anisodamrot::elmatstiff (matrix &d, long ipp, long ido)
{
  double e, e0;
  long idem = Mm->ip[ipp].gemid();

  if (Mm->ip[ipp].tm[idem] != elisomat)
    print_err("invalid type of elastic material is required", __FILE__, __LINE__, __func__);

  Mm->elmatstiff (d, ipp, ido);
  e = Mm->give_actual_ym(ipp);
  e0 = Mm->give_initial_ym(ipp);
  cmulm(e/e0, d);
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  Parameters:
  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created  by Tomas Koudelka,  9.2006
  Modified by Tomas Koudelka, 10.2009
*/
void anisodamrot::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr;

  // actual material parameters a 
  double aa, aat, aac;
  // previous values of volumetric damage parameter
  double dvo;
  // actual value of volumetric damage parameter
  double dv;
  // tensile and compressive damage tensors 
  matrix dt(3,3), dc(3,3);
  // decomposed strains to tensional and compressive parts
  matrix eps_t(3,3), eps_c(3,3);
  // new values of driving forces
  double yn;
  vector pyt(3), pyc(3);
  // actual principal values of damage tensors
  vector pdt(3), pdc(3);
  // actual strains and principal strains
  vector epsn(ncompstr);
  vector peps(3);
  // tensor form of strains and transformation matrix for principal direction of strains
  matrix epst(3,3), t(3,3);
  // auxiliary matrices for damage
  matrix omegat(3,3), omegac(3,3);
  // auxiliary matrices
  matrix tmp1(3,3), tmp2(3,3);
  // auxiliary transformation matrix for damage
  matrix tomega(3,3);
  // stress tensor and stress vector
  matrix sigt(3,3);
  vector sigma(ncompstr);
  // bulk modulus*3, shear modulus, volumetric strain, Young modulus, Poissons ratio
  double k0, g0, epsv, e, nu;

  e = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  // initialize values of material parameter a
  give_actual_param_a(ipp, ido, aa, aat, aac);

  //  total new strains
  for (i=0;i<ncompstr;i++)
    epsn[i] = Mm->ip[ipp].strain[i];

  // previous damage parameters
  dvo = Mm->ip[ipp].eqother[ido+3]; // volumetric damage parameter
  // deviatoric damage parameters for tension
  dt[0][0] = Mm->ip[ipp].eqother[ido+3+1+0];
  dt[1][1] = Mm->ip[ipp].eqother[ido+3+1+1];
  dt[2][2] = Mm->ip[ipp].eqother[ido+3+1+2];
  dt[1][2] = dt[2][1] = Mm->ip[ipp].eqother[ido+3+1+3];
  dt[0][2] = dt[2][0] = Mm->ip[ipp].eqother[ido+3+1+4];
  dt[0][1] = dt[1][0] = Mm->ip[ipp].eqother[ido+3+1+5];
  // deviatoric damage parameters for compression
  dc[0][0] = Mm->ip[ipp].eqother[ido+3+1+6+0];
  dc[1][1] = Mm->ip[ipp].eqother[ido+3+1+6+1];
  dc[2][2] = Mm->ip[ipp].eqother[ido+3+1+6+2];
  dc[1][2] = dc[2][1] = Mm->ip[ipp].eqother[ido+3+1+6+3];
  dc[0][2] = dc[2][0] = Mm->ip[ipp].eqother[ido+3+1+6+4];
  dc[0][1] = dc[1][0] = Mm->ip[ipp].eqother[ido+3+1+6+5];
  
  // new driving forces
  vector_tensor(epsn, epst, Mm->ip[ipp].ssst, strain);
  princ_val (epst,peps,t,nijac,limit,Mp->zero,3,1);
  yn = damdrvforce_vol(nu, peps);
  damdrvforce_dev(nu, peps, pyt, pyc);

  // transformation of attained damage form global coordinate system to new principal directions
  glmatrixtransf(dt, t);
  glmatrixtransf(dc, t);

  // new damage parameters
  dv = dam_vol(ipp, yn, aa, dvo);
  // initializing principal damage with old values
  for (i=0; i<3; i++)
  {
    pdt[i] = dt[i][i];
    pdc[i] = dc[i][i];
  }
  pdam_dev(ipp, pyc, pyt, aat, aac, pdt, pdc);

  // damage parameters can be increased only
  for (i = 0; i < 3; i++)
  {
    if (dt[i][i] < pdt[i])
      dt[i][i] = pdt[i];
    if (dc[i][i] < pdc[i])
      dc[i][i] = pdc[i];
  }
  // transformation of principal damage to the global coordinate system
  lgmatrixtransf(dt, t);
  lgmatrixtransf(dc, t);

  // damage parameters can be increased only
  if (dv < dvo)
    dv = dvo;    

  // compute Omega_t=(1-Dt)^(1/2)
  princ_val(dt, pdt, tomega, nijac, limit, Mp->zero, 3, 1);
  fillm(0.0, omegat);
  for(i=0; i<3; i++)
  {
    if ((pdt[i] < 0.0) || (pdt[i] > 1.0))
    {
      print_err("principal value of damage tensor for tension is out of range <0,1>.", __FILE__, __LINE__, __func__);
      if (pdt[i] < 0.0)
        pdt[i] = 0.0;
      else
        pdt[i] = 1.0;
    }
    omegat[i][i] = sqrt(1.0-pdt[i]);
  }
  lgmatrixtransf(omegat, tomega);

  // compute Omega_c=(1-Dc)^(1/2)
  princ_val(dc, pdc, tomega, nijac, limit, Mp->zero, 3, 1);
  fillm(0.0, omegac);
  for(i=0; i<3; i++)
  {
    if ((pdc[i] < 0.0) || (pdc[i] > 1.0))
    {
      print_err("principal value of damage tensor for compression is out of range <0,1>.", __FILE__, __LINE__, __func__);
      if (pdc[i] < 0.0)
        pdc[i] = 0.0;
      else
        pdc[i] = 1.0;
    }
    omegac[i][i] = sqrt(1.0-pdc[i]);
  }
  lgmatrixtransf(omegac, tomega);
      
  // new stresses
  k0 = e;
  k0 /= (1-2.0*nu);
  g0 = e;
  g0 /= 2.0*(1.0+nu);
  // compute volumetric part: 3/2*(ko-2*Go)*(1-d)*eps_v
  epsv = peps[0]+peps[1]+peps[2];
  epsv /= 3.0;
  // strain tensor due to positive principal strains
  fillm(0.0, eps_t);
  for (i = 0; i < 3; i++)
  {
    if (peps[i] > 0.0)
      eps_t[i][i] = peps[i];
  }
  lgmatrixtransf(eps_t, t);
  // strain tensor due to negative principal strains
  fillm(0.0, eps_c);
  for (i = 0; i < 3; i++)
  {
    if (peps[i] < 0.0)
      eps_c[i][i] = peps[i];
  }
  lgmatrixtransf(eps_c, t);

  sigt[0][0] = k0 - 2.0*g0;
  if (epsv > 0.0)
    epsv -= dv*epsv;
  sigt[0][0] *= epsv;
  sigt[1][1] = sigt[2][2] = sigt[0][0];
  // compute deviatoric part for tension: 2*Go*Omega_t*eps_t*Omega_t
  mxm(omegat, eps_t, tmp1);
  mxm(tmp1, omegat, tmp2);
  cmulm(2.0*g0, tmp2);
  addm(sigt, tmp2, sigt);
  // compute deviatoric part for compression: 2*Go*Omega_c*eps_c*Omega_c
  mxm(omegac, eps_c, tmp1);
  mxm(tmp1, omegac, tmp2);
  cmulm(2.0*g0, tmp2);
  addm(sigt, tmp2, sigt);
  tensor_vector(sigma, sigt, Mm->ip[ipp].ssst, stress);
  for (i=0;i<ncompstr;i++)
    Mm->ip[ipp].stress[i]=sigma[i];

  //  new internal data storage
  Mm->ip[ipp].other[ido+3] = dv; // volumetric damage parameter
  // deviatoric damage parameters for tension
  Mm->ip[ipp].other[ido+3+1+0] = dt[0][0];
  Mm->ip[ipp].other[ido+3+1+1] = dt[1][1];
  Mm->ip[ipp].other[ido+3+1+2] = dt[2][2];
  Mm->ip[ipp].other[ido+3+1+3] = dt[1][2];
  Mm->ip[ipp].other[ido+3+1+4] = dt[0][2];
  Mm->ip[ipp].other[ido+3+1+5] = dt[0][1];
  // deviatoric damage parameters for compression
  Mm->ip[ipp].other[ido+3+1+6+0] = dc[0][0];
  Mm->ip[ipp].other[ido+3+1+6+1] = dc[1][1];
  Mm->ip[ipp].other[ido+3+1+6+2] = dc[2][2];
  Mm->ip[ipp].other[ido+3+1+6+3] = dc[1][2];
  Mm->ip[ipp].other[ido+3+1+6+4] = dc[0][2];
  Mm->ip[ipp].other[ido+3+1+6+5] = dc[0][1];
  // principal values of damage tensor for tension
  Mm->ip[ipp].other[ido+3+1+6+6+0] = pdt[0];
  Mm->ip[ipp].other[ido+3+1+6+6+1] = pdt[1];
  Mm->ip[ipp].other[ido+3+1+6+6+2] = pdt[2];
  // principal values of damage tensor for compression
  Mm->ip[ipp].other[ido+3+1+6+6+3+0] = pdc[0];
  Mm->ip[ipp].other[ido+3+1+6+6+3+1] = pdc[1];
  Mm->ip[ipp].other[ido+3+1+6+6+3+2] = pdc[2];
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
void anisodamrot::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}
