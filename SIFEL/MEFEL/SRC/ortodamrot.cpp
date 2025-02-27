#include "ortodamrot.h"
#include "bmatrix.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
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
  Modified by Tomas Koudelka, 10.2008
*/
ortodamrot::ortodamrot (void)
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
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka, 9.2006
*/
ortodamrot::~ortodamrot (void)
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
void ortodamrot::read (XFILE *in)
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
          print_err("unknown type of correction of disipated energy is required.", __FILE__, __LINE__,__func__);
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
      print_err("unknown type of damage evolution function is required.", __FILE__, __LINE__,__func__);
  }
  xfscanf (in,"%k%m", "fatigue", &fatigue_flag_kwdset, (int *)(&fat));
  if (fat == fatigue_on)
    xfscanf (in,"%k%lf %k%lf", "beta_t", &betat, "beta_c", &betac);
}



/**
  Function computes value of load function 

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
*/
void ortodamrot::loadfunc(long ipp, double /*nu*/, vector &peps, vector &damt, vector &damc, 
                       double aat, double aac, vector &lft, vector &lfc)
{
  long i;
  double e, af, auf, aul, kappa, tmp1, tmp2, a, b, h;
  vector pyt(3), pyc(3);

  switch(damevf)
  {
    case brittle:
      for (i=0; i<3; i++)
      {
        // for tension
        if (peps(i) > 0.0)
	{
          kappa = peps(i);
          if (kappa <= y0t) 
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
	  kappa = fabs(peps(i));
          if (kappa <= y0c) 
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
        if (peps(i) > 0.0)
	{
          kappa = peps(i);
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
          kappa = fabs(peps(i));
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
        if (peps(i) > 0.0)
	{
          kappa = peps(i);
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
	  kappa = fabs(peps(i));
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
      print_err("unknown type of damage evolution function is required.", __FILE__, __LINE__,__func__);
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
double ortodamrot::brittle_damage(long ipp, double y, double e, double f, double uf, double omegao)
{
  double err,omega = 0.0, dtmp, tmp, h, indam, kappa;
  long i,ni, eid;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  indam = f/e;
  
  kappa = y;
  if (kappa < indam)
    // elastic loading/unloading
    return omegao;
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
      print_err("unknown dimension of element is required.", __FILE__, __LINE__, __func__);
      abort();
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
double ortodamrot::qbezier_damage(long ipp, double y, double e, double f, double uf, double ul,double omegao,long pq)
{
  double err,omega = 0.0, dtmp, tmp, h, indam, kappa, a, b;
  long i,ni, eid;
  
  ni=sra.give_ni ();
  err=sra.give_err ();
  
  indam = f/e;
  
  kappa = y;
  if (kappa-indam <= Mp->zero)
    // elastic loading
    return omegao;  
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
      print_err("unknown dimension of element is required.", __FILE__, __LINE__, __func__);
      abort();
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
    print_err("previous iteration of omega returned value < 0.", __FILE__, __LINE__, __func__);

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
      print_err("during iteration, omega reaches value < 0.", __FILE__, __LINE__, __func__);
      abort();
    }
    if (fabs(tmp) <= ft*err)
      break;
  }
  if (fabs(tmp) > ft*err)
    print_err("iteration of omega doesn't converge with required error.", __FILE__, __LINE__, __func__);
  
  if (omega > 1.0)
  { 
    if ((omega-1.0) > err)
      print_err("iteration of omega returns value > 1.", __FILE__, __LINE__, __func__);
    else
      omega = 1.0;
  }
  if (omega < 0.0)
  {
    if (fabs(omega) > err)
      print_err("iteration of omega returns value < 0.", __FILE__, __LINE__, __func__);
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
void ortodamrot::princ_dam(long ipp, vector &peps, double aat, double aac, vector &pdamt, vector &pdamc)
{
  long i;
  double tmp, ltmp;
  double e, af, auf, aul;


  switch(damevf)
  {
    case brittle:
      e = Mm->give_actual_ym(ipp);
      for(i=0; i<3; i++)
      {           
        // tension 
        if (peps(i) > 0.0)
	{
          // actual tensile strength
          af = Mm->give_actual_ft(ipp);
          // actual value of initial crack opening
          auf = uft*(af/ft);
          pdamt(i) = brittle_damage(ipp, peps(i), e, af, auf, pdamt(i));
	}
        // compression 
        if (peps(i) < 0.0)
	{
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc*(af/fc);
          pdamc(i) = brittle_damage(ipp, fabs(peps(i)), e, af, auf, pdamc(i));
	}
      }
      break;
    case quasi_brittle:
      for(i=0; i<3; i++)
      {    
        // tension 
        if (peps(i) > 0.0)
	{
          if (peps(i) > y0t)
          {    
            ltmp = log(peps(i) - y0t);
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
          if (fabs(peps(i)) > y0c)
          {    
            ltmp = log(fabs(peps(i)) - y0c);
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
      for(i=0; i<3; i++)
      {           
        // tension 
        if (peps(i) > 0.0)
	  {
	    // actual tensile strength
	    af = Mm->give_actual_ft(ipp);
	    // actual value of initial crack opening
	    auf = uft*(af/ft);
	    // actual value of limit crack opening
	    aul = ult*(af/ft);
	    pdamt(i) = qbezier_damage(ipp, peps(i), e, af, auf, aul, pdamt(i),pqt);
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
          pdamc(i) = qbezier_damage(ipp, fabs(peps(i)), e, af, auf, aul, pdamc(i),pqc);
	}
      }
      break;
    default:
      print_err("unknown type of damage evolution function is required.", __FILE__, __LINE__,__func__);
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
void ortodamrot::give_actual_param_a(long ipp, long ido, double &aat, double &aac)
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
        print_err("cannot determine dimension of element.", __FILE__, __LINE__, __func__);
        abort();
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

  Created by Tomas Koudelka, 10.2008 
*/
double ortodamrot::give_actual_ft(long /*ipp*/)
{
  return ft;
}



/**
  The function returns the value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns compressive strength. 

  Created by Tomas Koudelka, 10.2008 
*/
double ortodamrot::give_actual_fc(long /*ipp*/)
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
void ortodamrot::initvalues(long ipp, long ido)
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
void ortodamrot::matstiff (matrix &d,long ipp,long /*ido*/)
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
void ortodamrot::elmatstiff (matrix &d,long ipp)
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
void ortodamrot::pelmatstiff (long ipp, matrix &d, double e, double nu)
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
      print_err("unknown type of stress/strain state is required.", __FILE__, __LINE__, __func__);
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
void ortodamrot::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr;

  // actual material parameters a 
  double aat, aac;
  // tensile and compressive damage tensors 
  matrix dt(3,3), dc(3,3);
  // stiffness matrix
  matrix d(ncompstr, ncompstr);
  // decomposed strains to tensional and compressive parts
  matrix eps_t(3,3), eps_c(3,3);
  // new values of driving forces
  vector pyt(3), pyc(3);
  // actual principal values of damage tensors
  vector pdt(3), pdc(3);
  // actual strains and principal strains
  vector epsn(ncompstr);
  vector peps(3);
  // vector strains without inelastic part
  vector eps(ncompstr);
  // tensor form of the above vector and transformation matrix for principal direction of strains
  matrix epst(3,3), t(3,3);
  // vector and tensorial form of inelastic strains due to damage  (actualized for fatigue only)
  vector epsp(ncompstr);
  matrix epspt(3,3);
  // auxiliary matrices for damage
  matrix omegat(3,3), omegac(3,3);
  // auxiliary matrix
  matrix tmp(3,3);
  // auxiliary transformation matrix for damage
  matrix tomega(3,3);
  // stress tensor and stress vector
  matrix sigt(3,3);
  // auxiliary vector for stress computation
  vector sig(ncompstr);
  // resulting stress vector
  vector sigma(ncompstr);
  // resulting principal stresses
  vector psig(3);
  // Young modulus, Poissons ratio, average size of element
  double e, nu, h;
  // element id for the given ipp
  long eid;

  if ((Mp->matmodel == nonlocal) && (Mm->ip[ipp].hmt & 2))
  // nonlocal model is used in the problem and the given ip point has nonlocal model
  {
    if (Mp->nonlocphase == 1)
    {
      // in case of nonlocal damage in the first phase nothing is need to compute
      // pure strain vector will be averaged
      return;
    }
    else
    {
      //  new total averaged strains
      for (i=0;i<ncompstr;i++)
        epsn[i] = Mm->ip[ipp].nonloc[i];
    }
  }
  else
  // local model is used in this ip
  {
    //  total new strains
    for (i=0;i<ncompstr;i++)
      epsn[i] = Mm->ip[ipp].strain[i];
  }

  // initialize values of material parameter a
  give_actual_param_a(ipp, ido, aat, aac);
  e  = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    epsn[3] = -nu/(1.0-nu)*(epsn[0]+epsn[1]);


  // intialize damage tensors with previous values
  // D_t for tension
  dt[0][0] = Mm->ip[ipp].eqother[ido+2+0];
  dt[1][1] = Mm->ip[ipp].eqother[ido+2+1];
  dt[2][2] = Mm->ip[ipp].eqother[ido+2+2];
  dt[1][2] = dt[2][1] = Mm->ip[ipp].eqother[ido+2+3];
  dt[0][2] = dt[2][0] = Mm->ip[ipp].eqother[ido+2+4];
  dt[0][1] = dt[1][0] = Mm->ip[ipp].eqother[ido+2+5];
  // D_c for compression
  dc[0][0] = Mm->ip[ipp].eqother[ido+2+6+0];
  dc[1][1] = Mm->ip[ipp].eqother[ido+2+6+1];
  dc[2][2] = Mm->ip[ipp].eqother[ido+2+6+2];
  dc[1][2] = dc[2][1] = Mm->ip[ipp].eqother[ido+2+6+3];
  dc[0][2] = dc[2][0] = Mm->ip[ipp].eqother[ido+2+6+4];
  dc[0][1] = dc[1][0] = Mm->ip[ipp].eqother[ido+2+6+5];
    
  // inelastic strains due to previous damage - fatigue model
  if (fat == fatigue_on)
  {
    epspt[0][0] = Mm->ip[ipp].eqother[ido+2+6+6+0];
    epspt[1][1] = Mm->ip[ipp].eqother[ido+2+6+6+1];
    epspt[2][2] = Mm->ip[ipp].eqother[ido+2+6+6+2];
    epspt[1][2] = epspt[2][1] = Mm->ip[ipp].eqother[ido+2+6+6+3];
    epspt[0][2] = epspt[2][0] = Mm->ip[ipp].eqother[ido+2+6+6+4];
    epspt[0][1] = epspt[1][0] = Mm->ip[ipp].eqother[ido+2+6+6+5];
    tensor_vector(epsp, epspt, Mm->ip[ipp].ssst, strain); 
    subv(epsn,epsp,eps);
  }
  else
    copyv(epsn, eps);

  // new driving forces
  vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
  princ_val(epst, peps, t, nijac, limit, Mp->zero, 3, 1);
  // transformation of damage tensors to coordinate system of principal strains
  glmatrixtransf(dt, t);
  glmatrixtransf(dc, t);
  // initializing principal damage with old values
  for (i=0; i<3; i++)
  {
    pdt[i] = dt[i][i];
    pdc[i] = dc[i][i];
  }

  // new principal damage parameters
  princ_dam(ipp, peps, aat, aac, pdt, pdc);
  for (i=0; i<3; i++)
  {
    if (pdt[i] > dt[i][i])
      dt[i][i] = pdt[i];
    if (pdc[i] > dc[i][i])
      dc[i][i] = pdc[i];
  }
  lgmatrixtransf(dt, t);
  lgmatrixtransf(dc, t);
  
  if ((Mp->matmodel == nonlocal) && (Mm->ip[ipp].hmt & 2))
  // nonlocal model is used in the problem and the given ip point has nonlocal model
  {
    //  stresses have to be computed from the nonaveraged total strains
    for (i=0;i<ncompstr;i++)
      epsn[i] = Mm->ip[ipp].strain[i];
  }

  // compute Omega_t=(1-Dt)^(1/2)
  princ_val(dt, pdt, tomega, nijac, limit, Mp->zero, 3, 1);
  fillm(0.0, omegat);
  for(i=0; i<3; i++)
  {
    if ((pdt[i] < -1.0e-8) || (pdt[i] > 1.00000001))
    {
//      print_err("principal value of damage tensor for tension is out of range <0,1>.", __FILE__, __LINE__, __func__);
      if (pdt[i] < 0.0)
        pdt[i] = 0.0;
      else
        pdt[i] = 1.0;
    }
    omegat[i][i] = sqrt(1.0-pdt[i]);
  }
  lgmatrixtransf(omegat, tomega);
  // inelastic strains due to actual damage - fatigue model
  if (fat == fatigue_on)
  {
    fillm(0.0, epspt);
    for (i=0; i<3; i++)
      epspt[i][i]  = 0.5*betat*pdt[i]*pdt[i]/(1.0-pdt[i]);
    lgmatrixtransf(epspt, tomega);
    tensor_vector(epsp, epspt, Mm->ip[ipp].ssst, strain);
    subv(epsn,epsp,eps);
  }
  else
    copyv(epsn, eps);

  // compute Omega_c=(1-Dc)^(1/2)
  princ_val(dc, pdc, tomega, nijac, limit, Mp->zero, 3, 1);
  fillm(0.0, omegac);
  for(i=0; i<3; i++)
  {
    if ((pdc[i] < -1.0e-8) || (pdc[i] > 1.00000001))
    {
//      print_err("principal value of damage tensor for compression is out of range <0,1>.", __FILE__, __LINE__, __func__);
      if (pdc[i] < 0.0)
        pdc[i] = 0.0;
      else
        pdc[i] = 1.0;
    }
    omegac[i][i] = sqrt(1.0-pdc[i]);
  }
  lgmatrixtransf(omegac, tomega);
  // inelastic strains due to actual damage - fatigue model
  if (fat == fatigue_on)
  {
    fillm(0.0, epspt);
    for (i=0; i<3; i++)
      epspt[i][i] -= 0.5*betac*pdc[i]*pdc[i]/(1.0-pdc[i]);
    lgmatrixtransf(epspt, tomega);
    tensor_vector(epsp, epspt, Mm->ip[ipp].ssst, strain);
    subv(epsn,epsp,eps);
  }
  else
    copyv(epsn, eps);

  // decomposition of the strain tensor to the tensile and compressive part
  // eps = eps_t + eps_c
  vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
  princ_val(epst, peps, t, nijac, limit, Mp->zero, 3, 1);

  // compute tensile part of strain tensor and transform it to vector form
  fillm(0.0, epst);
  for(i=0; i<3; i++)
  {
    if (peps[i] > 0.0)
      epst[i][i] = peps[i];   
  }
  lgmatrixtransf(epst, t);
  tensor_vector(eps, epst, Mm->ip[ipp].ssst, strain);
  // eps contains eps_t vector now
  // compute effective tensile stress sig_t = D*eps_t
  Mm->elmatstiff(d, ipp);
  mxv(d, eps, sig);
  vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
  // compute sig_t = Omega_t*sig_t*Omega_t
  mxm(omegat, sigt, tmp);
  mxm(tmp, omegat, sigt);
  tensor_vector(sigma, sigt, Mm->ip[ipp].ssst, stress);

  // compute compressive part of strain tensor and transform it to vector form
  fillm(0.0, epst);
  for(i=0; i<3; i++)
  {
    if (peps[i] < 0.0)
      epst[i][i] = peps[i];
  }
  lgmatrixtransf(epst, t);
  tensor_vector(eps, epst, Mm->ip[ipp].ssst, strain);
  // eps contains eps_c vector now
  // compute effective compressive stress sig_c = D*eps_c
  Mm->elmatstiff(d, ipp);
  mxv(d, eps, sig);
  vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
  // compute sig_c = Omega_c*sig_c*Omega_c
  mxm(omegat, sigt, tmp);
  mxm(tmp, omegat, sigt);
  tensor_vector(sig, sigt, Mm->ip[ipp].ssst, stress);
  
  // compute resulting stress tensor  sigma = sig_t+sig_c
  addv(sigma, sig, sigma);

  //  new data storage
  for (i=0;i<ncompstr;i++)
    Mm->ip[ipp].stress[i]=sigma[i];

  Mm->ip[ipp].other[ido+2+0] = dt[0][0];
  Mm->ip[ipp].other[ido+2+1] = dt[1][1];
  Mm->ip[ipp].other[ido+2+2] = dt[2][2];
  Mm->ip[ipp].other[ido+2+3] = dt[1][2];
  Mm->ip[ipp].other[ido+2+4] = dt[0][2];
  Mm->ip[ipp].other[ido+2+5] = dt[0][1];
  // for compression
  Mm->ip[ipp].other[ido+2+6+0] = dc[0][0];
  Mm->ip[ipp].other[ido+2+6+1] = dc[1][1];
  Mm->ip[ipp].other[ido+2+6+2] = dc[2][2];
  Mm->ip[ipp].other[ido+2+6+3] = dc[1][2];
  Mm->ip[ipp].other[ido+2+6+4] = dc[0][2];
  Mm->ip[ipp].other[ido+2+6+5] = dc[0][1];
    
  // inelastic strains due to previous damage - fatigue model
  if (fat == fatigue_on)
  {
    subv(epsn, eps, epsp);
    Mm->ip[ipp].other[ido+2+6+6+0] = epsp[0];
    Mm->ip[ipp].other[ido+2+6+6+1] = epsp[1];
    Mm->ip[ipp].other[ido+2+6+6+2] = epsp[2];
    Mm->ip[ipp].other[ido+2+6+6+3] = epsp[3];
    Mm->ip[ipp].other[ido+2+6+6+4] = epsp[4];
    Mm->ip[ipp].other[ido+2+6+6+5] = epsp[5];
  }

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
      abort();
  }
  vector_tensor(sigma, sigt, Mm->ip[ipp].ssst, stress); 
  princ_val(sigt, psig, t, nijac, limit, Mp->zero, 3, 1);
  for (i = 0; i < 3; i++)
  {
    if (peps[i] > 0.0)
      Mm->ip[ipp].other[ido+2+6+6+6+i] =  pdt[i]*(peps[i]-psig[i]/e)*h;
    else
      Mm->ip[ipp].other[ido+2+6+6+6+i] = -pdc[i]*(peps[i]-psig[i]/e)*h;
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
void ortodamrot::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}
