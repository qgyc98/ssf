#include "ortodam2.h"
#include "bmatrix.h"
#include "matrix.h"
#include "vector.h"
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
  
  Created by Tomas Koudelka, 05.2019
*/
ortodam2::ortodam2 (void)
{
  cde = corr_off;
  fat = fatigue_off;
  ft = uft = ult = 0.0;
  fc = ufc = ulc = 0.0;
  pqt = pqc = 0;
  betac = betat = 0.0;
}



/**
  The destructor is only for the formal purposes.
  
  Created by Tomas Koudelka, 05.2019
*/
ortodam2::~ortodam2 (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 05.2019
*/
void ortodam2::read (XFILE *in)
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

  @param i[in]    - principal direction id
  @param peps[in] - %vector of principal strains
  @param psign[in] - nominal principal stresses (De*peps)
  @param de[in]   - elastic stiffness %matrix in principal directions

  @return The function returns calculated equivalent strain norm.

  Created by Tomas Koudelka 05.2019
*/
double ortodam2::dam_eq_strain(long i, vector &peps, vector &psig, matrix &de)
{
  double kappa = 0.0;
  switch (dameqstr)
  {
    case norstrain:
      kappa = fabs(peps(i));
      break;
    case norenergy:
      kappa = psig(i)*peps(i)/de(i,i);
      kappa = sqrt(kappa);
      break;
    case norrankine:
      kappa = psig(i)/de(i,i);
      kappa = fabs(kappa);
      break;
    default:
      print_err(" unknown type of equivalent strain norm is required", __FILE__, __LINE__,__func__);
  }
  return kappa;
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
  
  Created by Tomas Koudelka, 05.2019
*/
double ortodam2::brittle_damage(long ipp, double y, double e, double f, double uf, double omegao)
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
  
  Created by Tomas Koudelka, 05.2019
*/
double ortodam2::qbezier_damage(long ipp, double y, double e, double f, double uf, double ul,double omegao,long pq)
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
    print_err("previous iteration of omega returned value < 0", __FILE__, __LINE__, __func__);

  if (kappa*omegao*h >= ul) // limit crack opening was reached => full damage, omega = 1.0 
    return 1.0;
  tmp = 0.0;
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

  @param ipp[in]    - integration point number
  @param peps[in]   - %vector of principal strains
  @param pde[in]    - elastic stiffness matrix in principal directions 
  @param pdamc[out] - %vector of principal damage parameters for tension (output)
  @param pdamt[out] - %vector of principal damage parameters for compression (output)
  
  @return The function returns principal damage parameters in the parameters pdamt and pdamc.

  Created by Tomas Koudelka, 05.2019
*/
void ortodam2::princ_dam(long ipp, vector &peps, matrix &pde, vector &pdamt, vector &pdamc)
{
  long i;
  double af, auf, aul;
  double kappa;
  vector psig(ASTCKVEC(3));

  // compute principal nominal stresses
  mxv(pde, peps, psig);
  switch(damevf)
  {
    case brittle:
      for(i=0; i<3; i++)
      {           
        kappa = dam_eq_strain(i, peps, psig, pde);
        // tension 
        if (peps(i) >= 0.0)
        {
          // actual tensile strength
          af = Mm->give_actual_ft(ipp);
          // actual value of initial crack opening
          auf = uft*(af/ft);
          pdamt(i) = brittle_damage(ipp, kappa, pde(i,i), af, auf, pdamt(i));
        }
        // compression 
        if (peps(i) < 0.0)
        {
          // actual compressive strength
          af = Mm->give_actual_fc(ipp);
          // actual value of initial crack opening
          auf = ufc*(af/fc);
          pdamc(i) = brittle_damage(ipp, kappa, pde(i,i), af, auf, pdamc(i));
        }
      }
      break;
    case quadbezier:
      for(i=0; i<3; i++)
      {           
        kappa = dam_eq_strain(i, peps, psig, pde);
        // tension 
        if (peps(i) >= 0.0)
          {
            // actual tensile strength
            af = Mm->give_actual_ft(ipp);
            // actual value of initial crack opening
            auf = uft*(af/ft);
            // actual value of limit crack opening
            aul = ult*(af/ft);
            pdamt(i) = qbezier_damage(ipp, kappa, pde(i,i), af, auf, aul, pdamt(i), pqt);
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
          pdamc(i) = qbezier_damage(ipp, kappa, pde(i,i), af, auf, aul, pdamc(i),pqc);
        }
      }
      break;
    default:
      print_err(" unknown type of damage evolution function is required", __FILE__, __LINE__,__func__);
  }
}



/**
  The function returns the value of tensile strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns tensile strength. 

  Created by Tomas Koudelka, 05.2019
*/
double ortodam2::give_actual_ft(long /*ipp*/)
{
  return ft;
}



/**
  The function returns the value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns compressive strength. 

  10.2008 Created by Tomas Koudelka, 05.2019
*/
double ortodam2::give_actual_fc(long /*ipp*/)
{
  return fc;
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 05.2019
*/
void ortodam2::matstiff (matrix &d,long ipp,long /*ido*/)
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

  Created by Tomas Koudelka, 05.2019
*/
void ortodam2::elmatstiff (matrix &d,long ipp)
{
  Mm->elmatstiff (d,ipp);
  //  long idem = Mm->ip[ipp].gemid();         // index of elastic material
  //  long tmp=Mm->givencompeqother(ipp,idem); // number of internal variables for elastic material and thermal material 
  //  long idoem=Mm->givencompeqother(ipp,0) - tmp;  // total number of internal variables
  //  double e = Mm->give_actual_ym(ipp);
  //  double e0 = Mm->give_actual_ym(ipp,idem,idoem); // Young modulus from elastic material
  //  cmulm(e/e0, d);
}




/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 05.2019
*/
void ortodam2::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr;

  // old values of damage parameters
  vector pdto(ASTCKVEC(3)), pdco(ASTCKVEC(3));
  // damage parameters
  vector pdt(ASTCKVEC(3)), pdc(ASTCKVEC(3));
  // strains and principal strains
  vector epsn(ASTCKVEC(ncompstr));
  // tensor form of strain and transformation matrix for principal direction of strains
  matrix epst(ASTCKMAT(3,3)), tmat(ASTCKMAT(3,3));
  // stress tensor and stress vector
  matrix sigt(ASTCKMAT(3,3));
  vector sigma(ASTCKVEC(ncompstr));
  // elastic stiffness matrix - reduced form
  matrix de(ASTCKMAT(ncompstr,ncompstr));
  // elastic stiffness matrix - full form 6x6
  matrix deg(ASTCKMAT(6,6)), pdeg(ASTCKMAT(6,6));
  // elastic stiffness matrix in principal directions
  matrix pde(ASTCKMAT(3,3));
  // principal strains, principal stresses and principal irreversible strains
  vector peps(ASTCKVEC(3)), psig(ASTCKVEC(3)), pepsp(ASTCKVEC(3));
  double h;
  long eid;
  
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
  princ_val (epst, peps, tmat, nijac, limit, Mp->zero, 3, 1);
  
  // inelastic strains due to previous damage - fatigue model
  if (fat == fatigue_on)
  {
    for (i=0; i<3; i++)
    {
      pepsp[i]  = 0.5*betat*pdto[i]*pdto[i]/(1.0-pdto[i]);
      pepsp[i] -= 0.5*betac*pdco[i]*pdco[i]/(1.0-pdco[i]);
    }
    subv(peps, pepsp, peps);
  }

  // give full elastic stiffness matrix
  Mm->elmatstiff(deg, ipp, spacestress);
  // transformation to the principal coordinate system 
  gl_tens4transf(deg, pdeg, tmat, Mm->ip[ipp].ssst);
  // extract normal component block from principal stiffness matrix(6,6)
  extractm(pde, pdeg, 0, 3);

  // compute principal damage parameters
  princ_dam(ipp, peps, pde, pdt, pdc);
  for (i=0; i<3; i++)
  {
    if (pdt[i] < pdto[i])
      pdt[i] = pdto[i];
    if (pdc[i] < pdco[i])
      pdc[i] = pdco[i];
  }

  // new principal stresses
  mxv(pde, peps, psig);
  for (i=0; i<3; i++)
  {
    if (psig(i) > 0.0)
      psig(i) *= 1.0-pdt(i);
    if (psig(i) < 0.0)
      psig(i) *= 1.0-pdc(i);
  }
   
  // transformation of principal stresses to the global coordinate system
  nullm(sigt);
  for (i = 0; i < 3; i++)
    sigt(i,i) = psig(i);
  lgmatrixtransf(sigt, tmat);
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
  nullm(sigt);
  for (i = 0; i < 3; i++)
    sigt[i][i] = pdt[i];
  lgmatrixtransf(sigt, tmat);
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
      abort();
  }
  // storage of cracks openings computed from tensile damage tensor
  for (i = 0; i < 3; i++)
  {
    if (peps(i) > 0.0)
      Mm->ip[ipp].other[ido+14+i] = pdt[i]*(peps[i]-psig[i]/pde(i,i))*h;
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
void ortodam2::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}
