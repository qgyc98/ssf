#include "shrinkmat.h"
#include "creep.h"
#include "creep_b3.h"
#include "creep_rspec.h"
#include "creepb.h"
#include "creep_dpl.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "galias.h"
#include <math.h>



/**
  This constructor inializes attributes to zero values.

  Created 21. 11. 2013 by JK+TKo
*/
shrinkmat::shrinkmat (void)
{
  tshr = shr_measured;
  thumdef = no;
}



/**
  This destructor is only for the formal purposes.

  Created 21. 11. 2013 by JK+TKo
*/
shrinkmat::~shrinkmat (void)
{

}



/**
   Function reads input material parameters
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::read (XFILE *in)
{
  xfscanf(in, "%m", &tshrlaw_kwdset, &tshr);
  switch (tshr)
  {
    case shr_measured:
      shmeas.read(in);
      break;
    case beta_val:
      beta.read(in);
      break;
    default:
      print_err("unknown shrinkage law is required", __FILE__, __LINE__, __func__);
  }
  xfscanf(in, "%m", &answertype_kwdset, &thumdef);
  switch (thumdef)
  {
    case no:  // Relative humidity is not defined in the matrial. 
              // It is taken from the nonmechq which is calculated externally out of MEFEL
      break;
    case yes: // humidity is defined in the material by time dependent function
      humdef.read(in);
      break;
    default:
      print_err("invalid definition of relative humidity", __FILE__, __LINE__, __func__);
  }  
}



/**
   Function reads input material parameters
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::print (FILE *out)
{
  fprintf(out, "%d\n", tshr);
  switch (tshr)
  {
    case shr_measured:
      shmeas.print(out);
      break;
    case beta_val:
      beta.print(out);
      break;
    default:
      print_err("unknown shrinkage law is required", __FILE__, __LINE__, __func__);
  }
  fprintf(out, "%d\n", thumdef);
  switch (thumdef)
  {
    case no:  // Relative humidity is not defined in the matrial. 
              // It is taken from the nonmechq which is calculated externally out of MEFEL
      break;
    case yes: // humidity is defined in the material by time dependent function
      humdef.print(out);
      break;
    default:
      print_err("invalid definition of relative humidity", __FILE__, __LINE__, __func__);
  }  
}



/**
   This function computes material stiffnes %matrix.
   
   @param d - allocated matrix structure for material stiffness %matrix
   @param ipp - integration point number
   @param im - material index
   @param ido - index of internal variables for given material in the ipp eqother array
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::matstiff (matrix &d,long ipp,long im, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  Mm->matstiff(d, ipp, im+1, ido+1+2*ncomp);
}



/**
   This function computes correct stresses in the integration point and stores
   them into ip stress array.
   
   @param ipp - integration point pointer
   @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
   @param ido - index of internal variables for given material in the ipp other array
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::nlstressesincr (long ipp, long im, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  long ncompo = 1+2*ncomp;
  double hum, hum_prev, deps_sh, dshdh;
  vector sig(ASTCKVEC(ncomp)), eps_sh(ASTCKVEC(ncomp));
  matrix epstens(ASTCKMAT(3,3)), d(ASTCKMAT(ncomp, ncomp));

  // computation of stress increment from other materials
  Mm->computenlstressesincr (ipp, im+1, ido+ncompo);
  
  
  //  computation of stress increment caused by increment of shrinkage strains
  
  // actual volumetric moisture content
  if (thumdef == no)
    hum=Mm->givenonmechq(vol_moist_cont, ipp);
  else
    hum = humdef.getval(Mp->time);
  // volumetric moisture content from previous time step
  hum_prev = Mm->ip[ipp].eqother[ido];
  
  // isotropic increment of shrinkage strain
  // ---------------------------------------
  // The increment of shrinkage strain must be taken with positive sign because forces due to 
  // the corresponding stress incerement will be used on the right hand side of the equilibrium 
  // condition. Shrinkage strains are stored in eqother array and they are subtracted from 
  // total strains in the nlstresses procedure
  dshdh = 0.0;
  switch (tshr)
  {
    case shr_measured:
      dshdh = shmeas.getderiv(hum);
      break;
    case beta_val:
      dshdh = beta.getval(hum);
      break;
    default:
      print_err("unknown shrinkage law is required", __FILE__, __LINE__, __func__);
  }
  deps_sh = dshdh*(hum-hum_prev);
  
  epstens[0][0]=deps_sh;
  epstens[1][1]=deps_sh;
  epstens[2][2]=deps_sh;

  // now, the eps_sh will represent strain increment due to shrinkage  
  tensor_vector (eps_sh,epstens,Mm->ip[ipp].ssst,strain);
  
  // compute stress increment due to increments of shrinkage strains
  Mm->elmatstiff(d, ipp);
  mxv(d, eps_sh, sig);
  // store stress increments
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].eqother[ido+1+ncomp+i] = sig[i];

  //  actualize volumetric moisture content
  Mm->ip[ipp].other[ido] = hum;

  //  actualize total shrinkage strains 
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].other[ido+1+i] = Mm->ip[ipp].eqother[ido+1+i] + eps_sh[i];
}



/**
   This function computes correct stresses in the integration point and stores
   them into ip stress array.
   
   @param ipp - integration point pointer
   @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
   @param ido - index of internal variables for given material in the ipp other array
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::nlstresses (long ipp, long im, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  long ncompo;
  double nu;
  vector epsback(ASTCKVEC(ncomp));
  vector eps_sh(ASTCKVEC(ncomp));

  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }

  //  backup of the total strains
  for (i=0;i<ncomp;i++)
    epsback[i] = Mm->ip[ipp].strain[i];
  
  
  //  actual total shrinkage strains 
  for (i=0;i<ncomp;i++){
    eps_sh[i] = Mm->ip[ipp].other[ido+1+i];
  }
  
  //  initial values of elastic strain = total strain - shrinkage(irreversible) strain
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] -= eps_sh[i];
  
  // number of components in eqother array for shrinkage
  ncompo = 1+2*ncomp;
  
  Mm->computenlstresses(ipp,Mm->ip[ipp], im+1, ido+ncompo);
  
  //  recovery of the total strains and total stress increment
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] = epsback[i];  
}



/**
   This function computes correct stresses in the integration point and stores
   them into ip stress array.
   
   @param ipp - integration point pointer
   @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
   @param ido - index of internal variables for given material in the ipp other array
   
   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::nonloc_nlstresses (long ipp, long im, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  long ncompo;
  double nu;
  vector epsback(ASTCKVEC(ncomp));
  vector eps_sh(ASTCKVEC(ncomp));

  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }

  //  backup of the total strains
  for (i=0;i<ncomp;i++)
    epsback[i] = Mm->ip[ipp].strain[i];
  
  
  //  actual total shrinkage strains 
  for (i=0;i<ncomp;i++){
    eps_sh[i] = Mm->ip[ipp].other[ido+1+i];
  }
  
  //  initial values of elastic strain = total strain - shrinkage(irreversible) strain
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] -= eps_sh[i];
  
  // number of components in eqother array for shrinkage
  ncompo = 1+2*ncomp;
  
  Mm->compnonloc_nlstresses(ipp, im+1, ido+ncompo);
  
  //  recovery of the total strains and total stress increment
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] = epsback[i];  
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array

  Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::updateval (long ipp, long im, long ido)
{
  long i;
  //  other state variables
  long ncomp = Mm->ip[ipp].ncompstr;
  //  actualize volumetric moisture content
  Mm->ip[ipp].eqother[ido] = Mm->ip[ipp].other[ido];

  //  actualize total shrinkage strains 
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].eqother[ido+1+i] = Mm->ip[ipp].other[ido+1+i];

  Mm->updateipvalmat (ipp,im+1,ido+1+2*ncomp);
}



/**
   Function initializes eqother array with initial values. Actual values of 'vol_moist_cont' from 
   array Mm->nonmechq are taken as initial moisture.
   
   @param lcid - load case id
   @param ipp - integration point pointer
   @param im - material index
   @param ido - index of internal variables for given material in the ipp other array
   @param[in] rinit - flag for initialization after restorage from hdbackup

   Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::initvalues (long lcid, long ipp, long im, long ido, bool rinit)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  double hum;

  // actual value of volumetric moisture content = intial humidity
  if (thumdef == no)
    hum=Mm->givenonmechq(vol_moist_cont, ipp);
  else
    hum = humdef.getval(Mp->time);

  Mm->ip[ipp].eqother[ido] = hum;

  Mm->initvalues (lcid,ipp,im+1,ido+1+2*ncomp, rinit);
}



/**
  This function returns the value of tensile strength

  @param lcid - load case id
  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param fi - first index of the required stress increment component
  @param sig - %vector containing stress increment components (output)
  
  Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::givestressincr (long lcid, long ipp, long im, long ido, long fi,vector &sig)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  long i, j;
  
  // restore stress increments for other models
  Mm->givestressincr(lcid, ipp, im+1, ido+1+2*ncomp, fi, sig);

  // add stored stress increments due to shrinkage
  for (i=fi,j=0; j<sig.n; i++,j++)
    sig[j] += Mm->ip[ipp].eqother[ido+1+ncomp+i];  
}



/**
  This function returns the vector of irreversible strains

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param epsirr - %vector containing irreversible strain components (output)
  
  Created 21. 11. 2013 by JK+TKo
*/
void shrinkmat::giveirrstrains(long ipp, long im, long ido, vector &epsirr)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  long i;

  // restore irreversible strains for other models
  Mm->giveirrstrains(ipp, im+1, ido+1+2*ncomp, epsirr);

  // add stored strain due to shrinkage
  for (i=0;i<ncomp;i++)
    epsirr[i] += Mm->ip[ipp].eqother[ido+1+i];  
}



/**
  This function returns the value of tensile strength

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  Created 21. 11. 2013 by JK+TKo
*/
double shrinkmat::give_actual_ft (long ipp, long im, long ido)
{
  double ft;
  long ncomp = Mm->ip[ipp].ncompstr;

  ft = Mm->give_actual_ft(ipp, im+1, ido+1+2*ncomp);

  return ft;
}



/**
  This function returns the value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  Created 21. 11. 2013 by JK+TKo
*/
double shrinkmat::give_actual_fc (long ipp, long im, long ido)
{
  double fc;  
  long ncomp = Mm->ip[ipp].ncompstr;

  // strength is obtained from damage material;
  fc = Mm->give_actual_fc(ipp, im+1, ido+1+2*ncomp);

  return fc;
}



/**
   The funtion marks required non-mechanical quantities in the array anmq.
   
   @param anmq - array with flags for used material types
                 anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                 anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

   @return The function does not return anything, but it may change content of anmq array.
*/
void shrinkmat::give_reqnmq(long *anmq)
{
  if (thumdef == no)
    anmq[vol_moist_cont-1] = 1;
}
