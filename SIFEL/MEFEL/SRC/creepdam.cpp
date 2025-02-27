#include "creepdam.h"
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
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>

/**
  This constructor inializes attributes to zero values.
*/
creepdam::creepdam (void)
{

}

/**
  This destructor is only for the formal purposes.
*/
creepdam::~creepdam (void)
{

}

/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array


*/
void creepdam::matstiff (matrix &d,long ipp,long im, long ido)
{
  long ncompo = Mm->givencompeqother(ipp, im+1); 
  Mm->matstiff(d, ipp, im+2, ido+ncompo);
}


/**
  This function computes correct stresses in the integration point and stores
  them into ip stress array. It was used for Hinkley computations and supposes creep 
  strains to be read from backup file.
  
  @param ipp - integration point pointer
  7.10.2001

void creepdam::nlstresses (long ipp)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  vector epsback(ncomp);


  //  backup of the total strains
  for (i=0;i<ncomp;i++){
    epsback[i] = Mm->ip[ipp].strain[i];
  }

  // compute creep - actually it is supposed the creep strains are read from backup file
  //                 to the eqother array.
  //  initial values of elastic strain for damage = total strain - creep strain
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].strain[i] -= Mm->ip[ipp].eqother[i];
  }
  
  // compute damage, it has index 1, eqother values for damage starts at ncomp index
  Mm->computenlstresses(ipp, 1, 2*ncomp);
  
  //  recovery of the total strains
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].strain[i] = epsback[i];
  }

}
*/




/**
  This function computes stresses increments fro creep model
  in the integration point and stores them into ip stress array.
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp other array

  Created by Tomas Koudelka  7.2008
*/
void creepdam::nlstressesincr(long ipp, long im, long ido)
{
  double nu;

  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }

  // compute creep , it has index1, creep eqother values starts from 0.
  Mm->computenlstressesincr(ipp, im+1, ido);
  
  return;
}



/**
  This function computes correct stresses in the integration point and stores
  them into ip stress array.
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp other array

  7.10.2001
*/
void creepdam::nlstresses (long ipp, long im, long ido)
{
  long i, idem;
  long ncomp = Mm->ip[ipp].ncompstr;
  long ncompo;
  double nu;
  vector epsback(ncomp);
  vector epsirr(ncomp);


  if (Mm->ip[ipp].ssst == planestress)
  {
    idem = Mm->ip[ipp].gemid();
    nu = Mm->eliso[Mm->ip[ipp].idm[idem]].nu;
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }

  //  backup of the total strains
  for (i=0;i<ncomp;i++)
    epsback[i] = Mm->ip[ipp].strain[i];
  //  initial values of elastic strain for damage = total strain - irreversible strain
  Mm->giveirrstrains(ipp, im+1, ido, epsirr);
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] -= epsirr[i];

  // number of components in eqother array for creep 
  ncompo = Mm->givencompeqother(ipp, im+1); 
  // compute damage, it has index 2, eqother values for damage starts at ncompo index
  if ((Mm->ip[ipp].hmt & 2) && (Mp->nonlocphase == 2))// nonlocal damage
    Mm->compnonloc_nlstresses(ipp, im+2, ido+ncompo);
  else
    Mm->computenlstresses(ipp, Mm->ip[ipp], im+2, ido+ncompo);

  //  recovery of the total strains and total stress increment
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] = epsback[i];
}



/**
  This function computes correct stresses for the nonlocal damage models 
  in the integration point and stores them into ip stress array.
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp other array

  7.10.2001
*/
void creepdam::nonloc_nlstresses (long ipp, long im, long ido)
{
  long i, idem;
  long ncomp = Mm->ip[ipp].ncompstr;
  long ncompo;
  double nu;
  vector epsback(ncomp);
  vector epsirr(ncomp);


  if (Mm->ip[ipp].ssst == planestress)
  {
    idem = Mm->ip[ipp].gemid();
    nu = Mm->eliso[Mm->ip[ipp].idm[idem]].nu;
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }

  //  backup of the total strains
  for (i=0;i<ncomp;i++)
    epsback[i] = Mm->ip[ipp].strain[i];
  //  initial values of elastic strain for damage = total strain - irreversible strain
  Mm->giveirrstrains(ipp, im+1, ido, epsirr);
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].strain[i] -= epsirr[i];

  // number of components in eqother array for creep 
  ncompo = Mm->givencompeqother(ipp, im+1); 
  // compute damage, it has index 2, eqother values for damage starts at ncompo index
  if ((Mm->ip[ipp].hmt & 2) && (Mp->nonlocphase == 2))// nonlocal damage
    Mm->compnonloc_nlstresses(ipp, im+2, ido+ncompo);
  else
    Mm->compnonloc_nlstresses(ipp, im+2, ido+ncompo);

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

*/
void creepdam::updateval (long ipp, long im, long ido)
{
  //  update values of plasticity model
  //  long ncompo = 2*Mm->ip[ipp].ncompstr; // for Hinkley
  //  Mm->updateipvalmat (ipp,1,ncompo); // for Hinkley
  long ncompo = Mm->givencompeqother(ipp, im+1); 
  Mm->updateipvalmat (ipp,im+1,ido); 
  Mm->updateipvalmat (ipp,im+2,ido+ncompo);  
}



/**
   Function initializes eqother array with initial values.
   Actual values of quantities 'temperature' and 'rel_hum' from array 
   Mm->nonmechq are taken as initial temperature and initial moisture.
   
   @param lcid - load case id
   @param ipp - integration point pointer
   @param im - material index
   @param ido - index of internal variables for given material in the ipp other array
  @param[in] rinit - flag for initialization after restorage from hdbackup

   7.6.2005, TKo
*/
void creepdam::initvalues (long lcid, long ipp, long im, long ido, bool rinit)
{
  long ncompo = Mm->givencompeqother(ipp, im+1); 
  Mm->initvalues(lcid,ipp,im+1,ido, rinit); 
  Mm->initvalues(lcid,ipp,im+2,ido+ncompo, rinit); 
}



/**
  This function returns the value of tensile strength

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  TKo
*/
double creepdam::give_actual_ft (long ipp, long im, long ido)
{
  double ft;
  ft = Mm->give_actual_ft(ipp, im+1, ido);
  return ft;
}



/**
  This function returns the value of compressive strength

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  10.2008 TKo
*/
double creepdam::give_actual_fc (long ipp, long im, long ido)
{
  double fc;  
  long ncompo_cr = Mm->givencompeqother(ipp, im+1);
  // strength is obtained from damage material;
  fc = Mm->give_actual_fc(ipp, im+2, ido+ncompo_cr);
  return fc;
}
