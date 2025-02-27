#include "damplast.h"
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
  The constructor inializes attributes to zero values.
*/
damplast::damplast (void)
{

}

/**
  The destructor is only for the formal purposes.
*/
damplast::~damplast (void)
{

}

/**
  The function computes material stiffnes %matrix.

  Parameters:
  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  Returns:
  @retval The function does not return anything.
*/
void damplast::matstiff (matrix &d, long ipp, long ido)
{
  double dp;
  long ncompo = Mm->givencompother(ipp, 1);

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      Mm->elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff:
      Mm->elmatstiff(d, ipp, ido);
      dp=Mm->ip[ipp].eqother[ncompo+1];
      if (dp > 0.999999)
        dp = 0.999999;
      cmulm (1.0-dp,d);
      break;
    default:
      print_err("unknown type of stifness matrix is required", __FILE__, __LINE__, __func__);
      abort();
  }
}


/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  Parameters:
  @param ipp - integration point pointer
  
  Returns:
  @retval The function does not return anything.
  7.10.2001
*/
void damplast::nlstresses (long ipp)
{
  long i;
  long ncompo = Mm->givencompother(ipp, 1);
  long ncomp = Mm->ip[ipp].ncompstr;
  vector epsback(ncomp);
  vector kappa(1);
  vector sigma(ncomp);

  // compute plasticity, it has index 1
  Mm->computenlstresses(ipp, Mm->ip[ipp], 1, 0);  
  
  //  backup of the total strains
  for (i=0;i<ncomp;i++){
    epsback[i] = Mm->ip[ipp].strain[i];
  }
  //  initial values of elastic strain for damage = total strain - plastic strain
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].strain[i] -= Mm->ip[ipp].other[i];
  }
  
  kappa[0] = Mm->ip[ipp].eqother[ncompo+0];

  // damage stress solver, damage material has index 2
//  dp = Mm->scal_dam_sol (ipp, 2, ncompo, epsback, kappa, sigma);
  Mm->computenlstresses(ipp, Mm->ip[ipp], 2, ncompo);  
  
  //  recovery of the total strains
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].strain[i] = epsback[i];
  }
  
/*  
  //  new data storage
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].stress[i]=sigma[i];
  
  Mm->ip[ipp].other[ncompo+0]=kappa[0];
  Mm->ip[ipp].other[ncompo+1]=dp;
*/  
}

/**
  The function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  Parameters:
  @param ipp - integration point number in the mechmat ip array.

  Returns:
  @retval The function does not return anything.
*/
void damplast::updateval (long ipp)
{
  //  update values of plasticity model
  Mm->updateipvalmat (ipp,1,0);
  //  update values of damage model
  Mm->updateipvalmat (ipp,2,Mm->givencompother(ipp,1));
}



/**
  The function returns the actual value of tensile strength 
  obtained from the damage material parameters.
  
  @param ipp - number of integration point
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.
   
  @return The function returns value of tensile strength.

  Created by Tomas Koudelka, 12.2011
*/
double damplast::give_actual_ft(long ipp, long im, long ido)
{
  double ft;
  long ncompo;
  ncompo = Mm->givencompeqother(ipp, im+1);
  ft = Mm->give_actual_ft(ipp, im+2, ido+ncompo);
  return ft;
}
