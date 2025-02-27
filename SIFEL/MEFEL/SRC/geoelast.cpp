#include "geoelast.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "tensor.h"



/**
  The constructor inializes attributes to zero values.
  
  Created by Tomas Koudelka,
*/
geoelastmat::geoelastmat (void)
{
  keu = 0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka,
*/
geoelastmat::~geoelastmat (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void geoelastmat::read (XFILE *in)
{
  xfscanf (in,"%lf", &keu);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns computed stiffness %matrix in the parameter d.
  
  Created by Tomas Koudelka,
*/
void geoelastmat::matstiff (matrix &d, long ipp,long ido)
{
  Mm->elmatstiff(d, ipp, ido);
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void geoelastmat::nlstresses (long ipp, long im, long ido)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n),sigt(6);
  matrix d(n,n);
  double i1s;
  
  // compute plasticity
  Mm->computenlstresses(ipp, Mm->ip[ipp], im+1, ido+1);
  // extract computed stresses
  Mm->givestress (0, ipp, sig);
  //vector_tensor (sig, sigt, Mm->ip[ipp].ssst, stress);
  give_full_vector(sigt,sig,Mm->ip[ipp].ssst);
  i1s = first_invar (sigt)/3;
  if (i1s <= Mm->ip[ipp].eqother[ido])
  // loading state
  {
    Mm->storestress(0, ipp, sig);
    // update reached pressure
    Mm->ip[ipp].other[ido]=i1s;
  }
  else
  // unloading state
  {
    // assembling unloading matrix
    cmulm(keu, d);
    // give total strains
    Mm->givestrain(0, ipp, eps);
    // elastic strain = total strain - plastic strain
    for (i=0; i<n; i++)
      eps[i] -= Mm->ip[ipp].eqother[ido+1+i];
    mxv(d, eps, sig);
    Mm->storestress(0, ipp, sig);
    // copy reached pressure to current state
    Mm->ip[ipp].other[ido] = Mm->ip[ipp].eqother[ido];
  }
}



/**
  The function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void geoelastmat::updateval (long ipp, long im, long ido)
{
  // update reached pressure
  Mm->ip[ipp].eqother[ido] = Mm->ip[ipp].other[ido];  
  // update values of plasticity model
  Mm->updateipvalmat (ipp,im+1,ido+1);
}
