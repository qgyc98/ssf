#include "nonlocdamg.h"
#include "global.h"
#include "mechtop.h"
#include "mechmat.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "vecttens.h"
#include "mathem.h"
#include <string.h>



/**
  The constructor inializes attributes to zero values.

  Created by Tomas Koudelka
*/
nonlocdamg::nonlocdamg(void)
{
  r = 0;
  af = avrgf(0);
}



/**
  The destructor is defined only for the formal purposes.

  Created by Tomas Koudelka
*/
nonlocdamg::~nonlocdamg(void)
{
  r = 0;
  af = avrgf(0);
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  Parameters:
  @param in - pointer to the opened text file

  Returns:
  @retval 0 - on success
  @retval 1 - on reading error

  Created by Tomas Koudelka
*/
long nonlocdamg::read(XFILE *in)
{
  if (xfscanf(in, "%k%le%k%m", "radius", &r, "aver_func", &avrgf_kwdset, &af) == 4)
    return 0;

  print_err("cannot read parameters of nonlocal damage",__FILE__,__LINE__,__func__);
  return 1;
}



/**
  The function prints material parameters from the opened text file given
  by the parameter out.

  Parameters:
  @param out - pointer to the opened text file

  Returns:
  @return - The function does not return anything.

  Created by Tomas Koudelka
*/
void nonlocdamg::print(FILE *out)
{
  fprintf(out, "%le %d\n", r, int(af));
}



/**
  The function returns number of averaged quantities.

  Parameters:
  @param ipp - integration point number
  @param im - material index

  Returns: 
  @retval The function returns number of averaged quantities

  Created by Tomas Koudelka
*/
long nonlocdamg::give_num_averq (long ipp, long im)
{
  switch(Mm->ip[ipp].tm[im+1])
  {
    case scaldamage:{
      return(1); // kappa
    }
    case scaldamagecc:{
      return(6); // tensor eps
    }      
    case ortodamage:{
      return(6); // tensor eps
    }
    case ortodamagerot:{
      return(6); // tensor eps
    }
    case anisodamage:{
      return(6); // tensor eps
    }
    case anisodamagerot:{
      return(6); // tensor eps
    }
    default:
      print_err("unknown damage material type is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  The function returns vector of averaged quantities.

  Parameters:
  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param qv - %vector of averaged quantities (output)

  Returns : 
  @retval The function returns %vector of averaged quantities in the parameter qv

  Created by Tomas Koudelka
*/
void nonlocdamg::give_aver_quantv(long ipp, long im, long ido, vector &qv)
{
  switch(Mm->ip[ipp].tm[im+1])
  {
    case scaldamage:
      Mm->give_kappa(ipp, im+1, ido, qv); // kappa
      return;
    case scaldamagecc:
    case ortodamage:
    case ortodamagerot:
    case anisodamage:
    case anisodamagerot:  // tensor eps
    {
      long nc = Mm->ip[ipp].ncompstr;
      long i;
      for (i=0; i<nc; i++)
        qv[i] = Mm->ip[ipp].strain[i];
      return;
    }
    default:
      print_err("unknown damage material type is required", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function averages values of the damage function parameter in the given integration point.
  Type of average function is given by the data member af.

  Parameters:
  @param ipp - integration point number
  @param im  - material index
  @param ido - index of internal variables for given material in the ipp other array

  Returns:
  @retval The function does not return anything.

  Created by Tomas Koudelka
*/
void nonlocdamg::average (long ipp, long im, long ido)
{
  double rr, alpha, sum = 0.0;
  double weight;
  long i, j;
  // number of adjacent integration points
  long nad;
  // id of an adjacent integration point
  long aip;
  // number of averaged quantities
  long naq = Mm->give_num_averq(ipp, im);
  // vector of averaged quantities at the given integration point ipp
  vector aq(naq);
  // auxiliary vector of averaged quantities for particular adjacent integration points
  vector aux(naq);

  //nad .... number of adjacent points
  nad=Mt->nadjip[ipp];

  //loop over all adjacent points
  for (i = 0; i < nad; i++)
  {
    //number of adj. int. point
    aip=Mt->adjip[ipp][i];
    //compute distance of int. points ipp<->aip
    rr=Mt->dist[ipp][i];
    //compute nonlocal weight
    switch (af)
    {
      case parab:
        alpha = (rr >= r) ? 0 : (1-rr*rr/(r*r));
        break;
      case cubic:
        alpha = (rr >= r) ? 0 : (1-rr*rr/(r*r))*(1-rr*rr/(r*r));
        break;
      case exponential:
        alpha = (rr >= r) ? 0 : exp(-sqr(2.0*rr/r));
        break;
      default:
        print_err("unknown type of the averaging function is required",__FILE__,__LINE__,__func__);
        break;
    }
    if (alpha < 1.0e-3)
      continue;
//    sum += alpha;
    weight = alpha * Mm->ipv[aip];
    sum += weight;
    Mm->give_aver_quantv(aip, im, ido, aux);
    for (j=0; j<naq; j++)
      aq[j] += weight*aux[j];
//    aeps += weight*Mm->ip[aip].other[ido+0];
  }//loop i
  for(i=0; i<naq; i++)
    Mm->ip[ipp].nonloc[i] = aq[i] / sum;
}



/**
  The function returns the length of fracture process zone i.e. the diameter of 
  the averaged neighbourhood.

  Parameters:
  @param ipp - integration point number
  @param im  - material index
  @param ido - index of internal variables for given material in the ipp other array

  Returns:
  @return The function returns the length of process zone.

  Created by Tomas Koudelka, 26.8.2014
*/
double nonlocdamg::give_proczonelength(long ipp, long im, long ido)
{
  return 2.0*r;
}
