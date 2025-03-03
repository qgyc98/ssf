#include "nonlocmicroM4.h"
#include "global.h"
#include "mechmat.h"
#include "mechtop.h"
#include "intpoints.h"
#include "microM4.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


nonlocmicroM4::nonlocmicroM4 (void)
{
}
nonlocmicroM4::~nonlocmicroM4 (void)
{
}

void nonlocmicroM4::read(XFILE *in)
{
  xfscanf(in,"%lf",&r);
}

/**
  This function returns number of averaged quantities.

  @param ipp - integration point number

  @return The function returns number of averaged quantities for microplane M4 model.
*/
long nonlocmicroM4::give_num_averq (long ipp)
{
  return Mm->ip[ipp].ncompstr;
}


/**
  This function averages values of the gamma or plastic strains in the given integration point.
  Which values are averaged is defined by the data member waf.

  @param ipp - integration point number

*/
void nonlocmicroM4::average (long ipp)
{
  matrix e(6,6);//elastic stiffness matrix
  vector S1(6);
  double sigma_elast,sigma_loc,sigma_inelast;
  double S2=0;
  double rr,alfa;
  long i,j,k,aip,nad;

  //long numberOfMicroplanes = Mm->mpM4[0].numberOfMicroplanes;

  //zero matrix and vector
  fillv(0.0,S1);
  //nad .... number of adjacent points
  nad=Mt->nadjip[ipp];

  //loop over all adjacent points
  for (k=0;k<nad;k++){


    //number of adj. int. point
    aip=Mt->adjip[ipp][k];
    //compute distance of int. points ipp<->aip
    rr=Mt->dist[ipp][k];
    //compute nonlocal weight
    alfa=(rr>=r)?0:(1-rr*rr/(r*r));
    // multiply by the ip volume
    alfa *= Mm->ipv[aip];
    //compute elas. stiff. matrix for aip
    Mm->elmatstiff(e,aip);

    //loop over stress components
    for (i=0;i<6;i++){
      //compute elastic part of stress in aip
      sigma_elast=0;
      for (j=0;j<6;j++){
	sigma_elast += e[i][j] * (Mm->ip[aip].strain[j]);
      }
      //take local stress component stored in other[..x..y...]
      //v seznamu je take ipp
      sigma_loc=Mm->ip[aip].stress[i];
      //inelastic part of stress
      sigma_inelast=sigma_loc-sigma_elast;
      //nonlocal stress contribution
      S1[i]+=alfa*sigma_inelast;
    }    //loop i

    S2+=alfa;
  }//loop k


  for(i=0;i<6;i++)
    Mm->ip[ipp].nonloc[i] = S1[i]/S2;


}
