#include "nonlocplast.h"
#include "global.h"
#include "mechmat.h"
#include "mechtop.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "vecttens.h"


/**
  The constructor inializes attributes to zero values.

  Created by Tomas Koudelka,
*/
nonlocplast::nonlocplast(void)
{
  r = 0;
  waf = wavrg(0);
}



/**
  The destructor is defined only for the formal purposes.

  Created by Tomas Koudelka,
*/
nonlocplast::~nonlocplast(void)
{
  r = 0;
  waf = wavrg(0);
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
long nonlocplast::read(XFILE *in)
{
  if (xfscanf(in, "%le %m", &r, &wavrg_kwdset, (int*)&waf) == 2)
    return 0;

  print_err("cannot read parameters of nonlocal plasticity", __FILE__,__LINE__,__func__);
  return 1;

}



/**
  The function returns number of averaged quantities.

  @param ipp - integration point number

  @return The function returns number of averaged quantities.

  Created by Tomas Koudelka,
*/
long nonlocplast::give_num_averq (long ipp)
{
  switch (waf)
  {
    case avggamma:
      return Mm->ip[ipp].ncompeqother;
    case avgepsp:
      return Mm->ip[ipp].ncompstr;
    default:
      print_err("unknown plasticity material type is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  The function returns %vector of averaged quantities.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param qv - %vector of averaged quantities (output)

  @return The function returns number of averaged quantities in the parameter qv.

  Created by Tomas Koudelka,
*/
void nonlocplast::give_aver_quantv(long ipp, long im, long ido, vector &qv)
{
  switch (waf)
  {
    case avggamma:
      qv[0] = Mm->give_consparam(ipp, im+1, ido); // gamma
      return;
    case avgepsp:
      Mm->giveirrstrains(ipp, im+1, ido, qv); // eps_p
      return;
    default:
      print_err("unknown damage material type is required", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function averages values of the gamma or plastic strains in the given integration point.
  Data member waf defines which values are averaged.

  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void nonlocplast::average (long ipp, long ido)
{
  double rr, alpha, sum = 0.0;
  long ncompo, ncomps = Mm->ip[ipp].ncompstr;
  matrix d(ncomps,ncomps);
  vector dgdst(6),sigt(6);
  double gamma, agamma = 0.0;
  vector dgds(ncomps),aepsp(ncomps), epsa(ncomps), epsn(ncomps), epsp(ncomps), sig(ncomps);
  vector q(1);
  long i, j, nad, aip;

  for (i=0;i<ncomps;i++){
    Mm->ip[ipp].nonloc[i]=0.0;
  }

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
    alpha = (rr >= r) ? 0 : (1-rr*rr/(r*r));
    // multiply by the ip volume
    alpha *= Mm->ipv[aip];
    sum += alpha;

    switch (waf)
    {
      case avggamma:
        ncomps = Mm->ip[aip].ncompstr;
        gamma = Mm->ip[aip].other[ido+ncomps];
        agamma += alpha * gamma;
        break;
      case avgepsp:
        for (j = 0; j < ncomps; j++)
          Mm->ip[ipp].nonloc[j] += alpha*Mm->ip[aip].other[ido+j];
        break;
      default:
        print_err("unknown type of the averageing variable is required.",__FILE__,__LINE__,__func__);
        break;
    }
  }//loop i
  ncomps = Mm->ip[ipp].ncompstr;
  ncompo = Mm->ip[ipp].ncompother;
  switch (waf)
  {
    case avggamma:
    {
      //  initial values
      for (i = 0; i < ncomps; i++)
      {
        epsn[i]=Mm->ip[ipp].strain[i];
        epsp[i]=Mm->ip[ipp].eqother[ido+i];
      }
      //  elastic strain
      subv (epsn,epsp,epsa);
      //  elastic stiffness matrix
      Mm->elmatstiff (d,ipp);
      //  stress computation
      mxv (d,epsa,sig);
      // values of averaged plastic strains are computed
      //vector_tensor (sig, sigt, Mm->ip[ipp].ssst, stress);
      give_full_vector(sigt,sig,Mm->ip[ipp].ssst);
      // obtaining hardening parameters
      q[0] = 0.0;
      if (ncompo-ncomps-1 > 0)
        q[0] = Mm->ip[ipp].nonloc[ncomps+1];
      Mm->dgdsigma (ipp, 0, sigt, q, dgds);
      //vector_tensor (dgds,dgdst,Mm->ip[ipp].ssst, strain);
      give_full_vector(dgdst,dgds,Mm->ip[ipp].ssst);
      cmulv((agamma-Mm->ip[ipp].eqother[ido+ncomps])/sum, dgdst);
      //tensor_vector(aepsp, dgdst, Mm->ip[ipp].ssst, strain);
      give_red_vector(dgdst,aepsp,Mm->ip[ipp].ssst);
      // storing averaged values
      for (i = 0; i < ncomps; i++)
        Mm->ip[ipp].nonloc[i] = aepsp[i];
      Mm->ip[ipp].nonloc[ncomps] = agamma/sum;
      break;
    }
    case avgepsp:
      // storing averaged values
      for (i = 0; i < ncomps; i++)
        Mm->ip[ipp].nonloc[i] /= sum;
      break;
    default:
      print_err("unknown type of the averageing variable is required.",__FILE__,__LINE__,__func__);
      break;
  }
}
