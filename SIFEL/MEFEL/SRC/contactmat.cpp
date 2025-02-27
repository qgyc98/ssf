#include "contactmat.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"



contactmat::contactmat (void)
{
  ks = 0.0;  kn = 0.0;
}



contactmat::~contactmat (void)
{

}



void contactmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf",&ks, &kn);
}



/**
  The function returns stiffness %matrix of material
   
  @param[out] d - stiffness %matrix of the material
  @param[in] ipp - integration point pointer

  JK, 11.6.2006
*/
void contactmat::matstiff(matrix &d, long ipp)
{
  long n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n);

  switch(n){
    case 2: // 2D interface element
      d[0][0] = ks;
      d[1][1] = kn;
      break;
    case 3: // 3D interface element
      d[0][0] = ks;
      d[1][1] = ks;
      d[2][2] = kn;
      break;
    default:
      print_err("unknown number of strain array components (%ld), only values in range [2;3] are allowed.",
                __FILE__, __LINE__,__func__, n);
      abort();
  }
}



/**
  The function computes nodal forces in contact problems.
  Nodal forces are computed from difference of nodal displacements
  at appropriate nodes.
   
  @param[in] ipp - integration point id
  @param[in] im  - index of material type
  @param[in] ido - index in array other
   
  JF, 04/2013
  TKo, 10/2023
  TKo, 04/2024
*/
void contactmat::nlstresses (long ipp, long /*im*/, long /*ido*/)
{
  long n = Mm->ip[ipp].ncompstr;
  vector eps, sig;
  matrix d(ASTCKMAT(n,n));
  
  // make references to the actual values of strains and stresses
  makerefv(eps, Mm->ip[ipp].strain, n);
  makerefv(sig, Mm->ip[ipp].stress, n);

  // material elastic stiffness matrix
  matstiff(d, ipp);

  // compute stresses according elastic law
  mxv(d, eps, sig);
}



/**
  The function modifies model parametres according to stochastic calculation setup.

  @param[in] atm - attribute selection
  @param[in] val - %vector of new values

  Cretaed by Tomas Koudelka, 04.2024
*/
void contactmat::changeparam (atsel &atm, vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
      case 0:{
        ks=val[i];
        break;
      }
      case 1:{
        kn=val[i];
        break;
      }
      default:{
        fprintf (stderr,"\n\n wrong number of atribute in function changeparam (%s, line %d).\n",__FILE__,__LINE__);
      }
    }
  }
}



/**
  The function copies array other to the array eqother.
   
  @param[in] ipp - number of integration point
  @param[in] im  - index of material type
  @param[in] ido - index in array other
     
*/
void contactmat::updateval (long /*ipp*/, long /*im*/, long /*ido*/)
{
}
