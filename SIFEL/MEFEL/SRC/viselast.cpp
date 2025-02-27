#include "viselast.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"

viselast::viselast (void)
{

}


viselast::~viselast (void)
{

}

/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

*/
void viselast::matstiff (matrix &d,long ipp,long /*im*/,long /*ido*/)
{
  //switch (Mp->nlman->stmat){
  //case initial_stiff:{
  Mm->elmatstiff (d,ipp);
    //break;
    //}
    //case tangent_stiff:{
    //Mm->elmatstiff (d,ipp);
    //break;
  //}
  //default:{
  //print_err("unknown type of stifness matrix is required",__FILE__,__LINE__,__func__);
  //}}
  
}

/**
  This function computes material damping %matrix.

  @param d - allocated matrix structure for material damping %matrix
  @param ipp - integration point number

*/
void viselast::matdamp (matrix &d,long ipp,long im,long ido)
{
  Mm->matdamp (d,ipp,im+1,ido);
}


/**
  Function computes true stresses.
   
  @param ipp - number of integration point
   
  @return The function does not return anything.

  Created by JK,
*/
void viselast::nlstresses (long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n);
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  Mm->elmatstiff (d,ipp);
  mxv (d,eps,sig);
  /*
    if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (eps[0]+eps[1]);
  */
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
  {
    for (i=0;i<n;i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }
  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
}
