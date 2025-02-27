#include "elastortomat.h"
#include "global.h"
#include "mechmat.h"
#include "intpoints.h"
#include "elemswitch.h"
#include "vecttens.h"
#include "matrix.h"



elastortomat::elastortomat (void)
{
  allocm(6, 6, de);
  allocm(3, 3, t);
}



elastortomat::~elastortomat (void)
{
}



/**
  Function reads material parameters from the opened text file.
   
  @param[in] in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::read (XFILE *in)
{
  xfscanf(in, "%k", "ortomat");
  readm(in, de);
  xfscanf(in, "%k", "transfmat");
  readm(in, t);
}



/**
  Function prints material parameters into the opened text file.
   
  @param[in] out - pointer to the opened FILE

  @return The function does not return anything.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::print (FILE *out)
{
  printm(de, out);
  printm(out, t);
}



/**
  The function assembles the stiffness %matrix in the global
  coordinate system.

  @param[out] d   - resulting stiffness %matrix
  @param[in]  ipp - integration point id

  @return The funtion returns resulting %matrix in the argument d.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::matstiff (matrix &d, long ipp)
{
  strastrestate ssst = Mm->ip[ipp].ssst;

  if (ssst == spacestress)
    matstiff (d, ipp, ssst);
  else{
    matrix dg(ASTCKMAT(6,6));
    matstiff (dg, ipp, ssst);
    // convert stiffness matrix(6,6) to reduced form depending on ssst
    tensor4_ematrix(d, dg, ssst);
  }
}



/**
  The function assembles the stiffness %matrix in the global
  coordinate system.

  @param[out] d    - resulting stiffness %matrix
  @param[in]  ipp  - integration point id
  @param[in]  ssst - stress/strain state indicator

  @return The funtion returns resulting %matrix in the argument d.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::matstiff (matrix &d, long ipp, strastrestate ssst)
{
  matrix tmat(ASTCKMAT(3,3));
  vector p(ASTCKVEC(3));
  const char *namevar[] = {"x", "y", "z"};

  // x,y,z coordinates of the integration point will be stored in the vector p
  ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
  // evaluate transformation matrix
  t.evaluate(p, namevar, tmat);
  // transformation to the global coordinate system
  lg_tens4transf(d, de, tmat, ssst);
}



/**
  The function assembles the stiffness %matrix in the local
  coordinate system.

  @param d - resulting stiffness %matrix
  @param ssst - stress/strain state indicator

  @return The funtion returns resulting %matrix in the argument d.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::loc_matstiff (matrix &d, strastrestate /*ssst*/)
{
  copym(de, d);
}



/** 
  The function returns transformation %matrix from the local material coordinate system 
  to the global one in the argument. Components of the transformation %matrix
  are evaluated with respect to space coordinates of the given integration point ipp.
  The argument tmat is supposed to be allocated to dimensions 3x3.

  @param tmat - evaluated transformation %matrix 
  @param ipp - integration point id

  @return The function returns transformation %matrix in the argument tmat

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::give_transf_mat (matrix &tmat, long ipp)
{
  vector p(3);
  const char *namevar[] = {"x", "y", "z"};

  // x,y,z coordinates of the integration point will be stored in the vector p
  ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
  // evaluate transformation matrix
  t.evaluate(p, namevar, tmat);
}



/**
  The function computes stresses for the given strains.
  Resulting stresses are given in the global coordinate system.

  @param ipp - integration point id

  @return The funtion stores the resulting stresses in the stress array 
          of the given integration point.

  Created by Tomas Koudelka, 19.12.2014
*/
void elastortomat::nlstresses (long ipp)
{
  long n = Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector eps(ASTCKVEC(n)), sig(ASTCKVEC(n));
  matrix dg(ASTCKMAT(6,6)), d(ASTCKMAT(n,n));;
  
  // actual strains
  Mm->givestrain(0, ipp, eps);

  // assemble full stiffness matrix
  matstiff(dg, ipp, spacestress);

  // convert stiffness matrix(6,6) to reduced form depending on ssst
  tensor4_ematrix(d, dg, ssst);
  
  // compute sig = D.eps
  mxv(d, eps, sig);

  Mm->storestress(0, ipp, sig);

  // compute eps_z component in the case of plane-stress problem
  if (Mm->ip[ipp].ssst == planestress){
    eps(3) = -dg(2,0)/dg(2,2)*eps(0) - dg(2,1)/dg(2,2)*eps(1);
    Mm->storestrain(0, ipp, 3, 1, eps);
  }
}
