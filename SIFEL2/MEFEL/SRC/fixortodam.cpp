#include "fixortodam.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "sequent.h"
#include "intpoints.h"
#include "matrix.h"
#include "vecttens.h"


double fixortodam::L_MAX = 0.0;
double fixortodam::L_MIN = 1000.0;


fixortodam::fixortodam (void) : xeqtf(3), xeqcf(3)
{
  x1t = x2t = x3t = 0.0;
  x1c = x2c = x3c = 0.0;
}



fixortodam::~fixortodam (void)
{
}



/**
  Function reads material parameters from the opened text file.
   
  @param[in] in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2014
*/
void fixortodam::read (XFILE *in)
{
  xfscanf(in, "%k %le", "x1t", &x1t);
  xfscanf(in, "%k %le", "x2t", &x2t);
  xfscanf(in, "%k %le", "x3t", &x3t);

  xfscanf(in, "%k %le", "x1c", &x1c);
  xfscanf(in, "%k %le", "x2c", &x2c);
  xfscanf(in, "%k %le", "x3c", &x3c);

  xfscanf(in, "%k %le", "xeq1tf", &xeqtf(0));
  xfscanf(in, "%k %le", "xeq2tf", &xeqtf(1));
  xfscanf(in, "%k %le", "xeq3tf", &xeqtf(2));

  xfscanf(in, "%k %le", "xeq1cf", &xeqcf(0));
  xfscanf(in, "%k %le", "xeq2cf", &xeqcf(1));
  xfscanf(in, "%k %le", "xeq3cf", &xeqcf(2));
}



/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param[in] out - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 20.1.2015
*/
void fixortodam::print (FILE *out)
{
  fprintf(out, "%le ", x1t);
  fprintf(out, "%le ", x2t);
  fprintf(out, "%le ", x3t);

  fprintf(out, "%le ", x1c);
  fprintf(out, "%le ", x2c);
  fprintf(out, "%le ", x3c);

  fprintf(out, "%le ", xeqtf(0));
  fprintf(out, "%le ", xeqtf(1));
  fprintf(out, "%le ", xeqtf(2));

  fprintf(out, "%le ", xeqcf(0));
  fprintf(out, "%le ", xeqcf(1));
  fprintf(out, "%le ", xeqcf(2));
}



/**
   This function initializes material model data

   @param[in] ipp - integration point number
   @param[in] ido - index of internal variables for given material in the ipp other array

   12/06/2012 TKo
*/
void fixortodam::initval(long ipp, long ido)
{

  matrix d(ASTCKMAT(6,6));

  // It is supposed that elastic stiffness matrix is taken from the elastic orthotropic material 
  // which is given in the local coord. system of the material
  Mm->loc_elmatstiff(d, ipp);

  for (long i = 0; i < 3; i++)
  {
    Mm->ip[ipp].eqother[ido + 3 + i] = 2 * xeqtf(i);
    Mm->ip[ipp].eqother[ido + 6 + i] = xeqtf(i);
    Mm->ip[ipp].eqother[ido + 9 + i] = xeqcf(i);
  }

  if (L_MAX<pow(Mt->give_volume(Mm->elip[ipp]), 1.0 / 3.0))
    L_MAX = pow(Mt->give_volume(Mm->elip[ipp]), 1.0 / 3.0);

  if (L_MIN>pow(Mt->give_volume(Mm->elip[ipp]), 1.0 / 3.0))
    L_MIN = pow(Mt->give_volume(Mm->elip[ipp]), 1.0 / 3.0);
}



/**
  The function assembles the stiffness %matrix in the global
  coordinate system.

  @param[in]  ipp  - integration point id
  @param[in] ido - index of internal variables for given material in the ipp other array
  @param[out] d    - resulting stiffness %matrix in the reduced format for element computations

  @return The funtion returns resulting %matrix in the argument d.

*/
void fixortodam::matstiff (long ipp,long ido,matrix &d)
{
  matrix dl(ASTCKMAT(6,6)), dg(ASTCKMAT(6,6));
  matrix tmat(ASTCKMAT(3,3));

  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      // It is supposed that elastic stiffness matrix is taken from the elastic orthotropic material 
      // which is given in the local coord. system of the material
      Mm->loc_elmatstiff (dl,ipp);
      break;
    case tangent_stiff:
      tmatstiff(ipp, ido, dl);
      break;
    default:
      print_err("unknown type of stifness matrix is required.", __FILE__, __LINE__, __func__);
  }

  // evaluate transformation matrix
  Mm->loc_transf_mat(tmat, ipp);

  // transformation to the global coordinate system and conversion 
  // to the reduced form according to stress/straint state
  lg_tens4transf(dg, dl, tmat, Mm->ip[ipp].ssst);
  // convert stiffness matrix(6,6) to reduced form depending on ssst
  tensor4_ematrix(d, dg, Mm->ip[ipp].ssst);
}



/**
  The function computes tangent material stiffnes %matrix in the material local coordinate system but
  actually, only secant matrix is available.

  @param[in]  ipp - integration point number
  @param[in]  ido - index of internal variables for given material in the ipp other array
  @param[out] d   - allocated matrix structure for material stiffness %matrix in the full format (6x6) 
  @return The function returns elastic stiffness %matrix for principal stress computation 
          in the parameter d.
*/
void fixortodam::tmatstiff (long ipp, long ido, matrix &d)
{
  long i;
  vector omega(ASTCKVEC(3));

  // attained values of damage parameters from the last equilibrium state
  for (i=0; i<3; i++)
    omega(i) = Mm->ip[ipp].eqother[ido+i];

  /*  
  vector xeq0(ASTCKVEC(3));
  for (i = 0; i<3; i++)
  {
    xeq0(i)  = Mm->ip[ipp].eqother[ido+3+i];
    xeqtf(i) = Mm->ip[ipp].eqother[ido + 6 + i];
    xeqcf(i) = Mm->ip[ipp].eqother[ido + 9 + i];
  }
  */
  
  secstiffmat(ipp, omega, d);
}



/**
  The function returns actual values of equivalent displacements
  in the xeq.

  @param[in]  epst - strain tensor in the local coordinate system of the material
  @param[out] xeq  - equivalent displacement %vector (output parameter)

*/
void fixortodam::compute_eqdispl(long ipp, matrix &epst, vector &xeq)
{   
  long i;
  double l;

  long eid = Mm->elip[ipp];
  switch (Mt->give_dimension(eid))
  {
    case 1:
      l = Mt->give_length(eid);
      break;
    case 2:
      l = sqrt(Mt->give_area(eid));
      break;
    case 3:
      l = pow(Mt->give_volume(eid), 1.0/3.0);
      break;
    default:
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
      abort();
  }

  // It is supposed that elastic stiffness matrix is taken from the elastic orthotropic material 
  // which is given in the local coord. system of the material

  for (i = 0; i<3; i++)
    xeq(i) = l*(epst(i, i));
    
}



/**
  The function returns actual values of damage parameters in the particular
  local material directions stored in the %vector omega

  @param[in]  ipp  - integration point id
  @param[in]  xeq  - equivalent displacement %vector
  @param[in]  xeq0 - equivalent displacement treshold
  @param[in]  omegao - %vector of damage parametrs from the last equilibrium state
  @param[out] omega  - %vector of calculated damage parametrs

*/
void fixortodam::compute_dam(long ipp, vector &xeq, vector &xeq0, vector &omegao, vector &omega)
{
  long i;
  
  long eid = Mm->elip[ipp];
  for (i = 0; i < 3; i++)
  {
    if (xeq0(i) == 2 * xeqtf(i) || fabs(xeq(i)) <= fabs(xeq0(i)))
      omega(i) = 0.0;
    else
    {
      if ((xeq(i) > 0.0 && xeqtf(i) < xeq0(i)) || (xeq(i) < 0.0 && xeqcf(i) > xeq0(i)))
      {
        print_err("fracture energy of the element number %ld is too small. Step number = %ld", __FILE__, __LINE__, __func__, eid, Mp->istep);
        abort();
      }

      if (xeq(i) > 0.0 && xeq(i) >= xeqtf(i))
      {
        omega(i) = 1.0;
        continue;
      }

      if (xeq(i) < 0.0 && xeq(i) <= xeqcf(i))
      {
        omega(i) = 1.0;
        continue;
      }

      if (xeq(i) > 0.0)
        omega(i) = xeqtf(i)*(xeq(i) - xeq0(i)) / (xeq(i)*(xeqtf(i) - xeq0(i)));
      if (xeq(i) < 0.0)
        omega(i) = xeqcf(i)*(xeq(i) - xeq0(i)) / (xeq(i)*(xeqcf(i) - xeq0(i)));
    }
  }
  for (i = 0; i < 3; i++)
  {
    if (omega(i) < omegao(i))
      omega(i) = omegao(i);
    if (omega(i) < 0.0)
      omega(i) = 0.0;
    if (omega(i) > 1.0)
      omega(i) = 1.0;
  }
}



/**
  The function returns actual values of damage parameters in the particular
  local material directions stored in the %vector omega


  @param[in]  omega - %vector of damage parameters 
  @param[in]  ipp   - integration point id
  @param[out] d     - secant stiffness %matrix in the full format (6x6)

*/
void fixortodam::secstiffmat(long ipp, vector &omega, matrix &d)
{ 
  // It is supposed that elastic stiffness matrix is taken from the elastic orthotropic material 
  // which is given in the local coord. system of the material
  Mm->loc_elmatstiff (d,ipp);

  d(0,0)=pow(1-omega(0),2)*d(0,0);
  d(0,1)=(1-omega(0))*(1-omega(1))*d(0,1);
  d(0,2)=(1-omega(0))*(1-omega(2))*d(0,2);
  d(1,0)=d(0,1);
  d(1,1)=pow(1-omega(1),2)*d(1,1);
  d(1,2)=(1-omega(1))*(1-omega(2))*d(1,2);
  d(2,0)=d(0,2);
  d(2,1)=d(1,2);
  d(2,2)=pow(1-omega(2),2)*d(2,2);
  if (omega(0) == 1.0 && omega(1) == 1.0)
    d(3, 3) = 0;
  else
    d(3, 3) = pow((2 * (1 - omega(0))*(1 - omega(1))) / (2 - omega(0) - omega(1)), 2)*d(3,3);
  if (omega(0) == 1.0 && omega(2) == 1.0)
    d(4, 4) = 0;
  else
    d(4, 4) = pow((2 * (1 - omega(0))*(1 - omega(2))) / (2 - omega(0) - omega(2)), 2)*d(4,4);
  if (omega(1) == 1.0 && omega(2) == 1.0)
    d(5, 5) = 0;
  else
    d(5, 5) = pow((2 * (1 - omega(1))*(1 - omega(2))) / (2 - omega(1) - omega(2)), 2)*d(5,5);
}



/**
  The function computes actual stresses in the integration point and stores
  them into ip stress array.

  @param[in] ipp - integration point pointer
  @param[in] im  - index of material type for given ip
  @param[in] ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything but actualizes stress array of ipp int. point.
*/
void fixortodam::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncompstr=Mm->ip[ipp].ncompstr; // number of stress/strain components in vector form

  vector omega(ASTCKVEC(3)), omegao(ASTCKVEC(3)), xeq(ASTCKVEC(3));
  vector epsg(ASTCKVEC(ncompstr)), eps(ASTCKVEC(ncompstr)), sig(ASTCKVEC(ncompstr)), sig_g(ASTCKVEC(ncompstr));
  matrix epst(ASTCKMAT(3,3)), sigt(ASTCKMAT(3,3)), d(ASTCKMAT(6,6));
  matrix tmat(ASTCKMAT(3,3));

  /** equivalent displacement treshold in longitudinal, the first transversal and the second 
      transversal directions (tension)*/
  vector xeq0(3);

  // actual strains in the global coord. system
  for (i=0; i<ncompstr; i++)
    epsg(i) = Mm->ip[ipp].strain[i];

  // assemble transformation matrix(3,3) for the given integration point
  Mm->loc_transf_mat(tmat, ipp);

  //  transformation of strains given in the global coordinate system to the material local coord. system.
  vector_tensor (epsg,epst,Mm->ip[ipp].ssst,strain);
  glmatrixtransf(epst, tmat);  
  tensor_vector(eps,epst,Mm->ip[ipp].ssst,strain);

  // eps contains vector representation of strains in local coord. system of the material
  // epst contains tensor representation of strains in local coord. system of the material

  // attained values of damage parameters from the last equilibrium state
  for (i=0; i<3; i++)
  {
    omegao(i) = Mm->ip[ipp].eqother[ido + i];
    xeq0(i)   = Mm->ip[ipp].eqother[ido + 3 + i];    
    xeqtf(i)  = Mm->ip[ipp].eqother[ido + 6 + i];
    xeqcf(i)  = Mm->ip[ipp].eqother[ido + 9 + i];
  }

  // compute equivalent displacements
  compute_eqdispl(ipp, epst, xeq);

  // compute actual damage parameters
  compute_dam(ipp, xeq, xeq0, omegao, omega);

  // secant stiffnes matrix
  secstiffmat(ipp, omega, d);
  
  // compute actual stress vector in the local coord. system of the material
  // sig = D(omega) * eps
  mxv(d, eps, sig);
  for (i = 0; i < 3; i++)
  {
    if (eps(i) > 0.0 && omega(i) > 0.0 && sig(i) < 0.0)
      omega(i) = 1.0;
    if (eps(i) < 0.0 && omega(i) > 0.0 && sig(i) > 0.0)
      omega(i) = 1.0;
  }

  double l;
  long eid = Mm->elip[ipp];
  switch (Mt->give_dimension(eid))
  {
    case 1:
      l = Mt->give_length(eid);
      break;
    case 2:
      l = sqrt(Mt->give_area(eid));
      break;
    case 3:
      l = pow(Mt->give_volume(eid), 1.0 / 3.0);
      break;
    default:
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
      abort();
  }
  
  if (xeq0[0] == 2 * xeqtf(0)) 
  {
    if ((eps[0] > 0.0 && sig[0] >= x1t))
    {
      xeq0(0) = l * eps[0];
      xeqtf(0)=xeqtf(0)*x1t/sig[0];
      Mm->ip[ipp].other[ido + 6] = xeqtf(0);
    }
    if ((eps[0] < 0.0 && sig[0] <= x1c))
    {
      xeq0(0) = l * eps[0];
      xeqcf(0)=xeqcf(0)*x1c/sig[0];
      Mm->ip[ipp].other[ido + 9] = xeqcf(0);
    }
    if (xeq0[0] != 2*xeqtf(0))
    {
      Mm->ip[ipp].other[ido + 3] = xeq0(0);
      Mm->ip[ipp].other[ido + 12] = eps(0);
      Mm->ip[ipp].other[ido + 15] = sig(0);
    }
  }

  if (xeq0[1] == 2 * xeqtf(1)) 
  {
    if ((eps[1] > 0.0 && sig[1] >= x2t))
    {
      xeq0(1) = l * eps[1];
      xeqtf(1) = xeqtf(1)*x1t / sig[1];
      Mm->ip[ipp].other[ido + 6] = xeqtf(1);
    }
    if ((eps[1] < 0.0 && sig[1] <= x2c))
    {
      xeq0(1) = l * eps[1];
      xeqcf(1) = xeqcf(1)*x1c / sig[1];
      Mm->ip[ipp].other[ido + 9] = xeqcf(1);
    }
    if (xeq0[1] != 2 * xeqtf(1))
    {
      Mm->ip[ipp].other[ido + 4] = xeq0(1);
      Mm->ip[ipp].other[ido + 13] = eps(1);
      Mm->ip[ipp].other[ido + 16] = sig(1);
    }
  }

  if (xeq0[2] == 2 * xeqtf(2))
  {
    if ((eps[2] > 0.0 && sig[2] >= x3t))
    {
      xeq0(2) = l * eps[2];
      xeqtf(2) = xeqtf(2)*x1t / sig[2];
      Mm->ip[ipp].other[ido + 6] = xeqtf(2);
    }
    if ((eps[2] < 0.0 && sig[2] <= x3c))
    {
      xeq0(2) = l * eps[2];
      xeqcf(2) = xeqcf(2)*x1c / sig[2];
      Mm->ip[ipp].other[ido + 9] = xeqcf(2);
    }
    if (xeq0[2] != 2 * xeqtf(2))
    {
      Mm->ip[ipp].other[ido + 5] = xeq0(2);
      Mm->ip[ipp].other[ido + 14] = eps(2);
      Mm->ip[ipp].other[ido + 17] = sig(2);
    }
  }

  // transforms stress vector in the local coord. system of the material to global coord. system
  vector_tensor (sig,sigt,Mm->ip[ipp].ssst,stress);
  lgmatrixtransf(sigt, tmat);  
  tensor_vector(sig_g,sigt,Mm->ip[ipp].ssst,stress);

  // store actual stress in the global coord. system to the integration point
  for (i=0; i<ncompstr; i++)
    Mm->ip[ipp].stress[i] = sig_g[i];

  
  // store actual damage parameters to the integration point
  for (i = 0; i < 3; i++)
    Mm->ip[ipp].other[ido + i] = omega(i);
}



/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param[in] ipp - integration point number in the mechmat ip array.
  @param[in] im  - index of material type for given ip
  @param[in] ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

*/
void fixortodam::updateval (long ipp, long /*im*/, long ido)
{
  long i;

  for (i=0; i<18; i++)
    Mm->ip[ipp].eqother[ido+i] = Mm->ip[ipp].other[ido+i]; 
}
