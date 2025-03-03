#include "varelastisomat.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"



/**
  Constructor initializes data members to zero or default values.

  Created by TKo, 02.2013
*/
varelastisomat::varelastisomat (void)
{
  //  Young's modulus of elasticity
  e = 0.0;
  //  Poisson's ratio
  nu = 0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by TKo, 02.2013
*/
varelastisomat::~varelastisomat (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by TKo, 02.2013
*/
void varelastisomat::read (XFILE *in)
{
  //  e - intial Young's modulus of elasticity
  //  nu - Poisson's ratio
  xfscanf (in,"%k%lf %k%lf", "e", &e, "nu", &nu);
}



/**
  Function prints material parameters to the opened text file.
   
  @param out - pointer to the opened FILE

  @return The function does not return anything.

  Created by TKo, 09.2020
*/
void varelastisomat::print (FILE *out)
{
  //  e - intial Young's modulus of elasticity
  //  nu - Poisson's ratio
  fprintf(out,"%le %le", e, nu);
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff (matrix &d, long ipp)
{
  elmatstiff(d, ipp);
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ipp - integration point pointer

  @return The function returns material stiffness matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::elmatstiff (matrix &d, long ipp)
{
  elmatstiff (d, ipp, Mm->ip[ipp].ssst);
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ipp - integration point pointer

  @return The function returns material stiffness matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::elmatstiff (matrix &d, long ipp, strastrestate ssst)
{
  double e_c = Mm->give_actual_ym(ipp);  // current Young modulus
  double nu_c = Mm->give_actual_nu(ipp);  // current Poisson's ratio

  switch (ssst){
  case bar:{
    matstiff_bar (d, e_c);
    break;
  }
  case plbeam:{
    matstiff_plbeam (d, e_c, nu_c);
    break;
  }
  case spacebeam:{
    matstiff_spacebeam (d, e_c, nu_c);
    break;
  }
  case planestress:{
    matstiff_plstress (d, e_c, nu_c);
    break;
  }
  case planestrain:{
    matstiff_plstrain (d, e_c, nu_c);
    break;
  }
  case platek:{
    matstiff_platek (d, e_c, nu_c);
    break;
  }
  case plates:{
    matstiff_plates (d, e_c, nu_c);
    break;
  }
  case axisymm:{
    matstiff_axi (d, e_c, nu_c);
    break;
  }
  case spacestress:{
    matstiff_spacestr (d, e_c, nu_c);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for bar elements.
   
  d   - stiffness %matrix of the material (output)
  e_c - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_bar (matrix &d, double e_c)
{
  d[0][0] = e_c;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plane beam elements.
   
  @param d    - stiffness %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_plbeam (matrix &d, double e_c, double nu_c)
{
  d[0][0] = e_c;
  d[1][1] = e_c/2.0/(1.0+nu_c);
  d[2][2] = e_c;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plane beam elements.
   
  @param d    - stiffness %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_spacebeam (matrix &d, double e_c, double nu_c)
{
  d[0][0] = e_c;
  d[1][1] = e_c/2.0/(1.0+nu_c);
  d[2][2] = e_c/2.0/(1.0+nu_c);
  d[3][3] = e_c/2.0/(1.0+nu_c);
  d[4][4] = e_c;
  d[5][5] = e_c;
}



/**
  Function creates stiffness matrix of the elastic
  isotropic material for 2D problems (plane stress).

  @param d - stiffness matrix of the material
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_plstress (matrix &d, double e_c, double nu_c)
{
  double c;
  
  nullm(d);

  c = e_c/(1.0-nu_c*nu_c);
  
  d[0][0] = c;       d[0][1] = c*nu_c;  d[0][2] = 0.0;
  d[1][0] = c*nu_c;  d[1][1] = c;       d[1][2] = 0.0;
  d[2][0] = 0.0;     d[2][1] = 0.0;     d[2][2] = e_c/2.0/(1.0+nu_c);
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param d    - stiffness %matrix of the material
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_plstrain (matrix &d, double e_c, double nu_c)
{
  double c;
  
  nullm(d);

  c = e_c/(1.0+nu_c)/(1.0-2.0*nu_c);
  
  d[0][0] = c*(1.0-nu_c);   d[0][1] = c*nu_c;         d[0][2] = 0.0;
  d[1][0] = c*nu_c;         d[1][1] = c*(1.0-nu_c);   d[1][2] = 0.0;
  d[2][0] = 0.0;            d[2][1] = 0.0;            d[2][2] = e_c/2.0/(1.0+nu_c);

  if (d.m > 3)
  {
    d[0][3] = d[0][1]; d[1][3] = d[1][0];
    d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D axisymmetric problems.

  @param d    - stiffness %matrix of the material
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_axi (matrix &d, double e_c, double nu_c)
{
  double g,s;
  
  nullm(d);
  
  g = e_c/2.0/(1.0+nu_c);
  s = e_c/(1.0+nu_c)/(1.0-2.0*nu_c);
  
  d[0][0]=s*(1-nu_c);  d[0][1]=s*nu_c;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];     d[1][1]=d[0][0];    d[1][2]=d[0][1];
  d[2][0]=d[0][1];     d[2][1]=d[0][1];    d[2][2]=d[0][0];

  d[3][3]=g;

}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plate elements based on
  Kirchhoff theory.
   
  @param d - stiffness %matrix
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_platek (matrix &d, double e_c, double nu_c)
{
  double c,g;
  
  nullm(d);

  c = e_c/12.0/(1.0-nu_c*nu_c);
  g = e_c/2.0/(1.0+nu_c);
  
  d[0][0]=c;        d[0][1]=c*nu_c;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;       d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;     d[2][2]=g/12.0;
  
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plate elements based on
  Mindlin-Reissner theory.
   
  @param d    - stiffness %matrix
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_plates (matrix &d, double e_c, double nu_c)
{
  double c,g;
  
  nullm(d);

  c = e_c/12.0/(1.0-nu_c*nu_c);
  g = e_c/2.0/(1.0+nu_c);
  
  d[0][0]=c;        d[0][1]=c*nu_c;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;       d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;     d[2][2]=g/12.0;
  
  d[3][3]=g;  d[4][4]=g;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 3D problems.
   
  @param d    - stiffness %matrix of the material
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  Created by TKo, 02.2013
*/
void varelastisomat::matstiff_spacestr (matrix &d, double e_c, double nu_c)
{
  double g,s;
  
  nullm(d);
  
  g = e_c/2.0/(1.0+nu_c);
  s = e_c/(1.0+nu_c)/(1.0-2.0*nu_c);
  
  d[0][0]=s*(1-nu_c);  d[0][1]=s*nu_c;     d[0][2]=s*nu_c;
  d[1][0]=d[0][1];     d[1][1]=d[0][0];    d[1][2]=d[0][1];
  d[2][0]=d[0][1];     d[2][1]=d[0][1];    d[2][2]=d[0][0];

  d[3][3]=g;           d[4][4]=g;          d[5][5]=g;
}




/**
  Function assembles complience %matrix of material.
   
  @param c   - complience %matrix of material (output)
  @param ipp - integration point pointer

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl (matrix &c, long ipp)
{
  strastrestate ssst = Mm->ip[ipp].ssst; // detect stress state
  double e_c = Mm->give_actual_ym(ipp);  // current Young modulus
  double nu_c = Mm->give_actual_nu(ipp);  // current Poisson's ratio

  switch (ssst){
  case bar:{
    matcompl_bar (c, e_c);
    break;
  }
  case plbeam:{
    matcompl_plbeam (c, e_c, nu_c);
    break;
  }
  case planestress:{
    matcompl_plstress (c, e_c, nu_c);
    break;
  }
  case planestrain:{
    matcompl_plstrain (c, e_c, nu_c);
    break;
  }
  case axisymm:{
    matcompl_axi (c, e_c, nu_c);
    break;
  }
  case spacestress:{
    matcompl_spacestr (c, e_c, nu_c);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for bar elements.
   
  @param c   - compliance %matrix of the material (output)
  @param e_c - current Young modulus
   
  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_bar (matrix &c, double e_c)
{
  c[0][0] = 1.0/e_c;
}



/**
  Function creates compliance matrix of the elastic
  isotropic material for plane beam elements
   
  @param c    - compliance matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio
   
  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_plbeam (matrix &c, double e_c, double nu_c)
{
  c[0][0] = 1.0/e_c;
  c[1][1] = 2.0*(1.0+nu_c)/e_c;
  c[2][2] = 1.0/e_c;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane stress)

  @param c    - compliance %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_plstress (matrix &c, double e_c, double nu_c)
{
  fillm(0.0,c);
  
  c[0][0] =  1.0/e_c;       c[0][1] = -1.0*nu_c/e_c;  c[0][2] = 0.0;
  c[1][0] = -1.0*nu_c/e_c;  c[1][1] =  1.0/e_c;       c[1][2] = 0.0;
  c[2][0] =  0.0;           c[2][1] =  0.0;           c[2][2] = 2.0*(1.0+nu_c)/e_c;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param c    - compliance %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_plstrain (matrix &c, double e_c, double nu_c)
{
  double g;
  
  fillm(0.0,c);
  
  g = (1.0+nu_c)/e_c;
  
  c[0][0] =  g*(1.0-nu_c);   c[0][1] = -1.0*g*nu_c;     c[0][2] = 0.0;
  c[1][0] = -1.0*g*nu_c;     c[1][1] =  g*(1.0-nu_c);   c[1][2] = 0.0;
  c[2][0] =  0.0;            c[2][1] =  0.0;            c[2][2] = 2.0*g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for axisymmetric problems
   
  @param c    - compliance %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_axi (matrix &c, double e_c, double nu_c)
{
  double g;
  fillm(0.0,c);
  
  g = 2.0*(1.0+nu_c)/e_c;
  
  c[0][0]=1.0/e_c;     c[0][1]=-1.0*nu_c/e_c;  c[0][2]=c[0][1];    c[0][3]=0.0;
  c[1][0]=c[0][1];     c[1][1]=c[0][0];        c[1][2]=c[0][1];    c[1][3]=0.0;
  c[2][0]=c[0][1];     c[2][1]=c[0][1];        c[2][2]=c[0][0];    c[2][3]=0.0;

  c[3][0]=c[0][3];     c[3][1]=c[1][3];        c[3][2]=c[2][3];    c[3][3]=g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 3D problems.
   
  @param c - compliance %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu_c - current Poisson's ratio

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 02.2013
*/
void varelastisomat::matcompl_spacestr (matrix &c, double e_c, double nu_c)
{
  double g;
  fillm(0.0,c);
  
  g = 2.0*(1.0+nu_c)/e_c;
  
  c[0][0]=1.0/e_c;   c[0][1]=-1.0*nu_c/e_c;  c[0][2]=c[0][1];
  c[1][0]=c[0][1];   c[1][1]=c[0][0];        c[1][2]=c[0][1];
  c[2][0]=c[0][1];   c[2][1]=c[0][1];        c[2][2]=c[0][0];

  c[3][3]=g;         c[4][4]=g;              c[5][5]=g;
}



/**
  Function computes true stresses.
   
  @param ipp - number of integration point
   
  @return The function does not return anything.

  Created by TKo, 02.2013
*/
void varelastisomat::nlstresses (long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  double nu_c = Mm->give_actual_nu(ipp);
  vector eps(n),sig(n);
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff (d,ssst);
  mxv (d,eps,sig);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu_c / (1.0 - nu_c) * (eps[0]+eps[1]);

  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
  
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by TKo, 02.2013
*/
void varelastisomat::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      e=val[i];
      break;
    }
    case 1:{
      nu=val[i];
      break;
    }
    default:{
      print_err("wrong number of atribute is required", __FILE__, __LINE__, __func__);
    }
    }
  }
}
