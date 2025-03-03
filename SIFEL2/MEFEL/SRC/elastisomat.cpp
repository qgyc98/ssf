#include "elastisomat.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
elastisomat::elastisomat (void)
{
  //  Young's modulus of elasticity
  e = 0.0;
  //  Poisson's ratio
  nu = 0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK,
*/
elastisomat::~elastisomat (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by JK,
*/
void elastisomat::read (XFILE *in)
{
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  xfscanf (in,"%k%lf %k%lf","e",&e,"nu",&nu);
}

/**
  Function prints material parameters into the opened text file.
   
  @param out - pointer to the opened FILE

  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void elastisomat::print (FILE *out)
{
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  fprintf (out,"%le %le",e,nu);
}


/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK,
*/
void elastisomat::matstiff (matrix &d,strastrestate ssst)
{
  switch (ssst){
  case bar:{
    matstiff_bar (d);
    break;
  }
  case plbeam:{
    matstiff_plbeam (d);
    break;
  }
  case spacebeam:{
    matstiff_spacebeam (d);
    break;
  }
  case planestress:{
    matstiff_plstress (d);
    break;
  }
  case planestrain:{
    matstiff_plstrain (d);
    break;
  }
  case platek:{
    matstiff_platek (d);
    break;
  }
  case plates:{
    matstiff_plates (d);
    break;
  }
  case axisymm:{
    matstiff_axi (d);
    break;
  }
  case spacestress:{
    matstiff_spacestr (d);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness matrix in the parameter d.

  Created by JK,
*/
void elastisomat::elmatstiff (matrix &d,strastrestate ssst)
{
  switch (ssst){
  case bar:{
    matstiff_bar (d);
    break;
  }
  case plbeam:{
    matstiff_plbeam (d);
    break;
  }
  case spacebeam:{
    matstiff_spacebeam (d);
    break;
  }
  case planestress:{
    matstiff_plstress (d);
    break;
  }
  case planestrain:{
    matstiff_plstrain (d);
    break;
  }
  case platek:{
    matstiff_platek (d);
    break;
  }
  case plates:{
    matstiff_plates (d);
    break;
  }
  case axisymm:{
    matstiff_axi (d);
    break;
  }
  case spacestress:{
    matstiff_spacestr (d);
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
   
  d - stiffness %matrix of the material (output)
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 11.9.2001
*/
void elastisomat::matstiff_bar (matrix &d)
{
  d[0][0] = e;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plane beam elements.
   
  @param d - stiffness %matrix of the material (output)
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 11.9.2001
*/
void elastisomat::matstiff_plbeam (matrix &d)
{
  d[0][0] = e;
  d[1][1] = e/2.0/(1.0+nu);
  d[2][2] = e;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plane beam elements.
   
  @param d - stiffness %matrix of the material
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 14.3.2003
*/
void elastisomat::matstiff_spacebeam (matrix &d)
{
  d[0][0] = e;
  d[1][1] = e/2.0/(1.0+nu);
  d[2][2] = e/2.0/(1.0+nu);
  d[3][3] = e/2.0/(1.0+nu);
  d[4][4] = e;
  d[5][5] = e;
}



/**
  Function creates stiffness matrix of the elastic
  isotropic material for 2D problems (plane stress).

  @param d - stiffness matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_plstress (matrix &d)
{
  double c;
  
  fillm(0.0,d);

  c = e/(1.0-nu*nu);
  
  d[0][0] = c;     d[0][1] = c*nu;  d[0][2] = 0.0;
  d[1][0] = c*nu;  d[1][1] = c;     d[1][2] = 0.0;
  d[2][0] = 0.0;   d[2][1] = 0.0;   d[2][2] = e/2.0/(1.0+nu);
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param d - stiffness %matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_plstrain (matrix &d)
{
  double c;
  
  fillm(0.0,d);

  c = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0] = c*(1.0-nu);   d[0][1] = c*nu;         d[0][2] = 0.0;
  d[1][0] = c*nu;         d[1][1] = c*(1.0-nu);   d[1][2] = 0.0;
  d[2][0] = 0.0;          d[2][1] = 0.0;          d[2][2] = e/2.0/(1.0+nu);

  if (d.m > 3)
    {
      d[0][3] = d[0][1]; d[1][3] = d[1][0];
      d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
    }
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D axisymmetric problems.

  @param d - stiffness %matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_axi (matrix &d)
{
  double g,s;
  
  fillm(0.0,d);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;

}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plate elements based on
  Kirchhoff theory.
   
  @param d - stiffness %matrix
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_platek (matrix &d)
{
  double c,g;
  
  fillm(0.0,d);

  c = e/12.0/(1.0-nu*nu);
  g = e/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
  
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for plate elements based on
  Mindlin-Reissner theory.
   
  @param d - stiffness %matrix
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_plates (matrix &d)
{
  double c,g;
  
  fillm(0.0,d);

  c = e/12.0/(1.0-nu*nu);
  g = e/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
  
  d[3][3]=g;  d[4][4]=g;
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 3D problems.
   
  @param d - stiffness %matrix of the material

  Created by JK, 19.7.2001
*/
void elastisomat::matstiff_spacestr (matrix &d)
{
  double g,s;
  
  fillm(0.0,d);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;

  ////////////////////////////////////////////////////////////
  // this printing below is only for debug:
  //long i,j;
  //fprintf(Out,"\n\n vysledna matice D^M\n");
  //for (i=0;i<6;i++){
  //for (j=0;j<6;j++)
  //  fprintf(Out,"%e  \n",d[i][j]);
  //fprintf(Out,"\n");
  //}
}




/**
  Function assembles complience %matrix of material.
   
  @param c - complience %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material complience %matrix in the parameter d.

  Created by JK,
*/
void elastisomat::matcompl (matrix &c,strastrestate ssst)
{
  switch (ssst){
  case bar:{
    matcompl_bar (c);
    break;
  }
  case plbeam:{
    matcompl_plbeam (c);
    break;
  }
  case planestress:{
    matcompl_plstress (c);
    break;
  }
  case planestrain:{
    matcompl_plstrain (c);
    break;
  }
  case axisymm:{
    matcompl_axi (c);
    break;
  }
  case spacestress:{
    matcompl_spacestr (c);
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
   
  @param c - compliance %matrix of the material (output)
   
  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 5.11.2002
*/
void elastisomat::matcompl_bar (matrix &c)
{
  c[0][0] = 1.0/e;
}



/**
  Function creates compliance matrix of the elastic
  isotropic material for plane beam elements
   
  @param c - compliance matrix of the material (output)
   
  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 5.11.2002
*/
void elastisomat::matcompl_plbeam (matrix &c)
{
  c[0][0] = 1.0/e;
  c[1][1] = 2.0*(1.0+nu)/e;
  c[2][2] = 1.0/e;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane stress)

  @param c - compliance %matrix of the material (output)

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 5.11.2002
*/
void elastisomat::matcompl_plstress (matrix &c)
{
  fillm(0.0,c);
  
  c[0][0] =  1.0/e;     c[0][1] = -1.0*nu/e;  c[0][2] = 0.0;
  c[1][0] = -1.0*nu/e;  c[1][1] =  1.0/e;     c[1][2] = 0.0;
  c[2][0] = 0.0;        c[2][1] = 0.0;        c[2][2] = 2.0*(1.0+nu)/e;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param c - compliance %matrix of the material (output)

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 5.11.2002
*/
void elastisomat::matcompl_plstrain (matrix &c)
{
  double g;
  
  fillm(0.0,c);
  
  g = (1.0+nu)/e;
  
  c[0][0] = g*(1.0-nu);   c[0][1] = -1.0*g*nu;    c[0][2] = 0.0;
  c[1][0] = -1.0*g*nu;    c[1][1] = g*(1.0-nu);   c[1][2] = 0.0;
  c[2][0] = 0.0;          c[2][1] = 0.0;          c[2][2] = 2.0*g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for axisymmetric problems
   
  @param c - compliance %matrix of the material (output)

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 8.4.2005
*/
void elastisomat::matcompl_axi (matrix &c)
{
  double g;
  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/e;
  
  c[0][0]=1.0/e;     c[0][1]=-1.0*nu/e;  c[0][2]=c[0][1];    c[0][3]=0.0;
  c[1][0]=c[0][1];   c[1][1]=c[0][0];    c[1][2]=c[0][1];    c[1][3]=0.0;
  c[2][0]=c[0][1];   c[2][1]=c[0][1];    c[2][2]=c[0][0];    c[2][3]=0.0;

  c[3][0]=c[0][3];   c[3][1]=c[1][3];    c[3][2]=c[2][3];    c[3][3]=g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 3D problems.
   
  @param c - compliance %matrix of the material (output)

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 5.11.2002
*/
void elastisomat::matcompl_spacestr (matrix &c)
{
  double g;
  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/e;
  
  c[0][0]=1.0/e;     c[0][1]=-1.0*nu/e;  c[0][2]=c[0][1];
  c[1][0]=c[0][1];   c[1][1]=c[0][0];    c[1][2]=c[0][1];
  c[2][0]=c[0][1];   c[2][1]=c[0][1];    c[2][2]=c[0][0];

  c[3][3]=g;         c[4][4]=g;          c[5][5]=g;
}


/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array

   16/04/2022 TKr
*/
void elastisomat::initval(long ipp, long /*im*/, long ido)
{
  Mm->ip[ipp].other[ido+0] = Mm->ip[ipp].eqother[ido+0] = 0.0;
  Mm->ip[ipp].other[ido+1] = Mm->ip[ipp].eqother[ido+1] = 0.0;
}


/**
  Function computes true stresses.
   
  @param ipp - number of integration point
   
  @return The function does not return anything.

  Created by JK,
*/
void elastisomat::nlstresses (long ipp, long /*im*/, long ido)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n),epst(ASTCKVEC(6));
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(n,n);
  double dtime = Mp->timecon.actualforwtimeincr(),epsv;
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff (d,ssst);
  mxv (d,eps,sig);
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (eps[0]+eps[1]);
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
  {
    for (i=0;i<n;i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }
  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig(i);
  }

  give_full_vector (epst,eps,Mm->ip[ipp].ssst);
  epsv = first_invar(epst);
  
  //debug??!! for 1D problem
  //epsv = eps[0];
  Mm->ip[ipp].other[ido+0] = epsv;
  // volumetric strain rate
  //Mm->ip[ipp].other[ido+1] = (epsv - Mm->ip[ipp].eqother[ido+0])/dtime;
  // the volumetric strain rate is computed via generalized trapesoidal rule
  if (dtime > 0.0)
    Mm->ip[ipp].other[ido+1] = 0.5*((epsv - Mm->ip[ipp].eqother[ido+0])/dtime + Mm->ip[ipp].eqother[ido+1]);
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

*/
void elastisomat::updateval (long ipp,long ido)
{
  Mm->ip[ipp].eqother[ido+0] = Mm->ip[ipp].other[ido+0];     //
  Mm->ip[ipp].eqother[ido+1] = Mm->ip[ipp].other[ido+1];     //
}



/**
  The function returns rate of the volumetric strain rate at the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by Tomas Koudelka 05.2018
*/
double elastisomat::give_strain_vol(long ipp, long ido)
{
  return Mm->ip[ipp].eqother[ido+1];
}


/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by JK,
*/
void elastisomat::changeparam (atsel &atm,vector &val)
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
