#include "elastisopdmat.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "mechmat.h"


/**
  Constructor initializes data members to zero or default values.

  Created by TKo, 01.2021
*/
elastisopdmat::elastisopdmat (void)
{
  //  Young's modulus of elasticity
  e = 0.0;
  //  Poisson's ratio
  nu = 0.0;
  // mean stress treshold
  p0 = 0.0;  
}



/**
  Destructor is defined only for the formal purposes.

  Created by TKo, 01.2021
*/
elastisopdmat::~elastisopdmat (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by TKo,
*/
void elastisopdmat::read (XFILE *in)
{
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  xfscanf (in,"%k%le %k%le","e",&e,"nu",&nu);
  // mean stress treshold
  xfscanf (in,"%k%le","p0",&p0);
  // read expression for pressure dependent modulus
  pf.read(in);
}



/**
  Function prints material parameters into the opened text file.
   
  @param out - pointer to the opened FILE

  @return The function does not return anything.

  Created by TKo, 01.2021
*/
void elastisopdmat::print (FILE *out)
{
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  fprintf (out,"%le %le %le ", e, nu, p0);
  pf.print(out);
}



/**
  The function initializes material model state variables, i.e. initial value of mean stress
*/
void elastisopdmat::initval(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;
  double i1s;
  vector epse(ASTCKVEC(ncompstr)), sig(ASTCKVEC(ncompstr)), sigt(ASTCKVEC(6));
  matrix d(ASTCKMAT(ncompstr, ncompstr));
 
  if (Mm->ip[ipp].eqother[ido] == 0.0){
    Mm->givestrain(0, ipp, epse);
    //  elastic stiffness matrix
    elmatstiff (d, ipp, ido);
    // initial value of stress
    mxv(d, epse, sig);
    // convert sigma to 6 component vector
    give_full_vector (sigt, sig, Mm->ip[ipp].ssst);
    // initial value of total mean stress
    i1s = first_invar(sigt)/3.0;
    Mm->ip[ipp].eqother[ido] = Mm->ip[ipp].other[ido] = i1s;

    Mm->ip[ipp].eqother[ido+1] = Mm->ip[ipp].other[ido+1] = give_actual_ym(ipp, ido);
    //fprintf (stderr,"\n Bod %ld ma invariant %le, modul %le",ipp,i1s, Mm->ip[ipp].eqother[ido+1]);
  }
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d[out]  - stiffness %matrix of material
  @param ipp[in] - integration point id
  @param ido[in] - index of internal variables for given material in the ipp other array

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK,
*/
void elastisopdmat::matstiff (matrix &d, long ipp, long ido)
{
  elmatstiff(d, ipp, ido);
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d[out]  - stiffness %matrix of material
  @param ipp[in] - integration point id
  @param ido[in] - index of internal variables for given material in the ipp other array

  @return The function returns material stiffness matrix in the parameter d.

  Created by JK,
*/
void elastisopdmat::elmatstiff (matrix &d, long ipp, long ido)
{
  strastrestate ssst = Mm->ip[ipp].ssst;  // retrieve stress state
  double e_c = give_actual_ym(ipp, ido);  // current Young modulus
//  fprintf (stderr,"\n E v ipp %ld je %le", ipp, e_c);

  switch (ssst){
    case bar:
      matstiff_bar (d, e_c);
      break;
    case plbeam:
      matstiff_plbeam (d, e_c);
      break;
    case spacebeam:
      matstiff_spacebeam (d, e_c);
      break;
    case planestress:
      matstiff_plstress (d, e_c);
      break;
    case planestrain:
      matstiff_plstrain (d, e_c);
      break;
    case platek:
      matstiff_platek (d, e_c);
      break;
    case plates:
      matstiff_plates (d, e_c);
      break;
    case axisymm:
      matstiff_axi (d, e_c);
      break;
    case spacestress:
      matstiff_spacestr (d, e_c);
      break;
    default:
      print_err("unknown number of components of stress tensor is required on element %ld, ip=%le", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
      abort();
  }
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for bar elements.
   
  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_bar (matrix &d, double e_c)
{
  d[0][0] = e_c;
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for plane beam elements.
   
  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_plbeam (matrix &d, double e_c)
{
  d[0][0] = e_c;
  d[1][1] = e_c/2.0/(1.0+nu);
  d[2][2] = e_c;
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for 
  plane beam elements.
   
  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_spacebeam (matrix &d, double e_c)
{
  d[0][0] = e_c;
  d[1][1] = e_c/2.0/(1.0+nu);
  d[2][2] = e_c/2.0/(1.0+nu);
  d[3][3] = e_c/2.0/(1.0+nu);
  d[4][4] = e_c;
  d[5][5] = e_c;
}



/**
  Function creates stiffness matrix of the elastic isotropic material for 
  2D problems (plane stress).

  @param d[out]  - stiffness matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_plstress (matrix &d, double e_c)
{
  double c;
  
  nullm(d);

  c = e_c/(1.0-nu*nu);
  
  d[0][0] = c;     d[0][1] = c*nu;  d[0][2] = 0.0;
  d[1][0] = c*nu;  d[1][1] = c;     d[1][2] = 0.0;
  d[2][0] = 0.0;   d[2][1] = 0.0;   d[2][2] = e_c/2.0/(1.0+nu);
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for 
  2D problems (plane strain).

  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_plstrain (matrix &d, double e_c)
{
  double c;
  
  nullm(d);

  c = e_c/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0] = c*(1.0-nu);   d[0][1] = c*nu;         d[0][2] = 0.0;
  d[1][0] = c*nu;         d[1][1] = c*(1.0-nu);   d[1][2] = 0.0;
  d[2][0] = 0.0;          d[2][1] = 0.0;          d[2][2] = e_c/2.0/(1.0+nu);

  if (d.m > 3){
    d[0][3] = d[0][1]; d[1][3] = d[1][0];
    d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for 
  2D axisymmetric problems.

  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_axi (matrix &d, double e_c)
{
  double g,s;
  
  nullm(d);
  
  g = e_c/2.0/(1.0+nu);
  s = e_c/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for plate 
  elements based on Kirchhoff theory.
   
  @param d[out]  - stiffness %matrix
  @param e_c[in] - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_platek (matrix &d, double e_c)
{
  double c,g;
  
  nullm(d);

  c = e_c/12.0/(1.0-nu*nu);
  g = e_c/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
}



/**
  Function creates stiffness %matrix of the elastic  isotropic material for plate 
  elements based on Mindlin-Reissner theory.
   
  @param d[out]  - stiffness %matrix
  @param e_c[in] - current Young modulus
   
  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_plates (matrix &d, double e_c)
{
  double c,g;
  
  nullm(d);

  c = e_c/12.0/(1.0-nu*nu);
  g = e_c/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
  
  d[3][3]=g;  d[4][4]=g;
}



/**
  Function creates stiffness %matrix of the elastic isotropic material for 3D problems.
   
  @param d[out]  - stiffness %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matstiff_spacestr (matrix &d, double e_c)
{
  double g, s;
  
  nullm(d);
  
  g = e_c/2.0/(1.0+nu);
  s = e_c/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;
}



/**
  Function assembles complience %matrix of material.
   
  @param c[out]  - complience %matrix of material
  @param ipp[in] - integration point id
  @param ido[in] - index of internal variables for given material in the ipp other array

  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl(matrix &c, long ipp, long ido)
{
  strastrestate ssst = Mm->ip[ipp].ssst; // detect stress state
  double e_c = give_actual_ym(ipp, ido);  // current Young modulus

  switch (ssst){
    case bar:
      matcompl_bar (c, e_c);
      break;
    case plbeam:
      matcompl_plbeam (c, e_c);
      break;
    case planestress:
      matcompl_plstress (c, e_c);
      break;
    case planestrain:
      matcompl_plstrain (c, e_c);
      break;
    case axisymm:
      matcompl_axi (c, e_c);
      break;
    case spacestress:
      matcompl_spacestr (c, e_c);
      break;
    default:
      print_err("unknown number of components of stress tensor is required on element %ld, ip=%le", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
      abort();
  }
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for bar elements.
   
  @param c[out]  - compliance %matrix of the material
  @param e_c[in] - current Young modulus
   
  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_bar (matrix &c, double e_c)
{
  c[0][0] = 1.0/e_c;
}



/**
  Function creates compliance matrix of the elastic
  isotropic material for plane beam elements
   
  @param c[out]  - compliance matrix of the material
  @param e_c[in] - current Young modulus
   
  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_plbeam (matrix &c, double e_c)
{
  c[0][0] = 1.0/e_c;
  c[1][1] = 2.0*(1.0+nu)/e_c;
  c[2][2] = 1.0/e_c;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane stress)

  @param c[out]  - compliance %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_plstress (matrix &c, double e_c)
{
  nullm(c);
  
  c[0][0] =  1.0/e_c;     c[0][1] = -1.0*nu/e_c;  c[0][2] = 0.0;
  c[1][0] = -1.0*nu/e_c;  c[1][1] =  1.0/e_c;     c[1][2] = 0.0;
  c[2][0] =  0.0;         c[2][1] =  0.0;         c[2][2] = 2.0*(1.0+nu)/e_c;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param c    - compliance %matrix of the material (output)
  @param e_c  - current Young modulus
  @param nu - current Poisson's ratio

  @return The function returns material complience %matrix in the parameter d.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_plstrain (matrix &c, double e_c)
{
  double g;
  
  nullm(c);
  
  g = (1.0+nu)/e_c;
  
  c[0][0] =  g*(1.0-nu);   c[0][1] = -1.0*g*nu;     c[0][2] = 0.0;
  c[1][0] = -1.0*g*nu;     c[1][1] =  g*(1.0-nu);   c[1][2] = 0.0;
  c[2][0] =  0.0;          c[2][1] =  0.0;          c[2][2] = 2.0*g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for axisymmetric problems
   
  @param c[out]  - compliance %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_axi (matrix &c, double e_c)
{
  double g;
  nullm(c);
  
  g = 2.0*(1.0+nu)/e_c;
  
  c[0][0]=1.0/e_c;     c[0][1]=-1.0*nu/e_c;  c[0][2]=c[0][1];    c[0][3]=0.0;
  c[1][0]=c[0][1];     c[1][1]=c[0][0];      c[1][2]=c[0][1];    c[1][3]=0.0;
  c[2][0]=c[0][1];     c[2][1]=c[0][1];      c[2][2]=c[0][0];    c[2][3]=0.0;

  c[3][0]=c[0][3];     c[3][1]=c[1][3];      c[3][2]=c[2][3];    c[3][3]=g;
}



/**
  Function creates compliance %matrix of the elastic
  isotropic material for 3D problems.
   
  @param c[out]  - compliance %matrix of the material
  @param e_c[in] - current Young modulus

  @return The function returns material complience %matrix in the parameter c.

  Created by TKo, 01.2021
*/
void elastisopdmat::matcompl_spacestr (matrix &c, double e_c)
{
  double g;
  nullm(c);
  
  g = 2.0*(1.0+nu)/e_c;
  
  c[0][0]=1.0/e_c;   c[0][1]=-1.0*nu/e_c;  c[0][2]=c[0][1];
  c[1][0]=c[0][1];   c[1][1]=c[0][0];      c[1][2]=c[0][1];
  c[2][0]=c[0][1];   c[2][1]=c[0][1];      c[2][2]=c[0][0];

  c[3][3]=g;         c[4][4]=g;            c[5][5]=g;
}



/**
  Function computes true stresses.
   
  @param ipp[in] - number of integration point
  @param ido[in] - index of internal variables for given material in the ipp other array
   
  @return The function does not return anything.

  Created by TKo, 01.2021
*/
void elastisopdmat::nlstresses (long ipp, long ido)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(n)), sig(ASTCKVEC(n)), epst(ASTCKVEC(6));
  matrix d(ASTCKMAT(n,n));
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff (d, ipp, ido);
  mxv (d, eps, sig);
  
  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (eps[0]+eps[1]);
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
  {
    for (i=0; i<n; i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }
  
  for (i=0; i<n; i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp[in] - integration point number in the mechmat ip array.
  @param ido[in] - index of internal variables for given material in the ipp other array

  Created by TKo, 01.2021
*/
void elastisopdmat::updateval (long ipp,long ido)
{
  Mm->ip[ipp].eqother[ido] = Mm->ip[ipp].other[ido]; // initial mean stress
}



/**
  The function returns actual value of elastic modulus with respect to
  initial mean stress value.

  @param ipp[in] - integration point number in the mechmat ip array.
  @param ido[in] - index of internal variables for given material in the ipp other array

  @return The function returns actual value of elastic modulus.

  Created by TKo, 01.2021  
*/
double elastisopdmat::give_actual_ym(long ipp,long ido)
{
  double p = Mm->ip[ipp].eqother[ido];
  double e_act = e; // set initial value of elastic modulus
  if (p < p0)
    e_act = pf.getval(p);

  return e_act;
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm[in] - selected material parameters (parameters which are changed)
  @param val[in] - array containing new values of parameters
   
  @return The function does not return anything.

  Created by TKo, 01.2021
*/
void elastisopdmat::changeparam (atsel &atm, vector &val)
{
  long i;
  
  for (i=0; i<atm.num; i++){
    switch (atm.atrib[i]){
      case 0:
        e=val[i];
        break;
      case 1:
	nu=val[i];
	break;
      case 2:
	p0=val[i];
	break;
      default:
	print_err("wrong number of atribute is required", __FILE__, __LINE__, __func__);
	abort();
    }
  }
}
