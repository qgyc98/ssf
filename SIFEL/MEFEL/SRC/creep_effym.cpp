#include "creep_effym.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>

/**
  This constructor inializes attributes to zero values.
*/
creep_effym::creep_effym (void)
{
  evf=b3law;
  e0 = 0.0;
  nu = 0.0;

  //b3:
  q1=q2=q3=q4=q5=0.0;
  type_e = 0;
  e28=30.0e9; //Pa
  fc=35.8e6; //Pa
  wc=0.43;
  sc=3.4;
  gc=1.98;
  cs=305.0;  //kg*m^-3
  a1=1.05;
  a2=1.2;
  kd=0.15;   //m
  tb_time = 0.0;//s 

  //double power law - je nutne opravit ??!!
  qs = 0.0;  psi = 0.0; n = 0.0; m = 0.0; alpha = 0.0;  
  t0 = 0.0;
}



/**
  This destructor is only for the formal purposes.
*/
creep_effym::~creep_effym (void)
{

}



/**
   Function reads material parameters.
   
   @param in - input stream

   7.2008 TKo, TKr
*/
void creep_effym::read (XFILE *in)
{
  xfscanf (in,"%k%lf","nu",&nu);
  xfscanf (in, "%k%m", "ym_evolfunc", &ym_evolfunc_kwdset, &evf);
  switch (evf)
  {
    case b3law:
      xfscanf(in, "%k%ld", "type_e", &type_e);
      
      if (type_e == 1){xfscanf (in,"%lf ",&e28);}
      //input in Pa
      e0 = e28;
      e28 = e28/6.89476*1.0e-3;//to psi
      
      xfscanf (in,"%k%lf %k%lf %k%lf %k%lf %k%lf %k%lf %k%lf %k%lf %k%lf %k%lf", 
               "fc", &fc, "wc", &wc, "sc", &sc, "gc", &gc, "cs", &cs, "a1", &a1, 
               "a2", &a2, "ks", &ks, "kd", &kd, "tb_time", &tb_time);
      
      //input in Pa
      fc = fc/6.89476*1.0e-3;//to psi
      //input in kg*m^-3
      cs = cs/16.03;//to lb*ft^-3
      //input in m
      kd=kd*1000.0/25.4;//to inch

      //pro vypocet zacinajici od casu 0.0
      tl = tb_time/86400;

      if(tl <= 0.0){
	print_err("Age of concrte is less than zero", __FILE__, __LINE__, __func__);
	exit(0);
      }
    
      break;
    case doublepwrlaw:
      xfscanf(in, "%k%lf %k%lf %k%lf %k%lf %k%lf %k%lf %k%lf", "e0", &e0, "t0", &t0, "qs", &qs, "psi", &psi, "m", &m, "alpha", &alpha, "n", &n);
      break;
    case userfunc:
      xfscanf(in, "%k%lf", "e0", &e0);
      gf.read(in);
      break;
    default:
      print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
  }
}



/**
   Function prints material parameters.
   
   @param out - output stream

   6.2016 TKo
*/
void creep_effym::print (FILE *out)
{
  fprintf (out,"%le ", nu);
  fprintf (out, "%d ", evf);
  switch (evf)
  {
    case b3law:
      fprintf(out, "%ld ", type_e);
      
      if (type_e == 1){fprintf (out,"%le ", e28*6.89476*1.0e-3);}
      
      fprintf (out, "%le %le %le %le %le\n     %le %le %le %le %le\n", 
                      fc*6.89476*1.0e-3, wc, sc, gc, cs*16.03, 
                      a1, a2, ks, kd*25.4/1000.0, tb_time);
      
      //input in Pa
      fc = fc/6.89476*1.0e-3;//to psi
      //input in kg*m^-3
      cs = cs/16.03;//to lb*ft^-3
      //input in m
      kd=kd*1000.0/25.4;//to inch

      //pro vypocet zacinajici od casu 0.0
      tl = tb_time/86400;

      if(tl <= 0.0){
	print_err("Age of concrte is less than zero", __FILE__, __LINE__, __func__);
	exit(0);
      }
    
      break;
    case doublepwrlaw:
      fprintf(out, "%le %le %le %le %le %le %le\n", e0, t0, qs, psi, m, alpha, n);
      break;
    case userfunc:
      fprintf(out, "%le ", e0);
      gf.print(out);
      break;
    default:
      print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
  }
}



/**
  Function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array


*/
void creep_effym::matstiff (matrix &d,long ipp,long /*im*/, long /*ido*/)
{
  matstiff(d, Mm->ip[ipp].ssst);
}



/**
  Function computes increment of stresses due to ageing in the integration point and stores
  them into ip stress array. 
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp other array

  6.2016 TKo
*/
void creep_effym::nlstressesincr (long ipp, long im, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  double e, ep, nu;
  vector epsn(ncomp),eps_ag(ncomp);


  if (Mm->ip[ipp].ssst == planestress)
  {
    nu = Mm->give_actual_nu(ipp, im, ido);
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
  }
  // Young's modulus E_p = E(t-dt) from the previous time step
  ep = Mm->ip[ipp].eqother[ido+0];

  // actual Young's modulus E(t)
  e = Mm->give_actual_ym(ipp, im, ido);

  // store actual Young modulus
  Mm->ip[ipp].other[ido] = e;

  //  initial values of total strains
  for (i=0;i<ncomp;i++)
  {
    epsn[i]   = Mm->ip[ipp].strain[i];
    eps_ag[i] = Mm->ip[ipp].other[ido+1+i];
  }
  
  //  computation of strains due to ageing from the effective Young modulus
  for (i=0;i<ncomp;i++)
  {
    eps_ag[i] = (eps_ag[i]*ep + (e-ep)*epsn[i])/e;
    Mm->ip[ipp].other[ido+1+i] = eps_ag[i];
  }
}



/**
  Function computes correct stresses in the integration point and stores
  them into ip stress array. 
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp other array

  7.2008 TKo
*/
void creep_effym::nlstresses (long ipp, long /*im*/, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  vector epsn(ncomp),eps_ag(ncomp),sigma(ncomp);
  matrix d(ncomp,ncomp);



  //  initial values of total strains
  for (i=0;i<ncomp;i++)
    epsn[i] = Mm->ip[ipp].strain[i];
  
  Mm->matstiff(d,ipp);
  
  // restore  "irreversible ageing strains" from the effective Young modulus calculated in nlstressincr
  for (i=0;i<ncomp;i++)
    eps_ag[i] = Mm->ip[ipp].other[ido+1+i];

  // compute actual stresses sig = D_e(t)*(eps-eps_{ag})
  subv(epsn, eps_ag, eps_ag);
  mxv(d, eps_ag, sigma);   
  //  new stress data storage
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].stress[i]=sigma[i];
}



/**
  Function creates stiffness %matrix depending on stress-strain state and 
  type of required stiffness %matrix in problem description Mp->nlman->stmat 
  (initial, tangent,...).

  @param d    - stiffness %matrix (output)
  @param ssst - stress/strain state indicator

  7.2008 TKo
*/
void creep_effym::matstiff (matrix &d,strastrestate ssst)
{
  double ea; // actual Young's modulus

  if (Mp->nlman->stmat == initial_stiff)
    ea  = e0;
  else
    ea = actual_modulus();
  matstiff(d, ssst, ea);
}



/**
  Function creates stiffness %matrix depending on stress-strain state and 
  actual value of Young's modulus E_a.

  @param d    - stiffness %matrix (output)
  @param ssst - stress/strain state indicator
  @param ea   - actual value of Young's modulus

  7.2008 TKo
*/
void creep_effym::matstiff (matrix &d,strastrestate ssst,  double ea)
{
  switch (ssst)
  {
    case bar:
      matstiff_bar (d,ea);
      break;
    case plbeam:
      matstiff_plbeam (d,ea);
      break;
    case spacebeam:
      matstiff_spacebeam (d,ea);
      break;
    case planestress:
     matstiff_plstress (d,ea);
     break;
    case planestrain:
      matstiff_plstrain (d,ea);
      break;
    case platek:
      matstiff_plate (d,ea);
      break;
    case plates:
      matstiff_plate (d,ea);
      break;
    case axisymm:
      matstiff_axi (d,ea);
      break;
    case spacestress:
      matstiff_spacestr (d,ea);
      break;
    default:
      print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function creates stiffness %matrix depending on stress-strain state and 
  initial value of Young modulus E_0.

  @param d    - stiffness %matrix (output)
  @param ssst - stress/strain state indicator

  7.2008 TKo
*/
void creep_effym::elmatstiff (matrix &d,strastrestate ssst)
{
  matstiff (d, ssst, e0);
}



/**
   Function creates stiffness %matrix for bar elements.
   
   @param d  - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_bar (matrix &d, double ea)
{
  d[0][0] = ea;
}



/**
   Function creates stiffness %matrix for plane beam elements.
   
   @param d   - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_plbeam (matrix &d, double ea)
{
  fillm(0.0,d);

  d[0][0] = ea;
  d[1][1] = ea/2.0/(1.0+nu);
  d[2][2] = ea;
}



/**
   Function creates stiffness %matrix for space beam elements.
   
   @param d  - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_spacebeam (matrix &d, double ea)
{
  fillm(0.0,d);

  d[0][0] = ea;
  d[1][1] = ea/2.0/(1.0+nu);
  d[2][2] = ea/2.0/(1.0+nu);
  d[3][3] = ea/2.0/(1.0+nu);
  d[4][4] = ea;
  d[5][5] = ea;
}



/**
   Function creates stiffness %matrix for 2D problems (plane stress).

   @param d   - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_plstress (matrix &d, double ea)
{
  double c;
  
  fillm(0.0,d);

  c = ea/(1.0-nu*nu);
  
  d[0][0] = c;     d[0][1] = c*nu;  d[0][2] = 0.0;
  d[1][0] = c*nu;  d[1][1] = c;     d[1][2] = 0.0;
  d[2][0] = 0.0;   d[2][1] = 0.0;   d[2][2] = ea/2.0/(1.0+nu);
}



/**
   Function creates stiffness %matrix for 2D problems (plane strain).
 
   @param d   - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_plstrain (matrix &d, double ea)
{
  double c;

  fillm(0.0,d);
  c = ea/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0] = c*(1.0-nu);   d[0][1] = c*nu;         d[0][2] = 0.0;
  d[1][0] = c*nu;         d[1][1] = c*(1.0-nu);   d[1][2] = 0.0;
  d[2][0] = 0.0;          d[2][1] = 0.0;          d[2][2] = ea/2.0/(1.0+nu);

  if (d.m > 3)
  {
    d[0][3] = d[0][1]; d[1][3] = d[1][0];
    d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
}



/**
   Function creates stiffness %matrix for 2D problems (axisymmetric problem).

   @param d  - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_axi (matrix &d, double ea)
{
  double g,s;
 
  fillm(0.0,d);
  
  g = ea/2.0/(1.0+nu);
  s = ea/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;

}



/**
   Function creates stiffness %matrix for plate elements.
   
   @param d  - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_plate (matrix &d, double ea)
{
  double c,g;

  fillm(0.0,d);

  c = ea/12.0/(1.0-nu*nu);
  g = ea/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
  
  d[3][3]=g;  d[4][4]=g;
}



/**
   Function creates stiffness %matrix for 3D problems.
   
   @param d  - stiffness %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   TKo 7.2008
*/
void creep_effym::matstiff_spacestr (matrix &d, double ea)
{
  double g,s;
  
  fillm(0.0,d);
  
  g = ea/2.0/(1.0+nu);
  s = ea/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;
}



/**
  Function creates complience %matrix depending on stress-strain state and 
  type of required stiffness %matrix in problem description Mp->nlman->stmat (initial, tangent,...).

  @param c    - compliance %matrix of the material (output)
  @param ssst - stress/strain state indicator

  7.2008 TKo
*/
void creep_effym::matcompl (matrix &c,strastrestate ssst)
{
  double ea;

  if (Mp->nlman->stmat == initial_stiff)
    ea  = e0;
  else
    ea = actual_modulus();

  switch (ssst)
  {
    case bar:
      matcompl_bar (c,ea);
      break;
    case plbeam:
      matcompl_plbeam (c,ea);
      break;
    case planestress:
      matcompl_plstress (c,ea);
      break;
    case planestrain:
      matcompl_plstrain (c,ea);
      break;
    case axisymm:
      matcompl_axi (c,ea);
      break;
    case spacestress:
      matcompl_spacestr (c,ea);
      break;
    default:
      print_err ("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
}



/**
   Function creates compliance %matrix for bar elements.
   
   @param c  - compliance %matrix of the material (output)
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_bar (matrix &c, double ea)
{
  c[0][0] = 1.0/ea;
}



/**
   Function creates compliance %matrix for plane beam elements.
   
   @param c  - compliance %matrix of the material
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_plbeam (matrix &c, double ea)
{
  fillm(0.0,c);

  c[0][0] = 1.0/ea;
  c[1][1] = 2.0*(1.0+nu)/ea;
  c[2][2] = 1.0/ea;
}



/**
   Function creates compliance %matrix for 2D problems (plane stress).

   @param c  - compliance %matrix of the material
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_plstress (matrix &c, double ea)
{
  fillm(0.0,c);
  
  c[0][0] =  1.0/ea;     c[0][1] = -1.0*nu/ea;  c[0][2] = 0.0;
  c[1][0] = -1.0*nu/ea;  c[1][1] =  1.0/ea;     c[1][2] = 0.0;
  c[2][0] =  0.0;        c[2][1] =  0.0;        c[2][2] = 2.0*(1.0+nu)/ea;
}



/**
   Function creates compliance %matrix for 2D problems (plane strain).

   @param c  - compliance %matrix of the material
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_plstrain (matrix &c, double ea)
{
  double g;

  fillm(0.0,c);
  
  g = (1.0+nu)/ea;
  
  c[0][0] =  g*(1.0-nu);  c[0][1] = -1.0*g*nu;     c[0][2] = 0.0;
  c[1][0] = -1.0*g*nu;    c[1][1] =  g*(1.0-nu);   c[1][2] = 0.0;
  c[2][0] =  0.0;         c[2][1] =  0.0;          c[2][2] = 2.0*g;
}



/**
   Function creates compliance %matrix for axisymmetric problems
   
   @param c  - compliance %matrix of the material
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_axi (matrix &c, double ea)
{
  double g;

  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/ea;
  
  c[0][0]=1.0/ea;    c[0][1]=-1.0*nu/ea;  c[0][2]=c[0][1];    c[0][3]=0.0;
  c[1][0]=c[0][1];   c[1][1]=c[0][0];     c[1][2]=c[0][1];    c[1][3]=0.0;
  c[2][0]=c[0][1];   c[2][1]=c[0][1];     c[2][2]=c[0][0];    c[2][3]=0.0;

  c[3][0]=c[0][3];   c[3][1]=c[1][3];     c[3][2]=c[2][3];    c[3][3]=g;
}



/**
   Function creates compliance matrix for 3D problems.
   
   @param c  - compliance matrix of the material
   @param ea - actual value of Young's modulus
   
   7.2008 TKo
*/
void creep_effym::matcompl_spacestr (matrix &c, double ea)
{
  double g;

  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/ea;
  
  c[0][0]=1.0/ea;     c[0][1]=-1.0*nu/ea;  c[0][2]=c[0][1];
  c[1][0]=c[0][1];    c[1][1]=c[0][0];     c[1][2]=c[0][1];
  c[2][0]=c[0][1];    c[2][1]=c[0][1];     c[2][2]=c[0][0];

  c[3][3]=g;          c[4][4]=g;           c[5][5]=g;
}



/**
  Function computes value of effective Young modulus from given time.
 
  @retval actual value of Young modulus

  TKr + TKo 7.2008
*/
double creep_effym::actual_modulus()
{
  double e, j, tmp;
  double c0,ac,ag,m,z,r,qf,q;
  double t;

  switch (evf)
  {
  case b3law://b3 model for concrete - basic creep       
    //pro vypocet zacinajici od casu 0.0
    t = tb_time/86400.0 + Mp->time/86400;

    ac=sc+gc;   ag=ac/gc;   m=0.5; //m=0.28+1.0/fc/fc;
    n = 0.1;
    
    if (type_e == 1)
      q1=600000.0/e28;//measured
    else
      q1=600000.0/57000.0/sqrt(fc);//empirical
    
    r=1.7*pow(tl,0.12)+8.0;
    qf=1.0/(0.086*pow(tl,(2.0/9.0))+1.21*pow(tl,(4.0/9.0)));
    z = 1.0/pow(tl,m)*log(1.0 + pow((t - tl),n));
    q = qf/pow((1.0 + pow((qf/z),r)),(1.0/r));
    q2=451.1*sqrt(cs)/pow(fc,0.9);
    q3=0.29*pow(wc,4.0)*q2;
    q4=0.14/pow(ac,0.7);
    
    c0=q2*q+q3*log(1.0+pow((t-tl),n))+q4*log(t/tl);

    j = (q1 + c0)/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1

    e = 1.0/j; //effective modulus
    if (e0 == 0.0)
      e0 = e;
    break;
  case doublepwrlaw://toto neni double power law - je nutne opravit
    //pro vypocet zacinajici od casu 0.0
    t = tb_time/86400.0 + Mp->time/86400;
    m=0.5;
    e=e0;
    tmp = 1+psi*(pow(t0,-m)+alpha)*pow(t,n);
    j = 1/e0 + qs*log(tmp);
    e = 1.0/j;
    break;
  case userfunc:
    gf.getval(Mp->time);
    break;
  default:
    print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
    abort();
  }
  return e;
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array

*/
void creep_effym::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param epscr - vector of irreversible strains

  Created by TKo 2008
*/
void creep_effym::giveirrstrains (long ipp, long im, long ido, vector &epscr)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    epscr[i]=Mm->ip[ipp].eqother[ido+1+i];
  }
}



/**
  This function returns the actual value of Young's modulus.

  Created by TKo 6.2016
*/
double creep_effym::give_actual_ym ()
{
  double e =  actual_modulus();

  return e;
}



/**
  This function returns the actual value of Poisson's ratio.

  Created by TKo 6.2016
*/
double creep_effym::give_actual_nu ()
{
  return nu;
}



/**
  This function returns the initial value of Young's modulus.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  Created by TKo 6.2016
*/
double creep_effym::give_initial_ym ()
{
  return e0;
}



/**
  This function returns the value of tensile strength

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  
  Created by TKo 2008
*/
double creep_effym::give_actual_ft (long ipp, long im, long ido)
{
  double ft;
  if (im > 0)
    ft = Mm->give_actual_ft(ipp, im+1, ido);
  else
  {
    print_err("creep_effym model does not contain tensile strength", __FILE__, __LINE__, __func__);
    abort();
  }
  return ft;
}
