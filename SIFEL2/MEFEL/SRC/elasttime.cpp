#include "elasttime.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include <math.h>


/**
  Constructor initializes internal data prameters

  7.2008 TKo, TKr
*/
elasttime::elasttime (void)
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
  udstiff = NULL;

  //double power law - je nutne opravit ??!!
  qs = 0.0;  psi = 0.0; n = 0.0; m = 0.0; alpha = 0.0;  
  t0 = 0.0;
}



elasttime::~elasttime (void)
{
  if (udstiff)
  {
    delete udstiff[0];
    delete udstiff[1];
  }
  delete [] udstiff;
}



/**
   Function reads material parameters.
   
   @param in - input stream

   7.2008 TKo, TKr
*/
void elasttime::read (XFILE *in)
{
  long readgft, readgfs;

  //xfscanf (in,"%lf",&nu);
  //xfscanf (in, "%ld",(int*)&evf);
  xfscanf (in,"%k%lf","nu",&nu);
  xfscanf (in, "%k%m", "ym_evolfunc", &ym_evolfunc_kwdset, &evf);
  switch (evf)
    {
    case b3law:
      xfscanf(in,"%ld",&type_e);
      
      if (type_e == 1){xfscanf (in,"%lf ",&e28);}
      //input in Pa
      e0 = e28;
      e28 = e28/6.89476*1.0e-3;//to psi
      
      xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &fc, &wc, &sc, &gc, &cs, &a1, &a2, &ks, &kd, &tb_time);
      
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
      xfscanf(in, "%k", "user_def_stiff_function_table");
      itab.read(in);  // indeces in itab.y have already been decremented in the read function
      readgft = 0;  
      readgfs = 0;
      for(long i=0; i<itab.asize; i++)
      {
        if (itab.y[i] == 1)  // time dependent stiffness
        {
          readgft = 1;
          continue;
        }
        if (itab.y[i] == 2)  // strain/displacement dependent stiffness
        {
          readgfs = 1;
          continue;
        }
        print_err("invalid index %d is given for time %le, id must be in 1 or 2", __FILE__, __LINE__, __func__, itab.y[i], itab.x[i]);
        abort();
      }
      udstiff = new gfunct*[2];
      memset(udstiff, 0, sizeof(*udstiff)*2);
      if (readgft)
      {
        xfscanf(in, "%k", "time_dep_stiff");
        udstiff[0] = new gfunct;
        udstiff[0]->read(in);
      }
      if (readgfs)
      {
        xfscanf(in, "%k", "strain_dep_stiff");
        udstiff[1] = new gfunct;
        udstiff[1]->read(in);
      }
      break;
    default:
      print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
      abort();
  }
}



/**
   Function reads material parameters.
   
   @param in - input stream

   7.2008 TKo, TKr
*/
void elasttime::print (FILE *out)
{
  fprintf (out, "%le ", nu);
  fprintf (out, "%d ", int(evf));
  switch (evf)
  {
    case b3law:
      fprintf(out, "%ld ", type_e);
      
      if (type_e == 1)
        fprintf (out, "%le ", e28);
      fprintf (out, "%le %le %le %le\n  %le %le %le %le %le %le", fc, wc, sc, gc, cs, a1, a2, ks, kd, tb_time);
      
      break;
    case doublepwrlaw:
      fprintf(out, "%le %le %le %le %le\n  %le %le", e0, t0, qs, psi, m, alpha, n);
      break;
    case userfunc:
      fprintf(out, "%le\n", e0);
      itab.print(out);
      if (udstiff[0])
        udstiff[0]->print(out);
      if (udstiff[1])
        udstiff[1]->print(out);
      break;
    default:
      print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
  }
}



/**
  Function creates stiffness matrix depending on stress-strain state and 
  type of required stiffness matrix in problem description Mp->nlman->stmat 
  (initial, tangent,...).

  @param d    - stiffness matrix
  @param ipp - number of integration point

  7.2008 TKo
*/
void elasttime::matstiff (matrix &d, long ipp)
{
  matstiff (d, ipp, Mp->time);
}



/**
  Function creates stiffness matrix depending on stress-strain state and 
  type of required stiffness matrix in problem description Mp->nlman->stmat 
  (initial, tangent,...).

  @param d    - stiffness matrix
  @param ipp - number of integration point
  @param time - time for the stiffness evaluation

  7.2008 TKo
*/
void elasttime::matstiff (matrix &d, long ipp, double time)
{
  strastrestate ssst = Mm->ip[ipp].ssst;
  switch (ssst)
  {
    case bar:
      matstiff_bar (d, ipp, time, secant_stiff);
      break;
    case plbeam:
      matstiff_plbeam (d, ipp, time, secant_stiff);
      break;
    case spacebeam:
      matstiff_spacebeam (d, ipp, time, secant_stiff);
      break;
    case planestress:
      matstiff_plstress (d, ipp, time, secant_stiff);
     break;
    case planestrain:
      matstiff_plstrain (d, ipp, time, secant_stiff);
      break;
    case platek:
      matstiff_plate (d, ipp, time, secant_stiff);
      break;
    case plates:
      matstiff_plate (d, ipp, time, secant_stiff);
      break;
    case axisymm:
      matstiff_axi (d, ipp, time, secant_stiff);
      break;
    case spacestress:
      matstiff_spacestr (d, ipp, time, secant_stiff);
      break;
    default:
      print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function creates intial stiffness matrix depending on stress-strain state.

  @param d    - stiffness matrix
  @param ipp - integration point pointer
  @param time - time for the stiffness evaluation

  7.2008 TKo
*/
void elasttime::elmatstiff (matrix &d, long ipp)
{
  elmatstiff(d, ipp, Mp->time);
}



/**
  Function creates intial stiffness matrix depending on stress-strain state.

  @param d    - stiffness matrix
  @param ipp - integration point pointer
  @param time - time for the stiffness evaluation

  7.2008 TKo
*/
void elasttime::elmatstiff (matrix &d, long ipp, double time)
{
  strastrestate ssst = Mm->ip[ipp].ssst;
  switch (ssst){
  case bar:{
    matstiff_bar (d, ipp, time, initial_stiff);
    break;
  }
  case plbeam:{
    matstiff_plbeam (d, ipp, time, initial_stiff);
    break;
  }
  case spacebeam:{
    matstiff_spacebeam (d, ipp, time, initial_stiff);
    break;
  }
  case planestress:{
    matstiff_plstress (d, ipp, time, initial_stiff);
    break;
  }
  case planestrain:{
    matstiff_plstrain (d, ipp, time, initial_stiff);
    break;
  }
  case platek:{
    matstiff_plate (d, ipp, time, initial_stiff);
    break;
  }
  case plates:{
    matstiff_plate (d, ipp, time, initial_stiff);
    break;
  }
  case axisymm:{
    matstiff_axi (d, ipp, time, initial_stiff);
    break;
  }
  case spacestress:{
    matstiff_spacestr (d, ipp, time, initial_stiff);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   Function creates stiffness matrix for bar elements.
   
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_bar (matrix &d, long ipp, double time, stiffmatrix smt)
{
  /*
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus();
  */
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  d[0][0] = e;
}



/**
   Function creates stiffness matrix for plane beam elements.
   
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_plbeam (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  d[0][0] = e;
  d[1][1] = e/2.0/(1.0+nu);
  d[2][2] = e;
}



/**
   Function creates stiffness matrix for space beam elements.
   
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_spacebeam (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  d[0][0] = e;
  d[1][1] = e/2.0/(1.0+nu);
  d[2][2] = e/2.0/(1.0+nu);
  d[3][3] = e/2.0/(1.0+nu);
  d[4][4] = e;
  d[5][5] = e;
}



/**
   Function creates stiffness matrix for 2D problems (plane stress).

   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_plstress (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double c;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
  fillm(0.0,d);

  c = e/(1.0-nu*nu);
  
  d[0][0] = c;     d[0][1] = c*nu;  d[0][2] = 0.0;
  d[1][0] = c*nu;  d[1][1] = c;     d[1][2] = 0.0;
  d[2][0] = 0.0;   d[2][1] = 0.0;   d[2][2] = e/2.0/(1.0+nu);
}



/**
   Function creates stiffness matrix for 2D problems (plane strain).
 
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_plstrain (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double c;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
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
   Function creates stiffness matrix for 2D problems (axisymmetric problem).

   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_axi (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double g,s;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
  fillm(0.0,d);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;

}



/**
   Function creates stiffness matrix for plate elements.
   
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_plate (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double c,g;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
  fillm(0.0,d);

  c = e/12.0/(1.0-nu*nu);
  g = e/2.0/(1.0+nu);
  
  d[0][0]=c;        d[0][1]=c*nu;  d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c;     d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;   d[2][2]=g/12.0;
  
  d[3][3]=g;  d[4][4]=g;
}



/**
   Function creates stiffness matrix for 3D problems.
   
   @param d   - stiffness matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   TKo 7.2008
*/
void elasttime::matstiff_spacestr (matrix &d, long ipp, double time, stiffmatrix smt)
{
  double g,s;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
  fillm(0.0,d);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;
}



/**
  Function creates complience matrix depending on stress-strain state and 
  type of required stiffness matrix in problem description Mp->nlman->stmat (initial, tangent,...).

  @param c    - compliance matrix of the material
  @param time - time for the stiffness evaluation
  @param ipp - number of integration point

  7.2008 TKo
*/
void elasttime::matcompl (matrix &c, long ipp)
{
  matcompl(c, ipp, Mp->time);
}



/**
  Function creates complience matrix depending on stress-strain state and 
  type of required stiffness matrix in problem description Mp->nlman->stmat (initial, tangent,...).

  @param c    - compliance matrix of the material
  @param time - time for the stiffness evaluation
  @param ipp - number of integration point

  7.2008 TKo
*/
void elasttime::matcompl (matrix &c, long ipp, double time)
{
  strastrestate ssst = Mm->ip[ipp].ssst;
  switch (ssst)
  {
    case bar:
      matcompl_bar (c, ipp, time, Mp->nlman->stmat);
      break;
    case plbeam:
      matcompl_plbeam (c, ipp, time, Mp->nlman->stmat);
      break;
    case planestress:
      matcompl_plstress (c, ipp, time, Mp->nlman->stmat);
      break;
    case planestrain:
      matcompl_plstrain (c, ipp, time, Mp->nlman->stmat);
      break;
    case axisymm:
      matcompl_axi (c, ipp, time, Mp->nlman->stmat);
      break;
    case spacestress:
      matcompl_spacestr (c, ipp, time, Mp->nlman->stmat);
      break;
    default:
      print_err ("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
}



/**
   Function creates compliance matrix for bar elements.
   
   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_bar (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  c[0][0] = 1.0/e;
}



/**
   Function creates compliance matrix for plane beam elements.
   
   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_plbeam (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  c[0][0] = 1.0/e;
  c[1][1] = 2.0*(1.0+nu)/e;
  c[2][2] = 1.0/e;
}



/**
   Function creates compliance matrix for 2D problems (plane stress).

   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_plstress (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  fillm(0.0,c);
  
  c[0][0] =  1.0/e;     c[0][1] = -1.0*nu/e;  c[0][2] = 0.0;
  c[1][0] = -1.0*nu/e;  c[1][1] =  1.0/e;     c[1][2] = 0.0;
  c[2][0] = 0.0;        c[2][1] = 0.0;        c[2][2] = 2.0*(1.0+nu)/e;
}



/**
   Function creates compliance matrix for 2D problems (plane strain).

   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_plstrain (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double g;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  
  fillm(0.0,c);
  
  g = (1.0+nu)/e;
  
  c[0][0] = g*(1.0-nu);   c[0][1] = -1.0*g*nu;    c[0][2] = 0.0;
  c[1][0] = -1.0*g*nu;    c[1][1] = g*(1.0-nu);   c[1][2] = 0.0;
  c[2][0] = 0.0;          c[2][1] = 0.0;          c[2][2] = 2.0*g;
}



/**
   Function creates compliance matrix for axisymmetric problems
   
   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_axi (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double g;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/e;
  
  c[0][0]=1.0/e;     c[0][1]=-1.0*nu/e;  c[0][2]=c[0][1];    c[0][3]=0.0;
  c[1][0]=c[0][1];   c[1][1]=c[0][0];    c[1][2]=c[0][1];    c[1][3]=0.0;
  c[2][0]=c[0][1];   c[2][1]=c[0][1];    c[2][2]=c[0][0];    c[2][3]=0.0;

  c[3][0]=c[0][3];   c[3][1]=c[1][3];    c[3][2]=c[2][3];    c[3][3]=g;
}



/**
   Function creates compliance matrix for 3D problems.
   
   @param c   - compliance matrix of the material
   @param ipp - number of integration point
   @param time - time for the stiffness evaluation
   @param smt - type of assembled stiffness matrix (initial, tangent, ...)
   
   7.2008 TKo
*/
void elasttime::matcompl_spacestr (matrix &c, long ipp, double time, stiffmatrix smt)
{
  double g;
  double e;
  if (smt == initial_stiff)
    e = e0;
  else
    e = actual_modulus(ipp, time);

  fillm(0.0,c);
  
  g = 2.0*(1.0+nu)/e;
  
  c[0][0]=1.0/e;     c[0][1]=-1.0*nu/e;  c[0][2]=c[0][1];
  c[1][0]=c[0][1];   c[1][1]=c[0][0];    c[1][2]=c[0][1];
  c[2][0]=c[0][1];   c[2][1]=c[0][1];    c[2][2]=c[0][0];

  c[3][3]=g;         c[4][4]=g;          c[5][5]=g;
}



/**
  Function computes value of effective Young modulus from given time.
 
  @param ipp - number of integration point

  @retval actual value of Young modulus
  TKo 7.2008
*/
double elasttime::actual_modulus(long ipp)
{
  return actual_modulus(ipp, Mp->time);
}



/**
  Function computes value of effective Young modulus from given time.
   
  @param ipp - number of integration point
  @param time - time for the stiffness evaluation

  @retval actual value of Young modulus
  TKo 7.2008
*/
double elasttime::actual_modulus(long ipp, double time)
{
  double e, j, tmp;
  double c0,ac,ag,m,z,r,qf,q;
  double t;
  long id;

  switch (evf){
    case b3law://b3 model for concrete - basic creep       
      //pro vypocet zacinajici od casu 0.0
      t = tb_time/86400.0 + time/86400;
      
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
      t = tb_time/86400.0 + time/86400;
      m=0.5;
      e=e0;
      tmp = 1+psi*(pow(t0,-m)+alpha)*pow(t,n);
      j = 1/e0 + qs*log(tmp);
      e = 1.0/j;
      break;
    case userfunc:
      id = itab.getval(time)-1;
      switch (id){
        case 0:
          e = udstiff[id]->getval(time);
          break;
        case 1:
          if (Mm->ip[ipp].ncompstr == 1)
            e = udstiff[id]->getderiv(Mm->ip[ipp].strain[0]);
          else{
            print_err("displacement control of stiffness must be used only for 1D stress/strain state (time %le),\n"
                      " i.e. ncomptr must be 1 but it is %ld", __FILE__, __LINE__, __func__, time, Mm->ip[ipp].ncompstr);
            abort();
          }
          break;
        default:
          print_err("invalid index %ld required for user defined stiffness", __FILE__, __LINE__, __func__, id);
          abort();
      }
      break;
    default:
      print_err("unknown type of function of Young modulus evolution", __FILE__, __LINE__, __func__);
      abort();
  }
  return e;
}



/**
  Function returns value of initial Young modulus.
 
  @retval actual value of Young modulus

  TKo 7.2008
*/
double elasttime::initial_modulus()
{
  return e0;
}



/**
  The function computes true stresses
   
  @param ipp - number of integration point
  @param ido - index of internal variables for given material in the ipp other array   

  Created by Tomas Koudelka, 05.2018
*/
void elasttime::nlstresses (long ipp, long /*ido*/)
{
  long i, j, n = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(n)), sig(ASTCKVEC(n));
  matrix d(ASTCKMAT(n,n));
  double t = Mp->time;
  long id;

  //  initial values
  Mm->givestrain(0, ipp, eps);

  if (evf == userfunc)
  {
    id = itab.getval(t)-1;
    if (id == 1)   // strain dependent stresses
    {
      for (j=0; j<n; j++)
        Mm->ip[ipp].stress[j] = udstiff[id]->getval(eps[j]);
      return;
    }
    else  // time dependent stiffness
    {
      matstiff(d, ipp, t);
      mxv(d, eps, sig);
      for (i=0;i<n;i++)
        Mm->ip[ipp].stress[i]=sig[i];
    }
  }
  else
  {
    print_err("stress computation for evf=%d has not been implemented yet", __FILE__, __LINE__, __func__, int(evf));
    abort();
  }
      
  
  
/*
  long i, k,  n=Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(n)), epso(ASTCKVEC(n)), deps(ASTCKVEC(n)), sigo(ASTCKVEC(n)), sig(ASTCKVEC(n));
  vector epsa(ASTCKVEC(n)), depsa(ASTCKVEC(n)), dsig(ASTCKVEC(n)), sig_h(ASTCKVEC(n)), sig_h2(ASTCKVEC(n));
  vector sig_err(ASTCKVEC(n));
  matrix d(ASTCKMAT(n,n));
  double h, dh, t, dt, to, norm_r;
  long ni = 10000, stopk;
  
  //  initial values
  Mm->givestrain(0, ipp, eps);
  for(i=0; i<n; i++)
  {
    epso[i] = Mm->ip[ipp].eqother[ido+i];
    sigo[i] = Mm->ip[ipp].eqother[ido+n+i];
  }
  to = Mp->time - Mp->timecon.actualbacktimeincr();
  dt = Mp->timecon.actualbacktimeincr();

  subv(eps, epso, deps);
  copyv(sigo, sig);
  
  dh = 1.0; // actual relative step length
  h  = 0.0; // accumulated relative step length
  stopk=0;

  for (k=0; k<ni && h<1.0; k++)
  {
    if (dh==0.0)
    {
      print_err("dh=0.0 - invalid state of integration", __FILE__, __LINE__, __func__);
      abort();
    }

    copyv(sig, sig_h);
    copyv(sig, sig_h2);
    cmulv(dh, deps, depsa);

    t = to+h*dt;
    matstiff (d, ipp, t);
    mxv (d,depsa,dsig);
    addv(sig_h, dsig, sig_h);
    addmultv(sig_h2, dsig, 0.5, sig_h2);

    t += 0.5*dh*dt;
    cmulv(0.5*dh, deps, depsa);
    matstiff (d, ipp, t);
    mxv (d,depsa,dsig);
    addv(sig_h2, dsig, sig_h2);

    subv(sig_h2, sig_h, sig_err);
    norm_r = normv(sig_err);
    if (norm_r < 1.0e-8)
    {
      h += dh;
      dh *= 2.0;
      if (dh > 1.0)
        dh = 1.0;
      // check the maximum RK substep length
      if (1.0-h < dh)
      {
        dh = 1.0-h;
        stopk = 1;
      }
      copyv(sig_h2, sig);
    }
    else
    {
      t -= dh*dt;
      if (dh == 1.0e-7)
      {
        print_err("cannot reduce integration step length in %ld step", __FILE__, __LINE__, __func__, i);
        abort();
      }
      dh *= 0.5;
      if (dh < 1.0e-7) 
        dh = 1.0e-7;
    }
  }

  if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (eps[0]+eps[1]);

  for (i=0;i<n;i++)
  {
    Mm->ip[ipp].stress[i]=sig[i];
    Mm->ip[ipp].other[i] = eps[i];
    Mm->ip[ipp].other[i+n] = sig[i];
  }
*/
}



/**
  The function updates values in the eqother array containing values from the previous equlibrium state 
  to the values attained in the actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void elasttime::updateval (long ipp, long im, long ido)
{
  long i, n = Mm->givencompeqother(ipp,im);

  for (i=0; i<n; i++)
    Mm->ip[ipp].eqother[ido+i] = Mm->ip[ipp].other[ido+i];
}



void elasttime::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++)
  {
    switch (atm.atrib[i])
    {
      case 0:
        e0=val[i];
        break;
      case 1:
        nu=val[i];
        break;
      case 2:
        qs = val[i];  
        break;
      case 3:
       psi = val[i]; 
        break;
      case 4:
        n = val[i]; 
        break;
      case 5:
        m = val[i];  
        break;
      case 6:
        alpha = val[i];  
        break;
      default:
        print_err("wrong number of atribute in function changeparam",__FILE__,__LINE__, __func__);
    }
  }
}
