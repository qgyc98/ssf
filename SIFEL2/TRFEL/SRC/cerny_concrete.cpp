/*
    File:             cerny_concrete.cpp
    Author:           Tomas Krejci, 1.12.2002
    Purpose:          computes conductivity and capacity coefficients for fiber concrete 
                      for heat transfer (one medium, t ... temperature [deg. C])
                      - data measured in the laboratory of the Department of Physics 
                      of the Faculty of Civil Engineering CTU Prague
    Source:           
*/ 

#include "cerny_concrete.h"
#include "globalt.h"

cernymat::cernymat (void)
{
  rho=0.0; k=0.0;  c=0.0;
}
cernymat::~cernymat (void)
{}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void cernymat::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void cernymat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double t;
  
  t = Tm->ip[ipp].av[0];

  kk = get_k(t);
  
  d[0][0] = kk;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void cernymat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double t;
  
  t = Tm->ip[ipp].av[0];

  kk = get_k(t);
  
  fillm(0.0,d);
  
  d[0][0] = kk;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = kk;
}

/**
   function creates conductivity matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void cernymat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double t;
  
  t = Tm->ip[ipp].av[0];

  kk = get_k(t);
  
  fillm(0.0,d);
  
  d[0][0]=kk;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=kk;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=kk;
}


/**
   function creates capacity matrix of the material

   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void cernymat::matcap (double &cc,long /*ri*/,long /*ci*/,long ipp)
{
  double t;
  cc = 0.0;
    
  t = Tm->ip[ipp].av[0];
  
  cc = get_c(t);
}

/**
   function reads parameters
   
   @param in - input file
*/
void cernymat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&rho,&c,&k);
}

/**
   function prints parameters
   
   @param out - outut file
*/
void cernymat::print (FILE *out)
{
  fprintf (out,"  %e %e %e",rho,c,k);
}

/**
   function creates effective heat conductivity of the concrete
   - data measured in the laboratory of the Department of Physics 
   of the Faculty of Civil Engineering CTU Prague   

   @param t - temperature in deg. C

   @retval keff - heat conductivity of concrete
*/
double cernymat::get_k(double t)
{
  double keff;

  //keff = k - 0.00088*t*t;//k - 0.0028*t;
  keff = k - 0.0028*t;

  return(keff);
}

/**
   function creates derivative of effective heat conductivity of the concrete with respect to temperature
   - data measured in the laboratory of the Department of Physics 
   of the Faculty of Civil Engineering CTU Prague   

   @param t - temperature in deg. C

   @retval d - derivative of effective heat conductivity with respect to temperature
*/
double cernymat::get_dk_dt(double /*t*/)
{
  double d;

  //d = -0.00088*t/2.0;//-0.0028;
  d = -0.0028;

  return(d);
}

/**
   function creates effective heat capacity of concrete
   - data measured in the laboratory of the Department of Physics 
   of the Faculty of Civil Engineering CTU Prague   

   @param t - temperature in deg. C

   @retval ceff - effective heat capacity of concrete
*/
double cernymat::get_c(double t)
{
  double ceff;
  
  ceff = rho*(c + 0.25*t);
  //ceff = rho*(c + 0.025*t*t);

  return(ceff);
}

/**
   function creates derivative of effective heat capacity with respect to temperature
   - data measured in the laboratory of the Department of Physics 
   of the Faculty of Civil Engineering CTU Prague   

   @param t - temperature in deg. C

   @retval d - derivative of effective heat capacity with respect to temperature
*/
double cernymat::get_dc_dt(double /*t*/)
{
  double d;
  
  d = rho*0.25;
  //d = rho*0.025*t/2.0;

  return(d);
}

