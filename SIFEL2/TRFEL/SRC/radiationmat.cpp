/*
  File:             radiationmat.cpp
  Author:           Jaroslav Kruis, 20.12.2002, modified 13.06.2011 by Jan Koci & Jiri Madera
  Purpose:          Calculates properties of general isotropic material for linear onemedium transfer with radiative exchange of energy 
*/ 

#include "radiationmat.h"
#include "stochdrivert.h"

radiationmat::radiationmat (void)
{
  k=0.0;//conductivity
  c=0.0;//capacity
}
radiationmat::~radiationmat (void)
{}




/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void radiationmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void radiationmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void radiationmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
  fillm(0.0,d);
  
  d[0][0] = kk;   d[0][1] = 0.0;
  d[1][0] = 0.0;  d[1][1] = kk;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void radiationmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
  fillm(0.0,d);
  
  d[0][0]=kk;   d[0][1]=0.0;  d[0][2]=0.0;
  d[1][0]=0.0;  d[1][1]=kk;   d[1][2]=0.0;
  d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=kk;
}


/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void radiationmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();
}


/**
   function reads parameters
   
   @param in - input file
*/
void radiationmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf",&c,&k); //c=capacity, k=conductivity
}


/**
   function prints parameters
   
   @param out - outut file
*/
void radiationmat::print (FILE *out)
{
  fprintf (out,"  %e %e",c,k);
}

/**
   function creates conductivity coefficient of the isotropic material

   @retval k - heat conductivity %matrix of the isotropic material
*/

double radiationmat::get_k()
{
  return(k);
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double radiationmat::get_c()
{
  return(c);
}


/**
   function changes parameters of conductivity and capacity from a table
   @ param
*/
void radiationmat::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      c=val[0];
      break;
    }
    case 1:{
      k=val[1];
      break;
    }
    default:{
      print_err("wrong number of atribute is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}
