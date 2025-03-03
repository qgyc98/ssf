/*
    Author:           Jakub Goringer, email@
    Purpose:          Computes carbonation of concrete (uncoupled from moisture)
    Source:           CCR
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "carbmat1.h"
#include "globalt.h"

carbmat1::carbmat1(){
  global_variable = 3;
}

carbmat1::~carbmat1(){

}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void carbmat1::matcond (matrix &d,long ri,long ci,long ipp)
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
void carbmat1::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_conductivity();
  
  fillm(0.0,d);

  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void carbmat1::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_conductivity();

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

void carbmat1::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_conductivity();
  
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
void carbmat1::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
}


/**
   function reads parameters from input file
   @param in - input file
*/
void carbmat1::read(XFILE *in){
  xfscanf (in,"%lf", &param);
}


/**
   function computes effective CO2 diffusion of cracked material
   @param A - area of cone
*/

double carbmat1::Equivalent_diffusion(double A)
{
  double K_CO2 = 1.e-12;
  double K_Cracked;
  
  K_Cracked = A * K_CO2 * global_variable;
  
  return K_Cracked;
}


double carbmat1::get_conductivity(void){
  return param;
}
