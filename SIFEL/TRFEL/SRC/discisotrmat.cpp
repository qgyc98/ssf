/*
  File:             discisotrmat.cpp
  Author:           Jaroslav Kruis, 15.8.2006
  Purpose:          Calculates properties of general isotropic material for linear onemedium transfer
                    with discontinuity, it serves for debugging
*/ 

#include "discisotrmat.h"
#include "stochdrivert.h"

discisotrmat::discisotrmat (void)
{
  k=0.0;  c=0.0;  ju=0.0;
}
discisotrmat::~discisotrmat (void)
{
}



/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void discisotrmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    fprintf (stderr,"\n in function (file %s, line %d).\n",__FILE__,__LINE__);
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
void discisotrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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
void discisotrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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

void discisotrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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
void discisotrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();
}


/**
   function reads parameters
   
   @param in - input file
*/
void discisotrmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&c,&k,&ju);
}


/**
   function prints parameters
   
   @param out - outut file
*/
void discisotrmat::print (FILE *out)
{
  fprintf (out,"  %e %e %e",c,k,ju);
}

/**
   function creates conductivity coefficient of the isotropic material

   @retval k - heat conductivity coefficient of the isotropic material
*/

double discisotrmat::get_k()
{
  return(k);
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double discisotrmat::get_c()
{
  return(c);
}


/**
   function changes parameters of conductivity and capacity from a table
   @ param
*/
void discisotrmat::changeparam (atsel &atm,vector &val)
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
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}

/**
   function prescribes jump in values
*/
double discisotrmat::correct_val (double */*in*/,double */*iin*/)
{
  //return in[0]+ju;
  return ju;
}

/**
   function computes relative 
*/
double discisotrmat::compute_rel (double in)
{
  double out;
  //out=in/ju;
  out=in+ju;
  return out;
}

/**
   function computes absolute
*/
double discisotrmat::compute_abs (double in)
{
  double out;
  //out=ju*in;
  out=ju+in;
  return out;
}
