/*
    File:             kunzel.cpp
    Author:           Jiri Madera
    Purpose:          c
    Source:           Kunzel, 19xx;
    Assumptions:


*/


#include "richardscontam.h"
#include "richardscontam.h"
#include "globalt.h"


richardscontam::richardscontam()
{
  //  parameters of van Genuchten model
  alpha=0.0;
  n=0.0;
  m=0.0;
  //  dimension of problem solved
  dim=0;
  
  //  saturated hydraulic conductivity
  kksxx=0.0;
  kksxy=0.0;
  kksxz=0.0;
  kksyy=0.0;
  kksyz=0.0;
  kkszz=0.0;
  //  saturated water content
  thetas=0.0;
  //  residual water content
  thetar=0.0;
  
}

richardscontam::~richardscontam()
{

}

/**
   function reads material parameters
   
   @param in - input file
   
   1.12.2011
*/
void richardscontam::read(XFILE *in)
{
  //  dimenison of problem solved
  //  the dimension is determined from the first element
  dim = Tt->give_dimension (0);
  if (dim<1 || dim>3){
    print_err("wrong dimension of problem solved (%ld)",__FILE__,__LINE__,__func__,dim);
  }

  //  parameters of van Genuchten model
  //  saturated hydraulic conductivity
  //  saturated water content
  //  residual water content
  //  specific storage
  rich.read (in);
  
  if (dim==1){
    xfscanf (in,"%lf",&kksxx);
  }
  if (dim==2){
    xfscanf (in,"%lf %lf %lf",&kksxx,&kksxy,&kksyy);
  }
  if (dim==3){
    xfscanf (in,"%lf %lf %lf %lf %lf %lf",&kksxx,&kksxy,&kksxz,&kksyy,&kksyz,&kkszz);
  }
  
  //  cteni parametru, ktere nejsou v Richardsovi
  xfscanf (in,"%lf %lf %lf %lf %lf %lf",&alpha,&n,&m,&thetas,&thetar,&storage);
}

/**
   function prints material parameters
   
   @param out - output file
   
   1.12.2011
*/
void richardscontam::print(FILE *out)
{
  rich.print (out);
  
  if (dim==1){
    fprintf (out," %lf\n",kksxx);
  }
  if (dim==2){
    fprintf (out," %lf %lf %lf\n",kksxx,kksxy,kksyy);
  }
  if (dim==3){
    fprintf (out," %lf %lf %lf %lf %lf %lf\n",kksxx,kksxy,kksxz,kksyy,kksyz,kkszz);
  }
  
  //  parameters of van Genuchten model
  //  saturated hydraulic conductivity
  //  saturated water content
  //  residual water content
  fprintf (out," %lf %lf %lf %lf %lf %lf\n",alpha,n,m,thetas,thetar,storage);
}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void richardscontam::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of flux components is required",__FILE__,__LINE__,__func__);
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
void richardscontam::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double h,k;
  
  if (ri==0 && ci==0){
    //  K_hh
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];

    k = rich.kk_value (h);
  }
  if (ri==0 && ci==1){
    //  K_hc
    k=0.0;
  }
  if (ri==1 && ci ==0){
    //  K_ch
    k=0.0;
  }
  if (ri==1 && ci==1){
    //  K_cc
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];
    
    //  hydraulic conductivity
    k = cc_value (h);
    
  }

  
  d[0][0] = k*kksxx;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richardscontam::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double h,k;
  
  if (ri==0 && ci==0){
    //  K_hh
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];

    k = rich.kk_value (h);
  }
  if (ri==0 && ci==1){
    //  K_hc
    k=0.0;
  }
  if (ri==1 && ci ==0){
    //  K_ch
    k=0.0;
  }
  if (ri==1 && ci==1){
    //  K_cc
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];
    
    //  hydraulic conductivity
    k = cc_value (h);
    
  }
  
  d[0][0] = k*kksxx;   d[0][1] = k*kksxy;
  d[1][0] = k*kksxy;   d[1][1] = k*kksyy;

}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void richardscontam::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double h,k;
  
  if (ri==0 && ci==0){
    //  K_hh
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];

    k = rich.kk_value (h);
  }
  if (ri==0 && ci==1){
    //  K_hc
    k=0.0;
  }
  if (ri==1 && ci ==0){
    //  K_ch
    k=0.0;
  }
  if (ri==1 && ci==1){
    //  K_cc
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];
    
    //  hydraulic conductivity
    k = cc_value (h);
    
  }
  
  d[0][0] = k*kksxx;   d[0][1] = k*kksxy;  d[0][2]=k*kksxz;
  d[1][0] = k*kksxy;   d[1][1] = k*kksyy;  d[1][2]=k*kksyz;
  d[2][0] = k*kksxz;   d[2][1] = k*kksyz;  d[2][2]=k*kkszz;
}

void richardscontam::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;
  
  switch (n){
  case 1:{
    matcond1d_2 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d_2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d_2 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function creates conductivity %matrix of the material for 1D problems (convective term)
   
   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richardscontam::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double h;
  
  fillm(0.0,d);
  if (ri==0 && ci==0){
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];
    
    d[0][0]=0.0-rich.dkkdh_value (h);
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems (convective term)

   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void richardscontam::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double h;
  
  fillm(0.0,d);
  if (ri==0 && ci==0){
    //  actual hydraulic head
    h = Tm->ip[ipp].av[0];
    
    d[0][0]=0.0;
    d[0][1]=0.0-rich.dkkdh_value (h);
  }
}

/**
   function creates conductivity %matrix of the material for 1D problems (convective term)

   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richardscontam::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double h;
  
  fillm(0.0,d);
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  if (ri==0 && ci==0){
    d[0][0]=0.0;
    d[0][1]=0.0;
    d[0][2]=0.0-rich.dkkdh_value (h);
  }
}



/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richardscontam::matcap (double &/*c*/,long /*ri*/,long /*ci*/,long ipp)
{
  double h;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  //  capacity
  //c = c_value (h);

}

double richardscontam::cc_value (double /*h*/)
{
  return 0.0;
}
