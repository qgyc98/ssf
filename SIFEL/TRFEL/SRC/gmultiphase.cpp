/*
  File:             gmultiphase.cpp
  Author:           Tomas Krejci, 7.3.2006
  Purpose:          computes conductivity and capacity matrices in a material point for multiphase porous media
                    of geomaterials
  
  FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
  --------------------------------------------------
  
  sources: THE FINITE ELEMENT METHOD IN THE STATIC AND DYNAMIC DEFORMATION AND CONSOLIDATION OF POROUS MEDIA
           R. W. LEWIS, B. A. SCHREFLER
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "gmultiphase.h"
#include "constrel.h"
#include "globalt.h"
#include "globmatt.h"

gmultiph::gmultiph()
{
  mw = 18.01528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  scale_pw = Tp->scale[0];
  scale_pg = Tp->scale[1];
  scale_t = Tp->scale[2];
}
gmultiph::~gmultiph()
{}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void gmultiph::matcond (matrix &d,long ri,long ci,long ipp)
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
void gmultiph::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void gmultiph::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
}

/**
   function creates conductivity matrix of the material for 3D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void gmultiph::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}


/**
   function creates capacity matrix of the material
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void gmultiph::matcap (double &c,long ri,long ci,long ipp)
{
  double pw,pg,t;
  c = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = get_capwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    c = get_capwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    c = get_capgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    c = get_capgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    c = get_capgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    c = get_captw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    c = get_captg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    c = get_captt(pw,pg,t,ipp);// *scale_t;//scaling
}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void gmultiph::matcond2 (matrix &d,long ri,long ci,long ipp)
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
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}




/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void gmultiph::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd,g;
  double pw,pg,t;
      
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }

  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling
    g = 0.0;//gravity acceleration//??

    //d[0][0] = a*(c*g - Tm->ip[ipp].grad[0][0]) + b*(dd*g - Tm->ip[ipp].grad[1][0]);
    d[0][0] = 0.0;
  }  
}


/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void gmultiph::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd,*g;
  double pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling

    g = new double [2];
    g[0] = 0.0;//??
    g[1] = 0.0;//??
    
    //d[0][0] = a*(c*g[0] - Tm->ip[ipp].grad[0][0]) + b*(dd*g[0] - Tm->ip[ipp].grad[1][0]);
    //d[0][1] = a*(c*g[1] - Tm->ip[ipp].grad[0][1]) + b*(dd*g[1] - Tm->ip[ipp].grad[1][1]);
    d[0][0] = 0.0;
    d[0][1] = 0.0;

    delete [] g;
  }
}



/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void gmultiph::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd,*g;
  double pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
    
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling

    g = new double [3];
    g[0] = 0.0;//??
    g[1] = 0.0;//??
    g[2] = 0.0;//??
    
    
    //d[0][0] = a*(c*g[0] - Tm->ip[ipp].grad[0][0]) + b*(dd*g[0] - Tm->ip[ipp].grad[1][0]);
    //d[0][1] = a*(c*g[1] - Tm->ip[ipp].grad[0][1]) + b*(dd*g[1] - Tm->ip[ipp].grad[1][1]);
    //d[0][2] = a*(c*g[2] - Tm->ip[ipp].grad[0][2]) + b*(dd*g[2] - Tm->ip[ipp].grad[1][2]);
    d[0][0] = 0.0;
    d[0][1] = 0.0;
    d[0][2] = 0.0;

    delete [] g;
  }
}



/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   ...
   @param ipp - number of integration point
   
   11.2.2003
*/
void gmultiph::rhs_volume (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;
  
  switch (m){
  case 1:{
    rhs1d1 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    rhs2d1 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    rhs3d1 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void gmultiph::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;
  double g;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fc1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g = 0.0;//temp
    
    fillm(0.0,d);
    d[0][0] = f*g;
  }
  if(ri == 1){
    f = get_fg(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g = 0.0;//temp

    fillm(0.0,d);
    d[0][0] = f*g;
  }
  if(ri == 2){
    f = get_ft1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g = 0.0;//temp
    
    fillm(0.0,d);
    d[0][0] = f*g;
  }
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void gmultiph::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;
  double *g;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    g = new double [2];
    f = get_fc1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp
    
    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
  }
  if(ri == 1){
    g = new double [2];
    f = get_fg(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp

    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
  }
  if(ri == 2){
    g = new double [2];
    f = get_ft1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp
    
    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
  }
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void gmultiph::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;
  double *g;

  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling

  if(ri == 0){
    g = new double [3];
    f = get_fc1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp
    g[2] = 0.0;//temp

    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
    d[2][0] = f*g[2];
  }
  if(ri == 1){
    g = new double [2];
    f = get_fg(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp
    g[2] = 0.0;//temp

    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
    d[2][0] = f*g[2];
  }
  if(ri == 2){
    g = new double [3];
    f = get_ft1(pw,pg,t,ipp);
    //g = Gt->g;//complete vector of gravity acceleration
    g[0] = 0.0;//temp
    g[1] = 0.0;//temp
    g[2] = 0.0;//temp

    fillm(0.0,d);
    d[0][0] = f*g[0];
    d[1][0] = f*g[1];
    d[2][0] = f*g[2];
  }
}




/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double gmultiph::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t  = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_transcoeff_tt(pw,pg,t,bc,ipp);// *scale_t;//scaling

  c = c*trc;

  return (c);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double gmultiph::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t  = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    c = get_transmission_nodval_ww(nodval,pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_nodval_tt(nodval,trc2,pw,pg,t,bc,ipp);// *scale_t;//scaling

  return (c);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double gmultiph::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t  = nodalval (nn,2);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_flux_ww(nodval,pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_flux_tt(nodval,trc2,pw,pg,t,bc,ipp);// *scale_t;//scaling

  return (c);
}


/**
   function checks if gass pressure is greater than vapour pressure
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point


   TKr, 18.5.2005
*/
void gmultiph::gaspress_check(double pw,double &pg,double t,long ipp)
{
  double pgw,pc;
  state_eq tt;
  
  pc = pg - pw;

  //general
  //gas pressure check
  pgw = tt.get_pgw(pc,t);
  if (pgw > pg) 
    pg = pgw;
  
  //general
  //fully saturated state
  if(pc < 100.0)
    pg = 101325.0;
  
  //storing into gauss point
  Tm->ip[ipp].av[1] = pg;
}


/**
   function checks if capillary pressure is higher than maximum
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void gmultiph::cappress_check(double &/*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  /*  double pgw,pwmax,rhow;
      state_eq tt;
      
      //only boundary nodes
      //maximum capillary pressure check
      pgw = tt.get_pgw(pw,t);
      rhow = tt.get_rhow(t);
      pwmax = 0.0;
      if (pg < pgw && pg > 0.0)
      pwmax = -1.0*rhow*gasr*t/mw*log(pg/pgw);
      else
      pwmax = 0.0;
      
      if (pg < pgw){
      if (pwmax > 0.0){
      pw = pwmax;
      }
      }
      
      
      if (pw > pwmax)
      pw = pwmax;
      
      //storing into gauss point
      Tm->ip[ipp].av[0] = pw;
  */
}


/**
   function fixed capillary pressure if on one element is temperature higher than critical temperature for water
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void gmultiph::cappress_stop(double &/*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  //storing into gauss point
  //Tm->ip[ipp].av[0] = pw;
}


void gmultiph::values_correction (vector &nv)
{
  //  capillary pressure
  cappress_check(nv[0],nv[1],nv[2],0);

  //  gas pressure
  gaspress_check(nv[0],nv[1],nv[2],0);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kww - conductivity coefficient
*/
double gmultiph::get_kww(double pw,double pg,double t,long ipp)
{
  double rhow,pc,kintr,krw,muw,rhog,dg,dpgw_dpc,kww,mg;
  state_eq tt;

  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  rhog = tt.get_rhog(pc,pg,t);
  dg = tt.get_dg(pc,pg,t,ipp);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  mg = tt.get_mg(pc,pg,t);

  //according to Lewis and Schrefler p.394
  kww = rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc - rhow*kintr*krw/muw;//kontrola znamenek??!!

  //printf("kww = %e\n",kww);

  return(kww);
}

/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capww - capacity coefficient
*/
double gmultiph::get_capww(double pw,double pg,double t,long ipp)
{
  double capww,alpha,phi,ks,sw,sg,kw,dpgw_dpc,rhogw,pc,dsw_dpc,rhow;
  state_eq tt;
  
  pc = pg - pw;
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  alpha = tt.get_alpha(pc,pg,t,ipp);
  phi = tt.get_phi(t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  sw = tt.get_s(pc,t,ipp);
  rhogw = tt.get_rhogw(pc,t);
  sg = 1.0 - sw;
  rhow = tt.get_rhow(t);
  kw = tt.get_ks(pc,pg,t,ipp);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);

  //according to Lewis and Schrefler p.394
  capww = (alpha-phi)/ks*sw*(rhogw*sg+rhow*sw)+rhow*sw*phi/kw;
  capww = capww + sg*phi*mw/t/gasr*dpgw_dpc;
  capww = capww + ((alpha-phi)/ks*(rhogw*sg*pc+rhow*sw*pw-rhow*sw*pc) + phi*(rhow-rhogw))*dsw_dpc;

  //printf("capww = %e\n",capww);

  return(capww);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kwg - conductivity coefficient
*/
double gmultiph::get_kwg(double pw,double pg,double t,long ipp)
{
  double kwg,pc,rhogw,kintr,krg,mug,dg,pgw,dpgw_dpc,rhog,mg;
  state_eq tt;
  
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhogw = tt.get_rhogw(pc,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  dg = tt.get_dg(pc,pg,t,ipp);
  pgw = tt.get_pgw(pc,t);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  rhog = tt.get_rhow(t);
  mg = tt.get_mg(pc,pg,t);

  //according to Lewis and Schrefler p.395
  kwg = -1.0*(-rhogw*kintr*krg/mug - rhog*ma*mw/mg/mg*dg*(1.0/pg*dpgw_dpc - pgw/pg/pg));

  //printf("kwg = %e\n",kwg);

  kwg = 0.0;

  return(kwg);
}

/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capwg - capacity coefficient
*/
double gmultiph::get_capwg(double pw,double pg,double t,long ipp)
{
  double capwg,pc,sw,phi,alpha,ks,sg,rhogw,rhow,dpgw_dpc,dsw_dpc;
  state_eq tt;
 
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  sw = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);
  alpha = tt.get_alpha(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  sg = 1.0 - sw;
  rhogw = tt.get_rhogw(pc,t);
  rhow = tt.get_rhow(t);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);


  //according to Lewis and Schrefler p.394
  capwg = (alpha - phi)/ks*sg*(rhogw*sg + rhow*sw);
  capwg = capwg + sg*phi*mw/t/gasr*dpgw_dpc;
  capwg = capwg + ((alpha - phi)/ks*sg*(rhogw*sg*pc + rhow*sw*pw - rhow*sw*pc) + phi*(rhow - rhogw))*dsw_dpc;

  //printf("capwg = %e\n",capwg);
  capwg = 0.0;

  return(capwg);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kwt - conductivity coefficient
*/
double gmultiph::get_kwt(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double kwt;
 
  kwt = 0.0;

  return(kwt);
}


/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capwt - capacity coefficient
*/
double gmultiph::get_capwt(double pw,double pg,double t,long ipp)
{
  double capwt,pc,betaswg,phi,sw,sg,dsw_dt,rhow,rhogw,dpgw_dt,pgw,alpha,ks;
  state_eq tt;
  
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  betaswg = tt.get_betaswg_c(pc,pg,t,ipp);
  sw = tt.get_s(pc,t,ipp);
  sg = 1.0 - sw;
  phi = tt.get_phi(t,ipp); 
  dpgw_dt = tt.get_dpgw_dt(pc,t); 
  pgw = tt.get_pgw(pc,t); 
  alpha = tt.get_alpha(pc,pg,t,ipp); 
  ks = tt.get_ks(pc,pg,t,ipp); 
  dsw_dt = tt.get_ds_dt(pc,t,ipp);
  rhow = tt.get_rhow(t); 
  rhogw = tt.get_rhogw(pc,t); 
 
  //according to Lewis and Schrefler p.394
  capwt = betaswg*sw + sg*phi*mw/t/gasr*(dpgw_dt - pgw/t);
  capwt = capwt + ((alpha-phi)/ks*(rhogw*sg*pc + rhow*sw*pw - rhow*sw*pc) + phi*(rhow - rhogw))*dsw_dt;

  capwt = 0.0;
  //printf("capwt = %e\n",capwt);

  return(capwt);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgg - conductivity coefficient
*/
double gmultiph::get_kgg(double pw,double pg,double t,long ipp)
{
  double pc,rhoga,kintr,krg,mug,rhog,dg,kgg,mg,pgw,dpgw_dpc;
  state_eq tt;
  
  pc = pg - pw;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  rhoga = tt.get_rhoga(pc,pg,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  rhog = tt.get_rhog(pc,pg,t);
  dg = tt.get_dg(pc,pg,t,ipp);
  mg = tt.get_mg(pc,pg,t);
  pgw = tt.get_pgw(pc,t);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);

  //according to Lewis and Schrefler p.394
  kgg = rhoga*kintr*krg/mug - rhog*ma*mw/mg/mg*dg*(1.0/pg*dpgw_dpc - pgw/pg/pg);

  //printf("kgg = %e\n",kgg);

  return(kgg);
}

/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgg - capacity coefficient
*/
double gmultiph::get_capgg(double pw,double pg,double t,long ipp)
{
  double capgg,pc,sw,phi,alpha,ks,sg,rhoga,dpgw_dpc,dsw_dpc;
  state_eq tt;
 
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  sw = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);
  alpha = tt.get_alpha(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  sg = 1.0 - sw;
  rhoga = tt.get_rhoga(pc,pg,t);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);


  //according to Lewis and Schrefler p.395
  capgg = (alpha - phi)/ks*sg*sg*rhoga - ((alpha - phi)/ks*sg*pc + phi)*rhoga*dsw_dpc;
  capgg = capgg + sg*phi*ma/t/gasr - sg*phi*mw/t/gasr*dpgw_dpc;

  //printf("capgg = %e\n",capgg);

  return(capgg);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgw - conductivity coefficient
*/
double gmultiph::get_kgw(double pw,double pg,double t,long ipp)
{
  double rhog,dg,dpgw_dpc,kgw,mg,pc;
  state_eq tt;
  
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);
  
  rhog = tt.get_rhog(pc,pg,t);
  dg = tt.get_dg(pc,pg,t,ipp);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  mg = tt.get_mg(pc,pg,t);
  
  //according to Lewis and Schrefler p.395
  kgw = rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc;

  kgw = 0.0;
  //printf("kgw = %e\n",kgw);

  return(kgw);
}


/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgw - capacity coefficient
*/
double gmultiph::get_capgw(double pw,double pg,double t,long ipp)
{
  double capgw,pc,sw,phi,alpha,ks,sg,rhoga,dpgw_dpc,dsw_dpc;
  state_eq tt;
 
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  sw = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);
  alpha = tt.get_alpha(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  sg = 1.0 - sw;
  rhoga = tt.get_rhoga(pc,pg,t);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);


  //according to Lewis and Schrefler p.395
  capgw = (alpha - phi)/ks*sw*sg*rhoga + ((alpha - phi)/ks*sg*pc + phi)*rhoga*dsw_dpc;
  capgw = capgw + sg*phi*mw/t/gasr*dpgw_dpc;

  capgw = 0.0;
  //printf("capgw = %e\n",capgw);

  return(capgw);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgt - conductivity coefficient
*/
double gmultiph::get_kgt(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double kgt;

  kgt = 0.0;

  return(kgt);
}


/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgt - capacity coefficient
*/
double gmultiph::get_capgt(double pw,double pg,double t,long ipp)
{
  double capgt,rhoga,betasg,alpha,phi,ks,sw,sg,pc,dsw_dt,dpgw_dt,pgw;
  state_eq tt;
  
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);
  
  rhoga = tt.get_rhoga(pc,pg,t);
  betasg = tt.get_betasg_c(pc,pg,t,ipp);
  alpha = tt.get_alpha(pc,pg,t,ipp);
  phi = tt.get_phi(t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  sw = tt.get_s(pc,t,ipp);
  sg = 1.0 - sw;
  dsw_dt = tt.get_ds_dt(pc,t,ipp);
  dpgw_dt = tt.get_dpgw_dt(pc,t);
  pgw = tt.get_pgw(pc,t);
    
  //according to Lewis and Schrefler p.395
  capgt = -1.0*rhoga*betasg - ((alpha - phi)/ks*sg*pc + phi)*rhoga*dsw_dt;
  capgt = capgt + sg*phi*ma/t/t/gasr + sg*phi*mw/t/gasr*(dpgw_dt - pgw/t);

  capgt = 0.0;
  //printf("capgt = %e\n",capgt);

  return(capgt);
}

/**
   function creates conductivity coefficient of a general material - first part
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt1 - conductivity coefficient
*/
double gmultiph::get_ktt1(double pw,double pg,double t,long ipp)
{
  double lambdaeff,pc,ktt1;
  state_eq tt;
   
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  lambdaeff = tt.get_lambdaeff(pc,pg,t,ipp);

  ktt1 = lambdaeff;

  //printf("ktt = %e\n",ktt1);

  return(ktt1);
}


/**
   function creates conductivity coefficient of a general material - second (A) part (convective term)
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2a - conductivity coefficient
*/
double gmultiph::get_ktt2a(double pw,double pg,double t,long ipp)
{
  double cpw,pc,rhow,kintr,krw,muw,ktt2a,phi,sw;
  state_eq tt;
    
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  cpw = tt.get_cpw();
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  sw = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);

  ktt2a = phi*sw*cpw*rhow*kintr*krw/muw;
  
  return(ktt2a);
}

/**
   function creates conductivity coefficient of a general material - second (B) part (convective term)
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2b - conductivity coefficient
*/
double gmultiph::get_ktt2b(double pw,double pg,double t,long ipp)
{
  double cpg,pc,rhog,kintr,krg,mug,ktt2b,phi,sg;
  state_eq tt;

  pc = pg - pw;
 
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);
  
  cpg = tt.get_cpg(pc,pg,t);
  rhog = tt.get_rhog(pc,pg,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  sg = 1.0 - tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);

  ktt2b = phi*sg*cpg*rhog*kintr*krg/mug;

  return(ktt2b);
}


/**
   function creates conductivity coefficient of a general material - second (C) part (convective term)
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2c - conductivity coefficient
*/
double gmultiph::get_ktt2c(double pw,double pg,double t,long ipp)
{
  double ktt2c,rhow,pc;
  state_eq tt;

  pc = pg - pw;
 
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);
  
  rhow = tt.get_rhow(t);
  ktt2c = rhow;

  return(ktt2c);
}


/**
   function creates conductivity coefficient of a general material - second (D) part (convective term)
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2d - conductivity coefficient
*/
double gmultiph::get_ktt2d(double pw,double pg,double t,long ipp)
{
  double ktt2d,rhog,pc;
  state_eq tt;

  pc = pg - pw;
 
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);
  
  rhog = tt.get_rhog(pc,pg,t);
  ktt2d = rhog;

  return(ktt2d);
}


/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captt - capacity coefficient
*/
double gmultiph::get_captt(double pw,double pg,double t,long ipp)
{
  double captt,pc,rhocp,alpha,betasw,sw,dhvap,rhow,phi,dsw_dt,ks;
  state_eq tt;
 
  pc = pg - pw;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  
  alpha = tt.get_alpha(pc,pg,t,ipp);
  betasw = tt.get_betasw_c(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  rhocp = tt.get_rhocp(pc,pg,t,ipp);
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  dsw_dt = tt.get_ds_dt(pc,t,ipp);
  sw = tt.get_s(pc,t,ipp);

  //according to Lewis and Schrefler p.396
  captt = rhocp;
  captt = captt;// + dhvap*(betasw - rhow*((alpha-phi)/ks*sw*pw - (alpha-phi)/ks*sw*pg + phi))*dsw_dt;

  //printf("captt = %e\n",captt);

  return(captt);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval ktw - conductivity coefficient
*/
double gmultiph::get_ktw(double pw,double pg,double t,long ipp)
{
  double dhvap,rhow,kintr,krw,muw,ktw,pc;
  state_eq tt;

  pc = pg -pw;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);
  
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  
  //according to Lewis and Schrefler p.396
  ktw = -1.0*dhvap*(rhow*kintr*krw/muw);

  ktw = 0.0;
  //printf("ktw = %e\n",ktw);

  return(ktw);
}

/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captw - capacity coefficient
*/
double gmultiph::get_captw(double pw,double pg,double t,long ipp)
{
  double captw,pc,alpha,sw,dhvap,rhow,phi,dsw_dpc,ks,kw;
  state_eq tt;
 
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);


  alpha = tt.get_alpha(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  kw = tt.get_kw();
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);
  sw = tt.get_s(pc,t,ipp);

  //according to Lewis and Schrefler p.396
  captw = -1.0*dhvap*rhow*((alpha-phi)/ks*sw*sw + sw*phi/kw);
  captw = captw + dhvap*rhow*((alpha-phi)/ks*sw*pw - (alpha-phi)/ks*sw*pg + phi)*dsw_dpc;

  captw = 0.0;
  //printf("captw = %e\n",captw);

  return(captw);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktg - conductivity coefficient
*/
double gmultiph::get_ktg(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double ktg;

  ktg = 0.0;

  return(ktg);
}

/**
   function creates capacity coefficient of a general material
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captg - capacity coefficient
*/
double gmultiph::get_captg(double pw,double pg,double t,long ipp)
{
  double captg,pc,alpha,sw,sg,dhvap,rhow,phi,dsw_dpc,ks,kw;
  state_eq tt;
 
  pc = pg - pw;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);


  alpha = tt.get_alpha(pc,pg,t,ipp);
  ks = tt.get_ks(pc,pg,t,ipp);
  kw = tt.get_kw();
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  dsw_dpc = tt.get_ds_dpc(pc,t,ipp);
  sw = tt.get_s(pc,t,ipp);
  sg = 1.0 - sw;

  //according to Lewis and Schrefler p.396
  captg = -1.0*dhvap*rhow*(alpha-phi)/ks*sg*sw;
  captg = captg - dhvap*rhow*((alpha-phi)/ks*sw*pw - (alpha-phi)/ks*sw*pg + phi)*dsw_dpc;

  captg = 0.0;
  //printf("captg = %e\n",captg);

  return(captg);
}

/**
   function creates right-hand side coefficient of a general material for c medium
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double gmultiph::get_fc1(double pw,double pg,double t,long ipp)
{
  double fc1;
  double rhow,rhogw,rhog,kintr,krw,muw,krg,mug;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  rhow = tt.get_rhow(t);
  rhogw = tt.get_rhogw(pw,t);
  rhog = tt.get_rhog(pw,pg,t);  
  kintr = tt.get_kintr(pw,pg,t,ipp);
  krw = tt.get_krw(pw,t,ipp);
  muw = tt.get_muw(t);
  krg = tt.get_krg(pw,t,ipp);
  mug = tt.get_mug(pw,pg,t);

  fc1 = rhogw*kintr*krg*rhog/mug + rhow*kintr*krw*rhow/muw;

  return(fc1);
}


/**
   function creates right-hand side coefficient of a general material for g medium
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double gmultiph::get_fg(double pw,double pg,double t,long ipp)
{
  double fg;
  double rhoga,rhog,kintr,krg,mug;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  rhoga = tt.get_rhoga(pw,pg,t);
  rhog = tt.get_rhog(pw,pg,t);
  kintr = tt.get_kintr(pw,pg,t,ipp);
  krg = tt.get_krg(pw,t,ipp);
  mug = tt.get_mug(pw,pg,t);

  fg = rhoga*kintr*krg*rhog/mug;

  return(fg);
}


/**
   function creates right-hand side coefficient of a general material for t medium
   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double gmultiph::get_ft1(double pw,double pg,double t,long ipp)
{
  double ft1;
  double dhvap,rhow,kintr,krw,muw;
  state_eq tt;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pw,pg,t,ipp);
  krw = tt.get_krw(pw,t,ipp);
  muw = tt.get_muw(t);

  ft1 = -dhvap*rhow*kintr*krw*rhow/muw;

  return(ft1);
}





/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp)
{
  double fc3,pgws,rhow;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //relative humidity - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    
    fc3=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 31:{
    //mass concentration - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    
    fc3=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 32:{
    //water vapour pressure - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    
    fc3=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 33:{//relative humidity
    fc3=0.0;
    break;
  }
  case 34:{//mass concentration
    fc3=0.0;
    break;
  } 
  case 35:{//water vapour pressure
    fc3=0.0;
    break;
  } 
  case 36:{//relative humidity - pokus
    fc3=0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(fc3);
}

/**
   function creates correct new nodal value on the boundary (transmission) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_nodval_ww(double bv,double pw,double pg,double t,long bc,long ipp)
{
  double new_nodval,fc3,pgws,rhow,rhogw_surf,alpha;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //for testing = Kelvin eq. approximated by Taylor's serie in point 0
    //relative humidity -> rhogw
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    alpha = 0.30;
    
    fc3 = bv*pgws*mw/gasr/t;//rel. hum. -> rhogw
    fc3 = fc3 - pgws*mw/gasr/t;//adding of 0th member of T. serie
    fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pw*mw/rhow/t/gasr))*pw*pw/2.0;//adding of residuum of T. serie

    break;
  }
  case 31:{
    //for testing = Kelvin eq. approximated by Taylor's serie in point 0
    // mass concentration (rhogw)
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    alpha = 0.30;

    fc3 = bv - pgws*mw/gasr/t;//adding of 0th member of T. serie
    fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pw*mw/rhow/t/gasr))*pw*pw/2.0;//adding of residuum of T. serie
    
    break;
  }  
  case 32:{
    //for testing = Kelvin eq. approximated by Taylor's serie in point 0
    // water vapour pressure -> rhogw    
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    alpha = 0.30;

    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - pgws*mw/gasr/t;//adding of 0th member of T. serie
    fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pw*mw/rhow/t/gasr))*pw*pw/2.0;//adding of residuum of T. serie
    
    break;
  } 
  case 33:{
    // relative humidity -> rhogw
    /* pgws = tt.get_pgws(t);
       rhogw_surf = tt.get_rhogw(pw,t);
       
       fc3 = bv*pgws*mw/gasr/t;
       
       fc3 = fc3 - rhogw_surf;
    */
    //pokus
    pgws = tt.get_pgws(heat_rate(0.0,Tp->time));
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*pgws*mw/gasr/(heat_rate(0.0,Tp->time));
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  case 34:{
    // mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 35:{
    // water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 36:{
    // relative humidity -> rhogw pokus
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  new_nodval = fc3;
  return(new_nodval);
}



/**
   function creates flux on the boundary (transmission - convective mass transfer) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_flux_ww(double bv,double pw,double pg,double t,long bc,long ipp)
{
  double flux,fc3,pgws,rhogw_surf;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//for testing - boundary flux
    //relative humidity -> rhogw
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  case 31:{
    // mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 32:{
    // water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 33:{//for testing - boundary flux
    //relative humidity -> rhogw
    /*  pgws = tt.get_pgws(t);
	rhogw_surf = tt.get_rhogw(pw,t);
	
	fc3 = bv*pgws*mw/gasr/t;
	
	fc3 = fc3 - rhogw_surf;
    */
    //pokus
    pgws = tt.get_pgws(heat_rate(0.0,Tp->time));
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*pgws*mw/gasr/(heat_rate(0.0,Tp->time));
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  case 34:{
    //mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 35:{
    //water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 36:{//pokus
    //pokus
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pw,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
 default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  flux = fc3;
  return(flux);
}


/**
   function creates correct transfer coefficient on the boundary (transmission) for t medium

   @param pw - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_transcoeff_tt(double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    ft3=1.0;
    break;
  }
  case 31:{//heat transmission for testing
    ft3=0.0;
    break;
  }
  case 32:{//heat transmission
    ft3=1.0;
    break;
  }
  case 90:{//heat transfer + radiation(fire)
    ft3=0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  
  return(ft3);
}




/**
   function creates correct new nodal value on the boundary (transmission) for t medium

   @param bv - value of prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_nodval_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    ft3 = bv;
    break;
  }
  case 31:{//heat transmission for testing
    ft3 = (bv - t);
    break;
  }
  case 32:{//heat transmission - heat rate
    ft3 = heat_rate(0.0,Tp->time);
    break;
  }    
  case 90:{//heat transfer + radiation(fire)
    ft3 = (bv - t) + trr*(bv*bv*bv*bv - t*t*t*t);//(trr = e_sigma0/beta_t)
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(ft3);
}


/**
   function creates flux on the boundary (transmission - convective mass transfer) for c medium

   @param bv - prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double gmultiph::get_transmission_flux_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //capillary pressure check
  cappress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission - boundary flux
    ft3 = (bv - t);
    break;
  }
  case 31:{//heat transmission for testing - boundary flux
    ft3 = (bv - t);
    break;
  }
  case 32:{//heat transmission - boundary flux
    ft3 = (heat_rate(0.0,Tp->time) - t);
    break;
  }
  case 90:{//heat transfer + radiation(fire)
    ft3 = (bv - t) + trr*(bv*bv*bv*bv - t*t*t*t);//(trr = e_sigma0/beta_t)
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }

  return(ft3);
}

/**
   function computes temperature curve of heating rate up to 300deg. C
   @param stime - starting time
   @param time - current time
*/

double gmultiph::heat_rate(double /*stime*/,double time)
{
  double tinf;

  tinf = 298.15 + time/60.0;

  if (tinf >= 573.15)
    tinf = 573.15;

  return(tinf);
}


/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pw - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double gmultiph::get_othervalue(long compother,long ipp,double *r)
{
  double other,pc;
  state_eq tt;

  pc = r [1] - r[0];

  switch (compother){
  case 0:{//water pressure
    other = r[0];
      break;
  }
  case 1:{//gas pressure
    other = r[1];
    break;
  }
  case 2:{//temperature
    other = r[2];
      break;
  }
  case 3:{//relative humidity
    other = tt.get_rh(pc,r[2]);
      break;
  }
  case 4:{//saturation
    other = tt.get_s(pc,r[2],ipp);
    break;
  }
  case 5:{//vapour pressure
    other = tt.get_pgw(pc,r[2]);
    break;
  }
  case 6:{//capillary pressure
    other = pc;
    break;
  }
  case 7:{//moisture content
    other = tt.get_w(pc,r[1],r[2],ipp);
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);

}

/**
     function prints names of all variables in nodes
     @param out - output file
     @param compother - number of other components
*/
void gmultiph::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  case 1:{//gas pressure
    fprintf (out,"Gas pressure (Pa)             ");
    break;
  }
  case 2:{//temperature
    fprintf (out,"Temperature (K)               ");
    break;
  }
  case 3:{//relative humidity
    fprintf (out,"Relative humidity ()          ");
    break;
  }
  case 4:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 5:{//vapour pressure
    fprintf (out,"Pore water vapor pressure (Pa)");
    break;
  }
  case 6:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
    break;
  }
  case 7:{//moisture content
    fprintf (out,"Moisture content (kg/kg)      ");
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}
