/*
  File:             multiphase.cpp
  Author:           Tomas Krejci, 30.3.2003
  Purpose:          computes conductivity and capacity matrices in a material point for multiphase porous media
  of rigid prous materials
  
  FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
  --------------------------------------------------
  
  sources: 
  1. FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS BY AN ALGEBRAIC MULTIGRID METHOD
  Wang Xicheng, B.A. Schrefler (there are few mistakes)
  
  2. NUMERICAL ANALYSIS OF HYGRO-THERMAL BEHAVIOUR AND DAMAGE OF CONCRETE AT HIGH TEMPERATURE
  D. Gawin, C.E. Majorana, B.A. Schrefler
  
  3. NONLINEAR MODELLING OF CONCRETE AS MULTIPHASE POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
  Francesco Pesavento - doctoral thesis                      
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "multiphase.h"
#include "constrel.h"
#include "globalt.h"

multiph::multiph()
{
  mw = 18.01528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  scale_pc = Tp->scale[0];
  scale_pg = Tp->scale[1];
  scale_t = Tp->scale[2];
//  scale_pc = Tp->scale1;
//  scale_pg = Tp->scale2;
//  scale_t = Tp->scale3;
}
multiph::~multiph()
{}


/**
   JK, 11.4.2019
*/
double multiph::give_pg (long nn)
{
  long k;
  double pg;
  
  k=Gtt->give_dof(nn,1);
  if (k>0)   {pg = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {pg = 0.0;}
  if (k<0)   {pg = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return pg;
}

/**
   JK, 11.4.2019
*/
double multiph::give_pc (long nn)
{
  long k;
  double pc;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {pc = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {pc = 0.0;}
  if (k<0)   {pc = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return pc;
}

/**
   JK, 11.4.2019
*/
double multiph::give_temp (long nn)
{
  long k;
  double t;
  
  k=Gtt->give_dof(nn,2);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return t;
}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void multiph::matcond (matrix &d,long ri,long ci,long ipp)
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
void multiph::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kcg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kct(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pc,pg,t,ipp);// *scale_t;//scaling
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void multiph::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kcg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kct(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pc,pg,t,ipp);// *scale_t;//scaling
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;

  //fprintf(Outt,"\n i=%ld j=%ld kij=%e",ri,ci,k);//debug information

}

/**
   function creates conductivity matrix of the material for 3D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void multiph::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kcg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kct(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pc,pg,t,ipp);// *scale_t;//scaling
  
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
void multiph::matcap (double &c,long ri,long ci,long ipp)
{
  double pc,pg,t;
  c = 0.0;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capcc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    c = get_capcg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    c = get_capct(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    c = get_capgc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 1) && (ci == 1))
    c = get_capgg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    c = get_capgt(pc,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    c = get_captc(pc,pg,t,ipp);// *scale_pc;//scaling
  if((ri == 2) && (ci == 1))
    c = get_captg(pc,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    c = get_captt(pc,pg,t,ipp);// *scale_t;//scaling
}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void multiph::matcond2 (matrix &d,long ri,long ci,long ipp)
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
void multiph::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pc,pg,t;
      
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }

  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pc,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pc,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pc,pg,t,ipp);// *scale_t;//scaling

    d[0][0] = -1.0*(a*(Tm->ip[ipp].grad[1][0] - Tm->ip[ipp].grad[0][0] - b*Tp->gr[0]) + c*Tm->ip[ipp].grad[1][0]);
  }  
}


/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void multiph::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pc,pg,t;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pc,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pc,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pc,pg,t,ipp);// *scale_t;//scaling

    d[0][0] = -1.0*(a*(Tm->ip[ipp].grad[1][0] - Tm->ip[ipp].grad[0][0] - b*Tp->gr[0]) + c*Tm->ip[ipp].grad[1][0]);
    d[0][1] = -1.0*(a*(Tm->ip[ipp].grad[1][1] - Tm->ip[ipp].grad[0][1] - b*Tp->gr[1]) + c*Tm->ip[ipp].grad[1][1]);
  }
}



/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void multiph::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pc,pg,t;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
    
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pc,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pc,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pc,pg,t,ipp);// *scale_t;//scaling

    d[0][0] = -1.0*(a*(Tm->ip[ipp].grad[1][0] - Tm->ip[ipp].grad[0][0] - b*Tp->gr[0]) + c*Tm->ip[ipp].grad[1][0]);
    d[0][1] = -1.0*(a*(Tm->ip[ipp].grad[1][1] - Tm->ip[ipp].grad[0][1] - b*Tp->gr[1]) + c*Tm->ip[ipp].grad[1][1]);
    d[0][2] = -1.0*(a*(Tm->ip[ipp].grad[1][2] - Tm->ip[ipp].grad[0][2] - b*Tp->gr[2]) + c*Tm->ip[ipp].grad[1][2]);
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
void multiph::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
void multiph::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pc,pg,t;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fc1(pc,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 1){
    f = get_fg(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 2){
    f = get_ft1(pc,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void multiph::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pc,pg,t;
  
  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fc1(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 1){
    f = get_fg(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 2){
    f = get_ft1(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void multiph::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pc,pg,t;

  pc = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling

  if(ri == 0){
    f = get_fc1(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 1){
    f = get_fg(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 2){
    f = get_ft1(pc,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
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
double multiph::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pc,pg,t;
  c = 0.0;
  
  pc = give_pc (nn);
  pg = give_pg (nn);
  t = give_temp (nn);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_cc(pc,pg,t,bc,ipp);// *scale_pc;//scaling
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
    c = get_transmission_transcoeff_tt(pc,pg,t,bc,ipp);// *scale_t;//scaling

  c = c*trc;

  return (c);
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
double multiph::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
{
  double c,pc,pg,t;
  c = 0.0;
  
  pc = give_pc (nn);
  pg = give_pg (nn);
  t = give_temp (nn);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_cc(pc,pg,t,bc,ipp,flag);// *scale_pc;//scaling
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
    c = get_transmission_transcoeff_tt(pc,pg,t,bc,ipp);// *scale_t;//scaling

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
double multiph::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pc,pg,t;
  c = 0.0;
  
  pc = give_pc (nn);
  pg = give_pg (nn);
  t = give_temp (nn);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_nodval_cc(nodval,pc,pg,t,bc,ipp);// *scale_pc;//scaling
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
    c = get_transmission_nodval_tt(nodval,trc2,pc,pg,t,bc,ipp);// *scale_t;//scaling

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
double multiph::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pc,pg,t;
  c = 0.0;
  
  pc = give_pc (nn);
  pg = give_pg (nn);
  t = give_temp (nn);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_flux_cc(nodval,pc,pg,t,bc,ipp);// *scale_pc;//scaling
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
    c = get_transmission_flux_tt(nodval,trc2,pc,pg,t,bc,ipp);// *scale_t;//scaling

  return (c);
}


/**
   function checks if gas pressure is greater than vapour pressure

   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr, 18.5.2005
*/
void multiph::gaspress_check(double pc,double &pg,double t,long /*ipp*/)
{
  double pgw;
  state_eq tt;
  
  //general
  //gas pressure check
  pgw = tt.get_pgw(pc,t);
  if (pgw > pg) 
    pg = pgw;

  
  //general
  //fully saturated state
  if(pc < 100.0)
    pg = 101325.0;
  
  // M. Starnoni 24-11-2010
  if (pg<0.0){
    pg = 0.0;
    print_err("gas pressure was modified",__FILE__,__LINE__,__func__);
  }
  
}


/**
   function checks if capillary pressure is higher than maximum

   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void multiph::cappress_check(double &pc,double /*pg*/,double /*t*/,long /*ipp*/)
{
      state_eq tt;

  /*  double pgw,pcmax,rhow;
      state_eq tt;
      
      //only boundary nodes
      //maximum capillary pressure check
      pgw = tt.get_pgw(pc,t);
      rhow = tt.get_rhow(t);
      pcmax = 0.0;
      if (pg < pgw && pg > 0.0)
      pcmax = -1.0*rhow*gasr*t/mw*log(pg/pgw);
      else
      pcmax = 0.0;
      
      if (pg < pgw){
      if (pcmax > 0.0){
      pc = pcmax;
      }
      }
      
      
      if (pc > pcmax)
      pc = pcmax;
      
      //storing into gauss point
      Tm->ip[ipp].av[0] = pc;
  */

  // M. Starnoni 24-11-2010
  //corrected by TKr 08/12/2022
  if (pc <= tt.pcmin){
    pc = tt.pcmin;
    //print_err("capillary pressure was modified",__FILE__,__LINE__,__func__);
  }
  
}


/**
   function fixed capillary pressure if on one element is temperature higher than critical temperature for water
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void multiph::cappress_stop(double &/*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  //storing into gauss point
  //Tm->ip[ipp].av[0] = pc;
}

/**
   function corrects values of variables
   
   @param nv - array with variables
   
*/
void multiph::values_correction (vector &nv)
{
  //  capillary pressure
  cappress_check(nv[0],nv[1],nv[2],0);
  
  //  gas pressure
  gaspress_check(nv[0],nv[1],nv[2],0);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kcc - conductivity coefficient
*/
double multiph::get_kcc(double pc,double pg,double t,long ipp)
{
  double rhow,kintr,krw,muw,rhog,deff,dpgw_dpc,kcc,mg,ddbw,ds_dpc;
  state_eq tt;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  rhog = tt.get_rhog(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  mg = tt.get_mg(pc,pg,t);
  ds_dpc = tt.get_ds_dpc(pc,t,ipp);
  ddbw = tt.get_ddbw(pc,pg,t,ipp);

  //kcc = rhow*kintr*krw/muw + ddbw*rhow;//water
  //kcc = kcc - deff*dpgw_dpc/t/mg;//water vapour
  //kcc = -1.0*kcc;//sign correction
  //correction by TKr 09/12/2022
  kcc = -rhow*kintr*krw/muw - ddbw*rhow + rhog*ma*mw/mg/mg*deff/pg*dpgw_dpc;//water + bound water + water vapour

  return(kcc);
}

/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capcc - capacity coefficient
*/
double multiph::get_capcc(double pc,double pg,double t,long ipp)
{
  double capcc,s,drhogw_dpc,rhogw,drhow_dpc,rhow,phi,ds_dpc;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  ds_dpc = tt.get_ds_dpc(pc,t,ipp);
  s = tt.get_s(pc,t,ipp);
  drhow_dpc = 0.0;
  drhogw_dpc = tt.get_drhogw_dpc(pc,t);
  rhogw = tt.get_rhogw(pc,t);

  capcc = rhow*phi*ds_dpc + phi*s*drhow_dpc;//water
  capcc = capcc + (1.0 - s)*phi*drhogw_dpc - phi*rhogw*ds_dpc;//water vapour

  return(capcc);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kcg - conductivity coefficient
*/
double multiph::get_kcg(double pc,double pg,double t,long ipp)
{
  double rhogw,kintr,krg,mug,deff,pgw,rhow,krw,muw,kcg,mg;
  double rhog;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhogw = tt.get_rhogw(pc,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  pgw = tt.get_pgw(pc,t);
  rhow = tt.get_rhow(t);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  mg = tt.get_mg(pc,pg,t);
  rhog = tt.get_rhog(pc,pg,t);

  //kcg = -1.0*rhow*kintr*krw/muw;//water
  //kcg = -1.0*rhogw*kintr*krg/mug + deff*pgw/pg/t/mg;//water vapour
  //kcg = -1.0*kcg;//sign correction
  //correction by TKr 09/12/2022
  kcg = rhow*kintr*krw/muw + rhogw*kintr*krg/mug - rhog*ma*mw/mg/mg*deff*pgw/pg/pg;//water + water vapour

  return(kcg);
}

/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capcg - capacity coefficient
*/
double multiph::get_capcg(double pc,double pg,double t,long ipp)
{
  double capcg,s,phi,drhow_dpg;
  state_eq tt;
 
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  s = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);
  drhow_dpg = 0.0;

  capcg = phi*s*drhow_dpg;//water

  return(capcg);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kct - conductivity coefficient
*/
double multiph::get_kct(double pc,double pg,double t,long ipp)
{
  double rhog,deff,dpgw_dt,kct,mg;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhog = tt.get_rhog(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  dpgw_dt = tt.get_dpgw_dt(pc,t);
  mg = tt.get_mg(pc,pg,t);

  //kct = 0.0;//water
  //kct = kct - deff*dpgw_dt/t/mg;//water vapour
  //  kct = -1.0*kct;//sign correction
  //correction by TKr 09/12/2022
  kct = rhog*ma*mw/mg/mg*deff/pg*dpgw_dt;

  return(kct);
}


/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capct - capacity coefficient
*/
double multiph::get_capct(double pc,double pg,double t,long ipp)
{
  double capct,betasw,drhow_dt,dphi_dt,dehydw_dt,phi,s,drhogw_dt,ds_dt,rhow,rhogw;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  phi = tt.get_phi(t,ipp); 
  s = tt.get_s(pc,t,ipp); 
  drhogw_dt = tt.get_drhogw_dt(pc,t); 
  ds_dt = tt.get_ds_dt(pc,t,ipp); 
  rhow = tt.get_rhow(t); 
  rhogw = tt.get_rhogw(pc,t); 
  betasw = 0.0;//provisionally
  drhow_dt = tt.get_drhow_dt(pc,t); 
  dphi_dt = tt.get_dphi_dt(pc,pg,t,ipp); 
  dehydw_dt = tt.get_dehydw_dt(pc,pg,t,ipp);
 
  capct = betasw*rhow + rhow*phi*ds_dt + phi*s*drhow_dt + dphi_dt*s*rhow - dehydw_dt;//water
  capct = capct + (1-s)*phi*drhogw_dt - phi*rhogw*ds_dt + dphi_dt*(1.0 -s)*rhogw;//water vapour

  return(capct);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgg - conductivity coefficient
*/
double multiph::get_kgg(double pc,double pg,double t,long ipp)
{
  double rhoga,kintr,krg,mug,rhog,deff,kgg,mg,pgw;
  double dpgw_dpc;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhoga = tt.get_rhoga(pc,pg,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  rhog = tt.get_rhog(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  mg = tt.get_mg(pc,pg,t);
  pgw = tt.get_pgw(pc,t);

  //kgg = -1.0*rhoga*kintr*krg/mug - deff*pgw/pg/t/mg;//dry air + water vapour
  //kgg = -1.0*kgg;//sign correction

  // correction 09/12/2022:
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  kgg = rhoga*kintr*krg/mug - rhog*ma*mw/mg/mg*deff*(dpgw_dpc/pg-pgw/pg/pg);

  return(kgg);
}

/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgg - capacity coefficient
*/
double multiph::get_capgg(double pc,double pg,double t,long ipp)
{
  double capgg,phi,s,drhoga_dpg;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  phi = tt.get_phi(t,ipp);
  s = tt.get_s(pc,t,ipp);
  drhoga_dpg = tt.get_drhoga_dpg(pc,pg,t);
  
  capgg = phi*(1-s)*drhoga_dpg;//dry air

  return(capgg);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgc - conductivity coefficient
*/
double multiph::get_kgc(double pc,double pg,double t,long ipp)
{
  double rhog,deff,dpgw_dpc,kgc,mg;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);
  
  rhog = tt.get_rhog(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  dpgw_dpc = tt.get_dpgw_dpc(pc,t);
  mg = tt.get_mg(pc,pg,t);
  
  //kgc = deff*dpgw_dpc/t/mg;//water vapour
  //kgc = -1.0*kgc;//sign correction
  //correction 09/12/2022 by TKr
  kgc = -rhog*ma*mw/mg/mg*deff/pg*dpgw_dpc;

  return(kgc);
}


/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgc - capacity coefficient
*/
double multiph::get_capgc(double pc,double pg,double t,long ipp)
{
  double capgc,phi,s,drhoga_dpc,rhoga,ds_dpc;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  phi = tt.get_phi(t,ipp);
  s = tt.get_s(pc,t,ipp);
  drhoga_dpc = tt.get_drhoga_dpc(pc,pg,t);
  rhoga = tt.get_rhoga(pc,pg,t);
  ds_dpc = tt.get_ds_dpc(pc,t,ipp);
  
  capgc = phi*(1-s)*drhoga_dpc - phi*rhoga*ds_dpc;//dry air
  
  return(capgc);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgt - conductivity coefficient
*/
double multiph::get_kgt(double pc,double pg,double t,long ipp)
{
  double rhog,deff,dpgw_dt,kgt,mg;
  state_eq tt;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhog = tt.get_rhog(pc,pg,t);
  deff = tt.get_deff(pc,pg,t,ipp);
  dpgw_dt = tt.get_dpgw_dt(pc,t);
  mg = tt.get_mg(pc,pg,t);

  
  //kgt = deff*dpgw_dt/t/mg;//water vapour
  //kgt = -1.0*kgt;//sign correction
  //correction by TKr 09/12/2022
  kgt = -rhog*ma*mw/mg/mg*deff/pg*dpgw_dt;

  return(kgt);
}


/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgt - capacity coefficient
*/
double multiph::get_capgt(double pc,double pg,double t,long ipp)
{
  double capgt,phi,s,drhoga_dt,rhoga,ds_dt,dphi_dt;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  phi = tt.get_phi(t,ipp);
  s = tt.get_s(pc,t,ipp);
  drhoga_dt = tt.get_drhoga_dt(pc,pg,t);
  rhoga = tt.get_rhoga(pc,pg,t);
  ds_dt = tt.get_ds_dt(pc,t,ipp);
  dphi_dt = tt.get_dphi_dt(pc,pg,t,ipp);

  capgt = phi*(1.0-s)*drhoga_dt - phi*rhoga*ds_dt + dphi_dt*(1.0 - s)*rhoga;//dry air

  return(capgt);
}

/**
   function creates conductivity coefficient of a general material - first part
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt1 - conductivity coefficient
*/
double multiph::get_ktt1(double pc,double pg,double t,long ipp)
{
  double lambdaeff,ktt1;
  state_eq tt;
   
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  lambdaeff = tt.get_lambdaeff(pc,pg,t,ipp);

  ktt1 = lambdaeff;

  return(ktt1);
}


/**
   function creates conductivity coefficient of a general material - second (A) part (convective term)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2a - conductivity coefficient
*/
double multiph::get_ktt2a(double pc,double pg,double t,long ipp)
{
  double cpw,rhow,kintr,krw,muw,ktt2a,phi,s;
  state_eq tt;
    
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  cpw = tt.get_cpw();
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  phi = tt.get_phi(t,ipp);
  s = tt.get_s(pc,t,ipp);

  //correction 09/12/2022
  ktt2a = phi*s*cpw*rhow*kintr*krw/muw;
  
  return(ktt2a);
}

/**
   function creates conductivity coefficient of a general material - second (B) part (convective term)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2b - conductivity coefficient
*/
double multiph::get_ktt2b(double pc,double pg,double t,long ipp)
{
  double rhow,ktt2b;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);
  
  rhow = tt.get_rhow(t);
  
  ktt2b = rhow;

  return(ktt2b);
}

/**
   function creates conductivity coefficient of a general material - second (C) part (convective term)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2c - conductivity coefficient
*/
double multiph::get_ktt2c(double pc,double pg,double t,long ipp)
{
  double rhocpg,kintr,krg,mug,ktt2c,phi,s;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhocpg = tt.get_rhocpg(pc,pg,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);
  phi = tt.get_phi(t,ipp);
  s = tt.get_s(pc,t,ipp);

  ktt2c = phi*s*rhocpg*kintr*krg/mug;

  return(ktt2c);
}

/**
   function creates conductivity coefficient of a general material - second (D) part (convective term)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2d - conductivity coefficient
*/
double multiph::get_ktt2d(double pc,double pg,double t,long ipp)
{
  double rhog,ktt2d;
  state_eq tt;
    
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhog = tt.get_rhog(pc,pg,t);
  
  ktt2d = rhog;

  return(ktt2d);
}


/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captt - capacity coefficient
*/
double multiph::get_captt(double pc,double pg,double t,long ipp)
{
  double captt,s,drhow_dt,dphi_dt,dehydw_dt,fste,hydren,rhocp,dhvap,rhow,phi,ds_dt,betasw;
  state_eq tt;
   
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhocp = tt.get_rhocp(pc,pg,t,ipp);
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  ds_dt = tt.get_ds_dt(pc,t,ipp);
  //betasw = tt.get_betasw(pc,pg,t,ipp);
  betasw = 0.0;//provisionally
  s = tt.get_s(pc,t,ipp);
  drhow_dt = tt.get_drhow_dt(pc,t);
  dphi_dt = tt.get_dphi_dt(pc,pg,t,ipp);
  dehydw_dt = tt.get_dehydw_dt(pc,pg,t,ipp);
  fste = tt.get_fste(pc,pg,t,ipp);
  hydren = tt.get_hydren(pc,pg,t,ipp);

  captt = rhocp - dhvap*(betasw*rhow + rhow*phi*ds_dt + phi*s*drhow_dt + dphi_dt*s*rhow - dehydw_dt) + dehydw_dt*hydren/fste;//water

  return(captt);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval ktc - conductivity coefficient
*/
double multiph::get_ktc(double pc,double pg,double t,long ipp)
{
  double dhvap,ddbw,rhow,kintr,krw,muw,ktc,ds_dpc;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);
  
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  ddbw = tt.get_ddbw(pc,pg,t,ipp);
  ds_dpc = tt.get_ds_dpc(pc,t,ipp);
  
  ktc = dhvap*(rhow*kintr*krw/muw + ddbw*rhow);//water

  return(ktc);
}

/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captc - capacity coefficient
*/
double multiph::get_captc(double pc,double pg,double t,long ipp)
{
  double captc,s,drhow_dpc,dhvap,rhow,phi,ds_dpc;
  state_eq tt;
 
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  phi = tt.get_phi(t,ipp);
  ds_dpc = tt.get_ds_dpc(pc,t,ipp);
  s = tt.get_s(pc,t,ipp);
  drhow_dpc = 0.0;

  captc = -1.0*dhvap*(rhow*phi*ds_dpc + phi*s*drhow_dpc);//water

  return(captc);
}

/**
   function creates conductivity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktg - conductivity coefficient
*/
double multiph::get_ktg(double pc,double pg,double t,long ipp)
{
  double dhvap,rhow,kintr,krw,muw,ktg;
  state_eq tt;
    
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);

  ktg = -1.0*(dhvap*rhow*kintr*krw/muw);//water

  return(ktg);
}

/**
   function creates capacity coefficient of a general material
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captg - capacity coefficient
*/
double multiph::get_captg(double pc,double pg,double t,long ipp)
{
  double captg,s,phi,drhow_dpg,dhvap;
  state_eq tt;
 
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  s = tt.get_s(pc,t,ipp);
  phi = tt.get_phi(t,ipp);
  drhow_dpg = 0.0;
  dhvap = tt.get_dhvap(t);

  captg = -1.0*dhvap*phi*s*drhow_dpg;//water

  return(captg);
}

/**
   function creates right-hand side coefficient of a general material for c medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double multiph::get_fc1(double pc,double pg,double t,long ipp)
{
  double fc1;
  double rhow,rhogw,rhog,kintr,krw,muw,krg,mug;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhow = tt.get_rhow(t);
  rhogw = tt.get_rhogw(pc,t);
  rhog = tt.get_rhog(pc,pg,t);  
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);

  fc1 = rhogw*kintr*krg*rhog/mug + rhow*kintr*krw*rhow/muw;

  return(fc1);
}


/**
   function creates right-hand side coefficient of a general material for g medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double multiph::get_fg(double pc,double pg,double t,long ipp)
{
  double fg;
  double rhoga,rhog,kintr,krg,mug;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  rhoga = tt.get_rhoga(pc,pg,t);
  rhog = tt.get_rhog(pc,pg,t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krg = tt.get_krg(pc,t,ipp);
  mug = tt.get_mug(pc,pg,t);

  fg = rhoga*kintr*krg*rhog/mug;

  return(fg);
}


/**
   function creates right-hand side coefficient of a general material for t medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double multiph::get_ft1(double pc,double pg,double t,long ipp)
{
  double ft1;
  double dhvap,rhow,kintr,krw,muw;
  state_eq tt;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = tt.get_kintr(pc,pg,t,ipp);
  krw = tt.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);

  ft1 = -dhvap*rhow*kintr*krw*rhow/muw;

  return(ft1);
}





/**
   function creates correct transfer coefficient on the boundary (transmission) for pc

   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_transcoeff_cc(double pc,double pg,double t,long bc,long ipp)
{
  double fc3,pgws,rhow;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

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
  case 33:{//relative humidity - flux
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
  case 90:{//fire = heat transfer + radiation
    fc3 = 1.0;
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
   function creates correct transfer coefficient on the boundary (transmission) for pc

   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_transcoeff_cc(double pc,double pg,double t,long bc,long ipp,int flag)
{
  double fc3,pgws,rhow;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //relative humidity - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    //pgws = tt.get_pgws(t);
    //rhow = tt.get_rhow(t);

    fc3 = 1.0;//tt.get_pcrh(bv,t);//tady??!!
    //fc3=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie

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
  case 33:{//relative humidity - flux
    if(flag == 1)
      fc3=0.0;//into right hand side (matrix)
    else{
      pgws = tt.get_pgws(t);
      fc3=1.0;//into left hand side (flux)
      //fc3=gasr*t/mw*pgws;//pokus
    }
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
  case 90:{//fire = heat transfer + radiation
    fc3 = 1.0;
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
   function creates correct new nodal value on the boundary (transmission) for pc

   @param bv - prescribed value near the boundary
   @param pc - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_nodval_cc(double bv,double pc,double pg,double t,long bc,long ipp)
{
  double new_nodval,fc3,pgws,rhow,rhogw_surf,alpha,temp;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //for testing = Kelvin eq. approximated by Taylor's serie in point 0
    //relative humidity -> rhogw
    //pgws = tt.get_pgws(t);
    //rhow = tt.get_rhow(t);
    //alpha = 0.30;
    
    fc3 = tt.get_pcrh(bv,t);

    //fc3 = bv*pgws*mw/gasr/t;//rel. hum. -> rhogw
    //fc3 = fc3 - pgws*mw/gasr/t;//adding of 0th member of T. serie
    //fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pc*mw/rhow/t/gasr))*pc*pc/2.0;//adding of residuum of T. serie

    break;
  }
  case 31:{
    //for testing = Kelvin eq. approximated by Taylor's serie in point 0
    // mass concentration (rhogw)
    pgws = tt.get_pgws(t);
    rhow = tt.get_rhow(t);
    alpha = 0.30;

    fc3 = bv - pgws*mw/gasr/t;//adding of 0th member of T. serie
    fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pc*mw/rhow/t/gasr))*pc*pc/2.0;//adding of residuum of T. serie
    
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
    fc3 = fc3 - mw*mw*mw*pgws/t/t/t/gasr/gasr/gasr/rhow/rhow*(exp(-alpha*pc*mw/rhow/t/gasr))*pc*pc/2.0;//adding of residuum of T. serie
    
    break;
  } 
  case 33:{//relative humidity - flux
    // relative humidity -> rhogw
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*pgws*mw/gasr/t;//surface temperature is temporarily concidered 
    
    fc3 = fc3 - rhogw_surf;
    
    break;
  }
  case 34:{
    // mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 35:{
    // water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 36:{
    // relative humidity -> rhogw pokus
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  case 90:{//prescribed pgw for fire = heat transfer + radiation
    temp = iso_fire(Tp->time,293.0,0.0);
    rhogw_surf = tt.get_rhogw(pc,t);

    fc3 = rhogw_surf - bv/temp/gasr*mw; //flux for rhogw is computing
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
   function creates flux on the boundary (transmission - convective mass transfer) for pc

   @param bv - prescribed value near the boundary
   @param pc - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_flux_cc(double bv,double pc,double pg,double t,long bc,long ipp)
{
  double flux,fc3,pgws,rhogw_surf,temp;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//for testing - boundary flux
    //relative humidity -> rhogw
    //pgws = tt.get_pgws(t);
    //rhogw_surf = tt.get_rhogw(pc,t);
    
    //fc3 = bv*pgws*mw/gasr/t;
    
    //fc3 = fc3 - rhogw_surf;

    fc3 = tt.get_pcrh(bv,t);
    fc3 = fc3 - pc;


    break;
  }
  case 31:{
    // mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 32:{
    // water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 33:{//relative humidity - flux
    //relative humidity -> rhogw
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;

    break;
  }
  case 34:{
    //mass concentration (rhogw)
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv - rhogw_surf;
    break;
  }
    
  case 35:{
    //water vapour pressure -> rhogw
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*mw/gasr/t;
    fc3 = fc3 - rhogw_surf;
    break;
  } 
  case 36:{//pokus
    //pokus
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pc,t);
    
    fc3 = bv*pgws*mw/gasr/t;
    
    fc3 = fc3 - rhogw_surf;
    break;
  }
  case 90:{//prescribed pgw for fire = heat transfer + radiation
    temp = iso_fire(Tp->time,293.0,0.0);
    pgws = tt.get_pgws(t);
    rhogw_surf = tt.get_rhogw(pc,t);

    fc3 = bv/temp/gasr*mw;
    fc3 = rhogw_surf - fc3;//flux for rhogw is computing
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

   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_transcoeff_tt(double pc,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

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
  case 90:{//fire = heat transfer + radiation
    ft3=1.0;
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
   @param pc - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_nodval_tt(double bv,double trr,double pc,double pg,double t,long bc,long ipp)
{
  double ft3,temp;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    ft3 = bv;
    break;
  }
  case 31:{//heat transmission for testing
    ft3 = (bv - t);
    ft3 = -1.0*ft3;
    break;
  }
  case 32:{//heat transmission - heat rate
    ft3 = heat_rate(0.0,Tp->time);
    break;
  }    
  case 90:{//fire = heat transfer + radiation
    temp = iso_fire(Tp->time,293.0,0.0);
    ft3 = (temp - t) + trr*(temp*temp*temp*temp - t*t*t*t);//(trr = e_sigma0/beta_t)
    ft3 = -1.0*ft3;
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
   @param pc - actual capillary pressure on the boundary
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double multiph::get_transmission_flux_tt(double bv,double trr,double pc,double pg,double t,long bc,long ipp)
{
  double ft3,temp;

  //gas pressure check
  gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  cappress_check(pc,pg,t,ipp);

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
    temp = iso_fire(Tp->time,293.0,0.0);
    ft3 = (temp - t) + trr*(temp*temp*temp*temp - t*t*t*t);//(trr = e_sigma0/beta_t)
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

double multiph::heat_rate(double /*stime*/,double time)
{
  double tinf;

  tinf = 298.15 + time/60.0;

  if (tinf >= 573.15)
    tinf = 573.15;

  return(tinf);
}

 
/**  function computes temperature of fire curve ISO 834
     @param time   - current time
     @param t0     - reference temperature
     @retval tfirestart  - fire starting time
*/
double multiph::iso_fire (double time, double t0, double tfirestart)
{
  double tinf;

  //ISO 834 Fire Curve
  tinf=((t0-273.15)+345.0*log10(8.0*(time-tfirestart)/60.0+1.0))+273.15;  
  
  return(tinf);
}
  

/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pc - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double multiph::get_othervalue(long compother,long ipp,double *r)
{
  double other;
  state_eq tt;

  switch (compother){
  case 0:{//capillary pressure
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
    other = tt.get_rh(r[0],r[2]);
      break;
  }
  case 4:{//saturation
    other = tt.get_s(r[0],r[2],ipp);
    break;
  }
  case 5:{//vapour pressure
    other = tt.get_pgw(r[0],r[2]);
    break;
  }
  case 6:{//liquid water pressure
    other = tt.get_pw(r[0],r[1],r[2]);
    break;
  }
  case 7:{//moisture content
    other = tt.get_w(r[0],r[1],r[2],ipp);
    break;
  }    
  case 8:{//temp??!!
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
void multiph::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
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
  case 6:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  case 7:{//moisture content
    fprintf (out,"Moisture content (kg/kg)      ");
    break;
  }    
  case 8:{//temp??!!
    fprintf (out,"temp      ");
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}
