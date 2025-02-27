#include "difcalc.h"
#include "basefun.h"
#include "ordering.h"
#include <math.h>
#include <stdlib.h>

/**
  The function evaluates derivates of the function with variable xi
  with respect to x

  @param dx - array containing derivatives with respect to x
  @param jac - Jacobian
  @param x - array containing node coordinates
  @param xi - natural coordinate
   
  1.4.2002
*/
void derivatives_1d (vector &dx, double &jac, vector &x, double xi)
{
  vector nx(ASTCKVEC(x.n));
  
  switch (x.n){
    case 2:
      dx_bf_lin_1d (nx.a);
      break;
    case 3:
      dx_bf_quad_1d (nx.a,xi);
      break;
    case 4:
      dksi_bf_cubic_1d (nx.a,xi);
      break;
    default:
      print_err("the case for %ld number of nodes on element has not yet been implemented.",
                __FILE__,__LINE__, __func__, x.n);
  }

  scprd(nx, x, jac);
  cmulv(1.0/jac, nx, dx);  
}



/**
  The function evaluates jacobian for 1D elements given in 1D space.

  @param jac[out] - Jacobian
  @param x[in] - array containing node coordinates
  @param xi[in] - natural coordinate
   
  1.4.2002
*/
void jac_1d (double &jac, vector &x, double xi)
{
  vector nx(ASTCKVEC(x.n));
  
  switch (x.n){
    case 2:
      dx_bf_lin_1d(nx.a);
      break;
    case 3:
      dx_bf_quad_1d(nx.a, xi);
      break;
    case 4:
      dksi_bf_cubic_1d(nx.a, xi);
      break;
    default:
      print_err("the case for %ld number of nodes on element has not yet been implemented.",
                __FILE__,__LINE__, __func__, x.n);
  }
  scprd(nx, x, jac);
}



/**
  The function evaluates derivates of the function with variables xi and eta
  with respect to x and y

  @param dx,dy - arrays containing derivatives with respect to x and y
  @param jac - Jacobian
  @param x,y - array containing node coordinates
  @param xi,eta - natural coordinates
   
  24.6.2001
*/
void derivatives_2d (vector &dx,vector &dy,double &jac,
		     vector &x,vector &y,double xi,double eta)
{
  long i, n=x.n;
  double d1,d2,d3,d4,ddx,ddy;
  vector nx(ASTCKVEC(n)), ny(ASTCKVEC(n));
  
  switch (n){
  case 3:{
    dx_bf_lin_3_2d (nx.a);
    dy_bf_lin_3_2d (ny.a);
    break;
  }
  case 4:{
    dx_bf_lin_4_2d (nx.a, eta);
    dy_bf_lin_4_2d (ny.a, xi);
    break;
  }
  case 6:{
    dx_bf_quad_3_2d (nx.a, xi, eta);
    dy_bf_quad_3_2d (ny.a, xi, eta);
    break;
  }
  case 8:{
    dx_bf_quad_4_2d (nx.a, xi, eta);
    dy_bf_quad_4_2d (ny.a, xi, eta);
    break;
  }
  case 12:{
    dksi_bf_cubic_4_2d (nx.a, xi, eta);
    deta_bf_cubic_4_2d (ny.a, xi, eta);
    break;
  }
  default:{
    print_err("wrong number of nodes on 2D element", __FILE__, __LINE__, __func__);
  }
  }
  
  d1=0.0;  d2=0.0;	d3=0.0;  d4=0.0;
  for (i=0;i<n;i++){
    d1+=nx[i]*x[i];
    d2+=ny[i]*x[i];
    d3+=nx[i]*y[i];
    d4+=ny[i]*y[i];
  }
  
  jac = d1*d4 - d2*d3;
  
  for (i=0;i<dx.n;i++){
    ddx=dx[i];  ddy=dy[i];
    dx[i]=(ddx*d4-ddy*d3)/jac;
    dy[i]=(ddy*d1-ddx*d2)/jac;
  }
}



/**
  The function evaluates Jacobian %matrix composed of derivates of the approximation functions 
  of the global coordinates evaluated at variables xi and eta.

  @param jac       - Jacobian %matrix of derivatives (output)
                     / \sum_{i=1}^{x.n} dN_i/dksi*x_i; \sum_{i=1}^{x.n} dN_i/deta * x_i \
                     \ \sum_{i=1}^{x.n} dN_i/dksi*y_i; \sum_{i=1}^{x.n} dN_i/deta * y_i /
  @param x,y       - vectors containing nodal coordinates (input)
  @param xi,eta    - natural coordinates (input)

  11.2016 by TKo
*/
void jac_2d (matrix &jac, vector &x,vector &y,double xi,double eta)
{
  long n = x.n;
  vector dxi(ASTCKVEC(n)), deta(ASTCKVEC(n));

  switch (n)
  {
    case 3:
      dx_bf_lin_3_2d (dxi.a);
      dy_bf_lin_3_2d (deta.a);
      break;
    case 4:
      dx_bf_lin_4_2d (dxi.a,eta);
      dy_bf_lin_4_2d (deta.a,xi);
      break;
    case 6:
      dx_bf_quad_3_2d (dxi.a,xi,eta);
      dy_bf_quad_3_2d (deta.a,xi,eta);
      break;
    case 8:
      dx_bf_quad_4_2d (dxi.a,xi,eta);
      dy_bf_quad_4_2d (deta.a,xi,eta);
      break;
    case 12:
      dksi_bf_cubic_4_2d (dxi.a,xi,eta);
      deta_bf_cubic_4_2d (deta.a,xi,eta);
      break;
    default:
      print_err("wrong number of nodes (nne=%ld) on 2D element",__FILE__,__LINE__,__func__,x.n);
  }

  if ((jac.m == 2) && (jac.n == 2))
  {
    scprd(dxi, x, jac(0,0));
    scprd(deta, x, jac(0,1));
    scprd(dxi, y, jac(1,0));
    scprd(deta, y, jac(1,1));
  }
  else
  {
    print_err("wrong dimensions of Jacobian matrix (%ldx%ld), it must be 2x2", 
              __FILE__, __LINE__, __func__, jac.m, jac.n);
    abort();
  }

  return;
}



/**
  The function evaluates Jacobian determinant for 2D problems.
   
  @param jac[out] - Jacobian
  @param x,y[in] - array containing node coordinates
  @param xi,eta[in] - natural coordinates
   
  9.7.2001
*/
void jac_2d (double &jac,vector &x,vector &y,double xi,double eta)
{
  long i, n = x.n;
  double d1,d2,d3,d4;
  vector nx(ASTCKVEC(n));
  vector ny(ASTCKVEC(n));
  
  switch (n){
  case 3:{
    dx_bf_lin_3_2d (nx.a);
    dy_bf_lin_3_2d (ny.a);
    break;
  }
  case 4:{
    dx_bf_lin_4_2d (nx.a,eta);
    dy_bf_lin_4_2d (ny.a,xi);
    break;
  }
  case 6:{
    dx_bf_quad_3_2d (nx.a,xi,eta);
    dy_bf_quad_3_2d (ny.a,xi,eta);
    break;
  }
  case 8:{
    dx_bf_quad_4_2d (nx.a,xi,eta);
    dy_bf_quad_4_2d (ny.a,xi,eta);
    break;
  }
  case 12:{
    dksi_bf_cubic_4_2d (nx.a,xi,eta);
    deta_bf_cubic_4_2d (ny.a,xi,eta);
    break;
  }
  default:{
    print_err("wrong number of nodes on 2D element", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  d1=0.0;  d2=0.0;	d3=0.0;  d4=0.0;
  for (i=0;i<n;i++){
    d1+=nx[i]*x[i];
    d2+=ny[i]*x[i];
    d3+=nx[i]*y[i];
    d4+=ny[i]*y[i];
  }
  
  jac = d1*d4 - d2*d3;
}



/**
   functions evaluates second derivates of the function with variables xi and eta
   with respect to x and y
   
   @param bftype - type of approximation functions and coordinate functions
   @param dxdx, dxdy, dydy - arrays containing derivatives with respect to x and y
   @param jac - Jacobian
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   
   JK, 8. 10. 2016, modification 15. 12. 2019
*/
void sec_derivatives_2d (gtypebf bftype,vector &dxdx,vector &dxdy,vector &dydy,
			 vector &x,vector &y,double xi,double eta)
{
  long i;
  double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,jac,det;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double b1,b2,b3;
  double *dxi,*deta,*dxixi,*dxieta,*detaeta;
  double *bfdxi,*bfdeta,*bfdxixi,*bfdxieta,*bfdetaeta;
  double *dx,*dy;

  dxi = new double [x.n];
  deta = new double [x.n];
  dxixi = new double [x.n];
  dxieta = new double [x.n];
  detaeta = new double [x.n];

  bfdxi = new double [dxdx.n];
  bfdeta = new double [dxdx.n];
  bfdxixi = new double [dxdx.n];
  bfdxieta = new double [dxdx.n];
  bfdetaeta = new double [dxdx.n];

  dx = new double [dxdx.n];
  dy = new double [dxdx.n];

  switch (bftype){
  case quadrilateral_bilinear:{
    //  first derivatives of coordinate functions with respect to xi and eta
    dx_bf_lin_4_2d (dxi,eta);
    dy_bf_lin_4_2d (deta,xi);

    //  second derivatives of coordinate functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_lin_4_2d (dxixi);
    dxdy_bf_lin_4_2d (dxieta);
    dydy_bf_lin_4_2d (detaeta);


    //  first derivatives of approximation functions with respect to xi and eta
    dx_bf_lin_4_2d (bfdxi,eta);
    dy_bf_lin_4_2d (bfdeta,xi);
    
    //  second derivatives of approximation functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_lin_4_2d (bfdxixi);
    dxdy_bf_lin_4_2d (bfdxieta);
    dydy_bf_lin_4_2d (bfdetaeta);

    break;
  }
  case quadrilateral_biquadratic:{
    //  first derivatives of coordinate functions with respect to xi and eta
    dx_bf_quad_4_2d (dxi,xi,eta);
    dy_bf_quad_4_2d (deta,xi,eta);

    //  second derivatives of coordinate functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_quad_4_2d (dxixi,eta);
    dxdy_bf_quad_4_2d (dxieta,xi,eta);
    dydy_bf_quad_4_2d (detaeta,xi);

    //  first derivatives of approximation functions with respect to xi and eta
    dx_bf_quad_4_2d (bfdxi,xi,eta);
    dy_bf_quad_4_2d (bfdeta,xi,eta);

    //  second derivatives of approximation functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_quad_4_2d (bfdxixi,eta);
    dxdy_bf_quad_4_2d (bfdxieta,xi,eta);
    dydy_bf_quad_4_2d (bfdetaeta,xi);

    break;
  }
  case batoz_plate:{
    //  first derivatives of coordinate functions with respect to xi and eta
    dx_bf_lin_4_2d (dxi,eta);
    dy_bf_lin_4_2d (deta,xi);

    //  second derivatives of coordinate functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_lin_4_2d (dxixi);
    dxdy_bf_lin_4_2d (dxieta);
    dydy_bf_lin_4_2d (detaeta);


    //  first derivatives of approximation functions with respect to xi and eta
    dx_bf_quad_4_2d (bfdxi,xi,eta);
    dy_bf_quad_4_2d (bfdeta,xi,eta);

    //  second derivatives of approximation functions with respect to xi^2, eta^2 and xi eta
    dxdx_bf_quad_4_2d (bfdxixi,eta);
    dxdy_bf_quad_4_2d (bfdxieta,xi,eta);
    dydy_bf_quad_4_2d (bfdetaeta,xi);
    break;
  }
  default:{
    print_err("wrong number of components of second derivatives on 2D element", __FILE__, __LINE__, __func__);
  }
  }
  
  //  dx/d xi
  d1=0.0;
  //  dx/d eta
  d2=0.0;

  //  dy/d xi
  d3=0.0;
  //  dy/ deta
  d4=0.0;
  
  //  dx^2/d xi^2
  d5=0.0;
  //  dx^2/d xi/d eta
  d6=0.0;
  //  dx^2/d eta^2
  d7=0.0;

  //  dy^2/d xi^2
  d8=0.0;
  //  dy^2/d xi/d eta
  d9=0.0;
  //  dy^2/d eta^2
  d10=0.0;
  

  for (i=0;i<x.n;i++){
    //  dx/d xi
    d1+=dxi[i]*x[i];
    //  dx/d eta
    d2+=deta[i]*x[i];

    //  dy/d xi
    d3+=dxi[i]*y[i];
    //  dy/d eta
    d4+=deta[i]*y[i];
    
    //  dx^2/d xi^2
    d5+=dxixi[i]*x[i];
    //  dx^2/d xi/d eta
    d6+=dxieta[i]*x[i];
    //  dx^2/deta^2
    d7+=detaeta[i]*x[i];
    
    //  dy^2/d xi^2
    d8+=dxixi[i]*y[i];
    //  dy^2/d xi/d eta
    d9+=dxieta[i]*y[i];
    //  dy^2/d eta^2
    d10+=detaeta[i]*y[i];
  }
  
  //  matrix of the system
  a11 = d1*d1;
  a12 = 2.0 * d1*d3;
  a13 = d3*d3;
  
  a21 = d1*d2;
  a22 = d1*d4+d2*d3;
  a23 = d3*d4;
  
  a31 = d2*d2;
  a32 = 2.0*d2*d4;
  a33 = d4*d4;
  
  det = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a23*a32*a11 - a33*a12*a21;
  
  //  computation of the first derivatives of aproximation functions
  jac = d1*d4 - d2*d3;
  
  for (i=0;i<dxdx.n;i++){
    dx[i]=(bfdxi[i]*d4-bfdeta[i]*d3)/jac;
    dy[i]=(bfdeta[i]*d1-bfdxi[i]*d2)/jac;
  }
  
  for (i=0;i<dxdx.n;i++){
    b1 = bfdxixi[i] - dx[i]*d5 - dy[i]*d8;
    b2 = bfdxieta[i] - dx[i]*d6 - dy[i]*d9;
    b3 = bfdetaeta[i] - dx[i]*d7 - dy[i]*d10;
    
    dxdx[i]=(b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a13*a22*b3 - a23*a32*b1 - a33*a12*b2)/det;
    dxdy[i]=(a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a13*b2*a31 - a23*b3*a11 - a33*b1*a21)/det;
    dydy[i]=(a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - b1*a22*a31 - b2*a32*a11 - b3*a12*a21)/det;
  }
  
  delete [] bfdetaeta;
  delete [] bfdxieta;
  delete [] bfdxixi;
  delete [] bfdeta;
  delete [] bfdxi;

  delete [] detaeta;
  delete [] dxieta;
  delete [] dxixi;
  delete [] deta;
  delete [] dxi;

  delete [] dy;
  delete [] dx;
}



/**
   function evaluates derivatives of the function with variables xi, eta and zeta
   with respect to x, y and z
   
   @param dx,dy,dz - arrays containing derivatives with respect to x, y and z
   @param jac - Jacobian
   @param x,y,z - arrays containing node coordinates
   @param xi,eta,zeta - natural coordinates

   24.6.2001
*/
void derivatives_3d (vector &dx,vector &dy,vector &dz,double &jac,
		     vector &x,vector &y,vector &z,
		     double xi,double eta,double zeta)
{
  long i;
  double d1,d2,d3,d4,d5,d6,d7,d8,d9,ddx,ddy,ddz,*nx,*ny,*nz;
  
  nx = new double [x.n];
  ny = new double [x.n];
  nz = new double [x.n];

  switch (x.n){
  case 4:{
    dx_bf_lin_tet (nx);
    dy_bf_lin_tet (ny);
    dz_bf_lin_tet (nz);
    break;
  }
  case 6:{
    dx_bf_lin_wed_3d (nx,zeta);
    dy_bf_lin_wed_3d (ny,zeta);
    dz_bf_lin_wed_3d (nz,xi,eta);
    break;
  }
  case 8:{
    dx_bf_lin_hex_3d (nx,eta,zeta);
    dy_bf_lin_hex_3d (ny,xi,zeta);
    dz_bf_lin_hex_3d (nz,xi,eta);
    break;
  }
  case 10:{
    dx_bf_quad_tet(nx,xi,eta,zeta);
    dy_bf_quad_tet(ny,xi,eta,zeta);
    dz_bf_quad_tet(nz,xi,eta,zeta);
    break;
  }
  case 15:{
    dx_bf_quad_wed_3d (nx,xi,eta,zeta);
    dy_bf_quad_wed_3d (ny,xi,eta,zeta);
    dz_bf_quad_wed_3d (nz,xi,eta,zeta);
    break;
  }
  case 20:{
    dx_bf_quad_hex_3d (nx,xi,eta,zeta);
    dy_bf_quad_hex_3d (ny,xi,eta,zeta);
    dz_bf_quad_hex_3d (nz,xi,eta,zeta);
    break;
  }
  default:{
    fprintf (stderr,"\n\n wrong number of nodes on 3D element in the function");
    fprintf (stderr,"\n derivatives_3d (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  d1=0.0;  d2=0.0;  d3=0.0;  d4=0.0;
  d5=0.0;  d6=0.0;  d7=0.0;  d8=0.0;  d9=0.0;
  for (i=0;i<x.n;i++){
    d1+=nx[i]*x[i];
    d2+=ny[i]*x[i];
    d3+=nz[i]*x[i];
    d4+=nx[i]*y[i];
    d5+=ny[i]*y[i];
    d6+=nz[i]*y[i];
    d7+=nx[i]*z[i];
    d8+=ny[i]*z[i];
    d9+=nz[i]*z[i];
  }
  
  jac = d1*d5*d9 + d4*d8*d3 + d7*d2*d6 - d7*d5*d3 - d8*d6*d1 - d9*d4*d2;
  
  for (i=0;i<dx.n;i++){
    ddx=dx[i];  ddy=dy[i];  ddz=dz[i];
    dx[i]=(ddx*(d5*d9-d8*d6)+ddy*(d7*d6-d9*d4)+ddz*(d4*d8-d7*d5))/jac;
    dy[i]=(ddx*(d8*d3-d9*d2)+ddy*(d1*d9-d7*d3)+ddz*(d7*d2-d8*d1))/jac;
    dz[i]=(ddx*(d2*d6-d5*d3)+ddy*(d4*d3-d6*d1)+ddz*(d1*d5-d4*d2))/jac;
  }
  
  delete [] nx;  delete [] ny;  delete [] nz;
}



/**
   The functions evaluates Jacobian %matrix composed of derivates of the approximation functions 
   of the global coordinates evaluated at variables xi and eta.

   @param jac       - Jacobian %matrix of derivatives (output)
                      / \sum_{i=1}^{x.n} dN_i/dksi*x_i; \sum_{i=1}^{x.n} dN_i/deta * x_i; \sum_{i=1}^{x.n} dN_i/zeta * x_i \
                      | \sum_{i=1}^{x.n} dN_i/dksi*y_i; \sum_{i=1}^{x.n} dN_i/deta * y_i; \sum_{i=1}^{x.n} dN_i/zeta * y_i |
                      \ \sum_{i=1}^{x.n} dN_i/dksi*z_i; \sum_{i=1}^{x.n} dN_i/deta * z_i; \sum_{i=1}^{x.n} dN_i/zeta * z_i /
   @param x,y,z     - vectors containing nodal coordinates (input)
   @param xi,eta    - natural coordinates (input)

   11.2016 by TKo
*/
void jac_3d (matrix &jac, vector &x, vector &y, vector &z, double xi, double eta, double zeta)
{
  vector dxi(ASTCKVEC(x.n)), deta(ASTCKVEC(y.n)), dzeta(ASTCKVEC(y.n));

  switch (x.n)
  {
    case 4:
      dx_bf_lin_tet (dxi.a);
      dy_bf_lin_tet (deta.a);
      dz_bf_lin_tet (dzeta.a);
      break;
    case 6:
      dx_bf_lin_wed_3d (dxi.a, zeta);
      dy_bf_lin_wed_3d (deta.a, zeta);
      dz_bf_lin_wed_3d (dzeta.a, xi, eta);
    break;
    case 8:
      dx_bf_lin_hex_3d (dxi.a, eta, zeta);
      dy_bf_lin_hex_3d (deta.a, xi, zeta);
      dz_bf_lin_hex_3d (dzeta.a, xi, eta);
      break;
    case 10:
      dx_bf_quad_tet(dxi.a, xi, eta, zeta);
      dy_bf_quad_tet(deta.a, xi, eta, zeta);
      dz_bf_quad_tet(dzeta.a, xi, eta, zeta);
      break;
    case 15:
      dx_bf_quad_wed_3d (dxi.a, xi, eta, zeta);
      dy_bf_quad_wed_3d (deta.a, xi, eta, zeta);
      dz_bf_quad_wed_3d (dzeta.a, xi, eta, zeta);
      break;
    case 20:
      dx_bf_quad_hex_3d (dxi.a, xi, eta, zeta);
      dy_bf_quad_hex_3d (deta.a, xi, eta, zeta);
      dz_bf_quad_hex_3d (dzeta.a, xi, eta, zeta);
      break;
    default:
      print_err("wrong number of nodes (nne=%ld) on 3D element",
                __FILE__, __LINE__, __func__, x.n);
  }
  if ((jac.m == 3) && (jac.n == 3))
  {
    scprd(dxi, x, jac(0,0));
    scprd(deta, x, jac(0,1));
    scprd(dzeta, x, jac(0,2));
    scprd(dxi, y, jac(1,0));
    scprd(deta, y, jac(1,1));
    scprd(dzeta, y, jac(1,2));
    scprd(dxi, z, jac(2,0));
    scprd(deta, z, jac(2,1));
    scprd(dzeta, z, jac(2,2));
  }
  else
    print_err("wrong dimension of Jacobian matrix (%ldx%ld), it must be 3x3",
              __FILE__, __LINE__, __func__, jac.m, jac.n);
}



/**
   function evaluates Jacobian determinant for 3D problems
   
   @param jac - Jacobian
   @param x,y,z - arrays containing node coordinates
   @param xi,eta,zeta - natural coordinates

   24.6.2001
*/
void jac_3d (double &jac,vector &x,vector &y,vector &z,
	     double xi,double eta,double zeta)
{
  long i;
  double d1,d2,d3,d4,d5,d6,d7,d8,d9;
  vector nx(ASTCKVEC(x.n)), ny(ASTCKVEC(x.n)), nz(ASTCKVEC(x.n));
  
  switch (x.n){
  case 4:{
    dx_bf_lin_tet (nx.a);
    dy_bf_lin_tet (ny.a);
    dz_bf_lin_tet (nz.a);
    break;
  }
  case 6:{
    dx_bf_lin_wed_3d (nx.a,zeta);
    dy_bf_lin_wed_3d (ny.a,zeta);
    dz_bf_lin_wed_3d (nz.a,xi,eta);
    break;
  }
  case 8:{
    dx_bf_lin_hex_3d (nx.a,eta,zeta);
    dy_bf_lin_hex_3d (ny.a,xi,zeta);
    dz_bf_lin_hex_3d (nz.a,xi,eta);
    break;
  }
  case 10:{
    dx_bf_quad_tet (nx.a,xi,eta,zeta);
    dy_bf_quad_tet (ny.a,xi,eta,zeta);
    dz_bf_quad_tet (nz.a,xi,eta,zeta);
    break;
  }
  case 15:{
    dx_bf_quad_wed_3d (nx.a,xi,eta,zeta);
    dy_bf_quad_wed_3d (ny.a,xi,eta,zeta);
    dz_bf_quad_wed_3d (nz.a,xi,eta,zeta);
    break;
  }
  case 20:{
    dx_bf_quad_hex_3d (nx.a,xi,eta,zeta);
    dy_bf_quad_hex_3d (ny.a,xi,eta,zeta);
    dz_bf_quad_hex_3d (nz.a,xi,eta,zeta);
    break;
  }
  default:{
    print_err(" wrong number of nodes on 3D element.", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  d1=0.0;  d2=0.0;  d3=0.0;  d4=0.0;
  d5=0.0;  d6=0.0;  d7=0.0;  d8=0.0;  d9=0.0;
  for (i=0;i<x.n;i++){
    d1+=nx[i]*x[i];
    d2+=ny[i]*x[i];
    d3+=nz[i]*x[i];
    d4+=nx[i]*y[i];
    d5+=ny[i]*y[i];
    d6+=nz[i]*y[i];
    d7+=nx[i]*z[i];
    d8+=ny[i]*z[i];
    d9+=nz[i]*z[i];
  }
  
  jac = d1*d5*d9 + d4*d8*d3 + d7*d2*d6 - d7*d5*d3 - d8*d6*d1 - d9*d4*d2;

}



/**
  The function computes Jacobian of mapping from 2D curve to <-1;1> segment.
   
  @param jac[out] - the resulting Jacobian
  @param x[in] - nodal x coordinates
  @param y[in] - nodal y coordinates
  @param xi[in] - natural coordinate of the point where the Jacobian should be computed

  1.4.2002
*/
void jac1d2d (double &jac, vector &x, vector &y, double xi)
{
  long i, n=x.n;
  double dx,dy;
  vector nx(ASTCKVEC(n));
  
  switch (n){
  case 2:{
    dx_bf_lin_1d (nx.a);
    break;
  }
  case 3:{
    dx_bf_quad_1d (nx.a,xi);
    break;
  }
  case 4:{
    dksi_bf_cubic_1d (nx.a,xi);
    break;
  }
  default:
    print_err("wrong number of nodes (%ld) on element", __FILE__, __LINE__, __func__, n);
  }
  
  dx=0.0;  dy=0.0;
  for (i=0;i<n;i++){
    dx+=nx[i]*x[i];
    dy+=nx[i]*y[i];
  }
  
  jac = sqrt (dx*dx+dy*dy);
}



/**
  The function computes Jacobian of mapping from 3D curve to <-1;1> segment.
   
  @param jac[out] - resulting Jacobian
  @param x[in] - nodal x coordinates
  @param y[in] - nodal y coordinates
  @param z[in] - nodal z coordinates
  @param xi[in] - natural coordinate of the point where the Jacobian should be computed
   
   1.4.2002
*/
void jac1d3d (double &jac,vector &x,vector &y,vector &z,double xi)
{
  long i, n=x.n;
  double dx, dy, dz;
  vector nx(ASTCKVEC(n));
  
  switch (n){
    case 2:{
      dx_bf_lin_1d (nx.a);
      break;
    }
    case 3:{
      dx_bf_quad_1d (nx.a,xi);
      break;
    }
    case 4:{
      dksi_bf_cubic_1d (nx.a,xi);
      break;
    }
    default:
      print_err("wrong number of nodes (%ld) on element in function", __FILE__, __LINE__, __func__, n);
  }
  
  dx=0.0;  dy=0.0;  dz=0.0;
  for (i=0;i<n;i++){
    dx+=nx[i]*x[i];
    dy+=nx[i]*y[i];
    dz+=nx[i]*z[i];
  }
  
  jac = sqrt (dx*dx+dy*dy+dz*dz);
}



/**
   function computes Jacobian of mapping from surface in 3D to unit square or triangle
   
   @param jac - jacobian
   @param x,y,z - arrays containing node coordinates
   @param xi,eta - natural coordinates
   
   1.4.2002
*/
void jac2d3d (double &jac,vector &x,vector &y,vector &z,double xi,double eta)
{
  long i,n=x.n;
  double dxu,dxv,dyu,dyv,dzu,dzv,g11,g22,g12;
  vector nu(ASTCKVEC(n));  
  vector nv(ASTCKVEC(n));
  
  switch (n){
    case 3:{
      dx_bf_lin_3_2d (nu.a);
      dy_bf_lin_3_2d (nv.a);
      break;
    }
    case 4:{
      dx_bf_lin_4_2d (nu.a,eta);
      dy_bf_lin_4_2d (nv.a,xi);
      break;
    }
    case 6:{
      dx_bf_quad_3_2d (nu.a,xi,eta);
      dy_bf_quad_3_2d (nv.a,xi,eta);
      break;
    }
    case 8:{
      dx_bf_quad_4_2d (nu.a,xi,eta);
      dy_bf_quad_4_2d (nv.a,xi,eta);
      break;
    }
    default:
      print_err(" wrong number of nodes (%ld) on 2D element", __FILE__, __LINE__, __func__, n);
  }
  
  dxu=0.0;  dyu=0.0;  dzu=0.0;  dxv=0.0;  dyv=0.0;  dzv=0.0;
  for (i=0;i<n;i++){
    dxu+=nu[i]*x[i];
    dyu+=nu[i]*y[i];
    dzu+=nu[i]*z[i];
    dxv+=nv[i]*x[i];
    dyv+=nv[i]*y[i];
    dzv+=nv[i]*z[i];
  }
  
  g11=dxu*dxu+dyu*dyu+dzu*dzu;
  g22=dxv*dxv+dyv*dyv+dzv*dzv;
  g12=dxu*dxv+dyu*dyv+dzu*dzv;
  
  jac = sqrt (g11*g22-g12*g12);
}

/**
   function computes jacobian of mapping from 2D curve to line <-1;1>
   curve is edge on 2D finite element
   
   @param jac - jacobian
   @param x,y - arrays containing node coordinates
   @param xi - natural coordinate
   @param edid - edge id
   
   1.4.2002
*/
void jac1d_2d (double &jac,vector &x,vector &y,double xi,long edid)
{
  vector lx, ly;
  
  if (x.n==3){
    reallocv(RSTCKVEC(2, lx));  reallocv(RSTCKVEC(2, ly));
    if (edid==0){
      lx[0]=x[0];  lx[1]=x[1];
      ly[0]=y[0];  ly[1]=y[1];
    }
    if (edid==1){
      lx[0]=x[1];  lx[1]=x[2];
      ly[0]=y[1];  ly[1]=y[2];
    }
    if (edid==2){
      lx[0]=x[2];  lx[1]=x[0];
      ly[0]=y[2];  ly[1]=y[0];
    }
  }

  if (x.n==4){
    reallocv(RSTCKVEC(2, lx));  reallocv(RSTCKVEC(2, ly));
    if (edid==0){
      lx[0]=x[0];  lx[1]=x[1];
      ly[0]=y[0];  ly[1]=y[1];
    }
    if (edid==1){
      lx[0]=x[1];  lx[1]=x[2];
      ly[0]=y[1];  ly[1]=y[2];
    }
    if (edid==2){
      lx[0]=x[2];  lx[1]=x[3];
      ly[0]=y[2];  ly[1]=y[3];
    }
    if (edid==3){
      lx[0]=x[3];  lx[1]=x[0];
      ly[0]=y[3];  ly[1]=y[0];
    }
  }
  
  if (x.n==6){
    reallocv(RSTCKVEC(3, lx));  reallocv(RSTCKVEC(3, ly));
    if (edid==0){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[3];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[3];
    }
    if (edid==1){
      lx[0]=x[1];  lx[1]=x[2];  lx[2]=x[4];
      ly[0]=y[1];  ly[1]=y[2];  ly[2]=y[4];
    }
    if (edid==2){
      lx[0]=x[2];  lx[1]=x[0];  lx[2]=x[5];
      ly[0]=y[2];  ly[1]=y[0];  ly[2]=y[5];
    }
  }

  if (x.n==8){
    reallocv(RSTCKVEC(3, lx));  reallocv(RSTCKVEC(3, ly));
    if (edid==0){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[4];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[4];
    }
    if (edid==1){
      lx[0]=x[1];  lx[1]=x[2];  lx[2]=x[5];
      ly[0]=y[1];  ly[1]=y[2];  ly[2]=y[5];
    }
    if (edid==2){
      lx[0]=x[2];  lx[1]=x[3];  lx[2]=x[6];
      ly[0]=y[2];  ly[1]=y[3];  ly[2]=y[6];
    }
    if (edid==3){
      lx[0]=x[3];  lx[1]=x[0];  lx[2]=x[7];
      ly[0]=y[3];  ly[1]=y[0];  ly[2]=y[7];
    }
  }

  if (x.n==12){
    reallocv(RSTCKVEC(4, lx));  reallocv(RSTCKVEC(4, ly));
    if (edid==0){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[4];  lx[3]=x[5];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[4];  ly[3]=y[5];
    }
    if (edid==1){
      lx[0]=x[1];  lx[1]=x[2];  lx[2]=x[6];  lx[3]=x[7];
      ly[0]=y[1];  ly[1]=y[2];  ly[2]=y[6];  ly[3]=y[7];
    }
    if (edid==2){
      lx[0]=x[2];  lx[1]=x[3];  lx[2]=x[8];  lx[3]=x[9];
      ly[0]=y[2];  ly[1]=y[3];  ly[2]=y[8];  ly[3]=y[9];
    }
    if (edid==3){
      lx[0]=x[3];  lx[1]=x[0];  lx[2]=x[10];  lx[3]=x[11];
      ly[0]=y[3];  ly[1]=y[0];  ly[2]=y[10];  ly[3]=y[11];
    }
  }

  jac1d2d (jac,lx,ly,xi);
}



/**
   function computes jacobian of mapping from surface to unit square or triangle
   surface is surface of 3D finite element
   
   @param jac - jacobian
   @param x,y,z - node coordinates
   @param xi,eta - natural coordinates
   @param sid - surface id
   
*/
void jac2d_3d (double &jac,vector &x,vector &y,vector &z,double xi,double eta,long sid)
{
  vector lx,ly,lz;
  ivector sn;
  long i;
  
  if (x.n==4){
    // linear tetraherdon
    reallocv (RSTCKVEC(3,lx));  reallocv (RSTCKVEC(3,ly));  reallocv (RSTCKVEC(3,lz));
    reallocv (RSTCKIVEC(3,sn));
    if (sid==0){
      /*
      lx[0]=x[1];  lx[1]=x[2];  lx[2]=x[3];
      ly[0]=y[1];  ly[1]=y[2];  ly[2]=y[3];
      lz[0]=z[1];  lz[1]=z[2];  lz[2]=z[3];
      */
      lx[0]=x[2];  lx[1]=x[1];  lx[2]=x[3];
      ly[0]=y[2];  ly[1]=y[1];  ly[2]=y[3];
      lz[0]=z[2];  lz[1]=z[1];  lz[2]=z[3];
    }
    if (sid==1){
      lx[0]=x[0];  lx[1]=x[2];  lx[2]=x[3];
      ly[0]=y[0];  ly[1]=y[2];  ly[2]=y[3];
      lz[0]=z[0];  lz[1]=z[2];  lz[2]=z[3];
    }
    if (sid==2){
      lx[0]=x[1];  lx[1]=x[0];  lx[2]=x[3];
      ly[0]=y[1];  ly[1]=y[0];  ly[2]=y[3];
      lz[0]=z[1];  lz[1]=z[0];  lz[2]=z[3];
    }
    if (sid==3){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[2];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[2];
      lz[0]=z[0];  lz[1]=z[1];  lz[2]=z[2];
    }
  }

  if (x.n==6){
    // linear wedge
    reallocv (RSTCKVEC(3,lx));  reallocv (RSTCKVEC(3,ly));  reallocv (RSTCKVEC(3,lz));
    reallocv (RSTCKIVEC(3,sn));
    if (sid==0){
      reallocv (RSTCKVEC(4,lx));  reallocv (RSTCKVEC(4,ly));  reallocv (RSTCKVEC(4,lz));
      lx[0]=x[2];  lx[1]=x[1];  lx[2]=x[4];  lx[3]=x[5];
      ly[0]=y[2];  ly[1]=y[1];  ly[2]=y[4];  ly[3]=y[5];
      lz[0]=z[2];  lz[1]=z[1];  lz[2]=z[4];  lz[3]=z[5];
    }
    if (sid==1){
      reallocv (RSTCKVEC(4,lx));  reallocv (RSTCKVEC(4,ly));  reallocv (RSTCKVEC(4,lz));
      lx[0]=x[0];  lx[1]=x[2];  lx[2]=x[5];  lx[3]=x[3];
      ly[0]=y[0];  ly[1]=y[2];  ly[2]=y[5];  ly[3]=y[3];
      lz[0]=z[0];  lz[1]=z[2];  lz[2]=z[5];  lz[3]=z[3];
    }
    if (sid==2){
      reallocv (RSTCKVEC(4,lx));  reallocv (RSTCKVEC(4,ly));  reallocv (RSTCKVEC(4,lz));
      lx[0]=x[1];  lx[1]=x[0];  lx[2]=x[3];  lx[3]=x[4];
      ly[0]=y[1];  ly[1]=y[0];  ly[2]=y[3];  ly[3]=y[4];
      lz[0]=z[1];  lz[1]=z[0];  lz[2]=z[3];  lz[3]=z[4];
    }
    if (sid==3){
      reallocv (RSTCKVEC(3,lx));  reallocv (RSTCKVEC(3,ly));  reallocv (RSTCKVEC(3,lz));
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[2];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[2];
      lz[0]=z[0];  lz[1]=z[1];  lz[2]=z[2];
    }
    if (sid==4){
      reallocv (RSTCKVEC(3,lx));  reallocv (RSTCKVEC(3,ly));  reallocv (RSTCKVEC(3,lz));
      lx[0]=x[4];  lx[1]=x[3];  lx[2]=x[5];
      ly[0]=y[4];  ly[1]=y[3];  ly[2]=y[5];
      lz[0]=z[4];  lz[1]=z[3];  lz[2]=z[5];
    }
  }



  if (x.n==8){
    // linear hexahedron
    reallocv (RSTCKVEC(4,lx));  reallocv (RSTCKVEC(4,ly));  reallocv (RSTCKVEC(4,lz));
    reallocv (RSTCKIVEC(4,sn));
    if (sid==0){
      lx[0]=x[0];  lx[1]=x[3];  lx[2]=x[7];  lx[3]=x[4];
      ly[0]=y[0];  ly[1]=y[3];  ly[2]=y[7];  ly[3]=y[4];
      lz[0]=z[0];  lz[1]=z[3];  lz[2]=z[7];  lz[3]=z[4];
    }
    if (sid==1){
      lx[0]=x[1];  lx[1]=x[0];  lx[2]=x[4];  lx[3]=x[5];
      ly[0]=y[1];  ly[1]=y[0];  ly[2]=y[4];  ly[3]=y[5];
      lz[0]=z[1];  lz[1]=z[0];  lz[2]=z[4];  lz[3]=z[5];
    }
    if (sid==2){
      lx[0]=x[2];  lx[1]=x[1];  lx[2]=x[5];  lx[3]=x[6];
      ly[0]=y[2];  ly[1]=y[1];  ly[2]=y[5];  ly[3]=y[6];
      lz[0]=z[2];  lz[1]=z[1];  lz[2]=z[5];  lz[3]=z[6];
    }
    if (sid==3){
      lx[0]=x[3];  lx[1]=x[2];  lx[2]=x[6];  lx[3]=x[7];
      ly[0]=y[3];  ly[1]=y[2];  ly[2]=y[6];  ly[3]=y[7];
      lz[0]=z[3];  lz[1]=z[2];  lz[2]=z[6];  lz[3]=z[7];
    }
    if (sid==4){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[2];  lx[3]=x[3];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[2];  ly[3]=y[3];
      lz[0]=z[0];  lz[1]=z[1];  lz[2]=z[2];  lz[3]=z[3];
    }
    if (sid==5){
      lx[0]=x[4];  lx[1]=x[7];  lx[2]=x[6];  lx[3]=x[5];
      ly[0]=y[4];  ly[1]=y[7];  ly[2]=y[6];  ly[3]=y[5];
      lz[0]=z[4];  lz[1]=z[7];  lz[2]=z[6];  lz[3]=z[5];
    }
  }
  
  if (x.n==10)
  {
    // quadratic tetrahedron
    reallocv (RSTCKVEC(6,lx));  reallocv (RSTCKVEC(6,ly));  reallocv (RSTCKVEC(6,lz));
    reallocv (RSTCKIVEC(6,sn));
    quadtetrahedral_surfnod (sn.a,sid);
    for (i=0; i<6; i++)
    {
      lx[i] = x[sn[i]];
      ly[i] = y[sn[i]];
      lz[i] = z[sn[i]];
    }
  }

  if (x.n==20){
    // quadratic hexahedron
    reallocv (RSTCKVEC(8,lx));  reallocv (RSTCKVEC(8,ly));  reallocv (RSTCKVEC(8,lz));
    reallocv (RSTCKIVEC(8,sn));
    if (sid==0){
      lx[0]=x[0];  lx[1]=x[3];  lx[2]=x[7];  lx[3]=x[4];   lx[4]=x[11];  lx[5]=x[15];  lx[6]=x[19];  lx[7]=x[12];
      ly[0]=y[0];  ly[1]=y[3];  ly[2]=y[7];  ly[3]=y[4];   ly[4]=y[11];  ly[5]=y[15];  ly[6]=y[19];  ly[7]=y[12];
      lz[0]=z[0];  lz[1]=z[3];  lz[2]=z[7];  lz[3]=z[4];   lz[4]=z[11];  lz[5]=z[15];  lz[6]=z[19];  lz[7]=z[12];
    }
    if (sid==1){
      lx[0]=x[1];  lx[1]=x[0];  lx[2]=x[4];  lx[3]=x[5];   lx[4]=x[8];   lx[5]=x[12];  lx[6]=x[16];  lx[7]=x[13];
      ly[0]=y[1];  ly[1]=y[0];  ly[2]=y[4];  ly[3]=y[5];   ly[4]=y[8];   ly[5]=y[12];  ly[6]=y[16];  ly[7]=y[13];
      lz[0]=z[1];  lz[1]=z[0];  lz[2]=z[4];  lz[3]=z[5];   lz[4]=z[8];   lz[5]=z[12];  lz[6]=z[16];  lz[7]=z[13];
    }
    if (sid==2){
      lx[0]=x[2];  lx[1]=x[1];  lx[2]=x[5];  lx[3]=x[6];   lx[4]=x[9];   lx[5]=x[13];  lx[6]=x[17];  lx[7]=x[14];
      ly[0]=y[2];  ly[1]=y[1];  ly[2]=y[5];  ly[3]=y[6];   ly[4]=y[9];   ly[5]=y[13];  ly[6]=y[17];  ly[7]=y[14];
      lz[0]=z[2];  lz[1]=z[1];  lz[2]=z[5];  lz[3]=z[6];   lz[4]=z[9];   lz[5]=z[13];  lz[6]=z[17];  lz[7]=z[14];
    }
    if (sid==3){
      lx[0]=x[3];  lx[1]=x[2];  lx[2]=x[6];  lx[3]=x[7];   lx[4]=x[10];  lx[5]=x[14];  lx[6]=x[18];  lx[7]=x[15];
      ly[0]=y[3];  ly[1]=y[2];  ly[2]=y[6];  ly[3]=y[7];   ly[4]=y[10];  ly[5]=y[14];  ly[6]=y[18];  ly[7]=y[15];
      lz[0]=z[3];  lz[1]=z[2];  lz[2]=z[6];  lz[3]=z[7];   lz[4]=z[10];  lz[5]=z[14];  lz[6]=z[18];  lz[7]=z[15];
    }
    if (sid==4){
      lx[0]=x[0];  lx[1]=x[1];  lx[2]=x[2];  lx[3]=x[3];   lx[4]=x[8];   lx[5]=x[9];   lx[6]=x[10];  lx[7]=x[11];
      ly[0]=y[0];  ly[1]=y[1];  ly[2]=y[2];  ly[3]=y[3];   ly[4]=y[8];   ly[5]=y[9];   ly[6]=y[10];  ly[7]=y[11];
      lz[0]=z[0];  lz[1]=z[1];  lz[2]=z[2];  lz[3]=z[3];   lz[4]=z[8];   lz[5]=z[9];   lz[6]=z[10];  lz[7]=z[11];
    }
    if (sid==5){
      lx[0]=x[4];  lx[1]=x[7];  lx[2]=x[6];  lx[3]=x[5];   lx[4]=x[19];  lx[5]=x[18];  lx[6]=x[17];  lx[7]=x[16];
      ly[0]=y[4];  ly[1]=y[7];  ly[2]=y[6];  ly[3]=y[5];   ly[4]=y[19];  ly[5]=y[18];  ly[6]=y[17];  ly[7]=y[16];
      lz[0]=z[4];  lz[1]=z[7];  lz[2]=z[6];  lz[3]=z[5];   lz[4]=z[19];  lz[5]=z[18];  lz[6]=z[17];  lz[7]=z[16];
    }
  }


  jac2d3d (jac,lx,ly,lz,xi,eta);
}



/**
  The function solves the problem of inverse mapping from global to natural coordinate system of one point given in global coordinates px,py 
  on the 1D element defined by vectors of nodal coordinates x, y. Newton-Raphson method is being used for the soultion of the system of nonlinear 
  equations obtained in the problem and the resulting natural coordinate is stored in arguments xi. The initial estimate of values xi 
  should be given in these arguments at the function call. 

  @param px, py - global coordinates of the given point (input)
  @param x, y - nodal coordinates of the given element   (input)
  @param ni - the maximum number of iteration steps
  @param err - required tolerance, i.e. required Euclidean norm of residual
  @param xi - natural coordinate of the given point (input/output)  
*/
long point_natcoord_1d(/*double px, vector &x, long ni, double err, double &xi*/)
{
  /*  long i;
  double normg;
  matrix jac(ASTCKMAT(1,1)); // Jacobian matrix
  vector nc(ASTCKVEC(1));    // actual value of natural coordinate
  vector dnc(ASTCKVEC(1));   // natural coordinate increment
  vector gc(ASTCKVEC(1));    // global coordinate vector
  vector res(ASTCKVEC(1));   // residual vector
  
  nc(0) = xi;
  gc(0) = px;
  
  bf_1d(res, x, nc(0)); // compute initial estimate of global coordinates
  subv(gc, res, res);                       // calculate residual vector
  normg = normv(res);                       // compute norm of residual
  if (normg < err)
  { // norm of residual is less then tolerance => return natrual coordinates
    xi = nc(0);
    return 0;
  }                                         
  for(i=0; i<ni; i++)
  {
    jac_1d (jac(0,0), x, nc(0));   // assemble Jacobian matrix
    gemp (jac.a, dnc.a, res.a, 1, 1, 1.0e-15, 1); // solve equation system
    addv(nc, dnc);                                // actualize natural coordinate vector
    bf_1d(res, x, nc(0));            // compute global coordinates from the attained natural coordinates 
    subv(gc, res, res);                           // calculate residual vector
    normg = normv(res);                           // calculate norm of residual
    if (normg < err)
    { // residual norm is in the tolarenace => return attained natural coordinates
      xi = nc(0);
      return 0;
    }                                         
  }
  */
  return 1;  // natural coordinates cannot be computed with the required error
}  
 


/**
  The function solves the problem of inverse mapping from global to natural coordinate system of one point given in global coordinates px,py,pz 
  on the 1D element defined by vectors of nodal coordinates x, y and z. Newton-Raphson method is being used for the soultion of the system of nonlinear 
  equations obtained in the problem and the resulting natural coordinates are stored in arguments xi, eta and zeta. The initial estimate of values xi, eta and zeta 
  should be given in these arguments at the function call. 

  @param[in] px, py, pz - global coordinates of the given point
  @param[in] x, y, z - nodal coordinates of the given element
  @param[in] ni - the maximum number of iteration steps
  @param[in] err - required tolerance, i.e. required Euclidean norm of residual
  @param[in/out] xi - natural coordinate of the given point
  @param[out] aerr - attained norm of the residual vector

  @retval  ni+1 - if the tolerance cannot be attained
  @retval [0;ni] - returns the number of iterations performed in the case of successful computation
*/
long point_natcoord_1d_3d(double px, double py, double pz, vector &x, vector &y, vector &z, long ni, double err, double &xi, double &aerr)
{
  long i, j, k;
  vector jac(ASTCKVEC(3));
  double max_jac;
  double nc;    // actual value of natural coordinate
  double dnc;   // natural coordinate increment
  vector gc(ASTCKVEC(3));    // global coordinate vector
  vector res(ASTCKVEC(3));   // residual vector
  
  nc = xi;
  gc(0) = px;
  gc(1) = py;
  gc(2) = pz;
  
  bf_1d_3d(res, x, y, z, nc); // compute initial estimate of global coordinates
  subv(gc, res, res);         // calculate residual vector
  aerr = normv(res);          // compute norm of residual
  if (aerr < err)
  { // norm of residual is less then tolerance => return natrual coordinates
    xi = nc;
    return 0;
  }                                         
  for(i=0; i<ni; i++)
  {
    jac_1d(jac(0), x, nc);
    jac_1d(jac(1), y, nc);
    jac_1d(jac(2), z, nc);
    max_jac = jac(0);
    k = 0;
    for(j=1; j<3; j++){
      if (fabs(jac(j)) > max_jac){
        max_jac = jac(j);
        k = j;
      }
    }      
    dnc = res(k)/jac(k);         // solve equation jac * dnc = |res| = aerr for unknown dnc
    nc += dnc;                   // actualize natural coordinate vector
    bf_1d_3d(res, x, y, z, nc);  // compute global coordinates from the attained natural coordinates 
    subv(gc, res, res);          // calculate residual vector
    aerr = normv(res);           // calculate norm of residual
    if (aerr < err)
    { // residual norm is in the tolarenace => return attained natural coordinates
      xi = nc;
      return i+1;
    }                                         
  }
  xi = nc;
  return ni+1;  // natural coordinates cannot be computed with the required error
}  



/**
  The function solves the problem of inverse mapping from global to natural coordinate system of one point given in global coordinates px,py 
  on the 2D element defined by vectors of nodal coordinates x, y. Newton-Raphson method is being used for the soultion of the system of nonlinear 
  equations obtained in the problem and the resulting natural coordinates are stored in arguments xi and eta. The initial estimate of values xi and eta 
  should be given in these arguments at the function call. 

  @param[in] px, py - global coordinates of the given point
  @param[in] x, y - nodal coordinates of the given element
  @param[in] ni - the maximum number of iteration steps
  @param[in] err - required tolerance, i.e. required Euclidean norm of residual
  @param[in/out] xi, eta - natural coordinates of the given point
  @param[out] aerr - attained norm of residual vector

  @retval  ni+1 - if the tolerance cannot be attained
  @retval [0;ni] - returns the number of iterations performed in the case of successful computation
 
  Created by Tomas Koudelka, 30.11.2016
*/
long point_natcoord_2d(double px, double py, vector &x, vector &y, long ni, double err, double &xi, double &eta, double &aerr)
{
  long i;
  matrix jac(ASTCKMAT(2,2)); // Jacobian matrix
  vector nc(ASTCKVEC(2));    // actual value of natural coordinate
  vector dnc(ASTCKVEC(2));   // natural coordinate increment
  vector gc(ASTCKVEC(2));    // global coordinate vector
  vector res(ASTCKVEC(2));   // residual vector
  double f;                  // normalization factor of the Jacobian matrix and rhs vector
  
  nc(0) = xi;
  nc(1) = eta;
  gc(0) = px;
  gc(1) = py;
  
  bf_2d(res, x, y, nc(0), nc(1)); // compute initial estimate of global coordinates
  subv(gc, res, res);             // calculate residual vector
  aerr = normv(res);             // compute norm of residual
  if (aerr < err)
  { // norm of residual is less then tolerance => return natrual coordinates
    xi = nc(0);
    eta = nc(1);
    return 0;
  }                                         
  for(i=0; i<ni; i++)
  {
    jac_2d (jac, x, y, nc(0), nc(1));             // assemble Jacobian matrix
    f = normalize(jac);                           // normalize system matrix
    cmulv(f, res);                                // normalize rhs vector
    gemp (jac.a, dnc.a, res.a, 2, 1, 1.0e-15, 1); // solve equation system
    addv(nc, dnc, nc);                            // actualize natural coordinate vector
    bf_2d(res, x, y, nc(0), nc(1));               // compute global coordinates from the attained natural coordinates 
    subv(gc, res, res);                           // calculate residual vector
    aerr = normv(res);                           // calculate norm of residual
    if (aerr < err)
    { // residual norm is in the tolarenace => return attained natural coordinates
      xi = nc(0);
      eta = nc(1);
      return i+1;
    }                                         
  }
  xi = nc(0);
  eta = nc(1);
  return ni+1;  // natural coordinates cannot be computed with the required error
}  
 


/**
  The function solves the problem of inverse mapping from global to natural coordinate system of one point given in global coordinates px,py,pz 
  on the 3D element defined by vectors of nodal coordinates x, y and z. Newton-Raphson method is being used for the soultion of the system of nonlinear 
  equations obtained in the problem and the resulting natural coordinates are stored in arguments xi, eta and zeta. The initial estimate of values xi, eta and zeta 
  should be given in these arguments at the function call. 

  @param[in] px, py, pz - global coordinates of the given point
  @param[in] x, y - nodal coordinates of the given element
  @param[in] ni - the maximum number of iteration steps
  @param[in] err - required tolerance, i.e. required Euclidean norm of residual
  @param[in/out] xi, eta, zeta - natural coordinates of the given point
  @param[out] aerr - attained norm of the residual

  @retval ni+1 - if the tolerance cannot be attained
  @retval [0;ni] - returns the number of iterations performed in the case of successful computation
 
  Created by Tomas Koudelka, 30.11.2016
*/
long point_natcoord_3d(double px, double py, double pz, vector &x, vector &y, vector &z, long ni, double err, double &xi, double &eta, double &zeta, double &aerr)
{
  long i;
  matrix jac(ASTCKMAT(3,3)); // Jacobian matrix
  vector nc(ASTCKVEC(3));    // actual value of natural coordinate
  vector dnc(ASTCKVEC(3));   // natural coordinate increment
  vector gc(ASTCKVEC(3));    // global coordinate vector
  vector res(ASTCKVEC(3));   // residual vector
  double f;                  // normalization factor of the Jacobian matrix and rhs vector
  
  nc(0) = xi;
  nc(1) = eta;
  nc(2) = zeta;
  gc(0) = px;
  gc(1) = py;
  gc(2) = pz;
  
  bf_3d(res, x, y, z, nc(0), nc(1), nc(2)); // compute initial estimate of global coordinates
  subv(gc, res, res);                       // calculate residual vector
  aerr = normv(res);                       // compute norm of residual
  if (aerr < err)
  { // norm of residual is less then tolerance => return natrual coordinates
    xi = nc(0);
    eta = nc(1);
    zeta = nc(2);
    return 0;
  }                                         
  for(i=0; i<ni; i++)
  {
    jac_3d (jac, x, y, z, nc(0), nc(1), nc(2));   // assemble Jacobian matrix
    f = normalize(jac);                           // normalize system matrix
    cmulv(f, res);                                // normalize rhs vector
    gemp (jac.a, dnc.a, res.a, 3, 1, 1.0e-15, 1); // solve equation system
    addv(nc, dnc, nc);                            // actualize natural coordinate vector
    bf_3d(res, x, y, z, nc(0), nc(1), nc(2));     // compute global coordinates from the attained natural coordinates 
    subv(gc, res, res);                           // calculate residual vector
    aerr = normv(res);                           // calculate norm of residual
    if (aerr < err)
    { // residual norm is in the tolarenace => return attained natural coordinates
      xi = nc(0);
      eta = nc(1);
      zeta = nc(2);
      return i+1;
    }                                         
  }
  xi = nc(0);
  eta = nc(1);
  zeta = nc(2);
  return ni+1;  // natural coordinates cannot be computed with the required error
}  
