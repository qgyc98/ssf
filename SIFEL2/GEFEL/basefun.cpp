#include "basefun.h"
#include "vector.h"
#include "difcalc.h"
#include <stdlib.h>

/**
   function computes linear base functions on bar element
   natural coordinate is from segment <-1;1>
   
   @param n - array of base functions
   @param x - natural coodinate
*/
void bf_lin_1d (double *n,double x)
{
  n[0]=0.5*(1.0-x);
  n[1]=0.5*(1.0+x);
}

/**
   function computes derivatives of linear base functions on bar element with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
*/
void dx_bf_lin_1d (double *n)
{
  n[0]=-0.5;
  n[1]= 0.5;
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes quadratic approximation functions on bar element
   
   @param n - array containing approximation functions
   @param x - natural coordinate
*/
void bf_quad_1d (double *n,double x)
{
  n[0]=x*(x-1.0)/2.0;
  n[1]=x*(1.0+x)/2.0;
  n[2]=1.0-x*x;
}

/**
   function computes derivatives of quadratic approximation functions 
   with respect of x on bar element
   
   @param n - array containing approximation functions
   @param x - natural coordinate
*/
void dx_bf_quad_1d (double *n,double x)
{
  n[0]=x-0.5;
  n[1]=x+0.5;
  n[2]=-2.0*x;
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes cubic approximation functions on bar element
   
   @param n - array containing approximation functions
   @param ksi - natural coordinate
*/
void bf_cubic_1d(double *n, double ksi)
{
  const double third = 1.0/3.0;
  n[0] =  9.0/16.0*(ksi+third)*(ksi-third)*(1.0-ksi);
  n[1] =  9.0/16.0*(1.0+ksi)  *(ksi+third)*(ksi-third);
  n[2] = 27.0/16.0*(1.0+ksi)  *(ksi-third)*(ksi-1.0);
  n[3] = 27.0/16.0*(1.0+ksi)  *(ksi+third)*(1.0-ksi);
}

/**
   function computes derivatives of quadratic approximation functions 
   with respect of natural coordinate ksi on bar element
   
   @param n - array containing approximation functions
   @param x - natural coordinate
*/
void dksi_bf_cubic_1d(double *n, double ksi)
{
  double ksi2 = ksi*ksi;
  n[0]=1.0/16.0*(-27.0*ksi2 +18.0*ksi +1.0);
  n[1]=1.0/16.0*( 27.0*ksi2 +18.0*ksi -1.0);
  n[2]=27.0/16.0*( 3.0*ksi2 -2.0/3.0*ksi -1.0);
  n[3]=27.0/16.0*(-3.0*ksi2 -2.0/3.0*ksi +1.0);
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
  The function returns approximated global coordinates of the point given by natural coordinates xi 
  for 1D element defined at 1D by nodal coordinates stored in vector x.
  
  @param p - %vector conatining approximated function values (output)
  @param x - vectors of nodal coordinates of the 1D bar element
  @param xi - natural coordinate
*/
void bf_1d(vector &p, vector &x, double xi)
{
  vector n(ASTCKVEC(x.n));
  switch (x.n)
  {
    case 2:
      bf_lin_1d(n.a, xi);
      scprd(n, x, p(0));
      break;
    case 3:
      bf_quad_1d(n.a, xi);
      scprd(n, x, p(0));
      break;
    case 4:
      bf_cubic_1d(n.a, xi);
      scprd(n, x, p(0));
      break;
    default:
      print_err("unknown number of approximation functions (nap=%ld) for bar element", __FILE__, __LINE__, __func__, n.n);
  }
} 

/**
  The function returns approximated global coordinates of the point given by natural coordinates xi 
  for 1D element defined at 2D by nodal coordinates stored in vector x, y.
  
  @param p - %vector conatining approximated function values (output)
  @param x - vectors of nodal coordinates of the 1D bar element
  @param y - vectors of nodal coordinates of the 1D bar element
  @param z - vectors of nodal coordinates of the 1D bar element
  @param xi - natural coordinate

  Created by TKo, 08.2019
*/
void bf_1d_2d(vector &p, vector &x, vector &y, double xi)
{
  vector n(ASTCKVEC(x.n));
  switch (x.n)
  {
    case 2:
      bf_lin_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 3:
      bf_quad_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 4:
      bf_cubic_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    default:
      print_err("unknown number of approximation functions (nap=%ld) for bar element", __FILE__, __LINE__, __func__, n.n);
  }
} 

/**
  The function returns approximated global coordinates of the point given by natural coordinates xi 
  for 1D element defined at 3D by nodal coordinates stored in vector x, y, z.
  
  @param p - %vector conatining approximated function values (output)
  @param x - vectors of nodal coordinates of the 1D bar element
  @param y - vectors of nodal coordinates of the 1D bar element
  @param z - vectors of nodal coordinates of the 1D bar element
  @param xi - natural coordinate

  Created by TKo, 08.2019
*/
void bf_1d_3d(vector &p, vector &x, vector &y, vector &z, double xi)
{
  vector n(ASTCKVEC(x.n));
  switch (x.n)
  {
    case 2:
      bf_lin_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2));
      break;
    case 3:
      bf_quad_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2));
      break;
    case 4:
      bf_cubic_1d(n.a, xi);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2));
      break;
    default:
      print_err("unknown number of approximation functions (nap=%ld) for bar element", __FILE__, __LINE__, __func__, n.n);
  }
} 

/**
   function computes approximation function of dilatation of 2D beam element
   
   @param n - array containing approximation functions
   @param xi - natural coordinate from segment <0;1>
   
   20.2.2002
*/
void dilat (double *n,double xi)
{
  n[0] = 1.0-xi;
  n[1] = xi;
}

/**
   function computes derivative of approximation function of dilatation of 2D beam element with respect of x
   
   @param n - array containing approximation functions
   @param l - length of the beam
   
   20.2.2002
*/
void der_dilat (double *n,double l)
{
  n[0] = -1.0/l;
  n[1] = 1.0/l;
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes approximation function of deflection of 2D beam element
   
   @param n - array containing approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   20.2.2002
*/
void defl2D_fun (double *n,double xi,double l,double kappa)
{
  n[0] = (2.0*xi*xi*xi - 3.0*xi*xi - 2.0*kappa*xi + 2.0*kappa+1.0)/(1.0+2.0*kappa);
  n[1] = (-1.0*xi*xi*xi + (2.0+kappa)*xi*xi - (kappa+1.0)*xi)*l/(1.0+2.0*kappa);
  n[2] = (-2.0*xi*xi*xi + 3.0*xi*xi + 2.0*kappa*xi)/(1.0+2.0*kappa);
  n[3] = (-1.0*xi*xi*xi + (1.0-kappa)*xi*xi + kappa*xi)*l/(1.0+2.0*kappa);
}

/**
   function computes derivative of deflection function with respect of x

   @param n - array containing approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   20.2.2002
*/
void der_defl2D_fun (double *n,double xi,double l,double kappa)
{
  n[0] = (6.0*xi*xi - 6.0*xi - 2.0*kappa)/l/(1.0+2.0*kappa);
  n[1] = (-3.0*xi*xi + 2.0*(2.0+kappa)*xi - (kappa+1.0))/(1.0+2.0*kappa);
  n[2] = (-6.0*xi*xi + 6.0*xi + 2.0*kappa)/l/(1.0+2.0*kappa);
  n[3] = (-3.0*xi*xi + 2.0*(1.0-kappa)*xi + kappa)/(1.0+2.0*kappa);
}

/**
   function computes approximation function of rotation of 2D beam element
   
   @param n - array containing approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   20.2.2002
*/
void roty_fun (double *n,double xi,double l,double kappa)
{
  n[0] = (-6.0*xi*xi + 6.0*xi)/l/(1.0+2.0*kappa);
  n[1] = (3.0*xi*xi - 2.0*(kappa+2.0)*xi + 1.0+2.0*kappa)/(1.0+2.0*kappa);
  n[2] = (6.0*xi*xi - 6.0*xi)/l/(1.0+2.0*kappa);
  n[3] = (3.0*xi*xi - 2.0*(1.0-kappa)*xi)/(1.0+2.0*kappa);
}

/**
   function computes derivative of function of rotation of 2D beam element with respect of x
   
   @param n - array containing approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   20.2.2002
*/
void der_roty_fun (double *n,double xi,double l,double kappa)
{
  n[0] = (-12.0*xi + 6.0)/l/l/(1.0+2.0*kappa);
  n[1] = (6.0*xi - 2.0*(kappa+2.0))/l/(1.0+2.0*kappa);
  n[2] = (12.0*xi - 6.0)/l/l/(1.0+2.0*kappa);
  n[3] = (6.0*xi - 2.0*(1.0-kappa))/l/(1.0+2.0*kappa);
}


// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes a_i coefficients of area-coordinates on triangles
   L(x,y)_i = (a_i + b_i.x + c_i.y)/area/2
   
   @param a - coefficients
   @param x,y - node coordinates
*/
void a_coeff (double *a,double *x,double *y)
{
  a[0]=x[1]*y[2]-x[2]*y[1];
  a[1]=x[2]*y[0]-x[0]*y[2];
  a[2]=x[0]*y[1]-x[1]*y[0];
}



/**
   function computes b_i coefficients of area-coordinates on triangles
   L(x,y)_i = (a_i + b_i.x + c_i.y)/area/2
   
   @param b - coefficients
   @param y - node coordinates
*/
void b_coeff (double *b,double *y)
{
  b[0]=y[1]-y[2];
  b[1]=y[2]-y[0];
  b[2]=y[0]-y[1];
}



/**
   function computes c_i coefficients of area-coordinates on triangles
   L(x,y)_i = (a_i + b_i.x + c_i.y)/area/2
   
   @param c - coefficients
   @param x - node coordinates
*/
void c_coeff (double *c,double *x)
{
  c[0]=x[2]-x[1];
  c[1]=x[0]-x[2];
  c[2]=x[1]-x[0];
}



/**
   function computes a_i coefficients of area-coordinates on triangles
   L(x,y)_i = a_i + b_i.x + c_i.y
   
   @param a - coefficients
   @param x,y - node coordinates
   @param det - determinant=2*area
*/
void plsa (double *a,double *x,double *y,double det)
{
  a[0]=(x[1]*y[2]-x[2]*y[1])/det;
  a[1]=(x[2]*y[0]-x[0]*y[2])/det;
  a[2]=(x[0]*y[1]-x[1]*y[0])/det;
}



/**
   function computes b_i coefficients of area-coordinates on triangles
   L(x,y)_i = a_i + b_i.x + c_i.y
   
   @param b - coefficients
   @param y - node coordinates
   @param det - determinant=2*area
*/
void plsb (double *b,double *y,double det)
{
  b[0]=(y[1]-y[2])/det;
  b[1]=(y[2]-y[0])/det;
  b[2]=(y[0]-y[1])/det;
}



/**
   function computes c_i coefficients of area-coordinates on triangles
   L(x,y)_i = a_i + b_i.x + c_i.y
   
   @param c - coefficients
   @param x - node coordinates
   @param det - determinant=2*area
*/
void plsc (double *c,double *x,double det)
{
  c[0]=(x[2]-x[1])/det;
  c[1]=(x[0]-x[2])/det;
  c[2]=(x[1]-x[0])/det;
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************



/**
   function computes quadratic base functions on triangle
   
   @param n - array of base functions
   @param l - area coordinates
*/
void ac_bf_quad_3_2d (double *n,double *l)
{
  n[0]=l[0]*(2.0*l[0]-1.0);
  n[1]=l[1]*(2.0*l[1]-1.0);
  n[2]=l[2]*(2.0*l[2]-1.0);

  n[3]=l[0]*l[1]*4.0;
  n[4]=l[1]*l[2]*4.0;
  n[5]=l[2]*l[0]*4.0;

  //n[3]=l[1]*l[2]*4.0;
  //n[4]=l[2]*l[0]*4.0;
  //n[5]=l[0]*l[1]*4.0;
}



/**
   function computes derivatives of quadratic base functions on triangle with respect of x
   
   @param n - array of derivatives of base functions
   @param b - array of coefficients of area coordinates
   @param l - area coordinates
*/
void ac_dx_bf_quad_3_2d (double *n,double *b,double *l)
{
  n[0]=b[0]*(2.0*l[0]-1.0)+l[0]*2.0*b[0];
  n[1]=b[1]*(2.0*l[1]-1.0)+l[1]*2.0*b[1];
  n[2]=b[2]*(2.0*l[2]-1.0)+l[2]*2.0*b[2];

  n[3]=(b[0]*l[1]+l[0]*b[1])*4.0;
  n[4]=(b[1]*l[2]+l[1]*b[2])*4.0;
  n[5]=(b[2]*l[0]+l[2]*b[0])*4.0;

  //n[3]=(b[1]*l[2]+l[1]*b[2])*4.0;
  //n[4]=(b[2]*l[0]+l[2]*b[0])*4.0;
  //n[5]=(b[0]*l[1]+l[0]*b[1])*4.0;
}



/**
   function computes derivatives of quadratic base functions on triangle with respect of y
   
   @param n - array of derivatives of base functions
   @param c - array of coefficients of area coordinates
   @param l - area coordinates
*/
void ac_dy_bf_quad_3_2d (double *n,double *c,double *l)
{
  n[0]=c[0]*(2.0*l[0]-1.0)+l[0]*2.0*c[0];
  n[1]=c[1]*(2.0*l[1]-1.0)+l[1]*2.0*c[1];
  n[2]=c[2]*(2.0*l[2]-1.0)+l[2]*2.0*c[2];

  n[3]=(c[0]*l[1]+l[0]*c[1])*4.0;
  n[4]=(c[1]*l[2]+l[1]*c[2])*4.0;
  n[5]=(c[2]*l[0]+l[2]*c[0])*4.0;

  //n[3]=(c[1]*l[2]+l[1]*c[2])*4.0;
  //n[4]=(c[2]*l[0]+l[2]*c[0])*4.0;
  //n[5]=(c[0]*l[1]+l[0]*c[1])*4.0;
}


/**
   function computes linear approximation functions on triangle
   
   @param n - array containing approximation functions
   @param x,y - natural coordinates
   
   1.4.2002
*/
void bf_lin_3_2d (double *n,double x,double y)
{
  n[0]=x;
  n[1]=y;
  n[2]=1.0-x-y;
}

/**
   function computes derivatives of linear approximation functions with respect to natural coordinate xi
   
   @param n - array containing derivatives
   
   1.4.2002
*/
void dx_bf_lin_3_2d (double *n)
{
  n[0] =  1.0;
  n[1] =  0.0;
  n[2] = -1.0;
}

/**
   function computes derivatives of linear approximation functions with respect to natural coordinate eta
   
   @param n - array containing derivatives
   
   1.4.2002
*/
void dy_bf_lin_3_2d (double *n)
{
  n[0] =  0.0;
  n[1] =  1.0;
  n[2] = -1.0;
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes quadratic approximation functions on triangle
   
   @param n - array containing approximation functions
   @param x,y - natural coordinates
   
   1.4.2002
*/
void bf_quad_3_2d (double *n,double x,double y)
{
  n[0]=x*(x-0.5)*2.0;
  n[1]=y*(y-0.5)*2.0;
  n[2]=(1.0-x-y)*(0.5-x-y)*2.0;
  n[3]=x*y*4.0;
  n[4]=y*(1.0-x-y)*4.0;
  n[5]=x*(1.0-x-y)*4.0;
}

/**
   function computes derivatives of quadratic approximation functions on triangle with respect to xi
   
   @param n - array containing approximation functions
   @param x,y - natural coordinates
   
   1.4.2002
*/
void dx_bf_quad_3_2d (double *n,double x,double y)
{
  n[0]=4.0*x-1.0;
  n[1]=0.0;
  n[2]=4.0*x+4.0*y-3.0;
  n[3]=4.0*y;
  n[4]=-4.0*y;
  n[5]=4.0-8.0*x-4.0*y;
}

/**
   function computes derivatives of quadratic approximation functions on triangle with respect to eta
   
   @param n - array containing approximation functions
   @param x,y - natural coordinates
   
   1.4.2002
*/
void dy_bf_quad_3_2d (double *n,double x,double y)
{
  n[0]=0.0;
  n[1]=4.0*y-1.0;
  n[2]=4.0*x+4.0*y-3.0;
  n[3]=4.0*x;
  n[4]=4.0-4.0*x-8.0*y;
  n[5]=-4.0*x;
}


// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes bi-linear base functions on quadrilateral
   
   @param n - array of base functions
   @param x,y - natural coodinates
*/
void bf_lin_4_2d (double *n,double x,double y)
{
  n[0]=0.25*(1.0+x)*(1.0+y);
  n[1]=0.25*(1.0-x)*(1.0+y);
  n[2]=0.25*(1.0-x)*(1.0-y);
  n[3]=0.25*(1.0+x)*(1.0-y);
}



/**
   function computes first derivatives of bi-linear base functions
   on quadrilateral with respect to natural coordinate x
   
   @param n - array of derivatives of base functions
   @param y - natural coordinate
*/
void dx_bf_lin_4_2d (double *n,double y)
{
  n[0]=0.25*(1.0+y);
  n[1]=0.25*(-1.0-y);
  n[2]=0.25*(-1.0+y);
  n[3]=0.25*(1.0-y);
}



/**
   function computes first derivatives of bi-linear base functions
   on quadrilateral with respect to natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x - natural coordinate
*/
void dy_bf_lin_4_2d (double *n,double x)
{
  n[0]=0.25*(1.0+x);
  n[1]=0.25*(1.0-x);
  n[2]=0.25*(-1.0+x);
  n[3]=0.25*(-1.0-x);
}

/**
   function computes second derivatives of bi-linear base functions
   on quadrilateral with respect to natural coordinate x
   
   @param n - array of derivatives of base functions
*/
void dxdx_bf_lin_4_2d (double *n)
{
  n[0]=0.0;
  n[1]=0.0;
  n[2]=0.0;
  n[3]=0.0;
}


/**
   function computes second derivatives of bi-linear base functions
   on quadrilateral with respect to natural coordinates x and y
   
   @param n - array of derivatives of base functions
*/
void dxdy_bf_lin_4_2d (double *n)
{
  n[0] = 0.25;
  n[1] = -0.25;
  n[2] = 0.25;
  n[3] = -0.25;
}

/**
   function computes second derivatives of bi-linear base functions
   on quadrilateral with respect to natural coordinate y
   
   @param n - array of derivatives of base functions
*/
void dydy_bf_lin_4_2d (double *n)
{
  n[0]=0.0;
  n[1]=0.0;
  n[2]=0.0;
  n[3]=0.0;
}


// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************



/**
   function computes bi-quadratic base functions on quadrilateral
   
   @param n - array of base functions
   @param x,y - natural coodinates
*/
void bf_quad_4_2d (double *n,double x,double y)
{
  n[0]=0.25*(1.0+x)*(1.0+y)*(x+y-1.0);
  n[1]=0.25*(1.0-x)*(1.0+y)*(y-x-1.0);
  n[2]=0.25*(1.0-x)*(1.0-y)*(-x-y-1.0);
  n[3]=0.25*(1.0+x)*(1.0-y)*(x-y-1.0);
  n[4]=0.5*(1.0-x*x)*(1.0+y);
  n[5]=0.5*(1.0-x)*(1.0-y*y);
  n[6]=0.5*(1.0-x*x)*(1.0-y);
  n[7]=0.5*(1.0+x)*(1.0-y*y);
}



/**
   function computes derivatives of bi-quadratic base functions on quadrilateral with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
   @param x,y - natural coodinates
*/
void dx_bf_quad_4_2d (double *n,double x,double y)
{
  n[0]=0.25*(1.0+y)*(x+y-1.0)+0.25*(1.0+x)*(1.0+y);
  n[1]=0.25*(-1.0-y)*(y-x-1.0)-0.25*(1.0-x)*(1.0+y);
  n[2]=0.25*(-1.0+y)*(-x-y-1.0)-0.25*(1.0-x)*(1.0-y);
  n[3]=0.25*(1.0-y)*(x-y-1.0)+0.25*(1.0+x)*(1.0-y);
  n[4]=(1.0+y)*(-1.0)*x;
  n[5]=0.5*(y*y-1.0);
  n[6]=(y-1.0)*x;
  n[7]=0.5*(1.0-y*y);
}



/**
   function computes derivatives of bi-quadratic base functions on quadrilateral with respect of natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x,y - natural coodinates
*/
void dy_bf_quad_4_2d (double *n,double x,double y)
{
  n[0]=0.25*(1.0+x)*(x+y-1.0)+0.25*(1.0+x)*(1.0+y);
  n[1]=0.25*(1.0-x)*(y-x-1.0)+0.25*(1.0-x)*(1.0+y);
  n[2]=0.25*(1.0-x)*(x+y+1.0)-0.25*(1.0-x)*(1.0-y);
  n[3]=0.25*(1.0+x)*(y-x+1.0)-0.25*(1.0+x)*(1.0-y);
  n[4]=0.5*(1.0-x*x);
  n[5]=(x-1.0)*y;
  n[6]=0.5*(x*x-1.0);
  n[7]=(1.0+x)*(-1.0)*y;
}


/**
   function computes second derivatives of bi-quadratic base functions on quadrilateral with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
   @param y - natural coodinate
   
   JK, 14. 12. 2019
*/
void dxdx_bf_quad_4_2d (double *n,double y)
{
  n[0]=0.5*(1.0+y);
  n[1]=0.5*(1.0+y);
  n[2]=0.5*(1.0-y);
  n[3]=0.5*(1.0-y);
  n[4]=-1.0-y;
  n[5]=0.0;
  n[6]=-1.0+y;
  n[7]=0.0;
}


/**
   function computes second derivatives of bi-quadratic base functions on quadrilateral with respect of natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x - natural coodinate
   
   JK, 14. 12. 2019
*/
void dydy_bf_quad_4_2d (double *n,double x)
{
  n[0]=0.5*(1.0+x);
  n[1]=0.5*(1.0-x);
  n[2]=0.5*(1.0-x);
  n[3]=0.5*(1.0+x);
  n[4]=0.0;
  n[5]=-1.0+x;
  n[6]=0.0;
  n[7]=-1.0-x;
}


/**
   function computes second mixed derivatives of bi-quadratic base functions on quadrilateral
   
   @param n - array of derivatives of base functions
   @param x,y - natural coodinates
   
   JK, 14. 12. 2019
*/
void dxdy_bf_quad_4_2d (double *n,double x,double y)
{
  n[0]=0.5*(x+y)+0.25;
  n[1]=0.5*(x-y)-0.25;
  n[2]=-0.5*(x+y)+0.25;
  n[3]=-0.5*(x-y)-0.25;
  n[4]=0.0-x;
  n[5]=y;
  n[6]=x;
  n[7]=0.0-y;
}





/**
  Function evaluates bicubic base function on quadrilateral elements at the given point 
  on quarilateral element.

  @param n - array of base function values
  @param ksi - natural coordinate in the horizontal direction
  @param eta - natural coordinate in the vertical direction

  Created by Tomas Koudelka, 13.3.2019
*/
void bf_cubic_4_2d (double *n, double ksi, double eta)
{
  double reci32 = 1.0/32.0*(9.0*(ksi*ksi+eta*eta)-10.0);
  //
  // corner nodes
  // 
  n[0]=reci32*(1.0+ksi)*(1.0+eta);
  n[1]=reci32*(1.0-ksi)*(1.0+eta);
  n[2]=reci32*(1.0-ksi)*(1.0-eta);
  n[3]=reci32*(1.0+ksi)*(1.0-eta);
  //
  // mid-side nodes
  //
  // edge 1 between nodes 1 and 2, ksi={1/3;-1/3}, eta=1
  reci32 = 9.0/32.0*(1-ksi*ksi);
  n[4] = reci32*(1+3.0*ksi)*(1+eta);
  n[5] = reci32*(1-3.0*ksi)*(1+eta);
  // edge 3 between nodes 3 and 4, ksi={-1/3;+1/3}, eta=-1
  n[8] = reci32*(1-3.0*ksi)*(1-eta);
  n[9] = reci32*(1+3.0*ksi)*(1-eta);
  // edge 2 between nodes 2 and 3, ksi=-1, eta={1/3;-1/3}
  reci32 = 9.0/32.0*(1-eta*eta);
  n[6] = reci32*(1-ksi)*(1+3.0*eta);
  n[7] = reci32*(1-ksi)*(1-3.0*eta);
  // edge 4 between nodes 4 and 1, ksi=1, eta={-1/3;1/3}
  n[10] = reci32*(1+ksi)*(1-3.0*eta);
  n[11] = reci32*(1+ksi)*(1+3.0*eta);
}



/**
  The function computes derivatives of bi-cubic base functions on quadrilateral 
  with respect of natural coordinate ksi.
   
  @param n[out] - array of derivatives of base functions
  @param x[in], y[in] - natural coodinates

  Created by Tomas Koudelka, 13.3.2019
*/
void dksi_bf_cubic_4_2d(double *n, double ksi, double eta)
{
  const double c32 = 1.0/32.0;
  const double c932 = 9.0/32.0;
  double ksi2 = ksi*ksi;
  double eta2 = eta*eta;

  n[0] =   c32*(1.0+eta)*( 27.0*ksi2 +9.0*eta2 +18.0*ksi -10.0);
  n[1] =   c32*(1.0+eta)*(-27.0*ksi2 -9.0*eta2 +18.0*ksi +10.0);
  n[2] =   c32*(1.0-eta)*(-27.0*ksi2 -9.0*eta2 +18.0*ksi +10.0);
  n[3] =   c32*(1.0-eta)*( 27.0*ksi2 +9.0*eta2 +18.0*ksi -10.0);

  n[4] = c932*(1.0+eta)*(-9.0*ksi2 -2.0*ksi +3.0);
  n[5] = c932*(1.0+eta)*( 9.0*ksi2 -2.0*ksi -3.0);

  n[6] = c932*(1.0+3.0*eta)*(eta2 -1.0);
  n[7] = c932*(1.0-3.0*eta)*(eta2 -1.0);

  n[8] = c932*(1.0-eta)*( 9.0*ksi2 -2.0*ksi -3.0);
  n[9] = c932*(1.0-eta)*(-9.0*ksi2 -2.0*ksi +3.0);

  n[10] = c932*(1.0-3.0*eta)*(1.0 -eta2);
  n[11] = c932*(1.0+3.0*eta)*(1.0 -eta2);
}



/**
  The function computes derivatives of bi-cubic base functions on quadrilateral with 
  respect of natural coordinate eta
   
  @param n[out] - array of derivatives of base functions
  @param ksi[in],eta[in] - natural coodinates
*/
void deta_bf_cubic_4_2d (double *n,double ksi, double eta)
{
  const double c32 = 1.0/32.0;
  const double c932 = 9.0/32.0;
  double ksi2 = ksi*ksi;
  double eta2 = eta*eta;

  n[0] =   c32*(1.0+ksi)*( 27.0*eta2 +9.0*ksi2 +18.0*eta -10.0);
  n[1] =   c32*(1.0-ksi)*( 27.0*eta2 +9.0*ksi2 +18.0*eta -10.0);
  n[2] =   c32*(1.0-ksi)*(-27.0*eta2 -9.0*ksi2 +18.0*eta +10.0);
  n[3] =   c32*(1.0+ksi)*(-27.0*eta2 -9.0*ksi2 +18.0*eta +10.0);

  n[4] = c932*(1.0 -ksi2)*(1.0+3.0*ksi);
  n[5] = c932*(1.0 -ksi2)*(1.0-3.0*ksi);

  n[6] = c932*(1.0-ksi)*(-9.0*eta2 -2.0*eta +3.0);
  n[7] = c932*(1.0-ksi)*( 9.0*eta2 -2.0*eta -3.0);

  n[8] = c932*(ksi2 -1.0)*(1.0-3.0*ksi);
  n[9] = c932*(ksi2 -1.0)*(1.0+3.0*ksi);

  n[10] = c932*(1.0+ksi)*( 9.0*eta2 -2.0*eta -3.0);
  n[11] = c932*(1.0+ksi)*(-9.0*eta2 -2.0*eta +3.0);
}



/**
   function computes modified base functions on triangle for element with rotational DOFs
   
   @param n - array of base functions
   @param l - area coordinates
   @param b,c - coefficients of area coordinates (not divided by 2*area, i.e. y-y,x-x only)
*/
void bf_rot_3_2d (double *n,double *l,double *b,double *c)
{
  n[0] = l[0];
  n[1] = l[1];
  n[2] = l[2];

  n[3] = (l[0]*l[1]*b[2] - l[2]*l[0]*b[1])/2.0;
  n[4] = (l[1]*l[2]*b[0] - l[0]*l[1]*b[2])/2.0;
  n[5] = (l[2]*l[0]*b[1] - l[1]*l[2]*b[0])/2.0;

  n[6] = (l[0]*l[1]*c[2] - l[2]*l[0]*c[1])/2.0;
  n[7] = (l[1]*l[2]*c[0] - l[0]*l[1]*c[2])/2.0;
  n[8] = (l[2]*l[0]*c[1] - l[1]*l[2]*c[0])/2.0;
}



/**
   function computes derivatives of modified base functions on triangle for element with rotational DOFs with respect of x
   
   @param n - array of base functions
   @param l - area coordinates
   @param b,c - coefficients of area coordinates (not divided by 2*area, i.e. y-y,x-x only)
   @param area - area
*/
void dx_bf_rot_3_2d (double *n,double *l,double *b,double *c,double area)
{
  //  auxiliary coefficient for speedup
  double coef;
  coef = 0.5/area;
  //  coef is 1/2.0/area now

  n[0] = b[0]*coef;
  n[1] = b[1]*coef;
  n[2] = b[2]*coef;
  
  coef*=0.5;
  //  coef is 1/4.0/area now
  
  n[3] = (b[2]*l[1] - b[1]*l[2])*b[0]*coef;
  n[4] = (b[0]*l[2] - b[2]*l[0])*b[1]*coef;
  n[5] = (b[1]*l[0] - b[0]*l[1])*b[2]*coef;

  n[6] = (l[0]*(b[1]*c[2]-b[2]*c[1]) + b[0]*(c[2]*l[1]-c[1]*l[2]))*coef;
  n[7] = (l[1]*(b[2]*c[0]-b[0]*c[2]) + b[1]*(c[0]*l[2]-c[2]*l[0]))*coef;
  n[8] = (l[2]*(b[0]*c[1]-b[1]*c[0]) + b[2]*(c[1]*l[0]-c[0]*l[1]))*coef;
}



/**
   function computes derivatives of modified base functions on triangle for element with rotational DOFs with respect of y
   
   @param n - array of base functions
   @param l - area coordinates
   @param b,c - coefficients of area coordinates (not divided by 2*area, i.e. y-y,x-x only)
   @param area - area
*/
void dy_bf_rot_3_2d (double *n,double *l,double *b,double *c,double area)
{
  //  auxiliary coefficient for speedup
  double coef;
  coef = 0.5/area;
  //  coef is 1/2.0/area now

  n[0] = c[0]*coef;
  n[1] = c[1]*coef;
  n[2] = c[2]*coef;

  coef*=0.5;
  //  coef is 1/4.0/area now
  
  n[3] = (l[0]*(c[1]*b[2]-c[2]*b[1]) + c[0]*(b[2]*l[1]-b[1]*l[2]))*coef;
  n[4] = (l[1]*(c[2]*b[0]-c[0]*b[2]) + c[1]*(b[0]*l[2]-b[2]*l[0]))*coef;
  n[5] = (l[2]*(c[0]*b[1]-c[1]*b[0]) + c[2]*(b[1]*l[0]-b[0]*l[1]))*coef;
  
  n[6] = (c[2]*l[1] - c[1]*l[2])*c[0]*coef;
  n[7] = (c[0]*l[2] - c[2]*l[0])*c[1]*coef;
  n[8] = (c[1]*l[0] - c[0]*l[1])*c[2]*coef;
}



/**
   function computes modified base functions on quadrilateral for element with rotational DOFs
   
   @param n - array of base functions
   @param x,y - natural coordinates
   @param nx,ny - array of x and y components of outer normal vectors
   @param l - array of lengths of element edges
*/
void bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l)
{
  n[0]  = (1.0+x)*(1.0+y)/4.0;
  n[1]  = (1.0-x)*(1.0+y)/4.0;
  n[2]  = (1.0-x)*(1.0-y)/4.0;
  n[3]  = (1.0+x)*(1.0-y)/4.0;

  n[4]  = ((1.0-y*y)*(1.0+x)*l[3]*nx[3]-(1.0-x*x)*(1.0+y)*l[0]*nx[0])/16.0;
  n[5]  = ((1.0-x*x)*(1.0+y)*l[0]*nx[0]-(1.0-y*y)*(1.0-x)*l[1]*nx[1])/16.0;
  n[6]  = ((1.0-y*y)*(1.0-x)*l[1]*nx[1]-(1.0-x*x)*(1.0-y)*l[2]*nx[2])/16.0;
  n[7]  = ((1.0-x*x)*(1.0-y)*l[2]*nx[2]-(1.0-y*y)*(1.0+x)*l[3]*nx[3])/16.0;

  n[8]  = ((1.0-y*y)*(1.0+x)*l[3]*ny[3]-(1.0-x*x)*(1.0+y)*l[0]*ny[0])/16.0;
  n[9]  = ((1.0-x*x)*(1.0+y)*l[0]*ny[0]-(1.0-y*y)*(1.0-x)*l[1]*ny[1])/16.0;
  n[10] = ((1.0-y*y)*(1.0-x)*l[1]*ny[1]-(1.0-x*x)*(1.0-y)*l[2]*ny[2])/16.0;
  n[11] = ((1.0-x*x)*(1.0-y)*l[2]*ny[2]-(1.0-y*y)*(1.0+x)*l[3]*ny[3])/16.0;
}



/**
   function computes derivatives of modified base functions on quadrilateral for
   element with rotational DOFs with respect of natural coordinate x
   
   @param n - array of base functions
   @param x,y - natural coordinates
   @param nx,ny - array of x and y components of outer normal vectors
   @param l - array of lengths of element edges
*/
void dx_bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l)
{
  n[0]  = (1.0+y)/4.0;
  n[1]  = (1.0+y)/(-4.0);
  n[2]  = (1.0-y)/(-4.0);
  n[3]  = (1.0-y)/4.0;

  n[4]  = (     (1.0-y*y)*l[3]*nx[3] + 2.0*x*(1.0+y)*l[0]*nx[0])/16.0;
  n[5]  = (-2.0*x*(1.0+y)*l[0]*nx[0] +     (1.0-y*y)*l[1]*nx[1])/16.0;
  n[6]  = (-1.0*(1.0-y*y)*l[1]*nx[1] + 2.0*x*(1.0-y)*l[2]*nx[2])/16.0;
  n[7]  = (-2.0*x*(1.0-y)*l[2]*nx[2] -     (1.0-y*y)*l[3]*nx[3])/16.0;

  n[8]  = (     (1.0-y*y)*l[3]*ny[3] + 2.0*x*(1.0+y)*l[0]*ny[0])/16.0;
  n[9]  = (-2.0*x*(1.0+y)*l[0]*ny[0] +     (1.0-y*y)*l[1]*ny[1])/16.0;
  n[10] = (-1.0*(1.0-y*y)*l[1]*ny[1] + 2.0*x*(1.0-y)*l[2]*ny[2])/16.0;
  n[11] = (-2.0*x*(1.0-y)*l[2]*ny[2] -     (1.0-y*y)*l[3]*ny[3])/16.0;
}



/**
   function computes derivatives of modified base functions on quadrilateral for
   element with rotational DOFs with respect of natural coordinate y
   
   @param n - array of base functions
   @param x,y - natural coordinates
   @param nx,ny - array of x and y components of outer normal vectors
   @param l - array of lengths of element edges
*/
void dy_bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l)
{
  n[0]  = (1.0+x)/4.0;
  n[1]  = (1.0-x)/4.0;
  n[2]  = (1.0-x)/(-4.0);
  n[3]  = (1.0+x)/(-4.0);

  n[4]  = (-2.0*y*(1.0+x)*l[3]*nx[3] -     (1.0-x*x)*l[0]*nx[0])/16.0;
  n[5]  = (     (1.0-x*x)*l[0]*nx[0] + 2.0*y*(1.0-x)*l[1]*nx[1])/16.0;
  n[6]  = (-2.0*y*(1.0-x)*l[1]*nx[1] +     (1.0-x*x)*l[2]*nx[2])/16.0;
  n[7]  = (-1.0*(1.0-x*x)*l[2]*nx[2] + 2.0*y*(1.0+x)*l[3]*nx[3])/16.0;

  n[8]  = (-2.0*y*(1.0+x)*l[3]*ny[3] -     (1.0-x*x)*l[0]*ny[0])/16.0;
  n[9]  = (     (1.0-x*x)*l[0]*ny[0] + 2.0*y*(1.0-x)*l[1]*ny[1])/16.0;
  n[10] = (-2.0*y*(1.0-x)*l[1]*ny[1] +     (1.0-x*x)*l[2]*ny[2])/16.0;
  n[11] = (-1.0*(1.0-x*x)*l[2]*ny[2] + 2.0*y*(1.0+x)*l[3]*ny[3])/16.0;
}



/**
  The function returns approximated global coordinates of the point given by natural coordinates (xi, eta) 
  for 2D elements given by nodal coordinates stored in vectors x, y.
  
  @param p - %vector conatining approximated function values (output)
  @param x, y - vectors of nodal coordinates of the 2D plane element
  @param xi, eta - natural coordinates
*/
void bf_2d(vector &p, vector &x, vector &y, double xi, double eta)
{
  vector n(ASTCKVEC(x.n));
  switch (x.n)
  {
    case 3:
      bf_lin_3_2d (n.a, xi, eta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 4:
      bf_lin_4_2d (n.a, xi, eta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 6:
      bf_quad_3_2d (n.a, xi, eta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 8:
      bf_quad_4_2d (n.a, xi, eta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    case 12:
      bf_cubic_4_2d (n.a, xi, eta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      break;
    default:
      print_err("unknown number of approximation functions (nap=%ld) for bar element", __FILE__, __LINE__, __func__, n.n);
  }
} 



/**
   function computes base functions of constant curvature triangle
   (finite element for plate problems)
   
   @param n - array of base functions
   @param l - area coordinates
   @param sx,sy - arrays of components
*/
void bf_cct (double *n,double *l,double *sx,double *sy)
{
  n[0]=l[0];
  n[1]=(l[0]*l[1]*sy[0]-l[2]*l[0]*sy[2])/2.0;
  n[2]=(l[2]*l[0]*sx[2]-l[0]*l[1]*sx[0])/2.0;

  n[3]=l[1];
  n[4]=(l[1]*l[2]*sy[1]-l[0]*l[1]*sy[0])/2.0;
  n[5]=(l[0]*l[1]*sx[0]-l[1]*l[2]*sx[1])/2.0;

  n[6]=l[2];
  n[7]=(l[2]*l[0]*sy[2]-l[1]*l[2]*sy[1])/2.0;
  n[8]=(l[1]*l[2]*sx[1]-l[2]*l[0]*sx[2])/2.0;
}



/**
   function computes derivates of base functions of constant curvature triangle with respect of x
   (finite element for plate problems)
   
   @param n - array of base functions
   @param l - area coordinates
   @param b - array of coefficients of area coordinates
   @param sx,sy - arrays of components
*/
void dx_cct (double *n,double *l,double *b,double *sx,double *sy)
{
  n[0]=b[0];
  n[1]=((b[0]*l[1]+l[0]*b[1])*sy[0]-(b[2]*l[0]+l[2]*b[0])*sy[2])/2.0;
  n[2]=((b[2]*l[0]+l[2]*b[0])*sx[2]-(b[0]*l[1]+l[0]*b[1])*sx[0])/2.0;

  n[3]=b[1];
  n[4]=((b[1]*l[2]+l[1]*b[2])*sy[1]-(b[0]*l[1]+l[0]*b[1])*sy[0])/2.0;
  n[5]=((b[0]*l[1]+l[0]*b[1])*sx[0]-(b[1]*l[2]+l[1]*b[2])*sx[1])/2.0;

  n[6]=b[2];
  n[7]=((b[2]*l[0]+l[2]*b[0])*sy[2]-(b[1]*l[2]+l[1]*b[2])*sy[1])/2.0;
  n[8]=((b[1]*l[2]+l[1]*b[2])*sx[1]-(b[2]*l[0]+l[2]*b[0])*sx[2])/2.0;
}



/**
   function computes derivates of base functions of constant curvature triangle with respect of y
   (finite element for plate problems)
   
   @param n - array of base functions
   @param l - area coordinates
   @param b - array of coefficients of area coordinates
   @param sx,sy - arrays of components
*/
void dy_cct (double *n,double *l,double *c,double *sx,double *sy)
{
  n[0]=c[0];
  n[1]=((c[0]*l[1]+l[0]*c[1])*sy[0]-(c[2]*l[0]+l[2]*c[0])*sy[2])/2.0;
  n[2]=((c[2]*l[0]+l[2]*c[0])*sx[2]-(c[0]*l[1]+l[0]*c[1])*sx[0])/2.0;

  n[3]=c[1];
  n[4]=((c[1]*l[2]+l[1]*c[2])*sy[1]-(c[0]*l[1]+l[0]*c[1])*sy[0])/2.0;
  n[5]=((c[0]*l[1]+l[0]*c[1])*sx[0]-(c[1]*l[2]+l[1]*c[2])*sx[1])/2.0;

  n[6]=c[2];
  n[7]=((c[2]*l[0]+l[2]*c[0])*sy[2]-(c[1]*l[2]+l[1]*c[2])*sy[1])/2.0;
  n[8]=((c[1]*l[2]+l[1]*c[2])*sx[1]-(c[2]*l[0]+l[2]*c[0])*sx[2])/2.0;
}


/**
   function computes base functions on quadrilateral plate element based
   on the discrete Kirchhoff theory
   
   nodal unknowns are ordered
   w, \phi_x = dw/dy, \phi_y = - dw/dx, ...

   @param n - array of base functions
   @param x,y - natural coordinates
   @param l - array of lengths of element edges
   @param sx,sy - array of x and y components of direction vectors
   @param nx,ny - array of x and y components of outer normal vectors
   
   JK, 20. 10. 2019
*/
void bf_quad_dkq (double *n,double x,double y,double *l,double *sx,double *sy,double *nx,double *ny)
{
  double *q;
  q = new double [8];
  
  //  bi-quadratic approximation functions on a quadrilateral element
  bf_quad_4_2d (q,x,y);
  
  // ********
  //  dw/dx = - phi_y
  // ********
  
  //  w1
  n[0] = 1.5*(q[4]*ny[0]/l[0]-q[7]*ny[3]/l[3]);
  //  dw1/dy = phi_x1
  n[1] = 0.75*q[4]*sy[0]*ny[0] + 0.75*q[7]*sy[3]*ny[3];
  //  dw1/dx = - phi_y1
  n[2] = 0.0 - q[0] - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]) - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]);

  //  w2
  n[3] = 1.5*(q[5]*ny[1]/l[1]-q[4]*ny[0]/l[0]);
  //  dw2/dy = phi_x2
  n[4] = 0.75*q[5]*sy[1]*ny[1] + 0.75*q[4]*sy[0]*ny[0];
  //  dw2/dx = - phi_y2
  n[5] = 0.0 - q[1] - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]) - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]);

  //  w3
  n[6] = 1.5*(q[6]*ny[2]/l[2]-q[5]*ny[1]/l[1]);
  //  dw3/dy = phi_x3
  n[7] = 0.75*q[6]*sy[2]*ny[2] + 0.75*q[5]*sy[1]*ny[1];
  //  dw3/dx = - phi_y3
  n[8] = 0.0 - q[2] - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]) - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]);

  //  w4
  n[9] = 1.5*(q[7]*ny[3]/l[3]-q[6]*ny[2]/l[2]);
  //  dw4/dy = phi_x4
  n[10] = 0.75*q[7]*sy[3]*ny[3] + 0.75*q[6]*sy[2]*ny[2];
  //  dw4/dx = - phi_y4
  n[11] = 0.0 - q[3] - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]) - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]);
  
  
  // ********
  //  dw/dy = phi_x
  // ********
  
  //  w1
  n[12] = 1.5*(q[7]*nx[3]/l[3]-q[4]*nx[0]/l[0]);
  //  dw1/dy = phi_x1
  n[13] = q[0] - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]) - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]);
  //  dw1/dx = - phi_y1
  n[14] = 0.75*q[4]*sx[0]*nx[0] + 0.75*q[7]*sx[3]*nx[3];

  //  w2
  n[15] = 1.5*(q[4]*nx[0]/l[0]-q[5]*nx[1]/l[1]);
  //  dw2/dy = phi_x2
  n[16] = q[1] - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]) - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]);
  //  dw2/dx = - phi_y2
  n[17] = 0.75*q[5]*sx[1]*nx[1] + 0.75*q[4]*sx[0]*nx[0];

  //  w3
  n[18] = 1.5*(q[5]*ny[1]/l[1]-q[6]*nx[2]/l[2]);
  //  dw3/dy = phi_x3
  n[19] = q[2] - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]) - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]);
  //  dw3/dx = - phi_y3
  n[20] = 0.75*q[6]*sx[2]*nx[2] + 0.75*q[5]*sx[1]*nx[1];

  //  w4
  n[21] = 1.5*(q[6]*nx[2]/l[2]-q[7]*nx[3]/l[3]);
  //  dw4/dy = phi_x4
  n[22] = q[3] - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]) - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]);
  //  dw4/dx = -phi_y4
  n[23] = 0.75*q[7]*sx[3]*nx[3] + 0.75*q[6]*sx[2]*nx[2];
  
  delete [] q;
}


/**
   function computes derivatives with respect to x of the base functions on quadrilateral plate element based
   on the discrete Kirchhoff theory
   
   nodal unknowns are ordered
   w, \phi_x = dw/dy, \phi_y = - dw/dx, ...
   
   @param n - array of base functions
   @param xi,eta - natural coordinates
   @param x,y - vectors of nodal coordinates
   @param l - array of lengths of element edges
   @param sx,sy - array of x and y components of direction vectors
   @param nx,ny - array of x and y components of outer normal vectors
   @param signs - array of signs
   
   JK, 20. 10. 2019
*/
void dx_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs)
{
  double jac;
  double *q;
  vector dx(ASTCKVEC(8)), dy(ASTCKVEC(8));
  q = new double [8];
  
  //  bi-quadratic approximation functions on a quadrilateral element
  dx_bf_quad_4_2d (dx.a,xi,eta);
  dy_bf_quad_4_2d (dy.a,xi,eta);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  for (long i=0;i<8;i++){
    q[i]=dx[i];
  }

  // ********
  //  d2w/dx2 = - d phi_y/dx
  // ********
  
  //  w1
  n[0] = 1.5*(signs[0]*q[4]*ny[0]/l[0]-signs[3]*q[7]*ny[3]/l[3]);
  //  dw1/dy = phi_x1
  n[1] = 0.75*q[4]*sy[0]*ny[0] + 0.75*q[7]*sy[3]*ny[3];
  //  dw1/dx = - phi_y1
  n[2] = 0.0 - q[0] - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]) - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]);

  //  w2
  n[3] = 1.5*(signs[1]*q[5]*ny[1]/l[1]-signs[0]*q[4]*ny[0]/l[0]);
  //  dw2/dy = phi_x2
  n[4] = 0.75*q[5]*sy[1]*ny[1] + 0.75*q[4]*sy[0]*ny[0];
  //  dw2/dx = - phi_y2
  n[5] = 0.0 - q[1] - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]) - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]);

  //  w3
  n[6] = 1.5*(signs[2]*q[6]*ny[2]/l[2]-signs[1]*q[5]*ny[1]/l[1]);
  //  dw3/dy = phi_x3
  n[7] = 0.75*q[6]*sy[2]*ny[2] + 0.75*q[5]*sy[1]*ny[1];
  //  dw3/dx = - phi_y3
  n[8] = 0.0 - q[2] - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]) - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]);

  //  w4
  n[9] = 1.5*(signs[3]*q[7]*ny[3]/l[3]-signs[2]*q[6]*ny[2]/l[2]);
  //  dw4/dy = phi_x4
  n[10] = 0.75*q[7]*sy[3]*ny[3] + 0.75*q[6]*sy[2]*ny[2];
  //  dw4/dx = - phi_y4
  n[11] = 0.0 - q[3] - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]) - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]);
  
  
  // ********
  //  d2w/dy/dx = d phi_x/dx
  // ********
  
  //  w1
  n[12] = 1.5*(signs[3]*q[7]*nx[3]/l[3]-signs[0]*q[4]*nx[0]/l[0]);
  //  dw1/dy = phi_x1
  n[13] = q[0] - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]) - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]);
  //  dw1/dx = - phi_y1
  n[14] = 0.75*q[4]*sx[0]*nx[0] + 0.75*q[7]*sx[3]*nx[3];

  //  w2
  n[15] = 1.5*(signs[0]*q[4]*nx[0]/l[0]-signs[1]*q[5]*nx[1]/l[1]);
  //  dw2/dy = phi_x2
  n[16] = q[1] - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]) - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]);
  //  dw2/dx = - phi_y2
  n[17] = 0.75*q[5]*sx[1]*nx[1] + 0.75*q[4]*sx[0]*nx[0];

  //  w3
  n[18] = 1.5*(signs[1]*q[5]*nx[1]/l[1]-signs[2]*q[6]*nx[2]/l[2]);
  //  dw3/dy = phi_x3
  n[19] = q[2] - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]) - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]);
  //  dw3/dx = - phi_y3
  n[20] = 0.75*q[6]*sx[2]*nx[2] + 0.75*q[5]*sx[1]*nx[1];

  //  w4
  n[21] = 1.5*(signs[2]*q[6]*nx[2]/l[2]-signs[3]*q[7]*nx[3]/l[3]);
  //  dw4/dy = phi_x4
  n[22] = q[3] - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]) - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]);
  //  dw4/dx = -phi_y4
  n[23] = 0.75*q[7]*sx[3]*nx[3] + 0.75*q[6]*sx[2]*nx[2];
  
  delete [] q;
}


/**
   function computes derivatives with respect to y of the base functions on quadrilateral plate element based
   on the discrete Kirchhoff theory
   
   nodal unknowns are ordered
   w, \phi_x = dw/dy, \phi_y = - dw/dx, ...
   
   @param n - array of base functions
   @param xi,eta - natural coordinates
   @param x,y - vectors of nodal coordinates
   @param l - array of lengths of element edges
   @param sx,sy - array of x and y components of direction vectors
   @param nx,ny - array of x and y components of outer normal vectors
   @param signs - array of signs
   
   JK, 20. 10. 2019
*/
void dy_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs)
{
  double jac;
  double *q;
  vector dx(ASTCKVEC(8)), dy(ASTCKVEC(8));
  q = new double [8];
  
  //  bi-quadratic approximation functions on a quadrilateral element
  dx_bf_quad_4_2d (dx.a,xi,eta);
  dy_bf_quad_4_2d (dy.a,xi,eta);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  for (long i=0;i<8;i++){
    q[i]=dy[i];
  }
  
  
  // ********
  //  d2w/dx/dy = - d phi_y/dy
  // ********
  
  //  w1
  n[0] = 1.5*(signs[0]*q[4]*ny[0]/l[0]-signs[3]*q[7]*ny[3]/l[3]);
  //  dw1/dy = phi_x1
  n[1] = 0.75*q[4]*sy[0]*ny[0] + 0.75*q[7]*sy[3]*ny[3];
  //  dw1/dx = - phi_y1
  n[2] = 0.0 - q[0] - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]) - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]);

  //  w2
  n[3] = 1.5*(signs[1]*q[5]*ny[1]/l[1]-signs[0]*q[4]*ny[0]/l[0]);
  //  dw2/dy = phi_x2
  n[4] = 0.75*q[5]*sy[1]*ny[1] + 0.75*q[4]*sy[0]*ny[0];
  //  dw2/dx = - phi_y2
  n[5] = 0.0 - q[1] - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]) - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]);

  //  w3
  n[6] = 1.5*(signs[2]*q[6]*ny[2]/l[2]-signs[1]*q[5]*ny[1]/l[1]);
  //  dw3/dy = phi_x3
  n[7] = 0.75*q[6]*sy[2]*ny[2] + 0.75*q[5]*sy[1]*ny[1];
  //  dw3/dx = - phi_y3
  n[8] = 0.0 - q[2] - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]) - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]);

  //  w4
  n[9] = 1.5*(signs[3]*q[7]*ny[3]/l[3]-signs[2]*q[6]*ny[2]/l[2]);
  //  dw4/dy = phi_x4
  n[10] = 0.75*q[7]*sy[3]*ny[3] + 0.75*q[6]*sy[2]*ny[2];
  //  dw4/dx = - phi_y4
  n[11] = 0.0 - q[3] - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]) - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]);
  
  
  // ********
  //  d2w/dy2 = d phi_x/dy
  // ********
  
  //  w1
  n[12] = 1.5*(signs[3]*q[7]*nx[3]/l[3]-signs[0]*q[4]*nx[0]/l[0]);
  //  dw1/dy = phi_x1
  n[13] = q[0] - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]) - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]);
  //  dw1/dx = - phi_y1
  n[14] = 0.75*q[4]*sx[0]*nx[0] + 0.75*q[7]*sx[3]*nx[3];

  //  w2
  n[15] = 1.5*(signs[0]*q[4]*nx[0]/l[0]-signs[1]*q[5]*nx[1]/l[1]);
  //  dw2/dy = phi_x2
  n[16] = q[1] - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]) - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]);
  //  dw2/dx = - phi_y2
  n[17] = 0.75*q[5]*sx[1]*nx[1] + 0.75*q[4]*sx[0]*nx[0];

  //  w3
  n[18] = 1.5*(signs[1]*q[5]*nx[1]/l[1]-signs[2]*q[6]*nx[2]/l[2]);
  //  dw3/dy = phi_x3
  n[19] = q[2] - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]) - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]);
  //  dw3/dx = - phi_y3
  n[20] = 0.75*q[6]*sx[2]*nx[2] + 0.75*q[5]*sx[1]*nx[1];

  //  w4
  n[21] = 1.5*(signs[2]*q[6]*nx[2]/l[2]-signs[3]*q[7]*nx[3]/l[3]);
  //  dw4/dy = phi_x4
  n[22] = q[3] - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]) - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]);
  //  dw4/dx = -phi_y4
  n[23] = 0.75*q[7]*sx[3]*nx[3] + 0.75*q[6]*sx[2]*nx[2];
  
  delete [] q;
}


/**
   function computes derivatives with respect to x of the base functions on quadrilateral plate element based
   on the discrete Kirchhoff theory
   
   nodal unknowns are ordered
   w, \phi_x = dw/dy, \phi_y = - dw/dx, ...
   
   @param n - array of base functions
   @param xi,eta - natural coordinates
   @param x,y - vectors of nodal coordinates
   @param l - array of lengths of element edges
   @param sx,sy - array of x and y components of direction vectors
   @param nx,ny - array of x and y components of outer normal vectors
   @param signs - array of signs
   
   JK, 20. 10. 2019
*/
void dxdx_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs)
{
  double *q;
  vector dxdx(ASTCKVEC(8)), dxdy(ASTCKVEC(8)), dydy(ASTCKVEC(8));
  q = new double [8];
  
  sec_derivatives_2d (batoz_plate,dxdx,dxdy,dydy,x,y,xi,eta);

  for (long i=0;i<8;i++){
    q[i]=dxdx[i];
  }
  

  // ********
  //  d3w/dx3 = - d2 phi_y/dx2
  // ********
  
  //  w1
  n[0] = 1.5*(signs[0]*q[4]*ny[0]/l[0]-signs[3]*q[7]*ny[3]/l[3]);
  //  dw1/dy = phi_x1
  n[1] = 0.75*q[4]*sy[0]*ny[0] + 0.75*q[7]*sy[3]*ny[3];
  //  dw1/dx = - phi_y1
  n[2] = 0.0 - q[0] - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]) - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]);

  //  w2
  n[3] = 1.5*(signs[1]*q[5]*ny[1]/l[1]-signs[0]*q[4]*ny[0]/l[0]);
  //  dw2/dy = phi_x2
  n[4] = 0.75*q[5]*sy[1]*ny[1] + 0.75*q[4]*sy[0]*ny[0];
  //  dw2/dx = - phi_y2
  n[5] = 0.0 - q[1] - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]) - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]);

  //  w3
  n[6] = 1.5*(signs[2]*q[6]*ny[2]/l[2]-signs[1]*q[5]*ny[1]/l[1]);
  //  dw3/dy = phi_x3
  n[7] = 0.75*q[6]*sy[2]*ny[2] + 0.75*q[5]*sy[1]*ny[1];
  //  dw3/dx = - phi_y3
  n[8] = 0.0 - q[2] - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]) - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]);

  //  w4
  n[9] = 1.5*(signs[3]*q[7]*ny[3]/l[3]-signs[2]*q[6]*ny[2]/l[2]);
  //  dw4/dy = phi_x4
  n[10] = 0.75*q[7]*sy[3]*ny[3] + 0.75*q[6]*sy[2]*ny[2];
  //  dw4/dx = - phi_y4
  n[11] = 0.0 - q[3] - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]) - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]);
  
  
  // ********
  //  d3w/dy/dx2 = d2 phi_x/dx2
  // ********

  //  w1
  n[12] = 1.5*(signs[3]*q[7]*nx[3]/l[3]-signs[0]*q[4]*nx[0]/l[0]);
  //  dw1/dy = phi_x1
  n[13] = q[0] - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]) - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]);
  //  dw1/dx = - phi_y1
  n[14] = 0.75*q[4]*sx[0]*nx[0] + 0.75*q[7]*sx[3]*nx[3];

  //  w2
  n[15] = 1.5*(signs[0]*q[4]*nx[0]/l[0]-signs[1]*q[5]*nx[1]/l[1]);
  //  dw2/dy = phi_x2
  n[16] = q[1] - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]) - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]);
  //  dw2/dx = - phi_y2
  n[17] = 0.75*q[5]*sx[1]*nx[1] + 0.75*q[4]*sx[0]*nx[0];

  //  w3
  n[18] = 1.5*(signs[1]*q[5]*nx[1]/l[1]-signs[2]*q[6]*nx[2]/l[2]);
  //  dw3/dy = phi_x3
  n[19] = q[2] - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]) - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]);
  //  dw3/dx = - phi_y3
  n[20] = 0.75*q[6]*sx[2]*nx[2] + 0.75*q[5]*sx[1]*nx[1];

  //  w4
  n[21] = 1.5*(signs[2]*q[6]*nx[2]/l[2]-signs[3]*q[7]*nx[3]/l[3]);
  //  dw4/dy = phi_x4
  n[22] = q[3] - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]) - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]);
  //  dw4/dx = -phi_y4
  n[23] = 0.75*q[7]*sx[3]*nx[3] + 0.75*q[6]*sx[2]*nx[2];
  
  delete [] q;
}

/**
   function computes derivatives with respect to y of the base functions on quadrilateral plate element based
   on the discrete Kirchhoff theory
   
   nodal unknowns are ordered
   w, \phi_x = dw/dy, \phi_y = - dw/dx, ...
   
   @param n - array of base functions
   @param xi,eta - natural coordinates
   @param x,y - vectors of nodal coordinates
   @param l - array of lengths of element edges
   @param sx,sy - array of x and y components of direction vectors
   @param nx,ny - array of x and y components of outer normal vectors
   @param signs - array of signs
   
   JK, 20. 10. 2019
*/
void dydy_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs)
{
  double *q;
  vector dxdx(ASTCKVEC(8)), dxdy(ASTCKVEC(8)), dydy(ASTCKVEC(8));
  q = new double [8];
  
  sec_derivatives_2d (batoz_plate,dxdx,dxdy,dydy,x,y,xi,eta);

  for (long i=0;i<8;i++){
    q[i]=dydy[i];
  }
  
  
  // ********
  //  d3w/dx/dy2 = - d2 phi_y/dy2
  // ********
  
  //  w1
  n[0] = 1.5*(signs[0]*q[4]*ny[0]/l[0]-signs[3]*q[7]*ny[3]/l[3]);
  //  dw1/dy = phi_x1
  n[1] = 0.75*q[4]*sy[0]*ny[0] + 0.75*q[7]*sy[3]*ny[3];
  //  dw1/dx = - phi_y1
  n[2] = 0.0 - q[0] - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]) - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]);

  //  w2
  n[3] = 1.5*(signs[1]*q[5]*ny[1]/l[1]-signs[0]*q[4]*ny[0]/l[0]);
  //  dw2/dy = phi_x2
  n[4] = 0.75*q[5]*sy[1]*ny[1] + 0.75*q[4]*sy[0]*ny[0];
  //  dw2/dx = - phi_y2
  n[5] = 0.0 - q[1] - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]) - q[4]*(0.5*sy[0]*nx[0] + 0.25*sx[0]*ny[0]);

  //  w3
  n[6] = 1.5*(signs[2]*q[6]*ny[2]/l[2]-signs[1]*q[5]*ny[1]/l[1]);
  //  dw3/dy = phi_x3
  n[7] = 0.75*q[6]*sy[2]*ny[2] + 0.75*q[5]*sy[1]*ny[1];
  //  dw3/dx = - phi_y3
  n[8] = 0.0 - q[2] - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]) - q[5]*(0.5*sy[1]*nx[1] + 0.25*sx[1]*ny[1]);

  //  w4
  n[9] = 1.5*(signs[3]*q[7]*ny[3]/l[3]-signs[2]*q[6]*ny[2]/l[2]);
  //  dw4/dy = phi_x4
  n[10] = 0.75*q[7]*sy[3]*ny[3] + 0.75*q[6]*sy[2]*ny[2];
  //  dw4/dx = - phi_y4
  n[11] = 0.0 - q[3] - q[7]*(0.5*sy[3]*nx[3] + 0.25*sx[3]*ny[3]) - q[6]*(0.5*sy[2]*nx[2] + 0.25*sx[2]*ny[2]);
  
  
  // ********
  //  d3w/dy3 = d2 phi_x/dy2
  // ********
  
  //  w1
  n[12] = 1.5*(signs[3]*q[7]*nx[3]/l[3]-signs[0]*q[4]*nx[0]/l[0]);
  //  dw1/dy = phi_x1
  n[13] = q[0] - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]) - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]);
  //  dw1/dx = - phi_y1
  n[14] = 0.75*q[4]*sx[0]*nx[0] + 0.75*q[7]*sx[3]*nx[3];

  //  w2
  n[15] = 1.5*(signs[0]*q[4]*nx[0]/l[0]-signs[1]*q[5]*nx[1]/l[1]);
  //  dw2/dy = phi_x2
  n[16] = q[1] - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]) - q[4]*(0.5*sx[0]*ny[0] + 0.25*sy[0]*nx[0]);
  //  dw2/dx = - phi_y2
  n[17] = 0.75*q[5]*sx[1]*nx[1] + 0.75*q[4]*sx[0]*nx[0];

  //  w3
  n[18] = 1.5*(signs[1]*q[5]*nx[1]/l[1]-signs[2]*q[6]*nx[2]/l[2]);
  //  dw3/dy = phi_x3
  n[19] = q[2] - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]) - q[5]*(0.5*sx[1]*ny[1] + 0.25*sy[1]*nx[1]);
  //  dw3/dx = - phi_y3
  n[20] = 0.75*q[6]*sx[2]*nx[2] + 0.75*q[5]*sx[1]*nx[1];

  //  w4
  n[21] = 1.5*(signs[2]*q[6]*nx[2]/l[2]-signs[3]*q[7]*nx[3]/l[3]);
  //  dw4/dy = phi_x4
  n[22] = q[3] - q[7]*(0.5*sx[3]*ny[3] + 0.25*sy[3]*nx[3]) - q[6]*(0.5*sx[2]*ny[2] + 0.25*sy[2]*nx[2]);
  //  dw4/dx = -phi_y4
  n[23] = 0.75*q[7]*sx[3]*nx[3] + 0.75*q[6]*sx[2]*nx[2];
  
  delete [] q;
}



/**
   function computes bi-cubic base functions on quadrilateral with 4 nodes
   x in (-1,1)
   y in (-1,1)
   the functions are used in DKQ element, where the displacement is defined only
   on element edges
   
   @param n - array of base functions
   @param x,y - natural coodinates
   
   JK, 30.11.2019
*/
void bf_cubic_plate_4_2d (double *n,double x,double y)
{
  //  node 1
  //  w1 = 1
  n[0]=(-0.25*x*x*x + 0.75*x + 0.5)*(-0.25*y*y*y + 0.75*y + 0.5);
  //  phi_x1 = 1
  n[1]=0.125*(y*y*y + y*y - y - 1.0)*(1.0 + x);
  //  phi_y1 = 1
  n[2]=0.125*(0.0-x*x*x - x*x + x + 1.0)*(1.0 + y);

  //  node 2
  //  w2 = 1
  n[3]=(0.25*x*x*x - 0.75*x + 0.5)*(-0.25*y*y*y + 0.75*y + 0.5);
  //  phi_x2 = 1
  n[4]=0.125*(y*y*y + y*y - y - 1.0)*(1.0 - x);
  //  phi_y2 = 1
  n[5]=0.125*(0.0-x*x*x + x*x + x - 1.0)*(1.0 + y);
  
  //  node 3
  //  w3 = 1
  n[6]=(0.25*x*x*x - 0.75*x + 0.5)*(0.25*y*y*y - 0.75*y + 0.5);
  //  phi_x3 = 1
  n[7]=0.125*(y*y*y - y*y - y + 1.0)*(1.0 - x);
  //  phi_y3 = 1
  n[8]=0.125*(0.0-x*x*x + x*x + x - 1.0)*(1.0 - y);

  //  node 4
  //  w4 = 1
  n[9]=(-0.25*x*x*x + 0.75*x + 0.5)*(0.25*y*y*y - 0.75*y + 0.5);
  //  phi_x4 = 1
  n[10]=0.125*(y*y*y - y*y - y + 1.0)*(1.0 + x);
  //  phi_y4 = 1
  n[11]=0.125*(0.0-x*x*x - x*x + x + 1.0)*(1.0 - y);
}



/**
   function computes coefficients a_i of volume coordinates
   L(x,y,z)_i = (a_i + b_i.x + c_i.y + d_i.z)/6/volume
   
   @param a - array of coefficients
   @param x,y,z - array of node coordinates
   @param v - volume of element
*/
void a_coeff_3d (double *a,double *x,double *y,double *z,double v)
{
  a[0] = x[1]*(y[1]*(z[3]-z[2])+y[2]*(z[1]-z[3])+y[3]*(z[2]-z[1]));
  a[0]+= y[1]*(z[1]*(x[3]-x[2])+z[2]*(x[1]-x[3])+z[3]*(x[2]-x[1]));
  a[0]+= z[1]*(x[1]*(y[3]-y[2])+x[2]*(y[1]-y[3])+x[3]*(y[2]-y[1]));
  a[0]/=6.0*v;

  a[1] = x[2]*(y[2]*(z[0]-z[3])+y[3]*(z[2]-z[0])+y[0]*(z[3]-z[2]));
  a[1]+= y[2]*(z[2]*(x[0]-x[3])+z[3]*(x[2]-x[0])+z[0]*(x[3]-x[2]));
  a[1]+= z[2]*(x[2]*(y[0]-y[3])+x[3]*(y[2]-y[0])+x[0]*(y[3]-y[2]));
  a[1]/=6.0*v;

  a[2] = x[3]*(y[3]*(z[1]-z[0])+y[0]*(z[3]-z[1])+y[1]*(z[0]-z[3]));
  a[2]+= y[3]*(z[3]*(x[1]-x[0])+z[0]*(x[3]-x[1])+z[1]*(x[0]-x[3]));
  a[2]+= z[3]*(x[3]*(y[1]-y[0])+x[0]*(y[3]-y[1])+x[1]*(y[0]-y[3]));
  a[2]/=6.0*v;

  a[3] = x[0]*(y[0]*(z[2]-z[1])+y[1]*(z[0]-z[2])+y[2]*(z[1]-z[0]));
  a[3]+= y[0]*(z[0]*(x[2]-x[1])+z[1]*(x[0]-x[2])+z[2]*(x[1]-x[0]));
  a[3]+= z[0]*(x[0]*(y[2]-y[1])+x[1]*(y[0]-y[2])+x[2]*(y[1]-y[0]));
  a[3]/=6.0*v;
}



/**
   function computes coefficients b_i of volume coordinates
   L(x,y,z)_i = (a_i + b_i.x + c_i.y + d_i.z)/6/volume
   
   @param b - array of coefficients
   @param y,z - array of node coordinates
   @param v - volume of element
*/
void b_coeff_3d (double *b,double *y,double *z,double v)
{

  b[0]=-((y[2]-y[1])*(z[3]-z[1])-(y[3]-y[1])*(z[2]-z[1]))/6.0/v;
  b[1]=((y[3]-y[2])*(z[0]-z[2])-(y[0]-y[2])*(z[3]-z[2]))/6.0/v;
  b[2]=-((y[0]-y[3])*(z[1]-z[3])-(y[1]-y[3])*(z[0]-z[3]))/6.0/v;
  b[3]=((y[1]-y[0])*(z[2]-z[0])-(y[2]-y[0])*(z[1]-z[0]))/6.0/v;
/*
  b[0]=((y[1]-y[3])*(z[2]-z[3])-(y[2]-y[3])*(z[1]-z[3]))/6.0/v;
  b[1]=((y[2]-y[0])*(z[3]-z[0])-(y[3]-y[0])*(z[2]-z[0]))/6.0/v;
  b[2]=((y[3]-y[1])*(z[0]-z[1])-(y[0]-y[1])*(z[3]-z[1]))/6.0/v;
  b[3]=((y[0]-y[2])*(z[1]-z[2])-(y[1]-y[2])*(z[0]-z[2]))/6.0/v;
*/
}



/**
   function computes coefficients c_i of volume coordinates
   L(x,y,z)_i = (a_i + b_i.x + c_i.y + d_i.z)/6/volume
   
   @param c - array of coefficients
   @param x,z - array of node coordinates
   @param v - volume of element
*/
void c_coeff_3d (double *c,double *x,double *z,double v)
{

  c[0]=-((x[3]-x[1])*(z[2]-z[1])-(x[2]-x[1])*(z[3]-z[1]))/6.0/v;
  c[1]=((x[0]-x[2])*(z[3]-z[2])-(x[3]-x[2])*(z[0]-z[2]))/6.0/v;
  c[2]=-((x[1]-x[3])*(z[0]-z[3])-(x[0]-x[3])*(z[1]-z[3]))/6.0/v;
  c[3]=((x[2]-x[0])*(z[1]-z[0])-(x[1]-x[0])*(z[2]-z[0]))/6.0/v;
/*
  c[0]=((x[2]-x[3])*(z[1]-z[3])-(x[1]-x[3])*(z[2]-z[3]))/6.0/v;
  c[1]=((x[3]-x[0])*(z[2]-z[0])-(x[2]-x[0])*(z[3]-z[0]))/6.0/v;
  c[2]=((x[0]-x[1])*(z[3]-z[1])-(x[3]-x[1])*(z[0]-z[1]))/6.0/v;
  c[3]=((x[1]-x[2])*(z[0]-z[2])-(x[0]-x[2])*(z[1]-z[2]))/6.0/v;
*/
}



/**
   function computes coefficients d_i of volume coordinates
   L(x,y,z)_i = (a_i + b_i.x + c_i.y + d_i.z)/6/volume
   
   @param d - array of coefficients
   @param x,y - array of node coordinates
   @param v - volume of element
*/
void d_coeff_3d (double *d,double *x,double *y,double v)
{

  d[0]=-((x[2]-x[1])*(y[3]-y[1])-(x[3]-x[1])*(y[2]-y[1]))/6.0/v;
  d[1]=((x[3]-x[2])*(y[0]-y[2])-(x[0]-x[2])*(y[3]-y[2]))/6.0/v;
  d[2]=-((x[0]-x[3])*(y[1]-y[3])-(x[1]-x[3])*(y[0]-y[3]))/6.0/v;
  d[3]=((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/6.0/v;
/*
  d[0]=((x[1]-x[3])*(y[2]-y[3])-(x[2]-x[3])*(y[1]-y[3]))/6.0/v;
  d[1]=((x[2]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[2]-y[0]))/6.0/v;
  d[2]=((x[3]-x[1])*(y[0]-y[1])-(x[0]-x[1])*(y[3]-y[1]))/6.0/v;
  d[3]=((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]))/6.0/v;
*/
}


/**
   function computes coefficients a_i of volume coordinates
   L(x,y,z)_i = a_i + b_i.x + c_i.y + d_i.z
   
   @param a - array of coefficients
   @param x,y,z - array of node coordinates
   @param det - determinat=6*volume
*/
void vola_3d (double *a,double *x,double *y,double *z,double det)
{
  a[0]=(x[1]*(y[2]*z[3]-y[3]*z[2]) + x[2]*(y[3]*z[1]-y[1]*z[3]) + x[3]*(y[1]*z[2]-y[2]*z[1]))/det;
  a[1]=(x[0]*(y[3]*z[2]-y[2]*z[3]) + x[2]*(y[0]*z[3]-y[3]*z[0]) + x[3]*(y[2]*z[0]-y[0]*z[2]))/det;
  a[2]=(x[0]*(y[1]*z[3]-y[3]*z[1]) + x[1]*(y[3]*z[0]-y[0]*z[3]) + x[3]*(y[0]*z[1]-y[1]*z[0]))/det;
  a[3]=(x[0]*(y[2]*z[1]-y[1]*z[2]) + x[1]*(y[0]*z[2]-y[2]*z[0]) + x[2]*(y[1]*z[0]-y[0]*z[1]))/det;
}



/**
   function computes coefficients b_i of volume coordinates
   L(x,y,z)_i = a_i + b_i.x + c_i.y + d_i.z
   
   @param b - array of coefficients
   @param y,z - array of node coordinates
   @param det - determinat=6*volume
*/
void volb_3d (double *b,double *y,double *z,double det)
{
  b[0]=(y[1]*(z[3]-z[2]) + y[2]*(z[1]-z[3]) + y[3]*(z[2]-z[1]))/det;
  b[1]=(y[0]*(z[2]-z[3]) + y[2]*(z[3]-z[0]) + y[3]*(z[0]-z[2]))/det;
  b[2]=(y[0]*(z[3]-z[1]) + y[1]*(z[0]-z[3]) + y[3]*(z[1]-z[0]))/det;
  b[3]=(y[0]*(z[1]-z[2]) + y[1]*(z[2]-z[0]) + y[2]*(z[0]-z[1]))/det;
}



/**
   function computes coefficients c_i of volume coordinates
   L(x,y,z)_i = a_i + b_i.x + c_i.y + d_i.z
   
   @param c - array of coefficients
   @param x,z - array of node coordinates
   @param det - determinat=6*volume
*/
void volc_3d (double *c,double *x,double *z,double det)
{
  c[0]=(x[1]*(z[2]-z[3]) + x[2]*(z[3]-z[1]) + x[3]*(z[1]-z[2]))/det;
  c[1]=(x[0]*(z[3]-z[2]) + x[2]*(z[0]-z[3]) + x[3]*(z[2]-z[0]))/det;
  c[2]=(x[0]*(z[1]-z[3]) + x[1]*(z[3]-z[0]) + x[3]*(z[0]-z[1]))/det;
  c[3]=(x[0]*(z[2]-z[1]) + x[1]*(z[0]-z[2]) + x[2]*(z[1]-z[0]))/det;
}

/**
   function computes coefficients d_i of volume coordinates
   L(x,y,z)_i = a_i + b_i.x + c_i.y + d_i.z
   
   @param d - array of coefficients
   @param x,y - array of node coordinates
   @param det - determinat=6*volume
*/
void vold_3d (double *d,double *x,double *y,double det)
{
  d[0]=(x[1]*(y[3]-y[2]) + x[2]*(y[1]-y[3]) + x[3]*(y[2]-y[1]))/det;
  d[1]=(x[0]*(y[2]-y[3]) + x[2]*(y[3]-y[0]) + x[3]*(y[0]-y[2]))/det;
  d[2]=(x[0]*(y[3]-y[1]) + x[1]*(y[0]-y[3]) + x[3]*(y[1]-y[0]))/det;
  d[3]=(x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]))/det;
}



/**
   function computes tri-linear base functions on tetrahedron
   
   @param n - array of base functions
   @param x,y,z - natural coodinates
*/
void bf_lin_tet (double *n,double x,double y,double z)
{
  n[0]=x;
  n[1]=y;
  n[2]=z;
  n[3]=1.0-x-y-z;
}



void dx_bf_lin_tet (double *n)
{
  n[0]=1.0;
  n[1]=0.0;
  n[2]=0.0;
  n[3]=-1.0;
}



void dy_bf_lin_tet (double *n)
{
  n[0]=0.0;
  n[1]=1.0;
  n[2]=0.0;
  n[3]=-1.0;
}



void dz_bf_lin_tet (double *n)
{
  n[0]=0.0;
  n[1]=0.0;
  n[2]=1.0;
  n[3]=-1.0;
}



/**
   function computes tri-quadratic base functions on tetrahedron
   
   @param n - array of base functions
   @param x,y,z - natural coodinates
   TKo - 11.2008
*/
void bf_quad_tet(double *n,double x,double y,double z)
{
  n[0]=x*(x-0.5)*2.0;
  n[1]=y*(y-0.5)*2.0;
  n[2]=z*(z-0.5)*2.0;
  n[3]=(1.0-x-y-z)*(0.5-x-y-z)*2.0;
  n[4]=x*y*4.0;
  n[5]=y*z*4.0;
  n[6]=x*z*4.0;
  n[7]=x*(1.0-x-y-z)*4.0;
  n[8]=y*(1.0-x-y-z)*4.0;
  n[9]=z*(1.0-x-y-z)*4.0;
}



/**
   function computes derivatives of tri-quadratic base functions on tetrahedron with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
   TKo - 11.2008
*/
void dx_bf_quad_tet(double *n,double x,double y,double z)
{
  n[0] = 4.0*x-1.0;
  n[1] = 0.0;
  n[2] = 0.0;
  n[3] = 4.0*(x+y+z)-3.0;

  n[4] = 4.0*y;
  n[5] = 0.0;
  n[6] = 4.0*z;

  n[7] =  4.0*(1.0-2.0*x-y-z);
  n[8] = -4.0*y;
  n[9] = -4.0*z;
}



/**
   function computes derivatives of tri-quadratic base functions on tetrahedron with respect of natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
   TKo - 11.2008
*/
void dy_bf_quad_tet(double *n,double x,double y,double z)
{
  n[0] = 0.0;
  n[1] = 4.0*y-1.0;
  n[2] = 0.0;
  n[3] = 4.0*(x+y+z)-3.0;

  n[4] = 4.0*x;
  n[5] = 4.0*z;
  n[6] = 0.0;

  n[7] = -4.0*x;
  n[8] =  4.0*(1.0-x-2.0*y-z);
  n[9] = -4.0*z;
}



/**
   function computes derivatives of tri-quadratic base functions on tetrahedron with respect of natural coordinate z
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
   TKo - 11.2008
*/
void dz_bf_quad_tet (double *n,double x,double y,double z)
{
  n[0] = 0.0;
  n[1] = 0.0;
  n[2] = 4.0*z-1.0;
  n[3] = 4.0*(x+y+z)-3.0;

  n[4] = 0.0;
  n[5] = 4.0*y;
  n[6] = 4.0*x;

  n[7] = -4.0*x;
  n[8] = -4.0*y;
  n[9] =  4.0*(1.0-x-y-2.0*z);
}



/**
   function computes linear base functions on wedge
   
   @param n - array of basis functions
   @param x,y,z - natural coordinates
   
   JK, 9.5.2004
*/
void bf_lin_wed_3d (double *n,double x,double y,double z)
{
  n[0]=x*(1.0+z)/2.0;
  n[1]=y*(1.0+z)/2.0;
  n[2]=(1.0-x-y)*(1.0+z)/2.0;
  n[3]=x*(1.0-z)/2.0;
  n[4]=y*(1.0-z)/2.0;
  n[5]=(1.0-x-y)*(1.0-z)/2.0;
}

/**
   function computes derivatives of linear base functions on wedge with respect to x
   
   @param n - array of basis functions
   @param z - natural coordinate
   
   JK, 9.5.2004
*/
void dx_bf_lin_wed_3d (double *n,double z)
{
  n[0]=(1.0+z)/2.0;
  n[1]=0.0;
  n[2]=-1.0*(1.0+z)/2.0;
  n[3]=(1.0-z)/2.0;
  n[4]=0.0;
  n[5]=-1.0*(1.0-z)/2.0;
}

/**
   function computes derivatives of linear base functions on wedge with respect to y
   
   @param n - array of basis functions
   @param z - natural coordinate
   
   JK, 9.5.2004
*/
void dy_bf_lin_wed_3d (double *n,double z)
{
  n[0]=0.0;
  n[1]=(1.0+z)/2.0;
  n[2]=-1.0*(1.0+z)/2.0;
  n[3]=0.0;
  n[4]=(1.0-z)/2.0;
  n[5]=-1.0*(1.0-z)/2.0;
}

/**
   function computes derivatives of linear base functions on wedge with respect to z
   
   @param n - array of basis functions
   @param x,y - natural coordinates
   
   JK, 9.5.2004
*/
void dz_bf_lin_wed_3d (double *n,double x,double y)
{
  n[0]=x/2.0;
  n[1]=y/2.0;
  n[2]=(1.0-x-y)/2.0;
  n[3]=x/(-2.0);
  n[4]=y/(-2.0);
  n[5]=(1.0-x-y)/(-2.0);
}






/**
   function computes quadratic approximation functions on wedge
   
   @param n - array of basis functions
   @param x,y,z - natural coordinates
   
   JK, 19.9.2004
*/
void bf_quad_wed_3d (double *n,double x,double y,double z)
{
  n[0]  = x*(x-0.5)*(1.0+z)*z;
  n[1]  = y*(y-0.5)*(1.0+z)*z;
  n[2]  = (1.0-x-y)*(0.5-x-y)*(1.0+z)*z;
  n[3]  = x*(x-0.5)*(z-1.0)*z;
  n[4]  = y*(y-0.5)*(z-1.0)*z;
  n[5]  = (1.0-x-y)*(0.5-x-y)*(z-1.0)*z;
  n[6]  = 2.0*x*y*(1.0+z);
  n[7]  = 2.0*y*(1.0-x-y)*(1.0+z);
  n[8]  = 2.0*x*(1.0-x-y)*(1.0+z);
  n[9]  = x*(1.0-z*z);
  n[10] = y*(1.0-z*z);
  n[11] = (1.0-x-y)*(1.0-z*z);
  n[12] = 2.0*x*y*(1.0-z);
  n[13] = 2.0*y*(1.0-x-y)*(1.0-z);
  n[14] = 2.0*x*(1.0-x-y)*(1.0-z);
}

/**
   function computes derivatives of linear approximation functions on wedge with respect to x
   
   @param n - array of basis functions
   @param x,y,z - natural coordinates
   
   JK, 19.9.2004
*/
void dx_bf_quad_wed_3d (double *n,double x,double y,double z)
{
  n[0]  = (2.0*x-0.5)*(1.0+z)*z;
  n[1]  = 0.0;
  n[2]  = (2.0*x+2.0*y-1.5)*(1.0+z)*z;
  n[3]  = (2.0*x-0.5)*(z-1.0)*z;
  n[4]  = 0.0;
  n[5]  = (2.0*x+2.0*y-1.5)*(z-1.0)*z;
  n[6]  = 2.0*y*(1.0+z);
  n[7]  = -2.0*y*(1.0+z);
  n[8]  = (2.0-4.0*x-2.0*y)*(1.0+z);
  n[9]  = (1.0-z*z);
  n[10] = 0.0;
  n[11] = z*z-1.0;
  n[12] = 2.0*y*(1.0-z);
  n[13] = 2.0*y*(z-1.0);
  n[14] = (2.0-4.0*x-2.0*y)*(1.0-z);
}

/**
   function computes derivatives of quadratic approximation functions on wedge with respect to y
   
   @param n - array of basis functions
   @param x,y,z - natural coordinates
   
   JK, 19.9.2004
*/
void dy_bf_quad_wed_3d (double *n,double x,double y,double z)
{
  n[0]  = 0.0;
  n[1]  = (2.0*y-0.5)*(1.0+z)*z;
  n[2]  = (2.0*x+2.0*y-1.5)*(1.0+z)*z;
  n[3]  = 0.0;
  n[4]  = (2.0*y-0.5)*(z-1.0)*z;
  n[5]  = (2.0*x+2.0*y-1.5)*(z-1.0)*z;
  n[6]  = 2.0*x*(1.0+z);
  n[7]  = (2.0-2.0*x-4.0*y)*(1.0+z);
  n[8]  = -2.0*x*(1.0+z);
  n[9]  = 0.0;
  n[10] = 1.0-z*z;
  n[11] = z*z-1.0;
  n[12] = 2.0*x*(1.0-z);
  n[13] = (2.0-2.0*x-4.0*y)*(1.0-z);
  n[14] = 2.0*x*(z-1.0);
}

/**
   function computes derivatives of quadratic approximation functions on wedge with respect to z
   
   @param n - array of basis functions
   @param x,y,z - natural coordinates
   
   JK, 19.9.2004
*/
void dz_bf_quad_wed_3d (double *n,double x,double y,double z)
{
  n[0]  = x*(x-0.5)*(1.0+2.0*z);
  n[1]  = y*(y-0.5)*(1.0+2.0*z);
  n[2]  = (1.0-x-y)*(0.5-x-y)*(1.0+2.0*z);
  n[3]  = x*(x-0.5)*(2.0*z-1.0);
  n[4]  = y*(y-0.5)*(2.0*z-1.0);
  n[5]  = (1.0-x-y)*(0.5-x-y)*(2.0*z-1.0);
  n[6]  = 2.0*x*y;
  n[7]  = 2.0*y*(1.0-x-y);
  n[8]  = 2.0*x*(1.0-x-y);
  n[9]  = -2.0*x*z;
  n[10] = -2.0*y*z;
  n[11] = -2.0*z*(1.0-x-y);
  n[12] = -2.0*x*y;
  n[13] = -2.0*y*(1.0-x-y);
  n[14] = -2.0*x*(1.0-x-y);
}







/**
   function computes tri-linear base functions on hexahedron
   
   @param n - array of base functions
   @param x,y,z - natural coodinates
*/
void bf_lin_hex_3d (double *n,double x,double y,double z)
{
  n[0]=(1.0+x)*(1.0+y)*(1.0+z)/8.0;
  n[1]=(1.0-x)*(1.0+y)*(1.0+z)/8.0;
  n[2]=(1.0-x)*(1.0-y)*(1.0+z)/8.0;
  n[3]=(1.0+x)*(1.0-y)*(1.0+z)/8.0;
  n[4]=(1.0+x)*(1.0+y)*(1.0-z)/8.0;
  n[5]=(1.0-x)*(1.0+y)*(1.0-z)/8.0;
  n[6]=(1.0-x)*(1.0-y)*(1.0-z)/8.0;
  n[7]=(1.0+x)*(1.0-y)*(1.0-z)/8.0;
}

/**
   function computes derivatives of tri-linear base functions on hexahedron with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
   @param y,z - natural coodinates
*/
void dx_bf_lin_hex_3d (double *n,double y,double z)
{
  n[0]=(1.0+y)*(1.0+z)/8.0;
  n[1]=(1.0+y)*(1.0+z)/(-8.0);
  n[2]=(1.0-y)*(1.0+z)/(-8.0);
  n[3]=(1.0-y)*(1.0+z)/8.0;
  n[4]=(1.0+y)*(1.0-z)/8.0;
  n[5]=(1.0+y)*(1.0-z)/(-8.0);
  n[6]=(1.0-y)*(1.0-z)/(-8.0);
  n[7]=(1.0-y)*(1.0-z)/8.0;
}



/**
   function computes derivatives of tri-linear base functions on hexahedron with respect of natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x,z - natural coodinates
*/
void dy_bf_lin_hex_3d (double *n,double x,double z)
{
  n[0]=(1.0+x)*(1.0+z)/8.0;     
  n[1]=(1.0-x)*(1.0+z)/8.0;
  n[2]=(1.0-x)*(1.0+z)/(-8.0);
  n[3]=(1.0+x)*(1.0+z)/(-8.0);
  n[4]=(1.0+x)*(1.0-z)/8.0;
  n[5]=(1.0-x)*(1.0-z)/8.0;
  n[6]=(1.0-x)*(1.0-z)/(-8.0);
  n[7]=(1.0+x)*(1.0-z)/(-8.0);
}



/**
   function computes derivatives of tri-linear base functions on hexahedron with respect of natural coordinate z
   
   @param n - array of derivatives of base functions
   @param x,y - natural coodinates
*/
void dz_bf_lin_hex_3d (double *n,double x,double y)
{
  n[0]=(1.0+x)*(1.0+y)/8.0;
  n[1]=(1.0-x)*(1.0+y)/8.0;
  n[2]=(1.0-x)*(1.0-y)/8.0;
  n[3]=(1.0+x)*(1.0-y)/8.0;
  n[4]=(1.0+x)*(1.0+y)/(-8.0);
  n[5]=(1.0-x)*(1.0+y)/(-8.0);
  n[6]=(1.0-x)*(1.0-y)/(-8.0);
  n[7]=(1.0+x)*(1.0-y)/(-8.0);
}



/**
   function computes tri-quadratic base functions on hexahedron
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
*/
void bf_quad_hex_3d (double *n,double x,double y,double z)
{
  n[0]  = (1.0+x)*(1.0+y)*(1.0+z)*(x+y+z-2.0)/8.0;
  n[1]  = (1.0-x)*(1.0+y)*(1.0+z)*(y-x+z-2.0)/8.0;
  n[2]  = (1.0-x)*(1.0-y)*(1.0+z)*(z-x-y-2.0)/8.0;
  n[3]  = (1.0+x)*(1.0-y)*(1.0+z)*(x-y+z-2.0)/8.0;

  n[4]  = (1.0+x)*(1.0+y)*(1.0-z)*(x+y-z-2.0)/8.0;
  n[5]  = (1.0-x)*(1.0+y)*(1.0-z)*(y-x-z-2.0)/8.0;
  n[6]  = (1.0-x)*(1.0-y)*(1.0-z)*(x+y+z+2.0)/(-8.0);
  n[7]  = (1.0+x)*(1.0-y)*(1.0-z)*(x-y-z-2.0)/8.0;

  n[8]  = (1.0-x*x)*(1.0+y)*(1.0+z)/4.0;
  n[9]  = (1.0-y*y)*(1.0-x)*(1.0+z)/4.0;
  n[10] = (1.0-x*x)*(1.0-y)*(1.0+z)/4.0;
  n[11] = (1.0-y*y)*(1.0+x)*(1.0+z)/4.0;

  n[12] = (1.0+x)*(1.0+y)*(1.0-z*z)/4.0;
  n[13] = (1.0-x)*(1.0+y)*(1.0-z*z)/4.0;
  n[14] = (1.0-x)*(1.0-y)*(1.0-z*z)/4.0;
  n[15] = (1.0+x)*(1.0-y)*(1.0-z*z)/4.0;

  n[16] = (1.0-x*x)*(1.0+y)*(1.0-z)/4.0;
  n[17] = (1.0-y*y)*(1.0-x)*(1.0-z)/4.0;
  n[18] = (1.0-x*x)*(1.0-y)*(1.0-z)/4.0;
  n[19] = (1.0-y*y)*(1.0+x)*(1.0-z)/4.0;
}



/**
   function computes derivatives of tri-quadratic base functions on hexahedron with respect of natural coordinate x
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
*/
void dx_bf_quad_hex_3d (double *n,double x,double y,double z)
{
  n[0]  = ((1.0+y)*(1.0+z)*(x+y+z-2.0) + (1.0+x)*(1.0+y)*(1.0+z))/8.0;
  n[1]  = ((1.0+y)*(1.0+z)*(y-x+z-2.0) + (1.0-x)*(1.0+y)*(1.0+z))/(-8.0);
  n[2]  = ((1.0-y)*(1.0+z)*(z-x-y-2.0) + (1.0-x)*(1.0-y)*(1.0+z))/(-8.0);
  n[3]  = ((1.0-y)*(1.0+z)*(x-y+z-2.0) + (1.0+x)*(1.0-y)*(1.0+z))/8.0;

  n[4]  = ((1.0+y)*(1.0-z)*(x+y-z-2.0) + (1.0+x)*(1.0+y)*(1.0-z))/8.0;
  n[5]  = ((1.0+y)*(1.0-z)*(y-x-z-2.0) + (1.0-x)*(1.0+y)*(1.0-z))/(-8.0);
  n[6]  = ((1.0-y)*(1.0-z)*(x+y+z+2.0) - (1.0-x)*(1.0-y)*(1.0-z))/8.0;
  n[7]  = ((1.0-y)*(1.0-z)*(x-y-z-2.0) + (1.0+x)*(1.0-y)*(1.0-z))/8.0;

  n[8]  = x*(1.0+y)*(1.0+z)/(-2.0);
  n[9]  = (1.0-y*y)*(1.0+z)/(-4.0);
  n[10] = x*(1.0-y)*(1.0+z)/(-2.0);
  n[11] = (1.0-y*y)*(1.0+z)/4.0;

  n[12] = (1.0+y)*(1.0-z*z)/4.0;
  n[13] = (1.0+y)*(1.0-z*z)/(-4.0);
  n[14] = (1.0-y)*(1.0-z*z)/(-4.0);
  n[15] = (1.0-y)*(1.0-z*z)/4.0;

  n[16] = x*(1.0+y)*(1.0-z)/(-2.0);
  n[17] = (1.0-y*y)*(1.0-z)/(-4.0);
  n[18] = x*(1.0-y)*(1.0-z)/(-2.0);
  n[19] = (1.0-y*y)*(1.0-z)/4.0;
}



/**
   function computes derivatives of tri-quadratic base functions on hexahedron with respect of natural coordinate y
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
*/
void dy_bf_quad_hex_3d (double *n,double x,double y,double z)
{
  n[0]  = ((1.0+x)*(1.0+z)*(x+y+z-2.0) + (1.0+x)*(1.0+y)*(1.0+z))/8.0;
  n[1]  = ((1.0-x)*(1.0+z)*(y-x+z-2.0) + (1.0-x)*(1.0+y)*(1.0+z))/8.0;
  n[2]  = ((1.0-x)*(1.0+z)*(z-x-y-2.0) + (1.0-x)*(1.0-y)*(1.0+z))/(-8.0);
  n[3]  = ((1.0+x)*(1.0+z)*(x-y+z-2.0) + (1.0+x)*(1.0-y)*(1.0+z))/(-8.0);

  n[4]  = ((1.0+x)*(1.0-z)*(x+y-z-2.0) + (1.0+x)*(1.0+y)*(1.0-z))/8.0;
  n[5]  = ((1.0-x)*(1.0-z)*(y-x-z-2.0) + (1.0-x)*(1.0+y)*(1.0-z))/8.0;
  n[6]  = ((1.0-x)*(1.0-z)*(x+y+z+2.0) - (1.0-x)*(1.0-y)*(1.0-z))/8.0;
  n[7]  = ((1.0+x)*(1.0-z)*(x-y-z-2.0) + (1.0+x)*(1.0-y)*(1.0-z))/(-8.0);

  n[8]  = (1.0-x*x)*(1.0+z)/4.0;
  n[9]  = y*(1.0-x)*(1.0+z)/(-2.0);
  n[10] = (1.0-x*x)*(1.0+z)/(-4.0);
  n[11] = y*(1.0+x)*(1.0+z)/(-2.0);

  n[12] = (1.0+x)*(1.0-z*z)/4.0;
  n[13] = (1.0-x)*(1.0-z*z)/4.0;
  n[14] = (1.0-x)*(1.0-z*z)/(-4.0);
  n[15] = (1.0+x)*(1.0-z*z)/(-4.0);

  n[16] = (1.0-x*x)*(1.0-z)/4.0;
  n[17] = y*(1.0-x)*(1.0-z)/(-2.0);
  n[18] = (1.0-x*x)*(1.0-z)/(-4.0);
  n[19] = y*(1.0+x)*(1.0-z)/(-2.0);
}

/**
   function computes derivatives of tri-quadratic base functions on hexahedron with respect of natural coordinate z
   
   @param n - array of derivatives of base functions
   @param x,y,z - natural coodinates
*/
void dz_bf_quad_hex_3d (double *n,double x,double y,double z)
{
  n[0]  = ((1.0+x)*(1.0+y)*(x+y+z-2.0) + (1.0+x)*(1.0+y)*(1.0+z))/8.0;
  n[1]  = ((1.0-x)*(1.0+y)*(y-x+z-2.0) + (1.0-x)*(1.0+y)*(1.0+z))/8.0;
  n[2]  = ((1.0-x)*(1.0-y)*(z-x-y-2.0) + (1.0-x)*(1.0-y)*(1.0+z))/8.0;
  n[3]  = ((1.0+x)*(1.0-y)*(x-y+z-2.0) + (1.0+x)*(1.0-y)*(1.0+z))/8.0;

  n[4]  = ((1.0+x)*(1.0+y)*(x+y-z-2.0) + (1.0+x)*(1.0+y)*(1.0-z))/(-8.0);
  n[5]  = ((1.0-x)*(1.0+y)*(y-x-z-2.0) + (1.0-x)*(1.0+y)*(1.0-z))/(-8.0);
  n[6]  = ((1.0-x)*(1.0-y)*(x+y+z+2.0) - (1.0-x)*(1.0-y)*(1.0-z))/8.0;
  n[7]  = ((1.0+x)*(1.0-y)*(x-y-z-2.0) + (1.0+x)*(1.0-y)*(1.0-z))/(-8.0);

  n[8]  = (1.0-x*x)*(1.0+y)/4.0;
  n[9]  = (1.0-y*y)*(1.0-x)/4.0;
  n[10] = (1.0-x*x)*(1.0-y)/4.0;
  n[11] = (1.0-y*y)*(1.0+x)/4.0;

  n[12] = (1.0+x)*(1.0+y)*z/(-2.0);
  n[13] = (1.0-x)*(1.0+y)*z/(-2.0);
  n[14] = (1.0-x)*(1.0-y)*z/(-2.0);
  n[15] = (1.0+x)*(1.0-y)*z/(-2.0);

  n[16] = (1.0-x*x)*(1.0+y)/(-4.0);
  n[17] = (1.0-y*y)*(1.0-x)/(-4.0);
  n[18] = (1.0-x*x)*(1.0-y)/(-4.0);
  n[19] = (1.0-y*y)*(1.0+x)/(-4.0);
}


void bf_quad_hexrot_3d (double *n,double x,double y,double z)
{
  
  n[0]  = (1.0-x*x)*(1.0+y)*(1.0+z)/4.0;
  n[1]  = (1.0-y*y)*(1.0-x)*(1.0+z)/4.0;
  n[2] = (1.0-x*x)*(1.0-y)*(1.0+z)/4.0;
  n[3] = (1.0-y*y)*(1.0+x)*(1.0+z)/4.0;

  n[4] = (1.0+x)*(1.0+y)*(1.0-z*z)/4.0;
  n[5] = (1.0-x)*(1.0+y)*(1.0-z*z)/4.0;
  n[6] = (1.0-x)*(1.0-y)*(1.0-z*z)/4.0;
  n[7] = (1.0+x)*(1.0-y)*(1.0-z*z)/4.0;

  n[8] = (1.0-x*x)*(1.0+y)*(1.0-z)/4.0;
  n[9] = (1.0-y*y)*(1.0-x)*(1.0-z)/4.0;
  n[10] = (1.0-x*x)*(1.0-y)*(1.0-z)/4.0;
  n[11] = (1.0-y*y)*(1.0+x)*(1.0-z)/4.0;
}

void dx_bf_quad_hexrot_3d (double *n,double x,double y,double z)
{
  //  side 1-2, 2-3, 3-4, 4-1,  1-5, 2-6, 3-7, 4-8,  5-6, 6-7, 7-8, 8-5,   

  n[0]  = x*(1.0+y)*(1.0+z)/(-2.0);
  n[1]  = (1.0-y*y)*(1.0+z)/(-4.0);
  n[2] = x*(1.0-y)*(1.0+z)/(-2.0);
  n[3] = (1.0-y*y)*(1.0+z)/4.0;

  n[4] = (1.0+y)*(1.0-z*z)/4.0;
  n[5] = (1.0+y)*(1.0-z*z)/(-4.0);
  n[6] = (1.0-y)*(1.0-z*z)/(-4.0);
  n[7] = (1.0-y)*(1.0-z*z)/4.0;

  n[8] = x*(1.0+y)*(1.0-z)/(-2.0);
  n[9] = (1.0-y*y)*(1.0-z)/(-4.0);
  n[10] = x*(1.0-y)*(1.0-z)/(-2.0);
  n[11] = (1.0-y*y)*(1.0-z)/4.0;
}
void dy_bf_quad_hexrot_3d (double *n,double x,double y,double z)
{
  //  side 1-2, 2-3, 3-4, 4-1,  1-5, 2-6, 3-7, 4-8,  5-6, 6-7, 7-8, 8-5,   

  n[0]  = (1.0-x*x)*(1.0+z)/4.0;
  n[1]  = y*(1.0-x)*(1.0+z)/(-2.0);
  n[2] = (1.0-x*x)*(1.0+z)/(-4.0);
  n[3] = y*(1.0+x)*(1.0+z)/(-2.0);

  n[4] = (1.0+x)*(1.0-z*z)/4.0;
  n[5] = (1.0-x)*(1.0-z*z)/4.0;
  n[6] = (1.0-x)*(1.0-z*z)/(-4.0);
  n[7] = (1.0+x)*(1.0-z*z)/(-4.0);

  n[8] = (1.0-x*x)*(1.0-z)/4.0;
  n[9] = y*(1.0-x)*(1.0-z)/(-2.0);
  n[10] = (1.0-x*x)*(1.0-z)/(-4.0);
  n[11] = y*(1.0+x)*(1.0-z)/(-2.0);

}
void dz_bf_quad_hexrot_3d (double *n,double x,double y,double z)
{
  //  side 1-2, 2-3, 3-4, 4-1,  1-5, 2-6, 3-7, 4-8,  5-6, 6-7, 7-8, 8-5,   

  n[0]  = (1.0-x*x)*(1.0+y)/4.0;
  n[1]  = (1.0-y*y)*(1.0-x)/4.0;
  n[2] = (1.0-x*x)*(1.0-y)/4.0;
  n[3] = (1.0-y*y)*(1.0+x)/4.0;

  n[4] = (1.0+x)*(1.0+y)*z/(-2.0);
  n[5] = (1.0-x)*(1.0+y)*z/(-2.0);
  n[6] = (1.0-x)*(1.0-y)*z/(-2.0);
  n[7] = (1.0+x)*(1.0-y)*z/(-2.0);

  n[8] = (1.0-x*x)*(1.0+y)/(-4.0);
  n[9] = (1.0-y*y)*(1.0-x)/(-4.0);
  n[10] = (1.0-x*x)*(1.0-y)/(-4.0);
  n[11] = (1.0-y*y)*(1.0+x)/(-4.0);
}



/**
  The function returns approximated global coordinates of the point given by natural coordinates (xi, eta) 
  for 2D elements given by nodal coordinates stored in vectors x, y.
  
  @param p - %vector conatining approximated function values (output)
  @param x, y, z - vectors of nodal coordinates of the 3D element
  @param xi, eta, zeta - natural coordinates
*/
void bf_3d(vector &p, vector &x, vector &y, vector &z, double xi, double eta, double zeta)
{
  vector n(ASTCKVEC(x.n));
  switch (x.n)
  {
    case 4:
      bf_lin_tet(n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    case 6:
      bf_lin_wed_3d(n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    case 8:
      bf_lin_hex_3d (n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    case 10:
      bf_quad_tet(n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    case 15:
      bf_quad_wed_3d (n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    case 20:
      bf_quad_hex_3d (n.a, xi, eta, zeta);
      scprd(n, x, p(0));
      scprd(n, y, p(1));
      scprd(n, z, p(2)); 
      break;
    default:
      print_err("unknown number of approximation functions (nap=%ld) for bar element", __FILE__, __LINE__, __func__, n.n);
  }
}



/**
  The function checks whether the given natural coordinates xi, eta, zeta are within the 
  admisible range of the given element type et. If coordinates do not match the conditions, 
  then they are set to the closest point on the given element type.


  @param[in] et - type of element 
  @param[in/out] xi, eta, zeta - checked/resulting natural coordinates of the given point

  Created by TKo 07.2020
*/
void corr_nat_coord_bounds(gtypel et, double &xi, double &eta, double &zeta)
{
  switch (et)
  {
    case isolinear1d:
    case isoquadratic1d:
      if (xi < -1.0)   xi = -1.0;
      if (xi >  1.0)   xi =  1.0;
      break;
    case trianglelinear:
    case trianglequadratic:
      corr_nat_coord_tr(xi, eta);
      break;
    case isolinear2d:
    case isoquadratic2d:
    case isocubic2d:
      if (xi < -1.0)   xi = -1.0;
      if (xi >  1.0)   xi =  1.0;
      if (eta < -1.0)  eta = -1.0;
      if (eta >  1.0)  eta =  1.0;
      break;
    case tetrahedronlinear:
    case tetrahedronquadratic:
      corr_nat_coord_tet(xi, eta, zeta);
      break;
    case isolinear3d:
    case isoquadratic3d:
      if (xi < -1.0)    xi = -1.0;
      if (xi >  1.0)    xi =  1.0;
      if (eta < -1.0)   eta = -1.0;
      if (eta >  1.0)   eta =  1.0;
      if (zeta < -1.0)  zeta = -1.0;
      if (zeta >  1.0)  zeta =  1.0;
      break;
    default:
      print_err("unknown type (%d) of element is required", __FILE__, __LINE__, __func__, et);
      abort();
  }
}



/**
  The function checks whether the given natural coordinates xi, eta, zeta are within the 
  admisible range on the unit triangular element. If coordinates do not match the conditions, 
  then they are set to the closest point on unit triangle.


  @param[in/out] xi, eta - checked/resulting natural coordinates of the given point

  Created by TKo 07.2020
*/
void corr_nat_coord_tr(double &xi, double &eta)
{
  if ((xi  >= 0.0) && (xi <= 1.0) &&
      (eta >= 0.0) && (eta <= 1.0) &&
      (xi+eta <= 1.0))
    return;  // correct natural coordinates
  else{
    if ((xi < 0.0) && (eta < 0.0)){
      // projection on the vertex in the natrual coordinate origin
      xi = 0.0; eta = 0.0;
      return;
    }
    if ((xi < 0.0) && (eta >= 0.0) && (eta <= 1.0)){
      // projection on the edge in eta axis
      xi = 0.0;
      return;
    }
    if ((eta < 0.0) && (xi >= 0.0) && (xi <= 1.0)){
      // projection on the edge in xi axis
      eta = 0.0;
      return;
    }
    if ((xi+eta > 1.0) && (eta-xi < 1.0+xi) && (eta > xi-1.0)){
      // coordinates are above the diagonal edge and in the strip formed
      // by normals to the diagonal edge in unit verteces =>
      // => projection on the diagonal edge
      xi  = 0.5*(xi - eta + 1.0);
      eta = 0.5*(eta - xi + 1.0);
      return;
    }
    if (xi > 1.0){
      // projection to the unit vertex on xi axis
      xi = 1.0; eta = 0.0;
      return;
    }
    if (eta > 1.0){
      // projection to the unit vertex on eta axis
      eta = 1.0; xi = 0.0;
      return;
    }

    // this means that something was not considered in the above list of conditions
    print_err("internal error: this case should not happen (xi=%le, eta=%le).",
              __FILE__, __LINE__, __func__, xi, eta);
    abort();
  }
}



/**
  The function checks whether the given natural coordinates xi, eta, zeta are within the 
  admisible range of the unit tehtrahedron element. If coordinates do not match the conditions, 
  then they are set to the closest point on the unit tetraherdon.


  @param[in/out] xi, eta, zeta - checked/resulting natural coordinates of the given point

  Created by TKo 07.2020
*/
void corr_nat_coord_tet(double &xi, double &eta, double &zeta)
{
  double pxi, peta, pzeta, t;

  if ((xi  >= 0.0)  && (xi   <= 1.0) &&
      (eta >= 0.0)  && (eta  <= 1.0) &&
      (zeta >= 0.0) && (zeta <= 1.0) &&
      (xi+eta+zeta <= 1.0))
    return; // correct natural coordinates
  else{
    if ((xi < 0.0) && (eta < 0.0) && (zeta < 0.0)){
      xi = 0.0; eta = 0.0; zeta = 0.0;  // projection to the vertex in the natural coordinate origin
      return;
    }
    if ((xi < 0.0) && (eta < 0.0)){
      xi = 0.0; eta = 0.0;  // projection to the edge on zeta axis
      if (zeta > 1.0){
        zeta = 1.0;  // projection to the unit vertex on zeta axis
      }
      return;
    }
    if ((xi < 0.0) && (zeta < 0.0)){
      xi = 0.0; zeta = 0.0;  // projection to the edge on eta axis
      if (eta > 1.0){ 
        eta = 1.0;  // projection to the unit vertex on eta axis
      }
      return;
    }
    if ((eta < 0.0) && (zeta < 0.0)){
      eta = 0.0; zeta = 0.0;  // projection to the edge on xi axis
      if (xi > 1.0){
        xi = 1.0;   // projection to the unit vertex on xi axis
      }
      return;
    }
    if (xi < 0.0){
      xi = 0.0;   // projection to the eta-zeta plane
      if (eta+zeta > 1.0){
        // projection on the diagonal edge in eta-zeta plane
        peta  = 0.5*(eta - zeta + 1.0);
        pzeta = 0.5*(zeta - eta + 1.0);
        eta  = peta;
        zeta = pzeta;
        if (peta > 1.0){
          eta = 1.0;   zeta = 0.0; // projection to the unit vertex on eta
        }
        if (pzeta > 1.0){
          zeta = 1.0;  eta = 0.0;  // projection to the unit vertex on zeta
        }
      }
      return;
    }
    if (eta < 0.0){
      eta = 0.0; // projection on xi-zeta plane
      if (xi+zeta > 1.0){
        // projection on the diagonal edge in xi-zeta plane
        pxi   = 0.5*(xi - zeta + 1.0);
        pzeta = 0.5*(zeta - xi + 1.0);
        xi   = pxi;
        zeta = pzeta;
        if (pxi > 1.0){
          xi = 1.0;    zeta = 0.0; // projection to the unit vertex on xi
        }
        if (pzeta > 1.0){
          zeta = 1.0;  xi = 0.0; // projection to the unit vertex on zeta
        }
      }
      return;
    }
    if (zeta < 0.0){
      zeta = 0.0; // projection on xi-eta plane
      if (xi+eta > 1.0){
        // projection on the diagonal edge in xi-eta plane
        pxi   = 0.5*(xi - eta + 1.0);
        peta  = 0.5*(eta - xi + 1.0);
        xi  = pxi;
        eta = peta;
        if (pxi > 1.0){
          xi = 1.0;   eta = 0.0; // projection to the unit vertex on xi
        }
        if (peta > 1.0){
          eta = 1.0;  xi = 0.0; // projection to the unit vertex on eta
        }
      }
      return;
    }
    if ((xi >= 0.0) && (eta >= 0.0) && (zeta >= 0.0)){
      // projection on the plane with normal vector (1,1,1) formed by equilateral triangle surface
      t = 1.0/3.0*(1.0 - xi - eta - zeta);
      xi   += t;
      eta  += t;
      zeta += t;
      if ((xi >= 0.0) && (eta >= 0.0) && (zeta >= 0.0)) // projected point lies on the face formed by equilateral triangle
        return;
      else
        // the case when at least one coordinate xi, eta or zeta is negative => one level of recursion
        corr_nat_coord_tet(xi, eta, zeta);
      return;
    }

    // this means that something was not considered in the above list of conditions
    print_err("internal error: this case should not happen (xi=%le, eta=%le, zeta=%le).",
              __FILE__, __LINE__, __func__, xi, eta, zeta);
    abort();
  }
}
