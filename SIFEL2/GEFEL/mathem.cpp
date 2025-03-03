#include "mathem.h"
#include "basefun.h"
#include "gtopology.h"
#include "vector.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
//#include <complex.h>
#define MATH_ZERO 1.0e-100

/**
   function computes radius (distance from axisis of symmetry)
   of required point with coordinates xi,eta
   it is used in axisymmetric problems
   
   @param x - %vector containing real coordinates
   @param natcoord - %vector of natural coordinates
          for quadrilaterals
          natcoord[0]=xi, natcoord[1]=eta;
          for triangles
          natcoord[0]=areacoord[0], natcoord[1]=areacoord[1], natcoord[2]=areacoord[2];
   
   JK, 8.12.2001
*/
double radius (const vector &x, const vector &natcoord)
{
  long nne;
  double r;
  nne=x.n;
  vector bf(nne);
  
  switch (nne){
  case 3:{
    bf[0]=natcoord[0];
    bf[1]=natcoord[1];
    bf[2]=natcoord[2];
    break;
  }
  case 4:{
    bf_lin_4_2d (bf.a,natcoord[0],natcoord[1]);
    break;
  }
  default:{
    print_err("unknown number of nodes is required", __FILE__, __LINE__, __func__);
  }
  }
  
  scprd (x,bf,r);
  return r;
}

/**
   function computes distance between points A and B

   @param coorda - coordinates of point A
   @param coordb - coordinates of point B
   
   30.12.2001
*/
double length (const vector &coorda, const vector &coordb)
{
  long i,m,n;
  double l;
  
  m=coorda.n;  n=coordb.n;
  if (m!=n)  
    print_err("non-compatible numbers of components", __FILE__, __LINE__, __func__);
  
  l=0.0;
  for (i=0;i<m;i++){
    l+=(coorda[i]-coordb[i])*(coorda[i]-coordb[i]);
  }
  l=sqrt (l);

  return l;
}

/**
   function computes distance between points A and B

   @param[in] coorda - coordinates of point A
   @param[in] xb,yb,zb - coordinates of point B
   
   Cretaed by TKo 07.2020
*/
double length (const vector &coorda, double xb, double yb, double zb)
{
  double l;
   
  l  = sqr(xb-coorda[0]);
  l += sqr(yb-coorda[1]);
  l += sqr(zb-coorda[2]);
  l  = sqrt(l);

  return l;
}


double sqr(double x)
{
  return (x*x);
}


double pow3(double x)
{
  return (x*x*x);
}



/**
  The function computes volume of a hexahedron element with linear approximation functions. The element 
  is defined by nodal coordinates given in arrays x, y and z.

  @param[in] x - array of x coordinates  of hexahedron vertices (nodes), its dimension must be 8.
  @param[in] y - array of y coordinates  of hexahedron vertices (nodes), its dimension must be 8.
  @param[in] z - array of z coordinates  of hexahedron vertices (nodes), its dimension must be 8.

  @return The function returns the element volume.

  Created by Tomas Koudelka, 01.2024
*/
double linhex_volume(vector &x, vector &y, vector &z)
{
  double vol  = ((z[2] + z[3] - z[4] - z[5]) * y[1] + (-z[1] + z[3]) * y[2] + (-z[1] - z[2] + z[4] + z[7]) * y[3] + (z[1] - z[3] + z[5] - z[7]) * y[4] + (z[1] - z[4]) * y[5] - y[7] * (z[3] - z[4])) * x[0];
  vol += ((-z[2] - z[3] + z[4] + z[5]) * y[0] + (z[0] + z[3] - z[5] - z[6]) * y[2] + (z[0] - z[2]) * y[3] + (-z[0] + z[5]) * y[4] + (-z[0] + z[2] - z[4] + z[6]) * y[5] + y[6] * (z[2] - z[5])) * x[1];
  vol += ((z[1] - z[3]) * y[0] + (-z[0] - z[3] + z[5] + z[6]) * y[1] + (z[0] + z[1] - z[6] - z[7]) * y[3] + (-z[1] + z[6]) * y[5] + (-z[1] + z[3] - z[5] + z[7]) * y[6] + y[7] * (z[3] - z[6])) * x[2];
  vol += ((z[1] + z[2] - z[4] - z[7]) * y[0] + (-z[0] + z[2]) * y[1] + (-z[0] - z[1] + z[6] + z[7]) * y[2] + (z[0] - z[7]) * y[4] + (-z[2] + z[7]) * y[6] + y[7] * (z[0] - z[2] + z[4] - z[6])) * x[3];
  vol += ((-z[1] + z[3] - z[5] + z[7]) * y[0] + (z[0] - z[5]) * y[1] + (-z[0] + z[7]) * y[3] + (z[0] + z[1] - z[6] - z[7]) * y[5] + (z[5] - z[7]) * y[6] - y[7] * (z[0] + z[3] - z[5] - z[6])) * x[4];
  vol += ((-z[1] + z[4]) * y[0] + (z[0] - z[2] + z[4] - z[6]) * y[1] + (z[1] - z[6]) * y[2] + (-z[0] - z[1] + z[6] + z[7]) * y[4] + (z[1] + z[2] - z[4] - z[7]) * y[6] - y[7] * (z[4] - z[6])) * x[5];
  vol += ((-z[2] + z[5]) * y[1] + (z[1] - z[3] + z[5] - z[7]) * y[2] + (z[2] - z[7]) * y[3] + (-z[5] + z[7]) * y[4] + (-z[1] - z[2] + z[4] + z[7]) * y[5] + y[7] * (z[2] + z[3] - z[4] - z[5])) * x[6];
  vol -= x[7] * ((-z[3] + z[4]) * y[0] + (z[3] - z[6]) * y[2] + (z[0] - z[2] + z[4] - z[6]) * y[3] + (-z[0] - z[3] + z[5] + z[6]) * y[4] + (-z[4] + z[6]) * y[5] + y[6] * (z[2] + z[3] - z[4] - z[5]));

  vol /= 12.0;

  return vol;
}



/**
   function computes nodal values with help of least square problem solution
   
   @param val - array containing nodal values on element
   @param nx,ny,nz - arrays containing natural coordinates of nodes
   @param lhs - array of coefficients of linear functions (solution of least square problem)
   @param dim - problem dimension
   @param fi - first index
   @param ncomp - number of stored components

   15.7.2002
*/
void nodal_values (double **val,vector &nx,vector &ny,vector &nz,
                   double *lhs,long dim,long fi,long ncomp)
{
  long i,j,jj;

  if (dim==1){
    for (i=0;i<nx.n;i++){
      jj=fi;
      for (j=0;j<ncomp;j++){
        val[i][jj] = lhs[j*2]+lhs[j*2+1]*nx[i];
        jj++;
      }
    }
  }
  if (dim==2){
    for (i=0;i<nx.n;i++){
      jj=fi;
      for (j=0;j<ncomp;j++){
        val[i][jj] = lhs[j*3]+lhs[j*3+1]*nx[i]+lhs[j*3+2]*ny[i];
        jj++;
      }
    }
  }
  if (dim==3){
    for (i=0;i<nx.n;i++){
      jj=fi;
      for (j=0;j<ncomp;j++){
        val[i][jj] = lhs[j*4]+lhs[j*4+1]*nx[i]+lhs[j*4+2]*ny[i]+lhs[j*4+3]*nz[i];
        jj++;
      }
    }
  }
  
}



/**
   function returns the signum of given value (if value is < 0 returns -1, otherwise returns 1)
 */
double sgn (double i)
{
  return (i< 0.0 ? -1.0 : 1.0);
}


/**
   Heaviside function
   
   for x<0, H(x)=0
   for x>=0, H(x)=1;
   
   @param x - variable
   
   3. 10. 2013, JK
 */
double heaviside (double x)
{
  if (x<0.0)
    return 0.0;
  else
    return 1.0;

}



/**
  The function decomposes real number arg into mantisa and exponent, i.e.
  arg = ret * 10^exp

  @param[in] arg - real number to be decomposed
  @param[out] exp - exponent in power of 10

  @return The function returns mantisa of the required number arg.

  Created by Tomas Koudelka, 01.2024
*/
double frexp10(double arg, long &exp)
{
  if (arg == 0.0)
    exp = 0;
  else
    exp = 1 + (long)floor(log10(fabs(arg)));

  return arg*pow(10, -exp);
}


/**
   generalized Heaviside function
   
   for x<0, H(x)=0
   for x>=0, H(x)=1;
   
   @param x - variable
   @param eps - radius of transition zone for the Heaviside function
   
   3. 10. 2013, JK
 */
double genheaviside (double x,double eps)
{
  double h;
  
  if (x<0.0-eps)
    h= 0.0;
  if (0.0-eps <= x && x <= eps)
    h = (x+eps)/2.0/eps;
  if (eps < x)
    h = 1.0;

  //fprintf (stdout,"\n Heaviside function  %lf %lf",x,h);
  
  return h;
  //return heaviside (x);  
}



/**
   function establishs to polynom of 4th order
 */
double polynom_4 (double x,double *a)
{
  return (a[0]*x*x*x*x + a[1]*x*x*x + a[2]*x*x + a[3]*x + a[4]);
}




/**
   Function sorts first two members of array with ascending order.
   
   @param x - sorted array
   
   created  27.3.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void sort_2 (long *x)
{
  if (x[0] < x[1])
    return;
  long a = x[0];  x[0] = x[1];  x[1] = a;
}



/**
   Function sorts first three members of array with ascending order.
   
   @param x - sorted array
   
   created  27.3.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void sort_3 (long *x)
{
  long a;
  
  if (x[0] > x[1])
  {
    if (x[0] > x[2])  
    {
      if (x[1] > x[2]) 
      { 
        a=x[0];  x[0]=x[2];  x[2]=a; 
      }
      else             
      {  
        a=x[0];  x[0]=x[1];  x[1]=x[2];  x[2]=a; 
      }
    }
    else
    {  
      a=x[0];  x[0]=x[1];  x[1]=a; 
    }
  }
  else
  {              
    if (x[1] > x[2])
    {  
      if (x[0] > x[2]) 
      {  
        a=x[0];  x[0]=x[2];  x[2]=x[1];  x[1]=a; 
      }
      else
      {  
        a=x[1];  x[1]=x[2];  x[2]=a; 
      }
    }
  }
}



/**
   Function sorts first four members of array with ascending order
   
   @param x - sorted array
   
   created  27.3.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void sort_4 (long *x)
{
  long a;
  
  if (x[0] < x[1])  
  {
    if (x[1] > x[2])  
    {
      if (x[1] > x[3]) 
      {  
        a=x[1];  x[1]=x[3];  x[3]=a;  sort_3(x);  
        return; 
      }
      else 
      {  
        a=x[1];  x[1]=x[2];  x[2]=a;  
      }
    }
    else
    { 
      if (x[2] > x[3]) 
      {  
        a=x[2];  x[2]=x[3];  x[3]=a;  sort_3(x);  
        return; 
      }
      else
      {  
        return; 
      }
    }
  }
  else
  {
    if (x[0] > x[2])  
    {
      if (x[0] > x[3]) 
      {  
        a=x[0];  x[0]=x[3];  x[3]=a;  sort_3(x);  
        return; 
      }
      else
      {  
        a=x[0];  x[0]=x[2];  x[2]=a;  
      }
    }
    else 
    {
      if (x[2] > x[3]) 
      {  
        a=x[2];  x[2]=x[3];  x[3]=a;  sort_3(x);  
        return; 
      }
    }
  }
  if (x[0]>x[1]) 
  { 
    a=x[0]; x[0]=x[1]; x[1]=a; 
  }
}



/**
   Function sorts first three members of array with ascending order.
   
   @param x - sorted array
   
   created  12.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void sort_3 (double *x)
{
  double a;
  
  if (x[0] > x[1])  
  {
    if (x[0] > x[2])  
    {
      if (x[1] > x[2]) 
      {  
        a=x[0];  x[0]=x[2];  x[2]=a; 
      }
      else
      {  
        a=x[0];  x[0]=x[1];  x[1]=x[2];  x[2]=a; 
      }
    }
    else                               
    {  
      a=x[0];  x[0]=x[1];  x[1]=a; 
    }
  }
  else
  { 
    if (x[1] > x[2])
    {
      if (x[0] > x[2]) 
      {  
        a=x[0];  x[0]=x[2];  x[2]=x[1];  x[1]=a; 
      }
      else
      {
        a=x[1];  x[1]=x[2];  x[2]=a; 
      }
    }
  }
}



/**
  The function solves roots of quadratic equation a*x^2 + b*x +c = 0.

  @param[in] a - coefficient of the second order term
  @param[in] b - coefficient of the first order term
  @param[in] c - coefficient of the zero order term
  @param[in] zero - zero treshold value
  @param[out] x - array for result storage 

  @retval -1 - infinite number of roots
  @retval  0 - no roots
  @retval  1 - one real root
  @retval -1 - two real roots

  Created by Tomas Koudelka, 30.10.2018
*/
long solve_quadratic(double a, double b, double c, double zero, double x[2])
{
  if (fabs(a) < zero){
    if (fabs(b) < zero){
      if (fabs(c) < zero)
        return -1;
      else
        return 0;
    }
    else{
      x[0] = -c/b;
      return 1;
    }
  }
  else{
    double d = sqr(b) - 4.0*a*c;
    if (d < 0.0)
      return 0;
    else{
      double aux = 0.5/a;
      x[0] = (-b + sqrt(d))*aux;
      x[1] = (-b - sqrt(d))*aux;
    }
  }
  return 2;
} 



/**
  The function solves roots of cubic equation a*x^3 + b*x^2 + c*x + d = 0.

  @param[in] a - coefficient of the third order term
  @param[in] b - coefficient of the second order term
  @param[in] c - coefficient of the first order term
  @param[in] d - coefficient of the zero order term
  @param[in] zero - zero treshold value
  @param[out] nreal - the number of real roots
  @param[out] x - array for result storage 
  @param[out] ncompl - the number of complex roots
  @param[out] z - array for result storage 
*/
/* void solve_cubic(double a, double b, double c, double d, double zero, long &nreal, double *x, long &ncompl, double complex *z)
{
  nreal = ncompl = 0L;

  if (fabs(a) < zero){
    if (fabs(b) < zero){
      if (fabs(c) < zero){
        if (fabs(d) < zero){
          nreal = -1;
          return;
        }
        else
          return;
      }
      else{
        x[0] = -c/b;
        nreal = 1;
      }
    }
    else{
      nreal = solve_quadratic(b, c, d, zero, x);
      return;
    }
  }
  else{
    // solve cubic equation in the form x^3 + aa*x^2 + bb*x + cc = 0
    double aa = b/a;
    double bb = c/a;
    double cc = d/a;
    // transform above equation by substitution x = t -aa/3 in the from t^3 + p*t + q = 0
    double p = bb - pow3(aa)/3.0;
    double q = cc + (2.0*pow3(aa) - 9.0*aa*bb)/27.0;

    double dd = 18.0*a*b*c*d - 4.0*pow3(b)*d + sqr(b)*sqr(c) - 4.0*a*pow3(c) - 27.0*sqr(a)*sqr(d);
    if (dd >= 0.0)
      nreal = 3;
    else{
      nreal = 1;
      ncompl = 2;
    }
    if (dd < 0.0)
      return;
    else{
      double aux = 0.5/a;
      x[0] = (-b + sqrt(dd))*aux;
      x[1] = (-b - sqrt(dd))*aux;
    }
  }
  return;
} 
*/


/**
   function searchs roots of polynom of 2nd order = quadratic equation
   a*x^2 + b*x + c = 0
   answer: -1 = infinite number of results
           0,1,2 = no,one,two results
   
   @param a,b,c - coeficients of quadratic equation
   @param r1,r2 - roots of quadratic equation
   
   created  11.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long solv_polynom_2 (double a,double b,double c,double &r1,double &r2)
{
  if (fabs(a) < MATH_ZERO)
  {
    if (fabs(b) < MATH_ZERO)
    {
      if (fabs(c) < MATH_ZERO)
        return -1;
      else
        return 0;
    }
    else
    {
      r1 = -c/b;
      return 1;
    }
  }
  double D = b*b - 4*a*c;
  
  if (D < -MATH_ZERO)
    return 0;
  
  if (-MATH_ZERO<D && D<MATH_ZERO)
  {
    r1 = -b/2.0/a;
    return 1;
  }
  
  r1 = (-b + sqrt(D))/2.0/a;
  r2 = (-b - sqrt(D))/2.0/a;
  return 2;
}

/**
   function searchs roots of polynom of 3rd order = cubic equation
   ax^3 + bx^2 + cx + d = 0
   answer: -1 = infinite number of results
           0,1,2,3 = no,one,two,three results
   
   @param a,b,c,d - coeficients of cubic equation
   @param r1,r2,r3 - roots of cubic equation
   
   created  xx.x.xxxx, Borek Pacak, bp@cml.fsv.cvut.cz
 */
long solv_polynom_3 (double a, double b, double c, double d, double &r1, double &r2, double &r3)
{
  if(fabs(a) < MATH_ZERO)
    return solv_polynom_2 (b,c,d,r1,r2);
  
  
  double aa, p, q, D, pq, u, v, phi;
  double help,pi;
  
  pi=3.141592653589793238;
  
  aa = a;
  a = b/aa;
  b = c/aa;
  c = d/aa;
  p = b - a * a / 3.0;
  q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  pq = p * q;
  D = q * q / 4.0 + p * p * p / 27.0;
  
  if (fabs(D) < MATH_ZERO){
    if (fabs(pq) < MATH_ZERO){
      r1 = 0.0 - a / 3.0;
      r2 = r1;
      r3 = r1;
      return 3;
    }
    
    if (q < 0.0)
      r2 = -exp(log(-q / 2.0) / 3.0);
    else
      r2 = exp(log(q / 2.0) / 3.0);
    
    r1 = -2.0 * r2 - a / 3.0;
    r2 -= a / 3.0;
    return 2;
  }
  
  if (D > 0.0){
    u = -q / 2.0 + sqrt(D);
    v = -q / 2.0 - sqrt(D);
    if (u < 0.0)
      u = -exp(log(-u) / 3.0);
    else
      u = exp(log(u) / 3.0);
    if (v < 0.0)
      v = -exp(log(-v) / 3.0);
    else
      v = exp(log(v) / 3.0);
    r1 = u + v - a / 3.0;
    return 1;
  }
  
  p = sqrt(fabs(p) / 3.0);
  help = (-q / (2.0 * p * p * p));
  if (fabs (help) > 1.0) help = sgn(help); // prevent rounding errors
  
  phi = acos(help) / 3.0;
  r1 = 2.0 * p * cos(phi) - a / 3.0;
  r2 = -2.0 * p * cos(phi - pi / 3.0) - a / 3.0;
  r3 = -2.0 * p * cos(phi + pi / 3.0) - a / 3.0;
  return 3;
}

/**
   function searchs roots of polynom of 4th order in interval <a,b>
   coeff[0]*x^4 + coeff[1]*x^3 + coeff[2]*x^2 + coeff[3]*x + coeff[4] = 0
   answer: -1 = infinite number of results
           0,1,2,3,4 = no,one,two,three,four results
   
   @param coeff - array of coefficients of polynom
   @param a,b - interval of searching
   @param acc - relative accuracy of roots => (exact_root - root)/exact_root is in (-acc;acc) ;;  acc > 1e-14
   @param roots - roots of polynom
   
   created  1.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long solv_polynom_4 (double *coeff,double a,double b,double acc,double *roots)
{
  long i,j,ninflps,nint,pol,maxit,nr;
  double inflp[5],x,y,ya,yb;
  
  maxit = 200;
  
  ninflps = solv_polynom_3 (4*coeff[0],3*coeff[1],2*coeff[2],coeff[3],inflp[1],inflp[2],inflp[3]);
  sort_3 (inflp+1);
  
  nint = 1;
  inflp[0] = a;
  for (i=0;i<ninflps;i++)
    if ( a < inflp[i+1]  &&  inflp[i+1] < b )  inflp[nint++] = inflp[i+1];
  
  inflp[nint] = b;
  
  nr = 0;
  for (i=0;i<nint;i++){
    a = inflp[i];
    b = inflp[i+1];
    ya = polynom_4 (a,coeff);
    yb = polynom_4 (b,coeff);
    
    if      (fabs(ya) < MATH_ZERO)  roots[nr++] = a; 
    else if (fabs(yb) < MATH_ZERO)  roots[nr++] = b;
    else if (ya*yb<0) {
      pol = (ya < 0.0 ? 1 : -1);
      
      for (j=0;j<maxit;j++){
        
        switch (j%3){
        case 0:{  x = (a+b)/2.0;                       break; }
        case 1:{  x =  a - ya*(b-a)/(yb-ya);           break; }
        case 2:{  x = (a - ya*(b-a)/(yb-ya))*2.0 - x;  break; }
        }
        
        y = polynom_4 (x,coeff);
        if (0.0 < pol*y) { b = x; yb = y; }
        else             { a = x; ya = y; }
        
        if ( 2.0*(b-a)/(b+a) < acc ){
          roots[nr++] = x;
          break;
        }
      }
    }
    else continue;
    
    if (j==maxit){
      fprintf (stderr,"\n\n\n  ERROR 1 in function solv_polynom_4, PLEASE CONTACT AUTHOR of this function(%s, line %d)\n\n\n",__FILE__,__LINE__);
      exit (0);
    }
    
    if ( nr!=1 && roots[nr-2]+acc > roots[nr-1])  nr--; 
  }
  
  return nr;
}



/**
   function solves linear equation:
   a*x + b = 0
   answer: -1 = infinite number of results
           0,1 = number of results
   
   @param a,b - coeficients of equation
   @param x - result
   
   created  4.4.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long solv_1le (double a,double b,double &x)
{
  if (fabs(a) < MATH_ZERO)
  {
    if (fabs(b) < MATH_ZERO)
      return -1;
    else
      return 0;
  }
  
  x = -b/a;
  
  return 1;
}



/**
   function solves system of two linear equations:
   a[0]*x + b[0]*y + c[0] = 0
   a[1]*x + b[1]*y + c[1] = 0
   answer: -1 = infinite number of results
           0,1 = number of results
   
   @param a,b,c - arrays of coeficients of equations, dimension is 2
   @param x,y - results
   
   created  4.4.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long solv_2le (double *a, double *b, double *c,double &x,double &y)
{
  // test of linear dependence
  if (MATH_ZERO < fabs(a[0]))
  {
    if (MATH_ZERO > fabs(a[1]*(a[0]+b[0])-a[0]*(a[1]+b[1])))
    {
      if (MATH_ZERO < fabs(a[1]*c[0]-a[0]*c[1]))
        return  0;
      else
        return -1;
    }
  }
  else
  {
    if (MATH_ZERO > fabs(b[1]*(a[0]+b[0])-b[0]*(a[1]+b[1])))
    {
      if (MATH_ZERO < fabs(b[1]*c[0]-b[0]*c[1]) || MATH_ZERO < fabs(b[1]*c[0]-b[0]*c[1]))
        return  0;
      else
      {
        if (MATH_ZERO > fabs(b[0]) && MATH_ZERO > fabs(b[1]) && MATH_ZERO < fabs(c[0]))
          return 0;
        else
          return -1;
      }
    }
  }
  
  if (MATH_ZERO < fabs(a[0]))
  {
    y = (c[0]*a[1]-c[1]*a[0]) / (a[0]*b[1]-a[1]*b[0]);
    x = (-c[0]-b[0]*y) / a[0];
  }
  else
  {
    y = -c[0]/b[0];
    x = (-c[1]-b[1]*y) / a[1];
  }
  
  return 1;
}



/**
   function solves system of two non-linear equations:
   a[0]*x*y + b[0]*x + c[0]*y + d[0] = 0
   a[1]*x*y + b[1]*x + c[1]*y + d[1] = 0
   answer: -1 = infinite number of results
           0,1,2 = number of results
   
   @param a,b,c,d - arrays of coeficients of equations, dimension is 2
   @param x,y - arrays results, dimension is 2
   
   created  4.4.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long solv_2nle (double *a, double *b, double *c, double *d,double *x,double *y)
{
  long ii,jj,nr;
  
  if (fabs(a[0]) < MATH_ZERO && fabs(a[1]) < MATH_ZERO)
    return (solv_2le (b,c,d,x[0],y[0]));
  
  if (fabs(a[0]) < MATH_ZERO || fabs(a[1]) < MATH_ZERO)
  {
    if (fabs(a[1]) > MATH_ZERO) 
    { 
      ii=0; jj=1; 
    }
    else                        
    { 
      ii=1; jj=0; 
    }
    if (fabs(b[ii]) < MATH_ZERO)
    {
      if (fabs(c[ii]) < MATH_ZERO)
      {
        if (fabs(d[ii]) < MATH_ZERO)
          return -1;
        else
          return  0;
      }
      else
      {
        y[0] = -d[ii]/c[ii];
        return solv_1le ( a[jj]*y[0]+b[jj] , c[jj]*y[0]+d[jj] ,x[0]);
      }
    }
    else
    {
      if (fabs(c[ii]) < MATH_ZERO)
      {
        x[0] = -d[ii]/b[ii];
        return solv_1le ( a[jj]*x[0]+c[jj] , b[jj]*x[0]+d[jj] ,y[0]);
      }
      else
      {
        nr = solv_polynom_2 ( -a[jj]*b[ii] , -a[jj]*d[ii] + c[ii]*b[jj]-c[jj]*b[ii] , c[ii]*d[jj]-c[jj]*d[ii] ,x[0],x[1]);
        if (nr>0)
          y[0] = (-b[ii]*x[0]-d[ii]) / c[ii];
        if (nr>1)
          y[1] = (-b[ii]*x[1]-d[ii]) / c[ii];
        return nr;
      }
    }
  }
  
  nr = solv_polynom_2 ( a[0]*b[1]-a[1]*b[0] , a[0]*d[1]-a[1]*d[0] + c[0]*b[1]-c[1]*b[0] , c[0]*d[1]-c[1]*d[0] ,x[0],x[1]);
  
  for (long i=0;i<nr;i++)
  {
    if (fabs(y[i] = a[0]*x[i]+c[0]) > MATH_ZERO)
      y[i] = (-b[0]*x[i]-d[0]) / y[i] ;
    else
    {
      if (fabs(b[0]*x[i]+d[0]) > MATH_ZERO)
        fprintf (stderr,"\n\n\n  ERROR 1 in function solv_3r, PLEASE CONTACT AUTHOR of this function(%s, line %d)\n\n\n",__FILE__,__LINE__);
      else
        return -1;
    }
  }
  return nr;
}



/**
   function computes natural coordinates of "point" in linear triangle
   
   @param xx,yy - cartesian coordinates of "point"
   @param x,y - arrays of cartesian coordinates of element nodes, dimension is 3
   @param xi,eta - answer - natural coordinates of "point"
   
   created  11.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void nc_lin_3_2d (double xx,double yy,double *x,double *y,double &xi,double &eta)
{
  double det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  xi  = ( (x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*xx + (x[2] - x[1])*yy )/det;
  eta = ( (x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*xx + (x[0] - x[2])*yy )/det;
}



/**
   function computes natural coordinates of "point" in quadratic triangle
   
   @param acc - relative accuracy of natural coordinates => (exact_xi - xi)/exact_xi is in (-acc;acc) ;;  acc > 1e-14
   @param xx,yy - cartesian coordinates of "point"
   @param x,y - arrays of cartesian coordinates of element nodes, dimension is 6
   @param xi,eta - answer - natural coordinates of "point"
   
   created  13.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void nc_quad_3_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta)
{
  long i,ne,nk;
  double aux,coeff[5],eeta[4],accabs;
  double ax,bx,cx,dx,ex,fx;
  double ay,by,cy,dy,ey,fy;
  double k,l,m,n,o,p,q,r,s;
  
  accabs = acc*( fabs(x[0]-x[1]) + fabs(x[1]-x[2]) + fabs(x[2]-x[0]) )/2.0;
  
  ax = 2.0*x[0] + 2.0*x[2] - 4.0*x[5];
  bx = 2.0*x[1] + 2.0*x[2] - 4.0*x[4];
  cx =    -x[0] - 3.0*x[2] + 4.0*x[5];
  dx =    -x[1] - 3.0*x[2] + 4.0*x[4];
  ex = 4.0*x[2] + 4.0*x[3] - 4.0*x[4] - 4.0*x[5];
  fx = x[2] - xx;
  
  ay = 2.0*y[0] + 2.0*y[2] - 4.0*y[5];
  by = 2.0*y[1] + 2.0*y[2] - 4.0*y[4];
  cy =    -y[0] - 3.0*y[2] + 4.0*y[5];
  dy =    -y[1] - 3.0*y[2] + 4.0*y[4];
  ey = 4.0*y[2] + 4.0*y[3] - 4.0*y[4] - 4.0*y[5];
  fy = y[2] - yy;
  
  //test whether triangle is quadratic
  if ( fabs(ax)+fabs(bx)+fabs(ex)+fabs(ay)+fabs(by)+fabs(ey) < accabs){
    nc_lin_3_2d (xx,yy,x,y,xi,eta);
    return;
  }
  
  if (fabs(ax) < accabs){
    coeff[0] = ay*bx*bx + by*ex*ex - ey*ex*bx ;
    coeff[1] = 2.0*ay*bx*dx + 2.0*by*cx*ex + dy*ex*ex - ey*cx*bx - cy*ex*bx - ey*ex*dx ;
    coeff[2] = 2.0*ay*bx*fx + 2.0*dy*cx*ex + ay*dx*dx + by*cx*cx + fy*ex*ex - cy*cx*bx - ey*cx*dx - cy*ex*dx - ey*ex*fx ;
    coeff[3] = 2.0*ay*fx*dx + 2.0*fy*cx*ex + dy*cx*cx - cy*cx*dx - ey*cx*fx - cy*ex*fx ;
    coeff[4] = ay*fx*fx + fy*cx*cx - cy*cx*fx ;
    
    ne = solv_polynom_4 (coeff,-acc,1+acc,1e-13,eeta);
    
    eta = 2.0;
    for (i=0;i<ne;i++)
    {
      xi = -( bx*eeta[i]*eeta[i] + dx*eeta[i] + fx ) / ( cx + ex*eeta[i] );
      if ( -acc < xi && xi < 1+acc )
      {
        if (eta==2.0)
          eta = eeta[i];
        else
          fprintf (stderr,"\n\n\n  error 0 in function nc_quad_3_2d (%s, line %d)\n\n\n",__FILE__,__LINE__);
      }
    }
    if (eta==2.0)
      fprintf (stderr,"\n\n\n  error 4 in function nc_quad_3_2d (%s, line %d)\n\n\n",__FILE__,__LINE__);
    
    return;
  }
  
  q = ex*ex-4.0*ax*bx ;
  r = 2.0*(cx*ex - 2.0*ax*dx) ;
  s = cx*cx - 4.0*ax*fx ;
  k = 4.0*by*ax*ax - 2.0*ey*ex*ax + q*ay + ex*ex*ay ;
  l = 4.0*dy*ax*ax - 2.0*ey*cx*ax + r*ay + 2.0*ex*(cx*ay - ax*cy) ;
  m = 4.0*fy*ax*ax - 2.0*cy*cx*ax + s*ay + cx*cx*ay ;
  n =       ( 2.0*ax*ey - ex*(ay + ay) ) * ( 2.0*ax*ey - ex*(ay + ay) ) ;
  o = 2.0 * ( 2.0*ax*ey - ex*(ay + ay) ) * ( 2.0*ax*cy - ay*(cx + cx) ) ;
  p =       ( 2.0*ax*cy - ay*(cx + cx) ) * ( 2.0*ax*cy - ay*(cx + cx) ) ;
  
  coeff[0] = k*k - n*q ;
  coeff[1] = 2.0*k*l - n*r - o*q ;
  coeff[2] = 2.0*m*k + l*l - n*s - o*r - p*q ;
  coeff[3] = 2.0*l*m - o*s - p*r ;
  coeff[4] = m*m - p*s ;
  
  ne = solv_polynom_4 (coeff,-acc,1+acc,1e-13,eeta);
  
  nk = 0;
  eta = 2.0;
  for (i=0;i<ne;i++){
    nk = solv_polynom_2 ( ax , cx+ex*eeta[i] , bx*eeta[i]*eeta[i]+dx*eeta[i]+fx , xi , aux );
    
    if ( nk>0 && -acc < xi && xi < 1+acc )
    {
      if (eta==2.0)
        eta = eeta[i];
      else
        fprintf (stderr,"\n\n\n  error 1 in function nc_quad_3_2d (%s, line %d)\n\n\n",__FILE__,__LINE__);
    }
    if ( nk>1 && -acc < aux && aux < 1+acc )
    {
      if (eta==2.0)
      {
        xi = aux;
        eta = eeta[i];
      }
      else
        fprintf (stderr,"\n\n\n  error 2 in function nc_quad_3_2d (%s, line %d)\n\n\n",__FILE__,__LINE__);
    }
  }
  if (eta==2.0)
    fprintf (stderr,"\n\n\n  error 3 in function nc_quad_3_2d (%s, line %d)\n\n\n",__FILE__,__LINE__);
  
}



/**
   function computes natural coordinates of "point" in linear rectangle
   
   @param acc - relative accuracy of natural coordinates => (exact_xi - xi)/exact_xi is in (-acc;acc) ;;  acc > 1e-14
   @param xx,yy - cartesian coordinates of "point"
   @param x,y - arrays of cartesian coordinates of element nodes, dimension is 4
   @param xi,eta - answer - natural coordinates of "point"
   
   created  11.1.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void nc_lin_4_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta)
{
  long nroots;
  double ax,bx,cx,dx,ay,by,cy,dy;
  double a,b,c,ksi1,ksi2,accabs;
  
  accabs = acc*( fabs(x[0]-x[1]) + fabs(x[1]-x[2]) + fabs(x[2]-x[3]) + fabs(x[3]-x[0]) )/2.0;
  
  ax = x[0] - x[1] + x[2] - x[3];
  bx = x[0] - x[1] - x[2] + x[3];
  cx = x[0] + x[1] - x[2] - x[3];
  dx = x[0] + x[1] + x[2] + x[3] - 4.0 * xx;
  
  ay = y[0] - y[1] + y[2] - y[3];
  by = y[0] - y[1] - y[2] + y[3];
  cy = y[0] + y[1] - y[2] - y[3];
  dy = y[0] + y[1] + y[2] + y[3] - 4.0 * yy;
  
  a = ax * by - ay * bx;
  b = ax * dy - ay * dx + by * cx - bx * cy;
  c = cx * dy - cy * dx;
  
  nroots = solv_polynom_2 (a,b,c,ksi1,ksi2);
  
  if (nroots == 0){
    print_err("no real roots in the solution of 2nd order polynom", __FILE__, __LINE__, __func__);
    abort();
  }
  if (nroots == 1)  xi = ksi1;
  if (nroots == 2)
  {
    if (( -1-acc <= ksi1  &&  ksi1 <= 1+acc ) && ( ksi2 < -1-acc  ||  1+acc < ksi2 ))
      xi = ksi1;
    else 
    {
      if (( -1-acc <= ksi2  &&  ksi2 <= 1+acc ) && ( ksi1 < -1-acc  ||  1+acc < ksi1 ))
        xi = ksi2;
      else{
        print_err("both roots of polynom are in the range [-1;1], cannot determine ksi", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
  
  if ( fabs(cx+ax*xi) < accabs )
  {
    if ( fabs(bx) < accabs ){
      print_err ("zero bx and denominator in eta determination", __FILE__, __LINE__, __func__);
      abort();
    }
    else{
      if ( fabs(ay*xi+cy) < accabs || fabs(xi+dx/bx) > acc ){
        print_err ("zero denominator in eta determination",__FILE__, __LINE__, __func__);
        abort();
      }
      else
        eta = (-dy-by*xi)/(cy+ay*xi);
    }
  }
  else
    eta = (-dx-bx*xi)/(cx+ax*xi);
  
}



/**
   function computes natural coordinates of "point" in quadratic rectangle
   
   @param acc - relative accuracy of natural coordinates => (exact_xi - xi)/exact_xi is in (-acc;acc) ;;  acc > 1e-14
   @param xx,yy - cartesian coordinates of "point"
   @param x,y - arrays of cartesian coordinates of element nodes, dimension is 8
   @param xi,eta - answer - natural coordinates of "point"
   
   created  xx.x.xxxx, xxxxxxxx xxxxxxx, xxxxxx@cml.fsv.cvut.cz
*/
void nc_quad_4_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta)
{
  nc_lin_4_2d (acc,xx,yy,x,y,xi,eta);
  /*
    staci vyresit tuto soustavu rovnic :-)

    a*x*x*y + b*x*y*y + c*x*x + d*y*y + e*x*y + f*x + g*y + i = 0
    k*x*x*y + l*x*y*y + m*x*x + n*y*y + o*x*y + p*x + q*y + r = 0
   */
}

/**
   Function (i) divides all "elements" into linear "subelements",
   because postprocesors are able to image only linear - line, triangle, tetrahedra, rectangle and cube.
   (ii) fills array 'lin' by lines (between two nodes), of which is composed mesh.
   (iii) computes coordinates of centre point for quadratic rectangle elements.
   -> linear triangle is divided into one lineat triangle
   -> linear rectangle is divided into one lineat rectangle
   -> quadratic triangle is divided into four lineat triangles and composed of 6 lines
   -> quadratic rectangle is divided into four lineat rectangles and composed of 8 lines
   -> linear tetrahedra is divided into one linear tetrahedra
   -> linear hexahedra(cube) is divided into one linear cube

   @param nli - number of lines
   @param lin - array of line nodes
   @param ntr - number of triangles
   @param trn - array of triangle nodes
   @param nqu - number of quadrilaterals (rectangles)
   @param qun - array of quadrilateral nodes
   @param nte - number of tetrahedras
   @param ten - array of tetrahedra nodes
   @param ncu - number of cubes
   @param cun - array of cube nodes
   @param nmp - number of midpoints
   @param mpn - array of midpoint coordinates and numbers of corresponding elements
   @param flag - switch: 'e' = data for postprocesor Elixir, 'e' = data for postprocesor OpenDx
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/

void interpolelem (gtopology *gt,long &nli,long **&lin,long &ntr,long **&trn,long &nqu,long **&qun,long &nte,long **&ten,long &ncu,long **&cun,long &nmp,double **&mpn,char flag)
{
  long i;
  long *nodes;
  
  nli = ntr = nqu = nte = ncu = nmp = 0;
  
  lin = new long* [gt->ne*8];
  trn = new long* [gt->ne*4];
  qun = new long* [gt->ne*4];
  ten = new long* [gt->ne*1];
  cun = new long* [gt->ne*1];
  mpn = new double* [gt->ne];
  
  for (i=0;i<gt->ne;i++){
    nodes = gt->gelements[i].nodes;
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      trn[ntr] = new long[3];  trn[ntr][0] = nodes[0];  trn[ntr][1] = nodes[1];  trn[ntr][2] = nodes[2];  ntr++;
      break;
    }
    case 412:{
      qun[nqu] = new long [4];
      if (flag == 'e')  { qun[nqu][0] = nodes[0]; qun[nqu][1] = nodes[1]; qun[nqu][2] = nodes[2]; qun[nqu][3] = nodes[3];  nqu++; }
      if (flag == 'd')  { qun[nqu][0] = nodes[0]; qun[nqu][1] = nodes[3]; qun[nqu][2] = nodes[1]; qun[nqu][3] = nodes[2];  nqu++; }
      break;
    }
    case 622:{
      trn[ntr] = new long [3];  trn[ntr][0] = nodes[0];  trn[ntr][1] = nodes[3];  trn[ntr][2] = nodes[5];   ntr++;
      trn[ntr] = new long [3];  trn[ntr][0] = nodes[3];  trn[ntr][1] = nodes[1];  trn[ntr][2] = nodes[4];   ntr++;
      trn[ntr] = new long [3];  trn[ntr][0] = nodes[4];  trn[ntr][1] = nodes[2];  trn[ntr][2] = nodes[5];   ntr++;
      trn[ntr] = new long [3];  trn[ntr][0] = nodes[3];  trn[ntr][1] = nodes[4];  trn[ntr][2] = nodes[5];   ntr++;
      if (flag == 'd'){
        lin[nli] = new long [2];  lin[nli][0] = nodes[0];  lin[nli][1] = nodes[3];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[3];  lin[nli][1] = nodes[1];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[1];  lin[nli][1] = nodes[4];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[4];  lin[nli][1] = nodes[2];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[2];  lin[nli][1] = nodes[5];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[5];  lin[nli][1] = nodes[0];   nli++;
      }
      break;
    }
    case 822:{
      if (flag == 'e'){
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[7]; qun[nqu][1] = nodes[0]; qun[nqu][2] = nodes[4]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[4]; qun[nqu][1] = nodes[1]; qun[nqu][2] = nodes[5]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[5]; qun[nqu][1] = nodes[2]; qun[nqu][2] = nodes[6]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[6]; qun[nqu][1] = nodes[3]; qun[nqu][2] = nodes[7]; qun[nqu][3] = gt->nn+nmp;  nqu++;
      }
      if (flag == 'd'){
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[0]; qun[nqu][1] = nodes[7]; qun[nqu][2] = nodes[4]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[1]; qun[nqu][1] = nodes[4]; qun[nqu][2] = nodes[5]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[2]; qun[nqu][1] = nodes[5]; qun[nqu][2] = nodes[6]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        qun[nqu] = new long [4];  qun[nqu][0] = nodes[3]; qun[nqu][1] = nodes[6]; qun[nqu][2] = nodes[7]; qun[nqu][3] = gt->nn+nmp;  nqu++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[0];  lin[nli][1] = nodes[4];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[4];  lin[nli][1] = nodes[1];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[1];  lin[nli][1] = nodes[5];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[5];  lin[nli][1] = nodes[2];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[2];  lin[nli][1] = nodes[6];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[6];  lin[nli][1] = nodes[3];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[3];  lin[nli][1] = nodes[7];   nli++;
        lin[nli] = new long [2];  lin[nli][0] = nodes[7];  lin[nli][1] = nodes[0];   nli++;
      }
      mpn[nmp] = new double [3];
      mpn[nmp][0] = (gt->gnodes[nodes[4]].x + gt->gnodes[nodes[5]].x + gt->gnodes[nodes[6]].x + gt->gnodes[nodes[7]].x) / 4.0 ;
      mpn[nmp][1] = (gt->gnodes[nodes[4]].y + gt->gnodes[nodes[5]].y + gt->gnodes[nodes[6]].y + gt->gnodes[nodes[7]].y) / 4.0 ;
      mpn[nmp][2] = i;
      nmp++;
      break;
    }
    case 413:{
      if (flag == 'd'){
        ten[nte] = new long[4];  ten[nte][0] = nodes[0];  ten[nte][1] = nodes[1];  ten[nte][2] = nodes[2];  ten[nte][3] = nodes[3];  nte++;
      }
      if (flag == 'e')
        fprintf (stderr,"\n\n function interpolelem does not transform tetrahedras for drawing in elixir \n\n");
      break;
    }
    case 813:{
      if (flag == 'd'){
        cun[ncu] = new long[8];  cun[ncu][0] = nodes[0];  cun[ncu][1] = nodes[1];  cun[ncu][2] = nodes[3];  cun[ncu][3] = nodes[2];
                                 cun[ncu][4] = nodes[4];  cun[ncu][5] = nodes[5];  cun[ncu][6] = nodes[7];  cun[ncu][7] = nodes[6];  ncu++;
      }
      if (flag == 'e')
        fprintf (stderr,"\n\n function interpolelem does not transform hexahedras for drawing in elixir \n\n");
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong nnedegdim in function");
      fprintf (stderr," interpolelem (%s, line %d)",__FILE__,__LINE__);
    }
    }
  }
}

/**
   function prints data-file for Elixir
   it is possible only one
   
   @param gt - gtopology
   @param file - name of printed file
   @param valnod - array of depicted values in nodes
   @param valel - array of depicted values in centre of element - only for quadratic rectangel, for other elements **valel==anything
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void print_ex (gtopology *gt,const char *file, double *valnod, double *valel)
{
  long i,nn;
  long nli,ntr,nqu,ntet,ncub,nmp;
  long **linodes,**trnodes,**qunodes,**tetnodes,**cubnodes;
  double **mpnodes;
  FILE *data_ex;
  
  nn = gt->nn;
  
  interpolelem (gt,nli,linodes,ntr,trnodes,nqu,qunodes,ntet,tetnodes,ncub,cubnodes,nmp,mpnodes,'e');
  
  data_ex = fopen (file,"w");  if (data_ex==NULL) fprintf (stderr,"\n .ex file has not been opened.\n\n");
  
  // ** heading **
  fprintf (data_ex,"%10d%10d\n",7,1);
  fprintf (data_ex,"%10ld%10d%10ld%10ld%10d%10ld%10d%10d\n\n",nn+nmp,0,ntr,nqu,0,ntet,0,0);
  
  // ** nodes **
  for (i=0;i<nn;i++)
    fprintf (data_ex,"%10ld  %.6e  %.6e  %.6e  %.6e\n",i+1,gt->gnodes[i].x,gt->gnodes[i].y,gt->gnodes[i].z,valnod[i]);
  
  if (nmp)
    for (i=0;i<nmp;i++)
      fprintf (data_ex,"%10ld  %.6e  %.6e  %.6e  %.6e\n",nn+i+1,mpnodes[i][0],mpnodes[i][1],0.0,valel[(long)mpnodes[i][2]]);
  
  fprintf (data_ex,"\n");
  
  // ** elements **
  for (i=0;i<ntr;i++)
    fprintf (data_ex,"%10ld%10ld%10ld%10ld\n",i+1,trnodes[i][0]+1,trnodes[i][1]+1,trnodes[i][2]+1);

  for (i=0;i<nqu;i++)
    fprintf (data_ex,"%10ld%10ld%10ld%10ld%10ld\n",i+1,qunodes[i][0]+1,qunodes[i][1]+1,qunodes[i][2]+1,qunodes[i][3]+1);
      
  for (i=0;i<ntet;i++)
    fprintf (data_ex,"%10ld%10ld%10ld%10ld%10ld\n",i+1,tetnodes[i][0]+1,tetnodes[i][1]+1,tetnodes[i][2]+1,tetnodes[i][3]+1);
  
  
  fclose (data_ex);
  
  for (i=0;i<nli;i++) delete [] linodes[i];
  delete [] linodes;
  for (i=0;i<ntr;i++) delete [] trnodes[i];
  delete [] trnodes;
  for (i=0;i<nqu;i++) delete [] qunodes[i];
  delete [] qunodes;
  for (i=0;i<ntet;i++) delete [] tetnodes[i];
  delete [] tetnodes;
  for (i=0;i<ncub;i++) delete [] cubnodes[i];
  delete [] cubnodes;
  for (i=0;i<nmp;i++) delete [] mpnodes[i];
  delete [] mpnodes;
}

/**
   function prints data-file for OpenDX
   
   @param gt - gtopology
   @param file - name of printed file
   @param valnod - array of depicted values in nodes
   @param valel - array of depicted values in centre of element - only for quadratic rectangel, for other elements **valel==anything
   @param dimindex - 
   @param caption  - 
   @param nindex - dimension of dimindex
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void print_dx (gtopology *gt, const char *file,double **valnod,double **valel,char tve, long *dimindex, char **caption,long nindex)
{
  long i,j,k,nn,dim;
  long nli,ntr,nqu,ntet,ncub,nmp,object,value,vali;
  long **linodes,**trnodes,**qunodes,**tetnodes,**cubnodes;
  double **mpnodes;
  FILE *data_dx;
  
  nn = gt->nn;
  dim = gt->gelements[0].auxinf%10;
  object = 0;
  
  interpolelem (gt,nli,linodes,ntr,trnodes,nqu,qunodes,ntet,tetnodes,ncub,cubnodes,nmp,mpnodes,'d');
  
  data_dx = fopen (file,"w");  if (data_dx==NULL) fprintf (stderr,"\n .dx file has not been opened.\n\n");
  
  // ** nodecoords **
  fprintf (data_dx,"object %ld class array type float rank 1 shape %ld items %ld data follows\n",++object,dim,nn+nmp);
  
  for (i=0;i<nn;i++){
    fprintf (data_dx," %20.10f",gt->gnodes[i].x);
    if (dim==3 || dim==2) 
      fprintf (data_dx," %20.10f",gt->gnodes[i].y);
    if (dim==3)
      fprintf (data_dx," %20.10f",gt->gnodes[i].z);
    fprintf (data_dx,"\n");
  }
  
  if (nmp)
    for (i=0;i<nmp;i++)
      fprintf (data_dx," %20.10f %20.10f\n",mpnodes[i][0],mpnodes[i][1]);
  
  // ** elements **
  if (ntr){
    fprintf (data_dx,"object %ld class array type int rank 1 shape %d items %ld data follows\n",++object,3,ntr);
    
    for (i=0;i<ntr;i++){
      for (j=0;j<3;j++)
        fprintf (data_dx," %10ld",trnodes[i][j]);
      
      fprintf (data_dx,"\n");
    }
    
    fprintf (data_dx,"attribute \"element type\" string \"%s\"\n","triangles");
    fprintf (data_dx,"attribute \"ref\" string \"positions\"\n");
  }
  
  if (nqu){
    fprintf (data_dx,"object %ld class array type int rank 1 shape %d items %ld data follows\n",++object,4,nqu);
    
    for (i=0;i<nqu;i++){
      for (j=0;j<4;j++)
        fprintf (data_dx," %10ld",qunodes[i][j]);
      
      fprintf (data_dx,"\n");
    }
    
    fprintf (data_dx,"attribute \"element type\" string \"%s\"\n","quads");
    fprintf (data_dx,"attribute \"ref\" string \"positions\"\n");
  }
  
  if (ntet){
    fprintf (data_dx,"object %ld class array type int rank 1 shape %d items %ld data follows\n",++object,4,ntet);
    
    for (i=0;i<ntet;i++){
      for (j=0;j<4;j++)
        fprintf (data_dx," %10ld",tetnodes[i][j]);
      
      fprintf (data_dx,"\n");
    }
    
    fprintf (data_dx,"attribute \"element type\" string \"%s\"\n","tetrahedra");
    fprintf (data_dx,"attribute \"ref\" string \"positions\"\n");
  }
  
  if (ncub){
    fprintf (data_dx,"object %ld class array type int rank 1 shape %d items %ld data follows\n",++object,8,ncub);
    
    for (i=0;i<ncub;i++){
      for (j=0;j<8;j++)
        fprintf (data_dx," %10ld",cubnodes[i][j]);
      
      fprintf (data_dx,"\n");
    }
    
    fprintf (data_dx,"attribute \"element type\" string \"%s\"\n","cubes");
    fprintf (data_dx,"attribute \"ref\" string \"positions\"\n");
  }
  
  if (nli){
    fprintf (data_dx,"object %ld class array type int rank 1 shape %d items %ld data follows\n",++object,2,nli);
    
    for (i=0;i<nli;i++){
      for (j=0;j<2;j++)
        fprintf (data_dx," %10ld",linodes[i][j]);
      
      fprintf (data_dx,"\n");
    }
    
    fprintf (data_dx,"attribute \"element type\" string \"%s\"\n","lines");
    fprintf (data_dx,"attribute \"ref\" string \"positions\"\n");
  }
  
  
  // ** valnod **
  vali = 0;
  for (i=0;i<nindex;i++){
    if (dimindex[i]==1)
      fprintf (data_dx,"object %ld class array type float rank 0 items %ld data follows\n",++object,nn+nmp);
    else
      fprintf (data_dx,"object %ld class array type float rank 1 shape %ld items %ld data follows\n",++object,dimindex[i],nn+nmp);
    
    for (j=0;j<nn;j++){
      for (k=0;k<dimindex[i];k++)
        fprintf (data_dx," %20.10f",valnod[vali+k][j]);
      
      fprintf (data_dx,"\n");
    }
    
    if (nmp)
      for (j=0;j<nmp;j++){
        for (k=0;k<dimindex[i];k++)
          if      (tve=='n') fprintf (data_dx," %20.10f",valel[vali+k][(long)mpnodes[j][2]]);
          else if (tve=='t') fprintf (data_dx," %20.10f",valel[(long)mpnodes[j][2]][vali+k]);
          else fprintf (stderr,"\n\n\n error in print_dx in mathem.cpp \n\n\n");
        
        fprintf (data_dx,"\n");
      }
    
    
    fprintf (data_dx,"attribute \"dep\" string \"positions\"\n");
    vali += dimindex[i];
  }
  
  // ** fields **
  value = 1;
  if (ntr) value++;
  if (nqu) value++;
  if (ntet) value++;
  if (ncub) value++;
  if (nli) value++;
  for (i=0;i<nindex;i++){
    fprintf (data_dx,"object %ld class field\n",++object);
    fprintf (data_dx,"component \"positions\" value 1\n");
    fprintf (data_dx,"component \"connections\" value 2\n");
    fprintf (data_dx,"component \"data\" value %ld\n",++value);
    if (ntr && nqu){
      fprintf (data_dx,"object %ld class field\n",++object);
      fprintf (data_dx,"component \"positions\" value 1\n");
      fprintf (data_dx,"component \"connections\" value 3\n");
      fprintf (data_dx,"component \"data\" value %ld\n",value);
    }
  }

  // ** compositefields **
  if (ntr && nqu){
    for (i=0;i<nindex;i++){
      fprintf (data_dx,"object %ld class compositefield\n",++object);
      fprintf (data_dx,"member 0 value %ld\n",++value);
      fprintf (data_dx,"member 1 value %ld\n",++value);
    }
  }

  // ** lincon field **
  if (nli){
    fprintf (data_dx,"object %ld class field\n",++object);
    fprintf (data_dx,"component \"positions\" value 1\n");
    if (nqu && ntr){
      fprintf (data_dx,"component \"connections\" value 4\n");
      fprintf (data_dx,"component \"data\" value 5\n");
    }
    else{
      fprintf (data_dx,"component \"connections\" value 3\n");
      fprintf (data_dx,"component \"data\" value 4\n");
    }
  }

  // ** groups **
  fprintf (data_dx,"object %ld class group\n",++object);

  for (i=0;i<nindex;i++)
    fprintf (data_dx,"member \"%s\"  value %ld\n",caption[i],++value);

  if (nli)
    fprintf (data_dx,"member \"%s\"  value %ld\n","lincon",++value);

  fclose (data_dx);

  for (i=0;i<nli;i++) delete [] linodes[i];
  delete [] linodes;
  for (i=0;i<ntr;i++) delete [] trnodes[i];
  delete [] trnodes;
  for (i=0;i<nqu;i++) delete [] qunodes[i];
  delete [] qunodes;
  for (i=0;i<ntet;i++) delete [] tetnodes[i];
  delete [] tetnodes;
  for (i=0;i<ncub;i++) delete [] cubnodes[i];
  delete [] cubnodes;
  for (i=0;i<nmp;i++) delete [] mpnodes[i];
  delete [] mpnodes;
}

/**
   Function finds maximal and minimal number in first three members of array x.
   
   @param x - 1D array of given numbers
   @param max,min - answers
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void maxmin_3 (double *x,double &max,double &min)
{
  if (x[0] > x[1])   if (x[0] > x[2])   if (x[1] > x[2]){ max = x[0]; min = x[2];}
                                        else            { max = x[0]; min = x[1];}
                     else                               { max = x[2]; min = x[1];}
  else               if (x[1] > x[2])   if (x[0] > x[2]){ max = x[1]; min = x[2];}
                                        else            { max = x[1]; min = x[0];}
                     else                               { max = x[2]; min = x[0];}
}

/**
   Function finds maximal and minimal number in first four members of array x.
   
   @param x - 1D array of given numbers
   @param max,min - answers
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void maxmin_4 (double *x,double &max,double &min)
{
  if (x[0] > x[1]) if (x[0] > x[2]) if (x[0] > x[3]) if (x[1] < x[2]) if (x[1] < x[3]){ max = x[0]; min = x[1];}
                                                                      else            { max = x[0]; min = x[3];}
                                                     else             if (x[2] < x[3]){ max = x[0]; min = x[2];}
                                                                      else            { max = x[0]; min = x[3];}
                                    else                              if (x[1] < x[2]){ max = x[3]; min = x[1];}
                                                                      else            { max = x[3]; min = x[2];}
                   else                              if (x[2] > x[3]) if (x[1] < x[3]){ max = x[2]; min = x[1];}
                                                                      else            { max = x[2]; min = x[3];}
                                                     else                             { max = x[3]; min = x[1];}

  else             if (x[1] > x[2]) if (x[1] > x[3]) if (x[0] < x[2]) if (x[0] < x[3]){ max = x[1]; min = x[0];}
                                                                      else            { max = x[1]; min = x[3];}
                                                     else             if (x[2] < x[3]){ max = x[1]; min = x[2];}
                                                                      else            { max = x[1]; min = x[3];}
                                    else                              if (x[0] < x[2]){ max = x[3]; min = x[0];}
                                                                      else            { max = x[3]; min = x[2];}
                   else                              if (x[2] > x[3]) if (x[0] < x[3]){ max = x[2]; min = x[0];}
                                                                      else            { max = x[2]; min = x[3];}
                                                     else                             { max = x[3]; min = x[0];}
}



/**
  Function assembles matrix for point approximation with quadratic function
  y = a*x^2 + b*x + c. The least square method is used. 
  One new point is added and coefficients of the matrix a and vector of right hand side r
  are updated and solved. One old point is subtracted from the matrix (it is used for moving 
  least square method). The solution is saved in the vector l and the function returns values of 
  derivatives in the point x. The function is used in iterative solvers (e.g. Newton-Raphson) 
  for divergency check.
  
  @param a - matrix of normal equations
  @param r - right hand side vector
  @param l - vector of solution(coefficients of quadratic function)
             l[0] = c, l[1] = b, l[2] = a.
  @param x - x coordinate of added point
  @param y - y coordinate of added point
  @param x_old - x coordinate of old point which will be subtracted
  @param y_old - y coordinate of old point which will be subtracted
  @param zero - limit value for zero determinant detection
  @param solc - flag for computation of coefficients of quadratic approximation
                solc=0  accumulation of values only, no computation of coefficients, 
                solc=1  accumulation of values and computation of coeffcients and derivatives
                solc=2  accumulation/removal of values and computation of coeffcients and derivatives

  @retval The function returns value of the derivative in the given point [x,y] (2*a*x+b)

  TKo 8.2008
*/
double lsm_quad(matrix &a, vector &r, vector &l, double x, double y, double x_old, double y_old, double zero, long solc)
{
  double det_a, det_b1, det_b2, det_b3;

  if ((a.m != 3) && (a.n != 3) && (r.n != 3) && (l.n != 3))
  {
    print_err("Wrong dimension of matrix or vectors in lest square method", __FILE__, __LINE__, __func__);
    abort();
  }
  // accumulation of values into matrix of normal equations
  a[0][0] += 1.0;  a[0][1] += x;      a[0][2] += x*x;
  r[0] += y;
  a[1][0] += x;    a[1][1] += x*x;    a[1][2] += x*x*x;
  r[1] += x*y;
  a[2][0] += x*x;  a[2][1] += x*x*x;  a[2][2] += x*x*x*x;
  r[2] += x*x*y;
  
  if (solc == 0)
    return 0.0;

  if (solc == 2)
  {
    // removing old values from matrix of normal equations
    a[0][0] -= 1.0;  a[0][1] -= x_old;      a[0][2] -= x_old*x_old;
    r[0] -= y_old;
    a[1][0] -= x_old;    a[1][1] -= x_old*x_old;    a[1][2] -= x_old*x_old*x_old;
    r[1] -= x_old*y_old;
    a[2][0] -= x_old*x_old;  a[2][1] -= x_old*x_old*x_old;  a[2][2] -= x_old*x_old*x_old*x_old;
    r[2] -= x_old*x_old*y_old;
  }

  // solution of system of equatuions by the Cramers rule
  det_a  = a[0][0]*a[1][1]*a[2][2]; 
  det_a += a[0][1]*a[1][2]*a[2][0];
  det_a += a[0][2]*a[1][0]*a[2][1];
  det_a -= a[0][2]*a[1][1]*a[2][0];
  det_a -= a[0][0]*a[1][2]*a[2][1];
  det_a -= a[0][1]*a[1][0]*a[2][2];

  if (fabs(det_a) < zero)
  {
    print_err("Determinant in least square method is close to zero", __FILE__, __LINE__, __func__);
    return 0.0;
  }

  det_b1  = r[0]*a[1][1]*a[2][2]; 
  det_b1 += a[0][1]*a[1][2]*r[2];
  det_b1 += a[0][2]*r[1]*a[2][1];
  det_b1 -= a[0][2]*a[1][1]*r[2];
  det_b1 -= r[0]*a[1][2]*a[2][1];
  det_b1 -= a[0][1]*r[1]*a[2][2];

  det_b2  = a[0][0]*r[1]*a[2][2]; 
  det_b2 += r[0]*a[1][2]*a[2][0];
  det_b2 += a[0][2]*a[1][0]*r[2];
  det_b2 -= a[0][2]*r[1]*a[2][0];
  det_b2 -= a[0][0]*a[1][2]*r[2];
  det_b2 -= r[0]*a[1][0]*a[2][2];

  det_b3  = a[0][0]*a[1][1]*r[2]; 
  det_b3 += a[0][1]*r[1]*a[2][0];
  det_b3 += r[0]*a[1][0]*a[2][1];
  det_b3 -= r[0]*a[1][1]*a[2][0];
  det_b3 -= a[0][0]*r[1]*a[2][1];
  det_b3 -= a[0][1]*a[1][0]*r[2];

  l[0] = det_b1/det_a;
  l[1] = det_b2/det_a;
  l[2] = det_b3/det_a; 

  return 2*l[2]*x + l[1];
}


int check_math_err()
{
  switch (errno)
  {
    case EDOM:
      print_err("math domain error (nan) detected\n", __FILE__, __LINE__, __func__);
      abort();
    case ERANGE:
      print_err("math range error (inf) detected\n", __FILE__, __LINE__, __func__);
      abort();
    default:
      if (errno > 0)
      {
        print_err("%s detected\n", __FILE__, __LINE__, __func__, strerror(errno));
        abort();
      }
  }
  return errno;
}



int check_math_errel(long eid)
{
  switch (errno)
  {
    case EDOM:
      print_err("math domain error (nan) detected on element %ld\n", __FILE__, __LINE__, __func__, eid+1);
      //abort();//here, only for debug??!!
      break;
    case ERANGE:
      print_err("math range error (inf) detected on element %ld\n", __FILE__, __LINE__, __func__, eid+1);
      break;
    default:
      if (errno > 0){
        print_err("%s detected on element %ld\n", __FILE__, __LINE__, __func__, strerror(errno), eid+1);
	exit(0);
      }
  }
  return errno;
}



int check_math_errip(long eid, long ipp)
{
  switch (errno)
  {
    case EDOM:
      print_err("math domain error (nan) detected on element %ld (ipp=%ld)\n", __FILE__, __LINE__, __func__, eid+1, ipp);
      abort();
      break;
    case ERANGE:
      print_err("math range error (inf) detected on element %ld (ipp=%ld)\n", __FILE__, __LINE__, __func__, eid+1, ipp);
      abort();
      break;
    default:
      if (errno > 0){
        print_err("%s detected on element %ld  (ipp=%ld)\n", __FILE__, __LINE__, __func__, strerror(errno), eid+1, ipp);
	exit(0);
      }
  }
  return errno;
}



long test_math_err()
{
  if ((errno == EDOM) || (errno == ERANGE))
    return 1;

  if (errno > 0)
    return 2;

  return 0;
}


int test_math_errip(long eid, long ipp)
{
  switch (errno)
  {
    case EDOM:
      print_err("math domain error (nan) detected on element %ld (ipp=%ld)\n", __FILE__, __LINE__, __func__, eid+1, ipp);
      break;
    case ERANGE:
      print_err("math range error (inf) detected on element %ld (ipp=%ld)\n", __FILE__, __LINE__, __func__, eid+1, ipp);
      break;
    default:
      if (errno > 0){
        print_err("%s detected on element %ld  (ipp=%ld)\n", __FILE__, __LINE__, __func__, strerror(errno), eid+1, ipp);
      }
  }
  return errno;
}



int test_math_errel(long eid)
{
  switch (errno)
  {
    case EDOM:
      print_err("math domain error (nan) detected on element %ld\n", __FILE__, __LINE__, __func__, eid+1);
      //abort();//here, only for debug??!!
      break;
    case ERANGE:
      print_err("math range error (inf) detected on element %ld\n", __FILE__, __LINE__, __func__, eid+1);
      break;
    default:
      if (errno > 0){
        print_err("%s detected on element %ld\n", __FILE__, __LINE__, __func__, strerror(errno), eid+1);
      }
  }
  return errno;
}



