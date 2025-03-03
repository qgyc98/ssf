/*
 * mtk - Maths Toolkit for X11
 *
 * Copyright 1994-1997   andrewr@chiark.greenend.org.uk (Andrew Ross)
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 ********/

#include <math.h>
#include <errno.h>
#include <float.h>
#include "moremath.h"

#ifdef ACC
 #define ACCOLD ACC
 #define ACC 40
#else
 #undef ACCOLD
 #define ACC 40
#endif

#ifdef BIGNO
 #define BIGNOOLD BIGNO
 #define BIGNO 1.0e10
#else
 #undef BIGNOOLD
 #define BIGNO 1.0e10
#endif

#ifdef BIGNI
 #define BIGNIOLD BIGNI
 #define BIGNI 1.0e-10
#else
 #undef BIGNIOLD
 #define BIGNI 1.0e-10
#endif

#ifdef ITMAX
 #define ITMAXOLD ITMAX
 #define ITMAX 100
#else
 #undef ITMAXOLD
 #define ITMAX 100
#endif

#ifdef EPS
 #define EPSOLD EPS
 #define EPS 3.0e-7
#else
 #undef EPSOLD
 #define EPS 3.0e-7
#endif

#ifdef FPMIN
 #define FPMINOLD FPMIN
 #define FPMIN 1.0e-30
#else
 #undef FPMINOLD
 #define FPMIN 1.0e-30
#endif

/* Return asinh */
double ansi_asinh(double arg)
{
    arg = log(arg+sqrt(arg*arg+1));
    return arg;
}

/* Return acosh */
double ansi_acosh(double arg)
{
    if (arg < 1)
    {
       errno = EDOM;
       return 0;
    };
    arg = log(arg+sqrt(arg*arg-1));
    return arg;
}

/* Return atanh */
double ansi_atanh(double arg)
{
    if ((arg <= -1) || (arg >= 1))
    {
       errno = EDOM;
       return 0;
    };
    arg = 0.5*log(1+arg)-0.5*log(1-arg);
    return arg;
}

/* Return log of gamma function */
double ansi_lgamma(double arg)
{
   double x,y,ser,tmp;
   static double coeff[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,
	-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   int j;

   y = x = arg;
   tmp = x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser = 1.000000000190015;
   for (j = 0; j <= 5;j++) ser += coeff[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}

/* Return Bessel function J0 */

double ansi_j0(double arg)
{
  double absarg,z,x,y,ans,ans1,ans2;

  if ( (absarg = fabs(arg)) < 8 ){
    y = arg*arg;
    ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7 + y*(-11214424.18 + y*(77392.33017 - y*184.9052456))));
    ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718 + y*(59272.64853 + y*(267.8532712+y))));
    ans = ans1/ans2;
  }
  else {
    z = 8/absarg;
    y = z*z;
    x = absarg - 0.785398164;
    ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4 + y*(-0.2073370639e-5 + y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1 + y*(0.1430488765e-3 + y*(-0.6911147651e-5 + y*(0.7621095161e-6 - y*0.934935152e-7)));
    ans = sqrt( 0.636619772 / absarg ) * (cos(x)*ans1 - z*sin(x)*ans2);
  }
  return ans;
}

/* Return Bessel function J1 */

double ansi_j1(double arg)
{
  double absarg,z,x,y,ans,ans1,ans2;

  if ( (absarg = fabs(arg)) < 8 ){
    y = arg*arg;
    ans1 = arg*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1 + y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
    ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74 + y*(99447.43394 + y*(376.9991397 + y))));
    ans = ans1/ans2;
  }
  else {
    z = 8/absarg;
    y = z*z;
    x = absarg - 2.356194491;
    ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4 + y*(0.2457520174e-5 + y*(-0.240337019e-6))));
    ans2 = 0.04687499995 + y*(-0.2002690873e-3 + y*(0.8449199096e-5 + y*(-0.88228987e-6 + y*0.105787412e-6)));
    ans = sqrt( 0.636619772 / absarg ) * (cos(x)*ans1 - z*sin(x)*ans2);
    if (x < 0) ans = -ans;
  }
  return ans;
}

/* Return Bessel function Jn for integer n */

double ansi_jn(int n, double arg)
{
  int i,m,sumsign;
  double absarg,toarg,bj1,bj2,bj3,sum,ans;

  if (n < 0) {
    errno = EDOM;
    return 0;
  }
  if (n == 0) return ansi_j0(arg);
  if (n == 1) return ansi_j1(arg);
  absarg = fabs(arg);
  if (absarg == 0) return 0;
  else {
    toarg = 2/absarg;
    if ( absarg > n ) {
      bj1 = ansi_j0(arg);
      bj2 = ansi_j1(arg);
      for (i = 1; i < n; i++) {
	bj3 = i*toarg*bj2-bj1;
	bj1 = bj2;
	bj2 = bj3;
      }
      ans = bj3;
    }
    else {
      m = 2*((n + (int) sqrt(double(ACC*n)))/2 );
      sumsign = 0;
      bj3 = ans = sum = 0;
      bj2 = 1;
      for(i = m; i >0; i--) {
	bj1 = i*toarg*bj2 - bj3;
	bj3 = bj2;
	bj2 = bj1;
	if (fabs(bj2) > BIGNO ) {
	  bj2 *= BIGNI;
	  bj3 *= BIGNI;
	  ans *= BIGNI;
	  sum *= BIGNI;
	}
	if (sumsign) sum += bj2;
	sumsign = !sumsign;
	if (i == n) ans = bj3;
      }
      sum = 2*sum - bj2;
      ans /= sum;
    }
  }
  return ((arg < 0) && (n & 1)) ? -ans : ans;
}

/* Return incomplete gamma function by series method */

double ansi_gamma_series(double a, double x)
{
  int i;
  double sum,del,ap;

  if (x <= 0) {
    if (x < 0) errno = EDOM;
    return 0;
  }
  else {
    ap = a;
    del = sum = 1/a;
    for (i = 1; i <= ITMAX; i++) {
      ap++;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) return ( sum*exp(-x + a*log(x) - ansi_lgamma(a) ) );
    }
    errno = EDOM;
    return 0;
  }
}

/* Return incomplete gamma function by continued fraction method */

double ansi_gamma_cf(double a, double x)
{
  int i;
  double an,b,c,d,del,h;

  b = x+1-a;
  c = 1/FPMIN;
  h = d = 1/b;
  for (i = 1; i<=ITMAX; i++) {
    an = -i*(i-a);
    b += 2;
    d = an*d + b;
    if (fabs(d) < FPMIN ) d = FPMIN;
    c = b + an/c;
    if (fabs(c) < FPMIN ) c = FPMIN;
    d = 1/d;
    del = d*c;
    h *= del;
    if (fabs(del-1) < EPS) break;
  }
  if (i > ITMAX) errno = EDOM;
  return ( exp( -x + a*log(x) - ansi_lgamma(a) ) * h );
}

/* Return incomplete gamma function P(a,x) */

double ansi_gammap(double a, double x)
{
  if ((x < 0) || (a <= 0)) {
    errno = EDOM;
    return 0;
  }
  if ( x < a+1 ) return ansi_gamma_series(a,x);
  else return (1 - ansi_gamma_cf(a,x));
}

/* Return the error function erf(x) */

double ansi_erf(double arg)
{
  return (arg<0) ? -ansi_gammap(0.5,arg*arg) : ansi_gammap(0.5,arg*arg);
}

/* Return the complementary error function erfc(x) */

double ansi_erfc(double arg)
{
  return (arg<0) ? 1+ansi_gammap(0.5,arg*arg) : ansi_gammap(0.5,arg*arg);
}

#ifdef ACCOLD
 #define ACC ACCOLD
#endif

#ifdef BIGNOOLD
 #define BIGNO BIGNOOLD
#endif

#ifdef BIGNIOLD
 #define BIGNI BIGNIOLD
#endif

#ifdef ITMAXOLD
 #define ITMAX ITMAXOLD
#endif

#ifdef EPSOLD
 #define EPS EPSOLD
#endif

#ifdef FPMINOLD
 #define FPMIN FPMINOLD
#endif

/* Return arg! factorial */
double factorial(double arg)
{
   double Ans,Temp;
   static double a[13] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,
	39916800,479001600};
   if ((arg < 0) || (floor(arg) != arg))
   {
      errno = EDOM;
      return 0;
   };
   if (arg > 170)
   {
      errno = ERANGE;
      return HUGE_VAL;
   };
   if (arg <= 12) return a[(int)arg];
   Ans = 1;
   for (Temp = 1;Temp <= arg;Temp++)
   {
      Ans = Ans*Temp;
   };
   return Ans;
}

/* Return nCr combination */
double combination(double n,double r)
{
   double ans;
   if ((n < 0) || (r < 0) || (n < r))
   {
      errno = ERANGE;
      return 0;
   };
   if (n > 170)
   {
      errno = EDOM;
      return HUGE_VAL;
   };
   ans = factorial(n);
   ans = ans/factorial(r);
   ans = ans/factorial(n-r);
   return ans;
}

/* Return nPr permutation */
double permutation(double n,double r)
{
   double ans;
   if ((n < 0) || (r < 0) || (n < r))
   {
      errno = ERANGE;
      return 0;
   };
   if (n > 170)
   {
      errno = EDOM;
      return HUGE_VAL;
   };
   ans = factorial(n);
   ans = ans/factorial(r);
   return ans;
}

/* Return integer part of double */
double integ(double arg)
{
  double ans;
  modf(arg,&ans);
  return ans;
}

/* Return fractional part of double */
double frac(double arg)
{
  double ans;
  return modf(arg,&ans);
}

/* Return Heaviside step function */
double step(double arg)
{
  return (arg<0) ? 0 : 1;
}

/* Return sec function */
double sec(double arg)
{
  return 1/cos(arg);
}

/* Return cosec function */
double cosec(double arg)
{
  return 1/sin(arg);
}

/* Return sec function */
double cot(double arg)
{
  return 1/tan(arg);
}

