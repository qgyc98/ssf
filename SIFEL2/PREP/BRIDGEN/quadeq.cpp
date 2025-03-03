#include <stdio.h>
#include <stdlib.h>
#include "quadeq.h"


/**
  Funkce zjisti parametry kvadraticke rovnice ze dvou bodu, pricemz jeden z nich je vrchol paraboly.
  Parametry se zapisi do promennych a,b,c.
  

  @param a - parametr kvadratickeho clenu
  @param b - parametr linearniho celnu
  @param c - absolutni clen
  @param v - X-ova souradnice prvniho bodu
  @param w - Y-ova souradnice prvniho bodu
  @param x - X-ova souradnice vrcholu paraboly
  @param y - Y-ova souradnice vrcholu paraboly
*/
void getparam (double &a, double &b, double &c, double v, double w, double x, double y)
{
  double **s;
  long i;
  double pom;

  s = new double*[3];
  for (i=0;i<3;i++)
    s[i] = new double[4];

//  if (y>=w)
  if (y<=w)
    s[0][0]=2*x;
//  if (w>y)
  if (w<y)
    s[0][0]=2*v;
  s[0][1]=1.0;
  s[0][2]=0.0;
  s[0][3]=0.0;

  s[1][0]=v*v;
  s[1][1]=v;
  s[1][2]=1.0;
  s[1][3]=w;

  s[2][0]=x*x;
  s[2][1]=x;
  s[2][2]=1.0;
  s[2][3]=y;

  // vynasobeni posledniho radku pomerem clenu (-10/20)
  pom = (-1)*(s[1][0])/(s[2][0]);
  for (i=0;i<4;i++)
    s[2][i] *= pom;

  // secteni dvou poslednich radku
  for (i=0;i<4;i++)
    s[2][i] += s[1][i];

  // vynasobeni druheho radku pomerem clenu (-00/10)
  pom = (-1)*(s[0][0])/(s[1][0]);
  for (i=0;i<4;i++)
    s[1][i] *= pom;

  // secteni dvou prvnich radku
  for (i=0;i<4;i++)
    s[1][i] += s[0][i];

  // vynasobeni posledniho radku pomerem clenu (-11/21)
  pom = (-1)*(s[1][1])/(s[2][1]);
  for (i=0;i<4;i++)
    s[2][i] *= pom;

  // secteni dvou poslednich radku
  for (i=0;i<4;i++)
    s[2][i] += s[1][i];

  // vynasobeni druheho radku pomerem clenu (-22/12)
  pom = (-1)*(s[2][2])/(s[1][2]);
  for (i=0;i<4;i++)
    s[1][i] *= pom;

  // secteni dvou poslednich radku
  for (i=0;i<4;i++)
    s[1][i] += s[2][i];

  // vynasobeni prvniho radku pomerem clenu (-11/01)
  pom = (-1)*(s[1][1])/(s[0][1]);
  for (i=0;i<4;i++)
    s[0][i] *= pom;

  // secteni prvnich dvou radku
  for (i=0;i<4;i++)
    s[0][i] += s[1][i];

 // uprava na jednotkovou matici
  for (i=0;i<3;i++)
  {
    s[i][3]/=s[i][i];
    s[i][i]/=s[i][i];
  }

  a = s[0][3];
  b = s[1][3];
  c = s[2][3];
}
