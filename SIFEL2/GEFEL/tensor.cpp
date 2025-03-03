#include "tensor.h"
#include <math.h>
#include <stdlib.h>



/**
  The function transforms second order tensor tens with the help of transformation %matrix tmat 
  from the global coordinate system to the local one. The tensor is given in Voigt notation. 
  
  tens_{loc}=tmat^T tens_{glob} tmat

  @param tens - the second order tensor (input/output)
  @param tmat - transformation %matrix 
                 
  plane stress
  vector tens contains 3 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_12,

  plane strains
  vector tens contains 4 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_33, eps[3]=eps_12,

  3D cases
  vector tens contains 6 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_33,
  eps[3]=eps_23, eps[4]=eps_31, eps[5]=eps_12,

 17.1.2002
*/
void glob_loc_tens_trans (vector &tens, matrix &tmat)
{
  int n;
  
  switch ( tens.n )
    {
    case 3:
      if ( tmat.n != 2 ) 
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      if ( tmat.m != 2 )
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      n=2;
      break;
      
    case 4:
      if ( tmat.n != 3 )
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      if ( tmat.m != 3 )
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      n=3;
      break;
      
    case 6:
      if ( tmat.n != 3 )
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      if ( tmat.m != 3 )
        print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).",
                  __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
      n=3;
      break;
      
    default :
      n = 0;
      print_err("multiplying matrix by vector - incompatible size(tens.n=%ld, tmat.m=%ld, tmat.n=%ld).", 
                __FILE__, __LINE__, __func__, tens.n, tmat.m, tmat.n);
    }
  
  matrix a(n,n);
  matrix b(n,n);
  
  
  switch ( tens.n )
    {
    case 3 :
      
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[0][1]=tens[2];
      a[1][0]=tens[2];
      mtxm(tmat,a,b);
      mxm(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[0][1];
      break;
      
    case 4 :
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[2][2]=tens[2];
      a[1][2]=tens[3];
      a[2][1]=tens[3];
      mtxm(tmat,a,b);
      mxm(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[2][2];
      tens[3]=a[1][2];
      break;
      
    case 6 :
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[2][2]=tens[2];
      a[1][2]=tens[3];
      a[2][1]=tens[3];
      a[0][2]=tens[4];
      a[2][0]=tens[4];
      a[0][1]=tens[5];
      a[1][0]=tens[5];
      mtxm(tmat,a,b);
      mxm(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[2][2];
      tens[3]=a[1][2];
      tens[4]=a[0][2];
      tens[5]=a[0][1];
      break;
    }
}



/**
  The function transforms second order tensor tens with the help of transformation %matrix tmat 
  from the local coordinate system to the global one. The tensor is given in Voigt notation. 
  
  tens_{glob}=tmat tens_{loc} tmat^T

  @param tens - the second order tensor (input/output)
  @param tmat - transformation %matrix 
                 
  plane stress
  vector tens contains 3 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_12,

  plane strains
  vector tens contains 4 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_33, eps[3]=eps_12,

  3D cases
  vector tens contains 6 components
  eps[0]=eps_11, eps[1]=eps_22, eps[2]=eps_33,
  eps[3]=eps_23, eps[4]=eps_31, eps[5]=eps_12,

 17.1.2002
*/

void loc_glob_tens_trans (vector &tens, matrix &tmat)
{
  int n;
  
  switch ( tens.n )
    {
    case 3:
      if ( tmat.n != 2 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      if ( tmat.m != 2 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      n=2;
      break;
      
    case 4:
      if ( tmat.n != 3 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      if ( tmat.m != 3 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      
      n=3;
      break;
      
    case 6:
      if ( tmat.n != 3 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      if ( tmat.m != 3 ) fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
      
      n=3;
      break;
      
    default :
      n = 0;
      fprintf(stderr,"\nError multiplying matrix by vector - incompatible size (tens_trans in file tensor.*)\n");
    }
  
  matrix a(n,n);
  matrix b(n,n);
  
  
  switch ( tens.n )
    {
    case 3 :
      
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[0][1]=tens[2];
      a[1][0]=tens[2];
      mxm(tmat,a,b);
      mxmt(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[0][1];
      break;
      
    case 4 :
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[2][2]=tens[2];
      a[1][2]=tens[3];
      a[2][1]=tens[3];
      mxm(tmat,a,b);
      mxmt(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[2][2];
      tens[3]=a[1][2];
      break;
      
    case 6 :
      a[0][0]=tens[0];
      a[1][1]=tens[1];
      a[2][2]=tens[2];
      a[1][2]=tens[3];
      a[2][1]=tens[3];
      a[0][2]=tens[4];
      a[2][0]=tens[4];
      a[0][1]=tens[5];
      a[1][0]=tens[5];
      mxm(tmat,a,b);
      mxmt(b,tmat,a);
      tens[0]=a[0][0];
      tens[1]=a[1][1];
      tens[2]=a[2][2];
      tens[3]=a[1][2];
      tens[4]=a[0][2];
      tens[5]=a[0][1];
      break;
    }
}



/** 
  Function computes f(Aij), where f(x) is function of real and 
  Aij is second order tensor.

  @param a - matrix containing second order tensor
  @param f - pointer to function with one double argument and returning double
  @param nijac - number of iterations in Jacobi method
  @param error - required error of principal values
  @param zero - computer zero
  @param af - matrix containing resulting second order tensor
  
  @return The function returns result of f(A_ij) in the argument af.

  Created by Tomas Koudelka
*/
void f_tensor(matrix &a, double (*f)(double), long nijac, double error, double zero, matrix &af)
{
  vector pa(3);
  matrix t(3,3);
  long i;
  princ_val (a, pa, t, nijac, error, zero, 3, 1);
  fillm(0.0, af);
  for(i=0; i<3; i++)
    af[i][i] = (*f)(pa[i]);
  lgmatrixtransf(af, t);  
}



/** 
  Function computes f(Aij), where f(x) is function of real and 
  Aij is second order tensor. The principal directions are known 
  @param a - matrix containing second order tensor
  @param f - pointer to function with one double argument and returning double
  @param t - transformation matrix to the principal directions
  @param af - matrix containing resulting second order tensor (output)

  @return The function returns result of f(A_ij) in the argument af.

  Created by Tomas Koudelka
*/
void f_tensor(matrix &a, double (*f)(double), matrix &t, matrix &af)
{
  long i;
  copym(a, af);
  glmatrixtransf(af, t);  
  fillm(0.0, af);
  for(i=0; i<3; i++)
    af[i][i] = (*f)(a[i][i]);
  lgmatrixtransf(af, t);  
}



/** 
    Function computes derivatives of i-th component of alpha-th principal direction vector n of tensor A
    with respect to components of tensor A ie. d(n_{i,alpha})/d(A_{kl}) = (dpd_da)_{kl}
    Aij is second order tensor. The principal directions are known 
    @param t - matrix containing principal direction vectors stored in columns ie. x_g = T x_l
    @param pa - vector of principal values of A
    @param alpha - index of principal direction (index of column in t)
    @param i - index of vector component
    @param dpd_da - matrix containing resulting second order tensor

*/
long dpdir_da(matrix &t, vector &pa, long alpha, long i, matrix &dpd_da)
{
  long k,l,beta,ret;

  ret = 0;
  fillm(0.0, dpd_da);
  for(k=0; k<3; k++)
  {
    for(l=0; l<3; l++)
    {
      for(beta=0; beta<3; beta++)
      {
        if (alpha == beta)
          continue;
        if (pa[alpha] == pa[beta])
	{
          ret++;  
          continue;
	}
/*        if (pa[alpha] == pa[beta])
	{
          fprintf(stderr, "\n\n Warning - tensor has 2 or 3 same eigenvalues\n");
          fprintf(stderr, " function : dpdir_da (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
	}*/
        dpd_da[k][l] += (t[k][alpha]*t[l][beta]+t[k][beta]*t[l][alpha])*t[i][beta]/(pa[alpha]-pa[beta]);
      }
      dpd_da[k][l] /= 2.0; 
    } 
  }
  return ret;
}



/**
  The funtion computes single contraction of two second order tensors given in Voigt notation 
  c_{ij} = a_{ik} b_{kj}. The result must be given in the %matrix form 3x3.

  @param a - the second order tensor given in Voigt notation, i.e. a(6)
  @param b - the second order tensor given in Voigt notation, i.e. b(6)
  @param c - the resulting second order tensor given in %matrix form, i.e. 3x3 (output)

  @return The function stores result in the argument c.

  Created by Tomas Koudelka, 4.12.2015  
*/
void tensor_dot_prod(vector &a, vector &b, matrix &c)
{
  c(0,0) = a(0)*b(0) + a(5)*b(5) + a(4)*b(4);
  c(0,1) = a(0)*b(5) + a(5)*b(1) + a(4)*b(3);
  c(0,2) = a(0)*b(4) + a(5)*b(3) + a(4)*b(2);

  c(1,0) = a(5)*b(0) + a(1)*b(5) + a(3)*b(4);
  c(1,1) = a(5)*b(5) + a(1)*b(1) + a(3)*b(3);
  c(1,2) = a(5)*b(4) + a(1)*b(3) + a(3)*b(2);

  c(2,0) = a(4)*b(0) + a(3)*b(5) + a(2)*b(4);
  c(2,1) = a(4)*b(5) + a(3)*b(1) + a(2)*b(3);
  c(2,2) = a(4)*b(4) + a(3)*b(3) + a(2)*b(2);
}



/**
  The funtion computes single contraction of the two same second order tensors given in Voigt notation 
  b_{ij} = a_{ik} a_{kj}. The result is given in Voigt notation.

  @param a - the second order tensor given in Voigt notation, i.e. a(6)
  @param b - the resulting second order tensor given in Voigt notation, i.e. b(6)

  @return The function stores result in the argument b.

  Created by Tomas Koudelka, 4.12.2015  
*/
void tensor_dot_prod(vector &a, vector &b)
{
  b(0) = a(0)*a(0) + a(5)*a(5) + a(4)*a(4);
  b(1) = a(5)*a(5) + a(1)*a(1) + a(3)*a(3);
  b(2) = a(4)*a(4) + a(3)*a(3) + a(2)*a(2);

  b(3) = a(5)*a(4) + a(1)*a(3) + a(3)*a(2);
  b(4) = a(0)*a(4) + a(5)*a(3) + a(4)*a(2);
  b(5) = a(0)*a(5) + a(5)*a(1) + a(4)*a(3);
}



/**
   The function computes dot product of two of second order tensor stored in Voigt 
   notation. Product of the last three vector a components are scaled 
   by coefficient k. Vectors  contains 6 tensor t components in the following order:
   a[0]=t_11, a[1]=t_22, a[2]=t_33,
   a[3]=t_23, a[4]=t_31, a[5]=t_12.
  
   @param a - second order tensor stored in Voigt notation (six component %vector representation)
   @param b - second order tensor stored in Voigt notation (six component %vector representation)
   @param k - scaling factor

   Created by Tomas Koudelka, 20.4.2015
*/
double tensor_dbldot_prod (vector &a, vector &b, double k)
{
  long i;
  double nor;
  
  nor=0.0;
  for (i=0;i<3;i++){
    nor+=a(i)*b(i);
  }
  for (i=3;i<6;i++){
    nor+=k*a(i)*b(i);
  }
  return nor;
}



/**
  The function computes double contraction of the two second order tensors, 
  i.e. computes a_ij b_ij.
   
  @param a - second order tensor given in %matrix form 3x3
  @param b - second order tensor given in %matrix form 3x3
   
  @return The function returns a_{ij} b_{ij}

  JK, 9.8.2005
  Renamed, modified by TKo, 4.12.2015
*/
double tensor_dbldot_prod (matrix &a, matrix &b)
{
  long i,j,n;
  double ip = 0.0;
  
  if ((a.m == 3) && (a.n == 3) && (b.m == 3) && (b.n == 3))
  {
    n=a.n;
    for (i=0;i<n;i++){
      for (j=0;j<n;j++)
        ip+=a[i][j]*b[i][j];
    }
  }
  else
  {
    print_err("wrong dimensions of tensor(s) - a(%ld,%ld), b(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n);
    abort();
  }
  
  return ip;
}




