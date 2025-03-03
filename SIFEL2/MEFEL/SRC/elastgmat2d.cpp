#include "elastgmat2d.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "matrix.h"


/**
  The constructor inializes attributes to zero values.
  
  Created by TKr according to JK,
*/
elastgmat2d::elastgmat2d (void)
{
  d11 = d12 = d13 = 0.0;
  d22 = d23 = 0.0;
  d33 = 0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by TKr according to JK,
*/
elastgmat2d::~elastgmat2d (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by TKr according to Tomas Koudelka,
*/
void elastgmat2d::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&d11,&d12,&d13);
  xfscanf (in,"%lf %lf",&d22,&d23);
  xfscanf (in,"%lf",&d33);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated %matrix structure for material stiffness %matrix (output)

  @return The function returns computed stiffness %matrix in the parameter d.
  
  Created by TKr according to JK,
*/
void elastgmat2d::matstiff (matrix &d,strastrestate ssst)
{
  switch (ssst){
  case planestress:{
    d[0][0]=d11;  d[0][1]=d12;  d[0][2]=d13;
    d[1][0]=d12;  d[1][1]=d22;  d[1][2]=d23;
    d[2][0]=d13;  d[2][1]=d23;  d[2][2]=d33;
    break;
  }
  case planestrain:{
    d[0][0]=d11;  d[0][1]=d12;  d[0][2]=d13;
    d[1][0]=d12;  d[1][1]=d22;  d[1][2]=d23;
    d[2][0]=d13;  d[2][1]=d23;  d[2][2]=d33;

    if (d.m > 3)//this below is not the ideal solution
      {
	d[0][3] = d[0][1]; d[1][3] = d[0][1];
	d[3][0] = d[0][3]; d[3][1] = d[1][3]; d[3][3] = d[1][1];
      }
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKr according to JK,
*/
void elastgmat2d::elmatstiff (matrix &d,strastrestate ssst)
{
  matstiff(d,ssst);
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer

  @return The function does not return anything.

  Created by TKr according to Tomas Koudelka, 1.2010
*/
void elastgmat2d::nlstresses (long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n);
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff(d,ssst);
  mxv (d,eps,sig);

  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
}




/**
  The function computes material compliance %matrix.

   
  @param c - complience %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material compliance %matrix in the parameter c.
  
  Created by TKr according to JK,
*/
void elastgmat2d::matcompl (matrix &c,strastrestate ssst)
{
  long n = c.m;
  matrix d(n,n);

  fillm(0.0,c);

  matstiff(d,ssst);

  invm(d,c,1.0e-20);
}

