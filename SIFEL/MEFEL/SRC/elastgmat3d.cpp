#include "elastgmat3d.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "matrix.h"


/**
  The constructor inializes attributes to zero values.
  
  Created by JK,
*/
elastgmat3d::elastgmat3d (void)
{
  d11 = d12 = d13 = d14 = d15 = d16 = 0.0;
  d22 = d23 = d24 = d25 = d26 = 0.0;
  d33 = d34 = d35 = d36 = 0.0;
  d44 = d45 = d46 = 0.0;
  d55 = d56 = 0.0;
  d66 = 0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by JK,
*/
elastgmat3d::~elastgmat3d (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void elastgmat3d::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf",&d11,&d12,&d13,&d14,&d15,&d16);
  xfscanf (in,"%lf %lf %lf %lf %lf",&d22,&d23,&d24,&d25,&d26);
  xfscanf (in,"%lf %lf %lf %lf",&d33,&d34,&d35,&d36);
  xfscanf (in,"%lf %lf %lf",&d44,&d45,&d46);
  xfscanf (in,"%lf %lf",&d55,&d56);
  xfscanf (in,"%lf",&d66);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated %matrix structure for material stiffness %matrix (output)

  @return The function returns computed stiffness %matrix in the parameter d.
  
  Created by JK,
*/
void elastgmat3d::matstiff (matrix &d)
{
  d[0][0]=d11;  d[0][1]=d12;  d[0][2]=d13;  d[0][3]=d14;  d[0][4]=d15;  d[0][5]=d16;
  d[1][0]=d12;  d[1][1]=d22;  d[1][2]=d23;  d[1][3]=d24;  d[1][4]=d25;  d[1][5]=d26;
  d[2][0]=d13;  d[2][1]=d23;  d[2][2]=d33;  d[2][3]=d34;  d[2][4]=d35;  d[2][5]=d36;
  d[3][0]=d14;  d[3][1]=d24;  d[3][2]=d34;  d[3][3]=d44;  d[3][4]=d45;  d[3][5]=d46;
  d[4][0]=d15;  d[4][1]=d25;  d[4][2]=d35;  d[4][3]=d45;  d[4][4]=d55;  d[4][5]=d56;
  d[5][0]=d16;  d[5][1]=d26;  d[5][2]=d36;  d[5][3]=d46;  d[5][4]=d56;  d[5][5]=d66;
}



/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK,
*/
void elastgmat3d::elmatstiff (matrix &d)
{
  matstiff(d);
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer

  @return The function does not return anything.

  Created by Tomas Koudelka, 1.2010
*/
void elastgmat3d::nlstresses (long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n);
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff(d);
  mxv (d,eps,sig);

  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
}
