#include <stdlib.h>
#include <math.h>
#include "linhexrot.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "globmat.h"
#include "genfile.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"




linhexrot::linhexrot ()
{
  long i,j;

  //  number of nodes on element
  nne=8;
  //  number of DOFs on element
  ndofe=48;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=2;
  //  number of edges on element
  ned=12;
  //  number of nodes on one edge
  nned=2;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;
  //  number of surfaces on element
  nsurf=6;
  //  number of nodes on one surface
  nnsurf=4;
  //  strain/stress state
  ssst=spacestress;
  
  //  number of blocks (parts of geometric matrix)
  nb=1;
  
  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=6;

  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  
  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=8;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=2;
}

linhexrot::~linhexrot ()
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;
  
  delete [] cncomp;
  delete [] ncomp;
}



/**
   function approximates function defined by nodal values
   
   @param xi,eta,zeta - natural coordinates
   @param nodval - nodal values
   
   JK, 20.8.2001
*/
double linhexrot::approx (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_hex_3d (bf.a,xi,eta,zeta);
  
  scprd (bf,nodval,f);

  return f;
}

  //coordinate corners - 1(1,1,1),  2(-1,1,1),  3(-1,-1,1),  4(1,-1,1)
  //                     5(1,1,-1), 6(-1,1,-1), 7(-1,-1,-1), 8(1,-1,-1)
  //coord. inside   between node 1-2,  2-3,  3-4,  4-1, 
  //                             1-5,  2-6,  3-7,  4-8,
  //                             5-6,  6-7,  7-8,  8-5,


void linhexrot::tran_side (matrix &t12,matrix &t23,matrix &t34,matrix &t41,  matrix &t15, matrix &t26,
						   matrix &t37,matrix &t48,  matrix &t56,matrix &t67,matrix &t78, matrix &t85,   vector &x,vector &y,vector &z)
{
  long i;
  double pom;
  vector an(3),am(3),dl(ned);
  matrix a(3,ned);
  
  a[0][0]=x[1]-x[0]; a[0][1]=x[2]-x[1]; a[0][2] =x[3]-x[2]; a[0][3] =x[0]-x[3]; 
  a[0][4]=x[4]-x[0]; a[0][5]=x[5]-x[1]; a[0][6] =x[6]-x[2]; a[0][7] =x[7]-x[3];
  a[0][8]=x[5]-x[4]; a[0][9]=x[6]-x[5]; a[0][10]=x[7]-x[6]; a[0][11]=x[4]-x[7];

  a[1][0]=y[1]-y[0]; a[1][1]=y[2]-y[1]; a[1][2] =y[3]-y[2]; a[1][3] =y[0]-y[3]; 
  a[1][4]=y[4]-y[0]; a[1][5]=y[5]-y[1]; a[1][6] =y[6]-y[2]; a[1][7] =y[7]-y[3];
  a[1][8]=y[5]-y[4]; a[1][9]=y[6]-y[5]; a[1][10]=y[7]-y[6]; a[1][11]=y[4]-y[7];

  a[2][0]=z[1]-z[0]; a[2][1]=z[2]-z[1]; a[2][2] =z[3]-z[2]; a[2][3] =z[0]-z[3]; 
  a[2][4]=z[4]-z[0]; a[2][5]=z[5]-z[1]; a[2][6] =z[6]-z[2]; a[2][7] =z[7]-z[3];
  a[2][8]=z[5]-z[4]; a[2][9]=z[6]-z[5]; a[2][10]=z[7]-z[6]; a[2][11]=z[4]-z[7];
  
  
  for (i=0;i<ned;i++){
    dl[i]=sqrt(a[0][i]*a[0][i]+a[1][i]*a[1][i]+a[2][i]*a[2][i]);
    if(fabs(a[0][i])<=1.e-6 && fabs(a[1][i])<=1.e-6)
      {am[0]=1.;  am[1]=0.; am[2]=0.;
      }
    else
      {pom=sqrt(a[0][i]*a[0][i]+a[1][i]*a[1][i]);
        am[0]=-a[1][i]/pom; am[1]= a[0][i]/pom; am[2]=0.;
      }
    
    an[0]=-am[1]*a[2][i]/dl[i];
    an[1]= am[0]*a[2][i]/dl[i];
    an[2]= am[1]*a[0][i]/dl[i]-am[0]*a[1][i]/dl[i];
    
    if (i==0){ 
      t12[0][0]=0.;          t12[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t12[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t12[1][0]=-t12[0][1];  t12[1][1]=0.;                                  t12[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t12[2][0]=-t12[0][2];  t12[2][1]=-t12[1][2];                          t12[2][2]=0.;
    }
    else if  (i==1){
      t23[0][0]=0.;          t23[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t23[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t23[1][0]=-t23[0][1];  t23[1][1]=0.;                                  t23[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t23[2][0]=-t23[0][2];  t23[2][1]=-t23[1][2];                          t23[2][2]=0.;
    }
    else if  (i==2){
      t34[0][0]=0.;          t34[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t34[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t34[1][0]=-t34[0][1];  t34[1][1]=0.;                                  t34[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t34[2][0]=-t34[0][2];  t34[2][1]=-t34[1][2];                          t34[2][2]=0.;
    }
    else if  (i==3){
      t41[0][0]=0.;          t41[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t41[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t41[1][0]=-t41[0][1];  t41[1][1]=0.;                                  t41[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t41[2][0]=-t41[0][2];  t41[2][1]=-t41[1][2];                          t41[2][2]=0.;
    }
    else if  (i==4){
      t15[0][0]=0.;          t15[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t15[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t15[1][0]=-t15[0][1];  t15[1][1]=0.;                                  t15[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t15[2][0]=-t15[0][2];  t15[2][1]=-t15[1][2];                          t15[2][2]=0.;
    }
    else if  (i==5){
      t26[0][0]=0.;          t26[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t26[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t26[1][0]=-t26[0][1];  t26[1][1]=0.;                                  t26[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t26[2][0]=-t26[0][2];  t26[2][1]=-t26[1][2];                          t26[2][2]=0.;
    }
    else if  (i==6){
      t37[0][0]=0.;          t37[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t37[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t37[1][0]=-t37[0][1];  t37[1][1]=0.;                                  t37[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t37[2][0]=-t37[0][2];  t37[2][1]=-t37[1][2];                          t37[2][2]=0.;
    }
    else if  (i==7){
      t48[0][0]=0.;          t48[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t48[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t48[1][0]=-t48[0][1];  t48[1][1]=0.;                                  t48[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t48[2][0]=-t48[0][2];  t48[2][1]=-t48[1][2];                          t48[2][2]=0.;
    }
     else if  (i==8){
      t56[0][0]=0.;          t56[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t56[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t56[1][0]=-t56[0][1];  t56[1][1]=0.;                                  t56[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t56[2][0]=-t56[0][2];  t56[2][1]=-t56[1][2];                          t56[2][2]=0.;
    }
    else if  (i==9){
      t67[0][0]=0.;          t67[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t67[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t67[1][0]=-t67[0][1];  t67[1][1]=0.;                                  t67[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t67[2][0]=-t67[0][2];  t67[2][1]=-t67[1][2];                          t67[2][2]=0.;
    }
    else if  (i==10){
      t78[0][0]=0.;          t78[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t78[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t78[1][0]=-t78[0][1];  t78[1][1]=0.;                                  t78[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t78[2][0]=-t78[0][2];  t78[2][1]=-t78[1][2];                          t78[2][2]=0.;
    }
    else if  (i==11){
      t85[0][0]=0.;          t85[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t85[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t85[1][0]=-t85[0][1];  t85[1][1]=0.;                                  t85[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t85[2][0]=-t85[0][2];  t85[2][1]=-t85[1][2];                          t85[2][2]=0.;
    }

  }
  
}
/**
   function assembles %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param xi,eta,zeta - natural coordinates

   PF, 19.10.2008
*/
void linhexrot::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  long i,j,k,l;
  vector bf(nne);

  fillm (0.0,n);

  bf_lin_hex_3d (bf.a,xi,eta,zeta);
  
  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];  j+=3;
    n[1][k]=bf[i];  k+=3;
    n[2][l]=bf[i];  l+=3;
  }
}



/**
   function assembles strain-displacement (geometric) %matrix

   @param gm - geometric %matrix
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian

   JK, 19.7.2001
*/
void linhexrot::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac)
{
  long i,i1,i2,i3,i4,i5,i6,i7,i8;
  vector dx(nne),dy(nne),dz(nne), dxr(ned),dyr(ned),dzr(ned);
  matrix t12(3,3),t23(3,3),t34(3,3),t41(3,3), t15(3,3),t26(3,3),t37(3,3),t48(3,3), t56(3,3),t67(3,3),t78(3,3),t85(3,3);

  dx_bf_lin_hex_3d (dx.a, eta,zeta);
  dy_bf_lin_hex_3d (dy.a, xi,zeta);
  dz_bf_lin_hex_3d (dz.a, xi,eta);
  derivatives_3d (dx,dy,dz,jac, x,y,z,xi,eta,zeta);

  // ds12[0] is derivation ( d(12)/dx ), ds12[1] is ( d(12)/dy ), ds12[2] is ( d(12)/dz )
  dx_bf_quad_hexrot_3d (dxr.a, xi,eta,zeta);
  dy_bf_quad_hexrot_3d (dyr.a, xi,eta,zeta);
  dz_bf_quad_hexrot_3d (dzr.a, xi,eta,zeta);
  derivatives_3d (dxr,dyr,dzr,jac, x,y,z,xi,eta,zeta);

  tran_side (t12,t23,t34,t41, t15,t26,t37,t48, t56,t67,t78,t85,   x,y,z);
  
  fillm (0.0,gm);

  i1=0;  i2=1;  i3=2;
  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[1][i2]=dy[i];
    gm[2][i3]=dz[i];

    gm[3][i2]=dz[i];
    gm[3][i3]=dy[i];
    gm[4][i1]=dz[i];
    gm[4][i3]=dx[i];
    gm[5][i1]=dy[i];
    gm[5][i2]=dx[i];
    i1+=6;  i2+=6;  i3+=6;
  }
  
  i1=3; i2=9; i3=15; i4=21; i5=27; i6=33; i7=39; i8=45;
  for (i=0;i<3;i++){
    // node 1  du/dx, dv/dy, dw/dz, 
    gm[0][i1]=-t12[0][i]*dxr[0]+t41[0][i]*dxr[3]-t15[0][i]*dxr[4];
    gm[1][i1]=-t12[1][i]*dyr[0]+t41[1][i]*dyr[3]-t15[1][i]*dyr[4];
    gm[2][i1]=-t12[2][i]*dzr[0]+t41[2][i]*dzr[3]-t15[2][i]*dzr[4];
	//node 1   dv/dz+dw/dy,  du/dz+dw/dx,  du/dy+dv/dx
    gm[3][i1]=-t12[1][i]*dzr[0]+t41[1][i]*dzr[3]-t15[1][i]*dzr[4]
	          -t12[2][i]*dyr[0]+t41[2][i]*dyr[3]-t15[2][i]*dyr[4];
    gm[4][i1]=-t12[2][i]*dxr[0]+t41[2][i]*dxr[3]-t15[2][i]*dxr[4]
	          -t12[0][i]*dzr[0]+t41[0][i]*dzr[3]-t15[0][i]*dzr[4];
    gm[5][i1]=-t12[0][i]*dyr[0]+t41[0][i]*dyr[3]-t15[0][i]*dyr[4]
	          -t12[1][i]*dxr[0]+t41[1][i]*dxr[3]-t15[1][i]*dxr[4];
    // node 2
    gm[0][i2]=t12[0][i]*dxr[0]-t23[0][i]*dxr[1]-t26[0][i]*dxr[5];
    gm[1][i2]=t12[1][i]*dyr[0]-t23[1][i]*dyr[1]-t26[1][i]*dyr[5];
    gm[2][i2]=t12[2][i]*dzr[0]-t23[2][i]*dzr[1]-t26[2][i]*dzr[5];
    gm[3][i2]=t12[1][i]*dzr[0]-t23[1][i]*dzr[1]-t26[1][i]*dzr[5]
	         +t12[2][i]*dyr[0]-t23[2][i]*dyr[1]-t26[2][i]*dyr[5];
    gm[4][i2]=t12[2][i]*dxr[0]-t23[2][i]*dxr[1]-t26[2][i]*dxr[5]
	         +t12[0][i]*dzr[0]-t23[0][i]*dzr[1]-t26[0][i]*dzr[5];
    gm[5][i2]=t12[0][i]*dyr[0]-t23[0][i]*dyr[1]-t26[0][i]*dyr[5]
	         +t12[1][i]*dxr[0]-t23[1][i]*dxr[1]-t26[1][i]*dxr[5];
    // node 3
    gm[0][i3]=t23[0][i]*dxr[1]-t34[0][i]*dxr[2]-t37[0][i]*dxr[6];
    gm[1][i3]=t23[1][i]*dyr[1]-t34[1][i]*dyr[2]-t37[1][i]*dyr[6];
    gm[2][i3]=t23[2][i]*dzr[1]-t34[2][i]*dzr[2]-t37[2][i]*dzr[6];
    gm[3][i3]=t23[1][i]*dzr[1]-t34[1][i]*dzr[2]-t37[1][i]*dzr[6]
	         +t23[2][i]*dyr[1]-t34[2][i]*dyr[2]-t37[2][i]*dyr[6];
    gm[4][i3]=t23[2][i]*dxr[1]-t34[2][i]*dxr[2]-t37[2][i]*dxr[6]
	         +t23[0][i]*dzr[1]-t34[0][i]*dzr[2]-t37[0][i]*dzr[6];
    gm[5][i3]=t23[0][i]*dyr[1]-t34[0][i]*dyr[2]-t37[0][i]*dyr[6]
	         +t23[1][i]*dxr[1]-t34[1][i]*dxr[2]-t37[1][i]*dxr[6];
    // node 4
    gm[0][i4]=t34[0][i]*dxr[2]-t41[0][i]*dxr[3]-t48[0][i]*dxr[7];
    gm[1][i4]=t34[1][i]*dyr[2]-t41[1][i]*dyr[3]-t48[1][i]*dyr[7];
    gm[2][i4]=t34[2][i]*dzr[2]-t41[2][i]*dzr[3]-t48[2][i]*dzr[7];
	gm[3][i4]=t34[1][i]*dzr[2]-t41[1][i]*dzr[3]-t48[1][i]*dzr[7]
	         +t34[2][i]*dyr[2]-t41[2][i]*dyr[3]-t48[2][i]*dyr[7];
    gm[4][i4]=t34[2][i]*dxr[2]-t41[2][i]*dxr[3]-t48[2][i]*dxr[7]
	         +t34[0][i]*dzr[2]-t41[0][i]*dzr[3]-t48[0][i]*dzr[7];
    gm[5][i4]=t34[0][i]*dyr[2]-t41[0][i]*dyr[3]-t48[0][i]*dyr[7]
	         +t34[1][i]*dxr[2]-t41[1][i]*dxr[3]-t48[1][i]*dxr[7];
    // node 5
    gm[0][i5]=t15[0][i]*dxr[4]-t56[0][i]*dxr[8]+t85[0][i]*dxr[11];
    gm[1][i5]=t15[1][i]*dyr[4]-t56[1][i]*dyr[8]+t85[1][i]*dyr[11];
    gm[2][i5]=t15[2][i]*dzr[4]-t56[2][i]*dzr[8]+t85[2][i]*dzr[11];
	gm[3][i5]=t15[1][i]*dzr[4]-t56[1][i]*dzr[8]+t85[1][i]*dzr[11]
	         +t15[2][i]*dyr[4]-t56[2][i]*dyr[8]+t85[2][i]*dyr[11];
    gm[4][i5]=t15[2][i]*dxr[4]-t56[2][i]*dxr[8]+t85[2][i]*dxr[11]
	         +t15[0][i]*dzr[4]-t56[0][i]*dzr[8]+t85[0][i]*dzr[11];
    gm[5][i5]=t15[0][i]*dyr[4]-t56[0][i]*dyr[8]+t85[0][i]*dyr[11]
	         +t15[1][i]*dxr[4]-t56[1][i]*dxr[8]+t85[1][i]*dxr[11];
    // node 6
    gm[0][i6]=t26[0][i]*dxr[5]+t56[0][i]*dxr[8]-t67[0][i]*dxr[9];
    gm[1][i6]=t26[1][i]*dyr[5]+t56[1][i]*dyr[8]-t67[1][i]*dyr[9];
    gm[2][i6]=t26[2][i]*dzr[5]+t56[2][i]*dzr[8]-t67[2][i]*dzr[9];
    gm[3][i6]=t26[1][i]*dzr[5]+t56[1][i]*dzr[8]-t67[1][i]*dzr[9]
	         +t26[2][i]*dyr[5]+t56[2][i]*dyr[8]-t67[2][i]*dyr[9];
    gm[4][i6]=t26[2][i]*dxr[5]+t56[2][i]*dxr[8]-t67[2][i]*dxr[9]
	         +t26[0][i]*dzr[5]+t56[0][i]*dzr[8]-t67[0][i]*dzr[9];
    gm[5][i6]=t26[0][i]*dyr[5]+t56[0][i]*dyr[8]-t67[0][i]*dyr[9]
	         +t26[1][i]*dxr[5]+t56[1][i]*dxr[8]-t67[1][i]*dxr[9];
    // node 7
    gm[0][i7]=t37[0][i]*dxr[6]+t67[0][i]*dxr[9]-t78[0][i]*dxr[10];
    gm[1][i7]=t37[1][i]*dyr[6]+t67[1][i]*dyr[9]-t78[1][i]*dyr[10];
    gm[2][i7]=t37[2][i]*dzr[6]+t67[2][i]*dzr[9]-t78[2][i]*dzr[10];
    gm[3][i7]=t37[1][i]*dzr[6]+t67[1][i]*dzr[9]-t78[1][i]*dzr[10]
	         +t37[2][i]*dyr[6]+t67[2][i]*dyr[9]-t78[2][i]*dyr[10];
    gm[4][i7]=t37[2][i]*dxr[6]+t67[2][i]*dxr[9]-t78[2][i]*dxr[10]
	         +t37[0][i]*dzr[6]+t67[0][i]*dzr[9]-t78[0][i]*dzr[10];
    gm[5][i7]=t37[0][i]*dyr[6]+t67[0][i]*dyr[9]-t78[0][i]*dyr[10]
	         +t37[1][i]*dxr[6]+t67[1][i]*dxr[9]-t78[1][i]*dxr[10];
    // node 8
    gm[0][i8]=t48[0][i]*dxr[7]+t78[0][i]*dxr[10]-t85[0][i]*dxr[11];
    gm[1][i8]=t48[1][i]*dyr[7]+t78[1][i]*dyr[10]-t85[1][i]*dyr[11];
    gm[2][i8]=t48[2][i]*dzr[7]+t78[2][i]*dzr[10]-t85[2][i]*dzr[11];
    gm[3][i8]=t48[1][i]*dzr[7]+t78[1][i]*dzr[10]-t85[1][i]*dzr[11]
	         +t48[2][i]*dyr[7]+t78[2][i]*dyr[10]-t85[2][i]*dyr[11];
    gm[4][i8]=t48[2][i]*dxr[7]+t78[2][i]*dxr[10]-t85[2][i]*dxr[11]
	         +t48[0][i]*dzr[7]+t78[0][i]*dzr[10]-t85[0][i]*dzr[11];
    gm[5][i8]=t48[0][i]*dyr[7]+t78[0][i]*dyr[10]-t85[0][i]*dyr[11]
	         +t48[1][i]*dxr[7]+t78[1][i]*dxr[10]-t85[1][i]*dxr[11];
    i1+=1; i2+=1; i3+=1; i4+=1;  i5+=1; i6+=1; i7+=1; i8+=1;
  }
  /*
  fprintf (Out,"\n\n gm");
    for (i1=0;i1<6;i1++){
      fprintf (Out,"\n");
      for (i2=0;i2<ndofe;i2++){
      fprintf (Out," %le",gm[i1][i2]);
	  }
      fprintf (Out,"\n");
	}
  */

}

void linhexrot::geom_matrix_shear (matrix &gm,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac)
{
  long i,i1,i2,i3,i4,i5,i6,i7,i8;
  vector dx(nne),dy(nne),dz(nne), dxr(ned),dyr(ned),dzr(ned), bf(nne);
  matrix t12(3,3),t23(3,3),t34(3,3),t41(3,3), t15(3,3),t26(3,3),t37(3,3),t48(3,3), t56(3,3),t67(3,3),t78(3,3),t85(3,3);

  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);
  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);

  // ds12[0] is derivation ( d(12)/dx ), ds12[1] is ( d(12)/dy ), ds12[2] is ( d(12)/dz )
  dx_bf_quad_hexrot_3d (dxr.a, xi,eta,zeta);
  dy_bf_quad_hexrot_3d (dyr.a, xi,eta,zeta);
  dz_bf_quad_hexrot_3d (dzr.a, xi,eta,zeta);
  derivatives_3d (dxr,dyr,dzr,jac, x,y,z,xi,eta,zeta);

  tran_side (t12,t23,t34,t41, t15,t26,t37,t48, t56,t67,t78,t85,   x,y,z);
  
  fillm (0.0,gm);

  
  i1=3; i2=9; i3=15; i4=21; i5=27; i6=33; i7=39; i8=45;
  for (i=0;i<3;i++){
	//node 1   -dv/dz+dw/dy,  du/dz-dw/dx,  -du/dy+dv/dx
	gm[0][i1]=(t12[1][i]*dzr[0]-t41[1][i]*dzr[3]+t15[1][i]*dzr[4]
	          -t12[2][i]*dyr[0]+t41[2][i]*dyr[3]-t15[2][i]*dyr[4])/2.;
	gm[1][i1]=(t12[2][i]*dxr[0]-t41[2][i]*dxr[3]+t15[2][i]*dxr[4]
	          -t12[0][i]*dzr[0]+t41[0][i]*dzr[3]-t15[0][i]*dzr[4])/2.;
    gm[2][i1]=(t12[0][i]*dyr[0]-t41[0][i]*dyr[3]+t15[0][i]*dyr[4]
	          -t12[1][i]*dxr[0]+t41[1][i]*dxr[3]-t15[1][i]*dxr[4])/2.;
    // node 2
    gm[0][i2]=(-t12[1][i]*dzr[0]+t23[1][i]*dzr[1]+t26[1][i]*dzr[5]
	           +t12[2][i]*dyr[0]-t23[2][i]*dyr[1]-t26[2][i]*dyr[5])/2.;
    gm[1][i2]=(-t12[2][i]*dxr[0]+t23[2][i]*dxr[1]+t26[2][i]*dxr[5]
	           +t12[0][i]*dzr[0]-t23[0][i]*dzr[1]-t26[0][i]*dzr[5])/2.;
    gm[2][i2]=(-t12[0][i]*dyr[0]+t23[0][i]*dyr[1]+t26[0][i]*dyr[5]
	           +t12[1][i]*dxr[0]-t23[1][i]*dxr[1]-t26[1][i]*dxr[5])/2.;
    // node 3
	gm[0][i3]=(-t23[1][i]*dzr[1]+t34[1][i]*dzr[2]+t37[1][i]*dzr[6]
	           +t23[2][i]*dyr[1]-t34[2][i]*dyr[2]-t37[2][i]*dyr[6])/2.;
    gm[1][i3]=(-t23[2][i]*dxr[1]+t34[2][i]*dxr[2]+t37[2][i]*dxr[6]
	           +t23[0][i]*dzr[1]-t34[0][i]*dzr[2]-t37[0][i]*dzr[6])/2.;
    gm[2][i3]=(-t23[0][i]*dyr[1]+t34[0][i]*dyr[2]+t37[0][i]*dyr[6]
	           +t23[1][i]*dxr[1]-t34[1][i]*dxr[2]-t37[1][i]*dxr[6])/2.;
    // node 4
    gm[0][i4]=(-t34[1][i]*dzr[2]+t41[1][i]*dzr[3]+t48[1][i]*dzr[7]
	           +t34[2][i]*dyr[2]-t41[2][i]*dyr[3]-t48[2][i]*dyr[7])/2.;
    gm[1][i4]=(-t34[2][i]*dxr[2]+t41[2][i]*dxr[3]+t48[2][i]*dxr[7]
	           +t34[0][i]*dzr[2]-t41[0][i]*dzr[3]-t48[0][i]*dzr[7])/2.;
    gm[2][i4]=(-t34[0][i]*dyr[2]+t41[0][i]*dyr[3]+t48[0][i]*dyr[7]
	           +t34[1][i]*dxr[2]-t41[1][i]*dxr[3]-t48[1][i]*dxr[7])/2.;
    // node 5
    gm[0][i5]=(-t15[1][i]*dzr[4]+t56[1][i]*dzr[8]-t85[1][i]*dzr[11]
	           +t15[2][i]*dyr[4]-t56[2][i]*dyr[8]+t85[2][i]*dyr[11])/2.;
    gm[1][i5]=(-t15[2][i]*dxr[4]+t56[2][i]*dxr[8]-t85[2][i]*dxr[11]
	           +t15[0][i]*dzr[4]-t56[0][i]*dzr[8]+t85[0][i]*dzr[11])/2.;
    gm[2][i5]=(-t15[0][i]*dyr[4]+t56[0][i]*dyr[8]-t85[0][i]*dyr[11]
	           +t15[1][i]*dxr[4]-t56[1][i]*dxr[8]+t85[1][i]*dxr[11])/2.;
    // node 6
    gm[0][i6]=(-t26[1][i]*dzr[5]-t56[1][i]*dzr[8]+t67[1][i]*dzr[9]
	           +t26[2][i]*dyr[5]+t56[2][i]*dyr[8]-t67[2][i]*dyr[9])/2.;
    gm[1][i6]=(-t26[2][i]*dxr[5]-t56[2][i]*dxr[8]+t67[2][i]*dxr[9]
	           +t26[0][i]*dzr[5]+t56[0][i]*dzr[8]-t67[0][i]*dzr[9])/2.;
    gm[2][i6]=(-t26[0][i]*dyr[5]-t56[0][i]*dyr[8]+t67[0][i]*dyr[9]
	           +t26[1][i]*dxr[5]+t56[1][i]*dxr[8]-t67[1][i]*dxr[9])/2.;
    // node 7
    gm[0][i7]=(-t37[1][i]*dzr[6]-t67[1][i]*dzr[9]+t78[1][i]*dzr[10]
	           +t37[2][i]*dyr[6]+t67[2][i]*dyr[9]-t78[2][i]*dyr[10])/2.;
    gm[1][i7]=(-t37[2][i]*dxr[6]-t67[2][i]*dxr[9]+t78[2][i]*dxr[10]
	           +t37[0][i]*dzr[6]+t67[0][i]*dzr[9]-t78[0][i]*dzr[10])/2.;
    gm[2][i7]=(-t37[0][i]*dyr[6]-t67[0][i]*dyr[9]+t78[0][i]*dyr[10]
	           +t37[1][i]*dxr[6]+t67[1][i]*dxr[9]-t78[1][i]*dxr[10])/2.;
    // node 8
    gm[0][i8]=(-t48[1][i]*dzr[7]-t78[1][i]*dzr[10]+t85[1][i]*dzr[11]
	           +t48[2][i]*dyr[7]+t78[2][i]*dyr[10]-t85[2][i]*dyr[11])/2.;
    gm[1][i8]=(-t48[2][i]*dxr[7]-t78[2][i]*dxr[10]+t85[2][i]*dxr[11]
	           +t48[0][i]*dzr[7]+t78[0][i]*dzr[10]-t85[0][i]*dzr[11])/2.;
    gm[2][i8]=(-t48[0][i]*dyr[7]-t78[0][i]*dyr[10]+t85[0][i]*dyr[11]
	           +t48[1][i]*dxr[7]+t78[1][i]*dxr[10]-t85[1][i]*dxr[11])/2.;
    i1+=1; i2+=1; i3+=1; i4+=1;  i5+=1; i6+=1; i7+=1; i8+=1;
  }
  // -dv/dz+dw/dy,  du/dz-dw/dx,  -du/dy+dv/dx
   
  bf_lin_hex_3d (bf.a,xi,eta,zeta);
  i1=0;  i2=1;  i3=2;
  for (i=0;i<nne;i++){
  //  gm[0][i2]   =-dz[i]/2.;
  //  gm[0][i3]   = dy[i]/2.;
    gm[0][i1+3]-= bf.a[i];

  //  gm[1][i1]   = dz[i]/2.;  
  //  gm[1][i3]   =-dx[i]/2.;
    gm[1][i1+4]-= bf.a[i];

  //  gm[2][i1]   =-dy[i]/2.; 
  //  gm[2][i2]   = dx[i]/2.;  
    gm[2][i1+5]-= bf.a[i];
    i1+=6;  i2+=6;  i3+=6;
  }
  


}

/**
   function assembles auxiliary vectors B for evaluation of stiffness %matrix
   in geometrically nonlinear problems
   
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian
   @param b11,b12,b13,b21,b22,b23,b31,b32,b33 - vectors of derivatives of shape functions
   
   JK, 24.9.2005
*/
void linhexrot::bvectors (vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac,
		       vector &b11,vector &b12,vector &b13,
		       vector &b21,vector &b22,vector &b23,
		       vector &b31,vector &b32,vector &b33)
{
  vector dx(nne),dy(nne),dz(nne);
  
  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);
  
  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);
  
  fillv (0.0,b11);
  fillv (0.0,b12);
  fillv (0.0,b13);
  fillv (0.0,b21);
  fillv (0.0,b22);
  fillv (0.0,b23);
  fillv (0.0,b31);
  fillv (0.0,b32);
  fillv (0.0,b33);

  //  du/dx
  b11[0]=dx[0];  b11[3]=dx[1];  b11[6]=dx[2];  b11[9]=dx[3];   b11[12]=dx[4];  b11[15]=dx[5];  b11[18]=dx[6];  b11[21]=dx[7];
  //  du/dy
  b12[0]=dy[0];  b12[3]=dy[1];  b12[6]=dy[2];  b12[9]=dy[3];   b12[12]=dy[4];  b12[15]=dy[5];  b12[18]=dy[6];  b12[21]=dy[7];
  //  du/dz
  b13[0]=dz[0];  b13[3]=dz[1];  b13[6]=dz[2];  b13[9]=dz[3];   b13[12]=dz[4];  b13[15]=dz[5];  b13[18]=dz[6];  b13[21]=dz[7];

  //  dv/dx
  b21[1]=dx[0];  b21[4]=dx[1];  b21[7]=dx[2];  b21[10]=dx[3];  b21[13]=dx[4];  b21[16]=dx[5];  b21[19]=dx[6];  b21[22]=dx[7];
  //  dv/dy
  b22[1]=dy[0];  b22[4]=dy[1];  b22[7]=dy[2];  b22[10]=dy[3];  b22[13]=dy[4];  b22[16]=dy[5];  b22[19]=dy[6];  b22[22]=dy[7];
  //  dv/dz
  b23[1]=dz[0];  b23[4]=dz[1];  b23[7]=dz[2];  b23[10]=dz[3];  b23[13]=dz[4];  b23[16]=dz[5];  b23[19]=dz[6];  b23[22]=dz[7];

  //  dw/dx
  b31[2]=dx[0];  b31[5]=dx[1];  b31[8]=dx[2];  b31[11]=dx[3];  b31[14]=dx[4];  b31[17]=dx[5];  b31[20]=dx[6];  b31[23]=dx[7];
  //  dw/dy
  b32[2]=dy[0];  b32[5]=dy[1];  b32[8]=dy[2];  b32[11]=dy[3];  b32[14]=dy[4];  b32[17]=dy[5];  b32[20]=dy[6];  b32[23]=dy[7];
  //  dw/dz
  b33[2]=dz[0];  b33[5]=dz[1];  b33[8]=dz[2];  b33[11]=dz[3];  b33[14]=dz[4];  b33[17]=dz[5];  b33[20]=dz[6];  b33[23]=dz[7];

}

/**
   function computes strain-displacement %matrix for geometrically nonlinear problems
   
   @param gm - strain-displacement %matrix
   @param r - array of nodal displacements
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian

   JK, 24.9.2005
*/
void linhexrot::gngeom_matrix (matrix &gm,vector &r,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  double b11r,b12r,b13r,b21r,b22r,b23r,b31r,b32r,b33r;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe),av(ndofe);
  
  fillm (0.0,gm);
  
  bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
  
  scprd (b11,r,b11r);
  scprd (b12,r,b12r);
  scprd (b13,r,b13r);
  scprd (b21,r,b21r);
  scprd (b22,r,b22r);
  scprd (b23,r,b23r);
  scprd (b31,r,b31r);
  scprd (b32,r,b32r);
  scprd (b33,r,b33r);
  
  
  // *******
  //  E_11
  // *******
  
  //  B11 dr
  for (i=0;i<ndofe;i++){
    gm[0][i]+=b11[i];
  }
  
  //  r B11 B11 dr
  cmulv(b11r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }
  
  //  r B21 B21 dr
  cmulv(b21r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }

  //  r B31 B31 dr
  cmulv(b31r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }

  
  // *******
  //  E_22
  // *******
  
  //  B22 dr
  for (i=0;i<ndofe;i++){
    gm[1][i]+=b22[i];
  }
  
  //  r B12 B12 dr
  cmulv(b12r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  //  r B22 B22 dr
  cmulv(b22r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  //  r B32 B32 dr
  cmulv(b32r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  // *******
  //  E_33
  // *******
  
  //  B33 dr
  for (i=0;i<ndofe;i++){
    gm[2][i]+=b22[i];
  }
  
  //  r B13 B13 dr
  cmulv(b13r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B23 B23 dr
  cmulv(b23r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B33 B33 dr
  cmulv(b33r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  

  // **************
  //  E_23 = E_32
  // **************
  
  //  (B23 + B32) dr
  for (i=0;i<ndofe;i++){
    gm[3][i]+=b23[i]+b32[i];
  }
  
  //  r B13 B12 dr
  cmulv(b13r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B12 B13 dr
  cmulv(b12r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }

  //  r B23 B22 dr
  cmulv(b23r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B22 B23 dr
  cmulv(b22r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B33 B32 dr
  cmulv(b33r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B32 B33 dr
  cmulv(b32r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  // **************
  //  E_31 = E_13
  // **************
  
  //  (B31 + B13) dr
  for (i=0;i<ndofe;i++){
    gm[4][i]+=b31[i]+b13[i];
  }
  
  //  r B11 B13 dr
  cmulv(b11r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B13 B11 dr
  cmulv(b13r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }

  //  r B21 B23 dr
  cmulv(b21r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B23 B21 dr
  cmulv(b23r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B31 B33 dr
  cmulv(b31r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B33 B31 dr
  cmulv(b33r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }


  // **************
  //  E_12 = E_21
  // **************
  
  //  (B12 + B21) dr
  for (i=0;i<ndofe;i++){
    gm[5][i]+=b12[i]+b21[i];
  }
  
  //  r B12 B11 dr
  cmulv(b12r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B11 B12 dr
  cmulv(b11r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }

  //  r B22 B21 dr
  cmulv(b22r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B21 B22 dr
  cmulv(b21r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B32 B31 dr
  cmulv(b32r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B31 B32 dr
  cmulv(b31r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  


}


/**
   function computes gradient %matrix for geometrically nonlinear problems
   
   @param grm - gradient %matrix
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian

   JK, 24.9.2005
*/
void linhexrot::gnl_grmatrix (matrix &grm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe);
  
  bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
  
  for (i=0;i<ndofe;i++){
    grm[0][i]=b11[i];
    grm[1][i]=b12[i];
    grm[2][i]=b13[i];
    grm[3][i]=b21[i];
    grm[4][i]=b22[i];
    grm[5][i]=b23[i];
    grm[6][i]=b31[i];
    grm[7][i]=b32[i];
    grm[8][i]=b33[i];
  }
}


void linhexrot::dmatblock (matrix &dd, matrix &d)
{
  dd[0][0] += d[3][3]/intordsm[0][0];  dd[1][1] += d[4][4]/intordsm[0][0];  dd[2][2] += d[5][5]/intordsm[0][0];

}


/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - nodes of element
   @param tmat - transformation %matrix
   
   JK
*/
void linhexrot::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,n,m;

  fillm (0.0,tmat);

  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*3+0][i*3]=Mt->nodes[nodes[i]].e1[0];
      tmat[i*3+1][i*3]=Mt->nodes[nodes[i]].e1[1];
      tmat[i*3+2][i*3]=Mt->nodes[nodes[i]].e1[2];

      tmat[i*3+0][i*3+1]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*3+1][i*3+1]=Mt->nodes[nodes[i]].e2[1];
      tmat[i*3+2][i*3+1]=Mt->nodes[nodes[i]].e2[2];

      tmat[i*3+0][i*3+2]=Mt->nodes[nodes[i]].e3[0];
      tmat[i*3+1][i*3+2]=Mt->nodes[nodes[i]].e3[1];
      tmat[i*3+2][i*3+2]=Mt->nodes[nodes[i]].e3[2];
    }
  }
}

/**
   function computes stiffness %matrix of one element

   function computes stiffness %matrix for geometrically linear problems

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   JK, 19.7.2001
*/
void linhexrot::gl_stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ipp, ipshear;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  matrix gm,d(tncomp,tncomp), dd(3,3);
  
  Mt->give_node_coord3d (x,y,z,eid);
  fillm (0.0,sm);
  fillm (0.0,dd);

  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  reallocm (ncomp[0],ndofe,gm);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	//  stiffness matrix of the material
	Mm->matstiff (d,ipp);
	ipp++;    
	//  geometric matrix
	geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	//  contribution to the stiffness matrix of the element
	jac*=w[i]*w[j]*w[k];
	if (jac<0.0)
          print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	bdbjac (sm,gm,d,gm,jac);

	// average for shear	
      dmatblock ( dd, d);
      }
    }
  }
 
  // Shear
  ipshear=2;
  reallocv (ipshear,w);
  reallocv (ipshear,gp);
  reallocm (3,ndofe,gm);
  gauss_points (gp.a,w.a,ipshear);
  for (i=0;i<ipshear;i++){
    xi=gp[i];
    for (j=0;j<ipshear;j++){
      eta=gp[j];
      for (k=0;k<ipshear;k++){
	zeta=gp[k];
	//  geometric matrix
	geom_matrix_shear (gm,x,y,z,xi,eta,zeta,jac);
	//  contribution to the stiffness matrix of the element
	jac*=w[i]*w[j]*w[k];
    bdbjac (sm,gm,dd,gm,jac);
	
      }
    }
  }
/*
  fprintf (Out,"\n\n SMloc");
    for (i=0;i<8;i++){
      fprintf (Out,"\n");
      for (jj=0;jj<6;jj++){
	  ii=i*6+jj;
      fprintf (Out," %le",sm[ii][ii]);
	  }
      fprintf (Out,"\n");
	}
 */
}

/**
   function computes stiffness %matrix of hexahedral finite element
   
   function computes stiffness %matrix for geometrically nonlinear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   
   JK, 21.9.2005
*/
void linhexrot::gnl_stiffness_matrix (long lcid,long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac,jac2;
  vector w,gp,sig(tncomp),r(ndofe),x(nne),y(nne),z(nne);
  matrix gm(tncomp,ndofe),grm(9,ndofe),d(tncomp,tncomp),s(9,9);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  component setting to zero
  fillm (0.0,sm);
  
  //  array for weights of integration points
  reallocv (intordsm[0][0],w);
  //  array for coordinates of integration points
  reallocv (intordsm[0][0],gp);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	//
	//  linear stiffness matrix and inital deformation matrix
	//
	
	//  strain-displacement matrix
	gngeom_matrix (gm,r,x,y,z,xi,eta,zeta,jac);
	
	//  stiffness matrix of the material
	Mm->matstiff (d,ipp);
	
	jac*=w[i]*w[j]*w[k];
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,gm,d,gm,jac);
	
	
	//
	//  initial stress matrix
	//
	
	//  gradient matrix
	gnl_grmatrix (grm,x,y,z,xi,eta,zeta,jac2);
	
	//  stresses
	Mm->givestress (lcid,ipp,sig);
	
	s[0][0]=sig[0];  s[0][1]=sig[5];  s[0][2]=sig[4];
	s[1][0]=sig[5];  s[1][1]=sig[1];  s[1][2]=sig[3];
	s[2][0]=sig[4];  s[2][1]=sig[3];  s[2][2]=sig[2];
	
	s[3][3]=sig[0];  s[3][4]=sig[5];  s[3][5]=sig[4];
	s[4][3]=sig[5];  s[4][4]=sig[1];  s[4][5]=sig[3];
	s[5][3]=sig[4];  s[5][4]=sig[3];  s[5][5]=sig[2];
	
	s[6][6]=sig[0];  s[6][7]=sig[5];  s[6][8]=sig[4];
	s[7][6]=sig[5];  s[7][7]=sig[1];  s[7][8]=sig[3];
	s[8][6]=sig[4];  s[8][7]=sig[3];  s[8][8]=sig[2];
	
	
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,grm,s,grm,jac);
	
	ipp++;
      }
    }
  }
}


/**
   function assembles resulting stiffness %matrix of the element
   
   @param lcid - load case id
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 9.5.2002
*/
void linhexrot::res_stiffness_matrix (long /*lcid*/,long eid,matrix &sm)
{
  gl_stiffness_matrix (eid,0,0,sm);
  //  transformation of stiffness matrix
  ivector nodes (nne);
  Mt->give_elemnodes (eid,nodes);
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  //gnl_stiffness_matrix (lcid,eid,0,0,sm);
}

/**
   function computes mass %matrix

   @param eid - number of element
   @param mm - mass %matrix

   JK, 19.7.2001
*/
void linhexrot::mass_matrix (long eid,matrix &mm)
{
  long i,j,k;
  double jac,xi,eta,zeta,rho;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);

  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);
	
	jac=fabs(jac);

	bf_matrix (n,xi,eta,zeta);

	rho = approx (xi,eta,zeta,dens);
	
	jac*=w[i]*w[j]*w[k]*rho;
	
	nnj (mm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}

/**
   function computes mass %matrix

   @param eid - number of element
   @param mm - mass %matrix

   JK, 19.7.2001
*/
void linhexrot::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  
  mass_matrix (eid,mm);

  if (Mp->diagmass==1){
    diagonalization (mm);
  }

  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  ivector nodes (nne);
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}

/**
   function computes load matrix

   @param eid - number of element
   @param lm - load matrix
   
   25.7.2001
*/
void linhexrot::load_matrix (long eid,matrix &lm)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];  w3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);

	bf_matrix (n,xi,eta,zeta);

	jac*=w1*w2*w3;
	
	nnj (lm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}

/**
   function computes load matrix

   @param eid - number of element
   @param lm - load matrix
   
   25.7.2001
*/
void linhexrot::res_load_matrix (long eid,matrix &lm)
{
  long transf;
  ivector nodes(nne);

  load_matrix (eid,lm);
  
  //  transformation of load matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }
}




/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, modified 23.11.2006
*/
void linhexrot::res_ip_strains (long lcid,long eid)
{
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),r(ndofe),aux;
  matrix tmat;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  gl_ip_strains (lcid,eid,0,0,x,y,z,r);
  //gnl_ip_strains (lcid,eid,0,0,x,y,z,r);
  
}


/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y,z - %vectors of nodal coordinates
   @param r - %vector of nodal displacements
   
   10.5.2002, JK
*/
void linhexrot::gl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,jac;
  vector gp,w,eps;
  matrix gm;
  
  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  
	  mxv (gm,r,eps);
	  
	  Mm->storestrain (lcid,ipp,cncomp[ii],eps);
	  
	  ipp++;
	}
      }
    }
  }
}

/**
   function computes strains at integration points of element
   function is used in geometrically linear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ii - number of block
   @param x,y,z - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 24.9.2005
*/
void linhexrot::gnl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac,b11r,b12r,b13r,b21r,b22r,b23r,b31r,b32r,b33r;
  vector gp,w,eps;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe);
  
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  reallocv (tncomp,eps);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
	
	scprd (b11,r,b11r);
	scprd (b12,r,b12r);
	scprd (b13,r,b13r);
	scprd (b21,r,b21r);
	scprd (b22,r,b22r);
	scprd (b23,r,b23r);
	scprd (b31,r,b31r);
	scprd (b32,r,b32r);
	scprd (b33,r,b33r);
	
	eps[0] = b11r + 0.5*b11r*b11r + 0.5*b21r*b21r + 0.5*b31r*b31r;
	eps[1] = b22r + 0.5*b12r*b12r + 0.5*b22r*b22r + 0.5*b32r*b32r;
	eps[2] = b33r + 0.5*b13r*b13r + 0.5*b23r*b23r + 0.5*b33r*b33r;

	eps[3] = b23r+b32r + b12r*b13r + b22r*b23r + b32r*b33r;
	eps[4] = b31r+b13r + b13r*b11r + b23r*b21r + b33r*b31r;
	eps[5] = b12r+b21r + b11r*b12r + b21r*b22r + b31r*b32r;
	
	
	Mm->storestrain (lcid,ipp,eps);
	ipp++;
      }
    }
  }
}



/**
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 10.5.2002
*/
void linhexrot::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
//  nodip_linhexrot (ipp,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain (lcid,ipnum[i],eps);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storestrain (lcid,0,eps);
  }
}

/**
   function computes nodal strains directly
   
   @param lcid - load case id
   @param eid - element id
   @param stra - array for strain components
   
   stra[i][j] - the j-th strain component at the i-th node
   
   JK, 26.9.2004
*/
void linhexrot::nod_strains_comp (long lcid,long eid,double **stra)
{
  long i,j;
  double jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),nxi(nne),neta(nne),nzeta(nne),r(ndofe),eps(tncomp),aux;
  matrix tmat,gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
//  nodcoord_linhexrot (nxi,neta,nzeta);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,z,nxi[i],neta[i],nzeta[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    
    for (j=0;j<eps.n;j++){
      stra[i][j]=eps[j];
    }
  }

}



/**
   function computes strains in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi, eta, zeta - natural coordinates of the point
   @param fi - first index
   @param ncomp - number of components
   @param eps - array containing strains
   
   11.5.2002
*/
void linhexrot::appstrain (long lcid,long eid,double xi,double eta,double zeta,long fi,long ncomp,vector &eps)
{
  long i,j,k;
  ivector nod(nne);
  vector nodval(nne);
  
  if (ncomp != eps.n){
    fprintf (stderr,"\n\n wrong interval of indices in function strain (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  Mt->give_elemnodes (eid,nod);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nod[j]].strain[lcid*tncomp+i];
    }
    eps[k]=approx (xi,eta,zeta,nodval);
    k++;
  }
}



void linhexrot::strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
  double **stra = NULL;
  vector coord,eps;
  
  if (Mp->strainaver==0){
    stra = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
    }
    //elem_strains (stra,lcid,eid,ri,ci);
  }
  
  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_strains (stra,lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (ncp,eps);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	appstrain (lcid,eid,coord[0],coord[1],coord[2],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlq::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  if (Mp->strainaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }
}

/**
   function computes stresses at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002, JK
*/
void linhexrot::res_ip_stresses (long lcid,long eid)
{
  ip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses at integration points of element
   stresses are computed by material models
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002, JK
*/
void linhexrot::ip_stresses (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	
	ipp++;
      }
    }
  }
}

/**
   function computes elastic stresses at integration points of element
   stresses are computed from strains with the help of elastic stiffness
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002, JK
*/
void linhexrot::ip_elast_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  vector gp,w,eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);

  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	Mm->matstiff (d,ipp);
	
	Mm->givestrain (lcid,ipp,eps);
	
	mxv (d,eps,sig);
	
	Mm->storestress (lcid,ipp,sig);
	
	ipp++;
      }
    }
  }
}

/**
   function computes stresses at nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002, JK
*/
void linhexrot::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
//  nodip_linhexrot (ipp,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    
    //  storage of stresses to the node
    j=nod[i];
    Mt->nodes[j].storestress (lcid,0,sig);
  }
}


/**
   function computes stresses in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi, eta, zeta - natural coordinates of the point
   @param fi,li - first and last indices
   @param sig - array containing stresses
   
   11.5.2002
*/
void linhexrot::appstress (long lcid,long eid,double xi,double eta,double zeta,long fi,long ncomp,vector &sig)
{
  long i,j,k;
  ivector nodes(nne);
  vector nodval(nne);
  
  if (ncomp != sig.n){
    fprintf (stderr,"\n\n wrong interval of indices in function stress (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].stress[lcid*tncomp+i];
    }
    sig[k]=approx (xi,eta,zeta,nodval);
    k++;
  }
}



void linhexrot::stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
//  double **stra,**stre;
  vector coord,sig;
  
  /*
  if (Mp->stressaver==0){
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    //elem_strains (stra,lcid,eid,ri,ci);
    elem_stresses (stra,stre,lcid,eid,ri,ci);
  }

  */

  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (stre,lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (ncp,sig);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],coord[2],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  /*
  if (Mp->stressaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
      delete [] stre[i];
    }
    delete [] stra;
    delete [] stre;
  }
  */

  //mainip_stresses (0,eid,0,0);
  
  
}

/**
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 24.10.2005
*/
void linhexrot::nod_eqother_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp,ncompo;
  ivector ipnum(nne),nod(nne);
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
//  nodip_linhexrot (ipp,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  number of components of array other
    ncompo = Mm->givencompeqother (ipnum[i],0);
    reallocv (ncompo,eqother);
    //  components of array other at integration points
    Mm->giveeqother (ipnum[i],0,ncompo,eqother.a);
    
    //  storage of components to the node
    j=nod[i];
    Mt->nodes[j].storeother (lcid,0,ncompo,eqother);
    
  }
}














/**
   function computes internal forces in the case of geometrical linear computation
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates
   
   JK, 28.7.2001
   TKo 7.2008
*/
void linhexrot::gl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x, vector &y, vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}


/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   JK, 24.9.2005
*/
void linhexrot::gnl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector w,gp,x(nne),y(nne),z(nne),sig(tncomp),contr(ndofe),r(ndofe);
  matrix gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  
  fillv (0.0,ifor);
  
  //  array for coordinates of integration points
  reallocv (intordsm[0][0],gp);
  //  array for weights of integration points
  reallocv (intordsm[0][0],w);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	//  computation of stress
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	
	if (Mp->nodeintfor==1){
	  Mm->givestress (lcid,ipp,sig);
	  
	  //  strain-displacement (geometric) matrix
	  gngeom_matrix (gm,r,x,y,z,xi,eta,zeta,jac);
	  
	  mtxv (gm,sig,contr);
	  
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  addv (contr,ifor,ifor);
	}

	ipp++;
      }
    }
  }
}


/**
   function computes nonlocal internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y,z - nodal coordinates
   
   JK, 28.7.2001
   TKo 7.2008
*/
void linhexrot::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of nonlocal stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);

}



/**
   function computes increments of internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void linhexrot::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes contributions from eigenstrains to the right hand side
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   JK, 17.8.2004
*/
void linhexrot::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=eigstress;
  
  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor,x,y,z);
}



/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 24.9.2005
*/
void linhexrot::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  gl_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  //gnl_internal_forces (lcid,eid,0,0,ifor);
  
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //globloctransf (ifor,v,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   function computes resulting internal forces for nonlocal models
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void linhexrot::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //globloctransf (ifor,v,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



void linhexrot::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  incr_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //globloctransf (ifor,v,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - %vector of internal forces

   JK, 17.8.2004
*/
void linhexrot::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);

  Mt->give_node_coord3d (x,y,z,eid);

  eigstrain_forces (lcid,eid,0,0,nfor,x,y,z);

  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //globloctransf (nfor,v,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void linhexrot::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	
	ipp++;
      }
    }
  }
}



/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo  7.2008
*/
void linhexrot::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstressesincr (ipp);
	
	ipp++;
      }
    }
  }
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void linhexrot::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	
	ipp++;
      }
    }
  }
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void linhexrot::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->compnonloc_nlstresses (ipp);
	ipp++;
      }
    }
  }
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void linhexrot::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
        //
        //  eigenstresses are computed from eigenstrains \sigma_0 = D (-\eps_0)
        //
        // restore eigenstrains
	Mm->giveeigstrain (ipp,eigstr);
        // change sign of eigenstrain vector
        chsgnv(eigstr);
	//  matrix of stiffness of the material
	Mm->matstiff (d,ipp);
        // calculate eigenstresses    
	mxv (d,eigstr,sig);
	Mm->storeeigstress (ipp,sig);
	ipp++;
      }
    }
  }
}



/**
   function integrates selected quantity over the finite element
   it results in nodal values
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y,z - node coordinates
   
   JK, 27.11.2006
   TKo 7.2008
*/
void linhexrot::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x, vector &y, vector &z)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector w,gp,ipv(tncomp),contr(ndofe);
  matrix gm(tncomp,ndofe);

  fillv (0.0,nv);
  
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
        //  function assembles required quantity at integration point
        Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);
	//  strain-displacement (geometric) matrix
	geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	//  contribution to the nodal values
	mtxv (gm,ipv,contr);
	cmulv (jac*w[i]*w[j]*w[k],contr);
	//  summation
	addv (contr,nv,nv);
	ipp++;
      }
    }
  }
}









/**
   The function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param coord - %vector with global coordinates of integration point (ouput)
   
   @return The function returns global coordinates in the argument coord.
   
   10.1.2002
*/
void linhexrot::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,k,ii;
  double xi,eta,zeta;
  vector x(nne),y(nne),z(nne),w(intordsm[ri][ci]),gp(intordsm[ri][ci]);
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordsm[ri][ci];k++){
	zeta=gp[k];
	if (ii==ipp){
	  coord[0]=approx (xi,eta,zeta,x);
	  coord[1]=approx (xi,eta,zeta,y);
	  coord[2]=approx (xi,eta,zeta,z);
	}
	ii++;
      }
    }
  }
}



/**
   The function assembles natural coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ncoord - %vector with natural coordinates of integration point (ouput)
   
   @return The function returns natural coordinates in the argument ncoord.

   Created by TKo, 12.2016
*/
void linhexrot::ipncoord (long eid,long ipp,vector &ncoord)
{
  long i, j, k, ii, ri, ci;
  double xi, eta, zeta;
  vector w, gp;
  
  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp));
      gauss_points (gp.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];
  
      for (i=0;i<intordsm[ri][ci];i++){
        xi=gp[i];
        for (j=0;j<intordsm[ri][ci];j++){
          eta=gp[j];
          for (k=0;k<intordsm[ri][ci];k++){
            zeta=gp[k];
            if (ii==ipp){
              ncoord[0]=xi;
              ncoord[1]=eta;
              ncoord[2]=zeta;
            }
            ii++;
          }
        }
      }
    }
  }
}



/**
  The function computes initial values of the given quantities at each integration point of the
  element from the nodal values given by the parameter nodval. Initial condition types must be 
  the same for all nodes of the element.

  @param eid - element id
  @param ri  - block row index
  @param ci  - block column index
  @param nodval - nodal values of particular initial conditions.
                  nodval[i][j] represents value of j-th initial condition at i-th node of the given element.
  @param ictn - array of types of initial condition for each node of element.
                The type of initial condition determines which values are being specified in the node. 
                (ictn[i] & inistrain) returns nonzero if nodal values of initial strains are specified
                (ictn[i] & inistress) returns nonzero if nodal values of initial stresses are specified
                (ictn[i] & iniother)  returns nonzero if nodal values of initial values of eqother array are specified
                (ictn[i] & inicond)   returns nonzero if nodal values of other initial conditions are specified

  @retur The function does not return anything.

  Created by Tomas Koudelka 2004
  Revised by Tomas Koudelka 03.2012
*/
void linhexrot::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, m, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, zeta, ipval;
  vector w, gp, anv(nne);
  long nstra, nstre, ncompstr, ncompeqother;
  long idstra, idstre, idoth, idic;
  inictype ict;
  int aux;

  nstra = idstra = nstre = idstre = idoth = idic = 0;

  ict = ictn[0];
  for (i=0; i<nne; i++)
  {
    aux = int(ictn[i])-int(ict);
    if (aux < 0)  aux = -aux;    
    aux &= ~(inidisp);
    aux &= ~(inidisp_x);
    aux &= ~(inidisp_y);
    aux &= ~(inidisp_z);
    if ((ictn[i] != ict) && aux)
    {
      print_err("Incompatible types of initial conditions on element %ld\n"
                " at %ld. and %ld. nodes", __FILE__, __LINE__, __func__, eid+1, 1, i+1);
      abort();
    }
  }
  for (j = 0; j < nv; j++) // for all initial values
  {
    for(i = 0; i < nne; i++) // for all nodes on element
      anv[i] = nodval[i][j];
    for (ii = 0; ii < nb; ii++)
    {
      for (jj = 0; jj < nb; jj++)
      {
        ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
        if (intordsm[ii][jj] == 0)
          continue;
        reallocv (intordsm[ii][jj],gp);
        reallocv (intordsm[ii][jj],w);
        gauss_points (gp.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp[k];
          for (l = 0; l < intordsm[ii][jj]; l++)
          {
            eta=gp[l];
            for (m = 0; m < intordsm[ii][jj]; m++)
            {
              zeta=gp[m];
              //  value in integration point
              ipval = approx (xi,eta,zeta,anv);
              ncompstr =  Mm->ip[ipp].ncompstr;
              ncompeqother = Mm->ip[ipp].ncompeqother;
              if ((ictn[0] & inistrain) && (j < ncompstr))
              {
                Mm->ip[ipp].strain[idstra] += ipval;
                ipp++;
                continue;
              }
              if ((ictn[0] & inistress) && (j < nstra + ncompstr))
              {
                Mm->ip[ipp].stress[idstre] += ipval;
                ipp++;
                continue;
              }
              if ((ictn[0] & iniother) && (j < nstra+nstre+ncompeqother))
              {
                Mm->ip[ipp].eqother[idoth] += ipval;
                ipp++;
                continue;
              }
              if ((ictn[0] & inicond) && (j < nv))
              {
                if (Mm->ic[ipp] == NULL)
                {
                  Mm->ic[ipp] = new double[nv-j];
                  memset(Mm->ic[ipp], 0, sizeof(*Mm->ic[ipp])*(nv-j)); 
                }
                Mm->ic[ipp][idic] += ipval;
                ipp++;
                continue;
              }
              ipp++;
            }
          }
        }
      }
    }
    ipp=Mt->elements[eid].ipp[ri][ci];
    ncompstr =  Mm->ip[ipp].ncompstr;
    ncompeqother = Mm->ip[ipp].ncompeqother;
    if ((ictn[0] & inistrain) && (j < ncompstr))
    {
      nstra++;
      idstra++;
      continue;
    }
    if ((ictn[0] & inistress) && (j < nstra + ncompstr))
    {      
      nstre++;
      idstre++;
      continue;
    }  
    if ((ictn[0] & iniother)  && (j < nstra + nstre + ncompeqother))
    {
      idoth++;
      continue;
    }
    if ((ictn[0] & inicond) && (j < nv))
    {
      idic++;
      continue;
    }
  }
}



/**
   function computes volume appropriate to integration point
   
   2.3.2004, JK
*/
void linhexrot::ipvolume (long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    jac_3d (jac,x,y,z,xi,eta,zeta);
	    jac=fabs(jac);
	    
	    jac*=w[i]*w[j]*w[k];
	    
	    Mm->storeipvol (ipp,jac);
	    ipp++;
	  }
	}
      }
    }
  }
  
}

/**
   function computes nodal forces caused by presure on surface
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - vector of presure 
   @param eis - surface id 
   
   27.1.2006
*/
void linhexrot::node_forces_surf (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i,j;
  double xi,eta,zeta,jac;
  double *tnv;
  vector x(nne),y(nne),z(nne),gp,w,av(ndofe),v(ndofe);
  matrix n(napfun,ndofe),an(napfun,ndofe),am(ndofe,ndofe),tran(3,3);
  
  tnv = new double [12];

  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  reallocv (intordb,w);
  reallocv (intordb,gp);
  gauss_points (gp.a,w.a,intordb);
  
  //  surface number 1

  if (is[0]>0 ){
    xi=1.0;
    for (i=0;i<intordb;i++){
      eta=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,eta,zeta,0);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[0]==1){
      av[0] = nv[0*3+0];
      av[1] = nv[0*3+1];
      av[2] = nv[0*3+2];
      
      av[9] = nv[1*3+0];
      av[10] = nv[1*3+1];
      av[11] = nv[1*3+2];
      
      av[21] = nv[2*3+0];
      av[22] = nv[2*3+1];
      av[23] = nv[2*3+2];

      av[12] = nv[3*3+0];
      av[13] = nv[3*3+1];
      av[14] = nv[3*3+2];
    }


    if (is[0]==2){
      
      av[0] = nv[0*3+0];
      av[1] = nv[0*3+1];
      av[2] = nv[0*3+2];
      
      av[3] = nv[1*3+0];
      av[4] = nv[1*3+1];
      av[5] = nv[1*3+2];
      
      av[6] = nv[2*3+0];
      av[7] = nv[2*3+1];
      av[8] = nv[2*3+2];

      av[9]  = nv[3*3+0];
      av[10] = nv[3*3+1];
      av[11] = nv[3*3+2];
      
      locglob_nodeval (0,av,tnv,x,y,z);
      
      av[0] = tnv[0*3+0];
      av[1] = tnv[0*3+1];
      av[2] = tnv[0*3+2];
      
      av[9] = tnv[1*3+0];
      av[10] = tnv[1*3+1];
      av[11] = tnv[1*3+2];
      
      av[21] = tnv[2*3+0];
      av[22] = tnv[2*3+1];
      av[23] = tnv[2*3+2];

      av[12] = tnv[3*3+0];
      av[13] = tnv[3*3+1];
      av[14] = tnv[3*3+2];
    }


    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 2

  if (is[1]>0 ){
    eta=1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,zeta,1);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[1]==1){
      av[3] = nv[12+0*3+0];
      av[4] = nv[12+0*3+1];
      av[5] = nv[12+0*3+2];
      
      av[0] = nv[12+1*3+0];
      av[1] = nv[12+1*3+1];
      av[2] = nv[12+1*3+2];
      
      av[12] = nv[12+2*3+0];
      av[13] = nv[12+2*3+1];
      av[14] = nv[12+2*3+2];

      av[15] = nv[12+3*3+0];
      av[16] = nv[12+3*3+1];
      av[17] = nv[12+3*3+2];
    }


    if (is[1]==2){
      
      av[0] = nv[12+0*3+0];
      av[1] = nv[12+0*3+1];
      av[2] = nv[12+0*3+2];
      
      av[3] = nv[12+1*3+0];
      av[4] = nv[12+1*3+1];
      av[5] = nv[12+1*3+2];
      
      av[6] = nv[12+2*3+0];
      av[7] = nv[12+2*3+1];
      av[8] = nv[12+2*3+2];

      av[9]  = nv[12+3*3+0];
      av[10] = nv[12+3*3+1];
      av[11] = nv[12+3*3+2];
      
      locglob_nodeval (1,av,tnv,x,y,z);
      
      av[3] = tnv[0*3+0];
      av[4] = tnv[0*3+1];
      av[5] = tnv[0*3+2];
      
      av[0] = tnv[1*3+0];
      av[1] = tnv[1*3+1];
      av[2] = tnv[1*3+2];
      
      av[12] = tnv[2*3+0];
      av[13] = tnv[2*3+1];
      av[14] = tnv[2*3+2];

      av[15] = tnv[3*3+0];
      av[16] = tnv[3*3+1];
      av[17] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }

  //  surface number 3

  if (is[2]>0 ){
    xi=-1.0;
    for (i=0;i<intordb;i++){
      eta=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,eta,zeta,2);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[2]==1){
      av[6] = nv[24+0*3+0];
      av[7] = nv[24+0*3+1];
      av[8] = nv[24+0*3+2];
      
      av[3] = nv[24+1*3+0];
      av[4] = nv[24+1*3+1];
      av[5] = nv[24+1*3+2];
      
      av[15] = nv[24+2*3+0];
      av[16] = nv[24+2*3+1];
      av[17] = nv[24+2*3+2];

      av[18] = nv[24+3*3+0];
      av[19] = nv[24+3*3+1];
      av[20] = nv[24+3*3+2];
    }


    if (is[2]==2){
      
      av[0] = nv[24+0*3+0];
      av[1] = nv[24+0*3+1];
      av[2] = nv[24+0*3+2];
      
      av[3] = nv[24+1*3+0];
      av[4] = nv[24+1*3+1];
      av[5] = nv[24+1*3+2];
      
      av[6] = nv[24+2*3+0];
      av[7] = nv[24+2*3+1];
      av[8] = nv[24+2*3+2];

      av[9]  = nv[24+3*3+0];
      av[10] = nv[24+3*3+1];
      av[11] = nv[24+3*3+2];
      
      locglob_nodeval (2,av,tnv,x,y,z);
      
      av[6] = tnv[0*3+0];
      av[7] = tnv[0*3+1];
      av[8] = tnv[0*3+2];
      
      av[3] = tnv[1*3+0];
      av[4] = tnv[1*3+1];
      av[5] = tnv[1*3+2];
      
      av[15] = tnv[2*3+0];
      av[16] = tnv[2*3+1];
      av[17] = tnv[2*3+2];

      av[18] = tnv[3*3+0];
      av[19] = tnv[3*3+1];
      av[20] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 4

  if (is[3]>0 ){
    eta=-1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,zeta,3);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[3]==1){
      av[9] = nv[36+0*3+0];
      av[10] = nv[36+0*3+1];
      av[11] = nv[36+0*3+2];
      
      av[6] = nv[36+1*3+0];
      av[7] = nv[36+1*3+1];
      av[8] = nv[36+1*3+2];
      
      av[18] = nv[36+2*3+0];
      av[19] = nv[36+2*3+1];
      av[20] = nv[36+2*3+2];

      av[21] = nv[36+3*3+0];
      av[22] = nv[36+3*3+1];
      av[23] = nv[36+3*3+2];
    }


    if (is[3]==2){
      
      av[0] = nv[36+0*3+0];
      av[1] = nv[36+0*3+1];
      av[2] = nv[36+0*3+2];
      
      av[3] = nv[36+1*3+0];
      av[4] = nv[36+1*3+1];
      av[5] = nv[36+1*3+2];
      
      av[6] = nv[36+2*3+0];
      av[7] = nv[36+2*3+1];
      av[8] = nv[36+2*3+2];

      av[9]  = nv[36+3*3+0];
      av[10] = nv[36+3*3+1];
      av[11] = nv[36+3*3+2];
      
      locglob_nodeval (3,av,tnv,x,y,z);
      
      av[9] = tnv[0*3+0];
      av[10]= tnv[0*3+1];
      av[11]= tnv[0*3+2];
      
      av[6] = tnv[1*3+0];
      av[7] = tnv[1*3+1];
      av[8] = tnv[1*3+2];
      
      av[18] = tnv[2*3+0];
      av[19] = tnv[2*3+1];
      av[20] = tnv[2*3+2];

      av[21] = tnv[3*3+0];
      av[22] = tnv[3*3+1];
      av[23] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }

  //  surface number 5

  if (is[4]>0 ){
    zeta=1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	eta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,eta,4);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[4]==1){
      av[0] = nv[48+0*3+0];
      av[1] = nv[48+0*3+1];
      av[2] = nv[48+0*3+2];
      
      av[3] = nv[48+1*3+0];
      av[4] = nv[48+1*3+1];
      av[5] = nv[48+1*3+2];
      
      av[6] = nv[48+2*3+0];
      av[7] = nv[48+2*3+1];
      av[8] = nv[48+2*3+2];

      av[9] = nv[48+3*3+0];
      av[10] = nv[48+3*3+1];
      av[11] = nv[48+3*3+2];
    }


    if (is[4]==2){
      
      av[0] = nv[48+0*3+0];
      av[1] = nv[48+0*3+1];
      av[2] = nv[48+0*3+2];
      
      av[3] = nv[48+1*3+0];
      av[4] = nv[48+1*3+1];
      av[5] = nv[48+1*3+2];
      
      av[6] = nv[48+2*3+0];
      av[7] = nv[48+2*3+1];
      av[8] = nv[48+2*3+2];

      av[9]  = nv[48+3*3+0];
      av[10] = nv[48+3*3+1];
      av[11] = nv[48+3*3+2];
      
      locglob_nodeval (4,av,tnv,x,y,z);
      
      av[0] = tnv[0*3+0];
      av[1] = tnv[0*3+1];
      av[2] = tnv[0*3+2];
      
      av[3] = tnv[1*3+0];
      av[4] = tnv[1*3+1];
      av[5] = tnv[1*3+2];
      
      av[6] = tnv[2*3+0];
      av[7] = tnv[2*3+1];
      av[8] = tnv[2*3+2];

      av[9] = tnv[3*3+0];
      av[10] = tnv[3*3+1];
      av[11] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 6

  if (is[5]>0 ){
    zeta=-1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	eta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,eta,5);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[5]==1){
      av[12] = nv[60+0*3+0];
      av[13] = nv[60+0*3+1];
      av[14] = nv[60+0*3+2];
      
      av[21] = nv[60+1*3+0];
      av[22] = nv[60+1*3+1];
      av[23] = nv[60+1*3+2];
      
      av[18] = nv[60+2*3+0];
      av[19] = nv[60+2*3+1];
      av[20] = nv[60+2*3+2];

      av[15] = nv[60+3*3+0];
      av[16] = nv[60+3*3+1];
      av[17] = nv[60+3*3+2];
    }


    if (is[5]==2){
      
      av[0] = nv[60+0*3+0];
      av[1] = nv[60+0*3+1];
      av[2] = nv[60+0*3+2];
      
      av[3] = nv[60+1*3+0];
      av[4] = nv[60+1*3+1];
      av[5] = nv[60+1*3+2];
      
      av[6] = nv[60+2*3+0];
      av[7] = nv[60+2*3+1];
      av[8] = nv[60+2*3+2];

      av[9]  = nv[60+3*3+0];
      av[10] = nv[60+3*3+1];
      av[11] = nv[60+3*3+2];
      
      locglob_nodeval (5,av,tnv,x,y,z);
      
      av[12] = tnv[0*3+0];
      av[13] = tnv[0*3+1];
      av[14] = tnv[0*3+2];
      
      av[21] = tnv[1*3+0];
      av[22] = tnv[1*3+1];
      av[23] = tnv[1*3+2];
      
      av[18] = tnv[2*3+0];
      av[19] = tnv[2*3+1];
      av[20] = tnv[2*3+2];

      av[15] = tnv[3*3+0];
      av[16] = tnv[3*3+1];
      av[17] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }






  delete [] tnv;
}

void linhexrot::locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z)
{
  double norm;
  vector g1(3),g2(3),g3(3),lv(3),gv(3);
  matrix t(3,3);
  
  if (is==0){
    g1[0]=x[4]-x[7];
    g1[1]=y[4]-y[7];
    g1[2]=z[4]-z[7];
    
    g2[0]=x[3]-x[7];
    g2[1]=y[3]-y[7];
    g2[2]=z[3]-z[7];
  }
  if (is==1){
    g1[0]=x[5]-x[4];
    g1[1]=y[5]-y[4];
    g1[2]=z[5]-z[4];
    
    g2[0]=x[0]-x[4];
    g2[1]=y[0]-y[4];
    g2[2]=z[0]-z[4];
  }
  if (is==2){
    g1[0]=x[6]-x[5];
    g1[1]=y[6]-y[5];
    g1[2]=z[6]-z[5];
    
    g2[0]=x[1]-x[5];
    g2[1]=y[1]-y[5];
    g2[2]=z[1]-z[5];
  }
  if (is==3){
    g1[0]=x[7]-x[6];
    g1[1]=y[7]-y[6];
    g1[2]=z[7]-z[6];
    
    g2[0]=x[2]-x[6];
    g2[1]=y[2]-y[6];
    g2[2]=z[2]-z[6];
  }
  if (is==4){
    g1[0]=x[3]-x[2];
    g1[1]=y[3]-y[2];
    g1[2]=z[3]-z[2];
    
    g2[0]=x[1]-x[2];
    g2[1]=y[1]-y[2];
    g2[2]=z[1]-z[2];
  }
  if (is==5){
    g1[0]=x[5]-x[6];
    g1[1]=y[5]-y[6];
    g1[2]=z[5]-z[6];
    
    g2[0]=x[7]-x[6];
    g2[1]=y[7]-y[6];
    g2[2]=z[7]-z[6];
  }
  
  
  scprd (g1,g1,norm);
  norm=1.0/sqrt(norm);
  cmulv (norm,g1,g1);
  
  scprd (g1,g2,norm);
  g2[0]=g2[0]-norm*g1[0];
  g2[1]=g2[1]-norm*g1[1];
  g2[2]=g2[2]-norm*g1[2];
  
  scprd (g2,g2,norm);
  norm=1.0/sqrt(norm);
  cmulv (norm,g2,g2);
  
  g3[0]=g1[1]*g2[2]-g1[2]*g2[1];
  g3[1]=g1[2]*g2[0]-g1[0]*g2[2];
  g3[2]=g1[0]*g2[1]-g1[1]*g2[0];
  
  t[0][0]=g1[0];
  t[1][0]=g1[1];
  t[2][0]=g1[2];

  t[0][1]=g2[0];
  t[1][1]=g2[1];
  t[2][1]=g2[2];

  t[0][2]=g3[0];
  t[1][2]=g3[1];
  t[2][2]=g3[2];
  
  mxv (t.a,nv.a,tnv,3,3);
  mxv (t.a,nv.a+3,tnv+3,3,3);
  mxv (t.a,nv.a+6,tnv+6,3,3);
  mxv (t.a,nv.a+9,tnv+9,3,3);
  
}

/**
   function computes nodal forces caused by presure on surface
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - vector of presure 
   @param eis - surface id 
   
   4.2002, PF
*/
void linhexrot::node_forces_surf_old (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i,j,i1,ii;
  double xi=0.0,eta=0.0,zeta=0.0,jac, w1,w2;
  ivector nodes(nne);
  vector gx(nne),gy(nne),gz(nne),x(nnsurf),y(nnsurf),z(nnsurf),gp(intordb),w(intordb),av(ndofe),v(ndofe),pom(3);
  matrix n(napfun,ndofe),an(napfun,ndofe),am(ndofe,ndofe), tran(3,3);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (gx,gy,gz,eid);
  gauss_points (gp.a,w.a,intordb);
  fillv (0.0,nf);
  fillm (0.0,an);
  for (ii=0;ii<6;ii++){
    fillv (0.0,av);
    // is=0 not loading
    if (is[ii] ==0){}
    else {
      tran_mat( tran, gx, gy, gz, ii);
      // is=1 surface node 1,4,8,5
      if (ii ==0){
	xi=1.0;
        x[0]=gx[0];x[1]=gx[3];x[2]=gx[7];x[3]=gx[4]; 
        y[0]=gy[0];y[1]=gy[3];y[2]=gy[7];y[3]=gy[4]; 
        z[0]=gz[0];z[1]=gz[3];z[2]=gz[7];z[3]=gz[4]; 
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;     
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[i];  av[9+i]=nv[3+i];  av[21+i]=nv[6+i];  av[12+i]=nv[9+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=2 surface node 2,1,5,6
      else if (ii ==1){
	eta=1.0;
        x[0]=gx[1];x[1]=gx[0];x[2]=gx[4];x[3]=gx[5];
        y[0]=gy[1];y[1]=gy[0];y[2]=gy[4];y[3]=gy[5]; 
        z[0]=gz[1];z[1]=gz[0];z[2]=gz[4];z[3]=gz[5]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[3+i]=nv[12+i];  av[i]=nv[15+i];  av[12+i]=nv[18+i];  av[15+i]=nv[21+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=3 surface node 3,2,6,7
      else if (ii ==2){
	xi=-1.0;
        x[0]=gx[2];x[1]=gx[1];x[2]=gx[5];x[3]=gx[6]; 
        y[0]=gy[2];y[1]=gy[1];y[2]=gy[5];y[3]=gy[6]; 
        z[0]=gz[2];z[1]=gz[1];z[2]=gz[5];z[3]=gz[6]; 
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[6+i]=nv[24+i];  av[3+i]=nv[27+i];  av[15+i]=nv[30+i];  av[18+i]=nv[33+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=4 surface node 4,3,7,8
      else if (ii ==3){
	eta=-1.0;
        x[0]=gx[3];x[1]=gx[2];x[2]=gx[6];x[3]=gx[7]; 
        y[0]=gy[3];y[1]=gy[2];y[2]=gy[6];y[3]=gy[7];
        z[0]=gz[3];z[1]=gz[2];z[2]=gz[6];z[3]=gz[7]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;
	    
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[9+i]=nv[36+i];  av[6+i]=nv[39+i];  av[18+i]=nv[42+i];  av[21+i]=nv[45+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=5  surface node 1,2,3,4
      else if (ii ==4) {
	zeta=1.0;
        x[0]=gx[0];x[1]=gx[1];x[2]=gx[2];x[3]=gx[3]; 
        y[0]=gy[0];y[1]=gy[1];y[2]=gy[2];y[3]=gy[3]; 
        z[0]=gz[0];z[1]=gz[1];z[2]=gz[2];z[3]=gz[3]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    //	for constant			for (i1=0;i1<ndofe;i1++){
	    //								for (j1=0;j1<3;j1++){
	    //									an[j1][i1]=an[j1][i1]+n[j1][i1]*jac;}	}
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[48+i];  av[3+i]=nv[51+i];  av[6+i]=nv[54+i];  av[9+i]=nv[57+i];
	}
	//		fprintf (Out,"\n\n nn");
	//		for (i1=0;i1<ndofe;i1++){
	//			fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i1,an[0][i1],an[1][i1],an[2][i1]);}       
	
	// for constant		for (i=0;i<ndofe;i++){
	//						for (j=0;j<pom.n;j++){
	//							v[i]=an[j][i]*pom[j];}	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=6  surface node 6,5,8,7
      else if (ii ==5) {
	zeta=-1.0;
        x[0]=gx[5];x[1]=gx[4];x[2]=gx[7];x[3]=gx[6]; 
        y[0]=gy[5];y[1]=gy[4];y[2]=gy[7];y[3]=gy[6]; 
        z[0]=gz[5];z[1]=gz[4];z[2]=gz[7];z[3]=gz[6]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[12+i]=nv[60+i];  av[15+i]=nv[63+i];  av[18+i]=nv[66+i];  av[21+i]=nv[69+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransfblock (av, tran);
	}
        mxv (am,av,v);  
      }
      
			fprintf (Out,"\n\n zatiz");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,av[i1],av[i1+1],av[i1+2]);}       
			fprintf (Out,"\n\n vloc");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,v[i1],v[i1+1],v[i1+2]);}       
//      lgvectortransfblock(v, tran);
			fprintf (Out,"\n\n vglob");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,v[i1],v[i1+1],v[i1+2]);}       
      addv (nf,v,nf);      
    }
  }

  //  transformation of stiffness matrix to nodesystem
  //  transf = Mt->locsystems (nodes);
  //  if (transf>0){
  //    matrix tmat (ndofe,ndofe);
  //    transf_matrix (nodes,tmat);
  //    lgvectortransfblock (v,tmat);
  //  }

  //  fprintf (Out,"\n\n zatizeni na prvku cislo %ld",eid);
  //  for (i=0;i<ndofe;i++){
  //	  fprintf (Out,"\n %4ld   %20.10e",i,nf[i]);}
}

/**
   function computes transformation matrix on surface
   
   @param x, y - local coordinate  xL=adge12
   @param tran - tranformation metrix to GCS
   @param gx,gy,gz - vector of node global coordinate 
   @param is - surface id 
   
   4.2002, PF
*/
void linhexrot::tran_mat(matrix &tran, vector &gx, vector &gy, vector &gz, long is)
{
  long i,i1;
  double dl;
  matrix a(3,3);
// is=5  surface node 1,2,3,4
  if (is ==4) {
    for (i=0; i<3; i++) {
      i1=i+1; if(i1>2)i1=i1-3;
      a[0][i]=gx[i1]-gx[i];
      a[1][i]=gy[i1]-gy[i];
      a[2][i]=gz[i1]-gz[i];
    }
  }
// is=6  surface node 6,5,8,7
  else if (is ==5) {
    for (i=4; i<7; i++) {
      i1=i+1; if(i1>6)i1=i1-3;
      a[0][i-4]=gx[i1]-gx[i];
      a[1][i-4]=gy[i1]-gy[i];
      a[2][i-4]=gz[i1]-gz[i];
    }
  }
// is=1 surface node 1,4,8,5
// is=2 surface node 2,1,5,6
// is=3 surface node 3,2,6,7
// is=4 surface node 4,3,7,8
  else {
      if(is>0)i1=is-4;else i1=is;
      a[0][0]=gx[i1+3]-gx[is];
      a[1][0]=gy[i1+3]-gy[is];
      a[2][0]=gz[i1+3]-gz[is];
      a[0][1]=gx[i1+7]-gx[i1+3];
      a[1][1]=gy[i1+7]-gy[i1+3];
      a[2][1]=gz[i1+7]-gz[i1+3];
      a[0][2]=gx[is]-gx[i1+7];
      a[1][2]=gy[is]-gy[i1+7];
      a[2][2]=gz[is]-gz[i1+7];
  }
  
  dl=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);

  tran[0][0]=a[0][0]/dl;
  tran[1][0]=a[1][0]/dl;
  tran[2][0]=a[2][0]/dl;
  
  tran[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
  tran[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
  tran[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

  dl=sqrt(tran[0][2]*tran[0][2]+tran[1][2]*tran[1][2]+tran[2][2]*tran[2][2]);

  tran[0][2]=tran[0][2]/dl;
  tran[1][2]=tran[1][2]/dl;
  tran[2][2]=tran[2][2]/dl;
  
  tran[0][1]=tran[1][2]*tran[2][0]-tran[2][2]*tran[1][0];
  tran[1][1]=tran[2][2]*tran[0][0]-tran[0][2]*tran[2][0];
  tran[2][1]=tran[0][2]*tran[1][0]-tran[1][2]*tran[0][0];

  dl=sqrt(tran[0][1]*tran[0][1]+tran[1][1]*tran[1][1]+tran[2][1]*tran[2][1]);

  tran[0][1]=tran[0][1]/dl;
  tran[1][1]=tran[1][1]/dl;
  tran[2][1]=tran[2][1]/dl;
/*  
// local coordinate x=s12
    sl1[0]=0.0;    sl1[1]=0.0;    sl1[2]=0.0;
    sl2[0]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
    sl2[1]=0.0;    sl2[2]=0.0;
    sl3[0]=(s3(0)-s1(0))*tran(0,0)+(s3(1)-s1(1))*tran(1,0)+(s3(2)-s1(2))*tran(2,0);
    sl3[1]=(s3(0)-s1(0))*tran(0,1)+(s3(1)-s1(1))*tran(1,1)+(s3(2)-s1(2))*tran(2,1);
    sl3[2]=0.0;
*/
// local coordinate x=s12
//    x[0]=0.0;
//    y[0]=0.0;
//    x[1]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
//    y[1]=0.0;
//    x[2]=(gx(2)-gx(0))*tran(0,0)+(gy(2)-gy(0))*tran(1,0)+(gz(2)-gz(0))*tran(2,0);
//    y[2]=(gx(2)-gx(0))*tran(0,1)+(gy(2)-gy(0))*tran(1,1)+(gz(2)-gz(0))*tran(2,1);
//    x[3]=(gx(3)-gx(0))*tran(0,0)+(gy(3)-gy(0))*tran(1,0)+(gz(3)-gz(0))*tran(2,0);
//    y[3]=(gx(3)-gx(0))*tran(0,1)+(gy(3)-gy(0))*tran(1,1)+(gz(3)-gz(0))*tran(2,1);
}


/**
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   17.8.2004, JK
*/
void linhexrot::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k,l;
  double xi,eta,zeta;
  vector w,gp;
  
  l=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    ipval[l]=approx (xi,eta,zeta,nodval);
	    l++;
	  }
	}
      }
    }
  }
}



void linhexrot::aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),gp,w,eps;

  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  //Mm->givestrain (lcid,ipp,eps);
	  Mm->givestress (lcid,ipp,eps);
	  
         jac_3d (jac,x,y,z,xi,eta,zeta);	
         jac=fabs(jac);
	  
	  volume+=w[i]*w[j]*w[k]*jac;
	  
	  for (l=0;l<averstra.n;l++){
	    averstra[l]+=eps[l]*w[i]*w[j]*w[k]*jac;
	  }
	  
	  ipp++;
	}
      }
    }
  }
}





/**
   function picks up nodal values on required surface
   
   @param surf - number of required surface
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all surfaces
   
   JK, 19.8.2004
*/
void linhexrot::surfnodeval (long surf,vector &nodval,double *list)
{
  long i,j,k;
  ivector surfnod(nnsurf);
  
  fillv (0.0,nodval);
  quadhexahedral_surfnod (surfnod.a,surf);
  
  k=0;
  for (i=0;i<nnsurf;i++){
    for (j=0;j<napfun;j++){
      nodval[k]=list[surf*nnsurf*napfun+k];
      k++;
    }
  }
}
