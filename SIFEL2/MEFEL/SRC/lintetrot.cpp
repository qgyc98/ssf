#include "lintetrot.h"
#include "gadaptivity.h"
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
#include <stdlib.h>
#include <math.h>

lintetrot::lintetrot (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=24;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=1;
  //  number of edges on element
  ned=6;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=4;
  //  number of nodes on one surface
  nnsurf=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;
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
  
  nip[0][0]=4;  

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=4;
}

lintetrot::~lintetrot (void)
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
   
   @param volcoord - volume coordinates
   @param nodval - nodal values
   
   PF, 2.10.2008
*/
double lintetrot::approx (vector &volcoord,vector &nodval)
{
  double f;

  scprd (volcoord,nodval,f);
  
  return f;
}

/**
   function approximates function defined by nodal values
   
   @param volcoord - volume coordinates
   @param nodval - nodal values
   
   20.8.2001
*/
double lintetrot::approx_nat (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector volcoord(4);
  
  volcoord[0]=xi;
  volcoord[1]=eta;
  volcoord[2]=zeta;
  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
  
  scprd (volcoord,nodval,f);
  
  return f;
}

/**
   Transformation
   
   @param volcoord - volume coordinates
   @param nodval - nodal values
   
   20.8.2008
*/
void lintetrot::tran_side (matrix &t12, matrix &t13, matrix &t14, matrix &t23, matrix &t34, matrix &t42, vector &x,vector &y,vector &z)
{
  long i;
  double pom;
  vector an(3),am(3),dl(ned);
  matrix a(3,ned);
  
  //  in x,y,z  are global coord. node 1 , 2 , 3, 4
  //  side 1-2, 1-3, 1-4, 2-3, 3-4, 4-2
  a[0][0]=x[1]-x[0]; a[0][1]=x[2]-x[0]; a[0][2]=x[3]-x[0]; 
  a[0][3]=x[2]-x[1]; a[0][4]=x[3]-x[2]; a[0][5]=x[1]-x[3];
  a[1][0]=y[1]-y[0]; a[1][1]=y[2]-y[0]; a[1][2]=y[3]-y[0]; 
  a[1][3]=y[2]-y[1]; a[1][4]=y[3]-y[2]; a[1][5]=y[1]-y[3];
  a[2][0]=z[1]-z[0]; a[2][1]=z[2]-z[0]; a[2][2]=z[3]-z[0]; 
  a[2][3]=z[2]-z[1]; a[2][4]=z[3]-z[2]; a[2][5]=z[1]-z[3];
  
  
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
      t13[0][0]=0.;          t13[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t13[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t13[1][0]=-t13[0][1];  t13[1][1]=0.;                                  t13[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t13[2][0]=-t13[0][2];  t13[2][1]=-t13[1][2];                          t13[2][2]=0.;
    }
    else if  (i==2){
      t14[0][0]=0.;          t14[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t14[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t14[1][0]=-t14[0][1];  t14[1][1]=0.;                                  t14[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t14[2][0]=-t14[0][2];  t14[2][1]=-t14[1][2];                          t14[2][2]=0.;
    }
    else if  (i==3){
      t23[0][0]=0.;          t23[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t23[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t23[1][0]=-t23[0][1];  t23[1][1]=0.;                                  t23[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t23[2][0]=-t23[0][2];  t23[2][1]=-t23[1][2];                          t23[2][2]=0.;
    }
    else if  (i==4){
      t34[0][0]=0.;          t34[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t34[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t34[1][0]=-t34[0][1];  t34[1][1]=0.;                                  t34[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t34[2][0]=-t34[0][2];  t34[2][1]=-t34[1][2];                          t34[2][2]=0.;
    }
    else if  (i==5){
      t42[0][0]=0.;          t42[0][1]=(an[0]*am[1]-am[0]*an[1])*dl[i]/2.;  t42[0][2]=(an[0]*am[2]-am[0]*an[2])*dl[i]/2.;
      t42[1][0]=-t42[0][1];  t42[1][1]=0.;                                  t42[1][2]=(an[1]*am[2]-am[1]*an[2])*dl[i]/2.;
      t42[2][0]=-t42[0][2];  t42[2][1]=-t42[1][2];                          t42[2][2]=0.;
    }
    
  }
  
}

/**
   function assembles matrix of base functions

   @param n - matrix of base functions
   @param volcoord - volume coordinates
   
   2.10.2008
*/
void lintetrot::bf_matrix (matrix &n,vector &volcoord, vector &x,vector &y,vector &z)
{
  long i,j;
  matrix t12(3,3),t13(3,3),t14(3,3),t23(3,3),t34(3,3),t42(3,3);

  tran_side (t12,t13,t14,t23,t34,t42, x,y,z);

  fillm (0.0,n);

  // Shape function 
  // side 12  f=L1*L2, 13, 14, 23, 34, 42

  n[0][0]=volcoord[0]; n[1][1]=volcoord[0]; n[2][2]=volcoord[0];
  //node 1  - side 12, 13, 14
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      n[j][i+3]= volcoord[0]*volcoord[1]*t12(j,i)+volcoord[0]*volcoord[2]*t13(j,i)+volcoord[0]*volcoord[3]*t14(j,i);
	}
  }
  
  n[0][6]=volcoord[1]; n[1][7]=volcoord[1]; n[2][8]=volcoord[1]; 
   //node 2  - side 12, -, -, 23, -, 42
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      n[j][i+9]= -volcoord[0]*volcoord[1]*t12(j,i)+volcoord[1]*volcoord[2]*t23(j,i)-volcoord[3]*volcoord[1]*t42(j,i);
	}
  }
 
  n[0][12]=volcoord[2]; n[1][13]=volcoord[2]; n[2][14]=volcoord[2]; 
   //node 3  - side -, 13, -, 23, 34, -
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      n[j][i+15]= -volcoord[0]*volcoord[2]*t13(j,i)-volcoord[1]*volcoord[2]*t23(j,i)-volcoord[2]*volcoord[3]*t34(j,i);
	}
  }
  
  n[0][18] =volcoord[3];n[1][19]=volcoord[3];n[2][20]=volcoord[3];
   //node 4  - side -, -, 14, -, 34, 42
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      n[j][i+21]= -volcoord[0]*volcoord[3]*t14(j,i)-volcoord[2]*volcoord[3]*t34(j,i)+volcoord[3]*volcoord[1]*t42(j,i);
	}
  }

}
void lintetrot::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  fillm (0.0,n);
  
  n[0][0]=xi;  n[0][6]=eta;  n[0][12]=zeta;  n[0][18] =1.0-xi-eta-zeta;
  n[1][1]=xi;  n[1][7]=eta;  n[1][13]=zeta;  n[1][19]=1.0-xi-eta-zeta;
  n[2][2]=xi;  n[2][8]=eta;  n[2][14]=zeta;  n[2][20]=1.0-xi-eta-zeta;
}

/**
   function assembles geometric (strain/displacement) %matrix
   vector of strains has following ordering
   eps=(e_xx, e_yy, e_zz, e_yz, e_zx, e_xy)
   
   @param gm - geometric %matrix
   @param x,y,z - vectors of node coordinates
   @param volcoord - array of volume coordinates
   
   PF, 2.10.2008
*/
void lintetrot::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,vector &volcoord)
{
  long i,i1,i2,i3,i4;
  double det;
  vector b(4),c(4),d(4),ds12(4),ds13(4),ds14(4),ds23(4),ds34(4),ds42(4);
  matrix t12(3,3),t13(3,3),t14(3,3),t23(3,3),t34(3,3),t42(3,3);
 
  det = det3d (x.a,y.a,z.a);
  
  volb_3d (b.a,y.a,z.a,det);
  volc_3d (c.a,x.a,z.a,det);
  vold_3d (d.a,x.a,y.a,det);

  // side 12  f=L1*L2, 13, 14, 23, 34, 42
  tran_side (t12,t13,t14,t23,t34,t42, x,y,z);
  
  // ds12[0] is derivation ( d(L1*L2)/dx ), ds12[1] is ( d(L1*L2)/dy )
  ds12[0]=b[0]*volcoord[1]+b[1]*volcoord[0];
  ds12[1]=c[0]*volcoord[1]+c[1]*volcoord[0];
  ds12[2]=d[0]*volcoord[1]+d[1]*volcoord[0];

  ds13[0]=b[0]*volcoord[2]+b[2]*volcoord[0];
  ds13[1]=c[0]*volcoord[2]+c[2]*volcoord[0];
  ds13[2]=d[0]*volcoord[2]+d[2]*volcoord[0];

  ds14[0]=b[0]*volcoord[3]+b[3]*volcoord[0];
  ds14[1]=c[0]*volcoord[3]+c[3]*volcoord[0];
  ds14[2]=d[0]*volcoord[3]+d[3]*volcoord[0];

  ds23[0]=b[1]*volcoord[2]+b[2]*volcoord[1];
  ds23[1]=c[1]*volcoord[2]+c[2]*volcoord[1];
  ds23[2]=d[1]*volcoord[2]+d[2]*volcoord[1];

  ds34[0]=b[2]*volcoord[3]+b[3]*volcoord[2];
  ds34[1]=c[2]*volcoord[3]+c[3]*volcoord[2];
  ds34[2]=d[2]*volcoord[3]+d[3]*volcoord[2];

  ds42[0]=b[1]*volcoord[3]+b[3]*volcoord[1];
  ds42[1]=c[1]*volcoord[3]+c[3]*volcoord[1];
  ds42[2]=d[1]*volcoord[3]+d[3]*volcoord[1];


    // du/dx, dv/dy, dw/dz, dv/dz+dw/dy, du/dz+dw/dx, du/dy+dv/dx
  i1=0;  i2=1;  i3=2;  
  for (i=0;i<nne;i++){
    gm[0][i1]=b[i];  
    gm[1][i2]=c[i]; 
    gm[2][i3]=d[i]; 
    
    gm[3][i2]=d[i];  
    gm[3][i3]=c[i];  
    gm[4][i1]=d[i];  
    gm[4][i3]=b[i];  
    gm[5][i1]=c[i];  
    gm[5][i2]=b[i]; 
    i1+=6;  i2+=6;  i3+=6;
  }

  i1=3; i2=9; i3=15; i4=21;
  for (i=0;i<3;i++){
    // node 1  du/dx, dv/dy, dw/dz, 
    gm[0][i1]=-t12[0][i]*ds12[0]-t13[0][i]*ds13[0]-t14[0][i]*ds14[0];  
    gm[1][i1]=-t12[1][i]*ds12[1]-t13[1][i]*ds13[1]-t14[1][i]*ds14[1];  
    gm[2][i1]=-t12[2][i]*ds12[2]-t13[2][i]*ds13[2]-t14[2][i]*ds14[2]; 

	//node 1   dv/dz+dw/dy,  du/dz+dw/dx,  du/dy+dv/dx
    gm[3][i1]=-t12[1][i]*ds12[2]-t13[1][i]*ds13[2]-t14[1][i]*ds14[2]
              -t12[2][i]*ds12[1]-t13[2][i]*ds13[1]-t14[2][i]*ds14[1];  
    gm[4][i1]=-t12[2][i]*ds12[0]-t13[2][i]*ds13[0]-t14[2][i]*ds14[0]
              -t12[0][i]*ds12[2]-t13[0][i]*ds13[2]-t14[0][i]*ds14[2];  
    gm[5][i1]=-t12[0][i]*ds12[1]-t13[0][i]*ds13[1]-t14[0][i]*ds14[1]
              -t12[1][i]*ds12[0]-t13[1][i]*ds13[0]-t14[1][i]*ds14[0];  
    
    // node 2
    gm[0][i2]= t12[0][i]*ds12[0]-t23[0][i]*ds23[0]+t42[0][i]*ds42[0]; 
    gm[1][i2]= t12[1][i]*ds12[1]-t23[1][i]*ds23[1]+t42[1][i]*ds42[1];  
    gm[2][i2]= t12[2][i]*ds12[2]-t23[2][i]*ds23[2]+t42[2][i]*ds42[2]; 
    
    
    gm[3][i2]= t12[1][i]*ds12[2]-t23[1][i]*ds23[2]+t42[1][i]*ds42[2]
              +t12[2][i]*ds12[1]-t23[2][i]*ds23[1]+t42[2][i]*ds42[1];  
    gm[4][i2]= t12[2][i]*ds12[0]-t23[2][i]*ds23[0]+t42[2][i]*ds42[0]
              +t12[0][i]*ds12[2]-t23[0][i]*ds23[2]+t42[0][i]*ds42[2];  
    gm[5][i2]= t12[0][i]*ds12[1]-t23[0][i]*ds23[1]+t42[0][i]*ds42[1]
              +t12[1][i]*ds12[0]-t23[1][i]*ds23[0]+t42[1][i]*ds42[0];  
    
    
    // node 3
    gm[0][i3]= t23[0][i]*ds23[0]+t13[0][i]*ds13[0]-t34[0][i]*ds34[0]; 
    gm[1][i3]= t23[1][i]*ds23[1]+t13[1][i]*ds13[1]-t34[1][i]*ds34[1]; 
    gm[2][i3]= t23[2][i]*ds23[2]+t13[2][i]*ds13[2]-t34[2][i]*ds34[2]; 
    
    gm[3][i3]= t23[1][i]*ds23[2]+t13[1][i]*ds13[2]-t34[1][i]*ds34[2]
              +t23[2][i]*ds23[1]+t13[2][i]*ds13[1]-t34[2][i]*ds34[1];  
    gm[4][i3]= t23[2][i]*ds23[0]+t13[2][i]*ds13[0]-t34[2][i]*ds34[0]
              +t23[0][i]*ds23[2]+t13[0][i]*ds13[2]-t34[0][i]*ds34[2];  
    gm[5][i3]= t23[0][i]*ds23[1]+t13[0][i]*ds13[1]-t34[0][i]*ds34[1]
              +t23[1][i]*ds23[0]+t13[1][i]*ds13[0]-t34[1][i]*ds34[0];  
    
    
    // node 4
    gm[0][i4]= t14[0][i]*ds14[0]-t42[0][i]*ds42[0]+t34[0][i]*ds34[0];  
    gm[1][i4]= t14[1][i]*ds14[1]-t42[1][i]*ds42[1]+t34[1][i]*ds34[1];  
    gm[2][i4]= t14[2][i]*ds14[2]-t42[2][i]*ds42[2]+t34[2][i]*ds34[2];  
    
    gm[3][i4]= t14[2][i]*ds14[0]-t42[2][i]*ds42[0]+t34[2][i]*ds34[0]
              +t14[0][i]*ds14[2]-t42[0][i]*ds42[2]+t34[0][i]*ds34[2];  
    gm[4][i4]= t14[0][i]*ds14[1]-t42[0][i]*ds42[1]+t34[0][i]*ds34[1]
              +t14[1][i]*ds14[0]-t42[1][i]*ds42[0]+t34[1][i]*ds34[0];  
    gm[5][i4]= t14[1][i]*ds14[2]-t42[1][i]*ds42[2]+t34[1][i]*ds34[2]
              +t14[2][i]*ds14[1]-t42[2][i]*ds42[1]+t34[2][i]*ds34[1];  
    i1+=1; i2+=1; i3+=1; i4+=1;
 }
  
}

/**
   function assembles shear part of geometric (strain/displacement) %matrix
   
   @param gm - geometric %matrix
   @param x,y,z - vectors of nodal coordinates
   @param volcoord - array of volume coordinates
   
   PF, 2.10.2008
*/
void lintetrot::geom_matrix_shear (matrix &gm,vector &x,vector &y,vector &z,vector &volcoord)
{
  long i,i1,i2,i3,i4;
  double det;
  vector b(4),c(4),d(4),ds12(4),ds13(4),ds14(4),ds23(4),ds34(4),ds42(4);
  matrix t12(3,3),t13(3,3),t14(3,3),t23(3,3),t34(3,3),t42(3,3);

  det = det3d (x.a,y.a,z.a);
  
  volb_3d (b.a,y.a,z.a,det);
  volc_3d (c.a,x.a,z.a,det);
  vold_3d (d.a,x.a,y.a,det);

  // side 12  f=L1*L2, 13, 14, 23, 34, 42
  tran_side (t12,t13,t14,t23,t34,t42, x,y,z);

  // ds12[0] is derivation ( d(L1*L2)/dx ), ds12[1] is ( d(L1*L2)/dy )
  ds12[0]=b[0]*volcoord[1]+b[1]*volcoord[0];
  ds12[1]=c[0]*volcoord[1]+c[1]*volcoord[0];
  ds12[2]=d[0]*volcoord[1]+d[1]*volcoord[0];

  ds13[0]=b[0]*volcoord[2]+b[2]*volcoord[0];
  ds13[1]=c[0]*volcoord[2]+c[2]*volcoord[0];
  ds13[2]=d[0]*volcoord[2]+d[2]*volcoord[0];

  ds14[0]=b[0]*volcoord[3]+b[3]*volcoord[0];
  ds14[1]=c[0]*volcoord[3]+c[3]*volcoord[0];
  ds14[2]=d[0]*volcoord[3]+d[3]*volcoord[0];

  ds23[0]=b[1]*volcoord[2]+b[2]*volcoord[1];
  ds23[1]=c[1]*volcoord[2]+c[2]*volcoord[1];
  ds23[2]=d[1]*volcoord[2]+d[2]*volcoord[1];

  ds34[0]=b[2]*volcoord[3]+b[3]*volcoord[2];
  ds34[1]=c[2]*volcoord[3]+c[3]*volcoord[2];
  ds34[2]=d[2]*volcoord[3]+d[3]*volcoord[2];

  ds42[0]=b[1]*volcoord[3]+b[3]*volcoord[1];
  ds42[1]=c[1]*volcoord[3]+c[3]*volcoord[1];
  ds42[2]=d[1]*volcoord[3]+d[3]*volcoord[1];


  i1=3; i2=9; i3=15; i4=21;
  for (i=0;i<3;i++){
	//node 1   -dv/dz+dw/dy,  du/dz-dw/dx,  -du/dy+dv/dx
    gm[0][i1]=(t12[1][i]*ds12[2]+t13[1][i]*ds13[2]+t14[1][i]*ds14[2]
	       -t12[2][i]*ds12[1]-t13[2][i]*ds13[1]-t14[2][i]*ds14[1])/2.;  
    gm[1][i1]=(t12[2][i]*ds12[0]+t13[2][i]*ds13[0]+t14[2][i]*ds14[0]
	       -t12[0][i]*ds12[2]-t13[0][i]*ds13[2]-t14[0][i]*ds14[2])/2.;  
    gm[2][i1]=(t12[0][i]*ds12[1]+t13[0][i]*ds13[1]+t14[0][i]*ds14[1]
	       -t12[1][i]*ds12[0]-t13[1][i]*ds13[0]-t14[1][i]*ds14[0])/2.;  
    // node 2
    gm[0][i2]=(-t12[1][i]*ds12[2]+t23[1][i]*ds23[2]-t42[1][i]*ds42[2]
	       +t12[2][i]*ds12[1]-t23[2][i]*ds23[1]+t42[2][i]*ds42[1])/2.;  
    gm[1][i2]=(-t12[2][i]*ds12[0]+t23[2][i]*ds23[0]-t42[2][i]*ds42[0]
	       +t12[0][i]*ds12[2]-t23[0][i]*ds23[2]+t42[0][i]*ds42[2])/2.;  
    gm[2][i2]=(-t12[0][i]*ds12[1]+t23[0][i]*ds23[1]-t42[0][i]*ds42[1]
	       +t12[1][i]*ds12[0]-t23[1][i]*ds23[0]+t42[1][i]*ds42[0])/2.;  
    // node 3
    gm[0][i3]=(-t23[1][i]*ds23[2]-t13[1][i]*ds13[2]+t34[1][i]*ds34[2]
	       +t23[2][i]*ds23[1]+t13[2][i]*ds13[1]-t34[2][i]*ds34[1])/2.;  
    gm[1][i3]=(-t23[2][i]*ds23[0]-t13[2][i]*ds13[0]+t34[2][i]*ds34[0]
	       +t23[0][i]*ds23[2]+t13[0][i]*ds13[2]-t34[0][i]*ds34[2])/2.;  
    gm[2][i3]=(-t23[0][i]*ds23[1]-t13[0][i]*ds13[1]+t34[0][i]*ds34[1]
	       +t23[1][i]*ds23[0]+t13[1][i]*ds13[0]-t34[1][i]*ds34[0])/2.;  
    // node 4
    gm[0][i4]=(-t14[1][i]*ds14[2]+t42[1][i]*ds42[1]-t34[1][i]*ds34[2]
	       +t14[2][i]*ds14[1]-t42[2][i]*ds42[0]+t34[2][i]*ds34[1])/2.;  
    gm[1][i4]=(-t14[2][i]*ds14[0]+t42[2][i]*ds42[2]-t34[2][i]*ds34[0]
	       +t14[0][i]*ds14[2]-t42[0][i]*ds42[1]+t34[0][i]*ds34[2])/2.;  
    gm[2][i4]=(-t14[0][i]*ds14[1]+t42[0][i]*ds42[0]-t34[0][i]*ds34[1]
	       +t14[1][i]*ds14[0]-t42[1][i]*ds42[2]+t34[1][i]*ds34[0])/2.;  
    i1+=1; i2+=1; i3+=1; i4+=1;
  }
  
  i1=0;  i2=1;  i3=2; 
  for (i=0;i<nne;i++){
  //  gm[0][i2]  =-d[i]/2.;  
  //  gm[0][i3]  = c[i]/2.;
    gm[0][i1+3]-=volcoord[i];
    
  //  gm[1][i1]  = d[i]/2.;  
  //  gm[1][i3]  =-b[i]/2.;
    gm[1][i1+4]-= volcoord[i];
    
  //  gm[2][i1]  =-c[i]/2.; 
  //  gm[2][i2]  = b[i]/2.;  
    gm[2][i1+5]-= volcoord[i];
    i1+=6;  i2+=6;  i3+=6;
  }
 
  
}
/**
   function assembles stiffness matrices for bending and shear

   function extracts components of the bending
   part of stiffness %matrix of the material to the matrix db

   10.8.2008
*/
void lintetrot::dmatblock (matrix &dd, matrix &d)
{
  dd[0][0] += d[3][3]/intordsm[0][0];  dd[1][1] += d[4][4]/intordsm[0][0];  dd[2][2] += d[5][5]/intordsm[0][0];

}
/**
   function assembles transformation %matrix from nodal
   to global coordinate system
   
   @param nodes - nodes of element
   @param tmat - transformation %matrix
   
   JK, 2.10.2008
*/
void lintetrot::transf_matrix (ivector &nodes,matrix &tmat)
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

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   PF, 2.10.2008
*/
void lintetrot::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long ipp,i,  ipshear;
  double det,jac;
  vector x(nne),y(nne),z(nne),gp1,gp2,gp3,w,volcoord(4);
  matrix gm,d(tncomp,tncomp), dd(3,3);
  
  Mt->give_node_coord3d (x,y,z,eid);
  fillm (0.0,sm);
  fillm (0.0,dd);
  det = det3d (x.a,y.a,z.a);
  det = fabs(det);
    
    reallocv (intordsm[0][0],w);
    reallocv (intordsm[0][0],gp1);
    reallocv (intordsm[0][0],gp2);
    reallocv (intordsm[0][0],gp3);
    reallocm (ncomp[0],ndofe,gm);

    ipp=Mt->elements[eid].ipp[ri][ci];
    gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
            
    for (i=0;i<intordsm[0][0];i++){
	    volcoord[0]=gp1[i];
	    volcoord[1]=gp2[i];
	    volcoord[2]=gp3[i];
	    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];
	//  stiffness matrix of the material
	    Mm->matstiff (d,ipp);  ipp++;
	    jac=w[i]*det;
	    if (jac<0.0)
              print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	//  geometric matrix
	    geom_matrix (gm,x,y,z,volcoord);
	//  contribution to the stiffness matrix of the element
	    bdbjac (sm,gm,d,gm,jac);

	// average for shear	
      dmatblock ( dd, d);
    }

// Shear
    ipshear=1;
    reallocv (ipshear,w);
    reallocv (ipshear,gp1);
    reallocv (ipshear,gp2);
    reallocv (ipshear,gp3);
    reallocm (3,ndofe,gm);
    gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,ipshear);
   for (i=0;i<ipshear;i++){
	    volcoord[0]=gp1[i];
	    volcoord[1]=gp2[i];
	    volcoord[2]=gp3[i];
	    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];
	    jac=w[i]*det;
	    geom_matrix_shear (gm,x,y,z,volcoord);
	    bdbjac (sm,gm,dd,gm,jac);
      }
      
/*
  fprintf (Out,"\n\n SMloc");
    for (i=0;i<4;i++){
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
   function assembles resulting stiffness %matrix of the element

   @param eid - element id
   @param sm - stiffness %matrix

   PF, 2.10.2008
*/
void lintetrot::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);

  //  transformation of stiffness matrix from global to nodal coordinate system
  ivector nodes (nne);
  long transf;
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
    
}

/**
   function computes mass %matrix
   
   @param eid - number of element
   @param mm - mass %matrix
   
   26.8.2001
*/
void lintetrot::mass_matrix (long eid,matrix &mm)
{
  long i;
  double det,jac,rho;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp1(intordmm),gp2(intordmm),gp3(intordmm),volcoord(4),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);
  Mt->give_node_coord3d (x,y,z,eid);

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordmm);

  det = det3d (x.a,y.a,z.a);
  det = fabs(det);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];
    
    bf_matrix (n,volcoord, x, y, z);
    rho = approx (volcoord,dens);
    
    jac=det*w[i]*rho;
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function computes load %matrix
   
   @param eid - number of element
   @param lm - load %matrix

   26.8.2001
*/
void lintetrot::load_matrix (long eid,matrix &lm)
{
  long i;
  double det,jac;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp1(intordmm),gp2(intordmm),gp3(intordmm),volcoord(4);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordmm);

  det = det3d (x.a,y.a,z.a);
  det = fabs(det);

  fillm (0.0,lm);

  for (i=0;i<intordmm;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];

    bf_matrix (n,volcoord, x, y, z);

    jac=det*w[i];

    nnj (lm.a,n.a,jac,n.m,n.n);
  }
}





void lintetrot::res_ip_strains (long lcid,long eid)
{
  vector x(nne),y(nne),z(nne),r(ndofe),eps(tncomp),aux;
  ivector nodes(nne);
  matrix gm(tncomp,ndofe),tmat;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector from nodal to global coordinate system
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  ip_strains (lcid,eid,0,0,x,y,z,r);
  
}

/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y,z - vectors of node coordinates
   @param r - %vector of node displacements
   
   JK, 2.10.2008
*/
void lintetrot::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,ipp;
  vector gp1,gp2,gp3,w,eps(tncomp),volcoord(4);
  matrix gm(tncomp,ndofe);
   
    ipp=Mt->elements[eid].ipp[ri][ci];
    reallocv (intordsm[0][0],w);
    reallocv (intordsm[0][0],gp1);
    reallocv (intordsm[0][0],gp2);
    reallocv (intordsm[0][0],gp3);
    
    gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
    
    for (i=0;i<intordsm[0][0];i++){
      volcoord[0]=gp1[i];
      volcoord[1]=gp2[i];
      volcoord[2]=gp3[i];
      volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];
      
      geom_matrix (gm,x,y,z,volcoord);
      mxv (gm,r,eps);
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
}

/**
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void lintetrot::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  double vol;
  ivector nod(nne);
  vector eps(tncomp);

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];

  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  Mm->givestrain (lcid,ipp,eps);
  
  for (i=0;i<nne;i++){
    //  storage of strains to the node
    j=nod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2){
      vol=volumeip (eid,ri,ci);
      cmulv (vol,eps,eps);
      Mt->nodes[j].storestrain (lcid,0,vol,eps);
    }
  }
}



void lintetrot::strains (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{

}


/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 2.10.2008
*/
void lintetrot::res_ip_stresses (long lcid,long eid)
{
  compute_nlstress (lcid,eid,0,0);
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void lintetrot::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  double vol;
  vector nxi(nne),neta(nne),nzeta(nne),eps(tncomp),epst,epstt,sig(tncomp),natcoord(3);
  ivector nod(nne);
  
  Mt->give_elemnodes (eid,nod);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  stresses at the closest integration point
  Mm->givestress (lcid,ipp,sig);
  
  for (i=0;i<nne;i++){
    
    //  storage of stresses to the node
    j=nod[i];
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2){
      vol=volumeip (eid,ri,ci);
      cmulv (vol,sig,sig);
      Mt->nodes[j].storestress (lcid,0,vol,sig);
    }
  }
  
  
}



void lintetrot::stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
}








/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void lintetrot::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes internal forces for nonlocal models

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void lintetrot::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
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
void lintetrot::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes nodal forces caused by temperature changes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y,z - nodal coordinates
   
   7.2008, TKo
*/
void lintetrot::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
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
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void lintetrot::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y,z);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
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
void lintetrot::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   function computes resulting increment of internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void lintetrot::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   function computes resulting contributions from eigenstrains to the right hand side
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - %vector of internal forces

   TKo, 7.2008
*/
void lintetrot::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void lintetrot::compute_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp, i;
  
    ipp=Mt->elements[eid].ipp[ri][ci];
    for (i=0;i<intordsm[0][0];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
	      Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
    }
}



/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo  7.2008
*/
void lintetrot::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstressesincr (ipp);
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void lintetrot::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void lintetrot::compute_nonloc_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->compnonloc_nlstresses (ipp);
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void lintetrot::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ipp=Mt->elements[eid].ipp[ri][ci];
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
   
   TKo 7.2008
*/
void lintetrot::elem_integration(integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z)
{
  long ipp;
  double det,jac;
  vector ipv(tncomp),contr(ndofe);
  matrix gm(tncomp,ndofe);

  ipp=Mt->elements[eid].ipp[ri][ci];
  /*
  xi=0.25;
  eta=0.25;
  zeta=0.25;*/
  fillv (0.0,nv);
  det = det3d (x.a,y.a,z.a);  
  det = fabs(det);

  //  function assembles required quantity at integration point
  Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);

  //  strain-displacement (geometric) matrix
//  geom_matrix (gm,x,y,z);
  //  contribution to the nodal values
  mtxv (gm,ipv,contr);
  jac=det/6.0;
  cmulv (jac,contr);
  //  summation
  addv(contr,nv,nv);
}



/**
   The function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param coord - %vector with global coordinates of integration point (ouput)
   
   @return The function returns global coordinates in the argument coord.

   19.1.2002
*/
void lintetrot::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  vector x(nne),y(nne),z(nne),volcoord(4),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]),gp3(intordsm[ri][ci]);

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];

    if (ii==ipp){
      coord[0]=approx (volcoord,x);
      coord[1]=approx (volcoord,y);
      coord[2]=approx (volcoord,y);
    }
    ii++;
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
void lintetrot::ipncoord (long eid,long ipp,vector &ncoord)
{
  long i, ii, ri, ci;
  vector w, gp1 ,gp2, gp3;

  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      if (intordsm[ri][ci] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp3));
      gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];

      for (i=0;i<intordsm[ri][ci];i++){
        if (ii==ipp){
          ncoord[0]=gp1[i];
          ncoord[1]=gp2[i];
          ncoord[2]=gp3[i];
        }
        ii++;
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
void lintetrot::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, zeta, ipval;
  vector w, gp1, gp2, gp3, anv(nne);
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
        reallocv (intordsm[ii][jj],gp1);
        reallocv (intordsm[ii][jj],gp2);
        reallocv (intordsm[ii][jj],gp3);
        reallocv (intordsm[ii][jj],w);
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][ii]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp1[k];
          eta=gp2[k];
          zeta=gp3[k];
          //  value in integration point
          ipval = approx_nat (xi,eta,zeta,anv);
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
void lintetrot::ipvolume (long eid,long ri,long ci)
{
  long ipp;
  double jac,det;
  vector x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  det = det3d (x.a,y.a,z.a);
  jac=fabs(det)/6.0;
  
  Mm->storeipvol (ipp,jac);
}
/**
   function computes volume appropriate to integration point
   
   2.3.2004, JK
*/
double lintetrot::volumeip (long eid,long ri,long ci)
{
  double jac,det;
  vector x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  det = det3d (x.a,y.a,z.a);
  jac=fabs(det)/6.0;
  
  return jac;
}


void lintetrot::nod_eqother_ip (long lcid,long eid)
{
  long i,j,ipp,ncompo;
  double vol;
  ivector nod(nne);
  vector eqother;
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  ncompo = Mm->givencompeqother (ipp,0);
  reallocv (ncompo,eqother);
  Mm->giveeqother (ipp,0,ncompo,eqother.a);
  
  for (i=0;i<nne;i++){
    //  storage of other to the node
    j=nod[i];
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother (lcid,0,ncompo,eqother);
    if (Mp->otheraver==2){
      vol=volumeip (eid,0,0);
      cmulv (vol,eqother,eqother);
      Mt->nodes[j].storeother (0,ncompo,vol,eqother);
    }
  }
}

/**
   function computes nodal forces caused by presure on surface
   
   @param lcid - number of load case
   @param eid - element id
   @param is - identification of surfaces
   @param nv - nodal values
   @param nf - nodal forces
   
   5.4.2005, JK
   pracuje pouze v globalnim souradnem systemu
*/
void lintetrot::node_forces_surf_old (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),v(ndofe),av(ndofe),gp1(intordb),gp2(intordb),w(intordb);
  matrix am(ndofe,ndofe),n(napfun,ndofe);
  
  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  if (is[0]==1){
    for (i=0;i<intordb;i++){
      xi=0.0; eta=gp1[i];  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,eta,zeta,0);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[6] = nv[0*3+0];
    av[7] = nv[0*3+1];
    av[8] = nv[0*3+2];
    
    av[3] = nv[1*3+0];
    av[4] = nv[1*3+1];
    av[5] = nv[1*3+2];
    
    av[9]  = nv[2*3+0];
    av[10] = nv[2*3+1];
    av[11] = nv[2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[1]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=0.0;  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,xi,zeta,1);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[0] = nv[9+0*3+0];
    av[1] = nv[9+0*3+1];
    av[2] = nv[9+0*3+2];
    
    av[6] = nv[9+1*3+0];
    av[7] = nv[9+1*3+1];
    av[8] = nv[9+1*3+2];
    
    av[9]  = nv[9+2*3+0];
    av[10] = nv[9+2*3+1];
    av[11] = nv[9+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[2]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=0.0;
      
      jac2d_3d (jac,x,y,z,xi,eta,2);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[3] = nv[18+0*3+0];
    av[4] = nv[18+0*3+1];
    av[5] = nv[18+0*3+2];
    
    av[0] = nv[18+1*3+0];
    av[1] = nv[18+1*3+1];
    av[2] = nv[18+1*3+2];
    
    av[9]  = nv[18+2*3+0];
    av[10] = nv[18+2*3+1];
    av[11] = nv[18+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[3]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=1.0-xi-eta;
      
      jac2d_3d (jac,x,y,z,xi,eta,3);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[3] = nv[27+0*3+0];
    av[4] = nv[27+0*3+1];
    av[5] = nv[27+0*3+2];
    
    av[6] = nv[27+1*3+0];
    av[7] = nv[27+1*3+1];
    av[8] = nv[27+1*3+2];
    
    av[0] = nv[27+2*3+0];
    av[1] = nv[27+2*3+1];
    av[2] = nv[27+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
}


/**
   function computes nodal forces caused by presure on surface
   
   @param lcid - number of load case
   @param eid - element id
   @param is - identification of surfaces
   @param nv - nodal values
   @param nf - nodal forces
   
   12.12.2005, JK
*/
void lintetrot::node_forces_surf (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i;
  double xi,eta,zeta,jac;
  double *tnv;
  vector x(nne),y(nne),z(nne),v(ndofe),av(ndofe),gp1(intordb),gp2(intordb),w(intordb);
  matrix am(ndofe,ndofe),n(napfun,ndofe);
  
  tnv = new double [9];
  
  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  if (is[0]>0 ){
    for (i=0;i<intordb;i++){
      xi=0.0; eta=gp1[i];  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,eta,zeta,0);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[0]==1){
      av[6] = nv[0*3+0];
      av[7] = nv[0*3+1];
      av[8] = nv[0*3+2];
      
      av[3] = nv[1*3+0];
      av[4] = nv[1*3+1];
      av[5] = nv[1*3+2];
      
      av[9] = nv[2*3+0];
      av[10] = nv[2*3+1];
      av[11] = nv[2*3+2];
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
      
      locglob_nodeval (0,av,tnv,x,y,z);
      
      av[6] = tnv[0];
      av[7] = tnv[1];
      av[8] = tnv[2];
      
      av[3] = tnv[3];
      av[4] = tnv[4];
      av[5] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }


    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[1]>0){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=0.0;  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,xi,zeta,1);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[1]==1){
      av[0] = nv[9+0*3+0];
      av[1] = nv[9+0*3+1];
      av[2] = nv[9+0*3+2];
      
      av[6] = nv[9+1*3+0];
      av[7] = nv[9+1*3+1];
      av[8] = nv[9+1*3+2];
      
      av[9] = nv[9+2*3+0];
      av[10] = nv[9+2*3+1];
      av[11] = nv[9+2*3+2];
    }
    
    if (is[1]==2){
      
      av[0] = nv[9+0*3+0];
      av[1] = nv[9+0*3+1];
      av[2] = nv[9+0*3+2];
      
      av[3] = nv[9+1*3+0];
      av[4] = nv[9+1*3+1];
      av[5] = nv[9+1*3+2];
      
      av[6] = nv[9+2*3+0];
      av[7] = nv[9+2*3+1];
      av[8] = nv[9+2*3+2];
      
      locglob_nodeval (1,av,tnv,x,y,z);
      
      av[0] = tnv[0];
      av[1] = tnv[1];
      av[2] = tnv[2];
      
      av[6] = tnv[3];
      av[7] = tnv[4];
      av[8] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[2]>0){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=0.0;
      
      jac2d_3d (jac,x,y,z,xi,eta,2);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[2]==1){
      av[3] = nv[18+0*3+0];
      av[4] = nv[18+0*3+1];
      av[5] = nv[18+0*3+2];
      
      av[0] = nv[18+1*3+0];
      av[1] = nv[18+1*3+1];
      av[2] = nv[18+1*3+2];
      
      av[9] = nv[18+2*3+0];
      av[10] = nv[18+2*3+1];
      av[11] = nv[18+2*3+2];

    }

    if (is[2]==2){
      
      av[0] = nv[18+0*3+0];
      av[1] = nv[18+0*3+1];
      av[2] = nv[18+0*3+2];
      
      av[3] = nv[18+1*3+0];
      av[4] = nv[18+1*3+1];
      av[5] = nv[18+1*3+2];
      
      av[6] = nv[18+2*3+0];
      av[7] = nv[18+2*3+1];
      av[8] = nv[18+2*3+2];
      
      locglob_nodeval (2,av,tnv,x,y,z);
      
      av[3] = tnv[0];
      av[4] = tnv[1];
      av[5] = tnv[2];
      
      av[0] = tnv[3];
      av[1] = tnv[4];
      av[2] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[3]>0){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=1.0-xi-eta;
      
      jac2d_3d (jac,x,y,z,xi,eta,3);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[3]==1){
      av[3] = nv[27+0*3+0];
      av[4] = nv[27+0*3+1];
      av[5] = nv[27+0*3+2];
      
      av[6] = nv[27+1*3+0];
      av[7] = nv[27+1*3+1];
      av[8] = nv[27+1*3+2];
      
      av[0] = nv[27+2*3+0];
      av[1] = nv[27+2*3+1];
      av[2] = nv[27+2*3+2];

    }

    if (is[3]==2){
      
      av[0] = nv[27+0*3+0];
      av[1] = nv[27+0*3+1];
      av[2] = nv[27+0*3+2];
      
      av[3] = nv[27+1*3+0];
      av[4] = nv[27+1*3+1];
      av[5] = nv[27+1*3+2];
      
      av[6] = nv[27+2*3+0];
      av[7] = nv[27+2*3+1];
      av[8] = nv[27+2*3+2];
      
      locglob_nodeval (3,av,tnv,x,y,z);
      
      av[3] = tnv[0];
      av[4] = tnv[1];
      av[5] = tnv[2];
      
      av[6] = tnv[3];
      av[7] = tnv[4];
      av[8] = tnv[5];
      
      av[0] = tnv[6];
      av[1] = tnv[7];
      av[2] = tnv[8];
    }
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  
  delete [] tnv;
}


void lintetrot::locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z)
{
  double norm;
  vector g1(3),g2(3),g3(3),lv(3),gv(3);
  matrix t(3,3);
  
  if (is==0){
    g1[0]=x[2]-x[3];
    g1[1]=y[2]-y[3];
    g1[2]=z[2]-z[3];
    
    g2[0]=x[1]-x[3];
    g2[1]=y[1]-y[3];
    g2[2]=z[1]-z[3];
  }
  if (is==1){
    g1[0]=x[0]-x[3];
    g1[1]=y[0]-y[3];
    g1[2]=z[0]-z[3];
    
    g2[0]=x[2]-x[3];
    g2[1]=y[2]-y[3];
    g2[2]=z[2]-z[3];
  }
  if (is==2){
    g1[0]=x[1]-x[3];
    g1[1]=y[1]-y[3];
    g1[2]=z[1]-z[3];
    
    g2[0]=x[0]-x[3];
    g2[1]=y[0]-y[3];
    g2[2]=z[0]-z[3];
  }
  if (is==3){
    g1[0]=x[1]-x[0];
    g1[1]=y[1]-y[0];
    g1[2]=z[1]-z[0];
    
    g2[0]=x[2]-x[0];
    g2[1]=y[2]-y[0];
    g2[2]=z[2]-z[0];
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
  
}

/**
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   JK, 22.4.2005
*/
void lintetrot::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  vector volcoord(4);
  
  volcoord[0]=0.25;
  volcoord[1]=0.25;
  volcoord[2]=0.25;
  volcoord[3]=0.25;
  
  ipval[0]=approx (volcoord,nodval);
}


/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 22.4.2005
*/
void lintetrot::res_eigstrain_forces (long eid,vector &nfor)
{
  eigstrain_forces (eid,nfor);
}

/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 22.4.2005
*/
void lintetrot::eigstrain_forces (long eid,vector &nfor)
{
  long i,ipp;
  double jac,det;
  vector x(nne),y(nne),z(nne),eigstr(tncomp),sig(tncomp),contr(ndofe);
  matrix d(tncomp,tncomp),gm(tncomp,ndofe);
  
  Mt->give_node_coord3d (x,y,z,eid);
  det = det3d (x.a,y.a,z.a);
  det = fabs(det);
  
  fillv (0.0,nfor);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  Mm->giveeigstrain (ipp,cncomp[0],ncomp[0],eigstr);
  
  //  matrix of stiffness of the material
  Mm->matstiff (d,ipp);
  
  mxv (d,eigstr,sig);
  
//  geom_matrix (gm,x,y,z);
  
  mtxv (gm,sig,contr);
  
  jac = det/6.0;

  cmulv (jac,contr);
  
  for (i=0;i<contr.n;i++){
    nfor[i]+=contr[i];
  }
  
}


