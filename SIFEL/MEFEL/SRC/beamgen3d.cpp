#include "beamgen3d.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include <math.h>


beamgen3d::beamgen3d (void)
{
  long i,j;
  
  nne=2;  ndofe=14;  tncomp=8;  napfun=8;
  nb=1;  intordmm=4;  intordism=2;  ssst=spacebeam;

  ncomp = new long [nb];
  ncomp[0]=1;
  
  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  //  there are two integration points
  //  each integration point stores nodal forces and moments at one node
  nip[0][0]=2;
  intordsm[0][0]=1;
  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
}

beamgen3d::~beamgen3d (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
  }
  delete [] nip;
}


/**
   function assembles transformation %matrix from global to nodal system x_n = T x_g
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   PF, 20.12.2002
*/
void beamgen3d::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,i6,n,m;

  fillm (0.0,tmat);
  
  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      i6=i*7;
      tmat[i6][i6]   = Mt->nodes[nodes[i]].e1[0];  tmat[i6][i6+1]   = Mt->nodes[nodes[i]].e2[0];  tmat[i6][i6+2]   = Mt->nodes[nodes[i]].e3[0];
      tmat[i6+1][i6] = Mt->nodes[nodes[i]].e1[1];  tmat[i6+1][i6+1] = Mt->nodes[nodes[i]].e2[1];  tmat[i6+1][i6+2] = Mt->nodes[nodes[i]].e3[1];
      tmat[i6+2][i6] = Mt->nodes[nodes[i]].e1[2];  tmat[i6+2][i6+1] = Mt->nodes[nodes[i]].e2[2];  tmat[i6+2][i6+2] = Mt->nodes[nodes[i]].e3[2];
      i6=i*7+3;
      tmat[i6][i6]   = Mt->nodes[nodes[i]].e1[0];  tmat[i6][i6+1]   = Mt->nodes[nodes[i]].e2[0];  tmat[i6][i6+2]   = Mt->nodes[nodes[i]].e3[0];
      tmat[i6+1][i6] = Mt->nodes[nodes[i]].e1[1];  tmat[i6+1][i6+1] = Mt->nodes[nodes[i]].e2[1];  tmat[i6+1][i6+2] = Mt->nodes[nodes[i]].e3[1];
      tmat[i6+2][i6] = Mt->nodes[nodes[i]].e1[2];  tmat[i6+2][i6+1] = Mt->nodes[nodes[i]].e2[2];  tmat[i6+2][i6+2] = Mt->nodes[nodes[i]].e3[2];
      tmat[i6+3][i6] = 1.;      
    }
  }
}


/**
   function assembles transformation %matrix from local to global system x_g = T x_l
   
   columns of the transformation %matrix are created by coordinates
   of local base vectors expressed in global system
   
   @param tmat - transformation %matrix
   @param dl - lenght of the beam
   @param vec - direction %vector of local axis z (in global coordinates)
   @param x,y,z - vectors of nodal coordinates in global system
   @param eid - element id
   
   JK, 3.2.2002, revised 1.9.2006
*/
void beamgen3d::beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid)
{
  double c;

  fillm (0.0,tmat);
  
  //  local vector x_l
  tmat[0][0]=x[1]-x[0];
  tmat[1][0]=y[1]-y[0];
  tmat[2][0]=z[1]-z[0];
  
  //  length of the beam, it is equal to the norm of x_l
  dl=sqrt((tmat[0][0]*tmat[0][0]+tmat[1][0]*tmat[1][0]+tmat[2][0]*tmat[2][0]));
  
  if (dl<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beamgen3d element",eid);
    fprintf (stderr,"\n in function beamgen3d::beam_transf_matrix (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  //  normed local vector x_l
  tmat[0][0]=tmat[0][0]/dl;
  tmat[1][0]=tmat[1][0]/dl;
  tmat[2][0]=tmat[2][0]/dl;
  
  //  local vector y_l
  tmat[0][1]=vec[1]*tmat[2][0]-vec[2]*tmat[1][0];
  tmat[1][1]=vec[2]*tmat[0][0]-vec[0]*tmat[2][0];
  tmat[2][1]=vec[0]*tmat[1][0]-vec[1]*tmat[0][0];
  
  //  norm of the local vector y_l
  c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));

  if (c<Mp->zero){
    vec[0]=0.0;
    vec[1]=1.0;
    vec[2]=0.0;
    
    //  local vector z_l (defined by vector product)
    tmat[0][2]=tmat[1][0]*vec[2]-tmat[2][0]*vec[1];
    tmat[1][2]=tmat[2][0]*vec[0]-tmat[0][0]*vec[2];
    tmat[2][2]=tmat[0][0]*vec[1]-tmat[1][0]*vec[0];

    //  norm of the local vector z_l
    c=sqrt((tmat[0][2]*tmat[0][2]+tmat[1][2]*tmat[1][2]+tmat[2][2]*tmat[2][2]));

    //  normed local vector z_l
    tmat[0][2]=tmat[0][2]/c;
    tmat[1][2]=tmat[1][2]/c;
    tmat[2][2]=tmat[2][2]/c;
    
    if (c<Mp->zero){
      fprintf (stderr,"\n\n zero z base vector of the %ld beamgen3d element",eid);
      fprintf (stderr,"\n in function beamel2d::beam_transf_matrix.\n");
    }
    
    //  local vector y_l
    tmat[0][1]=tmat[1][2]*tmat[2][0]-tmat[2][2]*tmat[1][0];
    tmat[1][1]=tmat[2][2]*tmat[0][0]-tmat[0][2]*tmat[2][0];
    tmat[2][1]=tmat[0][2]*tmat[1][0]-tmat[1][2]*tmat[0][0];

    //  norm of the local vector y_l
    c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));

    if (c<Mp->zero){
      fprintf (stderr,"\n\n zero y base vector of the %ld beamgen3d element",eid);
      fprintf (stderr,"\n in function beamgen3d::beam_transf_matrix (file %s, line %d).\n",__FILE__,__LINE__);
    }
    
    //  normed local vector y_l
    tmat[0][1]=tmat[0][1]/c;
    tmat[1][1]=tmat[1][1]/c;
    tmat[2][1]=tmat[2][1]/c;
    
  }
  else{
    
    //  normed local vector y_l
    tmat[0][1]=tmat[0][1]/c;
    tmat[1][1]=tmat[1][1]/c;
    tmat[2][1]=tmat[2][1]/c;
    
    //  local vector z_l (defined by vector product)
    tmat[0][2]=tmat[1][0]*tmat[2][1]-tmat[2][0]*tmat[1][1];
    tmat[1][2]=tmat[2][0]*tmat[0][1]-tmat[0][0]*tmat[2][1];
    tmat[2][2]=tmat[0][0]*tmat[1][1]-tmat[1][0]*tmat[0][1];
    
    //  norm of the local vector z_l
    c=sqrt((tmat[0][2]*tmat[0][2]+tmat[1][2]*tmat[1][2]+tmat[2][2]*tmat[2][2]));
    
    if (c<Mp->zero){
      fprintf (stderr,"\n\n zero z base vector of the %ld beamgen3d element",eid);
      fprintf (stderr,"\n in function beamel2d::beam_transf_matrix.\n");
    }
    
    //  normed local vector z_l
    tmat[0][2]=tmat[0][2]/c;
    tmat[1][2]=tmat[1][2]/c;
    tmat[2][2]=tmat[2][2]/c;
  }
}
void beamgen3d::ck_matrix (matrix &ck,  double s,double c,double eh,double dl)
{
  double eh1,c1,cjm;

  fillm (0.0,ck);
  eh1=eh*dl;
  
// SESTAVENI MATICE Ck konstant
      c1=c-1.;
      cjm=c1*c1-(s-eh1)*s;

      ck[0][3]=s/cjm;
      ck[0][2]=-ck[0][3]*c1/s;
      ck[0][1]=-eh*ck[0][3];
      ck[0][0]=1.-ck[0][2];

      ck[1][3]= (eh1*s-c1)/eh/cjm;
      ck[1][2]= -(1./eh+ck[1][3]*c1)/s;
      ck[1][1]= 1.-eh*ck[1][3];
      ck[1][0]=-ck[1][2];

      ck[2][3]=-ck[0][3];
      ck[2][2]=-ck[2][3]*c1/s;
      ck[2][1]=-eh*ck[2][3];
      ck[2][0]=-ck[2][2];

      ck[3][3]= c1/eh/cjm;
      ck[3][2]=-ck[3][3]*(s-eh1)/c1;
      ck[3][1]=-eh*ck[3][3];
      ck[3][0]=-ck[3][2];
 }

/**
   function assembles geometric %matrix
   (function assembles %matrix of derivatives of approximation functions)
   
   
   ordering of approximation functions

   du/dx - strain e_xx
   phi+dv/dx - strain e_yz
   phi+dw/dx - strain e_xz
   dphix/dx - curvature
   dphiy/dx - curvature
   dphiz/dx - curvature

   @param n - array containing geometric %matrix
   @param s - natural coordinate from segment <0;1>
   @param dl - length of the beam
   @param gy - 6.0*E*I_y/k/G/A/l/l
   @param gz - 6.0*E*I_z/k/G/A/l/l
   
   PF, 20.9.2006
*/
void beamgen3d::geom_matrix (matrix &n,double s,double dl,double gy,double gz)
{
  double aj1,ll;

  fillm (0.0,n);
  ll=dl*dl;
  
//    Vnitrni sily   N, Vy, Vz, Mx, My, Mz
     n[0][0]=-1./dl;
     n[0][6]= 1./dl;
//  My
     aj1=1./(1.+2*gy);
     n[4][2] = aj1*(6.-12.*s)/ll;
     n[4][4] = aj1*(-2.*(2.+gy)+6.*s)/dl;
     n[4][8] =-n[4][2];
     n[4][10]= aj1*(-2.*(1.-gy)+6.*s)/dl;
// Qz
     n[2][2] =-2.*gy/dl/(1.+2.*gy);
     n[2][4] = gy/(1.+2.*gy);
     n[2][8] =-n[2][2];
     n[2][10]= n[2][4];
// Mz
     aj1=1./(1.+2.*gz);
     n[5][1] =-aj1*( 6.-12.*s)/ll;
     n[5][5] = aj1*(-2.*(2.+gz)+6.*s)/dl;
     n[5][7] =-n[5][1];
     n[5][11]= aj1*(-2.*(1.-gz)+6.*s)/dl;
// Qy
     n[1][1] =-2.*gz/dl/(1.+2.*gz);
     n[1][5] =-gz/(1.+2.*gz);
     n[1][7] =-n[1][1];
     n[1][11]= n[1][5];
// Mx
     n[3][3] =-1./dl;
     n[3][9] = 1./dl;
 }

/**
   function assembles %matrix of approximation functions
   
   ordering of approximation functions
   u v w- displacement along x,y,z
   fix,fiy,fiz - rotation around x,y,z
   
   @param n - %matrix of approximation functions
   @param s - natural coordinate from segment <0;1>
   @param dl - length of the beam
   @param gy - 6.0*E*I_y/k/G/A/l/l
   @param gz - 6.0*E*I_z/k/G/A/l/l
   
   PF, 20.9.2006
*/
void beamgen3d::bf_matrix (matrix &n,double s,double dl,double gy,double gz)
{
  double ll,aj1;
  ll=dl*dl;
//  u
        n[0][0] = 1.-s;
        n[0][6] = s;
//  w
        aj1=1./(1.+2.*gy);
        n[2][2] =aj1*(1.+2.*gy-2.*gy*s-3.*s*s+2*s*s*s);
        n[2][4] =aj1*(-(1.+gy)*s+(2.+gy)*s*s-s*s*s)*dl;
        n[2][8] =aj1*(2.*gy*s+3.*s*s-2.*s*s*s);
        n[2][10]=aj1*(gy*s+(1.-gy)*s*s-s*s*s)*dl;
// fy
        n[4][2] = aj1*(6.*s-6.*s*s)/dl;
        n[4][4] = aj1*(1.+2.*gy-2.*(2.+gy)*s+3.*s*s);
        n[4][8] =-n[4][2];
        n[4][10]= aj1*(-2.*(1.-gy)*s+3.*s*s);

//  v
        aj1=1./(1.+2.*gz);
        n[1][1] = aj1*(1.+2.*gz-2.*gz*s-3.*s*s+2.*s*s*s);
        n[1][5] =-aj1*(-(1.+gz)*s+(2.+gz)*s*s-s*s*s)*dl;
        n[1][7] = aj1*(2.*gz*s+3.*s*s-2.*s*s*s);
        n[1][11]=-aj1*(gz*s+(1.-gz)*s*s-s*s*s)*dl;
// fz
        n[5][1] =-aj1*(6.*s-6.*s*s)/dl;
        n[5][5] = aj1*(1.+2.*gz-2.*(2.+gz)*s+3.*s*s);
        n[5][7] =-n[5][1];
        n[5][11]= aj1*(-2.*(1.-gz)*s+3.*s*s);
// fx
        n[3][3] = 1.-s;
        n[3][9] = s;
  
}






/**
   function computes stiffness %matrix of 3D beam element
   
   @param eid - number of element
   @param sm - stiffness %matrix
   
   PF 20.9.2006
*/
void beamgen3d::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
  
  long i;
  for (i=0;i<sm.m;i++){
    if (sm[i][i]<Mp->zero){
      fprintf (stderr,"\n\n nulovy diagonalni prvek v matici tuhosti jednoho prvku (file %s, line %d)",__FILE__,__LINE__);
    }
  }  
}

/**
   function computes stiffness %matrix of 3D beam element
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   
   PF, 20.9.2006
*/
void beamgen3d::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long ipp,i,transf;
  double e,g,a,*ixyz,*ioyz,*beta,dl,ll,gy,gz,aj1, integn4;
  ivector nodes(nne);
  vector vec(3),x(nne),y(nne),z(nne),integn1(4),integn2(3),integn3(2);
  matrix d(tncomp,tncomp),tran(3,3);
  
  ixyz = new double [3];
  ioyz = new double [3];
  beta = new double [2];

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  ll=dl*dl;
  
  fillm (0.0,sm);
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  // vysecovy Io, Ioy, Ioz
//!!!  Mc->give_iomega (eid,ioyz);
   ioyz[0]=4.2692e-6 ; ioyz[1]=-2.6687e-5;  ioyz[2]=0.0; ixyz[0]=2.646e-7; ixyz[1]=2.135e-4; ixyz[2]=3.333e-5; 
 ioyz[0]=1.1924e-6 ; ioyz[1]=0.0;  ioyz[2]=7.5625e-6; ixyz[0]=2.333e-7; ixyz[1]=3.048e-5; ixyz[2]=6.25e-5; 
//  ioyz[0]=1.1924e-6 ; ioyz[1]=7.5625e-6;  ioyz[2]=0.0;  ixyz[0]=2.333e-7; ixyz[1]=6.25e-5; ixyz[2]=3.048e-5; 
//  ioyz[0]=0.2773e-6 ; ioyz[1]=0.0;  ioyz[2]=0.0;
//  ioyz[0]=1. ; ioyz[1]=0.;  ioyz[2]=0.0; ixyz[0]=1.; ixyz[1]=1.; ixyz[2]=1.;

  e=d[0][0];  g=d[1][1];

//  stiffness matrix in local coordinate system
// axial forces and torsion
  sm[0][0]=  e*a/dl; sm[0][7]= -sm[0][0]; sm[7][7]= sm[0][0];
  //sm[3][3]=  g*ixyz[0]/dl; sm[3][10]= -sm[3][3]; sm[10][10]= sm[3][3];
//  direct s Iy   w, fy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
  gy=0.;gz=0.;  
  aj1=ixyz[1]*e/(1.+2*gy);
  integn1[0]=12./ll/dl*aj1;      integn1[1]= 6./ll*aj1; integn1[2]=-12./ll/dl*aj1; integn1[3]= 6./ll*aj1;
  integn2[0]=2.*(2.+gy)/dl*aj1;  integn2[1]=-6./ll*aj1; integn2[2]=2.*(1.-gy)/dl*aj1;
  integn3[0]=12./ll/dl*aj1;      integn3[1]=-6./ll*aj1;
  integn4   =2.*(2.+gy)/dl*aj1;
  sm[2][2] = integn1[0];
  sm[2][4] = -integn1[1];
  sm[2][9] = integn1[2];
  sm[2][11]= -integn1[3];
  sm[4][4] = integn2[0];
  sm[4][9] = integn2[1];
  sm[4][11]= -integn2[2];
  sm[9][9] = integn3[0];
  sm[9][11]= -integn3[1];
  sm[11][11]= integn4;
//  direct s Iz   v, fz
  aj1=ixyz[2]*e/(1.+2.*gz);
  integn1[0]=12./ll/dl*aj1;      integn1[1]= 6./ll*aj1; integn1[2]=-12./ll/dl*aj1; integn1[3]=6./ll*aj1;
  integn2[0]=2.*(2.+gz)/dl*aj1;  integn2[1]=-6./ll*aj1; integn2[2]=2.*(1.-gz)/dl*aj1;
  integn3[0]=12./ll/dl*aj1;      integn3[1]=-6./ll*aj1;
  integn4   =2.*(2.+gz)/dl*aj1;
  sm[1][1] = integn1[0];
  sm[1][5] = integn1[1];
  sm[1][8] = integn1[2];
  sm[1][12]= integn1[3];
  sm[5][5] = integn2[0];
  sm[5][8] = integn2[1];
  sm[5][12]= integn2[2];
  sm[8][8] = integn3[0];
  sm[8][12]= integn3[1];
  sm[12][12]= integn4;
  
  stiffness_matrixtor (sm, dl,e,g,gy,gz,ixyz,ioyz);

  for (i=0;i<ndofe;i++){for (int j=i+1;j<ndofe;j++){sm[j][i]=sm[i][j];}}

//  transformation of stiffness matrix to the global system
// !!!zmenit   lgmatrixtransf3dblock (sm,tran);

  //  transformation of stiffness matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
  fprintf (Out,"\n\n MT prvek beam cislo %ld",eid);
  for (i=0;i<ndofe;i++){
     fprintf (Out,"\n %4ld",i);
     for (int j=0;j<ndofe;j++){
		fprintf (Out," %15.5e",sm[i][j]);
     }
  }

  delete [] ixyz;
  delete [] ioyz;
  delete [] beta;
}

/**
   function computes consistent mass %matrix of 3D beam element
   influence of inertial forces from rotations is accounted
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param mm - mass %matrix
   
   PF, 20.9.2006
*/
// *******************************************************
//    TORSION -   exact solver

  void beamgen3d::stiffness_matrixtor (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz)
{
  int i,i1,i2,j;
  double  eh,eh1,s,c,c1,s1,s2,c2,cs,pom1,pom2,pom3,xc,xs;
  vector sioyz(2),sm1(16),sm2(16);
  matrix ck(4,4);
  

//     ioyz[0] - I omega - vys. moment setrvac,   ioyz[1] = omegay,  ioyz[2] = omegaz
//     ixyz[0] - I k  - moment tuhosti ve volnem krouceni

   eh=pow(g*ixyz[0]/e/ioyz[0],0.5); 
//   eh=1.; 
   eh1=eh*dl;
   s=sinh(eh1); c=cosh(eh1); 
   ck_matrix (ck,  s,c,eh,dl);

//    Integraly  cosh, sinh, x.cosh, x.sinh, cosh^2, cosh*sinh, sinh^2
   c1=s/eh; s1=(c-1.)/eh; xc = (dl*s-s1)/eh; xs = (dl*c-c1)/eh; c2=s*c/2/eh+dl/2; cs=s*s/2/eh;  s2=s*c/2/eh-dl/2;
   pom1=e*ioyz[0]*eh*eh*eh*eh;
   pom2=g*ixyz[0];
//   pom2=pom1/eh/eh;
   pom3=pom1/eh/eh;
   for  (i=0;i<4;i++){
       i1=4*i-i*(i+1)/2; //1=4*(i-1)-i*(i-1)/2;
       for  (j=i;j<4;j++){
         i2=j+i1;
         sm1[i2]=(ck[i][2]*ck[j][2]*c2+ck[i][3]*ck[j][3]*s2+
              (ck[i][3]*ck[j][2]+ck[i][2]*ck[j][3])*cs)*pom1;
         sm1[i2]=sm1[i2]+(ck[i][1]*ck[j][1]*dl+eh*(ck[i][1]*ck[j][2]+
               ck[i][2]*ck[j][1])*s1+eh*(ck[i][1]*ck[j][3]+ck[i][3]*ck[j][1])*c1)*pom2;
         sm1[i2]=sm1[i2]+((ck[i][2]*ck[j][3]+ck[i][3]*ck[j][2])*cs+
			   ck[i][2]*ck[j][2]*s2+ck[i][3]*ck[j][3]*c2)*eh*eh*pom2;
		}
   }

   sm[3][3]  = sm1[0] ;
   sm[3][6]  = sm1[1] ;
   sm[3][10] = sm1[2] ;
   sm[3][13] = sm1[3] ;
   sm[6][6]  = sm1[4] ;
   sm[6][10] = sm1[5] ;
   sm[6][13] = sm1[6] ;
   sm[10][10]= sm1[7] ;
   sm[10][13]= sm1[8] ;
   sm[13][13]= sm1[9] ;

//   gy=0.;    gz=0.;
//   druhe der w
//   n[0] = (6.-12.*s)/ll; n[1]= (-2.*(2.+gy)+6.*s)/dl; n[2]=-(6.-12.*s)/ll; n[3]= (-2.*(1.-gy)+6.*s)/dl;
   
   for  (i=0;i<4;i++){
        sm1[i]   = 6.*eh*eh/dl/dl*(-(ck[i][2]*c1+ck[i][3]*s1)+
                 2./dl*(ck[i][2]*xc+ck[i][3]*xs) );
        sm1[i+4] =-2.*eh*eh/dl*(-(2.+gy)*(ck[i][2]*c1+ck[i][3]*s1)+
                 3./dl*(ck[i][2]*xc+ck[i][3]*xs) );
        sm2[i+4] =-2.*eh*eh/dl*(-(2.+gz)*(ck[i][2]*c1+ck[i][3]*s1)+
                 3./dl*(ck[i][2]*xc+ck[i][3]*xs) );
        sm1[i+8] =-sm1[i];
        sm1[i+12]=-2.*eh*eh/dl*(-(1.-gy)*(ck[i][2]*c1+ck[i][3]*s1)+
                 3./dl*(ck[i][2]*xc+ck[i][3]*xs) );
        sm2[i+12]=-2.*eh*eh/dl*(-(1.-gz)*(ck[i][2]*c1+ck[i][3]*s1)+
                 3./dl*(ck[i][2]*xc+ck[i][3]*xs) );
   }

   sioyz[0]=ioyz[1]*e/(1.+2.*gy);
   sioyz[1]=ioyz[2]*e/(1.+2.*gz);
   for  (i=0;i<2;i++){
        i1=1-i;
        i2=i+1;
        sm[i2][3] = sm1[0]*sioyz[i1];
        sm[i2][6] = sm1[1]*sioyz[i1];
        sm[i2][10]= sm1[2]*sioyz[i1];
        sm[i2][13]= sm1[3]*sioyz[i1];
        sm[3][4]  = sm1[4]*sioyz[0];
        sm[3][5]  =-sm2[4]*sioyz[1];
        sm[3][i+8]= sm1[8]*sioyz[i1];
        sm[3][11] = sm1[12]*sioyz[0];
        sm[3][12] =-sm2[12]*sioyz[1];

        sm[4][6]  = sm1[5]*sioyz[0];
        sm[5][6]  =-sm2[5]*sioyz[1];
        sm[4][10] = sm1[6]*sioyz[0];
        sm[5][10] =-sm2[6]*sioyz[1];
        sm[4][13] = sm1[7]*sioyz[0];
        sm[5][13] =-sm2[7]*sioyz[1];
        sm[6][i+8]= sm1[9]*sioyz[i1];
        sm[6][11] = sm1[13]*sioyz[0];
        sm[6][12] =-sm2[13]*sioyz[1];
        i2=i+8;
        sm[i2][10]= sm1[10]*sioyz[i1];
        sm[i2][13]= sm1[11]*sioyz[i1];
        sm[10][11]= sm1[14]*sioyz[0];
        sm[10][12]=-sm2[14]*sioyz[1];

        sm[11][13]= sm1[15]*sioyz[0];
        sm[12][13]=-sm2[15]*sioyz[1];
   }

}

// *******************************************************
//     TORSION      1        cubic function

void beamgen3d::stiffness_matrixtor1 (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz)
{
  double aj1,ll,integn4;
  vector sioyz(2),integn1(4),integn2(3),integn3(2);

//     ioyz[0] - I omega - vys. moment setrvac,   ioyz[1] = omegay,  ioyz[2] = omegaz
//     ixyz[0] - I k  - moment tuhosti ve volnem krouceni

//  gy=0.; gz=0.;
  ll=dl*dl; aj1=e*ioyz[0];
  integn1[0]=12./ll/dl*aj1;      integn1[1]= 6./ll*aj1; integn1[2]=-12./ll/dl*aj1; integn1[3]= 6./ll*aj1;
  integn2[0]=2.*(2.+0.)/dl*aj1;  integn2[1]=-6./ll*aj1; integn2[2]=2.*(1.-0.)/dl*aj1;
  integn3[0]=12./ll/dl*aj1;      integn3[1]=-6./ll*aj1;
  integn4=2.*(2.+0.)/dl*aj1;
//first deriv. +, +, -, +      +,- ,-,   +,-,  +, 
  sm[3][3]  = integn1[0]+g*1.2*ixyz[0]/dl;
  sm[3][6]  = integn1[1]+g*0.1*ixyz[0];
  sm[3][10] = integn1[2]-g*1.2*ixyz[0]/dl;
  sm[3][13] = integn1[3]+g*0.1*ixyz[0];
  sm[6][6]  = integn2[0]+g*4./30.*ixyz[0]*dl;
  sm[6][10] = integn2[1]-g*0.1*ixyz[0];
  sm[6][13] = integn2[2]-g*1./30.*ixyz[0]*dl;
  sm[10][10]= integn3[0]+g*1.2*ixyz[0]/dl;
  sm[10][13]= integn3[1]-g*0.1*ixyz[0];
  sm[13][13]= integn4   +g*4./30.*ixyz[0]*dl;

//  direct s Ioy - w (index 2,9), fy (index 4,11)
  aj1=ioyz[1]*e/(1.+2.*gy);
  integn1[0]=12./ll/dl*aj1;      integn1[1]= 6./ll*aj1;  integn1[2]=-12./ll/dl*aj1; integn1[3]= 6./ll*aj1;
  integn2[0]=2.*(2.+gy)/dl*aj1;  integn2[1]=-6./ll*aj1;  integn2[2]=2.*(1.-gy)/dl*aj1;
  integn3[0]=12./ll/dl*aj1;      integn3[1]=-6./ll*aj1;
  integn4=2.*(2.+gy)/dl*aj1;

   sm[2][3] = integn1[0];
   sm[2][6] = integn1[1];
   sm[2][10]= integn1[2];
   sm[2][13]= integn1[3];
    sm[3][4] = -integn1[1];
   sm[3][8]=integn1[2];
    sm[3][11] = -integn1[3];
    sm[4][6]  = -integn2[0];
    sm[4][10] = -integn2[1];
    sm[4][13] = -integn2[2];
   sm[6][9]= integn2[1];
    sm[6][11] = -integn2[2];
   sm[9][10]= integn3[0];
   sm[9][13]= integn3[1];
    sm[10][11]= -integn3[1];
    sm[11][13]= -integn4;

//  direct s Ioz - v (index 1,8), fz (index 5,12)
   aj1=ioyz[2]*e/(1.+2.*gz);
   integn1[0]=12./ll/dl*aj1;      integn1[1]= 6./ll*aj1; integn1[2]=-12./ll/dl*aj1; integn1[3]= 6./ll*aj1;
   integn2[0]=2.*(2.+gz)/dl*aj1;  integn2[1]=-6./ll*aj1; integn2[2]=2.*(1.-gz)/dl*aj1;
   integn3[0]=12./ll/dl*aj1;      integn3[1]=-6./ll*aj1;
   integn4=2.*(2.+gz)/dl*aj1;
   sm[1][3] = integn1[0];
   sm[1][6] = integn1[1];
   sm[1][10]= integn1[2];
   sm[1][13]= integn1[3];
    sm[3][5] = integn1[1];
   sm[3][8]=integn1[2];
    sm[3][12] = integn1[3];
    sm[5][6]  = integn2[0];
    sm[5][10] = integn2[1];
    sm[5][13] = integn2[2];
   sm[6][8]= integn2[1];
    sm[6][12] = integn2[2];
   sm[8][10]= integn3[0];
   sm[8][13]= integn3[1];
    sm[10][12]= integn2[1];
    sm[12][13]= integn4;

 /* 
  for (i=0;i<2;i++){
      i1=1-i;
      i2=i+1;
      sm[i2][3] = 12.*sioyz[i1]/dl;
      sm[i2][6] = 6.*sioyz[i1];
      sm[i2][10]= -sm[i2][3];
      sm[i2][13]=  sm[i2][6];
      sm[3][4]  =-6.*sioyz[0];
      sm[3][5]  = 6.*sioyz[1];
      sm[3][i+8]=-12.*sioyz[i1]/dl;
      sm[3][11] =-6.*sioyz[0];
      sm[3][12] = 6.*sioyz[1];
      sm[4][6]  =-4.*sioyz[0]*dl;
      sm[5][6]  = 4.*sioyz[1]*dl;
      sm[4][10] = 6.*sioyz[0];
      sm[5][10] =-6.*sioyz[1];
      sm[4][13] =-2.*sioyz[0]*dl;
      sm[5][13] = 2.*sioyz[1]*dl;
      sm[6][i+8]=-6.*sioyz[i1];
      sm[6][11] =-2.*sioyz[0]*dl;
      sm[6][12] = 2.*sioyz[1]*dl;
      i2=i+8;
      sm[i2][10]= 12.*sioyz[i1]/dl;
      sm[i2][13]= -6.*sioyz[i1];
      sm[10][11]= 6.*sioyz[0];
      sm[10][12]=-6.*sioyz[1];
      sm[11][13]=-4.*sioyz[0]*dl;
      sm[12][13]= 4.*sioyz[1]*dl;
  }
  */
}

// *******************************************************
//     TORSION      2        kubic function Mindlin element
//
void beamgen3d::stiffness_matrixtor2 (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz,double *iro)
{
  double  kakr,ajy,ajz,ajk;

//  Razeni v matici tuhosti

//  N1, Vy1, Vz1, Mx1, My1, Mz1, B1
//  N2, Vy2, Vz2, Mx2, My2, Mz2, B2

//  K R U T   -  Mx, Qz, My
//     ioyz[0] - I omega - vys. moment setrvac,   ioyz[1] = omegay,  ioyz[2] = omegaz
//     ixyz[0] - I k  - moment tuhosti ve volnem krouceni
//     iro[0] -  Iro   Irosin  Irocos
    kakr  =6.*ioyz[0]*e/iro[0]/g/dl/dl;
    
	ajy= (1.+2.*gy);
    ajk= (1.+2.*kakr);
    ajy=1.;
    sm[2][3] = (12.*ioyz[1]*e/dl/dl/dl+4.*kakr*gy*iro[1]*g/dl)/ajk/ajy;
    sm[3][3] = (12.*ioyz[0]*e/dl/dl/dl+(4.*kakr*kakr*iro[0] 
		      +ixyz[0]*(4.*kakr*(kakr+1.)+1.2) )*g/dl )/ajk/ajk;
    sm[3][4] = (-6.*ioyz[1]*e/dl/dl-2*kakr*gy*iro[1]*g)/ajk/ajy;
    sm[2][6] = sm[4][5];
    sm[3][6] = (-6.*ioyz[0]*e/dl/dl-(2.*kakr*kakr*iro[0]+iro[0]/10.)*g )/ajk/ajk;
    sm[4][6] = ((2.+kakr+gy+2.*kakr*gy)*2.*ioyz[1]*e/dl
		       +kakr*gy*iro[1]*g*dl)/ajk/ajy;
    sm[6][6] = ((1.+(1.+kakr)*kakr)*4.*ioyz[0]*e/dl+(kakr*kakr*iro[0]
		       +(kakr*(kakr+1)+0.4)/3.*ixyz[0])*g*dl)/ajk/ajk;

    sm[3][9] = -sm[2][3]; sm[6][9] = -sm[3][4];
    sm[2][10] = -sm[2][3]; sm[3][10] = -sm[3][3];
    sm[4][10] = -sm[3][4]; sm[6][10] = -sm[3][6];
    sm[9][10] = sm[2][3]; sm[10][10]=  sm[3][3];

    sm[3][11] =  sm[3][4];
    sm[6][11] = ((1.-kakr-gy-2.*kakr*gy)*2.*ioyz[1]*e/dl
                 -kakr*gy*iro[1]*g*dl)/ajk/ajy;
    sm[10][11]= -sm[3][4]; sm[2][13] =  sm[2][6];
    sm[3][13] =  sm[3][6];
    sm[4][13] = ((1.-kakr-gy-2.*kakr*gy)*2.*ioyz[1]*e/dl
                +kakr*gy*iro[1]*g*dl)/ajk/ajy;
    sm[6][13] = ((2.-(4.+4.*kakr)*kakr)*ioyz[0]*e/dl+(kakr*kakr*iro[0]
                +(kakr*kakr-0.1)/3.*ixyz[0] )*g*dl)/ajk/ajk;
    sm[9][13] = -sm[2][6];sm[10][13] = -sm[3][6];
    sm[11][13] =  sm[4][6];sm[13][13] =  sm[6][6];

//  K R U T   -  Mx, Qy, Mz
      ajz= (1.+2.*gz);
      ajz=1.;
      sm[1][3] = (12.*ioyz[2]*e/dl/dl/dl+4.*kakr*gz*iro[2]*g/dl)/ajk/ajz;
      sm[3][5] = (-6.*ioyz[2]*e/dl/dl-2.*kakr*gz*iro[2]*g)/ajk/ajz;
      sm[1][6] = sm[3][5];
      sm[5][6] = ((2.+kakr+gz+2.*kakr*gz)*2.*ioyz[2]*e/dl
                +kakr*gz*iro[2]*g*dl)/ajk/ajz;

      sm[3][8]  = -sm[1][3];sm[6][8]  = -sm[3][5];
      sm[1][10] = -sm[1][3];sm[5][10] = -sm[3][5];
      sm[8][10] =  sm[1][3];

      sm[3][12] =  sm[3][5];
      sm[6][12] = ((1.-kakr-gz-2.*kakr*gz)*2.*ioyz[2]*e/dl
                -kakr*gz*iro[2]*g*dl)/ajk/ajz;
      sm[10][12]= -sm[4][6];sm[2][14] =  sm[2][7];
      sm[5][13] = ((1.-kakr-gz-2.*kakr*gz)*2.*ioyz[2]*e/dl
                +kakr*gz*iro[2]*g*dl)/ajk/ajz;
      sm[8][13] = -sm[1][6];
      sm[12][13]=  sm[5][6];

}


/**
   function serves for computation of nodal forces and moments caused by body load
   
   @param eid - element id
   @param lm - load %matrix
   
   PF, 1.9.2006
*/
void beamgen3d::load_matrix (long eid,matrix &lm)
{
  long i,ii,transf;
  double jac,a,ixyz[3],beta[2],e,g,gy,gz,l;
  ivector nodes (nne);
  vector vec(3),x(nne),y(nne),z(nne),w(intordmm),gp(intordmm),xl(2);
  matrix d(tncomp,tncomp),n(napfun,ndofe),tran(3,3);
  
  //  nodal coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  coordinates of base vector z of the beam in global system
  Mc->give_vectorlcs (eid,vec);
  //  transformation matrix from local to global coordinates
  beam_transf_matrix (tran,l,vec,x,y,z,eid);
  //  area of cross section
  Mc->give_area (eid,a);
  //  shear coefficients for particular cross section
  Mc->give_shearcoeff (eid,beta);
  //  moments of inertia (I_x, I_y, I_z)
  Mc->give_mominer (eid,ixyz);
  //  number of integration point
  ii=Mt->elements[eid].ipp[0][0];
  //  stiffness matrix of material
  Mm->matstiff (d,ii);
  //  Young's modulus of leasticity
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  //  generalized shear modulus G_y
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/l;
  //  generalized shear modulus G_z
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/l;
  
  
  //  local coordinates of nodes for evaluation of jacobian
  xl[0]=0.0;
  xl[1]=l;
  
  //  element node numbers
  Mt->give_elemnodes (eid,nodes);
  
  //  integration points assembling
  gauss_points (gp.a,w.a,intordmm);
  
  //  matrix cleaning
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    //  computation of jacobian
    jac_1d (jac,xl,gp[i]);
    
    //  assembling of matrix of approximation functions
    bf_matrix (n,(gp[i]+1.0)/2.0,l,gy,gz);
    
    jac*=a*w[i];
    
    //  multiplication N^T N jac
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  
  //  transformation of load matrix to the global system
//   !!!zmenit   lgmatrixtransf3dblock (lm,tran);

  //  transformation of stiffness matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    //  assembling of transformation matrix from global to nodal system
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }
}


/**
   function stores nodal displacements and rotations in local coordinate system
   this is equivalent to functions computing strains
   nodal displacements and rotations are better quantities in case of beams
   
   @param eid - element id
   @param lcid - load case id
   
   PF, 20.9.2006
*/
void beamgen3d::nodal_displ (long eid,long lcid)
{
  long ii,transf;
  double l;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),d(ndofe),f(ndofe);
  matrix tm(3,3);
  
  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  element node numbers
  Mt->give_elemnodes (eid,nodes);

  //  transformation from nodal to global coordinate system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    //  assembling of transformation matrix
    transf_matrix (nodes,tmat);
    //  transformation of nodal displacements and rotations
    mtxv (tmat,d,f);
    copyv (f,d);
  }

  //  nodal coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  vector of z-axis
  Mc->give_vectorlcs (eid,vec);
  
  //  assembling of transformation matrix from global to local coordinate system
  beam_transf_matrix (tm,l,vec,x,y,z,eid);
  
  //  transformation from global to local coordinate system
//   !!!zmenit   glvectortransf3dblock (d,tm);

  //  number of integration point
  ii=Mt->elements[eid].ipp[0][0];
  
  
  Mm->storestrain (lcid,ii,0,7,d);
  Mm->storestrain (lcid,ii+1,7,7,d);
 
}

/**
   function computes nodal forces and moments expressed in local coodinate system
   this is equivalent to functions computing stresses
   nodal forces and moments are better quantities in case of beams

   x-axis goes through the beam
   z-axis is defined in cross section
   y-axis is defined by right hand orientation
   
   @param eid - element id
   @param lcid - load case id
   
   JK, 20.2.2002, revised 2.9.2006
*/
void beamgen3d::nodal_forces (long eid,long lcid)
{
  long ii,j,transf;
  double l;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),d(ndofe),f(ndofe);
  matrix sm(ndofe,ndofe),tm(3,3);
  
  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  assembling of stiffness matrix
  stiffness_matrix (eid,0,0,sm);
  
  //  computation of nodal forces and moments
   mxv (sm,d,f);
 	fprintf (Out,"\n\n MT prvek beam cislo %ld",lcid);
	fprintf (Out," \n u, v, w, fx, fy, fz, o  \n ");
	for (j=0;j<ndofe;j++){
		fprintf (Out," %15.5e",d[j]);
    }
	fprintf (Out," \n N,  Vy,  Vz, Mx, My, Mz, B,\n ");
	for (j=0;j<ndofe;j++){
		fprintf (Out," %15.5e",f[j]);
    }
 
  //  element node numbers
  Mt->give_elemnodes (eid,nodes);

  //  transformation from nodal to global coordinate system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    //  assembling of transformation matrix
    transf_matrix (nodes,tmat);
    //  transformation of nodal forces and moments
    mtxv (tmat,f,d);
    copyv (d,f);
  }
  
  //  nodal coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  vector of z-axis
  Mc->give_vectorlcs (eid,vec);
  
  //  assembling of transformation matrix from global to local coordinate system
  beam_transf_matrix (tm,l,vec,x,y,z,eid);

  internal_forces (lcid, eid, 0, 0, f);

  //  transformation from global to local coordinate system
  //glvectortransf3dblock (f,tm);
  
  //  number of integration point
  ii=Mt->elements[eid].ipp[0][0];
  
  //  storage of nodal forces and moments
  //  nodal forces and moments are expressed in local coordinate system
  Mm->storestress (lcid,ii,0,6,f);
  Mm->storestress (lcid,ii+1,6,6,f);
  
  
}

/**
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   PF 20.9.2006
*/
void beamgen3d::res_internal_forces (long lcid,long eid,vector &ifor)
{
   internal_forces1 (lcid, eid, 0, 0, ifor);
}
// *******************************************************
//     TORSION - Forces
void beamgen3d::internal_forces (long lcid,long eid,long /*ri*/,long /*ci*/,vector &/*ifor*/)
{
  long i,j,integr=2,ipp,transf;
  double dl,ll,a,ixyz[3],ioyz[3],beta[2],gy,gz,kakr,e,g;
  double eh,ehx,ajy,ajz,ajk,xs, s,c;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),rl(ndofe),rg(ndofe),f(ndofe),fx(tncomp),w(integr),gp(integr);
  vector dderno(4),dderny(4),ddernz(4);
  matrix n(tncomp,ndofe),d(tncomp,tncomp),tran(3,3),r(tncomp,2);
  matrix bb(tncomp,ndofe), ck(4,4);
  
// r  .... end forces last time
// m  .... dl(i)/dl(i-1)
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,rl.a);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  ll=dl*dl;
// r  .... end forces last time
// m  .... dl(i)/dl(i-1)
  fillm (0.0,r);
  
//  transformation of displacement vector
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //locglobtransf (rg,rl,tmat);
    lgvectortransf (rg,rl,tmat);
  }
  else{
    copyv (rl,rg);
  }
// *** transf. from gcs to lcs
//  globloctransf (rg,rl,tran);

//  stiffness_matrix (eid,ri,ci,sm);
//  mxv (sm,r,contr);

  ipp=Mt->elements[eid].ipp[0][0];
  Mm->matstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  // vysecovy Io, Ioy, Ioz
//!!!  Mc->give_iomega (eid,ioyz);
  ioyz[0]=4.2692e-6 ; ioyz[1]=-2.6687e-5;  ioyz[2]=0.0; ixyz[0]=2.646e-7; ixyz[1]=2.135e-4; ixyz[2]=3.333e-5; 
  ioyz[0]=1.1924e-6 ; ioyz[1]=0.0;  ioyz[2]=7.5625e-6; ixyz[0]=2.333e-7; ixyz[1]=3.048e-5; ixyz[2]=6.25e-5; 
//  ioyz[0]=1.1924e-6 ; ioyz[1]=7.5625e-6;  ioyz[2]=0.0;  ixyz[0]=2.333e-7; ixyz[1]=6.25e-5; ixyz[2]=3.048e-5; 
//  ioyz[0]=0.2773e-6 ; ioyz[1]=0.0;  ioyz[2]=0.0;
//  ioyz[0]=1. ; ioyz[1]=0.;  ioyz[2]=10.0; ixyz[0]=1.; ixyz[1]=10.; ixyz[2]=1.;

  e=d[0][0];  g=d[1][1];

//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
    gy=0.; gz=0.; 
	kakr=0.;

  eh=pow(g*ixyz[0]/e/ioyz[0],0.5); 
//  eh=1.; 
  eh=eh*dl; 
  s=sinh(eh); c=cosh(eh);
  eh=eh/dl; 
  ck_matrix (ck,  s,c,eh,dl);
  xs=0.0;
  for (i=0;i<2;i++){

    ehx=eh*xs*dl;
    s=sinh(ehx); c=cosh(ehx);


//    Vnitrni sily   N, My, Mz, B , Mvaz,  Mvol
    bb[0][0]=-e*a/dl;
    bb[0][7]= e*a/dl;
//        BB(i,4/7/11/14)=  e*pru(10+i)*eh*eh*(CK(i,3)*C+CK(i,4)*S)
    ajy= (1.+2.*gy); ajz= (1.+2.*gz); ajk= (1.+2.*kakr);
    dderno[0]=(-6.+12.*xs)/ll;  dderno[1]=(-2.*(2.+kakr)+6.*xs)/dl; dderno[2]=( 6.-12.*xs)/ll;  dderno[3]=(-2.*(1.-kakr)+6.*xs)/dl;
//    dderno[0]=eh*eh*(ck[0][2]*s+ck[0][3]*c);  dderno[1]=eh*eh*(ck[1][2]*s+ck[1][3]*c); dderno[2]=eh*eh*(ck[2][2]*s+ck[2][3]*c);  dderno[3]=eh*eh*(ck[3][2]*s+ck[3][3]*c);
    dderny[0]=(-6.+12.*xs)/ll;  dderny[1]=(-2.*(2.+gy)+6.*xs)/dl;   dderny[2]=( 6.-12.*xs)/ll;  dderny[3]=(-2.*(1.-gy)+6.*xs)/dl;
//  ex= - o O'' - y v'' - z w''
//  My = Int. z E ex = - E Iy w''- E Ioy O''=  E Iy fiy'- E Ioy O''
	bb[1][2]=  -e*ixyz[1]/ajy*dderny[0];
    bb[1][3]=  -e*ioyz[1]/ajk*dderno[0];
    bb[1][4]=   e*ixyz[1]/ajy*dderny[1];
    bb[1][6]=  -e*ioyz[1]/ajk*dderno[1];
    bb[1][9]=  -e*ixyz[1]/ajy*dderny[2];
    bb[1][10]= -e*ioyz[1]/ajk*dderno[2];
    bb[1][11]=  e*ixyz[1]/ajy*dderny[3];
    bb[1][13]= -e*ioyz[1]/ajk*dderno[3];
// Qz
    bb[5][2]=  e*ixyz[1]/ajy*12./ll/dl;
    bb[5][3]=  e*ioyz[1]/ajk*12./ll/dl;
	bb[5][4]= -e*ixyz[1]/ajy*6./ll;
    bb[5][6]=  e*ioyz[1]/ajk*6./ll;
    bb[5][9]=  e*ixyz[1]/ajy*(-12.)/ll/dl;
    bb[5][10]= e*ioyz[1]/ajk*(-12.)/ll/dl;
    bb[5][11]=-e*ixyz[1]/ajy*6./ll;
    bb[5][13]= e*ioyz[1]/ajk*6./ll;

// Mz=  -Int. y E ex = E Iz v''+ E Ioz O''=  E Iz fiz'+ E Ioy O''
    ddernz[0]=(-6.+12.*xs)/ll;  ddernz[1]=(-2.*(2.+gz)+6.*xs)/dl; ddernz[2]=( 6.-12.*xs)/ll;  ddernz[3]=(-2.*(1.-gz)+6.*xs)/dl;
    bb[2][1]=  e*ixyz[2]/ajz*ddernz[0];
    bb[2][3]=  e*ioyz[2]/ajk*dderno[0];
    bb[2][5]=  e*ixyz[2]/ajz*ddernz[1];
    bb[2][6]=  e*ioyz[2]/ajk*dderno[1];
    bb[2][ 8]= e*ixyz[2]/ajz*ddernz[2];
    bb[2][10]= e*ioyz[2]/ajk*dderno[2];
    bb[2][12]= e*ixyz[2]/ajz*ddernz[3];
    bb[2][13]= e*ioyz[2]/ajk*dderno[3];
// Qy
    bb[4][1]=  e*ixyz[2]/ajz*12./ll/dl;
    bb[4][3]=  e*ioyz[2]/ajk*12./ll/dl;
    bb[4][5]=  e*ixyz[2]/ajz*6./ll;
    bb[4][6]=  e*ioyz[2]/ajk*6./ll;
    bb[4][8] = e*ixyz[2]/ajz*(-12.)/ll/dl;
    bb[4][10]= e*ioyz[2]/ajk*(-12.)/ll/dl;
    bb[4][12]= e*ixyz[2]/ajz*6./ll;
    bb[4][13]= e*ioyz[2]/ajk*6./ll;

// B =  Int. o E ex = -E Io O'' - E Ioz v'' - E Ioy w''
    bb[3][1]= -e*ioyz[2]/ajz*ddernz[0];
//    bb[3][1]= -e*ioyz[2]/ajz*dderno[0];
    bb[3][2]= -e*ioyz[1]/ajy*dderny[0];
//    bb[3][2]= -e*ioyz[1]/ajy*dderno[0];
    bb[3][3]= -e*ioyz[0]*eh*eh*(ck[0][2]*c+ck[0][3]*s);
    bb[3][4]=  e*ioyz[1]/ajy*dderny[1];
//    bb[3][4]=  e*ioyz[1]/ajy*dderno[1];
    bb[3][5]= -e*ioyz[2]/ajz*ddernz[1];
//    bb[3][5]= -e*ioyz[2]/ajz*dderno[1];
    bb[3][6]= -e*ioyz[0]*eh*eh*(ck[1][2]*c+ck[1][3]*s);
    bb[3][8]= -e*ioyz[2]/ajz*ddernz[2];
//    bb[3][8]= -e*ioyz[2]/ajz*dderno[2];
    bb[3][9]= -e*ioyz[1]/ajy*dderny[2];
//    bb[3][9]= -e*ioyz[1]/ajy*dderno[2];
    bb[3][10]=-e*ioyz[0]*eh*eh*(ck[2][2]*c+ck[2][3]*s);
    bb[3][11]= e*ioyz[1]/ajy*dderny[3];
//    bb[3][11]= e*ioyz[1]/ajy*dderno[3];
    bb[3][12]=-e*ioyz[2]/ajz*ddernz[3];
//    bb[3][12]=-e*ioyz[2]/ajz*dderno[3];
    bb[3][13]=-e*ioyz[0]*eh*eh*(ck[3][2]*c+ck[3][3]*s);
   
//third deriv. +, +, -, + 
// Mvazane
    bb[6][1]= -e*ioyz[2]/ajz*12./ll/dl;
    bb[6][2]= -e*ioyz[1]/ajy*12./ll/dl;
    bb[6][3]= -e*ioyz[0]*eh*eh*eh*(ck[0][2]*s+ck[0][3]*c);
    bb[6][4]=  e*ioyz[1]/ajy*6./ll;
    bb[6][5]= -e*ioyz[2]/ajz*6./ll;
    bb[6][6]= -e*ioyz[0]*eh*eh*eh*(ck[1][2]*s+ck[1][3]*c);
    bb[6][8]= -e*ioyz[2]/ajz*(-12.)/ll/dl;
    bb[6][9]= -e*ioyz[1]/ajy*(-12.)/ll/dl;
    bb[6][10]=-e*ioyz[0]*eh*eh*eh*(ck[2][2]*s+ck[2][3]*c);
    bb[6][11]= e*ioyz[1]/ajy*6./ll;
    bb[6][12]=-e*ioyz[2]/ajz*6./ll;
    bb[6][13]=-e*ioyz[0]*eh*eh*eh*(ck[3][2]*s+ck[3][3]*c);
// Mvolne
    bb[7][3]=  g*ixyz[0]*(ck[0][1]+eh*(ck[0][2]*s+ck[0][3]*c));
    bb[7][6]=  g*ixyz[0]*(ck[1][1]+eh*(ck[1][2]*s+ck[1][3]*c));
    bb[7][10]= g*ixyz[0]*(ck[2][1]+eh*(ck[2][2]*s+ck[2][3]*c));
    bb[7][13]= g*ixyz[0]*(ck[3][1]+eh*(ck[3][2]*s+ck[3][3]*c));


    mxv (bb,rl,fx);

    xs=1.0;
	fprintf (Out,"\n\n MT prvek beam cislo %ld",lcid);
	fprintf (Out," \n u, v, w, fx, fy, fz, o  \n ");
	for (j=0;j<ndofe;j++){
		fprintf (Out," %15.5e",rl[j]);
    }
	fprintf (Out," \n N, My, Mz, B, Vy,  Vz, Mo  Mt\n ");
	for (j=0;j<tncomp;j++){
		fprintf (Out," %15.5e",fx[j]);
    }
  }

}
// *******************************************************
//     TORSION - Forces

void beamgen3d::internal_forces1 (long lcid,long eid,long /*ri*/,long /*ci*/,vector &/*ifor*/)
{
  long i,j,integr=2,ipp,transf;
  double dl,ll,a,ixyz[3],ioyz[3],beta[2],gy,gz,e,g;
  double  ajy,ajz,ajk, xs, kakr ;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),rl(ndofe),rg(ndofe),f(ndofe),w(integr),gp(integr), ddern(4),dderno(4);
  vector fx(4);
  matrix n(tncomp,ndofe),d(tncomp,tncomp),tran(3,3),r(tncomp,2);
  matrix bb(4,ndofe);
  
// r  .... end forces last time
// m  .... dl(i)/dl(i-1)
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,rl.a);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  ll=dl*dl;
// r  .... end forces last time
// m  .... dl(i)/dl(i-1)
  fillm (0.0,r);
  
//  transformation of displacement vector
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //locglobtransf (rg,rl,tmat);
    lgvectortransf (rg,rl,tmat);
  }
  else{
    copyv (rl,rg);
  }
// *** transf. from gcs to lcs
//  globloctransf (rg,rl,tran);

//  stiffness_matrix (eid,ri,ci,sm);
//  mxv (sm,r,contr);

  ipp=Mt->elements[eid].ipp[0][0];
  Mm->matstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  ioyz[0]=4.2692e-6 ; ioyz[1]=-2.6687e-5;  ioyz[2]=0.0; ixyz[0]=2.646e-7; ixyz[1]=2.135e-4; ixyz[2]=3.333e-5; 
//  ioyz[0]=1.1924e-6 ; ioyz[1]=0.0;  ioyz[2]=7.5625e-6; ixyz[0]=2.333e-7; ixyz[1]=3.048e-5; ixyz[2]=6.25e-5; 
//  ioyz[0]=1.1924e-6 ; ioyz[1]=7.5625e-6;  ioyz[2]=0.0;  ixyz[0]=2.333e-7; ixyz[1]=6.25e-5; ixyz[2]=3.048e-5; 
//  ioyz[0]=0.2773e-6 ; ioyz[1]=0.0;  ioyz[2]=0.0;
//  ioyz[0]=1. ; ioyz[1]=0.;  ioyz[2]=10.0; ixyz[0]=1.; ixyz[1]=10.; ixyz[2]=1.;

  e=d[0][0];  g=d[1][1];

//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;

  gy=0.; gz=0.; kakr=0.;

  xs=0.0;
  for (i=0;i<2;i++){
//    Vnitrni sily   N, My, B, Mz, B
    bb[0][0]=-e*a/dl;
    bb[0][7]= e*a/dl;

    ajy= (1.+2.*gy); ajz= (1.+2.*gz); ajk= (1.+2.*kakr);
    dderno[0]=(-6.+12.*xs)/ll; dderno[1]=(-2.*(2.+kakr)+6.*xs)/dl; dderno[2]=( 6.-12.*xs)/ll;  dderno[3]=(-2.*(1.-kakr)+6.*xs)/dl;
    ddern[0]=(-6.+12.*xs)/ll;  ddern[1]=(-2.*(2.+gy)+6.*xs)/dl;    ddern[2]=( 6.-12.*xs)/ll;   ddern[3]=(-2.*(1.-gy)+6.*xs)/dl;
//My
	bb[1][2]=  -e*ixyz[1]/ajy*ddern[0];
    bb[1][3]=  -e*ioyz[1]/ajk*dderno[0];
    bb[1][4]=   e*ixyz[1]/ajy*ddern[1];
    bb[1][6]=  -e*ioyz[1]/ajk*dderno[1];
    bb[1][9]=  -e*ixyz[1]/ajy*ddern[2];
    bb[1][10]= -e*ioyz[1]/ajk*dderno[2];
    bb[1][11]=  e*ixyz[1]/ajy*ddern[3];
    bb[1][13]= -e*ioyz[1]/ajk*dderno[3];
// B
	bb[3][2]=  -e*ioyz[1]/ajy*ddern[0];
    bb[3][3]=  -e*ioyz[0]/ajk*dderno[0];
    bb[3][4]=   e*ioyz[1]/ajy*ddern[1];
    bb[3][6]=  -e*ioyz[0]/ajk*dderno[1];
    bb[3][9]=  -e*ioyz[1]/ajy*ddern[2];
    bb[3][10]= -e*ioyz[0]/ajk*dderno[2];
    bb[3][11]=  e*ioyz[1]/ajy*ddern[3];
    bb[3][13]= -e*ioyz[0]/ajk*dderno[3];

//Mz
    ddern[0]=(-6.+12.*xs)/ll;  ddern[1]=(-2.*(2.+gz)+6.*xs)/dl; ddern[2]=( 6.-12.*xs)/ll;  ddern[3]=(-2.*(1.-gz)+6.*xs)/dl;
    bb[2][1]=  e*ixyz[2]/ajz*ddern[0];
    bb[2][3]=  e*ioyz[2]/ajk*dderno[0];
    bb[2][5]=  e*ixyz[2]/ajz*ddern[1];
    bb[2][6]=  e*ioyz[2]/ajk*dderno[1];
    bb[2][ 8]= e*ixyz[2]/ajz*ddern[2];
    bb[2][10]= e*ioyz[2]/ajk*dderno[2];
    bb[2][12]= e*ixyz[2]/ajz*ddern[3];
    bb[2][13]= e*ioyz[2]/ajk*dderno[3];
// B
    bb[3][1]=  -e*ioyz[2]/ajz*ddern[0];
    bb[3][5]=  -e*ioyz[2]/ajz*ddern[1];
    bb[3][8]=  -e*ioyz[2]/ajz*ddern[2];
    bb[3][12]= -e*ioyz[2]/ajz*ddern[3];

    mxv (bb,rl,fx);

    xs=1.0;
	fprintf (Out,"\n\n MT prvek beam cislo %ld",lcid);
	fprintf (Out," \n u, v, w, fx, fy, fz, o  \n ");
	for (j=0;j<ndofe;j++){
		fprintf (Out," %15.5e",rl[j]);
    }
	fprintf (Out," \n N, My, Mz, B\n ");
	for (j=0;j<4;j++){
		fprintf (Out," %15.5e",fx[j]);
    }
  }
}
