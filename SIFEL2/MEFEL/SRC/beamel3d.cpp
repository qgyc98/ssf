#include "beamel3d.h"
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


beamel3d::beamel3d (void)
{
  long i,j;
  
  //  number of nodes on element
  nne=2;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  //  in this case, displacements/rotations/forces/moments are used
  tncomp=6;
  //  number of functions approximated
  napfun=6;
  //  order of numerical integration of mass matrix
  intordmm=4;

  intordism=2;
  ssst=spacebeam;

  //  number of blocks (parts of geometric matrix)
  nb=1;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=1;
  
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

beamel3d::~beamel3d (void)
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
   Function approximates function defined by end nodal values.

   @param xi - coordinate on element
   @param nodval - nodal values
   
   @retval Function returns approximated value for the given natural coordinate xi.

   JK, 27.9.2005
*/
double beamel3d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));

  bf_lin_1d (bf.a,xi);

  scprd (bf,nodval,f);

  return f;
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
   
   PF, 20.12.2002, revised by JK 2.9.2006
*/
void beamel3d::geom_matrix (matrix &n,double s,double dl,double gy,double gz)
{
  double aj1,ll;

  fillm (0.0,n);
  ll=dl*dl;
  
  //    Vnitrni sily   N, Vy, Vz, Mx, My, Mz
  
  //  N
  n[0][0]=-1./dl;
  n[0][6]= 1./dl;

  //  Qy
  n[1][1] =-2.*gz/dl/(1.+2.*gz);
  n[1][5] =-gz/(1.+2.*gz);
  n[1][7] =-n[1][1];
  n[1][11]= n[1][5];

  //  Qz
  n[2][2] =-2.*gy/dl/(1.+2.*gy);
  n[2][4] = gy/(1.+2.*gy);
  n[2][8] =-n[2][2];
  n[2][10]= n[2][4];

  //  Mx
  n[3][3] =-1./dl;
  n[3][9] = 1./dl;
  
  //  My
  aj1=1./(1.+2*gy);
  n[4][2] = aj1*(6.-12.*s)/ll;
  n[4][4] = aj1*(-2.*(2.+gy)+6.*s)/dl;
  n[4][8] =-n[4][2];
  n[4][10]= aj1*(-2.*(1.-gy)+6.*s)/dl;

  //  Mz
  aj1=1./(1.+2.*gz);
  n[5][1] =-aj1*( 6.-12.*s)/ll;
  n[5][5] = aj1*(-2.*(2.+gz)+6.*s)/dl;
  n[5][7] =-n[5][1];
  n[5][11]= aj1*(-2.*(1.-gz)+6.*s)/dl;
}

/**
   function assembles %matrix of approximation functions
   
   ordering of approximation functions
   u v w - displacement along x,y,z
   fix,fiy,fiz - rotation around x,y,z
   
   @param n - %matrix of approximation functions
   @param s - natural coordinate from segment <0;1>
   @param dl - length of the beam
   @param gy - 6.0*E*I_y/k/G/A/l/l
   @param gz - 6.0*E*I_z/k/G/A/l/l
   
   PF, 20.12.2002, revised by JK 2.9.2006
*/
void beamel3d::bf_matrix (matrix &n,double s,double dl,double gy,double gz)
{
  double aj1;
  
  //  u
  n[0][0] = 1.-s;
  n[0][6] = s;

  //  v
  aj1=1./(1.+2.*gz);
  n[1][1] = aj1*(1.+2.*gz-2.*gz*s-3.*s*s+2.*s*s*s);
  n[1][5] =-aj1*(-(1.+gz)*s+(2.+gz)*s*s-s*s*s)*dl;
  n[1][7] = aj1*(2.*gz*s+3.*s*s-2.*s*s*s);
  n[1][11]=-aj1*(gz*s+(1.-gz)*s*s-s*s*s)*dl;

  //  w
  aj1=1./(1.+2.*gy);
  n[2][2] =aj1*(1.+2.*gy-2.*gy*s-3.*s*s+2*s*s*s);
  n[2][4] =aj1*(-(1.+gy)*s+(2.+gy)*s*s-s*s*s)*dl;
  n[2][8] =aj1*(2.*gy*s+3.*s*s-2.*s*s*s);
  n[2][10]=aj1*(gy*s+(1.-gy)*s*s-s*s*s)*dl;

  // fx
  n[3][3] = 1.-s;
  n[3][9] = s;

  // fy
  n[4][2] = aj1*(6.*s-6.*s*s)/dl;
  n[4][4] = aj1*(1.+2.*gy-2.*(2.+gy)*s+3.*s*s);
  n[4][8] =-n[4][2];
  n[4][10]= aj1*(-2.*(1.-gy)*s+3.*s*s);
  
  // fz
  n[5][1] =-aj1*(6.*s-6.*s*s)/dl;
  n[5][5] = aj1*(1.+2.*gz-2.*(2.+gz)*s+3.*s*s);
  n[5][7] =-n[5][1];
  n[5][11]= aj1*(-2.*(1.-gz)*s+3.*s*s);
}


/**
   function assembles %matrix of approximation functions for translational DOFs
   
   ordering of approximation functions
   u v w - displacement along x,y,z
   fix,fiy,fiz - rotation around x,y,z are neglected
   
   ordering of nodal values
   u_1, v_1, w_1, phi_x1, phi_y1, phi_z1, u_2, v_2, w_2, phi_x2, phi_y2, phi_z2

   @param n - %matrix of approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param gy - 6.0*E*I_y/k/G/A/l/l
   @param gz - 6.0*E*I_z/k/G/A/l/l
   
   JK, 16.11.2006
*/
void beamel3d::bf_matrix_transl (matrix &n,double xi,double l,double gy,double gz)
{
  fillm (0.0,n);

  //  displacement u (along x)
  n[0][0] = 1.0-xi;
  n[0][6] = xi;
  
  //  displacement v (along y)
  n[1][1] = (1.0 + 2.0*gz - 2.0*gz*xi - 3.0*xi*xi + 2.0*xi*xi*xi)/(1.0 + 2.0*gz);
  n[1][5] = ((1.0+gz)*xi - (2.0+gz)*xi*xi + xi*xi*xi)*l/(1.0+2.0*gz);
  n[1][7] = (2.0*gz*xi + 3.0*xi*xi - 2.0*xi*xi*xi)/(1.0+2.0*gz);
  n[1][11]= (-1.0*gz*xi - (1.0-gz)*xi*xi - xi*xi*xi)*l/(1.0+2.0*gz);
  
  //  displacement w (along z)
  n[2][2] = (1.0+2.0*gy - 2.0*gy*xi - 3.0*xi*xi + 2.0*xi*xi*xi)/(1.0+2.0*gy);
  n[2][4] = (-1.0*(1.0+gy)*xi + (2.0+gy)*xi*xi - xi*xi*xi)*l/(1.0+2.0*gy);
  n[2][8] = (2.0*gy*xi + 3.0*xi*xi - 2.0*xi*xi*xi)/(1.0+2.0*gy);
  n[2][10]= (gy*xi + (1.0-gy)*xi*xi - xi*xi*xi)*l/(1.0+2.0*gy);
}



/**
  The function assembles transformation %matrix from local nodal coordinate
  system to the global coordinate system x_g = T x_l.
   
  @param nodes - array containing node numbers
  @param tmat - transformation %matrix
  
  @return The function returns resulting transformation matrix in the parameter tmat. 
   
  PF, 20.12.2002, revised by JK, 28.11.2006
*/
void beamel3d::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*6+0][i*6] = Mt->nodes[nodes[i]].e1[0];
      tmat[i*6+1][i*6] = Mt->nodes[nodes[i]].e1[1];
      tmat[i*6+2][i*6] = Mt->nodes[nodes[i]].e1[2];

      tmat[i*6+0][i*6+1] = Mt->nodes[nodes[i]].e2[0];
      tmat[i*6+1][i*6+1] = Mt->nodes[nodes[i]].e2[1];
      tmat[i*6+2][i*6+1] = Mt->nodes[nodes[i]].e2[2];

      tmat[i*6+0][i*6+2] = Mt->nodes[nodes[i]].e3[0];
      tmat[i*6+1][i*6+2] = Mt->nodes[nodes[i]].e3[1];
      tmat[i*6+2][i*6+2] = Mt->nodes[nodes[i]].e3[2];
      
      
      tmat[i*6+3][i*6+3] = Mt->nodes[nodes[i]].e1[0];
      tmat[i*6+4][i*6+3] = Mt->nodes[nodes[i]].e1[1];
      tmat[i*6+5][i*6+3] = Mt->nodes[nodes[i]].e1[2];
      
      tmat[i*6+3][i*6+4] = Mt->nodes[nodes[i]].e2[0];
      tmat[i*6+4][i*6+4] = Mt->nodes[nodes[i]].e2[1];
      tmat[i*6+5][i*6+4] = Mt->nodes[nodes[i]].e2[2];

      tmat[i*6+3][i*6+5] = Mt->nodes[nodes[i]].e3[0];
      tmat[i*6+4][i*6+5] = Mt->nodes[nodes[i]].e3[1];
      tmat[i*6+5][i*6+5] = Mt->nodes[nodes[i]].e3[2];
    }
  }
}



/**
  The function assembles transformation %matrix from local element coordinate system 
  to global coordinate system x_g = T x_l.
   
  Columns of the transformation %matrix are created by coordinates
  of local base vectors expressed in global coordinate system.
   
  @param tmat - transformation %matrix
  @param dl - lenght of the beam
  @param vec - direction %vector of local axis z (in global coordinates)
  @param x,y,z - vectors of nodal coordinates in global system
  @param eid - element id
   
  JK, 3.2.2002, revised 1.9.2006
*/
void beamel3d::beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid)
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
    print_err("zero length of the %ld beamel3d element",__FILE__,__LINE__,__func__,eid);
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
      print_err("zero z base vector of the %ld beamel3d element",__FILE__,__LINE__,__func__,eid);
    }
    
    //  local vector y_l
    tmat[0][1]=tmat[1][2]*tmat[2][0]-tmat[2][2]*tmat[1][0];
    tmat[1][1]=tmat[2][2]*tmat[0][0]-tmat[0][2]*tmat[2][0];
    tmat[2][1]=tmat[0][2]*tmat[1][0]-tmat[1][2]*tmat[0][0];

    //  norm of the local vector y_l
    c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));

    if (c<Mp->zero){
      print_err("zero y base vector of the %ld beamel3d element",__FILE__,__LINE__,__func__,eid);
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
      print_err("zero z base vector of the %ld beamel3d element",__FILE__,__LINE__,__func__,eid);
    }
    
    //  normed local vector z_l
    tmat[0][2]=tmat[0][2]/c;
    tmat[1][2]=tmat[1][2]/c;
    tmat[2][2]=tmat[2][2]/c;
  }
  
  tmat[3][3]  = tmat[0][0];  tmat[3][4]   = tmat[0][1];  tmat[3][5]   = tmat[0][2];
  tmat[4][3]  = tmat[1][0];  tmat[4][4]   = tmat[1][1];  tmat[4][5]   = tmat[1][2];
  tmat[5][3]  = tmat[2][0];  tmat[5][4]   = tmat[2][1];  tmat[5][5]   = tmat[2][2];

  tmat[6][6]  = tmat[0][0];  tmat[6][7]   = tmat[0][1];  tmat[6][8]   = tmat[0][2];
  tmat[7][6]  = tmat[1][0];  tmat[7][7]   = tmat[1][1];  tmat[7][8]   = tmat[1][2];
  tmat[8][6]  = tmat[2][0];  tmat[8][7]   = tmat[2][1];  tmat[8][8]   = tmat[2][2];

  tmat[9][9]  = tmat[0][0];  tmat[9][10]  = tmat[0][1];  tmat[9][11]  = tmat[0][2];
  tmat[10][9] = tmat[1][0];  tmat[10][10] = tmat[1][1];  tmat[10][11] = tmat[1][2];
  tmat[11][9] = tmat[2][0];  tmat[11][10] = tmat[2][1];  tmat[11][11] = tmat[2][2];
}




/**
   function computes stiffness %matrix of 3D beam element
   
   @param eid - number of element
   @param sm - stiffness %matrix
   
   JK, 3.2.2002
*/
void beamel3d::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
  
}



/**
   The function computes stiffness %matrix of 3D beam element expressed
   in the local element coordinate system.
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   
   @return The function returns resulting %matrix in the parameter sm.

   PF, 20.12.2002
   Jan Zitny, 2012
*/
void beamel3d::stiffness_matrix_local (long eid,long ri,long ci,matrix &sm)
{
  long ii,i,i1;
  double e,g,a,ixyz[3],beta[2],dl,ll,gy,gz,aj1;
  ivector nodes(ASTCKIVEC(nne));
  vector vec(ASTCKVEC(3)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp)),tran(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid,vec);
  beam_transf_matrix (tran,dl,vec,x,y,z,eid);
  ll=dl*dl;
  
  fillm (0.0,sm);
  ii=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ii);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  e=d[0][0];  g=d[1][1];
  
  //  stiffness matrix in local coordinate system
  //  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
  //  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
  // axial forces and torsion
  sm[0][0]=  e*a/dl;
  sm[0][6]= -sm[0][0];
  sm[3][3]=  g*ixyz[0]/dl;
  sm[3][9]= -sm[3][3];
  //  direct s Iy
  aj1=ixyz[1]*e/(1.+2*gy);
  sm[2][2] = 12.*aj1/ll/dl;
  sm[2][4] =-6*aj1/ll;
  sm[2][8] =-sm[2][2];
  sm[2][10]= sm[2][4];
  sm[4][4] = 2.*(2.+gy)*aj1/dl;
  sm[4][8] = 6*aj1/ll;
  sm[4][10]= 2.*(1.-gy)*aj1/dl;
  //  direct s Iz
  aj1=ixyz[2]*e/(1.+2.*gz);
  sm[1][1] = 12.*aj1/ll/dl;
  sm[1][5] = 6.*aj1/ll;
  sm[1][7] =-sm[1][1];
  sm[1][11]= sm[1][5];
  sm[5][5] = 2.*(2.+gz)*aj1/dl;
  sm[5][7] =-6.*aj1/ll;
  sm[5][11]= 2.*(1.-gz)*aj1/dl;
  
  for (i=0;i<6;i++){
    i1=i+6;
    sm[i1][i1]=sm[i][i];
    sm[i1][i]=sm[i][i1];
  }
  
  sm[4][2] = sm[2][4];
  sm[5][1] = sm[1][5];
  sm[7][5] = sm[5][7];
  sm[8][4] = sm[4][8];
  sm[11][1]= sm[1][11];
  sm[10][2]= sm[2][10];
  sm[7][11]=-sm[1][5];
  sm[11][7]=-sm[1][5];
  sm[8][10]=-sm[2][4];
  sm[10][8]=-sm[2][4];
  
}

void beamel3d::stiffness_matrix_local2 (long eid,long ri,long ci,matrix &sm)
{
  long i,j,ipp;
  double e,g,a,ixyz[3],beta[2],l;
  ivector nodes(ASTCKIVEC(nne));
  vector vec(ASTCKVEC(3)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp)),tran(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid,vec);
  beam_transf_matrix (tran,l,vec,x,y,z,eid);
  
  fillm (0.0,sm);
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  e=d[0][0];  g=d[1][1];
  
  
  sm[0][0]=e*a/l;
  sm[0][5]=-1.0*e*a/l;

  sm[1][1]=12.0*e*ixyz[2]/l/l/l;
  sm[1][5]=6.0*e*ixyz[2]/l/l;
  sm[1][7]=-12.0*e*ixyz[2]/l/l/l;
  //sm[1][9]=6.0*e*ixyz[2]/l/l;
  sm[1][11]=6.0*e*ixyz[2]/l/l;
  
  sm[2][2]=12.0*e*ixyz[1]/l/l/l;
  sm[2][4]=-6.0*e*ixyz[1]/l/l;
  sm[2][8]=-12.0*e*ixyz[1]/l/l/l;
  sm[2][10]=-6.0*e*ixyz[1]/l/l;
  
  sm[3][3]=g*ixyz[0]/l;
  sm[3][9]=-1.0*g*ixyz[0]/l;
  
  sm[4][4]=4.0*e*ixyz[1]/l;
  sm[4][8]=6.0*e*ixyz[1]/l/l;
  sm[4][10]=2.0*e*ixyz[1]/l;

  sm[5][5]=4.0*e*ixyz[2]/l;
  sm[5][7]=-6.0*e*ixyz[2]/l/l;
  sm[5][11]=2.0*e*ixyz[2]/l;

  sm[6][6]=e*a/l;

  sm[7][7]=12.0*e*ixyz[2]/l/l/l;
  sm[7][11]=-6.0*e*ixyz[2]/l/l;

  sm[8][8]=12.0*e*ixyz[1]/l/l/l;
  sm[8][10]=6.0*e*ixyz[1]/l/l;
  
  sm[9][9]=g*ixyz[0]/l;
  
  sm[10][10]=4.0*e*ixyz[1]/l;

  sm[11][11]=4.0*e*ixyz[2]/l;
  
  for (i=0;i<12;i++){
    for (j=0;j<i;j++){
      sm[i][j]=sm[j][i];
    }
  }
  
}

/**
   function assembles local stiffness %matrix of 3D beam
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param sm - stiffness %matrix
   
   JK, 20. 7. 2019
*/
void beamel3d::stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm)
{
  long ipp;
  double e,g,a,ix,iy,iz,ixyz[3],l;
  vector vec(3),x(nne),y(nne),z(nne);
  matrix d(tncomp,tncomp),tran(ndofe,ndofe);
  
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid,vec);
  beam_transf_matrix (tran,l,vec,x,y,z,eid);
  
  fillm (0.0,sm);
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  
  e=d[0][0];  g=d[1][1];
  ix = ixyz[0];
  iy = ixyz[1];
  iz = ixyz[2];
  
  
  sm[0][0] =      e*a/l;
  sm[0][5] = -1.0*e*a/l;

  sm[1][1]  =  12.0*e*iz/l/l/l;
  sm[1][5]  =   6.0*e*iz/l/l;
  sm[1][7]  = -12.0*e*iz/l/l/l;
  sm[1][11] =   6.0*e*iz/l/l;
  
  sm[2][2]  =  12.0*e*iy/l/l/l;
  sm[2][4]  =  -6.0*e*iy/l/l;
  sm[2][8]  = -12.0*e*iy/l/l/l;
  sm[2][10] =  -6.0*e*iy/l/l;
  
  sm[3][3] =      g*ix/l;
  sm[3][9] = -1.0*g*ix/l;

  sm[4][2]  = -6.0*e*iy/l/l;
  sm[4][4]  =  4.0*e*iy/l;
  sm[4][8]  =  6.0*e*iy/l/l;
  sm[4][10] =  2.0*e*iy/l;
  
  sm[5][1]  =  6.0*e*iz/l/l;
  sm[5][5]  =  4.0*e*iz/l;
  sm[5][7]  = -6.0*e*iz/l/l;
  sm[5][11] =  2.0*e*iz/l;
  
  sm[6][0] = -1.0*e*a/l;
  sm[6][6] =      e*a/l;
  
  sm[7][1]  = -12.0*e*iz/l/l/l;
  sm[7][5]  =  -6.0*e*iz/l/l;
  sm[7][7]  =  12.0*e*iz/l/l/l;
  sm[7][11] =  -6.0*e*iz/l/l;
  
  sm[8][2]  = -12.0*e*iy/l/l/l;
  sm[8][4]  =   6.0*e*iy/l/l;
  sm[8][8]  =  12.0*e*iy/l/l/l;
  sm[8][10] =   6.0*e*iy/l/l;
  
  sm[9][3] = -1.0*g*ix/l;
  sm[9][9] =      g*ix/l;
  
  sm[10][2]  = -6.0*e*iy/l/l;
  sm[10][4]  =  2.0*e*iy/l;
  sm[10][8]  =  6.0*e*iy/l/l;
  sm[10][10] =  4.0*e*iy/l;
  
  sm[11][1]  =  6.0*e*iz/l/l;
  sm[11][5]  =  2.0*e*iz/l;
  sm[11][7]  = -6.0*e*iz/l/l;
  sm[11][11] =  4.0*e*iz/l;
  
}


/**
   The function computes stiffness %matrix of 3D beam element expressed
   in the global coordinate system or local nodal coordinate systems.
   Additionally, the condensation of the stiffness %matrix is performed.
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   @return The function returns resulting %matrix in the parameter sm.
   
   PF, 20.12.2002
   Jan Zitny 2012
*/
void beamel3d::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long ii,transf;
  double a,ixyz[3],beta[2],dl;
  ivector nodes(ASTCKIVEC(nne));
  vector vec(ASTCKVEC(3)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp)),tran(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid,vec);
  beam_transf_matrix (tran,dl,vec,x,y,z,eid);
  
  fillm (0.0,sm);
  ii=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ii);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  
  //stiffness_matrix_local2 (eid,ri,ci,sm);
  stiffness_matrix_local(eid,ri,ci,sm);
  
  
  //  transformation of stiffness matrix from local element coordinate
  //  system to the global coordinate system
  lgmatrixtransf (sm,tran);
  if(Gtm->gelements[eid].cne==2){
    ivector cu(ASTCKIVEC(ndofe));
  Gtm->give_cn(eid, cu.a);
  condense_matrix(sm, cu);
  }
  
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tran);
    glmatrixtransf (sm,tran);
  }
}



/**
  The function calculates member forces due to continuous load of element.

  @param eid - element id
  @param le  - 
  @param nv  - nodal values of continuous load
  @param nf  - resulting %vector of member forces

  @return The function returns resulting member forces in the parameter nf.

  Created by Jan Zitny 10.2012
*/
void beamel3d::nodeforces(long eid, long *le, double*nv, vector &nf)
{
  if(le[0]==1){
    double fxa=nv[0];
    double fya=nv[1];  
    double fza=nv[2];
    double mxa=nv[3];
    double mya=nv[4];
    double mza=nv[5];  
    double fxb=nv[6];
    double fyb=nv[7];  
    double fzb=nv[8];
    double mxb=nv[9];
    double myb=nv[10];
    double mzb=nv[11];

    double l;
    vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
    Mt->give_node_coord3d (x,y,z,eid);
    //  length of the element
    l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0])+(z[1]-z[0])*(z[1]-z[0]));
    if (l<Mp->zero)
      print_err("zero length of the %ld beamel3d element",__FILE__, __LINE__, __func__, eid);

    nf[0]= 1*(((fxa*l)/2) + (((fxb-fxa)*l)/6));
    nf[6]= 1*(((fxa*l)/2) + (((fxb-fxa)*l)/3));
    nf[1]= -1*(((fya*l)/2) + (((fyb-fya)*3*l)/20)) -1*(((18*mza)/l) + (((mzb-mza)*9)/l));
    nf[7]= -1*(((fya*l)/2) + (((fyb-fya)*7*l)/20)) + ((18*mza)/l) + (((mzb-mza)*9)/l);
    nf[2]= -1*(((fza*l)/2) + (((fzb-fza)*3*l)/20)) -1*(((18*mya)/l) + (((myb-mya)*9)/l));
    nf[8]= -1*(((fza*l)/2) + (((fzb-fza)*7*l)/20)) + ((18*mya)/l) + (((myb-mya)*9)/l);
    nf[3]= (((2*mxa+mxb)*l)/6);
    nf[9]= (((mxa+2*mxb)*l)/6);
    nf[4]= ((fza*l*l)/12)+(((fzb-fza)*l*l)/30) + 6*mya +3*(myb-mya);
    nf[10]= -1*(((fza*l*l)/12)+(((fzb-fza)*l*l)/20)) + 12*mya +6*(myb-mya);
    nf[5]= ((fya*l*l)/12)+(((fyb-fya)*l*l)/30) + 6*mza +3*(mzb-mza);
    nf[11]= -1*(((fya*l*l)/12)+(((fyb-fya)*l*l)/20)) + 12*mza +6*(mzb-mza);
  }    

  if(Gtm->gelements[eid].cne==2){
      
    matrix sm(ASTCKMAT(ndofe,ndofe));
    stiffness_matrix_local (eid,0,0,sm);
    
    ivector cu(ASTCKIVEC(ndofe));
    Gtm->give_cn(eid, cu.a);
    condense_vector(sm,nf,cu);
  }
}



/**
   function computes consistent mass %matrix of 3D beam element
   rotational DOFs are neglected
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param mm - mass %matrix of one element
   
   JK, 3.2.200
*/
void beamel3d::mass_matrix (long eid,long /*ri*/,long /*ci*/,matrix &mm)
{
  long i,transf;
  double a,l,gy,gz,xi,jac,rho,beta[2];
  ivector nodes(ASTCKIVEC(nne));
  vector vec(ASTCKVEC(3)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix n(ASTCKMAT(3,ndofe)),tran(ASTCKMAT(12,12));
  
  fillm (0.0,mm);
  
  // element nodes
  Mt->give_elemnodes (eid,nodes);
  //  area of cross section
  Mc->give_area (eid,a);
  //  shear coefficients
  Mc->give_shearcoeff (eid,beta);
  //  density of the material
  Mc->give_densitye (eid,rho);
  //  vector defining local coordinate system
  Mc->give_vectorlcs (eid, vec);

  gy=beta[0];
  gz=beta[1];

  Mt->give_node_coord3d (x,y,z,eid);

  beam_transf_matrix (tran,l,vec,x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordmm);
  
  
  for (i=0;i<intordmm;i++){
    //  xi has to be between 0.0 and 1.0 while gp[i] are between -1.0 and 1.0
    xi=0.5*(1.0+gp[i]);
    //  matrix of approximation functions
    bf_matrix_transl (n,xi,l,gy,gz);
    jac=w[i]*l/2.0*rho*a;
    
    //  N^T rho N
    nnj (mm.a,n.a,jac,n.m,n.n);
    
  }
  
  //  transformation of mass matrix from local element coordinate
  //  system to the global coordinate system
  lgmatrixtransf (mm,tran);

  //  transformation of mass matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tran);
    glmatrixtransf (mm,tran);
  }

}



/**
   function computes consistent mass %matrix of 3D beam element
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK, 2.9.2006
*/
void beamel3d::res_mass_matrix (long eid,matrix &mm)
{
  mass_matrix (eid,0,0,mm);
}



/**
   @param eid - element id
   @param ri,ci - row and column indices
   @param ism - initial stress %matrix
   
*/
void beamel3d::initstr_matrix (long eid,long ri,long ci,matrix &ism)
{
  long i,i1,ii,transf;
  double nn,  e,g,a,dl,ll,ixyz[3],beta[2],gy,gz,g2,a4;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),vec(ASTCKVEC(3));
  matrix d(ASTCKMAT(tncomp,tncomp)),tran(ASTCKMAT(12,12));

// nn....... Axial forces

  // !!!!!!!!!
  nn = 0.0; // temporary value, it has to be set to correct value
  // !!!!!!!!!

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);

//  Mc->give_vecbeam (eid,&vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  ll=dl*dl;
  
  fillm (0.0,ism);
  ii=Mt->elements[eid].ipp[ri][ci];
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  Mm->matstiff (d,ii);
  
  e=d[0][0];  g=d[1][1];


//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
  
  g2=gy*gy;
  a4=nn/(1.+2*gy)/(1.+2*gy);
  ism[2][2] = a4*4*(g2+gy+0.3)/dl;
  ism[4][4] = a4*(g2+gy+0.4)/3.*dl;
  ism[2][4] = - a4/10.;
  ism[2][8] = - a4*4*(g2+gy+0.3)/dl;
  ism[4][10]= a4*(-g2-gy-0.1)/3.*dl;
  ism[4][8] = a4/10.;
  ism[2][10]=- a4/10.;
//  direct s Iz
  g2=gz*gz;
  a4=nn/(1.+2.*gz)/(1.+2.*gz);
  ism[1][1] = a4*4.*(g2+gz+0.3)/dl;
  ism[5][5] = a4*(g2+gz+0.4)/3.*dl;
  ism[1][5] = a4/10.;
  ism[1][7] =- a4*4.*(g2+gz+0.3)/dl;
  ism[5][11]= a4*(-g2-gz-0.1)/3.*dl;
  ism[5][7] =- a4/10.;
  ism[1][11]= a4/10.;

  for (i=0;i<6;i++){
	  i1=i+6;
      ism[i1][i1]=ism[i][i];
      ism[i1][i]=ism[i][i1];
  }

  ism[4][2] = ism[2][4];
  ism[5][1] = ism[1][5];
  ism[7][5] = ism[5][7];
  ism[8][4] = ism[4][8];
  ism[11][1]= ism[1][11];
  ism[10][2]= ism[2][10];
  ism[7][11]=-ism[1][5];
  ism[11][7]=-ism[1][5];
  ism[8][10]=-ism[2][4];
  ism[10][8]=-ism[2][4];

  //  transformation of initial stress matrix from local element coordinate
  //  system to the global coordinate system
  lgmatrixtransf (ism,tran);

  //  transformation of initial stress matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tran);
    glmatrixtransf (ism,tran);
  }
}



/**
   function serves for computation of nodal forces and moments caused by body load
   
   @param eid - element id
   @param lm - load %matrix
   
   JK, 1.9.2006
*/
void beamel3d::load_matrix (long eid,matrix &lm)
{
  long i,ii,transf;
  double jac,a,ixyz[3],beta[2],e,g,gy,gz,l;
  ivector nodes (ASTCKIVEC(nne));
  vector vec(ASTCKVEC(3)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),xl(ASTCKVEC(2));
  matrix d(ASTCKMAT(tncomp,tncomp)),n(ASTCKMAT(napfun,ndofe)),tran(ASTCKMAT(12,12));
  
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
  lgmatrixtransf (lm,tran);

  //  transformation of stiffness matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    //  assembling of transformation matrix from global to nodal system
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }
}



/**
   The function stores nodal displacements and rotations in the
   local element coordinate system. This is equivalent to functions computing strains.
   Nodal displacements and rotations are better quantities in case of beams.
   
   @param eid - element id
   @param lcid - load case id
   
   JK, 20.2.2002
*/
void beamel3d::nodal_displ (long lcid,long eid)
{
  long i,j,ip,transf;
  double l;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),vec(ASTCKVEC(3)),d(ASTCKVEC(ndofe));
  vector dd(ASTCKVEC(6)),f(ASTCKVEC(ndofe));
  matrix tm(ASTCKMAT(3,3)),tmat (ASTCKMAT(ndofe,ndofe));
  
  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  element node numbers
  Mt->give_elemnodes (eid,nodes);

  //  transformation from nodal to global coordinate system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    //  assembling of transformation matrix
    transf_matrix (nodes,tmat);
    //  transformation of nodal displacements and rotations
    lgvectortransf (f,d,tmat);
    copyv (f,d);
  }
  
  //  nodal coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  vector of z-axis
  Mc->give_vectorlcs (eid,vec);
  
  //  assembling of transformation matrix from global to local coordinate system
  beam_transf_matrix (tmat,l,vec,x,y,z,eid);
  
  //  transformation from global problem coordinate system
  //  to local element coordinate system
  glvectortransf (d,f,tmat);
  copyv (f,d);
  
  j=6;
  for (i=0;i<6;i++){
    dd[i]=d[j];
    j++;
  }

  //  number of integration point
  ip=Mt->elements[eid].ipp[0][0];
  
  
  Mm->storestrain (lcid,ip,0,6,d);
  Mm->storestrain (lcid,ip+1,0,6,dd);
}



/**
   The function computes nodal forces and moments expressed in local 
   coodinate system. This is equivalent to functions computing stresses on 
   plane or space elements. Nodal forces and moments are better quantities 
   in case of beams.

   x-axis goes through the beam
   z-axis is defined in cross section
   y-axis is defined by right hand orientation
   
   @param eid - element id
   @param lcid - load case id
   
   JK, 20.2.2002, revised 2.9.2006
*/
void beamel3d::nodal_forces (long lcid,long eid)
{
  long i,j,ii,transf;
  double l;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector vec(ASTCKVEC(3)),d(ASTCKVEC(ndofe)),f(ASTCKVEC(ndofe)),ff(ASTCKVEC(6));
  matrix sm(ASTCKMAT(ndofe,ndofe)),tmat(ASTCKMAT(ndofe,ndofe));
  
  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  assembling of stiffness matrix
  stiffness_matrix (eid,0,0,sm);
  
  //  computation of nodal forces and moments
  mxv (sm,d,f);


  //  transformation of displacement vector from local nodal systems
  //  to global coordinate system
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    //  assembling of transformation matrix
    transf_matrix (nodes,tmat);
    //  transformation of nodal displacements and rotations
    lgvectortransf (d,f,tmat);
    copyv (d,f);
  }


  //  nodal coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  vector of z-axis
  Mc->give_vectorlcs (eid,vec);
  
  //  assembling of transformation matrix from global to local coordinate system
  beam_transf_matrix (tmat,l,vec,x,y,z,eid);
  
  //  transformation from global problem coordinate system
  //  to local element coordinate system
  glvectortransf (f,d,tmat);
  copyv (d,f);
  
  j=6;
  for (i=0;i<6;i++){
    ff[i]=f[j];
    j++;
  }
  
  //  number of integration point
  ii=Mt->elements[eid].ipp[0][0];
  
  //  storage of nodal forces and moments
  //  nodal forces and moments are expressed in local coordinate system
  Mm->storestress (lcid,ii,0,6,f);
  Mm->storestress (lcid,ii+1,0,6,ff);
}



/**
   The function computes internal forces, i.e. member forces of the element
   expressed in the global coordinate system or local nodal coordinate systems.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   20.12.2002
*/
void beamel3d::res_internal_forces (long lcid,long eid,vector &ifor)
{
//  internal_forces (lcid, eid, 0, 0, ifor);
  matrix sm(ASTCKMAT(ndofe, ndofe));
  vector r(ASTCKVEC(ndofe));

  stiffness_matrix (eid, 0, 0, sm);
  eldispl (lcid, eid, r.a);  

  mxv(sm, r, ifor);
}




void beamel3d::internal_forces (long lcid,long eid,long /*ri*/,long /*ci*/,vector &ifor)
{
  long i,integr=2,ipp,k,transf;
  double dl,ll,a,ixyz[3],beta[2],gy,gz,e,g,s,dj,   m;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),vec(ASTCKVEC(3));
  vector rl(ASTCKVEC(ndofe)),rg(ASTCKVEC(ndofe)),f(ASTCKVEC(ndofe));
  vector fx(ASTCKVEC(tncomp)),w(ASTCKVEC(integr)),gp(ASTCKVEC(integr));
  matrix n(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp)),tran(ASTCKMAT(12,12)),r(ASTCKMAT(tncomp,2));
  
  // r  .... end forces last time
  // m  .... dl(i)/dl(i-1)
  Mt->give_elemnodes(eid,nodes);
  Mt->give_node_coord3d(x,y,z,eid);
  eldispl (lcid,eid,rl.a);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran, dl, vec, x, y, z, eid);
  ll=dl*dl;

  // r  .... end forces last time
  // m  .... dl(i)/dl(i-1)
  m=dl;
  fillm (0.0,r);
  
  //  transformation of displacement vector
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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
  
  e=d[0][0];  g=d[1][1];
  
  //  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
  //  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
  gauss_points (gp.a,w.a,integr);
  for (i=0;i<integr;i++){
    //	 s=x/dl
    s=(gp[i]+1.)/2.; // gaus point in (-1,1) function on beam (0,1)
    dj=w[i]/2.*dl;   // weight divide by 2
    geom_matrix (n,s, dl, gy, gz);
    fx[0]=( n[0][0]*rl[0]+n[0][6]*rl[6] )*a*e;
    fx[3]=( n[3][3]*rl[3]+n[3][9]*rl[9] )*ixyz[0]*g;
    fx[1]=( n[1][1]*rl[1]+n[1][5]*rl[5]+n[1][7]*rl[7]+n[1][11]*rl[11] )*beta[1]*a*g;
    fx[5]=( n[5][1]*rl[1]+n[5][5]*rl[5]+n[5][7]*rl[7]+n[5][11]*rl[11] )*ixyz[2]*e;
    fx[2]=( n[2][2]*rl[2]+n[2][4]*rl[4]+n[2][8]*rl[8]+n[2][10]*rl[10] )*beta[0]*a*g;
    fx[4]=( n[4][2]*rl[2]+n[4][4]*rl[4]+n[4][8]*rl[8]+n[4][10]*rl[10] )*ixyz[1]*e;
    
    //  change of the II piola-Kirchhoff to Cauchyho real stress (through scale m)
    f[0]=f[0]+( r[0][0]*(1-s) + r[0][1]*s+fx[0] )*n[0][0]*dj;
    f[6]=f[6]+( r[0][0]*(1-s) + r[0][1]*s-fx[0] )*n[0][6]*dj;
    f[3]=f[3]+( r[3][0]*(1-s) + r[3][1]*s+fx[3] )*n[3][3]*dj;
    f[9]=f[9]+( r[3][0]*(1-s) + r[3][1]*s-fx[3] )*n[3][9]*dj;
    // Vy, Mz
    //     f[1]=f[1]+( (r[1][0]*(1-s)+r[1][1]*s)*n[1][1] + ((r[5][0]*(1-s)+r[5][1]*s)*dl/m)*n[1][5] )*dj;
    f[1]=f[1]  +(fx[1]*n[1][1] + fx[5]*n[1][5])*dj;
    
    //     f[7]=f[7]+( (r[1][0]*(1-s)+r[1][1]*s)*n[1][1] + ((r[5][0]*(1-s)+r[5][1]*s)*dl/m)*n[5][7] )*dj;
    f[7]=f[7]  +(fx[1]*n[1][7] + fx[5]*n[5][7])*dj;
    
    //     f[5]=f[5]+( ((r[5][0]*(1-s)+r[5][1]*s)*dl/m)*n[5][5] + (r[1][0]*(1-s)+r[1][1]*s)*n[1][5] )*dj;
    f[5]=f[5]  +(fx[5]*n[5][5] + fx[1]*n[1][5])*dj;
    
    //     f[11]=f[11]+( ((r[5][0]*(1-s)+r[5][1]*s)*dl/m)*n[11][5] + (r[1][0]*(1-s)+r[1][1]*s)*n[11][1] )*dj;
    f[11]=f[11]+(fx[5]*n[5][11] + fx[1]*n[1][11])*dj;
    // Vz, My
    //     f[2]=f[2]+( (r[2][0]*(1-s)+r[2][1]*s)*n[2][2] + ((r[4][0]*(1-s)+r[4][1]*s)*dl/m)*n[4][2] )*dj;
    f[2]=f[2]  +(fx[2]*n[2][2] + fx[4]*n[4][2])*dj;
    
    //     f[8]=f[8]+( (r[2][0]*(1-s)+r[2][1]*s)*n[2][8] + ((r[4][0]*(1-s)+r[4][1]*s)*dl/m)*n[4][8] )*dj;
    f[8]=f[8]  +(fx[2]*n[2][8] + fx[4]*n[4][8])*dj;
    
    //     f[4]=f[4]+( ((r[4][0]*(1-s)+r[4][1]*s)*dl/m)*n[4][4] + (r[2][0]*(1-s)+r[2][1]*s)*n[2][4] )*dj;
    f[4]=f[4]  +(fx[4]*n[4][4] + fx[2]*n[2][4] )*dj;
    
    //     f[10]=f[10]+( ((r[4][0]*(1-s)+r[4][1]*s)*dl/m)*n[4][10] + (r[2][0]*(1-s)+r[2][1]*s)*n[2][10] )*dj;
    f[10]=f[10]+(fx[4]*n[4][10] + fx[2]*n[2][10])*dj;
  }
  
  for (k=0;k<ndofe;k++){
    ifor[k]+=f[k];
    ifor[k]=0.0;                     //nekonsoliduje 
  }
  //  fprintf (Out,"\n\n R prvek beam cislo %ld\n",eid);
  //  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",rl[k]); }
  //  fprintf (Out,"\n\n Fint prvek cislo %ld\n",eid);
  //  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",f[k]);}
}



/**
   The function defines meaning of DOFs at nodes.
   
   @param eid - element id
   
   1.2.2005, JK
*/
void beamel3d::define_meaning (long eid)
{
  ivector cn(ASTCKIVEC(ndofe)),nod(ASTCKIVEC(nne));
  
  Mt->give_elemnodes (eid,nod);
  Mt->give_code_numbers (eid,cn.a);

  //  displacement in x direction
  if (cn[0]>0)  Mt->nodes[nod[0]].meaning[0]=1;
  //  displacement in y direction
  if (cn[1]>0)  Mt->nodes[nod[0]].meaning[1]=2;
  //  displacement in z direction
  if (cn[2]>0)  Mt->nodes[nod[0]].meaning[1]=3;
  //  displacement in x direction
  if (cn[6]>0)  Mt->nodes[nod[1]].meaning[0]=1;
  //  displacement in y direction
  if (cn[7]>0)  Mt->nodes[nod[1]].meaning[1]=2;
  //  displacement in z direction
  if (cn[8]>0)  Mt->nodes[nod[1]].meaning[1]=3;
}


