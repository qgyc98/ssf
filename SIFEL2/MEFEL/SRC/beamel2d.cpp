#include "beamel2d.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "loadcase.h"
#include "loadel.h"
#include <math.h>


beamel2d::beamel2d (void)
{
  long i;
  
  //  number nodes on element
  nne=2;
  //  number of DOFs on element
  ndofe=6;
  //  number of strain/stress components
  tncomp=3;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=4;
  //  order of numerical integration of initial stiffness matrix
  intordism=3;
  ssst=plbeam;

  //  coordinate system
  //  1 - xz, default version
  //  2 - xy, it is used for coupling with plane elements
  coordsys=2;
  
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
  nip[0][0]=3;
  intordsm[0][0]=3;

  //  total number of integration points
  tnip=3;
}

beamel2d::~beamel2d (void)
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
   Function approximates function defined by end node values.

   @param xi - natural coordinate
   @param nodval - nodal values
   
   @retval Function returns approximated value for the given natural coordinate xi.

   JK
*/
double beamel2d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));

  bf_lin_1d (bf.a,xi);

  scprd (bf,nodval,f);

  return f;
}



/**
   function assembles %matrix of approximation functions

   ordering of approximation functions
   u - displacement along x
   w - deflection - displacement along z
   phi - rotation around y

   @param n - array containing %matrix of approximation functions
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   JK, 20.2.2002
*/
void beamel2d::bf_matrix (matrix &n,double xi,double l,double kappa)
{
  vector u(2),w(4),r(4);

  dilat (u.a,xi);
  defl2D_fun (w.a,xi,l,kappa);
  roty_fun (r.a,xi,l,kappa);

  fillm (0.0,n);
  
  n[0][0]=u[0];
  n[0][3]=u[1];
  
  n[1][1]=w[0];
  n[1][2]=w[1];
  n[1][4]=w[2];
  n[1][5]=w[3];
  
  n[2][1]=r[0];
  n[2][2]=r[1];
  n[2][4]=r[2];
  n[2][5]=r[3];
}

/**
   function assembles %matrix of derivatives of approximation function of deflection
   
   @param n - %matrix of derivatives
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - shear coefficient 6.0*E*I/k/G/A/l/l
   
   JK
*/
void beamel2d::dbf_matrix (matrix &n,double xi,double l,double kappa)
{
  vector dw(4);
  
  der_defl2D_fun (dw.a,xi,l,kappa);
  
  fillm (0.0,n);
  
  n[0][1]=dw[0];
  n[0][2]=dw[1];
  n[0][4]=dw[2];
  n[0][5]=dw[3];
  
}

/**
   function assembles strain-displacement (geometric) %matrix
   
   ordering of strain components
   du/dx - strain e_xx
   phi+dw/dx - strain e_xz
   dphi/dx - curvature
   
   @param gm - geometric %matrix
   @param xi - natural coordinate from segment <0;1>
   @param l - length of the beam
   @param kappa - 6.0*E*I/k/G/A/l/l
   
   JK, 20.2.2002
*/
void beamel2d::geom_matrix (matrix &gm,double xi,double l,double kappa)
{
  vector du(2),dw(4),r(4),dr(4);
  
  der_dilat (du.a,l);
  der_defl2D_fun (dw.a,xi,l,kappa);
  roty_fun (r.a,xi,l,kappa);
  der_roty_fun (dr.a,xi,l,kappa);
  
  fillm (0.0,gm);
  
  gm[0][0] = du[0];
  gm[0][3] = du[1];
  
  gm[1][1] = r[0]+dw[0];
  gm[1][2] = r[1]+dw[1];
  gm[1][4] = r[2]+dw[2];
  gm[1][5] = r[3]+dw[3];
  
  gm[2][1] = dr[0];
  gm[2][2] = dr[1];
  gm[2][4] = dr[2];
  gm[2][5] = dr[3];
}





/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   JK, 3.2.2002
*/
void beamel2d::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*3+0][i*3] = Mt->nodes[nodes[i]].e1[0];    tmat[i*3+0][i*3+1] = Mt->nodes[nodes[i]].e2[0];    tmat[i*3+0][i*3+2] = 0.0;
      tmat[i*3+1][i*3] = Mt->nodes[nodes[i]].e1[1];    tmat[i*3+1][i*3+1] = Mt->nodes[nodes[i]].e2[1];    tmat[i*3+1][i*3+2] = 0.0;
      tmat[i*3+2][i*3] = 0.0;                          tmat[i*3+2][i*3+1] = 0.0;                          tmat[i*3+2][i*3+2] = 1.0;
      
    }
  }
}


/**
   function assembles transformation %matrix from local element coordinate
   system to the global problem coordinate system x_g = T x_l
   
   the transformation %matrix is the same for the coordinate system xy and xz
   therefore, only the notation xz is used
   if the system xy is used, one has to send into the function the y-coordinates
   stored in the %vector z
   
   xz coordinate system
   the beam element is formulated in x-z plane, x axis is identical with beam axis
   and it is oriented horizontally, z axis is oriented vertically (downwards)

   xy coordinate system
   the beam element is formulated in x-y plane, x axis is identical with beam axis
   and it is oriented horizontally, y axis is oriented vertically (upwards)
   
   @param x,z - node coordinates, the %vector z contains z coordinates in the case xz
                and y coordinates in the case xy
   @param tmat - transformation %matrix
   
   JK, 3.2.2002
*/
void beamel2d::beam_transf_matrix (long eid,vector &x,vector &z,matrix &tmat)
{
  double c,s,l;
  
  fillm (0.0,tmat);
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__,__LINE__,__func__,eid);
  }
  
  //  sine and cosine of element angle
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;
  
  //  transformation matrix
  tmat[0][0] = c;    tmat[0][1] = -1.0*s;    tmat[0][2] = 0.0;
  tmat[1][0] = s;    tmat[1][1] = c;         tmat[1][2] = 0.0;
  tmat[2][0] = 0.0;  tmat[2][1] = 0.0;       tmat[2][2] = 1.0;

  tmat[3][3] = c;    tmat[3][4] = -1.0*s;    tmat[3][5] = 0.0;
  tmat[4][3] = s;    tmat[4][4] = c;         tmat[4][5] = 0.0;
  tmat[5][3] = 0.0;  tmat[5][4] = 0.0;       tmat[5][5] = 1.0;
}



/**
   function computes stiffness %matrix of two-dimensional beam finite element
   element is based on Mindlin theory (shear stresses are taken into account)
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   JK, 3.2.2002
*/
void beamel2d::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,ipid,transf;
  double e,g,a,iy,shearcoeff,kappa,l,xi,jac,ww;
  ivector nodes(nne);
  vector x(nne),z(nne),w(intordsm[0][0]),gp(intordsm[0][0]);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp),tmat(ndofe,ndofe);
  
  //  node coordinates in the x-z plane
  Mt->give_node_coord2dxz (x,z,eid);
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__,__LINE__,__func__,eid);
  }
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  
  fillm (0.0,sm);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  area of element cross section
  Mc->give_area (eid,a);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero){  kappa=0.0;  shearcoeff=1.0;        }
  else{                      kappa=6.0*e*iy/shearcoeff/g/a/l/l; }
  
  d[0][0]*=a;
  d[1][1]*=a*shearcoeff;
  d[2][2]*=iy;
  
  //  stiffness matrix in local coordinate system
  for (i=0;i<intordsm[0][0];i++){
    xi=(1.0+gp[i])/2.0;  ww=w[i];
    geom_matrix (gm,xi,l,kappa);
    jac=l/2.0*ww;
    bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
  }
  
  //  transformation of stiffness matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (sm,tmat);

  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
}

/**
   function computes stiffness %matrix of two-dimensional beam finite element
   element is based on Mindlin theory (shear stresses are taken into account)
   
   @param eid - number of element
   @param sm - stiffness %matrix

   JK, 3.2.2002
*/
void beamel2d::res_stiffness_matrix (long eid,matrix &sm)
{
  //stiffness_matrix (eid,0,0,sm);
  stiffness_matrix_expl (eid,0,0,sm);
}

/**
   function computes stiffness %matrix of two-dimensional beam finite element
   element is based on Mindlin theory (shear stresses are taken into account)
   stiffness %matrix is formulated explicitly
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   JK, 3.2.2002
*/
void beamel2d::stiffness_matrix_expl_local (long eid,long ri,long ci,double l,matrix &sm)
{
  long ipid;
  double e,g,shearcoeff,a,iy,kappa;
  ivector nodes(nne);
  matrix d(tncomp,tncomp),tmat (ndofe,ndofe);
  
  fillm (0.0,sm);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  area of element cross section
  Mc->give_area (eid,a);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  if (coordsys==1){
    //  x-z plane
    
    //  stiffness matrix in local coordinate system
    sm[0][0]=      e*a/l;
    sm[0][3]= -1.0*e*a/l;
    
    sm[1][1] =  12.0*e*iy/l/l/l;
    sm[1][2] =  -6.0*e*iy/l/l;
    sm[1][4] = -12.0*e*iy/l/l/l;
    sm[1][5] =  -6.0*e*iy/l/l;
    
    sm[2][1] = -6.0*e*iy/l/l;
    sm[2][2] =  4.0*e*iy/l;
    sm[2][4] =  6.0*e*iy/l/l;
    sm[2][5] =  2.0*e*iy/l;
    
    sm[3][0] = -1.0*e*a/l;
    sm[3][3] =      e*a/l;
    
    sm[4][1] = -12.0*e*iy/l/l/l;
    sm[4][2] =   6.0*e*iy/l/l;
    sm[4][4] =  12.0*e*iy/l/l/l;
    sm[4][5] =   6.0*e*iy/l/l;
    
    sm[5][1] = -6.0*e*iy/l/l;
    sm[5][2] =  2.0*e*iy/l;
    sm[5][4] =  6.0*e*iy/l/l;
    sm[5][5] =  4.0*e*iy/l;
  }

  if (coordsys==2){
    //  x-y plane
    
    //  stiffness matrix in local coordinate system
    sm[0][0]=      e*a/l;
    sm[0][3]= -1.0*e*a/l;
    
    sm[1][1] =  12.0*e*iy/l/l/l;
    sm[1][2] =   6.0*e*iy/l/l;
    sm[1][4] = -12.0*e*iy/l/l/l;
    sm[1][5] =   6.0*e*iy/l/l;
    
    sm[2][1] =  6.0*e*iy/l/l;
    sm[2][2] =  4.0*e*iy/l;
    sm[2][4] = -6.0*e*iy/l/l;
    sm[2][5] =  2.0*e*iy/l;
    
    sm[3][0] = -1.0*e*a/l;
    sm[3][3] =      e*a/l;
    
    sm[4][1] = -12.0*e*iy/l/l/l;
    sm[4][2] =  -6.0*e*iy/l/l;
    sm[4][4] =  12.0*e*iy/l/l/l;
    sm[4][5] =  -6.0*e*iy/l/l;
    
    sm[5][1] =  6.0*e*iy/l/l;
    sm[5][2] =  2.0*e*iy/l;
    sm[5][4] = -6.0*e*iy/l/l;
    sm[5][5] =  4.0*e*iy/l;
  }


}

  

void beamel2d::stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm)
{
  long ipid,transf;
  double e,g,shearcoeff,a,iy,kappa,l;
  ivector nodes(nne);
  vector x(nne),z(nne);
  matrix d(tncomp,tncomp),tmat (ndofe,ndofe);
  
  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element", __FILE__, __LINE__, __func__,eid);
  }
  
  fillm (0.0,sm);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  area of element cross section
  Mc->give_area (eid,a);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  stiffness_matrix_expl_local(eid,ri,ci,l,sm);
  
  
  if(Gtm->gelements[eid].cne==2){
    ivector cu(ndofe);
    Gtm->give_cn(eid, cu.a);
    condense_matrix(sm, cu);
  }
  
  
  //  transformation of stiffness matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (sm,tmat);
  
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
}



/**
   function computes consistent mass %matrix of 2D beam element
   influence of inertial forces from rotations is accounted
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param mm - mass %matrix
   
   JK, 3.2.2006
*/
void beamel2d::mass_matrix (long eid,long ri,long ci,matrix &mm)
{
  long i,ipid,transf;
  double xi,ww,jac,e,g,shearcoeff,iy,a,l,kappa,rho;
  ivector nodes(nne);
  vector x(nne),z(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe),d(tncomp,tncomp),tmat (ndofe,ndofe);


  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
  }

  gauss_points (gp.a,w.a,intordmm);

  fillm (0.0,mm);
  ipid=Mt->elements[eid].ipp[ri][ci];

  //  area of element cross section
  Mc->give_area (eid,a);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  density of the material
  Mc->give_densitye (eid,rho);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  fillm (0.0,d);
  d[0][0]=1.0;  d[1][1]=1.0;  d[2][2]=iy;
  
  for (i=0;i<intordmm;i++){
    xi=(1.0+gp[i])/2.0;  ww=w[i];
    bf_matrix (n,xi,l,kappa);
    jac=a*l/2.0*ww*rho;
    bdbj (mm.a,n.a,d.a,jac,n.m,n.n);
  }


  //  transformation of mass matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (mm,tmat);

  //  transformation of mass matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}

/**
   function computes consistent mass %matrix of 2D beam element
   influence of inertial forces from rotations is accounted
   
   @param eid - number of element
   @param mm - mass %matrix
   
   JK, 3.2.2006
*/
void beamel2d::res_mass_matrix (long eid,matrix &mm)
{
  //mass_matrix (eid,0,0,mm);
  mass_matrix_expl (eid,0,0,mm);
}

/**
   function computes consistent mass %matrix of 2D beam element
   influence of inertial forces from rotations is accounted
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param mm - mass %matrix
   
   JK, 3.2.2006
*/
void beamel2d::mass_matrix_expl (long eid,long ri,long ci,matrix &mm)
{
  long ipid,transf;
  double e,g,shearcoeff,j,iy,a,l,kappa,rho;
  ivector nodes(nne);
  vector x(nne),z(nne);
  matrix d(tncomp,tncomp),tmat(ndofe,ndofe);

  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
  }
  
  fillm (0.0,mm);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  area of element cross section
  Mc->give_area (eid,a);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  density of the material
  Mc->give_densitye (eid,rho);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  j = (1.0+2.0*kappa)*(1.0+2.0*kappa);
  
  
  mm[0][0]=rho*a*l/3.0;
  mm[0][3]=rho*a*l/6.0;
  
  mm[1][1]=rho*a*l*(13.0/35.0+7.0/5.0*kappa+4.0/3.0*kappa*kappa)/j;
  mm[1][2]=rho*a*l*l*(-11.0/210.0-11.0/60*kappa-1.0/6.0*kappa*kappa)/j;
  mm[1][4]=rho*a*l*(9.0/70+3.0/5.0*kappa+2.0/3.0*kappa*kappa)/j;
  mm[1][5]=rho*a*l*l*(13.0/420+3.0/20.0*kappa+1.0/6.0*kappa*kappa)/j;
  
  mm[2][1]=mm[1][2];
  mm[2][2]=rho*a*l*l*l*(1.0/105.0+1.0/30.0*kappa+1.0/30.0*kappa*kappa)/j;
  mm[2][4]=rho*a*l*l*(-13.0/420.0-3.0/20.0*kappa-1.0/6.0*kappa*kappa)/j;
  mm[2][5]=rho*a*l*l*l*(-1.0/140.0-1.0/30.0*kappa-1.0/30.0*kappa*kappa)/j;
  
  mm[3][0]=mm[0][3];
  mm[3][3]=mm[0][0];
  
  mm[4][1]=mm[1][4];
  mm[4][2]=mm[2][4];
  mm[4][4]=rho*a*l*(13.0/35.0+7.0/5.0*kappa+4.0/3.0*kappa*kappa)/j;
  mm[4][5]=rho*a*l*l*(11.0/210.0+11.0/60.0*kappa+1.0/6.0*kappa*kappa)/j;
  
  mm[5][1]=mm[1][5];
  mm[5][2]=mm[2][5];
  mm[5][4]=mm[4][5];
  mm[5][5]=rho*a*l*l*l*(1.0/105+1.0/30.0*kappa+1.0/30.0*kappa*kappa)/j;

  
  /*
  mm[0][0]=rho*a*l/2.0;
  mm[1][1]=rho*a*l/24.0*12.0;
  mm[2][2]=rho*a*l*l*l/24.0;
  mm[3][3]=rho*a*l/2.0;
  mm[4][4]=rho*a*l/24.0*12.0;
  mm[5][5]=rho*a*l*l*l/24.0;
  */
  

  //  transformation of mass matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (mm,tmat);

  //  transformation of mass matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
  
}


/**
   function computes consistent geometric-stiffness %matrix
   notation based on R. Clough, J. Penzien: Dynamic of Structures,
   McGraw-Hill, 2nd edition, 1993, pages 195-196
   
   notation initial stress %matrix is also used
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param ism - inital stress %matrix
   
   JK
*/
void beamel2d::initstr_matrix (long eid,long ri,long ci,matrix &ism)
{
  long i,ipid,transf;
  double xi,ww,jac,e,g,shearcoeff,iy,a,l,kappa;
  ivector nodes(nne);
  vector x(nne),z(nne),w(intordism),gp(intordism);
  matrix n(1,ndofe),d(tncomp,tncomp),tmat(ndofe,ndofe);

  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
  }

  gauss_points (gp.a,w.a,intordism);

  
  fillm (0.0,ism);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  area of element cross section
  Mc->give_area (eid,a);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];

  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  fillm (0.0,d);
  d[0][0]=1.0;  d[1][1]=1.0;  d[2][2]=iy;

  for (i=0;i<intordism;i++){
    xi=(1.0+gp[i])/2.0;  ww=w[i];
    dbf_matrix (n,xi,l,kappa);
    jac=a*l/2.0*ww;
    nnj (ism.a,n.a,jac,n.m,n.n);
  }

  //  transformation of initial stress matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (ism,tmat);

  //  transformation of initial stress matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (ism,tmat);
  }
}


/**
   function computes consistent geometric-stiffness %matrix
   notation based on R. Clough, J. Penzien: Dynamic of Structures,
   McGraw-Hill, 2nd edition, 1993, pages 195-196
   
   notation initial stress %matrix is also used
   
   initial stress %matrix is formulated explicitly
   
   @param lcid - load case id
   @param eid - number of element
   @param ri,ci - row and column indices
   @param ism - inital stress %matrix
   
   JK
*/
void beamel2d::initstr_matrix_expl (long lcid,long eid,long ri,long ci,matrix &ism)
{
  long ipid,transf;
  double e,g,shearcoeff,iy,a,l,kappa,denom,nforce;
  ivector nodes(nne);
  vector x(nne),z(nne),w(intordism),gp(intordism),f(3);
  matrix d(tncomp,tncomp),tmat (ndofe,ndofe);


  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
  }

  gauss_points (gp.a,w.a,intordism);

  fillm (0.0,ism);
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  area of element cross section
  Mc->give_area (eid,a);
  //  shear coefficient
  Mc->give_shearcoeff (eid,&shearcoeff);
  //  moment of inertia of element cross section
  Mc->give_mominer (eid,&iy);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipid);
  //  Young's modulus
  e=d[0][0];
  //  shear modulus
  g=d[1][1];
  
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;
  
  denom=(1.0+2.0*kappa)*(1.0+2.0*kappa);
  
  fillm (0.0,ism);
  
  ism[1][1]=(6.0/5.0+4.0*kappa+4.0*kappa*kappa)/denom/l;
  ism[1][2]=-1.0/10.0/denom;
  ism[1][4]=(-6.0/5.0-4.0*kappa-4.0*kappa*kappa)/denom/l;
  ism[1][5]=-1.0/10.0/denom;
  
  ism[2][1]=ism[1][2];
  ism[2][2]=(2.0/15.0+kappa/3.0+kappa*kappa/3.0)*l/denom;
  ism[2][4]=1.0/10.0/denom;
  ism[2][5]=(-1.0/30.0-kappa/3.0-kappa*kappa/3.0)*l/denom;
  
  ism[4][1]=ism[1][4];
  ism[4][2]=ism[2][4];
  ism[4][4]=(6.0/5.0+4.0*kappa+4.0*kappa)/denom/l;
  ism[4][5]=1.0/10.0/denom;
  
  ism[5][1]=ism[1][5];
  ism[5][2]=ism[2][5];
  ism[5][4]=ism[4][5];
  ism[5][5]=(2.0/15.0+kappa/3.0+kappa*kappa/3.0)*l/denom;

  
  //  number of integration point
  ipid=Mt->elements[eid].ipp[0][0];
  //  storage of nodal forces and moments
  //  nodal forces and moments are expressed in local element coordinate system
  Mm->givestress (lcid,ipid,0,3,f);
  //  normal force
  nforce = f[0];
  
  cmulm (nforce,ism);
  
  
  
  //  transformation of initial stress matrix from local element coordinate
  //  system to the global coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  lgmatrixtransf (ism,tmat);

  //  transformation of initial stress matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (ism,tmat);
  }
}



/**
   function stores end displacements and rotations to integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 24.5.2006
*/
void beamel2d::nodal_displ (long lcid,long eid)
{
  long i,j,ipid,transf;
  ivector nodes(nne);
  vector x(nne),z(nne),d(ndofe),dd(3),f(ndofe);
  matrix tmat(ndofe,ndofe);
  
  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  
  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  transformation of displacement vector from local nodal systems
  //  to global coordinate system
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    //  assembling of transformation matrix
    transf_matrix (nodes,tmat);
    //  transformation of nodal displacements and rotations
    lgvectortransf (f,d,tmat);
    copyv (f,d);
  }
  
  //  assembling of transformation matrix from global problem coordinate system
  //  to local element coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  //  transformation from global problem coordinate system
  //  to local element coordinate system
  glvectortransf (d,f,tmat);
  copyv (f,d);
  
  j=3;
  for (i=0;i<3;i++){
    dd[i]=d[j];
    j++;
  }
  
  //  number of integration point
  ipid=Mt->elements[eid].ipp[0][0];
  
  //  storage of components
  Mm->storestrain (lcid,ipid,0,3,d);
  Mm->storestrain (lcid,ipid+1,0,3,dd);
}


/**
   function computes nodal forces and moments expressed in local coordinate system
   this is equivalent to functions computing stresses
   nodal forces and moments are better quantities in case of beams

   @param lcid - load case id
   @param eid - element id
   
   JK, 20.2.2002, revised 2.9.2006
*/
void beamel2d::nodal_forces (long lcid,long eid)
{
  long i,j,ipid,transf;
  double l;
  ivector nodes(nne);
  vector x(nne),z(nne),d(ndofe),f(ndofe),ff(3),ifor(ndofe);
  matrix sm(ndofe,ndofe),tmat(ndofe,ndofe);
  
  if (coordsys==1){
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
  }
  if (coordsys==2){
    //  node coordinates in the x-y plane
    Mt->give_node_coord2d (x,z,eid);
  }
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
  }

  //  nodal displacements and rotations
  eldispl (lcid,eid,d.a);
  
  //  assembling of stiffness matrix
  stiffness_matrix_expl (eid,0,0,sm);
  
  //  computation of nodal forces and moments
  mxv (sm,d,f);

  for(i=0; i<Mb->lc[lcid].nle; i++)
  {
      if(Mb->lc[lcid].loe[i].eid != eid)
      continue;

      for(j=0;j<ndofe;j++){
          f[j] -= Mb->lc[lcid].loe[i].nf[j];
      }
  }

  if (Mp->tprob == forced_dynamics || Mp->tprob == eigen_dynamics){
    
    //  mass matrix
    mass_matrix_expl (eid,0,0,sm);
    //  contribution from inertial forces
    mxv (sm,d,ifor);
    cmulv (Lsrs->w[lcid],ifor,ifor);
    
    f[0]-=ifor[0];
    f[1]-=ifor[1];
    f[2]-=ifor[2];

    f[3]-=ifor[3];
    f[4]-=ifor[4];
    f[5]-=ifor[5];
  }


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
  

  //  assembling of transformation matrix from global problem coordinate system
  //  to local element coordinate system
  beam_transf_matrix (eid,x,z,tmat);
  //  transformation from global problem coordinate system
  //  to local element coordinate system
  glvectortransf (f,d,tmat);
  copyv (d,f);
  
  
  j=3;
  for (i=0;i<3;i++){
    ff[i]=-f[j]; // make member force convention on the first node 
    j++;
  }
  
  //  number of integration point
  ipid=Mt->elements[eid].ipp[0][0];
  
  //  storage of nodal forces and moments
  //  nodal forces and moments are expressed in local element coordinate system
  Mm->storestress (lcid,ipid,0,3,f);
  Mm->storestress (lcid,ipid+1,0,3,ff);

}

/**
   Function computes internal forces (nodal forces and moments) 
   that are expressed in global problem or local nodal coordinate 
   system.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   JK, 12.8.2001
*/
void beamel2d::internal_forces (long lcid,long eid,vector &ifor)
{
  ivector nodes(nne);
  vector r(ndofe),contr(ndofe);
  matrix sm(ndofe,ndofe),tmat(ndofe,ndofe);
  
  fillv (0.0,ifor);

  //  nodal displacements and rotations
  eldispl (lcid,eid,r.a);
  //  assembling of stiffness matrix
  stiffness_matrix_expl (eid,0,0,sm);
  //  computation of nodal forces and moments
  mxv (sm,r,contr);
  
  //  nodal forces and moments are expressed in global problem or local nodal coordinate system
  addv (contr,ifor,ifor);


  
  
  //  ???!!! vypnuto ukladani koncovych sil vlivem posunu a natoceni - 
  //  ukladani provadi pouze fce nodal_forces, ktera
  //  zahrnuje i prispevek od mimostycnikoveho zatizeni + transformaci do 
  // lokalniho souradneho systemu prutu.
  // TKo
/*
  long i,j,ipid;
  vector ff(3);
  //  number of integration point
  ipid=Mt->elements[eid].ipp[0][0];
  j=3;
  for (i=0;i<3;i++){
    ff[i]=ifor[j];
    j++;
  }
  //  storage of nodal forces and moments
  Mm->storestress (lcid,ipid,0,3,ifor);
  Mm->storestress (lcid,ipid+1,0,3,ff);
*/
}

/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 12.8.2001
*/
void beamel2d::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,ifor);
}

/**
   function defines meaning of DOFs at nodes
   
   @param eid - element id
   
   1.2.2005, JK
*/
void beamel2d::define_meaning (long eid)
{
  ivector cn(ndofe),nod(nne);
  
  Mt->give_elemnodes (eid,nod);
  Mt->give_code_numbers (eid,cn.a);
  
  if (coordsys==1){
    //  displacement in x direction
    if (cn[0]>0)  Mt->nodes[nod[0]].meaning[0]=1;
    //  displacement in z direction
    if (cn[1]>0)  Mt->nodes[nod[0]].meaning[1]=3;
    //  displacement in x direction
    if (cn[3]>0)  Mt->nodes[nod[1]].meaning[0]=1;
    //  displacement in z direction
    if (cn[4]>0)  Mt->nodes[nod[1]].meaning[1]=3;
  }
  if (coordsys==2){
    //  displacement in x direction
    if (cn[0]>0)  Mt->nodes[nod[0]].meaning[0]=1;
    //  displacement in y direction
    if (cn[1]>0)  Mt->nodes[nod[0]].meaning[1]=2;
    //  displacement in x direction
    if (cn[3]>0)  Mt->nodes[nod[1]].meaning[0]=1;
    //  displacement in y direction
    if (cn[4]>0)  Mt->nodes[nod[1]].meaning[1]=2;
  }

}



/**
   function assembles nodal forces and moments caused by a load acting on element
   
   @param eid - element id
   @param le - type of coordinate system
   @param nv - values defining load
   @param nf - %vector of end forces and moments
   
   30. 9. 2012
*/
void beamel2d::nodeforces(long eid, long *le, double*nv, vector &nf)
{
  //  this is mechanical load
  /*
  if(le[0]==1){
    double fxa=nv[0];  
    double fza=nv[1];
    double mya=nv[2];  
    double fxb=nv[3];
    double fzb=nv[4];  
    double myb=nv[5];
    
    double l;
    vector x(nne),z(nne);
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
    
    //  length of the element
    l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
    if (l<Mp->zero){
      print_err("zero length of the %ld beamel2d element",__FILE__,__LINE__,__func__,eid);
    }
    
    nf[0]= 1*(((fxa*l)/2) + (((fxb-fxa)*l)/6));
    nf[3]= 1*(((fxa*l)/2) + (((fxb-fxa)*l)/3));
    nf[2]= -1*(((fza*l)/2) + (((fzb-fza)*3*l)/20)) -1*(((18*mya)/l) + (((myb-mya)*9)/l));
    nf[4]= -1*(((fza*l)/2) + (((fzb-fza)*7*l)/20)) + ((18*mya)/l) + (((myb-mya)*9)/l);
    nf[2]= ((fza*l*l)/12)+(((fzb-fza)*l*l)/30) + 6*mya +3*(myb-mya);
    nf[5]= -1*(((fza*l*l)/12)+(((fzb-fza)*l*l)/20)) + 12*mya +6*(myb-mya);
  }
  */
  
  //  this thermal load
  if (le[0]==1){
    long ipid;
    double dts,dt,e,a,iy,alpha,h;
    matrix d(tncomp,tncomp);
   
    //  integration point id
    ipid=Mt->elements[eid].ipp[0][0];
    //  area of element cross section
    Mc->give_area (eid,a);
    //  moment of inertia of element cross section
    Mc->give_mominer (eid,&iy);
    //  stiffness matrix of the material
    Mm->matstiff (d,ipid);
    //  Young's modulus
    e=d[0][0];
    
    //  change of temperature in the center of cross section
    dts = nv[0];
    //  difference of temperature changes of lower and upper surface of the beam
    dt  = nv[1];
    alpha = nv[2];
    h = nv[3];
    
    nf[0]=0.0-e*a*alpha*dts;
    nf[1]=0.0;
    nf[2]=0.0-e*iy*alpha/h*dt;
    nf[3]=e*a*alpha*dts;
    nf[4]=0.0;
    nf[5]=e*iy*alpha/h*dt;
    
    vector x(nne),z(nne);
    matrix tmat (ndofe,ndofe);
    //  node coordinates in the x-z plane
    Mt->give_node_coord2dxz (x,z,eid);
    //  system to the global coordinate system
    beam_transf_matrix (eid,x,z,tmat);
    lgvectortransfblock (nf,tmat);
  }
  
  if(Gtm->gelements[eid].cne==2){
    //
    matrix sm(ndofe,ndofe);
    stiffness_matrix_expl (eid,0,0,sm);
    
    ivector cu(ndofe);
    Gtm->give_cn(eid, cu.a);
    condense_vector(sm,nf,cu);
  }
}



/**
   function assembles end forces and moments cuased by a beam load type acting on element
   
   @param eid - element id
   @param le - type of coordinate system
   @param nv - values defining load
   @param nf - %vector of end forces and moments
   
   30. 9. 2012
*/
void beamel2d::beamnodeforces(long eid, elloadmeaning elm, double /*nnv*/, double *la, double *lf, double*nv, vector &nf)
{
  double l;

  switch(elm){
  case load_contmech_glob:{      
    //  this is mechanical load defined in the global coordinate system
    vector x(ASTCKVEC(nne)), z(ASTCKVEC(nne));
    
    if (coordsys==1){
      //  node coordinates in the x-z plane
      Mt->give_node_coord2dxz (x,z,eid);
    }
    if (coordsys==2){
      //  node coordinates in the x-y plane
      Mt->give_node_coord2d (x,z,eid);
      print_err("beam in x-y is not implemented",__FILE__, __LINE__, __func__);
      abort();
    }
    
    //  length of the element
    l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
    if (l<Mp->zero){
      print_err("zero length of the %ld beamel2d element",__FILE__, __LINE__, __func__,eid);
      abort();
    }
    
    //  transformation matrix from local coordinate system to the global coordinate system      
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    beam_transf_matrix (eid,x,z,tmat);
    // nodal values of load in local coordinate system
    vector nvl(ASTCKVEC(nf.n));      
    copyv(nv, nvl); // initialize with global values
    // compute nodal forces in the local coord. system of the beam
    double fxa=nvl[0];  
    double fza=nvl[1];
    double mya=nvl[2];  
    double fxb=nvl[3];
    double fzb=nvl[4];  
    double myb=nvl[5];
    double lb, d11, d12, d22, d1f, d2f, d1m, d2m, det_a, det_b1, det_b2, det_b1m, det_b2m, x1, x2, x1m, x2m;
    // end forces due to continuous load in x direction
    lb = l - lf[0] - la[0];
    nf[0] = -sqr(lf[0])*(2*fxa+fxb)/(6.0*l) - lf[0]*lb*(fxa+fxb)/(2.0*l);
    nf[3] = -sqr(lf[0])*(2.0*fxa+fxb)/(6.0*l) - lf[0]*lb*(fxa+fxb)/(2.0*l) - (fxa+fxb)*lf[0]/2.0;
    // end forces due to continuous load in z direction
    lb = l - lf[1] - la[1];
    d11 = (lf[1]+la[1])*sqr(lf[1]+la[1]) + lb*lb*lb;
    d22 = (lf[1]+lb)*sqr(lf[1]+lb) + la[1]*la[1]*la[1];
    d12 = -sqr(la[1])*(3.0*lf[1]+2.0*la[1])/(6.0*sqr(lf[1])) + lf[1]/6.0 - sqr(lb)*(3.0*lf[1]+2.0*lb)/(sqr(lf[1])*6.0);
    d1f = -sqr(la[1])*(3.0*lf[1]+2.0*la[1])*(2.0*fza+fzb)/36.0 + lf[1]*sqr(lf[1])*(8.0*fza+7.0*fzb)/360.0 + lb*sqr(lb)*(fza+2.0*fzb)/18.0;
    d2f = la[1]*sqr(la[1])*(2.0*fza+fzb)/18.0 + lf[1]*sqr(lf[1])*(7.0*fza+8.0*fzb)/360.0 - sqr(lb)*(3.0*lf[1]+2.0*lb)*(fza+2.0*fzb)/360.0;
    det_a = d11*d22-sqr(d12);
    det_b1 = d12*d2f-d22*d1f;
    det_b2 = d12*d1f-d11*d2f;
    x1 = det_b1/det_a;
    x2 = det_b2/det_a;
    nf[2] = -(1.0+la[1]/lf[1])*x1 + la[1]/lf[1]*x2 + la[1]*lf[1]*(2.0*fza+fzb)/6.0;
    nf[5] = -lb/lf[1]*x1 + (1.0+lb/lf[1])*x2 - lb*lf[1]*(fza+2.0*fzb)/6.0;
    nf[1] = (x1-x2)/lf[1] - 0.5*fza*lf[1] - (fzb-fza)*lf[1]/6.0;
    nf[4] = (x2-x1)/lf[1] - 0.5*fza*lf[1] - (fzb-fza)*lf[1]/3.0;
    // end forces due to continuous moment load in y direction
    lb = l - lf[2] - la[2];
    d1m = -(mya+myb)*sqr(la[2])*(3.0*lf[2]+2.0*la[2])/(12.0*lf[2]) + sqr(lf[2])*(myb-mya)/24.0 - (mya+myb)*sqr(lb)*lb/(6.0*lf[2]);
    d2m = (mya+myb)*la[2]*sqr(la[2])/(6.0*lf[2])+(myb-mya)*sqr(lf[2])/24.0+sqr(lb)*(mya+myb)*(3.0*lf[2]+2.0*lb)/(12.0*lf[2]);
    det_a = d11*d22-sqr(d12);
    det_b1m = d12*d2f-d22*d1f;
    det_b2m = d12*d1f-d11*d2f;
    x1m = det_b1m/det_a;
    x2m = det_b2m/det_a;
    nf[2] += -(1.0+la[2]/lf[2])*x1m + la[2]/lf[2]*x2m + (mya+myb)*la[2]/2.0;
    nf[5] += -lb/lf[2]*x1m + (1.0+lb/lf[2])*x2m + (mya+myb)*lb/2.0;
    nf[1] += (x1m-x2m)/lf[2] - (mya+myb)/2.0;
    nf[4] += (x2m-x1m)/lf[2] + (mya+myb)/2.0;
    // transform nodal forces to the global coord. system
    lgvectortransfblock (nf,tmat);
    break;
  }
  case load_contmech_loc:{
    //  this is mechanical load defined in the local coordinate system of the beam
    //  this is mechanical load defined in the global coordinate system
    vector x(ASTCKVEC(nne)), z(ASTCKVEC(nne));

    if (coordsys==1){
      //  node coordinates in the x-z plane
      Mt->give_node_coord2dxz (x,z,eid);
    }
    if (coordsys==2){
      //  node coordinates in the x-y plane
      Mt->give_node_coord2d (x,z,eid);
      print_err("beam in x-y is not implemented",__FILE__, __LINE__, __func__);
      abort();
    }
    
    //  length of the element
    double l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
    if (l < Mp->zero){
      print_err("zero length of the %ld beamel2d element", __FILE__, __LINE__, __func__, eid+1);
      abort();
    }
    // compute nodal forces in the local coord. system of the beam
    double fxa=nf[0];  
    double fza=nf[1];
    double mya=nf[2];  
    double fxb=nf[3];
    double fzb=nf[4];  
    double myb=nf[5];
    double lb, d11, d12, d22, d1f, d2f, d1m, d2m, det_a, det_b1, det_b2, det_b1m, det_b2m, x1, x2, x1m, x2m;
    // end forces due to continuous load in x direction
    lb = l - lf[0] - la[0];
    nf[0] = -sqr(lf[0])*(2*fxa+fxb)/(6.0*l) - lf[0]*lb*(fxa+fxb)/(2.0*l);
    nf[3] = -sqr(lf[0])*(2.0*fxa+fxb)/(6.0*l) - lf[0]*lb*(fxa+fxb)/(2.0*l) - (fxa+fxb)*lf[0]/2.0;
    // end forces due to continuous load in z direction
    lb = l - lf[1] - la[1];
    d11 = (lf[1]+la[1])*sqr(lf[1]+la[1]) + lb*lb*lb;
    d22 = (lf[1]+lb)*sqr(lf[1]+lb) + la[1]*la[1]*la[1];
    d12 = -sqr(la[1])*(3.0*lf[1]+2.0*la[1])/(6.0*sqr(lf[1])) + lf[1]/6.0 - sqr(lb)*(3.0*lf[1]+2.0*lb)/(sqr(lf[1])*6.0);
    d1f = -sqr(la[1])*(3.0*lf[1]+2.0*la[1])*(2.0*fza+fzb)/36.0 + lf[1]*sqr(lf[1])*(8.0*fza+7.0*fzb)/360.0 + lb*sqr(lb)*(fza+2.0*fzb)/18.0;
    d2f = la[1]*sqr(la[1])*(2.0*fza+fzb)/18.0 + lf[1]*sqr(lf[1])*(7.0*fza+8.0*fzb)/360.0 - sqr(lb)*(3.0*lf[1]+2.0*lb)*(fza+2.0*fzb)/360.0;
    det_a = d11*d22-sqr(d12);
    det_b1 = d12*d2f-d22*d1f;
    det_b2 = d12*d1f-d11*d2f;
    x1 = det_b1/det_a;
    x2 = det_b2/det_a;
    nf[2] = -(1.0+la[1]/lf[1])*x1 + la[1]/lf[1]*x2 + la[1]*lf[1]*(2.0*fza+fzb)/6.0;
    nf[5] = -lb/lf[1]*x1 + (1.0+lb/lf[1])*x2 - lb*lf[1]*(fza+2.0*fzb)/6.0;
    nf[1] = (x1-x2)/lf[1] - 0.5*fza*lf[1] - (fzb-fza)*lf[1]/6.0;
    nf[4] = (x2-x1)/lf[1] - 0.5*fza*lf[1] - (fzb-fza)*lf[1]/3.0;
    // end forces due to continuous moment load in y direction
    lb = l - lf[2] - la[2];
    d1m = -(mya+myb)*sqr(la[2])*(3.0*lf[2]+2.0*la[2])/(12.0*lf[2]) + sqr(lf[2])*(myb-mya)/24.0 - (mya+myb)*sqr(lb)*lb/(6.0*lf[2]);
    d2m = (mya+myb)*la[2]*sqr(la[2])/(6.0*lf[2])+(myb-mya)*sqr(lf[2])/24.0+sqr(lb)*(mya+myb)*(3.0*lf[2]+2.0*lb)/(12.0*lf[2]);
    det_a = d11*d22-sqr(d12);
    det_b1m = d12*d2f-d22*d1f;
    det_b2m = d12*d1f-d11*d2f;
    x1m = det_b1m/det_a;
    x2m = det_b2m/det_a;
    nf[2] += -(1.0+la[2]/lf[2])*x1m + la[2]/lf[2]*x2m + (mya+myb)*la[2]/2.0;
    nf[5] += -lb/lf[2]*x1m + (1.0+lb/lf[2])*x2m + (mya+myb)*lb/2.0;
    nf[1] += (x1m-x2m)/lf[2] - (mya+myb)/2.0;
    nf[4] += (x2m-x1m)/lf[2] + (mya+myb)/2.0;
    // transform nodal forces to the global coord. system
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    beam_transf_matrix (eid,x,z,tmat);
    lgvectortransfblock (nf,tmat);
    break;
  }
  case load_thermgrad_z:{
    //  this thermal load
    long ipid;
    double dts,dt,e,a,iy,alpha,h;
    matrix d(tncomp,tncomp);
    
    //  integration point id
    ipid=Mt->elements[eid].ipp[0][0];
    //  area of element cross section
    Mc->give_area (eid,a);
    //  moment of inertia of element cross section
    Mc->give_mominer (eid,&iy);
    //  stiffness matrix of the material
    Mm->matstiff (d,ipid);
    //  Young's modulus
    e=d[0][0];
    
    //  change of temperature in the center of cross section
    dts = nv[0];
    //  difference of temperature changes of lower and upper surface of the beam
    dt  = nv[1];
    alpha = nv[2];
    h = nv[3];
    
    nf[0]=0.0-e*a*alpha*dts;
    nf[1]=0.0;
    nf[2]=0.0-e*iy*alpha/h*dt;
    nf[3]=e*a*alpha*dts;
    nf[4]=0.0;
    nf[5]=e*iy*alpha/h*dt;
    
    vector x(ASTCKVEC(nne)),z(ASTCKVEC(nne));
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    
    if (coordsys==1){
      //  node coordinates in the x-z plane
      Mt->give_node_coord2dxz (x,z,eid);
    }
    if (coordsys==2){
      //  node coordinates in the x-y plane
      Mt->give_node_coord2d (x,z,eid);
      print_err("beam in x-y is not implemented",__FILE__, __LINE__, __func__);
      abort();
    }
    
    //  system to the global coordinate system
    beam_transf_matrix (eid,x,z,tmat);
    lgvectortransfblock (nf,tmat);
    break;
  }
  default:
    print_err("unknown meaning %d of the element load is required on element %ld", __FILE__, __LINE__, __func__, (int)elm, eid+1);
    abort();
  }
  if(Gtm->gelements[eid].cne==2){
    //
    matrix sm(ndofe,ndofe);
    stiffness_matrix_expl (eid,0,0,sm);
    
    ivector cu(ndofe);
    Gtm->give_cn(eid, cu.a);
    condense_vector(sm,nf,cu);
  }
}
