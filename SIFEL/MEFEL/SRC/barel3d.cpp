#include <math.h>
#include "barel3d.h"
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
#include "loadcase.h"



/**
   JK, 27.9.2005
*/
barel3d::barel3d (void)
{
  long i;
  
  //  number nodes on element
  nne=2;
  //  number of DOFs on element
  ndofe=6;
  //  number of strain/stress components
  tncomp=1;
  //  number of functions approximated
  napfun=3;
  ssst = bar;
  
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
  
  nip[0][0]=1;

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    tnip+=nip[i][i];
  }
  
  intordsm[0][0]=1;
  
  intordmm=2;
}



barel3d::~barel3d (void)
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
   Function computes direction %vector.
   
   @param s - direction %vector
   @param x,y,z - vectors of nodal coordinates
   
   @retval Function returns direction %vector s.

   JK, 28.9.2005
*/
void barel3d::dirvect (vector &s,vector &x,vector &y,vector &z)
{
  double l;

  s[0]=x[1]-x[0];
  s[1]=y[1]-y[0];
  s[2]=z[1]-z[0];
  
  l=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
  
  if (l<Mp->zero){
    print_err("zero norm of direction vector",__FILE__,__LINE__,__func__);
  }
  
  s[0]/=l;
  s[1]/=l;
  s[2]/=l;
}

/**
  Function returns coordinates of nodes at local element coordinate system.

  @param x  - %vector of nodal x-coordinates in global coordinate system
  @param y  - %vector of nodal y-coordinates in global coordinate system
  @param z  - %vector of nodal z-coordinates in global coordinate system
  @param lx - %vector of nodal coordinates at local element coordinate system

  @retval Function returns %vector of nodal coordinates at parameter lx.
*/
void barel3d::giveloccoord(vector &x, vector &y, vector &z, vector &lx)
{
  double l;

  lx[0] = 0.0;
  l = sqr(x[1] - x[0]) + sqr(y[1] - y[0]) + sqr(z[1] - z[0]);
  l = sqrt(l);
  lx[1] = l;
}

/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param xi - natural coordinate
   
   JK, 25.8.2001
*/
void barel3d::bf_matrix (matrix &n,vector &s,double xi)
{
  vector bf(ASTCKVEC(nne));
  
  bf_lin_1d (bf.a,xi);
  
  n[0][0]=bf[0]*s[0];
  n[0][1]=bf[0]*s[1];
  n[0][2]=bf[0]*s[2];
  n[0][3]=bf[1]*s[0];
  n[0][4]=bf[1]*s[1];
  n[0][5]=bf[1]*s[2];

}


/**
   function assembles strain-displacement (geometric) %matrix

   @param gm - strain-displacement (geometric) %matrix
   @param x,y,z - global coordinates of nodes
   @param jac - jacobian

   @retval Function returns geometric %matrix in parameter gm. 

   JK, 27.9.2005
*/
void barel3d::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,double &jac)
{
  double xi;
  vector s(ASTCKVEC(3)),lx(ASTCKVEC(nne)),dx(ASTCKVEC(nne));
  
  //  direction vector
  dirvect (s,x,y,z);
  //  local coordinates of nodes
  giveloccoord (x,y,z,lx);
  
  //  derivatives of the approximation functions
  //  the derivatives are independent on xi in this case
  xi=0.0;
  derivatives_1d (dx,jac,lx,xi);

  gm[0][0]=dx[0]*s[0];
  gm[0][1]=dx[0]*s[1];
  gm[0][2]=dx[0]*s[2];
  gm[0][3]=dx[1]*s[0];
  gm[0][4]=dx[1]*s[1];
  gm[0][5]=dx[1]*s[2];
}



/**
   function assembles strain-displacement (geometric) %matrix

   @param gm - strain-displacement (geometric) %matrix [out]
   @param xi - natural coordinate [in]

   @retval Function returns geometric %matrix at point with natural coordinate xi in parameter gm. 

   TKo, 22.11.2017
*/
/*
void barel3d::geom_matrix (matrix &gm,double xi)
{
  gm[0][0]=-xi;
  gm[0][1]=-xi;
  gm[0][2]=-xi;
  gm[0][3]=xi;
  gm[0][4]=xi;
  gm[0][5]=xi;
}
*/


/**
   Function approximates function defined by nodal values.

   @param xi - coordinate on element
   @param nodval - nodal values
   
   @retval Function returns approximated value for the given natural coordinate xi.

   JK, 27.9.2005
*/
double barel3d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));

  bf_lin_1d (bf.a,xi);

  scprd (bf,nodval,f);

  return f;
}



/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   @retval Function returns transformation %matrix in parameter tmat.

   JK, 3.2.2002
*/
void barel3d::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i;
  fillm (0.0,tmat);
  
  for (i=0;i<tmat.m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<nodes.n;i++){
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
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   T x_local = x_global
   x_local = T^T x_global
   
   direction %vector denotes %vector defining local coordinate axis
   
   @param eid[in] - element id
   @param x[in] - %vector of nodal coordinates in element coordinate system
   @param vec[in]   - direction %vector of local axis z (in global coordinate system)
   @param tmat[out] - tranformation %matrix T
   @param gx[in] - vector of nodal x coordinates in global coordinate system
   @param gy[in] - vector of nodal y coordinates in global coordinate system
   @param gz[in] - vector of nodal z coordinates in global coordinate system
   
*/
void barel3d::bar_transf_mat(long eid, vector &x, vector &vec, matrix &tmat, vector &gx, vector &gy, vector &gz)
{
  double c, l;

  nullm(tmat);
  
  //  local vector x_l
  tmat[0][0]=gx[1]-gx[0];
  tmat[1][0]=gy[1]-gy[0];
  tmat[2][0]=gz[1]-gz[0];
  
  //  length of the beam, it is equal to the norm of x_l
  l=sqrt((tmat[0][0]*tmat[0][0]+tmat[1][0]*tmat[1][0]+tmat[2][0]*tmat[2][0]));

  if (l<Mp->zero){
    print_err("zero length of the %ld barel3d element", __FILE__, __LINE__, __func__, eid);
  }

  x(0) = 0.0; x(1) = l;
  
  //  normed local vector x_l
  tmat[0][0]=tmat[0][0]/l;
  tmat[1][0]=tmat[1][0]/l;
  tmat[2][0]=tmat[2][0]/l;
  
  //  local vector y_l
  tmat[0][1]=vec[1]*tmat[2][0]-vec[2]*tmat[1][0];
  tmat[1][1]=vec[2]*tmat[0][0]-vec[0]*tmat[2][0];
  tmat[2][1]=vec[0]*tmat[1][0]-vec[1]*tmat[0][0];
  
  //  norm of the local vector y_l
  c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));

  if (c < Mp->zero){
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
    
    if (c < Mp->zero){
      print_err("zero z base vector of the %ld barel3d element",__FILE__,__LINE__,__func__,eid);
    }
    
    //  local vector y_l
    tmat[0][1]=tmat[1][2]*tmat[2][0]-tmat[2][2]*tmat[1][0];
    tmat[1][1]=tmat[2][2]*tmat[0][0]-tmat[0][2]*tmat[2][0];
    tmat[2][1]=tmat[0][2]*tmat[1][0]-tmat[1][2]*tmat[0][0];

    //  norm of the local vector y_l
    c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));

    if (c < Mp->zero){
      print_err("zero y base vector of the %ld barel3d element",__FILE__,__LINE__,__func__,eid);
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
    
    if (c < Mp->zero){
      print_err("zero z base vector of the %ld barel3d element",__FILE__,__LINE__,__func__,eid);
    }
    
    //  normed local vector z_l
    tmat[0][2]=tmat[0][2]/c;
    tmat[1][2]=tmat[1][2]/c;
    tmat[2][2]=tmat[2][2]/c;
  }
  
  tmat[3][3]  = tmat[0][0];  tmat[3][4]   = tmat[0][1];  tmat[3][5]   = tmat[0][2];
  tmat[4][3]  = tmat[1][0];  tmat[4][4]   = tmat[1][1];  tmat[4][5]   = tmat[1][2];
  tmat[5][3]  = tmat[2][0];  tmat[5][4]   = tmat[2][1];  tmat[5][5]   = tmat[2][2];
}



/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   T x_local = x_global
   x_local = T^T x_global
   
   direction %vector denotes %vector defining local coordinate axis
   
   @param x - %vector of nodal coordinates in element coordinate system
   @param tmat - tranformation %matrix T
   @param gx,gy - vectors of nodal coordinates in global coordinate systems
   
*/
/*
void barel3d::tran_mat (vector &x, matrix &tmat,vector &gx,vector &gy)
{
  long i,j;
  vector gv(2*gx.n),lv(2*gx.n);
  
  j=0;
  for (i=0;i<gx.n;i++){
    gv[j]=gx[i];  j++;
    gv[j]=gy[i];  j++;
  }
  
  //globloctransf (gv,lv,tmat);
  glvectortransf (gv,lv,tmat);

  for (i=0;i<gx.n;i++){
    x[i]=lv[i*2];
  }
}
*/


/**
   Function computes stiffness %matrix of threedimensonal bar element.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y,z - coordinates of nodes
   
   JK, 27.9.2005
*/
void barel3d::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y,vector &z)
{
  long ipp;
  double w,a,jac;
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));
  
  fillm (0.0,sm);
  
  //  the strain-displacement matrix is a constant matrix
  //  therfore, one point integration is enough
  
  //  strain-displacement matrix
  geom_matrix (gm,x,y,z,jac);
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of material
  Mm->matstiff (d,ipp);
  //  are of cross section
  Mc->give_area (eid,a);
  
  w=2.0;
  jac*=a*w;
  
  //  B^T D B
  bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
}



/**
   Function computes stiffness %matrix of one element. If it is required, nodal values are transformed to 
   the local coordinate systems.
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   @retval Function returns resulting stiffness %matrix in parameter sm

   JK, 27.9.2005
*/
void barel3d::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  // node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  stiffness matrix of element
  stiffness_matrix (eid,0,0,sm,x,y,z);


  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}



/**
   Function computes mass %matrix of threedimensional bar element.
   
   @param eid - element id
   @param mm - mass %matrix
   @param x,y,z - coordinates of nodes
   
   @retval Function returns resulting mass %matrix in parameter mm

   JK, 27.9.2005
*/
void barel3d::mass_matrix (long eid,matrix &mm,vector &x,vector &y,vector &z)
{
  long i;
  double a,rho,jac,xi;
  ivector nodes(ASTCKIVEC(nne));
  vector lx(ASTCKVEC(nne)),s(ASTCKVEC(3));
  vector dens(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm));
  matrix n(ASTCKMAT(1,ndofe));
  
  fillm (0.0,mm);
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);

  //  local coordinates of nodes
  //  only x coordinates are nonzero
  giveloccoord(x,y,z,lx);
  //  direction vector
  dirvect (s,x,y,z);
  
  //  area of bar cross-section
  Mc->give_area (eid,a);
  //  density of the material
  Mc->give_density (eid,nodes,dens);
  
  //  Gauss points
  gauss_points (gp.a,w.a,intordmm);
  
  
  for (i=0;i<intordmm;i++){
    xi = gp[i];
    
    //  jacobian
    jac_1d (jac,lx,xi);
 
    //  matrix of approximation functions
    bf_matrix (n,s,xi);
    //  density in integration point
    rho = approx (xi,dens);
    
    jac*=rho*a*w[i];
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
}



/**
   Function computes mass %matrix of one element. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param eid - number of element
   @param mm  - mass %matrix

   @retval Function returns resulting mass %matrix in parameter mm
   
   TKo 3.2009
   JK, 25. 7. 2018
*/
void barel3d::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix tmat (ASTCKMAT(ndofe,ndofe));
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);  
  
  //  mass matrix
  mass_matrix (eid,mm,x,y,z);
  
  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}





/**
   Function computes load %matrix.
   
   @param eid[in] - element id
   @param lm[out] - load %matrix
   @param x[in] - %vector of node coordinates
   
   @retval Function returns load %matrix in parameter lm.

   TKo according to JK, 1.5.2023
*/
void barel3d::load_matrix (long eid, matrix &lm, vector &x)
{
  double l,a;
  
  fillm (0.0,lm);
  
  //  length of the element
  l=sqrt((x[1]-x[0])*(x[1]-x[0]));
  if (l < Mp->zero){
    print_err("zero length of the %ld truss element",__FILE__, __LINE__, __func__, eid+1);
  }
  
  //  area of bar cross-section
  Mc->give_area(eid, a);

  // load is considered to be approximated by linear functions
  lm[0][0]=a*l*1.0/3.0;
  lm[0][1]=0.0;
  lm[0][2]=0.0;
  lm[0][3]=a*l*1.0/6.0;
  lm[0][4]=0.0;
  lm[0][4]=0.0;
  
  lm[1][0]=0.0;
  lm[1][1]=a*l*1.0/3.0;
  lm[1][2]=0.0;
  lm[1][3]=0.0;
  lm[1][4]=a*l*1.0/6.0;
  lm[1][5]=0.0;
  
  lm[2][0]=0.0;
  lm[2][1]=0.0;
  lm[2][2]=a*l*1.0/3.0;
  lm[2][3]=0.0;
  lm[2][4]=0.0;
  lm[2][5]=a*l*1.0/6.0;
  
  lm[3][0]=a*l*1.0/3.0;
  lm[3][1]=0.0;
  lm[3][2]=0.0;
  lm[3][3]=a*l*1.0/6.0;
  lm[3][4]=0.0;
  lm[3][4]=0.0;
  
  lm[4][0]=0.0;
  lm[4][1]=a*l*1.0/3.0;
  lm[4][2]=0.0;
  lm[4][3]=0.0;
  lm[4][4]=a*l*1.0/6.0;
  lm[4][5]=0.0;
  
  lm[5][0]=0.0;
  lm[5][1]=0.0;
  lm[5][2]=a*l*1.0/3.0;
  lm[5][3]=0.0;
  lm[5][4]=0.0;
  lm[5][5]=a*l*1.0/6.0;
}



/**
  The function computes load %matrix of the r
   finite element with bilinear approximation functions.
   Load vector is obtained after premultiplying load %matrix
   by nodal load values.
   
   @param eid[in] - number of element
   @param lm[out] - load %matrix
   
   TKo according to JK, 1.5.2023
*/
void barel3d::res_load_matrix(long eid, matrix &lm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), gx(ASTCKVEC(nne)), gy(ASTCKVEC(nne)), gz(ASTCKVEC(nne));
  vector lz(ASTCKVEC(3));
  matrix tran(ASTCKMAT(3,3));

  // ???!!! provisional direction vector of local z axis
  nullv(lz);
  lz(2) = 1.0;
  //  node coordinates in the global coordinate system
  Mt->give_node_coord3d(gx, gy, gz, eid);
  //  transformation matrix from global to element coordinate system
  bar_transf_mat(eid, x, lz, tran, gx, gy, gz);
  
  //  load matrix
  load_matrix(eid, lm, x);
  
  //  transformation of load matrix to global coordinate system
  lgmatrixtransfblock(lm, tran);
  
  //  transformation of load matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix(nodes, tmat);
    glmatrixtransf(lm, tmat);
  }
}




/**
   Function computes strains at integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 4.10.2005
*/
void barel3d::res_ip_strains (long lcid,long eid)
{
  long transf;
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne)),gz(ASTCKVEC(nne));
  vector s(ASTCKVEC(3)),r(ASTCKVEC(ndofe)),aux(ASTCKVEC(ndofe));
  ivector cn(ASTCKIVEC(ndofe)),nodes(ASTCKIVEC(nne));
  matrix tmat(ASTCKMAT(ndofe,ndofe));
  
  //  node coordinates
  Mt->give_node_coord3d (gx,gy,gz,eid);
  //  direction vector
  dirvect (s,gx,gy,gz);
  
  //  assembling of displacements defined on the element
  eldispl (lcid,eid,r.a);

  
  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    //  transformation matrix
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  ip_strains (lcid,eid,0,0,s,r);
}



/**
   Function computes strains at integration points of element.

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param s - direction %vector
   @param r - vector of nodal displacements
   
   JK, 3.10.2005
*/
void barel3d::ip_strains (long lcid,long eid,long ri,long ci,vector &/*s*/,vector &r)
{
  long i,ipid;
  double jac;
  vector eps(ASTCKVEC(tncomp));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  //  coordinates of nodes
  Mt->give_node_coord3d (x,y,z,eid);

  ipid=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    
    //  geometrix matrix
    geom_matrix (gm,x,y,z,jac);
    //  strain computation
    mxv (gm,r,eps);
    
    //  strain storage
    Mm->storestrain (lcid,ipid,eps);
    ipid++;
  }
}



/**
   Function computes strains in nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

   JK, 10.5.2002
*/
void barel3d::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipid;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector eps(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[ri][ci];
  nodip_bar (ipid,intordsm[0][0],ipnum);
  
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
   Function computes strains at strain points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void barel3d::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  vector coord,eps;

  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_strains (lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    nod_strains_ip (lcid,eid,ri,ci);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (RSTCKVEC(ncp,eps));
    reallocv (RSTCKVEC(1,coord));
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	//appstrain (lcid,eid,coord[0],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:
    print_err("unknown strain point is required",__FILE__,__LINE__,__func__);
  }

}



/**
   Function computes stresses at integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barel3d::res_ip_stresses (long lcid,long eid)
{
  ip_stresses (lcid,eid,0,0);
}



/**
   Function computes stresses at integration points of element.

   @param eid - element id
   @param ri - row index
   @param ci - column index

   JK, 10.5.2002
*/
void barel3d::ip_stresses (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipid,Mm->ip[ipid]);
}



/**
   Function computes stresses at integration points of element.

   @param eid - element id
   @param ri - row index
   @param ci - column index

   JK, 10.5.2002
*/
void barel3d::ip_elast_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ipid;
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ipid=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    //  stiffness matrix of material
    Mm->matstiff (d,ipid);
    //  strain
    Mm->givestrain (lcid,ipid,eps);
    //  stress computation
    mxv (d,eps,sig);
    //  storage of stress
    Mm->storestress (lcid,ipid,sig);
    ipid++;
  }
}



/**
   Function computes stresses at nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void barel3d::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipid;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[ri][ci];
  nodip_bar (ipid,intordsm[0][0],ipnum);
  
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
   Function computes stresses at required positions.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

*/
void barel3d::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  vector coord,sig;

  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    nod_stresses_ip (lcid,eid,ri,ci);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (RSTCKVEC(ncp,sig));
    reallocv (RSTCKVEC(1,coord));
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],0,ncp,sig);

      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:
    print_err("unknown stress point is required",__FILE__,__LINE__,__func__);
  }
}



/**
   Function computes other values in nodes of element.

   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void barel3d::nod_eqother_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipid,ncompo;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[ri][ci];
  nodip_bar (ipid,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Mm->givencompeqother (ipnum[i],0);
    reallocv (RSTCKVEC(ncompo,eqother));
    Mm->giveeqother (ipnum[i],0,ncompo,eqother.a);
    
    //  storage of other values to the node
    j=nod[i];
    Mt->nodes[j].storeother (lcid,0,ncompo,eqother);
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

   TKo, 7.2008
*/
void barel3d::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long ii,jj,i,aipp;
  double xi;
  vector w,gp;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));

  Mt->give_node_coord3d (x,y,z,eid);
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      //  integration point id
      aipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        if (aipp==ipp){
          coord[0]=approx (xi,x);
          coord[1]=approx (xi,y);
          coord[2]=approx (xi,z);
        }
        aipp++;
      }
      destrv (gp);  destrv (w);
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
void barel3d::ipncoord (long eid, long ipp, vector &ncoord)
{
  long ii, jj, i, aipp;
  double xi;
  vector w, gp;

  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      //  integration point id
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ii][jj],gp));
      reallocv(RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      aipp=Mt->elements[eid].ipp[ii][jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        if (aipp==ipp){
          ncoord[0]=xi;
          ncoord[1]=0.0;
          ncoord[2]=0.0;
        }
        aipp++;
      }
    }
  }
}



/**
   Function computes volume appropriate to integration point.

   @param eid - element id
   @param ri  - row index of the given integration point block
   @param ci  - column index of the given integration point block

   TKo, 7.2008
*/
void barel3d::ipvolume (long eid,long ri,long ci)
{
  long ipp;
  double l,a,vol;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  // node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + (z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    char emsg[100];
    sprintf(emsg, "zero length of the %ld truss element", eid);
    print_err(emsg, __FILE__, __LINE__, __func__);
  }
  //  are of cross section
  Mc->give_area (eid,a);
  // computation of volume
  vol=a*l;
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->storeipvol (ipp,vol);
}



/**
   Function interpolates the nodal values to the integration points on the element
   
   @param eid    - element id
   @param nodval - nodal values
   @param ipval  - value at integration points
   
   @retval Function returns approximated values at integration points in the %vector ipval.
   
   JK, 3.10.2005
*/
void barel3d::intpointval (long eid,vector &nodval,vector &ipval)
{
  long ii,jj,i,k;
  double xi;
  vector w,gp;
  
  k=0;
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        //  value at integration point
        ipval[k] = approx (xi,nodval);
        k++;
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
void barel3d::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, ipval;
  vector w, gp, anv(ASTCKVEC(nne));
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
        reallocv (RSTCKVEC(intordsm[ii][jj],gp));
        reallocv (RSTCKVEC(intordsm[ii][jj],w));
        gauss_points (gp.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp[k];
          //  value in integration point
          ipval = approx (xi,anv);
          if ((ictn[i] & inistrain) && (j < Mm->ip[ipp].ncompstr))
          {
            Mm->ip[ipp].strain[idstra] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[i] & inistress) && (j < nstra + Mm->ip[ipp].ncompstr))
          {
            Mm->ip[ipp].stress[idstre] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & iniother) && (j < nstra+nstre+Mm->ip[ipp].ncompeqother))
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
   Function computes internal forces.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param s - direction %vector
   @param ifor - %vector of internal forces
   
   @retval Function returns nodal values of internal forces in %vector ifor.

   JK, 3.10.2005
*/
void barel3d::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq;
  iq=locstress;


  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
    

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor);
}



/**
   Function computes internal forces for nonlocal models


   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   @retval Function returns nodal values of nonlocal internal forces in %vector ifor.

   TKo, 7.2008
*/
void barel3d::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq=nonlocstress;
  
  //  computation of nonlocal stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor);
}



/**
   Function computes increment of internal forces.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   @retval Function returns nodal values of increments of internal forces in %vector ifor.

   JK, 29.4.2008
*/
void barel3d::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor);

}



/**
   Function computes nodal forces caused by temperature changes.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   
   @retval Function returns nodal values of forces caused by eigenstrains in %vector nfor.

   30.11.2002, JK
       7.2008, TKo
*/
void barel3d::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor)
{
  integratedquant iq;
  iq=eigstress;
  
  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor);
}



/**
   Function computes resulting internal forces. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting internal forces in parameter ifor.
   
   JK, 3.10.2005
*/
void barel3d::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  internal_forces (lcid,eid,0,0,ifor);

  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting increment of nonlocal internal forces.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting nonlocal internal forces in parameter ifor.
   
   TKo 7.2008
*/
void barel3d::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  nonloc_internal_forces (lcid,eid,0,0,ifor);

  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting increments of internal forces.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting increments of internal forces in parameter ifor.
   
   29.4.2008
*/
void barel3d::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  incr_internal_forces (lcid,eid,0,0,ifor);

  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting internal forces caused by eigenstrains.

   @param lcid - load case id
   @param eid  - element id
   @param nfor - %vector of nodal forces caused by eigenstrains

   @retval Function returns resulting increments of internal forces in parameter nfor.

   12.8.2001
*/
void barel3d::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  eigstrain_forces (lcid,eid,0,0,nfor);

  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   Function computes correct stresses at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 1.12.2006
*/
void barel3d::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipid,Mm->ip[ipid]);
  
}



/**
   Function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barel3d::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
  
}



/**
   Function computes nonlocal correct stresses at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 1.12.2006
*/
void barel3d::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->compnonloc_nlstresses (ipid);
  
}



/**
   Function computes correct stress increments at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barel3d::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstressesincr (ipid);
  
}



/**
   Function computes correct stresses caused by eigenstrains at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 1.12.2006
*/
void barel3d::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid,rmid;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  //  integration point number
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  index of relaxation model in tm array
  rmid = Mm->ip[ipid].grmid ();
  
  if (rmid>-1){
    //  eigenstresses are computed from a given formula, not from eigstrains
    Mm->eigenstresses (sig,ipid,rmid);
  }else{
    //
    //  eigenstresses are computed from eigenstrains \sigma_0 = D (-\eps_0)
    //
    // restore eigenstrains
    Mm->giveeigstrain (ipid,eigstr);
    // change sign of eigenstrain vector
    chsgnv(eigstr);
    //  matrix of stiffness of the material
    Mm->matstiff (d,ipid);
    // calculate eigenstresses    
    mxv (d,eigstr,sig);
  }
  
  Mm->storeeigstress (ipid,sig);
  
}

/**
   Function integrates selected quantity over the finite element.
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   
   @retval Function returns integrated values at nodes in %vector nv.

   JK, 1.12.2006
*/
void barel3d::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv)
{
  long ipp;
  double a,jac,w;
  vector ipv(ASTCKVEC(tncomp)),contr(ASTCKVEC(ndofe));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  //  coordinates of nodes
  Mt->give_node_coord3d (x,y,z,eid);
  //  id of integration point
  ipp=Mt->elements[eid].ipp[ri][ci];

  //  function assembles required quantity at integration point
  Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);
  //  strain-displacement (geometric) matrix
  geom_matrix (gm,x,y,z,jac);
  //  contribution to the nodal values
  mtxv (gm,ipv,contr);
  //  area of bar cross-section
  Mc->give_area (eid,a);
  w=2.0;
  cmulv (a*jac*w,contr);
  //  summation
  addv(contr,nv,nv);
}

/**
   Function searches the minimum and maximum strain in integration points
   
   @param min - the minimum strain
   @param max - the maximum strain
   @param lcid - load case id
   @param eid - element id
   
   @retval Function returns integrated values at nodes in %vector nv.

   JK, 16. 11. 2012
*/
void barel3d::find_extreme_strains (vector &min,vector &max,long lcid,long eid)
{
  long i,ipp;
  vector eps(ASTCKVEC(tncomp));
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[0][0];
  
  //  strains
  Mm->givestrain (lcid,ipp,eps);
  
  //  loop over the number of strain components
  for (i=0;i<tncomp;i++){
    if (min[i]>eps[i])
      min[i]=eps[i];
    if (max[i]<eps[i])
      max[i]=eps[i];
  }//  end of the loop over the number of strain components
  
}

/**
   Function searches the minimum and maximum stresses in integration points
   
   @param min - the minimum stress
   @param max - the maximum stress
   @param lcid - load case id
   @param eid - element id
   
   JK, 16. 11. 2012
*/
void barel3d::find_extreme_stresses (vector &min,vector &max,long lcid,long eid)
{
  long i,ipp;
  vector sig(ASTCKVEC(tncomp));
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[0][0];
  
  //  strains
  Mm->givestress (lcid,ipp,sig);
  
  //  loop over the number of stress components
  for (i=0;i<tncomp;i++){
    if (min[i]>sig[i])
      min[i]=sig[i];
    if (max[i]<sig[i])
      max[i]=sig[i];
  }//  end of the loop over the number of stress components
  
}
