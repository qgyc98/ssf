#include <math.h>
#include "barel2d.h"
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



barel2d::barel2d (void)
{
  long i;
  
  //  number nodes on element
  nne=2;
  //  number of DOFs on element
  ndofe=4;
  //  number of strain/stress components
  tncomp=1;
  //  number of functions approximated
  napfun=2;
  //  strain/stress state
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

  zero=Mp->zero;
}



barel2d::~barel2d (void)
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
   @param x,y - vectors of nodal coordinates at global coordinate system

   @retval Function returns direction %vector s.
   
   JK, 31.10.2006
*/
void barel2d::dirvect (vector &s,vector &x,vector &y)
{
  double d;
  
  s[0]=x[1]-x[0];
  s[1]=y[1]-y[0];
  
  d=sqrt(s[0]*s[0]+s[1]*s[1]);
  
  if (d<zero){
    print_err("zero norm of direction vector in function dirvect",__FILE__,__LINE__,__func__);
  }
  
  s[0]/=d;
  s[1]/=d;
}

/**
  Function returns coordinates of nodes at local element coordinate system.

  @param x  - %vector of nodal x-coordinates in global coordinate system
  @param y  - %vector of nodal y-coordinates in global coordinate system
  @param lx - %vector of nodal coordinates at local element coordinate system

  @retval Function returns %vector of nodal coordinates at parameter lx.
*/
void barel2d::giveloccoord(vector &x, vector &y,vector &lx)
{
  double l;

  lx[0] = 0.0;
  l = sqr(x[1] - x[0]) + sqr(y[1] - y[0]);
  l = sqrt(l);
  lx[1] = l;
}

/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param xi - natural coordinate
   
   JK, 25.8.2001
*/
void barel2d::bf_matrix (matrix &n,vector &s,double xi)
{
  vector bf(ASTCKVEC(nne));
  
  bf_lin_1d (bf.a,xi);
  
  n[0][0]=bf[0]*s[0];
  n[0][1]=bf[0]*s[1];
  n[0][2]=bf[1]*s[0];
  n[0][3]=bf[1]*s[1];

}

/**
   function assembles strain-displacement (geometric) %matrix
   geometric %matrix is in global coordinate system

   @param gm - geometric %matrix
   @param x,y - node coordinates in the global coordinate system
   @param jac - jacobian
   
   @retval Function returns geometric %matrix in parameter gm. 

   JK, 10.8.2001
*/
void barel2d::geom_matrix (matrix &gm,vector &x,vector &y,double &jac)
{
  double xi;
  vector s(ASTCKVEC(2)),lx(ASTCKVEC(nne)),dx(ASTCKVEC(nne));
  
  //  direction vector
  dirvect (s,x,y);
  //  local coordinates of nodes
  giveloccoord (x,y,lx);
  
  //  derivatives of the approximation functions
  //  the derivatives are independent on xi in this case
  xi=0.0;
  derivatives_1d (dx,jac,lx,xi);
  

  gm[0][0]=dx[0]*s[0];
  gm[0][1]=dx[0]*s[1];
  gm[0][2]=dx[1]*s[0];
  gm[0][3]=dx[1]*s[1];
}



/**
   Function approximates function defined by nodal values.

   @param xi - natural coordinate
   @param nodval - nodal values
   
   @retval Function returns approximated value for the given natural coordinate xi.

   JK
*/
double barel2d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));

  bf_lin_1d (bf.a,xi);

  scprd (bf,nodval,f);

  return f;
}



/**
   Function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l.
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   JK, 3.2.2002
*/
void barel2d::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i;
  fillm (0.0,tmat);
  
  for (i=0;i<tmat.m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<nodes.n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*2+0][i*2]=Mt->nodes[nodes[i]].e1[0];   tmat[i*2+0][i*2+1]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2]=Mt->nodes[nodes[i]].e1[1];   tmat[i*2+1][i*2+1]=Mt->nodes[nodes[i]].e2[1];
    }
  }

}


/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   T x_local = x_global
   x_local = T^T x_global
   
   direction %vector denotes %vector defining local coordinate axis
   
   @param x - %vector of nodal coordinates in element coordinate system
   @param tran - tranformation %matrix T
   @param gx,gy - vectors of nodal coordinates in global coordinate systems
   
   JK, 12.10.2008
*/
void barel2d::tran_mat (vector &x, matrix &tran,vector &gx,vector &gy)
{
  double c,l;
  
  fillm(0.0, tran);
  //  first direction vector is located in the first column
  tran[0][0]=gx[1]-gx[0];
  tran[1][0]=gy[1]-gy[0];

  //  length of the first direction vector
  l=sqrt(tran[0][0]*tran[0][0]+tran[1][0]*tran[1][0]);
  
  //  normalization of the first direction vector
  tran[0][0]/=l;
  tran[1][0]/=l;

  //  second direction vector is located in the second column
  tran[0][1]=0.0-tran[1][0];
  tran[1][1]=tran[0][0];
  
  //  transformation of coordinates of the first node
  x[0]=tran[0][0]*gx[0]+tran[1][0]*gy[0];
  c   =tran[0][1]*gx[0]+tran[1][1]*gy[0];

  //  transformation of coordinates of the second node
  x[1]=tran[0][0]*gx[1]+tran[1][0]*gy[1];
  c   =tran[0][1]*gx[1]+tran[1][1]*gy[1];
}



/**
   Function computes stiffness %matrix of twodimensional bar element.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   JK, 10.8.2001
*/
void barel2d::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long ipid;
  double a,jac,w;
  ivector nodes (ASTCKIVEC(nne));
  vector s(ASTCKVEC(2));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  fillm (0.0,sm);
  
  //  geometric matrix
  geom_matrix (gm,x,y,jac);
  //  number of first integration point on element
  ipid=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of material
  Mm->matstiff (d,ipid);
  //  area of bar cross-section
  Mc->give_area (eid,a);
  
  //  weight for one point integration
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

   JK, 10.8.2001
*/
void barel2d::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  // node coordinates
  Mt->give_node_coord2d (x,y,eid);

  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix from global coordinate system
  //  to the local nodal systems
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}



/**
   Function computes mass %matrix of twodimensional bar element.
   
   @param eid - element id
   @param mm - mass %matrix

   @retval Function returns mass %matrix at parameter mm.

   JK, 10.8.2001
*/
void barel2d::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double a,rho,jac,xi;
  ivector nodes(ASTCKIVEC(nne));
  vector lx(ASTCKVEC(nne)),s(ASTCKVEC(2));
  vector dens(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm));
  matrix n(ASTCKMAT(1,ndofe));
  
  fillm (0.0,mm);
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  only x coordinates are nonzero
  giveloccoord(x,y,lx);
  //  direction vector
  dirvect (s,x,y);
  
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
*/
void barel2d::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);  
  
  //  mass matrix
  mass_matrix (eid,mm,x,y);
  
  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}



/**
   Function computes load %matrix.
   
   @param eid - element id
   @param lm - load %matrix
   @param x - %vector of node coordinates
   
   @retval Function returns load %matrix in parameter lm.

   JK, 12.10.2008
*/
void barel2d::load_matrix (long eid,matrix &lm,vector &x)
{
  double l,a;
  
  fillm (0.0,lm);
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0]));
  if (l<Mp->zero){
    print_err("zero length of the %ld truss element",__FILE__, __LINE__, __func__, eid+1);
  }
  
  //  area of bar cross-section
  Mc->give_area (eid,a);
  
  lm[0][0]=a*l*1.0/3.0;
  lm[0][1]=0.0;
  lm[0][2]=a*l*1.0/6.0;
  lm[0][3]=0.0;
  
  lm[1][0]=0.0;
  lm[1][1]=a*l*1.0/3.0;
  lm[1][2]=0.0;
  lm[1][3]=a*l*1.0/6.0;
  
  lm[2][0]=a*l*1.0/6.0;
  lm[2][1]=0.0;
  lm[2][2]=a*l*1.0/3.0;
  lm[2][3]=0.0;

  lm[3][0]=0.0;
  lm[3][1]=a*l*1.0/6.0;
  lm[3][2]=0.0;
  lm[3][3]=a*l*1.0/3.0;
}



/**
   Function computes load %matrix of bar element with linear approximation.
   Load vector is obtained after premultiplying load %matrix by nodal load values.
   
   @param eid[in] - number of element
   @param lm[out] - load %matrix
   
   JK, 25.7.2001
*/
void barel2d::res_load_matrix (long eid,matrix &lm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  matrix tran(ASTCKMAT(2,2));
  
  //  node coordinates in the global coordinate system
  Mt->give_node_coord2d (gx,gy,eid);
  //  transformation matrix from global to element coordinate system
  tran_mat (x,tran,gx,gy);
  
  //  load matrix
  load_matrix (eid,lm,x);
  
  //  transformation of load matrix to global coordinate system
  lgmatrixtransfblock(lm, tran);
  
  //  transformation of load matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }

}






/**
   Function computes strain at integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barel2d::res_ip_strains (long lcid,long eid)
{
  long transf;
  vector x(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),aux(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat(ASTCKMAT(ndofe,ndofe));
  
  //  displacements at element nodes
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  ip_strains (lcid,eid,0,0,r);
  
}



/**
   Function computes strains at integration point.

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param r - nodal displacements
   
   10.5.2002
*/
void barel2d::ip_strains (long lcid,long eid,long ri,long ci,vector &r)
{
  long ipid;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),eps(ASTCKVEC(tncomp));
  matrix gm(ASTCKMAT(tncomp,ndofe));

  //  node cooridnates
  Mt->give_node_coord2d (x,y,eid);
  
  //  geometrix matrix
  geom_matrix (gm,x,y,jac);
  //  strain computation
  mxv (gm,r,eps);
  //  strain storage
  ipid=Mt->elements[eid].ipp[ri][ci];
  Mm->storestrain (lcid,ipid,eps);
}



/**
   Function computes strains in nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

   10.5.2002
*/
void barel2d::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipid;
  double l;
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
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      l = 0.5*Mt->give_length(eid);
      cmulv (l, eps);
      Mt->nodes[j].storestrain(lcid,0,l,eps);
    }
  }
}



/**
  The function computes nodal strains directly, averageing of nodal strain values is performed according to setup.
   
  @param lcid - load case id
  @param eid - element id
   
  Tomas Koudelka, 14.11.2013
*/
void barel2d::nod_strains_comp(long lcid,long eid)
{
  long i,j;
  double jac;
  ivector enod(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),eps(ASTCKVEC(tncomp));
  vector r(ASTCKVEC(ndofe)), aux;
  matrix tmat, gm(ASTCKMAT(tncomp,ndofe));

  //  displacements at element nodes
  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }

  //  node cooridnates
  Mt->give_node_coord2d (x,y,eid);
  //  node numbers
  Mt->give_elemnodes (eid, enod);
  
  //  geometrix matrix
  geom_matrix (gm,x,y,jac);
  //  strain computation
  mxv (gm,r,eps);
  //  storage of strains to the node
  for (i=0;i<nne;i++)
  {
    //  storage of strains to the node
    j=enod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      cmulv (jac,eps);
      Mt->nodes[j].storestrain(lcid,0,2.0*jac,eps);
    }
  }
}



/**
  The function computes nodal strains directly.
   
  @param lcid - load case id
  @param eid - element id
   
  Tomas Koudelka, 14.11.2013
*/
void barel2d::nod_strains(long lcid,long eid)
{
  long k,m;
  double jac;
  ivector enod(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), eps(ASTCKVEC(tncomp)), r(ASTCKVEC(ndofe)), aux;
  matrix tmat, gm(ASTCKMAT(tncomp,ndofe));

  //  displacements at element nodes
  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }

  //  node cooridnates
  Mt->give_node_coord2d (x,y,eid);
  //  node numbers
  Mt->give_elemnodes (eid, enod);
  
  //  geometrix matrix
  geom_matrix (gm,x,y,jac);
  //  strain computation
  mxv (gm,r,eps);
  //  storage of strains to the node
  m = lcid*tncomp;
  for (k=0;k<tncomp;k++)
  {
    Mt->nodes[enod[0]].strain[m+k] = eps(k);
    Mt->nodes[enod[1]].strain[m+k] = eps(k);
  }
}



/**
   Function computes strains at strain points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void barel2d::strains (long lcid,long eid,long ri,long ci)
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
   Function computes stress at integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barel2d::res_ip_stresses (long lcid,long eid)
{
  ip_stresses (lcid,eid,0,0);
}



/**
   Function computes stress at integration points.

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void barel2d::ip_stresses (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipid,Mm->ip[ipid]);
}



/**
   Function computes stress at integration points.

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void barel2d::ip_elast_stresses (long lcid,long eid,long ri,long ci)
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
void barel2d::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipid;
  double l;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[ri][ci];
  nodip_bar (ipid,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    
    //  storage of strains to the node
    j=nod[i];
    //  storage of strains to the node
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2)
    {
      l = 0.5*Mt->give_length(eid);
      cmulv (l, sig);
      Mt->nodes[j].storestress(lcid,0,l,sig);
    }
  }
}



/**
   Function computes stresses at required positions.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

*/
void barel2d::stresses (long lcid,long eid,long ri,long ci)
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
  The function computes other values at nodes of eid-th element.

  @param eid[in] - element id
   
  10.5.2002
*/
void barel2d::nod_other_ip (long eid,long ri,long ci)
{
  long i, j, ipid, ncompo;
  double l;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[ri][ci];
  nodip_bar (ipid,intordsm[0][0],ipnum);

  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Mm->givencompother (ipnum[i],0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    
    j=nod[i];
    //  storage of other values to the node
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother(0, ncompo, other);
    if (Mp->otheraver==2)
    {
      l = 0.5*Mt->give_length(eid);
      cmulv (l, other);
      Mt->nodes[j].storeother(0, ncompo, l, other);
    }
  }
}



/**
   Function computes internal forces.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   @retval Function returns nodal values of internal forces in %vector ifor.

   12.8.2001
*/
void barel2d::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq;
  iq=locstress;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  Mt->give_node_coord2d (x,y,eid);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);

}



/**
   Function computes nonlocal internal forces.


   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   @retval Function returns nodal values of nonlocal internal forces in %vector ifor.

   TKo, 7.2008
*/
void barel2d::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq=nonlocstress;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  Mt->give_node_coord2d (x,y,eid);
  
  //  computation of nonlocal stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
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
void barel2d::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  integratedquant iq;
  iq=stressincr;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  Mt->give_node_coord2d (x,y,eid);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);

}



/**
   Function computes nodal forces caused by temperature changes.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   @retval Function returns nodal values of forces caused by eigenstrains in %vector nfor.

   30.11.2002, JK
*/
void barel2d::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor)
{
  integratedquant iq;
  iq=eigstress;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);

  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
      
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor,x,y);
}



/**
   Function computes resulting internal forces. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting internal forces in parameter ifor.
   
   12.8.2001
*/
void barel2d::res_internal_forces (long lcid,long eid,vector &ifor)
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
void barel2d::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void barel2d::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void barel2d::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
   Function interpolates the nodal values to the integration points on the element
   quadratic approximation functions are used
   
   @param eid    - element id
   @param nodval - nodal values
   @param ipval  - value at integration points
   
   @retval Function returns approximated values at integration points in the %vector ipval.
*/
void barel2d::intpointval (long eid,vector &nodval,vector &ipval)
{
  long ii,jj,i,k;
  double xi;
  vector w,gp;

  k=0;
  for (ii = 0; ii < nb; ii++)
  {
    for (jj = 0; jj < nb; jj++)
    {
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ii][jj],gp));
      reallocv(RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      for (i = 0; i < intordsm[ii][jj]; i++)
      {
        xi=gp[i];
        //  value in integration point
        ipval[k] = approx (xi,nodval);
        k++;
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

   TKo, 7.2008
*/
void barel2d::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,ii;
  double xi;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordsm[ri][ci])), gp(ASTCKVEC(intordsm[ri][ci]));

  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++){
    for (j=0;j<intordsm[ri][ci];j++){
      xi=gp[i];
      if (ii==ipp){
	coord[0]=approx (xi,x);
	coord[1]=approx (xi,y);
	coord[2]=0.0;
      }
      ii++;
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
void barel2d::ipncoord (long eid, long ipp, vector &ncoord)
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
   
*/
void barel2d::ipvolume (long eid,long ri,long ci)
{
  long ipp;
  double a,l,vol;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

  Mt->give_node_coord2d (x,y,eid);
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  if (l<Mp->zero)
    print_err("zero length of the %ld truss element", __FILE__, __LINE__, __func__, eid+1);

  //  area of bar cross-section
  Mc->give_area (eid,a);
  // computation of volume
  vol=a*l;
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->storeipvol (ipp,vol);
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
void barel2d::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
   Function computes correct stresses at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 30.11.2006
*/
void barel2d::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipid,Mm->ip[ipid]);
  
}



/**
   Function computes correct stress increments at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 30.11.2006
*/
void barel2d::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstressesincr (ipid);
  
}



/**
   Function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barel2d::local_values (long /*lcid*/,long eid,long ri,long ci)
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
   
   JK, 30.11.2006
*/
void barel2d::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipid;
  
  ipid=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->compnonloc_nlstresses (ipid);
  
}



/**
   Function computes correct stresses caused by eigenstrains at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 30.11.2006
*/
void barel2d::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
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
   Function integrates selected quantity over the finite element
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case (used for certain quantity types)
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y - node coordinates
   
   @retval Function returns integrated values at nodes in %vector nv.

   JK, 30.11.2006
*/
void barel2d::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long ipid;
  double a,w,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector ipv(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  ipid=Mt->elements[eid].ipp[ri][ci];

  //  function assembles required quantity at integration point
  Mm->givequantity (iq,lcid,ipid,0,ipv);

  //  strain-displacement (geometric) matrix
  geom_matrix (gm,x,y,jac);
  
  //  contribution to the nodal values
  mtxv (gm,ipv,contr);
  
  //  area of bar cross-section
  Mc->give_area (eid,a);
  
  w=2.0;
  cmulv (w*a*jac,contr);

  //  summation
  addv(contr,nv,nv);
  
  ipid++;

}



/**
   Function approximates nonmechanical quantities in nodes of element.
   The nodal values are taken from the closest integration point.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything.
   
   12/06/2012 TKr
*/
void barel2d::mechq_nodval (long eid,vector &nodval,nontransquant qt)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[0][0];
  nodip_bar (ipid,intordsm[0][0],ipnum);

  for (i=0;i<nne;i++){

    //copy mechanical quantity from closest int. point
    nodval[i] = Mm->givemechq(qt, ipnum[i]);

  }
}



/**
   Function computes mechanical quantities in nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncnv - number of computed nodes on element (only first ncnv of nodes is calculated)
   @param nq - number of required mechanical quantities
   @param qt - array of types of required mechanical quantities
   
   @return The function does not return anything.
   
   Created by Tomas Koudelka, 29.11.2013
*/
void barel2d::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
{
  long i, j, ncompstr, ncompo, ncompnl, ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  intpoints ipb; // backup of the first integration point of element

  // id of the first integration point on element
  ipid=Mt->elements[eid].ipp[0][0];

  // element nodes
  Mt->give_elemnodes(eid, enod);

  //  numbers of integration points closest to element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodip_bar (ipid,intordsm[0][0],ipnum);

  // compute strains at nodes
  nod_strains(0, eid);
  
  // number of nonloc array components
  ncompnl = Mm->give_num_averq(ipid, Mm->givenonlocid(ipid));

  // store original content of the first integration point on element because 
  // it becomes working int. point for nodal values calculations on the given element
  ipb.copy(Mm->ip[ipid], Mb->nlc, ncompnl, 1);

  // The first integration point will be used for computation of nodal values temporarily
  // then the original content of the first integration point will be restored
  for (i=0;i<ncnv;i++)
  {
    // number of strain components
    ncompstr = Mm->ip[ipnum[i]].ncompstr;
    // take nodal strain and store them to the first (working) integration point
    Mm->storestrain(0, ipid, 0, ncompstr, Mt->nodes[enod[i]].strain);

    ncompo = Mm->ip[ipnum[i]].ncompeqother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point the eqother values may be rewritten -> take them from the backup
        Mm->storeeqother(ipid, 0, ncompo, ipb.eqother);
      else
        Mm->storeeqother(ipid, 0, ncompo, Mm->ip[ipnum[i]].eqother);
      // take other values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point the, other values may be rewritten -> take them from the backup
        Mm->storeother(ipid, 0, ncompo, ipb.other);
      else
        Mm->storeother(ipid, 0, ncompo, Mm->ip[ipnum[i]].other);
    }

    if (ncompnl){
      // take nonloc values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point the nonloc values may be rewritten -> take them from the backup
        Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc);
      else
        Mm->storenonloc(ipid, 0, ncompnl, Mm->ip[ipnum[i]].nonloc);
    }

    // compute nodal stresses and internal variables of the material at node
    switch(Mp->matmodel)
    {
      case local:
        Mm->computenlstresses(ipid,Mm->ip[ipid]);
        break;
      case nonlocal:
        Mm->compnonloc_nlstresses (ipid);
        break;
      default:
        print_err("unknown approach in material model is required (Mp->matmodel=%ld)", __FILE__, __LINE__, __func__, Mp->matmodel);
    }
   
    // the internal material variables are stored in ip[ipid].other array and
    // they must be send to ip[ipid].eqother array
    Mm->updateipvalmat(ipid, 0, 0);
    
    //give calculated mechanical quantity from the first int. point
    for (j=0; j<nq; j++)
      nodval[j*ncnv+i] = Mm->givemechq(qt[j], ipid);
  }

  // restore original integration point content of strain/stress/other/eqother arrays
  Mm->storestrain(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.strain);
  Mm->storestress(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.stress);
  Mm->storeother(ipid, 0, ipb.ncompother, ipb.other);
  Mm->storeeqother(ipid, 0, ipb.ncompeqother, ipb.eqother);
  Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc);
}
