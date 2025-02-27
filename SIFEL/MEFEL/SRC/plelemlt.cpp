#include "plelemlt.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "adaptivity.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"
#include "gadaptivity.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>



planeelemlt::planeelemlt (void)
{
  long i,j;

  //  number nodes on element  
  nne=3;
  //  number of DOFs on element
  ndofe=6;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=2;
  //  order of numerical integration of mass matrix
  intordmm=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;
  
  //  number of blocks (parts of geometric matrix)
  nb=1;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=3;

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
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=1;
}

planeelemlt::~planeelemlt (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;
  
  delete [] ncomp;
  delete [] cncomp;
}



/**
   function approximates function defined by nodal values
   area coordinates are used
   
   @param areacoord - vector containing area coordinates
   @param nodval - nodal values
   
   JK
*/
double planeelemlt::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function approximates function defined by nodal values
   natural coordinates are used
   
   @param xi,eta - natural coordinates
   @param nodval - nodal values
   
   JK
*/
double planeelemlt::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector areacoord(3);

  //  conversion of natural coordinates to area coordinates
  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];

  scprd (areacoord,nodval,f);
  return f;
}

/**
   function returns %matrix of base function

   @param n - %matrix of approximation functions
   @param xi,eta - natural coordinates
   
   JK, 17.8.2001
*/
void planeelemlt::bf_matrix (matrix &n,double xi,double eta)
{
  vector bf(nne);

  bf_lin_3_2d (bf.a,xi,eta);

  fillm (0.0,n);

  n[0][0]=bf[0];
  n[0][2]=bf[1];
  n[0][4]=bf[2];
  
  n[1][1]=bf[0];
  n[1][3]=bf[1];
  n[1][5]=bf[2];
}

/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param x,y - node coordinates
   
   JK, 17.8.2001
*/
void planeelemlt::geom_matrix (matrix &gm,vector &x,vector &y)
{
  double det;
  vector b(3),c(3);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  fillm (0.0,gm);

  gm[0][0]=b[0];
  gm[0][2]=b[1];
  gm[0][4]=b[2];
  
  gm[1][1]=c[0];
  gm[1][3]=c[1];
  gm[1][5]=c[2];
  
  gm[2][0]=c[0];
  gm[2][1]=b[0];
  gm[2][2]=c[1];
  gm[2][3]=b[1];
  gm[2][4]=c[2];
  gm[2][5]=b[2];
}

/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param inodes - array containing node numbers
   @param tmat - transfomation %matrix
   
   JK, 17.8.2001
*/
void planeelemlt::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*2][i*2]   = Mt->nodes[nodes[i]].e1[0];  tmat[i*2][i*2+1]   = Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2] = Mt->nodes[nodes[i]].e1[1];  tmat[i*2+1][i*2+1] = Mt->nodes[nodes[i]].e2[1];
    }
  }
}

/**
   function computes stiffness %matrix of plane stress triangular
   finite element with linear approximation functions

   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - node coordinates
   
   JK, 17.8.2001
*/
void planeelemlt::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long ipp;
  double xi,eta,jac,det,thick;
  ivector nodes(nne);
  vector t(nne);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);

  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  xi=1.0/3.0;
  eta=1.0/3.0;

  fillm (0.0,sm);

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  // geometric matrix
  geom_matrix (gm,x,y);
  
  //  stiffness matrix of material
  Mm->matstiff (d,ipp);
  
  thick = approx_nat (xi,eta,t);
  
  //  det is equal to double area of the element
  jac=thick*det/2.0;
  
  //  contribution to the stiffness matrix of the element
  bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);

}


/**
   function computes B^T D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - B^T D %matrix

   53/10/2019 TKr according to TKo lintet.cpp
*/
void planeelemlt::bd_matrix (long eid,long ri,long ci,matrix &bd)
{
  long ipp;
  double xi,eta,jac,det,thick;
  ivector nodes(nne);
  vector t(nne);
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  xi=1.0/3.0;
  eta=1.0/3.0;

  fillm (0.0,bd);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  geometric matrix
  geom_matrix (gm,x,y);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  thick = approx_nat (xi,eta,t);
  
  //  det is equal to double area of the element
  jac=thick*det/2.0;
  
  //  B^T D
  mtxm (gm,d,bd);
  //  B^T D Jac

  cmulm (jac,bd,bd);
}

/**
   function computes integral of D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - integral of D %matrix

   25/10/2019 TKr according to TKo lintet.cpp
*/
void planeelemlt::dd_matrix (long eid,long ri,long ci,matrix &dd)
{
  long ipp;
  double xi,eta,jac,det,thick;
  ivector nodes(nne);
  vector t(nne);
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp));


  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  xi=1.0/3.0;
  eta=1.0/3.0;

  fillm (0.0,dd);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  thick = approx_nat (xi,eta,t);
  
  //  det is equal to double area of the element
  jac=thick*det/2.0;
  
  //  D Jac
  cmulm (jac,d,dd);
}


/**
   function computes stiffness %matrix of plane stress triangular
   finite element with linear approximation functions

   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 17.8.2001
*/
void planeelemlt::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);

  Mt->give_node_coord2d (x,y,eid);

  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}


/**
   function computes mass %matrix of the plane stress triangular
   finite element with linear approximation functions
   
   @param eid - number of element
   @param mm - mass %matrix
   @param x,y - node coordinates
   
   JK, 17.6.2001
*/
void planeelemlt::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double jac,det,thick,rho;
  ivector nodes(nne);
  vector w(intordmm),gp1(intordmm),gp2(intordmm),t(nne),dens(nne);
  matrix n(napfun,ndofe);
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  //  density of the material
  Mc->give_density (eid,nodes,dens);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    bf_matrix (n,gp1[i],gp2[i]);
    
    thick = approx_nat (gp1[i],gp2[i],t);
    rho = approx_nat (gp1[i],gp2[i],dens);
    
    jac=w[i]*thick*rho*det;
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function computes mass %matrix of the plane stress triangular
   finite element with linear approximation functions
   
   @param eid - number of element
   @param mm - mass %matrix
   
   JK, 17.6.2001
*/
void planeelemlt::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);

  Mt->give_node_coord2d (x,y,eid);

  mass_matrix (eid,mm,x,y);
  
  if (Mp->diagmass==1){
    diagonalization (mm);
  }

  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }

}


/**
   function computes load %matrix of the plane stress triangular
   finite element with linear approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   @param x,y - node coordinates

   JK, 25.7.2001
*/
void planeelemlt::load_matrix (long eid,matrix &lm,vector &x,vector &y)
{
  long i;
  double jac,det,thick;
  ivector nodes(nne);
  vector w(intordmm),gp1(intordmm),gp2(intordmm),b(3),c(3),t(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    bf_matrix (n,gp1[i],gp2[i]);
    
    thick = approx_nat (gp1[i],gp2[i],t);
    
    //  zkontrolovat deleni dvema
    jac=w[i]*thick*det;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function computes load %matrix of the plane stress triangular
   finite element with linear approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix

   JK, 25.7.2001
*/
void planeelemlt::res_load_matrix (long eid,matrix &lm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);

  load_matrix (eid,lm,x,y);

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
   
   JK
*/
void planeelemlt::res_ip_strains (long lcid,long eid)
{
  vector aux,x(nne),y(nne),r(ndofe);
  ivector nodes(nne);
  matrix tmat;

  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
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

  ip_strains (lcid,eid,0,0,x,y,r);
}

/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y - node coordinates
   @param r - nodal displacements
   
   JK, 10.5.2002
*/
void planeelemlt::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long ipp;
  vector eps(tncomp);
  matrix gm(tncomp,ndofe);

  geom_matrix (gm,x,y);
  mxv (gm,r,eps);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->storestrain (lcid,ipp,cncomp[0],ncomp[0],eps);
}

/**
   function computes strains at nodes of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   JK, 10.5.2002
*/
void planeelemlt::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(gncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelt (ipp,intordsm[0][0],ipnum);
  
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
   
   JK, 25.9.2004
*/
void planeelemlt::nod_strains_comp (long lcid,long eid)
{
  long i;
  ivector nodes(nne);
  vector x(nne),y(nne),r(ndofe),eps(tncomp),aux;
  matrix tmat,gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
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
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y);
    //  strain computation
    mxv (gm,r,eps);
    
    //  storage of strains to the node
    Mt->nodes[i].storestrain (lcid,0,eps);
  }

}



void planeelemlt::strains (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{

}




/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void planeelemlt::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval;
  
  k=0;
  reallocv (nne,nodval);
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx_nat (xi,eta,nodval);
    k++;
  }
}






/**
   function computes stresses at integration points of element
   stresses are computed by material models
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 10.5.2002
*/
void planeelemlt::res_ip_stresses (long /*lcid*/,long eid)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[0][0];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
}


/**
   function computes stresses at integration points of element
   stresses are computed from strains with the help of elastic stiffness
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   JK, 27.11.2006
*/
void planeelemlt::ip_elast_stresses (long lcid,long eid,long ri,long ci)
{
  long ipp;
  vector eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  ipp=Mt->elements[eid].ipp[ri][ci];

  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  //  strains
  Mm->givestrain (lcid,ipp,eps);
  
  //  elastic stresses
  mxv (d,eps,sig);
  
  Mm->storestress (lcid,ipp,sig);
}


/**
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 10.5.2002
*/
void planeelemlt::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(gncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelt (ipp,intordsm[0][0],ipnum);
  
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



void planeelemlt::stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
}



void planeelemlt::nod_others (long eid,long ri,long ci)
{
  long ipp, i, ncomp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),other,aux,natcoord(2);
  ivector nodes(nne);
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planelt (nxi,neta);

  Mt->give_elemnodes (eid,nodes);
  ipp=Mt->elements[eid].ipp[ri][ci];
  ncomp = Mm->ip[ipp].ncompother;

  reallocv (ncomp,other);
  lhs = new double [ncomp*3];
  rhs = new double [ncomp*3];
  lsm = new double [9];
  
  nullv (lsm,9);
  nullv (rhs,ncomp*3);
  
  for (i = 0;i < ncomp; i++)
    other[i] = Mm->ip[ipp].eqother[i];

  natcoord[0]=1.0/3.0;  natcoord[1]=1.0/3.0;
  matassem_lsm (lsm,natcoord);
  rhsassem_lsm (rhs,natcoord,other);
  
  solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp);
  Mt->other_nodal_values (nodes,nxi,neta,nxi,lhs,2,0,ncomp);
  
  delete [] lsm;  delete [] lhs;  delete [] rhs;
}



/**
   function computes internal forces

   this function is used in plane stress/strain elements (function is called
   by function res_internal_forces) and shell elements

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - node coordinates

   JK, 27.11.2006
*/
void planeelemlt::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);

}



/**
   function computes internal forces for nonlocal models

   this function is used in plane stress/strain elements (function is called
   by function res_nonloc_internal_forces) and shell elements

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y - node coordinates
   
   JK, 27.11.2006
*/
void planeelemlt::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of nonlocal stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes increments of internal forces

   this function is used in plane stress/strain elements (function is called
   by function res_internal_forces) and shell elements

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - node coordinates

   JK, 27.11.2006
*/
void planeelemlt::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=stressincr;
  
  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);

}



/**
   function computes nodal forces caused by temperature changes
   
   this function is used in plane stress/strain elements (function is called
   by function res_eigstrain_forces) and shell elements

   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   30.11.2002, JK
*/
void planeelemlt::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
{
  integratedquant iq;
  iq=eigstress;
  
  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor,x,y);
}



/**
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - internal forces
   
   JK, 27.11.2006
*/
void planeelemlt::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  
  internal_forces (lcid,eid,0,0,ifor,x,y);

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
   function computes nonlocal internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   JK, 27.11.2006
*/
void planeelemlt::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);

  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y);

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
   function computes resulting increments of internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - internal forces
   
   TKo 7.2008
*/
void planeelemlt::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  
  incr_internal_forces (lcid,eid,0,0,ifor,x,y);

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
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   TKr, 05/06/2014
*/
void planeelemlt::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  vector areacoord(3);
  
  fillv (1.0/3.0,areacoord);
  
  ipval[0]=approx (areacoord,nodval);
}


/**
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - array containing nodal forces
   
   JK, 27.11.2006
*/
void planeelemlt::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);

  Mt->give_node_coord2d (x,y,eid);

  eigstrain_forces (lcid,eid,0,0,nfor,x,y);

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
   
   TKo, 7.2008
*/
void planeelemlt::compute_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemlt::local_values (long /*lcid*/,long eid,long ri,long ci)
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
   
   JK, 27.11.2006
*/
void planeelemlt::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->compnonloc_nlstresses (ipp);
}



/**
   function computes correct stresses increment at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemlt::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstressesincr (ipp);
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemlt::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
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
   @param x,y - node coordinates
   
   JK, 27.11.2006
*/
void planeelemlt::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long ipp;
  double xi,eta,det,thick;
  ivector nodes(nne);
  vector t(nne),ipv(tncomp),contr(ndofe);
  matrix gm(tncomp,ndofe);
  
  Mc->give_thickness (eid,nodes,t);
  
  fillv (0.0,nv);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  xi=1.0/3.0;
  eta=1.0/3.0;

  thick = approx_nat (xi,eta,t);
  
  
  //  function assembles required quantity at integration point
  Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);
  //  strain-displacement (geometric) matrix
  geom_matrix (gm,x,y);
  
  //  contribution to the nodal values
  mtxv (gm,ipv,contr);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  cmulv (det*thick/2.0,contr);
  
  //  summation
  addv(contr,nv,nv);
}



/**
   The function integrates arbitrary selected quantity over the finite element, e.g.
   it performs \int_{\Omega} \mbf{\sigma} d\Omega which results in integrated values that can 
   be used in the homogenization problems.

   
   @param eid - element id
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param iv - integrated values (output)
   
   @return The function returns nodal values calculated in the %vector iv
   
   25/10/2019 TKr according to TKo
*/
void planeelemlt::elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv)
{
  long ipp;
  double xi,eta,jac,det,thick;
  ivector nodes(nne);
  vector t(nne);
  vector ipv(ASTCKVEC(iv.n));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);

  ipp=Mt->elements[eid].ipp[0][0];
  Mt->give_node_coord2d (x,y,eid);

  nullv(iv);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  xi=1.0/3.0;
  eta=1.0/3.0;
  
  thick = approx_nat (xi,eta,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;

  //  function assembles required quantity at integration point
  Mm->givequantity (iq, lcid, ipp, cncomp[0], ipv);

  cmulv (jac, ipv, iv);
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
void planeelemlt::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  vector x(nne),y(nne),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    if (ii==ipp){
      coord[0]=approx_nat (gp1[i],gp2[i],x);
      coord[1]=approx_nat (gp1[i],gp2[i],y);
      coord[2]=0.0;
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
void planeelemlt::ipncoord (long eid,long ipp,vector &ncoord)
{
  long i, ii, ri, ci;
  vector w, gp1, gp2;
  
  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      if (intordsm[ri][ci] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp2));
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];
  
      for (i=0;i<intordsm[ri][ci];i++){
        if (ii==ipp){
          ncoord[0]=gp1[i];
          ncoord[1]=gp2[i];
          ncoord[2]=0.0;
        }
        ii++;
      }
    }
  }
}



/**
   function returns coordinates of integration points

   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param coord - array of coordinates

   19.1.2002
*/
void planeelemlt::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  vector x(nne),y(nne),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  
  for (long i=0;i<intordsm[ri][ci];i++){
    coord[i][0]=approx_nat (gp1[i],gp2[i],x);
    coord[i][1]=approx_nat (gp1[i],gp2[i],y);
    coord[i][2]=0.0;
  }
}

/**
   function computes nodal forces generated by continuous edge load
   
   @param eid - element id
   @param le - loaded edge indicator
   @param nv - array of nodal values
   @param nf - vector of nodal forces
   
   JK, 5.4.2011
*/
void planeelemlt::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double xi,eta,jac;
  vector x(nne),y(nne),gp(intordb),w(intordb),av(ndofe),v(ndofe),tnv(ndofe);
  matrix n(napfun,ndofe),am(ndofe,ndofe);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);
  
  if (le[0]>0){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=(1.0-gp[i])/2.0;  eta=(1.0+gp[i])/2.0;
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,gp[i],0);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[0]=nv[0];
    av[1]=nv[1];
    av[2]=nv[2];
    av[3]=nv[3];
    
    if (le[0]==2){
      locglob_nodeval (0,av,tnv,x,y);
      fillv(0.0,av);
      
      av[0]=tnv[0];
      av[1]=tnv[1];
      av[2]=tnv[2];
      av[3]=tnv[3];
    }
    
    mxv (am,av,v);  addv (nf,v,nf);
  }

  if (le[1]>0){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=0.0;  eta=(1.0-gp[i])/2.0;
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,gp[i],1);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[2]=nv[4];
    av[3]=nv[5];
    av[4]=nv[6];
    av[5]=nv[7];
    
    if (le[1]==2){
      av[0]=nv[4];
      av[1]=nv[5];
      av[2]=nv[6];
      av[3]=nv[7];
      
      locglob_nodeval (1,av,tnv,x,y);
      fillv(0.0,av);
      
      av[2]=tnv[0];
      av[3]=tnv[1];
      av[4]=tnv[2];
      av[5]=tnv[3];
    }
    mxv (am,av,v);  addv (nf,v,nf);
  }

  if (le[2]>0 ){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=(1.0+gp[i])/2.0;  eta=0.0;
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,gp[i],2);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[4]=nv[8];
    av[5]=nv[9];
    av[0]=nv[10];
    av[1]=nv[11];
    
    if (le[2]==2){
      av[0]=nv[8];
      av[1]=nv[9];
      av[2]=nv[10];
      av[3]=nv[11];
      
      locglob_nodeval (2,av,tnv,x,y);
      fillv(0.0,av);
      
      av[4]=tnv[0];
      av[5]=tnv[1];
      av[0]=tnv[2];
      av[1]=tnv[3];
    }
    
    mxv (am,av,v);  addv (nf,v,nf);
  }
}

/**
   function transforms nodal values from the local system to the global system
   
   @param edid - edge id
   @param nv - array of nodal values
   @param tnv - array of transformed nodal values
   @param x,y - arrays of nodal coordinates
   
   JK, 5.4.2011
*/
void planeelemlt::locglob_nodeval (long edid,vector &nv,vector &tnv,vector &x,vector &y)
{
  double norm;
  vector g1(2),g2(2),lv(2),gv(2);
  matrix t(2,2);
  
  if (edid==0){
    //  first edge of the triangle
    //  edge between nodes 1 and 2
    
    //  vector parallel with the edge
    g1[0]=x[1]-x[0];
    g1[1]=y[1]-y[0];
    
    norm=sqrt(g1[0]*g1[0]+g1[1]*g1[1]);
    if (norm<Mp->zero){
      print_err("negative norm of direction vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit vector parallel with the edge
    g1[0]/=norm;
    g1[1]/=norm;
    
    //  vector parallel with the edge 3-1
    g2[0]=x[0]-x[2];
    g2[1]=y[0]-y[2];
    
    //  orthogonalization of the g2 with respect to g1
    norm=1.0/(g1[0]*g2[0]+g1[1]*g2[1]);
    g2[0]= g1[0]-norm*g2[0];
    g2[1]= g1[1]-norm*g2[1];
    
    norm=sqrt(g2[0]*g2[0]+g2[1]*g2[1]);
    if (norm<Mp->zero){
      print_err("negative norm of normal vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit normal vector to the edge
    g2[0]/=norm;
    g2[1]/=norm;
    
  }
  if (edid==1){
    //  second edge of the triangle
    //  edge between nodes 2 and 3
    
    //  vector parallel with the edge
    g1[0]=x[2]-x[1];
    g1[1]=y[2]-y[1];
    
    norm=sqrt(g1[0]*g1[0]+g1[1]*g1[1]);
    if (norm<Mp->zero){
      print_err("negative norm of direction vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit vector parallel with the edge
    g1[0]/=norm;
    g1[1]/=norm;
    
    //  vector parallel with the edge 1-2
    g2[0]=x[1]-x[0];
    g2[1]=y[1]-y[0];
    
    //  orthogonalization of the g2 with respect to g1
    norm=1.0/(g1[0]*g2[0]+g1[1]*g2[1]);
    g2[0]= g1[0]-norm*g2[0];
    g2[1]= g1[1]-norm*g2[1];
    
    norm=sqrt(g2[0]*g2[0]+g2[1]*g2[1]);
    if (norm<Mp->zero){
      print_err("negative norm of normal vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit normal vector to the edge
    g2[0]/=norm;
    g2[1]/=norm;
    
  }
  if (edid==2){
    //  third edge of the triangle
    //  edge between nodes 3 and 1
    
    //  vector parallel with the edge
    g1[0]=x[0]-x[2];
    g1[1]=y[0]-y[2];
    
    norm=sqrt(g1[0]*g1[0]+g1[1]*g1[1]);
    if (norm<Mp->zero){
      print_err("negative norm of direction vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit vector parallel with the edge
    g1[0]/=norm;
    g1[1]/=norm;
    
    //  vector parallel with the edge 3-2
    g2[0]=x[2]-x[1];
    g2[1]=y[2]-y[1];
    
    //  orthogonalization of the g2 with respect to g1
    norm=1.0/(g1[0]*g2[0]+g1[1]*g2[1]);
    g2[0]= g1[0]-norm*g2[0];
    g2[1]= g1[1]-norm*g2[1];
    
    norm=sqrt(g2[0]*g2[0]+g2[1]*g2[1]);
    if (norm<Mp->zero){
      print_err("negative norm of normal vector", __FILE__, __LINE__, __func__);
    }
    
    //  unit normal vector to the edge
    g2[0]/=norm;
    g2[1]/=norm;
    
  }

  t[0][0]=g1[0];
  t[1][0]=g1[1];

  t[0][1]=g2[0];
  t[1][1]=g2[1];
  
  mxv (t.a,nv.a+0,tnv.a+0,2,2);
  mxv (t.a,nv.a+2,tnv.a+2,2,2);
  mxv (t.a,nv.a+4,tnv.a+4,2,2);
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
void planeelemlt::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, ipval;
  vector w, gp1, gp2, anv(nne);
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
        reallocv (intordsm[ii][jj],w);
        gauss_points_tr (gp1.a, gp2.a, w.a, intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp1[k];
          eta=gp2[k];
          //  value in integration point
          ipval = approx_nat (xi,eta,anv);
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
   07.2008 TKo - multiplictaion by thickness added
*/
void planeelemlt::ipvolume (long eid,long ri,long ci)
{
  long ipp;
  double xi,eta,jac,thick;
  ivector nodes(nne);
  vector t(nne);
  vector x(nne),y(nne);
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);
  ipp=Mt->elements[eid].ipp[ri][ci];
  xi=1.0/3.0;  eta=1.0/3.0;
  
  thick = approx_nat (xi,eta,t);

  jac_2d (jac,x,y,xi,eta);
  jac *= thick/2.0;
  
  Mm->storeipvol (ipp,jac);
}


////////////////////       /* termitovo */       ////////////////////////////////////

/**
   Function integrates function "NT*D*B*r" (= NT*Stress) over whole element.
   N is matrix of basic functions.
   !!! Values in 'ntdbr' are stored unusually. Not after this manner {(val[0]; ... ;val[tncomp])[0] ; ... ; (val[0]; ... ;val[tncomp])[nne]}
   but after this manner {(val[0]; ... ;val[nne])[0] ; ... ; (val[0]; ... ;val[nne])[ntncomp]}
   
   @param eid   - element id
   @param ntdbr - empty(returned) array, dimension is tncomp*nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemlt::ntdbr_vector (long eid,vector &ntdbr)
{
  long i,j,ipp,ri,ci;
  double det,thick,jac,w;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),t(nne),eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  fillv (1.0/3.0,areacoord);
  w = 0.5;
  
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  thick = approx (areacoord,t);
  jac = w*thick*det;
  
  Mm->matstiff (d,ipp);
  Mm->givestrain (0,ipp,eps);
  mxv (d,eps,sig);
  
  for (i=0;i<tncomp;i++)
    for (j=0;j<nne;j++)
      ntdbr[j+i*nne] = jac * sig[i] * areacoord[j];
}

/**
   Function integrates function "NT*N" over whole element.
   N is matrix of basic functions.
   !!! Matrix N is composed of three same blocks, it is computed only one third.
   
   @param eid - element id
   @param ntn - empty(returned) 2Darray, dimension is nne x nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemlt::ntn_matrix (long eid,matrix &ntn)
{
  long i;
  double det,thick,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),t(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  fillv (1.0/3.0,areacoord);

  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  thick = approx (areacoord,t);
  jac = thick*det/24.0;

  fillm (jac,ntn);

  for (i=0;i<nne;i++)
    ntn[i][i] *= 2;
}

/**
   Function
   1. computes "e2" - square of energy norm of error of solution at element
      e2 = integral{ (der_fine - der)T * D * (der_fine - der) }
      der_fine = recovered strains or stresses, values in nodes are defined by "rderfull" and interpolated by base functions over element
      der = strains/stresses obtained by FEM
      D_inv = inverse stiffness matrix of material
   2. computes "u2" - square of energy norm of strain at element
      u2 = integral{ epsT * D * eps }
      eps = strain obtained by FEM
      D = stiffness matrix of material
   3. computes area of element and returns
   4. computes "sizel" (characteristic size of element)
   
   @param eid - element id
   @param e2 - empty(returned) value; 
   @param u2 - empty(returned) value; 
   @param sizel - empty(returned) value; 
   @param rderfull - 1Darray of strains/stresses in nodes; dimmension is gt->nn x ncomp
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
double planeelemlt :: compute_error (long eid, double &e2, double &u2, double &sizel, vector *rderfull,long flags)
{
  long i, ipp;
  double det, contr;
  ivector nodes(nne);
  vector x(nne), y(nne), areacoord(3), t(nne);
  vector der1(tncomp), der2(tncomp), der_fine(tncomp), der_err(tncomp);
  matrix d(tncomp,tncomp), dinv(tncomp,tncomp);
  
  Mt->give_elemnodes (eid, nodes);
  Mt->give_node_coord2d (x, y, eid);
  Mc->give_thickness (eid, nodes, t);
  
  ipp = Mt->elements[eid].ipp[0][0];
  
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  // stiffness matrix
  Mm->matstiff (d,ipp);
  if ( flags & 2 )  invm (d, dinv, 1.0e-20);
  
  // compute u2
  fillv (1.0/3.0, areacoord);
  
  if (!(flags & 2)){          // strain - der1 == strain
    Mm->givestrain (0, ipp, der1);
    mxv (d, der1, der2);
  }
  else                        // stress - der1 == stress
    if (!(flags & 4)){          // own
      Mm->givestrain (0,ipp,der2);
      mxv (d, der2, der1);
    }
    else{                       // other
      Mm->givestress (0,ipp,der1);
      mxv (dinv,der1,der2);
    }
  
  scprd (der1,der2,contr);
  u2 = contr * 0.5*det*approx (areacoord,t);
  
  // compute e2
  long intord = 3; // =intordcm
  vector gp1(intord), gp2(intord), w(intord);
  
  gauss_points_tr (gp1.a, gp2.a, w.a, intord);
  
  e2 = 0;
  for (i=0; i<intord; i++) {
    areacoord[0] = gp1[i];
    areacoord[1] = gp2[i];
    areacoord[2] = 1.0-gp1[i]-gp2[i];
    
    // der above
    give_der_star (areacoord, rderfull, nodes, der_fine, Mt->nn);
    subv (der_fine, der1, der_err);
    
    if (!(flags & 2))          // strain - der1 == strain
      vxmxv (der_err,d,contr);
    else                       // stress - der1 == stress
      vxmxv (der_err,dinv,contr);
    
    e2 += contr * w[i]*det*approx(areacoord,t);
  }
  
  sizel = sqrt(1.1547005384*det);
  return  det/2.0;
}







/**
   This function prepare array 'spder' for using in function `least_square::spr_default`.
   It allocates and fills up array 'spder' by derivatives(strain ro stress) in sampling points.
   For planeelemlt sampling points == main integration point == Gauss's point at the centroid of element.
   
   @param eid   - element id
   @param spder - allocated array 
   @param flags - determines by which derivate is spder filled
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemlt::elchar (long eid,double *&spder,long flags)
{
  long ipp,ri,ci;
  vector eps;
  matrix d;
  
  spder = new double[tncomp*1];
  
  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  
  if (!(flags & 2)){
    eps.n = tncomp;
    eps.a = spder;
    Mm->givestrain (0,ipp,eps);
    eps.a = NULL;
  }
  else if (!(flags & 4)){
    reallocv (tncomp,eps);
    reallocm (tncomp,tncomp,d);
    Mm->matstiff (d,ipp);
    Mm->givestrain (0,ipp,eps);
    mxv (d.a,eps.a,spder,d.m,d.n);
  }
  else{
    eps.n = tncomp;
    eps.a = spder;
    Mm->givestress (0,ipp,eps);
    eps.a = NULL;
  }
}

/**
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
double planeelemlt::error (long eid,vector &n,double &a)
{
  long i;
  double answer;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),nodval(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  
  fillv (1.0/3.0,areacoord);
  
  a = 0.5 * (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  for (i=0;i<nne;i++)
    nodval[i] = n[nodes[i]];
  
  answer = approx (areacoord,nodval);
  answer = a*answer*answer;
  
  return answer;
}
////////////////////       /* termitovo */       ////////////////////////////////////
