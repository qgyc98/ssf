#include "axisymqt.h"
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
#include "loadcase.h"
#include "intpoints.h"
#include <stdlib.h>

axisymqt::axisymqt (void)
{
  long i,j;
  
  //  number nodes on element
  nne=6;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=4;
  //  number of approximated functions
  napfun=2;
  //  number of edges
  ned=3;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration of mass matrix
  intordmm=3;
  //  order of numerical integration on edges
  intordb=2;
  //  strain/stress state
  ssst=axisymm;
  
  
  
  //  number of blocks
  //  implementation enables one block computation as well as three blocks
  //  three blocks are useful due to different integration orders
  //nb=3;
  nb=1;
  
  //  number of components in blocks
  ncomp = new long [nb];
  //  cumulative number of components in blocks
  cncomp = new long [nb];
  
  //  number of integration points in blocks
  //  order of numerical integration in blocks
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }

  
  if (nb==1){
    //
    //  one block is used, it means, geometric matric is not divided
    //
    //  number of components in block
    ncomp[0]=4;
    
    //  cumulative number of components in blocks
    cncomp[0]=0;
    
    nip[0][0]=3;
    intordsm[0][0]=3;
  }
  
  if (nb==3){
    //  three blocks are used, it means, geometric matrix is divided into
    //  three blocks
    //  first block contains epsilon_x and epsilon_y
    //  second block contains epsilon_fi
    //  third block contains epsilon_xy
    
    //  number of components in blocks
    ncomp[0]=2;
    ncomp[1]=1;
    ncomp[2]=1;
    
    //  cumulative number of components in blocks
    cncomp[0]=0;
    cncomp[1]=2;
    cncomp[2]=3;
    
    nip[0][0]=1;  nip[0][1]=3;  nip[0][2]=0;
    nip[1][0]=3;  nip[1][1]=3;  nip[1][2]=0;
    nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=1;
    
    intordsm[0][0]=1;  intordsm[0][1]=3;  intordsm[0][2]=0;
    intordsm[1][0]=3;  intordsm[1][1]=3;  intordsm[1][2]=0;
    intordsm[2][0]=0;  intordsm[2][1]=0;  intordsm[2][2]=1;
    
  }
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
}

axisymqt::~axisymqt (void)
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

   @param areacoord - %vector containing area coordinates
   @param nodval - nodal values
   
   JK
*/
double axisymqt::approx (vector &areacoord,vector &nodval)
{
  double f,xi,eta;
  vector bf(nne);

  xi=areacoord[0];
  eta=areacoord[1];

  bf_quad_3_2d (bf.a,xi,eta);

  //  scalar/dot product of vector of approximation
  //  functions and vector of nodal values
  scprd (bf,nodval,f);

  return f;
}

/**
   function approximates function defined by nodal values

   @param xi,eta - natural coordinates
   @param nodval - nodal values
   
   JK
*/
double axisymqt::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_3_2d (bf.a,xi,eta);
  scprd (bf,nodval,f);

  return f;
}

/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param xi,eta - area or natural coordinates
   
   JK, 29.3.2002
*/
void axisymqt::bf_matrix (matrix &n,double xi,double eta)
{
  vector bf(nne);
  
  bf_quad_3_2d (bf.a,xi,eta);

  fillm (0.0,n);
  
  n[0][0]=bf[0];
  n[0][2]=bf[1];
  n[0][4]=bf[2];
  n[0][6]=bf[3];
  n[0][8]=bf[4];
  n[0][10]=bf[5];
  
  n[1][1]=bf[0];
  n[1][3]=bf[1];
  n[1][5]=bf[2];
  n[1][7]=bf[3];
  n[1][9]=bf[4];
  n[1][11]=bf[5];
}

/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param areacoord - area coordinates
   
   JK, 3. 11. 2020
*/
void axisymqt::bf_matrix (matrix &n,vector &areacoord)
{
  double xi,eta;
  vector bf(nne);

  xi=areacoord[0];
  eta=areacoord[1];
  
  bf_quad_3_2d (bf.a,xi,eta);

  fillm (0.0,n);
  
  n[0][0]=bf[0];
  n[0][2]=bf[1];
  n[0][4]=bf[2];
  n[0][6]=bf[3];
  n[0][8]=bf[4];
  n[0][10]=bf[5];
  
  n[1][1]=bf[0];
  n[1][3]=bf[1];
  n[1][5]=bf[2];
  n[1][7]=bf[3];
  n[1][9]=bf[4];
  n[1][11]=bf[5];
}

/**
   function assembles geometric %matrix
   
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx

   @param gm - geometric %matrix
   @param xi, eta - natural coordinates
   @param x,y - node coordinates
   
   JK, 3. 11. 2020
*/
void axisymqt::geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y,double &jac)
{
  double r;
  vector n(nne),dx(nne),dy(nne);
  
  //  radius
  r = approx_nat (xi,eta,x);
  if (fabs(r)<Mp->zero){
    fprintf (stderr,"\n\n radius is equal %e in function axisymqt::geom_matrix (file %s, line %d)",r,__FILE__,__LINE__);
  }
  
  bf_quad_3_2d (n.a,xi,eta);
  dx_bf_quad_3_2d (dx.a,xi,eta);
  dy_bf_quad_3_2d (dy.a,xi,eta);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);

  fillm (0.0,gm);
  
  gm[0][0]=dx[0];
  gm[0][2]=dx[1];
  gm[0][4]=dx[2];
  gm[0][6]=dx[3];
  gm[0][8]=dx[4];
  gm[0][10]=dx[5];
  
  gm[1][1]=dy[0];
  gm[1][3]=dy[1];
  gm[1][5]=dy[2];
  gm[1][7]=dy[3];
  gm[1][9]=dy[4];
  gm[1][11]=dy[5];
  
  gm[2][0]=n[0]/r;
  gm[2][2]=n[1]/r;
  gm[2][4]=n[2]/r;
  gm[2][6]=n[3]/r;
  gm[2][8]=n[4]/r;
  gm[2][10]=n[5]/r;
  
  gm[3][0]=dy[0];
  gm[3][1]=dx[0];
  gm[3][2]=dy[1];
  gm[3][3]=dx[1];
  gm[3][4]=dy[2];
  gm[3][5]=dx[2];
  gm[3][6]=dy[3];
  gm[3][7]=dx[3];
  gm[3][8]=dy[4];
  gm[3][9]=dx[4];
  gm[3][10]=dy[5];
  gm[3][11]=dx[5];
}



/**
   nutno otestovat! pak je mozne smazat tuto hlasku
   
   transformation %matrix x_g = T x_l

   17.8.2001
*/
void axisymqt::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of triangular axisymmetric
   finite element with quadratic approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   
   JK, 3. 11. 2020
*/
void axisymqt::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,ipp;
  double jac,r;
  vector gp1,gp2,w,x(nne),y(nne);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  
  //  cleaning of stiffness matrix
  fillm (0.0,sm);
  
  //  coordinates and weights of integration points
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp1);
  reallocv (intordsm[0][0],gp2);
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);
  
  //  id of first integration point
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){

    //  geometric matrix
    geom_matrix (gm,gp1[i],gp2[i],x,y,jac);
    
    //  matrix of stiffness of the material
    Mm->matstiff (d,ipp);
    
    //  radius
    r = approx_nat (gp1[i],gp2[i],x);
    jac=w[i]*r*jac;
    
    //  contribution to the stiffness matrix of the element
    bdbjac (sm,gm,d,gm,jac);
    
    ipp++;
  }
}

/**
   function computes stiffness %matrix of triangular axisymmetric
   finite element with linear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   
   JK, 17.8.2001
*/
/*
void axisymqt::stiffness_matrix_block (long eid,long ri,long ci,matrix &sm)
{
  long i,ii,jj,ipp;
  double jac,det,r;
  vector gp1,gp2,w,x(nne),y(nne),areacoord(3);
  matrix gm(tncomp,ndofe),gmr,gmc,dd,d(tncomp,tncomp);
  
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  //  cleaning of stiffness matrix
  fillm (0.0,sm);
  
  
  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gmr);
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      //  coordinates and weights of integration points
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      
      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	areacoord[0]=gp1[i];
	areacoord[1]=gp2[i];
	areacoord[2]=1.0-areacoord[0]-areacoord[1];
	
	//  geometric matrices
	geom_matrix_block (gmr,ii,areacoord,x,y);
	geom_matrix_block (gmc,jj,areacoord,x,y);
	
	//  matrix of stiffness of the material
	Mm->matstiff (d,ipp);
	dmatblock (ii,jj,d,dd);
	
	//  radius
	r = approx (areacoord,x);
	jac=w[i]*r*det;
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,gm,d,gm,jac);
	
	ipp++;
      }
    }
  }
}
*/

/**
   function computes resulting stiffness %matrix of element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 5. 11. 2020
*/
void axisymqt::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(nne);

  stiffness_matrix (eid,0,0,sm);
  
  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  //
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}


/**
   function computes mass %matrix of the triangular axisymmetric
   finite element with linear approximation functions
   
   @param eid - number of element
   @param mm - mass %matrix

   JK, 5. 11. 2020
*/
void axisymqt::mass_matrix (long eid,matrix &mm)
{
  long i;
  double jac,det,rho,r;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp1(intordmm),gp2(intordmm),dens(nne);
  matrix n(napfun,ndofe);
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  coordinates of nodes
  Mt->give_node_coord2d (x,y,eid);
  //  nodal densities
  Mc->give_density (eid,nodes,dens);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  //  coordinates and weights of integration points
  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);
  
  //  cleaning of mass matrix
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    //  matrix of approximation functions
    bf_matrix (n,gp1[i],gp2[i]);
    //  density at integration point
    rho = approx_nat (gp1[i],gp2[i],dens);
    //  radius
    r = approx_nat (gp1[i],gp2[i],x);
    
    jac=r*w[i]*rho*det;
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
  
}


/**
   function computes resulting mass %matrix of element
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK, 5. 11. 2020
*/
void axisymqt::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  
  mass_matrix (eid,mm);
  
  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  //
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}


/**
   function integrates N^T c N over edges
   
   @param edg - edge id (number of edge)
   @param x, y - coordinates of element nodes
   @param intord - order of numerical integration
   @param gp, w  - coordinates and weights of integration points
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK, 5. 11. 2020
*/
void axisymqt::edge_integral (long edg,vector &x,vector &y,long /*intord*/,vector &gp,vector &w,
			      vector &/*coef*/,matrix &km)
{
  long i;
  double r,ww,jac;
  vector areacoord(3);
  matrix n(napfun,ndofe);
  
  //  cleaning of axiliary matrix
  fillm (0.0,km);
  
  if (edg==0){
    //  set up of area coordinates
    areacoord[2]=0.0;
    for (i=0;i<intordb;i++){
      areacoord[1]=(1.0+gp[i])/2.0;  areacoord[0]=1.0-areacoord[1];
      ww=w[i];
      
      //  matrix of approximation functions
      bf_matrix (n,areacoord);
      //  computation of jacobian
      jac1d_2d (jac,x,y,areacoord[1],0);
      //  radius
      r=approx (areacoord,x);

      jac*=ww*r;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  if (edg==1){
    //  set up of area coordinates
    areacoord[0]=0.0;
    for (i=0;i<intordb;i++){
      areacoord[1]=(1.0+gp[i])/2.0;  areacoord[2]=1.0-areacoord[1];
      ww=w[i];
      
      //  matrix of approximation functions
      bf_matrix (n,areacoord);
      //  computation of jacobian
      jac1d_2d (jac,x,y,areacoord[0],1);
      //  radius
      r=approx (areacoord,x);

      jac*=ww*r;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  if (edg==2){
    //  set up of area coordinates
    areacoord[1]=0.0;
    for (i=0;i<intordb;i++){
      areacoord[0]=(1.0+gp[i])/2.0;  areacoord[2]=1.0-areacoord[0];
      ww=w[i];
      
      //  matrix of approximation functions
      bf_matrix (n,areacoord);
      //  computation of jacobian
      jac1d_2d (jac,x,y,areacoord[0],2);
      //  radius
      r=approx (areacoord,x);

      jac*=ww*r;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  
}


/**
   function picks up nodal values on required edges
   
   @param edg - number of required edge
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all edges
   
   JK, 5. 11. 2020
*/
void axisymqt::edgenodeval (long edg,vector &nodval,double *list)
{
  long i,j,k,l;
  ivector edgenod(nned);
  
  fillv (0.0,nodval);
  //  nodes on required edge
  quadtriangle_edgnod (edgenod.a,edg);
  
  k=0;
  for (i=0;i<nned;i++){
    l=edgenod[i]*napfun;
    for (j=0;j<napfun;j++){
      nodval[l+j]=list[edg*nned*napfun+k];
      k++;
    }
  }
}

/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - indicators of loaded edges
   @param list - list of prescribed nodal values
   @param nf - %vector of nodal forces caused by edge load
   
   JK, 5. 11. 2020
*/
void axisymqt::edgeload (long eid,long *le,double *list,vector &nf)
{
  long i;
  vector x(nne),y(nne),gp(intordb),w(intordb),nodval(ndofe),coeff(nne),av(ndofe);
  matrix km(ndofe,ndofe);

  //  coefficient c in N^T c N
  for (i=0;i<nne;i++){
    coeff[i]=1.0;
  }
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordb);
  
  //  cleaning of nf
  fillv (0.0,nf);

  for (i=0;i<ned;i++){
    if (le[i]>0){
      fillm(0.0, km);
      //  selection of appropriate nodal values
      edgenodeval (i,nodval,list);
      //  integration over element edge
      edge_integral (i,x,y,intordb,gp,w,coeff,km);
      
      mxv (km,nodval,av);
      
      addv (nf,av,nf);
    }
  }
  
}
































































void axisymqt::res_mainip_strains (long lcid,long eid)
{
  vector aux,x(nne),y(nne),r(ndofe);
  ivector nodes(nne);
  matrix tmat;

  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
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

  mainip_strains (lcid,eid,0,0,x,y,r);  
}


/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void axisymqt::mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  double jac;
  vector gp1,gp2,w,eps(tncomp),aux;
  matrix gm(tncomp,ndofe);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    //reallocv (ncomp[ii],eps);
    //reallocm (ncomp[ii],ndofe,gm);

    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){

      geom_matrix (gm,gp1[i],gp2[i],x,y,jac);
      
      mxv (gm,r,eps);
      
      //Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
  }
}

/**
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void axisymqt::nod_strains_ip (long lcid,long eid)
{
  long i,j;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
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
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 23.9.2004
*/
void axisymqt::nod_strains_comp (long lcid,long eid)
{
  long i,j;
  double jac;
  vector x(nne),y(nne),nxi(nne),neta(nne),r(ndofe),eps(tncomp),aux;
  ivector nodes(nne);
  matrix gm(tncomp,ndofe),tmat;
  
  //  natural coordinates of nodes of element
  nodecoord (nxi,neta);
  //  node numbers of element
  Mt->give_elemnodes (eid,nodes);
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  
  for (i=0;i<nne;i++){
    //  block of geometric matrix
    geom_matrix (gm,nxi[i],neta[i],x,y,jac);
    //  strain computation
    mxv (gm,r,eps);
    
    //  storage of nodal strains
    j=nodes[i];
    Mt->nodes[j].storestrain (lcid,0,eps);
  }
}



/**
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqt::res_allip_strains (long lcid,long eid)
{
  allip_strains (lcid,eid,0,0);
}
/**
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqt::allip_strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
}


void axisymqt::strains (long /*lcid*/,long eid,long /*ri*/,long /*ci*/)
{
  vector coord,eps;
  
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
    /*
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (ncp,eps);
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	appval (coord[0],coord[1],0,ncp,eps,stra);
      if (Mp->strainaver==1)
      appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlt::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}

/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   
   10.5.2002
*/
void axisymqt::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 0.0;  eta[0] = 0.0;
  xi[1] = 1.0;  eta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 1.0;
}


/**
   function returns numbers of integration point closest to element nodes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipnum - array of numbers
   
   JK, 25.9.2004
*/
void axisymqt::nodipnum (long eid,ivector &ipnum)
{
  long i,j;
  
  j=intordsm[0][0];
  i=Mt->elements[eid].ipp[0][0];
  
  /*
  ipnum[0]=i+j*(j-1)+j-1;
  ipnum[1]=i+j-1;
  ipnum[2]=i;
  ipnum[3]=i+j*(j-1);
  */
  
  ipnum[0]=i+0;
  ipnum[1]=i+0;
  ipnum[2]=i+0;
  
}


/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void axisymqt::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
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
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void axisymqt::res_mainip_stresses (long lcid,long eid)
{
  mainip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void axisymqt::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ipp;
  vector gp,w,eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    Mm->matstiff (d,ipp);
    
    Mm->givestrain (lcid,ipp,eps);
    
    mxv (d,eps,sig);
    
    Mm->storestress (lcid,ipp,sig);
    
    ipp++;
  }
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqt::nod_stresses_ip (long lcid,long eid)
{
  long i,j;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
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


void axisymqt::nod_stresses_comp (long /*lcid*/,long /*eid*/)
{
/*
  long i,j,ipp;
  double jac;
  vector x(nne),y(nne),nxi(nne),neta(nne),r(ndofe),eps(tncomp),sig(tncomp),aux;
  ivector nodes(nne);
  matrix gm(tncomp,ndofe),tmat,d(tncomp,tncomp);
  
  //  natural coordinates of nodes of element
  nodecoord (nxi,neta);
  //  node numbers of element
  Mt->give_elemnodes (eid,nodes);
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    locglobtransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  
  ipp=Mt->elements[eid].ipp[0][0]+4;
  for (i=0;i<nne;i++){
    
    //Mm->givestrain (lcid,ipp,eps);
    Mm->givestress (lcid,ipp,sig);
    
    //  number of actual node
    j=nodes[i];
    //  storage of nodal strains
    //Mt->nodes[j].storestrain (lcid,0,tncomp,eps);
    //  storage of nodal strains
    Mt->nodes[j].storestress (lcid,0,tncomp,sig);
  }
*/  
}


/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqt::res_allip_stresses (long lcid,long eid)
{
  allip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqt::allip_stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  res_mainip_stresses (lcid,eid);
}


/**
   function computes stresses in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param areacoord - area coordinates of the point
   @param fi,li - first and last indices
   @param sig - array containing stresses
   
   11.5.2002
*/
void axisymqt::appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig)
{
  long i,j,k;
  ivector nodes;
  vector areacoord(3),nodval;
  
  if (ncomp != sig.n){
    fprintf (stderr,"\n\n wrong interval of indices in function stress (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  
  reallocv (nne,nodes);
  reallocv (nne,nodval);
  
  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].stress[lcid*tncomp+i];
    }
    sig[k]=approx (areacoord,nodval);
    k++;
  }
}



void axisymqt::stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  vector coord,sig;
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (stre,lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    nod_stresses_ip (lcid,eid);
    break;
  }
  case userdefined:{
    /*
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (ncp,sig);
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	appval (coord[0],coord[1],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}



/**
   function computes eqother components at nodes of element
   
   @param eid - element id
   
   10.5.2002
*/
void axisymqt::nod_other_ip (long eid)
{
  long i,j,ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other;
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Mm->givencompother (ipnum[i],0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    //  storage of other values to the node
    j=nod[i];
    Mt->nodes[j].storeother(0, ncompo, other);
  }
}



/**
   function computes load %matrix of the triangular axisymmetric
   finite element with linear approximation functions
   load %vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   
   25.7.2001
*/
void axisymqt::load_matrix (long eid,matrix &lm)
{
  long i;
  double jac,det;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp1(intordmm),gp2(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    bf_matrix (n,gp1[i],gp2[i]);
    
    //  zkontrolovat deleni dvema
    jac=w[i]*det;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  
}



/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void axisymqt::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void axisymqt::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes increments of internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void axisymqt::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   7.2008, TKo
*/
void axisymqt::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void axisymqt::res_internal_forces (long lcid,long eid,vector &ifor)
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
void axisymqt::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void axisymqt::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void axisymqt::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
void axisymqt::compute_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void axisymqt::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}



/**
   function computes nonlocal correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void axisymqt::compute_nonloc_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->compnonloc_nlstresses (ipp);
      ipp++;
    }
  }
}



/**
   function computes correct stress increments at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void axisymqt::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct increments of stresses
      if (Mp->strcomp==1)
        Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
}



void axisymqt::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
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
   TKo 7.2008
*/
void axisymqt::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double rad,jac;
  vector w,gp1,gp2;
  vector r(ndofe),ipv(tncomp),contr(ndofe),auxcontr(ndofe);
  matrix gm(tncomp,ndofe);
  
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
      //  strain-displacement (geometric) matrix
      geom_matrix (gm,gp1[i],gp2[i],x,y,jac);
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      rad = approx_nat (gp1[i],gp2[i],x);
      cmulv (rad*w[i]*jac,contr);
      //  summation
      addv(contr,nv,nv);
      ipp++;
    }
  }
}



/**
   function returns coordinates of integration points

   @param eid - element id
   @param ipp - integration point pointer
   @param ri - row index
   @param ci - column index
   @param coord - %vector of coordinates

   19.1.2002
*/
void axisymqt::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii,jj,aipp;
  vector x(nne),y(nne);
  vector areacoord(3),w,gp1,gp2;
  
  Mt->give_node_coord2d (x,y,eid);
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      //  integration point id
      aipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],w);
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      for (i=0;i<intordsm[ii][jj];i++){
        areacoord[0]=gp1[i];
        areacoord[1]=gp2[i];
        areacoord[2]=1.0-areacoord[0]-areacoord[1];
        if (aipp==ipp){
          coord[0]=approx (areacoord,x);
          coord[1]=approx (areacoord,y);
          coord[2]=0.0;
        }
        aipp++;
      }
    }
  }
}



void axisymqt::ipncoord(long eid, long ipp, vector &coord)
{
  long i, ii, jj, aipp;
  vector w, gp1, gp2;
  
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      //  integration point id
      aipp=Mt->elements[eid].ipp[ii][jj];
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv (RSTCKVEC(intordsm[ii][jj],gp1));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp2));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      for (i=0;i<intordsm[ii][jj];i++){
        if (aipp==ipp){
          coord[0]=gp1[i];
          coord[1]=gp2[i];
          coord[2]=0.0;
        }
        aipp++;
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
void axisymqt::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
    for(i = 0; i < nne; i++)
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
          //  value in the integration point
          ipval = approx_nat(xi, eta, anv);
          ncompstr =  Mm->ip[ipp].ncompstr;
          ncompeqother = Mm->ip[ipp].ncompeqother;
          if ((ict & inistrain) && (j < ncompstr))
          {
            Mm->ip[ipp].strain[idstra] += ipval;
            ipp++;
            continue;
          }
          if ((ict & inistress) && (j < nstra + ncompstr))
          {
            Mm->ip[ipp].stress[idstre] += ipval;
            ipp++;
            continue;
          }
          if ((ict & iniother) && (j < nstra+nstre+ncompeqother))
          {
            Mm->ip[ipp].eqother[idoth] += ipval;
            ipp++;
            continue;
          }
          if ((ict & inicond) && (j < nv))
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
    if ((ict & inistrain) && (j < ncompstr))
    {
      nstra++;
      idstra++;
      continue;
    }
    if ((ict & inistress) && (j < nstra + ncompstr))
    {      
      nstre++;
      idstre++;
      continue;
    }  
    if ((ict & iniother)  && (j < nstra + nstre + ncompeqother))
    {
      idoth++;
      continue;
    }
    if ((ict & inicond) && (j < nv))
    {
      idic++;
      continue;
    }
  }
}
