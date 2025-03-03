#include "dkq.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>



dkq::dkq (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=5;
  //  number of strain/stress components connected with bending
  bncomp=3;
  //  number of strain/stress components connected with shear
  sncomp=2;
  //  number of components for graphic purposes
  //gncomp=4;
  //  number of functions approximated
  napfun=1;
  //  order of numerical integration of mass matrix
  intordmm=3;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=1;
  //  number of nodes on surface
  nnsurf=4;
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
  
  intordsm[0][0]=2;
  
  nip[0][0]=4;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  //  strain/stress state
  ssst = platek;

}

dkq::~dkq (void)
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
   function prepares auxiliary data for quadrilateral plate element
   
   @param x,y - vectors of node coordinates
   @param l - %vector of lengths
   @param sx,sy - vectors of direction vectors components
   @param nx,ny - vectors of normal vectors components
   @param sign - %vector of signs
   @param nodes - element nodes
   
   JK, 20. 10. 2019
*/
void dkq::auxdata (vector &x,vector &y,vector &l,vector &sx,vector &sy,vector &nx,vector &ny,vector &signs,ivector &nodes)
{
  vector lsx(ASTCKVEC(nne)),lsy(ASTCKVEC(nne));
  
  if (nodes[1]>nodes[0]){
    sx[0]=x[1]-x[0];  sy[0]=y[1]-y[0];
  }else{
    sx[0]=x[0]-x[1];  sy[0]=y[0]-y[1];
  }
  lsx[0]=x[1]-x[0];  lsy[0]=y[1]-y[0];

  l[0]=sqrt(sx[0]*sx[0]+sy[0]*sy[0]);
  lsx[0]/=l[0];     lsy[0]/=l[0];
  sx[0]/=l[0];      sy[0]/=l[0];
  nx[0]=sy[0];      ny[0]=0.0-sx[0];
  signs[0]=sx[0]*lsx[0]+sy[0]*lsy[0];

  
  if (nodes[2]>nodes[1]){
    sx[1]=x[2]-x[1];  sy[1]=y[2]-y[1];
  }else{
    sx[1]=x[1]-x[2];  sy[1]=y[1]-y[2];
  }
  lsx[1]=x[2]-x[1];  lsy[1]=y[2]-y[1];

  l[1]=sqrt(sx[1]*sx[1]+sy[1]*sy[1]);
  lsx[1]/=l[1];     lsy[1]/=l[1];
  sx[1]/=l[1];      sy[1]/=l[1];
  nx[1]=sy[1];      ny[1]=0.0-sx[1];
  signs[1]=sx[1]*lsx[1]+sy[1]*lsy[1];


  if (nodes[3]>nodes[2]){
    sx[2]=x[3]-x[2];  sy[2]=y[3]-y[2];
  }else{
    sx[2]=x[2]-x[3];  sy[2]=y[2]-y[3];
  }
  lsx[2]=x[3]-x[2];  lsy[2]=y[3]-y[2];
  
  l[2]=sqrt(sx[2]*sx[2]+sy[2]*sy[2]);
  lsx[2]/=l[2];     lsy[2]/=l[2];
  sx[2]/=l[2];      sy[2]/=l[2];
  nx[2]=sy[2];      ny[2]=0.0-sx[2];
  signs[2]=sx[2]*lsx[2]+sy[2]*lsy[2];


  if (nodes[0]>nodes[3]){
    sx[3]=x[0]-x[3];  sy[3]=y[0]-y[3];
  }else{
    sx[3]=x[3]-x[0];  sy[3]=y[3]-y[0];
  }
  lsx[3]=x[0]-x[3];  lsy[3]=y[0]-y[3];
  
  l[3]=sqrt(sx[3]*sx[3]+sy[3]*sy[3]);
  lsx[3]/=l[3];     lsy[3]/=l[3];
  sx[3]/=l[3];      sy[3]/=l[3];
  nx[3]=sy[3];      ny[3]=0.0-sx[3];
  signs[3]=sx[3]*lsx[3]+sy[3]*lsy[3];

}

/**
   function approximates function defined by nodal values

   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
   JK, 20. 10. 2019
*/
double dkq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function assembles %matrix of base (approximation) functions
   there are bi-cubic shape functions

   the functions are used inside the element because the functions
   used for the stiffness are defined only on edges
   
   @param n - base function %matrix
   @param xi, eta - natural coordinates
   
   JK, 30. 11. 2019
*/
void dkq::bf_deflection_matrix (matrix &n,double xi,double eta)
{
  bf_cubic_plate_4_2d (n.a,xi,eta);
}

/**
   function assembles %matrix of base (approximation) functions
   there are bi-linear shape functions

   the functions are used for description of surface load
   
   @param n - base function %matrix
   @param xi, eta - natural coordinates
   
   JK, 30. 11. 2019
*/
void dkq::bf_linear_matrix (matrix &n,double xi,double eta)
{
  bf_lin_4_2d (n.a,xi,eta);
}


/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param xi,eta - natural coordinates
   @param x,y - coordinates of nodes on the element
   @param sx,sy - coordinates of direction vectors
   @param nx,ny - coordinates of normal vectors
   @param l - %vector of edge lengths
   @param signs - %vector of signs
   @param jac - Jacobian
   
   JK, 20. 10. 2019
*/
void dkq::geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y,vector &sx,vector &sy,vector &nx,vector &ny,vector &l,vector &signs,double &jac)
{
  long i;
  vector dx(ASTCKVEC(ndofe*2)),dy(ASTCKVEC(ndofe*2));
  
  dx_bf_quad_dkq (dx.a,xi,eta,x,y,l.a,sx.a,sy.a,nx.a,ny.a,signs.a);
  dy_bf_quad_dkq (dy.a,xi,eta,x,y,l.a,sx.a,sy.a,nx.a,ny.a,signs.a);
  
  //  Jacobian
  jac_2d (jac,x,y,xi,eta);

  for (i=0;i<ndofe;i++){
    //  d2w/dx2
    gm[0][i] = dx[i];
    //  d2w/dy2
    gm[1][i] = dy[ndofe+i];
    //  d2w/dx/dy
    gm[2][i] = dx[ndofe+i] + dy[i];
  }
  
}


/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param xi,eta - natural coordinates
   @param x,y - coordinates of nodes on the element
   @param sx,sy - coordinates of direction vectors
   @param nx,ny - coordinates of normal vectors
   @param l - %vector of edge lengths
   @param signs - %vector of signs
   
   JK, 15. 12. 2019
*/
void dkq::shear_geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y,vector &sx,vector &sy,vector &nx,vector &ny,vector &l,vector &signs)
{
  long i;
  vector dxdx(ASTCKVEC(ndofe*2)),dydy(ASTCKVEC(ndofe*2));
  
  dxdx_bf_quad_dkq (dxdx.a,xi,eta,x,y,l.a,sx.a,sy.a,nx.a,ny.a,signs.a);
  dydy_bf_quad_dkq (dydy.a,xi,eta,x,y,l.a,sx.a,sy.a,nx.a,ny.a,signs.a);
  
  for (i=0;i<ndofe;i++){
    //  d3w/dx3
    gm[0][i] = dxdx[i];
    //  d3w/dx2dy
    gm[1][i] = dxdx[ndofe+i];
    //  d3w/dx/dy2
    gm[2][i] = dydy[i];
    //  d3w/dy3
    gm[3][i] = dydy[ndofe+i];
  }
  
}

/**
   nutna kontrola

   function assembles transformation %matrix 
*/
void dkq::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,n,m;
  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*3+1][i*3+1]=Mt->nodes[nodes[i]].e1[0];  tmat[i*3+1][i*3+2]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*3+2][i*3+1]=Mt->nodes[nodes[i]].e1[1];  tmat[i*3+2][i*3+2]=Mt->nodes[nodes[i]].e2[1];
    }
  }
}


/**
   function computes stiffness %matrix of plate quadrilateral
   finite element based on discrete Kirchhoff theory
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of nodal coordinates

   JK, 20. 10. 2019
*/
void dkq::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,ww1,ww2,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector l(ASTCKVEC(nne)),sx(ASTCKVEC(nne)),sy(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),w,gp,t(ASTCKVEC(nne)),signs(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(bncomp,ndofe)),d(ASTCKMAT(bncomp,bncomp));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  thickness in nodes
  Mc->give_thickness (eid,nodes,t);
  //  auxiliary vectors (direction vectors, normal vectors and edge lengths
  auxdata (x,y,l,sx,sy,nx,ny,signs,nodes);
  
  fillm (0.0,sm);
  
  //  coordinates of integration points and their weights
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];  ww2=w[j];
      
      //  geometric matrix
      geom_matrix (gm,xi,eta,x,y,sx,sy,nx,ny,l,signs,jac);
      
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      //  thickness in integration point
      thick = approx (xi,eta,t);
      
      jac*=thick*thick*thick*ww1*ww2;
      
      //  contribution to the stiffness matrix of the element
      bdbjac (sm,gm,d,gm,jac);
      
      ipp++;
    }
  }
  
  /*
  fprintf (Out,"\n\n matice tuhosti \n");
  for (i=0;i<ndofe;i++){
    for (j=0;j<ndofe;j++){
      fprintf (Out,"%2ld %2ld  % 15.12le\n",i,j,sm[i][j]);
    }
    fprintf (Out,"\n");
  }
  fprintf (Out,"\n");
  */

}

/**
   function computes stiffness %matrix of quadrilateral plate
   element based on the discrete Kirchhoff theory
   
   @param lcid - load case id
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 20. 10. 2019
*/
void dkq::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  
  //  stiffness matrix
  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness %matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
}

/**
   function computes nodal forces and moments caused by
   surface load on element
   
   @param eid - element id
   @param nodvals - array of load nodal values
   @param x,y - vectors of nodal coordinates
   @param nf - nodal forces and moments
   
   JK, 30. 11. 2019
*/
void dkq::surfload (long /*eid*/, double *nodvals, vector &x, vector &y, vector &nf)
{
  long i,j;
  double xi,eta,w1,w2,jac;
  vector gp(ASTCKVEC(intordmm)),w(ASTCKVEC(intordmm));
  matrix nlin(ASTCKMAT(napfun,nne)),ncub(ASTCKMAT(napfun,ndofe)),am(ASTCKMAT(ndofe,nne));
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordmm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      
      jac_2d (jac,x,y,xi,eta);
      
      bf_deflection_matrix (ncub,xi,eta);
      bf_linear_matrix (nlin,xi,eta);
      
      jac*=w1*w2;
      
      nnjac (am,ncub,nlin,jac);
    }
  }
  
  mxv (am.a,nodvals,nf.a,ndofe,nne);
}

/**
   function coputes nodal forces caused by surface load
   
   @param eid - element id
   @param nodvals - array of nodal values of the surface load
   @param nf - %vector of nodal forces
   
   JK, 1. 12. 2019
*/
void dkq::node_forces_surf (long eid,double *nodvals,vector &nf)
{
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  
  //  nodal force caused by surface load
  surfload (eid,nodvals,x,y,nf);  
}


/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - id of loaded edge
   @param nv - nodal values
   @param nf - %vector of nodal forces
   
   JK, 9. 3. 2020
*/
void dkq::nodeforces (long /*eid*/,long *le,double *nv,vector &x,vector &y,vector &nf)
{
  long i;
  double ww,jac,xi,eta;
  vector gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),av(ASTCKVEC(nne)),v(ASTCKVEC(ndofe));
  matrix n1(ASTCKMAT(1,nne)),n3(ASTCKMAT(1,ndofe)),am(ASTCKMAT(ndofe,nne));
  
  //  integration point coordinates and weights
  gauss_points (gp.a,w.a,intordb);
  
  
  if (le[0]==1){
    fillm (0.0,am);
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_deflection_matrix (n3,xi,eta);
      bf_linear_matrix (n1,xi,eta);
      
      jac1d_2d (jac,x,y,xi,0);
      jac *= ww;
      
      nnjac (am,n3,n1,jac);
    }
    fillv (0.0,av);
    av[0]=nv[0];  av[1]=nv[1];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    fillm (0.0,am);
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];

      bf_deflection_matrix (n3,xi,eta);
      bf_linear_matrix (n1,xi,eta);
      
      jac1d_2d (jac,x,y,eta,1);
      jac *= ww;
      
      nnjac (am,n3,n1,jac);
    }
    fillv (0.0,av);
    av[1]=nv[2];  av[2]=nv[3];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    fillm (0.0,am);
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_deflection_matrix (n3,xi,eta);
      bf_linear_matrix (n1,xi,eta);
      
      jac1d_2d (jac,x,y,xi,2);
      jac *= ww;
      
      nnjac (am,n3,n1,jac);
    }
    fillv (0.0,av);
    av[2]=nv[4];  av[3]=nv[5];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    fillm (0.0,am);
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];
      
      bf_deflection_matrix (n3,xi,eta);
      bf_linear_matrix (n1,xi,eta);
      
      jac1d_2d (jac,x,y,eta,3);
      jac *= ww;

      nnjac (am,n3,n1,jac);
    }
    fillv (0.0,av);
    av[3]=nv[6];  av[0]=nv[7];
    mxv (am,av,v);  addv (nf,v,nf);
  }
}

/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - id of loaded edge
   @param nv - nodal values
   @param nf - %vector of nodal forces
   
   JK, 9. 3. 2020
*/
void dkq::node_forces_edge (long eid,long *le,double *nv,vector &nf)
{
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  
  nodeforces (eid,le,nv,x,y,nf);
}



/**
   function computes curvatures
      
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column indices
   @param x,y - vectors of nodal coordinates
   @param r - nodal values
   

   JK, 30. 11. 2019
*/
void dkq::ip_curvatures (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ipp;
  double xi,eta,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector l(ASTCKVEC(nne)),sx(ASTCKVEC(nne)),sy(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),w,gp,t(ASTCKVEC(nne)),signs(ASTCKVEC(nne)),eps(ASTCKVEC(bncomp));
  matrix gm(ASTCKMAT(bncomp,ndofe));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  thickness in nodes
  Mc->give_thickness (eid,nodes,t);
  //  auxiliary vectors (direction vectors, normal vectors and edge lengths
  auxdata (x,y,l,sx,sy,nx,ny,signs,nodes);
  
  //  coordinates of integration points and their weights
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  geometric matrix
      geom_matrix (gm,xi,eta,x,y,sx,sy,nx,ny,l,signs,jac);
      
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,cncomp[0],eps);
      
      ipp++;
    }
  }
}

/**
   function computes curvatures
      
   @param lcid - load case id
   @param eid - element id
   
   JK, 1. 3. 2020
*/
void dkq::res_ip_curvatures (long lcid,long eid)
{
  vector aux,x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  
  //  nodal values
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
  
  
  ip_curvatures (lcid,eid,0,0,x,y,r);
}


/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 30. 1. 2019
*/
void dkq::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double thick,xi,eta;
  ivector nodes(nne);
  vector sig(3),gp,w,t(ASTCKVEC(nne));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp);

    gauss_points (gp.a,w.a,intordsm[0][0]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	//  computation of correct stresses
	if (Mp->strcomp==1){
	  
	  //  thickness in integration point
	  thick = approx (xi,eta,t);
	  
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	  
	  Mm->givestress (lcid,ipp,sig);
	  
	  sig[0]*=thick*thick*thick;
	  sig[1]*=thick*thick*thick;
	  sig[2]*=thick*thick*thick;
	  
	  Mm->storestress (lcid,ipp,sig);
	  
	  ipp++;
	}
      }
    }
  }
}


/**
   function computes internal moments
   
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column indices
   
   JK, 30. 11. 2019
*/
void dkq::moments (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  double xi,eta,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  vector eps(ASTCKVEC(bncomp)),sig(ASTCKVEC(bncomp));
  matrix d(ASTCKMAT(bncomp,bncomp));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  thickness in nodes
  Mc->give_thickness (eid,nodes,t);
  
  //  coordinates of integration points and their weights
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      Mm->givestrain (0,ipp,eps);
      
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      mxv (d,eps,sig);
      
      //  thickness in integration point
      thick = approx (xi,eta,t);
      
      cmulv (thick*thick*thick,sig);
      
      Mm->storestress (0,ipp,sig);
      
      ipp++;
    }
  }
  

}

/**
   function computes internal moments
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 1. 3. 2020
*/
void dkq::res_moments (long lcid,long eid)
{
  //  internal moments
  moments (lcid,eid,0,0);
  
}


/**
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column indices
   @param x,y - vectors of nodal coordinates
   @param r - nodal values
   
   JK, 15. 12. 2019
*/
void dkq::forces (long /*lcid*/,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ipp;
  double xi,eta,thick,e,nu;
  ivector nodes(ASTCKIVEC(nne));
  vector l(ASTCKVEC(nne)),sx(ASTCKVEC(nne)),sy(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),w,gp,t(ASTCKVEC(nne)),signs(ASTCKVEC(nne));
  vector eps(ASTCKVEC(4)),sig(ASTCKVEC(2));
  matrix gm(ASTCKMAT(4,ndofe)),d(ASTCKMAT(2,4));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  thickness in nodes
  Mc->give_thickness (eid,nodes,t);
  //  auxiliary vectors (direction vectors, normal vectors and edge lengths
  auxdata (x,y,l,sx,sy,nx,ny,signs,nodes);
  
  //  coordinates of integration points and their weights
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  geometric matrix
      shear_geom_matrix (gm,xi,eta,x,y,sx,sy,nx,ny,l,signs);
      
      mxv (gm,r,eps);

      //  Young's modulus
      e=Mm->give_actual_ym (ipp,0,0);
      //  Poisson's number
      nu=Mm->give_actual_nu (ipp,0,0);

      d[0][0]=1.0;
      d[0][1]=0.0;
      d[0][2]=1.0+nu;
      d[0][3]=0.0;
      
      d[1][0]=0.0;
      d[1][1]=1.0+nu;
      d[1][2]=0.0;
      d[1][3]=1.0;
      
      //  thickness in integration point
      thick = approx (xi,eta,t);

      cmulm (e*thick*thick*thick/12.0/(1.0-nu*nu),d);
      
      mxv (d,eps,sig);
      
      Mm->storestress (0,ipp,bncomp,sig);
      
      ipp++;
    }
  }
  
}

/**
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 1. 3. 2020
*/
void dkq::res_forces (long lcid,long eid)
{
  vector aux,x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  
  //  nodal values
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
  
  //  internal shear forces
  forces (lcid,eid,0,0,x,y,r);
  
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
   
   TKo 26. 2. 2020
*/
void dkq::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector l(ASTCKVEC(nne)),sx(ASTCKVEC(nne)),sy(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),w,gp,signs(ASTCKVEC(nne));
  vector ipv,contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(bncomp,ndofe));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  auxiliary vectors (direction vectors, normal vectors and edge lengths
  auxdata (x,y,l,sx,sy,nx,ny,signs,nodes);
  
  //  coordinates of integration points and their weights
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];  ww2=w[j];
      
      //  geometric matrix
      geom_matrix (gm,xi,eta,x,y,sx,sy,nx,ny,l,signs,jac);
      
      jac*=ww1*ww2;
      
      reallocv (ncomp[0],ipv);
      
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);
      
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      cmulv (jac,contr);
      
      //  summation
      addv(contr,nv,nv);
      ipp++;
    }
  }
}



/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 26. 2. 2020
*/
void dkq::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}


/**
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 1. 3. 2020
*/
void dkq::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  
  //  internal forces
  internal_forces (lcid,eid,0,0,ifor,x,y);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}
