#include "plelemqt.h"
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
#include "intpoints.h"
#include "plelemlt.h"
#include "plelemsubqt.h"
#include "plelemlq.h"
#include "plelemqq.h"
#include "loadcase.h"
#include "gadaptivity.h"
#include <stdlib.h>
#include <math.h>


planeelemqt::planeelemqt (void)
{
  long i,j;
  
  //  number nodes on element
  nne=6;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=2;
  //  order of numerical integration of mass matrix
  intordmm=6;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;

  //  number of blocks (parts of geometric matrix)
  nb=2;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=2;
  ncomp[1]=1;

  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=2;

  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=3;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=3;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=3;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=3;

}

planeelemqt::~planeelemqt (void)
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

   @param xi,eta - natural coordinates
   @param nodval - nodal values
   
   1.4.2002
*/
double planeelemqt::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_3_2d (bf.a,xi,eta);
  scprd (bf,nodval,f);
  return f;
}

/**
   function assembles %matrix of approximation function
   
   @param n - %matrix of approximation functions
   @param xi,eta - natural coordinates
   
   17.8.2001
*/
void planeelemqt::bf_matrix (matrix &n,double xi,double eta)
{
  long i,i1,i2;
  vector bf(nne);
  
  bf_quad_3_2d (bf.a,xi,eta);
  fillm (0.0,n);
  
  i1=0;  i2=1;
  for (i=0;i<nne;i++){
    n[0][i1]=bf[i];  i1+=2;
    n[1][i2]=bf[i];  i2+=2;
  }
}

/**
   function assembles geometric %matrix
   
   @param gm - geometric %matrix
   @param x,y - node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   
   1.4.2002
*/
void planeelemqt::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  vector dx(nne),dy(nne);
  
  dx_bf_quad_3_2d (dx.a,xi,eta);
  dy_bf_quad_3_2d (dy.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);

  fillm (0.0,gm);

  i1=0;  i2=1;
  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[1][i2]=dy[i];
    gm[2][i1]=dy[i];  i1+=2;
    gm[2][i2]=dx[i];  i2+=2;
  }
}

/**
   function assembles geometric %matrix
   
   @param gm - geometric %matrix
   @param x,y - node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   
   1.4.2002
*/
void planeelemqt::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  vector dx(nne),dy(nne);
  
  dx_bf_quad_3_2d (dx.a,xi,eta);
  dy_bf_quad_3_2d (dy.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);

  fillm (0.0,gm);

  if (ri==0){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dx[i];  i1+=2;
      gm[1][i2]=dy[i];  i2+=2;
    }
  }
  
  if (ri==1){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dy[i];  i1+=2;
      gm[0][i2]=dx[i];  i2+=2;
    }
  }
}

/**
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
*/
void planeelemqt::dmatblock (long ri,long ci,matrix &d, matrix &dd)
{
  fillm (0.0,dd);
  
  if (ri==0 && ci==0){
    dd[0][0]=d[0][0];  dd[0][1]=d[0][1];
    dd[1][0]=d[1][0];  dd[1][1]=d[1][1];
  }
  if (ri==0 && ci==1){
    dd[0][0]=d[0][2];
    dd[1][0]=d[1][2];
  }
  if (ri==1 && ci==0){
    dd[0][0]=d[2][0];  dd[0][1]=d[2][1];
  }
  if (ri==1 && ci==1){
    dd[0][0]=d[2][2];
  }
}

/**
   function assembles transformation %matrix x_g = T x_l
   
   17.8.2001
*/
void planeelemqt::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of triangular
   finite element with quadratic approximation functions

   @param eid - element id
   @param sm - stiffness %matrix

   25.8.2001
*/
void planeelemqt::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii,jj,ipp;
  double jac,thick;
  ivector nodes(nne);
  vector t(nne),gp1,gp2,w;
  matrix gmr,gmc,dd,d(tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  fillm (0.0,sm);

  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gmr);
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],w);

      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	// geometric matrix
	geom_matrix_block (gmr,ii,x,y,gp1[i],gp2[i],jac);
	geom_matrix_block (gmc,jj,x,y,gp1[i],gp2[i],jac);
	
	//  stiffness matrix of material
	Mm->matstiff (d,ipp);
	dmatblock (ii,jj,d,dd);
	
	//  thickness
	thick = approx (gp1[i],gp2[i],t);
	
	jac*=w[i]*thick;
	
	//fprintf (stdout,"\n jakobian  %lf",jac);
	
	//  contribution to the stiffness matrix of the element
	//bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
	bdbjac (sm,gmr,dd,gmc,jac);
	
	ipp++;
      }
    }
  }
}

void planeelemqt::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  vector x(nne),y(nne);
  ivector nodes(nne);
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);

  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes mass %matrix of triangular
   finite element with quadratic approximation functions

   @param eid - element id
   @param mm - mass %matrix
   @param x,y - vectors of nodal coordinates

   25.8.2001
*/
void planeelemqt::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double jac,thick,rho;
  ivector nodes(nne);
  vector w(intordmm),gp1(intordmm),gp2(intordmm),t(nne),dens(nne);
  matrix n(napfun,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,mm);

  for (i=0;i<intordmm;i++){
    //  matrix of approximation functions
    bf_matrix (n,gp1[i],gp2[i]);

    //  thickness
    thick = approx (gp1[i],gp2[i],t);

    //  density
    rho = approx (gp1[i],gp2[i],dens);

    //  Jacobian
    jac_2d (jac,x,y,gp1[i],gp2[i]);

    jac*=w[i]*thick*rho;

    //  N^T.N multiplication
    nnj (mm.a,n.a,jac,n.m,n.n);
  }

}


/**
   function computes mass %matrix of the plane stress triangular
   finite element with linear approximation functions
   
   @param eid - number of element
   @param mm - mass %matrix
   
   JK, 14. 6. 2019
*/
void planeelemqt::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
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
   function computes load %matrix of triangular
   finite element with quadratic approximation functions
   load %vector is obtained after premultiplying load %matrix
   by nodal load values

   @param eid - number of element
   @param lm - load %matrix

   25.8.2001
*/
void planeelemqt::load_matrix (long eid,matrix &lm)
{
  long i;
  double jac,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp1(intordmm),gp2(intordmm),t(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    bf_matrix (n,gp1[i],gp2[i]);
    
    thick = approx (gp1[i],gp2[i],t);
    
    jac_2d (jac,x,y,gp1[i],gp2[i]);
    
    jac*=w[i]*thick;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  
}



/**
   function computes stresses at integration points of element
   stresses are computed by material models
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 8. 9. 2017
*/
void planeelemqt::res_ip_stresses (long /*lcid*/,long eid)
{
  long i,ii,ipp;
  long ri=0,ci=0;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
	Mm->computenlstresses (ipp,Mm->ip[ipp]);
      
      ipp++;
    }
  }
}



void planeelemqt::res_ip_strains (long lcid,long eid)
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

  ip_strains (lcid,eid,0,0,x,y,r);  
}

/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void planeelemqt::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  double jac;
  vector gp1,gp2,w,eps;
  matrix gm;

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);

    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      geom_matrix_block (gm,ii,x,y,gp1[i],gp2[i],jac);
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,cncomp[ii],eps);
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
void planeelemqt::nod_strains (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,aux,natcoord(2);
  ivector nodes(nne);

  lsm = new double [9];

  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      Mm->givestrain (lcid,ipp,cncomp[ii],eps);
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,eps);

      ipp++;
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->strain_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);
    
    delete [] lhs;  delete [] rhs;
  }
  delete [] lsm;
}

/**
   function computes strains on element
   
   @param lcid - load case id
   @param eid - element id
   
   18.7.2002
*/
void planeelemqt::elem_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,aux,natcoord(2);

  lsm = new double [9];

  nodecoord (nxi,neta);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      Mm->givestrain (lcid,ipp,cncomp[ii],eps);
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,eps);
      
      ipp++;
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    nodal_values (stra,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii]);
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes strains in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi,eta - natural coordinates
   @param fi,li - first and last indices
   @param eps - array containing strains
   
   11.5.2002
*/
void planeelemqt::appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps)
{
  long i,j,k;
  ivector nodes;
  vector nodval;
  
  if (ncomp != eps.n){
    fprintf (stderr,"\n\n wrong interval of indices in function strain (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  reallocv (nne,nodes);
  reallocv (nne,nodval);
  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].strain[lcid*tncomp+i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
}

/**
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void planeelemqt::allip_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector eps,gp1,gp2,w;
  
  reallocv (tncomp,eps);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],w);
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	
	if (Mp->strainaver==0)
	  appval (gp1[i],gp2[i],0,tncomp,eps,stra);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],0,tncomp,eps);
	
	Mm->storestrain (lcid,ipp,eps);
	ipp++;
      }
    }
  }
}

void planeelemqt::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra=NULL;
  vector coord,eps;
  
  if (Mp->strainaver==0){
    stra = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
    }
    elem_strains (stra,lcid,eid,ri,ci);
  }

  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    allip_strains (stra,lcid,eid,ri,ci);
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
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	appval (coord[0],coord[1],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlt::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  if (stra){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }
}




/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta

   10.5.2002
*/
void planeelemqt::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 1.0;  eta[0] = 0.0;
  xi[1] = 0.0;  eta[1] = 1.0;
  xi[2] = 0.0;  eta[2] = 0.0;
  xi[3] = 0.5;  eta[3] = 0.5;
  xi[4] = 0.0;  eta[4] = 0.5;
  xi[5] = 0.5;  eta[5] = 0.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void planeelemqt::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval;
  
  k=0;
  reallocv (nne,nodval);
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
}


/**
   function computes stresses in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void planeelemqt::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector gp1,gp2,w,eps,epst,epstt,sig,auxsig;
  matrix gm,d(tncomp,tncomp),dd;

  reallocm (tncomp,tncomp,d);

  for (ii=0;ii<nb;ii++){
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      Mm->matstiff (d,ipp);
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (ncomp[jj],eps);
	reallocm (ncomp[ii],ncomp[jj],dd);
	
	if (Mp->strainaver==0)
	  Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps);
	
	/*
	if (Mt->elements[eid].presctemp==1){
	  reallocv (tncomp,epstt);
	  tempstrains (lcid,eid,ipp,gp1[i],gp2[i],epstt);
	  reallocv (ncomp[jj],epst);
	  extract (epst,epstt,cncomp[jj],ncomp[jj]);
	  subv (eps,epst,eps);
	}
	*/

	dmatblock (ii,jj,d,dd);
	mxv (dd,eps,auxsig);
	addv (auxsig,sig,sig);
      }
      
      Mm->storestress (lcid,ipp,sig);
      
      ipp++;
    }
  }

}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void planeelemqt::nod_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,epst,epstt,sig,auxsig,natcoord(2);
  ivector nodes(nne);
  matrix d(tncomp,tncomp),dd;
  
  lsm = new double [9];

  Mt->give_elemnodes (eid,nodes);
  
  nodecoord (nxi,neta);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      
      Mm->matstiff (d,ipp);
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (ncomp[jj],eps);
	reallocm (ncomp[ii],ncomp[jj],dd);
	
	if (Mp->strainaver==0)
	  Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps);
	
	/*
	if (Mt->elements[eid].presctemp==1){
	  reallocv (tncomp,epstt);
	  tempstrains (lcid,eid,ipp,gp1[i],gp2[i],epstt);
	  reallocv (ncomp[jj],epst);
	  extract (epst,epstt,cncomp[jj],ncomp[jj]);
	  subv (eps,epst,eps);
	}
	*/
	
	dmatblock (ii,jj,d,dd);
	mxv (dd,eps,auxsig);
	addv (auxsig,sig,sig);
      }
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,sig);

      ipp++;
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->stress_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);
    
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void planeelemqt::elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,epst,epstt,sig,auxsig,natcoord(2);
  matrix d(tncomp,tncomp),dd;
  
  lsm = new double [9];

  nodecoord (nxi,neta);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      
      Mm->matstiff (d,ipp);
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (ncomp[jj],eps);
	reallocm (ncomp[ii],ncomp[jj],dd);

	if (Mp->strainaver==0)
	  appval (gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps,stra);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps);
	
	/*
	if (Mt->elements[eid].presctemp==1){
	  reallocv (tncomp,epstt);
	  tempstrains (lcid,eid,ipp,gp1[i],gp2[i],epstt);
	  reallocv (ncomp[jj],epst);
	  extract (epst,epstt,cncomp[jj],ncomp[jj]);
	  subv (eps,epst,eps);
	}
	*/

	dmatblock (ii,jj,d,dd);
	mxv (dd,eps,auxsig);
	addv (auxsig,sig,sig);
      }
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,sig);

      ipp++;
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    nodal_values (stre,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii]);
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes stresses in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi,eta - natural coordinates
   @param fi,li - first and last indices
   @param sig - array containing stresses

   11.5.2002
*/
void planeelemqt::appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig)
{
  long i,j,k;
  ivector nodes;
  vector nodval;
  
  if (ncomp != sig.n){
    fprintf (stderr,"\n\n wrong interval of indices in function stress (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  reallocv (nne,nodes);
  reallocv (nne,nodval);
  
  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].stress[lcid*tncomp+i];
    }
    sig[k]=approx (xi,eta,nodval);
    k++;
  }
}

/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

   10.5.2002
*/
void planeelemqt::allip_stresses (double **stre,long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector sig(tncomp),gp1,gp2,w;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],w);
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	
	if (Mp->stressaver==0)
	  appval (gp1[i],gp2[i],0,tncomp,sig,stre);
	if (Mp->stressaver==1)
	  appstress (lcid,eid,gp1[i],gp2[i],0,tncomp,sig);
	
	Mm->storestress (lcid,ipp,sig);
	ipp++;
      }
    }
  }
}

void planeelemqt::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra=NULL,**stre = NULL;
  vector coord,sig;
  
  if (Mp->stressaver==0){
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    elem_strains (stra,lcid,eid,ri,ci);
    elem_stresses (stra,stre,lcid,eid,ri,ci);
  }
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    allip_stresses (stre,lcid,eid,ri,ci);
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
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	appval (coord[0],coord[1],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  if (stra){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }
  if (stre){
    for (i=0;i<nne;i++){
      delete [] stre[i];
    }
    delete [] stre;
  }

}









/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void planeelemqt::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void planeelemqt::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void planeelemqt::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemqt::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void planeelemqt::res_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemqt::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemqt::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemqt::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
   function integrates selected quantity over the finite element
   it results in nodal values
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y - node coordinates
   
   TKo, 7.2008
*/
void planeelemqt::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double jac,thick;
  ivector nodes(nne);
  vector w,gp1,gp2,t(nne),ipv,contr(ndofe);
  matrix gm;
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],ipv);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      thick = approx (gp1[i],gp2[i],t);
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
	
      //  strain-displacement (geometric) matrix
      geom_matrix_block (gm,ii,x,y,gp1[i],gp2[i],jac);

      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      cmulv (jac*w[i]*thick,contr);
      
      //  summation
      addv(contr,nv,nv);
      
      ipp++;
    }
  }
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemqt::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  vector w,gp1,gp2,t(nne),eps(tncomp);
  
  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
      
    for (i=0;i<intordsm[ii][ii];i++){
      appval (gp1[i],gp2[i],0,tncomp,eps,stra);
      Mm->storestrain (lcid,ipp,eps);
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemqt::local_values (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  vector w,gp1,gp2,t(nne),eps(tncomp);
  
  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
      
    for (i=0;i<intordsm[ii][ii];i++){
      appval (gp1[i],gp2[i],0,tncomp,eps,stra);
      Mm->storestrain (lcid,ipp,eps);
      //  computation of correct stresses
      if (Mp->strcomp==1)
	Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}



/**
   function computes nonlocal correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemqt::compute_nonloc_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  vector w,gp1,gp2,t(nne),eps(tncomp);
  
  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
      
    for (i=0;i<intordsm[ii][ii];i++){
      appval (gp1[i],gp2[i],0,tncomp,eps,stra);
      Mm->storestrain (lcid,ipp,eps);
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->compnonloc_nlstresses (ipp);
      ipp++;
    }
  }
  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}



/**
   function computes correct stress increments at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemqt::compute_nlstressincr (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  vector w,gp1,gp2,t(nne),eps(tncomp);
  
  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
      
    for (i=0;i<intordsm[ii][ii];i++){
      appval (gp1[i],gp2[i],0,tncomp,eps,stra);
      Mm->storestrain (lcid,ipp,eps);
      //  computation of correct increments of stresses
      if (Mp->strcomp==1)
        Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void planeelemqt::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  for (ii=0;ii<nb;ii++){
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
   The function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param coord - %vector with global coordinates of integration point (ouput)
   
   @return The function returns global coordinates in the argument coord.

   19.1.2002
*/
void planeelemqt::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  vector x(nne),y(nne),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    if (ii==ipp){
      coord[0]=approx (gp1[i],gp2[i],x);
      coord[1]=approx (gp1[i],gp2[i],y);
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
void planeelemqt::ipncoord (long eid,long ipp,vector &ncoord)
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
void planeelemqt::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  vector x(nne),y(nne),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  
  for (long i=0;i<intordsm[ri][ci];i++){
    coord[i][0]=approx (gp1[i],gp2[i],x);
    coord[i][1]=approx (gp1[i],gp2[i],y);
    coord[i][2]=0.0;
  }
}

void planeelemqt::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double xi,eta,jac;
  vector x(nne),y(nne),gp(intordb),w(intordb),av(ndofe),v(ndofe);
  matrix n(napfun,ndofe),am(ndofe,ndofe);
  
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);
  
  if (le[0]==1){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      //xi=(1.0+gp[i])/2.0;  eta=1.0-xi;
      xi=(1.0-gp[i])/2.0;  eta=(1.0+gp[i])/2.0;
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,gp[i],0);
      
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[0]=nv[0]; av[1]=nv[1]; av[6]=nv[2]; av[7]=nv[3]; av[2]=nv[4]; av[3]=nv[5];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      //xi=0.0;  eta=(1.0+gp[i])/2.0;
      xi=0.0;  eta=(1.0-gp[i])/2.0;
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,gp[i],1);
      
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[2]=nv[6]; av[3]=nv[7]; av[8]=nv[8]; av[9]=nv[9]; av[4]=nv[10]; av[5]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=(1.0+gp[i])/2.0;  eta=0.0;

      bf_matrix (n,xi,eta);

      jac1d_2d (jac,x,y,gp[i],2);

      jac*=w[i];

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[4]=nv[12]; av[5]=nv[13]; av[10]=nv[14]; av[11]=nv[15]; av[0]=nv[16]; av[1]=nv[17];
    mxv (am,av,v);  addv (nf,v,nf);
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
void planeelemqt::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
          ipval = approx (xi,eta,anv);
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
void planeelemqt::ipvolume (long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double jac,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),t(nne),gp1,gp2,w;

  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],w);

      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	//  thickness
	thick = approx (gp1[i],gp2[i],t);
	jac_2d (jac,x,y,gp1[i],gp2[i]);
	jac*=w[i]*thick;
	
	Mm->storeipvol (ipp,jac);
	
	ipp++;
      }
    }
  }
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
void planeelemqt::ntdbr_vector (long eid,vector &ntdbr)
{
  long intord = 4;
  long i,k,l,ipp,ri,ci,lcid;
  double thick,jac,**stra;
  ivector nodes(nne);
  vector x(nne),y(nne),gp1(intord),gp2(intord),w(intord);
  vector t(nne),sig(tncomp),eps(tncomp),bf(nne);
  matrix d(tncomp,tncomp);

  lcid = ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];

  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  Mm->matstiff (d,ipp);

  gauss_points_tr (gp1.a,gp2.a,w.a,intord);

  fillv(0.0,ntdbr);

  for (i=0;i<intord;i++){
    appval (gp1[i],gp2[i],0,tncomp,eps,stra);
    mxv (d,eps,sig);

    bf_quad_3_2d (bf.a,gp1[i],gp2[i]);

    thick = approx (gp1[i],gp2[i],t);
    jac_2d (jac,x,y,gp1[i],gp2[i]);
    jac *= w[i]*thick;

    for (k=0;k<tncomp;k++)
      for (l=0;l<nne;l++)
	ntdbr[k*nne+l] += jac * bf[l] * sig[k];

  }

  for (i=0;i<nne;i++)
    delete [] stra[i];
  delete [] stra;
}

/**
   Function integrates function "NT*N" over whole element.
   N is matrix of basic functions.
   !!! Matrix N is composed of three same blocks, it is computed only one third.
   
   @param eid - element id
   @param ntn - empty(returned) 2Darray, dimension is nne x nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemqt::ntn_matrix (long eid,matrix &ntn)
{
  long intord = 6;
  long i,k,l;
  double thick,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),t(nne),gp1(intord),gp2(intord),w(intord),bf(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intord);
  
  fillm(0.0,ntn);
  
  for (i=0;i<intord;i++){
    bf_quad_3_2d (bf.a,gp1[i],gp2[i]);
    
    thick = approx (gp1[i],gp2[i],t);
    jac_2d (jac,x,y,gp1[i],gp2[i]);
    jac *= w[i]*thick;
    
    for (k=0;k<nne;k++)
      for (l=0;l<nne;l++)
        ntn[k][l] += jac * bf[k] * bf[l]; 
    
  }
}

/**
   Function 1.
   1. computes "e2" - square of energy norm of error of solution on element
      e2 = integral{ (sig_star - sig_roof)T * D_inv * (sig_star - sig_roof) }
      sig_star = recovered stress, values in nodes are defined by "rsigfull" and over element are interpolated by base functions
      sig_roof = stress obtained by FEM
      D_inv = inverse stiffness matrix of material
   2. computes "u2" - square of energy norm of strain on element
      u2 = integral{ epsT * D * eps }
      eps = strain obtained by FEM
      D = stiffness matrix of material
   3. computes area of element and adds to "volume" (sum of areas of previous elements)
   4. computes "sizel" (characteristic size of element)
   
   @param eid - element id
   @param volume - sum of areas of previous elements
   @param e2 - empty(returned) value; 
   @param u2 - empty(returned) value; 
   @param sizel - empty(returned) value; 
   @param rsigfull - 1Darray of stress in nodes; dimmension is ncomp*gt->nn after this manner {(val[0]; ... ;val[gt->nn])[0] ; ... ; (val[0]; ... ;val[gt->nn])[ncomp]}
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
double planeelemqt :: compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull)
{
  long intorde = 7;
  long i,ipp,ri,ci,lcid;
  double area,thick,jac,contr,**stra;
  double zero=1.0e-20;
  ivector nodes(nne);
  vector x(nne),y(nne),t(nne),gp1(intorde),gp2(intorde),w(intorde),bf(nne);
  vector sig_star(tncomp),sig_roof(tncomp),sig_err(tncomp),eps(tncomp);
  matrix d(tncomp,tncomp),dinv(tncomp,tncomp);
  
  ri = ci = lcid = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  Mm->matstiff (d,ipp);
  invm (d,dinv,zero);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intorde);
  
  e2 = u2 = 0;
  for (i=0;i<intorde;i++){
    thick = approx (gp1[i],gp2[i],t);
    jac_2d (jac,x,y,gp1[i],gp2[i]);
    jac *= w[i]*thick;
    
    bf_quad_3_2d (bf.a,gp1[i],gp2[i]);
    give_der_star (bf,rsigfull,nodes,sig_star,Mt->nn);
    
    appval (gp1[i],gp2[i],0,tncomp,eps,stra);
    mxv (d,eps,sig_roof);
    
    subv (sig_star,sig_roof,sig_err);
    
    vxmxv (sig_err,dinv,contr);
    e2 += jac*contr;
    
    vxmxv (eps,d,contr);
    u2 += jac*contr;
  }
  
  area = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  sizel = sqrt(1.1547005384*area);
  
  for (i=0;i<nne;i++)  delete [] stra[i];
  delete [] stra;
  
  return area/2.0;
}

/**
   Only for Termit
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemqt::give_sig_roof (matrix &d,double *areacoord,double *x,double *y,long *etnodes,double **xyrr,vector &sig_roof,vector &sig_star)
{                                   //pozor na zmenu poradi ip
  long i,nnecm;
  double xx,yy,detcm,shear,xicm,etacm,auxd;
  vector areacoordcm(3),eps(3),xcm,ycm,rcm,nodsig_x,nodsig_y,nodsig_xy;
  matrix gm;

  nnecm = etnodes[1];

  reallocv (nnecm,xcm);
  reallocv (nnecm,ycm);
  reallocv (nnecm,nodsig_x);
  reallocv (nnecm,nodsig_y);
  reallocv (nnecm,nodsig_xy);
  reallocv (2*nnecm,rcm);
  reallocm (3,2*nnecm,gm);

  for (i=0;i<nnecm;i++){
    xcm[i]       = xyrr[etnodes[i+2]][0];
    ycm[i]       = xyrr[etnodes[i+2]][1];
    rcm[2*i]     = xyrr[etnodes[i+2]][2];
    rcm[2*i+1]   = xyrr[etnodes[i+2]][3];
    nodsig_x[i]  = xyrr[etnodes[i+2]][4];
    nodsig_y[i]  = xyrr[etnodes[i+2]][5];
    nodsig_xy[i] = xyrr[etnodes[i+2]][6];
  }

  xx = yy = 0.0;
  for (i=0;i<3;i++){
    xx += areacoord[i] * x[i];
    yy += areacoord[i] * y[i];
  }

  switch (etnodes[0]){
  case planeelementlt:{
    detcm = (xcm[1]-xcm[0])*(ycm[2]-ycm[0])-(xcm[2]-xcm[0])*(ycm[1]-ycm[0]);
    areacoordcm[0] = ( (xcm[1]*ycm[2] - xcm[2]*ycm[1]) + (ycm[1] - ycm[2])*xx + (xcm[2] - xcm[1])*yy)/detcm;
    areacoordcm[1] = ( (xcm[2]*ycm[0] - xcm[0]*ycm[2]) + (ycm[2] - ycm[0])*xx + (xcm[0] - xcm[2])*yy)/detcm;
    areacoordcm[2] = ( (xcm[0]*ycm[1] - xcm[1]*ycm[0]) + (ycm[0] - ycm[1])*xx + (xcm[1] - xcm[0])*yy)/detcm;  //= 1-ar[0]-ar[1]

    Pelt->geom_matrix (gm,xcm,ycm);
    mxv (gm,rcm,eps);
    mxv (d,eps,sig_roof);

    sig_star[0] = Pelt->approx (areacoordcm,nodsig_x);
    sig_star[1] = Pelt->approx (areacoordcm,nodsig_y);
    sig_star[2] = Pelt->approx (areacoordcm,nodsig_xy);
    break;
  }
  case planeelementlq:{
    nc_lin_4_2d (1e-6,xx,yy,xcm.a,ycm.a,xicm,etacm);
    
    Pelq->geom_matrix (gm,xcm,ycm,0.0,0.0,auxd);
    mxv (gm,rcm,eps);
    mxv (d,eps,sig_roof);
    shear=sig_roof[2];
  
    Pelq->geom_matrix (gm,xcm,ycm,xicm,etacm,auxd);
    mxv (gm,rcm,eps);
    mxv (d,eps,sig_roof);
    sig_roof[2]=shear;

    sig_star[0] = Pelq->approx (xicm,etacm,nodsig_x);
    sig_star[1] = Pelq->approx (xicm,etacm,nodsig_y);
    sig_star[2] = Pelq->approx (xicm,etacm,nodsig_xy);
    break;
  }
  case planeelementqq:{
    nc_lin_4_2d (1e-6,xx,yy,xcm.a,ycm.a,xicm,etacm);   // funguje jen pro subparametriclou sit = rovne strany ctyruhelnika

    Peqq->geom_matrix (gm,xcm,ycm,xicm,etacm,auxd);
    mxv (gm,rcm,eps);
    mxv (d,eps,sig_roof);

    sig_star[0] = Peqq->approx (xicm,etacm,nodsig_x);
    sig_star[1] = Peqq->approx (xicm,etacm,nodsig_y);
    sig_star[2] = Peqq->approx (xicm,etacm,nodsig_xy);
    break;
  }
  case planeelementqt:{
    detcm = (xcm[1]-xcm[0])*(ycm[2]-ycm[0])-(xcm[2]-xcm[0])*(ycm[1]-ycm[0]);
    areacoordcm[0] = ( (xcm[1]*ycm[2] - xcm[2]*ycm[1]) + (ycm[1] - ycm[2])*xx + (xcm[2] - xcm[1])*yy)/detcm;
    areacoordcm[1] = ( (xcm[2]*ycm[0] - xcm[0]*ycm[2]) + (ycm[2] - ycm[0])*xx + (xcm[0] - xcm[2])*yy)/detcm;
    areacoordcm[2] = ( (xcm[0]*ycm[1] - xcm[1]*ycm[0]) + (ycm[0] - ycm[1])*xx + (xcm[1] - xcm[0])*yy)/detcm;

    //Pesqt->geom_matrix (gm,xcm,ycm,areacoordcm);
    mxv (gm,rcm,eps);
    mxv (d,eps,sig_roof);

    sig_star[0] = approx (xicm,etacm,nodsig_x);
    sig_star[1] = approx (xicm,etacm,nodsig_y);
    sig_star[2] = approx (xicm,etacm,nodsig_xy);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function");
    fprintf (stderr," planeelemqt::give_sig_roof (%s, line %d)",__FILE__,__LINE__);
  }
  }
}

///**
//   Only for Termit
//   
//   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
// */
//void planeelemqt::compute_error_fine (long eid,long *etnodes,double **xyrr,vector *rsigfull,double &e2)
//{
//  long intord = 7;
//  long i,ipp,ri,ci;
//  double thick,jac,contr;
//  double zero=1.0e-20;
//  ivector nodes(nne);
//  vector x(nne),y(nne),areacoord(3),t(nne),gp1(intord),gp2(intord),w(intord);
//  vector sig_fine(tncomp),sig_star(tncomp),sig_roof(tncomp),sig_err(tncomp),bf(nne);
//  matrix dinv(tncomp,tncomp),d(tncomp,tncomp);
//
//  ri = ci = 0;
//  ipp = Mt->elements[eid].ipp[ri][ci];
//
//  Mt->give_elemnodes (eid,nodes);
//  Mt->give_node_coord2d (x,y,eid);
//  Mc->give_thickness (eid,nodes,t);
//  
//  gauss_points_tr (gp1.a,gp2.a,w.a,intord);
//
//  Mm->matstiff (d,ipp);
//  invm (d,dinv,zero);
//
//  e2 = 0;
//
//  for (i=0;i<intord;i++){
//    areacoord[0]=gp1[i];  areacoord[1]=gp2[i];  areacoord[2]=1.0-gp1[i]-gp2[i];
//
//    thick = approx (gp1[i],gp2[i],t);
//    jac_2d (jac,x,y,gp1[i],gp2[i]);
//    jac *= w[i]*thick;
//
//    //Mm->givestress (ipp,sig_fine);
//    //ipp++;
//
//
//    //bf_quad_3_2d (bf.a,gp1[i],gp2[i]);
//    ac_bf_quad_3_2d (bf.a,areacoord.a);
//    give_der_star (bf,rsigfull,nodes,sig_fine,Mt->nn);
//    give_sig_roof (d,areacoord.a,x.a,y.a,etnodes,xyrr,sig_roof,sig_star);
//
//    subv (sig_fine,sig_roof,sig_err);
//
//    vxmxv (sig_err,dinv,contr);
//    e2 += jac*contr;
//  }
//}

/**
   This function prepare array 'spder' for using in function `least_square::spr_default`.
   It allocates and fills up array 'spder' by derivatives(strain ro stress) in sampling points.
   For planeelemlt sampling points == main integration point == Gauss's point at the centroid of element.
   
   @param eid   - element id
   @param spder - allocated array 
   @param flags - determines by which derivate is spder filled
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemqt :: elchar (long eid, matrix &spsig)
{
  long intord = 3;
  long i,ipp,ri,ci;
  vector eps(tncomp);
  matrix d(tncomp,tncomp);
  
  reallocm(intord, tncomp, spsig);
  
  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ipp);
  
  //body = Ada->body[2];
  
  //if (body){
    for (i=0; i<intord; i++) {
      eps[0] = Mm->ip[ipp  ].strain[0];
      eps[1] = Mm->ip[ipp  ].strain[1];
      eps[2] = Mm->ip[ipp+3].strain[2];
      ipp++;
      mxv (d.a, eps.a, spsig.a + i*tncomp, d.m, d.n);
    }
  //}
  //else{
  //  long lcid = 0;
  //  double **stra;
  //  vector gp1(intord),gp2(intord),w(intord);
  //  
  //  stra = new double* [nne];
  //  for (i=0;i<nne;i++){
  //    stra[i] = new double [tncomp];
  //  }
  //  elem_strains (stra,lcid,eid,ri,ci);
  //  
  //  gauss_points_tr (gp1.a,gp2.a,w.a,intord);
  //  
  //  for (i=0;i<intord;i++){
  //    appval (gp1[i],gp2[i],0,tncomp,eps,stra);
  //    sig = spsig + i*tncomp;
  //    mxv (d.a,eps.a,sig,d.m,d.n);
  //  }
  //
  //  for (i=0;i<nne;i++)  delete [] stra[i];
  //                       delete [] stra;
  //}
}
////////////////////       /* termitovo */       ////////////////////////////////////
