#include "axisymlq.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "loadcase.h"
#include "intpoints.h"
#include "mechtop.h"
#include "mechmat.h"
#include <math.h>
#include <stdlib.h>


axisymlq::axisymlq (void)
{
  long i,j;
  
  nne=4;  ndofe=8;  tncomp=4;  napfun=2;  ned=4;  nned=2;
  intordmm=3;  intordb=2;
  ssst=axisymm;
  nb=3;

  ncomp = new long [nb];
  ncomp[0]=2;
  ncomp[1]=1;
  ncomp[2]=1;

  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=2;
  cncomp[2]=3;
  
  
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=4;  nip[0][1]=4;  nip[0][2]=0;
  nip[1][0]=4;  nip[1][1]=4;  nip[1][2]=0;
  nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=1;
  
  intordsm[0][0]=2;  intordsm[0][1]=2;  intordsm[0][2]=0;
  intordsm[1][0]=2;  intordsm[1][1]=2;  intordsm[1][2]=0;
  intordsm[2][0]=0;  intordsm[2][1]=0;  intordsm[2][2]=1;

  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

}

axisymlq::~axisymlq (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] intordsm[i];
  }
  delete intordsm;
  
  delete [] cncomp;
  delete [] ncomp;
}

void axisymlq::eleminit (long eid)
{
  long ii,jj;
  Mt->elements[eid].nb=nb;
  Mt->elements[eid].intordsm = new long* [nb];
  Mt->elements[eid].nip = new long* [nb];

  for (ii=0;ii<nb;ii++){
    Mt->elements[eid].intordsm[ii] = new long [nb];
    Mt->elements[eid].nip[ii] = new long [nb];
    for (jj=0;jj<nb;jj++){
      Mt->elements[eid].intordsm[ii][jj]=intordsm[ii][jj];
      Mt->elements[eid].nip[ii][jj]=nip[ii][jj];
    }
  }
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
*/
double axisymlq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function returns matrix of approximation functions
   
   @param n - matrix of approximation functions
   @param xi,eta - natural coordinates
   
   9.7.2001
*/
void axisymlq::bf_matrix (matrix &n,double xi,double eta)
{
  long i,j,k;
  vector bf(nne);
  
  fillm (0.0,n);

  bf_lin_4_2d (bf.a,xi,eta);
  
  j=0;  k=1;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];
    n[1][k]=bf[i];
    j+=2;  k+=2;
  }
}

/**
   function assembles geometric matrix
   
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx
   
   @param gm - geometric matrix
   @param ri - block index
   @param x,y - arrays of node coordinates
   @param xi,eta - natural coordinates
   @param jac - jacobian
   
   8.12.2001
*/
void axisymlq::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  double r;
  vector bf(nne),dx(nne),dy(nne);
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  bf_lin_4_2d (bf.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  r = approx (xi,eta,x);
  if (fabs(r)<Mp->zero){
    fprintf (stderr,"\n\n radius is equal %e in function axisymlq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
  }
  
  fillm (0.0,gm);
  
  i1=0;  i2=1;
  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[1][i2]=dy[i];
    gm[2][i1]=bf[i]/r;
    gm[3][i1]=dy[i];
    gm[3][i2]=dx[i];
    i1+=2;  i2+=2;
  }
}

/**
   function assembles part of geometric matrix
   
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx
   
   @param gm - geometric matrix
   @param ri - block index
   @param x,y - arrays of node coordinates
   @param xi,eta - natural coordinates
   @param jac - jacobian
   
   8.12.2001
*/
void axisymlq::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  double r;
  vector bf(nne),dx(nne),dy(nne);
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  bf_lin_4_2d (bf.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  r = approx (xi,eta,x);
  if (fabs(r)<Mp->zero){
    fprintf (stderr,"\n\n radius is equal %e in function axisymlq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
  }
  
  fillm (0.0,gm);
  
  if (ri==0){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dx[i];
      gm[1][i2]=dy[i];
      i1+=2;  i2+=2;
    }
  }
  if (ri==1){
    i1=0;
    for (i=0;i<nne;i++){
      gm[0][i1]=bf[i]/r;
      i1+=2;
    }
  }
  if (ri==2){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dy[i];
      gm[0][i2]=dx[i];
      i1+=2;  i2+=2;
    }
  }
  
}

/**
   function extracts blocks from stiffness matrix of the material
   
   @param ri,ci - row and column indices
   @param d - stiffness matrix of material
   @param dd - required block from stiffness matrix of material
   
   10.5.2002
*/
void axisymlq::dmatblock (long ri,long ci,matrix &d, matrix &dd)
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
  if (ri==0 && ci==2){
    dd[0][0]=d[0][3];
    dd[1][0]=d[1][3];
  }
  
  if (ri==1 && ci==0){
    dd[0][0]=d[2][0];  dd[0][1]=d[2][1];
  }
  if (ri==1 && ci==1){
    dd[0][0]=d[2][2];
  }
  if (ri==1 && ci==2){
    dd[0][0]=d[2][3];
  }
  
  if (ri==2 && ci==0){
    dd[0][0]=d[3][0];  dd[0][1]=d[3][1];
  }
  if (ri==2 && ci==1){
    dd[0][0]=d[3][2];
  }
  if (ri==2 && ci==2){
    dd[0][0]=d[3][3];
  }
}


/**
   nutno otestovat! pak je mozne smazat tuto hlasku
   
   transformation matrix x_g = T x_l
*/
void axisymlq::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*2][i*2]   = Mt->nodes[nodes[i]].e1[0];    tmat[i*2][i*2+1]   = Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2] = Mt->nodes[nodes[i]].e1[1];    tmat[i*2+1][i*2+1] = Mt->nodes[nodes[i]].e2[1];
    }
  }
}

/**
   function computes stiffness matrix of axisymmetric quadrilateral
   finite element with bilinear approximation functions
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness matrix

   8.12.2001
*/
void axisymlq::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,ii,jj,ipp,transf;
  double xi,eta,jac,r;
  ivector nodes(nne);
  vector x(nne),y(nne),w,gp;
  matrix gmr,gmc,d(tncomp,tncomp),dd;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  
  fillm (0.0,sm);
  

  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gmr);
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  
	  //  geometric matrix
	  geom_matrix_block (gmr,ii,x,y,xi,eta,jac);
	  geom_matrix_block (gmc,jj,x,y,xi,eta,jac);
	  
	  //  matrix of stiffness of the material
	  Mm->matstiff (d,ipp);
	  dmatblock (ii,jj,d,dd);
	  
	  r = approx (xi,eta,x);
	  jac*=w[i]*w[j]*r;
	  
	  //  contribution to the stiffness matrix of the element
	  bdbjac (sm,gmr,dd,gmc,jac);
	  
	  ipp++;
	}
      }
    }
  }
  
  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes resulting stiffness matrix of element
   
   @param eid - element id
   @param sm - stiffness matrix
   
   10.5.2002
*/
void axisymlq::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

/**
   function computes mass matrix of the rectangular axisymmetric
   finite element with bilinear approximation functions
   
   @param eid - number of element
   @param mm - mass matrix

   24.6.2001
*/
void axisymlq::mass_matrix (long eid,matrix &mm)
{
  long i,j;
  double jac,xi,eta,rho,r;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp(intordmm),t(nne),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,mm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      rho = approx (xi,eta,dens);
      r = approx (xi,eta,x);
      jac*=w[i]*w[j]*rho*r;
      
      nnj (mm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

void axisymlq::res_mainip_strains (long lcid,long eid)
{
  long i;
  vector x(nne),y(nne),r(ndofe),aux;
  ivector nodes(nne);
  matrix tmat;

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
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
  
  for (i=0;i<nb;i++){
    mainip_strains (lcid,eid,0,0,i,x,y,r);
  }
}

/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void axisymlq::mainip_strains (long lcid,long eid,long ri,long ci,long ii,vector &x,vector &y,vector &r)
{
  long i,j,ipp;
  double xi,eta,jac;
  vector gp,w,eps;
  matrix gm;
  
  reallocv (intordsm[ii][ii],gp);
  reallocv (intordsm[ii][ii],w);
  reallocv (ncomp[ii],eps);
  reallocm (ncomp[ii],ndofe,gm);
  
  gauss_points (gp.a,w.a,intordsm[ii][ii]);
  
  ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
  for (i=0;i<intordsm[ii][ii];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ii][ii];j++){
      eta=gp[j];
      
      geom_matrix_block (gm,ii,x,y,xi,eta,jac);
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
void axisymlq::nod_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double xi,eta,*lsm,*lhs,*rhs;
  vector x(nne),y(nne),nxi(nne),neta(nne),r(ndofe),gp,w,eps,aux,natcoord(2);
  ivector nodes(nne);
  matrix tmat;

  lsm = new double [9];
  
  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
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

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
	
	natcoord[0]=xi;  natcoord[1]=eta;
	matassem_lsm (lsm,natcoord);
	rhsassem_lsm (rhs,natcoord,eps);
	
	ipp++;
      }
    }
        
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->strain_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);

    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes strains on element
   
   @param val - array containing strains on element
   @param lcid - load case id
   @param eid - element id
   
   15.7.2002
*/
void axisymlq::elem_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double xi,eta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp,w,eps,aux,natcoord(2);

  lsm = new double [9];

  nodecoord (nxi,neta);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
	
	natcoord[0]=xi;  natcoord[1]=eta;
	matassem_lsm (lsm,natcoord);
	rhsassem_lsm (rhs,natcoord,eps);
	
	ipp++;
      }
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
   @param xi, eta - natural coordinates of the point
   @param fi,li - first and last indices
   @param eps - array containing strains
   
   11.5.2002
*/
void axisymlq::appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps)
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
   function computes all strain components at all integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2004
*/
void axisymlq::res_allip_strains (long lcid,long eid)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
  //  all strain components at all integration points
  allip_strains (lcid,eid,0,0);
}

/**
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymlq::allip_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp1,ipp2,ipp3;
  vector eps(tncomp),aux(tncomp);
  
  ipp1=Mt->elements[eid].ipp[ri+0][ci+0];
  ipp2=Mt->elements[eid].ipp[ri+1][ci+1];
  ipp3=Mt->elements[eid].ipp[ri+2][ci+2];

  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      Mm->givestrain (lcid,ipp1,cncomp[0],ncomp[0],eps);
      Mm->givestrain (lcid,ipp2,cncomp[1],ncomp[1],eps);
      Mm->givestrain (lcid,ipp3,cncomp[2],ncomp[2],eps);
      
      Mm->storestrain (lcid,ipp1,eps);
      Mm->storestrain (lcid,ipp2,eps);
      
      addv (aux,eps,aux);

      ipp1++;  ipp2++;
    }
  }

  cmulv(0.25,aux,eps);
  Mm->storestrain (lcid,ipp3,eps);
  
}



void axisymlq::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
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
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlq::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}

/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   
   10.5.2002
*/
void axisymlq::nodecoord (vector &xi,vector &eta)
{
  xi[0] =  1.0;  eta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;
  xi[3] =  1.0;  eta[3] = -1.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void axisymlq::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval(nne);
  
  k=0;
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
}

/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void axisymlq::res_mainip_stresses (long lcid,long eid)
{
  long i;
  
  //  loop over blocks
  for (i=0;i<nb;i++){
    mainip_stresses (lcid,eid,0,0,i);
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
void axisymlq::mainip_stresses (long lcid,long eid,long ri,long ci,long ii)
{
  long i,j,jj,ipp;
  double xi,eta;
  vector gp,w,eps,epst,epstt,sig,auxsig;
  matrix d(tncomp,tncomp),dd;

  reallocv (ncomp[ii],sig);
  reallocv (ncomp[ii],auxsig);
  reallocv (intordsm[ii][ii],gp);
  reallocv (intordsm[ii][ii],w);
  
  gauss_points (gp.a,w.a,intordsm[ii][ii]);
  ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
  
  for (i=0;i<intordsm[ii][ii];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ii][ii];j++){
      eta=gp[j];
      
      Mm->matstiff (d,ipp);
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (ncomp[jj],eps);
	reallocm (ncomp[ii],ncomp[jj],dd);
	
	
	//  block of strains
	Mm->givestrain (lcid,ipp,cncomp[jj],eps);
	
	//  block of stiffness matrix of material
	dmatblock (ii,jj,d,dd);
	//  stress contributions
	mxv (dd,eps,auxsig);
	//  summation of contributions
	addv (auxsig,sig,sig);
      }
      
      //  storage of block of stress
      Mm->storestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
      
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
void axisymlq::nod_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  double xi,eta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),r(ndofe),gp,w,eps,epst,epstt,sig,auxsig,natcoord(2);
  ivector nodes(nne);
  matrix d(tncomp,tncomp),dd;

  lsm = new double [9];

  Mt->give_elemnodes (eid,nodes);
  nodecoord (nxi,neta);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	Mm->matstiff (d,ipp);
	ipp++;
	
	fillv (0.0,sig);
	for (jj=0;jj<nb;jj++){
	  reallocv (ncomp[jj],eps);
	  reallocm (ncomp[ii],ncomp[jj],dd);

	  if (Mp->strainaver==0)
	    Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	  if (Mp->strainaver==1)
	    appstrain (lcid,eid,xi,eta,cncomp[jj],ncomp[jj],eps);
	  
	  /*
	  if (Mt->elements[eid].presctemp==1){
	    reallocv (tncomp,epstt);
	    tempstrains (lcid,eid,ipp,xi,eta,epstt);
	    reallocv (ncomp[jj],epst);
	    extract (epst,epstt,cncomp[jj],ncomp[jj]);
	    subv (eps,epst,eps);
	  }
	  */
	  
	  dmatblock (ii,jj,d,dd);
	  mxv (dd,eps,auxsig);
	  addv (auxsig,sig,sig);
	}
	
	natcoord[0]=xi;  natcoord[1]=eta;
	matassem_lsm (lsm,natcoord);
	rhsassem_lsm (rhs,natcoord,sig);
	
      }
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->stress_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);
    
    
    delete [] lhs;  delete [] rhs;
  }
  delete [] lsm;
}

void axisymlq::elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  double xi,eta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp,w,eps,epst,epstt,sig,auxsig,natcoord(2);
  matrix d(tncomp,tncomp),dd;

  lsm = new double [9];

  nodecoord (nxi,neta);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	Mm->matstiff (d,ipp);
	ipp++;
	
	fillv (0.0,sig);
	for (jj=0;jj<nb;jj++){
	  reallocv (ncomp[jj],eps);
	  reallocm (ncomp[ii],ncomp[jj],dd);

	  if (Mp->strainaver==0)
	    appval (xi,eta,cncomp[jj],ncomp[jj],eps,stra);
	  if (Mp->strainaver==1)
	    appstrain (lcid,eid,xi,eta,cncomp[jj],ncomp[jj],eps);

	  /*
	  if (Mt->elements[eid].presctemp==1){
	    reallocv (tncomp,epstt);
	    tempstrains (lcid,eid,ipp,xi,eta,epstt);
	    reallocv (ncomp[jj],epst);
	    extract (epst,epstt,cncomp[jj],ncomp[jj]);
	    subv (eps,epst,eps);
	  }
	  */

	  dmatblock (ii,jj,d,dd);
	  mxv (dd,eps,auxsig);
	  addv (auxsig,sig,sig);
	}

	natcoord[0]=xi;  natcoord[1]=eta;
	matassem_lsm (lsm,natcoord);
	rhsassem_lsm (rhs,natcoord,sig);
	
      }
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
   @param xi, eta - natural coordinates of the point
   @param fi,li - first and last indices
   @param sig - array containing stresses
   
   11.5.2002
*/
void axisymlq::appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig)
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
   function computes all stress components at all integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void axisymlq::res_allip_stresses (long lcid,long eid)
{
  //  blocks of stress components at integration points
  res_mainip_stresses (lcid,eid);
  //  all stress components at all intergation points
  allip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymlq::allip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp1,ipp2,ipp3;
  vector sig(tncomp),aux(tncomp);
  
  ipp1=Mt->elements[eid].ipp[ri+0][ci+0];
  ipp2=Mt->elements[eid].ipp[ri+1][ci+1];
  ipp3=Mt->elements[eid].ipp[ri+2][ci+2];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      Mm->givestress (lcid,ipp1,cncomp[0],ncomp[0],sig);
      Mm->givestress (lcid,ipp2,cncomp[1],ncomp[1],sig);
      Mm->givestress (lcid,ipp3,cncomp[2],ncomp[2],sig);

      
      Mm->storestress (lcid,ipp1,sig);
      Mm->storestress (lcid,ipp2,sig);
      
      addv (aux,sig,aux);
      ipp1++;  ipp2++;
    }
  }
  
  cmulv(0.25,aux,sig);
  Mm->storestress (lcid,ipp3,sig);
}


void axisymlq::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
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
   function computes load matrix of the axisymmetric quadrilateral
   finite element with bilinear approximation functions
   load vector is obtained after premultiplying load matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load matrix

   8.12.2001
*/
void axisymlq::load_matrix (long eid,matrix &lm)
{
  long i,j;
  double jac,xi,eta,r;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      r = approx (xi,eta,x);
      jac*=r*w[i]*w[j];
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
}





void axisymlq::res_temp_forces (long lcid,long eid,vector &nfor)
{
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  temp_forces (lcid,eid,0,0,nfor,x,y);
}

/**
   function computes nodal forces caused by temperature changes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   22.11.2002, JK
*/
void axisymlq::temp_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
{
  long i,j,k,ii,ipp;
  double xi,eta,jac;
  vector eps,sig,contr(ndofe),epst(tncomp),gp,w;
  matrix d(tncomp,tncomp),dd,gm;
  
  fillv (0.0,nfor);

  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp);
    
    reallocm (ncomp[ii],ndofe,gm);
    reallocm (ncomp[ii],ncomp[ii],dd);
    reallocv (ncomp[ii],eps);
    reallocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	tempstrains (lcid,eid,ipp,xi,eta,epst);
	extract (eps,epst,cncomp[ii],ncomp[ii]);
	
	Mm->matstiff (d,ipp);
	ipp++;
	dmatblock (ii,ii,d,dd);
	mxv (dd,eps,sig);
	geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	mtxv (gm,sig,contr);
	cmulv (jac*w[i]*w[j],contr);
	
	for (k=0;k<contr.n;k++){
	  nfor[k]+=contr[k];
	}
	
      }
    }
  }
}

/**
   function computes strains caused by temperature changes
   
   @param lcid - load case id
   @param eid - element id
   @param ipp - integration point pointer
   @param xi,eta - natural coordinates
   @param eps - array containing strains
   
   22.12.2002, JK
*/
void axisymlq::tempstrains (long lcid,long eid,long ipp,double xi,double eta,vector &eps)
{
  double temp;
  ivector nodes(nne);
  vector dt(nne),tvect(tncomp);
  matrix d(tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mb->lc[lcid].tempchanges (dt.a,nodes);

  temp = approx (xi,eta,dt);
  fillv (temp,tvect);
  
  Mm->matdilat (d,ipp);
  mxv (d,tvect,eps);
  
  Mm->storeeigstrain (ipp,eps);
}







/**
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   8.12.2001
*/
void axisymlq::internal_forces (long lcid,long eid,vector &ifor)
{
  long i,j,k,ii,ipp;
  double xi,eta,jac,rad;
  ivector nodes(nne);
  vector x(nne),y(nne),w,gp;
  vector r(ndofe),eps(tncomp),sig,contr(ndofe),auxcontr(ndofe);
  matrix gm;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  eldispl (0,eid,r.a);
  
  fillv (0.0,ifor);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ii][ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	Mm->computenlstresses (ipp);
	
	Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);

	geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	mtxv (gm,sig,contr);
	
	rad = approx (xi,eta,x);
	cmulv (rad*jac*w[i]*w[j],contr);
	
	for (k=0;k<contr.n;k++){
	  ifor[k]+=contr[k];
	}
	
	ipp++;
	
      }
    }
  }
}

void axisymlq::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,ifor);
}


void axisymlq::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double ww,jac,xi,eta;
  vector x(nne),y(nne),gp(intordb),w(intordb),av(ndofe),v(ndofe);
  matrix n(napfun,ndofe),am(ndofe,ndofe);
  
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);

  if (le[0]==1){
    fillm (0.0,am);
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,0);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[0]=nv[0];  av[1]=nv[1];  av[2]=nv[2];  av[3]=nv[3];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    fillm (0.0,am);
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];

      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,1);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[2]=nv[4];  av[3]=nv[5];  av[4]=nv[6];  av[5]=nv[7];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    fillm (0.0,am);
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,2);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[4]=nv[8];  av[5]=nv[9];  av[6]=nv[10];  av[7]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    fillm (0.0,am);
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,3);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[6]=nv[12];  av[7]=nv[13];  av[0]=nv[14];  av[1]=nv[15];
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
void axisymlq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, ipval;
  vector w, gp, anv(nne);
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
        reallocv (intordsm[ii][jj],gp);
        reallocv (intordsm[ii][jj],w);
        gauss_points (gp.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp[k];
          for (l = 0; l < intordsm[ii][jj]; l++)
          {
            eta=gp[l];
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
