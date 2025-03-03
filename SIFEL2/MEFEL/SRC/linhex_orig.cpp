#include <stdlib.h>
#include <math.h>
#include "linhex.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"




linhex::linhex (void)
{
  long i;

  nne=8;  ndofe=24;  tnip=9;  tncomp=6;  napfun=3;
  intordmm=2;  ned=12;  nned=2;
  ssst=spacestress;
  
  nb=2;
  
  ncomp = new long [nb];
  ncomp[0]=3;
  ncomp[1]=3;

  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=3;
  
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=8;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=8;
  //nip[1][0]=0;  nip[1][1]=1;

  intordsm[0][0]=2;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=2;
  //intordsm[1][0]=0;  intordsm[1][1]=1;

}

linhex::~linhex (void)
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

void linhex::eleminit (long eid)
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
   
   @param xi,eta,zeta - natural coordinates
   @param nodval - nodal values
   
   20.8.2001
*/
double linhex::approx (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_hex_3d (bf.a,xi,eta,zeta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function assembles matrix of base functions
   
   @param n - matrix of base functions
   @param xi,eta,zeta - natural coordinates

   19.7.2001
*/
void linhex::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  long i,j,k,l;
  vector bf(nne);

  fillm (0.0,n);

  bf_lin_hex_3d (bf.a,xi,eta,zeta);
  
  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];  j+=3;
    n[1][k]=bf[i];  k+=3;
    n[2][l]=bf[i];  l+=3;
  }
}

/**
   function assembles geometric matrix

   @param gm - geometric matrix
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - naturalcoordinates
   @param jac - Jacobian

   19.7.2001
*/
void linhex::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac)
{
  long i,j,k,l;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);

  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);

  fillm (0.0,gm);

  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    gm[0][j]=dx[i];
    gm[1][k]=dy[i];
    gm[2][l]=dz[i];
    
    gm[3][k]=dz[i];
    gm[3][l]=dy[i];
    
    gm[4][j]=dz[i];
    gm[4][l]=dx[i];
    
    gm[5][j]=dy[i];
    gm[5][k]=dx[i];
    
    j+=3;  k+=3;  l+=3;
  }
}

/**
   function assembles block of geometric matrix

   @param gm - geometric matrix
   @param bi - block index
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - naturalcoordinates
   @param jac - Jacobian

   19.7.2001
*/
void linhex::geom_matrix_block (matrix &gm,long bi,vector &x,vector &y,vector &z,
				double xi,double eta,double zeta,double &jac)
{
  long i,j,k,l;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);

  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);
  
  fillm (0.0,gm);
  
  if (bi==0){
    j=0;  k=1;  l=2;
    for (i=0;i<nne;i++){
      gm[0][j]=dx[i];
      gm[1][k]=dy[i];
      gm[2][l]=dz[i];
      j+=3;  k+=3;  l+=3;
    }
  }
  if (bi==1){
    j=0;  k=1;  l=2;
    for (i=0;i<nne;i++){
      gm[0][k]=dz[i];
      gm[0][l]=dy[i];
      
      gm[1][j]=dz[i];
      gm[1][l]=dx[i];
      
      gm[2][j]=dy[i];
      gm[2][k]=dx[i];
      
      j+=3;  k+=3;  l+=3;
    }
  }

}

/**
   function extracts blocks from stiffness matrix of the material
   
   @param ri,ci - row and column indices
   @param d - stiffness matrix of the material
   @param dd - required block from stiffness matrix of material
   
   JK, 19.7.2001
*/
void linhex::dmatblock (long ri,long ci,matrix &d,matrix &dd)
{
  fillm (0.0,dd);

  if (ri==0 && ci==0){
    dd[0][0] = d[0][0];  dd[0][1] = d[0][1];  dd[0][2] = d[0][2];
    dd[1][0] = d[1][0];  dd[1][1] = d[1][1];  dd[1][2] = d[1][2];
    dd[2][0] = d[2][0];  dd[2][1] = d[2][1];  dd[2][2] = d[2][2];
  }
  if (ri==1 && ci==1){
    dd[0][0] = d[3][3];
    dd[1][1] = d[4][4];
    dd[2][2] = d[5][5];
  }

}

/**
   function assembles transformation matrix
   
   @param nodes - nodes of element
   @param tmat - transformation matrix
   
*/
void linhex::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness matrix of one element

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness matrix

   19.7.2001
*/
void linhex::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ii,jj,ipp,transf;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  matrix gmr,gmc,dd,d(tncomp,tncomp);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
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
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    //  geometric matrices
	    geom_matrix_block (gmr,ii,x,y,z,xi,eta,zeta,jac);
	    geom_matrix_block (gmc,jj,x,y,z,xi,eta,zeta,jac);
	    
	    Mm->matstiff (d,ipp);  ipp++;
	    dmatblock (ii,jj,d,dd);
	    
	    jac=fabs(jac);
	    
	    jac*=w[i]*w[j]*w[k];
	    
	    //  contribution to the stiffness matrix of the element
	    bdbjac (sm,gmr,dd,gmc,jac);
	    
	  }
	}
      }
    }
  }
  
  //  transformation of stiffness matrix
  ivector nodes (nne);
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
}

/**
   function assembles resulting stiffness matrix of the element
   
   @param eid - element id
   @param sm - stiffness matrix
   
   JK, 9.5.2002
*/
void linhex::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

/**
   function computes mass matrix

   @param eid - number of element
   @param mm - mass matrix

   19.7.2001
*/
void linhex::mass_matrix (long eid,matrix &mm)
{
  long i,j,k;
  double jac,xi,eta,zeta,rho;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);

  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);
	
	jac=fabs(jac);

	bf_matrix (n,xi,eta,zeta);

	rho = approx (xi,eta,zeta,dens);
	
	jac*=w[i]*w[j]*w[k]*rho;
	
	nnj (mm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}

/**
   function computes load matrix

   @param eid - number of element
   @param lm - load matrix
   
   25.7.2001
*/
void linhex::load_matrix (long eid,matrix &lm)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];  w3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);

	bf_matrix (n,xi,eta,zeta);

	jac*=w1*w2*w3;
	
	nnj (lm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}







void linhex::res_mainip_strains (long lcid,long eid)
{
  mainip_strains (lcid,eid,0,0);
}


/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void linhex::mainip_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),r(ndofe),gp,w,eps1,eps2,aux;
  ivector nodes(nne);
  matrix gm1,gm2,tmat;

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
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
    reallocv (ncomp[0],eps1);
    reallocm (ncomp[0],ndofe,gm1);
    reallocv (ncomp[1],eps2);
    reallocm (ncomp[1],ndofe,gm2);

    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  geom_matrix_block (gm1,0,x,y,z,xi,eta,zeta,jac);
	  geom_matrix_block (gm2,1,x,y,z,0.0,0.0,0.0,jac);
	  mxv (gm1,r,eps1);
	  mxv (gm2,r,eps2);
	  
	  Mm->storestrain (lcid,ipp,cncomp[0],ncomp[0],eps1);
	  Mm->storestrain (lcid,ipp,cncomp[1],ncomp[1],eps2);
	  ipp++;
	}
      }
    }
  }
}

/**
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void linhex::nod_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),gp,w,eps,natcoord(3);
  ivector nodes(nne);

  lsm = new double [16];

  nodecoord (nxi,neta,nzeta);
  Mt->give_elemnodes (eid,nodes);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*4];
    rhs = new double [ncomp[ii]*4];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,16);
    nullv (rhs,ncomp[ii]*4);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);

	  natcoord[0]=xi;  natcoord[1]=eta;  natcoord[2]=zeta;
	  matassem_lsm (lsm,natcoord);
	  rhsassem_lsm (rhs,natcoord,eps);
	  
	  ipp++;
	}
      }
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,4,ncomp[ii]);
    Mt->strain_nodal_values (nodes,nxi,neta,nzeta,lhs,3,cncomp[ii],ncomp[ii],lcid);
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes strains on element
   
   @param val - array containing strains on element
   @param lcid - load case id
   @param eid - element id
   
   19.9.2002
*/
void linhex::elem_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),gp,w,eps,natcoord(3);

  lsm = new double [16];

  nodecoord (nxi,neta,nzeta);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    lhs = new double [ncomp[ii]*4];
    rhs = new double [ncomp[ii]*4];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,16);
    nullv (rhs,ncomp[ii]*4);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
	  
	  natcoord[0]=xi;  natcoord[1]=eta;  natcoord[2]=zeta;
	  matassem_lsm (lsm,natcoord);
	  rhsassem_lsm (rhs,natcoord,eps);
	  
	  ipp++;
	}
      }
    }
        
    solve_lsm (lsm,lhs,rhs,Mp->zero,4,ncomp[ii]);
    nodal_values (stra,nxi,neta,nzeta,lhs,3,cncomp[ii],ncomp[ii]);

    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}


/**
   function computes strains in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi, eta, zeta - natural coordinates of the point
   @param fi - first index
   @param ncomp - number of components
   @param eps - array containing strains
   
   11.5.2002
*/
void linhex::appstrain (long lcid,long eid,double xi,double eta,double zeta,long fi,long ncomp,vector &eps)
{
  long i,j,k;
  ivector nod(nne);
  vector nodval(nne);
  
  if (ncomp != eps.n){
    fprintf (stderr,"\n\n wrong interval of indices in function strain (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  Mt->give_elemnodes (eid,nod);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nod[j]].strain[lcid*tncomp+i];
    }
    eps[k]=approx (xi,eta,zeta,nodval);
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
void linhex::allip_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta;
  vector eps(tncomp),gp,w;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    if (Mp->strainaver==0)
	      appval (xi,eta,zeta,0,tncomp,eps,stra);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,xi,eta,zeta,0,tncomp,eps);
	    
	    Mm->storestrain (lcid,ipp,eps);
	    ipp++;
	  }
	}
      }
    }
  }
}

void linhex::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra;
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
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	appstrain (lcid,eid,coord[0],coord[1],coord[2],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlq::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  if (Mp->strainaver==0){
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
   @param zeta - array containing natrual coordinates zeta
   
   10.5.2002
*/
void linhex::nodecoord (vector &xi,vector &eta,vector &zeta)
{
  xi[0] =  1.0;  eta[0] =  1.0;  zeta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;  zeta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;  zeta[2] =  1.0;
  xi[3] =  1.0;  eta[3] = -1.0;  zeta[3] =  1.0;
  xi[4] =  1.0;  eta[4] =  1.0;  zeta[4] = -1.0;
  xi[5] = -1.0;  eta[5] =  1.0;  zeta[5] = -1.0;
  xi[6] = -1.0;  eta[6] = -1.0;  zeta[6] = -1.0;
  xi[7] =  1.0;  eta[7] = -1.0;  zeta[7] = -1.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void linhex::appval (double xi,double eta,double zeta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval(nne);
  
  k=0;
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (xi,eta,zeta,nodval);
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
void linhex::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta;
  vector gp,w,eps,epst,epstt,sig,auxsig;
  matrix d(tncomp,tncomp),dd;

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
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
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->matstiff (d,ipp);
	  
	  fillv (0.0,sig);
	  for (jj=0;jj<nb;jj++){
	    reallocv (ncomp[jj],eps);
	    reallocm (ncomp[ii],ncomp[jj],dd);

	    if (Mp->strainaver==0)
	      Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,xi,eta,zeta,cncomp[jj],ncomp[jj],eps);
	    
	    /*
	    if (Mt->elements[eid].presctemp==1){
	      reallocv (tncomp,epstt);
	      tempstrains (lcid,eid,ipp,xi,eta,zeta,epstt);
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
  }
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void linhex::nod_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),r(ndofe),gp,w,eps,epst,epstt,sig,auxsig,natcoord(3);
  ivector nodes(nne);
  matrix d(tncomp,tncomp),dd;

  lsm = new double [16];

  nodecoord (nxi,neta,nzeta);
  Mt->give_elemnodes (eid,nodes);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*4];
    rhs = new double [ncomp[ii]*4];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,16);
    nullv (rhs,ncomp[ii]*4);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->matstiff (d,ipp);
	  
	  fillv (0.0,sig);
	  for (jj=0;jj<nb;jj++){
	    reallocv (ncomp[jj],eps);
	    reallocm (ncomp[ii],ncomp[jj],dd);
	    
	    if (Mp->strainaver==0)
	      Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,xi,eta,zeta,cncomp[jj],ncomp[jj],eps);
	    
	    /*
	    if (Mt->elements[eid].presctemp==1){
	      reallocv (tncomp,epstt);
	      tempstrains (lcid,eid,ipp,xi,eta,zeta,epstt);
	      reallocv (ncomp[jj],epst);
	      extract (epst,epstt,cncomp[jj],ncomp[jj]);
	      subv (eps,epst,eps);
	    }
	    */

	    dmatblock (ii,jj,d,dd);
	    mxv (dd,eps,auxsig);
	    addv (auxsig,sig,sig);
	  }
	  
	  natcoord[0]=xi;  natcoord[1]=eta;  natcoord[2]=zeta;
	  matassem_lsm (lsm,natcoord);
	  rhsassem_lsm (rhs,natcoord,sig);

	  ipp++;
	}
      }
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,4,ncomp[ii]);
    Mt->stress_nodal_values (nodes,nxi,neta,nzeta,lhs,3,cncomp[ii],ncomp[ii],lcid);
    
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

void linhex::elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),r,gp,w,eps,epst,epstt,sig,auxsig,natcoord(3);
  matrix d(tncomp,tncomp),dd;

  lsm = new double [16];

  nodecoord (nxi,neta,nzeta);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    lhs = new double [ncomp[ii]*4];
    rhs = new double [ncomp[ii]*4];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,16);
    nullv (rhs,ncomp[ii]*4);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->matstiff (d,ipp);
	  
	  fillv (0.0,sig);
	  for (jj=0;jj<nb;jj++){
	    reallocv (ncomp[jj],eps);
	    reallocm (ncomp[ii],ncomp[jj],dd);
	    
	    if (Mp->strainaver==0)
	      appval (xi,eta,zeta,cncomp[jj],ncomp[jj],eps,stra);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,xi,eta,zeta,cncomp[jj],ncomp[jj],eps);
	    
	    /*
	    if (Mt->elements[eid].presctemp==1){
	      reallocv (tncomp,epstt);
	      tempstrains (lcid,eid,ipp,xi,eta,zeta,epstt);
	      reallocv (ncomp[jj],epst);
	      extract (epst,epstt,cncomp[jj],ncomp[jj]);
	      subv (eps,epst,eps);
	    }
	    */

	    dmatblock (ii,jj,d,dd);
	    mxv (dd,eps,auxsig);
	    addv (auxsig,sig,sig);
	  }
	  
	  natcoord[0]=xi;  natcoord[1]=eta;  natcoord[2]=zeta;
	  matassem_lsm (lsm,natcoord);
	  rhsassem_lsm (rhs,natcoord,sig);

	  ipp++;
	}
      }
    }
        
    solve_lsm (lsm,lhs,rhs,Mp->zero,4,ncomp[ii]);
    nodal_values (stre,nxi,neta,nzeta,lhs,3,cncomp[ii],ncomp[ii]);


    delete [] lhs;  delete [] rhs;
  }

  delete [] lsm;
}


/**
   function computes stresses in arbitrary point on element
   
   @param lcid - load case id
   @param eid - element id
   @param xi, eta, zeta - natural coordinates of the point
   @param fi,li - first and last indices
   @param sig - array containing stresses
   
   11.5.2002
*/
void linhex::appstress (long lcid,long eid,double xi,double eta,double zeta,long fi,long ncomp,vector &sig)
{
  long i,j,k;
  ivector nodes(nne);
  vector nodval(nne);
  
  if (ncomp != sig.n){
    fprintf (stderr,"\n\n wrong interval of indices in function stress (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].stress[lcid*tncomp+i];
    }
    sig[k]=approx (xi,eta,zeta,nodval);
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
void linhex::allip_stresses (double **stre,long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta;
  vector sig(tncomp),gp,w;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    if (Mp->stressaver==0)
	      appval (xi,eta,zeta,0,tncomp,sig,stre);
	    if (Mp->stressaver==1)
	      appstress (lcid,eid,xi,eta,zeta,0,tncomp,sig);

	    Mm->storestress (lcid,ipp,sig);
	    ipp++;
	  }
	}
      }
    }
  }
}







void linhex::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra,**stre;
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
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],coord[2],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  if (Mp->stressaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
      delete [] stre[i];
    }
    delete [] stra;
    delete [] stre;
  }
}











void linhex::res_temp_forces (long lcid,long eid,vector &nfor)
{
  temp_forces (lcid,eid,0,0,nfor);
}

/**
   function computes nodal forces caused by temperature changes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   
   30.11.2002, JK
*/
void linhex::temp_forces (long lcid,long eid,long ri,long ci,vector &nfor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),eps,sig,contr(ndofe),t(nne),epst(tncomp),gp,w;
  matrix d(tncomp,tncomp),dd,gm;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  
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
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  tempstrains (lcid,eid,ipp,xi,eta,zeta,epst);
	  extract (eps,epst,cncomp[ii],ncomp[ii]);
	  
	  Mm->matstiff (d,ipp);
	  ipp++;
	  dmatblock (ii,ii,d,dd);
	  mxv (dd,eps,sig);
	  geom_matrix_block (gm,ii,x,y,z,xi,eta,zeta,jac);
	  mtxv (gm,sig,contr);
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  for (l=0;l<contr.n;l++){
	    nfor[l]+=contr[l];
	  }
	}
      }
    }
  }
}

/**
   function computes nodal forces caused by temperature changes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   
   30.11.2002, JK
*/
void linhex::temperaturestrains (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta;
  vector epst(tncomp),gp,w;
  
  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];

	  tempstrains (lcid,eid,ipp,xi,eta,zeta,epst);

	  ipp++;
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
   @param xi,eta,zeta - natural coordinates
   @param eps - array containing strains
   
   30.11.2002, JK
*/
void linhex::tempstrains (long lcid,long eid,long ipp,double xi,double eta,double zeta,vector &eps)
{
  double temp;
  ivector nodes(nne);
  vector dt(nne),tvect(tncomp);
  matrix d(tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mb->lc[lcid].tempchanges (dt.a,nodes);

  temp = approx (xi,eta,zeta,dt);
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
   
   28.7.2001
*/
void linhex::internal_forces2 (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  double **stra;
  vector x(nne),y(nne),z(nne),w,gp,r(ndofe),eps(tncomp),sig,contr(ndofe);
  matrix gm;

  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  Mt->give_node_coord3d (x,y,z,eid);
  
  fillv (0.0,ifor);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  appval (xi,eta,zeta,0,tncomp,eps,stra);
	  
	  Mm->storestrain (lcid,ipp,eps);
	  
	  Mm->computenlstresses (ipp);
	  
	  Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
	  
	  geom_matrix_block (gm,ii,x,y,z,xi,eta,zeta,jac);
	  mtxv (gm,sig,contr);
	  
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  for (l=0;l<contr.n;l++){
	    ifor[l]+=contr[l];
	  }
	  
	  ipp++;
	}
      }
    }
  }

  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}

/**
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   28.7.2001
*/
void linhex::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w,gp,r(ndofe),eps1(3),eps2(3),eps(tncomp),sig,contr(ndofe);
  matrix gm,gm1(3,ndofe),gm2(3,ndofe);

  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,r.a);
  
  fillv (0.0,ifor);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  geom_matrix_block (gm1,0,x,y,z,xi,eta,zeta,jac);
	  //geom_matrix_block (gm2,1,x,y,z,0.0,0.0,0.0,jac);
	  geom_matrix_block (gm2,1,x,y,z,xi,eta,zeta,jac);
	  
	  mxv (gm1,r,eps1);
	  mxv (gm2,r,eps2);
	  
	  eps[0]=eps1[0];
	  eps[1]=eps1[1];
	  eps[2]=eps1[2];
	  eps[3]=eps2[0];
	  eps[4]=eps2[1];
	  eps[5]=eps2[2];
	  
	  Mm->storestrain (lcid,ipp,eps);
	  
	  Mm->computenlstresses (ipp);
	  
	  Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
	  
	  geom_matrix_block (gm,ii,x,y,z,xi,eta,zeta,jac);
	  mtxv (gm,sig,contr);
	  
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  for (l=0;l<contr.n;l++){
	    ifor[l]+=contr[l];
	  }
	  
	  ipp++;
	}
      }
    }
  }

}

void linhex::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,0,0,ifor);
}

/**
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   28.7.2001
*/
void linhex::local_values (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta;
  double **stra;
  vector w,gp,eps(tncomp);
  matrix gm;

  stra = new double* [nne];
  for (i=0;i<nne;i++){
    stra[i] = new double [tncomp];
  }
  elem_strains (stra,lcid,eid,ri,ci);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  appval (xi,eta,zeta,0,tncomp,eps,stra);
	  
	  Mm->storestrain (lcid,ipp,eps);
	  
	  Mm->computenlstresses (ipp);
	  
	  ipp++;
	}
      }
    }
  }

  for (i=0;i<nne;i++){
    delete [] stra[i];
  }
  delete [] stra;
}


/**
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   28.7.2001
*/
void linhex::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp,sig,contr(ndofe);
  matrix gm;


  Mt->give_node_coord3d (x,y,z,eid);
  
  fillv (0.0,ifor);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  Mm->compnonloc_nlstresses (ipp);
	  
	  Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
	  
	  geom_matrix_block (gm,ii,x,y,z,xi,eta,zeta,jac);
	  mtxv (gm,sig,contr);
	  
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  for (l=0;l<contr.n;l++){
	    ifor[l]+=contr[l];
	  }
	  
	  ipp++;
	}
      }
    }
  }
}


/**
   function returns coordinates of integration points
   
   @param eid - element id
   @param ipp - integration point pointer
   @param coord - vector of coordinates
   
   10.1.2002
*/
void linhex::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,k,ii;
  double xi,eta,zeta;
  vector x(nne),y(nne),z(nne),w(intordsm[ri][ci]),gp(intordsm[ri][ci]);
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordsm[ri][ci];k++){
	zeta=gp[k];
	if (ii==ipp){
	  coord[0]=approx (xi,eta,zeta,x);
	  coord[1]=approx (xi,eta,zeta,y);
	  coord[2]=approx (xi,eta,zeta,z);
	}
	ii++;
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
void linhex::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, m, ipp;
  long ii, jj, nv = nodval.n;
  long nstra;
  double xi, eta, zeta, ipval;
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
    for(i = 0; i < nne; i++) // for all nodes on element
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
            for (m = 0; m < intordsm[ii][jj]; m++)
            {
              zeta=gp[m];
              //  value in integration point
              ipval = approx (xi,eta,zeta,anv);
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
*/
void linhex::ipvolume (long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    jac_3d (jac,x,y,z,xi,eta,zeta);
	    jac=fabs(jac);
	    
	    jac*=w[i]*w[j]*w[k];
	    
	    Mm->storeipvol (ipp,jac);
	    ipp++;
	  }
	}
      }
    }
  }
}
