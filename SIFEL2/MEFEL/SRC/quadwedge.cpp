#include <stdlib.h>
#include <math.h>
#include "quadwedge.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "gadaptivity.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"




quadwedge::quadwedge (void)
{
  long i;

  //  number of nodes on element
  nne=15;
  //  number of DOFs on element
  ndofe=45;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=1;
  //  number of edges on element
  ned=9;
  //  number of nodes on one edge
  nned=3;
  //  number of surfaces on element
  nsurf=5;
  //  number of nodes on one surface
  nnsurf=8;
  //  strain/stress state
  ssst=spacestress;
  
  //  number of blocks (parts of geometric matrix)
  nb=1;
  
  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=6;
  
  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  
  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsmt = new long* [nb];
  intordsmz = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsmt[i] = new long [nb];
    intordsmz[i] = new long [nb];
  }
  
  nip[0][0]=9;

  //  total number of integration points
  tnip=9;
  

  intordsmt[0][0]=3;
  intordsmz[0][0]=3;

}

quadwedge::~quadwedge (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsmt[i];
    delete [] intordsmz[i];
  }
  delete [] nip;
  delete [] intordsmt;
  delete [] intordsmz;
  
  delete [] cncomp;
  delete [] ncomp;
}


/**
   function approximates function defined by nodal values
   
   @param xi,eta,zeta - natural coordinates
   @param nodval - nodal values
   
   JK, 19.9.2004
*/
double quadwedge::approx (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_wed_3d (bf.a,xi,eta,zeta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function assembles matrix of base functions

   @param n - matrix of base functions
   @param xi,eta,zeta - natural coordinates
   
   JK, 19.9.2004
*/
void quadwedge::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  long i,j,k,l;
  vector bf(nne);

  fillm (0.0,n);
  
  bf_quad_wed_3d (bf.a,xi,eta,zeta);

  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];  j+=3;
    n[1][k]=bf[i];  k+=3;
    n[2][l]=bf[i];  l+=3;
  }

}

/**
   function assembles geometric matrix
   vector of strains has following ordering
   eps=(e_xx, e_yy, e_zz, e_yz, e_zx, e_xy)
   
   @param gm - geometric matrix
   @param x,y,z - vectors of node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian
   
   JK, 19.9.2004
*/
void quadwedge::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
			    double xi,double eta,double zeta,double &jac)
{
  long i,j,k,l;
  vector dx(nne),dy(nne),dz(nne);
  
  dx_bf_quad_wed_3d (dx.a,xi,eta,zeta);
  dy_bf_quad_wed_3d (dy.a,xi,eta,zeta);
  dz_bf_quad_wed_3d (dz.a,xi,eta,zeta);
  
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
   function assembles transformation matrix
   
   @param nodes - nodes of element
   @param tmat - transformation matrix
   
   NOT PROVED
*/
void quadwedge::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of one element

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   JK, 19.9.2004
*/
void quadwedge::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,k,ii,jj,ipp,transf;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,wt,gp,gp1,gp2;
  matrix gm,d(tncomp,tncomp);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  fillm (0.0,sm);

  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gm);
    for (jj=0;jj<nb;jj++){
      if (intordsmt[ii][jj]==0)  continue;

      reallocv (intordsmz[ii][jj],w);
      reallocv (intordsmz[ii][jj],gp);
      reallocv (intordsmt[ii][jj],wt);
      reallocv (intordsmt[ii][jj],gp1);
      reallocv (intordsmt[ii][jj],gp2);
      
      gauss_points (gp.a,w.a,intordsmz[ii][jj]);
      gauss_points_tr (gp1.a,gp2.a,wt.a,intordsmt[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsmt[ii][jj];i++){
	xi=gp1[i];  eta=gp2[i];
	for (k=0;k<intordsmz[ii][jj];k++){
	  zeta=gp[k];
	  
	  //  geometric matrices
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  
	  Mm->matstiff (d,ipp);  ipp++;
	  
	  jac*=wt[i]*w[k];
	  
	  //  contribution to the stiffness matrix of the element
	  bdbjac (sm,gm,d,gm,jac);
	  
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

   JK, 19.9.2004
*/
void quadwedge::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

