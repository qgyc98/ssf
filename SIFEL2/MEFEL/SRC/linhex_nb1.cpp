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
  long i,j;

  nne=8;  ndofe=24;  tncomp=6;  napfun=3;
  intordmm=2;  ned=12;  nned=2; intordb=2;  nsurf=6;  nnsurf=4;
  ssst=spacestress;
  
  nb=1;
  
  ncomp = new long [nb];
  ncomp[0]=6;

  cncomp = new long [nb];
  cncomp[0]=0;

  
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=8;

  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=2;

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
   function assembles strain-displacement (geometric) %matrix

   @param gm - geometric %matrix
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - natural coordinates
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
   function assembles geometric matrix

   @param gm - geometric matrix
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - naturalcoordinates
   @param jac - Jacobian

   19.7.2001
*/
void linhex::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,vector &z,
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
   function assembles auxiliary vectors B for evaluation of stiffness %matrix
   in geometrically nonlinear problems
   
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian
   @param b11,b12,b13,b21,b22,b23,b31,b32,b33 - vectors of derivatives of shape functions
   
   JK, 24.9.2005
*/
void linhex::bvectors (vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac,
		       vector &b11,vector &b12,vector &b13,
		       vector &b21,vector &b22,vector &b23,
		       vector &b31,vector &b32,vector &b33)
{
  vector dx(nne),dy(nne),dz(nne);
  
  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);
  
  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);
  
  fillv (0.0,b11);
  fillv (0.0,b12);
  fillv (0.0,b13);
  fillv (0.0,b21);
  fillv (0.0,b22);
  fillv (0.0,b23);
  fillv (0.0,b31);
  fillv (0.0,b32);
  fillv (0.0,b33);

  //  du/dx
  b11[0]=dx[0];  b11[3]=dx[1];  b11[6]=dx[2];  b11[9]=dx[3];   b11[12]=dx[4];  b11[15]=dx[5];  b11[18]=dx[6];  b11[21]=dx[7];
  //  du/dy
  b12[0]=dy[0];  b12[3]=dy[1];  b12[6]=dy[2];  b12[9]=dy[3];   b12[12]=dy[4];  b12[15]=dy[5];  b12[18]=dy[6];  b12[21]=dy[7];
  //  du/dz
  b13[0]=dz[0];  b13[3]=dz[1];  b13[6]=dz[2];  b13[9]=dz[3];   b13[12]=dz[4];  b13[15]=dz[5];  b13[18]=dz[6];  b13[21]=dz[7];

  //  dv/dx
  b21[1]=dx[0];  b21[4]=dx[1];  b21[7]=dx[2];  b21[10]=dx[3];  b21[13]=dx[4];  b21[16]=dx[5];  b21[19]=dx[6];  b21[22]=dx[7];
  //  dv/dy
  b22[1]=dy[0];  b22[4]=dy[1];  b22[7]=dy[2];  b22[10]=dy[3];  b22[13]=dy[4];  b22[16]=dy[5];  b22[19]=dy[6];  b22[22]=dy[7];
  //  dv/dz
  b23[1]=dz[0];  b23[4]=dz[1];  b23[7]=dz[2];  b23[10]=dz[3];  b23[13]=dz[4];  b23[16]=dz[5];  b23[19]=dz[6];  b23[22]=dz[7];

  //  dw/dx
  b31[2]=dx[0];  b31[5]=dx[1];  b31[8]=dx[2];  b31[11]=dx[3];  b31[14]=dx[4];  b31[17]=dx[5];  b31[20]=dx[6];  b31[23]=dx[7];
  //  dw/dy
  b32[2]=dy[0];  b32[5]=dy[1];  b32[8]=dy[2];  b32[11]=dy[3];  b32[14]=dy[4];  b32[17]=dy[5];  b32[20]=dy[6];  b32[23]=dy[7];
  //  dw/dz
  b33[2]=dz[0];  b33[5]=dz[1];  b33[8]=dz[2];  b33[11]=dz[3];  b33[14]=dz[4];  b33[17]=dz[5];  b33[20]=dz[6];  b33[23]=dz[7];

}

/**
   function computes strain-displacement %matrix for geometrically nonlinear problems
   
   @param gm - strain-displacement %matrix
   @param r - array of nodal displacements
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian

   JK, 24.9.2005
*/
void linhex::gngeom_matrix (matrix &gm,vector &r,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  double b11r,b12r,b13r,b21r,b22r,b23r,b31r,b32r,b33r;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe),av(ndofe);
  
  fillm (0.0,gm);
  
  bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
  
  scprd (b11,r,b11r);
  scprd (b12,r,b12r);
  scprd (b13,r,b13r);
  scprd (b21,r,b21r);
  scprd (b22,r,b22r);
  scprd (b23,r,b23r);
  scprd (b31,r,b31r);
  scprd (b32,r,b32r);
  scprd (b33,r,b33r);
  
  
  // *******
  //  E_11
  // *******
  
  //  B11 dr
  for (i=0;i<ndofe;i++){
    gm[0][i]+=b11[i];
  }
  
  //  r B11 B11 dr
  cmulv(b11r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }
  
  //  r B21 B21 dr
  cmulv(b21r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }

  //  r B31 B31 dr
  cmulv(b31r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }

  
  // *******
  //  E_22
  // *******
  
  //  B22 dr
  for (i=0;i<ndofe;i++){
    gm[1][i]+=b22[i];
  }
  
  //  r B12 B12 dr
  cmulv(b12r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  //  r B22 B22 dr
  cmulv(b22r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  //  r B32 B32 dr
  cmulv(b32r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  // *******
  //  E_33
  // *******
  
  //  B33 dr
  for (i=0;i<ndofe;i++){
    gm[2][i]+=b22[i];
  }
  
  //  r B13 B13 dr
  cmulv(b13r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B23 B23 dr
  cmulv(b23r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B33 B33 dr
  cmulv(b33r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  

  // **************
  //  E_23 = E_32
  // **************
  
  //  (B23 + B32) dr
  for (i=0;i<ndofe;i++){
    gm[3][i]+=b23[i]+b32[i];
  }
  
  //  r B13 B12 dr
  cmulv(b13r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B12 B13 dr
  cmulv(b12r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }

  //  r B23 B22 dr
  cmulv(b23r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B22 B23 dr
  cmulv(b22r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B33 B32 dr
  cmulv(b33r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  //  r B32 B33 dr
  cmulv(b32r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[3][i]+=av[i];
  }
  
  // **************
  //  E_31 = E_13
  // **************
  
  //  (B31 + B13) dr
  for (i=0;i<ndofe;i++){
    gm[4][i]+=b31[i]+b13[i];
  }
  
  //  r B11 B13 dr
  cmulv(b11r,b13,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B13 B11 dr
  cmulv(b13r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }

  //  r B21 B23 dr
  cmulv(b21r,b23,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B23 B21 dr
  cmulv(b23r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B31 B33 dr
  cmulv(b31r,b33,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }
  
  //  r B33 B31 dr
  cmulv(b33r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[4][i]+=av[i];
  }


  // **************
  //  E_12 = E_21
  // **************
  
  //  (B12 + B21) dr
  for (i=0;i<ndofe;i++){
    gm[5][i]+=b12[i]+b21[i];
  }
  
  //  r B12 B11 dr
  cmulv(b12r,b11,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B11 B12 dr
  cmulv(b11r,b12,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }

  //  r B22 B21 dr
  cmulv(b22r,b21,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B21 B22 dr
  cmulv(b21r,b22,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B32 B31 dr
  cmulv(b32r,b31,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  
  //  r B31 B32 dr
  cmulv(b31r,b32,av);
  for (i=0;i<ndofe;i++){
    gm[5][i]+=av[i];
  }
  


}


/**
   function computes gradient %matrix for geometrically nonlinear problems
   
   @param grm - gradient %matrix
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coordinates
   @param jac - Jacobian

   JK, 24.9.2005
*/
void linhex::gnl_grmatrix (matrix &grm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe);
  
  bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
  
  for (i=0;i<ndofe;i++){
    grm[0][i]=b11[i];
    grm[1][i]=b12[i];
    grm[2][i]=b13[i];
    grm[3][i]=b21[i];
    grm[4][i]=b22[i];
    grm[5][i]=b23[i];
    grm[6][i]=b31[i];
    grm[7][i]=b32[i];
    grm[8][i]=b33[i];
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
   function computes stiffness %matrix of one element

   function computes stiffness %matrix for geometrically linear problems

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness matrix

   JK, 19.7.2001
*/
void linhex::gl_stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ii,jj,ipp,transf;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  matrix gm,d(tncomp,tncomp);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  fillm (0.0,sm);

  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gm);
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
	    
	    //  geometric matrices
	    geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	    
	    Mm->matstiff (d,ipp);  ipp++;
	    
	    jac*=w[i]*w[j]*w[k];
	    
	    //  contribution to the stiffness matrix of the element
	    bdbjac (sm,gm,d,gm,jac);
	    
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
   function computes stiffness %matrix of quadrilateral finite element
   
   function computes stiffness %matrix for geometrically nonlinear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of nodal coordinates
   
   JK, 21.9.2005
*/
void linhex::gnl_stiffness_matrix (long lcid,long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac,jac2;
  vector w,gp,sig(tncomp),r(ndofe),x(nne),y(nne),z(nne);
  matrix gm(tncomp,ndofe),grm(9,ndofe),d(tncomp,tncomp),s(9,9);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a,cn.a,ndofe);
  
  //  component setting to zero
  fillm (0.0,sm);
  
  //  array for weights of integration points
  reallocv (intordsm[0][0],w);
  //  array for coordinates of integration points
  reallocv (intordsm[0][0],gp);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	//
	//  linear stiffness matrix and inital deformation matrix
	//
	
	//  strain-displacement matrix
	gngeom_matrix (gm,r,x,y,z,xi,eta,zeta,jac);
	
	//  stiffness matrix of the material
	Mm->matstiff (d,ipp);
	
	jac*=w[i]*w[j]*w[k];
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,gm,d,gm,jac);
	
	
	//
	//  initial stress matrix
	//
	
	//  gradient matrix
	gnl_grmatrix (grm,x,y,z,xi,eta,zeta,jac2);
	
	//  stresses
	Mm->givestress (lcid,ipp,sig);
	
	s[0][0]=sig[0];  s[0][1]=sig[5];  s[0][2]=sig[4];
	s[1][0]=sig[5];  s[1][1]=sig[1];  s[1][2]=sig[3];
	s[2][0]=sig[4];  s[2][1]=sig[3];  s[2][2]=sig[2];
	
	s[3][3]=sig[0];  s[3][4]=sig[5];  s[3][5]=sig[4];
	s[4][3]=sig[5];  s[4][4]=sig[1];  s[4][5]=sig[3];
	s[5][3]=sig[4];  s[5][4]=sig[3];  s[5][5]=sig[2];
	
	s[6][6]=sig[0];  s[6][7]=sig[5];  s[6][8]=sig[4];
	s[7][6]=sig[5];  s[7][7]=sig[1];  s[7][8]=sig[3];
	s[8][6]=sig[4];  s[8][7]=sig[3];  s[8][8]=sig[2];
	
	
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,grm,s,grm,jac);
	
	ipp++;
      }
    }
  }
}


/**
   function assembles resulting stiffness matrix of the element
   
   @param lcid - load case id
   @param eid - element id
   @param sm - stiffness matrix
   
   JK, 9.5.2002
*/
void linhex::res_stiffness_matrix (long lcid,long eid,matrix &sm)
{
  gl_stiffness_matrix (eid,0,0,sm);
  //gnl_stiffness_matrix (lcid,eid,0,0,sm);
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
	
	//jac=fabs(jac);
	//if (jac<0.0)  fprintf (stderr,"\n zaporny jakobian ve funkci load_matrix");

	nnj (lm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}





/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void linhex::res_mainip_strains (long lcid,long eid)
{
  long i;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),r(ndofe),aux;
  matrix tmat;
  
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
  
  gl_mainip_strains (lcid,eid,0,0,x,y,z,r);
  //gnl_mainip_strains (lcid,eid,0,0,x,y,z,r);
  
}


/**
   function computes block of strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y,z - %vectors of nodal coordinates
   @param r - %vector of nodal displacements
   
   10.5.2002
*/
void linhex::gl_mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,jac;
  vector gp,w,eps;
  matrix gm;
  
  for (ii=0;ii<nb;ii++){
    
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
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  
	  mxv (gm,r,eps);
	  
	  Mm->storestrain (lcid,ipp,cncomp[ii],eps);
	  
	  ipp++;
	}
      }
    }
  }
}

/**
   function computes strains at integration points of element
   function is used in geometrically linear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ii - number of block
   @param x,y,z - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 24.9.2005
*/
void linhex::gnl_mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac,b11r,b12r,b13r,b21r,b22r,b23r,b31r,b32r,b33r;
  vector gp,w,eps;
  vector b11(ndofe),b12(ndofe),b13(ndofe),b21(ndofe),b22(ndofe),b23(ndofe),b31(ndofe),b32(ndofe),b33(ndofe);
  
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  reallocv (tncomp,eps);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	bvectors (x,y,z,xi,eta,zeta,jac,b11,b12,b13,b21,b22,b23,b31,b32,b33);
	
	scprd (b11,r,b11r);
	scprd (b12,r,b12r);
	scprd (b13,r,b13r);
	scprd (b21,r,b21r);
	scprd (b22,r,b22r);
	scprd (b23,r,b23r);
	scprd (b31,r,b31r);
	scprd (b32,r,b32r);
	scprd (b33,r,b33r);
	
	eps[0] = b11r + 0.5*b11r*b11r + 0.5*b21r*b21r + 0.5*b31r*b31r;
	eps[1] = b22r + 0.5*b12r*b12r + 0.5*b22r*b22r + 0.5*b32r*b32r;
	eps[2] = b33r + 0.5*b13r*b13r + 0.5*b23r*b23r + 0.5*b33r*b33r;

	eps[3] = b23r+b32r + b12r*b13r + b22r*b23r + b32r*b33r;
	eps[4] = b31r+b13r + b13r*b11r + b23r*b21r + b33r*b31r;
	eps[5] = b12r+b21r + b11r*b12r + b21r*b22r + b31r*b32r;
	
	
	Mm->storestrain (lcid,ipp,eps);
	ipp++;
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
void linhex::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ri,ci,ipnum);
  
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
   @param stra - array for strain components
   
   JK, 26.9.2004
*/
void linhex::nod_strains_comp (long lcid,long eid,double **stra)
{
  long i,j;
  double jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),nxi(nne),neta(nne),nzeta(nne),r(ndofe),eps(tncomp),aux;
  matrix tmat,gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
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
    locglobtransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  nodecoord (nxi,neta,nzeta);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,z,nxi[i],neta[i],nzeta[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    
    for (j=0;j<eps.n;j++){
      stra[i][j]=eps[j];
    }
  }

}

/**
   function computes strains on element
   
   @param val - array containing strains on element
   @param lcid - load case id
   @param eid - element id
   
   19.9.2002
*/
/*  zruseno 26.9.2004
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
*/


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
   function computes strains at all integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2004
*/
void linhex::res_allip_strains (long lcid,long eid)
{
  allip_strains (lcid,eid,0,0);
}

/**
   function computes strains at all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void linhex::allip_strains (long lcid,long eid,long ri,long ci)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
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
    //elem_strains (stra,lcid,eid,ri,ci);
  }
  
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
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (ncp,eps);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
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
   function returns numbers of integration point closest to element nodes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipnum - array of numbers
   
   JK, 26.9.2004
*/
void linhex::nodipnum (long eid,long ri,long ci,ivector &ipnum)
{
  long i,j;
  
  j=intordsm[0][0];
  i=Mt->elements[eid].ipp[ri][ci];

  ipnum[0]=i+j*j*j-1;
  ipnum[1]=i+j*j-1;
  ipnum[2]=i+j-1;
  ipnum[3]=i+j*j*(j-1)+j-1;
  ipnum[4]=i+j*j*(j-1)+j*(j-1);
  ipnum[5]=i+j*(j-1);
  ipnum[6]=i;
  ipnum[7]=i+j*j*(j-1);
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
/*  zruseno 26.9.2004
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
*/

/**
   function computes stresses in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void linhex::res_mainip_stresses (long lcid,long eid)
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
void linhex::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta;
  vector gp,w,eps,epst,epstt,sig,auxsig;
  matrix d(tncomp,tncomp);

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


	    Mm->givestrain (lcid,ipp,eps);
	    

	    mxv (d,eps,auxsig);
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
   function computes stresses at nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002, JK
*/
void linhex::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ri,ci,ipnum);
  
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


void linhex::elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),r,gp,w,eps,epst,epstt,sig,auxsig,natcoord(3);
  matrix d(tncomp,tncomp);

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
	    
	    if (Mp->strainaver==0)
	      //appval (xi,eta,zeta,cncomp[jj],ncomp[jj],eps,stra);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,xi,eta,zeta,cncomp[jj],ncomp[jj],eps);
	    
	    mxv (d,eps,auxsig);
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
void linhex::res_allip_stresses (long lcid,long eid)
{
  res_allip_stresses (lcid,eid);
}

/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void linhex::allip_stresses (long lcid,long eid,long ri,long ci)
{
  res_mainip_stresses (lcid,eid);
}







void linhex::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra,**stre;
  vector coord,sig;
  
  /*
  if (Mp->stressaver==0){
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    //elem_strains (stra,lcid,eid,ri,ci);
    elem_stresses (stra,stre,lcid,eid,ri,ci);
  }

  */

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
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (ncp,sig);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
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
  

  //mainip_stresses (0,eid,0,0);
  
  
}




/**
   function computes other values in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void linhex::nod_others (long lcid,long eid,long ri,long ci)
{
  long i,j,k,l,ii,ipp,ncomp;
  double xi,eta,zeta,*lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),nzeta(nne),r(ndofe),gp,w,other,natcoord(3);
  ivector nodes(nne);

  lsm = new double [16];

  nodecoord (nxi,neta,nzeta);
  Mt->give_elemnodes (eid,nodes);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    ncomp = Mm->ip[ipp].ncompother;
    reallocv (ncomp,other);
    lhs = new double [ncomp*4];
    rhs = new double [ncomp*4];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,16);
    nullv (rhs,ncomp*4);
    
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
        eta=gp[j];
        for (k=0;k<intordsm[ii][ii];k++){
          zeta=gp[k];
          
          for (l=0; l<ncomp; l++)
            other[l] = Mm->ip[ipp].eqother[l];

          natcoord[0]=xi;  natcoord[1]=eta;  natcoord[2]=zeta;
          matassem_lsm (lsm,natcoord);
          rhsassem_lsm (rhs,natcoord,other);

          ipp++;
        }
      }
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,4,ncomp);
    Mt->other_nodal_values (nodes,nxi,neta,nzeta,lhs,3,0,ncomp,lcid);
    
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}



/**
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 24.10.2005
*/
void linhex::nod_eqother_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ncompo;
  ivector ipnum(nne),nod(nne);
  vector eqother;
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ri,ci,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    //Mm->givestrain (lcid,ipnum[i],eps);
    
    ncompo = Mm->givencompeqother (ipnum[i],0);
    reallocv (ncompo,eqother);
    Mm->giveeqother (ipnum[i],0,ncompo,eqother.a);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storeother (lcid,0,ncompo,eqother);
    
  }
}











/**
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   28.7.2001
*/
void linhex::gl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,l,ii,ipp,transf;
  double xi,eta,zeta,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w,gp,r(ndofe),eps(tncomp),sig,contr(ndofe),v(ndofe);
  matrix gm,tmat (ndofe,ndofe);

  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,r.a);
  
  
  //  transformation of nodal displacements
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    copyv (r,v);
    locglobtransf (r,v,tmat);
  }
  
  
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
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  
	  mxv (gm,r,eps);
	  
	  Mm->storestrain (lcid,ipp,eps);
	  
	  Mm->computenlstresses (ipp);
	  
	  Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
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
  
  //  transformation of nodal forces
  if (transf>0){
    transf_matrix (nodes,tmat);
    globloctransf (ifor,v,tmat);
    copyv (v,ifor);
  }
  
}

/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   
   JK, 24.9.2005
*/
void linhex::gnl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long i,j,k,l,ipp;
  double xi,eta,zeta,jac;
  vector w,gp,x(nne),y(nne),z(nne),sig(tncomp),contr(ndofe),r(ndofe);
  matrix gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  
  fillv (0.0,ifor);
  
  //  array for coordinates of integration points
  reallocv (intordsm[0][0],gp);
  //  array for weights of integration points
  reallocv (intordsm[0][0],w);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	//  computation of stress
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp);
	
	Mm->givestress (lcid,ipp,sig);
	
	//  strain-displacement (geometric) matrix
	gngeom_matrix (gm,r,x,y,z,xi,eta,zeta,jac);
	
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


void linhex::res_internal_forces (long lcid,long eid,vector &ifor)
{
  gl_internal_forces (lcid,eid,0,0,ifor);
  //gnl_internal_forces (lcid,eid,0,0,ifor);
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
  //elem_strains (stra,lcid,eid,ri,ci);

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
	  
	  //appval (xi,eta,zeta,0,tncomp,eps,stra);
	  
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
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
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

/**
   function computes nodal forces caused by presure on surface
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - vector of presure 
   @param eis - surface id 
   
   27.1.2006
*/
void linhex::node_forces_surf (long lcid,long eid,long *is,double *nv,vector &nf)
{
  long i,j;
  double xi,eta,zeta,jac;
  double *tnv;
  vector x(nne),y(nne),z(nne),gp,w,av(ndofe),v(ndofe);
  matrix n(napfun,ndofe),an(napfun,ndofe),am(ndofe,ndofe),tran(3,3);
  
  tnv = new double [12];

  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  reallocv (intordb,w);
  reallocv (intordb,gp);
  gauss_points (gp.a,w.a,intordb);
  
  //  surface number 1

  if (is[0]>0 ){
    xi=1.0;
    for (i=0;i<intordb;i++){
      eta=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,eta,zeta,0);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[0]==1){
      av[0] = nv[0*3+0];
      av[1] = nv[0*3+1];
      av[2] = nv[0*3+2];
      
      av[9] = nv[1*3+0];
      av[10] = nv[1*3+1];
      av[11] = nv[1*3+2];
      
      av[21] = nv[2*3+0];
      av[22] = nv[2*3+1];
      av[23] = nv[2*3+2];

      av[12] = nv[3*3+0];
      av[13] = nv[3*3+1];
      av[14] = nv[3*3+2];
    }


    if (is[0]==2){
      
      av[0] = nv[0*3+0];
      av[1] = nv[0*3+1];
      av[2] = nv[0*3+2];
      
      av[3] = nv[1*3+0];
      av[4] = nv[1*3+1];
      av[5] = nv[1*3+2];
      
      av[6] = nv[2*3+0];
      av[7] = nv[2*3+1];
      av[8] = nv[2*3+2];

      av[9]  = nv[3*3+0];
      av[10] = nv[3*3+1];
      av[11] = nv[3*3+2];
      
      locglob_nodeval (0,av,tnv,x,y,z);
      
      av[0] = tnv[0*3+0];
      av[1] = tnv[0*3+1];
      av[2] = tnv[0*3+2];
      
      av[9] = tnv[1*3+0];
      av[10] = tnv[1*3+1];
      av[11] = tnv[1*3+2];
      
      av[21] = tnv[2*3+0];
      av[22] = tnv[2*3+1];
      av[23] = tnv[2*3+2];

      av[12] = tnv[3*3+0];
      av[13] = tnv[3*3+1];
      av[14] = tnv[3*3+2];
    }


    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 2

  if (is[1]>0 ){
    eta=1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,zeta,1);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[1]==1){
      av[3] = nv[12+0*3+0];
      av[4] = nv[12+0*3+1];
      av[5] = nv[12+0*3+2];
      
      av[0] = nv[12+1*3+0];
      av[1] = nv[12+1*3+1];
      av[2] = nv[12+1*3+2];
      
      av[12] = nv[12+2*3+0];
      av[13] = nv[12+2*3+1];
      av[14] = nv[12+2*3+2];

      av[15] = nv[12+3*3+0];
      av[16] = nv[12+3*3+1];
      av[17] = nv[12+3*3+2];
    }


    if (is[1]==2){
      
      av[0] = nv[12+0*3+0];
      av[1] = nv[12+0*3+1];
      av[2] = nv[12+0*3+2];
      
      av[3] = nv[12+1*3+0];
      av[4] = nv[12+1*3+1];
      av[5] = nv[12+1*3+2];
      
      av[6] = nv[12+2*3+0];
      av[7] = nv[12+2*3+1];
      av[8] = nv[12+2*3+2];

      av[9]  = nv[12+3*3+0];
      av[10] = nv[12+3*3+1];
      av[11] = nv[12+3*3+2];
      
      locglob_nodeval (1,av,tnv,x,y,z);
      
      av[3] = tnv[0*3+0];
      av[4] = tnv[0*3+1];
      av[5] = tnv[0*3+2];
      
      av[0] = tnv[1*3+0];
      av[1] = tnv[1*3+1];
      av[2] = tnv[1*3+2];
      
      av[12] = tnv[2*3+0];
      av[13] = tnv[2*3+1];
      av[14] = tnv[2*3+2];

      av[15] = tnv[3*3+0];
      av[16] = tnv[3*3+1];
      av[17] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }

  //  surface number 3

  if (is[2]>0 ){
    xi=-1.0;
    for (i=0;i<intordb;i++){
      eta=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,eta,zeta,2);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[2]==1){
      av[6] = nv[24+0*3+0];
      av[7] = nv[24+0*3+1];
      av[8] = nv[24+0*3+2];
      
      av[3] = nv[24+1*3+0];
      av[4] = nv[24+1*3+1];
      av[5] = nv[24+1*3+2];
      
      av[15] = nv[24+2*3+0];
      av[16] = nv[24+2*3+1];
      av[17] = nv[24+2*3+2];

      av[18] = nv[24+3*3+0];
      av[19] = nv[24+3*3+1];
      av[20] = nv[24+3*3+2];
    }


    if (is[2]==2){
      
      av[0] = nv[24+0*3+0];
      av[1] = nv[24+0*3+1];
      av[2] = nv[24+0*3+2];
      
      av[3] = nv[24+1*3+0];
      av[4] = nv[24+1*3+1];
      av[5] = nv[24+1*3+2];
      
      av[6] = nv[24+2*3+0];
      av[7] = nv[24+2*3+1];
      av[8] = nv[24+2*3+2];

      av[9]  = nv[24+3*3+0];
      av[10] = nv[24+3*3+1];
      av[11] = nv[24+3*3+2];
      
      locglob_nodeval (2,av,tnv,x,y,z);
      
      av[6] = tnv[0*3+0];
      av[7] = tnv[0*3+1];
      av[8] = tnv[0*3+2];
      
      av[3] = tnv[1*3+0];
      av[4] = tnv[1*3+1];
      av[5] = tnv[1*3+2];
      
      av[15] = tnv[2*3+0];
      av[16] = tnv[2*3+1];
      av[17] = tnv[2*3+2];

      av[18] = tnv[3*3+0];
      av[19] = tnv[3*3+1];
      av[20] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 4

  if (is[3]>0 ){
    eta=-1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	zeta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,zeta,3);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[3]==1){
      av[9] = nv[36+0*3+0];
      av[10] = nv[36+0*3+1];
      av[11] = nv[36+0*3+2];
      
      av[6] = nv[36+1*3+0];
      av[7] = nv[36+1*3+1];
      av[8] = nv[36+1*3+2];
      
      av[18] = nv[36+2*3+0];
      av[19] = nv[36+2*3+1];
      av[20] = nv[36+2*3+2];

      av[21] = nv[36+3*3+0];
      av[22] = nv[36+3*3+1];
      av[23] = nv[36+3*3+2];
    }


    if (is[3]==2){
      
      av[0] = nv[36+0*3+0];
      av[1] = nv[36+0*3+1];
      av[2] = nv[36+0*3+2];
      
      av[3] = nv[36+1*3+0];
      av[4] = nv[36+1*3+1];
      av[5] = nv[36+1*3+2];
      
      av[6] = nv[36+2*3+0];
      av[7] = nv[36+2*3+1];
      av[8] = nv[36+2*3+2];

      av[9]  = nv[36+3*3+0];
      av[10] = nv[36+3*3+1];
      av[11] = nv[36+3*3+2];
      
      locglob_nodeval (3,av,tnv,x,y,z);
      
      av[9] = tnv[0*3+0];
      av[10]= tnv[0*3+1];
      av[11]= tnv[0*3+2];
      
      av[6] = tnv[1*3+0];
      av[7] = tnv[1*3+1];
      av[8] = tnv[1*3+2];
      
      av[18] = tnv[2*3+0];
      av[19] = tnv[2*3+1];
      av[20] = tnv[2*3+2];

      av[21] = tnv[3*3+0];
      av[22] = tnv[3*3+1];
      av[23] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }

  //  surface number 5

  if (is[4]>0 ){
    zeta=1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	eta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,eta,4);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[4]==1){
      av[0] = nv[48+0*3+0];
      av[1] = nv[48+0*3+1];
      av[2] = nv[48+0*3+2];
      
      av[3] = nv[48+1*3+0];
      av[4] = nv[48+1*3+1];
      av[5] = nv[48+1*3+2];
      
      av[6] = nv[48+2*3+0];
      av[7] = nv[48+2*3+1];
      av[8] = nv[48+2*3+2];

      av[9] = nv[48+3*3+0];
      av[10] = nv[48+3*3+1];
      av[11] = nv[48+3*3+2];
    }


    if (is[4]==2){
      
      av[0] = nv[48+0*3+0];
      av[1] = nv[48+0*3+1];
      av[2] = nv[48+0*3+2];
      
      av[3] = nv[48+1*3+0];
      av[4] = nv[48+1*3+1];
      av[5] = nv[48+1*3+2];
      
      av[6] = nv[48+2*3+0];
      av[7] = nv[48+2*3+1];
      av[8] = nv[48+2*3+2];

      av[9]  = nv[48+3*3+0];
      av[10] = nv[48+3*3+1];
      av[11] = nv[48+3*3+2];
      
      locglob_nodeval (4,av,tnv,x,y,z);
      
      av[0] = tnv[0*3+0];
      av[1] = tnv[0*3+1];
      av[2] = tnv[0*3+2];
      
      av[3] = tnv[1*3+0];
      av[4] = tnv[1*3+1];
      av[5] = tnv[1*3+2];
      
      av[6] = tnv[2*3+0];
      av[7] = tnv[2*3+1];
      av[8] = tnv[2*3+2];

      av[9] = tnv[3*3+0];
      av[10] = tnv[3*3+1];
      av[11] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }


  //  surface number 6

  if (is[5]>0 ){
    zeta=-1.0;
    for (i=0;i<intordb;i++){
      xi=gp[i];
      for (j=0;j<intordb;j++){
	eta=gp[j];
	
	jac2d_3d (jac,x,y,z,xi,eta,5);
	bf_matrix (n,xi,eta,zeta);
	jac = jac*w[i]*w[j];
	nnj (am.a,n.a,jac,n.m,n.n);
      }
    }
    
    if (is[5]==1){
      av[12] = nv[60+0*3+0];
      av[13] = nv[60+0*3+1];
      av[14] = nv[60+0*3+2];
      
      av[21] = nv[60+1*3+0];
      av[22] = nv[60+1*3+1];
      av[23] = nv[60+1*3+2];
      
      av[18] = nv[60+2*3+0];
      av[19] = nv[60+2*3+1];
      av[20] = nv[60+2*3+2];

      av[15] = nv[60+3*3+0];
      av[16] = nv[60+3*3+1];
      av[17] = nv[60+3*3+2];
    }


    if (is[5]==2){
      
      av[0] = nv[60+0*3+0];
      av[1] = nv[60+0*3+1];
      av[2] = nv[60+0*3+2];
      
      av[3] = nv[60+1*3+0];
      av[4] = nv[60+1*3+1];
      av[5] = nv[60+1*3+2];
      
      av[6] = nv[60+2*3+0];
      av[7] = nv[60+2*3+1];
      av[8] = nv[60+2*3+2];

      av[9]  = nv[60+3*3+0];
      av[10] = nv[60+3*3+1];
      av[11] = nv[60+3*3+2];
      
      locglob_nodeval (5,av,tnv,x,y,z);
      
      av[12] = tnv[0*3+0];
      av[13] = tnv[0*3+1];
      av[14] = tnv[0*3+2];
      
      av[21] = tnv[1*3+0];
      av[22] = tnv[1*3+1];
      av[23] = tnv[1*3+2];
      
      av[18] = tnv[2*3+0];
      av[19] = tnv[2*3+1];
      av[20] = tnv[2*3+2];

      av[15] = tnv[3*3+0];
      av[16] = tnv[3*3+1];
      av[17] = tnv[3*3+2];

    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }






  delete [] tnv;
}

void linhex::locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z)
{
  double norm;
  vector g1(3),g2(3),g3(3),lv(3),gv(3);
  matrix t(3,3);
  
  if (is==0){
    g1[0]=x[4]-x[7];
    g1[1]=y[4]-y[7];
    g1[2]=z[4]-z[7];
    
    g2[0]=x[3]-x[7];
    g2[1]=y[3]-y[7];
    g2[2]=z[3]-z[7];
  }
  if (is==1){
    g1[0]=x[5]-x[4];
    g1[1]=y[5]-y[4];
    g1[2]=z[5]-z[4];
    
    g2[0]=x[0]-x[4];
    g2[1]=y[0]-y[4];
    g2[2]=z[0]-z[4];
  }
  if (is==2){
    g1[0]=x[6]-x[5];
    g1[1]=y[6]-y[5];
    g1[2]=z[6]-z[5];
    
    g2[0]=x[1]-x[5];
    g2[1]=y[1]-y[5];
    g2[2]=z[1]-z[5];
  }
  if (is==3){
    g1[0]=x[7]-x[6];
    g1[1]=y[7]-y[6];
    g1[2]=z[7]-z[6];
    
    g2[0]=x[2]-x[6];
    g2[1]=y[2]-y[6];
    g2[2]=z[2]-z[6];
  }
  if (is==4){
    g1[0]=x[3]-x[2];
    g1[1]=y[3]-y[2];
    g1[2]=z[3]-z[2];
    
    g2[0]=x[1]-x[2];
    g2[1]=y[1]-y[2];
    g2[2]=z[1]-z[2];
  }
  if (is==5){
    g1[0]=x[5]-x[6];
    g1[1]=y[5]-y[6];
    g1[2]=z[5]-z[6];
    
    g2[0]=x[7]-x[6];
    g2[1]=y[7]-y[6];
    g2[2]=z[7]-z[6];
  }
  
  
  scprd (g1,g1,norm);
  norm=1.0/sqrt(norm);
  cmulv (norm,g1,g1);
  
  scprd (g1,g2,norm);
  g2[0]=g2[0]-norm*g1[0];
  g2[1]=g2[1]-norm*g1[1];
  g2[2]=g2[2]-norm*g1[2];
  
  scprd (g2,g2,norm);
  norm=1.0/sqrt(norm);
  cmulv (norm,g2,g2);
  
  g3[0]=g1[1]*g2[2]-g1[2]*g2[1];
  g3[1]=g1[2]*g2[0]-g1[0]*g2[2];
  g3[2]=g1[0]*g2[1]-g1[1]*g2[0];
  
  t[0][0]=g1[0];
  t[1][0]=g1[1];
  t[2][0]=g1[2];

  t[0][1]=g2[0];
  t[1][1]=g2[1];
  t[2][1]=g2[2];

  t[0][2]=g3[0];
  t[1][2]=g3[1];
  t[2][2]=g3[2];
  
  mxv (t.a,nv.a,tnv,3,3);
  mxv (t.a,nv.a+3,tnv+3,3,3);
  mxv (t.a,nv.a+6,tnv+6,3,3);
  mxv (t.a,nv.a+9,tnv+9,3,3);
  
}

/**
   function computes nodal forces caused by presure on surface
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - vector of presure 
   @param eis - surface id 
   
   4.2002, PF
*/
void linhex::node_forces_surf_old (long lcid,long eid,long *is,double *nv,vector &nf)
{
  long i,j,i1,ii;
  double xi=0.0,eta=0.0,zeta=0.0,jac, w1,w2;
  ivector nodes(nne);
  vector gx(nne),gy(nne),gz(nne),x(nnsurf),y(nnsurf),z(nnsurf),gp(intordb),w(intordb),av(ndofe),v(ndofe),pom(3);
  matrix n(napfun,ndofe),an(napfun,ndofe),am(ndofe,ndofe), tran(3,3);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (gx,gy,gz,eid);
  gauss_points (gp.a,w.a,intordb);
  fillv (0.0,nf);
  fillm (0.0,an);
  for (ii=0;ii<6;ii++){
    fillv (0.0,av);
    // is=0 not loading
    if (is[ii] ==0){}
    else {
      tran_mat( tran, gx, gy, gz, ii);
      // is=1 surface node 1,4,8,5
      if (ii ==0){
	xi=1.0;
        x[0]=gx[0];x[1]=gx[3];x[2]=gx[7];x[3]=gx[4]; 
        y[0]=gy[0];y[1]=gy[3];y[2]=gy[7];y[3]=gy[4]; 
        z[0]=gz[0];z[1]=gz[3];z[2]=gz[7];z[3]=gz[4]; 
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;     
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[i];  av[9+i]=nv[3+i];  av[21+i]=nv[6+i];  av[12+i]=nv[9+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=2 surface node 2,1,5,6
      else if (ii ==1){
	eta=1.0;
        x[0]=gx[1];x[1]=gx[0];x[2]=gx[4];x[3]=gx[5];
        y[0]=gy[1];y[1]=gy[0];y[2]=gy[4];y[3]=gy[5]; 
        z[0]=gz[1];z[1]=gz[0];z[2]=gz[4];z[3]=gz[5]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[3+i]=nv[12+i];  av[i]=nv[15+i];  av[12+i]=nv[18+i];  av[15+i]=nv[21+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=3 surface node 3,2,6,7
      else if (ii ==2){
	xi=-1.0;
        x[0]=gx[2];x[1]=gx[1];x[2]=gx[5];x[3]=gx[6]; 
        y[0]=gy[2];y[1]=gy[1];y[2]=gy[5];y[3]=gy[6]; 
        z[0]=gz[2];z[1]=gz[1];z[2]=gz[5];z[3]=gz[6]; 
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[6+i]=nv[24+i];  av[3+i]=nv[27+i];  av[15+i]=nv[30+i];  av[18+i]=nv[33+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=4 surface node 4,3,7,8
      else if (ii ==3){
	eta=-1.0;
        x[0]=gx[3];x[1]=gx[2];x[2]=gx[6];x[3]=gx[7]; 
        y[0]=gy[3];y[1]=gy[2];y[2]=gy[6];y[3]=gy[7];
        z[0]=gz[3];z[1]=gz[2];z[2]=gz[6];z[3]=gz[7]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;
	    
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[9+i]=nv[36+i];  av[6+i]=nv[39+i];  av[18+i]=nv[42+i];  av[21+i]=nv[45+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=5  surface node 1,2,3,4
      else if (ii ==4) {
	zeta=1.0;
        x[0]=gx[0];x[1]=gx[1];x[2]=gx[2];x[3]=gx[3]; 
        y[0]=gy[0];y[1]=gy[1];y[2]=gy[2];y[3]=gy[3]; 
        z[0]=gz[0];z[1]=gz[1];z[2]=gz[2];z[3]=gz[3]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    //	for constant			for (i1=0;i1<ndofe;i1++){
	    //								for (j1=0;j1<3;j1++){
	    //									an[j1][i1]=an[j1][i1]+n[j1][i1]*jac;}	}
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[48+i];  av[3+i]=nv[51+i];  av[6+i]=nv[54+i];  av[9+i]=nv[57+i];
	}
	//		fprintf (Out,"\n\n nn");
	//		for (i1=0;i1<ndofe;i1++){
	//			fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i1,an[0][i1],an[1][i1],an[2][i1]);}       
	
	// for constant		for (i=0;i<ndofe;i++){
	//						for (j=0;j<pom.n;j++){
	//							v[i]=an[j][i]*pom[j];}	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=6  surface node 6,5,8,7
      else if (ii ==5) {
	zeta=-1.0;
        x[0]=gx[5];x[1]=gx[4];x[2]=gx[7];x[3]=gx[6]; 
        y[0]=gy[5];y[1]=gy[4];y[2]=gy[7];y[3]=gy[6]; 
        z[0]=gz[5];z[1]=gz[4];z[2]=gz[7];z[3]=gz[6]; 
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[12+i]=nv[60+i];  av[15+i]=nv[63+i];  av[18+i]=nv[66+i];  av[21+i]=nv[69+i];
	}
	// Load in LCS and transformation to GCS
	if (is[ii] ==2){
	  lgvectortransf3dblock (av, tran);
	}
        mxv (am,av,v);  
      }
      
			fprintf (Out,"\n\n zatiz");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,av[i1],av[i1+1],av[i1+2]);}       
			fprintf (Out,"\n\n vloc");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,v[i1],v[i1+1],v[i1+2]);}       
//      lgvectortransf3dblock(v, tran);
			fprintf (Out,"\n\n vglob");
			for (i=0;i<nne;i++){
                i1=3*i;
				fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i,v[i1],v[i1+1],v[i1+2]);}       
      addv (nf,v,nf);      
    }
  }

  //  transformation of stiffness matrix to nodesystem
  //  transf = Mt->locsystems (nodes);
  //  if (transf>0){
  //    matrix tmat (ndofe,ndofe);
  //    transf_matrix (nodes,tmat);
  //    lgvectortransf3dblock (v,tmat);
  //  }

  //  fprintf (Out,"\n\n zatizeni na prvku cislo %ld",eid);
  //  for (i=0;i<ndofe;i++){
  //	  fprintf (Out,"\n %4ld   %20.10e",i,nf[i]);}
}

/**
   function computes transformation matrix on surface
   
   @param x, y - local coordinate  xL=adge12
   @param tran - tranformation metrix to GCS
   @param gx,gy,gz - vector of node global coordinate 
   @param is - surface id 
   
   4.2002, PF
*/
void linhex::tran_mat(matrix &tran, vector &gx, vector &gy, vector &gz, long is)
{
  long i,i1;
  double dl;
  matrix a(3,3);
// is=5  surface node 1,2,3,4
  if (is ==4) {
    for (i=0; i<3; i++) {
      i1=i+1; if(i1>2)i1=i1-3;
      a[0][i]=gx[i1]-gx[i];
      a[1][i]=gy[i1]-gy[i];
      a[2][i]=gz[i1]-gz[i];
    }
  }
// is=6  surface node 6,5,8,7
  else if (is ==5) {
    for (i=4; i<7; i++) {
      i1=i+1; if(i1>6)i1=i1-3;
      a[0][i-4]=gx[i1]-gx[i];
      a[1][i-4]=gy[i1]-gy[i];
      a[2][i-4]=gz[i1]-gz[i];
    }
  }
// is=1 surface node 1,4,8,5
// is=2 surface node 2,1,5,6
// is=3 surface node 3,2,6,7
// is=4 surface node 4,3,7,8
  else {
      if(is>0)i1=is-4;else i1=is;
      a[0][0]=gx[i1+3]-gx[is];
      a[1][0]=gy[i1+3]-gy[is];
      a[2][0]=gz[i1+3]-gz[is];
      a[0][1]=gx[i1+7]-gx[i1+3];
      a[1][1]=gy[i1+7]-gy[i1+3];
      a[2][1]=gz[i1+7]-gz[i1+3];
      a[0][2]=gx[is]-gx[i1+7];
      a[1][2]=gy[is]-gy[i1+7];
      a[2][2]=gz[is]-gz[i1+7];
  }
  
  dl=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);

  tran[0][0]=a[0][0]/dl;
  tran[1][0]=a[1][0]/dl;
  tran[2][0]=a[2][0]/dl;
  
  tran[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
  tran[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
  tran[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

  dl=sqrt(tran[0][2]*tran[0][2]+tran[1][2]*tran[1][2]+tran[2][2]*tran[2][2]);

  tran[0][2]=tran[0][2]/dl;
  tran[1][2]=tran[1][2]/dl;
  tran[2][2]=tran[2][2]/dl;
  
  tran[0][1]=tran[1][2]*tran[2][0]-tran[2][2]*tran[1][0];
  tran[1][1]=tran[2][2]*tran[0][0]-tran[0][2]*tran[2][0];
  tran[2][1]=tran[0][2]*tran[1][0]-tran[1][2]*tran[0][0];

  dl=sqrt(tran[0][1]*tran[0][1]+tran[1][1]*tran[1][1]+tran[2][1]*tran[2][1]);

  tran[0][1]=tran[0][1]/dl;
  tran[1][1]=tran[1][1]/dl;
  tran[2][1]=tran[2][1]/dl;
/*  
// local coordinate x=s12
    sl1[0]=0.0;    sl1[1]=0.0;    sl1[2]=0.0;
    sl2[0]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
    sl2[1]=0.0;    sl2[2]=0.0;
    sl3[0]=(s3(0)-s1(0))*tran(0,0)+(s3(1)-s1(1))*tran(1,0)+(s3(2)-s1(2))*tran(2,0);
    sl3[1]=(s3(0)-s1(0))*tran(0,1)+(s3(1)-s1(1))*tran(1,1)+(s3(2)-s1(2))*tran(2,1);
    sl3[2]=0.0;
*/
// local coordinate x=s12
//    x[0]=0.0;
//    y[0]=0.0;
//    x[1]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
//    y[1]=0.0;
//    x[2]=(gx(2)-gx(0))*tran(0,0)+(gy(2)-gy(0))*tran(1,0)+(gz(2)-gz(0))*tran(2,0);
//    y[2]=(gx(2)-gx(0))*tran(0,1)+(gy(2)-gy(0))*tran(1,1)+(gz(2)-gz(0))*tran(2,1);
//    x[3]=(gx(3)-gx(0))*tran(0,0)+(gy(3)-gy(0))*tran(1,0)+(gz(3)-gz(0))*tran(2,0);
//    y[3]=(gx(3)-gx(0))*tran(0,1)+(gy(3)-gy(0))*tran(1,1)+(gz(3)-gz(0))*tran(2,1);
}

/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 17.8.2004
*/
void linhex::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  eigstrain_forces (lcid,eid,nfor);
}

/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 17.8.2004
*/
void linhex::eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),gp,w,eigstr,sig,contr(ndofe);
  matrix d(tncomp,tncomp),gm;
  
  Mt->give_node_coord3d (x,y,z,eid);
  fillv (0.0,nfor);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eigstr);
    reallocv (ncomp[ii],sig);
    reallocm (ncomp[ii],ndofe,gm);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ii][ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	
	  Mm->giveeigstrain (ipp,cncomp[ii],ncomp[ii],eigstr);
	  
	  //  matrix of stiffness of the material
	  Mm->matstiff (d,ipp);
	  
	  mxv (d,eigstr,sig);
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  
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
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   17.8.2004, JK
*/
void linhex::intpointval (long eid,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k,l;
  double xi,eta,zeta;
  vector w,gp;
  
  l=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    ipval[l]=approx (xi,eta,zeta,nodval);
	    l++;
	  }
	}
      }
    }
  }
}



void linhex::aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),gp,w,eps;

  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  //Mm->givestrain (lcid,ipp,eps);
	  Mm->givestress (lcid,ipp,eps);
	  
	  
	  volume+=w[i]*w[j]*w[k]*jac;
	  
	  for (l=0;l<averstra.n;l++){
	    averstra[l]+=eps[l]*w[i]*w[j]*w[k]*jac;
	  }
	  
	  ipp++;
	}
      }
    }
  }
  
}
