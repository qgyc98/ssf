#include "plelemrotlq.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>



planeelemrotlq::planeelemrotlq (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=2;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  number of surface
  nsurf=1;
  //  number of nodes on surface
  nnsurf=4;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;

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
  
  intordsm[0][0]=3;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=3;
  
  nip[0][0]=9;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=9;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

}

planeelemrotlq::~planeelemrotlq (void)
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
   function prepares auxiliary data for element with rotational degrees of freedom
   
   @param x,y - vectors of node coordinates
   @param l - %vector of lengths
   @param nx,ny - vectors of normal vectors components

   8.12.2001
*/
void planeelemrotlq::auxdata (vector &x,vector &y,vector &l,vector &nx,vector &ny)
{
  double sx,sy;
  sx=x[1]-x[0];  sy=y[1]-y[0];
  l[0]=sqrt(sx*sx+sy*sy);
  nx[0]=sy/l[0];  ny[0]=sx*(-1.0)/l[0];
  
  sx=x[2]-x[1];  sy=y[2]-y[1];
  l[1]=sqrt(sx*sx+sy*sy);
  nx[1]=sy/l[1];  ny[1]=sx*(-1.0)/l[1];

  sx=x[3]-x[2];  sy=y[3]-y[2];
  l[2]=sqrt(sx*sx+sy*sy);
  nx[2]=sy/l[2];  ny[2]=sx*(-1.0)/l[2];

  sx=x[0]-x[3];  sy=y[0]-y[3];
  l[3]=sqrt(sx*sx+sy*sy);
  nx[3]=sy/l[3];  ny[3]=sx*(-1.0)/l[3];
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
   JK
*/
double planeelemrotlq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_4_2d (bf.a,xi,eta);
  scprd (bf,nodval,f);
  return f;
}

/**
   function assembles %matrix of base (approximation, shape) function
   
   @param n - %matrix of base functions
   @param xi,eta - natural coordinates
   @param l - %vector of edge lengths
   @param nx,ny - coordinates of normal vectors
   
   JK, 9.7.2001
*/
void planeelemrotlq::bf_matrix (matrix &n,double xi,double eta,vector &l,vector &nx,vector &ny)
{
  long i,j,k,m,ii,jj;
  vector bf(12);
  
  nullm (n);

  bf_rot_4_2d (bf.a,xi,eta,nx.a,ny.a,l.a);
  
  j=0;  k=1;  m=2;  ii=4;  jj=8;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];
    n[0][m]=bf[ii];
    n[1][k]=bf[i];
    n[1][m]=bf[jj];
    j+=3;  k+=3;  m+=3;  ii++;  jj++;
  }
}

/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param l - %vector of edge lengths
   @param nx,ny - coordinates of normal vectors
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemrotlq::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,
				  vector &l,vector &nx,vector &ny,double &jac)
{
  long i,i1,i2,i3,ii,jj;
  vector dx(3*nne),dy(3*nne);
  
  dx_bf_rot_4_2d (dx.a,xi,eta,nx.a,ny.a,l.a);
  dy_bf_rot_4_2d (dy.a,xi,eta,nx.a,ny.a,l.a);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullm (gm);
  
  i1=0;  i2=1;  i3=2;  ii=4;  jj=8;

  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[0][i3]=dx[ii];
    gm[1][i2]=dy[i];
    gm[1][i3]=dy[jj];
    gm[2][i1]=dy[i];
    gm[2][i2]=dx[i];
    gm[2][i3]=dy[ii]+dx[jj];
    i1+=3;  i2+=3;  i3+=3;  ii++;  jj++;
  }

  /*
  gm[0][0]=dx[0];
  gm[0][2]=dx[4];
  gm[0][3]=dx[1];
  gm[0][5]=dx[5];
  gm[0][6]=dx[2];
  gm[0][8]=dx[6];
  gm[0][9]=dx[3];
  gm[0][11]=dx[7];

  gm[1][1]=dy[0];
  gm[1][2]=dy[8];
  gm[1][4]=dy[1];
  gm[1][5]=dy[9];
  gm[1][7]=dy[2];
  gm[1][8]=dy[10];
  gm[1][10]=dy[3];
  gm[1][11]=dy[11];
  
  gm[2][0]=dy[0];
  gm[2][1]=dx[0];
  gm[2][2]=dy[4]+dx[8];
  gm[2][3]=dy[1];
  gm[2][4]=dx[1];
  gm[2][5]=dy[5]+dx[9];
  gm[2][6]=dy[2];
  gm[2][7]=dx[2];
  gm[2][8]=dy[6]+dx[10];
  gm[2][9]=dy[3];
  gm[2][10]=dx[3];
  gm[2][11]=dy[7]+dx[11];
  */
}

/**
   function assembles block of strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param ri - row index
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param l - %vector of edge lengths
   @param nx,ny - coordinates of normal vectors
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemrotlq::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,
					vector &l,vector &nx,vector &ny,double &jac)
{
  long i,i1,i2,i3,ii,jj;
  vector dx(3*nne),dy(3*nne);
  
  dx_bf_rot_4_2d (dx.a,xi,eta,nx.a,ny.a,l.a);
  dy_bf_rot_4_2d (dy.a,xi,eta,nx.a,ny.a,l.a);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullm (gm);
  
  if (ri==0){
    i1=0;  i2=1;  i3=2;  ii=4;  jj=8;

    for (i=0;i<nne;i++){
      gm[0][i1]=dx[i];
      gm[0][i3]=dx[ii];
      gm[1][i2]=dy[i];
      gm[1][i3]=dy[jj];
      i1+=3;  i2+=3;  i3+=3;  ii++;  jj++;
    }

    /*
    gm[0][0]=dx[0];
    gm[0][2]=dx[4];
    gm[0][3]=dx[1];
    gm[0][5]=dx[5];
    gm[0][6]=dx[2];
    gm[0][8]=dx[6];
    gm[0][9]=dx[3];
    gm[0][11]=dx[7];

    gm[1][1]=dy[0];
    gm[1][2]=dy[8];
    gm[1][4]=dy[1];
    gm[1][5]=dy[9];
    gm[1][7]=dy[2];
    gm[1][8]=dy[10];
    gm[1][10]=dy[3];
    gm[1][11]=dy[11];
    */
  }
  if (ri==1){
    i1=0;  i2=1;  i3=2;  ii=4;  jj=8;

    for (i=0;i<nne;i++){
      gm[0][i1]=dy[i];
      gm[0][i2]=dx[i];
      gm[0][i3]=dy[ii]+dx[jj];
      i1+=3;  i2+=3;  i3+=3;  ii++;  jj++;
    }

    /*
    gm[0][0]=dy[0];
    gm[0][1]=dx[0];
    gm[0][2]=dy[4]+dx[8];

    gm[0][3]=dy[1];
    gm[0][4]=dx[1];
    gm[0][5]=dy[5]+dx[9];

    gm[0][6]=dy[2];
    gm[0][7]=dx[2];
    gm[0][8]=dy[6]+dx[10];

    gm[0][9]=dy[3];
    gm[0][10]=dx[3];
    gm[0][11]=dy[7]+dx[11];
    */
  }
  
}

/**
   function assembles geometric %matrix used for rotation
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param l - %vector of edge lengths
   @param nx,ny - coordinates of normal vectors
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemrotlq::addgeommat (matrix &gm,vector &x,vector &y,double xi,double eta,
				 vector &l,vector &nx,vector &ny,double &jac)
{
  long i,i1,i2,i3,ii,jj;
  vector bf(ndofe),dx(ndofe),dy(ndofe);
  
  bf_rot_4_2d (bf.a,xi,eta,nx.a,ny.a,l.a);
  dx_bf_rot_4_2d (dx.a,xi,eta,nx.a,ny.a,l.a);
  dy_bf_rot_4_2d (dy.a,xi,eta,nx.a,ny.a,l.a);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullm (gm);
  
  i1=0;  i2=1;  i3=2;  ii=4;  jj=8;

  for (i=0;i<nne;i++){
    gm[0][i1]=0.0-dy[i]/2.0;
    gm[0][i2]=dx[i]/2.0;
    gm[0][i3]=dx[jj]/2.0-dy[ii]/2.0-bf[i];
    i1+=3;  i2+=3;  i3+=3;  ii++;  jj++;
  }

  /*
  gm[0][0]=0.0-dy[0]/2.0;
  gm[0][1]=dx[0]/2.0;
  gm[0][2]=dx[8]/2.0-dy[4]/2.0-bf[0];

  gm[0][3]=0.0-dy[1]/2.0;
  gm[0][4]=dx[1]/2.0;
  gm[0][5]=dx[9]/2.0-dy[5]/2.0-bf[1];

  gm[0][6]=0.0-dy[2]/2.0;
  gm[0][7]=dx[2]/2.0;
  gm[0][8]=dx[10]/2.0-dy[6]/2.0-bf[2];

  gm[0][9]=0.0-dy[3]/2.0;
  gm[0][10]=dx[3]/2.0;
  gm[0][11]=dx[11]/2.0-dy[7]/2.0-bf[3];
  */
}


/**
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
   
   JK
*/
void planeelemrotlq::dmatblock (long ri,long ci,matrix &d, matrix &dd)
{
  nullm (dd);
  
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
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   JK, 9.7.2001
*/
void planeelemrotlq::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,n,m;
  
  nullm (tmat);

  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*3][i*3]   = Mt->nodes[nodes[i]].e1[0];    tmat[i*3][i*3+1]   = Mt->nodes[nodes[i]].e2[0];  tmat[i*3][i*3+2]   = 0.0;
      tmat[i*3+1][i*3] = Mt->nodes[nodes[i]].e1[1];    tmat[i*3+1][i*3+1] = Mt->nodes[nodes[i]].e2[1];  tmat[i*3+1][i*3+2] = 0.0;
      tmat[i*3+2][i*3] = 0.0;                          tmat[i*3+2][i*3+1] = 0.0;                        tmat[i*3+2][i*3+2] = 1.0;
    }
  }
}

/**
   function computes stiffness %matrix of plane stress quadrilateral
   finite element with rotational degrees of freedom with
   bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of nodal coordinates

   8.12.2001
*/
void planeelemrotlq::stiffness_matrix (long eid,long /*ri*/,long /*ci*/,matrix &sm,vector &x,vector &y)
{
  long i,j,ii,jj,ipp;
  double xi,eta,ww1,ww2,jac,thick;
  ivector nodes(nne);
  vector l(nne),nx(nne),ny(nne),w,gp,t(nne);
  matrix gmr,gmc,dd,d(tncomp,tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  auxdata (x,y,l,nx,ny);
  
  nullm (sm);
  
  for (ii=0;ii<nb;ii++){
    
    reallocm (ncomp[ii],ndofe,gmr);
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      ipp=Mt->elements[eid].ipp[ii][jj];
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];  ww1=w[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];  ww2=w[j];
	  
	  //  geometrical matrix
	  geom_matrix_block (gmr,ii,x,y,xi,eta,l,nx,ny,jac);
	  geom_matrix_block (gmc,jj,x,y,xi,eta,l,nx,ny,jac);
	  
	  //  matrix of stiffness of the material
	  Mm->matstiff (d,ipp);
	  dmatblock (ii,jj,d,dd);
	  
	  thick = approx (xi,eta,t);
	  
	  jac*=thick*ww1*ww2;
	  
	  //  contribution to the stiffness matrix of the element
	  //bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
	  bdbjac (sm,gmr,dd,gmc,jac);
	  
	  ipp++;
	}
      }
    }
  }
  
  long intord = 2;
  reallocv (intord,w);
  reallocv (intord,gp);
  gauss_points (gp.a,w.a,intord);
  reallocm (intord,ndofe,gmr);
  reallocm (ndofe,ndofe,gmc);
  for (i=0;i<intord;i++){
    for (j=0;j<intord;j++){
      xi=gp[i];  eta=gp[j];
      
      addgeommat (gmr,x,y,xi,eta,l,nx,ny,jac);
      mtxm (gmr,gmr,gmc);
      thick = approx (xi,eta,t);
      cmulm (thick*jac*w[i]*w[j],gmc);
      addm (sm,gmc,sm);
    }
  }

}

/**
   function computes stiffness %matrix of quadrilateral element
   with rotational degrees of freedom
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK
*/
void planeelemrotlq::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);

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
   function computes mass %matrix of the plane stress quadrilateral
   finite element with rotationale degrees of freedom wtih
   bilinear approximation functions

   @param eid - number of element
   @param mm - mass %matrix
   @param x,y - vectors of nodal coordinates

   8.12.2001
*/
void planeelemrotlq::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick,rho;
  ivector nodes(nne);
  vector l(nne),nx(nne),ny(nne),w(intordmm),gp(intordmm),t(nne),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);
  auxdata (x,y,l,nx,ny);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (mm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta,l,nx,ny);
      
      thick = approx (xi,eta,t);
      rho = approx (xi,eta,dens);
      jac*=w1*w2*thick*rho;
      
      nnj (mm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function assembles mass %matrix of plane stress rectangular
   finite element with rotational degrees of freedom
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK
*/
void planeelemrotlq::res_mass_matrix (long eid,matrix &mm)
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
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}


/**
   function computes load %matrix of the plane stress rectangular
   finite element with rotational degrees of freedom
   load %vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   @param x,y - node coordinates
   
   25.7.2001
*/
void planeelemrotlq::load_matrix (long eid,matrix &lm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector l(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  //Mt->give_node_coord2d (x,y,eid);
  auxdata (x,y,l,nx,ny);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta,l,nx,ny);
      
      thick = approx (xi,eta,t);
      jac*=w1*w2*thick;
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function computes load %matrix of the plane stress rectangular
   finite element with rotational degrees of freedom
   load %vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 28. 8. 2018
*/
void planeelemrotlq::res_load_matrix (long eid,matrix &lm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);

  load_matrix (eid,lm,x,y);

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


void planeelemrotlq::res_ip_strains (long lcid,long eid)
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
   function computes strains in integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param ii - number of block
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   10.5.2002
*/
void planeelemrotlq::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ii,ipp;
  double xi,eta,jac;
  vector gp,w,eps,nx(ned),ny(ned),l(ned);
  matrix gm;
  
  auxdata (x,y,l,nx,ny);
  
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
	
	geom_matrix_block (gm,ii,x,y,xi,eta,l,nx,ny,jac);
	mxv (gm,r,eps);
	
	Mm->storestrain (lcid,ipp,cncomp[ii],eps);
	ipp++;
      }
    }
    
  }
}


/**
   function computes strains in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void planeelemrotlq::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(gncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq (ipp,intordsm[0][0],ipnum);
  
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
   
   10.5.2002
*/
void planeelemrotlq::nod_strains_comp (long lcid,long eid,double **stra)
{
  long i,j;
  double jac;
  ivector nodes(nne);
  vector x(nne),y(nne),l(nne),nx(nne),ny(nne),nxi(nne),neta(nne),r(ndofe),eps(tncomp),natcoord(2),aux;
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
  
  //  natural coordinates of element nodes
  nodecoord (nxi,neta);

  auxdata (x,y,l,nx,ny);

  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,nxi[i],neta[i],l,nx,ny,jac);
    //  strain computation
    mxv (gm,r,eps);
    
    for (j=0;j<eps.n;j++){
      stra[i][j]=eps[j];
    }
  }

}


/**
   function computes strains at strain points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
/*
void planeelemrotlq::strains (long lcid,long eid,long ri,long ci)
{
}
*/

/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   
   10.5.2002
*/
void planeelemrotlq::nodecoord (vector &xi,vector &eta)
{
  xi[0] =  1.0;  eta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;
  xi[3] =  1.0;  eta[3] = -1.0;
}

void planeelemrotlq::res_ip_stresses (long lcid,long eid)
{
  long ri,ci;
  ri=0;
  ci=0;
  compute_nlstress (lcid,eid,ri,ci);
}


/**
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 10.5.2002
*/
void planeelemrotlq::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(gncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq (ipp,intordsm[0][0],ipnum);
  
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

void planeelemrotlq::stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{

}



/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void planeelemrotlq::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes internal forces  for nonlocal models

   this function is used in plane stress/strain elements (function is called
   by function res_nonloc_internal_forces) and shell elements
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void planeelemrotlq::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes increment of  internal forces (from correct stresses increment)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void planeelemrotlq::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=stressincr;
  
  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of nodal forces
   @param x,y - vectors of nodal coordinates

   JK, 28.7.2001
*/
void planeelemrotlq::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void planeelemrotlq::res_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes resulting internal forces for nonlocal models
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo 7.2008
*/
void planeelemrotlq::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void planeelemrotlq::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes resulting contributions from eigenstrains to the right hand side
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - %vector of internal forces

   TKo 7.2008
*/
void planeelemrotlq::res_eigstrain_forces(long lcid,long eid,vector &nfor)
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
    //globloctransf (ifor,v,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void planeelemrotlq::compute_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
        //  computation of correct stresses
	if (Mp->strcomp==1)
          Mm->computenlstresses (ipp,Mm->ip[ipp]);
        ipp++;
      }
    }
  }
}



/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void planeelemrotlq::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
        //  computation of correct stresses
	if (Mp->strcomp==1)
          Mm->computenlstressesincr (ipp);
        ipp++;
      }
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
void planeelemrotlq::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
        //  computation of correct stresses
	if (Mp->strcomp==1)
          Mm->computenlstresses (ipp,Mm->ip[ipp]);
        ipp++;
      }
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
void planeelemrotlq::compute_nonloc_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
        //  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->compnonloc_nlstresses (ipp);
        ipp++;
      }
    }
  }
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemrotlq::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
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
void planeelemrotlq::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv, vector &x, vector &y)
{
  long i,j,ii,ipp;
  double xi,eta,jac,thick;
  ivector nodes(nne),cn(ndofe);
  vector w, gp, t(nne),ipv,contr(ndofe),l(ned),nx(ned),ny(ned);
  matrix gm;

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  auxdata (x,y,l,nx,ny);
  nullv (nv);

  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],ipv);
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	thick = approx (xi,eta,t);
	
        //  function assembles required quantity at integration point
        Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);

        //  strain-displacement (geometric) matrix
	geom_matrix_block (gm,ii,x,y,xi,eta,l,nx,ny,jac);
        //  contribution to the nodal values
	mtxv (gm,ipv,contr);
	cmulv (jac*w[i]*w[j]*thick,contr);

        //  summation
        addv(contr,nv,nv);      
	ipp++;
      }
    }
  }
}

/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - id of loaded edge
   @param nv - nodal values
   @param nf - %vector of nodal forces
   
   11. 8. 2015, JK
*/
void planeelemrotlq::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double ww,jac,xi,eta;
  vector l(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),av(ASTCKVEC(ndofe)),v(ASTCKVEC(ndofe));
  matrix n(napfun,ndofe),am(ndofe,ndofe);
  
  Mt->give_node_coord2d (x,y,eid);
  auxdata (x,y,l,nx,ny);
  gauss_points (gp.a,w.a,intordb);

  if (le[0]==1){
    nullm (am);
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta,l,nx,ny);
      
      jac1d_2d (jac,x,y,xi,0);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[0]=nv[0];  av[1]=nv[1];  av[2]=nv[2];  av[3]=nv[3];  av[4]=nv[4];  av[5]=nv[5];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    nullm (am);
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];

      bf_matrix (n,xi,eta,l,nx,ny);
      
      jac1d_2d (jac,x,y,eta,1);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[3]=nv[6];  av[4]=nv[7];  av[5]=nv[8];  av[6]=nv[9];  av[7]=nv[10];  av[8]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    nullm (am);
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta,l,nx,ny);
      
      jac1d_2d (jac,x,y,xi,2);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[6]=nv[12];  av[7]=nv[13];  av[8]=nv[14];  av[9]=nv[15];  av[10]=nv[16];  av[11]=nv[17];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    nullm (am);
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta,l,nx,ny);
      
      jac1d_2d (jac,x,y,eta,3);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[9]=nv[18];  av[10]=nv[19];  av[11]=nv[20];  av[0]=nv[21];  av[1]=nv[22];  av[2]=nv[23];
    mxv (am,av,v);  addv (nf,v,nf);
  }
}

/**
   function computes nodal forces caused by surface load
   
   @param eid - element id
   @param nodvals - array of load nodal values
   @param x,y - vectors of nodal coordinates
   @param nf - %vector of nodal forces
   
   JK, 9. 3. 2020
*/
void planeelemrotlq::surfload (long /*eid*/,double *nodvals, vector &x, vector &y,vector &nf)
{
  long i,j;
  double w1,w2,jac,xi,eta;
  vector gp(ASTCKVEC(intordmm)),w(ASTCKVEC(intordmm));
  vector l(ASTCKVEC(nne)),nx(ASTCKVEC(nne)),ny(ASTCKVEC(nne));
  matrix n(napfun,ndofe),am(ndofe,ndofe);
  
  //  auxiliary values
  auxdata (x,y,l,nx,ny);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (am);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      
      jac_2d (jac,x,y,xi,eta);
      
      bf_matrix (n,xi,eta,l,nx,ny);

      jac*=w1*w2;
      
      nnjac (am,n,n,jac);
    }
  }
  
  mxv (am.a,nodvals,nf.a,ndofe,ndofe);
}

/**
   function computes nodal forces caused by surface load
   
   @param eid - element id
   @param nodvals - array of load nodal values
   @param x,y - vectors of nodal coordinates
   @param nf - %vector of nodal forces
   
   JK, 9. 3. 2020
*/
void planeelemrotlq::node_forces_surf (long eid,double *nodvals,vector &nf)
{
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  
  //  nodal forces
  surfload (eid,nodvals,x,y,nf);
}

/**
   The function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param coord - %vector with global coordinates of integration point (ouput)
   
   @return The function returns global coordinates in the argument coord.

   10.1.2002
*/
void planeelemrotlq::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,ii;
  double xi,eta;
  vector x(nne),y(nne),w(intordsm[ri][ci]),gp(intordsm[ri][ci]);
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ri][ci];j++){
      eta=gp[j];
      if (ii==ipp){
	coord[0]=approx (xi,eta,x);
	coord[1]=approx (xi,eta,y);
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
void planeelemrotlq::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, j, ii, ri, ci;
  double xi, eta;
  vector w, gp;
  
  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      if (intordsm[ri][ci] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp));
      gauss_points (gp.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];
  
      for (i=0;i<intordsm[ri][ci];i++){
        xi=gp[i];
        for (j=0;j<intordsm[ri][ci];j++){
          eta=gp[j];
          if (ii==ipp){
            ncoord[0]=xi;
            ncoord[1]=eta;
            ncoord[2]=0.0;
          }
          ii++;
        }
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
void planeelemrotlq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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



/**
   function computes volume appropriate to integration point
   
   2.3.2004, JK
   07.2008 TKo - multiplictaion by thickness added
*/
void planeelemrotlq::ipvolume (long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  double xi,eta,jac,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),t(nne),w,gp;
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);
  
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

          thick = approx (xi,eta,t);
	  jac_2d (jac,x,y,xi,eta);
	  jac*=w[i]*w[j]*thick;
	  Mm->storeipvol (ipp,jac);
	    
	  ipp++;
	}
      }
    }
  }
  reallocv (1,w);
  reallocv (1,gp);
  gauss_points (gp.a,w.a,1);
  for (i=0;i<1;i++){
    for (j=0;j<1;j++){
      xi=gp[i];  eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      thick = approx (xi,eta,t);
      jac*=w[i]*w[j]*thick;
    }
  }
}
