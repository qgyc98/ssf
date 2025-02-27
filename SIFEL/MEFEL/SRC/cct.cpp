#include "cct.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "gadaptivity.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>

cctelem::cctelem (void)
{
  long i,j;

  //  number nodes on element
  nne=3;
  //  number of DOFs on element
  ndofe=9;
  //  number of strain/stress components
  tncomp=5;
  //  number of functions approximated
  napfun=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
  //  order of numerical integration of mass matrix
  intordmm=3;
  //  strain/stress state
  ssst=plates;

  //  number of blocks (parts of geometric matrix)
  nb=2;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=3;
  ncomp[1]=2;

  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=3;
  
  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=1;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=1;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=1;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=1;
}

cctelem::~cctelem (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] intordsm[i];
    delete [] nip[i];
  }
  delete [] intordsm;
  delete [] nip;

  delete [] cncomp;
  delete [] ncomp;
}


/**
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

   28.3.2002
*/
double cctelem::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function approximates function defined by nodal values

   @param xi,eta - natural coordinates
   @param nodval - nodal values

*/
double cctelem::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector areacoord(3);
  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  scprd (areacoord,nodval,f);
  return f;
}


/**
   function assembles %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param x,y - array containing node coordinates
   @param areacoord - array of area coordinates
   
   JK, 19.7.2001
*/
void cctelem::bf_matrix (matrix &n,vector &x,vector &y,vector &areacoord)
{
  vector bf(9),sx(3),sy(3);
  
  sx[0]=x[1]-x[0];  sx[1]=x[2]-x[1];  sx[2]=x[0]-x[2];
  sy[0]=y[1]-y[0];  sy[1]=y[2]-y[1];  sy[2]=y[0]-y[2];
  
  bf_cct (bf.a,areacoord.a,sx.a,sy.a);
  
  fillm (0.0,n);
  
  n[0][0]=bf[0];
  n[0][1]=bf[1];
  n[0][2]=bf[2];
  n[0][3]=bf[3];
  n[0][4]=bf[4];
  n[0][5]=bf[5];
  n[0][6]=bf[6];
  n[0][7]=bf[7];
  n[0][8]=bf[8];
  
  n[1][1]=bf[0];
  n[1][4]=bf[3];
  n[1][7]=bf[6];

  n[2][2]=bf[0];
  n[2][5]=bf[3];
  n[2][8]=bf[6];

}

/**
   function assembles part of geometric %matrix
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param areacoord - area coordinates
   
   JK, 19.7.2001
*/
void cctelem::geom_matrix_block (matrix &gm,long bi,vector &x,vector &y,vector &areacoord)
{
  double det;
  vector sx(3),sy(3),b(3),c(3),bf(9),dx(9),dy(9);

  sx[0]=x[1]-x[0];  sx[1]=x[2]-x[1];  sx[2]=x[0]-x[2];
  sy[0]=y[1]-y[0];  sy[1]=y[2]-y[1];  sy[2]=y[0]-y[2];
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);  
  
  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  bf_cct (bf.a,areacoord.a,sx.a,sy.a);
  dx_cct (dx.a,areacoord.a,b.a,sx.a,sy.a);
  dy_cct (dy.a,areacoord.a,c.a,sx.a,sy.a);
  
  fillm (0.0,gm);
  
  if (bi==0){
    gm[0][2]=dx[0];
    gm[0][5]=dx[3];
    gm[0][8]=dx[6];
    
    gm[1][1]=0.0-dy[0];
    gm[1][4]=0.0-dy[3];
    gm[1][7]=0.0-dy[6];
    
    gm[2][1]=0.0-dx[0];
    gm[2][2]=dy[0];
    gm[2][4]=0.0-dx[3];
    gm[2][5]=dy[3];
    gm[2][7]=0.0-dx[6];
    gm[2][8]=dy[6];
  }
  if (bi==1){
    gm[0][0]=dy[0];
    gm[0][1]=dy[1]-bf[0];
    gm[0][2]=dy[2];
    gm[0][3]=dy[3];
    gm[0][4]=dy[4]-bf[3];
    gm[0][5]=dy[5];
    gm[0][6]=dy[6];
    gm[0][7]=dy[7]-bf[6];
    gm[0][8]=dy[8];
    
    gm[1][0]=dx[0];
    gm[1][1]=dx[1];
    gm[1][2]=dx[2]+bf[0];
    gm[1][3]=dx[3];
    gm[1][4]=dx[4];
    gm[1][5]=dx[5]+bf[3];
    gm[1][6]=dx[6];
    gm[1][7]=dx[7];
    gm[1][8]=dx[8]+bf[6];
  }
  
}

/**
   function assembles geometric %matrix of the cct element
   
   @param gm - geometric %matrix
   @param x,y - array of node coordinates
   @param areacoord - area coordinates
   
   JK, 19.7.2001
*/
void cctelem::geom_matrix (matrix &gm,vector &x,vector &y,vector &areacoord)
{
  double det;
  vector sx(3),sy(3),b(3),c(3),bf(9),dx(9),dy(9);
  
  sx[0]=x[1]-x[0];  sx[1]=x[2]-x[1];  sx[2]=x[0]-x[2];
  sy[0]=y[1]-y[0];  sy[1]=y[2]-y[1];  sy[2]=y[0]-y[2];

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);  

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  bf_cct (bf.a,areacoord.a,sx.a,sy.a);
  dx_cct (dx.a,areacoord.a,b.a,sx.a,sy.a);
  dy_cct (dy.a,areacoord.a,c.a,sx.a,sy.a);

  fillm (0.0,gm);

  // expected ordering of nodal unknowns: d^T = {w_1, phi_{x1}, phi_{y1}, w_2, phi_{x2}, phi_{y2}, w_3, phi_{x3}, phi_{y3}}
  // expected ordering of gradient vector: eps^T = {kappa_x, kappa_y, kappa_{xy}, gamma_{yz}, gamma_{xy}}

  // kappa_x = kappa_{x1} + kappa_{x3} + kappa_{x3}
  gm[0][2]=dx[0];  // kappa_{x1} = dphi_{y1}/dx
  gm[0][5]=dx[3];  // kappa_{x2} = dphi_{y2}/dx
  gm[0][8]=dx[6];  // kappa_{x3} = dphi_{y3}/dx
  
  // kappa_y = kappa_{y1} + kappa_{y3} + kappa_{y3}
  gm[1][1]=0.0-dy[0];  // kappa_{y1} = -dphi_{x1}/dy
  gm[1][4]=0.0-dy[3];  // kappa_{y2} = -dphi_{x2}/dy
  gm[1][7]=0.0-dy[6];  // kappa_{y3} = -dphi_{x3}/dy

  // kappa_{xy} = kappa_{xy1} + kappa_{xy2} + kappa_{xy3}
  // kappa_{xy1} = dphi_{y1}/dy - dphi_{x1}/dx
  gm[2][1]=0.0-dx[0]; // -dphi_{x1}/dx
  gm[2][2]=dy[0];     //  dphi_{y1}/dy
  // kappa_{xy2} = dphi_{y2}/dy - dphi_{x2}/dx
  gm[2][4]=0.0-dx[3]; // -dphi_{x2}/dx
  gm[2][5]=dy[3];     //  dphi_{y2}/dy
  // kappa_{xy3} = dphi_{y3}/dy - dphi_{x3}/dx
  gm[2][7]=0.0-dx[6]; // -dphi_{x3}/dx
  gm[2][8]=dy[6];     //  dphi_{y3}/dy

  // gamma_{yz} = dw/dy - phi_{x} + dw_{sym}/dy
  gm[3][0]=dy[0];      
  gm[3][1]=dy[1]-bf[0];
  gm[3][2]=dy[2];
  gm[3][3]=dy[3];
  gm[3][4]=dy[4]-bf[3];
  gm[3][5]=dy[5];
  gm[3][6]=dy[6];
  gm[3][7]=dy[7]-bf[6];
  gm[3][8]=dy[8];
  
  // gamma_{xz} = dw/dx + phi_{y} + dw_{sym}/dx
  gm[4][0]=dx[0];
  gm[4][1]=dx[1];
  gm[4][2]=dx[2]+bf[0];
  gm[4][3]=dx[3];
  gm[4][4]=dx[4];
  gm[4][5]=dx[5]+bf[3];
  gm[4][6]=dx[6];
  gm[4][7]=dx[7];
  gm[4][8]=dx[8]+bf[6];
}

/**
   function extracts blocks from stiffness %matrix of the material
   
   @param ri,ci - row and column indices
   @param d - stiffness %matrix of the material
   @param dd - required block from stiffness %matrix of material
   @param t - thickness of plate
   
   JK, 19.7.2001
*/
void cctelem::dmatblock (long ri,long ci,matrix &d,matrix &dd,double t)
{
  double c;
  
  c=t*t*t;
  fillm (0.0,dd);
  
  if (ri==0 && ci==0){
    dd[0][0] = c*d[0][0];  dd[0][1] = c*d[0][1];  dd[0][2] = c*d[0][2];
    dd[1][0] = c*d[1][0];  dd[1][1] = c*d[1][1];  dd[1][2] = c*d[1][2];
    dd[2][0] = c*d[2][0];  dd[2][1] = c*d[2][1];  dd[2][2] = c*d[2][2];
  }
  if (ri==1 && ci==1){
    dd[0][0] = t*d[3][3]*5.0/6.0;  dd[0][1] = 0.0;
    dd[1][0] = 0.0;                dd[1][1] = t*d[4][4]*5.0/6.0;
  }

}

void cctelem::dmat (matrix &d,double t)
{
  double c;
  
  c=t*t*t;
  
  d[0][0] = c*d[0][0];  d[0][1] = c*d[0][1];  d[0][2] = c*d[0][2];
  d[1][0] = c*d[1][0];  d[1][1] = c*d[1][1];  d[1][2] = c*d[1][2];
  d[2][0] = c*d[2][0];  d[2][1] = c*d[2][1];  d[2][2] = c*d[2][2];
  
  d[3][3] = t*d[3][3]*5.0/6.0;
  d[4][4] = t*d[4][4]*5.0/6.0;
}

/**
   nutna kontrola

   function assembles transformation %matrix 
*/
void cctelem::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of cct element
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of node coordinates
   
   JK, 19.7.2001
*/
void cctelem::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii,jj,ipp;
  double jac,det,thick;
  ivector nodes(nne);
  vector l(3),t(nne),gp1,gp2,w,areacoord(3);
  matrix gmr,gmc,dd,d(tncomp,tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);  

  fillm (0.0,sm);
  

  for (ii=0;ii<nb;ii++){
    reallocm (ncomp[ii],ndofe,gmr);
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      
      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	areacoord[0]=gp1[i];
	areacoord[1]=gp2[i];
	areacoord[2]=1.0-areacoord[0]-areacoord[1];
	
	thick = approx (areacoord,t);

	//  geometric matrix
	geom_matrix_block (gmr,ii,x,y,areacoord);
	geom_matrix_block (gmc,jj,x,y,areacoord);
	
	//  matrix of stiffness of the material
	Mm->matstiff (d,ipp);
	dmatblock (ii,jj,d,dd,thick);
	
	jac=w[i]*det;
	
	//  contribution to the stiffness matrix of the element
	bdbjac (sm,gmr,dd,gmc,jac);
	
	ipp++;
      }
    }
  }

}

/**
   function assembles resulting stiffness %matrix of the element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 9.5.2002
*/
void cctelem::res_stiffness_matrix (long eid,matrix &sm)
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
   function computes mass %matrix of the cct element
   
   @param eid - element id
   @param mm - mass %matrix
   @param x,y - vectors of node coordinates
   
   JK, 19.7.2001
*/
void cctelem::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double det,ww,thick,rho;
  ivector nodes(nne);
  vector l(3),gp1(intordmm),gp2(intordmm),w(intordmm),t(nne),dens(nne);
  matrix n(3,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);
  
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,mm);

  for (i=0;i<intordmm;i++){
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
    ww=w[i];
    
    bf_matrix (n,x,y,l);
    
    thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];
    rho = dens[0]*l[0]+dens[1]*l[1]+dens[2]*l[2];
    
    ww*=det*thick*rho;

    nnj (mm.a,n.a,ww,3,ndofe);
  }
}

void cctelem::res_mass_matrix (long eid,matrix &mm)
{
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  mass_matrix (eid,mm,x,y);
}


/**
   function computes load %matrix of the cct element
   
   @param eid - element id
   @param lm - load %matrix
   
   25.7.2001
*/
void cctelem::load_matrix (long eid,matrix &lm)
{
  long i;
  double det,ww,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),l(3),gp1(intordmm),gp2(intordmm),w(intordmm),t(nne);
  matrix n(3,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);

  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
    ww=w[i];
    
    bf_matrix (n,x,y,l);
    
    thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];

    ww*=det*thick;
    
    nnj (lm.a,n.a,ww,3,ndofe);
  }
}

/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 26.9.2008
*/
void cctelem::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  vector gp1,gp2,w,eps(tncomp),areacoord(3);
  matrix gm(tncomp,ndofe);
  
  //  loop over blocks
  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      //  geometric matrix
      geom_matrix (gm,x,y,areacoord);
      
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
  }
}

/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void cctelem::res_ip_strains (long lcid,long eid)
{
  vector aux,x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
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
   function computes strains in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 26.9.2008
*/
void cctelem::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
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
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void cctelem::res_ip_stresses (long lcid,long eid)
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
   
   JK, 26.9.2008
*/
void cctelem::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
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


void cctelem::strains (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
  /*
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

  if (Mp->strainaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }
*/
}

/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   
   10.5.2002
*/
void cctelem::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 0.0;  eta[0] = 0.0;
  xi[1] = 1.0;  eta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 1.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void cctelem::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval(nne),areacoord(3);
  
  areacoord[0]=xi;  areacoord[1]=eta;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  k=0;
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (areacoord,nodval);
    k++;
  }
}


void cctelem::stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
/*
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

  if (Mp->stressaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
      delete [] stre[i];
    }
    delete [] stra;
    delete [] stre;
  }
*/
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
void cctelem::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void cctelem::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void cctelem::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void cctelem::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void cctelem::res_internal_forces (long lcid,long eid,vector &ifor)
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
void cctelem::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void cctelem::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void cctelem::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
void cctelem::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double thick;
  ivector nodes(nne);
  vector sig(5),gp1,gp2,w,l(3),t(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];

    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1){
	l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
	thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];

        Mm->computenlstresses (ipp,Mm->ip[ipp]);

	Mm->givestress (lcid,ipp,sig);
	sig[0]*=thick*thick*thick;
	sig[1]*=thick*thick*thick;
	sig[2]*=thick*thick*thick;
	sig[3]*=thick*5.0/6.0;
	sig[4]*=thick*5.0/6.0;
	Mm->storestress (lcid,ipp,sig);
	
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
void cctelem::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of local values
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
void cctelem::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;

  for (ii=0;ii<nb;ii++){
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
void cctelem::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct increments of stresses
      if (Mp->strcomp==1)
        Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void cctelem::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
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
   function integrates selected quantity over the finite element
   it results in nodal values
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y - node coordinates
   
   TKo 7.2008
*/
void cctelem::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double jac;
  vector w,gp1,gp2,areacoord(3);
  vector ipv,contr(ndofe);
  matrix gm;
  
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);

    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);

    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];

      jac_2d (jac,x,y,areacoord[0],areacoord[1]);
      reallocv (ncomp[ii],ipv);

      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);

      //  strain-displacement (geometric) matrix
      reallocm (ncomp[ii],ndofe,gm);
      geom_matrix_block (gm,ii,x,y,areacoord);
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      cmulv (jac*w[i],contr);

      //  summation
      addv(contr,nv,nv);
      ipp++;
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
void cctelem::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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


