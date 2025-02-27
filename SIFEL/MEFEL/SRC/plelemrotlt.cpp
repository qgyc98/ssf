#include "plelemrotlt.h"
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

planeelemrotlt::planeelemrotlt (void)
{
  long i,j;

  //  number nodes on element
  nne=3;
  //  number of DOFs on element
  ndofe=9;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
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

planeelemrotlt::~planeelemrotlt (void)
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
   
   @param areacoord - vector containing area coordinates
   @param nodval - nodal values
   
   6.1.2002
*/
double planeelemrotlt::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function approximates function defined by nodal values

   @param xi,eta - natural coordinates
   @param nodval - nodal values
   
   1.4.2002
*/
double planeelemrotlt::approx_nat (double xi,double eta,vector &nodval)
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
   function returns %matrix of base function
   
   @param n - array containing %matrix
   @param l - array of area coordinates
   @param b,c - arrays of coefficients of area coordinates

   6.1.2002
*/
void planeelemrotlt::bf_matrix (matrix &n,vector &l,vector &b,vector &c)
{
  vector bf(9);
  
  fillm (0.0,n);

  bf_rot_3_2d (bf.a,l.a,b.a,c.a);
  
  n[0][0]=bf[0];
  n[0][2]=bf[3];
  n[0][3]=bf[1];
  n[0][5]=bf[4];
  n[0][6]=bf[2];
  n[0][8]=bf[5];

  n[1][1]=bf[0];
  n[1][2]=bf[6];
  n[1][4]=bf[1];
  n[1][5]=bf[7];
  n[1][7]=bf[2];
  n[1][8]=bf[8];
}

/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param areacoord - area coordinates
   
   JK, 9.7.2001
*/
void planeelemrotlt::geom_matrix (matrix &gm,vector &x,vector &y,vector &areacoord)
{
  double area,det;
  vector b(3),c(3),dx(9),dy(9);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  area=det/2.0;
  
  b_coeff (b.a,y.a);
  c_coeff (c.a,x.a);

  dx_bf_rot_3_2d (dx.a,areacoord.a,b.a,c.a,area);
  dy_bf_rot_3_2d (dy.a,areacoord.a,b.a,c.a,area);

  fillm (0.0,gm);
  
  // expected ordering of nodal unknowns:  d^T = {u_1, v_1, phi_{z1}, u_2, v_2, phi_{z2}, u_3, v_3, phi_{z3}}
  // expected ordering of strain vector: eps^T = {eps_x, eps_y, gamma_{xy}}

  gm[0][0]=dx[0];
  gm[0][2]=dx[3];
  gm[0][3]=dx[1];
  gm[0][5]=dx[4];
  gm[0][6]=dx[2];
  gm[0][8]=dx[5];

  gm[1][1]=dy[0];
  gm[1][2]=dy[6];
  gm[1][4]=dy[1];
  gm[1][5]=dy[7];
  gm[1][7]=dy[2];
  gm[1][8]=dy[8];
  
  gm[2][0]=dy[0];
  gm[2][1]=dx[0];
  gm[2][2]=dy[3]+dx[6];
  gm[2][3]=dy[1];
  gm[2][4]=dx[1];
  gm[2][5]=dy[4]+dx[7];
  gm[2][6]=dy[2];
  gm[2][7]=dx[2];
  gm[2][8]=dy[5]+dx[8];
  
}


/**
   function assembles blocks of strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param ri - row index (number of required block)
   @param x,y - array containing node coordinates
   @param areacoord - area coordinates
   
   JK, 9.7.2001
*/
void planeelemrotlt::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,vector &areacoord)
{
  double area,det;
  vector b(3),c(3),dx(9),dy(9);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  area=det/2.0;
  
  b_coeff (b.a,y.a);
  c_coeff (c.a,x.a);
  
  dx_bf_rot_3_2d (dx.a,areacoord.a,b.a,c.a,area);
  dy_bf_rot_3_2d (dy.a,areacoord.a,b.a,c.a,area);

  fillm (0.0,gm);
  
  if (ri==0){
    gm[0][0]=dx[0];
    gm[0][2]=dx[3];
    gm[0][3]=dx[1];
    gm[0][5]=dx[4];
    gm[0][6]=dx[2];
    gm[0][8]=dx[5];
    
    gm[1][1]=dy[0];
    gm[1][2]=dy[6];
    gm[1][4]=dy[1];
    gm[1][5]=dy[7];
    gm[1][7]=dy[2];
    gm[1][8]=dy[8];
  }
  
  if (ri==1){
    gm[0][0]=dy[0];
    gm[0][1]=dx[0];
    gm[0][2]=dy[3]+dx[6];
    gm[0][3]=dy[1];
    gm[0][4]=dx[1];
    gm[0][5]=dy[4]+dx[7];
    gm[0][6]=dy[2];
    gm[0][7]=dx[2];
    gm[0][8]=dy[5]+dx[8];
  }
}

void planeelemrotlt::addgeommat (matrix &gm,vector &x,vector &y,vector &areacoord,double &jac)
{
  double det,area;
  vector b(3),c(3),bf(ndofe),dx(ndofe),dy(ndofe);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  jac = det;
  area=det/2.0;

  b_coeff (b.a,y.a);
  c_coeff (c.a,x.a);
  
  bf_rot_3_2d (bf.a,areacoord.a,b.a,c.a);
  dx_bf_rot_3_2d (dx.a,areacoord.a,b.a,c.a,area);
  dy_bf_rot_3_2d (dy.a,areacoord.a,b.a,c.a,area);
  
  
  fillm (0.0,gm);
  
  gm[0][0]=0.0-dy[0];
  gm[0][1]=dx[0];
  gm[0][2]=dx[6]-dy[3]-bf[0];
  gm[0][3]=0.0-dy[1];
  gm[0][4]=dx[1];
  gm[0][5]=dx[7]-dy[4]-bf[1];
  gm[0][6]=0.0-dy[2];
  gm[0][7]=dx[2];
  gm[0][8]=dx[8]-dy[5]-bf[2];

}

/**
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
   
   JK
*/
void planeelemrotlt::dmatblock (long ri,long ci,matrix &d, matrix &dd)
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
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix
   
   JK, 9.7.2001
*/
void planeelemrotlt::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*3][i*3]   = Mt->nodes[nodes[i]].e1[0];   tmat[i*3][i*3+1]   = Mt->nodes[nodes[i]].e2[0];  tmat[i*3][i*3+2]   = 0.0;
      tmat[i*3+1][i*3] = Mt->nodes[nodes[i]].e1[1];   tmat[i*3+1][i*3+1] = Mt->nodes[nodes[i]].e2[1];  tmat[i*3+1][i*3+2] = 0.0;
      tmat[i*3+2][i*3] = 0.0;                         tmat[i*3+2][i*3+1] = 0.0;                        tmat[i*3+2][i*3+2] = 1.0;
    }
  }
}

/**
   function computes stiffness %matrix of plane stress triangular
   finite element with rotational degrees of freedom with
   linear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of node coordinates

   8.12.2001
*/
void planeelemrotlt::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii,jj,ipp;
  double jac,det,thick;
  ivector nodes(nne);
  vector areacoord(3),t(nne),gp1,gp2,w;
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
      
      reallocm (ncomp[jj],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      reallocv (nip[ii][jj],gp1);
      reallocv (nip[ii][jj],gp2);
      reallocv (nip[ii][jj],w);
      
      
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];

      for (i=0;i<intordsm[ii][jj];i++){
        areacoord[0]=gp1[i];
        areacoord[1]=gp2[i];
        areacoord[2]=1.0-gp1[i]-gp2[i];
        
        // geometric matrix
        geom_matrix_block (gmr,ii,x,y,areacoord);
        geom_matrix_block (gmc,jj,x,y,areacoord);
	
        //  stiffness matrix of material
        Mm->matstiff (d,ipp);
        dmatblock (ii,jj,d,dd);
        
        thick = approx (areacoord,t);
        
        jac=thick*det*w[i];
        
        //  contribution to the stiffness matrix of the element
        bdbjac (sm,gmr,dd,gmc,jac);
	
        ipp++;
      }
    }
  }
  
  
  //long into = intordsm[0][0];
  long into = 1;
  reallocv (into,w);
  reallocv (into,gp1);
  reallocv (into,gp2);
  gauss_points_tr (gp1.a,gp2.a,w.a,into);
  
  for (i=0;i<into;i++){
    
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    reallocm (1,ndofe,gmr);
    reallocm (ndofe,ndofe,gmc);
    
    addgeommat (gmr,x,y,areacoord,jac);
    mtxm (gmr,gmr,gmc);
    thick = approx (areacoord,t);
    cmulm (thick*jac*w[i],gmc);
    addm (sm,gmc,sm);
  }

}

/**
   function computes stiffness %matrix of plane stress triangular
   finite element with rotational degrees of freedom with
   linear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of node coordinates

   8.12.2001
*/
void planeelemrotlt::stiffness_matrix2 (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii,jj,ipp;
  double jac,det,thick;
  ivector nodes(nne);
  vector areacoord(3),t(nne),gp1,gp2,w;
  matrix gm,gmr,gmc,dd,d(tncomp,tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  fillm (0.0,sm);

  reallocm (tncomp,ndofe,gm);

  ii=0;  jj=0;
  
  reallocv (nip[ii][jj],gp1);
  reallocv (nip[ii][jj],gp2);
  reallocv (nip[ii][jj],w);
  
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
  ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
  
  for (i=0;i<intordsm[ii][jj];i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    // geometric matrix
    geom_matrix (gm,x,y,areacoord);
    
    //  stiffness matrix of material
    Mm->matstiff (d,ipp);
    
    thick = approx (areacoord,t);
    
    jac=thick*det*w[i];
    
    //  contribution to the stiffness matrix of the element
    bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
    
    ipp++;
    
  }
  
  
  //long into = intordsm[0][0];
  long into = 1;
  reallocv (into,w);
  reallocv (into,gp1);
  reallocv (into,gp2);
  gauss_points_tr (gp1.a,gp2.a,w.a,into);
  
  for (i=0;i<into;i++){
    
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    reallocm (1,ndofe,gmr);
    reallocm (ndofe,ndofe,gmc);
    
    addgeommat (gmr,x,y,areacoord,jac);
    mtxm (gmr,gmr,gmc);
    thick = approx (areacoord,t);
    cmulm (thick*jac*w[i],gmc);
    addm (sm,gmc,sm);
  }

}


/**
   function assembles stiffness %matrix of element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK
*/
void planeelemrotlt::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  vector x(nne),y(nne);
  ivector nodes(nne);
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
   function computes mass %matrix of the plane stress triangular
   finite element with rotationale degrees of freedom wtih
   linear approximation functions
   
   @param eid - number of element
   @param mm - mass %matrix
   @param x,y - vectors of node coordinates
   
   8.12.2001
*/
void planeelemrotlt::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double jac,area,thick,rho;
  ivector nodes(nne);
  vector b(3),c(3),areacoord(3),w(intordmm),gp1(intordmm),gp2(intordmm),t(nne),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);

  //  det is equal to double area of the element
  area = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
  
  b_coeff (b.a,y.a);
  c_coeff (c.a,x.a);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++){
    
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    bf_matrix (n,areacoord,b,c);
    
    thick = approx (areacoord,t);
    rho = approx (areacoord,dens);
    
    jac=w[i]*thick*rho*area*2.0;
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function assembles mass %matrix of plane stress triangular
   finite element with rotational degrees of freedom
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK
*/

void planeelemrotlt::res_mass_matrix (long eid,matrix &mm)
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
   function computes load %matrix of the plane stress triangular
   finite element with linear approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values

   @param eid - number of element
   @param lm - load %matrix

   25.7.2001
*/
void planeelemrotlt::load_matrix (long eid,matrix &lm)
{
  long i;
  double jac,area,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp1(intordmm),gp2(intordmm),b(3),c(3),areacoord(3),t(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  area = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
  b_coeff (b.a,y.a);
  c_coeff (c.a,x.a);

  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    bf_matrix (n,areacoord,b,c);
    
    thick = approx (areacoord,t);
    
    //  zkontrolovat deleni dvema
    jac=w[i]*thick*area*2.0;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void planeelemrotlt::res_ip_strains (long lcid,long eid)
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
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 10.5.2002
*/
void planeelemrotlt::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  vector gp1,gp2,w,eps,epscum,areacoord(3);
  matrix gm;

  reallocm (tncomp,ndofe,gm);
  reallocv (tncomp,eps);
  
  //  loop over blocks
  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    //reallocv (ncomp[ii],eps);
    //reallocm (ncomp[ii],ndofe,gm);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-gp1[i]-gp2[i];
      
      // geometric matrix
      //geom_matrix_block (gmr,ii,x,y,areacoord);
      geom_matrix (gm,x,y,areacoord);
      
      mxv (gm,r,eps);
      
      //Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[0],eps);
      Mm->storestrain (lcid,ipp,0,tncomp,eps);
      ipp++;
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
void planeelemrotlt::nod_strains_ip (long lcid,long eid,long ri,long ci)
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
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void planeelemrotlt::ip_stresses (long lcid,long eid,long ri,long ci)
{
  compute_nlstress (lcid,eid,ri,ci);
}

/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void planeelemrotlt::res_ip_stresses (long lcid,long eid)
{
  long ri,ci;
  ri=0;
  ci=0;
  compute_nlstress (lcid,eid,ri,ci);
}

/**
   function computes stresses in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void planeelemrotlt::nod_stresses_ip (long lcid,long eid,long ri,long ci)
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
    //  strains at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storestress (lcid,0,sig);
  }
  
}

/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemrotlt::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
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

/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemrotlt::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
	Mm->computenlstressesincr (ipp);
      
      ipp++;
    }
  }
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
void planeelemrotlt::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);

}
























void planeelemrotlt::res_mainip_strains (long lcid,long eid)
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
   function computes strains in integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void planeelemrotlt::mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  vector gp1,gp2,w,eps,areacoord(3);
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
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-gp1[i]-gp2[i];
      
      geom_matrix_block (gm,ii,x,y,areacoord);
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
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
void planeelemrotlt::nod_strains (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,aux,natcoord(2),areacoord(3);
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
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-gp1[i]-gp2[i];
      
      Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
      
      natcoord[0]=areacoord[0];  natcoord[1]=areacoord[1];
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
   
   @param val - array containing strains on element
   @param lcid - load case id
   @param eid - element id
   
   15.7.2002
*/
void planeelemrotlt::elem_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,aux,natcoord(2),areacoord(3);
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
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
      
      natcoord[0]=areacoord[0];  natcoord[1]=areacoord[1];
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
   
   @param eid - element id
   @param xi, eta - natural coordinates of the point
   @param fi,li - first and last indices
   @param eps - array containing strains
   
   11.5.2002
*/
void planeelemrotlt::appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps)
{
  long i,j,k;
  ivector nodes(nne);
  vector nodval(nne),areacoord(3);
  
  if (ncomp != eps.n){
    fprintf (stderr,"\n\n wrong interval of indices in function strain (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];

  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].strain[lcid*tncomp+i];
    }
    eps[k]=approx (areacoord,nodval);
    k++;
  }
  
}

void planeelemrotlt::allip_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector eps(tncomp),gp1,gp2,w;
  
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

void planeelemrotlt::strains (long lcid,long eid,long ri,long ci)
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
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemrotlt::strains (%s, line %d).\n",__FILE__,__LINE__);
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
   
   10.5.2002
*/
void planeelemrotlt::nodecoord (vector &xi,vector &eta)
{
  xi[0] =  0.0;  eta[0] =  0.0;
  xi[1] =  1.0;  eta[1] =  0.0;
  xi[2] =  0.0;  eta[2] =  1.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void planeelemrotlt::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
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
   function computes stresses in integration points of element
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void planeelemrotlt::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector gp1,gp2,w,eps,sig,auxsig;
  matrix d(tncomp,tncomp),dd;

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
	
        dmatblock (ii,jj,d,dd);
        mxv (dd,eps,auxsig);
        addv (auxsig,sig,sig);
      }
      
      Mm->storestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
      
      ipp++;
    }
  }
}

void planeelemrotlt::nod_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,eps,sig,auxsig,natcoord(2);
  ivector nodes(nne);
  matrix d(tncomp,tncomp),dd;
  
  lsm = new double [9];

  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  
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
      ipp++;
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
        reallocv (ncomp[jj],eps);
        reallocm (ncomp[ii],ncomp[jj],dd);
	
	if (Mp->strainaver==0)
	  Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps);
	
        dmatblock (ii,jj,d,dd);
        mxv (dd,eps,auxsig);
        addv (auxsig,sig,sig);
      }
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,sig);
      
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->stress_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);
    
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

void planeelemrotlt::elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double *lsm,*lhs,*rhs;
  vector gp1,gp2,w,nxi(nne),neta(nne),eps,sig,auxsig,natcoord(2);
  ivector nodes(nne);
  matrix d(tncomp,tncomp),dd;

  lsm = new double [9];
  
  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  
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
      ipp++;
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (ncomp[jj],eps);
	reallocm (ncomp[ii],ncomp[jj],dd);
	
	if (Mp->strainaver==0)
	  appval (gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps,stra);
	if (Mp->strainaver==1)
	  appstrain (lcid,eid,gp1[i],gp2[i],cncomp[jj],ncomp[jj],eps);
	
	dmatblock (ii,jj,d,dd);
	mxv (dd,eps,auxsig);
	addv (auxsig,sig,sig);
      }
      
      natcoord[0]=gp1[i];  natcoord[1]=gp2[i];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,sig);
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    nodal_values (stre,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii]);
    
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

/**
   function computes stresses in arbitrary point on element
   
   @param eid - element id
   @param xi, eta - natural coordinates of the point
   @param fi,li - first and last indices
   @param sig - array containing stresses
   
   11.5.2002
*/
void planeelemrotlt::appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig)
{
  long i,j,k;
  ivector nodes(nne);
  vector areacoord(3),nodval(nne);
  
  if (ncomp != sig.n){
    fprintf (stderr,"\n\n wrong interval of indices in function stress (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }

  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];

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

void planeelemrotlt::allip_stresses (double **stre,long lcid,long eid,long ri,long ci)
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

void planeelemrotlt::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra=NULL,**stre=NULL;
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
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemrotlt::stresses (%s, line %d).\n",__FILE__,__LINE__);
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


/**
   function computes other values in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void planeelemrotlt::nod_others (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,j,ncomp,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp1,gp2,w,other,aux,natcoord(2),areacoord(3);
  ivector nodes(nne);

  lsm = new double [9];

  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);


  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    ncomp = Mm->ip[ipp].ncompother;
    reallocv (ncomp,other);
    lhs = new double [ncomp*3];
    rhs = new double [ncomp*3];
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    nullv (lsm,9);
    nullv (rhs,ncomp*3);
    
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-gp1[i]-gp2[i];
      
      for (j = 0; j < ncomp; j++)
        other[j] = Mm->ip[ipp].eqother[j];
      
      natcoord[0]=areacoord[0];  natcoord[1]=areacoord[1];
      matassem_lsm (lsm,natcoord);
      rhsassem_lsm (rhs,natcoord,other);
      
      ipp++;
    }
    
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp);
    Mt->other_nodal_values (nodes,nxi,neta,nxi,lhs,2,0,ncomp);
    
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}




/**
   function computes internal forces for nonlocal models

   this function is used in plane stress/strain elements (function is called
   by function res_nonloc_internal_forces) and shell elements

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void planeelemrotlt::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemrotlt::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void planeelemrotlt::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void planeelemrotlt::res_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemrotlt::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemrotlt::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemrotlt::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);

  Mt->give_node_coord2d (x,y,eid);
  eigstrain_forces (lcid,eid,0,0,nfor,x,y);
  Mt->give_node_coord2d (x,y,eid);
  internal_forces (lcid,eid,0,0,nfor,x,y);
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
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void planeelemrotlt::local_values (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  ivector nodes(nne);
  vector w,gp1,gp2,eps(tncomp);

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
void planeelemrotlt::compute_nonloc_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double **stra;
  ivector nodes(nne);
  vector w,gp1,gp2,eps(tncomp);

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
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void planeelemrotlt::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  vector eigstr(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
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
   
   TKo, 7.2008
*/
void planeelemrotlt::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double det,thick;
  ivector nodes(nne);
  vector areacoord(3),t(nne),ipv,contr(ndofe);
  vector w,gp1,gp2;
  matrix gm;

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));

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
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-gp1[i]-gp2[i];
      
      thick = approx (areacoord,t);
      
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);

      //  strain-displacement (geometric) matrix
      geom_matrix_block (gm,ii,x,y,areacoord);
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      cmulv (det*w[i]*thick,contr);
      
      //  summation
      addv(contr,nv,nv);      
      ipp++;
    }
  }
}


/**
   function computes nodal forces generated by continuous edge load
   
   @param eid - element id
   @param le - loaded edge indicator
   @param nv - array of nodal values
   @param nf - vector of nodal forces
   
   JK, 8.8.2015
*/
void planeelemrotlt::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double xi,eta,jac,det,area;
  vector x(nne),y(nne),gp(intordb),w(intordb),av(ndofe),v(ndofe),tnv(ndofe);
  vector b(3),c(3),l(3);
  matrix n(napfun,ndofe),am(ndofe,ndofe);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  area=det/2.0;
  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  if (le[0]>0){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=(1.0-gp[i])/2.0;  eta=(1.0+gp[i])/2.0;
      //bf_matrix (n,xi,eta);
      l[0]=xi;
      l[1]=eta;
      l[2]=1.0-xi-eta;

      bf_matrix (n,l,b,c);
    
      jac1d_2d (jac,x,y,gp[i],0);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[0]=nv[0];
    av[1]=nv[1];
    av[2]=nv[2];
    av[3]=nv[3];
    av[4]=nv[4];
    av[5]=nv[5];
    
    if (le[0]==2){
      //locglob_nodeval (0,av,tnv,x,y);
      fillv(0.0,av);
      
      av[0]=tnv[0];
      av[1]=tnv[1];
      av[2]=tnv[2];
      av[3]=tnv[3];
      av[4]=tnv[4];
      av[5]=tnv[5];
    }
    
    mxv (am,av,v);  addv (nf,v,nf);
  }

  if (le[1]>0){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=0.0;  eta=(1.0-gp[i])/2.0;
      //bf_matrix (n,xi,eta);
      l[0]=xi;
      l[1]=eta;
      l[2]=1.0-xi-eta;

      bf_matrix (n,l,b,c);
      
      jac1d_2d (jac,x,y,gp[i],1);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[3]=nv[6];
    av[4]=nv[7];
    av[5]=nv[8];
    av[6]=nv[9];
    av[7]=nv[10];
    av[8]=nv[11];
    
    if (le[1]==2){
      av[0]=nv[6];
      av[1]=nv[7];
      av[2]=nv[8];
      av[3]=nv[9];
      av[4]=nv[10];
      av[5]=nv[11];
      
      //locglob_nodeval (1,av,tnv,x,y);
      fillv(0.0,av);
      
      av[3]=tnv[0];
      av[4]=tnv[1];
      av[5]=tnv[2];
      av[6]=tnv[3];
      av[7]=tnv[4];
      av[8]=tnv[5];
    }
    mxv (am,av,v);  addv (nf,v,nf);
  }

  if (le[2]>0 ){
    fillm (0.0,am);
    for (i=0;i<intordb;i++){
      xi=(1.0+gp[i])/2.0;  eta=0.0;
      //bf_matrix (n,xi,eta);
      l[0]=xi;
      l[1]=eta;
      l[2]=1.0-xi-eta;

      bf_matrix (n,l,b,c);
      
      jac1d_2d (jac,x,y,gp[i],2);
      jac*=w[i];
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[6]=nv[12];
    av[7]=nv[13];
    av[8]=nv[14];
    av[0]=nv[15];
    av[1]=nv[16];
    av[2]=nv[17];
    
    if (le[2]==2){
      av[0]=nv[12];
      av[1]=nv[13];
      av[2]=nv[14];
      av[3]=nv[15];
      av[4]=nv[16];
      av[5]=nv[17];
      
      //locglob_nodeval (2,av,tnv,x,y);
      fillv(0.0,av);
      
      av[6]=tnv[0];
      av[7]=tnv[1];
      av[8]=tnv[2];
      av[0]=tnv[3];
      av[1]=tnv[4];
      av[2]=tnv[5];
    }
    
    mxv (am,av,v);  addv (nf,v,nf);
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
void planeelemrotlt::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  vector x(nne),y(nne),areacoord(3),w(intordsm[ri][ci]),gp1(intordsm[ri][ci]),gp2(intordsm[ri][ci]);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-areacoord[0]-areacoord[1];

    if (ii==ipp){
      coord[0]=approx (areacoord,x);
      coord[1]=approx (areacoord,y);
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
void planeelemrotlt::ipncoord (long eid, long ipp, vector &ncoord)
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
void planeelemrotlt::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
   07.2008 TKo - multiplictaion by thick added, correction of last ip point 
*/
void planeelemrotlt::ipvolume (long eid,long ri,long ci)
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
      
      reallocv (nip[ii][jj],gp1);
      reallocv (nip[ii][jj],gp2);
      reallocv (nip[ii][jj],w);
      
      
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];

      for (i=0;i<intordsm[ii][jj];i++){
	thick = approx_nat (gp1[i],gp2[i],t);
	jac_2d (jac,x,y,gp1[i],gp2[i]);
        jac*=w[i]*thick;
	
	Mm->storeipvol (ipp,jac);
	
        ipp++;
      }
    }
  }
  
  reallocv (1,w);
  reallocv (1,gp1);
  reallocv (1,gp2);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,1);
  
  thick = approx_nat (gp1[0],gp2[0],t);
  jac_2d (jac,x,y,gp1[0],gp2[0]);
  jac*=w[0]*thick;
  
  Mm->storeipvol (ipp,jac);
}
