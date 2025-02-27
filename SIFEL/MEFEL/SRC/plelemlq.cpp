#include "plelemlq.h"
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
#include "loadcase.h"
#include "gadaptivity.h"
#include <stdlib.h>
#include <math.h>



planeelemlq::planeelemlq (void)
{
  long i,j;
  
  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=8;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=2;
  //  order of numerical integration of mass matrix
  intordmm=2;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;
  
  //  number of blocks (parts of geometric matrix)
  //nb=2;
  nb=1;
  
  //  number of strain/stress components
  ncomp = new long [nb];
  if (nb==1){
    ncomp[0]=3;
  }
  if (nb==2){
    ncomp[0]=2;
    ncomp[1]=1;
  }
  
  //  cumulative number of components approximated
  cncomp = new long [nb];
  if (nb==1){
    cncomp[0]=0;
  }
  if (nb==2){
    cncomp[0]=0;
    cncomp[1]=2;
  }
  
  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  if (nb==1){
    nip[0][0]=4;
  }
  if (nb==2){
    nip[0][0]=4;  nip[0][1]=0;
    nip[1][0]=0;  nip[1][1]=1;
  }
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  if (nb==1){
    intordsm[0][0]=2;
  }
  if (nb==2){
    intordsm[0][0]=2;  intordsm[0][1]=0;
    intordsm[1][0]=0;  intordsm[1][1]=1;
  }
  
}

planeelemlq::~planeelemlq (void)
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

   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
   JK
*/
double planeelemlq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function assembles coordinates of selected integration point
   the integration point is defined by the parameter ipp
   
   @param eid - element id
   @param ipp - integration point id
   @param ri - row index
   @param ci - column index
   @param coord - array containing coordinates of selected integration point
   
   JK, 8.5.2002
*/
void planeelemlq::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,ii;
  double xi,eta;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordsm[ri][ci])),gp(ASTCKVEC(intordsm[ri][ci]));

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
void planeelemlq::ipncoord (long eid,long ipp,vector &ncoord)
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
   function assembles coordinates of integration points in block [ri][ci]
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param ipcoord - array containing coordinates of integration points
   
   JK, 8.5.2002
*/
void planeelemlq::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  long i,j,k;
  double xi,eta;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordsm[ri][ci])),gp(ASTCKVEC(intordsm[ri][ci]));
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  
  k=0;
  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ri][ci];j++){
      eta=gp[j];
      
      coord[k][0]=approx (xi,eta,x);
      coord[k][1]=approx (xi,eta,y);
      coord[k][2]=0.0;
      
      k++;
    }
  }
}

/**
   function assembles %matrix of base (approximation, shape) functions
   
   @param n - array containing %matrix
   @param xi, eta - natural coordinates
   
   JK, 9.7.2001
*/
void planeelemlq::bf_matrix (matrix &n,double xi,double eta)
{
  long i,j,k;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  nullm (n);

  j=0;  k=1;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];
    n[1][k]=bf[i];
    j+=2;  k+=2;
  }
}



/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemlq::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullm (gm);
  
  i1=0;  i2=1;
  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[1][i2]=dy[i];
    gm[2][i1]=dy[i];
    gm[2][i2]=dx[i];
    i1+=2;  i2+=2;
  }

}

/**
   function assembles blocks of strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param ri - row index (number of required block)
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemlq::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullm (gm);
  
  if (ri==0){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dx[i];
      gm[1][i2]=dy[i];
      i1+=2;  i2+=2;
    }
  }
  if (ri==1){
    i1=0;  i2=1;
    for (i=0;i<nne;i++){
      gm[0][i1]=dy[i];
      gm[0][i2]=dx[i];
      i1+=2;  i2+=2;
    }
  }
}

/**
   function assembles auxiliary vectors B for evaluation of stiffness %matrix
   in geometrically nonlinear problems
   
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   @param b11,b12,b21,b22 - vectors of derivatives of shape functions
   
   JK, 21.9.2005
*/
void planeelemlq::bvectors (vector &x,vector &y,double xi,double eta,double &jac,
			    vector &b11,vector &b12,vector &b21,vector &b22)
{
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  nullv (b11);
  nullv (b12);
  nullv (b21);
  nullv (b22);

  //  du/dx
  b11[0]=dx[0];  b11[2]=dx[1];  b11[4]=dx[2];  b11[6]=dx[3];
  //  du/dy
  b12[0]=dy[0];  b12[2]=dy[1];  b12[4]=dy[2];  b12[6]=dy[3];
  //  dv/dx
  b21[1]=dx[0];  b21[3]=dx[1];  b21[5]=dx[2];  b21[7]=dx[3];
  //  dv/dy
  b22[1]=dy[0];  b22[3]=dy[1];  b22[5]=dy[2];  b22[7]=dy[3];
}

/**
   function computes strain-displacement %matrix for geometrically nonlinear problems
   
   @param gm - strain-displacement %matrix
   @param r - array of nodal displacements
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian

   JK, 21.9.2005
*/
void planeelemlq::gngeom_matrix (matrix &gm,vector &r,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i;
  vector b11(ASTCKVEC(ndofe)),b12(ASTCKVEC(ndofe)),b21(ASTCKVEC(ndofe)),b22(ASTCKVEC(ndofe)),av(ASTCKVEC(ndofe));
  matrix am(ASTCKMAT(ndofe,ndofe));

  nullm (gm);
  
  bvectors (x,y,xi,eta,jac,b11,b12,b21,b22);
  
  // *******
  //  E_11
  // *******
  
  //  B11 dr
  for (i=0;i<ndofe;i++){
    gm[0][i]+=b11[i];
  }
  
  //  r B11 B11 dr
  vxv (b11,b11,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[0][i]+=av[i];
  }
  
  //  r B21 B21 dr
  vxv (b21,b21,am);
  vxm (r,am,av);
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
  
  //  r B22 B22 dr
  vxv (b22,b22,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  
  //  r B12 B12 dr
  vxv (b12,b12,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[1][i]+=av[i];
  }
  

  // **************
  //  E_12 = E_21
  // **************
  
  //  (B12 + B21) dr
  for (i=0;i<ndofe;i++){
    gm[2][i]+=b12[i]+b21[i];
  }
  
  //  r B11 B12 dr
  vxv (b11,b12,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B12 B11 dr
  vxv (b12,b11,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B22 B21 dr
  vxv (b22,b21,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  
  //  r B21 B22 dr
  vxv (b21,b22,am);
  vxm (r,am,av);
  for (i=0;i<ndofe;i++){
    gm[2][i]+=av[i];
  }
  


}


/**
   function computes gradient %matrix for geometrically nonlinear problems
   
   @param grm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian

   JK, 21.9.2005
*/
void planeelemlq::gnl_grmatrix (matrix &grm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i;
  vector b11(ASTCKVEC(ndofe)),b12(ASTCKVEC(ndofe)),b21(ASTCKVEC(ndofe)),b22(ASTCKVEC(ndofe));
  
  bvectors (x,y,xi,eta,jac,b11,b12,b21,b22);
  
  for (i=0;i<ndofe;i++){
    grm[0][i]=b11[i];
    grm[1][i]=b12[i];
    grm[2][i]=b21[i];
    grm[3][i]=b22[i];
  }
}





/**
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
   
   JK
*/
void planeelemlq::dmatblock (long ri,long ci,matrix &d, matrix &dd)
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
void planeelemlq::transf_matrix (ivector &nod,matrix &tmat)
{
  long i,n,m;
  
  nullm (tmat);

  n=nod.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nod[i]].transf>0){
      tmat[i*2][i*2]   = Mt->nodes[nod[i]].e1[0];    tmat[i*2][i*2+1]   = Mt->nodes[nod[i]].e2[0];
      tmat[i*2+1][i*2] = Mt->nodes[nod[i]].e1[1];    tmat[i*2+1][i*2+1] = Mt->nodes[nod[i]].e2[1];
    }
  }
}





/**
   function computes stiffness %matrix of plane rectangular
   finite element with bilinear approximation functions
   
   function computes stiffness %matrix for geometrically linear problems
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of nodal coordinates
   
   JK, 10.7.2001
*/
void planeelemlq::gl_stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ii,jj,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gmr,gmc,dd,d(ASTCKMAT(tncomp,tncomp));

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  nullm (sm);

  for (ii=0;ii<nb;ii++){
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gmr));
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));

      reallocm (RSTCKMAT(ncomp[jj],ndofe,gmc));
      reallocm (RSTCKMAT(ncomp[ii],ncomp[jj],dd));

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  
	  //  blocks of geometric matrices
	  if (nb==1){
	    geom_matrix (gmr,x,y,xi,eta,jac);
	    geom_matrix (gmc,x,y,xi,eta,jac);
	  }
	  else{
	    geom_matrix_block (gmr,ii,x,y,xi,eta,jac);
	    geom_matrix_block (gmc,jj,x,y,xi,eta,jac);
	  }
	  
	  //  matrix of stiffness of the material
	  Mm->matstiff (d,ipp);
	  if (nb==2){
	    //  block of the stiffness matrix of the material
	    dmatblock (ii,jj,d,dd);
	  }
	  
	  //  thickness in integration point
	  thick = approx (xi,eta,t);
	  
	  jac*=thick*w[i]*w[j];
	  
	  if (jac<0.0){
            print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	    jac=fabs(jac);
	  }
	  
	  //  contribution to the stiffness matrix of the element
	  if (nb==1){
	    bdbjac (sm,gmr,d,gmc,jac);
	  }
	  if (nb==2){
	    bdbjac (sm,gmr,dd,gmc,jac);
	  }
	  
	  ipp++;
	}
      }
    }
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
void planeelemlq::gnl_stiffness_matrix (long lcid,long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,jac,jac2,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne)),sig(ASTCKVEC(tncomp)),r(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe)),grm(ASTCKMAT(4,ndofe)),d(ASTCKMAT(tncomp,tncomp)),s(ASTCKMAT(4,4));
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of element
  Mc->give_thickness (eid,nodes,t);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  component setting to zero
  nullm (sm);

  //  array for weights of integration points
  reallocv (RSTCKVEC(intordsm[0][0],w));
  //  array for coordinates of integration points
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //
      //  linear stiffness matrix and inital deformation matrix
      //
      
      //  strain-displacement matrix
      gngeom_matrix (gm,r,x,y,xi,eta,jac);
      
      //  stiffness matrix of the material
      Mm->matstiff (d,ipp);
      
      //  thickness in integration point
      thick = approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
      
      //  contribution to the stiffness matrix of the element
      bdbjac (sm,gm,d,gm,jac);
      
      
      //
      //  initial stress matrix
      //
      
      //  gradient matrix
      gnl_grmatrix (grm,x,y,xi,eta,jac2);
      
      //  stresses
      Mm->givestress (lcid,ipp,sig);
      
      s[0][0]=sig[0];  s[0][1]=sig[2];
      s[1][0]=sig[2];  s[1][1]=sig[1];

      s[2][2]=sig[0];  s[2][3]=sig[2];
      s[3][2]=sig[2];  s[3][3]=sig[1];

      //  contribution to the stiffness matrix of the element
      bdbjac (sm,grm,s,grm,jac);
      
      ipp++;
    }
  }
}

/**
   function computes stiffness %matrix of quadrilateral element
   
   @param lcid - load case id
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK
*/
void planeelemlq::res_stiffness_matrix (long /*lcid*/,long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));

  Mt->give_node_coord2d (x,y,eid);

  
  gl_stiffness_matrix (eid,0,0,sm,x,y);
  //gnl_stiffness_matrix (lcid,eid,0,0,sm,x,y);
  
  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }

}



/**
   function computes B^T D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - B^T D %matrix

   23/10/2019 TKr according to TKo lintet.cpp
*/
void planeelemlq::bd_matrix (long eid,long ri,long ci,matrix &mbd)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp)),bd;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  nullm (mbd);
  reallocm (RSTCKMAT(mbd.m,mbd.n,bd));
  
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  only for one block of int. points
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  only one block of geometric matrices
      geom_matrix (gm,x,y,xi,eta,jac);
      
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      //  thickness in integration point
      thick = approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
      
      if (jac<0.0){
	print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	jac=fabs(jac);
      }
      
      //  contribution to the stiffness matrix of the element
      //  B^T D
      mtxm (gm,d,bd);
      //  B^T D Jac
      cmulm (jac,bd,bd);
      addm(mbd,bd,mbd);
      ipp++;
    }
  }
}

/**
   function computes integral of D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - integral of D %matrix

   23/10/2019 TKr according to TKo lintet.cpp
*/
void planeelemlq::dd_matrix (long eid,long ri,long ci,matrix &mdd)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  nullm (mdd);
  
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  only for one block of int. points
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  only one block of geometric matrices
      geom_matrix (gm,x,y,xi,eta,jac);

      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      //  thickness in integration point
      thick = approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
      
      if (jac<0.0){
	print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	jac=fabs(jac);
      }
      
      //  contribution to the stiffness matrix of the element
      //  D Jac
      cmulm (jac,d,d);
      addm(mdd,d,mdd);
      ipp++;
    }
  }
}



/**
   function computes mass %matrix of the plane stress rectangular
   finite element with bilinear approximation functions
   
   this function is used in plane stress/strain elements (function is called
   by function res_mass_matrix) and shell elements

   @param eid - number of element
   @param mm - mass %matrix
   @param x,y - vectors of nodal coordinates
   
   JK, 24.6.2001
*/
void planeelemlq::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,thick,rho;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne)),dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (mm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      thick = approx (xi,eta,t);
      rho = approx (xi,eta,dens);
      jac*=w[i]*w[j]*thick*rho;
      
      nnj (mm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function assembles mass %matrix of plane stress rectangular
   finite element with bilinear approximation functions
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK
*/
void planeelemlq::res_mass_matrix (long eid,matrix &mm)
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
   function computes damping %matrix of plane rectangular
   finite element with bilinear approximation functions
   
   function computes stiffness %matrix for geometrically linear problems
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of nodal coordinates
   
   JK, 20.11.2017
*/
void planeelemlq::damping_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ii,jj,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gmr,gmc,dd,d(ASTCKMAT(tncomp,tncomp));

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  nullm (sm);

  for (ii=0;ii<nb;ii++){
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gmr));
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));

      reallocm (RSTCKMAT(ncomp[jj],ndofe,gmc));
      reallocm (RSTCKMAT(ncomp[ii],ncomp[jj],dd));

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];

	  //  blocks of geometric matrices
	  geom_matrix_block (gmr,ii,x,y,xi,eta,jac);
	  geom_matrix_block (gmc,jj,x,y,xi,eta,jac);

	  //  matrix of damping of the material
	  Mm->matdamp (d,ipp);

	  //  block of the damping matrix of the material
	  dmatblock (ii,jj,d,dd);

	  //  thickness in integration point
	  thick = approx (xi,eta,t);

	  jac*=thick*w[i]*w[j];
	  
	  jac=fabs(jac);

	  //  contribution to the damping matrix of the element
	  bdbjac (sm,gmr,dd,gmc,jac);
	  
	  ipp++;
	}
      }
    }
  }
}

/**
   function assembles damping %matrix of plane stress rectangular
   finite element with bilinear approximation functions
   
   @param eid - element id
   @param mm - damping %matrix
   
   JK, 20.11.2017
*/
void planeelemlq::res_damping_matrix (long eid,matrix &dm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);

  damping_matrix (eid,0,0,dm,x,y);
  
  //  transformation of damping matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (dm,tmat);
  }
}



/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void planeelemlq::res_ip_strains (long lcid,long eid)
{
  vector aux,x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;

  Mt->give_elemnodes (eid,nodes);
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  Mt->give_node_coord2d (x,y,eid);
  gl_ip_strains (lcid,eid,0,0,x,y,r);
  //gnl_ip_strains (lcid,eid,0,0,x,y,r);
  
}

/**
   function computes strains at integration points of element
   function is used in geometrically linear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 10.5.2002
*/
void planeelemlq::gl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ii,ipp;
  double xi,eta,jac;
  vector gp,w,eps,epscum;
  matrix gm;
  
  //  loop over blocks
  for (ii=0;ii<nb;ii++){
    
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(ncomp[ii],eps));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	if (nb==1){
	  geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (nb==2){
	  geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	}
	
	mxv (gm,r,eps);
	
	Mm->storestrain (lcid,ipp,cncomp[ii],eps);
	ipp++;
      }
    }
  }
  
  if (nb==2){
    // ********************************************
    //  exchange of values at integration points
    // ********************************************
    //  values from 4 points to 1 point
    ipp=Mt->elements[eid].ipp[ri][ci];
    reallocv (RSTCKVEC(ncomp[0],epscum));
    reallocv (RSTCKVEC(ncomp[0],eps));
    for (i=0;i<intordsm[0][0];i++){
      for (j=0;j<intordsm[0][0];j++){
	Mm->givestrain (lcid,ipp,cncomp[0],eps);
	addv (eps,epscum,epscum);
	ipp++;
      }
    }
    
    cmulv (1.0/4.0,epscum,epscum);
    ipp=Mt->elements[eid].ipp[ri+1][ci+1];
    Mm->storestrain (lcid,ipp,epscum);
    
    //  value from 1 point to 4 points
    reallocv (RSTCKVEC(ncomp[1],eps));
    Mm->givestrain (lcid,ipp,cncomp[1],eps);
    ipp=Mt->elements[eid].ipp[ri][ci];
    Mm->storestrain (lcid,ipp,cncomp[1],eps);  ipp++;
    Mm->storestrain (lcid,ipp,cncomp[1],eps);  ipp++;
    Mm->storestrain (lcid,ipp,cncomp[1],eps);  ipp++;
    Mm->storestrain (lcid,ipp,cncomp[1],eps);  ipp++;
  }
  
}

/**
   function computes strains at integration points of element
   function is used in geometrically linear problems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 24.9.2005
*/
void planeelemlq::gnl_ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ipp;
  double xi,eta,jac,b11r,b12r,b21r,b22r;
  vector gp,w,eps,b11(ASTCKVEC(ndofe)),b12(ASTCKVEC(ndofe)),b21(ASTCKVEC(ndofe)),b22(ASTCKVEC(ndofe));
  
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(tncomp,eps));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      bvectors (x,y,xi,eta,jac,b11,b12,b21,b22);

      scprd (b11,r,b11r);
      scprd (b12,r,b12r);
      scprd (b21,r,b21r);
      scprd (b22,r,b22r);
      
      eps[0]=b11r+0.5*b11r*b11r+0.5*b21r*b21r;
      eps[1]=b22r+0.5*b12r*b12r+0.5*b22r*b22r;
      eps[2]=b12r+b21r+b11r*b12r+b21r*b22r;
      
      
      Mm->storestrain (lcid,ipp,eps);
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
void planeelemlq::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i, j, ipp;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double area;
  vector eps(ASTCKVEC(gncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes(eid, nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain(lcid, ipnum[i], eps);
    
    j=nod[i];
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.25;
      cmulv(area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}

/**
  The function computes nodal strains directly, averageing of nodal strains is performed according to setup.
   
  @param lcid - load case id
  @param eid - element id
   
  JK, 25.9.2004
*/
void planeelemlq::nod_strains_comp (long lcid, long eid)
{
  long i, j;
  double jac, area;
  ivector enod(ASTCKIVEC(nne)), ipnum(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  matrix tmat, gm(ASTCKMAT(tncomp, ndofe));
  
  //  node coordinates
  Mt->give_node_coord2d(x, y, eid);
  //  node numbers
  Mt->give_elemnodes(eid, enod);
  //  nodal displacements
  eldispl(lcid, eid, r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems(enod);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(enod, tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planelq(nxi, neta);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix(gm, x, y, nxi[i], neta[i], jac);
    //  strain computation
    mxv(gm, r, eps);
    
    j=enod[i];
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.25;
      cmulv(area,eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes nodal strains directly with no averageing.
   
  @param lcid - load case id
  @param eid - element id
   
  JK, 25.9.2004
*/
void planeelemlq::nod_strains (long lcid,long eid)
{
  long i, j, k, m;
  double jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  ivector enod(ASTCKIVEC(nne));
  matrix tmat, gm(ASTCKMAT(tncomp, ndofe));
  
  //  node coordinates
  Mt->give_node_coord2d (x, y, eid);
  //  node numbers
  Mt->give_elemnodes (eid, enod);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe,aux));
    reallocm(RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planelq (nxi,neta);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,nxi[i],neta[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    
    //  storage of strains to the node
    j=enod[i];
    m = lcid*tncomp;
    for (k=0; k<tncomp; k++)
      Mt->nodes[j].strain[m+k] = eps(k);
  }
}



/**
   function computes strains at strain points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void planeelemlq::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  vector coord,eps;

  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //  strains are computed at integration points
    res_ip_strains (lcid,eid);
    break;
  }
  case enodes:{
    //  strains are copied to nodes from the closest integration points
    nod_strains_ip (lcid,eid,ri,ci);
    break;
  }
  case cenodes:{
    //  strains are computed at nodes
    nod_strains_comp (lcid,eid);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (RSTCKVEC(ncp,eps));
    reallocv (RSTCKVEC(2,coord));
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);
      
      if (Mp->strainaver==0)
	//appval (coord[0],coord[1],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	//appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);
      
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    print_err("unknown strain point is required", __FILE__, __LINE__, __func__);
  }
  }
  
}


/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void planeelemlq::res_ip_stresses (long lcid, long eid)
{
  long ri,ci;
  ri=0;
  ci=0;
  compute_nlstress(lcid, eid, ri, ci);
}



/**
   function computes stresses at nodes of element

   @param lcid[in] - load case id
   @param eid[in] - element id
   @param ri,ci[in] - row and column indices
   
   JK, 10.5.2002
*/
void planeelemlq::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i, j, ipp;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double area;
  vector sig(ASTCKVEC(gncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
      Mm->givestress(lcid, ipnum[i], sig);
    
    //  storage of stresses to the node
    j=nod[i];
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress(lcid, 0, sig);
    if (Mp->stressaver==2)
    {
      area = Mt->give_area(eid)*0.25;
      cmulv(area, sig);
      Mt->nodes[j].storestress(lcid, 0, area, sig);
    }
  }
}



/**
  The function computes nodal stresses directly. Nodal stress values are avreaged according
  to the setup.
   
  @param lcid[in] - load case id
  @param eid[in] - element id
  @param stra[in] - array for strain components, stra[i]=pointer to the array of strain components at the i-th node
  @param stre[out] - array for stress components, stre[i]=pointer to the array of stress components at the i-th node
   
  JK, 25.9.2004
*/
void planeelemlq::nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre)
{
  long i, j, ipp;
  double area;
  vector eps(ASTCKVEC(tncomp)), sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  
  //  number of the first integration point on the element
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);

  //  loop over nodes
  for (i=0;i<nne;i++){
    for (j=0;j<eps.n;j++){
      eps[j]=stra[i][j];
    }
    //  stiffness matrix of the material
    Mm->matstiff(d, ipnum[i]);
    //  stress computation
    mxv (d,eps,sig);
    for (j=0;j<eps.n;j++){
      stre[i][j]=sig[j];
    }
    j=nod[i];
    //  storage of strains to the node
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid, 0, sig);
    if (Mp->stressaver==2)
    {
      area = Mt->give_area(eid)*0.25;
      cmulv(area, sig);
      Mt->nodes[j].storestress(lcid, 0, area, sig);
    }
  }
}



void planeelemlq::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  vector coord,sig;
  
  if (Mp->stressaver==0){
    /*
    double **stra,**stre;
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    elem_strains (stra,lcid,eid,ri,ci);
    elem_stresses (stra,stre,lcid,eid,ri,ci);
    */
  }
  //  computation of blocks of stresses at integration points
  res_ip_stresses (lcid,eid);
  
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    res_ip_stresses (lcid,eid);
    break;
  }
  case enodes:{
    nod_stresses_ip (lcid,eid,ri,ci);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (RSTCKVEC(ncp,sig));
    reallocv (RSTCKVEC(2,coord));
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],coord[1],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);

      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    print_err("unknown stress point is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}




/**
   function computes other values in nodes of element

   @param eid - element id
   
   10.5.2002
*/
void planeelemlq::nod_other_ip (long eid,long ri,long ci)
{
  long i, j, ncompo, ipp;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double area;
  vector other;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq (ipp, intordsm[0][0], ipnum);
 
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++)
  {
    ncompo = Mm->givencompother(ipnum[i], 0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    
    j=nod[i];
    //  storage of other values to the node
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother(0, ncompo, other);
    if (Mp->otheraver==2)
    {
      area = Mt->give_area(eid)*0.25;
      cmulv(area, other);
      Mt->nodes[j].storeother(0, ncompo, area, other);
    }    
  }
}





/**
   function computes load %matrix of the plane stress rectangular
   finite element with bilinear approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   @param x,y - node coordinates
   
   JK, 25.7.2001
*/
void planeelemlq::load_matrix (long eid,matrix &lm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  gauss_points (gp.a,w.a,intordmm);
  
  nullm (lm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      thick = approx (xi,eta,t);
      jac*=w1*w2*thick;
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function computes load %matrix of the plane stress rectangular
   finite element with bilinear approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 25.7.2001
*/
void planeelemlq::res_load_matrix (long eid,matrix &lm)
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













/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   JK, 28.7.2001
*/
void planeelemlq::gl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   JK, 22.9.2005
*/
void planeelemlq::gnl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  long i,j,k,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne)),sig(ASTCKVEC(tncomp)),contr(ASTCKVEC(ndofe)),r(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  //  element nodeds
  Mt->give_elemnodes (eid,nodes);
  //  thickness of element
  Mc->give_thickness (eid,nodes,t);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  nullv (ifor);
  
  //  array for coordinates of integration points
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  //  array for weights of integration points
  reallocv (RSTCKVEC(intordsm[0][0],w));
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  thickness at integration point
      thick = approx (xi,eta,t);
      
      //  computation of stress
      if (Mp->strcomp==1)
	Mm->computenlstresses (ipp,Mm->ip[ipp]);
      
      Mm->givestress (lcid,ipp,sig);
      
      //  strain-displacement (geometric) matrix
      gngeom_matrix (gm,r,x,y,xi,eta,jac);
      
      mtxv (gm,sig,contr);
      
      cmulv (jac*w[i]*w[j]*thick,contr);
      
      for (k=0;k<contr.n;k++){
	ifor[k]+=contr[k];
      }
      
      ipp++;
    }
  }
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
void planeelemlq::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemlq::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemlq::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
   
   JK, 22.9.2005
*/
void planeelemlq::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  
  gl_internal_forces (lcid,eid,0,0,ifor,x,y);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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
   
   JK, 22.9.2005
*/
void planeelemlq::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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
void planeelemlq::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  
  incr_internal_forces (lcid,eid,0,0,ifor,x,y);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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

   JK, 28.7.2001
*/
void planeelemlq::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  
  eigstrain_forces (lcid,eid,0,0,nfor,x,y);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    //globloctransf (nfor,v,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   Function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 10.6.2013
*/
void planeelemlq::local_values(long /*lcid*/,long eid,long ri,long ci)
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
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemlq::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
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
   
   TKo  7.2008
*/
void planeelemlq::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      for (j=0;j<intordsm[ii][ii];j++){
	//  computation of correct increments of stresses
	if (Mp->strcomp==1)
	  Mm->computenlstressesincr (ipp);
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
   
   JK, 27.11.2006
*/
void planeelemlq::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
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
   function computes correct eigen stresses caused by eigenstrains at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemlq::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
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
   
   JK, 27.11.2006
*/
void planeelemlq::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,j,ii,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne)),ipv,contr(ASTCKVEC(ndofe));
  matrix gm;
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  nullv (nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));
    reallocv (RSTCKVEC(ncomp[ii],ipv));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	thick = approx (xi,eta,t);
	
	//  function assembles required quantity at integration point
	Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
	
	if (nb==1){
	  //  strain-displacement (geometric) matrix
	  geom_matrix (gm,x,y,xi,eta,jac);
	}
	else{
	  //  strain-displacement (geometric) matrix
	  geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	}
	
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
   The function integrates arbitrary selected quantity over the finite element, e.g.
   it performs \int_{\Omega} \mbf{\sigma} d\Omega which results in integrated values that can 
   be used in the homogenization problems.

   
   @param eid - element id
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param iv - integrated values (output)
   
   @return The function returns nodal values calculated in the %vector iv
   
   23/10/2019 TKr according to TKo
*/
void planeelemlq::elem_volintegration_quant (long eid, integratedquant iq, long lcid, vector &iv)
{
  long i,j,ii,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gm;
  vector ipv(ASTCKVEC(iv.n)), contr(ASTCKVEC(iv.n));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d(x, y, eid);
  nullv(iv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));
    reallocv (RSTCKVEC(ncomp[ii],ipv));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ii][ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	thick = approx (xi,eta,t);
	
	//  function assembles required quantity at integration point
	Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
	
	if (nb==1){
	  //  strain-displacement (geometric) matrix
	  geom_matrix (gm,x,y,xi,eta,jac);
	}
	else{
	  //  strain-displacement (geometric) matrix
	  geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	}

	//  summation of contributions to the volume integral of the given quantity
	cmulv (jac*w[i]*w[j]*thick,ipv,contr);
	
	//  summation
	addv(contr,iv,iv);
	
	ipp++;
      }
    }
  }
}


void planeelemlq::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double ww,jac,xi,eta;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),av(ASTCKVEC(ndofe)),v(ASTCKVEC(ndofe));
  matrix n(ASTCKMAT(napfun,ndofe)),am(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);

  if (le[0]==1){
    nullm (am);
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,0);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[0]=nv[0];  av[1]=nv[1];  av[2]=nv[2];  av[3]=nv[3];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    nullm (am);
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];

      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,1);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[2]=nv[4];  av[3]=nv[5];  av[4]=nv[6];  av[5]=nv[7];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    nullm (am);
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,2);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[4]=nv[8];  av[5]=nv[9];  av[6]=nv[10];  av[7]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    nullm (am);
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,3);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    av[6]=nv[12];  av[7]=nv[13];  av[0]=nv[14];  av[1]=nv[15];
    mxv (am,av,v);  addv (nf,v,nf);
  }
}


/**
   function integrates N^T c N over edges
   
   @param edg - edge id (number of edge)
   @param x, y - coordinates of element nodes
   @param intord - order of numerical integration
   @param gp, w  - coordinates and weights of integration points
   @param t - nodal thicknesses
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK, 20.10.2007
*/
void planeelemlq::edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
				 vector &t,vector &coef,matrix &km)
{
  long i;
  double xi,eta,jac,ipval,thick;
  matrix n(ASTCKMAT(napfun,ndofe));
  
  if (edg==0){
    eta=1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,xi,edg);
      
      ipval=approx (xi,eta,coef);
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==1){
    xi=-1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,eta,edg);
      
      ipval=approx (xi,eta,coef);
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==2){
    eta=-1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,xi,edg);
      
      ipval=approx (xi,eta,coef);
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==3){
    xi=1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,eta,edg);
      
      ipval=approx (xi,eta,coef);
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}




/**
   function picks up nodal values on required edges
   
   @param edg - number of required edge
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all edges
   
   JK, 19.8.2004
*/
void planeelemlq::edgenodeval (long edg,vector &nodval,double *list)
{
  long i,j,k;
  ivector edgenod(ASTCKIVEC(nned));
  
  nullv (nodval);
  linquadrilat_edgnod (edgenod.a,edg);
  
  k=0;
  for (i=0;i<nned;i++){
    for (j=0;j<napfun;j++){
      nodval[k]=list[edg*nned*napfun+k];
      k++;
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
void planeelemlq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, ipval;
  vector w, gp, anv(ASTCKVEC(nne));
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
        reallocv (RSTCKVEC(intordsm[ii][jj],gp));
        reallocv (RSTCKVEC(intordsm[ii][jj],w));
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
   
   @param eid - element id
   @param ri,ci - row and column indices
   
   2.3.2004, JK
   07.2008 TKo - multiplictaion by thickness added
*/
void planeelemlq::ipvolume (long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),w,gp;
  
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
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

/**
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   @param nodval - array of nodal values
   @param ipval - array of values at integration points
   
   21.6.2004, JK
*/
void planeelemlq::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,k;
  double xi,eta;
  vector w,gp;
  
  k=0;
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	
	ipval[k]=approx (xi,eta,nodval);
	k++;
      }
    }
  }
}





/**
   function defines meaning of DOFs at nodes
   
   @param eid - element id
   
   JK, 21.8.2005
*/
void planeelemlq::define_meaning (long eid)
{
  ivector cn(ASTCKIVEC(ndofe)),nod(ASTCKIVEC(nne));
  
  Mt->give_elemnodes (eid,nod);
  Mt->give_code_numbers (eid,cn.a);

  //  displacement in x direction
  if (cn[0]>0)  Mt->nodes[nod[0]].meaning[0]=1;
  //  displacement in y direction
  if (cn[1]>0)  Mt->nodes[nod[0]].meaning[1]=2;
  //  displacement in x direction
  if (cn[2]>0)  Mt->nodes[nod[1]].meaning[0]=1;
  //  displacement in y direction
  if (cn[3]>0)  Mt->nodes[nod[1]].meaning[1]=2;
  //  displacement in x direction
  if (cn[4]>0)  Mt->nodes[nod[2]].meaning[0]=1;
  //  displacement in y direction
  if (cn[5]>0)  Mt->nodes[nod[2]].meaning[1]=2;
  //  displacement in x direction
  if (cn[6]>0)  Mt->nodes[nod[3]].meaning[0]=1;
  //  displacement in y direction
  if (cn[7]>0)  Mt->nodes[nod[3]].meaning[1]=2;

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
void planeelemlq::ntdbr_vector (long eid,vector &ntdbr)
{
  long intord = 2;
  long i,j,k,l,ipp,ippsh,ri,ci,lcid,body;
  double thick,xi,eta,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intord)),w(ASTCKVEC(intord)),t(ASTCKVEC(nne));
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp)),bf(ASTCKVEC(nne)),r;
  matrix d(ASTCKMAT(tncomp,tncomp)),gm;

  ri = ci = lcid = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  body = 1;
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);

  Mm->matstiff (d,ipp);

  //if (body){
    ippsh = ipp+4;
  // }
  // else{
  //   lcid = 0;
  //   reallocv (RSTCKIVEC(ndofe,cn));
  //   reallocv (RSTCKVEC(ndofe,r));
  //   reallocm (RSTCKMAT(tncomp,ndofe,gm));
  //   
  //   Mt->give_code_numbers (eid,cn.a);
  //   eldispl (lcid,eid,r.a,cn.a,ndofe);
  //   
  //   geom_matrix (gm,x,y,0.0,0.0,jac);
  //   mxv (gm,r,eps);
  //   sh = eps[2];
  // }
  
  gauss_points (gp.a,w.a,intord);
  
  nullv (ntdbr);
  
  for (i=0;i<intord;i++){
    xi=gp[i];
    for (j=0;j<intord;j++){
      eta=gp[j];
      
      eps[0] = Mm->ip[ipp  ].strain[0];
      eps[1] = Mm->ip[ipp  ].strain[1];
      eps[2] = Mm->ip[ippsh].strain[2];
      ipp++;
      jac_2d (jac,x,y,xi,eta);
      
      mxv (d,eps,sig);
      
      bf_lin_4_2d (bf.a,xi,eta);
      
      thick = approx (xi,eta,t);
      jac*=thick*w[i]*w[j];

      for (k=0;k<tncomp;k++)
        for (l=0;l<nne;l++)
          ntdbr[k*nne+l] += jac * bf[l] * sig[k];
      
    }
  }
}

/**
   Function integrates function "NT*N" over whole element.
   N is matrix of basic functions.
   !!! Matrix N is composed of three same blocks, it is computed only one third.
   
   @param eid - element id
   @param ntn - empty(returned) 2Darray, dimension is nne x nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemlq::ntn_matrix (long eid,matrix &ntn)
{
  long intord = 2; //intordmm;
  long i,j,k,l;
  double thick,xi,eta,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),gp(ASTCKVEC(intord)),w(ASTCKVEC(intord)),bf(ASTCKVEC(nne));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  gauss_points (gp.a,w.a,intord);
  
  nullm (ntn);
  
  for (i=0;i<intord;i++){
    xi=gp[i];
    for (j=0;j<intord;j++){
      eta=gp[j];
      
      thick = approx (xi,eta,t);
      jac_2d (jac,x,y,xi,eta);
      jac*=thick*w[i]*w[j];
      
      bf_lin_4_2d (bf.a,xi,eta);
      
      for (k=0;k<nne;k++)
        for (l=0;l<nne;l++)
          ntn[k][l] += jac * bf[k] * bf[l]; 
      
    }
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
double planeelemlq :: compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull)
{
  long intord = 2;
  long i,j,ipp,ippsh,ri,ci;
  double thick,area,xi,eta,jac,contr;
  double zero=1.0e-20;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),gp(ASTCKVEC(intord)),w(ASTCKVEC(intord)),bf(ASTCKVEC(nne)),r;
  vector sig_star(ASTCKVEC(tncomp)),sig_roof(ASTCKVEC(tncomp)),sig_err(ASTCKVEC(tncomp)),eps(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp)),dinv(ASTCKMAT(tncomp,tncomp)),gm;
  
  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  //body = Ada->body[5];
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  
  Mm->matstiff (d,ipp);
  invm (d,dinv,zero);

  //if (body){
    ippsh = ipp+4;
    //}
  //else{
  //  lcid = 0;
  //  reallocv (RSTCKIVEC(ndofe,cn));
  //  reallocv (RSTCKVEC(ndofe,r));
  //  reallocm (RSTCKMAT(tncomp,ndofe,gm));
  //  
  //  Mt->give_code_numbers (eid,cn.a);
  //  eldispl (lcid,eid,r.a,cn.a,ndofe);
  //  
  //  geom_matrix (gm,x,y,0.0,0.0,jac);
  //  mxv (gm,r,eps);
  //  area = eps[2];
  //}
  
  gauss_points (gp.a,w.a,intord);
  
  e2 = u2 = 0;
  
  for (i=0;i<intord;i++){
    xi=gp[i];
    for (j=0;j<intord;j++){
      eta=gp[j];
      
      //if (body){
	eps[0] = Mm->ip[ipp  ].strain[0];
	eps[1] = Mm->ip[ipp  ].strain[1];
	eps[2] = Mm->ip[ippsh].strain[2];
	ipp++;
	jac_2d (jac,x,y,xi,eta);
      //}
      //else{
      //	geom_matrix (gm,x,y,xi,eta,jac);
      //	mxv (gm,r,eps);
      //	eps[2] = area;
      //}
      mxv (d,eps,sig_roof);
      
      bf_lin_4_2d (bf.a,xi,eta);
      give_der_star (bf,rsigfull,nodes,sig_star,Mt->nn);
      
      subv (sig_star,sig_roof,sig_err);
      
      thick = approx (xi,eta,t);
      jac*=thick*w[i]*w[j];
      
      vxmxv (sig_err,dinv,contr);
      e2 += jac*contr;
      
      vxmxv (eps,d,contr);
      u2 += jac*contr;
    }
  }
  
  area = ( (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]) + (x[2]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[2]-y[0]) ) / 2.0;
  sizel = sqrt(area);
  return area;
}

/**
   This function prepare array 'spder' for using in function `least_square::spr_default`.
   It allocates and fills up array 'spder' by derivatives(strain ro stress) in sampling points.
   For planeelemlt sampling points == main integration point == Gauss's point at the centroid of element.
   
   @param eid   - element id
   @param spder - allocated array 
   @param flags - determines by which derivate is spder filled
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemlq :: elchar (long eid, matrix &spsig)
{
  long ipp,ri,ci,lcid;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),eps(ASTCKVEC(tncomp));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));
  
  reallocm(RSTCKMAT(1, tncomp, spsig));
  
  lcid = ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  Mm->matstiff (d,ipp);
  
  Mt->give_node_coord2d (x,y,eid);
  eldispl (lcid,eid,r.a);
  
  geom_matrix (gm,x,y,0.0,0.0,jac);
  mxv (gm,r,eps);
  mxv (d.a, eps.a, spsig.a, d.m, d.n);
  
  vector aux(ASTCKVEC(tncomp));
  Mm->givestress (0,ipp,aux);
  ;
  vector aux2(ASTCKVEC(tncomp));
  Mm->givestrain (0,ipp,aux2);
  ;
}
////////////////////       /* termitovo */       ////////////////////////////////////


/**
   function extracts variables on element
   it is used in stochastic or fuzzy computations
   
   @param at - array of attributes (it describes which variables will be extracted)
   @param val - %vector of extracted values
   
   JK, 17.4.2007
*/
void planeelemlq::extract (atsel &at,vector &/*val*/)
{
  long i;
  
  for (i=0;i<at.num;i++){
    switch (at.atrib[i]){
    case 0:{
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function extract (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
}


/**
   Function computes required mechanical quantity at nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything.
   
   12/06/2012 TKr
*/
void planeelemlq::mechq_nodval (long eid,vector &nodval,nontransquant qt)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  // id of the first integration point on element
  ipid=Mt->elements[eid].ipp[0][0];
  nodip_planelq (ipid,intordsm[0][0],ipnum);

  for (i=0;i<nne;i++){
    //copy nonmechanical quantity from closest int. point
    nodval[i] = Mm->givemechq(qt, ipnum[i]);
  }
}



/**
   Function computes mechanical quantities in nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncnv - number of computed nodes on element (only first ncnv of nodes is calculated)
   @param nq - number of required mechanical quantities
   @param qt - array of types of required mechanical quantities
   
   @return The function does not return anything.
   
   Created by Tomas Koudelka, 29.11.2013
*/
void planeelemlq::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
{
  long i, j, ncompstr, ncompo, ncompnl, ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  intpoints ipb; // backup of the first integration point of element

  // id of the first integration point on element
  ipid=Mt->elements[eid].ipp[0][0];

  // element nodes
  Mt->give_elemnodes(eid, enod);

  //  numbers of integration points closest to element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodip_planelq (ipid,intordsm[0][0],ipnum);

  // compute strains at nodes
  nod_strains(0, eid);
  
  // number of nonloc array components
  ncompnl = Mm->give_num_averq(ipid, Mm->givenonlocid(ipid));

  // store original content of the first integration point on element because 
  // it becomes working int. point for nodal values calculations on the given element
  ipb.copy(Mm->ip[ipid], Mb->nlc, ncompnl, 1);


  // The first integration point will be used for computation of nodal values temporarily
  // then the original content of the first integration point will be restored
  for (i=0;i<ncnv;i++)
  {     
    // number of strain components
    ncompstr = Mm->ip[ipnum[i]].ncompstr;
    // take nodal strain and store them to the first (working) integration point
    Mm->storestrain(0, ipid, 0, ncompstr, Mt->nodes[enod[i]].strain);

    ncompo = Mm->ip[ipnum[i]].ncompeqother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point, the eqother values may be rewritten -> take them from the backup
        Mm->storeeqother(ipid, 0, ncompo, ipb.eqother);
      else
        Mm->storeeqother(ipid, 0, ncompo, Mm->ip[ipnum[i]].eqother);
      // take other values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point the, other values may be rewritten -> take them from the backup
        Mm->storeother(ipid, 0, ncompo, ipb.other);
      else
        Mm->storeother(ipid, 0, ncompo, Mm->ip[ipnum[i]].other);
    }

    if (ncompnl){
      // take nonloc values for the given node from the closest integration point 
      if (ipid == ipnum[i]) // for the first int. point the nonloc values may be rewritten -> take them from the backup
        Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc);
      else
        Mm->storenonloc(ipid, 0, ncompnl, Mm->ip[ipnum[i]].nonloc);
    }

    // compute nodal stresses and internal variables of the material at node
    switch(Mp->matmodel)
    {
      case local:
        Mm->computenlstresses(ipid,Mm->ip[ipid]);
        break;
      case nonlocal:
        Mm->compnonloc_nlstresses (ipid);
        break;
      default:
        print_err("unknown approach in material model is required (Mp->matmodel=%ld)", __FILE__, __LINE__, __func__, Mp->matmodel);
    }
   
    // the internal material variables are stored in ip[ipid].other array and
    // they must be send to ip[ipid].eqother array
    Mm->updateipvalmat(ipid, 0, 0);
    
    //give calculated mechanical quantity from the first int. point
    for (j=0; j<nq; j++)
      nodval[j*ncnv+i] = Mm->givemechq(qt[j], ipid);
  }

  // restore original integration point content of strain/stress/other/eqother/nonloc arrays
  Mm->storestrain(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.strain);
  Mm->storestress(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.stress);
  Mm->storeother(ipid, 0, ipb.ncompother, ipb.other);
  Mm->storeeqother(ipid, 0, ipb.ncompeqother, ipb.eqother);
  Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc); // maybe this would not be necessary
}

