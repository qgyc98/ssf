#include "argyrisplate.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mathem.h"
#include "intp.h"
#include "global.h"
#include "globmat.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"
#include "loadcase.h"
#include <stdlib.h>
#include <math.h>



argyrisplate::argyrisplate (void)
{
  long i,j;
  
  //  number nodes on element
  nne=6;
  //  number of DOFs on element
  ndofe=21;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=1;
  //  order of numerical integration of mass matrix
  intordmm=2;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;
  //  strain/stress state
  ssst=platek;

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
  
  nip[0][0]=7;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=7;
}


argyrisplate::~argyrisplate (void)
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
   function assembles array of shape functions evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::fx (double x,double y,vector &shapef)
{
  shapef[0]  = 1.0;

  shapef[1]  = x;
  shapef[2]  = y;

  shapef[3]  = x*x;
  shapef[4]  = x*y;
  shapef[5]  = y*y;

  shapef[6]  = x*x*x;
  shapef[7]  = x*x*y;
  shapef[8]  = x*y*y;
  shapef[9]  = y*y*y;

  shapef[10] = x*x*x*x;
  shapef[11] = x*x*x*y;
  shapef[12] = x*x*y*y;
  shapef[13] = x*y*y*y;
  shapef[14] = y*y*y*y;

  shapef[15] = x*x*x*x*x;
  shapef[16] = x*x*x*x*y;
  shapef[17] = x*x*x*y*y;
  shapef[18] = x*x*y*y*y;
  shapef[19] = x*y*y*y*y;
  shapef[20] = y*y*y*y*y;
}

/**
   function assembles array of first derivatives of shape functions with respect to x evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::dfdx (double x,double y,vector &shapef)
{
  shapef[0]  = 0.0;

  shapef[1]  = 1.0;
  shapef[2]  = 0.0;

  shapef[3]  = 2.0*x;
  shapef[4]  = y;
  shapef[5]  = 0.0;

  shapef[6]  = 3.0*x*x;
  shapef[7]  = 2.0*x*y;
  shapef[8]  = y*y;
  shapef[9]  = 0.0;

  shapef[10] = 4.0*x*x*x;
  shapef[11] = 3.0*x*x*y;
  shapef[12] = 2.0*x*y*y;
  shapef[13] = y*y*y;
  shapef[14] = 0.0;

  shapef[15] = 5.0*x*x*x*x;
  shapef[16] = 4.0*x*x*x*y;
  shapef[17] = 3.0*x*x*y*y;
  shapef[18] = 2.0*x*y*y*y;
  shapef[19] = y*y*y*y;
  shapef[20] = 0.0;
}

/**
   function assembles array of second derivatives of shape functions with respect to x evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::dfdxdx (double x,double y,vector &shapef)
{
  shapef[0]  = 0.0;

  shapef[1]  = 0.0;
  shapef[2]  = 0.0;

  shapef[3]  = 2.0;
  shapef[4]  = 0.0;
  shapef[5]  = 0.0;

  shapef[6]  = 6.0*x;
  shapef[7]  = 2.0*y;
  shapef[8]  = 0.0;
  shapef[9]  = 0.0;

  shapef[10] = 12.0*x*x;
  shapef[11] = 6.0*x*y;
  shapef[12] = 2.0*y*y;
  shapef[13] = 0.0;
  shapef[14] = 0.0;

  shapef[15] = 20.0*x*x*x;
  shapef[16] = 12.0*x*x*y;
  shapef[17] = 6.0*x*y*y;
  shapef[18] = 2.0*y*y*y;
  shapef[19] = 0.0;
  shapef[20] = 0.0;
}

/**
   function assembles array of first derivatives of shape functions with respect to y evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::dfdy (double x,double y,vector &shapef)
{
  shapef[0]  = 0.0;

  shapef[1]  = 0.0;
  shapef[2]  = 1.0;

  shapef[3]  = 0.0;
  shapef[4]  = x;
  shapef[5]  = 2.0*y;

  shapef[6]  = 0.0;
  shapef[7]  = x*x;
  shapef[8]  = 2.0*x*y;
  shapef[9]  = 3.0*y*y;

  shapef[10] = 0.0;
  shapef[11] = x*x*x;
  shapef[12] = 2.0*x*x*y;
  shapef[13] = 3.0*x*y*y;
  shapef[14] = 4.0*y*y*y;

  shapef[15] = 0.0;
  shapef[16] = x*x*x*x;
  shapef[17] = 2.0*x*x*x*y;
  shapef[18] = 3.0*x*x*y*y;
  shapef[19] = 4.0*x*y*y*y;
  shapef[20] = 5.0*y*y*y*y;
}

/**
   function assembles array of second derivatives of shape functions with respect to y evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::dfdydy (double x,double y,vector &shapef)
{
  shapef[0]  = 0.0;

  shapef[1]  = 0.0;
  shapef[2]  = 0.0;

  shapef[3]  = 0.0;
  shapef[4]  = 0.0;
  shapef[5]  = 2.0;

  shapef[6]  = 0.0;
  shapef[7]  = 0.0;;
  shapef[8]  = 2.0*x;
  shapef[9]  = 6.0*y;

  shapef[10] = 0.0;
  shapef[11] = 0.0;
  shapef[12] = 2.0*x*x;
  shapef[13] = 6.0*x*y;
  shapef[14] = 12.0*y*y;

  shapef[15] = 0.0;
  shapef[16] = 0.0;
  shapef[17] = 2.0*x*x*x;
  shapef[18] = 6.0*x*x*y;
  shapef[19] = 12.0*x*y*y;
  shapef[20] = 20.0*y*y*y;
}

/**
   function assembles array of second derivatives of shape functions with respect to x and y evaluated in given point
   
   @param x, y - coordinates of the given points
   @param shapef - array of shape function values
   
   7.7.2012, JK
*/
void argyrisplate::dfdxdy (double x,double y,vector &shapef)
{
  shapef[0]  = 0.0;

  shapef[1]  = 0.0;
  shapef[2]  = 0.0;

  shapef[3]  = 0.0;
  shapef[4]  = 1.0;
  shapef[5]  = 0.0;

  shapef[6]  = 0.0;
  shapef[7]  = 2.0*x;
  shapef[8]  = 2.0*y;
  shapef[9]  = 0.0;

  shapef[10] = 0.0;
  shapef[11] = 3.0*x*x;
  shapef[12] = 4.0*x*y;
  shapef[13] = 3.0*y*y;
  shapef[14] = 0.0;

  shapef[15] = 0.0;
  shapef[16] = 4.0*x*x*x;
  shapef[17] = 6.0*x*x*y;
  shapef[18] = 6.0*x*y*y;
  shapef[19] = 4.0*y*y*y;
  shapef[20] = 0.0;
}

/**
   function computes coefficients of shape functions
   
   @param eid - element id
   
   7.7.2012, JK
*/
void argyrisplate::shapefunctions (long eid)
{
  long i;
  double nor;
  ivector nodes(nne);
  vector x(nne),y(nne),sf(ndofe),s(2),n(2);
  matrix mat(ndofe,ndofe),unitm(ndofe,ndofe);
  
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  
  
  //  function values in the first node
  fx (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[0][i]=sf[i];
  }
  //  first derivatives with respect to x in the first node
  dfdx (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[1][i]=sf[i];
  }
  //  first derivatives with respect to y in the first node
  dfdy (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[2][i]=sf[i];
  }
  //  second derivatives with respect to x in the first node
  dfdxdx (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[3][i]=sf[i];
  }
  //  second derivatives with respect to x and y in the first node
  dfdxdy (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[4][i]=sf[i];
  }
  //  second derivatives with respect to y in the first node
  dfdydy (x[0],y[0],sf);
  for (i=0;i<ndofe;i++){
    mat[5][i]=sf[i];
  }
  

  //  function values in the second node
  fx (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[6][i]=sf[i];
  }
  //  first derivatives with respect to x in the second node
  dfdx (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[7][i]=sf[i];
  }
  //  first derivatives with respect to y in the second node
  dfdy (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[8][i]=sf[i];
  }
  //  second derivatives with respect to x in the second node
  dfdxdx (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[9][i]=sf[i];
  }
  //  second derivatives with respect to x and y in the second node
  dfdxdy (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[10][i]=sf[i];
  }
  //  second derivatives with respect to y in the second node
  dfdydy (x[1],y[1],sf);
  for (i=0;i<ndofe;i++){
    mat[11][i]=sf[i];
  }
  
  

  //  function values in the third node
  fx (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[12][i]=sf[i];
  }
  //  first derivatives with respect to x in the third node
  dfdx (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[13][i]=sf[i];
  }
  //  first derivatives with respect to y in the third node
  dfdy (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[14][i]=sf[i];
  }
  //  second derivatives with respect to x in the third node
  dfdxdx (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[15][i]=sf[i];
  }
  //  second derivatives with respect to x and y in the third node
  dfdxdy (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[16][i]=sf[i];
  }
  //  second derivatives with respect to y in the third node
  dfdydy (x[2],y[2],sf);
  for (i=0;i<ndofe;i++){
    mat[17][i]=sf[i];
  }


  
  //  first edge
  if (nodes[0]<nodes[1]){
    //  direction vector
    s[0]=x[1]-x[0];
    s[1]=y[1]-y[0];
  }else{
    //  direction vector
    s[0]=x[0]-x[1];
    s[1]=y[0]-y[1];
  }
  //  normal vector
  n[0]=0.0-s[1];
  n[1]=s[0];
  
  nor = n[0]*n[0]+n[1]*n[1];
  n[0]/=sqrt(nor);
  n[1]/=sqrt(nor);
  
  //  first derivatives with respect to x in the first mid-side node
  dfdx (x[3],y[3],sf);
  for (i=0;i<ndofe;i++){
    mat[18][i]=sf[i]*n[0];
  }
  //  first derivatives with respect to y in the first mid-side node
  dfdy (x[3],y[3],sf);
  for (i=0;i<ndofe;i++){
    mat[18][i]+=sf[i]*n[1];
  }
  


  
  //  second edge
  if (nodes[1]<nodes[2]){
    //  direction vector
    s[0]=x[2]-x[1];
    s[1]=y[2]-y[1];
  }else{
    //  direction vector
    s[0]=x[1]-x[2];
    s[1]=y[1]-y[2];
  }
  //  normal vector
  n[0]=0.0-s[1];
  n[1]=s[0];
  
  nor = n[0]*n[0]+n[1]*n[1];
  n[0]/=sqrt(nor);
  n[1]/=sqrt(nor);
  
  //  first derivatives with respect to x in the second mid-side node
  dfdx (x[4],y[4],sf);
  for (i=0;i<ndofe;i++){
    mat[19][i]=sf[i]*n[0];
  }
  //  first derivatives with respect to y in the second mid-side node
  dfdy (x[4],y[4],sf);
  for (i=0;i<ndofe;i++){
    mat[19][i]+=sf[i]*n[1];
  }
  



  
  //  third edge
  if (nodes[2]<nodes[0]){
    //  direction vector
    s[0]=x[0]-x[2];
    s[1]=y[0]-y[2];
  }else{
    //  direction vector
    s[0]=x[2]-x[0];
    s[1]=y[2]-y[0];
  }
  //  normal vector
  n[0]=0.0-s[1];
  n[1]=s[0];
  
  nor = n[0]*n[0]+n[1]*n[1];
  n[0]/=sqrt(nor);
  n[1]/=sqrt(nor);
  
  //  first derivatives with respect to x in the third mid-side node
  dfdx (x[5],y[5],sf);
  for (i=0;i<ndofe;i++){
    mat[20][i]=sf[i]*n[0];
  }
  //  first derivatives with respect to y in the third mid-side node
  dfdy (x[5],y[5],sf);
  for (i=0;i<ndofe;i++){
    mat[20][i]+=sf[i]*n[1];
  }
  
  
  //  the unit matrix
  for (i=0;i<ndofe;i++){
    unitm[i][i]=1.0;
  }
  if (Mt->elements[eid].tmat==NULL){
    Mt->elements[eid].tmat = new matrix [1];
    Mt->elements[eid].tmat[0].m=ndofe;
    Mt->elements[eid].tmat[0].n=ndofe;
    Mt->elements[eid].tmat[0].a=new double [ndofe*ndofe];
  }
  //  solution of system of equations
  gemp (mat.a,Mt->elements[eid].tmat[0].a,unitm.a,ndofe,ndofe,1.0e-15,1);
  
}


/**
   function assembles strain-displacement (geometric) %matrix
   
   @param gm - geometric %matrix
   @param x,y - coordinates of given point
   
   JK, 7.7.2012
*/
void argyrisplate::geom_matrix (matrix &gm,double x,double y,long eid)
{
  long i;
  vector sfdxdx(ndofe),sfdydy(ndofe),sfdxdy(ndofe),v1(ndofe),v2(ndofe),v3(ndofe);
  
  fillm (0.0,gm);
  
  dfdxdx (x,y,sfdxdx);
  dfdydy (x,y,sfdydy);
  dfdxdy (x,y,sfdxdy);

  vxm (sfdxdx,Mt->elements[eid].tmat[0],v1);
  vxm (sfdydy,Mt->elements[eid].tmat[0],v2);
  vxm (sfdxdy,Mt->elements[eid].tmat[0],v3);
  
  for (i=0;i<ndofe;i++){
    gm[0][i]+=v1(i);
    gm[1][i]+=v2(i);
    gm[2][i]+=v3(i);
  }
  
}



/**
   function computes stiffness %matrix of Argyris triangular
   plate finite element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 8.7.2012
*/
void argyrisplate::stiffness_matrix (long eid,long /*ri*/,long /*ci*/,matrix &sm,vector &x,vector &y)
{
  //long i,ii,jj,ipp;
  long i,ipp;
  //double xi,eta,jac,thick;
  double xi,eta,jac;
  ivector nodes(nne);
  vector t(nne),gp1,gp2,w,xx(3),yy(3);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  xx[0]=x[0];
  xx[1]=x[1];
  xx[2]=x[2];
  yy[0]=y[0];
  yy[1]=y[1];
  yy[2]=y[2];
  
  shapefunctions (eid);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  fillm (0.0,sm);
  
  reallocm (ncomp[0],ndofe,gm);
  reallocv (intordsm[0][0],gp1);
  reallocv (intordsm[0][0],gp2);
  reallocv (intordsm[0][0],w);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);
  //  integration point id
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordsm[0][0];i++){
    //  real coordinates of the integration points
    xi = gp1[i]*x[0]+gp2[i]*x[1]+(1.0-gp1[i]-gp2[i])*x[2];
    eta = gp1[i]*y[0]+gp2[i]*y[1]+(1.0-gp1[i]-gp2[i])*y[2];
    
    // geometric matrix
    geom_matrix (gm,xi,eta,eid);
    
    //  stiffness matrix of material
    Mm->matstiff (d,ipp);
    
    //  thickness
    //thick = approx (gp1[i],gp2[i],t);
    
    jac=w[i]*det2d (xx.a,yy.a);
    
    //  contribution to the stiffness matrix of the element
    bdbjac (sm,gm,d,gm,jac);
    
    ipp++;
  }
}

/**
   function assembles
   
   @param eid - element id
   @param sm - the stiffness %matrix
   
   JK, 8.7.2012
*/
void argyrisplate::res_stiffness_matrix (long eid,matrix &sm)
{
  //long transf;
  vector x(nne),y(nne);
  ivector nodes(nne);
  
  //Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  
  stiffness_matrix (eid,0,0,sm,x,y);
  
  //  transformation of stiffness matrix
  //transf = Mt->locsystems (nodes);
  //if (transf>0){
  //matrix tmat (ndofe,ndofe);
  //transf_matrix (nodes,tmat);
  //glmatrixtransf (sm,tmat);
  //}
}
