#include "quadrilatkirch.h"
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


quadrilatkirch::quadrilatkirch (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=16;
  //  number of strain/stress components
  tncomp=3;
  //  number of functions approximated
  napfun=1;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=1;
  //  number of nodes on surface
  nnsurf=4;


  //  order of numerical integration of mass and load matrix
  intordmm=6;
  //intordb=2;
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
  nip[0][0]=16;

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=4;
}



quadrilatkirch::~quadrilatkirch (void)
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
   function approximates nodal values with the help of bi-linear functions
   
   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
   20. 5. 2018, JK
*/
double quadrilatkirch::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function assembles %matrix of polynomials,
   
   @param p - %matrix of polynomials
   @param x,y - coordinates where the %matrix has to be assembled
   
   5. 7. 2018
*/
void quadrilatkirch::pol_matrix (matrix &p,double &x,double &y)
{
  p[0][0]  = x*x*x*y*y*y;
  p[0][1]  = x*x*x*y*y;
  p[0][2]  = x*x*x*y;
  p[0][3]  = x*x*x;
  p[0][4]  = x*x*y*y*y;
  p[0][5]  = x*x*y*y;
  p[0][6]  = x*x*y;
  p[0][7]  = x*x;
  p[0][8]  = x*y*y*y;
  p[0][9]  = x*y*y;
  p[0][10] = x*y;
  p[0][11] = x;
  p[0][12] = y*y*y;
  p[0][13] = y*y;
  p[0][14] = y;
  p[0][15] = 1.0;
}

/**
   function assembles %matrix of polynomials,
   
   @param p - %matrix of polynomials
   @param x,y - coordinates where the %matrix has to be assembled
   
   5. 7. 2018
*/
void quadrilatkirch::bf_matrix (matrix &n,matrix &c,double &x,double &y)
{
  matrix p(ASTCKMAT(1,ndofe));
  
  //  matrix of polynomials
  pol_matrix (p,x,y);
  //  computation of the basis functions
  mxm (p,c,n);
  
}

/**
   function assembles %matrix of linear polynomials,
   
   @param p - %matrix of polynomials
   @param x,y - coordinates where the %matrix has to be assembled
   
   5. 7. 2018
*/
void quadrilatkirch::lbf_matrix (matrix &nn,matrix &cc,double &x,double &y)
{
  matrix pp(ASTCKMAT(1,nne));
  
  //  matrix of polynomials
  pp_matrix (pp,x,y);
  //  computation of the basis functions
  mxm (pp,cc,nn);
  
}

/**
   function assembles %matrix of polynomials, their first derivatives
   and second mixed derivative
   such %matrix is needed for determination of basis functions
   it contains full bi-cubic polynomials
   
   @param p - %matrix of polynomials
   @param x,y - coordinates where the %matrix has to be assembled
   
   19. 5. 2018, JK
*/
void quadrilatkirch::p_matrix (matrix &p,double &x,double &y)
{
  //  functional values
  p[0][0]  = x*x*x*y*y*y;
  p[0][1]  = x*x*x*y*y;
  p[0][2]  = x*x*x*y;
  p[0][3]  = x*x*x;
  p[0][4]  = x*x*y*y*y;
  p[0][5]  = x*x*y*y;
  p[0][6]  = x*x*y;
  p[0][7]  = x*x;
  p[0][8]  = x*y*y*y;
  p[0][9]  = x*y*y;
  p[0][10] = x*y;
  p[0][11] = x;
  p[0][12] = y*y*y;
  p[0][13] = y*y;
  p[0][14] = y;
  p[0][15] = 1.0;
  
  //  derivatve with respect to y which is rotation around x
  p[1][0]  = 3.0*x*x*x*y*y;
  p[1][1]  = 2.0*x*x*x*y;
  p[1][2]  = x*x*x;
  p[1][3]  = 0.0;
  p[1][4]  = 3.0*x*x*y*y;
  p[1][5]  = 2.0*x*x*y;
  p[1][6]  = x*x;
  p[1][7]  = 0.0;
  p[1][8]  = 3.0*x*y*y;
  p[1][9]  = 2.0*x*y;
  p[1][10] = x;
  p[1][11] = 0.0;
  p[1][12] = 3.0*y*y;
  p[1][13] = 2.0*y;
  p[1][14] = 1.0;
  p[1][15] = 0.0;
  
  //  negative derivative with respect to x which is rotation around y
  p[2][0]  = -3.0*x*x*y*y*y;
  p[2][1]  = -3.0*x*x*y*y;
  p[2][2]  = -3.0*x*x*y;
  p[2][3]  = -3.0*x*x;
  p[2][4]  = -2.0*x*y*y*y;
  p[2][5]  = -2.0*x*y*y;
  p[2][6]  = -2.0*x*y;
  p[2][7]  = -2.0*x;
  p[2][8]  = -1.0*y*y*y;
  p[2][9]  = -1.0*y*y;
  p[2][10] = -1.0*y;
  p[2][11] = -1.0;
  p[2][12] = 0.0;
  p[2][13] = 0.0;
  p[2][14] = 0.0;
  p[2][15] = 0.0;

  //  mixed second derivative
  p[3][0]  = 9.0*x*x*y*y;
  p[3][1]  = 6.0*x*x*y;
  p[3][2]  = 3.0*x*x;
  p[3][3]  = 0.0;
  p[3][4]  = 6.0*x*y*y;
  p[3][5]  = 4.0*x*y;
  p[3][6]  = 2.0*x;
  p[3][7]  = 0.0;
  p[3][8]  = 3.0*y*y;
  p[3][9]  = 2.0*y;
  p[3][10] = 1.0;
  p[3][11] = 0.0;
  p[3][12] = 0.0;
  p[3][13] = 0.0;
  p[3][14] = 0.0;
  p[3][15] = 0.0;
}


/**
   function assembles the S %matrix
   it is used for determination of basis functions
   
   @param s - the %matrix
   @param x,y - nodal coordinates
   
   19. 5. 2018, JK
*/
void quadrilatkirch::s_matrix (matrix &s,vector &x,vector &y)
{
  long i,j;
  double xx,yy;
  matrix p(4,16);
  
  xx = x[0];  yy = y[0];
  p_matrix (p,xx,yy);
  for (i=0;i<4;i++){
    for (j=0;j<16;j++){
      s[i][j]=p[i][j];
    }
  }

  xx = x[1];  yy = y[1];
  p_matrix (p,xx,yy);
  for (i=0;i<4;i++){
    for (j=0;j<16;j++){
      s[i+4][j]=p[i][j];
    }
  }

  xx = x[2];  yy = y[2];
  p_matrix (p,xx,yy);
  for (i=0;i<4;i++){
    for (j=0;j<16;j++){
      s[i+8][j]=p[i][j];
    }
  }

  xx = x[3];  yy = y[3];
  p_matrix (p,xx,yy);
  for (i=0;i<4;i++){
    for (j=0;j<16;j++){
      s[i+12][j]=p[i][j];
    }
  }

}

/**
   function computes coefficients of basis functions
   the C %matrix is the inverse of the S %matrix
   
   @param c - %matrix of coefficients
   @param x,y - nodal coordinates
   
   19. 5. 2018, JK
*/
void quadrilatkirch::c_matrix (matrix &c,vector &x, vector &y)
{
  long i;
  matrix s(ndofe,ndofe),id(ndofe,ndofe);
  
  s_matrix (s,x,y);
  nullm (id);
  for (i=0;i<ndofe;i++){
    id[i][i]=1.0;
  }
  //id[7][7]=1.0;
  
  gemp (s.a,c.a,id.a,ndofe,ndofe,1.0e-12,1);

  s_matrix (s,x,y);
  
  mxm (s,c,id);

  /*
  FILE *ma;
  ma=fopen ("coef.txt","w");
  for (i=0;i<16;i++){
    fprintf (ma,"\n");
    for (long j=0;j<16;j++){
      fprintf (ma,"%2ld %2ld   % 16.12le\n",i+1,j+1,id[i][j]);
      //fprintf (ma,"%2ld %2ld   % 16.12le\n",i+1,j+1,c[i][j]);
    }
  }
  fclose (ma);
  */

}

/**
   function computes coefficients of linear basis functions
   
   @param cc - %matrix of coefficients
   @param x,y - nodal coordinates
   
   19. 5. 2018, JK
*/
void quadrilatkirch::cc_matrix (matrix &cc,vector &x, vector &y)
{
  long i;
  matrix s(nne,nne),id(nne,nne);
  
  s[0][0] = x[0]*y[0];
  s[0][1] = x[0];
  s[0][2] = y[0];
  s[0][3] = 1.0;

  s[1][0] = x[1]*y[1];
  s[1][1] = x[1];
  s[1][2] = y[1];
  s[1][3] = 1.0;

  s[2][0] = x[2]*y[2];
  s[2][1] = x[2];
  s[2][2] = y[2];
  s[2][3] = 1.0;

  s[3][0] = x[3]*y[3];
  s[3][1] = x[3];
  s[3][2] = y[3];
  s[3][3] = 1.0;
  
  nullm (id);
  for (i=0;i<nne;i++){
    id[i][i]=1.0;
  }
  
  gemp (s.a,cc.a,id.a,nne,nne,1.0e-16,1);


}

/**
   function assembles %matrix of bi-linear polynomials
   
   @param cc - %matrix of coefficients
   @param x,y - nodal coordinates
   
   20. 5. 2018, JK
*/
void quadrilatkirch::pp_matrix (matrix &pp,double &x, double &y)
{
  pp[0][0] = x*y;
  pp[0][1] = x;
  pp[0][2] = y;
  pp[0][3] = 1.0;
}

/**
   function assembles %matrix of second partial derivatives of the deflection function
   it is used for computation of curvatures
   
   @param g - %matrix of polynomials
   @param x,y - coordinates where the %matrix has to be assembled
   
   19. 5. 2018, JK
*/
void quadrilatkirch::g_matrix (matrix &g,double &x,double &y)
{
  g[0][0]  = -6.0*x*y*y*y;
  g[0][1]  = -6.0*x*y*y;
  g[0][2]  = -6.0*x*y;
  g[0][3]  = -6.0*x;
  g[0][4]  = -2.0*y*y*y;
  g[0][5]  = -2.0*y*y;
  g[0][6]  = -2.0*y;
  g[0][7]  = -2.0;
  g[0][8]  = 0.0;
  g[0][9]  = 0.0;
  g[0][10] = 0.0;
  g[0][11] = 0.0;
  g[0][12] = 0.0;
  g[0][13] = 0.0;
  g[0][14] = 0.0;
  g[0][15] = 0.0;

  g[1][0]  = -6.0*x*x*x*y;
  g[1][1]  = -2.0*x*x*x;
  g[1][2]  = 0.0;
  g[1][3]  = 0.0;
  g[1][4]  = -6.0*x*x*y;
  g[1][5]  = -2.0*x*x;
  g[1][6]  = 0.0;
  g[1][7]  = 0.0;
  g[1][8]  = -6.0*x*y;
  g[1][9]  = -2.0*x;
  g[1][10] = 0.0;
  g[1][11] = 0.0;
  g[1][12] = -6.0*y;
  g[1][13] = -2.0;
  g[1][14] = 0.0;
  g[1][15] = 0.0;

  g[2][0]  = -18.0*x*x*y*y;
  g[2][1]  = -12.0*x*x*y;
  g[2][2]  = -6.0*x*x;
  g[2][3]  = 0.0;
  g[2][4]  = -12.0*x*y*y;
  g[2][5]  = -8.0*x*y;
  g[2][6]  = -4.0*x;
  g[2][7]  = 0.0;
  g[2][8]  = -6.0*y*y;
  g[2][9]  = -4.0*y;
  g[2][10] = -2.0;
  g[2][11] = 0.0;
  g[2][12] = 0.0;
  g[2][13] = 0.0;
  g[2][14] = 0.0;
  g[2][15] = 0.0;

}

/**
   function assembles geometric %matrix
   i.e. matrix between curvatures and nodal DOFs
   
   @param gm - the geometric %matrix
   @param c - %matrix of polynomial coefficients
   @param x,y - coordinates where the %matrix has to be assembled
   
   19. 5. 2018, JK
*/
void quadrilatkirch::geom_matrix (matrix &gm,matrix &c,double &x,double &y)
{
  matrix g(tncomp,ndofe);
  
  g_matrix (g,x,y);
  
  mxm (g,c,gm);
}

/**
   function assembles matrix of first partial derivatives of deflection function
   
   @param h - %matrix of first derivatives
   @param x,y - coordinates where the %matrix has to be assembled
   
   19. 5. 2018, JK
*/
void quadrilatkirch::h_matrix (matrix &h,double &x,double &y)
{
  h[0][0]  = 3.0*x*x*y*y*y;
  h[0][1]  = 3.0*x*x*y*y;
  h[0][2]  = 3.0*x*x*y;
  h[0][3]  = 3.0*x*x;
  h[0][4]  = 2.0*x*y*y*y;
  h[0][5]  = 2.0*x*y*y;
  h[0][6]  = 2.0*x*y;
  h[0][7]  = 2.0*x;
  h[0][8]  = y*y*y;
  h[0][9]  = y*y;
  h[0][10] = y;
  h[0][11] = 1.0;
  h[0][12] = 0.0;
  h[0][13] = 0.0;
  h[0][14] = 0.0;
  h[0][15] = 0.0;

  h[1][0]  = 3.0*x*x*x*y*y;
  h[1][1]  = 2.0*x*x*x*y;
  h[1][2]  = x*x*x;
  h[1][3]  = 0.0;
  h[1][4]  = 3.0*x*x*y*y;
  h[1][5]  = 2.0*x*x*y;
  h[1][6]  = x*x;
  h[1][7]  = 0.0;
  h[1][8]  = 3.0*x*y*y;
  h[1][9]  = 2.0*x*y;
  h[1][10] = x;
  h[1][11] = 0.0;
  h[1][12] = 3.0*y*y;
  h[1][13] = 2.0*y;
  h[1][14] = 1.0;
  h[1][15] = 0.0;
}

/**
   function assembles geometric %matrix for evaluation of the initial stress %matrix
   
   @param gm - geometric %matrix
   @param c - %matrix of polynomial coefficients
   @param x,y - coordinates where the %matrix has to be assembled
   
   19. 5. 2018, JK
*/
void quadrilatkirch::geom_initstiff_matrix (matrix &gm,matrix &c,double &x,double &y)
{
  matrix h(2,ndofe);
  
  h_matrix (h,x,y);
  
  mxm(h,c,gm);
}


/**
   function assembles the stiffness %matrix
   
*/
void quadrilatkirch::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,xx,yy,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix gm,d(ASTCKMAT(tncomp,tncomp)),c(ASTCKMAT(ndofe,ndofe)),cc(ASTCKMAT(nne,nne));

  //Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  fillm (0.0,sm);
  
  //  matrix of bi-cubic polynomial coefficients
  //  the coefficients are computed on the original mesh
  //  there is no transformation to unit element
  c_matrix (c,x,y);

  //  matrix of bi-linear polynomial coefficients
  //  the coefficients are computed on the original mesh
  //  there is no transformation to unit element
  cc_matrix (cc,x,y);

  
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  
  reallocm (RSTCKMAT(ncomp[0],ndofe,gm));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];

      jac_2d (jac,x,y,xi,eta);

      //  thickness in integration point
      thick = approx (xi,eta,t);

      xx = approx (xi,eta,x);
      yy = approx (xi,eta,y);
      
      
      //  geometric matrix
      geom_matrix (gm,c,xx,yy);
      
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      /*
      fprintf (Out,"\n\n ipp %3ld",ipp);
      fprintf (Out,"\n xx  %15.12lf",xx);
      fprintf (Out,"\n yy  %15.12lf",yy);
      fprintf (Out,"\n jac %15.12lf",jac);
      fprintf (Out,"\n t   %15.12lf",thick);
      */
      
      jac*=thick*thick*thick*w[i]*w[j];
      
      jac=fabs(jac);
      
      //  contribution to the stiffness matrix of the element
      bdbjac (sm,gm,d,gm,jac);
      
      ipp++;
    }
  }
}








/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   base vectors e1 and e2 have to contain 2 components
   the components are in the midplane
   
   @param nodes - element nodes
   @param tmat - transformation %matrix
   
   JK,
*/
void quadrilatkirch::transfmat (ivector &nodes,matrix &tmat)
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
   function assembles stiffness %matrix of Kirchhoff rectangular
   finite element with bi-cubic approximation functions
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   20. 5. 2018, JK
*/
void quadrilatkirch::res_stiffness_matrix (long eid,matrix &sm)
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
    transfmat (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}














/**
   function computes initial stress stiffness %matrix of Q4 element
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ism - i-s matrix
   
   15.11.2003
*/
void quadrilatkirch::initstr_matrix (long /*eid*/,long /*ri*/,long /*ci*/,matrix &/*ism*/)
{
  /*
  long i,i1,i2,j,j1,j2,k,ii,jj;
  double jac,a,a0,a1,thick,xx,yy;
  ivector nodes(nne);
  vector w,gp,l(2),x(nne),y(nne),z(nne),sig(3),sx(3),sy(3),t(nne),dwx(12),dwy(12);
  matrix tran(3,3),atd(16,12);

// sig....... Axial stress

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

// translation to center it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  
  //     jakobian
  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);

  atd_matrix (atd,x,y,eid);

  fillm (0.0,ism);

//  ii=Mt->elements[eid].ipp[sip][0];
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (intordsm[ii][jj],w); reallocv (intordsm[ii][jj],gp);

	  for (i=0;i<intordsm[ii][jj];i++){
		 l[0]=gp[i];  
		for (j=0;j<intordsm[ii][jj];j++){
		  l[1]=gp[j];
  	      xx=(x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
	      yy=(y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;
          for (k=0;k<dwx.n;k++){
           dwx[i]=atd[1][k]+yy*atd[3][k]+2.*xx*atd[4][k]+2.*xx*yy*atd[6][k]+
                  yy*yy*atd[7][k]+3.*xx*xx*atd[8][k]+3.*xx*xx*yy*atd[10][k]+yy*yy*yy*atd[11][k];
           dwy[i]=atd[2][k]+xx*atd[3][k]+2.*yy*atd[5][k]+xx*xx*atd[6][k]+2.*yy*xx*     
                  atd[7][k]+3.*yy*yy*atd[9][k]+xx*xx*xx*atd[10][k]+3.*xx*yy*yy*atd[11][k];
		  }  
  
	      thick=approx (l[0],l[1],t);
		  jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j]*thick/24.;
          for (i2=0;i2<sig.n;i2++){
            i1=3*i2-3;
            for (j2=0;j2<nne;j2++){
             j1=3*j2-3;
             ism[i1][j1]=ism[i1][j1]+( dwx[i]*dwx[j]*sig[0]
                +(dwx[i]*dwy[j]+dwx[j]*dwy[i])*sig[2]+dwy[i]*dwy[j]*sig[1] )*jac;
			}
		  }
		}
	  }
	}
  }
       

  */	
	
}



/**
   function computes load %matrix of the plate rectangular finite element
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 5. 7. 2018
*/
void quadrilatkirch::load_matrix (long /*eid*/,matrix &lm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,xx,yy;
  vector w,gp;
  matrix n(ASTCKMAT(napfun,ndofe)),nn(ASTCKMAT(napfun,nne)),c(ASTCKMAT(ndofe,ndofe)),cc(ASTCKMAT(nne,nne));
  matrix alm(ASTCKMAT(ndofe,nne));
  
  fillm (0.0,lm);

  //  matrix of bi-cubic polynomial coefficients
  //  the coefficients are computed on the original mesh
  //  there is no transformation to unit element
  c_matrix (c,x,y);

  //  matrix of bi-linear polynomial coefficients
  //  the coefficients are computed on the original mesh
  //  there is no transformation to unit element
  cc_matrix (cc,x,y);


  reallocv (RSTCKVEC(intordmm,w));
  reallocv (RSTCKVEC(intordmm,gp));
  
  gauss_points (gp.a,w.a,intordmm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      
      jac_2d (jac,x,y,xi,eta);
      
      xx = approx (xi,eta,x);
      yy = approx (xi,eta,y);
      
      bf_matrix (n,c,xx,yy);
      lbf_matrix (nn,cc,xx,yy);
      
      jac*=w[i]*w[j];
      jac=fabs(jac);
      
      //  contribution to the load matrix of the element
      mtxm (n,nn,alm);
      cmulm (jac,alm,alm);
      addm (alm,lm,lm);
      
    }
  }
}

/**
   @param nv - array of nodal values
   @param nf - array of nodal forces and moments
   
   JK, 5. 7. 2018
*/
void quadrilatkirch::surfload (long /*lcid*/,long eid,double *nv,vector &nf)
{
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  matrix lm(ASTCKMAT(ndofe,nne));
  
  //  node cordinates
  Mt->give_node_coord2d (x,y,eid);

  //  load matrix
  load_matrix (eid,lm,x,y);
  
  mxv (lm.a,nv,nf.a,ndofe,nne);
}
