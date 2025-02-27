#include "barelq2d.h"
#include "barel2d.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "loadcase.h"
#include <math.h>



barelq2d::barelq2d (void)
{
  long i;

  nne=3;  ndofe=6;  tncomp=1;  napfun=2;
  tnip=3;
  ssst = bar;

  nb=1;

  ncomp = new long [nb];
  ncomp[0]=1;

  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }

  nip[0][0]=3;

  intordsm[0][0]=3;
  intordmm=3;
  
  zero=Mp->zero;
}

barelq2d::~barelq2d (void)
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
   Function computes direction %vector.
   
   @param s - direction %vector
   @param x,y - vectors of nodal coordinates at global coordinate system

   @retval Function returns direction %vector s.
   
   TKR podle JK, 31.10.2006
*/
void barelq2d::dirvect (vector &s,vector &x,vector &y)
{
  double d;
  
  s[0]=x[1]-x[0];
  s[1]=y[1]-y[0];
  
  d=sqrt(s[0]*s[0]+s[1]*s[1]);
  
  if (d<zero){
    print_err("zero norm of direction vector in function dirvect",__FILE__,__LINE__,__func__);
  }
  
  s[0]/=d;
  s[1]/=d;
}

/**
  Function computes local coordinates.
  Coordinates of nodes after transformation of the bar to the x axis.

  @param x  - %vector of nodal x-coordinates in global coordinate system
  @param y  - %vector of nodal y-coordinates in global coordinate system
  @param lx - %vector of nodal coordinates at local element coordinate system

  @retval Function returns %vector of nodal coordinates at parameter lx.
*/
void barelq2d::giveloccoord(vector &x, vector &y,vector &lx)
{
  double l;

  lx[0] = 0.0;
  l = sqr(x[1] - x[0]) + sqr(y[1] - y[0]);
  l = sqrt(l);
  lx[1] = l;
  l = sqr(x[2] - x[0]) + sqr(y[2] - y[0]);
  l = sqrt(l);
  lx[2] = l;
}



/**
   Function approximates function defined by nodal values.

   @param xi,eta - coordinates on element
   @param nodval - nodal values

   @retval Function returns approximated value for the given natural coordinate xi.
*/
double barelq2d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));

  bf_quad_1d (bf.a,xi);

  scprd (bf,nodval,f);

  return f;
}


/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param s - direction %vector
   @param xi - natural coordinate
   
   JK, 25.8.2001
*/
void barelq2d::bf_matrix (matrix &n,vector &s,double xi)
{
  vector bf(ASTCKVEC(nne));
  
  bf_quad_1d (bf.a,xi);
  
  n[0][0]=bf[0]*s[0];
  n[0][1]=bf[0]*s[1];
  n[0][2]=bf[1]*s[0];
  n[0][3]=bf[1]*s[1];
  n[0][4]=bf[2]*s[0];
  n[0][5]=bf[2]*s[1];

}


/**
   Function assembles strain-displacement (geometric) %matrix.

   @param gm - geometric %matrix
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system
   @param xi - natural coordinate
   @param jac - Jacobian
   
   @retval Function returns Jacobian in parameter jac and geometric 
           %matrix in parameter gm. 

   JK, 10.8.2001
*/
void barelq2d::geom_matrix (matrix &gm, vector &gx,vector &gy,double xi,double &jac)
{
  vector s(ASTCKVEC(2)),lx(ASTCKVEC(nne)),dx(ASTCKVEC(nne));
  
  //  direction vector
  dirvect (s,gx,gy);
  //  local coordinates of nodes
  giveloccoord (gx,gy,lx);
  
  dx_bf_quad_1d (dx.a,xi);
  derivatives_1d (dx, jac, lx, xi);
  
  fillm(0.0, gm);
  
  gm[0][0]= dx[0]*s[0];
  gm[0][1]= dx[0]*s[1];
  gm[0][2]= dx[1]*s[0];
  gm[0][3]= dx[1]*s[1];
  gm[0][4]= dx[2]*s[0];
  gm[0][5]= dx[2]*s[1];
}





/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   x_local = T x_global
   
   @param x - %vector of nodal x-coordinates in the global coordinate system
   @param y - %vector of nodal y-coordinates in the global coordinate system
   @param tran - tranformation %matrix T
   
   @retval Function returns transformation %matrix in parameter tmat.

   JK, 12.10.2008
*/
/*
void barelq2d::give_glob_loc_tmat(vector &x, vector &y, matrix &tmat)
{
  double c, s, l;

  fillm(0.0, tmat);
  l = sqrt(sqr(x[1]-x[0]) + sqr(y[1]-y[0]));
  c = (x[1] - x[0]) / l;
  s = (y[1] - y[0]) / l;
  tmat[0][0] = tmat[2][2] = tmat[4][4] =  c;
  tmat[0][1] = tmat[2][3] = tmat[4][5] =  s;
  tmat[1][0] = tmat[3][2] = tmat[5][4] = -s;
  tmat[1][1] = tmat[3][3] = tmat[5][5] =  c;
}
*/


/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   T x_local = x_global
   
   @param x - %vector of nodal x-coordinates in the global coordinate system
   @param y - %vector of nodal y-coordinates in the global coordinate system
   @param tran - tranformation %matrix T
   
   @retval Function returns transformation %matrix in parameter tmat.

   JK, 12.10.2008
*/
 /*
void barelq2d::give_loc_glob_tmat(vector &x, vector &y, matrix &tmat)
{
  double c, s, l;

  fillm(0.0, tmat);
  l = sqrt(sqr(x[1]-x[0]) + sqr(y[1]-y[0]));
  c = (x[1] - x[0]) / l;
  s = (y[1] - y[0]) / l;
  tmat[0][0] = tmat[2][2] = tmat[4][4] =  c;
  tmat[0][1] = tmat[2][3] = tmat[4][5] = -s;
  tmat[1][0] = tmat[3][2] = tmat[5][4] =  s;
  tmat[1][1] = tmat[3][3] = tmat[5][5] =  c;
}
 */


/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   T x_local = x_global
   x_local = T^T x_global
   
   @param tran - tranformation %matrix T
   @param gx,gy - vectors of nodal coordinates in global coordinate systems
   
   JK, 12.10.2008
*/
void barelq2d::tran_mat (const vector &gx, const vector &gy, matrix &tran)
{
  double l;
  
  fillm(0.0, tran);
  //  first direction vector is located in the first column
  tran[0][0]=gx[1]-gx[0];
  tran[1][0]=gy[1]-gy[0];

  //  length of the first direction vector
  l=sqrt(tran[0][0]*tran[0][0]+tran[1][0]*tran[1][0]);
  
  //  normalization of the first direction vector
  tran[0][0]/=l;
  tran[1][0]/=l;

  //  second direction vector is located in the second column
  tran[0][1]=0.0-tran[1][0];
  tran[1][1]=tran[0][0]; 
}



/**
   Function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l.
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix

   @retval Function returns transformation %matrix in parameter tmat.
*/
void barelq2d::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i;
  fillm (0.0,tmat);

  for (i=0;i<tmat.m;i++){
    tmat[i][i]=1.0;
  }

  for (i=0;i<nodes.n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*2][i*2]=Mt->nodes[nodes[i]].e1[0];    tmat[i*2][i*2+1]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2]=Mt->nodes[nodes[i]].e1[1];  tmat[i*2+1][i*2+1]=Mt->nodes[nodes[i]].e2[1];
    }
  }

}



/**
   Function computes stiffness matrix of one element.

   @param eid - element id
   @param ri - row index
   @param sm - stiffness matrix
   @param gx,gy - node coordinates in global system
   
   @retval Function returns stiffness %matrix at parameter sm.

   10.8.2001
*/
void barelq2d::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &gx,vector &gy)
{
  long i,ipp;
  double a,xi,jac;
  vector w,gp;
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  fillm (0.0,sm);


  Mc->give_area (eid,a);

  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  reallocm (RSTCKMAT(ncomp[0],ndofe,gm));
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    
    //  stran-displacement matrix
    geom_matrix (gm,gx,gy,xi,jac);
    //  matrix of stiffness of the material
    Mm->matstiff (d,ipp);
    jac*=a*w[i];
    //  contribution to the stiffness matrix of the element
    bdbjac (sm,gm,d,gm,jac);

    ipp++;
  }
  
}


/**
   Function computes stiffness %matrix of one element. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param eid - number of element
   @param sm  - stiffness %matrix

   @retval Function returns resulting stiffness %matrix in parameter sm
*/
void barelq2d::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  //  node coordinates in global coordinate system
  Mt->give_node_coord2d (x,y,eid);
  //  node id
  Mt->give_elemnodes (eid,nodes);
  
  //  stiffness matrix
  stiffness_matrix (eid,0,0,sm,x,y);
  

  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}



/**
   Function computes mass %matrix of twodimensional bar element.
   
   @param eid - element id
   @param mm - mass %matrix
   @param x,y - node coordinates in global system
   
   @retval Function returns mass %matrix at parameter mm.

   JK, 10.8.2001
*/
void barelq2d::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double a,rho,jac,xi;
  ivector nodes(ASTCKIVEC(nne));
  vector lx(ASTCKVEC(nne)),s(ASTCKVEC(2));
  vector dens(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm));
  matrix n(ASTCKMAT(1,ndofe));
  
  fillm (0.0,mm);
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  
  //  local coordinates of nodes
  //  only x coordinates are nonzero
  giveloccoord(x,y,lx);
  //  direction vector
  dirvect (s,x,y);
  
  //  area of bar cross-section
  Mc->give_area (eid,a);
  //  density of the material
  Mc->give_density (eid,nodes,dens);
  
  //  Gauss points
  gauss_points (gp.a,w.a,intordmm);
  
  
  for (i=0;i<intordmm;i++){
    xi = gp[i];
    
    //  jacobian
    jac_1d (jac,lx,xi);
 
    //  matrix of approximation functions
    bf_matrix (n,s,xi);
    //  density in integration point
    rho = approx (xi,dens);
    
    jac*=rho*a*w[i];
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
}



/**
   Function computes mass %matrix of one element. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param eid - number of element
   @param mm  - mass %matrix

   @retval Function returns resulting mass %matrix in parameter mm
*/
void barelq2d::res_mass_matrix (long eid,matrix &mm)
{
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  long transf;
  
  //  node coordinates in global system
  Mt->give_node_coord2d (x,y,eid);
  
  //  mas matrix
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
   Function computes strain at main integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barelq2d::res_mainip_strains (long lcid,long eid)
{
  vector aux,gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;

  Mt->give_node_coord2d (gx,gy,eid);
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
  
  mainip_strains (lcid,eid,0,0,gx,gy,r);
}



/**
   Function computes strains at main integration points of element.

   @param lcid - load case id
   @param eid  - element id
   @param ri - row index
   @param ci - column index
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system
   @param r - %vector of nodal displacements

   10.5.2002
*/
void barelq2d::mainip_strains (long lcid,long eid,long ri,long ci,vector &gx,vector &gy,vector &r)
{
  long i,ii,ipp;
  double xi,jac;
  vector gp,w,eps;
  matrix gm;

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(ncomp[ii],eps));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));

    gauss_points (gp.a,w.a,intordsm[ii][ii]);

    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];

      geom_matrix (gm,gx,gy,xi,jac);
      mxv (gm,r,eps);

      Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
      ipp++;
    }
  }

}



/**
   Function computes strains in nodes of element.

   @param lcid - load case id
   @param eid  - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void barelq2d::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i, j, ii, jj, ipp, lipp;
  double l;
  ivector enod(ASTCKIVEC(nne)), ipnum(ASTCKIVEC(nne));
  vector gx(ASTCKVEC(nne)), gy(ASTCKVEC(nne)),nxi(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe)),eps(ASTCKVEC(tncomp)),aux;
  vector w(ASTCKVEC(tnip)),auxw(ASTCKVEC(tnip)),gp(ASTCKVEC(tnip));
  matrix tmat,gm(ASTCKMAT(tncomp,ndofe));
  
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodipnum (eid,0,0,ipnum);

  //  node numbers
  Mt->give_elemnodes (eid, enod);

  if (Mp->strainaver==2)
  {
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
        if (intordsm[ii][jj]==0)
          continue;
        gauss_points (gp.a, auxw.a, intordsm[ii][jj]);
        lipp = 0;
        for (i=0;i<intordsm[ii][jj];i++){
          w[lipp] = 0.5*auxw[i];
          lipp++;
        }
      }
    }
  }

  //  loop over nodes
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain (lcid,ipnum[i],eps);
    
    j=enod[i];
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      l = Mt->give_length(eid)*w[ipnum[i]-ipp];
      cmulv (l,eps);
      Mt->nodes[j].storestrain(lcid,0,l,eps);
    }
  }
}



/**
  The function computes nodal strains directly, averageing of nodal strain values is performed according to setup.
   
  @param lcid - load case id
  @param eid - element id

  @return The function does not return anything but it changes arrays of strains at 
          corresponding nodes of the element.

  Created by Tomas Koudelka, 14.11.2013
*/
void barelq2d::nod_strains_comp (long lcid,long eid)
{
  long i, j, ii, jj, ipp, lipp;
  double l, jac;
  ivector enod(ASTCKIVEC(nne)), ipnum(ASTCKIVEC(nne));
  vector gx(ASTCKVEC(nne)), gy(ASTCKVEC(nne)), nxi(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe)), eps(ASTCKVEC(tncomp)), aux;
  vector w(ASTCKVEC(tnip)),auxw(ASTCKVEC(tnip)),gp(ASTCKVEC(tnip));
  matrix tmat,gm(ASTCKMAT(tncomp,ndofe));
  
  Mt->give_node_coord2d (gx,gy,eid);
  //giveloccoord (gx,gy,x);
  //  direction vector
  //dirvect (s,gx,gy);
  Mt->give_elemnodes (eid,enod);
  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  if (Mp->strainaver==2)
  {
    ipp=Mt->elements[eid].ipp[0][0];
    nodipnum (eid,0,0,ipnum);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
        if (intordsm[ii][jj]==0)
          continue;
        gauss_points (gp.a, auxw.a, intordsm[ii][jj]);
        lipp = 0;
        for (i=0;i<intordsm[ii][jj];i++){
          w[lipp] = 0.5*auxw[i];
          lipp++;
        }
      }
    }
  }

  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_barq(nxi);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,gx,gy,nxi[i],jac);

    //  strain computation
    mxv (gm,r,eps);
    
    //  storage of strains to the node
    j=enod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      l = 2.0*jac*w[ipnum[i]-ipp];
      cmulv (l,eps);
      Mt->nodes[j].storestrain(lcid,0,l,eps);
    }
  }
}



/**
  The function computes nodal strains directly.
   
  @param lcid - load case id
  @param eid - element id

  @return The function does not return anything but it changes arrays of strains at 
          corresponding nodes of the element.

  Created by Tomas Koudelka, 23.5.2018
*/
void barelq2d::nod_strains (long lcid,long eid)
{
  long i,j,k,m;
  double jac;
  ivector enod(ASTCKIVEC(nne));
  vector gx(ASTCKVEC(nne)), gy(ASTCKVEC(nne)),nxi(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe)),eps(ASTCKVEC(tncomp)),aux;
  matrix tmat,gm(ASTCKMAT(tncomp,ndofe));
  
  Mt->give_node_coord2d (gx,gy,eid);
  Mt->give_elemnodes (eid,enod);
  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_barq(nxi);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,gx,gy,nxi[i],jac);

    //  strain computation
    mxv (gm,r,eps);
    
    //  storage of nodal strains
    j=enod[i];
    m = lcid*tncomp;
    for (k=0; k<tncomp; k++)
      Mt->nodes[j].strain[m+k] = eps(k);
  }
}



/**
   Function computes all strain components at all integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2004
*/
void barelq2d::res_allip_strains (long lcid,long eid)
{
  //  all strain components at all integration points
  allip_strains (lcid,eid,0,0);
}



/**
   Function assembles all values at all integration points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 25.9.2004
*/
void barelq2d::allip_strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
}



/**
   Function computes strains at strain points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void barelq2d::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
//  double **stra;
  vector coord,eps;

  switch (Mm->stra.tape[eid]){
    case nowhere:{
      break;
    }
    case intpts:{
      allip_strains (lcid,eid,ri,ci);
      break;
    }
    case enodes:{
      nod_strains_ip (lcid,eid,ri,ci);
      break;
    }
    case userdefined:{
      //  number of auxiliary element points
      naep = Mm->stra.give_naep (eid);
      ncp = Mm->stra.give_ncomp (eid);
      sid = Mm->stra.give_sid (eid);
      reallocv (RSTCKVEC(ncp,eps));
      reallocv (RSTCKVEC(1,coord));
      for (i=0;i<naep;i++){
        Mm->stra.give_aepcoord (sid,i,coord);

        if (Mp->strainaver==0)
          //appval (coord[0],0,ncp,eps,stra);
        if (Mp->strainaver==1)
	  //appstrain (lcid,eid,coord[0],0,ncp,eps);

        Mm->stra.storevalues(lcid,eid,i,eps);
      }
      /*
      if (Mp->strainaver==0)
      {
        for (i=0;i<nne;i++)
          delete [] stra[i];

        delete [] stra;
      }
      */
      break;
    }
    default:{
      print_err("unknown strain point is required on element (eid=%ld)",__FILE__,__LINE__, __func__, eid+1);
    }
  }
}



/**
   Function assembles natural coordinates of nodes of element.

   @param xi - array containing natural coordinates xi

   @retval Function returns natural coordinates in %vector xi.

   10.5.2002
*/
void barelq2d::nodecoord (vector &xi)
{
  xi[0] = -1.0;
  xi[1] =  1.0;
  xi[2] =  0.0;
}



/**
   Function returns numbers of integration point closest to element nodes.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipnum - array of local indeces of integration points closest to element nodes
   
   @retval Indices of integration points closest to element nodes are stored in ipnum.

   JK, 25.9.2004
*/
void barelq2d::nodipnum (long eid,long ri,long ci,ivector &ipnum)
{
  long i;
  
  i=Mt->elements[eid].ipp[ri][ci];

  ipnum[0]=i;
  ipnum[1]=i+2;
  ipnum[2]=i+1;
}



/**
   Function computes stresses at main integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barelq2d::res_mainip_stresses (long lcid,long eid)
{
  long i;
  
  //  loop over blocks
  for (i=0;i<nb;i++){
    mainip_stresses (lcid,eid,0,0,i);
  }
}



/**
   Function computes stresses in main integration points of element.

   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param ii - number of actual block

   10.5.2002
*/
void barelq2d::mainip_stresses (long lcid,long eid,long ri,long ci,long ii)
{
  long i,jj,ipp;
  vector eps,epst,epstt,sig,auxsig;
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  reallocv (RSTCKVEC(ncomp[ii],sig));
  reallocv (RSTCKVEC(ncomp[ii],auxsig));
  
  ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
  
  for (i=0;i<intordsm[ii][ii];i++){
    
    Mm->matstiff (d,ipp);
    
    fillv (0.0,sig);
    for (jj=0;jj<nb;jj++){
      reallocv (RSTCKVEC(ncomp[jj],eps));
      
	//  block of strains
      Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
      
      //  stress contributions
      mxv (d,eps,auxsig);
      //  summation of contributions
      addv (auxsig,sig,sig);

    }
    
    //  storage of block of stress
    Mm->storestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
    
    ipp++;
  }
  
}



/**
   Function computes stresses at nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void barelq2d::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i, j, ii, jj, ipp, lipp;;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double l;
  vector w(ASTCKVEC(tnip)),auxw(ASTCKVEC(tnip)),gp(ASTCKVEC(tnip));
  vector sig(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodipnum (eid,ri,ci,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);

  if (Mp->stressaver==2)
  {
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
        if (intordsm[ii][jj]==0)
          continue;
        gauss_points (gp.a, auxw.a, intordsm[ii][jj]);
        lipp = 0;
        for (i=0;i<intordsm[ii][jj];i++){
          w[lipp] = 0.5*auxw[i];
          lipp++;
        }
      }
    }
  }
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress(lcid, ipnum[i], sig);
    
    j=nod[i];
    //  storage of stresses to the node
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress(lcid, 0, sig);
    if (Mp->stressaver==2)
    {
      l = Mt->give_length(eid)*w[ipnum[i]-ipp];
      cmulv(l, sig);
      Mt->nodes[j].storestress(lcid, 0, sig);
    }
  }
}



/**
   Function computes stresses in nodes directly.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param stra - strains at nodes, stra[i][j] means eps[j] at node i
   @param stre - stresses at nodes - output parameter, stre[i][j] means sig[j] at node i

   @retval Function returns computed stresses at parameter stre.

   JK, 25.9.2004
*/
void barelq2d::nod_stresses_comp (long /*lcid*/,long eid,long ri,long ci,double **stra,double **stre)
{
  long i,j,ipp;
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  //  number of the first integration point on the element
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    for (j=0;j<eps.n;j++){
      eps[j]=stra[i][j];
    }
    mxv (d,eps,sig);
    for (j=0;j<eps.n;j++){
      stre[i][j]=sig[j];
    }
  }
}



/**
   Function computes all stress components at all integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void barelq2d::res_allip_stresses (long lcid,long eid)
{
  //  all stress components at all intergation points
  allip_stresses (lcid,eid,0,0);
}



/**
   Function computes all stress components at all integration points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK
*/
void barelq2d::allip_stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  //  blocks of stress components at integration points
  res_mainip_stresses (lcid,eid);
}



/**
   Function computes stresses at required positions.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices

*/
void barelq2d::stresses (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  vector coord,sig;

  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    allip_stresses (lcid,eid,ri,ci);
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
    reallocv (RSTCKVEC(1,coord));
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],0,ncp,sig);

      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    print_err("unknown stress point is required on element (eid=%ld)", __FILE__, __LINE__, __func__, eid+1);
  }
  }

}



/**
   Function computes other values in nodes of element.

   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void barelq2d::nod_other_ip (long eid,long ri,long ci)
{
  long i,j,ncompo,ii,jj,ipp,lipp;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  double l;
  vector other;
  vector w(ASTCKVEC(tnip)),auxw(ASTCKVEC(tnip)),gp(ASTCKVEC(tnip));
  
  //  numbers of integration points closest to nodes
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodipnum(eid,ri,ci,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes(eid,nod);
  
  if (Mp->otheraver==2)
  {
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
        if (intordsm[ii][jj]==0)
          continue;
        gauss_points(gp.a, auxw.a, intordsm[ii][jj]);
        lipp = 0;
        for (i=0;i<intordsm[ii][jj];i++){
          w[lipp] = 0.5*auxw[i];
          lipp++;
        }
      }
    }
  }

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
      l = Mt->give_length(eid)*w[ipnum[i]-ipp];
      cmulv (l, other);
      Mt->nodes[j].storeother(0, ncompo, l, other);
    }    
  }
}



/**
   Load vector is obtained after premultiplying load %matrix
   by nodal load values.
   
   @param eid - number of element
   @param lm - load %matrix
   
   @retval Function returns load %matrix in parameter lm.

   JK, 12.4.2003
*/
void barelq2d::load_matrix (long eid,matrix &lm)
{
  long i;
  double jac,xi,area;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm));
  vector lx(ASTCKVEC(nne)),s(ASTCKVEC(2));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);

  // global coordinates of nodes
  Mt->give_node_coord2d (x,y,eid);  

  //  local coordinates of nodes
  //  only x coordinates are nonzero
  giveloccoord(x,y,lx);
  //  direction vector
  dirvect (s,x,y);

  //  area of bar cross-section
  Mc->give_area (eid,area);
  
  //  Gauss points
  gauss_points (gp.a,w.a,intordmm);
  
  nullm(lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    
    jac_1d (jac, x, xi);
    
    bf_matrix (n, s, xi);
    
    jac*=w[i]*area;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
}



/**
   Function computes load %matrix of the plane stress rectangular
   finite element with bilinear approximation functions.
   Load vector is obtained after premultiplying load %matrix
   by nodal load values.
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 25.7.2001
*/
void barelq2d::res_load_matrix (long eid,matrix &lm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  matrix tran(ASTCKMAT(2,2));
  
  //  node coordinates in the global coordinate system
  Mt->give_node_coord2d (gx,gy,eid);
  //  transformation matrix from global to element coordinate system (2x2)
  tran_mat (gx, gy, tran);
  
  //  load matrix
  load_matrix (eid, lm);
  
  //  transformation of load matrix to global coordinate system
  lgmatrixtransfblock(lm, tran);
  
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
   Function computes internal forces.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system

   @retval Function returns nodal values of internal forces in %vector ifor.

   12.8.2001
*/
void barelq2d::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &gx,vector &gy)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,gx,gy);
}



/**
   Function computes internal forces for nonlocal models.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system

   @retval Function returns nodal values of nonlocal internal forces in %vector ifor.

   12.8.2001
*/
void barelq2d::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &gx,vector &gy)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,gx,gy);
}



/**
   Function computes increment of internal forces.
   
   @param lcid - load case id
   @param eid  - element id
   @param ri,ci - row and column indices
   @param ifor  - vector of internal forces
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system
   
   @retval Function returns nodal values of increments of internal forces in %vector ifor.

   JK, 29.4.2008
*/
void barelq2d::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &gx,vector &gy)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,gx,gy);
}



/**
   Function computes nodal forces caused by temperature changes or eigenstrains.
   
   @param lcid - load case id
   @param eid  - element id
   @param ri,ci - row and column indices
   @param nfor  - array containing nodal forces
   @param gx - x-coordinates in global system
   @param gy - y-coordinates in global system
   
   @retval Function returns nodal values of forces caused by eigenstrains in %vector nfor.

   30.11.2002, JK
*/
void barelq2d::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &gx,vector &gy)
{
  integratedquant iq;
  iq=eigstress;

  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor,gx,gy);
}





/**
   Function computes resulting internal forces. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting internal forces in parameter ifor.

   12.8.2001
*/
void barelq2d::res_internal_forces (long lcid,long eid,vector &ifor)
{
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  long transf;
  ivector nodes (ASTCKIVEC(nne));

  Mt->give_node_coord2d (gx,gy,eid);  

  internal_forces (lcid,eid,0,0,ifor,gx,gy);

  //  transformation of internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting increment of nonlocal internal forces.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting nonlocal internal forces in parameter ifor.

   TKo 7.2008
*/
void barelq2d::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  Mt->give_node_coord2d (gx,gy,eid);  

  incr_internal_forces (lcid,eid,0,0,ifor,gx,gy);

  //  transformation of nonlocal internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting increments of internal forces.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting increments of internal forces in parameter ifor.

   TKo 7.2008
*/
void barelq2d::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  Mt->give_node_coord2d (gx,gy,eid);  

  incr_internal_forces (lcid,eid,0,0,ifor,gx,gy);

  //  transformation of increments of internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function computes resulting internal forces caused by eigenstrains.

   @param lcid - load case id
   @param eid  - element id
   @param nfor - %vector of nodal forces caused by eigenstrains

   @retval Function returns resulting increments of internal forces in parameter nfor.

   12.8.2001
*/
void barelq2d::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne));
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  Mt->give_node_coord2d (gx,gy,eid);  

  eigstrain_forces (lcid,eid,0,0,nfor,gx,gy);

  //  transformation of forces caused by eigenstrains
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    vector v(ASTCKVEC(ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
   Function integrates selected quantity over the finite element.
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param gx - global x coordinates
   @param gy - global y coordinates
      
   @retval Function returns integrated values at nodes in %vector nv.

   TKo, 7.2008
*/
void barelq2d::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv,vector &gx,vector &gy)
{
  long i,ii,ipp;
  double xi,jac,area;
  vector w,gp,ipv,contr(ASTCKVEC(ndofe));
  matrix gm;
  
  Mc->give_area (eid,area);
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],gp));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));
    reallocv (RSTCKVEC(ncomp[ii],ipv));
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
      
      //  strain-displacement (geometric) matrix
      geom_matrix (gm,gx,gy,xi,jac);

      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      cmulv (jac*w[i]*area,contr);
      
      //  summation
      addv(contr, nv, nv);
      
      ipp++;
    }
  }
}



/**
   Function computes correct stresses at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void barelq2d::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->computenlstresses (ipp,Mm->ip[ipp]);
    
    ipp++;
  }
}


/**
   Function computes correct stresses at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void barelq2d::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->computenlstressesincr (ipp);
    
    ipp++;
  }
}


/**
   Function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barelq2d::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  
	  //  computation of correct stresses
	  if (Mp->strcomp==1)
	    Mm->computenlstresses (ipp,Mm->ip[ipp]);
	  
	  ipp++;
	}
      }
    }
  }
}



/**
   Function computes nonlocal correct stresses at integration points on element.
   
   @param lcid - number of load case
   @param eid  - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void barelq2d::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  
	  //  computation of correct stresses
	  if (Mp->strcomp==1)
	    Mm->compnonloc_nlstresses (ipp);
	  
	  ipp++;
	}
      }
    }
  }
}



/**
   Function computes nonlocal correct stresses at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barelq2d::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
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
void barelq2d::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, ipval;
  vector w, gp, anv(ASTCKVEC(3));
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
          //  value in integration point
          ipval = approx (xi,anv);
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
   Function interpolates the nodal values to the integration points on the element
   quadratic approximation functions are used
   
   @param eid    - element id
   @param nodval - nodal values
   @param ipval  - value at integration points
   
   @retval Function returns approximated values at integration points in the %vector ipval.
*/
void barelq2d::intpointval (long eid,vector &nodval,vector &ipval)
{
  long ii,jj,i,k;
  double xi;
  vector w,gp;

  k=0;
  for (ii = 0; ii < nb; ii++)
  {
    for (jj = 0; jj < nb; jj++)
    {
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      for (i = 0; i < intordsm[ii][jj]; i++)
      {
        xi=gp[i];
        //  value in integration point
        ipval[k] = approx (xi,nodval);
        k++;
      }
    }
  }
}



/**
   Function interpolates the nodal values to the integration points on the element
   linear approximation functions are used
   
   @param eid    - element id
   @param nodval - nodal values
   @param ipval  - value at integration points
   
   @retval Function returns approximated values at integration points in the %vector ipval.
*/
void barelq2d::intpointval2 (long /*eid*/,vector &nodval,vector &ipval)
{
  long ii,jj,i,k;
  double xi;
  vector w,gp;
  vector modnodval(ASTCKVEC(Bar2d->nne));
  
  for (i=0;i<Bar2d->nne;i++){
    modnodval[i]=nodval[i];
  }

  k=0;
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;

      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i = 0; i < intordsm[ii][jj]; i++){
        xi=gp[i];

        //  value in integration point
        ipval[k] = Bar2d->approx (xi,modnodval);

        k++;
      }
    }
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

   TKo, 7.2008
*/
void barelq2d::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long ii,jj,i,aipp;
  double xi;
  vector w,gp;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));

  Mt->give_node_coord2d (x,y,eid);
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      //  integration point id
      aipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        if (aipp==ipp){
          coord[0]=approx (xi,x);
          coord[1]=approx (xi,y);
          coord[2]=0.0;
        }
        aipp++;
      }
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
void barelq2d::ipncoord (long eid,long ipp,vector &ncoord)
{
  long ii,jj,i,aipp;
  double xi;
  vector w,gp;

  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      //  integration point id
      aipp=Mt->elements[eid].ipp[ii][jj];
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        if (aipp==ipp){
          ncoord[0]=xi;
          ncoord[1]=0.0;
          ncoord[2]=0.0;
        }
        aipp++;
      }
    }
  }
}



/**
   Function approximates nonmechanical quantities in nodes of element.
   The nodal values are taken from the closest integration point.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything.
   
   12/06/2012 TKr
*/
void barelq2d::mechq_nodval (long eid,vector &nodval,nontransquant qt)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Mt->elements[eid].ipp[0][0];
  nodip_barq (ipid,intordsm[0][0],ipnum);

  for (i=0;i<nne;i++){

    //copy mechanical quantity from closest int. point
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
void barelq2d::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
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
  nodip_barq (ipid,intordsm[0][0],ipnum);

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
      if (ipid == ipnum[i]) // for the first int. point the eqother values may be rewritten -> take them from the backup
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

  // restore original integration point content of strain/stress/other/eqother arrays
  Mm->storestrain(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.strain);
  Mm->storestress(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.stress);
  Mm->storeother(ipid, 0, ipb.ncompother, ipb.other);
  Mm->storeeqother(ipid, 0, ipb.ncompeqother, ipb.eqother);
  Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc);
}
