#include <math.h>
#include "barelq3d.h"
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



barelq3d::barelq3d (void)
{
  long i;

  nne=3;  ndofe=9;  tncomp=1;  napfun=1;
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



barelq3d::~barelq3d (void)
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
   @param x,y,z - vectors of nodal coordinates
   
   @retval Function returns direction %vector s.

   JK, 28.9.2005
*/
void barelq3d::dirvect (vector &s,vector &x,vector &y,vector &z)
{
  double d;
  
  s[0]=x[1]-x[0];
  s[1]=y[1]-y[0];
  s[2]=z[1]-z[0];
  
  d=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
  
  if (d<zero)
    print_err("zero norm of direction vector",__FILE__,__LINE__,__func__);
  
  s[0]/=d;
  s[1]/=d;
  s[2]/=d;
}

/**
   Function computes local coordinates.
   Coordinates of nodes after transformation of the bar to the x axis.
   
   @param x,y,z - vectors of nodal coordinates
   @param lx - %vector of the x coordinates after transformation
   
  @retval Function returns %vector of nodal coordinates at parameter lx.

   JK, 28.9.2005
*/
void barelq3d::giveloccoord(vector &x,vector &y,vector &z,vector &lx)
{
  double l;
  
  lx[0] = 0.0;
  l = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + (z[1]-z[0])*(z[1]-z[0]));
  lx[1] = l;
  l = sqrt((x[2]-x[0])*(x[2]-x[0]) + (y[2]-y[0])*(y[2]-y[0]) + (z[2]-z[0])*(z[2]-z[0]));
  lx[2] = l;
}


/**
   Function approximates function defined by nodal values.

   @param xi,eta - coordinates on element
   @param nodval - nodal values
   
   @retval Function returns approximated value for the given natural coordinate xi.

   JK, 28.9.2005
*/
double barelq3d::approx (double xi,vector &nodval)
{
  double f;
  vector bf(nne);

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
void barelq3d::bf_matrix (matrix &n,vector &s,double xi)
{
  vector bf(ASTCKVEC(nne));
  
  bf_quad_1d (bf.a,xi);
  
  n[0][0]=bf[0]*s[0];
  n[0][1]=bf[0]*s[1];
  n[0][2]=bf[0]*s[2];
  n[0][3]=bf[1]*s[0];
  n[0][4]=bf[1]*s[1];
  n[0][5]=bf[1]*s[2];
  n[0][6]=bf[2]*s[0];
  n[0][7]=bf[2]*s[1];
  n[0][8]=bf[2]*s[2];

}

/**
   Function returns strain-displacement (geometric) %matrix.

   @param gm - strain-displacement (geometric) %matrix
   @param x,y,z - node coordinates
   @param xi - natural coordinate
   @param jac - Jacobian

   @retval Function returns Jacobian in parameter jac and geometric 
           %matrix in parameter gm. 

   JK, 28.9.2005
*/
void barelq3d::geom_matrix (matrix &gm, vector &x,vector &y,vector &z,double xi,double &jac)
{
  vector s(ASTCKVEC(3)),lx(ASTCKVEC(nne)),dx(ASTCKVEC(nne));
  
  //  direction vector
  dirvect (s,x,y,z);
  //  local coordinates of nodes
  giveloccoord (x,y,z,lx);
  
  dx_bf_quad_1d (dx.a,xi);
  derivatives_1d (dx, jac, lx, xi);
  
  fillm (0.0, gm);
  
  gm[0][0]= dx[0]*s[0];
  gm[0][1]= dx[0]*s[1];
  gm[0][2]= dx[0]*s[2];
  gm[0][3]= dx[1]*s[0];
  gm[0][4]= dx[1]*s[1];
  gm[0][5]= dx[1]*s[2];
  gm[0][6]= dx[2]*s[0];
  gm[0][7]= dx[2]*s[1];
  gm[0][8]= dx[2]*s[2];
}






/**
   Function assembles transformation %matrix T between coordinate system
   defined on element and global coordinate system.
      
   x_local = T x_global
   
   @param x - %vector of nodal x-coordinates in the global coordinate system
   @param y - %vector of nodal y-coordinates in the global coordinate system
   @param tran - tranformation %matrix T
   
   @retval Function returns transformation %matrix in parameter tmat.

   JK    
*/
/*
void barelq3d::give_glob_loc_tmat(vector &x, vector &y, matrix &tmat)
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

   JK
*/
 /*
void barelq3d::give_loc_glob_tmat(vector &x, vector &y, matrix &tmat)
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
   Function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l.
   
   @param nodes - array containing node numbers
   @param tmat - transformation %matrix

   @retval Function returns transformation %matrix in parameter tmat.
*/
void barelq3d::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i;
  fillm (0.0,tmat);

  for (i=0;i<tmat.m;i++){
    tmat[i][i]=1.0;
  }
  
  
  for (i=0;i<nodes.n;i++){
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
   Function computes stiffness %matrix of one element.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y,z - node coordinates in global system
   
   @retval Function returns stiffness %matrix at parameter sm.

   JK, 28.9.2005
*/
void barelq3d::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y,vector &z)
{
  long i,ipp;
  double a,xi,jac;
  vector w,gp;
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  //  setting of stiffness matrix
  fillm (0.0,sm);
  //  area of cross section
  Mc->give_area (eid,a);
  
  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    
    //  strain-displacement matrix
    geom_matrix (gm,x,y,z,xi,jac);
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

   JK, 28.9.2005
*/
void barelq3d::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  nodal coordinates in global system
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  stiffness matrix
  stiffness_matrix (eid,0,0,sm,x,y,z);
  
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
   Function computes mass %matrix of twodimensional bar element.
   
   @param eid - element id
   @param mm - mass %matrix
   @param x,y,z - node coordinates in global system
   
   @retval Function returns mass %matrix at parameter mm.

*/
void barelq3d::mass_matrix (long eid,matrix &mm,vector &x,vector &y,vector &z)
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
  giveloccoord(x,y,z,lx);
  //  direction vector
  dirvect (s,x,y,z);
  
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
void barelq3d::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);  
  
  //  mass matrix
  mass_matrix (eid,mm,x,y,z);
  
  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  ivector nodes (ASTCKIVEC(nne));
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
void barelq3d::res_mainip_strains (long lcid,long eid)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector aux,x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  
  //  coordinates of nodes in global system
  Mt->give_node_coord3d (x,y,z,eid);
  //  nodes of element
  Mt->give_elemnodes (eid,nodes);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  computation of strains
  mainip_strains (lcid,eid,0,0,x,y,z,r);
}



/**
   Function computes strains at integration points of element.
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param lx - %vector of x-coordinates in local system
   @param s - direction %vector
   @param r - %vector of nodal displacements
   
   JK, 29.9.2005
*/
void barelq3d::mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,ii,ipp;
  double xi,jac;
  vector gp,w,eps;
  matrix gm;

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);

    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    //  integration point id
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      
      //  strain-displacement matrix
      geom_matrix (gm,x,y,z,xi,jac);
      mxv (gm,r,eps);

      Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
      ipp++;
    }
  }
}

/**
   Function computes strains at nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 29.9.2005
*/
void barelq3d::nod_strains_ip (long lcid,long eid,long ri,long ci)
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
   Function computes nodal strains directly.
   
   @param lcid - load case id
   @param eid - element id
   @param stra - strains at nodes, stra[i][j] means eps[j] at node i
   
   @retval Function returns resulting strains at parameter stra.

   JK, 29.9.2005
*/
void barelq3d::nod_strains_comp (long lcid,long eid,double **stra)
{
  long i,j;
  double jac;
  ivector nodes(nne);
  //vector gx(nne),gy(nne),gz(nne),x(nne),s(3),nxi(nne),r(ndofe),eps(tncomp),aux;
  vector x(nne),y(nne),z(nne),nxi(nne),r(ndofe),eps(tncomp),aux;
  matrix tmat,gm(tncomp,ndofe);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  x-coordinates in local system
  //giveloccoord (gx,gy,gz,x);
  //  computation of direction vector
  //dirvect (s,gx,gy,gz);
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
  nodecoord (nxi);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  strain-displacement (geometric) matrix
    geom_matrix (gm,x,y,z,nxi[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    
    for (j=0;j<eps.n;j++){
      stra[i][j]=eps[j];
    }
  }
}



/**
   Function computes all strain components at all integration points.
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 29.9.2005
*/
void barelq3d::res_allip_strains (long lcid,long eid)
{
  //  all strain components at all integration points
  allip_strains (lcid,eid,0,0);
}



/**
   function assembles all values at all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 29.9.2005
*/
void barelq3d::allip_strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
}



/**
   Function computes strains at strain points.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 29.9.2005
*/
void barelq3d::strains (long lcid,long eid,long ri,long ci)
{
//  long i;
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
    /*
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (ncp,eps);
    reallocv (1,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	//appstrain (lcid,eid,coord[0],0,ncp,eps);

      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    if (Mp->strainaver==0){
      for (i=0;i<nne;i++)
        delete [] stra[i];
      delete [] stra;
    }
    */
    break;
  }
  default:
    print_err("unknown strain point is required",__FILE__,__LINE__,__func__);
  }

}



/**
   Function assembles natural coordinates of nodes of element.

   @param xi - array containing natural coordinates xi

   @retval Function returns natural coordinates in %vector xi.

   JK, 29.9.2005
*/
void barelq3d::nodecoord (vector &xi)
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

   JK, 29.9.2005
*/
void barelq3d::nodipnum (long eid,long ri,long ci,ivector &ipnum)
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
void barelq3d::res_mainip_stresses (long lcid,long eid)
{
  mainip_stresses (lcid,eid,0,0);
}



/**
   function computes stresses at main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   JK, 29.9.2005
*/
void barelq3d::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  vector eps,epst,epstt,sig,auxsig;
  matrix d(tncomp,tncomp);
  
  for (ii=0;ii<nb;ii++){
    
    reallocv (ncomp[ii],sig);
    reallocv (ncomp[ii],auxsig);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      Mm->matstiff (d,ipp);  
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
        reallocv (ncomp[jj],eps);
	
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
}



/**
   Function computes stresses at nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 29.9.2005
*/
void barelq3d::nod_stresses_ip (long lcid,long eid,long ri,long ci)
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



/**
   Function computes nodal stresses directly.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param stra - strains at nodes, stra[i][j] means eps[j] at node i
   @param stre - stresses at nodes - output parameter, stre[i][j] means sig[j] at node i
   
   @retval Function returns computed stresses at parameter stre.

   JK, 25.9.2004
*/
void barelq3d::nod_stresses_comp (long /*lcid*/,long eid,long ri,long ci,double **stra,double **stre)
{
  long i,j,ipp;
  vector eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
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
void barelq3d::res_allip_stresses (long lcid,long eid)
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
void barelq3d::allip_stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
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
void barelq3d::stresses (long lcid,long eid,long ri,long ci)
{
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
    /*
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (ncp,sig);
    reallocv (1,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	//appval (coord[0],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],0,ncp,sig);

      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    */
    break;
  }
  default:
    print_err("unknown stress point is required",__FILE__,__LINE__,__func__);
  }

}



/**
   Function computes other values in nodes of element.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void barelq3d::nod_eqother_ip (long lcid,long eid,long ri,long ci)
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
    
    //  storage of other values to the node
    j=nod[i];
    Mt->nodes[j].storeother (lcid,0,ncompo,eqother);
    
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
void barelq3d::load_matrix (long /*eid*/,matrix &/*lm*/)
{
  /*
  long i;
  double jac,xi,area;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_area (eid,area);

  Mt->give_node_coord2d (x,y,eid);
  

  gauss_points (gp.a,w.a,intordmm);
  
  fillm (0.0,lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    
    jac_1d (jac,x,xi);
    bf_matrix (n,xi);
    
    jac*=w[i]*area;
    
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
  */
}



/**
   Function computes internal forces.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x - x-coordinates in local system
   @param s - direction %vector

   @retval Function returns nodal values of internal forces in %vector ifor.

   TKo 7.2008
*/
void barelq3d::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   Function computes internal forces for nonlocal models.

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x - x-coordinates in local system
   @param s - direction %vector

   @retval Function returns nodal values of nonlocal internal forces in %vector ifor.

   TKo 7.2008
*/
void barelq3d::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);

  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   Function computes increment of internal forces.
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x - x-coordinates in local system
   @param s - direction %vector
   
   @retval Function returns nodal values of increments of internal forces in %vector ifor.

   TKo 7.2008
*/
void barelq3d::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   Function computes nodal forces caused by temperature changes.
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x - x-coordinates in local system
   @param s - direction %vector
   
   @retval Function returns nodal values of forces caused by eigenstrains in %vector nfor.

   TKo 7.2008
*/
void barelq3d::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=eigstress;

  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,nfor,x,y,z);
}



/**
   Function computes resulting internal forces. If it is required, nodal values are transformed to 
   the local coordinate systems.

   @param lcid - load case id
   @param eid  - element id
   @param ifor - %vector of internal forces at nodes

   @retval Function returns resulting internal forces in parameter ifor.

   JK, 29.9.2005
*/
void barelq3d::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  nodal coordinates in global system
  Mt->give_node_coord3d (x,y,z,eid);
  
  internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
  //  transformation of internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    vector v(ASTCKVEC(ndofe));
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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

   JK, 29.9.2005
*/
void barelq3d::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  nodal coordinates in global system
  Mt->give_node_coord3d (x,y,z,eid);
  
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
  //  transformation of nonlocal internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    vector v(ASTCKVEC(ndofe));
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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

   JK, 29.9.2005
*/
void barelq3d::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));

  //  nodal coordinates in global system
  Mt->give_node_coord3d (x,y,z,eid);
  
  incr_internal_forces (lcid,eid,0,0,ifor,x,y,z);

  //  transformation of increments of internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    vector v(ASTCKVEC(ndofe));
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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

   JK, 29.9.2005
*/
void barelq3d::res_eigstrain_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));

  //  nodal coordinates in global system
  Mt->give_node_coord3d (x,y,z,eid);
  
  eigstrain_forces (lcid,eid,0,0,ifor,x,y,z);
  
  //  transformation of increments of internal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    vector v(ASTCKVEC(ndofe));
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}



/**
   Function integrates selected quantity over the finite element
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x - local x coordinate
   @param s - directional vector
   
   @retval Function returns integrated values at nodes in %vector nv.

   TKo, 7.2008
*/
void barelq3d::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z)
{
  long i,ii,ipp;
  double xi,jac,area;
  vector w,gp,ipv,contr(ndofe);
  matrix gm;
  
  
  Mc->give_area (eid,area);
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],ipv);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);

      //  strain-displacement (geometric) matrix
      geom_matrix (gm,x,y,z,xi,jac);

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
   
   TKo, 7.2008
*/
void barelq3d::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
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
   Function computes correct stresses at integration points on element.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void barelq3d::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
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
   Function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void barelq3d::local_values (long /*lcid*/,long eid,long ri,long ci)
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
   Function computes nonlocal correct stresses at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 7.2008
*/
void barelq3d::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
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
   Function computes nonlocal correct stresses at integration points on element.
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void barelq3d::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
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
void barelq3d::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, ipval;
  vector w, gp, anv(3);
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
void barelq3d::intpointval (long eid,vector &nodval,vector &ipval)
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
      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);
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
void barelq3d::intpointval2 (long /*eid*/,vector &nodval,vector &ipval)
{
  long ii,jj,i,k;
  double xi;
  vector w,gp;
  vector modnodval(Bar2d->nne);
  
  for (i=0;i<Bar2d->nne;i++){
    modnodval[i]=nodval[i];
  }

  k=0;
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;

      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);

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
void barelq3d::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long ii,jj,i,aipp;
  double xi;
  vector w,gp;
  vector x(nne),y(nne),z(nne);

  Mt->give_node_coord3d (x,y,z,eid);
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      //  integration point id
      aipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        if (aipp==ipp){
          coord[0]=approx (xi,x);
          coord[1]=approx (xi,y);
          coord[2]=approx (xi,z);
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
void barelq3d::ipncoord (long eid,long ipp,vector &ncoord)
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
      reallocv (intordsm[ii][jj],gp);
      reallocv (intordsm[ii][jj],w);
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
