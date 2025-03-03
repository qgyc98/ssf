#include "quadtet.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "lintet.h"
#include "linhex.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"
#include "intpoints.h"
#include <stdlib.h>



/**
  The constructor inializes attributes to zero values.
  
  Created by Tomas Koudelka, 11.2008
*/
quadtet::quadtet(void)
{
  long i,j;

  //  number of nodes on element
  nne=10;
  //  number of DOFs on element
  ndofe=30;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=4;
  //  number of edges on element
  ned=6;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;
  //  number of surfaces on element
  nsurf=4;
  //  number of nodes on one surface
  nnsurf=6;
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
  intordsm = new long* [nb];
  for (i=0;i<nb;i++)
  {
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=4;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++)
  {
    for (j=0;j<nb;j++)
      tnip+=nip[i][j];
  }

  intordsm[0][0]=4;

}



/**
  Destructor releases allocated memory of the quadtet object.

  Created by Tomas Koudelka, 11.2008
*/
quadtet::~quadtet(void)
{
  long i;
  
  for (i=0;i<nb;i++)
  {
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;
  
  delete [] cncomp;
  delete [] ncomp;
}





/**
  Function approximates quantity defined by nodal values to the given point.
   
  @param xi,eta,zeta - natural coordinates of required point
  @param nodval - nodal values

  @return The function returns approximated quantity in the given point.
   
  Created by Tomas Koudelka, 11.2008
*/
double quadtet::approx(double xi, double eta, double zeta, vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_tet(bf.a,xi,eta,zeta);
  
  scprd(bf,nodval,f);
  
  return f;
}



/**
  Function assembles %matrix of base functions in the required point.
   
  @param n - %matrix of base functions (output)
  @param xi,eta,zeta - natural coordinates of required point

  @return The function return required %matrix in the parametr n

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::bf_matrix(matrix &n, double xi, double eta, double zeta)
{
  long i,j,k,l;
  vector bf(nne);
  
  fillm(0.0,n);

  bf_quad_tet(bf.a,xi,eta,zeta);
  
  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];  j+=3;
    n[1][k]=bf[i];  k+=3;
    n[2][l]=bf[i];  l+=3;
  }
}



/**
  Function computes strain-displacement (geometrical) %matrix in the given point.

  @param gm - geometric %matrix
  @param x,y,z - vectors containing element node coordinates
  @param xi,eta,zeta - natural coordinates of required point
  @param jac - Jacobian (output)

  @return The function returns required geometrical %matrix in the parameter gm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::geom_matrix(matrix &gm, vector &x, vector &y, vector &z,
			  double xi, double eta, double zeta, double &jac)
{
  long i,j,k,l;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_quad_tet(dx.a,xi,eta,zeta);
  dy_bf_quad_tet(dy.a,xi,eta,zeta);
  dz_bf_quad_tet(dz.a,xi,eta,zeta);

  derivatives_3d(dx,dy,dz,jac,x,y,z,xi,eta,zeta);

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
  Function assembles transformation %matrix from local nodal coordinate
  system to the global coordinate system x_g = T x_l
   
  @param nodes - nodes of element
  @param tmat - transformation %matrix (output)

  @return The function returns required transformation %matrix in the parameter tmat.
   
  Created by Tomas Koudelka, 11.2008
*/
void quadtet::transf_matrix(ivector &nodes, matrix &tmat)
{
  long i,n,m;

  fillm(0.0,tmat);

  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++)
    tmat[i][i]=1.0;
  
  for (i=0;i<n;i++)
  {
    if (Mt->nodes[nodes[i]].transf>0)
    {
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

  @param eid - number of element
  @param ri,ci - row and column indices
  @param sm - stiffness %matrix (output)

  @return The function returns required stiffness %matrix in the parameter sm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::stiffness_matrix(long eid, long ri, long ci, matrix &sm)
{
  long i,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp1,gp2,gp3;
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  Mt->give_node_coord3d(x,y,z,eid);
  fillm(0.0,sm);

  reallocv(intordsm[0][0],w);
  reallocv(intordsm[0][0],gp1);
  reallocv(intordsm[0][0],gp2);
  reallocv(intordsm[0][0],gp3);
  
  gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    xi   = gp1[i];
    eta  = gp2[i];
    zeta = gp3[i];

    //  geometric matrix
    geom_matrix(gm,x,y,z,xi,eta,zeta,jac);
	
    //  stiffness matrix of the material
    Mm->matstiff(d,ipp);  ipp++;
    jac*=w[i];
    bdbjac(sm,gm,d,gm,jac);
  }
}



/**
  Function computes stiffness %matrix of one element. If it is required, nodal values are transformed to 
  the local coordinate systems.

  @param eid - number of element
  @param sm - stiffness %matrix (output)

  @return The function returns required stiffness %matrix in the parameter sm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_stiffness_matrix(long eid, matrix &sm)
{
  long transf;
  ivector nodes(nne);
  
  stiffness_matrix(eid,0,0,sm);
  
  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid,nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ndofe,ndofe);
    transf_matrix(nodes,tmat);
    glmatrixtransf(sm,tmat);
  }
}



/**
  Function computes mass %matrix. If it is required, nodal values are transformed to 
  the local coordinate systems.
   
  @param eid - number of element
  @param mm - mass %matrix

  @return The function returns required mass %matrix in the parameter mm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::mass_matrix(long eid, matrix &mm)
{
  long i;
  double jac,xi,eta,zeta,rho;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp1(intordmm);
  vector gp2(intordmm),gp3(intordmm),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes(eid,nodes);
  Mc->give_density(eid,nodes,dens);

  Mt->give_node_coord3d(x,y,z,eid);
  gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordmm);
  
  fillm (0.0,mm);
  
  for (i=0;i<intordmm;i++)
  {
    xi  =gp1[i];
    eta =gp2[i];
    zeta=gp3[i];

    jac_3d(jac,x,y,z,xi,eta,zeta);
    bf_matrix(n,xi,eta,zeta);
    rho = approx(xi,eta,zeta,dens);
    jac *= w[i]*rho;
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
}



/**
  Function computes mass %matrix without possible transformations at nodes.
   
  @param eid - number of element
  @param mm - mass %matrix

  @return The function returns required mass %matrix in the parameter mm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_mass_matrix(long eid, matrix &mm)
{
  long transf;
  ivector nodes(nne);
  
  mass_matrix(eid,mm);

  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat(ndofe,ndofe);
    transf_matrix(nodes,tmat);
    glmatrixtransf(mm,tmat);
  }
}



/**
  Function computes load %matrix without possible transformations at nodes.
   
  @param eid - number of element
  @param lm - load %matrix
   
  @return The function returns required load %matrix in the parameter lm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::load_matrix(long eid, matrix &lm)
{
  long i;
  double jac,xi,eta,zeta;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne);
  vector w(intordmm),gp1(intordmm),gp2(intordmm),gp3(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes(eid,nodes);
  Mt->give_node_coord3d(x,y,z,eid);
  gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordmm);
  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++)
  {
    xi  =gp1[i];  
    eta =gp2[i];
    zeta=gp3[i];	

    jac_3d (jac,x,y,z,xi,eta,zeta);
    bf_matrix (n,xi,eta,zeta);
    jac *= w[i];
    nnj (lm.a,n.a,jac,n.m,n.n);
  }
}



/**
  Function computes load %matrix. If it is required, nodal values are transformed to 
  the local coordinate systems.
   
  @param eid - number of element
  @param lm - load %matrix
   
  @return The function returns required load %matrix in the parameter lm.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_load_matrix(long eid, matrix &lm)
{
  long transf;
  ivector nodes(nne);

  load_matrix (eid,lm);
  
  //  transformation of load matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }

}



/**
  Function computes volume appropriate to integration point.
   
  @param eid - element id
  @param w   - weight of given integration point

  @return Function returns volume appropriate to integration point.

  Created by Tomas Koudelka, 3.2009
*/
double quadtet::volumeip (long eid, double w)
{
  double vol;
  
  vol = Mt->give_volume(eid)*w;
  
  return vol;
}



/**
  Function computes strains at integration points of given element.
   
  @param lcid - load case id
  @param eid - element id
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_ip_strains(long lcid, long eid)
{
  vector x(nne),y(nne),z(nne),r(ndofe),aux;
  ivector nodes(nne);
  matrix gm,tmat;

  Mt->give_elemnodes(eid,nodes);
  Mt->give_node_coord3d(x,y,z,eid);
  eldispl(lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems(nodes);
  if (transf>0){
    reallocv(ndofe,aux);
    reallocm(ndofe,ndofe,tmat);
    transf_matrix(nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf(aux,r,tmat);
    copyv(aux,r);
  }

  ip_strains(lcid,eid,0,0,x,y,z,r);
}



/**
  Function computes strains at integration points of given element.
  Used displacements are passed over parameter r.
   
  @param lcid - load case id
  @param eid - element id
  @param ri - row index
  @param ci - column index
  @param x,y,z - %vectors of nodal coordinates
  @param r - %vector of nodal displacements
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::ip_strains(long lcid, long eid, long ri, long ci, vector &x, vector &y, vector &z, vector &r)
{
  long i,ii,ipp;
  double xi,eta,zeta,jac;
  vector gp1,gp2,gp3,w,eps,aux;
  ivector nodes(nne),cn(ndofe);
  matrix gm,tmat;

  for (ii=0;ii<nb;ii++)
  {
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],gp3);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);
    gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];

    for(i=0;i<intordsm[ii][ii];i++)
    {
      xi  =gp1[i];
      eta =gp2[i];
      zeta=gp3[i];
	  
      geom_matrix(gm,x,y,z,xi,eta,zeta,jac);
      mxv(gm,r,eps);
      Mm->storestrain(lcid,ipp,cncomp[ii],ncomp[ii],eps);
      ipp++;
    }
  }
}



/**
  Function computes strains at nodes of given element. 
  The strains are copied to the node from the nearest
  integration point and they are averaged if it is 
  prescribed.

  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nod_strains_ip(long lcid, long eid, long ri, long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  double vol;
  vector eps(tncomp),w,gp1,gp2,gp3;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadtet(ipp,ipnum);
  
  if (Mp->strainaver==2)
  {
    reallocv(RSTCKVEC(intordsm[0][0],w));
    reallocv(RSTCKVEC(intordsm[0][0],gp1));
    reallocv(RSTCKVEC(intordsm[0][0],gp2));
    reallocv(RSTCKVEC(intordsm[0][0],gp3));
    gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
  }

  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain(lcid,ipnum[i],eps);
    
    j=nod[i];
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid,0,eps);
    if (Mp->strainaver==2)
    {
      vol=volumeip (eid,w[ipnum[i]-ipp]);
      cmulv (vol,eps,eps);
      Mt->nodes[j].storestrain(lcid,0,vol,eps);
    }
  }
}



/**
  Function computes strains at nodes of element.
   
  @param lcid - load case id
  @param eid - element id
  @param stra - array for strain components (output)
   
  stra[i][j] - the j-th strain component at the i-th node

  @return The function returns computed nodal strains in the parameter stra

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nod_strains_comp(long lcid, long eid, double **stra)
{
  long i,j;
  double jac;
  vector x(nne),y(nne),z(nne),nxi(nne),neta(nne),nzeta(nne),eps(tncomp),r(ndofe),aux;
  ivector nodes(nne);
  matrix gm(tncomp,ndofe),tmat;

  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0)
  {
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_quadtet (nxi,neta,nzeta);

  //  loop over nodes
  for (i=0;i<nne;i++)
  {
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
  Function computes strains at required positions.
   
  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::strains(long lcid, long eid, long ri, long ci)
{
  long i,naep,ncp,sid;
  vector coord,eps;
  double **stra=NULL;
  
  if (Mp->strainaver==0)
  {
    stra = new double* [nne];
    for (i=0;i<nne;i++)
      stra[i] = new double [tncomp];
    //elem_strains (stra,lcid,eid,ri,ci);
  }

  switch (Mm->stra.tape[eid])
  {
    case nowhere:
      break;
    case intpts:
      //  strains are computed at integration points
      res_ip_strains (lcid,eid);
      break;
    case enodes:
      //  strains are copied to nodes from the closest integration points
      nod_strains_ip (lcid,eid,ri,ci);
      break;
    case cenodes:
      //  strains are computed at nodes
      nod_strains_comp (lcid,eid,stra);
      break;
    case userdefined:
    {
      //  number of auxiliary element points
      naep = Mm->stra.give_naep (eid);
      ncp = Mm->stra.give_ncomp (eid);
      sid = Mm->stra.give_sid (eid);
      reallocv(ncp,eps);
      reallocv(3,coord);
      for (i=0;i<naep;i++)
      {
        Mm->stra.give_aepcoord (sid,i,coord);

        if (Mp->strainaver==0)
	  //appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
        if (Mp->strainaver==1)
	  //appstrain (lcid,eid,coord[0],coord[1],coord[2],0,ncp,eps);
      
        Mm->stra.storevalues(lcid,eid,i,eps);
      }
      break;
    }
    default:
      print_err("unknown strain point is required",__FILE__,__LINE__,__func__);
    
  }
  if (Mp->strainaver==0)
  {
    for (i=0;i<nne;i++)
      delete [] stra[i];

    delete [] stra;
  }
}



/**
  Function computes stresses at integration points of element.
   
  @param lcid - load case id
  @param eid - element id
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_ip_stresses(long lcid, long eid)
{
  ip_stresses (lcid,eid,0,0);
}



/**
  Function computes stresses at integration points of element
  stresses are computed by material models.
   
  @param lcid - load case id
  @param eid - element id
  @param ri - row index
  @param ci - column index
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::ip_stresses(long /*lcid*/, long eid, long ri, long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->computenlstresses (ipp,Mm->ip[ipp]);
    ipp++;
  }
}



/**
  Function computes stresses at integration points of element
  stresses are computed from strains with the help of elastic stiffness.
   
  @param lcid - load case id
  @param eid - element id
  @param ri - row index
  @param ci - column index
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::ip_elast_stresses(long lcid, long eid, long ri, long ci)
{
  long i,ipp;
  vector eps,sig(tncomp);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    //  stiffness matrix of the material
    Mm->matstiff(d,ipp);
    //  strains
    Mm->givestrain(lcid,ipp,eps);
    //  elastic stresses
    mxv(d,eps,sig);
    Mm->storestress(lcid,ipp,sig);
    ipp++;
  }
}



/**
  Function computes stresses at nodes of element.

  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nod_stresses_ip(long lcid, long eid, long ri, long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  double vol;
  vector sig(tncomp),w,gp1,gp2,gp3;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadtet (ipp,ipnum);
  
  if (Mp->stressaver==2)
  {
    reallocv(intordsm[0][0],w);
    reallocv(intordsm[0][0],gp1);
    reallocv(intordsm[0][0],gp2);
    reallocv(intordsm[0][0],gp3);
    gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
  }

  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++)
  {
    //  stresses at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    //  storage of stresses to the node
    j=nod[i];
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2)
    {
      vol=volumeip (eid,w[ipnum[i]-ipp]);
      cmulv (vol,sig,sig);
      Mt->nodes[j].storestress (lcid,0,vol,sig);
    }
  }
}



/**
  Function computes stresses in nodes.
   
  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices
  @param stra - strains at nodes, stra[i][j] means eps[j] at node i
  @param stre - stresses at nodes - output parameter, stre[i][j] means sig[j] at node i (output)

  @return The function returns nodal stresses in the parameter stre.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nod_stresses_comp(long /*lcid*/, long eid, long ri, long ci, double **stra, double **stre)
{
  long i,j,ipp;
  ivector ipnum(nne);
  vector eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  //  number of the first integration point on the element
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodip_quadtet (ipp,ipnum);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  stiffness matrix of the material is taken 
    // from the nearest integration point (not exact, it is valid only for elastic material)
    Mm->matstiff (d,ipnum[i]);
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
  Function computes stresses at required positions.
   
  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::stresses(long lcid, long eid, long ri, long ci)
{
  long i,naep,ncp,sid;
  vector coord,sig;
  
  switch (Mm->stre.tape[eid])
  {
    case nowhere:
      break;
    case intpts:
      //  stresses are computed at integration points
      res_ip_stresses (lcid,eid);
      break;
    case enodes:
      //  stresses are copied to nodes from the closest integration points
      nod_stresses_ip (lcid,eid,ri,ci);
      break;
    case userdefined:
    {
      //  number of auxiliary element points
      naep = Mm->stre.give_naep (eid);
      ncp = Mm->stre.give_ncomp (eid);
      sid = Mm->stre.give_sid (eid);
      reallocv (ncp,sig);
      reallocv (3,coord);
      for (i=0;i<naep;i++)
      {
        Mm->stre.give_aepcoord (sid,i,coord);

        if (Mp->stressaver==0)
          //appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
        if (Mp->stressaver==1)
         //appstress (lcid,eid,coord[0],coord[1],coord[2],0,ncp,sig);
        Mm->stre.storevalues(lcid,eid,i,sig);
      }
      break;
    }
    default:
      print_err("unknown stress point is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function computes other values in nodes of element.

  @param eid - element id
  @param ri,ci - row and column indices
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nod_other_ip (long eid,long ri,long ci)
{
  long i, j, ipp, ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double vol;
  vector other, w, gp1, gp2, gp3;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadtet(ipp,ipnum);

  if (Mp->otheraver==2)
  {
    reallocv(RSTCKVEC(intordsm[0][0], w));
    reallocv(RSTCKVEC(intordsm[0][0],gp1));
    reallocv(RSTCKVEC(intordsm[0][0],gp2));
    reallocv(RSTCKVEC(intordsm[0][0],gp3));
    gauss_points_tet(gp1.a, gp2.a, gp3.a, w.a, intordsm[0][0]);
  }

  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    //Mm->givestrain (lcid,ipnum[i],eps);
    
    ncompo = Mm->givencompother(ipnum[i], 0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    
    j=nod[i];
    //  storage of other values to the node
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother(0, ncompo, other);
    if (Mp->otheraver==2)
    {
      vol=volumeip(eid, w[ipnum[i]-ipp]);
      cmulv(vol, other, other);
      Mt->nodes[j].storeother(0, ncompo, vol, other);
    }
  }
}



/**
  Function computes internal forces in the case of geometrical linear computation.
   
  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
  @param ifor - %vector of internal forces (output)
  @param x,y,z - %vectors of nodal coordinates
   
  @return The function returns %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x, vector &y, vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
  Function computes nonlocal internal forces.

  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
  @param ifor - %vector of internal forces (output)
  @param x,y,z - %vectors of nodal coordinates
   
  @return The function returns %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of nonlocal stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);

}



/**
  Function computes increments of internal forces.

  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices
  @param ifor - %vector of increments of internal forces (output)
  @param x,y,z - %vectors of nodal coordinates

  @return The function returns %vector of increments of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
  Function computes contributions from eigenstrains to the right hand side
   
  @param lcid - load case id
  @param eid - element id
  @param ri,ci - row and column indices
  @param ifor - %vector of nodal forces due to eigenstrains
  @param x,y,z - %vectors of nodal coordinates

  @return The function returns %vector of nodal forces in the parameter nfor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
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
  Function computes resulting internal forces (from correct stresses)
   
  @param lcid - number of load case
  @param eid - element id
  @param ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
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
  Function computes resulting internal forces for nonlocal models
   
  @param lcid - load case id
  @param eid - element id
  @param ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
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
  Function computes increments of internal forces.

  @param lcid - load case id
  @param eid - element id
  @param ifor - %vector of increments of internal forces (output)

  @return The function returns resulting %vector of increments of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);

  incr_internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
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
  Function computes resulting contributions from eigenstrains to the right hand side
   
  @param eid - element id
  @param ifor - %vector of internal forces

  @return The function returns resulting %vector of nodal forces in the parameter nfor.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);

  Mt->give_node_coord3d (x,y,z,eid);

  eigstrain_forces (lcid,eid,0,0,nfor,x,y,z);

  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //globloctransf (nfor,v,tmat);
    glvectortransf (nfor,v,tmat);
    copyv (v,nfor);
  }
}



/**
  Function computes correct stresses at integration points on element.

  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for(i=0;i<intordsm[0][0];i++)
  {
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->computenlstresses(ipp,Mm->ip[ipp]);
    ipp++;
  }
}



/**
  Function computes correct increments of stresses at integration points on element.

  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
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
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    //  computation of local values
    if (Mp->strcomp==1)
      Mm->computenlstresses (ipp,Mm->ip[ipp]);
    ipp++;
  }
}



/**
  Function computes nonlocal correct stresses at integration points on element
   
  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->compnonloc_nlstresses (ipp);
    ipp++;
  }
}



/**
  Function computes nonlocal correct stresses at integration points on element.
   
  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
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
  Function integrates selected quantity over the finite element.
   
  @param iq - type of integrated quantity (see alias.h)
  @param lcid - number of load case
  @param eid - element id
  @param ri,ci - row and column indices
  @param nv - nodal values
  @param x,y,z - node coordinates
   
  @return Function returns integrated values at nodes in %vector nv.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x, vector&y, vector &z)
{
  long i,ipp;
  double xi,eta,zeta,jac;
  vector w, gp1, gp2, gp3;
  vector ipv(tncomp),contr(ndofe);
  matrix gm(tncomp,ndofe);

  fillv(0.0,nv);
  reallocv(intordsm[0][0],w);
  reallocv(intordsm[0][0],gp1);
  reallocv(intordsm[0][0],gp2);
  reallocv(intordsm[0][0],gp3);
  gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    xi  =gp1[i];
    eta =gp2[i];
    zeta=gp3[i];
    //  function assembles required quantity at integration point
    Mm->givequantity(iq,lcid,ipp,cncomp[0],ipv);
    //  strain-displacement (geometric) matrix
    geom_matrix(gm,x,y,z,xi,eta,zeta,jac);
    //  contribution to the internal forces
    mtxv(gm,ipv,contr);
    cmulv(jac*w[i],contr);
    //  summation
    addv(contr,nv,nv);
    ipp++;
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

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  double xi,eta,zeta;
  vector x(nne),y(nne),z(nne),w(intordsm[0][0]);
  vector gp1(intordsm[0][0]),gp2(intordsm[0][0]),gp3(intordsm[0][0]);
  
  gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++)
  {
    xi  =gp1[i];
    eta =gp2[i];
    zeta=gp3[i];
    if (ii==ipp)
    {
      coord[0]=approx (xi,eta,zeta,x);
      coord[1]=approx (xi,eta,zeta,y);
      coord[2]=approx (xi,eta,zeta,z);
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
void quadtet::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ii, ri, ci;
  vector w, gp1, gp2, gp3;
  
  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      if (intordsm[ri][ci] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp3));
      gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[0][0]);
      ii=Mt->elements[eid].ipp[ri][ci];

      for (i=0;i<intordsm[ri][ci];i++)
      {
        if (ii==ipp)
        {
          ncoord[0]=gp1[i];
          ncoord[1]=gp2[i];
          ncoord[2]=gp3[i];
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
void quadtet::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, zeta, ipval;
  vector w, gp1, gp2, gp3, anv(nne);
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
        reallocv (intordsm[ii][jj],gp3);
        reallocv (intordsm[ii][jj],w);
        gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi  = gp1[k];
          eta = gp2[k];
          zeta= gp3[k];
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
  Function computes volume appropriate to integration point.

  @param eid - element id
  @param ri  - row index of the given integration point block
  @param ci  - column index of the given integration point block
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::ipvolume (long eid, long ri, long ci)
{
  long i,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp1,gp2,gp3;
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++)
  {
    for (jj=0;jj<nb;jj++)
    {
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],gp3);
      gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++)
      {
	xi  =gp1[i];
        eta =gp2[i];
        zeta=gp3[i];
        jac_3d (jac,x,y,z,xi,eta,zeta);
        jac=fabs(jac);
        jac*=w[i];
        Mm->storeipvol (ipp,jac);
        ipp++;
      }
    }
  }
}



/**
  Function computes nodal forces caused by surface load.

  @param lcid - load case id   
  @param eid - element id
  @param is - array with surface id 
  @param nv - nodal values of surface load
  @param nf - nodal forces due to surface load (output)
   
  @return The function returns nodal forces caused by surface load in the parameter nf.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::node_forces_surf (long /*lcid*/, long eid, long *is, double *nv, vector &nf)
{
  long i,j,k,ii,transf;
  double xi=0.0,eta=0.0,zeta=0.0,jac;
  ivector nodes(nne), snodes(nnsurf);
  vector x(nne),y(nne),z(nne);
  vector gp1(intordb),gp2(intordb),w(intordb);
  vector v(ndofe),av(ndofe),tnv(nnsurf*napfun);
  
  matrix n(napfun,ndofe);
  matrix am(ndofe,ndofe), tran(3,3);
  matrix tmat;

  Mt->give_elemnodes(eid, nodes);
  Mt->give_node_coord3d(x, y, z, eid);
  gauss_points_tr(gp1.a, gp2.a, w.a, intordb);
  fillv(0.0, nf);
  for(ii=0; ii<4; ii++)
  {
    // is=0 not loading
    if (is[ii]==0)
      continue;
    fillv(0.0, av);
    fillm(0.0, am);
    quadtetrahedral_surfnod(snodes.a, ii);
    switch(ii+1)
    {
      case 1:
      {
        // surface 1 - node 3,2,4,  6,9,10
	xi=0.0;
	for(i=0; i<intordb; i++)
        {
	  eta=gp1[i];  zeta=gp2[i];   
          bf_matrix(n, xi, eta, zeta);
          jac2d_3d(jac, x, y, z, eta, zeta,ii);
	  jac = jac*w[i];     
	  nnj(am.a, n.a, jac, n.m, n.n);
	}
        break;
      }
      case 2:
      {
        // surface 2 - node 1,3,4,  7,10,8
  	eta=0.0;
	for (i=0; i<intordb; i++)
        {
          xi=gp1[i];  zeta=gp2[i];    
          bf_matrix(n,xi,eta,zeta);
          jac2d_3d(jac,x,y,z,xi,zeta,ii);
          jac = jac*w[i];      
	  nnj(am.a,n.a,jac,n.m,n.n);
	}
        break;
      }
      case 3:
      {
        // surface 3 - node 2,1,4,  5,8,9
	zeta=0.0;
	for (i=0;i<intordb;i++)
        {
          xi=gp1[i];  eta=gp2[i];      
          bf_matrix(n,xi,eta,zeta);
          jac2d_3d(jac,x,y,z,xi,eta,ii);
          jac = jac*w[i];      
          nnj (am.a,n.a,jac,n.m,n.n);
	}
        break;
      }
      case 4:
      {
        // surface 4 - node 1,2,3,  5,6,7
	for(i=0;i<intordb;i++)
        {
	  xi=gp1[i];  eta=gp2[i];  zeta=1.0-xi-eta;
          bf_matrix(n,xi,eta,zeta);
          jac2d_3d(jac,x,y,z,xi,zeta,ii);
          jac = jac*w[i];
          nnj(am.a,n.a,jac,n.m,n.n);
	}
        break;
      }
    }
    // Load in GCS
    if (is[ii]==1)
    {
      for(i=0; i<napfun; i++)
      {
        for (j=0; j<nnsurf; j++)
          av[snodes[j]*napfun+i] = nv[ii*napfun*nnsurf+j*napfun+i];
      }
    }
    // Load in LCS
    if (is[ii]==2)
    {
      k = 0;
      for(i=0; i<napfun; i++)
      {
        for (j=0; j<nnsurf; j++)
	{
          av[k] = nv[ii*napfun*nnsurf+j*napfun+i];
          k++;
	}
      }
      locglob_nodeval (ii,av,tnv.a,x,y,z);
      fillv(0.0, av);
      for(i=0; i<napfun; i++)
      {
        for (j=0; j<nnsurf; j++)
          av[snodes[j]*napfun+i] = tnv[ii*napfun*nnsurf+j*napfun+i];
      }
    }
    fillv(0.0,v);
    mxv(am,av,v);
    addv (v,nf,nf);
  }
  // transformation to loacal coordinate systems of nodes 
  transf = Mt->locsystems (nodes);
  if (transf>0)
  {
    reallocm(ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    lgvectortransfblock (nf,tmat);
  }
}



/**
  The function transforms nodal values of surface load expressed in local coordinate system 
  of the given surface is to the global coordinate system.

  @param is - index of surface
  @param nv - %vector of nodal values of surface load (3 components of load are assumed at ech node)
  @param tnv - resulting surface load in the globla coordinate system (output)
  @param x - nodal x-coordinates
  @param y - nodal y-coordinates
  @param z - nodal z-coordinates

  @return The function returns transformed surface load in the parameter tnv.

  Created by Tomas Koudelka,
*/
void quadtet::locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z)
{
  double norm;
  vector g1(3),g2(3),g3(3),lv(3),gv(3);
  matrix t(3,3);
  
  if (is==0){
    g1[0]=x[2]-x[3];
    g1[1]=y[2]-y[3];
    g1[2]=z[2]-z[3];
    
    g2[0]=x[1]-x[3];
    g2[1]=y[1]-y[3];
    g2[2]=z[1]-z[3];
  }
  if (is==1){
    g1[0]=x[0]-x[3];
    g1[1]=y[0]-y[3];
    g1[2]=z[0]-z[3];
    
    g2[0]=x[2]-x[3];
    g2[1]=y[2]-y[3];
    g2[2]=z[2]-z[3];
  }
  if (is==2){
    g1[0]=x[1]-x[3];
    g1[1]=y[1]-y[3];
    g1[2]=z[1]-z[3];
    
    g2[0]=x[0]-x[3];
    g2[1]=y[0]-y[3];
    g2[2]=z[0]-z[3];
  }
  if (is==3){
    g1[0]=x[1]-x[0];
    g1[1]=y[1]-y[0];
    g1[2]=z[1]-z[0];
    
    g2[0]=x[2]-x[0];
    g2[1]=y[2]-y[0];
    g2[2]=z[2]-z[0];
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
  mxv (t.a,nv.a+12,tnv+12,3,3);
  mxv (t.a,nv.a+15,tnv+15,3,3);
}



/**
  Function interpolates the nodal values to the integration points on the element
  quadratic approximation functions are used
   
  @param eid - element id
  @param nodval - nodal values
  @param ipval - value at integration points
   
  @return Function returns approximated values at integration points in the parameter ipval.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::intpointval(long /*eid*/,vector &nodval,vector &ipval)
{
  long i,ii,jj,l;
  double xi,eta,zeta;
  vector w,gp1,gp2,gp3;
  
  l=0;
  for (ii=0;ii<nb;ii++)
  {
    for (jj=0;jj<nb;jj++)
    {
      if (intordsm[ii][jj]==0)  
        continue;
      
      reallocv(intordsm[ii][jj],w);
      reallocv(intordsm[ii][jj],gp1);
      reallocv(intordsm[ii][jj],gp2);
      reallocv(intordsm[ii][jj],gp3);
      gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++)
      {
	xi  =gp1[i];
        eta =gp2[i];
        zeta=gp3[i];
        ipval[l]=approx(xi,eta,zeta,nodval);
        l++;
      }
    }
  }
}



/**
  Function interpolates the nodal values to the integration points on the element
  linear approximation functions are used
   
  @param eid - element id
  @param nodval - nodal values
  @param ipval - value at integration points
  
  @return Function returns approximated values at integration points in the parameter ipval.

  Created by Tomas Koudelka, 11.2008
*/
void quadtet::intpointval2(long /*eid*/,vector &nodval,vector &ipval)
{
  long i,ii,jj,l;
  double xi,eta,zeta;
  vector w,gp1,gp2,gp3;
  vector modnodval(Ltet->nne);
  
  for (i=0;i<Ltet->nne;i++)
    modnodval[i]=nodval[i];
  
  l=0;
  for (ii=0;ii<nb;ii++)
  {
    for (jj=0;jj<nb;jj++)
    {
      if (intordsm[ii][jj]==0)  
        continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      reallocv (intordsm[ii][jj],gp3);
      
      gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++)
      {
	xi  =gp1[i];
        eta =gp2[i];
        zeta=gp3[i];
        ipval[l]=Ltet->approx_nat (xi,eta,zeta,modnodval);
        l++;
      }
    }
  }
}



/**
  Functiion computes averaged strains at the given element.

  @param lcid - load case id   
  @param eid - element id
  @param ri  - row index of the given integration point block
  @param ci  - column index of the given integration point block
  @param averstra - %vector of avreaged values of strains (have to be divided by element volume) (output)
  @param volume - element volume (output)

  @return The function returns averaged strains in the parameter averstra and 
          element volume in the parameter volume.
   
  Created by Tomas Koudelka, 1.2009
*/
void quadtet::aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume)
{
  long i,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),gp1,gp2,gp3,w,eps;

  Mt->give_node_coord3d (x,y,z,eid);
  for (ii=0;ii<nb;ii++)
  {
    if (intordsm[ii][ii]==0)  
      continue;
    
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],gp3);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    gauss_points_tet(gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++)
    {
      xi  =gp1[i];
      eta =gp2[i];
      zeta=gp3[i];
      jac_3d (jac,x,y,z,xi,eta,zeta);
      jac=fabs(jac);
      Mm->givestrain (lcid,ipp,eps);
      volume+=w[i]*jac;
      for(l=0;l<averstra.n;l++)
        averstra[l]+=eps[l]*w[i]*jac;
      ipp++;
    }
  }
}
