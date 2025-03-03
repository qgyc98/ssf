#include "axisymqq.h"
#include "axisymlq.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "loadcase.h"
#include "intpoints.h"
#include <math.h>
#include <stdlib.h>


axisymqq::axisymqq (void)
{
  long i,j;
  
  nne=8;  ndofe=16;  tncomp=4;  napfun=2;  ned=4;  nned=3;
  intordmm=3;  intordb=2;
  ssst=axisymm;

  nb=1;
  
  ncomp = new long [nb];
  ncomp[0]=4;
  
  cncomp = new long [nb];
  cncomp[0]=0;
  
  
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=4;
  
  intordsm[0][0]=2;
  
  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

}

axisymqq::~axisymqq (void)
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
   
*/
double axisymqq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_quad_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function returns matrix of approximation functions
   
   @param n - matrix of approximation functions
   @param xi,eta - natural coordinates
   
   9.7.2001
*/
void axisymqq::bf_matrix (matrix &n,double xi,double eta)
{
  long i,j,k;
  vector bf(ASTCKVEC(nne));
  
  nullm (n);

  bf_quad_4_2d (bf.a,xi,eta);
  
  j=0;  k=1;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];
    n[1][k]=bf[i];
    j+=2;  k+=2;
  }
}

/**
   function assembles geometric matrix
   
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx
   
   @param gm - geometric matrix
   @param ri - block index
   @param x,y - arrays of node coordinates
   @param xi,eta - natural coordinates
   @param jac - jacobian
   
   8.12.2001
*/
void axisymqq::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  double r;
  vector bf(ASTCKVEC(nne)),dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  
  dx_bf_quad_4_2d (dx.a,xi,eta);
  dy_bf_quad_4_2d (dy.a,xi,eta);
  bf_quad_4_2d (bf.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  r = approx (xi,eta,x);
  if (fabs(r)<Mp->zero){
    //fprintf (stderr,"\n\n radius is equal %e in function axisymqq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
    r=0.00001;
  }
  
  nullm (gm);
  
  i1=0;  i2=1;
  for (i=0;i<nne;i++){
    gm[0][i1]=dx[i];
    gm[1][i2]=dy[i];
    gm[2][i1]=bf[i]/r;
    gm[3][i1]=dy[i];
    gm[3][i2]=dx[i];
    i1+=2;  i2+=2;
  }
}

/**
   function assembles part of geometric matrix
   
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx
   
   @param gm - geometric matrix
   @param ri - block index
   @param x,y - arrays of node coordinates
   @param xi,eta - natural coordinates
   @param jac - jacobian
   
   8.12.2001
*/
void axisymqq::geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac)
{
  if (nb==1){
    geom_matrix (gm,x,y,xi,eta,jac);
  }
  else{
    long i,i1,i2;
    double r;
    vector bf(ASTCKVEC(nne)),dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
    
    dx_bf_quad_4_2d (dx.a,xi,eta);
    dy_bf_quad_4_2d (dy.a,xi,eta);
    bf_quad_4_2d (bf.a,xi,eta);
    
    derivatives_2d (dx,dy,jac,x,y,xi,eta);
    
    r = approx (xi,eta,x);
    if (fabs(r)<Mp->zero){
      //fprintf (stderr,"\n\n radius is equal %e in function axisymqq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
      r=0.00001;
    }
    
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
      i1=0;
      for (i=0;i<nne;i++){
	gm[0][i1]=bf[i]/r;
	i1+=2;
      }
    }
    if (ri==2){
      i1=0;  i2=1;
      for (i=0;i<nne;i++){
	gm[0][i1]=dy[i];
	gm[0][i2]=dx[i];
	i1+=2;  i2+=2;
      }
    }
  }
}

/**
   nutno otestovat! pak je mozne smazat tuto hlasku
   
   transformation matrix x_g = T x_l
*/
void axisymqq::transf_matrix (ivector &nodes,matrix &tmat)
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
      tmat[i*2][i*2]   = Mt->nodes[nodes[i]].e1[0];    tmat[i*2][i*2+1]   = Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2] = Mt->nodes[nodes[i]].e1[1];    tmat[i*2+1][i*2+1] = Mt->nodes[nodes[i]].e2[1];
    }
  }
}

/**
   function computes stiffness matrix of axisymmetric quadrilateral
   finite element with bilinear approximation functions
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness matrix

   8.12.2001
*/
void axisymqq::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,ipp,transf;
  double xi,eta,jac,r;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w,gp;
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  
  nullm (sm);
  
  
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  geometric matrix
      geom_matrix (gm,x,y,xi,eta,jac);
      
      if(jac < 0.0){
	//jac = fabs (jac);
	print_err("wrong numbering of nodes on element number %ld, negative volume! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	abort();
      }
  
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      
      r = approx (xi,eta,x);
      jac*=w[i]*w[j]*r;
      
      //  contribution to the stiffness matrix of the element
      bdbjac (sm,gm,d,gm,jac);
      
      ipp++;
    }
  }
  
  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes resulting stiffness matrix of element
   
   @param eid - element id
   @param sm - stiffness matrix
   
   10.5.2002
*/
void axisymqq::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

/**
   function computes mass matrix of the rectangular axisymmetric
   finite element with bilinear approximation functions
   
   @param eid - number of element
   @param mm - mass matrix

   24.6.2001
*/
void axisymqq::mass_matrix (long eid,matrix &mm)
{
  long i,j;
  double jac,xi,eta,rho,r;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne)),dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (mm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      rho = approx (xi,eta,dens);
      r = approx (xi,eta,x);
      jac*=w[i]*w[j]*rho*r;
      
      nnj (mm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function assembles mass %matrix of axisymmetric rectangular
   finite element with biquadratic approximation functions
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK, 16. 6. 2019
*/
void axisymqq::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));

  mass_matrix (eid,mm);
  
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


void axisymqq::res_mainip_strains (long lcid,long eid)
{
  mainip_strains (lcid,eid,0,0);
}



/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void axisymqq::mainip_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  double xi,eta,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),gp,w,eps(ASTCKVEC(tncomp)),aux;
  ivector nodes(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),tmat;

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
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
  
  
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  reallocv (RSTCKVEC(intordsm[0][0],w));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      geom_matrix (gm,x,y,xi,eta,jac);
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
  }
}



/**
  The function computes strains in nodes of element. Strain values are taken from 
  int. points which are the closest to the particular nodes. Nodal strain values are 
  averaged according to setup.

  @param lcid[in] - load case id
  @param eid[in] - element id
   
  JK 10.5.2002
  Modified by TKo, 6.7.2018 
*/
void axisymqq::nod_strains_ip (long lcid,long eid)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double area;
  vector eps(ASTCKVEC(tncomp)), aux(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain(lcid,ipnum[i],eps);
    
    // compute average strain for midside nodes
    if ((intordsm[0][0]%2 == 0) && (i > 3)) // even number of int. points and midside nodes
    {
      // strain averaging from the two closest int. points
      if (i%2 == 0) // 5-th and 7-th local nodes
        Mm->givestrain(lcid, ipnum[i]+intordsm[0][0], aux);
      else  // 6-th and 8-th local nodes
        Mm->givestrain(lcid, ipnum[i]+1, aux);
      addv(eps, aux, eps);
      cmulv(0.5, eps);
    }

    j=nod[i];
    //  storage of strains to the node
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv (area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes strains at nodes of element, nodal strain values are 
  averaged according to setup.
   
  @param lcid[in] - load case id
  @param eid[in] - element id
   
  JK, 23.9.2004
*/
void axisymqq::nod_strains_comp(long lcid,long eid)
{
  long i, j;
  double jac, area;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  ivector enod(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;

  //  node numbers of element
  Mt->give_elemnodes(eid, enod);
  //  coordinates of element nodes
  Mt->give_node_coord2d(x, y, eid);
  //  nodal displacements
  eldispl(lcid, eid, r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (enod);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe,aux));
    reallocm(RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix(enod,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf(aux,r,tmat);
    copyv(aux,r);
  }
  
  //  natural coordinates of nodes of element
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planeqq (nxi,neta);

  for (i=0;i<nne;i++){
    //  block of geometric matrix
    geom_matrix(gm,x,y,nxi[i],neta[i],jac);
    //  strain computation
    mxv(gm,r,eps);
    
    //  storage of nodal strains
    j=enod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv(area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes strains at nodes of element, no averaging is performed
   
  @param lcid - load case id
  @param eid - element id
   
  TKo, 18.05.2018
*/
void axisymqq::nod_strains(long lcid,long eid)
{
  long i, j, k, m;
  double jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  ivector enod(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;

  //  natural coordinates of nodes of element
  nodcoord_planeqq(nxi,neta);
  //  node numbers of element
  Mt->give_elemnodes (eid,enod);
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  //  nodal displacements
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
  
  for (i=0;i<nne;i++){
    //  block of geometric matrix
    geom_matrix (gm,x,y,nxi[i],neta[i],jac);
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
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqq::res_allip_strains (long lcid,long eid)
{
  allip_strains (lcid,eid,0,0);
}
/**
   function computes strains in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqq::allip_strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid,eid);
}



void axisymqq::strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  vector coord,eps;
  
  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_strains (stra,lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    nod_strains_ip (lcid,eid);
    break;
  }
  case userdefined:{
    /*
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (RSTCKVEC(ncp,eps));
    reallocv (RSTCKVEC(2,coord));
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);
      
      if (Mp->strainaver==0)
	appval (coord[0],coord[1],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);
      
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function axisymqq::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}

/**
   function returns numbers of integration point closest to element nodes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipnum - array of numbers
   
   JK, 25.9.2004
*/
void axisymqq::nodipnum (long eid,ivector &ipnum)
{
  long i,j;
  
  j=intordsm[0][0];
  i=Mt->elements[eid].ipp[0][0];
  
  ipnum[0]=i+j*(j-1)+j-1;
  ipnum[1]=i+j-1;
  ipnum[2]=i;
  ipnum[3]=i+j*(j-1);
  if (j%2) // odd order of integration
  {
    ipnum[4] = (ipnum[0] + ipnum[1])/2;
    ipnum[5] = (ipnum[1] + ipnum[2])/2;
    ipnum[6] = (ipnum[2] + ipnum[3])/2;
    ipnum[7] = (ipnum[3] + ipnum[0])/2;
  }
  else // even order of integration
  {
    ipnum[4] = ipnum[0] - j*j/2;
    ipnum[5] = (ipnum[1] + ipnum[2])/2;
    ipnum[6] = ipnum[3] - j*j/2;
    ipnum[7] = (ipnum[3] + ipnum[0])/2;
  }
  
  /*  ipnum[0]=i+8;
  ipnum[1]=i+2;
  ipnum[2]=i+0;
  ipnum[3]=i+6;
  ipnum[4]=i+5;
  ipnum[5]=i+1;
  ipnum[6]=i+3;
  ipnum[7]=i+7;
  */
}



/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void axisymqq::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval;
  
  k=0;
  reallocv (RSTCKVEC(nne,nodval));
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
}



/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void axisymqq::res_mainip_stresses (long lcid,long eid)
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
void axisymqq::mainip_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      Mm->matstiff (d,ipp);
      
      Mm->givestrain (lcid,ipp,eps);
      
      mxv (d,eps,sig);
      
      Mm->storestress (lcid,ipp,sig);
      
      ipp++;
    }
  }
}



/**
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqq::nod_stresses_ip (long lcid,long eid)
{
  long i, j;
  double area;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(tncomp)), aux(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    
    if ((intordsm[0][0]%2 == 0) && (i > 3)) // even number of int. points and midside nodes
    {
      if (i%2 == 0) // 5-th and 7-th local nodes
        Mm->givestress(lcid, ipnum[i]+intordsm[0][0], aux);
      else  // 6-th and 8-th local nodes
        Mm->givestress(lcid, ipnum[i]+1, aux);
      addv(sig, aux, sig);
      cmulv(0.5, sig);
    }

    //  storage of stresses to the node
    j=nod[i];
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv (area,sig);
      Mt->nodes[j].storestress(lcid,0,area,sig);
    }
  }
  
}



/**
  The function computes stresses at nodes of the eid-th element. Nodal stress values are averaged 
  according to the setup.
   
  @param lcid[in] - load case id
  @param eid - element id
   
  JK, 23.9.2004
  Modified by TKo, 4.7.2018
*/
void axisymqq::nod_stresses_comp (long lcid,long eid)
{
  long i, j;
  double jac, area;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),nxi(ASTCKVEC(nne)),neta(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp)), aux;
  ivector ipnum(ASTCKIVEC(nne)),nodes(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;
  matrix d(ASTCKMAT(tncomp,tncomp)), auxd(ASTCKMAT(tncomp,tncomp));
  
  //  natural coordinates of nodes of element
  nodcoord_planeqq(nxi,neta);
  //  node numbers of element
  Mt->give_elemnodes (eid,nodes);
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (RSTCKVEC(ndofe,aux));
    reallocm (RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    //locglobtransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  nodipnum (eid,ipnum);

  for (i=0;i<nne;i++){
    //  block of geometric matrix
    geom_matrix (gm,x,y,nxi[i],neta[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    //  stiffness matrix of the material
    Mm->matstiff (d,ipnum[i]);

    // stiffness averaging for midside nodes
    if ((intordsm[0][0]%2 == 0) && (i > 3)) // even number of int. points and midside nodes
    {
      // compute average stiffness from two closest int. points
      if (i%2 == 0) // 5-th and 7-th local nodes
        Mm->matstiff (auxd,ipnum[i]+intordsm[0][0]);
      else  // 6-th and 8-th local nodes
        Mm->matstiff (auxd,ipnum[i]+1);
      addm(d, auxd, d);
      cmulm(0.5, d);
    }

    //  stress computation
    mxv (d,eps,sig);
    
    //  number of actual node
    j=nodes[i];
    //  storage of nodal strains
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv (area,sig);
      Mt->nodes[j].storestress(lcid,0,area,sig);
    }
  }
}



/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqq::res_allip_stresses (long lcid,long eid)
{
  allip_stresses (lcid,eid,0,0);
}



/**
   function computes stresses in all integration points
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymqq::allip_stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  res_mainip_stresses (lcid,eid);
}




void axisymqq::stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  vector coord,sig;

  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (stre,lcid,eid,ri,ci);
    //mainip_stresses (lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    nod_stresses_ip (lcid,eid);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    /*
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (RSTCKVEC(ncp,sig));
    reallocv (RSTCKVEC(2,coord));
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);
      
      if (Mp->stressaver==0)
	appval (coord[0],coord[1],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function axisymqq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}



/**
  The function computes other components at nodes of the eid-th element. Nodal values are averaged 
  according to the setup.

  @param eid[in] - element id
   
  JK, 10.5.2002
*/
void axisymqq::nod_other_ip (long eid)
{
  long i, j, ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other, aux;
  double area;
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  

  for (i=0;i<nne;i++)
  {
    ncompo = Mm->givencompother(ipnum[i], 0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    
    if ((intordsm[0][0]%2 == 0) && (i > 3)) // even number of int. points and midside nodes
    {
      reallocv(RSTCKVEC(ncompo, aux));
      if (i%2 == 0) // 5-th and 7-th local nodes
        Mm->giveother(ipnum[i]+intordsm[0][0], 0, ncompo, aux.a);
      else  // 6-th and 8-th local nodes
        Mm->giveother(ipnum[i]+1, 0, ncompo, aux.a);
      addv(other, aux, other);
      cmulv(0.5, other);
    }


    j=nod[i];
    //  storage of other values to the node
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother (0, ncompo, other);
    if (Mp->otheraver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv (area, other);
      Mt->nodes[j].storeother(0, ncompo, area, other);
    }    
  }
}



/**
   function computes load matrix of the axisymmetric quadrilateral
   finite element with bilinear approximation functions
   load vector is obtained after premultiplying load matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load matrix

   8.12.2001
*/
void axisymqq::load_matrix (long eid,matrix &lm)
{
  long i,j;
  double jac,xi,eta,r;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      r = approx (xi,eta,x);
      jac*=r*w[i]*w[j];
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
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
void axisymqq::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void axisymqq::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void axisymqq::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void axisymqq::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void axisymqq::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
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
void axisymqq::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void axisymqq::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void axisymqq::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
void axisymqq::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      if (Mp->strcomp==1)
	Mm->computenlstresses (ipp,Mm->ip[ipp]);
      ipp++;
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
void axisymqq::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      if (Mp->strcomp==1)
	Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
}



/**
   function computes local values which will be used for averaging 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void axisymqq::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
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
   
   TKo 7.2008
*/
void axisymqq::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      if (Mp->strcomp==1)
	Mm->computenlstresses (ipp,Mm->ip[ipp]);
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
void axisymqq::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
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
void axisymqq::elem_integration (integratedquant iq,long lcid,long eid,long /*ri*/,long /*ci*/,vector &nv,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,jac,rad;
  vector w,gp;
  vector ipv(ASTCKVEC(tncomp)),contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  nullv (nv);
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  reallocv (RSTCKVEC(intordsm[0][0],w));
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,0,ipv);
      //  strain-displacement (geometric) matrix
      geom_matrix (gm,x,y,xi,eta,jac);
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      rad = approx (xi,eta,x);
      cmulv (rad*jac*w[i]*w[j],contr);
      //  summation
      addv(contr,nv,nv);
      ipp++;
    }
  }
}



/**
   The function integrates arbitrary selected quantity over finite element volume, e.g.
   it performs \int_{\Omega} \mbf{\sigma} d\Omega which results in integrated values that can 
   be used in the homogenization problems.

   
   @param eid - element id
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param iv - integrated values (output)
   
   @return The function returns integrated values in the %vector iv

   Created by TKo, 21.2.2019
*/
void axisymqq::elem_volintegration_quant (long eid, integratedquant iq, long lcid, vector &iv)
{
  long i, j, ipp;
  double xi, eta, jac, rad;
  vector w(ASTCKVEC(intordsm[0][0])), gp(ASTCKVEC(intordsm[0][0]));
  vector ipv(ASTCKVEC(iv.n)), contr(ASTCKVEC(iv.n));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  gauss_points(gp.a, w.a, intordsm[0][0]);
  Mt->give_node_coord2d(x, y, eid);
  nullv(iv);
  ipp = Mt->elements[eid].ipp[0][0];
  
  for (i=0; i<intordsm[0][0]; i++){
    xi=gp[i];
    for (j=0; j<intordsm[0][0]; j++){
      eta=gp[j];
      //  function assembles required quantity at integration point
      Mm->givequantity(iq, lcid, ipp, 0, ipv);
      //  jacobian determinant
      jac_2d(jac, x, y, xi, eta);
      // compute radius at the given int. point
      rad = approx(xi, eta, x);
      //  contribution to the integrated values
      cmulv(rad*jac*w[i]*w[j], ipv, contr);
      //  summation of contributions to the volume integral of the given quantity
      addv(contr, iv, iv);
      ipp++;
    }
  }
}



/**
   function integrates N^T c N over edges
   
   @param edg - edge id (number of edge)
   @param x, y - coordinates of element nodes
   @param intord - order of numerical integration
   @param gp, w  - coordinates and weights of integration points
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK, 7. 7. 2018
*/
void axisymqq::edge_integral (long edg,vector &x,vector &y,long /*intord*/,vector &gp,vector &w,
			      vector &coef,matrix &km)
{
  long i;
  double xi,eta,jac,ipval,r;
  matrix n(ASTCKMAT(napfun,ndofe));
  
  if (edg==0){
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      
      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      //  computation of jacobian      
      jac1d_2d (jac,x,y,xi,edg);
      //  interpolation of c to integration point
      ipval=approx (xi,eta,coef);
      //  radius
      r=approx (xi,eta,x);
      
      jac *= w[i]*ipval*r;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==1){
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];

      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      //  computation of jacobian      
      jac1d_2d (jac,x,y,xi,edg);
      //  interpolation of c to integration point
      ipval=approx (xi,eta,coef);
      //  radius
      r=approx (xi,eta,x);
      
      jac *= w[i]*ipval*r;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==2){
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      
      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      //  computation of jacobian      
      jac1d_2d (jac,x,y,xi,edg);
      //  interpolation of c to integration point
      ipval=approx (xi,eta,coef);
      //  radius
      r=approx (xi,eta,x);
      
      jac *= w[i]*ipval*r;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==3){
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      
      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      //  computation of jacobian      
      jac1d_2d (jac,x,y,xi,edg);
      //  interpolation of c to integration point
      ipval=approx (xi,eta,coef);
      //  radius
      r=approx (xi,eta,x);
      
      jac *= w[i]*ipval*r;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}

/**
   function picks up nodal values on required edges
   
   @param edg - number of required edge
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all edges
   
   JK, 7. 7. 2018
*/
void axisymqq::edgenodeval (long edg,vector &nodval,double *list)
{
  long i,j,k,l;
  ivector edgenod(ASTCKIVEC(nned));
  
  nullv (nodval);
  //  nodes on required edge
  quadquadrilat_edgnod (edgenod.a,edg);
  
  k=0;
  for (i=0;i<nned;i++){
    l=edgenod[i]*napfun;
    for (j=0;j<napfun;j++){
      nodval[l+j]=list[edg*nned*napfun+k];
      k++;
    }
  }
}



/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - indicators of loaded edges
   @param list - list of prescribed nodal values
   @param nf - %vector of nodal forces caused by edge load
   
   JK, 7. 7. 2018
*/
void axisymqq::edgeload (long eid,long *le,double *list,vector &nf)
{
  long i;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),nodval(ASTCKVEC(ndofe)),coeff(ASTCKVEC(nne)),av(ASTCKVEC(ndofe));
  matrix km(ASTCKMAT(ndofe,ndofe));
  
  //  coefficient c in N^T c N
  for (i=0;i<nne;i++){
    coeff[i]=1.0;
  }
  //  coordinates of element nodes
  Mt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordb);
  
  //  cleaning of nf
  nullv (nf);
  
  for (i=0;i<ned;i++){
    if (le[i]>0){
      nullm (km);
      //  selection of appropriate nodal values
      edgenodeval (i,nodval,list);
      //  integration over element edge
      edge_integral (i,x,y,intordb,gp,w,coeff,km);
      
      mxv (km,nodval,av);
      
      addv (nf,av,nf);
    }
  }

}

/*
void axisymqq::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  long i;
  double ww,jac,xi,eta;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),av(ASTCKVEC(ndofe)),v(ASTCKVEC(ndofe));
  matrix n(ASTCKMAT(napfun,ndofe)),am(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordb);

  if (le[0]==1){
    fillm (0.0,am);
    eta = 1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,0);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[0]=nv[0];  av[1]=nv[1];  av[2]=nv[2];  av[3]=nv[3];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    fillm (0.0,am);
    xi = -1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];

      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,1);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[2]=nv[4];  av[3]=nv[5];  av[4]=nv[6];  av[5]=nv[7];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    fillm (0.0,am);
    eta = -1.0;
    for (i=0;i<intordb;i++){
      xi = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,xi,2);
      jac *= ww;
      
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[4]=nv[8];  av[5]=nv[9];  av[6]=nv[10];  av[7]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    fillm (0.0,am);
    xi = 1.0;
    for (i=0;i<intordb;i++){
      eta = gp[i];
      ww = w[i];
      
      bf_matrix (n,xi,eta);
      
      jac1d_2d (jac,x,y,eta,3);
      jac *= ww;

      nnj (am.a,n.a,jac,n.m,n.n);
    }
    fillv (0.0,av);
    av[6]=nv[12];  av[7]=nv[13];  av[0]=nv[14];  av[1]=nv[15];
    mxv (am,av,v);  addv (nf,v,nf);
  }
}
*/


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
void axisymqq::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
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
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param ncoord - %vector with natural coordinates of integration point (ouput)
   
   @return The function returns natural coordinates in the argument ncoord.

   Created by TKo, 12.2016
*/
void axisymqq::ipncoord (long eid, long ipp, vector &ncoord)
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
void axisymqq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
    for(i = 0; i < nne; i++)
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
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   21.6.2004, JK
*/
void axisymqq::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k;
  double xi,eta;
  vector w,gp;
  
  k=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  
	  ipval[k]=approx (xi,eta,nodval);
	  k++;
	}
      }
    }
  }
}

/**
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   21.6.2004, JK
*/
void axisymqq::intpointval2 (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k;
  double xi,eta;
  vector w,gp;
  vector modnodval(ASTCKVEC(Asymlq->nne));
  
  for (i=0;i<Asymlq->nne;i++){
    modnodval[i]=nodval[i];
  }
  
  k=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  
	  ipval[k]=Asymlq->approx (xi,eta,modnodval);
	  k++;
	}
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
   
   01/11/2016 TKr according to axisymlq
*/
void axisymqq::mechq_nodval (long eid,vector &nodval,nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,ipnum);

  for (i=0;i<nne;i++){
    //copy nonmechanical quantity from closest int. point
    nodval[i] = Mm->givemechq(qt, ipnum[i]);
  }
}

/**
   Function computes required mechanical quantity at nodes of element.
   linear approximation functions are used

   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything.
   
   01/11/2016 TKr according to axisymlq
*/
void axisymqq::mechq_nodval2 (long eid,vector &nodval,nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,ipnum);

  for (i=0;i<Asymlq->nne;i++){
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
   
   01/11/2016 TKr accordintg to Tomas Koudelka and axisymlq
*/
void axisymqq::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
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
  nodipnum (eid,ipnum);

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

  // restore original integration point content of strain/stress/other/eqother/nonloc arrays
  Mm->storestrain(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.strain);
  Mm->storestress(0, ipid, 0, Mb->nlc*ipb.ncompstr, ipb.stress);
  Mm->storeother(ipid, 0, ipb.ncompother, ipb.other);
  Mm->storeeqother(ipid, 0, ipb.ncompeqother, ipb.eqother);
  Mm->storenonloc(ipid, 0, ncompnl, ipb.nonloc); // maybe this would not be necessary
}
