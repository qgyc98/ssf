#include "plelemqq.h"
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
#include "gadaptivity.h"
#include "node.h"
#include "intpoints.h"
#include "element.h"
#include "loadcase.h"
#include <stdlib.h>
#include <math.h>



planeelemqq::planeelemqq (void)
{
  long i,j;
  
  //  number nodes on element
  nne=8;
  //  number of DOFs on element
  ndofe=16;
  //  number of strain/stress components
  tncomp=3;
  //  number of components for graphic purposes
  gncomp=4;
  //  number of functions approximated
  napfun=2;
  //  order of numerical integration of mass matrix
  intordmm=4;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;
  
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
  
  nip[0][0]=9;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=3;
}

planeelemqq::~planeelemqq (void)
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
   procedure approximates function defined by nodal values
   
   @param xi,eta - natural coordinates on element
   @param nodval - nodal values
   
   JK
*/
double planeelemqq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_quad_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function returns %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param xi,eta - natural coordinates
   
   JK, 25.8.2001
*/
void planeelemqq::bf_matrix (matrix &n,double xi,double eta)
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
   function assembles geometric %matrix (strain-displacement %matrix)
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian
   
   JK, 9.7.2001
*/
void planeelemqq::geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i,i1,i2;
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  
  dx_bf_quad_4_2d (dx.a,xi,eta);
  dy_bf_quad_4_2d (dy.a,xi,eta);
  
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
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - element nodes
   @param tmat - transformation %matrix
   
   JK,
*/
void planeelemqq::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of plane stress rectangular
   finite element with biquadratic approximation functions
   
   this function is used in plane stress/strain elements (function is called
   by function res_stiffness_matrix) and shell elements

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - node coordinates
   
   JK, 25.8.2001
*/
void planeelemqq::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,j,ii,jj,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp)),gm(ASTCKMAT(tncomp,ndofe));
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  
  nullm (sm);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  
	  //  geometric matrix
	  geom_matrix (gm,x,y,xi,eta,jac);
	    
	  if(jac < 0.0){
	    //det = fabs(det);
	    print_err("wrong numbering of nodes on 2D element number %ld, negative jacobian jac = %e", __FILE__, __LINE__, __func__, eid+1,jac);
	    abort();
	  }

	  //  matrix of material stiffness
	  Mm->matstiff (d,ipp);
	  
	  thick = approx (xi,eta,t);
	  
	  jac*=thick*w[i]*w[j];
	  
	  //  contribution to the stiffness matrix of the element
	  bdbj (sm.a,gm.a,d.a,jac,gm.m,gm.n);
	  
	  ipp++;
	}
      }
    }
  }
}

/**
   function assembles stiffness %matrix of plane stress rectangular
   finite element with biquadratic approximation functions
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK
*/
void planeelemqq::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  matrix tmat (ASTCKMAT(ndofe,ndofe));
  
  Mt->give_node_coord2d (x,y,eid);
  
  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes mass %matrix of the plane stress rectangular
   finite element with biquadratic approximation functions
   
   this function is used in plane stress/strain elements (function is called
   by function res_mass_matrix) and shell elements

   @param eid - number of element
   @param mm - mass %matrix
   @param x,y - node coordinates
   
   JK, 25.8.2001
*/
void planeelemqq::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick,rho;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne)),dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  //  density of material (defined at nodes or on element)
  Mc->give_density (eid,nodes,dens);

  gauss_points (gp.a,w.a,intordmm);
  
  nullm (mm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      
      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      
      //  thickness at integration point
      thick = approx (xi,eta,t);
      //  density at integration point
      rho = approx (xi,eta,dens);
      jac*=w1*w2*thick*rho;
      
      nnj (mm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function assembles mass %matrix of plane stress rectangular
   finite element with biquadratic approximation functions
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK
*/
void planeelemqq::res_mass_matrix (long eid,matrix &mm)
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
   finite element with biquadratic approximation functions
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   this function is used in plane stress/strain elements (function is called
   by function res_load_matrix) and shell elements

   @param eid - number of element
   @param lm - load %matrix
   @param x,y - node coordinates

   JK, 25.8.2001
*/
void planeelemqq::load_matrix (long eid,matrix &lm,vector &x,vector &y)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordmm)),gp(ASTCKVEC(intordmm)),t(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);

  gauss_points (gp.a,w.a,intordmm);
  
  nullm (lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];

      jac_2d (jac,x,y,xi,eta);

      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      
      //  thickness at integration point
      thick = approx (xi,eta,t);
      jac*=w1*w2*thick;
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

/**
   function assembles load %matrix of plane stress rectangular
   finite element with biquadratic approximation functions
   
   @param eid - element id
   @param lm - load %matrix
   
   JK
*/
void planeelemqq::res_load_matrix (long eid,matrix &lm)
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
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (lm,tmat);
  }

}

/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, modified 23.11.2006
*/
void planeelemqq::res_ip_strains (long lcid,long eid)
{
  long transf;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),aux(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat(ASTCKMAT(ndofe,ndofe));
  
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  ip_strains (lcid,eid,0,0,x,y,r);
}

/**
   function computes strains at integration points of element
   
   this function is used in plane stress/strain elements (function is called
   by function res_ip_strains) and shell elements

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y - node coordinates
   @param r - nodal displacements
   
   10.5.2002, JK, modified 23.11.2006
*/
void planeelemqq::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ipp;
  double xi,eta,jac;
  vector gp,w,eps;
  matrix gm;

  reallocv (RSTCKVEC(intordsm[0][0],gp));
  reallocv (RSTCKVEC(intordsm[0][0],w));
  reallocv (RSTCKVEC(ncomp[0],eps));
  reallocm (RSTCKMAT(ncomp[0],ndofe,gm));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      
      //  geometric matrix (strain-displacement matrix)
      geom_matrix (gm,x,y,xi,eta,jac);
      //  strain computation
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
  }
}



/**
  The function assembles strains at nodes of element.
  Strain values are obtained from the nearest integration points.

  @param lcid[in] - load case id
  @param eid[in] - element id
  @param ri,ci[in] - row and column indices (default value is 0, nonzero values are used in shell elements)

  10.5.2002
*/
void planeelemqq::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  double area;
  vector eps(ASTCKVEC(gncomp)), aux(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,ri,ci,ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain (lcid,ipnum[i],eps);
    
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
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv (area,eps);
      Mt->nodes[j].storestrain(lcid,0,area,eps);
    }
  }
}



/**
  The function computes nodal strains directly, averageing of nodal strains is performed according to setup.
   
  @param lcid - load case id
  @param eid - element id
   
  01/11/2016 by TKr according to JK
*/
void planeelemqq::nod_strains_comp (long lcid,long eid)
{
  long i, j;
  double jac, area;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)),aux;
  ivector enod(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;
  
  //  node coordinates
  Mt->give_node_coord2d (x, y, eid);
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
    //locglobtransf (aux, r, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planeqq(nxi, neta);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix(gm, x, y, nxi[i], neta[i], jac);
    //  strain computation
    mxv(gm, r, eps);
    
    //  storage of strains to the node
    j=enod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv(area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes nodal strains directly with no averageing.
   
  @param lcid - load case id
  @param eid - element id
   
  23/5/2018 by TKo
*/
void planeelemqq::nod_strains(long lcid,long eid)
{
  long i,j,k,m;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),nxi(ASTCKVEC(nne)),neta(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),eps(ASTCKVEC(tncomp)),aux;
  ivector enod(ASTCKIVEC(nne));
  matrix tmat,gm(ASTCKMAT(tncomp,ndofe));
  
  //  node coordinates
  Mt->give_node_coord2d (x, y, eid);
  //  node numbers
  Mt->give_elemnodes (eid, enod);
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
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planeqq (nxi,neta);
  
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
   function computes strains at required position on elements
   function possibly transforms strains to required coordinate systems
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 23.11.2006
*/
void planeelemqq::strains (long lcid,long eid,long ri,long ci)
{
  vector coord,eps;
  
  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //res_ip_strains (lcid,eid,ri,ci);
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
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemqq::strains (file %s, line %d).\n",__FILE__,__LINE__);
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
void planeelemqq::nodipnum (long eid,long ri,long ci,ivector &ipnum)
{
  // long j=intordsm[0][0];
  long i=Mt->elements[eid].ipp[ri][ci];

  ipnum[0]=i+8;
  ipnum[1]=i+2;
  ipnum[2]=i+0;
  ipnum[3]=i+6;
  ipnum[4]=i+5;
  ipnum[5]=i+1;
  ipnum[6]=i+3;
  ipnum[7]=i+7;

}

/**
   function computes stresses at integration points
   
   
*/
void planeelemqq::res_ip_stresses (long /*lcid*/,long /*eid*/)
{
  //mainip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses in integration points of element
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void planeelemqq::ip_stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
  /*
  long i,j,jj,ipp;
  double xi,eta;
  vector gp,w,eps,epst,epstt,sig,auxsig;
  matrix d(ASTCKMAT(tncomp,tncomp));

  long ii=0;

  reallocv (RSTCKVEC(ncomp[ii],sig));
  reallocv (RSTCKVEC(ncomp[ii],auxsig));
  reallocv (RSTCKVEC(intordsm[ii][ii],gp));
  reallocv (RSTCKVEC(intordsm[ii][ii],w));
  
  gauss_points (gp.a,w.a,intordsm[ii][ii]);
  ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
  
  for (i=0;i<intordsm[ii][ii];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ii][ii];j++){
      eta=gp[j];
      
      Mm->matstiff (d,ipp);
      
      fillv (0.0,sig);
      for (jj=0;jj<nb;jj++){
	reallocv (RSTCKVEC(ncomp[jj],eps));
	reallocm (RSTCKMAT(ncomp[ii],ncomp[jj],dd));
	
	Mm->givestrain (lcid,ipp,cncomp[jj],ncomp[jj],eps);
	
	dmatblock (ii,jj,d,dd);
	mxv (dd,eps,auxsig);
	addv (auxsig,sig,sig);
      }
      
      Mm->storestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
      
      ipp++;
    }
  }
  
  */
}


/**
   function computes stresses at nodes of element

   @param lcid[in] - load case id
   @param eid[in] - element id
   @param ri,ci[in] - row and column indices
   
   10.5.2002
*/
void planeelemqq::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i, j;
  double area;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(gncomp)), aux(ASTCKVEC(gncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ri, ci, ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress(lcid, ipnum[i], sig);
    
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
      Mt->nodes[j].storestress(lcid, 0, sig);
    if (Mp->stressaver==2)
    {
      area = Mt->give_area(eid)*0.125;
      cmulv(area,sig);
      Mt->nodes[j].storestress(lcid, 0, area, sig);
    }
  }
}


/**
  The function computes nodal stresses directly. Nodal stress values are averaged according to the setup.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id
  @param ri,ci[in] - row and column indices
  @param stra[in] - array for strain components (stra[i] = pointer to the array of strain components at the i-th node)
  @param stre[out] - array for computed nodal stresses (stre[i] = pointer to the array of stress components at the i-th node)
   
  01/11/2016 by TKr according to JK and planeelemlq
  Modified 06/07/2018 by TKo
*/
void planeelemqq::nod_stresses_comp (long /*lcid*/, long eid, long ri, long ci, double **stra, double **stre)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne));
  vector eps(ASTCKVEC(tncomp)), sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp, tncomp)), auxd(ASTCKMAT(tncomp, tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ri, ci, ipnum);

  //  loop over nodes
  for (i=0;i<nne;i++){
    for (j=0;j<eps.n;j++)
      eps[j]=stra[i][j];

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

    mxv (d,eps,sig);

    for (j=0;j<eps.n;j++)
      stre[i][j]=sig[j];
  }
}


void planeelemqq::stresses (long lcid,long eid,long ri,long ci)
{
  vector coord,sig;
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (stre,lcid,eid,ri,ci);
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
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemqq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}



/**
  The function computes other values in nodes at the eid-th element.

  @param eid[in] - element id
  @param ri,ci[in] - row and column indices (default value is 0, nonzero values are used in shell elements)
   
  10.5.2002
*/
void planeelemqq::nod_other_ip (long eid,long ri,long ci)
{
  long i, j, ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other, aux;
  double area;
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,ri,ci,ipnum);
  
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
      cmulv(area, other);
      Mt->nodes[j].storeother(0, ncompo, area, other);
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
   
   25.8.2001, JK modified 23.11.2006
*/
void planeelemqq::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemqq::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void planeelemqq::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=stressincr;
  
  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   this function is used in plane stress/strain elements (function is called
   by function res_eigstrain_forces) and shell elements

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   30.11.2002, JK
*/
void planeelemqq::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - internal forces
   
   JK, modified 23.11.2006
*/
void planeelemqq::res_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemqq::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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
void planeelemqq::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - array containing nodal forces
   
   30.11.2002, JK
*/
void planeelemqq::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemqq::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      //  computation of correct stresses
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
void planeelemqq::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      //  computation of correct stresses
      if (Mp->strcomp==1)
	Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void planeelemqq::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      //  computation of correct stresses
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
   
   JK, 27.11.2006
*/
void planeelemqq::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      
      //  computation of correct stresses
      if (Mp->strcomp==1)
	Mm->compnonloc_nlstresses (ipp);
      
      ipp++;
    }
  }
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void planeelemqq::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
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
   
   JK, 27.11.2006
   TKo 7.2008
*/
void planeelemqq::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp,t(ASTCKVEC(nne)),ipv(ASTCKVEC(tncomp)),contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  Mc->give_thickness (eid,nodes,t);
  
  nullv (nv);
  
  reallocv (RSTCKVEC(intordsm[0][0],gp));
  reallocv (RSTCKVEC(intordsm[0][0],w));
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      thick = approx (xi,eta,t);
      
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);

      //  strain-displacement (geometric) matrix
      geom_matrix (gm,x,y,xi,eta,jac);
      
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      cmulv (jac*w[i]*w[j]*thick,contr);
      
      //  summation
      addv(contr,nv,nv);
      
      ipp++;
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

   10.1.2002
*/
void planeelemqq::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
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
void planeelemqq::ipncoord (long eid,long ipp,vector &ncoord)
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
   
   8.5.2002
*/
void planeelemqq::ipcoordblock (long eid,long ri,long ci,double **coord)
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
      coord[k++][2]=0.0;
    }
  }
}

/**
   function computes nodal forces caused by edge load
   
   @param eid - element id
   @param le - list of loaded edges
   @param nv - nodal values of edge loads, it contains 24 components
   @param nf - %vector of nodal forces

   nv contains 2 values at each node on each edge, each edge contains 3 nodes,
   element contains 4 edges, therefore, nv contains 4*3*2 = 24 components

   JK, modification 23.11.2006
*/
void planeelemqq::res_nodeforces (long eid,long *le,double *nv,vector &nf)
{
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  Mt->give_node_coord2d (x,y,eid);
  nodeforces (eid,le,nv,nf,x,y);
}

/**
   function computes nodal forces caused by edge load
   
   this function is used in plane stress/strain elements (function is called
   by function res_nodeforces) and shell elements

   @param eid - element id
   @param le - list of loaded edges
   @param nv - nodal values of edge loads, it contains 24 components
   @param nf - %vector of nodal forces
   @param x,y - node coordinates

   nv contains 2 values at each node on each edge, each edge contains 3 nodes,
   element contains 4 edges, therefore, nv contains 4*3*2 = 24 components

   JK, modification 23.11.2006
*/
void planeelemqq::nodeforces (long /*eid*/,long *le,double *nv,vector &nf,vector &x,vector &y)
{
  long i;
  double ww,jac,xi,eta;
  vector gp(ASTCKVEC(intordb)),w(ASTCKVEC(intordb)),av(ASTCKVEC(ndofe)),v(ASTCKVEC(ndofe));
  matrix n(ASTCKMAT(napfun,ndofe)),am(ASTCKMAT(ndofe,ndofe));
  

  gauss_points (gp.a,w.a,intordb);

  if (le[0]==1){
    //  first edge is loaded
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
    av[0]=nv[0];  av[1]=nv[1];  av[8]=nv[2];  av[9]=nv[3];  av[2]=nv[4];  av[3]=nv[5];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[1]==1){
    //  second edge is loaded
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
    av[2]=nv[6];  av[3]=nv[7];  av[10]=nv[8];  av[11]=nv[9];  av[4]=nv[10];  av[5]=nv[11];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[2]==1){
    //  third edge is loaded
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
    av[4]=nv[12];  av[5]=nv[13];  av[12]=nv[14];  av[13]=nv[15];  av[6]=nv[16];  av[7]=nv[17];
    mxv (am,av,v);  addv (nf,v,nf);
  }
  if (le[3]==1){
    //  fourth edge is loaded
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
    av[6]=nv[18];  av[7]=nv[19];  av[14]=nv[20];  av[15]=nv[21];  av[0]=nv[22];  av[1]=nv[23];
    mxv (am,av,v);  addv (nf,v,nf);
  }
}


/**
   function computes displacement, strain and stress in the centre=middle of element
   
   @param eid - number of element
   @paran dd  - decomposed deformation
   @param valel - array of displacements, strains and stresses in the middle of all elements
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void planeelemqq::midpoints (long eid,long dd,double *valel)
{
  long i,j,lcid;
  double r[2];
  vector nodval(ASTCKVEC(nne)),der(ASTCKVEC(tncomp));
  ivector nodes(ASTCKIVEC(nne));
  
  lcid = 0;
  Mt->give_elemnodes (eid,nodes);  
  
  for (i=0;i<2;i++){
    for (j=0;j<nne;j++){
      noddispl (lcid,r,nodes[j]);
      nodval[j] = r[i];
    }
    valel[i] = approx (0.0,0.0,nodval);
    if (dd) valel[i+2] = valel[i];
  }

  //appstrain (lcid,eid,0.0,0.0,0,tncomp,der);
  for (i=0;i<tncomp;i++)
    valel[i+2*(1+dd)] = der[i];
  
  //appstress (lcid,eid,0.0,0.0,0,tncomp,der);
  for (i=0;i<tncomp;i++)
    valel[i+2*(1+dd)+tncomp] = der[i];
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
void planeelemqq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
   
   2.3.2004, JK
   07.2008 TKo - multiplictaion by thickness added
*/
void planeelemqq::ipvolume (long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),w,gp;
  
  //  element nodes
  Mt->give_elemnodes (eid,nodes);
  //  thickness of the element
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (RSTCKVEC(intordsm[ii][jj],w));
      reallocv (RSTCKVEC(intordsm[ii][jj],gp));

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
}

/**
   function interpolates the nodal values to the integration points on the element
   quadratic approximation functions are used

   @param eid - element id
   @param nodval - nodal values
   @param ipval - value at integration points
   
   21.6.2004, JK
*/
void planeelemqq::intpointval (long /*eid*/,vector &nodval,vector &ipval)
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
   linear approximation functions are used
   
   @param eid - element id
   @param nodval - nodal values
   @param ipval - value at integration points
   
   21.6.2004, JK
*/
void planeelemqq::intpointval2 (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k;
  double xi,eta;
  vector w,gp;
  vector modnodval(ASTCKVEC(Pelq->nne));
  
  for (i=0;i<Pelq->nne;i++){
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
	  
	  ipval[k]=Pelq->approx (xi,eta,modnodval);
	  k++;
	}
      }
    }
  }
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
void planeelemqq::ntdbr_vector (long eid,vector &ntdbr)
{
  long intord = 3;
  long i,j,k,l,ipp,ri,ci,lcid;
  double thick,xi,eta,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),gp(ASTCKVEC(intord)),w(ASTCKVEC(intord)),t(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp)),bf(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  ri = ci = lcid = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  eldispl (lcid,eid,r.a);

  Mm->matstiff (d,ipp);

  gauss_points (gp.a,w.a,intord);

  nullv (ntdbr);

  for (i=0;i<intord;i++){
    xi=gp[i];
    for (j=0;j<intord;j++){
      eta=gp[j];

      geom_matrix (gm,x,y,xi,eta,jac);
      mxv (gm,r,eps);
      mxv (d,eps,sig);

      bf_quad_4_2d (bf.a,xi,eta);
      
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
void planeelemqq::ntn_matrix (long eid,matrix &ntn)
{
  long i,j,k,l,intord = 3;
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
      
      bf_quad_4_2d (bf.a,xi,eta);

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
double planeelemqq :: compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull)
{
  long intord = 3;
  long i,j,ipp,ri,ci,lcid;
  double area,thick,xi,eta,jac,contr;
  double zero=1.0e-20;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),gp(ASTCKVEC(intord)),w(ASTCKVEC(intord)),bf(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  vector sig_star(ASTCKVEC(tncomp)),sig_roof(ASTCKVEC(tncomp)),sig_err(ASTCKVEC(tncomp)),eps(ASTCKVEC(tncomp));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp)),dinv(ASTCKMAT(tncomp,tncomp));

  ri = ci = lcid = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  Mc->give_thickness (eid,nodes,t);
  eldispl (lcid,eid,r.a);

  Mm->matstiff (d,ipp);
  invm (d,dinv,zero);

  gauss_points (gp.a,w.a,intord);

  e2 = u2 = 0;
  for (i=0;i<intord;i++){
    xi=gp[i];
    for (j=0;j<intord;j++){
      eta=gp[j];

      bf_quad_4_2d (bf.a,xi,eta);
      give_der_star (bf,rsigfull,nodes,sig_star,Mt->nn);
      
      geom_matrix (gm,x,y,xi,eta,jac);
      mxv (gm,r,eps);
      mxv (d,eps,sig_roof);

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
void planeelemqq :: elchar (long eid, matrix &spsig)
{
  long intord = 2;
  long i,j,ipp;
  vector eps(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  reallocm(RSTCKMAT(4, tncomp, spsig));
  
  ipp = Mt->elements[eid].ipp[0][0];
  Mm->matstiff (d,ipp);
  
  //body = Ada->body[7];
  
  //if (body){
    for (i=0; i<intord; i++)
      for (j=0;j<intord;j++){
	eps[0] = Mm->ip[ipp  ].strain[0];
	eps[1] = Mm->ip[ipp  ].strain[1];
	eps[2] = Mm->ip[ipp+4].strain[2];
	ipp++;
	mxv (d.a, eps.a, spsig.a + (i*intord+j)*tncomp, d.m, d.n);
      }
  // }
  // else{
  //   long lcid = 0;
  //   double jac;
  //   vector gp(intord),w(intord),x(nne),y(nne),r(ndofe);
  //   matrix gm(tncomp,ndofe);
  //   
  //   Mt->give_node_coord2d (x,y,eid);
  //   eldispl (lcid,eid,r.a);
  //   gauss_points (gp.a,w.a,intord);
  //   
  //   for (i=0;i<intord;i++)
  //     for (j=0;j<intord;j++){
  // 	geom_matrix (gm,x,y,gp[i],gp[j],jac);
  // 	mxv (gm,r,eps);
  // 	sig = spsig + (i*intord+j)*tncomp;
  // 	mxv (d.a,eps.a,sig,d.m,d.n);
  //     }
  // }
}
////////////////////       /* termitovo */       ////////////////////////////////////



/**
   Function computes required mechanical quantity at nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything.
   
   01/11/2016 TKr accordintg to Tomas Koudelka and planeelemlq
*/
void planeelemqq::mechq_nodval (long eid,vector &nodval,nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,0,0,ipnum);

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
   
   01/11/2016 TKr
*/
void planeelemqq::mechq_nodval2 (long eid,vector &nodval,nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,0,0,ipnum);

  for (i=0;i<Pelq->nne;i++){
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
   
   01/11/2016 TKr accordintg to Tomas Koudelka and planeelemlq
*/
void planeelemqq::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
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
  nodipnum (eid,0,0,ipnum);

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
