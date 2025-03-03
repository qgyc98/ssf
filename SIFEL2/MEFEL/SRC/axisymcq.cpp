#include "axisymcq.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "mechtop.h"
#include "axisymlq.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "loadcase.h"
#include "intpoints.h"
#include <math.h>
#include <stdlib.h>


axisymcq::axisymcq(void)
{
  long i,j;
  
  nne=12;  ndofe=24;  tncomp=4;  napfun=2;  ned=4;  nned=4;
  intordmm=3;  intordb=3;
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
  nip[0][0]=9;
  
  intordsm[0][0]=3;
  
  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

}

axisymcq::~axisymcq(void)
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
  The function approximates function defined by nodal values at the given point on element.
   
  @param xi[in], eta[in] - natural coordinates of point on element
  @param nodval[out]     - nodal values
   
  @retval The function returns approximated value from the given nodal values.

  Created by Tomas Koudelka, 19.3.2019
*/
double axisymcq::approx(double xi, double eta, vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_cubic_4_2d(bf.a, xi, eta);
  scprd(bf, nodval, f);

  return f;
}



/**
  The function returns matrix of approximation functions values evaluated 
  at the given point on element.
   
  @param n[out]          - matrix of approximation functions values
  @param xi[in], eta[in] - natural coordinates of the given point on element
   
  Created by Tomas Koudelka, 19.3.2019
*/
void axisymcq::bf_matrix(matrix &n, double xi, double eta)
{
  long i, j, k;
  vector bf(ASTCKVEC(nne));
  
  nullm(n);

  bf_cubic_4_2d(bf.a, xi, eta);
  
  j=0;  k=1;
  for (i=0; i<nne; i++){
    n[0][j]=bf[i];
    n[1][k]=bf[i];
    j+=2;  k+=2;
  }
}



/**
  The function assembles geometric %matrix at the given point on element.
  The order of starins is being assumed as follows: 
  epsilon_x = du/dx
  epsilon_y = dv/dy
  epsilon_fi = u/r
  epsilon_xy = du/dy + dv/dx
   
  @param gm[out]         - geometric matrix
  @param ri[in]          - block index
  @param x[in], y[in]    - arrays of nodal coordinates
  @param xi[in], eta[in] - natural coordinates of the given point on element
  @param jac[in]         - jacobian
   
  @return The function returns the resulting %matrix in the argument gm.

  Created by Tomas Koudelka, 19.3.2019
*/
void axisymcq::geom_matrix(matrix &gm, vector &x, vector &y, double xi, double eta, double &jac)
{
  long i, i1, i2;
  double r;
  vector bf(ASTCKVEC(nne)), dx(ASTCKVEC(nne)), dy(ASTCKVEC(nne));
  
  // compute derivatives with respect to natural coordinates and store them 
  // to the vectors dx and dy temporarily
  dksi_bf_cubic_4_2d(dx.a, xi, eta);
  deta_bf_cubic_4_2d(dy.a, xi, eta);
  // compute base function matrix
  bf_cubic_4_2d(bf.a, xi, eta);

  // compute derivatives with respect to x and y and jacobian of the transformation
  derivatives_2d(dx, dy, jac, x, y, xi, eta);

  // compute radius at the given point
  r = approx (xi, eta, x);
  if (fabs(r)<Mp->zero){
    //fprintf (stderr,"\n\n radius is equal %e in function axisymcq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
    r=0.00001;
  }
  
  nullm(gm);
  
  // assemble geometric matrix from the derivatives df/dx and df/dy
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
  The function assembles part of geometric matrix at the given point on element.
  The order of starins is being assumed as follows: 
  epsilon_x = du/dx
  epsilon_y = dv/dy
  epsilon_fi = u/r
  epsilon_xy = du/dy + dv/dx
   
  @param gm[out]         - geometric matrix
  @param ri[in]          - block index
  @param x[in], y[in]    - arrays of node coordinates
  @param xi[in], eta[in] - natural coordinates of the given point on element
  @param jac[out]        - jacobian
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::geom_matrix_block(matrix &gm, long ri, vector &x, vector &y,
                                 double xi, double eta, double &jac)
{
  if (nb==1){
    geom_matrix(gm, x, y, xi, eta, jac);
  }
  else{
    long i, i1, i2;
    double r;
    vector bf(ASTCKVEC(nne)), dx(ASTCKVEC(nne)), dy(ASTCKVEC(nne));
    
    // compute derivatives with respect to natural coordinates and store them 
    // to the vectors dx and dy temporarily
    dksi_bf_cubic_4_2d(dx.a, xi, eta);
    deta_bf_cubic_4_2d(dy.a, xi, eta);
    // compute base function matrix
    bf_cubic_4_2d(bf.a, xi, eta);
    
    // compute derivatives with respect to x and y and jacobian of the transformation
    derivatives_2d(dx, dy, jac, x, y, xi, eta);
    
    // compute radius at the given point
    r = approx (xi, eta, x);
    if (fabs(r)<Mp->zero){
      //fprintf (stderr,"\n\n radius is equal %e in function axisymcq::geom_matrix_block (%s, line %d)",r,__FILE__,__LINE__);
      r=0.00001;
    }
    
    nullm(gm);

    // assemble required block of geometric matrix    
    if (ri==0){
      i1=0;  i2=1;
      for (i=0; i<nne; i++){
	gm[0][i1]=dx[i];
	gm[1][i2]=dy[i];
	i1+=2;  i2+=2;
      }
    }
    if (ri==1){
      i1=0;
      for (i=0; i<nne; i++){
	gm[0][i1]=bf[i]/r;
	i1+=2;
      }
    }
    if (ri==2){
      i1=0;  i2=1;
      for (i=0; i<nne; i++){
	gm[0][i1]=dy[i];
	gm[0][i2]=dx[i];
	i1+=2;  i2+=2;
      }
    }
  }
}



/**
  nutno otestovat! pak je mozne smazat tuto hlasku
   
  The function assembles transformation matrix x_g = T x_l on the element
  from local coordinate systems of the element nodes to the global coordinate system.

  @param nodes[in] - array of element node id
  @param tmat[out] - resulting transformation %matrix

  @return The function returns transformation matrix in the argument tmat.

  Created by Tomas Koudelka, 19.3.2019
*/
void axisymcq::transf_matrix(ivector &nodes, matrix &tmat)
{
  long i, n;

  identm(tmat);

  n=nodes.n;

  for (i=0; i<n; i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*2][i*2]   = Mt->nodes[nodes[i]].e1[0];    tmat[i*2][i*2+1]   = Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2] = Mt->nodes[nodes[i]].e1[1];    tmat[i*2+1][i*2+1] = Mt->nodes[nodes[i]].e2[1];
    }
  }
}



/**
  The function computes stiffness %matrix of axisymmetric quadrilateral
  finite element with bicubic approximation functions.
   
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices
  @param sm[out]        - stiffness %matrix

  @return The function returns resulting stiffness %matrix in the argumnet sm.
  
  Created by Tomas Koudelka, 19.3.2019
*/
void axisymcq::stiffness_matrix(long eid, long ri, long ci, matrix &sm)
{
  long i, j;
  long ipp;
  double xi, eta, jac, r;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp;
  matrix gm(ASTCKMAT(tncomp,ndofe)), d(ASTCKMAT(tncomp,tncomp));
  
  Mt->give_node_coord2d(x, y, eid);
  
  nullm(sm);
  reallocv(RSTCKVEC(intordsm[0][0],w));
  reallocv(RSTCKVEC(intordsm[0][0],gp));
  
  gauss_points(gp.a, w.a, intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0; i<intordsm[0][0]; i++){
    xi=gp[i];
    for (j=0; j<intordsm[0][0]; j++){
      eta=gp[j];
      
      //  geometric matrix
      geom_matrix (gm, x, y, xi, eta, jac);
      if(jac < 0.0){
	//jac = fabs (jac);
	print_err("wrong numbering of nodes on element number %ld, negative jacobian! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	abort();
      }
  
      //  matrix of stiffness of the material
      Mm->matstiff (d,ipp);
      // compute radius at the given point
      r = approx(xi, eta, x);
      jac *= w[i]*w[j]*r;
      //  contribution to the stiffness matrix of the element
      bdbjac(sm, gm, d, gm, jac);
      
      ipp++;
    }
  }
}



/**
  The function computes resulting stiffness %matrix of element with respect to
  local coordinate systems defined at element nodes.
   
  @param eid[in] - element id
  @param sm[out] - stiffness %matrix
  
  @return The function returns resulting stiffness %matrix in the argument sm. 

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_stiffness_matrix(long eid, matrix &sm)
{
  ivector nodes(ASTCKIVEC(nne));

  // assemble stiffness matrix in the global coordinate system
  stiffness_matrix (eid,0,0,sm);

  //  transformation of stiffness matrix to the local nodal coord. systems
  Mt->give_elemnodes (eid,nodes);
  long transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe, ndofe));
    transf_matrix(nodes, tmat);
    glmatrixtransf(sm, tmat);
  }
}



/**
  The function computes mass %matrix of the rectangular axisymmetric
  finite element with bicubic approximation functions.
   
  @param eid[in] - number of element
  @param mm[out] - mass %matrix

  @return The function returns resulting mass %matrix in the argument mm.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::mass_matrix(long eid,matrix &mm)
{
  long i, j;
  double jac, xi, eta, rho, r;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordmm)), gp(ASTCKVEC(intordmm));
  vector t(ASTCKVEC(nne)), dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes(eid, nodes);
  Mc->give_density(eid, nodes, dens);
  Mt->give_node_coord2d(x, y, eid);
  gauss_points(gp.a, w.a, intordmm);
  
  nullm(mm);

  for (i=0; i<intordmm; i++){
    xi=gp[i];
    for (j=0; j<intordmm; j++){
      eta=gp[j];
      jac_2d(jac, x, y, xi, eta);
      bf_matrix(n, xi, eta);
      
      rho = approx(xi, eta, dens);
      r = approx(xi, eta, x);
      jac *= w[i]*w[j]*rho*r;
      
      nnj(mm.a, n.a, jac, n.m, n.n);
    }
  }
}



/**
*/
void axisymcq::res_mainip_strains(long lcid,long eid)
{
  mainip_strains (lcid,eid,0,0);
}



/**
  The function computes strains in main integration points of element.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id
  @param ri[in]   - row index
  @param ci[in]   - column index
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::mainip_strains(long lcid, long eid, long ri, long ci)
{
  long i, j, ipp;
  double xi, eta, jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector gp, w, eps(ASTCKVEC(tncomp)), aux;
  ivector nodes(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;

  eldispl(lcid, eid, r.a);
  
  //  transformation of displacement vector
  Mt->give_elemnodes(eid, nodes);
  long transf = Mt->locsystems(nodes);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  
  // retrieve Gauss point coordinates and weights
  reallocv(RSTCKVEC(intordsm[0][0], gp));
  reallocv(RSTCKVEC(intordsm[0][0], w));    
  gauss_points(gp.a, w.a, intordsm[0][0]);
  // assemble nodal coordinate vectors
  Mt->give_node_coord2d(x, y, eid);

  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordsm[0][0]; i++){
    xi=gp[i];
    for (j=0; j<intordsm[0][0]; j++){
      eta=gp[j];
      // compute geometric matrix B
      geom_matrix(gm, x, y, xi, eta, jac);
      // compute strains eps = B.d
      mxv (gm, r, eps);
      // store strains at the given integration point of element
      Mm->storestrain(lcid, ipp, eps);
      ipp++;
    }
  }
}



/**
  The function computes strains in nodes of element. Strain values are taken from 
  int. points which are the closest to the particular nodes. Nodal strain values are 
  averaged according to setup.

  @param lcid[in] - load case id
  @param eid[in]  - element id
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nod_strains_ip(long lcid, long eid)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double area;
  vector eps(ASTCKVEC(tncomp)), aux(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes(eid, nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain(lcid, ipnum[i], eps);
    // number of i-th node of element
    j=nod[i];
    // nodal strain values averageing
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2){
      area = Mt->give_area(eid)*1.0/12.0;
      cmulv (area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes strains at nodes of element, nodal strain values are 
  averaged according to setup.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nod_strains_comp(long lcid, long eid)
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
  long transf = Mt->locsystems(enod);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(enod, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux,r);
  }
  
  //  natural coordinates of nodes of element
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planecq (nxi, neta);

  for (i=0;i<nne;i++){
    //  block of geometric matrix
    geom_matrix(gm,x,y,nxi[i],neta[i],jac);
    //  strain computation
    mxv(gm,r,eps);
    //  storage of nodal strains
    j=enod[i];
    // nodal strain values averageing
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2){
      area = Mt->give_area(eid)*1.0/12.0;
      cmulv(area, eps);
      Mt->nodes[j].storestrain(lcid, 0, area, eps);
    }
  }
}



/**
  The function computes strains at nodes of element, no averaging is performed
   
  @param lcid - load case id
  @param eid  - element id
   
  Created by Tomas Koudelka, 13.03.2019
*/
void axisymcq::nod_strains(long lcid, long eid)
{
  long i, j, k, m;
  double jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  ivector enod(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)), tmat;

  //  natural coordinates of nodes of element
  nodcoord_planecq(nxi, neta);
  //  node numbers of element
  Mt->give_elemnodes(eid, enod);
  //  coordinates of element nodes
  Mt->give_node_coord2d(x, y, eid);
  //  nodal displacements
  eldispl(lcid, eid, r.a);
  
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
    geom_matrix(gm, x, y, nxi[i], neta[i], jac);
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
  The function computes strains in all integration points.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id
   
  @return The function does not return anything but actualizes strain arrays of 
          all integration points of the given element.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_allip_strains(long lcid, long eid)
{
  // all strain components at all integration points of all blocks on element
  allip_strains(lcid, eid, 0, 0);
}



/**
  The function computes strains in all integration points of the given block.
   
  @param lcid[in]       - load case id
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices
   
  @return The function does not return anything but actualizes strain arrays of 
          all integration points of the given element block.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::allip_strains(long lcid, long eid, long /*ri*/, long /*ci*/)
{
  //  blocks of strain components at integration points
  res_mainip_strains (lcid, eid);
}



/**
  The function computes starins at element points according to the setup in mechmat.
  
  @param lcid[in]       - load case id
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices

  @return The function does not return anything but it stores computed strains at the 
          given element block points.
*/
void axisymcq::strains(long lcid, long eid, long ri, long ci)
{
  vector coord, eps;
  
  switch (Mm->stra.tape[eid]){
    case nowhere:
      break;
    case intpts:
      allip_strains(lcid, eid, ri, ci);
      break;
    case enodes:
      nod_strains_ip (lcid,eid);
      break;
    case userdefined:
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
    default:
      print_err("unknown strain point type %d is required", __FILE__, __LINE__, __func__, Mm->stra.tape[eid]);
  }
}



/**
  The function returns the integration point numbers closest to element nodes.
   
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices
  @param ipnum[out]     - array of numbers
   
  @return The function returns the integration point numbers in the argument ipnum.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nodipnum(long eid, ivector &ipnum)
{
  long i, j;
  
  j=intordsm[0][0];
  i=Mt->elements[eid].ipp[0][0];
  
  ipnum[0]=i+j*(j-1)+j-1;
  ipnum[1]=i+j-1;
  ipnum[2]=i;
  ipnum[3]=i+j*(j-1);
  switch(j){
    case 3:
      ipnum[4]  = i+5;
      ipnum[5]  = i+5;
      ipnum[6]  = i+1;
      ipnum[7]  = i+1;
      ipnum[8]  = i+3;
      ipnum[9]  = i+3;
      ipnum[10] = i+7;
      ipnum[11] = i+7;
      break;    
    case 4:
      ipnum[4]  = i+7;
      ipnum[5]  = i+11;
      ipnum[6]  = i+2;
      ipnum[7]  = i+1;
      ipnum[8]  = i+4;
      ipnum[9]  = i+8;
      ipnum[10] = i+13;
      ipnum[11] = i+14;
      break;
    default:
      print_err("order of integration %ld has not yet been implemented", __FILE__, __LINE__, __func__, j);
      break;
  }
}



/**
  The function computes strains in arbitrary point on element by approximation 
  of nodal strain values.
   
  @param xi[in], eta[in] - natural coordinates of the point
  @param fi[in]          - index of the first required strain component
  @param nc[in]          - number of required strain components
  @param eps[out]        - array containing strains
  @param val[in]         - array containing strain values at element nodes
                           val[i][j] represents j-th strain component at 
                           i-th element node
   
  @return The function returns computed values in the argument eps.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::appval(double xi, double eta, long fi, long nc, vector &eps, double **val)
{
  long i, j, k;
  vector nodval;
  
  k=0;
  reallocv(RSTCKVEC(nne, nodval));
  for (i=fi; i<fi+nc; i++){
    for (j=0; j<nne; j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx(xi, eta, nodval);
    k++;
  }
}



/**
  The function computes stresses at main integration points at all element blocks.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id

  @return The function does not return anything but it stores computed stresses at the 
          main integration points of the given element.
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_mainip_stresses(long lcid, long eid)
{
  mainip_stresses(lcid, eid, 0, 0);
}



/**
  The function computes stresses in main integration points of element of the given block.
   
  @param lcid[in] - load case id
  @param eid[in]  - element id
  @param ri[in]   - row index of integration point block
  @param ci[in]   - column index of integration point block
   
  @return The function does not return anything but it stores computed stresses at the 
          given element block of main integration points.
*/
void axisymcq::mainip_stresses(long lcid, long eid, long ri, long ci)
{
  long i, j, ipp;
  vector eps(ASTCKVEC(tncomp)), sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp, tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for(i=0; i<intordsm[0][0]; i++){
    for (j=0; j<intordsm[0][0]; j++){
      // compute material stiffness matrix at the given ip
      Mm->matstiff(d, ipp);
      // retrieve strains at the given ip
      Mm->givestrain(lcid, ipp, eps);
      // compute stresses: sig = D.eps
      mxv(d, eps, sig);
      // store stresses to the given ip
      Mm->storestress(lcid, ipp, sig);
      // increase ip number
      ipp++;
    }
  }
}



/**
  The function computes stresses at nodes of element , the values a copied 
  from the closest integration point to node. 

  @param lcid[in] - load case id
  @param eid[in] - element id
   
  @return The function does not return anything but it stores computed stresses at the 
          element nodes.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nod_stresses_ip (long lcid, long eid)
{
  long i, j;
  double area;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(tncomp)), aux(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes(eid, nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress(lcid, ipnum[i], sig);    
    //  storage of stresses to the node
    j=nod[i];
    // nodal stress values averageing
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2){
      area = Mt->give_area(eid)*1.0/12.0;
      cmulv(area, sig);
      Mt->nodes[j].storestress(lcid,0,area,sig);
    }
  }
}



/**
  The function computes stresses at nodes of the eid-th element. Nodal stress values are averaged 
  according to the setup.
   
  @param lcid[in] - load case id
  @param eid - element id
   
  @return The function does not return anything but it stores computed stresses at the 
          element nodes.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nod_stresses_comp (long lcid,long eid)
{
  long i, j;
  double jac, area;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp)), aux;
  ivector ipnum(ASTCKIVEC(nne)), nodes(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp, ndofe)), tmat;
  matrix d(ASTCKMAT(tncomp, tncomp)), auxd(ASTCKMAT(tncomp, tncomp));
  
  //  natural coordinates of nodes of element
  nodcoord_planecq(nxi, neta);
  //  node numbers of element
  Mt->give_elemnodes(eid, nodes);
  //  coordinates of element nodes
  Mt->give_node_coord2d(x, y, eid);
  //  nodal displacements
  eldispl(lcid, eid, r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems(nodes);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv (aux, r);
  }
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ipnum);

  for (i=0; i<nne; i++){
    //  block of geometric matrix
    geom_matrix(gm, x, y, nxi[i], neta[i], jac);
    //  strain computation
    mxv(gm, r, eps);
    //  stiffness matrix of the material
    Mm->matstiff(d, ipnum[i]);
    //  stress computation
    mxv(d, eps, sig);    
    //  number of actual node
    j=nodes[i];
    //  storage of nodal strains
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress(lcid, 0, sig);
    if (Mp->stressaver==2){
      area = Mt->give_area(eid)*1.0/12.0;
      cmulv(area, sig);
      Mt->nodes[j].storestress(lcid, 0, area, sig);
    }
  }
}



/**
  The function computes stresses in all integration points of all element ip blocks.
   
  @param lcid[in]       - load case id
  @param eid[in]        - element id
   
  @return The function does not return anything but it stores computed stresses at the 
          all element integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_allip_stresses (long lcid, long eid)
{
  allip_stresses(lcid, eid, 0, 0);
}



/**
  The function computes stresses in all integration points of the given ip block of element.
   
  @param lcid[in] - load case id
  @param eid[in] - element id
  @param ri[in], ci[in] - row and column indices of the given ip block
   
  @return The function does not return anything but it stores computed stresses at the 
          the given element integration point block.
*/
void axisymcq::allip_stresses(long lcid, long eid, long /*ri*/, long /*ci*/)
{
  res_mainip_stresses(lcid, eid);
}




void axisymcq::stresses(long lcid, long eid, long ri, long ci)
{
  vector coord,sig;

  switch (Mm->stre.tape[eid]){
    case nowhere:
      break;
    case intpts:
      allip_stresses(lcid, eid, ri, ci);
      break;
    case enodes:
      nod_stresses_ip(lcid, eid);
      break;
    case userdefined:
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
    default:
      print_err("unknown stress point %d is required", __FILE__, __LINE__, __func__, Mm->stre.tape[eid]);
      break;
  }
}



/**
  The function computes other components at nodes of the eid-th element. Nodal values are averaged 
  according to the setup.

  @param eid[in] - element id
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nod_other_ip(long eid)
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
    // the node number of i-th element node    
    j=nod[i];
    // nodal values averageing
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother (0, ncompo, other);
    if (Mp->otheraver==2){
      area = Mt->give_area(eid)*1.0/12.0;
      cmulv (area, other);
      Mt->nodes[j].storeother(0, ncompo, area, other);
    }    
  }
}



/**
  The function computes load %matrix of the axisymmetric quadrilateral
  finite element with bicubic approximation functions. Load vector is obtained 
  after premultiplying load %matrix by nodal load values.
   
  @param eid[in] - number of element
  @param lm[out] - load matrix

  @return The function returns computed load matrix in the argument lm.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::load_matrix(long eid, matrix &lm)
{
  long i, j;
  double jac, xi, eta, r;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordmm)), gp(ASTCKVEC(intordmm));
  matrix n(ASTCKMAT(napfun, ndofe));
  
  Mt->give_elemnodes(eid, nodes);
  Mt->give_node_coord2d(x, y, eid);
  gauss_points(gp.a, w.a, intordmm);
  
  nullm(lm);

  for (i=0; i<intordmm; i++){
    xi=gp[i];
    for (j=0; j<intordmm; j++){
      eta=gp[j];
      jac_2d(jac, x, y, xi, eta);
      bf_matrix (n, xi, eta);
      r = approx(xi, eta, x);
      jac *= r*w[i]*w[j];
      nnj(lm.a, n.a, jac, n.m, n.n);
    }
  }
}



/**
  The function computes %vector internal forces on the given element.

  @param lcid[in]       - load case id
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
  @param ifor[out]      - resulting %vector of internal forces
  @param x[in], y[in]   - %vectors of nodal coordinates

  @return The function returns nodal values of internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::internal_forces(long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y)
{
  integratedquant iq=locstress;

  //  computation of stresses
  compute_nlstress(lcid, eid, ri, ci);
  //  integration of stresses over the element
  elem_integration(iq, lcid, eid, ri, ci, ifor, x, y);
}



/**
  The function computes %vector internal forces for nonlocal models on the given element.

  @param lcid[in]       - load case id
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
  @param ifor[out]      - resulting %vector of internal forces
  @param x[in], y[in]   - vectors of nodal coordinates

  @return The function returns computed internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::nonloc_internal_forces(long lcid, long eid, long ri, long ci, vector &ifor,
                                      vector &x, vector &y)
{
  integratedquant iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress(lcid, eid, ri, ci);
  //  integration of stresses over the element
  elem_integration(iq, lcid, eid, ri, ci, ifor, x, y);
}



/**
  The function computes increments of internal forces going to the right hand side due 
  to known values of strains or stresses in the material model 
  (e.g. creep strains, pore pressures, ...).

  @param lcid[in]       - load case id
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
  @param ifor[out]      - vector of internal forces
  @param x[in], y[in]   - vectors of nodal coordinates

  @return The function returns computed increments of internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::incr_internal_forces(long lcid, long eid, long ri, long ci, vector &ifor, 
                                    vector &x, vector &y)
{
  integratedquant iq=stressincr;

  //  computation of stresses
  compute_nlstressincr(lcid, eid, ri, ci);
  //  integration of stresses over the element
  elem_integration (iq, lcid, eid, ri, ci, ifor, x, y);
}



/**
  The function computes nodal forces caused by temperature changes or prescribed eigenstrains.
   
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
  @param nfor[out]      - array containing nodal forces
  @param x[in], y[in]   - nodal coordinates
   
  @return The function returns computed nodal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::eigstrain_forces(long lcid, long eid, long ri, long ci, vector &nfor,
                                vector &x, vector &y)
{
  integratedquant iq=eigstress;
  
  //  computation of eigenstresses
  if (Mp->eigstrcomp)
    compute_eigstress(lcid, eid, ri, ci);
  //  integration of stresses over the element
  elem_integration(iq, lcid, eid, ri, ci, nfor, x, y);
}



/**
  The function computes resulting internal forces with respect to local coordinate systems 
  defined at element nodes.
   
  @param lcid[in]  - load case id
  @param eid[in]   - element id
  @param ifor[out] - %vector of internal forces
   
  @return The function returns computed internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_internal_forces(long lcid, long eid, vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  // retrieve nodal coordinates
  Mt->give_node_coord2d(x, y, eid);
  // compute internal forces at nodes given in the global coordinate system
  internal_forces(lcid, eid, 0, 0, ifor, x, y);

  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe, ndofe));
    transf_matrix(nodes, tmat);
    glvectortransf(ifor, v, tmat);
    copyv(v, ifor);
  }
}



/**
  The function computes resulting internal forces for nonlocal models with respect to 
  local coordinate systems defined at element nodes.
   
  @param lcid[in]  - load case id
  @param eid[in]   - element id
  @param ifor[out] - resulting %vector of internal forces
   
  @return The function returns computed internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_nonloc_internal_forces (long lcid, long eid, vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  // retrieve nodal coordinates  
  Mt->give_node_coord2d(x, y, eid);
  // compute internal forces for nonlocal material models in the global coordinate system
  nonloc_internal_forces(lcid, eid, 0, 0, ifor, x, y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe, ndofe));
    transf_matrix(nodes, tmat);
    glvectortransf(ifor, v, tmat);
    copyv(v, ifor);
  }
}



/**
  The function computes resulting increment of internal forces going to the right hand side due 
  to known values of strains or stresses in the material model (e.g. creep strains, pore pressures, ...).
  The %vector is computed with respect to local coordinates systems defined at element nodes.
   
  @param lcid[in]  - load case id
  @param eid[in]   - element id
  @param ifor[out] - resulting %vector of internal forces increments
   
  @return The function returns computed internal forces in the argument ifor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_incr_internal_forces (long lcid, long eid, vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

  // retrieve nodal coordinates  
  Mt->give_node_coord2d(x, y, eid);
  // compute increments of internal forces
  incr_internal_forces(lcid, eid, 0, 0, ifor, x, y);

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
  The function computes resulting contributions from eigenstrains to the right hand side with
  respect to local coordinate systems defined at element nodes.
   
  @param lcid[in] - load case id
  @param eid[in] - element id
  @param nfor[out] - resulting %vector of internal forces

  @return The function returns computed nodal forces in the argument nfor.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::res_eigstrain_forces (long lcid, long eid, vector &nfor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

  // retrieve nodal coordinates  
  Mt->give_node_coord2d(x, y, eid);
  // compute nodal forces due to eigenstrains in the global coordinate system
  eigstrain_forces(lcid, eid, 0, 0, nfor, x, y);

  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe, ndofe));
    transf_matrix(nodes, tmat);
    glvectortransf(nfor, v, tmat);
    copyv(v, nfor);
  }
}



/**
  The function computes stresses at the given block of integration points on element.

  @param lcid[in]       - number of load case
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
   
  @return The function does not return anything but it stores computed values in 
          the integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::compute_nlstress(long /*lcid*/, long eid, long ri, long ci)
{
  long i, j, ipp;
  
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  for (i=0; i<intordsm[0][0]; i++){
    for (j=0; j<intordsm[0][0]; j++){
      if (Mp->strcomp==1)
	Mm->computenlstresses(ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}



/**
  The function computes increments of stresses (e.g. due to creep, pore pressure) at 
  the given integration point block on element.  

  @param lcid[in]       - number of load case
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
   
  @return The function does not return anything but it stores computed values in 
          the integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::compute_nlstressincr(long /*lcid*/, long eid, long ri, long ci)
{
  long i, j, ipp;
  
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  for (i=0; i<intordsm[0][0]; i++){
    for (j=0; j<intordsm[0][0]; j++){
      if (Mp->strcomp==1)
	Mm->computenlstressesincr(ipp);
      ipp++;
    }
  }
}



/**
  The function computes local values which will be used for averaging 
  in nonlocal models at integration points. Mp->nonlocphase have to be 1.

  @param lcid[in]      - number of load case
  @param eid[in]       - element id
  @param ri[in],ci[in] - row and column indices of the required integration point block
   
  @return The function does not return anything but it stores computed values in 
          the integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::local_values(long /*lcid*/, long eid, long ri, long ci)
{
  long i,j,ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      if (Mp->strcomp==1)
	Mm->computenlstresses(ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}



/**
  The function computes nonlocal stresses at the given integration point block on element.
   
  @param lcid[in]       - number of load case
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
   
  @return The function does not return anything but it stores computed values in 
          the integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::compute_nonloc_nlstress(long /*lcid*/, long eid, long ri, long ci)
{
  long i, j, ipp;
  
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  for (i=0; i<intordsm[0][0]; i++){
    for (j=0; j<intordsm[0][0]; j++){
      if (Mp->strcomp==1)
	Mm->computenlstresses(ipp,Mm->ip[ipp]);
      ipp++;
    }
  }
}



/**
  The function computes eigen stresses caused by temperature at the given integration 
  point block on element
   
  @param lcid[in] - number of load case
  @param eid[in] - element id
  @param ri[in], ci[in] - row and column indices of the required integration point block
   
  @return The function does not return anything but it stores computed values in 
          the integration points.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
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
  The function integrates selected stress quantity over the finite element volume
  which results in nodal force values. The \int_V B^T sig dV is being computed 
  by this function.
   
  @param iq[in]         - type of integrated quantity (see alias.h)
  @param lcid[in]       - number of load case
  @param eid[in]        - element id
  @param ri[in], ci[in] - row and column indices
  @param nv[in]         - nodal values
  @param x[in], y[in]   - nodal coordinates
   
  @return The function returns nodal values of integrated quantity 

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::elem_integration(integratedquant iq, long lcid, long eid, long /*ri*/, long /*ci*/, 
                                vector &nv, vector &x, vector &y)
{
  long i, j, ipp;
  double xi, eta, jac, rad;
  vector w, gp;
  vector ipv(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp, ndofe));
  
  nullv(nv);
  reallocv(RSTCKVEC(intordsm[0][0], gp));
  reallocv(RSTCKVEC(intordsm[0][0], w));
  gauss_points(gp.a, w.a, intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (i=0; i<intordsm[0][0]; i++){
    xi=gp[i];
    for (j=0; j<intordsm[0][0]; j++){
      eta=gp[j];
      //  function assembles required quantity at integration point
      Mm->givequantity(iq, lcid, ipp, 0, ipv);
      //  strain-displacement (geometric) matrix
      geom_matrix(gm, x, y, xi, eta, jac);
      //  contribution to the nodal values
      mtxv(gm, ipv, contr);
      rad = approx(xi, eta, x);
      cmulv(rad*jac*w[i]*w[j], contr);
      //  summation
      addv(contr, nv, nv);
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

   Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv)
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
  The function integrates N^T coef N over edges
   
  @param edg[in] - edge id (number of edge)
  @param x[in], y[in] - coordinates of element nodes
  @param intord[in] - order of numerical integration
  @param gp[in], w[in]  - coordinates and weights of integration points
  @param coef[in] - array of nodal values of coefficient
  @param km[out] - output %matrix
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::edge_integral (long edg, vector &x, vector &y, long /*intord*/, vector &gp, vector &w,
			      vector &coef, matrix &km)
{
  long i;
  double xi, eta, jac, ipval, r;
  matrix n(ASTCKMAT(napfun,ndofe));
  
  if (edg==0){
    eta = 1.0;
    for (i=0; i<intordb; i++){
      xi = gp[i];
      // matrix of approximation functions
      bf_matrix(n, xi, eta);
      // computation of jacobian      
      jac1d_2d(jac, x, y, xi, edg);
      // interpolation of c to integration point
      ipval=approx(xi, eta, coef);
      // radius
      r=approx(xi, eta, x);
    
      jac *= w[i]*ipval*r;
      nnj(km.a, n.a, jac, n.m, n.n);
    }
  }

  if (edg==1){
    xi = -1.0;
    for (i=0; i<intordb; i++){
      eta = gp[i];
      // matrix of approximation functions
      bf_matrix(n, xi, eta);
      // computation of jacobian      
      jac1d_2d(jac, x, y, xi, edg);
      // interpolation of coef to integration point
      ipval = approx(xi, eta, coef);
      // radius
      r=approx(xi, eta, x);

      jac *= w[i]*ipval*r;
      nnj(km.a, n.a, jac, n.m, n.n);
    }
  }

  if (edg==2){
    eta = -1.0;
    for (i=0; i<intordb; i++){
      xi = gp[i];
      // matrix of approximation functions
      bf_matrix(n, xi, eta);
      // computation of jacobian      
      jac1d_2d(jac, x, y, xi, edg);
      // interpolation of c to integration point
      ipval=approx(xi, eta, coef);
      //  radius
      r=approx(xi, eta, x);
      
      jac *= w[i]*ipval*r;
      nnj(km.a, n.a, jac, n.m, n.n);
    }
  }

  if (edg==3){
    xi = 1.0;
    for (i=0; i<intordb; i++){
      eta = gp[i];
      
      //  matrix of approximation functions
      bf_matrix(n, xi, eta);
      //  computation of jacobian      
      jac1d_2d(jac, x, y, xi, edg);
      //  interpolation of c to integration point
      ipval=approx(xi, eta, coef);
      //  radius
      r=approx(xi, eta, x);
      
      jac *= w[i]*ipval*r;
      nnj(km.a, n.a, jac, n.m, n.n);
    }
  }
}



/**
  The function picks up nodal values on required edge.
   
  @param edg[in] - the number of required edge
  @param nodval[out] - array of nodal values
  @param list[in] - array of nodal values defined on all edges
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::edgenodeval(long edg, vector &nodval, double *list)
{
  long i, j, k, l;
  ivector edgenod(ASTCKIVEC(nned));
  
  nullv(nodval);
  //  nodes on required edge
  cubicquadrilat_edgnod(edgenod.a, edg);
  
  k=0;
  for(i=0; i<nned; i++){
    l = edgenod[i]*napfun;
    for(j=0; j<napfun; j++){
      nodval[l+j] = list[edg*nned*napfun+k];
      k++;
    }
  }
}



/**
  The function computes nodal forces caused by edge load.
   
  @param eid[in] - element id
  @param le[in] - array of indicators of loaded edges
  @param list[in] - list of prescribed nodal values
  @param nf[out] - resulting %vector of nodal forces caused by edge load
   
  @return The function returns %vector of computed nodal forces in the argument nf.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::edgeload(long eid, long *le, double *list, vector &nf)
{
  long i;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), gp(ASTCKVEC(intordb)), w(ASTCKVEC(intordb));
  vector nodval(ASTCKVEC(ndofe)), coeff(ASTCKVEC(nne)), av(ASTCKVEC(ndofe));
  matrix km(ASTCKMAT(ndofe,ndofe));
  
  //  coefficient c in N^T c N
  for (i=0; i<nne; i++){
    coeff[i]=1.0;
  }
  //  coordinates of element nodes
  Mt->give_node_coord2d(x, y, eid);
  //  coordinates and weights of integration points
  gauss_points(gp.a, w.a, intordb);
  
  //  cleaning of nf
  nullv(nf);
  
  for (i=0; i<ned; i++){
    if (le[i]>0){
      nullm(km);
      //  selection of appropriate nodal values
      edgenodeval(i, nodval, list);
      //  integration over element edge
      edge_integral(i, x, y, intordb, gp, w, coeff, km);
      // compute contribution from i-th edge
      mxv(km, nodval, av);
      addv(nf, av, nf);
    }
  }
}



/**
  The function assembles global coordinates of integration point.
   
  @param eid[in]    - element id
  @param ipp[in]    - integration point pointer
  @param ri[in]     - row index of integration point block
  @param ci[in]     - column index of integration point block
  @param coord[out] - %vector with global coordinates of integration point
   
  @return The function returns global coordinates in the argument coord.

  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, j, ii;
  double xi, eta;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordsm[ri][ci])), gp(ASTCKVEC(intordsm[ri][ci]));

  gauss_points(gp.a, w.a, intordsm[ri][ci]);
  Mt->give_node_coord2d(x, y, eid);
  ii = Mt->elements[eid].ipp[ri][ci];

  for(i=0; i<intordsm[ri][ci]; i++){
    xi=gp[i];
    for(j=0; j<intordsm[ri][ci]; j++){
      eta=gp[j];
      if (ii==ipp){
	coord[0] = approx(xi, eta, x);
	coord[1] = approx(xi, eta, y);
	coord[2] = 0.0;
      }
      ii++;
    }
  }
}



/**
   The function assembles natural coordinates of integration point.
   
   @param eid[in]     - element id
   @param ipp[in]     - integration point pointer
   @param ri[in]      - row index of integration point block
   @param ci[in]      - column index of integration point block
   @param ncoord[out] - %vector with natural coordinates of integration point
   
   @return The function returns natural coordinates in the argument ncoord.

   Created by TKo, 12.2016
*/
void axisymcq::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, j, ii, ri, ci;
  double xi, eta;
  vector w, gp;

  for (ri=0; ri<nb; ri++){
    for (ci=0; ci<nb; ci++){
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
void axisymcq::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
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
  for (i=0; i<nne; i++){
    aux = int(ictn[i])-int(ict);
    if (aux < 0)  aux = -aux;    
    aux &= ~(inidisp);
    aux &= ~(inidisp_x);
    aux &= ~(inidisp_y);
    aux &= ~(inidisp_z);
    if ((ictn[i] != ict) && aux){
      print_err("Incompatible types of initial conditions on element %ld\n"
                " at %ld. and %ld. nodes", __FILE__, __LINE__, __func__, eid+1, 1, i+1);
      abort();
    }
  }
  for (j = 0; j < nv; j++){ // for all initial values
    for(i = 0; i < nne; i++)
      anv[i] = nodval[i][j];
    for (ii = 0; ii < nb; ii++){
      for (jj = 0; jj < nb; jj++){
        ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
        if (intordsm[ii][jj] == 0)
          continue;
        reallocv (RSTCKVEC(intordsm[ii][jj],gp));
        reallocv (RSTCKVEC(intordsm[ii][jj],w));
        gauss_points (gp.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++){
          xi=gp[k];
          for (l = 0; l < intordsm[ii][jj]; l++){
            eta=gp[l];
            //  value in integration point
            ipval = approx (xi,eta,anv);
            ncompstr =  Mm->ip[ipp].ncompstr;
            ncompeqother = Mm->ip[ipp].ncompeqother;
            if ((ictn[0] & inistrain) && (j < ncompstr)){
              Mm->ip[ipp].strain[idstra] += ipval;
              ipp++;
              continue;
            }
            if ((ictn[0] & inistress) && (j < nstra + ncompstr)){
              Mm->ip[ipp].stress[idstre] += ipval;
              ipp++;
              continue;
            }
            if ((ictn[0] & iniother) && (j < nstra+nstre+ncompeqother)){
              Mm->ip[ipp].eqother[idoth] += ipval;
              ipp++;
              continue;
            }
            if ((ictn[0] & inicond) && (j < nv)){
              if (Mm->ic[ipp] == NULL){
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
    if ((ictn[0] & inistrain) && (j < ncompstr)){
      nstra++;
      idstra++;
      continue;
    }
    if ((ictn[0] & inistress) && (j < nstra + ncompstr)){
      nstre++;
      idstre++;
      continue;
    }  
    if ((ictn[0] & iniother)  && (j < nstra + nstre + ncompeqother)){
      idoth++;
      continue;
    }
    if ((ictn[0] & inicond) && (j < nv)){
      idic++;
      continue;
    }
  }
}



/**
  The function interpolates the nodal values to the integration points on the element
   
  @param eid[in]    - element id
  @param nodval[in] - %vector of nodal values of quantity interpolated
  @param ipval[out] - array of quantity values in all element integration points

  @return The function returns computed values in the argument ipval.
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::intpointval(long /*eid*/, vector &nodval, vector &ipval)
{
  long i, j, ii, jj, k;
  double xi, eta;
  vector w, gp;
  
  k=0;
  for (ii=0; ii<nb; ii++){
    for (jj=0; jj<nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      gauss_points(gp.a, w.a, intordsm[ii][jj]);

      for (i=0; i<intordsm[ii][jj]; i++){
	xi = gp[i];
	for (j=0; j<intordsm[ii][jj]; j++){
	  eta = gp[j];
	  ipval[k] = approx(xi, eta, nodval);
	  k++;
	}
      }
    }
  }
}



/**
  The function interpolates the nodal values given on linear element 2D
  to the integration points on the element. The function uses linear interpolation function 
  defined on the corresponding linear element type.
   
  @param eid[in]    - element id
  @param nodval[in] - %vector of nodal values of quantity interpolated defined on the linear element type
  @param ipval[out] - array of quantity values in all element integration points

  @return The function returns computed values in the argument ipval.
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::intpointval2(long /*eid*/, vector &nodval, vector &ipval)
{
  long i,j,ii,jj,k;
  double xi,eta;
  vector w,gp;
  vector modnodval(ASTCKVEC(Asymlq->nne));
  
  for (i=0; i<Asymlq->nne; i++){
    modnodval[i]=nodval[i];
  }
  
  k=0;
  for (ii=0; ii<nb; ii++){
    for (jj=0; jj<nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      gauss_points(gp.a, w.a, intordsm[ii][jj]);

      for(i=0; i<intordsm[ii][jj]; i++){
	xi = gp[i];
	for(j=0; j<intordsm[ii][jj]; j++){
	  eta = gp[j];
	  ipval[k] = Asymlq->approx(xi, eta, modnodval);
	  k++;
	}
      }
    }
  }
}



/**
  The function computes required mechanical quantity at nodes of element.

  @param eid[in] - element id
  @param nodval[out] - %vector of nodal values
  @param qt[in] - type of mechanical quantity
   
  @return The function returns computed quantity values in the nodval argument.
  
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::mechq_nodval(long eid, vector &nodval, nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ipnum);

  for (i=0; i<nne; i++){
    //copy nonmechanical quantity from closest int. point
    nodval[i] = Mm->givemechq(qt, ipnum[i]);
  }
}



/**
  The function computes required mechanical quantity at nodes of element.
  linear approximation functions are used.

  @param eid[in] - element id
  @param nodval[out] - %vector of nodal values
  @param qt[in] - type of mechanical quantity
   
  @return The function returns computed quantity values in the nodval argument.
   
  Created by Tomas Koudelka, 13.3.2019
*/
void axisymcq::mechq_nodval2(long eid, vector &nodval, nontransquant qt)
{
  long i;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ipnum);

  for (i=0; i<Asymlq->nne; i++){
    //copy nonmechanical quantity from closest int. point
    nodval[i] = Mm->givemechq(qt, ipnum[i]);
  }
}



/**
  The function computes mechanical quantities in nodes of element.

  @param eid[in]     - element id
  @param nodval[out] - %vector of nodal values of all required quantities, i.e., 
                       nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                       is the number of calculated nodes on eid-th element.
  @param ncnv[in] - number of computed nodes on element (only first ncnv of nodes is calculated)
  @param nq[in]   - number of required mechanical quantities
  @param qt[in]   - array of types of required mechanical quantities
   
  @return The function does not return anything.
   
   01/11/2016 TKr accordintg to Tomas Koudelka and axisymlq
*/
void axisymcq::mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt)
{
  long i, j, ncompstr, ncompo, ncompnl, ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  intpoints ipb; // backup of the first integration point of element
  
  // id of the first integration point on element
  ipid=Mt->elements[eid].ipp[0][0];

  // element nodes
  Mt->give_elemnodes(eid, enod);

  //  numbers of integration points closest to element nodes
  nodipnum(eid,ipnum);

  // compute strains at nodes
  nod_strains(0, eid);
  
  // number of nonloc array components
  ncompnl = Mm->give_num_averq(ipid, Mm->givenonlocid(ipid));

  // store original content of the first integration point on element because 
  // it becomes working int. point for nodal values calculations on the given element
  ipb.copy(Mm->ip[ipid], Mb->nlc, ncompnl, 1);


  // The first integration point will be used for computation of nodal values temporarily
  // then the original content of the first integration point will be restored
  for (i=0;i<ncnv;i++){
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
    switch(Mp->matmodel){
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
