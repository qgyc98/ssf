#include "axisymlqiface.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "element.h"
#include "node.h"
#include "globmat.h"
#include "intp.h"
#include "basefun.h"

axisymlqinterface::axisymlqinterface(void)
{
  long i;

  nne=4;  ndofe=8;
  nb=1;  tncomp=3;
  ned=2; nned=2;
  napfun=4;

  ssst = planecontact;

  ncomp = new long [nb];
  ncomp[0]=3;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=2;
  intordsm[0][0]=2;
  tnip=2;
}



axisymlqinterface::~axisymlqinterface(void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;

  delete [] ncomp;
  
}



/**
   The function initializes basic data on elements.
   
   @param[in] eid - element id      
   
   JK, 11.6.2006

void axisymlqinterface::eleminit (long eid)
{

  long ii,jj;

  Mt->elements[eid].nb=nb;
  Mt->elements[eid].intordsm = new long* [nb];
  Mt->elements[eid].nip = new long* [nb];

  for (ii=0;ii<nb;ii++){
    Mt->elements[eid].intordsm[ii] = new long [nb];
    Mt->elements[eid].nip[ii] = new long [nb];
    for (jj=0;jj<nb;jj++){
      Mt->elements[eid].intordsm[ii][jj]=intordsm[ii][jj];
      Mt->elements[eid].nip[ii][jj]=nip[ii][jj];
    }
  }
}
*/


/**
   The function returns value of function in the given point.
   The function is given by nodal values and it is approximated 
   by shape functions of the element. The point is given in the 
   natural coordinate system of the element.

   @param[in] xi - natural coordinate on element
   @param[in] nodval - nodal values
   
   @return The function returns the value approximated.

   JK
*/
double axisymlqinterface::approx(double xi, vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_1d (bf.a,xi);
  // we need approximation as on 1D element with 2 nodes
  bf[2] = 0.0;
  bf[3] = 0.0;
  
  scprd (bf,nodval,f);

  return f;
}



/**
  The function computes integration point values from the given nodal 
  values of selected quantity.
 
  @param[in]  eid    - number of element
  @param[in]  nodval - %vector of nodal values of the given quantity
  @param[out] ipval  - %vector of integration point values approximated form the nodval.
 
  Created by Tomas Koudelka, 13.1.2016
*/
void axisymlqinterface::intpointval(long eid, vector &nodval, vector &ipval)
{
  long ii,jj,i,k,ipp;
  double xi;
  vector w,gp;

  k=0;

  for (ii = 0; ii < nb; ii++)
  {
    for (jj = 0; jj < nb; jj++)
    {
      ipp=Mt->elements[eid].ipp[ii][jj];
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      gauss_points (gp.a, w.a, intordsm[ii][jj]);
      for (i = 0; i < intordsm[ii][jj]; i++)
      {
        xi=gp[i];
        //  value in integration point
        ipval[k] = approx (xi, nodval);
        k++;
        ipp++;
      }
    }
  }
}



/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   
   10.5.2002
*/
void axisymlqinterface::nodecoord (vector &xi)
{
  xi[0] = -1.0;
  xi[1] =  1.0;
  xi[2] =  1.0;
  xi[3] = -1.0;
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
void axisymlqinterface::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, ii;
  double xi;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordsm[ri][ci])), gp(ASTCKVEC(intordsm[ri][ci]));

  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord2d (x,y,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    if (ii==ipp){
      coord[0] = approx(xi, x);
      coord[1] = approx(xi, y);
      coord[2] = 0.0;
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
void axisymlqinterface::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ii, ri, ci;
  double xi;
  vector w, gp;

  for (ri=0; ri<nb; ri++)
  {
    for (ci=0; ci<nb; ci++)
    {
      if (intordsm[ri][ci] == 0)
        continue;

      reallocv(RSTCKVEC(intordsm[ri][ci], w));
      reallocv(RSTCKVEC(intordsm[ri][ci], gp));
      gauss_points(gp.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];

      for (i=0;i<intordsm[ri][ci];i++){
        xi=gp[i];
        if (ii==ipp){
          ncoord[0]=xi;
          ncoord[1]=0.0;
          ncoord[2]=0.0;
        }
        ii++;
      }
    }
  }
}



/**
  The function assembles transformation %matrix x_g = T x_l
  for transformation from the local nodal coordinate systems.
   
  @param[in] nodes - array containing node numbers
  @param[out] tmat - transformation %matrix
   
  JK, 11.6.2006
*/
void axisymlqinterface::transf_matrix(ivector &nod, matrix &tmat)
{
  long i,n,m;

  nullm(tmat);

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
  The function computes stiffness %matrix of plane rectangular
  finite element for contact problems.
   
  @param[in] eid - number of element
  @param[in] ri,ci - row and column indices
  @param[out] sm - resulting stiffness %matrix

  @return The function returns resulting stiffnes %matrix in the parameter sm.
   
  JF, 29.10.2012
*/
void axisymlqinterface::stiffness_matrix(long eid, long ri, long ci, matrix &sm)
{
  long i,ipp;
  double xi,jac,r;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordsm[0][0])), gp(ASTCKVEC(intordsm[0][0]));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp[0], ndofe)), d(ASTCKMAT(tncomp,tncomp));

  nullm(sm);

  Mt->give_node_coord2d(x,y,eid);
  Mt->give_elemnodes (eid,nodes);

  gauss_points (gp.a,w.a,intordsm[0][0]);

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    xi=gp[i];

    geom_matrix (gm, x, y, xi, jac);
    //  matrix of stiffness of the material
    Mm->matstiff (d,ipp);
    //  radius in the given integration point
    r = approx (xi,x);
    jac*=w[i]*r;
    //  contribution to the stiffness matrix of the element
    bdbjac (sm,gm,d,gm,jac);

    ipp++;
  }
}



/**
  Function computes stiffness %matrix of one element. If it is required, nodal values 
  are transformed to the local nodal coordinate systems.

  @param[in] eid - number of element
  @param[out] sm - stiffness %matrix

  @return The function returns required stiffness %matrix in the parameter sm.

  Created by Tomas Koudelka, 11.2008
*/
void axisymlqinterface::res_stiffness_matrix(long eid, matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  
  stiffness_matrix(eid,0,0,sm);
  
  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid,nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix(nodes,tmat);
    glmatrixtransf(sm,tmat);
  }
}



/**
   The function computes strains at integration points. If there are 
   defined loacal nodal coordinate systems then the corresponding 
   components of the dislacement %vector are transformed.
   
   @param[in] lcid - load case id
   @param[in] eid - element id
   
   @return The function stores resulting strains in the array strains 
           of the element integration points.

   JK, 11.6.2006
*/
void axisymlqinterface::res_mainip_strains (long lcid, long eid)
{
 
  vector aux, r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;

  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems(nodes);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  mainip_strains(lcid, eid, 0, 0, r);
}



/**
  The function computes strains at integration points of element.
  It is used in geometrically linear problems.
   
  @param[in] lcid - load case id
  @param[in] eid - element id
  @param[in] ri,ci - row and column indices
  @param[in] r - %vector of nodal displacements
   
  @return The function stores resulting strains in the array strains 
          of the element integration points.

  JK, 11.6.2006
*/
void axisymlqinterface::mainip_strains (long lcid, long eid, long ri, long ci, vector &r)
{
  long i,ipp;
  double xi,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordsm[ri][ci])), gp(ASTCKVEC(intordsm[ri][ci])), t(ASTCKVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp[0],ndofe)), d(ASTCKMAT(tncomp,tncomp));
  vector eps(ASTCKVEC(tncomp));

  nullm (gm);

  Mt->give_node_coord2d(x,y,eid);
  Mt->give_elemnodes (eid,nodes);

  gauss_points (gp.a,w.a,intordsm[ri][ci]);

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++)
  {
    xi=gp[i];

    geom_matrix (gm, x, y, xi, jac);

    mxv (gm, r, eps);

    Mm->storestrain (lcid,ipp,eps);

    ipp++;
  }  
}



/**
  The function returns numbers of integration point closest to element nodes.
   
  @param[in] eid - element id
  @param[in] ri  - row index of the integration point block
  @param[in] ci  - column indices of the integration point block
  @param[out] ipnum - array of numbers
  
  TKo, 2.4.2024
*/
void axisymlqinterface::nodipnum (long eid, long ri, long ci, ivector &ipnum)
{
  long i,j;
  
  j=intordsm[0][0];
  i=Mt->elements[eid].ipp[ri][ci];

  ipnum[0] = ipnum[3] = i+0;
  ipnum[1] = ipnum[2] = i+j-1;
}

/**
   function computes strains in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 2.4.2024
*/
void axisymlqinterface::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double l;
  vector eps(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  nodipnum(eid, ri, ci, ipnum);
  
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
      l = Mt->give_length(eid)*0.5;
      cmulv(l, eps);
      Mt->nodes[j].storestrain(lcid, 0, l, eps);
    }
  }
}



/**
  The function computes nodal strains directly and average them according to setup.
   
  @param lcid - load case id
  @param eid - element id
   
  25/05/2016 TKr according to JK
*/
void axisymlqinterface::nod_strains_comp (long lcid,long eid)
{
  long i, j;
  double jac, l;
  ivector enod(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  matrix tmat,gm(ASTCKMAT(tncomp,ndofe));
  
  //  natural coordinates of nodes of element
  //  (function is from the file GEFEL/ordering.cpp)
  nodecoord (nxi);
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
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix(gm, x, y, nxi[i], jac);
    //  strain computation
    mxv (gm, r, eps);
    
    //  storage of strains to the node
    j=enod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain(lcid, 0, eps);
    if (Mp->strainaver==2)
    {
      l = Mt->give_length(eid)*0.5;
      cmulv(l, eps);
      Mt->nodes[j].storestrain(lcid, 0, l, eps);
    }
  }
}



/**
  The function computes nodal strains directly, no aveargeing is applied.
   
  @param lcid - load case id
  @param eid - element id
   
  25/05/2016 TKr according to JK
*/
void axisymlqinterface::nod_strains (long lcid,long eid)
{
  long i, j, k, m;
  double jac;
  ivector enod(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector nxi(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector eps(ASTCKVEC(tncomp)), aux;
  matrix tmat, gm(ASTCKMAT(tncomp,ndofe));
  
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
  nodecoord (nxi);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,nxi[i],jac);
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
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void axisymlqinterface::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i, j;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double l;
  vector sig(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodipnum (eid,ri,ci,ipnum);
  
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
      l = Mt->give_length(eid)*0.5;
      cmulv(l,sig);
      Mt->nodes[j].storestress(lcid, 0, l, sig);
    }
  }
}



/**
  The function computes other components at nodes of element. Nodal other values components are
  averaged according to the setup.

  @param eid[in] - element id
   
  JK, 10.5.2002
*/
void axisymlqinterface::nod_other_ip (long eid)
{
  long i, j, ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  double l;
  vector other;
  
  //  numbers of integration points closest to nodes
  nodipnum (eid,0,0,ipnum);
  
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
      Mt->nodes[j].storeother (0, ncompo, other);
    if (Mp->otheraver==2)
    {
      l = Mt->give_length(eid)*0.5;
      cmulv(l, other);
      Mt->nodes[j].storeother(0, ncompo, l, other);
    }    
  }
}



/**
  The function computes stresses at element integration points
   
  @param[in] lcid - load case id
  @param[in] eid - element id
   
  @return The function stores resulting stresses in the stress array 
          of the element integration points.
           
  Created by TKo, 11.2012
*/
void axisymlqinterface::res_mainip_stresses (long lcid, long eid)
{
  long ri, ci;
  ri = 0;
  ci = 0;
  compute_nlstress(lcid, eid, ri, ci);
}



/**
   The function assembles strain-displacement (geometric) %matrix.
   
   @param[out] gm - resulting geometric %matrix
   @param[in] x,y - nodal coordinates
   @param[in] xi -  natural coordinate of required integration point
   @param[out] jac - Jacobian of transformation
   
   @return The function returns resulting geometric %matrix in the parameter gm.

   JF, 15.10.2012
   Modified by TKo, 16.12.2015
*/

void axisymlqinterface::geom_matrix (matrix &gm, vector &x, vector &y, double xi, double &jac)
{
  matrix transf(ASTCKMAT(3,2));
  matrix aux(ASTCKMAT(3,2));
  double l,dx,dy;

  dx = x[1] - x[0];
  dy = y[1] - y[0];

  nullm(gm);

  double n1 = 0.5*(1-xi);
  double n2 = 0.5*(1+xi);

  l=sqrt(dx*dx+dy*dy);

  transf(0,0) =  dx/l; 
  transf(0,1) =  dy/l;
  transf(1,0) =  1.0;  // u_phi need not to be transformed
  transf(2,0) = -dy/l; 
  transf(2,1) =  dx/l;

  // node 1 coincides with node 4
  // node 2 coincides with node 3
  // local x axis corresponds to the direction from node 1 to node 2
  // The local y-axis is normal to the local x-axis, i.e. local axis x rotated counterclockwise by pi/2.
  // Nodes 1 and 4 in local element ordering have natural coordinate ksi=-1.
  // Nodes 2 and 3 in local element ordering have natural coordinate ksi=+1.
  // Ordering of the relative displacements: [[u^l]], [[u_phi]], [[v^l]], i.e.:
  // slip in the local x direction,
  // slip in circumferal direction,
  // relative displacement in the normal direction (local y direction).

  cmulm(n1, transf, aux); // B_1 . T = N_1 * I . T = N_1 * T

  gm(0,0) = -aux(0,0);
  gm(0,1) = -aux(0,1);
  gm(1,0) = -aux(1,0);
  gm(1,1) = -aux(1,1);
  gm(2,0) = -aux(2,0);
  gm(2,1) = -aux(2,1);

  gm(0,6) = aux(0,0);
  gm(0,7) = aux(0,1);
  gm(1,6) = aux(1,0);
  gm(1,7) = aux(1,1);
  gm(2,6) = aux(2,0);
  gm(2,7) = aux(2,1);

  cmulm(n2, transf, aux); // B_2 . T = N_2 * I . T = N_2 * T

  gm(0,2) = -aux(0,0);
  gm(0,3) = -aux(0,1);
  gm(1,2) = -aux(1,0);
  gm(1,3) = -aux(1,1);
  gm(2,2) = -aux(2,0);
  gm(2,3) = -aux(2,1);

  gm(0,4) = aux(0,0);
  gm(0,5) = aux(0,1);
  gm(1,4) = aux(1,0);
  gm(1,5) = aux(1,1);
  gm(2,4) = aux(2,0);
  gm(2,5) = aux(2,1);

  jac=l/2;
}



/**
  The function computes actual stresses at integration points on element.

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri,ci - row and column indices

  @return The function stores resulting stresses in the stress array 
          of the element integration points.

  Created by Tomas Koudelka, 11.2008
*/
void axisymlqinterface::compute_nlstress (long /*lcid*/, long eid, long ri, long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for(i=0; i<intordsm[0][0]; i++)
  {
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->computenlstresses(ipp, Mm->ip[ipp]);
    ipp++;
  }
}



/**
  The function computes actual stresses for the nonlocal material models at integration points on element.

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of the required integration point block on the given element
  @param[in] ci - column index of the required integration point block on the given element

  @return The function stores resulting stresses in the stress array 
          of the element integration points.

  Created by Tomas Koudelka, 1.2024
*/
void axisymlqinterface::compute_nonloc_nlstress (long /*lcid*/, long eid, long ri, long ci)
{
  long i,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for(i=0; i<intordsm[0][0]; i++)
  {
    //  computation of correct stresses
    if (Mp->strcomp==1)
      Mm->compnonloc_nlstresses(ipp);
    ipp++;
  }
}



/**
  The function computes increments of stresses at element integration points.

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of integration point block
  @param[in] ci - column indices of integration point block
   
  TKo  2.4.2024
*/
void axisymlqinterface::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long i, j ,ipp;

  ipp = Mt->elements[eid].ipp[ri+0][ci+0];
  for (i=0; i<intordsm[0][0]; i++){
    for (j=0; j<intordsm[0][0]; j++){
      //  computation of correct increments of stresses
      if (Mp->strcomp == 1)
        Mm->computenlstressesincr(ipp);
      ipp++;
    }
  }
}



/**
  The function computes eigen stresses caused by eigenstrains at all element integration points.
   
  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of the integration point block
  @param[in] ci - column index of the integration point block
   
  TKo 2.4.2024
*/
void axisymlqinterface::compute_eigstress(long /*lcid*/, long eid, long ri, long ci)
{
  long i,j,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ipp = Mt->elements[eid].ipp[ri+0][ci+0];
  for(i=0; i<intordsm[0][0]; i++){
    for(j=0; j<intordsm[0][0]; j++){
      //
      //  eigenstresses are computed from eigenstrains \sigma_0 = D (-\eps_0)
      //
      // restore eigenstrains
      Mm->giveeigstrain(ipp, eigstr);
      // change sign of eigenstrain vector
      chsgnv(eigstr);
      //  matrix of stiffness of the material
      Mm->matstiff(d, ipp);
      // calculate eigenstresses    
      mxv(d, eigstr, sig);
      Mm->storeeigstress(ipp, sig);
      ipp++;
    }
  }
}



/**
  The function computes internal forces due to actual stresses in the global coordinate system.
   
  @param[in]  lcid  - number of load case
  @param[in]  eid   - element id
  @param[in]  ri    - row index of the required integration point block on the given element
  @param[in]  ci    - column index of the required integration point block on the given element
  @param[out] ifor  - %vector of internal forces
  @param[in]  x     - %vector of x nodal coordinates
  @param[in]  y     - %vector of y nodal coordinates
   
  @return The function returns %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void axisymlqinterface::internal_forces (long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
  The function computes internal forces due to actual stresses in the global coordinate system.
  Nonlocal material models are considered in the stress computation.
   
  @param[in]  lcid  - number of load case
  @param[in]  eid   - element id
  @param[in]  ri    - row index of the required integration point block on the given element
  @param[in]  ci    - column index of the required integration point block on the given element
  @param[out] ifor  - %vector of internal forces
  @param[in]  x     - %vector of x nodal coordinates
  @param[in]  y     - %vector of y nodal coordinates
   
  @return The function returns %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 1.2024
*/
void axisymlqinterface::nonloc_internal_forces (long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
  The function computes resulting internal forces from actual stresses. If required, the transformation 
  to the nodal coordinate system is performed.
   
  @param[in]  lcid - number of load case
  @param[in]  eid  - element id
  @param[out] ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void axisymlqinterface::res_internal_forces (long lcid, long eid, vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d(x, y, eid);

  internal_forces(lcid, eid, 0, 0, ifor, x, y);
  
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
  The function computes resulting internal forces from actual stresses in the case of nonlocal material models. 
  If required, the transformation to the nodal coordinate system is performed.
   
  @param[in]  lcid - number of load case
  @param[in]  eid  - element id
  @param[out] ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 11.2008
*/
void axisymlqinterface::res_nonloc_internal_forces (long lcid, long eid, vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d(x, y, eid);

  nonloc_internal_forces(lcid, eid, 0, 0, ifor, x, y);
  
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
   function computes increment of  internal forces (from correct stresses increment)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void axisymlqinterface::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=stressincr;
  
  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes resulting increments of internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void axisymlqinterface::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
   
   this function is used in plane stress/strain elements (function is called
   by function res_eigstrain_forces) and shell elements

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   TKo 7.2008
*/
void axisymlqinterface::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - array containing nodal forces
   
   TKo 7.2008
*/
void axisymlqinterface::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
   The function integrates selected quantity on the selected finite element.
   It results in nodal values.
   
   @param[in] iq - type of integrated quantity (see alias.h)
   @param[in] lcid - number of load case
   @param[in] eid - element id
   @param[in] ri,ci - row and column indices
   @param[out] nv - nodal values
   @param[in] x,y - nodal coordinates

   @return The function returns nodal values of the integrated quantity in the parameter nv.
   
   JF, 29.10.2012
*/
void axisymlqinterface::elem_integration (integratedquant iq, long lcid, long eid, long ri, long ci,
                                      vector &nv, vector &x, vector &y)
{
  long ipp,i;
  double xi,jac,r;
  vector ipv(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordsm[0][0]));
  matrix gm(ASTCKMAT(tncomp, ndofe));
  vector gp(ASTCKVEC(intordsm[0][0]));

  
  nullv(nv);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[0][0];i++)
  {
    xi  =gp[i];
    //  function assembles required quantity at integration point
    Mm->givequantity(iq,lcid,ipp,0,ipv);
    //  strain-displacement (geometric) matrix
    geom_matrix(gm,x,y,xi,jac);
    //  contribution to the internal forces
    mtxv(gm,ipv,contr);
    //  radius in the given integration point
    r = approx (xi,x);
    cmulv(jac*w[i]*r,contr);
    //  summation
    addv(contr,nv,nv);
    ipp++;
  }
}
