#include "plquadiface.h"
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

plquadinterface::plquadinterface(void)
{
  long i;

  nne=4;  ndofe=8;
  nb=1;  tncomp=2;

  ssst = planecontact;

  ncomp = new long [nb];
  ncomp[0]=2;

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



plquadinterface::~plquadinterface(void)
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

void plquadinterface::eleminit (long eid)
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
double plquadinterface::approx(double xi, vector &nodval)
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
void plquadinterface::intpointval(long eid, vector &nodval, vector &ipval)
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
  The function assembles transformation %matrix x_g = T x_l
  for transformation from the local nodal coordinate systems.
   
  @param[in] nodes - array containing node numbers
  @param[out] tmat - transformation %matrix
   
  JK, 11.6.2006
*/
void plquadinterface::transf_matrix(ivector &nod, matrix &tmat)
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
void plquadinterface::stiffness_matrix(long eid, long ri, long ci, matrix &sm)
{
  long i,ipp;
  double xi,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector w(ASTCKVEC(intordsm[0][0])), gp(ASTCKVEC(intordsm[0][0]));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp[0], ndofe)), d(ASTCKMAT(tncomp,tncomp));

  nullm(sm);

  Mt->give_node_coord2d(x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  gauss_points (gp.a,w.a,intordsm[0][0]);

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++)
  {
    xi=gp[i];

    geom_matrix (gm, x, y, xi, jac);
    //  matrix of stiffness of the material
    Mm->matstiff (d,ipp);
    //  thickness in integration point
    thick = approx (xi,t);
    jac*=w[i]*thick;
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
void plquadinterface::res_stiffness_matrix(long eid, matrix &sm)
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
void plquadinterface::res_mainip_strains (long lcid, long eid)
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
void plquadinterface::mainip_strains (long lcid, long eid, long ri, long ci, vector &r)
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
   The function assembles strain-displacement (geometric) %matrix.
   
   @param[out] gm - resulting geometric %matrix
   @param[in] x,y - nodal coordinates
   @param[in] xi -  natural coordinate of required integration point
   @param[out] jac - Jacobian of transformation
   
   @return The function returns resulting geometric %matrix in the parameter gm.

   JF, 15.10.2012
   Modified by TKo, 16.12.2015
*/

void plquadinterface::geom_matrix (matrix &gm, vector &x, vector &y, double xi, double &jac)
{
  matrix transf(ASTCKMAT(2,2));
  matrix pom(ASTCKMAT(2,2));
  double l,dx,dy;

  dx = x[1] - x[0];
  dy = y[1] - y[0];

  nullm(gm);
  nullm(transf);
  nullm(pom);

  double n1 = 0.5*(1-xi);
  double n2 = 0.5*(1+xi);

  l=sqrt(dx*dx+dy*dy);

  transf[0][0]=dx/l;
  transf[0][1]=dy/l;
  transf[1][0]=-dy/l;
  transf[1][1]=dx/l;

  // node 1 coincides with node 4
  // node 2 coincides with node 3
  // The local x-axis corresponds to the direction from node 1 to node 2.
  // The local y-axis is normal to the local x-axis, i.e. local axis x rotated counterclockwise by pi/2.
  // Nodes 1 and 4 in local element ordering have natural coordinate ksi=-1.
  // Nodes 2 and 3 in local element ordering have natural coordinate ksi=+1.
  // Ordering of the relative displacements: [[u^l]], [[v^l]], i.e.:
  // slip in the local x direction,
  // relative displacement in the normal direction (local y direction).

  cmulm(n1, transf, pom); // B_1 . T = N_1 * I . T = N_1 * T

  gm[0][0]=-pom[0][0];
  gm[0][1]=-pom[0][1];
  gm[1][0]=-pom[1][0];
  gm[1][1]=-pom[1][1];

  gm[0][6]=pom[0][0];
  gm[0][7]=pom[0][1];
  gm[1][6]=pom[1][0];
  gm[1][7]=pom[1][1];


  cmulm(n2, transf, pom); // B_2 . T = N_2 * I . T = N_2 * T

  gm[0][2]=-pom[0][0];
  gm[0][3]=-pom[0][1];
  gm[1][2]=-pom[1][0];
  gm[1][3]=-pom[1][1];

  gm[0][4]=pom[0][0];
  gm[0][5]=pom[0][1];
  gm[1][4]=pom[1][0];
  gm[1][5]=pom[1][1];


  jac=l/2;
}



/**
  The function computes stresses at element integration points
   
  @param[in] lcid - load case id
  @param[in] eid - element id
   
  @return The function stores resulting stresses in the stress array 
          of the element integration points.
           
  Created by TKo, 11.2012
*/
void plquadinterface::res_mainip_stresses (long lcid, long eid)
{
  long ri, ci;
  ri = 0;
  ci = 0;
  compute_nlstress(lcid, eid, ri, ci);
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
void plquadinterface::compute_nlstress (long /*lcid*/, long eid, long ri, long ci)
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
void plquadinterface::compute_nonloc_nlstress (long /*lcid*/, long eid, long ri, long ci)
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
void plquadinterface::internal_forces (long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y)
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
void plquadinterface::nonloc_internal_forces (long lcid, long eid, long ri, long ci, vector &ifor, vector &x, vector &y)
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
void plquadinterface::res_internal_forces (long lcid, long eid, vector &ifor)
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
void plquadinterface::res_nonloc_internal_forces (long lcid, long eid, vector &ifor)
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
void plquadinterface::elem_integration (integratedquant iq, long lcid, long eid, long ri, long ci,
                                      vector &nv, vector &x, vector &y)
{
  long ipp,i;
  double xi,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector ipv(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordsm[0][0]));
  matrix gm(ASTCKMAT(tncomp, ndofe));
  vector gp(ASTCKVEC(intordsm[0][0]));

  
  nullv(nv);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);

  ipp=Mt->elements[eid].ipp[ri][ci];
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  for (i=0;i<intordsm[0][0];i++)
  {
    xi  =gp[i];
    //  function assembles required quantity at integration point
    Mm->givequantity(iq,lcid,ipp,0,ipv);
    //  strain-displacement (geometric) matrix
    geom_matrix(gm,x,y,xi,jac);
    //  contribution to the internal forces
    mtxv(gm,ipv,contr);
    //  thickness in integration point
    thick = approx (xi,t);
    cmulv(jac*w[i]*thick,contr);
    //  summation
    addv(contr,nv,nv);
    ipp++;
  }
}
