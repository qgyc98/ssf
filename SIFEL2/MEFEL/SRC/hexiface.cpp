#include "hexiface.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globmat.h"
#include "element.h"
#include "node.h"
#include "intp.h"
#include "basefun.h"
#include "difcalc.h"
#include "vector.h"
#include "matrix.h"
#include "ordering.h"
#include "math.h"



// initialization of the array with the order/names of the strain components
const mechquant hexinterface::stra_comp_ord[] = {eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy};



hexinterface::hexinterface(void)
{
  long i;
  
  //  the numner of nodes
  nne=8;
  //  the number of DOFs
  ndofe=24;
  //  the number of blocks
  nb=1;
  //  the number of components
  tncomp=3;
  
  //  stress-strain state
  ssst = planecontact;

  ncomp = new long [nb];
  ncomp[0]=3;

  //  cumulative number of components approximated
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
  tnip=4;
}



hexinterface::~hexinterface(void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;

  delete [] ncomp;
  delete [] cncomp;
}



/**
  The function approximates value at the given point with the help of 
  its nodal values.

  @param[in] xi, eta - natrual coordinates of the required point
  @param[in] nodval - nodal values of the approximated quantity
  
  @return The function returns approximated value at the given point.

  Created by JK. 10.2023
 */
double hexinterface::approx(double xi, double eta, vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  nodval[4]= 0.0;
  nodval[5]= 0.0;
  nodval[6]= 0.0;
  nodval[7]= 0.0;

  scprd (bf,nodval,f);

  return f;
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
void hexinterface::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,ii,jj,kk;
  double xi, eta;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), w(ASTCKVEC(intordsm[ri][ci])), gp(ASTCKVEC(intordsm[ri][ci]));
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  kk=Mt->elements[eid].ipp[ri][ci];
  
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      gauss_points (gp.a, w.a, intordsm[ii][jj]);
      for (i = 0; i < intordsm[ii][jj]; i++){
        xi = gp[i];
        for (j = 0; j < intordsm[ii][jj]; j++){
          eta = gp[j];
          if (kk==ipp){
            coord[0]=approx (xi,eta,x);
            coord[1]=approx (xi,eta,y);
            coord[2]=approx (xi,eta,z);
            return;
          }
          kk++;
        }
      }
    }
  }
}



/**
  The function computes integration point values from the given nodal 
  values of selected quantity.
 
  @param[in] eid    - number of element
  @param[in] nodval - %vector of nodal values of the given quantity
  @param[out] ipval - %vector of integration point values approximated form 
                      the nodval.

  @return The function returns calculated int. point values in the argument ipval. 

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::intpointval(long eid, vector &nodval, vector &ipval)
{
  long ii, jj, i, j, k;
  double xi, eta;
  vector w,gp;

  k=0;

  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      gauss_points (gp.a, w.a, intordsm[ii][jj]);
      for (i = 0; i < intordsm[ii][jj]; i++){
        xi = gp[i];
        for (j = 0; j < intordsm[ii][jj]; j++){
          eta = gp[j];
          //  value in integration point
          ipval[k] = approx (xi, eta, nodval);
          k++;
        }
      }
    }
  }
}




/**
  The function assembles transformation %matrix between coordinate system
  defined on element and global coordinate system.
  Transformation deals with coordinates.
   
  The transformation %matrix is from local element system to global system x_g = T x_e.
   
  @param[out] tran - transformation %matrix from local element system to the global coordinate system
  @param[in] gx,gy,gz - coordinates of nodes in the global coordinate system
   
  JK, 8. 7. 2023
*/
void hexinterface::coord_transf_matrix(matrix &tran, vector &gx, vector &gy, vector &gz)
{
  double norm,alpha,zero=1.0e-8;
  vector bv1(ASTCKVEC(3)),bv2(ASTCKVEC(3)),bv3(ASTCKVEC(3));
  
  //  first nonscaled basis vector (in the xi direction)
  bv1[0] = 0.5*(gx[7]+gx[3]-gx[2]-gx[6]);
  bv1[1] = 0.5*(gy[7]+gy[3]-gy[2]-gy[6]);
  bv1[2] = 0.5*(gz[7]+gz[3]-gz[2]-gz[6]);
  
  //  norm of the first basis vector
  norm = sqrt(bv1[0]*bv1[0]+bv1[1]*bv1[1]+bv1[2]*bv1[2]);
  
  if (norm<zero){
    print_err("nonpositive length of the first basis vector in hexahedral interface element", __FILE__, __LINE__, __func__);
  }
  
  //  first scaled basis vector (in the xi direction)
  bv1[0] = bv1[0]/norm;
  bv1[1] = bv1[1]/norm;
  bv1[2] = bv1[2]/norm;
  
  
  
  //  second nonscaled basis vector (in the eta direction)
  bv2[0] = 0.5*(gx[1]+gx[5]-gx[2]-gx[6]);
  bv2[1] = 0.5*(gy[1]+gy[5]-gy[2]-gy[6]);
  bv2[2] = 0.5*(gz[1]+gz[5]-gz[2]-gz[6]);
  
  alpha = bv1[0]*bv2[0] + bv1[1]*bv2[1] + bv1[2]*bv2[2];
  //  second vector orthogonal to the first vector
  bv2[0] = bv2[0] - alpha*bv1[0];
  bv2[1] = bv2[1] - alpha*bv1[1];
  bv2[2] = bv2[2] - alpha*bv1[2];

  //  norm of the second basis vector
  norm = sqrt(bv2[0]*bv2[0]+bv2[1]*bv2[1]+bv2[2]*bv2[2]);
  
  if (norm<zero){
    print_err("nonpositive length of the second basis vector in hexahedral interface element", __FILE__, __LINE__, __func__);
  }
  
  //  second scaled basis vector (in the eta direction)
  bv2[0] = bv2[0]/norm;
  bv2[1] = bv2[1]/norm;
  bv2[2] = bv2[2]/norm;
  

  //  third basis vector
  bv3[0] = bv1[1]*bv2[2] - bv1[2]*bv2[1];
  bv3[1] = bv1[2]*bv2[0] - bv1[0]*bv2[2];
  bv3[2] = bv1[0]*bv2[1] - bv1[1]*bv2[0];
  
  
  //  transformation matrix from the element coordinate system to the global coordinate system
  tran[0][0] = bv1[0];
  tran[1][0] = bv1[1];
  tran[2][0] = bv1[2];

  tran[0][1] = bv2[0];
  tran[1][1] = bv2[1];
  tran[2][1] = bv2[2];

  tran[0][2] = bv3[0];
  tran[1][2] = bv3[1];
  tran[2][2] = bv3[2];
 
}



/**
  The function assembles transformation %matrix from local nodal coordinate
  system to the global coordinate system x_g = T x_l
   
  @param nodes - nodes of element
  @param tmat - transformation %matrix
   
  JK, 8. 7. 2023
*/
void hexinterface::transf_matrix(ivector &nodes, matrix &tmat)
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
   The function transforms node coordinates from the global coordinate system to
   the element coordinate system.
   
   Z coordinates have to be equal to zero for all element nodes after transformation.
   
   @param[in] gx, gy, gz - arrays with node coordinates in the global system
   @param[out] lx, ly - arrays with node coordinates in the element coordinate system
   @param[in] zero - computer zero

   JK, 22. 1. 2024
*/
void hexinterface::local_coordinates (vector &gx,vector &gy,vector &gz,vector &lx,vector &ly,double zero)
{
  vector lv(ASTCKVEC(3)),gv(ASTCKVEC(3));
  matrix tmat (ASTCKMAT(3,3));
  
  //  assembling of transformation matrix
  coord_transf_matrix(tmat,gx,gy,gz);

  //  transformation of coordinates of the first node
  gv[0]=gx[0]-gx[2];  gv[1]=gy[0]-gy[2];  gv[2]=gz[0]-gz[2];
  mtxv (tmat,gv,lv);
  lx[0]=lv[0];  ly[0]=lv[1];
  if (fabs(lv[2])>1.0e-8)
    print_err("Nonzero z coordinate of the first node in the element coordinate system.",__FILE__,__LINE__,__func__);

  //  transformation of coordinates of the second node
  gv[0]=gx[1]-gx[2];  gv[1]=gy[1]-gy[2];  gv[2]=gz[1]-gz[2];
  mtxv (tmat,gv,lv);
  lx[1]=lv[0];  ly[1]=lv[1];
  if (fabs(lv[2])>1.0e-8)
    print_err("Nonzero z coordinate of the second node in the element coordinate system.",__FILE__,__LINE__,__func__);

  //  transformation of coordinates of the third node
  //  the third node will be the origin of the coordinate system
  //  the transformation matrix is therefore multiplied by zero vector
  //  no transformation is needed
  lx[2]=0.0;  ly[2]=0.0;
  
  //  transformation of coordinates of the fourth node
  gv[0]=gx[3]-gx[2];  gv[1]=gy[3]-gy[2];  gv[2]=gz[3]-gz[2];
  mtxv (tmat,gv,lv);
  lx[3]=lv[0];  ly[3]=lv[1];
  if (fabs(lv[2])>1.0e-8){
    print_err("Nonzero z coordinate of the fourth node in the element coordinate system.", __FILE__, __LINE__, __func__);
    abort();
  }
}



/**
   The function assembles strain-displacement (geometric) %matrix.
   
   @param[out] gm - resulting geometric %matrix (output)
   @param[in] gx,gy,gz - nodal coordinates in the global system
   @param[in] xi,eta -  natural coordinate of required integration point
   @param[out] jac - Jacobian of transformation (output)
   
   @return The function returns resulting geometric %matrix in the parameter gm.

   JK, 8. 7. 2023
*/
void hexinterface::geom_matrix(matrix &gm, vector &gx, vector &gy, vector &gz,
                               double xi, double eta, double &jac)
{
  long i;
  matrix transf(ASTCKMAT(3,3));
  matrix aux(ASTCKMAT(3,3));
  vector bf(ASTCKVEC(nne));
  vector lx(ASTCKVEC(4)),ly(ASTCKVEC(4));
  
  vector dx(ASTCKVEC(4)),dy(ASTCKVEC(4));
  double zero = 1.0e-8;
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  
  local_coordinates (gx,gy,gz,lx,ly,zero);

  derivatives_2d (dx,dy,jac,lx,ly,xi,eta);

  fillm(0.0,gm);
  fillm(0.0,transf);
  fillm(0.0,aux);
  
  //  transformation matrix from element coordinate system to the global system
  coord_transf_matrix (transf,gx,gy,gz);
  
  //  basis functions on a quadrilateral element
  bf_lin_4_2d (bf.a,xi,eta);
  
  for (i=0;i<4;i++){
    
    //  B_i T
    cmulm(bf[i], transf, aux);
    
    //  fisrt half of the geometric matrix
    gm[0][0+i*3] = aux[0][0];
    gm[0][1+i*3] = aux[1][0];
    gm[0][2+i*3] = aux[2][0];

    gm[1][0+i*3] = aux[0][1];
    gm[1][1+i*3] = aux[1][1];
    gm[1][2+i*3] = aux[2][1];

    gm[2][0+i*3] = aux[0][2];
    gm[2][1+i*3] = aux[1][2];
    gm[2][2+i*3] = aux[2][2];
    
    //  second half of the geometrix matrix
    gm[0][0+i*3+12] = 0.0 - aux[0][0];
    gm[0][1+i*3+12] = 0.0 - aux[1][0];
    gm[0][2+i*3+12] = 0.0 - aux[2][0];

    gm[1][0+i*3+12] = 0.0 - aux[0][1];
    gm[1][1+i*3+12] = 0.0 - aux[1][1];
    gm[1][2+i*3+12] = 0.0 - aux[2][1];

    gm[2][0+i*3+12] = 0.0 - aux[0][2];
    gm[2][1+i*3+12] = 0.0 - aux[1][2];
    gm[2][2+i*3+12] = 0.0 - aux[2][2];
  }
    
}



/**
   The function computes stiffness %matrix of hexahedral
   interface finite element
   
   @param[in] eid - number of element
   @param[out] sm - resulting stiffness %matrix

   @return The function returns resulting stiffnes %matrix in the parameter sm.
   
   JK, 8. 7. 2023
*/
void hexinterface::stiffness_matrix(long eid, matrix &sm)
{
  long ii, jj, i, j, ipp;
  double xi, eta, jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp[0], ndofe)), d(ASTCKMAT(tncomp,tncomp));
  vector gp, w;

  nullm (sm);

  Mt->give_node_coord3d (x,y,z,eid);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (ii = 0; ii < nb; ii++){
    for (jj = 0; jj < nb; jj++){
      if (intordsm[ii][jj] == 0)
        continue;
      
      reallocv(RSTCKVEC(intordsm[ii][jj], gp));
      reallocv(RSTCKVEC(intordsm[ii][jj], w));
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ii][jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        for (j=0;j<intordsm[ii][jj];j++){
          eta=gp[j];
          //  geometric matrix
          geom_matrix (gm,x,y,z,xi,eta,jac);
          //  stiffness matrix of the material
          Mm->matstiff (d,ipp);
          // wieghts of the integration point
          jac *= w[i]*w[j];      
          //  contribution to the stiffness matrix of the element due to one integration point
          bdbjac (sm,gm,d,gm,jac);
          ipp++;
        }
      }
    }
  }
}



/**
  Function computes stiffness %matrix of one element. If it is required, nodal values 
  are transformed to the local nodal coordinate systems.

  @param eid - number of element
  @param sm - stiffness %matrix (output)

  @return The function returns required stiffness %matrix in the parameter sm.

  JK, 8. 7. 2023
*/
void hexinterface::res_stiffness_matrix(long eid, matrix &sm)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  
  
  stiffness_matrix(eid, sm);
  
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

   TKo, 10.2023
*/
void hexinterface::res_mainip_strains(long lcid, long eid)
{
  vector aux, r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;

  //  nodal displacements
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems(nodes);
  if (transf>0) {
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  mainip_strains(lcid, eid, r);
}



/**
  The function computes strains at integration points of element.
  It is used in geometrically linear problems.
   
  @param[in] lcid - load case id
  @param[in] eid - element id
  @param[in] r - %vector of nodal displacements
   
  @return The function stores resulting strains in the array strains 
          of the element integration points.
           
   TKo, 10.2023
*/
void hexinterface::mainip_strains(long lcid, long eid, vector &r)
{
  long ii, i, j, ipp;
  double xi, eta, jac;
  ivector nodes(ASTCKIVEC(nne));
  vector w, gp;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  matrix gm;
  vector eps;

  nullm(gm);

  Mt->give_node_coord3d(x, y, z, eid);
  Mt->give_elemnodes (eid, nodes);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (ii=0;ii<nb;ii++){
    reallocv(RSTCKVEC(intordsm[ii][ii], gp));
    reallocv(RSTCKVEC(intordsm[ii][ii], w));
    reallocv(RSTCKVEC(ncomp[ii], eps));
    reallocm(RSTCKMAT(ncomp[ii], ndofe, gm));
    
    gauss_points(gp.a, w.a, intordsm[ii][ii]);
    
    ipp = Mt->elements[eid].ipp[ii][ii];
    for(i=0; i<intordsm[ii][ii]; i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
        geom_matrix(gm, x, y, z, xi, eta, jac);
        mxv(gm, r, eps);
        Mm->storestrain(lcid, ipp, cncomp[ii], eps);
        ipp++;
      }
    }
  }
}



/**
  The function computes stresses at element integration points
   
  @param[in] lcid - load case id
  @param[in] eid - element id
   
  @return The function stores resulting stresses in the stress array 
          of the element integration points.
           
  Created by TKo, 10.2023
*/
void hexinterface::res_mainip_stresses(long lcid, long eid)
{
  long ri,ci;
  ri=0;
  ci=0;
  compute_nlstress (lcid,eid,ri,ci);
}



/**
  The function computes strains in nodes of element, nodal values are taken 
  from the closest integration point.
   
  @param[in] lcid - load case id
  @param[in] eid - element id
   
  JK, 10.5.2002
*/
void hexinterface::nod_strains_ip (long lcid, long eid)
{
  long i,j,ipp;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne)), auxiv;
  vector eps(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[0][0];
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  // repeat ip num extraction for remaining element nodes
  makerefv(auxiv, ipnum.a+4, 4);
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  
  
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
  The function computes nodal strains directly.
   
  @param[in] lcid - load case id
  @param[in] eid - element id
  @param[in] stra - array for strain components
   
  stra[i][j] - the j-th strain component at the i-th node
   
  JK, 26.9.2004
*/
void hexinterface::nod_strains_comp (long lcid,long eid,double **stra)
{
  long i,j;
  double jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), nxi(ASTCKVEC(nne)), neta(ASTCKVEC(nne));
  vector r(ASTCKVEC(ndofe)), eps(ASTCKVEC(tncomp)), aux;
  matrix tmat, gm(ASTCKMAT(tncomp,ndofe));
  
  //  node coordinates
  Mt->give_node_coord3d(x, y, z, eid);
  //  node numbers
  Mt->give_elemnodes(eid, nodes);
  //  nodal displacements
  eldispl(lcid, eid, r.a);
  
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv(RSTCKVEC(ndofe, aux));
    reallocm(RSTCKMAT(ndofe, ndofe, tmat));
    transf_matrix(nodes, tmat);
    lgvectortransf(aux, r, tmat);
    copyv(aux, r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_planelq(nxi, neta);
  
  //  loop over nodes
  for (i=0; i<nne; i++){
    //  geometric matrix
    geom_matrix(gm, x, y, z, nxi[i], neta[i], jac);
    //  strain computation
    mxv(gm, r, eps);
    
    for (j=0; j<eps.n; j++){
      stra[i][j] = eps[j];
    }
  }

}



/**
  The function computes stresses at nodes, nodal values are taken from the closest integration point.
   
  @param[in] lcid - load case id
  @param[in] eid - element id
   
  10.5.2002, JK
*/
void hexinterface::nod_stresses_ip (long lcid,long eid)
{
  long i,j,ipp;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne)), auxiv;
  vector sig(ASTCKVEC(tncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[0][0];
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  // repeat ip num extraction for remaining element nodes
  makerefv(auxiv, ipnum.a+4, 4);
  nodip_planelq(ipp, intordsm[0][0], ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid, nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress(lcid, ipnum[i], sig);
    //  storage of stresses to the node
    j=nod[i];
    Mt->nodes[j].storestress(lcid, 0, sig);
  }
}



/**
  Function computes actual stresses at integration points on element.

  @param[in] lcid - number of load case
  @param[in] eid  - element id
  @param[in] ri - row index of the required integration point block on the given element
  @param[in] ci - column index of the required integration point block on the given element

  @return The function stores resulting stresses in the stress array 
          of the element integration points.

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::compute_nlstress(long lcid, long eid, long ri, long ci)
{
  long ii, jj, i, j, ipp;
  
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  for (ii=0; ii<nb; ii++){
    for (jj=0; jj<nb; jj++){
      for(i=0; i<intordsm[ii][jj]; i++){
        for (j=0; j<intordsm[ii][jj]; j++){
          //  computation of actual stresses
          if (Mp->strcomp==1)
            Mm->computenlstresses(ipp, Mm->ip[ipp]);
          
          ipp++;
        }
      }
    }
  }
}



/**
  Function computes actual stresses for nonlocal material models at integration points on element.

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of the required integration point block on the given element
  @param[in] ci - column index of the required integration point block on the given element

  @return The function stores resulting stresses in the stress array 
          of the element integration points.

  Created by Tomas Koudelka, 1.2024
*/
void hexinterface::compute_nonloc_nlstress(long lcid, long eid, long ri, long ci)
{
  long ii, jj, i, j, ipp;
  
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  for (ii=0; ii<nb; ii++){
    for (jj=0; jj<nb; jj++){
      for(i=0; i<intordsm[ii][jj]; i++){
        for (j=0; j<intordsm[ii][jj]; j++){
          //  computation of actual stresses
          if (Mp->strcomp==1)
            Mm->compnonloc_nlstresses(ipp);
          
          ipp++;
        }
      }
    }
  }
}



/**
  The function computes resulting internal forces due to actual stresses.

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of the required integration point block on the given element
  @param[in] ci - column index of the required integration point block on the given element
  @param[out] ifor - %vector of internal forces
  @param[in] x - array of the x coordinates of element nodes
  @param[in] y - array of the y coordinates of element nodes
  @param[in] z - array of the z coordinates of element nodes
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::internal_forces(long lcid, long eid, long ri, long ci,
                                   vector &ifor, vector &x, vector &y, vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress(lcid, eid, ri, ci);
  
  //  integration of stresses over the element
  elem_integration(iq, lcid, eid, ri, ci, ifor, x, y, z);
}



/**
  The function computes resulting internal forces due to actual stresses in the case of nonlocal material models. 

  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri - row index of the required integration point block on the given element
  @param[in] ci - column index of the required integration point block on the given element
  @param[out] ifor - %vector of internal forces
  @param[in] x - array of the x coordinates of element nodes
  @param[in] y - array of the y coordinates of element nodes
  @param[in] z - array of the z coordinates of element nodes
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::nonloc_internal_forces(long lcid, long eid, long ri, long ci,
                                   vector &ifor, vector &x, vector &y, vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nonloc_nlstress(lcid, eid, ri, ci);
  
  //  integration of stresses over the element
  elem_integration(iq, lcid, eid, ri, ci, ifor, x, y, z);
}



/**
  The function computes resulting internal forces from actual stresses. If required, the transformation 
  to the nodal coordinate system is performed.
   
  @param[in]  lcid - number of load case
  @param[in]  eid  - element id
  @param[out] ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d(x, y, z, eid);

  internal_forces(lcid, eid, 0, 0, ifor, x, y, z);
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ndofe, ndofe);
    transf_matrix(nodes, tmat);
    glvectortransf(ifor, v, tmat);
    copyv(v, ifor);
  }
}



/**
  The function computes resulting internal forces from actual stresses in the case of nonlocal material models. 
  If required, the transformation to the nodal coordinate system is performed.
   
  @param[in]  lcid - number of load case
  @param[in]  eid  - element id
  @param[out] ifor - %vector of internal forces
   
  @return The function returns resulting %vector of internal forces in the parameter ifor.

  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes(ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d(x, y, z, eid);

  nonloc_internal_forces(lcid, eid, 0, 0, ifor, x, y, z);

  
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes(eid, nodes);
  transf = Mt->locsystems(nodes);
  if (transf>0){
    matrix tmat(ndofe, ndofe);
    transf_matrix(nodes, tmat);
    glvectortransf(ifor, v, tmat);
    copyv(v, ifor);
  }

  matrix sm(ASTCKMAT(ndofe, ndofe));
  vector r(ASTCKVEC(ndofe)), cifor(ASTCKVEC(ndofe)), difor(ASTCKVEC(ndofe));
  eldispl(lcid, eid, r.a);
  res_stiffness_matrix(eid, sm);
  mxv(sm, r, cifor);
  subv(ifor, cifor, difor);
  double nif = normv(ifor);
  double ncif = normv(cifor);
  double ndif = normv(difor);
  double iferr = ndif;
  if (nif > 1.0e-15)
    iferr /= nif;
  if (iferr > 1.0e-7){
    //fprintf(stdout, "\ndifference in internal force computation detected on element %ld, ndif=%le, ncif=%le, nif=%le\n", eid+1, ndif, ncif, nif);
    fprintf(Out, "\ndifference in internal force computation detected on element %ld, ndif=%le, ncif=%le, nif=%le, iferr=%le", eid+1, ndif, ncif, nif, iferr);
    fflush(Out);
  }
}



/**
  The function integrates selected quantity on the selected finite element, i.e. \int B^T q dV.
  It results in nodal values.
   
  @param[in] iq - type of integrated quantity (see alias.h)
  @param[in] lcid - number of load case
  @param[in] eid - element id
  @param[in] ri,ci - row and column indices
  @param[out] nv - nodal values
  @param[in] x,y,z - nodal coordinates

  @return The function returns nodal values of the integrated quantity in the parameter nv.
   
  Created by Tomas Koudelka, 10.2023
*/
void hexinterface::elem_integration(integratedquant iq, long lcid, long eid, long ri, long ci,
                                    vector &nv, vector &x, vector &y, vector &z)
{
  long ii, jj, i, j, ipp;
  double xi, eta, jac;
  vector w, gp, ipv(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));
  
  nullv (nv);
  
  ipp=Mt->elements[eid].ipp[ri][ci];

  for (ii=0; ii<nb; ii++){
    for (jj=0; jj<nb; jj++){
    
      if (intordsm[ii][jj] == 0)
        continue;
      
      reallocv(RSTCKVEC(intordsm[ii][jj],gp));
      reallocv(RSTCKVEC(intordsm[ii][jj],w));  
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
  
      for (i=0;i<intordsm[ii][jj];i++){
        xi=gp[i];
        for (j=0;j<intordsm[ii][jj];j++){
          eta=gp[j];
          //  function assembles required quantity at integration point
          Mm->givequantity(iq, lcid, ipp, cncomp[ii], ipv);
          //  strain-displacement (geometric) matrix
          geom_matrix(gm, x, y, z, xi, eta, jac);
          //  contribution to the nodal values
          mtxv(gm, ipv, contr);
          cmulv(jac*w[i]*w[j], contr);
          //  summation
          addv(contr, nv, nv);
          ipp++;
        }
      }
    }
  }
}


