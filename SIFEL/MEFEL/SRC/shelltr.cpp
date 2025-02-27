#include "shelltr.h"
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
#include "plelemrotlt.h"
#include "dkt.h"
#include "intpoints.h"
#include <math.h>

shelltr::shelltr (void)
{
  long i,j;
  
  if (Perlt==NULL)   Perlt = new planeelemrotlt;
  if (Dkt==NULL)     Dkt = new dktelem;

  // ******
  //  NNE
  // ******
  //  number nodes on element
  nne=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=1;
  //  number of nodes on surface
  nnsurf=4;

  // **********************
  //  strain/stress state
  // **********************
  ssst=shell;

  // ********
  //  NDOFE
  // ********
  //  number of DOFs on plane element
  ndofes = Perlt->ndofe;
  //  number of DOFs on plate element
  ndofep = Dkt->ndofe;
  
  //  number of DOFs on shell element
  ndofe=ndofes+ndofep;
  
  
  // *********
  //  NAPFUN
  // *********
  //  number of functions approximated
  napfuns = Perlt->napfun;
  napfunp = Dkt->napfun;
  
  napfun=napfuns+napfunp;
  
  // *******************
  //  number of blocks
  // *******************
  nbs = Perlt->nb;
  nbp = Dkt->nb;
  
  nb = nbs + nbp;
  
  // ***********************************************
  //  number of components in strain/stress blocks
  // ***********************************************

  ncomps = new long [nbs];
  ncomps[0]=Perlt->ncomp[0];  //  =2
  ncomps[1]=Perlt->ncomp[1];  //  =1
 
  ncompp = new long [nbp];
  ncompp[0]=Dkt->ncomp[0];  //  =3


  // ************************************************
  //  the number of integration points
  // ************************************************
  nip = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i]=new long [nb];
  }
  
  for (i=0;i<nbs;i++){
    for (j=0;j<nbs;j++){
      nip[i][j]=Perlt->nip[i][j];
    }
  }
  for (i=0;i<nbp;i++){
    for (j=0;j<nbp;j++){
      nip[i+nbs][j+nbs]=Dkt->nip[i][j];
    }
  }

  // *****************************************
  //  the total number of integration points
  // *****************************************
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }


  // ************************************************
  //  the order of integration of stiffness matrix
  // ************************************************
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    intordsm[i] = new long [nb];
  }
  
  for (i=0;i<nbs;i++){
    for (j=0;j<nbs;j++){
      intordsm[i][j]=Perlt->intordsm[i][j];
    }
  }
  for (i=0;i<nbp;i++){
    for (j=0;j<nbp;j++){
      intordsm[i+nbs][j+nbs]=Dkt->intordsm[i][j];
    }
  }
  
  // *******************************
  //  the number of all components
  // *******************************
  tncomps=Perlt->tncomp;
  tncompp=Dkt->tncomp;
  
  //tncomp = tncomps+tncompp;
  tncomp = 6;

  // ********************************************
  //  ordering of unknowns on the shell element
  // ********************************************
  
  ordering = new long* [2];
  ordering[0] = new long [ndofes];
  ordering[1] = new long [ndofep];

  //  unknowns in the first node
  //  displacement in the x direction
  ordering[0][0]=1;
  //  displacement in the y direction
  ordering[0][1]=2;
  //  rotation around the z direction
  ordering[0][2]=6;
  
  //  displacement in the z direction
  ordering[1][0]=3;
  //  rotation around the x axis
  ordering[1][1]=4;
  //  rotation around the y axis
  ordering[1][2]=5;
  

  //  unknowns in the second node
  //  displacement in the x direction
  ordering[0][3]=7;
  //  displacement in the y direction
  ordering[0][4]=8;
  //  rotation around the z direction
  ordering[0][5]=12;
  
  //  displacement in the z direction
  ordering[1][3]=9;
  //  rotation around the x axis
  ordering[1][4]=10;
  //  rotation around the y axis
  ordering[1][5]=11;
  
  
  //  unknowns in the third node
  //  displacement in the x direction
  ordering[0][6]=13;
  //  displacement in the y direction
  ordering[0][7]=14;
  //  rotation around the z direction
  ordering[0][8]=18;
  
  //  displacement in the z direction
  ordering[1][6]=15;
  //  rotation around the x axis
  ordering[1][7]=16;
  //  rotation around the y axis
  ordering[1][8]=17;

}

shelltr::~shelltr (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;

  delete [] ncomps;
  delete [] ncompp;

  for (i=0;i<2;i++){
    delete [] ordering[i];
  }
  delete [] ordering;
}

/**
   function assembles element code numbers
   they are used for localization between all element unknowns and unknowns related to one subelement
   
   @param cn - array of code numbers
   @param ri - subelement id
   
   JK, 24. 8. 2018
*/
void shelltr::codnum (long *cn,long ri)
{
  long i;
  
  if (ri==0){
    for (i=0;i<ndofes;i++){
      cn[i]=ordering[0][i];
    }
  }
  if (ri==1){
    for (i=0;i<ndofep;i++){
      cn[i]=ordering[1][i];
    }
  }
}

/**
   function assembles transformation %matrix between coordinate system
   defined on element and global coordinate system
   transformation deals with coordinates
   
   transformation %matrix from local element system to global system
   x_g = T x_e
   
   @param tran - tranformation %matrix T
   @param gx,gy,gz - vectors of nodal coordinates in the global coordinate systems
   @param zero - computer zero
   
   JK, 1.10.2008
*/
void shelltr::coord_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz,double zero)
{
  double c,p,l;
  
  //  first direction vector is located in the first column
  tran[0][0]=gx[1]-gx[0];
  tran[1][0]=gy[1]-gy[0];
  tran[2][0]=gz[1]-gz[0];

  //  length of the first direction vector
  l=sqrt(tran[0][0]*tran[0][0]+tran[1][0]*tran[1][0]+tran[2][0]*tran[2][0]);
  if (l<zero){
    print_err("nonpositive length of the first basis vector in triangular shell element", __FILE__, __LINE__, __func__);
  }
  
  //  normalization of the first direction vector
  tran[0][0]/=l;
  tran[1][0]/=l;
  tran[2][0]/=l;


  //  second direction vector is located in the second column
  tran[0][1]=gx[2]-gx[0];
  tran[1][1]=gy[2]-gy[0];
  tran[2][1]=gz[2]-gz[0];
  
  //  dot product of the first and second direction vector
  p=tran[0][0]*tran[0][1]+tran[1][0]*tran[1][1]+tran[2][0]*tran[2][1];
  //  coefficient in Gramm-Schmidt orthonormalization process
  c=0.0-p;
  
  //  orthogonalization of the second vector with respect to the first vector
  tran[0][1]=tran[0][1]+tran[0][0]*c;
  tran[1][1]=tran[1][1]+tran[1][0]*c;
  tran[2][1]=tran[2][1]+tran[2][0]*c;
  
  //  length of the second direction vector
  l=sqrt(tran[0][1]*tran[0][1]+tran[1][1]*tran[1][1]+tran[2][1]*tran[2][1]);
  if (l<zero){
    print_err("nonpositive length of the second basis vector in triangular shell element", __FILE__, __LINE__, __func__);
  }

  //  normalization of the second direction vector
  tran[0][1]/=l;
  tran[1][1]/=l;
  tran[2][1]/=l;

  
  //  third direction vector is obtained from vector product of the first and second vector
  tran[0][2]=tran[1][0]*tran[2][1]-tran[2][0]*tran[1][1];
  tran[1][2]=tran[2][0]*tran[0][1]-tran[0][0]*tran[2][1];
  tran[2][2]=tran[0][0]*tran[1][1]-tran[1][0]*tran[0][1];
  
  //  length of the third vector (it should be one, it is calculated in order to improve numerical stability)
  l=sqrt(tran[0][2]*tran[0][2]+tran[1][2]*tran[1][2]+tran[2][2]*tran[2][2]);
  if (fabs(l-1.0)>zero){
    print_err("the size of the third vector %le is not equal to 1 triangular shell element", __FILE__, __LINE__, __func__,l);
  }
  
  //  normalization of the third vector
  tran[0][2]/=l;
  tran[1][2]/=l;
  tran[2][2]/=l;
  
}


/**
   function transforms node coordinates from the global coordinate system to
   the element coordinate system
   
   z coordinates have to be equal to zero after transformation
   
   @param gx, gy, gz - arrays with node coordinates in the global system
   @param lx, ly - arrays with node coordinates in the element coordinate system
   @param zero - computer zero

   JK, 17. 3. 2020
*/
void shelltr::local_coordinates (vector &gx,vector &gy,vector &gz,vector &lx,vector &ly,double zero)
{
  vector lv(ASTCKVEC(3)),gv(ASTCKVEC(3));
  matrix tmat (ASTCKMAT(3,3));
  
  //  assembling of transformation matrix
  coord_transf_matrix (tmat,gx,gy,gz,zero);

  //  transformation of coordinates of the first node
  //  the first node will be the origin of the coordinate system
  //  the transformation matrix is therefore multiplied by zero vector
  //  no transformation is needed

  //  transformation of coordinates of the second node
  gv[0]=gx[1]-gx[0];  gv[1]=gy[1]-gy[0];  gv[2]=gz[1]-gz[0];
  mtxv (tmat,gv,lv);
  lx[1]=lv[0];  ly[1]=lv[1];
  if (fabs(lv[2])>1.0e-8)
    print_err("Nonzero z coordinate of the second node in the element coordinate system.",__FILE__,__LINE__,__func__);

  //  transformation of coordinates of the third node
  gv[0]=gx[2]-gx[0];  gv[1]=gy[2]-gy[0];  gv[2]=gz[2]-gz[0];
  mtxv (tmat,gv,lv);
  lx[2]=lv[0];  ly[2]=lv[1];
  if (fabs(lv[2])>1.0e-8)
    print_err("Nonzero z coordinate of the third node in the element coordinate system.",__FILE__,__LINE__,__func__);
  
}


/**
   function assembles transformation %matrix between global coordinate system
   and local nodal coordinate system
   transformation deals with nodal forces

   transformation %matrix x_g = T x_l
   
   @param nodes - array of nodes on element
   @param tmat - transformation %matrix
   
   9.5.2002
*/
void shelltr::node_transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,i6,n,m;

  fillm (0.0,tmat);

  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      i6=i*6;
      tmat[i6][i6]   = Mt->nodes[nodes[i]].e1[0];  tmat[i6][i6+1]   = Mt->nodes[nodes[i]].e2[0]; tmat[i6][i6+2]     = Mt->nodes[nodes[i]].e3[0];
      tmat[i6+1][i6] = Mt->nodes[nodes[i]].e1[1];  tmat[i6+1][i6+1] = Mt->nodes[nodes[i]].e2[1]; tmat[i6+1][i6+2]   = Mt->nodes[nodes[i]].e3[1];
      tmat[i6+2][i6] = Mt->nodes[nodes[i]].e1[2];  tmat[i6+2][i6+1] = Mt->nodes[nodes[i]].e2[2]; tmat[i6+2][i6+2]   = Mt->nodes[nodes[i]].e3[2];
      i6=i*6+3;
      tmat[i6][i6]   = Mt->nodes[nodes[i]].e1[0];  tmat[i6][i6+1]   = Mt->nodes[nodes[i]].e2[0]; tmat[i6][i6+2]     = Mt->nodes[nodes[i]].e3[0];
      tmat[i6+1][i6] = Mt->nodes[nodes[i]].e1[1];  tmat[i6+1][i6+1] = Mt->nodes[nodes[i]].e2[1]; tmat[i6+1][i6+2]   = Mt->nodes[nodes[i]].e3[1];
      tmat[i6+2][i6] = Mt->nodes[nodes[i]].e1[2];  tmat[i6+2][i6+1] = Mt->nodes[nodes[i]].e2[2]; tmat[i6+2][i6+2]   = Mt->nodes[nodes[i]].e3[2];
    }
  }
}

/**
   function assembles transformation %matrix for the whole vectors/matrices
   from the local element coordinate system to the global coordinate system
   
   nodal unknowns - u, v, w, dw/dx, dw/dy, phi_z,
   
   @param gx, gy, gz - arrays of node coordinates in the global coordinate system
   @param tran - transformation %matrix (it contains 32 components)
   @param zero - computer zero
   
   JK, 17. 3. 2020
*/
void shelltr::elem_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz,double zero)
{
  long i,ri,ci;
  matrix tm (ASTCKMAT(3,3));
  
  //  transformation matrix with 3 rows and columns
  coord_transf_matrix (tm,gx,gy,gz,zero);
  
  //  all entries of the transformation matrix are set to zero
  fillm (0.0,tran);
  
  //  diagonal entries of the transformation matrix are set to one
  for (i=0;i<ndofe;i++){
    tran[i][i]=1.0;
  }
  
  //  the following transformation is connected with Perlq and Qkirch elements
  for (i=0;i<nne;i++){
    //  transformation of the displacements
    ri=i*6;  ci=i*6;
    tran[ri+0][ci+0]=tm[0][0];  tran[ri+0][ci+1]=tm[0][1];  tran[ri+0][ci+2]=tm[0][2];
    tran[ri+1][ci+0]=tm[1][0];  tran[ri+1][ci+1]=tm[1][1];  tran[ri+1][ci+2]=tm[1][2];
    tran[ri+2][ci+0]=tm[2][0];  tran[ri+2][ci+1]=tm[2][1];  tran[ri+2][ci+2]=tm[2][2];
    
    //  transformation of the rotations
    ri=i*6+3;  ci=i*6+3;
    tran[ri+0][ci+0]=tm[0][0];  tran[ri+0][ci+1]=tm[0][1];  tran[ri+0][ci+2]=tm[0][2];
    tran[ri+1][ci+0]=tm[1][0];  tran[ri+1][ci+1]=tm[1][1];  tran[ri+1][ci+2]=tm[1][2];
    tran[ri+2][ci+0]=tm[2][0];  tran[ri+2][ci+1]=tm[2][1];  tran[ri+2][ci+2]=tm[2][2];
    
  }
}

/**
   function assembles stiffness %matrix of triangular shell element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   4.5.2002
*/
void shelltr::res_stiffness_matrix (long eid,matrix &sm)
{
  long *cn,transf;
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne)),gz(ASTCKVEC(nne)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  ivector nodes(nne);
  matrix lsm,tmat(ASTCKMAT(ndofe,ndofe)),auxsm(ASTCKMAT(ndofe,ndofe));
  
  //  nodes on element
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates in the global system
  Mt->give_node_coord3d (gx,gy,gz,eid);
  //  transformation of the node coordinates in the global system gx, gy, gz to the element coordinate system x, y
  local_coordinates (gx,gy,gz,x,y,Mp->zero);
  
  fillm (0.0,sm);
  
  cn = new long [ndofes];
  reallocm (ndofes,ndofes,lsm);
  
  //  contribution from triangular plane element with rotational degrees of freedom
  Perlt->stiffness_matrix (eid,0,0,lsm,x,y);
  
  //  code numbers for plane stress matrix localization to the shell matrix
  codnum (cn,0);
  
  //  localization of the plane stress matrix into the shell matrix
  mat_localize (sm,lsm,cn,cn);
  
  delete [] cn;
  
  cn = new long [ndofep];
  reallocm (ndofep,ndofep,lsm);
  
  //  contribution from constant curve triangular plate element
  Dkt->stiffness_matrix (eid,2,2,lsm,x,y);
  
  //  code numbers for plate matrix into the shell matrix
  codnum (cn,1);

  //  localization of the plate matrix into the shell matrix
  mat_localize (sm,lsm,cn,cn);

  delete [] cn;
  
  //  transformation of stiffness matrix from element coordinate system to the global system
  //  transformation matrix
  elem_transf_matrix (tmat,gx,gy,gz,Mp->zero);
  mxm (tmat,sm,auxsm);
  mxmt (auxsm,tmat,sm);
  
  
  //  transformation of stiffness matrix to nodal coordinate system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    node_transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  
  
}



/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 23.9.2008
*/
void shelltr::res_ip_strains (long lcid,long eid)
{
  vector aux(ASTCKVEC(ndofe));
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne)),gz(ASTCKVEC(nne)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),rs(ASTCKVEC(ndofes)),rp(ASTCKVEC(ndofep));
  ivector cn,nodes(ASTCKIVEC(nne));
  matrix tmat(ASTCKMAT(ndofe,ndofe));
  
  //  node coordinates in global system
  Mt->give_node_coord3d (gx,gy,gz,eid);
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  transformation of the node coordinates in the global system gx, gy, gz to the element coordinate system x, y
  local_coordinates (gx,gy,gz,x,y,Mp->zero);

  //  nodal displacements and rotations
  eldispl (lcid,eid,r.a);
  
  //  transformation of displacement vector from nodal system to global system
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    node_transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  transformation of nodal values in global system to element system
  elem_transf_matrix (tmat,gx,gy,gz,Mp->zero);
  glvectortransf (r, aux, tmat);
  copyv (aux,r);
  
  reallocv (RSTCKIVEC(ndofes,cn));
  //  code numbers for plane stress components
  codnum (cn.a,0);

  globloc (r.a,rs.a,cn.a,ndofes);
  
  Perlt->ip_strains (lcid,eid,0,0,x,y,rs);
  
  /*
  vector eps(6);
  long ri=0;
  long ci=0;
  long ipp=Mt->elements[eid].ipp[ri][ci];
  for (long i=0;i<6;i++){
    eps[0]=Mm->ip[ipp].strain[0];
    eps[1]=Mm->ip[ipp].strain[1];
    eps[5]=Mm->ip[ipp].strain[2];
    loc_glob_tens_trans (eps,tran);
    Mm->ip[ipp].strain[0]=eps[0];
    Mm->ip[ipp].strain[1]=eps[1];
    Mm->ip[ipp].strain[2]=eps[5];
    ipp++;
  }
  */

  reallocv (RSTCKIVEC(ndofep,cn));
  codnum (cn.a,1);
  
  globloc (r.a,rp.a,cn.a,ndofep);
  
  Dkt->ip_strains (lcid,eid,2,2,x,y,rp);

}


/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void shelltr::res_ip_stresses (long lcid,long eid)
{
  Perlt->compute_nlstress (lcid,eid,0,0);
  Dkt->compute_nlstress (lcid,eid,2,2);
}

/**
   function computes stresses
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void shelltr::compute_nlstress (long lcid,long eid)
{
  Perlt->compute_nlstress (lcid,eid,0,0);
  Dkt->compute_nlstress (lcid,eid,2,2);
}


/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   
   JK, 23.9.2008
*/
void shelltr::compute_nlstressincr (long lcid,long eid)
{
  Perlt->compute_nlstressincr (lcid,eid,0,0);
  Dkt->compute_nlstressincr (lcid,eid,2,2);
}


/**
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 26.9.2008
*/
void shelltr::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long *cn,transf;
  ivector nodes (nne);
  vector aux(ASTCKVEC(ndofe));
  vector gx(ASTCKVEC(nne)),gy(ASTCKVEC(nne)),gz(ASTCKVEC(nne)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),v;
  matrix tmat(ASTCKMAT(ndofe,ndofe));
  
  //  node coordinates in the global system
  Mt->give_node_coord3d (gx,gy,gz,eid);
  //  transformation of the node coordinates in the global system gx, gy, gz to the element coordinate system x, y
  local_coordinates (gx,gy,gz,x,y,Mp->zero);
  
  reallocv (ndofes,v);
  
  Perlt->internal_forces (lcid,eid,0,0,v,x,y);
  
  cn = new long [ndofes];
  //  code numbers for plane stress matrix localization to the shell matrix
  codnum (cn,0);
  locglob (ifor.a,v.a,cn,ndofes);
  
  delete [] cn;
  cn = new long [ndofep];
  reallocv (ndofep,v);
  
  Dkt->internal_forces (lcid,eid,2,2,v,x,y);
  
  //  code numbers for plane stress matrix localization to the shell matrix
  codnum (cn,1);
  locglob (ifor.a,v.a,cn,ndofep);
  
  delete [] cn;
  
  //  transformation of forces in element coordinate system into the global system
  elem_transf_matrix (tmat,gx,gy,gz,Mp->zero);
  lgvectortransf (aux, ifor, tmat);
  copyv (aux,ifor);
  
  
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    node_transf_matrix (nodes,tmat);
    glvectortransf (ifor,v,tmat);
    copyv (v,ifor);
  }
}





void shelltr::nod_strains_ip (long lcid,long eid)
{
  
  //  contribution from triangular plane element with rotational degrees of freedom
  //Perlt->nod_strains_ip (lcid,eid,0,0);
  //  contribution from triangular plate element
  //Dkt->nod_strains_ip (lcid,eid,2,2);

  long i,j,ipp;
  ivector ipnums(nne),ipnump(nne),nod(nne);
  vector gx(nne),gy(nne),gz(nne),x(nne),y(nne);
  vector epss(tncomps),epsp(tncomps),eps(tncomp);
  matrix tran(3,3);

  //  node coordinates in global coordinate system
  Mt->give_node_coord3d (gx,gy,gz,eid);
  //  transformation from global to element coordinate system
  //coord_transf_matrix (x,y, tran, gx,  gy,  gz);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[0][0];
  nodip_planelt (ipp,intordsm[0][0],ipnums);
  ipp=Mt->elements[eid].ipp[2][2];
  nodip_planelt (ipp,intordsm[2][2],ipnump);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain (lcid,ipnums[i],epss);
    Mm->givestrain (lcid,ipnump[i],epsp);
    
    /*
    eps[0]=epss[0];
    eps[1]=epss[1];
    eps[2]=0.0;
    eps[3]=epsp[0];
    eps[4]=epsp[0];
    eps[5]=epss[0];
    */

    loc_glob_tens_trans (eps,tran);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storestrain (lcid,0,eps.n,eps);
  }
  
  /*
  reallocv (tncompd,eps);
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[2][2];
  nodip_planelt (ipp,intordsm[2][2],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    Mm->givestrain (lcid,ipnum[i],eps);
    
    loc_glob_tens_trans (eps,tran);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storestrain (lcid,tncompd,eps.n,eps);
  }
  */
}
  
void shelltr::strains (long /*lcid*/,long /*eid*/)
{

}

void shelltr::nod_stresses_ip (long lcid,long eid)
{
  
  //  contribution from triangular plane element with rotational degrees of freedom
  Perlt->nod_stresses_ip (lcid,eid,0,0);
    
  //  contribution from triangular plate element
  Dkt->nod_stresses_ip (lcid,eid,2,2);
}


  
void shelltr::stresses (long /*lcid*/,long /*eid*/)
{
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
void shelltr::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta;
//  double ipval;
  vector w, gp1, gp2, anv(nne);
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
        reallocv (intordsm[ii][jj],w);
        gauss_points_tr (gp1.a, gp2.a, w.a, intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp1[k];
          eta=gp2[k];
          //  value in integration point
//          ipval = approx_nat (xi,eta,anv);
          ncompstr =  Mm->ip[ipp].ncompstr;
          ncompeqother = Mm->ip[ipp].ncompeqother;
          if ((ictn[0] & inistrain) && (j < ncompstr))
          {
            //Mm->ip[ipp].strain[idstra] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & inistress) && (j < nstra + ncompstr))
          {
            //Mm->ip[ipp].stress[idstre] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & iniother) && (j < nstra+nstre+ncompeqother))
          {
            //Mm->ip[ipp].eqother[idoth] += ipval;
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
            //Mm->ic[ipp][idic] += ipval;
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
