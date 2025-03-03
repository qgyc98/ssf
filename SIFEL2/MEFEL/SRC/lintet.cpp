#include <stdlib.h>
#include <math.h>
#include "lintet.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "gadaptivity.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "loadcase.h"
#include "mathem.h"




lintet::lintet (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=1;
  //  number of edges on element
  ned=6;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=4;
  //  number of nodes on one surface
  nnsurf=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;
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
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=1;

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=1;

}

lintet::~lintet (void)
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
   
   @param volcoord - volume coordinates
   @param nodval - nodal values
   
   20.8.2001
*/
double lintet::approx (vector &volcoord,vector &nodval)
{
  double f;

  scprd (volcoord,nodval,f);
  
  return f;
}

/**
   The function computes approximated function value with the help of nodal values.
   
   @param xi, eta, zeta - natural coordinates
   @param nodval - nodal values
   
   20.8.2001
*/
double lintet::approx_nat (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector volcoord(ASTCKVEC(4));
  
  volcoord[0]=xi;
  volcoord[1]=eta;
  volcoord[2]=zeta;
  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
  
  scprd (volcoord,nodval,f);
  
  return f;
}

/**
   function assembles %matrix of base functions

   @param n - %matrix of base functions
   @param volcoord - volume coordinates
   
   24.8.2001
*/
void lintet::bf_matrix (matrix &n,vector &volcoord)
{
  nullm (n);
  
  n[0][0]=volcoord[0];  n[0][3]=volcoord[1];  n[0][6]=volcoord[2];  n[0][9] =volcoord[3];
  n[1][1]=volcoord[0];  n[1][4]=volcoord[1];  n[1][7]=volcoord[2];  n[1][10]=volcoord[3];
  n[2][2]=volcoord[0];  n[2][5]=volcoord[1];  n[2][8]=volcoord[2];  n[2][11]=volcoord[3];
}
void lintet::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  nullm (n);
  
  n[0][0]=xi;  n[0][3]=eta;  n[0][6]=zeta;  n[0][9] =1.0-xi-eta-zeta;
  n[1][1]=xi;  n[1][4]=eta;  n[1][7]=zeta;  n[1][10]=1.0-xi-eta-zeta;
  n[2][2]=xi;  n[2][5]=eta;  n[2][8]=zeta;  n[2][11]=1.0-xi-eta-zeta;
}

/**
   function assembles geometric %matrix
   vector of strains has following ordering
   eps=(e_xx, e_yy, e_zz, e_yz, e_zx, e_xy)
   
   @param gm - geometric %matrix
   @param x,y,z - vectors of node coordinates
   
   24.8.2001
*/
void lintet::geom_matrix_old (matrix &gm,vector &x,vector &y,vector &z)
{
  long i,i1,i2,i3,i4,i5,i6,i7,i8,i9;
  double det;
  vector b(ASTCKVEC(4)),c(ASTCKVEC(4)),d(ASTCKVEC(4));
  
  det = det3d (x.a,y.a,z.a);
  
  volb_3d (b.a,y.a,z.a,det);
  volc_3d (c.a,x.a,z.a,det);
  vold_3d (d.a,x.a,y.a,det);
  i1=0;  i2=1;  i3=2;  i4=1;  i5=2;  i6=0;  i7=2;  i8=0;  i9=1;
  for (i=0;i<4;i++){
    gm[0][i1]=b[i];  i1+=3;
    gm[1][i2]=c[i];  i2+=3;
    gm[2][i3]=d[i];  i3+=3;
    
    gm[3][i4]=d[i];  i4+=3;
    gm[3][i5]=c[i];  i5+=3;
    
    gm[4][i6]=d[i];  i6+=3;
    gm[4][i7]=b[i];  i7+=3;
    
    gm[5][i8]=c[i];  i8+=3;
    gm[5][i9]=b[i];  i9+=3;
  }
  
}

/**
   function assembles geometric %matrix
   vector of strains has following ordering
   eps=(e_xx, e_yy, e_zz, e_yz, e_zx, e_xy)
   
   @param gm - geometric %matrix
   @param x,y,z - vectors of node coordinates
   
   24.8.2001
*/
void lintet::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,double &jac)
{
  long i,i1,i2,i3,i4,i5,i6,i7,i8,i9;
  vector dx(ASTCKVEC(4)),dy(ASTCKVEC(4)),dz(ASTCKVEC(4));
  
  dx_bf_lin_tet (dx.a);
  dy_bf_lin_tet (dy.a);
  dz_bf_lin_tet (dz.a);

  derivatives_3d (dx,dy,dz,jac,x,y,z,0.25,0.25,0.25);

  i1=0;  i2=1;  i3=2;  i4=1;  i5=2;  i6=0;  i7=2;  i8=0;  i9=1;
  for (i=0;i<4;i++){
    gm[0][i1]=dx[i];  i1+=3;
    gm[1][i2]=dy[i];  i2+=3;
    gm[2][i3]=dz[i];  i3+=3;
    
    gm[3][i4]=dz[i];  i4+=3;
    gm[3][i5]=dy[i];  i5+=3;
    
    gm[4][i6]=dz[i];  i6+=3;
    gm[4][i7]=dx[i];  i7+=3;
    
    gm[5][i8]=dy[i];  i8+=3;
    gm[5][i9]=dx[i];  i9+=3;
  }
  
}



/**
   function assembles transformation %matrix
   
   @param nodes - nodes of element
   @param tmat - transformation %matrix
   
*/
void lintet::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of one element

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   19.7.2001
*/
void lintet::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long ipp;
  double det,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  Mt->give_node_coord3d (x,y,z,eid);
  nullm (sm);
  det = det3d (x.a,y.a,z.a);
  //det = fabs(det);
  
  if(det < 0.0){
    //det = fabs(det);
    print_err("wrong numbering of nodes on 3D element number %ld, negative volume! det = %e", __FILE__, __LINE__, __func__, eid+1,det);
    abort();
  }

  ipp=Mt->elements[eid].ipp[ri][ci];

  geom_matrix (gm,x,y,z,jac);
  Mm->matstiff (d,ipp);

  jac=det/6.0;
  bdbjac (sm,gm,d,gm,jac);
}

/**
   function computes B^T D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - B^T D %matrix

   17. 2. 2014
*/
void lintet::bd_matrix (long eid,long ri,long ci,matrix &bd)
{
  long ipp;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),d(ASTCKMAT(tncomp,tncomp));

  Mt->give_node_coord3d (x,y,z,eid);
  nullm (bd);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  geometric matrix
  geom_matrix (gm,x,y,z,jac);
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  //  B^T D
  mtxm (gm,d,bd);
  //  B^T D Jac
  jac/=6.0;
  cmulm (jac,bd,bd);
}

/**
   function computes integral of D %matrix of one element
   this %matrix is used in homogenization methods with prescribed stresses
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param bd - integral of D %matrix

   17. 2. 2014
*/
void lintet::dd_matrix (long eid,long ri,long ci,matrix &dd)
{
  long ipp;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  matrix d(ASTCKMAT(tncomp,tncomp));

  Mt->give_node_coord3d (x,y,z,eid);
  jac = det3d (x.a,y.a,z.a)/6.0;
  nullm (dd);
  
  //  integration point id
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  //  D Jac
  cmulm (jac,d,dd);
}

/**
   function assembles resulting stiffness %matrix of the element

   @param eid - element id
   @param sm - stiffness %matrix

   JK, 9.5.2002
*/
void lintet::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);

  //  transformation of stiffness matrix
  ivector nodes (ASTCKIVEC(nne));
  long transf;
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes mass %matrix
   
   @param eid - number of element
   @param mm - mass %matrix
   
   26.8.2001
*/
void lintet::mass_matrix (long eid,matrix &mm)
{
  long i;
  double det,jac,rho;
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp1(ASTCKVEC(intordmm)),gp2(ASTCKVEC(intordmm));
  vector gp3(ASTCKVEC(intordmm)),volcoord(ASTCKVEC(4)),dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);
  Mt->give_node_coord3d (x,y,z,eid);

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordmm);

  det = det3d (x.a,y.a,z.a);
  //det = fabs(det);
  
  nullm (mm);
  
  for (i=0;i<intordmm;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];
    
    bf_matrix (n,volcoord);
    rho = approx (volcoord,dens);
    
    jac=det*w[i]*rho;
    
    nnj (mm.a,n.a,jac,n.m,n.n);
  }
  
}

/**
   function assembles mass %matrix of tetrahedral
   finite element with bilinear approximation functions
   
   @param eid - element id
   @param mm - mass %matrix
   
   JK, 16. 6. 2019
*/
void lintet::res_mass_matrix (long eid,matrix &mm)
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


/**
   function computes load %matrix
   
   @param eid - number of element
   @param lm - load %matrix

   26.8.2001
*/
void lintet::load_matrix (long eid,matrix &lm)
{
  long i;
  double det,jac;
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),w(ASTCKVEC(intordmm)),gp1(ASTCKVEC(intordmm));
  vector gp2(ASTCKVEC(intordmm)),gp3(ASTCKVEC(intordmm)),volcoord(ASTCKVEC(4));
  matrix n(ASTCKMAT(napfun,ndofe));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordmm);

  det = det3d (x.a,y.a,z.a);
  //det = fabs(det);
  
  nullm (lm);

  for (i=0;i<intordmm;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];

    bf_matrix (n,volcoord);

    jac=det*w[i];

    nnj (lm.a,n.a,jac,n.m,n.n);
  }
}



/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   TKo - 1.2009
*/
void lintet::res_ip_strains (long lcid,long eid)
{
  ip_strains (lcid,eid,0,0);
}



/**
   function computes strains in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void lintet::ip_strains (long lcid,long eid,long ri,long ci)
{
  long ipp;
  double jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),eps(ASTCKVEC(tncomp)),aux;
  ivector nodes(ASTCKIVEC(nne));
  matrix gm(ASTCKMAT(tncomp,ndofe)),tmat;

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
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

  ipp=Mt->elements[eid].ipp[ri][ci];
  
  geom_matrix (gm,x,y,z,jac);
  mxv (gm,r,eps);
  Mm->storestrain (lcid,ipp,eps);
  
}

/**
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   
   10.5.2002
*/
void lintet::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  double vol;
  vector eps(ASTCKVEC(tncomp));
  ivector nod(ASTCKIVEC(nne));
  

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];

  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  Mm->givestrain (lcid,ipp,eps);
  
  for (i=0;i<nne;i++){
    //  storage of strains to the node
    j=nod[i];
    if (Mp->strainaver==1)
      Mt->nodes[j].storestrain (lcid,0,eps);
    if (Mp->strainaver==2){
      vol=volumeip (eid,ri,ci);
      cmulv (vol,eps,eps);
      Mt->nodes[j].storestrain (lcid,0,vol,eps);
    }
  }
}



void lintet::strains (long lcid,long eid,long ri,long ci)
{
  long i;
  double **stra = NULL;
  vector coord,eps;
  
  if (Mp->strainaver==0){
    stra = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
    }
    //elem_strains (stra,lcid,eid,ri,ci);
  }

  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //  strains are computed at integration points
    res_ip_strains (lcid,eid);
    break;
  }
  case enodes:{
    //  strains are copied to nodes from the closest integration points
    nod_strains_ip (lcid,eid,ri,ci);
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
/*    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (RSTCKVEC(ncp,eps));
    reallocv (RSTCKVEC(3,coord));
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	//appstrain (lcid,eid,coord[0],coord[1],coord[2],0,ncp,eps);
      
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlt::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  if (Mp->strainaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }
}


void lintet::res_ip_stresses (long lcid,long eid)
{
  ip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses in main integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002
*/
void lintet::ip_stresses (long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
  
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   10.5.2002
*/
void lintet::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  double vol;
  vector nxi(ASTCKVEC(nne)),neta(ASTCKVEC(nne)),nzeta(ASTCKVEC(nne)),eps(ASTCKVEC(tncomp));
  vector epst, epstt, sig(ASTCKVEC(tncomp)), natcoord(ASTCKVEC(3));
  ivector nod(ASTCKIVEC(nne));
  
  Mt->give_elemnodes (eid,nod);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  stresses at the closest integration point
  Mm->givestress (lcid,ipp,sig);
  
  for (i=0;i<nne;i++){
    
    //  storage of stresses to the node
    j=nod[i];
    if (Mp->stressaver==1)
      Mt->nodes[j].storestress (lcid,0,sig);
    if (Mp->stressaver==2){
      vol=volumeip (eid,ri,ci);
      cmulv (vol,sig,sig);
      Mt->nodes[j].storestress (lcid,0,vol,sig);
    }
  }
  
  
}



void lintet::stresses (long lcid,long eid,long ri,long ci)
{
  long i;
  double **stra=NULL,**stre=NULL;
  vector coord,sig;
  
  if (Mp->stressaver==0){
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    //elem_strains (stra,lcid,eid,ri,ci);
    //elem_stresses (stra,stre,lcid,eid,ri,ci);
  }

  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    res_ip_stresses (lcid,eid);
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
    reallocv (RSTCKVEC(3,coord));
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);
      
      if (Mp->stressaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],coord[1],coord[2],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    */
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  if (Mp->stressaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
      delete [] stre[i];
    }
    delete [] stra;
    delete [] stre;
  }
}








/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void lintet::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes internal forces for nonlocal models

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void lintet::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes increments of internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   TKo 7.2008
*/
void lintet::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes nodal forces caused by temperature changes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y,z - nodal coordinates
   
   7.2008, TKo
*/
void lintet::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
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
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void lintet::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y,z);
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
void lintet::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y,z);
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
void lintet::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  incr_internal_forces (lcid,eid,0,0,ifor,x,y,z);
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
void lintet::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)),x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  eigstrain_forces (lcid,eid,0,0,nfor,x,y,z);
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
void lintet::compute_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
}



/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo  7.2008
*/
void lintet::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstressesincr (ipp);
}



/**
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void lintet::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void lintet::compute_nonloc_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;

  ipp=Mt->elements[eid].ipp[ri][ci];

  //  computation of correct stresses
  if (Mp->strcomp==1)
    Mm->compnonloc_nlstresses (ipp);
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void lintet::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ipp=Mt->elements[eid].ipp[ri][ci];
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
}



/**
   The function integrates selected stress type quantity over the finite element, i.e.
   it performs \int_{\Omega} \mbf{B}^T \mbf{\sigma} d\Omega. It results in nodal values.
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values (output)
   @param x,y,z - node coordinates
   
   @return The function returns nodal values calculated in the %vector nv

   TKo 7.2008
*/
void lintet::elem_integration(integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z)
{
  long ipp;
  double det,jac;
  vector ipv(ASTCKVEC(tncomp)),contr(ASTCKVEC(ndofe));
  matrix gm(ASTCKMAT(tncomp,ndofe));

  ipp=Mt->elements[eid].ipp[ri][ci];
  nullv (nv);
  det = det3d (x.a,y.a,z.a);  
  //det = fabs(det);

  //  function assembles required quantity at integration point
  Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);

  //  strain-displacement (geometric) matrix
  geom_matrix (gm,x,y,z,jac);
  //  contribution to the nodal values
  mtxv (gm,ipv,contr);
  jac=det/6.0;
  cmulv (jac,contr);
  //  summation
  addv(contr,nv,nv);
}



/**
   The function integrates arbitrary selected quantity over the finite element, e.g.
   it performs \int_{\Omega} \mbf{\sigma} d\Omega which results in integrated values that can 
   be used in the homogenization problems.

   
   @param eid - element id
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param iv - integrated values (output)
   
   @return The function returns nodal values calculated in the %vector nv

   Created by TKo, 23.1.2015
*/
void lintet::elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv)
{
  long ipp;
  double jac;
  vector ipv(ASTCKVEC(iv.n));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));

  ipp=Mt->elements[eid].ipp[0][0];
  Mt->give_node_coord3d (x,y,z,eid);

  nullv(iv);

  jac = det3d (x.a,y.a,z.a)/6.0;  
  //det = fabs(det);

  //  function assembles required quantity at integration point
  Mm->givequantity (iq, lcid, ipp, cncomp[0], ipv);

  cmulv (jac, ipv, iv);
}



/**
   The function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of integration point block (input)
   @param ci - column index of integration point block (input)
   @param coord - %vector with global coordinates of integration point (ouput)
   
   @return The function returns global coordinates in the argument coord.

   19.1.2002
*/
void lintet::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,ii;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), volcoord(ASTCKVEC(4)), w(ASTCKVEC(intordsm[ri][ci]));
  vector gp1(ASTCKVEC(intordsm[ri][ci])), gp2(ASTCKVEC(intordsm[ri][ci])), gp3(ASTCKVEC(intordsm[ri][ci]));

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];

  for (i=0;i<intordsm[ri][ci];i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];

    if (ii==ipp){
      coord[0]=approx (volcoord,x);
      coord[1]=approx (volcoord,y);
      coord[2]=approx (volcoord,z);
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
void lintet::ipncoord (long eid,long ipp,vector &ncoord)
{
  long i, ii, ri, ci;
  vector gp1, gp2, gp3, w;

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
      gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ri][ci]);
      ii=Mt->elements[eid].ipp[ri][ci];

      for (i=0;i<intordsm[ri][ci];i++){
        if (ii==ipp){
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
void lintet::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, zeta, ipval;
  vector w, gp1, gp2, gp3, anv(ASTCKVEC(nne));
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
        if (intordsm[ii][jj] == 0)
          continue;
        ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
        ncompstr = Mm->ip[ipp].ncompstr;
        ncompeqother = Mm->ip[ipp].ncompeqother;
        reallocv (RSTCKVEC(intordsm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordsm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordsm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordsm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordsm[ii][ii]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp1[k];
          eta=gp2[k];
          zeta=gp3[k];
          //  value in integration point
          ipval = approx_nat (xi,eta,zeta,anv);
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
   function computes volume appropriate to integration point
   
   2.3.2004, JK
*/
void lintet::ipvolume (long eid,long ri,long ci)
{
  long ipp;
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  det = det3d (x.a,y.a,z.a);
  //jac=fabs(det)/6.0;
  jac=det/6.0;
  
  Mm->storeipvol (ipp,jac);
}
/**
   function computes volume appropriate to integration point
   
   2.3.2004, JK
*/
double lintet::volumeip (long eid,long ri,long ci)
{
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  det = det3d (x.a,y.a,z.a);
  //jac=fabs(det)/6.0;
  jac=det/6.0;
  
  return jac;
}



void lintet::nod_other_ip (long eid)
{
  long i,j,ipp,ncompo;
  double vol;
  ivector nod(ASTCKIVEC(nne));
  vector other;
  
  //  node numbers of the element
  Mt->give_elemnodes(eid,nod);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  ncompo = Mm->givencompeqother(ipp,0);
  reallocv(RSTCKVEC(ncompo, other));
  Mm->giveother(ipp, 0, ncompo, other.a);
  
  for (i=0;i<nne;i++){
    //  storage of other to the node
    j=nod[i];
    if (Mp->otheraver==1)
      Mt->nodes[j].storeother(0, ncompo, other);
    if (Mp->otheraver==2){
      vol=volumeip (eid,0,0);
      cmulv(vol, other);
      Mt->nodes[j].storeother(0, ncompo, vol, other);
    }
  }
}



/**
   function computes nodal forces caused by presure on surface
   
   @param lcid - number of load case
   @param eid - element id
   @param is - identification of surfaces
   @param nv - nodal values
   @param nf - nodal forces
   
   5.4.2005, JK
   pracuje pouze v globalnim souradnem systemu
*/
void lintet::node_forces_surf_old (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i;
  double xi,eta,zeta,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), v(ASTCKVEC(ndofe)), av(ASTCKVEC(ndofe));
  vector gp1(ASTCKVEC(intordb)), gp2(ASTCKVEC(intordb)), w(ASTCKVEC(intordb));
  matrix am(ASTCKMAT(ndofe,ndofe)),n(ASTCKMAT(napfun,ndofe));
  
  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  if (is[0]==1){
    for (i=0;i<intordb;i++){
      xi=0.0; eta=gp1[i];  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,eta,zeta,0);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[6] = nv[0*3+0];
    av[7] = nv[0*3+1];
    av[8] = nv[0*3+2];
    
    av[3] = nv[1*3+0];
    av[4] = nv[1*3+1];
    av[5] = nv[1*3+2];
    
    av[9]  = nv[2*3+0];
    av[10] = nv[2*3+1];
    av[11] = nv[2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[1]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=0.0;  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,xi,zeta,1);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[0] = nv[9+0*3+0];
    av[1] = nv[9+0*3+1];
    av[2] = nv[9+0*3+2];
    
    av[6] = nv[9+1*3+0];
    av[7] = nv[9+1*3+1];
    av[8] = nv[9+1*3+2];
    
    av[9]  = nv[9+2*3+0];
    av[10] = nv[9+2*3+1];
    av[11] = nv[9+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[2]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=0.0;
      
      jac2d_3d (jac,x,y,z,xi,eta,2);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[3] = nv[18+0*3+0];
    av[4] = nv[18+0*3+1];
    av[5] = nv[18+0*3+2];
    
    av[0] = nv[18+1*3+0];
    av[1] = nv[18+1*3+1];
    av[2] = nv[18+1*3+2];
    
    av[9]  = nv[18+2*3+0];
    av[10] = nv[18+2*3+1];
    av[11] = nv[18+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[3]==1){
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=1.0-xi-eta;
      
      jac2d_3d (jac,x,y,z,xi,eta,3);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    av[3] = nv[27+0*3+0];
    av[4] = nv[27+0*3+1];
    av[5] = nv[27+0*3+2];
    
    av[6] = nv[27+1*3+0];
    av[7] = nv[27+1*3+1];
    av[8] = nv[27+1*3+2];
    
    av[0] = nv[27+2*3+0];
    av[1] = nv[27+2*3+1];
    av[2] = nv[27+2*3+2];
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
}


/**
   function computes nodal forces caused by presure on surface
   
   @param lcid - number of load case
   @param eid - element id
   @param is - identification of surfaces
   @param nv - nodal values
   @param nf - nodal forces
   
   12.12.2005, JK
*/
void lintet::node_forces_surf (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i,transf;
  double xi,eta,zeta,jac;
  double *tnv;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), v(ASTCKVEC(ndofe)), av(ASTCKVEC(ndofe));
  vector gp1(ASTCKVEC(intordb)), gp2(ASTCKVEC(intordb)), w(ASTCKVEC(intordb));
  matrix am(ASTCKMAT(ndofe,ndofe)),n(ASTCKMAT(napfun,ndofe)),tmat;
  
  tnv = new double [9];
  
  // nodes on element
  Mt->give_elemnodes(eid, nodes);
  //  coordinates of element nodes
  Mt->give_node_coord3d (x,y,z,eid);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  if (is[0]>0 ){
    nullm (am);
    for (i=0;i<intordb;i++){
      xi=0.0; eta=gp1[i];  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,eta,zeta,0);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    nullv (av);
    
    if (is[0]==1){
      av[6] = nv[0*3+0];
      av[7] = nv[0*3+1];
      av[8] = nv[0*3+2];
      
      av[3] = nv[1*3+0];
      av[4] = nv[1*3+1];
      av[5] = nv[1*3+2];
      
      av[9] = nv[2*3+0];
      av[10] = nv[2*3+1];
      av[11] = nv[2*3+2];
    }


    if (is[0]==2){
      
      av[0] = nv[0*3+0];
      av[1] = nv[0*3+1];
      av[2] = nv[0*3+2];
      
      av[3] = nv[1*3+0];
      av[4] = nv[1*3+1];
      av[5] = nv[1*3+2];
      
      av[6] = nv[2*3+0];
      av[7] = nv[2*3+1];
      av[8] = nv[2*3+2];
      
      locglob_nodeval (0,av,tnv,x,y,z);
      nullv (av);
      
      av[6] = tnv[0];
      av[7] = tnv[1];
      av[8] = tnv[2];
      
      av[3] = tnv[3];
      av[4] = tnv[4];
      av[5] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }


    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[1]>0){
    nullm (am);
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=0.0;  zeta=gp2[i];
      
      jac2d_3d (jac,x,y,z,xi,zeta,1);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[1]==1){
      av[0] = nv[9+0*3+0];
      av[1] = nv[9+0*3+1];
      av[2] = nv[9+0*3+2];
      
      av[6] = nv[9+1*3+0];
      av[7] = nv[9+1*3+1];
      av[8] = nv[9+1*3+2];
      
      av[9] = nv[9+2*3+0];
      av[10] = nv[9+2*3+1];
      av[11] = nv[9+2*3+2];
    }
    
    if (is[1]==2){
      
      av[0] = nv[9+0*3+0];
      av[1] = nv[9+0*3+1];
      av[2] = nv[9+0*3+2];
      
      av[3] = nv[9+1*3+0];
      av[4] = nv[9+1*3+1];
      av[5] = nv[9+1*3+2];
      
      av[6] = nv[9+2*3+0];
      av[7] = nv[9+2*3+1];
      av[8] = nv[9+2*3+2];
      
      locglob_nodeval (1,av,tnv,x,y,z);
      nullv (av);
      
      av[0] = tnv[0];
      av[1] = tnv[1];
      av[2] = tnv[2];
      
      av[6] = tnv[3];
      av[7] = tnv[4];
      av[8] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }
    
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[2]>0){
    nullm (am);
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=0.0;
      
      jac2d_3d (jac,x,y,z,xi,eta,2);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[2]==1){
      av[3] = nv[18+0*3+0];
      av[4] = nv[18+0*3+1];
      av[5] = nv[18+0*3+2];
      
      av[0] = nv[18+1*3+0];
      av[1] = nv[18+1*3+1];
      av[2] = nv[18+1*3+2];
      
      av[9] = nv[18+2*3+0];
      av[10] = nv[18+2*3+1];
      av[11] = nv[18+2*3+2];

    }

    if (is[2]==2){
      
      av[0] = nv[18+0*3+0];
      av[1] = nv[18+0*3+1];
      av[2] = nv[18+0*3+2];
      
      av[3] = nv[18+1*3+0];
      av[4] = nv[18+1*3+1];
      av[5] = nv[18+1*3+2];
      
      av[6] = nv[18+2*3+0];
      av[7] = nv[18+2*3+1];
      av[8] = nv[18+2*3+2];
      
      locglob_nodeval (2,av,tnv,x,y,z);
      nullv (av);
      
      av[3] = tnv[0];
      av[4] = tnv[1];
      av[5] = tnv[2];
      
      av[0] = tnv[3];
      av[1] = tnv[4];
      av[2] = tnv[5];
      
      av[9] = tnv[6];
      av[10] = tnv[7];
      av[11] = tnv[8];
    }

    mxv (am,av,v);
    addv (v,nf,nf);
  }
  if (is[3]>0){
    nullm (am);
    for (i=0;i<intordb;i++){
      xi=gp1[i]; eta=gp2[i];  zeta=1.0-xi-eta;
      
      jac2d_3d (jac,x,y,z,xi,eta,3);
      bf_matrix (n,xi,eta,zeta);
      jac = jac*w[i];
      nnj (am.a,n.a,jac,n.m,n.n);
    }
    
    if (is[3]==1){
      av[0] = nv[27+0*3+0];
      av[1] = nv[27+0*3+1];
      av[2] = nv[27+0*3+2];
      
      av[3] = nv[27+1*3+0];
      av[4] = nv[27+1*3+1];
      av[5] = nv[27+1*3+2];
      
      av[6] = nv[27+2*3+0];
      av[7] = nv[27+2*3+1];
      av[8] = nv[27+2*3+2];

    }

    if (is[3]==2){
      
      av[0] = nv[27+0*3+0];
      av[1] = nv[27+0*3+1];
      av[2] = nv[27+0*3+2];
      
      av[3] = nv[27+1*3+0];
      av[4] = nv[27+1*3+1];
      av[5] = nv[27+1*3+2];
      
      av[6] = nv[27+2*3+0];
      av[7] = nv[27+2*3+1];
      av[8] = nv[27+2*3+2];
      
      locglob_nodeval (3,av,tnv,x,y,z);
      nullv (av);
      
      av[0] = tnv[0];
      av[1] = tnv[1];
      av[2] = tnv[2];
      
      av[3] = tnv[3];
      av[4] = tnv[4];
      av[5] = tnv[5];
      
      av[6] = tnv[6];
      av[7] = tnv[7];
      av[8] = tnv[8];
    }
    
    nullv (v);
    mxv (am,av,v);
    addv (v,nf,nf);
  }
  // transformation to loacal coordinate systems of nodes 
  transf = Mt->locsystems (nodes);
  if (transf>0)
  {
    reallocm(RSTCKMAT(ndofe,ndofe,tmat));
    transf_matrix (nodes,tmat);
    lgvectortransfblock (nf,tmat);
  }
  
  delete [] tnv;
}


void lintet::locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z)
{
  double norm;
  vector g1(ASTCKVEC(3)), g2(ASTCKVEC(3)), g3(ASTCKVEC(3)), lv(ASTCKVEC(3)), gv(ASTCKVEC(3));
  matrix t(ASTCKMAT(3,3));
  
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
  
}

/**
   function interpolates the nodal values to the integration points on the element
   
   @param eid - element id
   
   JK, 22.4.2005
*/
void lintet::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  vector volcoord(ASTCKVEC(4));
  
  volcoord[0]=0.25;
  volcoord[1]=0.25;
  volcoord[2]=0.25;
  volcoord[3]=0.25;
  
  ipval[0]=approx (volcoord,nodval);
}


/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 22.4.2005
*/
void lintet::res_eigstrain_forces (long eid,vector &nfor)
{
  eigstrain_forces (eid,nfor);
}

/**
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - vector of internal forces

   JK, 22.4.2005
*/
void lintet::eigstrain_forces (long eid,vector &nfor)
{
  long i,ipp;
  double jac,det;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), eigstr(ASTCKVEC(tncomp)), sig(ASTCKVEC(tncomp)), contr(ASTCKVEC(ndofe));
  matrix d(ASTCKMAT(tncomp,tncomp)),gm(ASTCKMAT(tncomp,ndofe));
  
  Mt->give_node_coord3d (x,y,z,eid);
  det = det3d (x.a,y.a,z.a);
  //det = fabs(det);

  nullv (nfor);
  
  ipp=Mt->elements[eid].ipp[0][0];
  
  Mm->giveeigstrain (ipp,cncomp[0],ncomp[0],eigstr);
  
  //  matrix of stiffness of the material
  Mm->matstiff (d,ipp);
  
  mxv (d,eigstr,sig);
  
  geom_matrix (gm,x,y,z,jac);
  
  mtxv (gm,sig,contr);
  
  jac = det/6.0;

  cmulv (jac,contr);
  
  for (i=0;i<contr.n;i++){
    nfor[i]+=contr[i];
  }
  
}



/**
   function defines meaning of DOFs at nodes
   
   @param eid - element id
   
   TKr 11.4.2024 according to JK
*/
void lintet::define_meaning (long eid)
{
  ivector cn(ASTCKIVEC(ndofe)),nod(ASTCKIVEC(nne));
  
  Mt->give_elemnodes (eid,nod);
  Mt->give_code_numbers (eid,cn.a);

  //  displacement in x direction
  if (cn[0]>0)  Mt->nodes[nod[0]].meaning[0]=1;
  //  displacement in y direction
  if (cn[1]>0)  Mt->nodes[nod[0]].meaning[1]=2;
  //  displacement in z direction
  if (cn[2]>0)  Mt->nodes[nod[0]].meaning[2]=3;

  //  displacement in x direction
  if (cn[3]>0)  Mt->nodes[nod[1]].meaning[0]=1;
  //  displacement in y direction
  if (cn[4]>0)  Mt->nodes[nod[1]].meaning[1]=2;
  //  displacement in z direction
  if (cn[5]>0)  Mt->nodes[nod[1]].meaning[2]=3;

  //  displacement in x direction
  if (cn[6]>0)  Mt->nodes[nod[2]].meaning[0]=1;
  //  displacement in y direction
  if (cn[7]>0)  Mt->nodes[nod[2]].meaning[1]=2;
  //  displacement in z direction
  if (cn[8]>0)  Mt->nodes[nod[2]].meaning[2]=3;

  //  displacement in x direction
  if (cn[9]>0)  Mt->nodes[nod[3]].meaning[0]=1;
  //  displacement in y direction
  if (cn[10]>0)  Mt->nodes[nod[3]].meaning[1]=2;
  //  displacement in z direction
  if (cn[11]>0)  Mt->nodes[nod[3]].meaning[2]=3;

}



////////////////////       /* termitovo */       ////////////////////////////////////

/**
   Function integrates function "NT*D*B*r" (= NT*Stress) over whole element.
   N is %matrix of basic functions.
   !!! Values in 'ntdbr' are stored unusually. Not after this manner {(val[0]; ... ;val[tncomp])[0] ; ... ; (val[0]; ... ;val[tncomp])[nne]}
   but after this manner {(val[0]; ... ;val[nne])[0] ; ... ; (val[0]; ... ;val[nne])[ntncomp]}
   
   @param eid   - element id
   @param ntdbr - empty(returned) array, dimension is tncomp*nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void lintet::ntdbr_vector (long eid,vector &ntdbr)
{
  long i,j,ipp,ri,ci;
  double det,jac,w;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), volcoord(ASTCKVEC(4)), eps(ASTCKVEC(tncomp)), sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));

  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];

  Mt->give_node_coord3d (x,y,z,eid);

  fillv (0.25,volcoord);
  w = 1.0/6.0;

  det = det3d (x.a,y.a,z.a);
  //jac = w*fabs(det);
  jac = w*det;

  Mm->matstiff (d,ipp);
  Mm->givestrain (0,ipp,eps);
  mxv (d,eps,sig);

  for (i=0;i<tncomp;i++)
    for (j=0;j<nne;j++)
      ntdbr[j+i*nne] = jac * sig[i] * volcoord[j];
}

/**
   Function integrates function "NT*N" over whole element.
   N is %matrix of basic functions.
   !!! Matrix N is composed of three same blocks, it is computed only one third.
   
   @param eid - element id
   @param ntn - empty(returned) 2Darray, dimension is nne x nne
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void lintet::ntn_matrix (long eid,matrix &ntn)
{
  long i;
  double det,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  
  Mt->give_node_coord3d (x,y,z,eid);
  det = det3d (x.a,y.a,z.a);
  //jac = fabs(det)/120.0;
  jac = det/120.0;

  fillm (jac,ntn);

  for (i=0;i<nne;i++)
    ntn[i][i] *= 2.0;
}

/**
   Function 1.
   1. computes "e2" - square of energy norm of error of solution on element
      e2 = integral{ (sig_star - sig_roof)T * D_inv * (sig_star - sig_roof) }
      sig_star = recovered stress, values in nodes are defined by "rsigfull" and over element are interpolated by base functions
      sig_roof = stress obtained by FEM
      D_inv = inverse stiffness %matrix of material
   2. computes "u2" - square of energy norm of strain on element
      u2 = integral{ epsT * D * eps }
      eps = strain obtained by FEM
      D = stiffness %matrix of material
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
double lintet :: compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull)
{
  long intord = 4;
  long i,ipp,ri,ci;
  double det,jac,contr,zero=1.0e-20;
  ivector nodes (ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), volcoord(ASTCKVEC(4));
  vector gp1(ASTCKVEC(intord)), gp2(ASTCKVEC(intord)), gp3(ASTCKVEC(intord)), w(ASTCKVEC(intord));
  vector sig_star(ASTCKVEC(tncomp)), sig_roof(ASTCKVEC(tncomp)), sig_err(ASTCKVEC(tncomp)), eps(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp)),dinv(ASTCKMAT(tncomp,tncomp));
  
  ri = ci = 0;
  ipp = Mt->elements[eid].ipp[ri][ci];
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);

  det = det3d (x.a,y.a,z.a);
  //det = fabs(det);

  Mm->matstiff (d,ipp);
  invm (d,dinv,zero);
  
  // compute u2
  fillv (0.25,volcoord);

  jac = det;

  Mm->givestrain (0,ipp,eps);

  mxv (d,eps,sig_roof);
  scprd (eps,sig_roof,contr);

  u2 = jac*contr;

  // compute e2
  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intord);

  e2 = 0;  
  for (i=0;i<intord;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-gp1[i]-gp2[i]-gp3[i];

    jac=w[i]*det;

    give_der_star (volcoord,rsigfull,nodes,sig_star,Mt->nn);
    subv (sig_star,sig_roof,sig_err);
    
    vxmxv (sig_err,dinv,contr);
    e2 += jac*contr;
  }
  
  sizel = pow (1.4142135624*det,1.0/3.0);
  
  return det/6.0;
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
void lintet :: elchar (long eid, matrix &spsig)
{
  long ipp;
  vector eps(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  reallocm(RSTCKMAT(1, tncomp, spsig));
  
  ipp = Mt->elements[eid].ipp[0][0];
  
  Mm->matstiff (d,ipp);
  Mm->givestrain (0,ipp,eps);
  mxv (d.a, eps.a, spsig.a, d.m, d.n);
}
////////////////////       /* termitovo */       ////////////////////////////////////
