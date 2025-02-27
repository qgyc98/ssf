#include "beam2dspec.h"
#include "global.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "normmat.h"
#include "mechtop.h"
#include "probdesc.h"
#include "mechcrsec.h"
#include <math.h>


beam2dspec::beam2dspec (void)
{
  long i;
  
  nne=2;  ndofe=6;  tncomp=6;  napfun=3;
  nb=1;  intordmm=0;  intordism=0;  ssst=plbeam;

  ncomp = new long [nb];
  ncomp[0]=6;
  
  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
  }
  nip[0][0]=2;
  tnip=2;

}

beam2dspec::~beam2dspec (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
  }
  delete [] nip;
}



/**
   function assembles transformation matrix x_g = T x_l
   
   @param nodes - array containing node numbers
   @param tmat - transformation matrix
   
   JK, 3.2.2002
*/
void beam2dspec::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,n,m;

  fillm (0.0,tmat);

  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*3+0][i*3] = Mt->nodes[nodes[i]].e1[0];    tmat[i*3+0][i*3+1] = Mt->nodes[nodes[i]].e2[0];    tmat[i*3+0][i*3+2] = 0.0;
      tmat[i*3+1][i*3] = Mt->nodes[nodes[i]].e1[1];    tmat[i*3+1][i*3+1] = Mt->nodes[nodes[i]].e2[1];    tmat[i*3+1][i*3+2] = 0.0;
      tmat[i*3+2][i*3] = 0.0;                          tmat[i*3+2][i*3+1] = 0.0;                          tmat[i*3+2][i*3+2] = 1.0;
      
    }
  }
}


/**
   function assembles transformation matrix x_g = T x_l
   
   @param c,s - cosine and sine of beam angle
   @param tmat - transformation matrix
   
   JK, 3.2.2002
*/
void beam2dspec::beam_transf_matrix (double c,double s,matrix &tmat)
{
  fillm (0.0,tmat);
  
  tmat[0][0] = c;    tmat[0][1] = -1.0*s;    tmat[0][2] = 0.0;
  tmat[1][0] = s;    tmat[1][1] = c;         tmat[1][2] = 0.0;
  tmat[2][0] = 0.0;  tmat[2][1] = 0.0;       tmat[2][2] = 1.0;

  tmat[3][3] = c;    tmat[3][4] = -1.0*s;    tmat[3][5] = 0.0;
  tmat[4][3] = s;    tmat[4][4] = c;         tmat[4][5] = 0.0;
  tmat[5][3] = 0.0;  tmat[5][4] = 0.0;       tmat[5][5] = 1.0;
}




void beam2dspec::res_stiffness_matrix (long eid,matrix &sm)
{
  //stiffness_matrix (eid,0,0,sm);
  stiffness_matrix_expl (eid,0,0,sm);
}

/**
   function computes stiffness matrix of plane stress rectangular
   finite element with bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indices
   @param sm - stiffness matrix

   JK, 3.2.2002
*/
void beam2dspec::stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm)
{
  long i,ipp,transf;
  double l,s,c;
  ivector nodes(nne);
  vector x(nne),z(nne);
  matrix d(tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2dxz (x,z,eid);
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beam2dspec element",eid);
    fprintf (stderr,"\n in function beam2dspec::stiffness_matrix (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;
  
  fillm (0.0,sm);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  
  i = Mm->ip[ipp].idm[0];
  //  stiffness matrix in local coordinate system
  Mm->normm[i].stiffness_matrix (sm,ipp,l);
  
  
  
  //  transformation of stiffness matrix to the global system
  matrix tmat (ndofe,ndofe);
  beam_transf_matrix (c,s,tmat);
  lgmatrixtransf (sm,tmat);
  //  transformation of stiffness matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }

}

void beam2dspec::res_mass_matrix (long eid,matrix &mm)
{
  mass_matrix_expl (eid,0,0,mm);
}

/**
   function computes mass matrix of 2D beam element
   
   @param eid - number of element
   @param ri,ci - row and column indices
   @param mm - mass matrix
   
   JK, 3.2.200
*/
void beam2dspec::mass_matrix_expl (long eid,long ri,long ci,matrix &mm)
{
  long ii,transf;
  double e,g,shearcoeff,j,iy,a,l,c,s,kappa,rho;
  ivector nodes(nne);
  vector x(nne),z(nne);
  matrix d (tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2dxz (x,z,eid);
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beam2dspec element",eid);
    fprintf (stderr,"\n in function beam2dspec::mass_matrix (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;
  
  fillm (0.0,mm);
  ii=Mt->elements[eid].ipp[ri][ci];
  Mc->give_area (eid,a);
  Mc->give_shearcoeff (eid,&shearcoeff);
  Mc->give_densitye (eid,rho);
  Mc->give_mominer (eid,&iy);
  Mm->matstiff (d,ii);
  
  e=d[0][0];  g=d[1][1];
  if (shearcoeff<Mp->zero)  kappa=0.0;
  else                      kappa=6.0*e*iy/shearcoeff/g/a/l/l;

  j = (1.0+2.0*kappa)*(1.0+2.0*kappa);

  mm[0][0]=rho*a*l/3.0;
  mm[0][3]=rho*a*l/6.0;
  
  mm[1][1]=rho*a*l*(13.0/35.0+7.0/5.0*kappa+4.0/3.0*kappa*kappa)/j;
  mm[1][2]=rho*a*l*l*(-11.0/210.0-11.0/60*kappa-1.0/6.0*kappa*kappa)/j;
  mm[1][4]=rho*a*l*(9.0/70+3.0/5.0*kappa+2.0/3.0*kappa*kappa)/j;
  mm[1][5]=rho*a*l*l*(13.0/420+3.0/20.0*kappa+1.0/6.0*kappa*kappa)/j;
  
  mm[2][1]=mm[1][2];
  mm[2][2]=rho*a*l*l*l*(1.0/105.0+1.0/30.0*kappa+1.0/30.0*kappa*kappa)/j;
  mm[2][4]=rho*a*l*l*(-13.0/420.0-3.0/20.0*kappa-1.0/6.0*kappa*kappa)/j;
  mm[2][5]=rho*a*l*l*l*(-1.0/140.0-1.0/30.0*kappa-1.0/30.0*kappa*kappa)/j;
  
  mm[3][0]=mm[0][3];
  mm[3][3]=mm[0][0];
  
  mm[4][1]=mm[1][4];
  mm[4][2]=mm[2][4];
  mm[4][4]=rho*a*l*(13.0/35.0+7.0/5.0*kappa+4.0/3.0*kappa*kappa)/j;
  mm[4][5]=rho*a*l*l*(11.0/210.0+11.0/60.0*kappa+1.0/6.0*kappa*kappa)/j;
  
  mm[5][1]=mm[1][5];
  mm[5][2]=mm[2][5];
  mm[5][4]=mm[4][5];
  mm[5][5]=rho*a*l*l*l*(1.0/105+1.0/30.0*kappa+1.0/30.0*kappa*kappa)/j;

  //  transformation of mass matrix to the global system
  matrix tmat (ndofe,ndofe);
  beam_transf_matrix (c,s,tmat);
  lgmatrixtransf (mm,tmat);
  //  transformation of mass matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}



/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ifor - vector of internal forces
   
   12.8.2001
*/
void beam2dspec::internal_forces (long lcid,long eid,vector &ifor)
{
  long i,ipp;
  double l,c,s;
  vector x(nne),z(nne),lifor(ndofe),r(ndofe);
  matrix d(tncomp,tncomp);
  
  //  node coordinates
  Mt->give_node_coord2dxz (x,z,eid);
  
  //  node displacements
  eldispl (lcid,eid,r.a);

  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beam2dspec element",eid);
    fprintf (stderr,"\n in function beam2dspec::internal_forces (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;

  //  deleting array
  fillv (0.0,ifor);
  
  //  number of integration point
  ipp=Mt->elements[eid].ipp[0][0];
  //  number of material
  i=Mm->ip[ipp].idm[0];
  

  //  internal forces
  Mm->normm[i].internal_forces (ifor,ipp,l);
  
  /*
  //  transformation matrix from local to global coordinate system
  matrix tmat (ndofe,ndofe);
  beam_transf_matrix (c,s,tmat);
  
  //  transformation of local forces into global forces
  mxv (tmat,lifor,ifor);
  */
}

void beam2dspec::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,ifor);
}


/**
   function computes internal forces

   @param lcid - load case id
   @param eid - element id
   @param ifor - vector of internal forces
   
   12.8.2001
*/
void beam2dspec::stresses (long eid,long lcid)
{
  long i,ipp;
  vector x(nne),z(nne),f(ndofe),r(ndofe);
  matrix d(tncomp,tncomp);
  
  //  node coordinates
  Mt->give_node_coord2dxz (x,z,eid);
  
  //  node displacements
  eldispl (lcid,eid,r.a);

  //  number of integration point
  ipp=Mt->elements[eid].ipp[0][0];
  //  number of material
  i=Mm->ip[ipp].idm[0];
  
  matrix sm(6,6);
  stiffness_matrix_expl (eid,0,0,sm);
  
  mxv (sm,r,f);
  
  for (i=0;i<6;i++){
    Mm->ip[ipp].stress[i]=f[i];
  }
  
}

/**
   function computes strains
   
   @param eid - element id
   @param lcid - load case id
   
   20.2.2002
*/
void beam2dspec::strains (long eid,long lcid)
{
  long ii,transf;
  double l,s,c;
  ivector nodes(nne);
  vector x(nne),z(nne),rl(ndofe),rg(ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2dxz (x,z,eid);

  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beamel2d element",eid);
    fprintf (stderr,"\n in function beamel2d::strains (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;

  eldispl (lcid,eid,rl.a);
  
  //  transformation of node displacements
  matrix tmat (ndofe,ndofe);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    transf_matrix (nodes,tmat);
    //locglobtransf (rg,rl,tmat);
    lgvectortransf (rg,rl,tmat);
  }
  else{
    copyv (rl,rg);
  }
  beam_transf_matrix (c,s,tmat);
  //globloctransf (rg,rl,tmat);
  glvectortransf (rg,rl,tmat);

  
  ii=Mt->elements[eid].ipp[0][0];
  Mm->storestrain (lcid,ii,rl);
}


/**
   function stores end displacements and rotations to integration points
   
   JK, 24.5.2006
*/
void beam2dspec::res_mainip_strains (long lcid,long eid)
{
  long ipp;
  vector aux,r(ndofe),eps(3);
  ivector nodes(nne);
  matrix tmat;

  //  nodal displacements and rotations
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
  
  eps[0]=r[0];
  eps[1]=r[1];
  eps[2]=r[2];
  
  //  number of integration point
  ipp=Mt->elements[eid].ipp[0][0];
  //  storage of end displacements and rotations
  Mm->storestrain (lcid,ipp,0,eps);

}

/**
   function computes 
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 10.10.2005
*/
void beam2dspec::res_mainip_stresses (long lcid,long eid)
{
  double s,c,l;
  ivector nodes(nne);
  vector x(nne),z(nne),r(ndofe),gf(ndofe),lf(ndofe),aux;
  matrix sm(ndofe,ndofe),tmat(ndofe,ndofe);

  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  Mt->give_elemnodes (eid,nodes);
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    //reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }

  
  //  stiffness matrix
  stiffness_matrix_expl (eid,0,0,sm);
  //  internal forces in global system
  mxv (sm,r,gf);

  

  //  node coordinates in x-z plane
  Mt->give_node_coord2dxz (x,z,eid);
  
  //  length of the element
  l=sqrt ((x[1]-x[0])*(x[1]-x[0])+(z[1]-z[0])*(z[1]-z[0]));
  if (l<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beamel2d element",eid);
    fprintf (stderr,"\n in function beamel2d::res_mainip_stresses (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
  //  sine and cosine
  s=(z[1]-z[0])/l;
  c=(x[1]-x[0])/l;
  
  //  transformation matrix
  //  local system - system connected with the element
  //  global system - system defined for all problem
  beam_transf_matrix (c,s,tmat);
  
  //  transformation of internal forces from global to local system
  //globloctransf (gf,lf,tmat);
  lgvectortransf (gf,lf,tmat);
  
  long ipp;
  ipp=Mt->elements[eid].ipp[0][0];
  Mm->ip[ipp].stress[0]=lf[0]*(-1.0);
  Mm->ip[ipp].stress[1]=lf[1]*(-1.0);
  Mm->ip[ipp].stress[2]=lf[2]*(-1.0);
  //fprintf (Out,"\n stiffness %le %le %le",Mm->ip[ipp].stress[0],Mm->ip[ipp].stress[1],Mm->ip[ipp].stress[2]);
  ipp++;


  Mm->ip[ipp].stress[0]=lf[3];
  Mm->ip[ipp].stress[1]=lf[4];
  Mm->ip[ipp].stress[2]=lf[5];
  //fprintf (Out,"     %le %le %le",Mm->ip[ipp].stress[0],Mm->ip[ipp].stress[1],Mm->ip[ipp].stress[2]);
  
}
