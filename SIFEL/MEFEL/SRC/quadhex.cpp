#include "quadhex.h"
#include "linhex.h"
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
#include "loadcase.h"
#include "intpoints.h"
#include <stdlib.h>


quadhex::quadhex (void)
{
  long i,j;

  //  number of nodes on element
  nne=20;
  //  number of DOFs on element
  ndofe=60;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=2;
  //  number of edges on element
  ned=12;
  //  number of nodes on one edge
  nned=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=3;
  //  number of surfaces on element
  nsurf=6;
  //  number of nodes on one surface
  nnsurf=8;
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
  
  nip[0][0]=27;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=3;

}

quadhex::~quadhex (void)
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
   
   @param xi,eta,zeta - natural coordinates
   @param nodval - nodal values
   
   JK, 20.8.2001
*/
double quadhex::approx (double xi,double eta,double zeta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_hex_3d (bf.a,xi,eta,zeta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function assembles %matrix of base functions
   
   @param n - %matrix of base functions
   @param xi,eta,zeta - coordinates

   JK, 16.8.2001
*/
void quadhex::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  long i,j,k,l;
  vector bf(nne);
  
  nullm (n);

  bf_quad_hex_3d (bf.a,xi,eta,zeta);
  
  j=0;  k=1;  l=2;
  for (i=0;i<nne;i++){
    n[0][j]=bf[i];  j+=3;
    n[1][k]=bf[i];  k+=3;
    n[2][l]=bf[i];  l+=3;
  }
}

/**
   function computes strain-displacement (geometric) %matrix

   @param gm - geometric %matrix
   @param x,y,z - vectors containing element node coordinates
   @param xi,eta,zeta - coordinates
   @param jac - Jacobian

   JK, 16.8.2001
*/
void quadhex::geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
			   double xi,double eta,double zeta,double &jac)
{
  long i,j,k,l;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_quad_hex_3d (dx.a,xi,eta,zeta);
  dy_bf_quad_hex_3d (dy.a,xi,eta,zeta);
  dz_bf_quad_hex_3d (dz.a,xi,eta,zeta);

  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);

  nullm (gm);

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
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   @param nodes - nodes of element
   @param tmat - transformation %matrix
   
   JK
*/
void quadhex::transf_matrix (ivector &nodes,matrix &tmat)
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

   JK, 16.8.2001
*/
void quadhex::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  Mt->give_node_coord3d (x,y,z,eid);
  nullm (sm);

  reallocv (intordsm[0][0],w);
  reallocv (intordsm[0][0],gp);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
	
	//  geometric matrix
	geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	
	//  stiffness matrix of the material
	Mm->matstiff (d,ipp);  ipp++;
	
	jac*=w[i]*w[j]*w[k];
	
	bdbjac (sm,gm,d,gm,jac);
      }
    }
  }
}

/**
   function computes stiffness %matrix of one element

   @param eid - number of element
   @param sm - stiffness %matrix

   JK, 16.8.2001
*/
void quadhex::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes (nne);
  
  stiffness_matrix (eid,0,0,sm);
  
  //  transformation of stiffness matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes mass %matrix
   
   @param eid - number of element
   @param mm - mass %matrix

   JK, 16.8.2001
*/
void quadhex::mass_matrix (long eid,matrix &mm)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3,rho;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm),dens(nne);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_density (eid,nodes,dens);

  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (mm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];  w3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);

	bf_matrix (n,xi,eta,zeta);

	rho = approx (xi,eta,zeta,dens);
	
	jac*=w1*w2*w3*rho;
	
	nnj (mm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}

/**
   function computes mass %matrix
   
   @param eid - number of element
   @param mm - mass %matrix

   JK, 16.8.2001
*/
void quadhex::res_mass_matrix (long eid,matrix &mm)
{
  long transf;
  ivector nodes(nne);
  
  mass_matrix (eid,mm);

  if (Mp->diagmass==1){
    diagonalization (mm);
  }
  
  //  transformation of mass matrix
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (mm,tmat);
  }
}

/**
   function computes load %matrix
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 16.8.2001
*/
void quadhex::load_matrix (long eid,matrix &lm)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3;
  ivector nodes (nne);
  vector x(nne),y(nne),z(nne),w(intordmm),gp(intordmm);
  matrix n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordmm);
  
  nullm (lm);
  
  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordmm;k++){
	zeta=gp[k];  w3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);

	bf_matrix (n,xi,eta,zeta);

	jac*=w1*w2*w3;
	
	nnj (lm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  
}

/**
   function computes load %matrix
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 16.8.2001
*/
void quadhex::res_load_matrix (long eid,matrix &lm)
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
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, modified 23.11.2006
*/
void quadhex::res_ip_strains (long lcid,long eid)
{
  vector x(nne),y(nne),z(nne),r(ndofe),gp,w,eps,aux;
  ivector nodes(nne);
  matrix gm,tmat;

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
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

  ip_strains (lcid,eid,0,0,x,y,z,r);
}

/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y,z - %vectors of nodal coordinates
   @param r - %vector of nodal displacements
   
   JK, 10.5.2002, modified 27.11.2006
*/
void quadhex::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long i,j,k,ii,ipp;
  double xi,eta,zeta,jac;
  vector gp,w,eps,aux;
  ivector nodes(nne),cn(ndofe);
  matrix gm,tmat;

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);

    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  geom_matrix (gm,x,y,z,xi,eta,zeta,jac);

	  mxv (gm,r,eps);
	  
	  Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
	  ipp++;
	}
      }
    }
  }
}

/**
   function computes strains in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.9.2005
*/
void quadhex::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadhex (ipp,intordsm[0][0],ipnum);
  
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
   function computes strains in nodes of element
   
   @param lcid - load case id
   @param eid - element id
   @param stra - array for strain components
   
   stra[i][j] - the j-th strain component at the i-th node

   JK, 10.5.2002
*/
void quadhex::nod_strains_comp (long lcid,long eid,double **stra)
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
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transf_matrix (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  //  natural coordinates of element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodcoord_quadhex (nxi,neta,nzeta);

  //  loop over nodes
  for (i=0;i<nne;i++){
    //  geometric matrix
    geom_matrix (gm,x,y,z,nxi[i],neta[i],nzeta[i],jac);
    //  strain computation
    mxv (gm,r,eps);
    
    for (j=0;j<eps.n;j++){
      stra[i][j]=eps[j];
    }
  }

}




void quadhex::strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
  vector coord,eps;
  
  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_strains (lcid,eid);
    break;
  }
  case enodes:{
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stra.give_naep (eid);
    ncp = Mm->stra.give_ncomp (eid);
    sid = Mm->stra.give_sid (eid);
    reallocv (ncp,eps);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	//appstrain (lcid,eid,coord[0],coord[1],coord[2],0,ncp,eps);
      
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function planeelemlq::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}




/**
   function computes stresses at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 10.5.2002
*/
void quadhex::res_ip_stresses (long lcid,long eid)
{
  ip_stresses (lcid,eid,0,0);
}

/**
   function computes stresses at integration points of element
   stresses are computed by material models
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   10.5.2002, JK
*/
void quadhex::ip_stresses (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	
	ipp++;
      }
    }
  }
}

/**
   function computes stresses at integration points of element
   stresses are computed from strains with the help of elastic stiffness
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   
   JK, 10.5.2002
*/
void quadhex::ip_elast_stresses (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  vector gp,w,eps,sig(tncomp);
  matrix gm(tncomp,ndofe),d(tncomp,tncomp);
  
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  stiffness matrix of the material
	Mm->matstiff (d,ipp);
	
	//  strains
	Mm->givestrain (lcid,ipp,eps);
	
	//  elastic stresses
	mxv (d,eps,sig);
	
	Mm->storestress (lcid,ipp,sig);
	
	ipp++;
      }
    }
  }
}

/**
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.9.2005
*/
void quadhex::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadhex (ipp,intordsm[0][0],ipnum);
  
  //  node numbers of the element
  Mt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    //  stresses at the closest integration point
    Mm->givestress (lcid,ipnum[i],sig);
    
    //  storage of stresses to the node
    j=nod[i];
    Mt->nodes[j].storestress (lcid,0,sig);
  }
  
}

/**
   function computes stresses in nodes
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param stra -
   @param stre -

   10.5.2002
*/
void quadhex::nod_stresses_comp (long /*lcid*/,long eid,long ri,long ci,double **stra,double **stre)
{
  long i,j,ipp;
  vector eps(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  //  number of the first integration point on the element
  ipp=Mt->elements[eid].ipp[ri][ci];
  //  stiffness matrix of the material
  Mm->matstiff (d,ipp);
  
  //  loop over nodes
  for (i=0;i<nne;i++){
    for (j=0;j<eps.n;j++){
      eps[j]=stra[i][j];
    }
    mxv (d,eps,sig);
    for (j=0;j<eps.n;j++){
      stre[i][j]=sig[j];
    }
  }
}







void quadhex::stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
  vector coord,sig;
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    //allip_stresses (lcid,eid,ri,ci);
    break;
  }
  case enodes:{
    break;
  }
  case userdefined:{
    //  number of auxiliary element points
    naep = Mm->stre.give_naep (eid);
    ncp = Mm->stre.give_ncomp (eid);
    sid = Mm->stre.give_sid (eid);
    reallocv (ncp,sig);
    reallocv (3,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);
      
      if (Mp->stressaver==0)
	//appval (coord[0],coord[1],coord[2],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	//appstress (lcid,eid,coord[0],coord[1],coord[2],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}



/**
   function computes other values in nodes of element

   @param eid[in] - element id
   @param ri,ci[in] - row and column indices
   
   JK, 24.10.2005
*/
void quadhex::nod_other_ip (long eid, long ri, long ci)
{
  long i, j, ipp, ncompo;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_quadhex(ipp, intordsm[0][0], ipnum);

  
  //  node numbers of the element
  Mt->give_elemnodes(eid, nod);
  
  for (i=0;i<nne;i++){
    //  strains at the closest integration point
    //Mm->givestrain (ipnum[i],eps);
    
    ncompo = Mm->givencompother(ipnum[i],0);
    reallocv(RSTCKVEC(ncompo, other));
    Mm->giveother(ipnum[i], 0, ncompo, other.a);
    
    //  storage of strains to the node
    j=nod[i];
    Mt->nodes[j].storeother(0, ncompo, other);
  }
}



/**
   function computes internal forces in the case of geometrical linear computation
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates
   
   JK, 28.7.2001
   TKo 7.2008
*/
void quadhex::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x, vector &y, vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes nonlocal internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - %vector of internal forces
   @param x,y,z - nodal coordinates
   
   JK, 28.7.2001
   TKo 7.2008
*/
void quadhex::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of nonlocal stresses
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
void quadhex::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes contributions from eigenstrains to the right hand side
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates

   JK, 17.8.2004
   TKo 7.2008
*/
void quadhex::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
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
   function computes resulting internal forces (from correct stresses)
   
   @param lcid - number of load case
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 24.9.2005
*/
void quadhex::res_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes resulting internal forces for nonlocal models
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void quadhex::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
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



void quadhex::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes contributions from eigenstrains to the right hand side
   
   @param eid - element id
   @param ifor - %vector of internal forces

   JK, 17.8.2004
*/
void quadhex::res_eigstrain_forces (long lcid,long eid,vector &nfor)
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
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.11.2006
*/
void quadhex::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	ipp++;
      }
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
void quadhex::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstressesincr (ipp);
	
	ipp++;
      }
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
void quadhex::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	ipp++;
      }
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
void quadhex::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
	//  computation of correct stresses
	if (Mp->strcomp==1)
	  Mm->compnonloc_nlstresses (ipp);
	ipp++;
      }
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
void quadhex::compute_eigstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<intordsm[0][0];j++){
      for (k=0;k<intordsm[0][0];k++){
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
}



/**
   function integrates selected quantity over the finite element
   it results in nodal values
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y,z - node coordinates
   
   JK, 27.11.2006
   TKo 7.2008
*/
void quadhex::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x, vector&y, vector &z)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector w,gp,ipv(tncomp),contr(ndofe);
  matrix gm(tncomp,ndofe);
  
  nullv (nv);
  reallocv (intordsm[0][0],gp);
  reallocv (intordsm[0][0],w);
  gauss_points (gp.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[0][0];i++){
    xi=gp[i];
    for (j=0;j<intordsm[0][0];j++){
      eta=gp[j];
      for (k=0;k<intordsm[0][0];k++){
	zeta=gp[k];
        //  function assembles required quantity at integration point
        Mm->givequantity (iq,lcid,ipp,cncomp[0],ipv);
	//  strain-displacement (geometric) matrix
	geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	//  contribution to the internal forces
	mtxv (gm,ipv,contr);
	cmulv (jac*w[i]*w[j]*w[k],contr);
	//  summation
	addv (contr,nv,nv);
	ipp++;
      }
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
void quadhex::ipcoord (long eid,long ipp,long ri,long ci,vector &coord)
{
  long i,j,k,ii;
  double xi,eta,zeta;
  vector x(nne),y(nne),z(nne),w(intordsm[ri][ci]),gp(intordsm[ri][ci]);
  
  gauss_points (gp.a,w.a,intordsm[ri][ci]);
  Mt->give_node_coord3d (x,y,z,eid);
  ii=Mt->elements[eid].ipp[ri][ci];
  
  for (i=0;i<intordsm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordsm[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordsm[ri][ci];k++){
	zeta=gp[k];
	if (ii==ipp){
	  coord[0]=approx (xi,eta,zeta,x);
	  coord[1]=approx (xi,eta,zeta,y);
	  coord[2]=approx (xi,eta,zeta,z);
	}
	ii++;
      }
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
void quadhex::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, j, k, ii, ri, ci;
  double xi,eta,zeta;
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
          for (k=0;k<intordsm[ri][ci];k++){
            zeta=gp[k];
            if (ii==ipp){
              ncoord[0]=xi;
              ncoord[1]=eta;
              ncoord[2]=zeta;
            }
            ii++;
          }
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
void quadhex::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, m, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, zeta, ipval;
  vector w, gp, anv(nne);
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
        reallocv (intordsm[ii][jj],gp);
        reallocv (intordsm[ii][jj],w);
        gauss_points (gp.a,w.a,intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          xi=gp[k];
          for (l = 0; l < intordsm[ii][jj]; l++)
          {
            eta=gp[l];
            for (m = 0; m < intordsm[ii][jj]; m++)
            {
              zeta=gp[m];
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
void quadhex::ipvolume (long eid,long ri,long ci)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),w,gp;
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    jac_3d (jac,x,y,z,xi,eta,zeta);
	    jac=fabs(jac);
	    
	    jac*=w[i]*w[j]*w[k];
	    
	    Mm->storeipvol (ipp,jac);
	    ipp++;
	  }
	}
      }
    }
  }
}

/**
   function computes nodal forces caused by presure on surface
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - vector of presure 
   @param eis - surface id 
   
   4.2002, PF
*/
void quadhex::node_forces_surf (long /*lcid*/,long eid,long *is,double *nv,vector &nf)
{
  long i,j,ii;
  double xi=0.0,eta=0.0,zeta=0.0,jac, w1,w2;
  ivector nodes(nne);
  vector gx(nne),gy(nne),gz(nne),x(nnsurf),y(nnsurf),z(nnsurf),gp(intordb),w(intordb),av(ndofe),v(ndofe),pom(3);
  matrix n(napfun,ndofe),an(napfun,ndofe),am(ndofe,ndofe), tran(3,3);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (gx,gy,gz,eid);
  gauss_points (gp.a,w.a,intordb);
  nullv (nf);
  nullm (an);
  for (ii=0;ii<6;ii++){
    nullv (av);
    // is=0 not loading
    if (is[ii] ==0){}
    else {
      tran_mat( tran, gx, gy, gz, ii);
      // is=1 surface node 1,4,8,5,  12,16,20,13
      if (ii ==0){
	xi=1.0;
        x[0]=gx[0];x[1]=gx[3];x[2]=gx[7];x[3]=gx[4]; x[4]=gx[11];x[5]=gx[15];x[6]=gx[18];x[7]=gx[12];
        y[0]=gy[0];y[1]=gy[3];y[2]=gy[7];y[3]=gy[4]; y[4]=gy[11];y[5]=gy[15];y[6]=gy[18];y[7]=gy[12];
        z[0]=gz[0];z[1]=gz[3];z[2]=gz[7];z[3]=gz[4]; z[4]=gz[11];z[5]=gz[15];z[6]=gz[18];z[7]=gz[12];
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;     
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[i];  av[9+i]=nv[3+i];  av[21+i]=nv[6+i];  av[15+i]=nv[9+i];
	  av[33+i]=nv[12+i];  av[45+i]=nv[15+i];  av[57+i]=nv[18+i];  av[36+i]=nv[21+i];
	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=2 surface node 2,1,5,6,  9,13,17,14
      else if (ii ==1){
	eta=1.0;
        x[0]=gx[1];x[1]=gx[0];x[2]=gx[4];x[3]=gx[5]; x[4]=gx[8];x[5]=gx[12];x[6]=gx[16];x[7]=gx[13];
        y[0]=gy[1];y[1]=gy[0];y[2]=gy[4];y[3]=gy[5]; y[4]=gy[8];y[5]=gy[12];y[6]=gy[16];y[7]=gy[13];
        z[0]=gz[1];z[1]=gz[0];z[2]=gz[4];z[3]=gz[5]; z[4]=gz[8];z[5]=gz[12];z[6]=gz[16];z[7]=gz[13];
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[3+i]=nv[24+i];  av[i]=nv[27+i];  av[12+i]=nv[30+i];  av[15+i]=nv[33+i];
	  av[24+i]=nv[36+i];  av[36+i]=nv[39+i];  av[48+i]=nv[42+i];  av[39+i]=nv[45+i];
	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=3 surface node 3,2,6,7, 10, 14,18,15
      else if (ii ==2){
	xi=-1.0;
        x[0]=gx[2];x[1]=gx[1];x[2]=gx[5];x[3]=gx[6]; x[4]=gx[9];x[5]=gx[13];x[6]=gx[17];x[7]=gx[14];
        y[0]=gy[2];y[1]=gy[1];y[2]=gy[5];y[3]=gy[6]; y[4]=gy[9];y[5]=gy[13];y[6]=gy[17];y[7]=gy[14];
        z[0]=gz[2];z[1]=gz[1];z[2]=gz[5];z[3]=gz[6]; z[4]=gz[9];z[5]=gz[13];z[6]=gz[17];z[7]=gz[14];
	for (i=0;i<intordb;i++){
	  eta=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,eta,zeta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[6+i]=nv[48+i];  av[3+i]=nv[51+i];  av[15+i]=nv[54+i];  av[18+i]=nv[57+i];
	  av[27+i]=nv[60+i];  av[39+i]=nv[63+i];  av[51+i]=nv[66+i];  av[42+i]=nv[69+i];
	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=4 surface node 4,3,7,8,  11,15,19,16
      else if (ii ==3){
	eta=-1.0;
        x[0]=gx[3];x[1]=gx[2];x[2]=gx[6];x[3]=gx[7]; x[4]=gx[10];x[5]=gx[14];x[6]=gx[18];x[7]=gx[15];
        y[0]=gx[3];y[1]=gy[2];y[2]=gy[6];y[3]=gy[7]; y[4]=gy[10];y[5]=gy[14];y[6]=gy[18];y[7]=gy[15];
        z[0]=gx[3];z[1]=gz[2];z[2]=gz[6];z[3]=gz[7]; z[4]=gz[10];z[5]=gz[14];z[6]=gz[18];z[7]=gz[15];
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    zeta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,zeta);
	    jac = jac*w1*w2;
	    
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[9+i]=nv[72+i];  av[6+i]=nv[75+i];  av[18+i]=nv[78+i];  av[21+i]=nv[81+i];
	  av[30+i]=nv[84+i];  av[42+i]=nv[87+i];  av[54+i]=nv[90+i];  av[45+i]=nv[93+i];
	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=5  surface node 1,2,3,4,  9,10,11,12
      else if (ii ==4) {
	zeta=1.0;
        x[0]=gx[0];x[1]=gx[1];x[2]=gx[2];x[3]=gx[3]; x[4]=gx[8];x[5]=gx[9];x[6]=gx[10];x[7]=gx[11];
        y[0]=gy[0];y[1]=gy[1];y[2]=gy[2];y[3]=gy[3]; y[4]=gy[8];y[5]=gy[9];y[6]=gy[10];y[7]=gy[11];
        z[0]=gz[0];z[1]=gz[1];z[2]=gz[2];z[3]=gz[3]; z[4]=gz[8];z[5]=gz[9];z[6]=gz[10];z[7]=gz[11];
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    //	for constant			for (i1=0;i1<ndofe;i1++){
	    //								for (j1=0;j1<3;j1++){
	    //									an[j1][i1]=an[j1][i1]+n[j1][i1]*jac;}	}
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[i]=nv[96+i];  av[3+i]=nv[99+i];  av[6+i]=nv[102+i];  av[9+i]=nv[105+i];
	  av[24+i]=nv[108+i];  av[27+i]=nv[111+i];  av[30+i]=nv[114+i];  av[33+i]=nv[117+i];
	}
	//		fprintf (Out,"\n\n nn");
	//		for (i1=0;i1<ndofe;i1++){
	//			fprintf (Out,"\n %4ld   %20.10e %20.10e %20.10e",i1,an[0][i1],an[1][i1],an[2][i1]);}       
	
	// for constant		for (i=0;i<ndofe;i++){
	//						for (j=0;j<pom.n;j++){
	//							v[i]=an[j][i]*pom[j];}	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);
      }
      // is=6  surface node 5,6,7,8, 17,18,19,20
      else if (ii ==5) {
	zeta=-1.0;
        x[0]=gx[4];x[1]=gx[5];x[2]=gx[6];x[3]=gx[7]; x[4]=gx[16];x[5]=gx[17];x[6]=gx[18];x[7]=gx[19];
        y[0]=gy[4];y[1]=gy[5];y[2]=gy[6];y[3]=gy[7]; y[4]=gy[16];y[5]=gy[17];y[6]=gy[18];y[7]=gy[19];
        z[0]=gz[4];z[1]=gz[5];z[2]=gz[6];z[3]=gz[7]; z[4]=gz[16];z[5]=gz[17];z[6]=gz[18];z[7]=gz[19];
	for (i=0;i<intordb;i++){
	  xi=gp[i];  w1=w[i];      
	  for (j=0;j<intordb;j++){
	    eta=gp[j];  w2=w[j];
	    bf_matrix (n,xi,eta,zeta);
	    jac2d3d (jac,x,y,z,xi,eta);
	    jac = jac*w1*w2;      
	    nnj (am.a,n.a,jac,n.m,n.n);
	  }
	}
	for (i=0;i<3;i++){
	  av[12+i]=nv[120+i];  av[15+i]=nv[123+i];  av[18+i]=nv[126+i];  av[21+i]=nv[129+i];
	  av[48+i]=nv[132+i];  av[51+i]=nv[135+i];  av[54+i]=nv[138+i];  av[57+i]=nv[141+i];
	}
	// Load in GCS
	if (is[ii] ==1){
	  glvectortransfblock (av, tran);
	}
        mxv (am,av,v);  
      }
      
      lgvectortransfblock(v, tran);
      addv (nf,v,nf);      
    }
  }

  //  transformation of stiffness matrix to nodesystem
  //  transf = Mt->locsystems (nodes);
  //  if (transf>0){
  //    matrix tmat (ndofe,ndofe);
  //    transf_matrix (nodes,tmat);
  //    lgvectortransfblock (v,tmat);
  //  }

  //  fprintf (Out,"\n\n zatizeni na prvku cislo %ld",eid);
  //  for (i=0;i<ndofe;i++){
  //	  fprintf (Out,"\n %4ld   %20.10e",i,nf[i]);}
}

/**
   function computes transformation matrix on surface
   
   @param x, y - local coordinate  xL=adge12
   @param tran - tranformation metrix to GCS
   @param gx,gy,gz - vector of node global coordinate 
   @param is - surface id 
   
   4.2002, PF
*/
void quadhex::tran_mat(matrix &tran, vector &gx, vector &gy, vector &gz, long is)
{
  long i,i1;
  double dl;
  matrix a(3,3);
// is=5  surface node 1,2,3,4
  if (is ==4) {
    for (i=0; i<3; i++) {
      i1=i+1; if(i1>2)i1=i1-3;
      a[0][i]=gx[i1]-gx[i];
      a[1][i]=gy[i1]-gy[i];
      a[2][i]=gz[i1]-gz[i];
    }
  }
// is=6  surface node 5,6,7,8
  else if (is ==5) {
    for (i=4; i<7; i++) {
      i1=i+1; if(i1>6)i1=i1-3;
      a[0][i-4]=gx[i1]-gx[i];
      a[1][i-4]=gy[i1]-gy[i];
      a[2][i-4]=gz[i1]-gz[i];
    }
  }
// is=1 surface node 1,4,8,5
// is=2 surface node 2,1,5,6
// is=3 surface node 3,2,6,7
// is=4 surface node 4,3,7,8
  else {
      if(is>0)i1=is-4;else i1=is;
      a[0][0]=gx[i1+3]-gx[is];
      a[1][0]=gy[i1+3]-gy[is];
      a[2][0]=gz[i1+3]-gz[is];
      a[0][1]=gx[i1+7]-gx[i1+3];
      a[1][1]=gy[i1+7]-gy[i1+3];
      a[2][1]=gz[i1+7]-gz[i1+3];
      a[0][2]=gx[is]-gx[i1+7];
      a[1][2]=gy[is]-gy[i1+7];
      a[2][2]=gz[is]-gz[i1+7];
  }
  
  dl=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);

  tran[0][0]=a[0][0]/dl;
  tran[1][0]=a[1][0]/dl;
  tran[2][0]=a[2][0]/dl;
  
  tran[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
  tran[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
  tran[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

  dl=sqrt(tran[0][2]*tran[0][2]+tran[1][2]*tran[1][2]+tran[2][2]*tran[2][2]);

  tran[0][2]=tran[0][2]/dl;
  tran[1][2]=tran[1][2]/dl;
  tran[2][2]=tran[2][2]/dl;
  
  tran[0][1]=tran[1][2]*tran[2][0]-tran[2][2]*tran[1][0];
  tran[1][1]=tran[2][2]*tran[0][0]-tran[0][2]*tran[2][0];
  tran[2][1]=tran[0][2]*tran[1][0]-tran[1][2]*tran[0][0];

  dl=sqrt(tran[0][1]*tran[0][1]+tran[1][1]*tran[1][1]+tran[2][1]*tran[2][1]);

  tran[0][1]=tran[0][1]/dl;
  tran[1][1]=tran[1][1]/dl;
  tran[2][1]=tran[2][1]/dl;
/*  
// local coordinate x=s12
    sl1[0]=0.0;    sl1[1]=0.0;    sl1[2]=0.0;
    sl2[0]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
    sl2[1]=0.0;    sl2[2]=0.0;
    sl3[0]=(s3(0)-s1(0))*tran(0,0)+(s3(1)-s1(1))*tran(1,0)+(s3(2)-s1(2))*tran(2,0);
    sl3[1]=(s3(0)-s1(0))*tran(0,1)+(s3(1)-s1(1))*tran(1,1)+(s3(2)-s1(2))*tran(2,1);
    sl3[2]=0.0;
*/
// local coordinate x=s12
//    x[0]=0.0;
//    y[0]=0.0;
//    x[1]=sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]);
//    y[1]=0.0;
//    x[2]=(gx(2)-gx(0))*tran(0,0)+(gy(2)-gy(0))*tran(1,0)+(gz(2)-gz(0))*tran(2,0);
//    y[2]=(gx(2)-gx(0))*tran(0,1)+(gy(2)-gy(0))*tran(1,1)+(gz(2)-gz(0))*tran(2,1);
//    x[3]=(gx(3)-gx(0))*tran(0,0)+(gy(3)-gy(0))*tran(1,0)+(gz(3)-gz(0))*tran(2,0);
//    y[3]=(gx(3)-gx(0))*tran(0,1)+(gy(3)-gy(0))*tran(1,1)+(gz(3)-gz(0))*tran(2,1);
}








/**
   function interpolates the nodal values to the integration points on the element
   quadratic approximation functions are used
   
   @param eid - element id
   @param nodval - nodal values
   @param ipval - value at integration points
   
   17.8.2004, JK
*/
void quadhex::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k,l;
  double xi,eta,zeta;
  vector w,gp;
  
  l=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);

      gauss_points (gp.a,w.a,intordsm[ii][jj]);

      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    ipval[l]=approx (xi,eta,zeta,nodval);
	    l++;
	  }
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
   
   17.8.2004, JK
*/
void quadhex::intpointval2 (long /*eid*/,vector &nodval,vector &ipval)
{
  long i,j,ii,jj,k,l;
  double xi,eta,zeta;
  vector w,gp;
  vector modnodval(Lhex->nne);
  
  for (i=0;i<Lhex->nne;i++){
    modnodval[i]=nodval[i];
  }
  
  l=0;
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
	xi=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  eta=gp[j];
	  for (k=0;k<intordsm[ii][jj];k++){
	    zeta=gp[k];
	    
	    ipval[l]=Lhex->approx (xi,eta,zeta,modnodval);
	    l++;
	  }
	}
      }
    }
  }
}



void quadhex::aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  vector x(nne),y(nne),z(nne),gp,w,eps;

  Mt->give_node_coord3d (x,y,z,eid);
  
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    
    reallocv (intordsm[ii][ii],gp);
    reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    for (i=0;i<intordsm[ii][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  //Mm->givestrain (lcid,ipp,eps);
	  Mm->givestress (lcid,ipp,eps);
	  
	  
          jac_3d (jac,x,y,z,xi,eta,zeta);	
          jac=fabs(jac);

	  volume+=w[i]*w[j]*w[k]*jac;
	  
	  for (l=0;l<averstra.n;l++){
	    averstra[l]+=eps[l]*w[i]*w[j]*w[k]*jac;
	  }
	  
	  ipp++;
	}
      }
    }
  }
}
