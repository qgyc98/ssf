#include "dkt.h"
#include "gadaptivity.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>

dktelem::dktelem (void)
{
  long i,j;

  //  number of nodes on element
  nne=3;
  //  number of DOFs on element
  ndofe=9;
  //  number of strain/stress components
  tncomp=3;
  //  number of functions approximated
  napfun=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
  // number of surfaces
  nsurf=1;
  // number of nodes on surface
  nnsurf = nne;

  intordmm=3;
  //  strain/stress state
  ssst=platek;

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

  nip[0][0]=3;

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=3;
}

dktelem::~dktelem (void)
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
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

   JK, 23.9.2008
*/
double dktelem::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function assembles %matrix of base functions
   
   @param gm - geometric %matrix
   @param x,y - array containing node coordinates
   @param l - areacoordinates
   
   15.3.2002
*/
void dktelem::geom_matrix (matrix &gm,vector &x,vector &y,vector &l)
{
  double det,dx1,dy1,dx12,dy12,dx31,dy31;
  vector sx(ASTCKVEC(3)),sy(ASTCKVEC(3)),q(ASTCKVEC(3)),pb(ASTCKVEC(3)),pc(ASTCKVEC(3)),dl(ASTCKVEC(3));
  
  //  det is equal to double area of the element
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
  dl[0] = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1] = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
  //  der   dL(1)/dy . det
  sx[0]=x[2]-x[1];  sx[1]=x[0]-x[2];  sx[2]=x[1]-x[0];
  //  der   dL(1)/dx . det
  sy[0]=y[1]-y[2];  sy[1]=y[2]-y[0];  sy[2]=y[0]-y[1];
  //  dl0 = between nodes 1,2
  q[0]=dl[0]*dl[0];
  q[1]=dl[1]*dl[1];
  q[2]=dl[2]*dl[2];
  pb[0]=(sy[0]*sy[0]/4.-sx[0]*sx[0]/2.)/q[1];
  pb[1]=(sy[1]*sy[1]/4.-sx[1]*sx[1]/2.)/q[2];
  pb[2]=(sy[2]*sy[2]/4.-sx[2]*sx[2]/2.)/q[0];
  pc[0]=(sx[0]*sx[0]/4.-sy[0]*sy[0]/2.)/q[1];
  pc[1]=(sx[1]*sx[1]/4.-sy[1]*sy[1]/2.)/q[2];
  pc[2]=(sx[2]*sx[2]/4.-sy[2]*sy[2]/2.)/q[0];
  
  dx1= sy[0]*(l[0]-1./4.)/det;
  dy1= sx[0]*(l[0]-1./4.)/det;
  dx12=(sy[0]*l[1]+sy[1]*l[0])/det;
  dy12=(sx[0]*l[1]+sx[1]*l[0])/det;
  dx31=(sy[2]*l[0]+sy[0]*l[2])/det;
  dy31=(sx[2]*l[0]+sx[0]*l[2])/det;       
  gm[0][0]= 6.*(sx[2]/q[0]*dx12-sx[1]/q[2]*dx31);
  gm[0][1]=-3.*(sx[1]*sy[1]/q[2]*dx31+sx[2]*sy[2]/q[0]*dx12);
  gm[0][2]= 4.*(dx1-pc[1]*dx31-pc[2]*dx12);
  gm[1][0]=-6.*(sy[2]/q[0]*dy12-sy[1]/q[2]*dy31);
  gm[1][1]=-4.*(dy1-pb[1]*dy31-pb[2]*dy12);
  gm[1][2]= 3.*(sx[1]*sy[1]/q[2]*dy31+sx[2]*sy[2]/q[0]*dy12);
  gm[2][0]= 6.*(sx[2]/q[0]*dy12-sx[1]/q[2]*dy31 -sy[2]/q[0]*dx12+sy[1]/q[2]*dx31);
  gm[2][1]=-3.*(sx[1]*sy[1]/q[2]*dy31+sx[2]*sy[2]/q[0]*dy12) -4.*(dx1-pb[1]*dx31-pb[2]*dx12);
  gm[2][2]= 4.*(dy1-pc[1]*dy31-pc[2]*dy12) +3.*(sx[1]*sy[1]/q[2]*dx31+sx[2]*sy[2]/q[0]*dx12);
  
  dx1= sy[1]*(l[1]-1./4.)/det;
  dy1= sx[1]*(l[1]-1./4.)/det;
  dx12=(sy[1]*l[2]+sy[2]*l[1])/det;
  dy12=(sx[1]*l[2]+sx[2]*l[1])/det;
  dx31=(sy[0]*l[1]+sy[1]*l[0])/det;
  dy31=(sx[0]*l[1]+sx[1]*l[0])/det;       
  gm[0][3]= 6.*(sx[0]/q[1]*dx12-sx[2]/q[0]*dx31);
  gm[0][4]=-3.*(sx[2]*sy[2]/q[0]*dx31+sx[0]*sy[0]/q[1]*dx12);
  gm[0][5]= 4.*(dx1-pc[2]*dx31-pc[0]*dx12);
  gm[1][3]=-6.*(sy[0]/q[1]*dy12-sy[2]/q[0]*dy31);
  gm[1][4]=-4.*(dy1-pb[2]*dy31-pb[0]*dy12);
  gm[1][5]= 3.*(sx[2]*sy[2]/q[0]*dy31+sx[0]*sy[0]/q[1]*dy12);
  gm[2][3]= 6.*(sx[0]/q[1]*dy12-sx[2]/q[0]*dy31 -sy[0]/q[1]*dx12+sy[2]/q[0]*dx31);
  gm[2][4]=-3.*(sx[2]*sy[2]/q[0]*dy31+sx[0]*sy[0]/q[1]*dy12) -4.*(dx1-pb[2]*dx31-pb[0]*dx12);
  gm[2][5]= 4.*(dy1-pc[2]*dy31-pc[0]*dy12) +3.*(sx[2]*sy[2]/q[0]*dx31+sx[0]*sy[0]/q[1]*dx12);
  
  dx1= sy[2]*(l[2]-1./4.)/det;
  dy1= sx[2]*(l[2]-1./4.)/det;
  dx12=(sy[2]*l[0]+sy[0]*l[2])/det;
  dy12=(sx[2]*l[0]+sx[0]*l[2])/det;
  dx31=(sy[1]*l[2]+sy[2]*l[1])/det;
  dy31=(sx[1]*l[2]+sx[2]*l[1])/det;       
  gm[0][6]= 6.*(sx[1]/q[2]*dx12-sx[0]/q[1]*dx31);
  gm[0][7]=-3.*(sx[0]*sy[0]/q[1]*dx31+sx[1]*sy[1]/q[2]*dx12);
  gm[0][8]= 4.*(dx1-pc[0]*dx31-pc[1]*dx12);
  gm[1][6]=-6.*(sy[1]/q[2]*dy12-sy[0]/q[1]*dy31);
  gm[1][7]=-4.*(dy1-pb[0]*dy31-pb[1]*dy12);
  gm[1][8]= 3.*(sx[0]*sy[0]/q[1]*dy31+sx[1]*sy[1]/q[2]*dy12);
  gm[2][6]= 6.*(sx[1]/q[2]*dy12-sx[0]/q[1]*dy31 -sy[1]/q[2]*dx12+sy[0]/q[1]*dx31);
  gm[2][7]=-3.*(sx[0]*sy[0]/q[1]*dy31+sx[1]*sy[1]/q[2]*dy12) -4.*(dx1-pb[0]*dx31-pb[1]*dx12);
  gm[2][8]= 4.*(dy1-pc[0]*dy31-pc[1]*dy12) +3.*(sx[0]*sy[0]/q[1]*dx31+sx[1]*sy[1]/q[2]*dx12);
  
}

void dktelem::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i,n,m;
  n=nodes.n;
  m=tmat.m;
  for (i=0;i<m;i++){
    tmat[i][i]=1.0;
  }
  
  for (i=0;i<n;i++){
    if (Mt->nodes[nodes[i]].transf>0){
      tmat[i*3+1][i*3+1]=Mt->nodes[nodes[i]].e1[0];  tmat[i*3+1][i*3+2]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*3+2][i*3+1]=Mt->nodes[nodes[i]].e1[1];  tmat[i*3+2][i*3+2]=Mt->nodes[nodes[i]].e2[1];
    }
  }
}

/**
   function extracts components of the shear
   part of stiffness %matrix of the material to the %matrix ds

   20.3.2002
*/
void dktelem::dbmat (matrix &d,matrix &db,double t)
{
  double c;

  c=t*t*t;

  db[0][0] = c*d[0][0];  db[0][1] = c*d[0][1];  db[0][2] = c*d[0][2];
  db[1][0] = c*d[1][0];  db[1][1] = c*d[1][1];  db[1][2] = c*d[1][2];
  db[2][0] = c*d[2][0];  db[2][1] = c*d[2][1];  db[2][2] = c*d[2][2];
}

/**
   function computes stiffness %matrix of dkt element

   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix

   15.3.2002
*/
void dktelem::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii;
  double jac,det,thick;
  ivector nodes(nne);
  vector gp1(ASTCKVEC(intordsm[0][0])), gp2(ASTCKVEC(intordsm[0][0]));
  vector w(ASTCKVEC(intordsm[0][0])), l(ASTCKVEC(3)), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(3,ndofe)), d(ASTCKMAT(tncomp,tncomp)), db(ASTCKMAT(tncomp,tncomp));

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  fillm (0.0,sm);

  ii=Mt->elements[eid].ipp[ri][ci];

  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);

  for (i=0;i<intordsm[0][0];i++){
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];

    geom_matrix (gm,x,y,l);

    Mm->matstiff (d,ii);  ii++;
    thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];
    dbmat (d,db,thick);
    jac=det*w[i];
    // bdbj (sm.a,gm.a,db.a,jac,gm.m,gm.n);
    bdbjac (sm, gm, db, gm, jac);
  }
}

void dktelem::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  ivector nodes(ASTCKIVEC(nne));
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);

  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}


void dktelem::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 0.0;  eta[0] = 0.0;
  xi[1] = 1.0;  eta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 1.0;
}


void dktelem::appval (vector &l, long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval;
  
  k=0;
  reallocv (nne,nodval);
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k] = nodval[0]*l[0]+nodval[1]*l[1]+nodval[2]*l[2];
    k++;
  }
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
   
   JK, 23.9.2008
*/
void dktelem::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  vector gp1,gp2,w,eps,l(ASTCKVEC(3));
  matrix gm;

  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],gp1));
    reallocv (RSTCKVEC(intordsm[ii][ii],gp2));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(ncomp[ii],eps));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));

    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];

    for (i=0;i<intordsm[ii][ii];i++){
      l[0]=gp1[i];
      l[1]=gp2[i];
      l[2]=1.0-l[0]-l[1];
      geom_matrix (gm,x,y,l);
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,cncomp[ii],eps);
      ipp++;
    }
  }
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
   
   JK, 23.9.2008
*/
void dktelem::res_ip_strains (long lcid,long eid)
{
  vector aux,x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  ivector nodes(ASTCKIVEC(nne));
  matrix tmat;

  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
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
  
  ip_strains (lcid,eid,0,0,x,y,r);

}

void dktelem::nod_strains_ip (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
}


/**
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void dktelem::res_ip_stresses (long lcid,long eid)
{
  compute_nlstress (lcid,eid,0,0);
}

void dktelem::nod_stresses_ip (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
}




void dktelem::strains (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
}



void dktelem::stresses (long /*lcid*/,long /*eid*/,long /*ri*/, long /*ci*/)
{
}





/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.9.2008
*/
void dktelem::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double thick;
  ivector nodes(ASTCKIVEC(nne));
  vector sig(ASTCKVEC(3)), gp1(ASTCKVEC(intordsm[0][0]));
  vector gp2(ASTCKVEC(intordsm[0][0])), w(ASTCKVEC(intordsm[0][0]));
  vector l(ASTCKVEC(3)), t(ASTCKVEC(nne));
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
	//  computation of correct stresses
	if (Mp->strcomp==1){
	  l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
	  thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];
	  
	  Mm->computenlstresses (ipp,Mm->ip[ipp]);
	  Mm->givestress (lcid,ipp,sig);
	  sig[0]*=thick*thick*thick;
	  sig[1]*=thick*thick*thick;
	  sig[2]*=thick*thick*thick;
	  Mm->storestress (lcid,ipp,sig);
 	  ipp++;
	}
      }
    }
  }
}


/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 23.9.2008
*/
void dktelem::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct increments of stresses
	  if (Mp->strcomp==1)
            Mm->computenlstressesincr (ipp);
	  ipp++;
	}
      }
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
void dktelem::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}

/**
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void dktelem::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
  Mt->give_node_coord2d (x,y,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat(ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
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
void dktelem::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
   
   JK, 23.9.2008
*/
void dktelem::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (ASTCKIVEC(nne));
  vector v(ASTCKVEC(ndofe)), x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  
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
   function integrates selected quantity over the finite element
   it results in nodal values
   
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param nv - nodal values
   @param x,y - node coordinates
   
   JK, 23.9.2008
*/
void dktelem::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double jac,thick,det;
  ivector nodes(ASTCKIVEC(nne));
  vector w,gp1,gp2,t(ASTCKVEC(nne)),ipv,contr(ASTCKVEC(ndofe)),l(ASTCKVEC(3));
  matrix gm;
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  fillv (0.0,nv);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  for (ii=0;ii<nb;ii++){
    reallocv (RSTCKVEC(intordsm[ii][ii],gp1));
    reallocv (RSTCKVEC(intordsm[ii][ii],gp2));
    reallocv (RSTCKVEC(intordsm[ii][ii],w));
    reallocv (RSTCKVEC(ncomp[ii],ipv));
    reallocm (RSTCKMAT(ncomp[ii],ndofe,gm));
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){

      l[0]=gp1[i];
      l[1]=gp2[i];
      l[2]=1.0-l[0]-l[1];

      thick = approx (l,t);
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
	
      //  strain-displacement (geometric) matrix
      geom_matrix (gm,x,y,l);
      
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      jac=det*w[i];
      cmulv (jac,contr);
      
      //  summation
      addv(contr,nv,nv);
      
      ipp++;
    }
  }
}

/**
  function computes load vector from edge forces fz of the DKT
  15.3.2002
*/
void dktelem::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  double cn0,sn0,cn1,sn1,cn2,sn2;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),dl(ASTCKVEC(3));
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);

  if (le[0]==1){
    cn0= (y[1]-y[0])*dl[0]/60.;
    sn0=-(x[1]-x[0])*dl[0]/60.;
    nf[1]+=-(3.*nv[0]+2.*nv[3])*cn0;
    nf[2]+=-(3.*nv[0]+2.*nv[3])*sn0;
    nf[4]+= (2.*nv[0]+3.*nv[3])*cn0;
    nf[5]+= (2.*nv[0]+3.*nv[3])*sn0;
    nf[0]+=((2.-0.1)*nv[0] +(1.+0.1)*nv[3])*dl[0]/6.;
    nf[3]+=((1.-0.1)*nv[0] +(2.+0.1)*nv[3])*dl[0]/6.;
  }
  if (le[1]==1){
    cn1= (y[2]-y[1])*dl[1]/60.;
    sn1=-(x[2]-x[1])*dl[1]/60.;
    nf[4]+=-(3.*nv[1]+2.*nv[0])*cn1;
    nf[5]+=-(3.*nv[1]+2.*nv[0])*sn1;
    nf[7]+= (2.*nv[1]+3.*nv[0])*cn1;
    nf[8]+= (2.*nv[1]+3.*nv[0])*sn1;
    nf[3]+=((2.-0.1)*nv[3] +(1.+0.1)*nv[6])*dl[1]/6.;
    nf[6]+=((1.-0.1)*nv[3] +(2.+0.1)*nv[6])*dl[1]/6.;
  }
  if (le[2]==1){
    cn2= (y[0]-y[2])*dl[2]/60.;
    sn2=-(x[0]-x[2])*dl[2]/60.;
    nf[7]+=-(3.*nv[2]+2.*nv[1])*cn2;
    nf[8]+=-(3.*nv[2]+2.*nv[1])*sn2;
    nf[1]+= (2.*nv[2]+3.*nv[1])*cn2;
    nf[2]+= (2.*nv[2]+3.*nv[1])*sn2;
    nf[6]+=((2.-0.1)*nv[6] +(1.+0.1)*nv[9])*dl[2]/6.;
    nf[0]+=((1.-0.1)*nv[6] +(2.+0.1)*nv[9])*dl[2]/6.;
  }
}


/**
  function computes load vector from area forces fz of the DKT
  15.3.2002
*/
void dktelem::areaforces (long eid,double *nv,vector &nf)
{
   double pl,cn0,cn1,cn2,sn0,sn1,sn2;
   ivector nodes(ASTCKIVEC(nne));
   vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));

   Mt->give_elemnodes (eid,nodes);
   Mt->give_node_coord2d (x,y,eid);
  
   pl = fabs(((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.);

   cn0= (y[1]-y[0])*pl/360.;
   cn1= (y[2]-y[1])*pl/360.;
   cn2= (y[0]-y[2])*pl/360.;
   sn0=-(x[1]-x[0])*pl/360.;
   sn1=-(x[2]-x[1])*pl/360.;
   sn2=-(x[0]-x[2])*pl/360.;
   nf[1]=-(7.*nv[0]+3.*nv[3]+5.*nv[6])*cn2+(7.*nv[0]+5.*nv[3]+3.*nv[6])*cn0;
   nf[2]=-(7.*nv[0]+3.*nv[3]+5.*nv[6])*sn2+(7.*nv[0]+5.*nv[3]+3.*nv[6])*sn0;
   nf[4]=-(5.*nv[0]+7.*nv[3]+3.*nv[6])*cn0+(3.*nv[0]+7.*nv[3]+5.*nv[6])*cn1;
   nf[5]=-(5.*nv[0]+7.*nv[3]+3.*nv[6])*sn0+(3.*nv[0]+7.*nv[3]+5.*nv[6])*sn1;
   nf[7]=-(3.*nv[0]+5.*nv[3]+7.*nv[6])*cn1+(5.*nv[0]+3.*nv[3]+7.*nv[6])*cn2;
   nf[8]=-(3.*nv[0]+5.*nv[3]+7.*nv[6])*sn1+(5.*nv[0]+3.*nv[3]+7.*nv[6])*sn2;
   nf[0]=(nv[0]/6. +nv[3]/12.+nv[6]/12.)*pl;
   nf[3]=(nv[0]/12.+nv[3]/6. +nv[6]/12.)*pl;
   nf[6]=(nv[0]/12.+nv[3]/12.+nv[6]/6. )*pl;


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
void dktelem::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double ipval;
  vector w, gp1, gp2, anv(ASTCKVEC(nne)), l(ASTCKVEC(3));
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
        reallocv (RSTCKVEC(intordsm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordsm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordsm[ii][jj],w));
        gauss_points_tr (gp1.a, gp2.a, w.a, intordsm[ii][jj]);
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          l[0]=gp1[k];
          l[1]=gp2[k];
          l[2]=1.0-l[0]-l[1];
          //  value in integration point
          ipval = approx(l, anv);
          if ((ictn[i] & inistrain) && (j < Mm->ip[ipp].ncompstr))
          {
            Mm->ip[ipp].strain[j] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[i] & inistress) && (j < nstra + Mm->ip[ipp].ncompstr))
          {
            Mm->ip[ipp].stress[j] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[i] & iniother) && (j < nv))
          {
            Mm->ip[ipp].other[j] += ipval;
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

