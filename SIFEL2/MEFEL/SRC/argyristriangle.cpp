#include "argyristriangle.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mathem.h"
#include "intp.h"
#include "difcalc.h"
#include "global.h"
#include "gadaptivity.h"
#include "globmat.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>

// ==================================================================


ArgyrisTriangle::ArgyrisTriangle (void)
{
  long i,j;

  //  number nodes on element
  nne=3;
  //  number of DOFs on element
  /*ndofe=9;*/
  ndofe = 21; /* ??? */
  //  number of strain/stress components
  tncomp=3; /* ??? */
  //  number of functions approximated
  napfun=21; /* ??? */
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;
  //  order of numerical integration of mass matrix
  /*intordmm=3;*/
  intordmm = 6;
  //  strain/stress state
  ssst=platek;

  //  number of blocks (parts of geometric matrix)
  /*nb=2;*/
  nb = 1;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0] = 3; /* ??? */

  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0] = 3; /* ??? */
  
  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  
  nip[0][0]=12; /* ??? */
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=12; /* ??? */
}

ArgyrisTriangle::~ArgyrisTriangle (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] intordsm[i];
    delete [] nip[i];
  }
  delete [] intordsm;
  delete [] nip;

  delete [] cncomp;
  delete [] ncomp;
}


/**
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

   28.3.2002
*/
double ArgyrisTriangle::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function approximates function defined by nodal values

   @param xi,eta - natural coordinates
   @param nodval - nodal values

*/
double ArgyrisTriangle::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector areacoord(3);
  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  scprd (areacoord,nodval,f);
  return f;
}

void ArgyrisTriangle::bf_matrix (matrix& n, vector& x, vector& y, vector& areacoord)
{
	double pointX = 0.0, pointY = 0.0;

	double* xCoords = x.a;
	double* yCoords = y.a;
	double* pointAreaCoords = areacoord.a;
	
	convertAreaCoordinatesToNaturalCoordinates(xCoords, yCoords, pointAreaCoords, &pointX, &pointY);

	for (int i = 0; i < 3; i++)
	{
		n[0][i*7 + 0] = get_Hi(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 1] = get_Hix(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 2] = get_Hiy(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 3] = get_Hixx(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 4] = get_Hiyy(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 5] = get_Hixy(i, xCoords, yCoords, pointX, pointY);
		n[0][i*7 + 6] = get_HjkNi(i, xCoords, yCoords, pointX, pointY);
	}
}

void ArgyrisTriangle::geom_matrix (matrix& gm, vector& x, vector& y, vector& areacoord)
{
	double pointX = 0.0, pointY = 0.0;
	
	double* xCoords = x.a;
	double* yCoords = y.a;
	double* pointAreaCoords = areacoord.a;
	
	convertAreaCoordinatesToNaturalCoordinates(xCoords, yCoords, pointAreaCoords, &pointX, &pointY);

	for (int d = 0; d < 3; d++)
	{
		Derivatives deriv = (Derivatives)d;
		for (int i = 0; i < 3; i++)
		{
			gm[d][i*7 + 0] = get_Hi(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 1] = get_Hix(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 2] = get_Hiy(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 3] = get_Hixx(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 4] = get_Hiyy(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 5] = get_Hixy(i, deriv, xCoords, yCoords, pointX, pointY);
			gm[d][i*7 + 6] = get_HjkNi(i, deriv, xCoords, yCoords, pointX, pointY);
		}
	}
}


/**
   function assembles %matrix of approximation functions
   
   @param n - %matrix of approximation functions
   @param x,y - array containing node coordinates
   @param areacoord - array of area coordinates
   
   JK, 19.7.2001
*/
/*
void ArgyrisTriangle::bf_matrix (matrix &n,vector &x,vector &y,vector &areacoord)
{
  vector bf(9),sx(3),sy(3);
  
  sx[0]=x[1]-x[0];  sx[1]=x[2]-x[1];  sx[2]=x[0]-x[2];
  sy[0]=y[1]-y[0];  sy[1]=y[2]-y[1];  sy[2]=y[0]-y[2];
  
  bf_cct (bf.a,areacoord.a,sx.a,sy.a);
  
  fillm (0.0,n);
  
  n[0][0]=bf[0];
  n[0][1]=bf[1];
  n[0][2]=bf[2];
  n[0][3]=bf[3];
  n[0][4]=bf[4];
  n[0][5]=bf[5];
  n[0][6]=bf[6];
  n[0][7]=bf[7];
  n[0][8]=bf[8];
  
  n[1][1]=bf[0];
  n[1][4]=bf[3];
  n[1][7]=bf[6];

  n[2][2]=bf[0];
  n[2][5]=bf[3];
  n[2][8]=bf[6];

}
*/

/**
   function assembles geometric %matrix of the cct element
   
   @param gm - geometric %matrix
   @param x,y - array of node coordinates
   @param areacoord - area coordinates
   
   JK, 19.7.2001
*/
/*
void ArgyrisTriangle::geom_matrix (matrix &gm,vector &x,vector &y,vector &areacoord)
{
  double det;
  vector sx(3),sy(3),b(3),c(3),bf(9),dx(9),dy(9);
  
  sx[0]=x[1]-x[0];  sx[1]=x[2]-x[1];  sx[2]=x[0]-x[2];
  sy[0]=y[1]-y[0];  sy[1]=y[2]-y[1];  sy[2]=y[0]-y[2];

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);  

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  bf_cct (bf.a,areacoord.a,sx.a,sy.a);
  dx_cct (dx.a,areacoord.a,b.a,sx.a,sy.a);
  dy_cct (dy.a,areacoord.a,c.a,sx.a,sy.a);

  fillm (0.0,gm);
  
  gm[0][2]=dx[0];
  gm[0][5]=dx[3];
  gm[0][8]=dx[6];
  
  gm[1][1]=0.0-dy[0];
  gm[1][4]=0.0-dy[3];
  gm[1][7]=0.0-dy[6];
  
  gm[2][1]=0.0-dx[0];
  gm[2][2]=dy[0];
  gm[2][4]=0.0-dx[3];
  gm[2][5]=dy[3];
  gm[2][7]=0.0-dx[6];
  gm[2][8]=dy[6];
  
  gm[3][0]=dy[0];
  gm[3][1]=dy[1]-bf[0];
  gm[3][2]=dy[2];
  gm[3][3]=dy[3];
  gm[3][4]=dy[4]-bf[3];
  gm[3][5]=dy[5];
  gm[3][6]=dy[6];
  gm[3][7]=dy[7]-bf[6];
  gm[3][8]=dy[8];
  
  gm[4][0]=dx[0];
  gm[4][1]=dx[1];
  gm[4][2]=dx[2]+bf[0];
  gm[4][3]=dx[3];
  gm[4][4]=dx[4];
  gm[4][5]=dx[5]+bf[3];
  gm[4][6]=dx[6];
  gm[4][7]=dx[7];
  gm[4][8]=dx[8]+bf[6];
}
*/

void ArgyrisTriangle::dmat (matrix& d, double t)
{
  double c;
  
  c = t*t*t;
  
  d[0][0] = c*d[0][0];  d[0][1] = c*d[0][1];  d[0][2] = c*d[0][2];
  d[1][0] = c*d[1][0];  d[1][1] = c*d[1][1];  d[1][2] = c*d[1][2];
  d[2][0] = c*d[2][0];  d[2][1] = c*d[2][1];  d[2][2] = c*d[2][2];
}

/**
   nutna kontrola

   function assembles transformation %matrix 
*/
void ArgyrisTriangle::transf_matrix (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of cct element
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - vectors of node coordinates
   
   JK, 19.7.2001
*/
void ArgyrisTriangle::stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y)
{
  long i,ii,jj,ipp;
  double jac,det,thick;
  ivector nodes(nne);
  vector l(3),t(nne),gp1,gp2,w,areacoord(3);
  matrix gm,d(tncomp,tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);  

  fillm (0.0,sm);  

  ii=0;
  {
    reallocm (ncomp[ii],ndofe,gm);
    jj=0;
    if (intordsm[ii][jj] != 0)
    {
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      
      gauss_points_tr(gp1.a, gp2.a, w.a, intordsm[ii][jj]);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++)
      {
				areacoord[0]=gp1[i];
				areacoord[1]=gp2[i];
				areacoord[2]=1.0-areacoord[0]-areacoord[1];
	
				thick = approx (areacoord,t);

				//  geometric matrix
				geom_matrix (gm,x,y,areacoord);
	
				//  matrix of stiffness of the material
				Mm->matstiff(d, ipp);
				dmat(d, thick);
	
				jac = w[i]*det;
	
				//  contribution to the stiffness matrix of the element
				bdbjac(sm, gm, d, gm, jac);
	
				ipp++;
      }
    }
  }
}

/**
   function assembles resulting stiffness %matrix of the element
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 9.5.2002
*/
void ArgyrisTriangle::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);

  stiffness_matrix (eid,0,0,sm,x,y);
  
  fprintf (Out,"\n\n matice tuhosti \n");
  for (long i=0;i<21;i++){
    fprintf (Out,"\n");
    for (long j=0;j<21;j++){
      fprintf (Out,"\n r %2ld %2ld   %le",i+1,j+1,sm[i][j]);
    }
  }
  fprintf (Out,"\n\n konec matice tuhosti \n");
  
  //  transformation of stiffness %matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}

/**
   function computes mass %matrix of the cct element
   
   @param eid - element id
   @param mm - mass %matrix
   @param x,y - vectors of node coordinates
   
   JK, 19.7.2001
*/
void ArgyrisTriangle::mass_matrix (long eid,matrix &mm,vector &x,vector &y)
{
  long i;
  double det,ww,thick,rho;
  ivector nodes(nne);
  vector l(3),gp1(intordmm),gp2(intordmm),w(intordmm),t(nne),dens(nne);
  matrix n(1,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mc->give_density (eid,nodes,dens);
  
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,mm);

  for (i=0;i<intordmm;i++){
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
    ww=w[i];
    
    bf_matrix (n,x,y,l);
    
    thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];
    rho = dens[0]*l[0]+dens[1]*l[1]+dens[2]*l[2];
    
    ww*=det*thick*rho;

    nnj (mm.a,n.a,ww,1,ndofe);
  }
}

void ArgyrisTriangle::res_mass_matrix (long eid,matrix &mm)
{
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  mass_matrix (eid,mm,x,y);
}


/**
   function computes load %matrix of the cct element
   
   @param eid - element id
   @param lm - load %matrix
   
   25.7.2001
*/
void ArgyrisTriangle::load_matrix (long eid,matrix &lm)
{
  long i;
  double det,ww,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),l(3),gp1(intordmm),gp2(intordmm),w(intordmm),t(nne);
  matrix n(1,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);

  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordmm);

  fillm (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
    ww=w[i];
    
    bf_matrix (n,x,y,l);
    
    thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];

    ww*=det*thick;
    
    nnj (lm.a,n.a,ww,1,ndofe);
  }
}

/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements
   
   JK, 26.9.2008
*/
void ArgyrisTriangle::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ii,ipp;
  vector gp1,gp2,w,eps(tncomp),areacoord(3);
  matrix gm(tncomp,ndofe);
  
  //  loop over blocks /*note: nb == 1*/
  for (ii=0;ii<nb;ii++){
    
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      //  geometric matrix
      geom_matrix (gm,x,y,areacoord);
      
      mxv (gm,r,eps);
      
      Mm->storestrain (lcid,ipp,eps);
      ipp++;
    }
  }
}

/**
   function computes strains at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void ArgyrisTriangle::res_ip_strains (long lcid,long eid)
{
  vector aux,x(nne),y(nne),r(ndofe);
  ivector nodes(nne);
  matrix tmat;
  
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
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
  
  ip_strains (lcid,eid,0,0,x,y,r);
  
}

/**
   function computes strains in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 26.9.2008
*/
void ArgyrisTriangle::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelt (ipp,intordsm[0][0],ipnum);
  
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
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK, 26.9.2008
*/
void ArgyrisTriangle::res_ip_stresses (long lcid,long eid)
{
  long ri,ci;
  ri=0;
  ci=0;

  compute_nlstress (lcid,eid,ri,ci);
}

/**
   function computes stresses at nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 26.9.2008
*/
void ArgyrisTriangle::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelt (ipp,intordsm[0][0],ipnum);
  
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


void ArgyrisTriangle::strains (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
  /*
  long i,naep,ncp,sid;
  double **stra;
  vector coord,eps;
  
  if (Mp->strainaver==0){
    stra = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
    }
    elem_strains (stra,lcid,eid,ri,ci);
  }

  switch (Mm->stra.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    allip_strains (stra,lcid,eid,ri,ci);
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
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stra.give_aepcoord (sid,i,coord);

      if (Mp->strainaver==0)
	appval (coord[0],coord[1],0,ncp,eps,stra);
      if (Mp->strainaver==1)
	appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);
      
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
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
*/
}

/**
   function assembles natural coordinates of nodes of element
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   
   10.5.2002
*/
void ArgyrisTriangle::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 0.0;  eta[0] = 0.0;
  xi[1] = 1.0;  eta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 1.0;
}

/**
   function computes strains in arbitrary point on element
   
   @param xi, eta - natural coordinates of the point
   @param eps - array containing strains
   @param val - array containing values on element
   
   11.5.2002
*/
void ArgyrisTriangle::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval(nne),areacoord(3);
  
  areacoord[0]=xi;  areacoord[1]=eta;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  k=0;
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (areacoord,nodval);
    k++;
  }
}


void ArgyrisTriangle::stresses (long /*lcid*/,long /*eid*/,long /*ri*/,long /*ci*/)
{
/*
  long i,naep,ncp,sid;
  double **stra,**stre;
  vector coord,sig;
  
  if (Mp->stressaver==0){
    stra = new double* [nne];
    stre = new double* [nne];
    for (i=0;i<nne;i++){
      stra[i] = new double [tncomp];
      stre[i] = new double [tncomp];
    }
    elem_strains (stra,lcid,eid,ri,ci);
    elem_stresses (stra,stre,lcid,eid,ri,ci);
  }
  
  switch (Mm->stre.tape[eid]){
  case nowhere:{
    break;
  }
  case intpts:{
    allip_stresses (stre,lcid,eid,ri,ci);
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
    reallocv (2,coord);
    for (i=0;i<naep;i++){
      Mm->stre.give_aepcoord (sid,i,coord);

      if (Mp->stressaver==0)
	appval (coord[0],coord[1],0,ncp,sig,stre);
      if (Mp->stressaver==1)
	appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);
      
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
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
*/
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
void ArgyrisTriangle::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;

  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes internal forces for nonlocal models

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void ArgyrisTriangle::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=nonlocstress;

  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes increments of internal forces

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void ArgyrisTriangle::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=stressincr;

  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y);
}



/**
   function computes nodal forces caused by temperature changes
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y - nodal coordinates
   
   7.2008, TKo
*/
void ArgyrisTriangle::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
   function computes resulting internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void ArgyrisTriangle::res_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
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
void ArgyrisTriangle::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
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
void ArgyrisTriangle::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  incr_internal_forces (lcid,eid,0,0,ifor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
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
void ArgyrisTriangle::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  eigstrain_forces (lcid,eid,0,0,nfor,x,y);
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  Mt->give_elemnodes (eid,nodes);
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
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
   
   TKo, 7.2008
*/
void ArgyrisTriangle::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,ipp;
  double thick;
  ivector nodes(nne);
  vector sig(5),gp1,gp2,w,l(3),t(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1){
	l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
	thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];

        Mm->computenlstresses (ipp,Mm->ip[ipp]);

	Mm->givestress (lcid,ipp,sig);
	sig[0]*=thick*thick*thick;
	sig[1]*=thick*thick*thick;
	sig[2]*=thick*thick*thick;
	sig[3]*=thick*5.0/6.0;
	sig[4]*=thick*5.0/6.0;
	Mm->storestress (lcid,ipp,sig);
	
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
   
   TKo, 7.2008
*/
void ArgyrisTriangle::local_values(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of local values
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
   
   TKo, 7.2008
*/
void ArgyrisTriangle::compute_nonloc_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct stresses
      if (Mp->strcomp==1)
        Mm->compnonloc_nlstresses (ipp);
      ipp++;
    }
  }
}



/**
   function computes correct stress increments at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void ArgyrisTriangle::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      //  computation of correct increments of stresses
      if (Mp->strcomp==1)
        Mm->computenlstressesincr (ipp);
      ipp++;
    }
  }
}



/**
   function computes correct eigen stresses caused by temperature at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void ArgyrisTriangle::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
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
   
   TKo 7.2008
*/
void ArgyrisTriangle::elem_integration (integratedquant iq,long lcid,long eid,long /*ri*/,long /*ci*/,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double jac;
  ivector nodes(nne);
  vector w,gp1,gp2,areacoord(3);
  vector ipv,contr(ndofe);
  matrix gm;
  
  Mt->give_elemnodes (eid,nodes);
  fillv (0.0,nv);
  
  ii=0;
  {
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ii][ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      areacoord[0]=gp1[i];
      areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[2];

      jac_2d (jac,x,y,areacoord[0],areacoord[1]);
      reallocv (ncomp[ii],ipv);

      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);

      //  strain-displacement (geometric) matrix
      reallocm (ncomp[ii],ndofe,gm);
      geom_matrix (gm,x,y,areacoord);
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      cmulv (jac*w[i],contr);

      //  summation
      addv(contr,nv,nv);
      ipp++;
    }
  }
}



void ArgyrisTriangle::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  long nstra;
  double xi, eta, ipval;
  vector w, gp1, gp2, anv(nne);

  nstra = 0;
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
          ipval = approx_nat (xi,eta,anv);
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
    if (ictn[i] & inistrain) nstra++;
  }
}



// ==================================================================
// ==================================================================
// ==================================================================


// ==================================================================
// Common functions

inline double cube(double x)
{
	return x * x * x;
}

// ==================================================================


void ArgyrisTriangle::computeTriangleConstants(/*in*/ double* xCoords, /*in*/ double* yCoords, /*out*/ double* A, /*out*/ double* L, /*out*/ double* c1, /*out*/ double* c2)
{
	// compute twice the area of triangle
	*A = xCoords[0]*yCoords[1] + xCoords[1]*yCoords[2] + xCoords[2]*yCoords[0] - xCoords[0]*yCoords[2] - xCoords[1]*yCoords[0] - xCoords[2]*yCoords[1];

	// compute lengths of edges
	for (int i = 0; i < 3; i++)
	{
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;
		L[i] = sqrt(sqr(xCoords[j] - xCoords[k]) + sqr(yCoords[j] - yCoords[k]));
	}

	// compute coefficients for area coordinates
	c1[0] = (yCoords[1] - yCoords[2]) / *A; // a
	c1[1] = (yCoords[2] - yCoords[0]) / *A; // b
	c1[2] = (yCoords[0] - yCoords[1]) / *A; // c

	c2[0] = (xCoords[2] - xCoords[1]) / *A; // f
	c2[1] = (xCoords[0] - xCoords[2]) / *A; // g
	c2[2] = (xCoords[1] - xCoords[0]) / *A; // h
}

double ArgyrisTriangle::getEij(int i, double* L)
{
	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	return (sqr(L[k]) - sqr(L[i]) - sqr(L[j])) / (2 * sqr(L[j]));
}

double ArgyrisTriangle::getEik(int i, double* L)
{
	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	return (sqr(L[j]) - sqr(L[i]) - sqr(L[k])) / (2 * sqr(L[k]));
}

void ArgyrisTriangle::computeAreaCoordinatesForPoint(/*in*/ double x, /*in*/ double y, /*in*/ double* xCoords, /*in*/ double* yCoords, /*out*/ double* omega)
{
	double A = xCoords[0]*yCoords[1] + xCoords[1]*yCoords[2] + xCoords[2]*yCoords[0] - xCoords[0]*yCoords[2] - xCoords[1]*yCoords[0] - xCoords[2]*yCoords[1];

	omega[0] = ((yCoords[1] - yCoords[2]) / A) * x + ((xCoords[2] - xCoords[1]) / A) * y + (xCoords[1] * yCoords[2] / A) - (xCoords[2] * yCoords[1] / A);
	omega[1] = ((yCoords[2] - yCoords[0]) / A) * x + ((xCoords[0] - xCoords[2]) / A) * y + (xCoords[2] * yCoords[0] / A) - (xCoords[0] * yCoords[2] / A);
	omega[2] = ((yCoords[0] - yCoords[1]) / A) * x + ((xCoords[1] - xCoords[0]) / A) * y + (xCoords[0] * yCoords[1] / A) - (xCoords[1] * yCoords[0] / A);
}

void ArgyrisTriangle::convertAreaCoordinatesToNaturalCoordinates(/*in*/ double* xCoords, /*in*/ double* yCoords, /*in*/ double* omega, /*out*/ double* pointX, /*out*/ double* pointY)
{
	double x = 0.0;
	double y = 0.0;
	for (int i = 0; i < 3; i++)
	{
		x += xCoords[i] * omega[i];
		y += yCoords[i] * omega[i];
	}
	*pointX = x;
	*pointY = y;
}

// ==================================================================================
// base functions
// ==================================================================================

double ArgyrisTriangle::get_Hi(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;

	double result = cube(omega[i]) * (10.0 - 15.0 * omega[i] + 6.0 * sqr(omega[i])) - 30.0 * (getEij(i, L) * omega[k] + getEik(i, L) * omega[j]) * sqr(omega[i]) * omega[j] * omega[k];

	return result;
}

double ArgyrisTriangle::get_Hix(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double alpha = 7 * (xCoords[i] - xCoords[k]) * getEij(i, L) - 5 * (xCoords[j] - xCoords[k]);
	double beta  = 7 * (xCoords[i] - xCoords[j]) * getEik(i, L) - 5 * (xCoords[k] - xCoords[j]);
	
	double result = cube(omega[i]) * (4.0 - 3.0 * omega[i]) * (pointX - xCoords[i])
	+ alpha * sqr(omega[i]) * omega[j] * sqr(omega[k])
	+ beta * sqr(omega[i]) * sqr(omega[j]) * omega[k];
	
	return result;
}

double ArgyrisTriangle::get_Hiy(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double gamma = 7 * (yCoords[i] - yCoords[k]) * getEij(i, L) - 5 * (yCoords[j] - yCoords[k]);
	double delta = 7 * (yCoords[i] - yCoords[j]) * getEik(i, L) - 5 * (yCoords[k] - yCoords[j]);
	
	double result = cube(omega[i]) * (4.0 - 3.0 * omega[i]) * (pointY - yCoords[i])
	+ gamma * sqr(omega[i]) * omega[j] * sqr(omega[k])
	+ delta * sqr(omega[i]) * sqr(omega[j]) * omega[k];
	
	return result;
}

double ArgyrisTriangle::get_Hixx(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double epsilon = -(xCoords[i] - xCoords[k]) * ((xCoords[i] - xCoords[k]) * getEij(i, L) - 2 * (xCoords[j] - xCoords[k]));
	double phi     = -(xCoords[i] - xCoords[j]) * ((xCoords[i] - xCoords[j]) * getEik(i, L) - 2 * (xCoords[k] - xCoords[j]));
	
	double result = 0.5 * ( cube(omega[i])*sqr(pointX-xCoords[i]) + epsilon*sqr(omega[i])*omega[j]*sqr(omega[k]) + phi*sqr(omega[i])*sqr(omega[j])*omega[k] );

	return result;
}

double ArgyrisTriangle::get_Hiyy(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double psi = -(yCoords[i] - yCoords[k]) * ((yCoords[i] - yCoords[k]) * getEij(i, L) - 2 * (yCoords[j] - yCoords[k]));
	double rho = -(yCoords[i] - yCoords[j]) * ((yCoords[i] - yCoords[j]) * getEik(i, L) - 2 * (yCoords[k] - yCoords[j]));
	
	double result = 0.5 * ( cube(omega[i])*sqr(pointY-yCoords[i]) + psi*sqr(omega[i])*omega[j]*sqr(omega[k]) + rho*sqr(omega[i])*sqr(omega[j])*omega[k] );
	
	return result;
}

double ArgyrisTriangle::get_Hixy(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double mju = -0.5 * ((yCoords[i]-yCoords[k])*((xCoords[i]-xCoords[k])*getEij(i, L) - 2*(xCoords[j]-xCoords[k])) + (xCoords[i]-xCoords[k])*((yCoords[i]-yCoords[k])*getEij(i, L) - 2*(yCoords[j]-yCoords[k])));
	double pi  = -0.5 * ((yCoords[i]-yCoords[j])*((xCoords[i]-xCoords[j])*getEik(i, L) - 2*(xCoords[k]-xCoords[j])) + (xCoords[i]-xCoords[j])*((yCoords[i]-yCoords[j])*getEik(i, L) - 2*(yCoords[k]-yCoords[j])));

	double result = cube(omega[i])*(pointX-xCoords[i])*(pointY-yCoords[i]) + mju*sqr(omega[i])*omega[j]*sqr(omega[k]) + pi*sqr(omega[i])*sqr(omega[j])*omega[k];
	
	return result;
}

double ArgyrisTriangle::get_HjkNi(int i, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	
	double result = (-16.0 * A / L[i]) * omega[i]*sqr(omega[j])*sqr(omega[k]);

	return result;
}

// ==================================================================================
// derivatives of base functions
// ==================================================================================

double ArgyrisTriangle::get_Hi(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY) // Second derivative of base function Hi by x at point [x;y]. i = 0, 1, 2
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	double result = -60 * getEij(i, L) * getFirstCommonExpressionValue(i, dType, omega, c1, c2) - 60 * getEik(i, L) * getSecondCommonExpressionValue(i, dType, omega, c1, c2);
	switch (dType)
	{
	case dX2:
		result += 60 * sqr(c1[i]) * (omega[i] - 3 * sqr(omega[i]) + 2 * cube(omega[i]));		
		break;
	case dY2:
		result += 60 * sqr(c2[i]) * (omega[i] - 3 * sqr(omega[i]) + 2 * cube(omega[i]));
		break;
	case dXdY:
		result += 60 * c1[i] * c2[i] * (omega[i] - 3 * sqr(omega[i]) + 2 * cube(omega[i]));
		break;
	}
	return result;
}

double ArgyrisTriangle::get_Hix(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	double a = c1[i];
	double f = c2[i];

	double alpha = 7 * (xCoords[i] - xCoords[k]) * getEij(i, L) - 5 * (xCoords[j] - xCoords[k]);
	double beta  = 7 * (xCoords[i] - xCoords[j]) * getEik(i, L) - 5 * (xCoords[k] - xCoords[j]);
	double result = 2 * alpha * getFirstCommonExpressionValue(i, dType, omega, c1, c2) + 2 * beta * getSecondCommonExpressionValue(i, dType, omega, c1, c2);
	
	switch (dType)
	{
	case dX2:
		result += 24*sqr(a)*omega[i]*pointX + 24*a*sqr(omega[i]) - 36*sqr(a)*sqr(omega[i])*pointX - 24*a*cube(omega[i]) - 24*sqr(a)*xCoords[i]*omega[i] + 36*sqr(a)*xCoords[i]*sqr(omega[i]);
		break;
	case dY2:
		result += 12*sqr(f) * (2*pointX*omega[i] - 3*pointX*sqr(omega[i]) - 2*xCoords[i]*omega[i] + 3*xCoords[i]*sqr(omega[i]));
		break;
	case dXdY:
		result += 24*a*f*omega[i]*pointX + 12*f*sqr(omega[i]) - 36*a*f*sqr(omega[i])*pointX - 12*f*cube(omega[i]) - 24*a*f*xCoords[i]*omega[i] + 36*a*f*xCoords[i]*sqr(omega[i]);
		break;
	}
	return result;
}

double ArgyrisTriangle::get_Hiy(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	double a = c1[i];
	double f = c2[i];

	double gamma = 7 * (yCoords[i] - yCoords[k]) * getEij(i, L) - 5 * (yCoords[j] - yCoords[k]);
	double delta = 7 * (yCoords[i] - yCoords[j]) * getEik(i, L) - 5 * (yCoords[k] - yCoords[j]);
	double result = 2 * gamma * getFirstCommonExpressionValue(i, dType, omega, c1, c2) + 2 * delta * getSecondCommonExpressionValue(i, dType, omega, c1, c2);
	
	switch (dType)
	{
	case dX2:
		result += 24*sqr(a)*omega[i]*pointY - 36*sqr(a)*sqr(omega[i])*pointY - 24*sqr(a)*yCoords[i]*omega[i] + 36*sqr(a)*yCoords[i]*sqr(omega[i]);
		break;
	case dY2:
		result += 24*sqr(f)*omega[i]*pointY + 24*f*sqr(omega[i]) - 36*sqr(f)*sqr(omega[i])*pointY - 24*f*cube(omega[i]) - 24*sqr(f)*yCoords[i]*omega[i] + 36*sqr(f)*yCoords[i]*sqr(omega[i]);
		break;
	case dXdY:
		result += 24*a*f*omega[i]*pointY + 12*a*sqr(omega[i]) - 36*a*f*sqr(omega[i])*pointY - 12*a*cube(omega[i]) - 24*a*f*yCoords[i]*omega[i] + 36*a*f*yCoords[i]*sqr(omega[i]);
		break;
	}
	return result;
}

double ArgyrisTriangle::get_Hixx(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	double a = c1[i];
	double f = c2[i];

	double epsilon = -(xCoords[i] - xCoords[k]) * ((xCoords[i] - xCoords[k]) * getEij(i, L) - 2 * (xCoords[j] - xCoords[k]));
	double phi     = -(xCoords[i] - xCoords[j]) * ((xCoords[i] - xCoords[j]) * getEik(i, L) - 2 * (xCoords[k] - xCoords[j]));
	double result = 2 * epsilon * getFirstCommonExpressionValue(i, dType, omega, c1, c2) + 2 * phi * getSecondCommonExpressionValue(i, dType, omega, c1, c2);

	switch (dType)
	{
	case dX2:
		result += 6*sqr(a)*omega[i]*sqr(pointX-xCoords[i]) + 6*a*sqr(omega[i])*(pointX-xCoords[i]) + 6*sqr(omega[i])*a*(pointX-xCoords[i]) + 2*cube(omega[i]);
		break;
	case dY2:
		result += 6*sqr(f)*omega[i]*sqr(pointX-xCoords[i]);
		break;
	case dXdY:
		result += 6*a*f*omega[i]*sqr(pointX-xCoords[i]) + 6*f*sqr(omega[i])*(pointX-xCoords[i]);
		break;
	}
	result *= 0.5;
	return result;
}

double ArgyrisTriangle::get_Hiyy(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	double a = c1[i];
	double f = c2[i];

	double psi = -(yCoords[i] - yCoords[k]) * ((yCoords[i] - yCoords[k]) * getEij(i, L) - 2 * (yCoords[j] - yCoords[k]));
	double rho = -(yCoords[i] - yCoords[j]) * ((yCoords[i] - yCoords[j]) * getEik(i, L) - 2 * (yCoords[k] - yCoords[j]));
	double result = 2 * psi * getFirstCommonExpressionValue(i, dType, omega, c1, c2) + 2 * rho * getSecondCommonExpressionValue(i, dType, omega, c1, c2);

	switch (dType)
	{
	case dX2:
		result += 6*sqr(a)*omega[i]*sqr(pointY-yCoords[i]);
		break;
	case dY2:
		result += 6*sqr(f)*omega[i]*sqr(pointY-yCoords[i]) + 6*f*sqr(omega[i])*(pointY-yCoords[i]) + 6*f*sqr(omega[i])*(pointY-yCoords[i]) + 2*cube(omega[i]);
		break;
	case dXdY:
		result += 6*a*f*omega[i]*sqr(pointY-yCoords[i]) + 6*a*sqr(omega[i])*(pointY-yCoords[i]);
		break;
	}
	result *= 0.5;
	return result;
}

double ArgyrisTriangle::get_Hixy(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	int j = (i + 1) % 3;
	int k = (i + 2) % 3;
	double a = c1[i];
	double f = c2[i];

	double mju = -0.5 * ((yCoords[i]-yCoords[k])*((xCoords[i]-xCoords[k])*getEij(i, L) - 2*(xCoords[j]-xCoords[k])) + (xCoords[i]-xCoords[k])*((yCoords[i]-yCoords[k])*getEij(i, L) - 2*(yCoords[j]-yCoords[k])));
	double pi  = -0.5 * ((yCoords[i]-yCoords[j])*((xCoords[i]-xCoords[j])*getEik(i, L) - 2*(xCoords[k]-xCoords[j])) + (xCoords[i]-xCoords[j])*((yCoords[i]-yCoords[j])*getEik(i, L) - 2*(yCoords[k]-yCoords[j])));
	double result = 2 * mju * getFirstCommonExpressionValue(i, dType, omega, c1, c2) + 2 * pi * getSecondCommonExpressionValue(i, dType, omega, c1, c2);

	switch (dType)
	{
	case dX2:
		result += 6*sqr(a)*omega[i]*(pointX-xCoords[i])*(pointY-yCoords[i]) + 6*a*sqr(omega[i])*(pointY-yCoords[i]);
		break;
	case dY2:
		result += 6*sqr(f)*omega[i]*(pointX-xCoords[i])*(pointY-yCoords[i]) + 6*f*sqr(omega[i])*(pointX-xCoords[i]);
		break;
	case dXdY:
		result += 6*a*f*omega[i]*(pointX-xCoords[i])*(pointY-yCoords[i]) + 3*a*sqr(omega[i])*(pointX-xCoords[i]) + 3*f*sqr(omega[i])*(pointY-yCoords[i]) + cube(omega[i]);
		break;
	}

	return result;
}

double ArgyrisTriangle::get_HjkNi(int i, Derivatives dType, double* xCoords, double* yCoords, double pointX, double pointY)
{
	double A;
	double L[3];
	double omega[3];
	double c1[3];
	double c2[3];

	computeAreaCoordinatesForPoint(/*in*/ pointX, /*in*/ pointY, /*in*/ xCoords, /*in*/ yCoords, /*out*/ omega);	
	computeTriangleConstants(/*in*/ xCoords, /*in*/ yCoords, /*out*/ &A, /*out*/ L, /*out*/ c1, /*out*/ c2);

	double theta = -32 * A / L[i];
	return theta * getThirdCommonExpressionValue(i, dType, omega, c1, c2);
}

double ArgyrisTriangle::getFirstCommonExpressionValue(int i, Derivatives dType, double* omega, double* c1, double* c2)
{
	int j = (i + 1) % 3;
	int k = (i + 2) % 3;

	double* ca=NULL;
	double* cb=NULL;

	switch (dType)
	{
	case dX2:
		ca = cb = c1;
		break;
	case dY2:
		ca = cb = c2;
		break;
	case dXdY:
		ca = c1;
		cb = c2;
		break;
	}

	return	ca[i] * (cb[i] * omega[j] * sqr(omega[k]) + cb[j] * omega[i] * sqr(omega[k]) + 2 * cb[k] * omega[i] * omega[j] * omega[k])
		  + ca[j] * (cb[i] * omega[i] * sqr(omega[k]) + cb[k] * sqr(omega[i]) * omega[k])
		  + ca[k] * (2 * cb[i] * omega[i] * omega[j] * omega[k] + cb[j] * sqr(omega[i]) * omega[k] + cb[k] * sqr(omega[i]) * omega[j]);
}

double ArgyrisTriangle::getSecondCommonExpressionValue(int i, Derivatives dType, double* omega, double* c1, double* c2)
{
	int j = (i + 1) % 3;
	int k = (i + 2) % 3;

	double* ca=NULL;
	double* cb=NULL;

	switch (dType)
	{
	case dX2:
		ca = cb = c1;
		break;
	case dY2:
		ca = cb = c2;
		break;
	case dXdY:
		ca = c1;
		cb = c2;
		break;
	}

	return	ca[i] * (cb[i] * sqr(omega[j]) * omega[k] + 2 * cb[j] * omega[i] * omega[j] * omega[k] + cb[k] * omega[i] * sqr(omega[j]))
		  + ca[j] * (2 * cb[i] * omega[i] * omega[j] * omega[k] + cb[j] * sqr(omega[i]) * omega[k] + cb[k] * sqr(omega[i]) * omega[j])
		  + ca[k] * (cb[i] * omega[i] * sqr(omega[j]) + cb[j] * sqr(omega[i]) * omega[j]);
}

double ArgyrisTriangle::getThirdCommonExpressionValue(int i, Derivatives dType, double* omega, double* c1, double* c2)
{
	int j = (i + 1) % 3;
	int k = (i + 2) % 3;

	double* ca=NULL;
	double* cb=NULL;

	switch (dType)
	{
	case dX2:
		ca = cb = c1;
		break;
	case dY2:
		ca = cb = c2;
		break;
	case dXdY:
		ca = c1;
		cb = c2;
		break;
	}

	return ca[i] * (cb[j]*omega[j]*sqr(omega[k]) + cb[k]*sqr(omega[j])*omega[k])
		 + ca[j] * (cb[i]*omega[j]*sqr(omega[k]) + cb[j]*omega[i]*sqr(omega[k]) + 2*cb[k]*omega[i]*omega[j]*omega[k])
		 + ca[k] * (cb[i]*sqr(omega[j])*omega[k] + 2*cb[j]*omega[i]*omega[j]*omega[k] + cb[k]*omega[i]*sqr(omega[j]));
}

