#include "dst.h"
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

dstelem::dstelem (void)
{
  long i,j;

  //  number of nodes on element
  nne=3;
  //  number of DOFs on element
  ndofe=9;
  //  number of strain/stress components
  tncomp=5;
  //  number of functions approximated
  napfun=3;
  //  number of edges on element
  ned=3;
  //  number of nodes on one edge
  nned=2;  
  intordmm=3;
  //  order of numerical integration on element edges (boundaries)
  intordb=2;
  //  strain/stress state
  ssst=plates;
  
  //  number of blocks (parts of geometric matrix)
  nb=3; 
  
  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=3;
  ncomp[1]=3;
  ncomp[2]=2;
  
  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=0;
  cncomp[2]=3;

  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=1;  nip[0][1]=0;  nip[0][2]=0;
  nip[1][0]=0;  nip[1][1]=3;  nip[1][2]=0;
  nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=1;
  
  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=1;  intordsm[0][1]=0; intordsm[0][2]=0;
  intordsm[1][0]=0;  intordsm[1][1]=3; intordsm[1][2]=0;
  intordsm[2][0]=0;  intordsm[2][1]=0; intordsm[2][2]=1;
}

dstelem::~dstelem (void)
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

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

   JK, 23.9.2008
*/
double dstelem::approx (vector &areacoord,vector &nodval)
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

void dstelem::tran_matrix (matrix &a,matrix &ct,long eid)
{
  double det,dj,gg,thick;
  long i,ipp;
  ivector nodes(nne);
  vector x(3),y(3),ax(3),ay(3),dl(3),t(nne);
  matrix dbx(3,3),dby(3,3),b(2,3),ba(3,3),dd(3,3),d(tncomp,tncomp);
  
  //  det is equal to double area of the element
  // dl(0)--  x(0)-x(1), 
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  ipp=Mt->elements[eid].ipp[0][0];
  Mm->matstiff (d,ipp);
  Mc->give_thickness (eid,nodes,t);
  thick = (t[0]+t[1]+t[2])/3.0;
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
  dl[0] = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1] = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
  // matrix vector of sides {cos,sin} 
  ct[0][0]=(x[1]-x[0])/dl[0];
  ct[0][1]=(y[1]-y[0])/dl[0];
  ct[1][0]=(x[2]-x[1])/dl[1];
  ct[1][1]=(y[2]-y[1])/dl[1];
  ct[2][0]=(x[0]-x[2])/dl[2];
  ct[2][1]=(y[0]-y[2])/dl[2];
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/2./det;  ax[1]=(x[0]-x[2])/2./det;  ax[2]=(x[1]-x[0])/2./det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/2./det;  ay[1]=(y[2]-y[0])/2./det;  ay[2]=(y[0]-y[1])/2./det;
  // dl0= betwen node 1,2
  dj=ay[1]*ax[0]+ay[0]*ax[1];
  dbx[0][0]= 8.*ct[0][0]*ay[0]*ay[1];
  dbx[1][0]= 4.*ct[0][1]*dj;
  dbx[2][0]= 4.*ct[0][0]*dj +8.*ct[0][1]*ay[0]*ay[1];
  dby[0][0]= 4.*ct[0][0]*dj;
  dby[1][0]= 8.*ct[0][1]*ax[0]*ax[1];
  dby[2][0]= 8.*ct[0][0]*ax[0]*ax[1] +4.*ct[0][1]*dj;

  dj=ay[2]*ax[1]+ay[1]*ax[2];
  dbx[0][1]= 8.*ct[1][0]*ay[1]*ay[2];
  dbx[1][1]= 4.*ct[1][1]*dj;
  dbx[2][1]= 4.*ct[1][0]*dj +8.*ct[1][1]*ay[1]*ay[2];
  dby[0][1]= 4.*ct[1][0]*dj;
  dby[1][1]= 8.*ct[1][1]*ax[1]*ax[2];
  dby[2][1]= 8.*ct[1][0]*ax[1]*ax[2] +4.*ct[1][1]*dj;

  dj=ay[0]*ax[2]+ay[2]*ax[0];
  dbx[0][2]= 8.*ct[2][0]*ay[2]*ay[0];
  dbx[1][2]= 4.*ct[2][1]*dj;
  dbx[2][2]= 4.*ct[2][0]*dj +8.*ct[2][1]*ay[2]*ay[0];
  dby[0][2]= 4.*ct[2][0]*dj;
  dby[1][2]= 8.*ct[2][1]*ax[2]*ax[0];
  dby[2][2]= 8.*ct[2][0]*ax[2]*ax[0] +4.*ct[2][1]*dj;

  dmatblock (dd, d,0,0,thick);

  gg=dd[2][2]/thick/thick*10.;
  for (i=0;i<3;i++){
     b[0][i]=(dd[0][0]*dbx[0][i]+dd[2][0]*dby[0][i]
             +dd[0][1]*dbx[1][i]+dd[2][1]*dby[1][i]
             +dd[0][2]*dbx[2][i]+dd[2][2]*dby[2][i])/gg;
     b[1][i]=(dd[1][0]*dby[0][i]+dd[2][0]*dbx[0][i]
             +dd[1][1]*dby[1][i]+dd[2][1]*dbx[1][i]
             +dd[1][2]*dby[2][i]+dd[2][2]*dbx[2][i])/gg;
  }

// vypocet matice transformace alfa-parametru do koncovych posunu
  dbx[0][0]=-.5*dl[0];
  dbx[1][0]=-dl[1]*(ct[0][0]*ct[1][0]+ct[0][1]*ct[1][1])/6.;
  dbx[2][0]=-dl[2]*(ct[0][0]*ct[2][0]+ct[0][1]*ct[2][1])/6.;
  dbx[0][1]=-dl[0]*(ct[1][0]*ct[0][0]+ct[1][1]*ct[0][1])/6.;
  dbx[1][1]=-.5*dl[1];
  dbx[2][1]=-dl[2]*(ct[1][0]*ct[2][0]+ct[1][1]*ct[2][1])/6.;
  dbx[0][2]=-dl[0]*(ct[2][0]*ct[0][0]+ct[2][1]*ct[0][1])/6.;
  dbx[1][2]=-dl[1]*(ct[2][0]*ct[1][0]+ct[2][1]*ct[1][1])/6.;
  dbx[2][2]=-.5*dl[2];

  for (i=0;i<3;i++){
     dbx[i][0]=dbx[i][0]+dl[i]*(ct[i][0]*b[0][0]+ct[i][1]*b[1][0]);
     dbx[i][1]=dbx[i][1]+dl[i]*(ct[i][0]*b[0][1]+ct[i][1]*b[1][1]);
     dbx[i][2]=dbx[i][2]+dl[i]*(ct[i][0]*b[0][2]+ct[i][1]*b[1][2]);
  }

   det= dbx[0][0]*dbx[1][1]*dbx[2][2]+dbx[0][1]*dbx[1][2]*dbx[2][0]+dbx[0][2]*dbx[1][0]*dbx[2][1]
       -dbx[0][1]*dbx[1][0]*dbx[2][2]-dbx[0][0]*dbx[1][2]*dbx[2][1]-dbx[0][2]*dbx[1][1]*dbx[2][0];
   ba[0][0]=(dbx[1][1]*dbx[2][2]-dbx[1][2]*dbx[2][1])/det;
   ba[0][1]=(dbx[0][2]*dbx[2][1]-dbx[0][1]*dbx[2][2])/det;
   ba[0][2]=(dbx[0][1]*dbx[1][2]-dbx[0][2]*dbx[1][1])/det;
   ba[1][0]=(dbx[1][2]*dbx[2][0]-dbx[1][0]*dbx[2][2])/det;
   ba[1][1]=(dbx[0][0]*dbx[2][2]-dbx[0][2]*dbx[2][0])/det;
   ba[1][2]=(dbx[0][2]*dbx[1][0]-dbx[0][0]*dbx[1][2])/det;
   ba[2][0]=(dbx[1][0]*dbx[2][1]-dbx[1][1]*dbx[2][0])/det;
   ba[2][1]=(dbx[0][1]*dbx[2][0]-dbx[0][0]*dbx[2][1])/det;
   ba[2][2]=(dbx[0][0]*dbx[1][1]-dbx[0][1]*dbx[1][0])/det;

  for (i=0;i<3;i++){
     a[i][0]=(-ba[i][0]+ba[i][2]);
     a[i][1]=-(ba[i][0]*ct[0][1]*dl[0]+ba[i][2]*ct[2][1]*dl[2])/2.;
     a[i][2]= (ba[i][0]*ct[0][0]*dl[0]+ba[i][2]*ct[2][0]*dl[2])/2.;
     a[i][3]=( ba[i][0]-ba[i][1]);
     a[i][4]=-(ba[i][0]*ct[0][1]*dl[0]+ba[i][1]*ct[1][1]*dl[1])/2.;
     a[i][5]= (ba[i][0]*ct[0][0]*dl[0]+ba[i][1]*ct[1][0]*dl[1])/2.;
     a[i][6]=( ba[i][1]-ba[i][2]);
     a[i][7]=-(ba[i][1]*ct[1][1]*dl[1]+ba[i][2]*ct[2][1]*dl[2])/2.;
     a[i][8]= (ba[i][1]*ct[1][0]*dl[1]+ba[i][2]*ct[2][0]*dl[2])/2.;
  }
/*
  fprintf (Out,"\n\n prvek cislo %ld",eid);
  for (j=0;j<3;j++){
		fprintf (Out,"\n");
	for (i=0;i<ndofe;i++){
		fprintf (Out," %20.10le",a[j][i]);
	}
  }
*/
}


void dstelem::geom_matrix_bconst (matrix &gm,vector &x,vector &y)
{
  double det;
  long i,i1,i2;
  vector ax(3),ay(3);
  
  //  det is equal to double area of the element
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.;
  //   der dL(i)/dy
  ax[0]=(x[2]-x[1])/2./det;  ax[1]=(x[0]-x[2])/2./det;  ax[2]=(x[1]-x[0])/2./det;
  //   der   dL(i)/dx
  ay[0]=(y[1]-y[2])/2./det;  ay[1]=(y[2]-y[0])/2./det;  ay[2]=(y[0]-y[1])/2./det;

  for (i=0;i<3;i++){
     i2=3*i+2;
     i1=i2-1;
     gm[0][i2]+=ay[i];
     gm[1][i1]+=-ax[i];
     gm[2][i2]+=ax[i];
     gm[2][i1]+=-ay[i];
  }
  
}

void dstelem::geom_matrix_bending (matrix &gm,matrix &a,matrix &ct,vector &x,vector &y,vector &l)
{
  double det;
  long i,i2,i3;
  vector ax(3),ay(3);
  matrix ba(3,3);
  
  //  det is equal to double area of the element
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.;
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/2./det;  ax[1]=(x[0]-x[2])/2./det;  ax[2]=(x[1]-x[0])/2./det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/2./det;  ay[1]=(y[2]-y[0])/2./det;  ay[2]=(y[0]-y[1])/2./det;
  
  for (i=0;i<3;i++){
     i2=i+1;
	 if (i2>2){i2=i2-3;}
     i3=i+2;
	 if (i3>2){i3=i3-3;}
	 ba[0][i]= 4.*ct[i][0]*( l[i2]*ay[i]+l[i]*ay[i2]+ay[i3]/3.);
     ba[1][i]= 4.*ct[i][1]*( l[i2]*ax[i]+l[i]*ax[i2]+ax[i3]/3.);
     ba[2][i]= 4.*ct[i][0]*( l[i2]*ax[i]+l[i]*ax[i2]+ax[i3]/3.)+4.*ct[i][1]*( l[i2]*ay[i]+l[i]*ay[i2]+ay[i3]/3.);
  }
  
  for (i=0;i<9;i++){
     gm[0][i]=ba[0][0]*a[0][i]+ba[0][1]*a[1][i]+ba[0][2]*a[2][i];
     gm[1][i]=ba[1][0]*a[0][i]+ba[1][1]*a[1][i]+ba[1][2]*a[2][i];
     gm[2][i]=ba[2][0]*a[0][i]+ba[2][1]*a[1][i]+ba[2][2]*a[2][i];

  }
  
}

void dstelem::geom_matrix_shear (matrix &gs,matrix &a,matrix &ct,long eid)
{
  double det,dj,gg,thick;
  long i,ipp;
  ivector nodes(nne);
  vector x(3),y(3),ax(3),ay(3),t(nne);
  matrix dbx(3,3),dby(3,3),dd(3,3),b(2,3),d(tncomp,tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
  ipp=Mt->elements[eid].ipp[0][0];
  Mm->matstiff (d,ipp);
  Mc->give_thickness (eid,nodes,t);
  thick = (t[0]+t[1]+t[2])/3.;
  //  det is equal to double area of the element
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.;
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/2./det;  ax[1]=(x[0]-x[2])/2./det;  ax[2]=(x[1]-x[0])/2./det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/2./det;  ay[1]=(y[2]-y[0])/2./det;  ay[2]=(y[0]-y[1])/2./det;
  dj=ay[1]*ax[0]+ay[0]*ax[1];
  dbx[0][0]= 8.*ct[0][0]*ay[0]*ay[1];
  dbx[1][0]= 4.*ct[0][1]*dj;
  dbx[2][0]= 4.*ct[0][0]*dj +8.*ct[0][1]*ay[0]*ay[1];
  dby[0][0]= 4.*ct[0][0]*dj;
  dby[1][0]= 8.*ct[0][1]*ax[0]*ax[1];
  dby[2][0]= 8.*ct[0][0]*ax[0]*ax[1] +4.*ct[0][1]*dj;

  dj=ay[2]*ax[1]+ay[1]*ax[2];
  dbx[0][1]= 8.*ct[1][0]*ay[1]*ay[2];
  dbx[1][1]= 4.*ct[1][1]*dj;
  dbx[2][1]= 4.*ct[1][0]*dj +8.*ct[1][1]*ay[1]*ay[2];
  dby[0][1]= 4.*ct[1][0]*dj;
  dby[1][1]= 8.*ct[1][1]*ax[1]*ax[2];
  dby[2][1]= 8.*ct[1][0]*ax[1]*ax[2] +4.*ct[1][1]*dj;

  dj=ay[0]*ax[2]+ay[2]*ax[0];
  dbx[0][2]= 8.*ct[2][0]*ay[2]*ay[0];
  dbx[1][2]= 4.*ct[2][1]*dj;
  dbx[2][2]= 4.*ct[2][0]*dj +8.*ct[2][1]*ay[2]*ay[0];
  dby[0][2]= 4.*ct[2][0]*dj;
  dby[1][2]= 8.*ct[2][1]*ax[2]*ax[0];
  dby[2][2]= 8.*ct[2][0]*ax[2]*ax[0] +4.*ct[2][1]*dj;

  dmatblock (dd, d,0,0,thick);
  gg=dd[2][2]/thick/thick*10.;
  for (i=0;i<3;i++){
     b[0][i]=(dd[0][0]*dbx[0][i]+dd[2][0]*dby[0][i]
              +dd[0][1]*dbx[1][i]+dd[2][1]*dby[1][i]
              +dd[0][2]*dbx[2][i]+dd[2][2]*dby[2][i])/gg;
     b[1][i]=(dd[1][0]*dby[0][i]+dd[2][0]*dbx[0][i]
              +dd[1][1]*dby[1][i]+dd[2][1]*dbx[1][i]
              +dd[1][2]*dby[2][i]+dd[2][2]*dbx[2][i])/gg;
  }
  for (i=0;i<9;i++){
     gs[0][i]=b[0][0]*a[0][i]+b[0][1]*a[1][i]+b[0][2]*a[2][i];
     gs[1][i]=b[1][0]*a[0][i]+b[1][1]*a[1][i]+b[1][2]*a[2][i];
  }

}

  
  
void dstelem::transf_matrix (ivector &nodes,matrix &tmat)
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
void dstelem::dmatblock (matrix &dd,matrix &d,long ri, long ci, double t)
{
  double c;
  
  c=t*t*t;
  fillm (0.0,dd);

  if ((ri==0 && ci==0)||(ri==1 && ci==1)){
  dd[0][0] = c*d[0][0];  dd[0][1] = c*d[0][1];  dd[0][2] = c*d[0][2];
  dd[1][0] = c*d[1][0];  dd[1][1] = c*d[1][1];  dd[1][2] = c*d[1][2];
  dd[2][0] = c*d[2][0];  dd[2][1] = c*d[2][1];  dd[2][2] = c*d[2][2];
  }
  if (ri==2 && ci==2){
    dd[0][0] = t*d[3][3]*5.0/6.0;  dd[0][1] = 0.0;
    dd[1][0] = 0.0;                dd[1][1] = t*d[4][4]*5.0/6.0;
  }
}

/**
   function computes stiffness %matrix of dst element
   
   @param eid - element id
   @param sm - stiffness %matrix

   15.3.2002
*/
void dstelem::stiffness_matrix (long eid,long ri,long ci, matrix &sm,vector &x, vector &y)
{
  long i,ii,jj,ipp;
  double jac,det,thick;
  ivector nodes(nne);
  vector w,gp1,gp2,l(3),t(nne);
  matrix gmc,dd,a(3,9),ct(3,2),d(tncomp,tncomp);

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  tran_matrix (a,ct,eid);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  fillm (0.0,sm);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1); reallocv (intordsm[ii][jj],gp2);
      reallocm (ncomp[ii],ndofe,gmc);
      reallocm (ncomp[ii],ncomp[jj],dd);

      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][jj]);

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];

      for (i=0;i<intordsm[ii][jj];i++){
	fillm (0.0,gmc);
	l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
	//	l[0]=0.5;  l[1]=0.5;  l[2]=1.0-l[0]-l[1];
	//	integr.3   l[0]=0.166; l[1]=0.166;  l[0]=0.66; l[1]=0.166;  l[0]=0.166; l[1]=0.666;

	if (ii==0){
	  geom_matrix_bconst (gmc,x,y);
	}
	if (ii==1){
	  geom_matrix_bending (gmc,a,ct,x,y,l);
	}
	if(ii==2){
	  geom_matrix_shear (gmc,a,ct,eid);
	}
	
	Mm->matstiff (d,ipp);
	ipp++;
	thick = t[0]*l[0]+t[1]*l[1]+t[2]*l[2];
	//  thick=0.7;
	dmatblock (dd, d,ii,jj,thick);
	jac=det*w[i];
	bdbjac (sm,gmc,dd,gmc,jac);
	
	//	   fprintf (Out,"\n %ld   %lf",ipp,sm[0][0]);
	
      }
    }
  }
}

void dstelem::res_stiffness_matrix (long eid,matrix &sm)
{
  long transf;
  ivector nodes(nne);
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);

  stiffness_matrix (eid,0,0,sm,x,y);

  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}
void dstelem::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 0.0;  eta[0] = 0.0;
  xi[1] = 1.0;  eta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 1.0;
}

/**
   function computes strains DST
*/
void dstelem::appval (vector &l, long fi,long nc,vector &eps,double **val)
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

*/
void dstelem::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long ii,i,ipp;
  vector eps,w,gp1,gp2,l(3);
  matrix gm,a(3,9),ct(3,2);
  
  tran_matrix (a,ct,eid);
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],eps);
    reallocv (intordsm[ii][ii],w);
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
      
      if (ii==0){
	//  first and second strain contributions have to be added together
	//  therefore, the first block is added to the second one
	geom_matrix_bending (gm,a,ct,x,y,l);
	geom_matrix_bconst (gm,x,y);
      }
      if (ii==1){
	geom_matrix_bending (gm,a,ct,x,y,l);
	geom_matrix_bconst (gm,x,y);
      }
      if(ii==2)
	geom_matrix_shear (gm,a,ct,eid);
      
      mxv (gm,r,eps);
      Mm->storestrain (lcid,ipp,cncomp[ii],eps);
      ipp++;
    }
  }
}


void dstelem::res_ip_strains (long lcid,long eid)
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
   function computes stresses at integration points
   
   @param lcid - load case id
   @param eid - element id
   
   JK
*/
void dstelem::res_ip_stresses (long lcid,long eid)
{
  compute_nlstress (lcid,eid,0,0);
}

void dstelem::stresses (long /*lcid*/,long /*eid*/,long /*ri*/, long /*ci*/)
{

}









/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.9.2008
*/
void dstelem::compute_nlstress (long lcid,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  double thick;
  ivector nodes(nne);
  vector sig(5),gp1,gp2,w,l(3),t(nne);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp1);
      reallocv (intordsm[ii][jj],gp2);
      gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);
      
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
	  sig[3]*=thick*5.0/6.0;
	  sig[4]*=thick*5.0/6.0;
	  Mm->storestress (lcid,ipp,sig);
	  
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

   JK, 24.9.2008
*/
void dstelem::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
   
   JK, 24.9.2008
*/
void dstelem::res_internal_forces (long lcid,long eid,vector &ifor)
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
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 23.9.2008
*/
void dstelem::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
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
   function computes increment of  internal forces (from correct stresses increment)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void dstelem::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void dstelem::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void dstelem::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double jac,thick,det;
  ivector nodes(nne);
  vector w,gp1,gp2,t(nne),ipv,contr(ndofe),l(3);
  matrix gm,a(3,9),ct(3,2);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  tran_matrix (a,ct,eid);
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],ipv);
    
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
//      if (ii==0){ 
//	geom_matrix_bending (gm,a,ct,x,y,l);
//	geom_matrix_bconst (gm,x,y);
//      }
      if(ii==1){
	geom_matrix_bending (gm,a,ct,x,y,l);
	geom_matrix_bconst (gm,x,y);
      }
      if(ii==2){
	geom_matrix_shear (gm,a,ct,eid);
      }
      
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
  function computes load vector from edge forces fz of the DST
  15.3.2002
*/
void dstelem::nodeforces (long eid,long *le,double *nv,vector &nf)
{
  double cn0,sn0,cn1,sn1,cn2,sn2;
  ivector nodes(nne);
  vector x(nne),y(nne),dl(3);
  
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
    nf[4]+=-(3.*nv[3]+2.*nv[6])*cn1;
    nf[5]+=-(3.*nv[3]+2.*nv[6])*sn1;
    nf[7]+= (2.*nv[3]+3.*nv[6])*cn1;
    nf[8]+= (2.*nv[3]+3.*nv[6])*sn1;
    nf[3]+=((2.-0.1)*nv[3] +(1.+0.1)*nv[6])*dl[1]/6.;
    nf[6]+=((1.-0.1)*nv[3] +(2.+0.1)*nv[6])*dl[1]/6.;
  }
  if (le[2]==1){
    cn2= (y[0]-y[2])*dl[2]/60.;
    sn2=-(x[0]-x[2])*dl[2]/60.;
    nf[7]+=-(3.*nv[6]+2.*nv[0])*cn2;
    nf[8]+=-(3.*nv[6]+2.*nv[0])*sn2;
    nf[1]+= (2.*nv[6]+3.*nv[0])*cn2;
    nf[2]+= (2.*nv[6]+3.*nv[0])*sn2;
    nf[6]+=((2.-0.1)*nv[6] +(1.+0.1)*nv[0])*dl[2]/6.;
    nf[0]+=((1.-0.1)*nv[6] +(2.+0.1)*nv[0])*dl[2]/6.;
  }
}


/**
  function computes load vector from area forces fz of the DST
  15.3.2002
*/
void dstelem::areaforces (long eid,double *nv,vector &nf)
{
   double pl,cn0,cn1,cn2,sn0,sn1,sn2;
   ivector nodes(nne);
   vector x(nne),y(nne);

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
void dstelem::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double ipval;
  vector w, gp1, gp2, anv(nne), l(ASTCKVEC(3));
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
          l[0]=gp1[i];
          l[1]=gp2[i];
          l[2]=1.0-l[0]-l[1];
          //  value in integration point
          ipval = approx(l,anv);
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
