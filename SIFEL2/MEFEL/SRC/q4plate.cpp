#include "q4plate.h"
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
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>


// w,fx,fy
q4plate::q4plate (void)
{
  long i,j;

  //  number nodes on element
  nne=4;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=5;
  //  number of functions approximated
  napfun=2;
  //  number of edges on element
  ned=4;
  //  number of nodes on one edge
  nned=2;
  intordmm=3;
  intordb=2;
  //  strain/stress state
  ssst=plates;
  
  //  number of blocks (parts of geometric matrix)
  nb=2;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=3;
  ncomp[1]=2;

  //  cumulative number of components approximated
  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=3;

  //  number of integration points
  //  order of numerical integration of stiffness matrix
  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=9;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=4;

  //  total number of integration points
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }

  intordsm[0][0]=3;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=2;
}



q4plate::~q4plate (void)
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

   @param xi,eta - coordinates on element
   @param nodval - nodal values
*/
double q4plate::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function assembles %matrix of transformation for Q4plate
   %matrix reduces 16 DOFs to 12 DOFs
   
   20.3.2002
*/
void q4plate::atd_matrix (matrix &atd,vector &x,vector &y,long /*eid*/)
{
  double sx,sy;
  matrix p(16,16),p1(16,16),ct(4,2);
  vector dl(4);
  long i,j,i1,i2,i3,i4,ii;

  dl[0]=sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1]=sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2]=sqrt((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
  dl[3]=sqrt((x[0]-x[3])*(x[0]-x[3])+(y[0]-y[3])*(y[0]-y[3]));
  // matrix vector of sides {cos,sin} 
  ct[0][0]=(x[1]-x[0])/dl[0];
  ct[0][1]=(y[1]-y[0])/dl[0];
  ct[1][0]=(x[2]-x[1])/dl[1];
  ct[1][1]=(y[2]-y[1])/dl[1];
  ct[2][0]=(x[3]-x[2])/dl[2];
  ct[2][1]=(y[3]-y[2])/dl[2];
  ct[3][0]=(x[0]-x[3])/dl[3];
  ct[3][1]=(y[0]-y[3])/dl[3];
  fillm (0.0,p);
  fillm (0.0,p1);
  fillm (0.0,atd);

  for (i=0;i<nne;i++){
    i1=3*i;
    i2=i1+1;
    i3=i1+2;
    i4=12+i;
    p[i1][0]=1.;
    p[i1][1]=x[i];
    p[i1][2]=y[i];
    p[i1][3]=x[i]*y[i];
    p[i1][4]=x[i]*x[i];
    p[i1][5]=y[i]*y[i];
    p[i1][6]=x[i]*x[i]*y[i];
    p[i1][7]=y[i]*y[i]*x[i];
    p[i1][8]=x[i]*x[i]*x[i];
    p[i1][9]=y[i]*y[i]*y[i];
    p[i1][10]=x[i]*x[i]*x[i]*y[i];
    p[i1][11]=y[i]*y[i]*y[i]*x[i];
    p[i2][5]=y[i]*2.;
    p[i2][6]=x[i]*x[i];
    p[i2][7]=y[i]*x[i]*2.;
    p[i2][9]=y[i]*y[i]*3.;
    p[i2][10]=x[i]*x[i]*x[i];
    p[i2][11]=y[i]*y[i]*x[i]*3.;
    p[i2][14]=1.;
    p[i2][15]=x[i];
    p[i3][4]=-x[i]*2.;
    p[i3][6]=-x[i]*y[i]*2.;
    p[i3][7]=-y[i]*y[i];
    p[i3][8]=-x[i]*x[i]*3.;
    p[i3][10]=-x[i]*x[i]*y[i]*3.;
    p[i3][11]=-y[i]*y[i]*y[i];
    p[i3][12]=-1.;
    p[i3][13]=-y[i];

    ii=i+1;

    if (ii > 3) ii=ii-4;
    
    sy=(y[i]+y[ii])/2.;
    sx=(x[i]+x[ii])/2.;
    p[i4][1]= ct[i][0];
    p[i4][2]= ct[i][1];
    p[i4][3]= sy*ct[i][0]+sx*ct[i][1];
    p[i4][12]=-ct[i][0];
    p[i4][13]=-sy*ct[i][0];
    p[i4][14]=-ct[i][1];
    p[i4][15]=-sx*ct[i][1];
  }

  double *ee;
  ee = new double [16*16];
  for (i=0;i<16;i++){
    for (j=0;j<16;j++){
      ee[i*16+j]=0.0;
    }
    ee[i*16+i]=1.0;
  }
  
  gemp (p.a,p1.a,ee,16,16,10.e-5,1);
//  invm (p,p1);  without move zerodiagonal
  delete [] ee;

/*  
  if(eid==122){
    fprintf (Out,"\n\n p1");
    for (i=0;i<16;i++){
      fprintf (Out,"\n");
      for (j=0;j<16;j++){
      fprintf (Out," %le",p1[i][j]);
	  }
    }
*/

  for (i=0;i<p1.n;i++){
    atd[i][0] =p1[i][0]  -p1[i][12]/dl[0]+p1[i][15]/dl[3];
    atd[i][1] =p1[i][1]  -p1[i][12]*ct[0][1]/2.-p1[i][15]*ct[3][1]/2.;
    atd[i][2] =p1[i][2]  +p1[i][12]*ct[0][0]/2.+p1[i][15]*ct[3][0]/2.;
    atd[i][3] =p1[i][3]  -p1[i][13]/dl[1]+p1[i][12]/dl[0];
    atd[i][4] =p1[i][4]  -p1[i][13]*ct[1][1]/2.-p1[i][12]*ct[0][1]/2.;
    atd[i][5] =p1[i][5]  +p1[i][13]*ct[1][0]/2.+p1[i][12]*ct[0][0]/2.;
    atd[i][6] =p1[i][6]  -p1[i][14]/dl[2]+p1[i][13]/dl[1];
    atd[i][7] =p1[i][7]  -p1[i][14]*ct[2][1]/2.-p1[i][13]*ct[1][1]/2.;
    atd[i][8] =p1[i][8]  +p1[i][14]*ct[2][0]/2.+p1[i][13]*ct[1][0]/2.;
    atd[i][9] =p1[i][9]  -p1[i][15]/dl[3]+p1[i][14]/dl[2];
    atd[i][10]=p1[i][10] -p1[i][15]*ct[3][1]/2.-p1[i][14]*ct[2][1]/2.;
    atd[i][11]=p1[i][11] +p1[i][15]*ct[3][0]/2.+p1[i][14]*ct[2][0]/2.;
  }
  
}

/**
  function assembles %matrix of base functions for w
  (approximation of deflection)

  20.3.2002
*/
void q4plate::bf_matrix (matrix &n,matrix &atd,vector &x,vector &y,vector &l)
{
  double sx,sy;
  int i;
  
  // first point is in I Q
  sx=(x[0]*(1+l[0])*(1+l[1])+x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
  sy=(y[0]*(1+l[0])*(1+l[1])+y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;
  // base w
  for (i=0;i<atd.n;i++){
    n[0][i]=atd[0][i]+sx*atd[1][i]+sy*atd[2][i]+sx*sy*atd[3][i]+sx*sx*atd[4][i]
      +sy*sy*atd[5][i]+sx*sx*sy*atd[6][i]+sy*sy*sx*atd[7][i]+sx*sx*sx*atd[8][i]
      +sy*sy*sy*atd[9][i]+sx*sx*sx*sy*atd[10][i]+sx*sy*sy*sy*atd[11][i];
  }
  
}

/**
   function assembles bending part of the geometrical %matrix of base functions Q4plate
   20.3.2002
*/
void q4plate::geom_matrix_bending (matrix &gm,matrix &atd,vector &x,vector &y,vector &l)
{
  double sx,sy;
  int i;
  // first point is in I Q
  sx=(x[0]*(1+l[0])*(1+l[1])+x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
  sy=(y[0]*(1+l[0])*(1+l[1])+y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;
  // original gm matrix, which would have to multiply by matrix atd
  //  gm[0][4]=-2.;  gm[0][6]=-y*2.;  gm[0][8]=-x*6.;  gm[0][10]=-x*y*6.;
  //  gm[1][5]=-2.;  gm[1][7]=-x*2.;  gm[1][9]=-y*6.;  gm[1][11]=-x*y*6.;
  //  gm[2][6]=-x*4.; gm[2][7]=-y*4.; gm[2][10]=-x*x*6.; gm[2][11]=-y*y*6.; gm[2][13]=-1.; gm[2][15]=-1.;
  
  for (i=0;i<atd.n;i++){
    gm[0][i]=-2.*atd[4][i]-2.*sy*atd[6][i]-6.*sx*atd[8][i]-6.*sx*sy*atd[10][i];
    gm[1][i]=-2.*atd[5][i]-2.*sx*atd[7][i]-6.*sy*atd[9][i]-6.*sx*sy*atd[11][i];
    gm[2][i]=-4.*sx*atd[6][i]-4.*sy*atd[7][i]-6.*sx*sx*atd[10][i]-6.*sy*sy*atd[11][i]-atd[13][i]-atd[15][i];
  }
  
}

/**
   function assembles shear part of the geometrical %matrix of base functions Q4plate
   20.3.2002
*/
void q4plate::geom_matrix_shear (matrix &gm,matrix &atd,vector &x,vector &y,vector &l)
{
  long i;
  double sx,sy;
  
  // first point is in I Q
  sx=(x[0]*(1+l[0])*(1+l[1])+x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
  sy=(y[0]*(1+l[0])*(1+l[1])+y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;
  
  for (i=0;i<atd.n;i++){
    gm[0][i]=atd[1][i]+sy*atd[3][i]-atd[12][i]-sy*atd[13][i];
    gm[1][i]=atd[2][i]+sx*atd[3][i]-atd[14][i]-sx*atd[15][i];
  }
}

/**
   function assembles stiffness matrices for bending and shear

   function extracts components of the bending
   part of stiffness %matrix of the material to the matrix db

   20.3.2002
*/
void q4plate::dmatblock (matrix &dd,matrix &d,long ri, long ci, double t)
{
  double c;
  
  fillm (0.0,dd);
  
  if (ri==0 && ci==0){
    c=t*t*t;
    dd[0][0] = c*d[0][0];  dd[0][1] = c*d[0][1];  dd[0][2] = c*d[0][2];
    dd[1][0] = c*d[1][0];  dd[1][1] = c*d[1][1];  dd[1][2] = c*d[1][2];
    dd[2][0] = c*d[2][0];  dd[2][1] = c*d[2][1];  dd[2][2] = c*d[2][2];
  }
  if (ri==1 && ci==1){
    dd[0][0] = t*d[3][3]*5.0/6.0;  dd[0][1] = 0.0;
    dd[1][0] = 0.0;                dd[1][1] = t*d[4][4]*5.0/6.0;
  }
}


/**
   function assembles transformation %matrix from local nodal coordinate
   system to the global coordinate system x_g = T x_l
   
   base vectors e1 and e2 have to contain 2 components
   the components are in the midplane
   
   @param nodes - element nodes
   @param tmat - transformation %matrix
   
   JK,
*/
void q4plate::transfmat (ivector &nodes,matrix &tmat)
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
   function computes stiffness %matrix of Q4plate element
   
   this function is used in plate elements (function is called
   by function res_stiffness_matrix) and shell elements
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param sm - stiffness %matrix
   @param x,y - nodal coordinates
   
   20.3.2002
*/
void q4plate::stiffness_matrix (long eid,long ri, long ci, matrix &sm,vector &x, vector &y)
{
  long i,j,ii,jj,ipp;
  double jac,a,a0,a1,thick;
  ivector nodes(nne);
  vector w,gp,l(2),t(nne);
  matrix gm,dd,atd(16,12),d(5,5);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  // translation to center it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  
  //     jakobian
  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);
  
  atd_matrix (atd,x,y,eid);
  
  fillm (0.0,sm);
  
  //  ii=Mt->elements[eid].ipp[sip][0];
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (intordsm[ii][jj],w); reallocv (intordsm[ii][jj],gp);
      reallocm (ncomp[jj],ndofe,gm);
      reallocm (ncomp[ii],ncomp[jj],dd);
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      
      for (i=0;i<intordsm[ii][jj];i++){
	l[0]=gp[i];  
	for (j=0;j<intordsm[ii][jj];j++){
	  l[1]=gp[j];
	  jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];
	  
	  if(ii==0) geom_matrix_bending (gm,atd,x,y,l);
	  else if(ii==1) geom_matrix_shear (gm,atd,x,y,l);
	  Mm->matstiff (d,ipp);
	  thick=approx (l[0],l[1],t);
          dmatblock (dd, d,ii,jj,thick);
          bdbjac (sm,gm,dd,gm,jac);
	  ipp++;
	}
      }
    }
  }
  
  
  /*
    fprintf (Out,"\n\n kontrola matice tuhosti");
    for (i=0;i<ndofe;i++){
    fprintf (Out,"\n");
    for (j=0;j<ndofe;j++){fprintf (Out," %le",sm[i][j]);}
    }
    if(eid<3)
    {
    fprintf (Out,"\n\n diagonalni prvky desky, %ld,\n",eid);
    for (i=0;i<ndofe;i++){fprintf (Out," %le",sm[i][i]);}
    }
  */  
  
}


/**
   function assembles stiffness %matrix of plane stress rectangular
   finite element with biquadratic approximation functions
   
   @param eid - element id
   @param sm - stiffness %matrix
   
*/
void q4plate::res_stiffness_matrix (long eid,matrix &sm)
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
    transfmat (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
}


/**
   function computes initial stress stiffness %matrix of Q4 element
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param ism - i-s matrix
   
   15.11.2003
*/
void q4plate::initstr_matrix (long eid,long /*ri*/,long /*ci*/,matrix &ism)
{
  long i,i1,i2,j,j1,j2,k,ii,jj;
  double jac,a,a0,a1,thick,xx,yy;
  ivector nodes(nne);
  vector w,gp,l(2),x(nne),y(nne),z(nne),sig(3),sx(3),sy(3),t(nne),dwx(12),dwy(12);
  matrix tran(3,3),atd(16,12);

// sig....... Axial stress

  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

// translation to center it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  
  //     jakobian
  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);

  atd_matrix (atd,x,y,eid);

  fillm (0.0,ism);

//  ii=Mt->elements[eid].ipp[sip][0];
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (intordsm[ii][jj],w); reallocv (intordsm[ii][jj],gp);

	  for (i=0;i<intordsm[ii][jj];i++){
		 l[0]=gp[i];  
		for (j=0;j<intordsm[ii][jj];j++){
		  l[1]=gp[j];
  	      xx=(x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
	      yy=(y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;
          for (k=0;k<dwx.n;k++){
           dwx[i]=atd[1][k]+yy*atd[3][k]+2.*xx*atd[4][k]+2.*xx*yy*atd[6][k]+
                  yy*yy*atd[7][k]+3.*xx*xx*atd[8][k]+3.*xx*xx*yy*atd[10][k]+yy*yy*yy*atd[11][k];
           dwy[i]=atd[2][k]+xx*atd[3][k]+2.*yy*atd[5][k]+xx*xx*atd[6][k]+2.*yy*xx*     
                  atd[7][k]+3.*yy*yy*atd[9][k]+xx*xx*xx*atd[10][k]+3.*xx*yy*yy*atd[11][k];
		  }  
  
	      thick=approx (l[0],l[1],t);
		  jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j]*thick/24.;
          for (i2=0;i2<sig.n;i2++){
            i1=3*i2-3;
            for (j2=0;j2<nne;j2++){
             j1=3*j2-3;
             ism[i1][j1]=ism[i1][j1]+( dwx[i]*dwx[j]*sig[0]
                +(dwx[i]*dwy[j]+dwx[j]*dwy[i])*sig[2]+dwy[i]*dwy[j]*sig[1] )*jac;
			}
		  }
		}
	  }
	}
  }
       

	
	
}



/**
   function computes load %matrix of the plate rectangular finite element
   load vector is obtained after premultiplying load %matrix
   by nodal load values
   
   @param eid - number of element
   @param lm - load %matrix
   
   JK, 3.10.2006
*/
void q4plate::load_matrix (long eid,matrix &lm)
{
  long i,j;
  double jac,xi,eta,w1,w2,a,a0,a1,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordmm),gp(intordmm),t(nne),l(2);
  matrix atd(16,12),n(napfun,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordmm);


  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;


  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);

  atd_matrix (atd,x,y,eid);


  fillm (0.0,lm);

  for (i=0;i<intordmm;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordmm;j++){
      eta=gp[j];  w2=w[j];
      
      bf_matrix (n,atd,x,y,gp);
      
      thick = approx (xi,eta,t);
      
      jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];

      jac*=thick;
      
      nnj (lm.a,n.a,jac,n.m,n.n);
    }
  }
  
}

void q4plate::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
{
  long i,j,k;
  vector nodval(nne);
  
  k=0;
  for (i=fi;i<fi+nc;i++){
    for (j=0;j<nne;j++){
      nodval[j]=val[j][i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
}
  






/**
   function computes strains Q4
   
   @param lcid - load case id
   @param eid - element id
*/
void q4plate::res_ip_strains (long lcid,long eid)
{
  vector aux,x(nne),y(nne),r(ndofe);
  ivector nodes(nne);
  matrix tmat;

  Mt->give_node_coord2d (x,y,eid);
  Mt->give_elemnodes (eid,nodes);
  eldispl (lcid,eid,r.a);

  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (ndofe,aux);
    reallocm (ndofe,ndofe,tmat);
    transfmat (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }

  ip_strains (lcid,eid,0,0,x,y,r);
}

/**
   function computes strains at integration points of element
   
   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y - %vectors of nodal coordinates
   @param r - %vector of nodal displacements
   
*/
void q4plate::ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ii,ipp;
  double a;
  vector gp,w,eps,aux,natcoord(2),l(3);
  matrix gm,tmat,atd(16,12);
  
  // translation to center it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  atd_matrix (atd,x,y,eid);
  
  //  fprintf (Out,"\n\n eps,Mainip prvek cislo %ld\n",eid);
  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;
    reallocv (intordsm[ii][ii],w); reallocv (intordsm[ii][ii],gp);
    reallocm (ncomp[ii],ndofe,gm); reallocv (ncomp[ii],eps);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];     
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    
    for (i=0;i<intordsm[ii][ii];i++){
      l[0]=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	l[1]=gp[j];
	if (ii==0)
	  geom_matrix_bending (gm,atd,x,y,l);
	if (ii==1)
	  geom_matrix_shear (gm,atd,x,y,l);

	mxv (gm,r,eps);
	Mm->storestrain (lcid,ipp,cncomp[ii],eps);
	ipp++;
	
	//          fprintf (Out,"\n");
	//		  for (ij=0;ij<ncomp[ii];ij++){
	//		    fprintf (Out,"%20.10e",eps[ij]);
	//			}
      }
    }
  }
  
}


/**
   function assembles strains at nodes of element
   strains are obtained from the nearest integration points

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices (default value is 0, nonzero values are used in shell elements)

   10.5.2002
*/
void q4plate::nod_strains_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector eps(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq (ipp,intordsm[0][0],ipnum);
  
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


void q4plate::strains (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
  double **stra=NULL;
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
    //allip_strains (stra,lcid,eid,ri,ci);
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
        //appstrain (lcid,eid,coord[0],coord[1],0,ncp,eps);
      Mm->stra.storevalues(lcid,eid,i,eps);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown strain point is required in function dstelem::strains (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  if (Mp->strainaver==0){
    for (i=0;i<nne;i++){
      delete [] stra[i];
    }
    delete [] stra;
  }

}



void q4plate::res_ip_stresses (long lcid,long eid)
{
  compute_nlstress (lcid,eid,0,0);
}



void q4plate::nod_stresses_ip (long lcid,long eid,long ri,long ci)
{
  long i,j,ipp;
  ivector ipnum(nne),nod(nne);
  vector sig(tncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Mt->elements[eid].ipp[ri][ci];
  nodip_planelq (ipp,intordsm[0][0],ipnum);
  
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

void q4plate::stresses (long lcid,long eid,long /*ri*/,long /*ci*/)
{
  long i,naep,ncp,sid;
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
    //allip_stresses (stre,lcid,eid,ri,ci);
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
	    //appstress (lcid,eid,coord[0],coord[1],0,ncp,sig);
      Mm->stre.storevalues(lcid,eid,i,sig);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown stress point is required in function planeelemrotlq::stresses (%s, line %d).\n",__FILE__,__LINE__);
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
   @param x,y - vectors of nodal coordinates

   TKo 7.2008
*/
void q4plate::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void q4plate::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void q4plate::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void q4plate::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y)
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
void q4plate::res_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  internal_forces (lcid,eid,0,0,ifor,x,y);
/*
  //  transformation of forces
  glvectortransf3dblock (ifor,tran); From LCS To GCS in shell
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transfmat (nodes,tmat);
    locglobtransf (ifor,rl,tmat);
  }
*/
}



/**
   function computes resulting internal forces for nonlocal models
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void q4plate::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y);
/*
  //  transformation of forces
  glvectortransf3dblock (ifor,tran); From LCS To GCS in shell
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transfmat (nodes,tmat);
    locglobtransf (ifor,rl,tmat);
  }
*/
}



/**
   function computes resulting increment of internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void q4plate::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  incr_internal_forces (lcid,eid,0,0,ifor,x,y);
/*
  //  transformation of forces
  glvectortransf3dblock (ifor,tran); From LCS To GCS in shell
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transfmat (nodes,tmat);
    locglobtransf (ifor,rl,tmat);
  }
*/
}



/**
   function computes resulting contributions from eigenstrains to the right hand side
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - %vector of internal forces

   TKo, 7.2008
*/
void q4plate::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
//  long transf;
  ivector nodes (nne);
  vector x(nne),y(nne);
  
  Mt->give_node_coord2d (x,y,eid);
  eigstrain_forces (lcid,eid,0,0,nfor,x,y);
/*
  //  transformation of forces
  glvectortransf3dblock (ifor,tran); From LCS To GCS in shell
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transfmat (nodes,tmat);
    locglobtransf (ifor,rl,tmat);
  }
*/
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void q4plate::compute_nlstress(long lcid,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  double thick;
  ivector nodes(nne);
  vector w,gp,t(nne),sig(tncomp);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],w);
      reallocv (intordsm[ii][jj],gp);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct stresses
	  if (Mp->strcomp==1){

	    thick=approx (gp[i],gp[j],t);

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
}



/**
   function computes correct increments of stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo  7.2008
*/
void q4plate::compute_nlstressincr(long /*lcid*/,long eid,long ri,long ci)
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
   function computes local values which will be used for averageing 
   in nonlocal models at integration points. Mp->nonlocphase have to be 1.

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void q4plate::local_values (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct stresses
	  if (Mp->strcomp==1)
	    Mm->computenlstresses (ipp,Mm->ip[ipp]);
	  ipp++;
	}
      }
    }
  }
}



/**
   function computes nonlocal correct stresses at integration points on element
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo 7.2008
*/
void q4plate::compute_nonloc_nlstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct nonlocal stresses
	  if (Mp->strcomp==1)
	    Mm->compnonloc_nlstresses (ipp);
	  ipp++;
	}
      }
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
void q4plate::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  vector eigstr(ASTCKVEC(tncomp)),sig(ASTCKVEC(tncomp));
  matrix d(ASTCKMAT(tncomp,tncomp));
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
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
void q4plate::elem_integration(integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,j,ii,jj,ipp;
  double jac,a,a0,a1,thick;
  ivector nodes(nne);
  vector w,gp,l(2),t(nne),ipv,f(ndofe);
  matrix gm,atd(16,12);
  
  // translation to center - it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  // jacobian
  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);
  atd_matrix (atd,x,y,eid);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  fillv (0.0,nv);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (intordsm[ii][jj],w); reallocv (intordsm[ii][jj],gp);
      reallocm (ncomp[jj],ndofe,gm); reallocv (ncomp[jj],ipv);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
	l[0]=gp[i];
	for (j=0;j<intordsm[ii][jj];j++){
	  l[1]=gp[j];
	  thick=approx (l[0],l[1],t);
          //  function assembles required quantity at integration point
          Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
          //  strain-displacement (geometric) matrix
	  if(ii==0) {
	    geom_matrix_bending (gm,atd,x,y,l);
	    //jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j]*thick*thick*thick;
	    jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];
	  }
	  else if(ii==1) {
	    geom_matrix_shear (gm,atd,x,y,l);
	    //jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j]*thick*5.0/6.0;
	    jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];
	  }
          //  contribution to the nodal values
	  mtxv (gm,ipv,f);
          cmulv (jac,f);

          //  summation
          addv(f,nv,nv);
/*	  for (k=0;k<f.n;k++){
	    ifor[k]+=f[k]*jac;
	    ifor[k]=0.0;                     //nekonsoliduje 
	  }	*/
	  ipp++;
	}
      }
    }
  }
}




/**
  function computes load vector from area forces fz of the Q4plate
  20.3.2002
*/
void q4plate::areaforces (long eid,double */*nv*/,vector &lm)
{
  long i,j;
  double jac,a,a0,a1,w1,w2;
  ivector nodes(nne);
  vector x(nne),y(nne),l(2),gp(intordmm),w(intordmm);
  matrix atd(16,12),n(1,ndofe),mm(ndofe,ndofe);

  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord2d (x,y,eid);
// translation to center it is necesary for w which in local coord. and no in relation
  a=(x[0]+x[1]+x[2]+x[3])/4.;
  x[0]=x[0]-a;
  x[1]=x[1]-a;
  x[2]=x[2]-a;
  x[3]=x[3]-a;
  a=(y[0]+y[1]+y[2]+y[3])/4.;
  y[0]=y[0]-a;
  y[1]=y[1]-a;
  y[2]=y[2]-a;
  y[3]=y[3]-a;
  
  //     jakobian
  a = (x[2]-x[0])*(y[3]-y[1])-(y[2]-y[0])*(x[3]-x[1]);
  a0=-(x[3]-x[2])*(y[1]-y[0])+(y[3]-y[2])*(x[1]-x[0]);
  a1=-(x[2]-x[1])*(y[3]-y[0])+(y[2]-y[1])*(x[3]-x[0]);
  
  gauss_points (gp.a,w.a,intordmm);
  atd_matrix (atd,x,y,eid);
  
  fillv (0.0,lm);
  
  for (i=0;i<intordmm;i++){
    for (j=0;j<intordmm;j++){
      l[0]=gp[i];  l[1]=gp[j];  w1=w[i];  w2=w[j];
      jac=(a+a0*l[0]+a1*l[1])*w1*w2;
      
      bf_matrix (n,atd,x,y,l);
      
      nnj (mm.a,n.a,jac,n.m,n.n);
      
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
void q4plate::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, l, ipp;
  long ii, jj, nv = nodval.n;
  double xi, eta, ipval;
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
            //  value in integration point
            ipval = approx (xi,eta,anv);
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
