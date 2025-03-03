#include "soilplateq.h"
#include "q4plate.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include <stdlib.h>
#include <math.h>

soilplateq::soilplateq (void)
{
  long i,j;

  nne=4;  ndofe=12;  tncomp=5;  napfun=2; ned=4;  nned=2;  
  intordmm=3; intordb=2; ssst=plates;
  nb=2; 
  
  ncomp = new long [nb];
  ncomp[0]=3;
  ncomp[1]=2;
  
  cncomp = new long [nb];
  cncomp[0]=0;
  cncomp[1]=3;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=9;  nip[0][1]=0;
  nip[1][0]=0;  nip[1][1]=4;
  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=3;  intordsm[0][1]=0;
  intordsm[1][0]=0;  intordsm[1][1]=2;
}

soilplateq::~soilplateq (void)
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
double soilplateq::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function assembles matrix of transformation for Q4plate
   20.3.2002
*/
void soilplateq::atd_matrix (matrix &atd,vector &x,vector &y)
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
  
  gemp (p.a,p1.a,ee,16,16,Mp->zero,1);
  delete [] ee;
  
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
   function extracts components of the soil
   part of stiffness matrix of the material to the matrix ds
   20.3.2002
*/
void soilplateq::dbmat (matrix &db,double c1, double c2, long ri, long ci)
{

  if (ri==0 && ci==0){
  db[0][0] = c1;   db[0][1] = 0.0;  db[0][2] = 0.0;
  db[1][0] = 0.0;  db[1][1] = c2;   db[1][2] = 0.0;
  db[2][0] = 0.0;  db[2][1] = 0.0;  db[2][2] = c2;
  }
  if (ri==1 && ci==1){
    db[0][0] = c2*5.0/6.0;		db[0][1] = 0.0;
    db[1][0] = 0.0;             db[1][1] = c2*5.0/6.0;
  }
}
void soilplateq::transfmat (ivector &nodes,matrix &tmat)
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

void soilplateq::appval (double xi,double eta,long fi,long nc,vector &eps,double **val)
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
  
void soilplateq::nodecoord (vector &xi,vector &eta)
{
  xi[0] = 1.0;  eta[0] = 1.0;
  xi[1] =-1.0;  eta[1] = 1.0;
  xi[2] =-1.0;  eta[2] =-1.0;
  xi[3] = 1.0;  eta[3] =-1.0;
}


void soilplateq::geom_matrix (matrix &gm,matrix &atd,vector &x,vector &y, vector &l)
{
  int i;
  double sx,sy;

// quadrilateral   w, dw/dx. dw/dy
  sx=(x[0]*(1+l[0])*(1+l[1])+x[1]*(1-l[0])*(1+l[1])+x[2]*(1-l[0])*(1-l[1])+x[3]*(1+l[0])*(1-l[1]))/4.;
  sy=(y[0]*(1+l[0])*(1+l[1])+y[1]*(1-l[0])*(1+l[1])+y[2]*(1-l[0])*(1-l[1])+y[3]*(1+l[0])*(1-l[1]))/4.;

  for (i=0;i<gm.n;i++){
    gm[0][i]=atd[0][i]+sx*atd[1][i]+sy*atd[2][i]+sx*sy*atd[3][i]+sx*sx*atd[4][i]+sy*sy*atd[5][i]+
      sx*sx*sy*atd[6][i]+sy*sy*sx*atd[7][i]+sx*sx*sx*atd[8][i]+sy*sy*sy*atd[9][i]+sx*sx*sx*sy*atd[10][i]+sx*sy*sy*sy*atd[11][i];
    gm[1][i]=atd[1][i]+sy*atd[3][i]+2.*sx*atd[4][i]+2.*sx*sy*atd[6][i]+sy*sy*atd[7][i]+3.*sx*sx*atd[8][i]+3.*sx*sx*sy*atd[10][i]+sy*sy*sy*atd[11][i];
//	gm[1][i]=0.0;
	gm[2][i]=atd[2][i]+sx*atd[3][i]+2.*sy*atd[5][i]+sx*sx*atd[6][i]+2.*sy*sx*atd[7][i]+3.*sy*sy*atd[9][i]+sx*sx*sx*atd[10][i]+3.*sx*sy*sy*atd[11][i];
//	gm[2][i]=0.0;
  }
}
void soilplateq::geom_matrix_shear (matrix &gm,matrix &atd,vector &x,vector &y,vector &l)
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

void soilplateq::res_stiffness_matrix (long eid,matrix &sm)
{
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  stiffness_matrix (eid,0,0,sm,x,y);
}

	  
void soilplateq::stiffness_matrix (long eid,long ri, long ci, matrix &sm,vector &x, vector &y)
{
  long ii,jj,i,j,ipp;
  double jac,a,a0,a1;
  ivector nodes(nne);
  vector w,gp,l(2);
  matrix gm,cc,atd(16,12),d(3,3);
//     jakobian
  Mt->give_elemnodes (eid,nodes);
  
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
  
  atd_matrix (atd,x,y);
  fillm (0.0,sm);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (d,ipp);
 
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      reallocv (intordsm[ii][jj],w); reallocv (intordsm[ii][jj],gp);
      reallocm (ncomp[jj],ndofe,gm);
      reallocm (ncomp[ii],ncomp[jj],cc);

      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      gauss_points (gp.a,w.a,intordsm[ii][jj]);

	  for (i=0;i<intordsm[ii][jj];i++){
		 l[0]=gp[i];  
		for (j=0;j<intordsm[ii][jj];j++){
		  l[1]=gp[j];
		  jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];

		  if(ii==0) geom_matrix (gm,atd,x,y,l);
		  else if(ii==1) geom_matrix_shear (gm,atd,x,y,l);
          dbmat (cc, d[0][0],d[1][1],ii,jj);
          bdbjac (sm,gm,cc,gm,jac);
		  ipp++;
		}
	  }
	}
 
  }


/*
  for (i=0;i<intordsm[0][0];i++){
	for (j=0;j<intordsm[0][0];j++){
	  l[0]=gp[i];  l[1]=gp[j];
	  jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];
	  geom_matrix (gm,atd,x,y,l);
      Mm->matstiff (cc,ipp);
//      imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
      bdbjac (sm,gm,cc,gm,jac);
	  ipp++;  
	}
  }
*/

//  fprintf (Out,"\n\n kontrola matice tuhosti");
//  for (i=0;i<ndofe;i++){fprintf (Out,"\n");
//    for (j=0;j<ndofe;j++){fprintf (Out," %le",sm[i][j]);}
//  }
//    fprintf (Out,"\n\n diagonalni prvky, %ld,\n",eid);
//    for (i=0;i<ndofe;i++){fprintf (Out," %le",sm[i][i]);}

}


void soilplateq::res_mainip_strains (long lcid,long eid)
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
    transfmat (nodes,tmat);
    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }

  mainip_strains (lcid,eid,0,0,x,y,r);
}


void soilplateq::mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,j,ii,ipp;
  double a;
  vector gp,w,eps,natcoord(2),l(3);
  matrix gm,atd(16,12);

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
  atd_matrix (atd,x,y);
  
//  fprintf (Out,"\n\n deformace");  fprintf (Out,"\n\n prvek cislo %ld",eid);
 for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocv (intordsm[ii][ii],gp);  reallocv (intordsm[ii][ii],w);
    reallocv (ncomp[ii],eps);
    reallocm (ncomp[ii],ndofe,gm);
    gauss_points (gp.a,w.a,intordsm[ii][ii]);

    for (i=0;i<intordsm[ii][ii];i++){
        l[0]=gp[i];
		for (j=0;j<intordsm[ii][ii];j++){
  		  l[1]=gp[j];
		  if(ii==0) geom_matrix (gm,atd,x,y,l);
		  else if(ii==1) geom_matrix_shear (gm,atd,x,y,l);
          mxv (gm,r,eps);
          Mm->storestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
          ipp++;
//          fprintf (Out,"\n"); for (ij=0;ij<ncomp[ii];ij++){fprintf (Out,"%20.10e",eps[ij]);}
		}
	}
  }
  
}

void soilplateq::elem_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp,w,eps,natcoord(2),l(3);
  ivector cn,nodes(nne);

  lsm = new double [9];
  nullv (lsm,9);

  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  

//  fprintf (Out,"\n\n deformace prvek cislo %ld",eid);
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocv (ncomp[ii],eps);
    reallocv (intordsm[ii][ii],gp);    reallocv (intordsm[ii][ii],w);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    nullv (rhs,ncomp[ii]*3);

    for (i=0;i<intordsm[ii][ii];i++){
        l[0]=gp[i];
		for (j=0;j<intordsm[ii][ii];j++){
  		  l[1]=gp[j];
	      Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
          natcoord[0]=l[0];  natcoord[1]=l[1];
          matassem_lsm (lsm,natcoord);
          rhsassem_lsm (rhs,natcoord,eps);
//          fprintf (Out,"\n");
//		  for (long ij=0;ij<ncomp[ii];ij++){
//		    fprintf (Out,"%20.10e",eps[ij]);
//			}
        ipp++;
		}
	}
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    nodal_values (stra,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii]);
    delete [] lhs;  delete [] rhs;
  }
  delete [] lsm;
}

void soilplateq::nod_strains (long lcid,long eid,long ri,long ci)
{
  long i,j,ii,ipp;
  double *lsm,*lhs,*rhs;
  vector nxi(nne),neta(nne),gp,w,eps,natcoord(2),l(3);
  ivector cn,nodes(nne);

  lsm = new double [9];

  nodecoord (nxi,neta);
  Mt->give_elemnodes (eid,nodes);
  

//  fprintf (Out,"\n\n deformace");
//  fprintf (Out,"\n\n prvek cislo %ld",eid);
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocv (ncomp[ii],eps);
    reallocv (intordsm[ii][ii],gp);    reallocv (intordsm[ii][ii],w);
    lhs = new double [ncomp[ii]*3];
    rhs = new double [ncomp[ii]*3];
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    nullv (lsm,9);
    nullv (rhs,ncomp[ii]*3);

    for (i=0;i<intordsm[ii][ii];i++){
        l[0]=gp[i];
		for (j=0;j<intordsm[ii][ii];j++){
  		  l[1]=gp[j];
	      Mm->givestrain (lcid,ipp,cncomp[ii],ncomp[ii],eps);
          natcoord[0]=l[0];  natcoord[1]=l[1];
          matassem_lsm (lsm,natcoord);
          rhsassem_lsm (rhs,natcoord,eps);
//          fprintf (Out,"\n");
//		  for (ij=0;ij<ncomp[ii];ij++){
//		    fprintf (Out,"%20.10e",eps[ij]);
//			}
        ipp++;
		}
	}
    solve_lsm (lsm,lhs,rhs,Mp->zero,3,ncomp[ii]);
    Mt->strain_nodal_values (nodes,nxi,neta,nxi,lhs,2,cncomp[ii],ncomp[ii],lcid);
    delete [] lhs;  delete [] rhs;
  }
  
  delete [] lsm;
}

void soilplateq::appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps)
{
  long i,j,k;
  ivector nodes;
  vector l(3),nodval;
  if (ncomp != eps.n){
    fprintf (stderr,"\n\n wrong interval of indices in function strain (%s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }


  reallocv (nne,nodes);
  reallocv (nne,nodval);
  Mt->give_elemnodes (eid,nodes);
  k=0;
  for (i=fi;i<fi+ncomp;i++){
    for (j=0;j<nne;j++){
      nodval[j]=Mt->nodes[nodes[j]].strain[lcid*tncomp+i];
    }
    eps[k]=approx (xi,eta,nodval);
    k++;
  }
  
}

void soilplateq::allip_strains (double **stra,long lcid,long eid,long ri,long ci)
{
  long i,j,ii,jj,ipp;
  vector eps(tncomp),gp,w,areacoord(3);
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      reallocv (intordsm[ii][jj],gp);  reallocv (intordsm[ii][jj],w);
      gauss_points (gp.a,w.a,intordsm[ii][jj]);
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
    	for (j=0;j<intordsm[ii][jj];j++){
    	if (Mp->strainaver==0)
		  appval (gp[i],gp[j],0,tncomp,eps,stra);
	    if (Mp->strainaver==1)
	      appstrain (lcid,eid,gp[i],gp[j],0,tncomp,eps);
          Mm->storestrain (lcid,ipp,eps);
	      ipp++;
		}
      }
    }
  }
}

void soilplateq::strains (long lcid,long eid,long ri,long ci)
{
  long i,naep,ncp,sid;
  double **stra=NULL;
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


void soilplateq::res_allip_stresses (long lcid,long eid)
{
   allip_stresses (lcid, eid, 0, 0);
}

/**
   function computes stresess
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   
   20.9.2006
*/
void soilplateq::allip_stresses (long lcid,long eid,long ri,long ci)
{
  long ipp,i,j;
  ivector nodes(nne);
  vector eps, sig;
  matrix cc;
  
  Mt->give_elemnodes (eid,nodes);
    
  // translation to center it is necesary for w which in local coord. and no in relation
//  elem_strains (stra,lcid,eid,ri,ci); 
  
  reallocv (ncomp[0],eps);  reallocv (ncomp[0],sig);
  reallocm (ncomp[0],ncomp[0],cc);
        
  ipp=Mt->elements[eid].ipp[ri][ci];
 
//  fprintf (Out,"\n\n Eps,sig soilprvek cislo %ld\n",eid);
  for (i=0;i<intordsm[0][0];i++){
	for (j=0;j<intordsm[0][0];j++){
      Mm->matstiff (cc,ipp);
	  Mm->givestrain (lcid,ipp,cncomp[0],ncomp[0],eps);
	  mxv (cc,eps,sig);	
	  Mm->storestress (lcid,ipp,cncomp[0],sig);      
//  for (k=0;k<3;k++){
//    printf (Out," %15.5e",eps[k]); fprintf (Out," %15.5e",sig[k]);  }
//    fprintf (Out,"\n");
	  ipp++;  
	}
  }

//  fprintf (Out,"\n\n Fint soil prvek cislo %ld\n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",ifor[k]);}
  
}



/**
   function computes internal forces
   
   this function is used in plane stress/strain elements (function is called
   by function res_internal_forces) and shell elements

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - node coordinates
   
   TKo, 7.2008
*/
void soilplateq::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}




/**
   function computes internal forces  for nonlocal models

   this function is used in plane stress/strain elements (function is called
   by function res_nonloc_internal_forces) and shell elements
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void soilplateq::nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=nonlocstress;
  
  //  computation of stresses
  compute_nonloc_nlstress (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes increment of  internal forces (from correct stresses increment)
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   @param x,y,z - vectors of nodal coordinates
   
   TKo, 7.2008
*/
void soilplateq::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z)
{
  integratedquant iq;
  iq=stressincr;
  
  //  computation of stresses
  compute_nlstressincr (lcid,eid,ri,ci);
  
  //  integration of stresses over the element
  elem_integration (iq,lcid,eid,ri,ci,ifor,x,y,z);
}



/**
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   this function is used in plane stress/strain elements (function is called
   by function res_eigstrain_forces) and shell elements

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param nfor - array containing nodal forces
   @param x,y,z - nodal coordinates
   
   TKo, 7.2008
*/
void soilplateq::eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z)
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
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - internal forces
   
   TKo, 7.2008
*/
void soilplateq::res_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  internal_forces (lcid,eid,0,0,ifor,x,y,z);
  
/*
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
*/
}



/**
   function computes resulting internal forces for nonlocal models
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void soilplateq::res_nonloc_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  nonloc_internal_forces (lcid,eid,0,0,ifor,x,y,z);

/*  
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
*/
}



/**
   function computes resulting increments of internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKo, 7.2008
*/
void soilplateq::res_incr_internal_forces (long lcid,long eid,vector &ifor)
{
//  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  incr_internal_forces (lcid,eid,0,0,ifor,x,y,z);

/*  
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
*/
}



/**
   function computes nodal forces caused by eigenstrains
   eigenstrain expresses e.g. temperature strains
   
   @param lcid - load case id
   @param eid - element id
   @param nfor - array containing nodal forces
   
   TKo, 7.2008
*/
void soilplateq::res_eigstrain_forces (long lcid,long eid,vector &nfor)
{
//  long transf;
  ivector nodes (nne);
  vector v(ndofe),x(nne),y(nne),z(nne);
  
  Mt->give_node_coord3d (x,y,z,eid);
  
  eigstrain_forces (lcid,eid,0,0,nfor,x,y,z);

/*  
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
*/
}



/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   TKo, 7.2008
*/
void soilplateq::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long ipp,ii,jj,i,j;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ii+ri][jj+ci];

      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct stresses
	  if (Mp->strcomp==1)
	    {
	      Mm->computenlstresses (ipp,Mm->ip[ipp]);
	      ipp++;  
	    }
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
void soilplateq::local_values (long lcid,long eid,long ri,long ci)
{
  long ipp,ii,i,j;
  vector eps(3),sig(3);
  matrix cc;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocm (ncomp[ii],ncomp[ii],cc); 
    for (i=0;i<intordsm[0][0];i++){
      for (j=0;j<intordsm[0][0];j++){
        //  computation of correct stresses
        if (Mp->strcomp==1)
	{
          Mm->computenlstresses (ipp,Mm->ip[ipp]);
          Mm->givestress (lcid,ipp,eps);  //strain for consolidation
          Mm->matstiff (cc,ipp);
	  mtxv (cc,eps,sig);	
          Mm->storestress (lcid,ipp,sig);  //stresses
	}
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
   
   TKo, 7.2008
*/
void soilplateq::compute_nonloc_nlstress (long lcid,long eid,long ri,long ci)
{
  long ipp,ii,i,j;
  vector eps(3),sig(3);
  matrix cc;
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocm (ncomp[ii],ncomp[ii],cc); 
    for (i=0;i<intordsm[0][0];i++){
      for (j=0;j<intordsm[0][0];j++){
        //  computation of correct stresses
        if (Mp->strcomp==1)
	{
          Mm->compnonloc_nlstresses (ipp);
          Mm->givestress (lcid,ipp,eps);  //strain for consolidation
          Mm->matstiff (cc,ipp);
	  mtxv (cc,eps,sig);	
          Mm->storestress (lcid,ipp,sig);  //stresses
	}
	ipp++;  
      }
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
void soilplateq::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long ipp,ii,jj,i,j;

  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;

      ipp=Mt->elements[eid].ipp[ii+ri][jj+ci];
      
      for (i=0;i<intordsm[ii][jj];i++){
	for (j=0;j<intordsm[ii][jj];j++){
	  //  computation of correct stresses
	  if (Mp->strcomp==1)
	    {
	      Mm->computenlstressesincr (ipp);
	      ipp++; 
	    }
	}
      }
    }
  }
}



void soilplateq::compute_eigstress(long /*lcid*/,long eid,long ri,long ci)
{
  long ipp,ii,i,j;
  vector eigstr(tncomp),sig(tncomp);
  matrix d(tncomp,tncomp);
  
  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    for (i=0;i<intordsm[0][0];i++){
      for (j=0;j<intordsm[0][0];j++){
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
   
   TKo, 7.2008
*/
void soilplateq::elem_integration (integratedquant iq, long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y, vector&z)
{
  long ipp,ii,i,j;
  double a,a0,a1,jac;
  ivector nodes(nne);
  vector f(ndofe),ipv(3);
  vector w,gp,l(2);
  matrix gm,cc,tran(3,3),atd(16,12);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
    
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
  atd_matrix (atd,x,y);
  fillv (0.0,nv);

  for (ii=0;ii<nb;ii++){
    ipp=Mt->elements[eid].ipp[ii+ri][ii+ci];
    reallocv (intordsm[ii][ii],gp);  reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    gauss_points (gp.a,w.a,intordsm[ii][ii]);

    for (i=0;i<intordsm[0][0];i++){
      for (j=0;j<intordsm[0][0];j++){
	l[0]=gp[i];  l[1]=gp[j];
	jac=(a+a0*l[0]+a1*l[1])/8.*w[i]*w[j];
        //  function assembles required quantity at integration point
        Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
        //  strain-displacement (geometric) matrix
	geom_matrix (gm,atd,x,y,l);
        //  contribution to the nodal values
	mtxv (gm,ipv,f);	
        cmulv (jac,f);
        //  summation
        addv(f,nv,nv);
	ipp++;  
    }
  }
 }
}
