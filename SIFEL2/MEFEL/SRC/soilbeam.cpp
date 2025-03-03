#include "soilbeam.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include "beamel3d.h"
#include "intpoints.h"
#include "elastisomat.h"
#include <stdlib.h>
#include <math.h>

soilbeam::soilbeam (void)
{
  long i;
  
  //  number nodes on element
  nne=2;
  //  number of DOFs on element
  ndofe=12;
  //  number of strain/stress components
  tncomp=6;
  //  number of functions approximated
  napfun=3;
  //  order of numerical integration of mass matrix
  intordmm=4;
  //  order of numerical integration of initial stiffness matrix
  intordism=2;
  //  strain/stress state
  ssst=spacebeam;

  //  number of blocks (parts of geometric matrix)
  nb=1;

  //  number of strain/stress components
  ncomp = new long [nb];
  ncomp[0]=1;
  
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
  nip[0][0]=2;
  intordsm[0][0]=2;

  c1= new double [3];
  c2= new double [3];

  
  bPod=1.;
  tnip=2;

}

soilbeam::~soilbeam (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
  }
  delete [] nip;
}


void soilbeam::transf_matrix (ivector &nodes,matrix &tmat)
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
void soilbeam::beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid)
{
  double c;
  fillm (0.0,tmat);
  
// tmat[0][0] = c;    tmat[0][1] = -1.0*s;    tmat[0][2] = 0.0;
// tmat[1][0] = s;    tmat[1][1] = c;         tmat[1][2] = 0.0;
// tmat[2][0] = 0.0;  tmat[2][1] = 0.0;       tmat[2][2] = 1.0;

// tmat[3][3] = c;    tmat[3][4] = -1.0*s;    tmat[3][5] = 0.0;
// tmat[4][3] = s;    tmat[4][4] = c;         tmat[4][5] = 0.0;
// tmat[5][3] = 0.0;  tmat[5][4] = 0.0;       tmat[5][5] = 1.0;
   tmat[0][0]=x[1]-x[0];
   tmat[1][0]=y[1]-y[0];
   tmat[2][0]=z[1]-z[0];
   dl=sqrt((tmat[0][0]*tmat[0][0]+tmat[1][0]*tmat[1][0]+tmat[2][0]*tmat[2][0]));
   if (dl<Mp->zero){
    fprintf (stderr,"\n\n zero length of the %ld beamel3d element",eid);
    fprintf (stderr,"\n in function beamel3d::beam_transf_matrix.\n");
   }
   tmat[0][0]=tmat[0][0]/dl;
   tmat[1][0]=tmat[1][0]/dl;
   tmat[2][0]=tmat[2][0]/dl;

   tmat[0][1] =vec[1]*tmat[2][0]-vec[2]*tmat[1][0];
   tmat[1][1] =vec[2]*tmat[0][0]-vec[0]*tmat[2][0];
   tmat[2][1] =vec[0]*tmat[1][0]-vec[1]*tmat[0][0];
   c=sqrt((tmat[0][1]*tmat[0][1]+tmat[1][1]*tmat[1][1]+tmat[2][1]*tmat[2][1]));
   if (c<Mp->zero){
    fprintf (stderr,"\n\n zero vec of the %ld beamel3d element",eid);
    fprintf (stderr,"\n in function beamel3d::beam_transf_matrix.\n");
   }
   tmat[0][1]=tmat[0][1]/c;
   tmat[1][1]=tmat[1][1]/c;
   tmat[2][1]=tmat[2][1]/c;

   tmat[0][2]=tmat[1][0]*tmat[2][1]-tmat[2][0]*tmat[1][1];
   tmat[1][2]=tmat[2][0]*tmat[0][1]-tmat[0][0]*tmat[2][1];
   tmat[2][2]=tmat[0][0]*tmat[1][1]-tmat[1][0]*tmat[0][1];
   c=sqrt((tmat[0][2]*tmat[0][2]+tmat[1][2]*tmat[1][2]+tmat[2][2]*tmat[2][2]));
   if (c<Mp->zero){
    fprintf (stderr,"\n\n zero vec of the %ld beamel3d element",eid);
    fprintf (stderr,"\n in function beamel2d::beam_transf_matrix.\n");
   }
   tmat[0][2]=tmat[0][2]/c;
   tmat[1][2]=tmat[1][2]/c;
   tmat[2][2]=tmat[2][2]/c;
}
void soilbeam::geom_matrix (matrix &n,double s,double dl,double gy,double gz)
{
  double ll,aj1;
  ll=dl*dl;
//  u
        n[0][0] = 1.-s;
        n[0][6] = s;
//  w
        aj1=1./(1.+2.*gy);
        n[2][2] =aj1*(1.+2.*gy-2.*gy*s-3.*s*s+2*s*s*s);
        n[2][4] =aj1*(-(1.+gy)*s+(2.+gy)*s*s-s*s*s)*dl;
        n[2][8] =aj1*(2.*gy*s+3.*s*s-2.*s*s*s);
        n[2][10]=aj1*(gy*s+(1.-gy)*s*s-s*s*s)*dl;
// fy
        n[4][2] = aj1*(6.*s-6.*s*s)/dl;
        n[4][4] = aj1*(1.+2.*gy-2.*(2.+gy)*s+3.*s*s);
        n[4][8] =-n[4][2];
        n[4][10]= aj1*(-2.*(1.-gy)*s+3.*s*s);

//  v
        aj1=1./(1.+2.*gz);
        n[1][1] = aj1*(1.+2.*gz-2.*gz*s-3.*s*s+2.*s*s*s);
        n[1][5] =-aj1*(-(1.+gz)*s+(2.+gz)*s*s-s*s*s)*dl;
        n[1][7] = aj1*(2.*gz*s+3.*s*s-2.*s*s*s);
        n[1][11]=-aj1*(gz*s+(1.-gz)*s*s-s*s*s)*dl;
// fz
        n[5][1] =-aj1*(6.*s-6.*s*s)/dl;
        n[5][5] = aj1*(1.+2.*gz-2.*(2.+gz)*s+3.*s*s);
        n[5][7] =-n[5][1];
        n[5][11]= aj1*(-2.*(1.-gz)*s+3.*s*s);
// fx
        n[3][3] = 1.-s;
        n[3][9] = s;
  
}


void soilbeam::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}

void soilbeam::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  long i,ipp,i1,transf, imat;
  double a2,a4,e,g,a,ixyz[3],beta[2],dl,ll,gy,gz,g2;
  ivector nodes(nne);
  vector vec(3),x(nne),y(nne),z(nne);
  matrix c(tncomp,tncomp),tran(3,3);
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  
  fillm (0.0,sm);
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (c,ipp);
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
//  Mc->give_widht (eid,&bPod);
  e=Mm->eliso[imat].e;  g=e/(1.+2.*Mm->eliso[imat].nu);
//  Mm->elmatstiff (d,ipp);
//  e=d[0][0];  g=d[1][1];
  c1[0]=c[0][0];
  c1[1]=c[1][1];
  c1[2]=c[2][2];
  c2[0]=c[3][3];
  c2[1]=c[4][4];
  c2[2]=c[5][5];
  ll=dl*dl;
//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;

// axial forces and torsion
  sm[0][0]= c1[0]*dl/3;
  sm[0][6]= c1[0]*dl/6; 
  g2=(c2[1]+c2[2])/2.;
  a2=bPod*bPod/2*( (c1[1]+c1[2])/2.*bPod/6.+sqrt((c1[1]+c1[2])/2.*g2) );
  a4=bPod*bPod/2*( (c1[1]+c1[2])/2.*bPod/6.+sqrt(g2*g2*g2/((c1[1]+c1[2])/2.)) );
  sm[3][3]=  a2/3.*dl+a4/dl;
  sm[3][9]=  a2/6.*dl-a4/dl;

//  direct s Iy
  g2=gy*gy;
  a2=bPod*(c1[2]+sqrt(c1[2]*c2[2])*2./bPod) *dl/(1.+2.*gy)/(1.+2.*gy);
  a4=bPod*(c2[2]+sqrt(c2[2]*c2[2]*c2[2]/c1[2])/bPod);
//  a4=0.0;
  sm[2][2] =  a2*(4.*g2/3.+1.4*gy+13./35.)+a4*4.*(g2+gy+.3);
  sm[2][4] =-(a2*(g2/6.+11.*gy/60.+11./210.)+ a4/10. )*dl;
  sm[2][8] =  a2*(2./3.*g2+0.6*gy+9./70.)-a4*4*(g2+gy+.3);
  sm[2][10]= (a2*(g2+0.9*gy+13./70.)/6.-a4/10. )*dl;
  sm[4][4] = (a2*(g2/30.+gy/30.+1./105.)+a4*(g2+gy+0.4)/3.)*ll;
  sm[4][8] =-(a2*(g2+0.9*gy+13./70.)/6.-a4/10. )*dl;
  sm[4][10]=-(a2*(g2/30.+gy/30.+1./140.)+a4*(g2+gy+0.1)/3.)*ll;
//  direct s Iz
  g2=gz*gz;
  a2=bPod*(c1[1]+sqrt(c1[1]*c2[1])*2./bPod)*dl/(1.+2.*gz)/(1.+2.*gz);
  a4=bPod*(c2[1]+sqrt(c2[1]*c2[1]*c2[1]/c1[1])/bPod);
//  a4=0.0;
  sm[1][1] =  a2*(4.*g2/3.+1.4*gz+13./35.)+a4*4.*(g2+gz+.3);
  sm[1][5] = (a2*(g2/6.+11.*gz/60.+11./210.)+a4/10. )*dl;
  sm[1][7] =  a2*(2./3.*g2+0.6*gz+9./70.)-a4*4.*(g2+gz+.3);
  sm[1][11]=-(a2*(g2+0.9*gz+13./70.)/6.-a4/10. )*dl;
  sm[5][5] = (a2*(g2/30.+gz/30.+1./105.)+a4*(g2+gz+0.4)/3.)*ll;
  sm[5][11]=-(a2*(g2/30.+gz/30+1./140.)+a4*(g2+gz+0.1)/3.)*ll;
  sm[5][7] = (a2*(g2+0.9*gz+13./70.)/6.-a4/10. )*dl;

  for (i=0;i<6;i++){
	  i1=i+6;
      sm[i1][i1]=sm[i][i];
      sm[i1][i]=sm[i][i1];
  }

  sm[4][2] = sm[2][4];
  sm[5][1] = sm[1][5];
  sm[7][5] = sm[5][7];
  sm[8][4] = sm[4][8];
  sm[11][1]= sm[1][11];
  sm[10][2]= sm[2][10];
  sm[7][11]=-sm[1][5];
  sm[11][7]=-sm[1][5];
  sm[8][10]=-sm[2][4];
  sm[10][8]=-sm[2][4];
  
  //  transformation of stiffness matrix to the global system
  lgmatrixtransfblock (sm,tran);
  //  transformation of stiffness matrix to the nodal system
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
//  fprintf (Out,"\n\n MT prvek soil cislo %ld",eid);
//  for (i=0;i<ndofe;i++){
//   fprintf (Out,"\n %4ld",i);
//   for (int j=0;j<ndofe;j++){
//    fprintf (Out," %15.5e",sm[i][i]);
//   }
//  }
//  fprintf (Out,"\n\n MT v MT Soilbeam cislo %ld \n",eid);
//  for (i=0;i<3;i++){fprintf (Out," %15.5e",c1[i]);}
//  for (i=0;i<3;i++){fprintf (Out," %15.5e",c2[i]);}
}



void soilbeam::strains (long lcid,long eid,long ri,long ci)
{
  long ipp,i,j,ij,transf;
  double dl;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),rl(ndofe),rg(ndofe),eps(6);
  matrix tran(3,3);
  
// r  .... end forces last time
// m  .... dl(i)/dl(i-1)
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  eldispl (lcid,eid,rl.a);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  
  //  transformation of displacement vector
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ndofe,ndofe);
    transf_matrix (nodes,tmat);
    //locglobtransf (rg,rl,tmat);
    lgvectortransf (rg,rl,tmat);
  }
  else{
    copyv (rl,rg);
  }
// *** transf. from gcs to lcs
//  globloctransf (rg,rl,tran);   

  ij=0;
  for (i=0;i<intordsm[0][0];i++){
    for (j=0;j<6;j++){  
      eps[j]=rl[ij];
	  ij++;
	}
//  fprintf (Out,"\n\n kontrola deformaci v soilbeam zapis,%ld \n",ipp);
//  for (long j=0;j<eps.n;j++){fprintf (Out," %15.5e",eps[j]);}
    Mm->storestrain(lcid,ipp,eps);
	ipp++;
  }
}
	  
void soilbeam::res_internal_forces (long lcid,long eid,vector &ifor)
{
   internal_forces (lcid, eid, 0, 0, ifor);
}

/**
   function computes internal forces
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param ifor - vector of internal forces
   
   20.12.2002
*/

void soilbeam::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long ipp,imat,i,j,k;
  double dl,ll,a2,a4,ixyz[3],beta[2],a,e,g,gy,gz,g2;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),rl(ndofe),rg(ndofe),f(ndofe),eps(6);
  matrix d(tncomp,tncomp),c(tncomp,tncomp),tran(3,3),r(2,ndofe);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  Mm->matstiff (c,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  e=Mm->eliso[imat].e;  g=e/(1.+2.*Mm->eliso[imat].nu);
//  Mm->elmatstiff (d,ipp);
//  e=d[0][0];  g=d[1][1]; 
  c1[0]=c[0][0];
  c1[1]=c[1][1];
  c1[2]=c[2][2];
  c2[0]=c[3][3];
  c2[1]=c[4][4];
  c2[2]=c[5][5];

  k=0;
  for (i=0;i<intordsm[0][0];i++){  
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
    Mm->givestress (lcid,ipp,eps);
    for (j=0;j<6;j++){  rl[k]=eps[j]; k++;}
	ipp++;
  }
  
  ll=dl*dl;
  f[0]= c1[0]*dl/3*rl[0]+ c1[0]*dl/6*rl[6]; 
  f[6]= c1[0]*dl/6*rl[0]+ c1[0]*dl/3*rl[6]; 
  g2=(c2[1]+c2[2])/2.;
  a2=bPod*bPod/2*( (c1[1]+c1[2])/2.*bPod/6.+sqrt((c1[1]+c1[2])/2.*g2) );
  a4=bPod*bPod/2*( (c1[1]+c1[2])/2.*bPod/6.+sqrt(g2*g2*g2/((c1[1]+c1[2])/2.)) );
  f[3]=  (a2/3.*dl+a4/dl)*rl[3]+(a2/6.*dl-a4/dl)*rl(9);
  f[9]=  (a2/6.*dl-a4/dl)*rl[3]+(a2/3.*dl+a4/dl)*rl(9);
  
//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;

//  direct s Iy
  g2=gy*gy;
  a2=bPod*(c1[2]+sqrt(c1[2]*c2[2])*2./bPod)*dl/(1.+2.*gy)/(1.+2.*gy);
  a4=bPod*(c2[2]+sqrt(c2[2]*c2[2]*c2[2]/c1[2])/bPod);
//  a4=0.0;
  f[2] =   (a2*(4.*g2/3.+1.4*gy+13./35.)+a4*4*(g2+gy+0.3))*rl[2]
	      +(a2*(2./3.*g2+0.6*gy+9./70.)-a4*4*(g2+gy+0.3))*rl[8]
	      -(a2*(g2/6.+11.*gy/60.+11./210.)+a4/10.)*dl*rl[4]
	      +(a2*(g2+0.9*gy+13./70.)/6.-a4/10. )*dl*rl[10];

//  f[4] = -a2*( (g2/6.+11.*gy/60.+11./210.)*rl[2]+(g2+0.9*gy+13./70.)/6.*rl[8] )*dl
//	     -a4/10.*(rl[4]-rl[10])*dl;
  f[4] = -(a2*(g2/6.+11.*gy/60.+11./210.)+ a4/10. )*dl*rl[2]
         +(a2*(g2/30.+gy/30.+1./105.)+a4*(g2+gy+0.4)/3.)*ll*rl[4]
         -(a2*(g2+0.9*gy+13./70.)/6.-a4/10. )*dl*rl[8]
         -(a2*(g2/30.+gy/30.+1./140.)+a4*(g2+gy+0.1)/3.)*ll*rl[10];
  f[8] =   (a2*(2./3.*g2+0.6*gy+9./70.)-a4*4*(g2+gy+0.3))*rl[2]
	      +(a2*(4.*g2/3.+1.4*gy+13./35.)+a4*4*(g2+gy+0.3))*rl[8]
	      -(a2*(g2+0.9*gy+13./70.)/6.-a4/10.)*dl*rl[4]
		  +(a2*(g2/6.+11.*gy/60.+11./210.)+ a4/10. )*dl*rl[10];

//  f[10]=  a2*( (g2+0.9*gy+13./70.)/6.*rl[2]+(g2/6.+11.*gy/60.+11./210.)*rl[8] )*dl
//	     +a4/10.*(rl[4]-rl[10])*dl;
  f[10]= +(a2*(g2+0.9*gy+13./70.)/6.-a4/10. )*dl*rl[2]
         -(a2*(g2/30.+gy/30.+1./140.)+a4*(g2+gy+0.1)/3.)*ll*rl[4]
		 +(a2*(g2/6.+11.*gy/60.+11./210.)+ a4/10. )*dl*rl[8]       
         +(a2*(g2/30.+gy/30.+1./105.)+a4*(g2+gy+0.4)/3.)*ll*rl[10];

//  direct s Iz
  g2=gz*gz;
  a2=bPod*(c1[1]+sqrt(c1[1]*c2[1])*2./bPod)*dl/(1.+2.*gz)/(1.+2.*gz);
  a4=bPod*(c2[1]+sqrt(c2[1]*c2[1]*c2[1]/c1[1])/bPod);
//  a4=0.0;
//  f[1] =  a2*( (4.*g2/3.+1.4*gz+13./35.)*rl[1]+(2./3.*g2+0.6*gz+9./70.)*rl[7] )
//	     +a4*4.*(g2+gz+.3)*(rl[5]-rl[11]);
  f[1] =   (a2*(4.*g2/3.+1.4*gz+13./35.)+a4*4*(g2+gz+0.3))*rl[1]
	      +(a2*(2./3.*g2+0.6*gz+9./70.)-a4*4*(g2+gz+0.3))*rl[7]
	      -(a2*(g2/6.+11.*gz/60.+11./210.)+a4/10.)*dl*rl[5]
	      +(a2*(g2+0.9*gz+13./70.)/6.-a4/10. )*dl*rl[11];
//  f[5] =  a2*( (g2/6.+11.*gz/60.+11./210.)*rl[1]+(g2+0.9*gz+13./70.)/6.*rl[7] )*dl
//	     +a4/10.*(rl[5]-rl[11])*dl;
  f[5] = -(a2*(g2/6.+11.*gz/60.+11./210.)+ a4/10. )*dl*rl[1]
         +(a2*(g2/30.+gz/30.+1./105.)+a4*(g2+gz+0.4)/3.)*ll*rl[5]
         -(a2*(g2+0.9*gz+13./70.)/6.-a4/10. )*dl*rl[7]
         -(a2*(g2/30.+gz/30.+1./140.)+a4*(g2+gz+0.1)/3.)*ll*rl[11];
//  f[7] =  a2*( (2./3.*g2+0.6*gz+9./70.)*rl[1]+(4.*g2/3.+1.4*gz+13./35.)*rl[7] )
//	     -a4*4.*(g2+gz+.3)*(rl[5]-rl[11]);
  f[7] =   (a2*(2./3.*g2+0.6*gz+9./70.)-a4*4*(g2+gz+0.3))*rl[1]
	      +(a2*(4.*g2/3.+1.4*gz+13./35.)+a4*4*(g2+gz+0.3))*rl[7]
	      -(a2*(g2+0.9*gz+13./70.)/6.-a4/10.)*dl*rl[5]
		  +(a2*(g2/6.+11.*gz/60.+11./210.)+ a4/10. )*dl*rl[11];

//		  -(a2*(g2/6.+11.*gz/60.+11./210.)+a4/10.)*(rl[5]-rl[11]);
//  f[11]= -a2*( (g2+0.9*gz+13./70.)/6.*rl[1]+(g2/6.+11.*gz/60.+11./210.)*rl[7] )*dl
//	     -a4/10.*(rl[5]-rl[11])*dl;
  f[11]= +(a2*(g2+0.9*gz+13./70.)/6.-a4/10. )*dl*rl[1]
         -(a2*(g2/30.+gz/30.+1./140.)+a4*(g2+gz+0.1)/3.)*ll*rl[5]
		 +(a2*(g2/6.+11.*gz/60.+11./210.)+ a4/10. )*dl*rl[7]       
         +(a2*(g2/30.+gz/30.+1./105.)+a4*(g2+gz+0.4)/3.)*ll*rl[11];

// do GCS z LCS  glvectortransfblock (f, tran);
  for (k=0;k<ndofe;k++){
    ifor[k]+=f[k];
  }
//  fprintf (Out,"\n\n R prvek Soilbeam cislo %ld \n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",rl[k]);}
//  fprintf (Out,"\n\n Fint prvek Soilbeam cislo %ld \n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",f[k]);}
//  fprintf (Out,"\n\n MT v Soilbeam cislo %ld \n",eid);
//  for (k=0;k<3;k++){fprintf (Out," %15.5e",c1[k]);}
//  for (k=0;k<3;k++){fprintf (Out," %15.5e",c2[k]);}
}
void soilbeam::internal_forces1 (long lcid,long eid,long ri,long ci,vector &/*ifor*/)
{
  long i,j,k,integr=2,ipp;
  double dl,ll,a,ixyz[3],beta[2],g2,gy,gz,a2x,a2y,a2z,e,g,s;
  ivector cn(ndofe),nodes(nne);
  vector x(nne),y(nne),z(nne),vec(3),rl(ndofe),rg(ndofe),f(ndofe),fx(tncomp),w(integr),gp(integr),eps(6);
  matrix n(tncomp,ndofe),c(tncomp,tncomp),d(tncomp,tncomp),tran(3,3),r(tncomp,2);
  
  Mt->give_elemnodes (eid,nodes);
  Mt->give_node_coord3d (x,y,z,eid);
  
  ipp=Mt->elements[eid].ipp[ri][ci];
  Mm->matstiff (c,ipp);
  Mm->elmatstiff (d,ipp);
  Mc->give_area (eid,a);
  Mc->give_mominer (eid,ixyz);
  Mc->give_shearcoeff (eid,beta);
  Mc->give_vectorlcs (eid, vec);
  beam_transf_matrix (tran,dl, vec,x,y,z,eid);
  e=d[0][0];  g=d[1][1]; ll=dl*dl;
  c1[0]=c[0][0];
  c1[1]=c[1][1];
  c1[2]=c[2][2];
  c2[0]=c[3][3];
  c2[1]=c[4][4];
  c2[2]=c[5][5];

  k=0;
  for (i=0;i<intordsm[0][0];i++){  
    Mm->computenlstresses (ipp,Mm->ip[ipp]);
    Mm->givestress (lcid,ipp,eps);
    for (j=0;j<6;j++){      rl[k]=eps[j];	 k++; }
	ipp++;
  }
  
//  direct s Iy
  if(beta[0]<Mp->zero) gy=0.0;
  else                 gy=6.0*e*ixyz[1]/beta[0]/g/a/ll;
//  direct s Iz
  if(beta[1]<Mp->zero) gz=0.;
  else                 gz=6.0*e*ixyz[2]/beta[1]/g/a/ll;
  g2=(c2[1]+c2[2])/2.;
  a2x=bPod*bPod/2*( (c1[1]+c1[2])/2.*bPod/6.+sqrt((c1[1]+c1[2])/2.*g2) );
  g2=gy*gy;
  a2y=bPod*(c1[2]+sqrt(c1[2]*c2[2])*2./bPod)*dl/(1.+2.*gy)/(1.+2.*gy);
  g2=gz*gz;
  a2z=bPod*(c1[1]+sqrt(c1[1]*c2[1])*2./bPod)*dl/(1.+2.*gz)/(1.+2.*gz);
//	 s=x/dl
     s=0.; 
  for (i=0;i<2;i++){
     geom_matrix (n,s, dl, gy, gz);
     fx[0]=( n[0][0]*rl[0]+n[0][6]*rl[6] )*c1[0];
     fx[3]=( n[3][3]*rl[3]+n[3][9]*rl[9] )*a2x;
     fx[1]=( n[1][1]*rl[1]+n[1][5]*rl[5]+n[1][7]*rl[7]+n[1][11]*rl[11] )*a2y;
     fx[5]=( n[5][1]*rl[1]+n[5][5]*rl[5]+n[5][7]*rl[7]+n[5][11]*rl[11] )*a2y;
     fx[2]=( n[2][2]*rl[2]+n[2][4]*rl[4]+n[2][8]*rl[8]+n[2][10]*rl[10] )*a2z;
     fx[4]=( n[4][2]*rl[2]+n[4][4]*rl[4]+n[4][8]*rl[8]+n[4][10]*rl[10] )*a2z;
	 s=dl;
  }

//  fprintf (Out,"\n\n R prvek beam cislo %ld\n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",rl[k]); }
//  fprintf (Out,"\n\n Fint prvek cislo %ld\n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",f[k]);}
}
