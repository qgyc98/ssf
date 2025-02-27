#include "soilplatetr.h"
#include "dkt.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "node.h"
#include "element.h"
#include <stdlib.h>
#include <math.h>

soilplatetr::soilplatetr (void)
{
  long i,j;

  nne=3;  ndofe=9;  tncomp=5;  napfun=3; ned=3;  nned=2;  
  intordmm=3; intordb=2; ssst=plates;
  nb=1; 
  
  ncomp = new long [nb];
  ncomp[0]=3;
  
  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }
  nip[0][0]=7;
  
  tnip=0;
  for (i=0;i<nb;i++){
    for (j=0;j<nb;j++){
      tnip+=nip[i][j];
    }
  }
  
  intordsm[0][0]=7;
}

soilplatetr::~soilplatetr (void)
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


double soilplatetr::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector areacoord(3);

  //  conversion of natural coordinates to area coordinates
  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  
  scprd (areacoord,nodval,f);
  
  return f;
}

/**
   function extracts components of the soil
   part of stiffness %matrix of the material to the matrix ds

   20.3.2002
*/
void soilplatetr::dbmat (matrix &db,double c1, double c2)
{

  db[0][0] = c1;   db[0][1] = 0.0;  db[0][2] = 0.0;
  db[1][0] = 0.0;  db[1][1] = c2;   db[1][2] = 0.0;
  db[2][0] = 0.0;  db[2][1] = 0.0;  db[2][2] = c2;
}


void soilplatetr::geom_matrix (matrix &gm,matrix &cn,vector &ax,vector &ay,vector &dl, vector &l)
{
  int i,i1,i2,i4;
  double lll;
//  vector ll(3);

//  ll[0]=l(0)*l(1)/2.;
//  ll[1]=l(1)*l(2)/2.;
//  ll[2]=l(2)*l(0)/2.;
  lll=l(0)*l(1)*l(2)/2.;
  for (i=0;i<nne;i++){
	 i1=i+1;
	 if (i1>=nne) i1=i1-nne;
     i2=i1+1;     
	 if (i2>=nne) i2=i2-nne;
     i4=i*3;     
     gm[0][i4]  =  l[i];
     gm[0][i4+1]=(l[i]*l[i]*l[i1]+lll)*dl[i]*cn[i][0]-(l[i]*l[i]*l[i2]+lll)*dl[i2]*cn[i2][0];
     gm[0][i4+2]=(l[i]*l[i]*l[i1]+lll)*dl[i]*cn[i][1]-(l[i]*l[i]*l[i2]+lll)*dl[i2]*cn[i2][1];
     gm[1][i4]  = ay[i];
//     gm[1][i4+1]=-ay[i] *( (2.*l[i]*l[i2]+ll[i1])*dl[i2]*cn[i2][1]+(2.*l[i]*l[i1]+ll[i1])*dl[i]*cn[i][1] )
//                 -ay[i1]*( (ll[i2])*dl[i2]*cn[i2][1]+(l[i]*l[1]+ll[i2])*dl[i]*cn[i][1] )
//                 -ay[i2]*( (l[i]*l[i]+ll[i])*dl[i2]*cn[i2][1]+(ll[i])*dl[i]*cn[i][1] )

//     gm[1][i4+2]=-ay[i] *( (2.*l[i]*l[i2]+ll[i1])*dl[i2]*cn[i2][0]+(2.*l[i]*l[i1]+ll[i1])*dl[i]*cn[i][0] )
//                 -ay[i1]*( (ll[i2])*dl[i2]*cn[i2][0]+(l[i]*l[1]+ll[i2])*dl[i]*cn[i][0] )
//                 -ay[i2]*( (l[i]*l[i]+ll[i])*dl[i2]*cn[i2][0]+(ll[i])*dl[i]*cn[i][0] )

     gm[2][i4]  = ax[i];
//     gm[2][i4+1]=-ax[i] *( (2.*l[i]*l[i2]+ll[i1])*dl[i2]*cn[i2][1]+(2.*l[i]*l[i1]+ll[i1])*dl[i]*cn[i][1] )
//                 -ax[i1]*( (ll[i2])*dl[i2]*cn[i2][1]+(l[i]*l[1]+ll[i2])*dl[i]*cn[i][1] )
//                 -ax[i2]*( (l[i]*l[i]+ll[i])*dl[i2]*cn[i2][1]+(ll[i])*dl[i]*cn[i][1] )

//     gm[2][i4+2]=-ax[i] *( (2.*l[i]*l[i2]+ll[i1])*dl[i2]*cn[i2][0]+(2.*l[i]*l[i1]+ll[i1])*dl[i]*cn[i][0] )
//                 -ax[i1]*( (ll[i2])*dl[i2]*cn[i2][0]+(l[i]*l[1]+ll[i2])*dl[i]*cn[i][0] )
//                 -ax[i2]*( (l[i]*l[i]+ll[i])*dl[i2]*cn[i2][0]+(ll[i])*dl[i]*cn[i][0] )

  }
}

void soilplatetr::transf_matrix (ivector &nodes,matrix &tmat)
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
void soilplatetr::stiffness_matrix (long eid,long ri, long ci, matrix &sm,vector &x, vector &y)
{
  long i,ipp;
  double det,jac,thick;
  ivector nodes(nne);
  vector w,gp1,gp2,c1(nne),c2(nne),t(nne),l(3),dl(3),ax(3),ay(3);
  matrix gm,cn(3,2),cc(3,3);
//     jakobian
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
//  Mc->give_soil (eid,nodes,c1,c2);
//  c1[0]=1000; c1[1]=1000;c1[2]=1000;
//  c2[0]=1000;c2[1]=1000;c2[2]=1000;

// triangular  
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
  dl[0] = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1] = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
  // matrix vector of normal of sides {cos,sin} 
  cn[0][0]= (y[1]-y[0])/dl[0];
  cn[0][1]=-(x[1]-x[0])/dl[0];
  cn[1][1]=-(x[2]-x[1])/dl[1];
  cn[1][0]= (y[2]-y[1])/dl[1];
  cn[2][1]=-(x[0]-x[2])/dl[2];
  cn[2][0]= (y[0]-y[2])/dl[2];
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/det;  ax[1]=(x[0]-x[2])/det;  ax[2]=(x[1]-x[0])/det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/det;  ay[1]=(y[2]-y[0])/det;  ay[2]=(y[0]-y[1])/det;

  fillm (0.0,sm);
  
//  ii=Mt->elements[eid].ipp[sip][0];
  reallocv (intordsm[0][0],w);  reallocv (intordsm[0][0],gp1);  reallocv (intordsm[0][0],gp2);
  reallocm (ncomp[0],ndofe,gm); reallocm (ncomp[0],ncomp[0],cc);
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);
  ipp=Mt->elements[eid].ipp[ri][ci];
  for (i=0;i<intordsm[0][0];i++){
    thick = approx_nat (gp1[i],gp2[i],t);
    l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
    jac=det*w[i];
    geom_matrix (gm,cn,ax,ay,dl,l);
    Mm->matstiff (cc,ipp);
    //     c11 = c1[0]*l[0]+c1[1]*l[1]+c1[2]*l[2];
    //     c22 = c2[0]*l[0]+c2[1]*l[1]+c2[2]*l[2];
    //    dbmat (cc,c11,c22);
    bdbjac (sm,gm,cc,gm,jac*thick);
    ipp++;  
  }
/*
  long ij; 
  fprintf (Out,"\n");
  fprintf (Out,"\n\n prvek cislo %ld \n",eid);
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[0][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[1][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[2][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[3][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[4][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[5][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[6][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[7][ij]);}fprintf (Out,"\n");
  for (ij=0;ij<ndofe;ij++){fprintf (Out,"%12.3e",sm[8][ij]);}fprintf (Out,"\n");*/
}

void soilplatetr::res_stiffness_matrix (long eid,matrix &sm)
{
  vector x(nne),y(nne);
  Mt->give_node_coord2d (x,y,eid);
  stiffness_matrix (eid,0,0,sm,x,y);
}

void soilplatetr::res_mainip_strains (long lcid,long eid)
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

  mainip_strains (lcid,eid,0,0,x,y,r);
}


void soilplatetr::mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long i,ipp;
  double det;
  vector ax(3),ay(3),dl(3),l(3),gp1,gp2,w,eps;
  matrix cn(3,2),gm;

  // translation to center it is necesary for w which in local coord. and no in relation
// triangular  
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
  dl[0] = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1] = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
  // matrix vector of normal of sides {cos,sin} 
  cn[0][0]= (y[1]-y[0])/dl[0];
  cn[0][1]=-(x[1]-x[0])/dl[0];
  cn[1][1]=-(x[2]-x[1])/dl[1];
  cn[1][0]= (y[2]-y[1])/dl[1];
  cn[2][1]=-(x[0]-x[2])/dl[2];
  cn[2][0]= (y[0]-y[2])/dl[2];
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/2./det;  ax[1]=(x[0]-x[2])/2./det;  ax[2]=(x[1]-x[0])/2./det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/2./det;  ay[1]=(y[2]-y[0])/2./det;  ay[2]=(y[0]-y[1])/2./det;
  
//  fprintf (Out,"\n\n deformace");
//  fprintf (Out,"\n\n prvek cislo %ld",eid);
  ipp=Mt->elements[eid].ipp[ri][ci];
  reallocv (intordsm[0][0],w);  reallocv (intordsm[0][0],gp1);  reallocv (intordsm[0][0],gp2);
  reallocm (ncomp[0],ndofe,gm);
  reallocv (ncomp[0],eps);
  gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[0][0]);

    for (i=0;i<intordsm[0][0];i++){
     l[0]=gp1[i];  l[1]=gp2[i];  l[2]=1.0-l[0]-l[1];
     geom_matrix (gm,cn,ax,ay,dl,l);
     mxv (gm,r,eps);
     Mm->storestrain (lcid,ipp,cncomp[0],ncomp[0],eps);
     ipp++;
//		fprintf (Out,"\n");
//		for (long ij=0;ij<ncomp[0];ij++){fprintf (Out,"%20.10e",eps[ij]);}
	}
}

void soilplatetr::res_allip_stresses (long lcid,long eid)
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
void soilplatetr::allip_stresses (long lcid,long eid,long ri,long ci)
{
  long ipp,i;
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
      Mm->matstiff (cc,ipp);
	  Mm->givestrain (lcid,ipp,cncomp[0],ncomp[0],eps);
	  mxv (cc,eps,sig);	
	  Mm->storestress (lcid,ipp,cncomp[0],sig);      
//  for (k=0;k<3;k++){
//    printf (Out," %15.5e",eps[k]); fprintf (Out," %15.5e",sig[k]);  }
//    fprintf (Out,"\n");
	  ipp++;  
  }

//  fprintf (Out,"\n\n Fint soil prvek cislo %ld\n",eid);
//  for (k=0;k<ndofe;k++){fprintf (Out," %15.5e",ifor[k]);}
  
}


/**
   function computes correct stresses at integration points on element

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 27.9.2008
*/
void soilplatetr::compute_nlstress (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++)
      {
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
   
   JK, 23.9.2008
*/
void soilplatetr::compute_nlstressincr (long /*lcid*/,long eid,long ri,long ci)
{
  long i,ii,jj,ipp;
  
  for (ii=0;ii<nb;ii++){
    for (jj=0;jj<nb;jj++){
      if (intordsm[ii][jj]==0)  continue;
      
      ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
      
      for (i=0;i<intordsm[ii][jj];i++){
        //  computation of correct increments of stresses
        if (Mp->strcomp==1)
          Mm->computenlstressesincr (ipp);
        ipp++;
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
void soilplatetr::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
{
  integratedquant iq;
  iq=locstress;
  
  //  computation of stresses
  compute_nlstress (lcid,eid,ri,ci);
  
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
void soilplatetr::res_internal_forces (long lcid,long eid,vector &ifor)
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
    //globloctransf (ifor,v,tmat);
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
void soilplatetr::incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y)
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
void soilplatetr::res_incr_internal_forces (long lcid,long eid,vector &ifor)
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
void soilplatetr::elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y)
{
  long i,ii,ipp;
  double thick,det;
  ivector nodes(nne);
  vector w,gp1,gp2,t(nne),ipv,contr(ndofe),l(3),dl(3),ax(3),ay(3);
  matrix gm,cn(3,2);
  
  Mt->give_elemnodes (eid,nodes);
  Mc->give_thickness (eid,nodes,t);
  
  fillv (0.0,nv);

  // triangular  
  det = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
  dl[0] = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  dl[1] = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
  dl[2] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
  // matrix vector of normal of sides {cos,sin} 
  cn[0][0]= (y[1]-y[0])/dl[0];
  cn[0][1]=-(x[1]-x[0])/dl[0];
  cn[1][1]=-(x[2]-x[1])/dl[1];
  cn[1][0]= (y[2]-y[1])/dl[1];
  cn[2][1]=-(x[0]-x[2])/dl[2];
  cn[2][0]= (y[0]-y[2])/dl[2];
  //   der   dL(1)/dy
  ax[0]=(x[2]-x[1])/det;  ax[1]=(x[0]-x[2])/det;  ax[2]=(x[1]-x[0])/det;
  //   der   dL(1)/dx
  ay[0]=(y[1]-y[2])/det;  ay[1]=(y[2]-y[0])/det;  ay[2]=(y[0]-y[1])/det;
  

  
  for (ii=0;ii<nb;ii++){
    reallocv (intordsm[ii][ii],gp1);
    reallocv (intordsm[ii][ii],gp2);
    reallocv (intordsm[ii][ii],w);
    reallocm (ncomp[ii],ndofe,gm);
    reallocv (ncomp[ii],ipv);
    
    gauss_points_tr (gp1.a,gp2.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ri+ii][ci+ii];
    
    for (i=0;i<intordsm[ii][ii];i++){
      thick = approx_nat (gp1[i],gp2[i],t);
      //  function assembles required quantity at integration point
      Mm->givequantity (iq,lcid,ipp,cncomp[ii],ipv);
	
      //  strain-displacement (geometric) matrix
      l[0]=gp1[i];
      l[1]=gp2[i];
      l[2]=1.0-l[0]-l[1];
      geom_matrix (gm,cn,ax,ay,dl,l);
      
      //  contribution to the nodal values
      mtxv (gm,ipv,contr);
      
      cmulv (det*w[i]*thick,contr);
      
      //  summation
      addv(contr,nv,nv);
      
      ipp++;
    }
  }
}
