#include "consol.h"
#include "matrix.h"
#include "vector.h"
#include "math.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "globmat.h"
#include "intpoints.h"
#include "winpast.h"
#include <stdlib.h>

consol::consol (void)
{
  nRetTime=7;  
}

consol::~consol (void)
{
  
}

void consol::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&cn,&gama,&vv,&ni,&e0,&m,&hpodl,&hpodlmax,&rr,&nf);
  
}

long consol::numberOfConsol ()
{
	return nRetTime;
}
/**
   function computes stress increment due to viscous strains
   
   @param ipp - integration point pointer
   
   7.2008
*/
void consol::nlstressesincr (long ipp)
{
  nc=Mm->ip[ipp].ncompstr;
  if(nc==6) ncc=2;
  if(nc==3) ncc=3;
 
  ddTime=Mp->timecon.actualforwtimeincr ();
  ccTime=Mp->timecon.actualtime ();
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  phase1(ipp);

}



/**
   function returns eps from history
   
   @param screep - vector of history
   @param epsscr - vector deformation of history
   
   10.10.2002
*/
void consol::nlstresses (long ipp)
{
  nc=Mm->ip[ipp].ncompstr;
  if(nc==6) ncc=2;
  if(nc==3) ncc=3;
  
  ddTime=Mp->timecon.actualforwtimeincr ();
  ccTime=Mp->timecon.actualtime ();
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  phase2(ipp);
}



void consol::phase1 (long ipp)
{
  long i,ii,liche;
  double delY,dj;
  vector epscr(nc),deps(nc),sig(nc);
  matrix d(nc,nc),screep(ncc,nRetTime);
  
  //  right hand side from time dependent value computation
  //  creep function for shrink od casu pomt do cctime
  for (i=0; i<nc; i++){epscr[i]=0.0;}
  get_hc ( nc, ipp);
  delY=2.4674*ncn*cn*ddTime/hpodl/nh/hpodl/nh;
  // for beam are  sig,epscr,(N,Vy,Vz,Mx,My,Mz) in int.point
  for (i=0;i<nRetTime;i++){
    ii=nc+ncc*i;
    liche=2*i+1;
    dj=delY*liche*liche;
    if(nc==6){
      screep[0][i]=Mm->ip[ipp].other[ii];	
      epscr[2]=epscr[2]+(1.-exp(-dj))*screep[0][i];
      screep[1][i]=Mm->ip[ipp].other[ii+1];	
      epscr[4]=epscr[4]+(1.-exp(-dj))*screep[1][i];
    }
    // for plate  ( w, fix, fiy )
    else{
      screep[0][i]=Mm->ip[ipp].other[ii];	
      epscr[0]=epscr[0]+(1.-exp(-dj))*screep[0][i];
      screep[1][i]=Mm->ip[ipp].other[ii+1];	
      epscr[1]=epscr[1]+(1.-exp(-dj))*screep[1][i];
      screep[2][i]=Mm->ip[ipp].other[ii+2];	
      epscr[2]=epscr[2]+(1.-exp(-dj))*screep[2][i];
    }
    
  }
  
  //        fprintf (Out,"\n\n *** phase 1 v case t  prirustek dt ");
  //        fprintf (Out," %e\t%e", ccTime, ddTime);
  //        fprintf (Out,"\n\n cteni screep z other");
  //        fprintf (Out,"\n %le %le ",screep[0][0],screep[1][0]);
  //        fprintf (Out,"\n\n F A Z E  1 epscr od skrytych = zatizeni \n");
  //        fprintf (Out," %le %le \n ",epscr[2],epscr[4]);
  
  Mm->storestress (0,ipp,epscr);
  // in arr. stress is displac.
  
}

void consol::phase2 (long ipp)
{
  long i,j,ii,liche;
  double delY,dj;
  vector epscr(nc),deps(nc),sig(nc);
  matrix d(nc,nc),screep(ncc,nRetTime);
  
  //  new total strain and stress
  //        fprintf (Out,"\n\n **** phase 2 v case t   prirustek dt");
  //        fprintf (Out," %e\t%e", ccTime, ddTime);
  //        fprintf (Out,"\n\n strain new");
  //        fprintf (Out,"\n %le %le ",Mm->ip[ipp].strain[2],Mm->ip[ipp].strain[4]);
  //        fprintf (Out,"\n\n strain old");
  //        fprintf (Out,"\n %le %le ",Mm->ip[ipp].other[2],Mm->ip[ipp].other[4]);
  // arr. strain = is total strain (including e from screep) from t=0
  for (i=0; i<nc; i++){deps[i]=Mm->ip[ipp].strain[i] - Mm->ip[ipp].other[i];}   //increment of str
  get_hc ( nc,ipp);
  delY=2.4674*ncn*cn*ddTime/hpodl/nh/hpodl/nh;
  //  increment of e=e-es  for  stress components
  for (i=0;i<nRetTime;i++){
    ii=nc+ncc*i;
    liche=2*i+1;
    dj=delY*liche*liche;
    if(nc==6){
      screep[0][i]=Mm->ip[ipp].other[ii];
      deps[2]=deps[2]-(1.-exp(-dj))*screep[0][i];
      screep[1][i]=Mm->ip[ipp].other[ii+1];
      deps[4]=deps[4]-(1.-exp(-dj))*screep[1][i];
    }
    else{
      screep[0][i]=Mm->ip[ipp].other[ii];
      deps[0]=deps[0]-(1.-exp(-dj))*screep[0][i];
      screep[1][i]=Mm->ip[ipp].other[ii+1];
      deps[1]=deps[1]-(1.-exp(-dj))*screep[1][i];
      screep[2][i]=Mm->ip[ipp].other[ii+2];
      deps[2]=deps[2]-(1.-exp(-dj))*screep[2][i];
    }
  }
  // save total strain (from t=0)
  for (i=0; i<nc; i++){
    Mm->ip[ipp].other[i]=Mm->ip[ipp].strain[i];     // it is last strain in arr. strain
    sig[i]=deps[i];
  }
  //  time history
  seps_time (screep,sig,nc,ncc,ipp);
  for (i=0; i<nRetTime; i++){
    ii=nc+ncc*i;
    for (j=0; j<ncc; j++){
      Mm->ip[ipp].other[ii+j]=screep[j][i];
    }
  }
  //        fprintf (Out,"\n\n nh, nc1,nc2");
  //        fprintf (Out," %e\t%e\t%e", nh,nc1,nc2);
  //        fprintf (Out,"\n\n ukladani screep=other");
  //        fprintf (Out,"\n %le %le ",screep[0][0],screep[1][0]);
  
  // save total strain only from consolidation(from t=0) 	position	ii=nc+ncc*nRetTime;
  // for get_hc
  ii=nc+ncc*nRetTime;
  if(nc==6){
    Mm->ip[ipp].other[ii]+=deps[2];
  }
  else{
    Mm->ip[ipp].other[ii]+=deps[0];
  }
  
  
  //        fprintf (Out,"\n\n screep=other");
  //        fprintf (Out,"\n %le %le ",screep[0][0],screep[1][0]);
  //        fprintf (Out,"\n\n eps celkove,   eps celkove minule=other, screep=other, delta eps");
  //        fprintf (Out,"\n %le %le %le %le %le %le",epscr[0],epscr[1],epscr[2],epscr[3],epscr[4],epscr[5]);
  //        fprintf (Out,"\n %le %le ",Mm->ip[ipp].other[0],Mm->ip[ipp].other[1]);
  //        fprintf (Out,"\n %le %le ",Mm->ip[ipp].other[2],Mm->ip[ipp].other[3]);
  //        fprintf (Out,"\n %le %le %le %le %le %le",deps[0],deps[1],deps[2],deps[3],deps[4],deps[5]);
  //  Stress data storage first position, Hide (strain) seconde position
  
}

/**
   function apdates material parameters
   
   @param d - matrix of soil 

   10.10.2002
*/
void consol::matstiff (matrix &d,long ipp)
//  function returns elastic stiffness matrix of soil
{
  long iwink,j,liche,nc=Mm->ip[ipp].ncompstr;
  double delY,dj;  
  if(Mm->ip[ipp].gemid()==2)  //consol is jointed with soil on 2 pozition
    {
      //    iwink=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
      iwink=Mm->ip[ipp].idm[1];
      Mm->wpast[iwink].matstiff (d,ipp);
    }
  //  d[0][0]=1000.;d[1][1]=1000.;d[2][2]=1000.;
  //  d[3][3]=0.001;d[4][4]=0.001;d[5][5]=0.001;
  
  //ddTime=Mp->timefun.getval(Mp->time);
  ddTime=Mp->timecon.actualforwtimeincr ();
  vlivTCSum=0.0;
  get_hc ( d.n,ipp);
  delY=2.4674*ncn*cn*ddTime/hpodl/nh/hpodl/nh;
  for (j=0;j<nRetTime;j++){
    liche=2*j+1;
    dj=delY*liche*liche;
    vlivTCSum=vlivTCSum+0.81057*(1.-(1.-exp(-dj))/dj)/liche/liche;
  }
  //        fprintf (Out,"\n vlivt v MT %le \n ",vlivTCSum);
  //		fprintf (Out,"\n\n eh v case t   od t0 prirustek dt");
  //		fprintf (Out," %e\t%e\t%e\t%e",et, pomt, ccTime, ddTime);
  
  //  New material modulus in Time
  for (j=0;j<d.n;j++){
    d[j][j]=d[j][j]/vlivTCSum;
  }
  //  New material modulus in change h
  if(nc==6){
    for (j=0;j<d.n/2;j++){
      d[j][j]=d[j][j]*nc1;
    }
    for (j=d.n/2;j<d.n;j++){
      d[j][j]=d[j][j]*nc2;
    }
  }
  else{
    d[0][0]=d[0][0]*nc1;
    for (j=1;j<d.n;j++){
      d[j][j]=d[j][j]*nc2;
    }
  }
  /*
    fprintf (Out,"\n\n D cons cislo %ld",ipp);
    for (j=0;j<d.n;j++){
    fprintf (Out," %15.5e",d[j][j]);
    }
  */
}

/**
   function returns history
   
   @param screep - vector of history
   @param sig - vector od add stress
   
   10.10.2002
*/
void consol::seps_time (matrix &screep,vector &sig,long nc,long ncc,long ipp)
{
  int j,i,liche;
  double dj,delY;
  vector rp(ncc);
  

  
  ddTime=Mp->timecon.actualforwtimeincr ();
  vlivTCSum=0.0;
  get_hc ( nc,ipp);
  delY=2.4674*ncn*cn*ddTime/hpodl/nh/hpodl/nh;
  for (j=0;j<nRetTime;j++){
    liche=2*j+1;
    dj=delY*liche*liche;
    vlivTCSum=vlivTCSum+0.81057*(1.-(1.-exp(-dj))/dj)/liche/liche;
  }
  
  if(nc==6){
	rp[0]=sig[2]/vlivTCSum;
	rp[1]=sig[4]/vlivTCSum;
//        fprintf (Out,"\n  posuny od prir. zat. pro screep");
//        fprintf (Out,"\n %le %le ",sig[2],sig[4]);
  }
  else {
	rp[0]=sig[0]/vlivTCSum;
	rp[1]=sig[1]/vlivTCSum;
	rp[2]=sig[2]/vlivTCSum;
  }

  
  
  delY=2.4674*ncn*cn*ddTime/hpodl/nh/hpodl/nh;
  for (j=0;j<nRetTime;j++){
	liche=2*j+1;
    dj=delY*liche*liche;
	for (i=0;i<ncc;i++){
		screep[i][j]=exp(-dj)*screep[i][j]+0.81057*rp[i]*(1.-exp(-dj))/dj/liche/liche;
	}
  }

//        fprintf (Out,"\n\n screep=other");
//        fprintf (Out,"\n %le %le ",screep[0][0],screep[1][0]);
}


/**
   function returns change of value
   
   @param vv - value
   @param sig - vector of presure
   
   10.10.2002
*/
void consol::get_hc (long nc,long ipp)
{
  int ii, ncc;
  double f,smax,ee,p,pp,hactual,*sumf,*sizexyz;
  vector deps(nc);
  
//Elements without change of parameters of consolidation
  if(nf<1.e-6)     //Elements without load
  {
     nh=1.;  nc1=1.;  nc2=1.;  ncn=1.;
	 return;
  }
  
    sizexyz = new double [3];
    sumf = new double [nc];
    Gtm->give_domain_sizes (sizexyz);
    Mb->give_comp_sum (sumf);  
  
    if(nc==6){
     if(rr<1.e-6)
		{if(sizexyz[0]>sizexyz[1])
			{smax=sizexyz[0]/2.;}
		else
			{smax=sizexyz[1]/2.;}
		}
	 else
		{smax=rr;}

     f=nf*sumf[2]/smax/2.;       //Beam forces z
	 ncc=2;
	}
    else {
     if(rr<1.e-6)
		{smax=(sizexyz[0]+sizexyz[1])/4.;
		}
	 else
		{smax=rr;}

     f=nf*sumf[0]/smax/smax/3.14;       //Plate forces z
     ncc=3;
	}

    p=0.5*(1.-gama*vv/(-f));
    if(p>=0.5)
	 pp=0.0001;
    else if(p>=0.4)
	 pp=(0.5*0.4/(0.5-0.4)+0.5)-0.5/(0.5-0.4)*p;
    else if(p>=0.24)
	 pp=(0.5*0.24/(0.4-0.24)+1.)-0.5/(0.4-0.24)*p;
    else if(p>=0.06)
	 pp=(0.06/(0.24-0.06)+2.)-1./(0.24-0.06)*p;
    else if(p>=0.02)
	 pp=(0.02/(0.06-0.02)+3.)-1./(0.06-0.02)*p;
    else if(p>=0.005)
	 pp=(0.005/(0.02-0.005)+4.)-1./(0.02-0.005)*p;
    else
	 pp=5.;

    hactual=smax/pp*sqrt((2.-2.*ni)/(1.-2.*ni));
    if(hpodlmax>=0.0)
	  if(hactual>hpodlmax) 
		  hactual=hpodlmax;

    nh=hactual/hpodl;
//from p=0.1   c1/c10 =0.15,    1=1, 2,5=2,6
    nc1=(0.051/nh/nh+0.888/nh+0.0606);   
//from p=0.1   c2/c20 =4,    1=1, 2,5=0,55
    nc2=(1.2639/nh/nh-4.7236/nh+4.46);   

// epsilon of soil from consolidation
  if(e0<1.e-6)  
  {
     ncn=1.;
  }
  else
  {
    ii=nc+ncc*nRetTime;
    ee =	Mm->ip[ipp].other[ii]/nh/hpodl;
    ee = -ee*(1.+e0)+e0;
    ncn=pow((ee/e0),m);
  }
    delete [] sizexyz;
    delete [] sumf;

}

void consol::updateval (long ipp, long im, long ido)
{
	long i,n=Mm->givencompeqother(ipp,im);
	for (i=0;i<n;i++){
		Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
	}
}


void consol::give_dstresses_eqother(vector &dsigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    dsigma[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc)+i];
  }
}
