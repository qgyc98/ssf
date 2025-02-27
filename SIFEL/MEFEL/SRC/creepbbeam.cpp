#include "creepbbeam.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "globmat.h"
#include "intpoints.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vector.h"
#include <stdlib.h>
#include <math.h>

creepbbeam::creepbbeam (void)
{
  nRetTime=7;  
  allocv (nRetTime,retTime);
  allocv (nRetTime,ert);
  retTime[0]=1.0e-9;
  retTime[1]=1.0e-3;
  retTime[2]=1.0;
  retTime[3]=14.0;
  retTime[4]=100.0;
  retTime[5]=800.0;
  retTime[6]=3000.0;
  type_h=1;
  type_temp=1;
  timeMax=Mp->timecon.endtime ();
  h_s = 0.4;   // end humidity
  temp_s=20;    // temperature
  esht=0.0;
  nc=1;
  k_s = 1.0;    //shape factor
  r_s = 0.8;
  k_d = 0.12;  //effective cross area 2V/S
  timemat=0.0;
}

creepbbeam::~creepbbeam (void)
{
}

void creepbbeam::creepinit (long mie)
{
}

void creepbbeam::read (FILE *in)
{
  fscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %ld %ld",
	          &tb,&t_w,&fc,&wc,&sc,&gc,&c_s,&a1,&type_h,&type_temp);
  
  if (type_h==1){fscanf (in,"%lf ",&h_s); }
  if (type_temp==1){fscanf (in,"%lf ",&temp_s); }

}
/**
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

*/
double creepbbeam::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function inversion sym. matrix

   @param a - matrix

*/
void creepbbeam::inv_sym (matrix &a)
{
  long i,j,k;
  double diag,an;
    for (k=0;k<a.n;k++){
       diag = a[k][k];
       if (diag<=0.0) { abort();}
       for (i=0;i<a.n;i++){
         an = a[i][k]/diag;
         if (i!=k) { 
           for (j=i;j<a.n;j++){
             if(j!=k) {
               a[i][j] = a[i][j] - an*a[k][j];
               a[j][i] = a[i][j];
			 }
		   }
		 }
         a[i][k] = an;
         a[k][i] = an;
       }
       a[k][k] = - 1.0/diag;
    }

	for (i=0;i<a.n;i++){
       for (j=0;j<a.n;j++){
         a[i][j] = - a[i][j];
	   }
	}
}

/**
   function returns eps from history
   
   @param screep - vector of history
   @param epsscr - vector deformation of history
   
   10.10.2002
*/
void creepbbeam::nlstresses (long ipp)
{
  //long nc=Mm->ip[ipp].ncompstr;
  //timeMax=Mp->end_time;
  //ddTime=Mp->timefun.getval(Mp->time);
  ddTime=Mp->timecon.actualforwtimeincr ();
  //ccTime=Mp->time;
  ccTime=Mp->timecon.actualtime ();
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //if(type_h!=1) {h_s=Mm->givenonmechq(rel_hum, ipp);}   //actual moisture
  //if(type_temp!=1) {temp=Mm->givenonmechq(temperature, ipp);}  //actual temperature
  get_h (ipp);     //initial  moisture
  get_temp (ipp);     //initial  temperature
  if (Mp->phase==1){
	phase1(ipp);
  }

  if (Mp->phase==2){
	phase2(ipp);
  }

}

void creepbbeam::phase1 (long ipp)
{
//  right hand side from time dependent value computation
//  creep function for shrink from pomt to cctime
  long i,j,ii;
  double qq,dj;
  vector epscr(nc),deps(nc),sig(nc);
  matrix d(nc,nc),screep(nc,nRetTime);

    Mm->matstiff(d,ipp);
    e0=Mm->eliso[imat].e;
    mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
	for (i=0; i<nc; i++){epscr[i]=0.0;}
	qq=2.0/3.0;
	for (i=0;i<nRetTime;i++){
		ii=nc*(i+1);
		dj=pow((ccTime+ddTime)/retTime[i],qq)-pow(ccTime/retTime[i],qq);
	    for (j=0; j<nc; j++){
		 screep[j][i]=Mm->ip[ipp].other[ii+j];	
	     epscr[j]=epscr[j]+(1.-exp(-dj))*screep[j][i];
		}
	}
//        fprintf (Out,"\n\n v case t  prirustek dt faze 1");
//        fprintf (Out," %e\t%e", ccTime, ddTime);
//        fprintf (Out,"\n\n eps faze1");
//        fprintf (Out," %le %le %le %le %le %le",epscr[0],epscr[1],epscr[2],epscr[3],epscr[4],epscr[5]);
//        fprintf (Out,"\n\n tau faze1");
//        fprintf (Out," %le %le %le %le %le %le",screep[0][0],screep[1][0],screep[2][0],screep[3][0],screep[4][0],screep[5][0]);
//  shrink and temp
	epscr[0]=epscr[0]+deps0;
	epscr[1]=epscr[1]+deps0;
	if(epscr.n==6) epscr[2]=epscr[2]+deps0;
//		fprintf (Out,"\n\n time, esht prir faze 1");
  
//		fprintf (Out,"\n\n phase 1 eps1, eps2 prir");
//		fprintf (Out," %e %e",epscr[0], epscr[1]);
    Mm->stiff_deps_creep (ipp, epscr, nc,nc);  // save sig
}



void creepbbeam::phase2 (long ipp)
{
//  new total strain and stress
  long i,j,ii;
  double qq,dj;
  vector epscr(nc),deps(nc),sig(nc);
  matrix d(nc,nc),screep(nc,nRetTime);

		fprintf (Out,"\n\n phase 2 v case t   prirustek dt");
		fprintf (Out," %e\t%e", ccTime, ddTime);
// arr. strain = is total strain (including e from screep) from t=0
	for (i=0; i<nc; i++){deps[i]=Mm->ip[ipp].strain[i]-Mm->ip[ipp].other[i];} //increment deps
	qq=2.0/3.0;
	for (i=0;i<nRetTime;i++){
		ii=nc*(i+1);
		dj=pow((ccTime+ddTime)/retTime[i],qq)-pow(ccTime/retTime[i],qq);
	    for (j=0; j<nc; j++){
		 screep[j][i]=Mm->ip[ipp].other[ii+j];
		 deps[j]=deps[j]-(1.-exp(-dj))*screep[j][i];     // deps=deps-e(screep)
		}
	}
//  shrink and temp
	deps[0]=deps[0]-deps0;
	deps[1]=deps[1]-deps0;
	if(deps.n==6) deps[2]=deps[2]-deps0;
//  increment of stress components
//        fprintf (Out,"\n\n eps celkove,   eps celkove minule,  delta eps");
//        fprintf (Out,"\n %le %le %le %le %le %le",epscr[0],epscr[1],epscr[2],epscr[3],epscr[4],epscr[5]);
//        fprintf (Out,"\n %le %le %le ",Mm->ip[ipp].other[0],Mm->ip[ipp].other[1],Mm->ip[ipp].other[2]);
//        fprintf (Out,"\n %le %le %le ",Mm->ip[ipp].other[3],Mm->ip[ipp].other[4],Mm->ip[ipp].other[5]);
        fprintf (Out,"\n\n delta eps");
        fprintf (Out,"\n %e %e %e",deps[0],deps[1],deps[2]);
//  stiffness matrix of material
    Mm->matstiff(d,ipp);
    e0=Mm->eliso[imat].e;
    mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
    mxv (d,deps,sig);
// save total strain (from t=0)
    for (i=0; i<nc; i++){Mm->ip[ipp].other[i]=Mm->ip[ipp].strain[i];}// it is last strain in arr. strain
//  time history
	seps_time (screep,sig);
	for (i=0; i<nRetTime; i++){
		ii=nc*(i+1);
	    for (j=0; j<nc; j++){
		 Mm->ip[ipp].other[ii+j]=screep[j][i];
		}
	}
//  Stress data storage first position, Hide (strain) seconde position
//    fprintf (Out,"\n\n *****prir sigma phase 2");
//    fprintf (Out,"\n %e %e %e ",sig[0],sig[1],sig[2]);
//    fprintf (Out,"\n %e %e %e %e %e %e",sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]);
	Mm->ip[ipp].other[nc*(nRetTime+1)+1]=esht;
	Mm->ip[ipp].other[nc*(nRetTime+1)+2]=h_s;
	Mm->ip[ipp].other[nc*(nRetTime+1)+3]=temp_s;
}


void creepbbeam::get_h (long ipp)
{
    h_slast=0.0;
	if(type_h != 1){h_slast=Mm->ip[ipp].other[nc*(nRetTime+1)+2]; }	
}


void creepbbeam::get_temp (long ipp)
{
    templast=0.0;
    if(type_temp != 1){templast=Mm->ip[ipp].other[nc*(nRetTime+1)+2]; }
}

/**
   function apdates material parameters
   
   @param E - function of creep

   10.10.2002
*/
void creepbbeam::matstiff (matrix &d, long ipp)
//  function returns elastic stiffness matrix
//  d - elastic stiffness matrix
{
  long i,j,k;
  double jt,et,pomt,pomt1,qq,delYi,delYj;
  vector bpom(nRetTime);
  matrix apom(nRetTime,nRetTime);
//  matrix c(d.m,d.n);
  
  
  //  initial elastic stiffness matrix
  Mm->elmatstiff (d,ipp);
  //  initial elastic compliance matrix    Mm->elmatcompl (c,ipp);

  e0=Mm->eliso[imat].e;
  mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
  
  //ddTime=Mp->timefun.getval(Mp->time);
  ddTime=Mp->timecon.actualforwtimeincr ();
  //ccTime=Mp->time;
  ccTime=Mp->timecon.actualtime ();
  pomt=ccTime+ddTime/2.;
  qq=2./3.;
  jt=0.0;
  //for (ii=0;ii<numMat;ii++){
  //  long ii;
  //     ii=imat;
  //   from concrete starts   tb 
  if(timemat<=pomt) {
    //        ert[nRetTime+1][imat]=0.;
    et=0.;
    for (i=0;i<nRetTime;i++){
      bpom[i]=0.;
      for (j=0;j<nRetTime;j++){
	apom[i][j]=0.;
      }
    }
    k=0;
    pomt1=pomt+0.0001;   
    while (pomt1<=timeMax*5.){
      //   creep function od casu pomt do casu pomt1
      b3_law(jt,esht, pomt, pomt1);
      //		fprintf (Out,"\n\n jt v case t   od t0 prirustek dt");
      //		fprintf (Out," %e\t%e\t%e\t%e",jt, pomt, ccTime, ddTime);
      for (i=0;i<nRetTime;i++){
	delYi=pow((pomt/retTime[i]),qq)-pow((pomt1/retTime[i]),qq);
	bpom[i]=bpom[i]+(1.-exp(delYi))*jt;
	for (j=0;j<nRetTime;j++){
	  delYj= pow((pomt/retTime[j]),qq)-pow((pomt1/retTime[j]),qq);
	  apom[i][j]=apom[i][j]+(1.-exp(delYi))*(1.-exp(delYj));
	}
      }
      //	      pomt1=pomt1+10**(k)
      k=k+1;
      pomt1=pomt1+20*(k);
    }
    
    inv_sym(apom);
    for (i=0;i<nRetTime;i++){
      //	      ert[i][imat]=0.;
      ert[i]=0.;
      for (j=0;j<nRetTime;j++){
	//	       ert[i][imat]=ert[i][imat]+apom[i][j]*bpom[j];
	ert[i]=ert[i]+apom[i][j]*bpom[j];
      }
      //	      ert[i][imat]=1./ert[i][imat];
      ert[i]=1./ert[i];
      if(ccTime<0.0001){
	delYi=pow((0.0001+ddTime)/retTime[i],qq)-pow(0.0001/retTime[i],qq);}
      else {
	delYi=pow((ccTime+ddTime)/retTime[i],qq)-pow(ccTime/retTime[i],qq);}
      //          lam[i]=1-(1.-(exp(-delYi)))/delYi;
      //          ert[(nRetTime+1][imat]=ert[nRetTime+1][imat]+(delYi-1.+exp(-delYi))/delYi/ert[i][imat];
      et=et+(delYi-1.+exp(-delYi))/delYi/ert[i];
    }
    //	    ert[nRetTime+1][imat]=1./ert[nRetTime+1][imat];
    et=1./et;
    //         call Jkontrola(imat)
    
    
    //      shrink v time    ccTime+ddTime/2.
    pomt1=pomt+0.0001;   
    b3_law(jt,esht, pomt, pomt1);
  }
  
  else {
    //        ert[nRetTime+1][imat]=e0/100000.;
    et=e0/100000.;
    for (i=0;i<nRetTime;i++){ ert[i]=0.;}
  }
  //}
  
  //		fprintf (Out,"\n\n eh v case t   od t0 prirustek dt");
  //		fprintf (Out," %e\t%e\t%e\t%e",et, pomt, ccTime, ddTime);
  
  //  New material modulus in Time
  qq=et/e0;
  
  cmulm( qq, d);
}

/**
   function returns history
   
   @param screep - vector of history
   @param sig - vector od add stress
   
   10.10.2002
*/
void creepbbeam::seps_time (matrix &screep,vector &sig)
{
  int i;
  double qq,dj,pomt;
  
  qq=2./3.;
  if((ccTime+ddTime)==ddTime) {
  }
  else{
    pomt=ccTime+ddTime;
	if(sig.n==4){
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],qq)-pow(ccTime/retTime[i],qq);
		screep[0][i]=exp(-dj)*screep[0][i]+(sig[0]-mi*sig[1])*(1.-exp(-dj))/dj/ert[i];
		screep[1][i]=exp(-dj)*screep[1][i]+(sig[1]-mi*sig[0])*(1.-exp(-dj))/dj/ert[i];
		screep[2][i]=exp(-dj)*screep[2][i]+ sig[2]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
		screep[3][i]=0.0;
	 }
	}
	else if(sig.n==6){
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],qq)-pow(ccTime/retTime[i],qq);
		screep[0][i]=exp(-dj)*screep[0][i]+(sig[0]-mi*sig[1]-mi*sig[2])*(1.-exp(-dj))/dj/ert[i];
		screep[1][i]=exp(-dj)*screep[1][i]+(sig[1]-mi*sig[0]-mi*sig[2])*(1.-exp(-dj))/dj/ert[i];
		screep[2][i]=exp(-dj)*screep[2][i]+(sig[2]-mi*sig[0]-mi*sig[1])*(1.-exp(-dj))/dj/ert[i];
		screep[3][i]=exp(-dj)*screep[3][i]+ sig[3]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
		screep[4][i]=exp(-dj)*screep[4][i]+ sig[4]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
		screep[5][i]=exp(-dj)*screep[5][i]+ sig[5]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
	 }
	}
  }

}


/**
   function computes J(t,t') DOUBLE POWER LAW by Prof. Bazant
   
   @param jt - function of creep

   10.10.2002
*/
void creepbbeam::b2_law (double &jt,double &des_hn, double t0, double t )
{
  double n, c, c0,cd,  ac,ag,m,x,fi,esn,esn1,tau,sd,st,st0, kt,kh,et1,et2;
  
//  K_s shape factor slab=1.0, cylinder=1.15, sguare prism.=1.25, sphere=1.3, cube=1.55
  
  //  stiffness matrix of material
//  Mm->matstiff (d,ipp);
/*
  from t0  to t 
  from concrete starts   tb 
  t_w  age when drying begins 
  (fc') is 28 day average cilinder strenght fc' [ksi] ksi=1000psi=6.895 MPa(f.e.6.454=44.5MPa)***6.381
  (w/c) is water-cement ratio of the mix by weight   ***0.43
  (s/c) is send-cement ratio of the mix by weight    ***3.4
  (g/c) is gravel-cement ratio of the mix by weight g/c=a/c-s/c     ***1.98
  (a/c) is aggregate-cement ratio of the mix by weight a/c=g/c+s/c   
  (a1)  is coef. for cements of type I,II a1=1.00, III a1=0.93, IV a1=1.05   ***1.05
  (ro)  is mass of concrete in [lb/ft3] =16.03 kg/m3 ***156
  (tl)  effective cross section thickness       D=2*vs_s
  cs  cement content in m3  .. kg/m3
     E0=(0.09+1/(1.7*(0.5*ro*ro*fc*1e-4)*(0.5*ro*ro*fc*1e-4)))
     Et=E0*sqrt(t/(4+0.85*t))          podle ACI Commite 209/II

*/
  
  
  if (tb<t0){
    if (t_w<tb) t_w = tb;
    if (t_w<=0.0) t_w = 1.0;
    ac=sc+gc; ag=ac/gc; m=-0.28-47.541025/fc/fc;
    
    x=(2.1*ac/pow(sc,1.4)+0.1*pow((fc/6895.),1.5)*pow(wc,(1./3.))*pow(ag,2.2))*a1-4.;
    if (x<=4){n=0.12;}
    else {n=0.12+(0.07*pow(x,6.0)/(5130.0+pow(x,6.0)));}
  
    fi=0.5*pow(10.,(3.*n))/(pow(28.,m)+0.025/wc);
    c0=fi*(pow((t0-tb),m)+0.025/wc)*pow((t-t0),n);
  
  //  shrink and drying
    x=( (1.25*sqrt(ac)+0.5*pow((gc/sc),2))*pow(((1.+sc)/wc),(1/3.))*sqrt(fc/6895.) )-12;
    if (x>=0.0){esn=0.0121;}
    else {esn=0.0121-0.0088/(390./pow(x,4)+1.);}

    c=(0.000000125*wc*c_s-0.000012);
    if (c<0.000007) c=0.000007;
    if (c>0.000021) c=0.000021;
    kt=(temp_s+273.)/296*exp(5000/296.-5000./(temp_s+273.));   //for temperature 23C, 296K
    c=c*kt*(0.05+sqrt(6.3/t_w));
    tau=0.267*(k_s*k_d)*(k_s*k_d)/c;
    if (5<=2){
//    if ( (t0>t_w) && (t<=(t_w+tau)) ){
  //  shrink
      et1=e0/(1.+pow(0.1,(3.*n))*fi*( pow(607.,m)    +0.025/wc));
      et2=e0/(1.+pow(0.1,(3.*n))*fi*(pow((t_w+tau),m)+0.025/wc));
      kh=1-pow(h_s,3.)+pow((1-h_s),5.);
      st= 1./(1.+pow(pow(tau/(t-t_w),r_s),0.5*r_s) );
      st0= 1./(1.+pow(pow(tau/(t0-t_w),r_s),0.5*r_s) );
      esn1=esn*et1/et2*st*kh/1000000.0;
  //  drying
      x=56000*pow((ac*fc/6895.),0.3)*pow(gc,1.3)*pow((wc/esn),1.5)-0.85;
      if (x>=0.0){fi=0.008+0.027*(1./(1.+0.7/pow(x,1.4)));  }
      else fi=0.008;

      fi = 1./sqrt(1.+(t-t0)/10./tau)*fi;
      sd = 1./pow((1.+10.*tau/(t-t_w)),(2.6*n - 6.*n*n));
      //cd = fi*pow(t_w,(m/2.))*(pow(h_0,1.5) - pow(h_s,1.5))*esn1*sd;
      cd = fi*pow(t_w,(m/2.))*(1.0 - pow(h_s,1.5))*esn1*sd;
	  des_hn=-esn*et1/et2*(st0-st)*kh/1000000.0;
	}
    else {
	  cd = 0.0;
      des_hn=0.0;
	}
  }
  else {
    cd = 0.0;
    des_hn=0.0;
  }
  //  measure static modulus - creep modulus  1/Ec=2/(3 Es)
  //  1/E=e-6/psi =e-6/6.895kPa
  jt=(1+c0+cd)/(1.5*e0);
//    fprintf (Out,"\n t0, t,  c0");
//    fprintf (Out," %e %e\t%e",t0,t,c0);
}


void creepbbeam::b3_law (double &jt,double &des_hn, double t0, double t)
{
  double n,c0,cd,  ac,ag,m,x,y,q,q1,q2,q3,q4,q5,esn,kh,tau,ht,ht0,dh_s,dtemp;
  short experiment=0;
//  call only from matstiff 
//  shape factor slab=1.0, cylinder=1.15, sguare prism.=1.25, sphere=1.3, cube=1.55
  
/*
  from t0  to t 
  from concrete starts   tb 
  t_w  age when drying begins 
  (fc') is 28 day average cilinder strenght fc' [ksi] ksi=1000psi=6.895 MPa(f.e.6.454=44.5MPa)***6.381
  (w/c) is water-cement ratio of the mix by weight   ***0.43
  (s/c) is send-cement ratio of the mix by weight    ***3.4
  (g/c) is gravel-cement ratio of the mix by weight g/c=a/c-s/c     ***1.98
  (a/c) is aggregate-cement ratio of the mix by weight a/c=g/c+s/c   
  (a1)  is coef. for cements of type I,II a1=1.00, III a1=0.93, IV a1=1.05   ***1.05
  (ro)  is mass of concrete in [lb/ft3] =16.03 kg/m3 ***156
  (k_d)  effective cross section thickness       D=2*V/S
  cs  cement content in m3  .. kg/m3
     E0=(0.09+1/(1.7*(0.5*ro*ro*fc*1e-4)*(0.5*ro*ro*fc*1e-4)))
     Et=E0*sqrt(t/(4+0.85*t))          podle ACI Commite 209/II

*/
  
  if (tb<t0){
   if (t_w<tb) t_w = tb;
   if (t_w<=0.0) t_w = 1.0;
   ac=sc+gc;   ag=ac/gc;   m=0.5; //m=0.28+47541025/fc/fc;
    
   x=(2.1*ac/pow(sc,1.4)+0.1*pow((fc/6895),1.5)*pow(wc,(1.0/3.0))*pow(ag,2.2))*a1-4;
   if (x<=4){n=0.12;}
   else {n=0.12+(0.07*pow(x,6.0)/(5130.0+pow(x,6.0)));}
   q1=600000.0/4734.0/sqrt(fc);
   x=1.7*pow(t0,0.12)+8.0;
   y=(0.086*pow(t0,(2.0/9.0))+1.21*pow(t0,(4.0/9.0)));
   q=1.0/y/pow( ( 1.0+ pow( (pow(t0,m)/log(1.+pow((t-t0),n))/y),x ) ),(1./x));
   q2=185.4*sqrt(c_s)/pow(fc,0.9);
   q3=0.29*pow(wc,4.0)*q2;
   q4=20.3/pow(ac,0.7);
  
   c0=q2*q+q3*log(1.0+pow((t-t0),n))+q4*log(t/t0);
 
  //  shrink and drying
   if ( type_h ==1){
      if ( (t0>=t_w) ){        //   if ( (t0>t_w) && (t<=(t_w+600.)) ){
  //  shrink from humidity=1 to h_s
        tau=k_s*k_s*k_d*k_d*85000./pow(t,0.08)/pow(fc,0.25);
        esn=-a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
	    ht0 =tanh(sqrt((t0-t_w)/tau));
	    ht  =tanh(sqrt((t -t_w)/tau));
		if(h_s<=0.98){kh=(1.0-pow(h_s,3.0));}
		else {kh=-10*(1.-h_s);}
        des_hn = esn*kh/1000000.0*(ht-ht0); //increment of eps
  //  drying from humidity=1 to h_s 
        q5 = 0.757/fc/pow(fabs(esn),0.6);
        cd = q5*sqrt( -exp(-8+8.*(1.0-h_s)*ht)+exp(-8.+8.*(1.0-h_s)*ht0) );
	  }
	  else {
		cd = 0.0;
		des_hn=0.0;
	  }
   }
   else {
  //  shrink and drying from computing
     dh_s = h_s-h_slast;   // humidity
     esn = -a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
	 if(h_s<=0.98){kh=pow(h_s,3.0)-pow(h_slast,3.0);}
	 else {kh=-10.0*dh_s;}
  //  shrink
     des_hn = esn*kh/1000000.0; //increment of eps
  //  drying
     q5 = 0.757/fc/pow(fabs(esn),0.6);
     cd = -q5*dh_s;
   }
   if ( type_temp ==1 ){
  //  temperature constant
   }
   else {
  //  temperature
     dtemp = temp_s-templast;   
     esn=alfa;
  //  shrink
     des_hn=des_hn+esn*dtemp;
  //  drying
	 q5= 1.5/fc;
     cd = cd+q5*esn*dtemp;
   }
  }

  else {
    q1 = 0.0;
    c0 = 0.0;
    cd = 0.0;
    des_hn=0.0;
  }
  //  measure static modulus - creep modulus  1/Ec=2/(3 Es)
  //  1/E=e-6/psi =e-6/6.895kPa
  jt=(q1+c0+cd)*1.e-6;
//    fprintf (Out,"\n t0, t,  q1,q2,q3,q4");
//    fprintf (Out," %e %e\t%e\t%e\t%e\t%e",t0,t,q1,q2*q,q3*log(1+pow((t-t0),n)),q4*log((t-tb)/(t0-tb)));
}
