#include "creepbs.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "globmat.h"
#include "intpoints.h"
#include "elastisomat.h"
#include "probdesc.h"
#include "mechmat.h"
#include <stdlib.h>
#include <math.h>


creepbs::creepbs (void)
{
  timeMax=Mp->timecon.endtime ();
  double m = log10(timeMax);
  nRetTime=7;  
  allocv (nRetTime,retTime);
  allocv (nRetTime,ert);
  retTime[0]=1.0e-9;
  retTime[1]=1.0e-3;
  retTime[2]=1.0;
                  //  retTime[3]=14.0;  //retTime[4]=100.0; //retTime[5]=800.0;
  retTime[3]=pow(10.,0.25*m);
  retTime[4]=pow(10.,0.5*m);
  retTime[5]=pow(10.,0.75*m);
  retTime[6]=timeMax;
  type_h=1;
  type_temp=1;
  h_s = 0.4;   // end humidity
  desht=0.0;
  k_s = 1.0;	   //shape factor
  k_d = 1.0;      // thicknes
  timemat=0.0;
  alfa=0.0;
  allocm(nRetTime, nRetTime, apom); 
  ccTime = -1.0e30;
}

creepbs::~creepbs (void)
{
}

void creepbs::creepinit (long ipp,double val,nonmechquant nmq)
{
  long nc=Mm->ip[ipp].ncompstr;

  if(nmq == rel_hum)       
    Mm->ip[ipp].eqother[nc*(nRetTime+1)+1]=val;

  if(nmq == temperature )
    Mm->ip[ipp].eqother[nc*(nRetTime+1)+2]=val;
}

void creepbs::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %ld %ld",
	          &tb,&t_w,&fc,&wc,&sc,&gc,&c_s,&a1,&k_d,&type_h,&type_temp);
  
  if (type_h==1){xfscanf (in,"%lf ",&h_s); }
  if (type_temp==1){xfscanf (in,"%lf ",&temp_s); }

}
/**
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

*/
double creepbs::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function inversion sym. matrix

   @param a - matrix

*/
void creepbs::inv_sym (matrix &a)
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

void creepbs::updateval()
{
 if(type_h!=1&&type_temp!=1) updatevalchange();
 else updatevalkons();
 
}
void creepbs::updatevalchange()
{
  long i,j,k;
  double pomt,pomt1,qq,delYi,delYj,m0,m;
  
  if (ccTime != Mp->timecon.actualtime ())
  {
    //ddTime=Mp->timefun.getval(Mp->time);
    ddTime=Mp->timecon.actualforwtimeincr ();
    //ccTime=Mp->time;
    ccTime=Mp->timecon.actualtime ();
    pomt=ccTime+ddTime/2.;
    qq=2./3.;
    //   from concrete starts   tb 
    fillm(0.0, apom); 
    if (timemat<=pomt) {
      pomt1=pomt+0.0001;   
      m0 = log10(pomt1);
      m = log10(2.*timeMax)*0.999;
      k=0;
      while (pomt1<=2.*timeMax)
      {
        for (i=0;i<nRetTime;i++)
        {  
	  delYi=pow((pomt/retTime[i]),qq)-pow((pomt1/retTime[i]),qq);
	  for (j=0;j<nRetTime;j++)
	  {
	    delYj= pow((pomt/retTime[j]),qq)-pow((pomt1/retTime[j]),qq);
	    apom[i][j]=apom[i][j]+(1.-exp(delYi))*(1.-exp(delYj));
	  }
        }
        k=k+1;
        pomt1=pow(10.,m0+(m-m0)/40*k);
      }
      inv_sym(apom);
    }
    fprintf(stdout, "\n ### inversion\n");
  }
}
void creepbs::updatevalkons()
{
  long i,j,k;
  double pomt,pomt1,qq,jt,delYi,delYj,m0,m;
  vector bpom(nRetTime);
  
  if (ccTime != Mp->timecon.actualtime ())
  {
    //ddTime=Mp->timefun.getval(Mp->time);
    ddTime=Mp->timecon.actualforwtimeincr ();
    //ccTime=Mp->time;
    ccTime=Mp->timecon.actualtime ();
    pomt=ccTime+ddTime/2.;
    qq=2./3.;
     jt=0.0;
  //   from concrete starts   tb 
   if(timemat<=pomt) {
    et=0.;
    for (i=0;i<nRetTime;i++){
      bpom[i]=0.;
      for (j=0;j<nRetTime;j++){
	    apom[i][j]=0.;
      }
    }

    pomt1=pomt+0.0001;   
    m0 = log10(pomt1);
    m = log10(2.*timeMax)*0.999;
    k=0;
    while (pomt1<=2.*timeMax){
      //   creep function from time pomt to time pomt1
      {b3_law(jt, pomt, pomt1);}
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
      k=k+1;
	  pomt1=pow(10.,m0+(m-m0)/40*k);
//      pomt1=pomt1+20*(k);
    }
    
    inv_sym(apom);
    for (i=0;i<nRetTime;i++){
      ert[i]=0.;
      for (j=0;j<nRetTime;j++){
	           ert[i]=ert[i]+apom[i][j]*bpom[j];
      }
      ert[i]=1./ert[i];
      if(ccTime<0.0001){
	    delYi=pow((0.0001+ddTime)/retTime[i],qq)-pow(0.0001/retTime[i],qq);}
      else {
	    delYi=pow((ccTime+ddTime)/retTime[i],qq)-pow(ccTime/retTime[i],qq);}
        et=et+(delYi-1.+exp(-delYi))/delYi/ert[i];
    }
    if(et>1.e-50)  et=1./et;
    else{fprintf (Out,"\n\n Stiff is zero"); abort();}
    //         call Jkontrola(imat)
    
   }
  
   else {
//    et=e0/100000.;
    et=1.;
    for (i=0;i<nRetTime;i++){ ert[i]=0.;}
   }
	
  }

}


/**
   function returns eps from history
   
   @param screep - vector of history
   @param epsscr - vector deformation of history
   
   10.10.2002
*/
void creepbs::nlstresses (long ipp)
{
  //  number of components od strain/stress arrays
  nc=Mm->ip[ipp].ncompstr;
  
  //  type of stress state
  //bar=1,plbeam=2,spacebeam=5,
  //planestress=10,planestrain=11,plate=15,
  //axisymm=20,shell=25,spacestress=30
  ss=Mm->ip[ipp].ssst;  
  
  //timeMax=Mp->end_time;
  //ddTime=Mp->timefun.getval(Mp->time);
  ddTime=Mp->timecon.actualforwtimeincr ();
  //ccTime=Mp->time;
  ccTime=Mp->timecon.actualtime ();
  
  //  id of elastic material
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  
  //  actual moisture
  if(type_h!=1){
    h_s=Mm->givenonmechq(rel_hum, ipp);
  }
  
  //  actual temperature
  if(type_temp!=1){
    temp_s=Mm->givenonmechq(temperature, ipp);
  }
  
  get_h (ipp);     //initial  moisture
  get_temp (ipp);     //initial  temperature
  
  if (Mp->phase==1){
    phase1(ipp);
  }
  
  if (Mp->phase==2){
    phase2(ipp);
  }
  
}


void creepbs::phase1 (long ipp)
{
  //  right hand side from time dependent value computation
  //  creep function for shrink from pomt to cctime
  long i,j,ii;
  double qq,dj,endTime;
  
  vector epscr(nc),depsshtem(nc),sig(nc);
  matrix d(nc,nc),screep(nc,nRetTime);
  
  //    Mm->matstiff(d,ipp);  //compute desht
  //    e0=Mm->eliso[imat].e;
  //    mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
  endTime=ccTime+ddTime;
  get_desht (desht, ccTime, endTime);
  for (i=0; i<nc; i++){
    epscr[i]=0.0;
  }
  qq=2.0/3.0;
  for (i=0;i<nRetTime;i++){
    ii=nc*(i+1);
    dj=pow(endTime/retTime[i],qq)-pow(ccTime/retTime[i],qq);
    for (j=0; j<nc; j++){
      screep[j][i]=Mm->ip[ipp].eqother[ii+j];	
      epscr[j]=epscr[j]+(1.-exp(-dj))*screep[j][i];
    }
  }
  //        fprintf (Out,"\n\n v case t  prirustek dt faze 1");
  //        fprintf (Out," %e\t%e", ccTime, ddTime);
  //        fprintf (Out,"\n\n eps faze1");
  //        fprintf (Out," %le %le %le %le %le %le",epscr[0],epscr[1],epscr[2],epscr[3],epscr[4],epscr[5]);
  //        fprintf (Out,"\n\n tau faze1");
  //        fprintf (Out," %le %le %le %le %le %le",screep[0][0],screep[1][0],screep[2][0],screep[3][0],screep[4][0],screep[5][0]);
  //  increment of shrink and temperature
  
  //  jen kvuli debuggovani
  //desht=0.0;

  /*  epscr[0]=epscr[0]+desht;
      epscr[1]=epscr[1]+desht;
      if(ss==planestrain) epscr[2]=epscr[2]+desht;
      else if(ss==axisymm) epscr[2]=epscr[2]+desht;
      else if(ss==spacestress) epscr[2]=epscr[2]+desht;

  */
  //new for steel!!
  epscr[0]= 0.0;
  epscr[1]= 0.0;
  if(ss==planestrain) epscr[2]=0.0;
  else if(ss==axisymm) epscr[2]=0.0;
  else if(ss==spacestress) epscr[2]=0.0;
  
  //	fprintf (Out,"\n\n desht prir faze 1"); fprintf (Out," %le",desht);
  //	fprintf (Out,"\n h_s prir faze 1"); fprintf (Out," %le %le",h_s, h_slast);
  //	fprintf (Out,"\n temp_s prir faze 1"); fprintf (Out," %le %le",temp_s, temp_slast);
  
  //		fprintf (Out,"\n\n phase 1 eps1, eps2 prir");
  //		fprintf (Out," %e %e",epscr[0], epscr[1]);
  
  Mm->stiff_deps_creep (ipp, 0,epscr, nc,nc);  // save sig
}



void creepbs::phase2 (long ipp)
{
//  new total strain and stress
  long i,j,ii;
  double qq,dj, endTime;
  vector epscr(nc),deps(nc),depsshtem(nc),sig(nc);
  matrix d(nc,nc),screep(nc,nRetTime);
  
  endTime=ccTime+ddTime;
    fprintf (Out,"\n\n phase 2 v case t   prirustek dt");
    fprintf (Out," %e\t%e", ccTime, ddTime);
  // arr. strain = is total strain (including e from screep) from t=0
  for (i=0; i<nc; i++){
    deps[i]=Mm->ip[ipp].strain[i]-Mm->ip[ipp].eqother[i];
  } //increment deps
  qq=2.0/3.0;
  for (i=0;i<nRetTime;i++){
    ii=nc*(i+1);
    dj=pow(endTime/retTime[i],qq)-pow(ccTime/retTime[i],qq);
    for (j=0; j<nc; j++){
      screep[j][i]=Mm->ip[ipp].eqother[ii+j];
      deps[j]=deps[j]-(1.-exp(-dj))*screep[j][i];     // deps=deps-e(screep)
    }
  }
  get_desht (desht, ccTime, endTime);
  //  fprintf (Out,"\n\n desht prir faze 1"); fprintf (Out," %le",desht);
  
  //  debuggovani
  //desht=0.0;

  /*  deps[0]=deps[0]-desht;
      deps[1]=deps[1]-desht;
      if(ss==planestrain) deps[2]=deps[2]-desht;
      else if(ss==axisymm) deps[2]=deps[2]-desht;
      else if(ss==spacestress) deps[2]=deps[2]-desht;
  */
  //new for steel!!
  deps[0]=0.0;
  deps[1]=0.0;
  if(ss==planestrain) deps[2]=0.0;
  else if(ss==axisymm) deps[2]=0.0;
  else if(ss==spacestress) deps[2]=0.0;
  


  //  increment of stress components
  //        fprintf (Out,"\n\n eps celkove,   eps celkove minule,  delta eps");
  //        fprintf (Out,"\n %le %le %le %le %le %le",epscr[0],epscr[1],epscr[2],epscr[3],epscr[4],epscr[5]);
  //        fprintf (Out,"\n %le %le %le ",Mm->ip[ipp].eqother[0],Mm->ip[ipp].eqother[1],Mm->ip[ipp].eqother[2]);
  //        fprintf (Out,"\n %le %le %le ",Mm->ip[ipp].eqother[3],Mm->ip[ipp].eqother[4],Mm->ip[ipp].eqother[5]);
  //        fprintf (Out,"\n\n delta eps");
  //        fprintf (Out,"\n %e %e %e",deps[0],deps[1],deps[2]);
  //  stiffness matrix of material
  Mm->matstiff(d,ipp);  
  e0=Mm->eliso[imat].e;
  mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
  mxv (d,deps,sig);
  // save total strain (from t=0)
  for (i=0; i<nc; i++){Mm->ip[ipp].eqother[i]=Mm->ip[ipp].strain[i];}// it is last strain in arr. strain
  //  time history
  seps_time (screep,sig);
  for (i=0; i<nRetTime; i++){
    ii=nc*(i+1);
    for (j=0; j<nc; j++){
      Mm->ip[ipp].eqother[ii+j]=screep[j][i];
    }
  }
  //       fprintf (Out,"\n\n tau faze2");
  //       fprintf (Out," %le %le %le %le %le %le",screep[0][0],screep[1][0],screep[2][0],screep[3][0],screep[4][0],screep[5][0]);
  //  Stress data storage first position, Hide (strain) seconde position
  //    fprintf (Out,"\n\n *****prir sigma phase 2");
  //    fprintf (Out,"\n %e %e %e ",sig[0],sig[1],sig[2]);
  //    fprintf (Out,"\n %e %e %e %e %e %e",sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]);
  Mm->ip[ipp].eqother[nc*(nRetTime+1)+0]+=desht;
  Mm->ip[ipp].eqother[nc*(nRetTime+1)+1]=h_s;
  Mm->ip[ipp].eqother[nc*(nRetTime+1)+2]=temp_s;
  //       fprintf (Out,"\n\n desht hs, temps faze2");
  //       fprintf (Out," %le  %le %le ",desht, h_s, temp_s);
}

void creepbs::get_h (long ipp)
{
  h_slast=1.0;
  if(type_h !=1 ){h_slast=Mm->ip[ipp].eqother[nc*(nRetTime+1)+1]; }
}


void creepbs::get_temp (long ipp)
{
  temp_slast=283.0;   //10 C
  if(type_temp != 1 ){temp_slast=Mm->ip[ipp].eqother[nc*(nRetTime+1)+2]; }
}



/**
   function apdates material parameters
   
   @param E - function of creep

   10.10.2002
*/
void creepbs::matstiff (matrix &d, long ipp)
//  function returns elastic stiffness matrix
//  d - elastic stiffness matrix
{
  double qq=1.;
  
  nc=Mm->ip[ipp].ncompstr;
  ss=Mm->ip[ipp].ssst;  //planestress=10,planestrain=11,plate=15,
                        //axisymm=20,shell=25,spacestress=30
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
   //  initial elastic stiffness matrix
  Mm->elmatstiff (d,ipp);
  //  initial elastic compliance matrix    Mm->elmatcompl (c,ipp);
  e0=Mm->eliso[imat].e;
  mi=Mm->eliso[imat].nu;    //    d[1][0] = nu*e/(1.0-nu*nu);
  
 if(type_h!=1 && type_temp!=1) matstiffchange(qq, ipp);
 else matstiffkons(qq);
  //  New material modulus in Time
 //new for steel!!
 qq = 1.0;

  cmulm( qq, d);
}
void creepbs::matstiffkons (double &qq)
//  function returns elastic stiffness matrix
{
  //  New material modulus in Time
  qq=et/e0;
}

void creepbs::matstiffchange (double &qq, long ipp)
//  function returns elastic stiffness matrix
//  d - elastic stiffness matrix
{
  long i,j,k;
  double jt,pomt,pomt1,q,delYi,m0,m;
  vector bpom(nRetTime);
  
  if(Mm->ip[ipp].hmt & 1) 
	{i=Mm->ip[ipp].idm[Mm->ip[ipp].nm-1]; alfa=Mm->tidilat[i].alpha;}
  // initial and actual moisture and temperature
  if(type_h!=1) {h_s=Mm->givenonmechq(rel_hum, ipp);}   //actual moisture
  if(type_temp!=1) {temp_s=Mm->givenonmechq(temperature, ipp);}  //actual temperature
  get_h (ipp);     //initial  moisture
  get_temp (ipp);     //initial  temperature

  pomt=ccTime+ddTime/2.;
  q=2./3.;
  jt=0.0;
  //   from concrete starts   tb 
  if(timemat<=pomt) {
    et=0.;
    for (i=0;i<nRetTime;i++){
      bpom[i]=0.;
/*      for (j=0;j<nRetTime;j++){
	    apom[i][j]=0.;
      }*/
    }

    pomt1=pomt+0.0001;   
    m0 = log10(pomt1);
    m = log10(2.*timeMax)*0.999;
    k=0;
    while (pomt1<=2.*timeMax){
      //   creep function from time pomt to time pomt1
      {b3_law(jt, pomt, pomt1);}
      //		fprintf (Out,"\n\n jt v case t   od t0 prirustek dt");
      //		fprintf (Out," %e\t%e\t%e\t%e",jt, pomt, ccTime, ddTime);
      for (i=0;i<nRetTime;i++){
	   delYi=pow((pomt/retTime[i]),q)-pow((pomt1/retTime[i]),q);
	   bpom[i]=bpom[i]+(1.-exp(delYi))*jt;
	   for (j=0;j<nRetTime;j++){
//	  delYj= pow((pomt/retTime[j]),q)-pow((pomt1/retTime[j]),q);
//	  apom[i][j]=apom[i][j]+(1.-exp(delYi))*(1.-exp(delYj));
	   }
      }
      k=k+1;
	  pomt1=pow(10.,m0+(m-m0)/40*k);
//      pomt1=pomt1+20*(k);
    }
    
//    inv_sym(apom);
    for (i=0;i<nRetTime;i++){
      ert[i]=0.;
      for (j=0;j<nRetTime;j++){
	   ert[i]=ert[i]+apom[i][j]*bpom[j];
      }
      ert[i]=1./ert[i];
      if(ccTime<0.0001){
	   delYi=pow((0.0001+ddTime)/retTime[i],q)-pow(0.0001/retTime[i],q);}
      else {
	   delYi=pow((ccTime+ddTime)/retTime[i],q)-pow(ccTime/retTime[i],q);}
       et=et+(delYi-1.+exp(-delYi))/delYi/ert[i];
    }
    if(et>1.e-50)  et=1./et;
    else{fprintf (Out,"\n\n Stiff is zero"); abort();}
  }
  
  else {
    et=1.;
    for (i=0;i<nRetTime;i++){ ert[i]=0.;}
  }
  
  //		fprintf (Out,"\n\n eh v case t   od t0 prirustek dt");
  //		fprintf (Out," %e\t%e\t%e\t%e",et, pomt, ccTime, ddTime);
  
  //  New material modulus in Time
  qq=et/e0;


}

/**
   function returns history
   
   @param screep - vector of history
   @param sig - vector od add stress
   
   10.10.2002
*/
void creepbs::seps_time (matrix &screep,vector &sig)
{
  int i;
  double q,dj,pomt,pp;
  
  q=2./3.;
  if((ccTime+ddTime)==ddTime) {
  }
  else{
    pomt=ccTime+ddTime;
	if(ss==planestress){
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],q)-pow(ccTime/retTime[i],q);
		screep[0][i]=exp(-dj)*screep[0][i]+(sig[0]-mi*sig[1])*(1.-exp(-dj))/dj/ert[i];
		screep[1][i]=exp(-dj)*screep[1][i]+(sig[1]-mi*sig[0])*(1.-exp(-dj))/dj/ert[i];
		screep[2][i]=exp(-dj)*screep[2][i]+ sig[2]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
		screep[3][i]=0.0;
	 }
	}
	else if(ss==planestrain){  //sx,sy,sy,txy
	 pp=1.+mi;
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],q)-pow(ccTime/retTime[i],q);
		screep[0][i]=exp(-dj)*screep[0][i]+ pp*(sig[0]*(1.-mi)-mi*sig[1])*(1.-exp(-dj))/dj/ert[i];
		screep[1][i]=exp(-dj)*screep[1][i]+ pp*(sig[1]*(1.-mi)-mi*sig[0])*(1.-exp(-dj))/dj/ert[i];
		screep[2][i]=exp(-dj)*screep[2][i]+ pp*sig[2]*2.   *(1.-exp(-dj))/dj/ert[i];
		screep[3][i]=0.0;
	 }
	}
	else if(ss==axisymm){
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],q)-pow(ccTime/retTime[i],q);
		screep[0][i]=exp(-dj)*screep[0][i]+(sig[0]-mi*sig[1]-mi*sig[2])*(1.-exp(-dj))/dj/ert[i];
		screep[1][i]=exp(-dj)*screep[1][i]+(sig[1]-mi*sig[0]-mi*sig[2])*(1.-exp(-dj))/dj/ert[i];
		screep[2][i]=exp(-dj)*screep[2][i]+(sig[2]-mi*sig[0]-mi*sig[1])*(1.-exp(-dj))/dj/ert[i];
		screep[3][i]=exp(-dj)*screep[3][i]+ sig[3]*(1.+mi)   *(1.-exp(-dj))/dj/ert[i];
	 }
	}
	else if(ss==spacestress){
	 for (i=0;i<nRetTime;i++){
		dj=pow(pomt/retTime[i],q)-pow(ccTime/retTime[i],q);
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
   function computes eps  by Prof. Bazant
   
   @param des_hn - function computing of shrink and temperature strain

   10.9.2004
*/
void creepbs::get_desht (double &des_hn, double t0, double t)
{
  double esn,kh,tau,ht,ht0,dh_s,dtemp;

/*
  from t0  to t 
  k_s   shape factor slab=1.0, cylinder=1.15, sguare prism.=1.25, sphere=1.3, cube=1.55
  tb    from concrete starts
  t_w   age when drying begins 
  (fc') is 28 day average cilinder strenght fc' [ksi] ksi=1000psi=6.895 MPa(f.e.6.454=44.5MPa)***6.381
  (k_d) effective cross section thickness       D=2*vs_s
   h_s  actual relative humidity 
   h_slast   relative humidity at time t0
   temp_s    actual temperature in K 
   temp_slast  temperature at time t0
*/
  
  if (tb<=t0){
 
  //  shrink
   if ( type_h ==1){
      if ( (t0>=t_w) ){        //   if ( (t0>t_w) && (t<=(t_w+600.)) ){
  //  shrink from huminidy   from  1 to h_s
        tau=k_s*k_s*k_d*k_d*85000./pow(t,0.08)/pow(fc,0.25);
        esn=-a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
	    ht0 =tanh(sqrt((t0-t_w)/tau));
	    ht  =tanh(sqrt((t -t_w)/tau));
		if(h_s<=0.98){kh=(1.0-pow(h_s,3.0));}
		else {kh=-10.0*(1.-h_s);}
        des_hn = esn*kh/1000000.0*(ht-ht0); //increment of eps
	  }
	  else {
		des_hn=0.0;
	  }
   }
   else {
  //  shrink from changeing humidity
     dh_s = h_s-h_slast;   // humidity
     esn = a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
     kh = 3.0*h_s*h_s;
     des_hn = esn*kh*dh_s*1.e-6;
   }
   if ( type_temp ==1 ){
  //  temperature constant
   }
   else {
  //  temperature
     dtemp = temp_s-temp_slast;   
     esn=alfa*1.e6;
     des_hn=des_hn+esn*dtemp*1.e-6;
   }
  }

  else {
    des_hn=0.0;
  }

  des_hn = 0.0;//tady ??!!

}


/**
   function computes J(t,t')  by Prof. Bazant    (Jt = strain)
   
   @param jt - function of creep

   10.10.2002
*/
void creepbs::b3_law (double &jt, double t0, double t)
{
  double n,c0,cd,  ac,ag,m,x,y,q,q1,q2,q3,q4,q5,esn,kh,tau,ht,ht0,dh_s,dtemp;
//  call only from matstiff 
/*
  from t0  to t 
  K_s  shape factor slab=1.0, cylinder=1.15, sguare prism.=1.25, sphere=1.3, cube=1.55 
  tb   from concrete starts
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
  //  shrink from huminidy=1 to h_s
        tau=k_s*k_s*k_d*k_d*85000./pow(t,0.08)/pow(fc,0.25);
        esn=a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
	    ht0 =tanh(sqrt((t0-t_w)/tau));
	    ht  =tanh(sqrt((t -t_w)/tau));
		if(h_s<=0.98){kh=(1.0-pow(h_s,3.0));}
		else {kh=-10.0*(1.-h_s);}
  //  drying from huminidy=1 to h_s 
        q5 = 0.757/fc/pow(esn,0.6);
        cd = q5*sqrt( exp(-8+8.*(1.0-h_s)*ht)-exp(-8.+8.*(1.0-h_s)*ht0) );
	  }
	  else {
		cd = 0.0;
	  }
   }
   else {
  //  shrink and drying from changeing humidity 
     dh_s = h_s-h_slast;   // humidity
     esn = a1*1.2*(0.019*pow((wc*c_s),2.1)/pow(fc,0.28)+270.);
     kh = 3.0*h_s*h_s;
	 q5 = 0.35/fc;
     cd = fabs(q5*esn*kh*dh_s);
   }
   if ( type_temp ==1 ){
  //  temperature constant
   }
   else {
  //  temperature
     dtemp = temp_s-temp_slast;   
     esn=alfa*1.e6;
	 q5= 1.5/fc;
     cd = cd+fabs(q5*esn*dtemp);
   }
  }

  else {
    q1 = 0.0;
    c0 = 0.0;
    cd = 0.0;
  }
  //  measure static modulus - creep modulus  1/Ec=2/(3 Es)
  //  1/E=e-6/psi =e-6/6.895kPa
  //Jt = strain
  jt=(q1+c0+cd)*1.e-6;     
//    fprintf (Out,"\n t0, t,  q1,q2,q3,q4");
//    fprintf (Out," %e %e\t%e\t%e\t%e\t%e",t0,t,q1,q2*q,q3*log(1+pow((t-t0),n)),q4*log((t-tb)/(t0-tb)));
}
