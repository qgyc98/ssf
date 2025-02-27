#include "eigsol.h"
#include "vector.h"
#include <math.h>



/**
   function computes all eigenvalues and eigenvectors of the matrix
   matrix is stored as dense matrix and will be overwritten!
   
   eigenvectors are columns of the matrix evec but matrix is stored
   in usual way, by rows.
   
   @param a - array containing matrix
   @param evec - array containing eigenvectors
   @param eval - array containing eigenvalues
   @param n - order of matrix a (number of rows or columns)
   @param ni - maximum number of iterations
   @param ani - number of performed iterations
   @param limit - maximum acceptable absolute value of offdiagonal element
   @param normalize - switch for normalizing of matrix (on>0, off=0)
   
   !!!! nutno predelat trideni vlastnich vektoru
   
   23.7.2001
   updated 17.8.2007 by TKo. - normalizing of matrix
*/
void jacobi_rot (double *a,double *evec,double *eval,long n,long ni,long &ani,double limit, long normalize)
{
  long i,j,k,l,ii,jj,en,nn;
  double c,s,q,r,t,ai,aj,maxa;
  
  if (normalize)
  {
    maxa = 0.0; 
    for (i=0; i<n*n; i++) // search any nonzero element for initial value of max
    {
      if (a[0] != 0.0)
      {
	maxa = a[0];
        break;;
      }
    }    
    if (maxa != 0.0) // if matrix is nonzero then normalize
    {
      for (i=0; i<n*n; i++) // find maximal nonzero element of matrix
      {
        if ((fabs(a[i]) > fabs(maxa)) && (a[i] != 0.0))
          maxa = a[i];
      }
      for (i=0; i<n*n; i++) // normalizing
        a[i] /= maxa;
    }
  }
  en=(n*n-n)/2;
  //  initial values
  k=0;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      evec[k]=0.0;  k++;
    }
  }
  for (i=0;i<n;i++){
    evec[i*n+i]=1.0;
  }

  //  main iteration loop
  for (k=0;k<ni;k++){
    
    nn=0;
    for (i=0;i<n;i++){
      for (j=i+1;j<n;j++){
	if (fabs(a[i*n+j])<limit){
	  nn++;  continue;
	}
	else{
	  q=(a[j*n+j]-a[i*n+i])/2.0/a[i*n+j];
	  r=sqrt(1.0+q*q);
	  if (q>=0.0)  t=-q+r;
	  else         t=-q-r;
	  r=sqrt(1.0+t*t);
	  c=1.0/r;
	  s=t*c;
	  
	  //  A.J
	  ii=i;  jj=j;
	  for (l=0;l<n;l++){
	    ai=a[ii];
	    aj=a[jj];
	    a[ii]=ai*c-aj*s;
	    a[jj]=ai*s+aj*c;
	    ii+=n;  jj+=n;
	  }

	  //  J.A.J
	  ii=i*n;  jj=j*n;
	  for (l=0;l<n;l++){
	    ai=a[ii];
	    aj=a[jj];
	    a[ii]=ai*c-aj*s;
	    a[jj]=ai*s+aj*c;
	    ii++;  jj++;
	  }

	  //  Q=Q.J
	  ii=i;  jj=j;
	  for (l=0;l<n;l++){
	    ai=evec[ii];
	    aj=evec[jj];
	    evec[ii]=ai*c-aj*s;
	    evec[jj]=ai*s+aj*c;
	    ii+=n;  jj+=n;
	  }

	}
      }
    }
    if (nn==en){
      break;
    }
  }

  //  number of performed iterations
  ani=k;

  for (i=0;i<n;i++){
    eval[i]=a[i*n+i];
  }
  if ((normalize) && (maxa != 0.0))
  {
    for (i=0;i<n;i++)
      eval[i]*=maxa;
  }
  /*
  //  sorting eigenvalues and eigenvectors according to magnitudes of eigenvalues
  for (i=n-1;i>=0;i--){
    c=0.0;
    for (j=i;j>=0;j--){
      s=fabs(a[j*n+j]);
      if (s>c){
	c=s;  k=j;
      }
    }
    eval[i]=c;
    a[k*n+k]=a[i*n+i];

    ii=i;  jj=k;
    for (j=0;j<n;j++){
      t=evec[ii];
      evec[ii]=evec[jj];
      evec[jj]=t;
      ii+=n;  jj+=n;
    }
  }
  */

}



/**
   generalized Jacobi iteration method for problems A.x = w.B.x
      
   output
   @param x - array containing eigenvectors (stored in columns)
   @param w - array containing eigenvalues
   
   input
   @param a - matrix a (dense storage)
   @param b - matrix b (dense storage)
   @param n - number of rows (columns) of matrices a and b
   @param ni - maximum number of iterations
   @param thresholds - array containing thresholds
   @param nthr - number of thresholds
   @param zero - computer zero
   
   20.1.2002
*/
void gen_jacobi (double *a,double *b,double *x,double *w,long n,long ni,
		 double *thresholds,long nthr,double zero)
{
  long i,j,k,l,nt,numel,nze,stop;
  double discr,sol1,sol2,alpha,beta,limit,c0,c1,c2,denom;
  
  numel=(n*n-n)/2;
  nt=0;
  limit=thresholds[nt];
  stop=0;

  //  initial values
  for (i=0;i<n*n;i++){
    x[i]=0.0;
  }
  for (i=0;i<n;i++){
    x[i*n+i]=1.0;
  }
  
  //  iteration process
  for (i=0;i<ni;i++){
    fprintf (stdout,"\n\n JM iterace cislo %ld",i);


    nze=0;
    for (j=0;j<n;j++){
      for (k=j+1;k<n;k++){
	if (fabs(a[j*n+k])<limit && fabs(b[j*n+k])<limit)  nze++;
	else{
	  //  coefficients of quadratic equation
	  c2=a[j*n+k]*b[j*n+j]-a[j*n+j]*b[j*n+k];
	  c1=a[k*n+k]*b[j*n+j]-a[j*n+j]*b[k*n+k];
	  c0=a[k*n+k]*b[j*n+k]-a[j*n+k]*b[k*n+k];
	  
	  if (fabs(c2)<zero){
	    //  degenerated equation
	    if (fabs(c1)<zero){
	      fprintf (stderr,"\n\n non-solvable equation in function jacobi.\n");
	    }
	    alpha=(0.0-c0)/c1;
	  }
	  else{
	    discr=c1*c1-4.0*c2*c0;
	    if (discr<0.0){
	      fprintf (stderr,"\n\n negative discriminant in quadratic equation in function jacobi.\n");
	    }
	    discr=sqrt(discr);
	    sol1=(-1.0*c1+discr)/2.0/c2;
	    sol2=(-1.0*c1-discr)/2.0/c2;
	    //if (fabs(sol1)<fabs(sol2))  alpha=sol1;
	    //else                        alpha=sol2;
	    if (fabs(sol1)<fabs(sol2))  alpha=sol2;
	    else                        alpha=sol1;
	    
	  }
	  denom=a[k*n+k]+alpha*a[j*n+k];
	  if (fabs(denom)<zero){
	    fprintf (stderr,"\n\n zero denominator in expression for beta in function jacobi.\n");
	  }
	  beta=(0.0-a[j*n+k]-a[j*n+j]*alpha)/denom;
	  
	  //  AA=T.A and BB=T.B computation
	  for (l=0;l<n;l++){
	    sol1=a[j*n+l];
	    a[j*n+l]=a[j*n+l]+beta*a[k*n+l];
	    a[k*n+l]=a[k*n+l]+alpha*sol1;
	    sol1=b[j*n+l];
	    b[j*n+l]=b[j*n+l]+beta*b[k*n+l];
	    b[k*n+l]=b[k*n+l]+alpha*sol1;
	  }
	  
	  //  AA.T and BB.T computation
	  for (l=0;l<n;l++){
	    sol1=a[l*n+j];
	    a[l*n+j]=a[l*n+j]+beta*a[l*n+k];
	    a[l*n+k]=a[l*n+k]+alpha*sol1;
	    sol1=b[l*n+j];
	    b[l*n+j]=b[l*n+j]+beta*b[l*n+k];
	    b[l*n+k]=b[l*n+k]+alpha*sol1;
	  }
	  
	  //  T_i.T_{i+1} computation
	  for (l=0;l<n;l++){
	    sol1=x[l*n+j];
	    x[l*n+j]=x[l*n+j]+beta*x[l*n+k];
	    x[l*n+k]=x[l*n+k]+alpha*sol1;
	  }

	}
      }
    }
    if (nze==numel){
      fprintf (stdout,"\n\n zmena zavory");
      nt++;
      if (nt==nthr)  stop=1;
      else  limit=thresholds[nt];
    }
    
    if (stop==1)  break;
  }
  
  //  eigenvalues
  for (i=0;i<n;i++){
    w[i]=a[i*n+i]/b[i*n+i];
  }
  
  //  eigenvalues and eigenvectors sorting
  for (i=0;i<n;i++){
    sol1=fabs(w[i]);  k=i;
    for (j=i+1;j<n;j++){
      if (sol1>fabs(w[j])){
	sol1=fabs(w[j]);  k=j;
      }
    }
    sol2=w[i];
    w[i]=w[k];
    w[k]=sol2;
    for (j=0;j<n;j++){
      sol2=x[j*n+i];
      x[j*n+i]=x[j*n+k];
      x[j*n+k]=sol2;
    }
  }
}



void gen_jacobi2 (double *k,double *m,double *x,double *w,long n,long ni,
		 double *gate,long ng,double limit)
/*
  generalized Jacobi iteration method for problems A.x = w.B.x

   
   vystupy
   x - pole obsahujici vlastni vektory
   w - pole obsahujici vlastni cisla
   
   vstupy
   k - pole obsahujici matici K
   m - pole obsahujici matici M
   n - rozmer matic K a M
   ni - maximalni pocet iteraci
   gate - pole zavor
   ng - pocet zavor
   limit - konstanta pro testovani na nulu
   
   3.6.1998
*/
{
  long i,j,l,ii,jj,ll,rr,ss,zero,numel,nz;
  double kii,kij,kjj,mii,mij,mjj,s,k0,k1,k2,d,min,alpha,beta,err;
  
  numel=(n*n-n)/2;
  
  /**********************************/
  /*  nastaveni pocatecnich hodnot  */
  /**********************************/
  /*  nulovani pole pro vlastni vektory  */
  nullv (x,n*n);
  /*  pocatecni hodnoty  */
  for (i=0;i<n;i++){
    x[i*n+i]=1.0;
  }
  
  /*************/
  /*  iterace  */
  /*************/
  nz=0;  err=gate[nz];
  for (l=0;l<ni;l++){
    zero=0;
    for (i=0;i<n;i++){
      ii=i*n+i;  jj=ii;  kii=k[ii];  mii=m[ii];
      for (j=i+1;j<n;j++){
	ii++;  kij=k[ii];  mij=m[ii];
	jj+=n+1;  kjj=k[jj];  mjj=m[jj];
	if (fabs(kii)<limit || fabs(kjj)<limit ||
	    fabs(mii)<limit || fabs(mjj)<limit){
	  if (fabs(kij)<err && fabs(mij)<err){
	    zero++;  continue;
	  }
	}
	else{
	  if (fabs(kij*kij/kii/kjj)<err && fabs(mij*mij/mii/mjj)<err){
	    zero++;  continue;
	  }
	}
	k2=kij*mjj-kjj*mij;
	k1=kii*mjj-kjj*mii;
	k0=kii*mij-kij*mii;
	if (fabs(k2)<limit){
	  if (fabs(k1)<limit){
	    fprintf (stderr,"\n V Jacobiove metode je zdegenerovana kvadraticka rovnice");
	    fprintf (stderr,"\n Vypocet konci. \n\n");
	  }
	  else{
	    beta=-1.0*k0/k1;
	  }
	}
	else{
	  d=k1*k1-4.0*k2*k0;
	  if (d<-1.0*limit){
	    fprintf (stderr,"\n V Jacobiove metode je zaporny diskriminant.");
	    fprintf (stderr,"\n Vypocet konci. \n\n");
	  }
	  s=-k1+sqrt(d);
	  min=-k1-sqrt(d);
	  if (fabs(s)>fabs(min))  beta=s/2.0/k2;
	  else              beta=min/2.0/k2;
/*	  beta=(sqrt(d)-k1)/2.0/k2;*/
	}
	if (fabs(kii+beta*kij)<limit){
	  beta=-1.0*kij/kjj;
	  if (fabs(mii+beta*mij)<limit){
	    fprintf (stderr,"\n Zmatky v Jacobiove metode rotaci.");
	    fprintf (stderr,"\n Vypocet konci. \n\n");
	  }
	  else{
	    alpha=(mij+beta*mjj)/(mii+beta*mij)*(-1.0);
	  }
	}
	else{
	  alpha=(kij+beta*kjj)/(kii+beta*kij)*(-1.0);
	}

	/**********************/
	/*  provedeni rotace  */
	/**********************/
	k[i*n+i]+=2.0*beta*kij+beta*beta*kjj;
	m[i*n+i]+=2.0*beta*mij+beta*beta*mjj;
	k[j*n+j]+=2.0*alpha*kij+alpha*alpha*kii;
	m[j*n+j]+=2.0*alpha*mij+alpha*alpha*mii;
	k[i*n+j]=0.0;  m[i*n+j]=0.0;
	rr=i;  ss=j;
	for (ll=0;ll<i;ll++){
	  s=k[rr];
	  k[rr]+=beta*k[ss];
	  k[ss]+=alpha*s;
	  s=m[rr];
	  m[rr]+=beta*m[ss];
	  m[ss]+=alpha*s;
	  rr+=n;  ss+=n;
	}
	rr=i*n+i+1;  ss=(i+1)*n+j;
	for (ll=i+1;ll<j;ll++){
	  s=k[rr];
	  k[rr]+=beta*k[ss];
	  k[ss]+=alpha*s;
	  s=m[rr];
	  m[rr]+=beta*m[ss];
	  m[ss]+=alpha*s;
	  rr++;  ss+=n;
	}
	rr=i*n+j+1;  ss=j*n+j+1;
	for (ll=j+1;ll<n;ll++){
	  s=k[rr];
	  k[rr]+=beta*k[ss];
	  k[ss]+=alpha*s;
	  s=m[rr];
	  m[rr]+=beta*m[ss];
	  m[ss]+=alpha*s;
	  rr++;  ss++;
	}

	/*****************************/
	/*  vypocet vlastnich tvaru  */
	/*****************************/
	rr=i;  ss=j;
	for (ll=0;ll<n;ll++){
	  s=x[rr];
	  x[rr]+=beta*x[ss];
	  x[ss]+=alpha*s;
	  rr+=n;  ss+=n;
	}
      }
    }
    /*  posouzeni konvergence  */
    fprintf (stdout,"\n iterace %ld      %ld",l,zero);
    if (zero==numel){
      nz++;
      if (nz<ng){
	err=gate[nz];
	fprintf (stdout,"\n Bude se pracovat s %ld. zavorou.",nz+1);
	fprintf (stdout,"\n Velikost zavory je  %e",err);
      }
      else{
	break;
      }
    }
  }
  
  /*****************************/
  /*  vypocet vlastnich cisel  */
  /*****************************/
/*
  ii=0;
  for (i=0;i<n;i++){
    if (fabs(m[ii])<limit){
      w[i]=1.0e44;
    }
    else{
      w[i]=k[ii]/m[ii];
    }
    ii+=n+1;
  }
*/
  for (i=0;i<n;i++){
    w[i]=k[i*n+i]/m[i*n+i];
  }

  /******************************************************/
  /*  trideni pole vlastnich cisel a vlastnich vektoru  */
  /******************************************************/
  for (i=0;i<n;i++){
    min=1.0e50;
    for (j=i;j<n;j++){
      if (w[j]<min){
	min=w[j];  ii=j;
      }
    }
    if (ii!=i){
      s=w[i];
      w[i]=w[ii];
      w[ii]=s;
    }
  }
}
