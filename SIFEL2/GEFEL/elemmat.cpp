#include "elemmat.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

elemmat::elemmat (void)
{
  //  number of degrees of freedom
  n=0;
  //  number of elements
  ne=0;
  //  decomposition indicator
  decompid=0;
  //  total number of stored entries
  mem=0;
  //  maximum number of components in local vector
  maxarr=0;
  
  //  array containing element matrices
  a=NULL;
  //  array containing code numbers of finite elements
  acn=NULL;
  //  array containing number of matrix entries on elements
  neme = NULL;
  //  array containing number of DOFs on elements
  andofe = NULL;
}

elemmat::~elemmat (void)
{
  long i;
  
  for (i=0;i<ne;i++){
    delete [] a[i];
    delete [] acn[i];
  }
  delete [] a;
  delete [] acn;
  
  delete [] neme;
  delete [] andofe;
}

void elemmat::alloc (gtopology *top)
{
  long i,ndofe;
  
  //  number of finite elements
  ne=top->ne;
  //  array containing number of entries of one element
  neme = new long [ne];
  //  array containing number of DOF on one element
  andofe = new long [ne];
  
  mem=0;  maxarr=0;
  for (i=0;i<ne;i++){
    ndofe = top->give_ndofe (i);
    andofe[i] = ndofe;
    neme[i] = ndofe*ndofe;
    mem += neme[i];
    if (ndofe>maxarr)  maxarr=ndofe;
  }
  
  //  array containing matrices of finite elements
  a = new double* [ne];
  //  array containing code numbers of finite elements
  acn = new long* [ne];
  for (i=0;i<ne;i++){
    a[i] = new double [neme[i]];
    acn[i] = new long [andofe[i]];
  }
}

double** elemmat::status ()
{
  return a;
}

long elemmat::decomp ()
{
  return decompid;
}

void elemmat::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}

void elemmat::nullmat ()
{
  long i,j;
  for (i=0;i<ne;i++){
    for (j=0;j<neme[i];j++){
      a[i][j]=0.0;
    }
  }
}

void elemmat::localize (matrix &b,long *cn,long eid)
{
  long i,j,k;
  
  k=0;
  for (i=0;i<b.m;i++){
    acn[eid][i]=cn[i];
    for (j=0;j<b.m;j++){
      a[eid][k]=b[i][j];
      k++;
    }
  }
}

void elemmat::localized (double *b,long *cn,long eid,long m)
{
  long i,j,k;
  
  k=0;
  for (i=0;i<m;i++){
    acn[eid][i]=cn[i];
    for (j=0;j<m;j++){
      a[eid][k]=b[i*m+j];
      k++;
    }
  }
}

void elemmat::initiate (gtopology *top,long ndof,long /*mespr*/)
{
  if (status ()==NULL){
    alloc (top);
  }
  else{
    nullmat ();
  }
  
  n=ndof;

  //if (mespr==1)  fprintf (stdout,"\n total number of matrix entries   %ld",mem);
}


void elemmat::mxv_em (double *b,double *c)
{
  long i,j,ii,uj;
  double *u,*v;
  
  u = new double [maxarr];
  v = new double [maxarr];
  
  nullv (c,n);

  for (i=0;i<ne;i++){
    uj=andofe[i];
    
    nullv (u,maxarr);
    for (j=0;j<uj;j++){
      ii=acn[i][j];
      if (ii>0)
	u[j]=b[ii-1];
    }

    mxv (a[i],u,v,uj,uj);

    for (j=0;j<uj;j++){
      ii=acn[i][j];
      if (ii>0)
	c[ii-1]+=v[j];
    }
  }
  
  delete [] u;
  delete [] v;
}

void elemmat::addmat_em (double c,elemmat &dm)
{
  long i,j;
  for (i=0;i<ne;i++){
    for (j=0;j<neme[i];j++){
      a[i][j]+=c*dm.a[i][j];
    }
  }
}

void elemmat::scalmat_em (double c)
{
  long i,j;
  for (i=0;i<ne;i++){
    for (j=0;j<neme[i];j++){
      a[i][j]*=c;
    }
  }
}

void elemmat::printmat (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n");
  for (i=0;i<ne;i++){
    fprintf (out,"\n%ld",i);
    for (j=0;j<neme[i];j++){
      fprintf (out," %5.2f",a[i][j]);
    }
  }
}

/**
   function solves system of linear algebraic equations
   by conjugate gradient method,
   
   @param x - vector of unknowns
   @param y - vector of right hand side
   @param ni - maximum number of iterations
   @param err - required error (norm of residual vector)
   @param ani - number of performed iterations
   @param ares - attained error (norm of residual vector)
   @param zero - computer zero
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero vector
   iv=1 - initial vector is taken from x array
   
   26.6.2003
*/
void elemmat::cg (double *x,double *y,
		  long ni,double err,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double nom,denom,nory,alpha,beta;
  double *d,*r,*p;
  
  d = new double [n];
  r = new double [n];
  p = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  mxv_em (x,p);

  nory=0.0;  nom=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=p[i]-y[i];
    nom+=r[i]*r[i];
    d[i]=-1.0*r[i];
  }
  
  if (nory<zero){
    fprintf (stderr,"\n\n norm of right hand side in conjugate gradient method is smaller than %e",zero);
    fprintf (stderr,"\n see file %s, line %d.\n",__FILE__,__LINE__);
    ares=nory;  ani=0;
    return;
  }
  

  //  iteration loop
  for (i=0;i<ni;i++){
    
    
    //  new coefficient alpha
    //mxv_cr (d,p);
    mxv_em (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      fprintf (stdout,"\n there is zero denominator in alpha computation in conjugate gradient method (cr.cpp)\n");
      break;
    }
    
    alpha = nom/denom;
    
    //  new approximation of x and r
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]+=alpha*p[j];
    }
    
    denom=nom;
    
    nom = ss (r,r,n);
    
    fprintf (stdout,"\n iteration   %ld  norres/norrhs %e",i,nom/nory);
    
    if (nom/nory<err)  break;
    //if (fabs(nom)<limit)  break;
    
    
    beta = nom/denom;
    
    //  new vector of direction
    for (j=0;j<n;j++){
      d[j]=beta*d[j]-r[j];
    }
  }
  
  ani=i;  ares=nom;
  
  delete [] p;  delete [] r;  delete [] d;
}

