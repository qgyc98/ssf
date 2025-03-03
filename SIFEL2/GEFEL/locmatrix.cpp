#include "locmatrix.h"
#include <stdlib.h>
#include <stdio.h>

locmatrix::locmatrix ()
{
  nr=0;  nc=0;  nnc=0;
  
  nncr=NULL;  adr=NULL;
  lm=NULL;  ci=NULL;
}

locmatrix::~locmatrix ()
{
  delete [] nncr;
  delete [] adr;
  delete [] lm;
  delete [] ci;
}

/**
   function sets up basic variables
   
   @param nrows - number of rows
   @param ncolumns - number of columns
   
   JK, 23.12.2006
*/
void locmatrix::initiate_var (long nrows,long ncolumns)
{
  if (nrows!=nr){
    if (nncr!=NULL){
      delete [] nncr;
      nncr=NULL;
    }
    if (adr!=NULL){
      delete [] adr;
      adr=NULL;
    }
  }
  
  nr=nrows;
  nc=ncolumns;
  
  if (nncr==NULL){
    nncr = new long [nr];
  }
}

/**
   function computes addresses of the first nonzero components in particular rows
   
   JK, 23.12.2006
*/
void locmatrix::addresses ()
{
  long i;
  
  if (adr==NULL){
    adr = new long [nr+1];
  }

  adr[0]=0;
  for (i=1;i<nr+1;i++){
    adr[i]=adr[i-1]+nncr[i-1];
  }
  
  //  number of nonzero components
  nnc=adr[nr];
}

/**
   function allocates array ci
   
   JK, 23.12.2006
*/
void locmatrix::allocate_ci ()
{
  if (ci!=NULL){
    delete [] ci;
  }
  ci = new long [nnc];
}

/**
   function allocates array lm
   
   JK, 23.12.2006
*/
void locmatrix::allocate_lm ()
{
  if (lm!=NULL){
    delete [] lm;
  }
  lm = new double [nnc];
}

/**
   function initiates array of numbers of nonzero components in rows
   
   @param a - array containing numbers of nonzero components in rows
   a[i]=j - the i-th row contains j nonzero components

   JK, 23.12.2006
*/
void locmatrix::initiate_nncr (long *a)
{
  long i;
  for (i=0;i<nr;i++){
    nncr[i]=a[i];
  }
}

/**
   function initiates array column indices ci
   all components are initiated
   
   @param colind - array of column indices
   
   JK, 23.12.2006
*/
void locmatrix::initiate_ci (long *colind)
{
  long i;
  for (i=0;i<adr[nr];i++){
    ci[i]=colind[i];
  }
}

/**
   function initiates array column indices ci
   only column indices of one row are initiated
   
   @param id - number of actual row
   @param colind - array of column indices

   JK, 23.12.2006
*/
void locmatrix::initiate_ci (long id,long *colind)
{
  long i,j;
  
  if (id<0 || id>nr){
    fprintf (stderr,"\n wrong number of row is used in localization matrix (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  j=0;
  for (i=adr[id];i<adr[id+1];i++){
    ci[i]=colind[i];
    j++;
  }
}

/**
   function initiates array localization %matrix
   all components are initiated
   
   @param a - array contains entries of localization %matrix
   
   JK, 23.12.2006
*/
void locmatrix::initiate_lm (double *a)
{
  long i;
  for (i=0;i<adr[nr];i++){
    lm[i]=a[i];
  }
}

/**
   function initiates array localization %matrix
   
   @param id - number of actual row
   @param a - array contains entries of localization %matrix
   
   JK, 23.12.2006
*/
void locmatrix::initiate_lm (long id,double *a)
{
  long i,j;

  if (id<0 || id>nr){
    fprintf (stderr,"\n wrong number of row is used in localization matrix (file %s, line %d)\n",__FILE__,__LINE__);
  }

  j=0;
  for (i=adr[id];i<adr[id+1];i++){
    lm[i]=a[j];
    j++;
  }
}



/**
   function multiplies localization %matrix by the %vector
   L a = b
   
   JK, 23.12.2006
*/
void locmatrix::lmxv (double *a,double *b)
{
  long i,j,k;
  double s;

  for (i=0;i<nr;i++){
    s=0.0;
    for (j=adr[i];j<adr[i+1];j++){
      k=ci[j];
      s+=lm[j]*a[k];
    }
    b[i]=s;
  }
}

/**
   function multiplies transposed localization %matrix by the %vector
   L^T a = b
   
   JK, 23.12.2006
*/
void locmatrix::lmtxv (double *a,double *b)
{
  long i,j,k;
  double s;

  for (i=0;i<nc;i++){
    b[i]=0.0;
  }
  
  for (i=0;i<nr;i++){
    s=a[i];
    for (j=adr[i];j<adr[i+1];j++){
      k=ci[j];
      b[k]+=lm[j]*s;
    }
  }
}

/**
   function multiplies localization %matrix by the %vector
   L a = b

   localization %matrix contains entries equal to 0 or 1
   in this special case, function selects appropriate components
   and no real multiplication is required

   JK, 23.12.2006
*/
void locmatrix::lm01xv (double *a,double *b)
{
  long i,j,k;
  double s;

  for (i=0;i<nr;i++){
    s=0.0;
    for (j=adr[i];j<adr[i+1];j++){
      k=ci[j];
      s+=a[k];
    }
    b[i]=s;
  }
}

/**
   function multiplies transposed localization %matrix by the %vector
   L^T a = b

   localization %matrix contains entries equal to 0 or 1
   in this special case, function selects appropriate components
   and no real multiplication is required

   JK, 23.12.2006
*/
void locmatrix::lmt01xv (double *a,double *b)
{
  long i,j,k;
  double s;
  
  for (i=0;i<nc;i++){
    b[i]=0.0;
  }
  
  for (i=0;i<nr;i++){
    s=a[i];
    for (j=adr[i];j<adr[i+1];j++){
      k=ci[j];
      b[k]+=s;
    }
  }
}



/**
   function multiplies localization %matrix by the %matrix
   L A = B
   %matrix A has m rows and n columns

   localization %matrix contains entries equal to 0 or 1
   in this special case, function selects appropriate components
   and no real multiplication is required

   JK, 23.12.2006
*/
void locmatrix::lm01xm (gmatrix &a,gmatrix &b)
{
  long i,j,k,l;
  double s;
  
  for (i=0;i<nr;i++){
    for (j=0;j<a.n;j++){
      s=0.0;
      for (k=adr[i];k<adr[i+1];k++){
	l=ci[k];
	s+=a.give_entry (i,l);
      }
    }
    b.add_entry (s,i,j);
  }
}


/**
   function reduces global %matrix to coarse grid %matrix
   localization %matrix must be Boolean and only one unit entry has to be in each row and column
   number of rows of localization %matrix must be less than number of columns
   
   lm(m,n).a(n,n).lm^T(n,m)=b(m,m)
   
   function stores selected %matrix entries in the skyline storage format
   
   @param a - global %matrix
   @param b - coarse grid %matrix

   JK, 27.12.2006
*/
void locmatrix::lmxmxlmt01 (gmatrix &a,gmatrix &b)
{
  long i,j,ii,jj,ncontr;
  long *rind,*cind,*aux;
  double s;
  double *entr;
  
  //  evaluation of number of contributions
  ncontr=0;
  for (i=0;i<nr;i++){
    ii=ci[i];
    for (j=0;j<nr;j++){
      jj=ci[j];
      s=a.give_entry (ii,jj);
      if (fabs(s)>threshold)
	ncontr++;
    }
  }
  
  //  row and column indices, matrix entries
  rind = new long [ncontr];
  cind = new long [ncontr];
  entr = new double [ncontr];
  
  ncontr=0;
  for (i=0;i<nr;i++){
    ii=ci[i];
    for (j=0;j<nr;j++){
      jj=ci[j];
      s=a.give_entry (ii,jj);
      if (fabs(s)>threshold){
	rind[ncontr]=ii;
	cind[ncontr]=jj;
	entr[ncontr]=s;
	ncontr++;
      }
    }
  }
  
  
  //  transformation of storage system from the compresed rows to skyline
  aux = new long [nr+1];
  for (i=0;i<nr+1;i++){
    aux[i]=-1;
  }
  //  lengths of columns in skyline
  for (i=0;i<ncontr;i++){
    if (rind[i]>cind[i]){
      //  lower part of the matrix
      continue;
    }
    if (cind[i]-rind[i]+1>aux[cind[i]])
      aux[cind[i]]=cind[i]-rind[i]+1;
  }
  //  addresses of diagonal entries in skyline
  ii=aux[0];
  aux[0]=0;
  for (i=1;i<nr+1;i++){
    jj=aux[i];
    aux[i]=aux[i-1]+ii;
    ii=jj;
  }
  
  //  allocation of array with addresses of diagonal entries
  b.sky->allocadr (nr);
  //
  for (i=0;i<nr+1;i++){
    b.sky->adr[i]=aux[i];
  }
  //  number of stored matrix entries
  b.sky->neglobmat ();
  //  allocation of array for matrix entries
  b.sky->allocglomat ();
  
  //  matrix assembling
  for (i=0;i<ncontr;i++){
    if (rind[i]>cind[i]){
      //  lower part of the matrix
      continue;
    }
    j=aux[cind[i]];
    b.sky->a[j+cind[i]-rind[i]]=entr[i];
  }
  
  //  array cleaning
  delete [] aux;
  delete [] rind;
  delete [] cind;
  delete [] entr;
}


/**
   function reduces global %matrix to coarse grid %matrix
   number of rows of localization %matrix must be less than number of columns
   
   lm(m,n).a(n,n).lm^T(n,m)=b(m,m)
   
   function stores selected %matrix entries in the skyline storage format
   
   @param a - global %matrix
   @param b - coarse grid %matrix

   JK, 6.1.2007
*/
void locmatrix::lmxmxlmt (gmatrix &a,gmatrix &b)
{
  long i,j,k,l,n,cind;
  long *aux,*nncaux,*adraux,**ciaux;
  double s;

  n=a.n;
  
  aux = new long [n];
  
  //  auxiliary array, number of nonzero entries in matrix LM x M
  nncaux = new long [nr];
  
  //  auxiliary array, addresses of first entries
  adraux = new long [nr+1];
  
  //  auxiliary array, column indices
  ciaux = new long* [nr];
  
  //  determination of nonzero entries in matrix LM x M
  adraux[0]=0;
  for (i=0;i<nr;i++){
    nncaux[i]=0;
    for (j=0;j<n;j++){
      s=0.0;  l=adr[i];  aux[j]=-1;
      for (k=0;k<nncr[i];k++){
	cind=ci[l];
	s+=lm[l]*a.give_entry (i,cind);
	l++;
      }
      if (fabs(s)>threshold){
	nncaux[i]++;
	aux[j]=j;
      }
    }
    adr[i+1]=adr[i]+nncaux[i];
    ciaux[i] = new long [nncaux[i]];
    k=0;
    for (j=0;j<n;j++){
      if (aux[j]>-1){
	ciaux[i][k]=aux[j];
	k++;
      }
    }
  }
  
  
}
