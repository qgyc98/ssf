#include "csv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

compvect::compvect ()
{
  n=0;  nz=0;
  a=NULL;
  ind=NULL;
  
  //  computer zero
  limit=1.0e-15;
  //  indices are not ordered
  ordering=0;
}

compvect::compvect (long *ii,long nonz)
{
  long i;
  
  n=0;  nz=nonz;
  a=new double [nz];
  ind=new long [nz];
  
  for (i=0;i<nz;i++){
    ind[i]=ii[i];
    a[i]=0.0;
  }
  
  //  computer zero
  limit=1.0e-15;
  //  indices are not ordered
  ordering=0;
}

compvect::~compvect ()
{
  delete [] a;
  delete [] ind;
}

/**
   function copies current %vector to back up %vector bc
   
   @param bc - back up %vector
   
   JK, 27.2.2007
*/
void compvect::copy (compvect *bc)
{
  long i;
  
  bc->~compvect ();
  
  bc->n = n;
  bc->nz = nz;
  bc->ind = new long [nz];
  bc->a = new double [nz];
  
  for (i=0;i<nz;i++){
    bc->ind[i] = ind[i];
    bc->a[i] = a[i];
  }

}

/**
   function sets up compressed sparse %vector from dense %vector
   
   @param b - dense %vector
   @param nc - number of components of the dense %vector
   
   JK, 25.2.2007
*/
void compvect::setup_vector (double *b,long nc)
{
  long i;
  
  //  number of all components
  n=nc;
  
  //  number of nonzero components
  nz=0;
  for (i=0;i<n;i++){
    if (fabs(b[i])>limit)
      nz++;
  }
  
  if (a!=NULL){
    delete [] a;
  }
  a = new double [nz];
  if (ind!=NULL){
    delete [] ind;
  }
  ind = new long [nz];
  
  nz=0;
  for (i=0;i<n;i++){
    if (fabs(b[i])>limit){
      ind[nz]=i;
      a[nz]=b[i];
      nz++;
    }
  }
  
  //  indices are ordered increasingly
  ordering=1;
}

/**
   function tests ordering of indices of the %vector
   
   JK, 27.2.2007
*/
void compvect::test_ordering ()
{
  long i;
  
  for (i=1;i<nz;i++){
    if (ind[i]<=ind[i-1])
      break;
  }
  
  ordering=0;
}

/**
   function reorders components
   indices are reordered increasingly
   
   JK, 25.2.2007
*/
void compvect::reorder_components ()
{
  long i,j,k,l,min;
  double s;
  
  for (i=0;i<nz;i++){
    min=ind[i];
    for (j=i;j<nz;j++){
      if (min>ind[j]){
	min=ind[j];
	k=j;
      }
    }
    l=ind[i];
    ind[i]=k;
    ind[k]=l;
    
    s=a[i];
    a[i]=min;
    a[k]=s;
  }
  
  //  indices are ordered increasingly
  ordering=1;
}

/**
   function computes y = a*x + y

   @param x - compressed sparse %vector
   @param aa - scalar
   
   JK, 27.2.2007
*/
void compvect::axpy (compvect *x,double aa)
{
  long i,j,k;
  compvect bc;
  
  //  back up of the current vector
  copy (&bc);
  
  i=0;  j=0;  k=0;
  do{
    if (ind[i] == x->ind[j]){
      i++;  j++;  k++;
    }
    else{
      if (ind[i] < x->ind[j]){
	i++;
      }
      else{
	j++;
      }
    }
  }while(i<nz && j<x->nz);
  
  nz=k;
  delete [] ind;
  delete [] a;
  ind = new long [nz];
  a = new double [nz];
  
  i=0;  j=0;  k=0;
  do{
    if (bc.ind[i] == x->ind[j]){
      a[k]=bc.a[i]+x->a[j]*aa;
      ind[k]=bc.ind[i];
      i++;  j++;  k++;
    }
    else{
      if (bc.ind[i] < x->ind[j]){
	i++;
      }
      else{
	j++;
      }
    }
  }while(i<bc.nz && j<x->nz);
  
}

/**
   function computes y = a*x + y
   
   @param x - compressed sparse %vector
   @param aa - scalar
   
   JK, 27.2.2007
*/
void compvect::axpy_known (compvect *x,double aa)
{
  long i;
  
  for (i=0;i<nz;i++){
    a[i] += aa*x->a[i];
  }
}



/**
   function computes dot product of two compressed sparse vectors
   
   @param 

   JK, 25.2.2007
*/
double compvect::dotprod (compvect *cv)
{
  long j,k;
  long min1,max1,min2,max2,nc;
  long *ind2;
  double s,*a2;
  
  if (cv->ordering == 0){
    fprintf (stderr,"\n\n compressed sparse vector is not ordered increasingly in function comprow::ss (file %s, line %d)\n",__FILE__,__LINE__);
  }
  if (ordering == 0){
    fprintf (stderr,"\n\n compressed sparse vector is not ordered increasingly in function comprow::ss (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  ind2=cv->ind;
  a2=cv->a;
  
  //  number of nonzero components in the cv vector
  nc=cv->nz;
  //  minimum index in the cv vector
  min2=ind2[0];
  //  maximum index in the cv vector
  max2=ind2[nc-1];
  
  //  minimum index in the actual vector
  min1=ind[0];
  //  maximum index in the actual vector
  max1=ind[nz-1];

  if (max1<min2)
    return s=0.0;
  if (max2<min1)
    return s=0.0;
  
  j=0;  k=0;
  s=0.0;
  do{
    if (ind[j]==ind2[k]){
      s+=a[j]*a2[k];
      j++;  k++;
    }
    else{
      if (ind[j]<ind[k]){
	j++;
      }
      else{
	k++;
      }
    }
  }while (j<nz && k<nc);
  
  return s;
}
