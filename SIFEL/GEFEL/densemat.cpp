#include "densemat.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "precond.h"

densemat::densemat (void)
{
  //  the number of rows in the matrix (it is equal to the number of unknowns)
  n=0;
  //  array containing matrix entries
  a=NULL;
  //  indicator of matrix factorization
  decompid=0;
}

densemat::~densemat (void)
{
  delete [] a;
}

/**
   function allocates array containing %matrix
   
   @param m - number of rows/columns
   
   JK
*/
void densemat::alloc (long m)
{
  n=m;
  negm=m*m;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));
}

void densemat::dealloc (void)
{
  delete [] a;
  a = NULL;
}

void densemat::copy (densemat *dm)
{
  long i,j;

  if (a != NULL)
    delete [] a;  

  n=dm->n;
  negm=n*n;
  a = new double [negm];
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]=dm->a[i*n+j];
    }
  }
}

/**
   function creates a copy of the %matrix
   
   @param dm - dense %matrix which will be set up by data
   
   JK, 15.3.2007
*/
void densemat::copy_dm (densemat &dm)
{
  long i,j;
  
  //  the number of rows/columns
  dm.n=n;
  //  the number of matrix entries
  dm.negm=n*n;
  
  if (dm.a != NULL)
    delete [] dm.a;  
  dm.a = new double [negm];
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      dm.a[i*n+j]=a[i*n+j];
    }
  }
}


/**
   function returns status of array a
   
   JK
*/
double* densemat::status ()
{
  return a;
}

/**
   function returns indicator of decomposition (factorization)
   
   JK
*/
long densemat::decomp ()
{
  return decompid;
}

/**
   function changes indicator of decomposition (factorization)
   
   JK
*/
void densemat::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}

void densemat::setfact ()
{
  decompid=1;
}
void densemat::setnotfact ()
{
  decompid=0;
}

/**
   function fills array by zeros
   
   JK
*/
void densemat::nullmat ()
{
  long i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]=0.0;
    }
  }
}

/**
   function localizes local %matrix b into global %matrix
   
   @param b - local %matrix
   @param cn - array containing code numbers
   
   JK
*/
void densemat::localize (matrix &b,long *cn)
{
  long i,j,ii,jj;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      a[ii*n+jj]+=b[i][j];
    }
  }
}

/**
   function localizes local %matrix stored in array into global %matrix
   
   @param b - array containing local %matrix
   @param cn - array containing code numbers
   @param m - number of rows/columns of the %matrix b
   
   JK
*/
void densemat::localized (double *b,long *cn,long m)
{
  long i,j,ii,jj;
  
  for (i=0;i<m;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      a[ii*n+jj]+=b[i*m+j];
    }
  }
}

/**
   function localizes general rectangular (non-square) %matrix into global %matrix
   
   @param b - local %matrix
   @param rcn - row code numbers
   @param ccn - column code numbers
   
   JK
*/
void densemat::glocalize (matrix &b,long *rcn,long *ccn)
{
  long i,j,ii,jj;
  
  for (i=0;i<b.m;i++){
    ii=rcn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<b.n;j++){
      jj=ccn[j]-1;
      if (jj<0)  continue;
      a[ii*n+jj]+=b[i][j];
    }
  }
}

/**
   function localizes contributions from Lagrange multipliers to the %skyline storage
   
   @param nm - number of Lagrange multipliers
   @param ncn1 - nodal code numbers of the first node (on one side of interface)
   @param ncn2 - nodal code numbers of the second node (on the other side of interface)
   @param mcn - code numbers of Lagrange multipliers defined between the previous nodes

   JK, 18.8.2008
*/
void densemat::mult_localize (long nm,long *ncn1,long *ncn2,long *mcn)
{
  long i,cn,m;
  
  for (i=0;i<nm;i++){
    cn=mcn[i]-1;
    if (cn>-1){
      m=ncn1[i]-1;
      if (m>-1){
	a[m*n+cn]=-1.0;
	a[cn*n+m]=-1.0;
      }
      m=ncn2[i]-1;
      if (m>-1){
	a[m*n+cn]=1.0;
	a[cn*n+m]=1.0;
      }
    }
  }
}

/**
   function initiates densemat storage
   
   @param ndof - number of rows/columns of the %matrix
   @param mespr - message printing indicator
   
   JK
*/
void densemat::initiate (long ndof,long mespr)
{
  if (status ()==NULL){
    alloc (ndof);
  }
  else{
    nullmat ();
  }
  
  if (mespr==1)  fprintf (stdout,"\n number of matrix entries   %ld",negm);
}

/**
   function multiplies %matrix stored as dense %matrix by %vector b,
   resulting %vector is c
   
   @param b - array containing %vector b
   @param c - array containing %vector c
   
   JK
*/
void densemat::mxv_dm (double *b,double *c)
{
  long i,j,aca;
  double s;

  aca=0;
  for (i=0;i<n;i++){
    s=0.0;
    for (j=0;j<n;j++){
      s+=a[aca]*b[j];
      aca++;
    }
    c[i]=s;
  }
}

/**
   function adds premultiplied %matrix stored in dm by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param dm - another dense %matrix storage
   
   JK
*/
void densemat::addmat_dm (double c,densemat &dm)
{
  long i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]+=c*dm.a[i*n+j];
    }
  }
}

/**
   function multiplies components of the %matrix by coefficient c
   
   @param c - multiplicative coefficient
   
   JK
*/
void densemat::scalmat_dm (double c)
{
  long i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]*=c;
    }
  }
}
/**
   function multiplies matrix a stored in dense format with matrix b stored in dense format
   ma - matrix a
   mb - matrix b
   nm - number of colums of matrix a and number of rows of matrix b
   mc - result of computation
   
   JB
   
*/
void densemat::maxmb(long nm, double *ma, double *mb, double *mc)
{
  long i,j,k;
  double s;
  
  for(i = 0; i < nm; i++){
    for(j = 0; j < nm; j++){
      s=0.0;
      for(k = 0; k< nm; k++){
	s+=ma[i*nm+k]*mb[k*nm+j];
      }
      mc[i*n+j]=s;
    }
  }
  
  
}


/**
   function solves equation system by Gauss elimination method
   %matrix is stored as dense %matrix
   x and y are also dense matrices, stored by rows
   right hand side vectors and vectors of unknowns are columns!

   @param x - array containing vectors of unknowns
   @param y - array containing right hand sides
   @param m - number of right hand sides
   @param zero - computer zero
   @param pivot - type of pivoting
   
   pivot=1 - pivot is searched only when there is zero on the diagonal
   pivot=2 - pivot is searched every time
   
   JK, 23.7.2001
*/
void densemat::gemp (double *x,double *y,long m,double zero,long pivot)
{
  long i,j,k,ii,jj,kk,ir,ic,li,ui,*order;
  double f,s;
  
  order = new long [n];
  
    
  //  initialization of order array
  for (i=0;i<n;i++){
    order[i]=i;
  }
  
  for (k=0;k<n-1;k++){
    ir=k;  ic=k;

    if (pivot==1){
      //  pivot is searched only when A(k,k)=0
      if (fabs(a[k*n+k])<zero){
        s=0.0;  ii=k*n+k;
        for (i=k;i<n;i++){
          for (j=k;j<n;j++){
            f=fabs(a[ii]);
            if (f>s){
              s=f;  ir=i;  ic=j;
            }
            ii++;
          }
          ii+=k;
        }
        if (s<zero){
          print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
          break;
        }
      }
    }
    
    if (pivot==2){
      //  pivot is searched in every step
      s=0.0;  ii=k*n+k;
      for (i=k;i<n;i++){
        for (j=k;j<n;j++){
          f=fabs(a[ii]);
          if (f>s){
            s=f;  ir=i;  ic=j;
          }
          ii++;
        }
        ii+=k;
      }
      if (s<zero){
        print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
        break;
      }
    }
    
    //  rows exchange
    if (ir!=k){
      li=k*n+k;  ui=(k+1)*n;  j=ir*n+k;
      for (i=li;i<ui;i++){
        s=a[i];
        a[i]=a[j];
        a[j]=s;
        j++;
      }

      li=k*m;  ui=li+m;  j=ir*m;
      for (i=li;i<ui;i++){
        s=y[i];
        y[i]=y[j];
        y[j]=s;
        j++;
      }
    }

    //  column exchange
    if (ic!=k){
      i=order[k];
      order[k]=order[ic];
      order[ic]=i;
      
      j=k;  ii=ic;
      for (i=0;i<n;i++){
        s=a[j];
        a[j]=a[ii];
        a[ii]=s;
        j+=n;  ii+=n;
      }                         
    }
    
    
    //  elimination
    
    for (i=k+1;i<n;i++){
      ii=i*n+k;  kk=k*n+k;
      s=a[ii]/a[kk];

      //  modification of matrix A
      for (j=k;j<n;j++){
        a[ii]-=s*a[kk];
        ii++;  kk++;
      }

      //  modification of right hand sides
      ii=i*m;  kk=k*m;
      for (j=0;j<m;j++){
        y[ii]-=s*y[kk];
        ii++;  kk++;
      }
    }
    
    
  }
  
  if (fabs(a[(n-1)*n+n-1])<zero){
    print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
  }
  

  //  back substitution
  
  for (k=n-1;k>-1;k--){
    f=a[k*n+k];  kk=k*m;
    for (i=0;i<m;i++){
      s=0.0;  ii=k*n+n-1;  jj=(n-1)*m+i;
      for (j=n-1;j>k;j--){
        s+=a[ii]*x[jj];
        ii--;  jj-=m;
      }
      x[kk]=(y[kk]-s)/f;
      kk++;
    }
  }
  
  //  reordering of unknowns into original positions
  for (i=0;i<n;i++){
    if (order[i]!=i){
      for (j=i;j<n;j++){
        if (order[j]==i){
          jj=j;  break;
        }
      }
      
      ii=order[i];
      order[i]=order[jj];
      order[jj]=ii;
      
      ii=i*m;  kk=jj*m;
      for (j=0;j<m;j++){
        s=x[ii];
        x[ii]=x[kk];
        x[kk]=s;
        ii++;  kk++;
      }
    }
  }
  
  delete [] order;
}


/**
   function reduces equation system by Gauss elimination method
   %matrix is stored as dense %matrix
   x and y are also dense matrices, stored by rows
   right hand side vectors and vectors of unknowns are columns!
   
   it is condensation strategy

   @param b - array containing condensed %matrix
   @param x - array containing left hand side
   @param y - array containing right hand side
   @param m - number of remaining unknowns
   @param zero - computer zero
   @param tc - type of computation
   
   tc=1 - function decomposes %matrix and creates condensed %matrix
   tc=2 - function computes backward substitution
   
   JK, 2.6.2002
*/
void densemat::gempkon (double *b,double *c,double *x,double *y,long m,double zero,long tc)
{
  long i,j,k,ii,kk;
  double s;

  if (tc==1){
    for (k=0;k<n-m;k++){
      
      if (fabs(a[k*n+k])<zero)  fprintf (stderr,"\n");
      
      for (i=k+1;i<n;i++){
        ii=i*n+k;  kk=k*n+k;
        s=a[ii]/a[kk];
        
        //  modification of matrix A
        for (j=k;j<n;j++){
          a[ii]-=s*a[kk];
          ii++;  kk++;
        }
        
        //  modification of right hand side
        y[i]-=s*y[k];
        
      }
      
    }
    
    //  reduced matrix assembling
    ii=n-m;
    for (i=0;i<m;i++){
      kk=n-m;
      for (j=0;j<m;j++){
        b[i*m+j]=a[ii*n+kk];
        kk++;
      }
      ii++;
    }
    
    //  reduced vector assembling
    ii=n-m;
    for (i=0;i<m;i++){
      c[i]=y[ii];
      ii++;
    }
    
    
  }
  
  if (tc==2){
    //  interface components are localized to the vector of subdomain
    j=0;
    for (k=n-m;k<n;k++){
      x[k]=c[j];  j++;
    }
    
    //  back substitution
    for (k=n-m-1;k>-1;k--){
      s=0.0;  ii=k*n+k+1;
      for (j=k+1;j<n;j++){
        s+=a[ii]*x[j];
        ii++;
      }
      x[k]=(y[k]-s)/a[k*n+k];
    }
  }
  
}

/**
   function decomposes %matrix A into L.U form
   
   @param x - left hand side
   @param y - right hand side
   @param zero - computer zero
   @param tc - type of computation

   tc = 1 - decomposition and back substitution
   tc = 2 - decomposition of %matrix
   tc = 3 - back substitution
   
   JK, 6.1.2000
*/
void densemat::lu (double *x,double *y,double zero,long tc)
{
  long i,j,k;
  double s;
  
  if (tc==1 || tc==2){
    if (decompid==0){
      for (i=0;i<n;i++){
	for (j=0;j<i;j++){
	  s=0.0;
	  for (k=0;k<j;k++){
	    s+=a[i*n+k]*a[k*n+j];
	  }
	  a[i*n+j]-=s;
	}
	
	s=0.0;
	for (k=0;k<i;k++){
	  s+=a[i*n+k]*a[k*n+i];
	}
	a[i*n+i]-=s;
	if (fabs(a[i*n+i])<zero){
          print_err("zero diagonal entry in lu factorization", __FILE__, __LINE__, __func__);
	}
	
	for (j=i+1;j<n;j++){
	  s=0.0;
	  for (k=0;k<i;k++){
	    s+=a[i*n+k]*a[k*n+j];
	  }
	  a[i*n+j]=(a[i*n+j]-s)/a[i*n+i];
	}
      }
      decompid=1;
    }
  }
  
  if (tc==1 || tc==3){
    if (decompid==0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      for (i=0;i<n;i++){
	s=0.0;
	for (j=0;j<i;j++){
	  s+=a[i*n+j]*y[j];
	}
	y[i]=(y[i]-s)/a[i*n+i];
      }
      for (i=n-1;i>-1;i--){
	s=0.0;
	for (j=i+1;j<n;j++){
	  s+=a[i*n+j]*x[j];
	}
	x[i]=y[i]-s;
      }
    }
  }
}

/**
   function computes Cholesky decomposition of %matrix and solves
   system of equations

   notice: this function can be used only for symmetric positive definite matrices
   function was implemented for testing other Cholesky decompositions performed
   on more complicated schemes of %matrix storage
   
   @param x - left hand side
   @param y - right hand side
   @param zero - computer zero
   @param tc - type of computation

   tc = 1 - decomposition and back substitution
   tc = 2 - decomposition of %matrix
   tc = 3 - back substitution
   
   JK, 17.12.2005
*/
void densemat::ll (double *x,double *y,double /*zero*/,long tc)
{
  long i,j,k;
  double s;
  double *z;
  
  z = new double [n];

  if (tc==1 || tc==2){
    if (decompid==0){
      for (i=0;i<n;i++){
	for (j=0;j<i;j++){
	  s=0.0;
	  for (k=0;k<j;k++){
	    s+=a[i*n+k]*a[j*n+k];
	  }
	  a[i*n+j]-=s;
	  a[i*n+j]/=a[j*n+j];
	}
	
	s=0.0;
	for (k=0;k<i;k++){
	  s+=a[i*n+k]*a[i*n+k];
	}
	a[i*n+i]-=s;
	if (a[i*n+i]<0.0){
          print_err("negative diagonal entry in Cholesky factorization", __FILE__, __LINE__, __func__);
	}
	a[i*n+i]=sqrt(a[i*n+i]);
	
      }
      
      decompid=1;
    }
  }
  
  if (tc==1 || tc==3){
    if (decompid==0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      for (i=0;i<n;i++){
	s=0.0;
	for (k=0;k<i;k++){
	  s+=a[i*n+k]*z[k];
	}
	z[i]=(y[i]-s)/a[i*n+i];
      }
      
      for (i=n-1;i>-1;i--){
	s=0.0;
	for (k=i+1;k<n;k++){
	  s+=a[k*n+i]*x[k];
	}
	x[i]=(z[i]-s)/a[i*n+i];
      }
    }
  }
  
  delete [] z;
}

/**
   function computes incomplete Cholesky decomposition of %matrix and solves
   system of equations

   notice: this function can be used only for symmetric positive definite matrices
   function was implemented for testing other Cholesky decompositions performed
   on more complicated schemes of %matrix storage
   
   @param x - left hand side
   @param y - right hand side
   @param zero - computer zero
   @param limit - threshold for fill-in
   @param tc - type of computation

   tc = 1 - decomposition and back substitution
   tc = 2 - decomposition of %matrix
   tc = 3 - back substitution
   
   JK, 17.12.2005
*/
void densemat::ill (double *x,double *y,double /*zero*/,double limit,long tc)
{
  long i,j,k;
  double s;
  double *z;
  
  z = new double [n];

  if (tc==1 || tc==2){
    if (decompid==0){
      for (i=0;i<n;i++){
	for (j=0;j<i;j++){
	  if (fabs(a[i*n+j])>limit){
	    s=0.0;
	    for (k=0;k<j;k++){
	      s+=a[i*n+k]*a[j*n+k];
	    }
	    a[i*n+j]-=s;
	    a[i*n+j]/=a[j*n+j];
	  }
	}
	
	s=0.0;
	for (k=0;k<i;k++){
	  s+=a[i*n+k]*a[i*n+k];
	}
	a[i*n+i]-=s;
	if (a[i*n+i]<0.0){
          print_err("negative diagonal entry in Cholesky decomposition", __FILE__, __LINE__, __func__);
	}
	a[i*n+i]=sqrt(a[i*n+i]);
	
      }
      decompid=1;
    }
  }
  
  if (tc==1 || tc==3){
    if (decompid==0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      for (i=0;i<n;i++){
	s=0.0;
	for (k=0;k<i;k++){
	  s+=a[i*n+k]*z[k];
	}
	z[i]=(y[i]-s)/a[i*n+i];
      }
      
      for (i=n-1;i>-1;i--){
	s=0.0;
	for (k=i+1;k<n;k++){
	  s+=a[k*n+i]*x[k];
	}
	x[i]=(z[i]-s)/a[i*n+i];
      }
    }
  }
  
  delete [] z;
}



/**
   function computes kernel of %matrix
   
   @param r - array containing base vectors of kernel
   @param dim - dimension of kernel
   @param se - array containing numbers (indices) of singular equations
   @param ense - estimated number of singular equations
   @param limit - defines zero on diagonal

   21.4.2003, JK
*/
void densemat::ker (double *r,long &dim,long *se,long ense,double limit)
{
  long i,j,k,stop;
  double s,*b;
  
  //  array for backup
  b = new double [ense*n];
  
  //  elimination of the matrix
  dim=0;
  for (i=0;i<n;i++){
    if (fabs(a[i*n+i])<limit){
      se[dim]=i;

      for (j=0;j<n;j++){
	a[i*n+j]=0.0;
      }

      a[i*n+i]=1.0;
      k=dim*n;
      for (j=0;j<n;j++){
	b[k]=a[j*n+i];  k++;
	a[j*n+i]=0.0;
      }
      a[i*n+i]=1.0;

      dim++;
      if (dim>ense){
        print_err("wrong number of singular equations, there are more singular equations than has been estimated", __FILE__, __LINE__, __func__);
	abort ();
      }
    }
    else{
      for (j=i+1;j<n;j++){
	s=a[j*n+i]/a[i*n+i];
	for (k=i;k<n;k++){
	  a[j*n+k]-=s*a[i*n+k];
	}
      }
    }
  }
  
  //  computation of base vectors of the kernel
  //  return of original entries
  k=0;
  for (i=0;i<dim;i++){
    for (j=0;j<n;j++){
      a[j*n+se[i]]=b[k];  k++;
    }
  }
  //  backward substitution
  nullv (r,dim*n);
  for (i=0;i<dim;i++){
    r[i*n+se[i]]=1.0;
    for (j=n-1;j>-1;j--){
      stop=0;
      for (k=0;k<dim;k++){
	if (j==se[k]){
	  stop=1;  break;
	}
      }
      if (stop==1)  continue;
      
      s=0.0;
      for (k=j+1;k<n;k++){
	s-=a[j*n+k]*r[i*n+k];
      }
      
      r[i*n+j]=s/a[j*n+j];
    }
  }
  delete [] b;
  
}


/**
   function solves system of linear algebraic equations
   by conjugate gradient method, %matrix is stored as dense
   notice: this function can be used only for symmetric positive definite matrices
   function was implemented for testing other conjugate gradient methods performed
   on more complicated schemes of %matrix storage

   @param x - %vector of unknowns
   @param y - %vector of right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   @param zero - computer zero
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero vector
   iv=1 - initial vector is taken from x array
   
   JK, 17.12.2005
*/
void densemat::cg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv)
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
  mxv_dm (x,p);

  nory=0.0;  nom=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=p[i]-y[i];
    nom+=r[i]*r[i];
    d[i]=-1.0*r[i];
  }
  
  if (nory<zero){
    print_err("norm of right hand side in conjugate gradient method is smaller than %e", __FILE__, __LINE__, __func__,zero);
    ares=nory;  ani=0;
    for (i=0;i<n;i++)
      x[i]=0.0;
    return;
  }
  

  //  iteration loop
  for (i=0;i<ni;i++){
    
    
    //  new coefficient alpha
    mxv_dm (d,p);
    //mv_cr15 (d,p);
    
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
    
    if (nom/nory<res)  break;
    
    beta = nom/denom;
    
    //  new vector of direction
    for (j=0;j<n;j++){
      d[j]=beta*d[j]-r[j];
    }
  }
  
  ani=i;  ares=nom;
  
  delete [] p;  delete [] r;  delete [] d;

}

/**
   @param x - %vector of unknowns
   @param y - %vector of right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residula %vector
   @param ani - number of computed iterations
   @param ares - reached residual norm
   @param zero - computer zero
   @param iv - initial values indicator
   @param tprec - type of preconditioner
   @param par - preconditioning parameter
   
   iv=0 - initial vector is assumed as zero vector
   iv=1 - initial vector is taken from x array
   
   tprec=1 - diagonal (Jacobi) preconditioner
   tprec=10 - incomplete decomposition
   
   par: relaxation parameter - for SSOR preconditioner
        weight coefficient - for incomplete decomposition
	
   JK, 17.12.2005
*/
void densemat::cg_prec (densemat &dm,double *x,double *y,long ni,double res,long &ani,double &ares,
			double zero,long iv,long tprec,double /*par*/)
{
  long i,j;
  double nom,denom,nory,alpha,beta;
  double *d,*r,*p,*h;
  
  d = new double [n];
  r = new double [n];
  p = new double [n];
  h = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  mxv_dm (x,p);
  nory=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=p[i]-y[i];
  }
  
  if (nory<zero){
    print_err("right hand side vector has zero norm", __FILE__, __LINE__, __func__);
    ares=nory;  ani=0;
  }
  
  switch (tprec){
    //case 1:{   prec_diag_scr (h,r);  break;  }
    //case 5:{   prec_ssor_scr (h,r,par);  break; }
    case 10:{
      //prec_ildl_scr (h,r);
      dm.ill (h,r,zero,zero,3);
      break;  
    }
    //case sparseindec:{
    //sdirect->Solve (h,r);
    //break;
    //}
  default:{
    print_err("wrong preconditioner type is required", __FILE__, __LINE__, __func__);
  }
  }

  nom = ss (r,h,n);
  
  for (i=0;i<n;i++){
    d[i]=0.0-h[i];
  }
  
  //  iteration loop
  for (i=0;i<ni;i++){

    //  matrix-vector multiplication
    mxv_dm (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      print_err("denominator in alpha expression is equal to zero", __FILE__, __LINE__, __func__);
      ares=nom;  ani=i;
      return;
    }
    
    //  new alpha
    alpha = nom/denom;
    
    //  new vectors
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]+=alpha*p[j];
    }
    
    
    //  preconditioning
    switch (tprec){
      //case 1:{   prec_diag_scr (h,r);  break;  }
      //case 5:{   prec_ssor_scr (h,r,par);  break; }
      case 10:{
	//prec_ildl_scr (h,r);  break;  }
	dm.ill (h,r,zero,zero,3);
	break;
      }
	//case sparseindec:{
	//sdirect->Solve (h,r);
	//break;
	//}
    default:{
      print_err("wrong preconditioner type is required", __FILE__, __LINE__, __func__);
    }
    }
    
    denom=nom;
    nom = ss (r,h,n);
    beta = nom/denom;
    
    fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/nory);
    
    if (fabs(nom/nory)<res){
      ani=i;  ares=nom/nory;
      break;
    }
    
    //  new direction vector
    for (j=0;j<n;j++){
      d[j]=beta*d[j]-h[j];
    }
  }
  
  delete [] h;
  delete [] p;
  delete [] r;
  delete [] d;
}








void densemat::cg_prec_new (precond &pr,double *x,double *y,long ni,double res,long &ani,double &ares,
			    double zero,long iv)
{
  long i,j;
  double nom,denom,nory,alpha,beta;
  double *d,*r,*p,*h;

  d = new double [n];
  r = new double [n];
  p = new double [n];
  h = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  
  //  A x_0 = p
  mxv_dm (x,p);
  
  //  norm of the right hand side
  nory=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=y[i]-p[i];
  }
  
  if (nory<zero){
    print_err("right hand side vector has zero norm", __FILE__, __LINE__, __func__);
    ares=nory;  ani=0;
  }
  
  //  h = M^-1 r
  pr.preconditioning (r,h,y);
  
  nom = ss (r,h,n);
  
  //  direction vector
  for (i=0;i<n;i++){
    d[i]=h[i];
  }
  
  
  //  iteration loop
  for (i=0;i<ni;i++){
    
    //  matrix-vector multiplication
    mxv_dm (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      print_err("denominator in alpha expression is equal to zero", __FILE__, __LINE__, __func__);
      ares=nom;  ani=i;
      return;
    }
    
    //  new alpha
    alpha = nom/denom;
    
    //  new vectors
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
    }
    
    //  h = M^-1 r
    pr.preconditioning (r,h,y);
    
    
    denom=nom;
    nom = ss (r,h,n);
    beta = nom/denom;
    
    fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/nory);
    
    if (fabs(nom/nory)<res){
      ani=i;  ares=nom/nory;
      break;
    }
    
    //  new direction vector
    for (j=0;j<n;j++){
      d[j]=beta*d[j]+h[j];
    }
  }
  
  delete [] h;
  delete [] p;
  delete [] r;
  delete [] d;
}

/**
   function solves system of linear equations by the Jacobi iterative method
   
   @param x - the %vector of solution
   @param y - the %vector of right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   
   17. 9. 2015, JK
*/
void densemat::jacobi (double *x, double *y,long ni,double res,long &ani,double &ares,FILE */*out*/)
{
  long i,j,k;
  double norm,rr,*r;
  
  r = new double [n];
  
  //  the vector of solution is set to zero vector
  for (i=0;i<n;i++){
    x[i]=0.0;
  } 
  
  //  iteration loop
  for (i=0;i<ni;i++){
    
    copyv (x,r,n);
    
    //  loop over the number of rows in the matrix
    for (j=0;j<n;j++){
      rr=y[j];
      for (k=0;k<n;k++){
	if (k==j)
	  continue;
	else
	  rr -= a[j*n+k]*r[k]; 
      }
      x[j] = rr/a[j*n+j];
    }
    
    //  computation of the residual vector
    mxv_dm (x,r);
    // adds 2 double arrays, 2nd array is multplied by constant, the result is stored in the third array
    addmultv (y,r,-1.0,r,n);
    //  the norm of the residual vector
    norm = sqrt(ss (r,r,n));
    
    fprintf (stdout,"\n iteration   %ld  norres %e",i,norm);
    if (norm<res){
      break;
    }
  }
  
  ani=i;
  ares=norm;

}


/**
   function solves system of linear equations by the Gauss-Seidel iterative method
   
   @param x - the %vector of solution
   @param y - the %vector of right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   
   17. 9. 2015, JK
*/
void densemat::gauss_seidel (double *x, double *y,long ni,double res,long &ani,double &ares,FILE */*out*/)
{
  long i,j,k;
  double norm,rr,*r;
  
  r = new double [n];
  
  //  the vector of solution is set to zero vector
  for (i=0;i<n;i++){
    x[i]=0.0;
  } 
  
  //  iteration loop
  for (i=0;i<ni;i++){
    
    //  loop over the number of rows in the matrix
    for (j=0;j<n;j++){
      rr=y[j];
      for (k=0;k<n;k++){
	if (k==j)
	  continue;
	else
	  rr -= a[j*n+k]*x[k];
      }
      x[j] = rr/a[j*n+j];
    }
    
    //  computation of the residual vector
    mxv_dm (x,r);
    // adds 2 double arrays, 2nd array is multplied by constant, the result is stored in the third array
    addmultv (y,r,-1.0,r,n);
    //  the norm of the residual vector
    norm = sqrt(ss (r,r,n));
    
    fprintf (stdout,"\n iteration   %ld  norres %e",i,norm);
    if (norm<res){
      break;
    }
  }
  
  ani=i;
  ares=norm;

}









/**
   function prints %matrix into output file
   
   @param out - output stream
   
   JK
*/
void densemat::printmat (FILE *out)
{
  long i,j;
  
  /*
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"\n%3ld",i);
    fprintf (out,"\n");
    for (j=0;j<n;j++){
      //fprintf (out," %5.2f",a[i*n+j]);
      fprintf (out,"%3ld %3ld  %20.15le\n",i,j,a[i*n+j]);
    }
    //fprintf (out,";");
  }
  */

  
  /*
  FILE *p;
  
  p = fopen ("matice.txt","w");
  for (i=0;i<n;i++){
    for (j=0;j<n-1;j++){
      fprintf (p," %lf,",a[i*n+j]);
    }
    j=n-1;
    fprintf (p," %lf",a[i*n+j]);
    fprintf (p," XX\n");
  }
  fclose (p);
  */
  fprintf(out, "\nMatrix in dense format:\n");
  fprintf(out, "%6c", ' ');
  for (i=0;i<n;i++)
    fprintf(out, " %16ld", i+1);
  fprintf(out, "\n");
  for (i=0;i<n;i++){
    fprintf(out, "%6ld", i+1);
    for (j=0;j<n;j++){
      fprintf(out, " % 16e", give_entry(i, j));
    }
    fprintf(out, "\n");
  }
  fprintf (out,"\n\n");  
}


/**
   function prints diagonal entries of the %matrix
   
   @param out - output stream
   
   JK
*/
void densemat::printdiag (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"%5ld  %.10e\n",i,a[i*n+i]);
  }
}

/**
   function returns required %matrix entry
   
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
double densemat::give_entry (long ri,long ci)
{
  double e;
  
  return e = a[ri*n+ci];
}

/**
   function adds required %matrix entry
   
   @param e - %matrix entry
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
void densemat::add_entry (double e,long ri,long ci)
{
  a[ri*n+ci]+=e;
}

/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long densemat::give_negm ()
{
  return negm;
}


/**
   function scales %matrix by its diagonal elements
   
   @param d - array, where scale factors from diagonal are stored
   
   JK, 28.5.2008
*/
void densemat::diag_scale (double *d)
{
  long i,j;
  
  for (i=0;i<n;i++){
    d[i]=1.0/sqrt(a[i*n+i]);
  }
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]*=d[i]*d[j];
    }
  }
}

/**
   function estimates spectral radius

   the estimates are based on the Gershgorin circle theorem
   for details see G.H. Golub, C.F. Van Loan: Matrix computations.
   The Johns Hopkins University Press, 3rd edition, page 320, 1996
   
   JK, 27.8.2008
*/
double densemat::estim_spect_radius ()
{
  long i,j;
  double s,sr;
  
  //  estimate of the spectral radius
  sr=0.0;
  for (i=0;i<n;i++){
    s=0.0;
    for (j=0;j<n;j++){
      s+=fabs(a[i*n+j]);
    }
    if (s>sr)
      sr=s;
  }
  
  return sr;
}

/**
   function assembles %matrix in dense format form scr format
   
   scradr - array with adr in scr format
   scrci - array with ci in scr format
   scra - array with a in scr format
   
   12.02.2009 JB
*/
void densemat::assemble_dense_from_scr(long *scradr,long *scrci,double *scra,long neq)
{
  long i,j;
  long address1,address2,index;
  
  // number of rows of matrix b (a)
  n = neq;
  // alocation of matrix b
  a = new double [n*n];
  
  for(i = 0; i < n*n; i++){
     a[i] = 0.0;
   }
  // rewritten of scr format to dense format
  for(i = 0; i < n; i++){
    address1 = scradr[i];
    address2 = scradr[i+1];
    for(j = address1; j < address2; j++){
      index = scrci[j];
      a[i*n+index]=scra[j];
      a[index*n+i]=scra[j];
    }
  }
  
}

