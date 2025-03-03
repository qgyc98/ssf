#include "cr.h"
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "precond.h"
#ifdef INC_OPENMP
#include <omp.h>
#endif


comprow::comprow (void)
{
  //  number of rows of the matrix
  n=0;
  //  number of entries in the global matrix
  negm=0;
  //  threshold for entries rejection
  limit = 0.0;
  //  decomposition indicator
  decompid=0;
  
  adr=NULL;
  ci=NULL;
  a=NULL;
  adra=NULL;
  aux=NULL;
}



comprow::~comprow ( void )
{
  delete [] adr;
  delete [] ci;
  delete [] a;
  
  delete [] adra;

  if (aux != NULL)
    delete [] aux;
}



/**
   function returns indicator of decomposition (factorization) of %matrix
   
   JK
*/
long comprow::decomp ()
{
  return decompid;
}



/**
   function changes indicator of decomposition (factorization) of %matrix
   
   JK
*/
void comprow::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}



/**
   function allocates array containing addresses of the first entries in rows
   function also allocates auxiliary array
   
   @param m - number of unknowns in solved problem
   
   JK
*/
void comprow::allocadr (long m)
{
  n=m;
  adr  = new long [n+1];
  adra = new long [n+1];
  
  memset (adr,0,(n+1)*sizeof(*adr));
  memset (adra,0,(n+1)*sizeof(*adra));
}



/**
   function returns status of array a
   
   JK
*/
double* comprow::status ()
{
  return a;
}



/**
   function fills the array of zeros
   
   JK
*/
void comprow::nullmat ()
{
  memset(a, 0, negm*sizeof(*a));
}



/**
   function evaluates number of contributions to the %matrix from one element
   
   @param cn - array containing code numbers of the element
   @param ndofe - number of DOFs on the element
   
   JK, 22.6.2001
*/
void comprow::numcontr_elem (long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<ndofe;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      adr[ii]++;
    }
  }
}



/**
   function evaluates number of contributions to the %matrix
   
   @param top - topology
   
   JK, 22.6.2001
*/
void comprow::numcontr (gtopology *top)
{
  long i,ne,ndofe;
  ivector cn;

  ne=top->ne;
  for (i=0;i<ne;i++){
    ndofe=top->give_ndofe (i);
    reallocv(RSTCKIVEC(ndofe, cn));
    top->give_code_numbers (i,cn.a);
    numcontr_elem (cn.a,ndofe);
  }
}



/**
   function fills auxiliary array by one element
   
   @param cn - array containg code numbers of the element
   @param ndofe - number of DOFs on the element
   
   aux - auxiliary array containing column indices

   JK, 22.6.2001
*/
void comprow::fillarray_elem (long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<ndofe;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      aux[adra[ii]]=jj;
      adra[ii]++;
    }
  }
}



/**
   function fills auxiliary array
   
   @param top - topology
   
   aux - auxiliary array containing column indices

   JK, 22.6.2001
*/
void comprow::fillarray (gtopology *top)
{
  long i,ne,ndofe;
  ivector cn;

  ne=top->ne;
  for (i=0;i<ne;i++){
    ndofe=top->give_ndofe (i);
    reallocv(RSTCKIVEC(ndofe, cn));
    top->give_code_numbers (i,cn.a);
    fillarray_elem (cn.a,ndofe);
  }
}

/**
   function computes addresses of the first entries in the rows
   
   JK, 22.6.2001
*/
void comprow::addresses (void)
{
  long i,j,k;
  
  j=adr[0];  adr[0]=0;
  for (i=1;i<n+1;i++){
    k=adr[i];
    adr[i]=j;
    j+=k;
  }
  
  for (i=0;i<n+1;i++){
    adra[i]=adr[i];
  }
  
  aux = new long [adr[n]];
  memset (aux,0,adr[n]*sizeof(*aux));
}



/**
   function sortes array aux
   function also allocates array for %matrix storage
   
   JK, 22.6.2001
*/
void comprow::sort_and_colindex (void)
{
  long i,j,k,ii,jj,lj,uj,min,prev;
  
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adra[i];  prev=-1;
    for (j=lj;j<uj;j++){
      min=LONG_MAX;
      for (k=j;k<uj;k++){
	if (aux[k]<min){
	  min=aux[k];  ii=k;
	}
      }
      if (min==prev){
	uj--;  j--;
	aux[ii]=aux[uj];
      }
      else{
	jj=aux[j];
	aux[j]=min;
	aux[ii]=jj;
	prev=min;
      }
    }
    adra[i]=uj;
  }

  j=0;
  for (i=0;i<n;i++){
    j+=adra[i]-adr[i];
  }

  //  number of non-zero elements in matrix
  negm=j;
  
  
  
  printf ("\n\n negm %ld\n\n",negm);
  //abort ();
  
  
  
  //  array containing resulting column indices
  delete [] ci;
  ci = new long [negm];
  memset (ci,0,negm*sizeof(*ci));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));

  ii=0;  lj=0;  adr[0]=0;
  for (i=0;i<n;i++){
    uj=adra[i];
    for (j=lj;j<uj;j++){
      ci[ii]=aux[j];  ii++;
    }
    lj=adr[i+1];
    adr[i+1]=ii;
  }

  for (i=0;i<=n;i++){
    adra[i]=adr[i];
  }

  delete [] aux;
  aux=NULL;
}


/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as dense %matrix
   
   @param b - local %matrix in row ordering
   @param cn - array containing code numbers of element
   
   JK, 25.6.2001
*/
void comprow::localize (matrix &b,long *cn)
{
  long i,j,k,ii,jj,kk,ll,lk,uk;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    ll=0;  lk=adr[ii];  uk=adr[ii+1];  kk=lk;
    for (j=0;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0) continue;
      if (ll<jj){
	ll=jj;
	for (k=kk;k<uk;k++){
	  if (ci[k]!=jj)  continue;
	  else{
	    a[k]+=b[i][j];
	    kk=k;
	    break;
	  }
	}
      }
      else{
	ll=jj;
	for (k=kk;k>=lk;k--){
	  if (ci[k]!=jj)  continue;
	  else{
	    a[k]+=b[i][j];
	    kk=k;
	    break;
	  } 
	}
      }
    }
  }
}

/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as dense %matrix
   
   @param b - local %matrix in row ordering
   @param cn - array containing code numbers of element
   @param n - order of %matrix b (number of rows or columns)
   
   JK, 25.6.2001
*/
void comprow::localized (double *b,long *cn,long nc)
{
  long i,j,k,ii,jj,kk,ll,lk,uk;
  
  for (i=0;i<nc;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    ll=0;  lk=adr[ii];  uk=adr[ii+1];  kk=lk;
    for (j=0;j<nc;j++){
      jj=cn[j]-1;
      if (jj<0) continue;
      if (ll<jj){
	ll=jj;
	for (k=kk;k<uk;k++){
	  if (ci[k]!=jj)  continue;
	  else{
	    a[k]+=b[i*nc+j];
	    kk=k;
	    break;
	  }
	}
      }
      else{
	ll=jj;
	for (k=kk;k>=lk;k--){
	  if (ci[k]!=jj)  continue;
	  else{
	    a[k]+=b[i*nc+j];
	    kk=k;
	    break;
	  } 
	}
      }
    }
  }
}

/**
   function rejects zero entries from the global %matrix
   function returns number of rejected entries
   
   @param limit - threshold for rejection
   
   JK, 25.6.2001
*/
long comprow::minimize (double limit)
{
  long i,j,n1,cor;
  
  cor=-1;  n1=0;
  for (i=0;i<n;i++){
    adr[i]=adra[i]-n1;
    for (j=adra[i];j<adra[i+1];j++){
      cor++;
      if (fabs(a[j])>limit){
	a[cor]=a[j];
	ci[cor]=ci[j];
      }
      else{
	cor--;  n1++;
      }
    }
  }
  adr[n]=adra[n]-n1;
  negm=adr[n];
  delete [] adra;
  adra=NULL;
  return n1;
}


/**
   function initiates compressed rows storage
   
   @param top - pointer to general topology
   @param ndof - number of rows/columns of the %matrix
   @param mespr - indicator of message printing
   
   JK
*/
void comprow::initiate (gtopology *top,long ndof,long mespr)
{
  if (status()==NULL){
    allocadr (ndof);
    numcontr (top);
    addresses ();
    fillarray (top);
    sort_and_colindex ();
  }
  else{
    nullmat ();
  }
  
  if (mespr==1)  fprintf (stdout,"\n number of matrix entries  %ld",negm);
}

/**
   function prints %matrix into output file
   
   @param out - output stream
   
   JK
*/
void comprow::printmat (FILE *out)
{
  /*
  long i,j;
  double *am;
  
  am = new double [n*n];
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      am[i*n+j]=0.0;
    }
  }
  
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      am[i*n+ci[j]]=a[j];
    }
  }
  
  for (i=0;i<n;i++){
    fprintf (out,"\n");
    for (j=0;j<n;j++){
      fprintf (out," %5.2f",am[i*n+j]);
    }
  }
  */

  //delete [] am;


  /*
    long k=0;
  for (i=0;i<n;i++){
    fprintf (out,"\n");
    for (j=0;j<n;j++){
      if (j==ci[k]){
	fprintf (out," %10.5f",a[k]);
	k++;
      }
      else{
	fprintf (out," %10.5f",0.0);
      }
    }
  }
  */
  
  
  /*
  long i,j;
  FILE *p;
  
  p = fopen ("matice.txt","w");
  fprintf (p,"%ld %ld\n",n,negm);
  //fprintf (p,"S=sparse(%ld,%ld);\n",n,n);
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      //fprintf (p,"S(%ld,%ld)=%e;\n",i+1,ci[j]+1,a[j]);
      //fprintf (p,"%8ld %8ld  % 12.9e\n",i+1,ci[j]+1,a[j]);
      fprintf (p,"%8ld %8ld  % 20.15e\n",i,ci[j],a[j]);
    }
  }
  fclose (p);
  */
  
  long i,j;
  fprintf (out,"\n\n JKKJ \n\n");
  fprintf (out,"%ld  %ld\n",n,negm);
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"%8ld %8ld  % 20.15e\n",i,ci[j],a[j]);
    }
  }
  fprintf (out,"\n\n JKKJ \n\n");

  
  /*
  long i,j;
  FILE *rif,*cif,*mat;
  rif = fopen ("ri.txt","w");
  cif = fopen ("ci.txt","w");
  mat = fopen ("mat.txt","w");
  fprintf (mat,"%ld\n",negm);
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (rif,"%ld\n",i+1);
      fprintf (cif,"%ld\n",ci[j]+1);
      fprintf (mat,"%16.12e\n",a[j]);
    }

  }
  fclose (rif);
  fclose (cif);
  fclose (mat);
  */


}


/**
   function prints diagonal entries of %matrix stored in the compressed row storage scheme
   
   @param out - output stream
   
   JK, 17.3.2007
*/
void comprow::printdiag (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n diagonal entries of matrix stored in the compressed row storage scheme\n");
  
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      if (ci[j]==i){
	fprintf (out,"%f\n",a[ci[j]]);
	break;
      }
    }
  }
}


/**
   function returns required %matrix entry
   
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
double comprow::give_entry (long ri,long rci)
{
  long i;
  double e=0.0;
  
  for (i=adr[ri];i<adr[ri+1];i++){
    if (rci==ci[i]){
      e=a[i];
      break;
    }
  }
  
  return e;
}


/**
   function adds required matrix entry
   
   @param e - matrix entry
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
void comprow::add_entry (double e,long ri,long rci)
{
  long i;
  
  for (i=adr[ri];i<adr[ri+1];i++){
    if (rci==ci[i]){
      a[i]+=e;
      break;
    }
  }
}


/**
   function multiplies %matrix by %vector
   
   @param b - array containing %vector b
   @param c - array containing resulting %vector c = A.b
   
   JK
*/
void comprow::mxv_cr (double *b,double *c)
{
  #ifdef INC_OPENMP
  //  pocet vlaken
  int nthreads=2;
  #pragma omp parallel num_threads(nthreads)
  {
    int id;
    long i,j,ii,lj,uj;
    double s;
    
    //  thread id
    id = omp_get_thread_num();

    for (i=id;i<n;i+=nthreads){
      lj=adr[i];  uj=adr[i+1];
      s=0.0;
      for (j=lj;j<uj;j++){
	ii=ci[j];
	s+=a[j]*b[ii];
      }
      c[i]=s;
    }
  }
  #else
  long i,j,ii,lj,uj;
  double s;
  
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];
    s=0.0;
    for (j=lj;j<uj;j++){
      ii=ci[j];
      s+=a[j]*b[ii];
    }
    c[i]=s;
  }
  #endif

}



/**
   function multiplies %matrix by %vector
   
   @param b - array containing %vector b
   @param c - array containing resulting %vector c = A.b
   interleaving of 2 adjacent rows
   2003 IS
*/

void comprow::mxv_cr_new (double *b,double *c)
{
  long i,j,ii,j2,lj,uj,uj2,p1,p2;
  double s1,s2;

  lj=adr[0];  
  for (i=0;i<(n-1);i+=2)
  {
    uj=adr[i+1];
    uj2=adr[i+2];
    p1=uj-lj;
    p2=uj2-uj;
    s1=0.0;
    s2=0.0;
    if (p2>=p1)
    {
      for (j=lj,j2=uj;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j2++;
      }
      for (;j2<uj2;j2++)
      {
        s2+=a[j2]*b[ci[j2]];
      }      
    }
    else
    {
      for (j=lj,j2=uj;j2<uj2;j2++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j++;
      }
      for (;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
      }             
    }   
    lj=uj2;
    c[i]=s1;
    c[i+1]=s2;    
  }
  if (n&1)
  {
    lj=adr[n-1];  uj=adr[n];
    s1=0;
    for (j=lj;j<uj;j++)
    {
      ii=ci[j];
      s1+=a[j]*b[ii];
    }
    c[n-1]=s1;
  }  
}


/**
   function multiplies %matrix by %vector =%vector2
   return %vector2 * %vector
   @param b - array containing %vector b
   @param c - array containing resulting %vector c = A.b
   interleaving of 2 adjacent rows
   2003 IS
*/

double comprow::mxv_cr_new2 (double *b,double *c)
{
  long i,j,ii,j2,lj,uj,uj2,p1,p2;
  double s1,s2,suma;

  suma=0.0;
  lj=adr[0];  
  for (i=0;i<(n-1);i+=2)
  {
    uj=adr[i+1];
    uj2=adr[i+2];
    p1=uj-lj;
    p2=uj2-uj;
    s1=0.0;
    s2=0.0;
    if (p2>=p1)
    {
      for (j=lj,j2=uj;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j2++;
      }
      for (;j2<uj2;j2++)
      {
        s2+=a[j2]*b[ci[j2]];
      }      
    }
    else
    {
      for (j=lj,j2=uj;j2<uj2;j2++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j++;
      }
      for (;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
      }             
    }   
    lj=uj2;
    c[i]=s1;
    c[i+1]=s2;    
    suma+=s1*b[i];
    suma+=s2*b[i+1];
  }
  if (n&1)
  {
    lj=adr[n-1];  uj=adr[n];
    s1=0;
    for (j=lj;j<uj;j++)
    {
      ii=ci[j];
      s1+=a[j]*b[ii];
    }
    c[n-1]=s1;
    suma+=s1*b[n-1];
  }  
  return suma;
}


/**
   function multiplies transposed %matrix by %vector b, result is stored in %vector c
   
   @param b - input %vector
   @param c - output %vector
   
   JK, 16.8.2000/4.5.2002
*/
void comprow::mtxv_cr (double *b,double *c)
{
  long i,j,ii,lj,uj;
  double s;
  
  nullv (c,n);

  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];  s=b[i];
    for (j=lj;j<uj;j++){
      ii=ci[j];
      c[ii]+=a[j]*s;
    }
  }
}

void comprow::mv_cr15 (double *b,double *c)
/*
  funkce nasobi matici a s vektorem b, vysledek je vektor c
  matice a je ulozena v compressed row
  
  prokladani dvou sousednich
  
  vystup
  c - vysledny vektor
  
  vstupy
  a - matice v compressed row
  b - vektor
  adr - pole adres prvnich prvku v radku
  ci - pole sloupcovych indexu
  n - rozmer matice a
  
  28.7.1997
*/
{
  long i,j,ii,j2,lj,uj,uj2,p1,p2;
  double s1,s2;

  lj=adr[0];  
  for (i=0;i<(n-1);i+=2)
  {
    uj=adr[i+1];
    uj2=adr[i+2];
    p1=uj-lj;
    p2=uj2-uj;
    s1=0.0;
    s2=0.0;
    if (p2>=p1)
    {
      for (j=lj,j2=uj;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j2++;
      }
      for (;j2<uj2;j2++)
      {
        s2+=a[j2]*b[ci[j2]];
      }      
    }
    else
    {
      for (j=lj,j2=uj;j2<uj2;j2++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j++;
      }
      for (;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
      }          	
    }	
    lj=uj2;
    c[i]=s1;
    c[i+1]=s2;    
  }
  if (n&1)
  {
    lj=adr[n-1];  uj=adr[n];
    s1=0;
    for (j=lj;j<uj;j++)
    {
      ii=ci[j];
      s1+=a[j]*b[ii];
    }
    c[n-1]=s1;
  }

  
}

/**
   function multiplies %matrix by %vector b1, result is stored in %vector c1
                  multiplies transposed %matrix by %vector b2, result is stored in %vector c2
                  return c1*b2   
   @param b1,b2 - input %vectors
   @param c1,c2 - output %vectors
   
   2003 IS
*/

double comprow::mxv_cr_pom (double *b1,double *c1,double *b2,double *c2)
{
  long i,j,ii,lj,uj;
  double s1,s2,pom,suma;
  
  suma=0.0;
  nullv (c2,n);
  
  for (i=0;i<n;i++)
  {
    lj=adr[i];
    uj=adr[i+1];
    s2=b2[i];
    s1=0.0;
    for (j=lj;j<uj;j++)
    {
      ii=ci[j];
      pom=a[j];
      s1+=pom*b1[ii];
      c2[ii]+=pom*s2;
    }
    c1[i]=s1;
    suma+=s1*s2;
  }
  return suma;
}

/**
   function adds multiplied %matrix stored in cr by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param cr - another compressed row storage
   
   JK
*/
void comprow::addmat_cr (double c,comprow &cr)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]+=c*cr.a[i];
  }
}

/**
   function multiplies %matrix by coefficient c
   
   @param c - multiplicative coefficient
   
   JK
*/
void comprow::scalmat_cr (double c)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]*=c;
  }
}

/**
   function solves system of linear algebraic equations
   by conjugate gradient method, %matrix is stored as compressed rows
   
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
   
   JK, 17.7.2001
*/
void comprow::cg (double *x,double *y,
		  long ni,double res,long &ani,double &ares,double zero,long iv)
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
  mxv_cr (x,p);

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
    mxv_cr (d,p);
    //mv_cr15 (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      fprintf (stdout,"\n there is zero denominator (%e) in alpha computation in conjugate gradient method (cr.cpp)\n",fabs(denom));
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
    
    if (i%100==0)
    fprintf (stdout,"\n iteration   %ld  norres  %e    norres/norrhs %e",i,nom,nom/nory);
    //fprintf (stdout,"\n iteration   %ld  norres/norrhs %e",i,sqrt(nom/nory));
    
    if (sqrt(nom/nory)<res)  break;
    //if (fabs(nom)<limit)  break;
    
    
    beta = nom/denom;
    
    //  new vector of direction
    for (j=0;j<n;j++){
      d[j]=beta*d[j]-r[j];
    }
  }
  
  fprintf (stdout,"\n iteration   %ld  norres  %e    norres/norrhs %e",i,nom,nom/nory);
  ani=i;  ares=nom;
  
  delete [] p;  delete [] r;  delete [] d;
}


/**
   function solves system of linear algebraic equations
   by conjugate gradient method, %matrix is stored as compressed rows
   
   @param x - %vector of unknowns
   @param y - %vector of right hand side
   @param ni - maximum number of iterations
   @param res - required error of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   @param zero - computer zero
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero vector
   iv=1 - initial vector is taken from x array
   
   JK, 17.7.2001
*/
void comprow::cg_prec (precond &pr,double *x,double *y,
		       long ni,double res,long &ani,double &ares,double zero,long iv)
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
  mxv_cr (x,p);

  //  norm of the right hand side
  nory=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=y[i]-p[i];
  }
  
  if (nory<zero){
    print_err("norm of right hand side in conjugate gradient method is smaller than %e", __FILE__, __LINE__, __func__,zero);
    ares=nory;  ani=0;
    return;
  }
  
  //  h = M^-1 r
  pr.preconditioning (r,h,r);

  nom = ss (r,h,n);
  
  //  direction vector
  for (i=0;i<n;i++){
    d[i]=h[i];
  }


  //  iteration loop
  for (i=0;i<ni;i++){
    
    
    //  new coefficient alpha
    mxv_cr (d,p);
    //mv_cr15 (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      fprintf (stdout,"\n there is zero denominator (%e) in alpha computation in conjugate gradient method (cr.cpp)\n",fabs(denom));
      break;
    }
    
    //  new alpha
    alpha = nom/denom;
    
    //  new approximation of x and r
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
    }
    
    //  h = M^-1 r
    //if (i<5)
      pr.preconditioning (r,h,r);
      /*
    else{
      for (j=0;j<n;j++){
	h[j]=r[j];
      }
    }
      */

    denom=nom;
    nom = ss (r,h,n);
    beta = nom/denom;

    //double norres = ss (r,r,n);
    
    //fprintf (stdout,"\n iteration   %ld   alpha %e  beta %e   norres/norrhs %e",i,alpha,beta,sqrt(nom/nory));
    //if (i%100==0)
    //fprintf (stdout,"\n iteration   %ld   norres  %e  norrhs %e   norres/norrhs %e",i,sqrt(nom),sqrt(nory),sqrt(nom/nory));
    fprintf (stdout,"\n iteration   %ld   norres  %e  norrhs %e   norres/norrhs %e",i,sqrt(nom),sqrt(nory),sqrt(nom/nory));


    //if (norres<res){

    if (sqrt(nom/nory)<res){
      ani=i;  ares=nom/nory;
      break;
    }
    
    
    //  new vector of direction
    for (j=0;j<n;j++){
      d[j]=beta*d[j]+h[j];
    }
  }
  
  ani=i;  ares=nom;
  
  delete [] h;  delete [] p;  delete [] r;  delete [] d;
}


/**
   function solves system of linear algebraic equations
   by conjugate gradient method, %matrix is stored as compressed rows
   
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
   2 fussions + new multiplication   
   2003 IS
*/
void comprow::cg_new (double *x,double *y,
          long ni,double res,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double nom,denom,nory,alpha,beta,nomx;
  double *d,*r,*p;
  
  d = new double [n];
  r = new double [n];
  p = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  mxv_cr_new (x,p);

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
    return;
  }
  

  //  iteration loop
  for (i=0;i<ni;i++){
    
    
    //  new coefficient alpha
    //mxv_cr (d,p);    
    //denom = ss (d,p,n);
    denom=mxv_cr_new2(d,p);
    if (fabs(denom)<zero){
      fprintf (stdout,"\n there is zero denominator in alpha computation in conjugate gradient method (cr.cpp)\n");
      break;
    }
    
    alpha = nom/denom;
    
    //  new approximation of x and r
    nomx=0.0;
    for (j=0;j<n;j++){
      r[j]+=alpha*p[j]; 
      x[j]+=alpha*d[j];
      nomx+=r[j]*r[j];
    }
    
    denom=nom;
    
    //nom = ss (r,r,n);
    nom = nomx;

    if (i%100==0)
    fprintf (stdout,"\n iteration   %ld  norres  %e    norres/norrhs %e",i,nom,nom/nory);
    //fprintf (stdout,"\n iteration   %ld  norres/norrhs %e",i,nom/nory);
    
    if (nom/nory<res)  break;
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

/**
   function solves system of linear algebraic equations
   by bi-conjugate gradient method
   
   @param x - array containing solution
   @param y - array containing right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   @param zero - computer zero
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero vector
   iv=1 - initial vector is taken from x array
   
   JK, 16.8.2000/4.5.2002
*/
void comprow::bicg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double normy,normr,normrr,nom,denom,alpha,beta;
  double *r,*rr,*d,*dd,*p,*pp;

  r = new double [n];  rr = new double [n];
  d = new double [n];  dd = new double [n];
  p = new double [n];  pp = new double [n];
  
  normy = ss (y,y,n);
  
  if (normy<zero){
    print_err("norm of right hand side in conjugate gradient method is smaller than %e", __FILE__, __LINE__, __func__,zero);
    ares=normy;  ani=0;
    return;
  }

  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  
  mxv_cr (x,p);
  for (i=0;i<n;i++){
    r[i]=y[i]-p[i];
    d[i]=r[i];
    rr[i]=r[i];
    dd[i]=rr[i];
  }
  
  nom = ss (r,rr,n);
  
  for (i=0;i<ni;i++){
    mxv_cr (d,p);
    mtxv_cr (dd,pp);
    
    denom = ss (dd,p,n);
    alpha = nom/denom;
    
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
      rr[j]-=alpha*pp[j];
    }

    normr = ss (r,r,n);
    normrr = ss (rr,rr,n);
    if (i%100==0)
    fprintf (stdout,"\n iteration %ld     %e     %e",i,normr/normy,normrr/normy);
    if (normr/normy<res && normrr/normy<res) break;
    
    
    denom = nom;
    nom = ss (r,rr,n);
    beta = nom/denom;

    for (j=0;j<n;j++){
      d[j]=r[j]+beta*d[j];
      dd[j]=rr[j]+beta*dd[j];
    }
  }
  
  ani=i;  ares=normr;
  delete [] pp;  delete [] p;
  delete [] dd;  delete [] d;
  delete [] rr;  delete [] r;
}

//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************
//**********************************************************

/**
  funkce nasobi matici a s vektorem b, vysledek je vektor c
  matice a je ulozena v compressed row
  
  prokladani dvou sousednich
  
  vystup
  c - vysledny vektor
  vraci c*b
  
  vstupy
  a - matice v compressed row
  b - vektor
  adr - pole adres prvnich prvku v radku
  ci - pole sloupcovych indexu
  n - rozmer matice a
  
  15.4.2003
*/
double comprow::mv_cr15_rev (double *b,double *c)
{
  long i,j,j2,lj,uj,uj2,p1,p2;
  double s1,s2,soucet=0;

  lj=adr[0];  
  for (i=0;i<n;i+=2)
  {
    uj=adr[i+1];
    uj2=adr[i+2];
    p1=uj-lj;
    p2=uj2-uj;
    s1=0.0;
    s2=0.0;
    if (p2>=p1)
    {
      for (j=lj,j2=uj;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j2++;
      }
      for (;j2<uj2;j2++)
      {
        s2+=a[j2]*b[ci[j2]];
      }      
    }
    else
    {
      for (j=lj,j2=uj;j2<uj2;j2++)
      {
        s1+=a[j]*b[ci[j]];
        s2+=a[j2]*b[ci[j2]];
        j++;
      }
      for (;j<uj;j++)
      {
        s1+=a[j]*b[ci[j]];
      }                 
    }   
    lj=uj2;
    c[i]=s1;    
    c[i+1]=s2;    
    soucet+=s1*b[i]+s2*b[i+1];
  }
  if(n&1)
  {
    s1=0.0;
    for (j=adr[n-1];j<adr[n];j++)
    {
      s1+=a[j]*b[ci[j]];
    }
    c[n-1]=s1;
  
    soucet+=s1*b[n-1];
  }
  return soucet;
}
 
/**
  funkce resi soustavu linearnich algebraickych rovnic
  metodou sdruzenych gradientu (podle Axelssona)
  
  matice soustavy je ulozena v compresovane forme
  
  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  ci - pole sloupcovych indexu
  adr - pole adres prvnich prvku
  n - pocet neznamych
  ni - maximalni pocet iteraci
  res - maximalni pripustna chyba
  ani - pocet skutecne provedenych operaci
  ares - norma vektoru rezidui
  limit - konstanta pro testovani na nulu
  
  IS,  10.5.2003
*/
void comprow::cg_cr_rev (double *x,double *y,
			 long ni,double res,long &ani,double &ares,double limit,long iv)
{
  long i,j;
  double nom,denom,nory,alpha,beta,ener;
  double *d,*r,*p;
  
  d = new double [n];
  r = new double [n];
  p = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  mv_cr15_rev (x,p);
  
  nory=0.0;  nom=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=0.0-y[i];
    d[i]=y[i];
    nom+=r[i]*r[i];
  }
  
  if (nory<limit){
    printf ("\n\n norma vektoru prave strany v metode sdruzenych gradientu");
    printf ("\n je mensi nez %e",limit);
    ares=nory;  ani=0;
    return;
  }
  
  for (i=0;i<ni;i++)
  {
    /*  vypocet nove alphy  */
    denom=mv_cr15_rev (d,p);
    
    //denom = ss (d,p,n);
    if (fabs(denom)<limit){
      printf ("\n V metode sdruzenych gradientu je nulovy jmenovatel u alfa");
      abort ();
    }
    
    alpha = nom/denom;
    ener=0;
    /*  nova aproximace x a r */
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]+=alpha*p[j];
      ener+=r[j]*r[j];
    }
    
    //ener = energie_cr (a,x,y,adr,ci,n);
    
    //printf ("\n iterace   %ld   %e    %e",i,nom/nory,ener);
    //printf ("\n iterace   %ld   %e",i,nom/nory);
    
    denom=nom;

    /*  vypocet beta  */
    nom=ener;
    //nom = ss (r,r,n);
    
    
    if (nom/nory<res){
      printf ("\n iterace   %ld   %e",i,nom/nory);
      
      break;
    }
    if (fabs(nom)<limit){
      printf ("\n iterace   %ld   %e",i,nom/nory);
      
      break;
    }
    
    beta = nom/denom;
    
    /*  vypocet noveho d */
    for (j=n-1;j>=0;j--){
      d[j]=beta*d[j]-r[j];
    }
  }
  
  ani=i;  ares=nom;
  
  delete [] p;  delete [] r;  delete [] d;
}

/**
   function solves system of linear algebraic equations
   by bi-conjugate gradient method
   
   @param x - array containing solution
   @param y - array containing right hand side
   @param ni - maximum number of iterations
   @param res - required norm of residual %vector
   @param ani - number of performed iterations
   @param ares - attained norm of residual %vector
   @param zero - computer zero
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero vector
   iv=1 - initial vector is taken from x array
   
   modified multiplication,3 fussions   
   2003 IS
*/
void comprow::bicg_new (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double normy,normr,normrr,nom,denom,alpha,beta,nomx;
  double *r,*rr,*d,*dd,*p,*pp;

  r = new double [n];  rr = new double [n];
  d = new double [n];  dd = new double [n];
  p = new double [n];  pp = new double [n];
  
  normy = ss (y,y,n);
  
  if (normy<zero){
    print_err("norm of right hand side in conjugate gradient method is smaller than %e", __FILE__, __LINE__, __func__,zero);
    ares=normy;  ani=0;
    return;
  }

  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  
  mxv_cr (x,p);
  for (i=0;i<n;i++){
    r[i]=y[i]-p[i];
    d[i]=r[i];
    rr[i]=r[i];
    dd[i]=rr[i];
  }
  
  nom = ss (r,rr,n);
  
  for (i=0;i<ni;i++){
    //mxv_cr (d,p);
    //mtxv_cr (dd,pp);    
    //denom = ss (dd,p,n);
    denom=mxv_cr_pom (d,p,dd,pp);
    
    alpha = nom/denom;
    normr = normrr = nomx = 0.0;
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
      normr += r[j]*r[j];
      rr[j]-=alpha*pp[j];
      normrr+= rr[j]*rr[j];
      nomx+=r[j]*rr[j];
    }

    //normr = ss (r,r,n);
    //normrr = ss (rr,rr,n);
    fprintf (stdout,"\n iteration %ld     %e     %e",i,normr/normy,normrr/normy);
    if (normr/normy<res && normrr/normy<res) break;
    
    
    denom = nom;
    //nom = ss (r,rr,n);
    nom=nomx;
    beta = nom/denom;

    for (j=0;j<n;j++){
      d[j]=r[j]+beta*d[j];
      dd[j]=rr[j]+beta*dd[j];
    }
  }
  
  ani=i;  ares=normr;
  delete [] pp;  delete [] p;
  delete [] dd;  delete [] d;
  delete [] rr;  delete [] r;
}

/**
   function copies data to additional storage
   
   @param cr - additional storage (it is filled)
   
   JK, 29.9.2006
*/
void comprow::copy_cr (comprow &cr)
{
  long i;
  
  cr.n=n;
  cr.negm=negm;
  
  if (cr.adr!=NULL){
    delete [] cr.adr;
  }
  cr.adr=new long [n+1];

  if (cr.ci!=NULL){
    delete [] cr.ci;
  }
  cr.ci=new long [negm];
  
  if (cr.a!=NULL){
    delete [] cr.a;
  }
  cr.a=new double [negm];
  
  for (i=0;i<=n;i++){
    cr.adr[i]=adr[i];
  }
  for (i=0;i<negm;i++){
    cr.a[i]=a[i];
    cr.ci[i]=ci[i];
  }
}


/**
   function fills new object by data
   
   @param nsub - number of rows of %matrix
   @param smadr - array of addresses of first nonzero entries in rows
   @param smci - array of column indices
   @param sma - array of nenzero %matrix entries
   
   JK, 19.2.2007
*/
void comprow::fill_data (long nsub,long *smadr,long *smci,double *sma)
{
  long i;

  //  number of rows
  n=nsub;
  
  //  addresses of first nonzero matrix entries
  if (adr!=NULL){
    delete [] adr;
  }
  adr = new long [n+1];
  
  for (i=0;i<n+1;i++){
    adr[i]=smadr[i];
  }
  
  //  number of nonzero matrix entries
  negm = adr[n];
  
  //  column indices
  if (ci!=NULL){
    delete [] ci;
  }
  ci = new long [negm];
  
  for (i=0;i<negm;i++){
    ci[i]=smci[i];
  }
  
  //  nonzero matrix entries
  if (a!=NULL){
    delete [] a;
  }
  a = new double [negm];
  
  for (i=0;i<negm;i++){
    a[i]=sma[i];
  }
}



/**
   function selects submatrix from the original %matrix stored in the compressed
   row storage scheme, the submatrix is defined by indices
   
   @param li - list of indices of submatrix
   @param nsub - number of rows and columns of submatrix (it is also length of the array li)
   @param smcr - pointer to object of the class comprow which stores submatrix
   
   JK, 19.2.2007
*/
void comprow::select_submatrix (long *li,long nsub,comprow *smcr)
{
  long i,j,k,m,ii,minm,maxm,mins,maxs;
  long *smadr,*smci,*newnumb;
  double *sma;
  
  smadr = new long [nsub+1];
  
  //  minimum index in the submatrix
  mins=li[0]-1;
  //  maximum index in the submatrix
  maxs=li[nsub-1]-1;

  //
  //  computation of numbers of nonzero matrix entries in particular rows defined by submatrix
  //
  for (i=0;i<nsub;i++){
    smadr[i]=0;
    //  number of actual row
    m=li[i]-1;
    //  minimum index in the m-th row of the matrix
    minm=ci[adr[m]];
    //  maximum index in the m-th row of the matrix
    maxm=ci[adr[m+1]-1];
    
    if (maxm<mins)
      continue;
    if (maxs<minm)
      continue;
    
    j=adr[m];  k=0;
    do{
      if (ci[j]==li[k]-1){
	smadr[i]++;
	j++;  k++;
      }
      else{
	if (ci[j]<li[k]-1){
	  j++;
	}
	else{
	  k++;
	}
      }
    }while (j<adr[m+1] && k<nsub);
  }
  
  //  array smadr contains numbers of nonzero matrix entries in particular rows at this moment
  
  j=0;
  for (i=0;i<nsub;i++){
    j+=smadr[i];
  }
  
  smci = new long [j];
  sma = new double [j];
  
  ii=0;
  for (i=0;i<nsub;i++){
    //  number of actual row
    m=li[i]-1;
    //  minimum index in the m-th row of the matrix
    minm=ci[adr[m]];
    //  maximum index in the m-th row of the matrix
    maxm=ci[adr[m+1]-1];
    
    if (maxm<mins)
      continue;
    if (maxs<minm)
      continue;
    
    j=adr[m];  k=0;
    do{
      if (ci[j]==li[k]-1){
	sma[ii]=a[j];
	smci[ii]=ci[j];
	ii++;
	j++;  k++;
      }
      else{
	if (ci[j]<li[k]-1){
	  j++;
	}
	else{
	  k++;
	}
      }
    }while (j<adr[m+1] && k<nsub);
    
  }

  m=smadr[0];  smadr[0]=0;
  for (j=1;j<nsub+1;j++){
    k=smadr[j];
    smadr[j]=m;
    m+=k;
  }
  
  //  array smadr contains addresses of first nonzero matrix entries in submatrix at this moment
  
  
  //  array of column indices must be reorganized
  //  original unknown numbers 3 4 7 9 must be converted to unknow numbers 0 1 2 3
  //  array ci must be converted with respect to this mapping
  
  newnumb = new long [maxs+1];
  for (i=0;i<maxs+1;i++){
    newnumb[i]=-1;
  }
  for (i=0;i<nsub;i++){
    j=li[i]-1;
    newnumb[j]=1;
  }
  j=0;
  for (i=0;i<maxs+1;i++){
    if (newnumb[i]==1){
      newnumb[i]=j;
      j++;
    }
  }
  
  for (i=0;i<smadr[nsub];i++){
    j=smci[i];
    smci[i]=newnumb[j];
  }
  
  
  //  submatrix is stored in new object of the class comprow
  smcr->fill_data (nsub,smadr,smci,sma);
  
  delete [] sma;
  delete [] smci;
  delete [] smadr;
  delete [] aux;
  aux = NULL;
}


/**
   function selects submatrix from the original %matrix stored in the compressed
   row storage scheme, the submatrix is defined by indices
   
   @param li - list of indices of submatrix
   @param nsub - number of rows and columns of submatrix (it is also length of the array li)
   @param smdm - pointer to object of the class densemat which stores submatrix
   
   JK, 14.3.2007
*/
void comprow::select_submatrix (long *li,long nsub,densemat *smdm)
{
  long i,j,k,ii,jj,kk;
  
  smdm->n = nsub;
  smdm->a = new double [nsub*nsub];
  for (i=0;i<nsub;i++){
    for (j=0;j<nsub;j++){
      smdm->a[i*nsub+j]=0.0;
    }
  }
  
  for (i=0;i<nsub;i++){
    //  row index in the original matrix
    ii=li[i]-1;
    for (j=adr[ii];j<adr[ii+1];j++){
      //  column index in the new matrix
      jj=ci[j];  kk=-1;
      for (k=0;k<nsub;k++){
	if (li[k]-1==jj){
	  kk=k;
	  break;
	}
      }
      if (kk>-1)
	smdm->a[i*nsub+k]=a[j];
    }
  }

}


/**
   function computes dot product of the ri-th row with sparse %vector sv
   
   @param cvi - input compressed sparse %vector
   @param cvo - output compressed sparse %vector
   
   JK, 25.2.2007
*/
void comprow::crxcv_cv (compvect *cvi,compvect *cvo)
{
  long i,j,k;
  long minm,maxm,minv,maxv,ncm,ncv;
  long *ind;
  double s,*b,*c;
  
  if (cvi->ordering == 0){
    print_err("compressed sparse vector is not ordered increasingly", __FILE__, __LINE__, __func__);
  }
  
  //  auxiliary array
  c = new double [n];

  //  number of nonzero components in the vector
  ncv=cvi->nz;
  //  minimum index in the vector
  minv=cvi->ind[0];
  //  maximum index in the vector
  maxv=cvi->ind[ncv-1];
  
  ind=cvi->ind;
  b=cvi->a;
  
  for (i=0;i<n;i++){
    c[i]=0.0;
    
    //  number of nonzero components in the i-th row of the matrix
    ncm=adr[i+1]-adr[i];
    //  minimum index in the i-th row of the matrix
    minm=ci[adr[i]];
    //  maximum index in the i-th row of the matrix
    maxm=ci[adr[i+1]-1];
    
    if (maxm<minv)
      continue;
    if (maxv<minm)
      continue;
    
    j=adr[i];  k=0;
    s=0.0;
    do{
      if (ci[j]==ind[k]){
	s+=a[j]*b[k];
	j++;  k++;
      }
      else{
	if (ci[j]<ind[k]){
	  j++;
	}
	else{
	  k++;
	}
      }
    }while (j<adr[i+1] && k<ncv);
    c[i]=s;
  }
  
  //  conversion of the dense vector b to the compressed sparse vector
  cvo->setup_vector (c,n);
  
  delete [] c;
}

/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long comprow::give_negm ()
{
  return negm;
}

/**
   function estimates spectral radius

   the estimates are based on the Gershgorin circle theorem
   for details see G.H. Golub, C.F. Van Loan: Matrix computations.
   The Johns Hopkins University Press, 3rd edition, page 320, 1996
   
   JK, 27.8.2008
*/
double comprow::estim_spect_radius ()
{
  long i,j;
  double s,sr;
  
  //  estimate of the spectral radius
  sr=0.0;
  for (i=0;i<n;i++){
    s=0.0;
    for (j=adr[i];j<adr[i+1];j++){
      s+=fabs(a[j]);
    }
    if (s>sr)
      sr=s;
  }
  
  return sr;
}
