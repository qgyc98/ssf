#include "cc.h"
#include "diagmat.h"
#include "umfpack.h"
#include "globalg.h"
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef INC_CHMD
 #include <cblas.h>
#endif


compcol::compcol (void)
{
  //  number of columns of the matrix
  n=0L;
  //  number of entries in the global matrix
  negm=0L;
  //  threshold for entries rejection
  limit = 0.0;
  //  decomposition indicator
  decompid=0L;
  
  adr=NULL;
  adra=NULL;
  aux=NULL;
  ri=NULL;
  a=NULL;

#ifdef INC_CHMD
  symbolic = numeric = NULL;
#endif  
}



compcol::~compcol ( void )
{
  delete [] adr;
  delete [] adra;
  if (aux != NULL)
    delete [] aux;
  delete [] ri;
  delete [] a;

#ifdef INC_CHMD  
  // deallocation of memory of UMFPACK structures
  umfpack_dl_free_symbolic(&symbolic);
  umfpack_dl_free_numeric(&numeric);
#endif  
}



/**
   function returns indicator of decomposition (factorization) of %matrix
   
   JK
*/
long compcol::decomp ()
{
  return decompid;
}



/**
   function changes indicator of decomposition (factorization) of %matrix
   
   JK
*/
void compcol::changedecomp ()
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
void compcol::allocadr (long m)
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
double* compcol::status ()
{
  return a;
}



/**
   function fills the array of zeros
   
   JK
*/
void compcol::nullmat ()
{
  memset(a, 0, negm*sizeof(*a));
}



/**
   function evaluates number of contributions to the %matrix from one element
   
   @param cn - array containing code numbers of the element
   @param ndofe - number of DOFs on the element
   
   JK, 22.6.2001
*/
void compcol::numcontr_elem (long *cn,long ndofe)
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
void compcol::numcontr (gtopology *top)
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
   
   aux - auxiliary array containing row indices

   JK, 22.6.2001
*/
void compcol::fillarray_elem (long *cn,long ndofe)
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
void compcol::fillarray (gtopology *top)
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
   function computes addresses of the first entries in the columns
   
   JK, 22.6.2001
*/
void compcol::addresses (void)
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
void compcol::sort_and_colindex (void)
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
  
  
  
  //  array containing resulting row indices
  delete [] ri;
  ri = new long [negm];
  memset (ri,0,negm*sizeof(*ri));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));

  ii=0;  lj=0;  adr[0]=0;
  for (i=0;i<n;i++){
    uj=adra[i];
    for (j=lj;j<uj;j++){
      ri[ii]=aux[j];  ii++;
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
void compcol::localize (matrix &b,long *cn)
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
	  if (ri[k]!=jj)  continue;
	  else{
	    a[k]+=b[j][i];
	    kk=k;
	    break;
	  }
	}
      }
      else{
	ll=jj;
	for (k=kk;k>=lk;k--){
	  if (ri[k]!=jj)  continue;
	  else{
	    a[k]+=b[j][i];
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
void compcol::localized (double *b,long *cn,long nc)
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
	  if (ri[k]!=jj)  continue;
	  else{
	    a[k]+=b[j*nc+i];
	    kk=k;
	    break;
	  }
	}
      }
      else{
	ll=jj;
	for (k=kk;k>=lk;k--){
	  if (ri[k]!=jj)  continue;
	  else{
	    a[k]+=b[j*nc+i];
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
long compcol::minimize (double limit)
{
  long i,j,n1,cor;
  
  cor=-1;  n1=0;
  for (i=0;i<n;i++){
    adr[i]=adra[i]-n1;
    for (j=adra[i];j<adra[i+1];j++){
      cor++;
      if (fabs(a[j])>limit){
	a[cor]=a[j];
	ri[cor]=ri[j];
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
void compcol::initiate (gtopology *top,long ndof,long mespr)
{
  if (status()==NULL){
    decompid = 0;
    allocadr (ndof);
    numcontr (top);
    addresses ();
    fillarray (top);
    sort_and_colindex ();    
#ifdef INC_CHMD
    openblas_set_num_threads(Numth);
    double *null = (double*)NULL;
    if (symbolic){
      // deallocation of memory of UMFPACK structure
      umfpack_dl_free_symbolic(&symbolic);
    }
    if (numeric){
      // deallocation of memory of UMFPACK structure
      umfpack_dl_free_numeric(&numeric);
    }
    umfpack_dl_symbolic (n, n, adr, ri, a, &symbolic, null, null);
#endif    
  }
  else{    
    nullmat ();
    decompid = 0;
#ifdef INC_CHMD
    if (numeric){
      // deallocation of memory of UMFPACK structure
      umfpack_dl_free_numeric(&numeric);
    }
#endif    
  }
  
  
  if (mespr==1){
#ifdef INC_CHMD
    fprintf(stdout, "\n setting OpenBLAS to use %d threads\n", Numth);
#endif
    fprintf (stdout," number of matrix entries  %ld\n",negm);
  }
}

/**
   function prints %matrix into output file
   
   @param out - output stream
   
   JK
*/
void compcol::printmat (FILE *out)
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
  fprintf (out,"\n\n compressed columns \n\n");
  fprintf (out,"%ld  %ld\n",n,negm);
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"%8ld %8ld  % 20.15e\n",i,ri[j],a[j]);
    }
    fprintf (out,"\n");
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
void compcol::printdiag (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n diagonal entries of matrix stored in the compressed row storage scheme\n");
  
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      if (ri[j]==i){
	fprintf(out, "%f\n", a[ri[j]]);
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
double compcol::give_entry (long rri,long ci)
{
  long i;
  double e=0.0;
  
  for (i=adr[ci];i<adr[ci+1];i++){
    if (rri==ri[i]){
      e=a[i];
      break;
    }
  }
  
  return e;
}



/**
  The function adds required matrix entry.
   
  @param[in] e - matrix entry
  @param[in] rri - row index of the entry
  @param[in] ci - column index of the entry
   
  JK, 24.7.2005
*/
void compcol::add_entry (double e,long rri,long ci)
{
  long i;
  
  for (i=adr[ci];i<adr[ci+1];i++){
    if (rri==ri[i]){
      a[i]+=e;
      break;
    }
  }
}



/**
  The function multiplies the given %matrix by %vector b.
   
  @param[in] b - array containing %vector b
  @param[out] c - array containing resulting %vector c = A.b
   
  Created by TKo, 04.2024
*/
void compcol::mxv_cc (double *b, double *c)
{
  long i,j,ii,lj,uj;

  nullv(c, n);
  for (i=0; i<n; i++){
    lj=adr[i];  uj=adr[i+1];
    for (j=lj; j<uj; j++){
      ii = ri[j];
      c[ii] += a[j]*b[i];
    }
  }
}



/**
  The function adds multiplied %matrix stored in cr by coefficient c to actual %matrix.
  
  @param[in] c - multiplicative coefficient
  @param[in] cc - summand in compressed column storage
   
  JK
*/
void compcol::addmat_cc (double c, compcol &cc)
{
  long i;
  for (i=0;i<negm;i++){
    a[i] += c*cc.a[i];
  }
}



/**
  The function adds diagonal %matrix multiplied by the coefficient c to actual %matrix.
  
  @param[in] c - multiplicative coefficient
  @param[in] d - diagonal %matrix summand  storage
   
  TKo
*/
void compcol::addmat_diag (double c, diagmat &d)
{
  long i, j, ncomp;
  bool diag_found;
  
  for (i=0; i<n; i++){
    ncomp = adr[i+1] - adr[i];
    diag_found = false;
    for (j=0; j<ncomp; j++){
      if (ri[j] == i){
        a[adr[i]+j] += c*d.a[i];
        diag_found = true;
      }      
    }
    if ((diag_found == false) && (d.a[i] != 0.0)){
      print_err("diagonal element on %ld row was not stored in compressed column format matrix",
                __FILE__, __LINE__, __func__, i+1);
      abort();
    }
  }
}



/**
  The function multiplies %matrix by coefficient c.
   
  @param[in] c - multiplicative coefficient
   
  JK
*/
void compcol::scalmat_cc (double c)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]*=c;
  }
}



#ifdef INC_CHMD  

/**
  The function solves system of linear algebraic equations
  by external slover CHOLMOD from SuiteSprase library, 
  %matrix is stored as compressed column format.
   
  @param x - %vector of unknowns
  @param y - %vector of right hand side
   
  TKo, 29.3.2023
*/
void compcol::solve (double *x, double *y)
{
  double *null = (double*)NULL;
  if (decompid == 0){
    umfpack_dl_numeric (adr, ri, a, symbolic, &numeric, null, null);
    decompid = 1;
  }
  umfpack_dl_solve(UMFPACK_A, adr, ri, a, x, y, numeric, null, null); 
  decompid = 1;
}
#else
void compcol::solve (double *x, double *y)
{
  print_err("\n Support for UMFPACK solver was not defined,"
            "\n Use different storage/solver type, e.g. cr with biconjugated gradients\n", __FILE__, __LINE__, __func__);  
  abort();
}
#endif



/**
   function fills new object by data
   
   @param nsub - number of columns of %matrix
   @param smadr - array of addresses of first nonzero entries in columns
   @param smri - array of row indices
   @param sma - array of nenzero %matrix entries
   
   JK, 19.2.2007
*/
void compcol::fill_data (long nsub,long *smadr,long *smri,double *sma)
{
  long i;

  //  number of columns
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
  
  //  row indices
  if (ri!=NULL){
    delete [] ri;
  }
  ri = new long [negm];
  
  for (i=0;i<negm;i++){
    ri[i]=smri[i];
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
   @param smcr - pointer to object of the class compcol which stores submatrix
   
   JK, 19.2.2007
*/
/*
void compcol::select_submatrix (long *li,long nsub,compcol *smcr)
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
  
  
  //  submatrix is stored in new object of the class compcol
  smcr->fill_data (nsub,smadr,smci,sma);
  
  delete [] sma;
  delete [] smci;
  delete [] smadr;
  delete [] aux;
  aux = NULL;
}
*/

/**
   function selects submatrix from the original %matrix stored in the compressed
   row storage scheme, the submatrix is defined by indices
   
   @param li - list of indices of submatrix
   @param nsub - number of rows and columns of submatrix (it is also length of the array li)
   @param smdm - pointer to object of the class densemat which stores submatrix
   
   JK, 14.3.2007
*/
/*
void compcol::select_submatrix (long *li,long nsub,densemat *smdm)
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
*/



/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long compcol::give_negm ()
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
double compcol::estim_spect_radius ()
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
