#include "scc.h"
#include "diagmat.h"
#include "globalg.h"
#include <limits.h>
#include <math.h>
#include <string.h>

#ifdef INC_CHMD
 // macros needed by CHOLMOD
 #define TRUE  1
 #define FALSE 0
#endif

symcompcol::symcompcol(void)
{
  limit = 0.0;
  adra=NULL;
  ri=NULL;
  a=NULL;
  aux=NULL;
  decompid=0;
  adr = NULL;
#ifdef INC_CHMD  
  chmd_spr = NULL;
  chmd_fact = NULL;
#endif
}



symcompcol::~symcompcol(void)
{
#ifdef INC_CHMD
  // deallocation of memory of CHOLMOD structures is not needed
  // they contain references to arrays managed by SIFEL
  cholmod_l_free_sparse(&chmd_spr, &chmd_com);

  // this structure is created by CHOLMOD -> deallocate memory
  if (chmd_fact)
    cholmod_l_free_factor (&chmd_fact, &chmd_com);
  cholmod_l_finish(&chmd_com);
#else
  delete [] adr;
  delete [] ri;
  delete [] a;
#endif  
  delete [] adra;
  delete [] aux;
}



/**
  The function allocates array containing addresses of the first entries in columns (adra)
  and the last column entry (adra)
  
   
  @param m[in] - number of unknowns in solved problem

  @return The function does not return anything but allocates arrays adr and adra.

  Created by JK
*/
void symcompcol::allocadr(long m)
{
  n=m;

  adr  = new long [n+1];
  adra = new long [n+1];
  
  memset (adr,0,(n+1)*sizeof(*adr));
  memset (adra,0,(n+1)*sizeof(*adra));
}



/**
  The function evaluates number of contributions to the %matrix from one element
   
  @param cn[in] - array containing code numbers of the element

  Created by JK
*/
void symcompcol::numcontr_elem(const ivector &cn)
{
  long i,j,ii,jj;

  long ndofe = cn.n;
  for (i=0; i<ndofe; i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<ndofe;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)  continue;
      adr[ii]++;
    }
  }
}



/**
  The function evaluates number of contributions to the %matrix due to elements in the array pelems.
   
  @param top[in] - pointer to general topology
   
  @return The function does not return anything but it stores results in the array adr.   

  JK, 22.6.2001
*/
void symcompcol::numcontr(gtopology *top)
{
  long i, ne, ndofe;
  ivector cn;

  ne=top->ne;
  for (i=0;i<ne;i++){
    ndofe=top->give_ndofe (i);
    reallocv(RSTCKIVEC(ndofe, cn));
    top->give_code_numbers(i, cn.a);
    numcontr_elem (cn);
  }
}



/**
  The function computes addresses of the first entries in the %matrix columns.
  Components of array adr has to be initialized by the number of contributions 
  from all considered elements in the problem solved.

  @returns The function stores computed addresses in the array adr and adra.   
*/
void symcompcol::addresses(void)
{
  long i,j,k;
  
  j=adr[0];  adr[0]=0L;
  for (i=1L; i<n+1L; i++){
    k=adr[i];
    adr[i]=j;
    j+=k;
  }
  
  for (i=0L; i<n+1L; i++){
    adra[i]=adr[i];
  }
  
  aux = new long[adr[n]];
  memset(aux, 0, adr[n]*sizeof(*aux));
}



/**
  The function fills auxiliary array by contributions from one element.
   
  @param cn[in] - array of element code numbers
   
  @return The function stores in the array aux row indices where
          the contributions from the given element will be stored.

  Created by JK
*/
void symcompcol::fillarray_elem(const ivector &cn)
{
  long i,j,ii,jj;

  long ndofe = cn.n;
  
  for (i=0L; i<ndofe; i++){
    ii = cn[i]-1L;
    if (ii<0)  continue;
    for (j=0L; j<ndofe; j++){
      jj = cn[j]-1L;
      if (jj < 0L)  continue;
      if (ii < jj)  continue;
      aux[adra[ii]] = jj;
      adra[ii]++;
    }
  }
}



/**
  The function fills auxiliary array with row indices.
   
  @param top[in] - pointer to general topology
   
  @return The function stores in the array aux row indices where
          the contributions all elements will be stored.
   
  Created by JK, 22.6.2001
*/
void symcompcol::fillarray(gtopology *top)
{
  long i, ne, ndofe;
  ivector cn;
  ne=top->ne;
  for (i=0L; i<ne; i++){
    ndofe=top->give_ndofe (i);
    reallocv(RSTCKIVEC(ndofe, cn));
    top->give_code_numbers(i, cn.a);
    fillarray_elem(cn);
  }
}



/**
   function sortes array aux
   function also allocates array for matrix storage
   
   JK, 8.7.2001
*/
void symcompcol::sort_and_colindex (void)
{
  long i,j,k,ii,jj,lj,uj,min,prev;

  // aux contains DOF code numbers of particular contributions from elements
  // adr contains addresses of column beginnings in the array aux
  // adra contains addresses of column ends in the array aux (adra[i]-adr[i] = line length of array aux)

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

#ifdef INC_CHMD    
  // initialize CHOLMOD sparse matrix
  if (cholmod_l_start(&chmd_com) == 0){
    print_err("CHOLMOD cannot be started", __FILE__, __LINE__, __func__);
    abort();
  }
  chmd_com.supernodal = -1;
  chmd_com.nthreads_max = Numth;
  chmd_spr = cholmod_l_allocate_sparse(n, // nrow
                                       n, // ncol
                                       negm, // nzmax
                                       TRUE, // sorted columns
                                       TRUE, // array nz is ignored (packed format)
                                       1,    // stype=1 -> upper triangular part is stored
                                       CHOLMOD_REAL, // xtype
                                       &chmd_com);
  fprintf(stdout, "\nSetup of sparse matrix:\n"
          "stype=%d, itype=%d, xtype=%d, dtype=%d, sorted=%d, packed=%d\n",
          chmd_spr->stype, chmd_spr->itype, chmd_spr->xtype, chmd_spr->dtype,
          chmd_spr->sorted, chmd_spr->packed);

  memcpy(chmd_spr->p, adr, sizeof(*adr)*(n+1));
  delete [] adr;
  adr = (decltype(adr))chmd_spr->p;
  
  memset(chmd_spr->i, 0,   sizeof(*ri)*(negm));
  ri = (decltype(ri))chmd_spr->i;
  memset(chmd_spr->x, 0,   sizeof(*a)*(negm));
  a  = (decltype(a))chmd_spr->x;
  
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
  aux = NULL;
#else
  //  array containing resulting row indices
  if (ri)
    delete [] ri;

  ri = new long [negm];
  memset (ri,0,negm*sizeof(*ri));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(*a));


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
  aux = NULL;
#endif
}



/**
  The function allocates CHOLMOD sparse matrix object and initializes
  its array of column beginnings and row indices.

  @retval 0 - on success
  @retval 1 - in the case of an error

  JK, 8.7.2001
*/
void symcompcol::sort_and_colindex_tko (void)
{
  long i,j,ii,nd,lj,uj;

  // aux contains DOF code numbers of particular contributions from elements
  // adr contains addresses of column beginnings in the array aux
  // adra contains addresses of column ends in the array aux (adra[i]-adr[i]) = column length of array aux)

  long *udof = new long[n]; // array of used dof indicators
  memset(udof, 0, sizeof(*udof)*n);

  for (i=0L; i<n; i++)
  {
    // remove duplicity of used DOFs
    ii = adr[i];
    lj = adr[i];
    uj = adra[i];
    for (j=lj; j<uj; j++)
    {
      if (udof[aux[j]] == 0L)
      {
        udof[aux[j]] = 1;
        aux[ii] = aux[j];
        ii++;
      }
    }
    uj = adra[i] = ii;
    nd = uj - lj;


    if (double(nd)*nd < double(n))
    {
      for (j=lj; j<uj; j++)
        udof[aux[j]] = 0L;
      // sort used DOFs
      shell_sort(aux+adr[i], nd);  
    }
    else
    {
      ii = lj; // = adr[i]
      for (j=0; j<n; j++)
      {
        if (udof[j] == 0L)
          continue;

        // udof[j] is nonzero
        aux[ii] = j;
        ii++;
        udof[j] = 0L;
      }
    }
  }

  delete [] udof;

  j=0;
  for (i=0;i<n;i++){
    j+=adra[i]-adr[i];
  }

  //  number of non-zero elements in matrix
  negm=j;

#ifdef INC_CHMD    
  // initialize CHOLMOD sparse matrix
  if (cholmod_l_start(&chmd_com) == 0){
    print_err("CHOLMOD cannot be started", __FILE__, __LINE__, __func__);
    abort();
  }
  // chmd_com.supernodal = CHOLMOD_SIMPLICIAL;
  // chmd_com.supernodal = CHOLMOD_AUTO;
  chmd_com.supernodal = CHOLMOD_SUPERNODAL;
  
  chmd_spr = cholmod_l_allocate_sparse(n, // nrow
                                       n, // ncol
                                       negm, // nzmax
                                       TRUE, // sorted columns
                                       TRUE, // array nz is ignored (packed format)
                                       1,    // stype=1 -> upper triangular part is stored
                                       CHOLMOD_REAL, // xtype
                                       &chmd_com);
  if (chmd_spr == NULL){
    print_err("CHOLMOD cannot be started", __FILE__, __LINE__, __func__);
    abort();
  }
  fprintf(stdout, "\nSetup of sparse matrix:\n"
          "stype=%d, itype=%d, xtype=%d, dtype=%d, sorted=%d, packed=%d\n",
          chmd_spr->stype, chmd_spr->itype, chmd_spr->xtype, chmd_spr->dtype,
          chmd_spr->sorted, chmd_spr->packed);
  memcpy(chmd_spr->p, adr, sizeof(*adr)*(n+1));
  delete [] adr;
  adr = (decltype(adr))chmd_spr->p;
  
  memset(chmd_spr->i, 0,   sizeof(*ri)*(negm));
  ri = (decltype(ri))chmd_spr->i;
  memset(chmd_spr->x, 0,   sizeof(*a)*(negm));
  a  = (decltype(a))chmd_spr->x;
  
  ii = 0L;   lj = 0L;   adr[0] = 0L;
  for (i=0L; i<n; i++){
    uj=adra[i];
    for (j=lj; j<uj; j++){
      ri[ii] = aux[j];
      ii++;
    }
    lj = adr[i+1];
    adr[i+1] = ii;
  }
  
  for (i=0L; i<=n; i++){
    adra[i] = adr[i];
  }

  delete [] aux;
  aux = NULL;
#else
  //  array containing resulting column indices
  if (ri)   delete [] ri;

  ri = new long [negm];
  memset (ri, 0, negm*sizeof(*ri));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double[negm];
  memset (a, 0, negm*sizeof(*a));

  ii=0L;  lj=0L;  adr[0]=0L;
  for (i=0L; i<n; i++){
    uj = adra[i];
    for(j=lj; j<uj; j++){
      ri[ii]=aux[j];
      ii++;
    }
    lj=adr[i+1];
    adr[i+1]=ii;
  }
  
  for (i=0;i<=n;i++){
    adra[i]=adr[i];
  }

  delete [] aux;
  aux = NULL;
#endif
}



/**
  The function localizes local %matrix b to the global %matrix.
  %Matrix b is stored as dense %matrix.
   
  @param b[in] - local %matrix in columns ordering
  @param cn[in] - array containing code numbers of element

  @return The function does not return anything but it modifies 
          matrix components.   

  Created by JK
*/
void symcompcol::localize(matrix &b, long *cn)
{
  long i,j,k,ii,jj,kk,lk,uk,ll;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    ll=0;  lk=adr[ii];  uk=adr[ii+1];  kk=lk;
    for (j=0;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0) continue;
      if (ii<jj)  continue;
      if (ll<jj){
	ll=jj;
	for (k=kk;k<uk;k++){
	  if (ri[k]!=jj)  continue;
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
	  if (ri[k]!=jj)  continue;
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
void symcompcol::localized(double *b, long *cn, long nc)
{
  long i,j,k,ii,jj,kk,lk,uk,ll;
  
  for (i=0;i<nc;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    ll=0;  lk=adr[ii];  uk=adr[ii+1];  kk=lk;
    for (j=0;j<nc;j++){
      jj=cn[j]-1;
      if (jj<0) continue;
      if (ii<jj)  continue;
      if (ll<jj){
	ll=jj;
	for (k=kk;k<uk;k++){
	  if (ri[k]!=jj)  continue;
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
	  if (ri[k]!=jj)  continue;
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
long symcompcol::minimize(double limit)
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
  adra = NULL;
  return n1;
}



/**
  The function initiates a %matrix stored in the symmetric compressed 
  column storage format. CHOLMOD sparse structure is exploited for the %matrix storage.

  @param top[in]   - pointer to general topology
  @param ndof[in]  - number of rows/columns of the %matrix
  @param mespr[in] - indicator of message printing
   
  Created by TKo 04.2023 according to JK
*/
long symcompcol::initiate(gtopology *top, long ndof, long mespr)
{
  if (status() == NULL){
    decompid = 0;
    allocadr(ndof);    
    numcontr(top);       // adr[i] obsahuje pocet prispevku z prvku do i-teho stupne volnosti
    addresses();         // adr[i] obsahuje zacatek radku v poli aux, pole aux je vynulovane, adra[i]=adr[i]
    fillarray(top);      // naplni se pole aux kodovymi cisly jednotlivych prispevku od prvku, adra[i]=index konce radku s kodovymi cisly+1(tj. zacatek nasl radku s kodovymi cisly)
    //sort_and_colindex(); // used for common problems
    sort_and_colindex_tko(); // used for the stress approach in homogenization problem

#ifdef INC_CHMD
    /*
    // initialize CHOLMOD sparse matrix
    chmd_spr.nrow  = n;
    chmd_spr.ncol  = n;
    chmd_spr.nzmax = negm;
    chmd_spr.p  = adr;
    chmd_spr.i  = ri;
    chmd_spr.nz = NULL; // for unpacked matrices - number of nonzeros in each column
    chmd_spr.x  = a;     // array of matrix elements in copressed column format
    chmd_spr.z  = NULL;
    //    chmd_spr.stype  = -1; // lower trianglular part is stored
    chmd_spr.stype  = 1; // upper trianglular part is stored
    chmd_spr.itype  = CHOLMOD_LONG; // arrays chmd_spr.i and chmd_spr.p are type of long
    chmd_spr.xtype  = CHOLMOD_REAL; // real number algebra type
    chmd_spr.dtype  = CHOLMOD_DOUBLE; // matrix entries x are type of double
    chmd_spr.sorted = TRUE; // sorted columns
    chmd_spr.packed = TRUE; // array nz is ignored
    fprintf(stdout, "\nSetup of sparse matrix:\n"
                     "stype=%d, itype=%d, xtype=%d, dtype=%d, sorted=%d, packed=%d\n",
            chmd_spr.stype, chmd_spr.itype, chmd_spr.xtype, chmd_spr.dtype, chmd_spr.sorted, chmd_spr.packed);
    if (cholmod_start(&chmd_com) == 0){
      print_err("CHOLMOD cannot be started", __FILE__, __LINE__, __func__);
      abort();
    }
    // allocate workspace of the problem
    cholmod_allocate_work(0, n+n, 0, &chmd_com);*/
    // make symbolic factorization
    if (mespr == 1)
      fprintf(stdout, " Starting matrix analysis ... ");
    //    chmd_com.supernodal = 2;
    chmd_fact = cholmod_l_analyze(chmd_spr, &chmd_com);
    if (mespr == 1){
      fprintf(stdout, " - reordering methods used (%d): ", chmd_com.nmethods);
      for (int i=0; i<chmd_com.nmethods; i++)  fprintf(stdout, "%s,", chmd_ord_method(chmd_com.method[i].ordering));
      fprintf(stdout, "\b\n");
    }
    if (chmd_fact == NULL){
      print_err("CHOLMOD matrix cannot be analyzed", __FILE__, __LINE__, __func__);
      abort();
    }
    if (mespr == 1){
      fprintf(stdout, " - selected reordering method: %s\n", chmd_ord_method(chmd_com.selected));      
      switch(chmd_com.supernodal){
        case 1:
          fprintf(stdout, " - automatic switch between supernodal and simplicial factorization\n");
          break;
        case 2:
          fprintf(stdout, " - supernodal factorization will be done\n");
          break;
        default:
          fprintf(stdout, " - simplicial factorization will be done\n");
          break;
      }
      fprintf(stdout, " - maximum number of threads used: %d\n", Numth);
    }
      
#endif    
  }
  else{
    nullmat ();
  }
  
  fprintf(stdout, "\n Number of non-zero matrix entries %ld", negm);
  return 0L;
}



/**
  Funtion returns name of reordering CHOLMOD methods for the given identifier.

  @param i[in] - method identifier

  @return The funtion returns string with the method name

  Created by Tomas Koudelka
*/
const char*symcompcol::chmd_ord_method(int i)
{
  switch (i){
    case 0:
      return "UserDef";
    case 1:
      return "ApproxMinDeg";
    case 2:
      return "METIS";
    case 3:
      return "NestDisect";
    case 4:
      return "Natural";
    default:
      return "Other";
  }
  return NULL;
}



/**
  The function returns required %matrix entry.
   
  @param ir[in] - row index of the required %matrix entry
  @param ic[in] - column index of the required %matrix entry

  @return It returns matrix entry at with coordinates (ir,ic)
   
  JK, 24.7.2005
*/
double symcompcol::give_entry(long ir, long ic)
{
  long i;
  double e = 0.0;
  
  if (ir>ic){
    i=ic;  ic=ir;  ir=i;    
  }
  for (i=adr[ic]; i<adr[ic+1]; i++){
    if (ir==ri[i]){
      e=a[i];
      break;
    }
  }

  return e;
}



/**
   function adds required matrix entry
   
   @param e[in] - matrix entry
   @param ir[in] - row index
   @param ic[in] - column index
   
   JK, 24.7.2005
*/
void symcompcol::add_entry(double e, long ir, long ic)
{
  long i;
  
  if (ir>=ic){
    i=ic;  ic=ir;  ir=i;    
  }
  for (i=adr[ic];i<adr[ic+1];i++){
    if (ir==ri[i]){
      a[i] += e;
      break;
    }
  }
}



#ifdef INC_CHMD  
void symcompcol::factorize(void)
{
  cholmod_l_factorize(chmd_spr, chmd_fact, &chmd_com);
  decompid = 1;
}
#else
void symcompcol::factorize(void)
{
  print_err("\n Support for CHOLMOD solver was not defined,"
            "\n Use different storage/solver type, e.g. skyline with LDL factorization.\n", __FILE__, __LINE__, __func__);  
  abort();
}
#endif



#ifdef INC_CHMD  
/**
  The function solves system of linear algebraic equations
  by external slover CHOLMOD from SuiteSprase library, 
  %matrix is stored as compressed column format, left and 
  right-hand side vectors are stored in dense storage format.
   
  @param x[out] - %vector of unknowns (left-hand side), 
                  it must be allocated at least to n components
  @param y[in/out] - %vector of right hand side - may be modified (right-hand side),
                  it must be allocated at least to n components

  @return The function stores solution in the array x.

  Created by TKo, 29.3.2023
*/
void symcompcol::solve(double *x, double *y)
{
  if (decompid == 0)   factorize();

  chmd_rhs.nrow  = n; // number of matrix rows
  chmd_rhs.ncol  = 1; // number of matrix columns
  chmd_rhs.nzmax = n; // total number of the densmatrix components
  chmd_rhs.d     = n; // leading dimension for the accessing element x(i,j) = chmd_rhs.x[i+j*chmd_rhs.d] -> rhs vectors are stored by rows
  chmd_rhs.x     = y; // array of matrix components
  chmd_rhs.z     = NULL; // array for imaginary components in the case of complex number algebra
  chmd_rhs.xtype = CHOLMOD_REAL;   // flag for real or complex number algebra
  chmd_rhs.dtype = CHOLMOD_DOUBLE; // flag for type of matrix components (double vs float)

  // solve equation system Ax = B, i.e. type CHOLMOD_A
  fflush(stdout);
  //fprintf(stdout, "\n Solution started ...");
  chmd_lhs = cholmod_l_solve(CHOLMOD_A, chmd_fact, &chmd_rhs, &chmd_com);
  //fprintf(stdout, " OK\n");
  memcpy(x, chmd_lhs->x, n*sizeof(*x));
  cholmod_l_free_dense(&chmd_lhs, &chmd_com);

  chmd_rhs.x = NULL;
  decompid = 1;
}
#else
void symcompcol::solve(double *x, double *y)
{
  print_err("\n Support for CHOLMOD solver was not defined,"
            "\n Use different storage/solver type, e.g. scr with conjugated gradients\n", __FILE__, __LINE__, __func__);  
  abort();
}
#endif



/**
   function multiplies %matrix by %vector
   
   @param b - array containing %vector b
   @param c - array containing resulting %vector c = A.b
   
   JK
*/
void symcompcol::mxv_scc(double *b, double *c)
{
  long i,j,ii,lj,uj;
  double s,d;
  
  nullv(c, n);
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];
    s=0.0;  d=b[i];
    for (j=lj;j<uj;j++){
      ii=ri[j];
      s+=a[j]*b[ii];
      c[ii]+=a[j]*d;
    }
    c[i]=s;
  }
}



/**
   function adds multiplied %matrix stored in scc by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param scc - another symmetric compressed columns storage
   
   JK
*/
void symcompcol::addmat_scc(double c, symcompcol &scc)
{
  long i;
  for (i=0;i<negm;i++){
    a[i] += c*scc.a[i];
  }
}



/**
  The function adds diagonal %matrix multiplied by the coefficient c to actual %matrix.
  
  @param[in] c - multiplicative coefficient
  @param[in] d - diagonal %matrix summand  storage
   
  TKo
*/
void symcompcol::addmat_diag (double c, diagmat &d)
{
  long i, ii, ncomp;
  bool diag_found;
  
  for (i=0; i<n; i++){
    ncomp = adr[i+1] - adr[i];
    diag_found = false;
    ii = adr[i+1]-1; // guessed index of diagonal entry
    if (ncomp && (ri[ii] == i)){
      a[ii] += c*d.a[i];
      diag_found = true;
    }
    if ((diag_found == false) && (d.a[i] != 0.0)){
      print_err("diagonal element on %ld row was not stored in matrix stored in the symmetric compressed column format",
                __FILE__, __LINE__, __func__, i+1);
      abort();
    }
  }
}



/**
   function multiplies %matrix by coefficient c
   
   @param c - multiplicative coefficient
   
   JK
*/
void symcompcol::scalmat_scc(double c)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]*=c;
  }
}



/**
  The function prints matrix in the compressed column storage.
   
  @param out[in/out] - pointer to the opened text file

  @return The function deos not return anything, but prints the %matrix content to the text file.
   
  Created by JK, 10.10.2003
*/
void symcompcol::printmat(FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n\nMatrix in the symmetric compressed column format:\n");
  fprintf (out,"\n Number of rows/columns    %ld",n);
  fprintf (out,"\n Number of stored entries  %ld",negm);
  fprintf (out,"\n\n Array ADR:\n");
  for (i=0;i<n+1;i++){
    fprintf (out,"\n%ld",adr[i]);
  }
  fprintf (out,"\n\n CI  A:\n\n");
  for (i=0;i<n;i++){
    fprintf(out, "Column %ld:\n", i+1);
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"%8ld   % 16.12e\n",ri[j],a[j]);
    }
  }
  fprintf (out,"\n\n");
    
  fprintf(out, "\nDense format generated from the symmetric compressed column format:\n");
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
  The function prints diagonal entries of %matrix stored in the symmetric compressed column 
  storage format.
   
  @param out[in/out] - pointer to the opned text file

  @return The function deos not return anything, but prints the %matrix diagonal to the text file.
   
  Created by JK, 17.3.2007
*/
void symcompcol::printdiag(FILE *out)
{
  fprintf (out,"\n\n diagonal entries of matrix stored in the symmetric compressed column storage scheme\n");
  
  for (long i=0; i<n; i++){
    fprintf(out,"%e\n", give_entry(i,i));
  }
}



/**
   function copies data to additional storage
   
   @param scc - additional storage (it is filled)
   
   JK, 29.9.2006
*/
void symcompcol::copy_scc(symcompcol &scc)
{
  long i;

  scc.n=n;
  scc.negm=negm;
  
  if (scc.adr!=NULL){
    delete [] scc.adr;
  }
  scc.adr=new long [n+1];

  if (scc.ri!=NULL){
    delete [] scc.ri;
  }
  scc.ri=new long [negm];
  
  if (scc.a!=NULL){
    delete [] scc.a;
  }
  scc.a=new double [negm];

  for (i=0;i<=n;i++){
    scc.adr[i]=adr[i];
  }
  for (i=0;i<negm;i++){
    scc.a[i]=a[i];
    scc.ri[i]=ri[i];
  }
}



/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long symcompcol::give_negm()
{
  return negm;
}
