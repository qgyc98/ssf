#include "scr.h"
#include "gmatrix.h"
#include <limits.h>
#include <math.h>
#include <string.h>


symcomprow::symcomprow (void)
{
  limit = 0.0;
  adra=NULL;
  ci=NULL;
  a=NULL;
  aux=NULL;
  incdec=NULL;
  decompid=0;
  adr = NULL;
}



symcomprow::~symcomprow (void)
{
  delete [] adr;
  delete [] ci;
  delete [] a;
  delete [] incdec;
  delete [] adra;
  delete [] aux;
}



/**
   function returns indicator of decomposition (factorization)
   
   JK
*/
long symcomprow::decomp ()
{
  return decompid;
}



/**
   function changes indicator of decomposition (factorization)
   
   JK
*/
void symcomprow::changedecomp ()
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
void symcomprow::allocadr (long m)
{
  n=m;

  adr  = new long [n+1];
  adra = new long [n+1];
  
  memset (adr,0,(n+1)*sizeof(*adr));
  memset (adra,0,(n+1)*sizeof(*adra));
}



/**
   function returns status of array s
   
   JK
*/
double* symcomprow::status ()
{
  return a;
}



/**
   function fills array a by zeros
   
   JK
*/
void symcomprow::nullmat ()
{
  memset(a, 0, negm*sizeof(*a));
}



/**
   function evaluates number of contributions to the %matrix from one element
   
   @param cn - array containing code numbers of the element
   @param ndofe - number of DOFs on the element
   
   JK, 22.6.2001
*/
void symcomprow::numcontr_elem (long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){
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
   function evaluates number of contributions to the %matrix
   
   @param top - pointer to general topology
   
   JK, 22.6.2001
*/
void symcomprow::numcontr (gtopology *top)
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
   
   param aux - auxiliary array containing column indexes

   JK, 22.6.2001
*/
void symcomprow::fillarray_elem (long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<ndofe;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)  continue;
      aux[adra[ii]]=jj;
      adra[ii]++;
    }
  }
}



/**
   function fills auxiliary array
   
   @param top - pointer to general topology
   
   param aux - auxiliary array containing column indices
   
   JK, 22.6.2001
*/
void symcomprow::fillarray (gtopology *top)
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
void symcomprow::addresses (void)
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
   function also allocates array for matrix storage
   
   JK, 8.7.2001
*/
void symcomprow::sort_and_colindex (void)
{
  long i,j,k,ii,jj,lj,uj,min,prev;

  // aux contains DOF code numbers of particular contributions from elements
  // adr contains addresses of row beginnings in the array aux
  // adra contains addresses of row ends in the array aux (adra[i]-adr[i] = line length of array aux)

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

  //  array containing resulting column indices
  delete [] ci;

  ci = new long [negm];
  memset (ci,0,negm*sizeof(*ci));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(*a));


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
  aux = NULL;
  
}



/**
   function sortes array aux
   function also allocates array for matrix storage
   
   TKo according JK, 2.2.2015
*/
void symcomprow::sort_and_colindex_tko (void)
{
  long i,j,ii,nd,lj,uj;

  // aux contains DOF code numbers of particular contributions from elements
  // adr contains addresses of row beginnings in the array aux
  // adra contains addresses of row ends in the array aux (adra[i]-adr[i] = line length of array aux)

  long *udof = new long[n]; // array of used dof indicators
  memset(udof, 0, sizeof(*udof)*n);

  for (i=0;i<n;i++)
  {
    // remove duplicity of used DOFs
    ii = adr[i];
    lj = adr[i];
    uj = adra[i];
    for (j=lj; j<uj; j++)
    {
      if (udof[aux[j]] == 0)
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
        udof[aux[j]] = 0;
      // sort used DOFs
      shell_sort(aux+adr[i], nd);  
    }
    else
    {
      ii = lj; // = adr[i]
      for (j=0; j<n; j++)
      {
        if (udof[j] == 0)
          continue;

        // udof[j] is nonzero
        aux[ii] = j;
        ii++;
        udof[j] = 0;
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

  //  array containing resulting column indices
  delete [] ci;

  ci = new long [negm];
  memset (ci,0,negm*sizeof(*ci));
  
  //  array containing non-zero entries of the matrix
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(*a));


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
  aux = NULL;  
}



/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as dense %matrix
   
   @param b - local %matrix in row ordering
   @param cn - array containing code numbers of element
   
   JK, 25.6.2001
*/
void symcomprow::localize (matrix &b,long *cn)
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
void symcomprow::localized (double *b,long *cn,long nc)
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
long symcomprow::minimize (double limit)
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
  adra = NULL;
  return n1;
}



/**
   function initiates symmetric compressed rows storage
   
   @param top - pointer to general topology
   @param ndof - number of rows/columns of the %matrix
   @param mespr - indicator of message printing
   
   JK
*/
void symcomprow::initiate (gtopology *top,long ndof,long mespr)
{
  if (status()==NULL){
    allocadr (ndof);    
    numcontr (top); // adr[i] obsahuje pocet prispevku z prvku do i-teho stupne volnosti
    addresses (); // adr[i] obsahuje zacatek radku v poli aux, pole aux je vynulovane, adra[i]=adr[i]
    fillarray (top); // naplni se pole aux kodovymi cisly jednotlivych prispevku od prvku, adra[i]=index konce radku s kodovymi cisly+1(tj. zacatek nasl radku s kodovymi cisly)
    sort_and_colindex_tko (); // use for the stress approach in homogenization problem
    //sort_and_colindex ();
  }
  else{
    nullmat ();
  }
  
  if (mespr==1)  fprintf (stdout,"\n number of matrix entries (kontrola) %ld",negm);
}



/**
   function returns required %matrix entry
   
   @param ir,ic - row and column indices
   
   JK, 24.7.2005
*/
double symcomprow::give_entry (long ir,long ic)
{
  long i;
  double e=0.0;
  
  if (ir>ic){
    for (i=adr[ir];i<adr[ir+1];i++){
      if (ic==ci[i]){
	e=a[i];
	break;
      }
    }
  }
  else{
    i=ir;
    ir=ic;
    ic=i;
    for (i=adr[ir];i<adr[ir+1];i++){
      if (ic==ci[i]){
	e=a[i];
	break;
      }
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
void symcomprow::add_entry (double e,long ri,long rci)
{
  long i;
  
  if (ri>rci){
    for (i=adr[ri];i<adr[ri+1];i++){
      if (rci==ci[i]){
	a[i]+=e;
	break;
      }
    }
  }
  else{
    i=ri;
    ri=rci;
    rci=i;
    for (i=adr[ri];i<adr[ri+1];i++){
      if (rci==ci[i]){
	a[i]+=e;
	break;
      }
    }
  }
}



/**
   function multiplies %matrix by %vector
   
   @param b - array containing %vector b
   @param c - array containing resulting %vector c = A.b
   
   JK
*/
void symcomprow::mxv_scr (double *b,double *c)
{
  long i,j,ii,lj,uj;
  double s,d;
  
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];
    s=0.0;  d=b[i];
    for (j=lj;j<uj;j++){
      ii=ci[j];
      s+=a[j]*b[ii];
      c[ii]+=a[j]*d;
    }
    c[i]=s;
  }
}



/**
   function adds multiplied %matrix stored in scr by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param scr - another symmetric compressed rows storage
   
   JK
*/
void symcomprow::addmat_scr (double c,symcomprow &scr)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]+=c*scr.a[i];
  }
}

/**
   function adds multiplied %matrix stored in scr by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param scr - another symmetric compressed rows storage
   
   JK, 22. 10. 2023
*/
void symcomprow::addmat_gm (double c,gmatrix &gm)
{
  long i;
  
  switch (gm.ts){
  case diag_mat:{
    if (gm.n != n){
      print_err("inconsistent number of rows/columns in matrices",__FILE__,__LINE__,__func__);
      abort ();
    }
    for (i=0;i<n;i++){
      a[adr[i+1]-1]+=c*gm.diagm->a[i];
    }
    break;
  }
    //case dense_matrix:{
    //dm -> addmat_dm (a,*gm.dm);
    //break;
    //}
    //case skyline_matrix:{
    //sky -> addmat_sky (a,*gm.sky);
    //break;
    //}
    //case double_skyline:{
    //dsky -> addmat_dsky (a,*gm.dsky);
    //break;
    //}
    //case compressed_rows:{
    //cr -> addmat_cr (a,*gm.cr);
    //break;
    //}
    // case symm_comp_rows:{
    //scr->addmat_gm (a,gm);
    //break;
    //}
    //case element_matrices:{
    //em -> addmat_em (a,*gm.em);
    //break;
    //}
    //case spdirect_stor_scr:{
    //scr -> addmat_scr (a,*gm.scr);
    //sdirect -> AddNumbers (a,gm.sdirect->GetSparseMatrix());
    //break;
    //}
    //case spdirect_stor_cr:{
    //cr -> addmat_cr (a,*gm.cr);
    //break;
    //}
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}

/**
   function multiplies %matrix by coefficient c
   
   @param c - multiplicative coefficient
   
   JK
*/
void symcomprow::scalmat_scr (double c)
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
void symcomprow::cg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double nom,denom,nory,alpha,beta;
  double *d,*r,*p;
  
  fprintf (stdout,"NDOF %ld\n",n);
  
  d = new double [n];
  r = new double [n];
  p = new double [n];
  
  //  initial values
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  
  mxv_scr (x,p);

  //for (i=0;i<n;i++){
  //fprintf (stdout," vektor p %6ld  %le\n",i,p[i]);
  //}

  nory=0.0;  nom=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=y[i]-p[i];
    nom+=r[i]*r[i];
    d[i]=r[i];
  }
  
  fprintf (stdout," SCR nory %le\n",nory);
  fprintf (stdout," SCR nom  %le\n",nom);

  if (nory<zero){
    print_warning("norm of right hand side in conjugate gradient method is smaller than %e", __FILE__, __LINE__, __func__,zero);
    ares=nory;  ani=0;
    for (i=0;i<n;i++)
      x[i]=0.0;
    
    delete [] p;  delete [] r;  delete [] d;
    
    return;
  }
  
  
  //  iteration loop
  for (i=0;i<ni;i++){

    //  new coefficient alpha
    mxv_scr (d,p);

    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      print_warning("there is zero denominator in alpha computation in conjugate gradient method", __FILE__, __LINE__, __func__);
      break;
    }
    
    alpha = nom/denom;
    
    //  new approximation of x and r
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
    }
    
    denom=nom;

    nom = ss (r,r,n);

    if (i%10==0)
      fprintf (stdout,"\n iteration  %6ld  norres  %le    norres/norrhs %e",i,nom,sqrt(nom/nory));
    //fprintf (stdout,"\n iteration   %ld  norres/norrhs %e",i,sqrt(nom/nory));
    
    
    if (sqrt(nom/nory)<res){//this printing added by TKr 30/07/2008
      fprintf (stdout,"\n\n Last iteration step   %ld  norres/norrhs %e\n\n",i,sqrt(nom/nory));
      break;
    }
    //if (nom<res)  break;
    if (fabs(nom)<zero){//this printing added by TKr 30/07/2008
      fprintf (stdout,"\n\n Last iteration step   %ld  norres/norrhs %e\n\n",i,sqrt(nom/nory));
      break;
    }

    beta = nom/denom;
    
    //  new vector of direction
    for (j=0;j<n;j++){
      d[j]=beta*d[j]+r[j];
    }
  }
  
  ani=i;  ares=nom;
  
  delete [] p;  delete [] r;  delete [] d;
}



/**
   function computes preconditioned %vector via Jacobi method
   new %vector is computed from equation A.x = y, where
   A is diagonal %matrix and contains diagonal entries of system %matrix
   
   @param x - %vector of unknowns
   @param y - %vector of right hand side

   JK, 2.12.2001
*/
void symcomprow::prec_diag_scr (double *x,double *y)
{
  long i;
  
  for (i=0;i<n;i++){
    x[i]=y[i]/a[adr[i+1]-1];
  }
}



/**
   function computes preconditioned %vector via SSOR method
   new %vector is computed approximately from equation A.x = y
   computation is done like one step of SSOR method
   
   @param x - %vector of unknowns
   @param y - %vector of right hand side
   @param omega - relaxation parameter
   
   JK, 2.12.2001
*/
void symcomprow::prec_ssor_scr (double *x,double *y,double omega)
{
  long i,j,uj;
  double s,*v;
  
  v = new double [n];
  memset (v,0,n*sizeof(double));

  for (i=0;i<n;i++){
    uj=adr[i+1]-1;  s=0.0;
    for (j=adr[i];j<uj;j++){
      s+=a[j]*v[ci[j]];
    }
    v[i]=(y[i]-s)/a[uj]*omega;
  }
  
  for (i=0;i<n;i++){
    v[i]=v[i]*a[adr[i+1]-1]/omega*(2.0-omega);
  }
  
  for (i=n-1;i>-1;i--){
    x[i]=v[i]*omega/a[adr[i+1]-1];
    uj=adr[i+1]-1;  s=x[i];
    for (j=adr[i];j<uj;j++){
      v[ci[j]]-=a[j]*s;
    }
  }
  
  delete [] v;
}



/**
   function computes multiplication of two %matrix rows in scr storage
   
   @param m,o - row numbers, equality m>=o must be satisfied
   
   JK, 2.12.2001
*/
double symcomprow::rows_mult_ldl_scr (long m,long o)
{
  long i,j,ii,jj;
  double s;
  
  if (m<o){
    print_err("wrong row numbers", __FILE__, __LINE__, __func__);
  }
  
  s=0.0;  ii=adr[m];  jj=adr[o];
  i=ci[ii];  j=ci[jj];
  while (j<o){
    if (i==j) {  s+=incdec[ii]*incdec[jj]*incdec[adr[j+1]-1];  ii++;  jj++;  }
    if (i<j)  {  ii++;  }
    if (i>j)  {  jj++;  }
    i=ci[ii];  j=ci[jj];
  }
  
  return s;
}



/**
   function decomposes %matrix in symmetric compressed rows storage
   by incomplete decomposition (factorization)
   
   @param gamma - weight coefficient
   
   JK, 2.12.2001
*/
void symcomprow::incomplete_ldl (double gamma)
{
  long i,j,ii,lj,uj;
  double s;
  
  //  allocation of incdec array
  incdec = new double [negm];
  
  for (i=0;i<n;i++){
    //  offdiagonal elements evaluation
    ii=adr[i];  lj=ci[ii];  uj=adr[i+1]-1;  incdec[uj]=a[uj];
    for (j=lj;j<i;j++){
      s = rows_mult_ldl_scr (i,j);
      if (j==ci[ii]) {  incdec[ii]=(a[ii]-s)/incdec[adr[j+1]-1];  ii++;  }
      else           {  incdec[uj]-=gamma*s;  }
    }

    //  diagonal element evaluation
    lj=adr[i];  s=0.0;
    for (j=lj;j<uj;j++){
      s+=incdec[j]*incdec[j]*incdec[adr[ci[j]+1]-1];
    }
    incdec[uj]-=s;
  }
  
}



/**
   function computes preconditioned vector via incomplete decomposition
   new %vector is computed from equation A.x = y
   %matrix A is incompletely decomposed in scr storage
   
   @param x - %vector of unknowns
   @param y - %vector of right hand side

   JK, 2.12.2001
*/
void symcomprow::prec_ildl_scr (double *x,double *y)
{
  long i,j,lj,uj;
  double s,*p;
  
  p = new double [n];
  
  //  Lp=y
  for (i=0;i<n;i++){
    s=0.0;  lj=adr[i];  uj=adr[i+1]-1;
    for (j=lj;j<uj;j++){
      s+=incdec[j]*p[ci[j]];
    }
    p[i]=y[i]-s;
  }
  
  //  D^-1
  for (i=0;i<n;i++){
    p[i]/=incdec[adr[i+1]-1];
  }
  
  //  Lx=p
  for (i=n-1;i>-1;i--){
    lj=adr[i];  uj=adr[i+1]-1;
    x[i]=p[i];  s=x[i];
    for (j=lj;j<uj;j++){
      p[ci[j]]-=incdec[j]*s;
    }
  }
  
  delete [] p;
}



/**
   function solves system of linear algebraic equations by
   preconditioned conjugate gradient method
   
   if incomplete decomposition preconditioner is used, 
   %matrix must be incompletely decomposed before this function is called
   
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
	
   JK, 2.12.2001
*/
void symcomprow::cg_prec (double *x,double *y,long ni,double res,long &ani,double &ares,
			  double zero,long iv,long tprec,double par,ISolver *sdirect)
{
  long i,j;
  double nom,denom,nory,alpha,beta,rr;
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
  mxv_scr (x,p);
  nory=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=y[i]-p[i];
  }
  
  if (nory<zero){
    print_err("right hand side vector has zero norm", __FILE__, __LINE__, __func__);
    ares=nory;  ani=0;
  }
  
  switch (tprec){
  case 1:{   prec_diag_scr (h,r);  break;  }
  case 5:{   prec_ssor_scr (h,r,par);  break; }
  case 10:{  prec_ildl_scr (h,r);  break;  }
  case sparseindec:{
    sdirect->Solve (h,r);
    break;
  }
  default:{
    print_err("wrong preconditioner type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  nom = ss (r,h,n);
  
  for (i=0;i<n;i++){
    d[i]=h[i];
  }
  
  //  iteration loop
  for (i=0;i<ni;i++){

    //  matrix-vector multiplication
    mxv_scr (d,p);
    
    denom = ss (d,p,n);
    if (fabs(denom)<zero){
      print_err("denominator in alpha expression is equal to zero", __FILE__, __LINE__, __func__);
      ares=nom;  ani=i;
      return;
    }
    
    //  new alpha
    alpha = nom/denom;
    
    //  new vectors
    rr=0.0;
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
      rr+=r[j]*r[j];
    }
    
    
    //  preconditioning
    switch (tprec){
    case 1:{   prec_diag_scr (h,r);  break;  }
    case 5:{   prec_ssor_scr (h,r,par);  break; }
    case 10:{  prec_ildl_scr (h,r);  break;  }
    case sparseindec:{
      sdirect->Solve (h,r);
      break;
    }
    default:{
      print_err("wrong preconditioner type is required", __FILE__, __LINE__, __func__);
    }
    }
    
    denom=nom;
    nom = ss (r,h,n);
    beta = nom/denom;
    
    if (i%100==0)
      fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,rr/nory);
    //fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,rr/nory);
    
    if (fabs(rr/nory)<res){
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
   function prints matrix in the compressed rows storage
   
   10.10.2003, JK
*/
void symcomprow::printmat (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n\nMatrix in the symmetric compressed row format:\n");
  fprintf (out,"\n Number of rows/columns    %ld",n);
  fprintf (out,"\n Number of stored entries  %ld",negm);
  fprintf (out,"\n\n Array ADR:\n");
  for (i=0;i<n+1;i++){
    fprintf (out,"\n%ld",adr[i]);
  }
  fprintf (out,"\n\n CI  A:\n\n");
  for (i=0;i<n;i++){
    fprintf(out, "Row %ld:\n", i+1);
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"%8ld   % 16.12le\n",ci[j],a[j]);
    }
  }
  fprintf (out,"\n\n");
    
  fprintf(out, "\nDense format generated from the symmetric compressed row format:\n");
  fprintf(out, "%6c", ' ');
  for (i=0;i<n;i++)
    fprintf(out, " %16ld", i+1);
  fprintf(out, "\n");
  for (i=0;i<n;i++){
    fprintf(out, "%6ld", i+1);
    for (j=0;j<n;j++){
      fprintf(out, "% 16le ", give_entry(i, j));
    }
    fprintf(out, "\n");
  }
  fprintf (out,"\n\n");  
}



/**
   function prints diagonal entries of %matrix stored in the symmetric compressed row storage scheme
   
   @param out - output stream
   
   JK, 17.3.2007
*/
void symcomprow::printdiag (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n diagonal entries of matrix stored in the symmetric compressed row storage scheme\n");
  
  for (i=0;i<n;i++){
    fprintf (out,"%lf\n",a[adr[i+1]-1]);
  }
}



/**
   function copies data to additional storage
   
   @param scr - additional storage (it is filled)
   
   JK, 29.9.2006
*/
void symcomprow::copy_scr (symcomprow &scr)
{
  long i;

  scr.n=n;
  scr.negm=negm;
  
  if (scr.adr!=NULL){
    delete [] scr.adr;
  }
  scr.adr=new long [n+1];

  if (scr.ci!=NULL){
    delete [] scr.ci;
  }
  scr.ci=new long [negm];
  
  if (scr.a!=NULL){
    delete [] scr.a;
  }
  scr.a=new double [negm];

  for (i=0;i<=n;i++){
    scr.adr[i]=adr[i];
  }
  for (i=0;i<negm;i++){
    scr.a[i]=a[i];
    scr.ci[i]=ci[i];
  }
}

/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long symcomprow::give_negm ()
{
  return negm;
}


/**
   Funkce vybere prvky z matice na zaklade pole sdof, kde jsou umisteny prvky,
   ktere se maji vybrat
   Podminka pouziti funkce - pole sdof je serazeno od nejmensiho kodoveho cisla
   po nejvetsi kodove cislo
   
   JB, 11.10.2007 

*/
void symcomprow::select_submatrix(symcomprow *smscr,long nsdof,long *sdof)
{
  
  long i,j,k,m,*auxci;
  long address1,address2;
     
  smscr->n = nsdof;
   
  // pro cislovani sdof od 0 jako u ci
  // alokace pole adres pro vybrane dof -> mala matice 
  // tvorba pomocneho pole sl. indexu
  smscr->adr = new long[nsdof+1];
  auxci = new long[nsdof];
  for(i = 0; i < nsdof; i++){
    sdof[i]--;
    smscr->adr[i] = -1;
    auxci[i] = i;
  }

  // negm - pocet prvku v male matici
  smscr->negm = 0;
  
  for(i = 0; i < nsdof; i++){
    address1 = adr[sdof[i]];
    address2 = adr[sdof[i]+1];
    m = 0;
    for(j = address1; j < address2; j++){
      // vypocet poctu prvku v male matici 
      for(k = 0; k < nsdof; k++){
	if(ci[j] == sdof[k]){
	  smscr->negm++;
	  m++;
	  break;
	}
      }
    }
    // vytvoreni pole adr pro malou matici
    smscr->adr[i] = smscr->negm-m;
  }
  smscr->adr[nsdof]= smscr->negm;
  
  
  
  smscr->ci  = new long[smscr->negm];
  smscr->a   = new double[smscr->negm];
  
  m = 0;
  for(i = 0; i < nsdof; i++){
    address1 = adr[sdof[i]];
    address2 = adr[sdof[i]+1];
    for(j = address1; j < address2; j++){
      for(k = 0; k < nsdof; k++){
	if(ci[j] == sdof[k]){
	  // pridani prvku z a do male matice
	  smscr->a[m] = a[j];
	  // ci 
	  smscr->ci[m] = auxci[k];
	  m++;
	}
      }

    }
  }
  
  delete []auxci;

}
