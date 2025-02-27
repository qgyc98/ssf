#include "dskyline.h"
#include <math.h>
#include <limits.h>
#include <string.h>

dskyline::dskyline (void)
{
  //  the number of rows in the matrix (it is equal to the number of unknowns)
  n=0;
  //  the number of stored matrix entries
  negm=0;
  
  //  array containing addresses of diagonal matrix entries
  adr = NULL;
  //  array containing matrix entries
  a = NULL;
  //  indicator of matrix factorization
  decompid=0;
}



dskyline::~dskyline (void)
{
  delete [] adr;
  delete [] a;
}

long dskyline::decomp ()
{
  return decompid;
}

void dskyline::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}

void dskyline::setfact ()
{
  decompid=1;
}
void dskyline::setnotfact ()
{
  decompid=0;
}

/**
   function allocates array containing addresses of diagonal entries

   @param m - number of unknowns in solved problems
   
   JK
*/
void dskyline::allocadr (long m)
{
  n=m;
  adr = new long [n+1];
  memset (adr,0,(n+1)*sizeof(long));
}

double* dskyline::status ()
{
  return a;
}

/**
   function computes contributions to the column lengths from one element
   
   @param cn - array of code numbers on actual element
   @param ndofe - number of DOFs on actual element
   
   21.2.2002
*/
void dskyline::column_lengths_elem (long *cn,long ndofe)
{
  long i,j,k,min;
  
  min=LONG_MAX;
  for (i=0;i<ndofe;i++){
    k=cn[i]-1;
    if (k>-1){
      if (k<min)  min=k;
    }
  }
  for (i=0;i<ndofe;i++){
    k=cn[i]-1;
    if (k>-1){
      j=k-min+1;
      if (j>adr[k])  adr[k]=j;
    }
  }
}

/**
   function computes contributions to the column lengths from Lagrange multipliers
   
   @param ncn1 - code numbers of first node (on one side of interface)
   @param ncn2 - code numbers of second node (on the other side of interface)
   @param mcn - code numbers of multipliers defined in between
   @param nm - number of multipliers
   
   JK, 6.8.2008
*/
void dskyline::column_lengths_mult (long *ncn1,long *ncn2,long *mcn,long nm)
{
  long i,j,k,min;
  
  for (i=0;i<nm;i++){
    //  minimum code number
    if (ncn1[i]<ncn2[i])
      min=ncn1[i];
    else
      min=ncn2[i];

    //  distance of the farthest nonzero offdiagonal entry
    j=mcn[i]-min+1;
    
    //  actual column id
    k=mcn[i]-1;

    //  comparison with existing value
    if (j>adr[k])
      adr[k]=j;
  }
}


/**
  function prepares array of the column lengths
  lengths will be used for addresses of the diagonal entries
  column lenghts are collected in the array address

  21.6.2001
*/
void dskyline::column_lengths(gtopology *top)
{
  long i,ne,ndofe,nm,*cn,*ncn1,*ncn2,*mcn;

  //  number of elements
  ne=top->ne;

  //  contributions from elements
  for (i=0;i<ne;i++){
    ndofe=top->give_ndofe (i);
    cn = new long [ndofe];
    top->give_code_numbers (i,cn);
    column_lengths_elem (cn,ndofe);
    delete [] cn;
  }

  //  this is used in cases, where Lagrange multipliers are
  //  used for saddle point problems
  //  in many cases, the number of end nodes and general edges is zero
  //  and this part is not executed
  for (i=0;i<top->nen;i++){
    //  number of multipliers
    nm=top->endnodes[i].ndofn;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    //  code numbers on nodes and appropriate multipliers
    top->give_endnode_code_numbers (i,ncn1,ncn2,mcn);
    //  computation of column lengts cased by general edges
    column_lengths_mult (ncn1,ncn2,mcn,nm);
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }

  //  contributions from general edges
  for (i=0;i<top->nged;i++){
    //  number of multipliers between first nodes
    nm=top->gedges[i].ndofnf;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    //  code numbers on first nodes and appropriate multipliers
    top->give_edge_code_numbers (i,1,ncn1,ncn2,mcn);
    //  computation of column lengts cased by general edges
    column_lengths_mult (ncn1,ncn2,mcn,nm);
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
    
    //  number of multipliers between last nodes
    nm=top->gedges[i].ndofnl;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    //  code numbers on first nodes and appropriate multipliers
    top->give_edge_code_numbers (i,2,ncn1,ncn2,mcn);
    //  computation of column lengts cased by general edges
    column_lengths_mult (ncn1,ncn2,mcn,nm);
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }
}


/**
   function computes addresses of diagonal entries in the %matrix
   
   21.2.2002
*/
void dskyline::addresses ()
{
  long i,j,k;
  
  j=adr[0];  adr[0]=0;
  for (i=1;i<n+1;i++){
    k=adr[i];
    adr[i]=j;
    j+=k;
  }

}



/**
   function determines number of entries in the dskyline which is
   equal to the size of the array containing the global %matrix
   
   21.2.2002
*/
void dskyline::neglobmat ()
{
  negm = 2*adr[n];
}



/**
   function allocates array containing global %matrix
*/
void dskyline::allocglomat ()
{
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));
}

/**
   function fills dskyline array by zero
   
   21.2.2002
*/
void dskyline::nullsky ()
{
  long i;
  for (i=0;i<negm;i++){
    a[i]=0.0;
  }
}

/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as a dense %matrix in row ordering
   b is object of class %matrix

   @param b - array containing local %matrix
   @param cn - array containing code numbers of finite element
   
   21.2.2002
*/
void dskyline::localize (matrix &b,long *cn)
{
  long i,j,ii,jj,m;
  
  m=adr[n];
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    for (j=0;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)
	a[adr[jj]+jj-ii+m] += b[i][j];
      else
	a[adr[ii]+ii-jj]   += b[i][j];
    }
  }
}

/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as a dense %matrix in row ordering
   b is object of class %matrix

   @param b - array containing local %matrix
   @param cn - array containing code numbers of finite element
   @param k - order of local %matrix (number of rows or columns)
   
   21.2.2002
*/
void dskyline::localized (double *b,long *cn,long k)
{
  long i,j,ii,jj,m;
  
  m=adr[n];
  
  for (i=0;i<k;i++){
    ii=cn[i]-1;
    if (ii==-1)  continue;
    for (j=0;j<k;j++){
      jj=cn[j]-1;
      if (jj==-1)  continue;
      if (ii<jj)
	a[adr[jj]+jj-ii+m] += b[i*k+j];
      else
	a[adr[ii]+ii-jj]   += b[i*k+j];
    }
  }
}

/**
   function localizes general rectangular (non-square) %matrix b
   into global %matrix
   
   @param b - array containing local %matrix
   @param rcn - row code numbers
   @param ccn - column code numbers
   
   JK
*/
void dskyline::glocalize (matrix &b,long *rcn,long *ccn)
{
  long i,j,ii,jj,m;
  
  m=adr[n];

  for (i=0;i<b.m;i++){
    ii=rcn[i]-1;
    if (ii<0) continue;
    for (j=0;j<b.n;j++){
      jj=ccn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)
	a[adr[jj]+jj-ii+m] += b[i][j];
      else
	a[adr[ii]+ii-jj]   += b[i][j];
      
    }
  }
  
}

/**
   function localizes contributions from Lagrange multipliers to the %dskyline storage
   
   @param nm - number of Lagrange multipliers
   @param ncn1 - nodal code numbers of the first node (on one side of interface)
   @param ncn2 - nodal code numbers of the second node (on the other side of interface)
   @param mcn - code numbers of Lagrange multipliers defined between the previous nodes

   JK, 8.8.2008
*/
void dskyline::mult_localize (long nm,long *ncn1,long *ncn2,long *mcn)
{
  long i,ai,cn,m,mm;
  
  m=adr[n];

  //  loop over number of multipliers
  for (i=0;i<nm;i++){
    cn=mcn[i]-1;
    if (cn>-1){
      //  position in the array a
      ai=adr[cn];
      
      mm=ncn1[i]-1;
      if (mm>-1){
	a[ai+cn-mm]=-1.0;
	a[ai+cn-mm+m]=-1.0;
      }
      mm=ncn2[i]-1;
      if (mm>-1){
	a[ai+cn-mm]=1.0;
	a[ai+cn-mm+m]=1.0;
      }
    }      
  }
}


void dskyline::initiate (gtopology *top,long ndof,long mespr)
{
  if (status()==NULL){
    allocadr (ndof);
    column_lengths (top);
    addresses ();
    neglobmat ();
    allocglomat ();
  }
  else{
    nullsky ();
  }
  
  if (mespr==1)  fprintf (stdout,"\n number of dskyline entries  negm %ld",negm);
}


/**
   function computes %matrix-vector product A.x=y
   
   @param x - input vector
   @param y - output vector
   
   JK, 22.2.2002
*/
void dskyline::mxv_dsky (double *x,double *y)
{
  long i,j,ii,m;
  double s;

  m=adr[n];

  for (i=0;i<n;i++){
    ii=i-(adr[i+1]-adr[i])+1;  s=0.0;
    for (j=adr[i+1]-1;j>=adr[i];j--){
      s+=a[j]*x[ii];  ii++;
    }
    y[i]=s;
  }
  for (i=0;i<n;i++){
    ii=i-(adr[i+1]-adr[i])+1;  s=x[i];
    for (j=adr[i+1]+m-1;j>adr[i]+m;j--){
      y[ii]+=a[j]*s;  ii++;
    }
  }
}

/**
   function computes %matrix-vector product A.x=y
   
   @param x - input vector
   @param y - output vector
   
   JK, 20.9.2002
*/
void dskyline::mtxv_dsky (double *x,double *y)
{
  long i,j,ii,m;
  double s;

  m=adr[n];

  for (i=0;i<n;i++){
    ii=i-(adr[i+1]-adr[i])+1;  s=0.0;
    for (j=adr[i+1]-1+m;j>=adr[i]+m;j--){
      s+=a[j]*x[ii];  ii++;
    }
    y[i]=s;
  }
  for (i=0;i<n;i++){
    ii=i-(adr[i+1]-adr[i])+1;  s=x[i];
    for (j=adr[i+1]-1;j>=adr[i];j--){
      y[ii]+=a[j]*s;  ii++;
    }
  }
}

void dskyline::addmat_dsky (double c,dskyline &dsky)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]+=c*dsky.a[i];
  }
}

void dskyline::scalmat_dsky (double c)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]*=c;
  }
}

/**
   function fills skyline array by zero
   
   25.1.2002
*/
void dskyline::copy_dsky (dskyline &dsky)
{
  long i;
  
  //  the number of rows/columns
  dsky.n=n;
  //  the number of matrix entries stored
  dsky.negm=negm;
  
  if (dsky.adr!=NULL){
    delete [] dsky.adr;
  }
  dsky.adr = new long [n+1];
  if (dsky.a!=NULL){
    delete [] dsky.a;
  }
  dsky.a = new double [negm];
  
  for (i=0;i<=n;i++){
    dsky.adr[i]=adr[i];
  }
  for (i=0;i<negm;i++){
    dsky.a[i]=a[i];
  }
}

/**
   function decomposes general real %matrix to LU form
   
       | l_11    0    0    0   |
       | l_21 l_22    0    0   |
   L = | l_31 l_32 l_33    0   |
       |  .    .    .      .   |
       |  .    .    .      .   |
       
       |  1 u_12 u_13 u_14 ... |
       |  0    1 u_23 u_24 ... |
   U = |  0    0    1 u_34 ... |
       |  0    0    0    1 ... |
       |  .    .    .    .     |

   @param x - array containing left hand side
   @param x - array containing left hand side
   @param zero - computer zero
   @param tc - type of computation

   tc=1 - LU decomposition and back substitution
   tc=2 - only LU decomposition
   tc=3 - only back substitution
   
   JK, 21.2.2002
*/
void dskyline::lu_dsky (double *x,double *y,double zero,long tc)
{
  long i,j,k,ii,jj,kk,lj,lk,uk,m;
  double s;
  
  m = adr[n];
  
  //
  //  matrix decomposition into L.U
  //
  if (tc==1 || tc==2){
    if (decompid == 0){
      
      for (i=0;i<n;i++){
	
	lj=i-adr[i+1]+adr[i]+1;
	
	for (j=lj;j<i;j++){
	  kk=j-adr[j+1]+adr[j]+1;
	  if (kk<lj)  ii=lj;
	  else        ii=kk;
	  
	  s=0.0;  uk=adr[j]+j-ii;  lk=adr[j];  jj=adr[i]+i-ii+m;
	  for (k=uk;k>lk;k--){
	    s+=a[k]*a[jj];  jj--;
	  }
	  a[jj]=(a[jj]-s)/a[lk];
	  
	  s=0.0;  uk+=m;  lk+=m;  jj=adr[i]+i-ii;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[jj]-=s;
	}
	
	s=0.0;  uk=adr[i+1]-1;  lk=adr[i];  jj=adr[i+1]-1+m;
	for (k=uk;k>lk;k--){
	  s+=a[jj]*a[k];  jj--;
	}
	a[lk]-=s;
	if (fabs(a[lk])<zero){
	  fprintf (stderr,"\n\n zero diagonal entry in the dskyline (%ld-th row)",i);
	}
      }
      
      decompid=1;
    }
  }
  
  //
  //  back substitution
  //
  if (tc==1 || tc==3){
    if (decompid == 0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      for (i=0;i<n;i++){
	ii=adr[i];
	s=0.0;  lj=i+adr[i]-adr[i+1]+1;  k=adr[i+1]-1;
	for (j=lj;j<i;j++){
	  s+=a[k]*y[j];  k--;
	}
	if (fabs(a[ii])<zero){
	  fprintf (stderr,"\n\n zero diagonal entry in the dskyline (%ld-th row)",i);
	}
	y[i]=(y[i]-s)/a[ii];
      }
      
      for (i=n-1;i>-1;i--){
	x[i]=y[i];  s=x[i];
	lj=i-adr[i+1]+adr[i]+1;  k=adr[i+1]-1+m;
	for (j=lj;j<i;j++){
	  y[j]-=a[k]*s;  k--;
	}
      }
    }
  }
}

/**
   function eliminates internal unknowns of non-symmetric system of equations
   function is based on L.U %matrix decomposition
   
   @param b - reduced %matrix (stored as dense %matrix)
   @param c - condensed right hand side %vector
   @param x - array containing unknowns (left hand side)
   @param y - array containing right hand side
   @param zero - computer zero
   @param nr - number of remaining unknowns
   @param tc - type of computation
   
   JK, 18.4.2002
*/
void dskyline::lukon_dsky (double *b,double *c,double *x,double *y,double zero,long nr,long tc)
{
  long i,j,k,l,ii,jj,kk,lj,lk,uk,uj,m;
  double s;
  
  m = adr[n];
  
  //
  //  matrix decomposition into L.U
  //
  if (tc==1){
    if (decompid==0){
      //  internal unknown elimination
      for (i=0;i<n-nr;i++){
	
	lj=i-adr[i+1]+adr[i]+1;
	
	for (j=lj;j<i;j++){
	  kk=j-adr[j+1]+adr[j]+1;
	  if (kk<lj)  ii=lj;
	  else        ii=kk;
	  
	  s=0.0;  uk=adr[j]+j-ii;  lk=adr[j];  jj=adr[i]+i-ii+m;
	  for (k=uk;k>lk;k--){
	    s+=a[k]*a[jj];  jj--;
	  }
	  a[jj]=(a[jj]-s)/a[lk];
	  
	  s=0.0;  uk+=m;  lk+=m;  jj=adr[i]+i-ii;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[jj]-=s;
	}
	
	s=0.0;  uk=adr[i+1]-1;  lk=adr[i];  jj=adr[i+1]-1+m;
	for (k=uk;k>lk;k--){
	  s+=a[jj]*a[k];  jj--;
	}
	a[lk]-=s;
	if (fabs(a[lk])<zero){
	  fprintf (stderr,"\n\n zero diagonal entry (%e < %e) in the dskyline (%ld-th row), line %d",a[lk],zero,i,__LINE__);
	}
	
	
      }
      
      //  modification of remaining part of matrix
      for (i=n-nr;i<n;i++){
	lj=i-adr[i+1]+adr[i]+1;
	
	for (j=lj;j<n-nr;j++){
	  kk=j-adr[j+1]+adr[j]+1;
	  if (kk<lj)  ii=lj;
	  else        ii=kk;
	  s=0.0;  uk=adr[j]+j-ii;  lk=adr[j];  jj=adr[i]+i-ii+m;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[jj]=(a[jj]-s)/a[lk];
	  
	  s=0.0;  uk+=m;  lk+=m;  jj=adr[i]+i-ii;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[jj]-=s;
	}
	
	for (j=n-nr;j<i;j++){
	  kk=j-adr[j+1]+adr[j]+1;
	  if (kk<lj)  ii=lj;
	  else        ii=kk;
	  s=0.0;  uk=adr[j]+j-ii;  lk=adr[j]+j-(n-nr);  jj=adr[i]+i-ii+m;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[adr[i]+i-j+m]-=s;
	  
	  s=0.0;  uk=adr[i]+i-ii;  lk=adr[i]+i-(n-nr);  jj=adr[j]+j-ii+m;
	  for (k=uk;k>lk;k--){
	    s+=a[jj]*a[k];  jj--;
	  }
	  a[adr[i]+i-j]-=s;
	  
	}
	
	s=0.0;  uk=adr[i+1]-1;  lk=adr[i]+i-(n-nr);  jj=uk+m;
	for (k=uk;k>lk;k--){
	  s+=a[k]*a[jj];  jj--;
	}
	a[adr[i]]-=s;
	
	
      }
      
      decompid=1;
    }
    
    //  assembling of condensed matrix
    /*
    for (i=n-nr;i<n;i++){
      for (j=n-nr;j<=i;j++){
	b[(i-n+nr)*nr+j-n+nr]=a[adr[i]+i-j];
      }
      for (j=i+1;j<n;j++){
	b[(i-n+nr)*nr+j-n+nr]=a[adr[j]+j-i+m];
      }
    }
    */
    
    long ack,ack1,acrk,ri,ci;
    for (i=n-nr;i<n;i++){
      ack=adr[i];  ack1=adr[i+1]-1;
      acrk=i-ack1+ack;
      if (acrk<n-nr){
	ri=i-n+nr;  ci=i-n+nr;
	j=ack;
	b[ri*nr+ci]=a[j];
	ri--;
	for (j=ack+1;j<=i-n+nr+ack;j++){
	  b[ri*nr+ci]=a[j+m];
	  b[ci*nr+ri]=a[j];
	  //b[ri*nr+ci]=a[j];
	  //b[ci*nr+ri]=a[j+m];
	  ri--;
	}
      }else{
	ri=i-n+nr;  ci=i-n+nr;
	j=ack;
	b[ri*nr+ci]=a[j];
	ri--;
	for (j=ack+1;j<=ack1;j++){
	  b[ri*nr+ci]=a[j+m];
	  b[ci*nr+ri]=a[j];
	  //b[ri*nr+ci]=a[j];
	  //b[ci*nr+ri]=a[j+m];
	  ri--;
	}
      }
    }
    
    //  modification of right hand side vector
    for (i=0;i<n-nr;i++){
      ii=adr[i];
      s=0.0;  lj=i+adr[i]-adr[i+1]+1;  k=adr[i+1]-1;
      for (j=lj;j<i;j++){
	s+=a[k]*y[j];  k--;
      }
      if (fabs(a[ii])<zero){
	print_err("zero diagonal entry (%e < %e) in the dskyline (%ld-th row)",__FILE__, __LINE__, __func__,a[ii],zero,i);
      }
      y[i]=(y[i]-s)/a[ii];
    }
    
    for (i=n-nr;i<n;i++){
      s=0.0;  lj=i+adr[i]-adr[i+1]+1;  k=adr[i+1]-1;
      for (j=lj;j<n-nr;j++){
	s+=a[k]*y[j];  k--;
      }
      y[i]-=s;
    }
    
    //  assembling of right hand side
    l=0;
    for (k=n-nr;k<n;k++){
      c[l]=y[k];
      l++;
    }
  }// end of the statement if (tc==1)
  
  //
  //  back substitution
  //
  if (tc==2){
    if (decompid == 0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{

      //  interface components are localized to the vector of subdomain
      j=0;
      for (k=n-nr;k<n;k++){
	x[k]=c[j];  j++;
	//fprintf (stdout,"\n x %ld %ld   %e",k,j,x[k]);
      }

      for (i=n-1;i>=n-nr;i--){
	s=x[i];
	lj=adr[i]+1+m;  uj=adr[i+1]+m;  k=i-1;
	for (j=lj;j<uj;j++){
	  y[k]-=a[j]*s;  k--;
	}
      }
      
      for (i=n-nr-1;i>-1;i--){
	x[i]=y[i];  s=x[i];
	lj=i-adr[i+1]+adr[i]+1;  k=adr[i+1]-1+m;
	for (j=lj;j<i;j++){
	  y[j]-=a[k]*s;  k--;
	}
      }
    }
  }
  
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
   
   16.8.2000/20.9.2002
*/
void dskyline::bicg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv)
{
  long i,j;
  double normy,normr,normrr,nom,denom,alpha,beta;
  double *r,*rr,*d,*dd,*p,*pp;

  r = new double [n];  rr = new double [n];
  d = new double [n];  dd = new double [n];
  p = new double [n];  pp = new double [n];
  
  normy = ss (y,y,n);
  
  if (normy<zero){
    fprintf (stderr,"\n\n norm of right hand side in bi-conjugate gradient method is smaller than %e",zero);
    fprintf (stderr,"\n see file %s, line %d.\n",__FILE__,__LINE__);
    ares=normy;  ani=0;
    return;
  }

  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  
  mxv_dsky (x,p);
  for (i=0;i<n;i++){
    r[i]=y[i]-p[i];
    d[i]=r[i];
    rr[i]=r[i];
    dd[i]=rr[i];
  }
  
  nom = ss (r,rr,n);
  
  for (i=0;i<ni;i++){
    mxv_dsky (d,p);
    mtxv_dsky (dd,pp);
    
    denom = ss (dd,p,n);
    alpha = nom/denom;
    
    for (j=0;j<n;j++){
      x[j]+=alpha*d[j];
      r[j]-=alpha*p[j];
      rr[j]-=alpha*pp[j];
    }

    normr = ss (r,r,n);
    normrr = ss (rr,rr,n);
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

void dskyline::printmat (FILE *out)
{
  long i,j,m;
  /*
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"\n");
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"\n %20.10f",a[j]);
    }
  }
  */

  m=adr[n];

  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"\n");
    for (j=0;j<=i;j++){
      if (adr[i]+i-j>=adr[i+1])  fprintf (out," %e",0.0);
      else fprintf (out," %e",a[adr[i]+i-j]);
    }
    for (j=i+1;j<n;j++){
      if (adr[j]+j-i>=adr[j+1])  fprintf (out," %e",0.0);
      else fprintf (out," %e",a[adr[j]+j-i+m]);
    }
  }
}

void dskyline::printdiag (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"%5ld  %20.10f\n",i,a[adr[i]]);
  }
}


/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long dskyline::give_negm ()
{
  return negm;
}

/**
   function checks diagonal entries
   
   the function is used in some nonlinear nostationary problems
   where high jumps in coefficients occur
   some of element matrices are zero matrices and this
   function puts nonzero values on the diagonal
   
   @param thr - prescribed threshold
   @param rhs - %vector of the right hand side
   
   JK, 14.7.2008
*/
void dskyline::diag_check (double thr,double *rhs)
{
  long i,j,nch;
  
  nch=0;
  for (i=0;i<n;i++){
    j=adr[i];
    if (fabs(a[j])<thr){
      a[j]=1.0;
      rhs[i]=0.0;
      nch++;
    }
  }
  fprintf (stdout,"\n number of changes in the function dskyline::diag_check %ld\n",nch);
}
