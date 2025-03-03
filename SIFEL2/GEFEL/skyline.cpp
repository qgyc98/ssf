#include "skyline.h"
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
//#include "block2.cpp"
#include "cr.h"
#include "scr.h"
#include "gmatrix.h"

skyline::skyline (void)
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


skyline::~skyline (void)
{
  delete [] adr;
  delete [] a;
}

/**
   function returns indicator of decomposition (factorization)
   
   decompid=0 - %matrix is not decomposed (factorized)
   decompid=1 - %matrix is decomposed (factorized)
   
   JK
*/
long skyline::decomp ()
{
  return decompid;
}

/**
   function changes indicator of decomposition (factorization)
   
   JK
*/
void skyline::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}

/**
   set up indicator to factorized
*/
void skyline::setfact ()
{
  decompid=1;
}

/**
   set up indicator to not factorized
*/
void skyline::setnotfact ()
{
  decompid=0;
}


/**
   function allocates array containing addresses of diagonal entries

   @param m - number of unknowns in solved problems
   
   JK
*/
void skyline::allocadr (long m)
{
  n=m;
  adr = new long [n+1];
  memset (adr,0,(n+1)*sizeof(long));
}

/**
   function returns status of array a
   
   JK
*/
double* skyline::status ()
{
  return a;
}

/**
   function computes contributions to the column lengths from one element
   
   @param cn - array of code numbers on actual element
   @param ndofe - number of DOFs on actual element
   
   JK, 21.6.2001
*/
void skyline::column_lengths_elem (long *cn,long ndofe)
{
  long i,j,k,min;
  
  min=LONG_MAX;
  for (i=0;i<ndofe;i++){
    k=cn[i]-1;
    if (k>-1){
      if (k<min)  min=k;
    }
    if (k > n){
      print_err("invalid DOF number %ld\n", __FILE__, __LINE__, __func__, k);
    }
  }
  if (min == LONG_MAX){
    print_err("column length failure\n", __FILE__, __LINE__, __func__);
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
void skyline::column_lengths_mult (long *ncn1,long *ncn2,long *mcn,long nm)
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
  
  @param top - pointer to general topology
  
  JK, 21.6.2001
*/
void skyline::column_lengths (gtopology *top)
{
  long i,ne,ndofe,nm,*ncn1,*ncn2,*mcn;
  ivector cn;

  //  number of elements
  ne=top->ne;
  
  //  contributions from elements
  for (i=0; i<ne; i++){
    if (top->leso[i] == 1){
      ndofe = top->give_gndofe(i);
      reallocv(RSTCKIVEC(ndofe, cn));
      top->give_gcode_numbers(i, cn.a);
      column_lengths_elem(cn.a, ndofe);
    }
  }
  
  //  this is used in cases, where Lagrange multipliers are
  //  used for saddle point problems
  //  in many cases, the number of end nodes and general edges is zero
  //  and this part is not executed
  
  //  contributions from end nodes
  for (i=0;i<top->nen;i++){
    //  number of multipliers
    nm=top->endnodes[i].ndofn;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    //  code numbers on nodes and appropriate multipliers
    top->give_endnode_code_numbers (i,ncn1,ncn2,mcn);
    //  computation of column lengts caused by general edges
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
    //  computation of column lengts caused by general edges
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
    //  computation of column lengts caused by general edges
    column_lengths_mult (ncn1,ncn2,mcn,nm);
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }
}




/**
   function computes addresses of diagonal entries in the %matrix
   
   JK, 21.6.2001
*/
void skyline::addresses ()
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
   function determines number of entries in the skyline which is
   equal to the size of the array containing the global %matrix
   
   JK, 21.6.2001
*/
void skyline::neglobmat ()
{
  negm = adr[n];
}



/**
   function allocates array containing global %matrix
   
   JK
*/
void skyline::allocglomat ()
{
  delete [] a;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));
}



/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as a dense %matrix in row ordering
   b is object of class %matrix

   @param b - array containing local %matrix
   @param cn - array containing code numbers of finite element
   
   JK, 25.6.2001
*/
void skyline::localize (matrix &b,long *cn)
{
  long i,j,ii,jj,k;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    for (j=0;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii>jj)  continue;
      k=adr[jj]+jj-ii;
      a[k]+=b[i][j];
    }
  }
  
  /*
    zmena kvuli spojenym stupnum volnosti

  long i,j,ii,jj,k;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    for (j=i;j<b.m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)  k=adr[jj]+jj-ii;
      else        k=adr[ii]+ii-jj;
      a[k]+=b[i][j];
    }
  }
  */
}




/**
   function localizes local %matrix b to the global %matrix
   %matrix b is stored as a dense %matrix in row ordering
   b is double array pointer
   
   @param b - array containing local %matrix
   @param cn - array containing code numbers
   @param m - order of local %matrix (number of rows or columns)

   JK, 25.6.2001
*/
void skyline::localized (double *b,long *cn,long m)
{
  long i,j,ii,jj,k;
  
  for (i=0;i<m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    for (j=0;j<m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii>jj)  continue;
      k=adr[jj]+jj-ii;
      a[k]+=b[i*m+j];
    }
  }

  /*
    zmena kvuli spojenym stupnum volnosti
  long i,j,ii,jj,k;
  
  for (i=0;i<m;i++){
    ii=cn[i]-1;
    if (ii<0) continue;
    for (j=i;j<m;j++){
      jj=cn[j]-1;
      if (jj<0)  continue;
      if (ii<jj)  k=adr[jj]+jj-ii;
      else        k=adr[ii]+ii-jj;
      a[k]+=b[i*m+j];
    }
  }
  */
}

/**
   function localizes general rectangular (non-square) %matrix b
   into global %matrix
   
   @param b - array containing local %matrix
   @param rcn - row code numbers
   @param ccn - column code numbers
   
   JK
*/
void skyline::glocalize (matrix &b,long *rcn,long *ccn)
{
  long i,j,ii,jj,k;
  
  for (i=0;i<b.m;i++){
    ii=rcn[i]-1;
    if (ii<0) continue;
    for (j=0;j<b.n;j++){
      jj=ccn[j]-1;
      if (jj<0)  continue;
      if (ii>jj)  continue;
      k=adr[jj]+jj-ii;
      a[k]+=b[i][j];
    }
  }
  
}

/**
   function localizes contributions from Lagrange multipliers to the %skyline storage
   
   @param nm - number of Lagrange multipliers
   @param ncn1 - nodal code numbers of the first node (on one side of interface)
   @param ncn2 - nodal code numbers of the second node (on the other side of interface)
   @param mcn - code numbers of Lagrange multipliers defined between the previous nodes

   JK, 8.8.2008
*/
void skyline::mult_localize (long nm,long *ncn1,long *ncn2,long *mcn)
{
  long i,ai,cn,m;
  
  //  loop over number of multipliers
  for (i=0;i<nm;i++){
    cn=mcn[i]-1;
    if (cn>-1){
      //  position in the array a
      ai=adr[cn];
      
      m=ncn1[i]-1;
      if (m>-1)
	a[ai+cn-m]=-1.0;
      m=ncn2[i]-1;
      if (m>-1)
	a[ai+cn-m]=1.0;
    }      
  }
}

/**
   function fills skyline array by zeros
   
   JK, 25.1.2002
*/
void skyline::nullsky ()
{
  long i;
  for (i=0;i<negm;i++){
    a[i]=0.0;
  }
}

/**
   function copies actual skyline into another one
   
   @param sky - another skyline
   
   JK, 25.1.2002
*/
void skyline::copy_sky (skyline &sky)
{
  long i;
  
  //  the number of rows/columns
  sky.n=n;
  // the number of matrix entries stored
  sky.negm=negm;
  
  if (sky.adr!=NULL){
    delete [] sky.adr;
  }
  sky.adr = new long [n+1];
  if (sky.a!=NULL){
    delete [] sky.a;
  }
  sky.a = new double [negm];
  
  for (i=0;i<=n;i++){
    sky.adr[i]=adr[i];
  }
  for (i=0;i<negm;i++){
    sky.a[i]=a[i];
  }
}

/**
   function initiates skyline storage
   
   @param top - pointer to general topology
   @param ndof - number of rows/columns of the %matrix
   @param mespr - indicator of message printing
   
   JK
*/
void skyline::initiate (gtopology *top,long ndof,long mespr)
{
  if (status()==NULL){
    allocadr (ndof);
    column_lengths (top);
    addresses ();
    neglobmat ();
    if (mespr==1){
      fprintf (stdout,"\n number of skyline entries  negm %ld",negm);
      fflush (stdout);
    }
    allocglomat ();
    return;
  }
  else{
    nullsky ();
  }
  if (mespr==1){
    fprintf (stdout,"\n number of skyline entries  negm %ld",negm);
    fflush (stdout);
  }
}





/**
   function computes column lengths
   the %matrix is stored in the format compressed rows
   this strategy is developed for the BOSS method
   
   @param cradr - array of the first entries in rows
   @param ci - array of column indices
   
   JK, 8.5.2007
*/
void skyline::column_lengths_cr (long *cradr,long *ci)
{
  long i,j;
  
  for (i=0;i<n;i++){
    j=cradr[i];
    adr[i]=i-ci[j]+1;
  }
}

/**
   function assembles array of %matrix entries
   
   @param cradr - addresses in compressed row format
   @param ci - column indices in compressed row format
   @param cra - %matrix entries stored in the compressed row format
   
   JK, 8.5.2007
*/
void skyline::mat_entries (long *cradr,long *ci,double *cra)
{
  long i,j,lj,uj,k;
  
  for (i=0;i<n;i++){
    lj=cradr[i];
    uj=cradr[i+1]-1;
    for (j=uj;j>=lj;j--){
      if (ci[j]<=i){
	k=adr[i]+i-ci[j];
	a[k]=cra[j];
      }
    }
  }
}

/**
   function assembles all data about %matrix stored in the %skyline
   which is defined in the compressed row format
   
   @param cr - pointer to the compressed row format
   
   JK, 8.5.2007
*/
void skyline::assemble_from_cr (long crn,long *cradr,long *crci,double *cra)
{
  //  allocation of the array adr and definition of the number of rows/columns
  allocadr (crn);
  
  //  lengths of columns in the skyline
  column_lengths_cr (cradr,crci);
  
  //  addresses of diagonal entries
  addresses ();
  
  //  number of matrix entries stored
  neglobmat ();
  
  //  allocation of the array for matrix entries
  allocglomat ();
  
  //  matrix entries assembling
  mat_entries (cradr,crci,cra);
  
  fprintf (stdout,"\n skyline for aggregate contains %8ld unknowns and stores %10ld entries",n,negm);
}




/**
   function computes lengths of columns
   
   @param scr - pointer to %matrix stored in symmetric compressed rows
   @param se - array of selected unknown numbers

   JK, 22.8.2007
*/
void skyline::column_lengths_scr (long *scradr,long *scrci,long *se)
{
  long i,j,k,lk,uk;
  long *aux;
  long min,max;
  
  max=0;
  for (i=0;i<n;i++){
    if (se[i]>max)
      max=se[i];
  }
  aux = new long [max];
  for (i=0;i<max;i++){
    aux[i]=-1;
  }
  for (i=0;i<n;i++){
    aux[se[i]-1]=i;
  }
  
  
  for (i=0;i<n;i++){
    //  real number of the unknown
    j=se[i]-1;
    //  address of the first nonzero entry in the row
    lk=scradr[j];
    uk=scradr[j+1];
    min=max;
    for (k=lk;k<uk;k++){
      if (aux[scrci[k]]!=-1){
	if (aux[scrci[k]]<min)
	  min=aux[scrci[k]];
      }
    }
    adr[i]=aux[j]-min+1;
  }
  
  delete [] aux;


  
  /*
  //  auxiliary array
  //  it contains gaps between two successive numbers of selected unknowns
  //  this number has to be subtracted from the column lengths
  aux = new long [n];
  aux[0]=0;
  for (i=0;i<n-1;i++){
    if (se[i+1]-se[i] == 0)
      aux[i+1]=aux[i];
    else
      aux[i+1]=aux[i]+se[i+1]-se[i]-1;
  }
  
  for (i=0;i<n;i++){
    j=se[i]-1;
    k=scradr[j];
    //if (adr[j-aux[i]] < j - scrci[k] + 1 - aux[i])
    adr[i]=j - scrci[k] + 1 - aux[i];
  }
  
  delete [] aux;
  */
}

/**
   
JK, 22.8.2007
*/
void skyline::mat_entries_scr (long *scradr,long *scrci,double *scra,long *se)
{
  long i,j,k,l,lj,uj,ii,jj;
  long *aux;
  
  //  auxiliary array
  //  it contains gaps between two successive numbers of selected unknowns
  //  this number has to be subtracted from the column lengths
  aux = new long [n];
  aux[0]=0;
  for (i=0;i<n-1;i++){
    aux[i+1]=aux[i]+se[i+1]-se[i]-1;
  }
   
  for (i=0;i<n;i++){
    k=se[i]-1;
    lj=scradr[k];
    uj=scradr[k+1]-1;
    for (j=uj;j>lj-1;j--){
      for (ii=0;ii<n;ii++){
	if (scrci[j]==se[ii]-1){
	  jj=scrci[j]-aux[ii];
	  break;
	}
      }
      //l=adr[i] + k - aux[i] - scrci[j];
      l=adr[i] + k - aux[i] - jj;
      a[l]=scra[j];
    }
  }
  
  delete [] aux;
}


/**
   function assembles %matrix of subdomain/aggregate etc. from
   %matrix of the whole problem stored in the symmetric compressed row storage
   
   @param scr - pointer to %matrix stored in symmetric compressed rows
   @param neq - number of selected equations (rows of the %matrix)
   @param se - array of selected unknown numbers
   
   JK, 22.8.2007
*/
void skyline::assemble_from_scr (long *scradr,long *scrci,double *scra,long neq,long *se)
{
  //  allocation of the array adr and definition of the number of rows/columns
  allocadr (neq);
  
  //  computation of column lengths
  column_lengths_scr (scradr,scrci,se);

  //  addresses of diagonal entries
  addresses ();

  //  number of matrix entries stored
  neglobmat ();
  
  //  allocation of the array for matrix entries
  allocglomat ();
  
  //  matrix entries assembling
  mat_entries_scr (scradr,scrci,scra,se);

  fprintf (stdout,"\n skyline for subdomain/aggregate contains %8ld unknowns and stores %10ld entries",n,negm);
}


/**
   either function decomposes %matrix into L.D.L form
   or computes solution of linear algebraic system
   nonsingular %matrix of system is supposed
   
   @param x - array containing system solution
   @param y - array containing right hand side
   @param zero - computer zero (for testing small numbers near zero)
   @param tc - computation indicator
   
   tc=1 - L.D.L decomposition and system solution
   tc=2 - only L.D.L. decomposition
   tc=3 - only back-substitution
   
   JK, 21.6.2001
*/
void skyline::ldl_sky (double *x,double *y,double zero,long tc)
{
  long i,j,k,ac,ac1,ac2,acs,ack,ack1,acrk,aci,aci1,acri,acj,acj1;
  double s,g;
  
  
  if (tc==1 || tc==2){
    if (decompid == 0){
      //
      //  matrix decomposition into L.D.L
      //
      for (k=1;k<n;k++){
	ack=adr[k];  ack1=adr[k+1];
	ac1=k+ack;   acrk=ac1-ack1+1;
	acj=ack1-2;
	for (i=acrk+1;i<k;i++){
	  aci=adr[i];  aci1=adr[i+1];
	  ac2=i+aci;   acri=ac2-aci1+1;
	  if (acri<acrk)  ac=acrk;
	  else            ac=acri;
	  acj1=ac1-ac;  acs=ac2-ac;
	  s=0.0;
	  for (j=acj1;j>acj;j--){
	    s+=a[j]*a[acs];
	    acs--;
	  }
	  a[acj]-=s;  acj--;
	}
	s=0.0;
	for (i=ack1-1;i>ack;i--){
	  ac1=adr[acrk];  acrk++;
	  g=a[i];
	  a[i]/=a[ac1];
	  s+=a[i]*g;
	}
	a[ack]-=s;
	if (fabs(a[ack])<zero){
	  printf ("\n\n zero diagonal entry in the skyline (%ld-th row)",k);
	}
      }
      
      decompid=1;
    }
  }

  if (tc==1 || tc==3){
    if (decompid == 0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      //
      //  back-substitution
      //
      //  computation of L.z = y  => z
      for (k=1;k<n;k++){
	ack=adr[k]+1;  ack1=adr[k+1];
	ac1=k-1;  s=0.0;
	for (i=ack;i<ack1;i++){
	  s+=a[i]*y[ac1];
	  ac1--;
	}
	y[k]-=s;
      }
      //  D.L.x = z  =>  L.x = D^{-1}.z
      for (k=0;k<n;k++){
	ac=adr[k];
	if (fabs(a[ac])<zero){
	  printf ("\n\n zero diagonal entry in the skyline (%ld-th row) (D^-1)",k);
	}
	y[k]/=a[ac];
      }
      //  L.x = D^{-1}.z => x
      for (k=n-1;k>-1;k--){
	ack=adr[k]+1;  ack1=adr[k+1];
	x[k]=y[k];  g=x[k];
	ac=k-1;
	for (i=ack;i<ack1;i++){
	  y[ac]-=a[i]*g;
	  ac--;
	}
      }
    }
  }
}












/**
   either function decomposes %matrix into L.D.L form
   or computes solution of linear algebraic system
   regular %matrix of system is supposed
   
   @param x - array containing system solution
   @param y - array containing right hand side
   @param zero - computer zero (for testing small numbers near zero)
   @param tc - computation indicator
   
   tc=1 - L.D.L decomposition and system solution
   tc=2 - only L.D.L. decomposition
   tc=3 - only back-substitution
   
   JK, IS, 21.6.2001
*/
/*
void skyline::ldl_sky_new (double *x,double *y,double zero,long tc)
{
  long i,j,k,ac,ac1,ac2,acs,ack,ack1,acrk,aci,aci1,acri,acj,acj1;
  double s,g;
  
  if (tc==1 || tc==2){
    if (decompid == 1){
      fprintf (stderr,"\n\n factorized matrix is once again factorized (file %s, line %d)\n",__FILE__,__LINE__);
    }
    //
    //  matrix decomposition into L.D.L
    //
    for (k=1;k<n;k++){
      ack=adr[k];  ack1=adr[k+1];
      ac1=k+ack;   acrk=ac1-ack1+1;
      acj=ack1-2;
      for (i=acrk+1;i<k;i++){
	aci=adr[i];  aci1=adr[i+1];
	ac2=i+aci;   acri=ac2-aci1+1;
	if (acri<acrk)  ac=acrk;
	else            ac=acri;
	acj1=ac1-ac;  acs=ac2-ac;
	s=0.0;
	for (j=acj1;j>acj;j--){
	  s+=a[j]*a[acs];
	  acs--;
	}
	a[acj]-=s;  acj--;
      }
      s=0.0;
      for (i=ack1-1;i>ack;i--){
	ac1=adr[acrk];  acrk++;
	g=a[i];
	a[i]/=a[ac1];
	s+=a[i]*g;
      }
      a[ack]-=s;
      if (fabs(a[ack])<zero){
	printf ("\n\n zero diagonal entry in the skyline (%ld-th row)",k);
      }
    }
    
    decompid=1;
  }
  if (tc==1 || tc==3){
    //
    //  back-substitution
    //
    //  computation of L.z = y  => z
    for (k=1;k<n;k++){
      ack=adr[k]+1;  ack1=adr[k+1];
      ac1=k-1;  s=0.0;
      for (i=ack;i<ack1;i++){
	s+=a[i]*y[ac1];
	ac1--;
      }
      y[k]-=s;
    }
    //  D.L.x = z  =>  L.x = D^{-1}.z
    for (k=0;k<n;k++){
      ac=adr[k];
      if (fabs(a[ac])<zero){
	printf ("\n\n zero diagonal entry in the skyline (%ld-th row) (D^-1)",k);
      }
      y[k]/=a[ac];
    }
    //  L.x = D^{-1}.z => x
    for (k=n-1;k>-1;k--){
      ack=adr[k]+1;  ack1=adr[k+1];
      x[k]=y[k];  g=x[k];
      ac=k-1;
      for (i=ack;i<ack1;i++){
	y[ac]-=a[i]*g;
	ac--;
      }
    }
  }
}
*/

long skyline::auxmax(long i,long j)
{
  if (i>j) return i;
  else     return j;
}


long skyline::auxmin(long i,long j)
{
  if (j>i) return i;
  else     return j;
}

/**
   nejnovejsi verze
*/
/*
void skyline::ldl_sky_10 (double *x,double *y,double zero,long tc)
{
  long band;
  
  if (tc==1 || tc==2){
    //  bandwidth of matrix
    band=bandwidth(adr,n);
    
    //  matrix factorization
    blokove2 (a,adr,n,band);
  }
  
  long i,k,ac,ac1,ack,ack1;
  double s,g;
  if (tc==1 || tc==3){
// ********************
// *  vypocet reseni  *
    // ********************
//  vypocet  Lz=y => z (y se prepisuji na z) 
    for (k=0;k<n;k++){
      //  smycka pres nezname 
      ack=adr[k]+1;  ack1=adr[k+1];
      ac1=k-1;  s=0.0;
      for (i=ack;i<ack1;i++){
    s+=a[i]*y[ac1];
    ac1--;
      }
      y[k]=(y[k]-s)/a[ack-1];
    }
    for (k=0;k<n;k++){
      x[k]=y[k];
    }
    //  vypocet Lx=1/Dz => x  
    for (k=n-1;k>-1;k--){
      //  smycka pres nezname 
      ack=adr[k]+1;  ack1=adr[k+1];
      x[k]/=a[ack-1];
      g=x[k];
      ac=k-1;
      for (i=ack;i<ack1;i++){
        x[ac]-=a[i]*g;
        ac--;
      }
    }
  }

}
*/

void skyline::ldl_sky3 (double *x,double *y,long tc)
/*
  funkce resi soustavu linearnich algebraickych rovnic
  reseni se provadi rozkladem LDL
  matice soustavy je ulozena ve skylinu
  
  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  adr - pole adres diagonalnich prvku
  n - pocet neznamych
  tc - typ vypoctu  tc=1 - provede se eliminace i zpetny chod
                    tc=2 - provede se pouze eliminace
                    tc=3 - provede se pouze zpetny chod
                                                        
  IS, 10.6.2003


  funkce je totozna s procedurou ve fortranu, ktera je stejne
  rychla jako colsol od Batheho
  funkce byla testovana s vysledky od Batheho
                                      ldl ve fortranu
                                      stare ldl v c
*/
{

  double sa,sb,s,g;	

  long i,j,k,ac1,ack,ack1,acrk,acj;

  long aci,aci1,ac2,acri,ac,acj1,acs,pom;  

  long acia,aci1a,ac2a,acria,aca,acj1a,acsa,poma;

  long acib,aci1b,ac2b,acrib,acb,acj1b,acsb,pomb,pombz;  

  long pocet;



  

  if (tc==1 || tc==2){

    /*****************************/

    /*  rozklad matice soustavy  */

    /*****************************/

    for (k=1;k<n;k++){

      /*  smycka pres vsechny radky matice  */

      if ((k%200)==0) 
      {
        //printf("%li ",k);

	//fflush(stdout);
      }
      ack=adr[k];

      ack1=adr[k+1];

      ac1=k+ack;  //jakoby konecny

      acrk=ac1-ack1+1;  //kolik do konce

      acj=ack1-2;  // ?

      /*  uprava mimodiagonalnich prvku k-teho sloupce  */

      for (i=acrk+1;i<(k-1);i+=2){

        /*  smycka pres prvky k-teho sloupce  */

        acia=adr[i];

        aci1a=acib=adr[i+1];

        aci1b=adr[i+2];

        ac2a=i+acia; //jakoby konecny

        ac2b=i+1+acib; //jakoby konecny        

        acria=ac2a-aci1a+1; //kolik do konce

        acrib=ac2b-aci1b+1; //kolik do konce        

        aca=auxmax(acrk,acria);

        acb=auxmax(acrk,acrib);        

        acj1a=ac1-aca;

        acj1b=ac1-acb;        

        acsa=ac2a-aca;

        acsb=ac2b-acb;        

        /* old

        for (j=acj1;j>acj;j--){

          s+=a[j]*a[acs];

          acs--;

        } */



        // new

        poma=acsa+acj+1-acj1a; // (ac2-ac)+acj+1-(ac1-ac) = ac2-ac1+acj+1 = (i+aci)-(k+ack)+(ack1-2)+1 = 

        pomb=pombz=acsb+acj+1-acj1b; // (ac2-ac)+acj+1-(ac1-ac) = ac2-ac1+acj+1 = (i+aci)-(k+ack)+(ack1-2)+1 =         

        pocet=acj1a-acj; // (ac1-ac)-

        #ifdef LADIC        

        if (pocet<0) pocet=0;

        zvys_op(2l*pocet+1);                            

        pocty[pocet]++;

        #endif         

        

        /*        

        for (j=acj+1;j<=acj1;j++){

          s+=a[j]*a[pom++];

        }*/

        

        //s=cblas_ddot(pocet, a+acj+1l, 1, a+pom, 1); 

        //a[acj]-=s;  acj--; 

        sa=0.0;

        sb=0.0;        



        if (acj1a<=acj1b) 

        {

          for (j=acj+1;j<=acj1a;j++){

            g=a[j]; 	

            sa+=g*a[poma++];

            sb+=g*a[pomb++];          

          }

          for (j=acj1a+1;j<=acj1b;j++){	

            sb+=a[j]*a[pomb++];          

          }          

        }  

        else

        {

          for (j=acj+1;j<=acj1b;j++){

            g=a[j]; 	

            sa+=g*a[poma++];

            sb+=g*a[pomb++];          

          }

          for (j=acj1b+1;j<=acj1a;j++){            	

            sa+=a[j]*a[poma++];             

          }          

        }  

        

        a[acj]-=sa;

        sb+=a[acj]*a[pombz-1];

        a[acj-1]-=sb;

        acj-=2; 

        

      }//end of i

      if (i<k)

      {

        aci=adr[i];

        aci1=adr[i+1];

        ac2=i+aci; //jakoby konecny

        acri=ac2-aci1+1; //kolik do konce

        ac=auxmax(acrk,acri);

        

        acj1=ac1-ac;

        acs=ac2-ac;

        s=0.0;

        pom=acs+acj+1-acj1; // (ac2-ac)+acj+1-(ac1-ac) = ac2-ac1+acj+1 = (i+aci)-(k+ack)+(ack1-2)+1 = 

        pocet=acj1-acj; // (ac1-ac)-

        #ifdef LADIC        

        if (pocet<0) pocet=0;

        zvys_op(2l*pocet+1);                            

        pocty[pocet]++;

        #endif         

        

                

        for (j=acj+1;j<=acj1;j++){

          s+=a[j]*a[pom++];

        }

        

        //s=cblas_ddot(pocet, a+acj+1l, 1, a+pom, 1); 

        a[acj]-=s;  acj--; 

      	

      } 	

      /*  uprava diagonalniho prvku v k-tem sloupci  */

      s=0.0;

      for (i=ack1-1;i>ack;i--){

        /*  smycka pres mimodiagonalni prvky k-teho sloupce  */

        ac1=adr[acrk];  acrk++;

        g=a[i];

        a[i]/=a[ac1];

        s+=a[i]*g;

      }//end of i

      

      #ifdef LADIC      

      pocet=ack1-1-ack;

      if (pocet<0) pocet=0;      

      zvys_op(3l*pocet+1l);

      #endif

      a[ack]-=s;

    }

  }

  if (tc==1 || tc==3){

    /********************/

    /*  vypocet reseni  */

    /********************/

    /*  vypocet  Lz=y => z (y se prepisuji na z) */

    for (k=1;k<n;k++){

      /*  smycka pres nezname  */

      ack=adr[k]+1;  ack1=adr[k+1];

      ac1=k-1;  s=0.0;

      for (i=ack;i<ack1;i++){

        s+=a[i]*y[ac1];

        ac1--;

      }

      #ifdef LADIC            

      pocet=ack1-ack;

      if (pocet<0) pocet=0;      

      zvys_op(2l*pocet+1l);

      #endif      

      y[k]-=s;

    }

    /*  deleni zbyle soustavy diagonalnimi prvky (DLx=z => Lx=1/Dz) */

    for (k=0;k<n;k++){

      ac=adr[k];  y[k]/=a[ac];

    }

    #ifdef LADIC                

    zvys_op(n);

    #endif    

    /*  vypocet Lx=1/Dz => x  */

    for (k=n-1;k>-1;k--){

      /*  smycka pres nezname  */

      ack=adr[k]+1;  ack1=adr[k+1];

      x[k]=y[k];  g=x[k];

      ac=k-1;

      for (i=ack;i<ack1;i++){

        y[ac]-=a[i]*g;

        ac--;

      }

      #ifdef LADIC                  

      pocet=ack1-ack;

      if (pocet<0) pocet=0;      

      zvys_op(2l*pocet);      

      #endif      

    }

  }

}

void skyline::eliminuj_4i_rev(double *a,long n,long s)
{
  double sum1,sum2,sum3,sum4,x;
  long i,j,k;
  
  for(j=0;j<n;j++)
  {
    sum1=sum2=0.0;      
    for(k=0;k<(j-1);k+=2)
    {
      sum1+=(a[j*s+k]*a[j*s+k]);
      sum2+=(a[j*s+k+1]*a[j*s+k+1]);
    }       
    if (j&1) sum1+=(a[j*s+j-1]*a[j*s+j-1]);
    //if ((a[j*s+j]-sum1-sum2)<=0) { printf ("Chyba \n"); }
    a[j*s+j]=sqrt(a[j*s+j]-sum1-sum2);
    
    if(j&1)
    { 
      for(i=j+1;i<(n-3);i+=4)
      {
        sum1=sum2=sum3=sum4=0.0;        
        for(k=0;k<j;k++)
        {
          x=a[j*s+k];     
          sum1+=a[i*s+k]*x;
          sum2+=a[(i+1)*s+k]*x;   
          sum3+=a[(i+2)*s+k]*x;   
          sum4+=a[(i+3)*s+k]*x;   
        } 
        x=a[j*s+j];
        a[i*s+j]=(a[i*s+j]-sum1)/x;
        a[(i+1)*s+j]=(a[(i+1)*s+j]-sum2)/x;     
        a[(i+2)*s+j]=(a[(i+2)*s+j]-sum3)/x;     
        a[(i+3)*s+j]=(a[(i+3)*s+j]-sum4)/x;     
      }
      for(;i<n;i++)
      {
        sum1=0.0;       
        for(k=0;k<j;k++)
        {
          sum1+=a[i*s+k]*a[j*s+k];    
        } 
        a[i*s+j]=(a[i*s+j]-sum1)/a[j*s+j];            
      }   
    }
    else
    {
      for(i=n-1;i>(j+3);i-=4)
      {
        sum1=sum2=sum3=sum4=0.0;        
        for(k=0;k<j;k++)
        {
          x=a[j*s+k];     
          sum1+=a[i*s+k]*x;
          sum2+=a[(i-1)*s+k]*x;   
          sum3+=a[(i-2)*s+k]*x;   
          sum4+=a[(i-3)*s+k]*x;   
        } 
        x=a[j*s+j];
        a[i*s+j]=(a[i*s+j]-sum1)/x;
        a[(i-1)*s+j]=(a[(i-1)*s+j]-sum2)/x;     
        a[(i-2)*s+j]=(a[(i-2)*s+j]-sum3)/x;     
        a[(i-3)*s+j]=(a[(i-3)*s+j]-sum4)/x;     
      }
      for(;i>j;i--)
      {
        sum1=0.0;       
        for(k=0;k<j;k++)
        {
          sum1+=a[i*s+k]*a[j*s+k];    
        } 
        a[i*s+j]=(a[i*s+j]-sum1)/a[j*s+j];            
      }   
    } 
  }       
}   
  

void skyline::napln_a(long i1,long i2,long band,double *pole,double *a,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;    
  for(i=i1;i<i2;i++)    
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=min2(a2-a1,i+1-i1);
    x=pole+a1;
    for(j=0;j<kolik;j++) {a[(i-i1)*band+(i-i1)-j]=*x;   x++; } 
    for(j=kolik;j<=(i-i1);j++) a[(i-i1)*band+(i-i1)-j]=0.0;
  } 

}   


void skyline::uloz_a(long i1,long i2,long band,double *pole,double *a,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;    
  for(i=i1;i<i2;i++)    
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=min2(a2-a1,i+1-i1);
    x=pole+adr[i];
    for(j=0;j<kolik;j++) { *x=a[(i-i1)*band+(i-i1)-j];    x++; }
  }     

}   


void skyline::napln_b(long i1,long i2,long band,double *pole,double *b,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;    
  for(i=i1;i<i2;i++)    
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=i+1-i1;    
    x=pole+adr[i]+kolik;
    for(j=0;j<((a2-a1)-kolik);j++) { b[(i-i1)*band+(band-1)-j]=*x;    x++; }
    for(j=((a2-a1)-kolik);j<band;j++) b[(i-i1)*band+(band-1)-j]=0.0;
  }     

}


void skyline::uloz_b(long i1,long i2,long band,double *pole,double *b,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;    
  for(i=i1;i<i2;i++)    
  {
    a1=adr[i];
    a2=adr[i+1];
    //kolik=min2(a2-a1,i+1-i1)+1;    
    kolik=i+1-i1;    
    x=pole+adr[i]+kolik;
    for(j=0;j<((a2-a1)-kolik);j++) { *x=b[(i-i1)*band+(band-1)-j];    x++; }
  } 

}

void skyline::faze1_block3(double *a,double *b,long n,long s,long n2,long s2)
{
   //lepsi cache
   //unroll 4i
   long i,j,k,i2,j2,stepi,stepj;   
   double x,sum1,sum2,sum3,sum4;  
   
   stepi=1000;
   stepj=min2(120000l/(16l*max2(n,n2)),64l);
   for(i2=0;i2<n2;i2+=stepi)
   {
     for(j2=0;j2<n;j2+=stepj)
     {    
       for(i=i2;i<(min2(i2+stepi,n2)-3);i+=4)
       {
         for(j=j2;j<min2(j2+stepj,n);j++)
         {
           sum1=sum2=sum3=sum4=0.0;
           for(k=0;k<j;k++)
           {
             x=a[j*s+k];
             sum1+=b[i*s2+k]*x;//L                   
             sum2+=b[(i+1)*s2+k]*x;//L                   
             sum3+=b[(i+2)*s2+k]*x;//L                   
             sum4+=b[(i+3)*s2+k]*x;//L                   
           }  
           x=a[j*s+j];
           b[i*s2+j]=(b[i*s2+j]-sum1)/x;//A,L   
           b[(i+1)*s2+j]=(b[(i+1)*s2+j]-sum2)/x;//A,L   
           b[(i+2)*s2+j]=(b[(i+2)*s2+j]-sum3)/x;//A,L   
           b[(i+3)*s2+j]=(b[(i+3)*s2+j]-sum4)/x;//A,L   
         }   
       }
       for(;i<min2(i2+stepi,n2);i++)
       {
         for(j=j2;j<min2(j2+stepj,n);j++)
         {
           sum1=0.0;
           for(k=0;k<j;k++) sum1+=b[i*s2+k]*a[j*s+k];//L                   
           b[i*s2+j]=(b[i*s2+j]-sum1)/a[j*s+j];//A,L   
         }   
       }         
     }
  } 
}

void skyline::faze2_block3(double *a,double *b,long n,long s,long n2,long s2)
{
  double x,sum1,sum2,sum3,sum4;    
  long i,j,k,i2,j2,step;   
  
  //step=32l;
  //step=min2(120000l/(16l*max2(n,n2)),64l);
  step=120000l/(16l*max2(n,n2));
  for(i2=0;i2<n2;i2+=step)
  {
    for(j2=0;j2<=i2;j2+=step)
    {   
      for(i=i2;i<min2(n2,i2+step);i++)
      {
        for(j=j2;j<=(min2(i,j2+step-1)-3);j+=4)
        {
          sum1=sum2=sum3=sum4=0.0;
          for(k=0;k<n;k++)
          {
            x=b[i*s2+k];
            sum1+=x*b[j*s2+k];//A
            sum2+=x*b[(j+1)*s2+k];//A
            sum3+=x*b[(j+2)*s2+k];//A
            sum4+=x*b[(j+3)*s2+k];//A
          }  
          a[i*s+j]-=sum1;//B
          a[i*s+j+1]-=sum2;
          a[i*s+j+2]-=sum3;
          a[i*s+j+3]-=sum4;          
        }
        for(;j<=min2(i,j2+step-1);j++)
        {
            sum1=0.0;
            for(k=0;k<n;k++)  sum1+=b[i*s2+k]*b[j*s2+k];//A
            a[i*s+j]-=sum1;//B
        }       
      }
    }
  }
}

void skyline::blokove2(double *pole,long *adr,long n, long band)
{
  //zmenena faze1

  long i,i1,i2;
  double *a,*b; 
  
  //printf("Band = %li \n",band);
  a=new double [band*band];
  b=new double [band*band];
  
  i=0;  
  //if ((n-i)<(2*band)) krok=n-i;
  i1=i2=band;
  napln_a(i,i1,band,pole,a,adr);
  /*
  printf("*************************A**********************\n");  
  for(k=0;k<band;k++)
  {
    for(j=0;j<band;j++)  printf("  %g  ",a[j+k*band]);
    printf("\n");
  } */   

  /*
  printf("Pred elim \n");
  fflush(stdout);*/  
  eliminuj_4i_rev(a,band,band);
  /*
  printf("*************************sqrt(A)***********************\n");  
  for(k=0;k<band;k++)
  {
    for(j=0;j<band;j++)  printf("  %g  ",a[j+k*band]);
    printf("\n");
  } */       

  while(i1<n)
  {
    //printf(" %li \n",i1);
    i2=min2(i1+band,n);        
    //printf("Cyklus   i= %li i1= %li i2= %li\n",i,i1,i2);
    
    //napln_b(i1,i2,band,pole,b,adr);
    napln_b(i1,i2,i1-i,pole,b,adr);        
    /*
    printf("**********************************\n");
    printf("Po napln B \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",b[j+k*(i2-i1)]);
      printf("\n");
    } */       
                
    faze1_block3(a,b,i1-i,i1-i,i2-i1,i1-i);
    //faze1_block(a,b,i1-i,i2-i1,4);
    //puv faze1_block(a,b,i1-i,i1-i,i2-i1,i1-i);
    /*
    printf("**********************************\n");
    printf("Po faze1 \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",b[j+k*(i2-i1)]);
      printf("\n");
    } */       
        
    //uloz_a(i,i1,band,pole,a,adr);
    uloz_a(i,i1,i1-i,pole,a,adr);
    
    //napln_a(i1,i2,band,pole,a,adr);
    napln_a(i1,i2,i2-i1,pole,a,adr);
    /*
    printf("**********************************\n");
    printf("Po napln A \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */       
    
    faze2_block3(a,b,i1-i,i2-i1,i2-i1,i1-i);
    //faze2(a,b,i1-i,i2-i1,4);
    //puv faze2_block(a,b,i1-i,i2-i1,i2-i1,i1-i);
    /*
    printf("**********************************\n");
    printf("Po faze2 \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */           
    //uloz_b(i1,i2,band,pole,b,adr);
    uloz_b(i1,i2,i1-i,pole,b,adr);
    
    eliminuj_4i_rev(a,i2-i1,i2-i1);
    /*
    printf("**********************************\n");
    printf("Po elim A \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */       
    
    i=i1;
    i1=i2;
  }
  //uloz_a(i,n,band,pole,a,adr);
  uloz_a(i,n,i1-i,pole,a,adr);
  //tisk_sky(adr,pole,n);
}

/**
  funkce resi soustavu linearnich algebraickych rovnic
  reseni se provadi rozkladem LDL
  matice soustavy je ulozena ve skylinu
  
  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  adr - pole adres diagonalnich prvku
  n - pocet neznamych
  tc - typ vypoctu  tc=1 - provede se eliminace i zpetny chod
                    tc=2 - provede se pouze eliminace
		    tc=3 - provede se pouze zpetny chod
							
  10.7.1996
  funkce je totozna s procedurou ve fortranu, ktera je stejne
  rychla jako colsol od Batheho
  funkce byla testovana s vysledky od Batheho
	                              ldl ve fortranu
				      stare ldl v c
*/
/*
void skyline::ldl_sky4 (double *x,double *y,long tc)
{
  long i,k,ac,ac1,ack,ack1;
  double s,g;
  
  if (tc==1 || tc==2)
  {
    i=bandwidth(adr,n);
    blokove(a,adr,n,i);
  }
  if (tc==1 || tc==3){
    // ********************
    // *  vypocet reseni  *
    // ********************
    //  vypocet  Lz=y => z (y se prepisuji na z) 
    for (k=1;k<n;k++){
      //  smycka pres nezname  
      ack=adr[k]+1;  ack1=adr[k+1];
      ac1=k-1;  s=0.0;
      for (i=ack;i<ack1;i++){
	s+=a[i]*y[ac1];
	ac1--;
      }
      y[k]-=s;
    }
    //  deleni zbyle soustavy diagonalnimi prvky (DLx=z => Lx=1/Dz) 
    for (k=0;k<n;k++){
      ac=adr[k];  y[k]/=a[ac];
    }
    //  vypocet Lx=1/Dz => x  
    for (k=n-1;k>-1;k--){
      //  smycka pres nezname 
      ack=adr[k]+1;  ack1=adr[k+1];
      x[k]=y[k];  g=x[k];
      ac=k-1;
      for (i=ack;i<ack1;i++){
	y[ac]-=a[i]*g;
	ac--;
      }
    }
  }
}
*/




















/**
   function multiplies %matrix in skyline storage by %vector b
   
   @param b - array containing %vector b
   @param c - array containing resulting %vector A.b
   
   JK, 21.10.2001
*/
void skyline::mxv_sky (double *b,double *c)
{
  long i,j,acb,aci,aci1;
  double s,g;
  
  for (i=0;i<n;i++){
    aci=adr[i];  aci1=adr[i+1]-1;
    g=b[i];  s=0.0;  acb=i-aci1+aci;
    for (j=aci1;j>aci;j--){
      s+=a[j]*b[acb];
      c[acb]+=a[j]*g;
      acb++;
    }
    c[i]=s+a[aci]*g;
  }
}

/**
   function computes L^T b = c
   
   @param b - input %vector
   @param c - output %vector
   
       | 1 x x x x x x |
       | 0 1 x x x x x |
       | 0 0 1 x x x x |
   A = | 0 0 0 1 x x x |
       | 0 0 0 0 1 x x |
       | 0 0 0 0 0 1 x |
       | 0 0 0 0 0 0 1 |
       
   JK
*/
void skyline::utv (double *b,double *c)
{
  long i,j,k,lj,uj;
  double s;
  
  nullv (c,n);
  
  for (i=n-1;i>-1;i--){
    s=b[i];  k=i-1;
    c[i]+=b[i];
    lj=adr[i]+1;  uj=adr[i+1];
    for (j=lj;j<uj;j++){
      c[k]+=a[j]*s;  k--;
    }
  }
}

/**
   function computes L b = c
   
   @param b - input %vector
   @param c - output %vector
   
       | 1 0 0 0 0 0 0 |
       | x 1 0 0 0 0 0 |
   A = | x x 1 0 0 0 0 |
       | x x x 1 0 0 0 |
       | x x x x 1 0 0 |
       | x x x x x 1 0 |
       | x x x x x x 1 |
       
   JK
*/
void skyline::ltv (double *b,double *c)
{
  long i,j,k,lj,uj;
  double s;
  
  nullv (c,n);
  
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1]-1;
    k=i-(uj-lj);
    s=0.0;
    for (j=uj;j>lj;j--){
      s+=a[j]*b[k];  k++;
    }
    c[i]=s+b[i];
  }
}

void skyline::ldlmxv_sky (double *b,double *c)
{
  long i;
  
  utv (b,c);
  for (i=0;i<n;i++){
    c[i]*=a[adr[i]];
  }
  ltv (c,b);
  
  for (i=0;i<n;i++){
    c[i]=b[i];
  }
}


/**
   function adds premultiplied components of %matrix stored in another skyline by c
   to components of actual skyline
   
   @param c - multiplicative coefficient
   @param sky - another skyline
   
   JK
*/
void skyline::addmat_sky (double c,skyline &sky)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]+=c*sky.a[i];
  }
}

void skyline::addmat_gm (double c,gmatrix &gm)
{
  long i;
  
  switch (gm.ts){
  case diag_mat:{
    if (gm.n != n){
      print_err("inconsistent number of rows/columns in matrices",__FILE__,__LINE__,__func__);
      abort ();
    }
    for (i=0;i<n;i++){
      a[adr[i]]+=c*gm.diagm->a[i];
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
   function multiplies skyline by c
   
   @param c - multiplicative coefficient
   
   JK
*/
void skyline::scalmat_sky (double c)
{
  long i;
  for (i=0;i<negm;i++){
    a[i]*=c;
  }
}



/**
   either function decomposes %matrix into L.D.L form
   or computes base vectors of %matrix kernel
   singular or non-singular matrices are accepted
   function is used especially in FETI method

   @param r - base vectors of %matrix kernel
   @param nse - kernel dimension
   @param se - array containing indices of linearly dependent equations
   @param ense - estimation of number of kernel dimension
   @param limit - linear dependency threshold
   @param tc - computation indicator
   
   tc=1 - only %matrix decomposition
   tc=2 - only base vectors computation
   tc=3 - %matrix decomposition and base vectors computation
   
   JK, 25.3.1999
*/
void skyline::ker (double *r,long &nse,long *se,
		   long ense,double limit,long tc)
{
  long i,j,k,ii,jj,kk,lj,uj,li,ui,lk,uk,mi,ise,ib,*adrb;
  double s,g,*b;
  
  b = new double [ense*n];
  adrb = new long [ense+1];
  
  //
  //  matrix decomposition
  //
  if (tc==1 || tc==3){
    if (decompid == 1){
      print_err("factorized matrix has to be factorized once again", __FILE__, __LINE__, __func__);
    }
    else{
      
      //  dependent equation counter
      ise=0;
      ib=0;
      adrb[0]=0;
      
      for (i=1;i<n;i++){
	lj=adr[i];  uj=adr[i+1]-2;
	
	mi=i-(uj-lj)-1;  j=mi+1;
	for (jj=uj;jj>lj;jj--){
	  li=adr[j];  ui=adr[j+1]-1;  k=j-(ui-li);
	  
	  if (k<mi)  {  uk=uj+1;     ii=li+j-mi;  }
	  else       {  uk=lj+i-k;   ii=ui;       }
	  
	  s=0.0;
	  for (kk=uk;kk>jj;kk--){
	    s+=a[kk]*a[ii];  ii--;
	  }
	  
	  a[jj]-=s;  j++;
	}
	
	//  off-diagonal element modification
	s=0.0;  j=mi;
	for (jj=uj+1;jj>lj;jj--){
	  g=a[jj];
	  a[jj]/=a[adr[j]];
	  s+=a[jj]*g;
	  j++;
	}
	a[lj]-=s;
	
	//  diagonal entry check
	if (fabs(a[lj])<limit){
	  //  array containing numbers of dependent equations
	  se[ise]=i;  ise++;
	  
	  if (ise>ense){
	    print_err("larger number of dependent rows than the estimated number", __FILE__, __LINE__, __func__);
	  }
	  
	  
	  //  storage of dependent row and column in array b
	  for (jj=uj+1;jj>lj;jj--){
	    b[ib]=a[jj];  ib++;
	    a[jj]=0.0;
	  }
	  a[lj]=1.0;
	  adrb[ise]=ib;
	  
	  for (j=i+1;j<n;j++){
	    if (j-(adr[j+1]-adr[j])<i)  a[adr[j]+j-i]=0.0;
	  }
	  
	}
	
      }
      
      nse=ise;
      
      decompid=1;
    }
  }
  
  
  //
  //  base vectors of kernel computation
  //
  if (tc==2 || tc==3){
    if (tc==2 && decompid == 0){
      print_err("matrix is not factorized, base vectors of kernel cannot be computed", __FILE__, __LINE__, __func__);
    }
    else{
      //  original entries return
      ise=nse;
      for (i=0;i<ise;i++){
	uj=adr[se[i]+1]-1;  lj=adr[se[i]];
	ib=adrb[i];
	for (jj=uj;jj>lj;jj--){
	  a[jj]=b[ib];  ib++;
	}
      }
      
      //  base vectors assembling
      nullv (r,ise*n);
      for (i=0;i<ise;i++){
	ib=i*n;
	for (j=n-1;j>-1;j--){
	  r[ib+se[i]]=1.0;  s=r[ib+j];
	  uk=adr[j+1]-1;  lk=adr[j];  k=j-(uk-lk);
	  for (kk=uk;kk>lk;kk--){
	    r[ib+k]-=a[kk]*s;  k++;
	  }
	}
      }
      
    }
  }
  
  delete [] adrb;  delete [] b;
}




/**
   function eliminates internal unknowns
   function is used in the Schur complement method and the FETI-DP method
   internal unknowns must be at the beginning
   boundary unknowns must be at the end
   
   @param b - condensed %matrix, it is stored as dense %matrix
   @param c - condensed right hand side %vector
   @param x - %vector of solution
   @param y - right hand side %vector
   @param m - number of boundary unknowns
   @param tc - computation indicator
   
   tc=1 - condensation of the %matrix and the right hand side %vector (forward reduction) only,
          if the matrix is factorized, only modification of the right hand side %vector is performed
   tc=2 - back-substitution only (matrix has to be factorized),
   tc=3 - partial back substitution only (used in the FETI-DP method)
   
   tc=4 - condensation of the right hand side %vector only (matrix has been factorized before)
   
   JK, 23.7.1996
*/
void skyline::ldlkon_sky (double *b,double *c,double *x,double *y,long m,long tc)
{
  long i,j,k,l,ac,ac1,aca,acc,acs,ack,ack1,acrk,ackk,aci,aci1,acri,acii;
  long acj,acj1,acui,acb1,acb2;
  double s,g;
  
  acui=n-m;
  if (tc==1){
    if (decompid == 0){
      //  elimination of rows belonging to internal unknowns
      for (k=1;k<acui;k++){
	ack=adr[k];  ack1=adr[k+1];
	ackk=k+ack;  acrk=ackk-ack1+1;
	acj=ack1-2;
	//  modification of entries from k-th column
	for (i=acrk+1;i<k;i++){
	  aci=adr[i];  aci1=adr[i+1];
	  acii=i+aci;  acri=acii-aci1+1;
	  if (acri<acrk)  ac=acrk;
	  else            ac=acri;
	  acj1=ackk-ac;  acs=acii-ac;
	  s=0.0;
	  for (j=acj1;j>acj;j--){
	    s+=a[j]*a[acs];
	    acs--;
	  }
	  a[acj]-=s;  acj--;
	}
	//  modification of diagonal entry
	s=0.0;
	for (i=ack1-1;i>ack;i--){
	  ackk=adr[acrk];  acrk++;
	  g=a[i];
	  a[i]/=a[ackk];
	  s+=a[i]*g;
	}
	a[ack]-=s;
      }
      
      //  modification of rows belonging to boundary unknowns
      for (k=acui;k<n;k++){
	ack=adr[k];  ack1=adr[k+1]-1;
	ackk=k+ack;  acrk=ackk-ack1;
	acj=ack1-1;
	if (acrk<acui){
	  for (i=acrk+1;i<acui;i++){
	    //  modification of entries in boundary columns and internal rows
	    aci=adr[i];  aci1=adr[i+1];
	    acii=i+aci;  acri=acii-aci1+1;
	    if (acri<acrk)  ac=acrk;
	    else            ac=acri;
	    acj1=ackk-ac;  acs=acii-ac;
	    s=0.0;
	    for (j=acj1;j>acj;j--){
	      s+=a[j]*a[acs];
	      acs--;
	    }
	    a[acj]-=s;  acj--;
	  }
	  aca=ackk-acui;
	  for (i=acui;i<k;i++){
	    //  modification of entries in boundary columns and rows
	    aci=adr[i];  aci1=adr[i+1];
	    acii=i+aci;  acri=acii-aci1+1;
	    if (acri<acrk)  ac=acrk;
	    else            ac=acri;
	    acj=ackk-ac;  acj1=ackk-acui;
	    acs=acii-ac;  s=0.0;
	    for (j=acj;j>acj1;j--){
	      s+=a[j]*a[acs];
	      acs--;
	    }
	    a[aca]-=s;  aca--;
	  }
	  s=0.0;
	  for (i=acrk;i<acui;i++){
	    aci=adr[i];
	    g=a[ack1];
	    a[ack1]/=a[aci];
	    s+=a[ack1]*g;
	    ack1--;
	  }
	  a[ack]-=s;
	}
      }
      
      decompid=1;
    }
    
    //  assembling of condensed matrix
    acc=0;
    for (k=acui;k<n;k++){
      ack=adr[k];  ack1=adr[k+1]-1;
      acrk=k+ack-ack1;
      if (acrk<acui)  aci=acui;
      else            aci=acrk;
      aci1=ack+k-aci;
      acb1=acc*m+acc;  acb2=acb1;
      for (i=ack;i<=aci1;i++){
	b[acb1]=a[i];  b[acb2]=a[i];
	acb1-=m;  acb2--;
      }
      acc++;
    }
    
    
    //  modification of right hand side (internal part)
    //  L.z = y
    for (k=0;k<acui;k++){
      ack=adr[k];  ack1=adr[k+1]-1;
      ackk=k+ack;  acrk=ackk-ack1;
      s=0.0;
      for (i=ack1;i>ack;i--){
	s+=a[i]*y[acrk];  acrk++;
      }
      y[k]-=s;
    }
    for (k=acui;k<n;k++){
      ack=adr[k];  ack1=adr[k+1]-1;
      ackk=k+ack;  acrk=ackk-ack1;
      s=0.0;  aci=ackk-acui;
      for (i=ack1;i>aci;i--){
	s+=a[i]*y[acrk];  acrk++;
      }
      y[k]-=s;
    }
    
    //  assembling of right hand side
    l=0;
    for (k=acui;k<n;k++){
      c[l]=y[k];
      l++;
    }
    
  }
  
  
  if (tc==2){
    if (decompid == 0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{

      //  interface components are localized to the vector of subdomain
      j=0;
      for (k=acui;k<n;k++){
	x[k]=c[j];  j++;
      }

      //  back-substitution
      //
      //  D.L.x = z  => L.x = D^{-1}.z 
      for (k=0;k<acui;k++){
	ack=adr[k];
	y[k]/=a[ack];
      }
      //  modification of right hand side
      for (k=n-1;k>=acui;k--){
	ack=adr[k];  ack1=adr[k+1]-1;
	ackk=k+ack;  acrk=ackk-ack1;
	s=x[k];
	for (i=ack1;i>ack;i--){
	  y[acrk]-=a[i]*s;
	  acrk++;
	}
      }
      //  computation of solution
      for (k=acui-1;k>-1;k--){
	ack=adr[k];  ack1=adr[k+1]-1;
	ackk=k+ack;  acrk=ackk-ack1;
	x[k]=y[k];  s=x[k];
	for (i=ack1;i>ack;i--){
	  y[acrk]-=s*a[i];  acrk++;
	}
      }
    }
  }
  
  if (tc==3){
    if (decompid == 0){
      print_err("matrix is not factorized, back substitution cannot be performed", __FILE__, __LINE__, __func__);
    }
    else{
      //
      //  partial back-substitution
      //
      //  computation of L.z = y  => z
      for (k=1;k<acui;k++){
	ack=adr[k]+1;  ack1=adr[k+1];
	ac1=k-1;  s=0.0;
	for (i=ack;i<ack1;i++){
	  s+=a[i]*y[ac1];
	  ac1--;
	}
	y[k]-=s;
      }
      
      
      //  D.L.x = z  =>  L.x = D^{-1}.z
      for (k=0;k<acui;k++){
	ac=adr[k];
	//if (fabs(a[ac])<zero){
	//printf ("\n\n zero diagonal entry in the skyline (%ld-th row) (D^-1)",k);
	//}
	y[k]/=a[ac];
      }
      //  L.x = D^{-1}.z => x
      for (k=acui-1;k>-1;k--){
	ack=adr[k]+1;  ack1=adr[k+1];
	x[k]=y[k];  g=x[k];
	ac=k-1;
	for (i=ack;i<ack1;i++){
	  y[ac]-=a[i]*g;
	  ac--;
	}
      }
    }
  }
  
  if (tc==4){
    if (decompid == 0){
      print_err ("matrix is not factorized, forward reduction of the right hand side cannot be performed", __FILE__, __LINE__, __func__);
    }
    
    //  modification of right hand side (internal part)
    //  L.z = y
    for (k=0;k<acui;k++){
      ack=adr[k];  ack1=adr[k+1]-1;
      ackk=k+ack;  acrk=ackk-ack1;
      s=0.0;
      for (i=ack1;i>ack;i--){
	s+=a[i]*y[acrk];  acrk++;
      }
      y[k]-=s;
    }
    for (k=acui;k<n;k++){
      ack=adr[k];  ack1=adr[k+1]-1;
      ackk=k+ack;  acrk=ackk-ack1;
      s=0.0;  aci=ackk-acui;
      for (i=ack1;i>aci;i--){
	s+=a[i]*y[acrk];  acrk++;
      }
      y[k]-=s;
    }
    
    //  assembling of right hand side
    l=0;
    for (k=acui;k<n;k++){
      c[l]=y[k];
      l++;
    }
    
  }

}



/**
   function eliminates internal unknowns
   function is used in the Schur complement method and the FETI-DP method
   internal unknowns must be at the beginning
   boundary unknowns must be at the end
   
   @param b - condensed %matrix, it is stored as dense %matrix
   @param c - condensed right hand side %vector
   @param x - %vector of solution
   @param y - right hand side %vector
   @param m - number of boundary unknowns
   @param tc - computation indicator
   
   tc=1 - condensation of the %matrix and the right hand side %vector (forward reduction) only,
          if the matrix is factorized, only modification of the right hand side %vector is performed
   tc=2 - back-substitution only (matrix has to be factorized),
   tc=3 - partial back substitution only (used in the FETI-DP method)
   
   tc=4 - condensation of the right hand side %vector only (matrix has been factorized before)
   
   JK, 23.7.1996
*/
double skyline::ldlkoncount_sky (double */*b*/,double */*c*/,double */*x*/,double */*y*/,long m,long tc)
{
  long i,j,k,ac,aca,acs,ack,ack1,acrk,ackk,aci,aci1,acri,acii;
  long acj,acj1,acui;
  double s;
  
  double no;
  
  no=0.0;
  
  acui=n-m;
  if (tc==1){
    if (decompid == 0){
      //  elimination of rows belonging to internal unknowns
      for (k=1;k<acui;k++){
	ack=adr[k];  ack1=adr[k+1];
	ackk=k+ack;  acrk=ackk-ack1+1;
	acj=ack1-2;
	//  modification of entries from k-th column
	for (i=acrk+1;i<k;i++){
	  aci=adr[i];  aci1=adr[i+1];
	  acii=i+aci;  acri=acii-aci1+1;
	  if (acri<acrk)  ac=acrk;
	  else            ac=acri;
	  acj1=ackk-ac;  acs=acii-ac;
	  s=0.0;
	  for (j=acj1;j>acj;j--){
	    no+=1.0;
	    acs--;
	  }
	  acj--;
	}
	//  modification of diagonal entry
	s=0.0;
	for (i=ack1-1;i>ack;i--){
	  ackk=adr[acrk];  acrk++;
	  no+=1.0;
	}
      }
      
      fprintf (stdout,"\n number of operations after inner unknowns elimination  %le",no);

      //  modification of rows belonging to boundary unknowns
      for (k=acui;k<n;k++){
	ack=adr[k];  ack1=adr[k+1]-1;
	ackk=k+ack;  acrk=ackk-ack1;
	acj=ack1-1;
	if (acrk<acui){
	  for (i=acrk+1;i<acui;i++){
	    //  modification of entries in boundary columns and internal rows
	    aci=adr[i];  aci1=adr[i+1];
	    acii=i+aci;  acri=acii-aci1+1;
	    if (acri<acrk)  ac=acrk;
	    else            ac=acri;
	    acj1=ackk-ac;  acs=acii-ac;
	    s=0.0;
	    for (j=acj1;j>acj;j--){
	      no+=1.0;
	      acs--;
	    }
	    acj--;
	  }
	  aca=ackk-acui;
	  for (i=acui;i<k;i++){
	    //  modification of entries in boundary columns and rows
	    aci=adr[i];  aci1=adr[i+1];
	    acii=i+aci;  acri=acii-aci1+1;
	    if (acri<acrk)  ac=acrk;
	    else            ac=acri;
	    acj=ackk-ac;  acj1=ackk-acui;
	    acs=acii-ac;  s=0.0;
	    for (j=acj;j>acj1;j--){
	      no+=1.0;
	      acs--;
	    }
	    aca--;
	  }
	  s=0.0;
	  for (i=acrk;i<acui;i++){
	    aci=adr[i];
	    no+=1.0;
	    ack1--;
	  }
	}
      }
      
      decompid=1;
    }
    
  }
  
  return no;
}



/**
   function computes %vector y = K^+ x   where K^+ is pseudoinverse %matrix
   it is used in FETI method
   %matrix must be decomposed into L.D.L^T form
   
   @param x - left hand side (solution of the problem)
   @param y - right hand side
   @param nse - number of linearly dependent rows
   @param se - indecies of linearly dependent rows
   @param zero - computer zero
   
   JK, 8.3.2002
*/
void skyline::ldl_feti_sky (double *x,double *y,long nse,long *se,double zero)
{
  long i,j,k,ii,lj,uj;
  double s;
  
  //  computation of z from equation L . z = y
  //  vector z overwrites vector y
  k=0;
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1]-1;  ii=i-uj+lj;  s=0.0;
    for (j=uj;j>lj;j--){
      s+=a[j]*y[ii];  ii++;
    }
    
    if (se[k]==i){
      y[i]=0.0;  k++;
    }
    else{
      y[i]-=s;
    }
  }

  //  vector modification D . (L.x) = z => L.x = D^{-1}.z
  for (i=0;i<n;i++){
    s=a[adr[i]];
    if (fabs(s)<zero)  fprintf (stderr,"\n\n zero diagonal entry (row number %ld) in function skyline::ldl_feti_sky\n",i);
    y[i]/=s;
  }
  
  //  computation of vector x from equation L . x = D^{-1}.z
  k=nse-1;
  for (i=n-1;i>=0;i--){
    lj=adr[i];  uj=adr[i+1];
    if (k>-1){
      if (se[k]==i){
	//y[i]=1.0;  k--;
	y[i]=0.0;  k--;
      }
    }
    x[i]=y[i];
    s=x[i];  ii=i-1;
    for (j=lj+1;j<uj;j++){
      y[ii]-=s*a[j];  ii--;
    }
  }

}

/**
   function assembles block 12 of %matrix
   
       | A_11 A_12 |
   A = |           |
       | A_21 A_22 |
       
   @param block - %matrix containing required block
   @param nrdof - number of uneliminated unknowns
   
   JK, 24.7.2002
*/
void skyline::ldl_a12block (double *block,long nrdof)
{
  long i,j,ai,ai1,li,lj,ri,ci,minri;
  
  nullv (block,nrdof*(n-nrdof));
  
  li=n-nrdof;
  
  for (i=li;i<n;i++){
    ai=adr[i];  ai1=adr[i+1];
    minri=ai1-ai-i+1;
    if (minri>=li)
      continue;
    lj=ai+i-li+1;
    ri=li-1;
    ci=i-li;
    for (j=lj;j<ai1;j++){
      block[ri*nrdof+ci]=a[j];
      ri--;
    }
  }
  
}

/**
   function prints %matrix into output file
   
   @param out - output stream
   
   JK
*/
void skyline::printmat (FILE *out)
{
  long i;
  
  FILE *mat;
  mat = fopen ("matice.txt","w");
  fprintf (mat,"%ld",n);
  fprintf (out,"\n\n");
  for (i=0;i<n+1;i++){
    fprintf (out,"\n %9ld",adr[i]);
  }
  fprintf (out,"\n\n");
  for (i=0;i<negm;i++){
    fprintf (out,"\n %20.10f",a[i]);
  }
  fclose (mat);
  
  
  /*
  fprintf (out,"%ld",n);
  for (i=0;i<n;i++){
    fprintf (out,"\n");
    for (j=0;j<i;j++){
      if (adr[i]+i-j>=adr[i+1])  fprintf (out," %10.5le",0.0);
      else fprintf (out," %10.5le",a[adr[i]+i-j]);
    }
    for (j=i;j<n;j++){
      if (adr[j]+j-i>=adr[j+1])  fprintf (out," %10.5le",0.0);
      else fprintf (out," %10.5le",a[adr[j]+j-i]);
    }
  }
  */
  
  /*
  FILE *mat;
  mat = fopen ("matice.txt","w");
  fprintf (mat,"%ld",n);
  for (i=0;i<n;i++){
    fprintf (mat,"\n");
    for (j=0;j<i;j++){
      if (adr[i]+i-j>=adr[i+1])  fprintf (mat," %12.10le",0.0);
      else fprintf (mat," %12.10le",a[adr[i]+i-j]);
    }
    for (j=i;j<n;j++){
      if (adr[j]+j-i>=adr[j+1])  fprintf (mat," %12.10le",0.0);
      else fprintf (mat," %12.10le",a[adr[j]+j-i]);
    }
  }
  fclose (mat);
  */
  
  
  /*
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (out,"%20.15e\n",a[j]);
    }
  }
  */
  
  
  /*
  FILE *rif,*cif,*mat,*ds;
  rif = fopen ("ri.txt","w");
  cif = fopen ("ci.txt","w");
  mat = fopen ("mat.txt","w");
  ds = fopen ("ds.txt","w");
  for (i=0;i<n;i++){
    for (j=adr[i];j<adr[i+1];j++){
      fprintf (rif,"%ld ",i-j+adr[i]+1);
      fprintf (cif,"%ld ",i+1);
      fprintf (mat,"%le ",a[i]);
    }
    fprintf (ds,"%ld %ld\n",i,adr[i+1]-adr[i]);
  }
  fclose (rif);
  fclose (cif);
  fclose (mat);
  fclose (ds);
  */
  
  
  
}

/**
   function prints diagonal components of the %matrix
   
   @param out - output stream
   
   JK
*/
void skyline::printdiag(FILE *out)
{
  fprintf (out,"\n\n");
  for (long i=0;i<n;i++)
    fprintf(out, "%5ld  %.10e\n", i,a[adr[i]]);
}

/**
   function returns required %matrix entry
   
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
double skyline::give_entry (long ri,long ci)
{
  long i;
  double e;
  
  if (ci>ri){
    i=adr[ci];
    e = a[i+ci-ri];
  }
  else{
    i=adr[ri];
    e = a[i+ri-ci];
  }
  
  return e;
}

/**
   function adds required %matrix entry
   
   @param e - %matrix entry
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
void skyline::add_entry (double e,long ri,long ci)
{
  long i;
  
  if (ci>ri){
    i=adr[ci]+ci-ri;
    a[i]+=e;
  }
  else{
    i=adr[ri]+ri-ci;
    a[i]+=e;
  }
}

/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long skyline::give_negm ()
{
  return negm;
}

/**
   function scales %matrix by its diagonal elements

   JK, 23.5.2008
*/
void skyline::diag_scale (double *d)
{
  long i,j,k;
  
  for (i=0;i<n;i++){
    d[i]=1.0/sqrt(a[adr[i]]);
  }
  
  for (i=0;i<n;i++){
    k=i;
    for (j=adr[i];j<adr[i+1];j++){
      a[j]*=d[i]*d[k];
      k--;
    }
  }
}

/*
  JB
*/

void skyline::select_submatrix(skyline *smsky,long nsdof,long *sdof,FILE *out)
{
  long i,j,k,l,address1,address2,upper,lower;
   
  smsky->n = nsdof;
  smsky->adr = new long[nsdof+1];
  fprintf(out,"kontrola v select_submatrix\n");
  
  for(i = 0; i < nsdof; i++){
    fprintf(out,"%ld\n",sdof[i]);
    smsky->adr[i]=-1;
  }
  
  
  // jen pro test
  // ->
  // velka A
  fprintf(out,"negm = %ld n = %ld\n",negm,n);
  for(i = 0; i < negm; i++){
    fprintf(out,"%le\n",a[i]);
  }
  // adr
  for(i = 0; i < n+1; i++){
    fprintf(out,"%ld\n",adr[i]);
  }
  // <-
  
  // smsky->negm - pocet prvku v male matici
  smsky->negm = 0;
  for(i = 0; i < nsdof; i++){
    address1 =  adr[sdof[i]];
    address2 =  adr[sdof[i]+1];
    k = address2 - address1;
    upper = sdof[i]+1;
    lower = upper-k;
    l = 0;
    //printf("sdof %ld a1 %ld a2 %ld k %ld l %ld u %ld\n",sdof[i],address1,address2,k,lower,upper);
    for(j = 0; j < nsdof; j++){
      if(sdof[j] < upper && sdof[j] >= lower){
	//printf("%ld\n",sdof[j]);
	smsky->negm++;
	l++;
      }
    }
    smsky->adr[i]=smsky->negm-l;
  }
  
  fprintf(out,"smsky->negm = %ld\n",smsky->negm); 
  smsky->adr[nsdof]= smsky->negm;
  for(i = 0; i < nsdof+1; i++){
    fprintf(out,"%ld\n",smsky->adr[i]); 
  }
  
  smsky->a = new double[smsky->negm];
  
  for(i = 0; i < nsdof; i++){
    address1 =  adr[sdof[i]];
    address2 =  adr[sdof[i]+1];
    k = address2 - address1;
    upper = sdof[i]+1;
    lower = upper-k;
    //printf("sdof %ld a1 %ld a2 %ld k %ld l %ld u %ld\n",sdof[i],address1,address2,k,lower,upper);
    for(j = 0; j < nsdof; j++){
      if(sdof[j] < upper && sdof[j] >= lower){
	// sdof[j] - vyznam cisla radku sdof[i] - vyznam cisla sloupce - velka matice
	// i cislo sloupce j - cislo radku - mala matice
	//printf("adr = %ld j = %ld i = %ld a = %ld\n",address1,sdof[j],sdof[i],address1+sdof[i]-sdof[j]);
	//printf("sadr = %ld j = %ld i = %ld a = %ld\n",sadr[i],j,i,sadr[i]+i-j);
	smsky->a[smsky->adr[i]+i-j]=a[address1+sdof[i]-sdof[j]];
	//printf("%lf\n",sa[sadr[i]+sdof[i]-sdof[j]]);
      }
    }
  }
  
  // jen pro test
  // ->
  // velka A
  fprintf(out,"smsky negm = %ld n = %ld\n",smsky->negm,smsky->n);
  for(i = 0; i < smsky->negm; i++){
    fprintf(out,"%le\n",smsky->a[i]);
  }
  // adr
  for(i = 0; i < smsky->n+1; i++){
    fprintf(out,"%ld\n",smsky->adr[i]);
  }
  // <-
  
  
}


/**
   function checks diagonal entries
   
   the function is used in some nonlinear nostationary problems
   where high jumps in coefficients occur
   some of element matrices are zero matrices and this
   function puts nonzero values on the diagonal
   
   @param thr - prescribed threshold
   
   JK, 25.8.2011
*/
void skyline::diag_check (double thr)
{
  long i,j;
  
  for (i=0;i<n;i++){
    j=adr[i];
    if (fabs(a[j])<thr)
      a[j]=1.0;
  }
}
