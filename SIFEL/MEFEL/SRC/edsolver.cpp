#include <math.h>
#include "edsolver.h"
#include "global.h"
#include "probdesc.h"
#include "globmat.h"
#include "mechprint.h"
#include "elemswitch.h"

#include "gmatrix.h"
#include "eigsol.h"
#include "vector.h"


void solve_eigen_dynamics (double *x,double *w)
{
  long i;
  
  //  stiffness matrix
  stiffness_matrix (0);
  //  mass matrix
  mass_matrix (0);
  
  switch (Mp->eigsol.teigsol){
  case inv_iteration:{
    inverse_iteration (x,w);
    break;
  }
  case subspace_it_jacobi:{
    subspace_iter_jac (x,w);
    break;
  }
  case subspace_it_gsortho:{
    subspace_iter_ortho (x,w);
    break;
  }
  case shifted_subspace_it_gsortho:{
    subspace_shift_iter_ortho (x,w);
    break;
  }  default:{
    print_err("unknown eigenvalue solver is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }
  
  print_eigenvalues (w);
  print_eigenvectors ();
  //print_eigenvect_martin (Out);
  
  
  print_init(-1, "wt");    
  for (i=0;i<Mp->eigsol.neigv;i++){
    //  computes and prints required quantities
    //compute_req_val(i);
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    //compute_req_val (i);
    print_step(i, i+1, w[i], NULL);    
  }
  print_close();

  //Mm->allipstresses (0);
  
  
}

/**
   function computes the smallest eigenvalue and appropriate eigenvector
   of the generalized problem Ax=wBx
   the eigenvector of the smallest eigenvalue is stored in Lsrs in lhs

   @param w - the smallest eigenvalue
   
   3.11.2001
*/
void inverse_iteration (double *x,double *w)
{
  long i,ni,n;
  double nom,denom,zero,rho,prevrho,err,error=0.0;
  double *p,*z;
  
  //  maximum number of iterations
  ni=Mp->eigsol.nies;
  //  required error
  err=Mp->eigsol.erres;
  //  number of unknowns
  n=Ndofm;
  //  computer zero
  zero=Mp->zero;
  
  //Smat->tlinsol=Mp->tlinsol;

  p = new double [n];
  z = new double [n];
  //x = Lsrs->give_lhs(0);
  
  //  LDL decomposition of matrix A
  //Sm_sky->ldl_sky (z,z,zero,2);
  
  //  initial value of vector x
  for (i=0;i<n;i++){
    x[i]=1.0;
  }
  //  initial vector z
  //massxvect (x,z);
  Mmat->gmxv (x,z);

  //  main iteration loop
  for (i=0;i<ni;i++){
    //  copy of vector z which will be overwriten in back substitution
    copyv (z,p,n);
    //  back substitution
    //Sm_sky->ldl_sky (x,p,zero,3);
    //Smat->solve_system (Gtm,x,p);
    Mp->ssle->solve_system (Gtm,Smat,x,p,Out);
    
    //  x.K.x
    nom=ss(x,z,n);
    //  new vector z
    //massxvect (x,z);
    Mmat->gmxv (x,z);

    //  x.M.x
    denom=ss(x,z,n);
    rho=nom/denom;
    denom=1.0/sqrt(denom);
    //  normed vector z
    cmulv (denom,z,n);

    if (Mespr==1)  printf ("\n iterace %4ld     rho %e",i,rho);
    if (i!=0){
      error=(prevrho-rho)/prevrho;
      if (error<err){
	cmulv (denom,x,n);  break;
      }
    }
    prevrho=rho;
  }
  Mp->eigsol.anies=i;
  Mp->eigsol.aerres=error;
  w[0]=rho;

  delete [] z;  delete [] p;
}

/**
   function computes eigenvalues and eigenvectors
   subspace iteration with Gram-Schmidt orthonormalization is used
   
   @param w - array containing eigenvalues
   
   4.11.2001
*/
void subspace_iter_ortho (double *x,double *w)
{
  long i,j,k,l,ii,ni,n,nv,nrv;
  double err,error,maxerror=0.0,nom,denom,alpha;
  double *z,*p,*wo;
  
  //  problem size
  n=Ndofm;
  //  maximum number of iterations
  ni=Mp->eigsol.nies;
  //  required error
  err=Mp->eigsol.erres;
  //  number of required eigenvectors and eigenvalues
  nrv=Mp->eigsol.neigv;
  //  number of vectors used for iteration process
  nv=Mp->eigsol.nev;
  
  z = new double [n];
  p = new double [n];
  wo = new double [nv];

  //  initial value of vector x
  nom=1.0;
  for (i=0;i<nv*n;i++){
    x[i]=nom;
    nom+=1.0;
  }
  
  
  //  main iteration loop
  for (i=0;i<ni;i++){
    for (j=0;j<nv;j++){
      //  right hand side vector
      nullv (z,n);
      Mmat->gmxv (x+j*n,z);

      copyv (z,p,n);
      //  back substitution
      Mp->ssle->solve_system (Gtm,Smat,x+j*n,z,Out);

      //  Rayleigh' quotients
      //  x.K.x
      nom=ss(x+j*n,p,n);
      //  x.M.x
      nullv (z,n);
      Mmat->gmxv (x+j*n,z);

      denom=ss(x+j*n,z,n);
      w[j]=nom/denom;
      //printf ("\n iterace  %4ld %ld   %30.20e",i,j,w[j]);
    }
    
    //  Gram-Schmidt orthonormalization
    for (j=0;j<nv;j++){
      nullv (z,n);
      Mmat->gmxv (x+j*n,z);
      
      for (k=0;k<j;k++){
	alpha=ss(x+k*n,z,n);
	for (l=0;l<n;l++){
	  x[j*n+l]-=alpha*x[k*n+l];
	}
      }
      nullv (z,n);
      Mmat->gmxv (x+j*n,z);
      
      alpha=ss(x+j*n,z,n);
      alpha=1.0/sqrt(alpha);
      cmulv (alpha,x+j*n,n);
    }
    
    //  convergence check
    if (i!=0){
      ii=0;  maxerror=0.0;
      for (j=0;j<nrv;j++){
	error=fabs(wo[j]-w[j])/w[j];
	if (error>maxerror)  maxerror=error;
	if (error<err)  ii++;
      }
      
      if (Mespr==1)  printf ("\n iterace  %4ld    pocet doiterovanych tvaru %4ld  max. odchylka %e",i,ii,maxerror);
      //if (Mespr==1)  printf ("\n iterace  %4ld    w[0] %30.20e",i,w[0]);
      
      if (ii==nrv){
	//  for fuzzy computation
	//for (j=0;j<nrv;j++){
	//max=0.0;
	//for (k=0;k<n;k++){
	//if (fabs(x[j*n+k])>max)  max = x[j*n+k];
	//}
	//if (max<0.0){
	//for (k=0;k<n;k++){
	//x[j*n+k]*=-1.0;
	//}
	//}
	//}
	
	
	//  modification for fuzzy computation
	/*
	for (j=0;j<nrv;j++){
	  if (x[j*n]<0.0){
	    for (k=0;k<n;k++){
	      x[j*n+k]*=-1.0;
	    }
	  }
	}
	*/
	
	break;
      }
    }
    
    for (j=0;j<nv;j++){
      wo[j]=w[j];
    }
  }

  Mp->eigsol.anies=i;
  Mp->eigsol.aerres=maxerror;
  
  /*
  //  kontrola
  //  stiffness matrix
  stiffness_matrix (0);
  //  mass matrix
  mass_matrix (0);
  
  double *jk,*kj;
  jk = new double [n];
  kj = new double [n];
  Smat->gmxv (x,jk);
  Mmat->gmxv (x,kj);
  cmulv(w[0],kj,n);
  
  fprintf (Out,"\n\n kontrola jk a kj \n");
  for (i=0;i<n;i++){
    fprintf (Out,"\n %le %le  %le",jk[i],kj[i],jk[i]-kj[i]);
  }
  fprintf (Out,"\n\n konec kontroly\n");
  */



  delete [] z;
  delete [] wo;
  delete [] p;
}

/**
   function is not ready for computation
   initial vectors must be assembled
   
   function computes eigenvalues and eigenvectors
   subspace iteration with Jacobi method of rotations is used

   @param w - array containing eigenvalues

   4.11.2001
*/
void subspace_iter_jac (double *x,double *w)
{
  long i,j,k,l,ii,ni,n,nv,nrv,nij,njacthr,stop;
  double err,error,maxerror,zero,alpha;
  double *jacthr;
  //double *z, *p, *wo, *redstiff, *redmass, *redeigv, *aux, *invec;
  vector z, p, wo, redstiff, redmass, redeigv, aux, invec;
  // long *ind
  ivector ind;
  
  //  problem size
  n=Ndofm;
  //  maximum number of iterations
  ni=Mp->eigsol.nies;
  //  required error
  err=Mp->eigsol.erres;
  //  number of required eigenvectors and eigenvalues
  nrv=Mp->eigsol.neigv;
  //  number of vectors used for iteration process
  nv=Mp->eigsol.nev;
  //  maximum number of iterations in Jacobi method
  nij=Mp->eigsol.nijmr;
  //  number of thresholds in Jacobi' method
  njacthr=Mp->eigsol.njacthr;
  //  computer zero
  zero=Mp->zero;
  
  stop=0;

  jacthr=Mp->eigsol.jacthr;
  //x=Lsrs->give_lhs(0);

  //z = new double [n*nv];
  reallocv(n*nv, z);
  //p = new double [n*nv];
  reallocv(n*nv, p);
  //wo = new double [nv];
  reallocv(nv, wo);
  //redstiff = new double [nv*nv];
  reallocv(nv*nv, redstiff);
  //redmass = new double [nv*nv];
  reallocv(nv*nv, redmass);
  //redeigv = new double [nv*nv];
  reallocv(nv*nv, redeigv);
  //aux = new double [nv];
  reallocv(nv, aux);
  //invec = new double [n];
  reallocv(n, invec);
  //ind = new long [n];
  reallocv(n, ind);

  for (i=0;i<n;i++){
    //invec[i]=Mm_sky->a[Mm_sky->adr[i]]/Sm_sky->a[Sm_sky->adr[i]];
  }
  for (i=0;i<n;i++){
    alpha=invec[i];  k=i;
    for (j=i+1;j<n;j++){
      if (alpha<invec[j]){
	alpha=invec[j];  k=j;
      }
    }
    ind[i]=k;
    alpha=invec[i];
    invec[i]=invec[k];
    invec[k]=alpha;
  }
  
  //  LDL decomposition of matrix A
  //Sm_sky->ldl_sky (z.a,z.a,zero,2);
  //Smat->solve_system (Gtm,z.a,z.a);
  //Mp->ssle->solve_system (Gtm,Smat,z.a,z.a,Out);

  //  initial value of vector x
  for (i=0;i<nv*n;i++){
    x[i]=0.0;
  }
  for (i=0;i<n;i++){
    x[i]=1.0;
  }
  for (i=1;i<nv;i++){
    x[i*n+ind[i]]=1.0;
  }

  //  right hand side vector
  for (j=0;j<nv;j++){
    //massxvect (x+j*n,z+j*n);
    Mmat->gmxv (x+j*n,z.a+j*n);
  }
  
  /*
  fprintf (Out,"\n\n\n Matice hmotnosti");
  for (i=0;i<n;i++){
    fprintf (Out,"\n");
    for (j=Mm_sky->adr[i];j<Mm_sky->adr[i+1];j++){
      fprintf (Out," %lf",Mm_sky->a[j]);
    }
  }
  fprintf (Out,"\n\n slozky soucinu M.x=z");
  for (i=0;i<nv;i++){
    fprintf (Out,"\n");
    for (j=0;j<n;j++){
      fprintf (Out," %lf",z[i*n+j]);
    }
  }
  */

  //  main iteration loop
  for (i=0;i<ni;i++){
    
    //  A.x_new = z = B.x_old => x_new
    for (j=0;j<nv;j++){
      //  backup of B.x_old (will be overwritten)
      copyv (z.a+j*n,p.a+j*n,n);
      //  back substitution - solution x_new
      //Sm_sky->ldl_sky (x+j*n,z.a+j*n,zero,3);
      //Smat->solve_system (Gtm,x+j*n,z.a+j*n);
      Mp->ssle->solve_system (Gtm,Smat,x+j*n,z.a+j*n,Out);
    }

    
    /*
    fprintf (Out,"\n\n slozky x");
    for (long kk=0;kk<nv;kk++){
      fprintf (Out,"\n");
      for (j=0;j<n;j++){
	fprintf (Out," %lf",x[kk*n+j]);
      }
    }
    */


    //  reduced stiffness matrix computaion
    //  K_red = x.K.x = x.z
    mtxmccr (x,p.a,redstiff.a,n,nv,nv);
    
    //  new right hand side - z = B.x
    for (j=0;j<nv;j++){
      //massxvect (x+j*n,z.a+j*n);
      Mmat->gmxv (x+j*n,z.a+j*n);
    }
    
    //  reduced mass matrix computation
    //  M_red = x.M.x = x.z
    mtxmccr (x,z.a,redmass.a,n,nv,nv);
    
    /*
    fprintf (Out,"\n\n\n iterace cislo %ld",i);
    for (long kk=0;kk<nv;kk++){
      fprintf (Out,"\n");
      for (long jj=0;jj<nv;jj++){
	fprintf (Out,"\n %4ld %4ld k %30.15lf     m %30.15lf",kk,jj,redstiff[kk*nv+jj],redmass[kk*nv+jj]);
      }
    }
    */

    
    //  solution of reduced system by Jacobi' method of rotations
    gen_jacobi (redstiff.a,redmass.a,redeigv.a,w,nv,nij,jacthr,njacthr,zero);
    
    /*
    for (long kk=0;kk<nv;kk++){
      fprintf (Out,"\n w %30.15lf",w[kk]);
      for (long jj=0;jj<nv;jj++){
	fprintf (Out,"\n %4ld %4ld x %30.15lf",kk,jj,redeigv[jj*nv+kk]);
      }
    }
    */

    //  updated right hand side - z = p.redeigv
    for (j=0;j<n;j++){
      for (k=0;k<nv;k++){
	aux[k]=0.0;
	for (l=0;l<nv;l++){
	  aux[k]+=z[l*n+j]*redeigv[l*nv+k];
	}
      }
      for (k=0;k<nv;k++){
	z[k*n+j]=aux[k];
      }
    }
    
    
    
    //  normalization process
    double norm = sqrt(ss (z.a,z.a,n*nv));
    cmulv (1.0/norm,z.a,n*nv);

    
    //  convergence check
    if (i!=0){
      ii=0;  maxerror=0.0;
      for (j=0;j<nrv;j++){
	error=(wo[j]-w[j])/w[j];
	if (error>maxerror)  maxerror=error;
	if (error<err)  ii++;
      }
      
      if (Mespr==1)  fprintf (stdout,"\n iterace  %4ld    pocet doiterovanych tvaru %4ld  max. odchylka %e",i,ii,maxerror);

      if (ii==nrv)  stop=1;
    }
    
    for (j=0;j<nv;j++){
      wo[j]=w[j];
    }
    
    if (stop==1)  break;
  }

  Mp->eigsol.anies=i;
  Mp->eigsol.aerres=maxerror;
}


/**
   function computes eigenvalues and eigenvectors
   subspace iteration with Gram-Schmidt orthonormalization is used
   shift is applied in order to solve systems with singular stiffness %matrix
   
   @param w - array containing eigenvalues
   
   6. 2. 2014
*/
void subspace_shift_iter_ortho (double *x,double *w)
{
  long i,j,k,l,ii,ni,n,nv,nrv;
  double err,error,maxerror=0.0,zero,nom,denom,alpha,shift;
  double *z,*p,*wo;
  
  //  problem size
  n=Ndofm;
  //  maximum number of iterations
  ni=Mp->eigsol.nies;
  //  required error
  err=Mp->eigsol.erres;
  //  number of required eigenvectors and eigenvalues
  nrv=Mp->eigsol.neigv;
  //  number of vectors used for iteration process
  nv=Mp->eigsol.nev;
  //  computer zero
  zero=Mp->zero;
  //  shift
  shift = Mp->eigsol.shift;
  
  z = new double [n];
  p = new double [n];
  wo = new double [nv];
  
  for (i=0;i<nv;i++){
    wo[i]=0.0;
  }

  //  initial value of vector x
  nom=1.0;
  for (i=0;i<nv*n;i++){
    x[i]=nom;
    nom+=1.0;
  }
  
  //  shift of the matrix
  Smat->addgm (shift,*Mmat);
  
  //  main iteration loop
  for (i=0;i<ni;i++){
    for (j=0;j<nv;j++){
      //  right hand side vector
      Mmat->gmxv (x+j*n,z);

      copyv (z,p,n);
      //  back substitution
      //Sm_sky->ldl_sky (x+j*n,z,zero,3);
      //Smat->solve_system (Gtm,x+j*n,z);
      Mp->ssle->solve_system (Gtm,Smat,x+j*n,z,Out);

      //  Rayleigh' quotients
      nom=ss(x+j*n,p,n);
      //massxvect (x+j*n,z);
      Mmat->gmxv (x+j*n,z);

      denom=ss(x+j*n,z,n);
      w[j]=nom/denom;
    }
    
    //  Gram-Schmidt orthonormalization
    for (j=0;j<nv;j++){
      //massxvect (x+j*n,z);
      Mmat->gmxv (x+j*n,z);
      
      for (k=0;k<j;k++){
	alpha=ss(x+k*n,z,n);
	for (l=0;l<n;l++){
	  x[j*n+l]-=alpha*x[k*n+l];
	}
      }
      //massxvect (x+j*n,z);
      Mmat->gmxv (x+j*n,z);

      alpha=ss(x+j*n,z,n);
      alpha=1.0/sqrt(alpha);
      cmulv (alpha,x+j*n,n);
    }
    
    //  convergence check
    //if (i>5){
    
    ii=0;  maxerror=0.0;
    for (j=0;j<nrv;j++){
      error=fabs(wo[j]-w[j])/w[j];
      //error=fabs(wo[j]-w[j]);
      if (error>maxerror)
        maxerror=error;
      if (error<err)
        ii++;
    }
    
    
    
    if (Mespr==1)  printf ("\n iterace  %4ld    pocet doiterovanych tvaru %4ld  max. odchylka %e",i,ii,maxerror);
    //if (Mespr==1)  printf ("\n iterace  %4ld    w[0] %30.20e",i,w[0]);
    
    if (maxerror<err && i>4)
      break;
    
      //if (ii==nrv){
	//  for fuzzy computation
	//for (j=0;j<nrv;j++){
	//max=0.0;
	//for (k=0;k<n;k++){
	//if (fabs(x[j*n+k])>max)  max = x[j*n+k];
	//}
	//if (max<0.0){
	//for (k=0;k<n;k++){
	//x[j*n+k]*=-1.0;
	//}
	//}
	//}
	
	
	//  modification for fuzzy computation
	/*
	for (j=0;j<nrv;j++){
	  if (x[j*n]<0.0){
	    for (k=0;k<n;k++){
	      x[j*n+k]*=-1.0;
	    }
	  }
	}
	*/
	
	//break;
     // }
   // }
    
    for (j=0;j<nv;j++){
      wo[j]=w[j];
    }
  }

  Mp->eigsol.anies=i;
  Mp->eigsol.aerres=maxerror;
  
  for (i=0;i<nv;i++){
    w[i]-=shift;
    if (w[i]<0.0){
      //  this can happen for zero eigenvalues
      //  because of cancellation errors, w[i] can be e.g. -1.0e-8
      w[i]=fabs(w[i]);
    }
  }


  delete [] z;
  delete [] wo;
  delete [] p;
}
