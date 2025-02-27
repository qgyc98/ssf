#include "slsolver.h"
#include "global.h"
#include "probdesc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"

/**
   function solves problems of linear stability
   
   JK, 12.1.2003
   modification 8. 5. 2018
*/
void solve_linear_stability (double *x,double *w)
{
  long i,n;
  double *rhs;
 
  //  number of mechanical degrees of freedom
  n=Ndofm;

  //  aray for the right hand sde vector
  rhs = Lsrs->give_rhs (0);
  
  //  assembling of the right hand side
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
  }
  
  //  computation of internal forces, mainly of the normal forces
  for (i=0;i<Lsrs->nlc;i++){
    compute_req_val (i); 
  }
  
  //  initial stiffness matrix assembling
  for (i=0;i<Lsrs->nlc;i++){
    initial_stiffness_matrix (i);
  }
  
  
  inverse_iteration2 (x,w);
  
  print_eigenvalues (w);
  print_eigenvectors ();
  
}




/**
   function computes the smallest eigenvalue and appropriate eigenvector
   of the generalized problem Ax=wBx
   the eigenvector of the smallest eigenvalue is stored in Lsrs in lhs

   @param w - the smallest eigenvalue
   
   3.11.2001
*/
void inverse_iteration2 (double *x,double *w)
{
  long i,ni,n;
  double nom,denom,zero,rho,prevrho,err,error;
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
  Ismat->gmxv (x,z);

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
    Ismat->gmxv (x,z);

    //  x.M.x
    denom=ss(x,z,n);
    rho=nom/denom;
    denom=1.0/sqrt(denom);
    //  normed vector z
    cmulv (denom,z,n);

    if (Mespr==1)  printf ("\n iterace %4ld     rho %15.12e",i,rho);
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
