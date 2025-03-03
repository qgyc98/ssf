#include "eigvalsol.h"
#include "iotools.h"
#include "global.h"
#include <string.h>

eigvalsol::eigvalsol (void)
{
  //  type of solver of eigenvalues and eigenvectors
  teigsol = no_eigsolver;
  //  the number of required eigenvectors
  neigv=0;
  //  the number of vectors used in computation
  nev=0;
  //  the maximum number of iterations
  nies=0;
  //  the number of performed iterations
  anies=0;
  //  required error
  erres=0.0;
  //  attained error
  aerres=0.0;
  //  the maximum number of iteration in Jacobi' method
  nijmr=0;
  //  the number of thresholds in Jacobi' method
  njacthr=0;
  //  array containing thresholds in Jacobi' method
  jacthr=NULL;
  //  shift
  shift=0.0;
}

eigvalsol::~eigvalsol (void)
{
}

/**
   function reads data from input file
   
   @param in - pointer to the input file
   @param mespr - message printing indicator
   
   JK, 20.8.2005
*/
void eigvalsol::read (XFILE *in)
{
  long i;

  //  type of eigenvalue problem solver
  //  0 - no_eigsolver
  //  1 - inv_iteration
  //  4 - subspace_it_jacobi
  //  5 - subspace_it_gsortho
  //  6 - shifted_subspace_it_gsortho
  xfscanf (in,"%k%m","type_of_eig_solver",&eigensolver_kwdset,(int*)&teigsol);
  
  switch (teigsol){
  case inv_iteration:{
    if (Mespr==1)  fprintf (stdout,"\n eigenvalue problem will be solved by inverse iteration method");
    
    //  number of eigenvectors
    nev=1;  neigv=1;
    
    xfscanf (in,"%ld %lf",&nies,&erres);
    break;
  }
  case subspace_it_jacobi:{
    if (Mespr==1){
      fprintf (stdout,"\n eigenvalue problem will be solved by subspace");
      fprintf (stdout,"\n iteration method with Jacobi method of rotations");
    }
    
    //  neigv - required number of eigenvectors
    //  nev - the number of vectors which will be used in the iteration
    //  nies - the maximum number of iterations
    //  erres - the required norm of residuals
    xfscanf (in,"%ld %ld %ld %lf",&neigv,&nev,&nies,&erres);
    if (neigv>nev){
      print_err("wrong number of vectors is used in subspace iteration method",__FILE__,__LINE__,__func__);
    }
    xfscanf (in,"%ld %ld",&nijmr,&njacthr);
    jacthr = new double [njacthr];
    for (i=0;i<njacthr;i++){
      xfscanf (in,"%lf",jacthr+i);
    }
    
    break;
  }
  case subspace_it_gsortho:{
    if (Mespr==1){
      fprintf (stdout,"\n eigenvalue problem will be solved by subspace");
      fprintf (stdout,"\n iteration method with Gram-Schmidt orthonormalization technique");
    }
    
    //  neigv - required number of eigenvectors
    //  nev - the number of vectors which will be used in the iteration
    //  nies - the maximum number of iterations
    //  erres - the required norm of residuals
    xfscanf (in,"%ld %ld %ld %le",&neigv,&nev,&nies,&erres);
    if (neigv>nev){
      print_err("wrong number of vectors is used in subspace iteration method",__FILE__,__LINE__,__func__);
      print_err("nev must be greater or equal to neigv",__FILE__,__LINE__,__func__);
    }
    
    break;
  }
  case shifted_subspace_it_gsortho:{
    if (Mespr==1){
      fprintf (stdout,"\n eigenvalue problem will be solved by subspace");
      fprintf (stdout,"\n iteration method with Gram-Schmidt orthonormalization technique");
    }
    
    //  neigv - required number of eigenvectors
    //  nev - the number of vectors which will be used in the iteration
    //  nies - the maximum number of iterations
    //  erres - the required norm of residuals
    //  shift
    xfscanf (in,"%ld %ld %ld %le %le",&neigv,&nev,&nies,&erres,&shift);
    if (neigv>nev){
      print_err("wrong number of vectors is used in subspace iteration method",__FILE__,__LINE__,__func__);
      print_err("nev must be greater or equal to neigv",__FILE__,__LINE__,__func__);
    }    
    break;
  }
  default:{
    print_err("wrong type of eigenvalue problem solver is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function prints data about solver of eigenvalues and eigenvectors
   
   @param out - pointer to output file
   
   JK, 20.8.2005
*/
void eigvalsol::print (FILE *out)
{
  long i;
  
  fprintf (out,"%d\n",(int)teigsol);
  
  switch (teigsol){
  case inv_iteration:{
    fprintf (out, "%ld %e", nies, erres);
    break;
  }
  case subspace_it_jacobi:{
    fprintf (out, "%ld %ld %ld %e", neigv, nev, nies, erres);
    fprintf (out," %ld %ld\n",nijmr,njacthr);
    for (i=0;i<njacthr;i++){
      fprintf (out," %e",jacthr[i]);
    }
    break;
  }
  case subspace_it_gsortho:{
    fprintf (out, "%ld %ld %ld %e", neigv, nev, nies, erres);
    break;
  }
  case shifted_subspace_it_gsortho:{
    fprintf (out, "%ld %ld %ld %e %e", neigv, nev, nies, erres, shift);
    break;
  }  default:{
    print_err("wrong type of eigenvalue problem solver is required",__FILE__,__LINE__,__func__);
  }
  }
  
  fprintf (out,"\n");
  
}

