#include "homogtrans.h"
#include "spsolvert.h"
#include "nspsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include <string.h>
#include "npsolvert.h"


/**
   function solves homogeniztion for transport problems on PUC

   TKr, 17/08/2010
*/
void transport_homogenization_old_old ()
{
  long i,j;
  long ntm,ncomp;
  double puc_area;
  double *rhs,*lhs;
  matrix ltm,lm,kpuc,km,kM,kinv,khelp,cM;

  //  initiation of transport material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //solving of fluctuation patr
  switch (Tp->tprob){
  case stationary_problem:{
    //stationary problem solution
    solve_stationary_problem ();
    break;
  }
  case nonlinear_stationary_problem:{
    solve_nonlinear_stationary_problem_pokus ();
    //solve_nonlinear_stationary_problem ();
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

  //---------------------------------------//
  //updating of unknown from master unknown, gradients and fluctuation part
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //homogenization
  //homogenization of conductivity matrix
  //computing of l matrix and l_t matrix
  ntm = Tp->ntm;
  ncomp = Tt->nodes[0].ncompgrad;
  allocm (Ndoft,ntm*ncomp,ltm);
  allocm (ntm*ncomp,Ndoft,lm);
  allocm (Ndoft,Ndoft,kpuc);
  allocm (Ndoft,Ndoft,kinv);
  allocm (Ndoft,ntm*ncomp,khelp);

  assemble_l_matrix (lm,ltm);

  
  if (Mesprt != 0){
    fprintf(Outt,"\n\n matice Lm\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<Ndoft;j++)
	fprintf(Outt,"%e  ",lm[i][j]);
      fprintf(Outt,"\n");
    }

    fprintf(Outt,"\n\n matice Ltm\n");
    for (i=0;i<Ndoft;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",ltm[i][j]);
      fprintf(Outt,"\n");
    }
  }
  

  //conductivity matrix
  //assemble_conductivity_matrix (kpuc);//conductivity matrix of PUC stored as an array
  //conductivity matrix of PUC
  conductivity_matrix (0);

  //computing of average matrix (D_m in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,km);

  assemble_average_d_matrix(km,puc_area);
  
  //computing of macroscopic matrix (D_M in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,kM);  

  //contribution from fluctuation part
  //inversion of matrix K
  // lhs
  lhs = new double [kinv.m];
  nullv (lhs,kinv.m);
  // rhs
  rhs = new double [kinv.m];
  nullv (rhs,kinv.m);
  
  for(i=0;i<kinv.m;i++){
    nullv (rhs,kinv.m);
    nullv (lhs,kinv.m);
    rhs[i] = 1.0;
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    for(j=0;j<kinv.n;j++){
      kinv[j][i]=lhs[j];
    }
  }

  //inversion of matrix kpuc stored as an array
  //invm(kpuc,kinv,Tp->zero);

  /* //kontrola:
     matrix aux(kinv.m, kinv.n);
     mxm(kpuc, kinv, aux);
     printm(aux, Outt, 3, 11);
     destrm(aux);
  */

  mxm(kinv,ltm,khelp);

  //pomocny tisk:
  fprintf(Outt,"\n\n matice khelp\n");
  for (j=0;j<ntm*ncomp;j++){
    for (i=0;i<Ndoft;i++){
      fprintf(Outt,"%e  ",khelp[i][j]);
    }
    fprintf(Outt,"\n");
  }
  
  mxm(lm,khelp,kM);
  cmulm(-1.0,kM);
  puc_area = 1.0/puc_area;
  cmulm(puc_area,kM);


  if (Mesprt != 0){
    fprintf(Outt,"\n\n matice D^m\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",km[i][j]);
      fprintf(Outt,"\n");
    }
    
    fprintf(Outt,"\n\n fluktuacni prispevky D\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",kM[i][j]);
      fprintf(Outt,"\n");
    }
  }

  //final macroscopic conductivity matrix (D_M in TRFEL notation)
  addm(km,kM,kM);

  //computing of macroscopic cpacity matrix (C_M in TRFEL notation)
  allocm (ntm,ntm,cM);
  //homogenization of capacity matrix
  //computing of average matrix
  assemble_average_c_matrix(cM);

  if (Mesprt != 0){
    fprintf(Outt,"\n\n vysledna matice D^M\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",kM[i][j]);
      fprintf(Outt,"\n");
    }
    
    fprintf(Outt,"\n\n vysledna matice C^M\n");
    for (i=0;i<ntm;i++){
      for (j=0;j<ntm;j++)
	fprintf(Outt,"%e  \n",cM[i][j]);
      fprintf(Outt,"\n");
    }
  }
  
  destrm (lm);
  destrm (ltm);
  destrm (kpuc);
  destrm (kinv);
  destrm (km);
  destrm (kM);
  destrm (cM);
  destrm (khelp);
  delete [] rhs;
  delete [] lhs;
}



/**
   function solves homogenization for transport problems on PUC

   TKr+JK, 20/04/2015
*/
void transport_homogenization ()
{
  long i,j;
  long ntm,ncomp;
  double puc_area;
  matrix lhs,lm,ltm,kpuc,km,kM,cM;
  double *vhelp;

  //  initiation of transport material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation_puc ();


  //---------------------------------------//
  //solving of fluctuation patr
  switch (Tp->tprob){
  case stationary_problem:{
    //stationary problem solution
    solve_stationary_problem ();
    break;
  }
  case nonlinear_stationary_problem:{
    solve_nonlinear_stationary_problem_pokus ();
    //solve_nonlinear_stationary_problem ();
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

  //---------------------------------------//
  //updating of unknown from master unknown, gradients and fluctuation part
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //homogenization
  //homogenization of conductivity matrix
  //computing of l and lt matrix
  ntm = Tp->ntm;
  ncomp = Tt->nodes[0].ncompgrad;
  allocm (ntm*ncomp,Ndoft,lm);
  allocm (Ndoft,ntm*ncomp,ltm);
  allocm (Ndoft,Ndoft,kpuc);

  assemble_l_matrix (lm);
  assemble_lt_matrix (ltm);
  
  if (Mesprt != 0){
    fprintf(Outt,"\n\n matice Lm\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<Ndoft;j++)
	fprintf(Outt,"%e  ",lm[i][j]);
      fprintf(Outt,"\n");
    }
    fprintf(Outt,"\n\n matice Ltm\n");
    for (i=0;i<Ndoft;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  ",ltm[i][j]);
      fprintf(Outt,"\n");
    }
  }

  //conductivity matrix of PUC
  conductivity_matrix (0);

  //computing of average matrix (D_m in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,km);

  assemble_average_d_matrix(km,puc_area);
  
  //computing of macroscopic matrix (D_M in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,kM);  

  //contribution from fluctuation part
  //inversion of matrix K
  allocm (ntm*ncomp,Ndoft,lhs);
  // vhelp
  vhelp = new double [Ndoft];

  //inversion of matrix kpuc stored as an array
  for (i=0;i<ntm*ncomp;i++){
    if (Mesprt != 0) 
      fprintf(Outt,"\n\n tisk vhelp pro i=%ld\n",i);
    for (j=0;j<Ndoft;j++){
      vhelp[j] = ltm[j][i];//correct indexes 
      if (Mesprt != 0) 
	fprintf(Outt,"%e  \n",vhelp[j]);
    }
    Tp->ssle->solve_system (Gtt,Kmat,lhs.a+i*Ndoft,vhelp,Outt);    
  }

  if (Mesprt != 0){
    fprintf(Outt,"\n\n matice lhs\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<Ndoft;j++)
	fprintf(Outt,"%e  ",lhs[i][j]);
      fprintf(Outt,"\n");
    }
  }
  
  mxmt (lm,lhs,kM);
  cmulm(-1.0,kM);
  puc_area = 1.0/puc_area;
  cmulm(puc_area,kM);
  
  if (Mesprt != 0){
    fprintf(Outt,"\nkontrolni tisk pro novou verzi homogenizace:\n");
    fprintf(Outt,"\n\n matice D^m\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",km[i][j]);
      fprintf(Outt,"\n");
    }
    
    fprintf(Outt,"\n\n fluktuacni prispevky D\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",kM[i][j]);
      fprintf(Outt,"\n");
    }
  }

  //final macroscopic conductivity matrix (D_M in TRFEL notation)
  addm (kM, km, kM);

  //computing of macroscopic cpacity matrix (C_M in TRFEL notation)
  allocm (ntm,ntm,cM);
  //homogenization of capacity matrix
  //computing of average matrix
  assemble_average_c_matrix(cM);
  
  
  if (Mesprt != 0){
    fprintf(Outt,"\n\n vysledna matice D^M\n");
    for (i=0;i<ntm*ncomp;i++){
      for (j=0;j<ntm*ncomp;j++)
	fprintf(Outt,"%e  \n",kM[i][j]);
      fprintf(Outt,"\n");
    }
    
    fprintf(Outt,"\n\n vysledna matice C^M\n");
    for (i=0;i<ntm;i++){
      for (j=0;j<ntm;j++)
	fprintf(Outt,"%e  \n",cM[i][j]);
      fprintf(Outt,"\n");
    }
  }
  
  destrm (lm);
  destrm (ltm);
  destrm (lhs);
  destrm (kpuc);
  destrm (km);
  destrm (kM);
  destrm (cM);
  delete [] vhelp;
}


/**
   function solves homogeniztion for transport problems on PUC and returns overall (homogenized) conductivity and capacity matrices 

   @param unkn_r - unknows and gradients of unknowns in integration point (in master node on PUC) stored as an array (one dimensional)
   @param matrix_k - homogenized conductivity matrix of PUC stored as an array (one dimensional)
   @param matrix_c - homogenized capacity matrix of PUC stored as an array (one dimensional)

   TKr, 17/08/2010
*/
void paral_transport_homogenization_old_old (double *unkn_r,double *matrix_k,double *matrix_c)
{
  long i,j,k;
  long ntm,ncomp;
  double puc_area;
  double *rhs,*lhs;
  matrix ltm,lm,kpuc,km,kM,kinv,khelp,cM;

  //  initiation of transport material models
  Tm->initmaterialmodels();

  //initiation of master node unknowns and gradients
  ntm = Tp->ntm;
  ncomp = Tt->nodes[0].ncompgrad;
  k = -1;
  for (i=0;i<ntm;i++){
    k++;
    Tb->lc[i].masterval = unkn_r[k];
    for (j=0;j<ncomp;j++){
      k++;
      Tb->lc[i].mastergrad[j] = unkn_r[k];
    }
  }
  
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //solving of fluctuation patr
  switch (Tp->tprob){
  case stationary_problem:{
    //stationary problem solution
    solve_stationary_problem ();
    break;
  }
  case nonlinear_stationary_problem:{
    solve_nonlinear_stationary_problem_pokus ();
    //solve_nonlinear_stationary_problem ();
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

  //---------------------------------------//
  //updating of unknown from master unknown, gradients and fluctuation part
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //homogenization
  //homogenization of conductivity matrix
  //computing of l matrix and l_t matrix
  allocm (Ndoft,ntm*ncomp,ltm);
  allocm (ntm*ncomp,Ndoft,lm);
  allocm (Ndoft,Ndoft,kpuc);
  allocm (Ndoft,Ndoft,kinv);
  allocm (Ndoft,ntm*ncomp,khelp);

  assemble_l_matrix (lm,ltm);

  //conductivity matrix
  //conductivity matrix of PUC
  conductivity_matrix (0);

  //computing of average matrix (D_m in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,km);

  assemble_average_d_matrix(km,puc_area);
  
  //computing of macroscopic matrix (D_M in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,kM);  

  //contribution from fluctuation part
  //inversion of matrix K
  // lhs
  lhs = new double [kinv.m];
  nullv (lhs,kinv.m);
  // rhs
  rhs = new double [kinv.m];
  nullv (rhs,kinv.m);
  
  for(i=0;i<kinv.m;i++){
    nullv (rhs,kinv.m);
    nullv (lhs,kinv.m);
    rhs[i] = 1.0;
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    for(j=0;j<kinv.n;j++){
      kinv[j][i]=lhs[j];
    }
  }

  mxm(kinv,ltm,khelp);
  mxm(lm,khelp,kM);
  cmulm(-1.0,kM);
  puc_area = 1.0/puc_area;
  cmulm(puc_area,kM);

  //final macroscopic conductivity matrix (D_M in TRFEL notation)
  addm(km,kM,kM);

  //computing of macroscopic cpacity matrix (C_M in TRFEL notation)
  allocm (ntm,ntm,cM);
  //homogenization of capacity matrix
  //computing of average matrix
  assemble_average_c_matrix(cM);


  //storing of overall properties
  //conductivity matrix
  k = 0;
  for (i=0;i<ntm*ncomp;i++){
    for (j=0;j<ntm*ncomp;j++){
      matrix_k[k] = kM[i][j];
      k++;
    }
  }

  //capacity matrix
  k = 0;
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      matrix_c[k] = cM[i][j];
      k++;
    }
  }
  
  destrm (lm);
  destrm (ltm);
  destrm (kpuc);
  destrm (kinv);
  destrm (km);
  destrm (kM);
  destrm (cM);
  destrm (khelp);
  delete [] rhs;
  delete [] lhs;
}


/**
   function solves homogenization for transport problems on PUC and returns overall (homogenized) conductivity and capacity matrices 

   @param unkn_r - unknows and gradients of unknowns in integration point (in master node on PUC) stored as an array (one dimensional)
   @param matrix_k - homogenized conductivity matrix of PUC stored as an array (one dimensional)
   @param matrix_c - homogenized capacity matrix of PUC stored as an array (one dimensional)

   TKr+JK, 20/04/2015
*/
void paral_transport_homogenization (double *unkn_r,double *matrix_k,double *matrix_c)
{
  long i,j,k;
  long ntm,ncomp;
  double puc_area;
  matrix lhs,lm,ltm,kpuc,km,kM,cM;
  double *vhelp;

  //  initiation of transport material models
  Tm->initmaterialmodels();

  //initiation of master node unknowns and gradients
  ntm = Tp->ntm;
  ncomp = Tt->nodes[0].ncompgrad;
  k = -1;
  for (i=0;i<ntm;i++){
    k++;
    Tb->lc[i].masterval = unkn_r[k];
    for (j=0;j<ncomp;j++){
      k++;
      Tb->lc[i].mastergrad[j] = unkn_r[k];
    }
  }
  
  //kontrolni tisk gradientu a neznamych:
  /* for (i=0;i<ntm;i++){
     fprintf(Outt,"\n\n neznama Tb->lc[%ld].masterval = %e\n",i,Tb->lc[i].masterval);
     for (j=0;j<ncomp;j++){
     fprintf(Outt,"\n\n gradient Tb->lc[%ld].mastergrad[%ld] = %e\n",i,j,Tb->lc[i].mastergrad[j]);
     }
     }
  */

  //  nodes - integration points interpolation
  approximation_puc ();
 
  //---------------------------------------//
  //solving of fluctuation part
  switch (Tp->tprob){
  case stationary_problem:{
    //stationary problem solution
    solve_stationary_problem ();
    break;
  }
  case nonlinear_stationary_problem:{
    solve_nonlinear_stationary_problem_pokus ();
    //solve_nonlinear_stationary_problem ();
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

  //---------------------------------------//
  //updating of unknown from master unknown, gradients and fluctuation part
  //  nodes - integration points interpolation
  approximation_puc ();

  //---------------------------------------//
  //homogenization
  //homogenization of conductivity matrix
  //computing of l matrix and l_t matrix
  reallocm (ntm*ncomp,Ndoft,lm);
  allocm (Ndoft,ntm*ncomp,ltm);
  reallocm (Ndoft,Ndoft,kpuc);

  assemble_l_matrix (lm);
  assemble_lt_matrix (ltm);
  
  //conductivity matrix
  //conductivity matrix of PUC
  conductivity_matrix (0);
  
  //computing of average matrix (D_m in TRFEL notation)
  reallocm (ntm*ncomp,ntm*ncomp,km);

  assemble_average_d_matrix(km,puc_area);
  
  //computing of macroscopic matrix (D_M in TRFEL notation)
  allocm (ntm*ncomp,ntm*ncomp,kM);  

  //contribution from fluctuation part
  //inversion of matrix K
  reallocm (ntm*ncomp,Ndoft,lhs);
  // vhelp
  vhelp = new double [Ndoft];

  //inversion of matrix kpuc stored as an array
  for (i=0;i<ntm*ncomp;i++){
    for (j=0;j<Ndoft;j++)
      vhelp[j] = ltm[j][i]; //correct indexes
    Tp->ssle->solve_system (Gtt,Kmat,lhs.a+i*Ndoft,vhelp,Outt);
  }
  
  mxmt (lm,lhs,kM);
  cmulm(-1.0,kM);
  puc_area = 1.0/puc_area;
  cmulm(puc_area,kM);
  
  // this printing below is only for debug 
  /* fprintf(Outt,"\n\n matice D^m\n");
     for (i=0;i<ntm*ncomp;i++){
     for (j=0;j<ntm*ncomp;j++)
     fprintf(Outt,"%e  \n",km[i][j]);
     fprintf(Outt,"\n");
     }
     
     fprintf(Outt,"\n\n fluktuacni prispevky D\n");
     for (i=0;i<ntm*ncomp;i++){
     for (j=0;j<ntm*ncomp;j++)
     fprintf(Outt,"%e  \n",kM[i][j]);
     fprintf(Outt,"\n");
     }
  */
  
  //final macroscopic conductivity matrix (D_M in TRFEL notation)
  addm (kM, km, kM);

  //computing of macroscopic cpacity matrix (C_M in TRFEL notation)
  reallocm (ntm,ntm,cM);
  //homogenization of capacity matrix
  //computing of average matrix
  assemble_average_c_matrix(cM);

  //storing of overall properties
  //conductivity matrix
  k = 0;
  for (i=0;i<ntm*ncomp;i++){
    for (j=0;j<ntm*ncomp;j++){
      matrix_k[k] = kM[i][j];
      k++;
    }
  }

  //capacity matrix
  k = 0;
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      matrix_c[k] = cM[i][j];
      k++;
    }
  }

  // this printing below is only for debug 
  /* fprintf(Outt,"\n\n vysledna matice D^M\n");
     for (i=0;i<ntm*ncomp;i++){
     for (j=0;j<ntm*ncomp;j++)
     fprintf(Outt,"%e  \n",kM[i][j]);
     fprintf(Outt,"\n");
     }
     
     fprintf(Outt,"\n\n vysledna matice C^M\n");
     for (i=0;i<ntm;i++){
     for (j=0;j<ntm;j++)
     fprintf(Outt,"%e  \n",cM[i][j]);
     fprintf(Outt,"\n");
     }
     fflush(Outt);
  */

  if (Kmat != NULL){
    //if (Mesprt != 0)
    //fprintf (stdout,"\n\n Cleaning of conductivity matrix \n\n ");
    delete Kmat;
    Kmat=NULL;
  }

  destrm (lm);
  destrm (ltm);
  destrm (lhs);
  destrm (kpuc);
  destrm (km);
  destrm (kM);
  destrm (cM);
  delete [] vhelp;
}
