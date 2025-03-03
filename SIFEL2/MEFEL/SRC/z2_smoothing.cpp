#include "z2_smoothing.h"
#include "gadaptivity.h"
#include "adaptivity.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "alias.h"
#include "elemhead.h"
#include "vector.h"
#include "matrix.h"
#include "skyline.h"




z2_smoothing::z2_smoothing (long ncomp)
{
  z2_smoothing::ncomp = ncomp;
  nn = Mt->nn;
  ne = Mt->ne;
  flags = Ada->give_adaptflag();
}

z2_smoothing::~z2_smoothing ()  { }



void z2_smoothing::run (vector *rsigfull)
{
  fprintf(stdout,"\n ==============  Z2_smoothing  ===============");

  vector ainv,ntdbr;
  skyline *ntn_sky;
  char label[50];

  allocv (ncomp*nn,ntdbr);
  compute_ntdbr (ntdbr);

  if (1) { //!Mp->diagonalization){
    ntn_sky = new skyline;
    
    alloc_ntn (ntn_sky);
    compute_ntn_sky (ntn_sky);
    compute_rsigfull (ntn_sky,ntdbr,rsigfull);
    
    delete ntn_sky;
  }
  else {
    fprintf (stdout,"\n used diagonalization");
    allocv (nn,ainv);
    compute_ainv (ainv);
    compute_rsigfull (ainv,ntdbr,rsigfull);
  }


  // *** PRINT ***

  sprintf(label, "ntdbr:");
  fprintf_vector(Out,ntdbr,label,9);
  //delete [] label;
  //fprintf_long_1D_1(stdout,ntn_sky->adr,ntn_sky->n+1,"ntn_sky.adr:");
  //fprintf_double_1D_1(stdout,ntn_sky->a,ntn_sky->negm,"ntn_sky:");
  //fprintf_vector(stdout,ainv,"ainv:");
}


void z2_smoothing::compute_ntdbr (vector &ntdbr)
{
  long i,j,te,nne;
  ivector cn;
  vector ntdbri;

  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    nne = Mt->give_nne (i);
    allocv (ncomp*nne,ntdbri);

    switch (te){
    case planeelementlt:{     Pelt->ntdbr_vector (i,ntdbri);   break; }
    case planeelementlq:{     Pelq->ntdbr_vector (i,ntdbri);   break; }
    case planeelementqt:{     Peqt->ntdbr_vector (i,ntdbri);   break; }
    case planeelementqq:{     Peqq->ntdbr_vector (i,ntdbri);   break; }
    case lineartet:{          Ltet->ntdbr_vector (i,ntdbri);   break; }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr," adaptivity::compute_ntdbr (%s, line %d)",__FILE__,__LINE__);
    }
    }

    allocv (nne,cn);
    give_adapt_code_numbers (i,cn);
    for (j=0;j<ncomp;j++){
      locglob (ntdbr.a + j*nn,ntdbri.a + j*nne,cn.a,nne);
    }
    destrv (ntdbri);  destrv (cn);
  }
}

void z2_smoothing::alloc_ntn (skyline *ntn_sky)
{
  ntn_sky->allocadr (nn);

  column_lengths_nn (ntn_sky);

  ntn_sky->addresses ();
  ntn_sky->neglobmat ();
  ntn_sky->allocglomat ();

  if (flags & 1)
    fprintf(stdout,"\n number of skyline entries negm %ld",ntn_sky->negm);
}


void z2_smoothing::column_lengths_nn (skyline *ntn_sky)
{
  long i,nne;
  ivector cn;

  for (i=0;i<ne;i++){
    nne = Mt->give_nne (i);
    allocv (nne,cn);
    give_adapt_code_numbers (i,cn);
    ntn_sky->column_lengths_elem (cn.a,nne);
    destrv (cn);
  }
}

void z2_smoothing::give_adapt_code_numbers (long eid,ivector &cn)
{                                                              // cn E <1;ncomp*nn>
  long i,nne;
  
  nne = cn.n;
  Mt->give_elemnodes (eid,cn);
  
  for (i=0;i<nne;i++)
    cn[i] += 1;
}

void z2_smoothing::compute_ntn_sky (skyline *ntn_sky)
{
  long i,te,nne;
  ivector cn;
  matrix ntni;
  
  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    nne = Mt->give_nne (i);
    allocm (nne,nne,ntni);
    
    switch (te){
    case planeelementlt:{    Pelt->ntn_matrix (i,ntni);   break; }
    case planeelementlq:{    Pelq->ntn_matrix (i,ntni);   break; }
    case planeelementqt:{    Peqt->ntn_matrix (i,ntni);   break; }
    case planeelementqq:{    Peqq->ntn_matrix (i,ntni);   break; }
    case lineartet :{        Ltet->ntn_matrix (i,ntni);   break; }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr," adaptivity::compute_ntn_sky (%s, line %d)",__FILE__,__LINE__);
    }
    }
    
    allocv (nne,cn);
    give_adapt_code_numbers (i,cn);
    ntn_sky->localized (ntni.a,cn.a,nne);
    destrm (ntni);  destrv(cn);
  }
}

void z2_smoothing :: compute_rsigfull (skyline *ntn_sky,vector &ntdbr, vector *rsigfull)
{
  long i, j;
  vector sol(nn);
  
  ntn_sky->ldl_sky (sol.a,ntdbr.a,Mp->zero,2);
  
  for (i=0;i<ncomp;i++) {
    ntn_sky->ldl_sky (sol.a, ntdbr.a + i*nn, Mp->zero, 3);
    
    for (j=0; j<nn; j++)
      rsigfull[j].a[i] = sol.a[j];
  }
}

void z2_smoothing::compute_ainv (vector &ainv)
{
  long i,te,nne;
  ivector cn;
  vector ntnvi,ntn(nn);
  matrix ntnmi;

  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    nne = Mt->give_nne (i);
    allocv (nne,ntnvi);
    allocm (nne,nne,ntnmi);

    switch (te){
    case planeelementlt:{      Pelt->ntn_matrix	(i,ntnmi);   break; }
    case planeelementlq:{      Pelq->ntn_matrix	(i,ntnmi);   break; }
    case planeelementqt:{      Peqt->ntn_matrix	(i,ntnmi);   break; }
    case planeelementqq:{      Peqq->ntn_matrix	(i,ntnmi);   break; }
    case lineartet:{           Ltet->ntn_matrix	(i,ntnmi);   break; }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr," adaptivity::compute_ainv (%s, line %d)",__FILE__,__LINE__);
    }
    }

    ntnmtov (ntnmi,ntnvi);

    allocv (nne,cn);
    give_adapt_code_numbers (i,cn);
    locglob (ntn.a,ntnvi.a,cn.a,nne);
    destrv (ntnvi); destrm (ntnmi);  destrv (cn);
  }

  for (i=0;i<nn;i++)
    ainv[i] = 1/ntn[i];
}


void z2_smoothing::compute_rsigfull (vector &ainv,vector &ntdbr, vector *rsigfull)
{
  long i,j;
  
  for (i=0; i<ncomp; i++)
    for (j=0; j<nn; j++)
      //rsigfull[j + i*nn] = ainv[j] * ntdbr[j + i*nn];
      rsigfull[j].a[i] = ainv[j] * ntdbr[j + i*nn];
}
