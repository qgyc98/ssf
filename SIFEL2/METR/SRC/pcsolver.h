#ifndef PCSOLVER_H
#define PCSOLVER_H

#include <stdio.h>

/**
   partially coupled problem
   transport problem influences mechanical problem
*/

void solve_pcouplprob ();

void newton_raphson_parcoupl (long lcid);
void newton_raphson_parcoupl_comp (long lcid);
void newton_raphson_parcoupl_common_dt (long lcid);
void newton_raphson_parcoupl_lin (long lcid);
void newton_raphson_parcoupl_nonlin (long lcid);
void newton_raphson_parcoupl_nonlin_old (long lcid);
void newton_raphson_parcoupl_nonlin_new (long lcid);

void vector_assemb (double *c,double *m,double *t);
void vector_decomp (double *c,double *m,double *t);

#endif
