#ifndef EDSOLVER_H
#define EDSOLVER_H

#include <stdio.h>

void solve_eigen_dynamics (double *x,double *w);

void inverse_iteration (double *x,double *w);
void subspace_iter_ortho (double *x,double *w);
void subspace_iter_jac (double *x,double *w);
void subspace_shift_iter_ortho (double *x,double *w);

#endif
