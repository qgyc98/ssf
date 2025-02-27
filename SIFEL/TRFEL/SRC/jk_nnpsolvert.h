#ifndef NNPSOLVERT_H
#define NNPSOLVERT_H

#include <stdio.h>

void solve_nonlinear_nonstationary_problem ();

void nonlinear_nonstat_solv ();
void nonlinear_nonstat_solv_old ();
void nonlinear_nonstat_solv_pokus ();
void nonlinear_nonstat_solv_new ();
void nonlinear_nonstat_solv_linesearch ();

long selectivenorm (double *f,double err,long ni,long nii);

#endif
