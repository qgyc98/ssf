#ifndef NNPSOLVERT_H
#define NNPSOLVERT_H

#include <stdio.h>

void solve_nonlinear_nonstationary_problem ();

void nonlinear_nonstat_solv ();
void nonlinear_nonstat_solv_oldd ();

void nonlinear_nonstat_solv_vform ();
void nonlinear_nonstat_solv_dform ();
void nonlinear_nonstat_solv_pokus ();
void nonlinear_nonstat_solv_new ();
void nonlinear_nonstat_solv_linesearch ();

long norm_computation (double *f,double *rhs,double err,long nt,long sc);
long norm_computation_vec (double *f,double *rhs,double *err,double *thresh,long nt,long sc, double &norfb);

void nonlinear_nonstat_solv_dform_dneska ();

void nonlinear_nonstat_solv_fnr_dform ();
void nonlinear_nonstat_solv_nr_dform ();
void nonlinear_nonstat_solv_fnr_dform_old ();

#endif
