#ifndef PNPSOLVERT_H
#define PNPSOLVERT_H

#include <stdio.h>
#include "timecontr.h"
#include "npglvec.h"

void par_solve_nonstationary_problem ();
void par_linear_nonstat_solv_vform ();
void par_nonstat_solv_vform_comp ();
void par_nonstat_solver_init (long lcid, long &rest_calc, np_glob_vec &np_gv);
long par_one_step_linear (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv);
long par_one_step_nonlinear (long lcid,double time, double dt, double prev_time, long rest_calc, long istep, long li, np_glob_vec &np_gv, double *gv, double *grhs);
long par_norm_computation_vec (double *f,double *rhs,double *err,double *thresh, long nt,long sc, double *gf, double *grhs, double &norfb);

#endif
