#ifndef DNNPSOLVERT_H
#define DNNPSOLVERT_H

#include <stdio.h>

void solve_discont_nonlin_nonstationary_problem ();

/// d-form solver of nonlinear nonstationar problems
void nonlin_nonstat_dform ();

/// initialization phase of d-form solver for one_step concept
void nonstat_solver_dform_init (long lcid, np_glob_vec &np_gv);

/// function performs one time step of d-form solver for nonlinear non-stationary problem
long one_step_nonlinear_dform (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv);

/// d-form solver of nonlinear nonstationar problems built on one step concept
void nonstat_solv_dform_comp ();
#endif
