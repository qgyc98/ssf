#ifndef NPSOLVERT_H
#define NPSOLVERT_H

#include <stdio.h>
#include "timecontr.h"
#include "npglvec.h"

void solve_nonstationary_problem ();

void nonstat_solv_vform_comp ();
void nonstat_solver_init (long lcid, long &rest_calc, np_glob_vec &np_gv);
long one_step_linear (long lcid,double time, double dt, double prev_time, long rest_calc, long istep, long li, np_glob_vec &np_gv);
long one_step_nonlinear (long lcid,double time, double dt, double prev_time, long rest_calc, long istep, long li, np_glob_vec &np_gv);

/// function updates status of elements and nodes after restoreg from backup in growing mechanical problem
void update_elnod_stat_after_hdbrestt(np_glob_vec &np_gv);

/// updates status of nodes and elements according to previous time step in growing mechanical problem
void update_elemnod_prev_timet(double prev_time, long ncd, long nce);

/// updates status of nodes, DOFs and elements in growing mechanical problem
void update_elnod_statt(long lcid, long istep, double prev_time, double *lhs, double *tdlhs, long &tnce, long &tncd);

/// calculates and hadles force %vector due to removed elements
//long forces_due2removed_elem(long lcid, long nrel, long istep, double prev_time, double *fi, double *fp, 
//                             double *flp, double *fb, double *r);///tady??!!

void linear_nonstat_solv_vform ();
void linear_nonstat_solv_dform ();

void linear_nonstat_radiation_solv_dform ();

void aux_nonlintime_print (FILE *aux, double *r, double l);

void linear_nonstat_solv_dform_subcycl ();


void explicit_difference_method (long lcid);
void optim_linear_nonstat_solv_dform ();

#endif
