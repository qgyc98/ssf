#ifndef PMTSOLVER_H
#define PMTSOLVER_H

#include "mtglvec.h"
#include <stdio.h>

void par_solve_timemech_prob ();

void par_solve_timemech_prob2 (long lcid);
void par_visco_solver_init(long lcid, long &rest_calc, mt_glob_vec &mt_gv);
long par_one_step (long lcid,double time, double dt, double &dtr, double prev_time, long rest_calc, long istep, long li, mt_glob_vec &mt_gv);
//long par_one_step (long lcid,double time, double dt, long istep, long li, mt_glob_vec &mt_gv);
void par_solve_timemech_prob (long lcid);

#endif
