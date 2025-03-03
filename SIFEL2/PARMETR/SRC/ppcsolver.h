#ifndef PPCSOLVER_H
#define PPCSOLVER_H

#include <stdio.h>
struct np_glob_vec;
struct mt_glob_vec;
struct nonlinman;

void par_solve_pcouplprob ();

void par_newton_raphson_parcoupl_comp (long lcid);
void par_newton_raphson_parcoupl_common_dt (long lcid);
void par_nonstat_trfel_init (long lcid, np_glob_vec &np_gv);
void par_visco_mefel_init(long lcid, mt_glob_vec &mt_gv);
long par_one_step_trfel_linear (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv);
long par_one_step_trfel_nonlinear (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv);
long par_one_step_mefel (long lcid,double time, double dt, long istep, long li, mt_glob_vec &mt_gv);
double par_gnewton_raphson_one_step_mefel(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, 
					  double *fi,long istep, long &j, long li, long ini, double ierr);

void par_newton_raphson_parcoupl (long lcid);
void par_newton_raphson_parcoupl_lin (long lcid);
void par_newton_raphson_parcoupl_nonlin (long lcid);


#endif
