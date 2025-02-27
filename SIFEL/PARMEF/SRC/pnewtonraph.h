#ifndef PNEWTONRAPH_H
#define PNEWTONRAPH_H



#include "nonlinman.h"



// The function solves the system of nonlinear equations by the Newton-Raphson method.
double par_gnewton_raphson_one_step(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, double *fi,
                                    double &dtr, long istep, long &j, long li, long ini, double ierr);

// function checks for the minimum time step reduction required by the material models
double par_check_dtr(void);

double par_compute_res_norm_nr(long lcid, resnormt normt, double *fb, double *fa, long n, double &norfa);

double par_compute_react_norm();
#endif
