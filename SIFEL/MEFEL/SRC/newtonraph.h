#ifndef NEWTONRAPH_H
#define NEWTONRAPH_H



#include "nonlinman.h"



// The function solves the system of nonlinear equations by the Newton-Raphson method.
double gnewton_raphson (long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp,
                        long li, double ilambda, answertype outres);


double gnewton_raphson2 (long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp, 
                         double *flp, double *flc, long li, double ilambda, answertype outres);

double gnewton_raphson_one_step(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, 
                                double *dr, double *fi, double &dtr, long istep, long &j, long li, answertype fusm);

double compute_res_norm_nr(long lcid, resnormt normt, double *fb, double *fa, long n, double &norfa);
#endif
