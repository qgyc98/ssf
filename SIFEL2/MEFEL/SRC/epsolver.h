#ifndef EPSOLVER_H
#define EPSOLVER_H

#include <stdio.h>

void solve_epressure ();  ///< main computation function for earth pressure problem type

void epressure_solver (long lcid); ///< function controling earth pressure solvers

void general_epressure (long lcid); ///< earth pressure solver with load iteration

void general_epressure_varsup (long lcid); ///< earth pressure solver with temporary supports

void femplast_epressure (long lcid); ///< earth pressure solver with temporary supports for plasticity models

void restore_displ(long lcid, double **bckr); ///< extraction of the reached displacement form the backup vector

double arclengthrv (long lcid, double rlambda, double rlerr); ///< the arc-length method which stops iteration on the given value of the lambda parameter

double newton_raphsonep (long lcid, double *forig); ///< the Newton-Raphsonmethop modified for starting from a reached state

#endif
