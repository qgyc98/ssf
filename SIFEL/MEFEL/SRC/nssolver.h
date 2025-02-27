#ifndef NSSOLVER_H
#define NSSOLVER_H

#include <stdio.h>
#include "galias.h"
struct matrix;
struct vector;
class nonlinman;

/// function calls appropriate solver for the nonlinear statics problem
void solve_nonlinear_statics ();

/// function calls appropriate solver of system of nonlinear equations for the given load case
void nonlinear_solver (long lcid, double *ra, double *fa, double *fc, double *fp, double *flc, double *flp, 
                       double ilambda, long li);

/// The function assembles the attained load %vector fa.
void assemble_attained_load_vector(double *fa, double *fc, double *fp, long n, double lambda, nonlinman *nlman);

/// The function performs the divergency check of the inner loop if it is required in the nlman setup.
long check_divergency(nonlinman *nlman, matrix &lsm_a, vector &lsm_r, vector &lsm_l, long j, double norf);

/*
/// Newton-Raphson method for systems of nonlinear equations
void newton_raphson (long lcid);

/// Newton-Raphson method for systems of nonlinear equations with direct displacement control
void displ_control (long lcid);

/// Newton-Raphson method for systems of nonlinear equations with direct displacement control (required lambda can be prescribed)
double displ_controlrv (long lcid, double inilambda, double dmplambda, double rlambda, double errl, long incrdispl);

/// Arclength method for systems of nonlinear equations
long arclength (long lcid);

long arclength_old (long lcid,long adaptcontrol);

/// function initializes selected dofs/nodes/distancies for the arclength method 
void seldofinit ();
*/

/// function assembles stiffness matrix if needed
long assemble_stiffness_matrix (long lcid,long i,long j,long li,answertype fusm);

/*
///  function determines increment of lambda parameter in the arc-length method
void determine_dlambda (double *ddr,double *v,double *ddrprev,long n,double ddlambda,double psi,double norfp,double dl,double &dlambda,long &stop);

/// function coputes coefficient of quadratic equation solved in the arclength method
void quadeqcoeff (double *ddr,double *v,long n,double ddlambda,double psi,double norfp,double dl,
		  double &a0,double &a1,double &a2);

/// function computes generalized norm of displacement increments
double displincr (double *displincr,long n);

/// function computes generalized norm of load increments
double loadincr (double *loadincr,long n);
*/

/// function saves actual step of arclength method
void arclsave (long fi,long n,double blambda,double dl,double *r,double *fp);

/// function restarts saved step of arclength method
void arclopen (long &fi,long &n,double &blambda,double &dl,double *r,double *fp);

/// function restarts Newton-Raphsom method from the file (used for analysis of Hinkley)
void newton_raphson_restart (long lcid);



#endif
