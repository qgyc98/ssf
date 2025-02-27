#ifndef FCSOLVER_H
#define FCSOLVER_H

#include <stdio.h>
struct vector;
struct ivector;

void solve_fcouplprob ();

void nonlinear_solver_coupl (long lcid);
void newton_raphson_coupl_tkr (long lcid);
void newton_raphson_coupl (long lcid);
void newton_raphson_coupl_vform (long lcid);
void newton_raphson_coupl_dform (long lcid);

void assemble_dof_block_id(ivector &dofbid);
void compute_vector_sel_norms(const double *c, const ivector &dofbid, vector &norms);
long check_tolerance(const vector &vnorfa, const vector &vnorfb, const vector &tol, vector &aterr);
void vector_assemb (double *c,double *m,double *t);
void vector_decomp (double *c,double *m,double *t);
void print_log();

void newton_raphson_fully_coupled (long lcid);

#endif
