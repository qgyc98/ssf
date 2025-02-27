#ifndef MTSOLVER_H
#define MTSOLVER_H

#include <stdio.h>
#include "selection.h"
#include "timecontr.h"
#include "mtglvec.h"


/// function invokes solver of time dependent problems
void solve_time_dep_prob ();

/// solver of time dependent problems, viscous models can be used
void visco_solver (long lcid);

/// solver of time dependent problems built on one step concept
void visco_solver2 (long lcid);

/// initialization phase of visco_solver for one_step concept
void visco_solver_init(long lcid, long &rest_calc, mt_glob_vec &mt_gv);

/// function performs one time step of visco_solver2
long one_step (long lcid,double time, double dt, double &dtr, double prev_time, long rest_calc, 
               long istep, long li, mt_glob_vec &mt_gv);

/// function updates status of elements and nodes after restoreg from backup in growing mechanical problem
void update_elnod_stat_after_hdbrest(long lcid);

/// updates status of nodes and elements according to previous time step in growing mechanical problem
void update_elemnod_prev_time(double prev_time, long ncd, long nce);

/// updates status of nodes, DOFs and elements in growing mechanical problem
void update_elnod_stat(long lcid, long istep, double prev_time, long *ifn, double *r, double *fb, double *fp,
                       long &mnce, long &mncd, long &mnae);

/// calculates and hadles force %vector due to removed elements
long forces_due2removed_elem(long lcid, long nrel, long istep, double prev_time, double *fi, double *fp, 
                             double *flp, double *fb, double *r);


// FUNCTIONS OF HARDDISK BACKUP

/// function saves actual step of the solver and internal data for the future restart
void solver_save (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a text format
void solver_save_text (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a binary format
void solver_save_binary (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a text file (single file)
void solver_save_text_single (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a text file (multiple files)
void solver_save_text_multiple (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a binary file (single file)
void solver_save_binary_single (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);
/// function saves actual step of the solver and internal data for the future restart in a binary file (multiple files)
void solver_save_binary_multiple (double *r, double *fp, long ni, double time, double dt, timecontr *tc, long n);

/// function restores saved step of the solver and internal data and restarts computation
void solver_restore (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation (text format)
void solver_restore_text (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation (binary format)
void solver_restore_binary (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation from a text file
void solver_restore_text_single (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation from text files
void solver_restore_text_multiple (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation from a binary file
void solver_restore_binary_single (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n);
/// function restores saved actual step of the solver and internal data and restarts computation from binary files
void solver_restore_binary_multiple (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n);


#endif
