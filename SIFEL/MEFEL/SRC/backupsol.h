#ifndef BACKUPOL_H
#define BACKUPOL_H

#include "timecontr.h"

//
// FUNCTIONS OF HARDDISK BACKUP
//

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
