#ifndef BACKUPSOLT_H
#define BACKUPSOLT_H

#include "timecontr.h"

void solvert_save (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_text (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_binary (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_text_single (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_text_multiple (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_binary_single (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_save_binary_multiple (double *r, double *dr, double *fp, long ni, double time, double dt, timecontr &tc, long n);
void solvert_restore (double *r, double *dr, double *fp, long &ni, double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_text (double *r, double *dr, double *fp, long &ni, double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_binary (double *r, double *dr, double *fp, long &ni, double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_text_single (double *r, double *dr, double *fp, long &ni, double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_text_multiple (double *r, double *dr, double *fp, long &ni, double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_binary_single (double *r,double *dr, double *fp,long &ni,double &time, double &dt, timecontr &tc, long &n);
void solvert_restore_binary_multiple (double *r,double *dr, double *fp,long &ni,double &time, double &dt, timecontr &tc, long &n);


#endif
