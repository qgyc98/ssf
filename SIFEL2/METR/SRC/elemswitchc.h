#ifndef ELEMSWITCHC_H
#define ELEMSWITCHC_H

#include "aliasc.h"
#include "matrix.h"

void upper_cond_coupl_mat (long eid,long lcid,matrix &vm);
void lower_cond_coupl_mat (long eid,long lcid,matrix &vm);
void upper_cap_coupl_mat (long eid,long lcid,matrix &vm);
void lower_cap_coupl_mat (long eid,long lcid,matrix &vm);

void zero_order_matrix (long eid,matrix &km);
void first_order_matrix (long eid,matrix &cm);

void intpointvaluesc (long eid);
void intpointstrainsc (long eid);
void intpointgradientsc (long eid);

void assemble_coup (long lcid,double *rhs,long n);

void volume_rhs_vectorc (vector &lv,long lcid,long eid);


void internal_coup_fluxes (long lcid,double *intflux);
void internal_coup_forces (long lcid,double *intflux);

#endif
