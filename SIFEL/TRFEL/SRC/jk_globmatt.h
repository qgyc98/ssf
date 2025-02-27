#ifndef GLOBMATT_H
#define GLOBMATT_H

#include "genfile.h"
#include "aliast.h"

void conductivity_matrix (long lcid);

void capacity_matrix (long lcid);

void residuum (double *r,double *p,double *v,double dt,long n,long lcid);

void nodalvalues (long lcid,double *r,long *cn,long ndofe);
void nodval (long lcid,double *r, long idn);
void nodalderivatives (double *r,long *cn,long ndofe);
void prescvalues (double *r,long *cn,long ndofe);
void initialvalues (double *r,long *cn,long ndofe);

void conductmat (long eid,long lcid,matrix &km);
void capacmat (long eid,long lcid,matrix &cm);

void internal_fluxes (double *intflux,long n);

void approximation ();
void intpointvalues (long eid);
void intpointgradients (long eid);
void assemble_init (double *rhs);
void trfel_right_hand_side (long lcid,double *rhs,long n);
void trfel_bound_flux (long lcid,double *rhs,long n);

void compute_req_valt (long lcid);

void source_vector (long lcid,long eid,vector &nodval,vector &lv);

void give_nodal_humid (double *gv,long mnt);

void solution_correction ();

#endif
