#ifndef GLOBMAT_H
#define GLOBMAT_H

#include <stdio.h>

struct matrix;
struct vector;
//class problem;


void stiffness_matrix (long lcid);
void mass_matrix (long lcid);
void initial_stiffness_matrix (long lcid);
void damping_matrix (long lcid);
void layered_stiffness_matrix ();

void internal_forces (long lcid,double *intfor);
void loc_internal_forces (long lcid,double *intfor);
void nonloc_internal_forces (long lcid,double *intfor);
void incr_internal_forces (long lcid,double *intfor);
void lagrmultcontr_intforces (long lcid,double *intfor);

void nodal_eigstrain_forces (long lcid,double *nodfor,double time);
void elem_eigstrain_forces (long lcid,long eid,vector &nfor);

void eldispl (long lcid,long eid,double *r);
void elprdispl (long lcid,long eid,double *r);
void elprevprdispl (long lcid,long eid,double *r);

void noddispl (long lcid,double *r, long nid);
void select_noddispl (long lcid,double *r,long nid);

/// extracts force compoennts of one node from the global force %vector
void nodforce (double *f, double *nf, long nid);

double macrostrains(long lcid, long cn);
void macrostrains(long lcid, vector &meps);
double macrostresses(long lcid, long cn);
void macrostresses(long lcid, vector &msig);


void nodforce (double *f, long nid, vector &nf, bool react=true);
void select_nodforce(double *f, long nid, vector &nf, bool react=true);


void constr_matrix (long nid,long cid,matrix &lcm);
void mefel_right_hand_side (long lcid,double *rhs,double *flv=NULL);
void compute_req_val (long lcid);
void compute_reactions(long lcid);

void local_global_displ_transf (long lcid);
void stress_initdispl(long lcid);

void print_nodval_vec(const double *v, const char *label, FILE *out);
void check_hex_iface_elem();

#endif
