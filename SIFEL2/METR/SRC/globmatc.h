#ifndef GLOBMATC_H
#define GLOBMATC_H

#include "global.h"
#include "globalt.h"
#include "globalc.h"

struct matrix;
struct ivector;

void zero_order_matrix (long lcid);
void first_order_matrix (long lcid);

void zero_order_matrix_fc ();
void first_order_matrix_fc ();
void residual_vector ();

void internal_gforces (long lcid,double *intfor);
void incr_internal_gforces (long lcid,double *intfor);
void elemvaluesc (long lcid,double *r,long *cn,long ndofe);

void metr_right_hand_side (long lcid,double *rhs,double *flv);
void right_hand_side (long lcid,double *rhs,long n);
void metr_right_hand_side_fc (double *rhs);

void updateval();
void approximationc ();
void approximationcoup ();
void approximationc_fc ();

void internal_coup_fluxes (double *intflux,long n);
void internal_coup_forces (long lcid,double *intflux);
void right_hand_side (long lcid,double *rhs,long n);
void assemble_init_coup (double *rhs);

void pass_coup_data(long lcid);

void metr_mefel ();

void mefel_trfel (long lcid);
void mefel_trfel_by_nodes(void);
void mefel_trfel_by_nodes_comp(void);
void mefel_trfel_by_aip(long lcid, long n, ipmap *ipm);

void trfel_mefel ();
void trfel_mefel_by_nodes(void);
void trfel_mefel_by_nodes_comp(void);
void trfel_mefel_by_aip(long n, ipmap *ipm);

void init_trfel_mefel();
void init_trfel_mefel_by_nodes(void);
void init_trfel_mefel_by_nodes_comp(void);
void init_trfel_mefel_by_aip(long n, ipmap *ipm);

void actualize_aip_nonmechq(long n, ipmap *ipm);
void actualize_aip_nontransq(long n, ipmap *ipm);

void metr_ip_mapping(long ni, double err);
void trfel_mefel_ip_mapping(long ni, double err);
void mefel_trfel_ip_mapping(long ni, double err);

void trfmef_give_transq_nodval (double *gv, long *nodmap, nonmechquant nmq);
void create_trfmef_nod_map(long *nodmap);

// functions for the direct mapping of passed values by int. point indeces
void mefel_trfel_ip_mapping_old(long *tm_ip, long *mt_ip);
void mefel_trfel_copyip(void);
void trfel_mefel_copyip(void);
void init_trfel_mefel_copyip(void);

// obsolete functions for the passing values, used in old solvers
void init_trfel_mefel_orig ();
void nodal_nodal_values (double *gv, long *nodmap, nonmechquant nmq);
void approximation_temper ();
void approximation_humid ();
void approximation_inittemper ();

void eldispl_fc (long lcid,long eid,double *r);

#endif
