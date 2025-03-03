#ifndef CREEP_H
#define CREEP_H

#include "alias.h"

struct matrix;
struct vector;

void creep_matstiff (matrix &d, long ipp,long im,long ido);
double creep_matstiffchange (long ipp,long im,long ido);
void creep_initmaterialmodel(long ipp,long im,long ido);
void creep_updateval(long ipp,long im,long ido);
double creep_compute_actual_ym (long ipp,long im,long ido);
double creep_give_actual_ym (long ipp,long im,long ido);
double creep_compute_inital_ym (long ipp,long im,long ido);
double creep_give_actual_ft (long ipp,long im,long ido);
double creep_give_actual_fc (long ipp,long im,long ido);
void creep_givestressincr (long ipp,long im,long ido,long fi,vector &sig);
void creep_nlstressesincr (long ipp,long im,long ido);
void creep_incrtotstresses (long ipp,long im,long ido,vector &dsigma);
void creep_aeging_strains (long ipp,long im,long ido,vector &eps_ag);
void creep_nlstresses (long ipp,long im,long ido);
void creep_hidden_strains (matrix &screep,vector &sig,vector &emu, double mi, long n_ret_times, vector &ret_times, double actualtime, double dt,strastrestate ss);
void creep_giveirrstrains (long ipp, long im, long ido, vector &epsir);
long creep_number_rettimes(long ipp,long im);
long creep_ncompo (long ipp,long im);
void unit_compl_matrix (matrix &c,double nu,strastrestate ssst);

#endif
