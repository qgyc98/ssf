#ifndef ELEMSWITCHT_H
#define ELEMSWITCHT_H

#include "aliast.h"
#include "galias.h"
#include "matrix.h"
class ipmap;

///  computes a conductivity matrix of an element
void conductmat (long eid,long lcid,matrix &km);
///  computes a capacity matrix of an element
void capacmat (long eid,long lcid,matrix &cm);
///  computes a reaction matrix of an element (in reaction-diffusion problems)
void reactmat (long eid,long lcid,matrix &rm);
///  computes an advection matrix of an element
void advectmat (long eid,long lcid,matrix &hm);
/// computes a gradient %matrix at the given point of element
void grad_matrix (long eid, vector &x, vector &y, vector &z, double xi, double eta, double zeta, matrix &gm);

void intpointvalues (long eid);
void intpointvalues_puc (long eid);
void initintpointvalues (long eid);
void intpointgradients (long eid);
void intpointfluxes (long eid);
void averageflux (long lcid,long eid,vector &fl);
void intpointothers (long eid);
void intpointeqothers (long eid);
/// function computes approximated value of the qunatity given by nodal values
double approx(long eid, double xi, double eta, double zeta, vector &nv);

void source_vector (long lcid,long eid,vector &nodval,vector &lv);

//void compute_ipfluxes ();
//void compute_ipgrads ();
//void compute_ipotherst ();
//void compute_ipeqotherst ();

//void compute_ipfluxes ();
//void compute_ipgrads ();
//void compute_ipotherst ();
//void compute_ipeqotherst ();

// compute quantities for auxiliary integration points

/// computes fluxes at auxiliary integration points
void compute_ipfluxes ();
/// computes gradients of primary unknowns at auxiliary integration points
void compute_ipgrads ();
/// computes other values/state variables at auxiliary integration points
void compute_ipotherst ();
/// computes eqother values/state variables at auxiliary integration points
void compute_ipeqotherst ();

void compute_nodegrads ();
void compute_nodefluxes ();
void compute_nodeotherst ();
void compute_nodeotherst_comp ();
void compute_nodeeqotherst ();
void compute_nodeeqotherst_comp ();

void internal_fluxes (double *intflux,long n);
void nodal_energy (double *nodener,long n,double dt);
void compute_average_fluxes (long lcid, vector &fl);

void elem_neumann_vector (vector &lv,long lcid,long eid,long i);
void elem_newton_vector (vector &lv,long lcid,long eid,long i);
void elem_transmission_flux (vector &lv,long lcid,long eid,long i);
void volume_rhs_vector (vector &lv,long lcid,long eid);
void volume_rhs_vector2 (vector &lv,long lcid,long eid);
void lmat (long eid,long lcid,matrix &lm);
void ltmat (long eid,long lcid,matrix &lm);
void averdmat (long eid,double &elemarea,matrix &km);
void avercmat (long eid,double &elemarea,matrix &cm);

double give_elemarea (long eid);

void higher_to_lower_level_elem (long eid,long *counter,double *buff);
void intpointvalt (double *gv, nontransquant ntq, double scale);
void elem_intpointvalt (long eid,vector &nodval,vector &ipval);
void give_transq_nodval (double *gv,nonmechquant nmq);
void elem_transq_nodval (long eid, vector &nodval, nonmechquant nmq);
void elem_transq_init_nodval (long eid, vector &nodval, nonmechquant nmq);
void elem_transq_nodval_comp (long eid, vector &nodval, long ncne, long nq, nonmechquant *qt);
void elem_transq_init_nodval_comp (long eid, vector &nodval, long ncne, long nq, nonmechquant *qt);

///  function returns integral over a single finite element
///  variable is stored in nodes
double elem_total_integral (long eid,vector &nodval);
///  function returns integral over a single finite element
///  variable is stored in integration points
double elem_total_integral_ip (long eid,long varid);
/// function computes global coordinates of the given integration point on element
void ipcoordt (long eid, long ipp, long ri, long ci, vector &coord);
/// function returns natural coordinates of the given integration point on element
void ipncoordt (long eid, long ipp, vector &ncoord);
/// function computes the centroid of the selected element
void centroidt(long eid, vector &coord);

/// function computes variable values at auxiliary integration points
void compute_aipvals (long n, ipmap *ipm);

/// function computes gradients at auxiliary integration points
void compute_aipgrads(long lcid, long n, ipmap *ipm);

/// function computes fluxes and eventual state variables at auxiliary integration points
void compute_aipfluxes();

/// function computes other values at auxiliary integration points
void compute_aipotherst(long n, ipmap *ipm);

/// function computes surface flux
void elem_surface_flux (long eid,long lcid,long beid,double *fluxes);

#endif
