#ifndef GLOBMATT_H
#define GLOBMATT_H

#include "genfile.h"
#include "aliast.h"
#include "ipmap.h"

void conductivity_matrix (long lcid);

void capacity_matrix (long lcid);

//void residuum (double *r,double *p,double *v,double dt,long n,long lcid);

///  computes/extracts all nodal value at the given PUC element
void elemvalues_puc (long eid, vector &r);
/// adds contribution from master node to nodal value %vector of the given PUC element
void nodalval_puc_mastercont(long nid, vector &r);
///  computes/extracts all nodal value at the given element
void elemvalues (long eid, vector &r);

///  computes/extracts all nodal value at the given node
void nodalval (long nid, vector &r);
/// extracts all nodal values at the given regular node
void select_nodalval (long nid, vector &r);

/// computes/extracts nodal value at the given DOF at the given node
double nodalval (long nid, long dofid);
/// extracts nodal value at the given DOF at the given regular node
double select_nodalval (long nid, long dofid);

/// computes/extracts initial nodal value at the given node
void initnodval (long idn, vector &r);
/// extracts initial initial nodal value at the given regular node
void select_initnodval (long idn, vector &r);

/// computes/extracts initial nodal value at the given node
void initnodval2(long idn, vector &r);
/// extracts initial initial nodal value at the given regular node
void select_initnodval2(long idn, vector &r);

/// computes/extracts nodal value from given global nodal value %vector at one node
void gvnodval (long idn, double *lhs, vector &r);
/// extracts nodal value from given global nodal value %vector at one regular node
void select_gvnodval (long idn, double *lhs, vector &r);

/// extracts all nodal values at one node from global value %vector
void  gen_gvnodval(double *ifor, long nid, vector &nf);
/// extracts given DOF value from flux resultant %vector at one node
double gen_gvnodval(double *ifor, long nid, long dofid);

/// extracts other values from nodes whose numbers are defined in the %vector nod
void nodalotherval (ivector &nod, vector &r);

/// computes/extracts first time derivative of nodal values for given DOF id
void nodalderivatives (long eid, vector &r);

void prescvalues (double *r,long *cn,long ndofe);
void prevprescvalues (double *r,long *cn,long ndofe);
void initialvalues (long eid, vector &r);


void approximation ();
void initapproximation ();
void aip_approximation (long n, ipmap *ipm);
void aip_initapproximation (long n, ipmap *ipm);
void approximation_puc ();
void actual_previous_change ();
void actual_previous_nodval ();

void assemble_init (double *rhs);
void trfel_right_hand_side (long lcid,double *rhs,long n);
void trfel_right_hand_side2 (long lcid,double *rhs,long n);
void trfel_bound_flux (long lcid,double *rhs,long n);

void compute_req_valt (long lcid);


void give_nodal_humid (double *gv, long *nodmap);

void solution_correction ();

void compute_cycles (double *plhs,double *pplhs);

void copy_nodval (long lcid);

void assemble_gradients_contrib(double *rhs);
void assemble_l_matrix (matrix &lm, matrix &ltm);
void assemble_l_matrix (matrix &lm);
void assemble_lt_matrix (matrix &ltm);
void assemble_conductivity_matrix (matrix &km);
void assemble_average_d_matrix (matrix &km, double &area);
void assemble_average_c_matrix (matrix &cm);

///  function computes integral of a variable over the whole domain
///  the variable is stored on nodes
double total_integral (long lcid);

///  function computes integral of a variable over the whole domain
///  the variable is stored in integration points
double total_integral_ip (long varid);

void dt_subdom ();

void lnso_leso_setting (long *lsso);

void surface_fluxes (FILE *out);
double surface_fluxes (long lcid, long nbf);
#endif
