#ifndef ELEMSWITCH_H
#define ELEMSWITCH_H

#include "alias.h"
#include "galias.h"
#include "vector.h"
#include "matrix.h"
#include "ipmap.h"

/// function computes strains at integration points
void compute_ipstrains (long lcid);
void compute_strains (long lcid);
/// function computes strains at nodes with respect to setup in problem description
void compute_nodestrains (long lcid);
/// general function for strain computation
void computestrains (long lcid);
/// function computes strains directly at nodes
void nodestrains_comp (long lcid);

/// function computes stresses at integration points
void compute_ipstresses (long lcid);
///
void find_extreme_ipstresses (long lcid);
/// function computes stresses at nodes
void compute_nodestresses (long lcid);
/// general function for stress computation
void computestresses (long lcid);

/// function computes values of other array at nodes
void compute_nodeothers ();

/// function interpolates nodal values to integration points on one element
void elem_intpointval (long eid,vector &nodval,vector &ipval);

/// function interpolates nodal values into integration points
void intpointval (double *gv,nonmechquant nmq,double scale);
/// function interpolates nodal values into integration points using lower order approximation
void intpointval2 (double *gv,nonmechquant nmq);
/// function interpolates nodal values into integration points using lower order approximation
void elem_intpointval2 (long eid, vector &nodval, vector &ipval);

/// function interpolates arbitrary nodal quantity to the required inner point
double interpolate (long eid, vector &nodval, vector &coord);
/// function computes volume appropriate to integration points
void ipvolume ();

/// function computes contributions to internal forces from one element
void elem_internal_forces (long i,long lcid,vector &ifor);
/// function computes stresses
void elem_nlstresses(long i,long lcid);
/// function computes local values for consecutive averaging
void elem_local_values (long i,long lcid);
/// function computes contributions to internal forces from one element nonlocally
void elem_nonloc_internal_forces (long i,long lcid,vector &ifor);
/// function computes contributions to increments of internal forces from one element
void elem_incr_internal_forces (long i,long lcid,vector &ifor);
/// function computes nodal forces caused by eigenstrains on one element
void elem_eigstrain_forces (long lcid,long eid,vector &nfor);

/// function integrates arbitrary selected quantity over element volume
void elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv);

/// function computes stiffness %matrix of required element
void stiffmat (long lcid,long eid,matrix &sm);
/// function computes bd %matrix of required element
void bdmatrix (long lcid,long eid,matrix &bd);
/// function computes dd %matrix of required element
void ddmatrix (long lcid,long eid,matrix &dd);
/// function computes db %matrix of required element
void dbmatrix (long lcid,long eid,matrix &db);
/// function computes mass %matrix of required element
void massmat (long lcid,long eid,matrix &mm);
/// function computes intial stiffness %matrix of required element
void initstiffmat (long lcid,long eid,matrix &sm);
/// function computes damping %matrix of required element
void dampmat (long lcid,long eid,matrix &dm);
/// function computes load %matrix of required element
void loadmat (long eid,matrix &lm);
/// function assembles transformation matrix of the element according to local coordinate systems defined in element nodes
void elem_transf_matrix (long eid, ivector &enodes, matrix &tmat);

/// function determines mechanical quantities at nodes on one element, they are used by TRFEL in coupled problems
void elem_mechq_nodval (long eid, vector &nodval, nontransquant ntq);
/// function determines mechanical quantities at nodes on one element, they are used by TRFEL in coupled problems
/// linear approx. functions are needed
void elem_mechq_nodval2 (long eid, vector &nodval, nontransquant ntq);

/// function computes mechanical quantities in nodes on one element, they are used by TRFEL in coupled problems
void elem_mechq_nodval_comp(long eid, vector &nodval, long ncne, long nq, nontransquant *qt);
/// function computes mechanical quantities in nodes on one element, they are used by TRFEL in coupled problems
/// linear approx. functions are needed
void elem_mechq_nodval_comp2(long eid, vector &nodval, long ncne, long nq, nontransquant *qt);

/// function computes the centroid of the selected element
void centroid(long eid, vector &coord);

/// function computes global coordinates of the given integration point
void ipcoord (long eid,long ipp,long ri,long ci,vector &ipcoord);

/// function returns natural coordinates of the given integration point
void ipncoord (long eid,long ipp,vector &ipncoord);

/// function returns strain-displacement (geometric) %matrix for the given point on element
void geom_matrix(long eid, double xi, double eta, double zeta, matrix &gm);

/// function computes actual values of eigenstrains/eigenstresses at auxliary integration points
void compute_aipeigstr(double time, long n, ipmap *ipm);

/// function computes strains at auxiliary integration points
void compute_aipstrains(long lcid, long n, ipmap *ipm);

/// function computes stresses and eventual state variables at auxiliary integration points
void compute_aipstresses(long lcid);

/// function selects components from the macro level to meso level for homogenization problems 
void higher_to_lower_level_elemm (long eid,long *counter,double *buff);

#endif
