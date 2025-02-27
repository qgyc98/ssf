#ifndef VECTTENS_H
#define VECTTENS_H

#include <stdio.h>
#include "alias.h"

struct matrix;
struct vector;
struct ivector;



/// converts %vector to second order tensor
void vector_tensor (vector &v,matrix &t,strastrestate ssst,strastre ss);
/// converts second order tensor to %vector
void tensor_vector (vector &v,matrix &t,strastrestate ssst,strastre ss);
/// converts second order tensor to %vector in the full form
void tensor_vector_full (vector &v, matrix &t, strastre ss);

/// creates auxiliary matrix (useful for thermal strains,...)
void tensor_vector_matrix (strastrestate ssst,matrix &m);

/// converts fourth order tensor to reduced material %matrix
void tensor4_matrix (matrix &m, const matrix &t, strastrestate ssst);
/// converts fourth order tensor to reduced element %matrix
void tensor4_ematrix (matrix &m, const matrix &t, strastrestate ssst);
/// converts %matrix to fourth order tensor
void matrix_tensor4 (const matrix &m, matrix &t, strastrestate ssst);

/// converts matrices with 6 row components to matrices with ncompstr number of rows, the number of columns is fixed
void red_rows_matrix (matrix &m, matrix &t, strastrestate ssst);
/// converts matrices with 6 column components to matrices with ncompstr number of columns, the number of rows is fixed
void red_cols_matrix (matrix &m, matrix &t, strastrestate ssst);

/// returns guess of stress/strain state depending on the number of components
strastrestate guess_ssst(long ncomp);

/// returns number of stress/strain array components depending on the stress/strain indicator
long give_ncompstr(strastrestate ssst);

// converts long %vector to reduced %vector
//void longvect_shortvect (double *lv,double *sv,strastrestate ssst);

/// converts full-length stress/strain %vector to reduced %vector
void give_red_vector(const vector &fv, vector &rv, strastrestate ssst);

/// converts reduced stress/strain %vector to full-length %vector
void give_full_vector(vector &fv, const vector &rv, strastrestate ssst);

/// returns vector of indeces of shear components in the reduced vector/matrix notation
void give_shear_indices(strastrestate ssst, ivector &id);

/// returns vector of indeces of of all shear components in the reduced vector/matrix notation, -1 is set for not used components
void give_all_shear_indices(strastrestate ssst, ivector &id);

/// returns vector of indeces of of all normal components in the reduced vector/matrix notation, -1 is set for not used components
void give_all_normal_indices(strastrestate ssst, ivector &id);

/// converts full-length stress/strain %vector to reduced %vector
void give_red_vector(double *fv, double *rv, strastrestate ssst);

/// converts reduced stress/strain %vector to full-length %vector
void give_full_vector(double *fv, double *rv, strastrestate ssst);

/// converts reduced stress/strain %vector to full-length %vector
void give_full_vector(vector &fv, double *rv, strastrestate ssst);

/// transforms stress/strain %vector l given in the local coordinate system to stress/strain %vector g in global coordinate system (g = T . l)
void lg_engvectortransf (vector &g, const vector &l, const matrix &tmat, strastrestate ssst, strastre str);

/// transforms stress/strain %vector g given in the global coordinate system to stress/strain %vector l in local coordinate system (l = T^T . g)
void gl_engvectortransf (const vector &g, vector &l, const matrix &tmat, strastrestate ssst, strastre stra);

/// transforms stress/strain %vector g given in the global coordinate system to stress/strain %vector l in local coordinate system (l = T^T . g) and returns i-th component of l
void gl_comp_engvectortransf (const vector &g, double &l, long i, const matrix &tmat, strastrestate ssst, strastre stra);

/// transforms 4-th order tensor given in the local coordinate system to global one 
void lg_tens4transf (matrix &g, const matrix &l, const matrix &tmat, strastrestate ssst);

/// transforms 4-th order tensor given in the global coordinate system to local one 
void gl_tens4transf (const matrix &g, matrix &l, const matrix &tmat, strastrestate ssst);

///  function computes second invariant of deviator of stress tensor
///  components of the stress tensor are used
double j2_stress_invar (vector &lv);
double j2_strain_invar (vector &lv);
void deviator (vector &tens, vector &dev);

double first_invar (vector &tens);
double second_stress_invar (vector &tens);
double third_stress_invar (vector &tens);
double second_strain_invar (vector &tens);
double third_strain_invar (vector &tens);
double tensor_stress_norm(vector &sig);
double tensor_strain_norm(vector &eps);
void normed_stress_tensor(vector &t, vector &nt);
void normed_strain_tensor(vector &eps, vector &neps);

double first_invar (matrix &t);
double second_invar (matrix &t);
double third_invar (matrix &t);

#endif
