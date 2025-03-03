#ifndef TENSOR_H
#define TENSOR_H

#include "vector.h"
#include "matrix.h"


/// transformation of the second order tensor from the global coordinate system to the local one
void glob_loc_tens_trans (vector &tens, matrix &tmat);
/// transformation of the second order tensor from the local coordinate system to the global one
void loc_glob_tens_trans (vector &tens, matrix &tmat);


/// computes function of the second order tensor
void f_tensor(matrix &a, double (*f)(double), long nijac, double error, double zero, matrix &af);
/// computes function of the second order tensor
void f_tensor(matrix &a, double (*f)(double), matrix &t, matrix &af);

/// computes derivatioves i-th component of alpha-th principal direction vector with respect to the second order tensor
long dpdir_da(matrix &t, vector &pa, long alpha, long i, matrix &dpd_da);

/// computes single contraction of the second order tensors
void tensor_dot_prod(vector &a, vector &b, matrix &c);
/// computes single contraction of the two same second order tensors
void tensor_dot_prod(vector &a, vector &b);
/// computes double contraction of the second order tensors
double tensor_dbldot_prod (vector &a, vector &b, double k);
/// computes double contarction of the two second order tensors
double tensor_dbldot_prod (matrix &a, matrix &b);

#endif
