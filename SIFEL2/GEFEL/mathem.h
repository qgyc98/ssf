#ifndef MATHEM_H
#define MATHEM_H

#include "vector.h"
#include "matrix.h"
#include "basefun.h"

#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif

class gtopology;
struct vector;

double radius (const vector &x, const vector &natcoord);
double length (const vector &coorda, double xb, double yb, double zb);
double length (const vector &coorda, const vector &coordb);
double sqr(double x);
double pow3(double x);
double linhex_volume(vector &x, vector &y, vector &z);

void nodal_values (double **val,vector &nx,vector &ny,vector &nz,
                   double *lhs,long dim,long fi,long ncomp);

#ifndef max2
 #define max2(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef min2
 #define min2(a,b) ((a)<(b)?(a):(b))
#endif

double sgn (double i);
double heaviside (double x);
double frexp10(double arg, int * exp);
double genheaviside (double x,double eps);
double polynom_4 (double x,double *a);

void sort_2 (long *x);
void sort_3 (long *x);
void sort_4 (long *x);
void sort_3 (double *x);

long solv_polynom_2 (double a,double b,double c,double &r1,double &r2);
long solv_polynom_3 (double a, double b, double c, double d, double &r1, double &r2, double &r3);
long solv_polynom_4 (double *coeff,double a,double b,double acc,double *roots);

long solv_1le (double a,double b,double &x);
long solv_2le (double *a, double *b, double *c,double &x,double &y);
long solv_2nle (double *a, double *b, double *c, double *d,double *x,double *y);

void nc_lin_3_2d (double xx,double yy,double *x,double *y,double &xi,double &eta);
void nc_quad_3_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta);
void nc_lin_4_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta);
void nc_quad_4_2d (double acc,double xx,double yy,double *x,double *y,double &xi,double &eta);

void interpolelem (gtopology *gt,long &nli,long **&lin,long &ntr,long **&trn,long &nqu,long **&qun,long &nte,long **&ten,long &ncu,long **&cun,long &nmp,double **&mpn,char flag);
void print_ex (gtopology *gt, const char *file, double *valnod, double *valel);
void print_dx (gtopology *gt, const char *file,double **valnod,double **valel,char tve, long *dimindex, char **caption,long nindex);

void maxmin_3 (double *x,double &max,double &min);
void maxmin_4 (double *x,double &max,double &min);


double lsm_quad(matrix &a, vector &r, vector &l, double x, double y, double x_old, double y_old, double zero, long solc);

int check_math_err();
int check_math_errel(long eid);
int check_math_errip(long eid, long ipp);
long test_math_err();
int test_math_errel(long eid);
int test_math_errip(long eid, long ipp);

#endif
