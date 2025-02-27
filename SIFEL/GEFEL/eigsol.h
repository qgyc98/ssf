#ifndef EIGSOL_H
#define EIGSOL_H

void jacobi_rot (double *a,double *evec,double *eval,long n,long ni,long &ani,double limit, long normalize);
void gen_jacobi (double *a,double *b,double *x,double *w,long n,long ni,
		 double *thresholds,long nthr,double zero);
void gen_jacobi2 (double *k,double *m,double *x,double *w,long n,long ni,
		 double *gate,long ng,double limit);

#endif
