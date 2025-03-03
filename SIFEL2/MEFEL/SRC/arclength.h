#ifndef ARCLENGTH_H
#define ARCLENGTH_H


#include "galias.h"
#include "nonlinman.h"
#include <stdio.h>



/**
  The macro is used for coefficient debugging of quadratic equations in the arclength method.
  It can be called only from the 
*/
#define quadeqcoef_log() \
{\
  fprintf (Out,"\n\n\n Check of quadratic equation");\
  fprintf (Out,"\n norfp     %15.10le",norfp);\
  fprintf (Out,"\n norv      %15.10le",norv);\
  fprintf (Out,"\n (ddr,v)   %15.10le",norddrv);\
  fprintf (Out,"\n (ddr,ddr) %15.10le",norddr);\
  fprintf (Out,"\n dl        %15.10le",dl);\
  fprintf (Out,"\n ddlambda  %15.10le",ddlambda);\
  fprintf (Out,"\n psi       %15.10le",psi);\
  fprintf (Out,"\n a2        %15.10le",a2);\
  fprintf (Out,"\n a1        %15.10le",a1);\
  fprintf (Out,"\n a0        %15.10le",a0);\
  fprintf (Out,"\n discrim   %15.10le",a1*a1-4.0*a2*a0);\
  fprintf (Out,"\n\n");\
}




/// The main procedure which solves the system of the nonlinear equations by the arclength method.
double garclength(long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp, 
                  double *flc, double *flp, long li, double ilambda, answertype outres);

/// The function assembles the stiffness matrix in dependence on the setup of the arclength method.
long assemble_stiffness_matrix (long lcid,long i,long j,long li);

/// The function initializes selected dofs used for the arclength method controling.
void seldofinit();

/// The function calculates generalized norm of the load vector.
double loadincr (long lcid, nonlinman *nlman, double lambda, double *fp, long n);

/// The function computes square of the generalized norm of displacement vector.
double displincr (long lcid, long istep, nonlinman *nlman, double *dr1, double *dr2, long n);

/// The function returns prescribed displacements for the i-th seleceted dof in nonlinman
double prdisplincr(long i, long lcid, long istep);

/// The function computes coefficiets of quadratic equation solved in the spherical arclength method.
void quadeqcoeff (long lcid, long istep, nonlinman *nlman, double *ddr,double *v,long n,double ddlambda,double psi,double norfp,double dl,
		  double &a0,double &a1,double &a2);

/// The function computes increment of load coefficient.
void determine_dlambda (double *ddr, double *u, double *v, double *ddrprev, long n, 
                        double ddlambda, double psi, double norfp, double dl, 
                        double &dlambda, long &stop, nonlinman *nlman, long lcid, long istep);

///  The function modifies the psi coefficient according to the  ddlambda increment adaptively.
void adapt_psi_coeff(long i, double ddlambda, double &ddlambda0, double &psi0, double &psi);

#endif
