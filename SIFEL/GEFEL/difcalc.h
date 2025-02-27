#ifndef DIFCALC_H
#define DIFCALC_H

#include "vector.h"
#include "matrix.h"
#include "galias.h"

void derivatives_1d (vector &dx,double &jac,vector &x,double xi);

void jac_1d (double &jac,vector &x,double xi);

void jac_1d_2d (double &jac, vector &x, vector &y, double xi);

void jac_1d_3d (double &jac, vector &x, vector &y, vector &z, double xi);

void derivatives_2d (vector &dx,vector &dy,double &jac,
		     vector &x,vector &y,double xi,double eta);

void jac_2d (matrix &jac,vector &x,vector &y,double xi,double eta);
void jac_2d (double &jac,vector &x,vector &y,double xi,double eta);

void sec_derivatives_2d (gtypebf bftype,vector &dxdx,vector &dxdy,vector &dydy,
			 vector &x,vector &y,double xi,double eta);

void derivatives_3d (vector &dx,vector &dy,vector &dz,double &jac,
		     vector &x,vector &y,vector &z,
		     double xi,double eta,double zeta);

void jac_3d (matrix &jac,vector &x,vector &y,vector &z,
	     double xi,double eta,double zeta);
void jac_3d (double &jac,vector &x,vector &y,vector &z,
	     double xi,double eta,double zeta);

void jac1d2d (double &jac,vector &x,vector &y,double xi);
void jac1d3d (double &jac,vector &x,vector &y,vector &z,double xi);
void jac2d3d (double &jac,vector &x,vector &y,vector &z,double xi,double eta);

void jac1d_2d (double &jac,vector &x,vector &y,double xi,long edid);
void jac2d_3d (double &jac,vector &x,vector &y,vector &z,double xi,double eta,long sid);

long point_natcoord_1d(double px, double py, vector &x, vector &y, long ni, double err, double &xi);
long point_natcoord_1d_3d(double px, double py, double pz, vector &x, vector &y, vector &z, long ni, double err, double &xi, double &aerr);
long point_natcoord_2d(double px, double py, vector &x, vector &y, long ni, double err, double &xi, double &eta, double &aerr);
long point_natcoord_3d(double px, double py, double pz, vector &x, vector &y, vector &z, long ni, double err, double &xi, double &eta, double &zeta, double &aerr);
#endif
