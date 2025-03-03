#ifndef BASEFUN_H
#define BASEFUN_H

#include "galias.h"
struct vector;

void bf_lin_1d (double *n,double x);
void dx_bf_lin_1d (double *n);

void bf_quad_1d (double *n,double x);
void dx_bf_quad_1d (double *n,double x);

void bf_cubic_1d (double *n, double ksi);
void dksi_bf_cubic_1d(double *n,double ksi);

void bf_1d_2d(vector &p, vector &x, vector &y, double xi);
void bf_1d_3d(vector &p, vector &x, vector &y, vector &z, double xi);


void dilat (double *n,double xi);
void der_dilat (double *n,double l);

void defl2D_fun (double *n,double xi,double l,double kappa);
void der_defl2D_fun (double *n,double xi,double l,double kappa);
void roty_fun (double *n,double xi,double l,double kappa);
void der_roty_fun (double *n,double xi,double l,double kappa);

void a_coeff (double *a,double *x,double *y);
void b_coeff (double *b,double *y);
void c_coeff (double *c,double *x);

void plsa (double *a,double *x,double *y,double det);
void plsb (double *b,double *y,double det);
void plsc (double *c,double *x,double det);

void ac_bf_quad_3_2d (double *n,double *l);
void ac_dx_bf_quad_3_2d (double *n,double *b,double *l);
void ac_dy_bf_quad_3_2d (double *n,double *c,double *l);


void bf_lin_3_2d (double *n,double x,double y);
void dx_bf_lin_3_2d (double *n);
void dy_bf_lin_3_2d (double *n);
void bf_quad_3_2d (double *n,double x,double y);
void dx_bf_quad_3_2d (double *n,double x,double y);
void dy_bf_quad_3_2d (double *n,double x,double y);


void bf_lin_4_2d (double *n,double x,double y);
void dx_bf_lin_4_2d (double *n,double y);
void dy_bf_lin_4_2d (double *n,double x);
void dxdx_bf_lin_4_2d (double *n);
void dxdy_bf_lin_4_2d (double *n);
void dydy_bf_lin_4_2d (double *n);


void bf_quad_4_2d (double *n,double x,double y);
void dx_bf_quad_4_2d (double *n,double x,double y);
void dy_bf_quad_4_2d (double *n,double x,double y);
void dxdx_bf_quad_4_2d (double *n,double y);
void dydy_bf_quad_4_2d (double *n,double x);
void dxdy_bf_quad_4_2d (double *n,double x,double y);


void bf_cubic_4_2d (double *n,double ksi,double eta);
void dksi_bf_cubic_4_2d (double *n,double ksi,double eta);
void deta_bf_cubic_4_2d (double *n,double ksi,double eta);

void bf_rot_3_2d (double *n,double *l,double *b,double *c);
void dx_bf_rot_3_2d (double *n,double *l,double *b,double *c,double area);
void dy_bf_rot_3_2d (double *n,double *l,double *b,double *c,double area);

void bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l);
void dx_bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l);
void dy_bf_rot_4_2d (double *n,double x,double y,double *nx,double *ny,double *l);

void bf_2d(vector &p, vector &x, vector &y, double xi, double eta);

void bf_cct (double *n,double *l,double *sx,double *sy);
void dx_cct (double *n,double *l,double *b,double *sx,double *sy);
void dy_cct (double *n,double *l,double *c,double *sx,double *sy);


void bf_quad_dkq (double *n,double x,double y,double *l,double *sx,double *sy,double *nx,double *ny);
void dx_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs);
void dy_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs);
void dxdx_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs);
void dydy_bf_quad_dkq (double *n,double xi,double eta,vector &x,vector &y,double *l,double *sx,double *sy,double *nx,double *ny,double *signs);
void bf_cubic_plate_4_2d (double *n,double x,double y);


void a_coeff_3d (double *a,double *x,double *y,double *z,double v);
void b_coeff_3d (double *b,double *y,double *z,double v);
void c_coeff_3d (double *c,double *x,double *z,double v);
void d_coeff_3d (double *d,double *x,double *y,double v);

void vola_3d (double *a,double *x,double *y,double *z,double det);
void volb_3d (double *b,double *y,double *z,double det);
void volc_3d (double *c,double *x,double *z,double det);
void vold_3d (double *d,double *x,double *y,double det);

void bf_lin_tet (double *n,double x,double y,double z);
void dx_bf_lin_tet (double *n);
void dy_bf_lin_tet (double *n);
void dz_bf_lin_tet (double *n);

void bf_quad_tet (double *n,double x,double y,double z);
void dx_bf_quad_tet (double *n,double x,double y,double z);
void dy_bf_quad_tet (double *n,double x,double y,double z);
void dz_bf_quad_tet (double *n,double x,double y,double z);

void bf_lin_wed_3d (double *n,double x,double y,double z);
void dx_bf_lin_wed_3d (double *n,double z);
void dy_bf_lin_wed_3d (double *n,double z);
void dz_bf_lin_wed_3d (double *n,double x,double y);

void bf_quad_wed_3d (double *n,double x,double y,double z);
void dx_bf_quad_wed_3d (double *n,double x,double y,double z);
void dy_bf_quad_wed_3d (double *n,double x,double y,double z);
void dz_bf_quad_wed_3d (double *n,double x,double y,double z);


void bf_lin_hex_3d (double *n,double x,double y,double z);
void dx_bf_lin_hex_3d (double *n,double y,double z);
void dy_bf_lin_hex_3d (double *n,double x,double z);
void dz_bf_lin_hex_3d (double *n,double x,double y);

void bf_quad_hex_3d (double *n,double x,double y,double z);
void dx_bf_quad_hex_3d (double *n,double x,double y,double z);
void dy_bf_quad_hex_3d (double *n,double x,double y,double z);
void dz_bf_quad_hex_3d (double *n,double x,double y,double z);

void bf_quad_hexrot_3d (double *n,double x,double y,double z);
void dx_bf_quad_hexrot_3d (double *n,double x,double y,double z);
void dy_bf_quad_hexrot_3d (double *n,double x,double y,double z);
void dz_bf_quad_hexrot_3d (double *n,double x,double y,double z);
void bf_3d(vector &p, vector &x, vector &y, vector &z, double xi, double eta, double zeta);

void corr_nat_coord_bounds(gtypel et, double &xi, double &eta, double &zeta);
void corr_nat_coord_tr(double &xi, double &eta);
void corr_nat_coord_tet(double &xi, double &eta, double &zeta);
#endif
