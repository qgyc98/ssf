#ifndef HERMTRF_H
#define HERMTRF_H

///  The function returns required conductivity coefficient from the material conductivity %matrix.
void condmat(double phi, double t, double *cm, long ri, long ci, long im);

///  The function returns required capacity coefficient from the material capacity %matrix.
double capcoeff(double phi, double t, long ri, long ci, long im);


/// The function initializes the TRFEL from the HERMES code.
void hermtrf_init(int argc, const char *argv[]);

void hermtrf_close();
#endif
