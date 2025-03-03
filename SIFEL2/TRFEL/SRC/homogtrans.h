#ifndef HOMOGTRANS_H
#define HOMOGTRANS_H

#include <stdio.h>

void transport_homogenization_old_old();
void transport_homogenization();
void paral_transport_homogenization_old_old(double *u,double *k,double *c);
void paral_transport_homogenization(double *u,double *k,double *c);

#endif
