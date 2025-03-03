#ifndef SLOPESOL_H
#define SLOPESOL_H

#include <stdio.h>

/// the arc-length method which reaches given value of the lambda parameter
double arclengthrv1 (long lcid, double rlambda, double rlerr);

/// the Newton-Raphson method which reaches given value of the lambda parameter
void newton_raphsonrv1 (long lcid, double rlambda, double rlerr);

#endif
