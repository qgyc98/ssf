#ifndef CPCSOLVER_H
#define CPCSOLVER_H

#include <stdio.h>

/**
   partially coupled problem
   transport problem influences mechanical problem
*/

void solve_gpcouplprob ();

void newton_raphson_gparcoupl (long lcid);
void newton_raphson_gparcoupl_lin (long lcid);
void newton_raphson_gparcoupl_nonlin (long lcid);

#endif
