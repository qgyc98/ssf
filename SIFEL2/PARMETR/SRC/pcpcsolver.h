#ifndef PCPCSOLVER_H
#define PCPCSOLVER_H

#include <stdio.h>

void par_solve_gpcouplprob ();

void par_newton_raphson_gparcoupl (long lcid);

void par_newton_raphson_gparcoupl_lin (long lcid);
void par_newton_raphson_gparcoupl_nonlin (long lcid);


#endif
