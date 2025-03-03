#ifndef FDSOLVER_H
#define FDSOLVER_H

#include <stdio.h>

void solve_forced_dynamics ();


void newmark_method (long lcid);
void optim_newmark_method (long lcid);

void difference_method (long lcid);

void verlet_method (long lcid);

void response_spectrum_method (long lcid);

void explicit_difference_method (long lcid);

void difference_method_2 (long lcid);

#endif
