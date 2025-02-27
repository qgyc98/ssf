#ifndef HPNPSOLVERT_H
#define HPNPSOLVERT_H

//#include <stdio.h>

void parallel_homogenization_lin_nonstat ();

void parallel_homogenization_new ();

void parallel_homogenization ();

void higher_to_lower_level (long llid,long ncbuff,double *buff);

#endif
