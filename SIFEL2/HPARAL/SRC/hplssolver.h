#ifndef HPLSSOLVER_H
#define HPLSSOLVER_H

//#include <stdio.h>

void parallel_homogenization_linear_statics_tiles ();
void parallel_homogenization_linear_statics ();
void higher_to_lower_levelm (long llid,long ncbuff,double *buff);

#endif
