#ifndef POINTSET_H
#define POINTSET_H

#include <stdio.h>

/**
  This class groups data about point setson the elements, where the user requires
  to compute stresses or strains. It is used for the mechprep preprocessor
*/
class pointset
{
  public:
	 pointset();
	 ~pointset();
	 long read(FILE *in);
	 long print(FILE *out);

	 long   n;      ///< number	of set
	 long   np;     ///< number of point
	 long   ncoord; ///< number of point coordinate
	 double *coord; ///< array with coordinates of each point
};

#endif
