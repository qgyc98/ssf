#ifndef ELEMENTT_H
#define ELEMENTT_H

#include <stdio.h>
#include "aliast.h"
#include "iotools.h"

struct matrix;
struct vector;

/**
   class elementt defines general element for transport problems
*/

class elementt
{
 public:
  elementt (void);
  ~elementt (void);
  void read (XFILE *in,long eid);
  void readmat (long m,XFILE *in);
  void print (FILE *out,long eid);
  void printmat (long m,FILE *out);

  void alloc_initnodval (long ndofe);
  void initnodvalues (vector &r);
  void subtrinitnodval (double *r,long ndofe);
 
  ///  type of element
  elemtypet te;

  ///  indicator of nodes with defined source of quantity
  ///  source = 0 - elements contains no node with source
  ///  source = 1 - elements contains nodes with source
  long source;
  ///  integration point pointer
  long **ipp;
  ///  type of cross section
  crsectypet crst;
  ///  number of appropriate cross section type
  long idcs;
  ///  transmission indicator
  ///  indicator of boundary condition on element
  ///  transi=0 - default value, no boundary conditions
  ///  transi=2 - prescribed fluxes (defined directly or due to climatic conditions)
  ///  transi=3 - prescribed transmission
  ///  transi=4 - element contains prescribed fluxes as well as prescribed transmission (e.g. there is one edge with Neumann condition and one edge with Newton condition)
  long *transi;

  ///  material type
  mattypet *tm;
  ///  material id
  long *idm;

  static long ntm;
  ///  array of initial nodal values
  ///  it is used in problems with changing number of elements
  double *initnodval;
  
  ///  transformation %matrix due to hanging nodes
  matrix *tmat;
};

#endif
