#ifndef ELEMENT_H
#define ELEMENT_H

#include <stdio.h>
#include "alias.h"
#include "iotools.h"
//#include "vector.h"
#include "matrix.h"


/**
  Class element defines general finite element.
   
  Created by JK,
*/

class element
{
 public:
  element (void);
  ~element (void);
  void read (XFILE *in,long eid);
  void readmat (XFILE *in);
  void print (FILE *out,long eid);
  void printmat (FILE *out);

  void alloc_initdispl (long ndofe);
  void initdisplacement (double *r,long ndofe);
  void subtrinitdispl (double *r,long ndofe);
  void alloc_growstr (long eid);

  ///  type of element
  elemtype te;
  ///  indicator of prescribed displacements on element
  long prescdispl;
  ///  indicator of temperature changes on element
  long presctemp;
  ///  computation of reactions
  long react;
  ///  array of integration point pointers
  long **ipp;
  ///  type of cross section
  crsectype crst;  
  ///  number of appropriate cross section type
  long idcs;
  ///  stress/strain state
  strastrestate ssst;
  ///  number of strain/stress components
  long ncomp;
  
  ///  array of initial displacements
  ///  it is used in problems with changing number of elements
  double *initdispl;
  
  ///  transformation %matrix due to hanging nodes
  matrix *tmat;

  /*
  priprava pro kodova cisla na prvku
  
  number of element only code numbers
  long necn;
  array containing element only code numbers
  long *cne;
  */

  ///  type of material
  mattype *tm;
  ///  number of appropriate material type
  long *idm;
  ///  number of material types
  long nm;
  
};

#endif
