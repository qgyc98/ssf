#ifndef STOCHDRIVERT_H
#define STOCHDRIVERT_H

#include <stdio.h>
#include "aliast.h"
struct vector;
struct atsel;
#include  "matrix.h"
#include  "iotools.h"

/**

*/
class stochdrivert
{
 public:
  stochdrivert (void);
  ~stochdrivert (void);
  
  void compute_nprunknowns ();
  
  void read (XFILE* in);
  void readtable (XFILE *in);
  void writetable ();
  
  void changevalues (long sampleid);
  void changematerials (long id,vector &val);
  void changecrsections (long id,vector &val);
  
  void extractor (long sampleid);

  void importvalues (double *dat);
  void exportvalues (double *val);
  
  //  name of auxiliary file
  char auxfile[1001];
  //  number of stochastic materials
  long nsmt;
  //  number of stochastic cross-sections
  long nscs;
  
  mattypet *mt;
  long *idm;
  atsel *atm;

  crsectypet *cst;
  long *idcs;
  atsel *atcs;

  //  stochastic loaded nodes
  long *idln;
  atsel *atln;
  
  long nsampl;
  long nstochvar;
  matrix stochtabin;
  long nprunknowns;
  matrix stochtabout;
  
  //  number of printed nodal displacements
  long npnv;
  //  array conatining numbers of nodes
  long *nna;
  //  numbers of particular DOFs
  atsel *pnv;
  
  
};

#endif
