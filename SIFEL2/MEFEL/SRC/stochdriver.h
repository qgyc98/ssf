#ifndef STOCHDRIVER_H
#define STOCHDRIVER_H

#include "xfile.h"
#include "alias.h"
#include "fuzzygen.h"
#include "fuzzynum.h"
#include <stdio.h>
#include "matrix.h"

struct vector;

/**
   class stochdriver
   
   this class controles stochastic computations
   
   JK, TK
*/

class stochdriver
{
 public:
  stochdriver ();
  ~stochdriver ();
  
  void compute_nprunknowns ();
  
  void read (XFILE* in);
  void readtable (XFILE *in);
  void writetable ();
  
  void changevalues (long sampleid);
  void assemble_new_values (long sampleid);
  void replace_values ();
  
  void changematerials (long id,vector &val);
  void changecrsections (long id,vector &val);
  void changenodloads (long id,vector &val);
  
  void extractor ();
  void save_results (long sampleid);
  
  void diagpostproc ();
  
  void update_auxparam ();
  
  /*
  void give_new_invalues (double *buff);
  void save_new_invalues (double *buff);
  void give_new_outvalues (double *buff);
  void save_new_outvalues (double *buff);
  */

  ///  name of auxiliary file
  char auxfilein[1001];
  char auxfileout[1001];

  ///  number of stochastic materials
  long nsmt;
  ///  number of stochastic cross-sections
  long nscs;
  ///  number of stochastic nodal loads
  long nsnl;
  
  ///  data about material types
  mattype *mt;
  long *idm;
  atsel *atm;
  
  ///  data about cross-section types
  crsectype *cst;
  long *idcs;
  atsel *atcs;
  
  ///  stochastic loaded nodes
  long *idln;
  atsel *atln;
  
  ///  number of samples
  long nsampl;
  ///  number of stochastic variables
  long nstochvar;
  ///  number of printed output variables
  long nprunknowns;

  matrix stochtabin;
  matrix stochtabout;
  
  
  ///  number of printed nodal displacements
  long npnd;
  ///  array conatining numbers of nodes
  long *nna;
  ///  numbers of particular DOFs
  atsel *pnd;

  ///  number of required eigenvectors
  long neigv;
  
  ///  number of elements with printed values
  long npev;
  ///  array containing numbers of elements
  long *ena;
  ///  description of printed values on elements
  atsel *ev;
  
  
  long ndispl,nelem;





  
  ///  array of actual input variables
  double *avi;
  ///  array of actual output variables
  double *avo;
  ///  input stream; contains particular samples
  XFILE *datin;
  ///  output stream; contains particular samples
  FILE *datout;
  
  
  ///  generator of fuzzy numbers
  fuzzygen fg;
  ///  output fuzzy numbers
  fuzzynum *fn;
};

#endif
