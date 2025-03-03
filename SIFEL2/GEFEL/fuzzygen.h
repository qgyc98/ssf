#ifndef FUZZYGEN_H
#define FUZZYGEN_H

#include "iotools.h"
#include "fuzzynum.h"


/**
   class fuzzygen
   
   
   JK, 18.8.2005
*/
class fuzzygen
{
 public:
  fuzzygen (void);
  ~fuzzygen (void);
  
  void read (XFILE *in);
  void totcombnumber ();
  void combnumber ();

  void gener_alphacuts (long alphid,double *avi);
  void gener_allcomb (double *avi);

  void give_new_values (long sampleid,double *avi,FILE *out);
  void save_values (fuzzynum *fnum,long nprunknowns,double *avo);
  
  
  ///  number of fuzzy variables
  long nfv;
  ///  number of alpha-cuts
  long nalph;
  ///  number of values
  long nval;

  ///  number of combinations in one alpha-cut
  long ncomb;
  ///  total number of combinations
  long tncomb;
  
  ///  number of actual alpha-cut
  long actalph;
  ///  number of actual combination in actual alpha-cut
  long actcomb;

  
  ///  array of required alpha-cuts
  double *alpha;

  long *aux;
  
  ///  fuzzy numbers
  fuzzynum *fn;

};

#endif
