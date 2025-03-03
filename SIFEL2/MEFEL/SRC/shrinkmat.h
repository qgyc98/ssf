#ifndef SHRINKMAT_H
#define SHRINKMAT_H

#include <stdio.h>
#include "alias.h"
#include "iotools.h"
#include "gfunct.h"
struct matrix;
struct vector;

/**
  This class defines material model which computes
  irreversible shrinkage strains
  
  21. 11. 2013
*/
class shrinkmat
{
 public:
  shrinkmat (void);
  ~shrinkmat (void);

  void read (XFILE *in);
  void print (FILE *out);
  
  void matstiff (matrix &d,long ipp,long im,long ido);
  void nlstressesincr (long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  void initvalues (long lcid, long ipp, long im, long ido, bool rinit);
  void givestressincr (long lcid, long ipp, long im, long ido, long fi, vector &sig);
  void giveirrstrains(long ipp, long im, long ido, vector &epsirr);
  double give_actual_ft (long ipp, long im, long ido);
  double give_actual_fc (long ipp, long im, long ido);  

  void give_reqnmq(long *anmq);
  
  tshrlaw tshr;  // type of shrinkage law
  gfunct beta;   // general function with derivatives of shrinkage with respect to humidity
  gfunct shmeas; // general function with measured curve of shrinkage in dependence on humidity
  answertype thumdef; // flag of relative humidity definition (1/yes = defined in the model, 0/no = taken from nonmechq)
  gfunct humdef; // time dependent course of relative humidity defined in the material in the case of thumdef=1/yes
};

#endif
