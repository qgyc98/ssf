#ifndef TIMESWMAT_H
#define TIMESWMAT_H

#include "iotools.h"
#include "alias.h"
#include "galias.h"
struct matrix;
struct vector;
class  gfunct;

/**
  This class defines artificial material model which combines
  several arbitrary material models for time dependent problems
  These models are switched on and off depending on their time functions.
  The time functions must be defined for every material model in 
  the model chain. The models in the model chain must be ascending ordered by 
  the time of activation and organized into the logical groups comprising
  standard definition of the MEFEL material model for the required time
  (i.e. the control material model must be on the frist position in the group,
  the elastic model should be the last and so on).
  Only those time functions, which are related with the concatenated active material models, 
  should be defined to switch on the given model at required time while the 
  remaining models in the model chain must be switched off.

  There is also switch for application of stress free state for newly switched on 
  material models. If the flag store_eigstr is set to 'yes' then actual strains
  will be stored and subtracted before stress calculation.
  
  15.1.2004
  Updated by Tomas Koudelka, 8.2.2016
*/
class timeswmat
{
 public:
  timeswmat (void);
  ~timeswmat (void);
  void read(XFILE* in);
  void print(FILE* out);
  double give_actual_ym(long ipp, long im, long ido);
  double give_initial_ym(long ipp, long im, long ido);
  double give_actual_nu(long ipp, long im, long ido);
  double give_actual_ft(long ipp, long im, long ido);
  //double give_pore_press(long ipp, long im, long ido);
  long givencompother (long ipp, long im);
  long givencompeqother (long ipp, long im);
  void initvalues(long lcid, long ipp, long im, long ido, bool rinit);
  void matstiff (matrix &d,long ipp, long im, long ido);
  void givestressincr (long lcid,long ipp,long im,long ido,long fi,vector &sig);
  void nlstressesincr (long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  long update_ami(long im);

  /// number of combined material models
  long ncm;
  /// flag for storage of eigenstrains due to stress free state of new activated material models
  answertype store_eigstr;
  ///  time functions for material switching
  gfunct *gf;
  /// index of the first active material, ami=0 means the first material AFTER the timeswmat
  long ami;
  /// number of active materials
  long nam;
  
};

#endif
