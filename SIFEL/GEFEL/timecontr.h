#ifndef TIMECONTR_H
#define TIMECONTR_H

#include <stdio.h>
#include "gfunct.h"
#include "iotools.h"

/**
   class timecontr controles time in problems with time
   
   JK
*/

class timecontr
{
 public:
  timecontr (void);
  ~timecontr (void);
  void read (XFILE *in);
  void print (FILE *out);
  double newtime ();
  double newtime (double &dt);
  void oldtime ();
  double starttime ();
  double endtime ();
  double actualtime ();
  double actualforwtimeincr ();
  double actualbacktimeincr ();
  double initialtimeincr ();
  void take_values (timecontr &tc);
  void seconds_days ();
  long isitimptime ();
  void initiate (timecontr &tc);
  void be_copy_of (timecontr &tc);
  void save_txt (FILE *out, int prec);
  void save_bin (FILE *out);
  void restore_txt (FILE *in);
  void restore_bin (FILE *in);

  //void reduce_step (double rf);
  
  ///  actual time
  double time;
  ///  actual backward time step
  double backwarddt;
  ///  actual forward time step
  double forwarddt;
  ///  backup of forward time step
  double bfdt;
  ///  backup of solver time step
  double bsdt;
  /// backup of user defined forward time step
  double bfdtu;
  ///  starting time
  double start_time;
  ///  end time
  double end_time;
  ///  function describing time increments
  gfunct timefun;
  ///  array of important times
  double *imptime;
  ///  number of important times
  long nit;
  ///  actual position in imptime array
  long apit;
  ///  is it important time now?
  ///  iiit=0 - no
  ///  iiit=1 - yes
  long iiit;
  ///  backup of iiit from previous time step
  long biiit;

  /** 
    flag for taking important times exactly, i.e. 
    exact_impt == true (default) -> time increment may be reduced in order the given important times can be attained exactly
    exact_impt == true, time increment will not be reduced with respect to given important times and the closest attained time 
                  step will be considered to be important one.
  */
  bool exact_impt;

  double it_step;

  /**  type of time controller
       tct = 0 - fixed - time increments are defined by user due the function timefun, solvers are not able to change these time steps
       tct = 1 - adaptive - there are minimum and maximum time steps prescribed by user, the minimum and maximum time steps are constants,
                            solver prescribes the time steps, the final time step is equal to the time step from solver except the case
                            where user time step differs for the actual time and the previous one
       tct = 2 - adaptivemin - the minimum time step is defined by a constant, smaller value of solver time step and user time is used
                               in the case of change of user time step,
       tct = 3 - adaptivemax
  */
  timecontrtype tct;
  ///  minimum time step
  double dtmin;
  ///  maximum time step
  double dtmax;

  /// actual minimum time step given by general time function 
  gfunct *dtminfun;

  /// actual maximum time step given by general time function 
  gfunct *dtmaxfun;
  
};

#endif
