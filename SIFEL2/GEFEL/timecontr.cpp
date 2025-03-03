#include "timecontr.h"
#include <stdlib.h>
#include <math.h>

timecontr::timecontr (void)
{
  //  actual time
  time=0.0;
  //  time of start of analysis
  start_time=0.0;
  //  time of end of analysis
  end_time=0.0;

  //  forward time step
  forwarddt=0.0;
  //  backup of forward time step
  bfdt=0.0;
  //  backup of changed user forward time step
  bfdtu=0.0;
  //  backup of solver time step
  bsdt=0.0;
  //  backward time step
  backwarddt=0.0;
  
  //  number of important times
  nit=0;
  //  actual position in important times
  apit=0;
  //  array containing important times
  imptime=NULL;
  // no automatic generation of important times is required
  it_step = -1.0;
  //  indicator of important times
  iiit=0;
  // backup of iiit, i.e. iiit from previous time step
  biiit = 0;
  // flag for taking important times exactly
  exact_impt=true;
  
  //  type of time controller
  tct=notct;
  // actual minimum time step
  dtmin=0.0;
  // actual maximum time step
  dtmax=0.0;
  // actual minimum time step given by general time function 
  // (used in connection with adaptive_minmax)
  dtminfun = NULL;
  // actual maximum time step given by general time function 
  // (used in connection with adaptive_minmax)
  dtmaxfun = NULL;
  
}

timecontr::~timecontr (void)
{
  delete [] imptime;
  delete dtminfun; 
  delete dtmaxfun; 
}

/**
   function reads necessary information about time measurement
   
   @param in - input stream
   
   6.4.2003, JK
   19.9.2011, Modified by Tomas Koudelka (it_step)
*/
void timecontr::read (XFILE *in)
{
  long i;
  
  //  type of time controller
  //  fixed=0
  //  adaptive=1
  //  adaptivemin=2
  //  adaptivemax=3
  xfscanf (in,"%k%m","time_contr_type",&timecontrtype_kwdset,(int*)&tct);
  
  //  starting and ending times
  xfscanf (in,"%k%lf%k%lf","start_time",&start_time,"end_time",&end_time);
  
  //  number of important times
  xfscanf (in,"%k%ld","num_imp_times", &nit);
  
  //  important times
  // pokus JM 29.5.2008
  if (nit < 0){
    if (nit == -2)
      exact_impt = false;
    
    xfscanf (in,"%lf",&it_step);
    nit = long((end_time - start_time)/it_step);
    if (nit < 2)
      nit = 2;
    
    imptime=new double [nit];
    imptime[0] = start_time + it_step;
    
    for (i=1; i<nit-1; i++)
      imptime[i] = imptime[i-1] + it_step;
    imptime[nit-1] = end_time;
  }
  else{
    if (nit==0){
      imptime = new double [1];
      imptime[0]=end_time;
    }
    else{
      imptime=new double [nit];
    }
    for (i=0;i<nit;i++){
      xfscanf (in,"%lf",imptime+i);
    }
  }

  //  function prescribing time increments/steps
  timefun.read (in);

  switch (tct){
  case fixed:{
    dtmin = dtmax = initialtimeincr ();
    break;
  }
  case adaptive:{ // solver has main control over the time increments
    xfscanf (in,"%k%lf %k%lf","dtmin",&dtmin,"dtmax",&dtmax);
    break;
  }
  case adaptivemin:{
    xfscanf (in,"%k%lf","dtmin",&dtmin);
    dtmax = initialtimeincr (); // maximum time increment dtmax is given by timefun
    bfdtu = initialtimeincr ();
    break;
  }
  case adaptivemax:{
    xfscanf (in,"%k%lf","dtmax",&dtmax);
    dtmin = initialtimeincr (); // minimum time increment dtmax is given by timefun
    break;
  }
  case adaptive_minmax:{ // solver has main control over the time increments, limits are given by function
    dtminfun = new gfunct();
    dtmaxfun = new gfunct();
    xfscanf(in, "%k", "dtmin_fun");
    dtminfun->read(in);
    xfscanf(in, "%k", "dtmax_fun");
    dtmaxfun->read(in);
    bfdtu = initialtimeincr ();
    break;
  }
  default:{
    print_err("unknown type of time controller is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  //  initialization of actual time
  time = start_time;
  //  initialization of the first forward time increment
  forwarddt = initialtimeincr ();
  //  backup of forward time increment
  bfdt = forwarddt;
}



/**
   function prints necessary information about time measurement
   
   @param out - output stream
   
   6.4.2003, JK
*/
void timecontr::print (FILE *out)
{
  long i;
  
  fprintf (out,"\n%d\n", tct);
  //  starting and ending times
  fprintf (out,"%e %e\n",start_time,end_time);
  
  //  number of important times
  if (it_step >= 0.0)
    fprintf (out,"-1 %e\n", it_step);
  else
  {  
    fprintf (out,"%ld\n",nit);
    //  important times
    for (i=0;i<nit;i++){
      fprintf (out," %e \n",imptime[i]);
    }
  }
  
  //  function prescribing time increments/steps
  timefun.print (out);

  switch (tct)
  {
    case fixed:
      break;
    case adaptive:
      fprintf(out,"%e %e\n", dtmin, dtmax);
      break;
    case adaptivemin:
      fprintf (out, "%e\n", dtmin);
      break;
    case adaptivemax:
      fprintf (out, "%e\n", dtmax);
      break;
    case adaptive_minmax:
      dtminfun->print(out);
      fprintf(out, "\n");
      dtmaxfun->print(out);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of time controller is required",__FILE__,__LINE__,__func__);
  }
}



/**
   function returns new time with respect to important times
   
   6.4.2003, JK
   
   ???!!! candidate for removal
*/
double timecontr::newtime ()
{
  double newtime;
  
  //  trial time increment
  forwarddt=timefun.getval(time);
  
  //  trial new time
  newtime=time+forwarddt;
  
  //  indicator of important time
  iiit=0;
  if (apit<nit){
    if (imptime[apit]<=newtime){
      newtime=imptime[apit];
      iiit=1;
      apit++;
    }
    
    //  actual time step
    forwarddt=newtime-time;
    //  new time
    time=newtime;
  }
  else{
    //  new time
    time=newtime;
  }
  
  //  new backward time step
  backwarddt=forwarddt;
  return time;
}



/**
   function returns new time with respect to important times
   and time step enforced by an analysis
   
   @param dt - required time step from an analysis
   
   30.5.2008, JK
   31.10.2013 Rewritten by TKo+JK
*/
double timecontr::newtime (double &dt)
{
  double newtime, forwarddtu;

  biiit = iiit;
  if (exact_impt == false){
    // important times are required approximately
    if (apit<nit){
      // time difference between attained time and actual important time
      double ptitdt = fabs(time - imptime[apit]);
      // time difference between new time and actual important time
      double atitdt = fabs(time + dt - imptime[apit]);
      // time difference between next possible new time and actual important time
      double ntitdt = fabs(time + dt + dt - imptime[apit]);
      // test if the solver enlarged the time increment comparing to previous one
      if ((atitdt < ntitdt) && (ptitdt < atitdt) && (iiit == 0))
        // the new actual time should be the closest one -> do not enlarge time increment
        // if solver proposed so and use the original time increment
        if (dt > forwarddt)
          dt = forwarddt;
    }
  }

  
  //  trial time step prescribed by user
  forwarddtu=timefun.getval(time);
  
  if (tct == fixed)
  { // time increments are given by user

    if ((dt < forwarddtu) && (dt < bsdt))
    {
      // user has prescribed greater time increment than the one came from solver
      // and solver requires shorter time increment then the previous one
      // => wrong time increment precribed by user
      print_err("time increment %le required by solver is less\n"
                " than fixed time increment %le specified by user.\n"
                " Use adaptive time control instead of fixed one", 
                __FILE__, __LINE__, __func__, dt, forwarddtu);
      time=start_time;
      return -1.0;
    }

    forwarddt = forwarddtu;
    dtmin = dtmax = forwarddt;
  }
  else
  {// time increments has to be set with respect to dt prescribed by solver


    forwarddt=dt;
    if (iiit==1){
      //  previous time has been important time
      //  theoretically, the time step may be extremely small
      //  therefore, the time step before important time is used
      //  if solver did not require shorter time increment
      if (dt >= bsdt)
        forwarddt=bfdt;
    }


    if (tct == adaptivemin) 
    {
      // maximum time increment dtmax is given by timefun
      dtmax = forwarddtu;

      if ((bfdtu < forwarddtu) && (dt >= bsdt)) 
      {
        // user has prescribed greater time increment just in the given time step
        // and solver did not require shorter time increment in the last performed time step
        forwarddt = forwarddtu;
        bfdtu = forwarddtu;
      }

    }

    if (tct == adaptivemax) // minimum time increment dtmin is given by timefun
      dtmin = forwarddtu;
    
    if (tct == adaptive_minmax)
      {
	dtmin = dtminfun->getval(time);
	dtmax = dtmaxfun->getval(time);
      }
    
    //if(adaptive || adaptive_minmax)
    if(tct == adaptive_minmax)
      {
	if ((bfdtu != forwarddtu) && (dt >= bsdt)) 
	  {
	    // user has prescribed change in time increment just for the given time step
	    // and solver did not require shorter time increment in the performed time step
	    forwarddt = forwarddtu;
	    bfdtu = forwarddtu;
	  }
      }
  }
  
  if (forwarddt>dtmax)//maximum time increment
    forwarddt = dtmax;
  

  if (forwarddt<dtmin)
  {
    // there is prescribed greater time increment than the one came from solver
    // and solver requires shorter time increment then the previous one
    // => 
    print_err("time increment %.15le is less than minumum time increment %.15le.\n"
              "time increment proposed by solver is %.15le.\n",
              __FILE__, __LINE__, __func__, forwarddt, dtmin, dt);
    //abort();
    time=start_time;
    return -1.0;
  }
  

  //  trial new time
  newtime=time+forwarddt;

  // check new time not to be greater than end time
  if (newtime > end_time)
  {
    newtime=end_time;
    forwarddt = newtime - time;
  }
  
  // backup of forward time step determined with respect to solver,
  // user given time increments and end time
  bfdt=forwarddt;

  //  indicator of important time
  iiit=0;
  if (exact_impt){
    if (apit<nit){
      if (imptime[apit]<=newtime){
        //fprintf (stdout,"\n\n jsme v important time %ld\n\n",apit);
        newtime=imptime[apit];
        //  change actual time step
        forwarddt=newtime-time;
        iiit=1;
        apit++;
      }
    }
  }
  else{
    // important times will be taken approximately,
    // the closest attained time step to the actual important time will be considered as the important one
    // no reduction of time step will be performed
    if (apit<nit){
      // time difference between actual time and actual important time
      double atitdt = fabs(newtime - imptime[apit]);
      // time difference between next possible time and actual important time
      double ntitdt = fabs(newtime + forwarddt - imptime[apit]);
      if (atitdt <= ntitdt){
        // the actual time is the closest one probably
        iiit=1;
        apit++;
      }
    }
  }
  //  new time
  time=newtime;
  
  //  new backward time step
  backwarddt=forwarddt;

  // time increment is passed to solver in the dt argument
  dt=forwarddt;

  // backup of time increment passed to solver
  bsdt = dt;

  return time;
}



/**
   function shifts time to the previous time instant
   
   29.5.2008, JK
*/
void timecontr::oldtime ()
{
  if (exact_impt){
    if (apit>0){
      if (imptime[apit-1]==time){
        apit--;
      }
    }
  }
  else{
    if ((apit>0) && iiit) apit--;
  }
  
  //  back shift in time to the previous value
  time-=backwarddt;

  iiit = biiit;
  if (exact_impt){
    if (imptime[apit]==time)
      iiit=1;
    else
      iiit=0;
  }
}



/**
   function returns starting time
   
   6.4.2003, JK
*/
double timecontr::starttime ()
{
  return start_time;
}



/**
   function returns ending time
   
   6.4.2003, JK
*/
double timecontr::endtime ()
{
  return end_time;
}



/**
   function returns actual time
   
   6.4.2003, JK
*/
double timecontr::actualtime ()
{
  return time;
}



/**
   function returns actual forward time increment
   
   6.4.2003, JK
*/
double timecontr::actualforwtimeincr ()
{
  return forwarddt;
}



/**
   function returns actual backward time increment
   
   6.4.2003, JK
*/
double timecontr::actualbacktimeincr ()
{
  return backwarddt;
}



/**
   function returns initial time increment
   
   6.4.2003, JK
*/
double timecontr::initialtimeincr ()
{
  double dt;

  dt = timefun.getval(start_time);

  if (nit==0){
    return dt;
  }
  else{
    if (start_time+dt<imptime[0]){
      return dt;
    }
    else{
      dt = imptime[0]-start_time;
      return dt;
    }
  }
}



/**
   function takes values from another time controller
   
   @param timecontr - pointer to another time controller
   
   JK, 20.9.2004
*/
void timecontr::take_values (timecontr &tc)
{
  time = tc.time;
  backwarddt = tc.backwarddt;
  forwarddt = tc.forwarddt;
  start_time = tc.start_time;
  end_time = tc.end_time;
  iiit = tc.iiit;
  exact_impt = tc.exact_impt;
}



/**
   function transforms time in seconds to time in days
   
   JK, 20.9.2004
*/
void timecontr::seconds_days ()
{
  time = time/86400.0;
  backwarddt = backwarddt/86400.0;
  forwarddt = forwarddt/86400.0;
  start_time = start_time/86400.0;
  end_time = end_time/86400.0;
}



/**
   function returns indicator of important time
   
   iiit=0 - it is not important time
   iiit=1 - it is important time
*/
long timecontr::isitimptime ()
{
  return iiit;
}



/**
   function sets up input data from other object of the class timecontr
   
   @param tc - another time controller
   
   JK
*/
void timecontr::initiate (timecontr &tc)
{
  long i;
  
  start_time = tc.start_time;
  end_time = tc.end_time;
  
  timefun.initiate (tc.timefun);

  nit = tc.nit;
  exact_impt = tc.exact_impt;
  imptime = new double [nit];
  for (i=0;i<nit;i++){
    imptime[i]=tc.imptime[i];
  }
  
  time = start_time;

  //  initialization of the first forward time increment
  forwarddt = initialtimeincr ();
}



/**
   function reduces time step with respect to external conditions
   
   this function is used in nonlinear dynamics in cases of
   nonequilibriated stages
   
   @param rf - reduction factor
   
   in the case of reduction three times, rf = 3;
   
   10.9.2007, JK
*/
/*
void timecontr::reduce_step (double rf)
{
  //  time shift back
  time -= forwarddt;

  //  reduction of time step
  if (rf<1.0e-3){
    fprintf (stderr,"\n nonpositive reduction factor in function reduce_step (file %s, line %d).\n",__FILE__,__LINE__);
  }
  forwarddt /= rf;
  
  //  new time
  time += forwarddt;
}
*/



/**
   function copy values to another time controller
   
   @param timecontr - pointer to another time controller
   
   LS, 29.8.2012
*/
void timecontr::be_copy_of (timecontr &tc)
{
  time = tc.time;
  backwarddt  = tc.backwarddt;
  forwarddt  = tc.forwarddt;
  bfdt = tc.bfdt;
  bsdt = tc.bsdt;
  bfdtu = tc.bfdtu;
  apit = tc.apit;
  iiit = tc.iiit;
}



/**
  Function saves internal data to the given backup text file.

  @param out - pointer to the opened text backup file
  @param prec - precision for output of the real values

  @retval The function does not return anything

  Created 3.2009 by TKo
*/
void timecontr::save_txt (FILE *out, int prec)
{
  fprintf(out,"%.*le\n",prec,time);
  fprintf(out,"%.*le\n",prec,backwarddt);
  fprintf(out,"%.*le\n",prec,forwarddt);
  fprintf(out,"%.*le\n",prec,bfdt);
  fprintf(out,"%ld\n",apit);
  fprintf(out,"%ld\n",iiit);
}



/**
  Function saves internal data to the given backup binary file.

  @param out - pointer to the opened binary backup file

  @retval The function does not return anything

  Created 3.2009 by TKo
*/
void timecontr::save_bin (FILE *out)
{
  fwrite(&time, sizeof(time), 1, out);
  fwrite(&backwarddt,  sizeof(backwarddt),  1, out);
  fwrite(&forwarddt,  sizeof(forwarddt),  1, out);
  fwrite(&bfdt, sizeof(bfdt), 1, out);
  fwrite(&apit, sizeof(apit), 1, out);
  fwrite(&iiit, sizeof(iiit), 1, out);
}



/**
  Function restores internal data from the given backup text file.

  @param out - pointer to the opened text backup file

  @retval The function does not return anything

  Created 3.2009 by TKo
*/
void timecontr::restore_txt (FILE *in)
{
  long i;

  fscanf(in,"%le",&time);
  fscanf(in,"%le",&backwarddt);
  fscanf(in,"%le",&forwarddt);
  fscanf(in,"%le",&bfdt);
  fscanf(in,"%ld",&apit);
  fscanf(in,"%ld",&iiit);
  for (i=0, apit=0; i<nit; i++, apit++)
  {
    if (imptime[i] > time)
      break;
  }
}



/**
  Function restores internal data from the given backup binary file.

  @param out - pointer to the opened binary backup file

  @retval The function does not return anything

  Created 3.2009 by TKo
*/
void timecontr::restore_bin (FILE *in)
{
  long i;
  fread(&time, sizeof(time), 1, in);
  fread(&backwarddt,  sizeof(backwarddt),  1, in);
  fread(&forwarddt,  sizeof(forwarddt),  1, in);
  fread(&bfdt, sizeof(bfdt), 1, in);
  fread(&apit, sizeof(apit), 1, in);
  fread(&iiit, sizeof(iiit), 1, in);
  for (i=0, apit=0; i<nit; i++, apit++)
  {
    if (imptime[i] > time)
      break;
  }
}



