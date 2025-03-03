/** 
    Continuos retardation spectrum for solidification theory of concrete creep for B3 model
    according to Bazant and Prasannan and V. Smilauer
*/
#include "creep_rspec.h"
#include "iotools.h"
#include "galias.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "globmat.h"
#include "intpoints.h"
#include "elastisomat.h"
#include "creep.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>



rspecmat::rspecmat (void)
{
  type_h=0;
  type_temp=0;

  type_e=0;

  hum_env = temp_env = 0.0;

  flag_drshr = flag_shr = flag_temp = 0;

  ft_flag = 0;

  //effect of temperature on creep
  kappa = 0.0;
  beta_t = beta_h = beta_s = 1.0;

  //basic creep parameters
  q1_ini = 0.20179785365853656;
  q2_ini = 1.2818023122513782;
  q3_ini = 0.0392024641578089227;
  q4_ini = 0.129358593798123332;
  q5_ini = 0;

  q1 = 0.20179785365853656;
  q2 = 1.2818023122513782;
  q3 = 0.0392024641578089227;
  q4 = 0.129358593798123332;
  q5 = 0;
  
  C_const = 0.0;
  
  n_param = 0.1;
  m_param = 0.5;
  
  min_ret_time = 1.e-2;//minimum retardation time used in calculation [days], truncated values will be hidden in C_{const}, should be about 1.0e-16 for exact q1 determination, otherwise can be like 1.e-3
  
  alpha = 12.0e-6; //K^-1
  e28=30.0e9; //Pa
  fc=35.8e6; //Pa
  ft_flag = 0;
  ft=1.5e6; //Pa
  ft_ratio = 1.0e-5;
  wc=0.43;
  sc=3.4;
  gc=1.98;
  cs=305.0;  //kg*m^-3
  a1=1.05;
  a2=1.2;
  kd=0.15;   //m
  ks=1.0;
  type_b3 = 1;
  type_e=0;
  type_h=0;
  type_temp=0;

  tb_time = 0.0;
  th_time = 0.0;
  napproxtime = 0;
  nRetTime = 7;
  type_rt = 0;
  e0 = 0.0;
  previoustime = 0.0;
  actualtime = 0.0;
  dtb = 0.0;
  tb_age_dt = tb_age = tbl_age = tbh_age = maxtimeb = 0.0;

  retTime=NULL;
  timeMax=(Mp->timecon.endtime ())/86400.0;//in days

  eps_ainf = 0.0;
}

rspecmat::~rspecmat (void)
{
  delete [] retTime;
}

/**
   function reads material parameters

   @param in - input file

   TKr, 3.6.2005
*/

void rspecmat::read (XFILE *in)
{
  long i;

  xfscanf(in,"%ld",&type_b3);
  //////////////////////////////////////////////
  if (type_b3 == 2){//reading of material parameters of B3 model
    xfscanf (in,"%lf %lf %lf %lf %lf", &q1_ini, &q2_ini, &q3_ini, &q4_ini, &q5_ini);
  }
  //////////////////////////////////////////////
  else{//reading of material parameters of B3 model for mixture composition
    xfscanf(in,"%ld",&type_e);
    
    if (type_e == 1){xfscanf (in,"%lf ", &e28);}
    //input in Pa
    e0 = 1.5*e28/6.89476*1.0e-3;//to psi
    e28 = e28/6.89476*1.0e-3;//to psi
    
    xfscanf (in,"%lf %ld", &fc, &ft_flag);
    
    switch(ft_flag){
    case 0:
      xfscanf (in,"%lf ", &ft);
      break;
    case 1:
      xfscanf (in,"%lf ", &ft_ratio);
      break;
    case 2:{//temperature effect on tensile streght
      xfscanf (in,"%lf %lf ", &ft, &ft_ratio);
      break;
    }
    default:{
      print_err("\n unknown type ft_flag is required \n",__FILE__, __LINE__, __func__);
    }
    }
    
    xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &alpha, &wc, &sc, &gc, &cs, &a1, &a2, &ks, &kd);
    
    //input in Pa
    fc = fc/6.89476*1.0e-3;//to psi
    //input in kg*m^-3
    cs = cs/16.03;//to lb*ft^-3
    //input in m
    kd=kd/0.0254;//to inch
  }
  
  ////////////////////////////////////////////
  //reading of temperature and humidity effect
  xfscanf (in,"%ld %ld", &type_h, &type_temp);  
  
  //environmental humidity and temperature
  if (type_h == 0){xfscanf (in,"%lf ", &hum_env); }
  if (type_temp == 0){xfscanf (in,"%lf ", &temp_env); }
  
  xfscanf (in,"%lf %lf %ld %ld %ld %d %d %d", &tb_time, &th_time, &napproxtime, &nRetTime, &type_rt, &flag_drshr, &flag_shr, &flag_temp);
  
  tb_time = tb_time/86400.0;//days
  th_time = th_time/86400.0;//days
  
  retTime = new double [nRetTime];//vector of ret. times
  for (i=0;i<nRetTime;i++){
    retTime[i]=0.0;
  }
  
  if (type_rt==1){//reading of retardation times (ages) [days]
    for (i=0;i<nRetTime;i++)
      xfscanf (in,"%lf",&retTime[i]);
  }
  
  xfscanf (in,"%lf",&eps_ainf);
  xfscanf (in,"%lf",&kappa);
  
}





/**
   function prints material parameters

   @param out - output file

   TKr, 28/8/2008 - revised//opravit??!!
*/
void rspecmat::print (FILE *out)
{
  long i;
  
  fprintf (out," %ld",type_b3);
  //////////////////////////////////////////////
  if (type_b3 == 2){//printing of material parameters of B3 model
    fprintf (out," %lf %lf %lf %lf %lf", q1_ini, q2_ini, q3_ini, q4_ini, q5_ini);
  }
  //////////////////////////////////////////////
  else{//reading of material parameters of B3 model for mixture composition
    fprintf(out," %ld",type_e);
    
    if (type_e == 1){
      e28 = e28*6.89476/1.0e-3;//to Pa
      
      fprintf (out," %lf",e28);
    }
    
    fc = fc*6.89476/1.0e-3;//to Pa
    
    fprintf (out," %lf %d", fc, ft_flag);
    
    if (ft_flag == 0)
      fprintf (out," %lf", ft);
    else
      fprintf (out," %lf", ft_ratio);
    
    cs = cs*16.03;//to kg*m^-3
    kd=kd*0.0254;//to m
    
    fprintf (out," %lf %lf %lf %lf %lf %lf %lf %lf %lf", alpha, wc, sc, gc, cs, a1, a2, ks, kd);
  }
  
  fprintf (out," %ld %ld", type_h, type_temp);
  
  if (type_h == 0){fprintf (out," %lf ", hum_env); }
  if (type_temp == 0){fprintf (out," %lf ", temp_env); }
  
  tb_time = tb_time*86400.0;//to seconds
  th_time = th_time*86400.0;//to seconds
  
  fprintf (out," %lf %lf %ld %ld %ld %d %d %d", tb_time, th_time, napproxtime, nRetTime, type_rt, flag_drshr, flag_shr, flag_temp);
  
  if (type_rt==1){
    for (i=0;i<nRetTime;i++)
      fprintf (out," %lf", retTime[i]);
  }
  
  fprintf (out," %lf",eps_ainf);
  fprintf (out," %lf",kappa);
}

/**
   function computes ages in days

   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2008 - revised
*/
void rspecmat::compute_ages (long ipp, long ido)
{
  //times and ages in days
  if(Mp->time == previoustime){
    //    dtb = Mp->timecon.initialtimeincr ()/86400.0;
    dtb = Mp->timecon.actualforwtimeincr()/86400.0;
  }
  else{
    dtb = (Mp->time - previoustime)/86400.0;
  }
  actualtime = Mp->time/86400.0;
  tbh_age = actualtime - th_time;
  tbl_age = actualtime - tb_time + dtb/2.0;
  tb_age = actualtime - tb_time;
  napptimeb = give_napproxtime(ipp);
  maxtimeb = Mp->timecon.endtime ()/86400.0;//days
  tb_age_dt = tb_age + dtb;

  if (tb_age <= 0.0){
    print_err("\n age of concrete must be greater than zero! (ipp=%ld, eid=%ld)"
              "\n Age of concrete = %e, Actual time = %e, Time of end of concrete casting = %e;"
              "\n (Age of concrete = Actual time - Time of end of concrete casting)",
              __FILE__, __LINE__, __func__, ipp, Mm->elip[ipp]+1, tb_age, actualtime, tb_time);
    exit(0);
  }
  if (tbh_age <= 0.0){
    fprintf (stderr,"\n Age of concrete when drying begins must be greater than zero!!!");
    fprintf (stderr,"\n Age of concrete when drying begins = %e, Actual time = %e, Time when drying begins = %e; (Age of concrete when drying begins = Actual time - Time when drying begins).",tbh_age,actualtime,th_time);
    print_err("\n In function compute_ages ", __FILE__, __LINE__, __func__);
    exit(1);
  }

  if (th_time <= tb_time){
    fprintf (stderr,"\n Time when drying begins must be greater than Time of end of concrete casting!!!");
    fprintf (stderr,"\n Time when drying begins = %e   Time of end of concrete casting = %e.",th_time,tb_time);
    print_err("\n In function compute_ages ", __FILE__, __LINE__, __func__);
    exit(2);
  }

  //toto opravit
  beta_t = give_beta_t_eqother(ipp,ido);

  dtb = beta_t*beta_h*beta_s*dtb;

  //fprintf (Out,"\n actualtime = %lf",actualtime);
  //fprintf (Out,"\n tbh_age =  %lf",tbh_age);
  //fprintf (Out,"\n tbl_age =  %lf",tbl_age);
  //fprintf (Out,"\n tb_age =  %lf",tb_age);
  //fprintf (Out,"\n napptimeb =  %lf",napptimeb);
  //fprintf (Out,"\n maxtimeb =  %lf",maxtimeb);
  //fprintf (Out,"\n tb_age_dt =  %lf",tb_age_dt);
  //fflush(Out);
}


/**
   function returns computed ages

   @param t_age_dt -    
   @param t_age    -    
   @param tl_age   -    
   @param th_age   -    
   @param dt       -    
   @param maxtime  -    
   @param napptime -    
   @param ipp      - number of integration point

   TKr, 28/8/2008 - revised
*/
void rspecmat::give_ages (double &t_age_dt,double &t_age,double &tl_age,double &th_age,double &dt,double &maxtime,long &napptime,long /*ipp*/)
{
  t_age_dt = tb_age_dt;
  t_age = tb_age;
  tl_age = tbl_age;
  th_age = tbh_age;
  dt = dtb;
  maxtime = maxtimeb;
  napptime = napptimeb;
}


/**
   function returns number of components of eqother array

   @param ipp - number of integration point

   TKr, 28/8/2008 - revised
*/
long rspecmat::give_nceqother (long ipp)
{
  long nc,nceqother;

  nc=Mm->ip[ipp].ncompstr;

  nceqother = nc + nc;                  //total strains, total stresses
  nceqother = nceqother + nc + nc;      //strain increments, stress increments
  nceqother = nceqother + 2;            //total humidity, total temperature
  nceqother = nceqother + nRetTime;     //coefficients of Dirichlet series
  nceqother = nceqother + nc*nRetTime;  //hidden strains increments
  nceqother = nceqother + nc;           //creep-irreversible strains increments
  nceqother = nceqother + nc;           //shrinkage-irreversible strains increments
  nceqother = nceqother + nc;           //stress-induced shrinkage-irreversible strains increments
  nceqother = nceqother + nc;           //total irreversible strains (creep-irreversible strains + shrinkage-irreversible strains + stress-induced shrinkage-irreversible strains)
  nceqother = nceqother + 1;            //previous free shrinkage
  nceqother = nceqother + 1;            //actual Young's modulus
  nceqother = nceqother + 1;            //actual beta_t
  nceqother = nceqother + nc;           //total "aeging" strains = creep-irreversible strains + shrinkage-irreversible strains + stress-induced shrinkage-irreversible strains + strains due to aegin of concrete - changing effective Young's modulus
  nceqother = nceqother + 1;            //actual ft
  nceqother = nceqother + 1;            //old Young's modulus from previous time step
  
  return(nceqother);
}


/**
   function returns number of components of other array

   @param ipp - number of integration point

   TKr, 07/08/2008, not used now
*/

long rspecmat::give_ncother (long ipp)
{
  long nc,ncother;

  nc=Mm->ip[ipp].ncompstr;

  //ncother = nc;
  //temporarily
  ncother = give_nceqother (ipp);

  return(ncother);
}


/**
   Function initializes eqother array with initial values of temperature.
   Actual values of quantity 'initial_temperature' from array Mm->nonmechq are
   taken as initial temperature.
   
   @param ipp - integration point pointer
   @param ido - index of internal variables for given material in the ipp other array
   
   TKo+TKr, 28/8/2008 - revised
*/
void rspecmat::initvalues (long ipp, long /*im*/, long ido)
{
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  //if(Mm->moist != NULL)
  //Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)]=Mm->initmoist[ipp];

  if(Mm->givestatusnmq(initial_temperature) == 1)
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+1)]=Mm->givenonmechq(initial_temperature, ipp);//initial temperature

  //storing of beta_t
  store_beta_t_eqother(1.0,ipp,ido);
  //storing of actual ft
  //ft = compute_actual_ft (ipp,im,ido);//toto promyslet
  store_ft_eqother(ipp,ido);
}


/**
   
   @param ipp - integration point pointer
   @param ido - index of internal variables for given material in the ipp eqother array
   
   TKo+TKr, 28/8/2008 - revised
*/
void rspecmat::updatevalues (long ipp, long im,long ido)
{
  double beta;

  //updating of beta_t
  //beta_t computation
  beta = compute_beta_t(ipp);
  //storing of beta_t
  store_beta_t_eqother(beta,ipp,ido);
  //storing of actual ft
  ft = compute_actual_ft (ipp,im,ido);
  store_ft_eqother(ipp,ido);
}


/**
   function returns number of retardation times

   TKr, 3/6/2005
*/
long rspecmat::give_nret_time (void)
{
  return(nRetTime);
}


/**
   function computes ... 

   @param q3  - 
   @param q   - 
   @param tau - 

   TKr, 10/9/2008 - revised
*/
double rspecmat::give_L(double qq3, double q, double tau){

  double L;

  L = -2.*n_param*n_param*pow(3.*tau, 2.*n_param-3.)*(n_param-1.-pow(3.*tau,n_param)) / pow((1.+pow(3.*tau,n_param)),3.);
  L += (n_param*(n_param-2.)*pow(3*tau, n_param-3.)*(n_param-1.-pow(3.*tau,n_param))-n_param*n_param*pow(3.*tau,2.*n_param-3.)) / ((1.+pow(3.*tau,n_param))*(1.+pow(3.*tau,n_param)));
  L *= qq3*pow(3.*tau,3.)/2.;
  L *= log (10.) * q;//now L is A_mu = 1/E_mu

  //   fprintf(Out,"\ntau = %e,   L = %e",tau,L);//debug??!!
  
  return L;
}


/**
   function computes creep compliacne function J(t,t') [1/Pa] in a material point B3 Bazant's model
   dependent on temerature and humidity changes

   @param t0        - age when drying begins [days]
   @param tl        - age at loading [days]
   @param t         - time representing age of concrete [days]
   @param ipp       - number of integration point
  
   TKr, 4/9/2006
   nonlinear creep effect will be added soon
*/
double rspecmat::give_J_E_mu(vector &e_mu,double t0, double tl, double t, long ipp,long ido){
  
  double jt,L,ret_time_hlp;
  double n,cd,ac,ag,m,z,r,qf,q,eps_shinf,tau,ht,htl,st,stl;
  double hum,temp,hum_prev,temp_prev;
  double kt;
  double aux;
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  if ( type_h == 0){
    hum = 0.0;
    hum_prev = 0.0;
  }
  else{
    // actual humidity
    hum=Mm->givenonmechq(rel_hum, ipp);
    // humidity from previous time step
    hum_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)];
  }
  if ( type_temp == 0){
    temp = 0.0;
    temp_prev = 0.0;
  }
  else{
    // actual temperature
    temp=Mm->givenonmechq(temperature, ipp);
    // temperature from previous time step
    temp_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+1];
  }

  if (t < t0){
    fprintf (stderr,"\n Age of concrete is lower than age when dyring begins!!!");
    fprintf (stderr,"\n Age of concrete = %e,   Age when dyring begins = %e",t,t0);
    print_err("\n In ", __FILE__, __LINE__, __func__);
    exit(1);
  }
  if (t < tl){
    fprintf (stderr,"\n Age of concrete is lower than age at loading!!!");
    fprintf (stderr,"\n Age of concrete = %e,   Age at loading = %e",t,tl);
    print_err("\n In ", __FILE__, __LINE__, __func__);
    exit(0);
  }

  switch(type_b3){
  case 1:
    //basic creep    
    ac=sc+gc;   
    ag=ac/gc;   
    m = m_param; //m=0.28+1.0/fc/fc;
    n = n_param;

    if (type_e == 1)
      q1=600000.0/e28;//measured
    else
      q1=600000.0/57000.0/sqrt(fc);//empirical
    
    r=1.7*pow(tl,0.12)+8.0;
    qf=1.0/(0.086*pow(tl,(2.0/9.0))+1.21*pow(tl,(4.0/9.0)));
    z = 1.0/pow(tl,m)*log(1.0 + pow((t - tl),n));
    q = qf/pow((1.0 + pow((qf/z),r)),(1.0/r));
    q2=451.1*sqrt(cs)/pow(fc,0.9);
    q3=0.29*pow(wc,4.0)*q2;
    q4=0.14/pow(ac,0.7);

    
    //drying creep
    cd = 0.0;
    //humidity effect
    if (type_h == 0){
      //constant huminidy
      //kt = 190.8/pow(t0,0.08)*fc;//this is wrong this is from literature (t0 = age when drying begins)
      kt = 190.8/pow(t,0.08)/pow(fc,0.25);//this is according to creepb.cpp (Fajman)
      // according to SI units
      // kt = 0.085/pow(t,0.08)/pow(fc,0.25);//Bazant& Jirasek red book
      tau = kt*ks*ks*kd*kd;
      st = tanh(sqrt((t - t0)/tau));
      stl = tanh(sqrt((tl - t0)/tau));
      ht = 1.0 - (1.0 - hum_env)*st;
      htl = 1.0 - (1.0 - hum_env)*stl;
      eps_shinf = a1*a2*(26.0*pow((wc*cs),2.1)/pow(fc,0.28)+270.);
      // according to SI units
      //Bazant& Jirasek red book
      // eps_shinf = a1*a2*(0.019*pow((wc*cs),2.1)/pow(fc,0.28)+270.);
      q5 = 7.57e5/fc/pow(eps_shinf,0.6);
      cd = q5*sqrt(exp(-8.0*ht) - exp(-8.0*htl)); 
      
      if (tl < t0)
	cd = 0.0;
    }
    
    
    /*  
	if (type_h == 1){
	//changing humidity
	//stress-induced effect is included as strain in function: give_deps_stressinduced, this is old:
	dhum = hum-hum_prev; 
	eps_shinf = a1*a2*(26.0*pow((wc*cs),2.1)/pow(fc,0.28)+270.);
	q5 = 0.35/fc;
	cd = fabs(q5*eps_shinf*3.0*hum*hum*dhum);
	}
	
	//temperature effect
	if ( type_temp == 0 ){
	}
	else {
	//changing temperature
	dtemp = temp - temp_prev;   
	q5= 1.5/fc;
	cd = cd+fabs(q5*alpha*dtemp);
	}
    */
    
    break;
  case 2:
    q1 = q1_ini;
    q2 = q2_ini;
    q3 = q3_ini;
    q4 = q4_ini;
    q5 = q5_ini;
    break;
  default:{
    print_err("\n unknown type B3 concrete is required \n",__FILE__, __LINE__, __func__);
  }
  }
  
  if (flag_drshr == 0)
    cd = 0.0;
  
  // nastaveni parametru q1 az q5 pro Temelin !!!
  //q1 = 0.13563462295081966;
  //q2 = 1.0319759698093798;
  //q3 = 0.010231549384924803;
  //q4 = 0.059358593798123332;
  //q5 = 2.0528081279457426;
  
  //units:
  q1 = q1/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1
  q2 = q2/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1
  q3 = q3/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1
  q4 = q4/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1
  q5 = q5/6.89476*1.0e-3*1.0e-6;//units 1e-6*psi^-1 -> Pa^-1
  
  q = log10(retTime[1]/retTime[0]);//log spacing of consequent retardation times  
  jt = 0.;
  C_const = 0;
  
  //calculate degenerated (already crept) Kelvin series from min_ret_time to 10^-20 days
  ret_time_hlp = retTime[0];
  while (ret_time_hlp >= 1.e-20){
    ret_time_hlp /= pow(10.,q);//calculate next retardation time which is smaller in log scale
    C_const += give_L(q3, q, ret_time_hlp);
  }
  
  for(int i=0; i<nRetTime; i++){
    L = give_L(q3, q, retTime[i]);
    e_mu[i] = 1./L;
    aux = -(t-tl)/retTime[i];
    if (aux < log(DBL_MIN)) // exponent out of range for type of double
      jt+=L;
    else
      jt+=L*(1.-exp(aux));
  }
  
  jt = jt + q1 + C_const;
  
  //returns only e_mu
  return (jt);
}


/**
   function returns immediate compliance q1

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_q1()
{
  return(q1);
}


/**
   function returns ...

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_C_const()
{
  return(C_const);
}


/**
   function returns q4 ... 

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_q4()
{ 
  return(q4);
}


/**
   function returns constant for nonlinear creep effect

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_nonlin_func()
{
  double nonlin_func;
  
  nonlin_func = 1.0;

  return(nonlin_func);
}


/**
   function returns ... 

   @param time_mid - 

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_inv_v(double time_mid){
  
  double inv_v;

  inv_v = q2/q3 * pow(1.0/time_mid, m_param) + 1.;

  return(inv_v);
}


/**
   function computes retardation times

   @param rettimes    - %vector of ret. times
   @param n_ret_times - number of ret. times
   @param ipp         - number of integration point

   TKr, 3/6/2005
*/
void rspecmat::give_rettimes (vector &rettimes,long /*n_ret_times*/,long /*ipp*/)
{
  long i;

  //retardation times on a spectrum of non-aging log-power law from Laplace transformation
  rettimes.a[0] = retTime[0] = min_ret_time;
  
  for (i=1;i<nRetTime;i++){
    rettimes.a[i] = retTime[i] = retTime[0] * pow(10.,i);//less spacing does not make any significant improvement 
  }
}


/**
   function stores E_mu stiffnesses of Kelvin chain units into eqother

   @param e_mu        - %vector of stiffnesses
   @param n_ret_times - number of ret. times
   @param ipp         - number of integration point
   @param ido         - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::store_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  if (n_ret_times != nRetTime){//check of number of ret. times
    fprintf (stderr,"\n Wrong number of ret. times in function store_emu (file %s, line %d).\n",__FILE__,__LINE__);
    abort();
  }
  
  for (i=0;i<nRetTime;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2)+ i]  = e_mu.a[i];
  }
}


/**
   function returns E_mu stiffnesses of Kelvin chains

   @param e_mu        - %vector of stiffnesses
   @param n_ret_times - number of ret. times
   @param ipp         - number of integration point
   @param ido         - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::give_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  if (n_ret_times != nRetTime){//check of number of ret. times
    fprintf (stderr,"\n Wrong number of ret. times in function store_emu (file %s, line %d).\n",__FILE__,__LINE__);
    abort();
  }
  
  for (i=0;i<nRetTime;i++){
    e_mu.a[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2)+ i];
  }
}


/**
   function stores actual Young's modulus into eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 09/08/2008
*/
void rspecmat::store_ym_eqother(double ym,long ipp,long ido)
{
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1)] = ym;
}


/**
   function returns actual Young's modulus from eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 09/08/2008
*/
double rspecmat::give_ym_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double ym;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  ym = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1)];

  return(ym);
}


/**
   function stores previous Young's modulus into eqother

   @param ym       - previous Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 19/09/2014
*/
void rspecmat::store_ym_old_eqother(double ym,long ipp,long ido)
{
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1+nc+1)] = ym;
}


/**
   function returns previous Young's modulus from eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 19/09/2014
*/
double rspecmat::give_ym_old_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double ym;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  ym = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1+nc+1)];

  return(ym);
}




/**
   function computes of beta_t which reduces or increases retardation times

   @param ipp      - number of integration point

   TKr, 09/08/2008
   TKo, 06/08/2013, added case without temperature influence
*/
double rspecmat::compute_beta_t(long ipp)
{
  long nc,n_ret_times;
  double beta,tempr,uh_R,temprinit;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();

  if ((Mm->givestatusnmq(temperature) == 0) ||  (Mm->givestatusnmq(initial_temperature) == 1))
    return 1.0;

  tempr = Mm->givenonmechq(temperature, ipp);
  //toto musime casem upravit??!!
  if (tempr < 273.15)
    tempr += 273.15;//pro stupne C
  
  temprinit = Mm->givenonmechq(initial_temperature, ipp);
  //toto musime casem upravit??!!
  if (temprinit < 273.15)
    temprinit += 273.15;//pro stupne C
  
  uh_R = 4600.0*pow((30.0/(tempr - 263.0)),0.39);
  
  //varianty vypoctu??!!
  beta = 1.0;
  if(fabs(tempr - temprinit) > 0.001){
    if(kappa > 0.0){
      beta = exp(uh_R*(1.0/temprinit - 1.0/tempr));
      //beta = exp(3000*(1.0/temprinit - 1.0/tempr))*(0.1 + (1.0 - 0.1)*0.6*0.6);//reduced time
      if(kappa != 1.0)
	beta = kappa;
    }
  }

  //   printf("\nbeta = %lf",beta);//debug??!!
  return beta;
}


/**
   function stores actual beta_t into eqother

   @param beta     - beta_t
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 09/08/2008
*/
void rspecmat::store_beta_t_eqother(double beta,long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double beta_last;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  beta_last = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1)];

  if(beta_last < beta)//toto musime casem upravit??!!
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1)] = beta;
}


/**
   function gives actual beta_t from eqother

   @param beta     - beta_t
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 09/08/2008
*/
double rspecmat::give_beta_t_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double beta;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  beta = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1)];

  return beta;
}



/**
   function stores actual tensile strenght ft into eqother

   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 17/04/2009
*/
void rspecmat::store_ft_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1)+nc] = ft;
}


/**
   function gives actual tensile strenght ft from eqother
   
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array
   
   TKr, 17/04/2009
*/
double rspecmat::creep_give_actual_ft(long ipp,long /*im*/,long ido)
{
  long n_ret_times;
  long nc;
  double ftt;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  ftt = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1)+nc];

  return ftt;
}



/**
   function returns hidden strains (gamma_mu)

   @param gamma_mu - %matrix of hidden strains
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::give_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido)
{
  long i,j,ii;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();

  ii = 0;
  for (i=0;i<n_ret_times;i++){
    for (j=0;j<nc;j++){
      gamma_mu[j][i]=Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times)+ ii];
    ii++;
    }
  }
}



/**
   function stores hidden strains (gamma_mu) into eqother

   @param gamma_mu - %matrix of hidden strains
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::store_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido)
{
  long i,j,ii;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();

  ii = 0;
  for (i=0;i<n_ret_times;i++){
    for (j=0;j<nc;j++){
      Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times)+ ii] = gamma_mu[j][i];
    ii++;
    }
  }
}


/**
   function returns total stresses from eqother

   @param sigma - %vector of total stresses
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::give_stresses_eqother(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    sigma[i] = Mm->ip[ipp].eqother[ido+(nc)+i];
  }
}


/**
   function stores total stresses into eqother

   @param sigma - %vector of total stresses
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void rspecmat::store_stresses_eqother(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc)+i] = sigma[i];
  }
}


/**
   function returns total stresses from other

   @param sigma - %vector of total stresses
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_stresses_other(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    sigma[i] = Mm->ip[ipp].other[ido+(nc)+i];
  }
}


/**
   function stores total stresses into other

   @param sigma - %vector of total stresses
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_stresses_other(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].other[ido+(nc)+i] = sigma[i];
  }
}


/**
   function returns increments of total stresses from eqother

   @param dsigma - %vector of increments of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_dstresses_eqother(vector &dsigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    dsigma[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc)+i];
    //fprintf (Out,"\n Vraci dsigma = %e\n",dsigma[i]);//debug??!!
  }
}


/**
   function stores increments of total stresses into eqother

   @param dsigma - %vector of increments of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_dstresses_eqother(vector &dsigma,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc)+i] = dsigma[i];
    //fprintf (Out,"\n Uklada dsigma = %e\n",dsigma[i]);//debug??!!
  }
}


/**
   function returns total strains from eqother

   @param eps - %vector of total strains
   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_strains_eqother(vector &eps,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    eps[i] = Mm->ip[ipp].eqother[ido+(0)+i];
    //fprintf (Out,"\n Vraci total strains eps = %e\n",eps[i]);//debug??!!
    //fprintf (Out,"\n Pro srovnani Mm->ip[0].strain[0] = %e\n",Mm->ip[ipp].strain[i]);//debug??!!
  }
}


/**
   function stores total strains from eqother

   @param eps - %vector of total strains
   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_strains_eqother(vector &eps,long ipp,long ido)
{
  long i;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(0)+i] = eps[i];
    //fprintf (Out,"\n Uklada total strains eps = %e\n",eps[i]);//debug??!!
    //fprintf (Out,"\n Pro srovnani Mm->ip[0].strain[0] = %e\n",Mm->ip[ipp].strain[i]);//debug??!!
  }
}


/**
   function returns increments of irreversible creep strains

   @param deps_cr - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_creepdstrains_eqother(vector &deps_cr,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    deps_cr[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times)+i];
    //fprintf (Out,"\n Vraci deps_cr = %e\n",deps_cr[i]);//debug??!!
  }
}


/**
   function stores increments of irreversible creep strains into eqother

   @param deps_cr - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_creepdstrains_eqother(vector &deps_cr,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times)+i] = deps_cr[i];
    //fprintf (Out,"\n Uklada deps_cr = %e\n",deps_cr[i]);//debug??!!
  }
}


/**
   function returns increments of irreversible shrinkage strains from eqother
\
   @param deps_sh - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array 

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_irrdstrains_eqother(vector &deps_sh,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    deps_sh[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc)+i];
  }
}


/**
   function storess increments of irreversible shrinkage strains into eqother

   @param deps_sh - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_irrdstrains_eqother(vector &deps_sh,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc)+i] = deps_sh[i];
  }
}


/**
   function returns increments of irreversible stress-induced shrinkage strains from eqother

   @param deps_ss - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::give_stressirrdstrains_eqother(vector &deps_ss,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    deps_ss[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc)+i];
  }
}


/**
   function stores increments of irreversible stress-induced shrinkage strains into eqother

   @param deps_ss - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_stressirrdstrains_eqother(vector &deps_ss,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc)+i] = deps_ss[i];
  }
}


/**
   function returns free shrinkage

   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
double rspecmat::give_shrinkage_eqother(long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;
  double eps_sh=0.0;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    eps_sh = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc)];
  }

  return eps_sh;
}


/**
   function stores free shrinkage into eqother

   @param eps_sh  - shrinkage strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_shrinkage_eqother(double eps_sh,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc)] = eps_sh;
  }
}


/**
   function stores actual humidity

   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_hum_eqother(long ipp,long ido)
{
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  if(Mm->givestatusnmq(rel_hum) == 1)
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)] = Mm->givenonmechq(rel_hum, ipp);
  else
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)] = 0.0;
}


/**
   function stores actual temperature

   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::store_temp_eqother(long ipp,long ido)
{
  long nc;

  nc=Mm->ip[ipp].ncompstr;

  if(Mm->givestatusnmq(temperature) == 1)
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+1)] = Mm->givenonmechq(temperature, ipp);
  else
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+1)] = 0.0;
}


/**
   function returns number of times

   @param ipp - number of integration point

   TKr, 3/6/2005
*/
long rspecmat::give_napproxtime(long /*ipp*/)
{
  return(napproxtime);
}


/**
   function computes increment of irreversible free strains in a material point by B3 Bazant's model
   from temerature and humidity changes (drying shrikage, shrinkage, stress induced strain etc.)
   free thermal strain is not included!!      

   @param deps_sh   - %vector of free shrinkage strains
   @param t0        - age when drying begins [days]
   @param t_dt      - time representing age of concrete + dt [days]
   @param t         - time representing age of concrete [days]
   @param ipp       - number of integration point
   @param ido       - index of internal variables for given material in the ipp other array
   
   TKr, 28/8/2008 - revised
*/
void rspecmat::give_deps_free (vector &deps_sh, double t0, double t_dt, double t, double /*dt*/, long ipp,long im,long ido)
{
  double eps_shinf,kh,kt,tau,st,dhum,dtemp;
  double hum,temp,hum_prev,temp_prev;
  double deps,deps_h,deps_h_old,deps_a,deps_t,e_607,e_t0_tau,e_t0,e_actual;
  long nc;
  //  stress/strain state
  strastrestate ss;

  nc=Mm->ip[ipp].ncompstr;
    
  if ( type_h == 0){
    hum = 0.0;
    hum_prev = 0.0;
  }
  else{
    // actual humidity
    hum=Mm->givenonmechq(rel_hum, ipp);
    // humidity from previous time step
    hum_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)];
  }
  if ( type_temp == 0){
    temp = 0.0;
    temp_prev = 0.0;
  }
  else{
    // actual temperature
    temp=Mm->givenonmechq(temperature, ipp);
    // temperature from previous time step
    temp_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+1];
  }
  
  if (t < t0){
    fprintf (stderr,"\n Age of concrete is lower than age when dyring begins!!!");
    fprintf (stderr,"\n Age of concrete = %10.15e,   Age when dyring begins = %10.15e",t,t0);
    print_err("\n In ", __FILE__, __LINE__, __func__);
    exit(1);
  }

  //humidity effect
  if (type_h == 0){
    //constant huminidy
    //kt = 190.8/pow(t0,0.08)*fc;//this is wrong this is from literature (t0 = age when drying begins)
    kt = 190.8/pow(t,0.08)/pow(fc,0.25);//this is according to creepb.cpp (Fajman)
    tau = kt*ks*ks*kd*kd;
    st = tanh(sqrt((t_dt - t0)/tau));
    //derivative of st with respect to time
    //dst = tanh(sqrt((t_dt - t0)/tau)) - tanh(sqrt((t - t0)/tau));
    eps_shinf = a1*a2*(26.0*pow((wc*cs),2.1)/pow(fc,0.28)+270.); // 26 should be replaced by 0.019 for SI units
    //time dependence of ultimate shrinkage = effect of aging
    e_607 = e28*sqrt(607.0/(4+0.85*607.0));// empirical formulation
    e_t0_tau = e28*sqrt((t0+tau)/(4+0.85*(t0+tau)));// empirical formulation
    eps_shinf = eps_shinf*e_607/e_t0_tau;//zatim??!!

    if(hum_env<=0.98)
      kh=(1.0-pow(hum_env,3.0));
    else 
      kh= (-0.2 - (1.0-pow(0.98,3.0)))/0.02*(hum_env - 0.98);
    if(hum_env>=1.0)
      kh = -0.2;

    deps_h = -eps_shinf*kh*st; //shrinkage strain

    deps_h_old = give_shrinkage_eqother(ipp,ido);
    store_shrinkage_eqother(deps_h,ipp,ido);
    deps_h = deps_h - deps_h_old; //shrinkage strain increment

    deps_h = deps_h*1.e-6;//
  }    
  else{
    //changing humidity
    dhum = hum - hum_prev;
    eps_shinf = a1*a2*(26.0*pow((wc*cs),2.1)/pow(fc,0.28)+270.);
    eps_shinf = eps_shinf*1.e-6;//units
    //eps_shinf = 0.0008;//new 24.2.2006 corresponds to Pa,m,N

    //time dependence of ultimate shrinkage
    e_t0 = e28*sqrt((t0)/(4+0.85*(t0)));// empirical formulation
    e_actual = creep_give_actual_ym (ipp,im,ido)/6.89476*1.0e-3;//to psi
    //eps_shinf = eps_shinf*e_t0/e_actual;

    deps_h = eps_shinf*3.0*hum*hum*dhum; //shrinkage strain increment

    //simple empirical formula for autogenous shrinkage
    //added total autogenous shrinkage
    deps_a = -eps_ainf*((1.0-exp(-.125*t_dt*24.0))-(1.0-exp(-.125*t*24.0)));//time in hours
    
    deps_h = deps_h + deps_a;
  }
  

  //temperature effect
  dtemp = temp - temp_prev;  
  deps_t = 0.0;//deps_t = alpha*dtemp;//free thermal strain is not included now
  
  if (flag_shr == 0)
    deps_h = 0.0;
  
  if (flag_temp == 0)
    deps_t = 0.0;
  
  deps = deps_h + deps_t;
   
  fillv(0.0,deps_sh);
  ss=Mm->ip[ipp].ssst;
  
  //only uniaxial
  if(ss==bar){
    deps_sh[0]=deps;
  }
  else{
    deps_sh[0]=deps;
    deps_sh[1]=deps;
  }
  if(ss==planestrain) deps_sh[3]=deps;
  if(ss==planestress) deps_sh[3]=deps;
  if(ss==axisymm) deps_sh[2]=deps;
  if(ss==spacestress) deps_sh[2]=deps;

}


/**
   function computes increment of irreversible stress induced strains in a material point B3 Bazant's model
   from temerature and humidity changes (drying shrikage, shrinkage, stress induced strain etc.)
   
   @param deps_ss   - %vector of stress-induced shrinkage strains
   @param t0        - age when drying begins [days]
   @param sigma     - actual uniaxial stress [Pa]
   @param t_dt      - time representing age of concrete + dt [days]
   @param t         - time representing age of concrete [days]
   @param ipp       - number of integration point
   @param ido       - index of internal variables for given material in the ipp other array

   TKr, 28/8/2008 - revised   
*/
void rspecmat::give_deps_stressinduced (vector &deps_ss, double t0, double /*t_dt*/, double t, vector &sigma, long ipp,long im,long ido)
{
  double eps_shinf,deps,deps_cs,deps_ts,r,gh,dhum,rho,gt,dtemp,k;
  double hum,temp,hum_prev,temp_prev,e_t0,e_actual;
  long i,nc;
  //  stress/strain state
  strastrestate ss;

  nc=Mm->ip[ipp].ncompstr;
  
  if ( type_h == 0){
    hum = 0.0;
    hum_prev = 0.0;
  }
  else{
    // actual humidity
    hum=Mm->givenonmechq(rel_hum, ipp);
    // humidity from previous time step
    hum_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)];
  }
  if ( type_temp == 0){
    temp = 0.0;
    temp_prev = 0.0;
  }
  else{
    // actual temperature
    temp=Mm->givenonmechq(temperature, ipp);
    // temperature from previous time step
    temp_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+1];
  }

  if (t < t0){
    fprintf (stderr,"\n Age of concrete is lower than age when dyring begins!!!");
    fprintf (stderr,"\n In function deps_stressinduced (file %s, line %d).\n",__FILE__,__LINE__);
    abort();
  }

  eps_shinf = a1*a2*(26.0*pow((wc*cs),2.1)/pow(fc,0.28)+270.);
  eps_shinf = eps_shinf*1.e-6;//units
  //eps_shinf = 0.0008;//new 24.2.2006 corresponds to Pa,m,N

  //time dependence of ultimate shrinkage
  e_t0 = e28*sqrt((t0)/(4+0.85*(t0)));// empirical formulation
  e_actual = creep_give_actual_ym (ipp,im,ido)/6.89476*1.0e-3;//to psi
  //eps_shinf = eps_shinf*e_t0/e_actual;//zatim??!!

  k = -eps_shinf;
  r = 0.3*ft*1.e-6;//in MPa^-1
  rho = 1.5*ft*1.e-6;//in MPa^-1

  dhum = hum-hum_prev;
  dtemp = temp - temp_prev;

  if(dhum >= 0.0)
    gh = 1.0;
  else
    gh = -1.0;

  if(dtemp >= 0.0)
    gt = 1.0;
  else
    gt = -1.0;

  
  fillv(0.0,deps_ss);
  ss=Mm->ip[ipp].ssst;  

  for(i=0;i<nc;i++){
    if (type_h == 0)
      deps_cs = 0.0;
    else
      deps_cs = -k*r*gh*sigma[i]*dhum;
    
    if (type_temp == 0)
      deps_ts = 0.0;
    else
      deps_ts = -alpha*rho*gt*sigma[i]*dtemp;
    
    deps = (deps_cs + deps_ts)*1.e-12;//units
      
    if (flag_drshr == 0)
      deps = 0.0;
    
    deps_ss[i]=deps;
  }

  //only uniaxial
  if(ss==planestress){
    deps_ss[2]=0.0;
  }
  if(ss==planestrain){
    deps_ss[2]=0.0;
  }
  if(ss==axisymm){
    deps_ss[3]=0.0;
  }
  if(ss==spacestress){
    deps_ss[3]=0.0;
    deps_ss[4]=0.0;
    deps_ss[5]=0.0;
  }

}


/** 
   function returns actual tensile strenght
    
   @param ipp - index of integration point
   @param im  - index of material
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double rspecmat::compute_actual_ft (long ipp,long im,long ido) 
{ 
  long imat,nc;
  double q,e_0,e_t,fft,fft_prev,temp,temp_prev;
  
  nc = Mm->ip[ipp].ncompstr;

  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //creep coefficient
  q = creep_matstiffchange (ipp,im,ido);
  //asympotic elastic modulus 
  e_0=Mm->eliso[imat].e;
  e_t = q*e_0;
  
  switch(ft_flag){
    case 0:
    fft = ft;
    break;
  case 1:
    fft = ft_ratio*e_t;
    break;
  case 2:{//temperature effect on tensile streght
    temp=Mm->givenonmechq(temperature, ipp);
    // temperature from previous time step
    temp_prev = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+1];
    fft_prev = creep_give_actual_ft(ipp,im,ido);
    if(temp > temp_prev)
      fft = fft_prev - ft_ratio*(temp - temp_prev);
    else
      fft = fft_prev;
    break;
  }
  default:{
    print_err("\n unknown type ft_flag is required \n",__FILE__, __LINE__, __func__);
    abort();
  }
  }

  return (fft);
} 


/**
   function adds increments of irreversible creep strains to total irreversible creep strains in eqother

   @param deps - %vector of increments of strains
   @param ipp  - number of integration point
   @param ido  - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::addirrstrains_eqother (vector &deps,long ipp, long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;
  
  n_ret_times = give_nret_time ();

  for (i=0;i<deps.n;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc)+i] += deps[i];
  }
}


/**
   function returns total irreversible creep strains from eqother

   @param epscr - %vector of increments of strains
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void rspecmat::giveirrstrains_eqother (long ipp, long ido, vector &epscr)
{
  long i;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  n_ret_times = give_nret_time ();
  
  for (i=0;i<epscr.n;i++){
    epscr[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc)+i];
  }
}



/**
   function stores total "aeging" strains into eqother

   @param eps_ag - %vector of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 17/3/2009
*/
void rspecmat::store_agstrains_eqother(vector &eps_ag,long ipp,long ido)
{
  long i;
  long n_ret_times,nc;


  nc=Mm->ip[ipp].ncompstr;

  n_ret_times = give_nret_time ();

  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1)+i] = eps_ag[i];
  }
}


/**
   function returns total "aeging" strains from eqother

   @param eps_ag - %vector of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 17/3/2009
*/
void rspecmat::give_agstrains_eqother(vector &eps_ag,long ipp,long ido)
{
  long i;
  long n_ret_times,nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    eps_ag[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+2+n_ret_times+nc*n_ret_times+nc+nc+nc+nc+1+1+1)+i];
  }
}



/**
   function returns other value

   @param compother - number of other components
   @param ipp       - index of the first integration point on element

   TKr, 28/8/2008 - revised
*/
double rspecmat::get_othervalue(long compother,long /*ipp*/)
{
  double other;

  //not finished
  switch (compother){ 
  case 0:{//eps_x
    other = 0.0;
      break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
    abort();
  }
  }
  return (other);

}


/**
   function prints name of other value

   @param out       - output file
   @param compother - number of other component

   TKr, 28/8/2005 - revised
*/
void rspecmat::print_othervalue_name(FILE *out,long compother)
{
  //not finished
  switch (compother){
  case 0:{//eps_x
    fprintf (out,"eps_x ()");
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
*/
void rspecmat::give_reqnmq(long *anmq)
{
  if (type_temp == 1)
  {
    anmq[temperature-1] = 1;
    anmq[initial_temperature-1] = 1;
  }
  if (type_h == 1)
    anmq[rel_hum-1] = 1;
}
