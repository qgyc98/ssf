/** 
    Concrete creep double-power-law model according to Bazant.
*/
#include "creep_dpl.h"
#include "iotools.h"
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


dplmat::dplmat (void)
{
  q1 = 0.0;

  e28=30.0e9; //Pa
  fc=35.8e6; //Pa
  ro = 2500.0; //kg*m-3
  wc=0.43;
  sc=3.4;
  gc=1.98;
  a1 = 1.05;

  tb_time = 0.0;
  th_time = 0.0;
  napproxtime = 0;
  nRetTime = 0;
  type_rt = 0;
  ft_flag = 0;
  ft = 0.0;
  ft_ratio = 1.0e-5;

  e0 = 0.0;

  previoustime = 0.0;
  actualtime = 0.0;
  dtb = 0.0;
  tb_age_dt = tb_age = tbl_age = tbh_age = maxtimeb = 0.0;
  napptimeb = 0;

  retTime=NULL;
  emu=NULL;
  timeMax=(Mp->timecon.endtime ())/86400.0;//in days
}


dplmat::~dplmat (void)
{
  delete [] retTime;
}


/**
   function reads material parameters

   @param in - input file
   
   TKr, 28/8/2008 - revised
*/
void dplmat::read (XFILE *in)
{
  long i;

  xfscanf(in,"%ld",&type_e);
  
  if (type_e == 1){xfscanf (in,"%lf ",&e28);}
  //input in Pa
  e0 = 1.5*e28/6.89476*1.0e-3;//to psi
  e28 = e28/6.89476*1.0e-3;//to psi
  
  xfscanf (in,"%lf %lf %lf %lf %lf %lf",&fc,&ro,&wc,&sc,&gc,&a1);
  
  //input in Pa
  fc = fc/6.89476*1.0e-3;//to psi
  //input in kg*m^-3
  ro = ro/16.03;//to lb*ft^-3
  
  xfscanf (in,"%lf %lf %ld %ld %ld",
	  &tb_time,&th_time,&napproxtime,&nRetTime,&type_rt);
  
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
}


/**
   function prints material parameters

   @param out - output file

   TKr, 28/8/2008 - revised
*/
void dplmat::print (FILE *out)
{
  long i;

  fprintf(out," %ld",type_e);
  
  if (type_e == 1){fprintf (out," %e ",e28);}
  
  fprintf (out," %lf %lf %lf %lf %lf %lf",fc,ro,wc,sc,gc,a1);
  
  fprintf (out," %lf %lf %ld %ld %ld",
	   tb_time,th_time,napproxtime,nRetTime,type_rt);
  
  if (type_rt==1){//retardation times (ages) [days]
    for (i=0;i<nRetTime;i++)
      fprintf (out," %lf",retTime[i]);
  }
}


/**
   function computes ages in days

   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2008 - revised
*/
void dplmat::compute_ages (long ipp,long /*ido*/)
{
  //times and ages in days
  if(Mp->time == previoustime)
    //    dtb = Mp->timecon.initialtimeincr ()/86400.0;
    dtb = Mp->timecon.actualforwtimeincr()/86400.0;
  else
    dtb = (Mp->time - previoustime)/86400.0;
  actualtime = Mp->time/86400.0;
  tbh_age = actualtime - th_time;
  tbl_age = actualtime - tb_time + dtb/2.0;
  tb_age = actualtime - tb_time;
  napptimeb = give_napproxtime(ipp);

  maxtimeb = Mp->timecon.endtime ()/86400.0;
  tb_age_dt = tb_age + dtb;

  if (tb_age <= 0.0){
    fprintf (stderr,"\n Age of concrete must be greater than zero!!!");
    fprintf (stderr,"\n Age of concrete = %e, Actual time = %e, Time of end of concrete casting = %e; (Age of concrete = Actual time - Time of end of concrete casting)",tb_age,actualtime,tb_time);
    print_err("\n In function compute_ages ", __FILE__, __LINE__, __func__);
    exit(0);
  }
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
void dplmat::give_ages (double &t_age_dt,double &t_age,double &tl_age,double &th_age,double &dt,double &maxtime,long &napptime,long /*ipp*/)
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
long dplmat::give_nceqother (long ipp)
{
  long nc,nceqother;

  nc=Mm->ip[ipp].ncompstr;

  nceqother = nc + nc;                  //total strains, total stresses
  nceqother = nceqother + nc + nc;      //strain increments, stress increments
  nceqother = nceqother + nRetTime;     //coefficients of Dirichlet series
  nceqother = nceqother + nc*nRetTime;  //hidden strains increments
  nceqother = nceqother + nc;           //creep-irreversible strains increments
  nceqother = nceqother + nc;           //total irreversible strains (creep-irreversible strains)
  nceqother = nceqother + 1;            //actual Young's modulus
  nceqother = nceqother + nc;           //total "aeging" strains = creep-irreversible strains + shrinkage-irreversible strains + stress-induced shrinkage-irreversible strains + strains due to aegin of concrete - changing effective Young's modulus
  nceqother = nceqother + 1;            //previous Young's modulus
  
  return(nceqother);
}


/**
   function returns number of components of other array

   @param ipp - number of integration point

   TKr, 07/08/2008, not used now
*/
long dplmat::give_ncother (long ipp)
{
  long nc,ncother;

  nc=Mm->ip[ipp].ncompstr;

  ncother = nc;

  return(ncother);
}


/**
   function returns number of retardation times

   TKr, 3/6/2005
*/
long dplmat::give_nret_time (void)
{
  return(nRetTime);
}


/**
   function computes retardation times

   @param rettimes    - %vector of ret. times
   @param n_ret_times - number of ret. times
   @param ipp         - number of integration point

   TKr, 3/6/2005
*/
void dplmat::give_rettimes (vector &rettimes,long n_ret_times,long /*ipp*/)
{
  long i;
  double m,mm;

  if (n_ret_times != nRetTime){//check of number of ret. times
    print_err("\n Wrong number of ret. times", __FILE__, __LINE__, __func__);
    exit(0);
  }

  if (type_rt==1){
  }
  else
    {
      m = log10(2.0*timeMax);
      
      retTime[0]=1.0e-9;
      retTime[1]=1.0;
      mm=1./(nRetTime-2);
      
      for (i=2;i<nRetTime-1;i++){
	retTime[i]=pow(10.,(i-1)*mm*m);
      }
      retTime[nRetTime -1]=2.0*timeMax;
    }
  
  for (i=0;i<nRetTime;i++){
    rettimes.a[i]=retTime[i];
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
void dplmat::store_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;

  if (n_ret_times != nRetTime){//check of number of ret. times
    print_err("\n Wrong number of ret. times", __FILE__, __LINE__, __func__);
    exit(0);
  }
  
  for (i=0;i<nRetTime;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+ i]  = e_mu.a[i];
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
void dplmat::give_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;

  if (n_ret_times != nRetTime){//check of number of ret. times
    print_err("\n Wrong number of ret. times", __FILE__, __LINE__, __func__);
    exit(0);
  }
  
  for (i=0;i<nRetTime;i++){
    e_mu.a[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc)+ i];
  }
}


/**
   function stores actual Young's modulus into eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 09/08/2008
*/
void dplmat::store_ym_eqother(double ym,long ipp,long ido)
{
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc)] = ym;
}


/**
   function returns actual Young's modulus from eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 17/03/2009
*/
double dplmat::give_ym_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double ym;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  ym = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc)];

  return(ym);
}


/**
   function stores previous Young's modulus into eqother

   @param ym       - previous Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 19/09/2014
*/
void dplmat::store_ym_old_eqother(double ym,long ipp,long ido)
{
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc+1+nc)] = ym;
}


/**
   function returns previous Young's modulus from eqother

   @param ym       - actual Young's modulus
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 19/09/2014
*/
double dplmat::give_ym_old_eqother(long ipp,long ido)
{
  long n_ret_times;
  long nc;
  double ym;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  ym = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc+1+nc)];

  return(ym);
}


/**
   function returns hidden strains (gamma_mu)

   @param gamma_mu - %matrix of hidden strains
   @param ipp      - number of integration point
   @param ido      - index of internal variables for given material in the ipp other array

   TKr, 3/6/2005
*/
void dplmat::give_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido)
{
  long i,j,ii;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();

  ii = 0;
  for (i=0;i<n_ret_times;i++){
    for (j=0;j<nc;j++){
      gamma_mu[j][i]=Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times)+ ii];
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
void dplmat::store_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido)
{
  long i,j,ii;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();

  ii = 0;
  for (i=0;i<n_ret_times;i++){
    for (j=0;j<nc;j++){
      Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times)+ ii] = gamma_mu[j][i];
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
void dplmat::give_stresses_eqother(vector &sigma,long ipp,long ido)
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
void dplmat::store_stresses_eqother(vector &sigma,long ipp,long ido)
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
void dplmat::give_stresses_other(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    sigma[i] = Mm->ip[ipp].eqother[ido+(nc)+i];
  }
}


/**
   function stores total stresses into other

   @param sigma - %vector of total stresses
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::store_stresses_other(vector &sigma,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc)+i] = sigma[i];
  }
}


/**
   function returns increments of total stresses from eqother

   @param dsigma - %vector of increments of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::give_dstresses_eqother(vector &dsigma,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    dsigma[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc)+i];
  }
}


/**
   function stores increments of total stresses into eqother

   @param dsigma - %vector of increments of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::store_dstresses_eqother(vector &dsigma,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc)+i] = dsigma[i];
  }
}


/**
   function returns total strains from eqother

   @param eps - %vector of total strains
   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::give_strains_eqother(vector &eps,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    eps[i] = Mm->ip[ipp].eqother[ido+(0)+i];
  }
}


/**
   function stores total strains from eqother

   @param eps - %vector of total strains
   @param ipp - number of integration point
   @param ido - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::store_strains_eqother(vector &eps,long ipp,long ido)
{
  long i;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(0)+i] = eps[i];
  }
}


/**
   function returns increments of irreversible creep strains

   @param deps_cr - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::give_creepdstrains_eqother(vector &deps_cr,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    deps_cr[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times)+i];
  }
}


/**
   function stores increments of irreversible creep strains into eqother

   @param deps_cr - %vector of increments of strains
   @param ipp     - number of integration point
   @param ido     - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::store_creepdstrains_eqother(vector &deps_cr,long ipp,long ido)
{
  long i;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times)+i] = deps_cr[i];
  }
}


/**
   function returns number of times

   @param ipp - number of integration point

   TKr, 3/6/2005
*/
long dplmat::give_napproxtime(long /*ipp*/)
{
  return(napproxtime);
}


/**
   function computes creep compliacne function from Double power law

   @param tt        - age at loading [s]
   @param t         - time representing age of concrete [s]
   @param ipp       - number of integration point
   @param ido       - index of internal variables for given material in the ipp other array
  
   TKr, 3/6/2005
*/
double dplmat::double_power_law(double tt,double t,long /*ipp*/,long /*ido*/)
{
  double ee,a,n,m,fi;
  double x,ag,ac,jt;

  if (type_e == 1){//measured
    ee = 600000.0/e28;//measured
  }
  else {//empirical
    ee = (0.09+1/(1.7*(0.00005*ro*ro*fc/1000.0)*(0.00005*ro*ro*fc/1000.0)));//this is for fc in psi
  }
  
  m=0.28+1/(fc*fc/1000.0/1000.0);
  a=1./(40.0*wc); 
  ag=(sc+gc)/gc;
  ac=sc+gc;
  x=(2.1*ac/pow(sc,1.4)+0.1*pow((fc/1000.0),1.5)*pow(wc,1./3.)*pow(ag,2.2))*a-4.0;//this is for fc in psi
  if (x <= 4.0)
    n=0.12;
  else
    n=.12+(0.07*pow(x,6.)/(5130.0+pow(x,6.)));
  
  fi=0.5*pow(10.0,(3.0*n))/(pow(28.0,-m)+a);
  
  q1 = ee;
  jt = fi*ee*(pow(tt,-m)+a)*pow((t-tt),n);

  q1 = q1/6.89476*1.0e-3*1.0e-6;//units 1e-6*ksi^-1 -> Pa^-1
  jt = jt/6.89476*1.0e-3*1.0e-6;//units 1e-6*ksi^-1 -> Pa^-1

  return(jt);
}


/**
   function returns immediate compliance q1

   TKr, 28/8/2005 - revised
*/
double  dplmat::give_q1()
{
  return(q1);
}


/** 
   function returns actual tensile strenght
    
   @param ipp - index of integration point
   @param im  - index of material
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double dplmat::creep_give_actual_ft (long ipp,long im,long ido) 
{ 
  long imat;
  double q,e_0,e_t,fft;
  
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //creep coefficient
  q = creep_matstiffchange (ipp,im,ido);
  //asympotic elastic modulus 
  e_0=Mm->eliso[imat].e;
  e_t = q*e_0;
  
  if(ft_flag == 0)
    fft = ft;
  else
    fft = ft_ratio*e_t;
  
  return (fft);
} 


/**
   function adds increments of irreversible creep strains to total irreversible creep strains in eqother

   @param deps - %vector of increments of strains
   @param ipp  - number of integration point
   @param ido  - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::addirrstrains_eqother (vector &deps,long ipp, long ido)
{
  long i;
  long n_ret_times;
  long nc;

  nc=Mm->ip[ipp].ncompstr;  
  n_ret_times = give_nret_time ();

  for (i=0;i<deps.n;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc)+i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc)+i] + deps[i];
  }
}


/**
   function returns total irreversible creep strains from eqother

   @param epscr - %vector of increments of strains
   @param ipp   - number of integration point
   @param ido   - index of internal variables for given material in the ipp other array

   TKr, 28/8/2005 - revised
*/
void dplmat::giveirrstrains_eqother (long ipp, long ido, vector &epscr)
{
  long i;
  long n_ret_times;
  long nc;
  
  nc=Mm->ip[ipp].ncompstr;
  
  n_ret_times = give_nret_time ();
  
  for (i=0;i<epscr.n;i++){
    epscr[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc)+i];
  }
}


/**
   function stores total "aeging" strains into eqother

   @param eps_ag - %vector of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 17/3/2009
*/
void dplmat::store_agstrains_eqother(long ipp,long ido,vector &eps_ag)
{
  long i;
  long n_ret_times,nc;


  nc=Mm->ip[ipp].ncompstr;

  n_ret_times = give_nret_time ();

  for (i=0;i<nc;i++){
    Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc+1)+i] = eps_ag[i];
  }
}


/**
   function returns total "aeging" strains from eqother

   @param eps_ag - %vector of total stresses
   @param ipp    - number of integration point
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 17/3/2009
*/
void dplmat::give_agstrains_eqother(long ipp,long ido,vector &eps_ag)
{
  long i;
  long n_ret_times,nc;

  nc=Mm->ip[ipp].ncompstr;
  n_ret_times = give_nret_time ();
  
  for (i=0;i<nc;i++){
    eps_ag[i] = Mm->ip[ipp].eqother[ido+(nc+nc+nc+nc+n_ret_times+nc*n_ret_times+nc+nc+1)+i];
  }
}



/**
   function returns other value

   @param compother - number of other components
   @param ipp       - index of the first integration point on element

   TKr, 28/8/2008 - revised
*/
double dplmat::get_othervalue(long compother,long /*ipp*/)
{
  double other;

  //not finished
  switch (compother){
  case 0:{//eps_x
    other = 0.0;
      break;
  }
  default:{
    print_err("\n Unknown type of component is required in ", __FILE__, __LINE__, __func__);
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
void dplmat::print_othervalue_name(FILE *out,long compother)
{
  //not finished
  switch (compother){
  case 0:{//eps_x
    fprintf (out,"eps_x ()");
    break;
  }
  default:{
    print_err("\n Unknown type of component is required in ", __FILE__, __LINE__, __func__);
  }
  }
}
