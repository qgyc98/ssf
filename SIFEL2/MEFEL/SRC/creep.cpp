#include "creep.h"
#include "iotools.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "elastisomat.h"
#include "creep_b3.h"
#include "creep_rspec.h"
#include "creepb.h"
#include "creep_dpl.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>



/**
   function returns actual stiffness matrix

   @param d   - stiffness %matrix
   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/
void creep_matstiff (matrix &d, long ipp, long im, long ido)
{
  double q;
  q = 1.0;
  // asymptotic elastic stiffness matrix
  Mm->elmatstiff(d, ipp, ido);  
  // actual stiffness  coefficient (modulus)
  q = creep_matstiffchange(ipp, im, ido);

  // actual stiffness matrix
  cmulm(q, d);
}


/**
   function returns actual stiffness (modulus)

   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/
double creep_matstiffchange (long ipp,long im,long ido)
{
  long i;
  long imat,napproxtime;
  double tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,q,delYi;
  double e0,et,qq,inv_v,q1,c_const,param_exp,exp_prec=-50;

  //number of retardation times
  long n_ret_times = creep_number_rettimes (ipp,im);
  vector ret_times(ASTCKVEC(n_ret_times)),emu(ASTCKVEC(n_ret_times));
  
  // memory clearing
  fillv(0.0, ret_times);
  fillv(0.0, emu);
  
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //asympotic elastic modulus
  e0=Mm->eliso[imat].e;

  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //times -> ages [days]
    //actuall times in days
    Mm->crb3[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //retardation times from ip
    Mm->crb3[i].give_rettimes(ret_times,n_ret_times,ipp);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crb3[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    //ageing function 1/v(t)    
    inv_v = 1.0;
    //inital compliance
    q1 = Mm->crb3[i].give_q1();
    // const
    c_const = 0.0;
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //times -> ages [days]
    //actuall times in days
    Mm->crrs[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //retardation times from ip
    Mm->crrs[i].give_rettimes(ret_times,n_ret_times,ipp);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crrs[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    //ageing function 1/v(t)    
    inv_v = Mm->crrs[i].give_inv_v(tl_age);//aritmetical mid_interval value
    //inital compliance
    q1 = Mm->crrs[i].give_q1();
    // const
    c_const = Mm->crrs[i].give_C_const();
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //times -> ages [days]
    //actuall times in days
    Mm->crdpl[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //retardation times from ip
    Mm->crdpl[i].give_rettimes(ret_times,n_ret_times,ipp);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crdpl[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    //ageing function 1/v(t)    
    inv_v = 1.0;
    //inital compliance
    q1 = Mm->crdpl[i].give_q1();
    // const
    c_const = 0.0;
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  //Computing of actual stiffness modulus from dirichlet series
  et=0.;
  q=2./3.;
  //Kelvin chain model
  for (i=0;i<n_ret_times;i++){
    if(tb_age<0.0001){
      delYi=pow((0.0001+dt)/ret_times[i],q)-pow(0.0001/ret_times[i],q);
    }
    else{
      delYi=pow((tb_age+dt)/ret_times[i],q)-pow(tb_age/ret_times[i],q);
    }

    param_exp = -delYi;
    if(param_exp <= exp_prec)//numerical correction
      param_exp = exp_prec;
    
    et=et+(delYi-1.+exp(param_exp))/delYi/emu[i];
  }

  if((Mm->ip[ipp].tm[im]) == creeprs){
    et = et + c_const;
    et = et * inv_v;
  }
  //added q1 term (inital compliance)
  et = et + q1;
  
  //Do not delete this part. This is here for testing creep models.
  /* double dd=0.0;
     long j;
     if (Mm->ip[ipp].tm[im] == creeprs){
     i=Mm->ip[ipp].idm[im];
     for(j=1;j<1000;j++){
     //et = Mm->crrs[i].give_J_E_mu(emu,th_age,tb_age,tb_age+dd,ipp,ido);
     //et = (et-q1)*inv_v;
     et = Mm->crb3[0].b3_law (th_age,tb_age,tb_age+dd,ipp,ido);
     et = et + q1;
     fprintf (Out,"%e   %e\n",tb_age+dd,-et);
     dd = dd + 1.0;
     }
     }
     fclose(Out);
  */
  
  //actual stiffness
  if(et>1.e-50)
    et=1./et;
  else{
    print_err("stiffness modulus is zero (et=%le) on ip=%ld, eid=%ld.\n",
              __FILE__, __LINE__, __func__, et, ipp, Mm->elip[ipp]+1);
    abort();
  }

  // actual material modulus
  qq=et/e0;

  return(qq);
}


/** 
   function computes actual Young's modulus
    
   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double creep_compute_actual_ym (long ipp,long im,long ido) 
{ 
  long i,imat;
  double q,e_0,e_t;

  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //creep coefficient
  q = creep_matstiffchange (ipp,im,ido);
  //asympotic elastic modulus 
  e_0=Mm->eliso[imat].e;
  e_t = q*e_0;
  
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    Mm->crb3[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    Mm->crrs[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    Mm->crdpl[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }

  return (e_t);
} 



/** 
   function returns actual Young's modulus
    
   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double creep_give_actual_ym (long ipp,long im,long ido) 
{ 
  long i;
  double e_t;

  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    e_t = Mm->crb3[i].give_ym_eqother(ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    e_t = Mm->crrs[i].give_ym_eqother(ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    e_t = Mm->crdpl[i].give_ym_eqother(ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return (e_t);
} 


/** 
   function computes inital Young's modulus    

   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double creep_compute_inital_ym (long ipp,long im,long ido) 
{ 
  long i,imat;
  double q,e_0,e_t;

  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //creep coefficient
  q = creep_matstiffchange (ipp,im,ido);
  //asympotic elastic modulus 
  e_0=Mm->eliso[imat].e;
  e_t = q*e_0;
  
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    Mm->crb3[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    Mm->crrs[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    Mm->crdpl[i].store_ym_eqother(e_t,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }

  return (e_t);
} 


/** 
    function returns actual tensile strenght
    
    @param ipp - integration point
    @param im  - index of material type for given ip
    @param ido - index of internal variables for given material in the ipp other array
    
    TKr, 07/08/2008 - revised
*/ 
double creep_give_actual_ft (long ipp,long im,long ido) 
{ 
  long i;
  double fft;

  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    fft = Mm->crb3[i].creep_give_actual_ft (ipp,im,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    fft = Mm->crrs[i].creep_give_actual_ft (ipp,im,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    fft = Mm->crdpl[i].creep_give_actual_ft (ipp,im,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return (fft);
} 


/** 
   function returns actual compression strenght

   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/ 
double creep_give_actual_fc (long ipp,long im,long ido) 
{ 
  long imat;
  double q,e_0,e_t,ffc;
  
  imat=Mm->ip[ipp].idm[Mm->ip[ipp].gemid()];
  //creep coefficient
  q = creep_matstiffchange (ipp,im,ido);
  //asympotic elastic modulus 
  e_0=Mm->eliso[imat].e;
  e_t = q*e_0;
  
  //not completed??!!
  ffc = e_t/1000.0;

  return (ffc);
} 


/**
   function initializes creep material model and retardation coefficient (stiffnesses of Kelvin chain units)

   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
*/
void creep_initmaterialmodel(long ipp,long im,long ido)
{
  long i,j,k,napproxtime;
  double t_age,tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,q,delYi,delYj,m0,m;
  double jt,ym0,exp_prec=-50;

  //number of retardation times
  long n_ret_times = creep_number_rettimes (ipp,im);
  vector bpom(ASTCKVEC(n_ret_times)), ret_times(ASTCKVEC(n_ret_times)), emu(ASTCKVEC(n_ret_times));
  matrix apom(ASTCKMAT(n_ret_times, n_ret_times)), cpom(ASTCKMAT(n_ret_times, n_ret_times));


  // memory allocation
  fillv(0.0, bpom);
  fillm(0.0, apom);
  fillm(0.0, cpom);
  fillv(0.0, ret_times);
  fillv(0.0, emu);  
  jt=0.0;

  //integration coffecient
  q=2./3.;

  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];    
    //time initialization
    Mm->crb3[i].previoustime = Mp->time;
    //compute initial values
    Mm->crb3[i].initvalues(ipp, ido);
    //times -> ages [days]
    //actuall times in days
    Mm->crb3[i].compute_ages (ipp,ido);
    Mm->crb3[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crb3[i].give_rettimes(ret_times,n_ret_times,ipp);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];    
    //time initialization
    Mm->crrs[i].previoustime = Mp->time;
    //compute initial values
    Mm->crrs[i].initvalues(ipp,im,ido);
    //times -> ages [days]
    //actuall times in days
    Mm->crrs[i].compute_ages (ipp,ido);
    Mm->crrs[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crrs[i].give_rettimes(ret_times,n_ret_times,ipp);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];    
    //time initialization
    Mm->crdpl[i].previoustime = Mp->time;
    //times -> ages [days]
    //actuall times in days
    Mm->crdpl[i].compute_ages (ipp,ido);
    Mm->crdpl[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crdpl[i].give_rettimes(ret_times,n_ret_times,ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  //Dirichlet series
  if (Mm->ip[ipp].tm[im] != creeprs){
    t_age=tl_age+0.0001;
    
    m0 = log10(tl_age+0.0001);
    m = log10(2.0*timeMax);
    k=0;
    //computing actual stiffness from compliance function
    //coverting creep function into Dirichlet-Prony series
    do{
      //creep function
      switch (Mm->ip[ipp].tm[im]){
      case creepb3:{
	i=Mm->ip[ipp].idm[im];
	//actuall compliance function
	jt = Mm->crb3[i].b3_law (th_age,tl_age,t_age,ipp,ido);
	break;
      }
      case creepdpl:{
	i=Mm->ip[ipp].idm[im];
	//actuall compliance function
	jt = Mm->crdpl[i].double_power_law (tl_age,t_age,ipp,ido);
	break;
      }
      default:{
	print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
        abort();
      }
      }
      
      for (i=0;i<n_ret_times;i++){
	delYi=pow((tl_age/ret_times[i]),q)-pow((t_age/ret_times[i]),q);
	if(delYi <= exp_prec)//numerical correction
	  delYi = exp_prec;
	bpom[i]=bpom[i]+(1.-exp(delYi))*jt;
	for (j=0;j<n_ret_times;j++){
	  delYj= pow((tl_age/ret_times[j]),q)-pow((t_age/ret_times[j]),q);
	  if(delYj <= exp_prec)//numerical correction
	    delYj = exp_prec;
	  apom[i][j]=apom[i][j]+(1.-exp(delYi))*(1.-exp(delYj));
	}
      }
      k=k+1;

      t_age=pow(10.,m0+(m-m0)/napproxtime*k);
    }while (t_age <= 2.0*timeMax);

    invm(apom,cpom,1.0e-20);
    
    //computing coefficients of Dirichlet series
    for (i=0;i<n_ret_times;i++){    
      for (j=0;j<n_ret_times;j++)
	emu[i]=emu[i]+cpom[i][j]*bpom[j];
      
      emu[i]=1./emu[i];
    }
  }
  
  //Continuous retardation spectrum
  if (Mm->ip[ipp].tm[im] == creeprs){
    jt = Mm->crrs[i].give_J_E_mu(emu,th_age,tb_age,tb_age_dt,ipp,ido);
  }
  
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //storing stiffness E_mu of Kelvin chains into itegration point
    Mm->crb3[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //storing stiffness E_mu of Kelvin chains into itegration point
    Mm->crrs[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //storing stiffness E_mu of Kelvin chains into itegration point
    Mm->crdpl[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  //inital effective Young's modulus
  ym0 = creep_compute_inital_ym (ipp,im,ido);
}


/**
   function updates creep material model
    
   @param ipp - integration point
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array
   
   TKr, 07/08/2008 - revised
   TKr, 17/09/2014 - new revision
*/
void creep_updateval(long ipp,long im,long ido)
{
  long i,j,napproxtime;
  double mu,tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,e;
  strastrestate ss;
  //  number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(nc)), sigma(ASTCKVEC(nc)), sigmao(ASTCKVEC(nc)), dsigma(ASTCKVEC(nc));

  //number of retardation times
  long n_ret_times = creep_number_rettimes (ipp,im);
  vector ret_times(ASTCKVEC(n_ret_times)), emu(ASTCKVEC(n_ret_times));
  matrix screep(ASTCKMAT(nc, n_ret_times));


  //  stress/strain state; planestress=10,planestrain=11,plate=15,axisymm=20,shell=25,spacestress=30
  ss=Mm->ip[ipp].ssst;  
  //Poisons ratio
  mu=Mm->give_actual_nu(ipp);

  // memory clearing
  fillv(0.0, ret_times);
  fillv(0.0, emu);
  fillm(0.0,screep);
  fillv(0.0,eps);
  fillv(0.0,sigma);
  fillv(0.0,sigmao);
  fillv(0.0,dsigma);

  //updating of strain increments
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    
    // stores total strains into eqother
    for (j=0; j<nc; j++){
      eps[j]=Mm->ip[ipp].strain[j];
    }
    Mm->crb3[i].store_strains_eqother(eps,ipp,ido);   
    //actuall times in days
    Mm->crb3[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //update of material model
    Mm->crb3[i].updatevalues(ipp, ido);
    //ret. times from ip
    Mm->crb3[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crb3[i].give_hidden_strains_eqother(screep,ipp,ido);
    //restoring of previous total overall stresses from eqother
    Mm->crb3[i].give_stresses_eqother(sigmao,ipp,ido);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crb3[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    
    //increment of overall stresses
    for (j=0;j<nc;j++){
      sigma[j] = Mm->ip[ipp].stress[j];
      dsigma[j] = sigma[j] - sigmao[j];
    }
    Mm->crb3[i].store_stresses_eqother(sigma,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    
    // stores total strains into eqother
    for (j=0; j<nc; j++){
      eps[j]=Mm->ip[ipp].strain[j];
    }
    Mm->crrs[i].store_strains_eqother(eps,ipp,ido);
    //actuall times in days
    Mm->crrs[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //update of material model
    Mm->crrs[i].updatevalues(ipp,im,ido);
    //ret. times from ip
    Mm->crrs[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crrs[i].give_hidden_strains_eqother(screep,ipp,ido);
    //restoring of previous total overall stresses from eqother
    Mm->crrs[i].give_stresses_eqother(sigmao,ipp,ido);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crrs[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    
    //increment of overall stresses
    for (j=0;j<nc;j++){
      sigma[j] = Mm->ip[ipp].stress[j];
      dsigma[j] = sigma[j] - sigmao[j];
    }
    Mm->crrs[i].store_stresses_eqother(sigma,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    
    // stores total strains into eqother
    for (j=0; j<nc; j++){
      eps[j]=Mm->ip[ipp].strain[j];
    }
    Mm->crdpl[i].store_strains_eqother(eps,ipp,ido);
    //actuall times in days
    Mm->crdpl[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crdpl[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crdpl[i].give_hidden_strains_eqother(screep,ipp,ido);
    //restoring of previous total overall stresses from eqother
    Mm->crdpl[i].give_stresses_eqother(sigmao,ipp,ido);
    //coefficients of Dirichlet series from ip (Kelvin chain model)
    Mm->crdpl[i].give_emu_eqother(emu,n_ret_times,ipp,ido);
    
    //increment of overall stresses
    for (j=0;j<nc;j++){
      sigma[j] = Mm->ip[ipp].stress[j];
      dsigma[j] = sigma[j] - sigmao[j];
    }
    Mm->crdpl[i].store_stresses_eqother(sigma,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  //computes new hidden strains
  creep_hidden_strains (screep,dsigma,emu,mu,n_ret_times,ret_times,tb_age,dt,ss);
  
  //stores new hidden strains to eqother array and previous young's modulus
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //stores hidden strains to other
    Mm->crb3[i].store_hidden_strains_eqother(screep,ipp,ido);
    e = Mm->crb3[i].give_ym_eqother(ipp,ido);
    Mm->crb3[i].store_ym_old_eqother(e,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //stores hidden strains to other
    Mm->crrs[i].store_hidden_strains_eqother(screep,ipp,ido);
    e = Mm->crrs[i].give_ym_eqother(ipp,ido);
    Mm->crrs[i].store_ym_old_eqother(e,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //stores hidden strains to other
    Mm->crdpl[i].store_hidden_strains_eqother(screep,ipp,ido);
    //stores previous Young's modulus
    e = Mm->crdpl[i].give_ym_eqother(ipp,ido);
    Mm->crdpl[i].store_ym_old_eqother(e,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function computes new increments of irreversible strains and stresses 
   and retardation coefficient (stiffnesses of Kelvin chain units)

   @param ipp - index of integration point
   @param im  - index of
   @param ido - index in array other  

   TKr, 07/08/2008 - revised
   TKr, 17/09/2014 - new revision
*/
void creep_nlstressesincr (long ipp,long im,long ido)
{
  long i,j,k,napproxtime;
  double t_age,tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,q,delYi,delYj,m0,m;
  double jt,ym,inv_v,q4,nonlin_func,dj,mu,param_exp,exp_prec=-50;
  strastrestate ss;
  vector deps_help;
  matrix c;

  // number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector deps(ASTCKVEC(nc)), deps_cr(ASTCKVEC(nc)), deps_sh(ASTCKVEC(nc)), deps_ss(ASTCKVEC(nc));
  vector dsigma(ASTCKVEC(nc)), sigma(ASTCKVEC(nc));
  // number of retardation times
  long n_ret_times = creep_number_rettimes (ipp,im);
  vector bpom(ASTCKVEC(n_ret_times)), ret_times(ASTCKVEC(n_ret_times)), emu(ASTCKVEC(n_ret_times));
  matrix apom(ASTCKMAT(n_ret_times,n_ret_times)), cpom(ASTCKMAT(n_ret_times,n_ret_times)), screep(ASTCKMAT(nc,n_ret_times)), d(ASTCKMAT(nc,nc));

  // stress/strain state; planestress=10,planestrain=11,plate=15,axisymm=20,shell=25,spacestress=30
  ss=Mm->ip[ipp].ssst;  
  // Poisons ratio
  mu=Mm->give_actual_nu(ipp);

  // memory allocation
  fillv(0.0,deps);
  fillv(0.0,deps_cr);
  fillv(0.0,deps_sh);
  fillv(0.0,deps_ss);
  fillv(0.0,dsigma);
  fillv(0.0,sigma);
  fillv(0.0, ret_times);
  fillv(0.0, emu);
  fillm(0.0,screep);
  fillv(0.0, bpom);
  fillm(0.0, apom);
  fillm(0.0, cpom);
  jt = 0.0;

  //integration coffecient
  q=2./3.;

  //overall stresses
  for (j=0;j<nc;j++){
    sigma[j] = Mm->ip[ipp].stress[j];
  }

  //new quantities
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //actual time
    Mm->crb3[i].previoustime = Mp->time-Mp->dtime;//this is necessary
    //new ages
    Mm->crb3[i].compute_ages (ipp,ido);
    Mm->crb3[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crb3[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crb3[i].give_hidden_strains_eqother(screep,ipp,ido);

    //ageing function 1/v(t)
    inv_v = 1.0;
    //flow creep is already in comliance function
    q4 = 0.0;
    nonlin_func = 1.0;

    //non-mechanical loading
    //computes increments of irreversible strains, free shrinkage and stress-induced shrinkage, in actual time step
    Mm->crb3[i].give_deps_free (deps_sh,th_age,tb_age_dt,tb_age,dt,ipp,im,ido);
    Mm->crb3[i].give_deps_stressinduced (deps_ss,th_age,tb_age_dt,tb_age,sigma,ipp,im,ido);
    //stores increments of shrinkage-irreversible strains to other
    Mm->crb3[i].store_irrdstrains_eqother(deps_sh,ipp,ido);
    //stores increments of stress-induced shrinkage-irreversible strains to other
    Mm->crb3[i].store_stressirrdstrains_eqother(deps_ss,ipp,ido);
    // stores humidity
    Mm->crb3[i].store_hum_eqother(ipp,ido);
    // stores temperature
    Mm->crb3[i].store_temp_eqother(ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //actual time
    Mm->crrs[i].previoustime = Mp->time-Mp->dtime;//this is necessary
    //new ages
    Mm->crrs[i].compute_ages (ipp,ido);
    Mm->crrs[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crrs[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crrs[i].give_hidden_strains_eqother(screep,ipp,ido);

    //ageing function 1/v(t)
    inv_v = Mm->crrs[i].give_inv_v(tl_age);//aritmetical mid_interval value
    //flow creep
    q4 = Mm->crrs[i].give_q4();
    nonlin_func = Mm->crrs[i].give_nonlin_func();

    //non-mechanical loading
    //computes increments of irreversible strains, free shrinkage and stress-induced shrinkage, in actual time step
    Mm->crrs[i].give_deps_free (deps_sh,th_age,tb_age_dt,tb_age,dt,ipp,im,ido);
    Mm->crrs[i].give_deps_stressinduced (deps_ss,th_age,tb_age_dt,tb_age,sigma,ipp,im,ido);
    //stores increments of shrinkage-irreversible strains to other
    Mm->crrs[i].store_irrdstrains_eqother(deps_sh,ipp,ido);
    //stores increments of stress-induced shrinkage-irreversible strains to other
    Mm->crrs[i].store_stressirrdstrains_eqother(deps_ss,ipp,ido);
    // stores humidity
    Mm->crrs[i].store_hum_eqother(ipp,ido);
    // stores temperature
    Mm->crrs[i].store_temp_eqother(ipp,ido);    
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //actual time
    Mm->crdpl[i].previoustime = Mp->time-Mp->dtime;//this is necessary
    //new ages
    Mm->crdpl[i].compute_ages (ipp,ido);
    Mm->crdpl[i].give_ages (tb_age_dt,tb_age,tl_age,th_age,dt,timeMax,napproxtime,ipp);
    //ret. times from ip
    Mm->crdpl[i].give_rettimes(ret_times,n_ret_times,ipp);
    //hidden strains from other
    Mm->crdpl[i].give_hidden_strains_eqother(screep,ipp,ido);

    //ageing function 1/v(t)
    inv_v = 1.0;
    //flow creep is already in compliance function
    q4 = 0.0;
    nonlin_func = 1.0;
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  //computing of new retardation coefficients
  //Classical Dirichlet series
  if (Mm->ip[ipp].tm[im] != creeprs){
    t_age=tl_age+0.0001;
    
    m0 = log10(tl_age+0.0001);
    m = log10(2.0*timeMax);
    k=0;
    
    //computing actual stiffness from compliance function
    //coverting creep function into Dirichlet-Prony series
    while (t_age <= 2.0*timeMax){
      
      //creep function
      switch (Mm->ip[ipp].tm[im]){
      case creepb3:{
	i=Mm->ip[ipp].idm[im];
	//actuall compliance function
	jt = Mm->crb3[i].b3_law (th_age,tl_age,t_age,ipp,ido);
	break;
      }
      case creepdpl:{
	i=Mm->ip[ipp].idm[im];
	//actuall compliance function
	jt = Mm->crdpl[i].double_power_law (tl_age,t_age,ipp,ido);
	break;
      }
      default:{
	print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
        abort();
      }
      }
      
      for (i=0;i<n_ret_times;i++){
	delYi=pow((tl_age/ret_times[i]),q)-pow((t_age/ret_times[i]),q);
	if(delYi <= exp_prec)//numerical correction
	  delYi = exp_prec;
	bpom[i]=bpom[i]+(1.-exp(delYi))*jt;
	for (j=0;j<n_ret_times;j++){
	  delYj= pow((tl_age/ret_times[j]),q)-pow((t_age/ret_times[j]),q);
	  if(delYj <= exp_prec)//numerical correction
	    delYj = exp_prec;
	  apom[i][j]=apom[i][j]+(1.-exp(delYi))*(1.-exp(delYj));
	}
      }
      k=k+1;
      
      t_age=pow(10.,m0+(m-m0)/napproxtime*k);
    }
    
    invm(apom,cpom,1.0e-20);
    
    //computing of new coefficients of Dirichlet series
    for (i=0;i<n_ret_times;i++){    
      for (j=0;j<n_ret_times;j++)
	emu[i]=emu[i]+cpom[i][j]*bpom[j];
      
      emu[i]=1./emu[i];
    }
  }
  //Retardation spectrum
  if (Mm->ip[ipp].tm[im] == creeprs){
    jt = Mm->crrs[i].give_J_E_mu(emu,th_age,tb_age,tb_age_dt,ipp,ido);
  }
  
  //storing retardation coefficients (E_mu of Kelvin chains)into integration point
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    Mm->crb3[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    Mm->crrs[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    Mm->crdpl[i].store_emu_eqother(emu,n_ret_times,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  //computes increments of irreversible strains (creep) from previous time step
  for (i=0;i<n_ret_times;i++){
    dj=pow((tb_age + dt)/ret_times[i],q)-pow(tb_age/ret_times[i],q);
    param_exp = -dj;
    if(param_exp < exp_prec)//numerical correction
      param_exp = exp_prec;
    for (j=0; j<nc; j++){
      //deps_cr[j]=deps_cr[j]+(1.-exp(-dj))*screep[j][i];//old
      //new: added effect of non-linear creep and aging function for continuous retardation spectrum
      deps_cr[j] = deps_cr[j] + nonlin_func*((1.-exp(param_exp))*(screep[j][i]*inv_v));
    }
  }
  
  //added flow term from continuous retardation spectrum
  if (Mm->ip[ipp].tm[im] == creeprs){
    reallocv(RSTCKVEC(nc,deps_help));
    fillv(0.0,deps_help);
    reallocm(RSTCKMAT(nc,nc,c));
    fillm(0.0,c);

    //commented
    //reallocv(RSTCKVEC(nc,sigma));
    //fillv(0.0,sigma);

    //overall stresses
    //for (j=0;j<nc;j++){
    //sigma[j] = Mm->ip[ipp].stress[j];
    //}

    unit_compl_matrix(c,mu,ss);
    mxv (c,sigma,deps_help);
    if (Mm->ip[ipp].ssst == planestress)
      {
	deps_help[3] = -mu / (1.0 - mu) * (deps_help[0]+deps_help[1]);
      }	
    
    for (j=0; j<nc; j++){
      deps_cr[j] = deps_cr[j] + nonlin_func*q4*dt/tl_age*deps_help[j];
    }
  }

  //computes increments of creep strains plus shrinkage and stress-induced shrinkage  
  for(i=0;i<nc;i++){
    deps[i]=deps_cr[i]+deps_sh[i]+deps_ss[i];
  }

  //computes actual Young's modulus
  ym = creep_compute_actual_ym (ipp,im,ido);
  //computes incremets of stresses from irreversible strain increments
  fillm(0.0,d);
  creep_matstiff (d, ipp, im, ido);
  mxv (d,deps,dsigma);

  //storing of actual increments of overall stresses into other
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //irreversible strains added
    Mm->crb3[i].addirrstrains_eqother (deps,ipp,ido);
    //stores increments of creep-irreversible strains to other
    Mm->crb3[i].store_creepdstrains_eqother(deps_cr,ipp,ido);
    Mm->crb3[i].store_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //irreversible strains added
    Mm->crrs[i].addirrstrains_eqother (deps,ipp,ido);
    //stores increments of creep-irreversible strains to other
    Mm->crrs[i].store_creepdstrains_eqother(deps_cr,ipp,ido);
    Mm->crrs[i].store_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //irreversible strains added
    Mm->crdpl[i].addirrstrains_eqother (deps,ipp,ido);
    //stores increments of creep-irreversible strains to other
    Mm->crdpl[i].store_creepdstrains_eqother(deps_cr,ipp,ido);
    Mm->crdpl[i].store_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function returns new increments of irreversible strains and stresses

   @param ipp   - integration point
   @param im    - index of material type for given ip
   @param ido   - index of internal variables for given material in the ipp other array
   @param fi  - first index of the required stress increment component
   @param sigma - %vector of new increments of stresses from irreversible strain increments

   TKr, 07/08/2008 - revised 
*/
void creep_givestressincr (long ipp,long im,long ido,long fi,vector &sig)
{
  long i,j;
  //  number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector dsigma(ASTCKVEC(nc));

  // memory clearing
  fillv(0.0,dsigma);
 
  //actual increments of overall stresses from eqother
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    Mm->crb3[i].give_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    Mm->crrs[i].give_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    Mm->crdpl[i].give_dstresses_eqother(dsigma,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }

  for (i=fi,j=0; j<sig.n; i++, j++){
    sig[j]=dsigma[i];
  }
}
     

/**
   function computes new increments of total stresses

   @param ipp    - integration point
   @param im     - index of material type for given ip
   @param ido    - index of internal variables for given material in the ipp other array
   @param dsigma - %vector of new increments of stresses from irreversible strain increments

   TKr, 07/08/2008 - revised 
*/
void creep_incrtotstresses (long ipp,long im,long ido,vector &dsigma)
{
  long i,j;
  double nu;

  //  number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector eps_old(ASTCKVEC(nc)), deps(ASTCKVEC(nc)), depstot(ASTCKVEC(nc)), deps_cr(ASTCKVEC(nc)), deps_sh(ASTCKVEC(nc)), deps_ss(ASTCKVEC(nc));
  matrix d(ASTCKMAT(nc,nc));

  // memory clearing
  fillv(0.0,eps_old);
  fillv(0.0,deps);
  fillv(0.0,depstot);
  fillv(0.0,deps_cr);
  fillv(0.0,deps_sh);
  fillv(0.0,deps_ss);
  
  //quantities from eqother array
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from - previous time step
    Mm->crb3[i].give_strains_eqother(eps_old,ipp,ido);
    //increments of creep-irreversible strains from other
    Mm->crb3[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    //increments of shrinkage-irreversible strains from other
    Mm->crb3[i].give_irrdstrains_eqother(deps_sh,ipp,ido);
    //increments of stress-induced shrinkage-irreversible strains from other
    Mm->crb3[i].give_stressirrdstrains_eqother(deps_ss,ipp,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from - previous time step
    Mm->crrs[i].give_strains_eqother(eps_old,ipp,ido);
    //increments of creep-irreversible strains from other
    Mm->crrs[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    //increments of shrinkage-irreversible strains from other
    Mm->crrs[i].give_irrdstrains_eqother(deps_sh,ipp,ido);
    //increments of stress-induced shrinkage-irreversible strains from other
    Mm->crrs[i].give_stressirrdstrains_eqother(deps_ss,ipp,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from - previous time step
    Mm->crdpl[i].give_strains_eqother(eps_old,ipp,ido);
    //increments of creep-irreversible strains from other
    Mm->crdpl[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  //computes increments of total reversible strains from previous time step
  nu = Mm->give_actual_nu(ipp);
  if (Mm->ip[ipp].ssst == planestress)
    {
      Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (Mm->ip[ipp].strain[0]+Mm->ip[ipp].strain[1]);
    }
  for (j=0; j<nc; j++){
    depstot[j] = Mm->ip[ipp].strain[j] - eps_old[j];
    deps[j] = depstot[j] - deps_cr[j] - deps_sh[j] - deps_ss[j];
  }

  //  stiffness matrix of material
  fillm(0.0,d);
  creep_matstiff (d, ipp, im, ido);
  //computes new increments of stresses
  mxv (d,deps,dsigma);
}


/**
   function computes new increments of aeging strains

   @param ipp    - integration point
   @param im     - index of material type for given ip
   @param ido    - index of internal variables for given material in the ipp other array

   TKr, 12/03/2009 
*/
void creep_aeging_strains (long ipp,long im,long ido)
{
  long i,j;  
  double ym_old,ym_act;

  //  number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector eps_old(ASTCKVEC(nc)), eps_ag_old(ASTCKVEC(nc)), eps_ir_old(ASTCKVEC(nc)), deps_cr(ASTCKVEC(nc)), deps_sh(ASTCKVEC(nc)), deps_ss(ASTCKVEC(nc));
  
  // memory allocation
  fillv(0.0,eps_old);
  fillv(0.0,eps_ir_old);
  fillv(0.0,deps_cr);
  fillv(0.0,deps_sh);
  fillv(0.0,deps_ss);
  fillv(0.0,eps_ag_old);


  //restoring from eqother
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from eqother - previous time step
    Mm->crb3[i].give_strains_eqother(eps_old,ipp,ido);
    //actual irreversible strains
    Mm->crb3[i].giveirrstrains_eqother(ipp,ido,eps_ir_old);
    //actual creep dstrains
    Mm->crb3[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    //actual shrinkage-irreversible and stress-induced shrinkage-irreversible strains
    Mm->crb3[i].give_irrdstrains_eqother(deps_sh,ipp,ido);
    Mm->crb3[i].give_stressirrdstrains_eqother(deps_ss,ipp,ido);
    for (j=0; j<nc; j++)
      eps_ir_old[j] = eps_ir_old[j] - deps_cr[j] - deps_sh[j] - deps_ss[j];
    //aeging strains from prevous time step 
    Mm->crb3[i].give_agstrains_eqother(ipp,ido,eps_ag_old);
    // previous Young's modulus from eqother - previous time step
    ym_old = Mm->crb3[i].give_ym_old_eqother(ipp,ido);
    // actual Young's modulus
    ym_act = creep_give_actual_ym (ipp,im,ido);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from eqother - previous time step
    Mm->crrs[i].give_strains_eqother(eps_old,ipp,ido);
    //actual irreversible strains
    Mm->crrs[i].giveirrstrains_eqother(ipp,ido,eps_ir_old);
    //actual creep dstrains
    Mm->crrs[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    //actual shrinkage-irreversible and stress-induced shrinkage-irreversible strains
    Mm->crrs[i].give_irrdstrains_eqother(deps_sh,ipp,ido);
    Mm->crrs[i].give_stressirrdstrains_eqother(deps_ss,ipp,ido);
    for (j=0; j<nc; j++)
      eps_ir_old[j] = eps_ir_old[j] - deps_cr[j] - deps_sh[j] - deps_ss[j];
    //aeging strains from prevous time step 
    Mm->crrs[i].give_agstrains_eqother(eps_ag_old,ipp,ido);
    // previous Young's modulus from eqother - previous time step
    ym_old = Mm->crrs[i].give_ym_old_eqother(ipp,ido);
    // actual Young's modulus
    ym_act = creep_give_actual_ym (ipp,im,ido);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];    
    //total strains from eqother - previous time step
    Mm->crdpl[i].give_strains_eqother(eps_old,ipp,ido);
    //actual irreversible strains
    Mm->crdpl[i].giveirrstrains_eqother(ipp,ido,eps_ir_old);
    //actual creep dstrains
    Mm->crdpl[i].give_creepdstrains_eqother(deps_cr,ipp,ido);
    //actual shrinkage-irreversible and stress-induced shrinkage-irreversible strains are not defined in dpl
    for (j=0; j<nc; j++)
      eps_ir_old[j] = eps_ir_old[j] - deps_cr[j];
    //aeging strains from prevous time step 
    Mm->crdpl[i].give_agstrains_eqother(ipp,ido,eps_ag_old);
    // previous Young's modulus from eqother - previous time step
    ym_old = Mm->crdpl[i].give_ym_old_eqother(ipp,ido);
    // actual Young's modulus
    ym_act = creep_give_actual_ym (ipp,im,ido);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  for (j=0; j<nc; j++){
    eps_ag_old[j] = (1.0 - ym_old/ym_act)*(eps_old[j] - eps_ir_old[j]) + ym_old/ym_act*eps_ag_old[j];
  }
  
  //storing into eqother
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    Mm->crb3[i].store_agstrains_eqother(ipp,ido,eps_ag_old);
    break;
  }
  case creeprs:{
    Mm->crrs[i].store_agstrains_eqother(eps_ag_old,ipp,ido);
    break;
  }
  case creepdpl:{
    Mm->crdpl[i].store_agstrains_eqother(ipp,ido,eps_ag_old);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function computes new total stresses

   @param ipp - index of integration point
   @param im  - index of
   @param ido - index in array other  

   TKr 07/08/2008 - revised
*/
void creep_nlstresses (long ipp,long im,long ido)
{
  long i,j;

  //  number of stres components on element
  long nc=Mm->ip[ipp].ncompstr;
  vector sigma(ASTCKVEC(nc)), sigmao(ASTCKVEC(nc)), dsigma(ASTCKVEC(nc)), eps(ASTCKVEC(nc)), epsag(ASTCKVEC(nc)), epscr(ASTCKVEC(nc));


  // memory allocation
  fillv(0.0,sigma);
  fillv(0.0,sigmao);
  fillv(0.0,dsigma);
  fillv(0.0,eps);
  fillv(0.0,epsag);
  fillv(0.0,epscr);

  creep_incrtotstresses(ipp, im, ido, dsigma);

  //stores new quantities to eqother array
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    //old stresses from eqother
    Mm->crb3[i].give_stresses_eqother(sigmao,ipp,ido);
    //new total stresses
    for (j=0; j<nc; j++)
      sigma[j] = sigmao[j] + dsigma[j];
    //store stresses into ip
    for (j=0;j<nc;j++){
      Mm->ip[ipp].stress[j] = sigma[j];
    }
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    //old stresses from eqother
    Mm->crrs[i].give_stresses_eqother(sigmao,ipp,ido);
    //new total stresses
    for (j=0; j<nc; j++)
      sigma[j] = sigmao[j] + dsigma[j];
    //store stresses into ip
    for (j=0;j<nc;j++){
      Mm->ip[ipp].stress[j] = sigma[j];
    }
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    //old stresses from eqother
    Mm->crdpl[i].give_stresses_eqother(sigmao,ipp,ido);
    //new total stresses
    for (j=0; j<nc; j++)
      sigma[j] = sigmao[j] + dsigma[j];
    //store stresses into ip
    for (j=0;j<nc;j++){
      Mm->ip[ipp].stress[j] = sigma[j];
    }
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function computes hidden strains
   
   @param screep - %vector of hidden strains
   @param sig    - %vector of stress increments
   @param ipp    - index of integration point
   @param im     - index of
   @param ido    - index in array other  

   TKr 07/08/2008 - revised
*/
void creep_hidden_strains (matrix &screep,vector &sig,vector &emu, double mi, long n_ret_times, vector &ret_times, double actualtime, double dt,strastrestate ss)
{
  int i;
  double q,dj,pomt,pp;
  double param_exp,exp_prec=-50;
  q=2./3.;
  if((actualtime+dt) == dt) {
  }
  else{
    pomt=actualtime+dt;
    
    if(ss==bar){
      for (i=0;i<n_ret_times;i++){
	dj=pow(pomt/ret_times[i],q)-pow(actualtime/ret_times[i],q);
	param_exp = -dj;
	if(param_exp < exp_prec)//numerical correction
	  param_exp = exp_prec;
	screep[0][i]=exp(param_exp)*screep[0][i]+sig[0]*(1.-exp(param_exp))/dj/emu[i];
      }
    }
    else if(ss==planestress){
      for (i=0;i<n_ret_times;i++){
	dj=pow(pomt/ret_times[i],q)-pow(actualtime/ret_times[i],q);
	param_exp = -dj;
	if(param_exp < exp_prec)//numerical correction
	  param_exp = exp_prec;
	screep[0][i]=exp(param_exp)*screep[0][i]+(sig[0]-mi*sig[1])*(1.-exp(param_exp))/dj/emu[i];
	screep[1][i]=exp(param_exp)*screep[1][i]+(sig[1]-mi*sig[0])*(1.-exp(param_exp))/dj/emu[i];
	screep[2][i]=exp(param_exp)*screep[2][i]+ sig[2]*2.*(1.+mi)*(1.-exp(param_exp))/dj/emu[i];
	screep[3][i]=exp(param_exp)*screep[3][i]+(-mi*sig[0]-mi*sig[1])*(1.-exp(param_exp))/dj/emu[i];
      }
    }
    else if(ss==planestrain){
      pp=1.+mi;
      for (i=0;i<n_ret_times;i++){
	dj=pow(pomt/ret_times[i],q)-pow(actualtime/ret_times[i],q);
	param_exp = -dj;
	if(param_exp < exp_prec)//numerical correction
	  param_exp = exp_prec;
	screep[0][i]=exp(param_exp)*screep[0][i]+ pp*(sig[0]*(1.-mi)-mi*sig[1])*(1.-exp(param_exp))/dj/emu[i];
	screep[1][i]=exp(param_exp)*screep[1][i]+ pp*(sig[1]*(1.-mi)-mi*sig[0])*(1.-exp(param_exp))/dj/emu[i];
	screep[2][i]=exp(param_exp)*screep[2][i]+ pp*sig[2]*2.   *(1.-exp(param_exp))/dj/emu[i];
	screep[3][i]=0.0;
      }
    }
    else if(ss==axisymm){
      for (i=0;i<n_ret_times;i++){
	dj=pow(pomt/ret_times[i],q)-pow(actualtime/ret_times[i],q);
	param_exp = -dj;
	if(param_exp < exp_prec)//numerical correction
	  param_exp = exp_prec;
	screep[0][i]=exp(param_exp)*screep[0][i]+(sig[0]-mi*sig[1]-mi*sig[2])*(1.-exp(param_exp))/dj/emu[i];
	screep[1][i]=exp(param_exp)*screep[1][i]+(sig[1]-mi*sig[0]-mi*sig[2])*(1.-exp(param_exp))/dj/emu[i];
	screep[2][i]=exp(param_exp)*screep[2][i]+(sig[2]-mi*sig[0]-mi*sig[1])*(1.-exp(param_exp))/dj/emu[i];
	screep[3][i]=exp(param_exp)*screep[3][i]+ (sig[3]*2.0*(1.0+mi))*(1.-exp(param_exp))/dj/emu[i];//??!!
      }
    }
    else if(ss==spacestress){
      for (i=0;i<n_ret_times;i++){
	dj=pow(pomt/ret_times[i],q)-pow(actualtime/ret_times[i],q);
	param_exp = -dj;
	if(param_exp < exp_prec)//numerical correction
	  param_exp = exp_prec;
	screep[0][i]=exp(param_exp)*screep[0][i]+(sig[0]-mi*sig[1]-mi*sig[2])*(1.-exp(param_exp))/dj/emu[i];
	screep[1][i]=exp(param_exp)*screep[1][i]+(sig[1]-mi*sig[0]-mi*sig[2])*(1.-exp(param_exp))/dj/emu[i];
	screep[2][i]=exp(param_exp)*screep[2][i]+(sig[2]-mi*sig[0]-mi*sig[1])*(1.-exp(param_exp))/dj/emu[i];
	screep[3][i]=exp(param_exp)*screep[3][i]+ sig[3]*(1.+mi)   *(1.-exp(param_exp))/dj/emu[i];
	screep[4][i]=exp(param_exp)*screep[4][i]+ sig[4]*(1.+mi)   *(1.-exp(param_exp))/dj/emu[i];
	screep[5][i]=exp(param_exp)*screep[5][i]+ sig[5]*(1.+mi)   *(1.-exp(param_exp))/dj/emu[i];
      }
    }
  }
}


/**
   function returns total irreversible strains (aeging is included)

   @param ipp   - integration point
   @param im    - index of material type for given ip
   @param ido   - index of internal variables for given material in the ipp other array
   @param epsir - %vector of irreversible strains

   TKr, 07/08/2008 - revised 

*/
void creep_giveirrstrains (long ipp, long im, long ido, vector &epsir)
{
  long i;
  vector epsag(ASTCKVEC(epsir.n));

  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    Mm->crb3[i].giveirrstrains_eqother (ipp,ido,epsir);
    //added aeging strains due to new Young's modulus
    //computing of total "aeging" strains
    creep_aeging_strains (ipp,im,ido);    
    Mm->crb3[i].give_agstrains_eqother(ipp,ido,epsag);
    addv(epsir,epsag,epsir);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    Mm->crrs[i].giveirrstrains_eqother (ipp,ido,epsir);
    //added aeging strains due to new Young's modulus
    //computing of total "aeging" strains
    creep_aeging_strains (ipp,im,ido);    
    Mm->crrs[i].give_agstrains_eqother(epsag,ipp,ido);
    addv(epsir,epsag,epsir);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    Mm->crdpl[i].giveirrstrains_eqother (ipp,ido,epsir);
    //added aeging strains due to new Young's modulus
    //computing of total "aeging" strains
    creep_aeging_strains (ipp,im,ido);    
    Mm->crdpl[i].give_agstrains_eqother(ipp,ido,epsag);
    addv(epsir,epsag,epsir);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function returns number of retardation times

   @param ipp    - integration point
   @param im     - index of material type for given ip

   TKr, 07/08/2008 - revised 
*/
long creep_number_rettimes (long ipp,long im)
{
  long i,n_ret_times;

  //creep function
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    
    //number of retardation times
    n_ret_times = Mm->crb3[i].give_nret_time();
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    
    //number of retardation times
    n_ret_times = Mm->crrs[i].give_nret_time();
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    
    //number of retardation times
    n_ret_times = Mm->crdpl[i].give_nret_time();
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(n_ret_times);
}


/**
  function returns number of eqother components

   @param ipp    - integration point
   @param im     - index of material type for given ip

   TKr, 07/08/2008 - revised 
*/
long creep_ncompo (long ipp,long im)
{
  long i,ncompo;

  //creep function
  switch (Mm->ip[ipp].tm[im]){
  case creepb3:{
    i=Mm->ip[ipp].idm[im];
    
    ncompo = Mm->crb3[i].give_nceqother (ipp);
    break;
  }
  case creeprs:{
    i=Mm->ip[ipp].idm[im];
    
    ncompo = Mm->crrs[i].give_nceqother (ipp);
    break;
  }
  case creepdpl:{
    i=Mm->ip[ipp].idm[im];
    
    ncompo = Mm->crdpl[i].give_nceqother (ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type required in ", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(ncompo);
}


/**
   function assembles compliance %matrix of material
   
   @param c - unit compliance %matrix of material
   @param ssst - strain/stress state
*/
void unit_compl_matrix (matrix &c,double nu,strastrestate ssst)
{
  fillm(0.0,c);
  
  switch (ssst){
  case bar:{
    c[0][0] = 1.0;
    break;
  }
  case plbeam:{
    c[0][0] = 1.0;
    c[1][1] = 1/(1/2.0/(1.0+nu));
    c[2][2] = 1.0;
    break;
  }
  case spacebeam:{
    c[0][0] = 1;
    c[1][1] = 1.0/(1/2.0/(1.0+nu));
    c[2][2] = 1.0/(1/2.0/(1.0+nu));
    c[3][3] = 1.0/(1/2.0/(1.0+nu));
    c[4][4] = 1;
    c[5][5] = 1;
    break;
  }
  case planestress:{
    double g;
    
    g = 1.0/(2*(1.0+nu))*2.0*(1+nu);
    
    c[0][0] = g;     c[0][1] = -nu*g; c[0][2] = 0.0;
    c[1][0] = -nu*g; c[1][1] = g;     c[1][2] = 0.0;
    c[2][0] = 0.0;   c[2][1] = 0.0;   c[2][2] = 1;

    break;
  }
  case planestrain:{
    double g;
    
    g = (1.0 - nu)/2.0*2.0*(1+nu);
    
    c[0][0] = g;       c[0][1] = -nu/2.0; c[0][2] = 0.0;
    c[1][0] = -nu/2.0; c[1][1] = g;       c[1][2] = 0.0;
    c[2][0] = 0.0;     c[2][1] = 0.0;     c[2][2] = 1;

    break;
  }
    
  case axisymm:{
    double g;
    g = 2.0*(1.0+nu);
    
    c[0][0]=1.0;       c[0][1]=-nu;        c[0][2]=c[0][1];    c[0][3]=0.0;
    c[1][0]=c[0][1];   c[1][1]=c[0][0];    c[1][2]=c[0][1];    c[1][3]=0.0;
    c[2][0]=c[0][1];   c[2][1]=c[0][1];    c[2][2]=c[0][0];    c[2][3]=0.0;
    c[3][0]=c[0][3];   c[3][1]=c[1][3];    c[3][2]=c[2][3];    c[3][3]=g;
    break;
  }

  case spacestress:{
    c[0][0]=1;         c[0][1]=-nu;      c[0][2]=-nu;
    c[1][0]=c[0][1];   c[1][1]=c[0][0];  c[1][2]=c[0][1];
    c[2][0]=c[0][1];   c[2][1]=c[0][1];  c[2][2]=c[0][0];
    c[3][3]=2.0*(1.0+nu);
    c[4][4]=2.0*(1.0+nu);
    c[5][5]=2.0*(1.0+nu);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}
