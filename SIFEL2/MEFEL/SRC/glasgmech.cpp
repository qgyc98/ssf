#include <math.h>
#include "glasgmech.h"
#include "vector.h"
#include "matrix.h"
#include "vecttens.h"
#include "tensor.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "elastisomat.h"


#define nijac 20
#define limit 1.0e-8

glasgmech::glasgmech (void)
{
  tkappa0=0.0;
  st=0.0;
  k=0.0;
  gf0=0.0;
  lc=0.0;
  a=0.0;
  b=0.0;
  c=0.0;
  ft = paramf_type(0);
}

glasgmech::~glasgmech (void)
{

}

/**
   function reads input data from the input text file
*/
void glasgmech::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %m", &tkappa0, &st, &gf0, &lc, &paramf_type_kwdset, (int *)(&ft));
  xfscanf(in, "%lf", &k);
  xfscanf (in,"%lf %lf %lf",&a,&b,&c);
}

/**
   function computes coefficient beta
   
   @param ipp - number of integration point
   @param t - temperature
   
   13.7.2004
*/
double glasgmech::betacoeff (long ipp,double t)
{
  double tstar,tbar,beta;
  
  tstar=(470.0-20.0)/100.0;
  tbar=(t-20.0)/100.0;
  
  if (0.0 > tbar){
    fprintf (stderr,"\n\n Warning: normalized temperature t_bar is less than 0.0 at integration point %ld",ipp);
    return 0.0;
  }
  if (0.0<=tbar && tbar<=tstar){
    beta=(2.0*a*tbar+b)*0.01;
  }
  else{
    beta=(2.0*c*(tbar-tstar)+2.0*a*tstar+b)*0.01;
  }
  
  return beta;
}


/**
   function computes norm of equivalent strains
   
   @param ipp - number of integration point
   @param  eps - %vector of the strains
   @param  kappa - %vector where the computed parameters are stored

*/
void glasgmech::damfuncpar(long /*ipp*/, vector &/*eps*/, vector &/*kappa*/)
{
  /*
  vector poseps;
  vector princeps;
  matrix pvect;
  matrix epst;
  matrix dev;
  double i1e, j2e, nu, e, tmp, a, b;
  long ncomp, idem;

  switch (ft)
  {
    case vonmises :
      ncomp=Mm->ip[ipp].ncompstr;
      reallocm(RSTCKMAT(3, 3, epst));
      reallocm(RSTCKMAT(3, 3, dev));
      vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
      i1e = first_invar(epst);
      deviator (epst,dev);
      j2e = second_invar (dev);
      idem = Mm->ip[ipp].gemid();
      if (Mm->ip[ipp].tm[idem] != elisomat)
      {
        fprintf(stderr, "\n\n Invalid type of elastic material is required");
        fprintf(stderr, "\n  in function glasgmech::damfuncpar, (file %s, line %d)\n", __FILE__, __LINE__);
      }
      e  = Mm->eliso[idem].e;
      nu = Mm->eliso[idem].nu;
      a   = (k-1.0)/(2.0*k)/(1.0-2.0*nu);
      b   = 3.0/k/(1+nu);
      tmp = a*a*i1e*i1e + b*j2e;
      kappa[0] = a*i1e+sqrt(tmp);
      break;
    case normazar :
      reallocm(RSTCKMAT(3, 3, epst));
      reallocm(RSTCKMAT(3, 3, pvect));
      reallocv(RSTCKVEC(3, princeps));
      reallocv(RSTCKVEC(3,poseps));
      nullv(princeps);
      vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
      princ_val (epst,princeps,pvect,nijac,limit,Mp->zero,3,1);
      extractposv (princeps, poseps);
      scprd(poseps, poseps, tmp);
      kappa[0] = sqrt(tmp);
      break;
    default :
    {
      fprintf(stderr, "\n\n Unknown damage parameter function type");
      fprintf(stderr, "\n  in function glasgmech::damfuncpar, (file %s, line %d)\n", __FILE__, __LINE__);
    }
  }
  return;
  */
}

/**
  This function computes damage parameter omega which is the result of the
  damage function.

  @param ipp - integration point number
  @param tempr - actual temperature
  @param chi - thermal damage parameter
  @param kappa - %vector of the parameters of damage function
  
*/
double glasgmech::damfunction(long ipp,double tempr,double chi,vector &kappa)
{
  double omega = 0.0, e0, kappa0, tn, gamma, gf,th;
  long idem = Mm->ip[ipp].gemid();
  
  if (Mm->ip[ipp].tm[idem] != elisomat)
  {
    fprintf(stderr, "\n\nError - ip %ld has wrong type of elastic material combined with", ipp);
    fprintf(stderr, "\n scalar damge. The elastic isotropic material is requried");
    fprintf(stderr, "\n Function glasgmech::damfunction (file %s, line %d)\n", __FILE__, __LINE__);
    return 0.0;
  }

  
  // Young modulus
  e0 = Mm->eliso[Mm->ip[ipp].idm[idem]].e;
  
  //  normalized temperature
  tn = (tempr-20.0)/100.0;
  //  maximum reached normalized temperature
  th = (kappa[1]-20.0)/100.0;
  if (tn>th){
    th=tn;
  }
  //  check of maximum reached normalized temperature
  if ((th < 0.0) || (th > 7.9)){
    fprintf(stderr, "\n\nWarning - largest value of normalised temperature\n");
    fprintf(stderr, "          is out of suggested range 0.0-7.9 (ipp %ld)\n",ipp);
  }
  
  //  threshold of the mechanical damage (depends on temperature)
  kappa0 = st/e0*(1.0 - 0.016*th*th)/(1.0 - 0.1*th)/(1.0 - 0.1*th);
  
  //  fracture energy
  gf = gf0*(1.0+0.39*th-0.07*th*th);
  
  //  softening parameter
  gamma = (1.0-chi)*st*lc/gf;
  
  if (kappa[0] < kappa0)
    // elastic loading
    return 0.0;
  
  //  mechanical damage parameter
  omega = 1.0-(kappa0/kappa[0])*exp(-gamma*(kappa[0]-kappa0));

  return omega;
}

/**
  Function computes derivative of mechanical damage parameter with respect to equivalent strain

  @param ipp - integration point number
  @param tempr - actual temperature
  @param chi - thermal damage parameter
  @param kappa - %vector of the parameters of damage function

*/
double glasgmech::domegadkmd (long ipp,double tempr,double chi, vector &kappa)
{
  double domega = 0.0, e0, kappa0, tn, gamma, gf,th;
  long idem = Mm->ip[ipp].gemid();
  
  if (Mm->ip[ipp].tm[idem] != elisomat)
  {
    fprintf(stderr, "\n\nError - ip %ld has wrong type of elastic material combined with", ipp);
    fprintf(stderr, "\n scalar damge. The elastic isotropic material is requried");
    fprintf(stderr, "\n Function glasgmech::damfunction (file %s, line %d)\n", __FILE__, __LINE__);
    return 0.0;
  }

  
  // Young modulus
  e0 = Mm->eliso[Mm->ip[ipp].idm[idem]].e;
  //  normalized temperature
  tn = (tempr-20.0)/100.0;
  //  maximum reached normalized temperature
  th = (kappa[1]-20.0)/100.0;
  if (tn>th){
    th=tn;
  }
  //  check of maximum reached normalized temperature
  if ((th < 0.0) || (th > 7.9)){
    fprintf(stderr, "\n\nWarning - largest value of normalised temperature\n");
    fprintf(stderr, "          is out of suggested range 0.0-7.9 (ipp %ld)\n",ipp);
  }
  
  //  threshold of the mechanical damage (depends on temperature)
  kappa0 = st/e0*(1.0 - 0.016*th*th)/(1.0 - 0.1*th)/(1.0 - 0.1*th);
  
  //  fracture energy
  gf = gf0*(1.0+0.39*th-0.07*th*th);
  
  //  softening parameter
  gamma = (1.0-chi)*st*lc/gf;

  if (kappa[0] < kappa0)
    // elastic loading
    return 0.0;
  
  domega  = -kappa0*exp(-gamma*(kappa[0]-kappa0))/kappa[0];
  domega += kappa0*kappa[0]*gamma*exp(-gamma*(kappa[0]-kappa0));
  domega /= kappa[0]*kappa[0];
  
  return domega;
}


/**
  Function computes derivative of mechanical damage parameter with respect to temperature

  @param ipp - integration point number
  @param tempr - actual temperature
  @param chi - thermal damage parameter
  @param kappa - %vector of the parameters of damage function

*/
double glasgmech::domegadt (long ipp,double tempr,double chi, vector &kappa)
{
  double domega = 0.0, e0, kappa0, tn, gamma, gf,th, dkmd0dt, dgdt;
  long idem = Mm->ip[ipp].gemid();
  
  if (Mm->ip[ipp].tm[idem] != elisomat)
  {
    fprintf(stderr, "\n\nError - ip %ld has wrong type of elastic material combined with", ipp);
    fprintf(stderr, "\n scalar damge. The elastic isotropic material is requried");
    fprintf(stderr, "\n Function glasgmech::damfunction (file %s, line %d)\n", __FILE__, __LINE__);
    return 0.0;
  }

  
  // Young modulus
  e0 = Mm->eliso[Mm->ip[ipp].idm[idem]].e;
  //  normalized temperature
  tn = (tempr-20.0)/100.0;
  //  maximum reached normalized temperature
  th = (kappa[1]-20.0)/100.0;
  if (tn>th){
    th=tn;
  }
  //  check of maximum reached normalized temperature
  if ((th < 0.0) || (th > 7.9)){
    fprintf(stderr, "\n\nWarning - largest value of normalised temperature\n");
    fprintf(stderr, "          is out of suggested range 0.0-7.9 (ipp %ld)\n",ipp);
  }
  
  //  threshold of the mechanical damage (depends on temperature)
  kappa0 = st/e0*(1.0 - 0.016*th*th)/(1.0 - 0.1*th)/(1.0 - 0.1*th);
  
  //  fracture energy
  gf = gf0*(1.0+0.39*th-0.07*th*th);
  
  //  softening parameter
  gamma = (1.0-chi)*st*lc/gf;

  if (kappa[0] < kappa0)
    // elastic loading
    return 0.0;
  
  // derivative of kappa_md with respect of temperature
  dkmd0dt = st/e0/100.0*(-0.032*th*(1.0-0.1*th)*(1.0-0.1*th)+(1.0-0.016*th*th)*2.0*(1.0-0.1*th)*0.1);
  dkmd0dt /= (1.0-0.1*th)*(1.0-0.1*th)*(1.0-0.1*th)*(1.0-0.1*th);
  // derivative of gamma with respect of temperature
  dgdt = -(1-chi)*st*lc/(100.0*gf*gf)*(0.39-0.14*th)*gf0;

  domega  = dkmd0dt*exp(-gamma*(kappa[0]-kappa0))/kappa[0];
  domega += kappa0/kappa[0]*exp(-gamma*(kappa[0]-kappa0))*(-(dgdt*(kappa[0]-kappa0)-gamma*dkmd0dt));
  domega  = -domega;
  
  return domega;
}

/**
  This function computes material stiffnes %matrix.

  @param d[out] - allocated matrix structure for material stiffness %matrix
  @param ipp[in] - integration point number
  @param ido[in] - index of internal variables for given material in the ipp other array

*/
void glasgmech::matstiff (matrix &d, long ipp, long ido)
{
  double dp;

  switch (Mp->nlman->stmat)
  {
    case initial_stiff :
      Mm->elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff :
      Mm->elmatstiff(d, ipp, ido);
      dp=Mm->ip[ipp].eqother[ido+1];
      if (dp > 0.999999)
        dp = 0.999999;
      cmulm (1.0-dp,d);
      break;
    default :
      fprintf(stderr, "\n\nError - unknown type of stifness matrix");
      fprintf(stderr, "\n in function glasgmech::matstiff (%s, line %d)\n", __FILE__, __LINE__);
  }
}

/**
   function computes 
   
   @param t - temeprature
   @param dt - incrment of temperature
   @param eps - thermal strain increment
   
*/
void glasgmech::compute_thermdilat (double t,double dt,matrix &eps)
{
  double auxt,alpha;
  
  auxt = (t-20.0)/100.0;
  if (auxt>=0.0 && auxt<=6.0){
    alpha = 6.0e-5/(7.0-auxt);
  }
  else{  alpha=0.0; }
  
  fillm (0.0,eps);
  
  eps[0][0]=alpha*dt;
  eps[1][1]=alpha*dt;
  eps[2][2]=alpha*dt;
  
}

/**
  This function computes thermal damage parameter chi which is the result of the
  thermal damage function.

  @param ipp - integration point number
  @param tempr - actual temperature
  @param kappa - %vector of the parameters of thermal damage function
                 it contains the maximum of either the largest value attained by temperature
		 or the reference temperature

*/
double glasgmech::thermdamfunction (long /*ipp*/,double tempr,vector &kappa)
{
  double chi = 0.0;
  
  if (tempr > kappa[0])
    kappa[0] = tempr;
  if (kappa[0] > tkappa0)
    chi = 20.0*(kappa[0] - tkappa0)*(1.0-5.0*(kappa[0]-tkappa0));
  else
    chi = 0.0;
  
  return chi;
}

/**
   function computes derivative of thermal damage function with respect to temperature
   kappa_td is notation in the notes
   
   @param ipp - integration point number
   @param tempr - actual temperature
   @param kappa - %vector of the parameters of thermal damage function
                  it contains the maximum of either the largest value attained by temperature
		  or the reference temperature

*/
double glasgmech::dtdfdt (long /*ipp*/,double tempr,vector &kappa)
{
  double htd;
  
  if (tempr > kappa[0])
    kappa[0] = tempr;
  
  if (kappa[0] > tkappa0)
    htd = 20.0*(1.0-5.0*(kappa[0]-tkappa0)) - 5.0*20.0*(kappa[0] - tkappa0);
  else
    htd = 0.0;
  
  return htd;
}

/**
   function computes increments of load induced thermal strains
   @param ipp - integration point number
   @param epstm - result vector of increments of load induced thermal strains
   @param sigt - stress tensor
   @param told - previous temperature
   @param tnew - actual temperature
   
*/
void glasgmech::compute_lits (long /*ipp*/,vector &/*epstm*/,matrix &/*sigt*/,double /*told*/,double /*tnew*/)
{
  /*
  long idem;
  double dt,i1,beta,nu,fc;
  vector princsig(3);
  matrix pvect(3,3);
  strastrestate ssst;
  
  idem=Mm->ip[ipp].gemid();
  if (Mm->ip[ipp].tm[idem] != elisomat)
    {
      fprintf(stderr, "\n\n Invalid type of elastic material is required");
      fprintf(stderr, "\n  in function glasgmech::compute_lits, (file %s, line %d)\n", __FILE__, __LINE__);
    }
  nu = Mm->eliso[idem].nu;
  
  fillv (0.0,epstm);
  
  //  increment of temperature
  dt=tnew-told;
  if (dt<=0.0){
    return;
  }
  else{
    princ_val (sigt,princsig,pvect,nijac,limit,Mp->zero,3,1);
    extractnegv (princsig,princsig);
    fillm (0.0,sigt);
    sigt[0][0]=princsig[0];
    sigt[1][1]=princsig[1];
    sigt[2][2]=princsig[2];
    lgmatrixtransf(sigt, pvect);
    
    i1=first_invar (sigt);

    beta = betacoeff (ipp,tnew);
    
    fc = st*k;
    
    cmulm ((1.0+nu),sigt);
    sigt[0][0]-=nu*i1;
    sigt[1][1]-=nu*i1;
    sigt[2][2]-=nu*i1;
    cmulm (beta/fc*dt,sigt);
    
    ssst=Mm->ip[ipp].ssst;
    tensor_vector (epstm,sigt,ssst,strain);
  }
  */
}

/**
   Function computes derivatives of equivalent strain with respect of elastic strains.
   @param ipp - integration point number
   @param epstm - result vector of increments of load induced thermal strains
   @param sigt - stress tensor
   @param told - previous temperature
   @param tnew - actual temperature
   
*/
void glasgmech::depseqdepsel(long /*ipp*/, vector &/*eps*/, vector &/*deeqdeel*/)
{
  /*
  double q1, q2, q3, e, nu, i1e, j2e, tmp;
  double r, af, bf;
  long i;
  long idem = Mm->ip[ipp].gemid();
  long ncomp = Mm->ip[ipp].ncompstr;
  vector pv(ncomp), tv(ncomp);
  matrix pm(ncomp,ncomp);
  matrix epst(3,3), dev(3,3);
  
  fillv(0.0, pv); 
  fillm(0.0, pm);
  fillv(0.0,deeqdeel);

  vector_tensor(eps, epst, Mm->ip[ipp].ssst, strain);
  i1e = first_invar(epst);
  deviator (epst,dev);
  j2e = second_invar (dev);
  idem = Mm->ip[ipp].gemid();
  if (Mm->ip[ipp].tm[idem] != elisomat)
  {
    fprintf(stderr, "\n\n Invalid type of elastic material is required");
    fprintf(stderr, "\n  in function glasgmech::depseqds, (file %s, line %d)\n", __FILE__, __LINE__);
  }
  e  = Mm->eliso[idem].e;
  nu = Mm->eliso[idem].nu;
  af  = (k-1.0)/(2.0*k)/(1.0-2.0*nu);
  bf  = 3.0/k/(1+nu);
  r   = af*af*i1e*i1e + bf*j2e;
  r   = sqrt(r);
  if (r > 0.0)
    r = 0.5 / r;
  else     
    r = 0.0;
  q3 = nu/(1.0-nu);
  q1 = 1.0 + q3 + q3*q3;
  q2 = 1.0 - 2.0*q3 - 2.0*q3*q3;
  switch (Mm->ip[ipp].ssst)
  {
    case planestrain:
    case axisymm:
      pv[0] = 1.0-q3;
      pv[1] = 1.0-q3;
      pv[2] = 0.0;
      pv[3] = 0.0;
      pm[0][0] = 2.0/3.0*q1;
      pm[1][1] = pm[0][0];
      pm[2][2] = 0.0;
      pm[3][3] = 2.0;
      pm[0][1] = -q2/3.0;
      pm[1][0] = pm[0][1];
      break;
    case spacestress:
      for (i = 0; i < 3; i++)
      {
        pv[i] = 1.0;
        pm[i][i] = 2.0/3.0;
        pm[i+3][i+3] = 2.0;
      }
      pm[0][1] = -q2/3.0; 
      pm[0][2] = pm[0][1]; 
      pm[1][2] = pm[0][1]; 
      pm[2][0] = pm[0][1]; 
      pm[2][1] = pm[0][1]; 
      break;
    default:
      fprintf(stderr, "\n\n Invalid stress/strain state is required in function glasgmech::depseqds\n");
      fprintf(stderr, " in integeration point %ld, (file %s, line %d)\n", ipp, __FILE__, __LINE__);           
  }
  for(i = ncomp / 2; i < ncomp; i++)
    eps[i] *= 0.5;
  tmp = af+2.0*af*af*i1e*r;
  cmulv(tmp, pv);
  mxv(pm, pv, tv);
  cmulv(bf*r, tv);
  addv(tv, pv, deeqdeel);
  for(i = ncomp / 2; i < ncomp; i++)
  {
    deeqdeel[i] *= 0.5;
    eps[i] *= 2.0;
  }
  */
}                



/**
   function computes correct stresses in the integration point and stores
   them into ip stress array.
   
   @param ipp - integration point pointer
   @param im - index of material
   @param ido - index of the viscous material in the array eqother
   
*/
void glasgmech::nlstresses (long /*ipp*/,long /*im*/,long /*ido*/)
{
  /*
  long i,ncompstr;
  
  //  number of strain/stress components in the problem
  ncompstr = Mm->ip[ipp].ncompstr;
  
  double newt,oldt,dt,chi,htd,omega,domega, omegao, kappao, depseq;
  
  vector epsn(ncompstr),epso(ncompstr),epse(ncompstr),epsft(ncompstr),epstm(ncompstr);
  vector epsc(ncompstr),epsirr(ncompstr),epstirr(ncompstr), deeqdeel(ncompstr), deps(ncompstr);
  vector sig(ncompstr),sigtrial(ncompstr);
  vector kappa(2),tkappa(1);
  vector deftdt(ncompstr), s(ncompstr); 
//  vector detmdt(ncompstr);
  matrix eps(3,3),d(ncompstr,ncompstr);
  matrix sigt(3,3);
  
  i = 0;
  //if (Mp->phase==1){
  if (i==1){
  */
    /* *******************************/
    //  right hand side computation  //
    /* *******************************/
    /*
    //  total new strains
    for (i=0;i<ncompstr;i++){
      epsn[i] = Mm->ip[ipp].strain[i];
    }
    //  previous total strains
    for (i=0;i<ncompstr;i++){
      epso[i] = Mm->ip[ipp].eqother[ncompstr+i];
    }
    */

    /* ************************/
    //  free thermal strains  //
    /* ************************/
    /*
    //  actual temperature
    newt = Mm->givenonmechq(temperature, ipp);
    //  temperature from the previous step
    oldt = Mm->ip[ipp].eqother[4*ncompstr];
    //  increment of temperature
    dt = newt-oldt;
    
    //  increments of free thermal strains
    compute_thermdilat (newt,dt,eps);
    //  conversion of tensor notation into vector one
    tensor_vector (epsft,eps,Mm->ip[ipp].ssst,strain);
    */
    /* ********************************/
    //  load induced thermal strains  //
    /* ********************************/
  /*
    vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
    compute_lits (ipp, epstm, sigt, oldt, newt);
   //compute_lits (epstm);
   */
    /* ******************/
    //  thermal damage  //
    /* ******************/
    /*
    //  history parameter of thermal damage
    tkappa[0]=Mm->ip[ipp].eqother[4*ncompstr+3];
    
    //  thermal damage parameter
    chi = thermdamfunction (ipp,newt,tkappa);
    */
    /* *********************/
    //  mechanical damage  //
    /* *********************/
    /*
    //  norm of equivalent strain
    damfuncpar(ipp,epsn,kappa);
    
    // previous equivalent strain norm
    kappao = Mm->ip[ipp].eqother[4*ncompstr+1];
    //  previous damage
    omegao = Mm->ip[ipp].eqother[4*ncompstr+2];
    // previous temperature
    kappa[1] = Mm->ip[ipp].eqother[4*ncompstr];
    //  mechanical damage parameter
    omega = damfunction(ipp,newt,chi,kappa);
    if (omega < kappao)
    {
      omega = omegao;
      domega = 0.0;
    }
    else
    {
      // depseq = kappa[0] - kappao;

      // increment of equvalent strain
      depseqdepsel(ipp, epsn, deeqdeel);
      subv (epsn, epso, deps);
      subv (deps, epsft, deps);
      subv (deps, epstm, deps);
      scprd(deps, deeqdeel, depseq);
      // mechanical damage increment
      domega = domegadt(ipp, newt, chi, kappa)*depseq;
      domega += domegadkmd(ipp, newt, chi, kappa)*(newt-oldt);
    }
    
    */
    /* *********/
    //  creep  //
    /* *********/
    
    //give_creep_eps (epsc);
    
    /* *************************************/
    //  increment of irreversible strains  //
    /* *************************************/
  /*
    for (i=0;i<ncompstr;i++){
      epsirr[i]=epsft[i]+epsc[i]+epstm[i];
    }
  */
    /* ******************************/
    //  storage of computed values  //
    /* ******************************/
  /*
    for (i=0;i<ncompstr;i++){
      Mm->ip[ipp].eqother[3*ncompstr+i]=epsirr[i];
    }
    Mm->ip[ipp].eqother[4*ncompstr+1]=kappa[0];
    Mm->ip[ipp].eqother[4*ncompstr+2]=omega;
    Mm->ip[ipp].eqother[4*ncompstr+3]=tkappa[0];
    Mm->ip[ipp].eqother[4*ncompstr+4]=chi;
    Mm->ip[ipp].eqother[4*ncompstr+5]=domega;
    
    //  stiffness matrix of the material
    Mm->elmatstiff(d,ipp);
    
    //  stress increments
    mxv (d,epsirr,sig);
    
    Mm->storestress (0,ipp,sig);
  }
  
  
  //if (Mp->phase==2){
  if (i==2){
    
    for (i=0;i<ncompstr;i++){
      //  new total strain
      epsn[i] = Mm->ip[ipp].strain[i];
      //  previous total strain
      epso[i] = Mm->ip[ipp].eqother[ido+ncompstr+i];
      //  previous total irreversible strains
      epstirr[i] = Mm->ip[ipp].eqother[ido+2*ncompstr+i];
      //  increment of irreversible strain
      epsirr[i] = Mm->ip[ipp].eqother[ido+3*ncompstr+i];
      
      //  total elastic strains
      epse[i]=epsn[i]-epsirr[i];
    }
    omega =Mm->ip[ipp].eqother[4*ncompstr+2];
    chi   =Mm->ip[ipp].eqother[4*ncompstr+4];
    domega=Mm->ip[ipp].eqother[4*ncompstr+5];
    //  actual temperature
    //newt = Mm->givenonmechq(temperature, ipp);
    newt = 20.0;
    //  temperature from the previous step
    oldt = Mm->ip[ipp].eqother[4*ncompstr];
    //  increment of temperature
    dt = newt-oldt;
    //  history parameter of thermal damage
    tkappa[0]=Mm->ip[ipp].eqother[4*ncompstr+3];
    //  derivative of thermal damage function with respect to temperature
    htd = dtdfdt (ipp,newt,tkappa);
    
    
    //  stiffness matrix of material
    Mm->matstiff(d,ipp);
    //  effective stress
    mxv (d,epse,sigtrial);

    //  total strain increment
    subv (epsn,epso,epso);
    //  elastic strain increment
    subv (epso,epsirr,epse);
    //  secant stiffness matrix
    cmulm ((1.0-omega)*(1.0-chi),d);
    //  stress increment
    mxv (d,epse,sig);

    cmulv (htd*dt*(1.0-omega)+domega*(1.0-chi),sigtrial);
    
    subv (sig,sigtrial,sig);
    

    //  new data storage
    Mm->ip[ipp].eqother[4*ncompstr+0]=newt;
    for (i=0;i<ncompstr;i++){
      //  new stress
      Mm->ip[ipp].eqother[ido+i] += sig[i];
      //  total strains
      Mm->ip[ipp].eqother[ido+ncompstr+i] = epsn[i];
      //  total irreversible strains
      Mm->ip[ipp].eqother[ido+2*ncompstr+i] += epsirr[i];
      //  stress increment
      Mm->ip[ipp].stress[i]=sig[i];
    }

  }
  */
}
