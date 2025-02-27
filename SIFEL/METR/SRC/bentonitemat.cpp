#include "bentonitemat.h"
#include "globalc.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"

bentonitemat::bentonitemat (void)
{
  phideg = lam_star = kap_star = n_star = 0.0;
  nu = n = l = nt = lt = m = alpha_s = kap_m = sm_star = 0.0;
  em_star = s_airentry0 = em0 = tr = at = bt = lambdap0 = aer = p_t = 0.0;

  // retrieve number of parameters, number of state variables and number of generalized vector components
  
  nparams=huet.nparms;         // number of parameters
  nstatev=huet.nstatev;        // number of state variables (9 for actual implementation, see constructor in generalmod.h)
  ngstrain=huet.ngstrain;      // number of generalised strain tensor components (8 for actual implementation, see constructor in generalmod.h)
  ngstress=huet.ngstress;      // number of generalised stress tensor components (7 for actual implementation, see constructor in generalmod.h)
  nrkfparams=huet.nrkf_parms;  // number of RKF method parameters
  nrkfstatev=huet.nrkf_statev; // number of RKF method parameters

  params    = new double[nparams];
  rkfparams = new double[nrkfparams];

  pr_suc_fl = no;
  pr_tempr_fl = no;

  // maximum value of 'zero' level of state variables for RKF time integration
  // it may be rewritten by RKF tolerance if its value is less than the value specified here 
  sv_rkf_zero = 1.0e-3;
  
  ncompstr=4;
  ncompeqother = 2*ncompstr+3+nstatev+12+ngstress*ngstrain;
  ncompother =   ncompeqother;
}

bentonitemat::~bentonitemat (void)
{
  delete [] params;
  delete [] rkfparams;
}

/**
  The function returns the number of eqother array componets of the material model.

  @param ipp - integration point number in the mechmat ip array.

  @retval The function returns value of 

  Created by Tomas Koudelka, 8.4.2019
*/
long bentonitemat::givencompeqother(long ipp)
{
  long ncompstr = Cm->ip[ipp].ncompstr;

  return 2*ncompstr+3+nstatev+12+ngstress*ngstrain;
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 5.4.2019
*/
void bentonitemat::read (XFILE *in)
{
  xfscanf(in, "%k %le %k %le %k %le %k %le",
              "phi", &phideg, "lam_star", &lam_star, 
              "kap_star", &kap_star, "N_star", &n_star);
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le %k %le",
              "nu", &nu, "n", &n, "l", &l, "nT", &nt, "lT", &lt, "m", &m);
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le %k %le %k %le",
              "alpha_s", &alpha_s, "kap_m", &kap_m, "sm_star", &sm_star, 
              "em_star", &em_star, "csh", &csh, "s_airentry0", &s_airentry0, "em0", &em0);
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le %k %le",
              "tr", &tr, "at", &at, "bt", &bt,
              "aer", &aer, "lambdap0", &lambdap0, "p_t", &p_t);

  //  xfscanf(in, "%k %le %k %le","einit", &einit, "ascaninit", &ascaninit);

  xfscanf(in, "%k%m", "presc_suction", &answertype_kwdset, &pr_suc_fl);
  if (pr_suc_fl == yes)
    suc_fn.read(in);

  xfscanf(in, "%k%m", "presc_tempr", &answertype_kwdset, &pr_tempr_fl);
  if (pr_tempr_fl == yes)
    tempr_fn.read(in);

  sra.read(in);

  if (sv_rkf_zero > sra.err)
    sv_rkf_zero = sra.err;
 
  if (phideg <= 0.0)
  {
    print_err("phi must be positive (phi = %le)", __FILE__, __LINE__, __func__, phideg);
    abort();
  }
  if (p_t < 0.0)
  {
    print_err("p_t must not be negative (p_t = %le)", __FILE__, __LINE__, __func__, p_t);
    abort();
  }
  if (lam_star <= 0.0)
  {
    print_err("lam_star must be positive (lam_star = %le)", __FILE__, __LINE__, __func__, lam_star);
    abort();
  }
  if (kap_star <= 0.0)
  {
    print_err("kap_star must be positive (kap_star = %le)", __FILE__, __LINE__, __func__, kap_star);
    abort();
  }
  if (n_star <= 0.0)
  {
    print_err("n_star must be positive (n_star = %le)", __FILE__, __LINE__, __func__, n_star);
    abort();
  }
  // temporary check for Mock-up test
  if (aer != 1.0)
  {
    print_err("temporary check for Mock-up experiment - aer is not set to 1.0 (aer = %le)", __FILE__, __LINE__, __func__, aer);
    abort();
  }
  // temporary check for Mock-up test
  if (aer != 1.0)
  {
    print_err("temporary check for Mock-up experiment - aer is not set to 1.0 (aer = %le)", __FILE__, __LINE__, __func__, aer);
    abort();
  }

  /*
  if (rr <= 0.0)
  {
    print_err("rr must be positive (rr = %le)", __FILE__, __LINE__, __func__, rr);
    abort();
  }
  if (s_e0 < 0.0)
  {
    print_err("s_e0 must not be negative (s_e0 = %le)", __FILE__, __LINE__, __func__, s_e0);
    abort();
  }
  */
  params[0] = phideg;
  params[1] = lam_star;
  params[2] = kap_star;
  params[3] = n_star;
  params[4] = nu;
  params[5] = n;
  params[6] = l;
  params[7] = nt;
  params[8] = lt;
  params[9] = m;
  params[10] = alpha_s;
  params[11] = kap_m;
  params[12] = sm_star;
  params[13] = em_star;
  params[14] = csh;
  params[15] = s_airentry0;
  params[16] = em0;
  params[17] = tr;
  params[18] = at;
  params[19] = bt;
  params[20] = aer;
  params[21] = lambdap0;
  params[22] = p_t;  

  //params[21] = einit;
  //params[22] = ascaninit;

  rkfparams[0] = sra.give_err_sig();
  rkfparams[1] = sra.give_hmin();
  rkfparams[2] = sra.give_ni();
  rkfparams[3] = sv_rkf_zero;
  rkfparams[4] = (int)sra.give_rkt();

}

/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened output text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 5.4.2019
*/
void bentonitemat::print(FILE *out)
{
  fprintf(out, "%le %le %le %le\n", 
               phideg, lam_star, kap_star, n_star);
  fprintf(out, "%le %le %le %le %le %le\n",
               nu, n, l, nt, lt, m);
  fprintf(out, "%le %le %le %le %le %le %le\n",
              alpha_s, kap_m, sm_star, em_star, csh, s_airentry0, em0);
  fprintf(out, "%le %le %le %le %le %le\n",tr, at, bt, aer, lambdap0, p_t);

  fprintf(out, "%d ", int(pr_suc_fl));
  if (pr_suc_fl == yes)
    suc_fn.print(out);

  fprintf(out, "%d ", int(pr_tempr_fl));
  if (pr_tempr_fl == yes)
    tempr_fn.print(out);

  fprintf(out, "\n");
  sra.print(out);
}


//double bentonitemat::give_saturation_degree (intpointsc *ip,long ipp)
double bentonitemat::give_saturation_degree (intpointsc &ip)
{
  return ip.eqother[2*ncompstr+3+2];
}

double bentonitemat::give_dSr_ds (intpointsc &ip)
{
  return ip.eqother[2*ncompstr+3+nstatev+8];
}

double bentonitemat::give_dSr_depsv (intpointsc &ip)
{
  return ip.eqother[2*ncompstr+3+nstatev+11];
}

double bentonitemat::give_porosity (intpointsc &ip)
{
  double e;
  e = ip.eqother[2*ncompstr+3];
  return e/(1.0+e);
}

double bentonitemat::give_pressure (intpointsc &ip)
{
  return ip.av[0];
}



/**
  The function initializes material model. Retrieves initial values of the following
  quantities

  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything but updates initial value of Sr and stress vector 
          in eqother array.

  Created by Tomas Koudelka, 5.4.2019
*/
void bentonitemat::initval (long ipp)
{
  // data handled by MEFEL
  long ncomp  =  Cm->ip[ipp].ncompstr; // number of reduced stress/strain components
  strastrestate ssst =  Cm->ip[ipp].ssst; // stress/strain state indicator
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain));
  vector deps(ASTCKVEC(ngstrain)), dstatev(ASTCKVEC(nstatev));  
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));
  double *statev = Cm->ip[ipp].eqother+2*ncomp+3;
  int herr;

  if (statev[0] == 0.0) // initial void ratio was not set => model has not been initialized yet
  {
    statev[0] = Cm->ic[ipp][0];  // initialize actual value of void ratio
    statev[7] = Cm->ic[ipp][1];  // initialize actual value of ascan

    // set initial value of effective stresses from eigenstresses
    // eigenstress is assumed to be net stress, i.e. sig_{net} = sig_{tot} - u_a.I
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
    {
      for (long i=0; i<ncomp; i++)
      {
        //Cm->ip[ipp].eqother[ncomp+i] = Cm->eigstresses[ipp][i];
        //Cm->ip[ipp].stress[i] = Cm->eigstresses[ipp][i];
      }
    }  
  }
  
  statev[1] = Cm->ip[ipp].av[0]; // initial suction pressure s (positive)

  statev[2] = 0.0;    // initial degree of saturation S_r

  statev[3] = Cm->ip[ipp].av[1];  // initial temperature T

  
  give_full_vector(sig.a, Cm->ip[ipp].stress, ssst);
  give_full_vector(eps.a, Cm->ip[ipp].strain, ssst);
  
  //  Proc zde neni konverze do ABAQUSu?
  
  // set remaining components of generalized strain vector
  eps(6) = statev[1];  // set suction to the corresponding component of generalized strain vector
  eps(7) = statev[3];  // set temperature to the corresponding component of generalized strain vector
  // initialize model
  herr = huet.soil_model(eps.a, sig.a, statev, deps.a, 1.0, NULL, params, NULL, rkfparams, 0, Cp->inner_step);
  if (herr)
  {
    print_err("hypoplasticity model cannot be initialized on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }

  statev[2] = sig(6); // store computed initial Sr

  // prepare working copy of array of internal variables
  // eqother values have been already set because statev array is a reference to the eqother array
  copyv(statev, Cm->ip[ipp].other+2*ncomp+3, nstatev);

  // set initial values for dSr/ds and dSr/dT according to generalized stiffness matrix
  herr = huet.soil_model(eps.a, sig.a, statev, NULL, 1.0, dd_gen.a, params, NULL, rkfparams, 3, Cp->inner_step);
  if (herr)
  {
    print_err("generalized stiffness matrix cannot be initialized on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  Cm->ip[ipp].eqother[2*ncomp+3+nstatev+8] = dd_gen(ngstress-1, ngstrain-2); // store dSr/ds matrix component
  Cm->ip[ipp].eqother[2*ncomp+3+nstatev+9] = dd_gen(ngstress-1, ngstrain-1); // store dSr/dT matrix component
  Cm->ip[ipp].eqother[2*ncomp+3+nstatev+11] = dd_gen(ngstress-1, 0);//+dd_gen(ngstress-1, 1)+dd_gen(ngstress-1, 2); // compute and store dSr/depsv from matrix components
  
  // actualize 'consistent' material stiffness matrix
  //copym(dd_gen, Mm->ip[ipp].eqother+2*ncomp+3+nstatev+12);

  // actualize other array by the computed values
  memcpy(Cm->ip[ipp].other, Cm->ip[ipp].eqother, sizeof(*Cm->ip[ipp].other)*Cm->ip[ipp].ncompother);
}





/**
  Function computes stiffness %matrix of the material
  in the required integration point.

  @param d - stiffness %matrix (output)
  @param ipp - number of integration point
   
  @return The function returns stiffness %matrix in the parametr d.

  31.3.2019, JK
*/
void bentonitemat::matstiff (matrix &d,long ipp)
{
  
  long i,j;
  //  pocet slozek strainu
  long ncomp = Cm->ip[ipp].ncompstr;
  //strastrestate ssst = Cm->ip[ipp].ssst;
  double aux;
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain)), dd(ASTCKMAT(6,6)), ddd(ASTCKMAT(ngstress,ngstrain));
  vector statev(ASTCKVEC(nstatev));  
 vector rkf_statev(ASTCKVEC(nrkfstatev));
 vector sig(ASTCKVEC(ngstress));
 vector eps(ASTCKVEC(ngstrain));
 vector deps(ASTCKVEC(ngstrain));

  //  strain
  //  stress
  //  prirustek strainu
  //  stavove promenne
  
  eps[0]=Cm->ip[ipp].strain[0];
  eps[1]=Cm->ip[ipp].strain[1];
  eps[2]=Cm->ip[ipp].strain[2];
  eps[3]=Cm->ip[ipp].strain[3];
  eps[4]=0.0;
  eps[5]=0.0;
  eps[6] = Cm->ip[ipp].eqother[2*ncomp+3+1];   // suction from the previous time step
  eps[7] = Cm->ip[ipp].eqother[2*ncomp+3+3];   // temperature from the previous time step

  sig[0] = Cm->ip[ipp].stress[0];
  sig[1] = Cm->ip[ipp].stress[1];
  sig[2] = Cm->ip[ipp].stress[2];
  sig[3] = Cm->ip[ipp].stress[3];
  sig[4] = 0.0;
  sig[5] = 0.0;
  sig[6] = Cm->ip[ipp].eqother[2*ncomp+3+2];   // degree of saturation Sr from the previous time step
  
  /*
  deps[0]=Cm->ip[ipp].dstrain[0];
  deps[1]=Cm->ip[ipp].dstrain[1];
  deps[2]=Cm->ip[ipp].dstrain[2];
  deps[3]=Cm->ip[ipp].dstrain[3];
  deps[4]=0.0;
  deps[5]=0.0;
  //  zde musi byt take prirustky
  deps[6] = Cm->ip[ipp].eqother[2*ncomp+3+1];   // suction from the previous time step
  deps[7] = Cm->ip[ipp].eqother[2*ncomp+3+3];   // temperature from the previous time step
  */
  
  j=0;
  for (i=2*ncomp+3;i<2*ncomp+3+nstatev;i++){
    statev[j]=Cm->ip[ipp].eqother[i];
    j++;
  }

  rkf_statev[0] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+5];
  rkf_statev[1] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+6];
  rkf_statev[2] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+7];
  rkf_statev[3] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev];

  //  Cp->jstep cislo kroku v Newton-Raphson
  //  treba zjistit dtime
  //long rkf_stat = huet.soil_model(eps.a, sig.a, statev.a, deps.a, dtime, dd_gen.a, params, rkf_statev.a, rkfparams, 1, Cp->nr_step);

  //  matice dd_gen ma slozky podle ABAQUSu
  // strain increment vector must be reoredred
  // SIFEL ordering is deps11, deps22, deps33, dgamma23, dgamma13, dgamma12
  // ABAQUS ordering is deps11, deps22, deps33, dgamma12, dgamma13, dgamma23
  
  //  dSr/dv = d 7,1
  //  dSr/dp = d 7,7
  //  Sr je v statev
  
  // retrieve actual generalized material stiffness matrix from the other array
  //makerefm(dd_gen, Mm->ip[ipp].other+2*ncomp+3+nstatev+12, ngstress, ngstrain);
  makerefm(dd_gen, ddd.a, ngstress, ngstrain);

  // reorder stiffness matrix in the full-length format
  //
  // swap the fourth and sixth columns except of the fourth and sixth elements
  for (i=0; i<ngstress; i++)
  {
    if ((i == 3) || (i==5))
      continue;
    aux = dd_gen(i,3);
    dd_gen(i,3) = dd_gen(i,5);
    dd_gen(i,5) = aux;
  }
  // swap the fourth and sixth rows except of the fourth and sixth elements
  for (i=0; i<ngstrain; i++)
  {
    if ((i == 3) || (i==5))
      continue;
    aux = dd_gen(3,i);
    dd_gen(3,i) = dd_gen(5,i);
    dd_gen(5,i) = aux;
  }
  // swap elements on intersections of the fourth and sixth columns and rows
  aux = dd_gen(3,3);
  dd_gen(3,3) = dd_gen(5,5);
  dd_gen(5,5) = aux;
  aux = dd_gen(3,5);
  dd_gen(3,5) = dd_gen(5,3);
  dd_gen(5,3) = aux;

  // extract the stiffness matrix for the mechanics
  //extractm(dd, dd_gen, 0, 6); 
  // convert stiffness matrix to the reduced format suitable for element stiffness matrix
  //tensor4_ematrix (d, dd, ssst);


}

/**
   function assembles conductivity/permeability %matrix
   
   @param d - conductivity/permeability %matrix
   @param ipp - number of integration point
   
   2.4.2019, JK
 */
void bentonitemat::matcond (matrix &d)
{
  d[0][0]=k;     d[0][1]=0.0;
  d[1][0]=0.0;   d[1][1]=k;
}


/**
   function returns c_pp coefficient
   
   @param ipp - number of integration point
   
   4.4.2019, JK
 */
double bentonitemat::c_pp_coeff (intpointsc &ip)
{
  double cpp,sr,p,n,dsdp,dsdv;
  
  //  saturation degree
  sr = give_saturation_degree (ip);
  //  porosity
  n = give_porosity (ip);
  //  dSr/dp
  dsdp = give_dSr_ds (ip);
  //  dSr/depsv
  dsdv = give_dSr_depsv (ip);
  //  pressure
  p = give_pressure (ip);
  
  cpp = (alpha-n)/ks*sr*sr + n*sr/kw + (alpha-n)/ks*p*dsdp + n*dsdv;

  return cpp;
}

/**
   function returns c_up coefficient
   
   @param ipp - number of integration point
   
   4.4.2019, JK
 */
double bentonitemat::c_up_coeff (intpointsc &ip)
{
  double cup,sr;
  
  //  saturation degree
  sr = give_saturation_degree (ip);
  
  cup = alpha*sr;
  
  return cup;
}

/**
   function returns c_pu coefficient
   
   @param ipp - number of integration point
   
   4.4.2019, JK
 */
double bentonitemat::c_pu_coeff (intpointsc &ip)
{
  double cpu,sr,p,n,dsdv;
  
  //  saturation degree
  sr = give_saturation_degree (ip);
  //  porosity
  n = give_porosity (ip);
  //  dSr/depsv
  dsdv = give_dSr_depsv (ip);
  //  pressure
  p = give_pressure (ip);
  
  cpu = sr*alpha + (alpha-n)/ks*sr*p*dsdv + n*dsdv;
  
  return cpu;
}

void bentonitemat::nlstresses (long ipp, long /*im*/)
{
  long i, j, rkf_stat;
  //  long herr;
  
  // data handled by MEFEL
  long ncomp  =  Cm->ip[ipp].ncompstr;
  long neval = long(Cm->ip[ipp].other[2*ncomp+3+nstatev+5]);
  strastrestate ssst =  Cm->ip[ipp].ssst;
  double aux, dtsub, sro;
  double pnewdt = 1.0; // reduction coefficient for new time step
  double dtime = Cp->timecon.actualforwtimeincr();
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain)), deps(ASTCKVEC(ngstrain));
  vector statev(ASTCKVEC(nstatev));
  vector epseq(ASTCKVEC(ngstrain));
  vector sigeq(ASTCKVEC(ngstress));
  vector stateveq(ASTCKVEC(nstatev));
  //  rktype rkt = sra.give_rkt(); // type of Runge-Kutta method
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));

  // counter of number of steps performed in RKF method within one time step 
  long nstep = long(Cm->ip[ipp].other[2*ncomp+3+nstatev+6]);
  // the minimum step length attained in RKF method within one time step 
  double mindt = Cm->ip[ipp].other[2*ncomp+3+nstatev+7];
  
  // stress array have to contain values from the last attained equilibrium state
  for (i=0; i<ncomp; i++)
    Cm->ip[ipp].stress[i] = Cm->ip[ipp].eqother[ncomp+i];
  
  // calculate strain increments and store them in eps temporarily
  for (i=0; i<ncomp; i++)
    eps(i) = Cm->ip[ipp].strain[i]-Cm->ip[ipp].eqother[i];

  // convert stress, strain and strain increment vectors to full 6 componet vectors
  //give_full_vector(sig.a, Cm->ip[ipp].stress, ssst);
  sig[0]=Cm->ip[ipp].stress[0];
  sig[1]=Cm->ip[ipp].stress[1];
  sig[2]=Cm->ip[ipp].stress[2];
  sig[3]=Cm->ip[ipp].stress[3];
  sig[4]=0.0;
  sig[5]=0.0;
  //give_full_vector(deps.a, eps.a, ssst);
  deps[0]=eps[0];
  deps[1]=eps[1];
  deps[2]=eps[2];
  deps[3]=eps[3];
  deps[4]=0.0;
  deps[5]=0.0;
  nullv(eps);
  //give_full_vector(eps.a, Cm->ip[ipp].eqother, ssst);
  eps[0]=Cm->ip[ipp].eqother[0];
  eps[1]=Cm->ip[ipp].eqother[1];
  eps[2]=Cm->ip[ipp].eqother[2];
  eps[3]=Cm->ip[ipp].eqother[3];
  eps[4]=0.0;
  eps[5]=0.0;

  // assemble generalized strain
  eps(6) = Cm->ip[ipp].eqother[2*ncomp+3+1];   // suction from the previous time step
  eps(7) = Cm->ip[ipp].eqother[2*ncomp+3+3];   // temperature from the previous time step

  // assemble generalized stress vector
  sro = sig(6) = Cm->ip[ipp].eqother[2*ncomp+3+2];   // degree of saturation Sr from the previous time step

  // assemble generalized strain increment vector
  //if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
  //deps(6) = suc_fn.getval(Mp->time) - Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction increment ds
  //else
  //deps(6) = Mm->givenonmechq(suction, ipp) - Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction pressure increment ds
  //if (pr_tempr_fl == yes) //prescribed temperature T
  //deps(7) = tempr_fn.getval(Mp->time) - Mm->ip[ipp].eqother[ido+2*ncomp+3+3]; // temperature increment dT
  //else
  //deps(7) = Mm->givenonmechq(temperature, ipp) - Mm->ip[ipp].eqother[ido+2*ncomp+3+3];   // temperature increment dT

  // prepare working copy of internal variables from the last equilibrium state
  for (i=2*ncomp+3, j=0; i<2*ncomp+3+nstatev; i++, j++)
    statev(j) = Cm->ip[ipp].eqother[i];

  dtsub = Cm->ip[ipp].eqother[2*ncomp+3+nstatev];

  // create copy of strain, stress and state variable vectors from the last equilibrium state
  copyv(eps, epseq);
  copyv(sig, sigeq);
  copyv(statev, stateveq);
  
  // compute new values of stress vector and state variable increments  
  vector rkf_statev(ASTCKVEC(nrkfstatev));

  rkf_statev[0] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+5];
  rkf_statev[1] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+6];
  rkf_statev[2] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev+7];
  rkf_statev[3] = Cm->ip[ipp].eqother[2*ncomp+3+nstatev];
  rkf_stat = huet.soil_model(eps.a, sig.a, statev.a, deps.a, dtime, dd_gen.a, params, rkf_statev.a, rkfparams, 1, Cp->inner_step);
  // actualize RKF state variables
  neval = long(rkf_statev[0]);
  nstep = long(rkf_statev[1]);
  mindt = rkf_statev[2];
  dtsub = rkf_statev[3];
  // store 'consistent' material stiffness matrix
  copym(dd_gen, Cm->ip[ipp].other+2*ncomp+3+nstatev+12);
  
  if (rkf_stat){
    fprintf(Out, "\nRKF failed - rkf_stat=%ld, element=%ld, ip=%ld, loc_ip=%ld, incr_step=%ld, inner_step=%ld, dt=%le", 
            rkf_stat, Cm->elip[ipp], ipp, ipp-Ct->elements[Cm->elip[ipp]].ipp+1, Cp->incr_step, Cp->inner_step, Cp->timecon.actualforwtimeincr());
    fprintf(Out, "\nsig\n");
    printv(Out, sig);
    fprintf(Out, "\neps\n");
    printv(Out, eps);
    fprintf(Out, "\ndeps\n");
    printv(Out, deps);
    fprintf(Out, "\nstatev\n");
    printv(Out, statev);
    fprintf(Out, "\n\n");
    pnewdt = 0.5;  // RKF was not able integrate quantities with required error => decrease time step in the solver
  }

  // stress vector must be reoredred to SIFEL notation
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig(3);
  sig(3) = sig(5);
  sig(5) = aux;

  // store stress in reduced notation to the integration point
  give_red_vector(sig.a, Cm->ip[ipp].stress, ssst);

  // store actual value of strain vector
  Cm->storeother(ipp, 0, ncomp, Cm->ip[ipp].strain);
  // store actual value of stress vector
  Cm->storeother(ipp, ncomp, ncomp, Cm->ip[ipp].stress);
  // store actual value of state variables
  Cm->storeother(ipp, 2*ncomp+3, nstatev, statev.a);

  // store actual value of suction
  Cm->ip[ipp].other[2*ncomp+3+1] = eps(6);
  // store actual value of temperature
  Cm->ip[ipp].other[2*ncomp+3+3] = eps(7);
  // store actual value of RK substep coefficient
  Cm->ip[ipp].other[2*ncomp+3+nstatev]= dtsub;
  // store actual value of time step coefficient
  Cm->ip[ipp].other[2*ncomp+3+nstatev+1]= pnewdt;
  // store actual number of model evaluation - only for testing purposes
  Cm->ip[ipp].other[2*ncomp+3+nstatev+5]= neval;
  // store actual number of performed RKF steps - only for testing purposes
  Cm->ip[ipp].other[2*ncomp+3+nstatev+6]= nstep;
  // store the minimum RKF step length - only for testing purposes
  Cm->ip[ipp].other[2*ncomp+3+nstatev+7]= mindt;
  // volumetric strain rate 
  Cm->ip[ipp].other[2*ncomp+3+nstatev+10]= (deps(0)+deps(1)+deps(2))/dtime;

  // calculate additional quantities used for output of results
  double epsv;
  // calculate volumetric strain
  Cm->ip[ipp].other[2*ncomp+3+nstatev+2]= epsv = eps(0)+eps(1)+eps(2);
  // suction increment
  Cm->ip[ipp].other[2*ncomp+3+nstatev+3]= eps(6)-Cm->ip[ipp].eqother[2*ncomp+3+1];
  // depsV
  Cm->ip[ipp].other[2*ncomp+3+nstatev+4]= deps(0)+deps(1)+deps(2);

  /*
  // calculate strain deviator
  epsv *= 1.0/3.0;
  eps(0) -= epsv;
  eps(1) -= epsv;
  eps(2) -= epsv;
  deviator(eps, eps);
  // vertical strain deviator component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+3]= eps(1);
  // calculate stress deviator
  deviator(sig, sig);
  // vertical stress deviator component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+4]= sig(1);
  */
}

/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 8.4.2019
*/
void bentonitemat::updateval (long ipp,long /*im*/,long /*ido*/)
{
  long i, j;
  //  long herr;
  long n = givencompeqother(ipp);
  long ncomp = Cm->ip[ipp].ncompstr;
  strastrestate ssst =  Cm->ip[ipp].ssst;
  double aux;
  vector sigeq(ASTCKVEC(ngstress));
  vector epseq(ASTCKVEC(ngstrain));
  vector stateveq(ASTCKVEC(nstatev));
  //  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));

  // convert stress, strain and strain increment vectors to full 6 componet vectors
  give_full_vector(epseq.a, Cm->ip[ipp].strain, ssst);
  give_full_vector(sigeq.a, Cm->ip[ipp].stress, ssst);

  // stress vector must be reoredred
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sigeq(3);
  sigeq(3) = sigeq(5);
  sigeq(5) = aux;

  // strain vector must be reoredred
  // SIFEL ordering is eps11, eps22, eps33, gamma23, gamma13, gamma12
  // ABAQUS ordering is eps11, eps22, eps33, gamma12, gamma13, gamma23
  aux = epseq(3);
  epseq(3) = epseq(5);
  epseq(5) = aux;

  // assemble generalized strain increment vector
  //if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    //epseq(6) = suc_fn.getval(Cp->time);  // suction s defined by user function
  //else
  //   epseq(6) = Mm->givenonmechq(suction, ipp);   // suction s from TRFEL

  //if (pr_tempr_fl == yes) //prescribed temperature T
  //epseq(7) = tempr_fn.getval(Mp->time); // temperature T defined by user function
  //else
  //epseq(7) = Mm->givenonmechq(temperature, ipp);   // temperature T from TRFEL

  // assemble generalized stress vector
  sigeq(6) = Cm->ip[ipp].other[2*ncomp+3+2];   // actual degree of saturation Sr

  // prepare working copy of internal variables from the last equilibrium state
  for (i=2*ncomp+3, j=0; i<2*ncomp+3+nstatev; i++, j++)
    stateveq(j) = Cm->ip[ipp].other[i];

  // copy internal variables from working array to the equilibrium array
  for (i=0;i<n;i++)
    Cm->ip[ipp].eqother[i]=Cm->ip[ipp].other[i];

  Cm->ip[ipp].other[2*ncomp+3+nstatev+5] = 0;   // zero neval at the end of actual time step
  Cm->ip[ipp].other[2*ncomp+3+nstatev+6] = 0;   // zero nstep at the end of actual time step
  Cm->ip[ipp].other[2*ncomp+3+nstatev+7] = 1.0; // set mindt to 1 at the end of actual time step
}


