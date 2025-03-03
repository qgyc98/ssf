#include <stdlib.h>
#include <float.h>
#include "generalmod.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "intpoints.h"
#include "element.h"
#include "tensor.h"
#include "vecttens.h"
#include "mathem.h"
#include "vector.h"
#include "vecttens.h"
#include "hypoplunsatexptherm2.h"


#define nijac 200
#define limit 1.0e-15  // zero limit for principal values of tensors (Jacobi method)



hypoplast_unsat_exp_therm::hypoplast_unsat_exp_therm()
{
  phideg = lam_star = kap_star = n_star = 0.0;
  nu = n = l = nt = lt = m = alpha_s = kap_m = sm_star = 0.0;
  em_star = s_airentry0 = em0 = tr = at = bt = lambdap0 = aer = p_t = 0.0;

  // retrieve number of parameters, number of state variables and number of generalized vector components
  
  nparams=huet.nparms;         // number of parameters
  nstatev=huet.nstatev;        // number of state variables (10 for actual implementation, see constructor in generalmod.h)
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
}



hypoplast_unsat_exp_therm::~hypoplast_unsat_exp_therm()
{
  delete [] params;
  delete [] rkfparams;
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 30.1.2014
*/
void hypoplast_unsat_exp_therm::read(XFILE *in)
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

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void hypoplast_unsat_exp_therm::print(FILE *out)
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



/**
  The function initializes material model. Retrieves initial values of the following
  quantities

  @param[in] ipp - integration point number
  @param[in] ido - index of internal variables for given material in the ipp other array
  @param[in] rinit - flag for initialization after restorage from hdbackup
  
  @return The function does not return anything but updates initial value of Sr and stress vector 
          in eqother array.

  Created by Tomas Koudelka, 1.10.2015
*/
void hypoplast_unsat_exp_therm::initval (long ipp,long ido, bool rinit)
{
  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr; // number of reduced stress/strain components
  strastrestate ssst =  Mm->ip[ipp].ssst; // stress/strain state indicator
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain));
  vector deps(ASTCKVEC(ngstrain)), dstatev(ASTCKVEC(nstatev));  
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));
  double *statev = Mm->ip[ipp].eqother+ido+2*ncomp+3;
  int herr;

  if (statev[0] == 0.0) // initial void ratio was not set => model has not been initialized yet
  {
    statev[0] = Mm->ic[ipp][0];  // initialize actual value of void ratio
    statev[7] = Mm->ic[ipp][1];  // initialize actual value of ascan

    // set initial value of effective stresses from eigenstresses
    // eigenstress is assumed to be net stress, i.e. sig_{net} = sig_{tot} - u_a.I
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
    {
      for (long i=0; i<ncomp; i++)
      {
        Mm->ip[ipp].eqother[ido+ncomp+i] = Mm->eigstresses[ipp][i];
        Mm->ip[ipp].stress[i] = Mm->eigstresses[ipp][i];
      }
    }  
  }
  

  if (pr_suc_fl == yes)
    statev[1] = suc_fn.getval(Mp->time); //initial prescribed suction pressure s (positive)
  else
    statev[1] = Mm->givenonmechq(suction, ipp); // initial suction pressure s (positive)

  statev[2] = 0.0;    // initial degree of saturation S_r

  if (pr_tempr_fl == yes)
    statev[3] = tempr_fn.getval(Mp->time);  //initial prescribed temperature
  else
    statev[3] = Mm->givenonmechq(temperature, ipp);  // initial temperature T

  
  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(eps.a, Mm->ip[ipp].strain, ssst);

  // set remaining components of generalized strain vector
  eps(6) = statev[1];  // set suction to the corresponding component of generalized strain vector
  eps(7) = statev[3];  // set temperature to the corresponding component of generalized strain vector
  // initialize model
  int flag = 0;
  if (rinit == true)  
    flag = -1;
  herr = huet.soil_model(eps.a, sig.a, statev, deps.a, 1.0, NULL, params, NULL, rkfparams, flag, Mp->jstep, ipp);
  if (herr)
  {
    print_err("hypoplasticity model cannot be initialized on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }

  statev[2] = sig(6); // store computed initial Sr

  // prepare working copy of array of internal variables
  // eqother values have been already set because statev array is a reference to the eqother array
  copyv(statev, Mm->ip[ipp].other+ido+2*ncomp+3, nstatev);

  // set initial values for dSr/ds and dSr/dT according to generalized stiffness matrix
  herr = huet.soil_model(eps.a, sig.a, statev, NULL, 1.0, dd_gen.a, params, NULL, rkfparams, 3, Mp->jstep, ipp);
  if (herr)
  {
    print_err("generalized stiffness matrix cannot be initialized on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  Mm->ip[ipp].eqother[2*ncomp+3+nstatev+8] = dd_gen(ngstress-1, ngstrain-2); // store dSr/ds matrix component
  Mm->ip[ipp].eqother[2*ncomp+3+nstatev+9] = dd_gen(ngstress-1, ngstrain-1); // store dSr/dT matrix component
  Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev+11] = dd_gen(ngstress-1, 0);//+dd_gen(ngstress-1, 1)+dd_gen(ngstress-1, 2); // compute and store dSr/depsv from matrix components
  
  // actualize 'consistent' material stiffness matrix
  copym(dd_gen, Mm->ip[ipp].eqother+ido+2*ncomp+3+nstatev+12);

  // actualize other array by the computed values
  memcpy(Mm->ip[ipp].other, Mm->ip[ipp].eqother, sizeof(*Mm->ip[ipp].other)*Mm->ip[ipp].ncompother);
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 1.10.2015
*/
void hypoplast_unsat_exp_therm::matstiff (matrix &d,long ipp,long /*ido*/)
{
  long i;

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double aux;
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain)), dd(ASTCKMAT(6,6));

  /*
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain)), deps(ASTCKVEC(ngstrain));
  vector statev(ASTCKVEC(nstatev));  
  int herr;
  long j;

  // convert stress and strain vectors to full 6 component vectors
  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(eps.a, Mm->ip[ipp].strain, ssst);

  // stress vector must be reoredred
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig(3);
  sig(3) = sig(5);
  sig(5) = aux;

  // strain vector must be reoredred
  // SIFEL ordering is eps11, eps22, eps33, gamma23, gamma13, gamma12
  // ABAQUS ordering is eps11, eps22, eps33, gamma12, gamma13, gamma23
  aux = eps(3);
  eps(3) = eps(5);
  eps(5) = aux;

  // assemble generalized strain, stress and strain increment vectors
  eps(6) = Mm->ip[ipp].other[ido+2*ncomp+3+1];   // suction from the actual time step
  eps(7) = Mm->ip[ipp].other[ido+2*ncomp+3+3];   // temperature from the actual time step
  sig(6) = Mm->ip[ipp].other[ido+2*ncomp+3+2];   // degree of saturation Sr from the actual time step

  // prepare working copy of internal variables
  for (i=ido+2*ncomp+3, j=0; i<ido+2*ncomp+3+nstatev; i++, j++)
    statev(j) = Mm->ip[ipp].other[i];

  herr = huet.soil_model(eps.a, sig.a, statev.a, NULL, NULL, NULL, dd_gen.a, NULL, 1, Mp->jstep);
  if (herr)
  {
    print_err("math error detected in hypoplastic model on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  */

  // retrieve actual generalized material stiffness matrix from the other array
  makerefm(dd_gen, Mm->ip[ipp].other+2*ncomp+3+nstatev+12, ngstress, ngstrain);

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
  extractm(dd, dd_gen, 0, 6); 
  // convert stiffness matrix to the reduced format suitable for element stiffness matrix
  tensor4_ematrix (d, dd, ssst);
}



/**
  The function computes stress increments due to pore pressure change in the integration point and stores
  them into eqother array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything but stores the resulting values in the eqother array.

  Created by Tomas Koudelka, 11.2015
*/
void hypoplast_unsat_exp_therm::nlstressesincr (long ipp, long /*im*/, long ido)
{
  long i, j;
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double aux;
  double ds;
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain));
  vector statev(ASTCKVEC(nstatev));  
  matrix dd_gen;
  //int herr;
  
  // This procedure is supposed to be called at the beginning of the actual time step
  // and therefore the ip.stress and ip.strain arrays must contain values 
  // from the last equilibrium step

  // convert stress and strain vectors to full 6 component vectors
  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(eps.a, Mm->ip[ipp].strain, ssst);

  // stress vector must be reoredred
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig(3);
  sig(3) = sig(5);
  sig(5) = aux;

  // strain vector must be reoredred
  // SIFEL ordering is eps11, eps22, eps33, gamma23, gamma13, gamma12
  // ABAQUS ordering is eps11, eps22, eps33, gamma12, gamma13, gamma23
  aux = eps(3);
  eps(3) = eps(5);
  eps(5) = aux;

  // assemble generalized strain, stress and strain increment vectors
  eps(6)  = Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction from the previous time step
  eps(7)  = Mm->ip[ipp].eqother[ido+2*ncomp+3+3];   // temperature from the previous time step
  sig(6) = Mm->ip[ipp].eqother[ido+2*ncomp+3+2];   // degree of saturation Sr from the previous time step

  // prepare working copy of internal variables
  for (i=ido+2*ncomp+3, j=0; i<ido+2*ncomp+3+nstatev; i++, j++)
    statev(j) = Mm->ip[ipp].eqother[i];
  
  // retrieve actual generalized material stiffness matrix from the other array
  makerefm(dd_gen, Mm->ip[ipp].other+2*ncomp+3+nstatev+12, ngstress, ngstrain);

  /*
  reallocm(RSTCKMAT(dd_gen, ngstress,ngstrain));
  herr = huet.soil_model(eps.a, sig.a, statev.a, NULL, NULL, NULL, dd_gen.a, NULL, 3, Mp->jstep);
  if (herr)
  {
    print_err("math error detected in hypoplastic model on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  */

  if (pr_suc_fl == yes) //prescribed suction pressure s
    ds = (suc_fn.getval(Mp->time)-Mm->ip[ipp].eqother[ido+2*ncomp+3+1]);  //suction pressure increment ds = du_a - du_w
  else // suction is obtained from METR
    ds = Mm->givenonmechq(suction, ipp) - Mm->ip[ipp].eqother[ido+2*ncomp+3+1];

  // calculate negative increment of stresses due to change in suction (-du = -ds*chi)
  // which will be used in calculation of nodal forces and added to the right hand side vector
  // K*dr = df_e - \int B^T*du dV
  for (i=0;i<3;i++)
    Mm->ip[ipp].eqother[ido+2*ncomp+i] = -ds*dd_gen(i,6);
  //Mm->ip[ipp].eqother[ido+2*ncomp+i] = 0.0;
}



/**
  The function computes stress increment due to pore pressure change  in the integration point and stores
  them into ip stress array.

  @param lcid - load case id
  @param ipp - integration point number
  @param im  - index of the material in the tm and idm arrays on integration point
  @param ido - index of internal variables in the ip's ipp eqother array
  @param fi  - first index of the required stress increment component
  @param sig - %vector containing stress increment components (output)
   
  @return The function returns %vector of stress increments in the parameter sig.

  Created by Tomas Koudelka, 11.2015
*/
void hypoplast_unsat_exp_therm::givestressincr (long /*lcid*/, long ipp, long /*im*/, long ido, long fi, vector &sig)
{
  nullv(sig);  

  long i, j, ncomp = Mm->ip[ipp].ncompstr;

  for (i=fi, j=0; i<3 && j<sig.n;  i++, j++)
    sig(j) = Mm->ip[ipp].eqother[ido+2*ncomp+i];
}




/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 1.10.2015
*/
void hypoplast_unsat_exp_therm::nlstresses (long ipp, long /*im*/, long ido)
{
  long i, j, rkf_stat;
  //  long herr;

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  long neval = long(Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+5]);
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double aux, dtsub, sro;
  double pnewdt = 1.0; // reduction coefficient for new time step
  double dtime = Mp->timecon.actualforwtimeincr();
  vector sig(ASTCKVEC(ngstress)), eps(ASTCKVEC(ngstrain)), deps(ASTCKVEC(ngstrain));
  vector statev(ASTCKVEC(nstatev));
  vector epseq(ASTCKVEC(ngstrain));
  vector sigeq(ASTCKVEC(ngstress));
  vector stateveq(ASTCKVEC(nstatev));
  //  rktype rkt = sra.give_rkt(); // type of Runge-Kutta method
  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));

  // counter of number of steps performed in RKF method within one time step 
  long nstep = long(Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+6]);
  // the minimum step length attained in RKF method within one time step 
  double mindt = Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+7];
  
  // stress array have to contain values from the last attained equilibrium state
  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].stress[i] = Mm->ip[ipp].eqother[ido+ncomp+i];

  // calculate strain increments and store them in eps temporarily
  for (i=0; i<ncomp; i++)
    eps(i) = Mm->ip[ipp].strain[i]-Mm->ip[ipp].eqother[ido+i];

  // convert stress, strain and strain increment vectors to full 6 componet vectors
  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(deps.a, eps.a, ssst);
  nullv(eps);
  give_full_vector(eps.a, Mm->ip[ipp].eqother+ido, ssst);

  // stress vector must be reoredred
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig(3);
  sig(3) = sig(5);
  sig(5) = aux;

  // strain vector must be reoredred
  // SIFEL ordering is eps11, eps22, eps33, gamma23, gamma13, gamma12
  // ABAQUS ordering is eps11, eps22, eps33, gamma12, gamma13, gamma23
  aux = eps(3);
  eps(3) = eps(5);
  eps(5) = aux;

  // strain increment vector must be reoredred
  // SIFEL ordering is deps11, deps22, deps33, dgamma23, dgamma13, dgamma12
  // ABAQUS ordering is deps11, deps22, deps33, dgamma12, dgamma13, dgamma23
  aux = deps(3);
  deps(3) = deps(5);
  deps(5) = aux;

  // assemble generalized strain
  eps(6) = Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction from the previous time step
  eps(7) = Mm->ip[ipp].eqother[ido+2*ncomp+3+3];   // temperature from the previous time step

  // assemble generalized stress vector
  sro = sig(6) = Mm->ip[ipp].eqother[ido+2*ncomp+3+2];   // degree of saturation Sr from the previous time step

  // assemble generalized strain increment vector
  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    deps(6) = suc_fn.getval(Mp->time) - Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction increment ds
  else
    deps(6) = Mm->givenonmechq(suction, ipp) - Mm->ip[ipp].eqother[ido+2*ncomp+3+1];   // suction pressure increment ds
  if (pr_tempr_fl == yes) //prescribed temperature T
    deps(7) = tempr_fn.getval(Mp->time) - Mm->ip[ipp].eqother[ido+2*ncomp+3+3]; // temperature increment dT
  else
    deps(7) = Mm->givenonmechq(temperature, ipp) - Mm->ip[ipp].eqother[ido+2*ncomp+3+3];   // temperature increment dT

  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp+3, j=0; i<ido+2*ncomp+3+nstatev; i++, j++)
    statev(j) = Mm->ip[ipp].eqother[i];

  dtsub = Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev];

  // create copy of strain, stress and state variable vectors from the last equilibrium state
  copyv(eps, epseq);
  copyv(sig, sigeq);
  copyv(statev, stateveq);

  // compute new values of stress vector and state variable increments  
  vector rkf_statev(ASTCKVEC(nrkfstatev));

  rkf_statev[0] = Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev+5];
  rkf_statev[1] = Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev+6];
  rkf_statev[2] = Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev+7];
  rkf_statev[3] = Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev];
  
  //if(Mp->time >= 90899 && ipp == 0)
  //printf("Time = %lf\n",Mp->time);

  rkf_stat = huet.soil_model(eps.a, sig.a, statev.a, deps.a, dtime, dd_gen.a, NULL, rkf_statev.a, rkfparams, 1, Mp->jstep, ipp);
  // actualize RKF state variables
  neval = long(rkf_statev[0]);
  nstep = long(rkf_statev[1]);
  mindt = rkf_statev[2];
  dtsub = rkf_statev[3];
  // store 'consistent' material stiffness matrix
  copym(dd_gen, Mm->ip[ipp].other+ido+2*ncomp+3+nstatev+12);
  /*
  switch (rkt)
  {
    case fwdeulert:
    {
      vector dsig(ASTCKVEC(ngstress));
      vector dstatev(ASTCKVEC(nstatev));
      rkf_stat = huet.soil_model(eps.a, sig.a, statev.a, deps.a, dsig.a, dstatev.a, NULL, NULL, 2, Mp->jstep);
      addv(sig, dsig, sig);
      addv(statev, dstatev, statev);
      neval += 1;
      nstep += 1;
      mindt = 1.0;
      break;
    }
    case heunt:
      rkf_stat = adfwdeuler(eps, sig, statev, deps, dtsub, neval, nstep, mindt);
      break;
    case rkf23t:
      rkf_stat = rkf23(ipp, eps, sig, statev, deps, dtsub, neval, nstep, mindt);
      break;
    case rkf23bst:
      rkf_stat = rkf23bs(eps, sig, statev, deps, dtsub, neval, nstep, mindt);
      break;
    case rkf34t:
      rkf_stat = rkf34(eps, sig, statev, deps, dtsub, neval, nstep, mindt);
      break;
    case rkf45t:
      rkf_stat = rkf45(eps, sig, statev, deps, dtsub, neval, nstep, mindt);
      break;
    default:
      print_err("unknown Runge-Kutta method is required", __FILE__, __LINE__, __func__);
  }
  */
  if (rkf_stat){
    fprintf(Out, "\nRKF failed - rkf_stat=%ld, element=%ld, ip=%ld, loc_ip=%ld, istep=%ld, jstep=%ld, dt=%le", 
            rkf_stat, Mm->elip[ipp], ipp, ipp-Mt->elements[Mm->elip[ipp]].ipp[0][0]+1, Mp->istep, Mp->jstep, Mp->timecon.actualforwtimeincr());
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
  give_red_vector(sig.a, Mm->ip[ipp].stress, ssst);

  // actual value of generalized strain vector
  // addv(eps,deps,eps);

  // store actual value of strain vector
  Mm->storeother(ipp, ido, ncomp, Mm->ip[ipp].strain);
  // store actual value of stress vector
  Mm->storeother(ipp, ido+ncomp, ncomp, Mm->ip[ipp].stress);
  // store actual value of state variables
  Mm->storeother(ipp, ido+2*ncomp+3, nstatev, statev.a);
  // store actual value of suction
  Mm->ip[ipp].other[ido+2*ncomp+3+1] = eps(6);
  // store actual value of temperature
  Mm->ip[ipp].other[ido+2*ncomp+3+3] = eps(7);
  // store actual value of RK substep coefficient
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev]= dtsub;
  // store actual value of time step coefficient
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+1]= pnewdt;
  // store actual number of model evaluation - only for testing purposes
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+5]= neval;
  // store actual number of performed RKF steps - only for testing purposes
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+6]= nstep;
  // store the minimum RKF step length - only for testing purposes
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+7]= mindt;
  // volumetric strain rate is computed via generalized trapesoidal rule
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+10]= 0.5*((deps(0)+deps(1)+deps(2))/dtime + Mm->ip[ipp].eqother[ido+2*ncomp+3+nstatev+10]);
  // volumetric strain rate
  //Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+10]= (deps(0)+deps(1)+deps(2))/dtime;
  /*
  // assemble non-mechanical blocks (all except of sig-eps block) of generalized stiffness matrix
  herr = huet.soil_model(epseq.a, sigeq.a, stateveq.a, NULL, NULL, NULL, dd_gen.a, NULL, 3, Mp->jstep);
  if (herr)
  {
    print_err("generalized stiffness matrix cannot be assembled on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  // store derivatives of degree of saturation Sr
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+8] = dd_gen(ngstress-1, ngstrain-2); // store dSr/ds matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+9] = dd_gen(ngstress-1, ngstrain-1); // store dSr/dT matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+11] = dd_gen(ngstress-1, 0)+dd_gen(ngstress-1, 1)+dd_gen(ngstress-1, 2); // compute and store dSr/depsv from matrix components
  */

  // calculate additional quantities used for output of results
  double epsv;
  // calculate volumetric strain
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+2]= epsv = eps(0)+eps(1)+eps(2);
  // suction increment
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+3]= eps(6)-Mm->ip[ipp].eqother[ido+2*ncomp+3+1];
  // depsV
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+4]= deps(0)+deps(1)+deps(2);
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
  The function reduces step length in Runge-Kutta-Fehlberg methods according to given 
  coefficient and minimum step length.

  @param rc   - reduction coefficient of step length (input)
  @param h    - actual/modified step length (input/output)
  @param hmin - minimum step length (input)

  @retval 0 - On sucessfull step reduction.
  @retval 1 - The minimum step length has been attained and the step cannot be reduced further.

  Created by Tomas Koudelka, 24.5.2016
*/
long hypoplast_unsat_exp_therm::rkf_redstep(double rc, double &h, double hmin)
{
  if (h > hmin)
  {
    h *= rc;
    if (h < hmin)
      h = hmin;
    return 0;
  }
  else
    return 1;
}



/**
  The function returns time step size required by the hypoplsaticity model. It is represented by ratio 
  of required time step size to the actual one or 1.0 for no change in time step size.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns required time step size provided by the hypoplasticity model

  Created by Tomas Koudelka, 1.10.2015
*/
double hypoplast_unsat_exp_therm::dstep_red(long ipp, long /*im*/, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  return (Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+1]);
}



/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 1.10.2015
*/
void hypoplast_unsat_exp_therm::updateval (long ipp,long im,long ido)
{
  long i, j;
  //  long herr;
  long n = Mm->givencompeqother(ipp,im);
  long ncomp = Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double aux;
  vector sigeq(ASTCKVEC(ngstress));
  vector epseq(ASTCKVEC(ngstrain));
  vector stateveq(ASTCKVEC(nstatev));
  //  matrix dd_gen(ASTCKMAT(ngstress,ngstrain));

  // convert stress, strain and strain increment vectors to full 6 componet vectors
  give_full_vector(epseq.a, Mm->ip[ipp].strain, ssst);
  give_full_vector(sigeq.a, Mm->ip[ipp].stress, ssst);

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
  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    epseq(6) = suc_fn.getval(Mp->time);  // suction s defined by user function
  else
    epseq(6) = Mm->givenonmechq(suction, ipp);   // suction s from TRFEL

  if (pr_tempr_fl == yes) //prescribed temperature T
    epseq(7) = tempr_fn.getval(Mp->time); // temperature T defined by user function
  else
    epseq(7) = Mm->givenonmechq(temperature, ipp);   // temperature T from TRFEL

  // assemble generalized stress vector
  sigeq(6) = Mm->ip[ipp].other[ido+2*ncomp+3+2];   // actual degree of saturation Sr

  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp+3, j=0; i<ido+2*ncomp+3+nstatev; i++, j++)
    stateveq(j) = Mm->ip[ipp].other[i];

  /*
  // assemble non-mechanical blocks (all except of sig-eps block) of generalized stiffness matrix
  herr = huet.soil_model(epseq.a, sigeq.a, stateveq.a, NULL, NULL, NULL, dd_gen.a, NULL, 3, Mp->jstep);
  if (herr)
  {
    print_err("generalized stiffness matrix cannot be assembled on element %ld (ip=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
    abort();
  }
  // store derivatives of degree of saturation Sr
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+8] = dd_gen(ngstress-1, ngstrain-2); // store dSr/ds matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+9] = dd_gen(ngstress-1, ngstrain-1); // store dSr/dT matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+11] = dd_gen(ngstress-1, 0)+dd_gen(ngstress-1, 1)+dd_gen(ngstress-1, 2); // compute dSr/depsv from matrix components
  */
  matrix dd_gen;
  makerefm(dd_gen, Mm->ip[ipp].other+ido+2*ncomp+3+nstatev+12, ngstress, ngstrain);

  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+8] = dd_gen(ngstress-1, ngstrain-2); // store dSr/ds matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+9] = dd_gen(ngstress-1, ngstrain-1); // store dSr/dT matrix component
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+11] = dd_gen(ngstress-1, 0);//+dd_gen(ngstress-1, 1)+dd_gen(ngstress-1, 2); // compute dSr/depsv from matrix components

  // copy internal variables from working array to the equilibrium array
  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];

  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+5] = 0;   // zero neval at the end of actual time step
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+6] = 0;   // zero nstep at the end of actual time step
  Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+7] = 1.0; // set mindt to 1 at the end of actual time step
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
  Created by Tomas Koudelka,  1.10.2015
*/
void hypoplast_unsat_exp_therm::give_reqnmq(long *anmq)
{
  //  anmq[saturation_degree-1] = 1;
  if (pr_suc_fl == no)
    anmq[suction-1] = 1;
  if (pr_tempr_fl == no)
    anmq[temperature-1] = 1;
}



/**
  The function returns the number of eqother array componets of the material model.

  @param ipp - integration point number in the mechmat ip array.

  @retval The function returns value of 

  Created by Tomas Koudelka, 4.3.2019
*/
long hypoplast_unsat_exp_therm::givencompeqother(long ipp)
{
  long ncompstr = Mm->ip[ipp].ncompstr;

  return 2*ncompstr+3+nstatev+12+ngstress*ngstrain;
}



/**
  The function extracts actual degree of saturation for the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of actual degree of saturation Sr.

  Created by Tomas Koudelka, 25.5.2016
*/
double hypoplast_unsat_exp_therm::give_sr(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;

  return Mm->ip[ipp].other[ido+2*ncompstr+3+2];
}



/**
  The function computes actual derivative of degree of saturation with respect to suction 
  for the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of derivative of actual degree of saturation with respect to suction dSr/ds.

  Created by Tomas Koudelka, 25.5.2016
*/
double hypoplast_unsat_exp_therm::give_dsr_ds(long ipp, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double dsr = Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+8];
  return dsr;
}



/**
  The function computes actual derivative of degree of saturation with respect to temperature 
  for the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of derivative of actual degree of saturation with respect to suction dSr/ds.

  Created by Tomas Koudelka, 29.6.2017
*/
double hypoplast_unsat_exp_therm::give_dsr_dtemp(long ipp, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double dsr = Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+9];
  return dsr;
}



/**
  The function computes actual derivative of degree of saturation with respect to volumetric strain 
  for the attained equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of derivative of actual degree of saturation with respect to suction dSr/ds.

  Created by Tomas Koudelka, 29.6.2018
*/
double hypoplast_unsat_exp_therm::give_dsr_depsv(long ipp, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double dsr = Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+11];
  return dsr;
}



/**
  The function extracts actual porosity for the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of actual porosity.

  Created by Tomas Koudelka, 25.5.2016
*/
double hypoplast_unsat_exp_therm::give_porosity(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;

  double e = Mm->ip[ipp].other[ido+2*ncompstr+3]; // actual void ratio
  return (e/(1.0+e));
}



/**
  The function extracts actual void ratio for the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of actual porosity e.

  Created by Tomas Koudelka, 25.5.2016
*/
double hypoplast_unsat_exp_therm::give_void_ratio(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;

  return Mm->ip[ipp].other[ido+2*ncompstr+3];
}



/**
  The function returns rate of the volumetric strain at the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by Tomas Koudelka 05.2018
*/
double hypoplast_unsat_exp_therm::give_strain_vol_rate(long ipp, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;

  return Mm->ip[ipp].other[ido+2*ncomp+3+nstatev+10];
}



/**
  The function returns approximation of the bulk modulus from the stored 
  tangent stiffness %matrix.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns approximated value of the bulk modulus.
  
  Created by Tomas Koudelka 09.2023
*/
double hypoplast_unsat_exp_therm::give_bulk_modulus(long ipp, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  matrix d;
  makerefm(d, Mm->ip[ipp].other+ido+2*ncomp+3+nstatev+12, ngstress, ngstrain);
  double k = 0.0;
  k = d(0,0) + d(1,1) + d(2,2);
  k /= 3.0;
  return k;  
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,  30.1.2014
*/
void hypoplast_unsat_exp_therm::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      phideg=val[i];
      break;
    }
    case 1:{
      lam_star=val[i];
      break;
    }
    case 2:{
      kap_star=val[i];
      break;
    }
    case 3:{
      n_star=val[i];
      break;
    }
    case 4:{
      nu=val[i];
      break;
    }
    case 5:{
      n=val[i];
      break;
    }
    case 6:{
      l=val[i];
      break;
    }
    case 7:{
      nt=val[i];
      break;
    }
    case 8:{
      lt=val[i];
      break;
    }
    case 9:{
      m=val[i];
      break;
    }
    case 10:{
      alpha_s=val[i];
      break;
    }
    case 11:{
      kap_m=val[i];
      break;
    }
    case 12:{
      sm_star=val[i];
      break;
    }
    case 13:{
      em_star=val[i];
      break;
    }
    case 14:{
      s_airentry0=val[i];
      break;
    }
    case 15:{
      em0=val[i];
      break;
    }
    case 16:{
      tr=val[i];
      break;
    }
    case 17:{
      at=val[i];
      break;
    }
    case 18:{
      bt=val[i];
      break;
    }
    case 19:{
      aer=val[i];
      break;
    }
    case 20:{
      p_t=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}
