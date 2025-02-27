#include "umatunsat2.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include "mathem.h"
#include "vector.h"
#include <math.h>

#define nijac 200
#define limit 1.0e-15  // zero limit for principal values of tensors (Jacobi method)

hypoplast::hypoplast()
{
  phideg = phi = p_t = lam_star = kap_star = n_star = 0.0;
  rr = n = l = m = s_e0 = e_0;
  nprops = 11;
  nstatv = 6;
  pr_suc_fl = no;
}



hypoplast::~hypoplast()
{
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 30.1.2014
*/
void hypoplast::read(XFILE *in)
{
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le",
              "phi", &phideg, "p_t", &p_t, "lam_star", &lam_star, 
              "kap_star", &kap_star, "N_star", &n_star);
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le %k %le",
              "rr", &rr, "n", &n, "l", &l, "m", &m, 
              "s_e0", &s_e0, "e_0", &e_0);
  xfscanf(in, "%k%m", "presc_suction", &answertype_kwdset, &pr_suc_fl);
  if (pr_suc_fl == yes)
    suc_fn.read(in);

  phi = phideg*M_PI/180.0; // convert phi to radians

  if (phi <= 0.0)
  {
    print_err("phi must be positive (phi = %le)", __FILE__, __LINE__, __func__, phi);
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
}



void hypoplast::initval (long ipp,long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double *statev = Mm->ip[ipp].eqother+ido+2*ncomp;

  if (statev[0] == 0.0) // initial void ratio was not set  model has not been initialized yet
  {
    statev[0] = Mm->ic[ipp][0];  // initialize actual value of void ratio

    // set initial value of stresses from eigenstresses
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
    {
      for (long i=0; i<ncomp; i++)
        Mm->ip[ipp].eqother[ido+ncomp+i] = Mm->eigstresses[ipp][i];
    }  
  }

  // statev[1] = Mm->givenonmechq(suction, ipp);            // initial suction pressure s (positive)
  // statev[2] = Mm->givenonmechq(saturation_degree, ipp);  // initial degree of saturation S_r

  // ???!! docasne
  if (pr_suc_fl == yes)
    statev[1] = suc_fn.getval(Mp->time);  //initial prescribed suction pressure s (positive)
  else
    statev[1] = 200.0;  // initial suction pressure s (positive)

  statev[2] = 0.0;    // initial degree of saturation S_r
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 11.2.2014
*/
void hypoplast::matstiff (matrix &d,long ipp,long ido)
{
  long i, j;

  // dummy variables
  double ddum=0.0;

  // should be updated according METR/TRFEL
  vector unsatvar(ASTCKVEC(5));

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double aux, time, dtime;
  vector sig(ASTCKVEC(6)), eps(ASTCKVEC(6)), deps(ASTCKVEC(6)), statev(ASTCKVEC(nstatv));  
  matrix dd(ASTCKMAT(6,6));
  
  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp, j=0; i<ido+2*ncomp+nstatv; i++, j++)
  {
    Mm->ip[ipp].other[i] = Mm->ip[ipp].eqother[i];
    statev(j) = Mm->ip[ipp].eqother[i];
  }
/*
  unsatvar[0]= Mm->givenonmechq(suction, ipp);                    //suction pressure s (positive)
  unsatvar[1]= unsatvar[0]-Mm->ip[ipp].eqother[ido+ncomp+nstatv]; //suction pressure increment ds (positive)
  unsatvar[2]= Mm->givenonmechq(saturation_degree, ipp);          // degree of saturation S_r
*/

  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
  {
    unsatvar[0] = Mm->ip[ipp].eqother[ido+2*ncomp+1];   // suction from the previous time step
    unsatvar[1] = suc_fn.getval(Mp->time)-unsatvar[0];  //suction pressure increment ds (positive)
  }
  else
  {
    unsatvar[0] = 200.0;  //constant suction pressure s (positive) (from the previous time step)
    unsatvar[1] = 0.0;    // suction pressure increment
  }
  unsatvar[2]= 0.0;   // degree of saturation S_r

  statev[1]=unsatvar[0]; // suction pressure s (positive)
  statev[2]=unsatvar[2]; // degree of saturation S_r

  // set stress array have to contain values from the last attained equilibrium state
  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].stress[i] = Mm->ip[ipp].eqother[ido+ncomp+i];

  // convert stress and strain vector to full 6 componet vectors
  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(eps.a, Mm->ip[ipp].strain, ssst);



  // stress vector must be reoredred
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig[3];
  sig[3] = sig[5];
  sig[5] = aux;

  // strain vector must be reoredred
  // SIFEL ordering is eps11, eps22, eps33, eps23, eps13, eps12
  // ABAQUS ordering is eps11, eps22, eps33, eps12, eps13, sig23
  aux = eps[3];
  eps[3] = eps[5];
  eps[5] = aux;

  if (Mp->timecon.tct == notct)
  {
    dtime = Mp->dlambda;
    time = Mp->lambda;
  }
  else
  {
    dtime = Mp->dtime;
    time = Mp->time;
  }

  umatunsat(ipp, sig, statev, dd, eps, deps, dtime, ddum, unsatvar);

  // reorder stiffness matrix in the full-length format
  //
  // swap the fourth and sixth columns except of the fourth and sixth elements
  for (i=0; i<5; i++)
  {
    if (i == 3)
      continue;
    aux = dd[i][3];
    dd[i][3] = dd[i][5];
    dd[i][5] = aux;
  }
  // swap the fourth and sixth rows except of the fourth and sixth elements
  for (i=0; i<5; i++)
  {
    if (i == 3)
      continue;
    aux = dd[3][i];
    dd[3][i] = dd[5][i];
    dd[5][i] = aux;
  }
  // swap elements on intersections of the fourth and sixth columns and rows
  aux = dd[3][3];
  dd[3][3] = dd[5][5];
  dd[5][5] = aux;
  aux = dd[3][5];
  dd[3][5] = dd[5][3];
  dd[5][3] = aux;

  // convert stiffness matrix to the reduced format
  tensor4_ematrix (d, dd, ssst);
}



/**
  The function computes stress increment due to pore pressure change  in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  @param sig - %vector of stress increment due to suction change (output) 

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2.2014
*/
void hypoplast::givestressincr (long ipp, long /*im*/, long ido, vector &sig)
{
  long i;

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double ds;
  matrix sigt(ASTCKMAT(3,3));

  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    ds = suc_fn.getval(Mp->time)-Mm->ip[ipp].eqother[ido+2*ncomp+1];  //suction pressure increment ds (positive)
  else
    ds = 0.0;

  // calculate stress increment tensor due to change in suction (dsig = -dpp = ds)
  for (i=0;i<3;i++)
    sigt[i][i] += ds;
  
  // convert total stress tensor back to vector form - sig
  tensor_vector(sig, sigt, Mm->ip[ipp].ssst,stress);
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2.2014
*/
void hypoplast::nlstresses (long ipp, long /*im*/, long ido)
{
  long i, j;

  // should be updated according METR/TRFEL
  vector unsatvar(ASTCKVEC(5));

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  vector statev(ASTCKVEC(nstatv));
  double aux, time, dtime;
  double pnewdt = 1.0; // reduction coefficient for new time step
  vector sig(ASTCKVEC(6)), eps(ASTCKVEC(6)), deps(ASTCKVEC(6));  
  matrix dd(ASTCKMAT(6,6));



  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp, j=0; i<ido+2*ncomp+nstatv; i++, j++)
  {
    Mm->ip[ipp].other[i] = Mm->ip[ipp].eqother[i];
    statev(j) = Mm->ip[ipp].eqother[i];
  }

/*  
  unsatvar[0]= Mm->givenonmechq(suction, ipp);                      //suction pressure s (positive)
  unsatvar[1]= unsatvar[0]-Mm->ip[ipp].eqother[ido+2*ncomp+nstatv]; //suction pressure increment ds (positive)
  unsatvar[2]= Mm->givenonmechq(saturation_degree, ipp);            // degree of saturation S_r
*/

  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
  {
    unsatvar[0] = Mm->ip[ipp].eqother[ido+2*ncomp+1];   // suction from the previous time step
    unsatvar[1] = suc_fn.getval(Mp->time)-unsatvar[0];  //suction pressure increment ds (positive)
  }
  else
  {
    unsatvar[0] = 200.0;  //constant suction pressure s (positive) (from the previous time step)
    unsatvar[1] = 0.0;    // suction pressure increment
  }
  unsatvar[2]= 0.0;   // degree of saturation S_r

  statev[1]=unsatvar[0]; // suction pressure s (positive)
  statev[2]=unsatvar[2]; // degree of saturation S_r


  // calculate strain increment = actual strain - previous strain and store it in the eps temporarily
  for (i=0; i<ncomp; i++)
    eps[i] = Mm->ip[ipp].strain[i] - Mm->ip[ipp].eqother[ido+i];

  // stress array have to contain values from the last attained equilibrium state
  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].stress[i] = Mm->ip[ipp].eqother[ido+ncomp+i];


  give_full_vector(sig.a, Mm->ip[ipp].stress, ssst);
  give_full_vector(deps.a, eps.a, ssst);
  nullv(eps);
  give_full_vector(eps.a, Mm->ip[ipp].eqother+ido, ssst);
  //  give_full_vector(eps.a, Mm->ip[ipp].strain, ssst);


  // stress vector must be reoredred to ABAQUS notation
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig[3];
  sig[3] = sig[5];
  sig[5] = aux;

  // strain vector must be reoredred to ABAQUS notation
  // SIFEL ordering is eps11, eps22, eps33, eps23, eps13, eps12
  // ABAQUS ordering is eps11, eps22, eps33, eps12, eps13, sig23
  aux = eps[3];
  eps[3] = eps[5];
  eps[5] = aux;

  // increment of strain vector must be reoredred to ABAQUS notation
  // SIFEL ordering is eps11, eps22, eps33, eps23, eps13, eps12
  // ABAQUS ordering is eps11, eps22, eps33, eps12, eps13, sig23
  aux = deps[3];
  deps[3] = deps[5];
  deps[5] = aux;

  if (Mp->timecon.tct == notct)
  {
    dtime = Mp->dlambda;
    time = Mp->lambda;
  }
  else
  {
    dtime = Mp->dtime;
    time = Mp->time;
  }

  umatunsat(ipp, sig, statev, dd, eps, deps, dtime, pnewdt, unsatvar);

  // stress vector must be reoredred to SIFEL notation
  // SIFEL ordering is sig11, sig22, sig33, sig23, sig13, sig12
  // ABAQUS ordering is sig11, sig22, sig33, sig12, sig13, sig23
  aux = sig[3];
  sig[3] = sig[5];
  sig[5] = aux;

  // store stress in reduced notation to the integration point
  give_red_vector(sig.a, Mm->ip[ipp].stress, ssst);

  // store actual value of strain vector
  Mm->storeother(ipp, ido, ncomp, Mm->ip[ipp].strain);
  // store actual value of stress vector
  Mm->storeother(ipp, ido+ncomp, ncomp, Mm->ip[ipp].stress);
  // store actual value of state variables
  Mm->storeother(ipp, ido+2*ncomp, nstatv, statev.a);
  // store actual value of suction
  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    Mm->ip[ipp].other[ido+2*ncomp+1] = suc_fn.getval(Mp->time);
  // store actual value of time step coefficient
  Mm->ip[ipp].other[ido+2*ncomp+nstatv]= pnewdt;
}



/**
  The function returns time step size required by the hypoplsaticity model. It is represented by ratio 
  of required time step size to the actual one or 1.0 for no change in time step size.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns required time step size provided by the hypoplasticity model

  Created by Tomas Koudelka, 10.3.2014
*/
double hypoplast::dstep_red(long ipp, long /*im*/, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  return (Mm->ip[ipp].other[ido+2*ncomp+nstatv]);
//  return 1.0;
}



/**
  Function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2.2014
*/
void hypoplast::updateval (long ipp,long im,long ido)
{
  long i;
  long n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
  Created by Tomas Koudelka,  12.2.2014
*/
void hypoplast::give_reqnmq(long *anmq)
{
  anmq[saturation_degree-1] = 1;
  anmq[suction-1] = 1;
}



/**
  The function extracts virgin porosity e_lambda1 for the attained equilibrium state
  from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of virgin porosity e_lambda1.

  Created by Tomas Koudelka, 13.2.2014  
*/
double hypoplast::give_virgporosity(long /*ipp*/, long /*ido*/)
{
//  return n_star;
  return exp(n_star)-1.0;
}



/**
  The function extracts initial porosity e_ini for the attained equilibrium state
  from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of initial porosity e_ini.

  Created by Tomas Koudelka, 13.2.2014  
*/
double hypoplast::give_iniporosity(long ipp, long /*ido*/)
{
  return Mm->ic[ipp][0];
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,  30.1.2014
*/
void hypoplast::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      phi=val[i];
      break;
    }
    case 1:{
      p_t=val[i];
      break;
    }
    case 2:{
      lam_star=val[i];
      break;
    }
    case 3:{
      kap_star=val[i];
      break;
    }
    case 4:{
      n_star=val[i];
      break;
    }
    case 5:{
      rr=val[i];
      break;
    }
    case 6:{
      n=val[i];
      break;
    }
    case 7:{
      l=val[i];
      break;
    }
    case 8:{
      m=val[i];
      break;
    }
    case 9:{
      s_e0=val[i];
      break;
    }
    case 10:{
      e_0=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}












void hypoplast::umatunsat(long ipp, vector &stress, vector &statev, matrix &ddsdde,
                          vector &/*stran*/, vector &dstran, double dtime, 
                          double &pnewdt, vector &unsatvar)
{
//    ----------------------------------------------------------------------
//
//    Solution dependent state variables (statev):
//    definition via sdvini
//
//      1 ... void    void ratio
//      2 ... s    suction
//      3 ... S_r       degree of saturation
//      4 ... nfev    number of function evaluation
//      5 ... phi_mob phi_mob in degrees
//      6 ... dtsub   suggested substep size

  const long maxnint = 10000;     /* maximum number of time substeps allowed - if the limit is exceeded abaqus is forced to reduce 
                                     the overall time step size (cut-back) */
  const long maxninttest = 1000;

  const double tolintt = 1.0e-3;     // prescribed error tolerance for the adaptive substepping scheme
  const double tolintttest = 1.0e-1;
  const double dtmin = 1.0e-17;   /* minimum substeps size allowed - if the limit is exceeded SIFEL 
                                     is forced to reduce the overall time step size (cut-back) */
  const double perturb = 1.0e-5;  // perturbation parameter for numerical computation of Jacobian matrices
  const double youngel=100.0; // Young modulus used in the stress states close to tension
  const double nuel=0.48;     // Poissons ratio used in the stress states close to tension
  const long nasv = 2;  /* define number of additional state variables
                             asv(0) ... void    void ratio
                             asv(1) ... suction suction */
  const double  nyact = 6 + nasv;
  const long nfasv = 0;   // number of first additional state variable in statev field 
  const long nasvdim = 6; // maximum number of additional state variables ???!!! may be removed
  const long nydim = 6+nasvdim;



  long i;
  long testing, inittension;
  long nfev; // number of function evaluation
  long error;
  double depsv_np1, dtsub, pore, dsuction, norm_deps, norm_deps2, theta;
  vector sig_n(ASTCKVEC(6)), deps_np1(ASTCKVEC(6)), y(ASTCKVEC(nydim)), y_n(ASTCKVEC(nydim)), asv(ASTCKVEC(nasv));
  matrix ddtan(ASTCKMAT(6,6));

  
  // suggested time substep size, and initial excess pore pressure
  dtsub = statev(nasv+3);
  pore = 0.0;

  // vector of additional state variables
  rcopyv(statev, nfasv, asv, 0, nasv);

  // compute volume strain increment and current net stress tensor

  // compute effective stress components sig' = sig_tot + pore
  // NOTE: stress = total stress tensor (tension positive)
  //       pore   = exc. pore pressure (undrained conds., compression positive)
  //       sig    = effective stress (tension positive) 
  copyv(stress, sig_n);
  sig_n(0) += pore;
  sig_n(1) += pore;
  sig_n(2) += pore;

  // Move strain increment dstran into deps and compute volumetric strain increment
  // NOTE: all strains negative in compression; deps has always 6 components
  copyv(dstran, deps_np1);
  depsv_np1 = deps_np1(0) + deps_np1(1) + deps_np1(2);

  dsuction = unsatvar(1);

  // --------------------
  //   Time integration
  // --------------------

  // the following code replaces procedure iniy_hu(y,nydim,nasv,ntens,sig_n,asv);
  rcopyv(sig_n, 0, y, 0, 6);
  rcopyv(asv, 0, y, 6, nasv);

  // the following code replaces procedure push_hu(y,y_n,nydim);
  copyv(y, y_n);

  // ... check whether the initial state is not tensile
  inittension = check_rkf(y);

  // ... parameter of the numerical differentiation: sqrt(macheps)*||deps||
  norm_deps2  = tensor_dbldot_prod(deps_np1, deps_np1, 0.5);
  norm_deps   = sqrt(norm_deps2);
  // ... check whether the strain rate from the ABAQUS is not NAN          
  // norm_D need not to be checked because there is no source of nan at the material 
  // check_isnan(norm_d);

  theta=-perturb*max2(norm_deps,1.0e-6);  // negative sign for compression

  // ... local integration using adaptive RKF-23 method, consistent Jacobian and error estimation
  if((dtsub <= 0.0) || (dtsub > dtime))
    dtsub = dtime;

  // if testing==1 PLAXIS is testing for the initial strain increment.
  testing=0;

  // For use in ABAQUS, comment out the following line
  // ???!!!
  if ((Mp->istep == 1) && (Mp->jstep == 1))     testing=1;

  if (norm_deps == 0.0)      testing = 2;
  //     FEM asking for ddsdde only

  nfev = 0; // initialisation

  if (inittension == 0)
  {
    if (testing == 1) // calculate stresses
    {
      error = rkf23_update(ipp, y, nasv, dtsub, tolintttest, maxninttest, 
                           dtmin, deps_np1, dsuction, nfev, dtime);
      // give original state if the model fails without substepping
      if (error == 3)
      {
        for(i=0; i<nyact; i++)
          y(i) = y_n(i);

        if(unsatvar(1) != 0.0)  // suction pressure increment
          y(6+1) = y_n(6+1) + unsatvar(1);
        error = 0;
      }
    }
    else
    { 
      if (testing == 2)
      {
        for(i=0; i<nyact; i++)
          y(i) = y_n(i);

        if(unsatvar(1) != 0.0) 
          y(6+1) = y_n(6+1)+unsatvar(1);
        
        error = 0;
      }
      // Normal RKF23 integration
      else // (inittension == 0) && (testing == 0)
      {
        error = rkf23_update(ipp, y, nasv, dtsub, tolintt, maxnint, dtmin,
                             deps_np1, dsuction, nfev, dtime);
      }
    }
    // error conditions (if any)

    if (error == 3)
    {
      // reduce abaqus load increment
      pnewdt= 0.25;
      // do not do anything, we are the most likely close to the tensile region
      // y(i)=y_n(i)
      rcopyv(y_n, 0, y, 0, long(nyact));
    }
    else
    {
      if (error == 10)
      {
        print_err("severe error in hypoplastic model (element %ld, ipp=%ld)", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
        abort();
      }
    }

    // update dtsub and nfev
    if(dtsub <= 0.0)
      dtsub = 0.0;
    else
    { 
      if (dtsub >= dtime)
        dtsub = dtime;
    }
    statev(nasv+3) = dtsub;
    statev(nasv+1) = double(nfev);

    // compute consistent tangent via numerical perturbation
    perturbate(y_n, nasv, deps_np1, ddtan);
  }
  else // inittension != 0
  {
    // we were initilly in the tensile stress, calc elastic
    calc_elasti(y, deps_np1, ddtan, youngel, nuel);

    if(unsatvar(1) != 0.0)   
      y(6+1) = y_n(6+1)+unsatvar(1);
  } 

  // convert solution (stress + cons. tangent) to abaqus format
  // update pore pressure and compute total stresses 
  solout(stress, asv, ddsdde, y, ddtan);
     

  // updated vector of additional state variables to abaqus statev vector
  // statev(i+nfasv) = asv(i);
  rcopyv(asv, 0, statev, nfasv, nasv);
  
  if (inittension == 0)
    calc_statev(stress, statev, nasv);

  // -----------------------
  // End of time integration
  // -----------------------

  // Temporary WRC for testing purposes

  unsatvar(0) = unsatvar(0) + unsatvar(1);
  unsatvar(2) = pow(10.0/unsatvar(0), 0.38);
  unsatvar(3) = unsatvar(2) - pow(10.0/(unsatvar(0)-unsatvar(1)), 0.38);
  unsatvar(4) = unsatvar(3)/unsatvar(1);

  return;
}



/**
  The function checks whether the solution vector y is OK for hypoplasticity.
  
  @param y - %vector of function values (it contains stress components)

  @retval 0 - on success
  @retval 1 - if the tensile stresses detected

  Created by Tomas Koudelka 20.4.2015 according to C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff code
*/
long hypoplast::check_rkf(vector &y)
{
  long iopt;

  double pmean;
  vector sig_star(ASTCKVEC(6));

  double tmin;
  matrix sig_t(ASTCKMAT(3,3)), pvect(ASTCKMAT(3,3));
  vector pval(ASTCKVEC(3));
  double minstress, suction, zero;

  zero=0.0;

  minstress=p_t/4.0;

  rcopyv(y, 0, sig_star, 0, 6);

  sig_star(0) -= p_t;
  sig_star(1) -= p_t;
  sig_star(2) -= p_t;

  suction=y(7);
  
  calc_khalili_stress(sig_star, sig_star, suction);
              
  pmean =-(1.0/3.0)*(sig_star(0)+sig_star(1)+sig_star(2));
            
  // check for positive mean stress
  if(pmean <= minstress)
    return 1;

  // calculate minimum principal stress

  iopt=0;
  vector_tensor(sig_star, sig_t, spacestress, stress);
  princ_val(sig_t, pval, pvect, nijac, limit, Mp->zero, 3, 1);
  tmin = -pval(0);
     
  // check for tension
  if(tmin <= minstress) 
    return 1;
        
  // check for NAN
  // y vector should be checked for NAN but is performed at the element level
  if (test_math_err())
    return 1;

  return 0;
}      



/**
  The function calculates effective stress %vector according to Khalili.

  @param sig_net - total stress %vector
  @param sig_ef  - effective stress %vector
  @param suction - suction pressure

  Created by Tomas Koudelka 20.4.2015, according to C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff code
*/
void hypoplast::calc_khalili_stress(vector &sig_net, vector &sig_ef, double suction)
{
  double khalfact;
  double gamma = 0.55;

  if(suction > s_e0)
    khalfact=-suction*pow(s_e0/suction, gamma);
  else
    khalfact=-suction;

  copyv(sig_net, sig_ef);
  sig_ef(0) += khalfact;
  sig_ef(1) += khalfact;
  sig_ef(2) += khalfact;

  return;
}      


/**
  The function performs numerical solution of y'=f(y) by
  explicit, adapive RKF23 scheme with local time step extrapolation

  @param y - (input/output)
  @param 
  @param 
  @param 

  @retval 0 - on success
  @retval 3 - on RKF error 
 
  Created by Tomas Koudelka 21.4.2015 according to Tamagnini, Sellari & Miriano code 
*/
long hypoplast::rkf23_update(long ipp, vector &y, long nasv, double &dtsub, double err_tol, 
                     long maxnint, double dtmin, vector &deps_np1,
                     double dsuction, long &nfev, double dtime)
{
  long i, ksubst, kreject;
  long error_rkf = 0;
  long ny = y.n;

  vector y_k(ASTCKVEC(ny)), y_2(ASTCKVEC(ny)), y_3(ASTCKVEC(ny)), y_til(ASTCKVEC(ny)), y_hat(ASTCKVEC(ny));
  vector krk_1(ASTCKVEC(ny)), krk_2(ASTCKVEC(ny)), krk_3(ASTCKVEC(ny));
  double t_k, dt_k;
  double norm_r, s_hull;
  double temp;


  // start of update process                
  t_k = 0.0;
  dt_k = dtsub/dtime;
  ksubst = 0;
  kreject = 0;
  nfev = 0;
  copyv(y, y_k);

  // start substepping 
  while(t_k < 1.0) 
  {
    ksubst=ksubst+1;
    // check for maximum number of substeps
    if (ksubst > maxnint)
    {
      print_warning("number of substeps %ld is too big, step rejected (elem=%ld, ipp=%ld)",
                    __FILE__, __LINE__, __func__, ksubst, Mm->elip[ipp]+1, ipp+1);
      return 3;
    }
    // build RK functions
    error_rkf = check_rkf(y_k);
    if (error_rkf == 1)
    {
      return 3;
    }
    else
      rhs(y_k, nasv, deps_np1, krk_1, nfev, dsuction);

    // find y_2
    temp=0.5*dt_k;

    // y_2(i)=y_k(i)+temp*krk_1(i)
    addmultv(y_k, krk_1, temp, y_2);

    error_rkf = check_rkf(y_2);
    if (error_rkf == 1)
    {
      return 3;
    }
    else
      rhs(y_2, nasv, deps_np1, krk_2, nfev, dsuction);

    // find y_3
    for(i=0; i<ny; i++)
      y_3(i)=y_k(i)-dt_k*krk_1(i) + 2.0*dt_k*krk_2(i);
   
    error_rkf = check_rkf(y_3);
    if (error_rkf == 1)
    {
      return 3;
    }
    else
      rhs(y_3, nasv, deps_np1, krk_3, nfev, dsuction);

    // approx. solutions of 2nd (y_til) and 3rd (y_hat) order
    for(i=0; i<y_til.n; i++)
    {
      y_til(i)=y_k(i)+dt_k*krk_2(i);
      y_hat(i)=y_k(i)+dt_k*((1.0/6.0)*krk_1(i) + (2.0/3.0)*krk_2(i) + (1.0/6.0)*krk_3(i));
    }

    // local error estimate
    norm_r = norm_res(y_til,y_hat, nasv);

    // check if output y_hat can be used as an input into the next step
    error_rkf = check_rkf(y_hat);
           
    if (error_rkf)
    {
      return 3;
    }

    // time step size estimator according to Hull
   
    if (norm_r != 0)
      s_hull = (0.9*dt_k)*pow((err_tol/norm_r), 1.0/3.0);
    else
      s_hull = 1.0;

    if (norm_r < err_tol)
    {
      // substep is accepted, update y_k and t_k and estimate new substep size dt_k
      copyv(y_hat, y_k);
      t_k  = t_k + dt_k;
      dt_k = min2(4.0*dt_k, s_hull);
      if (dt_k > 1.0)
        dt_k = 1.0;
      dtsub = dt_k*dtime;
      dt_k = min2((1.0 - t_k), dt_k);
    }
    else
    {
      // substep is not accepted, recompute with new (smaller) substep size DT
      dt_k = max2(dt_k/4.0, s_hull);
      // check for minimum step size
      if (dt_k < dtmin) 
      {
        return 3;
      }
    }
  } // end of while loop
        
  // recover final state
  copyv(y_k, y);

  return 0;
}



/**
  The function calculates coefficient krk from the current state y 
  and strain increment deps for Masin hypoplastic model for clays 
  with intergranular strains.

  @param y - 
  @param nasv - 
  @param deps - %vector of strain increment
  @param krk - %vector of coefficients for Runge-Kutta method (output)
  @param nfev - the number of function evaluations
  @param dsuction - suction increment

  @return

  Created by Tomas Koudelka, 21.4.2015 according to Tamagnini & Sellari
*/
void hypoplast::rhs(vector &y, long nasv, vector &deps, vector &krk, 
                    long &nfev, double dsuction)
{
  vector sig(ASTCKVEC(6)), q(ASTCKVEC(nasv));
  vector f_sig(ASTCKVEC(6)), f_q(ASTCKVEC(nasv));

  // update counter for the number of function f(y) evaluations
  nfev=nfev+1;

  // initialize kRK
  fillv(0.0, krk);

  // recover current state variables (sig,q)                   
  rcopyv(y, 0, sig, 0, 6);
  rcopyv(y, 6, q, 0, nasv);

  // build f_sig(6) and f_q(nasv) vectors and move them into kRK
  get_f_sig_q(sig, q, deps, f_sig, f_q, dsuction);

  rcopyv(f_sig, 0, krk, 0, 6);
  rcopyv(f_q, 0, krk, 6, nasv);

  return;
}



/**
  The function finds vectors f_sigma and f_q in f(y)

  @param sig - stress %vector
  @param q -
  @param deps - %vector of strain increment
  @param f_sig -  (output)
  @param f_q -  (output)
  @param dsuction - suction increment

  @retval 0 - on success
  @retval 1 - 
  @retval 1 - 
  @retval 1 - 
  @retval 10 - 

  Created by Tomas Koudelka, 21.4.2015, according to (Tamagnini, Sellari & Miriano)
*/
void hypoplast::get_f_sig_q(vector &sig, vector &q, vector &deps,
                            vector &f_sig, vector &f_q, double dsuction)
{
  long nasv = q.n;

  matrix hh(ASTCKMAT(nasv,6)), ll(ASTCKMAT(6,6));
  vector nn(ASTCKVEC(6)), hunsat(ASTCKVEC(6));
  double norm_d,norm_d2;
  double suction, zero, gamma;
  double khalratefact;
      
  suction = q(1);
  zero=0.0;
      
  // compute tangent operators
  get_tan(deps, sig, q, hh, ll, nn, hunsat);

  // compute F_sig=MM*deps
  mxv(ll, deps, f_sig);
  norm_d2 = tensor_dbldot_prod(deps, deps, 0.5);
  norm_d  = sqrt(norm_d2);
  addmultv(f_sig, nn, norm_d);

  khalratefact=-dsuction;
  if (suction > s_e0)
  {
    gamma=0.55;
    khalratefact=-(1-gamma)*dsuction*pow((s_e0/suction),gamma);
  }

  f_sig(0) -= khalratefact;
  f_sig(1) -= khalratefact;
  f_sig(2) -= khalratefact;

  // wetting collapse contribution                      
  if (suction > s_e0)
  {
    if (dsuction < zero)
      addmultv(f_sig, hunsat, dsuction);
  }
 
  // compute F_q=HH*deps
  mxv(hh, deps, f_q);
  // f_q(0)=HH(0,0)*(deps(0)+deps(1)+deps(2))
  f_q(1)=dsuction;

  return;
}



/** 
  The function computes matrices M and H for Masin and Khalili (2008) model for unsaturated soils.

  @param deps - %vector(6) of strain increments
  @param sig - %vector(6) of
  @param q - %vector(nasv) of
  @param hh - %matrix(nasv,6) of 
  @param ll - %matrix(6,6) of (input)
  @param nn - %vector(6) of
  @param hunsat - %vector(6) of

  @return

  Created by Tomas Koudelka 22.4.2015 according to 
*/
void hypoplast::get_tan(vector &/*deps*/, vector &sig, vector &q, matrix &hh, matrix &ll, vector &nn, vector &hunsat)
{
  const double sqrt3=sqrt(3.0), twosqrt2=2.0*sqrt(2.0), sqrt2=sqrt(2.0), ln2m1=1.0/log(2.0);

  long i, j;
  vector eta(ASTCKVEC(6)), eta_dev(ASTCKVEC(6)), sig_star(ASTCKVEC(6));

  matrix ii(ASTCKMAT(6,6)), aa(ASTCKMAT(6,6)), aainv(ASTCKMAT(6,6));
  vector m(ASTCKVEC(6)), ainvn(ASTCKVEC(6)), nnpure(ASTCKVEC(6));
  vector m_dir(ASTCKVEC(6)), m_dir1(ASTCKVEC(6)), h_e(ASTCKVEC(6));
  vector ikron(ASTCKVEC(6));

  double voidr;
  double fdsbs;
  double pp, qq, cos3t, i1, i2, i3, tanpsi;
  double a, a2, ff, alpha, fd, fdi, fs, c_1, c_2, yi, yy;
  double num, den, af, fa2, eta_dn2, eta_n2, norm_m, norm_m2;
  double dpeds_divpe;
  double logsse, fufact, pehvor;
  double temp1, temp2, suction;
  double lam_star_n, n_star_n, r_lc, p_ref;
  double sinphi, sinphi2;

  fillrow(0.0, 0, hh);
     
  // fourth order identity tensors in Voigt notation
  ii(0,0)=1.0;
  ii(1,1)=1.0;
  ii(2,2)=1.0;
  ii(3,3)=0.5;
  ii(4,4)=0.5;
  ii(5,5)=0.5;

  ikron(0)=1.0;
  ikron(1)=1.0;
  ikron(2)=1.0;

  p_ref=1.0;

  sinphi  = sin(phi);
  sinphi2 = sinphi*sinphi;
      
  // recover internal state variables
  voidr=q(0);
  suction=q(1);

  // axis translation due to cohesion (p_t>0)
  copyv(sig, sig_star);
  sig_star(0) -= p_t;
  sig_star(1) -= p_t;
  sig_star(2) -= p_t;
      
  // calculate Khalili stress
  calc_khalili_stress(sig_star, sig_star, suction);
      
  // calculate n_star_n and lambda_star_n using current suction
  r_lc = rr;
      
  if (suction > s_e0)
    logsse = log(suction/s_e0);
  else
    logsse = 0.0;

  n_star_n   = n_star   + this->n*logsse;
  lam_star_n = lam_star + this->l*logsse;

  // auxiliary stress tensors
  inv_sig(sig_star, pp, qq, cos3t, i1, i2, i3);

  cmulv(1/i1, sig_star, eta);
  copyv(eta, eta_dev);
  eta_dev(0) -= 1.0/3.0;
  eta_dev(1) -= 1.0/3.0;
  eta_dev(2) -= 1.0/3.0;

  // functions a and F
  eta_dn2 = tensor_dbldot_prod(eta_dev, eta_dev, 2.0);
  tanpsi  = sqrt3*sqrt(eta_dn2);
  temp1   = (1.0/8.0*tanpsi*tanpsi) + 
            (2.0-tanpsi*tanpsi)/(2.0+sqrt2*tanpsi*cos3t);
  temp2   = tanpsi/twosqrt2;

  a  = sqrt3*(3.0-sinphi)/(twosqrt2*sinphi);
  a2 = a*a;
  ff = sqrt(temp1)-temp2;

  // barotropy and pyknotropy functions
  temp1 = (lam_star_n - kap_star) / (lam_star_n + kap_star);
  temp2 = (3.0 + a2) / (sqrt3 * a);
  alpha = ln2m1 * log(temp1*temp2);

  temp1 = (2.0*i1)/(3.0*p_ref);
  pehvor = exp((n_star_n - log(1.0+voidr))/lam_star_n);
  fd = pow(-temp1/pehvor, alpha);

  fdi = pow(2.0, alpha);
  temp1 = -i1/lam_star_n;
  temp2 = 3.0+a2-fdi*sqrt3*a;
  fs = temp1/temp2;
      
  // tensor L
  temp1 = 2.0/(9.0*r_lc);
  c_1 = temp1*temp2;
  c_2 = 1.0+(1.0-c_1)*3.0/a2;

  for(i = 0; i < 6; i++)
  {
    for(j = 0; j < 6; j++)
      ll(i,j)=(3.0*c_1)*ii(i,j) + (3.0*c_2*a2)*(eta(i)*eta(j));
  }

  // function YY
  yi = (sqrt3*a)/(3.0+a2);
  num = (i1*i2+9.0*i3)*(1.0-sinphi2);
  den = 8.0*i3*sinphi2;
  yy = (yi-1.0)*(num/den)+yi;

  // tensor m and NN
  af = a/ff;
  fa2 = ff*ff/a2;
  eta_n2 = tensor_dbldot_prod(eta, eta, 2.0);
  temp1 = (1.0/3.0)*(6.0*eta_n2-1.0)/(fa2+eta_n2);

  for(i=0; i<6; i++)
    m(i) = -af*((eta(i)+eta_dev(i)) - (temp1*eta(i)));

  norm_m2 = tensor_dbldot_prod(m, m, 2.0);
  norm_m  = sqrt(norm_m2);
  cmulv(1.0/norm_m, m, m_dir);

  cmulv(-yy, m_dir, m_dir1);
  m_dir1(3) *= 2.0;
  m_dir1(4) *= 2.0;
  m_dir1(5) *= 2.0;

  mxv(ll, m_dir1, nnpure);
  cmulm(fs, ll);
  cmulv(fs*fd,nnpure, nn);

  fillv(0.0, hunsat);        
      
  if (suction > s_e0)
  {
    dpeds_divpe = (this->n - log(pehvor)*this->l)/(lam_star_n*suction);
    cmulv(dpeds_divpe, sig_star, hunsat);

    for(i=0; i<6; i++)
    {
      for(j=0; j<6; j++)
        aa(i,j) = ll(i,j)+sig_star(i)*ikron(j)/lam_star_n;
    }
         
    invm(aa, aainv, Mp->zero);
    mxv(aainv, nnpure, ainvn);
    cmulv(fs, ainvn);

    norm_m2 = tensor_dbldot_prod(ainvn, ainvn, 2.0);
    norm_m  = sqrt(norm_m2);
    fdsbs   = 1.0/norm_m;        
    fufact  = pow(fd/fdsbs, this->m/alpha);
    cmulv(fufact, hunsat);         
  }

  // void ratio evolution function (tension positive)
  h_e(0) = 1.0 + voidr;
  h_e(1) = 1.0 + voidr;
  h_e(2) = 1.0 + voidr;

        
  for(j=0; j<6; j++)
    hh(0,j) = h_e(j);
      
  return;
}



/**
  The function calculates stress invariants from the given stress %vector sig
  where the Voigt notation is used with the following index conversion
  11 -> 1
  22 -> 2
  33 -> 3
  12 -> 4
  13 -> 5
  23 -> 6

  @param sig - stress %vector in Voigt notation (input)
  @param pp - mean stress
  @param qq - sqrt(3.0*J2s) where J2s is the second invariant of stress deviator
  @param cos3t - Lodes angle
  @param i1 - the first invariant of stress tensor
  @param i2 - the second invariant of stress tensor
  @param i3 - the third invariant of stress tensor
*/
void hypoplast::inv_sig(vector &sig, double &pp, double &qq, double &cos3t, double &i1, double &i2, double &i3)
{
  vector sdev(ASTCKVEC(6)), eta(ASTCKVEC(6)), eta_d(ASTCKVEC(6)), eta_d2(ASTCKVEC(6));
  double xmin1,xmin2,xmin3;
  double tretadev3;
  double norm2, norm2sig, norm2eta, numer, denom;

  const double sqrt6=sqrt(6.0);

  // trace and mean stress
  i1 = first_invar(sig);
  pp = (1.0/3.0)*i1;
  deviator(sig, sdev);

  // normalized stress and dev. normalized stress
  cmulv(1.0/i1, sig, eta);
  copyv(eta, eta_d);
  eta_d(0) -= 1.0/3.0;
  eta_d(1) -= 1.0/3.0;
  eta_d(2) -= 1.0/3.0;

  // second invariants
  norm2    = tensor_dbldot_prod(sdev, sdev, 2.0);
  norm2sig = tensor_dbldot_prod(sig, sig, 2.0);
  norm2eta = tensor_dbldot_prod(eta_d, eta_d, 2.0);

  qq = sqrt(1.5*norm2);
  i2 = 0.5 * (norm2sig - i1*i1);

  // components of (eta_d_ij)(eta_d_jk)
  eta_d2(0) = eta_d(0)*eta_d(0) + eta_d(3)*eta_d(3) + eta_d(4)*eta_d(4);
  eta_d2(1) = eta_d(3)*eta_d(3) + eta_d(1)*eta_d(1) + eta_d(5)*eta_d(5);
  eta_d2(2) = eta_d(5)*eta_d(5) + eta_d(4)*eta_d(4) + eta_d(2)*eta_d(2);
  eta_d2(3) = eta_d(0)*eta_d(3) + eta_d(3)*eta_d(1) + eta_d(5)*eta_d(4);
  eta_d2(4) = eta_d(4)*eta_d(0) + eta_d(5)*eta_d(3) + eta_d(2)*eta_d(4);
  eta_d2(5) = eta_d(3)*eta_d(4) + eta_d(1)*eta_d(5) + eta_d(5)*eta_d(2);

  // Lode angle
  if (norm2eta < 1.0e-18)
    cos3t = -1.0;
  else
  {
    tretadev3 =  tensor_dbldot_prod(eta_d, eta_d2, 2.0);
    numer     = -sqrt6 * tretadev3;
    denom     = pow(sqrt(norm2eta), 3.0);
    cos3t     = numer/denom;
    if(fabs(cos3t) > 1.0) 
      cos3t = cos3t/fabs(cos3t);
  }

  // determinant
  xmin1 = sig(1)*sig(2) - sig(5)*sig(5);
  xmin2 = sig(3)*sig(2) - sig(5)*sig(4);
  xmin3 = sig(3)*sig(5) - sig(4)*sig(1);

  i3 = sig(0)*xmin1 - sig(3)*xmin2 + sig(4)*xmin3;
}



/**
  The function evaluates norm of residual vector Res=||y_hat-y_til||

  Created by Tomas Koudelka, 22.4.2015 according to Tamagnini, Sellari & Miriano
*/
double hypoplast::norm_res(vector &y_til, vector &y_hat, long nasv)
{
  long i, testnan;
  long ny = y_til.n;

  double void_til, void_hat, del_void;
  double norm_r2, norm_r;
  double norm_sig2, norm_sig;

  vector sig_hat(ASTCKVEC(6)), sig_til(ASTCKVEC(6)), del_sig(ASTCKVEC(6));
  vector err(ASTCKVEC(ny));
  
  //recover stress tensor and internal variables
  rcopyv(y_hat, 0, sig_hat, 0, 6);
  rcopyv(y_til, 0, sig_til, 0, 6);
  for(i=0; i<6; i++)
    del_sig(i) = fabs(sig_hat(i)-sig_til(i));
  
  void_hat = y_hat(6);
  void_til = y_til(6);
  del_void = fabs(void_hat-void_til);

  // relative error norms

  norm_sig2 = tensor_dbldot_prod(sig_hat, sig_hat, 2.0);
  norm_sig  = sqrt(norm_sig2);

  if (norm_sig > 0.0)
  {
    double aux = 1.0/norm_sig;
    for(i=0; i<6; i++)
      err(i) = del_sig(i)*aux;
  }

  err(6+nasv-2) = del_void/void_hat;
      
  // global relative error norm
  scprd(err, err, norm_r2);
  norm_r=sqrt(norm_r2);

  testnan=0;
  //testnan += check_isnan(norm_sig);
  //testnan += check_isnan(void_hat);
  if (testnan)
    norm_r = 1.0e20;
      
  return norm_r;
}



/**
  The function computes consistent tangent stiffness numerically.

  Created by Tomas Koudelka, 22.4.2015, according to Tamagnini
*/
void hypoplast::perturbate(vector &y_n, long nasv, vector &deps_np1, matrix &dd)
{
  long error_rkf=0;
  vector sig(ASTCKVEC(6)), q(ASTCKVEC(nasv));
  vector nn(ASTCKVEC(6)), hunsattmp(ASTCKVEC(6));
  matrix hhtmp(ASTCKMAT(nasv,6)), ll(ASTCKMAT(6,6));

  // initialize dd and y_star
  fillm(0.0, dd);      
  rcopyv(y_n, 0, sig, 0, 6);
  rcopyv(y_n, 6, q, 0, nasv);

  if (error_rkf != 10) 
    get_tan(deps_np1, sig, q, hhtmp, dd, nn, hunsattmp);

  return;
}



/**
   The function performs numerical solution of y'=f(y)
   explicit, adapive RKF23 scheme with local time step extrapolation.

   Created by Tomas Koudelka, 22.4.2015, according to Tamagnini, Sellari & Miriano
*/
void hypoplast::calc_elasti(vector &y, vector &deps_np1, matrix &ddtan, 
                            double youngel, double nuel)
{
  long i, j;

  vector krondelta(ASTCKVEC(6)), f_sig(ASTCKVEC(6));
  matrix ii(ASTCKMAT(6,6));

  // initialize variables

  // fourth order identity tensors in Voigt notation
  ii(0,0) = 1.0;
  ii(1,1) = 1.0;
  ii(2,2) = 1.0;
  ii(3,3) = 0.5;
  ii(4,4) = 0.5;
  ii(5,5) = 0.5;

  krondelta(1) = 1.0;
  krondelta(2) = 1.0;
  krondelta(3) = 1.0;

  // Elastic stiffness tensor 
  for(i=0; i<6; i++)
  {
    for(j=0; j<6; j++)
      ddtan(i,j) = (youngel/(1+nuel))*
                   (ii(i,j) + (nuel/(1-2*nuel))*(krondelta(i)*krondelta(j)));
  }
      
  mxv(ddtan, deps_np1, f_sig);

  for (i=0; i<6; i++)
    y(i) = y(i) + f_sig(i);

  return;
}



/**
  The function copies the vector of state variables to umat output

  NOTE: solid mechanics convention for stress and strain components
  pore is always positive in compression

  Created by Tomas Koudelka, 22.4.2015 according to  Tamagnini, Sellari
*/
void hypoplast::solout(vector &stress, vector &asv, matrix &ddsdde, vector &y, matrix &dd)
{
  // updated total stresses (effective stresses stored in y(0:5))
  rcopyv(y, 0, stress, 0, 6);

  // additional state variables (void ratio, suction)
  rcopyv(y, 6, asv, 0, asv.n);

  // consistent tangent stiffness
  copym(dd, ddsdde);

  return;
}



/**
  The function computes additional state variables for postprocessing.
*/
void hypoplast::calc_statev(vector &stress, vector &statev, long nasv)
{
  double i1, i2, i3, cos3t, pp, qq;
  double sin2phi, sinphi;
  vector sig_star(ASTCKVEC(6));

  // calc phimob (statev(4)) from Matsuoka-Nakai YS
  copyv(stress, sig_star);
  sig_star(0) -= p_t;
  sig_star(1) -= p_t;
  sig_star(2) -= p_t;

  calc_khalili_stress(sig_star, sig_star, statev(1));

  inv_sig(sig_star, pp, qq, cos3t, i1,i2,i3);
  if (i3 != 0.0)
    sin2phi=(9.0 + i1*i2/i3)/(1.0+i1*i2/i3);
  else 
    sin2phi=0.0;

  if (sin2phi < 0.0)
    sin2phi=0.0;

  if (sin2phi > 1.0)
    sin2phi = 1.0;

  sinphi = sqrt(sin2phi);
  statev(nasv+2) = asin(sinphi)*(180.0/M_PI);
  return;
}
