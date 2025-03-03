#include "hypoplast.h"
#include "global.h"
#include "intpoints.h"
#include "vecttens.h"
#include "mechmat.h"
#include "probdesc.h.h"

#include "iotools.h"
#include "matrix.h"
#include "vector.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef FCC_DEF
extern "C" void umatunsat_(double &, double& , double &, double &,
                            double &, double &, double &, double &, 
                            double &, double &, double &, double &,
                            double &, double &, double &, double &,
                            double &, double &, char &,   int &,  
                            int &,    int &,    int &,    double &,
                            int &,    double &, double &, double &,
                            double &, double &, double &, int &, 
                            int &,    int &,    int &,    int &,
                            int &,    double &);
#else
void umatunsat_(double &, double& , double &, double &,
                double &, double &, double &, double &, 
                double &, double &, double &, double &,
                double &, double &, double &, double &,
                double &, double &,   char &,    int &,  
                   int &,    int &,    int &, double &,
                   int &, double &, double &, double &,
                double &, double &, double &,    int &, 
                   int &,    int &,    int &,    int &,
                   int &, double &)
{
  print_err("Hypoplastic material from ABAQUS was not compiled\n"
            "(Fortran compiler was not found)", __FILE__, __LINE__, __func__);
}
#endif


hypoplast::hypoplast()
{
  phi = p_t = lam_star = kap_star = n_star = 0.0;
  rr = n = l = m = s_e0 = e_0;
  nprops = 11;
  nstatv = 6;
  pr_suc_fl = no;
  props = new double[nprops];
  memset(props, 0, sizeof(*props)*nprops);
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
              "phi", &phi, "p_t", &p_t, "lam_star", &lam_star, 
              "kap_star", &kap_star, "N_star", &n_star);
  xfscanf(in, "%k %le %k %le %k %le %k %le %k %le %k %le",
              "rr", &rr, "n", &n, "l", &l, "m", &m, 
              "s_e0", &s_e0, "e_0", &e_0);
  xfscanf(in, "%k%m", "presc_suction", &answertype_kwdset, &pr_suc_fl);
  if (pr_suc_fl == yes)
    suc_fn.read(in);

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

  props[0]  = phi;
  props[1]  = p_t;
  props[2]  = lam_star;
  props[3]  = kap_star;
  props[4]  = n_star;
  props[5]  = rr;
  props[6]  = n;
  props[7]  = l;
  props[8]  = m;
  props[9]  = s_e0 ;
  props[10] = e_0;
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
  long i;
  int noel = int(Mm->elip[ipp]); // element id
  int npt = int(ipp);
  int ndi = 3; // number of normal stress components
  int nsh = 3; // number of shear stress components
  int ntens = ndi+nsh; // total number of stress components

  // dummy variables
  double ddum=0.0;
  double ddum3[3]={0.0, 0.0, 0.0};
  double ddum6[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dfgrd[3*3]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  int    idum;
  char   cmname []="unsatm";

  // should be updated according METR/TRFEL
  double unsatvar[5]={0.0, 0.0, 0.0, 0.0, 0.0};

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double *statev = Mm->ip[ipp].other+ido+2*ncomp;
  double aux, time, dtime;
  vector sig(6), eps(6);  
  matrix dd(6,6);
  
  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp; i<ido+2*ncomp+nstatv; i++)
    Mm->ip[ipp].other[i] = Mm->ip[ipp].eqother[i];
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

  umatunsat_(*sig.a, *statev, *dd.a,    ddum,   ddum, 
              ddum,   ddum,   *ddum6,  *ddum6,  ddum,
             *eps.a, *ddum6,   time,    dtime,  ddum, 
              ddum,  *ddum6,  *ddum6,  *cmname, ndi, 
              nsh,    ntens,   nstatv, *props,  nprops, 
             *ddum3, *dfgrd,   ddum,    ddum,  *dfgrd, 
             *dfgrd,  noel,    npt,     idum,   idum, 
             (int&)Mp->istep, (int&)Mp->jstep, *unsatvar);

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
void hypoplast::givestressincr (long ipp, long im, long ido, vector &sig)
{
  long i;

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  double ds;
  matrix sigt(3,3);

  if (pr_suc_fl == yes) //prescribed suction pressure s (positive)
    ds = suc_fn.getval(Mp->time)-Mm->ip[ipp].eqother[ido+2*ncomp+1];  //suction pressure increment ds (positive)
  else
    ds = 0.0;

  // calculate stress increment tensor due to change in suction (dsig = -dpp = ds)
  for (i=0;i<3;i++)
    sigt[i][i] += ds;
  
  // convert total stress tensor back to vector form - sig
  //  tensor_vector(sig, sigt, Mm->ip[ipp].ssst,stress);
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
void hypoplast::nlstresses (long ipp, long im, long ido)
{
  long i;
  int noel = int(Mm->elip[ipp]); // element id
  int npt = int(ipp);
  int ndi = 3; // number of normal stress components
  int nsh = 3; // number of shear stress components
  int ntens = ndi+nsh; // total number of stress components

  // dummy/unused variables
  double ddum=0.0;
  double ddum3[3]={0.0, 0.0, 0.0};
  double ddum6[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddum36[6*6];
  memset(ddum36, 0, sizeof(*ddum36)*36);
  double dfgrd[3*3]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  int    idum;
  char cmname []="unsatm";

  // should be updated according METR/TRFEL
  double unsatvar[5]={0.0, 0.0, 0.0, 0.0, 0.0};

  // data handled by MEFEL
  long ncomp  =  Mm->ip[ipp].ncompstr;
  strastrestate ssst =  Mm->ip[ipp].ssst;
  double *statev = Mm->ip[ipp].other+ido+2*ncomp;
  double aux, time, dtime;
  double pnewdt = 1.0; // reduction coefficient for new time step
  vector sig(6), eps(6), deps(6);  
  matrix dd(6,6);



  // prepare working copy of internal variables from the last equilibrium state
  for (i=ido+2*ncomp; i<ido+2*ncomp+nstatv; i++)
    Mm->ip[ipp].other[i] = Mm->ip[ipp].eqother[i];

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

  umatunsat_(*sig.a, *statev, *ddum36,  ddum,   ddum, 
              ddum,   ddum,   *ddum6,  *ddum6,  ddum,
             *eps.a, *deps.a,  time,    dtime,  ddum, 
              ddum,  *ddum6,  *ddum6,  *cmname, ndi, 
             nsh,    ntens,    nstatv, *props,  nprops, 
             *ddum3, *dfgrd,   pnewdt,  ddum,  *dfgrd, 
             *dfgrd,  noel,    npt, idum,   idum, 
             (int&)Mp->istep, (int&)Mp->jstep, *unsatvar);

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
double hypoplast::dstep_red(long ipp, long im, long ido)
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
double hypoplast::give_virgporosity(long ipp, long ido)
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
double hypoplast::give_iniporosity(long ipp, long ido)
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
