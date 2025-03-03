#include "damplastifacemat.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "iotools.h"
#include "vector.h"
#include "matrix.h"
#include "mathem.h"
#include "intpoints.h"
#include <math.h>
#include <stdlib.h>


/**
  Constructor initializes model parameters to the zero values.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
damplastifacemat::damplastifacemat()
{
  ks = kn = 0.0;
  ft = wf = phideg = tanphi = c = 0.0;
  coref = no;
}



/**
  Destructor deallocates used heap memory.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
damplastifacemat::~damplastifacemat()
{
}



/**
  The function reads the model parameter sform the opened text file
  given by the pointer in.

  @param in - pointer to the opened XFILE structure
  
  @return - The function does not return anything but stores new values in the model parameters.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::read(XFILE *in)
{
  xfscanf(in, "%k%le %k%le", "ks", &ks, "kn", &kn);
  xfscanf(in, "%k%le %k%le %k%le %k%le", 
          "ft", &ft, "wf", &wf, "phi_deg", &phideg, "cohesion", &c);
  sra.read(in);

  // read an optional treatment of corrosion effect
  xfscanf(in, "%k%m", "corrosion_effect", &answertype_kwdset, &coref);
  if (coref == yes){
    // evolution function of the corrossion expansion factor alpha_c(t)
    xfscanf(in, "%k", "alpha_c");
    alpha_c.read(in);  

    // evolution function of the corrossion current density i_cor_func(t)
    xfscanf(in, "%k", "i_cor_func");
    icor_func.read(in);

    // evolution function of the flux of corrosion products flux_cor_func(t)
    xfscanf(in, "%k", "flux_cor_func");
    flux_cor_func.read(in);
  }
  tanphi = tan(phideg*M_PI/180);
}



/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened output text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::print(FILE *out)
{
  fprintf(out, "%le %le ", ks, kn);
  fprintf(out, "%le %le %le %le ", ft, wf, phideg, c);
  sra.print(out);

  fprintf(out, " %d", coref);
  if (coref == yes){    
    alpha_c.print(out);  
    icor_func.print(out);
    flux_cor_func.print(out);    
  }
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  @return The function returns computed elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::elmatstiff(matrix &d, long ipp)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  for (i=0; i<nc-1; i++){
    d(i,i) = ks;
  }
  d(nc-1,nc-1) = kn;
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns computed stiffness %matrix in the parameter d.
  
  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::matstiff(matrix &d,long ipp,long ido)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  double omega;

  if (Mp->nlman == NULL){
    elmatstiff (d,ipp);
  }
  else{  
    switch (Mp->nlman->stmat){
      case initial_stiff:
        elmatstiff (d,ipp);
        break;
      case tangent_stiff:
      case secant_stiff:
      case incr_tangent_stiff:
      case ijth_tangent_stiff:     
        elmatstiff (d,ipp);
        omega=Mm->ip[ipp].eqother[ido+1];
        if (omega > 0.999999)
        omega = 0.999999;
        d(nc-1,nc-1) *= (1-omega);
        // if plasticity has evolved (actual gamma - previous gamma > 0) then reduce shear stiffness
        if (Mm->ip[ipp].other[ido+3] - Mm->ip[ipp].eqother[ido+3] > 0.0){
        for (i=0; i<nc-1; i++)
          d(i,i) *= 0.001;
        }
        break;
      default:
        print_err("unknown type of stifness matrix is required", __FILE__, __LINE__, __func__);
    }
  }
}



/**
  The function computes actual stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::nlstresses (long ipp, long /*im*/, long ido)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  double omega, w, kappa, aft;
  double gamma, dup, epst, tau;
  double u_cor=0.0;
  vector sigma(ASTCKVEC(nc)), dd(ASTCKVEC(nc)), dds(ASTCKVEC(nc-1));
  vector sepst(ASTCKVEC(nc-1)), sepsp(ASTCKVEC(nc-1)), dsepsp(ASTCKVEC(nc-1)), aux;

  // restore actual relative displacements and irreversible relative shear displacements sepsp
  for(i=0; i<nc; i++){
    dd(i) = Mm->ip[ipp].strain[i];
    if (i<nc-1){
      dds(i) = dd(i);
      sepsp(i) = Mm->ip[ipp].eqother[ido+4+i];
    }
  }

  // calculate trial increment of relative displacements in slip directions, sepst
  subv(dds, sepsp, sepst);
  epst = normv(sepst);
    
  // restore state variables
  kappa = Mm->ip[ipp].eqother[ido+0];
  omega = Mm->ip[ipp].eqother[ido+1];
  w     = Mm->ip[ipp].eqother[ido+2]; 
  gamma = Mm->ip[ipp].eqother[ido+3];
  if (coref == yes){
    u_cor  = Mm->ip[ipp].eqother[ido+nc-1+4];
  }

  // actual tensile strength
  aft = Mm->give_actual_ft(ipp);

  if (coref == yes){
    // corrosion effect is taken into account
    // numerical quadrature of the corroded thickness
    double dt = Mp->timecon.actualforwtimeincr();
    double tn = Mp->time;
    double to = Mp->time - dt;
    const long order = 2;
    vector gp(ASTCKVEC(order)), w(ASTCKVEC(order));
    double du_cor = 0.0;
    gauss_points (gp.a, w.a, order);
    for (i=0; i<order; i++){
      // shift Gauss points and their weights due to integration limits \int_{to}^{tn} dtime
      gp[i] = 0.5*(tn+to)+0.5*dt*gp[i];
      w[i] = 0.5*dt*w[i];
      double icor = icor_func.getval(gp[i]);     // current density at the shifted Gauss point
      double alpha = alpha_c.getval(gp[i]);      // corrosion expansion factor at the shifted Gauss point
      double flux = flux_cor_func.getval(gp[i]); // flux contribution at the shifted Gauss point
      double dxcor = icor*0.0315*1e-6/3600.0;    // increment of corroded thickness [m/s]
      // add contribution to the increment of thickness of corroded layer at the shifted Gauss point
      du_cor += w[i]*((alpha-1.0)*dxcor - flux);
    }
    // add integrated increment of the thickness of corroded layer
    u_cor += du_cor;
    // change the normal displacement difference
    dd(nc-1) -= u_cor;
  }
  // first calculate normal stress component that will be used in the yield condition for the shear stress component
  sigma(nc-1) = compute_normal_stress(ipp, dd(nc-1), aft, kappa, omega, w);

  

  // calculate shear stress in the actual shear direction
  dup = 0.0;
  tau = compute_shear_stress(ipp, epst, aft, omega, kappa, sigma(nc-1), dup, gamma);
  // calculate base vector of the actual shear direction -> sepst
  if (fabs(epst) > 0.0) 
    cmulv(1.0/epst, sepst);
  // calculate increments of plastic shear components
  cmulv(dup, sepst, dsepsp);
  addv(sepsp, dsepsp, sepsp);

  // calcualte shear stress components
  makerefv(aux, sigma.a, nc-1);
  cmulv(tau, sepst, aux);

  // save actual stresses
  Mm->storestress(0, ipp, sigma);

  // save new values of the state variables
  Mm->ip[ipp].other[ido+0] = kappa;
  Mm->ip[ipp].other[ido+1] = omega;
  Mm->ip[ipp].other[ido+2] = w; 
  Mm->ip[ipp].other[ido+3] = gamma;
  for(i=0; i<nc-1; i++)
    Mm->ip[ipp].other[ido+4+i] = sepsp(i);
  if (coref == yes){
    Mm->ip[ipp].other[ido+nc-1+4] = u_cor;
  }
}



/**
  The function updates values in the eqother array, containing values from the previous equlibrium state, 
  to the values attained in the actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::updateval (long ipp,long im,long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}


/**
  The function computes actual value of the normal stress component on the interface 
  element according to the simple damage model.

  @param ipp - integration point pointer
  @param dv - actual value of the relative displacements in the normal direction (input)
  @param aft - actual tensile strength in normal direction (input)
  @param kappa - history variable with the maximum attained crack width (input/ouptut)
  @param omega - damage parameter (input/output)
  @param w - crack width (input/output)

  @return The function returns actual normal stress component corresponding to the actual dv.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::compute_normal_stress(long ipp, double dv, double aft, double &kappa, double &omega, double &w)
{
  long i, ni  = sra.give_ni();
  double err = sra.give_err();
  double indam, ift, awf, omegan, tmp, dtmp;
  double sig = 0.0;

  // initial value of the tensile strength is taken from this material model
  ift = ft;

  if (dv > kappa) // loading
  {
    kappa = dv;
    // actual value of initial crack opening
    awf = wf*(aft/ift);
    // initial damage treshold
    indam = aft/kn;

    if (kappa <= indam) // elastic loading
      omegan = 0.0;
    else // cracking => compute damage parameter omega
    {
      // it is better to start from omega=1 in order to avoid of convergention problems
      omegan = 1.0;
      tmp = 0.0;
      for (i=0; i<ni; i++)
      {
        tmp = (1.0-omegan)*kn*kappa - aft*exp(-omegan*kappa/awf);
        if (fabs(tmp) < err*aft)
          break;
        dtmp = -kn*kappa + aft*kappa/awf*exp(-omegan*kappa/awf);
        omegan += -tmp/dtmp;
      }
      if (fabs(tmp) > err*aft)
      {
        print_err("iteration procedure cannot compute damage parameter omega with required error\n"
                  "(dsigma=%le, err=%le)", __FILE__, __LINE__, __func__, tmp, err*aft);
        abort();
      }
    }
    // caution - increment of omega must be positive
    if (omegan > omega){
      omega = omegan;
      fprintf(Out, "\npositive omega increment, eid=%ld, ipp=%ld", Mm->elip[ipp]+1, ipp);
    }

    // actual values of stress and crack width
    sig = (1.0-omega)*kn*dv;
    w   = omega*dv; 
  }
  else{ // unloading/reloading/compression
    if (dv >= 0.0) // unloading/reloading in tension
      sig = (1.0-omega)*kn*dv;
    else // compression => no stiffness reduction
      sig = kn*dv;
  }
  return sig;
}



/**
  The function computes actual value of the shear stress component on the interface 
  element according to the simple plasticity model with Mohr-Coulomb criterion.

  @param ipp - integration point pointer
  @param du - actual value of the relative displacements in the tangential direction (input)
  @param aft - actual tensile strength in normal direction (input)
  @param omega - actual value of damage parameter (input)
  @param kappa - maximum relative displacement in the normal direction ever attained in history (input)
  @param sigma - actual value of the normal stress component (input)
  @param dup - irreversible component of the relative displacements in the tangential direction (input/output)
  @param gamma - cumulative value of consistency parameter (input/output)

  @return The function returns actual shear stress component corresponding to the actual du.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::compute_shear_stress(long ipp, double du, double /*aft*/, double omega, double /*kappa*/, double sigma, double &dup, double &gamma)
{
  long i, ni = sra.give_ni();
  double ac, tau0;
  double err = sra.give_err(); // required relative error of residual
  double tau_tr = 0.0;
  double f, dgamma;

  // actual cohesion according to attained damage parameter
  ac = (1.0-omega)*c;

  // actual yield stress
  tau0 = tanphi*sigma - ac;
  if (tau0 > 0.0)
    tau0 = 0.0;

  // absolute value of residual error
  if (c != 0.0) 
    err *= c; 

  i  = 0;
  do
  {
    // trial shear stress - elastic predictor
    tau_tr = ks*(du-dup);
    // yield function evaluation - Mohr-Coulomb criterion with given normal stress => simple plasticity
    f = fabs(tau_tr) + tau0;
    if (f <= err) // elastic step
      break;
    else // compute plastic corrector
    {
      // increment of consistency parameter
      dgamma = f/ks;
      // increment of irreversible relative displacement
      dup += dgamma*sgn(tau_tr);
      // cumulative consistency parameter
      gamma += dgamma;
      i++;
    }
  } while((f > err) && (i < ni));

  if (f > err)
  {
    print_err("convergence problem of stress return algorithm (f=%le, err=%le)", 
              __FILE__, __LINE__, __func__, f, err);
    abort();
  }
  if(gamma > 0.0)
    fprintf(Out, "\nconsitency parameter gamma > 0, eid=%ld, ipp=%ld", Mm->elip[ipp]+1, ipp);


  return tau_tr;
}



/**
  The function returns actual value of the tensile strength.

  @return The function returns actual value of the tensile strength.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::give_actual_ft()
{
  return ft;
}



/**
  The function returns the value of the limit elastic relative displacement 
  which is dependent on actual value of the normal stiffness coefficient 
  and actual value of tensile strength.

  @param ipp - integration point number in the mechmat ip array.

  @return Function returns actual value of limit elastic relative displacement.
  
  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::epsefunction(long ipp)
{
  double ft = Mm->give_actual_ft(ipp);
  return (ft/kn) ;
}



/**
  This function returns the value of damage parameter obtained from the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the value of damage parameter omega.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::givedamage(long ipp, long ido)
{
  return Mm->ip[ipp].other[ido+1];
}



/**
  This function returns the value of damage parameter obtained from the actual attained state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the value of damage parameter omega.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::give_crackwidth(long ipp, long ido)
{
  return Mm->ip[ipp].other[ido+2];
}



/**
  This function returns the value of cumulative consistency parameter gamma obtained from 
  the last converged state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the value of cumulative consistency parameter gamma.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
double damplastifacemat::give_consparam(long ipp, long ido)
{
  return Mm->ip[ipp].eqother[ido+3];
}



/**
  This function returns the value of irreversible relative displacements 
  obtained from the last converged state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  @param ddp - vector of irreversible displacements (output)
  
  @return Function returns the value of irreversible relative displacements in the argument ddp.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::giveirrstrains(long ipp, long ido, vector &ddp)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  
  for(i=0; i<nc-1; i++){
    ddp(i) = Mm->ip[ipp].eqother[ido+4+i];
  }
  ddp(nc-1) = 0.0;
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.


  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
void damplastifacemat::changeparam(atsel &atm,vector &val)
{
  long i;
  
  for (i=0; i<atm.num; i++)
  {
    switch (atm.atrib[i])
    {
      case 0:
        ks=val[i];
        break;
      case 1:
        kn=val[i];
        break;
      case 2:
        ft=val[i];
        break;
      case 3:
        wf=val[i];
        break;
      case 4:
        phideg=val[i];
        break;
      default:
        print_err("wrong index of changed atribute",__FILE__,__LINE__,__func__);
        abort();
    }
  }
}
  

