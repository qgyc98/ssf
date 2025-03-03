#include "plastifacemat.h"
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

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
plastifacemat::plastifacemat()
{
  ks = kn = 0.0;
  sig0t = sig0c = phideg = tanphi = c = 0.0;
  rdn_lt = rdn_lc = 0.0;
}



/**
  Destructor deallocates used heap memory.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
plastifacemat::~plastifacemat()
{
}



/**
  The function reads the model parameter sform the opened text file
  given by the pointer in.

  @param in - pointer to the opened XFILE structure
  
  @return - The function does not return anything but stores new values in the model parameters.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::read(XFILE *in)
{
  xfscanf(in, "%k%le %k%le", "ks", &ks, "kn", &kn);
  xfscanf(in, "%k%le %k%le %k%le %k%le", 
          "sig_0t", &sig0t, "sig_0c", &sig0c, "phi_deg", &phideg, "cohesion", &c);
  if (sig0t < 0.0){
    print_err("tensile yield stress sig_0t=%le is negative", __FILE__, __LINE__, __func__, sig0t);
    abort();
  }
  if (sig0c > 0.0){
    print_err("compressive yield stress sig_0c=%le is positive", __FILE__, __LINE__, __func__, sig0c);
    abort();
  }
  sra.read(in);
  tanphi = tan(phideg*M_PI/180);
  rdn_lt = sig0t/kn;
  rdn_lc = sig0c/kn;
}



/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened output text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::print(FILE *out)
{
  fprintf(out, "%le %le ", ks, kn);
  fprintf(out, "%le %le %le %le ", sig0t, sig0c, phideg, c);
  sra.print(out);
}



/**
  The function computes elastic material stiffnes %matrix from actual Young modulus.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

  @return The function returns computed elastic stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::elmatstiff(matrix &d, long ipp)
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
  
  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::matstiff(matrix &d,long ipp,long /*ido*/)
{
  if (Mp->nlman == NULL){
    elmatstiff (d,ipp);
  }
  else{  
    switch (Mp->nlman->stmat){
      case initial_stiff:
      case tangent_stiff:
      case secant_stiff:
      case incr_tangent_stiff:
      case ijth_tangent_stiff:     
        elmatstiff (d,ipp);
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

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::nlstresses (long ipp, long /*im*/, long ido)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  double rdn, rdnp, sign;
  double trdst, trdsp, gammas, tau;
  vector sigma(ASTCKVEC(nc));
  vector rd(ASTCKVEC(nc)), rds(ASTCKVEC(nc-1));
  vector rdst(ASTCKVEC(nc-1)), rdsp(ASTCKVEC(nc-1)), drdsp(ASTCKVEC(nc-1)), aux;

  // restore actual relative displacements and irreversible relative shear displacements rdsp
  for(i=0; i<nc; i++){
    rd(i) = Mm->ip[ipp].strain[i];
    if (i<nc-1){
      rds(i) = rd(i);
      rdsp(i) = Mm->ip[ipp].eqother[ido+2+i];
    }
  }
  // realtive displacement in the normal direction
  rdn = rd(nc-1);

  // calculate trial increment of relative displacements in slip directions, rdst
  subv(rds, rdsp, rdst);
  trdst = normv(rdst);
    
  // restore state variables
  rdnp   = Mm->ip[ipp].eqother[ido+0];
  gammas = Mm->ip[ipp].eqother[ido+1];

  // first calculate normal stress component that will be used in the yield condition for the shear stress component
  sign = compute_normal_stress(ipp, rdn, rdnp);

  // calculate shear stress in the actual shear direction
  trdsp = 0.0;
  tau = compute_shear_stress(ipp, trdst, sign, trdsp, gammas);
  // calculate base vector of the actual shear direction -> rdst
  if (fabs(trdst) > 0.0) 
    cmulv(1.0/trdst, rdst);
  // calculate increments of plastic shear components
  cmulv(trdsp, rdst, drdsp);
  addv(rdsp, drdsp, rdsp);

  // calcualte shear stress components
  makerefv(aux, sigma.a, nc-1);
  cmulv(tau, rdst, aux);

  // store normal stress
  sigma(nc-1) = sign;
  // save actual stresses
  Mm->storestress(0, ipp, sigma);

  // save new values of the state variables
  Mm->ip[ipp].other[ido+0] = rdnp;
  Mm->ip[ipp].other[ido+1] = gammas;
  for(i=0; i<nc-1; i++)
    Mm->ip[ipp].other[ido+2+i] = rdsp(i);
}



/**
  The function updates values in the eqother array, containing values from the previous equlibrium state, 
  to the values attained in the actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function does not return anything.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::updateval (long ipp,long im,long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}


/**
  The function computes actual value of the normal stress component on the interface 
  element according to the simple damage model.

  @param[in]     ipp - integration point pointer
  @param[in]     rdn - tottal attained relative displacement in the normal direction
  @param[in,out] rdnp - actual tensile strength in the normal direction

  @return The function returns actual normal stress component corresponding to the actual relative displacement rdn.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
double plastifacemat::compute_normal_stress(long ipp, double rdn, double &rdnp)
{
  double drdnp = 0.0;
  double rdnt = 0.0;
  double sig  = 0.0;

  // compute trial normal strain increment
  rdnt = fabs(rdn - rdnp);
  if (rdnt > rdn_lt){ // limit tensile relative displacement is exceeded
    // compute plastic strain increment in the normal direction
    drdnp = sgn(rdn-rdnp)*(rdnt - rdn_lt);
  }
  if (rdnt < rdn_lc){ // limit compressive relative displacement is exceeded
    // compute plastic strain increment in the normal direction
    drdnp = sgn(rdn-rdnp)*(rdnt - rdn_lc);
  }

  // actualize plastic strains and compute stresses
  rdnp += drdnp;
  sig   = kn*(rdn-rdnp);

  return sig;
}



/**
  The function computes actual value of the shear stress component on the interface 
  element according to the simple plasticity model with Mohr-Coulomb criterion.

  @param[in]     ipp - integration point pointer
  @param[in]     rds - actual value of the relative displacements in the slip direction (input)
  @param[in]     sigma - actual value of the normal stress component (input)
  @param[in,out] rdsp - irreversible component of the relative displacements in the slip direction (input/output)
  @param[in,out] gamma - cumulative value of consistency parameter (input/output)

  @return The function returns actual shear stress component corresponding to the actual rds.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
double plastifacemat::compute_shear_stress(long ipp, double rds, double sigma, double &rdsp, double &gamma)
{
  long i, ni = sra.give_ni();
  double tau0;
  double err = sra.give_err(); // required relative error of residual
  double tau_tr = 0.0;
  double f, dgamma;

  // actual yield stress
  tau0 = tanphi*sigma - c;
  if (tau0 > 0.0)
    tau0 = 0.0;

  // absolute value of residual error
  if (c != 0.0)
    err *= c; 

  i  = 0;
  do
  {
    // trial shear stress - elastic predictor
    tau_tr = ks*(rds-rdsp);
    // yield function evaluation - Mohr-Coulomb criterion with given normal stress => simple plasticity
    f = fabs(tau_tr) + tau0;
    if (f <= err) // elastic step
      break;
    else // compute plastic corrector
    {
      // increment of consistency parameter
      dgamma = f/ks;
      // increment of irreversible relative displacement
      rdsp += dgamma*sgn(tau_tr);
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
  This function returns the value of cumulative consistency parameter gamma obtained from 
  the last converged state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  
  @return Function returns the value of cumulative consistency parameter gamma.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
double plastifacemat::give_consparam(long ipp, long ido)
{
  return Mm->ip[ipp].eqother[ido+3];
}



/**
  This function returns the value of irreversible relative displacements 
  obtained from the last converged state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array
  @param rdp - vector of irreversible relative displacements (output)
  
  @return Function returns the value of irreversible relative displacements in the argument rdp.

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::giveirrstrains(long ipp, long ido, vector &rdp)
{
  long i, nc = Mm->ip[ipp].ncompstr;
  
  for(i=0; i<nc-1; i++){
    rdp(i) = Mm->ip[ipp].eqother[ido+4+i];
  }
  rdp(nc-1) = 0.0;
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.


  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
void plastifacemat::changeparam(atsel &atm,vector &val)
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
        sig0t=val[i];
        break;
      case 3:
        sig0c=val[i];
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
  

