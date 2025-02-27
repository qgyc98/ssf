#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "drprag2.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"



#define nijac 20
#define limit 1.0e-8




/**
  This constructor inializes attributes to zero values.
*/
drprag2::drprag2 (void)
{
  phi=0.0;  c=0.0;  psi=0.0;
  theta = 0.0; clim = 0.0;
}



/**
  This destructor is only for the formal purposes.
*/
drprag2::~drprag2 (void)
{

}



/**
  This function reads material parameters from the opened text file given
  by the parameter in. Then it computes material constants alpha, alpha1 and
  beta

  @param in - pointer to the opned text file

  25.3.2002  by Tomas Koudelka
*/
void drprag2::read (XFILE *in)
{
  xfscanf (in,"%k%lf%k%lf%k%lf%k%lf%k%lf","phi",&phi,"coh",&c,"psi",&psi,"theta",&theta,"clim",&clim);
  sra.read (in);

  // outter coincidence with Mohr-Coulomb - compression cone
  // alpha=6.0*sin(phi)/(sqrt(3.0)*(3.0-sin(phi)));
  // inner coincidence with Mohr-Coulomb - extension cone
  // alpha=6.0*sin(phi)/(sqrt(3.0)*(3.0+sin(phi)));
  // identical collapse as Mohr-Coulomb in plain strain conditions
  alpha=3.0*tan(phi)/sqrt(9.0 + 12.0*tan(phi)*tan(phi));

  // alpha1=6.0*sin(psi)/(sqrt(3.0)*(3.0-sin(psi)));
  // alpha1=6.0*sin(psi)/(sqrt(3.0)*(3.0+sin(psi)));
  alpha1=3.0*tan(psi)/sqrt(9.0 + 12.0*tan(psi)*tan(psi));

  // beta=6.0*cos(phi)/(sqrt(3.0)*(3.0-sin(phi)));
  // beta=6.0*cos(phi)/(sqrt(3.0)*(3.0+sin(phi)));
  beta=3.0/sqrt(9.0 + 12.0*tan(phi)*tan(phi));
}



/**
  This function prints material parameters to the opened text file given
  by the parameter in.

  @param out - pointer to the opned text file

  25.3.2015  by Tomas Koudelka
*/
void drprag2::print (FILE *out)
{
  fprintf (out,"%le %le %le %le %le ", phi, c, psi, theta, clim);
  sra.print (out);
}



void drprag2::giveIota(strastrestate ssst, vector &iota)
{
  fillv(0.0, iota);
  switch (ssst)
  {
    case bar:
      iota[0]=1.0;
      break;
    case planestrain:
      iota[0]=iota[1]=iota[3]=1.0;
      break;
    case axisymm:
    case spacestress:
      iota[0]=iota[1]=iota[2]=1.0;
      break;
    default:
      print_err("stress/strain state %d is not supported in this model", __FILE__, __LINE__, __func__, ssst);
  }
}



void drprag2::giveIdev(strastrestate ssst, matrix &idev)
{
  fillm(0.0, idev);
  switch (ssst)
  {
    case bar:
      idev[0][0]=1.0;
      break;
    case planestrain:
      idev[0][0] =  2.0/3.0; idev[0][1] = -1.0/3.0; idev[0][3] = -1.0/3.0;
      idev[1][0] = -1.0/3.0; idev[1][1] =  2.0/3.0; idev[1][3] = -1.0/3.0;
      idev[2][2] = 0.5;
      idev[3][0] = -1.0/3.0; idev[3][1] = -1.0/3.0; idev[3][3] =  2.0/3.0;
      break;
    case axisymm:
      idev[0][0] =  2.0/3.0; idev[0][1] = -1.0/3.0; idev[0][2] = -1.0/3.0;
      idev[1][0] = -1.0/3.0; idev[1][1] =  2.0/3.0; idev[1][2] = -1.0/3.0;
      idev[2][0] = -1.0/3.0; idev[2][1] = -1.0/3.0; idev[2][2] =  2.0/3.0;
      idev[3][3] = 0.5;
      break;
    case spacestress:
      idev[0][0] =  2.0/3.0; idev[0][1] = -1.0/3.0; idev[0][2] = -1.0/3.0;
      idev[1][0] = -1.0/3.0; idev[1][1] =  2.0/3.0; idev[1][2] = -1.0/3.0;
      idev[2][0] = -1.0/3.0; idev[2][1] = -1.0/3.0; idev[2][2] =  2.0/3.0;
      idev[3][3] = 0.5;
      idev[4][4] = 0.5;
      idev[5][5] = 0.5;
      break;
    default:
      print_err("stress/strain state %d is not supported in this model", __FILE__, __LINE__, __func__, ssst);
  }
}



void drprag2::getImportantValues(long ipp, long ido, double &K_con, double &G_con, double &rho_tr, double &p_tr)
{
  //double K_con, G_con, rho_tr, p_tr;
  double E, nu;
  //double over_epsilon_ptr;
  double scalar_prod1, scalar_prod2;
  //double H_value, Hder_value;
  //double part_q1, part_q2, part_q3;
  long n=Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epsilon(n), epsilon_p(n), epsilon_etr(n), aux_epsilon_etr1(n); 
  vector iota(n);
  matrix Idev(n,n);
  vector aux_epsilon_etr2(n);
  
  E = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  K_con = E/(3*(1-2*nu));
  G_con = E/(2*(1+nu));
  
  giveIdev(ssst, Idev);
  
  //compute trial epsilon^e
  giveirrstrains (ipp, ido, epsilon_p);
  Mm->givestrain(0, ipp, epsilon);
  subv(epsilon, epsilon_p, epsilon_etr);
  
  //compute rho^{tr}
  mxv(Idev, epsilon_etr, aux_epsilon_etr1);
  scprd(epsilon_etr, aux_epsilon_etr1, scalar_prod1);
  if (scalar_prod1 < 0.0)
    scalar_prod1 = 0.0;
  rho_tr = 2.0*G_con*sqrt(scalar_prod1);
  
  //compute p^{tr}
  giveIota(ssst, iota);
  scprd(epsilon_etr, iota, scalar_prod2);
  p_tr = K_con*scalar_prod2;
  
}


void drprag2::getImportantValues(long ipp, long ido, double &K_con, double &G_con, double &rho_tr, double &p_tr, 
                                 vector &iota, vector &epsilon_etr, matrix &Idev)
{
  //double K_con, G_con, rho_tr, p_tr;
  double E, nu;
  //double over_epsilon_ptr;
  double scalar_prod1, scalar_prod2;
  //double H_value, Hder_value;
  //double part_q1, part_q2, part_q3;
  long n=Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epsilon(n), epsilon_p(n), aux_epsilon_etr1(n); //epsilon_etr,
  
  E = Mm->give_initial_ym(ipp);
  nu = Mm->give_initial_nu(ipp);
  K_con = E/(3*(1-2*nu));
  G_con = E/(2*(1+nu));
  
  //define matrix Idev
  giveIdev(ssst, Idev);
  
  giveirrstrains (ipp, ido, epsilon_p);
  Mm->givestrain(0, ipp, epsilon);
  subv(epsilon, epsilon_p, epsilon_etr);
  
  //compute rho^{tr}
  mxv(Idev, epsilon_etr, aux_epsilon_etr1);
  scprd(epsilon_etr, aux_epsilon_etr1, scalar_prod1);
  if (scalar_prod1 < 0.0)
    scalar_prod1 = 0.0;
  rho_tr = 2.0*G_con*sqrt(scalar_prod1);
    
  
  //compute p^{tr}
  giveIota(ssst, iota);
  scprd(epsilon_etr, iota, scalar_prod2);
  p_tr = K_con*scalar_prod2;
}


/**
 * This function get the value of function q which is define in the paper prepare
 * by Standa Sysala
 * 
 * @param gamma      - single value, on this value depends the function q
 * @param ipp        - integration point pointer
 * @param ido        - index of internal variables for given material in the ipp other array
 * @return q_value   - return value for this function
 * 
 * 22.8.2014 Martin Cermak
 */

double drprag2::getQvalue(double gamma, long ipp, long ido)
{
  double K_con, G_con, rho_tr, p_tr, q_value;
  //double E, nu;
  //double over_epsilon_ptr;
  //double scalar_prod1, scalar_prod2;
  double H_value, Hder_value;
  double part_q1, part_q2, part_q3;
  getImportantValues(ipp, ido, K_con, G_con, rho_tr, p_tr);
  
  //get value of \overline{\varepsilon}^{p,tr} (hardening value)
  //over_epsilon_ptr=Mm->ip[ipp].eqother[ido+ncomp];
  
  getValueOfHardening(ipp, ido, H_value, Hder_value);
  
  part_q1 = sqrt(0.5)*max2((rho_tr-gamma*G_con*sqrt(2.0)),0);
  part_q2 = alpha*(p_tr - gamma*K_con*alpha1);
  part_q3 = beta*(c+H_value);       //suppose that beta is ksi and c is c_0 in Standa paper
  
  q_value = part_q1 + part_q2 - part_q3;
  
  
  return q_value;
  
}


/**
 * This function get the value of function H which describe in the paper prepare
 * by Standa Sysala 
 * @param ipp          - integration point pointer
 * @param ido          - index of internal variables for given material in the ipp other array
 * @param H_value      - return value for hardening
 * @param Hder_value   - return value for derivation of function H
 * @return 
 * 
 * 22.8.2014 Martin Cermak
 */

void drprag2::getValueOfHardening(long ipp, long ido, double &H_value, double &Hder_value)
{
  double eps_eqpl = 0.0;
  long ncomp=Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epsp(ncomp);
  matrix epst(ncomp,ncomp);
  double tgt = tan(theta);
  double cq;
  
  giveirrstrains (ipp, ido, epsp);
  // converts vector form of the epsp to tensor form epst 
  vector_tensor (epsp, epst, ssst, strain);
  // calculates equivalent plastic strain sqrt(2/3)*|eps_p|

  // linear dependence of cohesion on hardening parameter with cutoff on the clim
  cq = c + tgt*eps_eqpl;
  if ((cq > clim) && (c < clim))
  {
    cq = clim;
    H_value = cq - c;
    Hder_value = 0;
  }
  if ((cq < clim) && (c > clim))
  {
    cq = clim;
    H_value = cq - c;
    Hder_value = 0;
  }
  
  H_value = cq - c;
  Hder_value = tgt * sqrt(1.0+2.0/9.0*alpha1*alpha1);
}

/**
 * This function compute sufficient value of \delta \lambda by the Newton method 
 * for one Gauss point for the second case in the function tangentstiffness 
 * matrix where we come back to the appex
 * 
 * @param ipp - integration point pointer
 * @param ido - index of internal variables for given material in the ipp other array
 * @return sufficient value of \delta\lambda
 * 
 * 26.8.2014 Martin Cermak
 */

double drprag2::newton_one_element_case2(long ipp, long ido)
{
  double K_con, G_con, rho_tr, p_tr, q_value;
  double H_value, Hder_value;
  double delta_lambda, delta_lambda_next, tau;
  double part_q2, part_q3, qDer_value;
  int counter_iter;
  
  getImportantValues(ipp, ido, K_con, G_con, rho_tr, p_tr);
  
  delta_lambda = rho_tr/(G_con*sqrt(2));
  delta_lambda_next = 0.0;
  tau = sra.err;
  counter_iter = 0;
  
  while(1) {
        
    getValueOfHardening(ipp, ido, H_value, Hder_value);
    part_q2 = alpha*(p_tr - delta_lambda*K_con*alpha1);
    part_q3 = beta*(c+H_value);       //suppose that beta is ksi and c is c_0 in Standa paper
    q_value = part_q2 - part_q3;
    qDer_value = -alpha*alpha1*K_con - beta*Hder_value;

    delta_lambda_next = delta_lambda - q_value/qDer_value;

    if ((delta_lambda_next-delta_lambda)/(delta_lambda_next+delta_lambda)<=tau || counter_iter >= sra.ni) {
      break;
    }
    else {
      delta_lambda = delta_lambda_next;
    }
    counter_iter= counter_iter+1;
  }
  
  return delta_lambda_next;
}

/**
 * This function compute sufficient value of \delta \lambda by the Newton method 
 * for one Gauss point for the third case in the function tangentstiffness 
 * matrix where we come back to the smooth surface
 * 
 * @param ipp - integration point pointer
 * @param ido - index of internal variables for given material in the ipp other array
 * @return sufficient value of \delta\lambda
 * 
 * 26.8.2014 Martin Cermak
 */

double drprag2::newton_one_element_case3(long ipp, long ido)
{
  double K_con, G_con, rho_tr, p_tr, q_value;
  double H_value, Hder_value;
  double delta_lambda, delta_lambda_next, tau;
  double part_q1, part_q2, part_q3, qDer_value;
  int counter_iter;
  
  getImportantValues(ipp, ido, K_con, G_con, rho_tr, p_tr);
  
  delta_lambda = 0.0;
  delta_lambda_next = 0.0;
  tau = sra.err;
  counter_iter = 0;
  
  while(1) {
        
    getValueOfHardening(ipp, ido, H_value, Hder_value);
    part_q1 = sqrt(0.5)*max2((rho_tr-delta_lambda*G_con*sqrt(2.0)),0);
    part_q2 = alpha*(p_tr - delta_lambda*K_con*alpha1);
    part_q3 = beta*(c+H_value);       //suppose that beta is ksi and c is c_0 in Standa paper
    q_value = part_q1 + part_q2 - part_q3;
    qDer_value = -G_con - alpha*alpha1*K_con - beta*Hder_value;

    delta_lambda_next = delta_lambda - q_value/qDer_value;

    if ((delta_lambda_next-delta_lambda)/(delta_lambda_next+delta_lambda)<=tau || counter_iter >= sra.ni) {
      break;
    }
    else {
      delta_lambda = delta_lambda_next;
    }
    counter_iter= counter_iter+1;
  }
  
  return delta_lambda_next;
}


/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns calculated stiffness %matrix in the parameter d.

  4.4.2002 by Tomas Koudelka
*/
void drprag2::matstiff (matrix &d, long ipp,long ido)
{
  if (Mp->nlman->stmat==initial_stiff)
  {
    //  initial elastic matrix
    Mm->elmatstiff(d, ipp, ido);
  }
  if (Mp->nlman->stmat==tangent_stiff)
  {
    //  tangent stiffness matrix
    tangentstiff(d, ipp, ido);
  }
}



/**
  This function computes elastic-plastic material stiffness %matrix.

  @param td - allocated matrix structure for material tangent stiffness %matrix 
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns calculated tangent stiffness %matrix in the parameter td.

  26.8.2014 by Martin Cermak
*/
void drprag2::tangentstiff (matrix &tde,long ipp,long ido)
{
  long n=Mm->ip[ipp].ncompstr;
  long i,j;

  matrix Idev(n,n);
  double K_con, G_con, rho_tr, p_tr, q_value0, q_value1, delta_lambda;
  double H_value, H1, gamma;
  vector iota(n); 
  vector epsilon_etr(n);
  double aux_con1, aux_con2;
  vector s_tr(n), n_tr(n), aux_s_tr(n);
  matrix aux_td2(n,n), aux_td3(n,n), aux_Ntr(n,n), td(n,n);
  vector vec_td31(n), vec_td32(n);
  matrix d(n,n);

  Mm->elmatstiff (d,ipp);
  getImportantValues(ipp, ido, K_con, G_con, rho_tr, p_tr, iota, epsilon_etr, Idev);
  
  gamma = 0.0;
  q_value0 = getQvalue(gamma, ipp, ido);
  gamma = rho_tr/(G_con*sqrt(2));
  q_value1 = getQvalue(gamma, ipp, ido);
  
  //elastic case q(0) <= 0
  if (q_value0 < 1e-10) {
    for(i=0; i<d.m; i++)
    {
      for(j=0; j<d.n; j++)
        td(i,j)=d(i,j);
    } 
  }
  //case 2: return to the apex
  else if (q_value1 >= 1e-10) {
    //satisfy equation 0 = q(\Delta\lambda) = \eta(p^{tr}-\Delta\lambda K\overline{\eta})
    // - \ksi(c_0 + H(\overline{\varepsilon}^{p,tr} + \Delta\lambda \ksi))
    delta_lambda = newton_one_element_case2(ipp, ido);
    //get value of H1 which is value of derivation of function H
    getValueOfHardening(ipp, ido, H_value, H1);
    
    //Tder_lok = (ksi^2*K_con*H_1)/(K_con*eta*eta_over+ksi^2*H_1)*iota*iota';
    aux_con1 = (beta*beta*K_con*H1)/(K_con*alpha*alpha1 + beta*beta*H1);       // suppose beta is ksi and alpha is eta
    vxv(iota, iota, td);
    cmulm(aux_con1, td);
  }
  else if ((q_value0 > 1e-10) && (q_value1 < 1e-10)) {
    //satisfy equation 0 = q(\Delta\lambda) = \sqrt{0.5}(rho^{tr}-\Delta\lambda*G_con*\sqrt{2}) + \eta(p^{tr}-\Delta\lambda K\overline{\eta})
    // - \ksi(c_0 + H(\overline{\varepsilon}^{p,tr} + \Delta\lambda \ksi))
    delta_lambda = newton_one_element_case3(ipp, ido);
    //get value of H1 which is value of derivation of function H
    getValueOfHardening(ipp, ido, H_value, H1);
    
    mxv(Idev, epsilon_etr, aux_s_tr);
    cmulv((2*G_con), aux_s_tr, s_tr);
    cmulv((1.0/rho_tr), s_tr, n_tr);
    
    //compute elmentary part of matrix td
    aux_con1 = delta_lambda * (G_con*sqrt(2)/rho_tr);
    vxv(n_tr, n_tr, aux_Ntr);
    subm(Idev, aux_Ntr, aux_td2);
    cmulm(aux_con1, aux_td2);
    
    aux_con2 = 1.0/(G_con + K_con*alpha*alpha1 + beta*beta*H1);   //suppose that alpha is eta and beta is ksi    
    fillv(0.0, vec_td31);
    addmultv (vec_td31, n_tr, (G_con*sqrt(2)));
    addmultv (vec_td31, iota, (K_con*alpha1));
    fillv(0.0, vec_td32);
    addmultv (vec_td32, n_tr, (G_con*sqrt(2)));
    addmultv (vec_td32, iota, (K_con*alpha));
    cmulv(aux_con2, vec_td32);
    vxv(vec_td31, vec_td32, aux_td3);
    
    //copym (d,td);
    subm(d, aux_td2, td);
    subm(td, aux_td3, td);
    
  }
  else {
    
  }
  for(i=0; i<tde.m; i++)
  {
    for(j=0; j<tde.n; j++)
      tde(i,j)=td(i,j);
  } 
  
}



/**
  This function computes stresses at given integration point ipp,
  depending on the reached strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the stress array of the given integration point.

  27.8.2014 by Martin Cermak
*/
void drprag2::nlstresses (long ipp, long /*im*/, long ido)
{
  long i,n=Mm->ip[ipp].ncompstr;  
  vector epsp(n),q(1),q_next(1);
  
  matrix Idev(n,n);
  matrix De(n,n);
  double K_con, G_con, rho_tr, p_tr, q_value0, q_value1, delta_lambda=0.0;
  double gamma, aux_con1, aux;
  vector iota(n), epsilon(n), epsilon_etr(n), epsilon_p_previous(n);
  vector epsilon_p_next(n), aux_vec(n);
  vector s_tr(n), n_tr(n), aux_s_tr(n), sigma_tr(n), T_lok(n);
  
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  
  getImportantValues(ipp, ido, K_con, G_con, rho_tr, p_tr, iota, epsilon_etr, Idev);
  
  gamma = 0.0;
  q_value0 = getQvalue(gamma, ipp, ido);
  gamma = rho_tr/(G_con*sqrt(2));
  q_value1 = getQvalue(gamma, ipp, ido);
  
  Mm->elmatstiff(De,ipp);
  mxv(De, epsilon_etr, sigma_tr);
  
  //get strain and plastic strain
  Mm->givestrain(0, ipp, epsilon);
  giveirrstrains (ipp, ido, epsilon_p_previous);
      
  //elastic case q(0) <= 0
  if (q_value0 < 1e-10) {
    copyv (sigma_tr,T_lok);
    copyv (epsilon_p_previous, epsilon_p_next);
    copyv (q, q_next);
    delta_lambda = 0.0;
  }
  //case 2: return to the apex
  else if (q_value1 >= 1e-10) {
    //satisfy equation 0 = q(\Delta\lambda) = \eta(p^{tr}-\Delta\lambda K\overline{\eta})
    // - \ksi(c_0 + H(\overline{\varepsilon}^{p,tr} + \Delta\lambda \ksi))
    delta_lambda = newton_one_element_case2(ipp, ido);        
    
    //stress
    aux_con1 = (p_tr - delta_lambda*K_con*alpha1);       // suppose  alpha1 is \overline{eta}    
    cmulv(aux_con1, iota, T_lok);
    
    //plastic strain
    aux_con1 = (p_tr - delta_lambda*K_con*alpha1)/(3.0*K_con);
    fillv(0.0, aux_vec);
    addmultv (aux_vec, iota, aux_con1);
    subv(epsilon, aux_vec, epsilon_p_next);
    
    //overline{\varepsilon}^p
    q_next[0] = q[0] + delta_lambda*beta;       //suppose that beta is ksi
    
  }
  else if ((q_value0 > 1e-10) && (q_value1 < 1e-10)) {
    //satisfy equation 0 = q(\Delta\lambda) = \sqrt{0.5}(rho^{tr}-\Delta\lambda*G_con*\sqrt{2}) + \eta(p^{tr}-\Delta\lambda K\overline{\eta})
    // - \ksi(c_0 + H(\overline{\varepsilon}^{p,tr} + \Delta\lambda \ksi))
    delta_lambda = newton_one_element_case3(ipp, ido);
    
    mxv(Idev, epsilon_etr, aux_s_tr);
    cmulv((2*G_con), aux_s_tr, s_tr);
    cmulv((1.0/rho_tr), s_tr, n_tr);
    
    //stress
    fillv(0.0, aux_vec);
    addmultv (aux_vec, n_tr, (delta_lambda*G_con*sqrt(2)));
    addmultv (aux_vec, iota, (delta_lambda*K_con*alpha1));
    subv(sigma_tr, aux_vec, T_lok);
    
    //plastic strain
    fillv(0.0, aux_vec);
    addmultv (aux_vec, n_tr, (delta_lambda*sqrt(0.5)));
    addmultv (aux_vec, iota, (delta_lambda*alpha1/3.0));
    addv(epsilon_p_previous, aux_vec, epsilon_p_next);
    
    //overline{\varepsilon}^p
    q_next[0] = q[0] + delta_lambda*beta;       //suppose that beta is ksi
    
  }
  else {
    
  }
  
  Mm->storestress(0, ipp, T_lok);

  //  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsilon_p_next[i];
  }
  Mm->ip[ipp].other[ido+n]+=delta_lambda;
  Mm->ip[ipp].other[ido+n+1]=q_next[0];
  aux = 0.0;
  for (i=0; i<n; i++)
  {
    aux += (epsilon_p_next(i) - Mm->ip[ipp].eqother[ido+i])*(epsilon_p_next(i) - Mm->ip[ipp].eqother[ido+i]);
  }
  Mm->ip[ipp].other[ido+n+2]=sqrt(aux);
  
}





/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the eqother array of internal variables in the given integration point.

  4.4.2002 by Tomas Koudelka
*/
void drprag2::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  @return The function returns the %vector of irreversible strains in the parameter epsp

  4.4.2002 by Tomas Koudelka
*/
void drprag2::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  This function extracts consistency parametr gamma for the attained equilibrium state
  from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of consistency parameter.

  4.4.2002 by Tomas Koudelka
*/
double drprag2::give_consparam (long ipp, long ido)
{
  long ncompstr;
  double gamma;

  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];

  return gamma;
}



/**
  This function changes required material parameters and updates
  material constants alpha, alpha1 and beta. It is used for stochastic 
  caluclations.

  @param atm - structure with description of selected parameters
  @param val - %vector of new values for the required material parameters.

  @retval The function returns value of consistency parameter.

  4.4.2002 by Tomas Koudelka
*/
void drprag2::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      phi=val[i];
      break;
    }
    case 1:{
      c=val[i];
      break;
    }
    case 2:{
      psi=val[i];
      break;
    }
    case 3:{
      theta=val[i];
      break;
    }
    case 4:{
      clim=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  alpha=6.0*sin(phi)/(3.0-sin(phi));
  alpha1=6.0*sin(psi)/(3.0-sin(psi));

  beta=6.0*cos(phi)/(3.0-sin(phi));
}
