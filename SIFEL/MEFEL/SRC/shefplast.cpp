#include "shefplast.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include <float.h>

#define nijac 20
#define limit 1.0e-8



/**
   This constructor inializes attributes to zero values.
*/
shefplast::shefplast (void)
{
  rc = 0.0;
  rt = 0.0;
  gamma = 0.99;
  alpha = 0.5;
  p = 0.0;
  a = -1.0;
  k0 = 0.1;
  ah = bh = ch = 0.0;
  kh0 = 0.001;
  reallocv(6, dev); // it must be allocated on the heap because it is passed across all methods
  khat = rhoc = r = m = xi = ft = 0.0;
  b0 = b1 = c0 = c1 = c2 = c3 = c4 = c5 = theta = d0 = d1 = d2 = 0.0;
}



/**
   This destructor is only for the formal purposes.
*/
shefplast::~shefplast (void)
{

}



/**
   This function reads material parameters from the opened text file given
   by the parameter in.

   @param in - pointer to the opned text file

*/
void shefplast::read (XFILE *in)
{
  xfscanf (in, "%lf %lf %lf %lf %lf %lf", &rc, &rt, &p, &ah, &bh, &ch);
  sra.read (in);
  
  ft = rt/rc;
  m = 3.0*(1.0-pow(ft,2.0/gamma));
  m /= ft+2.0*pow(ft,1.0/gamma);
}



void shefplast::compute_khat (vector &q)
{
  double k;

  k = k0 + (1.0 - k0)*sqrt(q[0]*(2.0-q[0]));
  khat = pow(k,2.0*p)*(1.0 - xi*xi*(1.0-k)*(1.0-k)/(a*a));
}



void shefplast::compute_rhoc ()
{
  double tmp;

  tmp = -m+sqrt(m*m-12.0*sqrt(3.0)*m*xi+36.0);
  rhoc = pow(1.0/6.0, gamma)*sqrt(2.0/3.0)*pow(tmp,gamma);
}
 


void shefplast::compute_r ()
{
  double tmp, rhoe;
  double pi=4.0*atan(1.0) ;

  tmp = -m+sqrt(m*m-3.0*sqrt(3.0)*m*xi+9.0);
  rhoe = pow(1.0/3.0, gamma)*sqrt(2.0/3.0)*pow(tmp, gamma);
  b0 = rhoe/rhoc;
  if(b0<=0.5) {
    b0 = 0.5 + 1.0e-6 ;
  }
  else if(b0>1.0) {
    b0 = 1.0 ;
  }
  b1 = sqrt(3.0)*(1.0-alpha)*b0/(1.0+b0)+2.0*alpha*b0/sqrt(3.0);
  c0 = (2.0-sqrt(3.0)*b1)*(2.0*b0-sqrt(3.0)*b1) / ((b1*(1.0+b0)-sqrt(3.0)*b0)*(b1*(1.0+b0)-sqrt(3.0)*b0));
  c1 = 3.0-c0*(1.0+b0)*(1.0+b0);
  c2 = 1.0+3.0*c0*(1.0-b0)*(1.0-b0);
  c3 = 2.0*c0*sqrt(3.0)*(1.0-b0*b0);
  c4 = (1.0+b0)*(1.0-b0*c0);
  c5 = (1.0-b0)*(1.0-3.0*b0*c0);
  if((-sqrt(3.0)/2.0*3.0*j3s/pow(j2s,1.5)) > 1.0) {
    theta = pi/6.0 ;
  }
  else if((-sqrt(3.0)/2.0*3.0*j3s/pow(j2s,1.5)) < -1.0) {
    theta = -pi/6.0 ;
  }
  else
    theta = 1.0/3.0*asin(-sqrt(3.0)/2.0*3.0*j3s/pow(j2s,1.5));
  d0 = c1*cos(theta)*cos(theta) - c2*sin(theta)*sin(theta) + c3*sin(theta)*cos(theta);
  d1 = 2.0*(c4*sqrt(3.0)*cos(theta)-c5*sin(theta));
  d2 = b0*(4.0-3.0*b0*c0);
  r = 2.0*d0/(d1-sqrt(d1*d1-4.0*d0*d2));
}



void shefplast::compute_classvar(vector &sig, vector &q)
{
  vector devv(ASTCKVEC(6));
  
  deviator(sig,devv);
  i1s = first_invar(sig);
  j2s = j2_stress_invar (sig);
  j3s = third_stress_invar(devv);  
  xi = i1s/(rc*sqrt(3.0));
  compute_khat (q);
  compute_rhoc ();
  compute_r ();

}



/**
   This function computes the value of yield functions.

   @param sig - stress tensor
   @param q   - %vector of hardening parameter

   @retval The function returns value of yield function for the given stress tensor

   14.1.2004
*/
double shefplast::yieldfunction (vector &sig, vector &q)
{
  double rho, f;

  compute_classvar (sig, q);
  rho = pow(2.0*j2s,0.5)/rc;
  f = rho*rho - khat*rhoc*rhoc /(r*r);

  return f;
}



/**
   This function computes derivatives of r parameter
   q = 1/r i.e. q in this case does not mean hardening parameters vector.
   with respect of tensor sigma.

   @param sig - stress tensor
   @param drds - %vector where the resulting derivatives are stored


   14.1.2004
*/
void shefplast::derqdsigma (vector &/*sig*/, vector &drds)
{
  long i;
  double dqdd0, dqdd1, dqdd2, dd0db0, dd1db0, dd2db0;
  double db1db0, dc0db0,dc1db0,dc2db0,dc3db0,dc4db0,dc5db0;
  double dqdb0, dqdth, db0dxi;
  double dd0dth, dd1dth;  
  double dthdj2, dthdj3;
  double tmp = sqrt(d1*d1-4.0*d0*d2);
  double tmp1, tmp2, tmp3;
  double drhocdxi,drhoedxi,rhoe;
  double pi=4.0*atan(1.0) ;
  double epsi=1.0e-6 ;
  vector dthds(ASTCKVEC(6));
  vector dj3ds(ASTCKVEC(6));
  vector dxids(ASTCKVEC(6));
  vector db0ds(ASTCKVEC(6));

  dqdd0 = 4*d2*d0/2.0/tmp - d1+tmp;
  dqdd0 /= (2.0*d0*d0);
  dqdd1 = (1.0-d1/tmp)/2.0/d0;
  dqdd2 = 1.0/tmp;
  db1db0 = sqrt(3.0)*(1.0-alpha)/(1.0+b0)/(1.0+b0)+2.0*alpha/sqrt(3.0);
 
  tmp1 = -sqrt(3.0)*db1db0*(2.0*b0-sqrt(3.0)*b1)+(2.0-sqrt(3.0)*b1)*(2.0-sqrt(3.0)*db1db0);
  tmp2 = b1*(1.0+b0)-sqrt(3.0)*b0;
  tmp3 = 2.0*tmp2*(db1db0*(1.0+b0)+b1-sqrt(3.0))*(2.0-sqrt(3.0)*b1)*(2.0*b0-sqrt(3.0)*b1);
  
  dc0db0 = (tmp1*tmp2*tmp2 - tmp3)/pow(tmp2, 4.0);
  dc1db0 = -dc0db0*(1.0+b0)*(1.0+b0) - c0*2.0*(1.0+b0);
  dc2db0 = 3.0*dc0db0*(1.0-b0)*(1.0-b0)-6.0*c0*(1.0-b0);
  dc3db0 = 2.0*sqrt(3.0)*dc0db0*(1.0-b0*b0)-4.0*sqrt(3.0)*c0*b0;
  dc4db0 = 1.0-b0*c0-(1.0+b0)*(c0+b0*dc0db0);
  dc5db0 = 3.0*b0*c0-1.0-3.0*(1.0-b0)*(c0+b0*dc0db0);

  dd0db0 = dc1db0*cos(theta)*cos(theta)-dc2db0*sin(theta)*sin(theta)+dc3db0*sin(theta)*cos(theta);
  dd1db0 = 2.0*(dc4db0*sqrt(3.0)*cos(theta)-dc5db0*sin(theta));
  dd2db0 = 4.0 - 3.0*b0*c0 - 3.0*b0*(c0+b0*dc0db0);
  
  dqdb0 = dqdd0*dd0db0 + dqdd1*dd1db0 + dqdd2*dd2db0;

  // Compute dB0ds
  tmp1 = sqrt(m*m-12.0*sqrt(3.0)*m*xi+36.0);
  drhocdxi = -gamma*pow((tmp1-m)/6.0, gamma-1.0);  
  drhocdxi *= m*sqrt(2.0)/tmp1;

  tmp1 = sqrt(m*m-3.0*sqrt(3.0)*m*xi+9.0);
  drhoedxi = -gamma*pow((tmp1-m)/3.0, gamma-1.0);  
  drhoedxi *= m/sqrt(2.0)/tmp1;
  
  tmp = -m+sqrt(m*m-3.0*sqrt(3.0)*m*xi+9.0);
  rhoe = pow(1.0/3.0, gamma)*sqrt(2.0/3.0)*pow(tmp, gamma);

  db0dxi = (drhoedxi*rhoc-drhocdxi*rhoe)/rhoc/rhoc;
  for(i=0; i < 3; i++)
    dxids(i) = 1.0/rc/sqrt(3.0);
  cmulv(db0dxi,dxids,db0ds);
  // End compute dB0ds
  
  dd0dth = -2.0*(c1+c2)*cos(theta)*sin(theta)+c3*(cos(theta)*cos(theta)-sin(theta)*sin(theta));
  dd1dth = 2.0*(-c4*sqrt(3.0)*sin(theta)-c5*cos(theta));
  dqdth = dqdd0*dd0dth+dqdd1*dd1dth;

  // Compute dthds
  /* originall version commented 3.12.2015
    for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
      {
        for (k = 0; k < 3; k++)
          dj3ds[i][j] += dev[i][k]*dev[k][j];
      }
      dj3ds[i][i] -= 2.0/3.0*j2s;
    }
  */
  tensor_dot_prod(dev, dj3ds);
  dj3ds(0) -= 2.0/3.0*j2s;
  dj3ds(1) -= 2.0/3.0*j2s;
  dj3ds(2) -= 2.0/3.0*j2s;
  
  if(fabs(fabs(theta/(pi/6.0))-1.0)<epsi) {
    nullv(dthds);
  }
  else {
    dthdj2 = 3.0 * sqrt(3.0)*j3s/4.0/cos(3.0*theta)/pow(j2s,2.5);
    dthdj3 = -sqrt(3.0)/2.0/cos(3.0*theta)/pow(j2s,1.5); 
    /* originall version commented 3.12.2015
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
	  dthds[i][j] = dthdj2*dev[i][j] + dthdj3*dj3ds[i][j];
      }
    */
    addmultv(dev, dthdj2, dj3ds, dthdj3, dthds);
  }
  // End compute dthds
  
  // The final addition
  /* originall version commented 3.12.2015
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        drds[i][j] = dqdb0*db0ds[i][j]+dqdth*dthds[i][j];
      }
    }
  */
  addmultv(db0ds, dqdb0, dthds, dqdth, drds);
  cmulv((2.0/r)*khat*rhoc*rhoc,drds);
}



/**
   This function computes derivatives of khat parameter
   with respect of tensor sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameter
   @param dkhatds - %matrix where the resulting derivatives are stored
   @param ido - index of internal variables for given material in the ipp other array


   4.1.2002
*/
void shefplast::derkhatds (long /*ipp*/, vector &/*sig*/, vector &q, vector &dkhatds, long /*ido*/)
{
  long i;
  double k, dkhatdxi;

  vector dxids(6);

  k = k0 + (1.0-k0)*sqrt(q[0]*(2.0-q[0]));

  dkhatdxi = -2.0*pow(k,2.0*p)*(1.0-k)*(1.0-k)*xi/a/a;

  for(i=0; i < 3; i++)
    dxids(i) = 1.0/rc/sqrt(3.0);

  // Compute the full term corresponding to dkhatds in dfds
  cmulv(dkhatdxi, dxids, dkhatds);
  cmulv(rhoc*rhoc/r/r, dkhatds);
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameter
   @param dfds - %vector where the resulting derivatives are stored in Voigt notation
   @param ido - index of internal variables for given material in the ipp other array


   4.1.2002
*/
void shefplast::deryieldfsigma (long ipp, vector &sig, vector &q, vector &dfds, long ido)
{
  long i;
  double tmp, tmp1, rhoc;
  vector drhods(ASTCKVEC(6));
  vector drhocds(ASTCKVEC(6));
  vector dqdst(ASTCKVEC(6));
  vector dkhatds(ASTCKVEC(6));
  vector drhoeds(ASTCKVEC(6));

  compute_classvar (sig, q);

  // Compute the full term corresponding to drhods in dfds
  cmulv(2.0/rc/rc, dev, drhods);
  // Compute the full term corresponding to drhocds in dfds
  tmp = sqrt(m*m-12.0*sqrt(3.0)*m*xi+36.0);
  tmp1 = -gamma*pow((tmp-m)/6.0,gamma-1.0);  
  tmp1 *= m*sqrt(2.0)/tmp;
  tmp1 /= (sqrt(3.0)*rc);
  for (i = 0; i < 3; i++)
    drhocds(i) = tmp1;
  rhoc = pow((1.0/6.0),gamma)*sqrt(2.0/3.0)*pow((tmp-m),gamma);
  cmulv(2.0*rhoc*khat/r/r, drhocds);
  
  // Compute the full term corresponding to dqds in dfds
  derqdsigma(sig, dqdst); // OK 2 

  derkhatds(ipp, sig, q, dkhatds, ido);

  addv(dkhatds, drhocds, dfds);
  addv(dfds, dqdst, dfds);
  subv(drhods, dfds, dfds);

  // conversion into Voigt notation
  //dfds(3) *= 2.0;
  //dfds(4) *= 2.0;
  //dfds(5) *= 2.0;
}



/**
   function computes the second derivatives of yield function
   with respect of vector sigma

   @param sig - stress components
   @param ssst - assumed stress/strain state at the given ip

   19.12.2002
*/
void shefplast::dderyieldfsigma (matrix &/*ddfds*/, strastrestate /*ssst*/)
{
}



/**
   This function computes derivatives of plastic potential function
   with respect of vector sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameter
   @param dgds - %vector where the resulting derivatives are stored in Voigt notation
   @param ido - index of internal variables for given material in the ipp other array
*/
void shefplast::derpotsigma (long ipp, vector &sig, vector &q, vector &dgds, long ido)
{
  deryieldfsigma (ipp, sig, q, dgds, ido);
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector of hradening parameters.

   @param sig - stress tensor
   @param dfds - %matrix where the resulting derivatives are stored


   4.1.2002
*/
void shefplast::deryieldfq(vector &/*sig*/, vector &q, vector &dfq)
{
  double k ;
  double dkdkh, dkhatdk ;

  k = k0 + (1.0-k0)*sqrt(q[0]*(2.0-q[0]));
  dkhatdk = 2.0*p*pow(k,2.0*p-1.0)*(1.0-xi*xi*(1.0-k)*(1.0-k)/a/a)
    +2.0*(1.0-k)*(pow(k,p)*xi/a)*(pow(k,p)*xi/a) ;

  dkdkh = (1.0-k0)*(1.0-q[0])/pow(q[0]*(2.0-q[0]),0.5) ;

  dfq[0] = -dkhatdk*dkdkh*rhoc*rhoc/r/r;
  return;
}



/**
   This function computes material stiffnes matrix.

   @param d - allocated matrix structure for material stiffness %matrix
   @param ipp - integration point number

*/
void shefplast::matstiff (matrix &d, long ipp,long ido)
{
  if (Mp->nlman->stmat==initial_stiff) {
    //  initial elastic matrix
    Mm->elmatstiff (d,ipp);
  }
  if (Mp->nlman->stmat==tangent_stiff) {
    //  tangent stiffness matrix

    matrix ad(d.m,d.n);
    Mm->elmatstiff (ad,ipp);
    tangentstiff (ad,d,ipp,ido);
  }
}



/**
   This function returns the tangent stiffness matrix.
   @param d - allocated matrix structure for material stiffness matrix

*/
void shefplast::tangentstiff (matrix &d,matrix &td,long ipp,long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  long i, j;
  double tmps, gamma, hf, dhdk;
  vector mxi, tmpv, str, dfdk, dhds, dfdsdk, dfds, av, q, sig;
  matrix td7;
  matrix a, ainv, dinv, dfdstens, dfdsds, am, tmp;

  gamma=Mm->ip[ipp].eqother[ido+ncomp];

  if (gamma<1.0e-10) {
    copym (d,td);
  }
  else {
    reallocv(RSTCKVEC(ncomp, str));
    reallocv(RSTCKVEC(7, mxi));
    reallocv(RSTCKVEC(7, tmpv));
    reallocv(RSTCKVEC(6, av));
    reallocv(RSTCKVEC(6, dhds));
    reallocv(RSTCKVEC(6, dfdsdk));
    reallocv(RSTCKVEC(6, dfds));
    reallocv(RSTCKVEC(1, dfdk));
    reallocv(RSTCKVEC(1, q));
    reallocv(RSTCKVEC(6, sig));

    reallocm(RSTCKMAT(6, 6, dinv));
    reallocm(RSTCKMAT(7, 7, a));
    reallocm(RSTCKMAT(7, 7, tmp));
    reallocm(RSTCKMAT(7, 7, ainv));
    reallocm(RSTCKMAT(3, 3, dfdstens));
    reallocm(RSTCKMAT(6, 6, dfdsds));
    reallocm(RSTCKMAT(6, 6, am));
    reallocm(RSTCKMAT(d.m+1, d.n+1, td7));

    invm(d, dinv, Mp->zero);

    q[0] = Mm->ip[ipp].eqother[ido+ncomp+1];

    Mm->givestress (0,ipp,str);
    give_full_vector(sig, str, Mm->ip[ipp].ssst);

    deryieldfsigma(ipp, sig, q, dfds, ido);
    numdiff_dfdsdkc(ipp, sig, q, dfdsdk, ido);
    numdiff_dfdsdsc(ipp, sig, q, dfdsds, ido);

    hf = hardening(ipp, sig, q, ido);

    numdiff_dhdkc(ipp, sig, q, dhdk, ido);
    numdiff_dhdsc(ipp, sig, q, dhds, ido);
    deryieldfq(sig, q, dfdk);

    //  block 1,1 (6*6)
    cmulm(gamma, dfdsds, am);
    addm(am, dinv, am);
    
    for (i=0;i<6;i++) {
      for (j=0;j<6;j++) {
	a[i][j]=am[i][j];
      }
    }
    
    //  block 1,2
    cmulv(gamma, dfdsdk, av);
    
    for (i=0;i<6;i++) {
      a[i][6] = av[i];
      mxi[i] = dfds[i] ;
    }
    mxi[6] = dfdk[0] ;
    
    //  block 2,1
    cmulv(-1.0*gamma,dhds,av);
    
    for (i=0;i<6;i++) {
      a[6][i]=av[i];
    }
    
    //  block 2,2
    a[6][6]=1.0-gamma*dhdk;

    invm(a,ainv,Mp->zero);   
    vxm(mxi,ainv,tmpv);

    // We change the last value of vector mxi to get the vector (m,-h)
    mxi[6] = -1.0* hf ; 
    scprd(tmpv,mxi,tmps) ;
    cmulv(1.0/tmps,tmpv,tmpv);
    vxv(mxi,tmpv,tmp);
    mxm(ainv,tmp,td7);
    subm(ainv,td7,td7);

    // We only keep the 6*6 (in 3D) first components for the tangent stiffness matrix
    for (i=0;i<d.m;i++) {
      for (j=0;j<d.n;j++) {
	td[i][j] = td7[i][j] ;
      }
    }
  }
}



/**
   This function computes stresses at given integration point ipp,
   depending on the reached strains.
   The cutting plane algorithm is used. The stress and the other attribute of
   given integration point is actualized.

   @param ipp - integration point number in the mechmat ip array.
   @param ido - index of internal variables for given material in the ipp other array
*/
void shefplast::nlstresses (long ipp, long ido)
//
{
  long i,n=Mm->ip[ipp].ncompstr;
  double gamma;
  vector epsn(n),epsp(n),q(1);

  //  initial values
  for (i=0; i<n; i++) {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  gamma=Mm->ip[ipp].eqother[ido+n];
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  if (q[0] == 0.0)
    q[0] = kh0;

  //  stress return algorithm
  //Mm->cutting_plane (ipp,gamma,epsn,epsp,q);
  stress_return (ipp,gamma,q,epsn,epsp,ido);
  
  //  new data storage
  for (i=0; i<n; i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;
  Mm->ip[ipp].other[ido+n+1]=q[0];
}



/**
   This function updates values in the other array reached in the previous equlibrium state to
   values reached in the new actual equilibrium state.

   @param ipp - integration point number in the mechmat ip array.
   @param ido - index of internal variables for given material in the ipp other array

*/
void shefplast::updateval (long ipp,long ido)
{
  long i,n = Mm->ip[ipp].ncompstr;

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }

  Mm->ip[ipp].eqother[ido+n]=Mm->ip[ipp].other[ido+n];
  Mm->ip[ipp].eqother[ido+n+1]=Mm->ip[ipp].other[ido+n+1];
}



void shefplast::stress_return (long ipp,double &lambda,vector &k,vector &eps,vector &epsp,long ido)
{
  long i,ii,jj,stop,ni;
  double f,hf,nlambda,dlambda,dhdk,nory,err,tmp1,tmp2;
  vector dfds(ASTCKVEC(6)), dfdk(ASTCKVEC(1)), dfdsdk(ASTCKVEC(6));
  vector dhds(ASTCKVEC(6)), nsig(ASTCKVEC(6)), dsig(ASTCKVEC(6));
  vector mxi(ASTCKVEC(7)), tmp(ASTCKVEC(7)), x(ASTCKVEC(7)), y(ASTCKVEC(7));
  vector av(ASTCKVEC(6)), sigtrial(ASTCKVEC(6)), nk(ASTCKVEC(1)), tmp6(ASTCKVEC(6));
  vector nsigtens(ASTCKVEC(6)), sigtens(ASTCKVEC(6)), dfdstens(ASTCKVEC(6));
  matrix d(ASTCKMAT(6,6)), dinv(ASTCKMAT(6,6)), dfdsds(ASTCKMAT(6,6));
  matrix a(ASTCKMAT(7,7)), ainv(ASTCKMAT(7,7)), am(ASTCKMAT(6,6));

  //ni = Mp->nicp;
  //err = Mp->errcp;

  ni=sra.give_ni ();
  err=sra.give_err ();  
  nlambda = 0.0;
  nk[0] = k[0];

  //  elastic stiffness matrix
  Mm->elmatstiff (d,ipp);
  // We try to get the same set up as in ASTER but does not change anything
  /*  for (i=3;i<6;i++) {
      d[i][i] *= 2.0 ;
      eps[i] /= 2.0 ;
      epsp[i] /= 2.0 ;
      }*/
  invm(d,dinv,Mp->zero);
  
  //  elastic strain
  subv(eps,epsp,av);

  //  trial stress
  mxv(d,av,sigtrial);
  
  give_full_vector(nsigtens, sigtrial, Mm->ip[ipp].ssst);
  copyv(sigtrial,nsig);
  // ***********************
  //  main iteration loop
  // ***********************
  stop=0;
  for (i=0;i<ni;i++) {
    //  f computing
    f=yieldfunction (nsigtens,nk);
    // if f is negative at the first step material is still elastic
    if (f<err && i==0) {
      stop =1 ;
      break ;
    }

    //  dfds assembling
    //deryieldfsigma (ipp,nsigtens,nk,dfdstens,ido);
    give_red_vector(dfdstens, dfds, Mm->ip[ipp].ssst);
    
    //  dfdk computing
    //deryieldfq (nsigtens,nk,dfdk);
    
    //  dfdsds assembling
    numdiff_dfdsdsc (ipp,nsigtens,nk,dfdsds,ido);
    
    //  dfdsdk assembling
    numdiff_dfdsdkc (ipp,nsigtens,nk,dfdsdk,ido);
    
    //  h computing
    hf = hardening(ipp,nsigtens,nk,ido);
    
    //  dhds assembling
    numdiff_dhdsc (ipp,nsigtens,nk,dhds,ido);

    //  dhdk computing
    numdiff_dhdkc (ipp,nsigtens,nk,dhdk,ido);
    
    //
    //  Jacobian assembling
    //

    //  block 1,1 (6*6)
    cmulm(nlambda,dfdsds,am);
    addm(am,dinv,am);
    
    for (ii=0;ii<6;ii++) {
      for (jj=0;jj<6;jj++) {
	a[ii][jj]=am[ii][jj];
      }
    }

    //  block 1,2
    cmulv(nlambda,dfdsdk,av);
    
    for (ii=0;ii<6;ii++) {
      a[ii][6]=av[ii];
      mxi[ii] = dfds[ii] ;
    }
    mxi[6] = dfdk[0] ;
    
    //  block 2,1
    cmulv(-1.0*nlambda,dhds,av);
    
    for (ii=0;ii<6;ii++) {
      a[6][ii]=av[ii];
    }
    
    //  block 2,2
    a[6][6]=1.0-nlambda*dhdk;
    // Compute the vector of residuals (7*1)

    //  block 1 (6*1)
    cmulv(nlambda,dfds,av);
    mxv(dinv,nsig,tmp6);
    addv(tmp6,av,av);

    mxv(dinv,sigtrial,tmp6);
    subv(av,tmp6,av);

    for (ii=0;ii<6;ii++) {
      y[ii]=av[ii];
    }
    
    //  block 2
    if(nk[0]==1.0)
      y[6] = 0.0 ;
    else
      y[6]=nk[0]-nlambda*hf-k[0];

    // Compute the Lambda increment dlambda
    invm(a,ainv,Mp->zero);
    vxm(mxi,ainv,tmp);
    scprd(tmp,y,tmp1) ;
    // We change the last value of vector mxi to get the vector (m,-h)
    mxi[6] = -1.0* hf ; 
    scprd(tmp,mxi,tmp2) ;
    dlambda = (f-tmp1)/tmp2 ;
    // Compute the sigma and kappa increment in the x vector
    cmulv(dlambda,mxi,mxi);
    addv(y,mxi,tmp);
    cmulv(-1.0,tmp,tmp);
    mxv(ainv,tmp,x) ;
    
    //
    //  new values
    //
    for (ii=0;ii<6;ii++) {
      dsig[ii]=x[ii];
    }
    
    //  new stresses
    addv(nsig,dsig,nsig);
    give_full_vector(nsigtens, nsig, Mm->ip[ipp].ssst);

    //  new k
    nk[0]+=x[6];
    if(nk[0] < 0.0) {
      nk[0] = 0.0 ;
    }
    if(nk[0] > 1.0) {
      nk[0] = 1.0 ;
    }
    //  new lambda
    nlambda+=dlambda;

      //  norm of the vector of residuals
      nory = fabs(y[0]) ;
      for (ii=1;ii<7;ii++) {
        nory = maxim(nory,fabs(y[ii]));
     }
      nory = maxim(nory,f);

    if(ipp==0)
      fprintf(stderr,"local error %e/%e \n",nory,err);
    if (nory<err) {
      stop=1;
      break;
    }
  }

  if (stop==1) {
    copyv (nk,k);
    Mm->storestress (0,ipp,nsig);
    lambda=nlambda;
    mxv(dinv,nsig,tmp6); 
    subv(eps,tmp6,epsp);
  }
  if(i==ni) {
    fprintf (stderr,"\n\n stress return algo in %s is not successfull after %ld steps for ipp = %ld.\n",__FILE__,ni,ipp);
  }
}



/**
   function computes zeta
   
   14.1.2004
*/
void shefplast::compzeta()
{
  if (xi>0.0)
    /*    zeta = -ah + sqrt(ah*ah+ch);*/
    zeta = 2.0e-5;
  else
    zeta = -ah + sqrt(ah*ah-bh*xi+ch);
}



double shefplast::hardening(long ipp, vector &sigtens, vector &q, long ido)
{

  double h;
  vector dfds;

  if (q[0]<1.0) {
    reallocv (RSTCKVEC(6, dfds));
    deryieldfsigma (ipp, sigtens, q, dfds, ido);
    compzeta ();
    // The tensor norm already computes the square root !!!
    h = sqrt(2.0/3.0)*tensor_stress_norm (dfds)/zeta;
  }
  else {
    h=0.0;
  }

  return h;
}



/**
   function computes numerically the second derivatives of yield function
   with respect to stress tensor
   the first derivatives are expressed explicitly, the second derivatives
   are computed numerically using a centred integration scheme
   
   14.1.2004
*/
void shefplast::numdiff_dfdsdsc(long ipp, vector &sigtens, vector &q, matrix &dfdsds, long ido)
{

  long i,j;
  double dh;
  long ncomp = Mm->ip[ipp].ncompstr;
  vector sig(ASTCKVEC(ncomp)), sigdh(ASTCKVEC(ncomp)), sigmdh(ASTCKVEC(ncomp));
  vector dfds(ASTCKVEC(6)), dfdsdh(ASTCKVEC(6)), sigtensdh(ASTCKVEC(6)), sigtensmdh(ASTCKVEC(6));
  
  give_red_vector(sigtens, sig, Mm->ip[ipp].ssst);
  
  for (i=0;i<6;i++) {
    copyv (sig, sigdh);
    copyv (sig, sigmdh);
    dh = normv(sig);
    dh *= sqrt(DBL_EPSILON);

    sigdh[i]+=dh;
    sigmdh[i]-=dh;
    give_full_vector(sigtensdh, sigdh, Mm->ip[ipp].ssst);
    give_full_vector(sigtensmdh, sigmdh, Mm->ip[ipp].ssst);
    
    deryieldfsigma(ipp, sigtensmdh, q, dfds, ido);
    deryieldfsigma(ipp, sigtensdh, q, dfdsdh, ido);
    
    subv(dfdsdh, dfds, dfds);
    cmulv(0.5/dh, dfds, dfds);
    
    //  originally, wrong limit value of j index used for vector sigdh and number of rows in dfdsds
    // in the below loop and therefor this command was commented out 3.12.2015
    // give_red_vector(dfds, sigdh, Mm->ip[ipp].ssst);
    
    for (j=0;j<6;j++) {
      //  originally, wrong limit value of j index used for vector sigdh and number of rows in dfdsds
      //  this command was commented out 3.12.2015
      //      dfdsds[j][i]=sigdh[j];
      dfdsds[j][i]=dfds[j];
    }
  }
}



/**
   function computes numerically the second derivatives of yield function
   with respect to stress tensor and hardening parameters
   the first derivatives are expressed explicitly, the second derivatives
   are computed numerically
   
   14.1.2004
*/
void shefplast::numdiff_dfdsdkc(long ipp, vector &sigtens, vector &q, vector &dfdsdk, long ido)
{
  double dh;
  vector qdh(ASTCKVEC(1)), qmdh(ASTCKVEC(1));
  vector dfds(ASTCKVEC(6)), dfdsdh(ASTCKVEC(6)), dfdsmdh(ASTCKVEC(6));
  
  dh = normv(q);
  dh *= sqrt(DBL_EPSILON);

  qdh[0]=q[0]+dh;
  qmdh[0]=q[0]-dh;
  
  deryieldfsigma (ipp, sigtens, qmdh, dfds, ido);
  deryieldfsigma (ipp, sigtens, qdh, dfdsdh, ido);
  
  subv(dfdsdh, dfds, dfds);
  cmulv(0.5/dh, dfds, dfds);
  copyv(dfds, dfdsdk);

  // convert to Voigt notation
  //dfdsdk(3) *= 2.0;
  //dfdsdk(4) *= 2.0;
  //dfdsdk(5) *= 2.0;
  
  //  tensor_vector (dfdsdk,dfds,Mm->ip[ipp].ssst,stress);
}



/**
   function computes numerically derivatives of hardening function
   with respect to stress tensor
   
   14.1.2004
*/
void shefplast::numdiff_dhdsc(long ipp, vector &sigtens, vector &q, vector &dhds, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  double h,hdh,dh;
  vector sig(ASTCKVEC(ncomp)), sigdh(ASTCKVEC(ncomp)), sigmdh(ASTCKVEC(ncomp));
  vector sigtensdh(6);
  
  give_red_vector(sigtens, sig, Mm->ip[ipp].ssst);
  
  for (i=0;i<6;i++) {
    copyv(sig, sigdh);
    copyv(sig, sigmdh);
    dh = normv(sig);
    dh *= sqrt(DBL_EPSILON);

    sigdh[i]+=dh;
    sigmdh[i]-=dh;

    give_full_vector(sigtensdh, sigdh, Mm->ip[ipp].ssst);
    give_full_vector(sigtens, sigmdh, Mm->ip[ipp].ssst);
    
    hdh = hardening (ipp, sigtensdh, q, ido);
    h = hardening (ipp, sigtens, q, ido);

    dhds[i]=(hdh-h)*0.5/dh;
  }
}



/**
   function computes numerically derivatives of hardening function
   with respect to hardening parameters
   
   14.1.2004
*/
void shefplast::numdiff_dhdkc(long ipp,vector &sigtens,vector &q,double &dhdk,long ido)
{
  double h,hdh,dh;
  vector qdh(ASTCKVEC(1)), qmdh(ASTCKVEC(1));
  
  dh = normv(q);
  dh *= sqrt(DBL_EPSILON);
  qdh[0]=q[0]+dh;
  qmdh[0]=q[0]-dh;

  hdh = hardening (ipp,sigtens,qdh,ido);
  h = hardening (ipp,sigtens,qmdh,ido);
  
  dhdk = (hdh-h)*0.5/dh;

}



inline double
shefplast::maxim (double a,double b)
{
  return(a>b ? a:b);
}
