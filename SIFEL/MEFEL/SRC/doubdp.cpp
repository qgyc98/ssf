#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "doubdp.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "elastisomat.h"

/**
  This constructor initializes attributes to zero values.
*/
doubdp::doubdp (void)
{
  fc = 0.0; ft = 0.0; fb = 0.0; gtf = 0.0;
  sft = 0;
}


/**
  The destructor
*/
doubdp::~doubdp (void)
{

}


/**
  This function reads material parameters from the opened text file given
  by the parameter in and evaluates starting parameters of the model.

  @param in - pointer to the opened text file

  06/2014   J. Fiedler
*/
void doubdp::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf",&fc,&ft,&fb,&gtf);
  sra.read (in);
  xfscanf (in,"%ld",&sft);
  
  alphat = 1.73205/3*(fc-ft)/(fc+ft);
  taut   = fc*(1.73205-3*alphat)/3;

  alphac = 1.73205/3*(fb-fc)/(2*fb+fc);
  tauc   = fc*(1.73205-3*alphac)/3;

  icr = (taut-tauc)/(alphat-alphac);
  jcr = taut - alphat*icr;

  ivrch = taut/alphat;
  jvrch = 0.0;
}

/**
   From the given parameter q (hardening parameters) function computes
   new tensile strength and appropriately reevaluates all parameters
   of the material model.

   @param ipp - integration point number
   @param q   - vector of hardening parametrs
                index 0 -> gammat
                index 1 -> q = ft(gammat)

   01/2015  J. Fiedler  
*/
void doubdp::params (long ipp, vector &q)
{
  long eid;
  double h;
  
  eid = Mm->elip[ipp];
  switch (Mt->give_dimension(eid))          // switch for characteristic length of an element
  {
    case 1:
      h = Mt->give_length(eid);
      break;
    case 2:
      h = sqrt(Mt->give_area(eid));
      break;
    case 3:
      h = pow(Mt->give_volume(eid), 1.0/3.0);
      break;
    default:
      print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
      abort();
  }
  
  gamult = gtf*2/ft/h;
  q(1) = ft*(1-q(0)/gamult);
  if (q(1)<0)
    q(1) = 0.0;

  alphat = 1.73205/3*(fc-q(1))/(fc+q(1));
  taut   = fc*(1.73205-3*alphat)/3;

  alphac = 1.73205/3*(fb-fc)/(2*fb+fc);
  tauc   = fc*(1.73205-3*alphac)/3;

  icr = (taut-tauc)/(alphat-alphac);
  jcr = taut - alphat*icr;

  ivrch = taut/alphat;
  jvrch = 0.0;
}

/**
   Function returns value of the yield function which is eventually returned to.
   It performs following tasks:
     - evaluation of the stress invariants
     - dividing the whole region to stress return areas
     - determination of the current return area
            v = 1 -> elastic behaviour
            v = 2 -> return to the tensile yield function
            v = 3 -> return to the compressive yield function
            v = 4 -> return to the intersection
            v = 5 -> return to the top

   @param sig - engineering components of stress, it contains 6 components
          sig[0] = sigma_x
          sig[1] = sigma_y
          sig[2] = sigma_z
          sig[3] = tau_yz
          sig[4] = tau_zx
          sig[5] = tau_xy
   @param q   - %vector of hardening parameters

   06/2014  J. Fiedler
*/
double doubdp::yieldfunction (vector &sig, vector &/*q*/)
{

  double i1, j2, f1, f2, n1, n2, n3;
  double tcr=0.0, tvrch=0.0;
  double a, b, c;
  double l, m, n;
  //  matrix dev(3,3);

  i1 = first_invar (sig);
  //  deviator (sig,dev);
  //  j2 = sqrt(second_invar (dev));
  j2 = sqrt(j2_stress_invar(sig));
  
  f1 = alphat*i1+j2-taut;                   //tensile yield function
  f2 = alphac*i1+j2-tauc;                   //compressive yield function
  n1 = alphac*j2-i1-alphac*jcr+icr;         //compressive normal line of the intersection
  n2 = alphat*j2-i1-alphat*jcr+icr;         //tensile normal line of the intersection
  n3 = alphat*j2-i1-alphat*jvrch+ivrch;     //normal line of the top

  a = i1-icr;
  b = j2-jcr;
  c = -a*icr-b*jcr;
  l = i1-ivrch;
  m = j2-jvrch;
  n = -l*ivrch-m*jvrch;

  if ((f1<=0.0)&&(f2<=0.0)){v = 1; return f1;}
  
  if ((f1>0.0) &&(f2<=0.0))
  {
    if (n3<0)
    {
      v = 5;
      tvrch = l*i1 + m*j2 + n;
      return tvrch;
    }
    else      {v = 2; return f1;}
  }
  
  if ((f1<=0.0)&&(f2>0.0)) {v = 3; return f2;}
  
  if ((f1>0.0) &&(f2>0.0))
  {
    if  (n1>=0)           {v = 3; return f2;}
    if ((n2<=0)&&(n3>=0)) {v = 2; return f1;}
    if ((n1<0) &&(n2>0))
    {
      v = 4;
      tcr = a*i1 + b*j2 + c;
      return tcr;
    }
    if  (n3<0)
    {
      v = 5;
      tvrch = l*i1 + m*j2 + n;
      return tvrch;
    }
  }

  return -1;
}


/**
   This function computes derivatives of the yield function (corresponding
   with the actual stress return area) with respect of vector sigma.

   @param sig - stress tensor
   @param dfds - %matrix where the resulting derivatives are stored
   @param q - %vector of hardening parameters

   @return The function returns resulting %matrix of 
           derivatives in the parameter dfds.

   01/2015 J. Fiedler
*/
void doubdp::deryieldfsigma (vector &sig, vector &dfds, vector &/*q*/)
{
  double i1, j2;
  double a, b, d;
  double l, m;
  //  matrix dev(3,3);

  nullv(dfds);
  i1 = first_invar(sig);
  //  deviator(sig,dev);
  //  j2 = sqrt(second_invar(dev));
  j2 = sqrt(j2_stress_invar(sig));

  a = i1-icr;
  b = j2-jcr;
  l = i1-ivrch;
  m = j2-jvrch;
  /*
  dfds[0][0] += (sig[0][0]-sig[1][1])/3;
  dfds[0][1] += 2*sig[0][1];
  dfds[1][0] += 2*sig[0][1];
  dfds[1][1] += (sig[1][1]-sig[2][2])/3;
  dfds[1][2] += 2*sig[1][2];
  dfds[2][1] += 2*sig[1][2];
  dfds[2][0] += 2*sig[2][0];
  dfds[0][2] += 2*sig[2][0];
  dfds[2][2] += (sig[2][2]-sig[0][0])/3;
  */
  dfds(0) += (sig(0)-sig(1))/3;
  dfds(1) += (sig(1)-sig(2))/3;
  dfds(2) += (sig(2)-sig(0))/3;
  dfds(3) += 2*sig(3);
  dfds(4) += 2*sig(4);
  dfds(5) += 2*sig(5);

  d = 1/(2*j2);
  //  cmulm(d, dfds);
  cmulv(d, dfds);

  // if (v==1)
  if (v==2)
  {
    //    dfds(0,0) += alphat;
    //    dfds(1,1) += alphat;
    //    dfds(2,2) += alphat;
    dfds(0) += alphat;
    dfds(1) += alphat;
    dfds(2) += alphat;
  }
  if (v==3)
  {
    //    dfds(0,0) += alphac;
    //    dfds(1,1) += alphac;
    //    dfds(2,2) += alphac;
    dfds(0) += alphac;
    dfds(1) += alphac;
    dfds(2) += alphac;
  }
  if (v==4)
  {
    //    cmulm(b, dfds);
    //    dfds(0,0) += a;
    //    dfds(1,1) += a;
    //    dfds(2,2) += a;
    cmulv(b, dfds);
    dfds(0) += a;
    dfds(1) += a;
    dfds(2) += a;
  }
  if (v==5)
  {
    //    cmulm(m, dfds);
    //    dfds(0,0) += l;
    //    dfds(1,1) += l;
    //    dfds(2,2) += l;
    cmulv(m, dfds);
    dfds(0) += l;
    dfds(1) += l;
    dfds(2) += l;
  }
}


/**
  This function computes material stiffness %matrix.

  @param d   - material stiffness matrix
  @param ipp - integration point number in the mechmat ip array
  @param ido - index of internal variables for given material in the ipp other array

  @return The function computes material stiffness matrix of the given integration point.

  07/2014 J. Fiedler
*/
void doubdp::matstiff (matrix &d, long ipp,long ido)
{
  if (Mp->nlman->stmat==initial_stiff)
  {
    //  initial elastic matrix
    Mm->elmatstiff(d, ipp, ido);
  }
  if (Mp->nlman->stmat==tangent_stiff)
  {
    //  tangent stiffness matrix
//    matrix ad(d.m,d.n);
    Mm->elmatstiff(d, ipp, ido);
//    tangentstiff (ad,d,ipp,ido);
  }
}


/**
  This function computes stresses at given integration point ipp,
  depending on the reached strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.
  A special approach is used for the plane stress state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the stress array and the other array of the given integration point.

  01/2015 J. Fiedler
*/
void doubdp::nlstresses (long ipp, long im, long ido)
{
  long i,ni,n=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(n),epsp(n),q(2),sig(n),pom(n);
  matrix d(n,n), sigt(3,3);
  
//  initial values
  for (i=0;i<n;i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];          // initial total strain
    epsp[i]=Mm->ip[ipp].eqother[ido+i];     // actual plastic strain
  }
  gamma=Mm->ip[ipp].eqother [ido+n];
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  ni=sra.give_ni ();
  err=sra.give_err ();

// the special approach for plane stress - transverse stress correction
/*  if ((Mm->ip[ipp].ssst == planestress))
  {
    for (i=0;i<n;i++)
      pom[i] = epsn[i];
    epsn[3] = pom[2];
      
    nu = Mm->give_actual_nu(ipp);
    epsn[2] = -nu/(1-nu)*(epsn[0]-epsp[0]+epsn[1]-epsp[1])+epsp[3];     // trial transverse strain
    
    for (j=0;j<ni;j++)      // iteration loop for transverse stress correction
    {
      Mm->cutting_plane2 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
      if (abs(Mm->ip[ipp].stress[3])<(err*100000))      // trasnverse stress testing
        break;
      if (j==ni-1 && abs(Mm->ip[ipp].stress[3])>err*100000)
      {
        print_err("transverse stress correction was not successful",__FILE__,__LINE__,__func__);
        fprintf (Out,"\n transverse stress correction was not successful");
      }
      str = Mm->ip[ipp].stress[3];
      Mm->ip[ipp].ssst = axisymm;
      Mm->elmatstiff (d,ipp);
      Mm->ip[ipp].ssst = planestress;
      epsn[2]+= -Mm->ip[ipp].stress[3]/d(2,2);      // adjustment of the transverse strain      
    }      
  }
// for remaining stress states
else*/
    Mm->cutting_plane3 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);

//  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;
  Mm->ip[ipp].other[ido+n+1]=q[0];  
}


/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the values from the other array to the eqother array.  
*/
void doubdp::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}

/**
  This function computes plastic modulus.

  @param sig - stress tensor
  @param q   - %vector of the hardening parameters

  @return The function returns value of the plastic modulus.

  01/2015 J. Fiedler
*/
double doubdp::plasmodscalar(vector &sig, vector &q)
{

  double ret;

  double dfdq;
  double dqdg;
  double i1 = first_invar(sig);

  dfdq = -(2*sqrt(3)/3)*fc/((fc+q(1))*(fc+q(1)))*(i1+fc);
  dqdg = -ft/gamult;
  ret  = dfdq * dqdg;

  if (sft==0)
    ret = 0;

  return -ret;
}

/**
  This function computes new value of the hardening parameter q
  and calls the "params" function to reevaluate the paramaters
  of the model.

  @param ipp    - integration point pointer
  @param dgamma - increment of the plastic multiplier
  @param q      - %vector of the hardening parameters

  @return The function updates components of the vector of the
          hardening parameters q.

  01/2015 J. Fiedler
*/
void doubdp::updateq(long ipp, double dgamma, vector &q)
{
  if (sft == 1)
  {
    q(0) += dgamma;
    params (ipp,q);
  }
}
