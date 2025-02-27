#include "chen.h"
#include "matrix.h"
#include "iotools.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "alias.h"
#include "stochdriver.h"
#include "hardsoft.h"
#include <math.h>

chen::chen (void)
{
  //  compression yield stress
  fyc=0.0;
  //  tension yield stress
  fyt=0.0;
  //  bi-axial compression yield stress
  fybc=0.0;
  //  compression ultimate stress
  fc=0.0;
  //  tension ultimate stress
  ft=0.0;
  //  bi-axial compression ultimate stress
  fbc=0.0;
  
  //  auxiliary constants
  ay=0.0;  au=0.0;
  ky=0.0;  ku=0.0;
  
  
  //  general strain/stress state
  state=1;
  //  compression
  //state=2;
  //  tension
  //state=3;
}

chen::~chen (void)
{

}

/**
   function reads material parameters of the model
   
   JK, 21.2.2005
*/
void chen::read (XFILE *in)
{
  //  compression yield stress
  //  tension yield stress
  //  bi-axial compression yield stress
  //  compression ultimate stress
  //  tension ultimate stress
  //  bi-axial compression ultimate stress
  xfscanf (in,"%lf %lf %lf  %lf %lf %lf",&fyc,&fyt,&fybc,&fc,&ft,&fbc);
  //  type of stress return algorithm
  sra.read (in);
  //  type of hardening/softening
  hs.read (in);
  
}



/**
  The function returns elastic stiffness %matrix.
   
  @param d[out]  - elastic stiffness %matrix
  @param ipp[in] - integration point id
  @param ido[in] - first index in the array other
   
  4.8.2001
*/
void chen::matstiff (matrix &d, long ipp, long ido)
{
  if (Mp->nlman->stmat==initial_stiff){
    //  initial elastic matrix
    Mm->elmatstiff(d, ipp);
  }
  if (Mp->nlman->stmat==tangent_stiff){
    //  tangent stiffness matrix
    matrix ad(d.m, d.n);
    Mm->elmatstiff(ad, ipp);
    tangentstiff(ad, d, ipp, ido);
  }
}



/**
   function assembles tangent stiffness %matrix
   
   @param d - elastic stiffness %matrix of the material
   @param td - tangent stiffness %matrix
   @param ipp - integration point id
   @param ido - first index in the array other
   
*/
void chen::tangentstiff (matrix &d,matrix &td,long ipp,long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  double denom1,denom2,denom,gamma;
  vector str(ncomp),av(d.m),q(1),dfdq(1);
  vector sig(6),dfds(6);
  matrix am(d.m,d.n);
  
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  if (gamma<1.0e-10){
    copym (d,td);
  }
  else{
    
    //  actual stress
    Mm->givestress (0,ipp,str);
    give_full_vector(sig,str,Mm->ip[ipp].ssst);
    //  actual hardening parameter
    q[0]=Mm->ip[ipp].eqother[ido+ncomp+1];
    
    deryieldfdsigma (sig,q,dfds);
    //tensor_vector (str,dfds,Mm->ip[ipp].ssst,stress);
    
    if (Mm->ip[ipp].ssst==planestress){
      vector auxstr(3);
      auxstr[0]=str[0];auxstr[1]=str[1];auxstr[2]=str[2];
      destrv (str);
      allocv (d.m,str);
      str[0]=auxstr[0];str[1]=auxstr[1];str[2]=auxstr[2];
    }
    
    mxv (d,str,av);
    scprd (av,str,denom1);
    
    scprd (str,str,denom2);
    
    deryieldfdq (sig,q,dfdq);
    
    denom2 = sqrt(denom2)*dfdq[0];
    
    denom=denom1-denom2;
    
    //if (fabs(denom)<1.0e-1){
    if (fabs(denom)<1.0e-1 || denom>1.0e5){
      copym (d,td);
    }
    else{
      vxv (str,str,am);
      mxm (d,am,td);
      mxm (td,d,am);
      
      cmulm (1.0/denom,am);
      
      subm (d,am,td);
    }
    
    destrv (str);
  }

}


/**
   function evaluates yield function for given stresses
   
   @param sig - stress components stored in 3x3 %matrix
   @param q - internal parameters (hardening)
   
   JK, 4.8.2001
*/
double chen::yieldfunction (vector &sig,vector &q)
{
  double f,invar,j2,zone;//kappa;
  matrix dev(3,3);
  
  //  computation of stress deviator
  //deviator (sig,dev);
  //  first invariant of the stress tensor
  invar=first_invar (sig);
  //  second invariant of the stress deviator
  //j2=second_stress_invar(dev);
  j2=second_stress_invar(sig);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  if (state==1){
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      
      if (q[0]>1.0)  q[0]=1.0;
      //  yield function with hardening
      f = j2 + ay/3.0*invar - ky + q[0]*((au-ay)*invar/3.0 - (ku-ky));
    }
    //  otherwise
    else{
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      
      if (q[0]>1.0)  q[0]=1.0;
      //  yield function with hardening
      f = j2 - (invar*invar)/6.0 + ay/3.0*invar - ky +q[0]*((au-ay)*invar/3.0 - (ku-ky));
    }
  }
  if (state==2){
    //  compression-compression
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]>1.0)  q[0]=1.0;
    //  yield function with hardening
    f = j2 + ay/3.0*invar - ky + q[0]*((au-ay)*invar/3.0 - (ku-ky));
  }
  if (state==3){
    //  non compression-compression zone
    ay = (fyc-fyt)/2.0;
    au = (fc-ft)/2.0;
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
    
    if (q[0]>1.0)  q[0]=1.0;
    //  yield function with hardening
    f = j2 - (invar*invar)/6.0 + ay/3.0*invar - ky +q[0]*((au-ay)*invar/3.0 - (ku-ky));
  }
  
  return f;
}


/**
   function evaluates derivatives of yield function
   with respect to stress components
   
   @param sig - stress tensor (stored in 3x3 %matrix)
   @param q - hardening parameters
   @param dfds - derivatives of yield function with respect to stress components (stored in 3x3 %matrix)
   
   JK, 20.4.2005
*/
void chen::deryieldfdsigma (vector &sig,vector &/*q*/,vector &dfds)
{
  double invar,j2,zone;//kappa;
  matrix dev(3,3);
  
  //  computation of stress deviator
  //deviator (sig,dev);
  //  first invariant of the stress tensor
  invar=first_invar (sig);
  //  second invariant of the stress deviator
  j2=second_stress_invar(sig);
  
  nullv (dfds);

  zone=sqrt(j2)+invar/sqrt(3.0);
  /*
  if (state==1){
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      
      if (q[0]>1.0)  q[0]=1.0;
      
      dfds[0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
      dfds[1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
      dfds[2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 + q[0]/3.0*(au-ay);
      
      dfds[3]=2*sig[3];
      dfds[4]=2*sig[4];
      dfds[5]=2*sig[5];
    }
    //  otherwise
    else{ 
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      
      if (q[0]>1.0)  q[0]=1.0;
      
      dfds[0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      dfds[1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      dfds[2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      
      dfds[3]=2*sig[3];
      dfds[4]=2*sig[4];
      dfds[5]=2*sig[5];
    }
  }
  if (state==2){
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]>1.0)  q[0]=1.0;
    
    dfds[0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
    dfds[1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
    dfds[2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 + q[0]/3.0*(au-ay);
    
    dfds[3]=2*sig[3];
    dfds[4]=2*sig[4];
    dfds[5]=2*sig[5];
  }
  if (state==3){
    ay = (fyc-fyt)/2.0;
    au = (fc-ft)/2.0;
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
    
    if (q[0]>1.0)  q[0]=1.0;
    
    dfds[0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    dfds[1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    dfds[2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    
    dfds[3]=2*sig[3];
    dfds[4]=2*sig[4];
    dfds[5]=2*sig[5];
    
  }
  */
}

/**
   function computes second derivatives of yiled function with
   respect to stresses
   
   @param sig - stress tensor (stored in 3x3 %matrix)
   @param dfdsds - second derivatives of yield function with respect to stress components (stored in 6x6 %matrix)

   JK, 20.4.2005
*/
void chen::deryieldfdsigmadsigma (vector &sig,matrix &dfdsds)
{
  double invar,j2,zone;
  matrix dev(3,3);

  //  computation of stress deviator
  //deviator (sig,dev);
  //  first invariant of the stress tensor
  invar=first_invar (sig);
  //  second invariant of the stress deviator
  j2=second_stress_invar(sig);
  
  nullm (dfdsds);
  
  zone = sqrt(j2)+invar/sqrt(3.0);

  if (state==1){
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      dfdsds[0][0] =  2.0/3.0;
      dfdsds[0][1] = -1.0/3.0;
      dfdsds[0][2] = -1.0/3.0;
      
      dfdsds[1][0] = -1.0/3.0;
      dfdsds[1][1] =  2.0/3.0;
      dfdsds[1][2] = -1.0/3.0;
      
      dfdsds[2][0] = -1.0/3.0;
      dfdsds[2][1] = -1.0/3.0;
      dfdsds[2][2] =  2.0/3.0;
      
      dfdsds[3][3] = 2.0;
      dfdsds[4][4] = 2.0;
      dfdsds[5][5] = 2.0;
    }
    //  otherwise
    else{ 
      dfdsds[0][0] =  1.0/3.0;
      dfdsds[0][1] = -2.0/3.0;
      dfdsds[0][2] = -2.0/3.0;
      
      dfdsds[1][0] = -2.0/3.0;
      dfdsds[1][1] =  1.0/3.0;
      dfdsds[1][2] = -2.0/3.0;
      
      dfdsds[2][0] = -2.0/3.0;
      dfdsds[2][1] = -2.0/3.0;
      dfdsds[2][2] =  1.0/3.0;
      
      dfdsds[3][3] = 2.0;
      dfdsds[4][4] = 2.0;
      dfdsds[5][5] = 2.0;
    }
  }
  if (state==2){
    dfdsds[0][0] =  2.0/3.0;
    dfdsds[0][1] = -1.0/3.0;
    dfdsds[0][2] = -1.0/3.0;
    
    dfdsds[1][0] = -1.0/3.0;
    dfdsds[1][1] =  2.0/3.0;
    dfdsds[1][2] = -1.0/3.0;
    
    dfdsds[2][0] = -1.0/3.0;
    dfdsds[2][1] = -1.0/3.0;
    dfdsds[2][2] =  2.0/3.0;
    
    dfdsds[3][3] = 2.0;
    dfdsds[4][4] = 2.0;
    dfdsds[5][5] = 2.0;
    
  }
  if (state==3){
    dfdsds[0][0] =  1.0/3.0;
    dfdsds[0][1] = -2.0/3.0;
    dfdsds[0][2] = -2.0/3.0;
    
    dfdsds[1][0] = -2.0/3.0;
    dfdsds[1][1] =  1.0/3.0;
    dfdsds[1][2] = -2.0/3.0;
    
    dfdsds[2][0] = -2.0/3.0;
    dfdsds[2][1] = -2.0/3.0;
    dfdsds[2][2] =  1.0/3.0;
    
    dfdsds[3][3] = 2.0;
    dfdsds[4][4] = 2.0;
    dfdsds[5][5] = 2.0;
  }
}

/**
   function evaluates derivatives of yield function
   with respect to internal variable q

   @param sig - stress tensor (stored in 3x3 %matrix)
   @param q - internal variables
   @param dfdq - derivatives of yield function with respect to internal variable q (stored in 1x1 %matrix)
   
   JK, 20.4.2005
*/
void chen::deryieldfdq (vector &sig,vector &q,vector &dfdq)
{
  double invar,j2,zone;
  matrix dev(3,3);
  
  //  computation of stress deviator
  //deviator (sig,dev);
  //  first invariant of the stress tensor
  invar=first_invar (sig);
  //  second invariant of the stress deviator
  j2=second_stress_invar(sig);

  nullv (dfdq);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  if (state==1){
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      
      if (q[0]<1.0)
	dfdq[0] = (au-ay)*invar/3.0 - (ku-ky);
      else
	dfdq[0]=0.0;
    }
    //  otherwise
    else{
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      
      if (q[0]<1.0)
	dfdq[0] = (au-ay)*invar/3.0 - (ku-ky);
      else
	dfdq[0]=0.0;
    }
  }
  if (state==2){
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]<1.0)
      dfdq[0] = (au-ay)*invar/3.0 - (ku-ky);
    else
      dfdq[0]=0.0;
  }
  if (state==3){
    ay = (fyc-fyt)/2.0;
    au = (fc-ft)/2.0;
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
    
    if (q[0]<1.0)
      dfdq[0] = (au-ay)*invar/3.0 - (ku-ky);
    else
      dfdq[0]=0.0;
  }
}

/**
   function evaluates second derivatives of yield function
   with respect to stress components and internal variables
   
   @param sig - stress tensor (stored in 3x3 %matrix)
   @param q - hardening parameters
   @param dfdsdq - derivatives of yield function with respect to stress components (stored in 6x1 %matrix)
   
   JK, 20.4.2005
*/
void chen::deryieldfdsigmadq (vector &sig,vector &q,matrix &dfdsdq)
{
  double invar,j2,zone;
  matrix dev(3,3);
  
  //  computation of stress deviator
  //deviator (sig,dev);
  //  first invariant of the stress tensor
  invar=first_invar (sig);
  //  second invariant of the stress deviator
  j2=second_stress_invar(sig);
  
  nullm (dfdsdq);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  if (state==1){
    //compression-compression
    if(invar<0.0 && zone<0.0){
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      
      if (q[0]<1.0){
	dfdsdq[0][0]=(au-ay)/3.0;
	dfdsdq[1][0]=(au-ay)/3.0;
	dfdsdq[2][0]=(au-ay)/3.0;
      }

      dfdsdq[3][0]=0.0;
      dfdsdq[4][0]=0.0;
      dfdsdq[5][0]=0.0;
    }
    //otherwise
    else{ 
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      
      if (q[0]<1.0){
	dfdsdq[0][0]=(au-ay)/3.0;
	dfdsdq[1][0]=(au-ay)/3.0;
	dfdsdq[2][0]=(au-ay)/3.0;
      }
      
      dfdsdq[3][0]=0.0;
      dfdsdq[4][0]=0.0;
      dfdsdq[5][0]=0.0;
    }
  }
  if (state==2){
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]<1.0){
      dfdsdq[0][0]=(au-ay)/3.0;
      dfdsdq[1][0]=(au-ay)/3.0;
      dfdsdq[2][0]=(au-ay)/3.0;
    }
    
    dfdsdq[3][0]=0.0;
    dfdsdq[4][0]=0.0;
    dfdsdq[5][0]=0.0;
  }
  if (state==3){
    ay = (fyc-fyt)/2.0;
    au = (fc-ft)/2.0;
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
    
    if (q[0]<1.0){
      dfdsdq[0][0]=(au-ay)/3.0;
      dfdsdq[1][0]=(au-ay)/3.0;
      dfdsdq[2][0]=(au-ay)/3.0;
    }
    
    dfdsdq[3][0]=0.0;
    dfdsdq[4][0]=0.0;
    dfdsdq[5][0]=0.0;
  }
  
}



/**
   function computes derivatives of hardening function with respect to stress components
   
   @param sig - stress components (stored in 3x3 %matrix)
   @param q - internal parameters
   @param dhds - derivatives of hardening functions with respect to stresses (stored in 6 x ncomphard %matrix)
   
   JK, 17.2.2007
*/
void chen::dhdsigma (vector &sigt,vector &q,vector &dhds)
{
  vector dfds(6);
  matrix dfdsds(6,6);

  //  first derivatives of yield function with respect to stresses
  deryieldfdsigma (sigt,q,dfds);
  //  second derivatives of yield function with respect to stresses
  deryieldfdsigmadsigma (sigt,dfdsds);
  //  derivatives of hardening function with respect to stresses
  hs.dhdsigma (sigt,dfds,dfdsds,dhds);
  
  if (q[0]>=1.0)
    nullv (dhds);
}

/**
   function computes derivatives of hardening function with respect to internal parameters

   @param sig - stress components (stored in 3x3 %matrix)
   @param q - internal parameters
   @param dhdq - derivatives of hardening functions with respect to internal parameters (stored in ncomphard x ncomphard %matrix)
   
   JK, 18.2.2007
*/
void chen::dhdqpar (vector &sigt,vector &q,vector &dhdq)
{
  vector dfds(6);
  matrix dfdsdq(6,1);

  //  first derivatives of yield function with respect to stresses
  deryieldfdsigma (sigt,q,dfds);
  //  second derivatives of yield function with respect to stresses
  deryieldfdsigmadq (sigt,q,dfdsdq);
  //  derivatives of hardening function with respect to hardening parameters
  hs.dhdqpar (sigt,dfds,dfdsdq,dhdq);

  if (q[0]>=1.0)
    nullv (dhdq);
}

/**
   function computes derivatives of hardening function with respect to consistency parameter

   @param dhdg - derivatives of hardening function with respect to consistency parameter stored in ncomphard x 1 %vector
   
   JK, 18.2.2007
*/
void chen::dhdgamma (vector &dhdg)
{
  //  derivatives of hardening function with respect to consistency parameter
  hs.dhdgamma (dhdg);
}

/**
   function computes values of hardening function
   
   @param sig - stress components (stored in 3x3 %matrix)
   @param q - internal parameters stored in %vector
   @param h - values of hardening function stored in %vector
   
   JK, 18.2.2007
*/
void chen::hvalues (vector &sigt,vector &q,vector &h)
{
  vector dfds(6);
  matrix dfdsds(6,6);

  //  first derivatives of yield function with respect to stresses
  deryieldfdsigma (sigt,q,dfds);
  //  values of hardening function
  hs.hvalues (sigt,dfds,h);
}



/**
   function computes true stresses from attained strains
   
   @param ipp - number of integration point
   @param im - type of material
   @param ido - index in the array other
   
   21.2.2005
*/
void chen::nlstresses (long ipp, long im, long ido)
{
  long i,ni,n;//nhard;
  double gamma,err;
  
  //  number of strain/stress components
  n = Mm->ip[ipp].ncompstr;
  
  vector epsn(n),epsp(n),q(1);
  
  //  initial values
  for (i=0;i<n;i++){
    //  new total strains
    epsn[i]=Mm->ip[ipp].strain[i];
    //  attained equilibriated plastic starins
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  //  consistency parameter
  gamma = Mm->ip[ipp].eqother[ido+n];
  
  //if (gamma>0.0){
  //printf ("\n ipp %ld  gamma %le",ipp,gamma);
  //}
  
  //  hardening parameter
  q[0]=Mm->ip[ipp].eqother[ido+n+1];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    //Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);

    //Mm->newton_stress_return_2 (ipp,im,ido,gamma,epsn,epsp,q,ni,err,epslim);
    Mm->newton_stress_return_2 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    
    //Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }
  
  //fprintf (Out,"\n q %le         gamma %le",q[0],gamma);
  
  //  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;
  Mm->ip[ipp].other[ido+n+1]=q[0];
  
}

void chen::nonloc_nlstresses (long ipp, long im, long ido)
{
  long i,ni, n = Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(n),epsp(n),q(1);

  //  initial values
  for (i=0;i<n;i++){
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].nonloc[i];
  }
  gamma = Mm->ip[ipp].eqother[ido+n];
  q[0] = Mm->ip[ipp].eqother[ido+n+1];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }
  
  //  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;
  Mm->ip[ipp].other[ido+n+1]=q[0];
}

/**
   function updates values of the array eqother
   
   @param ipp - number of integration point
*/
void chen::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}

/**
   function changes material parameters in stochastic/fuzzy computations
   
   21.2.2005
*/
void chen::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.nba;i++){
    switch (atm.atrib[i]){
    case 0:{
      fyc=val[i];
      break;
    }
    case 1:{
      fyt=val[i];
      break;
    }
    case 2:{
      fybc=val[i];
      break;
    }
    case 3:{
      fc=val[i];
      break;
    }
    case 4:{
      ft=val[i];
      break;
    }
    case 5:{
      fbc=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
  //  hardening/softening
  hs.changeparam (atm,val);
  
}

double chen::give_consparam (long ipp,long ido)
{
  long ncompstr;
  double gamma;
  
  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];
  
  return gamma;
}

long chen::give_num_interparam ()
{
  return 1;
}

void chen::give_interparam (long ipp,long ido,vector &q)
{
  long ncompstr=Mm->ip[ipp].ncompstr;
  
  q[0]=Mm->ip[ipp].eqother[ido+ncompstr+1];
}

