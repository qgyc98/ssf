#include <math.h>
#include "chen.h"
#include "global.h"
#include "intpoints.h"
#include "vecttens.h"
#include "alias.h"
#include "stochdriver.h"
#include "iotools.h"



chen::chen (void)
{
  fyc=0.0;
  fyt=0.0;
  fybc=0.0;
  fc=0.0;
  ft=0.0;
  fbc=0.0;
  epsu=0.0;
  
  ay=0.0;  au=0.0;
  ky=0.0;  ku=0.0;
  alpha=0.0;  beta=0.0;
  
  //  hardening is switched off
  //hard=0;
  //  hardening is switched on
  hard=1;
}

chen::~chen (void)
{

}

/**
   function reads material parameters of the model
   
   21.2.2005
*/
void chen::read (FILE *in)
{
  fscanf (in,"%lf %lf %lf %lf %lf %lf %lf",&fyc,&fyt,&fybc,&fc,&ft,&fbc,&epsu);
  sra.read (in);
}


/**
   function returns elastic stiffness matrix
   
   @param d - elastic stiffness matrix
   
   4.8.2001
*/
void chen::matstiff (matrix &d,long ipp,long ido)
{
  if (Mp->nlman->stmat==initial_stiff){
    //  initial elastic matrix
    Mm->elmatstiff (d,ipp);
  }
  if (Mp->nlman->stmat==tangent_stiff){
    //  tangent stiffness matrix
    //fprintf (stderr,"\n\n tangent stiffness matrix is not implemented yet\n");
    
    matrix ad(d.m,d.n);
    Mm->elmatstiff (ad,ipp);
    tangentstiff (ad,d,ipp,ido);
  }
}


void chen::tangentstiff (matrix &d,matrix &td,long ipp,long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  double denom,gamma;
  vector str,av(d.m),q(1);
  matrix sig(3,3),am(d.m,d.n);
  
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  if (gamma<1.0e-10){
    copym (d,td);
  }
  else{
    
    allocv (ncomp,str);
    
    Mm->givestress (0,ipp,str);
    vector_tensor (str,sig,Mm->ip[ipp].ssst,stress);
    deryieldfdsigma (sig,q,sig);
    tensor_vector (str,sig,Mm->ip[ipp].ssst,stress);
    
    if (Mm->ip[ipp].ssst==planestress){
      vector auxstr(3);
      auxstr[0]=str[0];auxstr[1]=str[1];auxstr[2]=str[2];
      destrv (str);
      allocv (d.m,str);
      str[0]=auxstr[0];str[1]=auxstr[1];str[2]=auxstr[2];
    }
    
    mxv (d,str,av);
    scprd (av,str,denom);
    
    
    q[0] = Mm->ip[ipp].eqother[ido+ncomp+1];
    denom+= plasmodscalar(q);
    
    if (fabs(denom)<1.0e-10){
      copym (d,td);
    }
    else{
      vxv (str,str,am);
      mxm (d,am,td);
      mxm (td,d,am);
      
      cmulm (1.0/denom,am);
      
      subm (d,am,td);
    }
  }
  
}

/**
   function evaluates yield function for given stresses
   
   @param sig - stresses
   @param q - internal parameters (hardening)
   
   4.8.2001
*/
double chen::yieldfunction (matrix &sig,vector &q)
{
  double f,invar,j2,zone,kappa;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  //  compression-compression
  if(invar<0.0 && zone<0.0)
    {
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);
      
      kappa = sqrt(ky + q[0]*(ku-ky));
      if (kappa>sqrt(ku))
	kappa=sqrt(ku);
      
      if (hard==0){
	//  yield function without hardening
	f = j2 + ay/3.0*invar-ky;
      }
      if (hard==1){
	//  yield function with hardening
	f = j2 + beta/3.0*invar-kappa*kappa*(1.0-alpha/3.0*invar);
      }
    }
  //  otherwise
  else
    {
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);
      
      kappa = sqrt(ky + q[0]*(ku-ky));
      if (kappa>sqrt(ku))
	kappa=sqrt(ku);
      
      if (hard==0){
	//  yield function without hardening
	f = j2 - (invar*invar)/6.0+ay/3.0*invar-ky;
      }
      if (hard==1){
	//  yield function with hardening
	f = j2 - (invar*invar)/6.0 + beta/3.0*invar-kappa*kappa*(1.0-alpha/3.0*invar);
      }
      
    }
  //fprintf (Out,"\n a0  %le",a0);
  //fprintf (Out,"\n tau02  %le",tau02);
  
  //  f = tensornorm (dev) - (sqrt(2.0/3.0)*fs+q[0]);
  return f;
}


/**
   function evaluates derivatives of yield function
   with respect of stress components
   4.8.2001
*/
void chen::deryieldfdsigma (matrix &sig,vector &q,matrix &dfds)
{
  double invar,j2,zone,kappa;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  //  compression-compression
  if(invar<0.0 && zone<0.0)
    {
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);
      
      kappa = sqrt(ky + q[0]*(ku-ky));
      if (kappa>sqrt(ku))
	kappa=sqrt(ku);
      
      if (hard==0){
	dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + beta/3.0+alpha/3.0*kappa*kappa;
	dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + beta/3.0+alpha/3.0*kappa*kappa;
	dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + beta/3.0+alpha/3.0*kappa*kappa;
      }
      if (hard==1){
	dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0;
	dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0;
	dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0;
      }

      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];
    }
  
  
  
  //  otherwise
  else
    { 
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);

      kappa = sqrt(ky + q[0]*(ku-ky));
      if (kappa>sqrt(ku))
	kappa=sqrt(ku);
      
      if (hard==0){
	dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + beta/3.0+alpha/3.0*kappa*kappa;
	dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + beta/3.0+alpha/3.0*kappa*kappa;
	dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + beta/3.0+alpha/3.0*kappa*kappa;
      }
      if (hard==1){
	dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0;
	dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0;
	dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0;
      }
      
      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];

      dfds[0][0]-=1.0/3.0*invar;
      dfds[1][1]-=1.0/3.0*invar;
      dfds[2][2]-=1.0/3.0*invar;

      //fprintf (Out,"\n derivace: jsme v otherwise");
    }
  
}

/**
   function evaluates derivatives of yield function
   with respect of stress components
   4.8.2001
*/
void chen::deryieldfdsigma_old (matrix &sig,matrix &dfds)
{
  double invar,j2,zone;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  //compression-compression
  if(invar<0.0 && zone<0.0)
    {
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);

      dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2])+ay/3.0;
      dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2])+ay/3.0;
      dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1])+ay/3.0;

      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];
    }
  
  
  
  //otherwise
  else
    { 
      ay = (fyc-fyt)/2.0;

      dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2])+ay/3.0;
      dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2])+ay/3.0;
      dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1])+ay/3.0;

      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];

      dfds[0][0]-=1.0/3.0*invar;
      dfds[1][1]-=1.0/3.0*invar;
      dfds[2][2]-=1.0/3.0*invar;

      //fprintf (Out,"\n derivace: jsme v otherwise");
    }
  
}


/**
   function computes second derivatives of yiled function with
   respect to stresses
*/
void chen::deryieldfdsigmadsigma (matrix &sig,matrix &dfdsds)
{
  double invar,j2,zone;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  fillm (0.0,dfdsds);
  
  zone = sqrt(j2)+invar/sqrt(3.0);

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

/**
   function evaluates derivatives of yield function
   with respect of stress components
   4.8.2001
*/
void chen::deryieldfdsigma_old_old (matrix &sig,matrix &dfds)
{
  double invar,j2,a0;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  normedtensor (dev,dfds);    
  j2=second_invar(dev);
  
  //compression-compression
  if(invar<0.0 && (sqrt(j2)+invar/sqrt(3.0))<0.0)
    {
      a0 = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      dfds[0][0]+=a0/3.0;
      dfds[1][1]+=a0/3.0;
      dfds[2][2]+=a0/3.0;
    }
  
  
  
  //otherwise
  else
    { 
      a0 = (fyc-fyt)/2.0;
      dfds[0][0]+=a0/3.0-1.0/3.0*invar;
      dfds[1][1]+=a0/3.0-1.0/3.0*invar;
      dfds[2][2]+=a0/3.0-1.0/3.0*invar;
    }
  
}

/**
   function evaluates derivatives of yield function
   with respect of internal variable q
   27.10.2001
*/
void chen::deryieldfdq (matrix &sig,vector &q,vector &dfdq)
{
  double invar,j2,zone,kappa;
  matrix dev(3,3);
  
  if (hard==0){
    nullv (dfdq.a,dfdq.n);
  }
  if (hard==1){
    deviator (sig,dev);
    invar=first_invar (sig);
    j2=second_invar(dev);
    
    zone=sqrt(j2)+invar/sqrt(3.0);
    
    //  compression-compression
    if(invar<0.0 && zone<0.0)
      {
	ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
	au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
	ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
	ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
	alpha = (au-ay)/(ku-ky);
	beta = (ay*ku-au*ky)/(ku-ky);
      }
    //  otherwise
    else{
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);
    }
    
    kappa = sqrt(ky + q[0]*(ku-ky));
    if (kappa>sqrt(ku))
      kappa=sqrt(ku);
    
    dfdq[0]=-2.0*kappa*(1.0-alpha/3.0*invar);
  }
}

/**
   function evaluates derivatives of yield function
   with respect of internal variable q
   27.10.2001
*/
void chen::deryieldfdqdq (matrix &sig,vector &q,matrix &dfdqdq)
{
  double invar,j2,zone,kappa;
  matrix dev(3,3);
  
  if (hard==0){
    dfdqdq[0][0]=0.0;
  }
  if (hard==1){
    deviator (sig,dev);
    invar=first_invar (sig);
    j2=second_invar(dev);
    
    zone=sqrt(j2)+invar/sqrt(3.0);
    
    //compression-compression
    if(invar<0.0 && zone<0.0)
      {
	ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
	au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
	ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
	ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
	alpha = (au-ay)/(ku-ky);
	beta = (ay*ku-au*ky)/(ku-ky);
      }
    else{
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      alpha = (au-ay)/(ku-ky);
      beta = (ay*ku-au*ky)/(ku-ky);
    }
    
    dfdqdq[0][0]=-2.0*(1.0-alpha/3.0*invar);
  }
}


void chen::deryieldfdsigmadq (matrix &sig,vector &q,matrix &dfdsdq)
{
  double invar,j2,zone,kappa;
  matrix dev(3,3);
  
  if (hard==0){
    fillm (0.0,dfdsdq);
  }
  if (hard==1){
    deviator (sig,dev);
    invar=first_invar (sig);
    j2=second_invar(dev);
    
    zone=sqrt(j2)+invar/sqrt(3.0);
    
    //compression-compression
    if(invar<0.0 && zone<0.0)
      {
	ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
	au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
	ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
	ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
	alpha = (au-ay)/(ku-ky);
	beta = (ay*ku-au*ky)/(ku-ky);
	
	kappa = sqrt(ky + q[0]*(ku-ky));
	if (kappa>sqrt(ku))
	  kappa=sqrt(ku);
	
	dfdsdq[0][0]=2.0*alpha/3.0*kappa;
	dfdsdq[1][0]=2.0*alpha/3.0*kappa;
	dfdsdq[2][0]=2.0*alpha/3.0*kappa;
	
	dfdsdq[3][0]=0.0;
	dfdsdq[4][0]=0.0;
	dfdsdq[5][0]=0.0;
	
      }
    
    
    
    //otherwise
    else
      { 
	ay = (fyc-fyt)/2.0;
	au = (fc-ft)/2.0;
	ky = fyc*fyt/6.0;
	ku = fc*ft/6.0;
	alpha = (au-ay)/(ku-ky);
	beta = (ay*ku-au*ky)/(ku-ky);
	
	kappa = sqrt(ky + q[0]*(ku-ky));
	if (kappa>sqrt(ku))
	  kappa=sqrt(ku);
	
	dfdsdq[0][0]=2.0*alpha/3.0*kappa;
	dfdsdq[1][0]=2.0*alpha/3.0*kappa;
	dfdsdq[2][0]=2.0*alpha/3.0*kappa;
	
	dfdsdq[3][0]=0.0;
	dfdsdq[4][0]=0.0;
	dfdsdq[5][0]=0.0;
      }
  }
  
}


void chen::plasmod (long ipp,vector &epsp,matrix &sig,matrix &h)
{
  long ncompstr=epsp.n;  
  double s,eq;
  matrix epspt(3,3);
  
  /*
  //  conversion from vector notation to tensor notation
  vector_tensor (epsp,epspt,Mm->ip[ipp].ssst,strain);
  
  //  equivalent plastic strain
  s  = (epspt[0][0]-epspt[1][1])*(epspt[0][0]-epspt[1][1]);
  s += (epspt[1][1]-epspt[2][2])*(epspt[1][1]-epspt[2][2]);
  s += (epspt[2][2]-epspt[0][0])*(epspt[2][2]-epspt[0][0]);
  s += 6.0*epspt[1][2]*epspt[1][2];
  s += 6.0*epspt[2][0]*epspt[2][0];
  s += 6.0*epspt[0][1]*epspt[0][1];

  eq = sqrt(2.0)/2.0*sqrt(s);
  
  h[0][0]=eq/100.0;
  */


  /*
  //  conversion from vector notation to tensor notation
  vector_tensor (epsp,epspt,Mm->ip[ipp].ssst,strain);
  
  s=Mm->ip[ipp].eqother[ncompstr+2];
  
  fprintf (Out,"\n\n Hodnota pred prirustkem  %le",s);
  cumulstrain (epspt,s);
  fprintf (Out,"\n hodnota po prirustkem    %le",s);
  
  Mm->ip[ipp].eqother[ncompstr+2]=s;
  
  h[0][0]=s;
  */
  
  double f,invar,j2,zone;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  //  compression-compression
  if(invar<0.0 && zone<0.0){
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
  }
  //  otherwise
  else{
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
  }
  
  h[0][0]=(ku-ky)/epsu*epsu;
  

  /*
  //s=Mm->ip[ipp].eqother[ncompstr+2]/10.0;
  s=Mm->ip[ipp].eqother[ncompstr];
  if (s>1.0)
    h[0][0]=50.0;
  else
    h[0][0]=s;
  */
  
  
}

double chen::plasmodscalar(vector &qtr)
{
  /*
  double ret;
  vector dfq(qtr.n), hp(qtr.n);
  matrix h(qtr.n, qtr.n);

  deryieldfdq(sig,q,dfq);
  plasmod (h);
  mxv (h,dfq,hp);
  scprd (hp,dfq,ret);
  return ret;
  */
}

void chen::updateq(long ipp,double dgamma,vector &epsp,matrix &sig,vector &q)
{
  /*
  long i;
  vector dfdq(1),hp(1);
  matrix h(1,1);
  
  if (hard==0){
    nullv (q.a,q.n);
  }
  if (hard==1){
    deryieldfdq (sig,q,dfdq);
    plasmod (ipp,epsp,h);
    mxv (h,dfdq,hp);
    for (i=0;i<1;i++)
      q[i]-=dgamma*hp[i];
  }
*/
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
  long i,ni,n,nhard;
  double gamma,err;
  
  //  number of strain/stress components
  n = Mm->ip[ipp].ncompstr;
  
  vector epsn(n),epsp(n);
  /*
  if (hard==0){
    nhard=0;
  }
  if (hard==1){
    nhard=1;
  }
  */
  vector q(1);
  
  //  initial values
  for (i=0;i<n;i++){
    //  new total strains
    epsn[i]=Mm->ip[ipp].strain[i];
    //  attained equilibriated plastic starins
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  //  consistency parameter
  gamma = Mm->ip[ipp].eqother[ido+n];
  //  hardening parameter
  
  if (hard==0)
    q[0]=0.0;
  if (hard==1){
    
    double f,invar,j2,zone,kappa;
    vector sig (epsp.n);
    matrix sigt(3,3),dev(3,3);
    
    
    vector_tensor (sig,sigt,Mm->ip[ipp].ssst,stress);
    deviator (sigt,dev);
    invar=first_invar (sigt);
    j2=second_invar(dev);
    
    zone=sqrt(j2)+invar/sqrt(3.0);
    
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    }
    //  otherwise
    else{
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
    }
    
    q[0] = Mm->ip[ipp].eqother[ido+n+1]+ky;

  }
  
  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
   
    //Mm->closest_point_proj (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    //Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    Mm->newton_stress_return_2 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    //Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
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
  if (hard==1){
    Mm->ip[ipp].other[ido+n+1]=q[0];
  }
  
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
  
  for (i=0;i<atm.num;i++){
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
