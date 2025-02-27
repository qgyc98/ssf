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
  epsuc=0.0;
  epsut=0.0;
  
  ay=0.0;  au=0.0;
  ky=0.0;  ku=0.0;
  
  hp=0.0;
  
  //  general strain/stress state
  //state=1;
  //  compression
  state=2;
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
void chen::read (FILE *in)
{
  fscanf (in,"%lf %lf %lf  %lf %lf %lf  %lf %lf  %lf",&fyc,&fyt,&fybc,&fc,&ft,&fbc,&epsuc,&epsut,&hp);
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
  double denom1,denom2,denom,gamma;
  vector str,av(d.m),q(1),dfdq(1);
  matrix sig(3,3),am(d.m,d.n),dfds(3,3);
  
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  if (gamma<1.0e-10){
    copym (d,td);
  }
  else{
    allocv (ncomp,str);
    
    //  actual stress
    Mm->givestress (0,ipp,str);
    vector_tensor (str,sig,Mm->ip[ipp].ssst,stress);
    //  actual hardening parameter
    q[0]=Mm->ip[ipp].eqother[ido+ncomp+1];
    
    deryieldfdsigma (sig,q,dfds);
    tensor_vector (str,dfds,Mm->ip[ipp].ssst,stress);
    
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
    
    destrv (str);
  }

}

void chen::plasmod (matrix &h)
{
  h[0][0]=-1.0*hp;
}

/**
   function evaluates yield function for given stresses
   
   @param sig - stresses
   @param q - internal parameters (hardening)
   
   JK, 4.8.2001
*/
double chen::yieldfunction (matrix &sig,vector &q)
{
  double f,invar,j2,zone,kappa;
  matrix dev(3,3);

  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
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
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]>1.0)  q[0]=1.0;
    //  yield function with hardening
    f = j2 + ay/3.0*invar - ky + q[0]*((au-ay)*invar/3.0 - (ku-ky));
  }
  if (state==3){
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
   
   JK, 20.4.2005
*/
void chen::deryieldfdsigma (matrix &sig,vector &q,matrix &dfds)
{
  double invar,j2,zone,kappa;
  matrix dev(3,3);
  
  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  zone=sqrt(j2)+invar/sqrt(3.0);
  
  if (state==1){
    //  compression-compression
    if(invar<0.0 && zone<0.0){
      ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
      au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
      ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
      ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
      
      if (q[0]>1.0)  q[0]=1.0;
      
      dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
      dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
      dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 + q[0]/3.0*(au-ay);
      
      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];
    }
    //  otherwise
    else{ 
      ay = (fyc-fyt)/2.0;
      au = (fc-ft)/2.0;
      ky = fyc*fyt/6.0;
      ku = fc*ft/6.0;
      
      if (q[0]>1.0)  q[0]=1.0;
      
      dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
      
      dfds[0][1]=2*sig[0][1];
      dfds[1][0]=2*sig[0][1];
      dfds[0][2]=2*sig[0][2];
      dfds[2][0]=2*sig[0][2];
      dfds[1][2]=2*sig[1][2];
      dfds[2][1]=2*sig[1][2];
    }
  }
  if (state==2){
    ay = (fybc*fybc-fyc*fyc)/(2.0*fybc-fyc);
    au = (fbc*fbc-fc*fc)/(2.0*fbc-fc);
    ky = fyc*fybc*(2*fyc-fybc)/3.0/(2*fybc-fyc);
    ku = fc*fbc*(2*fc-fbc)/3.0/(2*fbc-fc);
    
    if (q[0]>1.0)  q[0]=1.0;
    
    dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
    dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 + q[0]/3.0*(au-ay);
    dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 + q[0]/3.0*(au-ay);
    
    dfds[0][1]=2*sig[0][1];
    dfds[1][0]=2*sig[0][1];
    dfds[0][2]=2*sig[0][2];
    dfds[2][0]=2*sig[0][2];
    dfds[1][2]=2*sig[1][2];
    dfds[2][1]=2*sig[1][2];
  }
  if (state==3){
    ay = (fyc-fyt)/2.0;
    au = (fc-ft)/2.0;
    ky = fyc*fyt/6.0;
    ku = fc*ft/6.0;
    
    if (q[0]>1.0)  q[0]=1.0;
    
    dfds[0][0]=1.0/3.0*(2.0*sig[0][0]-sig[1][1]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    dfds[1][1]=1.0/3.0*(2.0*sig[1][1]-sig[0][0]-sig[2][2]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    dfds[2][2]=1.0/3.0*(2.0*sig[2][2]-sig[0][0]-sig[1][1]) + ay/3.0 - 1.0/3.0*invar + q[0]/3.0*(au-ay);
    
    dfds[0][1]=2*sig[0][1];
    dfds[1][0]=2*sig[0][1];
    dfds[0][2]=2*sig[0][2];
    dfds[2][0]=2*sig[0][2];
    dfds[1][2]=2*sig[1][2];
    dfds[2][1]=2*sig[1][2];
  }
}

/**
   function computes second derivatives of yiled function with
   respect to stresses
   
   JK, 20.4.2005
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
   
   JK, 20.4.2005
*/
void chen::deryieldfdq (matrix &sig,vector &q,vector &dfdq)
{
  double invar,j2,zone;
  matrix dev(3,3);
  
  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
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

*/
void chen::deryieldfdsigmadq (matrix &sig,vector &q,matrix &dfdsdq)
{
  double invar,j2,zone;
  matrix dev(3,3);
  
  deviator (sig,dev);
  invar=first_invar (sig);
  j2=second_invar(dev);
  
  fillm (0.0,dfdsdq);
  
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

void chen::deryieldfdqdq (matrix &dfdqdq)
{
  dfdqdq[0][0]=0.0;
}

void chen::updateq(long ipp,double dgamma,vector &epsp,matrix &sig,vector &q)
{

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
  //  hardening parameter
  q[0]=Mm->ip[ipp].eqother[ido+n+1];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    //Mm->closest_point_proj (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    //Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    Mm->newton_stress_return_2 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    //Mm->newton_stress_return_3 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
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
    case 8:{
      hp=val[i];
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
