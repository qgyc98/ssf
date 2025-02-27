#include <math.h>
#include "j2flow.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include "alias.h"
#include "stochdriver.h"
#include "matrix.h"

j2flow::j2flow (void)
{
  //  yiled stress
  fs=0.0;
  //  plastic modulus
  //  Simo, Hughes: Computational inelasticity, page 9
  k=0.0;
}

j2flow::~j2flow (void)
{

}

/**
   function reads model parametrs
   
   @param in - input file
   
   modified 17. 4. 2015, JK
*/
void j2flow::read (XFILE *in)
{
  //  fs - yield stress
  //  k - plastic modulus
  xfscanf (in,"%k%lf %k%lf","fs",&fs,"k",&k);
  //  stress return algorithm
  sra.read (in);
  
  //  magnitude of shear stress at yielding in pure shear
  //  Jirasek, Bazant
}


/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by T. Koudelka 2.7.2015
*/
void j2flow::print (FILE *out)
{
  fprintf (out,"%lf %lf ", fs, k);
  sra.print (out);
}

/**
   function evaluates yield function for given stresses
   
   @param sig - engineering components of stress, it contains 6 components
          sig[0] = sigma_x
          sig[1] = sigma_y
          sig[2] = sigma_z
          sig[3] = tau_yz
          sig[4] = tau_zx
          sig[5] = tau_xy
   @param q - internal variables (hardening)
   
   4.8.2001
*/
double j2flow::yieldfunction (vector &sig,vector &q)
{
  double j2,f;
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  //  components of the deviator are not needed
  j2 = j2_stress_invar (sig);
  
  //  yield function
  //f = sqrt(j2) - (fs+k*q[0]);
  f = sqrt(j2) - (fs+q[0]);
  
  return f;
}

/**
   function evaluates derivatives of the yield function
   with respect to stress components
   
   @param sigma - stress components
   @param dfds - %vector of derivatives of the yield function
   
   4.8.2001, modified 17. 4. 2015, JK
*/
void j2flow::dfdsigma (vector &sig,vector &dfds)
{
  double j2;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j2 = sqrt(j2);
  
  //  deviator of stress tensor
  deviator (sig,dev);
  
  dfds[0] = dev[0]/j2/2.0;
  dfds[1] = dev[1]/j2/2.0;
  dfds[2] = dev[2]/j2/2.0;
  dfds[3] = dev[3]/j2;
  dfds[4] = dev[4]/j2;
  dfds[5] = dev[5]/j2;
}

/**
   function evaluates derivatives of the yield function
   with respect to internal variable q
   
   @param dq - derivatives of the yield function with respect to the internal variable
   
   27.10.2001
*/
void j2flow::dfdqpar (vector &dq)
{
  dq[0]=-1.0;
}


/**
   function assembles second derivatives of the yield function
   with respect to the stress tensor
   
   @param sig - %vector of stress components
   @param dfdsds - %matrix of the derivatives
   
   20. 4. 2015, JK
*/
void j2flow::dfdsigmadsigma (vector &sig,matrix &dfdsds)
{
  double j2,j23;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j23 = j2*j2*j2;
  j2 = sqrt(j2);
  j23 = sqrt(j23);
  
  //  deviator of stress tensor
  deviator (sig,dev);
  
  //  only one half of components is evaluated
  //  due to symmetry

  dfdsds[0][0] = 0.0 - dev[0]*dev[0]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[0][1] = 0.0 - dev[0]*dev[1]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[0][2] = 0.0 - dev[0]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[0][3] = 0.0 - dev[0]*sig[3]/2.0/j23;
  dfdsds[0][4] = 0.0 - dev[0]*sig[4]/2.0/j23;
  dfdsds[0][5] = 0.0 - dev[0]*sig[5]/2.0/j23;

  dfdsds[1][1] = 0.0 - dev[1]*dev[1]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[1][2] = 0.0 - dev[1]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[1][3] = 0.0 - dev[1]*sig[3]/2.0/j23;
  dfdsds[1][4] = 0.0 - dev[1]*sig[4]/2.0/j23;
  dfdsds[1][5] = 0.0 - dev[1]*sig[5]/2.0/j23;

  dfdsds[2][2] = 0.0 - dev[2]*dev[2]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[2][3] = 0.0 - dev[2]*sig[3]/2.0/j23;
  dfdsds[2][4] = 0.0 - dev[2]*sig[4]/2.0/j23;
  dfdsds[2][5] = 0.0 - dev[2]*sig[5]/2.0/j23;

  dfdsds[3][3] = 0.0 - sig[3]*sig[3]/j23 + 1.0/j2;
  dfdsds[3][4] = 0.0 - sig[3]*sig[4]/j23;
  dfdsds[3][5] = 0.0 - sig[3]*sig[5]/j23;

  dfdsds[4][4] = 0.0 - sig[4]*sig[4]/j23 + 1.0/j2;
  dfdsds[4][5] = 0.0 - sig[4]*sig[5]/j23;

  dfdsds[5][5] = 0.0 - sig[5]*sig[5]/j23 + 1.0/j2;
  

  //  second half of components are copied from
  //  the first half

  dfdsds[1][0] = dfdsds[0][1];

  dfdsds[2][0] = dfdsds[0][2];
  dfdsds[2][1] = dfdsds[1][2];

  dfdsds[3][0] = dfdsds[0][3];
  dfdsds[3][1] = dfdsds[1][3];
  dfdsds[3][2] = dfdsds[2][3];

  dfdsds[4][0] = dfdsds[0][4];
  dfdsds[4][1] = dfdsds[1][4];
  dfdsds[4][2] = dfdsds[2][4];
  dfdsds[4][3] = dfdsds[3][4];

  dfdsds[5][0] = dfdsds[0][5];
  dfdsds[5][1] = dfdsds[1][5];
  dfdsds[5][2] = dfdsds[2][5];
  dfdsds[5][3] = dfdsds[3][5];
  dfdsds[5][4] = dfdsds[4][5];

}

/**
   function evaluates derivatives of the yield function
   with respect to internal variable q
   
   @param dq - derivatives of the yield function with respect to the internal variable
   
   27.10.2001
*/
void j2flow::dfdsigmadq (matrix &dfdsdq)
{
  nullm (dfdsdq);
}

/**
   function evaluates derivates of the yield function with respect to stresses and
   internal variables (hardening/softening parameters)
   
   @param dfdsdq - derivatives of yield function with respect to stresses and internal variables
   
   JK, 8.8.2005
*/
void j2flow::dfdqpardqpar (matrix &dfdqdq)
{
  nullm (dfdqdq);
}

/**
   function assembles the %vector of hardening h which governs
   the evolution of internal parameters
   
   17. 6. 2015
*/
void j2flow::hardvect (vector &hv)
{
  //hv[0]=1.0;
  hv[0]=(0.0-k)/sqrt(3.0);
  //hv[0]=0.0-k;
}

/**
   function returns elastic stiffness %matrix
   
   @param d - elastic stiffness %matrix
   @param ipp - integration point id
   @param ido - 
   
   4.8.2001
*/
void j2flow::matstiff (matrix &d,long ipp,long ido)
{
  if (Mp->nlman->stmat==initial_stiff){
    //  initial elastic matrix
    Mm->elmatstiff(d, ipp, ido);
  }
  else{
    //  tangent stiffness matrix
    
    matrix ad(d.m,d.n);
    Mm->elmatstiff(ad, ipp, ido);
    tangentstiff(ad, d, ipp, ido);
  }
}

/**
  The function assembles the elasto-plastic stiffness %matrix
   
  @param d[out] - elastic stiffness %matrix
  @param td[out] - elasto-plastic stiffness %matrix
  @param ipp[in] - integration point id
  @param ido[in] - index of internal variables for given material in the ipp other array
     
  20. 4. 2015
*/
void j2flow::tangentstiff (matrix &d,matrix &td,long ipp,long ido)
{
  long ncomp, i, j;

  //  the number of stress components
  ncomp=Mm->ip[ipp].ncompstr;

  double denom,gamma;
  vector dfds(ASTCKVEC(ncomp)),sig(ASTCKVEC(ncomp)),av(ASTCKVEC(ncomp)),q(ASTCKVEC(1));
  matrix de(ASTCKMAT(ncomp, ncomp));
  matrix auxm1(ASTCKMAT(ncomp, ncomp));
  matrix auxm2(ASTCKMAT(ncomp, ncomp));

  
  
  //  the consistency parameter
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  if (gamma<1.0e-10){
    //  there is no plastic strain
    //  elastic matrix is used
    copym (d,td);
  }
  else{
    strastrestate ssst = Mm->ip[ipp].ssst;
    if ((ssst == planestrain) || (ssst == planestress))
    {
      for(i=0; i<3; i++)
        for (j=0; j<3; j++)
          de(i,j) = d(i,j);
    }
    else
      copym(d, de);
   
    //  actual stress components
    Mm->givestress (0,ipp,sig);
    //  derivative of the yiled function with respect to the stress
    dfdsigma (sig,dfds);
    
    //  df/ds D df/ds
    mxv (de,dfds,av);
    scprd (av,dfds,denom);
    
    //  hardening variable
    q[0] = Mm->ip[ipp].eqother[ido+ncomp+1];
    denom+= plasmodscalar ();
    
    if (fabs(denom)<1.0e-10){
      //  there is no plastic strain
      //  elastic matrix is used
      copym (d,td);
    }
    else{
      vxv (dfds,dfds,auxm1);
      mxm (de,auxm1,auxm2);
      mxm (auxm2,de,auxm1);
      
      if ((de(0,0)-auxm1(0,0)/denom)/de(0,0) > 0.01)
      {  
        cmulm (1.0/denom,auxm1);
        subm (de,auxm1,auxm2);
      }
      if ((ssst == planestrain) || (ssst == planestress))
      {
        for(i=0; i<3; i++)
          for (j=0; j<3; j++)
            td(i,j) = auxm2(i,j);
      }
      else
        copym(auxm2, td);
    } 
  }
}

void j2flow::updateq(double dgamma, vector &q)
{
  vector dfdq(1);
  vector hv(1);
  
  //dfdqpar (dfdq);
  hardvect (hv);
  
  //q[0]-=dgamma*dfdq[0];
  q[0]-=dgamma*hv[0]/sqrt(3.0);
}


/**
   function computes df/dq h
   this number if used e.g. in the cutting plane method
   
   20. 4. 2015, JK
*/
double j2flow::plasmodscalar ()
{
  double ret;
  //vector dfdq(1);
  
  // derivative of the yield function with respect to the hardening variable
  //dfdqpar (dfdq);
  // scprd (dfdq,dfdq,ret);

  vector dfq(1),hv(1);
  
  hardvect (hv);
  dfdqpar (dfq);
  ret = ss(dfq.a,hv.a,hv.n);

  return ret;
}

/**
   function computes stresses
*/
void j2flow::nlstresses (long ipp, long im, long ido)
  //
{
  long i,ni, n = Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(ASTCKVEC(n)),epsp(ASTCKVEC(n)),q(ASTCKVEC(1));

  //  initial values
  for (i=0;i<n;i++){
    //  total strains
    epsn[i]=Mm->ip[ipp].strain[i];
    //  plastic strains
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  //  consistency parameter
  gamma = Mm->ip[ipp].eqother[ido+n];
  //  hardening variable
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  
  //  stress return algorithm
  switch (sra.tsra){
  case cp:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  case gsra:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  default:{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  }
   
  //  new data storage
  for (i=0;i<n;i++){
    //  plastic strains
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  //  consistency parameter
  Mm->ip[ipp].other[ido+n]=gamma;
  //  hardening variable
  Mm->ip[ipp].other[ido+n+1]=q[0];

}









void j2flow::nonloc_nlstresses (long ipp, long im, long ido)
{
  long i,ni, n = Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(ASTCKVEC(n)),epsp(ASTCKVEC(n)),q(ASTCKVEC(1));

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
   function copies array other to the array eqother
   
   @param ipp - number of integration point
   @param im - index of material type
   @param ido - index in array other
   
*/
void j2flow::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  Returns vector of irreversible strains via parameter epsp
*/
void j2flow::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



void j2flow::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      fs=val[i];
      break;
    }
    case 1:{
      k=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}

double j2flow::give_consparam (long ipp,long ido)
{
  long ncompstr;
  double gamma;
  
  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];
  
  return gamma;
}

long j2flow::give_num_interparam ()
{
  return 1;
}

void j2flow::give_interparam (long ipp,long ido,vector &q)
{
  long ncompstr=Mm->ip[ipp].ncompstr;
  
  q[0]=Mm->ip[ipp].eqother[ido+ncompstr+1];
}
