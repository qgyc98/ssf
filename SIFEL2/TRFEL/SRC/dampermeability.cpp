#include <errno.h>
#include "dampermeability.h"
#include "globalt.h"

dampermeability::dampermeability (void)
{
  //  kinematic viscosity of water with respect to temperature in [deg C]
  kinvis = new gfunct (tab,12);
  
  kinvis->tabf->x[0]=0.0;
  kinvis->tabf->y[0]=1.787e-6;

  kinvis->tabf->x[1]=5.0;
  kinvis->tabf->y[1]=1.519e-6;

  kinvis->tabf->x[2]=10.0;
  kinvis->tabf->y[2]=1.307e-6;

  kinvis->tabf->x[3]=20.0;
  kinvis->tabf->y[3]=1.004e-6;

  kinvis->tabf->x[4]=30.0;
  kinvis->tabf->y[4]=0.801e-6;

  kinvis->tabf->x[5]=40.0;
  kinvis->tabf->y[5]=0.658e-6;

  kinvis->tabf->x[6]=50.0;
  kinvis->tabf->y[6]=0.553e-6;

  kinvis->tabf->x[7]=60.0;
  kinvis->tabf->y[7]=0.475e-6;

  kinvis->tabf->x[8]=70.0;
  kinvis->tabf->y[8]=0.413e-6;

  kinvis->tabf->x[9]=80.0;
  kinvis->tabf->y[9]=0.365e-6;

  kinvis->tabf->x[10]=90.0;
  kinvis->tabf->y[10]=0.326e-6;

  kinvis->tabf->x[11]=100.0;
  kinvis->tabf->y[11]=0.294e-6;


  /*
  rhow = new gfunct (tab,40);
  rhow->tabf->x[0]= 0.0;
  rhow->tabf->y[0]=   999.84;
  rhow->tabf->x[1]= 1.0;
  rhow->tabf->y[1]=   999.90;
  rhow->tabf->x[2]= 2.0;
  rhow->tabf->y[2]=   999.94;
  rhow->tabf->x[3]= 3.0;
  rhow->tabf->y[3]=   999.96;
  rhow->tabf->x[4]= 4.0;
  rhow->tabf->y[4]=   999.97;
  rhow->tabf->x[5]= 5.0;
  rhow->tabf->y[5]=   999.96;
  rhow->tabf->x[6]= 6.0;
  rhow->tabf->y[6]=   999.94;
  rhow->tabf->x[7]= 7.0;
  rhow->tabf->y[7]=   999.90;
  rhow->tabf->x[8]= 8.0;
  rhow->tabf->y[8]=   999.85;
  rhow->tabf->x[9]= 9.0;
  rhow->tabf->y[9]=   999.78;
  rhow->tabf->x[10]=10.0;
  rhow->tabf->y[10]=   999.70;
  rhow->tabf->x[11]=11.0;
  rhow->tabf->y[11]=   999.60;
  rhow->tabf->x[12]=12.0;
  rhow->tabf->y[12]=   999.50;
  rhow->tabf->x[13]=13.0;
  rhow->tabf->y[13]=   999.38;
  rhow->tabf->x[14]=14.0;
  rhow->tabf->y[14]=   999.24;
  rhow->tabf->x[15]=15.0;
  rhow->tabf->y[15]=   999.10;
  rhow->tabf->x[16]=16.0;
  rhow->tabf->y[16]=   998.94;
  rhow->tabf->x[17]=17.0;
  rhow->tabf->y[17]=   998.77;
  rhow->tabf->x[18]=18.0;
  rhow->tabf->y[18]=   998.59;
  rhow->tabf->x[19]=19.0;
  rhow->tabf->y[19]=   998.40;
  rhow->tabf->x[20]=20.0;
  rhow->tabf->y[20]=   998.20;
  rhow->tabf->x[21]=21.0;
  rhow->tabf->y[21]=   997.99;
  rhow->tabf->x[22]=22.0;
  rhow->tabf->y[22]=   997.77;
  rhow->tabf->x[23]=23.0;
  rhow->tabf->y[23]=   997.54;
  rhow->tabf->x[24]=24.0;
  rhow->tabf->y[24]=   997.30;
  rhow->tabf->x[25]=25.0;
  rhow->tabf->y[25]=   997.04;
  rhow->tabf->x[26]=26.0;
  rhow->tabf->y[26]=   996.78;
  rhow->tabf->x[27]=27.0;
  rhow->tabf->y[27]=   996.51;
  rhow->tabf->x[28]=28.0;
  rhow->tabf->y[28]=   996.23;
  rhow->tabf->x[29]=29.0;
  rhow->tabf->y[29]=   995.94;
  rhow->tabf->x[30]=30.0;
  rhow->tabf->y[30]=   995.65;
  rhow->tabf->x[31]=31.0;
  rhow->tabf->y[31]=   995.34;
  rhow->tabf->x[32]=32.0;
  rhow->tabf->y[32]=   995.08;
  rhow->tabf->x[33]=33.0;
  rhow->tabf->y[33]=   994.70;
  rhow->tabf->x[34]=34.0;
  rhow->tabf->y[34]=   994.37;
  rhow->tabf->x[35]=35.0;
  rhow->tabf->y[35]=   994.03;
  rhow->tabf->x[36]=40.0;
  rhow->tabf->y[36]=   992.36;
  rhow->tabf->x[37]=50.0;
  rhow->tabf->y[37]=   988.24;
  rhow->tabf->x[38]=60.0;
  rhow->tabf->y[38]=   983.38;
  rhow->tabf->x[39]=70.0;
  rhow->tabf->y[39]=   977.99;
  rhow->tabf->x[40]=80.0;
  rhow->tabf->y[40]=   972.01;
*/
}

dampermeability::~dampermeability (void)
{
  delete kinvis;
}


/**
   function reads parameters
   
   @param in - input file
*/
void dampermeability::read (XFILE */*in*/)
{
}


/**
   function prints parameters
   
   @param out - outut file
*/
void dampermeability::print (FILE */*out*/)
{
}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point

*/
void dampermeability::matcond (matrix &d,long ipp)
{
  long i,n;
  double omega,l,wc,alpha,beta,k,kd,kl,km,t,logk;
  double tmp; 
  
  //  according to Picandet and Dufour
  alpha = 11.3;
  beta = 1.64;
  
  //  size of D matrix
  n = d.n;

  //  permeability
  k = d[0][0];
  
  //  temperature
//  t = Tm->givetransq (temperature,ipp);
  t = 300.0;
  
  //  damage parameter
  omega = Tm->givenontransq (scal_iso_damage,ipp);

  if (omega == 0.0)  // no damage of the material => no modification of the permeability
    return;

  //  process zone length
  l = Tm->givenontransq (proc_zone_length,ipp);
  //  crack width
  wc = Tm->givenontransq (crack_width,ipp);
  
  if (omega<0.15){
    //  zone I
    errno = 0;
    km = k*(1.0+pow((alpha*omega),beta)+pow((alpha*omega),2.0*beta)/2.0);
    check_math_errel(Tm->elip[ipp]);      
  }
/*  
  if (omega>0.9){
    //  zone III
    errno = 0;
    km = wc*wc*wc/12.0/l/kinvis->getval (t-273.15);
    check_math_errel(Tm->elip[ipp]);      
    if (km < k)
      km = k;
  }
  
  if ((0.15 <= omega) && (omega <= 0.9)){
    //  zone II
    
    errno = 0;
    kd = k*(1.0+pow((alpha*omega),beta)+pow((alpha*omega),2.0*beta)/2.0);
    check_math_errel(Tm->elip[ipp]);      
    errno = 0;
    kl = wc*wc*wc/12.0/l/kinvis->getval (t-273.15);
    check_math_errel(Tm->elip[ipp]);      
    
    if (kd>kl)
      km = kd+kl;
    else
    {
      logk = (1.0-omega)*log (kd) + omega *log (kl);    
      errno = 0;
      if (test_math_err())
        km = kd+kl;
        
      errno = 0;
      km = exp (logk);
      if (test_math_err())
        km = kd+kl;
    }
  }
*/
  if (0.15 <= omega){
    //  zone II
    
    tmp = omega;
    if (omega > 0.9)
      omega = 0.9;
    errno = 0;
    kd = k*(1.0+pow((alpha*omega),beta)+pow((alpha*omega),2.0*beta)/2.0);
    check_math_errel(Tm->elip[ipp]);      
    errno = 0;
    kl = wc*wc*wc/12.0/l/kinvis->getval (t-273.15);
    check_math_errel(Tm->elip[ipp]);      
    
    if (kd>kl)
      km = kd+kl;
    else
    {
      logk = (1.0-omega)*log (kd) + omega *log (kl);    
      errno = 0;
      if (test_math_err())
        km = kd+kl;
        
      errno = 0;
      km = exp (logk);
      if (test_math_err())
        km = kd+kl;
    }
    omega = tmp;
  }
  if (omega > 0.9)
  {
    kd = wc*wc*wc/12.0/l/kinvis->getval (t-273.15);
    if (kd > km)
      km = kd;
  }
  for (i=0;i<n;i++){
    d[i][i]=km;
  }
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   4. 2. 2014 by TKo
*/
void dampermeability::give_reqntq(long *antq)
{
  //  damage parameter
  antq[scal_iso_damage-1] = 1;
}



