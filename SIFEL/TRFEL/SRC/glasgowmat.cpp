#include "globalt.h"
#include "glasgowmat.h"
#include "globmatt.h"

/****************************
*                           *
* Tenchev Model - 8/3/04    *
*                           *
****************************/

glasgowmat::glasgowmat ()
{
  model=2;              //  1=original tenchev, 2=modified tenchev
  ra=287.0;
  rv=461.5;
  stef=5.67e-8;
  
  rhoc = 0.0;
  rhos = 0.0;
  por0 = 0.0;
  k0 = 0.0;
  rhol0 = 0.0;
  sssp = 0.0;
  t0 = 0.0;
  rhov0 = 0.0;
  pginf = 0.0;
  emmi = 0.0;
  alph = 0.0;
  hq = 0.0;
  crhoair = 0.0;
  tfirestart = 0.0;
}

glasgowmat::~glasgowmat ()
{}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void glasgowmat::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void glasgowmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pg,t,rhov;
  k = 0.0;
  
  rhov = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = kmv (t,pg,rhov);
  if((ri == 0) && (ci == 1))
    k = kmp (t,pg,rhov);
  if((ri == 0) && (ci == 2))
    k = kmt (t,pg,rhov);

  if((ri == 1) && (ci == 0))
    k = kav (t,pg,rhov);
  if((ri == 1) && (ci == 1))
    k = kap (t,pg,rhov);
  if((ri == 1) && (ci == 2))
    k = kat (t,pg,rhov);

  if((ri == 2) && (ci == 0))
    k = ktv (t,pg,rhov);
  if((ri == 2) && (ci == 1))
    k = ktp (t,pg,rhov);
  if((ri == 2) && (ci == 2))
    k = ktt (t,pg,rhov);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void glasgowmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pg,t,rhov;
  k = 0.0;
  
  rhov = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = kmv (t,pg,rhov);
  if((ri == 0) && (ci == 1))
    k = kmp (t,pg,rhov);
  if((ri == 0) && (ci == 2))
    k = kmt (t,pg,rhov);

  if((ri == 1) && (ci == 0))
    k = kav (t,pg,rhov);
  if((ri == 1) && (ci == 1))
    k = kap (t,pg,rhov);
  if((ri == 1) && (ci == 2))
    k = kat (t,pg,rhov);

  if((ri == 2) && (ci == 0))
    k = ktv (t,pg,rhov);
  if((ri == 2) && (ci == 1))
    k = ktp (t,pg,rhov);
  if((ri == 2) && (ci == 2))
    k = ktt (t,pg,rhov);
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
}

/**
   function creates conductivity matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void glasgowmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pg,t,rhov;
  k = 0.0;
  
  rhov = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = kmv (t,pg,rhov);
  if((ri == 0) && (ci == 1))
    k = kmp (t,pg,rhov);
  if((ri == 0) && (ci == 2))
    k = kmt (t,pg,rhov);

  if((ri == 1) && (ci == 0))
    k = kav (t,pg,rhov);
  if((ri == 1) && (ci == 1))
    k = kap (t,pg,rhov);
  if((ri == 1) && (ci == 2))
    k = kat (t,pg,rhov);

  if((ri == 2) && (ci == 0))
    k = ktv (t,pg,rhov);
  if((ri == 2) && (ci == 1))
    k = ktp (t,pg,rhov);
  if((ri == 2) && (ci == 2))
    k = ktt (t,pg,rhov);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}


/**
   function creates capacity matrix of the material

   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void glasgowmat::matcap (double &c,long ri,long ci,long ipp)
{
  double pg,t,rhov;
  c = 0.0;
  
  rhov = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    c = cmv (t,pg,rhov);
  if((ri == 0) && (ci == 1))
    c = cmp ();
  if((ri == 0) && (ci == 2))
    c = cmt (t,pg,rhov);

  if((ri == 1) && (ci == 0))
    c = cav (t,pg,rhov);
  if((ri == 1) && (ci == 1))
    c = cap (t,pg,rhov);
  if((ri == 1) && (ci == 2))
    c = cat (t,pg,rhov);

  if((ri == 2) && (ci == 0))
    c = ctv (t,pg,rhov);
  if((ri == 2) && (ci == 1))
    c = ctp ();
  if((ri == 2) && (ci == 2))
    c = ctt (t,pg,rhov);
}




void glasgowmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&rhoc,&rhos,&por0,&k0,&rhol0,&sssp,&t0,&rhov0,&pginf);
  xfscanf (in,"%lf %lf %lf %lf %lf",&emmi,&alph,&hq,&crhoair,&tfirestart);
}



void glasgowmat::print (FILE *out)
{
  fprintf (out,"  %e %e %e %e %e %e %e %e %e",rhoc,rhos,por0,k0,rhol0,sssp,t0,rhov0,pginf);
  fprintf (out,"%e %e %e %e %e",emmi,alph,hq,crhoair,tfirestart);
}

/**  Calculate Far Field Vapour Content
     @param rhov0   - 
     @retval rhoinf - 
*/
double glasgowmat::f_rhovinf ()
{
  double rhovinf;

  rhovinf=rhov0*0.8;
  
  return(rhovinf);
}
  
/**  Calculate Far Field Temperature
     @param time   - 
     @param t0     - 
     @retval tinf  - 
*/
double glasgowmat::f_tinf (double time)
{
  double tinf;

  //Linear Temperature increase
  //tinf=(t0-273.15)+(time*20.0/60.0)+273.15;  
  
  //ISO 834 Fire Curve
  tinf=((t0-273.15)+345.0*log10(8.0*(time-tfirestart)/60.0+1.0))+273.15;  
  
  return(tinf);
}
  
/**  Calculate Vapour Pressure
     @param    - 
     @retval   - 
*/
double glasgowmat::f_pv (double rhov, double t)
{
  double pv;
 
  pv=rhov*rv*t;

  return(pv);
}
    
/**  Calculate Saturation Vapour Pressure
     @param    - 
     @retval   - 
*/
double glasgowmat::f_psat(double t)
{
  double psat;

  psat  = -0.000000001437422219446870*t*t*t*t*t*t +
    0.000004424390583021230*t*t*t*t*t -
    0.003928080821257910*t*t*t*t +
    1.591032529443030*t*t*t -
    325.8874385048470*t*t +
    32147.79527519750*t -
    1154663.603250870;
  
  return(psat);
}

/**  Calculate Derivative of Saturation Vapour Pressure wrt Temperature
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dpsatdt (double t)
{
  double dpsatdt;
  
  dpsatdt=32147.7952751975 - 
    651.774877009694*t + 
    4.77309758832909*t*t - 
    0.01571232328503164*t*t*t + 
    0.00002212195291510615*t*t*t*t - 
    8.62453331668122e-9*t*t*t*t*t;
  
  return(dpsatdt);
}

/**  Calculate Density of Liquid Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_rhow (double /*t*/)
{
  double rhow;

  //Tenchev
  //rhow0=rhow;
  rhow=1000.0;

  /*   //Gawin et al.
       
       //Temperature in Centigrade
       TC=T-273.15;
       TC0=T0-273.15;

       pw0 = 1.0e7;  // Water pressure in range 0 - 20MPa
       
       rhow0= (4.8863e-7 - 1.6528e-9*TC0 + 1.8621e-12*(TC0*TC0) +
       2.4266e-13*(TC0*TC0*TC0) - 1.5996e-15*(TC0*TC0*TC0*TC0) +
       3.3703e-18*(TC0*TC0*TC0*TC0*TC0))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC0 + 8.7696e-3*(TC0*TC0)-
       9.2118e-5*(TC0*TC0*TC0) + 3.3534e-7*(TC0*TC0*TC0*TC0) -
       4.4034e-10*(TC0*TC0*TC0*TC0*TC0));
       
       
       if(TC <= 373.946){
       rhow= (4.8863e-7 - 1.6528e-9*TC + 1.8621e-12*(TC*TC) +
       2.4266e-13*(TC*TC*TC) - 1.5996e-15*(TC*TC*TC*TC) +
       3.3703e-18*(TC*TC*TC*TC*TC))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC + 8.7696e-3*(TC*TC)-
       9.2118e-5*(TC*TC*TC) + 3.3534e-7*(TC*TC*TC*TC) -
       4.4034e-10*(TC*TC*TC*TC*TC));
       
       //Calculate Derivative of Density of Liquid Water wrt Temperature
       drhowdT=-0.7737700000000001 + (-2.0e7 + pw0)*(-1.6528000000000001e-9 + 
       3.7242e-12*(TC) + 7.2798e-13*(TC*TC) - 6.3984e-15*(TC*TC*TC) + 
       1.68515e-17*(TC*TC*TC*TC)) + 0.0175392*(TC) - 
       0.00027635400000000003*(TC*TC) + 1.3413600000000001e-6*(TC*TC*TC) - 
       2.2017000000000004e-9*(TC*TC*TC*TC);
       }
       else{
       rhow= 6.18545405822237e-06*(pw0 - 2.0e7) + 478.725113592986;
       drhowdT=0.0;
       }
  */
  
  return(rhow);
}
  

/**  Calculate Density of Liquid Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_rhow0 (double /*t*/)
{
  double rhow0;

  //Tenchev
  //rhow0=rhow;
  rhow0=1000.0;

  /*   //Gawin et al.
       
       //Temperature in Centigrade
       TC=T-273.15;
       TC0=T0-273.15;

       pw0 = 1.0e7;  // Water pressure in range 0 - 20MPa
       
       rhow0= (4.8863e-7 - 1.6528e-9*TC0 + 1.8621e-12*(TC0*TC0) +
       2.4266e-13*(TC0*TC0*TC0) - 1.5996e-15*(TC0*TC0*TC0*TC0) +
       3.3703e-18*(TC0*TC0*TC0*TC0*TC0))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC0 + 8.7696e-3*(TC0*TC0)-
       9.2118e-5*(TC0*TC0*TC0) + 3.3534e-7*(TC0*TC0*TC0*TC0) -
       4.4034e-10*(TC0*TC0*TC0*TC0*TC0));
       
       
       if(TC <= 373.946){
       rhow= (4.8863e-7 - 1.6528e-9*TC + 1.8621e-12*(TC*TC) +
       2.4266e-13*(TC*TC*TC) - 1.5996e-15*(TC*TC*TC*TC) +
       3.3703e-18*(TC*TC*TC*TC*TC))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC + 8.7696e-3*(TC*TC)-
       9.2118e-5*(TC*TC*TC) + 3.3534e-7*(TC*TC*TC*TC) -
       4.4034e-10*(TC*TC*TC*TC*TC));
       
       //Calculate Derivative of Density of Liquid Water wrt Temperature
       drhowdT=-0.7737700000000001 + (-2.0e7 + pw0)*(-1.6528000000000001e-9 + 
       3.7242e-12*(TC) + 7.2798e-13*(TC*TC) - 6.3984e-15*(TC*TC*TC) + 
       1.68515e-17*(TC*TC*TC*TC)) + 0.0175392*(TC) - 
       0.00027635400000000003*(TC*TC) + 1.3413600000000001e-6*(TC*TC*TC) - 
       2.2017000000000004e-9*(TC*TC*TC*TC);
       }
       else{
       rhow= 6.18545405822237e-06*(pw0 - 2.0e7) + 478.725113592986;
       drhowdT=0.0;
       }
  */
  
  return(rhow0);
}

/**  Calculate  Derivative of Density of Liquid Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_drhowdt (double /*t*/)
{
  double drhowdt;

  //Tenchev
  drhowdt=0.0;
  
  /*   //Gawin et al.
  
       pw0 = 1.0e7;  // Water pressure in range 0 - 20MPa
       
       rhow0= (4.8863e-7 - 1.6528e-9*TC0 + 1.8621e-12*(TC0*TC0) +
       2.4266e-13*(TC0*TC0*TC0) - 1.5996e-15*(TC0*TC0*TC0*TC0) +
       3.3703e-18*(TC0*TC0*TC0*TC0*TC0))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC0 + 8.7696e-3*(TC0*TC0)-
       9.2118e-5*(TC0*TC0*TC0) + 3.3534e-7*(TC0*TC0*TC0*TC0) -
       4.4034e-10*(TC0*TC0*TC0*TC0*TC0));
       
       
       if(TC <= 373.946){
       rhow= (4.8863e-7 - 1.6528e-9*TC + 1.8621e-12*(TC*TC) +
       2.4266e-13*(TC*TC*TC) - 1.5996e-15*(TC*TC*TC*TC) +
       3.3703e-18*(TC*TC*TC*TC*TC))*(pw0 - 2.0e7) +
       (1.0213e3 - 7.7377e-1*TC + 8.7696e-3*(TC*TC)-
       9.2118e-5*(TC*TC*TC) + 3.3534e-7*(TC*TC*TC*TC) -
       4.4034e-10*(TC*TC*TC*TC*TC));
       
       //Calculate Derivative of Density of Liquid Water wrt Temperature
       drhowdT=-0.7737700000000001 + (-2.0e7 + pw0)*(-1.6528000000000001e-9 + 
       3.7242e-12*(TC) + 7.2798e-13*(TC*TC) - 6.3984e-15*(TC*TC*TC) + 
       1.68515e-17*(TC*TC*TC*TC)) + 0.0175392*(TC) - 
       0.00027635400000000003*(TC*TC) + 1.3413600000000001e-6*(TC*TC*TC) - 
       2.2017000000000004e-9*(TC*TC*TC*TC);
       }
       else{
       rhow= 6.18545405822237e-06*(pw0 - 2.0e7) + 478.725113592986;
       drhowdT=0.0;
       }
  */
  return(drhowdt);
}
  
/**  Calculate Porosity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_por (double t)
{
  double tc,por,pora,porb,porc,pord;
  
  tc=t-273.15;

  //Constant
  //  por=por0;
  
  //Gawin et al.
  //  Apor=0.00017;               //  Material constant
  //  por=por0+Apor*(T-293.15);
  
  //Tenchev
  pora= (-1.0/85750000.0);   //  Constants
  porb= (27.0/1715000.0);    //
  porc= (-24.0/8575.0);      //
  pord= (389.0/343.0);       //
  
  if(tc < 100.0){
    por=por0;
  }
  if((tc >= 100.0) && (tc <= 800.0)){
    por=por0*((pora*tc*tc*tc)+(porb*tc*tc)+(porc*tc)+pord);
  }
  if(tc > 800.0){
    por=por0*3.0;
  }
  
  return(por);
}
  

  
/**  Calculate Derivative of Porosity wrt Temperature
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dpordt (double t)
{
  double dpordt,tc,pora,porb,porc,pord;

  tc=t-273.15;

  //Constant
  //  por=por0;
  
  //Gawin et al.
  //  Apor=0.00017;               //  Material constant
  //  dpordT=Apor;                //  Derivative of Porosity wrt Temperature
  
  //Tenchev
  pora=-(1.0/85750000.0);   //  Constants
  porb=(27.0/1715000.0);    //
  porc=-(24.0/8575.0);      //
  pord=(389.0/343.0);       //
  
  if(tc < 100.0){
    dpordt=0.0;                //  Derivative of Porosity wrt Temperature
  }
  if((tc >= 100.0) && (tc <= 800.0)){
    dpordt=por0*((3*pora*tc*tc)+(2*porb*tc)+porc);          //  Derivative of Porosity wrt Temperature
  }
  if(tc > 800.0){
    dpordt=0.0;              //  Derivative of Porosity wrt Temperature
  }
  
  return(dpordt);
}
  
  
/**  Calculate Volume Fraction of Free Water - from Sorption Isotherms
     @param    - 
     @retval   - 
*/
double glasgowmat::f_fracl0 (double t)
{
  double rhow0,fracl0;

  rhow0 = f_rhow0(t);

  fracl0=rhol0/rhow0;
  if(fracl0 > por0) 
    fracl0=por0;        //  CORRECTION FOR ORIGINALLY OVERSATURATED MATERIAL

  return(fracl0);
}
  
  
/**  Calculate relative humidity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_h (double rhov, double t)
{
  double pv,psat,h;

  pv = f_pv(rhov,t);
  psat = f_psat(t);

  h=pv/psat;            //  Relative Humidity

  return(h);
}
  
  
/**  Calculate factor m
     @param    - 
     @retval   - 
*/
double glasgowmat::f_m (double t)
{
  double m;
  double tc;

  tc=t-273.15;

  m=1.04-(((tc+10.0)*(tc+10.0))/((27317.5+((tc+10.0)*(tc+10.0)))));
  
  if (m > 1.0){
    m = 1.0;
  }

  return(m);
}
   
  
  
/**  Calculate factor mi
     @param    - 
     @retval   - 
*/
double glasgowmat::f_mi (double t)
{
  double m,mi;
  double tc;

  tc=t-273.15;

  m=1.04-(((tc+10.0)*(tc+10.0))/((27317.5+((tc+10.0)*(tc+10.0)))));
  mi=1.0/m;
  
  if (m > 1.0){
    m = 1.0;
    mi = 1.0;
  }

  return(mi);
}
   
  
/**  Calculate Derivative of m factor wrt Temperature
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dmdt (double t)
{
  double m,dmdt;
  double tc;

  tc=t-273.15;
  
  m=1.04-(((tc+10.0)*(tc+10.0))/((27317.5+((tc+10.0)*(tc+10.0)))));
  dmdt= -1.0*((54635.0*(10.0+tc)))/(27417.5 + 20.0*tc + tc*tc)/(27417.5 + 20.0*tc + tc*tc);
  
  
  if (m > 1.0){
    dmdt = 0.0;
  }

  return(dmdt);
}

  
/**  Calculate Derivative of mi factor wrt Temperature
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dmidt (double t)
{
  double m,dmdt,dmidt;
  double tc;

  tc=t-273.15;
  
  m=1.04-(((tc+10.0)*(tc+10.0))/((27317.5+((tc+10.0)*(tc+10.0)))));
  dmdt= -1.0*((54635.0*(10.0+tc)))/(27417.5 + 20.0*tc + tc*tc)/(27417.5 + 20.0*tc + tc*tc);
  
  dmidt = -1.0/m/m*dmdt;
  

  if (m > 1.0){
    dmidt = 0.0;
  }

  return(dmidt);
}

/**  Calculate Volume Fraction of Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_fracl (double rhov, double t)
{
  double fracl0,rhow0,fracl,h;
  double aaa,bbb,ccc,ddd;
  double mi,rhow,m,pv,psat,por;

  h = f_h(rhov,t);
  fracl0 = f_fracl0(t);
  rhow0 = f_rhow0(t);
  mi = f_mi(t);
  rhow = f_rhow(t);
  m = f_m(t);
  pv = f_pv(rhov,t);
  psat = f_psat(t);
  por = f_por(t);


  //Sorption Isotherms
  if(h <= 0.96){
    //conversion pow(a,b) -> exp(b*log(a))
    fracl=exp(mi*log((fracl0*rhow0/rhoc)*h))*rhoc/rhow;
    //fracl=(pow(((fracl0*rhow0/rhoc)*h), mi))*rhoc/rhow;
  }
  if((h > 0.96) && (h < 1.04)){

    /*  //////////////////////////////////////////////////////////////////////////////////////
    //Linear Intermediate
    double rhol96,rhol104;
    rhol96=pow((fracl0*rhow0/rhoc)*0.96,mi)*rhoc;//free water per m3 concrete (kg/m3)
    rhol104=fracl0*rhow0;//free water per m3 concrete (kg/m3)
    
    fracl=(rhol96 +((h-0.96)*((rhol104-rhol96)/0.08)))/rhow;
    //////////////////////////////////////////////////////////////////////////////////////
    */

    //////////////////////////////////////////////////////////////////////////////////////
    //Cubic Intermediate Curve Coefficients
    //conversion pow(a,b) -> exp(b*log(a))
    /*   aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.04166666666666362 + m))/(m*rhow);
	 
	 bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.042222222222221294 + m))/(m*rhow);
	 
	 ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.042824074074075444 + m))/(m*rhow);
	 
	 ddd=(3870.02879998521*fracl0*rhow0*m - 3886.9999999851375*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.043478260869569095 + m))/(m*rhow);
    */
    aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.04166666666666362 + m))/(m*rhow);
    
    bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.042222222222221294 + m))/(m*rhow);
    
    ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.042824074074075444 + m))/(m*rhow);
    
    ddd=(3870.02879998521*fracl0*rhow0*m - 3886.9999999851375*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.043478260869569095 + m))/(m*rhow);
    //fracl=ddd + (pv*(ccc*pow(psat,2) + pv*(bbb*psat + aaa*pv)))/pow(psat,3);
    fracl = aaa*(pv/psat)*(pv/psat)*(pv/psat) + bbb*(pv/psat)*(pv/psat) + ccc*(pv/psat) + ddd;
    //////////////////////////////////////////////////////////////////////////////////////

  }
  if(h >= 1.04){
    fracl=(fracl0*rhow0*(1.0+(0.12*(h-1.04))))/rhow;
  }

  //  CORRECTION FOR OVERSATURATED MATERIAL
  if(fracl > por) 
    fracl=por;
    
  return(fracl);
}

/**  Calculate Derivative of Volume Fraction of Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dfracldt (double rhov, double t)
{
  double dfracldt;
  double fracl0,rhow0,ff,ffc,h;
  double aaa,bbb,ccc,daaadt,dbbbdt,dcccdt,dddddt;
  double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3;
  double rhow,m,mi,pv,psat;
  double drhowdt,dmdt,dmidt,dpsatdt;
  double derf,derff,help;
    
  h = f_h(rhov,t);
  fracl0 = f_fracl0(t);
  rhow0 = f_rhow0(t);
  rhow = f_rhow(t);
  m = f_m(t);
  mi = f_mi(t);
  pv = f_pv(rhov,t);
  psat = f_psat(t);
  
  drhowdt = f_drhowdt(t); 
  dmdt = f_dmdt(t);
  dmidt = f_dmidt(t);
  dpsatdt = f_dpsatdt(t);

  a1 = -3887.499999985135;
  a2 = 3906.249999985064;
  a3 = 0.04166666666666362;

  b1 = 11663.249999955411;
  b2 = -11718.749999955196;
  b3 = 0.042222222222221294;

  c1 = -11645.279999955485;
  c2 = 11699.999999955271;
  c3 = 0.042824074074075444;

  d1 = 3870.02879998521;
  d2 = -3886.9999999851375;
  d3 = 0.043478260869569095;

  ff = fracl0*rhow0;
  ffc= fracl0*rhow0/rhoc;
  

  //Sorption Isotherms
  if(h <= 0.96){
    //temperature Derivative
    /*    dfracldt=-((rhoc*pow((fracl0*rhow0*pv)/(rhoc*psat),1.0/m)*drhowdt)/pow(rhow,2))+
	  (rhoc*pow((fracl0*rhow0*pv)/(rhoc*psat),1.0/m)*(-((log((fracl0*rhow0*pv)/
	  (rhoc*psat))*dmdt)/pow(m,2)) + 
	  (rhoc*psat*(-((fracl0*rhow0*pv*dpsatdt)/
	  (rhoc*pow(psat,2))) + (fracl0*rhow0*rv*rhov)/(rhoc*psat)))/
	  (fracl0*rhow0*m*pv)))/rhow;
    */
    
    //conversion pow(a,b) -> exp(b*log(a))
    /*    derf = ffc*rhov*rv/psat - ffc*pv/psat/psat*dpsatdt;
	  derff = (pow((ffc*pv/psat), mi))*(dmidt*log(ffc*pv/psat) + 
	  mi*derf/(ffc*pv/psat));
	  dfracldt = derff*rhoc/rhow - pow((ffc*pv/psat), mi)*rhoc/rhow/rhow*drhowdt;
    */
    derf = ffc*rhov*rv/psat - ffc*pv/psat/psat*dpsatdt;
    derff = (exp(mi*log(ffc*pv/psat)))*(dmidt*log(ffc*pv/psat) + 
				      mi*derf/(ffc*pv/psat));
    dfracldt = derff*rhoc/rhow - exp(mi*log(ffc*pv/psat))*rhoc/rhow/rhow*drhowdt;
  }
  if((h > 0.96) && (h < 1.04)){
    /*  //////////////////////////////////////////////////////////////////////////////////////
    //Linear Intermediate
    dfracldt=(-12.5*fracl0*rhow0*m*m*(pv*(1.0*rhow*dpsatdt + 1.0*psat*drhowdt) 
    + psat*(-0.96*psat*drhowdt 
    - 1.0*rhow*rv*rhov)) 
    - 13.0*pow(0.96,mi)*rhoc*pow((fracl0*rhow0)/rhoc,mi)
    *(-0.9615384615384615*m*m*pv*rhow*dpsatdt 
    + psat*psat*((-0.040821994520255166 + 1.0*log((fracl0*rhow0)/rhoc))*rhow*dmdt 
    + 1.0*m*m*drhowdt) + psat*(pv*((0.03925191780793766 
    -  0.9615384615384615*log((fracl0*rhow0)/rhoc))*rhow*dmdt 
    - 0.9615384615384615*m*m*drhowdt) 
    + 0.9615384615384615*m*m*rhow*rv*rhov)))/(m*m*psat*psat*rhow*rhow);
    //////////////////////////////////////////////////////////////////////////////////////
    */

    //////////////////////////////////////////////////////////////////////////////////////
    //Cubic Intermediate Curve Coefficients
    //conversion pow(a,b) -> exp(b*log(a))
    /*   aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.04166666666666362 + m))/(m*rhow);
	 
	 bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.042222222222221294 + m))/(m*rhow);
	 
	 ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*pow(0.96,mi)*rhoc*
	 pow((fracl0*rhow0)/rhoc,mi)*(0.042824074074075444 + m))/(m*rhow);
    */
    aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.04166666666666362 + m))/(m*rhow);
    
    bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.042222222222221294 + m))/(m*rhow);
    
    ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*exp(mi*log(0.96))*rhoc*
	 exp(mi*log((fracl0*rhow0)/rhoc))*(0.042824074074075444 + m))/(m*rhow);

    //temperature Derivatives of Polynomial Coefficients
    //conversion pow(a,b) -> exp(b*log(a))
    /*    daaadt = -1.0*a1*ff/rhow/rhow*drhowdt + a2*a3*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/m/rhow +
	  pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/m/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/m/rhow*dmdt +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/rhow/rhow*drhowdt) +
	  a2*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/rhow + pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/rhow/rhow*drhowdt);
	  
	  
	  dbbbdt = -1.0*b1*ff/rhow/rhow*drhowdt + b2*b3*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/m/rhow +
	  pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/m/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/m/rhow*dmdt +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/rhow/rhow*drhowdt) +
	  b2*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/rhow + pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/rhow/rhow*drhowdt);
	  
	  
	  dcccdt = -1.0*c1*ff/rhow/rhow*drhowdt + c2*c3*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/m/rhow +
	  pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/m/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/m/rhow*dmdt +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/rhow/rhow*drhowdt) +
	  c2*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/rhow + pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/rhow/rhow*drhowdt);
	  
	  
	  dddddt = -1.0*d1*ff/rhow/rhow*drhowdt + d2*d3*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/m/rhow +
	  pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/m/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/m/rhow*dmdt +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/m/rhow/rhow*drhowdt) +
	  d2*rhoc*(pow(0.96,mi)*log(0.96)*dmidt*pow(ffc,mi)/rhow + pow(0.96,mi)*pow(ffc,mi)*log(ffc)*dmidt/rhow +
	  -1.0*pow(0.96,mi)*pow(ffc,mi)/rhow/rhow*drhowdt);
    */

    daaadt = -1.0*a1*ff/rhow/rhow*drhowdt + a2*a3*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/m/rhow +
							exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)*dmidt/m/rhow +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/m/rhow*dmdt +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/rhow/rhow*drhowdt) +
      a2*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/rhow + exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)
	       *dmidt/rhow +
	       -1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/rhow/rhow*drhowdt);
    
    
    dbbbdt = -1.0*b1*ff/rhow/rhow*drhowdt + b2*b3*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/m/rhow +
							exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)*dmidt/m/rhow +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/m/rhow*dmdt +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/rhow/rhow*drhowdt) +
      b2*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/rhow + exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)
	       *dmidt/rhow +
	       -1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/rhow/rhow*drhowdt);
    
    
    dcccdt = -1.0*c1*ff/rhow/rhow*drhowdt + c2*c3*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/m/rhow +
							exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)*dmidt/m/rhow +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/m/rhow*dmdt +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/rhow/rhow*drhowdt) +
      c2*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/rhow + exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)
	       *dmidt/rhow +
	       -1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/rhow/rhow*drhowdt);
    
    
    dddddt = -1.0*d1*ff/rhow/rhow*drhowdt + d2*d3*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/m/rhow +
							exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)*dmidt/m/rhow +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/m/rhow*dmdt +
							-1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/m/rhow/rhow*drhowdt) +
      d2*rhoc*(exp(mi*log(0.96))*log(0.96)*dmidt*exp(mi*log(ffc))/rhow + exp(mi*log(0.96))*exp(mi*log(ffc))*log(ffc)
	       *dmidt/rhow +
	       -1.0*exp(mi*log(0.96))*exp(mi*log(ffc))/rhow/rhow*drhowdt);
    
    
    /*    daaadt=(3887.499999985135*fracl0*rhow0*pow(m,3)*drhowdt - 
	  3906.249999985064*pow(0.96,1.0/m)*rhoc*pow((fracl0*rhow0)/rhoc,1.0/m)*
	  (1.0*(-0.0017009164383438408 + 0.04166666666666362*log((fracl0*rhow0)/rhoc) + 
	  (0.0008446721464084541 + log((fracl0*rhow0)/rhoc))*m)*rhow*dmdt + 
	  pow(m,2)*(0.04166666666666362 + m)*drhowdt))/(pow(m,3)*pow(rhow,2));
	  
	  dbbbdt=(-11663.249999955411*fracl0*rhow0*pow(m,3)*drhowdt + 
	  11718.749999955196*pow(0.96,1.0/m)*rhoc*pow((fracl0*rhow0)/rhoc,1.0/m)*
	  ((-0.0017235953241885138 + 0.042222222222221294*log((fracl0*rhow0)/rhoc) + 
	  (0.0014002277019661302 + log((fracl0*rhow0)/rhoc))*m)*rhow*dmdt + 
	  pow(m,2)*(0.042222222222221294 + m)*drhowdt))/(pow(m,3)*pow(rhow,2));
	  
	  dcccdt=(11645.279999955485*fracl0*rhow0*pow(m,3)*drhowdt - 
	  11699.999999955271*pow(0.96,1.0/m)*rhoc*pow((fracl0*rhow0)/rhoc,1.0/m)*
	  ((-0.001748164117186909 + 0.042824074074075444*log((fracl0*rhow0)/rhoc) + 
	  (0.002002079553820279 + log((fracl0*rhow0)/rhoc))*m)*rhow*dmdt + 
	  pow(m,2)*(0.042824074074075444 + m)*drhowdt))/(pow(m,3)*pow(rhow,2));
	  
	  dddddt=(-3870.028799985211*fracl0*rhow0*pow(m,3)*drhowdt + 
	  3886.9999999851375*pow(0.96,1.0/m)*rhoc*pow((fracl0*rhow0)/rhoc,1.0/m)*
	  ((-0.0017748693269677752 + 0.043478260869569095*log((fracl0*rhow0)/rhoc) + 
	  (0.0026562663493139294 + log((fracl0*rhow0)/rhoc))*m)*rhow*dmdt + 
	  pow(m,2)*(0.043478260869569095 + m)*drhowdt))/(pow(m,3)*pow(rhow,2));
    */
    //temperature Derivative
    /*  dfracldt=(pow(psat,4)*dddddt - 3*aaa*pow(pv,3)*dpsatdt + psat*pow(pv,2)*
	(pv*daaadt - 2*bbb*dpsatdt + 3*aaa*rv*rhov) + pow(psat,2)*pv*
	(pv*dbbbdt - ccc*dpsatdt + 2*bbb*rv*rhov) + pow(psat,3)*(pv*
	dcccdt + ccc*rv*rhov))/pow(psat,4);
    */

    help = (rhov*rv*psat - pv*dpsatdt)/psat/psat;
    dfracldt = daaadt*(pv/psat)*(pv/psat)*(pv/psat) + dbbbdt*(pv/psat)*(pv/psat) + dcccdt*(pv/psat) + dddddt;
    dfracldt = dfracldt + 3.0*aaa*(pv/psat)*(pv/psat)*help + 2.0*bbb*(pv/psat)*help + ccc*help; 
    //////////////////////////////////////////////////////////////////////////////////////
    
  }
  if(h >= 1.04){
    //temperature Derivative
    /*   dfracldt=-((fracl0*rhow0*(1.0 + 0.12*(-1.04 + pv/psat))*drhowdt)/pow(rhow,2))
	 + (0.12*fracl0*rhow0*(-((pv*dpsatdt)/pow(psat,2)) + rv*rhov/psat))/rhow;
    */
    help = 0.12*(rhov*rv*psat - pv*dpsatdt)/psat/psat;
    dfracldt = fracl0*rhow0*(help*rhow - (1.0+(0.12*(h+1.04)))*drhowdt)/rhow/rhow;
  }
  
  /*  if (Tp->time >= 305.0){
      fprintf (Outt,"\nh        = %e",h);
      fprintf (Outt,"\ndfracldt = %e",dfracldt);
      }
      //dfracldt = -1.437836e-03;
      
      if (fabs(dfracldt) > 9.0e-03)
      dfracldt = -1.437836e-03;
  */

  return(dfracldt);
}

/**  Calculate Vapour Content Derivative
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dfracldrhov (double rhov, double t)
{
  double dfracldrhov;
  double fracl0,rhow0,h;
  double rhow,m,pv,psat;
  double aaa,bbb,ccc;
  
  h = f_h(rhov,t);
  fracl0 = f_fracl0(t);
  rhow0 = f_rhow0(t);
  rhow = f_rhow(t);
  m = f_m(t);
  pv = f_pv(rhov,t);
  psat = f_psat(t);

  if(h <= 0.96){
    //Vapour Content Derivative
    //conversion pow(a,b) -> exp(b*log(a))
    //dfracldrhov=(fracl0*rhow0*pow((fracl0*rhow0*pv)/(rhoc*psat),-1.0 + 1.0/m)*rv*t)/(m*psat*rhow);
    dfracldrhov=(fracl0*rhow0*exp((-1.0 + 1.0/m)*log((fracl0*rhow0*pv)/(rhoc*psat)))*rv*t)/(m*psat*rhow);
  }
  if((h > 0.96) && (h < 1.04)){

    /*   //////////////////////////////////////////////////////////////////////////////////////
    //Linear Intermediate
    double mi;
    mi = f_mi(t);
    dfracldrhov=(12.5*fracl0*rhow0*rv*t - 12.5*pow(0.96,mi)*rhoc*pow((fracl0*rhow0)/rhoc,mi)*rv*t)/(psat*rhow);
    //////////////////////////////////////////////////////////////////////////////////////
    */
    
    //////////////////////////////////////////////////////////////////////////////////////
    //Cubic Intermediate Curve Coefficients
    //conversion pow(a,b) -> exp(b*log(a))
    /*    aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*pow(0.96,1.0/m)*rhoc*
	  pow((fracl0*rhow0)/rhoc,1.0/m)*(0.04166666666666362 + m))/(m*rhow);
	  
	  bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*pow(0.96,1.0/m)*rhoc*
	  pow((fracl0*rhow0)/rhoc,1.0/m)*(0.042222222222221294 + m))/(m*rhow);
	  
	  ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*pow(0.96,1.0/m)*rhoc*
	  pow((fracl0*rhow0)/rhoc,1.0/m)*(0.042824074074075444 + m))/(m*rhow);
    */
    aaa=(-3887.499999985135*fracl0*rhow0*m + 3906.249999985064*exp((1.0/m)*log(0.96))*rhoc*
	 exp((1.0/m)*log((fracl0*rhow0)/rhoc))*(0.04166666666666362 + m))/(m*rhow);
    
    bbb=(11663.249999955411*fracl0*rhow0*m - 11718.749999955196*exp((1.0/m)*log(0.96))*rhoc*
	 exp((1.0/m)*log((fracl0*rhow0)/rhoc))*(0.042222222222221294 + m))/(m*rhow);
    
    ccc=(-11645.279999955485*fracl0*rhow0*m + 11699.999999955271*exp((1.0/m)*log(0.96))*rhoc*
	 exp((1.0/m)*log((fracl0*rhow0)/rhoc))*(0.042824074074075444 + m))/(m*rhow);
    
    //Vapour Content Derivative
    //dfracldrhov=((ccc*pow(psat,2) + pv*(2*bbb*psat + 3*aaa*pv))*rv*t)/pow(psat,3);
    dfracldrhov = 3.0*aaa*(pv/psat)*(pv/psat/psat)*(rv*t) + 2.0*bbb*(pv/psat/psat)*(rv*t) + ccc*(rv*t/psat);
    //////////////////////////////////////////////////////////////////////////////////////

  }
  if(h >= 1.04){
    //Vapour Content Derivative
    dfracldrhov=(0.12*fracl0*rhow0*rv*t)/(psat*rhow);
  }
  return(dfracldrhov);
}


/**  Calculate Mass of Dehydrated Water Released per m3
     @param    - 
     @retval   - 
*/
double glasgowmat::f_mcbwrel (double t)
{
  double tc,mcbwrel,mcbw0;
    
  tc=t-273.15;

  //Constant
  //  mcbw0=0.0;
  //  mcbwrel=0.0;
  //  dmcbwreldt=0.0;      
  
  //Tenchev 1
  mcbw0=0.09*rhoc;
  if(tc <= 200.0){
    mcbwrel=0.0;
  }
  if((tc > 200.0) && (tc <= 300.0)){
    mcbwrel=rhoc*(7.0e-4*(tc-200.0));
  }
  if((tc > 300.0) && (tc <= 800.0)){
    mcbwrel=rhoc*(0.4e-4*(tc-300.0)+0.07);
  }
  if(tc > 800.0){
    mcbwrel=rhoc*0.09;
  }
  
  //Tenchev 2
  //  mcbw0=0.24*rhoc;
  //  if(tc < 100.0){
  //    mcbwrel=0.0;
  //    dmcbwreldT=0.0;
  //  }
  //  if(tc >= 100.0 && tc <= 700.0){
  //    mcbwrel=0.04*(tc/100.0-1.0)*rhoc;
  //    dmcbwreldT=4.0e-4*rhoc;
  //  }
  //  if(tc > 700.0){
  //    mcbwrel=0.24*rhoc;
  //    dmcbwreldT=0.0;
  //  }
  
  //Ichikawa
  //  hdeg=1.0;                   // Degree of Hydration
  //  mcbw0=0.23*hdeg*rhoc;
  //  if(tc <= 105.0){
  //    mcbwrel=0.0;
  //    dmcbwreldT=0.0;
  //  }
  //  if(tc > 105.0 && tc <= 850.0){
  //    mcbwrel=mcbw0*(tc-105.0)/745.0;
  //    dmcbwreldT=mcbw0/745.0;
  //  }
  //  if(tc > 850.0){
  //    mcbwrel=mcbw0;
  //    dmcbwreldT=0.0;
  //  }
  
  return(mcbwrel);
}

/**  Calculate Derivative of Mass of Dehydrated Water Released per m3
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dmcbwreldt (double t)
{
  double tc,dmcbwreldt,mcbw0;
    
  tc=t-273.15;

  //Constant
  //  mcbw0=0.0;
  //  dmcbwreldt=0.0;      
  
  //Tenchev 1
  mcbw0=0.09*rhoc;
  if(tc <= 200.0){
    dmcbwreldt=0.0;
  }
  if((tc > 200.0) && (tc <= 300.0)){
    dmcbwreldt=rhoc*7.0e-4;
  }
  if((tc > 300.0) && (tc <= 800.0)){
    dmcbwreldt=rhoc*0.4e-4;
  }
  if(tc > 800.0){
    dmcbwreldt=0.0;
  }
  
  //Tenchev 2
  //  mcbw0=0.24*rhoc;
  //  if(tc < 100.0){
  //    mcbwrel=0.0;
  //    dmcbwreldt=0.0;
  //  }
  //  if(tc >= 100.0 && tc <= 700.0){
  //    mcbwrel=0.04*(tc/100.0-1.0)*rhoc;
  //    dmcbwreldt=4.0e-4*rhoc;
  //  }
  //  if(tc > 700.0){
  //    mcbwrel=0.24*rhoc;
  //    dmcbwreldt=0.0;
  //  }
  
  //Ichikawa
  //  hdeg=1.0;                   // Degree of Hydration
  //  mcbw0=0.23*hdeg*rhoc;
  //  if(tc <= 105.0){
  //    mcbwrel=0.0;
  //    dmcbwreldt=0.0;
  //  }
  //  if(tc > 105.0 && tc <= 850.0){
  //    mcbwrel=mcbw0*(tc-105.0)/745.0;
  //    dmcbwreldt=mcbw0/745.0;
  //  }
  //  if(tc > 850.0){
  //    mcbwrel=mcbw0;
  //    dmcbwreldt=0.0;
  //  }
  
  return(dmcbwreldt);
}
  
/**  Calculate Volume Fraction of Dehydrated Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_fracd (double t)
{
  double fracd,mcbwrel,rhow;

  mcbwrel = f_mcbwrel(t);
  rhow = f_rhow(t);
  fracd=mcbwrel/rhow;
  
  return(fracd);
}
 
/**  Calculate Derivative of Volume Fraction of Dehydrated Water wrt Temperature
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dfracddt(double t)
{
  double dfracddt,dmcbwreldt,rhow;

  dmcbwreldt = f_dmcbwreldt(t);
  rhow = f_rhow(t);

  dfracddt=dmcbwreldt/rhow;
  
  return(dfracddt);
}
  
/**  Calculate Volume Fraction of Gas Phase 
     @param    - 
     @retval   - 
*/
double glasgowmat::f_fracg(double rhov, double t)
{
  double fracg,por,fracl;

  por = f_por(t);
  fracl = f_fracl(rhov,t);

  fracg=por-fracl;

  return(fracg);
}
  
/**  Calculate Saturation 
     @param    - 
     @retval   - 
*/
double glasgowmat::f_s(double rhov, double t)
{
  double s,fracl,por;  
 
  por = f_por(t);
  fracl = f_fracl(rhov,t);
  
  s=fracl/por;

  return(s);
}

/**  Calculate pl
     @param    - 
     @retval   - 
*/
double glasgowmat::f_pl(double rhov, double pg, double t)
{
  double pl,s,pc;
  
  s = f_s(rhov,t);
  pc = f_pc(rhov,t);

  switch (model){
  case 1:{
    pl=pg;
    break;
  }
  case 2:{
    if(s <= sssp){
      pl = 0.0;
    }
    else{
      pl = pg - pc;
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }  
  return(pl);
}
  
/**  Calculate pc
     @param    - 
     @retval   - 
*/
double glasgowmat::f_pc (double rhov, double t)
{
  double pc,s,h,rhow;

  s = f_s(rhov,t);
  h = f_h(rhov,t);
  rhow = f_rhow(t);

  switch (model){
  case 1:{
    pc=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      //pc=0.0;//tady pokus??!!
      pc=-1.0*rv*t*rhow*log(h);
    }
    else{
      pc=-1.0*rv*t*rhow*log(h);
      if(pc < 0.0) { //changed by Tomas Krejci (added braces for "if")
	pc=0.0;              //  CORRECTION FOR RH>100%
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return(pc);
}
  
/**  Calculate dpcdt
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dpcdt (double rhov, double t)
{
  double dpcdt,pc,s,h,rhow,drhowdt,psat,dpsatdt;

  s = f_s(rhov,t);
  h = f_h(rhov,t);
  rhow = f_rhow(t);
  drhowdt = f_drhowdt(t);
  psat = f_psat(t);
  dpsatdt = f_dpsatdt(t);

  switch (model){
  case 1:{
    dpcdt=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      dpcdt=0.0;
    }
    else{
      pc=-1.0*rv*t*rhow*log(h);
      if(pc < 0.0) { //changed by Tomas Krejci (added braces for "if")
	dpcdt=0.0;
      }
      else{
	dpcdt=rv*rhow*((t*dpsatdt/psat)-((1.0+t*drhowdt/rhow)*log(rhov*rv*t/psat))-1.0);
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(dpcdt);
}
  
/**  Calculate dpcdrhov
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dpcdrhov (double rhov, double t)
{
  double dpcdrhov,pc,s,h,rhow;

  s = f_s(rhov,t);
  h = f_h(rhov,t);
  rhow = f_rhow(t);
  
  switch (model){
  case 1:{
    dpcdrhov=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      dpcdrhov=0.0;
    }
    else{
      pc=-1.0*rv*t*rhow*log(h);
      if(pc < 0.0) { //changed by Tomas Krejci (added braces for "if")
	dpcdrhov=0.0;
      }
      else{
	dpcdrhov=-rv*t*rhow/rhov;
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(dpcdrhov);
}
  
/**  Calculate sln
     @param    - 
     @retval   - 
*/
double glasgowmat::f_sln (double rhov, double t)
{
  double sln,s;

  s = f_s(rhov,t);

  switch (model){
  case 1:{
    sln=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      sln=0.0;
    }
    else{
      sln=(s-sssp)/(1.0-sssp);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(sln);
}
  
/**  Calculate ppore
     @param    - 
     @retval   - 
*/
double glasgowmat::f_ppore (double rhov, double pg, double t, double pginf)
{
  double ppore;
  double sln,pc,s,h,rhow;

  s = f_s(rhov,t);
  h = f_h(rhov,t);
  rhow = f_rhow(t);

  switch (model){
  case 1:{
    ppore=pg-pginf;
    pc=0.0;
    sln=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      pc=0.0;
      sln=0.0;
      ppore=pg-pginf;
    }
    else{
      pc=-1.0*rv*t*rhow*log(h);
      if(pc < 0.0) { //changed by Tomas Krejci (added braces for "if")
	pc=0.0;              //  CORRECTION FOR RH>100%
      }
      else{}
      sln=(s-sssp)/(1.0-sssp);
      ppore=pg-(sln*pc)-pginf;
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(ppore);
}
  
/**  Calculate sb
     @param    - 
     @retval   - 
*/
double glasgowmat::f_sb (double rhov, double t)
{
  double sb,s;

  s = f_s(rhov,t);

  switch (model){
  case 1:{
    sb=0.0;
    break;
  }
  case 2:{
    if(s <= sssp){
      sb=s;                             //  Physically bound water saturation
    }
    else{
      sb=sssp;                          //  Physically bound water saturation
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return(sb);
}

/**  Calculate Dry Air Pressure
     @param    - 
     @retval   - 
*/
double glasgowmat::f_pa (double rhov, double pg, double t)
{
  double pa,pv;

  pv = f_pv(rhov,t);

  pa=pg-pv;
  
  return(pa);
}
 
/**  Calculate Air Content
     @param    - 
     @retval   - 
*/
double glasgowmat::f_rhoa (double rhov, double pg, double t)
{
  double rhoa,pa;

  pa = f_pa(rhov,pg,t);

  rhoa=pa/ra/t;

  return(rhoa);
}
  
/**  Calculate Gas Content
     @param    - 
     @retval   - 
*/
double glasgowmat::f_rhog (double rhov, double pg, double t)
{
  double rhog,rhoa;

  rhoa = f_rhoa(rhov,pg,t);

  rhog=rhoa+rhov;

  //fprintf (Outt,"\nrhog = %e",rhog);

  return(rhog);
}
  
/**  Calculate dav - Calculate Gas Diffusivity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_dav (double /*rhov*/, double pg, double t)
{
  double dif,delta,tor,dav;
    
  //Constant
  //  dav=2.42e-5;                                    //   From CENGEL BOOK
    
  //Tenchev (Cengel)
  //Internal
  //conversion pow(a,b) -> exp(b*log(a))
  dif=1.87*exp(2.072*log(t))/pg*1.0e-5;
  //dif=1.87*(pow(t,2.072)/pg)*1.0e-5;
  delta=0.5;                                        // Constrictivity
  tor=3.0;                                          // Tortuosity
  dav=dif*delta/tor/tor;

  //fprintf(Outt,"\n------------------\n");
  //fprintf(Outt,"deff = %e\n",dav);
  //fprintf(Outt,"------------------\n");
  
  return(dav);
}

/**  Calculate DavEx - Calculate Gas Diffusivity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_davex (double pg, double tinf)
{
  double davex,tor,delta,difex;
    
  //Constant
  //  dav=2.42e-5;                                    //   From CENGEL BOOK
  //  davex=dav;
    
  //Tenchev (Cengel)
  //Internal
  delta=0.5;                                        // Constrictivity
  tor=3.0;                                          // Tortuosity

  //CENGEL BOOK
  //conversion pow(a,b) -> exp(b*log(a))
  difex=1.87*exp(2.072*log(tinf))/pg*1.0e-5;
  //difex=1.87*(pow(tinf,2.072)/pg)*1.0e-5;
  davex=difex*delta/tor/tor;
  
  return(davex);
}

/** Calculate Diffusivity of Physically Bound Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_db(double rhov, double t)
{
  double db,db0,dbtref,s;

  s = f_s(rhov,t);

  //Gawin et al.
  db0=1.57e-11;
  dbtref=295.0;
  if(s <= sssp)
    db=db0*exp(-2.08*(s/sssp)*(t/dbtref));
  else
    db=0.0;
  
  return(db);
}
  
  
/**  Calculate Intrinsic Permeability 
     @param    - 
     @retval   - 
*/
double glasgowmat::f_kk (double t)
{
  double kk,por;
 
  //Constant
  //  KK=k0
  
  //Gawin et al.
  //  TK0=293.15;         //  Damage not included at present
  //  AKT=0.005;          //  Constants AKT, AKP and AKD should be changed for damage
  //  AKP=0.368;
  //  AKD=
  //  P0=0.10e6;          //  Atmospheric pressure
  //  KK=k0*pow(10,(AKT*(T-TK0)))*pow((pg/P0),AKP)   //*pow(10,(AKD*D))
  
  
  //Tenchev
  por = f_por(t);

  //conversion pow(a,b) -> exp(b*log(a))
  kk=exp((2.0/3.0)*log(por/por0))*k0;
  //kk=pow((por/por0),(2.0/3.0))*k0;
  //tady??!!
  kk = k0;//pokus??!!
  

  return(kk);
}
  
/**  Calculate Calculate Relative Permeability
     @param    - 
     @retval   - 
*/
double glasgowmat::f_kg (double rhov, double t)
{
  double kg,kb,km,s;

  //Calculate Relative Permeabilities========================================================***************OPTIONS***************
  //  Sir=0.0;                  //  Irreducible saturation
  //  Scr=1.0;                  //  Critical saturation for no gas flow
  //  Se=(S-Sir)/(1-Sir);       //  Effective Saturation
  
  //Brooks and Corey
  //  PSDI=1.0;                 //  Pore size distribution index
  //  Kg=pow((1.0-Se),2)*pow(1.0-Se,((2.0+PSDI)/PSDI));
  //  Kl=pow(Se,((2.0+3.0*PSDI)/PSDI));
  
  //Gawin et al.
  //  KgAg=2.0;                  //  Constant between 1 and 3
  //  KlAw=2.0;                  //  Constant between 1 and 3
  //  Kg=pow(1.0-(S/Scr),KgAg);
  //  Kl=pow(Se,KlAw);
  
  //Alternative Gawin et al.
  //  KlBw=6.0;                        //  Constant equal to 6 or 16
  //  Kl=pow(pow(1.0+((1.0-h)/0.25),KlBw),-1.0)*pow(S,KlAw);
  
  s = f_s(rhov,t);
  
  //Baroghel-Bouny et al.
  kb=2.2748;               //  For ordinary concrete (See Baroghel-Bouny et al.)
  km=1.0/kb;
  //conversion pow(a,b) -> exp(b*log(a))
  kg=(sqrt(1.0-s))*exp((2.0*km)*log(1.0-exp((1.0/km)*log(s))));
  //kg=(sqrt(1.0-s))*pow((1.0-pow(s,(1.0/km))),(2.0*km));
  
  //Tenchev
  //  Kg=1.0-S;
  //  Kl=0.01;
  
  return(kg);
}
  
/**  Calculate Calculate Relative Permeability
     @param    - 
     @retval   - 
*/
double glasgowmat::f_kl (double rhov, double t)
{
  double kl,km,kb,s;

  //Calculate Relative Permeabilities========================================================***************OPTIONS***************
  //  Sir=0.0;                  //  Irreducible saturation
  //  Scr=1.0;                  //  Critical saturation for no gas flow
  //  Se=(S-Sir)/(1-Sir);       //  Effective Saturation
  
  //Brooks and Corey
  //  PSDI=1.0;                 //  Pore size distribution index
  //  Kg=pow((1.0-Se),2)*pow(1.0-Se,((2.0+PSDI)/PSDI));
  //  Kl=pow(Se,((2.0+3.0*PSDI)/PSDI));
  
  //Gawin et al.
  //  KgAg=2.0;                  //  Constant between 1 and 3
  //  KlAw=2.0;                  //  Constant between 1 and 3
  //  Kg=pow(1.0-(S/Scr),KgAg);
  //  Kl=pow(Se,KlAw);
  
  //Alternative Gawin et al.
  //  KlBw=6.0;                        //  Constant equal to 6 or 16
  //  Kl=pow(pow(1.0+((1.0-h)/0.25),KlBw),-1.0)*pow(S,KlAw);
  
  s = f_s(rhov,t);
  
  //Baroghel-Bouny et al.
  kb=2.2748;               //  For ordinary concrete (See Baroghel-Bouny et al.)
  km=1.0/kb;
  //conversion pow(a,b) -> exp(b*log(a))
  kl=sqrt(s)*exp(2.0*log(1.0-exp(km*log(1.0-exp((1.0/km)*log(s))))));
  //kl=sqrt(s)*pow((1.0-pow((1.0-pow(s,(1.0/km))),km)),2.0);
  
  //Tenchev
  //  Kg=1.0-S;
  //  Kl=0.01;
  
  return(kl);
}
  
/**  Calculate Dynamic Viscosity of Liquid Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_mul (double t)
{
  double mul;

  //Constant
  //  mul=1.002e-3
  
  //  Gawin et al.  
  if(t <= 647.3)
    //conversion pow(a,b) -> exp(b*log(a))
    //mul=0.6612*(pow(t-229.0,-1.562));
    mul=0.6612*exp((-1.562)*log(t-229.0));
  else
    mul=0.000043;
  
  return(mul);
}
  
  
/**  Calculate Dynamic Viscosity of Water Vapour
     @param    - 
     @retval   - 
*/
double glasgowmat::f_muv (double t)
{
  double muv,t0muv,amuv,muv0;

  //Constant
  //  muv=0.973e-5
  
  //Gawin et al.
  if(t <= 647.3){
    t0muv=273.15;                 // I think
    muv0=8.85e-6;
    amuv=3.53e-8;
    muv=muv0+amuv*(t-t0muv);      // DOESN'T MATCH CENGEL'S DATA
  }
  else
    muv=0.00004313;
  
  return(muv);
}
  
  
/**  Calculate Dynamic Viscosity of Dry Air
     @param    - 
     @retval   - 
*/
double glasgowmat::f_mua (double t)
{
  double mua,t0mua,amua,bmua,mua0;

  //Constant
  //  mua=1.825e-5;
  
  //Gawin et al.
  t0mua=273.15;                     // I think
  mua0=17.17e-6;
  amua=4.73e-8;
  bmua=2.22e-11;
  //conversion pow(a,b) -> exp(b*log(a))
  //mua=mua0+(amua*(t-t0mua))-(bmua*pow(t-t0mua,2.0));
  mua=mua0+(amua*(t-t0mua))-(bmua*exp(2.0*log(t-t0mua)));
  
  return(mua);
}
  
  
/**  Calculate Dynamic Viscosity of Gas
     @param    - 
     @retval   - 
*/
double glasgowmat::f_mug (double rhov, double pg, double t)
{
  double mug,mua,muv,rhoa;

  mua = f_mua(t);
  muv = f_muv(t);
  rhoa = f_rhoa(rhov,pg,t);

  //Tenchev
  if((rhoa == 0.0) && (rhov == 0.0))
    mug=0.0;
  else
    mug=((rhoa*mua)+(rhov*muv))/(rhoa+rhov);
  
  //Gawin et al.
  //  mug=muv+((mua-muv)*pow(pa/pg,0.608));
  
  return(mug);
}
  
   
/**  Calculate Specific Heat of Liquid Water
     @param    - 
     @retval   - 
*/
double glasgowmat::f_cl (double t)
{
  double cl;

  //Constant
  //  cl=4182.0;
  
  //Tenchev (Cengel)
  if(t <= 647.3)
    //FOR SATURATED STATE - ONLY VALID UP TO CRITICAL POINT
    //conversion pow(a,b) -> exp(b*log(a))
    //cl=((2.4768*t)+3368.2)+pow(((1.08542631988638*t)/513.15),31.4447657616636);
    cl=((2.4768*t)+3368.2)+exp(31.4447657616636*log((1.08542631988638*t)/513.15));
  else
    cl=24515.0;
  
  return(cl);
}
  
  
/**  Calculate Specific Heat of Water Vapour
     @param    - 
     @retval   - 
*/
double glasgowmat::f_cv (double t)
{
  double cv;

  //Constant
  //cv=1867.0;
  
  //Tenchev (Cengel)
  if(t <= 647.3)
    //FOR SATURATED STATE - ONLY VALID UP TO CRITICAL POINT
    //conversion pow(a,b) -> exp(b*log(a))
    //cv=((7.1399*t)-443.0)+pow(((1.13771502228162*t)/513.15),29.4435287521143);
    cv=((7.1399*t)-443.0)+exp(29.4435287521143*log((1.13771502228162*t)/513.15));
  else
    cv=45821.04;
  
  return(cv);
}
  
  
/**  Calculate Specific Heat of Dry Air
     @param    - 
     @retval   - 
*/
double glasgowmat::f_ca (double  t)
{
  double ca;

  //Constant
  //  Ca=1007.0;
  
  ca = -0.00000009849367018147350*t*t*t + 0.0003564362577698610*t*t - 0.1216179239877570*t +1012.502552163240;
  
  return(ca);
}
  
   
/**  Calculate Specific Heat of Solid Skeleton
     @param    - 
     @retval   - 
*/
double glasgowmat::f_cs (double t)
{
  double cs,tc;

  tc=t-273.15;

  //Constant
  //  cs=1071.562;
  
  //Tenchev
  cs=900.0+80.0*(tc/120.0)-(4.0*tc*tc/120.0/120.0);
    
  //Gawin et al.
  
  return(cs);
}
  
  
/**  Calculate Thermal Conductivity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_keff (double t)
{
  double keff,tc;

  tc=t-273.15;

  //Constant
  //  keff=1.642218;
  
  //Tenchev
  keff=2.0-(0.24*(tc/120.0))+(0.012*tc*tc/120.0/120.0);
  
  //Averaged from Eurocode - Anderberg
  //  keff=1.68-(0.19055*(TC/100.0))+(0.0082*TC*TC/100.0/100.0);
  
  return(keff);
}
  
  
/**  Calculate Volume Fraction of Solid Skeleton
     @param    - 
     @retval   - 
*/
double glasgowmat::f_fracs (double rhov, double t)
{
  double fracs,fracg,fracl;

  fracg = f_fracg(rhov,t);
  fracl = f_fracl(rhov,t);

  fracs=1.0-fracg-fracl;
  
  return(fracs);
}
  
  
/**  Calculate Effective Heat Capacity
     @param    - 
     @retval   - 
*/
double glasgowmat::f_crho (double rhov, double pg, double t)
{
  double crho,fracs,cs,fracl,rhow,cl,cv,fracg,rhoa,ca;

  fracs = f_fracs(rhov,t);
  cs = f_cs(t);
  fracl = f_fracl(rhov,t);
  rhow = f_rhow(t);
  cl = f_cl(t);
  fracg = f_fracg(rhov,t);
  cv = f_cv(t);
  rhoa = f_rhoa(rhov,pg,t);
  ca = f_ca(t);

  crho=(rhos*fracs*cs)+(fracl*rhow*cl)+(fracg*rhov*cv)+(fracg*rhoa*ca);
  
  return(crho);
}
  
  
/**  Calculate Specific Heat of Evaporation
     @param    - 
     @retval   - 
*/
double glasgowmat::f_le (double t)
{
  double le;

  //Constant
  /*  if(t <= 647.3)
      le=2454000.0;
      else
      le=0.0;
  */

  //Gawin et al.
  if(t <= 647.3)
    //conversion pow(a,b) -> exp(b*log(a))
    //le=2.672e5*pow((647.3-t),0.38);
    le=2.672e5*exp(0.38*log(647.3-t));
  else
    le=0.0;

  return(le);
}
  
   
/**  Calculate Specific Heat of Dehydration
     @param    - 
     @retval   - 
*/
double glasgowmat::f_ld ()
{
  double ld;

  //Tenchev
  ld=2400.0e3;
  
  return(ld);
}
  
   
/**  Calculate Radiative Heat Transfer Coefficient
     @param    - 
     @retval   - 
*/
double glasgowmat::f_hrad(double t, double tinf)
{
  double hrad;

  //Constant
  //hrad=3.428182994;
  
  //Tenchev
  hrad=emmi*stef*((t*t)+(tinf*tinf))*(t+tinf);

  return(hrad);
}
  
/**  Calculate Combined Heat Transfer Coefficient
     @param    - 
     @retval   - 
*/
double glasgowmat::f_hqr(double t, double tinf)
{
  double hrad,hqr;

  hrad=f_hrad(t,tinf);
  
  //Constant
  //hqr=1.0e6;
  
  hqr=hq+hrad;
  
  return(hqr);
}
  
  
/**  Calculate Water Vapour Transfer Coefficient
     @param    - 
     @retval   - 
*/
double glasgowmat::f_beta (double pg, double tinf)
{
  double beta,davex;

  davex = f_davex(pg,tinf);

  //  Should this be DavEx or DifEx I.e. no tortuosity or constrictivity effects
  //conversion pow(a,b) -> exp(b*log(a))
  //beta=(hq/crhoair)*pow((davex/alph),(2.0/3.0));
  beta=(hq/crhoair)*exp((2.0/3.0)*log(davex/alph));
  
  return(beta);
}
  

/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

double glasgowmat::ktt (double t,double /*pg*/,double rhov)
{
  double ktt,keff;

  keff = f_keff(t);

  switch (model){
  case 1:{
    ktt = keff;
    break;
  }
  case 2:{
    double klt,fracl,rhow,sb,db,s,por,dfracldt,dpordt,kk,kl,dpcdt,mul,le;
    
    fracl = f_fracl(rhov,t);
    rhow = f_rhow(t);
    sb = f_sb(rhov,t);
    db = f_db(rhov,t);
    s = f_s(rhov,t);
    por = f_por(t);
    dfracldt = f_dfracldt(rhov,t);
    dpordt = f_dpordt(t);
    kk = f_kk(t);
    kl = f_kl(rhov,t);
    dpcdt = f_dpcdt(rhov,t);
    mul = f_mul(t);
    le = f_le(t);

    klt=fracl*rhow*(((sb*db/s/por)*(dfracldt-(fracl*dpordt/por)))-(1.0-(sb/s))*(kk*kl*dpcdt/mul));
    ktt = keff-(le*klt);
    
    /*    fprintf (Outt,"\n");
	  fprintf (Outt,"\nt = %e",t);
	  fprintf (Outt,"\npg = %e",pg);
	  fprintf (Outt,"\nrhov = %e",rhov);
	  fprintf (Outt,"\nklt = %e",klt);
	  fprintf (Outt,"\nkeff = %e",keff);
	  fprintf (Outt,"\nktt = %e",ktt);
    */
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return (ktt);
}

double glasgowmat::ktp (double t,double /*pg*/,double rhov)
{
  double ktp,klp;
  double sb,s,fracl,rhow,kk,kl,mul,le;
  
  sb = f_sb(rhov,t);
  s = f_s(rhov,t);
  fracl = f_fracl(rhov,t);
  rhow = f_rhow(t);
  kk = f_kk(t);
  kl = f_kl(rhov,t);
  mul = f_mul(t);
  le = f_le(t);

  switch (model){
  case 1:{
    ktp= -le*fracl*rhow*kk*kl/mul;
    break;
  }
  case 2:{
    klp=(1.0-(sb/s))*fracl*rhow*kk*kl/mul;
    ktp= -le*klp;
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return (ktp);
}

double glasgowmat::ktv (double t,double /*pg*/,double rhov)
{
  double ktv;
  double klv,db,sb,s,fracl,rhow,kk,kl,mul,le,dfracldrhov,dpcdrhov,por;
  
  sb = f_sb(rhov,t);
  db = f_db(rhov,t);
  s = f_s(rhov,t);
  fracl = f_fracl(rhov,t);
  dfracldrhov = f_dfracldrhov(rhov,t);
  dpcdrhov = f_dpcdrhov(rhov,t);
  rhow = f_rhow(t);
  kk = f_kk(t);
  por = f_por(t);
  kl = f_kl(rhov,t);
  mul = f_mul(t);
  le = f_le(t);

  switch (model){
  case 1:{
    ktv=0.0;
    break;
  }
  case 2:{
    klv=fracl*rhow*((sb*db*dfracldrhov/s/por)-((1.0-(sb/s))*kk*kl*dpcdrhov/mul));
    ktv= -le*klv;
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return (ktv);
}

double glasgowmat::kat (double t,double pg,double rhov)
{
  double kat,dav,fracg,rhog;
   
  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhog = f_rhog(rhov,pg,t);
 
  /*  fprintf (Outt,"\n");
      fprintf (Outt,"\nt = %e",t);
      fprintf (Outt,"\npg = %e",pg);
      fprintf (Outt,"\nrhov = %e",rhov);
      fprintf (Outt,"\ndav = %e",dav);
      fprintf (Outt,"\nfracg = %e",fracg);
      fprintf (Outt,"\nrhog = %e",rhog);
  */
  kat= -1.0*dav*fracg*rhov*pg/rhog/ra/t/t;
  
  return (kat);
}

double glasgowmat::kap (double t,double pg,double rhov)
{
  double kap;
  double dav,fracg,rhoa,rhog,kk,kg,mug;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);
  kk = f_kk(t);
  kg = f_kg(rhov,t);
  mug = f_mug(rhov,pg,t);

  kap= (kk*kg*fracg*rhoa/mug)+(dav*fracg*rhov/rhog/ra/t);

  return (kap);
}

double glasgowmat::kav (double t,double pg,double rhov)
{
  double kav;
  double dav,fracg,rhoa,rhog;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);
  
  kav= -1.0*(rhoa+(rv*rhov/ra))*dav*fracg/rhog;
  
  return (kav);
}

double glasgowmat::kmt (double t,double pg,double rhov)
{
  double kmt,klt,kvt;
  double dav,fracg,rhoa,rhog;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);

  switch (model){
  case 1:{
    kvt=fracg*dav*rhov*pg/rhog/ra/t/t;
    klt=0.0;
    kmt=kvt+klt;
    break;
  }
  case 2:{
    double db,sb,s,fracl,rhow,kk,kl,mul,le,dfracldt,dpcdt,por,dpordt;

    sb = f_sb(rhov,t);
    db = f_db(rhov,t);
    s = f_s(rhov,t);
    fracl = f_fracl(rhov,t);
    dfracldt = f_dfracldt(rhov,t);
    dpcdt = f_dpcdt(rhov,t);
    rhow = f_rhow(t);
    kk = f_kk(t);
    por = f_por(t);
    dpordt = f_dpordt(t);
    kl = f_kl(rhov,t);
    mul = f_mul(t);
    le = f_le(t);

    kvt=fracg*dav*rhov*pg/rhog/ra/t/t;
    klt=fracl*rhow*(((sb*db/s/por)*(dfracldt-(fracl*dpordt/por)))-(1.0-(sb/s))*(kk*kl*dpcdt/mul));
    kmt=kvt+klt;
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return (kmt);
}

double glasgowmat::kmp (double t,double pg,double rhov)
{
  double kmp,kvp,klp;
  double dav,fracg,rhoa,rhog;
  double sb,s,fracl,rhow,kk,kl,mul,le;
  double kg,mug;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);
  sb = f_sb(rhov,t);
  s = f_s(rhov,t);
  fracl = f_fracl(rhov,t);
  rhow = f_rhow(t);
  kk = f_kk(t);
  kl = f_kl(rhov,t);
  kg = f_kg(rhov,t);
  mul = f_mul(t);
  mug = f_mug(rhov,pg,t);
  le = f_le(t);
  
  switch (model){
  case 1:{
    kvp=fracg*rhov*((kk*kg/mug)-(dav/rhog/ra/t));
    klp=fracl*rhow*(kk*kl/mul);
    kmp=kvp+klp;
    break;
  }
  case 2:{
    kvp=fracg*rhov*((kk*kg/mug)-(dav/rhog/ra/t));
    klp=(1.0-(sb/s))*fracl*rhow*kk*kl/mul;
    kmp=kvp+klp;
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return(kmp);
}

double glasgowmat::kmv (double t,double pg,double rhov)
{
  double kmv,kvv,klv;
  double dav,fracg,fracl,rhoa,rhog,sb,db,s,dfracldrhov,rhow,kk,kl,mul,por,dpcdrhov;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);
  sb = f_sb(rhov,t);
  db = f_db(rhov,t);
  s = f_s(rhov,t);
  fracl = f_fracl(rhov,t);
  dfracldrhov = f_dfracldrhov(rhov,t);
  rhow = f_rhow(t);
  kk = f_kk(t);
  kl = f_kl(rhov,t);
  mul = f_mul(t);
  por = f_por(t);
  dpcdrhov = f_dpcdrhov(rhov,t);
  
  switch (model){
  case 1:{
    kvv=(fracg*dav/rhog)*(rhoa+(rv*rhov/ra));
    klv=0.0;
    kmv=kvv+klv;
    break;
  }
  case 2:{
    kvv=(fracg*dav/rhog)*(rhoa+(rv*rhov/ra));
    klv=fracl*rhow*((sb*db*dfracldrhov/s/por)-((1.0-(sb/s))*kk*kl*dpcdrhov/mul));
    kmv=kvv+klv;
    break;
  }
  default:{
      fprintf (stderr,"\n\n unknown model type in glasgowmat.cpp (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (kmv);
}

double glasgowmat::ctt (double t,double pg,double rhov)
{
  double ctt;
  double crho,le,ld,rhow,drhowdt,fracl,fracd,dfracldt,dfracddt;

  crho = f_crho(rhov,pg,t);
  rhow = f_rhow(t);
  drhowdt = f_drhowdt(t);
  fracl = f_fracl(rhov,t);
  fracd = f_fracd(t);
  dfracldt = f_dfracldt(rhov,t);
  dfracddt = f_dfracddt(t);
  ld = f_ld();
  le = f_le(t);

  ctt =  (crho+((ld+le)*(fracd*drhowdt+rhow*dfracddt))-(le*(fracl*drhowdt+rhow*dfracldt)));

  return(ctt);
}
	
double glasgowmat::ctp ()
{
  return 0.0;
}

double glasgowmat::ctv (double t,double /*pg*/,double rhov)
{
  double ctv;
  double le,rhow,dfracldrhov;

  le = f_le(t);
  rhow = f_rhow(t);
  dfracldrhov = f_dfracldrhov(rhov,t);

  ctv = (-le*rhow*dfracldrhov);

  return(ctv);
}

double glasgowmat::cat (double t,double pg,double rhov)
{
  double cat;
  double dpordt,rhoa,fracg,dfracldt;
  
  dpordt = f_dpordt(t);
  rhoa = f_rhoa(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  dfracldt = f_dfracldt(rhov,t);

  cat = ((rhoa*(dpordt-dfracldt))-(fracg*pg/ra/t/t));

  return(cat);
}

double glasgowmat::cap (double t,double /*pg*/,double rhov)
{
  double cap;
  double fracg;
  
  fracg = f_fracg(rhov,t);  

  cap = (fracg/ra/t);

  return(cap);
}

double glasgowmat::cav (double t,double pg,double rhov)
{
  double cav;
  double fracg,rhoa,dfracldrhov;

  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  dfracldrhov = f_dfracldrhov(rhov,t);

  cav = (-1.0*(fracg*rv/ra)-(rhoa*dfracldrhov));

  return(cav);
}

double glasgowmat::cmt (double t,double /*pg*/,double rhov)
{
  double cmt;
  double dpordt,rhow,drhowdt,fracl,fracd,dfracldt,dfracddt;
  
  dpordt = f_dpordt(t);
  rhow = f_rhow(t);
  drhowdt = f_drhowdt(t);
  fracl = f_fracl(rhov,t);
  fracd = f_fracd(t);
  dfracldt = f_dfracldt(rhov,t);
  dfracddt = f_dfracddt(t);
  
  cmt = ((rhov*dpordt)+((fracl-fracd)*drhowdt)+((rhow-rhov)*dfracldt)-(rhow*dfracddt));
  
  return(cmt);
}

double glasgowmat::cmp ()
{
  return 0.0;
}

double glasgowmat::cmv (double t,double /*pg*/,double rhov)
{
  double cmv;
  double fracg,rhow,dfracldrhov;
  
  fracg = f_fracg(rhov,t);
  rhow = f_rhow(t);
  dfracldrhov = f_dfracldrhov(rhov,t);

  cmv = (fracg+((rhow-rhov)*dfracldrhov));
  
  return(cmv);
}

/////////////////////////////////
//added by Tomas Krejci 28.1.2005

/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double glasgowmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,rhov,pg,t;
  new_trc = 0.0;
  
  rhov = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_mv(t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  if((ri == 0) && (ci == 2))
    new_trc = get_transmission_transcoeff_mt(t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  
  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
    new_trc = 0.0;
  if((ri == 1) && (ci == 2))
    new_trc = 0.0;
  
  if((ri == 2) && (ci == 0))
    new_trc = 0.0;
  if((ri == 2) && (ci == 1))
    new_trc = 0.0;
  if((ri == 2) && (ci == 2))
    new_trc = get_transmission_transcoeff_tt(t,pg,rhov,bc,ipp);// *scale_t;//scaling

  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double glasgowmat::transmission_nodval(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval,rhov,pg,t;
  new_nodval = 0.0;
  
  rhov = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_mv(nodval,t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 0) && (ci == 2))
    new_nodval = get_transmission_nodval_mt(nodval,t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  
  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 2))
    new_nodval = 0.0;
  
  if((ri == 2) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 2) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 2) && (ci == 2))
    new_nodval = get_transmission_nodval_tt(nodval,t,pg,rhov,bc,ipp);// *scale_t;//scaling

  return (new_nodval);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double glasgowmat::transmission_flux(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double flux,rhov,pg,t;
  flux = 0.0;
  
  rhov = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_mv(nodval,t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  if((ri == 0) && (ci == 1))
    flux = 0.0;
  if((ri == 0) && (ci == 2))
    flux = get_transmission_flux_mt(nodval,t,pg,rhov,bc,ipp);// *scale_pc;//scaling
  
  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
    flux = 0.0;
  if((ri == 1) && (ci == 2))
    flux = 0.0;
  
  if((ri == 2) && (ci == 0))
    flux = 0.0;
  if((ri == 2) && (ci == 1))
    flux = 0.0;
  if((ri == 2) && (ci == 2))
    flux = get_transmission_flux_tt(nodval,t,pg,rhov,bc,ipp);// *scale_t;//scaling

  return (flux);
}


/**
   function creates correct transfer coefficient on the boundary (transmission) for rhov

   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure
   @param t - actual temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_transcoeff_mv(double t,double pg,double /*rhov*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//mass transmission
    trc=1.0;
    break;
  }
  case 31:{//mass transmission
    double beta;

    beta = f_beta(pg,t);
    
    trc = beta;
    break;
  }
  case 32:{//mass transmission
    double beta,tinf;

    tinf = f_tinf(Tp->time); 
    beta = f_beta(pg,tinf);

    trc=beta;
    break;
  }
  case 36:{//rel.hum.
    trc=0.0;
    break;
  }    
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }

  return(trc);
}


/**
   function creates correct new nodal value on the boundary (transmission) for rhov

   @param bv - value of prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_nodval_mv(double bv,double t,double pg,double rhov,long bc,long /*ipp*/)
{
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//mass transmission
    nodval = bv;
    break;
  } 
  case 31:{//mass transmission
    nodval = bv;
    break;
  }
  case 32:{//mass transmission
    double beta,tinf;

    tinf = f_tinf(Tp->time); 
    beta = f_beta(pg,tinf);

    nodval = beta*bv;
    break;
  } 
  case 36:{//rel.hum
    double psat;

    psat = f_psat(t);
    
    nodval = bv*psat/rv/t;
    
    nodval = nodval - rhov;
    break;
  }   
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates flux on the boundary (transmission - convective mass transfer) for rhov

   @param bv - prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_flux_mv(double bv,double t,double pg,double rhov,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//mass transmission - boundary flux
    flux = (bv - rhov);
    break;
  }
  case 31:{//mass transmission - boundary flux
    double beta;

    beta = f_beta(pg,t);
    
    flux = beta*(bv - rhov);
    break;
  }
  case 32:{//mass transmission - boundary flux
    double beta,tinf;

    tinf = f_tinf(Tp->time); 
    beta = f_beta(pg,tinf);

    flux = beta*(bv - rhov);
    break;
  }    
  case 36:{//rel.hum
    double psat;

    psat = f_psat(t);
    
    flux = bv*psat/rv/t;
    
    flux = flux - rhov;
    break;
  }   
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(flux);
}

/**
   function creates correct transfer coefficient on the boundary (transmission) for t

   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure
   @param t - actual temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_transcoeff_tt(double t,double /*pg*/,double /*rhov*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    trc=1.0;
    break;
  }
  case 31:{//heat transmission
    double hqr;

    hqr = f_hqr(t,t);
        
    trc=hqr;
    break;
  }
  case 32:{//heat transmission - ISO 834 Fire Curve
    double hqr,tinf;

    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);
        
    trc=hqr;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(trc);
}



/**
   function creates correct new nodal value on the boundary (transmission) for t

   @param bv - value of prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_nodval_tt(double bv,double t,double /*pg*/,double /*rhov*/,long bc,long /*ipp*/)
{
  double nodval;
  
  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    nodval = bv;
    break;
  } 
  case 31:{//heat transmission
    nodval = bv;
    break;
  } 
  case 32:{//heat transmission - ISO 834 Fire Curve
    double hqr,tinf;

    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);

    nodval = hqr*tinf;
    break;
  } 
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(nodval);
}



/**
   function creates flux on the boundary (transmission - convective mass transfer) for t

   @param bv - prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_flux_tt(double bv,double t,double /*pg*/,double /*rhov*/,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission - boundary flux
    flux = (bv - t);
    break;
  }
  case 31:{//heat transmission - boundary flux
    double hqr;

    hqr = f_hqr(t,bv);
    
    flux = hqr*(bv - t);
    break;
  }
  case 32:{//heat transmission - boundary flux - ISO 834 Fire Curve
    double hqr,tinf;
    
    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);

    flux = hqr*(tinf - t);
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(flux);
}

/**
   coupled boundary condition - function creates correct new nodal value on the boundary (transmission) for rhov from t

   @param bv - value of prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_transcoeff_mt(double t,double pg,double rhov,long bc,long /*ipp*/)
{
  double trc;
  double kvt,keff,hqr,tinf;
  double dav,fracg,rhoa,rhog;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);

  kvt=fracg*dav*rhov*pg/rhog/ra/t/t;
  keff = f_keff (t);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission

    trc=-1.0*kvt/keff;
    break;
  }
  case 31:{//heat transmission
    hqr = f_hqr(t,t);

    trc=-1.0*kvt*hqr/keff;
    break;
  }
  case 32:{//heat transmission - ISO 834 Fire Curve  
    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);
    
    trc=-1.0*kvt*hqr/keff;
    break;
  }
  case 36:{//rel.hum
    trc=0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }

  return(trc);
}


/**
   coupled boundary condition - function creates correct new nodal value on the boundary (transmission) for rhov from t

   @param bv - value of prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_nodval_mt(double bv,double t,double pg,double rhov,long bc,long /*ipp*/)
{
  double nodval;
  double kvt,keff,hqr,tinf;
  double dav,fracg,rhoa,rhog;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);

  kvt=fracg*dav*rhov*pg/rhog/ra/t/t;
  
  keff = f_keff (t);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission

    nodval = -1.0*kvt*bv/keff;
    break;
  } 
  case 31:{//heat transmission
    nodval = -1.0*kvt*bv/keff;
    break;
  } 
  case 32:{//heat transmission - ISO 834 Fire Curve
    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);

    nodval = -1.0*kvt*hqr*tinf/keff;
    break;
  } 
  case 36:{//rel.hum
    nodval=0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(nodval);
}

/**
   coupled boundary condition - function creates correct flux on the boundary (transmission) for rhov from t

   @param bv - value of prescribed value near the boundary
   @param rhov - actual water vapour concentration
   @param pg - actual capillary gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double glasgowmat::get_transmission_flux_mt(double bv,double t,double pg,double rhov,long bc,long /*ipp*/)
{
  double flux;
  double kvt,keff,hqr,tinf;
  double dav,fracg,rhoa,rhog;

  dav = f_dav(rhov,pg,t);
  fracg = f_fracg(rhov,t);
  rhoa = f_rhoa(rhov,pg,t);
  rhog = f_rhog(rhov,pg,t);

  kvt=fracg*dav*rhov*pg/rhog/ra/t/t;
  
  keff = f_keff (t);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission - boundary flux

    flux = -1.0*kvt*(bv - t)/keff;
    break;
  }
  case 31:{//heat transmission - boundary flux
    hqr = f_hqr(t,bv);

    flux = -1.0*kvt*hqr*(bv - t)/keff;
    break;
  }
  case 32:{//heat transmission - boundary flux - ISO 834 Fire Curve
    tinf = f_tinf(Tp->time);
    hqr = f_hqr(t,tinf);

    flux = -1.0*kvt*hqr*(tinf - t)/keff;
    break;
  }
  case 36:{//rel.hum
    flux=0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param t - temperature on actual node
   @param pg - gas pressure on actual node
   @param rhov - concentration of water vapour
*/

double glasgowmat::get_othervalue(long compother,long /*ipp*/, double *r)
{
  double other;

  switch (compother){
  case 0:{//concentration of water vapour
    other = r[0];
    break;
  }
  case 1:{//gas pressure
    other = r[1];
    break;
  }
  case 2:{//temperature
    other = r[2];
    break;
  }
  case 3:{//relative humidity
    double h;
    h = f_h(r[0],r[2]);

    other = h;
    break;
  }
  case 4:{//saturation
    double s;
    s = f_s(r[0],r[2]);

    other = s;
    break;
  }
  case 5:{//vapour pressure
    double pv;
    pv = f_pv(r[0],r[2]);

    other = pv;
    break;
  }
  case 6:{//liquid water pressure
    double pl;
    pl = f_pl(r[0],r[1],r[2]);

    other = pl;
    break;
  }
  case 7:{//capillary pressure
    double pc;
    pc = f_pc(r[0],r[2]);

    other = pc;
    break;
  }    
  case 8:{//capillary pressure
    double pc;
    pc = f_pc(r[0],r[2]);

    other = pc;
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);
}

/**
     function prints names of all variables in nodes
     @param out - output file
     @param compother - number of other components
*/
void glasgowmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//concentration of water vapour
    fprintf (out,"Vapour content (kg/m3)        ");
    break;
  }
  case 1:{//gas pressure
    fprintf (out,"Gas pressure (Pa)             ");
    break;
  }
  case 2:{//temperature
    fprintf (out,"Temperature (K)               ");
    break;
  }
  case 3:{//relative humidity
    fprintf (out,"Relative humidity ()          ");
    break;
  }
  case 4:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 5:{//vapour pressure
    fprintf (out,"Pore water vapor pressure (Pa)");
    break;
  }
  case 6:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  case 7:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
    break;
  }    
  case 8:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 16/07/2014
*/
void glasgowmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_vapor_content;
  dofname[1] = trf_gas_press;
  dofname[2] = trf_temperature;
}
