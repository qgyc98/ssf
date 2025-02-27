/*
  File:             xxx.cpp
  Author:           Jiri Madera and Jaroslav kruis
  Purpose:          c
  Source:           modified Kunzel model
  Assumptions:
  1.9.2014
  
  
 */



#include "moistheat.h"
#include "globalt.h"
#include "math.h"
#include "globmatt.h"

dampermeability moistheatmat::damper;

moistheatmat::moistheatmat()
{

  
  //  density of solid material
  rho_m = 0.0;
  //  water density
  rho_w = 1000.0;
  //  water molecular weight
  M = 0.01801528;
  // universal gas constant
  R = 8.314472;
  
  // permeability influenced by damage - default option is off
  daminfl = off;
}

moistheatmat::~moistheatmat()
{
}

/**
   function reads material characteristics
   
   @param in - input file

   JM, 12.10.2011
*/
void moistheatmat::read (XFILE *in)
{
  //  density
  rho.read (in);

  //  permeability
  dpp.read (in);

  //  transport coefficient
  transpar.read (in);

  //  sorption isotherm - akumulation curfe w(pv)
  sorpiso.read (in);

  //  saturated moisture
  sm.read (in);

  //  specific heat capacity
  c.read (in);

  //  thermal conductivity
  lambda.read (in);

  // temper of transport parametr
  teptr.read(in);

  // porozimetrie of transport parametr
  poroz.read(in);

  xfscanf (in,"%m", &answertype_kwdset, &daminfl);
}


/**
   function prints material characteristics
   
   @param out - output file

   JM, 12.10.2011
*/
void moistheatmat::print (FILE *out)
{
  //  density
  rho.print (out);

  //  permeability
  dpp.print (out);

  //  transport coefficient
  transpar.print (out);

  //  sorption isotherm - akumulation curfe w(pv)
  sorpiso.print (out);

  //  saturated moisture
  sm.print (out);

  //  specific heat capacity
  c.print (out);

  //  thermal conductivity
  lambda.print (out);

  // temper of transport parametr
  teptr.print(out);

  // poroz of transport parametr
  poroz.print(out);

  fprintf(out, "%d\n", daminfl);
}

/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void moistheatmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
  }
  }
  
  if (daminfl == on){
    if (ri==0 && ci==0){
      //  conductivity matrix is modified due to damage
      damper.matcond (d,ipp);
      Tm->ip[ipp].eqother[3] = d[0][0];
    }
    if (ri==1 && ci==0){
      //  conductivity matrix is modified due to damage
      damper.matcond (d,ipp);
      Tm->ip[ipp].eqother[4] = d[0][0];
    }
    
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

    JM, 12.10.2011
*/
void moistheatmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  k = 0.0;

  if((ri == 0) && (ci == 0))
    k = kmm(ipp);
  if((ri == 0) && (ci == 1))
    k = kmt(ipp);
  if((ri == 1) && (ci == 0))
    k = khm(ipp);
  if((ri == 1) && (ci == 1))
    k = kht(ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

    JM, 12.10.2011
*/
void moistheatmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  k = 0.0;


  if((ri == 0) && (ci == 0))
    k = kmm(ipp);
  if((ri == 0) && (ci == 1))
    k = kmt(ipp);
  if((ri == 1) && (ci == 0))
    k = khm(ipp);
  if((ri == 1) && (ci == 1))
    k = kht(ipp);
  
  fillm(0.0,d);

 /* if (ipp==2)
    fprintf(Outt,"\n ri %ld ci %ld  k  %20.15le\n",ri,ci,k);
  */
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JM, 12.10.2011
*/

void moistheatmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  k = 0.0;


  if((ri == 0) && (ci == 0))
    k = kmm(ipp);
  if((ri == 0) && (ci == 1))
    k = kmt(ipp);
  if((ri == 1) && (ci == 0))
    k = khm(ipp);
  if((ri == 1) && (ci == 1))
    k = kht(ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}



/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
    JM, 12.10.2011
*/
void moistheatmat::matcap (double &c,long ri,long ci,long ipp)
{
 c=0.0;
  
  if((ri == 0) && (ci == 0))
    c =  cmm(ipp);
  if((ri == 0) && (ci == 1))
    c = cmt(ipp);
  if((ri == 1) && (ci == 0))
    c = chm(ipp);
  if((ri == 1) && (ci == 1)){
    c = cht(ipp);
  }

  //fprintf(Outt,"\n ipp %ld ri %ld ci %ld  c  %20.15le\n",ipp,ri,ci,c);
}

/**
   function checks values of relative humidity and temperature
   if the values are out of range, the limit values are assigned
   
   @param nv - array with relative humidity and temperature
   
   JK, 7.10.2011
*/
void moistheatmat::values_correction (vector &nv)
{
  double tk, pv;
  
  // partial pressure of water vapor
  pv=nv[0];
  //  temperature
  tk=nv[1];
  // tady to udelat lepe....
  if(pv >= 3.16995416e+03)
    pv = 3.1699e+03;
  if(pv <= 3.17014677e-06)
    pv = 3.18e-06;
  
  nv[0]=pv;
  nv[1]=tk;
  
}

/**
   function checks values of water vapour pressure and temperature on integration point
   if the values are out of range, the limit values are assigned
   
   TKr, 02/08/2017
*/
void moistheatmat::values_correction_ipp (long ipp)
{
  double tk, pv;
  
  // partial pressure of water vapor
  pv = Tm->ip[ipp].av[0];
  //  temperature
  tk = Tm->ip[ipp].av[1];
  // tady to udelat lepe....
  if(pv >= 3.1699e+03){
    fprintf(Outt,"Vapour pressure corrected ipp = %ld pv = %e\n",ipp,pv);
    fprintf(Outt,"Tb->lc[0].masterval = %e,Tb->lc[0].mastergrad[0] = %e\n", Tb->lc[0].masterval,Tb->lc[0].mastergrad[0]);
    fflush(Outt);
    pv = 3.1699e+03;
    abort();
  }
  if(pv <= 3.17e-06){
    fprintf(Outt,"Vapour pressure corrected ipp = %ld pv = %e\n",ipp,pv);
    fprintf(Outt,"Tb->lc[0].masterval = %e,Tb->lc[0].mastergrad[0] = %e\n", Tb->lc[0].masterval,Tb->lc[0].mastergrad[0]);
    fflush(Outt);
    pv = 3.18e-06;
  }
  Tm->ip[ipp].av[0] = pv;
}

//---------------------------------------------------------------------------
/**
   conductivity coefficient between moisture gradient and moisture flux
   
   @param ipp - id of integration point

   JM, 12.10.2011
*/
double moistheatmat::kmm (long ipp)
{
  double kk, transp, pv;//w=0.0, ttra=0.0, pv2=0.0;
  
  //  partial pressure of water vapor
  pv = Tm->ip[ipp].av[0];
  
  //  volume of unfrozen water
  //  the value is used for determination of permeability
  //w=Tm->ip[ipp].eqother[7];
  // w=0.0; 
  /*
  if (w>0.0){
    //  there is ice in pores
    //  the permeability is calculated with respect to the smallest radius of pores where no frozen water occurs
    ttra = transpar.getval (pv);
    pv2 = sorpiso.inverse_isotherm_value (w);
    if (pv<pv2)
	pv2=pv;
    transp = transpar.getval (pv2);
    //transp = 1e-14;
    
//    if (ipp < 30)
//	fprintf(Outt,"ipp = %ld, pv = %8.6e, dg = %8.6e, pv2 = %8.6e, dg2= %8.6e  \n",ipp, pv,ttra,pv2, transp);
    
  }else{
    // there is no ice
    transp = transpar.getval (pv);
  }
  */
  
  transp = transpar.getval (pv);
  kk = transp;
  
  if (ipp==0)
    fprintf (Outt,"%le  % 16.12le   % 16.12le",Tp->time,pv,kk);
  
  return (kk);
}

//---------------------------------------------------------------------------

/**
   conductivity coefficient between temperature gradient and moisture flux


   @param ipp - id of integration point

	JM, 12.10.2011
*/
double moistheatmat::kmt(long /*ipp*/)
{
  double kk;
  
     kk = 0.0;
  return (kk);
}

//---------------------------------------------------------------------------
/**
   conductivity coefficient between temperature gradient and heat flux

   @param ipp - id of integration point

    JM, 12.10.2011

 */
double moistheatmat::kht(long ipp)
{
  double moist, lamb, kk;
	// volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];
  
	lamb = lambda.getval (moist);

	//lambda+ L_v delta_p R Rov/M
	kk = lamb;

  if (ipp==0)
    fprintf (Outt,"    % 16.12le   % 16.12le\n",moist,lamb);

  return (kk);
}

//---------------------------------------------------------------------------

/**
   conductivity coefficient between moisture gradient and heat flux
   
   @param ipp - id of integration point
   
   JM, 12.10.2011
*/
double moistheatmat::khm (long ipp)
{
  double kk,tkk,lv,dp,pv;
  
  // temperature
	tkk = Tm->ip[ipp].av[1];

	// water vapor pressure
	pv = Tm->ip[ipp].av[0];

	//CorD(10,kd,moist,lambda,a2,a3);
	dp = dpp.getval (pv);

	//dp = water_vapour_permeability(tk,ipp);
	lv = latent_heat_of_evaporation_of_water (tkk);

	kk = lv*dp;

	if (ipp==0)
	  fprintf (Outt,"    % 16.12le",kk);
  	
  return (kk);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between moisture gradient and moisture flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double moistheatmat::cmm(long ipp)
{
  double cc,tk;
  double pv, smc, w;
	// partial pressure of water vapor
	 pv = Tm->ip[ipp].av[0];
	// temperature
	tk = Tm->ip[ipp].av[1];

	// akumulacni parametr - derivative
	//akumpar = akump.getval (moist);
	cc = Tm->ip[ipp].eqother[1]*rho_w;
	w = sorpiso.isotherm_value (pv);

	smc = sm.getval(pv);
// 	fprintf(Outt,"cc je %3.2e pro pv= %lf, w = %lf, smc = %lf\n",cc,pv, w, smc);

	cc = cc + (smc-w)*M/R/tk;
//	fprintf(Outt,"cc++ je %3.2e pro pv= %lf, w = %lf, smc = %lf\n",cc,pv, w, smc);

	return (cc);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between temperature gradient and moisture flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double moistheatmat::cmt(long /*ipp*/)
{
  double cc;

  cc = 0.0;
	return (cc);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between moisture gradient and heat flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double moistheatmat::chm(long /*ipp*/)
{
  double cc;
  
      cc = 0.0;
  
  return (cc);
}

 //---------------------------------------------------------------------------
/**
   capacity coefficient between temperature gradient and heat flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double moistheatmat::cht(long ipp)
{
  double  cc1;

  //cc1 = derivative_of_the_enthalpy_density(ipp);
  cc1 = derivative_of_the_enthalpy_density_h(ipp);

  return (cc1);
}

//---------------------------------------------------------------------------
/*
   function computes the latent heat of evaporation
   
   @param tk - temperature
*/
double moistheatmat::latent_heat_of_evaporation_of_water(double tk)
{
  return (2.5008e6)*pow((273.15/tk),(0.167+tk*3.67e-4));
}


//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/**
   Function calculates saturation water vapor pressure
   saturation water vapor pressure by xxxx
   tk ... temperature [K]

   
   @param tk - temperature [K]
*/
double  moistheatmat::saturated_water_vapor_pressure(double tk)
{
     return (exp(23.5771 - 4042.9/(tk - 37.58)));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/**
   Function calculates derivative of the enthalpy density
   
   @param tk - temperature [K]
   @param ipp - id of integration point
   
   JM, 12.10.2011
*/
double moistheatmat::derivative_of_the_enthalpy_density (long ipp)
{
  double cap, moist,rho_m;

  // volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];
  
  // cap - specific heat capacity
  cap = c.getval (moist);
  //  density
  rho_m = rho.getval (0.0);
  
  return (rho_m * cap);
}


double moistheatmat::derivative_of_the_enthalpy_density_h (long ipp)
{
  double c_eff, cap, moist,rho_m, l_ice, c_w, tk, c_effmax, f;
  int t_kod,e;

  l_ice = 334000.0; ///[J/kg]

  c_w = 4180.0; //[J/kg/K]

  ro_ice = 934.0; // [kg/m3]

  c_ice = 2110.0; // [J/kg/K]

  beta = 1.11; // objemova zmena
  
  f = 0.375;
  e = 2; 
  // temperature
  tk = Tm->ip[ipp].av[1];
  
// volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];
  
  // cap - specific heat capacity
  cap = c.getval (moist);

    //  density
  rho_m = rho.getval (0.0);
  
  del_l = l_ice*moist*rho_w/rho_m; 
  
  cl = (rho_m*cap+rho_w*c_w*moist)/(rho_m+rho_w*moist);

  cs = (rho_m*cap+ro_ice*c_ice*moist*beta)/(rho_m+ro_ice*moist*beta);

 gibbs_thomson(ipp);

// t_1 = 0.0; // doplnit
// t_2 = 10.0; // doplnit
 
 if (tk<=t_1){
    t_kod = 0;
 }
 else{
     if (tk<=((t_1+t_2)/2)){
	t_kod = 1;
     }
     else{
	if (tk<t_2){
	    t_kod = 2;
	}
	else{
	    t_kod = 3;
	}
     }
} 
 // fprintf(Outt,"del_l = %lf, t_kod =  %ld,t= %lf  t_1= %lf,  t_2= %lf  \n",del_l,t_kod,tk, t_1, t_2);

switch (t_kod){
    case 0:{ c_eff = cs;
    break;
    }
    case 1:{ 
    c_effmax = del_l/(f*(t_2-t_1))+(cs+cl)/2; 
    c_eff =pow((1-cos(2*M_PI*(tk-t_1)/(t_2-t_1))),e)*((c_effmax-cs)/(pow(2,e))) +cs;
/*    if (ipp == 800){
    fprintf(Outt,"del_l = %lf, c_effmax= %3.2e, c_eff  = %3.2e,  t_1= %lf,  t_2= %lf  \n",del_l,c_effmax,c_eff, t_1, t_2);
    }*/
    break;
    } 
    case 2:{ 
    c_effmax = del_l/(f*(t_2-t_1))+(cs+cl)/2; 
    c_eff =pow((1-cos(2*M_PI*(t_2-tk)/(t_2-t_1))),e)*((c_effmax-cl)/(pow(2,e))) +cl;
/*    if (ipp == 800){
    fprintf(Outt,"del_l = %lf, c_effmax= %3.2e, c_eff  = %3.2e,  t_1= %lf,  t_2= %lf  \n",del_l,c_effmax,c_eff, t_1, t_2);
    }*/
    break;
    }
    case 3:{ c_eff = cl;
    break;
    } 
}

    return (rho_m * c_eff);
    
    //fprintf(Outt,"c_effmax je %3.2e , c_eff= %3.2e \n\n\n",c_effmax,c_eff);
}

/**
   Gibbs-Thomson equation for evaluation of freezing point depression
   
   @param ipp - integration point id
*/
void moistheatmat::gibbs_thomson (long ipp)
{
  double por, moist, t_cr;

  // volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];
 
  // porimezimetricka krivka
  //  for a given volumetric moisture content, the function returns the pore radius
  por = poroz.getval (moist);

  //  freezing point depression
  t_cr = (-2*273.15*0.0317*0.000018)/(6010*por/2);
  t_1 = t_cr-1.5+273.15;
  t_2 = t_cr+1.5+273.15;
  //fprintf(Outt,"t_cr je %3.2e , t_1= %3.2e t_2 = %lf  \n",t_cr,t_1,t_2);
}

/** 
    function computes volume change caused by water freezing
    
    @param w - volumetric moisture content
    @param t - temperature
    @param pt - temperature from the previous time step
    @param rl - lower radius of ice filling a pore
    @param ru - upper radius of ice filling a pore
    @param icevol - volume of unfrozen water which is used for determination of actual permeability
    
    14. 3. 2019
*/
double moistheatmat::freezing_volume_increment (double w, double t, double pt, double &rl, double &ru, double &icevol)
{
  double r,d,t_cr,wb,wt,dv;
  double threshold=1.0e-20;
  
  if (rl<threshold && ru<threshold){
    //  there is no ice from the previous steps
    
    //  for a given volumetric moisture content, the function returns the pore diameter
    d = poroz.getval (w);
    
    //  freezing point depression (Gibbs-Thomson equation)
    //  T_0 = 273.15
    //  gamma_sl = 0.0317
    //  nu = 18e-6
    //  delta H = 6010
    t_cr = 2.0*273.15*0.0317*0.000018/(6010.0*d/2.0);
    
    if (t<273.15-t_cr){
      //  the actual temperature is less than freezing point
      
      //  radius of pores with unfrozen water (Gibbs-Thomson equation)
      r = 2.0*273.15*0.0317*0.000018/(6010.0*(273.15-t));
      //  volumetric moisture equal to pore space with unfrozen water
      wb = poroz.getinvval (2.0*r);
      
      //  volume increment caused by frozen water
      //  volume increment is positive
      dv = (w-wb)*0.0905;

      rl=r;
      ru=d/2.0;
      icevol=wb;
    }else{
      //  the actual temperature is greater than freezing point
      icevol=0.0;
      dv = 0.0;
    }
  }else{
    //  there is ice from previous steps
    
    if (t<pt){
      //  the actual temperature is less than the previous one

      //  radius of pores with unfrozen water (Gibbs-Thomson equation)
      r = 2.0*273.15*0.0317*0.000018/(6010.0*(273.15-t));
      //  volumetric moisture equal to pore space with unfrozen water
      wb = poroz.getinvval (2.0*r);
      //  volumetric moisture equal to pore space with unfrozen water from previous step
      wt = poroz.getinvval (2.0*rl);
      //  volume increment caused by frozen water
      //  the volume increment is positive
      dv = (wt-wb)*0.0905;
      
      rl=r;
      icevol=wb;
    }else{
      //  the actual temperature is greater than the previous one
      
      if (t<273.15){
	
	//  radius of pores with unfrozen water (Gibbs-Thomson equation)
	r = 2.0*273.15*0.0317*0.000018/(6010.0*(273.15-t));
	
	if (r<ru){
	  //  the radius for the actual temperature is less than the upper radius attained
	  //  the temperature does not cause total thawing
	  
	  //  volumetric moisture equal to pore space with unfrozen water from previous step
	  wb = poroz.getinvval (2.0*rl);
	  //  volumetric moisture equal to pore space for the actual radius
	  wt = poroz.getinvval (2.0*r);
	  
	  //  volume increment caused by frozen water
	  //  the volume increment is negative
	  dv = (wb-wt)*0.0905;
	  
	  rl=r;
	  icevol=wt;
	}else{
	  //  the radius for the actual temperature is greater than the upper radius attained
	  //  the temperature causes total thawing
	  
	  //  volumetric moisture equal to pore space with unfrozen water from previous step
	  wb = poroz.getinvval (2.0*rl);
	  //  volumetric moisture equal to pore space under upper radius
	  wt = poroz.getinvval (2.0*ru);
	  
	  //  volume increment caused by frozen water
	  //  the volume increment is negative
	  dv = (wb-wt)*0.0905;
	  
	  rl=0.0;
	  ru=0.0;
	  icevol=w;
	}
      }
      else{
	//  t is greater than 273.15 and therefore, there is no ice without any condition
	
	//  volumetric moisture equal to pore space with unfrozen water from previous step
	wb = poroz.getinvval (2.0*rl);
	//  volumetric moisture equal to pore space under upper radius
	wt = poroz.getinvval (2.0*ru);
	
	//  volume increment caused by frozen water
	//  the volume increment is negative
	dv = (wb-wt)*0.0905;
	
	rl=0.0;
	ru=0.0;
	icevol=w;
      }
    }
  }
  
  return dv;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



/**
   function computes new nodal value (for transmission_vector)
   for boundary condition (third kind of boundary condition)
   
   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of the boundary condition
   
*/
double moistheatmat::transmission_nodval (double nodval,long ri,long ci,long nid,long bc)
{

  double pv,new_nodval,t;
  new_nodval = 0.0;
  
  //partial pressure of water vapor
  pv = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);

  if((ri == 0) && (ci == 0))
	new_nodval = get_transmission_nodval_hh (nodval,pv,t,bc);
  if((ri == 0) && (ci == 1))
	new_nodval = 0.0;

  if((ri == 1) && (ci == 0))
	new_nodval = get_transmission_nodval_th (nodval,bc);
  if((ri == 1) && (ci == 1))
	new_nodval = get_transmission_nodval_tt (nodval,bc);

  return (new_nodval);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (partial pressure of water vapor)

   @param bv - prescribed value near the boundary
   @param rh - relative humidity in domain solved
   @param t - temperature in domain solved
   @param bc - type of boundary condition
   
*/
double moistheatmat::get_transmission_nodval_hh (double bv,double rh,double t,long bc)
{
  double nodval,pgws,moist;

  switch (bc){
  case 5:{
    // water vapour pressure transmission from climatic conditions type 2
    nodval = bv;
    break;
  }
  case 30:{
	// relative humidity -> pgw

   //	pgws = saturated_water_vapor_pressure(t);
   //	nodval = bv*pgws;
	nodval = bv;
	break;
  }
  case 31:{
	//pgw
	nodval = bv;
	break;
  }
  case 32:{
	//relative humidity -> pgw

	pgws = saturated_water_vapor_pressure(t);
	nodval = pgws*rh;
	bv = pgws*bv;
	//inverse eq.
	nodval = bv - nodval;

	break;
  }
  case 33:{
	//h -measurements

	//  density
	rho_m = rho.getval (0.0);
	// volumetric moisture content
	moist = bv/100.0*rho_m/rho_w;
	//relative humidity from inverse sorption isotherm;
	//inverse_sorption_izoterms_giva_data(0,nodval,t,moist,a3,ipp);
	nodval = sorpiso.inverse_isotherm_value (moist);
	break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return nodval;
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (relative humidity)

   @param bv - prescribed value near the boundary
   @param rh - relative humidity in domain solved
   @param t - temperature in domain solved
   @param bc - type of boundary condition
   
*/
double moistheatmat::get_transmission_nodval_th (double bv,long bc)
{
  double nodval;

  switch (bc){
  case 5:{
    nodval = bv;
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return nodval;
}

/**
   function returns prescribed temperature of external environment for the transmisson boundary condition

   @param bv - prescribed external temperature
   @param bc - type of boundary condition
   
   TKr, modified by JK
*/
double moistheatmat::get_transmission_nodval_tt (double bv,long bc)
{
  double nodval;
  
  switch (bc){
  case 5:{
    nodval = bv;
    break;
  }
  case 30:{
    nodval = bv;
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return (nodval);
}



/**
   function computes new transmission coefficient
   
   @param trc - prescribed transmission coefrhcient on the boundary
   @param ri  - row index
   @param ci  - column index
   @param nid - number of node
   @param bc  - type of boundary condition
   
*/
double moistheatmat::transmission_transcoeff (double trc,long ri,long ci,long /*nid*/,long bc)
{
  //long k;
  double new_trc;
  //double pv,t;
  new_trc = 0.0;
  /*
  // partial pressure of water vapor
  pv = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);
  */
  if((ri == 0) && (ci == 0))
	new_trc = get_transmission_transcoeff_hh (trc,bc);
  if((ri == 0) && (ci == 1))
	new_trc = 0.0;

  if((ri == 1) && (ci == 0))
	new_trc = get_transmission_transcoeff_th (trc,bc);
  if((ri == 1) && (ci == 1))
	new_trc = get_transmission_transcoeff_tt (trc,bc);

  return (new_trc);
}

/**
   function creates transfer coefrhcient on the boundary for prescribed condition (relative humidity)

   @param t - temperature
   @param bc - type of boundary condition

*/
double moistheatmat::get_transmission_transcoeff_hh (double t,long bc)
{
  double trc,pgws;

  switch (bc){
  case 5:{
    trc = t;
    break;
  }
  case 30:{
	// relative humidity
	//pgws = saturated_water_vapor_pressure (t);
	//trc = pgws;
	trc = t;
	break;
  }
  case 31:{
	//water vapour pressure pgw
	pgws = saturated_water_vapor_pressure (t);
	trc = pgws;
	trc = 1.0;
	break;
  }
  case 32:{
	//water vapour pressure pgw
	trc = 0.0;
	break;
  }
  case 33:{
	//h - measurements
	trc = 1.0;
	break;
  }
  default:{
	print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
	exit(0);
  }
  }

  return trc;
}

/**
   function creates transfer coefficient on the boundary for prescribed condition (relative humidity)

   @param t - temperature
   @param bc - type of boundary condition

*/
double moistheatmat::get_transmission_transcoeff_th (double t,long bc)
{
  double trc;

  switch (bc){
  case 5:{
    trc = t;
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return trc;
}

/**
   function creates transfer coefrhcient on the boundary for prescribed condition (temperature)

   @param trc - prescribed transmission coefficent
   @param bc  - type of boundary condition

*/
double moistheatmat::get_transmission_transcoeff_tt (double trcp,long bc)
{
  double trc;

  switch (bc){
  case 5:{
    trc = trcp;
    break;
  }
  case 30:{
	//heat transmission
	trc = trcp;
	break;
  }
  default:{
	print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
	exit(0);
  }
  }

  return trc;
}



/**
   function computes flux (for transmission_vector)
   for boundary condition (third kind of boundary condition)

   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition

*/
double moistheatmat::transmission_flux (double nodval,long ri,long ci,long nid,long bc)
{
  
  double flux,pv,t;

  flux = 0.0;

  // partial pressure of water vapor
  pv = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);

  if((ri == 0) && (ci == 0))
	flux = get_transmission_flux_hh (nodval,pv,t,bc);
  if((ri == 0) && (ci == 1))
	flux = 0.0;

  if((ri == 1) && (ci == 0))
	flux = 0.0;
  if((ri == 1) && (ci == 1))
	flux = get_transmission_flux_tt (nodval,t,bc);

  return flux;
}

/**
   function creates flux on the boundary (convective mass transfer)

   @param bv - prescribed value near the boundary
   @param pv - relative humidity in domain solved
   @param t - temperature in domain solved
   @param bc - type of boundary condition

*/
double moistheatmat::get_transmission_flux_hh (double /*bv*/,double /*pv*/,double /*t*/,long /*bc*/)
{
  double flux;
  //double flux,pgws,moist,nodval;

  /* tady to udelat znova
  double flux,pgws,moist,nodval;

  switch (bc){//type of prescribed variable
  case 30:{
	//relative humidity -> pgw
	pgws = saturated_water_vapor_pressure (t);

	flux = pgws*rh;
	bv = pgws*bv;
	flux = bv - flux;//inverse eq.

	break;
  }
  case 31:{
	//water vapour pressure pgw
	pgws = saturated_water_vapor_pressure (t);
	flux = pgws*rh;
	flux = bv - flux;
	break;
  }
  case 32:{
	//relative humidity -> pgw
	pgws = saturated_water_vapor_pressure(t);
	flux = pgws*rh;
	bv = pgws*bv;
	flux = bv - flux;//inverse eq.

	break;
  }
  case 33:{
	//h -measurements

	//  density
	rho_m = rho.getval (0.0);
	// volumetric moisture content
	moist = bv/100*rho_m/rho_w;

	//relative humidity from inverse sorption isotherm;
	//inverse_sorption_izoterms_giva_data(0,nodval,t,moist,a3,ipp);
	nodval = sorpiso.inverse_isotherm_value (moist);

	flux = nodval - rh;
	break;
  }
  
  default:{
	print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
	exit(0);
  }
  }
  */
    flux = 0.0;

  return flux;
}




/**
   function creates heat flux on the boundary

   @param bv - prescribed value near the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

*/
double moistheatmat::get_transmission_flux_tt (double bv,double t,long bc)
{
  double flux_tt;

  switch (bc){//type of prescribed variable
  case 30:{
	//heat transmission
	flux_tt = (bv - t);
	break;
  }
  default:{
	print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
	exit(0);
  }
  }
  return(flux_tt);
}


/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - rhrst integration point on element
   @param pc - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double moistheatmat::get_othervalue(long compother,double rh,double t, long ipp)
{
  double other,w;

  switch (compother){
  case 0:{//relative humidity
    other = rh/saturated_water_vapor_pressure (teptr.getval(0));
	break;
  }
  case 1:{//temperature
	other = Tm->ip[ipp].av[1]-273.15;
	  break;
  }
  case 2:{//moisture content w = m3/m3
	//sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
	w = sorpiso.isotherm_value (rh);
	other = w;
	break;
  }
  case 3:{//  volume increment of ice
    other = Tm->ip[ipp].eqother[4];
    break;
  }
  case 4:{//moisture content u = kg/kg (u.rho_m/rho_w = w)
	//  sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
	rho_m = rho.getval (0.0);
	//sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
	w = sorpiso.isotherm_value (rh);
		other = w*rho_w/rho_m;
	break;
  }
  case 5:{
  	// t = 298.15;
	double lamb;
	 w = sorpiso.isotherm_value (rh);
	lamb = lambda.getval (w);
	 other = lamb;
	  break;
  }
  case 6:{//Moisture concentration
	  other = rh*saturated_water_vapor_pressure (t)*M/(R*t);
	  break;
  }
  case 7:{//Pore water vapor pressure
	  other = R*t/(saturated_water_vapor_pressure (t)*M)/rh;
	  break;
  }
  case 8:{//temperature
	  other = 0.0;
	  break;
  }
  case 9:{//temperature
	  other = 0.0;
	  break;
  }
  case 10:{//temperature
	  other = 0.0;
	  break;
  }
  default:{
	print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
  return (other);

}

/**
     function prints names of all variables in nodes
     @param out - output rhle
     @param compother - number of other components
*/
void moistheatmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//Relative humidity
    fprintf (out,"Relative humidity ()       ");
    break;
  }
  case 1:{//Temperature
    fprintf (out,"Temperature (K)             ");
	break;
  }
  case 2:{//Moisture content w
   fprintf (out,"Moisture content w (m3/m3)   ");
    break;
  }
  case 3:{//water vapour pressure pgw
    fprintf (out,"Water vapour pressure (Pa)  ");
    break;
  }
  case 4:{//Moisture content u 
   fprintf (out,"Moisture content u (kg/kg)   ");
    break;
  }
  case 5:{// Thermal conductivity L
   fprintf (out,"Thermal conductivity L (W/m/K)   ");
    break;
  }
  case 6:{// Moisture concentration
   fprintf (out,"Moisture concentration rho_gw (kg/m3)   ");
    break;
  }
  case 7:{// Pore water vapor pressure
   fprintf (out,"Pore water vapor pressure p^gw (Pa)   ");
    break;
  }
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
}


//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

/**
   function computes auxilary values
   all input values have to be imported as parameters of the function
   because it is used in connection with integration points and nodes
   therefore, the input values cannot be read in this function
   
   @param ipp - integration point id
   @param inv - actual values
   @param inp - values from the previous time step
   @param ine - values of the array other
   
   JM, 12.10.2011
*/
void moistheatmat::aux_values (long /*ipp*/,double *inv,double *inp,double *ine,double *out)
{
  double tk,tp,rl,ru,w,dwdf, smc, tep;
  double pv;
  
  // partial pressure of water vapour
  pv = inv[0];
  // temperature 
  tk = inv[1];
  
  //  temperature at previous time step
  tp = inp[1];
  
  //  lower radius of ice filling a pore
  rl = ine[5];
  //  upper radius of ice filling a pore
  ru = ine[6];
  
  //  saturated volumetric moisture content
  smc = sm.getval (0.0);
  //  w - volumetric moisture content
  w = sorpiso.isotherm_value (pv);
  //  derivative of sorption isotherm with respect to relative humidity
  dwdf = sorpiso.derivative_isotherm_value (pv);
  // temperature of transport parameter
  tep = teptr.getval(0.0);
  //  volume increment caused by frozen water
  //dv = freezing_volume_increment (w,tk,tp,rl,ru);

  //if (ipp==1525)
  //fprintf (Outt,"%15.12le %15.12le  %15.12le   %15.12le\n",Tp->time,tk,tp,dv);
  
  if (w>smc)
    w = smc;
  if (w<0.0)
    w = 0.0;
  
  // volumetric moisture content
  out[0] = w;
  // derivative of sorption isoterm
  out[1] = dwdf;
  // saturated volumetric moisture content
  out[2] = smc;
  // temperature of transport parameter
  out[3] = tep;
  //  volume increment caused by frozen water
  //out[4] = dv;
  //  lower radius of ice filling a pore
  //out[5] = rl;
  //  upper radius of ice filling a pore
  //out[6] = ru;
}

/**
   function stores values to eqother arrays in integration points

   @param ipp - integration point id
   @param out - array with values

   7.10.2011
*/
void moistheatmat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  //Tm->ip[ipp].eqother[4]=out[4];
  //Tm->ip[ipp].eqother[5]=out[5];
  //Tm->ip[ipp].eqother[6]=out[6];
  
}

/**
   function reads values from arrays in integration points

   @param ipp - integration point id
   @param av - actual values
   @param pv - previous values
   @param eq - array with values

   7.10.2011
*/
void moistheatmat::give_values (long ipp,double *av,double *pv,double *eq)
{
  av[0] = Tm->ip[ipp].av[0];
  av[1] = Tm->ip[ipp].av[1];

  pv[0] = Tm->ip[ipp].pv[0];
  pv[1] = Tm->ip[ipp].pv[1];

  eq[0] = Tm->ip[ipp].eqother[0];
  eq[1] = Tm->ip[ipp].eqother[1];
  eq[2] = Tm->ip[ipp].eqother[2];
  eq[3] = Tm->ip[ipp].eqother[3];
  eq[4] = Tm->ip[ipp].eqother[4];
  eq[5] = Tm->ip[ipp].eqother[5];
  eq[6] = Tm->ip[ipp].eqother[6];

}

/*
  function defines initial values
  
  @param ipp - integration point id
  
  11.7.2012
*/
void moistheatmat::initvalues (long /*ipp*/,long /*ido*/)
{
}
/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  //pv - partial pressure of water vapor
  Created by Tomas Koudelka, 4.2.2014
*/
void moistheatmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_press_water_vapor;
  dofname[1] = trf_temperature;
}



/**
   Function returns temperature in the given integration point.
   
   @param ipp - integration point id

   @return Funtion returns value of temperature stored in the integartion point.

   24. 10. 2013, JK
*/
double moistheatmat::give_temperature (long ipp)
{
  return Tm->ip[ipp].av[1];
}



/**  
  Function returns initial temperature in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of initial temperature stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double moistheatmat::give_inittemperature (long /*ipp*/)
{
  print_err("Not yet implemented", __FILE__, __LINE__, __func__);
  abort();
}



/**
  Function returns relative humidity in the given integration point.

  @param ipp - integration point id

  @return Funtion returns value of relative humidity stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double moistheatmat::give_rel_hum (long ipp)
{
  double pv, pvs, relhum, t;
  pv = Tm->ip[ipp].av[0];
  t =  teptr.getval(0.0);
  pvs = saturated_water_vapor_pressure (t);
  relhum = pv/pvs;

  return relhum;               // opravit
}

/**
  Function returns partial pressure of water vapor in the given integration point.

  @param ipp - integration point id

  @return Funtion returns value of partial pressure of water vapor stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double moistheatmat::give_press_water_vapor (long ipp)
{
  return Tm->ip[ipp].av[0];
}



/**  
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integration point.

  Created by TKo, 25.11.2013
*/
double moistheatmat::give_vol_moist (long ipp)
{
  return Tm->ip[ipp].eqother[0];
}


/**  
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integartion point.

  Created by TKo, 12.2.2019
*/
double moistheatmat::give_volume_change (long ipp)
{
  return Tm->ip[ipp].eqother[4];
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   23. 4. 2014, by TKo
*/
void moistheatmat::give_reqntq(long *antq)
{
  if (daminfl == on){
	//  damage parameter
    antq[scal_iso_damage-1] = 1;
    //  process zone length
    antq[proc_zone_length-1] = 1;
    //  crack width
    antq[crack_width-1] = 1;
  }
}
