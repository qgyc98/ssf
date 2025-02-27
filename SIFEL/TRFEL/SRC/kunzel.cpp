/*
  File:             kunzel.cpp
  Author:           Jiri Madera
  Purpose:          c
  Source:           Kunzel, 19xx;
  Assumptions:
  27.9.2003
  
  
  position CORD:  2 - density
		  3 - porosity
		  4 - water vapour diffusion resistance factor
		  5 - moisture diffusivity
		  6 - sorption isoterm
		  7 - saturated moisture
		  8 - none
		  9 - specific heat capacity
		 10 - thermal conductivity
		 11 - 13 - none
		 14 Dcoef
		 15 - binding isotherm
		 16 - cfmax
		 17 ws
		 18 - kmm type
		 19 - kunzeltype
*/



#include "kunzel.h"
#include "globalt.h"
#include "math.h"
#include "globmatt.h"

dampermeability kunmat::damper;

kunmat::kunmat()
{
  long i;
  
  //  type of the model is not defiend
  kunzeltype = 0;
  
  //  type of the moisture function is not defiend
  kmmtype = 0;
  
  //  density of solid material
  rho_m = 0.0;
  //  water density
  rho_w = 1000.0;
  //  water molecular weight
  M = 0.01801528;
  // universal gas constant
  R = 8.314472;
  
  //  correction coefficient for the sorption isotherm
  asi = 0.0;
  //  correction coefficient for the desorption isotherm
  adsi = 0.0;
  //  correction coefficient for the moisture diffusivity
  akappa = 0.0;
  //  correction coefficient for the moisture diffusivity of drying
  adkappa = 0.0;
  
  //  minimum time before hysteresis
  time_hyst =0.0;
  
  //  influence of damage on permeability
  daminfl=off;

  for (i=0;i<20;i++){
    MatChar[i]=0;
    data[i] = NULL;
  }
}

kunmat::~kunmat()
{
  for (int i=0; i<20; i++)
    if (data[i])
      delete data[i];
}

/**
   function reads material characteristics
   
   @param in - input file

   JM, 12.10.2011
*/
void kunmat::read (XFILE *in)
{
  
  xfscanf (in,"%ld %ld",&kunzeltype,&kmmtype);
  
  if (kunzeltype==8) {
    //  density
    rho.read (in);
    //  permeabilita...
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
  }
  else{
    //  density
    rho.read (in);
    
    //  porosity
    por.read (in);
    
    //  water vapour diffusion resistance factor
    mu.read (in);

    //  moisture diffusivity
    kappa.read (in);

    if (kunzeltype == 11){
      //  minimum time before hysteresis
      xfscanf (in,"%lf",&time_hyst);
    
      //  correction coefficient for the moisture diffusivity
      xfscanf (in,"%lf",&akappa);
    
      //  moisture diffusivity for drying
      kappadry.read (in);
    
      //  correction coefficient for the moisture diffusivity of drying
      xfscanf (in,"%lf",&adkappa);
    }
  
    //  sorption isotherm
    sorpiso.read (in);

    if (kunzeltype == 11){
      //  correction coefficient for the sorption isotherm
      xfscanf (in,"%lf",&asi);
    
      //  desorption isotherm
      desorpiso.read (in);
    
      //  correction coefficient for the desorption isotherm
      xfscanf (in,"%lf",&adsi);
    }
  
    //  saturated moisture
    sm.read (in);
    
    //  specific heat capacity
    c.read (in);

    //  thermal conductivity
    lambda.read (in);
  }
  //  damage influence
  xfscanf (in,"%m",&flagsw_kwdset,&daminfl);
}


/**
   function prints material characteristics
   
   @param out - output file

   JM, 12.10.2011
*/
void kunmat::print (FILE *out)
{
  
  fprintf (out,"\n %ld %ld\n",kunzeltype,kmmtype);
  
  if (kunzeltype==8) 
  {
    //  density
    rho.print (out);
    //  permeabilita...
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
  }
  else
  {
    //  density
    rho.print (out);
    //  porosity
    por.print (out);
    //  water vapour diffusion resistance factor
    mu.print (out);
    //  moisture diffusivity
    kappa.print (out);

    if (kunzeltype == 11){
      //  minimum time before hysteresis
      fprintf (out,"\n %lf",time_hyst);
    
      //  correction coefficient for the moisture diffusivity
      fprintf (out,"\n %lf",akappa);
    
      //  moisture diffusivity for drying
      kappadry.print (out);

      fprintf (out,"\n %lf",adkappa);
    }
  
  
    //  sorption isotherm
    sorpiso.print (out);

    if (kunzeltype == 11){
      fprintf (out,"\n %lf",asi);
    
      //  desorption isotherm
      desorpiso.print (out);
    
      fprintf (out,"\n %lf",adsi);
    }
  
    //  saturated moisture
    sm.print (out);
    //  specific heat capacity
    c.print (out);
    //  thermal conductivity
    lambda.print (out);
  }
  
  //  damage influence
  fprintf (out," %d\n",daminfl);
}









/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void kunmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    //  conductivity matrix is modified due to damage
    if((ri == 0) && (ci == 0))
      damper.matcond (d,ipp);
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
void kunmat::matcond1d (matrix &d,long ri,long ci,long ipp)
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
void kunmat::matcond2d (matrix &d,long ri,long ci,long ipp)
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

  //if (ipp==0)
  //fprintf(Outt,"\n ri %ld ci %ld  k  %20.15le\n",ri,ci,k);
  
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

void kunmat::matcond3d (matrix &d,long ri,long ci,long ipp)
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
void kunmat::matcap (double &c,long ri,long ci,long ipp)
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

  //if (ipp==0)
  //fprintf(Outt,"\n ri %ld ci %ld  c  %20.15le\n",ri,ci,c);
}

/**
   function checks values of relative humidity and temperature
   if the values are out of range, the limit values are assigned
   
   @param nv - array with relative humidity and temperature
   
   JK, 7.10.2011
*/
void kunmat::values_correction (vector &nv)
{
  double rh, tk, pv, pvs;
  
  switch (kunzeltype){
  case 1:{
    //  relative humidity
    rh=nv[0];
    //  temperature
    tk=nv[1];
    
    if (rh >= 1.0)
      rh = 1.0;
    if (rh <= 0.0)
      rh = 0.0;
    
    if (tk >= 350.0)
      tk = 350.0;
    if (tk <= 240.0)
      tk = 240.0;
    
    nv[0]=rh;
	nv[1]=tk;
    break;
  }
  case 2:
  case 3:
  case 4:{
    // partial pressure of water vapor
    pv=nv[0];
    //  temperature
    tk=nv[1];
    
    pvs = saturated_water_vapor_pressure(tk);
    // relative humidity
    rh = pv/pvs;
    
    if(rh >= 1.0)
      pv = pvs;
    if(rh <= 0.0)
      //	pv = pvs*0.01;
      
      if (tk >= 350.0)
	tk = 350.0;
    if (tk <= 240.0)
      tk = 240.0;
    
    nv[0]=pv;
    nv[1]=tk;
    break;
  }
  case 5:
  case 6:
  case 7:{
    // partial pressure of water vapor
    pv=nv[0];
    //  temperature
    tk=nv[1];
    
    
    nv[0]=pv;
    nv[1]=tk;
    break;
  }
  case 8:{
	// partial pressure of water vapor
	pv=nv[0];
	//  temperature
	tk=nv[1];

	if(pv >= 3.16995416e+03)

	  pv = 3.1699e+03;
	if(pv <= 3.17014677e-06)
	  pv = 3.18e-06;

	nv[0]=pv;
	nv[1]=tk;
	break;
  }
  default:{
    print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
  }
  }

}

//---------------------------------------------------------------------------
/**
   conductivity coefficient between moisture gradient and moisture flux
   
   @param ipp - id of integration point

   JM, 12.10.2011
*/
double kunmat::kmm (long ipp)
{
  double rh, dp, dmretc, ps, dphi, moist, kk, tk, kapa, pv;
  double A,B,r1,r2;

  switch (kunzeltype) {
  case 1:
  case 11:{

	// relative humidity
	rh = Tm->ip[ipp].av[0];

	// temperature
	tk = Tm->ip[ipp].av[1];

	if (rh>1.0)
	  rh = 1.0;
	if (rh<0.0)
	  rh = 0.0;

	// tady je zmena v transportu plyne a kapalne faze...
	// od 0 - cca90% jen plynna a od 90 - 97.6 castecne obe a od 907.6 - 1.0 jen kapalna

	r1 = 0.90;
	r2 = 0.976;

	// volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];
	A = 0.0;
	B = 0.0;

	if (rh>=0.0){
	  if (rh<=r1){
	switch (kmmtype){
	case 0:{
	  A=1.0;
	  B=1.0;
	  break;
	}
	case 1:
	case 2:
	case 3:{
	  A = 1.0;
	  B = 0.0;
	  break;
	}
	default:{
	  print_err("unknown type of Kunzel model is required",__FILE__,__LINE__,__func__);
	}
	}//  end of switch (kmmtype){
	  }//  end of if (rh<=r1){
	  else{
	if (rh > r2){
	  switch (kmmtype){
	  case 0:{
		A=1.0;
		B=1.0;
		break;
	  }
	  case 1:
	  case 2:
	  case 3:{
		A = 0.0;
		B = 1.0;
		break;
	  }
	  default:{
		print_err("unknown type of Kunzel model is required",__FILE__,__LINE__,__func__);
	  }
	  }//  end of switch (kmmtype){
	}//  end of if (rh > r2){
	else{
	  switch (kmmtype){
	  case 0:{
		A=1.0;
		B=1.0;
		break;
	  }
	  case 1:{
		B = (1/(r2-r1))*(rh-r1);
		A = 1-B;
		break;
	  }
	  case 2:{
		B= 0.5*sin(M_PI/(r2-r1)*(rh-r1)-M_PI/2)+0.5;
		A=1-B;
		break;
	  }
	  case 3:{
		if (rh<0.938){
		  B=32*pow((1/(r2-r1))*(rh-r1),6);
		  A=1-B;
		}
		else{
		  B=1-32*pow((1/(r2-r1))*(r2-rh),6);
		  A=1-B;
		}
		break;
	  }
	  default:{
	    print_err("unknown definition of Material KUNZEL is required",__FILE__,__LINE__,__func__);
	  }
	  }//  end of switch (kmmtype){
	}//  end of else
      }
	}
    else{
      print_err("unknown error",__FILE__,__LINE__,__func__);
      printf("rh = %lf\n\n",rh);
      abort();
    }
    
    
    // volumetric moisture content
   // moist =  Tm->ip[ipp].eqother[0];
    // moisture diffusivity
    kapa =  Tm->ip[ipp].eqother[3];
    
    dp = water_vapour_permeability(tk,ipp);
    ps = saturated_water_vapor_pressure( tk);

    //derivative of sorption isotherms
	dmretc =  Tm->ip[ipp].eqother[1];
    //    dmretc = 455.0;
	dphi = kapa * dmretc*rho_w;
    kk = dphi*B + A*dp*ps;
    
    break;
  }
    /*
  case 2:{
    
    // partial pressure of water vapor
   // pv = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    // moisture diffusivity
    dw =  Tm->ip[ipp].eqother[3];

    dp = water_vapour_permeability(tk,ipp);

    kk = dp+dw*M/(R*tk);
	break;
 }
 case 3:{

    // partial pressure of water vapor
    //pv = Tm->ip[ipp].av[0];

    // temperature
    //tk = Tm->ip[ipp].av[1];

    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];

    // mi - water vapour diffusion resistance factor
    //CorD(18,kd,moist,D,a2,a3);
    //kk = D;
    kk = mu.getval (moist);
    break;
 }
 case 4:{

    // partial pressure of water vapor
    //pv = Tm->ip[ipp].av[0];

	// temperature
	tk = Tm->ip[ipp].av[1];

    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];

    // mi - water vapour diffusion resistance factor
    CorD(18,kd,moist,D,a2,a3);
    kk = D*M/(R*tk);;
    break;
 }
 case 5:
 case 6:{

	// partial pressure of water vapor
	//pv = Tm->ip[ipp].av[0];

	// temperature
	tk = Tm->ip[ipp].av[1];

	dp = water_vapour_permeability(tk,ipp);

	 kk = dp;
	break;
 }
 case 7:{


	// partial pressure of water vapor
	//pv = Tm->ip[ipp].av[0];
	moist= 0.05;
	CorD(18,kd,moist,D,a2,a3);

	// testovani
	dp = 2.217e-11;
	D = dp*R*293.15/M;

	kk = D;
	break;
 }
*/
 case 8:{

	// partial pressure of water vapor
	pv = Tm->ip[ipp].av[0];

	// temperature
	//tk = Tm->ip[ipp].av[1];
	 double transp;
	 transp = transpar.getval (pv);

	 //dp = water_vapour_permeability(tk,ipp);

	 kk = transp;
	break;
 }
 default:{
	  print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }
 }

return (kk);
}

//---------------------------------------------------------------------------

/**
   conductivity coefficient between temperature gradient and moisture flux


   @param ipp - id of integration point

	JM, 12.10.2011
*/
double kunmat::kmt(long ipp)
{
  double kk, tk, rh;
  double dp, dpvs;
  
  switch (kunzeltype) {
  case 1:
  case 11:{
    
    // relative humidity
    rh = Tm->ip[ipp].av[0];
    
	// temperature
    tk = Tm->ip[ipp].av[1];
    
    if (rh>1.0)
      rh = 1.0;
    if (rh<0.0)
      rh = 0.0;
    
	dpvs = derivative_of_saturated_water_vapor_pressure(tk);
    dp = water_vapour_permeability(tk,ipp);
    
    // delta_p phi dp_vs/dT
    kk = dp * rh * dpvs;
    break;
  }
  case 2:{
      double pv,tk, dw;
      // partial pressure of water vapor
      pv = Tm->ip[ipp].av[0];
      // temperature
      tk = Tm->ip[ipp].av[1];

      // moisture diffusivity
      dw =  Tm->ip[ipp].eqother[3];

      kk = -dw*M*pv/(R*tk*tk);

      break;
  }

  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:{
	   kk = 0.0;

      break;
  }
  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }

  }

  //	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);
// fprintf(Outt,"Cas je %d sekund \n",Tp->time);
//fprintf(Outt,"kmt je %3.2e pro  tk= %lf \n",kk,tk);
  return (kk);
}

//---------------------------------------------------------------------------
/**
   conductivity coefficient between temperature gradient and heat flux

   @param ipp - id of integration point

    JM, 12.10.2011

 */
double kunmat::kht(long ipp)
{
  double kk,tk,rh;
  double lv1, lamb, dpvs, dp, moist, rov;
  
  switch (kunzeltype){
  case 1:
  case 11:{
    
    // relative humidity
    rh = Tm->ip[ipp].av[0];
    
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    if (rh>1.0)
      rh = 1.0;
    if (rh<0.0)
      rh = 0.0;

	// volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];
	// lambda - thermal conductivity
	lamb = lambda.getval (moist);

	dpvs =  derivative_of_saturated_water_vapor_pressure(tk);
	dp =  water_vapour_permeability(tk,ipp);

	lv1 = latent_heat_of_evaporation_of_water (tk);

	//lambda+ L_v delta_p phi dp_vs/dT
	kk = lamb + lv1 * dp * rh * dpvs;
	break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:{
	double  moist;
	// volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];

	// lambda - thermal conductivity
	//CorD(10,kd,moist,lambda,a2,a3);
	kk  = lambda.getval (moist);
	break;
  }
  case 7:{

	// relative humidity
	rov = Tm->ip[ipp].av[0];

	// temperature
	tk = Tm->ip[ipp].av[1];

	// volumetric moisture content
	moist =  0.01;
	// lambda - thermal conductivity
	//CorD(10,kd,moist,lambda,a2,a3);
	lamb = lambda.getval (moist);
	dp =  water_vapour_permeability(tk,ipp);

	lv1 = latent_heat_of_evaporation_of_water (tk);

	//lambda+ L_v delta_p R Rov/M
	kk = lamb + lv1 * dp * R * rov/M;

	break;
 }
  case 8:{

	// relative humidity
	//pv = Tm->ip[ipp].av[0];

	// temperature
	//tk = Tm->ip[ipp].av[1];

	// lambda - thermal conductivity
	//CorD(10,kd,moist,lambda,a2,a3);
	// volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];

	lamb = lambda.getval (moist);

	//lambda+ L_v delta_p R Rov/M
	kk = lamb;

	break;
 }

 default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }
 }

  return (kk);
}

//---------------------------------------------------------------------------

/**
   conductivity coefficient between moisture gradient and heat flux
   
   @param ipp - id of integration point
   
   JM, 12.10.2011
*/
double kunmat::khm (long ipp)
{
  double kk,tk;
  double lv,dp,ps,moist,D, pv;
  
  switch (kunzeltype){
  case 1:
  case 11:{
    
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    dp = water_vapour_permeability(tk,ipp);
    ps = saturated_water_vapor_pressure (tk);
    lv = latent_heat_of_evaporation_of_water (tk);
    
    //L_v delta_p p_vs
    kk = lv * dp * ps;
    
    break;
  }
  case 2:
  case 3:
  case 4:
  case 5:{
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    dp = water_vapour_permeability(tk,ipp);
    lv = latent_heat_of_evaporation_of_water (tk);
    
    kk = lv*dp;
    
    break;
  }
  case 6:{
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    dp = water_vapour_permeability(tk,ipp);
    lv = latent_heat_of_evaporation_of_water (tk);
    
    kk = lv*dp;
    
    break;
  }
  case 7:{

    // relative humidity
    // rh = Tm->ip[ipp].av[0];
    
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    // if (rh>1.0) rh = 1.0;
    // if (rh<0.0) rh = 0.0;
    
    moist= 0.05;
    CorD(18,kd,moist,D,a2,a3);
    lv = latent_heat_of_evaporation_of_water (tk);
    dp = water_vapour_permeability(tk,ipp);
    //L_v delta_p p_vs
    kk = lv * dp*R*tk/M;
    
    break;
  }
  case 8:{
  //double dpp;
	// temperature
	tk = Tm->ip[ipp].av[1];

	// water vapor pressure
	pv = Tm->ip[ipp].av[0];

	//CorD(10,kd,moist,lambda,a2,a3);
	dp = dpp.getval (pv);

	//dp = water_vapour_permeability(tk,ipp);
	lv = latent_heat_of_evaporation_of_water (tk);

	kk = lv*dp;

	break;
  }
  default:{
    print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
    
  }
  }
  
  return (kk);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between moisture gradient and moisture flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::cmm(long ipp)
{
  double cc,tk;
  double moist=0.0, pi;

  switch (kunzeltype){
  case 1:
  case 11:{
    // derivative sorption isotherms * density of water
	cc = Tm->ip[ipp].eqother[1]*rho_w;
    break;
  }
  case 2:{
    double tk;
    // partial pressure of water vapor
    // pv = Tm->ip[ipp].av[0];
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    cc = M/(R*tk);
    
    break;
  }
  case 3:{
    
    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];
    
    // pi - porosity
    pi = por.getval (moist);
    
    cc = (pi-moist);
    
    break;
  }
  case 4:{
    // partial pressure of water vapor
    // pv = Tm->ip[ipp].av[0];
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];
    
    // pi - porosity
    pi = por.getval (moist);
    
    cc = (pi-moist)*M/(R*tk);;
    
    break;
  }
  case 5:{
	// partial pressure of water vapor
    // pv = Tm->ip[ipp].av[0];
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];
    
    // pi - porosity
    pi = por.getval (moist);
    
    cc = (pi)*M/(R*tk);;
    
    break;
  }
  case 6:{
    // partial pressure of water vapor
    // pv = Tm->ip[ipp].av[0];
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];
    
    // pi - porosity
    pi = por.getval (moist);
    
    cc = (pi)*M/(R*tk);
    
    break;
  }
  case 7:{
	// partial pressure of water vapor
	// pv = Tm->ip[ipp].av[0];
	// temperature
	// tk = Tm->ip[ipp].av[1];
    
        // volumetric moisture content
	moist =  Tm->ip[ipp].eqother[0];
    

	// pi - porosity
	pi = por.getval (moist);

	cc = pi;
	// testovani
	  cc = 1.0;

	  break;
  }
  case 8:{
	double pv;
	// partial pressure of water vapor
	 pv = Tm->ip[ipp].av[0];
	// temperature
	// tk = Tm->ip[ipp].av[1];

	// akumulacni parametr - derivative
	//akumpar = akump.getval (moist);
	cc = Tm->ip[ipp].eqother[1]*rho_w;
 /*	t = Tp->time;
	if (t>1.1e5) {
	   fprintf(Outt,"cc je %3.2e pro pv= %lf\n",cc,pv);
	}
   */
	//cc = akumpar;
	// testovani
	//  cc = 1.0;

	  break;
  }
  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }


 }

  return (cc);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between temperature gradient and moisture flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::cmt(long ipp)
{
  double cc,tk,pv;

  switch (kunzeltype){
  case 1:
  case 11:{
	cc = 0.0;

	break;
  }
  case 2:{
	// partial pressure of water vapor
	pv = Tm->ip[ipp].av[0];
	// temperature
	tk = Tm->ip[ipp].av[1];

	cc = -M*pv/(R*tk*tk);
	break;
  }
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:{
	cc = 0.0;
	break;
  }

  default:{
	print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
  }

  }

  return (cc);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between moisture gradient and heat flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::chm(long /*ipp*/)
{
  double cc;

  switch (kunzeltype){
  case 1:
  case 11:{
	cc = 0.0;
	break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:{
    cc = 0.0;
    
    break;
  }
    
  default:{
    print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
    
  }
    
  }
  
  return (cc);
}

 //---------------------------------------------------------------------------
/**
   capacity coefficient between temperature gradient and heat flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::cht(long ipp)
{
  double cc;

  switch (kunzeltype) {
  case 1:
  case 11:{
	cc = derivative_of_the_enthalpy_density(ipp);
    
    break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:{
    cc = derivative_of_the_enthalpy_density(ipp);
    
    break;
  }
    
  default:{
    print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
  }
    
  }
  
  return (cc);
}

//---------------------------------------------------------------------------
 /**
   function computes derivative of the saturated water vapor pressure

      @param tk - temperature
 */

double kunmat::derivative_of_saturated_water_vapor_pressure(double tk)
{

   return ((4042.9 * exp(23.5771-4042.9/(tk - 37.58)))/((tk-37.58)*(tk-37.58)));
}

//---------------------------------------------------------------------------

/**
   function computes water vapour permeability delta_p = delta_a/mu
   delta_a - vapour permeability of air
   mu - water vapour diffusion resistance factor

   @param rh - relative humidity
   @param tk - temperature
   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::water_vapour_permeability (double tk,long ipp)
{
  double P, Rv, Pa;
  double da, dp, mi;

    // volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];

  //mi -  water vapour diffusion resistance factor
  //CorD(4,kd,moist,mi,a2,a3);
  mi = mu.getval (moist);

  // Rv = R/M , where R - universal gas constant, M - water molecular weight
  Rv = R/M;

  // air pressure
  Pa = 101325.0;
  P = 101325.0;

  da = (2.306e-5 * Pa)/(Rv * tk * P)*pow((tk/273.15),1.81);

  if (mi==0.0)
    dp=0.0;
  else
    dp = da/mi;

  return (dp);
}
//---------------------------------------------------------------------------
/*
   function computes the latent heat of evaporation
   
   @param tk - temperature
*/
double kunmat::latent_heat_of_evaporation_of_water(double tk)
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
double  kunmat::saturated_water_vapor_pressure(double tk)
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
double kunmat::derivative_of_the_enthalpy_density (long ipp)
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
double kunmat::transmission_nodval (double nodval,long ri,long ci,long nid,long bc)
{
  double new_nodval,h,t;
  new_nodval = 0.0;
  
  //  humidity
  h = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);
  
  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_hh (nodval,h,t,bc);
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  
  if((ri == 1) && (ci == 0))
    //new_nodval = get_transmission_nodval_th (nodval,bc);
    new_nodval = 0.0;//must be corrected??!!
  if((ri == 1) && (ci == 1))
    new_nodval = get_transmission_nodval_tt (nodval,bc);
  
  return (new_nodval);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (relative humidity)
   
   @param bv - prescribed value near the boundary
   @param rh - relative humidity in domain solved
   @param t - temperature in domain solved
   @param bc - type of boundary condition
   
*/
double kunmat::get_transmission_nodval_hh (double bv,double rh,double t,long bc)
{
  double nodval,pgws,moist;
  
  switch (bc){
  case 5:{
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
double kunmat::get_transmission_nodval_th (double bv,long bc)
{
  double nodval;
  
  switch (bc){
  case 4:{
    nodval = bv;
    break;
  }
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
  return nodval;
}

/**
   function returns prescribed temperature of external environment for the transmisson boundary condition
   
   @param bv - prescribed external temperature
   @param bc - type of boundary condition
   
   TKr, modified by JK
*/
double kunmat::get_transmission_nodval_tt (double bv,long bc)
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
double kunmat::transmission_transcoeff (double trc,long ri,long ci,long nid,long bc)
{
  double new_trc,h,t;
  new_trc = 0.0;
  
  //  humidity
  h = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);
  
  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_hh (trc,bc);
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  
  if((ri == 1) && (ci == 0))
    //new_trc = get_transmission_transcoeff_th (trc,bc);
    new_trc = 0.0;//must be corrected??!!
  if((ri == 1) && (ci == 1))
    new_trc = get_transmission_transcoeff_tt (trc,bc);
  
  return (new_trc);
}

/**
   function creates transfer coefrhcient on the boundary for prescribed condition (relative humidity)
   
   @param t - temperature
   @param bc - type of boundary condition
   
*/
double kunmat::get_transmission_transcoeff_hh (double t,long bc)
{
  double trc;
  
  switch (bc){
  case 5:{
    trc = t;
    break;
  }
  case 30:{
    // relative humidity
    // double pgws;
    //pgws = saturated_water_vapor_pressure (t);
    //trc = pgws;
    trc = t;
    trc = 1.0;//debug??!!

    break;
  }
  case 31:{
    //water vapour pressure pgw
    //pgws = saturated_water_vapor_pressure (t);
    //trc = pgws;
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
double kunmat::get_transmission_transcoeff_th (double t,long bc)
{
  double trc;
  
  switch (bc){
  case 4:{
    trc = t;
    break;
  }
  case 5:{
    trc = t;
    break;
  }
  case 30:{
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
double kunmat::get_transmission_transcoeff_tt (double trcp,long bc)
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
double kunmat::transmission_flux (double nodval,long ri,long ci,long nid,long bc)
{
  double flux,h,t;
  
  flux = 0.0;
  
  //  humidity
  h = nodalval (nid,0);
  //  temperature
  t = nodalval (nid,1);
  
  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_hh (nodval,h,t,bc);
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
   @param rh - relative humidity in domain solved
   @param t - temperature in domain solved
   @param bc - type of boundary condition
   
*/
double kunmat::get_transmission_flux_hh (double bv,double rh,double t,long bc)
{
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
  
  return flux;
}




/**
   function creates heat flux on the boundary

   @param bv - prescribed value near the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

*/
double kunmat::get_transmission_flux_tt (double bv,double t,long bc)
{
  double flux;
  
  switch (bc){//type of prescribed variable
  case 30:{
    //heat transmission
    flux = (bv - t);
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - rhrst integration point on element
   @param pc - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double kunmat::get_othervalue(long compother,double rh,double t, long /*ipp*/)
{
  double other,w=0.0;
  
  switch (compother){
  case 0:{//relative humidity
    other = rh*100;
    break;
  }
  case 1:{//temperature
    other = t-273.15;
    break;
  }
  case 2:{//moisture content w = m3/m3
    //sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    w = sorpiso.isotherm_value (rh);
    other = w;
    break;
  }
  case 3:{//water vapour pressure pgw
    other = rh * saturated_water_vapor_pressure (t);
    break;
  }
  case 4:{//moisture content u = kg/kg (u.rho_m/rho_w = w)
    //  sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    rho_m = rho.getval (0.0);
    w = sorpiso.isotherm_value (rh);
    other = w*rho_w/rho_m;
    break;
  }
  case 5:{
    t = 298.15;
    other = rh/saturated_water_vapor_pressure (t);
    break;
  }
  case 6:{//temperature
    other = rh*saturated_water_vapor_pressure (t)*M/(R*t);
    break;
  }
  case 7:{//temperature
    other = rh*R*t/(saturated_water_vapor_pressure (t)*M);
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
void kunmat::print_othervalue_name(FILE *out,long compother)
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
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
}


//---------------------------------------------------------------------------

/**
   function computes exponential coefficient of moisture diffzusivity
   
   @param a,b - function coefficients
   @param rh - relative humidity
   
   7.10.2011
*/
double kunmat::kapa_exp(double a, double b, double rh)
{
  return (a * exp(b * rh));
}

//---------------------------------------------------------------------------

/**
   function determines required material parameter
   
    position CORD:  2 - density
		    3 - porosity
		    4 - water vapour diffusion resistance factor
		    5 - moisture diffusivity
		    6 - sorption isoterm
		    7 - saturated moisture
		    8 - none
		    9 - specific heat capacity
		    10 - thermal conductivity
		    11 - 13 - none
		    14 Dcoef
		    15 - binding isotherm
		    16 - cfmax
		    17 ws
		    18 - none
		    19 - kunzeltype
		     
   @param charid - 
   @param kvyhl -
   @param in - 
   @param x - independent variable
   @param y - dependent variable
   @param
   @param

   JM
*/
void kunmat::CorD(long charid,long &kvyhl, double x, double &y, double &z, double &z2)
{

  switch (MatChar[charid]){
  case 0:{
    kvyhl = 0;
    break;
  }
  case 1:{
    //  material parameter is given by a constant
    //  the function getval requires an argument, in this case, it is omitted
    y = data[charid]->getval(0.0);
    kvyhl = 1;
    break;
  }
  case 2:{
    //  material parameter is given by a table
    
    y = data[charid]->getval(x);
    kvyhl = 2;
    
    //  number of rows in the table
    long nr=data[charid]->tabf[0].asize;
    //  the last value in the list
    //  it contains hygroscopic moisture
    z2 = data[charid]->tabf[0].y[nr-1];

    break;
  }
  case 30:{
    y = MatFunce[charid][0];
    z = MatFunce[charid][1];
    kvyhl = 30;
    break;
  }
  case 31:{
    y = MatFunce[charid][0];
    z = MatFunce[charid][1];
    if (charid == 6)
      {
	z2 = MatFunce[charid][2];
      }
    kvyhl = 31;
    break;
  }
  case 32:{
	y = MatFunce[charid][0];
    z = MatFunce[charid][1];
    kvyhl = 32;
    break;
  }
  case 33:{
    kvyhl = 33;
    //zatim nic
    break;
  }
  default:{
    print_err("unknown definition of material parameter is required",__FILE__,__LINE__,__func__);
  }
  }
}

/*
void kunmat::sorption_izoterms_giva_data(long kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp)
{

//kod ..... 0- jen vlhkost, 1 - jen derivaci...

  long s;
  double hmrh, mhmc, smc, u, a, n;
  
  if (rh < 0) rh = 0;
  
  CorD(6,s,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN
  
  switch (s){
  case 2:{
    if(rh>0.976)
      {
	
	CorD(7,kd,rh,smc,a2,a3);
	CorD(6,kd,0.976,mhmc,a2,mhmc);
	give_data(rh,mhmc, smc, 0.976,moist);
	dmoistdrh  =1000 * (smc - mhmc)/(1 - 0.976);
	// toto pak smazat,,,,,,
	//dmoistdrh = 455.0;
      }
    else
      {
	CorD(6,s,rh,moist,dmoistdrh,a3);
	//dmoistdrh =1000* derivative_sorption_izoterm_data(rh);
	dmoistdrh =1000* sorpiso. derivative_isotherm_value (rh);
      }
    break;
  }
  case 30:{
    CorD(6,s,rh,a1,hmrh,a3);
    CorD(7,s,rh,smc,a2,a3);
    CorD(6,s,rh,mhmc,hmrh,a3);
    give_data(rh,mhmc, smc, hmrh,moist);
    if (rh < hmrh)
      {
	dmoistdrh = sorpiso.derivative_isotherm_value (double in);
derivative_of_the_sorption_isotherm_root (rh, tk);
      }
    else {
      dmoistdrh = derivative_of_retention_curve_root(rh, tk, ipp);
    }
    
    break;
  }
  case 31:{
    switch (kod){
    case 0:{
      CorD(6,kod,0,u,a,n);
      //moist = si_kk_hansen(rh, 0,u,a,n);
      moist = sorpiso.sorption_isotherm_value (rh);
      //moist = sihans.hansen_sorption_isotherm (rh);
      //dmoistdrh = sihans.derivative_relhum(rh);
      break;
    }
    case 1:{
      //dmoistdrh = sorption_izotherm_derivative_kk_hansen(rh, tk);
      dmoistdrh = sorpiso.derivative_sorption_isotherm_value (rh);
      break;
    }
    default:{
      print_err("unknown definition of Sorption isotherm HANSEN is required",__FILE__,__LINE__,__func__);
    }
    }
    break ;
  }
  default:{
    print_err("unknown definition of Sorption isotherm  is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
}
*/

/**
   Function calculates relative humidity from water content (inverse relation form sorption isotherm).

   @param kod       - code (0 or 1)
   @param rh        - water content
   @param tk        - temperature
   @param moist     - water content
   @param dmoistdrh - dearivative of water content with respect to relative humidity
   @param ipp       - intpoint number
   
   01/08/2013,TKr
*/
/*
void kunmat::inverse_sorption_izoterms_giva_data(long kod,double & rh, double tk, double moist, double & dmoistdrh, long ipp)
{
  long s;
  
  if (rh < 0)
    rh = 0;
  
  // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN
  CorD(6,s,0,a1,a2,a3);
  
  switch (s){
  case 2:{
    
    //  material parameter is given by a table, it returns inverse value
    rh = data[6]->tabf[0].getval3(moist,dmoistdrh);
    
    break;
  }
  default:{
    print_err("unknown definition of Sorption isotherm  is required",__FILE__,__LINE__,__func__);
  }
  }
}
*/

/*
double kunmat::si_kk_hansen(double rh, double tk,double u, double a, double n)
{
    return (u*pow((1-(log(rh))/a),(-1/n)));
}
*/


/**
   function computes derivative of the sorption isotherm

   @param x1 - moisture content

   7.10.2011
*/
/*
double kunmat::derivative_sorption_izoterm_data(double x1)
{
  long charid;
  double derbi;
  
  //  sorption isotherm
  charid=6;
  derbi = data[charid]->getderiv (x1);

  return derbi;
}
*/

/**
   function computes auxilary values
   all input values have to be imported as parameters of the function
   because it is used in connection with integration points and nodes
   therefore, the input values cannot be read in this function
   
   @param ipp - integration point id

   JM, 12.10.2011
*/
void kunmat::aux_values (long ipp,double *inv,double *inp,double *ine,double *out)
{
  double rh,tk,w,dwdf, kapaA, smc;
  
  switch (kunzeltype){
  case 1:{
	//  Kunzel model is formulated in relative humidity and temperature


	// relative humidity
	rh = inv[0];
	// temperature
	tk = inv[1];

	if (rh>1.0)
	  rh = 1.0;
	if (rh<0.0)
	  rh = 0.0;

	//  saturated volumetric moisture content
	smc = sm.getval (0.0);
	//  w - volumetric moisture content
	w = sorpiso.isotherm_value (rh);
	//  derivative of sorption isotherm with respect to relative humidity
	dwdf = sorpiso.derivative_isotherm_value (rh);

	if (w>smc)
	  w = smc;
	if (w<0.0)
	  w = 0.0;

	//  moisture diffusivity
	kapaA = kappa.getval (w);

	// volumetric moisture content
	out[0] = w;
	// derivative of sorption isoterm
	out[1] = dwdf;
    // saturated volumetric moisture content
	out[2] = smc;
    // moisture diffusivity
	out[3] = kapaA;
    //  type of Kunzel model
	out[4] = kunzeltype;

	break;
  }
	/*
  case 2:
  case 3:
  case 4:{

	// partial pressure of water vapor
	pv = Tm->ip[ipp].av[0];
	// temperature
	tk = Tm->ip[ipp].av[1];

	pvs = saturated_water_vapor_pressure(tk);
	// relative humidity
	rh = pv/pvs;
	if (rh>1) {
      rh = 1;
    }
    if (rh<0) {
      rh = 0.0;
    }
    // saturated volumetric moisture content
    smc = data[7]->getval(0.0);
    // w - volumetric moisture content
    sorption_izoterms_giva_data(0,rh,tk,w,dwdf, ipp);
    if (w>smc) w = smc;
    if(w<0) w = 0.0;
    
    kapa_values(MatChar[5],ipp,w, kapaA);
    
    // volumetric moisture content
    out[0] = w;
    // derivative of sorption isoterms
    out[1] = rh;
    // saturated volumetric moisture content
    out[2] = smc;
    // moisture diffusivity
    out[3] = kapaA;
    
    break;
  }
  case 5:
  case 6:
  case 7:{
    double pv,tk, pvs;
    double w,dwdf, kapaA, smc, rh;
    
    // partial pressure of water vapor
    pv = Tm->ip[ipp].av[0];
    // temperature
    tk = Tm->ip[ipp].av[1];
    
    // saturated volumetric moisture content
    smc = data[7]->getval(0.0);
    // w - volumetric moisture content
    
    // volumetric moisture content
    out[0] =smc;
    // derivative of sorption isoterms
    
    out[1] = 0.1;
    // saturated volumetric moisture content
    out[2] = smc;
    // moisture diffusivity
    out[3] = -10000.0;
    
    break;
  }
*/
  case 8:{
	//  Kunzel model is formulated in relative humidity and temperature

	double pv;
	// bla bla bla
	pv = inv[0];
	// temperature
	tk = inv[1];

	//  saturated volumetric moisture content
	smc = sm.getval (0.0);
	//  w - volumetric moisture content
	w = sorpiso.isotherm_value (pv);
	//  derivative of sorption isotherm with respect to relative humidity
	dwdf = sorpiso.derivative_isotherm_value (pv);

	if (w>smc)
	  w = smc;
	if (w<0.0)
	  w = 0.0;

	//  moisture diffusivity
   //	kapaA = kappa.getval (w);

	// volumetric moisture content
	out[0] = w;
	// derivative of sorption isoterm
	out[1] = dwdf;
	// saturated volumetric moisture content
	out[2] = smc;
	// moisture diffusivity
	out[3] = 10000.0;
    //  type of Kunzel model
	out[4] = kunzeltype;

	break;
  }
  case 11:{
	//  Kunzel model is formulated in relative humidity and temperature with hysteresis

	//  saturated volumetric moisture content
	smc = sm.getval (0.0);

	hystereze (ipp,inv,inp,ine,out);

	out[2] = smc;
	out[4] = kunzeltype;

	break;
  }
  default:
    print_err("unknown type of Kunzel model is required",__FILE__,__LINE__,__func__);
  }
}

/**
   function stores values to eqother arrays in integration points

   @param ipp - integration point id
   @param out - array with values

   7.10.2011
*/
void kunmat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
}

/**
   function reads values from arrays in integration points

   @param ipp - integration point id
   @param av - actual values
   @param eq - array with values

   7.10.2011
*/
void kunmat::give_values (long ipp,double *av,double *pv,double *eq)
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
}

/*
  function defines initial values
  
  @param ipp - integration point id
  
  11.7.2012
*/
void kunmat::initvalues (long ipp,long /*ido*/)
{
  Tm->ip[ipp].eqother[4] = kunzeltype;
}


/**
   function evaluates the sorption isotherm
   
   @param in - input value
   
   10. 10. 2013
*/
/*
double kunmat::sorption_isotherm_value (double in)
{
  double si;
  
  si = isoth.isotherm_value (in);
  
  return si;
}
*/


/**
   function solves the hysteresis
   
   @param ipp - integration point id
   
   17. 10. 2013
*/
void kunmat::hystereze (long /*ipp*/,double *inv,double *inp,double *ine,double *out)
{
  double rh,rhp,w,wa,wp,wd,dsi,ddsi,xi,xikappa,kapa,dwdf,kapaa,kapap,kapad,dkdwa,dkdwd;
  double logkapap,logkapaa,logkapad,logdkdwa,logdkdwd;
  double time,act_time,start_time;
    xi = 0.0;  
  //  actual time
  act_time = Tp->time;
  //  start time (of the analysis
  start_time = Tp->timecont.start_time;
  //  time from the start
  time = act_time - start_time;
  
  //  actual relative humidity
  rh = inv[0];
  
  if (time < time_hyst){
    //  hysteresis is not taken into account in this time period
    
    //  w - volumetric moisture content
    w = sorpiso.isotherm_value (rh);
    //  moisture diffusivity
    kapa = kappa.getval (w);
    //  derivative of sorption isotherm with respect to relative humidity
    dwdf = sorpiso.derivative_isotherm_value (rh);
    
  }
  else{
    //  hysteresis is taken into account in this time period
    
    //  volumetric moisture content from the previous time
    wp = ine[0];
    //  relative humidity from the previous time
    rhp = inp[0];
    
    //  wa - volumetric moisture content for adsorption cycles
    wa = sorpiso.isotherm_value (rh);
    //  wd - volumetric moisture content for desorption cycles
    wd = desorpiso.isotherm_value (rh);
    //  derivative of the sorption isotherm with respect to the relative humidity
    dsi = sorpiso.derivative_isotherm_value (rh);
    //  derivative of the desorption isotherm with respect to the relative humidity
    ddsi = desorpiso.derivative_isotherm_value (rh);
    
    //  slope of the hysteretic parameter
    dwdf = ((wp-wa)*(wp-wa)*adsi*ddsi + (wp-wd)*(wp-wd)*asi*dsi)/((wd-wa)*(wd-wa));
    
    //  volumetric moisture content
    w = wp + xi*(rh-rhp);
    
    
    //  moisture diffusivity from the previous time
    kapap = ine[3];
    //  moisture diffusivity for adsorption cycles at the actual time
    kapaa = kappa.getval (w);
    //  moisture diffusivity for desorption cycles at the actual time
    kapad = kappadry.getval (w);
    //  derivative of the moisture diffusivity for adsorption cycles with respect to the moisture content
    dkdwa = kappa.getderiv (w);
    //  derivative of the moisture diffusivity for desorption cycles with respect to the moisture content
    dkdwd = kappadry.getderiv (w);
    
    logkapap = log (kapap);
    logkapaa = log (kapaa);
    logkapad = log (kapad);
    logdkdwa = log (dkdwa);
    logdkdwd = log (dkdwd);
    
    //  slope of the hysteretic parameter
    xikappa = ((logkapap-logkapaa)*(logkapap-logkapaa)*adkappa*logdkdwd + (logkapap-logkapad)*(logkapap-logkapad)*akappa*logdkdwa)/((logkapad-logkapaa)*(logkapad-logkapaa));
    
    //kapa = kapap + xikappa
    

  }

  //  w - volumetric moisture content
  out[0] = w;
  //  derivative of sorption isotherm with respect to relative humidity
  out[1] = dwdf;
  //  moisture diffusivity
  out[3] = kapa;

  
/*
    if (w > wp){
      //  navlhani
      hyst=0;
    }else{
      //  vysychani
      hyst=1;
    }
    
    if (hyst==0){
      //  navlhani
      
    case 0:{			// zvysovani promene....
	if(logkapap  == ypvANorm)
	  {
	    yavNorm = yavANorm;	
	  }
	else {
	  newder = (k1[matchar]*derpvDNorm*(ypvNorm-ypvANorm)*(ypvNorm-ypvANorm)+k2[matchar]*derpvANorm*(ypvNorm -ypvDNorm)*(ypvNorm - ypvDNorm))/((ypvDNorm-ypvANorm)*(ypvDNorm-ypvANorm));
	  //yavNorm = (xav-xpv)*newder*k3[matchar]  + ypvNorm;
	  p = -1*newder; 
	  logP = log(p);
	  yavNorm =ypvNorm + logP*(xav-xpv)*k3[matchar];
	  
	  
	  if(yavNorm > yavANorm) yavNorm = yavANorm;
	  if(yavNorm < yavDNorm) yavNorm = yavDNorm;
	}
	break;
      }
      case 1:{			// klesani zakl promene
	if(ypvNorm == ypvDNorm)
	  {
	    yavNorm = yavDNorm;
	  }
	else {
	  newder  = (k4[matchar]*derpvDNorm*(ypvNorm-ypvANorm)*(ypvNorm-ypvANorm)+k5[matchar]*derpvANorm*(ypvNorm -ypvDNorm)*(ypvNorm - ypvDNorm))/((ypvDNorm-ypvANorm)*(ypvDNorm-ypvANorm));
	  
	  
	  p = -1*newder; 
	  logP = log(p);
	  yavNorm =ypvNorm + logP*(xav-xpv)*k6[matchar];
	  //yavNorm = (xav-xpv)*newder*k6[matchar]  + ypvNorm;
	  
	  if(yavNorm < yavDNorm) yavNorm = yavDNorm;
	  if(yavNorm > yavANorm) yavNorm = yavANorm;
	}
	break;
      }
      }

    yav = exp(yavNorm);
    outvalue = yav;
    outvalue2 = -10000.0;
    if (t>15559200)
      {
	if(t< 15724800)
	  {
	    if(ipp == 80)
	      fprintf(Outt2,"%lf %e  %e  %e  %e\n",t, xav,yavA,yav, yavD);
	  }
	else{
	  if(t> 31536000 && t< 31622400)
	    {
	      if(ipp == 80)
		fprintf(Outt2,"%lf %e  %e  %e  %e\n",t, xav,yavA,yav, yavD);
	    }
	}
	
      }
    
    //if(ipp == 0)
    //	  fprintf(Outt,"cas t %lf ipp je %ld   xav %e  a ypv %e  A KAPA JE %e \n",t, ipp,xav,ypv, outvalue);
    
    
    if(xav > xpv)            
      {
	hyst = 0; // navlhani
      }
    else	{
      hyst = 1; // vysouseni
    }
    switch(hyst)
      {
      case 0:{			// zvysovani promene....
	if(ypv == ypvA)                    
	  {
	    yav = yavA;
	  }
	else {
	  newder = (k1[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k2[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k3[matchar]  + ypv;
	  if(yav > yavA) yav = yavA;
	  if(yav < yavD) yav = yavD;
	}
	break;
      }
      case 1:{			// klesani zakl promene
	if(ypv == ypvD)
		{
		  yav = yavD;
		}
	else {
	  newder  = (k4[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k5[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k6[matchar]  + ypv;
	  if(yav < yavD) yav = yavD;
	  if(yav > yavA) yav = yavA;
	}
	break;
      }
      }
    outvalue = yav;
    outvalue2 = -10000.0;
    //if(ipp == 0)
    //	  fprintf(Outt,"cas t %lf ipp je %ld   xav %e  a ypv %e  A KAPA JE %e \n",t, ipp,xav,ypv, outvalue);
    
    break;
  }
 case 6:
   {
     xav = x;
     ypv = ineq1;		  
	    CorD(6,other,1,xav,yavA,yavD,other1);
	    CorD(6,other,1,xpv,ypvA,ypvD,other1);
	    
	    
	    der_value_hyst(matchar,40,xpv,derpvA,derpvD,ipp);
	    if(xav >= xpv)
	      {
			hyst = 0;
	      }
	    else{
			hyst = 1;
	    }
	    
	    switch(hyst)
	      {
	      case 0:{			// zvysovani promene....
			if(ypv == ypvA)                    
				{
					yav = yavA;
				}
			else {
					newder = (k1[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k2[matchar]*derpvA*(ypv - ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
					yav = (xav-xpv)*newder*k3[matchar]  + ypv;
					if(yav < yavA) yav = yavA;
					if(yav > yavD) yav = yavD;
				}
			break;
			}
	      case 1:{			// klesani zakl promene
			if(ypv == ypvD) 
				{
					yav = yavD;
				}
			else {
					newder  = (k4[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k5[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
					yav = (xav-xpv)*newder*k6[matchar] + ypv;
					if(yav < yavA) yav = yavA;
					if(yav > yavD) yav = yavD;
		  
				}
			break;
			}
	      }
	    
		der_value_hyst(matchar,40,xav,deravA,deravD,ipp);
	    
	    if (yav == yavA)
			newder = deravA;
	    if (yav == yavD)
			newder = deravD;
	    
	    double newder2;

		newder2 = (yav-ypv)/(xav-xpv);
		outvalue = yav;
	    outvalue2 = newder*1000;
		if (t>15559200)
		{
			if(t< 15724800)
			{
			  if(ipp == 80){
			    //fprintf(Outt3,"%lf %e  %e  %e  %e\n",t, xav,yavA,yav, yavD);
			  }
			}
			else{
			  if(t> 31536000 && t< 31622400)
			    {
			      if(ipp == 80){
				//fprintf(Outt3,"%lf %e  %e  %e  %e\n",t, xav,yavA,yav, yavD);
			      }
			    }
			}

		}
	    
	  }
	  break;
	}
    }
*/

}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Koudelka, 4.2.2014
*/
void kunmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_rel_humidity;
  dofname[1] = trf_temperature;
}



/**
   Function returns temperature in the given integration point.
   
   @param ipp - integration point id

   @return Funtion returns value of temperature stored in the integartion point.

   24. 10. 2013, JK
*/
double kunmat::give_temperature (long ipp)
{
  return Tm->ip[ipp].av[1];
}



/**  
  Function returns initial temperature in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of initial temperature stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double kunmat::give_inittemperature (long /*ipp*/)
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
double kunmat::give_rel_hum (long ipp)
{
  return Tm->ip[ipp].av[0];
}



/**  
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integartion point.

  Created by TKo, 25.11.2013
*/
double kunmat::give_vol_moist (long ipp)
{
  return Tm->ip[ipp].eqother[0];
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   23. 4. 2014, by TKo
*/
void kunmat::give_reqntq(long *antq)
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
