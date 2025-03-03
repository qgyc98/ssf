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
		    18 - diffusion coefficient
		    19 - kunzeltype
*/


#include "kunzel.h"
#include "globalt.h"


kunmat::kunmat()
{
  long i;
  
  //  type of the model is not defiend
  kunzeltype = 0;
  
  //  density of solid material
  rho_m = 0.0;
  //  water density
  rho_w = 1000.0;
  //  water molecular weight
  M = 0.01801528;
  // universal gas constant
  R = 8.314472;
  
  for (i=0;i<20;i++){
    MatChar[i]=0;
  }
}

kunmat::~kunmat()
{

}

/**
   function reads material characteristics
   
   @param in - input file

   JM, 12.10.2011
*/
void kunmat::read(XFILE *in)
{
  long i,nr,nc;
  
  for (i = 2; i<20; i++){
    xfscanf(in, "%d",&MatChar[i]);
  }
  
  for (i = 2;i < 19; i++){
    
    switch (MatChar[i]){
    case 0:{
      //  material parameter is not present
      break;
    }
    case 1:{
      //  material parameter is given by a constant
      
      data[i]=new gfunct ();
      
      //  function is defined by a constant
      data[i]->tfunc=stat;
      
      xfscanf(in, "%lf",&data[i]->f);
      
      break;
    }
    case 2:{
      //  material parameter is given by a table
      
      //  nr - number of rows
      //  nc - number of columns, it must be 2
      xfscanf(in, "%ld %ld",&nr,&nc);
      
      //  function is defined by a table
      data[i]=new gfunct (tab,nr);
      //  piecewise linear approximation
      data[i]->tabf[0].itype = piecewiselin;
      //  reading of the table
      data[i]->tabf[0].readval (in);
      
      //  the number of rows is read again
      xfscanf(in, "%d",&nr);
      
      break;
    }
    case 3:{
      //edit_chechbox(i,3);
      break;
    }
    case 30:{
      if( MatChar[i-1] ==1)
	{
	  xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	}
      else
	{
	  xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	}
      break;
    }
    case 31:{
      if( MatChar[i-1] ==1)
	{
	  if (i == 6)
	    {
	      xfscanf(in, "%lf  %lf  %lf",&MatFunce[i][0],&MatFunce[i][1], &MatFunce[i][2]);
	    }
	  else
	    {
	      xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	}
      else
	{
	  if (i == 6)
	    {
	      xfscanf(in, "%lf  %lf  %lf",&MatFunce[i][0],&MatFunce[i][1], &MatFunce[i][2]);
	    }
	  else
	    {
	      xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	}
      
      break;
    }
    case 32:{
      if( MatChar[i-1] ==1)
	{
	  xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	}
      else
	{
	  xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	}
      break;
    }
    default:{
      print_err("unknown definition of Material KUNZEL is required",__FILE__,__LINE__,__func__);
    }
    }
  }

  //  determination of the type of Kunzel model
  kunzeltype = MatChar[19];
  
  CorD(2,kd,0,rho_m,a2,a3);
}


/**
   function prints material characteristics
   
   @param out - output file
   JM, 12.10.2011
*/
void kunmat::print(FILE *out)
{
  long i,j;
  
  fprintf (out,"\n");
  for (i = 2; i<20; i++)
    {
      fprintf(out, " %ld ",MatChar[i]);
    }
  fprintf(out, " %ld ",MatChar[19]);
  for (i = 2;i < 19; i++)
    {

      switch (MatChar[i]){
      case 0:{
	//  material parameter is not present
	break;
      }
      case 1:{
	//  material parameter is given by a constant
	fprintf(out, "\n %e ",data[i]->f);
	break;
      }
      case 2:{
	//  material parameter is given by a table
	
	fprintf(out, "\n %ld %d ",data[i]->tabf[0].asize,2);
	
	for (j=0;j<data[i]->tabf[0].asize;j++){
	  fprintf(out, "\n %e %e ",data[i]->tabf[0].x[j],data[i]->tabf[0].y[j]);
	}
	
	fprintf(out, "\n %ld ",data[i]->tabf[0].asize);
	break;
      }
      case 3:{
	//edit_chechbox(i,3);
	break;
      }
      case 30:{
	if( MatChar[i-1] ==1)
	  {
	    fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	  }
	else
	  {
	    fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	  }
	break;
      }
      case 31:{
	if( MatChar[i-1] ==1)
	  {
	    if (i == 6)
	      {
		fprintf(out, "\n %e  %e  %e ",MatFunce[i][0],MatFunce[i][1], MatFunce[i][2]);
	      }
	    else
	      {
		fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	      }
	  }
	else
	  {
	    if (i == 6)
	      {
		fprintf(out, "\n %e  %e  %e",MatFunce[i][0],MatFunce[i][1], MatFunce[i][2]);
	      }
	    else
	      {
		fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	      }
	  }
	
	break;
      }
      case 32:{
	if( MatChar[i-1] ==1)
	  {
	    fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	  }
	else
	  {
	    fprintf(out, "\n %e  %e ",MatFunce[i][0],MatFunce[i][1]);
	  }
	break;
      }
      default:{
	print_err("unknown definition of Material KUNZEL is required",__FILE__,__LINE__,__func__);
	fprintf (out,"\n");
      }
      }
    }
  fprintf (out,"\n");
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

 // fprintf(Outt,"t je %lf phi je %lf k je %lf pro ipp %d\n",t, phi,k, ipp);

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
}

/**
   function checks values of relative humidity and temperature
   if the values are out of range, the limit values are assigned
   
   @param nv - array with relative humidity and temperature
   
   JK, 7.10.2011
*/
void kunmat::values_correction (vector &nv)
{

 switch (kunzeltype) {
 case 1:{
    double rh, tk;
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
 case 4:
 {
     double pv, tk, pvs, rh;
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
 double pv, tk, pvs, rh;
    // partial pressure of water vapor
    pv=nv[0];
    //  temperature
    tk=nv[1];


    nv[0]=pv;
    nv[1]=tk;
 }
 break;
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

double kunmat::kmm(long ipp)
{
  double kk, tk;


 switch (kunzeltype) {
 case 1:{
    //double tk;
    double dw, mi, dp, dmretc, ps, dphi, moist, smc;

    // relative humidity
   // rh = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    //if (rh>1.0) rh = 1.0;
    //if (rh<0.0) rh = 0.0;

    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];

    // mi - water vapour diffusion resistance factor
    CorD(4,kd,moist,mi,a2,a3);

    // moisture diffusivity
    dw =  Tm->ip[ipp].eqother[3];

    dp = water_vapour_permeability(tk,ipp);
    ps = saturated_water_vapor_pressure( tk);

    //derivative of sorption isotherms
    dmretc =  Tm->ip[ipp].eqother[1];

    dphi = dw * dmretc;
    kk = dphi + dp * ps;

    break;
 }
 case 2:{
    double dw, dp;

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
    double D, moist;

    // partial pressure of water vapor
    //pv = Tm->ip[ipp].av[0];

    // temperature
    //tk = Tm->ip[ipp].av[1];

    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];

    // mi - water vapour diffusion resistance factor
    CorD(18,kd,moist,D,a2,a3);
    kk = D;
    break;
 }
 case 4:{
    double D, moist;

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
    double dp;

    // partial pressure of water vapor
    //pv = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    dp = water_vapour_permeability(tk,ipp);

     kk = dp;
    break;
 }
 case 7:{
    double moist, D, dp;

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
 default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }
 }
 //  tk = Tm->ip[ipp].av[1];

//fprintf(Outt,"Cas je %lf sekund \n",Tp->time);
//fprintf(Outt,"kmm je %3.2e pro  tk= %lf ipp = %ld\n",kk,tk,ipp);
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
  double kk, tk;

  switch (kunzeltype) {
  case 1:{
       double rh,tk;
       double dp, dpvs;

       // relative humidity
       rh = Tm->ip[ipp].av[0];

       // temperature
       tk = Tm->ip[ipp].av[1];

       if (rh>1.0) rh = 1.0;
       if (rh<0.0) rh = 0.0;

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
  case 7:{
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
 double kk,tk;

 switch (kunzeltype) {
 case 1:{
    double rh;
    double lv1, lambda, dpvs, dp, moist;

    // relative humidity
    rh = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    if (rh>1.0) rh = 1.0;
    if (rh<0.0) rh = 0.0;

    // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];
    // lambda - thermal conductivity
    CorD(10,kd,moist,lambda,a2,a3);

    dpvs =  derivative_of_saturated_water_vapor_pressure(tk);
    dp =  water_vapour_permeability(tk,ipp);

    lv1 = latent_heat_of_evaporation_of_water (tk);

    //lambda+ L_v delta_p phi dp_vs/dT
    kk = lambda + lv1 * dp * rh * dpvs;
     break;
 }
 case 2:
 case 3:
 case 4:
 case 5:
 case 6:{
    double  moist, lambda;
     // volumetric moisture content
    moist =  Tm->ip[ipp].eqother[0];

    // lambda - thermal conductivity
     CorD(10,kd,moist,lambda,a2,a3);

     kk = lambda;
    break;
 }
 case 7:{
    double rov;
    double lv1, lambda, pvs, dp, moist;

    // relative humidity
    rov = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    // volumetric moisture content
    moist =  0.01;
    // lambda - thermal conductivity
    CorD(10,kd,moist,lambda,a2,a3);
    dp =  water_vapour_permeability(tk,ipp);

    lv1 = latent_heat_of_evaporation_of_water (tk);

    //lambda+ L_v delta_p R Rov/M
    kk = lambda + lv1 * dp * R * rov/M;

    break;
 }
 default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }
 }
   //tk = Tm->ip[ipp].av[1];
 //fprintf(Outt,"Cas je %lf sekund \n",Tp->time);
//fprintf(Outt,"kht je %3.2e pro  tk= %lf ipp= %ld \n",kk,tk,ipp);
   //	fprintf(Outt,"tokJ4 je %3.2e pro rh= %lf a t= %lf \n",J4,rh,tk);
  return (kk);
}

//---------------------------------------------------------------------------

/**
   conductivity coefficient between moisture gradient and heat flux

   @param ipp - id of integration point

    JM, 12.10.2011

 */
double kunmat::khm(long ipp)
{
  double kk,tk;

  switch (kunzeltype) {
 case 1:{
    double  tk;

    double lv, dp, ps;
    // relative humidity
   // rh = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

   // if (rh>1.0) rh = 1.0;
   // if (rh<0.0) rh = 0.0;

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
    double tk;
    double lv, dp, ps;

    // temperature
    tk = Tm->ip[ipp].av[1];

    dp = water_vapour_permeability(tk,ipp);
    lv = latent_heat_of_evaporation_of_water (tk);

    kk = lv*dp;

    break;
 }
 case 6:{
    //double tk;
    double lv, dp, ps;

    // partial pressure of water vapor
   // pv = Tm->ip[ipp].av[0];

    // temperature
    tk = Tm->ip[ipp].av[1];

    dp = water_vapour_permeability(tk,ipp);
    lv = latent_heat_of_evaporation_of_water (tk);

    kk = lv*dp;

    break;
 }
 case 7:{
    double  tk, moist,dp;

    double lv,D;
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
 default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }
 }
  //tk = Tm->ip[ipp].av[1];
// fprintf(Outt,"Cas je %lf sekund \n",Tp->time);
//fprintf(Outt,"khm je %3.2e pro  tk= %lf ipp = %ld\n",kk,tk,ipp);

  //	fprintf(Outt,"tokJ3 je  %3.2e pro rh= %lf a t= %lf \n",J3,rh,tk);
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

  switch (kunzeltype) {
  case 1:{
       // derivative sorptin isotherms * density of water
       cc = Tm->ip[ipp].eqother[1];
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
      double moist, pi;

      // volumetric moisture content
     moist =  Tm->ip[ipp].eqother[0];

    // pi - porosity
     CorD(3,kd,moist,pi,a2,a3);

      cc = (pi-moist);

      break;
  }
   case 4:{
      double moist, pi;
      // partial pressure of water vapor
     // pv = Tm->ip[ipp].av[0];
      // temperature
      tk = Tm->ip[ipp].av[1];

      // volumetric moisture content
     moist =  Tm->ip[ipp].eqother[0];

    // pi - porosity
     CorD(3,kd,moist,pi,a2,a3);

      cc = (pi-moist)*M/(R*tk);;

      break;
  }
  case 5:{
      double moist, pi;
      // partial pressure of water vapor
     // pv = Tm->ip[ipp].av[0];
      // temperature
      tk = Tm->ip[ipp].av[1];

    // pi - porosity
     CorD(3,kd,1,pi,a2,a3);

      cc = (pi)*M/(R*tk);;

      break;
  }
  case 6:{
      double moist, pi;
      // partial pressure of water vapor
     // pv = Tm->ip[ipp].av[0];
      // temperature
      tk = Tm->ip[ipp].av[1];

    // pi - porosity
     CorD(3,kd,1,pi,a2,a3);

      cc = (pi)*M/(R*tk);

      break;
  }
  case 7:{
      double moist, pi;
      // partial pressure of water vapor
     // pv = Tm->ip[ipp].av[0];
      // temperature
     // tk = Tm->ip[ipp].av[1];

    // pi - porosity
     CorD(3,kd,1,pi,a2,a3);

      cc = pi;
      ///testovani
      cc = 1.0;

      break;
  }
  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }


 }
//   tk = Tm->ip[ipp].av[1];
//  fprintf(Outt,"Cas je %lf sekund \n",Tp->time);
//  fprintf(Outt,"cmm je %3.2e pro  tk= %lf ipp=%ld\n",cc,tk,ipp);
  //	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);

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
  double cc,tk ;

  switch (kunzeltype) {
  case 1:{
     //  double rh,tk;

       // relative humidity
     //  rh = Tm->ip[ipp].av[0];

       // temperature
     //  tk = Tm->ip[ipp].av[1];

       cc = 0.0;

      break;
  }
  case 2:{
      double pv,tk, dw;
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
  case 7:{
      //double pv,tk, dw;
      // partial pressure of water vapor
      //pv = Tm->ip[ipp].av[0];
      // temperature
      //tk = Tm->ip[ipp].av[1];

      cc = 0.0;
      break;
  }

  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
  }

  }
 // fprintf(Outt,"Cas je %d sekund \n",Tp->time);
 // fprintf(Outt,"cmt je %3.2e pro  tk= %lf \n",cc,tk);

  //	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);

  return (cc);
}

//---------------------------------------------------------------------------
/**
   capacity coefficient between moisture gradient and heat flux


   @param ipp - id of integration point

    JM, 12.10.2011
*/
double kunmat::chm(long ipp)
{
  double cc,tk;

  switch (kunzeltype) {
  case 1:{
    //   double rh,tk;

       // relative humidity
      // rh = Tm->ip[ipp].av[0];

       // temperature
      // tk = Tm->ip[ipp].av[1];

       cc = 0.0;

      break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:{
    //  double pv,tk;
      // partial pressure of water vapor
    //  pv = Tm->ip[ipp].av[0];
      // temperature
    //  tk = Tm->ip[ipp].av[1];

      cc = 0.0;

      break;
  }

  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);

  }

  }

  //	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);
 // fprintf(Outt,"Cas je %d sekund \n",Tp->time);
 // fprintf(Outt,"chm je %3.2e pro  tk= %lf \n",cc,tk);

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
  double cc,tk;

  switch (kunzeltype) {
  case 1:{
    //   double rh,tk;

       // relative humidity
     //  rh = Tm->ip[ipp].av[0];

       // temperature
    //   tk = Tm->ip[ipp].av[1];

     //  if (rh>1.0) rh = 1.0;
     //  if (rh<0.0) rh = 0.0;

       cc = derivative_of_the_enthalpy_density(ipp);

      break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:{
     // double pv,tk;
      // partial pressure of water vapor
    //  pv = Tm->ip[ipp].av[0];
      // temperature
     // tk = Tm->ip[ipp].av[1];

      cc = derivative_of_the_enthalpy_density(ipp);

      break;
  }

  default:{
      print_err("unknown type of Kunzel is required",__FILE__,__LINE__,__func__);
  }

  }

  //	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);
 // fprintf(Outt,"Cas je %lf sekund \n",Tp->time);
//   tk = Tm->ip[ipp].av[1];
 //fprintf(Outt,"cht je %3.2e pro  tk= %lf ipp = %ld\n",cc,tk,ipp);

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
double kunmat::water_vapour_permeability( double tk,long ipp)
{
  double P, Rv, Pa;
  double da, dp, mi;

    // volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];

  //mi -  water vapour diffusion resistance factor
  CorD(4,kd,moist,mi,a2,a3);
  // Rv = R/M , where R - universal gas constant, M - water molecular weight
  Rv = R/M;

  // air pressure
  Pa = 101325;
  P = 101325;

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

double kunmat::derivative_of_the_enthalpy_density(long ipp)
{
  double rom, c, moist;

 // volumetric moisture content
  moist =  Tm->ip[ipp].eqother[0];
  
  // rom - density of solid material
  CorD(2,kd,0,rom,a2,a3);

  // c- specific heat capacity
  CorD(9,kd,moist,c,a2,a3);

  rho_m = rom;

  return (rom * c);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat::hygroscopic_moisture(double rh, double tk)
{
  // pro funkci ROOT
  double mhmc, mhrh,a3;
  CorD(6,kd,rh,mhmc,mhrh,a3);
  
  return (mhmc/(1-sqrt(1-mhrh)));
  
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat::derivative_of_the_sorption_isotherm_root(double rh, double tk)
{
// pro funkci ROOT
 double  mhv;
        if (rh >1) rh = 1;
        if (rh < 0) rh = 0;

        mhv = hygroscopic_moisture (rh,tk);

        return (1000 * mhv/(2*sqrt(1-rh)));

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double  kunmat::derivative_of_retention_curve_root(double rh, double tk, long ipp)
{
// pro funkci ROOT
        double smc, mhmc, mhrh;
       // CorD(7,kd,0,rh,smc,a2,a3);
        smc = Tm->ip[ipp].eqother[2];
        CorD(6,kd,rh,mhmc,mhrh,a3);
       // mhrh = 0.95;

        return (1000 * (smc - mhmc)/(1 - mhrh));
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat::sorption_izotherm_derivative_kk_hansen(double rh, double tk)
{
  double a, b, n;
  
  CorD(6,kd,rh,a,b,n);
  
  
  if (rh >1) rh = 1;
  if (rh < 0) rh = 0;
  
  
  return (1000 * (a/(b*rh*n))*pow((1-(log(rh))/b),(-1-(1/n))));
  
}




/**
   function computes new transmission coefrhcient (for transmission_vector)
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefrhcient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of rhrst integration point on element
*/
double kunmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_trc,h,t;
  new_trc = 0.0;

  k=Gtt->give_dof(nn,0);
  if (k>0)   {h = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {h = 0.0;}
  if (k<0)   {h = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval();}

  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_hh(h,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;

  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
    new_trc = get_transmission_transcoeff_tt(h,t,bc,ipp);

  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector)
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefrhcient on the boundary,
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of rhrst integration point on element
*/
double kunmat::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_nodval,h,t;
  new_nodval = 0.0;

  k=Gtt->give_dof(nn,0);
  if (k>0)   {h = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {h = 0.0;}
  if (k<0)   {h = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval();}

  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_hh(nodval,h,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;

  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
    new_nodval = get_transmission_nodval_tt(nodval,h,t,bc,ipp);

  return (new_nodval);
}


/**
   function computes flux (for transmission_vector)
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefrhcient on the boundary,
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of rhrst integration point on element
*/
double kunmat::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double flux,h,t;
  flux = 0.0;

  k=Gtt->give_dof(nn,0);
  if (k>0)   {h = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {h = 0.0;}
  if (k<0)   {h = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval();}

  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_hh(nodval,h,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    flux = 0.0;

  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
    flux = get_transmission_flux_tt(nodval,h,t,bc,ipp);

  return (flux);
}


/**
   function creates transfer coefrhcient on the boundary for prescribed condition (relative humidity)

   @param rh - relative humidity
   @param t - temperature
   @param bc - type of boundary condition
*/
double kunmat::get_transmission_transcoeff_hh(double rh,double t,long bc,long ipp)
{
  double trc,pgws;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity
    pgws = saturated_water_vapor_pressure(t);
    trc = pgws;
    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = saturated_water_vapor_pressure(t);
    trc = pgws;
    break;
  }
  case 32:{//water vapour pressure pgw
    trc = 0.0;
    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (relative humidity)

   @param bv - prescribed value near the boundary
   @param rh - actual replative humidity on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double kunmat::get_transmission_nodval_hh(double bv,double rh,double t,long bc,long ipp)
{
  double nodval,pgws, pgwspred;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity -> pgw

    pgws = saturated_water_vapor_pressure(t);
    pgwspred = saturated_water_vapor_pressure(t);

    nodval = bv*pgwspred;

    break;
  }
  case 31:{//pgw
    nodval = bv;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = saturated_water_vapor_pressure(t);

    nodval = pgws*rh;
    bv = pgws*bv;
    nodval = bv - nodval;//inverse eq.

    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates flux on the boundary (convective mass transfer)

   @param bv - prescribed value near the boundary
   @param rh - rel. hum.
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double kunmat::get_transmission_flux_hh(double bv,double rh,double t,long bc,long ipp)
{
  double flux,pgws;

  switch (bc){//type of prescribed variable
  case 30:{//relative humidity -> pgw
    pgws = saturated_water_vapor_pressure( t);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;//inverse eq.

    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = saturated_water_vapor_pressure(t);
    flux = pgws*rh;

    flux = bv - flux;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = saturated_water_vapor_pressure(t);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;//inverse eq.

    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function creates transfer coefrhcient on the boundary for prescribed condition (temperature)

   @param h   - rel. hum
   @param t   - temperature
   @param bc  - type of boundary condition
*/
double kunmat::get_transmission_transcoeff_tt(double h,double t,long bc,long ipp)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    trc = 1.0;
    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }

  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (temperature)

   @param bv - prescribed value near the boundary
   @param h - actual rel. hum. on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double kunmat::get_transmission_nodval_tt(double bv,double h,double t,long bc,long ipp)
{
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    nodval = bv;
    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates heat flux on the boundary

   @param bv - prescribed value near the boundary
   @param h - rel. hum.
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double kunmat::get_transmission_flux_tt(double bv,double h,double t,long bc,long ipp)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    flux = (bv - t);
    break;
  }
  default:{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
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

double kunmat::get_othervalue(long compother,double rh,double t, long ipp)
{
  double other;

  switch (compother){
  case 0:{//relative humidity
    other = rh;
    break;
  }
  case 1:{//temperature
    other = t;
      break;
  }
  case 2:{//moisture content w = m3/m3
    double w, a3;
    //sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    other = w;
    break;
  }
  case 3:{//water vapour pressure pgw
    other = rh * saturated_water_vapor_pressure (t);
    break;
  }
  case 4:{//moisture content u = kg/kg (u.rho_m/rho_w = w)
    double w, a3;
  //  sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    other = w*rho_w/rho_m;
    break;
  }    
  case 5:{
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


void kunmat::give_data(double rh,double Mhmc, double Smc, double Mhrh, double & moistakt)
{
  double hvezdicka;
  hvezdicka = Mhmc/(1-sqrt(1-Mhrh));

  if (rh < Mhrh) 
  {                  
    moistakt = (1-sqrt(1-rh)) * hvezdicka ;
  }
  else {
    if (rh > 1)
      {
	moistakt = Smc;
      }
    else
      {
	moistakt = Mhmc + (rh - Mhrh )/(1-Mhrh )* (Smc - Mhmc );
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

void kunmat::sorption_izoterms_giva_data(long kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp)
{
  /*
    kod ..... 0- jen vlhkost, 1 - jen derivaci...
  */
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
      }
    else
      {
	CorD(6,s,rh,moist,dmoistdrh,a3);
	dmoistdrh =1000* derivative_sorption_izoterm_data(rh);
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
	dmoistdrh = derivative_of_the_sorption_isotherm_root (rh, tk);
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
      moist = si_kk_hansen(rh, 0,u,a,n);
      dmoistdrh = sorption_izotherm_derivative_kk_hansen(rh, tk);
      break;
    }
    case 1:{
      dmoistdrh = sorption_izotherm_derivative_kk_hansen(rh, tk);
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

double kunmat::si_kk_hansen(double rh, double tk,double u, double a, double n)
{
    return (u*pow((1-(log(rh))/a),(-1/n)));
}

/**
   function computes derivative of the sorption isotherm

   @param x1 - moisture content

   7.10.2011
*/
double kunmat::derivative_sorption_izoterm_data(double x1)
{
  long charid;
  double derbi;
  
  //  sorption isotherm
  charid=6;
  derbi = data[charid]->getderiv (x1);

  return derbi;
}

/**
   function computes auxilary values to eqother arrays in integration points

   @param ipp - integration point id
   @param out - array with values

   JM, 12.10.2011
*/
void kunmat::aux_values (long ipp,double *in,double *inp, double *ineq,double *out)
{


  switch (kunzeltype) {
  case 1:{
       double rh,tk;
       double w,dwdf, kapaA, smc;

       // relative humidity
       rh = in[0];
       // temperature
       tk = in[1];
       // saturated volumetric moisture content
       smc = data[7]->getval(0.0);
       // w - volumetric moisture content
       sorption_izoterms_giva_data(0,rh,tk,w,dwdf, ipp);
       if (w>smc) w = smc;
       if(w<0) w = 0.0;

       kapa_values(MatChar[5],ipp,w, kapaA);

       // moisture diffusivity
       out[3] = kapaA;
       // volumetric moisture content
       out[0] = w;
       // derivative of sorption isoterms
       out[1] = dwdf;
       // saturated volumetric moisture content
       out[2] = smc;
       // moisture diffusivity
       out[4] = kunzeltype;

       break;
  }
  case 2:
  case 3:
  case 4:{
	double pv,tk, pvs;
	double w,dwdf, kapaA, smc, rh;
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

       // moisture diffusivity
      out[3] = kapaA;
       // volumetric moisture content
      out[0] = w;
       // derivative of sorption isoterms

      out[1] = rh;
       // saturated volumetric moisture content
      out[2] = smc;

      out[4] = kunzeltype;
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

      // moisture diffusivity
      out[3] = -10000.0;
       // volumetric moisture content
      out[0] =smc;
       // derivative of sorption isoterms

      out[1] = 0.1;
       // saturated volumetric moisture content
      out[2] = smc;

      out[4] = kunzeltype;
      break;
  }
  default:
      ;
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
void kunmat::give_values (long ipp,double *av,double *eq)
{
  av[0] = Tm->ip[ipp].av[0];
  av[1] = Tm->ip[ipp].av[1];

  eq[0] = Tm->ip[ipp].eqother[0];
  eq[1] = Tm->ip[ipp].eqother[1];
  eq[2] = Tm->ip[ipp].eqother[2];
  eq[3] = Tm->ip[ipp].eqother[3];
  eq[4] = Tm->ip[ipp].eqother[4];
}
/* function computes moisture diffusivity



   JM, 12.10.2011
*/
void kunmat::kapa_values (long kod, long ipp,double x1, double &kapa)
{
  double kapa1, kapa2, kapa3;
  long kd1;

  switch (kod){
  case 0: kapa = 0.0;
    break;
  case 1: CorD(5,kod,x1,kapa,kapa1,kapa2);
    break;
  case 2: CorD(5,kod,x1,kapa,kapa1,kapa2); //CorD(5,kod,1,x1w,kapa,kapa1,kapa2);
    break;
  case 30:{
    CorD(5,kod,x1,kapa1,kapa2,kapa3);
    kapa = kapa_exp(kapa1,kapa2,x1);
    break;
  }
  default:{
    print_err("unknown kod is required",__FILE__,__LINE__,__func__);
  }
  }
}

void kunmat::initvalues (long ipp,long ido)
{
  double x1,x2;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
}


