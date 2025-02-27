/*
    File:             kunzel.cpp
    Author:           Jiri Madera
    Purpose:          c
    Source:           Kunzel, 19xx;
    Assumptions:
    27.9.2003


    position CORD:  2 - density
                    3 - porosity
                    4 - faktor difusniho odporu
                    5 - kapa
                    6 - sorption izoterm
                    7 - saturated moisture
                    8 - none
                    9 - Cecko
                    10 - Lambda
                    11 - 13 - none
                    14 Dcoef
                    15 - binding isotherm
                    16 - cfmax
                    17 ws
                    18 - 19 - none
*/


#include "kunzel.h"
#include "globalt.h"


kunmat::kunmat()
{
  rho_m = 0.0;  //density
  rho_w = 1000.0;  //water density
}

kunmat::~kunmat()
{

  
}



/**
   function computes conductivity matrix of the material
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
void kunmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double phi,t;
  k = 0.0;
  
  phi = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(phi,t,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(phi,t,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(phi,t,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(phi,t,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void kunmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double phi,t;
  k = 0.0;
  
  phi = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(phi,t,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(phi,t,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(phi,t,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(phi,t,ipp);
  
  fillm(0.0,d);

 // fprintf(Outt,"t je %lf phi je %lf k je %lf pro ipp %d\n",t, phi,k, ipp);

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

void kunmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double phi,t;
  k = 0.0;
  
  phi = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(phi,t,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(phi,t,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(phi,t,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(phi,t,ipp);
  
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
*/
void kunmat::matcap (double &c,long ri,long ci,long ipp)
{
  double phi,t;
  c = 0.0;
  
  phi = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    //c = DerivativeOfTheMoistureRetentionCharacteristik(phi,t);
	c =  Tm->ip[ipp].eqother[1];
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = DerivaceHustotyEntalpiePodleTeploty(phi,t,ipp);
}


void kunmat::values_correction (vector &nv, long ipp)
{
  //  relative humidity
  relhum_check(nv[0],nv[1],ipp);
}


/**
   function checks if relative humidity is grater than 1.0
   @param rh - relative humidity
   @param t  - temperature
   @param ipp - number of integration point


   TKr, 5.8.2005
*/
void kunmat::relhum_check(double &rh,double &t,long ipp)
{
  if (rh >= 1.0) 
    rh = 1.0;
  if (rh <= 0.0) 
    rh = 0.0;
  

  if (t >= 350.0) 
    t = 350.0;
  if (t <= 240.0) 
    t = 240.0;

  //storing into gauss point
 // Tm->ip[ipp].av[0] = rh;
 // Tm->ip[ipp].av[1] = t;
}


//---------------------------------------------------------------------------
 /**
   popis

 */

double kunmat::DerivaceTlakuNasycenychParNaTeplote(double rh, double tk)
{

   return ((4042.9 * exp(23.5771-4042.9/(tk - 37.58)))/((tk-37.58)*(tk-37.58)));
}

//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat::PermeabilitaVodniPary(double rh, double tk,long ipp)
{
  double P, Rv, Pa;
  //double da, dp, mi, moist, smc;
  double da, dp, mi, smc;
  
  if (rh >1) rh = 1.0;
  if (rh < 0) rh = 0.0;
  
  //sorption_izoterms_giva_data(0,rh,tk,moist,a3);
  
  moist =  Tm->ip[ipp].eqother[0];
  smc = Tm->ip[ipp].eqother[2];
  //CorD(7,kd,0,rh,smc,a2,a3);
  CorD(4,kd,smc,moist,mi,a2,a3);
  
  Rv = 461.5;
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
/**
   popis

 */

double kunmat::tokJ2(double rh, double tk, long ipp)
{

double dp, dpsdt;


        if (rh >1) rh = 1;
        if (rh < 0) rh = 0;

       // nutno opravit znaceni ... vymenit dp za dpsdt
        dp = DerivaceTlakuNasycenychParNaTeplote(rh,tk);
        dpsdt = PermeabilitaVodniPary(rh, tk,ipp);

        J2 = dp * rh * dpsdt;
	
//	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,rh,tk);
        
	return (J2);

}

//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat::tokJ4(double rh, double tk, long ipp)
{
double lv1, lambdaakt, dpsdt, dp, moist,smc;


        if (rh >1) rh = 1;
        if (rh < 0) rh = 0;

       // sorption_izoterms_giva_data(0,rh,tk,moist,a3);
        
		moist =  Tm->ip[ipp].eqother[0];
		smc = Tm->ip[ipp].eqother[2];
	       //	CorD(7,kd,0,rh,smc,a2,a3);
        CorD(10,kd,smc,moist,lambdaakt,a2,a3);

        dpsdt =  DerivaceTlakuNasycenychParNaTeplote(rh, tk);
        dp =  PermeabilitaVodniPary(rh, tk,ipp);

        lv1 = LatentHeatofEvaporationOfWater (rh, tk);

        J4 = lambdaakt + lv1 * dp * rh * dpsdt;
//	fprintf(Outt,"tokJ4 je %3.2e pro rh= %lf a t= %lf \n",J4,rh,tk);
        return (J4);

}

//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat::tokJ3(double rh, double tk, long ipp)
{
double lv1;
double dp, ps;

       if (rh >1) rh = 1;
       if (rh < 0) rh = 0;

       dp = PermeabilitaVodniPary(rh, tk,ipp);
       ps = TlakNasycenychVodnichParNaTeplote (rh, tk);
       lv1 = LatentHeatofEvaporationOfWater (rh, tk);

       J3 = lv1 * dp * ps;
//	fprintf(Outt,"tokJ3 je  %3.2e pro rh= %lf a t= %lf \n",J3,rh,tk);
       return (J3);
}
 //---------------------------------------------------------------------------

double kunmat::LatentHeatofEvaporationOfWater(double rh, double tk)
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
double  kunmat::TlakNasycenychVodnichParNaTeplote(double rh, double tk)
{
     return (exp(23.5771 - 4042.9/(tk - 37.58)));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat::DerivaceHustotyEntalpiePodleTeploty(double rh, double tk, long ipp)
{
  double rom, c, moist, smc;
  if (rh >1) rh = 1;
  if (rh < 0) rh = 0;
  
 // sorption_izoterms_giva_data(0,rh,tk,moist,a3);
  moist =  Tm->ip[ipp].eqother[0];
  smc = Tm->ip[ipp].eqother[2];

  //CorD(7,kd,0,rh,smc,a2,a3);
  CorD(2,kd,0,0,rom,a2,a3);
   CorD(9,kd,smc,moist,c,a2,a3);
 rho_m = rom;

 return (rom * c);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat::HygroscopicMoisture(double rh, double tk)
{
// pro funkci ROOT
double mhmc, mhrh;
        CorD(6,kd,0,rh,mhmc,mhrh,a3);

        return (mhmc/(1-sqrt(1-mhrh)));

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat::DerivativeOfTheSorptionIsotherm(double rh, double tk)
{
// pro funkci ROOT
 double  mhv;
        if (rh >1) rh = 1;
        if (rh < 0) rh = 0;

        mhv = HygroscopicMoisture (rh,tk);

        return (1000 * mhv/(2*sqrt(1-rh)));

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double  kunmat::DerivativeOfTheRetentionCurve(double rh, double tk, long ipp)
{
// pro funkci ROOT
        double smc, mhmc, mhrh;
       // CorD(7,kd,0,rh,smc,a2,a3);
        smc = Tm->ip[ipp].eqother[2];
        CorD(6,kd,0,rh,mhmc,mhrh,a3);
       // mhrh = 0.95;

        return (1000 * (smc - mhmc)/(1 - mhrh));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double  kunmat::DerivativeOfTheMoistureRetentionCharacteristik(double rh,double tk, long ipp)
{

  double drfdrh;

  sorption_izoterms_giva_data(1,rh,tk,a1,drfdrh, ipp);

  return (drfdrh);

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat::tokJ1(double rh, double tk, long ipp)
{
//double P;
double kapa1,b, mi, dp, dmretc, ps, drh, moist, smc;
int s;


        if (rh < 0) rh = 0;
        if (rh >1) rh = 1;

       // sorption_izoterms_giva_data(0,rh,tk,moist,a3);
        moist =  Tm->ip[ipp].eqother[0];
		//CorD(7,kd,0,rh,smc,a2,a3);
        smc = Tm->ip[ipp].eqother[2];
        CorD(4,kd,smc,moist,mi,a2,a3);
        kapa1 =  Tm->ip[ipp].eqother[3];
     /*   CorD(5,s,smc,moist,kapa1,b,a3);

        switch (s){
	case 0:
	case 1:
	case 2:{
	  break;
	}
        case 30:{
	  kapa1 = kapa_exp(kapa1,b,rh, tk);
	  break;
        }
	case 32:
	case 33:{
	  break;
	}
        default:{
              fprintf (stderr,"\n\n unknown definition of kapa is required");
              fprintf (stderr,"\n in function tokJ1 (file %s, line %d)\n",__FILE__,__LINE__);
        }
        }*/

       // P = 101325;

        dp = PermeabilitaVodniPary(rh, tk,ipp);
        ps = TlakNasycenychVodnichParNaTeplote (rh, tk);
       // dmretc = DerivativeOfTheMoistureRetentionCharacteristik (rh, tk);

       //sorption_izoterms_giva_data(1,rh,tk,a1,dmretc);
       dmretc =  Tm->ip[ipp].eqother[1];

       drh = kapa1 * dmretc;

        J1 = drh + dp * ps;
	
//	fprintf(Outt,"Cas je %d sekund \n",Tp->time);
//	fprintf(Outt,"tokJ1 je %3.2e pro rh= %lf a t= %lf \n",J1,rh,tk);
        return (J1);
}
//---------------------------------------------------------------------------
void kunmat::read(XFILE *in)
{
  int i, j;
  
  for (i = 2; i<19; i++)
    {
      xfscanf(in, "%d",&MatChar[i]);
    }
  xfscanf(in, "%d",&MatChar[19]);
  for (i = 2;i <= 19; i++)
    {
      
      switch (MatChar[i]){
      case 0:{
	break;
      }
      case 1:{
	xfscanf(in, "%lf",&MatConst[i]);
	//edit_chechbox(i,1);
	break;
      }
      case 2:{
	
	
	if( MatChar[i-1] ==2)
	  {
	    xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
	  }
	else
	  {
	    xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
	  }
	
	switch (int(MatData[i][1][0])){
	case 2:{
	  for (j=1; j<=MatData[i][0][0];j++)
	    {
	      xfscanf(in, "%lf %lf",&MatData[i][0][j],&MatData[i][1][j]);
	    }
	  break;
	}
	case 3:{
	  for (j=1; j<=MatData[i][0][0];j++)
	    {
	      xfscanf(in, "%lf %lf %lf",&MatData[i][0][j],&MatData[i][1][j],&MatData[i][2][j]);
	    }
	  break;
	}
	default:{
	  fprintf (stderr,"\n\n unknown definition of Material KUNZEL-DATA  is required");
	  fprintf (stderr,"\n in function read (file %s, line %d)\n",__FILE__,__LINE__);
	}
	}
	xfscanf(in, "%d",&MatData[i][0][0]);
	//   edit_chechbox(i,2);
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
	fprintf (stderr,"\n\n unknown definition of Material KUNZEL is required");
	fprintf (stderr,"\n in function read (file %s, line %d)\n",__FILE__,__LINE__);
      }
      }
    }
  
  madripom = 0;
  CorD(2,kd,0,0,rho_m,a2,a3);
}


void kunmat::print(FILE *out)
{
  //fprintf (out,"\n kunzel material \n");
  int i, j;

  fprintf (out,"\n");
  for (i = 2; i<19; i++)
    {
      fprintf(out, " %d ",MatChar[i]);
    }
  fprintf(out, " %d ",MatChar[19]);
  for (i = 2;i <= 19; i++)
    {
      
      switch (MatChar[i]){
      case 0:{
	break;
      }
      case 1:{
	fprintf(out, "\n %e ",MatConst[i]);
	//edit_chechbox(i,1);
	break;
      }
      case 2:{
	if( MatChar[i-1] ==2)
	  {
	    fprintf(out, "\n %d %d ",int(MatData[i][0][0]),int(MatData[i][1][0]));
	  }
	else
	  {
	    fprintf(out, "\n %d %d ",int(MatData[i][0][0]),int(MatData[i][1][0]));
	  }
	
	switch (int(MatData[i][1][0])){
	case 2:{
	  for (j=1; j<=MatData[i][0][0];j++)
	    {
	      fprintf(out, "\n %e %e ",MatData[i][0][j],MatData[i][1][j] );
	    }
	  break;
	}
	case 3:{
	  for (j=1; j<=MatData[i][0][0];j++)
	    {
	      fprintf(out, "\n %e %e %e ",MatData[i][0][j],MatData[i][1][j],MatData[i][2][j] );
	    }
	  break;
	}
	default:{
	  fprintf (out,"\n");
	}
	}
	
	fprintf(out, "\n %d ",int(MatData[i][0][0]));
	//   edit_chechbox(i,2);
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
	fprintf (out,"\n");
      }
      }
    }
  fprintf (out,"\n");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat::soptionizothermDerivation(double rh, double tk)
{
double a, b, n;

     CorD(6,kd,0,rh,a,b,n);


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
    pgws = TlakNasycenychVodnichParNaTeplote(rh,t);
    trc = pgws;
    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = TlakNasycenychVodnichParNaTeplote(rh,t);
    trc = pgws;
    break;
  }
  case 32:{//water vapour pressure pgw
    trc = 0.0;
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

    pgws = TlakNasycenychVodnichParNaTeplote(rh,t);
    pgwspred = TlakNasycenychVodnichParNaTeplote(rh,t);

    nodval = bv*pgwspred;

    break;
  }
  case 31:{//pgw
    nodval = bv;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = TlakNasycenychVodnichParNaTeplote(rh,t);

    nodval = pgws*rh;
    bv = pgws*bv;
    nodval = bv - nodval;//inverse eq.

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
    pgws = TlakNasycenychVodnichParNaTeplote(rh, t);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;//inverse eq.

    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = TlakNasycenychVodnichParNaTeplote(rh, t);
    flux = pgws*rh;

    flux = bv - flux;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = TlakNasycenychVodnichParNaTeplote(rh, t);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;//inverse eq.

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
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
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
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
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
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
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
    sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    other = w;
    break;
  }
  case 3:{//water vapour pressure pgw
    other = rh * TlakNasycenychVodnichParNaTeplote (rh, t);
    break;
  }
  case 4:{//moisture content u = kg/kg (u.rho_m/rho_w = w)
    double w, a3;
    sorption_izoterms_giva_data(0,rh,t,w,a3, ipp);
    other = w*rho_w/rho_m;
    break;
  }    
  case 5:{
	  other = 0.0;
      break;
  }
  case 6:{//temperature
	  other = 0.0;
      break;
  }
  case 7:{//temperature
	  other = 0.0;
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
    //fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
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
   // fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
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

double kunmat::kapa_exp(double a, double b, double rh, double tk, long ipp)
{
 double moist;
 //if (rh <= 0) rh = 0;

 //sorption_izoterms_giva_data(0,rh,tk,moist,a3);
// moist =  Tm->ip[ipp].eqother[0];
 //return (a * exp(b * moist));
 return (a * exp(b * rh));

}
//---------------------------------------------------------------------------
/*position CORD:     2 - density
                    3 - porosity
                    4 - faktor difusniho odporu
                    5 - kapa
                    6 - sorption izoterm
                    7 - saturated moisture
                    8 - hydraulic conductivity
                    9 - Cecko
                    10 - Lambda
                    11 - water retention curve
                    12 -13 - none
                    14 Dcoef
                    15 - binding isotherm
                    16 - cfmax
                    17 ws
                    18 - 19 - none    */

void kunmat::CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2)
{

  switch (MatChar[cislochar]){
  case 0:{
    kvyhl = 0;
    break;
  }
  case 1:{
    y = MatConst[cislochar];
    kvyhl = 1;
    break;
  }
  case 2:{
    int i = 1;

    int psl, prad;
    prad  = int( MatData[cislochar][0][0]);
    psl = int (MatData[cislochar][1][0]);
    if (x <=0)
      {
	x = 0;
	y = MatData[cislochar][1][1];
	if (psl == 3)
	  {
	    z = MatData[cislochar][2][1];
	  }

      }
    else {
      if (x > in)
	{
	  x = in;
	  y = MatData[cislochar][1][prad];
	  if (psl == 3)
	    {
	      z = MatData[cislochar][2][prad];
	    }
	}
      else
	{
	  while (x > MatData[cislochar][0][i]){
	    i++;
	  }

	  y = MatData[cislochar][1][i-1] + ((x - MatData[cislochar][0][i-1])*(MatData[cislochar][1][i] - MatData[cislochar][1][i-1]))/(MatData[cislochar][0][i] - MatData[cislochar][0][i-1]);
	  if (psl == 3)
	    {
	      z = MatData[cislochar][2][i-1] + ((x - MatData[cislochar][0][i-1])*(MatData[cislochar][2][i] - MatData[cislochar][2][i-1]))/(MatData[cislochar][0][i] - MatData[cislochar][0][i-1]);
	    }
	}
    }
    kvyhl = 2;
	z2=MatData[cislochar][1][prad];

    break;
  }
  case 30:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    kvyhl = 30;
    break;
  }
  case 31:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    if (cislochar == 6)
      {
	z2 = MatFunce[cislochar][2];
      }
    kvyhl = 31;
    break;
  }
  case 32:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    kvyhl = 32;
    break;
  }
  case 33:{
    kvyhl = 33;
    //zatim nic
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown definition of material parameter is required");
    fprintf (stderr,"\n in function CorD (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
}

void kunmat::sorption_izoterms_giva_data(int kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp)
{
  /*
    kod ..... 0- jen vlhkost, 1 - jen derivaci...
  */
  int s;
  double hmrh, mhmc, smc, u, a, n;
  
  if (rh < 0) rh = 0;
  
  CorD(6,s,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN
  
  switch (s){
  case 2:{
    if(rh>0.976)
      {
	
	CorD(7,kd,0,rh,smc,a2,a3);
	CorD(6,kd,0,rh,mhmc,a2,mhmc);
	give_data(rh,mhmc, smc, 0.976,moist);
	dmoistdrh  =1000 * (smc - mhmc)/(1 - 0.976);
      }
    else
      {
	CorD(6,s,1,rh,moist,dmoistdrh,a3);
	dmoistdrh =1000* derivation_sorption_izoterm_data(rh,tk, 0);
      }
    break;
  }
  case 30:{
    CorD(6,s,1,rh,a1,hmrh,a3);
    CorD(7,s,0,rh,smc,a2,a3);
    CorD(6,s,0,rh,mhmc,hmrh,a3);
    give_data(rh,mhmc, smc, hmrh,moist);
    if (rh < hmrh)
      {
	dmoistdrh = DerivativeOfTheSorptionIsotherm (rh, tk);
      }
    else {
      dmoistdrh = DerivativeOfTheRetentionCurve(rh, tk, ipp);
    }
    
    break;
  }
  case 31:{
    switch (kod){
    case 0:{
      CorD(6,kod,0,0,u,a,n);
      moist = si_kk_hansen(rh, 0,u,a,n);
      dmoistdrh = soptionizothermDerivation(rh, tk);
      break;
    }
    case 1:{
      dmoistdrh = soptionizothermDerivation(rh, tk);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown definition of Sorption isotherm HANSEN is required");
      fprintf (stderr,"\n in function sorption_izoterms_giva_data (file %s, line %d)\n",__FILE__,__LINE__);
    }
    }
    break ;
  }
  default:{
    fprintf (stderr,"\n\n unknown definition of Sorption isotherm  is required");
    fprintf (stderr,"\n in function sorption_izoterms_giva_data (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  
  
}

double kunmat::si_kk_hansen(double rh, double tk,double u, double a, double n)
{
    return (u*pow((1-(log(rh))/a),(-1/n)));
}

double kunmat::derivation_sorption_izoterm_data(double x1, double x2, double x3)
{
  int i = 1;
  double derbi;
     if ( x1 < 0)
          {
          x1 = 0;
          }

     while (x1 > MatData[6][0][i]){
               i++;
               }
 derbi = (MatData [6][1][i] - MatData [6][1][i-1])/(MatData [6][0][i] - MatData [6][0][i-1]);   //TODO: Add your source code here
 return(derbi);
 }

/*void kunmat::aux_values (long ipp)
{
  double rh,tk;
  
  rh = Tm->ip[ipp].av[0];
  tk = Tm->ip[ipp].av[1];
  
  //  computation of moisture and derivative of moisture with respect to relative humidity
  sorption_izoterms_giva_data(0,rh,tk,moist,dmoistdrh);
  
  
}*/
void kunmat::aux_values (long ipp,double *in,double *inp, double *ineq,double *out)
{
  double x1,x2,x1p, x2p, ypv1,ypv2,ypv3,ypv4,ypv5, ypv0;
  
  double w,dwdf, kapaA;
  double t;
  t = Tp->time;
  x1 = in[0];
  x2 = in[1];
  
  x1p = inp[0];
  x2p = inp[1];
  

  ypv0 =ineq[0];
  ypv1 =ineq[1];
  ypv2 =ineq[2];
  ypv3 =ineq[3];
  ypv4 =ineq[4];
  ypv5 =ineq[5];

  sorption_izoterms_giva_data(0,x1,x2,w,dwdf, ipp);
  // sorption_izoterms_values(MatChar[6],ipp,x1,x1p,ypv1,w,dwdf);
  kapa_values(MatChar[5],ipp,w,ypv0,ypv2, kapaA);
  out[3] = kapaA;

  
  out[5] = -10000.0;
 
  out[0] = w;
  out[1] = dwdf;
  //out[3] = -10000.0;
  out[4] = -10000.0;
  out[2] = MatConst[7];

 /* if (ipp == 10)
  fprintf(Outt,"t= %lf je w == %3.2e dwdrh je %3.2e a kapa je %3.2e \n",t,w,dwdf,kapaA);
 */
}

void kunmat::sorption_izoterms_values(int kod, long ipp, double x1,double xpv, double ineq1, double & w, double & dwdf)
{
 
  //int s;
  //double hmrh, mhmc, smc, u, a, n;
  //double x2;

/* x2 = -10000.0;

       if (x1 < 0) x1 = 0;

      // CorD(6,s,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (kod){
       case 2:{
             CorD(6,s,1,x1,w,dwdf,a3);
             dwdf =derivation_sorption_izoterm_data(x1,x2,0,ipp)*1000;
             break;
       }
       case 30:{
                    CorD(7,s,0,x1,smc,a2,a3);
                    CorD(6,s,0,x1,mhmc,hmrh,a3);
                    give_data(x1,mhmc, smc, hmrh,w);
                    if (x1 < hmrh)
                        {
                        dwdf = DerivativeOfTheSorptionIsotherm (x1, x2);
                        }
                    else {
                        dwdf = DerivativeOfTheRetentionCurve(x1, x2);
                        }
             break;
       }
       case 31:{
                         CorD(6,kod,0,0,u,a,n);
                         w = si_kk_hansen(x1, 0,u,a,n);
                         dwdf = soptionizothermDerivation(x1, x2);
             break ;
	   }
	        
       default:{
             fprintf (stderr,"\n\n unknown definition of Sorption isotherm  is required");
             fprintf (stderr,"\n in function sorption_izoterms_giva_data (file %s, line %d)\n",__FILE__,__LINE__);
			   }
			   */
}

void kunmat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
  Tm->ip[ipp].eqother[5]=out[5];

}

void kunmat::give_values (long ipp,double *av, double *pv, double *eq)
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

    
}

void kunmat::kapa_values (int kod, long ipp,double x1, double xpv, double ineq1, double &kapa)
{
	double kapa1, kapa2, kapa3, x1w;
	double smc, a, n, x2;
	int kd1;
//	x2 = 0.0;

	CorD(7,kd1,1,x1,smc,a,n);

    if( x1 > smc)
    {
     x1 = smc;
    }

    switch (kod){
    case 0: kapa = 0.0;
	    break;
    case 1: CorD(5,kod,1,x1,kapa,kapa1,kapa2);
	    break;
    case 2: CorD(5,kod,1,x1,kapa,kapa1,kapa2); //CorD(5,kod,1,x1w,kapa,kapa1,kapa2);
	    break;
    case 30:{
	    CorD(5,kod,1,x1,kapa1,kapa2,kapa3);
	    kapa = kapa_exp(kapa1,kapa2,x1, x2, ipp);
	    break;
		  }
/*    case 40: hystereze2 (5,x1,xpv,ineq1,kapa,kapa1,ipp);
	    break;*/
  }
}

void kunmat::initvalues (long ipp,long ido)
{
  double x1,x2;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
}


