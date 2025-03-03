/*
    File:             kunzel2.cpp
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
                    14 - Dcoef
                    15 - binding isotherm
                    16 - cfmax
                    17 - ws
                    18 - 19 - none
*/


#include "kunzel2.h"
#include "globalt.h"


kunmat2::kunmat2()
{
  rho_m = 0.0;  //density
  rho_w = 1000.0;  //water density
}

kunmat2::~kunmat2()
{

  
}

/**
   JK, 11.4.2019
*/
double kunmat2::give_hum (long nn)
{
  long k;
  double h;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {h = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {h = 0.0;}
  if (k<0)   {h = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return h;
}
/**
   JK, 11.4.2019
*/
double kunmat2::give_temp (long nn)
{
  long k;
  double t;
  
  k=Gtt->give_dof(nn,1);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return t;
}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void kunmat2::matcond (matrix &d,long ri,long ci,long ipp)
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
void kunmat2::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(x1,x2,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(x1,x2,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void kunmat2::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(x1,x2,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(x1,x2,ipp);
  
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

void kunmat2::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = tokJ1(x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = tokJ2(x1,x2,ipp);
  if((ri == 1) && (ci == 0))
    k = tokJ3(x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = tokJ4(x1,x2,ipp);
  
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
void kunmat2::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2;
  c = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
  //  c = DerivativeOfTheMoistureRetentionCharacteristik(x1,x2,ipp);
   c =  Tm->ip[ipp].eqother[1];
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = DerivaceHustotyEntalpiePodleTeploty(x1,x2, ipp);
}


void kunmat2::values_correction (vector &nv, long ipp)
{
  //  relative humidity
  relhum_check(nv[0],nv[1],ipp);
}


/**
   function checks if relative humidity is grater than 1.0
   @param rh - relative humidity ..x1
   @param t  - temperature ..x2
   @param ipp - number of integration point


   TKr, 5.8.2005
*/
void kunmat2::relhum_check(double &x1,double /*x2*/,long ipp)
{
  if (x1 >= 1.0) 
    x1 = 1.0;
  if (x1 <= 0.0) 
    x1 = 0.0;
  
  //storing into gauss point
  Tm->ip[ipp].av[0] = x1;
}


//---------------------------------------------------------------------------
 /**
   popis

 */

double kunmat2::DerivaceTlakuNasycenychParNaTeplote(double /*x1*/, double x2)
{

   return ((4042.9 * exp(23.5771-4042.9/(x2 - 37.58)))/((x2-37.58)*(x2-37.58)));
}

//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat2::PermeabilitaVodniPary(double x1, double x2, long ipp)
{
       double P, Rv, Pa;
       //double da, dp, mi, moist, smc;
       double da, dp, mi, smc, moist1;

        if (x1 > 1) x1 = 1;
        if (x2 < 0) x1 = 0;

	//sorption_izoterms_giva_data(0,x1,x2,moist,a3);
       //CorD(7,kd,0,x1,smc,a2,a3);
       moist1 = Tm->ip[ipp].eqother[0];
       smc = Tm->ip[ipp].eqother[2];
       CorD(4,kd,smc,moist1,mi,a2,a3);

       Rv = 461.5;
       Pa = 101325;
       P = 101325;

       da = (2.306e-5 * Pa)/(Rv * x2 * P)*pow((x2/273.15),1.81);
       dp = da/mi;
       return (dp);
 }
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat2::tokJ2(double x1, double x2, long ipp)
{

double dp, dpsdt;


        if (x1 > 1) x1 = 1;
        if (x1 < 0) x1 = 0;

       // nutno opravit znaceni ... vymenit dp za dpsdt
        dp = DerivaceTlakuNasycenychParNaTeplote(x1,x2);
        dpsdt = PermeabilitaVodniPary(x1,x2, ipp);

        J2 = dp * x1 * dpsdt;
	
//	fprintf(Outt,"tokJ2 je %3.2e pro rh= %lf a t= %lf \n",J2,x1,x2);
        
	return (J2);

}

//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat2::tokJ4(double x1, double x2, long ipp)
{
double lv1, lambdaakt, dpsdt, dp, moist1,smc;


        if (x1 > 1) x1 = 1;
        if (x1 < 0) x1 = 0;

        //sorption_izoterms_giva_data(0,x1,x2,moist,a3,ipp);
        moist1 = Tm->ip[ipp].eqother[0];
        smc = Tm->ip[ipp].eqother[2];
        //CorD(7,kd,0,x1,smc,a2,a3);
         CorD(10,kd,smc,moist1,lambdaakt,a2,a3);

        dpsdt =  DerivaceTlakuNasycenychParNaTeplote(x1, x2);
        dp =  PermeabilitaVodniPary(x1, x2, ipp);

        lv1 = LatentHeatofEvaporationOfWater (x1, x2);

        J4 = lambdaakt + lv1 * dp * x1 * dpsdt;
//  fprintf(Outt,"tokJ4 je %3.2e pro rh= %lf a t= %lf vlhkost  %lf  labda  %lf   dp  %lf  Dpsdt  %lf\n",J4,x1,x2, moist1, lambdaakt, dp, dpsdt);
        return (J4);

}

//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat2::tokJ3(double x1, double x2, long ipp)
{
double lv1;
double dp, ps;

       if (x1 > 1) x1 = 1;
       if (x1 < 0) x1 = 0;

       dp = PermeabilitaVodniPary(x1, x2, ipp);
       ps = TlakNasycenychVodnichParNaTeplote (x1, x2);
       lv1 = LatentHeatofEvaporationOfWater (x1, x2);

       J3 = lv1 * dp * ps;
//	fprintf(Outt,"tokJ3 je  %3.2e pro rh= %lf a t= %lf \n",J3,x1,x2);
       return (J3);
}
 //---------------------------------------------------------------------------

double kunmat2::LatentHeatofEvaporationOfWater(double /*x1*/, double x2)
{
        return (2.5008e6)*pow((273.15/x2),(0.167+x2*3.67e-4));
}


//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/**
   Function calculates saturation water vapor pressure
   saturation water vapor pressure by xxxx
   tk ... temperature [K]

   
   @param tk - temperature [K]
*/
double  kunmat2::TlakNasycenychVodnichParNaTeplote(double /*x1*/, double x2)
{
     return (exp(23.5771 - 4042.9/(x2 - 37.58)));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat2::DerivaceHustotyEntalpiePodleTeploty(double x1, double /*x2*/, long ipp)
{
  double rom, c, moist1, smc;
  if (x1 > 1) x1 = 1;
  if (x1 < 0) x1 = 0;
  
  moist1 = Tm->ip[ipp].eqother[0];
  //sorption_izoterms_giva_data(0,x1,x2,moist,a3,ipp);
  //CorD(7,kd,0,x1,smc,a2,a3);
  smc = Tm->ip[ipp].eqother[2];
  CorD(2,kd,0,0,rom,a2,a3);
  
  rho_m = rom;
  
  CorD(9,kd,smc,moist1,c,a2,a3);
  
  return (rom * c);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat2::HygroscopicMoisture(double x1, double /*x2*/)
{
// pro funkci ROOT
double mhmc, mhrh;
        CorD(6,kd,0,x1,mhmc,mhrh,a3);

        return (mhmc/(1-sqrt(1-mhrh)));

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */


double kunmat2::DerivativeOfTheSorptionIsotherm(double x1, double x2)
{
// pro funkci ROOT
 double  mhv;
        if (x1 > 1) x1 = 1;
        if (x1 < 0) x1 = 0;

        mhv = HygroscopicMoisture (x1,x2);

        return (1000 * mhv/(2*sqrt(1-x1)));

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double  kunmat2::DerivativeOfTheRetentionCurve(double x1, double x2)
{
// pro funkci ROOT
        double smc, mhmc, mhrh;
        CorD(7,kd,0,x1,smc,a2,a3);
        CorD(6,kd,0,x2,mhmc,mhrh,a3);
       // mhrh = 0.95;

        return (1000 * (smc - mhmc)/(1 - mhrh));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double  kunmat2::DerivativeOfTheMoistureRetentionCharacteristik(double /*x1*/,double /*x2*/, long /*ipp*/)
{

  double drfdrh;
                drfdrh = 1e20;
//  sorption_izoterms_giva_data(1,x1,x2,a1,drfdrh,ipp);

  return (drfdrh);

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
   popis

 */

double kunmat2::tokJ1(double x1, double x2, long ipp)
{
//double P;
double kapa1, mi, dp, dmretc, ps, drh, moist1, smc;



        if (x1 < 0) x1 = 0;
        if (x1 > 1) x1 = 1;

      //  sorption_izoterms_giva_data(0,x1,x2,moist,a3,ipp);
	moist1 = Tm->ip[ipp].eqother[0];
        smc = Tm->ip[ipp].eqother[2];
        //CorD(7,kd,0,x1,smc,a2,a3);
        CorD(4,kd,smc,moist1,mi,a2,a3);


	//kapa_values(MatChar[5],ipp,kapa1,x1);
	 kapa1 = Tm->ip[ipp].eqother[3];
 

	dp = PermeabilitaVodniPary(x1, x2, ipp);
        ps = TlakNasycenychVodnichParNaTeplote (x1, x2);

	dmretc = Tm->ip[ipp].eqother[1];
        drh = kapa1 * dmretc;
        J1 = drh + dp * ps;
	
        return (J1);
}
//---------------------------------------------------------------------------
void kunmat2::read(XFILE *in)
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
					xfscanf(in, "%lf %lf",&MatData[i][0][j],&MatData[i][1][j] );
				}
				break;
				   }
			case 3:{
				for (j=1; j<=MatData[i][0][0];j++)
				{
					xfscanf(in, "%lf %lf %lf",&MatData[i][0][j],&MatData[i][1][j],&MatData[i][2][j] );
				}
				break;
				   }
			default:{
				fprintf (stderr,"\n\n unknown definition of Material KUNZEL-DATA  is required");
				fprintf (stderr,"\n in function read (file %s, line %d)\n",__FILE__,__LINE__);
					}
			}
			xfscanf(in, "%d",&MatData[i][0][0]);
			break;
			   }
		case 3:{
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
		case 40:{
			//nacitani hystreze doplnit
			xfscanf(in, "%lf %lf",&MatData[i][0][0],&Init[i]);
		    xfscanf(in, "%lf %lf %lf %lf %lf %lf",&k1[i], &k2[i], &k3[i],&k4[i], &k5[i], &k6[i]);
		    for (j=1; j<=MatData[i][0][0];j++)
			{
				xfscanf(in, "%lf %lf",&MatData[i][0][j],&MatData[i][1][j]);
			}
			xfscanf(in, "%lf",&MatData[i][0][0]);
			xfscanf(in, "%lf",&MatData[i][2][0]);
		    for (j=1; j<=MatData[i][2][0];j++)
			{
				xfscanf(in, "%lf %lf",&MatData[i][2][j] ,&MatData[i][3][j] );
			}
			xfscanf(in, "%lf",&MatData[i][2][0]);
			break;
				}
		default:{
            fprintf (stderr,"\n\n unknown definition of Material KUNZEL is required");
            fprintf (stderr,"\n in function read (file %s, line %d)\n",__FILE__,__LINE__);
				}
    }
   }
}

/*
void kunmat2::print(FILE *out)
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
}*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/**
   popis

 */

double kunmat2::soptionizothermDerivation(double x1, double /*x2*/)
{
double a, b, n;

     CorD(6,kd,0,x1,a,b,n);


     if (x1 > 1) x1 = 1;
     if (x1 < 0) x1 = 0;


     return (1000 * (a/(b*x1*n))*pow((1-(log(x1))/b),(-1-(1/n))));

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
double kunmat2::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,h,t;
  new_trc = 0.0;
  
  h = give_hum (nn);
  t = give_temp (nn);

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
double kunmat2::transmission_nodval(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval,h,t;
  new_nodval = 0.0;
  
  h = give_hum (nn);
  t = give_temp (nn);

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
double kunmat2::transmission_flux(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double flux,h,t;
  flux = 0.0;
  
  h = give_hum (nn);
  t = give_temp (nn);

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
double kunmat2::get_transmission_transcoeff_hh(double x1,double x2,long bc,long /*ipp*/)
{
  double trc,pgws;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity
    pgws = TlakNasycenychVodnichParNaTeplote(x1,x2);
    trc = pgws;
    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = TlakNasycenychVodnichParNaTeplote(x1,x2);
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
double kunmat2::get_transmission_nodval_hh(double bv,double x1,double x2,long bc,long /*ipp*/)
{
  double nodval,pgws, pgwspred;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity -> pgw

    pgws = TlakNasycenychVodnichParNaTeplote(x1,x2);
    pgwspred = TlakNasycenychVodnichParNaTeplote(x1,x2);

    nodval = bv*pgwspred;

    break;
  }
  case 31:{//pgw
    nodval = bv;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = TlakNasycenychVodnichParNaTeplote(x1,x2);

    nodval = pgws*x1;
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
double kunmat2::get_transmission_flux_hh(double bv,double x1,double x2,long bc,long /*ipp*/)
{
  double flux,pgws;

  switch (bc){//type of prescribed variable
  case 30:{//relative humidity -> pgw
    pgws = TlakNasycenychVodnichParNaTeplote(x1, x2);

    flux = pgws*x1;
    bv = pgws*bv;
    flux = bv - flux;//inverse eq.

    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = TlakNasycenychVodnichParNaTeplote(x1, x2);
    flux = pgws*x1;

    flux = bv - flux;
    break;
  }
  case 32:{//relative humidity -> pgw
    pgws = TlakNasycenychVodnichParNaTeplote(x1, x2);

    flux = pgws*x1;
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
double kunmat2::get_transmission_transcoeff_tt(double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double kunmat2::get_transmission_nodval_tt(double bv,double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double kunmat2::get_transmission_flux_tt(double bv,double /*x1*/,double x2,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    flux = (bv - x2);
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

double kunmat2::get_othervalue(long compother, double x1,double x2, long ipp)
{
  double other;
  double t;
  t = Tp->time;
 

  switch (compother){
  case 0:{//relative humidity
    other = x1;
    break;
  }
  case 1:{//temperature
    other = x2;
      break;
  }
  case 2:{//moisture content w = m3/m3
    //double w, a3;
    other = Tm->ip[ipp].eqother[0];
     // sorption_izoterms_giva_data(0,x1,x2,x1w,x1wder,ipp);
	//other = x1w;
/*	if(ipp > 38 && ipp < 43){
		double rh1;
		rh1 = Tm->ip[ipp].eqother[5];
		fprintf(Outt,"Cas t %lf ipp je %ld  w je  %3.15e pro rh= %3.15e \n  rhipp= %3.15e \n",t, ipp,w,x1, rh1);
	}*/
    //other= -100000.0;
    break;
  }
  case 3:{//water vapour pressure pgw
    other = x1 * TlakNasycenychVodnichParNaTeplote (x1, x2);
    break;
  }
  case 4:{//moisture content u = kg/kg (u.rho_m/rho_w = w)
    double w;
    w = Tm->ip[ipp].eqother[0];
    other = w*rho_w/rho_m;
    break;
  } 
  case 5:other = Tm->ip[ipp].eqother[2];
  case 6:
  case 7:
	  break;
	
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);

}

/**
     function prints names of all variables in nodes
     @param out - output rhle
     @param compother - number of other components
*/
void kunmat2::print_othervalue_name(FILE *out,long compother)
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
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


void kunmat2::give_data(double x1,double Mhmc, double Smc, double Mhrh, double & moistakt)
{
  double hvezdicka;
  hvezdicka = Mhmc/(1-sqrt(1-Mhrh));

  if (x1 < Mhrh) {
                  moistakt = (1-sqrt(1-x1)) * hvezdicka ;
                     }
                   else {
                      if (x1 > 1)
                        {
                         moistakt = Smc;
                        }
                      else
                        {
                      moistakt = Mhmc + (x1 - Mhrh )/(1-Mhrh )* (Smc - Mhmc );
                         }
                     }
   
}
//---------------------------------------------------------------------------

double kunmat2::kapa_exp(double a, double b, double x1w, double /*x2*/, long /*ipp*/)
{
  	if (x1w <= 0) x1w = 0;
 	return (a * exp(b * x1w));
}
//---------------------------------------------------------------------------

void kunmat2::CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2)
{
  int i;
  switch (MatChar[cislochar]){
  case 0:{
    kvyhl = 0;
    break;
  }
  case 1:{
    y = MatConst[cislochar];
	z = -10000;
	z2 = -10000;
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
    z2 = MatData[cislochar][1][prad];
    break;
  }
  case 30:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    kvyhl = 30;
	z2 = -10000;
    break;
  }
  case 31:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    if (cislochar == 6)
      {
	z2 = MatFunce[cislochar][2];
      }
	else {
		z2 = -10000;
	}
    kvyhl = 31;
    break;
  }
  case 32:{
    y = MatFunce[cislochar][0];
    z = MatFunce[cislochar][1];
    kvyhl = 32;
	z2 = -10000;
    break;
  }
  case 33:{
    kvyhl = 33;
    //zatim nic
	z2 = -10000;
	z = -10000;
	y = -10000;
    break;
  }
  case 40:{
	int pradA, pradD;
	pradA  = int(MatData[cislochar][0][0]);
	pradD  = int(MatData[cislochar][2][0]);
	if (x <=0)
		{
		x = 0;
		y = MatData[cislochar][1][1];
		z = MatData[cislochar][3][1];
		}
	else
	{
		if (x > in)
		{
			y = MatData[cislochar][1][pradA];
			z = MatData[cislochar][3][pradD];
		}
		else
		{
			i = 1;
			while (x > MatData[cislochar][0][i]){
				i++;
			}
			y = MatData[cislochar][1][i-1] + ((x - MatData[cislochar][0][i-1])*(MatData[cislochar][1][i] - MatData[cislochar][1][i-1]))/(MatData[cislochar][0][i] - MatData[cislochar][0][i-1]); 
			i = 1;
			while (x > MatData[cislochar][2][i]){
				i++;
			}
			z = MatData[cislochar][3][i-1] + ((x - MatData[cislochar][2][i-1])*(MatData[cislochar][3][i] - MatData[cislochar][3][i-1]))/(MatData[cislochar][2][i] - MatData[cislochar][2][i-1]);
		}
	}
	z2 = MatData[cislochar][1][pradA];
	kvyhl = 40;
		 }
	break;
  }
}

void kunmat2::sorption_izoterms_values(int kod, long ipp, double x1,double xpv, double ineq1, double & w, double & dwdf)
{
  int s;
  double hmrh, mhmc, smc, u, a, n;
  double x2;
  
  x2 = -10000.0;
  
  if (x1 < 0) x1 = 0;
  
  // CorD(6,s,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN
  
  switch (kod){
  case 2:{
    if(x1>0.976)
      {
	
	CorD(7,kd,0,x1,smc,a2,a3);
	CorD(6,kd,0,x1,a3,a2,mhmc);
	give_data(x1,mhmc, smc, 0.976,w);
	dwdf  =1000 * (smc - mhmc)/(1 - 0.976);
      }
    else
      {
	CorD(6,s,1,x1,w,dwdf,a3);
	dwdf =derivation_sorption_izoterm_data(x1,x2,0,ipp)*1000;
      }
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
	 
  case 40:{
    
    if(x1>=0.976)
      {
	CorD(7,kd,0,x1,smc,a2,a3);
	CorD(6,kd,0,x1,a1,a2,mhmc);
	give_data(x1,mhmc, smc, 0.976,w);
	dwdf  =1000 * (smc - mhmc)/(1 - 0.976);				
      }
    else{
      hystereze2 (6,x1, xpv, ineq1,w,dwdf, ipp);
    }
    break;
  }
    
  default:{
    fprintf (stderr,"\n\n unknown definition of Sorption isotherm  is required");
    fprintf (stderr,"\n in function sorption_izoterms_giva_data (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  
}

double kunmat2::si_kk_hansen(double x1, double /*x2*/,double u, double a, double n)
{
    return (u*pow((1-(log(x1))/a),(-1/n)));
}

double kunmat2::derivation_sorption_izoterm_data(double x1, double /*x2*/, double /*x3*/, long /*ipp*/)
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


void kunmat2::hystereze2 (int matchar, double x, double xpv, double ineq1,double & outvalue,double & outvalue2, long ipp)
{
  double xav,yav,newder,ypv,yavA,yavD,ypvA,ypvD,derpvA,derpvD,other1,deravA,deravD;
  int hyst, other;
   //double  k1,k2, k3, k4, k5, k6;
  double t;
  t = Tp->time;
   
   
  if (t < Init[matchar])
    {
      switch (matchar)
	{
	case 5: {
	  CorD(5,other,1,x,outvalue,outvalue2,other1);
	  
	}
	  break;
	case 6: {
	  CorD(6,other,1,x,outvalue,outvalue2,other1);
	  outvalue2 = derivation_dy_dx(6,x,0,1)*1000;
	}
	  break;
	}
    }
  else
    {
      switch(matchar)
	{
	case 5:{
	  double yavANorm,yavDNorm,ypvANorm,ypvDNorm,derpvANorm,derpvDNorm,ypvNorm,yavNorm;
	  double p, logP;
	  
	  xav = x; //chyba
	  //  kapa z predchoziho kroku
	  ypv = ineq1;
	  //  
	  ypvNorm = log(ineq1);
	  CorD(5,other,1,xav,yavA,yavD,other1);
	  CorD(5,other,1,xpv,ypvA,ypvD,other1);
	  // yavA - kappa pro navlhani pri aktualni vlhkosti
	  yavANorm = log(yavA);
	  //  yavD - kappa pro vysychani pri aktualni vlhkosti
	  yavDNorm = log(yavD);	
	  // ypvA - kappa pro navlhani pro predchozi vlhkosti
	  ypvANorm = log(ypvA);
	  //  ypvD - kappa pro vysychani pri predchozi vlhkosti
	  ypvDNorm = log(ypvD);	
	  
	  der_value_hyst(matchar,40,xpv,derpvA,derpvD,ipp); 
	  //  derpvA - derivace kapa pro navlhani
	  //  derpvD - derivace kapa pro vysychani
	  derpvANorm = log(derpvA);
	  derpvDNorm = log(derpvD);
	  // xav je aktualni vlhkost
	  //  xpv je predchozi vlhkost
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
	      if(ypvNorm == ypvANorm)
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
	  

	 /* if(xav > xpv)            
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
	  */
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
}

double kunmat2::derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2)

{
	int i = 1;
	double outvalue;
	
	while (prom > MatData[matchar][pomk1][i]){
			i++;
		}
	outvalue = (MatData [matchar][pomk2][i] - MatData [matchar][pomk2][i-1])/(MatData [matchar][pomk1][i] - MatData [matchar][pomk1][i-1]);
	
	return(outvalue);
}

void kunmat2::der_value_hyst (int matchar,int kod, double pv, double & outvalue,double & outvalue2, long /*ipp*/)
{

	switch (kod){
	case 0:{
		break;
		   }
	case 1:{
		break;
		   }
	case 2:{
		break;
		   }
	case 30:{
		break;
			}
	case 31:{
		break;
			}
	case 40:{

		outvalue = derivation_dy_dx(matchar,pv,0,1);
		outvalue2 = derivation_dy_dx(matchar,pv,2,3);

		}
		break;
	}
}

void kunmat2::kapa_values (int kod, long ipp,double x1, double xpv, double ineq1, double &kapa)
{
	double kapa1, kapa2, kapa3;
	double smc, a, n, x2;
	int kd1;
	x2 = 0.0;

	CorD(7,kd1,1,x1,smc,a,n);
  // x1w = Tm->ip[ipp].eqother[0];
    if( x1 > smc) x1 = smc;

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
    case 40: hystereze2 (5,x1,xpv,ineq1,kapa,kapa1,ipp);
	    break;
  }
}

/**
   function computes auxiliary values which are necessary for future computation
   
   JM, 29.5.2007
*/

void kunmat2::aux_values (long ipp,double *in,double *inp, double *ineq,double *out)
{
  double x1,x2,x1p, x2p, ypv1,ypv2,ypv3,ypv4,ypv5, ypv0;
  
  double kapaA,w,dwdf;
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

  sorption_izoterms_values(MatChar[6],ipp,x1,x1p,ypv0,w,dwdf); 
  
  kapa_values(MatChar[5],ipp,w,ypv0,ypv2, kapaA); 
  out[3] = kapaA;

  out[0] = w;
  out[1] = dwdf;
  out[2] = MatConst[7];
  out[4] = -10000.0;
  out[5] = -10000.0;

  if (ipp == 1)
  fprintf(Outt,"t= %lf je w == %3.2e dwdrh je %3.2e a kapa je %3.2e \n",t,w,dwdf,kapaA);
 
}

void kunmat2::initvalues (long ipp,long /*ido*/)
{
  double x1,x2;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
}

/**
   function selects auxiliary values
   
   JK, 7.1.2008
*/
void kunmat2::give_values (long ipp,double *av, double *pv, double *eq)
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


/**
   function saves auxiliary values
   
   JK, 7.1.2008
*/
void kunmat2::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
  Tm->ip[ipp].eqother[5]=out[5];

}
