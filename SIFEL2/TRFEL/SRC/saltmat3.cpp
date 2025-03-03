#include "globalt.h"
#include "saltmat3.h"
#include "globmatt.h"



kunmat *kuns;
grunewaldmat *grun;
/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
***********************************************/
/* x1 (w) is the volumetric water content [m3/m3]
   x2  (Cf) is the concentration of free ions in water [kg/m3]
   x3  (Cfyz)
   D the ion diffusivity [m2/s]
   Kapa = the moisture diffusivity [m2/s]

 position CORD:     2 - density
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
FILE *in1,*in2,*in3;
saltmat3::saltmat3 ()
{
  mw = 0.01801528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8.31441; //universal gas constant J.mol-1.K-1

}

saltmat3::~saltmat3 ()
{}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void saltmat3::matcond (matrix &d,long ri,long ci,long ipp)
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
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void saltmat3::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    //matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    //matcond3d (d,ri,ci,ipp);//3D
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
void saltmat3::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
  
  //  check of salt concentration
  //CorD(16,kd,0,x1,cfmax,a2,a3);
  cfmax=Tm->ip[ipp].eqother[3];

  if (x2 > cfmax){
    x2 = cfmax;
  }
  
    
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3,ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3,ipp);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3,ipp);
  
  d[0][0] = k;
}

void saltmat3::matcond2d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  d[0][0] = 0.0;  d[0][1] = 0.0;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void saltmat3::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  //  check of salt concentration
  //CorD(16,kd,0,x1,cfmax,a2,a3);
  cfmax=Tm->ip[ipp].eqother[3];

  if (x2 > cfmax){
    x2 = cfmax;
  }

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3,ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3,ipp);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3,ipp);

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

void saltmat3::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  //  check of salt concentration
  //CorD(16,kd,0,x1,cfmax,a2,a3);
  cfmax=Tm->ip[ipp].eqother[3];

  if (x2 > cfmax){
    x2 = cfmax;
  }

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3,ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3,ipp);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3,ipp);
  
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
void saltmat3::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2,x3;
  c = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    c = c11 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 1))
    c = c12 (x1,x2,x3,ipp);
  if((ri == 0) && (ci == 2))
    c = c13 (x1,x2,x3,ipp);

  if((ri == 1) && (ci == 0))
    c = c21 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 1))
    c = c22 (x1,x2,x3,ipp);
  if((ri == 1) && (ci == 2))
    c = c23 (x1,x2,x3,ipp);

  if((ri == 2) && (ci == 0))
    c = c31 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 1))
    c = c32 (x1,x2,x3,ipp);
  if((ri == 2) && (ci == 2))
    c = c33 (x1,x2,x3,ipp);
}


void saltmat3::values_correction (vector &nv,long ipp)
{
  double cfmax;
  
  cfmax=Tm->ip[ipp].eqother[3];

  if (nv[1] > cfmax){
    nv[1] = cfmax;
  }
}



/**
   function reads data and material parameters

   @param in  - input file
*/
void saltmat3::read (XFILE *in)
{
  int i,j;
  for (i = 2; i<19; i++)
    {
      xfscanf(in, "%d",&MatChar[i]);
    }
  xfscanf(in, "%d",&MatChar[19]);
  for (i = 2;i <= 19; i++)
    {
      
      switch (MatChar[i])
	{
	case 0:
	  break;
	case 1:xfscanf(in, "%lf",&MatConst[i]);
	  //edit_chechbox(i,1);
	  break;
	case 2:   if( MatChar[i-1] ==2)
			  {
				  xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
			  }
		else
		{
			xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
		}
		switch (int(MatData[i][1][0])){
			case 2: for (j=1; j<=MatData[i][0][0];j++)
					{
						xfscanf(in, "%lf %lf",&MatData[i][0][j],&MatData[i][1][j] );
					}
				break;
			case 3: for (j=1; j<=MatData[i][0][0];j++)
					{
						xfscanf(in, "%lf %lf %lf",&MatData[i][0][j],&MatData[i][1][j],&MatData[i][2][j] );
					}
				break;
		}
		xfscanf(in, "%d",&MatData[i][0][0]);
		break;
	case 3://edit_chechbox(i,3);
	  break;
	case 30:   if( MatChar[i-1] ==1)
	    {
	      xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	  else
	    {
	      xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	  
	  break;
	case 31: if( MatChar[i-1] ==1)
	    {
	      if (i == 6)
		{
		  xfscanf(in, "%lf %lf %lf",&MatFunce[i][0],&MatFunce[i][1], &MatFunce[i][2]);
		}
	      else
		{
		  xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
		}
	    }
	  else
	    {
	      if (i == 6)
		{
		  xfscanf(in, "%lf %lf %lf",&MatFunce[i][0],&MatFunce[i][1], &MatFunce[i][2]);
		}
	      else
		{
		  xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
		}
	    }
	  
	  break;
	case 32:  if( MatChar[i-1] ==1){
		xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
			  }
		else
		{
			xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
		}
		break;
	case 40:{
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
	  }
  }
  
}


/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

double saltmat3::k11(double x1,double /*x2*/,double /*x3*/,long ipp)
{
  double k11;
  double kapak, ps, delta, dfdw, w_hyg ,rh_hyg, hvezda;
  double p = 101325;
  double a,u,n;
  int kod;
  double smc;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

  kapak=Tm->ip[ipp].eqother[0];
 
  // ps = pgws(t);
  ps = pgws(294.15);
  //delta = permeabilitavodnipary(t,p);
  delta = permeabilitavodnipary(294.15,p);
  CorD(6,kod,1,x1,u,a,n);

  switch (kod)
  {
  case 2:   // CorD(6,kod,1,0,x1,fi,a2,a3);
			 dfdw = derivation_dy_dx(6,x1,0,1);
	  break;
  case 30:
			CorD(6,kod,0,x1,w_hyg,rh_hyg,a3);
            hvezda =grun ->sortpion_isotherm_root_shifted(x1,w_hyg ,rh_hyg);
            CorD(7,kod,0,x1,smc,a2,a3);
		//	grun->give_data_si_root_fi(w,0.0,0.0,w_hyg, smc, rh_hyg,hvezda,relhum);
			grun->give_data_si_root_dfidw(x1,0.0,0.0,w_hyg, smc, rh_hyg,hvezda,dfdw);
	  break;
  case 31:  {
	  dfdw = (u/(a*x1*n))*pow((1-(log(x1))/a),(-1-(1/n)));
			}
	  break;
  }
  
  k11 = 1000*kapak + ps*delta*dfdw;
  return (k11);
}

double saltmat3::k12(double x1,double /*x2*/,double /*x3*/,long ipp)
{
  double k12;
   double dcoef;
   double smc,a,n;
   int kod;
   CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

   dcoef=Tm->ip[ipp].eqother[1];
   k12 = -dcoef*x1;

  return (k12);
}

double saltmat3::k13(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k13;

  k13 = 0.0;

  return (k13);
}

double saltmat3::k21(double /*x1*/,double x2,double /*x3*/,long ipp)
{
  double k21;
  double kapak;

  kapak=Tm->ip[ipp].eqother[0];
  k21 = kapak*x2*1;
  
  return (k21);
}

double saltmat3::k22(double x1,double /*x2*/,double /*x3*/,long ipp)
{
  double k22;
  double dcoef;

   double smc,a,n;
   int kod;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

  dcoef=Tm->ip[ipp].eqother[1];
  k22 = dcoef*x1;

  return (k22);
}

double saltmat3::k23(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k23;

  k23 = 0.0;
  
  return (k23);
}

double saltmat3::k31(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k31;

  k31 = 0.0;

  return (k31);
}

double saltmat3::k32(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k32;

  k32= 0.0;

  return (k32);
}

double saltmat3::k33(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k33;

  k33= 0.0;

  return (k33);
}

double saltmat3::c11(double x1,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c11;
  double ps;
  double fi,dfdw, wsoli, por, w_hyg ,rh_hyg, hvezda;
  double u, a, n;
  int kod;
  double smc;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;


  ps = pgws(294.15);
  CorD(3,kd,1,x1,por,a2,a3);
  CorD(6,kod,1,x1,u,a,n);

  switch (kod)
  {
  case 2:
	  CorD(6,kod,1,x1,fi,a2,a3);
	  dfdw = derivation_dy_dx(6,x1,0,1);
	  break;
  case 30:
	  CorD(6,kod,0,x1,w_hyg,rh_hyg,a3);
            hvezda =grun ->sortpion_isotherm_root_shifted(x1,w_hyg ,rh_hyg);
            CorD(7,kod,0,x1,smc,a2,a3);
			grun->give_data_si_root_fi(x1,0.0,0.0,w_hyg, smc, rh_hyg,hvezda,fi);
			grun->give_data_si_root_dfidw(x1,0.0,0.0,w_hyg, smc, rh_hyg,hvezda,dfdw);
	  break;
  case 31:  {
	  dfdw = (u/(a*x1*n))*pow((1-(log(x1))/a),(-1-(1/n)));
	  fi = u*pow((1-(log(x1))/a),(-1/n));
	  //   printf ("\n w je %lf  fi je %lf a derivace je %lf",x1,fi, dfdw);
			}
	  break;
  }

 // if (ipp == 2) fprintf (Outt,"\n w je %lf  fi je %lf a derivace je %lf",x1,fi, dfdw);

  CorD(17,kd,0,x1,wsoli,a2,a3);

  c11 = 1000 + mw/(gasr*294.15)*ps*(por-x1-wsoli)*dfdw-mw*ps*fi/(gasr*294.15);

  return(c11);
}

double saltmat3::c12(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c12;

  c12 = 0.0;

  return(c12);
}

double saltmat3::c13(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c13;

  c13 = 0.0;

  return(c13);
}

double saltmat3::c21(double /*x1*/,double x2,double /*x3*/,long ipp)
{
  double c21,cfmax;

  cfmax=Tm->ip[ipp].eqother[3];

  if (x2 > cfmax)
  {
	  c21 = cfmax;
  }
  else
  {
    c21 = x2;
  }
  
  return(c21);
}

double saltmat3::c22(double x1,double /*x2*/,double /*x3*/,long ipp)
{
  double c22;
  double dcbdcf;
   double smc,a,n;
   int kod;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

 
  dcbdcf=Tm->ip[ipp].eqother[2];
  c22 = x1 + dcbdcf;
 
  return(c22);
}

double saltmat3::c23(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c23;

  c23 = 1.0;

  return(c23);
}

double saltmat3::c31(double /*x1*/,double x2,double /*x3*/,long ipp)
{
  double c31;
  double cfmax;

  cfmax=Tm->ip[ipp].eqother[3];
  
  if (x2 < cfmax){
    c31 = 0.0;
  }
  else{
    c31 = cfmax-x2;
  }
  
  return(c31);
}

double saltmat3::c32(double x1,double x2,double /*x3*/,long ipp)
{
  double c32;
  double cfmax;
   double smc,a,n;
   int kod;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;


  cfmax=Tm->ip[ipp].eqother[3];
  
  if (x2 < cfmax){
    c32 = 0.0;
  }
  else{
    c32 = -x1;
  }
  
  return(c32);
}

double saltmat3::c33(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c33;

  c33 = 1.0;

  return(c33);
}


/* Function calculates all auxiliary data - optional */
void saltmat3::auxiliarydata (double /*x1*/,double /*x2*/,double /*x3*/)
{}


/******************
Boundary conditions
*******************/
	
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
double saltmat3::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,x1,x2,x3;
  new_trc = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);
  x3 = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_11(x1,x2,x3,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  if((ri == 0) && (ci == 2))
    new_trc = 0.0;
  
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
    new_trc = 0.0;

  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double saltmat3::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval,x1,x2,x3;
  new_nodval = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);
  x3 = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_11(nodval,x1,x2,x3,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 0) && (ci == 2))
    new_nodval = 0.0;
  
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
    new_nodval = 0.0;

  return (new_nodval);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double saltmat3::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  double flux,x1,x2,x3;
  flux = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);
  x3 = nodalval (nn,2);
  
  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_11(nodval,x1,x2,x3,bc,ipp);
  if((ri == 0) && (ci == 1))
    flux = 0.0;
  if((ri == 0) && (ci == 2))
    flux = 0.0;
  
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
    flux = 0.0;

  return (flux);
}


/**
   function creates correct new nodal value on the boundary (transmission) for 1st medium
   @param new_nodval - new prescribed value near the boundary
   @param bv         - value of prescribed value near the boundary
   @param x1 ... x3  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double saltmat3::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    nodval = bv;//should be changed
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
   function creates correct transfer coefficient on the boundary (transmission) for 1st medium
   @param f11        - correct transfer coefficient
   @param x1 ... x3  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double saltmat3::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    trc=1.0;//should be changed
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
   function creates flux on the boundary (transmission - convective mass transfer) for 1st medium
   @param new_nodval - flux on the boundary
   @param bv         - prescribed value near the boundary
   @param x1 ... x3  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double saltmat3::get_transmission_flux_11(double bv,double x1,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//transmission - boundary flux
    flux = (bv - x1);//should be changed
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
   @param compother  - number of other components
   @param ipp        - first integration point on element
   @param x1 ... x3  - actual unknowns on the boundary
*/

double saltmat3::get_othervalue(long compother,long /*ipp*/, double x1,double /*x2*/,double /*x3*/)
{
  double other;

  switch (compother){
  case 0:{//first unknown
    other = x1;//should be changed
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
     @param out       - output file
     @param compother - number of other components
*/
void saltmat3::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//first unknown
    fprintf (out,"First unknown (Units)        ");//should be changed
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}
/** tlak nasycenych vodnich par
    @param t  - temperature K
*/

double  saltmat3::pgws(double t)
{
    double pgw;

	pgw = exp(23.5771 - 4042.9/(t - 37.58));

	return (pgw);
}
 /** Permeabilita Vodni Pary
    @param t  - temperature K
    @param p  -  presure Pa

*/
double saltmat3::permeabilitavodnipary(double t, double  p)
{

	double Da, Dp,rv,pa, mi;
	
	CorD(4,kd,0,0,mi,a2,a3);

	rv = 461.5;
	pa = 101325;

	Da = (2.306e-5 * pa)/(rv * t * p)*pow((t/294.15),1.81);
	Dp = Da/mi;
	return (Dp);
 }

void saltmat3::CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2)
{
int i;
  switch (MatChar[cislochar])
     {
     case 0:
     break;
     case 1: y = MatConst[cislochar];
     break;
     case 2:{
		 i = 1;
		 int psl, prad;
		 prad  = int(MatData[cislochar][0][0]);
		 psl = int(MatData[cislochar][1][0]);
		 if (x <=0)
		 {
			 x = 0;
			 y = MatData[cislochar][1][1];
			 if (psl == 3)
			 {
				 z = MatData[cislochar][2][1];
			 }
		 }
		 else 
		 {
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
			}
		 break;
     case 30:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
               kvyhl = 30;
			   break;
     case 31:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
               if (cislochar == 6)
			   {
				   z2 = MatFunce[cislochar][2];
			   }
               kvyhl = 31;
			   break;
     case 32:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
               kvyhl = 32;
			   break;
     case 33: //zatim nic
		 break;
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
		 kvyhl = 40;
			 }
		 break;
  }

}

void saltmat3::binding_izoterm_derivation(double x2, double & derbi)
{
  int i = 0;
  if ( x2 < 0)
  {
      x2 = 0;
  }
  
  while (x2 > MatData[15][0][i]){
	  i++;
  }
  
  derbi = (MatData [15][1][i] - MatData [15][1][i-1])/(MatData [15][0][i] - MatData [15][0][i-1]);
  //TODO: Add your source code here
}

void saltmat3::kapa_values (int kod, long ipp, double &kapa, double x1)
{
  double kapa1, kapa2, kapa3;
  double smc,a,n;
  CorD(7,kod,1,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

  switch (kod){
  case 0: kapa = 0.0;
	  break;
  case 1: CorD(5,kod,1,x1,kapa,kapa1,kapa2);
	  break;
  case 2: CorD(5,kod,1,x1,kapa,kapa1,kapa2);
	  break;
  case 30:{
	  CorD(5,kod,1,x1,kapa1,kapa2,kapa3);
	  kapa = kapa1*exp(kapa2*x1);
	  break;
		  }
  case 40: hystereze (5,0, kapa, ipp);
	  break;
  }
 /* if (ipp == 160)
  {
	  in1 = fopen("kapavyp.txt","a+");
	  fprintf(in1,"%e %e\n",x1,kapa);
	  fclose(in1);
  }*/

}

void saltmat3::salt_diffusivity_values (int kod, long ipp, double &diff, double x2)
{
  double diff1, diff2, diff3;
  
  
  switch (kod){
  case 0: diff = 0.0;
	  break;
  case 1: CorD(14,kod,1000,x2,diff,diff1,diff2);
	  break;
  case 2: CorD(14,kod,1000,x2,diff,diff1,diff2);
	  break;
  case 30:{
	  CorD(14,kod,1000,x2,diff1,diff2,diff3);
	  diff = diff1*exp(diff2*x2);
	  break;
		  }
  case 40: hystereze (14,0,diff, ipp);
	  break;
  }
  /*if (ipp == 160)
  {
	  in3 = fopen("Deckovyp.txt","a+");
	  fprintf(in3,"%e %e\n",x2,diff);
	  fclose(in3);
  }*/

}

void saltmat3::binding_izoterm_derivation_values (int kod, long ipp, double &dcbdcf, double x2)
{
  double a,b,a3,Cb,x1;
  
  x1 = Tm->ip[ipp].av[0];
  
  if (x2 <= 0) x2 = 1e-10;
  if (x1 <0) x1 = 0.0;

  switch (kod)
     {
  case 0: dcbdcf = 0.0;
	  break;
  case 2: binding_izoterm_derivation(x2,dcbdcf);
	  break;
  case 30:{
	  CorD(15,kod,1000,x2,a,b,a3);
	  Cb =  1/(1/(a*b*x2)+1/b);
	  dcbdcf = Cb*Cb/(a*b*x2*x2);
	  break;
		  }
  case 31:{
	  CorD(15,kod,1000,x2,a,b,a3);
	  dcbdcf = b*a*pow(x2,(b-1));
	  break;
		  }
  case 40: hystereze (15,0,dcbdcf, ipp);
	  break;
  }
/*  if (ipp == 160)
  {
	  in2 = fopen("VIder.txt","a+");
	  fprintf(in2,"%e %e\n",x2,dcbdcf);
	  fclose(in2);
  }
*/
}



/**
   function computes auxiliary values which are necessary for future computation
   
   JM, 29.5.2007
*/
void saltmat3::aux_values (long ipp)
{
  double x1,x2,x3,a1,a2,cfmax,jm1, d, dcbdcf;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];

  CorD(16,kd,0,x1,cfmax,a1,a2);
  Tm->ip[ipp].eqother[3]=cfmax;
  
  kapa_values(MatChar[5],ipp,jm1, x1); 
  Tm->ip[ipp].eqother[0]=jm1;
  
  salt_diffusivity_values(MatChar[14],ipp,d, x2);
  Tm->ip[ipp].eqother[1]=d;

  binding_izoterm_derivation_values(MatChar[15],ipp,dcbdcf, x2);
  Tm->ip[ipp].eqother[2]=dcbdcf;
}


/**
   function defines inital values of quantities at integration points
   
   @param ipp - integration point pointer
   @param ido - index in array eq_other   

   JM, 29.5.2007
*/
void saltmat3::initvalues (long ipp,long /*ido*/)
{
  if (MatChar[5] == 40)  Tm->ip[ipp].eqother[0]=Init[5];
  if (MatChar[14] == 40) Tm->ip[ipp].eqother[1]=Init[14];
  if (MatChar[15] == 40) Tm->ip[ipp].eqother[2]=Init[15];
  Tm->ip[ipp].eqother[4] = 3800;
  
}

void saltmat3::hystereze (int matchar, int /*matchar2*/, double & outvalue, long ipp)
{
  double x1,x2,x3,xpv,xav,yav,rhpv, newder,ypv,yavA,yavD, ypvA, ypvD, derpvA,derpvD, other1;
  int hyst, k, other;
  //double  k1,k2, k3, k4, k5, k6;
  double t;
  t = Tp->time;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
  
 
  switch(matchar){
  case 5: k = 0;
	  xav = x1;
	  xpv = Tm->ip[ipp].pv[0];
	  rhpv = Tm->ip[ipp].pv[0];
//	  k1 = k1[matchar];
	  /*k2[5] = 10;
	  k3[5] = 1000;
	  k4[5] = 5;
	  k6[5] = 100;
	  k1[5] = 0.1;
	  k4 =0.1*/;
	  CorD(5,other,1,xav,yavA,yavD,other1);
	  CorD(5,other,1,xpv,ypvA,ypvD,other1);
	  ypv = Tm->ip[ipp].eqother[0];
	  break;
  case 14: k = 1;
	  xav = x2;
	  xpv = Tm->ip[ipp].pv[1];
	  rhpv = Tm->ip[ipp].pv[0];
	  CorD(14,other,1000,xav,yavA,yavD,other1);
	  CorD(14,other,1000,xpv,ypvA,ypvD,other1);
	  /*k2 = 0.1;
	  k1 = 10;
	  k3 = 100;
	  k4 = 2000;
	  k5 = 0.1;
	  k6 = 800;*/
	  ypv = Tm->ip[ipp].eqother[1];
	  break;
  case 15: k=2;
	  xav = x2;
	  if(t == 3800)
	  {
		  xpv = 100.0;
	  }
	  else
	  {
		  xpv = Tm->ip[ipp].pv[1];
	  }

	  rhpv = Tm->ip[ipp].pv[0];
	  CorD(15,other,1000,xav,yavA,yavD,other1);
	  CorD(15,other,1000,xpv,ypvA,ypvD,other1);
	  /*k2 = 1;
	  k1 = 0.5;
	  k3 = 300;
	  k4 = 1;
	  k5 = 0.5;
	  k6 = 300;*/
	  if( t== 3800)
	  {
		  ypv =Tm->ip[ipp].eqother[2]= 110.0;
	  }
	  else
	  {
		  ypv = Tm->ip[ipp].eqother[5];
	  }
	  break;
  }

  
  
  der_value_hyst(matchar,40,xpv,derpvA,derpvD,ipp); 

  if (x1 > rhpv){
	  hyst = 0;
  }
  else{
	  hyst = 1;
  }

  if(yavA>yavD)
  {
	  if (x1 > rhpv){
		  hyst = 2;
	  }
	  else{
		  hyst = 3;
	  }
  }
  
  switch(hyst){
  case 0:{			// navlhani
	  newder = (k1[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k2[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k3[matchar]  + ypv;
	  
	  if (yav < yavA)
		  {
			  yav = yavA;
		  }
	  if (yav > yavD) yav = yavD;
	  break;
		 }
  case 1:{			// vysychani
	  newder  = (k4[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k5[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k6[matchar] + ypv;
	  if (yav > yavD)
		  {
			  yav = yavD;
		  }
	  if(yav < yavA) yav = yavA;
	  break;
		 }
  case 2:{			// navlhani
	  newder = (k1[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k2[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k3[matchar]+ ypv;
	  if (yav > yavA)
		  {
			  yav = yavA;
		  }
	  if(yav < yavD) yav = yavD;
	  break;
		 }
  case 3:{			// vysychani
	  newder  = (k4[matchar]*derpvD*(ypv-ypvA)*(ypv-ypvA)+k5[matchar]*derpvA*(ypv -ypvD)*(ypv - ypvD))/((ypvD-ypvA)*(ypvD-ypvA));
	  yav = (xav-xpv)*newder*k6[matchar] + ypv;
	  if (yav < yavD)
		  {
			  yav = yavD;
		  }
	  if(yav > yavA)yav = yavA;
	  break;
		 }
  }

  if (matchar == 15) 
    {
      if (yav == yavA)
	{
	  newder = (yav-ypv)/(xav - xpv);
	  if ((xav - xpv)== 0)
	    newder = 0;
	  
	  Tm->ip[ipp].eqother[k]=newder;
	}
      else
	{
	  if (yav == yavD)
	    {
	      newder = (yav-ypv)/(xav - xpv);
	      if ((xav - xpv)== 0)
		newder = 0;
	      Tm->ip[ipp].eqother[k]=newder;
	    }
	  else
	    {
	      Tm->ip[ipp].eqother[k]=newder;
	      
	    }
	}
      
      Tm->ip[ipp].eqother[5]=yav;
    }
  else
  {
	  Tm->ip[ipp].eqother[k]=yav;
  }
  
  //if (t== Tm->ip[ipp].eqother[4])
  //{
	  if (ipp ==160)
	  {
		  switch (matchar){
		  case 5: {
			  in1 = fopen("kapa.txt","a+");
			  fprintf(in1,"%lf %e %e %e %e\n",t,xav, yavA, yavD,yav);
			  fclose(in1);
				  }
			  break;
		  case 14: {
			  in2 = fopen("decko.txt","a+");
			  fprintf(in2,"%lf %e %e %e %e\n",t,xav, yavA, yavD,yav);
			  fclose(in2);
				  }
			  break;
		  case 15:{
			  in3 = fopen("vi.txt","a+");
			  fprintf(in3,"%lf %e %e %e %e\n",t,xav, yavA, yavD,yav);
			  fclose(in3);
				  }
			  break;
		  }
	//  }
	 // Tm->ip[ipp].eqother[4] = t+ 60;

	  }
	  if(matchar == 15)
	  {
		  outvalue = newder;
	  }
	  else
	  {
		  outvalue = yav;
	  }
}

void saltmat3::der_value_hyst (int matchar,int kod, double pv, double & outvalue,double & outvalue2, long /*ipp*/)
{
	int i;


	switch (kod){
	case 0:
		break;
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
		i = 1;
		while (pv > MatData[matchar][0][i]){
			i++;
		}
		outvalue = (MatData [matchar][1][i] - MatData [matchar][1][i-1])/(MatData [matchar][0][i] - MatData [matchar][0][i-1]);
		i = 1;
		while (pv > MatData[matchar][2][i]){
			i++;
		}
		outvalue2 = (MatData [matchar][3][i] - MatData [matchar][3][i-1])/(MatData [matchar][2][i] - MatData [matchar][2][i-1]);
			}
		break;
	}
}

double saltmat3::get_moisture(double rh)
{
  int s;
 double moist,hmrh, mhmc, smc, u, a, n;

       if (rh < 0) rh = 0;
	  
	   //CorD(6,s,1,rh,u,a,n);
	   s = MatChar[6];

       //CorD(6,s,0,0,0.0,a1,a2);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
          case 2:   CorD(6,s,1,rh,moist,a1,a3);
          break;
          case 30:  CorD(6,s,1,rh,a1,hmrh,a3);
                    CorD(7,s,0,rh,smc,a2,a3);
                    CorD(6,s,0,rh,mhmc,hmrh,a3);
                    //kuns->give_data(rh,mhmc, smc, hmrh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,1,rh,u,a,n);
	    //moist = kuns->si_kk_hansen(rh, 0,u,a,n);
          break ;
          }
 return(moist);
    //TODO: Add your source code here
}

double saltmat3::get_rel_hum(double w)
{
  int kod;
  double relhum;
  double	u, a, n, a2, a3,w_hyg,rh_hyg, hvezda, smc;
  
  //   sorption_izoterms_giva_data(0,w,0,0,relhum,a1);    //TODO: Add your source code here
  
  //  this variable has to be defined
  relhum=0.0;

  CorD(6,kod,1,relhum,u,a,n);
  
  switch (kod)
    {
    case 2:   CorD(6,kod,1,w,relhum,a2,a3);
      break;
    case 30:  //CorD(6,s,1,x1,a1,hmrh,a3);
      CorD(6,kod,0,w,w_hyg,rh_hyg,a3);
      hvezda =grun ->sortpion_isotherm_root_shifted(w,w_hyg ,rh_hyg);
      CorD(7,kod,0,w,smc,a2,a3);
      grun->give_data_si_root_fi(w,0.0,0.0,w_hyg, smc, rh_hyg,hvezda,relhum);
      
      break;
    case 31: 
      relhum  = a1*pow((1-(log(w))/a2),(-1/a3));
      break ; 
    }
  
  
  return (relhum);
}

double saltmat3::derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2)
{
  int i = 1;
  double outvalue;
  
  while (prom > MatData[matchar][pomk1][i]){
    i++;
  }
  outvalue = (MatData [matchar][pomk2][i] - MatData [matchar][pomk2][i-1])/(MatData [matchar][pomk1][i] - MatData [matchar][pomk1][i-1]);
  return(outvalue);
  
}
