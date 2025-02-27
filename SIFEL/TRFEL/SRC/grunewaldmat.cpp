#include "globalt.h"
#include "grunewaldmat.h"
#include "kunzel.h"
#include "globmatt.h"

static double Tref = 293.15; // reference temperature
static double Rvap = 462.0;  // specific gas constant of water vapour
static double cwat = 4180;   // specific heat capacity of pure liquid water
static double cvap = 1617;   // specific heat capacity of water vapour
static double hvap = 2256000; // specific enthalpy of water vapour at 0?C
static double alfal = 0.00; // thermal expansion coefficient of the liquid phase
static double rowref = 1000;


/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
***********************************************/

grunewaldmat::grunewaldmat ()
{
}

grunewaldmat::~grunewaldmat ()
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
void grunewaldmat::matcond (matrix &d,long ri,long ci,long ipp)
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
void grunewaldmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = -1000000.0;
  if(x1<0)
{
	x3 = -1000000.0;
}
 // x1 = Tm->ip[ipp].pv[0];
    
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

  
  //fprintf (Outt2,"Time je %lf ipp %ld  k %20.15lf \n",Tp->time,ipp,k);
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void grunewaldmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = -1000000.0;


 // x1 = Tm->ip[ipp].pv[0];
    
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

void grunewaldmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = -1000000.0;

    
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
void grunewaldmat::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2,x3;
  c = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = -1000000.0;

    
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

  //fprintf (Outt2,"Time je %lf ipp %ld  c %20.15lf \n",Tp->time,ipp,c);
  
}



/**
   function reads data and material parameters

   @param in  - input file
*/
void grunewaldmat::read (XFILE *in)
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
                              break;
                              case 2:   if( MatChar[i-1] ==2)
                                             {
                                              xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
                                              }
                                        else
                                             {
                                             xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
                                             }

                                        switch (int(MatData[i][1][0]))
                                             {
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
                                     //   edit_chechbox(i,2);
                              break;
                              case 3://edit_chechbox(i,3);
                              break;
                              case 30:   if( MatChar[i-1] ==1)
                                             {
                                               xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
                                              }
                                        else
                                             {
                                             xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
                                             }

                              break;
                              case 31: if( MatChar[i-1] ==1)
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
                              case 32:  if( MatChar[i-1] ==1)
                                             {
                                                xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
                                              }
                                        else
                                             {
                                               xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
                                             }
                              break;
                              }
                         }

  madripom = 0;

}

void grunewaldmat::print(FILE *out)
{
  fprintf (out,"\n grunewald material \n");
}



/*position CORD:     2 - density
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
                    18 - 19 - none    */




/******************************************************
Conductivity and capacity terms for C and K matrices
*******************************************************/




double grunewaldmat::k11(double x1,double x2,double x3,long ipp)
{
  double k11;
  rol = density_lql(x1,x2,x3);
 // kapa = give_kapa(x1,x2);// kapa
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  dvn = diffusion_number_function(x1,x2,x3)/mi;
  dp = dvn/(Rvap*x2);

  //sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfidw); // jeste doladit
  
  dfidw =  Tm->ip[ipp].eqother[1];
  kapa =  Tm->ip[ipp].eqother[3];


  k11 = rol*kapa + pvs *dp *dfidw;


  return (k11);
}

double grunewaldmat::k12(double x1,double x2,double x3,long ipp)
{
  double k12;
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  dpvsdt =  derivation_saturation_water_vapour_pressure_temperature(x1,x2,x3);
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  dvn = diffusion_number_function(x1,x2,x3)/mi;
  dp = dvn/(Rvap*x2);
 // sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1); // jeste doladit

  dfidt = 0.0; // neni namereno
  
  relhum =  Tm->ip[ipp].eqother[0];
  

  k12 = dp*dpvsdt*relhum + pvs*dp*dfidt;

  return (k12);
}

double grunewaldmat::k13(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k13;

  //here will be something??!!
  k13 = 0.0;
  return (k13);
}

double grunewaldmat::k21(double x1,double x2,double x3,long ipp)
{
  double k21;

  rol = density_lql(x1,x2,x3);
   //kapa = give_kapa(x1,x2);// kapa
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  dvn = diffusion_number_function(x1,x2,x3)/mi;
  dp = dvn/(Rvap*x2);
  ul = specific_internal_energy_of_the_liquid_phase(x1,x2,x3);
  hv = specific_enthalpy_of_water_vapour_hv(x1,x2,x3);
  //sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfidw); // jeste doladit


  
  dfidw =  Tm->ip[ipp].eqother[1];
  kapa =  Tm->ip[ipp].eqother[3];

  k21 = rol*kapa*ul + hv* pvs* dp*dfidw;
  return (k21);
}

double grunewaldmat::k22(double x1,double x2,double x3,long ipp)
{
  double k22;

  CorD(10,kd,1.0,0,x1,lambda,a2,a3);//lambda
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  dvn = diffusion_number_function(x1,x2,x3)/mi;
  dp = dvn/(Rvap*x2);
  //sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1); // jeste doladit

  hv = specific_enthalpy_of_water_vapour_hv(x1,x2,x3);
  dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2,x3);
  dfidt = 0.0;// neni namereno



  relhum =  Tm->ip[ipp].eqother[0];
  
  k22 = lambda+dp*relhum*hv*dpvsdt + dp*pvs*hv*dfidt;

  return (k22);
}

double grunewaldmat::k23(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k23;

  k23 = 0.0;
  //here will be something??!!

  return (k23);
}

double grunewaldmat::k31(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k31;

  //here will be something??!!

  k31 = 0.0;
  return (k31);
}

double grunewaldmat::k32(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k32;

  //here will be something??!!
  k32 = 0.0;

  return (k32);
}

double grunewaldmat::k33(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double k33;

  //here will be something??!!
  k33 = 0.0;

  return (k33);
}

double grunewaldmat::c11(double x1,double x2,double x3,long ipp)
{
  double c11;

  row = density_lqw(x1,x2,x3);
  rov = ((partial_water_vapour_pressure_function(x1,x2,x3))/(Rvap*x2));
  CorD(3,kd,0,0,x1,pir,a2,a3);//pir porovitost
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
//  sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfidw); // jeste doladit
  //here will be something??!!


  dfidw =  Tm->ip[ipp].eqother[1];


  c11 = row - rov + (pir - x1)/(Rvap * x2)*pvs*dfidw;

  return(c11);
}

double grunewaldmat::c12(double x1,double x2,double x3,long ipp)
{
  double c12;

  drowdt = derivation_density_drowdt(x1,x2,x3);
  CorD(3,kd,0,0,x1,pir,a2,a3);//pir porovitost
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  dfidt = 0.0; // neni namereno
  //sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1); // jeste doladit

  dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2,x3);
  pv = partial_water_vapour_pressure_function(x1,x2,x3);
 
  relhum =  Tm->ip[ipp].eqother[0];

  //here will be something??!!
  c12 = (x1*drowdt + (pir - x1)/(Rvap*x2)*(pvs *dfidt + relhum* dpvsdt-pv/x2));

  return(c12);
}

double grunewaldmat::c13(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c13;

  //here will be something??!!
  c13  = 0.0;
  return(c13);
}

double grunewaldmat::c21(double x1,double x2,double x3,long ipp)
{
  double c21;

  rol = density_lql(x1,x2,x3);
  ul = specific_internal_energy_of_the_liquid_phase(x1,x2,x3);
  uv = specific_internal_energy_of_water_vapour(x1,x2,x3);
  CorD(3,kd,0,0,x1,pir,a2,a3);//pir porovitost
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
//  sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfidw); // jeste doladit
  rov = ((partial_water_vapour_pressure_function(x1,x2,x3))/(Rvap*x2));
   pv = partial_water_vapour_pressure_function(x1,x2,x3);


  dfidw =  Tm->ip[ipp].eqother[1];

  c21 = rol*ul + uv*(pir - x1)/(Rvap*x2)*pvs * dfidw - uv*rov;
  return(c21);
}

double grunewaldmat::c22(double x1,double x2,double x3,long ipp)
{
  double c22;

  CorD(2,kd,0,0,x1,rom,a2,a3);// ro matrice
  dumdt = derivation_of_specific_internal_energy_of_the_solid_material_dependence_temperature(x1,x2,x3);
  CorD(3,kd,0,0,x1,pir,a2,a3);//pir porovitost
  rol = density_lql(x1,x2,x3);
  duldt = derivation_of_specific_internal_energy_of_the_liquid_phase_dependence_temperature(x1,x2,x3);
  rov = ((partial_water_vapour_pressure_function(x1,x2,x3))/(Rvap*x2));
  duvdt = derivation_specific_internal_energy_of_water_vapour(x1,x2,x3) ; // dodelat
  ul = specific_internal_energy_of_the_liquid_phase(x1,x2,x3);
  droldt = derivation_density_droldt(x1,x2,x3);
  uv = specific_internal_energy_of_water_vapour(x1,x2,x3);
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  dfidt = 0.0; // neni namereno
  //sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1); // jeste doladit
  /*double pokus;
  pokus = get_moisture(relhum); */
  dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2,x3);
  pv = partial_water_vapour_pressure_function(x1,x2,x3);
  
  relhum =  Tm->ip[ipp].eqother[0];

  c22 = rom *dumdt + rol * x1*duldt + rov *duvdt*(pir - x1)+ul * x1*droldt + uv *(pir - x1)/(Rvap *x2)*(pvs *dfidt+relhum*dpvsdt-pv/x2);


  return(c22);
}

double grunewaldmat::c23(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c23;

  //here will be something??!!
  c23 = 0.0;
  return(c23);
}

double grunewaldmat::c31(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c31;

  //here will be something??!!
  c31  = 0.0;
  return(c31);
}

double grunewaldmat::c32(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c32;

  //here will be something??!!
  c32  = 0.0;
  return(c32);
}

double grunewaldmat::c33(double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
{
  double c33;

  //here will be something??!!
  c33 = 0.0;
  return(c33);
}


/* Function calculates all auxiliary data - optional */
void grunewaldmat::auxiliarydata (double /*x1*/,double /*x2*/,double /*x3*/,long /*ipp*/)
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
double grunewaldmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
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
  
  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
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
double grunewaldmat::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
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
  
  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
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
double grunewaldmat::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
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
  
  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
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
double grunewaldmat::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double nodval;

  if (bc>10){
    nodval = bv;//should be changed
  } 
  else{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
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
double grunewaldmat::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double trc;
  
  if (bc>10){
    trc=1.0;//should be changed
  }
  else{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
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
double grunewaldmat::get_transmission_flux_11(double bv,double x1,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
{
  double flux;

  if (bc>10){
    flux = (bv - x1);//should be changed
  }
  else{
    print_err("no real boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  
  return(flux);
}

/**
   function computes all variables in nodes
   @param compother  - number of other components
   @param ipp        - first integration point on element
   @param x1 ... x3  - actual unknowns on the boundary
*/

double grunewaldmat::get_othervalue(long compother,long /*ipp*/, double x1,double x2,double x3)
{
  double other;

  switch (compother){
  case 0:{//moisture
    other = x1;//should be changed
      break;
  }
  case 1:{//temperature
    other = x2;//should be changed
      break;
  }
  case 2:{//relative humidity
    //other = x1;//should be changed
    double relhum, a1;
    sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1);
    other = relhum;
      break;
  }
  case 3:{//relative humidity
    //other = x1;//should be changed
    other = 0.0;
      break;
  }
  case 5:
  case 6:
  case 7:
  case 4:{//relative humidity
    other = 0.0;
      break;
  }
  
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);
}

/**
     function prints names of all variables in nodes
     @param out       - output file
     @param compother - number of other components
*/
void grunewaldmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//first unknown
    fprintf (out,"Moisture (m3/m3)        ");//should be changed
    break;
  }
  case 1:{//temperature
    fprintf (out,"Temperature (K)        ");//should be changed
    break;
  }
  case 2:{//relative humidity
    fprintf (out,"Relative humidity (--)        ");//should be changed
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}



double grunewaldmat::density_lql(double /*x1*/,double x2, double /*x3*/)
{
 /* density of the liquid phase rol
 Tr = T - Tref
 without salt transport
 rowref ...... mass density of pure liquid water for T = Tref
 alfa1  ...... thermal expansion coefficient of the liquid phase
 Tref   ...... reference temperature
 */
 double tr;
 tr = x2 - Tref;

 rol = rowref * exp(alfal * tr);
 return(rol);
 }

void grunewaldmat::give_data_si_root_fi(double x1, double /*x2*/, double /*x3*/, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh)
{

  if (x1 < w_hyg)
     {
     relh = 1 - (1-x1/shift_w)*(1-x1/shift_w);
     }
  else
     {
     if (x1 > w_sat)
          {
          relh = 1;
          }
     else
          {
          relh = rh_hyg+(1-rh_hyg)*(x1-w_hyg)/(w_sat-w_hyg);
          }
     }

 // moissat = Smc ;
    //TODO: Add your source code here
}
//---------------------------------------------------------------------------

double grunewaldmat::kapa_exp(double a, double b, double x1, double /*x2*/, double /*x3*/)
{
 if (x1 <= 0) x1 = 0;

// sorption_izoterms_giva_data(0,x1,x2,x3,moist,a3);

 return (a * exp(b * x1));

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

void grunewaldmat::CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2)
{

  switch (MatChar[cislochar])
     {
     case 0:
     break;
     case 1: y = MatConst[cislochar];
		 kvyhl = 1;
		 z= -10000.0;
		 z2= -10000.0;
     break;
     case 2:   {
                int i = 1;

                int psl, prad;
                prad  = int (MatData[cislochar][0][0]);
                psl = int (MatData[cislochar][1][0]);
                if(rhw == 1)
                        {
                        if (x <=0)
                                {
                                x = 0;
                                y = MatData[cislochar][0][1];
                                if (psl == 3)
                                        {
                                        z = MatData[cislochar][2][1];
                                        }

                                }
                        else {
                                if (x > in)
                                        {
                                        x = in;
                                        y = MatData[cislochar][0][prad];
                                        if (psl == 3)
                                                {
                                                z = MatData[cislochar][2][prad];
                                                }
                                        }
                                else
                                        {
                                        while (x > MatData[cislochar][1][i]){
                                        i++;
                                        }

                                        y = MatData[cislochar][0][i-1] + ((x - MatData[cislochar][1][i-1])*(MatData[cislochar][0][i] - MatData[cislochar][0][i-1]))/(MatData[cislochar][1][i] - MatData[cislochar][1][i-1]);
                                        if (psl == 3)
                                                {
                                                z = MatData[cislochar][2][i-1] + ((x - MatData[cislochar][1][i-1])*(MatData[cislochar][2][i] - MatData[cislochar][2][i-1]))/(MatData[cislochar][1][i] - MatData[cislochar][1][i-1]);
                                                }
                                        }
                                }
                        }
                else    {
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
                }
                kvyhl = 2;
				z2=MatData[cislochar][1][prad];
                }
     break;
     case 30:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
			   z2= -10000.0;
               kvyhl = 30;
     break;
     case 31:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
               if (cislochar == 6)
                    {
                    z2 = MatFunce[cislochar][2];
                    }
			   else {
				    z2= -10000.0;
			   }
               kvyhl = 31;
     break;
     case 32:  y = MatFunce[cislochar][0];
               z = MatFunce[cislochar][1];
			   z2= -10000.0;
               kvyhl = 32;
     break;
     case 33: //zatim nic
     break;
     }
}

void grunewaldmat::sorption_izoterms_giva_data(int /*kod*/,double x1, double x2, double x3, double & fi, double & dfdw)
{
 /*
 kod ..... 0- jen vlhkost, 1 - jen derivaci...
 */
 int s;
 double w_hyg,rh_hyg, smc, hvezda;
       
       if (x1 < 0) x1 = 0;

       CorD(6,s,0,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
	   case 2:{
		   CorD(7,s,0,0,x1,smc,a2,a3);
		   CorD(6,s,smc,1,x1,fi,a2,a3);
                if(x1>a3)
				{
					CorD(7,s,0,0,x1,smc,a2,a3);
					CorD(6,s,smc,0,x1,a2,rh_hyg,w_hyg);
					rh_hyg = 0.976;
					hvezda =sortpion_isotherm_root_shifted(x1,w_hyg ,0.976);
					
					give_data_si_root_fi(x1,x2,x3,w_hyg, smc, rh_hyg,hvezda,fi);
					give_data_si_root_dfidw(x1,x2,x3,rh_hyg,smc,w_hyg,hvezda,dfdw);
				}
				else{
					sorption_izotherm_derivation(x1,x2,x3,dfdw);
				}
			  }
          break;
	   case 30:{  //CorD(6,s,1,x1,a1,hmrh,a3);
                      CorD(6,s,0,0,x1,w_hyg,rh_hyg,a3);
                      hvezda =sortpion_isotherm_root_shifted(x1,w_hyg ,rh_hyg);
                      CorD(7,s,0,0,x1,smc,a2,a3);
                      give_data_si_root_fi(x1,x2,x3,w_hyg, smc, rh_hyg,hvezda,fi);
					  give_data_si_root_dfidw(x1,x2,x3,rh_hyg,smc,w_hyg,hvezda,dfdw);
			   }
          break;
         /* case 31:  switch (kod)
                         {
                         case 0:   CorD(6,kod,0,0,u,a,n);
                                   moist = si_kk_hansen(x1,x2,x3,u,a,n);
                         break;
                         case 1:   //dmoistdrh = soptionizothermDerivation(x1, x2);
                         break;
                         }
          break ; */
        }

   // printf ("\n rh je %lf derivace je %lf",rh, DrvDrh);
     
}

double grunewaldmat::si_kk_hansen(double x1, double /*x2*/, double /*x3*/,double u, double a, double n)
{
    return (u*pow((1-(log(x1))/a),(-1/n)));
}

double grunewaldmat::derivation_density_droldt(double x1, double x2, double x3)
{
 return  (alfal *  density_lql(x1,x2,x3));  //TODO: Add your source code here
}

double grunewaldmat::density_lqw(double x1, double x2, double x3)
{
    return(density_lql(x1,x2,x3)); //TODO: Add your source code here
}
 
double grunewaldmat::derivation_density_drowdt(double x1, double x2, double x3)
{
   return (derivation_density_droldt(x1,x2,x3));  //TODO: Add your source code here
}

double grunewaldmat::saturation_water_vapour_pressure(double /*x1*/, double x2, double /*x3*/)
{
  return (exp(23.5771 - 4042.9/(x2 - 37.58)));   //TODO: Add your source code here
}

double grunewaldmat::derivation_saturation_water_vapour_pressure_temperature(double /*x1*/, double x2, double /*x3*/)
{
  return ((4042.9 * exp(23.5771-4042.9/(x2 - 37.58)))/((x2-37.58)*(x2-37.58)));   //TODO: Add your source code here
}

double grunewaldmat::specific_internal_energy_of_the_liquid_phase(double x1, double x2, double x3)
{
 return (specific_internal_energy_of_the_liquid_water(x1,x2,x3));       //TODO: Add your source code here
}

double grunewaldmat::specific_internal_energy_of_the_liquid_water(double /*x1*/, double x2, double /*x3*/)
{
  return(cwat*(x2 - 273.15));      //TODO: Add your source code here
}

double grunewaldmat::specific_internal_energy_of_water_vapour(double /*x1*/, double x2, double /*x3*/)
{
 return(cvap*(x2 - 273.15));       //TODO: Add your source code here
}

double grunewaldmat::derivation_of_specific_internal_energy_of_the_solid_material_dependence_temperature(double x1, double /*x2*/, double /*x3*/)
{
   double cmat;
   CorD(9,kd,0,0,x1,cmat,a2,a3);;
   return(cmat);      //TODO: Add your source code here
}

double grunewaldmat::derivation_of_specific_internal_energy_of_the_liquid_phase_dependence_temperature(double /*x1*/, double /*x2*/, double /*x3*/)
{
  // specific internal energy of the liquid water-- derivation equal .... liquid phase
  return(cwat);      //TODO: Add your source code here
}

double grunewaldmat::partial_water_vapour_pressure_function(double x1, double x2, double x3)
{
  double relhum;
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
    sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1);
  return(pvs*relhum);
        //TODO: Add your source code here
}

double grunewaldmat::diffusion_number_function(double /*x1*/, double x2, double /*x3*/)
{
 return(2.3e-5 * pow(x2/273,1.81));       //TODO: Add your source code here
}

double grunewaldmat::specific_enthalpy_of_water_vapour_hv(double /*x1*/, double x2, double /*x3*/)
{
   return(cvap*(x2-273.15)+hvap);  //TODO: Add your source code here
}

void grunewaldmat::sorption_izotherm_derivation(double x1, double /*x2*/, double /*x3*/, double & derfi)
{
 int i = 1;
if ( x1 < 0)
     {
     x1 = 0;
      }

 while (x1 > MatData[6][1][i]){
                              i++;
                              }
 derfi = (MatData [6][0][i] - MatData [6][0][i-1])/(MatData [6][1][i] - MatData [6][1][i-1]);
    //TODO: Add your source code here       //TODO: Add your source code here
}



double grunewaldmat::derivation_specific_internal_energy_of_water_vapour(double /*x1*/, double /*x2*/, double /*x3*/)
{
 return(cvap );    //TODO: Add your source code here
}

void grunewaldmat::get_rel_hum(double w, double &fi, double &dfdw)
{
   sorption_izoterms_giva_data(0,w,0,0,fi,dfdw);    //TODO: Add your source code here
   
}

double grunewaldmat::sortpion_isotherm_root_shifted(double /*x1*/, double w_hyg, double rh_hyg)
{
  return(w_hyg/(1-sqrt(1-rh_hyg)));
}

void grunewaldmat::give_data_si_root_dfidw(double x1, double /*x2*/, double /*x3*/, double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw)
{

     if ( x1< 0) x1 = 0;
     if (x1>w_sat) x1 = w_sat;

     if(x1<w_hyg)
          {
          dfdw = 2*(shift_w-x1)/(shift_w*shift_w);
          }
     else
          {
          dfdw = (1-rh_hyg)/(w_sat-w_hyg);
          }
}

double grunewaldmat::get_moisture(double rh)
{
  int s;
 double moist,hmrh, mhmc, smc, u, a, n;

       if (rh < 0) rh = 0;

       CorD(6,s,0,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
	   case 2:   {
			if(rh > 0.976)
			{
				CorD(7,s,0,0,rh,smc,a2,a3);
				CorD(6,s,0,0,rh,a2,a3,mhmc);
				hmrh = 0.976;
				//kun->give_data(rh,mhmc, smc, hmrh,moist);
			}
			else
			{
				CorD(6,s,1,0,rh,moist,a1,a3);
			}
				 }

          break; 
          case 30:  CorD(6,s,1,0,rh,a1,hmrh,a3);
                    CorD(7,s,0,0,rh,smc,a2,a3);
                    CorD(6,s,0,0,rh,mhmc,hmrh,a3);
                    //kun->give_data(rh,mhmc, smc, hmrh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,0,0,0,u,a,n);
	    //moist = kun->si_kk_hansen(rh, 0,u,a,n);
          break ;
          }
 return(moist);
    //TODO: Add your source code here
}

double grunewaldmat::give_kapa(double x1, double x2)
{
	int s;
	double kapa,u, a2, a3, smc,x3;
	x3 = 0.0;
	CorD(7,s,1,0,x1,smc,a2,a3);

	if (x1 < 0) x1 = 0;
	if (x1> smc) x1 = smc;

	CorD(5,s,1,0,x1,u,a2,a3);// kapa
	
	switch (s){
	case 1:
	case 2: CorD(5,s,1,0,x1,kapa,a2,a3);// kapa
		break;
	case 30: kapa = kapa_exp(u,a2,x1,x2,x3);
		break;
	}
 return(kapa);
    //TODO: Add your source code here
}

/*


inicializace.... 

*/
void grunewaldmat::initvalues (long ipp,long /*ido*/)
{
  double x1,x2;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
}

void grunewaldmat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
  Tm->ip[ipp].eqother[5]=out[5];

}

void grunewaldmat::give_values (long ipp,double *av, double */*pv*/, double *eq)
{
	av[0] = Tm->ip[ipp].av[0];
	av[1] = Tm->ip[ipp].av[1];

//	pv[0] = Tm->ip[ipp].pv[0];
//	pv[1] = Tm->ip[ipp].pv[1];
  

	eq[0] = Tm->ip[ipp].eqother[0];
	eq[1] = Tm->ip[ipp].eqother[1];
	eq[2] = Tm->ip[ipp].eqother[2];
	eq[3] = Tm->ip[ipp].eqother[3];
	eq[4] = Tm->ip[ipp].eqother[4];
	eq[5] = Tm->ip[ipp].eqother[5];

    
}

void grunewaldmat::aux_values (long /*ipp*/,double *in,double *inp, double *ineq,double *out)
{
  double x1,x2,x3,x1p, x2p, ypv1,ypv2,ypv3,ypv4,ypv5, ypv0;
  
  double rh,dfidw, kapa;
  double t;
  t = Tp->time;
  x1 = in[0];
  x2 = in[1];
  x3 = -100000;


  x1p = inp[0];
  x2p = inp[1];
  

  ypv0 =ineq[0];
  ypv1 =ineq[1];
  ypv2 =ineq[2];
  ypv3 =ineq[3];
  ypv4 =ineq[4];
  ypv5 =ineq[5];



 sorption_izoterms_giva_data(0,x1,x2,x3,rh,dfidw); // jeste doladit
 // sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1); // jeste doladit



  //sorption_izoterms_giva_data(0,x1,x2,w,dwdf);
 // rh = get_rel_hum(x1);
  kapa = give_kapa(x1,x2);
  // sorption_izoterms_values(MatChar[6],ipp,x1,x1p,ypv1,w,dwdf); 
  //kapa_values(MatChar[5],ipp,w,ypv0,ypv2, kapaA); 
  
  
  out[0] = rh;
  out[1] = dfidw;
  out[2] = MatConst[7];
  out[3] = kapa;
  out[4] = -10000.0;
  out[5] = -10000.0;
 

 /* if (ipp == 10)
  fprintf(Outt,"t= %lf je w == %3.2e dwdrh je %3.2e a kapa je %3.2e \n",t,w,dwdf,kapaA);
 */
}
