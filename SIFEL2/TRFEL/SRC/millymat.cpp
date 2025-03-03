#include "globalt.h"
#include "millymat.h"
#include "math.h"
#include "globmatt.h"

/**********************************************
*                                             *
* General material model for 2 media transfer *
*      4 independent parameters               *
***********************************************/

static double T0 = 293.15; // reference temperature
static double M = 0.018; // doplnit
static double R = 8.314; // univerzalni plynnova konstanta.. kontrola
static double rl = 1000;  // density of water
static double p = 101325; // tlak
static double g = 9.81; // gravitacni konstanta


//static double Rvapd = 462.0;  // specific gas constant of water vapour
static double cwatd = 4180;   // specific heat capacity of pure liquid water
static double cvapd = 1617;   // specific heat capacity of water vapour
//double rold = 1000;  // density of water

//static kunmat *kun;

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


millymat::millymat ()
{}

millymat::~millymat ()
{}


 

void millymat::values_correction (vector &nv)
{
  //  relative humidity
  nv[0] = nv[0];
  nv[1] = nv[1];

 // relhum_check(nv[0],nv[1],0);
}

/**
   function computes conductivity matrix of the material
   in the required integration point

   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void millymat::matcond (matrix &d,long ri,long ci,long ipp)
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
void millymat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;

  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2, ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2, ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2, ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2, ipp);

  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void millymat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;

  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2, ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2, ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2, ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2, ipp);

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

void millymat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;

  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2, ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2, ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2, ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2, ipp);

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
void millymat::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2;
  c = 0.0;

  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];

  if((ri == 0) && (ci == 0))
    c = c11 (x1,x2, ipp);
  if((ri == 0) && (ci == 1))
    c = c12 (x1,x2, ipp);

  if((ri == 1) && (ci == 0))
    c = c21 (x1,x2, ipp);
  if((ri == 1) && (ci == 1))
    c = c22 (x1,x2, ipp);

}



/**
   function reads data and material parameters

   @param in  - input file
*/
void millymat::read (XFILE *in)
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

  //madripom = 0;
  //CorD(2,kd,0,0,rho_m,a2,a3);
}


/******************************************************
Conductivity and capacity terms for C and K matrices
*******************************************************/

double millymat::k11(double x1,double x2,long ipp)
{
  double k11, rs, Kl, Kv;

  CorD(2,kd,0,0,x1,rs,a2,a3);// ros

  Kl =  Tm->ip[ipp].eqother[5];
  Kv =coeff_Kv(x1,x2,ipp);

 // k11 = rl/rs*(Kl+Kv);
  k11 = rl*(Kl+Kv);
  return (k11);
}

double millymat::k12(double x1,double x2,long ipp)
{
  double k12, rs,Dtv;

  CorD(2,kd,0,0,x1,rs,a2,a3);// ros
  Dtv =coeff_Dtv(x1,x2,ipp);

 // k12 = rl/rs*Dtv;
  k12 = rl*Dtv;
 // k12  = 0.0;
  return (k12);
}


double millymat::k21(double x1,double x2,long ipp)
{
  double k21, rs,Kl, Kv, hl, hv;

  CorD(2,kd,0,0,x1,rs,a2,a3);// ros
  Kl =  Tm->ip[ipp].eqother[5];
  Kv =coeff_Kv(x1,x2,ipp);
  hl =specific_enthalpy_of_liquid_phase(x1,x2);
  hv = specific_enthalpy_of_gaseous_phase(x1,x2);

//  k21 = 1/rs*(-hl*rl*Kl-hv*rl*Kv);
  k21 = (-hl*rl*Kl-hv*rl*Kv);
  //k21 = 0.0;
  return (k21);
}

double millymat::k22(double x1,double x2,long ipp)
{
  double k22, rs, lambda, Dtv, hv;
  CorD(2,kd,0,0,x1,rs,a2,a3);// ros
  lambda =  Tm->ip[ipp].eqother[6];
  hv = specific_enthalpy_of_gaseous_phase(x1,x2);
  Dtv =coeff_Dtv(x1,x2,ipp);

 // k22 = 1/rs*(lambda-Dtv*hv*rl);
  k22 = (lambda-Dtv*hv*rl);
 //   k22 = (lambda+Dtv*hv*rl);
  fprintf(Outt,"k22 je %3.6e pro rh= %lf a t= %lf \n",k22,x1,x2);
  return (k22);
}

double millymat::c11(double x1,double x2,long ipp)
{
  double c11, rs,A, B, duldpsi;
   CorD(2,kd,0,0,x1,rs,a2,a3);// ros
  A =coeef_zaporneA(x1,x2, ipp);
  B = coeef_B(x1,x2, ipp);
  duldpsi =  Tm->ip[ipp].eqother[8];

  c11 = (1+A)*duldpsi + B;
  c11 = c11 *rs;
  return(c11);
}

double millymat::c12(double x1,double x2,long ipp)
{
  double c12,rs, A, C, duldt;

  CorD(2,kd,0,0,x1,rs,a2,a3);// ros
   A =coeef_zaporneA(x1,x2, ipp);
  C = coeef_C(x1,x2, ipp);
  duldt =  Tm->ip[ipp].eqother[4];

  c12 = (1+A)*duldt + C;
  c12 = c12 *rs;
//  c12 = 0.0;
  return(c12);
}


double millymat::c21(double x1,double x2,long ipp)
{
  double c21, A, B, duldpsi, hl, hv, rs;
  A =coeef_zaporneA(x1,x2, ipp);
  B = coeef_B(x1,x2, ipp);
  duldpsi =  Tm->ip[ipp].eqother[8];
  hv = specific_enthalpy_of_gaseous_phase(x1,x2);
  hl = specific_enthalpy_of_liquid_phase(x1,x2);

  c21 = (hl+hv*A)*duldpsi+B*hv;
 //  c21 = 0.0;
  CorD(2,kd,0,0,x1,rs,a2,a3);//rs .. hustota materialu
  c21 = c21 * rs;
  return(c21);
}

double millymat::c22(double x1,double x2,long ipp)
{
  double c22, A, C, duldt, hl, hv, cs, ul, uv, rs;

  A =coeef_zaporneA(x1,x2, ipp);
  C = coeef_C(x1,x2, ipp);
  duldt =  Tm->ip[ipp].eqother[4];
  ul =  Tm->ip[ipp].eqother[0];
  hv = specific_enthalpy_of_gaseous_phase(x1,x2);
  hl = specific_enthalpy_of_liquid_phase(x1,x2);
  CorD(9,kd,0,0,x1,cs,a2,a3);// rcs - opravit
  uv = gaseous_moisture_content_by_mass(x1,x2,ipp);

  c22 = (hl+A*hv)*duldt+(cs+cwatd*ul+cvapd*uv+hv*C);
   CorD(2,kd,0,0,x1,rs,a2,a3);//rs .. hustota materialu
   c22 = c22*rs;
  fprintf(Outt,"c22 je %3.6e pro rh= %lf a t= %lf \n",c22,x1,x2);
  return(c22);
}



/* Function calculates all auxiliary data - optional */
void millymat::auxiliarydata (double /*x1*/,double /*x2*/)
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
double millymat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,x1,x2;
  new_trc = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_11(x1,x2,bc,ipp);
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
double millymat::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
{

  double new_nodval,x1,x2;
  new_nodval = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_11(nodval,x1,x2,bc,ipp);
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
double millymat::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
{

  double flux,x1,x2;
  flux = 0.0;
  
  x1 = nodalval (nn,0);
  x2 = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_11(nodval,x1,x2,bc,ipp);
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
double millymat::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double millymat::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double millymat::get_transmission_flux_11(double bv,double x1,double /*x2*/,long bc,long /*ipp*/)
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

double millymat::get_othervalue(long compother,long /*ipp*/, double x1,double /*x2*/)
{
  double other;

  switch (compother){
  case 0:{//first unknown
    other = relative_humidity_psi(x1,293.15);//should be changed
      break;
  }
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
  case 11:
  case 12:
  break;
  default:{
   // fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);
}

/**
     function prints names of all variables in nodes
     @param out       - output file
     @param compother - number of other components
*/
void millymat::print_othervalue_name(FILE *out,long compother)
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

double millymat::saturation_water_vapour_pressure(double /*x1*/, double x2)
{
  return (exp(23.5771 - 4042.9/(x2 - 37.58)));   //TODO: Add your source code here
}

double millymat::derivation_saturation_water_vapour_pressure_temperature(double /*x1*/, double x2)
{
  return ((4042.9 * exp(23.5771-4042.9/(x2 - 37.58)))/((x2-37.58)*(x2-37.58)));   //TODO: Add your source code here
}

double millymat::relative_volume_ration_a(double x1, double /*x2*/, double ul)
{
   double pi, rs;
   CorD(3,kd,0,0,x1,pi,a2,a3);//pir porovitost
   CorD(2,kd,0,0,x1,rs,a2,a3);//rs .. hustota materialu
  // rl = 1000; // udelat zavislost na teplote
   // z retencni krivky urcime ul, zkontrolovat jednotky

   return (pi-rs/rl*ul);  //TODO: Add your source code here
}
double millymat::relative_humidity_pc(double /*x1*/, double x2)
{
	double pc, fi;
       // x2 = T0;
        // pc .. saci tlak z retencni krivky
        //rl ...hustotat vody v zavislosti nateplote
        pc = 1000000000000.0;// nevim,kdo to tam dal, ale je potreba vyresit, odkuz vezmu tuto hodnotu
        rl = 1000;
        fi = exp(-pc*M/(rl*R*x2));
        //fi = exp(-pc*M/(rl*R*293.15));
        if (fi>1.0) {
            fi = 1.0;
        }
        return (fi);
}
  double millymat::relative_humidity_psi(double x1, double x2)
{
	double fi;
        //rl ...hustotat vody v zavislosti nateplote
        //x2 = T0;

        fi = exp(x1*M*g/(R*x2));
        return (fi);
}
double millymat::coeef_zaporneA (double x1, double x2, long ipp)
{
	double fi, pvs;
         fi = Tm->ip[ipp].eqother[2];
        pvs = saturation_water_vapour_pressure(x1,x2);
        return (-1/rl*fi*pvs*M/(R*x2));

}

double millymat::coeef_B (double x1, double x2, long ipp)
{
	double fi, pvs, a, rs;
        fi = Tm->ip[ipp].eqother[2];
        pvs = saturation_water_vapour_pressure(x1,x2);
        a = Tm->ip[ipp].eqother[7];
        CorD(2,kd,0,0,x1,rs,a2,a3);//rs .. hustota materialu
        return (fi*pvs*a/rs*(M/(R*x2))*(M/(R*x2)));

}
double millymat::coeef_C (double x1, double x2, long ipp)
{
	double fi, pvs, a, rs, dpvsdt;

        fi = Tm->ip[ipp].eqother[2];
        pvs = saturation_water_vapour_pressure(x1,x2);
        dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2);
        a = Tm->ip[ipp].eqother[7];

        CorD(2,kd,0,0,x1,rs,a2,a3);//rs .. hustota materialu


       return (1/rs*M/(R*x2)*a*pvs*fi*(-x1*g*M/(R*x2*x2)+1/pvs*dpvsdt-1/x2));
   //      return (M/(R*x2)*a*pvs*fi*(-x1*g*M/(R*x2*x2)+1/pvs*dpvsdt-1/x2));

}
double millymat::diffusion_coefficient_of_water_vapor_in_air(double x1, double x2)
{
   double mi;
   CorD(4,kd,0,0,x1,mi,a2,a3);// mi
   return(2.3e-5 * pow(x2/273.15,1.81));  //TODO: Add your source code here
}

double millymat::partial_water_vapour_pressure_function(double x1, double x2, long ipp)
{
  double fi, pvs;

  pvs = saturation_water_vapour_pressure(x1,x2);
  //sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1);
   fi = Tm->ip[ipp].eqother[2];
  return(pvs*fi);
        //TODO: Add your source code here
}

double millymat::coeff_Kv (double x1, double x2, long ipp)
{
   double  alf,a, mi, D, pv, fi,pvs;

   CorD(4,kd,0,0,x1,mi,a2,a3);// mi
   alf =1/mi; // alfa neboli tortuozita

   a = Tm->ip[ipp].eqother[7];
   D = diffusion_coefficient_of_water_vapor_in_air(x1,x2);// pozor dodelat ->kontrola....
   fi = Tm->ip[ipp].eqother[2]; // relativni vlhkost z nezname psi (x1)
   pv = partial_water_vapour_pressure_function(x1,x2,ipp);
   pvs = saturation_water_vapour_pressure(x1,x2);

   return(1/rl*alf*a*D*p/(p-pv)*(M/(R*x2))*(M/(R*x2))*g*pvs*fi );
}

double millymat::coeff_Dtv (double x1, double x2, long ipp)
{
   double  alf,a, mi, D, pv, fi,pvs,dpvsdt;

   CorD(4,kd,0,0,x1,mi,a2,a3);// mi
   alf =1/mi; // alfa neboli tortuozita

   a = Tm->ip[ipp].eqother[7];
   D = diffusion_coefficient_of_water_vapor_in_air(x1,x2);// pozor dodelat ->kontrola....
   fi = Tm->ip[ipp].eqother[2]; // relativni vlhkost z nezname psi (x1)
   pv = partial_water_vapour_pressure_function(x1,x2,ipp);
   pvs = saturation_water_vapour_pressure(x1,x2);
   dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2);

   return(1/rl*alf*a*D*p/(p-pv)*(M/(R*x2))*fi*(-pvs/x2+dpvsdt-M/(R*x2*x2)*g*x1*pvs));
}
double  millymat::gaseous_moisture_content_by_mass (double x1, double x2, long ipp)
{
// v modelu uv
	double pvs, ul, pi, rs;

         pvs = saturation_water_vapour_pressure(x1,x2);
         CorD(3,kd,0,0,x1,pi,a2,a3);//pir porovitost
         CorD(2,kd,0,0,x1,rs,a2,a3);//ros
         ul = Tm->ip[ipp].eqother[0]; // z retencni krivky

         return( (pvs*M/(R*x2))*(pi/rs-ul/rl)*exp(x1*g*M/(R*x2)));

}

double millymat:: specific_enthalpy_of_porous_matrix (double x1, double x2)
{
	double cs;
        CorD(9,kd,0,0,x1,cs,a2,a3);//cs - opravit

        return(cs*(x2-T0));
}

double millymat:: specific_enthalpy_of_liquid_phase (double /*x1*/, double x2)
{
	return(cwatd*(x2-T0));
}
double millymat:: specific_enthalpy_of_gaseous_phase (double x1, double x2)
{
	double Lv;
        Lv =latent_heat_of_evaporation_of_water(x1,x2);
	return(cvapd*(x2-T0)+Lv);
}
double millymat::latent_heat_of_evaporation_of_water(double /*x1*/, double x2)
{
  return (2.5008e6)*pow((273.15/x2),(0.167+x2*3.67e-4));
}

double millymat:: derivation_water_retention_curve_by_pressure_head (double /*x1*/, double /*x2*/)
{
  double dulpsi;
   dulpsi = 0.0; // tady taky nevim, tak nastaveno na 0.
  return (dulpsi);
}

double millymat:: derivation_water_retention_curve_by_temperature (double /*x1*/, double /*x2*/)
{
  double dut;
  dut = 0.0;
  return (dut);
}
double millymat::give_ul (double x1, double x2, double &fi, double &dmoistdrh )
{
    double ul;
    // tady prepocitavam psi na relativni vlhkost z duvodu rozdeleni a zapojeni SI a RK


    fi = relative_humidity_psi(x1,x2);
    if (fi>1.0) {
        fi = 1.0;
    }
    if (fi>0.976) {
       ul =give_ul_retention_curve(x1,x2,fi, dmoistdrh);

    }
    else {
    	ul = give_ul_sorption_izothemrs(x1,x2,fi, dmoistdrh);
      //  ul = 0.0;
    }

  return(ul);
}
double millymat::give_ul_sorption_izothemrs (double x1, double x2,double rh, double &dmoistdrh )
{
 	double  rs;
 	// tady musi byt klasicke hledani SI

        int s;
 	double moist,hmrh, mhmc, smc, u, a, n;

        if (rh < 0) rh = 0;

       	CorD(6,s,0,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       	switch (s)
          {
	   case 2:   {
                      CorD(6,s,1,0,rh,moist,a1,a3);
                      dmoistdrh =derivation_sorption_izoterm_data(rh,x2, 6);
             }
          break;
          case 30:  CorD(6,s,1,0,rh,a1,hmrh,a3);
                    CorD(7,s,0,0,rh,smc,a2,a3);
                    CorD(6,s,0,0,rh,mhmc,hmrh,a3);
                    //kun->give_data(rh,mhmc, smc, hmrh,moist);
                    dmoistdrh = DerivativeOfTheSorptionIsotherm (x1,x2,rh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,0,0,0,u,a,n);
	    //moist = kun->si_kk_hansen(rh, 0,u,a,n);
                    //dmoistdrh = kun->sorptionizothermDerivation(rh, x2);
          break ;
          }
        CorD(2,kd,0,0,moist,rs,a2,a3);//ros

        return (moist*rl/rs);
}

double millymat::give_ul_retention_curve (double x1, double x2, double /*rh*/, double &dmoistdrh )
{
        double moist, rs, pc;
        int s;
 	//fi = relative_humidity_psi(x1,x2);
        pc = capilar_pressure (x1,x2);

        // tady musi byt klasicke hledani SI
        CorD(11,s,0,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN
         // tady je to blbe....


       	switch (s)
          {
	   case 2:   {
                      CorD(11,s,1,0,pc,moist,a1,a3);
                      dmoistdrh =derivation_sorption_izoterm_data(pc,x2, 11);
             }
          break;
        }

        CorD(2,kd,0,0,moist,rs,a2,a3);//ros

        return (moist*rl/rs);

}
double millymat::capilar_pressure (double x1, double /*x2*/)
{
	return(-x1*g*rl);
}


double millymat::ul_to_w (double ul)
{
	double rs;
        CorD(2,kd,0,0,ul,rs,a2,a3);//ros
        return(ul*rs/rl);
}

double millymat::derivation_sorption_izoterm_data(double x1, double /*x2*/, int kod)
{
  // kod ... rozdeluje na SI = 6 a Ret krivku= 11
  int i = 1;
  double derbi;
     if ( x1 < 0)
          {
          x1 = 0;
          }

     while (x1 > MatData[kod][0][i]){
               i++;
               }
 derbi = (MatData [kod][1][i] - MatData [kod][1][i-1])/(MatData [kod][0][i] - MatData [kod][0][i-1]);   //TODO: Add your source code here
 return(derbi);
}

double millymat::DerivativeOfTheSorptionIsotherm(double x1, double x2,double rh,double &moistakt1)
{
// pro funkci ROOT
 double mhmc, mhrh, rs,rh2,moistakt2, ul1, ul2, psi1, psi2;
 double hvezdicka;

        if (rh < 0) rh = 0;

         CorD(6,kd,0,0,rh,mhmc,mhrh,a3);

  	 hvezdicka = mhmc/(1-sqrt(1-mhrh));
         moistakt1 = (1-sqrt(1-rh)) * hvezdicka ;
         rh2 = rh +.001;
         moistakt2 = (1-sqrt(1-rh2)) * hvezdicka ;

         CorD(2,kd,0,0,moistakt1,rs,a2,a3);//ros

         ul1 = (moistakt1*rl/rs);
         ul2 = (moistakt2*rl/rs);

         //x2 = T0;

         psi1 = x1;
         psi2 = R*x2/(g*M)*log(rh2);

        return ((ul2-ul1)/(psi2-psi1));

}

void millymat::CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2)
{

  switch (MatChar[cislochar])
     {
     case 0:
     break;
     case 1: y = MatConst[cislochar];
     break;
     case 2:   {
                int i = 1;

                int psl, prad;
                prad  = int(MatData[cislochar][0][0]);
                psl = int(MatData[cislochar][1][0]);
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
     }
}


void millymat::initvalues (long ipp,long /*ido*/)
{
  double x1,x2;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
}

void millymat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
  Tm->ip[ipp].eqother[5]=out[5];
  Tm->ip[ipp].eqother[6]=out[6];
  Tm->ip[ipp].eqother[7]=out[7];
  Tm->ip[ipp].eqother[8]=out[8];

}

void millymat::give_values (long ipp,double *av, double */*pv*/, double *eq)
{
  av[0] = Tm->ip[ipp].av[0];
  av[1] = Tm->ip[ipp].av[1];
  
  
  eq[0] = Tm->ip[ipp].eqother[0];
  eq[1] = Tm->ip[ipp].eqother[1];
  eq[2] = Tm->ip[ipp].eqother[2];
  eq[3] = Tm->ip[ipp].eqother[3];
  eq[4] = Tm->ip[ipp].eqother[4];
  eq[5] = Tm->ip[ipp].eqother[5];
  eq[6] = Tm->ip[ipp].eqother[6];
  eq[7] = Tm->ip[ipp].eqother[7];
  eq[8] = Tm->ip[ipp].eqother[8];
  
}

void millymat::aux_values (long /*ipp*/,double *in,double */*inp*/, double */*ineq*/,double *out)
{
  //  double x1,x2,x1p, x2p, ypv1,ypv2,ypv3,ypv4,ypv5, ypv0;
  double x1,x2;
  
  double ul, w, rh, dudpsi, dudT, kl, lambda,a;
  
  x1 = in[0];
  x2 = in[1];
  /*
    x1p = inp[0];
    x2p = inp[1];
  */
  /* ypv0 =ineq[0];
     ypv1 =ineq[1];
     ypv2 =ineq[2];
     ypv3 =ineq[3];
     ypv4 =ineq[4];
     ypv5 =ineq[5];
     //ypv5 =ineq[6];
     */
  ul = give_ul(x1,x2, rh, dudpsi);
  w = ul_to_w(ul);
  dudT = derivation_water_retention_curve_by_temperature (x1, x2);
  // rh = relative_humidity_psi(x1,x2);
  CorD(8,kd,0,0,w,kl,a2,a3);//hydraulicka vodivost
  CorD(10,kd,0,0,w,lambda,a2,a3);//lambda
  a = relative_volume_ration_a(x1,x2,ul);
  if(rh<0.97)
    {
      ul = 0.0;
    }
  
  out[0] = ul;
  out[1] = w;
  out[2] = rh;
  out[8] = dudpsi;
  out[4] = dudT;
  out[5] = kl;
  out[6] = lambda;
  out[7] = a;
  out[3] = MatConst[7];
  
  /*
    double t;
    t= Tp->time;
    if (ipp == 10)
    fprintf(Outt,"t= %lf je rh == %3.2e ul je %3.2e a dudpsi je %3.2e \n",t,rh,ul,dudpsi);
  */
  
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
    18 - 19 - none
  */
}
