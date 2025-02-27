#include "globalt.h"
#include "devriesmat.h"
#include "kunzel.h"
//#include "Math.h"

/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
***********************************************/


static double T0 = 273.15; // reference temperature
static double Rvapd = 462.0;  // specific gas constant of water vapour
static double cwatd = 4180;   // specific heat capacity of pure liquid water
static double cvapd = 1617;   // specific heat capacity of water vapour
static double rold = 1000;  // density of water

devriesmat::devriesmat ()
{}

devriesmat::~devriesmat ()
{}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void devriesmat::matcond (matrix &d,long ri,long ci,long ipp)
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
void devriesmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void devriesmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3);

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

void devriesmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,x3);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,x3);
  if((ri == 0) && (ci == 2))
    k = k13 (x1,x2,x3);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,x3);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,x3);
  if((ri == 1) && (ci == 2))
    k = k23 (x1,x2,x3);

  if((ri == 2) && (ci == 0))
    k = k31 (x1,x2,x3);
  if((ri == 2) && (ci == 1))
    k = k32 (x1,x2,x3);
  if((ri == 2) && (ci == 2))
    k = k33 (x1,x2,x3);
  
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
void devriesmat::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2,x3;
  c = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  if((ri == 0) && (ci == 0))
    c = c11 (x1,x2,x3);
  if((ri == 0) && (ci == 1))
    c = c12 (x1,x2,x3);
  if((ri == 0) && (ci == 2))
    c = c13 (x1,x2,x3);

  if((ri == 1) && (ci == 0))
    c = c21 (x1,x2,x3);
  if((ri == 1) && (ci == 1))
    c = c22 (x1,x2,x3);
  if((ri == 1) && (ci == 2))
    c = c23 (x1,x2,x3);

  if((ri == 2) && (ci == 0))
    c = c31 (x1,x2,x3);
  if((ri == 2) && (ci == 1))
    c = c32 (x1,x2,x3);
  if((ri == 2) && (ci == 2))
    c = c33 (x1,x2,x3);
}



/**
   function reads data and material parameters

   @param in  - input file
*/
void devriesmat::read(XFILE *in)
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

  madripom = 0;
}

void devriesmat::print(FILE *out)
{
  fprintf (out,"\n devries material \n");
}

/*position CORD:     2 - density
                    3 - porosity
                    4 - faktor difusniho odporu
                    5 - kapa
                    6 - sorption izoterm
                    7 - saturated moisture
                    8 - hydraulic conductivity
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

double devriesmat::k11(double x1,double x2,double x3)
{
  double k11;
  double dwl, dwv;
  dwl = coeff_of_water_diffusion_by_moisture_grad(x1,x2,x3);
  dwv = coeff_dwv(x1,x2,x3);
  k11 = dwl + dwv;

  return (k11*rold);
}

double devriesmat::k12(double x1,double x2,double x3)
{
  double k12;

  double dtv, dtl;
  dtl = coeff_of_water_diffusion_by_temperature_grad(x1,x2,x3);
  dtv =coeff_dtv(x1,x2,x3);

  k12 = dtv + dtl;

  return (k12*rold);
}

double devriesmat::k13(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k13;
  k13 = 0.0;
  return (k13);
}

double devriesmat::k21(double x1,double x2,double x3)
{
  double k21;
  double dwl, dwv, lv;
  dwl =coeff_of_water_diffusion_by_moisture_grad(x1,x2,x3);
  dwv = coeff_dwv(x1,x2,x3);
  lv = latent_heat_of_evaporation_of_water(x1,x2,x3);



  k21 =  rold*cwatd*dwl *(x2-T0)+rold*dwv*(cvapd*(x2- T0)+lv);//to je ok
// k21 =  rold*(x2-273.15)*(cwatd*dwl + cvapd*dwv)+rold*lv*dwv ; 
// k21 =  rol*x2*(cwat*dwl + cvap*dwv)+rol*lv*dwv ;
  return (-k21);
}

double devriesmat::k22(double x1,double x2,double x3)
{
  double k22;
  double lambda, dtl, dtv, lv;
  dtl =  coeff_of_water_diffusion_by_temperature_grad(x1,x2,x3);
  dtv = coeff_dtv(x1,x2,x3);
  CorD(10,kd,1.0,0,x1,lambda,a2,a3);//lambda
  lv = latent_heat_of_evaporation_of_water(x1,x2,x3);

// k22 = lambda-rold*(x2-273.15)*(cwatd*dtl+cvapd*dtv)-rold*lv*dtv;
  k22 = lambda-rold*cwatd*dtl*(x2-T0)-rold*dtv*(cvapd*(x2-T0)+lv); //to je ok
  // k22 = lambda-rol*x2*(cwat*dtl+cvap*dtv)-rol*lv*dtv;

  return (k22);
}

double devriesmat::k23(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k23;

 k23 = 0.0;

  return (k23);
}

double devriesmat::k31(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k31;

  k31 = 0.0;

  return (k31);
}

double devriesmat::k32(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k32;

  k32 = 0.0;

  return (k32);
}

double devriesmat::k33(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k33;

  k33 = 0.0;

  return (k33);
}

double devriesmat::c11(double x1,double x2,double x3)
{
  double c11;
  double a, drovdw,rov;

  a = relative_volume_ration_a(x1,x2,x3);
  drovdw = derivation_drov_dw(x1,x2,x3);
  rov =partial_density_of_water_vapor(x1,x2,x3);

 c11 = 1+a*drovdw/rold-rov/rold; //to je ok
// c11 = 1+1/rold*a*drovdw-1/rold*rov;
  return(c11*rold);
 // return(c11);
}

double devriesmat::c12(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c12;
  c12 = 0.0;
  return(c12);
}

double devriesmat::c13(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c13;
  c13 = 0.0;

  return(c13);
}

double devriesmat::c21(double x1,double x2,double x3)
{
  double c21;
  double a,lv,rov,drovdw;

  a = relative_volume_ration_a(x1,x2,x3);
  drovdw = derivation_drov_dw(x1,x2,x3);
  rov =partial_density_of_water_vapor(x1,x2,x3);

  //CorD(9,kd,0,0,x1,cs,a2,a3);// c

  lv =latent_heat_of_evaporation_of_water(x1,x2,x3);

 // c21 = lv*a*drovdw-lv*rov+rol*cwat*(x2)+cvap*(x2)*(a*drovdw-rov);
  c21 = rold*cwatd*(x2-T0)+a*drovdw*(cvapd+lv)-rov*(cvapd*(x2-T0)+lv);//toto je OK
  
//c21 = lv*a*drovdw-lv*rov-rold*cwatd*(x2-273.15)+cvapd*(x2-283.15)*(a*drovdw-rov);

  return(c21);
}

double devriesmat::c22(double x1,double x2,double x3)
{
  double c22;
  double ros,cstar, cs,a,rov;

a = relative_volume_ration_a(x1,x2,x3);
  rov =partial_density_of_water_vapor(x1,x2,x3);
  CorD(2,kd,0,0,x1,ros,a2,a3);// bulk density
  CorD(9,kd,0,0,x1,cs,a2,a3);// bulk density
  cstar = specific_heat_capacities_star(x1,x2,x3);

// c22 =ros*cstar;
  c22 =ros*cs+rold*x1*cwatd+a*rov*cvapd;//toto je OK
  //here will be something??!!

  return(c22);
}

double devriesmat::c23(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c23;
  c23 = 0.0;

  return(c23);
}

double devriesmat::c31(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c31;
  c31 = 0.0;
  return(c31);
}

double devriesmat::c32(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c32;
  c32 = 0.0;
  return(c32);
}

double devriesmat::c33(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c33;
  c33 = 0.0;

  return(c33);
}


/* Function calculates all auxiliary data - optional */
void devriesmat::auxiliarydata (double /*x1*/,double /*x2*/,double /*x3*/)
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
double devriesmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_trc,x1,x2,x3;
  new_trc = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
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
double devriesmat::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_nodval,x1,x2,x3;
  new_nodval = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
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
double devriesmat::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double flux,x1,x2,x3;
  flux = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
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
double devriesmat::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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
double devriesmat::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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
double devriesmat::get_transmission_flux_11(double bv,double x1,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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

double devriesmat::get_othervalue(long compother,long /*ipp*/, double x1,double x2,double x3)
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
void devriesmat::print_othervalue_name(FILE *out,long compother)
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

void devriesmat::CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2)
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

double devriesmat::saturation_water_vapour_pressure(double /*x1*/, double x2, double /*x3*/)
{
  return (exp(23.5771 - 4042.9/(x2 - 37.58)));   //TODO: Add your source code here
}

double devriesmat::derivation_saturation_water_vapour_pressure_temperature(double /*x1*/, double x2, double /*x3*/)
{
  return ((4042.9 * exp(23.5771-4042.9/(x2 - 37.58)))/((x2-37.58)*(x2-37.58)));   //TODO: Add your source code here
}

double devriesmat::partial_water_vapour_pressure_function(double x1, double x2, double x3)
{
  double relhum, pvs;

  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  sorption_izoterms_giva_data(0,x1,x2,x3,relhum,a1);

  return(pvs*relhum);
        //TODO: Add your source code here
}


void devriesmat::sorption_izotherm_derivation(double x1, double /*x2*/, double /*x3*/, double & derfi)
{
 int i = 0;
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
 
double devriesmat::derivation_specific_internal_energy_of_water_vapour(double /*x1*/, double /*x2*/, double /*x3*/)
{
 return(cvapd );    //TODO: Add your source code here
}

double devriesmat::get_rel_hum(double w)
{
   double relhum;
   sorption_izoterms_giva_data(0,w,0,0,relhum,a1);    //TODO: Add your source code here
   return (relhum);
}

double devriesmat::sortpion_isotherm_root_shifted(double /*x1*/, double w_hyg, double rh_hyg)
{
  return(w_hyg/(1-sqrt(1-rh_hyg)));
}

void devriesmat::give_data_si_root_dfidw(double x1, double /*x2*/, double /*x3*/, double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw)
{
     //double dfdw,
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

double devriesmat::get_moisture(double rh)
{
  int s;
 double moist,hmrh, mhmc, smc, u, a, n;

       if (rh < 0) rh = 0;

       CorD(6,s,0,0,0,a1,a2,a3);      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
          case 2:   CorD(6,s,1,0,rh,moist,a1,a3);
          break;
          case 30:  CorD(6,s,1,0,rh,a1,hmrh,a3);
                    CorD(7,s,0,0,rh,smc,a2,a3);
                    CorD(6,s,0,0,rh,mhmc,hmrh,a3);
                  //  kun->give_data(rh,mhmc, smc, hmrh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,0,0,0,u,a,n);
                  //  moist = kun->si_kk_hansen(rh, 0,u,a,n);
          break ;
          }
 return(moist);
    //TODO: Add your source code here
}
void devriesmat::sorption_izoterms_giva_data(int kod,double x1, double x2, double x3, double & fi, double & dfdw)
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
          case 2:   CorD(6,s,1,1,x1,fi,a2,a3);
                    sorption_izotherm_derivation(x1,x2,x3,dfdw);
          break;
          case 30:  //CorD(6,s,1,x1,a1,hmrh,a3);
                      CorD(6,s,0,0,x1,w_hyg,rh_hyg,a3);
                      hvezda =sortpion_isotherm_root_shifted(x1,w_hyg ,rh_hyg);
                      CorD(7,s,0,0,x1,smc,a2,a3);
                    switch (kod)
                         {
                         case 0:
                                   give_data_si_root_fi(x1,x2,x3,w_hyg, smc, rh_hyg,hvezda,fi);
                         break;
                         case 1:   give_data_si_root_dfidw(x1,x2,x3,rh_hyg,smc,w_hyg,hvezda,dfdw);
                         break;
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
void devriesmat::give_data_si_root_fi(double x1, double /*x2*/, double /*x3*/, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh)
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

double devriesmat::pressure_head(double x1, double x2, double /*x3*/)
{
  double relhum;
   sorption_izoterms_giva_data(0,x1,0,0,relhum,a1);
 if (relhum<=0.1)
          {
           return(1e6);
  }
     else {
          return (Rvapd*x2*log(relhum)/9.80665);   //TODO: Add your source code here
          }



// return (Rvapd*x2*log(relhum)/9.80665); 
}

double devriesmat::derivation_pressure_head(double x1, double x2, double x3)
{
     double relhum, dfdw;
     sorption_izoterms_giva_data(0,x1,0,0,relhum,a1);
     sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfdw);
if (relhum<=0.1)
          {
           return(Rvapd*x2/(9.80665*0.1)*1/dfdw);
          }
     else {
       return (Rvapd*x2/(9.80665*relhum)*1/dfdw);   //TODO: Add your source code here

       }
 //return (Rvapd*x2/(9.80665*relhum)*1/dfdw);   //TODO: Add your source co
}


double devriesmat::relative_volume_ration_a(double x1, double /*x2*/, double /*x3*/)
{
   double pi;
   CorD(3,kd,0,0,x1,pi,a2,a3);//pir porovitost
   return (pi-x1);  //TODO: Add your source code here
}


double devriesmat::surface_tension(double /*x1*/, double x2, double /*x3*/)
{
double teplota;
teplota = x2-273.15;

     return ((-0.1647*teplota + 75.882)/1000);//TODO: Add your source code here
}

double devriesmat::derivation_surface_tension_on_temperature(double /*x1*/, double /*x2*/, double /*x3*/)
{
 return (-0.0001647);    //TODO: Add your source code here
}


double devriesmat::coeff_of_water_diffusion_by_moisture_grad(double x1, double x2, double x3)
{
 double dhdw, kl;

 CorD(8,kd,1,0,x1,kl,a2,a3);//hydraulic conductivity
 dhdw = derivation_pressure_head(x1,x2,x3);
 return(kl*dhdw);
      //TODO: Add your source code here
}

double devriesmat::coeff_of_water_diffusion_by_temperature_grad(double x1, double x2, double x3)
{
 double kl, sigm, h,dsigmdt;
 CorD(8,kd,1,0,x1,kl,a2,a3);//hydraulic conductivity
 h = pressure_head(x1,x2,x3);
 sigm = surface_tension(x1,x2,x3);
 dsigmdt = derivation_surface_tension_on_temperature(x1,x2,x3);
 return (kl*h*1/sigm*dsigmdt);

     //TODO: Add your source code here
}

double devriesmat::coeff_dwv(double x1, double x2, double x3)
{
  double rol, mi, tor, a,d,p,pv,pvs,dfdw;
  rol = 1000;
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  tor = 1/mi;
  a =relative_volume_ration_a(x1,x2,x3);
  d = diffusion_coefficient_of_water_vapor_in_air(x1,x2,x3);// pozor dodelat....
  p = 101325;
  pv = partial_water_vapour_pressure_function(x1,x2,x3);
  pvs = saturation_water_vapour_pressure(x1,x2,x3);

  sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfdw);

  return (1/rol*tor*a*d*p/(p-pv)*1/(Rvapd*x2)*pvs *dfdw);

     //TODO: Add your source code here
}

double devriesmat::coeff_dtv(double x1, double x2, double x3)
{
  double rold, mi, tor, a,d,p,pv,pvs,dfdw, relhum, dpvsdt, dfdt;
  rold = 1000;
  CorD(4,kd,0,0,x1,mi,a2,a3);// mi
  tor = 1/mi;
  a =relative_volume_ration_a(x1,x2,x3);
  d = diffusion_coefficient_of_water_vapor_in_air(x1,x2,x3);// pozor dodelat....
  p = 101325;
  pv = partial_water_vapour_pressure_function(x1,x2,x3);
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  dpvsdt = derivation_saturation_water_vapour_pressure_temperature(x1,x2,x3);
  sorption_izoterms_giva_data(0,x1,0,0,relhum,a1);
  sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfdw);
  dfdt = 0.0;

  return(1/rold*tor*a*d*p/(p-pv)*1/(Rvapd*x2)*(-relhum*pvs/x2+relhum*dpvsdt+pvs *dfdt));
   //TODO: Add your source code here
}

double devriesmat::specific_heat_capacities_star(double x1, double x2, double x3)
{
  double cs, cl, cv, ros, rov, a;
  CorD(9,kd,0,0,x1,cs,a2,a3);// c
  CorD(2,kd,0,0,x1,ros,a2,a3);// bulk density
  a =relative_volume_ration_a(x1,x2,x3);
  cl = cwatd;
  cv = cvapd;
  rold = 1000;
  rov = partial_density_of_water_vapor(x1,x2,x3);

  return (cs +rold/ros*cl*x1+rov/ros*a*cv);

     //TODO: Add your source code here
}

double devriesmat::partial_density_of_water_vapor(double x1, double x2, double x3)
{
 double pv;
 pv = partial_water_vapour_pressure_function(x1,x2,x3);
 return (pv/(Rvapd*x2));    //TODO: Add your source code here
}

double devriesmat::latent_heat_of_evaporation_of_water(double /*x1*/, double x2, double /*x3*/)
{
  return (2.5008e6)*pow((273.15/x2),(0.167+x2*3.67e-4));
}


double devriesmat::derivation_drov_dw(double x1, double x2, double x3)
{
  double pvs,dfdw;
  pvs = saturation_water_vapour_pressure(x1,x2,x3);
  sorption_izoterms_giva_data(1,x1,x2,x3,a1,dfdw);

  return (1/(Rvapd*x2)*pvs*dfdw);
     //TODO: Add your source code here
}

double devriesmat::diffusion_coefficient_of_water_vapor_in_air(double x1, double x2, double /*x3*/)
{
   double mi;
   CorD(4,kd,0,0,x1,mi,a2,a3);// mi
   return(2.3e-5 * pow(x2/273,1.81));  //TODO: Add your source code here
}
