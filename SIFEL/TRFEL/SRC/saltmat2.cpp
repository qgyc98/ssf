#include "globalt.h"
#include "saltmat2.h"
#include "globmatt.h"

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
                    14 Dcoef
                    15 - binding isotherm
                    16 - cfmax
                    17 ws
                    18 - 19 - none

*/

saltmat2::saltmat2 ()
{
  mw = 0.01801528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8.31441; //universal gas constant J.mol-1.K-1

}

saltmat2::~saltmat2 ()
{}



 
/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void saltmat2::matcond (matrix &d,long ri,long ci,long ipp)
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
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void saltmat2::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    //matcond1d (d,ri,ci,ipp);//1D
    //break;
  }
  case 2:{
    matcond2d2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    //matcond3d (d,ri,ci,ipp);//3D
    //break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void saltmat2::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
  
  //  check of salt concentration
  CorD(16,kd,0,x1,cfmax,a2,a3);
  if (x2 > cfmax){
    x2 = cfmax;
  }
  
    
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

void saltmat2::matcond2d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  d[0][0] = 0.0;   d[0][1] = 0.0;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void saltmat2::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  //  check of salt concentration
  CorD(16,kd,0,x1,cfmax,a2,a3);
  if (x2 > cfmax){
    x2 = cfmax;
  }

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

void saltmat2::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2,x3,cfmax;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  x3 = Tm->ip[ipp].av[2];
    
  //  check of salt concentration
  CorD(16,kd,0,x1,cfmax,a2,a3);
  if (x2 > cfmax){
    x2 = cfmax;
  }

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
void saltmat2::matcap (double &c,long ri,long ci,long ipp)
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


void saltmat2::values_correction (vector &nv)
{
  double cfmax;
  
  CorD(16,kd,0,1,cfmax,a2,a3);
  if (nv[1] > cfmax){
    nv[1] = cfmax;
  }
}



/**
   function reads data and material parameters

   @param in  - input file
*/
void saltmat2::read (XFILE *in)
{
int i, j;
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
	case 32:  if( MatChar[i-1] ==1)
	    {
	      xfscanf(in, "%lf %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	  else
	    {
	      xfscanf(in, "%lf  %lf",&MatFunce[i][0],&MatFunce[i][1]);
	    }
	  break;
	}
    }
  
}


/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

double saltmat2::k11(double x1,double /*x2*/,double /*x3*/)
{
  double k11;
  double kapak, ps, delta, dfdw;
  double p = 101325;
  double a,u,n;
  int kod;
  CorD(5,kd,1,x1,kapak,a2,a3);
 // ps = pgws(t);
  ps = pgws(294.15);
  //delta = permeabilitavodnipary(t,p);
  delta = permeabilitavodnipary(294.15,p);
  //sisotherm(isi1, x1,fi, dfdw);
  //CorD(6,0,1,x1,0,dfdw,0);
  CorD(6,kod,1,x1,u,a,n);

  switch (kod)
          {
          case 2:   //fi = u;
                    dfdw = a;
          break;
          case 30:
          break;
          case 31:  {
                    dfdw = (u/(a*x1*n))*pow((1-(log(x1))/a),(-1-(1/n)));
                   // fi = u*pow((1-(log(x1))/a),(-1/n));
                  //   printf ("\n w je %lf  fi je %lf a derivace je %lf",x1,fi, dfdw);
                    }
          break;
          }
  k11 = 1000*kapak + ps*delta*dfdw;
  return (k11);
}

double saltmat2::k12(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k12;

  k12 =0.0; //here will be something??!!

  return (k12);
}

double saltmat2::k13(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k13;

  k13 = 0.0;

  //here will be something??!!

  return (k13);
}

double saltmat2::k21(double x1,double x2,double /*x3*/)
{
  double k21;
  double kapak;

  CorD(5,kd,1,x1,kapak,a2,a3);

  k21 = kapak*x2*1;
  
  return (k21);
}

double saltmat2::k22(double x1,double x2,double /*x3*/)
{
  double k22;
  double dcoef;

   // musim doladit a upresnit
  CorD(14,kd,1000,x2,dcoef,a2,a3);

  k22 = dcoef*x1;


  return (k22);
}

double saltmat2::k23(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k23;

  k23 = 0.0;
  //here will be something??!!

  return (k23);
}

double saltmat2::k31(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k31;

  k31 = 0.0;

  //here will be something??!!

  return (k31);
}

double saltmat2::k32(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k32;

  k32= 0.0;

  //here will be something??!!

  return (k32);
}

double saltmat2::k33(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double k33;

  k33= 0.0;

  //here will be something??!!

  return (k33);
}

double saltmat2::c11(double x1,double /*x2*/,double /*x3*/)
{
  double c11;
  double ps;
  double fi,dfdw, wsoli, por;
  double u, a, n;
  int kod;

  ps = pgws(294.15);
  CorD(3,kd,1,x1,por,a2,a3);
  CorD(6,kod,1,x1,u,a,n);

  switch (kod)
          {
          case 2:   fi = u;
                    dfdw = a;
          break;
          case 30:
          break;
          case 31:  {
                    dfdw = (u/(a*x1*n))*pow((1-(log(x1))/a),(-1-(1/n)));
                    fi = u*pow((1-(log(x1))/a),(-1/n));
                  //   printf ("\n w je %lf  fi je %lf a derivace je %lf",x1,fi, dfdw);
                    }
          break;
          }

  CorD(17,kd,0,x1,wsoli,a2,a3);

  c11 = 1000 + mw/(gasr*294.15)*ps*(por-x1-wsoli)*dfdw-mw*ps*fi/(gasr*294.15);

  return(c11);
}

double saltmat2::c12(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c12;

 c12 = 0.0;

  return(c12);
}

double saltmat2::c13(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c13;

  c13 = 0.0;

  return(c13);
}

double saltmat2::c21(double x1,double x2,double /*x3*/)
{
  double c21,cfmax;

  //  check of salt concentration
  CorD(16,kd,0,x1,cfmax,a2,a3);
  if (x2 > cfmax){
    c21 = cfmax;

  }
  else{
    c21 = x2;
  }
  
  return(c21);
}

double saltmat2::c22(double x1,double x2,double /*x3*/)
{
  double c22;
    double dCbDcf;
  int k;
  double a,b, Cb;

  if (x2 <= 0) x2 = 1e-10;
  if (x1 <0) x1 = 0;

  CorD(15,k,1000,x2,a,b,a3);

  switch (k)
     {
     case 31: dCbDcf = b*a*pow(x2,(b-1));
     break;
     case 30:   Cb =  1/(1/(a*b*x2)+1/b);
                dCbDcf = Cb*Cb/(a*b*x2*x2);
     break;
     case 2:   binding_izoterm_derivation(x2,dCbDcf);

     break;
     case 0: dCbDcf =0;// pouziti dat
     break;
     }

  c22 = x1 + dCbDcf;
 // c22 = x1;

  return(c22);
}

double saltmat2::c23(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c23;


  c23 = 1.0;

  return(c23);
}

double saltmat2::c31(double x1,double x2,double /*x3*/)
{
  double c31;
  double cfmax;

  CorD(16,kd,0,x1,cfmax,a2,a3);
  
  if (x2 < cfmax){
    c31 = 0.0;
  }
  else{
    c31 = cfmax-x2;
  }
  
  return(c31);
}

double saltmat2::c32(double x1,double x2,double /*x3*/)
{
  double c32;
  double cfmax;

  CorD(16,kd,0,x1,cfmax,a2,a3);
  
  if (x2 < cfmax){
    c32 = 0.0;
  }
  else{
    c32 = -x1;
  }
  
  return(c32);
}

double saltmat2::c33(double /*x1*/,double /*x2*/,double /*x3*/)
{
  double c33;

  c33 = 1.0;

  return(c33);
}


/* Function calculates all auxiliary data - optional */
void saltmat2::auxiliarydata (double /*x1*/,double /*x2*/,double /*x3*/)
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
double saltmat2::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
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
double saltmat2::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
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
double saltmat2::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
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
double saltmat2::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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
double saltmat2::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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
double saltmat2::get_transmission_flux_11(double bv,double x1,double /*x2*/,double /*x3*/,long bc,long /*ipp*/)
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

double saltmat2::get_othervalue(long compother,long /*ipp*/, double x1,double /*x2*/,double /*x3*/)
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
void saltmat2::print_othervalue_name(FILE *out,long compother)
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

double  saltmat2::pgws(double t)
{
        double pgw;

        pgw = exp(23.5771 - 4042.9/(t - 37.58));

        return (pgw);
}
 /** Permeabilita Vodni Pary
    @param t  - temperature K
    @param p  -  presure Pa

*/
double saltmat2::permeabilitavodnipary(double t, double  p)
{

       double Da, Dp,rv,pa;
       double mi;

        CorD(4,kd,0,0,mi,a2,a3);

       rv = 461.5;
       pa = 101325;



       Da = (2.306e-5 * pa)/(rv * t * p)*pow((t/294.15),1.81);
       Dp = Da/mi;
       return (Dp);
 }

void saltmat2::CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2)
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
                prad  = int (MatData[cislochar][0][0]);
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
	 default:{
			fprintf (stderr,"\n\n unknown definition of material parameter is required");
			fprintf (stderr,"\n in function CorD (file %s, line %d)\n",__FILE__,__LINE__);
			}
     }

}



void saltmat2::binding_izoterm_derivation(double x2, double & derbi)
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
