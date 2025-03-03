#include "globalt.h"
#include "saltmat4.h"
#include "kunzel.h"

kunmat *kun3; 

/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
***********************************************/
/* x1 (w) is the volumetric water content [m3/m3]
   x2  (Cf) is the concentration of free ions in water [kg/m3]


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

discmat::discmat ()
{
  mw = 0.01801528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8.31441; //universal gas constant J.mol-1.K-1

}

discmat::~discmat ()
{}

/**
   JK, 11.4.2019
*/
double discmat::give_x1 (long nn)
{
  long k;
  double x1;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval (Tp->time);}
  
  return x1;
}

/**
   JK, 11.4.2019
*/
double discmat::give_x2 (long nn)
{
  long k;
  double x2;
  
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval (Tp->time);}
  
  return x2;
}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void discmat::matcond (matrix &d,long ri,long ci,long ipp)
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
  case 4:{
    matcond4d (d,ri,ci,ipp);//4D
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
void discmat::matcond2 (matrix &d,long ri,long ci,long ipp)
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
  case 4:{
    //matcond4d (d,ri,ci,ipp);//4D
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
void discmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
 
    
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,ipp);

  d[0][0] = k;
}

void discmat::matcond2d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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
void discmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];

  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,ipp);

 

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

void discmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,ipp);
  
  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}


void discmat::matcond4d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double x1,x2;
  k = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = k11 (x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    k = k12 (x1,x2,ipp);
  
  if((ri == 1) && (ci == 0))
    k = k21 (x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (x1,x2,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0; d[0][3]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0; d[1][3]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;   d[2][3]=0.0;
  d[3][0]=0.0; d[3][1]=0.0; d[3][2]=0.0; d[3][3]=k;
}

/**
   function creates capacity matrix of the material

   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void discmat::matcap (double &c,long ri,long ci,long ipp)
{
  double x1,x2;
  c = 0.0;
  
  x1 = Tm->ip[ipp].av[0];
  x2 = Tm->ip[ipp].av[1];
  
    
  if((ri == 0) && (ci == 0))
    c = c11 (x1,x2,ipp);
  if((ri == 0) && (ci == 1))
    c = c12 (x1,x2,ipp);
  
  if((ri == 1) && (ci == 0))
    c = c21 (x1,x2,ipp);
  if((ri == 1) && (ci == 1))
    c = c22 (x1,x2,ipp);
}


void discmat::values_correction (vector &/*nv*/,long /*ipp*/)
{
  // values correction
  
}



/**
   function reads data and material parameters

   @param in  - input file
*/
void discmat::read (XFILE *in)
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
				  if(i == 6)
				  {
				  xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
				  }
				  else
				  {
				  xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
				  }
			  }
		else
		{
			if(i == 6)
				  {
				  xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
				  }
				  else
				  {
				  xfscanf(in, "%lf %lf",&MatData[i][0][0],&MatData[i][1][0]);
				  }
		}
		switch (int(MatData[i][1][0])){
			case 2: for (j=1; j<=MatData[i][0][0];j++)
					{
						if(i ==6)
						{
							xfscanf(in, "%lf %lf",&MatData[i][1][j],&MatData[i][0][j] );
						}
						else
						{
							xfscanf(in, "%lf %lf",&MatData[i][0][j],&MatData[i][1][j] );
						}
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
	  }
  }
}


void discmat::print(FILE *out)
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



/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

double discmat::k11(double x1,double x2, long ipp)
{
  double k11;
  double kapak, ps, delta, dfdw;
  double a,n;
  int kod;
  double smc;
  CorD(7,kod,1,0,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

  kapak=Tm->ip[ipp].eqother[0];
 
  ps = pgws(x1,x2);
  delta = permeabilitavodnipary(x1,x2);
  
  dfdw = Tm->ip[ipp].eqother[2];
 
  k11 = 1000*kapak + ps*delta*dfdw;
  return (k11);
}

double discmat::k12(double x1,double x2, long ipp)
{
 double k12, delta,dpvsdt;
  double fi;

  delta = permeabilitavodnipary(x1,x2);
  fi = Tm->ip[ipp].eqother[1];
  dpvsdt =  derivation_saturation_water_vapour_pressure_temperature(x1,x2);
  
  k12 = delta*fi*dpvsdt;

  return (k12);
}


double discmat::k21(double x1,double x2, long ipp)
{
  double k21;
  double smc, lv1, ps, dp,dfdw , a,n;
  int kod;
  
  CorD(7,kod,1,0,x1,smc,a,n);

  if( x1 > smc) x1 = smc;
  
  
  lv1 = latent_heat_of_evaporation_of_water (x1,x2);
  ps = pgws (x1,x2);
  dp = permeabilitavodnipary(x1,x2);
   dfdw = Tm->ip[ipp].eqother[2];
     
  k21= lv1*ps*dp*dfdw;
  
  return (k21);
}

double discmat::k22(double x1,double x2, long /*ipp*/)
{
   double k22;
   double lv1, lambda, dpvsdt, dp, a, n;
  double smc;
  int kod;

  CorD(7,kod,1,0,x1,smc,a,n);

  if( x1 > smc) x1 = smc;
  lv1 = latent_heat_of_evaporation_of_water (x1,x2);
  dp = permeabilitavodnipary(x1,x2);
  CorD(10,kd,1,0,x1,lambda,a2,a3);//lambda
  dpvsdt =  derivation_saturation_water_vapour_pressure_temperature(x1,x2);

  k22= lambda + lv1 * dp * dpvsdt;

  return (k22);
}

double discmat::c11(double x1,double x2, long ipp)
{
  double c11;
  double ps;
  double fi,dfdw, wsoli, por;
  double a, n;
  int kod;
  double smc;
  CorD(7,kod,1,0,x1,smc,a,n);

  if( x1 > smc) x1 = smc;


  ps = pgws(x1,x2);
  CorD(3,kd,1,0,x1,por,a2,a3);
   
  fi = Tm->ip[ipp].eqother[1];
  dfdw = Tm->ip[ipp].eqother[2];
  CorD(17,kd,0,0,x1,wsoli,a2,a3);

  c11 = 1000 + mw/(gasr*294.15)*ps*(por-x1)*dfdw-mw*ps*fi/(gasr*294.15);

  return(c11);
}

double discmat::c12 (double /*x1*/, double /*x2*/, long /*ipp*/)
{
  return 0;
}

double discmat::c21 (double /*x1*/, double /*x2*/, long /*ipp*/)
{
  return 0;
}

double discmat::c22(double x1,double /*x2*/, long /*ipp*/)
{
  double c22;
  double rom, c,a2,a3, smc;
  int kod;
  
  CorD(7,kod,1,0,x1,smc,a2,a3);

  if( x1 > smc) x1 = smc;

  if (x1 < 0) x1 = 0;

  CorD(2,kd,0,0,0,rom,a2,a3);
  CorD(9,kd,smc,0,x1,c,a2,a3);
  
  c22 = rom * c;

  return(c22);
}

/* Function calculates all auxiliary data - optional */
void discmat::auxiliarydata (double /*x1*/,double /*x2*/)
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
double discmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,x1,x2;
  new_trc = 0.0;
  
  x1 = give_x1 (nn);
  x2 = give_x2 (nn);

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
double discmat::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval,x1,x2;
  new_nodval = 0.0;
  
  x1 = give_x1 (nn);
  x2 = give_x2 (nn);

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
double discmat::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  double flux,x1,x2;
  flux = 0.0;
  
  x1 = give_x1 (nn);
  x2 = give_x2 (nn);

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
double discmat::get_transmission_nodval_11(double bv,double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double discmat::get_transmission_transcoeff_11(double /*x1*/,double /*x2*/,long bc,long /*ipp*/)
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
double discmat::get_transmission_flux_11(double bv,double x1,double /*x2*/,long bc,long /*ipp*/)
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

double discmat::get_othervalue(long compother,long ipp, double x1,double /*x2*/)
{
  double other, a;

  
  switch (compother){
  case 1:{
	    other = -10000;
		 }
	  break;
  case 2:get_rel_hum2(x1,other,a);
	  break;
  case 3:other  = Tm->ip[ipp].eqother[10];
	  break;
  case 4:other  = Tm->ip[ipp].eqother[11];
	  break;
  case 5: get_rel_hum2(x1,other,a); // Tm->ip[ipp].eqother[5];
	  break;
  case 6:other  = 0.0; //Tm->ip[ipp].eqother[5];
	  break;
  case 7: other = 0.0;
	  break;
  case 0:{//first unknown

	  other = 0.0;
      break; 
  }
  default:{
	  print_err("unknown type of component is required in function ",__FILE__,__LINE__,__func__);
		  }
  }
  return (other);
}

/**
     function prints names of all variables in nodes
     @param out       - output file
     @param compother - number of other components
*/
void discmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//first unknown
    fprintf (out,"First unknown (Units)        ");//should be cha4nged
    break;
  }
  default:{
    print_err("unknown type of component is required in function ",__FILE__,__LINE__,__func__);
  }
  }
}
/** tlak nasycenych vodnich par
    @param x2  - temperature K
*/

double  discmat::pgws(double /*x1*/, double x2)
{
    double pgw;

	pgw = exp(23.5771 - 4042.9/(x2 - 37.58));

	return (pgw);
}
 /** Permeabilita Vodni Pary 
    @param x2  - temperature K
    @param p  -  presure Pa

*/
double discmat::permeabilitavodnipary(double x1, double x2)
{

	double Da, Dp,rv,pa, mi, p;
	
	CorD(4,kd,1,0,x1,mi,a2,a3);

	rv = 461.5;
	pa = 101325;
	p = 101325;

	Da = (2.306e-5 * pa)/(rv * x2 * p)*pow((x2/294.15),1.81);
	Dp = Da/mi;
	return (Dp);

 }

void discmat::CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2)
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
	 default:{
		 print_err("unknown type of component is required in function ",__FILE__,__LINE__,__func__);
			 }
  }

}



void discmat::kapa_values (int kod, double x1, double /*xpv*/, double /*ineq1*/, double &kapa)
{
  double kapa1, kapa2, kapa3;
  double smc,a,n;
  CorD(7,kod,1,0,x1,smc,a,n);

  if( x1 > smc) x1 = smc;

  switch (kod){
  case 0: kapa = 0.0;
	  break;
  case 1: CorD(5,kod,1,0,x1,kapa,kapa1,kapa2);
	  break;
  case 2: CorD(5,kod,1,0,x1,kapa,kapa1,kapa2);
	  break;
  case 30:{
	  CorD(5,kod,1,0,x1,kapa1,kapa2,kapa3);
	  kapa = kapa1*exp(kapa2*x1);
	  break;
		  }
   }

}


/**
   function selects auxiliary values
   
   JK, 7.1.2008
*/
void discmat::give_values (long ipp,double *av, double *pv, double *eq)
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
  eq[7] = Tm->ip[ipp].eqother[7];

}

/**
   function computes auxiliary values which are necessary for future computation
   
   JM, 29.5.2007
*/
void discmat::aux_values (double *in,double *inp, double *ineq,double *out)
{
  double x1,x2,x1p, x2p,  ypv1,ypv2,ypv3,ypv4,ypv5,ypv6,ypv7,ypv0;
 double jm1,fi,dfdw;
 int kod;
  
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
  ypv6 =ineq[6];
  ypv7 =ineq[7];
  
  double smc,a,n;
  CorD(7,kod,1,0,x1,smc,a,n);
	
	kapa_values(MatChar[5],x1,x1p,ypv0, jm1); 
	out[0]=jm1;
  
	sorption_izoterms_values(MatChar[6],x1,x1p,ypv1,fi,dfdw); 
	
	out[1]=fi;
	out[2]=dfdw;
	out[3]=smc;
    out[5]=-10000;
    out[6]=-10000;
	out[7]=-10000.0;


}
/**
   function saves auxiliary values
   
   JK, 7.1.2008
*/
void discmat::save_values (long ipp,double *out)
{
  Tm->ip[ipp].eqother[0]=out[0];
  Tm->ip[ipp].eqother[1]=out[1];
  Tm->ip[ipp].eqother[2]=out[2];
  Tm->ip[ipp].eqother[3]=out[3];
  Tm->ip[ipp].eqother[4]=out[4];
  Tm->ip[ipp].eqother[5]=out[5];
  Tm->ip[ipp].eqother[6]=out[6];
  Tm->ip[ipp].eqother[7]=out[7];

}

/**
   function defines inital values of quantities at integration points
   
   @param ipp - integration point pointer
   @param ido - index in array eq_other   

   JM, 29.5.2007
*/
void discmat::initvalues (long ipp,long /*ido*/)
{
	double x1,x2;

	x1 = Tm->ip[ipp].av[0];
	x2 = Tm->ip[ipp].av[1];
}

double discmat::derivation_saturation_water_vapour_pressure_temperature(double /*x1*/, double x2)
{
  return ((4042.9 * exp(23.5771-4042.9/(x2 - 37.58)))/((x2-37.58)*(x2-37.58)));   //TODO: Add your source code here
}

double discmat::latent_heat_of_evaporation_of_water(double /*x1*/, double x2)
{
        return (2.5008e6)*pow((273.15/x2),(0.167+x2*3.67e-4));
}

void discmat::sorption_izotherm_derivation(double x1, double /*x2*/, double & derfi)
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
    
}

void discmat::sorption_izoterms_values(int kod, double x1,double /*xpv*/, double /*ineq1*/, double & fi, double & dfdw)
{
 /*
 kod ..... 0- jen vlhkost, 1 - jen derivaci...
 */
 int s;
 double w_hyg,rh_hyg, smc, hvezda;
// double vyh1, vyh2, vyh3;
 double x2;

 x2 = 0.0;
       
       if (x1 < 0) x1 = 0;

       switch (kod)
          {
	   case 2: {
		   CorD(6,s,1,0,x1,fi,a2,a3);
		   dfdw = derivation_dy_dx(6,x1,0,1);
		   if(fi>0.976)
					{
						 
						CorD(7,s,0,0,x1,smc,a2,a3);
						CorD(6,s,0,0,x1,w_hyg,rh_hyg,w_hyg);
						hvezda =sortpion_isotherm_root_shifted(x1,w_hyg ,0.976);
						give_data_si_root_fi(x1,x2,w_hyg, smc, rh_hyg,hvezda,fi);
						give_data_si_root_dfidw(x1,x2,rh_hyg,smc,w_hyg,hvezda,dfdw);
					}
					else
					{
						dfdw = derivation_dy_dx(6,x1,0,1);
					}
			   }
            break;
          case 30:  //CorD(6,s,1,x1,a1,hmrh,a3);
                      CorD(6,s,0,0,x1,w_hyg,rh_hyg,a3);
                      hvezda =sortpion_isotherm_root_shifted(x1,w_hyg ,rh_hyg);
                      CorD(7,s,0,0,x1,smc,a2,a3);
					  give_data_si_root_fi(x1,x2,w_hyg, smc, rh_hyg,hvezda,fi);
					  give_data_si_root_dfidw(x1,x2,rh_hyg,smc,w_hyg,hvezda,dfdw);
                    
          break;
          case 31: 
			  CorD(6,s,1,0,x1,a1,a2,a3); 
			  fi = exp(a2*(1.0-pow((a1/x1),(a3))));
			  //fi  = a1*pow((1-(log(x1))/a2),(-1/a3));
			  dfdw = exp(a2*(1.0-pow((a1/x1),a3)))*a2*a3*pow(a1,a3)*pow(x1,(-1.0-a3));
                       //dfdw = (a1/(a2*x1*a3))*pow((1-(log(x1))/a2),(-1-(1/a3)));
		      break ; 
		      }
}
 

void discmat::give_data_si_root_fi(double x1, double /*x2*/, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh)
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

double discmat::sortpion_isotherm_root_shifted(double /*x1*/, double w_hyg, double rh_hyg)
{
  return(w_hyg/(1-sqrt(1-rh_hyg)));
}

void discmat::give_data_si_root_dfidw(double x1, double /*x2*/, double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw)
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

//double discmat::get_moisture(double rh)

void discmat::get_moisture (long /*nid*/,double in,double */*inp*/, double &out)
{
  int s;
 double moist,hmrh, mhmc, smc, u, a, n, rh;
 double t;

  t = Tp->time;
        rh = in;
       if (rh < 0) rh = 0;

       s = MatChar[6];      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
          case 2:   CorD(6,s,1,1,rh,moist,a1,a3);
		  break;
          case 30:  CorD(6,s,1,0,rh,a1,hmrh,a3);
                    CorD(7,s,0,0,rh,smc,a2,a3);
                    CorD(6,s,0,0,rh,mhmc,hmrh,a3);
                    //kun3->give_data(rh,mhmc, smc, hmrh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,0,0,0,u,a,n);
	    //moist = kun3->si_kk_hansen(rh, 0,u,a,n);
          break ;
	   }
	   out = moist;

    //TODO: Add your source code here
} 


double discmat::get_moisture2 ( double rh)
{
  int s;
 double moist,hmrh, mhmc, smc, u, a, n;


     if (rh < 0) rh = 0;

       s = MatChar[6];      // zjisteni zpusobu zadani SI ... 2-data, 30-ROOT, 31-HANSEN

       switch (s)
          {
          case 2:   CorD(6,s,1,1,rh,moist,a1,a3);
		  break;
          case 30:  CorD(6,s,1,0,rh,a1,hmrh,a3);
                    CorD(7,s,0,0,rh,smc,a2,a3);
                    CorD(6,s,0,0,rh,mhmc,hmrh,a3);
                    //kun3->give_data(rh,mhmc, smc, hmrh,moist);
                    //kunmat->give_data(fi,mhmc, smc, hmrh,moist);
          break;
          case 31:  CorD(6,s,0,0,0,u,a,n);
	    //moist = kun3->si_kk_hansen(rh, 0,u,a,n);
          break ;
	   }
	   return( moist);

    //TODO: Add your source code here
} 

double discmat::derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2)

{
	int i = 1;
	double outvalue;
	
	while (prom > MatData[matchar][pomk1][i]){
			i++;
		}
	outvalue = (MatData [matchar][pomk2][i] - MatData [matchar][pomk2][i-1])/(MatData [matchar][pomk1][i] - MatData [matchar][pomk1][i-1]);
	return(outvalue);

}



void discmat::water_content_relhum (long /*nid*/,double *in,double *inp,double *ineq, double *out)
//void discmat::water_content_relhum (long nid,double *in,double *inp,double *ineq, double out)
{
  double x1,x1p, ypv1;
 double fi,dfdw;
  
	x1 = in[0];
	x1p = inp[0];
	ypv1 =ineq[11];

	sorption_izoterms_values(MatChar[6],x1,x1p,ypv1,fi,dfdw); 
	out[10] = fi;

}
void discmat::get_rel_hum2 (double w,double &fi, double &dfdw)
{
	//double fi, dfdw;
  	sorption_izoterms_values(MatChar[6],w,0.0,0.0,fi,dfdw); 
	//return (fi);
}

void discmat::aux_values_elements (double *in,double *inp, double *ineq,double *out)
{
  double x1,x1p,ypv8;
 double fi,dfdw;
  
	x1 = in[0];
  	x1p = inp[0];

	ypv8 =ineq[8];
  
	sorption_izoterms_values(MatChar[6],x1,x1p,ypv8,fi,dfdw); 
	out[8] = fi;
	out[9] = x1;

}
