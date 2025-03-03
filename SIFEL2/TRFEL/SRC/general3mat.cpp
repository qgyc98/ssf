#include "globalt.h"
#include "general3mat.h"

/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
***********************************************/

general3mat::general3mat ()
{}

general3mat::~general3mat ()
{}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void general3mat::matcond (matrix &d,long ri,long ci,long ipp)
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
void general3mat::matcond1d (matrix &d,long ri,long ci,long ipp)
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
void general3mat::matcond2d (matrix &d,long ri,long ci,long ipp)
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

void general3mat::matcond3d (matrix &d,long ri,long ci,long ipp)
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
void general3mat::matcap (double &c,long ri,long ci,long ipp)
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
void general3mat::read (FILE *in)
{
  //fscanf (in,"%lf %lf",&par1,&par2);
}


/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

double general3mat::k11(double x1,double x2,double x3)
{
  double k11;

  //here will be something??!!

  return (k11);
}

double general3mat::k12(double x1,double x2,double x3)
{
  double k12;

  //here will be something??!!

  return (k12);
}

double general3mat::k13(double x1,double x2,double x3)
{
  double k13;

  //here will be something??!!

  return (k13);
}

double general3mat::k21(double x1,double x2,double x3)
{
  double k21;

  //here will be something??!!

  return (k21);
}

double general3mat::k22(double x1,double x2,double x3)
{
  double k22;

  //here will be something??!!

  return (k22);
}

double general3mat::k23(double x1,double x2,double x3)
{
  double k23;

  //here will be something??!!

  return (k23);
}

double general3mat::k31(double x1,double x2,double x3)
{
  double k31;

  //here will be something??!!

  return (k31);
}

double general3mat::k32(double x1,double x2,double x3)
{
  double k32;

  //here will be something??!!

  return (k32);
}

double general3mat::k33(double x1,double x2,double x3)
{
  double k33;

  //here will be something??!!

  return (k33);
}

double general3mat::c11(double x1,double x2,double x3)
{
  double c11;

  //here will be something??!!

  return(c11);
}

double general3mat::c12(double x1,double x2,double x3)
{
  double c12;

  //here will be something??!!

  return(c12);
}

double general3mat::c13(double x1,double x2,double x3)
{
  double c13;

  //here will be something??!!

  return(c13);
}

double general3mat::c21(double x1,double x2,double x3)
{
  double c21;

  //here will be something??!!

  return(c21);
}

double general3mat::c22(double x1,double x2,double x3)
{
  double c22;

  //here will be something??!!

  return(c22);
}

double general3mat::c23(double x1,double x2,double x3)
{
  double c23;

  //here will be something??!!

  return(c23);
}

double general3mat::c31(double x1,double x2,double x3)
{
  double c31;

  //here will be something??!!

  return(c31);
}

double general3mat::c32(double x1,double x2,double x3)
{
  double c32;

  //here will be something??!!

  return(c32);
}

double general3mat::c33(double x1,double x2,double x3)
{
  double c33;

  //here will be something??!!

  return(c33);
}


/* Function calculates all auxiliary data - optional */
void general3mat::auxiliarydata (double x1,double x2,double x3)
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
double general3mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_trc,x1,x2,x3;
  new_trc = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval();}
  
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
double general3mat::transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_nodval,x1,x2,x3;
  new_nodval = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval();}
  
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
double general3mat::transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double flux,x1,x2,x3;
  flux = 0.0;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {x1 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x1 = 0.0;}
  if (k<0)   {x1 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,1);
  if (k>0)   {x2 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x2 = 0.0;}
  if (k<0)   {x2 = Tb->lc[0].pv[0-k-1].getval();}
  k=Gtt->give_dof(nn,2);
  if (k>0)   {x3 = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {x3 = 0.0;}
  if (k<0)   {x3 = Tb->lc[0].pv[0-k-1].getval();}
  
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
double general3mat::get_transmission_nodval_11(double bv,double x1,double x2,double x3,long bc,long ipp)
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
double general3mat::get_transmission_transcoeff_11(double x1,double x2,double x3,long bc,long ipp)
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
double general3mat::get_transmission_flux_11(double bv,double x1,double x2,double x3,long bc,long ipp)
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

double general3mat::get_othervalue(long compother,long ipp, double x1,double x2,double x3)
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
void general3mat::print_othervalue_name(FILE *out,long compother)
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
