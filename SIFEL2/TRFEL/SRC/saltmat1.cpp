#include "globalt.h"
#include "saltmat1.h"
#include "globmatt.h"

dampermeability saltmat1::damper;

saltmat1::saltmat1 ()
{
  //  permeability influenced by damage - default option is off
  daminfl = off;
}

saltmat1::~saltmat1 ()
{}

/**
   function reads material characteristics
   
   @param in - input file

   JM, 18. 11. 2013
*/
void saltmat1::read (XFILE *in)
{
  //  moisture diffusivity
  kappa.read (in);
  //  saturated moisture
  sm.read (in);
  //  Dcoef
  dcoef.read (in);
  //  binding isotherm
  bindiso.read (in);
}

/**
   function prints material characteristics
   
   @param out - output file

   JM, 18. 11.2013
*/
void saltmat1::print (FILE *out)
{
  //  moisture diffusivity
  kappa.print (out);
  //  saturated moisture
  sm.print (out);
  //  Dcoef
  dcoef.print (out);
  //  binding isotherm
  bindiso.print (out);
}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 30. 5. 2014
*/
void saltmat1::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of space dimension is required",__FILE__,__LINE__,__func__);
  }
  }

  if (daminfl == on){
    //  conductivity matrix is modified due to damage
    damper.matcond (d,ipp);
  }
  Tm->ip[ipp].eqother[0] = d[0][0];
}


/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - integration point id 
*/
void saltmat1::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,cf;
  k = 0.0;

  w = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = k11 (w);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 1) && (ci == 0))
    k = k21 (w,cf);
  if((ri == 1) && (ci == 1))
    k = k22 (w,cf);
  
  d[0][0] = k;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void saltmat1::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,cf;
  k = 0.0;

  w = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = k11 (w);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 1) && (ci == 0))
    k = k21 (w,cf);
  if((ri == 1) && (ci == 1))
    k = k22 (w,cf);
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - integration point id
*/

void saltmat1::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,cf;
  k = 0.0;

  w = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = k11 (w);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 1) && (ci == 0))
    k = k21 (w,cf);
  if((ri == 1) && (ci == 1))
    k = k22 (w,cf);
  
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
   @param ipp - integration point id
*/
void saltmat1::matcap (double &c,long ri,long ci,long ipp)
{
  double w,cf;
  
  if((ri == 0) && (ci == 0))
    c = c11 ();
  if((ri == 0) && (ci == 1))
    c = c12 ();
  if((ri == 1) && (ci == 0)){
    cf = Tm->ip[ipp].av[1];
    c = c21 (cf);
  }
  if((ri == 1) && (ci == 1)){
    w = Tm->ip[ipp].av[0];
    cf = Tm->ip[ipp].av[1];
    c = c22 (w,cf);
  }
}


/******************************************************
Conductivity and capacity terms for C and K matrices
*******************************************************/

double saltmat1::c11 ()
{
  double c11;
  c11 = 1.0;
  return c11;
}

double saltmat1::c12 ()
{
  double c12;
  c12 = 0.0;
  return c12;
}

double saltmat1::c21 (double cf)
{
  double c21;
  c21 = cf;
  return c21;
}

double saltmat1::c22 (double w,double cf)
{
  double c22,dCbDcf;
  
  if (cf <= 0.0)
    cf = 1.0e-10;
  if (w < 0.0)
    w = 0.0;
  
  dCbDcf = bindiso.derivative_isotherm_value (cf);
  
  c22 = w + dCbDcf;
  
  return c22;
}



double saltmat1::k11 (double w)
{
  double k11,wsat;
  
  if (w < 0.0)
    w = 0.0;

  wsat = sm.getval (w);
  if (w > wsat)
    w = wsat;
  
  k11 = kappa.getval (w);

  return k11;
}

double saltmat1::k12 ()
{
  double k12;
  
  k12 = 0.0;

  return k12;
}

double saltmat1::k21 (double w,double cf)
{
  double k21,wsat;

  if (w < 0.0)
    w = 0.0;
  wsat = sm.getval (w);
  if (w > wsat)
    w = wsat;
  
  k21 = w*dcoef.getval (cf);

  return k21;
}

double saltmat1::k22 (double w,double cf)
{
  double k22,wsat;
  
  if (w < 0.0)
    w = 0.0;
  wsat = sm.getval (w);
  if (w > wsat)
    w = wsat;


  k22 = cf*kappa.getval (w);
  
  return k22;
}

/******************
Boundary conditions
*******************/
	

/**
   function determines transmission coefficient
   in some cases, the boundary conditions are prescribed in different
   variables than variables used in the problem
   for example, moisture content is used in the problem but boundary
   condition is prescribed with the help of pressures

   @param trc - prescribed transmission coefficient on the boundary
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat1::transmission_transcoeff (double trc,long ri,long ci,long nid,long /*bc*/)
{
  double new_trc,w,cf;
  new_trc = 0.0;
  
  //  moisture content
  w = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  
  if((ri == 0) && (ci == 0))
    new_trc = trc;
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  
  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
    new_trc = trc;
  
  return new_trc;
}

/**
   function creates correct transfer coefficient on the boundary (transmission) for 1st medium
   @param f11        - correct transfer coefficient
   @param w ... cc  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
/*
double saltmat4::get_transmission_transcoeff_11 (double w,double cf,double cc,double t,long bc,long ipp)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    trc=1.0;//should be changed
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }

  return(trc);
}
*/


/**
   function determines nodal value
   in some cases, the boundary conditions are prescribed in different
   variables than variables used in the problem
   for example, moisture content is used in the problem but boundary
   condition is prescribed with the help of pressures

   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat1::transmission_nodval (double nodval,long ri,long ci,long nid,long /*bc*/)
{
  double new_nodval,w,cf;
  new_nodval = 0.0;
  
  //  moisture content
  w = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  
  if((ri == 0) && (ci == 0))
    //new_nodval = get_transmission_nodval_11(nodval,w,cf,cc,t,bc,ipp);
    new_nodval = nodval;
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  
  
  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
    new_nodval = nodval;
  
  return new_nodval;
}


/**
   function creates correct new nodal value on the boundary (transmission) for 1st medium
   @param new_nodval - new prescribed value near the boundary
   @param bv         - value of prescribed value near the boundary
   @param w ... cc  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
/*
double saltmat4::get_transmission_nodval_11 (double bv,double w,double cf,double cc,double t,long bc,long ipp)
{
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    nodval = bv;//should be changed
    break;
  } 
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(nodval);
}
*/

/**
   function computes flux through boundary
   
   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat1::transmission_flux (double nodval,long ri,long ci,long nid,long bc)
{
  double flux,w,cf;
  flux = 0.0;
  
  //  moisture content
  w = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  
  
  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_ww (nodval,w,bc);
  if((ri == 0) && (ci == 1))
    flux = 0.0;
  
  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
    flux = 0.0;
  
  return (flux);
}


/**
   function computes flux through the boundary (transmission - convective mass transfer) for the first medium
   
   @param bv - prescribed value on the boundary
   @param w - actual moisture content on the boundary
   @param bc - type of boundary condition
*/
double saltmat1::get_transmission_flux_ww (double bv,double w,long bc)
{
  double flux;
  
  switch (bc){//type of prescribed variable
  case 30:{//transmission - boundary flux
    flux = (bv - w);//should be changed
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
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

double saltmat1::get_othervalue(long compother,long /*ipp*/, double x1,double /*x2*/,double /*x3*/)
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
void saltmat1::print_othervalue_name(FILE *out,long compother)
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


/**
   The function returns ordered dof names of primary unknowns 
   required by the model.
   
   @param dofname   - array of uknown name for particular nodal dofs (output)
                      dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
   @param ntm       - number of transported media = number of nodal dof = length of array dofname
   
   JK 29. 5. 2014
*/
void saltmat1::give_dof_names (namevart *dofname, long ntm)
{
  if (ntm < 1){
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_moisture;
  dofname[1] = trf_salt_conc;
}
