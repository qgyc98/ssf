/*
  File:             isotrmat.cpp
  Author:           Jaroslav Kruis, 20.12.2002
  Purpose:          Calculates properties of general isotropic material for linear onemedium transfer
*/ 
#include "isotrmat.h"
#include "stochdrivert.h"
#include "globalt.h"

isotrmat::isotrmat (void)
{
  //  thermal conductivity  
  k=0.0;
  //  capacity
  c=0.0;
  //  reaction coefficient
  r=0.0;
}



isotrmat::~isotrmat (void)
{
}



/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void isotrmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
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
void isotrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  //double cf = Tm->ip[ipp].av[0];
  //kk=kk/(1.0+4.9/(1.0+3335.0*cf)/(1.0+3335.0*cf));
  //if (cf<1.0e-9)
  //cf=1.0e-9;
  
  //kk=kk/(1.0+0.013096544*0.36*pow(cf,-0.64));
  
  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void isotrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  //double cf = Tm->ip[ipp].av[0];
  //kk=kk/(1.0+4.9/(1.0+3335.0*cf)/(1.0+3335.0*cf));
  //double cf = Tm->ip[ipp].av[0];
  //if (cf<1.0e-9)
  //cf=1.0e-9;
  
  //kk=kk/(1.0+1.03*0.36/0.08*pow(cf,-0.64));
  //kk=kk/(1.0+0.013096544*0.36*pow(cf,-0.64));


  fillm(0.0,d);
  
  d[0][0] = kk;   d[0][1] = 0.0;
  d[1][0] = 0.0;  d[1][1] = kk;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void isotrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
  fillm(0.0,d);
  
  d[0][0]=kk;   d[0][1]=0.0;  d[0][2]=0.0;
  d[1][0]=0.0;  d[1][1]=kk;   d[1][2]=0.0;
  d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=kk;
}


/**
   function creates capacity %matrix of the material

   @param cc   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void isotrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();

  //double cf = Tm->ip[ipp].av[0];
  //cc = 1.0 + 9.8/(1.0+7250.0*cf)/(1.0+7250.0*cf);
  //cc = 1.0 + 4.9/(1.0+3335.0*cf)/(1.0+3335.0*cf);
  //cc=1.0;
}

/**
   function creates reaction %matrix of the material

   @param rr  - reaction %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 12. 6. 2019
*/
void isotrmat::matreact (double &rr,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  rr = get_r ();
}

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
double isotrmat::transmission_transcoeff(double trc,long /*ri*/,long /*ci*/,long /*nn*/,long bc,long /*ipp*/)
{
  double new_trc=0.0;
  
  switch (bc){//type of prescribed variable
  case 5:{
    new_trc=1.0;
    break;
  }
  case 11:{
    new_trc=trc;
    break;
  }
  case 30:{//heat transmission
    new_trc=1.0;
    break;
  }
  case 31:{//heat transmission for testing
    new_trc=0.0;
    break;
  }
  case 50:{//heat transmission - boundary condition computed from the actual boundary flux
    new_trc=1.0;
    break;
  }
  case 51:{//heat transmission - boundary condition computed from the actual boundary flux
    new_trc=1.0;
    break;
  }
  case 90:{//radiation
    new_trc=1.0;
    break;
  }
  default:{
    print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  
  new_trc = new_trc*trc;

  return (new_trc);
}


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
double isotrmat::transmission_transcoeff(double trc,long /*ri*/,long /*ci*/,long /*nn*/,long bc,long /*ipp*/, int flag)
{
  double new_trc=0.0;
  
  switch (bc){//type of prescribed variable
  case 5:{
    new_trc=1.0;
    break;
  }
  case 11:{
    new_trc=trc;
    break;
  }
  case 30:{//heat transmission
    new_trc=1.0;
    break;
  }
  case 31:{//heat transmission for testing
    new_trc=0.0;
    break;
  }
  case 50:{//heat transmission - boundary condition computed from the actual boundary flux
    new_trc=1.0;
    break;
  }
  case 51:{//heat transmission => heat flux - this boundary condition is created provisionally for salt diffusion problem
    // the flux of salt ions is calculated for inflow flux from external and interal values of concentration 
    if(flag == 1)//for rhs vector
      new_trc=1.0;
    else
      new_trc=0.0;//for matrix
    break;
  }
  case 90:{//radiation
    if(flag == 1)//for rhs vector
      new_trc=1.0;
    else
      new_trc=0.0;//for matrix
    break;
  }
  default:{
    print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  
  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double isotrmat::transmission_nodval(double nodval,double trc2,long /*ri*/,long /*ci*/,long nn,long bc,long ipp)
{
  long k;
  double dt=0.0,new_nodval=0.0,t;
  //provisionally for salt diffusion:
  double por,rhoc,rhow,ecl,c0,kcl;


  k=Gtt->give_dof(nn,0);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  switch (bc){//type of prescribed variable
  case 5:{
    new_nodval = nodval;
    break;
  }
  case 11:{
    new_nodval = nodval;
    break;
  }
  case 30:{//heat transmission
    new_nodval = nodval;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    new_nodval = (nodval - t);
    break;
  }
  case 50:{//heat transmission - boundary condition computed from the actual boundary flux
    dt = Tp->timecont.actualbacktimeincr();
    new_nodval=Tm->ip[ipp].av[0] + Tb->fluxes[0]*nodval*dt;//cumulated value from the first integration point of the element
    //fprintf(Outt,"new_nodval = %e\n",new_nodval);
    //fflush(Outt);
    break;
  }
  case 51:{//heat transmission => heat flux - this boundary condition is created provisionally for salt diffusion problem
    // the flux of salt ions is calculated for inflow flux from external and interal values of concentration
    por = 0.122;
    rhoc = 2233.0;
    rhow = 998;
    //ecl = 1.16e-12;
    ecl = 1.0e-7;
    c0 = 0.51;
    //kcl = 1.74e-12;
    kcl = 1.5e-7;
    new_nodval = -0.035453*(ecl*(1/por*rhoc/rhow*Tm->ip[ipp].av[0] - nodval)+kcl*(1/por*rhoc/rhow*Tm->ip[ipp].av[0]/c0)*
			    (1/por*rhoc/rhow*Tm->ip[ipp].av[0]/c0)*exp(-1.15/por*rhoc/rhow*Tm->ip[ipp].av[0]));
    break;
  }
  case 90:{//radiation
    //new_nodval = (nodval - t) + trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
    new_nodval = trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
    break;
  }
  default:{
    print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  
  return (new_nodval);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double isotrmat::transmission_flux(double nodval,double trc2,long /*ri*/,long /*ci*/,long nn,long bc,long ipp)
{
  long k;
  double dt=0.0,flux=0.0,t;
  //provisionally for salt diffusion:
  double por,rhoc,rhow,ecl,c0,kcl;

  k=Gtt->give_dof(nn,0);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    flux = (nodval - t);
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    flux = (nodval - t);
    break;
  }
  case 50:{//heat transmission - boundary condition computed from the actual boundary flux
    dt = Tp->timecont.actualbacktimeincr();
    flux = (Tm->ip[ipp].av[0] + Tb->fluxes[0]*nodval*dt - t);//cumulated value from the first integration point of the element
    break;
  }
  case 51:{//heat transmission => heat flux - this boundary condition is created provisionally for salt diffusion problem
    // the flux of salt ions is calculated for inflow flux from external and interal values of concentration
    por = 0.122;
    rhoc = 2233.0;
    rhow = 998;
    //ecl = 1.16e-12;
    ecl = 1.0e-7;
    c0 = 0.51;
    //kcl = 1.74e-12;
    kcl = 1.5e-7;
    flux = -0.035453*(ecl*(1/por*rhoc/rhow*Tm->ip[ipp].av[0] - nodval)+kcl*(1/por*rhoc/rhow*Tm->ip[ipp].av[0]/c0)*
		      (1/por*rhoc/rhow*Tm->ip[ipp].av[0]/c0)*exp(-1.15/por*rhoc/rhow*Tm->ip[ipp].av[0]));
    break;
  }
  case 90:{//radiation
    //flux = (nodval - t) + trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
    flux = trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
    break;
  }
  default:{
    print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }

  return (flux);
}



/**
   function reads parameters
   
   @param in - input file
*/
void isotrmat::read (XFILE *in)
{
  //  c - capacity
  //  k - conductivity
  xfscanf (in,"%lf %lf",&c,&k);
  
  if (Tp->react>0){
    //  r - reaction coefficient
    xfscanf (in,"%lf",&r);
  }
  
}


/**
   function prints parameters
   
   @param out - outut file
*/
void isotrmat::print (FILE *out)
{
  fprintf (out,"  %e %e",c,k);
  if (Tp->react>0){
    fprintf (out,"  %e\n",r);
  }else{
    fprintf (out,"\n");
  }
}

/**
   function creates conductivity coefficient of the isotropic material

   @retval k - heat conductivity %matrix of the isotropic material
*/

double isotrmat::get_k()
{
  return(k);
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double isotrmat::get_c()
{
  return(c);
}

/**
   function returns reaction coefficnet of the isotropic material
   
   @retval r - reaction coefficient of the isotropic material
   
   JK, 12. 6. 2019
*/
double isotrmat::get_r()
{
  return (r);
}


/**  
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integartion point.

  Created by Tomas Koudelka, 13.1.2016
*/
double isotrmat::give_vol_moist(long ipp)
{
  if (Tp->mednam == moisture)
    return  Tm->ip[ipp].av[0];
  else
    print_err("invalid moisture quantity is required in isotrmat", __FILE__, __LINE__, __func__);

  return 0.0;
}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Koudelka, 6.12.2013
*/
void isotrmat::give_dof_names (namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  if (Tp->mednam == moisture)
    dofname[0] = trf_moisture;
  else
    dofname[0] = trf_temperature;
}


/**
   function returns temperature

   @param ipp - integration point number

   @retval t - temperature

   17/07/2018, TKr
*/
double isotrmat::give_temperature(long ipp)
{
  double t;

  t = Tm->ip[ipp].av[0];

  return(t);
}


/**
   function changes parameters of conductivity and capacity from a table
   @ param
*/
void isotrmat::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      c=val[0];
      break;
    }
    case 1:{
      k=val[1];
      break;
    }
    default:{
      print_err("wrong number of atribute is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   25. 4. 2014, by TKo
*/
void isotrmat::give_reqntq(long */*antq*/)
{
}


/**
   function computes other variables in nodes
   @param compother - number of other components
   @param ipp - rhrst integration point on element
   @param t - temperature on actual node
*/

double isotrmat::get_othervalue(long compother,double t, long /*ipp*/)
{
  double other;
  
  switch (compother){
  case 0:{//temperature in deg. C
    other = t-273.15;
    break;
  }
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
  return (other);
  
}


/**
     function prints names of all variables in nodes
     @param out - output rhle
     @param compother - number of other components
*/
void isotrmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//Temperature in deg. C
    fprintf (out,"Temperature (C)             ");
    break;
  }
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
}
