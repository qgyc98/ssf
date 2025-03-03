/*
  File:             tdisotrmat.cpp
  Author:           Jaroslav Kruis, 7.12.2017
  Purpose:          Calculates properties of general time dependent isotropic material for onemedium transfer
*/ 
#include "tdisotrmat.h"
#include "stochdrivert.h"
#include "globalt.h"

tdisotrmat::tdisotrmat (void)
{
  //  capacity
  c=0.0;
}



tdisotrmat::~tdisotrmat (void)
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
void tdisotrmat::matcond (matrix &d,long ri,long ci,long ipp)
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
void tdisotrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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
void tdisotrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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

void tdisotrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void tdisotrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();

  //double cf = Tm->ip[ipp].av[0];
  //cc = 1.0 + 9.8/(1.0+7250.0*cf)/(1.0+7250.0*cf);
  //cc = 1.0 + 4.9/(1.0+3335.0*cf)/(1.0+3335.0*cf);
  //cc=1.0;
}


/**
   function reads parameters
   
   @param in - input file
*/
void tdisotrmat::read (XFILE *in)
{
  //  c - capacity
  xfscanf (in,"%lf",&c);
  //  k - conductivity
  k.read (in);
}


/**
   function prints parameters
   
   @param out - outut file
*/
void tdisotrmat::print (FILE *out)
{
  fprintf (out,"  %e\n",c);
  k.print (out);
}

/**
   function creates conductivity coefficient of the isotropic material

   @retval k - heat conductivity %matrix of the isotropic material
*/

double tdisotrmat::get_k()
{
  return k.getval(Tp->time);
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double tdisotrmat::get_c()
{
  return(c);
}



/**  
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integartion point.

  Created by Tomas Koudelka, 13.1.2016
*/
double tdisotrmat::give_vol_moist(long ipp)
{
  if (Tp->mednam == moisture)
    return  Tm->ip[ipp].av[0];
  else
    print_err("invalid moisture quantity is required in tdisotrmat", __FILE__, __LINE__, __func__);

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
void tdisotrmat::give_dof_names (namevart *dofname, long ntm)
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
   function changes parameters of conductivity and capacity from a table
   @ param
*/
void tdisotrmat::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      c=val[0];
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
void tdisotrmat::give_reqntq(long */*antq*/)
{
}
