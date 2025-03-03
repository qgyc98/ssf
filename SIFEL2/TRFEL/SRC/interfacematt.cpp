/*
  File:             interfacematt.cpp
  Author:           Jaroslav Kruis, 5. 1. 2016
  Purpose:          Calculates properties of interface edge for onemedium transfer
*/ 
#include "interfacematt.h"
#include "stochdrivert.h"
#include "globalt.h"

interfacematt::interfacematt (void)
{
  //  coefficient of conductivity in the direction of element
  ks=0.0;
  //  coefficient of conductivity in the normal direction
  kn=0.0;
  //  fictitious width of the element
  h=0.0;
}



interfacematt::~interfacematt (void)
{
}



/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   Created by JK, 2016
*/
void interfacematt::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  default:{
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   Created by JK, 2016
*/
void interfacematt::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  fillm(0.0,d);
  
  d[0][0] = ks;   d[0][1] = 0.0;
  d[1][0] = 0.0;  d[1][1] = kn;
}



/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   Created by JK, 2016
*/
void interfacematt::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
}


/**
   function reads parameters
   
   @param in - input file

   Created by JK, 2016
*/
void interfacematt::read (XFILE *in)
{
  //  ks - conductivity in the element direction
  //  kn - conductivity in the normal direction
  //  h - fictitious width of the element
  //  h and kn are correlated and they should be selected properly
  xfscanf (in,"%lf %lf %lf",&ks,&kn,&h);
}


/**
   function prints parameters
   
   @param out - outut file

   Created by JK, 2016
*/
void interfacematt::print (FILE *out)
{
  fprintf (out,"  %e %e %e",ks,kn,h);
}

/**
   function returns conductivity coefficient in the direction of the element

   @retval ks - conductivity coefficient in the direction of the element

   Created by JK, 2016
*/
double interfacematt::get_ks ()
{
  return ks;
}

/**
   function returns conductivity coefficient in the normal direction to the element

   @retval kn - conductivity coefficient in the normal direction to the element

   Created by JK, 2016
*/
double interfacematt::get_kn ()
{
  return kn;
}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Koudelka, 13.1.2016
*/
void interfacematt::give_dof_names (namevart *dofname, long ntm)
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
  Function returns volumetric moisture content in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of volumetric moisture content stored in the integartion point.

  Created by Tomas Koudelka, 13.1.2016
*/
double interfacematt::give_vol_moist(long ipp)
{
  if (Tp->mednam == moisture)
    return  Tm->ip[ipp].av[0];
  else
  return 0.0;
  //    print_err("invalid moisture quantity is required in interfacemat", __FILE__, __LINE__, __func__);

  return 0.0;
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
  Created by Tomas Koudelka, 13.1.2016
*/
void interfacematt::give_reqntq(long */*antq*/)
{
}


