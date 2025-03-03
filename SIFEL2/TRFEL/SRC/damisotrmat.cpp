/*
  File:             damisotrmat.cpp
  Author:           Tomas Koudelka, Jaroslav Kruis, 6.5.2014
  Purpose:          Calculates properties of general isotropic material for linear onemedium transfer
*/ 

#include "damisotrmat.h"
#include "stochdrivert.h"
#include "globalt.h"

dampermeability damisotrmat::damper;

damisotrmat::damisotrmat (void)
{
  //  thermal conductivity  
  k=0.0;
  //  capacity
  c=0.0;
  // permeability influenced by damage - default option is off
  daminfl = off;
}



damisotrmat::~damisotrmat (void)
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
void damisotrmat::matcond (matrix &d,long ri,long ci,long ipp)
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
   @param ipp - number of integration point
*/
void damisotrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void damisotrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double kk;

  kk = get_k();
  
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

void damisotrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
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
void damisotrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();
}


/**
   function reads parameters
   
   @param in - input file
*/
void damisotrmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf",&c,&k); //c=capacity, k=conductivity
  xfscanf (in,"%m", &answertype_kwdset, &daminfl);
}


/**
   function prints parameters
   
   @param out - outut file
*/
void damisotrmat::print (FILE *out)
{
  fprintf (out,"  %e %e",c,k);
}

/**
   function creates conductivity coefficient of the isotropic material

   @retval k - heat conductivity %matrix of the isotropic material
*/

double damisotrmat::get_k()
{
  return(k);
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double damisotrmat::get_c()
{
  return(c);
}



void damisotrmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//actual permeability
    fprintf (out,"Actual permeability");
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

  Created by Tomas Koudelka, 6.12.2013
*/
void damisotrmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_moisture;
}



/**
   function changes parameters of conductivity and capacity from a table
   @ param
*/
void damisotrmat::changeparam (atsel &atm,vector &val)
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
void damisotrmat::give_reqntq(long *antq)
{
  if (daminfl == on){
    //  damage parameter
    antq[scal_iso_damage-1] = 1;
    //  process zone length
    antq[proc_zone_length-1] = 1;
    //  crack width
    antq[crack_width-1] = 1;
  }
}
