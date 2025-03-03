#include <math.h>
#include "isoviscous.h"
#include "global.h"
#include "matrix.h"
#include "vector.h"
#include "intpoints.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
isoviscous::isoviscous (void)
{
  //  coefficient of volumetric viscosity
  xi=0.0;
  //  coefficient of dynamic viscosity
  eta=0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK,
*/
isoviscous::~isoviscous (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by JK, 20.11.2017
*/
void isoviscous::read (XFILE *in)
{
  //  coefficient of volumetric viscosity, coefficient of dynamic viscosity
  xfscanf (in,"%k %lf %k %lf", "xi", &xi, "eta", &eta);
}



/**
  The function prints model parameters in to the opened text file.

  @param out - pointer to the opened text file for output

  @return The function does not return anything.

  Created by JK, 20.11.2017
*/
void isoviscous::print (FILE *out)
{
  fprintf (out, "%le %le", xi,eta);
}


/**
  Function assembles damping %matrix of material.
   
  @param d - damping %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material damping %matrix in the parameter d.

  Created by JK, 20.11.2017
*/
void isoviscous::matdamp (matrix &d,strastrestate ssst)
{
  switch (ssst){
    /*
  case bar:{
    matstiff_bar (d);
    break;
  }
  case plbeam:{
    matstiff_plbeam (d);
    break;
  }
  case spacebeam:{
    matstiff_spacebeam (d);
    break;
  }
    */
  case planestress:{
    matdamp_plstress (d);
    break;
  }
    /*
  case planestrain:{
    matstiff_plstrain (d);
    break;
  }
  case platek:{
    matstiff_platek (d);
    break;
  }
  case plates:{
    matstiff_plates (d);
    break;
  }
  case axisymm:{
    matstiff_axi (d);
    break;
  }
    */
  case spacestress:{
    matdamp_spacestr (d);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
  Function creates damping matrix of the viscous
  isotropic material for 2D problems (plane stress).

  @param d - damping matrix of the material

  @return The function returns material damping %matrix in the parameter d.

  Created by JK, 20.11.2017
*/
void isoviscous::matdamp_plstress (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0] = xi+eta/3.0;      d[0][1] = xi-2.0/3.0*eta;  d[0][2] = 0.0;
  d[1][0] = xi-2.0/3.0*eta;  d[1][1] = xi+eta/3.0;      d[1][2] = 0.0;
  d[2][0] = 0.0;             d[2][1] = 0.0;             d[2][2] = eta;
}

/**
  Function creates damping %matrix of the viscous
  isotropic material for 3D problems.
   
  @param d - damping %matrix of the material

  Created by JK, 20.11.2017
*/
void isoviscous::matdamp_spacestr (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0] = xi+eta/3.0;      d[0][1] = xi-2.0/3.0*eta;  d[0][2] = xi-2.0/3.0*eta;
  d[1][0] = xi-2.0/3.0*eta;  d[1][1] = xi+eta/3.0;      d[1][2] = xi-2.0/3.0*eta;
  d[2][0] = xi-2.0/3.0*eta;  d[2][1] = xi-2.0/3.0*eta;  d[2][2] = xi+eta/3.0;

  d[3][3] = eta;             d[4][4]=eta;               d[5][5]=eta;
}
