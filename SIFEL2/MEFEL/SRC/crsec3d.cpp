#include "crsec3d.h"
#include "iotools.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"
#include <stdio.h>



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
crsec3d::crsec3d (void)
{
  rho=0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK,
*/
crsec3d::~crsec3d (void)
{
  
}



/**
  Function reads cross-section parameters of threedimensional object from the 
  opened text file.
   
  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void crsec3d::read (XFILE *in)
{
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%k%lf","rho",&rho);
    break;
  }
  case eigen_dynamics:{
    xfscanf (in,"%k%lf","rho",&rho);
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%k%lf","rho",&rho);
    break;
  }
  case mat_nonlinear_statics:{
    xfscanf (in,"%k%lf","rho",&rho);
    break;
  }
  case growing_mech_structure:
  case mech_timedependent_prob:{
    xfscanf (in,"%k%lf","rho",&rho);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
  Function prints cross-section parameters of threedimensional object into the 
  opened text file.
   
  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by TKr, 02/03/2013
*/
void crsec3d::print (FILE *out)
{
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n %le",rho);
    break;
  }
  case eigen_dynamics:{
    fprintf (out,"\n %le",rho);
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n %le",rho);
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n %le",rho);
    break;
  }
  case growing_mech_structure:
  case mech_timedependent_prob:{
    fprintf (out,"\n %le",rho);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function changes class parameters.
  The function is used in stochastic or fuzzy computations.
   
  @param atcs - selected cross-section parameters (parameters which are changed)
  @param val - %vector containing new values of parameters
   
  @return The function does not return anything.

  Created by JK,
*/
void crsec3d::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      rho=val[i];
      break;
    }
    default:{
      print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
    }
  }
}
