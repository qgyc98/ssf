#include "crseclayer.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"



/**
  Constructor initializes data members to zero or default values.
*/
crseclayer::crseclayer (void)
{
  nl = 0;
  th = 0.0; rho = 0.0; m = 0.0;
  layth = NULL; layz  = NULL;
}



/**
  Destructor is defined only for the formal purposes. 
*/
crseclayer::~crseclayer (void)
{
  delete [] layth;
  delete [] layz;
}



/**
  Function reads parameters of layered cross-section from the opened 
  text file and eventually evaluates z-coordinates and thickness.
   
  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JF 10/2014,
*/
void crseclayer::read (XFILE *in)
{
  long i;
  
  switch (Mp->tprob){
  case linear_statics:
  case mat_nonlinear_statics:{
    xfscanf (in,"%ld",&nl);
    layth = new double[nl];
    layz  = new double[nl];
      
    for (i=0;i<nl;i++){
      xfscanf (in,"%le",&layth[i]);
    }
    break;
  }
  case eigen_dynamics:{
    
    break;
  }
  case forced_dynamics:{
    break;
  }
  case geom_nonlinear_statics:{
    break;
  }
  case earth_pressure:{
    break;
  }
  case mech_timedependent_prob:
  case growing_mech_structure:{
    break;
  }
  case layered_linear_statics:{
    break;
  }
  case lin_floating_subdomain:{
    break;
  }
  case nonlin_floating_subdomain:{
    break;
  }
  case hemivar_inequal:{
    break;
  }
  case load_balancing:{
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
  zcoordinates();
}

/**
  Function computes thickness and z-coordintes of the layers.

  @return The function does not return anything.

  Created by JF 10/2014,
*/
void crseclayer::zcoordinates ()
{
  long i;
  double p=0.0;
  
  for (i=0;i<nl;i++){
        th += layth[i];
  }
  
  for (i=0;i<nl;i++){
        layz[i] = th/2-p-layth[i]/2;
        p += layth[i];
  }
}


/**
  Function changes class parameters.
  The function is used in stochastic or fuzzy computations.
   
  @param atcs - selected cross-section parameters (parameters which are changed)
  @param val - %vector containing new values of parameters
   
  @return The function does not return anything.

*/
void crseclayer::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 1:{    
      rho=val[i];
      break;
    }
    default:{
      print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
    }
  }
}
