#include "crsec2dbar.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
crsec2dbar::crsec2dbar (void)
{
  a=0.0;  rho=0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK,
*/
crsec2dbar::~crsec2dbar (void)
{
  
}



/**
  Function reads data of cross-sections from the opened text file.
   
  @param in - pointer to the opened XFILE
   
  @return The function does not return anything.

  Created by JK,
*/
void crsec2dbar::read (XFILE *in)
{
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case eigen_dynamics:{
    xfscanf (in,"%k%lf %k%lf %k%lf","a",&a,"rho",&rho,"m",&m);
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%k%lf %k%lf %k%lf","a",&a,"rho",&rho,"m",&m);
    break;
  }
  case mat_nonlinear_statics:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case mech_timedependent_prob:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case growing_mech_structure:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case lin_floating_subdomain:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case nonlin_floating_subdomain:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  case load_balancing:{
    xfscanf (in,"%k%lf","a",&a);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function prints data of cross-sections into the opened text file.
   
  @param out - pointer to the opened FILE
   
  @return The function does not return anything.

  TKr, 11/02/2013, according to read (XFILE *in)
*/
void crsec2dbar::print (FILE *out)
{
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n %le",a);
    break;
  }
  case eigen_dynamics:{
    fprintf (out,"\n %le %le",a,rho);
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n %le %le",a,rho);
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n %le",a);
    break;
  }
  case mech_timedependent_prob:{
    fprintf (out,"\n %le",a);
    break;
  }
  case growing_mech_structure:{
    fprintf (out,"\n %le",a);
    break;
  }
  case lin_floating_subdomain:{
    fprintf (out,"\n %le",a);
    break;
  }
  case nonlin_floating_subdomain:{
    fprintf (out,"\n %le",a);
    break;
  }
  case load_balancing:{
    fprintf (out,"\n %le",a);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function changes corss-section parameters. 
  It is used for stochastic calculations.

  @param atcs - selected cross-section parameters (parameters which are changed)
  @param val - %vector containing new values of parameters
   
  @return The function does not return anything.

  Created by JK,
*/
void crsec2dbar::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      a=val[i];
      break;
    }
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
