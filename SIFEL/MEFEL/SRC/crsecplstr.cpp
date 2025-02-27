#include "crsecplstr.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
crsecplstr::crsecplstr (void)
{
  t=0.0;  rho=0.0;  m=0.0;
}



/**
  Destructor is defined only for the formal purposes. 

  Created by JK,
*/
crsecplstr::~crsecplstr (void)
{
  
}



/**
  Function reads parameters of cross-section of twodimensional plane problem from the
  opened text file.
   
  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void crsecplstr::read (XFILE *in)
{
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case eigen_dynamics:{
    xfscanf (in,"%k%lf %k%lf %k%lf","thickness",&t,"rho",&rho,"m",&m);
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%k%lf %k%lf %k%lf","thickness",&t,"rho",&rho,"m",&m);
    break;
  }
  case mat_nonlinear_statics:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case geom_nonlinear_statics:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case earth_pressure:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case mech_timedependent_prob:
  case growing_mech_structure:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case layered_linear_statics:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case lin_floating_subdomain:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case nonlin_floating_subdomain:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case hemivar_inequal:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  case load_balancing:{
    xfscanf (in,"%k%lf","thickness",&t);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function prints parameters of cross-section of twodimensional plane problem into the
  opened text file.
   
  @param out - pointer to the opened text file

  @return The function does not return anything.

  TKr, 11/02/2013 according to read (XFILE *in)
*/
void crsecplstr::print (FILE *out)
{
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n%le",t);
    break;
  }
  case eigen_dynamics:{
    fprintf (out,"\n%le %le %le",t,rho,m);
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n%le %le %le",t,rho,m);
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n%le",t);
    break;
  }
  case geom_nonlinear_statics:{
    fprintf (out,"\n%le",t);
    break;
  }
  case earth_pressure:{
    fprintf (out,"\n%le",t);
    break;
  }
  case mech_timedependent_prob:
  case growing_mech_structure:{
    fprintf (out,"\n%le",t);
    break;
  }
  case layered_linear_statics:{
    fprintf (out,"\n%le",t);
    break;
  }
  case lin_floating_subdomain:{
    fprintf (out,"\n%le",t);
    break;
  }
  case nonlin_floating_subdomain:{
    fprintf (out,"\n%le",t);
    break;
  }
  case hemivar_inequal:{
    fprintf (out,"\n%le",t);
    break;
  }
  case load_balancing:{
    fprintf (out,"\n%le",t);
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
void crsecplstr::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      t=val[i];
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
