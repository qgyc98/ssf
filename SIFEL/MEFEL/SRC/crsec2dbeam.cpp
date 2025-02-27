#include "crsec2dbeam.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
crsec2dbeam::crsec2dbeam (void)
{
  a=0.0;  iy=0.0;  shearcoeff=0.0;  rho=0.0;  t=0.0;
}



/**
  Destructor releases allocated memory of the crsec2dbeam object.

  Created by JK,
*/
crsec2dbeam::~crsec2dbeam (void)
{
  
}



/**
  Function reads parameters of cross-section of twodimensional beams from the 
  opened text file.
   
  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void crsec2dbeam::read (XFILE *in)
{
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  case eigen_dynamics:{
    xfscanf (in,"%lf %lf %lf %lf",&a,&iy,&shearcoeff,&rho);
    break;
  }
  case linear_stability:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  case mech_timedependent_prob:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  case mat_nonlinear_statics:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%lf %lf %lf %lf",&a,&iy,&shearcoeff,&rho);
    break;
  }
  case earth_pressure:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  case layered_linear_statics:{
    xfscanf (in,"%lf %lf %lf %lf",&a,&iy,&shearcoeff,&t);
    break;
  }
  case growing_mech_structure:{
    xfscanf (in,"%lf %lf %lf",&a,&iy,&shearcoeff);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
  Function prints parameters of cross-section of twodimensional beams into the 
  opened text file.
   
  @param out - pointer to the opened text file

  @return The function does not return anything.

  TKr, 11/02/2013, according to read (XFILE *in)
*/
void crsec2dbeam::print (FILE *out)
{
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n %le %le %le",a,iy,shearcoeff);
    break;
  }
  case mech_timedependent_prob:{
    fprintf (out,"\n %le %le %le",a,iy,shearcoeff);
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n %le %le %le",a,iy,shearcoeff);
    break;
  }
  case eigen_dynamics:{
    fprintf (out,"\n %le %le %le %le",a,iy,shearcoeff,rho);
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n %le %le %le %le",a,iy,shearcoeff,rho);
    break;
  }
  case earth_pressure:{
    fprintf (out,"\n %le %le %le",a,iy,shearcoeff);
    break;
  }
  case layered_linear_statics:{
    fprintf (out,"\n %le %le %le %le",a,iy,shearcoeff,t);
    break;
  }
  case growing_mech_structure:{
    fprintf (out,"\n %le %le %le",a,iy,shearcoeff);
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function returns area of the cross-section.
 
  @return The function returns corss-section area.

  Created by JK,
*/
double crsec2dbeam::give_area ()
{
  return a;
}



/**
  The function returns moments of inertia of the cross-section.

  @param i - array for storage of moments of inertia (output)
 
  @return The function returns moments of inertia in the parameter i.

  Created by JK,
*/
void crsec2dbeam::give_moments (double *i)
{
  i[0]=iy;
}



/**
  The function returns shear coefficients of the cross-section.

  @param sc - array for storage of shear coefficients (output)
 
  @return The function returns shear coefficients in the parameter sc.

  Created by JK,
*/
void crsec2dbeam::give_shearcoeff (double *sc)
{
  sc[0]=shearcoeff;
}



/**
  The function returns density of the cross-section.
 
  @return The function returns corss-section density.

  Created by JK,
*/
double crsec2dbeam::give_density ()
{
  return rho;
}



/**
  Function changes class parameters.
  The function is used in stochastic or fuzzy computations.
   
  @param atcs - selected cross-section parameters (parameters which are changed)
  @param val - %vector containing new values of parameters
   
  @return The function does not return anything.

  Created by JK,
*/
void crsec2dbeam::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      a=val[i];
      break;
    }
    case 1:{
      iy=val[i];
      break;
    }
    case 2:{
      shearcoeff=val[i];
      break;
    }
    case 3:{
      rho=val[i];
      break;
    }
    default:
      print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
  }
}
