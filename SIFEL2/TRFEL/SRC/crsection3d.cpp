#include "crsection3d.h"
#include "globalt.h"
#include "stochdrivert.h"

crsection3d::crsection3d (void)
{
  rho=0.0;
}
crsection3d::~crsection3d (void)
{
  
}

/**
   function reads variables of cross section
   
   @param in - input file
   
*/
void crsection3d::read (XFILE *in)
{
  switch (Tp->tprob){
  case stationary_problem:{
    xfscanf (in,"%lf",&rho);
    break;
  }
  case nonlinear_stationary_problem:{
    break;
  }
  case nonstationary_problem:{
    xfscanf (in,"%lf",&rho);
    break;
  }
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    xfscanf (in,"%lf",&rho);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function prints variables of cross section
   
   @param out - output file
   
*/
void crsection3d::print (FILE *out)
{
  switch (Tp->tprob){
  case stationary_problem:{
    break;
  }
  case nonstationary_problem:{
    fprintf (out," %lf\n",rho);
    break;
  }
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    fprintf (out," %lf\n",rho);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void crsection3d::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      rho=val[1];
      break;
    }
    default:{
      print_err("wrong number of atribute is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}
