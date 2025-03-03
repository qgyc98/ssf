#include "crsection1d.h"
#include "globalt.h"
#include "stochdrivert.h"

crsection1d::crsection1d (void)
{
  a=0.0;  rho=0.0;
}
crsection1d::~crsection1d (void)
{
  
}

void crsection1d::read (XFILE *in)
{
  switch (Tp->tprob){
  case stationary_problem:{
    if (Tp->homogt!=0)//in case of homogenization
      xfscanf (in,"%lf",&rho);
    xfscanf (in,"%lf",&a);
    break;
  }
  case nonlinear_stationary_problem:{
    if (Tp->homogt!=0)//in case of homogenization
      xfscanf (in,"%lf",&rho);
    xfscanf (in,"%lf",&a);
    break;
  }
  case nonstationary_problem:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  case nonlinear_nonstationary_problem:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  case discont_nonstat_problem:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  case discont_nonlin_nonstat_problem:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  case growing_np_problem:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  case growing_np_problem_nonlin:{
    xfscanf (in,"%lf %lf",&rho,&a);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void crsection1d::print (FILE *out)
{
  switch (Tp->tprob){
  case stationary_problem:{
    fprintf (out," %lf\n",a);
    break;
  }
  case nonstationary_problem:{
    fprintf (out," %lf %lf\n",rho,a);
    break;
  }
  case nonlinear_nonstationary_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    fprintf (out," %lf %lf",rho,a);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void crsection1d::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      a=val[0];
      break;
    }
    case 1:{
      rho=val[1];
      break;
    }
    default:{
      print_err("wrong number of atribute is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}
