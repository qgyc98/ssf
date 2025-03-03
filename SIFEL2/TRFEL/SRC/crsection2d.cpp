#include "crsection2d.h"
#include "globalt.h"
#include "stochdrivert.h"

crsection2d::crsection2d (void)
{
  t=0.0;  rho=0.0;
}
crsection2d::~crsection2d (void)
{
  
}

void crsection2d::read (XFILE *in)
{
  switch (Tp->tprob){
  case stationary_problem:{
    if (Tp->homogt!=0)//in case of homogenization
      xfscanf (in,"%lf",&rho);
    xfscanf (in,"%lf",&t);
    break;
  }
  case nonlinear_stationary_problem:{
    if (Tp->homogt!=0)//in case of homogenization
      xfscanf (in,"%lf",&rho);
    xfscanf (in,"%lf",&t);
    break;
  }
  case nonstationary_problem:{
    xfscanf (in,"%lf %lf",&rho,&t);
    break;
  }
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    xfscanf (in,"%lf %lf",&rho,&t);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void crsection2d::print (FILE *out)
{
  switch (Tp->tprob){
  case stationary_problem:{
    fprintf (out," %lf\n",t);
    break;
  }
  case nonstationary_problem:{
    fprintf (out," %lf %lf\n",rho,t);
    break;
  }
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    fprintf (out," %lf %lf\n",rho,t);
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void crsection2d::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      t=val[0];
      break;
    }
    case 1:{
      rho=val[1];
      break;
    }
    default:{
      print_err("wrong number of atributes",__FILE__,__LINE__,__func__);
    }
    }
  }
}
