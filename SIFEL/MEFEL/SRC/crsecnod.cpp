#include "crsecnod.h"
#include "global.h"
#include "probdesc.h"


crsecnod::crsecnod (void)
{
  t=0.0;  m=0.0;
}
crsecnod::~crsecnod (void)
{
  
}

void crsecnod::read (XFILE *in)
{
  switch (Mp->tprob){
  case layered_linear_statics:{
    xfscanf (in,"%lf",&t);
    break;
  }
  case nonlinear_dynamics:{
    xfscanf (in,"%lf",&m);
    break;
  }
  default:
    print_err("unknown problem type %d is required.", __FILE__, __LINE__, __func__, Mp->tprob);
  }
}


void crsecnod::print (FILE *out)
{
  switch (Mp->tprob){
  case layered_linear_statics:{
    fprintf (out,"\n %le",t);
    break;
  }
  case nonlinear_dynamics:{
    fprintf (out,"\n %le",m);
    break;
  }    
  default:
    print_err("unknown problem type %d is required.", __FILE__, __LINE__, __func__, Mp->tprob);
  }
}

