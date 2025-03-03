#include "multpvalt.h"
#include "globalt.h"

multpvalt::multpvalt()
{
  v = 0.0;
}

multpvalt::~multpvalt()
{
}

long multpvalt::read(FILE *in)
{
  if (Tp->tprob == nonstationary_problem || Tp->tprob == nonlinear_nonstationary_problem ||  Tp->tprob == growing_np_problem ||  Tp->tprob == growing_np_problem_nonlin)
{
    fscanf (in,"%lf",&v);
  }
  return(0);
}

double multpvalt::getval(void)
{
  double ret,t;
  
  switch (Tp->tprob){
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    t=Tp->time;
    ret=v*Tb->lc[0].give_fact (t);
    return(ret);
  }
  default:{
    fprintf (stderr,"\n\n");
  }
  }
  
  return ret;
}

long multpvalt::read_prop(FILE *in, long lc)
{
  /*
  fscanf(in, "%ld", &nlc);
  if ((nlc < 1) || (nlc > lc))
    return(2);
  fscanf(in, "%le", &f);
  //fgets(func, 255, in);
  xfscanf(in, "%255s", func);

  */
  return(0);

}
