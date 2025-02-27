#include "pvalt.h"
#include "globalt.h"

pvalt::pvalt()
{
  v = 0.0;
  ipv = 0.0;
  func[0] = '\0';
  eq = NULL;
  var = NULL;
  nmultpv = 0;
}

pvalt::~pvalt()
{
  delete eq;
}

/**
   function reads data from input file
   
   @param in - input file
*/
long pvalt::read(XFILE *in)
{
  Parser par;
  
  //  initial prescribed value
  xfscanf (in,"%lf",&ipv);
  
  switch (Tp->tprob){
  case stationary_problem:
  case nonlinear_stationary_problem:{
    xfscanf (in,"%lf",&v);
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:{
    
    //  type of function
    xfscanf(in,"%d",(int*)&tfunc);
    
    switch (tfunc){
    case constant:{
      xfscanf (in,"%lf",&v);
      break;
    }
    case pars:{
      xfscanf(in, "%255s", func);
      //      fgets(func, 255, in->file);
      eq = par.TextToTree(func);
      if (eq == NULL)
	{
	  print_err("Error parsing expression",__FILE__,__LINE__,__func__);
	  return(1);
	}
      var = eq->Variables.at(0);
      break;
    }
    case tab_file:{
      tabf.read_data_file(in);
      break;
    }
    case tab_file2:{
      tabf.read_data_file2(in);
      break;
    }
    case tab:{
      tabf.read(in);
      break;
    }
    case multpv:{//prescribed values multiplied by time function
      xfscanf (in,"%ld",&nmultpv);
      xfscanf (in,"%lf",&v);
      break;
    }
    default:{
      print_err("unknown function is required",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  default:{
    print_err("unknown problem type",__FILE__,__LINE__,__func__);
  }
  }
  
  return(0);
}

long pvalt::print(FILE *out)
{
  Parser par;
  
  fprintf (out,"\n\n %lf",ipv);
  
  switch (Tp->tprob){
  case stationary_problem:
  case nonlinear_stationary_problem:{
    fprintf (out," %lf",v);
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:{
    
    fprintf(out,"  %d",tfunc);
    
    switch (tfunc){
    case constant:{
      fprintf (out," %lf",v);
      break;
    }
    case pars:{//dodelat??
      //xfscanf(in, "%255s", func);
      ////fgets(func, 255, in);
      //eq = par.TextToTree(func);
      //if (eq == NULL)
      //{
      //  fprintf(stderr, "\nError parsing expression in function pvalt::read (file %s, line %d)\n",__FILE__,__LINE__);
      //  return(1);
      //}
      //var = eq->Variables.at(0);
      break;
    }
    case tab_file:{
      tabf.print_data_file(out);
      break;
    }
    case tab_file2:{
      tabf.print_data_file2(out);
      break;
    }
    case tab:{
      tabf.print(out);
      break;
    }
    case multpv:{//prescribed values multiplied by time function
      fprintf (out,"\n %ld",nmultpv);
      fprintf (out,"\n %lf",v);
      break;
    }
    default:{
    }
    }
    
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }

  return(0);
}

/**
   @param time - actual time
   
   JK, 15.4.2019
*/
double pvalt::getval (double time)
{
  double ret;
  
  switch (Tp->tprob){
  case stationary_problem:
  case nonlinear_stationary_problem:{
    ret=v;
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    switch (tfunc){
    case constant:{
      ret=v;
      break;
    }
    case pars:{
      if (var)
	var->Value = time;
      ret = eq->Evaluate();
      break;
    }
    case tab_file2:
    case tab_file:
    case tab:{
      ret=tabf.getval(time);
      break;
    }
    case multpv:{
      //ret=v*Tb->lc[0].give_fact (t);//old
      ret=ipv + (v - ipv)*Tb->lc[0].give_fact (time,nmultpv);//new
      break;
    }
    default:{
      print_err("unknown time function is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(ret);
}


/**
   @param t - previous time
   
   JK, 15. 4. 2019
*/
double pvalt::getprevval (double t)
{
  double ret;
  
  switch (Tp->tprob){
  case stationary_problem:
  case nonlinear_stationary_problem:{
    ret=v;
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    //t=Tp->time-Tp->timecont.backwarddt;
    switch (tfunc){
    case constant:{
      ret=v;
      break;
    }
    case pars:{
      if (var)
	var->Value = t;
      ret = eq->Evaluate();
      break;
    }
    case tab_file2:
    case tab_file:
    case tab:{
      ret=tabf.getval(t);
      break;
    }
    case multpv:{
      //t=Tp->time-Tp->timecont.backwarddt;
      //ret=v*Tb->lc[0].give_fact (t);//old
      ret=ipv + (v - ipv)*Tb->lc[0].give_fact (t,nmultpv);//new
      break;
    }
    default:{
      print_err("unknown time function is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(ret);
}


long pvalt::read_prop(FILE */*in*/, long /*lc*/)
{
  /*
  fscanf(in, "%ld", &nlc);
  if ((nlc < 1) || (nlc > lc))
    return(2);
  fscanf(in, "%le", &f);
  xfscanf(in, "%255s", func);
  //fgets(func, 255, in);

  */
  return(0);

}

double pvalt::getipv(void)
{
  return (ipv);
}
