#include "sourcet.h"
#include "globalt.h"
#include "globmatt.h"

sourcet::sourcet()
{
  //  type of model of source
  sourtype = (sourcetype) 0;
  
  //  general function
  gf=NULL;
  
  //  constant value for stationary problems only
  v = 0.0;
  
  //  heat caused by cement hydration described by a special function
  hydrh=NULL;

  //nmultpv = 0;
  //cemh = NULL;//zatim??

  seebh=NULL;
}

sourcet::~sourcet()
{
  delete gf;
  delete hydrh;
  //  delete cemh;
  delete seebh;
}

/**
   function reads input data
   
   @param in - pointer to input file
   
   JK, revised 25.6.2005
*/
long sourcet::read (XFILE *in)
{
  //char *files = new char[1000];
  
  switch (Tp->tprob){
  case stationary_problem:{
    xfscanf (in, "%k%lf", "src_val", &v);
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    //  type of model of source
    //  matfunction=1
    //  cement_hydration=5
    xfscanf (in,"%k%m","src_type",  &sourcetype_kwdset, &sourtype);
    
    switch (sourtype){
    case matfunction:{
      gf = new gfunct;
      gf->read (in);
      break;
    }
    case concrete_heat:{
      hydrh = new hydrationheat;
      hydrh->read (in);
      break;
    }
    case cement_hydration:{
      //zatim??!!
      //cemh = new cemhyd (NULL);
      //getnexttxt (in);
      //inputln (in,files,1000);
      //cemh->read (files);
      break;
    }
    case seebeck:{
      seebh = new seebeckheat;
      seebh->read (in);
      break;
    }
    default:{
      print_err("unknown type of model of quantity source is required",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(0);
}

/**
   function prints output data
   
   @param out - pointer to output file
   
   TKr, 3.1.2006
*/
long sourcet::print(FILE *out)
{
  switch (Tp->tprob){
  case stationary_problem:{
    fprintf (out,"\n %lf ",v);
    break;
  }
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    //  type of model of source
    fprintf (out,"\n%d \n",sourtype);
    
    switch (sourtype){
    case matfunction:{
      gf->print (out);
      break;
    }
    case concrete_heat:{
      hydrh->print (out);
      break;
    }
    case cement_hydration:{
      //cemh->print (out);
      break;
    }
    case seebeck:{
      //
    break;
    }

    default:{
      print_err("unknown type of model of quantity source is required",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(0);
}



/**
   function returns value of source of quantity at required time
   
   @param eid - element id
   
   JK, revised 25.6.2005
*/
double sourcet::giveval (long eid)
{
  long i,ipp;
  double ret,t;
  
  switch (Tp->tprob){
  case stationary_problem:{
    ret=v;
    break;
  }
    
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:{
    
    switch (sourtype){
    case matfunction:{
      t = Tp->time;
      ret = gf->getval (t);
      break;
    }
    case concrete_heat:{
      t = Tp->time;
      ret = hydrh->give_value (t);
      break;
    }
    case cement_hydration:{
      ipp=Tt->elements[eid].ipp[0][0];
      i = Tm->ip[ipp].idm;
      ret = Tm->cemhydr[i].get_source(ipp,0,0);
      break;
    }
    case seebeck:{
      // heat conduction set as a source (seebeck effect)
      ipp=Tt->elements[eid].ipp[0][0];
      ret = seebh->give_value ();
      
      if(Tm->ip[ipp].ncompgrad == 1)
	ret = ret*Tm->ip[ipp].grad[0][0]*Tm->ip[ipp].grad[0][0];

      if(Tm->ip[ipp].ncompgrad == 2)
	ret = ret*(Tm->ip[ipp].grad[0][0]*Tm->ip[ipp].grad[0][0] 
		   + Tm->ip[ipp].grad[0][1]*Tm->ip[ipp].grad[0][1]);

      if(Tm->ip[ipp].ncompgrad == 3)
	ret = ret*(Tm->ip[ipp].grad[0][0]*Tm->ip[ipp].grad[0][0] 
		   + Tm->ip[ipp].grad[0][1]*Tm->ip[ipp].grad[0][1]
		   + Tm->ip[ipp].grad[0][2]*Tm->ip[ipp].grad[0][2]);

      break;
    }
    default:{
      print_err("unknown type of model of quantity source is required",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(ret);
}




long sourcet::read_prop(FILE */*in*/, long /*lc*/)
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



/**
  The function compares content of the actual object with the 
  content of src.

  @param src - compared object

  Returns:
   @retval 0 - in the case that objects are identical
   @retval 1 - in the case that objects differs

  Created by TKo 09.2010
*/
long sourcet::compare(sourcet &src)
{
  if (sourtype != src.sourtype)
    return 1;

  if (v != src.v)
    return 1;

  if (gf->compare(*src.gf))
    return 1;

  if (hydrh->compare(*src.hydrh))
    return 1;

//  if (cemh->compare(*src.hydrh))
//    return 1;
  return 0;
}
