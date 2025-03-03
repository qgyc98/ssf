#include "strretalg.h"
#include "global.h"
#include "alias.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
strretalg::strretalg ()
{
  zero=0.0;
  ni=0;
  err=0.0;
  err_stv_req=0.0;
  err_stv_max=0.0;
  h_min = 0.0;
  rkt = nork;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by JK,
*/
strretalg::~strretalg ()
{
}



/**
  The function reads parameters for stress-return algorithm from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by JK,
*/
void strretalg::read (XFILE *in)
{
  //  type of stress return algorithm
  //  0 - no stress return algorithm
  //  1 - cutting plane algorithm
  //  10 - general nonlinear stress return algorithm
  //  11 - Runge-Kutta-Fehlsberg algorithm
  xfscanf (in,"%m", &stressretalgtype_kwdset,(int*)&tsra);
  
  switch (tsra){
  case nostressretalg:{
    if (Mespr==1)
      fprintf (stdout,"\n no stress return algorithm is required");
    break;
  }
  case cp:{
    xfscanf (in, "%k%ld %k%lf", "ni", &ni, "err", &err);
    break;
  }
  case gsra:{
    xfscanf (in, "%k%ld %k%lf", "ni", &ni, "err", &err);
    break;
  }
  case rkfsra:{
    // xfscanf (in, "%k%m %k%ld %k%lf %k%lf %k%lf %k%lf", "rktype", &rktype_kwdset, &rkt, "ni", &ni, "err_sig", &err, "err_stv_req", &err_stv, "err_stv_max", &err_stv_max, "h_min", &h_min);
    xfscanf (in, "%k%m %k%ld %k%lf %k%lf", "rktype", &rktype_kwdset, &rkt, "ni", &ni, "err", &err, "h_min", &h_min);
    break;
  }
  default:
    print_err("unknown stress return algorithm is required", __FILE__, __LINE__, __func__);
  }
  if (err < 0.0)
    print_err("negative error threshold of stress integration algorithm is required\n (err=%le < 0)\n",
              __FILE__, __LINE__, __func__, err);
  if (ni < 0)
    print_err("maximum number of stress integration iterations is negative\n (ni=%ld < 0)\n",
              __FILE__, __LINE__, __func__, ni);
}

/**
  The function prints parameters for stress-return algorithm from the opened text file given
  by the parameter out.

  @param out - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void strretalg::print (FILE *out)
{
  //  type of stress return algorithm
  fprintf (out,"%d",(int)tsra);
  
  switch (tsra){
  case nostressretalg:{
    break;
  }
  case cp:{
    fprintf (out, " %ld %e", ni, err);
    break;
  }
  case gsra:{
    fprintf (out, " %ld %e", ni, err);
    break;
  }
  case rkfsra:{
    //fprintf (out, " %d %ld %e %e %e %e", int(rkt), ni, err, err_stv_req, err_stv_max, h_min);
    fprintf (out, " %d %ld %e %e", int(rkt), ni, err, h_min);
    break;
  }
  default:
    print_err("unknown stress return algorithm is required", __FILE__, __LINE__, __func__);
  }
  
}


/**
  The function returns required type of stress-return algorithm.

  @return The function returns required type of stress-return algorithm.

  Created by JK,
*/
stressretalgtype strretalg::give_tsra ()
{
  return tsra;
}



/**
  The function returns required maximum number of iterations of reuired stress-return algorithm.

  @return The function returns maximum number of iterations.

  Created by JK,
*/
long strretalg::give_ni ()
{
  return ni;
}



/**
  The function returns required error of stress-return algorithm.

  @return The function returns required error of stress-return algorithm.

  Created by JK,
*/
double strretalg::give_err ()
{
  return err;
}



/**
  The function returns type of Runge-Kutta stress-return algorithm.

  @return The function returns required type of Runge-Kutta stress-return algorithm.

  Created by TKo, 12.2015
*/
rktype strretalg::give_rkt()
{
  if (tsra == rkfsra)
    return rkt;

  print_err("unknown type of stress return algorithm is required", __FILE__, __LINE__, __func__);
  return nork;
}



/**
  The function returns required error for stresses of stress-return algorithm.

  @return The function returns required error of stress-return algorithm.

  Created by TKo, 08.2017
*/
double strretalg::give_err_sig ()
{
  return err;
}



/**
  The function returns required error and maximum error for state variables of stress-return algorithm.


  @param

  @return The function returns required error of stress-return algorithm.

  Created by TKo, 08.2017,
*/
void strretalg::give_err_stv (double &err_req, double &err_max)
{
  err_req = err_stv_req;
  err_max = err_stv_max;
}



/**
  The function returns minimum step length for Runge-Kutta-Fehlsberg stress-return algorithm.

  @return The function returns required error of stress-return algorithm.

  Created by TKo, 10.2015
*/
double strretalg::give_hmin()
{
  if (tsra == rkfsra)
    return h_min;

  print_err("unknown type of stress return algorithm is required", __FILE__, __LINE__, __func__);
  return 0;
}
