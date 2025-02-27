/*
    File:             twomediac.cpp
    Author:           Tomas Krejci
    Purpose:          computes conductivity and capacity matrices in a material point for coupled two media trasport with mechanics
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "twomediac.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "globalc.h"
#include "multiphasec.h"
#include "consol_awf2.h"
#include "consol_wf2.h"

medc2::medc2()
{}
medc2::~medc2()
{}


/**
   function computes conductivity matrix D in a material point for three media transfer
   for upper block
   
   @param d   - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc2::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolwf2c:{
      Cmu->consol_wf2c[i].matcond_u(d,ri,ci,ipp);
      break;
    } 
    case consolawf2c:{
      Cmu->consol_awf2c[i].matcond_u(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].matcond_u(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }    
}


/**
   function computes capacity matrix D in a material point for three media transfer
   for upper block
   
   @param d   - capacity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc2::matcap_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolwf2c:{
      Cmu->consol_wf2c[i].matcap_u(d,ri,ci,ipp);
      break;
    }       
    case consolawf2c:{
      Cmu->consol_awf2c[i].matcap_u(d,ri,ci,ipp);
      break;
    }       
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].matcap_u(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }    
}



/**
   function computes conductivity matrix D in a material point for three media transfer
   for lower block
   
   @param d   - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc2::matcond_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cml->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cml->ip[ipp].tm){//material type
    case consolwf2c:{
      Cml->consol_wf2c[i].matcond_l(d,ri,ci,ipp);
      break;
    } 
    case consolawf2c:{
      Cml->consol_awf2c[i].matcond_l(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].matcond_l(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }    
}


/**
   function computes capacity matrix D in a material point for three media transfer
   for lower block
   
   @param d   - capacity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc2::matcap_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cml->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cml->ip[ipp].tm){//material type
    case consolwf2c:{
      Cml->consol_wf2c[i].matcap_l(d,ri,ci,ipp);
      break;
    }       
    case consolawf2c:{
      Cml->consol_awf2c[i].matcap_l(d,ri,ci,ipp);
      break;
    }             
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].matcap_l(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }    
}


/**
   function creates first part right-hand side matrix of the general material - influence of gravity accel.
   
   @param d   - right-hand side  %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void medc2::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolwf2c:{
      Cmu->consol_wf2c[i].rhs_u1(d,ri,ci,ipp);
      break;
    }      
    case consolawf2c:{
      Cmu->consol_awf2c[i].rhs_u1(d,ri,ci,ipp);
      break;
    }      
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].rhs_u1(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  } 
}


/**
   function creates first part right-hand side matrix of the general material - influence of initial values (temperature)
   
   @param d   - right-hand side  %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void medc2::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_air_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolwf2c:{
      //Cmu->consol_wf2c[i].rhs_u2(d,ri,ci,ipp);//no influence
      break;
    }      
    case consolawf2c:{
      //Cmu->consol_awf2c[i].rhs_u2(d,ri,ci,ipp);//no influence
      break;
    }      
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  case mech_heat_water:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case consolhwf2c:{
      Cmu->consol_hwf2c[i].rhs_u2(d,ri,ci,ipp);
      break;
    } 
    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
      abort();
    }
    }  
    break;
  }
  default:{
    print_err("unknown media name is required",__FILE__,__LINE__,__func__);
    abort();
  }
  } 
}
