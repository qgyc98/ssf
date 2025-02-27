/*
    File:             threemediac.cpp
    Author:           Tomas Krejci
    Purpose:          computes conductivity and capacity matrices in a material point for coupled three media trasport with mechanics
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "threemediac.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "globalc.h"
#include "multiphasec.h"
#include "glasgowmatc.h"
#include "consol_hawf3.h"

medc3::medc3()
{}
medc3::~medc3()
{}

/**
   function computes conductivity matrix D in a material point for three media transfer
   for upper block
   
   @param d   - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc3::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.matcond_u(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[i].matcond_u(d,ri,ci,ipp);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].matcond_u(d,ri,ci,ipp);
      break;
    } 
      
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  case mech_moisture_salt:{
    fillm(0.0, d);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
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
void medc3::matcap_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.matcap_u(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[i].matcap_u(d,ri,ci,ipp);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].matcap_u(d,ri,ci,ipp);
      break;
    }

    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  case mech_moisture_salt:{
    fillm(0.0, d);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
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
void medc3::matcond_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cml->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cml->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.matcond_l(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[i].matcond_l(d,ri,ci,ipp);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].matcond_l(d,ri,ci,ipp);
      break;
    }
      
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  case mech_moisture_salt:{
    fillm(0.0, d);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
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
void medc3::matcap_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cml->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cml->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.matcap_l(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[i].matcap_l(d,ri,ci,ipp);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].matcap_l(d,ri,ci,ipp);
      break;
    }
      
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  case mech_moisture_salt:{
    fillm(0.0, d);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }    
}

/**
   function creates first part right-hand side matrix of the general material
   
   @param d   - right-hand side  %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc3::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.rhs_u1(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[Tm->ip[ipp].idm].fa1(t,pg,rhov);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].rhs_u1(d,ri,ci,ipp);//temporarilly commented
      break;
    }
      
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
}

/**
   function creates volume (second) part right-hand side matrix of the general material
   
   @param d   - right-hand side  %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void medc3::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cp->mednam){//names of transported media
  case mech_heat_moisture:{
    
    switch (Cmu->ip[ipp].tm){//material type
    case concreteBc:
    case baroghelBc:
    case C60baroghelBc:
    case C30baroghelBc:
    case o30bazantBc:
    case C60bazantBc:{
      multiphc mtphc;
      mtphc.rhs_u2(d,ri,ci,ipp);
      break;
    }
    case glasgowc:{
      //Cmu->tenchc[Tm->ip[ipp].idm].fa1(t,pg,rhov);
      break;
    }
    case consolhawf3c:{
      Cmu->consol_hawf3c[i].rhs_u2(d,ri,ci,ipp);
      break;
    }
     
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
    }
    }  
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown media name is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
}
