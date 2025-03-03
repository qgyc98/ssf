#include "hardsoft.h"
#include "alias.h"
#include "global.h"
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>

hardsoft::hardsoft ()
{
  //  no hardening/softening
  ths = nohs;
  //  norm of plastic strain tensor
  ivhs = plstrnorm;
  
  //  limit for norm of plastic strain
  epspu=0.0;
  //  limit for norm of plastic strain in tension
  epsput=0.0;
  //  limit for norm of plastic strain in compression
  epspuc=0.0;
  
  //  computer zero
  zero=1.0e-8;
  
}

hardsoft::~hardsoft ()
{
}

/**
   function reads input data
   
   @param in - pointer to input file
   
   JK, 9.8.2005, modified 17.2.2007
*/
void hardsoft::read (XFILE *in)
{
  //  type of hardening/softening
  xfscanf (in,"%m", &hardensoften_kwdset, (int*)&ths);
  
  switch (ths){
  case nohs:{
    //  no hardening/softening
    if (Mespr==1)
      fprintf (stdout,"\n no hardening/softening is required");
    break;
  }
  case plstrainnorm:{
    //  hardening/softening is proportional to the norm of the plastic strain tensor
    if (Mespr==1)
      fprintf (stdout,"\n hardening/softening is proportional to the norm of the plastic strain tensor");
    break;
  }
  case limplstrainnorm:{
    //  hardening/softening is proportional to the norm of the plastic strain tensor, norm is limited by the parameter
    if (Mespr==1)
      fprintf (stdout,"\n hardening/softening is proportional to the limited norm of the plastic strain tensor");
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of hardening/softening is required in function read (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  if (ths!=nohs){
    //  input variable to hardening/softening
    xfscanf (in,"%m", &hsinputvar_kwdset, (int*)&ivhs);
    
    if (ivhs==limplstrnorm){
      xfscanf (in,"%lf",&epspu);
    }
    if (ivhs==limanplstrnorm){
      xfscanf (in,"%lf %lf",&epsput,&epspuc);
    }
  }
  
}



/**
   The function prints data about hardening to the opened text file.
   
   @param out - pointer to the output text file
   
   TKo 4.2016
*/
void hardsoft::print(FILE *out)
{
  //  type of hardening/softening
  fprintf(out, "%d", int(ths));
  
  switch (ths){
    case nohs:
      //  no hardening/softening
      break;
    case plstrainnorm:
      break;
    case limplstrainnorm:
      break;
    default:
      print_err("unknown type %d of hardening/softening is required", __FILE__, __LINE__, __func__, int(ths));
  }
  
  if (ths != nohs){
    //  input variable to hardening/softening
    fprintf(out, " %d", int(ivhs));
    
    if (ivhs==limplstrnorm){
      fprintf(out, " %lf", epspu);
    }
    if (ivhs==limanplstrnorm){
      fprintf(out, " %lf %lf", epsput, epspuc);
    }
  }
}



/**
   function computes values of hardening/softening function
   
   @param sigt - stress components stored in 3x3 %matrix
   @param dgds - first derivatives of plastic potential with respect to stress components (stored in 3x3 %matrix)
   @param dgdsds - second derivatives of plastic potential with respect to stress components (stored in 6x6 %matrix)
   @param dhds - derivatives of hardening/softening function with respect to stress components (stored in %matrix 6 x ncomphard)
   
   JK, 17.2.2007
*/
void hardsoft::hvalues (vector &/*sigt*/,vector &/*dgds*/,vector &/*h*/)
{
  /*
  double norm,inv;
  
  //  cleaning output vector
  fillv (0.0,h);
  
  switch (ths){
  case nohs:{
    break;
  }
  case plstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    h[0]=norm;
    break;
  }
  case limplstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    //  first invariant of the stress tensor
    inv = first_invar (sigt);
    
    switch (ivhs){
    case limplstrnorm:{
      h[0]=norm/epspu;
      break;
    }
    case limanplstrnorm:{
      if (inv<0.0)
	h[0]=norm/epspuc;
      else
	h[0]=norm/epsput;
      break;
    }
    default:{
      fprintf (stderr,"\n unknown type of limit plastic strain is required in function hvalues (file %s, line %d)",__FILE__,__LINE__);
    }
    }
    
    break;
  }
  default:{
    fprintf (stderr,"\n unknown type of hardening/softening is required in function hvalues (file %s, line %d)",__FILE__,__LINE__);
  }
  }
  */
}


/**
   function computes derivatives of hardening/softening function with respect
   to stress components
   
   @param sigt - stress components stored in 3x3 %matrix
   @param dgds - first derivatives of plastic potential with respect to stress components (stored in 3x3 %matrix)
   @param dgdsds - second derivatives of plastic potential with respect to stress components (stored in 6x6 %matrix)
   @param dhds - derivatives of hardening/softening function with respect to stress components (stored in %matrix 6 x ncomphard)
   
   JK, 17.2.2007
*/
void hardsoft::dhdsigma (vector &/*sigt*/,vector &/*dgds*/,matrix &/*dgdsds*/,vector &/*dhds*/)
{
  /*
  long i,j;
  double norm,s,inv;

  //  cleaning output matrix
  fillv (0.0,dhds);
  
  switch (ths){
  case nohs:{
    break;
  }
  case plstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    
    for (i=0;i<6;i++){
      s=0.0;
      for (j=0;j<6;j++){
	s+=dgdsds[i][j];
      }
      dhds[i][0]=s/norm;
    }
    
    break;
  }
  case limplstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    //  first invariant of the stress tensor
    inv = first_invar (sigt);
    
    for (i=0;i<6;i++){
      s=0.0;
      for (j=0;j<6;j++){
	s+=dgdsds[i][j];
      }

      switch (ivhs){
      case limplstrnorm:{
	dhds[i][0]=s/norm/epspu;
	break;
      }
      case limanplstrnorm:{
	if (inv<0.0)
	  dhds[i][0]=s/norm/epspuc;
	else
	  dhds[i][0]=s/norm/epsput;
	break;
      }
      default:{
	fprintf (stderr,"\n unknown type of limit plastic strain is required in function dhdsigma (file %s, line %d)",__FILE__,__LINE__);
      }
      }

    }
    
    break;
  }
  default:{
    fprintf (stderr,"\n unknown type of hardening/softening is required in function dhdsigma (file %s, line %d)",__FILE__,__LINE__);
  }
  }
  */
}

/**
   function computes derivatives of hardening/softening function with respect to internal parameters
   
   @param sigt - stress components stored in 3x3 %matrix
   @param dgds - first derivatives of plastic potential with respect to stress components (stored in 3x3 %matrix)
   @param dgdsdq - second derivatives of plastic potential with respect to stress components and hardening parameters (stored in 6 x ncomphard %matrix)
   @param dhdq - derivatives of hardening/softening function with respect to internal parameters (stored in %matrix ncomphard x ncomphard)
   
   JK, 17.2.2007
*/
void hardsoft::dhdqpar (vector &/*sigt*/,vector &/*dgds*/,matrix &/*dgdsdq*/,vector &/*dhdq*/)
{
  /*
  long j;
  double norm,s,inv;
  

  //  cleaning output matrix
  fillm (0.0,dhdq);
  
  switch (ths){
  case nohs:{
    break;
  }
  case plstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    
    s=0.0;
    for (j=0;j<6;j++){
      s+=dgdsdq[j][0];
    }
    dhdq[0][0]=s/norm;
    
    break;
  }
  case limplstrainnorm:{
    //  norm of first derivatives of plastic potential with respect to stress components
    norm = tensornorm (dgds);
    //  first invariant of the stress tensor
    inv = first_invar (sigt);
    
    s=0.0;
    for (j=0;j<6;j++){
      s+=dgdsdq[j][0];
    }
    
    switch (ivhs){
    case limplstrnorm:{
      dhdq[0][0]=s/norm/epspu;
      break;
    }
    case limanplstrnorm:{
      if (inv<0.0)
	dhdq[0][0]=s/norm/epspuc;
      else
	dhdq[0][0]=s/norm/epsput;
      break;
    }
    default:{
      fprintf (stderr,"\n unknown type of limit plastic strain is required in function dhdqpar (file %s, line %d)",__FILE__,__LINE__);
    }
    }
    
    break;
  }
  default:{
    fprintf (stderr,"\n unknown type of hardening/softening is required in function dhdqpar (file %s, line %d)",__FILE__,__LINE__);
  }
  }
  */
}

/**
   function computes derivatives of hardening/softening function with respect
   to consistency parameter gamma
   
   @param dhdg - derivatives of hardening functions with respect to consistency parameter stored in ncomphard x 1 %vector
   
   JK, 17.2.2007
*/
void hardsoft::dhdgamma (vector &dhdg)
{
  fillv (0.0,dhdg);
}

/**
   function changes material parameters in stochastic or fuzzy computations

   particular material parameters are order in the following way

   limit plastic strain                  epspu
   limit plastic strain in tension       epsput
   limit plastic strain in compression   epspuc

   JK, 23.8.2005
*/
void hardsoft::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=atm.nba;i<atm.num;i++){
    switch (atm.atrib[i]){
      
      //  limit plastic strains
    case 0:{
      epspu=val[i];
      break;
    }
    case 1:{
      epsput=val[i];
      break;
    }
    case 2:{
      epspuc=val[i];
      break;
    }
      
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

}
