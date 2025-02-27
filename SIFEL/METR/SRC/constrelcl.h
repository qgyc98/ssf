#ifndef CONSTRELCL_H
#define CONSTRELCL_H

#include "genfile.h"

class state_eqcl
{
 public:
  state_eqcl();//constructor
  ~state_eqcl();   //destructor
  
  //DEGREE OF SATURATION AND DERIVATIVES
  double get_s(double pc,double pg,double t,long ipp);
  
  // BIOT'S CONSTANT
  double get_alpha(double pc,double pg,double t,long ipp);
  //MATERIAL PROPERTIES
  double get_rhos(double pc,double pg,double t,long ipp);
  double get_betas(double pc,double pg,double t,long ipp);
  double get_krw(double pc,double t,long ipp);
  double get_kintr(double pc,double pg,double t,long ipp);
  
 private: 

};  

#endif
