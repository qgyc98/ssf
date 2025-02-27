#ifndef CONSTRELCU_H
#define CONSTRELCU_H

#include "genfile.h"
#include "alias.h"

class state_eqcu
{
 public:
  state_eqcu();//constructor
  ~state_eqcu();   //destructor
  
  //DEGREE OF SATURATION AND DERIVATIVES
  double get_s(double pc,double pg,double t,long ipp);
  // BIOT'S CONSTANT
  double get_alpha(double pc,double pg,double t,long ipp);
  //MATERIAL PROPERTIES
  double get_rhos(double pc,double pg,double t,long ipp);
  double get_betas(double pc,double pg,double t,long ipp);
  //STIFFNESS MATRIX
  void matstiff (matrix &d,strastrestate mssst, long ipp);
  void matstiff_bar (matrix &d,long ipp);
  void matstiff_plstrain (matrix &d,long ipp);
  void matstiff_plstress (matrix &d,long ipp);
  void matstiff_spacestress (matrix &d,long ipp);
  void matstiff_axi (matrix &d,long ipp);
  double give_e(long ipp);
  double give_nu(long ipp);

 private: 
  double scale_pc,scale_pg,scale_t,scale_u;

};  

#endif
