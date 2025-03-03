#ifndef MULTIPHASEC_H
#define MULTIPHASEC_H

#include "genfile.h"

class multiphc
{
 public:
  multiphc();    //constructor
  ~multiphc();   //destructor

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond1d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_u (matrix &d,long ri,long ci,long ipp);
  void matcond3d_u (matrix &d,long ri,long ci,long ipp);
  void matcap1d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_u (matrix &d,long ri,long ci,long ipp);
  void matcap3d_u (matrix &d,long ri,long ci,long ipp);

  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void matcond1d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_l (matrix &d,long ri,long ci,long ipp);
  void matcond3d_l (matrix &d,long ri,long ci,long ipp);
  void matcap1d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_l (matrix &d,long ri,long ci,long ipp);
  void matcap3d_l (matrix &d,long ri,long ci,long ipp);
  
  void gaspress_check(double pc,double &pg,double t,long ipp);
  void cappress_check(double &pc,double pg,double t,long ipp);
  void cappress_stop(double &pc,double pg,double t,long ipp);

  void values_correction (vector &nv);
  double get_kcu(double pc,double pg,double t,long ipp);
  double get_capcu(double pc,double pg,double t,long ipp);
  double get_kgu(double pc,double pg,double t,long ipp);
  double get_capgu(double pc,double pg,double t,long ipp);
  double get_ktu(double pc,double pg,double t,long ipp);
  double get_captu(double pc,double pg,double t,long ipp);

  double get_kug(double pc,double pg,double t,long ipp);
  double get_capug(double pc,double pg,double t,long ipp);
  double get_kuc(double pc,double pg,double t,long ipp);
  double get_capuc(double pc,double pg,double t,long ipp);
  double get_kut(double pc,double pg,double t,long ipp);
  double get_caput(double pc,double pg,double t,long ipp);

  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs1d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d_u1 (matrix &d,long ri,long ci,long ipp);

  double get_fuc1(double pc,double pg,double t,long ipp);
  double get_fug1(double pc,double pg,double t,long ipp);
  double get_fut1(double pc,double pg,double t,long ipp);

  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs1d_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs3d_u2 (matrix &d,long ri,long ci,long ipp);

  double get_fu2(double pc,double pg,double t,long ipp);

 private: 

  double scale_pc,scale_pg,scale_t;

};  

#endif
