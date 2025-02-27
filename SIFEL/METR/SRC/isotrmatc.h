#ifndef ISOTRMATC_H
#define ISOTRMATC_H

#include <stdio.h>
#include "genfile.h"
struct vector;
struct atsel;

class isotrmatc
{
 public:
  isotrmatc (void);
  ~isotrmatc (void);
  void read (XFILE *in);
  void print(FILE *out);
  
  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond1d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_ax_u (matrix &d,long ri,long ci,long ipp);
  void matcond3d_u (matrix &d,long ri,long ci,long ipp);
  void matcap1d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_ax_u (matrix &d,long ri,long ci,long ipp);
  void matcap3d_u (matrix &d,long ri,long ci,long ipp);

  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void matcond1d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_ax_l (matrix &d,long ri,long ci,long ipp);
  void matcond3d_l (matrix &d,long ri,long ci,long ipp);
  void matcap1d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_ax_l (matrix &d,long ri,long ci,long ipp);
  void matcap3d_l (matrix &d,long ri,long ci,long ipp);

  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs1d2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_ax_2 (matrix &d,long ri,long ci,long ipp);
  void rhs3d2 (matrix &d,long ri,long ci,long ipp);

  double get_k();
  double get_c();
  double get_e();
  double get_nu();
  double get_alpha();
  void changeparam (atsel &atm,vector &val);

  //  coefficient of conductivity
  double k;
  //  coefficient of capacity
  double c;
  //Young's modulus
  double e;
  //Poisson's constant
  double nu;
  //  coefficient of thermal dilatation
  double alpha;
  


};

#endif
