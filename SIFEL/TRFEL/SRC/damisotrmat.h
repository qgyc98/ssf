#ifndef DAMISOTRMAT_H
#define DAMISOTRMAT_H

#include <stdio.h>
#include "genfile.h"
#include "dampermeability.h"
#include "aliast.h"


struct vector;
struct atsel;

class damisotrmat
{
 public:
  damisotrmat (void);
  ~damisotrmat (void);
 
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);
  void print (FILE *out);
  double get_k();
  double get_c();
  void print_othervalue_name(FILE *out,long compother);
  void give_dof_names(namevart *dofname, long ntm);
  void changeparam (atsel &atm,vector &val);
  void give_reqntq(long *antq);

  //  coefficient of permeability
  double k;
  //  coefficient of capacity
  double c;

  ///  influence of damage on permeability
  static dampermeability damper;

  ///  flag for influence of damage on permeability
  flagsw daminfl;
  
};

#endif
