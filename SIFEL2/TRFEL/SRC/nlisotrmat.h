#ifndef NLISOTRMAT_H
#define NLISOTRMAT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "aliast.h"
#include "dampermeability.h" 
struct vector;
struct atsel;

class nlisotrmat
{
 public:
  nlisotrmat (void);
  ~nlisotrmat (void);
 
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void give_dof_names(namevart *dofname, long ntm);
  void read (XFILE *in);
  void print (FILE *out);
  double get_k (double v);
  double get_c();
  //void changeparam (atsel &atm,vector &val);
 
  
  ///  coefficient of conductivity has the form a*v+b
  ///  v - unknown variable, in the case of heat transfer, v denotes the temperature
  double a,b;

  //  coefficient of capacity
  double c;
  
 /* double omega;
  double lambda;
  double epsilon;
  double sigma;
  double w; */
  //  coefficient of conductivity
  gfunct  k;

  
};

#endif
