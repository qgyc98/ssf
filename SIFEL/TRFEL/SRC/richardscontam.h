#ifndef RICHARDSCONTAM_H
#define RICHARDSCONTAM_H

#include "genfile.h"
#include "richards.h"

class richardscontam
{
 public:
  richardscontam();    //constructor
  ~richardscontam();   //destructor
  
  void read(XFILE *in);
  void print(FILE *out);
  
  double cc_value (double h);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcond1d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond2d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond3d_2 (matrix &d,long ri,long ci,long ipp);
  

  ///  parameters of van Genuchten model
  double alpha;
  double n;
  double m;
  
  ///  dimension of problem solved
  long dim;

  ///  saturated hydraulic conductivities
  double kksxx,kksxy,kksxz,kksyy,kksyz,kkszz;
  ///  saturated water content
  double thetas;
  ///  residual water content
  double thetar;
  ///  specific storage
  double storage;
  
  ///
  richards rich;
};

#endif
