#ifndef HOMOGMAT_H
#define HOMOGMAT_H

#include <stdio.h>
#include "genfile.h"
#include "aliast.h"
struct vector;

class homogmat
{
 public:
  homogmat (void);
  ~homogmat (void);
 
  /// reads material parameters from the file
  void read (XFILE *in);
  /// prints material parameters to the file
  void print (FILE *out);

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void assemble_matrices (double *d,long ntm,long dim);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  
  double get_transmission_transcoeff_ww(double w,double t,long bc,long ipp);
  double get_transmission_nodval_ww(double bv,double w,double t,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double w,double t,long bc,long ipp);
  
  double get_transmission_transcoeff_tt(double w,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double w,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double w,double t,long bc,long ipp);
  
  double get_transmission_transcoeff_tw(double w,double t,long bc,long ipp);
  double get_transmission_nodval_tw(double bv,double w,double t,long bc,long ipp);
  double get_transmission_flux_tw(double bv,double w,double t,long bc,long ipp);

  double get_othervalue(long compother,double rh,double t,long ipp);
  void print_othervalue_name(FILE *out,long compother);

  void give_dof_names(namevart *dofname, long ntm);

   /// fills array with non-transport quantities required by the model
  void give_reqntq(long *antq); 

  matrix dd;
  matrix cc;

  long hom_mattype; //material type for mogenization
  long hom_mattype_number; //number of material for homogenization
};

#endif
