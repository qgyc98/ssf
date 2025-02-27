#ifndef INTPOINTS_H
#define INTPOINTS_H

#include <stdio.h>
#include "alias.h"
#include "selection.h"
#include "iotools.h"

/**
  Class intpoints:
   The class defines structure of integration point which is used in the mechmat class.
   Integeration point holds data about used materials and it also 
   contains arrays of strains, stresses, internal variables of materialmodels and nonlocal values.
   Other array is used for internal variables and its size is defined by the used model.
   
  Created by JK, Tomas Koudelka
*/

class intpoints
{
 public:
  intpoints (void);
  ~intpoints (void);
  void read (XFILE *in);

  void alloc (long nlc,long ipp,long ncompnl);
  void alloc_strains (long nlc);
  void alloc_stresses (long nlc);
  //void alloc_other (long ipp);
  void copy(intpoints &ip, long nlc, long ncompnl, long realloc);

  long gemid(void);
  long grmid(void);

  void clean (long nlc);
  void clean_strains (long nlc);
  void clean_stresses (long nlc);
  void clean_other ();

  void save_data_txt    (FILE *aux,long nlc, sel &selother);
  void restore_data_txt (FILE *aux,long nlc, long ncompo, sel &selother, long *selid);
  void save_data_bin    (FILE *aux,long nlc, sel &selother);
  void restore_data_bin (FILE *aux,long nlc, long ncompo, sel &selother, long *selid);

  long give_ncompstr ();


  ///  type of material
  mattype *tm;
  ///  number of appropriate material type
  long *idm;
  ///  number of material types defined at integration point
  long nm;
  /// bit array with material type presence flags (i.e nonlocal model, thermal dilatancy)
  long hmt;
  ///  number of component of stress/strain array
  long ncompstr;
  ///  number of component of eqother array
  long ncompeqother;
  ///  number of component of other array
  long ncompother;
  ///  stress-strain state
  strastrestate ssst;

  ///  stress components
  double *stress;
  ///  strain components
  double *strain;
  ///  other components
  double *other;
  ///  equilibriated components of other array
  double *eqother;
  ///  nonlocal values
  double *nonloc;

  /*
    PLANE STRESS PROBLEMS
    strain components   e_xx, e_yy, e_xy
    stress components   s_xx, s_yy, s_xy
                   
    PLANE STRAIN PROBLEMS
    strain components   e_xx, e_yy, e_xy
    stress components   s_xx, s_yy, s_zz, s_xy
    
    3D PROBLEMS
    strain components   e_xx, e_yy, e_zz, e_yz, e_zx, e_xy
    stress components   s_xx, s_yy, s_zz, s_yz, s_zx, s_xy
    

    ARRAY OTHER FOR VISCO-PLASTIC PROBLEMS
    other[0*ncompstr - (1*ncompstr-1)] - real stresses
    other[1*ncompstr - (2*ncompstr-1)] - irreversible strain increments
    other[2*ncompstr - (3*ncompstr-1)] - previous total strains
    other[3*ncompstr - ncompother] - hardening parameters
    
    
    
    if thermal dilatancy is necessary, it is always the last component in the tm array and
    the variable hmt is equal to 1, otherwise is equal to 0
    
    position of elastic material in the arrays tm or idm is defined by function gemid ()
    ip[ipp].gemid ()
    position of nonlocal material in the arrays tm or idm is defined by function gnlmid ()
    ip[ipp].gnlmid ()
    
    ELASTIC MODELS
    tm[0] - type of elastic material
    (tm[1] - type of thermal dilatancy)
    
    PLASTIC MODELS
    tm[0] - type of plastic model (type of yield function)
    tm[1] - type of elastic model
    (tm[2] - type of thermal dilatancy)
    
    VISCO-PLASTIC MODELS
    tm[0] - type of viscous material model
    tm[1] - type of plastic model (type of yield function)
    tm[2] - type of elastic model
    (tm[3] - type of thermal dilatancy)
    
    SCALAR DAMAGE MODELS
    tm[0] - type of damage model
    tm[1] - type of elastic model
    (tm[2] - type of thermal dilatancy)

  */

};

#endif
