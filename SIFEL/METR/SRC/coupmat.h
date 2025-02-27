#ifndef COUPMAT_H
#define COUPMAT_H

#include <stdio.h>
#include "aliasc.h"
#include "intpointsc.h"
#include "isotrmatc.h"
#include "sejtkrmatc.h"
#include "baroghelBc.h"
#include "concreteBc.h"
#include "C60baroghelc.h"
#include "C30baroghelc.h"
#include "o30bazantc.h"
#include "C60bazantc.h"
#include "glasgowmatc.h"
#include "genfile.h"
#include "glasgowcoup.h"
#include "consol_awf1c.h"
#include "consol_awf2c.h"
#include "consol_hwf2c.h"
#include "consol_hawf3c.h"
#include "bentonitemat.h"
#include "xfile.h"

class coupmat
{
 public:
  coupmat (void);
  ~coupmat (void);

  void readmatchar (XFILE *in);
  void intpointalloc ();
  void intpointinit ();
  void read (XFILE *in);

  long givencompother (long ipp,long im);
  long givencompeqother (long ipp,long im);

 /*
  long intpnum (void);
  void readip (XFILE *in);
  void matcond (matrix &d,long ipp,long ri,long ci);
  void matcap (matrix &d,long ipp,long ri,long ci);

  void computenlfluxes (matrix &d,long lcid,long ipp);
  void storestrain_cml (long ipp,long fi,vector &eps);
  void storegrad_cml (long lcid,long ipp,vector &gr);
  void givegrad_cml (long lcid,long ipp,vector &gr);
  void storeflux_cml (long lcid,long ipp,vector &fl);
  void givefluxes_cml (long lcid,long ipp,vector &fl);
  */
  
  ///  stiffness matrix
  void matstiff (matrix &d,long ipp);
  ///  conductivity/permeability matrix
  void matcond (matrix &d,long ipp);
  ///  c_pp coefficient
  double c_pp_coeff (long ipp);
  ///  c_up coefficient
  double c_up_coeff (long ipp);
  ///  c_pu coefficient
  double c_pu_coeff (long ipp);

  void initvalues (long ipp);

  void storeother (long ipp,long fi,long ncomp,double *comp);
  void storestrain (long ipp,vector &eps);

  long intpnum_fc (void);  
  
  ///  number of material types
  long nmt;
  ///  total number of all integration points
  long tnip;
  ///  pointer to integration points
  intpointsc *ip;

  /// array with initial conditions at regular integration points (allocated at mechbclc::readinic)
  double **ic;
  ///  array containing eigenstrains id
  ///  the array is allocated in function mechbclc::read
  long **eigstrid;
  ///  array containing eigenstrains
  ///  the array is allocated only for several problems
  ///  eigstrains[number of int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **eigstrains;
  ///  array containing eigenstresses
  ///  the array is allocated only for several problems
  ///  eigstresses[number of int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **eigstresses;





  //  COUPLING PROBLEMS
  //  one-phase medium - heat transfer - isotropic heat-elast material
  isotrmatc *itrmc;
  //  saturate-unsaturated one-phase flow in deforming meduim
  sejtkrmatc *sejtkrmc;
  //three-phase medium
  baroghelmatc *baroghelc;
  concreteBmatc *concretec;
  C60barmatc  *C60baroghelc;
  C30barmatc  *C30baroghelc;
  o30bazmatc  *o30bazantc;
  C60bazmatc  *C60bazantc;
  glasgowmatc *tenchc;
  //one-phase flow
  con_awf1matc *consol_awf1c;  
  //two-phase flow
  con_awf2matc *consol_awf2c;  
  con_hwf2matc *consol_hwf2c;  
  //three-phase flow
  con_hawf3matc *consol_hawf3c;
  glasgowcoup *gcm;
  
  ///  fully coupled material model for bentonite
  ///  Masin hypoplastic model is completed by ingredients of Schrefler model
  bentonitemat *bento;

  /**  integration point - element correspondation
       elip[number of int point] = number of element */
  long *elip;
 
};

#endif
