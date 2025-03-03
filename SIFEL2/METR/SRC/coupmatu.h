#ifndef COUPMATU_H
#define COUPMATU_H

#include <stdio.h>
#include "alias.h"
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
#include "glasgowcoup.h"
#include "consol_awf1c.h"
#include "consol_wf1c.h"
#include "consol_wf2c.h"
#include "consol_awf2c.h"
#include "consol_hawf3c.h"
#include "consol_hwf2c.h"
#include "xfile.h"

class coupmatu
{
 public:
  coupmatu (void);
  ~coupmatu (void);
  void ipalloc (void);
  long intpnum (void);
  void readip (XFILE *in);
  void readmatchar (XFILE *in);
  void intpointalloc ();
  void intpointinit ();
  void initmaterialmodels ();
  void initvalues (long ipp,long im,long ido);
  void updateipval (void);
  void updateipvalmat (long ipp,long im,long ido);
  void read (XFILE *in);
  void matstiff (matrix &d,long ipp);
  void matcond (matrix &d,long ipp,long ri,long ci);
  void matcap (matrix &d,long ipp,long ri,long ci);
  void volume_rhs1 (matrix &d,long ipp,long ri,long ci,long ncomp);
  void volume_rhs2 (matrix &d,long ipp,long ri,long ci);

  void computenlstresses (matrix &d,long ri, long ci, long ipp);
  void storestresses_cmu (long ipp,long fi,vector &gr);
  void givestresses_cmu (long ipp,long fi,vector &fl);
  void storestrain_cmu (long ipp,long fi,vector &eps);
  void givestrain_cmu (long ipp,long fi,vector &eps);
  void storegrad_cmu (long lcid,long ipp,vector &gr);
  void givegrad_cmu (long lcid,long ipp,vector &gr);
  void storeflux_cmu (long lcid,long ipp,vector &fl);
  void givefluxes_cmu (long lcid,long ipp,vector &fl);

  //  number of material types
  long nmt;
  //  total number of all integration points
  long tnip;
  //  pointer to integration points
  intpointsc *ip;

  //  COUPLING PROBLEMS
  //  one-phase medium - heat transfer - isotropic heat-elast material
  isotrmatc *itrmc;
  //  saturated-unsaturated one-phase flow in deforming meduim
  sejtkrmatc *sejtkrmc;
  con_awf1matc *consol_awf1c;
  con_wf1matc *consol_wf1c;
  //two-phase flow
  con_wf2matc *consol_wf2c;
  con_awf2matc *consol_awf2c;
  con_hwf2matc *consol_hwf2c;
  //three-phase flow
  con_hawf3matc *consol_hawf3c;
  baroghelmatc *baroghelc;
  concreteBmatc *concretec;
  C60barmatc  *C60baroghelc;
  C30barmatc  *C30baroghelc;
  o30bazmatc  *o30bazantc;
  C60bazmatc  *C60bazantc;
  glasgowmatc *tenchc;
  glasgowcoup *gcm;

};

#endif
