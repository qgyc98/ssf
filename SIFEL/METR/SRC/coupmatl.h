#ifndef COUPMATL_H
#define COUPMATL_H

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
#include "consol_wf1c.h"
#include "consol_wf2c.h"
#include "consol_awf2c.h"
#include "consol_hwf2c.h"
#include "consol_hawf3c.h"
#include "xfile.h"

class coupmatl
{
 public:
  coupmatl (void);
  ~coupmatl (void);
  void ipalloc (void);
  long intpnum (void);
  void readip (XFILE *in);
  void readmatchar (XFILE *in);
  void intpointalloc ();
  void intpointinit ();
  void read (XFILE *in);
  void matcond (matrix &d,long ipp,long ri,long ci);
  void matcap (matrix &d,long ipp,long ri,long ci);

  void computenlfluxes (matrix &d,long lcid,long ipp);
  void storestrain_cml (long ipp,long fi,vector &eps);
  void storegrad_cml (long lcid,long ipp,vector &gr);
  void givegrad_cml (long lcid,long ipp,vector &gr);
  void storeflux_cml (long lcid,long ipp,vector &fl);
  void givefluxes_cml (long lcid,long ipp,vector &fl);

  //  number of material types
  long nmt;
  //  total number of all integration points
  long tnip;
  //  pointer to integration points
  intpointsc *ip;

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
  con_wf1matc *consol_wf1c;  
  //two-phase flow
  con_wf2matc *consol_wf2c;  
  con_awf2matc *consol_awf2c;  
  con_hwf2matc *consol_hwf2c;  
  //three-phase flow
  con_hawf3matc *consol_hawf3c;
  glasgowcoup *gcm;
};

#endif
