#ifndef MECHCRSEC_H
#define MECHCRSEC_H

#include "iotools.h"
#include "alias.h"
class crsec2dbar;
class crsec2dbeam;
class crsec3dbeam;
class crsecplstr;
class crsec3d;
class crsecnod;
class crseclayer;
struct ivector;
struct vector;



/**
  Class mechcrsec:
   
  It is one of the 5 most important classes of the program.
  (probdesc, mechtop, mechmat, mechbclc, mechcrsec)
   
  Class mechcrsec contains data about cross sections.
   
  Created by JK, TKo
*/
class mechcrsec
{
 public:
  mechcrsec ();
  ~mechcrsec ();
  void read (XFILE *in);
  void readcrsectype(XFILE *in, crsectype cstype, long numtype);
  void print (FILE *out);
  void printcrschar (FILE *out, crsectype ct, long numinst);

  void give_thickn (ivector &nod,vector &t);
  void give_thicke (long eid,double &t);
  void give_thickness (long eid,ivector &nodes,vector &th);
  double give_onethickness (crsectype crst,long idcs);
  
  void give_area (long ipp,double &a);
  void give_mominer (long ipp,double *i);
  void give_shearcoeff (long ipp,double *shearcoeff);
  void give_vectorlcs (long eid,vector &vec);

  void give_densityn (ivector &nod,vector &rho);
  void give_densitye (long eid,double &rho);
  void give_density (long eid,ivector &nodes,vector &dens);
  
  double give_weight (crsectype cr,long idcs);

  double* give_layer_thicke (long eid);
  double* give_layer_zcoord (long eid);
  long    give_num_lay (long eid);
  
  ///  number of cross section types
  long ncst;
  ///  type of cross sections
  crsectype *cstype;
  ///  number of instances of particular cross sections
  long *numtype;
  
  ///  cross section of 2D bar elements
  crsec2dbar *cs2dbar;
  ///  cross section of 2D beam elements
  crsec2dbeam *cs2dbeam;
  ///  cross section of 3D beam elements
  crsec3dbeam *cs3dbeam;
  ///  cross section for plane stress and plate problems
  crsecplstr *csplstr;
  ///  cross section for 3D problems
  crsec3d *cs3d;
  ///  cross section defined in nodes
  ///  used especially in layered static analysis
  crsecnod *csnod;
  ///  cross section for layered plate problems
  crseclayer *cslay;

};

#endif
