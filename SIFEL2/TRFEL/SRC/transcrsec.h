#ifndef TRANSCRSEC_H
#define TRANSCRSEC_H

#include <stdio.h>
#include "aliast.h"
#include "crsection1d.h"
#include "crsection2d.h"
#include "crsection3d.h"
#include "genfile.h"

class transcrsec
{
 public:
  transcrsec (void);
  ~transcrsec (void);
  void read (XFILE *in);
  void readcrsectype(XFILE *in, crsectypet ct, long numt);
  void print (FILE *out);
  void printcrschar (FILE *out, crsectypet ct, long numinst);
  void give_thickn (ivector &nod,vector &t);
  void give_thicke (long eid,double &t);
  void give_thickness (long eid,ivector &nodes,vector &th);
  void give_areae (long eid,double &a);
  void give_arean (ivector &nod,vector &a);
  void give_area (long eid,ivector &nodes,vector &a);
  void give_densityn (ivector &nod,vector &rho);
  void give_densitye (long eid,double &rho);
  void give_density (long eid,ivector &nodes,vector &dens);

  
  ///  number of cross section types
  long ncst;
  ///  type of cross sections
  crsectypet *cstype;
  ///  number of instances of particular cross sections
  long *numtype;
  
  ///  cross section of 1D elements
  crsection1d *cs1d;
  ///  cross section of 2D elements
  crsection2d *cs2d;
  ///  cross section of 3D elements
  crsection3d *cs3d;

};

#endif
