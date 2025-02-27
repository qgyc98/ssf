#ifndef AEPOINTS_H
#define AEPOINTS_H

#include <stdio.h>
#include "alias.h"
struct vector;


class aepoints
{
 public:
  aepoints ();
  ~aepoints ();
  void read (FILE *in);
  void init(strastre ssf);
  long give_naep (long eid);
  long give_ncomp (long eid);
  long give_sid (long eid);
  void give_aepcoord (long sid,long pid,vector &coord);
  void alloc (long nlc);
  void storevalues (long lcid,long eid,long pid,vector &val);
  void transformvalues (long tt);

  
  //  array containing type of auxiliary points on elements
  long *tape;
  //  array containing pointers to type of auxiliary points on elements
  long *ptape;
  //  array containing number of auxiliary points on elements
  long *nape;
  //  number of user defined sets of points
  long nudsets;
  //  auxiliary array containing number of points and components in one set
  long **udpa;
  //  array containing auxiliary points coordinates
  double ***udpc;
  
  //  array containing evaluated variables
  double ***ev;
  
  //  number of elements with transformation
  long net;
  //  array containing element numbers
  long *ent;
  //  array containing numbers of points with transformation
  long *npt;
  //  array containing point numbers
  long **pt;
  //  array containing pointers to local bases
  long **plcs;
  
  //  number of local coordinate systems
  long nlcs;
  //  array containing numbers of components of local coordinates systems
  long *nclcs;
  //  array containing base vectors of local coordinate systems
  double **lcs;
  
  
};
#endif
