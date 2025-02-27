#ifndef COUPTOP_H
#define COUPTOP_H

#include <stdio.h>
#include "aliasc.h"
#include "alias.h"
#include "nodec.h"
#include "elementc.h"
#include "iotools.h"

class node;
class element;
struct vector;
struct ivector;


/**
   CLASS COUPTOP
   
   it is one of the 5 most important classes of the program

   class couptop contains topology data of problem
   
   basic data
   nn (long number) - stands for number of nodes in problem
   ncn (long) - stands for number of constrained nodes
   ne (long) - stands for number of all finite elements in problem
   
   class mechtop creates
   nn objects of class node
   ne objects of class element
   
   
   JK, TK
*/
class couptop
{
 public:
  couptop (void);
  ~couptop (void);
  void read (XFILE *in);
  long mesh_check(void);
  
  elemtypec give_elem_type (long eid);
  long give_ndofn (long nid);
  long give_ndofe (long eid);
  long give_nne (long eid);
  //long give_nne_inner (elemtypec te);
  long give_upper_tnip (long eid);
  long give_lower_tnip (long eid);
  long give_upper_nip (long eid,long ri,long ci);
  long give_lower_nip (long eid,long ri,long ci);
  void give_elemnodes (long eid,ivector &nodes);
  void give_code_numbers (long eid,long *cn);
  void give_node_coord2d (vector &x,vector &y,long eid);
  void give_node_coord2dxz (vector &x,vector &z,long eid);
  void give_node_coord3d (vector &x,vector &y,vector &z,long eid);
  
  long give_mnb (long eid);
  strastrestate give_ssst (long eid,long bi);
  //  function determines the number of boundary objects (edges, surfaces) and the number of nodes on a boundary object
  void give_nbobjects (long eid,long &bco,long &ncompbco);

  ///  function returns the number of strain/stress components
  long give_ncompstr (long eid);
  ///  function returns the number of displacement components
  long give_ncompdispl (long eid);
  ///  function returns the number of transported media
  long give_ntm (long eid);
  ///  function returns the number of gradient/flux components
  long give_ncompgrad (long eid);
  
  long give_nip (long eid);
  
  
  ///  the number of nodes
  long nn;
  ///  nodes
  nodec *nodes;
  
  ///  number of nodes with prescribed values
  ///  the number of nodes with Dirichlet boundary condition
  long nnd;

  ///  number of elements
  long ne;
  ///  finite elements
  elementc *elements;
  
};

#endif
