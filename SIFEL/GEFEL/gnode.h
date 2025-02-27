#ifndef GNODE_H
#define GNODE_H

#include <stdio.h>
#include "iotools.h"
#include "galias.h"
class gtopology;
class gfunct;

/**
   class gnode
   
   this class contains general topological informations about node
   
   JK
*/

class gnode
{
 public:
  gnode (void);
  ~gnode (void);
  void read (XFILE *in);
  void print (FILE *out);
  void constr (long dofcontr,XFILE *in);
  void print_constr (long dofcontr,FILE *out);
  long give_ndofn ();
  long give_dof (long m);
  void save_dof (long m,long num);
  double distance2 (long dim,const double *c);

  void update_dofs (gfunct *gf,double time, long lnso);
  /// function searches for changed DOFs
  long search_changed_dofs (gfunct *gf,double time,double prev_time, long lnso, long plnso);
  /// function sets DOF number array to zero values
  void clear_dof ();
  // function clears selected DOF number
  void clear_dof (long j);

  ///  number of DOFs on node
  long ndofn;
  ///  array containing code numbers
  long *cn;
  ///  coordinates of the node
  double x,y,z;
  ///  auxiliary idicator
  long ai;
  ///  numbers (pointers) of general functions
  long *tgf;
  
  ///  array of the master nodes (if the node is hanging node)
  long *mnodes;
  ///  natural coordinates of the hanging node
  double *natcoord;
  ///  type of master entity
  ///  it describes to which quantity (edge, surface, etc.) is the hanging node attached
  gtypel masentity;
};

/*
  v kazdem uzlu je treba vzhledem k obecnosti softwaru
  uchovavat ndofn, protoze existuji prvky s rozdilnymi
  pocty DOF v ruznych uzlech, napr. trojuhelnikova deska
  s polynomem 5. stupne ma ve vrcholech 6 neznamych, ale
  v uzlech uprostred stran ma pouze jeden stupen volnosti

  podpory jsou tim padem dany ve vstupnim souboru, protoze
  kazdy uzel obsahuje pole cn, kde 0 znamena braneno,
  1 znamena volno a zaporna cisla jsou po provedeni absolutni
  hodnoty pointry na predepsana posunuti
  
*/

#endif
