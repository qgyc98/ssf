#ifndef LGNODE_H
#define LGNODE_H

#include <stdio.h>

class lgnode
{
 public:
  lgnode (void);
  ~lgnode (void);
  
  //  number of layers
  long nl;
  //  number of multipliers
  long nmult;
  //  number of DOFs of one-layer node
  long ndofn;
  //  list of one-layer node numbers
  long *nodes;
  //  array containing code numbers
  long **cn;

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
