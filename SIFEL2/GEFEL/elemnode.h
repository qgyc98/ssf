#ifndef ELEMNODE_H
#define ELEMNODE_H

#include "gtopology.h"

/**
   class deals with correspondence between nodes and elements
   it assembles and contains list of nodes shared by selected elements
   it assembles and contains list of elements which contain selected nodes
   
   JK, 20.10.2007
*/

class elemnode
{
 public:
  elemnode ();
  ~elemnode ();
  
  void selelem (long nselelem,long *lselem);
  void selnode (long nselnod,long *lselnod);
  
  void elemnodes (gtopology *gt);
  
  ///  number of selected elements
  long nse;
  
  ///  list of selected elements
  ///  lse[i]=j - the i-th selected element has number j
  long *lse;
      
  ///  number of selected nodes
  long nsn;
  
  ///  list of selected nodes
  ///  lsn[i]=j - the i-th selected node has number j
  long *lsn;
  
  
  ///  number of influenced elements
  long nie;

  ///  list of elements influenced by selected nodes, it contains number of selected node
  ///  in other words, it contains positions in array lsn
  ///  elnod[i][j] = k - the j-th node on the i-th influenced element has number k
  ///  k=-1 - the j-th node is not selected, k>-1 - number of selected node
  ///  example: list of selected nodes lsn = 1, 4, 7, 12
  ///  elnod[3][1]=2 - second node on the fourth element has number 7, because lsn[2]=7
  ///  array elnod has nie x nne[i] components
  long **elnod;

};

#endif
