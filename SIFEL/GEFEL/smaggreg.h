#ifndef SMAGGREG_H
#define SMAGGREG_H

#include "locmatrix.h"
#include "aggregator.h"
#include "gmatrix.h"
#include "gtopology.h"

/**
   class contains solver of system of linear algebraic equations based
   on smoothed aggregation concept
   
   JK, 23.12.2006
*/
class smaggreg
{
 public:
  smaggreg ();
  ~smaggreg ();
  
  //void aggregates ();
  //void tentative_prolongator ();
  void localization_matrices (gtopology &gt);
  void local_aggregated_matrices (gmatrix &sm);
  
  void decompose_aggreg_matrix (gmatrix &sm);
  void smoothed_prolong ();

  ///  number of unknowns
  long nu;
  ///  number of aggregates
  long na;
  ///  kernel dimension for each aggregate
  long dimker;
  
  ///  list of numbers of nodes in aggregates
  ///  nnagr[i]=j - the i-th aggregate contains j nodes
    //long *nnagr;

  ///  list of nodes in aggregates
  ///  lnnagr[i][j]=k - the j-th node in the i-th aggregate has number k
    //long **lnnagr;

    ///  list of number of unknowns on aggregates
    ///  nuagr[i]=j - the i-th aggregate contains j unknowns
    long *nuagr;
    
    ///  list of unknown id's on aggregates
    ///  uagr[i][j]=k - the j-th unknown on the i-th aggregate has number k
      ///  k is number if global ordering
      long **uagr;

  ///  tentative prolongator
  ///  (Brezina, page 68, plongator is defined columnwise there)
  ///  this formulation assembles transposed prolongator because of
  ///  structure of class localization %matrix (locmatrix.h)
    //locmatrix tp;
	
  ///  smoothed prolongator
  /// (Brezina, page 69, plongator is defined columnwise there)
  ///  this formulation assembles transposed prolongator because of
  ///  structure of class localization %matrix (locmatrix.h)
  locmatrix sp;

  ///  localization matrices
  locmatrix *lm;
    

  ///  local aggregated matrices
  gmatrix *gm;
  
  ///  aggregator
  aggregator aggreg;
};

#endif
