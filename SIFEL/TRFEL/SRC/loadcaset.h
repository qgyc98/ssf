#ifndef LOADCASET_H
#define LOADCASET_H

#include <stdio.h>
#include "iotools.h"
#include "pvalt.h"
#include "loadelt.h"
#include "sourcet.h"
#include "climatcond.h"
#include "climatcond2.h"
#include "bnodvalt.h"
#include "elemnode.h"

class loadcaset
{
 public:
  loadcaset (void);
  ~loadcaset (void);
  void read (XFILE *in,long lcid);
  void print (FILE *out,long lcid);
  void assemble (long lcid,double *rhs);
  void assemble_flux (long lcid,double *rhs,long n);
  double give_fact (double t, long i);
  void elemsource ();
  void sourcenodalvalues (long idse,vector &nodval);

  void source_contrib (long lcid,double *rhs);

  // *************************************************
  //  DESCRIPTION OF BOUNDARY CONDITIONS ON ELEMENTS
  // *************************************************
  ///  number of element with boundary conditions
  long neb;
  ///  number of objects describing nodal values
  long nnv;
  ///  nodal values of variable fluxes or external values
  bnodvalt *nodval;
  ///  number of objects describing climatic conditions
    long ncc,ncc2;
  ///  detailed climatic conditions
  climatcond *climcond;
  climatcond2 *climcond2;
  ///  type of reading of climatic conditions
  ///  trcc=1 - climatice conditions are described in the input file
  ///  trcc=2 - climatic conditions are desrcibed in additional files
  long trcc;
  long trcc2;

  
  // *************************
  //  DESCRIPTION OF SOURCES
  // *************************
  ///  number of quantity sources
  long nqs;
  
  ///  number of nodes with defined quantity source
  long nnqs;
  ///  list of nodes with defined quantity source
  ///  lnqs[nid][qsid] - at node nid (node id) is defined source of quantity with number qsid (quantity source id)
  long **lnqs;

  ///  number of elements with defined quantity source
  long neqs;
  ///  list of elements with defined quantity source
  ///  leqs[eid][qsid] - at element eid (element id) is defined source of quantity with number qsid (quantity source id)
  long **leqs;

  ///  number of elements influenced by nodes with prescribed source
  long neins;
  
  ///  number of nodes with defined point quantity source
  long nnpqs;
  ///  list of nodes with defined point quantity source
  ///  lnpqs[nid][qsid] - at node nid (node id) is defined point source of quantity with number qsid (quantity source id)
  long **lnpqs;

  ///  objects describing sources
  sourcet *sour;
  
  ///  class contains map between selected nodes and elements
  ///  it is used for sources defined at nodes and at elements
  ///  JK, 20.10.2007
  elemnode *en;

  ///  preprocessing for sources
  ///  sip = 0 - no preprocessing is required
  ///  sip = 1 - preprocessing is required
  long sip;
  ///  list of nodes in cluster where the same sources are defined
  long **lnc;
  ///  array containing numbers of nodes in clusters
  long *nc;
  /// number of nodes
  long ne;
  
  
  
  // *******************************
  //  prescribed values in nodes
  // *******************************
  ///  number of prescribed values
  long npv;
  ///  prescribed values
  pvalt *pv;
  
  // *******************************************************
  //  elements with Neumann and Newton boundary conditions
  // *******************************************************
  loadelt *elemload;
  
  
  
  //  type of time function
  generalfunct tfunc;
  //  function defined by table
  tablefunct *tabf;
  long nmultfunc;//total number of functions
  
  
  /// values in master node
  double masterval;    //unknown value
  double *mastergrad;  //unknown gradient
  int ndim;            //dimension of problem
};

#endif
