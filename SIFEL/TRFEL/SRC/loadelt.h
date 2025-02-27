#ifndef LOADELT_H
#define LOADELT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "climatcond.h"
#include "bnodvalt.h"

class loadelt
{
 public:
  loadelt ();
  loadelt (long id, long numbo);
  ~loadelt ();

  void read (XFILE *in,long lcid,long elemid,long inbo,long innbo);
  void read_comp (XFILE *in,long lcid,long elemid,long inbo,long innbo);
  void bc_indicators (long &transi);
  void updatevalues (long lcid, long &transi);
  void print(FILE *out,long lcid);

  long merge(loadelt &ebc);
  void renumber_id(long nnv, long ntrc);
  void give_bc (bocontypet *av);
  void give_nodval (long lcid,vector &fe);
  void give_external_nodval (long lcid,long cid,vector &fe);
  //void give_trc (long lcid,long cid,vector &trc);
  void give_trc (double time,bnodvalt *nodval,long lcid,long cid,vector &trc);
  //void give_trr (long lcid,vector &trr);
  void give_trr (double time,bnodvalt *nodval,vector &trr);

  ///  loaded element id (number of loaded element)
  long eid;
  ///  number of boundary objects (end nodes, edges, surfaces)
  long nbo;
  ///  number of nodes on boundary object
  long nnbo;
    
  ///  boundary conditions (for nongrowing problems)
  ///  bc=0 - inner object (end node, edge, surface)
  ///  bc=1 - object with prescribed values (Dirichlet boundary condition)
  ///  bc=2 - object with prescribed flux (Neumann boundary condition)
  ///  bc=3 - object with detailed climatic conditions
  ///  bc=4 - object with generalized climatic conditions
  ///  bc>10 - object with transmission boundary condition (Newton b.c.)
  bocontypet *bc;
  ///  boundary conditions (for growing problems)
  ///  bcf=0 - inner object (end node, edge, surface)
  ///  bcf=1 - object with prescribed values (Dirichlet boundary condition)
  ///  bcf=2 - object with prescribed flux (Neumann boundary condition)
  ///  bcf=3 - object with detailed climatic conditions
  ///  bcf=4 - object with generalized climatic conditions
  ///  bcf>10 - object with transmission boundary condition (Newton b.c.)
  long *bcf;

  ///  id of nodal values or climatic conditions (for nongrowing problems)
  ///  one id for each end node, edge or surface
  long *nvid;
  ///  id of nodal values or climatic conditions (for growing problems)
  ///  one id for each end node, edge or surface
  long *nvidf;
  ///  id of nodal values describing transmission coefficients (for nongrowing problems)
  ///  one id for each end node, edge or surface
  long *trcid;
  ///  id of nodal values describing transmission coefficients (for growing problems)
  ///  one id for each end node, edge or surface
  long *trcidf;
  ///  id of nodal values describing transmission/radiation coefficients (for nongrowing problems)
  ///  one id for each end node, edge or surface
  long *trrid;
  ///  id of nodal values describing transmission/radiation coefficients (for growing problems)
  ///  one id for each end node, edge or surface
  long *trridf;
  
};

#endif
