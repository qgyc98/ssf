#ifndef NODE_H
#define NODE_H

#include <stdio.h>
#include "alias.h"
#include "xfile.h"
struct vector;
struct matrix;

/**
  Class node:
   
  It contains all necessary informations which are not collected in the class gnode from GEFEL.
   
  Created by JK,
*/

class node
{
 public:
  node (void);
  ~node (void);
  void read (XFILE *in);
  void print (FILE *out);
  void alloc (long ncomp,long ncompo,long nlc);
  void alloc_strain (long ncomp,long nlc);
  void alloc_stress (long ncomp,long nlc);
  void alloc_other (long ncompo);
  void realloc (long ncompo,long nlc);

  void storestrain (long lcid,long fi,vector &eps);
  void storestrain (long lcid,long fi,double vol,vector &eps);
  void storestrain (long lcid,long fi,long ncomp,vector &eps);
  void storestress (long lcid,long fi,vector &sig);
  void storestress (long lcid,long fi,double vol,vector &sig);
  void storestress (long lcid,long fi,long ncomp,vector &sig);
  void storeother  (long fi,long ncomp,vector &otherv);
  void storeother (long fi,long ncomp,double vol,vector &otherv);
  void strain_averageval (long lcid);
  void stress_averageval (long lcid);
  void other_averageval  ();
  void nullstrain (long lcid);
  void nullstress (long lcid);
  double giveother(long compid);
  void nullother  ();
  void clean (long nlc);
  void alloc_meaning (long nid);
  void alloc_growstr (long nid);
  void give_transfmat(long ndofn, matrix &tmat);
  
  ///  type of cross section
  crsectype crst;
  ///  number of appropriate cross section type
  long idcs;
  ///  indicator of special coordinate system in node
  long transf;
  ///  base vectors of local coordinate system
  double *e1,*e2,*e3;
  /// principal values of strains
  double *pstra;
  /// principal values of stresses
  double *pstre;
  ///  presence of reactions
  long react;
  ///  reactions
  double *r;
  ///  number of components of strain/stress array
  long ncompstr;
  ///  number of components of other array
  long ncompother;
  ///  number of contributions to the node
  long *ncontr_strain,*ncontr_stress, ncontr_other;
  ///  volumes appropriate to contributions to the node
  double *vol_strain, *vol_stress, vol_other;
  ///  array containing strains
  double *strain;
  ///  array containing stresses
  double *stress;
  ///  array containing other values
  ///  it serves only for temporary purposes,
  ///  it is often rewritten
  double *other;
  
  ///  array of nodal values
  double *nodval;
  

  ///  meaning of DOFs
  ///  meaning=1 - displacement in x direction
  ///  meaning=2 - displacement in y direction
  ///  meaning=3 - displacement in z direction
  ///  meaning=4 - rotation about x direction
  ///  meaning=5 - rotation about y direction
  ///  meaning=6 - rotation about z direction
  long *meaning;

};

/*
  v kazdem uzlu je treba vzhledem k obecnosti softwaru
  uchovavat ndofn, protoze existuji prvky s rozdilnymi
  pocty DOF v ruznych uzlech, napr. trojuhelnikova deska
  s polynomem 5. stupne ma ve vrcholech 6 neznamych, ale
  v uzlech uprostred stran ma pouze jeden stupen volnosti

  podpory jsou tim padem dany ve vstupnim souboru, protoze
  kazdy uzel obsahuje pole cn, kde 0 znamena braneno,
  1 znamena volno a zaporna cisla jsou po porovedeni absolutni
  hodnoty pointry na predepsana posunuti
  
*/

#endif
