#ifndef DLOADCASE_H
#define DLOADCASE_H

#include "xfile.h"
#include "alias.h"
#include "seismtool.h"
#include <stdio.h>

struct vector;
class dloadn;
class dloadel;
class dloadpd;
class loadcase;
class ipmap;

/**
  Class dloadcase defines load cases in dynamic problems 
  or time dependent problems.
   
  JK, 27.9.2004
*/

class dloadcase
{
 public:
  dloadcase (void);
  ~dloadcase (void);
  /// reads definition of load case from the opened text file
  void read (XFILE *in);
  
  /// check consistency of prescribed macro-stress/strain components at subload cases
  long check_consistency_macrostr();
  
  /// prints definition of load case to the opened text file
  void print (FILE *out);
  
  /// assembles right hand side of the given load case, all types of load are taken into account
  void assemble (long lcid,double *rhs,double *flv,long n, double t);
  
  /// assembles right hand side of the given load case, all types of load are taken into account
  void assemble (double *rhs, double *flv,double *lhs);
  
  /// compute reactions for the given load case
  void compute_reactions (long lcid);
  
  /// reallocates array of inital values of prescribed displacements at gradual construction process problems
  void realloc_pid_array();
  
  /// computes the maximum number of different values of prescribed displacements
  void compute_max_npd();
  
  /// returns total value of prescribed displacement for the given DOF and time
  double get_pd(double time, long dof);
  
  /// returns total value of macro-strain component at the given time
  double get_macrostre(double time, long dof);
  
  /// returns total value of macro-strain %vector at the given time
  void get_macrostre(double time, vector &meps);

  /// returns total value of macro-stress component at the given time
  double get_macrostra(double time, long dof);
  
  /// returns total value of macro-stress %vector at the given time
  void get_macrostra(double time, vector &msig);

  /// assembles right hand side of the given load case, only temperature loads are taken into account
  void tempercontrib (long lcid,double *rhs,long n,double t);
  
  /// compute temperature strains at auxiliary integration points for the given load case
  void aip_temperstrains(long lcid, long n, ipmap *ipm);
  
  /// returns the number of prescribed macro-stress components
  long give_num_mstress_comp();
  
  /// returns the number of prescribed macro-strain components
  long give_num_mstrain_comp();
  
  /// returns the pointer to the array of types of prescribed macro-value components
  strastre* give_mstrastre();

  /// returns pointer to the array of code (DOF) numbers of macro-stress components
  long* give_mstress_cn();

  /// returns pointer to the prescribed macro-strain array of the given load case
  double* give_mstrains();

  /// the function returns actual values of prescribed macro-strain components in the argument mstra
  void give_mstrains(double time, vector &mstra);

  /// the function returns actual values of prescribed macro-stress components in the argument mstre
  void give_mstresses(double time, vector &mstre);
  
  
  ///  type of dynamic load
  dynload tdl;

  ///  the number of subload cases
  long nslc;  

  ///  the number of time dependent loads at nodes
  long nln;

  ///  the number of time dependent loads defined on elements 
  long nle;

  ///  the number of time dependent prescribed displacements
  long npd;

  /// the maximum nuber of prescribed displacements at all subloadcases (used in gradual construction problems)
  long max_npd;

  /// the number of prescribed initial displacements at initial conditions
  long npid;

  /// the number of prescribed initial displacements by rotations
  long nrpid;

  /// the number of all prescribed initial displacements: tnpid = npid+nrpid, it is dimension of array pid
  long tnpid;
  
  /// array of time dependent loads defined at nodes
  dloadn  *lon;

  /// array of time dependent loads defined on elements
  dloadel *loe;

  /// array of time dependent prescribed displacements defined at nodes
  dloadpd *pd;

  /// array of prescribed initial displacements at nodes for gradual contruction problem
  double  *pid; 
  
  /// array of subload cases
  loadcase *slc;

  /// array of time functions
  gfunct *gf;

  ///  seismic tool
  seismtool stool;
};

#endif
