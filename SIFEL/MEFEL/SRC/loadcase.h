#ifndef LOADCASE_H
#define LOADCASE_H

#include "xfile.h"
#include "galias.h"
#include "alias.h"
#include <stdio.h>

#ifndef FNAMELEN
 #define FNAMELEN 1001
#endif

class loadn;
class loadel;
struct ivector;
struct vector;
class ipmap;

/**
  Class loadcase defines load cases in static problems.
   
  Created by JK,
*/

class loadcase
{
 public:
  loadcase (void);
  ~loadcase (void);

  /// reads definition of load case from the opened text file
  void read (XFILE *in);

  /// prints definition of load case to the opened text file
  void print (FILE *out);

  /// assembles right hand side of the given load case, all types of load are taken into account
  void assemble (long lcid,double *rhs,double *flv, double scale);

  /// assembles right hand side of the given load case, only force and temperature loads are taken into account (without prescribed displacements)
  void assemblewopd (long lcid,double *rhs,double *flv,double scale);

  /// assembles right hand side of the given load case, only force load contributions are taken into account
  void forcecontrib (long lcid,double *rhs,double scale);

  /// assembles right hand side of the given load case, only prescribed displacement contributions are taken into account
  void prdisplcontrib (long lcid, double *rhs);

  /// assembles right hand side of the given load case, only temperature loads are taken into account
  void tempercontrib (long lcid,double *rhs,double scale);

  /// accumulate temperature strains at auxiliary integration points for the given load case
  void aip_cumultemperstrains_contrib (long lcid, double scale, long n, ipmap *ipm);

  /// compute reactions for the given load case
  void compute_reactions (long lcid, double scale, answertype comp_intfc);

  /// returns the number of prescribed macro-stress components
  long give_num_mstress_comp();

  /// returns the number of prescribed macro-strain components
  long give_num_mstrain_comp();

  /// returns pointer to the array of types of precribed macro-value components
  strastre* give_mstrastre();

  /// returns pointer to the array of code (DOF) numbers of macro-stress components
  long* give_mstress_cn();

  /// returns prescribed macro-strain components in the argument mstra
  void give_mstrains(vector &mstra);
  
  /// returns prescribed macro-stress components in the argument mstre
  void give_mstresses(vector &mstre);

  ///  number of loaded nodes
  long nln;

  ///  number of loaded elements
  long nle;

  ///  number of prescribed displacements
  long npd;

  ///  number of prescribed temperature changes
  long npt;
  
  ///  type of temperature changes
  long tempchang;

  ///  loaded nodes
  loadn *lon;

  ///  loaded elements
  loadel *loe;

  ///  prescribed displacements
  double *pd;

  ///  prescribed temperature changes
  double *pt;
  
  /** 
    the total number of prescribed macro-values (strains/stresses) 
    in homogenization problems (Mp->homog > 2), 
    it is the maximum number of stress/strain components from all elements and
    it also represents dimension of the below arrays
  */
  long ncompmv;
  
  /** 
    array of component types of prescribed macro-values in homogenization problems (Mp->homog > 2), mstrastre[ncompmv]
    mstrastre[i] = strain => the i-th macro-strain value is being prescribed as i-th component,
    mstrastre[i] = stress => the i-th macro-stress value is being prescribed as i-th component
  */
  strastre *mstrastre;

  /// array of macro stresses values in homogenization problems (Mp->homog = 3,5,7,9), mstress[ncompmv]
  double *mstress;

  /// the number of prescribed macro-stress components (Mp->homog > 2)
  long nmstrecomp;

  /** 
    array of code numbers of prescribed macro-stress components, mstress_cn[ncompmv]
    mstress_cn[i] > 0 <=> code (DOF) number of i-th macro-value component which represents
                          prescribed macro-stress at the i-th direction
    mstress_cn[i] = 0 <=> no macro-stress is being prescribed in the i-th direction
  */
  long *mstress_cn;

  /// array of macro strain values in homogenization problems (Mp->homog = 4,6,8,9), mstrain[ncompmv]
  double *mstrain;

  /// the number of prescribed macro-strain components (Mp->homog > 2)
  long nmstracomp;

  ///  data file for reading
  char temp_file[FNAMELEN];

  /// number of multiple data files
  long ntemp_file;

  ///  prescribed temperature changes as multiple values according to time
  double **pts;

  /// time of individual temperature settings (file)
  double *temp_time;
};

#endif
