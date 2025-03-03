#ifndef MECHBCLC_H
#define MECHBCLC_H

#include "iotools.h"
#include "alias.h"
#include "galias.h"

class inicd;
class loadcase;
class dloadcase;
class gfunct;
class axisrotrec;
class ipmap;
struct vector;

/**
  Class mechbclc:
   
  It is one of the 5 most important classes of the program.
  (probdesc, mechtop, mechmat, mechbclc, mechcrsec)
   
  Class mechbclc contains data about boundary conditions, load cases.
   
  Created by JK, TKo
*/
class mechbclc
{
 public:
  mechbclc ();
  ~mechbclc ();
  void read (XFILE *in);
  void print (FILE *out);
  void read_eigenstrains (XFILE *in);
  void print_eigenstrains (FILE *out);
  void eigstrain_computation (double time);
  void aip_eigstrain_computation (long n, ipmap *ipm, double time);
  long readinic (XFILE *in);
  long printinic (FILE *out);
  void inicipval(void);
  void read_rotinidispl(XFILE *in);
  void clear_rotinidispl(long *ifn, double time);
  void store_rotinidispl(long *ifn, double time, double *r);
  long num_dofs_rotinidispl(long *ifn, double time);
  void actualize_displ_def_node(long lcid, double time);
  void apply_rotinidispl(long lcid, long *ifn, double time);
  void clear_inidispl(long *ifn);
  void store_inidispl(long *ifn, double *r);
  long num_dofs_inidispl(long *ifn);
  void apply_inidispl(long lcid, long *ifn);
  void alloc_sumcomp ();
  void comp_sum (double *rhs);
  void give_comp_sum (double *sum);
  void comp_sum_react();
  void comp_sum_pdreact();

  /// computes temperature strains at auxiliary integration points
  void aip_temperstrains(long lcid, long n, ipmap *ipm);

  /// returns the number of prescribed macro-stress components
  long give_num_mstress_comp(long lcid);

  /// returns the number of prescribed macro-strain components
  long give_num_mstrain_comp(long lcid);

  /// returns pointer to the array of types of prescribed macro-value components
  strastre* give_mstrastre(long lcid);

  /// returns pointer to the array of code (DOF) numbers of macro-stress components
  long* give_mstress_cn(long lcid);

  /// the function returns actual values of prescribed macro-strain components in the argument mstra
  void give_mstrains(long lcid, double time, vector &mstra);
  
  /// the function returns actual values of prescribed macro-stress components in the argument mstre
  void give_mstresses(long lcid, double time, vector &mstre);

  /// number of load cases
  long nlc;
  /// number of initcond
  long nico;
  /// the number of axis for calculation prescribed initial displacements
  long naxis;
  
  ///  number of components of array sumcomp
  long ncsum;
  ///  array containing sums of components of load %vector in particular directions 
  double *sumcomp;
  ///  array containing sums of components of reactions at supports in particular directions 
  double *reactsumcomp;
  ///  array containing sums of components of reactions at prescribed displacements in particular directions 
  double *pd_reactsumcomp;
  
  /// array of static load cases
  loadcase  *lc;
  /// array of time dependent load cases
  dloadcase *dlc;
  /// array of initial conditions
  inicd     *ico;
  /// array of records with definition of axis and rotations for calculation prescribed initial displacements
  axisrotrec *arotrec;

  ///  the number of general functions describing the eigenstrains
  long ngfes;
  ///  array of general functions describing the eigenstrains
  gfunct *eigstrfun;
  ///  stress-strain state
  strastrestate ssst;
};

#endif
