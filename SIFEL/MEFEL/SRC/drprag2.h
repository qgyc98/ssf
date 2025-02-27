#ifndef DRPRAG2_H
#define DRPRAG2_H

#include "iotools.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   This class defines Drucker-Prager plasticity model

   structure of stored data
   plastic strains, consistency parameter, hardening parametr

   28.3.2002
*/

class drprag2
{
 public:
  drprag2 (void);
  ~drprag2 (void);
  void read (XFILE *in);
  void print (FILE *out);
  void giveIota(strastrestate ssst, vector &iota);
  void giveIdev(strastrestate ssst, matrix &idev);
  void getImportantValues(long ipp, long ido, double &K_con, double &G_con, double &rho_tr, double &p_tr);
  void getImportantValues(long ipp, long ido, double &K_con, double &G_con, double &rho_tr, double &p_tr, vector &iota, vector &epsilon_etr, matrix &Idev);
  double getQvalue(double gamma, long ipp, long ido);
  void getValueOfHardening(long ipp, long ido, double &H_value, double &Hder_value);
  double newton_one_element_case2(long ipp, long ido);
  double newton_one_element_case3(long ipp, long ido);
  void matstiff (matrix &d, long ipp, long ido);
  void tangentstiff (matrix &td,long ipp,long ido);
  void nlstresses (long ipp,long im,long ido);
  //void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im,long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  double give_consparam (long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  

  ///  friction angle
  double phi;
  ///  cohesion
  double c;
  ///  dilatation
  double psi;
  /// angle of linear hardening/softening
  double theta;
  /// limit cohesion
  double clim;


  double alpha; ///< material constant alpha
  double alpha1; ///< material constant alpha1
  double beta; ///< material constant beta

  ///  stress return algorithm
  strretalg sra;
};

#endif
