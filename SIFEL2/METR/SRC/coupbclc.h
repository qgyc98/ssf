#ifndef COUPBCLC_H
#define COUPBCLC_H

#include <stdio.h>
#include "iotools.h"
#include "alias.h"
class inicd;
class gfunct;
//class loadcase;
class dloadcase;


class coupbclc
{
 public:
  coupbclc (void);
  ~coupbclc (void);
  void read (FILE *in);
  long readinic (XFILE *in);

  void read_eigenstrains (XFILE *in);
  void eigstrain_computation (double time);
  
  /// number of load cases
  long nlc;
  /// number of initcond
  long nico;

  //loadcase  *lc;
  dloadcase *dlc;
  inicd     *ico;

  ///  the number of general functions describing the eigenstrains
  long ngfes;
  ///  array of general functions describing the eigenstrains
  gfunct *eigstrfun;
  ///  stress-strain state
  strastrestate strastre;

};

#endif
