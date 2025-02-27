#ifndef HARDSOFT_H
#define HARDSOFT_H

#include "iotools.h"
#include "alias.h"
struct matrix;
struct vector;
struct atsel;

/**
   class defines type of hardening/softening
   
   JK, 8.8.2005
*/

class hardsoft
{
 public:
  hardsoft ();
  ~hardsoft ();

  void read (XFILE *in);
  void print (FILE *out);

  void hvalues (vector &sigt,vector &dgds,vector &h);
  void dhdsigma (vector &sigt,vector &dgds,matrix &dgdsds,vector &dhds);
  void dhdqpar (vector &sigt,vector &dgds,matrix &dgdsdq,vector &dhdq);
  void dhdgamma (vector &dhdg);
  
  void changeparam (atsel &atm,vector &val);

  
  ///  type of hardening/softening
  hardensoften ths;
  ///  input variable to hardening/softening
  hsinputvar ivhs;
  
  ///  limit for norm of plastic strain
  double epspu;
  ///  limit for norm of plastic strain in tension
  double epsput;
  ///  limit for norm of plastic strain in compression
  double epspuc;

  ///  computer zero
  double zero;
  
};

#endif
