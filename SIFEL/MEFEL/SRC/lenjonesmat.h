#ifndef LENJONESMAT_H
#define LENJONESMAT_H

#include "iotools.h"
#include "alias.h"
struct matrix;
struct vector;
struct atsel;

/**
   class lenjonesmat defines Lennard-Jones 6-12 interatomic potential
   
   JK, 19.6.2005
*/

class lenjonesmat
{
 public:
  lenjonesmat (void);
  ~lenjonesmat (void);
  void read (XFILE *in);
  double compute_force (double r);
  double equilib_distance ();
  
  double first_derivative (double r);
  double second_derivative (double r);
  
  ///  constants of potential
  double eps;
  double sig;
  
  ///  minimum accaptable distance between atoms
  double mindist;
  ///  equilibrium distance between two atoms and for given constants of potential
  double eqdist;

};

#endif
