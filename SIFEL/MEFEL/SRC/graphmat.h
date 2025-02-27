#ifndef GRAPHMAT_H
#define GRAPHMAT_H

#include "iotools.h"
#include "alias.h"

struct matrix;
class Equation;
class tablefunct;


/**
  The class defines special type of material model which is
  used for the spring support with nonlinear stiffness. The
  constitutive equation is described with force/displacement
  diagram. The diagram can be given by the table or by the
  expression written to the string which is parsed. The several
  equations for the different ranges can be given.
*/
class graphmat
{
  public :
   graphmat(void);
   ~graphmat(void);
   long read (XFILE *in);
   long print (FILE *out);
   void matstiff (matrix &d,long ipp);
   void nlstresses (long ipp);

   /// type of graph
   graphtype gt;
   /// number of function
   long numf;
   /// stiffness for linear elastic material
   double k;
   /// array of parsed expression strings
   char **func;
   /// array of functions of graph (displacement-force graph)
   Equation **eq;
   /// array of functions of graph (displacement-stiffness graph)
   Equation **deq;
   /// array of limit values for each function
   double   *limval;
   /// table of values of graph
   tablefunct *tab;
};

#endif
