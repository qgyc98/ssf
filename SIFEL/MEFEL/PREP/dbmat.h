#ifndef DBMAT_H
#define DBMAT_H

#include <stdio.h>
#include "alias.h"
#include "iotools.h"



/**
   Material class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   Modified 04.2008 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
class matr
{
  public :
   mattype type; ///< material type number
   long ninst;   ///< number of indeces
   char **inst;  ///< array of strings with material properties (material indeces)
   long *instu;  ///< array indicators of indeces which are used
   long *ridx;   ///< array with renumbered used indeces
   long ninstu;  ///< counter used indeces this material type

   matr(void);
   ~matr(void);

   /// allocates memory for material paramaters stored as strings
   void alloc(long ni);

   /// allocates memory for material paramaters read/stored by mechmat procedures
   void allocmm(long ni);

   /// reads material parameters from the opened text file as strings
   long read(XFILE *in);
};

class mechmat;
struct descrip;

/**
   Material database class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
class dbmat
{
  public :
   long numt;    ///< number of material entries
   matr *mat;    ///< array with materials
   long *matu;   ///< array with indicators of material usage in the given problem
   long nmatu;   ///< number of used materials in the given problem

   dbmat(void);
   ~dbmat(void);

   /// reads materials from the text file as strings
   long read(XFILE *in);

   /// reads materials from the text file via mechmat procedures
   long readmm(XFILE *in, mechmat *mechm, descrip *d);

   /// searches material database for the given material type and index
   long search_mat(mattype ce, long ci);

   /// marks used materials
   void mark_used(long ic, long csti);

   /// renumbers instance indeces of particular material types
   void renumber_id(void);
};

#endif
