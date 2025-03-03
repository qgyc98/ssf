#ifndef DBMATT_H
#define DBMATT_H

#include <stdio.h>
#include "aliast.h"
#include "iotools.h"



/**
   Material class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   Modified 04.2008 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
class matrt
{
  public :
   mattypet type; ///< material type number
   long ninst;   ///< number of indeces
   char **inst;  ///< array of strings with material properties (material indeces)
   long *instu;  ///< array indicators of indeces which are used
   long *ridx;   ///< array with renumbered used indeces
   long ninstu;  ///< counter used indeces this material type

   matrt(void);
   ~matrt(void);

   /// allocates memory for material paramaters stored as strings
   void alloc(long ni);

   /// allocates memory for material paramaters read/stored by transmat procedures
   void alloctm(long ni);

   /// reads material parameters from the opened text file as strings
   long read(XFILE *in);
};

class transmat;
struct descript;


/**
   Material database class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
class dbmatt
{
  public :
   long   numt;   ///< number of material entries
   matrt *mat;    ///< array with materials
   long  *matu;   ///< array with indicators of material usage in the given problem
   long   nmatu;  ///< number of used materials in the given problem

   dbmatt(void);
   ~dbmatt(void);

   /// reads materials from the text file as strings
   long read(XFILE *in);

   /// reads materials from the text file via mechmat procedures
   long readtm(XFILE *in, transmat *mechm, descript *d);

   /// searches for given material type and index
   long search_mat(mattypet ce, long ci);

   /// marks used materials
   void mark_used(long ic, long csti);

   /// renumbers instance indeces of particular material types
   void renumber_id(void);
};

#endif
