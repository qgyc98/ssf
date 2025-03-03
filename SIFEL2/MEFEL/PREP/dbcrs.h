#ifndef DBCRS_H
#define DBCRS_H

#include <stdio.h>
#include "iotools.h"
#include "alias.h"

class csec
/**
   Cross-section class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
{
  public :
   crsectype type;
   long ninst;   ///< number of indeces
   char **inst;  ///< array of strings with cross-section properties (cross-section indeces)
   long *instu;  ///< array indicators of indeces which are used
   long *ridx;   ///< array with renumbered used indeces
   long ninstu;  ///< counter used indeces this cross-section type

   csec(void);
   ~csec(void);

   /// alllocates memory for storage of cross section parameters as strings
   void alloc(long ni);

   /// alllocates memory for storage of cross section parameters in mechcrsec
   void allocmc(long ni);

   /// reading cross-section
   long read(XFILE *in);
};

class mechcrsec;
struct descrip;

class dbcrs
/**
   Cross-section database class

   Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
{
  public :
   long numt;  ///< number of cross-section entries
   csec *crs;  ///< array with cross-sections
   long *crsu; ///< array with indicators of cross-section usage in the given problem
   long ncrsu; ///< number of cross-sections used in the given problem

   dbcrs(void);
   ~dbcrs(void);

   /// reads cross section parameters from the opened text file  as strings
   long read(XFILE *in);

   /// reads cross section parameters from the opened text file with help of mechcrsec procedures
   long readmc(XFILE *in, mechcrsec *mechc, descrip *d);

   /// searching cross-section with given type and index
   long search_crs(crsectype ce, long ci);

   /// marks ic-th corss-section type and its csti-th parameter set as used
   void mark_used(long ic, long csti);

   /// renumbers database indeces according with used cross-sections
   void renumber_id(void);
};

#endif
