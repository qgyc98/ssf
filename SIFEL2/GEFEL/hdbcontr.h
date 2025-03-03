#ifndef HDBCONTR_H
#define HDBCONTR_H
#include <stdio.h>
#include "iotools.h"
#include "selection.h"

#ifndef FNAMELEN
 #define FNAMELEN 1001
#endif



/**
  The class stores the setup of backup for particular 
  SIFEL modules. 
*/
class hdbcontr
{
 public:
  /// type of backup on harddisk
  hdbackuptype hdbtype;
  /// format of backup file for restoring
  hdbackupfmttype hdbfmtr;
  /// format of backup file for saving
  hdbackupfmttype hdbfmts;
  /// flag for removing old previous backup files (rmold=1 - i.e. only one set of backup file will be hold)
  answertype rmold;
  /// id of previous backup files 
  long rmold_id;
  /// precision of real numbers in output files
  long prec;
  /// backup filename for restoring
  char hdbnamer[FNAMELEN];
  /// backup filename for saving
  char hdbnames[FNAMELEN];
  sel selelemr;
  sel selelems;
  /// selection of other components for saving
  sel *selother_s;
  /// selection of other components for restoring
  sel *selother_r;
  /// array of indices of starting positions for restoring of other array
  long **selother_id;
  
  hdbcontr();
  ~hdbcontr();
 
  /// reads steup from file
  void read(XFILE *in);

  /// prints setup to file
  void print(FILE *out);

  /// returns whether restorage is required
  long restore_stat(void);

  /// returns whether saving is required
  long save_stat(void);
};

#endif
