#ifndef PGLOBAL_H
#define PGLOBAL_H

#ifndef EXTERN
#define EXTERN extern
#endif

#include "pprobdesc.h"
#include "genfile.h"
#include "seqfilesm.h"
#include "pglobalg.h"


//*****************************
//*****************************
//  CLASSES CONTAINING DATA
//*****************************
//*****************************
//  problem description
EXTERN pprobdesc *Pmp;

//  parallel solver
EXTERN psolver *Psolm;

/// initializes all parallel global variables to null values
void pinitnull_glob(void);

/// deallocates memory of all parallel global variables and set them to null values
void pdelete_glob(void);
#endif
