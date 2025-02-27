#ifndef PGLOBALT_H
#define PGLOBALT_H

#ifndef EXTERN
#define EXTERN extern
#endif

#include "pprobdesct.h"
#include "genfile.h"
#include "seqfilest.h"

// Myrank, Ndom and Nproc defined in the below header file due to PMETR compilation conflicts
#include "pglobalg.h"


//*****************************
//*****************************
//  CLASSES CONTAINING DATA
//*****************************
//*****************************

//  problem description
EXTERN pprobdesct *Ptp;

//  parallel solver
EXTERN psolver *Psolt;

/// initializes all parallel global variables to null values
void pinitnull_globt(void);

/// deallocates memory of all parallel global variables and set them to null values
void pdelete_globt(void);
#endif
