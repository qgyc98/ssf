#ifndef PGLOBAL_H
#define PGLOBAL_H

#ifndef EXTERN
#define EXTERN extern
#endif

#include "genfile.h"
#include "seqfilesm.h"



//**************************************
//**************************************
//  GLOBAL CONSTANTS AND VARIABLES
//**************************************
//**************************************

//  number of processors
EXTERN int Nproc;

//  rank of process
EXTERN int Myrank;

//  number of subdomain
EXTERN int Ndom;

#endif
