#ifndef HPGLOBAL_H
#define HPGLOBAL_H

#ifndef EXTERN
#define EXTERN extern
#endif

#include "mpi.h"
#include "genfile.h"
#include "seqfilesm.h"
#include "seqfilest.h"
#include "global.h"
#include "globalt.h"


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

//  the number of tile types
EXTERN long Ntiletypes;

//  the number of tiles in tiling
EXTERN long Ntiles;

//  the number of internal DOFs on tiles
//  it contains Ntiletypes components
EXTERN long* Nidof;

//  the number of interface DOFs on tiles
//  it contains Ntiletypes components
EXTERN long* Nbdof;

//  the number of boundary/interface DOFs on tiles in tilings
//  it contains Ntiles components
EXTERN long* Nbdoftiling;

//  the code numbers of Schur complements on tiles
//  it contains Ntiles rows and Nbdof[i] columns
EXTERN long** Cnbn;

//  the number of tilings
EXTERN long Ncomb;

//  the list of tiles in tilings
//  it contains Ncomb rows and Ntiles columns
EXTERN long** Tiling;


//  type of problem on the MEFEL master processor
//  it is needed also on slave processors
EXTERN problemtype Mtprobm;

//  type of problem on the TRFEL master processor
//  it is needed also on slave processors
EXTERN problemtypet Mtprob;

// name of processors
EXTERN char proc_namet[MPI_MAX_PROCESSOR_NAME];
EXTERN int nameLengtht;

//  general matrices for stiffness matrices of tiles
//  general matrices for conductivity matrices of tiles
EXTERN gmatrix *SSmat;

//  arrays for right hand sides of tiles
EXTERN double **Srhs;

//  arrays for left hand sides of tiles
EXTERN double **Slhs;

//  array for prescribed DOFs on tiling
EXTERN double *Pdofval;

#endif
