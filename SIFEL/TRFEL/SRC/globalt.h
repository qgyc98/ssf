#ifndef GLOBALT_H
#define GLOBALT_H

#include "genfile.h"
#include "probdesct.h"
#include "transtop.h"
#include "transmat.h"
#include "transcrsec.h"
#include "transbclc.h"
#include "lhsrhst.h"
#include "linbart.h"
#include "linbart3d.h"
#include "linbartax.h"
#include "quadbart.h"
#include "quadbartax.h"
#include "trlineart.h"
#include "trlinaxisym.h"
#include "quadlineart.h"
#include "quadquadrilatt.h"
#include "quadquadrilattax.h"
#include "quadlinaxisym.h"
#include "interfacequadrilat.h"
#include "lintett.h"
#include "linhext.h"
#include "quadhext.h"
#include "linwedget.h"
#include "gen2delem.h"
#include "outdrivert.h"
#include "stochdrivert.h"
#include "adaptivityt.h"


//**************************************
//**************************************
//  GLOBAL CONSTANTS AND VARIABLES
//**************************************
//**************************************

#ifndef EXTERN
#define EXTERN extern
#endif

//  number of DOF in the problem
EXTERN long Ndoft;

//  print detail messages
EXTERN long Mesprt;

EXTERN FILE *Outt;
EXTERN FILE *Outt1;
EXTERN FILE *Outt2;

//*****************************
//*****************************
//  CLASSES CONTAINING DATA
//*****************************
//*****************************
//  problem description
EXTERN probdesct *Tp;

EXTERN gtopology *Gtt;

//  especially topological data
EXTERN transtop *Tt;

//  especially material data
EXTERN transmat *Tm;

//  especially cross-section data
EXTERN transcrsec *Tc;

//  especially boundary condition and load case data
EXTERN transbclc *Tb;

//  right and left side of the system
EXTERN lhsrhst *Lsrst;

//  adaptivity
EXTERN adaptivityt *Adat;

//  stochastic driver
EXTERN stochdrivert *Stt;

// output driver
EXTERN outdrivert *Outdt;

//*************************
//*************************
//  MATRIX STORAGES
//*************************
//*************************

EXTERN gmatrix *Kmat,*Cmat, *Jmat, *Bmat;

//*************************
//*************************
//  FINITE ELEMENTS
//*************************
//*************************

//  1D linear bar element
EXTERN linbart *Lbt;

//  3D linear bar element
EXTERN linbart3d *Lbt3d;

//  1D linear bar element axisymmetric
EXTERN linbartax *Lbat;

//  1D quadratic bar element
EXTERN quadbart *Qbt;

//  1D quadratic axisymmetric bar element
EXTERN quadbartax *Qbat;

//  2D triangular element
EXTERN trlineart *Ltt;

//  axisymmetric triangular element
EXTERN trlinaxisym *Ltat;

//  2D quadrilateral element with linear approximation functions
EXTERN quadlineart *Lqt;

//  2D quadrilateral element with quadratic approximation functions
EXTERN quadquadrilatt *Qqt;

//  2D quadrilateral element with quadratic approximation functions for axisymmetric problems
EXTERN quadquadrilattax *Qqat;

//  axisymmetric  quadrilateral element
EXTERN quadlinaxisym *Lqat;

//  quadrilateral contact element
EXTERN interfacequadrilat *Ifcquadt;

//  3D tetrahedral element
EXTERN lintett *Ltett;

//  3D hexahedral element with linear approximation functions
EXTERN linhext *Lht;

//  3D hexahedral element with quadratic approximation functions
EXTERN quadhext *Qht;

//  3D wedge element with linear approximation functions
EXTERN linwedget *Lwt;

// general 2D element used in HERMES
EXTERN gen2delem *G2d;

/// function initializes all global variables to null values
void initnull_globt (void);

/// deletes all allocated global variables and initializes all global variables to null values
void delete_globt (void);

#endif
