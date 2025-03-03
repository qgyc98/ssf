#ifndef GLOBALC_H
#define GLOBALC_H

#include "probdescc.h"
#include "lhsrhsc.h"
#include "global.h"
#include "globalt.h"
#include "coupmat.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "couptop.h"
#include "coupbclc.h"
#include "barelc.h"
#include "quadrilatc.h"
#include "axiquadc.h"
#include "axisymc.h"
#include "hexahedc.h"
#include "outdriverc.h"
#include "ipmap.h"

class gtopology;


#ifndef EXTERN
#define EXTERN extern
#endif


//**************************************
//**************************************
//  GLOBAL CONSTANTS AND VARIABLES
//**************************************
//**************************************
//  number of DOF in the problem
EXTERN long Ndofc;

//  print detail messages
EXTERN long Mesprc;

EXTERN FILE *Outc;

//*****************************
//*****************************
//  CLASSES CONTAINING DATA
//*****************************
//*****************************
//  problem description
EXTERN probdescc *Cp;

//  general topology (problem independent)
EXTERN gtopology *Gtu;

//  especially topological data
EXTERN couptop *Ct;

//  especially material data
EXTERN coupmat *Cm;
EXTERN coupmatu *Cmu;
EXTERN coupmatl *Cml;

// array of mapping from TRFEL integration points to MEFEL integration points
// TM_ip[ipp_trf] = ipp_mef
EXTERN long *TM_ip;

// array of mapping from MEFEL integration points to TRFEL integration points
// MT_ip[ipp_mef] = ipp_trf
EXTERN long *MT_ip;

// array of mapping from TRFEL integration points to MEFEL integration points
EXTERN ipmap *TMipmap;

// array of mapping from MEFEL integration points to TRFEL integration points
EXTERN ipmap *MTipmap;

// array of TRFEL->MEFEL correspondence map TM_nod_map[i]=MEFEL node id of TRFEL i-th node 
EXTERN long *TM_nod_map;

//  especially cross section data

//  especially boundary condition and load case data
EXTERN coupbclc *Cb;

//  right and left side of the system
EXTERN lhsrhsc *Lsrsc;


// output driver
EXTERN outdriverc *Outdc;

//*************************
//*************************
//  MATRIX STORAGES
//*************************
//*************************

//  zero-order matrix
EXTERN gmatrix *D0mat;

//  first-order matrix
EXTERN gmatrix *D1mat;

//*************************
//*************************
//  FINITE ELEMENTS
//*************************
//*************************

//  onedimensional coupling finite element
EXTERN barelc *Cbar;

//  quadrilateral coupling finite element
EXTERN quadrilatc *Cquad;

//  quadrilateral coupling finite element
EXTERN axiquadc *Caxiq;

//  quadrilateral axisymmetric element for fully coupled material
EXTERN axisymc *Caxifc;

//  hexahedral coupling finite element
EXTERN hexahedc *Chex;

/// function initializes all global variables to null values
void initnull_globc (void);

/// deletes all allocated global variables and initializes all global variables to null values
void delete_globc (void);
#endif
