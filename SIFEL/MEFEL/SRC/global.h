#ifndef GLOBAL_H
#define GLOBAL_H

//extern "C" { int PARDISO (void *, int *, int *, int *, int *, int *,
//		    double *, int *, int *, int *, int *, int *,
//		    int *, double *, double *, int *);
//}


//#include "probdesc.h"
//#include "mechtop.h"
//#include "mechmat.h"
//#include "mechcrsec.h"
//#include "mechbclc.h"
//#include "lhsrhs.h"
//#include "adaptivity.h"
//#include "stochdriver.h"
//#include "outdriverm.h"
//#include "flsubdom.h"

#include <stdio.h>

#ifdef EXTENDED_GLOBINC
 #include "probdesc.h"
 #include "mechtop.h"
 #include "mechmat.h"
 #include "mechcrsec.h"
 #include "mechbclc.h"
 #include "lhsrhs.h"
 #include "adaptivity.h"
 #include "stochdriver.h"
 #include "outdriverm.h"
 #include "flsubdom.h"
 #include "gtopology.h"
 #include "gmatrix.h"
 #include "elemhead.h"
 #include "vector.h"
#else
 class gtopology;
 class gmatrix;
 class probdesc;
 class mechtop;
 class mechmat;
 class mechcrsec;
 class mechbclc;
 class lhsrhs;
 class adaptivity;
 class stochdriver;
 class outdriverm;
 class flsubdom;

 class barel2d;
 class barel3d;
 class barelq2d;
 class barelq3d;
 class beamel2d;
 class beamel3d;
 class beamgen3d;
 class soilbeam;
 class beam2dspec;
 class springel;
 class planeelemlt;
 class planeelemqt;
 class planeelemrotlt;
 class planeelemlq;
 class planeelemqq;
 class planeelemrotlq;
 class planeelemsubqt; 
 class cctelem;
 class dktelem;
 class dstelem;
 class q4plate;
 class ArgyrisTriangle;
 class argyrisplate;
 class soilplatetr;
 class soilplateq;
 class axisymlq;
 class axisymlt;
 class axisymqq;
 class axisymqt;
 class axisymcq;
 class axisymlqinterface;
 class shelltr;
 class shelltrm;
 class shellq;
 class lintet;
 class quadtet;
 class linhex;
 class quadhex;
 class lintetrot;
 class linhexrot;
 class linwedge;
 class quadwedge;
 class hexinterface;
 class elemparticle;
 class plquadinterface;
 class ArgyrisTriangle;
 class argyrisplate;
 class quadrilatkirch;
 class dkq;
 class tetralattice;
 struct vector;
#endif

#ifndef EXTERN
#define EXTERN extern
#endif


//**************************************
//**************************************
//  GLOBAL CONSTANTS AND VARIABLES
//**************************************
//**************************************
//  number of DOF in the problem
EXTERN long Ndofm;

//  print detail messages
EXTERN long Mespr;

EXTERN FILE *Out;

//*****************************
//*****************************
//  CLASSES CONTAINING DATA
//*****************************
//*****************************
//  problem description
EXTERN probdesc *Mp;

//  general topology (problem independent)
EXTERN gtopology *Gtm;

//  especially topological data
EXTERN mechtop *Mt;
EXTERN mechtop *Mtt;

//  especially material data
EXTERN mechmat *Mm;
EXTERN mechmat *Mmm;

//  especially cross section data
EXTERN mechcrsec *Mc;

//  especially boundary condition and load case data
EXTERN mechbclc *Mb;

//  right and left side of the system
EXTERN lhsrhs *Lsrs;

//  adaptivity
EXTERN adaptivity *Ada;

//  stochastic driver
EXTERN stochdriver *St;

// output driver
EXTERN outdriverm *Outdm;

// floating subdomains
EXTERN flsubdom *Fsd;

//*************************
//*************************
//  MATRIX STORAGES
//*************************
//*************************


//  stiffness matrix
EXTERN gmatrix *Smat;

//  mass matrix
EXTERN gmatrix *Mmat;

//  damping matrix
EXTERN gmatrix *Dmat;

//  initial stiffness matrix
EXTERN gmatrix *Ismat;

//  array of auxiliary matrices
EXTERN gmatrix **Amat;

//*************************
//*************************
//  FINITE ELEMENTS
//*************************
//*************************

//  bar element with linear functions in 2D
EXTERN barel2d *Bar2d;

//  bar element with linear functions in 2D
EXTERN barel3d *Bar3d;

//  bar element with quadratic functions in 2D
EXTERN barelq2d *Barq2d;

//  bar element with quadratic functions in 3D
EXTERN barelq3d *Barq3d;

//  beam element in 2D
EXTERN beamel2d *Beam2d;

//  beam element in 3D
EXTERN beamel3d *Beam3d;

//  generalized beam element in 3D
EXTERN beamgen3d *Beam3dg;

//  subsoil beam element
EXTERN soilbeam *Sbeam;

//  special beam in 2D
EXTERN beam2dspec *Spbeam2d;

//  spring element
EXTERN springel *Spring;

//  linear triangular plane finite element
EXTERN planeelemlt *Pelt;

//  quadratic triangular plane finite element
EXTERN planeelemqt *Peqt;

//  linear triangular plane finite element with rotational degree of freedom
EXTERN planeelemrotlt *Perlt;

//  linear quadrilateral plane finite element
EXTERN planeelemlq *Pelq;

//  quadratic quadrilateral plane finite element
EXTERN planeelemqq *Peqq;

//  linear quadrilateral plane finite element with rotational degree of freedom
EXTERN planeelemrotlq *Perlq;

//  quadratic triangular plane finite element with subparametric geometry description
EXTERN planeelemsubqt *Pesqt;

//  constant curve triangular plate finite element
EXTERN cctelem *Cct;

//  Kirchhoff discrete plate finite element
EXTERN dktelem *Dkt;

//  thick triangular plate element
EXTERN dstelem *Dst;

//  quadrilateral discrete plate finite element
EXTERN q4plate *Q4pl;

//  triangular Argyris's element for Kirchhoff plates
EXTERN ArgyrisTriangle *Argtr;
//  triangular Argyris's element for Kirchhoff plates
EXTERN argyrisplate *Argtrpl;

//  quadrilateral plate element based on the Kirchhoff theory
EXTERN quadrilatkirch *Qkirch;

//  quadrilateral plate element based on the discrete Kirchhoff theory
EXTERN dkq *Dkqelem;


//  subsoil plate triangular element
EXTERN soilplatetr *Spltr;

//  subsoil plate triangular element
EXTERN soilplateq *Splq;

//  axisymmetric quadrilateral finite element with bilinear shape functions
EXTERN axisymlq *Asymlq;

//  axisymmetric triangular finite element with linear shape functions
EXTERN axisymlt *Asymlt;

//  axisymmetric quadrilateral finite element with quadratic shape functions
EXTERN axisymqq *Asymqq;

//  axisymmetric triangular finite element with quadratic shape functions
EXTERN axisymqt *Asymqt;

//  axisymmetric quadrilateral finite element with cubic shape functions
EXTERN axisymcq *Asymcq;

//  axisymmetric quadrilateral finite element with linear shape functions
EXTERN axisymlqinterface *Asymlqifc;

//  shell triangular element
EXTERN shelltr *Shtr;

//  shell triangular element based on Mindlin theory
EXTERN shelltrm *Shtrm;

//  shell quadrilateral element
EXTERN shellq *Shq;

//  threedimensional linear tetrahedral finite element
EXTERN lintet *Ltet;

//  threedimensional quadratic tetrahedral finite element
EXTERN quadtet *Qtet;

//  threedimensional linear hexagonal finite element
EXTERN linhex *Lhex;

//  threedimensional quadratic hexagonal finite element
EXTERN quadhex *Qhex;

//  threedimensional linear tetrahedral finite element with rotational degrees of freedom
EXTERN lintetrot *Ltetrot;

//  threedimensional linear hexahedral finite element with rotational degrees of freedom
EXTERN linhexrot *Lhexrot;

//  threedimensional linear wedge finite element
EXTERN linwedge *Lwed;

//  threedimensional quadratic wedge finite element
EXTERN quadwedge *Qwed;

EXTERN hexinterface *Hexifc;

//  element for particle computation
EXTERN elemparticle *Pelem;

//  plane quadrilateral element for interface between edges of 2D elements
EXTERN plquadinterface *Pqifc;

//  element for lattice model
EXTERN tetralattice *Tlatt;

// The total number of hypoplasticty model evaluation
EXTERN double Neval;

EXTERN double Omp_wtime;

// array of element internal force vectors Elemfrc[i][j] means j-th component of internal force vector of i-th element
EXTERN vector *Elemfrc;

// following variables are just for the debugging OpenMP
//#ifndef INC_OPENMP
//EXTERN FILE *out;
//EXTERN FILE *out2;
//EXTERN FILE *in;
//#endif
//#ifdef INC_OPENMP
//EXTERN FILE *out;
//EXTERN FILE *out2;
//#endif

/// function initializes all global variables to null values
void initnull_glob (void);

/// deletes all allocated global variables and initializes all global variables to null values
void delete_glob (void);
#endif
