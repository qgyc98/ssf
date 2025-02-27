#ifndef OUTDIAGM_H
#define OUTDIAGM_H
#include <stdio.h>
#include "galias.h"
#include "xfile.h"
#include "alias.h"
#include "selection.h"




/**
  The class manages output of diagrams to the text file.
  It is used for diagram output in the nonlinear or time dependent 
  problems. Selected quantities are printed at each seleted time step.
  
  Created by Tomas Koudelka,
*/
class outdiagm
{
  public :
  /// constructor
  outdiagm();
  /// destructor
  ~outdiagm();
  /// reads data from the input file
  long read(XFILE *in);
  /// prints data to the input file
  long print(FILE *out);
  /// prints header to the output diagram file
  long print_header(FILE *out);
  /// prints diagram to output diagram file
  long printval(FILE *out, long lcid, double lambda, long istep, double *fi);
  /// prints diagram to output diagram file
  long printval(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr);
  /// forced print of diagram to output diagram file
  long printval_forced(FILE *out, long lcid, double lambda, long istep, double *fi);
  /// forced print of diagram to output diagram file
  long printval_forced(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr);
  /// prints displacements
  long print_displacements(FILE *out, long lcid, long idp);
  /// prints strains
  long print_strains(FILE *out, long lcid, long idp);
  /// prints stresses
  long print_stresses(FILE *out, long lcid, long idp);
  /// prints macro-strain
  long print_macrostrain(FILE *out, long lcid, long idp);
  /// prints macro-stress
  long print_macrostress(FILE *out, long lcid, long idp);
  /// prints forces
  long print_forces(FILE *out, long idp, double *fi);
  /// prints reactions
  long print_reactions(FILE *out, long lcid, long idp);
  /// prints others
  long print_others(FILE *out, long lcid, long idp);

  /// number of printed unknowns
  long npun;
  /// selection of printed steps
  sel dstep;
  /// pid is node or ipp flag
  nodip *nif;
  /// point id
  long *pid;
  /// element id
  long *eid;
  /// coordinates of selected points 
  double *x, *y, *z;
  /// order of integration point on element
  long *ipeid;
  /// array of type of printed unknowns
  prunk *pu;
  /// indeces of pu components
  long *ipu;
};

#endif
