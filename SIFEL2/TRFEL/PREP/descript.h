#ifndef DESCRIP_H
#define DESCRIP_H

#include "galias.h"
#include <stdio.h>


/**
  This structure holds data about input files of TRFEL preprocessor
*/
struct descript
{
  char topf[1025];   ///< topology file name
  char matf[1025];   ///< material file name
  char crf[1025];    ///< cross-section file name
  char icf[1025];    ///< initial condition file name
  char hangnf[1025]; ///< hanging nodes file name

  meshform meshfmt;  ///< format of topology file
  long paral;        ///< indicator whether sequential version of preprocessor is used (= 0) or paralell (= 1)
  long redgn;        ///< indicator whether edge numbers of element should be read in topology file (=1)
  answertype inicdf; ///< indicator of initial conditions given in a file
  answertype matsec; ///< indicator of material section in the preprocessor file
  answertype crssec; ///< indicator of cross-section section in the preprocessor file
  answertype matstr; ///< indicator whether materials are read as strings (yes) or by transmat procedures (no)
  answertype crsstr; ///< indicator whether cross section parameters are read as strings (yes) or by transcrsec procedures (no)
  answertype matkwd; ///< indicator whether the keywords are required (yes) in materials read by transmat procedures
  answertype crskwd; ///< indicator whether the keywords are required (yes) in cross sections read by transcrsec procedures

  descript();
  ~descript();

  long print(FILE *out); ///< the function prints the preprocessor description to the text file
};

#endif
