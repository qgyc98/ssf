#ifndef OUTRESFILE
#define OUTRESFILE

#include "galias.h"
#include "selection.h"
#include <stdio.h>

class lcoordsys;
class outquantm;
class stochdriver;
class gnodvalvm;
struct XFILE;


/**
  The class represents description of the result ouput file
*/
class outresfilem
{
  public:
   outresfilem();
   ~outresfilem();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// printf header to the output file
   void print_header(long istep, double time);
   /// creates modified file name (default format extension, stepid and parameter set id are added optionally)
   void get_modified_fname(long stochid, long stepid);
   /// prints required quantities to the given output file
   void print_quants(long lcid, double time, long istep, lcoordsys *lcsarr, gnodvalvm &gnv, long idn1, long ide1);
   /// returns default extension of output file format
   static const char* give_fmt_dext(resfilefmt fmt);

   /// maximum length of file name strings
   static const long fnl = 1025;
   /// base name of the output file
   char outfbn[fnl];
   /// length of outfn string - it is set in outgrfile::read 
   long outfbnl;
   /// format of the given result file
   resfilefmt rffmt;
   /// selection of steps which the quantities will be printed in
   sel selstep;
   /// selection of load cases which the quantities will be printed in
   sel sellc;
   /// number of quantities to be printed out
   long nqnt;
   /// array of records with quantity output desprition
   outquantm *qnt;
   /// pointer to the associated output file of the given quantity
   FILE *outf;
   /// modified output file name which corresponds to the opened file in outf
   char outfmn[fnl];
};





#endif
