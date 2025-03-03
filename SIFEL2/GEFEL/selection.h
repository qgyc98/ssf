#ifndef SELECTION_H
#define SELECTION_H

#include <stdio.h>
#include "galias.h"
#include "xfile.h"
#include "timecontr.h"
#include "siftop.h"



/**
  Class manages the selections of steps, nodes, elements, strains, stresses and
  other values. It is devoted especially to outdrivers, but generally 
  it should be used in whatever case where the above quantities are need to select.
  
  created 03.2004 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz 
  modified 10.10.2007 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz 
*/
class sel
{
  public:
   /// constructor
   sel(); 
   /// destructor
   ~sel(); 
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// copies given instance to the destination instance
   void copy(sel &dest);
   /// checks presence of id at given selection
   long presence_id(double id);
   /// checks presence of id at given selection
   long presence_id(long id);
   /// checks presence of ids at given selection
   long presence_id(long id, double rid);
   /// checks presence of ids at given selection
   long presence_id(long id, double rid, timecontr &tc);
   /// checks presence of id at given selection
   long presence_id(double id, long &ir);
   /// checks presence of id at given selection
   long presence_id(long id, long &ir);
   /// checks presence of id at given selection with index greater than ir
   long presence_idgt(long id, long &ir);
   /// checks presence of id at given selection with index sid
   long presence_id(long id, long sid, long &ir);
   /// checks presence of ids at given selection
   long presence_id(long id, double rid, long &ir);
   /// checks presence of ids at given selection and conjugated selections
   long presence_id(sel *consel, long id, long idc);
   /// checks presence of ids at given selection and conjugated selections
   long presence_id(sel *consel, long id, long idc, long &ir);
   /// counts number of selected components in array with length tncomp
   long give_nselcomp(long tncomp);
   /// returns number of selected objects
   long give_num_lst_items(long nobj, long *selobj);
   /// returns number of detected ranges in the array selobj
   long give_num_range_items(long nobj, long *selobj);
   /// convert the given list of objects to the range selection type
   void conv2range(long nit, long nobj, long **selobj, sel *conselent, sel *conselr, long *transent, long *transr);
   /// convert the given list of objects to the list selection type
   void conv2lst(long nit, long nobj, long **selobj, sel *conselent, sel *consell, long *transent, long *transl);
   /// convert the given list of objects to the range selection type
   void conv2range(long nit, long nobj, long **selobj);
   /// convert the given list of objects to the list selection type
   void conv2lst(long nit, long nobj, long **selobj);
   /// converts selection type of property to selection type of list
   long conv_selprop(siftop *top, objtype ot, sel *conselent, sel *&consel, long *transent, long *&trans);
   /// converts selection type of property to selection type of list
   long conv_selprop(const siftop *top, objtype ot);
   /// converts list of the selected objects to the range type of selection
   void conv2range(long nit, long nobj, ivector &selobj);
   /// converts list of the selected objects to the list type of selection
   void conv2lst(long nit, long nobj, ivector &selobj);

   /// type of selection
   seltype st;
   /// number of selections   
   long n;
   /// array with the first indices
   long *id1;
   /// array with lengthes of selections
   long *ncomp;
   /// array with begining real values of real range, real list or real period
   double *rid1;
   /// array with ending real values of real range
   double *rid2;
   /// array of selections
   sel    *sarray;
   /// required error of real range, real list or real period
   double err;
   /// initial time for selection of real period
   double initime;
   /// final time for selection of real period
   double fintime;
   /// real period
   double r;
   /// type of entity for property selection of nodes/elements 
   gentity *ent;
};

#endif
