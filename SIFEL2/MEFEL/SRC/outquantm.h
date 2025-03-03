#ifndef OUTQUANTM_H
#define OUTQUANTM_H

#include "iotools.h"
#include "selection.h"
#include "alias.h"
#include "galias.h"
#include <stdio.h>

class siftop;
class lcoordsys;
class gnodvalvm;

class outquantm
{
  public:
   /// constructor
   outquantm();
   /// destructor
   ~outquantm();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);
   /// prints values for selected items to grahics file
   void print_quant(long lcid, double time, long stepid, lcoordsys *lcsarray, gnodvalvm &nv, resfilefmt fmt);

   /// collects required quantity values at all points
   void collect_quant_val(long lcid, lcoordsys *lcsarray, gnodvalvm &nv);
   /// collects the given quantity values at all nodes and store them in the auxiliary arrays qvs or qvv
   void collect_nodal_quant_val(long lcid, gnodvalvm &nv);
   /// collects the given quantity values on all elements and store them in the auxiliary arrays qvs or qvv
   void collect_elem_quant_val(long lcid);
   /// transforms values in auxiliary array qvv to the required coordinate system
   void make_coord_transf(lcoordsys *lcsarray);
   /// makes quantity value scaling or other transformation
   void make_scale_transf();
   /// returns position %vector of the given point
   void give_pnt_coord(long pid, vector &pc);
   /// returns stress/strain state indicator in the given point
   strastrestate give_pnt_ssst(long pid);
   /// zeroes values in auxiliary array qvs
   void clean_qvs();
   /// zeroes values in auxiliary array qvv
   void clean_qvv();
   /// returns pointer to the string with default label of the required quantity
   static void give_quant_deflabel(mechquant mq, mechquant rmq, long maxll, char *label);   
   /// returns pointer to the string with default label of the required quantity component
   static void give_quant_comp_deflabel(mechquant mq, mechquant rmq, strastrestate ssst,
                                        long compid, long maxll, char *label);

   /// name of the required mechanical quantity
   mechquant mqn;
   /// name of the referenced quantity if the mqn represents an ordinary tensor
   mechquant reftensq;
   /// identifier of tensor type (stress/strain/other) for tensors stored in other or nonmech arrays 
   strastre other_tens_strastre;
   /// maximum length of quantity label string
   static const long maxl_qlabel = 256;
   /// quantity label (if '@' is the first character then automatic label will be generated)
   char qlabel[maxl_qlabel];
   /// indicator whether to use default quantity name
   bool defqlabel;
   /// type of point where the quantity is to be printed out (node/elem ipp/elem UDP)
   nodip tpnt;
   /// total number of points (nodes or elems according to tpnt) defined in the problem solved, i.e. Mt->nn or Mm->tnip
   long tnpnt;
   /// selection nodal/element ids where the quantity will be printed out
   sel selid;
   /// selected integration points on elements given by selid
   sel seleip;
   /// number of user defined points on selected elements
   long nudp;
   /// natural coordinates of selected user defined points
   double *ksi, *eta, *zeta;
   /// the given quantity representation (scalar/vector/tensor2d) in the output format
   quantrep qrf;
   /// indices of selected quantity components (displacement comps, strain/stress comps, ...)
   sel seliq;
   /// quantity scaling/tranformation flag
   answertype qstf;
   /// function for scaling/tranformation of the quantity
   gfunct qst;
   /// local coordinate system for vector and tensor2d quantities
   long lcs;

   // auxiliary members determined in mechmat 

   /// natural mathematical representation of the given quantity (scalar/vector/tensor2) determined in mechmat
   quantrep qr;
   /// number of quantity components in the case of vector or tensor quantity in the natural representation 
   long ncmp;
   /// reduced number of quantity components in the case of vector or tensor quantity according to the component selection
   long ncmpr;   
   /// tensor quantity storage notation (Voigt reduced/Voigt full/Full matrix 3x3)
   tensqnot tsn;
   /// tensor quantity type identifier (strain/stress/other)
   strastre tti;

   // auxiliary members not read in the quantity record, read/handled in outresfile, outdriverm
  
   /// pointer to the output file name string
   char *outfname;
   /// pointer to the associated output file of the given quantity
   FILE *outf;
   /// index shift of the first element
   long ide1;
   /// index shift of the first node
   long idn1;
   /** Auxiliary array of quantity values, qvv[i] represents pointer to array of the given vector/tensor2d quantity components 
       at i-th node/element ip/element UDP. In the case of elements ip or UDP, qvv[i][j] is the value of j-th selected 
       quantity component at k-th element ip or UDP and ncomp is the maximum number of vector/tensor2d components.
   */
   static double **qvv;
   /** Auxiliary array of quantity values, qvv[i] represents value of the selected quantity cpomponent
       at i-th node/element ip/element UDP in the case of scalar output format */
   static double *qvs;
   /// maximum number of points where the quantities will be printed out, it is the first dimension of arrays qvv and qvs
   static long maxnpnt;
   /// maximum number of components of quantities that will be printed out, it is the second dimension of array qvv
   static long maxncmp;
   /// parameter names of the general function transformation matrix
   static const char *namepar[];
   /// maximum number of iterations in Jacobi's method for the determination of the tensor's principal values
   static long ni_jac;
   /// tolerance of Jacobi's method for the determination of the tensor's principal values
   static double err_jac;
   /// zero treshold for the distinguishing of zero values in Jacobi's method
   static double zero;
};

#endif
