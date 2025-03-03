#ifndef OUTDRIVERM_H
#define OUTDRIVERM_H

#include "galias.h"
#include "selection.h"
#include "xfile.h"
#include "alias.h"
#include "outdiagm.h"
#include "outresfilem.h"
#include "siftop.h"
#include "lcoordsys.h"

#include <stdio.h>

#ifndef FNAMELEN
 #define FNAMELEN 1001
#endif

//class nodeoutm;
//class nodeoutgm;
//class elemoutm;
//class elemoutgm;
//class pointoutm;
//class outdriverm;
struct matrix;
struct gfmatrix;
class outquantm;


class nodeoutm
{
  public :
   /// constructor
   nodeoutm();
   /// destructor
   ~nodeoutm();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected nodes
   void print_out(FILE *out, long lcid);
   /// prints displacements at selected nodes
   void print_disp(FILE *out, long lcid);
   /// prints averaged strains at selected nodes
   void print_stra(FILE *out, long lcid);
   /// prints averaged stresses at selected nodes
   void print_stre(FILE *out, long lcid);
   /// prints averaged other values at selected nodes
   void print_other(FILE *out);
   /// prints reactions at selected nodes
   void print_react(FILE *out);
   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);

   /// prints nodal values at selected steps
   sel dstep;
   /// selections for load cases 
   sel  sellc;

   /// selection of nodes for displacements
   sel selndisp;
   /// selections for displacements
   sel  *seldisp;

   /// selection of nodes for strains
   sel selnstra;
   /// selections for strains
   sel  *selstra;
   /// indices of local coordinate systems for strains
   long *transtra;

   /// selection of nodes for stresses
   sel selnstre;
   /// selections for stresses
   sel  *selstre;
   /// indices of local coordinate systems for stresses
   long *transtre;

   /// selection of nodes for other values
   sel selnoth;
   /// selections for other values
   sel  *seloth;
   
   /// reactions output indicator
   long react;    
};



class nodeoutgm
{
  public :
   /// constructor
   nodeoutgm ();
   /// destructor
   ~nodeoutgm ();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected items to grahics file
   void print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, double *ifor);
   /// prints values for selected items to several grahics files
   void print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf, double *ifor);
   /// prints values for selected items to grahics file
   void print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, double *ifor, double *fr);
   /// prints values for selected items to several grahics files
   void print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf, double *ifor, double *fr);

   /// prints values of selected strains as scalars to graphics file
   void print_gr_stra_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as scalars to graphics file
   void print_gr_stre_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as scalars to graphics file
   void print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as vectors to graphics file
   void print_gr_stra_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as vectors to graphics file
   void print_gr_stre_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as vectors to graphics file
   void print_gr_oth_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as tensors to graphics file
   void print_gr_stra_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as tensors to graphics file
   void print_gr_stre_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as tensors to graphics file
   void print_gr_oth_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as scalars to several graphics files
   void print_gr_stra_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as scalars to several graphics files
   void print_gr_stre_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as scalars to several graphics files
   void print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// prints values of selected strains as vectors to several graphics files
   void print_gr_stra_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as vectors to several graphics files
   void print_gr_stre_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as vectors to several graphics files
   void print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// prints values of selected strains as tensors to several graphics files
   void print_gr_stra_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as tensors to sevral graphics files
   void print_gr_stre_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as tensors to several graphics files
   void print_gr_oth_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);

   /// prints nodal values at selected steps
   sel dstep;

   /// selections for load cases 
   sel  sellc;

   /// selection of nodes for displacements
   sel selndisp;
   /// selections for displacements
   sel  *seldisp;

   /// selection of items for strains
   sel selnstra;
   /// selections for strains
   sel *selstra;
   /// indices of local coordinate systems for strains
   long *transtra;

   /// selection of items for stresses
   sel selnstre;
   /// selections for stresses
   sel *selstre;
   /// indices of local coordinate systems for stresses
   long *transtre;

   /// selection of items for other values
   sel selnoth;
   /// selections for other values
   sel  *seloth;
   /// selection of items for forces
   sel selnforce;
   /// selections for components of forces
   sel *selforce;
};



class elemoutm
{
  public :
   /// constructor
   elemoutm();
   /// destructor
   ~elemoutm();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected elements
   void print_out(FILE *out, long lcid);
   /// prints strains at selected elements
   void print_stra(FILE *out, long lcid);
   /// prints stresses at selected elements
   void print_stre(FILE *out, long stre);
   /// prints other values at selected elements
   void print_other(FILE *out);
   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);

   /// prints element values at selected steps
   sel dstep;

   /// selections for load cases 
   sel sellc;

   /// selection of elements for strains
   sel selestra;
   /// selections for strains
   sel *selstra;
   /// indices of local coordinate systems for strains
   long *transtra;

   /// selection of elements for stresses
   sel selestre;
   /// selections for stresses
   sel *selstre;
   /// indices of local coordinate systems for stresses
   long *transtre;

   /// selection of elements for other values
   sel seleoth;
   /// selections for other values
   sel  *seloth;
};



class elemoutgm
{
  public :
   /// constructor
   elemoutgm ();
   /// destructor
   ~elemoutgm ();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected items to grahics file
   void print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, long idelem1);
   /// prints values for selected items to to several grahics file
   void print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf, long idelem1);

   /// prints values of selected strains as scalars to graphics file
   void print_gr_stra_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as scalars to graphics file
   void print_gr_stre_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as scalars to graphics file
   void print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as vectors to graphics file
   void print_gr_stra_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as vectors to graphics file
   void print_gr_stre_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as vectors to graphics file
   void print_gr_oth_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as tensors to graphics file
   void print_gr_stra_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected stresses as tensors to graphics file
   void print_gr_stre_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);
   /// prints values of selected other items as tensors to graphics file
   void print_gr_oth_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf);

   /// prints values of selected strains as scalars to several graphics files
   void print_gr_stra_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as scalars to several graphics files
   void print_gr_stre_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as scalars to several graphics files
   void print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// prints values of selected strains as vectors to several graphics files
   void print_gr_stra_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as vectors to several graphics files
   void print_gr_stre_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as vectors to several graphics files
   void print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// prints values of selected strains as tensors to several graphics files
   void print_gr_stra_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected stresses as tensors to sevral graphics files
   void print_gr_stre_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);
   /// prints values of selected other items as tensors to several graphics files
   void print_gr_oth_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid);

   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);

   /// prints element values at selected steps
   sel dstep;

   /// selections for load cases 
   sel  sellc;

   /// selection of items for strains
   sel selestra;
   /// selections for strains
   sel *selstra;
   /// indices of local coordinate systems for strains
   long *transtra;

   /// selection of items for stresses
   sel selestre;
   /// selections for stresses
   sel *selstre;
   /// indices of local coordinate systems for stresses
   long *transtre;

   /// selection of items for other values
   sel seleoth;
   /// selections for other values
   sel  *seloth;
   
   /// number of the first element for GiD output (normally should be 1), set by print_graphics function
   long ide1;

   /// number of local coordinate systems used at output procedures
   long nlcs;
   /// array of transformation matrices for individual local coordinate systems
   matrix *lcs;
};



class pointoutm
{
  public :
   /// constructor
   pointoutm();
   /// destructor
   ~pointoutm();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for given user defined point
   void print_out(FILE *out, long lcid);
   /// prints strains at given user defined point
   void print_stra(FILE *out, long lcid);
   /// prints stresses at given user defined point
   void print_stre(FILE *out, long lcid);
   /// prints other values at given user defined point
   void print_other(FILE *out, long pid);

   /// prints UDP values at selected steps
   sel dstep;

   /// number of user defined points
   long npnt;

   /// x-coordinate
   double *ksi;
   /// y-coordinate
   double *eta;
   /// z-coordinate
   double *zeta;

   ///  selection of elements
   sel  selelem;

   /// selections for load cases 
   sel  sellc;

   ///  selections of points
   sel  *selpnt;

   /// selections for strains
   sel  *selstra;
    /// indices of local coordinate systems for strains
   long *transtra;

   /// selections for stresses
   sel  *selstre;
   /// indices of local coordinate systems for stresses
   long *transtre;

   /// selections for other values
   sel  *seloth;
};



class outdriverm
{
  public :
   /// constructor
   outdriverm();
   /// destructor
   ~outdriverm();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints header for output file
   void print_header(FILE *out);
   /// prints header for new step
   void print_newstep(FILE *out, long lcid, long istep, double time);
   /// prints values to the text file
   void print_out(FILE *out, long lcid, long istep, double time);
   /// forced print of values to the text file
   void print_out_forced(FILE *out, long lcid, long istep, double time);
   /// prints diagrams
   void print_diags(long lcid, double lambda, long istep, double *fi);
   /// forced print of diagrams
   void print_diags_forced(long lcid, double lambda, long istep, double *fi);
   /// prints diagrams
   void print_diags(long lcid, double lambda, long istep, double *fi, double *fr);
   /// forced print of diagrams
   void print_diags_forced(long lcid, double lambda, long istep, double *fi, double *fr);
   /// prints graphics
   void print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi);
   /// forced print of graphics
   void print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi);
   /// prints graphics
   void print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr);
   /// forced print of graphics
   void print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr);
   /// creates output graphics files and their headers for GiD separated format (grfmt_gidsp)
   void create_files_gidsp(const char *mode);
   /// creates output graphics files and their headers for general quantity output (grfmt_quant)
   void create_files_quant(const char *mode, long lcid, long istep, const char *dlcid);
   /// prints result files for the given load case and (time) step
   void print_sep_files_quant(long lcid, double time, long istep, double *r, double *fl, double *fi, double *fr);
   /// prints header to the result file in plain text output format
   void print_header_plain_out(FILE *out);

   /// tests the output setting whether it requires some local coordinate system
   long testlcs();
   /// converts selections given by property id to list or range type
   void conv_sel_prop(siftop *top);

   /// output text filename
   char outfn[FNAMELEN];
   /// output diagrams filename
   char outdiagfn[FNAMELEN];
   /// output graphics filename
   char outgrfn[FNAMELEN];
   /// generated output graphic file name for grfmt_gid_sep format
   char outgrfngs[FNAMELEN+50];

   /// text output flag
   flagsw textout;
   /// output text file
   FILE *outf;
   /// output diagram files
   FILE **outdiagf;
   /// output graphics file
   FILE *outgr;

   /// description what will be printed from nodal values
   nodeoutm  no;
   /// description what will be printed from element values
   elemoutm  eo;
   /// description what will be printed for each user def. point
   pointoutm po;

   /// description what will be printed from nodal values for graphics
   nodeoutgm  nog;
   /// description what will be printed from element values for graphics
   elemoutgm  eog;
   

   /// number of local coordinate systems  
   long nlcs;
   /**
       array of transformation matrices of the individual local coordinate systems
       rows of the matrix contains base vectors of the given lcs in the global 
       coordinate system  {e'} = [T]*{e} = lcs[i]*{e}
   */
   lcoordsys *lcs;

   /// format of grahics file - 0 - none, 1 - openDX, 2 - FemCAD, 3 - GiD
   graphfmt gf;
   /// number of cuts exported from 3d mesh to GiD
   long ncut;
   /// number of the first node for GiD output (normally should be 1), set by print_init function
   long idn1;
   /// number of the first element for GiD output (normally should be 1), set by print_init function
   long ide1;

   /// number of output diagrams
   long ndiag;
   /// array with description of printed diagrams
   outdiagm *odiag;
   
   ///number of file sequence in VTK graphic output
  long vtk_num;
  /// number of output files with given quantities
  long nofq;
  /// array of file and quantity description for result output in the case gf == grfmt_quant
  outresfilem *ofq;
  /** 
     array of indicators for creation of output file with quantities: cf[i] = 1 -> create i-th file, cf[i] = j < 0 -> i-th file is the same as -(j+1)-th,
     dimension of cf is nofq.
  */
  ivector cf;
};

#endif
