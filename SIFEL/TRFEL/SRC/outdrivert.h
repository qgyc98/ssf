#ifndef OUTDRIVERT_H
#define OUTDRIVERT_H

#include <stdio.h>
#include "galias.h"
#include "selection.h"
#include "xfile.h"
#include "aliast.h"
#include "outdiagt.h"

#define FNAMELEN 1001

//class nodeoutt;
//class nodeoutgt;
//class elemoutt;
//class elemoutgt;
//class outdrivert;

class nodeoutt
{
  public :
   /// constructor
   nodeoutt();
   /// destructor
   ~nodeoutt();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected nodes
   void print_out(FILE *out, long lcid);
   /// prints unknowns at selected nodes
   void print_unkn(FILE *out, long lcid);
   /// prints gradients at selected nodes
   void print_grad(FILE *out, long lcid);
   /// prints fluxes at selected nodes
   void print_flux(FILE *out, long lcid);
   /// prints other values at selected nodes
   void print_other(FILE *out);
   /// prints eq_other values at selected nodes
   void print_eqother(FILE *out);
   /// prints optional fluxes at selected nodes
   void print_forcet(FILE *out); 
   /// converts selections according to the givne property to the standard types (list or range)
   void conv_sel_prop(siftop *top);

   /// prints nodal values at selected steps
   sel dstep;
   
   /// selection of nodes for unknowns
   sel selnunkn;
   /// selections for displacements
   sel *selunkn;
   
   /// selection of nodes for gradients
   sel selngrad;
   /// selections for strains
   sel  *selgrad;
   
   /// selection of nodes for stresses
   sel selnflux;
   /// selections for stresses
   sel  *selflux; 

   /// selection of nodes for other values
   sel selnoth;
   /// selections for other values
   sel  *seloth; 

   /// selection of nodes for eq_other values
   sel selneqoth;
   /// selections for eq_other values
   sel  *seleqoth; 
 
   /// selection of items for forces=fluxes
   sel selnforcet;
   /// selections for components of forces=fluxes
   sel *selforcet;
};


class nodeoutgt
{
  public :
   /// constructor
   nodeoutgt();
   /// destructor
   ~nodeoutgt();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected nodes to graphics file
   void print_graphics(FILE *out, long lcid, const char *desclcid, graphftt gf, double *ifor);
   /// prints values for selected nodes to to several grahics file
   void print_graphics(char *outfn, const char *mode, long lcid, const char *desclcid, graphftt gf, double *ifor, hflagsw hf);
   /// converts selections according to the givne property to the standard types (list or range)
   void conv_sel_prop(siftop *top);
   
   /// prints nodal values at selected steps
   sel dstep;
   
   /// selection of nodes for unknowns
   sel selnunkn;
   /// selections for unknowns
   sel *selunkn;
   
   /// selection of nodes for gradients
   sel selngrad;
   /// selections for gradients
   sel  *selgrad;
   
   /// selection of nodes for fluxes
   sel selnflux;
   /// selections for fluxes
   sel  *selflux;  

   /// selection of nodes for other values
   sel selnoth;
   /// selections for other values
   sel  *seloth; 

   /// selection of nodes for eq_other values
   sel selneqoth;
   /// selections for eq_other values
   sel  *seleqoth;

   /// selection of items for forces=fluxes
   sel selnforcet;
   /// selections for components of forces=fluxes
   sel *selforcet;
};


class elemoutt
{
  public :
   /// constructor
   elemoutt();
   /// destructor
   ~elemoutt();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected elements
   void print_out(FILE *out, long lcid);
   /// prints unknowns at selected elements
   void print_unkn(FILE *out, long lcid);
   /// prints gradients at selected elements
   void print_grad(FILE *out, long lcid);
   /// prints fluxes at selected elements
   void print_flux(FILE *out, long lcid);
   /// prints averaged other values at selected elements
   void print_other(FILE *out);
   /// prints averaged eq_other values at selected elements
   void print_eqother(FILE *out);
   /// converts selections according to the givne property to the standard types (list or range)
   void conv_sel_prop(siftop *top);
   
   /// prints element values selected steps
   sel dstep;
   
   /// selection of elements for gradients
   sel selegrad;
   /// selections for gradients
   sel *selgrad;
   
   /// selection of elements for fluxes
   sel seleflux;
   /// selections for fluxes
   sel *selflux;
   
   /// selection of elements for other values
   sel seleoth;
   /// selections for other values
   sel *seloth; 

   /// selection of elements for eq_other values
   sel seleeqoth;
   /// selections for eq_other values
   sel *seleqoth;
};



class elemoutgt
{
  public :
   /// constructor
   elemoutgt();
   /// destructor
   ~elemoutgt();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for selected elements to grahics file
   void print_graphics(FILE *out, long lcid, const char *desclcid, graphftt gf, long idelem1);
   /// prints gradients for selected elements to grahics file as individual scalars
   void print_gr_grad_scal(FILE *out, const char *desclcid, graphftt gf);
   /// prints fluxes for selected elements to grahics file as individual scalars
   void print_gr_flux_scal(FILE *out, const char *desclcid, graphftt gf);
   ///  prints other values for selected elements to grahics file as individual scalars
   void print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphftt gf);
   ///  prints eqother values for selected elements to grahics file as individual scalars
   void print_gr_eqoth_scal(FILE *out, long lcid, const char *desclcid, graphftt gf);
   /// prints gradients for selected elements to grahics file as vectors
   void print_gr_grad_vec(FILE *out, long lcid, const char *desclcid, graphftt gf);
   /// prints fluxes for selected elements to grahics file as vectors 
   void print_gr_flux_vec (FILE *out, long lcid, const char *desclcid, graphftt gf);
   ///  prints other values for selected elements to grahics file as vectors
   void print_gr_oth_vec (FILE *out, long lcid, const char *desclcid, graphftt gf);
   ///  prints eqother values for selected elements to grahics file as vectors
   void print_gr_eqoth_vec (FILE *out, long lcid, const char *desclcid, graphftt gf);  

   /// prints values for selected elements to several grahics file
   void print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphftt gf, hflagsw hf, long idelem1);
   /// prints gradients for selected elements to several grahics file as individual scalars
   void print_gr_grad_scal(const char *outfn, const char *mode, const char *desclcid, hflagsw hf);
   /// prints fluxes for selected elements to several grahics file as individual scalars
   void print_gr_flux_scal(const char *outfn, const char *mode, const char *desclcid, hflagsw hf);
   /// prints other values for selected elements to several grahics file as individual scalars
   void print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// prints eqother values for selected elements to several grahics file as individual scalars
   void print_gr_eqoth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// prints gradients for selected elements to several grahics file as vectors
   void print_gr_grad_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// prints fluxes for selected elements to several grahics file as vectors
   void print_gr_flux_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// prints other values for selected elements to several grahics file as vectors
   void print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// prints eqother values for selected elements to several grahics file as vectors
   void print_gr_eqoth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf);
   /// converts selections according to the givne property to the standard types (list or range)
   void conv_sel_prop(siftop *top);
  
   /// prints element values selected steps
   sel dstep;
   
   /// selection of elements for gradients
   sel selegrad;
   /// selections for gradients
   sel *selgrad;
   
   /// selection of elements for fluxes
   sel seleflux;
   /// selections for fluxes
   sel *selflux;
   
   /// selection of elements for other values
   sel seleoth;
   /// selections for other values
   sel *seloth;  

   /// selection of elements for eq_other values
   sel seleeqoth;
   /// selections for eq_other values
   sel *seleqoth; 

   /// number of the first element for GiD output (normally should be 1), set by print_graphics function
   long ide1;
};


class pointoutt
{
  public :
   /// constructor
   pointoutt();
   /// destructor
   ~pointoutt();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints values for given user defined point
   void print_out(FILE *out, long lcid);
   /// prints gradients at given user defined point
   void print_grad(FILE *out, long lcid);
   /// prints fluxes at given user defined point
   void print_flux(FILE *out, long lcid);
   /// prints other values at given user defined point
   void print_other(FILE *out, long pid);
   /// prints eq_other values at given user defined point
   void print_eqother(FILE *out, long pid);

   /// prints element values selected steps
   sel dstep;

   long npnt;

   /// x-coordinate
   double *ksi;
   /// y-coordinate
   double *eta;
   /// z-coordinate
   double *zeta;

   ///  selection of elements
   sel  selelem;

   /// selections for unknowns 
   sel  sellc;

   ///  selections of points
   sel  *selpnt;

   /// selections for strains
   sel  *selgrad;

   /// selections for stresses
   sel  *selflux;

   /// selections for other values
   sel  *seloth;

   /// selections for eq_other values
   sel  *seleqoth;
};


class outdrivert
{
  public :
   /// constructor
   outdrivert();
   /// destructor
   ~outdrivert();
   /// reads data from file
   long read(XFILE *in);
   /// prints description data to file
   void print(FILE *out);
   /// prints header for output file
   void print_header(FILE *out);
   /// prints header for new step
   void print_newstep(FILE *out, long lcid, long istep, double time);
   /// prints values at nodes
   void print_out(FILE *out, long lcid, long istep, double time);
   /// the forced print of values to the text file
   void print_out_forced(FILE *out, long lcid, long istep, double time);
   /// prints diagrams
   void print_diags(long lcid, double lambda, long istep, double *fi);
   /// the forced print of diagrams
   void print_diags_forced(long lcid, double lambda, long istep, double *fi);
   /// prints graphics
   void print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi);
   /// the forced print of graphics
   void print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi);
   /// creates output graphics files and their headers for GiD separated format (grftt_gidsp)
   void create_files_gidsp(const char *mode);
   /// prints footers for GiD separated format
   void close_files_gidsp();
   /// converts selections of quantities at nodes and elements according to the givne property to the standard types (list or range)
   void conv_sel_prop(siftop *top);
   /// output text filename
   char outfn[FNAMELEN];
   /// output diagrams filename
   char outdiagfn[FNAMELEN];
   /// output graphics filename
   char outgrfn[FNAMELEN];
   /// generated output graphic file name for grftt_gid_sep format
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
   nodeoutt  no;
   /// description what will be printed from element values
   elemoutt  eo;
   /// description what will be printed for each user def. point
   pointoutt po;
   
   /// description what will be printed from nodal values for graphics
   nodeoutgt  nog;
   /// description what will be printed from element values for graphics
   elemoutgt  eog;  

   /// number of local coordinate systems  
   long nlcs;
   /// number of components of base vectors
   long *nclcs;
   /// base vectors each lcs has base vectors stored in one row
   double **bvlcs;
   
   /// format of grahics file - 0 - none, 1 - openDX, 2 - FemCAD
   graphftt gf;
   /// number of cuts exported from 3d mesh to GiD
   long ncut;
   /// number of the first node for GiD output (normally should be 1), set by print_init function
   long idn1;
   /// number of the first element for GiD output (normally should be 1), set by print_init function
   long ide1; 

   ///number of file sequence in VTK graphic output
   long vtk_num;
  
   /// number of output diagrams
   long ndiag;
   /// array with description of printed diagrams
   outdiagt *odiag;
};

#endif
