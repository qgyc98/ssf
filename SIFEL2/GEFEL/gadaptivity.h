#ifndef GADAPTIVITY_H
#define GADAPTIVITY_H

#include <stdio.h>
class gtopology;
struct ivector;
struct matrix;
struct vector;

void give_der_star (vector &bf, vector *rderfull,ivector &nodes,vector &der_star,long nn);
void ntnmtov (matrix &ntnm,vector &ntnv);
void print_contoures (gtopology *gt,char *file,double **midpoints,double *valel);
void print_confile (const char *file,long nnod,long nel,double **auxxyzv,long **auxnod,long dim);

void fprintf_matrix (FILE *stream, matrix &mx, char s[]);
void fprintf_vector (FILE *stream,vector &v,char s[],long c);
void fprintf_ivector (FILE *stream,ivector &v,char s[],long c);

void fprintf_d_1D (FILE *stream,double *p,long p_n,char s[],long c);
void fprintf_l_1D (FILE *stream,long *p,long p_n,char s[],long c);
void fprintf_l_2D (FILE *stream,long **p,long p_n,long p_m,char s[],long c);

long edge_position (gtopology *gt,long node1,long node2);
long surface_position (gtopology *gt,long node1,long node2,long node3);

long adjelem2edge (gtopology *gt,long node1,long node2,long eid);

long opposite_node (gtopology *gt,long *nod,long eid);

void print_valel (FILE *stream, gtopology *gt, const char *path, const char *file, const char *caption, double *valel, char flag);

//long common_numbers (long n1,long *p1,long n2,long *p2);

/// SPR Superconvergent patch recovery


class patch_averaging
{
 private:
  ///  actual gtopology
  gtopology *gt;
  
  ///  FOLLOWING CHARACTERISTICS MUST BE IDENTICAL FOR EVERY ELEMENT ALL OVER THE DOMAIN
  ///  dim,ncomp(ncompother),deg(polynom degree of base functions),material
  
  long dim;     ///  dimension (2 or 3)
  long nvals;   ///  number of smoothed ~ aproximated values
  long nn;      ///  number nodes at domain
  long ne;      ///  number elements at domain
  // ...        ///  number of coeficients of patch polynom
  
  
  /////  array of natural coordinates of sampling points, dimension (gt->ne ; dim*nsp[i]), for 3D, nsp[i]==2  =>  spcoord[i]==(x1; y1; z1; x2; y2; z2)
  double **spcoord;
  /////  array of rough values in sampling points, dimension (gt->ne ; ncomp)
  //double **spvalue;
  /////  array of extreme element natural coordinates, dimension (gt->ne ; 4), maxcoord[i]==(max x; min x; max y; min y)
  double **maxcoord;
  /////  array indicating position of nodes on domain, dimension (gt->nn ; ???)
  long **insidenod;
  /////  array indicating position of elements on domain, dimension (gt->ne ; ???)
  //long *insidelem;
  //
  //
  /////
  long flag;
  
 public:
  /// CONSTRUCTOR
  patch_averaging (gtopology *gt, long dim, long nvals, long flag);
  /// DESTRUCTOR
  ~patch_averaging (void);
  
  /// solve - main function
  void solve (FILE *out, const matrix *spvalue, vector *nodvalue);
  
 private:
  ///
  void nsp_spcoord_maxcoord_assembling (const matrix *spvalue, char nsma_flag);
  void insidenod_assembling            (void);
  void compute_patches_spr             (const matrix *spvalue, vector *nodvalue, long cut_flag);
  ///
  void normal_coordinates_ae   (long nadjel, const long *adjel,       vector &coef_normcoord);
  void polynom_coefficients_ae (long nadjel, const long *adjel, const vector &coef_normcoord,       vector *a, const matrix *spvalue);
  void nodvalue_assembling_ae  (long nadjel, const long *adjel, const vector &coef_normcoord, const vector *a, vector *nodvalue, long nid, ivector &magnitude);
  
  ///
  void polynom (long ncoef, const double *normcoord, double *p);
  ///
  void sigma (vector *nodvalue, long ncoef,long nid,const vector &coef_normcoord,const vector *a,ivector &magnitude);
  
};


#endif
