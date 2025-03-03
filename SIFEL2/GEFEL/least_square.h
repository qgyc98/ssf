#ifndef LEAST_SQUARE_H
#define LEAST_SQUARE_H

class gtopology;
struct matrix;
struct vector;
struct ivector;


class least_square
{
 public:
  least_square (long dim,long ncomp,gtopology *gt,long flag);
  ~least_square (void);
  
  void spr_default (FILE *out,double **spvalue,double *nodvalue);
  void spr_optional (FILE *out,long *nsp,double **spcoord,double **spvalue,double *nodvalue,long cut_flag);
  void L2_nod2sp (FILE *out,long *nsp,double **spcoord,double *nodvalue,double **spvalue,char typ_patch);
  void L2_sp2sp (FILE *out,double **spvalue,long nap,const long *parentel,const double **apcoord,double **apvalue);
  
 private:
  void nsp_spcoord_maxcoord_assembling (char nsma_flag);
  void insidenod_assembling (void);
  void adjap_assembling (long nap,const long *parentel,const double **apcoord,long *nadjap,long **adjap);
  
  void compute_patches_spr (long cut_flag);
  void compute_patches_nod2sp (char typ_patch);
  void compute_patches_sp2sp (long nap,const double **apcoord,const long *nadjap,const long **adjap,double **apvalue);
  
  void normal_coordinates_ae (vector &coef_normcoord);
  void polynom_coefficients_ae (const vector &coef_normcoord,vector *a);
  void polynom_coefficients_inv_ae (const vector &coef_normcoord,vector *a,char typ_patch);
  
  void nodvalue_assembling_ae (const vector &coef_normcoord,const vector *a,long nid,ivector &magnitude);
  void spvalue_assembling_patch_ae (const vector &coef_normcoord,const vector *a,long id,double **magnitude);
  void apvalue_assembling (const vector &coef_normcoord,const vector *a,long nadjap,const long *adjap,const double **apcoord,double **apvalue);
  
  void sigma (long ncoef,long nid,const vector &coef_normcoord,const vector *a,ivector &magnitude);
  void polynom (long ncoef,const double *normcoord,double *p);
  
  
  ///  actual gtopology
  gtopology *gt;
  
  ///  FOLLOWING CHARACTERISTICS MUST BE IDENTICAL FOR EVERY ELEMENT ALL OVER THE DOMAIN
  ///  dim,ncomp(ncompother),deg(polynom degree of base functions),material
  
  ///  dimension (2 or 3)
  long dim;
  ///  number of smoothed ~ aproximated values
  long nvals;
  ///  number nodes on domain
  long nn;
  ///  number elements on domain
  long ne;
  ///  number of coeficients of patch polynom
  
  
  ///  number of sampling points on constituent element, dimension == gt->ne
  long *nsp;
  ///  array of natural coordinates of sampling points, dimension (gt->ne ; dim*nsp[i]), for 3D, nsp[i]==2  =>  spcoord[i]==(x1; y1; z1; x2; y2; z2)
  double **spcoord;
  ///  array of rough values in sampling points, dimension (gt->ne ; ncomp)
  double **spvalue;
  ///  
  double *nodvalue;
  
  ///  array of extreme element natural coordinates, dimension (gt->ne ; 4), maxcoord[i]==(max x; min x; max y; min y)
  double **maxcoord;
  ///  array indicating position of nodes on domain, dimension (gt->nn ; ???)
  long **insidenod;
  ///  array indicating position of elements on domain, dimension (gt->ne ; ???)
  long *insidelem;
  
  ///
  long nadjel;
  ///
  long *adjel;
  
  ///
  long flag;
};

#endif
