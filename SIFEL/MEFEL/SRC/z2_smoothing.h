#ifndef Z2_ESTIMATOR_H
#define Z2_ESTIMATOR_H

struct vector;
struct ivector;
class skyline;

class z2_smoothing
{
 public:
  z2_smoothing (long ncomp);
  ~z2_smoothing (void);

  void run (vector *rsigfull);

  void give_adapt_code_numbers (long eid,ivector &cn);
  void compute_ntdbr (vector &ntdbr);

  void alloc_ntn (skyline *ntn_sky);
  void column_lengths_nn (skyline *ntn_sky);
  void compute_ntn_sky (skyline *ntn_sky);
  void compute_rsigfull (skyline *ntn_sky,vector &ntdbr, vector *rsigfull);

  void compute_ainv (vector &ainv);
  void compute_rsigfull (vector &ainv,vector &ntdbr, vector *rsigfull);

  long ncomp;
  long nn;
  long ne;
  long flags;
};

#endif





