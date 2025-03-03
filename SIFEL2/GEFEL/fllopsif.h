#ifndef FLLOPSIF_H
#define FLLOPSIF_H
#ifdef INC_PERMON
 #include "permonaif.h"
#endif

class fllopsif
{
 public:
  fllopsif();
  ~fllopsif();
  int solver_fllop (long n, long *kci, long *kadr, double *ka, double *rhs, double *lhs, 
                    long d, long bm, long *bci, long *badr, double *ba, double *ker);
#ifdef INC_PERMON
  PetscInt *kci_p;
  PetscInt *kadr_p;
  long nkci;  // number of components of kci_p array
  long nkadr; // number of components of kadr_p array
#endif
};

#endif
