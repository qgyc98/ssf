#ifndef GLOBALG_H
#define GLOBALG_H

#ifndef EXTERN
 #define EXTERN extern
#endif

#if defined(INC_OPENMP) || defined(INC_CHMD)
  /// number of threads used in OpenMP
  EXTERN int Numth;
#endif

#endif
