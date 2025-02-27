#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
// PETSC
class PETSC_CONT{
 public:
  // PETSC matrix
  Mat fact;
  // PETSC vector
  Vec r_rhs;
  // PETSC vector
  Vec z_lhs;
  PetscInt *ix ;
  PetscScalar *val;
  
  KSP                solver;
  PC                 prec;
  Mat                A;
  Vec                x,b,d;
  PetscErrorCode     ierr;
  PetscReal rtol;
  PetscReal abstol;
  PetscReal dtol;
  PetscInt maxits;
};

