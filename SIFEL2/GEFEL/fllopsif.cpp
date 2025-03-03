#ifdef INC_PERMON
 #include "permonaif.h"
#endif

#include "fllopsif.h"
#include "iotools.h"


fllopsif::fllopsif()
{
#ifdef INC_PERMON
  kci_p = NULL;
  kadr_p = NULL;
  nkci = 0L;
  nkadr = 0L;
#endif
}



fllopsif::~fllopsif()
{
#ifdef INC_PERMON
  delete [] kci_p;
  delete [] kadr_p;
#endif
}


/**
  @param n - number of unknowns on the domain for stiffness %matrix
  @param kci - array of of column indeces for stiffness %matrix
  @aparam kadr - array of of indeces of row beginnings for stiffness %matrix
  @param ka - array of nonzero stiffness %matrix components
  @param rhs - array of %vector of right hand side
  @param lhs - array of solution %vector
  @param d - number of defects
  @param bm - number of rows for constrained %matrix B
  @param bci - array of of column indeces for constrained %matrix B 
  @aparam badr - array of of indeces of the row beginnings for constrained %matrix B 
  @param ba - array of nonzero constrained %matrix B components
  @param ker - array of kernel %matrix components (d x n)

  Created by Tomas Koudelka & Jaroslav Kruis, according to V. Hapla hints
*/
int fllopsif::solver_fllop (long /*n*/, long */*kci*/, long */*kadr*/, double */*ka*/, double */*rhs*/, double */*lhs*/, 
                            long /*d*/, long /*bm*/, long */*bci*/, long */*badr*/, double */*ba*/, double */*ker*/)
{
 #ifdef INC_PERMON
  int err;
  int argc = 1;
  char *argv = "mefel";
  char **pargv = &argv;
  long i;
  
  if (kci_p == NULL)
  {
    nkci = kadr[n];
    kci_p = new PetscInt[kadr[n]];
  }
  else
    if (nkci != kadr[n])
    {      
      delete [] kci_p;
      nkci = kadr[n];
      kci_p = new PetscInt[kadr[n]];
    }
  if (kadr_p == NULL)
  {
    nkadr = n+1;
    kadr_p = new PetscInt[n+1];
  }
  else
    if (nkadr != n+1)
    {
      delete [] kci_p;
      nkadr = n+1;
      kadr_p = new PetscInt[n+1];
    }

  for (i=0; i< kadr[n]; i++)
    kci_p[i] = PetscInt(kci[i]);
  for (i=0; i< n+1; i++)
    kadr_p[i] = PetscInt(kadr[i]);


  FllopAIFInitialize(&argc, &pargv, "permon-options.txt");

  // zde se muze nastavit typ solveru, kterym to bude PETSc resit
  //
  // sdruzene gradienty 
  // PetscOptionsInsertString(NULL, "-qps_view_convergence -qps_ksp_type cg -qps_ksp_monitor -qps_rtol 1e-10 -pc_type jacobi -options_left -options_table -log_summary");
  // primy resic
  // PetscOptionsInsertString(NULL, "-qps_view_convergence -qps_ksp_type cg -qps_ksp_monitor -qps_rtol 1e-10 -pc_factor_mat_solver_package superlu_dist -pc_type lu -options_left -options_table -log_summary");

  // primy resic, kontakt
  // PetscOptionsInsertString(NULL, "-fllop_info -dual -qp_R_orth_type none -qps_view_convergence -qps_ksp_type cg -qps_ksp_monitor -qps_rtol 1e-10 -pc_factor_mat_solver_package superlu_dist -pc_type lu -options_left -options_table -log_summary");

  //  je-li matice v symmetric compressed rows
  //  je treba poskytnout horni trojuhelnik, ne nase symcomprow
  //  err = FllopAIFSetFETIOperator (n, adr, ci, a, QP_SYM_UPPER_TRIANGULAR);
  //  err = FllopAIFSetFETIOperator (n, kadr, kci, ka, QP_SYM_UNDEF, "K");
  // sdruzene gradienty, priprava na FETI
  // err = FllopAIFSetFETIOperator (n, kadr_p, kci_p, ka, AIF_MAT_SYM_UNDEF, "K");
  // primy resic, jednoprocesorove
  err = FllopAIFSetOperatorByStripes (n, n, n, kadr_p, kci_p, ka, AIF_MAT_SYM_UNDEF, "K");
  if (err)
  {
    print_err("initialization of FLLOP failed", __FILE__, __LINE__, __func__);
    abort();
  }

  //  rhs - nase prava strana na jedne domene
  FllopAIFSetRhs(n, rhs, "f");
	
  //  m dualni dimenze - pocet vsech multiplikatoru v uloze
  //  je treba sestavit matici vazeb B do compressed rows
  //  matice B ma tolik radku, kolik je multiplikatoru v cele uloze a tolik
  //  sloupcu, kolik je posunu na konkretni domene
  //  pokud chceme TFETI, nedefinujeme Dirichletovy podminky, pridame prvky do matice B

  //  vektor prave strany rovnostnich vazeb - skoky nebo slip mezi domenami
  //  v pripade klasicke dekompozice, sl je nulovy vektor
  // double *sl;
  // long m = n;
  //FllopAIFAddEq(bm, n, PETSC_FALSE, PETSC_TRUE, badr, bci, ba, "B", sl, "s");
  PetscInt br[16] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
  PetscInt bc[16] = {55, 63, 56, 64,   57, 65, 58, 66,   59, 67, 60, 68,   61, 69, 62, 70};
  double bb[16] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  FllopAIFSetEqCOO(8, n, PETSC_FALSE, PETSC_TRUE, br, bc, bb, 16, "Be", NULL, "ce");


  // matice vazeb pro kontakt
  PetscInt bri_[10] = {0,0,1,1,2,2,3,3,4,4}; 

  // toto bci_ je pro NEpodeprenou 2. podoblast
  // PetscInt bci_[8] = {3,34,9,42,15,50,21,58}; 

  // toto bci_ je pro z boku podeprenou 2. podoblast
  PetscInt bci_[10] = {3,33,9,40,15,48,21,56,27,72}; 

  double ba_[10] = {1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0};
  // vektor vzdalenosti mezi podoblastmi
  double gap[5] = {1.0,1.0,1.0,1.0,1.0};

  // matice nepronikani se nastavuje nasl.
  // gap je vektor mezery
  // 4 = pocet vazeb
  // 8 = pocet nenulovych prvku v matici b (delka pole bri)
  FllopAIFSetIneqCOO(5, n, PETSC_FALSE, PETSC_TRUE, bri_, bci_, ba_, 10, "Bi", gap, "g");

  //  lhs je vektor reseni
  FllopAIFSetSolutionVector(n, lhs, "r");

  FllopAIFFromOptions();

  FllopAIFSolve();

  // nasl. prikaz dat pouze po updateu matice, jinak se uchovava rozlozena matice
  //  FllopAIFReset();
 #endif

  return 0;
}
