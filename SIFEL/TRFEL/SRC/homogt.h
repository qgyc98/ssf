#include "vector.h"
#include "matrix.h"

class HomogData 
{
 public:

  HomogData();
  HomogData(double *PhaseMatrix);
  HomogData(const double iE_m, const double inu_m, const double if_m,
	    const double iE_i, const double inu_i, const double iHirsch);

  ~HomogData(void){};

  void Voigt(void);
  void Reuss(void);
  void Hirsch(void);
  void Hansen(void);
  void Counto(void);
  void MoriTanaka(void);
  void MT_mtrx(double *PhaseMatrix, int NumRows, int MatrixRow);
  void SelfConsistent(void);
  void SCS(double *PhaseMatrix, int NumRows);
  void Voigt(double *PhaseMatrix, int NumRows);
  void Reuss(double *PhaseMatrix, int NumRows);
  void HashinShtrikman(void);
  void WalpoleMulti(double *PhaseMatrix, int NumRows);
  void KusterToksoz(void);
  void HerveZaoui(double *PhaseMatrix, int NumRows);

  double E_m;  // Young's modulus of matrix phase
  double nu_m; // Poisson's ration of matrix phase
  double f_m;  // Volume fraction of matrix phase
  double k_m;  // Bulk modulus of materix
  double mu_m; // Effective shear modulus of matrix

  double E_i;  // Young's modulus of inclusion
  double nu_i; // Poisson's ration of inclusion 
  double f_i;  // Volume fraction of inclusion
  double k_i;  // Bulk modulus of inclusion
  double mu_i; // Effective shear modulus of inclusionnction to call to select a buffer from the buffers tab. *

  double E_hmg;
  double nu_hmg;
  double k_hmg, k_hmg_2;
  double mu_hmg, mu_hmg_2;

  double E_hmg_2;  //auxiliary for lower and upper values
  double nu_hmg_2; //auxiliary for lower and upper values
 
  double H_chi;  //Hirsch factor for model

  double temp1, temp2, temp3, temp4; //help variables 

 private:
  double alpha_m; // Auxiliary factor for bulk response
  double beta_m;  // Auxiliary factor for shear response
 
  void ENuToKMu( const double oE, const double onu, double &ok, double &omu ); 
  void KMuToENu( const double ok, const double omu, double &oE, double &onu ); 
  void Swap(double &A, double &B);
  double Lambda(double *PhaseMatrixKMu, int NumRows, double mu);
  double Gamma(double *PhaseMatrixKMu, int NumRows, double k);
  double Zeta(double k, double mu);
  void fillJ(matrix &J, double r, const vector &mu, const vector &k, int phase);
  void fillL(matrix &L, double r, const vector &mu, const vector &k, int phase);
};

