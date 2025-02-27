#include "cusatismaterial.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#ifndef MAX
 #define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
 #define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif



/**
  Constructor initializes data members to zero or default values.

  Created by JK, 7. 3. 2021
*/
cusatismaterial::cusatismaterial(void) {
  //  Young's modulus of elasticity
  e = 0.0;
  //  Poisson's ratio
  nu = 0.0;
  // Normal modulus
  eN = 0.0;
  // Normal-Shear coupling
  alpha = 0.0;
  // Shear modulus
  eT = 0.0;
  // Tensile strength
  ft = 0.0;
  // Tensile fracture energy 
  gft = 0.0;
  // Tensile characteristic length
  lt = 0.0;
  // Shear strength
  fs = 0.0;
  // Shear fracture energy
  gfs = 0.0;
  // Shear characteristic length
  ls = 0.0;
  // Tension-shear interaction exponent
  nt = 0.0;
  // Initial hardening moduls
  hc0 = 0.0;
  // Transitional strain ratio
  kc0 = 0.0;
  // Deviatoric-volumetric strain ratio treshold
  kc1 = 0.0;
  // Deviatoric damage parameter - reduction of hc for increasing rdv
  kc2 = 0.0;
  // Volumetric-deviatoric coupling coefficient
  beta = 0.0;
  // Yielding compressive strength
  sigC0 = 0.0;
  // Densified normal modulus
  ed = 0.0;
  // Initial friction coefficient 
  m0 = 0.0;
  // Final friction coefficient
  mInf = 0.0;
  // Transitional normal stress for friction coefficient
  sigN0 = 0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK, 7. 3. 2021
*/
cusatismaterial::~cusatismaterial(void) {

}



/**
  Function reads material parameters from the opened text file.

  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::read(XFILE* in) {
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  xfscanf(in, "%k%lf %k%lf", "e", &e, "nu", &nu);
  
  // tady pridej vsechny parametry, ktere se ctou, dal jsem jich tam jen par
  //  zda se mi, ze nektere se pocitaji z jinych a neni treba je cist
  xfscanf (in,"%k%lf %k%lf %k%lf %k%lf","ft", &ft, "lt", &lt, "nt", &nt, "sigC0", &sigC0);

  init();
}

/**
  Function initializes material parameters after reading.

  @param Nothing

  @return The function does not return anything.

  11. 3. 2021, JV
*/
void cusatismaterial::init() {
  eN = e / (1. - 2. * nu);
  alpha = (1. - 4. * nu) / (1. + nu);

  if (alpha < 0.0) {
    print_err("Negtive alpha = %e, Poisson's number = %e is probably greater than 0.25, this is not allowed for Cusatis model.", __FILE__, __LINE__, __func__, alpha, nu);
    abort();
  }

  eT = alpha * eN;

  // characteristic length
  lt = 2. * eN * gft / (ft * ft); // lt = 2 * characteristic size (length)
  ls = 2. * alpha * eN * gfs / (fs * fs); // see Cusatis & Cedolin "Two-scale study of concrete fracturing behavior" 


  // @todo check for the material prameters


  // testing parameters
  ft = 2.0e6; // TensileStrength
  lt = 0.8;   // TensileCharacteristicLength
  nt = 1.5;   // SofteningExponent
  fs = 3. * ft; // ShearStrengthRatio
  sigC0 = 190.0e6; // CompressiveStrength or CompressiveYieldingStrength
  hc0 = 0.4 * eN;  // InitialHardeningModulusRatio
  kc0 = 4.0;       // TransitionalStrainRatio
  ed = 1. * eN;    // DensificationRatio
  m0 = 0.2;        // InitialFriction
  mInf = 0.0;      // AsymptoticFriction
  sigN0 = 600.0e6; // TransitionalStress
  beta = 0.0;      // VolumetricDeviatoricCoupling
  kc1 = 1.0;       // DeviatoricStrainThresholdRatio
  kc2 = 5.0;       // DeviatoricDamageParameter

  ls = 1e20;
}

/**
  Function prints material parameters into the opened text file.

  @param out - pointer to the opened FILE

  @return The function does not return anything.

  7. 3. 2021, JK
*/
void cusatismaterial::print(FILE* out) {
  //  e - Young's modulus of elasticity
  //  nu - Poisson's ratio
  fprintf(out, "%le %le", e, nu);
}


/**
  Function assembles stiffness %matrix of material.

  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness %matrix in the parameter d.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::matstiff(matrix& d, strastrestate ssst) {
  switch (ssst) {
  case spacestress: {
    matstiff_spacestr(d);
    break;
  }
  default: {
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function assembles stiffness %matrix of material.

  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness matrix in the parameter d.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::elmatstiff(matrix& d, strastrestate ssst) {
  switch (ssst) {
  case spacestress: {
    matstiff_spacestr(d);
    break;
  }
  default: {
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 3D problems.

  @param d - stiffness %matrix of the material

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::matstiff_spacestr(matrix& d) {
  double g, s;

  nullm(d);

  g = e / 2.0 / (1.0 + nu);
  s = e / (1.0 + nu) / (1.0 - 2.0 * nu);

  d[0][0] = s * (1 - nu);  d[0][1] = s * nu;     d[0][2] = s * nu;
  d[1][0] = d[0][1];   d[1][1] = d[0][0];  d[1][2] = d[0][1];
  d[2][0] = d[0][1];   d[2][1] = d[0][1];  d[2][2] = d[0][0];

  d[3][3] = g;         d[4][4] = g;        d[5][5] = g;
}


/**
  Function assembles complience %matrix of material.

  @param c - complience %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::matcompl(matrix& c, strastrestate ssst) {
  switch (ssst) {
  case spacestress: {
    matcompl_spacestr(c);
    break;
  }
  default: {
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function creates compliance %matrix of the elastic
  isotropic material for 3D problems.

  @param c - compliance %matrix of the material (output)

  @return The function returns material complience %matrix in the parameter d.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::matcompl_spacestr(matrix& c) {
  double g;
  nullm(c);

  g = 2.0 * (1.0 + nu) / e;

  c[0][0] = 1.0 / e;     c[0][1] = -1.0 * nu / e;  c[0][2] = c[0][1];
  c[1][0] = c[0][1];   c[1][1] = c[0][0];    c[1][2] = c[0][1];
  c[2][0] = c[0][1];   c[2][1] = c[0][1];    c[2][2] = c[0][0];

  c[3][3] = g;         c[4][4] = g;          c[5][5] = g;
}

/* ELASTIC
void cusatismaterial::nlstresses(long ipp) {
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(n)), sig(ASTCKVEC(n)), epst(ASTCKVEC(6));
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(ASTCKMAT(n, n));

  //  initial values
  for (i = 0; i < n; i++) {
    eps[i] = Mm->ip[ipp].strain[i];
  }

  //matstiff (d,ssst);
  //mxv (d,eps,sig);

  double eT = alpha * eN;

  sig[0] = eN * eps[0];
  sig[1] = eT * eps[1];
  sig[2] = eT * eps[2];


  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)) {
    for (i = 0; i < n; i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }

  for (i = 0; i < n; i++) {
    Mm->ip[ipp].stress[i] = sig[i];
  }
}
*/


/**
  Function computes true stresses.

  @param ipp - number of integration point

  @return The function does not return anything.

  Created by JK, 7. 3. 2021
*/
void cusatismaterial::nlstresses(long ipp) {
  long i, n = Mm->ip[ipp].ncompstr;
  vector sig(ASTCKVEC(n)), epst(ASTCKVEC(6));

  // Inelastic behaviour (4.2)
  // Implementation according to Cusatis et al. "Lattice Discrete Particle Model (LDPM) for failure behavior of concrete. I: Theory"
  double dum, *a;
  double ht, hs, ftc, fsc, omega, h0, sig0, hc, eNc;
  double epsV, epsD, epsDV, epsC0, epsC1, sigC1, ele;

  double epsT, epsNMax, epsTMax, epsEqMax;
  double sigbt, sigbc, sigbs;

  double epsEq, epsEqOld;
  double epsN, epsL, epsM;
  double epsNOld, epsLOld, epsMOld;

  double sigN, sigL, sigM, sigEq;
  double sigNOld, sigLOld, sigMOld, sigEqOld;

  //  init values
  epsN = Mm->ip[ipp].strain[0];
  epsL = Mm->ip[ipp].strain[1];
  epsM = Mm->ip[ipp].strain[2];

  ele = Mm->ip[ipp].other[13]; // edge length
  epsV = Mm->ip[ipp].other[14]; // volumetric strain

  a = Mm->ip[ipp].other + 15;
  sigNOld = a[0];
  sigLOld = a[1];
  sigMOld = a[2];

  epsNOld = a[3];
  epsLOld = a[4];
  epsMOld = a[5];

  epsNMax = a[6];
  epsTMax = a[7];

  /* --- stress calculation --- */
  epsT = sqrt(epsL * epsL + epsM * epsM);
  epsEq = sqrt(epsN * epsN + alpha * epsT * epsT);
  epsEqOld = sqrt(epsNOld * epsNOld + alpha * (epsLOld * epsLOld + epsMOld * epsMOld));
  sigEqOld = sqrt(sigNOld * sigNOld + (sigLOld * sigLOld + sigMOld * sigMOld) / alpha);

  // max strain values
  if (epsNMax < epsN) {
    epsNMax = epsN;
  }
  if (epsTMax < epsT) {
    epsTMax = epsT;
  }
  epsEqMax = sqrt(epsNMax * epsNMax + alpha * epsTMax * epsTMax);
  
  // coupling variable
  dum = sqrt(alpha) * epsT;
  if (dum != 0.0) {
    omega = atan(epsN / dum);
  } else if (epsN < 0.0) { // pore collapse
    omega = -0.5 * M_PI;
  } else {                 // fracture
    omega = 0.5 * M_PI;
  }

  if (epsN > 0.0) { // fracturing behavior (Cusatis 4.2.1)
    double so = sin(omega);
    double co = cos(omega);
    double rstc;

    // post-peak slope (softening modulus) and tesile strength
    ht = 2. * eN; // modulus in pure tension
    ftc = ft;     // current tensile strength
    if (lt / ele <= 1.001) {
      ht *= 1000.0;
      ftc *= sqrt(lt / ele);  // according Bazant's: "Fracture and Size Effect in Concrete and Other Quasibrittle Materials", section 8.6.4 (linear sofening) 
    } else {
      ht /= (lt / ele - 1.);
    }
    hs = 2. * eT; // modulus in pure shear
    fsc = fs;     // current shear strength
    if (ls / ele <= 1.001) {
      hs *= 1000.0;
      fsc *= sqrt(ls / ele);  // according Bazant's: "Fracture and Size Effect in Concrete and Other Quasibrittle Materials", section 8.6.4 (linear sofening) 
    } else {
      hs /= (ls / ele - 1.);
    }
    rstc = fsc / ftc;

    // h0 = ht * pow(2. * omega / M_PI, nt);  // original formulation without softening in shear
    dum = hs / alpha;
    h0 = dum + (ht - dum) * pow(2. * omega / M_PI, nt);
    
    // strain-dependent boundary
    if (fabs(co) > Mp->zero) {
      dum = co * co / (rstc * rstc);
      sig0 = ftc * (-so + sqrt(so * so + 4. * alpha * dum)) / (2. * alpha * dum);
    } else {
      sig0 = ftc;
    }
    sigbt = sig0 * exp(-h0 * MAX(0.0, epsEqMax - sig0 / eN) / sig0);

    sigEq = MIN(MAX(0.0, sigEqOld + eN * (epsEq - epsEqOld)), sigbt);

    if (epsEq > Mp->zero) {
      dum = sigEq / epsEq;
      sigN = dum * epsN;
      sigL = dum * alpha * epsL;
      sigM = dum * alpha * epsM;
    } else {
      sigN = sigL = sigM = 0.0;
    }

  } else { // behaviour in compression
    
    // pore collapse and material compaction (Cusatis 4.2.2)
    epsD = epsN - epsV;
    epsDV = epsV + beta * epsD;
    eNc = eN;
    if (sigNOld < -sigC0) {
      eNc = ed;
    }

    hc = hc0;
    dum = epsD / epsV - kc1;
    if (dum > 0.0) {
      hc /= 1. + kc2 * dum;
    } 

    epsC0 = sigC0 / eN;
    epsC1 = kc0 * epsC0;
    sigC1 = sigC0 + (epsC1 - epsC0) * hc;
    if (-epsDV <= epsC0) {
      sigbc = -sigC0;
    } else if (-epsDV <= epsC1) {
      sigbc = -sigC0 + (epsDV + epsC0) * hc;
    } else {
      sigbc = -sigC1 * exp(-(epsDV + epsC1) * hc / sigC1);
    }

    sigN = MAX(MIN(0.0, sigNOld + eNc * (epsN - epsNOld)), sigbc);

    // frictional behaviour
    if (fabs(sigN0) > Mp->zero) {
      dum = (m0 - mInf) * sigN0;
      sigbs = fs + dum - mInf * sigN - dum * exp(sigN / sigN0);
    } else {
      sigbs = fs - mInf * sigN;
    }

    sigL = sigLOld + eT * (epsL - epsLOld); // elastic predictor
    sigM = sigMOld + eT * (epsM - epsMOld); // elastic predictor

    dum = sqrt(sigL * sigL + sigM * sigM);
    if (fabs(dum) > Mp->zero) {
      dum = MIN(dum, sigbs) / dum; // ratio sigT / sigT_elast

      sigL *= dum;
      sigM *= dum;
    }
  }


  // FINAL STRESS ADJUSTMENT AND SAVING
  sig[0] = sigN;
  sig[1] = sigL;
  sig[2] = sigM;

  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)) {
    for (i = 0; i < n; i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }

  for (i = 0; i < n; i++) {
    Mm->ip[ipp].stress[i] = sig[i];
  }

  a[0] = sigN;
  a[1] = sigL;
  a[2] = sigM;
  
  a[3] = epsN;
  a[4] = epsL;
  a[5] = epsM;
  
  a[6] = epsNMax;
  a[7] = epsTMax;
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  7. 3. 2021, JK
*/
void cusatismaterial::updateval(long ipp, long ido) {
  Mm->ip[ipp].eqother[ido+0] = Mm->ip[ipp].other[ido+0];     //
  Mm->ip[ipp].eqother[ido+1] = Mm->ip[ipp].other[ido+1];     //
}


/**
  The function changes material parameters for stochastic analysis.

  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters

  @return The function does not return anything.

  Created by JK,
*/
void cusatismaterial::changeparam(atsel& atm, vector& val) {
  long i;

  for (i = 0; i < atm.num; i++) {
    switch (atm.atrib[i]) {
    case 0: {
      e = val[i];
      break;
    }
    case 1: {
      nu = val[i];
      break;
    }
    default: {
      print_err("wrong number of atribute is required", __FILE__, __LINE__, __func__);
    }
    }
  }
}
