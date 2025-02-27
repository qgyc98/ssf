#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "layplate.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "elastisomat.h"
#include "strretalg.h"
#include "alias.h"
#include "element.h"

/**
  This constructor inializes attributes to zero values.
*/
layplate::layplate (void)
{
  nl = 0; nli = 0; err = 0.0;
  nm = NULL; idm = NULL;
}


/**
  The destructor
*/
layplate::~layplate (void)
{
  bcup.tm = NULL;
  bcup.idm = NULL;
  bcup.stress = NULL;
  bcup.strain = NULL;
  bcup.other = NULL;
  bcup.eqother = NULL;

  for (long i=0;i<nl;i++)
  {
    delete [] tm[i];
    delete [] idm[i];
  }
  delete [] tm;
  delete [] idm;
  delete [] nm;
}


/**
  This function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened text file

  10/2014   J. Fiedler
*/
void layplate::read (XFILE *in)
{
  long i,j;
//  long eid = Mm->elip[ipp];
  
  xfscanf (in,"%k%ld", "num_lay", &nl);
/*  if (Mc->give_num_lay(eid)!=nl){
    print_err("wrong number of layers defined", __FILE__, __LINE__, __func__);
  }
*/    
  tm   = new mattype*[nl];
  idm  = new long*[nl];
  nm   = new long[nl];
  
  for (i=0;i<nl;i++)
  {
    xfscanf (in,"%ld", &nm[i]);
    tm[i]  = new mattype[nm[i]];
    idm[i] = new long[nm[i]];

    for (j=0;j<nm[i];j++)
    {
      xfscanf(in, "%k%m%k%ld", "mattype", &mattype_kwdset, (int*)(tm[i]+j), "num_inst", idm[i]+j);
      idm[i][j]--;
    }
  }
  xfscanf (in,"%k%ld", "num_iter", &nli);
  xfscanf (in,"%k%lf", "error", &err);

  allocm(3,3,c);
  allocm(3,3,cn);
}



/**
  This function computes material stiffness %matrix considering layer stiffnesses.

  @param d[out]  - material stiffness %matrix
  @param ipp[in] - integration point number in the mechmat ip array.
  @param ido[in] - index of internal variables for given material in the ipp other array

  11/2014   J. Fiedler
*/
void layplate::matstiff (matrix &d, long ipp, long ido)
{
  long i,j,jj;
  double *th;
  double *zet;
  matrix e(3,3);
  matrix dn(3,3);
  matrix dnm(3,3);
  matrix dm(3,3);
  matrix pomat(3,3);
  matrix pommat(3,3);
  long eid = Mm->elip[ipp];
  double *k=NULL;

  th   = Mc->give_layer_thicke(eid);
  zet  = Mc->give_layer_zcoord(eid);
  
  fillm(0,e);
  fillm(0,d);
  backup(ipp,k,0);
  for (i=0;i<nl;i++){
    Mm->ip[ipp].tm  = tm [i];     // mattype aktualni vrstvy
    Mm->ip[ipp].idm = idm[i];     // index sady mat. vl. aktualni sady
    Mm->ip[ipp].nm  = nm [i];
  
    Mm->elmatstiff(e, ipp, ido);
    
    for (j=0;j<3;j++){
      for (jj=0;jj<3;jj++){
        dn[j][jj] += th[i]*e[j][jj];
      }
    }
    
    for (j=0;j<3;j++){
      for (jj=0;jj<3;jj++){
        dnm[j][jj] += zet[i]*th[i]*e[j][jj];
      }
    }
    
    for (j=0;j<3;j++){
      for (jj=0;jj<3;jj++){
        dm[j][jj] += th[i]*zet[i]*zet[i]*e[j][jj];
      }
    }
  }
  restore_values(ipp,k,0);
// stiffness matrix calculation
// D = -D_{nm}D^{-1}_{n}D_{nm}+D_{m}
  invm(dn,c,1.0e-20);
  cmulm(-1.0,dnm,pomat);
// matrix "c" is used for relation between eps and curavatures
// C = -D^{-1}_{n}D_{nm}
  mxm(c,pomat,cn);   //
  mxm(pomat,c,pommat);
  mxm(pommat,dnm,pomat);
  addm(pomat,dm,d);  // resulting stiffness matrix stored in "d"

}


/**
  This function computes internal forces of the plate at given integration point ipp,
  depending on the reached strains. The internal forces and the other attribute of
  given integration point are actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip.
  @param ido - index of internal variables for given material in the ipp other array.

  @return The function updates the stress array of the given integration point.

  1/2015 J. Fiedler
*/
void layplate::nlstresses (long ipp, long /*im*/, long ido)
{
  long i,j,nc = Mm->ip[ipp].ncompstr;
  long eid = Mm->elip[ipp];
  double *k;          // vector of curvatures
  double *layth;      // pointer to array of thicnkesses of layers
  double *layz;       // pointer to array of z-coordinates of layers
  vector intfor(nc);  // internal forces - bending and torque moments
  vector df(3);       // vector of unbalanced axial forces
  vector eps(3);      // deformation of middle plane of the plate
  vector dpeps(3);
  vector deps(3);
  double normf;

  
//  initial values
  Mm->ip[ipp].other[3]=1.0;
  layth = Mc->give_layer_thicke(eid);
  layz  = Mc->give_layer_zcoord(eid);
  backup (ipp,k,1);                      // backuping original values at given integration point

// initial deformations of middle plane of the plate
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      eps[i] += cn[i][j]*k[j];

// iteration loop for axial forces eqilibrium
  for (i=0;i<nli;i++)
  {
    stress_calc(df,eps,k,intfor,layz,layth,ipp,ido);
    normf = normv(df);
    if (normf<err)
      break;
    mxv(c,df,deps);
    for (j=0;j<3;j++)
      eps[j] -= deps[j];
    if (i==nli-1 && normf>err)
      {
        Mm->ip[ipp].other[3]=0.5;   // reducing the main loading step
      }
  }
    
  Mm->storestress(0, ipp, intfor);  // saving the internal forces
  restore_values (ipp,k,1);         // restoring original values at given integration point
}


/**
  This function evaluates the size of eqother array according to number of layers
  and material models used.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns the size of eqother array depending of number of layers
          and their material type.

  10/2014 J. Fiedler
*/
long layplate::compeqother(long ipp)
{
  long sizeeqother=0;
  long i;
  mattype *btm = Mm->ip[ipp].tm;
  Mm->ip[ipp].ncompstr = 4;

  // loop going through layers
  for (i=0;i<nl;i++){
    Mm->ip[ipp].tm  = tm [i];
    sizeeqother += Mm->givencompeqother(ipp,0)+3;
  }

  Mm->ip[ipp].tm  = btm;  // returning to original "layplate" mattype
  Mm->ip[ipp].ncompstr = 3;
  
  return sizeeqother+5;     // first 5 spots in other array is used for computation purposes
}


/**
  This function evaluates the size of other array according to number of layers
  and material models used.

  @param ipp - integration point number in the mechmat ip array.

  @return The function returns the size of other array depending of number of layers
          and their material type.

  10/2014 J. Fiedler
*/
long layplate::compother(long ipp)
{
  long sizeother=0;
  long i;
  mattype *btm = Mm->ip[ipp].tm;
  Mm->ip[ipp].ncompstr = 4;

  // loop going through layers
  for (i=0;i<nl;i++)
  {
    Mm->ip[ipp].tm  = tm[i];
    sizeother += Mm->givencompother(ipp,0)+3;
  }

  Mm->ip[ipp].tm  = btm;  // returning to original "layplate" mattype
  Mm->ip[ipp].ncompstr = 3;
  
  return sizeother+5;     // first 5 spots in other array is used for computation purposes
}


/**
  Function backs up some informations stored on integration point and prepares arrays for computation.

  @param ipp - integration point number in the mechmat ip array.

  10/2014 J. Fiedler
*/
void layplate::backup(long ipp, double *&k, long j)
{
  bcup.tm   = Mm->ip[ipp].tm;
  bcup.idm  = Mm->ip[ipp].idm;
  bcup.nm   = Mm->ip[ipp].nm;
  bcup.ssst = Mm->ip[ipp].ssst;
  bcup.ncompstr = Mm->ip[ipp].ncompstr;
  Mm->ip[ipp].ssst = planestress;
  Mm->ip[ipp].ncompstr = 4;

  if (j==1)
  {
    bcup.stress = Mm->ip[ipp].stress;
    k = Mm->ip[ipp].strain;
    Mm->ip[ipp].alloc_strains(1);
    Mm->ip[ipp].alloc_stresses(1);
    Mm->ip[ipp].clean_strains(1);
    Mm->ip[ipp].clean_stresses(1);
  }
}

/**
  Function restores original informations stored on integration point.

  @param ipp - integration point number in the mechmat ip array.

  10/2014 J. Fiedler
*/
void layplate::restore_values(long ipp, double *k, long j)
{
  long i;
  
  Mm->ip[ipp].tm       = bcup.tm;
  Mm->ip[ipp].idm      = bcup.idm;
  Mm->ip[ipp].nm       = bcup.nm;
  Mm->ip[ipp].ssst     = bcup.ssst;
  Mm->ip[ipp].ncompstr = bcup.ncompstr;

  if (j==1)
  {
    delete [] Mm->ip[ipp].strain;
    Mm->ip[ipp].strain = k;
    for (i=0; i<3; i++)
      bcup.stress[i] = Mm->ip[ipp].stress[i];
    delete [] Mm->ip[ipp].stress;
    Mm->ip[ipp].stress = bcup.stress;
  }
}

/**
  The function calculates stresses on all layers and resulting internal forces
  according to eps and curvatures of a plate.

  @param df     - vector of unbalanced normal forces.
  @param eps    - deformation of middle plane of the plate.
  @param k      - curvatures of the plate.
  @param intfor - vector of the internal forces.
  @param layz   - z-coordinates of the layers.
  @param layth  - thicnkesses of the layers.
  @param ipp    - integration point number in the mechmat ip array.
  @param ido    - index of internal variables for given material in the ipp other array

  @return The function returns calculated internal forces and vector of unbalanced normal forces.

  01/2015  J. Fiedler
*/
void layplate::stress_calc(vector &df, vector &eps, double *k, vector &intfor, double *layz, double *layth, long ipp, long ido)
{
  long i,j;
  long neq;
  
  for (i=0;i<3;i++){
    df[i] = 0.0;
    intfor[i] = 0.0;
  }
  
  ido = 5;
  Mm->ip[ipp].other[4]=0.0;

// main loop going through layers
  for (i=0;i<nl;i++)
  {
    Mm->ip[ipp].tm  = tm [i];     // mattype of the current layer
    Mm->ip[ipp].idm = idm[i];     // index of material atributes for the current layer
    Mm->ip[ipp].nm  = nm [i];

    // calculating deformation of the current layer
    for (j=0;j<3;j++)
      Mm->ip[ipp].strain[j] = k[j]*layz[i]+eps[j];

    // calculating stresses of the layer according to the defined material model
    Mm->computenlstresses(ipp,Mm->ip[ipp],0,ido);

    // adding up the internal forces
    for (j=0;j<3;j++)
      intfor[j] += Mm->ip[ipp].stress[j]*layth[i]*layz[i];

    // storing layer stresses to the "other" array    
    neq = Mm->givencompeqother(ipp,0);
    for (j=0;j<3;j++)
      Mm->ip[ipp].other[ido+neq+j] = Mm->ip[ipp].stress[j];

    // adding up the axial forces
    for (j=0;j<3;j++)
      df[j] += Mm->ip[ipp].stress[j]*layth[i];

    // adding up the gamma parameter
    for (j=0;j<3;j++)
      Mm->ip[ipp].other[4] += Mm->ip[ipp].other[ido+4];

    ido += neq+3;      // recalculation of "ido" for next layer
  }
}


/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip.
  @param ido - index of internal variables for given material in the ipp other array.

  @return The function does not return anything.
*/
void layplate::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
}


/**
  This function is used for reducing a step of the main nonlinear loading process.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip.
  @param ido - index of internal variables for given material in the ipp other array.

  @return The function returns a ratio of the step reduction.
*/
double layplate::dstep_red(long ipp, long /*im*/, long /*ido*/)
{
  return (Mm->ip[ipp].other[3]);
}

