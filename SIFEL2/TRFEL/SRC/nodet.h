#ifndef NODET_H
#define NODET_H

#include <stdio.h>
#include "iotools.h"
#include "aliast.h"
#include "vector.h"

/**
   class nodet defines node for transport problems
   it contains informations about nodes which are not stored in gnode
*/

class nodet
{
 public:
  nodet (void);
  ~nodet (void);
  void read (XFILE *in,long ndof);
  void print (FILE *out);

  void alloc_grad (long ncomp);
  void alloc_flux (long ncomp);
  void alloc_other (long ncomp0);
  void alloc_eqother (long ncomp0);
  void alloc_nodval ();
  void alloc_nodvali ();
  void alloc_nodvalp ();
  void alloc_nodvalt ();

  void nullgrad ();
  void nullflux ();
  void nullother ();
  void nulleqother ();
  void nullnodval ();
  void nullnodvalp ();
  void nullnodvali ();
  void nullnodvalt ();

  void storegrad (long fi,vector &gradv);
  void storegrad (long lcid,double vol,vector &fluxv);
  void storeflux (long fi,vector &fluxv);
  void storeflux (long lcid,double vol,vector &fluxv);
  void storeother (long ncomp,double otherv);
  void storeother (long ncomp,double vol,double otherv);
  void storeeqother (long ncomp,double eqotherv);
  void storeeqother (long ncompeq,double vol,double eqotherv);

  void   givegrad (long lcid, vector &gv);
  double givegrad (long lcid, long compid);
  void   giveflux (long lcid, vector &fv);
  double giveflux (long lcid, long compid);
  double giveother (long compid);
  void   giveother (vector &othv);
  double giveeqother (long compid);
  void   giveeqother (vector &eqothv);

  void grad_averageval ();
  void flux_averageval ();
  void other_averageval  ();
  void eqother_averageval  ();

  void save_nodval (double *nv);
  void give_values (double *in, double *inp, double *ineq );
  void save_values (double *out);
  void actual_previous_change ();

  void give_flux (long lcid,double *nv);

  ///  number of DOFs on node
  long ndofn;
  ///  number of components of gradients array
  long ncompgrad;
  ///  number of components of other array
  long ncompother;
  ///  number of components of eq_other array
  long ncompeqother;

  ///  type of cross section
  crsectypet crst;
  ///  number of appropriate cross section type
  long idcs;

  ///  array of contributions to the gradient array
  long *ncontr_grad;
  double *vol_grad;
  ///  array of contributions to the flux array
  long *ncontr_flux;
  double *vol_flux;
  ///  number of contributions to the other array
  long *ncontr_other;
  double *vol_other;
  ///  number of contributions to the eqother array
  long *ncontr_eqother;
  double *vol_eqother;

  ///  array containing gradients
  double **gradient;
  ///  array containing fluxes
  ///  flux[i][j]=k - the j-th component of the i-th flux is k
  double **flux;
  ///  array containing other values
  double *other;
  ///  array containing eqother values
  double *eqother;

  ///  array containing actual nodal values
  double *nodval;
  ///  array containing nodal values from the previous step
  double *nodvalp;
  ///  array containing initial nodal values
  double *nodvali;
  ///  array containing time derivative of nodal values
  double *nodvalt;
  /// master node weights for the scaling of contributions at hanging node
  vector *mnw;
};

#endif
