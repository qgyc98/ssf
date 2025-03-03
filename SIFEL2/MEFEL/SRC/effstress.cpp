#include "effstress.h"

#include "galias.h"
#include "matrix.h"
#include "vector.h"
#include "gfunct.h"
#include "iotools.h"

#include "alias.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>



effstress::effstress()
{
  ppv = 0.0;
  srv = 0.0;
  varppv = NULL;
  varsrv = NULL;
}



effstress::~effstress()
{
  delete varppv;
  delete varsrv;
}



/**
  The function reads type of pore pressure computation and in the 
  case of contant for all problem solved, the value is read from 
  the input file.
 
  @param in - pointer to the opened XFILE

  @return The funtion does not return anything.

  Created by Tomas Koudelka, 2.10.2013
*/
void effstress::read(XFILE *in)
{
  xfscanf(in, "%k%m", "pore_press_comp", &pore_press_comp_kwdset, &ppct);
  switch (ppct)
  {
    case noporep:
      break;
    case const_all:
      xfscanf(in, "%le%le", &ppv, &srv);
      if (Mp->pore_press < 2)
        Mp->pore_press = 1; // prescribed pore presures are constant in time, no data transfer from TRFEL
      else
      {
        print_err("uncoupled and partially coupled approaches cannot be used together", __FILE__, __LINE__, __func__);
        abort();
      }
      break;
    case nonmechq:
      if (Mp->pore_press == 1)
      {
        print_err("uncoupled and partially coupled approaches cannot be used together", __FILE__, __LINE__, __func__);
        abort();
      }
      else
      {
        // partially coupled approach is supposed
        // it can be rewritten in METR according to METR input file
        // pore pressure values are taken from TRFEL
        Mp->pore_press = 2; 
      }
      break;
    case var_press:
      varppv = new gfunct;
      varsrv = new gfunct;
      varppv->read(in);
      varsrv->read(in);
      if (Mp->pore_press < 2)
        Mp->pore_press = 1; // prescribed pore presures are given by general function (variation in time and space can be given),
                            // no data transfer from TRFEL
      else
      {
        print_err("uncoupled and partially coupled approaches cannot be used together", __FILE__, __LINE__, __func__);
        abort();
      }
      break;
    default:
      print_err("Unknown type of pore pressure computation is required", __FILE__, __LINE__, __func__);
  }

  return;
}



/**
  The function calculates stresses according to effective stress concept.

  @param ipp - number of integration point
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables in the ip's ipp eqother array. The standard value is zero.

  @return The function does not return anything but fills stress array on the given integration point  
          with the actual stress values.

  Created by TKo & TKr, 09.2012
  Modified by TKo, 10.2013
*/
void effstress::nlstresses (long ipp, long im, long ido)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;
  vector uw(ASTCKVEC(ncomp));
  matrix uwt(ASTCKMAT(3,3));
  vector sig(ASTCKVEC(ncomp));
  matrix sigt(ASTCKMAT(3,3));
  double pp, sr;
  
  // restore initial value of effective pore pressure
  pp = Mm->ip[ipp].eqother[ido];

  // create tensor form of effective pore pressure  -> uwt
  for (i=0; i<3; i++)
    uwt(i,i) = pp;
  
  // convert effective pore pressure tensor back to vector form -> uw
  tensor_vector(uw, uwt, Mm->ip[ipp].ssst,stress);
  
  // compute effective eigenstress values because the total values have to be given in the input file
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (i=0; i<ncomp; i++)
      Mm->eigstresses[ipp][i] -= uw(i);
  }
  
  // eventually, actualize pore pressure values 
  if (ppct == var_press)
  {
    pp = varppv->getval(Mp->time);
    sr = varsrv->getval(Mp->time);
    // set effective_pore_pressure_component of nonmechq array to constant value, saturation degree is assued as the effective stress factor chi
    Mm->storenonmechq(eff_pore_pressure, ipp, pp*sr);
  }

  // compute actual values of effective stresses
  Mm->computenlstresses(ipp, Mm->ip[ipp], im+1, ido+1);//eqother values for other materials starts at 2 index

  // retrieve actual effective stress vector -> sig
  Mm->givestress(0, ipp, sig);

  // create tensor form of effective stress vector -> sigt
  vector_tensor(sig, sigt, Mm->ip[ipp].ssst,stress);
  
  // actual effective pore pressure
  pp = Mm->givenonmechq(eff_pore_pressure,ipp);

  // calculate total stress tensor
  // sig^{tot}_{ij} = sig^{'}_{ij} + p*\delta_{ij}
  for (i=0;i<3;i++)
    sigt(i,i) += pp;

  // convert total stress tensor back to vector form -> sig
  tensor_vector(sig, sigt, Mm->ip[ipp].ssst,stress);
  
  // store the total stresses 
  for (i=0;i<ncomp;i++)
    Mm->ip[ipp].stress[i] = sig[i];    

  // recover total eigenstress values
  // (it may not be necessary because total values of eigenstresses should be recovered in
  //  the course of right hand side assembling where the original values given in the input file are used)
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (i=0; i<ncomp; i++)
      Mm->eigstresses[ipp][i] += uw(i);
  }
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
*/
void effstress::give_reqnmq(long *anmq)
{
  switch (ppct)
  {
    case noporep:
      break;
    case const_all:
    case nonmechq:
      anmq[eff_pore_pressure-1] = 1;
      break;
    default:
      print_err("Unknown type of pore pressure computation is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function initializes corresponding component of the mechmat::nonmechq array
  in the case of constant pore pressure value defined in the material.
  Additionally, it sets flag probdesc::eigstrcomp for computation of nodal force 
  contribution due to pore pressure.

  @param[in] lcid - load case id
  @param[in] ipp  - integration point pointer
  @param[in] im   - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param[in] ido  - index of internal variables in the ip's ipp eqother array. The standard value is zero.
  @param[in] rinit - flag for initialization after restorage from hdbackup

  @return The function does not return anything but it changes Mm->nonmechq array.

  Created by Tomas Koudelka, 14.10.2013 
*/
void effstress::initvalues(long lcid, long ipp, long im, long ido, bool rinit)
{
  long i,ncomp = Mm->ip[ipp].ncompstr;
  vector uw(ASTCKVEC(ncomp));
  matrix uwt(ASTCKMAT(3,3));
  double pp, sr;

  if (ppct == const_all)
  {
    // set pore_pressure_component of nonmechq array to constant value
    Mm->storenonmechq(eff_pore_pressure, ipp, ppv*srv);   
  }
  if (ppct == var_press)
  {
    pp = varppv->getval(Mp->time);
    sr = varsrv->getval(Mp->time);
    // set pore_pressure_component of nonmechq array to constant value
    Mm->storenonmechq(eff_pore_pressure, ipp, pp*sr);
  }

  // actual effective pore pressure  
  pp = Mm->givenonmechq(eff_pore_pressure,ipp);

  // store initial value of effective pore pressure
  Mm->ip[ipp].other[ido] = Mm->ip[ipp].eqother[ido] = pp;


  // compute effective eigenstress values because the total values have to be given in the input file
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    // calculate pore pressure tensor -> uwt
    for (i=0; i<3; i++)
      uwt(i,i) += pp;
    
    // convert total pore pressure tensor back to vector form -> uw
    tensor_vector(uw, uwt, Mm->ip[ipp].ssst,stress);

    // calculate effective eigenstresses
    for (i=0; i<ncomp; i++)
      Mm->eigstresses[ipp][i] -= uw(i);
  }
    
  // initialization of other material models from the model sequence -> stress array on integration points
  // may be changed
  Mm->initvalues(lcid, ipp, im+1, ido+1, rinit);//eqother values for other materials starts at 2 index

  //initial value of stresses
  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].stress[i] += uw(i);

  // recover total eigenstress values
  // (it may not be necessary because total values of eigenstresses should be recovered in
  //  the course of right hand side assembling where the original values given in the input file are used)
  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
    for (i=0; i<ncomp; i++)
      Mm->eigstresses[ipp][i] += uw(i);
  }
}
