#include "timeswmat.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "mechmat.h"



/**
  The constructor inializes attributes to zero values.

  Created by Tomas Koudelka, 21.11.2006
*/
timeswmat::timeswmat (void)
{
  ncm = 0;
  gf = NULL;
  ami = -1;
  nam = 0;
  store_eigstr = no;
}



/**
  The destructor is defined only for the formal purposes.

  Created by Tomas Koudelka, 21.11.2006
*/
timeswmat::~timeswmat (void)
{
  delete [] gf;
}



/**
  The function reads material data from the opened file in.

  @param in - pointer to the opened XFILE

  @return The function does not return anything.
   
  Created by Tomas Koudelka, 21.11.2006
*/
void timeswmat::read(XFILE* in)
{
  long i;
  xfscanf(in, "%k%ld", "ncm", &ncm);
  xfscanf(in, "%k%m", "store_eigstr", &answertype_kwdset, &store_eigstr);
  if (ncm < 2)
  {
    print_err("insufficient number of materials (%ld < 2) is required", __FILE__, __LINE__, __func__, ncm);
    abort();
  }
  gf = new gfunct[ncm];
  for (i=0; i < ncm; i++)
  {
    gf[i].read(in);
  } 
}



/**
  The function prints material data from the opened file in.

  @param out - pointer to the opened FILE

  @return The function does not return anything.
   
  Created by Tomas Koudelka, 21.11.2006
*/
void timeswmat::print(FILE* out)
{
  long i;
  fprintf(out, "%ld ", ncm);
  fprintf(out, "%d ", int(store_eigstr));
  for (i=0; i < ncm; i++)
    gf[i].print(out);
}



/**
  The function returns actual value of Young modulus depending on active material model.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns actual value of Young modulus.
   
  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
double timeswmat::give_actual_ym(long ipp, long /*im*/, long ido)
{
  double e = 0.0;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    e = Mm->give_actual_ym(ipp, ami, ido+ncompstr);      
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
  return e;
}
 


/**
  The function returns initial value of Young modulus depeneding on the active material model.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns initial value of Young modulus.
   
  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
double timeswmat::give_initial_ym(long ipp, long /*im*/, long ido)
{
  double e = 0.0;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    e = Mm->give_initial_ym(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
  return e;
}
 


/**
  The function returns actual value of Poissons ratio depeneding on the active material model.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns actual value of Poissons ratio.

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
double timeswmat::give_actual_nu(long ipp, long /*im*/, long ido)
{
  double nu = 0.0;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    nu = Mm->give_actual_nu(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
  return nu;
}
 


/**
  The function returns actual value of tensile strength depeneding on the active material model.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns actual value tensile strength.

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
double timeswmat::give_actual_ft(long ipp, long /*im*/, long ido)
{
  double ft = 0.0;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    ft = Mm->give_actual_ft(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
  return ft;
}



/**
  The function returns actual value of pore pressure depeneding on the active material model.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns actual value tensile strength.

  Created by Tomas Koudelka, 2.10.2013
  Rewritten by Tomas Koudelka, 9.2.2016
*/
/*
double timeswmat::give_pore_press(long ipp, long im, long ido)
{
  double p;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    p = Mm->give_pore_press(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on", __FILE__, __LINE__, __func__);
    abort();
  }
  return p;
}
*/


/**
  The function returns total number of components of other array.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns the total number of other components for the active material models
          or the maximum number of other components from all models in the model chain in the case
          that no model is active (initialization phase).

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
long timeswmat::givencompother (long ipp, long im)
{
  long ncompo = 0;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    ncompo = Mm->givencompother(ipp, ami);
  else 
  {   
    // return the maximum number of other components
    // for the initial phase when ami < 0
    long maxncompo = 0;
    for (long i=0; i<ncm; i++)
    {
      ncompo = Mm->givencompother(ipp, im+1+i);
      if (ncompo > maxncompo)
        maxncompo = ncompo;    
    }
    ncompo = maxncompo;
  }

  // add space for eigenstrains due to the stress free switching of the material models
  ncompo += ncompstr;

  return ncompo;
}



/**
  The function returns total number of components of eqother array.

  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns the total number of eqother components for the active material models
          or the maximum number of eqother components from all models in the model chain in the case
          that no model is active (initialization phase).

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
long timeswmat::givencompeqother (long ipp, long im)
{
  long ncompo = 0L;
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    ncompo = Mm->givencompeqother(ipp, ami); 
  else 
  {   
    // return the maximum number of eqother components
    // for the initial phase when ami < 0
    long maxncompo = 0;
    for (long i=0; i<ncm; i++)
    {
      ncompo = Mm->givencompeqother(ipp, im+1+i);
      if (ncompo > maxncompo)
        maxncompo = ncompo;    
    }
    ncompo = maxncompo;
  }

  // add space for eigenstrains due to the stress free switching of the material models
  ncompo += ncompstr;

  return ncompo;
}



/**
  The function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number
  @param im - material index
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function returns material stiffnes matrix in the parameter d

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
void timeswmat::matstiff (matrix &d,long ipp,long /*im*/, long ido)
{
  long ncompstr = 0L;
  vector epsb;

  if (store_eigstr == yes)
  {
    ncompstr = Mm->ip[ipp].ncompstr;
    reallocv(RSTCKVEC(ncompstr, epsb));
    // backup attained strains
    Mm->givestrain(0, ipp, epsb);
    // subtract stored eigenstrains attained at the material model activation
    for (long i=0; i<ncompstr; i++)
      Mm->ip[ipp].strain[i] -= Mm->ip[ipp].eqother[ido+i];
  }

  // compute the material stiffness matrix 
  if (ami > 0)
    Mm->matstiff(d, ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on  (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }

  // restore the total strains
  if (store_eigstr == yes)
    Mm->storestrain(0, ipp, epsb);
}


/**
  The function restores part of stress increments at the integration point ipp 
  to the %vector sig.
  
  @param lcid - load case id
  @param ipp - integration point number
  @param im  - index of the material in the tm and idm arrays on integration point
  @param ido - index of internal variables in the ip's ipp eqother array
  @param fi  - first index of the required stress increment component
  @param sig - %vector containing stress increment components (output)

  @return The function returns %vector of stress increments in the parameter sig.

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
void timeswmat::givestressincr (long lcid,long ipp,long /*im*/,long ido,long fi,vector &sig)
{
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    Mm->givestressincr(lcid, ipp, ami, ido+ncompstr, fi, sig);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
}


/**
  The function computes correct increments of stresses in the integration point
  for computed strains.
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function deos not return anything.

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
void timeswmat::nlstressesincr (long ipp, long /*im*/, long ido)
{
  long ncompstr = 0L;
  vector epsb;

  if (store_eigstr == yes)
  {
    ncompstr = Mm->ip[ipp].ncompstr;
    reallocv(RSTCKVEC(ncompstr, epsb));
    // backup attained strains
    Mm->givestrain(0, ipp, epsb);
    // subtract stored eigenstrains attained at the material model activation
    for (long i=0; i<ncompstr; i++)
      Mm->ip[ipp].strain[i] -= Mm->ip[ipp].eqother[ido+i];
  }

  // compute stress increments for the active material model
  if (ami > 0)
    Mm->computenlstressesincr(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }

  // restore the total strains
  if (store_eigstr == yes)
    Mm->storestrain(0, ipp, epsb);
}


/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.
  
  @param ipp - integration point pointer
  @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
  @param ido - index of internal variables for given material in the ipp eqother array

  @return The function deos not return anything.

  Created by Tomas Koudelka, 21.11.2006
  Rewritten by Tomas Koudelka, 9.2.2016
*/
void timeswmat::nlstresses (long ipp, long /*im*/, long ido)
{
  long ncompstr = 0L;
  vector epsb;

  if (store_eigstr == yes)
  {
    ncompstr = Mm->ip[ipp].ncompstr;
    reallocv(RSTCKVEC(ncompstr, epsb));
    // backup attained strains
    Mm->givestrain(0, ipp, epsb);
    // subtract stored eigenstrains attained at the material model activation
    for (long i=0; i<ncompstr; i++)
      Mm->ip[ipp].strain[i] -= Mm->ip[ipp].eqother[ido+i];
  }

  // compute stresses for the active material model
  if (ami > 0)
    Mm->computenlstresses(ipp,Mm->ip[ipp], ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }

  // restore the total strains
  if (store_eigstr == yes)
    Mm->storestrain(0, ipp, epsb);
}


/**
  The function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array

  @return The function deos not return anything.

  Created by Tomas Koudelka, 21.11.2006
*/
void timeswmat::updateval (long ipp, long /*im*/, long ido)
{
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    Mm->updateipvalmat(ipp, ami, ido+ncompstr);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp], ipp, Mp->time);
    abort();
  }
}



/**
  The function updates index of the first actual material, i.e.
  index of the actual control material model in the model chain.
  It is supposed to be called at time solver at the beginning of 
  each time step.

  @param im - index of material (i.e. index of given timeswmat in ip[].tm[] array)

  @retval  0 - in the case of no change in the model.
  @retval >0 - in the case of change of the active material model (ami is returned)
  @retval -1 - in the case that no material model is active.

  Created by Tomas Koudelka, 9.2.2016
*/
long timeswmat::update_ami(long im)
{
  long i, ret = -1;
  long j; // total index in the material chain

  for (i = 0, j=im+1; i < ncm; i++, j++)
  {
    if (gf[i].getval_long(Mp->time))
    {
      if (ami == j) // the same index of the first active model is being detected
        return 0;
      else 
      { // index of the first active model is different
        if (ret < 0)
        { // the first model is being detected
          ami = im+i+1; // im is the index of timeswmat in ip[].tm[] array, all governed material models follows => i+1
          nam = 1;
          ret = ami;
        }
        else
          nam++;
      }
    }  
  }
  if (ret < 0)
    ami = -1;

  return ret;
}



/**
  Function initializes eqother array with initial values.
  
  @param lcid - load case id   
  @param ipp - integration point pointer
  @param im - material index
  @param ido - index of internal variables for given material in the ipp other array
  @param[in] rinit - flag for initialization after restorage from hdbackup

  @return The function deos not return anything.

  Created by Tomas Koudelka, 21.11.2006
*/
void timeswmat::initvalues (long lcid, long ipp, long /*im*/, long ido, bool rinit)
{
  long ncompstr = 0L;

  if (store_eigstr == yes)
    ncompstr = Mm->ip[ipp].ncompstr;

  if (ami > 0)
    Mm->initvalues(lcid, ipp, ami, ido+ncompstr, rinit);
  else
  {
    print_err("no material is switched on (elem=%ld, ipp=%ld, time=%le)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp, Mp->time);
    abort();
  }
}
