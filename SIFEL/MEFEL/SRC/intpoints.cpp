#include <string.h>
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include "microM4.h"
#include "microSIM.h"
#include "microfiber.h"
#include "timeswmat.h"

#ifndef FNAMELEN
 #define FNAMELEN 1001
#endif



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
  Modified by Tomas Koudelka,
*/
intpoints::intpoints ()
{
  tm = NULL;  idm = NULL;
  nm = 0; hmt = 0;
  ncompstr=0;  ncompeqother=0;  ncompother=0;  ssst = (strastrestate) 0;
  strain = NULL;  stress = NULL;  other = NULL;  eqother = NULL;   nonloc = NULL;
}



/**
  Destructor releases allocated memory of the mechbclc object.

  Created by JK,
  Modified by Tomas Koudelka,
*/
intpoints::~intpoints (void)
{
  delete [] tm;      delete [] idm;
  delete [] strain;  delete [] stress;
  delete [] other;   delete [] eqother;
  delete [] nonloc;
}



/**
  The function reads material types and indices of used sets of material parameters on
  integration points from the opened text file. Actually, the function is not used and 
  the same material models on the integration points of one finite element are assumed.

  @param in - pointer to the opened XFILE
   
  @return The function does not return anything.

  Created by JK,
  Modified by Tomas Koudelka,
*/
void intpoints::read (XFILE *in)
{
  long i, idx, idxt;
  
  xfscanf (in, "%ld", &nm);
  tm  = new mattype[nm];
  idm = new long[nm];
  memset (tm,  0, sizeof(*tm)*nm);
  memset (idm, 0, sizeof(*idm)*nm);
  idx = 0;
  idxt = nm-1;
  for (i = 0; i < nm; i++)
  {
    xfscanf(in, "%k%m%k%ld", "mattype", &mattype_kwdset, (int*)(tm+idx), "num_inst", idm+idx);
    idm[idx]--;
    /* commented out due to timeswmat
      if (tm[idx] == therisodilat){
      tm[idxt] = tm[idx];
      idm[idxt] = idm[idx];
      hmt ^= 1;
      idxt--;
      idx--;
      }*/
    idx++;
  }
}



/**
  The function returns index of elastic material in the tm and idm arrays.
  
  @return The function returns index of elastic material model.

  Created by Tomas Koudelka, 28.5.2003
*/
long intpoints::gemid(void)
{
  long ret = nm - 1; // normally the elastic material is on the last position

  if (tm[0] == time_switchmat) // for time switched materials, the elastic model is the last one from the actual model chain
    ret = Mm->tswmat[idm[0]].ami + Mm->tswmat[idm[0]].nam - 1;

  if (hmt & 1) // the material has thermal dilatancy
    ret--;
  if (hmt & 4) // there is model of stress relaxation
    ret--;

  return ret;
}

/**
  The function returns index of relaxation material in the tm and idm arrays.
  
  @return The function returns index of material model of relaxation.

  Created by Tomas Koudelka, 13. 6. 2013
*/
long intpoints::grmid(void)
{
  long ret = -1;  //  usually, the model of relaxation is not present
  
  if (hmt & 4){  //  there is model of stress relaxation
    ret=nm-1;    //  the model of relaxation is on the last position for standard models

    // for time switched materials, the model of relaxation 
    // is on the last position of the active material chain
    if (tm[0] == time_switchmat) // for time switched materials, the elastic model is the last one from the actual model chain
      ret = Mm->tswmat[idm[0]].ami + Mm->tswmat[idm[0]].nam - 1;

    if (hmt & 1) //  if the model of thermal dilatation is defined, the relaxation is next to last position
      ret--;
  }
  
  return ret;
}



/**
  The function returns index of the nonlocal material in the tm and idm arrays.
  If no nonlocal material is present it returns -1.

  @return The function returns index of elstic material or -1.

  Created by Tomas Koudelka, 28.5.2003
*/
/*
long intpoints::gnlmid(void)
{
  long ret;

  if (hmt & 2) // the material is nonlocal
  {
    ret = nm - 1; // normally it is on the last position
    if (hmt & 1)  // the material has thermal dilatancy
      ret--;
  }
  else
    ret = -1;

  return ret;
}
*/



/**
  Function allocates memory for stress, strain, eqother, other and nonloc arrays
  for material with the index 0.
   
  @param nlc - number of load cases
  @param ipp - number of integration point
  @param ncompnl - number of nonloc array components

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void intpoints::alloc (long nlc,long ipp,long ncompnl)
{
  long ncompo = 0,ncompeqo = 0;
  
  //fprintf (Out,"\n ncompstr %ld  nlc %ld",ncompstr,nlc);

  //  number of components of array eqother
  ncompeqo = Mm->givencompeqother(ipp, 0);
  ncompeqother = ncompeqo;
  
  //  number of components of array other
  ncompo = Mm->givencompother(ipp, 0);
  ncompother = ncompo;
  
  strain = new double [nlc*ncompstr];
  stress = new double [nlc*ncompstr];
  memset(strain, 0, sizeof(*strain)*nlc*ncompstr);
  memset(stress, 0, sizeof(*stress)*nlc*ncompstr);
  
  
  if (ncompeqo > 0)
  {
    eqother = new double[ncompeqo];
    memset(eqother, 0, sizeof(*eqother)*ncompeqo);
  }
  if (ncompo > 0)
  {
    other   = new double[ncompo];
    memset(other, 0, sizeof(*other)*ncompo);
  }
  if (ncompnl > 0)
  {
    nonloc   = new double[ncompnl];
    memset(nonloc, 0, sizeof(*nonloc)*ncompnl);
  }
}



/**
  Function allocates array for strains.
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::alloc_strains (long nlc)
{
  strain = new double [nlc*ncompstr];
}



/**
  Function allocates array for stresses.
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::alloc_stresses (long nlc)
{
  stress = new double [nlc*ncompstr];
}



/**
  Function allocates array for other values.
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
/*
void intpoints::alloc_other (long ipp)
{
  
  //  number of components of array eqother
  ncompeqother = Mm->givencompeqother(ipp, 0);
  
  //  number of components of array other
  ncompother = Mm->givencompother(ipp, 0);
  
  if (ncompeqother > 0){
    eqother = new double[ncompeqother];
    memset(eqother, 0, sizeof(*eqother)*ncompeqother);
  }
  if (ncompother > 0){
    other   = new double[ncompother];
    memset(other, 0, sizeof(*other)*ncompother);
  }
  
}
*/


/**
  Copies values of stress/strain/other/eqother/nonloc arrays from the given integration point.

  @param ip - integration point whose values will be copied to the current instance
  @param nlc - number of load cases taken into account for allocation of stress/strain arrays
  @param ncompnl - number of components of nonloc array which will be copied
  @param realloc - inidicator for int. point arrays allocation
                   realloc=1 all arrays are reallocated
                   realloc=0 data are just copied, array must be allocated on the right dimensions

  @return The function does not return anything.

  Created by Tomas Koudelka, 29.11.2013 
*/
void intpoints::copy(intpoints &ip, long nlc, long ncompnl, long realloc)
{
  //  number of material types defined at integration point
  nm = ip.nm;

  // bit array with material type presence flags (i.e nonlocal model, thermal dilatancy)
  hmt = ip.hmt; 

  //  number of component of stress/strain array
  ncompstr = ip.ncompstr;

  //  number of component of eqother array
  ncompeqother = ip.ncompeqother;

  //  number of component of other array
  ncompother = ip.ncompother;

  //  stress-strain state
  ssst = ip.ssst;

  //  type of material
  if (realloc)
  {
    if (tm)
      delete [] tm;
    tm  = new mattype[nm];
  }
  memcpy (tm, ip.tm, sizeof(*tm)*nm);

  //  number of appropriate material type
  if (realloc)
  {
    if (idm)
      delete [] idm;
    idm = new long[nm];
  }
  memcpy (idm, ip.idm, sizeof(*idm)*nm);
  
  //  stress components
  if (realloc)
  {
    if (stress)
      delete [] stress;
    stress = new double[nlc*ncompstr];
  }
  memcpy (stress, ip.stress, sizeof(*stress)*ncompstr);

  //  strain components
  if (realloc)
  {
    if (strain)
      delete [] strain;
    strain = new double[nlc*ncompstr];
  }
  memcpy (strain, ip.strain, sizeof(*strain)*ncompstr);

  //  other components
  if (realloc)
  {
    if (other)
      delete [] other;
    if (ip.other)
      other = new double[ncompother];
    else
      other = NULL;
  }
  if(ip.other)
    memcpy (other, ip.other, sizeof(*other)*ncompother);

  //  equilibriated components of other array
  if (realloc)
  {
    if (eqother)
      delete [] eqother;
    if (ip.eqother)
      eqother = new double[ncompeqother];
    else
      eqother = NULL;
  }
  if (ip.eqother)
    memcpy (eqother, ip.eqother, sizeof(*eqother)*ncompeqother);

  //  nonlocal values
  if (realloc)
  {
    if (nonloc)
      delete [] nonloc;
    if (ip.nonloc)
      nonloc = new double[ncompnl];
    else
      nonloc = NULL;
  }
  if (ip.nonloc)
    memcpy (nonloc, ip.nonloc, sizeof(*nonloc)*ncompnl);
}



/**
  Function cleans all arrays defined at integration point
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::clean (long nlc)
{
  long i;
  for (i=0;i<nlc*ncompstr;i++){
    strain[i]=0.0;
    stress[i]=0.0;
  }
  for (i=0;i<ncompeqother;i++){
    eqother[i]=0.0;
  }
  for (i=0;i<ncompother;i++){
    other[i]=0.0;
  }
}



/**
  Function cleans strain array defined at integration point.
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::clean_strains (long nlc)
{
  long i;
  for (i=0;i<nlc*ncompstr;i++){
    strain[i]=0.0;
  }
}



/**
  Function cleans stress array defined at integration point
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::clean_stresses (long nlc)
{
  long i;
  for (i=0;i<nlc*ncompstr;i++){
    stress[i]=0.0;
  }
}



/**
  Function cleans arrays other and eqother defined at integration point.
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void intpoints::clean_other ()
{
  long i;
  for (i=0;i<ncompeqother;i++){
    eqother[i]=0.0;
  }
  for (i=0;i<ncompother;i++){
    other[i]=0.0;
  }
}



/**
  Function saves all arrays into the auxiliary file.
   
  @param aux - pointer to auxiliary file
  @param nlc - number of load cases
  @param selother - selection of components of eqother array which will be saved
   
  @return The function does not return anything.

  Created by JK, 19.9.2004
  Modified by Tomas Koudelka, 9.2008
*/
void intpoints::save_data_txt (FILE *aux,long nlc, sel &selother)
{
  long i;
  int prec = (int)Mp->hdbcont.prec;
  
  for (i=0;i<nlc*ncompstr;i++)
    fprintf (aux,"%.*le\n",prec,strain[i]);

  for (i=0;i<nlc*ncompstr;i++)
    fprintf (aux,"%.*le\n",prec,stress[i]);

  for (i=0;i<ncompeqother;i++)
  {
    if (selother.presence_id(i))
      fprintf (aux,"%.*le\n",prec,eqother[i]);
  }
}


/**
  Function restores all arrays from the auxiliary file.
   
  @param aux - pointer to auxiliary file
  @param nlc - number of load cases
  @param ncompo - number of saved eqother components
  @param selother - selection of components of eqother array to which will be other values restored
  @param selid - array of indices 
   
  @return The function does not return anything.

  Created by JK, 19.9.2004
  Modified by Tomas Koudelka, 9.2008
*/
void intpoints::restore_data_txt (FILE *aux,long nlc, long ncompo, sel &selother, long *selid)
{
  long i, ik, is;
  double tmp;
  
  for (i=0;i<nlc*ncompstr;i++)
    fscanf (aux,"%le",strain+i);

  for (i=0;i<nlc*ncompstr;i++)
    fscanf (aux,"%le",stress+i);

  for (i=0;i<ncompo;i++)
  {
    fscanf (aux,"%le",&tmp);
    if (selother.presence_id(i,ik))
    {
      is = selid[ik]+i-selother.id1[ik]; // storage index in eqother array 
      if (is >= ncompeqother)
        print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
      else
        eqother[is] = tmp;
    }
  }
}



/**
  Function saves all arrays into the auxiliary file in binary format.
   
  @param aux - pointer to auxiliary file
  @param nlc - number of load cases
  @param selother - selection of components of eqother array which will be saved
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void intpoints::save_data_bin (FILE *aux,long nlc, sel &selother)
{
  long i;
  
  fwrite (strain, sizeof(*strain),nlc*ncompstr, aux);
  fwrite (stress, sizeof(*stress),nlc*ncompstr, aux);
  for (i=0;i<ncompeqother;i++)
  {
    if (selother.presence_id(i))
      fwrite (eqother+i, sizeof(*eqother), 1, aux);
  }
}



/**
  Function restores all arrays from the auxiliary file in binary format.
   
  @param aux - pointer to auxiliary file
  @param nlc - number of load cases
  @param selother - selection of components of eqother array to which will be other values restored
  @param selid    - array of indices of positions in eqother array to which will be 
                    restored selected saved eqother components
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void intpoints::restore_data_bin (FILE *aux,long nlc, long ncompo, sel &selother, long *selid)
{
  long i, ik, is;
  double tmp;
  
  fread (strain, sizeof(*strain), nlc*ncompstr, aux);
  fread (stress, sizeof(*stress), nlc*ncompstr, aux);
  for (i=0;i<ncompo;i++)
  {
    fread (&tmp, sizeof(tmp), 1, aux);
    if (selother.presence_id(i,ik))
    {
      is = selid[ik]+i-selother.id1[ik]; // storage index in eqother array 
      if (is >= ncompeqother)
        print_err("invalid index for eqother restoring is required", __FILE__, __LINE__, __func__);
      else
        eqother[is] = tmp;
    }
  }
}

/**
   function returns the number of components of strain/stress vectors
   
   this function is important for plane stress and plane strain
   the arrays strain and stress contains 4 components in the code
   but user expects 3 components only
   
   11. 7. 2014, JK
*/
long intpoints::give_ncompstr ()
{
  if (ssst==planestress || ssst==planestrain)
    return 3;
  else
    return ncompstr;
}
