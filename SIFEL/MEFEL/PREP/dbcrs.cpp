#include "dbcrs.h"
#include "descrip.h"
#include "alias.h"
#include "iotools.h"
#include "mechcrsec.h"
#include <string.h>



/**
   Constructor which initializes data members with null or zero values.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001

*/
csec :: csec(void)
{
  type   = crsectype(-1);
  ninst  = 0;
  inst   = NULL;
  instu  = NULL;
  ridx   = NULL;
  ninstu = 0;
}



/**
   Destructor which deletes each string with cross-section properties.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
csec :: ~csec(void)
{
  long i;
  if (inst)
  {
    for (i=0; i < ninst; i++)
      delete [] inst[i];
  }
  delete [] inst;
  delete [] instu;
  delete [] ridx;
}



/**
  The function allocates memory for ni sets of cross section parameters that
  will be read as strings.

  @param ni - number of allocated instances

  @return The function does not return anything.

  Created by Tomas Koudelka, 1.7.2014
*/
void csec::alloc(long ni)
{
  ninst = ni;
  inst = new char*[ninst];
  memset(inst, 0, sizeof(*inst)*ninst);
  instu = new long[ninst];
  memset(instu, 0, sizeof(*instu)*ninst);
  ridx = new long[ninst];
  memset(ridx, 0, sizeof(*ridx)*ninst);
}



/**
  The function allocates memory for ni sets of cross section parameters that
  will be read/stored via mechcrsec procedures.

  @param ni - number of allocated instances

  @return The function does not return anything.

  Created by Tomas Koudelka, 1.7.2014
*/
void csec::allocmc(long ni)
{
  ninst = ni;
  instu = new long[ninst];
  memset(instu, 0, sizeof(*instu)*ninst);
  ridx = new long[ninst];
  memset(ridx, 0, sizeof(*ridx)*ninst);
}



/**
   This method reads one type crosssection with different parameters from the opened text file in.
   The parameters are stored as strings.

   @param in - opened text file with cross-section data.

   @retval 0 - on succes
   @retval 1 - error reading cross-section type and number of different parametre sets of given cross-section type
   @retval 2 - error reading cross-section properties string

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
long csec::read(XFILE *in)
{
  long i, id;
  char emsg[1001];

  xfscanf(in, "%k%m%k%ld", "crstype", &crsectype_kwdset, &type, "num_inst", &ninst);
  // allocating apropriate ammount of memory
  alloc(ninst);
  for (i=0; i < ninst; i++)
  //  Reading strings with cross section properties for each cross section index
  //  each line in the input file corresponds to one string with index of 
  //  instance and cross section properties for this instance
  {
    inst[i] = new char[1025];
    memset(inst[i], 0, sizeof(*inst[i])*1025);
    xfscanf(in, "%ld", &id);
    if ((id < 1) || (id > ninst))
    {
      sprintf(emsg, "index of the cross-section type instance is out of range <1,%ld>\n"
                    "input file : %s, line=%ld, column=%ld", ninst, in->fname, in->line, in->col);
      print_err(emsg, __FILE__, __LINE__, __func__);
      return 1;
    }
    else
    {
      xfscanf(in,"%a", inst[id-1]);
      skipline(in);
    }
  }
  return(0);
}



/**
   Constructor which initializes data members with null or zero values.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001

*/
dbcrs::dbcrs(void)
{
  numt = 0L;
  crs = NULL;
  crsu = NULL;
  ncrsu = 0L;
}



/**
   Destructor which deletes cross-section database.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
dbcrs::~dbcrs(void)
{
  delete [] crs;
  delete [] crsu;
}



/**
  The method reads cross section characteristics from the opened text file.
  The parameters are stored as strings

  @param in - pointer to the opened file

  @retval 0 - on success
  @retval i>0 - in case of failure in the reading of i-th cross section type

  Created by TKo, 6.2014
*/
long dbcrs::read(XFILE *in)
{
  // reading number of cross section types
  xfscanf(in, "%k%ld", "num_crsec_types", &numt);
 
  crs = new csec[numt]; // allocating memory for each cross section type
  crsu = new long[numt]; 
  memset(crsu, 0, sizeof(*crsu)*numt);
  for (long i = 0; i < numt; i++)
  {
    // reading each index of cross section types
    if (crs[i].read(in))
    {
      return(i+1);
    }
  }
  return 0;
}



/**
  The method reads cross section characteristics from the opened text file
  with help of mechcrsec procedures.

  @param in - pointer to the opened file
  @param mechc - pointer to the allocated mechcrsec structure
  @param d - pointer to structure with description of preprocessor setup 

  @retval 0 - on success
  @retval i>0 - in case of failure in the reading of i-th cross section type

  Created by TKo, 6.2014
*/
long dbcrs::readmc(XFILE *in, mechcrsec *mechc, descrip *d)
{
  kwd_handling bmode;
  // reading number of cross section types
  xfscanf(in, "%k%ld", "num_crsec_types", &numt);
 
  crs = new csec[numt]; // allocating memory for each cross section type
  crsu = new long[numt]; 
  memset(crsu, 0, sizeof(*crsu)*numt);
  mechc->ncst = numt;
  mechc->cstype  = new crsectype [numt];  
  mechc->numtype = new long [numt];
  for (long i = 0; i < numt; i++)
  {
    xfscanf(in, "%k%m %k%ld", "crstype", &crsectype_kwdset, &crs[i].type, "num_inst", &crs[i].ninst);
    crs[i].allocmc(crs[i].ninst);
    mechc->cstype[i] = crs[i].type;
    mechc->numtype[i] = crs[i].ninst;
    if (d->crskwd == no)
    {
      bmode = in->kwdmode;
      if (in->index_created)
        in->kwdmode = sect_mode_ignore;
      else
        in->kwdmode = ignore_kwd;
    }
    // reading each index of cross section types
    mechc->readcrsectype(in, mechc->cstype[i], mechc->numtype[i]);
    if (d->crskwd == no)
      in->kwdmode = bmode;
  }
  return 0;
}



/**
   This method searches database whether contains cross-section with
    type number ce and type index number ci.

   @param ce - cross-section type number
   @param ci - cross-section type index

   Returns :
    @retval If the cross-section with given parameters exists in the database function returns
     index of given cross-section in the array crs else returns -1.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
long dbcrs::search_crs(crsectype ce, long ci)
{
  long i,ret = -1;

  for (i = 0; i < numt; i++)
  {
    if ((crs[i].type == ce) && ((ci >= 0) && (ci < crs[i].ninst)))
    {
      ret = i;
      break;
    }
  }
  return(ret);
}



/**
  The method marks of used cross-section type and its instance(parameters).
  Eventually, it increases the number of used cs types and instances(property) of
  the given cs type.

  @param ic - index of cross-section type in the database
  @param csti - index of cross-section type instance(prarameters)

  Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 7. 2009
*/
void dbcrs::mark_used(long ic, long csti)
{
  if (crs[ic].instu[csti] == 0)
  {        
    crs[ic].ninstu++; // number of used instances of cs type
    crs[ic].instu[csti]++; // indicator whether instance of cs is used
    if (crsu[ic] == 0)
    {
      crsu[ic]++;
      ncrsu++;
    }
  }
}



/**
  The method renumbers indeces of used instances of cross-section types.
  
  Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 7. 2009
*/
void dbcrs::renumber_id(void)
{
  long i, k, rindex;

  for (i = 0; i < numt; i++)
  { // for all cross-section types in the database
    if (crs[i].ninstu) 
    { // some instances were marked as used in this cs type
      rindex = 1;
      for (k = 0; k < crs[i].ninst; k++)
      {
        if (crs[i].instu[k])
        { // instance was used => assign new index starting from 1
          crs[i].ridx[k] = rindex;
          rindex++;
        }
      }
    }
  }
}
