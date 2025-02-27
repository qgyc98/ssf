#include "dbmatt.h"
#include "descript.h"
#include "iotools.h"
#include "transmat.h"
#include <string.h>



/**
   Constructor which initializes data members with null or zero values.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001

*/
matrt :: matrt(void)
{
  type = mattypet(-1);
  ninst = 0;
  inst = NULL;
  instu = NULL;
  ridx = NULL;
  ninstu = 0;
}



/**
   Destructor which deletes each string with material properties.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
matrt :: ~matrt(void)
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
  The function allocates memory for ni sets of material parameters that
  will be read as strings.

  @param ni - number of allocated instances

  @return The function does not return anything.

  Created by Tomas Koudelka, 30.6.2014
*/
void matrt::alloc(long ni)
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
  The function allocates memory for ni sets of material parameters that
  will be read/stored via transmat procedures.

  @param ni - number of allocated instances

  @return The function does not return anything.

  Created by Tomas Koudelka, 30.6.2014
*/
void matrt::alloctm(long ni)
{
  ninst = ni;
  instu = new long[ninst];
  memset(instu, 0, sizeof(*instu)*ninst);
  ridx = new long[ninst];
  memset(ridx, 0, sizeof(*ridx)*ninst);
}



/**
   This method reads one type of material with different parameters from the opened text file in.

   @param in - opened text file with material data.

   Returns:
     @retval 0 - on succes
     @retval 1 - error reading index of material type instance

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 7. 2009
*/
long matrt :: read(XFILE *in)
{
  long i, id;
  char emsg[200];

  xfscanf(in, "%k%m%k%ld", "mattype", &mattypet_kwdset, &type, "num_inst", &ninst);
  // allocating appropriate ammount of memory
  alloc(ninst);
  for (i=0; i < ninst; i++)
  //  Reading strings with material properties for each material index
  //  each line in the input file corresponds to one string with index of 
  //  instance and material properties for this instance
  {
    inst[i] = new char[20025];
    memset(inst[i], 0, sizeof(*inst[i])*20025);
    xfscanf(in, "%ld", &id);
    if ((id < 1) || (id > ninst))
    {
      sprintf(emsg, "index of the material type instance is out of range <1,%ld>\n"
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
dbmatt :: dbmatt(void)
{
  numt = 0;
  mat = NULL;
  matu = NULL;
  nmatu = 0L;
}



/**
   Destructor which deletes material database.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
dbmatt :: ~dbmatt(void)
{
  delete [] mat;
  delete [] matu;
}



/**
  The method reads material characteristics from the opened text file.
  The characteristic are stored as strings.

  @param in - pointer to the opened file

  @retval 0 - on success
  @retval i>0 - in case of failure in the reading of i-th material type

  Created by TKo, 6.2014
*/
long dbmatt::read(XFILE *in)
{
  // reading number of material types
  xfscanf(in, "%k%ld", "num_mat_types", &numt);
 
  mat = new matrt[numt]; // allocating memory for each material type
  matu = new long[numt]; 
  memset(matu, 0, sizeof(*matu)*numt);
  for (long i = 0; i < numt; i++)
  {
    // reading each index of material types
    if (mat[i].read(in))
    {
      return(i+1);
    }
  }
  return 0;
}



/**
  The method reads material characteristics from the opened text file
  with help of transmat procedures.

  @param in - pointer to the opened file
  @param mechm - pointer to the allocated mechmat structure
  @param d - pointer to structure with description of preprocessor setup 

  @retval 0 - on success
  @retval i>0 - in case of failure in the reading of i-th material type

  Created by TKo, 6.2014
*/
long dbmatt::readtm(XFILE *in, transmat *transm, descript *d)
{
  kwd_handling bmode;
  // reading number of material types
  xfscanf(in, "%k%ld", "num_mat_types", &numt);
 
  mat = new matrt[numt]; // allocating memory for each material type
  matu = new long[numt]; 
  memset(matu, 0, sizeof(*matu)*numt);
  transm->nmt = numt;
  transm->mattype = new mattypet [numt];  
  transm->numtype = new long [numt];
  for (long i = 0; i < numt; i++)
  {
    xfscanf(in, "%k%m %k%ld", "mattype", &mattypet_kwdset, &mat[i].type, "num_inst", &mat[i].ninst);
    mat[i].alloctm(mat[i].ninst);
    transm->mattype[i] = mat[i].type;
    transm->numtype[i] = mat[i].ninst;
    if (d->matkwd == no)
    {
      bmode = in->kwdmode;
      if (in->index_created)
        in->kwdmode = sect_mode_ignore;
      else
        in->kwdmode = ignore_kwd;
    }
    // reading each index of material types
    transm->readmattype(in, transm->mattype[i], transm->numtype[i]);
    if (d->matkwd == no)
      in->kwdmode = bmode;
  }
  return 0;
}



/**
   This method searches database whether contains material with
    type number ce and type index number ci.

   @param ce - material type number
   @param ci - material type index

   Returns :
    @retval If the material with given parameters exists in the database function returns
     index of given material in the array crs else returns -1.

   Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 11. 2001
*/
long dbmatt :: search_mat(mattypet ce, long ci)
{
  long i,ret = -1;
  for (i = 0; i < numt; i++)
  {
    if ((mat[i].type == ce) && ((ci >= 0) && (ci < mat[i].ninst)))
    {
      ret = i;
      break;
    }
  }
  return(ret);
}



/**
  The method marks of used material type and its instance(parameters).
  Eventually, it increases the number of used material types and instances(properties) of
  the given material type.

  @param im - index of material type in the database
  @param mti - index of material type instance(prarameters)

  Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 7. 2009
*/
void dbmatt::mark_used(long im, long mti)
{
  if (mat[im].instu[mti] == 0)
  {        
    mat[im].ninstu++; // number of used instances of material type
    mat[im].instu[mti]++; // indicator whether instance of material is used
    if (matu[im] == 0)
    {
      matu[im]++;
      nmatu++;
    }
  }
}



/**
  The method renumbers indeces of used instances of material types.
  
  Created by Tomas Koudelka koudelka@cml.fsv.cvut.cz, 7. 2009
*/
void dbmatt::renumber_id(void)
{
  long i, k, rindex;

  for (i = 0; i < numt; i++)
  { // for all material types in the database
    if (mat[i].ninstu) 
    { // some instances were marked as used in this material type
      rindex = 1;
      for (k = 0; k < mat[i].ninst; k++)
      {
        if (mat[i].instu[k])
        { // instance was used => assign new index starting from 1
          mat[i].ridx[k] = rindex;
          rindex++;
        }
      }
    }
  }
}
