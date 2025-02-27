#include <string.h>
#include <stdlib.h>
#include "bocon.h"
#include "parser.h"
#include "globprep.h"


/**
  Constructor initializes data members
  Parameters :
    none
  Returns :
    nothing
*/
bocon::bocon()
{
  nn = 0L;
  ndir = 0L;
  dir = NULL;
  con = NULL;
  nspd = NULL;
  ndpd = NULL;
  gf = NULL;
}



/**
  Method reads data from the text file specified by parameter in, data are
    in preprocessor format
  Parameters :
    @param in    - pointer to the opened text file, where the data will be read
    @param ndofn - number of dofs in given node with prescribed condition
  Returns :
    @retval 0 - on success
    @retval 1 - number of boundary conditions is not in range <1, ndofn>
    @retval 2 - direction of a boundary condition is not in range <1, ndofn>
    @retval 3 - load case of a boundary condition is not in range <1, Nlc>
    @retval 4 - subload case of a boundary condition is not in range <1, Nslc>
    @retval 5 - error in parsing string for dynamic boundary condition
*/
long bocon::read(XFILE *in, long ndofn)
{
  long nbc, td, tlc, tslc, tnlc, i, sc, dc, gc, ret;
  double tc;
  gfunct tgf;

  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  if (tnlc == 0)
    tnlc = 1L;

  ndir = ndofn;

  // deallocate eventually used memory
  if (dir)
    delete [] dir;
  if (con)
    delete [] con;
  if (nspd)
    delete [] nspd;
  if (ndpd)
    delete [] ndpd;
  if (gf)
  {
    for (i = 0; i < tnlc*ndir; i++)
      delete gf[i];
    delete [] gf;
  }

  xfscanf (in, "%k%ld", "num_bc", &nbc);
  // reading number of conditions
  if ((nbc < 1) || (nbc > ndir))
  // which should be in range <1,ndofn>
    return(1);
  // allocating memory for arrays
  dir = new bctype[ndir];
  memset(dir, 0, sizeof(*dir)*ndir);
  nspd = new long[ndir];
  memset(nspd, 0, sizeof(*nspd)*ndir);
  ndpd = new long[ndir];
  memset(ndpd, 0, sizeof(*ndpd)*ndir);
  con = new double[ndir*tnlc];
  memset(con, 0, sizeof(*con)*ndir*tnlc);
  gf = new gfunct*[ndir*tnlc];
  memset(gf, 0, sizeof(*gf)*ndir*tnlc);

  // loop which reads conditions in particular directions
  for (i = 0; i < nbc; i++)
  {
    xfscanf(in, "%k%ld", "dir", &td);  // reading direction
    if ((td < 1) || (td > ndir))
    // direction should be in range <1, ndir>
      return(2);

    sc = dc = gc = 0L;
    sc = xfscanf(in, "%+k", "cond");
    if (sc == 0)
    {
      dc = xfscanf(in, "%+k", "tdcond");
      if (dc == 0)
        gc = xfscanf(in, "%+k", "gcond");
    }
    if (sc+dc+gc == 0)
    {
      print_err("cannot locate keyword 'cond', 'tdcond' nor 'gcond' in file %s, line %ld, col %ld.", 
                __FILE__, __LINE__, __func__, in->fname, in->line, in->col);
      abort();
    }
    if (sc+dc+gc > 1)
    {
      print_err("both keywords 'cond' and 'tdcond' or 'gcond' located in file %s, line %ld, col %ld.\n"
                "Correct type of boundary condition cannot be determined",
                __FILE__, __LINE__, __func__, in->fname, in->line, in->col);
      abort();
    }

    if (sc > 0) // static Dirichlet's condition
    {
      td--;
      dir[td] = sbc;
      tlc = 0;
      xfscanf(in, "%le", &tc); // reading static condition
      if (tc != 0.0) 
      // load case id is required for nonzero prescribed displacements
      {
        xfscanf (in, "%k%ld", "lc_id", &tlc);
        if ((tlc < 1) || (tlc > Nlc))
          return (3);
        tlc--;
        if (Tnslc)
        {
          xfscanf(in, "%k%ld", "slc_id", &tslc);
          if ((tslc < 1) || (tslc > Nslc[tlc]))
            return (4);
          tslc--;
          tlc = Nslc_cum[tlc]+tslc;
        }
        nspd[td]++; 
      }
      con[tlc*ndir+td] = tc;
    }
    if (gc > 0) // general static Dirichlet's condition
    {
      td--;
      dir[td] = gbc;
      tlc = 0;
      tgf.read(in); // reading static condition
      xfscanf (in, "%k%ld", "lc_id", &tlc);
      if ((tlc < 1) || (tlc > Nlc))
        return (3);
      tlc--;
      if (Tnslc)
      {
        xfscanf(in, "%k%ld", "slc_id", &tslc);
        if ((tslc < 1) || (tslc > Nslc[tlc]))
          return (4);
        tslc--;
        tlc = Nslc_cum[tlc]+tslc;
      }
      nspd[td]++;
      gf[tlc*ndir+td] = new gfunct; 
      gf[tlc*ndir+td]->copy(tgf);
    }

    // time dependent Dirichlet's condition
    if (dc > 0)
    {
      td--;
      dir[td] = tdbc;
      tlc = 0;
      ret = tgf.read(in); // reading time dependent condition
      if (ret)
        return (5);

      xfscanf (in, "%k%ld", "lc_id", &tlc);
      if ((tlc < 1) || (tlc > Nlc))
        return (3);

      tlc--;
      if (Tnslc)
      {
        xfscanf(in, "%k%ld", "slc_id", &tslc);
        if ((tslc < 1) || (tslc > Nslc[tlc]))
          return (4);
        tslc--;
        tlc = Nslc_cum[tlc]+tslc;
      }
      ndpd[td]++; 
      gf[tlc*ndir+td] = new gfunct; 
      gf[tlc*ndir+td]->copy(tgf);
      con[tlc*ndir+td] = 0.0;
    }
  }
  return(0);
}



/**
  The function initializes bocon structure for the given number of nodal DOFs.
  It is used for the periodic boundary condition for the transformation of
  master/slave condtions at corners of PUC.

  @param ndofn[in] - the number of DOFs  

  Created by Tomas Koudelka, 06.2021
*/
void bocon::init(long ndofn)
{
  long tnlc, i;
  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  if (tnlc == 0)
    tnlc = 1L;

  nn = 0L;
  ndir = ndofn;
  // allocating memory for arrays
  dir = new bctype[ndir];
  memset(dir, 0, sizeof(*dir)*ndir);
  nspd = new long[ndir];
  memset(nspd, 0, sizeof(*nspd)*ndir);
  ndpd = new long[ndir];
  memset(ndpd, 0, sizeof(*ndpd)*ndir);
  con = new double[ndir*tnlc];
  memset(con, 0, sizeof(*con)*ndir*tnlc);
  gf = new gfunct*[ndir*tnlc];
  memset(gf, 0, sizeof(*gf)*ndir*tnlc);
  for (i=0; i<ndir; i++){
    dir[i] = sbc;
    con[i] = 0.0;
  }
}



/**
  Method tries to merge two boundary conditions

  @param tbc  - object of boundary condition which will be copied

  @retval Method returns pointer to new object of boundary condition
*/
bocon *bocon::merge(bocon *tbc)
{
  bocon *ret = new bocon;
  long i, k, nlcid;

  ret->nn = nn;
  ret->ndir = ndir;
  nlcid = Nlc;
  if (Tnslc)
    nlcid = Tnslc;

  ret->dir = new bctype[ndir];
  memset(ret->dir, 0, sizeof(*ret->dir)*ndir);
  ret->nspd = new long[ndir];
  memset(ret->nspd, 0, sizeof(*ret->nspd)*ndir);
  ret->ndpd = new long[ndir];
  memset(ret->ndpd, 0, sizeof(*ret->ndpd)*ndir);
  ret->con = new double[ndir*nlcid];
  memset(ret->con, 0, sizeof(*ret->con)*ndir*nlcid);
  ret->gf = new gfunct*[ndir*nlcid];
  memset(ret->gf, 0, sizeof(*ret->gf)*ndir*nlcid);

  for (k=0; k < nlcid; k++)
  {
    for (i=0; i<ndir; i++)
    {
      // compare BC direction indicators
      if (dir[i] == tbc->dir[i])
      {      
        // both BC have prescribed the same type of condition in this direction
        ret->dir[i] = dir[i];
        // compare BC condition values
        switch (dir[i]){
          case nobc:
            // no BC prescribed => no other values have to be merged
            break;
          case sbc:
            if (con[ndir*k+i] == tbc->con[ndir*k+i]){
              // values of conditions are the same
              ret->con[ndir*k+i] = con[ndir*k+i];
              if (con[ndir*k+i] != 0.0)
                ret->nspd[i]++;
              break;
            }
            else{
              // conflict in merged BC
              // different values of boundary conditions in the same direction cannot be merged
              delete ret;
              return NULL;
            }
          case gbc:
            if (ret->gf[ndir*k+i]->compare(*gf[ndir*k+i]) == 0){
              ret->gf[ndir*k+i]->copy(*gf[ndir*k+i]);
              ret->nspd[i]++;
              break;
            }
            else{
              // conflict in merged BC
              delete ret;
              return NULL;
            }
          case tdbc:
            if (con[ndir*k+i] == tbc->con[ndir*k+i]){
              // values of conditions are the same
              ret->con[ndir*k+i] = con[ndir*k+i];
            }
            else{
              // conflict in merged BC
              delete ret;
              return NULL;
            }
            if (ret->gf[ndir*k+i]->compare(*gf[ndir*k+i]) == 0){
              // time functions are indentical
              ret->gf[ndir*k+i]->copy(*gf[ndir*k+i]);
              ret->ndpd[i]++;
              break;
            }
            else{
              // conflict in merged conditions - time functions are different
              delete ret;
              return NULL;
            }
          default:
            print_err("unknown type of boundary condition is required", __FILE__, __LINE__, __func__);
        }
      }
      else{
        // BC types for the given dir differ
        
        if((dir[i] == nobc) || (tbc->dir[i] == nobc)){
          // only one of the BC has been assigned to this direction
          if (tbc->dir[i] == nobc){
            // set dir of returned BC to this BC
            ret->dir[i] = dir[i];
            ret->con[ndir*k+i] = con[ndir*k+i];
            ret->nspd[i] = nspd[i];
            ret->ndpd[i] = ndpd[i];
            if (gf[ndir*k+i]){
              ret->gf[ndir*k+i] = new gfunct;
              ret->gf[ndir*k+i]->copy(*gf[ndir*k+i]);
            }
          }
          if (dir[i] == nobc){
            // set dir of returned BC to tbc
            ret->dir[i] = tbc->dir[i];
            ret->con[ndir*k+i] = tbc->con[ndir*k+i];
            ret->nspd[i] = tbc->nspd[i];
            ret->ndpd[i] = tbc->ndpd[i];
            if (tbc->gf[ndir*k+i]){
              ret->gf[ndir*k+i] = new gfunct;
              ret->gf[ndir*k+i]->copy(*tbc->gf[ndir*k+i]);
            }
          }
        }
        else{
          // conflict in merged BC types
          delete ret;
          return NULL;        
        }
      }
    }
  }
  return ret;
}



/**
  Method creates copy of given object 

  Returns:
  Method returns pointer to newly alocated object of boundary condition.
*/
bocon *bocon::copy(void)
{
  bocon *ret = new bocon;
  long i, k, nlcid;

  ret->nn = nn;
  ret->ndir = ndir;
  nlcid = Nlc;
  if (Tnslc)
    nlcid = Tnslc;
  if (nlcid == 0)
    nlcid = 1L;
  ret->dir = new bctype[ndir];
  memcpy(ret->dir, dir, sizeof(*dir)*ndir);
  ret->nspd = new long[ndir];
  memcpy(ret->nspd, nspd, sizeof(*nspd)*ndir);
  ret->ndpd = new long[ndir];
  memcpy(ret->ndpd, ndpd, sizeof(*ndpd)*ndir);
  ret->con = new double[ndir*nlcid];
  memcpy(ret->con, con, sizeof(*con)*ndir*nlcid);
  ret->gf = new gfunct*[ndir*nlcid];
  memset(ret->gf, 0, sizeof(*ret->gf)*ndir*nlcid);

  for (k=0; k < nlcid; k++)
  {
    for(i=0; i<ndir; i++)
    {
      if (gf[ndir*k+i])
      {
        ret->gf[ndir*k+i] = new gfunct;
        ret->gf[ndir*k+i]->copy(*gf[ndir*k+i]);
      }
    }
  }
  return ret;
}



bocon::~bocon()
{
  long tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;

  delete [] dir;
  delete [] con;
  delete [] nspd; 
  delete [] ndpd; 
  if (gf)
  {
    for (long i = 0; i < tnlc*ndir; i++)
      delete gf[i];
    delete [] gf;
  }
}
