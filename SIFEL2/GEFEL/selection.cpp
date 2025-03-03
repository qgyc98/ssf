#include "selection.h"
#include "iotools.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



/**
  Constructor initializes data to zero values.
 */
sel::sel()
{
  st = sel_no;
  n = 0;
  ncomp = id1 = NULL;
  rid1 = rid2 = NULL;
  sarray = NULL;
  err = 0.0;
  initime = fintime = r = 0.0;
  ent = NULL;
}


/**
  Destructor deallocates used memory.
 */
sel::~sel()
{
  delete [] ncomp;
  delete [] id1;
  delete [] rid1;
  delete [] rid2;
  delete [] sarray;
  delete [] ent;
}



/**
  Function reads data of selections from a text file given by parameter in.
  
  @param in - pointer to opened text file

  @retval 0 - on success
  @retval 1 - error selection type
  @retval 2 - error reading selection values
  @retval 3 - error reading range or list
  @retval 4 - error reading real range or real list

  Created by Tomas Koudelka, 2004
 */
long sel::read(XFILE *in)
{
  long i;
  const char *kwdstr = NULL;

  if (xfscanf(in, "%m", &seltype_kwdset, (int *)&st) != 1)
  {
    print_err("cannot read type of selection", __FILE__, __LINE__, __func__);
    return 1;
  }
  switch (st)
  {
    case sel_no:
      n = 0;
      return 0;
    case sel_all:
    case sel_impvalues:
      n = 1;
      id1 = new long [1];
      ncomp = new long [1];
      id1[0] = 0;
      ncomp[0] = 0;
      return 0;
    case sel_period:
      n = 1;
      id1 = new long [1];
      ncomp = new long [1];
      ncomp[0] = 1;
      if (xfscanf(in, "%ld", id1) != 1)
      {
        print_err("cannot read period", __FILE__, __LINE__, __func__, seltype_kwdset.get_str(st));
        return 2;
      }
      return 0;
    case sel_single:
      n = 1;
      id1 = new long [1];
      ncomp = new long [1];
      ncomp[0] = 1;
      if (xfscanf(in, "%ld", id1) != 1)
      {
        print_err("cannot read single selection", __FILE__, __LINE__, __func__, seltype_kwdset.get_str(st));
        return 2;
      }
      id1[0]--;
      return 0;
    case sel_realperiod:
      if (xfscanf(in, "%k%le%k%le%k%le%k%le", "ini_time", &initime,
                  "fin_time", &fintime, "period", &r, "err", &err) != 8)
      {
        print_err("Cannot read real period selection", __FILE__, __LINE__, __func__);
        return 2;
      }
      if ((fintime - initime) < 0.0)
      {
        print_err("Final time is less than initial time", __FILE__, __LINE__, __func__);
        return 2;
      }
      n = long((fintime - initime)/r) + 1;
      if (n < 2)
        n = 2;
      rid1 = new double[n];
      memset(rid1, 0, sizeof(*rid1)*n);
      rid1[0] = initime;
      for (i=1; i<n-1; i++)
        rid1[i] += rid1[i-1] + r;
      rid1[n-1] = fintime; 
      return 0;
    case sel_mtx:  
      n = 1;
      id1 = new long [1];
      ncomp = new long [1];
      id1[0] = 0;
      ncomp[0] = 0;
      return 0;
    case sel_vec:
    case sel_range_vec:
    case sel_range_mtx:
      n = 1;
      break;
    case sel_range:
    case sel_realrange:
      kwdstr = "numranges";
      break;
    case sel_list:
    case sel_reallist:
    case sel_impvallst:
      kwdstr = "numlist_items";
      break;
    case sel_prop:
      if (xfscanf(in, "%k%ld", "numprop", &n) != 2)
      {
        print_err("cannot read number of selected property id", __FILE__, __LINE__, __func__);
        return 2;
      }
      id1 = new long [n];
      ncomp = new long [n];
      ent = new gentity [n];
      for (i=0; i<n; i++) 
      {
        ncomp[i] = 1;
        if (xfscanf(in, "%k%ld", "prop", id1+i) != 2)
        {
          print_err("cannot read property number", __FILE__, __LINE__, __func__);
          return 2;
        }
        if (xfscanf(in, "%k%m", "ent", &gentity_kwdset, ent+i) != 2)
        {
          print_err("cannot read entity type", __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      return 0;
    case sel_sarray:
      if (xfscanf(in, "%k%ld", "num_sel", &n) != 2)
      {
        print_err("cannot read number of selections", __FILE__, __LINE__, __func__);
        return 2;
      }
      sarray = new sel[n];
      for (i=0; i<n; i++)
        sarray[i].read(in);
      return 0;
    default:
      break;
  }
  if ((st != sel_range_vec) && (st != sel_range_mtx) && (st != sel_vec))
  {
    if (xfscanf(in, "%k%ld", kwdstr, &n) != 2)
    {
      print_err("cannot read number of ranges or items", __FILE__, __LINE__, __func__);
      return 2;
    }
  }
  if ((st == sel_range) || (st == sel_range_mtx) || (st == sel_range_vec) || (st == sel_list) || 
      (st == sel_vec)   || (st == sel_impvallst))
  {
    id1 = new long[n];
    memset(id1, 0, sizeof(*id1)*n);
    if ((st ==  sel_range)|| (st == sel_range_mtx) || (st == sel_range_vec))
    {
      ncomp = new long[n];
      memset(ncomp, 0, sizeof(*ncomp)*n);
    }
/*    switch(st)
    {
      case sel_range:
      case sel_range_vec:
      case sel_range_mtx:
        kwdstr = "first_id";
        break;
      case sel_list:
        kwdstr = "item_id";
        break;
      default:
        break;
    }*/
    for (i=0; i<n; i++)
    {
//      if (xfscanf(in, "%k%ld", kwdstr, id1+i) != 2)
      if (xfscanf(in, "%ld", id1+i) != 1)
      {
        print_err("cannot read range or list of selection", __FILE__, __LINE__, __func__);
        return 3;
      }
      id1[i]--;
      if ((st == sel_list) || (st == sel_vec))
        continue;
//      if (xfscanf(in, "%k%ld", "ncomp", ncomp+i) != 2)
      if (xfscanf(in, "%ld", ncomp+i) != 1)
      {
        print_err("cannot read range of selection", __FILE__, __LINE__, __func__);
        return 3;
      }
    }
  }
  if ((st == sel_realrange) || (st == sel_reallist))
  {
    rid1 = new double[n];
    memset(rid1, 0, sizeof(*rid1)*n);
    if (st ==  sel_realrange)
    {
      rid2 = new double[n];
      memset(rid2, 0, sizeof(*rid2)*n);
    }
/*    switch(st)
    {
      case sel_realrange:
        kwdstr = "first_val";
        break;
      case sel_reallist:
        kwdstr = "item_val";
        break;
      default:
        break;
    }*/
    for (i=0; i<n; i++)
    {
    
//      if (xfscanf(in, "%k%lf", kwdstr, rid1+i) != 2)
      if (xfscanf(in, "%lf", rid1+i) != 1)
      {
        print_err("cannot read real range or list", __FILE__, __LINE__, __func__);
        return 4;
      }
      if (st == sel_reallist)
        continue;
//      if (xfscanf(in, "%k%lf", "sec_val", rid2+i) != 2)
      if (xfscanf(in, "%lf", rid2+i) != 1)
      {
        print_err("cannot read real range in function", __FILE__, __LINE__, __func__);
        return 4;
      }
    }
    if (st == sel_reallist)
    {
//      if (xfscanf(in, "%k%lf", "err", &err) != 2)
      if (xfscanf(in, "%lf", &err) != 1)
      {
        print_err("cannot read real range or list", __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }
  return 0;
}



/**
  Function prints data of selections to the text file given by the parameter out.
  
  @param in - pointer to opened text file

  @retval 0 - on success
  @retval 1 - error reading number of ranges
  @retval 2 - error reading range
 */
void sel::print(FILE *out)
{
  long i;
  fprintf(out, "%d", (int)st);
  switch (st)
  {
    case sel_no:
    case sel_all:
    case sel_impvalues:
    case sel_mtx:
      fprintf(out, "\n");
      break;
    case sel_vec:
      fprintf(out, " %ld\n", id1[0]+1);
      break;
    case sel_range:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%ld %ld\n", id1[i]+1, ncomp[i]);
      break;
    case sel_range_vec:
      for (i=0; i<n; i++)
        fprintf(out, " %ld %ld\n", id1[i]+1, ncomp[i]);
      break;
    case sel_range_mtx:
      for (i=0; i<n; i++)
        fprintf(out, " %ld %ld\n", id1[i]+1, ncomp[i]);
      break;
    case sel_list:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%ld\n", id1[i]+1);
      break;
    case sel_period:
      fprintf(out, " %ld\n", id1[0]);
      break;
    case sel_realperiod:
      fprintf(out, " %le %le %le %le\n", initime, fintime, r, err);
      break;
    case sel_realrange:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%e %e\n", rid1[i], rid2[i]);
      break;
    case sel_reallist:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%e\n", rid1[i]);
      fprintf(out, "%e\n", err);
      break;
    case sel_prop:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%ld %d\n", id1[i], int(ent[i]));
      break;
    case sel_impvallst:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
        fprintf(out, "%ld\n", id1[i]+1);
      break;
    case sel_sarray:
      fprintf(out, " %ld\n", n);
      for (i=0; i<n; i++)
      {
        sarray[i].print(out);
        fprintf(out, "\n");
      }
      break;
    case sel_single:
      fprintf(out, " %ld\n", id1[0]+1);
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__); 
  }
}


/**
  Function copies data of selection to the destination selection instance.
  
  @param dest - destination selection instance where the content of the actual object 
                will be copied to

  @retval 0 - on success
 */
void sel::copy(sel &dest)
{
  long i;
  dest.st = st;

  switch (st)
  {
    case sel_no:
      dest.n = 0;
      break;
    case sel_all:
    case sel_range:
    case sel_period:
    case sel_single:
    case sel_impvalues:
    case sel_mtx:
    case sel_range_vec:
    case sel_range_mtx:
      dest.n = n;
      dest.id1 = new long[n];
      dest.ncomp = new long[n];
      memcpy(dest.id1, id1, n*sizeof(*id1));
      memcpy(dest.ncomp, ncomp, n*sizeof(*ncomp));
      break;
    case sel_vec:
    case sel_list:
    case sel_impvallst:
      dest.n = n;
      dest.id1 = new long[n];
      memcpy(dest.id1, id1, n*sizeof(*id1));
      break;
    case sel_realperiod:
      dest.n = n;
      dest.initime = initime;
      dest.fintime = fintime;
      dest.r = r;
      dest.err = err;
      break;
    case sel_realrange:
      dest.n = n;
      dest.rid1 = new double[n];
      dest.rid2 = new double[n];
      memcpy(dest.rid1, rid1, n*sizeof(*rid1));
      memcpy(dest.rid2, rid2, n*sizeof(*rid2));
      break;
    case sel_reallist:
      dest.n = n;
      dest.rid1 = new double[n];
      memcpy(dest.rid1, rid1, n*sizeof(*rid1));
      dest.err = err;
      break;
    case sel_prop:
      dest.n = n;
      dest.id1 = new long[n];
      dest.ent = new gentity[n];
      memcpy(dest.id1, id1, n*sizeof(*id1));
      memcpy(dest.ent, ent, n*sizeof(*ent));
      break;
    case sel_sarray:
      dest.n = n;
      dest.sarray = new sel[n];
      for (i=0; i<n; i++)
        sarray[i].copy(dest.sarray[i]);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__); 
  }
}


/**
  Function checks presence of the value id in ranges of the selection.
  
  
  @param id - checked value

  @retval 0 - value was not found
  @retval 1 - value was found
 */
long sel::presence_id(double id)
{
  long i;
  switch (st)
  {
    case sel_realrange:
      for (i=0; i<n; i++)
        if ((id >= rid1[i]) && (id <= rid2[i]))
          return 1;
      break;
    case sel_reallist:
    case sel_realperiod:
      for (i=0; i<n; i++)
        if (fabs(id - rid1[i]) <= err)
          return 1;
      break;
    case sel_sarray:
      for (i=0; i<n; i++)
      {
        if ((sarray[i].st == sel_reallist) || 
            (sarray[i].st == sel_realperiod) || 
            (sarray[i].st == sel_realrange))
        {
          if (sarray[i].presence_id(id))
            return 1;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id in ranges of the selection.
  
  
  @param id - checked index

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_id(long id)
{
  long i;
  switch (st)
  {
    case sel_no:
      return 0;
    case sel_all:
      return 1;
    case sel_mtx:
      if (id < 6)
        return 1;
      break;
    case sel_vec:
      if (id == id1[0])
        return 1;
      break;
    case sel_range:
    case sel_range_vec:
    case sel_range_mtx:
      for (i=0; i<n; i++)
        if ((id >= id1[i]) && (id < id1[i]+ncomp[i]))
          return 1;
      break;
    case sel_list:
      for (i=0; i<n; i++)
        if (id == id1[i])
          return 1;
      break;
    case sel_period:
      if (id%id1[0] == 0)
        return 1;
      break;
    case sel_sarray:
      for (i=0; i<n; i++)
      {
        if ((sarray[i].st != sel_reallist) && 
            (sarray[i].st != sel_realperiod) && 
            (sarray[i].st != sel_realrange))
        {
          if (sarray[i].presence_id(id))
            return 1;
        }
      }
      break;
    case sel_single:
      if (id1[0] == id)
        return 1;
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id or value rid in ranges of the selection.
  
  
  @param id  - checked index
  @param rid - checked value
  

  @retval 0 - value was not found
  @retval 1 - value was found
 */
long sel::presence_id(long id, double rid)
{
  if (st == sel_sarray)
  {
    // check for both real and integer values
    return presence_id(id) + presence_id(rid);
  }
  if ((st == sel_realrange) || (st == sel_reallist) || (st == sel_realperiod))
    return presence_id(rid);
  else
    return presence_id(id);
}



/**
  Function checks presence of the index id or value rid in ranges of the selection
  or checks important times in time controler tc. The function is used
  especially for step selection.
  
  
  @param id  - checked index
  @param rid - checked value
  @param tc  - time controler from the problem description
  

  @retval 0 - value was not found
  @retval 1 - value was found
 */
long sel::presence_id(long id, double rid, timecontr &tc)
{
  long i;
  if (st == sel_sarray)
  {
    for (i=0; i<n; i++)
    {
      if(sarray[i].presence_id(id, rid, tc))
        return 1;
    }
    return 0;
  }
  if (st == sel_impvalues)
    return tc.iiit;
  if (st == sel_impvallst)
  {
    for (long i=0; i<n; i++)
    {
      if (tc.apit == id1[i])
        return tc.iiit;
    }
    return 0;
  }
  if ((st == sel_realrange) || (st == sel_reallist) || (st == sel_realperiod))
    return presence_id(rid);
  else
    return presence_id(id);
}



/**
  Function checks presence of the value id in ranges of the selection.
  
  
  @param id - checked value
  @param ir - output parameter which contains index of range where
              the id was found.

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_id(double id, long &ir)
{
  long i;
  ir = -1;
  switch (st)
  {
    case sel_realrange:
      for (i=0; i<n; i++)
      {
        if ((id >= rid1[i]) && (id <= rid2[i]))
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_reallist:
    case sel_realperiod:
      for (i=0; i<n; i++)
      {
        if (fabs(id - rid1[i]) <= err)
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_sarray:
      for (i=0; i<n; i++)
      {
        if(sarray[i].presence_id(id))
        {
          ir = i;
          return 1;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id in ranges of the selection.
  
  
  @param id - checked index
  @param ir - output parameter which contains index of range where
              the id was found.

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_id(long id, long &ir)
{
  long i;
  ir = -1;
  switch (st)
  {
    case sel_no:
      return 0;
    case sel_all:
      ir = 0;
      return 1;
    case sel_mtx:
      if (id < 6)
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_vec:
      if (id == id1[0])
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_range:
    case sel_range_vec:
    case sel_range_mtx:
      for (i=0; i<n; i++)
      {
        if ((id >= id1[i]) && (id < id1[i]+ncomp[i]))
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_list:
      for (i=0; i<n; i++)
      {
        if (id == id1[i])
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_period:
      if (id%id1[0] == 0)
      {
        ir = 0;
        return 1;
      }       
      break;
    case sel_sarray:
      for (i=0; i<n; i++)
      {
        if (sarray[i].presence_id(id))
        {
          ir = i;
          return 1;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id in ranges of the selection.
  Index of ranges have to be greater than ir. It allows for the checking 
  of all overlapping ranges of main selection and corresponding conjugated 
  selections.
  
  
  @param id - checked index
  @param ir - input/output parameter which contains 
              befor the call:
              minimum value of required index of range
              after call :
              index of range where the id was found.

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_idgt(long id, long &ir)
{
  long i;
  switch (st)
  {
    case sel_no:
      return 0;
    case sel_all:
      if (ir < 0)
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_mtx:
      if ((id < 6) && (ir < 0))
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_vec:
      if ((id == id1[0]) && (ir < 0))
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_range:
    case sel_range_vec:
    case sel_range_mtx:
      for (i=0; i<n; i++)
      {
        if ((id >= id1[i]) && (id < id1[i]+ncomp[i]) && (ir < i))
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_list:
      for (i=0; i<n; i++)
      {
        if ((id == id1[i]) && (ir < i))
        {
          ir = i;
          return 1;
        }
      }
      break;
    case sel_period:
      if ((id%id1[0] == 0) && (ir < 0))
      {
        ir = 0;
        return 1;
      }       
      break;
    case sel_sarray:
      for (i=0; i<n; i++)
      {
        if (sarray[i].presence_id(id) && (ir < i))
        {
          ir = i;
          return 1;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id in sid-th range of the selection.
  
  
  @param id - checked index
  @param ir - output parameter which contains index of range where
              the id was found.
  @param sid - index of range where the index should be found

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_id(long id, long sid, long &ir)
{
  ir = -1;
  switch (st)
  {
    case sel_no:
      return 0;
    case sel_all:
      if (sid == 0)
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_mtx:
      if ((id < 6) && (sid == 0))
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_vec:
      if ((id == id1[0]) && (sid == 0))
      {
        ir = 0;
        return 1;
      }
      break;
    case sel_range:
    case sel_range_vec:
    case sel_range_mtx:
      if ((sid < 0) && (sid >= n))
        break;
      if ((id >= id1[sid]) && (id < id1[sid]+ncomp[sid]))
      {
        ir = sid;
        return 1;
      }
      break;
    case sel_list:
      if ((sid < 0) && (sid >= n))
        break;
      if (id == id1[sid])
      {
        ir = sid;
        return 1;
      }
      break;
    case sel_period:
      if ((id%id1[0] == 0) && (sid == 0))
      {
        ir = 0;
        return 1;
      }       
      break;
    case sel_sarray:
      if ((sid < 0) && (sid >= n))
        break;
      if (sarray[sid].presence_id(id))
      {
        ir = sid;
        return 1;
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return 0;
}



/**
  Function checks presence of the index id or value rid in ranges of the selection.
  
  
  @param id  - checked index [in]
  @param rid - checked value [in]
  @param ir - output parameter which contains index of range where
              the id was found. [out]

  @retval 0 - index was not found
  @retval 1 - index was found
 */
long sel::presence_id(long id, double rid, long &ir)
{
  if (st == sel_sarray)
  {
    // check for both real and integer values
    return presence_id(id, ir) + presence_id(rid, ir);
  }
  if ((st == sel_realrange) || (st == sel_reallist) || (st == sel_realperiod))
    return presence_id(rid, ir);
  else
    return presence_id(id, ir);
}



/**
  Function checks presence of the index id in ranges of the selection and
  also checks presence of the index idc in conjugated selections consel.
  Conjugated selections means array of selections whose size is number of 
  ranges of given selection. Each element (i.e. selection) of this array 
  corresponds to one range of given selection.
  
  @param consel - array of conjugated selections [in]
  @param id  - checked index of given selection  [in]
  @param idc - checked index of of conjugated selections [in]

  @retval 0 - indices were not found
  @retval 1 - indices were found
 */
long sel::presence_id(sel *consel, long id, long idc)
{
  long ir=-1;
  long i;

  for (i=0; i < n; i++)
  {
    if (presence_idgt(id, ir))
    {
      if (consel[ir].presence_id(idc))
        return 1;
    }
    else
      return 0;
  }
  
  return 0;
}



/**
  Function checks presence of the index id in ranges of the selection and
  also checks presence of the index idc in conjugated selections consel.
  Conjugated selections means array of selections whose size is number of 
  ranges of given selection. Each element (i.e. selection) of this array 
  corresponds to one range of given selection.
  
  @param consel - array of conjugated selections [in]
  @param id  - checked index of given selection  [in]
  @param idc - checked index of of conjugated selections [in]
  @param ir  - index of conjugated selection where the idc was found [out]

  @retval 0 - indices were not found
  @retval 1 - indices were found
 */
long sel::presence_id(sel *consel, long id, long idc, long &ir)
{
  ir=-1;
  long i;

  for (i=0; i < n; i++)
  {
    if (presence_idgt(id, ir))
    {
      if (consel[ir].presence_id(idc))
        return 1;
    }
    else
      return 0;
  }
  
  return 0;
}



/** 
  Function computes number of selected items.
  @param tncomp - total number of components which could be selected

  @retval number of selected components

  TKo 9.2008
*/
long sel::give_nselcomp(long tncomp)
{
  long ret=0;
  long i, di;

  switch(st)
  {
    case sel_no:
      ret=0;
      break;
    case sel_all:
      ret=tncomp;
      break;
    case sel_range:
      for(i=0; i<n; i++)
      {
        if (id1[i] >= tncomp) // id is completely out of maximum value
          continue;
        di = id1[i]+ncomp[i]-1; // id of last component in the given range
        if (di < tncomp)
	{
          ret += ncomp[i];  // all components of the given range can be selected
          continue;         
	}
        else
          ret += di+1-tncomp; // only part of range can be selected
      }
      break;
    case sel_list:
      for(i=0; i<n; i++)
      {
        if (id1[i] < tncomp)    ret++;
      }  
      break;
    case sel_single:
      ret = 1;
      break;
    case sel_period:
      ret = 0;
      for (i=0; i<tncomp; i+=id1[0])
        ret++;
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
  }
  return ret;
}



/**
  The function converts selection of nodes/elements according to given property
  id and entity type to selection types list or range.

  @param top - mesh topology with defined properties of nodes and elements
  @param ot  - object type indicator for the given property selection
               it specifies whether the given selection by property id is related  
               to nodes or elements
  @param conselent - array of conjugated selections with particular entities
  @param consel - resulting array of conjugated selections with created ranges or list items (output)
  @param transent - array of conjugated transformation indicators with particular entities
  @param trans - array of conjugated transformation indicators with created ranges or list items (output)

  @retval 0 - on success
  @retval 1 - in case of unknown type of object
  @retval 2 - no edge or surface property ids were specified on elements in the given topology

  Created by Tomas Koudelka, 9.1.2015
*/
long sel::conv_selprop(siftop *top, objtype ot, sel *conselent, sel *&consel, long *transent, long *&trans)
{
  long i,j;
  long nlit, nrit, nold;
  long *nlitent, *nritent;
  long ret; 
  long **selobjent;
  long nobj;
  long *entid;
  long nentid;

  switch (ot)
  {
    case gelem:
      nobj = top->ne;
      break;
    case gnod:
      nobj = top->nn;
      break;
    default:
      print_err("unknown type of object(%d) is required", __FILE__, __LINE__, __func__, ot);
      return 1;
  }
  nlitent = new long[n];
  memset(nlitent, 0, sizeof(*nlitent)*n);
  nritent = new long[n];
  memset(nritent, 0, sizeof(*nritent)*n);
  nold = n;
  
  selobjent = new long*[n];
  memset(selobjent, 0, sizeof(*selobjent)*n);

  for (i=0; i<n; i++)  // loop for particular specified entities
  {
    selobjent[i] = new long[nobj];
    memset(selobjent[i], 0, sizeof(*selobjent[i])*nobj);
    switch (ot)
    {
      case gelem:
        for (j=0; j<top->ne; j++)
        {
          entid = NULL;
          nentid = 0;
          ret = top->elements[j].searchprop(id1[i], ent[i], NULL, j, top->edges, top->surfaces, entid, nentid);
          if (ret < 0)
          {
            print_err("no edge or surface properties are given on elements", __FILE__, __LINE__, __func__);
            return 2;
          }
          selobjent[i][j] += ret;
        }
        break;
      case gnod:
        entid = NULL;
        nentid = 0;
        top->get_propent_nodes(id1[i], ent[i], selobjent[i]);
        for(j=0; j<nobj; ++j)
          ++selobjent[i][j];
        break;
      default:
        print_err("unknown type of object(%d) is required", __FILE__, __LINE__, __func__, ot);
        return 1;
    }
    nlitent[i] = give_num_lst_items(nobj, selobjent[i]);   // number of  selected objects of the given entity
    nritent[i] = give_num_range_items(nobj, selobjent[i]); // number of ranges in the selection of objects of the given entity
  }

  // compute total number of records for the resulting list and range selection types
  nlit = nrit = 0;
  for (j=0; j<n; j++)
  {
    nlit += nlitent[j];
    nrit += nritent[j];
  }

  delete [] id1;
  delete [] ncomp;
  delete [] ent;
  ent = NULL;

  if (2*nrit < nlit)  // range needs two numbers for the specification
  {
    // conversion to the range type selection is more favourable    
    consel = new sel[nrit];
    if (transent)
    {
      trans = new long[nrit];
      memset(trans, 0, sizeof(*trans)*nrit);
    }
    else
      trans = NULL;
    conv2range(nrit, nobj, selobjent, conselent, consel, transent, trans); 
  }
  else 
  {
    // conversion to the list type selection is more favourable
    consel = new sel[nlit];
    if (transent)
    {
      trans = new long[nlit];
      memset(trans, 0, sizeof(*trans)*nlit);
    }
    else
      trans = NULL;
    conv2lst(nlit, nobj, selobjent, conselent, consel, transent, trans); 
  }

  for (i=0; i < nold; i++) // n was changed in course of conversion to new range or list selection => nold must be used
    delete [] selobjent[i];
  delete [] selobjent;
  delete [] nlitent;
  delete [] nritent;

  return 0;
}



/**
  The function converts selection of nodes/elements according to given property
  id and entity type to selection types list or range.

  @param[in] top - mesh topology with defined properties of nodes and elements
  @param[in] ot  - object type indicator for the given property selection
                   it specifies whether the given selection by property id is related  
                   to nodes or elements

  @retval 0 - on success
  @retval 1 - in case of unknown type of object
  @retval 2 - no edge or surface property ids were specified on elements in the given topology

  Created by Tomas Koudelka, 09.2023
*/
long sel::conv_selprop(const siftop *top, objtype ot)
{
  long i,j;
  long nlit, nrit, nold;
  long *nlitent, *nritent;
  long ret; 
  long **selobjent;
  long nobj;
  long *entid;
  long nentid;

  switch (ot)
  {
    case gelem:
      nobj = top->ne;
      break;
    case gnod:
      nobj = top->nn;
      break;
    default:
      print_err("unknown type of object(%d) is required", __FILE__, __LINE__, __func__, ot);
      return 1;
  }
  nlitent = new long[n];
  memset(nlitent, 0, sizeof(*nlitent)*n);
  nritent = new long[n];
  memset(nritent, 0, sizeof(*nritent)*n);
  nold = n;
  
  selobjent = new long*[n];
  memset(selobjent, 0, sizeof(*selobjent)*n);

  for (i=0; i<n; i++)  // loop for particular specified entities
  {
    selobjent[i] = new long[nobj];
    memset(selobjent[i], 0, sizeof(*selobjent[i])*nobj);
    switch (ot)
    {
      case gelem:
        for (j=0; j<top->ne; j++){
          entid = NULL;
          nentid = 0;
          ret = top->elements[j].searchprop(id1[i], ent[i], NULL, j, top->edges, top->surfaces, entid, nentid);
          if (ret < 0){
            print_err("no edge or surface properties are given on elements", __FILE__, __LINE__, __func__);
            return 2;
          }
          selobjent[i][j] += ret;
        }
        break;
      case gnod:
        entid = NULL;
        nentid = 0;
        top->get_propent_nodes(id1[i], ent[i], selobjent[i]);
        for(j=0; j<nobj; ++j)
          ++selobjent[i][j];
        break;
      default:
        print_err("unknown type of object(%d) is required", __FILE__, __LINE__, __func__, ot);
        return 1;
    }
    nlitent[i] = give_num_lst_items(nobj, selobjent[i]);   // number of  selected objects of the given entity
    nritent[i] = give_num_range_items(nobj, selobjent[i]); // number of ranges in the selection of objects of the given entity
  }

  // compute total number of records for the resulting list and range selection types
  nlit = nrit = 0;
  for (j=0; j<n; j++){
    nlit += nlitent[j];
    nrit += nritent[j];
  }

  delete [] id1;
  delete [] ncomp;
  delete [] ent;
  ent = NULL;

  if (2*nrit < nlit){
    // range needs two numbers for the specification
    // conversion to the range type selection is more favourable    
    conv2range(nrit, nobj, selobjent); 
  }
  else{ 
    // conversion to the list type selection is more favourable
    conv2lst(nlit, nobj, selobjent); 
  }

  for (i=0; i < nold; i++) // n was changed in course of conversion to new range or list selection => nold must be used
    delete [] selobjent[i];
  delete [] selobjent;
  delete [] nlitent;
  delete [] nritent;

  return 0;
}



/**
  The function returns number of selected objects given in the array selobj.
  
  @param nobj - length of the array selobj, i.e. the total number of selectable objects
  @param selobj - array of indicators whether the given object is selected (selobj[i] != 0) 
                  or not (selobj[i] = 0)

  @return The function returns number of selected objects.

  Created by Tomas Koudelka, 9.1.2015
*/
long sel::give_num_lst_items(long nobj, long *selobj)
{
  long i, j;

  j = 0;
  for (i=0; i<nobj; i++)
  {
    if (selobj[i])   j++;
  }

  return j;
}



/**
  The function returns number of contiuous index ranges of selected objects given by selobj 
  array.
  
  @param nobj - length of the array selobj, i.e. the total number of selectable objects
  @param selobj - array of indicators whether the given object is selected (selobj[i] != 0) 
                  or not (selobj[i] = 0)

  @return The number of detected index ranges in the selobj array.

  Created by Tomas Koudelka, 9.1.2015
*/
long sel::give_num_range_items(long nobj, long *selobj)
{
  long i, j, start_r;
 
  j = 0;
  start_r = 0;
  for(i=0; i<nobj; i++)
  {
    if (selobj[i]) // i-th object is selected
    {
      if (start_r == 0) // no range is started => start new range
        start_r = 1;  // set indicator that we have started new range

      continue;
    }
    else // i-th object is not selected 
    {
      if (start_r) // a range has been started => finish the range
      {
        j++;  // increase range index
        start_r = 0; // switch off the new range indicator
      }
    }
  }
  if (start_r) // a range has been started => finish the range
  {
    j++;  // increase range index
    start_r = 0; // switch off the new range indicator
  }

  return j;
}



/**
  The function converts list of selected objects to the range type of selection.
  It is supposed that the number of selected ranges was detected by call of give_num_range_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining  components are zero.
  The function is used for the selection of nodes or elements according to property id in the 
  outdriver class.
  
  @param nit - number of detected ranges
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - 2D array of indicators for particular entities whether the given j-th object 
                  is selected for i-th entity (selobj[i][j] != 0) or not (selobj[i][j] = 0)
  @param conselent - array of conjugated selections with particular entities
  @param conselr - resulting array of conjugated selections with created ranges
  @param transent - array of conjugated transformation indicators with particular entities
  @param transr - array of conjugated transformation indicators with created ranges

  @return The function does not return anything but it changes internal parameters of the class 
          to the range selection type.

  Created by Tomas Koudelka, 08.2016
*/
void sel::conv2range(long nit, long nobj, long **selobj, sel *conselent, sel *conselr,
                     long *transent, long *transr)
{
  long i, j, k, start_r, nold = n;
 
  st = sel_range;
  n = nit;
  id1 = new long[n];
  ncomp = new long[n];
  memset(id1, 0, sizeof(*id1)*n);
  memset(ncomp, 0, sizeof(*ncomp)*n);
  

  j = 0;
  for (i=0; i < nold; i++)
  {
    start_r = 0;
    for(k=0; k<nobj; k++)
    {
      if (selobj[i][k]) // k-th object is selected
      {
        if (start_r == 0) // no range is started => start new range
        {
          id1[j] = k;
          start_r = 1;  // set indicator that we have started new range
        }
        continue;
      }
      else // k-th object is not selected 
      {
        if (start_r) // a range has been started => finish the range
        {
          ncomp[j] = k-id1[j]; // number of components in the range
          conselent[i].copy(conselr[j]);
          if (transent)
            transr[j] = transent[i];
          j++;  // increase range index
          start_r = 0; // switch off the new range indicator
        }
      }
    }
    if (start_r) // a range has been started => finish the range
    {
      ncomp[j] = k-id1[j]; // number of components in the range
      conselent[i].copy(conselr[j]);
      if (transent)
        transr[j] = transent[i];
    }
  }
}



/**
  The function converts list of selected objects to the list type of selection.
  It is supposed that the number of selected items was detected by call of give_num_lst_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining components are zero.
  The function is used for the selection of nodes or elements according to property id in the 
  outdriver class.
  
  @param nit - number of list items
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - 2D array of indicators for particular entities whether the given j-th object 
                  is selected for i-th entity (selobj[i][j] != 0) or not (selobj[i][j] = 0)
  @param conselent - array of conjugated selections with particular entities
  @param consell - resulting array of conjugated selections with created list items
  @param transent - array of conjugated transformation indicators with particular entities
  @param transr - array of conjugated transformation indicators with created ranges

  @return The function does not return anything but it changes internal parameters of the class 
          to the list selection type.

  Created by Tomas Koudelka, 9.1.2015
*/
void sel::conv2lst(long nit, long nobj, long **selobj, sel *conselent, sel *consell, 
                   long *transent, long *transl)
{
  long i, j, k, nold = n;

  st = sel_list;
  n = nit;
  id1 = new long[n];
  ncomp = NULL;
  memset(id1, 0, sizeof(*id1)*n);

  j = 0;
  for (k=0; k<nold; k++)
  {
    for (i=0; i<nobj; i++)
    {
      if (selobj[k][i])
      {
        id1[j]=i;
        conselent[k].copy(consell[j]);
        if (transent)
          transl[j] = transent[k];
        j++;
      }
    }
  }
}



/**
  The function converts list of selected objects to the range type of selection.
  It is supposed that the number of selected ranges was detected by call of give_num_range_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining  components are zero.
  The function is used for the selection of nodes or elements according to property id in the 
  outdriver class.
  
  @param nit - number of detected ranges
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - 2D array of indicators for particular entities whether the given j-th object 
                  is selected for i-th entity (selobj[i][j] != 0) or not (selobj[i][j] = 0)

  @return The function does not return anything but it changes internal parameters of the class 
          to the range selection type.

  Created by Tomas Koudelka, 09.2023
*/
void sel::conv2range(long nit, long nobj, long **selobj)
{
  long i, j, k, start_r, nold = n;
 
  st = sel_range;
  n = nit;
  id1 = new long[n];
  ncomp = new long[n];
  memset(id1, 0, sizeof(*id1)*n);
  memset(ncomp, 0, sizeof(*ncomp)*n);
  

  j = 0;
  for (i=0; i < nold; i++){
    start_r = 0;
    for(k=0; k<nobj; k++){
      if (selobj[i][k]){
        // k-th object is selected
        if (start_r == 0){
          // no range is started => start new range
          id1[j] = k;
          start_r = 1;  // set indicator that we have started new range
        }
        continue;
      }
      else{
        // k-th object is not selected 
        if (start_r){
          // a range has been started => finish the range
          ncomp[j] = k-id1[j]; // number of components in the range
          j++;  // increase range index
          start_r = 0; // switch off the new range indicator
        }
      }
    }
    if (start_r) // a range has been started => finish the range
      ncomp[j] = k-id1[j]; // number of components in the range
  }
}



/**
  The function converts list of selected objects to the list type of selection.
  It is supposed that the number of selected items was detected by call of give_num_lst_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining components are zero.
  The function is used for the selection of nodes or elements according to property id in the 
  outdriver class.
  
  @param nit - number of list items
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - 2D array of indicators for particular entities whether the given j-th object 
                  is selected for i-th entity (selobj[i][j] != 0) or not (selobj[i][j] = 0)

  @return The function does not return anything but it changes internal parameters of the class 
          to the list selection type.

  Created by Tomas Koudelka, 9.1.2015
*/
void sel::conv2lst(long nit, long nobj, long **selobj)
{
  long i, j, k, nold = n;

  st = sel_list;
  n = nit;
  id1 = new long[n];
  ncomp = NULL;
  memset(id1, 0, sizeof(*id1)*n);

  j = 0;
  for (k=0; k<nold; k++){
    for (i=0; i<nobj; i++){
      if (selobj[k][i]){
        id1[j]=i;
        j++;
      }
    }
  }
}



/**
  The function converts list of selected objects to the range type of selection.
  It is supposed that the number of selected ranges was detected by call of give_num_range_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining  components are zero.
  The function is used for the selection of nodes with prescribed initial displacements by the given rotation
  in gradual construction problems.
  
  @param nit - number of detected ranges
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - integer %vector of indicators whether the given i-th object is selected (selobj[i] != 0) 
                  or not (selobj[i] = 0)

  @return The function does not return anything but it changes internal parameters of the class 
          to the range selection type.

  Created by Tomas Koudelka, 08.2016
*/
void sel::conv2range(long nit, long nobj, ivector &selobj)
{
  long j, k, start_r;
 
  st = sel_range;
  n = nit;
  id1 = new long[n];
  ncomp = new long[n];
  memset(id1, 0, sizeof(*id1)*n);
  memset(ncomp, 0, sizeof(*ncomp)*n);
  

  j = 0;
  start_r = 0;
  for(k=0; k<nobj; k++)
  {
    if (selobj[k]) // k-th object is selected
    {
      if (start_r == 0) // no range is started => start new range
      {
        id1[j] = k;
        start_r = 1;  // set indicator that we have started new range
      }
      continue;
    }
    else // k-th object is not selected 
    {
      if (start_r) // a range has been started => finish the range
      {
        ncomp[j] = k-id1[j]; // number of components in the range
        j++;  // increase range index
        start_r = 0; // switch off the new range indicator
      }
    }
  }
  if (start_r) // a range has been started => finish the range
    ncomp[j] = k-id1[j]; // number of components in the range
}



/**
  The function converts list of selected objects to the list type of selection.
  It is supposed that the number of selected items was detected by call of give_num_lst_items.
  The list of selected objects is represented by the array selobj with length nobj. Selected
  object is indicated by nonzero selobj component while the remaining components are zero.
  The function is used for the selection of nodes with prescribed initial displacements by the given rotation
  in gradual construction problems.
  
  @param nit - number of list items
  @param nobj - length of array selobj, i.e. the total number of selectable objects
  @param selobj - integer %vector of indicators whether the given i-th object is selected (selobj[i] != 0) 
                  or not (selobj[i] = 0)

  @return The function does not return anything but it changes internal parameters of the class 
          to the list selection type.

  Created by Tomas Koudelka, 08.2016
*/
void sel::conv2lst(long nit, long nobj, ivector &selobj)
{
  long i, j;

  st = sel_list;
  n = nit;
  id1 = new long[n];
  ncomp = NULL;
  memset(id1, 0, sizeof(*id1)*n);

  j = 0;
  for (i=0; i<nobj; i++)
  {
    if (selobj[i])
    {
      id1[j]=i;
      j++;
    }
  }
}



