#include <stdlib.h>
#include <string.h>
#include "hdbcontr.h"
#include "galias.h"



/**
  The constructor initializes data to the default values - 
  no backup is set by default

  Created by Tomas Koudelka, 06.2008  
*/
hdbcontr::hdbcontr()
{
  hdbtype = nohdb;  // no backup
  hdbfmtr = hdbfmts = text;
  rmold = no;
  rmold_id = -1;
  hdbnamer[0] = 0;
  hdbnames[0] = 0;
  selother_r = selother_s = NULL;
  selother_id = NULL;
}



/**
  The destructor releases allocated memory.

  Created by Tomas Koudelka, 06.2008  
*/
hdbcontr::~hdbcontr()
{
  delete [] selother_r;
  delete [] selother_s;
  for (long i=0; i<selelemr.n; i++)
    delete [] selother_id[i];
  delete [] selother_id;
  
}



/**
  The function reads the setup of the backup controller form the opened text file.

  @param in - pointer to the opened text file

  @return The function does not return anything.
 
  Created by Tomas Koudelka, 06.2008
*/
void hdbcontr::read(XFILE *in)
{
  long i,j;
  // indicator of backup on harddisk
  xfscanf (in, "%k%m", "hdbackup", &hdbackuptype_kwdset, (int *)&hdbtype);
  if (hdbtype)
  {
    switch(hdbtype)
    {
      case hdbr_single:
      case hdbr_multiple:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmtr", &hdbackupfmttype_kwdset, (int *)&hdbfmtr);
        // name of backup file
        xfscanf(in, " %a", hdbnamer);
        // selection of other components for backup restoring
        xfscanf(in, "%k", "elems_saved");
        selelemr.read(in);
        switch (selelemr.st)
        {
          case sel_no:
            break;
          case sel_all:
          case sel_range:
          case sel_list:
            xfscanf(in, "%k", "other_comp");
            selother_r = new sel[selelemr.n];
            selother_id = new long*[selelemr.n];
            memset(selother_id, 0, sizeof(*selother_id)*selelemr.n);
            for (i=0; i<selelemr.n; i++)
	    {
              selother_r[i].read(in);
              selother_id[i] = new long[selother_r[i].n];
              memset(selother_id[i], 0, sizeof(*selother_id[i])*selother_r[i].n);
              for (j=0; j<selother_r[i].n; j++)
	      {
                xfscanf(in, "%k%ld", "other_id", selother_id[i]+j);
                selother_id[i][j]--;
	      }
	    }
            break;
          default:
	    print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
            abort();
        }
        break;
      case hdbs_single:
      case hdbs_multiple:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmts", &hdbackupfmttype_kwdset, (int *)&hdbfmts);
        xfscanf(in, "%k%m","hdbackup_rm_old", &answertype_kwdset, (int *)&rmold);
        if (hdbfmts==text)
          xfscanf(in,"%k%ld","precision",&prec);
        // name of backup file
        xfscanf(in, " %a", hdbnames);
        // selection of other components for backup saving
        xfscanf(in, "%k", "elems_other");
        selelems.read(in);
        switch (selelems.st)
        {
          case sel_no:
            break;
          case sel_all:
          case sel_range:
          case sel_list:
            xfscanf(in, "%k", "other_comp");
            selother_s = new sel[selelems.n];
            for (i=0; i<selelems.n; i++)
              selother_s[i].read(in);
            break;
          default:
	    print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
            abort();
        }
        break;
      case hdbrs_single:
      case hdbrs_multiple:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmtr", &hdbackupfmttype_kwdset, (int *)&hdbfmtr);
        // name of backup file
        xfscanf(in, " %a", hdbnamer);
        // selection of other components for backup restoring
        xfscanf(in, "%k", "elems_saved");
        selelemr.read(in);
        switch (selelemr.st)
        {
          case sel_no:
            break;
          case sel_all:
          case sel_range:
          case sel_list:
            xfscanf(in, "%k", "other_comp");
            selother_r = new sel[selelemr.n];
            selother_id = new long*[selelemr.n];
            memset(selother_id, 0, sizeof(*selother_id)*selelemr.n);
            for (i=0; i<selelemr.n; i++)
	    {
              selother_r[i].read(in);
              selother_id[i] = new long[selother_r[i].n];
              memset(selother_id[i], 0, sizeof(*selother_id[i])*selother_r[i].n);
              for (j=0; j<selother_r[i].n; j++)
	      {
                xfscanf(in, "%k%ld", "other_id", selother_id[i]+j);
                selother_id[i][j]--;
	      }
	    }
            break;
          default:
	    print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
            abort();
        }
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmts", &hdbackupfmttype_kwdset, (int *)&hdbfmts);
        xfscanf(in, "%k%m","hdbackup_rm_old", &answertype_kwdset, (int *)&rmold);
        if (hdbfmts==text)
          xfscanf(in,"%k%ld","precision",&prec);
        // name of backup file
        xfscanf(in, " %a", hdbnames);
        // selection of other components for backup saving
        xfscanf(in, "%k", "elems_other");
        selelems.read(in);
        switch (selelems.st)
        {
          case sel_no:
            break;
          case sel_all:
          case sel_range:
          case sel_list:
            xfscanf(in, "%k", "other_comp");
            selother_s = new sel[selelems.n];
            for (i=0; i<selelems.n; i++)
              selother_s[i].read(in);
            break;
          default:
	    print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
            abort();
        }
        break;
      case hdbr_nonloc:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmtr", &hdbackupfmttype_kwdset, (int *)&hdbfmtr);
        // name of backup file
        xfscanf(in, " %a", hdbnamer);
        break;
      case hdbs_nonloc:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmts", &hdbackupfmttype_kwdset, (int *)&hdbfmts);
        // name of backup file
        xfscanf(in, " %a", hdbnames);
        break;
      case hdbrs_nonloc:
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmtr", &hdbackupfmttype_kwdset, (int *)&hdbfmtr);
        // name of backup file
        xfscanf(in, " %a", hdbnamer);
        // format of backup file
        xfscanf(in, "%k%m", "hdbackup_fmts", &hdbackupfmttype_kwdset, (int *)&hdbfmts);
        // name of backup file
        xfscanf(in, " %a", hdbnames);
        break;
      default:
        print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
    }
  }
}



/**
  The function prints the setup of the backup controller to the opened text file.
  The main usage is in the preprocessor.

  @param out - pointer to the opened text file

  @return The function does not return anything.
 
  Created by Tomas Koudelka, 06.2008
*/
void hdbcontr::print(FILE *out)
{
  long i, j;
  // indicator of backup on harddisk
  fprintf(out, "%d", hdbtype);
  switch(hdbtype)
  {
    case nohdb:
      fprintf(out, "\n");
      break;
    case hdbr_single:
    case hdbr_multiple:
      // format of backup file
      fprintf(out, " %d\n", hdbfmtr);
      // name of backup file
      fprintf(out,"%s\n", hdbnamer);
      // selection of other components for backup restoring
      selelemr.print(out);
      switch (selelemr.st)
      {
        case sel_no:
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          for (i=0; i<selelemr.n; i++)
          {
            selother_r[i].print(out);
            for(j=0; j<selother_r[i].n; j++)
              fprintf(out, " %ld\n", selother_id[i][j]+1);
          }
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          abort();
      }
      break;
    case hdbs_single:
    case hdbs_multiple:
      // format of backup file
      fprintf(out, " %d", hdbfmts);
      fprintf(out, " %d", (int)rmold);
      if (hdbfmts==text)
        fprintf(out, " %ld", prec);
      fprintf(out, "\n");
      // name of backup file
      fprintf(out, "%s\n", hdbnames);
      // selection of other components for backup saving
      selelems.print(out);
      switch (selelems.st)
      {
        case sel_no:
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          for (i=0; i<selelems.n; i++)
            selother_s[i].print(out);
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          abort();
      }
      break;
    case hdbrs_single:
    case hdbrs_multiple:
      // format of backup file
      fprintf(out, " %d\n", hdbfmtr);
      // name of backup file
      fprintf(out, "%s\n", hdbnamer);
      // selection of other components for backup restoring
      selelemr.print(out);
      switch (selelemr.st)
      {
        case sel_no:
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          for (i=0; i<selelemr.n; i++)
          {
            selother_r[i].print(out);
            for(j=0; j<selother_r[i].n; j++)
              fprintf(out, "%ld\n", selother_id[i][j]+1);
          }
          break;
          default:
	    print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
            abort();
      }
      // format of backup file
      fprintf(out, " %d", hdbfmts);
      fprintf(out, " %d", (int)rmold);
      if (hdbfmts==text)
        fprintf(out, " %ld", prec);
      fprintf(out, "\n");
      // name of backup file
      fprintf(out, "%s\n", hdbnames);
      // selection of other components for backup saving
      selelems.print(out);
      switch (selelems.st)
      {
        case sel_no:
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          for (i=0; i<selelems.n; i++)
            selother_s[i].print(out);
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          abort();
      }
      break;
    case hdbr_nonloc:
      // format of backup file
      fprintf(out, " %d\n", hdbfmtr);
      // name of backup file
      fprintf(out, "%s\n", hdbnamer);
      break;
    case hdbs_nonloc:
      // format of backup file
      fprintf(out, "%d", hdbfmts);
      // name of backup file
      fprintf(out, "%s\n", hdbnames);
      break;
    case hdbrs_nonloc:
      // format of backup file
      fprintf(out, " %d\n", hdbfmtr);
      // name of backup file
      fprintf(out, "%s\n", hdbnamer);
      // format of backup file
      fprintf(out, " %d\n", hdbfmts);
      // name of backup file
      fprintf(out, "%s\n", hdbnames);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function returns status of restorage

  Returns:
   @retval 0 = no restorage is required
           1 = restorage is required

  Created by Tomas Koudelka, 06.2008
*/
long hdbcontr::restore_stat(void)
{
  long ret;

  switch(hdbtype)
  {
    case hdbr_single:
    case hdbrs_single:
    case hdbr_multiple:
    case hdbrs_multiple:
      ret=1;
      break;
    default:
      ret=0;
  }
  return ret;
}   



/**
  Function returns status of saving

  Returns:
   @retval 0 = no backup required
           1 = backup required

  Created by Tomas Koudelka, 06.2008
*/
long hdbcontr::save_stat(void)
{
  long ret;

  switch(hdbtype)
  {
    case hdbs_single:
    case hdbrs_single:
    case hdbs_multiple:
    case hdbrs_multiple:
      ret=1;
      break;
    default:
      ret=0;
  }
  return ret;
}   
