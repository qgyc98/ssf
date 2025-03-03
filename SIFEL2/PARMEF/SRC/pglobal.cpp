#include "pglobal.h"



/**
  The function initializes global variables to null values.

  Created by Tomas Koudelka, 08.2021
*/
void pinitnull_glob(void)
{
  Nproc  = 0;
  Myrank = 0;
  Ndom   = 0;
  proc_name[0] = '\0';
  nameLength = 0;

  Pmp = NULL;
  Psolm = NULL;
}



/**
  The function deallocates global variables and sets them to null values.
  Deallocation should be in reverse order of allocation in mefelinit procedure,
  otherwise there could be memory errors in cases that some global variables depend 
  on other.

  Modified by Tomas Koudelka, 08.2021
*/
void pdelete_glob(void)
{
  delete Pmp;    Pmp = NULL;
  delete Psolm;  Psolm = NULL;

  Nproc  = 0;
  Myrank = 0;
  Ndom   = 0;
  proc_name[0] = '\0';
  nameLength = 0;
}

