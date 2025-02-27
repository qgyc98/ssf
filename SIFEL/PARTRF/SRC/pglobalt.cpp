#include "pglobalt.h"



/**
  The function initializes global variables to null values.

  Created by Tomas Koudelka, 08.2021
*/
void pinitnull_globt(void)
{
  Nproc  = 0;
  Myrank = 0;
  Ndom   = 0;
  proc_name[0] = '\0';
  nameLength = 0;

  Ptp = NULL;
  Psolt = NULL;
}



/**
  The function deallocates global variables and sets them to null values.
  Deallocation should be in reverse order of allocation in ptrfelinit procedure,
  otherwise there could be memory errors in cases that some global variables depend 
  on other.

  Modified by Tomas Koudelka, 08.2021
*/
void pdelete_globt(void)
{
  delete Ptp;    Ptp = NULL;
  delete Psolt;  Psolt = NULL;

  Nproc  = 0;
  Myrank = 0;
  Ndom   = 0;
  proc_name[0] = '\0';
  nameLength = 0;
}

