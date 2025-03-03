#include "globalc.h"
#include "mechmat.h"
#include "mechtop.h"

/**
  The function initializes global variables to null values.

  Created by Tomas Koudelka, 08.2021
*/
void initnull_globc(void)
{
  Ndofc=0;
  Mesprc=0;
  Outc = NULL;

  Cp = NULL;
  Gtu = NULL;
  Ct = NULL;
  Cm = NULL;
  Cmu = NULL;
  Cml = NULL;
  Cb = NULL;

  TM_ip = NULL;
  MT_ip = NULL;
  TMipmap = NULL;
  MTipmap = NULL;
  TM_nod_map = NULL;

  Lsrsc = NULL;
  Outdc = NULL;

  D0mat = NULL;
  D1mat = NULL;

  Cbar = NULL;
  Cquad = NULL;
  Caxiq = NULL;
  Caxifc = NULL;
  Chex = NULL;
}



/**
  The function deallocates global variables and sets them to null values.
  Deallocation should be in reverse order of allocation in metrinit procedure,
  otherwise there could be memory errors in cases that some global variables depend 
  on other.

  Modified by Tomas Koudelka, 08.2021
*/
void delete_globc(void)
{
  delete D0mat;  D0mat = NULL;
  delete D1mat;  D1mat = NULL;

  delete Cbar;    Cbar = NULL;
  delete Cquad;   Cquad = NULL;
  delete Caxiq;   Caxiq = NULL;
  delete Caxifc;  Caxifc = NULL;
  delete Chex;    Chex = NULL;

  delete Lsrsc;  Lsrsc = NULL;
  delete Outdc;  Outdc = NULL;

  delete [] TM_ip;       TM_ip = NULL;
  delete [] MT_ip;       MT_ip = NULL;
  delete [] TMipmap;     TMipmap = NULL;
  delete [] MTipmap;     MTipmap = NULL;
  delete [] TM_nod_map;  TM_nod_map = NULL;

  delete Cm;   Cm = NULL;
  delete Cmu;  Cmu = NULL;
  delete Cml;  Cml = NULL;
  delete Cb;   Cb = NULL;
  delete Ct;   Ct = NULL;
  delete Gtu;  Gtu = NULL;
  delete Cp;   Cp = NULL;

  Ndofc=0;
  Mesprc=0;
}
