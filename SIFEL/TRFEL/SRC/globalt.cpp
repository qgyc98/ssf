#include "globalt.h"
#include "elemheadt.h"

/**
  The function initializes global variables to null values.

  Created by Tomas Koudelka, 08.2021
*/
void initnull_globt (void)
{
  Ndoft = Mesprt = 0;
  Outt = Outt1 = Outt2 = NULL;

  Tp = NULL;
  Gtt = NULL;
  Tt = NULL;
  Tm = NULL;
  Tc = NULL;
  Tb = NULL;

  Lbt  = NULL;
  Lbt3d = NULL;
  Lbat = NULL;
  Qbt   = NULL;
  Qbat  = NULL;
  Ltt   = NULL;
  Ltat  = NULL;
  Lqt   = NULL;
  Qqt   = NULL;
  Qqat  = NULL;
  Lqat  = NULL;
  Ltett = NULL;
  Lht   = NULL;
  Qht   = NULL;
  Lwt   = NULL;
  G2d   = NULL;
  Outdt = NULL;

  Lsrst = NULL;
  Stt  = NULL;
  Adat = NULL;

  Kmat = NULL;
  Cmat = NULL;
  Jmat = NULL;
  Bmat = NULL;
}



/**
  The function deallocates global variables and sets them to null values.
  Deallocation should be in reverse order of allocation in trfelinit procedure,
  otherwise there could be memory errors in cases that some global variables depend 
  on other.

  Created by Ladislav Svoboda?
  Modified by Tomas Koudelka, 08.2021
*/
void delete_globt (void)
{
  //fprintf(Outt, "\ndelete_globt():\n Tm->ip = %p", Tm->ip);
  delete Kmat;     Kmat = NULL;
  delete Cmat;	   Cmat = NULL;
  delete Jmat;	   Jmat = NULL;
  delete Bmat;	   Bmat = NULL;

  delete Lbt;	   Lbt = NULL;
  delete Lbt3d;    Lbt3d = NULL;
  delete Lbat;	   Lbat = NULL;
  delete Qbt;	   Qbt = NULL;
  delete Qbat;	   Qbat = NULL;
  delete Ltt;	   Ltt = NULL;
  delete Ltat;	   Ltat = NULL;
  delete Lqt;	   Lqt = NULL;
  delete Qqt;	   Qqt = NULL;
  delete Qqat;	   Qqat = NULL;
  delete Lqat;	   Lqat = NULL;
  delete Ltett;	   Ltett = NULL;
  delete Lht;	   Lht = NULL;
  delete Qht;	   Qht = NULL;
  delete Lwt;      Lwt = NULL;
  delete G2d;	   G2d = NULL;

  delete Stt;      Stt = NULL;
  delete Lsrst;	   Lsrst = NULL;
  delete Adat;	   Adat = NULL;

  delete Outdt;	   Outdt = NULL;
  delete Tb;	   Tb = NULL;
  delete Tc;	   Tc = NULL;
  delete Tm;	   Tm = NULL;
  delete Tt;	   Tt = NULL;
  delete Gtt;	   Gtt = NULL;
  delete Tp;       Tp = NULL;
  
  Ndoft = Mesprt = 0;
}

