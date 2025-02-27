#include "global.h"
#include "gmatrix.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "mechcrsec.h"
#include "lhsrhs.h"
#include "outdriverm.h"
#include "adaptivity.h"
#include "stochdriver.h"
#include "flsubdom.h"
#include "elemhead.h"

/**
  The function initializes global variables to null values.

  Created by Tomas Koudelka, 08.2021
*/
void initnull_glob(void)
{
  Ndofm = 0;
  Mespr = 0;
  Out = NULL;

  Mp = NULL;
  Gtm = NULL;
  Mt = NULL;
  Mtt = NULL;
  Mm = NULL;
  Mmm = NULL;
  Mc = NULL;
  Mb = NULL;
  Outdm = NULL;
  
  Lsrs = NULL;
  Ada = NULL;
  St = NULL;
  Fsd = NULL;
  
  Smat = NULL;
  Mmat = NULL;
  Dmat = NULL;
  Ismat = NULL;
  Amat = NULL;
  
  Bar2d = NULL;
  Bar3d = NULL;
  Barq2d = NULL;
  Barq3d = NULL;
  Beam2d = NULL;
  Beam3d = NULL;
  Beam3dg = NULL;
  Sbeam = NULL;
  Spbeam2d = NULL;
  Spring = NULL;
  Pelt = NULL;
  Peqt = NULL;
  Perlt = NULL;
  Pelq = NULL;
  Peqq = NULL;
  Perlq = NULL;
  Pesqt = NULL;
  Cct = NULL;
  Dkt = NULL;
  Dst = NULL;
  Q4pl = NULL;
  Argtr = NULL;
  Argtrpl = NULL;
  Qkirch = NULL;
  Dkqelem = NULL;
  Spltr = NULL;
  Splq = NULL;
  Asymlq = NULL;
  Asymlt = NULL;
  Asymqq = NULL;
  Asymqt = NULL;
  Asymcq = NULL;
  Asymlqifc = NULL;
  Shtr = NULL;
  Shtrm = NULL;
  Shq = NULL;
  Ltet = NULL;
  Qtet = NULL;
  Lhex = NULL;
  Qhex = NULL;
  Ltetrot = NULL;
  Lhexrot = NULL;
  Lwed = NULL;
  Qwed = NULL;
  Pelem = NULL;
  Pqifc = NULL;
  Hexifc = NULL;
  Tlatt = NULL;
  
  Neval = 0.0;
  Omp_wtime = 0.0;
  Elemfrc = NULL;
//#ifndef INC_OPENMP
//  out = NULL;
//  out2 = NULL;
//  in = NULL;
//#endif
//#ifdef INC_OPENMP
//  out = NULL;
//  out2 = NULL;
//#endif
}



/**
  The function deallocates global variables and sets them to null values.
  Deallocation should be in reverse order of allocation in mefelinit procedure,
  otherwise there could be memory errors in cases that some global variables depend 
  on other.

  Modified by Tomas Koudelka, 08.2021
*/
void delete_glob(void)
{
  delete Smat;   Smat = NULL;
  delete Mmat;   Mmat = NULL;
  delete Dmat;   Dmat = NULL;
  delete Ismat;  Ismat = NULL;
  delete Amat;   Amat = NULL;

  delete Ada;    Ada = NULL;
  delete Fsd;    Fsd = NULL;
  delete St;     St = NULL;
  delete Lsrs;   Lsrs = NULL;

  delete Outdm;  Outdm = NULL;
  delete Mb;     Mb = NULL;
  delete Mc;     Mc = NULL;
  delete Mmm;    Mmm = NULL;
  delete Mm;     Mm = NULL;
  delete Mtt;    Mtt = NULL;
  delete Mt;     Mt = NULL;
  delete Gtm;    Gtm = NULL;
  delete Mp;     Mp = NULL;
  
  delete Bar2d;     Bar2d = NULL;
  delete Bar3d;     Bar3d = NULL;
  delete Barq2d;    Barq2d = NULL;
  delete Barq3d;    Barq3d = NULL;
  delete Beam2d;    Beam2d = NULL;
  delete Beam3d;    Beam3d = NULL;
  delete Beam3dg;   Beam3dg = NULL;
  delete Sbeam;     Sbeam = NULL;
  //delete Spbeam2d;  Spbeam2d = NULL;
  delete Spring;    Spring = NULL;
  delete Pelt;      Pelt = NULL;
  delete Peqt;      Peqt = NULL;
  delete Perlt;     Perlt = NULL;
  delete Pelq;      Pelq = NULL;
  delete Peqq;      Peqq = NULL;
  delete Perlq;     Perlq = NULL;
  delete Pesqt;     Pesqt = NULL;
  delete Cct;       Cct = NULL;
  delete Dkt;       Dkt = NULL;
  delete Dst;       Dst = NULL;
  delete Q4pl;      Q4pl = NULL;
  delete Argtr;     Argtr = NULL;
  delete Argtrpl;   Argtrpl = NULL;
  delete Qkirch;    Qkirch = NULL;
  delete Dkqelem;   Dkqelem = NULL;
  delete Spltr;     Spltr = NULL;
  delete Splq;      Splq = NULL;
  delete Asymlq;    Asymlq = NULL;
  delete Asymlt;    Asymlt = NULL;
  delete Asymqq;    Asymqq = NULL;
  delete Asymqt;    Asymqt = NULL;
  delete Asymcq;    Asymcq = NULL;
  delete Asymlqifc; Asymlqifc = NULL;
  delete Shtr;      Shtr = NULL;
  delete Shtrm;     Shtrm = NULL;
  delete Shq;       Shq = NULL;
  delete Ltet;      Ltet = NULL;
  delete Qtet;      Qtet = NULL;
  delete Lhex;      Lhex = NULL;
  delete Qhex;      Qhex = NULL;
  delete Ltetrot;   Ltetrot = NULL;
  delete Lhexrot;   Lhexrot = NULL;
  delete Lwed;      Lwed = NULL;
  delete Qwed;      Qwed = NULL;
  delete Pelem;     Pelem = NULL;
  delete Pqifc;     Pqifc = NULL;
  delete Hexifc;    Hexifc = NULL;
  delete Tlatt;     Tlatt = NULL;

  delete [] Elemfrc;

// #ifndef INC_OPENMP
//   fclose(in);
//   fclose(out);
//   fclose(out2);
// #endif
// #ifdef INC_OPENMP
//   fclose(out);
//   fclose(out2);
// #endif
  Ndofm = Mespr = 0;
}

