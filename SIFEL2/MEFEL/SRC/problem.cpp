#include <stdio.h>

#include "problem.h"
#include "global.h"
#include "probdesc.h"
#include "gtopology.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "gmatrix.h"



problem::problem (void)
{
  mp = NULL;    gt = NULL;
  mt = NULL;    mm = NULL;
  mc = NULL;    mb = NULL;
  lsrs = NULL;
  smat = NULL;
}

problem::~problem (void)
{

}

void problem::globinic (void)
{
  mp   = Mp;
  gt   = Gtm;
  mt   = Mt;
  mm   = Mm;
  mc   = Mc;
  mb   = Mb;
  lsrs = Lsrs;
  smat = Smat;
}

void problem::deinic (void)
{
  mp = NULL;    gt = NULL;
  mt = NULL;    mm = NULL;
  mc = NULL;    mb = NULL;
  lsrs = NULL;
  smat = NULL;
}

void problem::dealoc (void)
{
  delete mp;
  delete gt;
  delete mt;
  delete mm;
  delete mc;
  delete mb;
  delete lsrs;
  delete smat;
}
