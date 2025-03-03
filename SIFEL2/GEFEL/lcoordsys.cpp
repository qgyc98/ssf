#include "lcoordsys.h"
#include "gfmatrix.h"
#include "matrix.h"
#include "vector.h"
#include "gnode.h"
#include "iotools.h"

#include <stdlib.h>


const double lcoordsys::zero = 1.0e-10;
const char *lcoordsys::namepar[] = {"x", "y", "z", "t"};

lcoordsys::lcoordsys() : dim{3L}, tt{no_tr}, caxis{NULL}, paxis{NULL}, tmat{NULL}, gftmat{NULL}
{  
}



lcoordsys::lcoordsys(long pdim) : dim{pdim}, tt{no_tr}, caxis{NULL}, paxis{NULL}, tmat{NULL}, gftmat{NULL}
{  
}



lcoordsys::~lcoordsys()
{
  delete paxis;
  delete caxis;  
  delete tmat;
}



/**
  The function reads definition of one local coordinate system from the opened text file.

  @param[in] in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error of the reading of coordinate system definition

  Created by Tomas Koudelka, 09.2023
*/
void lcoordsys::read (XFILE *in)
{  
  xfscanf(in, "%k%m", "transf_type", &tmat_type_kwdset, &tt);
  switch (tt){
    case gen_tr:
      tmat = new matrix(dim, dim);
      readm(in, *tmat);
      break;
    case cyl_tr:
      if (dim != 3)
        dim = 3;
      paxis = new vector(dim);
      caxis = new vector(dim);     
      readv(in, *paxis);
      readv(in, *caxis);
      normalize(*caxis);
      break;
    case gfgen_tr:
      gftmat = new gfmatrix(dim, dim);
      readm(in, *tmat);
      break;
    default:
      print_err("unknown type %d of a local coordinate system\n",
                __FILE__, __LINE__, __func__, tt);
      abort();
  }
}



/**
  The function reads definition of one local coordinate system from the opened text file.

  @param[in] in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error of the reading of coordinate system definition

  Created by Tomas Koudelka, 09.2023
*/
void lcoordsys::print (FILE *out)
{  
  fprintf(out, "%d\n", int(tt));
  switch (tt){
    case gen_tr:
      printm(out, *tmat);
      break;
    case cyl_tr:
      printv(out, *paxis);
      printv(out, *caxis);
      break;
    case gfgen_tr:
      printm(out, *tmat);
      break;
    default:
      print_err("unknown type %d of a local coordinate system\n",
                __FILE__, __LINE__, __func__, tt);
      abort();
  }
  fprintf(out, "\n");
}



/**
  The function returns transformation matrix in the argument t.
  
  @param[in,out] t - transformation matrix, it must be allocated to dimensions (dim,dim).
  @param[in]     pc - point coordinates at which the transformation matrix will be assembled.
  @param[in]     time - actual time for the transformation with the help of a general function matrix

  Created by Tomas Koudelka, 09.2023
*/
void lcoordsys::give_transfmat(matrix &t, vector &pc, double time)
{
  vector auxv;
  switch (tt){
    case no_tr:
      identm(t);
      break;
    case gen_tr:
      if (copym(*tmat, t)){
        print_err("cannot assemble of the lcs transformation matrix\n",
                  __FILE__, __LINE__, __func__);
        abort();
      }
      break;
    case cyl_tr:
      give_cyl_transfmat(t, pc);
      break;
    case gfgen_tr:
      reallocv(RSTCKVEC(4, auxv));
      auxv(0) = pc(0);
      auxv(1) = pc(1);
      auxv(2) = pc(2);
      auxv(3) = time;
      gftmat->evaluate(auxv, namepar, t);      
      break;
    default:
      print_err("unknown type %d of a local coordinate system\n",
                __FILE__, __LINE__, __func__, tt);
      abort();
  }
}



/**
  The function assembles transformation matrix for a cylindrical domain. The matrix is returned in 
  the argument t.
  
  @param[in,out] t - transformation matrix, it must be allocated to dimensions (dim,dim)
  @param[in]     pc - point coordinates at which the transformation matrix will be assembled.

  Created by Tomas Koudelka, 09.2023
*/
void lcoordsys::give_cyl_transfmat(matrix &t, vector &pc)
{
  vector lx(ASTCKVEC(3));
  vector ly(ASTCKVEC(3));
  vector lz(ASTCKVEC(3));

  copyv(*caxis, lz);
  lx(0L) = pc(0) - (*paxis)(0L);
  lx(1L) = pc(1) - (*paxis)(1L);
  lx(2L) = pc(2) - (*paxis)(2L);
  double nlx = normv(lx);
  if (nlx < zero){
    lx(0L) = 1.0;
    lx(1L) = lx(2L) = 0.0;
  }
  double aux;
  scprd(lx, lz, aux);
  lx(0L) -= aux*lz(0L);
  lx(1L) -= aux*lz(1L);
  lx(2L) -= aux*lz(2L);
  normalize(lx);

  crprd(lx, lz, ly);
  normalize(ly);

  // assemble transformation matrix from the base vectors
  for (long i=0L; i<3L; i++){
    t(i, 0L) = lx(i);
    t(i, 1L) = ly(i);
    t(i, 2L) = lz(i);
  }
}
