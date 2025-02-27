#include "mtglvec.h"
#include "probdesc.h"
#include "vector.h"
#include "global.h"



/**
  Default constructor initializes data members
  to the default values.

  Created by Tomas Koudelka, 11.2011
*/
mt_glob_vec::mt_glob_vec()
{
  istep = ncd = nce = 0;
  r = f = dr = fp = fl = fi = fb = lhsb = flp = NULL;
  ifn = NULL;
}



/**
  The function allocates vectors dr, fp, fl, fi, fb and lhsb
  with length n. Vectors r and f remains unallocated because
  they are used as pointers to the right and left hand sides 
  that are allocated separately.

  @param n - length of vectors allocated i.e. total number of unknowns

  @return The function does not return anything.

  Created by Tomas Koudelka, 11.2011
*/
void mt_glob_vec::alloc(long n)
{
  //  vector of increments of nodal displacements
  dr   = new double [n];
  nullv (dr,n);
  //  vector of prescribed forces from the previous step
  fp   = new double [n];
  nullv (fp,n);
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc
  fl   = new double [n];
  nullv (fl,n);
  //  vector of internal forces
  fi   = new double [n];
  nullv (fi,n);
  //  auxiliary force vector
  fb   = new double [n];
  nullv (fb,n);
  //  backup of the nodal displacements
  lhsb = new double [n];
  nullv (lhsb,n);
  
  if (Mp->tprob == growing_mech_structure)
  {
    //  vector of prescribed forces from the previous step due to prescribed force load
    flp   = new double [n];
    nullv (flp,n);
    // vector of indicators for nodes on the interface between new and old parts of the structure
    ifn   = new long [n];
    nullv (ifn,n);
  }
}



/**
  The function deallocates vectors dr, fp, flp, fl, fi, fb, lhsb and ifn.
*/
void mt_glob_vec::dealloc()
{
  delete [] fi;  
  delete [] fb;  
  delete [] fp;  
  delete [] dr; 
  delete [] fl;  
  delete [] lhsb;
  
  delete [] flp;
  delete [] ifn;
}



/**
  Destructor deallocates vectors dr, fp, flp, fl, fi, fb, lhsb and ifn.

  Created by Tomas Koudelka, 11.2011
*/
mt_glob_vec::~mt_glob_vec()
{
  dealloc();
}
