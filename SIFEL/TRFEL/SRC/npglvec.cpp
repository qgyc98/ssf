#include "npglvec.h"
#include "vector.h"
#include "globalt.h"

/**
  Default constructor initializes data members
  to the default values.

  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
*/
np_glob_vec::np_glob_vec()
{
  istep = ncd = nce = 0;
  lhs = tdlhs = rhs = f = d = p = fb = fi = NULL;
  v = z = lhsb = tdlhsb = NULL;
  ifn = NULL;
}



/**
  The function allocates vectors tdlhs, f, d, p
  with length n. Vectors lhs and rhs remains unallocated because
  they are used as pointers to the right and left hand sides 
  that are allocated separately.

  @param n - length of vectors allocated i.e. total number of unknowns

  @return The function does not return anything.

  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
*/
void np_glob_vec::alloc(long n)
{
  //  vector of prescribed fluxes (right hand side)
  f   = new double [n];
  nullv (f,n);
  //  predictor
  d   = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p   = new double [n];
  nullv (p,n);

  //  backup of nodal values
  lhsb = new double [n];
  nullv (lhsb,n);
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];
  nullv (tdlhsb,n);
}

/**
  The function allocates auxiliary vectors fb, fi
  
  @param n - length of vectors allocated i.e. total number of unknowns
  
  @return The function does not return anything.
  
  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
*/
void np_glob_vec::alloc_aux(long n)
{
  //  auxiliary vector
  fb = new double [n];
  nullv (fb,n);
  fi = new double [n];
  nullv (fi,n);
}



/**
  The function allocates auxiliary vectors v, z, lhsb and tdlhsb for 
  the dform solver type.
  
  @param n - length of vectors allocated i.e. total number of unknowns
  
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 6/2014
*/
void np_glob_vec::alloc_daux(long n)
{
  //  auxiliary vector
  v = new double [n];
  nullv (v,n);
  //  auxiliary vector
  z = new double [n];
  nullv (z,n);
}



/**
  The function deallocates vectors lhs, tdlhs, rhs, f, d, p.

  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
*/
void np_glob_vec::dealloc()
{
  delete [] f; 
  delete [] d;  
  delete [] p;
  
  if(fb != NULL)
    delete [] fb;
  if(fi != NULL)
    delete [] fi;

  if(v != NULL)
    delete [] v;
  if(z != NULL)
    delete [] z;
  if(lhsb != NULL)
    delete [] lhsb;
  if(tdlhsb != NULL)
    delete [] tdlhsb;

  //delete [] flp;
  if(ifn != NULL)
    delete [] ifn;
}



/**
  Destructor deallocates vectors lhs, tdlhs, rhs, f, d, p.

  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
*/
np_glob_vec::~np_glob_vec()
{
  dealloc();
}
