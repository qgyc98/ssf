#include "endnodem.h"
#include "global.h"
#include "gtopology.h"
#include "globmat.h"
#include <string.h>

/**
   constructor
   
   JK, 1.3.2009
*/
endnodem::endnodem (void)
{
  //  number of approximated functions
  napfun = 0;
  if (napfun<1){
    print_err("wrong number of approximated functions is required",__FILE__,__LINE__,__func__);
  }
  
  //  number of assigned general endnode from GEFEL
  nen = 0;
  
  //  nodal displacements in the global coordinate system
  u1=0.0;  u2=0.0;
  v1=0.0;  v2=0.0;

  //  jumps in the displacement field
  ju=0.0;  jv=0.0;

  //  type of material model
  tm = NULL;
  //  id of material model
  idm = NULL;

  //  auxiliary array for nodal displacements
  //  this class serves only for 2D problems
  //  it means that nodes contain 2 DOFs
  r = new double [2];
}

endnodem::~endnodem (void)
{
  delete [] tm;
  delete [] idm;
  delete [] r;
}


/**
   function assembles nodal displacements
   
   @param lcid - load case id
   
   JK, 1.3.2009
*/
void endnodem::nodal_displacements (long lcid)
{
  long nid;
  
  //  first node at the end node
  nid=Gtm->endnodes[nen].fn;
  //  nodal displacements
  noddispl (lcid,r,nid);
  u1=r[0];  v1=r[1];
  
  //  second node on the end node
  nid=Gtm->gedges[nen].ln;
  //  nodal displacements
  noddispl (lcid,r,nid);
  u2=r[0];  v2=r[1];
  
} 

/**
   computes jumps in displacement field at the end node
   
   JK, 1.3.2009
*/
void endnodem::compute_jumps ()
{
  ju = u2-u1;
  jv = v2-v1;
}

/**
   function computes derivative of the end node functional with respect
   to the nodal displacements
   the end node functional is called lipschitz continuous perturbation or
   local lipschitz functional
   
   @param lcid - load case id
   @param v - array containing the %vector of derivative
   
   JK, 1.3.2009
*/
void endnodem::compute_endnode_functional_derivative (long lcid,double */*v*/)
{
  //  nodal displacements at nodes on edge
  nodal_displacements (lcid);
  //  magnitudes of jumps between two subdomains
  compute_jumps ();
  
}

