#include "edgem.h"
#include "global.h"
#include "gtopology.h"
#include "globmat.h"

#include <math.h>

/**
   constructor
   
   JK, 28.8.2007
*/
edgem::edgem (void)
{
  //long fn,ln;
  
  /*
  //  number of nodes on edge
  nn = i;
  if (nn<2){
    fprintf (stderr,"\n\n wrong number of nodes on edge (number of nodes %ld)",nn);
    fprintf (stderr,"\n in constructor edgem (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  //  number of approximated functions
  napfun = j;
  if (napfun<1){
    fprintf (stderr,"\n\n wrong number of approximated functions on edge (number of functions %ld)",napfun);
    fprintf (stderr,"\n in constructor edgem (file %s, line %d)\n",__FILE__,__LINE__);
  }
  */

  //  number of assigned general edge from GEFEL
  ned = 0;
  
  //  nodal displacements in the global coordinate system
  u1=0.0;  u2=0.0;  u3=0.0;  u4=0.0;
  v1=0.0;  v2=0.0;  v3=0.0;  v4=0.0;

  //  tangentional displacements (displacements in the direction defined by the edge)
  td1=0.0;  td2=0.0;  td3=0.0;  td4=0.0;
  
  //  normal displacements (displacements normal to the direction defined by the edge)
  nd1=0.0;  nd2=0.0;  nd3=0.0;  nd4=0.0;
  
  //  jumps in the tangential and normal directions
  jt1=0.0;  jt2=0.0;  jn1=0.0;  jn2=0.0;
  
  //  auxiliary array for nodal displacements
  //  this class serves only for 2D problems
  //  it means that nodes contain 2 DOFs
  r = new double [2];

  //  type of material model
  tm = NULL;
  //  id of material model
  idm = NULL;
}

edgem::~edgem ()
{
  delete [] tm;
  delete [] idm;
  delete [] r;
}

/**
   JK, 8.8.2007
*/
void edgem::read (FILE */*in*/)
{

}

/**
   function initializes edge
   
   @param edid - edge id
   
   JK, 1.3.2009
*/
void edgem::init (long /*edid*/)
{

}


/**
   function assembles nodal displacements
   
   example:
   ---3------4---
   ---7------8---
   nlist[0]=3, nlist[1]=7, nlist[2]=4, nlist[3]=8
   
   @param lcid - load case id
   
   JK, 28.8.2007
*/
void edgem::nodal_displacements (long lcid)
{
  long nid;
  
  //  first node on the edge
  nid=Gtm->gedges[ned].nlist[0];
  //  nodal displacements
  noddispl (lcid,r,nid);
  u1=r[0];  v1=r[1];
  
  //  second node on the edge
  nid=Gtm->gedges[ned].nlist[1];
  //  nodal displacements
  noddispl (lcid,r,nid);
  u2=r[0];  v2=r[1];
  
  //  third node on the edge
  nid=Gtm->gedges[ned].nlist[2];
  //  nodal displacements
  noddispl (lcid,r,nid);
  u3=r[0];  v3=r[1];
  
  //  fourth node on the edge
  nid=Gtm->gedges[ned].nlist[3];
  //  nodal displacements
  noddispl (lcid,r,nid);
  u4=r[0];  v4=r[1];
  
} 

/**
   function computes tangential and normal displacements
   
   JK, 28.8.2007
*/
void edgem::tan_nor_displacements ()
{
  //  tangent direction
  Gtm->gedges[ned].give_dirvect (r);
  
  td1=u1*r[0]+v1*r[1];
  td2=u2*r[0]+v2*r[1];
  td3=u3*r[0]+v3*r[1];
  td4=u4*r[0]+v4*r[1];
  
  //  normal direction
  Gtm->gedges[ned].give_norvect (r);
  
  nd1=u1*r[0]+v1*r[1];
  nd2=u2*r[0]+v2*r[1];
  nd3=u3*r[0]+v3*r[1];
  nd4=u4*r[0]+v4*r[1];
  
}

/**
   computes jumps in displacement field along the edge
   
   JK, 28.8.2007
*/
void edgem::compute_jumps ()
{
  //  jump in tangential direction
  jt1 = td2-td1;
  jt2 = td4-td3;

  //  jump in normal direction
  jn1 = nd2-nd1;
  jn2 = nd4-nd3;
}

/**
   function computes derivative of the edge functional with respect
   to the nodal displacements
   the edge functional is called lipschitz continuous perturbation or
   local lipschitz functional
   
   @param lcid - load case id
   @param v - array containing the %vector of derivative
   
   JK, 1.3.2009
*/
void edgem::compute_edge_functional_derivative (long lcid,double */*v*/)
{
  //  nodal displacements at nodes on edge
  nodal_displacements (lcid);
  //  displacements in the tangential and normal directions
  tan_nor_displacements ();
  //  magnitudes of jumps between two subdomains
  compute_jumps ();
  
}

