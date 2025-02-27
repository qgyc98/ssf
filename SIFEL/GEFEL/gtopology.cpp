#include "gtopology.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <utility>
#include <algorithm>
#include "galias.h"
#include "ordering.h"
#include "gmatrix.h"
#include "basefun.h"
#include "dofrange.h"
#include "gfunct.h"
#include "iotools.h"
#include "boundbox.h"

gtopology::gtopology ()
{
  //  number of nodes
  nn=0;
  //  number of elements
  ne=0;
  //  number of layered nodes
  nln=0;
  //  number of degrees of freedom in the problem
  ndof=0;
  //  number of internal DOFs
  nidof=0;
  //  number of boundary/interface DOFs
  nbdof=0;
  //  number of unknowns in saddle point problem
  nsad=0;
  //  number of nodes which are switched on
  //  it is used in problems with changing number of nodes and elements
  nnso=0;
  //  number of elements which are switched on
  //  it is used in problems with changing number of nodes and elements
  neso=0;
  //  number of subdomains/aggregates
  ns=0;
  //  number of phases
  nph=0;
  //  indicator of sequential topology reading
  //  sequential topology is not read
  rst=0;
  //  hanging nodes (default - no hanging nodes)
  hangnod=0;
  
  
  //  general nodes
  gnodes=NULL;
  //  general elements
  gelements=NULL;
  //  layered general nodes
  lgnodes=NULL;
  
  domsizes = NULL;
  
  //  array of numbers of adjacent elements to nodes
  nadjelnod=NULL;
  //  array of adjacent elements to nodes
  adjelnod=NULL;
  //  array of numbers of adjacent elements to elements
  nadjelel=NULL;
  //  array of adjacent elements to elements
  adjelel=NULL;
  //  numbers of adjacent nodes to nodes
  nadjnodnod=NULL;
  //  array of adjacent nodes to nodes
  adjnodnod=NULL;


  bckcn=NULL;  unln=NULL;  unnl=NULL;
  cngtopcorr=NULL;

  gphases=NULL;
  
  
  
  //  DOF control
  //  default value - DOFs are independent on time, they do not change
  dofcontr=0;
  //  type of generator of code numbers
  //  default value - the classical generation of code numbers
  cngen=1;
  //  reordered node numbers
  ordering=NULL;
  //  type of renumbering
  nodren = no_renumbering;
  

  //
  //  TOOLS FOR PROBLEMS WITH CHANGING NUMBER OF ELEMENTS
  //
  //  state of code numbers
  cnstate=0;
  //  general functions
  gf=NULL;
  //  list of nodes which are switched on
  lnso=NULL;
  //  list of elements which are switched on
  leso=NULL;
  
  
  //  type of egdes
  edtype = noedge;
  //  number of general edges
  nged=0;
  //  number of edge series
  nser=0;
  //  series of edges
  edgeser=NULL;

  nedser=NULL;

  edgelist = NULL;
  
  //  general edges
  gedges=NULL;
  
  //  number of end nodes
  nen=0;
  //  end nodes
  endnodes=NULL;

  //  sequential topology
  stop=NULL;
  
  nadjnodnod = NULL;
  
  adjnodnod = NULL;  
  
  
  //  arrays for hanging node determination
  //  arrays of coordinates of element centers
  xc = NULL;
  yc = NULL;
  zc = NULL;
  //  array of radii of circumscribed balls
  circumrad = NULL;

  eldofr = NULL;
  ompelord = NULL;
}

gtopology::~gtopology ()
{
  long i;
  delete [] gnodes;  
  delete [] gelements;
  delete [] lgnodes;
  
  delete [] domsizes;
  
  delete [] nadjelnod;  
  if (adjelnod)
  {
    for (i=0;i<nn;i++)
      delete [] adjelnod[i];
  }
  delete [] adjelnod;
  
  delete [] nadjelel;
  if (adjelel)
  {
    for (i=0;i<ne;i++)
      delete [] adjelel[i];
  }
  delete [] adjelel;
  
  delete [] bckcn;
  delete [] unln;  
  delete [] unnl;
  delete [] cngtopcorr;

  delete [] gphases;
  delete [] ordering;

  //
  //  TOOLS FOR PROBLEMS WITH CHANGING NUMBER OF ELEMENTS
  //
  delete [] gf;
  delete [] lnso;
  delete [] leso;
  
  delete [] nedser;
  delete [] edgeser;
  delete [] gedges;

  if (edgelist)
  {
    for (i=0;i<nser;i++)
      delete [] edgelist[i];
    delete [] edgelist;
  }

  delete [] endnodes;
  
  //  sequential toplogy
  delete stop;

  delete [] nadjnodnod;

  if (adjnodnod)
  {
    for (i=0;i<nn;i++)
      delete [] adjnodnod[i];
    delete [] adjnodnod;
  }

  delete [] eldofr;
  delete [] ompelord;
}

/**
   function initiates class gtopology
   
   @param icn - array containing code numbers
   @param nc - array containing numbers of components
   @param m - number of generalized elements
   
   24.7.2002
*/
void gtopology::initiate (long **icn,long *nc,long m)
{
  long i;
  ne=m;
  cnstate=1;
  alloc_elements (ne);
  for (i=0;i<ne;i++){
    gelements[i].initiate (icn[i],nc[i]);
  }
}

/**
   function allocates nodes
   
   @param m - the number of nodes
*/
void gtopology::alloc_nodes (long m)
{
  nn=m;
  if (m<0){
    print_err("negative number of nodes is required", __FILE__, __LINE__, __func__);
  }
  gnodes = new gnode [nn];
}

/**
   function allocates elements
   
   @param m - the number of elements
*/
void gtopology::alloc_elements (long m)
{
  ne=m;
  if (m<0){
    print_err("negative number of elements is required", __FILE__, __LINE__, __func__);
  }
  gelements = new gelement [ne];
}

/**
   function allocates layered nodes
   
   @param m - the number of layered nodes
*/
void gtopology::alloc_lnodes (long m)
{
  nln=m;
  if (m<0){
    print_err("negative number of nodes is required", __FILE__, __LINE__, __func__);
  }
  lgnodes = new lgnode [nln];
}



/**
   function allocates phases
   
   @param m - number of nodes
   6.4.2005 TKr
*/
void gtopology::alloc_phases (long m)
{
  nph=m;
  if (m<0){
    print_err("negative number of phases is required", __FILE__, __LINE__, __func__);
  }
  gphases = new gphase [nph];
}

/**
   function returns number of degrees of freedom of required node
   
   @param nid - node id (number of node)
*/
long gtopology::give_ndofn (long nid)
{
  long ndofn;
  ndofn = gnodes[nid].give_ndofn ();
  return ndofn;
}

/**
   function returns number of multipliers defined between two layers in one node
   
   @param lnid - layerd node id (number of layered node)
*/
long gtopology::give_nmult (long lnid)
{
  return lgnodes[lnid].nmult;
}

/**
   function returns number of degrees of freedom of required element
   
   @param eid - element id (number of element)
 */
long gtopology::give_ndofe (long eid)
{
  return gelements[eid].give_ndofe ();
}

/**
   function returns number of all degrees of freedom of required element
   it means with Lagrange multipliers defined on element
   
   @param eid - element id (number of element)
*/
long gtopology::give_gndofe (long eid)
{
  long ndofe,nne;
  nne = gelements[eid].give_nne ();
  ndofe = gelements[eid].give_ndofe ();
  ndofe+= nne*2*gelements[eid].give_nmult ();
  return ndofe;
}

/**
   function returns number of nodes on selected edge
   
   @param edid - edge id
*/
long gtopology::give_nnedge (long edid)
{
  return gedges[edid].nn;
}

/**
   function returns the number of nodes on required element
   
   @param eid - element id (number of element)
*/
long gtopology::give_original_nne (long eid)
{
  return gelements[eid].give_nne ();
}

/**
   function returns the number of nodes on required element
   
   @param eid - element id (number of element)
*/
long gtopology::give_nne (long eid)
{
  long nne;
  
  if (gelements[eid].give_nmne () == 0){
    //  there are no hanging nodes
    nne = gelements[eid].give_nne ();
  }else{
    //  there are hanging nodes
    nne = gelements[eid].give_nmne ();
  }
  return nne;
}

/**
   function returns the number of master nodes on required element
   
   @param eid - element id (number of element)
*/
long gtopology::give_nmne (long eid)
{
  return gelements[eid].give_nmne ();
}

/**
   function returns indicator of code numbers defined on element
   
   @param eid - element id (number of element)
*/
long gtopology::give_cne (long eid)
{
  return gelements[eid].give_cne ();
}

/**
   function returns required degree of freedom of required node
   
   @param nid - node id
   @param n - number of degree of freedom
*/
long gtopology::give_dof (long nid,long n)
{
  return gnodes[nid].give_dof (n);
}

/**
   function sets one code number on node
   
   @param nid - node id
   @param n - number of degree of freedom
   @param num - value of code number
*/
void gtopology::save_dof (long nid,long n,long num)
{
  gnodes[nid].save_dof (n,num);
}

/**
   function returns numbers of nodes on element
   
   @param eid - element id
   @param nod - array containing node numbers
*/
void gtopology::give_nodes (long eid,ivector &nod)
{
  if (gelements[eid].nmne==0){
    //  there are no hanging nodes
    gelements[eid].give_nodes (nod);
  }else{
    //  there are hanging nodes
    gelements[eid].give_master_nodes (nod);
  }
}

/**
   function returns numbers of nodes on element
   
   @param eid - element id
   @param nod - array containing node numbers
*/
void gtopology::give_original_nodes (long eid,ivector &nod)
{
  gelements[eid].give_nodes (nod);
}

/**
   function returns numbers of nodes on element
   
   @param eid - element id
   @param nod - array containing node numbers
   
   JK, 17.6.2012
*/
void gtopology::give_master_nodes (long eid,ivector &nod)
{
  gelements[eid].give_master_nodes (nod);
}


/**
   function returns coordinates of the appropriate node
   
   @param nid - node id (input)
   @param coord - vector containing nodal coordinates (output)
   
   Created by Tomas Koudelka, 07.2016
*/
void gtopology::give_node_coord(long nid, vector &coord)
{
  if (coord.n!=3)
    print_err("wrong size of array (%ld<3) is required", __FILE__, __LINE__, __func__, coord.n);

  if ((nid < 0) && (nid >= nn))
    print_err("invalid node id nid=%ld is required, nid must be from <0,%ld>", __FILE__, __LINE__, __func__, nid, nn-1); 

  coord[0] = gnodes[nid].x;  
  coord[1] = gnodes[nid].y;  
  coord[2] = gnodes[nid].z;  
}


/**
   function returns node coordinates of the appropriate element
   
   @param x - %vector containing node coordinates
   @param eid - element id
   
   16.7.2002
*/
void gtopology::give_node_coord1d (vector &x,long eid)
{
  long i,nne;
  nne = gelements[eid].give_nne ();
  if (x.n!=nne){
    print_err("wrong size of array is required", __FILE__, __LINE__, __func__);
  }
  ivector nodes(ASTCKIVEC(nne));
  gelements[eid].give_nodes (nodes);
  
  for (i=0;i<nne;i++){
    x[i]=gnodes[nodes[i]].x;
  }
}


/**
   function returns node coordinates of the appropriate element
   
   @param x,y - vectors containing node coordinates
   @param eid - element id
   
   16.7.2002
*/
void gtopology::give_node_coord2d (vector &x,vector &y,long eid)
{
  long i,nne;
  nne = gelements[eid].give_nne ();
  if (x.n!=nne || y.n!=nne){
    print_err("wrong size of array is required", __FILE__, __LINE__, __func__);
  }
  ivector nodes(ASTCKIVEC(nne));
  gelements[eid].give_nodes (nodes);
  
  for (i=0;i<nne;i++){
    x[i]=gnodes[nodes[i]].x;
    y[i]=gnodes[nodes[i]].y;
  }
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,z - vectors containing coordinates
   @param eid - element id
   
   16.7.2002
*/
void gtopology::give_node_coord2dxz (vector &x,vector &z,long eid)
{
  long i,nne;
  nne = gelements[eid].give_nne ();
  if (x.n!=nne || z.n!=nne){
    print_err("wrong size of array is required", __FILE__, __LINE__, __func__);
  }
  ivector nodes(ASTCKIVEC(nne));
  gelements[eid].give_nodes (nodes);
  
  for (i=0;i<nne;i++){
    x[i]=gnodes[nodes[i]].x;
    z[i]=gnodes[nodes[i]].z;
  }
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,y,z - vectors containing coordinates
   @param eid - element id
   
   16.7.2002
*/
void gtopology::give_node_coord3d (vector &x,vector &y,vector &z,long eid)
{
  long i,nne;
  nne = gelements[eid].give_nne ();
  if (x.n!=nne || y.n!=nne || z.n!=nne){
    print_err("wrong size of array is assumed.",__FILE__,__LINE__,__func__);
  }
  ivector nodes(ASTCKIVEC(nne));
  gelements[eid].give_nodes (nodes);
  
  for (i=0;i<nne;i++){
    x[i]=gnodes[nodes[i]].x;
    y[i]=gnodes[nodes[i]].y;
    z[i]=gnodes[nodes[i]].z;
  }
}

/**
   function extracts all code numbers of actual element
   
   @param eid - element id
   @param cn - array containing all code numbers on element
   
   16.7.2002
*/
void gtopology::give_code_numbers (long eid,long *cn)
{
  long i,j,k,nne,ndofn;
  
  if (cnstate==0){
    print_err("code numbers are required but not available",__FILE__,__LINE__,__func__);
  }
  
  if (gelements[eid].cne==1){
    //  code numbers are stored on elements
    for (i=0;i<give_ndofe (eid);i++){
      cn[i]=gelements[eid].cn[i];
    }
  }
  else{
    //  code numbers are obtained from nodes
    k=0;
    nne=give_nne (eid);
    ivector nod(ASTCKIVEC(nne));
    give_nodes (eid,nod);
    for (i=0;i<nne;i++){
      ndofn=give_ndofn (nod[i]);
      for (j=0;j<ndofn;j++){
	cn[k]=give_dof (nod[i],j);
	k++;
      }
    }
  }
  
}

/**
   function returns array cn stored on the eid-th element
   
   @param eid - element id
   @param cn - array containing all code numbers on element

   16.7.2012
*/
void gtopology::give_cn (long eid,long *cn)
{
  long i;

  if (cnstate==0){
    print_err("code numbers are required but not available",__FILE__,__LINE__,__func__);
  }

  if (gelements[eid].cne==1 || gelements[eid].cne==2){
    //  code numbers are stored on elements
    for (i=0;i<give_ndofe (eid);i++){
      cn[i]=gelements[eid].cn[i];
    }
  }
}

/**
   function extracts all code numbers of actual element
   including code numbers of Lagrange multipliers
   
   @param eid - element id
   @param cn - array containing all code numbers on element
   
   16.7.2002
*/
void gtopology::give_gcode_numbers (long eid,long *cn)
{
  long i,j,k,nne,ndofe,ndofn,lnid,lid,nmult;
  ivector nod;

  if (cnstate==0){
    //fprintf (stderr,"\n\n code numbers are required but not available (file %s, line %d)\n",__FILE__,__LINE__);
  }

  ndofe = give_ndofe (eid);
  
  if (gelements[eid].cne==1){
    //  code numbers are stored on elements
    for (i=0;i<ndofe;i++){
      cn[i]=gelements[eid].cn[i];
    }
  }
  else{
    //  code numbers are obtained from nodes
    k=0;
    nne=give_nne (eid);
    reallocv (RSTCKIVEC(nne,nod));
    give_nodes (eid,nod);
    for (i=0;i<nne;i++){
      ndofn=give_ndofn (nod[i]);
      for (j=0;j<ndofn;j++){
	cn[k]=give_dof (nod[i],j);
	k++;
      }
    }
  }
  
  if (gelements[eid].nmult>0){
    k=ndofe;
    nne=give_nne (eid);
    reallocv (RSTCKIVEC(nne,nod));
    give_nodes (eid,nod);
    for (i=0;i<nne;i++){
      lnid=unln[nod[i]];
      lid=unnl[nod[i]];
      nmult=give_nmult (lnid);
      give_mult_code_numbers (lnid,lid,cn+k);
      k+=nmult;
    }
    for (i=0;i<nne;i++){
      lnid=unln[nod[i]];
      lid=unnl[nod[i]]+1;
      nmult=give_nmult (lnid);
      give_mult_code_numbers (lnid,lid,cn+k);
      k+=nmult;
    }
  }
}

/**
   function extracts code numbers of actual node
   
   @param nid - node id
   @param cn - array containing code numbers on element
   
   7.12.2002
*/
void gtopology::give_node_code_numbers (long nid,long *cn)
{
  long i,ndofn;
  
  if (cnstate==0){
    print_err("code numbers are required but not available", __FILE__, __LINE__, __func__);
  }

  ndofn = give_ndofn (nid);
  for (i=0;i<ndofn;i++){
    cn[i]=gnodes[nid].cn[i];
  }
}



/**
  The function extracts code numbers of the given node, hanging nodes are taken into account.
   
  @param nid - id of required node
  @param cn - array containing code numbers on element, its dimension should be either
              ndofn for ordinary nodes or total number of ndofn collected from all master 
              nodes in the case of the hanging node. The array is reallocated and thus it 
              NEED NOT to be allocated on input. If not allocated or determined dimension 
              is greater then the cn size, the DYNAMIC memory must be allocated in order 
              to prevent stack deallocation after return from the function.
   
  Created by Tomas Koudelka, 05.2023
*/
void gtopology::give_gnode_code_numbers (long nid, ivector &cn)
{
  if (cnstate==0){
    print_err("code numbers are required but not available", __FILE__, __LINE__, __func__);
  }
  
  long ndofn = give_ndofn (nid);
  if (ndofn > 0){
    reallocv(ndofn, cn);
    for (long i=0; i<ndofn; i++){
      cn[i]=gnodes[nid].cn[i];
    }
  }
  else{
    long tndofn = 0;
    long nmn = -ndofn;
    for (long i=0; i<nmn; i++){
      tndofn += give_ndofn(gnodes[nid].mnodes[i]);
    }
    reallocv(tndofn, cn);
    long k=0;
    for (long i=0; i<nmn; i++){
      long nmid = gnodes[nid].mnodes[i];
      ndofn = give_ndofn(nmid);
      for(long j=0; j<ndofn; j++){
        cn[k]=gnodes[nmid].cn[j];
        k++;
      }
    }
  }
}



/**
   function extracts code numbers of multipliers of actual node
   
   @param nid - node id
   @param lid - layer id
   @param cn - array containing code numbers on element
   
   7.12.2002
*/
void gtopology::give_mult_code_numbers (long nid,long lid,long *cn)
{
  long i,nmult;
  
  if (cnstate==0){
    print_err("code numbers are required but not available", __FILE__, __LINE__, __func__);
  }

  nmult = give_nmult (nid);
  for (i=0;i<nmult;i++){
    cn[i]=lgnodes[nid].cn[lid][i];
  }
}

/**
   function assembles code numbers of nodes and multipliers defined on selected end node
   
   @param enid - end node id
   @param ncn1 - array containing code numbers of the first node (on one side of edge)
   @param ncn2 - array containing code numbers of the second node (on the other side of edge)
   @param mcn - array containing code numbers of multipliers
   
   JK, 6.8.2008
*/
void gtopology::give_endnode_code_numbers (long enid,long *ncn1,long *ncn2,long *mcn)
{ 
  ivector nid(ASTCKIVEC(2));

  //  node number
  endnodes[enid].give_node_numbers (nid.a);
  
  //  number of DOFs on node
  // long ndofn = give_ndofn (nid[0]);
  
  //  code numbers of nodes
  give_node_code_numbers (nid[0],ncn1);
  give_node_code_numbers (nid[1],ncn2);
  
  //  code numbers of Lagrange multipliers
  endnodes[enid].give_mult_code_numbers (mcn);
}

/**
   function assembles code numbers of nodes and multipliers defined on selected edge
   
   @param edid - edge id
   @param fln - first or last nodes indicator
          fln=1 - first nodes are used
	  fln=2 - last nodes are used
   @param ncn1 - array containing code numbers of the first node (on one side of edge)
   @param ncn2 - array containing code numbers of the second node (on the other isde of edge)
   @param mcn - array containing code numbers of multipliers
   
   JK, 6.8.2008
*/
void gtopology::give_edge_code_numbers (long edid,long fln,long *ncn1,long *ncn2,long *mcn)
{
  ivector nid(ASTCKIVEC(2));

  if (fln==1){
    gedges[edid].give_first_node_numbers (nid.a);
  }
  if (fln==2){
    gedges[edid].give_last_node_numbers (nid.a);
  }
  
  //  number of DOFs on node
  // long ndofn = give_ndofn (nid[0]);
  
  //  code numbers of nodes
  give_node_code_numbers (nid[0],ncn1);
  give_node_code_numbers (nid[1],ncn2);
  
  gedges[edid].give_mult_code_numbers (fln,mcn);
}



/**
   function returns type of finite element
   
   @param eid - element id
   
   JK, 13.4.2007
*/
gelemtype gtopology::give_elem_type (long eid)
{
  return gelements[eid].get;
}

/**
   function returns type of finite element according to siftop class
   
   @param eid - element id
   
   TKo, 07.2020
*/
gtypel gtopology::give_siftopelem_type (long eid)
{
  switch(gelements[eid].get){
    case linbar:
      return isolinear1d;
    case quadbar:
      return isoquadratic1d;
    case lintriag:
      return trianglelinear;
    case quadtriag:
      return trianglequadratic;
    case linquad:
      return isolinear2d;
    case quadquad:
      return isoquadratic2d;
    case cubicquad:
      return isocubic2d;
    case lintetra:
      return tetrahedronlinear;
    case quadtetra:
      return tetrahedronquadratic;
    case linhexa:
      return isolinear3d;
    case quadhexa:
      return isoquadratic3d;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);

  }
  return noel;
}



/**
   function selects end nodes on required 1D element
   
   @param eid - element id
   @param nodes - end node numbers
   
   JK, 12.10.2008
*/
void gtopology::give_end_nodes (long eid,long *nodes)
{
  long nne;
  gelemtype get;
  ivector nod;
  ivector enod(ASTCKIVEC(2));
  
  //  number nodes on element
  nne = give_original_nne (eid);
  //  allocation of the array of node numbers
  reallocv (RSTCKIVEC(nne,nod));
  //  node numbers
  give_original_nodes (eid,nod);
  //  type of finite element
  get = give_elem_type (eid);

  switch (get){
  case linbar:{
    linbar_endpoints (enod.a);
    nodes[0]=nod[enod[0]];
    nodes[1]=nod[enod[1]];
    break;
  }
  case quadbar:{
    quadbar_endpoints (enod.a);
    nodes[0]=nod[enod[0]];
    nodes[1]=nod[enod[1]];
    break;
  }

  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function returns local node numbers for the given edge of element.

  @param eid - element id
  @param edid - edge id 
  @param edgenod - %vector of local node id

  @return The function returns output in parameter nodes.

  Created by Tomas Koudelka, 3.12.2013
*/
void gtopology::give_edge_loc_nodes (long eid,long edid,ivector &edgenod)
{
  gelemtype get;
  
  //  type of finite element
  get = give_elem_type (eid);
  
  switch (get)
  {
    case linbar:
      edgenod[0] = 0;
      edgenod[1] = 1;
      break;
    case quadbar:
      edgenod[0] = 0;
      edgenod[1] = 1;
      edgenod[2] = 2;
      break;
    case lintriag:
      lintriangle_edgnod(edgenod.a,edid);
      break;
    case quadtriag:
      quadtriangle_edgnod(edgenod.a,edid);
      break;
    case linquad:
      linquadrilat_edgnod (edgenod.a,edid);
      break;
    case quadquad:
      quadquadrilat_edgnod (edgenod.a,edid);
      break;
    case lintetra:
      lintetrahedral_edgnod(edgenod.a,edid);
      break;
    case quadtetra:
      quadtetrahedral_edgnod (edgenod.a,edid);
      break;
    case linhexa:
      linhexahedral_edgnod(edgenod.a,edid);
      break;
    case quadhexa:
      quadhexahedral_edgnod (edgenod.a,edid);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
   function selects nodes on required edge
   
   @param eid - element id
   @param edid - edge id
   @param nodes - edge node numbers
   
   JK, 13.4.2007
*/
void gtopology::give_edge_nodes (long eid,long edid,long *nodes)
{
  long i,nne,nned;
  gelemtype get;
  ivector nod;
  ivector edgenod;
  
  //  number nodes on element
  nne = give_original_nne (eid);
  //  allocation of the array of node numbers
  reallocv (RSTCKIVEC(nne,nod));
  //  node numbers
  give_original_nodes (eid,nod);
  //  type of finite element
  get = give_elem_type (eid);
  // number of nodes on element edge
  nned = give_nned (eid);
  
  reallocv(RSTCKIVEC(nned, edgenod));
  switch (get){
  case lintriag:{
    lintriangle_edgnod(edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case quadtriag:{
    quadtriangle_edgnod(edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case linquad:{
    linquadrilat_edgnod (edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case quadquad:{
    quadquadrilat_edgnod (edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case lintetra:{
    lintetrahedral_edgnod(edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case quadtetra:{
    quadtetrahedral_edgnod (edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case linhexa:{
    linhexahedral_edgnod(edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  case quadhexa:{
    quadhexahedral_edgnod (edgenod.a,edid);
    for(i = 0; i < nned; i++){
      nodes[i]=nod[edgenod[i]];
    }
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function returns node numbers on required element surface
   
   @param eid - element id
   @param edid - edge id
   @param nodes - array containing edge nodes

   JB, 18.4.2007
*/
void gtopology::give_surf_nodes (long eid,long surfid,long *nodes)
{
  long i,nne,nnsurf;
  gelemtype get;
  ivector nod;
  ivector surfnod;
  
  //  number nodes on element
  nne = give_nne (eid);
  //  allocation of the array of node numbers
  reallocv (RSTCKIVEC(nne,nod));
  //  node numbers
  give_nodes (eid,nod);
  //  type of finite element
  get = give_elem_type (eid);
  //  number of nodes on surface
  nnsurf = give_nnsurf (eid);
  
  reallocv(RSTCKIVEC(nnsurf, surfnod));
  switch(get){
  case lintetra:{
    lintetrahedral_surfnod (surfnod.a,surfid);
    for(i = 0; i < nnsurf; i++){
      nodes[i]=nod[surfnod[i]];
    }
    break;
  }
  case quadtetra:{
    quadtetrahedral_surfnod (surfnod.a,surfid);
    for(i = 0; i < nnsurf; i++){
      nodes[i]=nod[surfnod[i]];
    }
    break;
  }
  case linhexa:{
    linhexahedral_surfnod (surfnod.a,surfid);
    for(i = 0; i < nnsurf; i++){
      nodes[i]=nod[surfnod[i]];
    }
    break;
  }
  case quadhexa:{
    quadhexahedral_surfnod (surfnod.a,surfid);
    for(i = 0; i < nnsurf; i++){
      nodes[i]=nod[surfnod[i]];
    }
    break;
  }
    
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function returns number of end nodes on the element
   
   @param eid - element id
   
   JK, 17.10.2008
*/
long gtopology::give_nen (long eid)
{
  long ned=0;
  gelemtype get;
  
  get = give_elem_type (eid);
  
  switch (get){
  case linbar:{     ned=2;    break;  }
  case quadbar:{    ned=2;    break;  }
  case lintriag:{   ned=0;    break;  }
  case quadtriag:{  ned=0;    break;  }
  case linquad:{    ned=0;    break;  }
  case quadquad:{   ned=0;    break;  }
  case lintetra:{   ned=0;    break;  }
  case quadtetra:{  ned=0;    break;  }
  case linhexa:{    ned=0;    break;  }
  case quadhexa:{   ned=0;    break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return ned;
}

/**
   function returns number of edges of the element
   
   @param eid - element id
   
   JB, 18.4.2007
*/
long gtopology::give_ned (long eid)
{
  long ned=0;
  gelemtype get;
  
  get = give_elem_type (eid);
  
  switch (get){
  case linbar:{     ned=0;    break;  }
  case quadbar:{    ned=0;    break;  }
  case lintriag:{   ned=3;    break;  }
  case quadtriag:{  ned=3;    break;  }
  case linquad:{    ned=4;    break;  }
  case quadquad:{   ned=4;    break;  }
  case lintetra:{   ned=6;    break;  }
  case quadtetra:{  ned=6;    break;  }
  case linhexa:{    ned=12;   break;  }
  case quadhexa:{   ned=12;   break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return ned;
}

/**
   function returns number of nodes on one element edge
   
   @param eid - number of element
   
   JB, 18.4.2007
*/
long gtopology::give_nned (long eid)
{
  long nned=0;
  gelemtype get;

  get = give_elem_type (eid);
  
  switch (get){
  case linbar:{     nned=2;    break;  }
  case quadbar:{    nned=3;    break;  }
  case lintriag:{   nned=2;    break;  }
  case quadtriag:{  nned=3;    break;  }
  case linquad:{    nned=2;    break;  }
  case quadquad:{   nned=3;    break;  }
  case lintetra:{   nned=2;    break;  }
  case quadtetra:{  nned=3;    break;  }
  case linhexa:{    nned=2;    break;  }
  case quadhexa:{   nned=3;    break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return nned;
}

/**
   function returns number of surfaces of the element
   
   @param eid - element id
   
   JB, 18.4.2004
*/
long gtopology::give_nsurf (long eid)
{
  long nsurf;
  gelemtype get;

  get = give_elem_type (eid);
  
  switch (get){
  case lintetra:{   nsurf=4;    break;  }
  case quadtetra:{  nsurf=4;    break;  }
  case linhexa:{    nsurf=6;    break;  }
  case quadhexa:{   nsurf=6;    break;  }
  default:{
    nsurf = 0;
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return nsurf;
}

/**
   function returns number of nodes on element surface
   
   @param eid - element id
   
   JB, 18.4.2004
*/
long gtopology::give_nnsurf (long eid)
{
  long nnsurf;
  gelemtype get;
  
  get = give_elem_type (eid);
  
  switch (get){
  case lintetra:{         nnsurf=3;    break;  }
  case linhexa:{          nnsurf=4;    break;  }
  case quadtetra:{        nnsurf=6;    break;  }
  case quadhexa:{         nnsurf=8;    break;  }
  default:{
    nnsurf = 0;
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return nnsurf;
}


/**
   function returns degree of polynomial approximation
   
   @param eid - number of element
   
   JB, 18.4.2007
*/
long gtopology::give_degree (long eid)
{
  long deg=0;
  gelemtype get;
  
  get = give_elem_type (eid);
  
  switch (get){
  case linbar:{     deg=1;    break;  }
  case quadbar:{    deg=2;    break;  }
  case lintriag:{   deg=1;    break;  }
  case quadtriag:{  deg=2;    break;  }
  case linquad:{    deg=1;    break;  }
  case quadquad:{   deg=2;    break;  }
  case lintetra:{   deg=1;    break;  }
  case quadtetra:{  deg=2;    break;  }
  case linhexa:{    deg=1;    break;  }
  case quadhexa:{   deg=2;    break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  return deg;
}

/**
   Function return spatial dimension of the element
   @param eid - number of element
   16.11.10 JB
**/
long gtopology::give_whole_dim(long eid)
{
  
  gelemtype get;
  long dim;
  get =  give_elem_type (eid);
  switch(get){
  case linbar:{     dim = 1;    break;  }
  case quadbar:{    dim = 1;    break;  }
  case lintriag:{   dim = 2;    break;  }
  case quadtriag:{  dim = 2;    break;  }
  case linquad:{    dim = 2;    break;  }
  case quadquad:{   dim = 2;    break;  }
  case lintetra:{   dim = 3;    break;  }
  case quadtetra:{  dim = 3;    break;  }
  case linhexa:{    dim = 3;    break;  }
  case quadhexa:{   dim = 3;    break;  }
  default:{
    dim = 0;
    print_err("unknown element type is required\n", __FILE__, __LINE__, __func__);
  }
  }
  return (dim);
}



/**
  The function determines the maximum dimension of elements used in the problem.

  @return The function returns the maximum dimension of elements in the problem.

  Created by Tomas Koudelka, 08.2016  
*/
long gtopology::give_maxdimension ()
{
  long i, aux, ret = 0L;

  for (i=0L; i<ne; i++)
  {
    aux = give_whole_dim(i);
    if (ret < aux)
      ret = aux;
  }
  return ret;
}



/**
   function generates code numbers (ordering of unknowns and equations)
   function possibly reorders node numbers
   
   @param out - output file
   
   JK, 28.2.2007
*/
long gtopology::codenum_generation (FILE *out)
{
  long i,ndof;
  
  ordering = new long [nn];
  
  //  switch over renumbering type
  switch (nodren){
  case no_renumbering:{
    for (i=0;i<nn;i++){
      ordering[i]=i;
    }
    break;
  }
  case cuthill_mckee:
  case rev_cuthill_mckee:{
    cuthill_mckee_renumb (out);
    break;
  }
  case sloan:{
    sloan_renumb(out);
    break;
  }
  default:{
    print_err("unknown type of renumbering is required",__FILE__,__LINE__,__func__);
  }
  }
  
  //  no code numbers generation
  if (cngen==0){
    ndof = 0;
  }
  //  generation of the code numbers
  if (cngen==1){
    //  the classical generation of code numbers
    ndof = gencodnum ();
  }
  if (cngen==2){
    //  function allocates and perform back up of the array cn on nodes
    //backup_cn ();
    
    //  internal unknowns are ordered first, interface unknowns are ordered last
    ndof = schur_ordering ();
  }
  if (cngen==3){
    //  primal unknowns are ordered first, dual unknowns are ordered last
    //  in more detail: internal unknowns are ordered first,
    //  interface unknowns are ordered then and finally, dual unknowns are ordered last
    ndof = saddlepoint_ordering (stop->ns,stop->nnsd,stop->ltg);
  }
  
  delete [] ordering;
  ordering = NULL;
  
  if (hangnod==1){
    //  there are hanging nodes
    //  DOFs have to be copied from nodes to elements
    dof_transfer_hangnod ();
  }

  //  codenum_elem_print(out);
  //  print_dof_diff(out);
  
  return ndof;
}



/**
  The function prints the content of code number arrays of nodes to the opened text file.

  @param out - pointer to the opned text file where the data will be printed out

  @return The funtion does not return anything.

  Created by Tomas Koudelka, 08.2016 
*/
void gtopology::codenum_print(FILE *out)
{
  long i, j;
  for(i=0L; i<nn; i++)
  {
    fprintf(out, "%ld : ", i+1);
    for(j=0; j<give_ndofn(i); j++)
      fprintf(out, "%ld ", give_dof(i, j));
    fprintf(out, "\n");
  }
}



/**
  The function prints code numbers of element to the opened text file.

  @param out - pointer to the opned text file where the data will be printed out

  @return The funtion does not return anything.

  Created by Tomas Koudelka, 10.2022 
*/
void gtopology::codenum_elem_print(FILE *out)
{
  long i, j, k, ndofn, ndofe;
  ivector dofe, enod;
  for(i=0L; i<ne; i++)
  {
    ndofe = gelements[i].give_ndofe();
    fprintf(out, "%ld : ", i+1);
    if (gelements[i].give_cne()){
      // pozor reference se muze premazat ve vetvi else pro nasledujic krok !!!??? 
      //      makerefv(dofe, gelements[i].cn, ndofe);
    }
    else{
      reallocv(ndofe, dofe);
      long nne = gelements[i].give_nne();
      makerefv(enod, gelements[i].nodes, nne);
      for(j=0; j<nne; j++){
        ndofn = give_ndofn(enod(j));
        for(k=0; k<ndofn; k++)
          fprintf(out, "%ld ", give_dof(enod(j), k));
      }
    }
    fprintf(out, "\n");
  }
}



/**
  The function prints code numbers of element to the opened text file.

  @param out - pointer to the opned text file where the data will be printed out

  @return The funtion does not return anything.

  Created by Tomas Koudelka, 10.2022 
*/
void gtopology::print_dof_diff(FILE *out)
{
  long i, j, k, l, ndofn, ndofe;
  ivector dofe, enod, dofmin(nn), dofmax(nn);
  fillv(LONG_MAX, dofmin);
  
  for(i=0L; i<ne; i++)
  {
    ndofe = gelements[i].give_ndofe();
    if (gelements[i].give_cne()){
      // pozor reference se muze premazat ve vetvi else pro nasledujic krok !!!??? 
      //      makerefv(dofe, gelements[i].cn, ndofe);
    }
    else{
      reallocv(ndofe, dofe);
      long nne = gelements[i].give_nne();
      makerefv(enod, gelements[i].nodes, nne);
      for(j=0; j<nne; j++){
        ndofn = give_ndofn(enod(j));
        for(k=0; k<ndofn; k++){
          long d = give_dof(enod(j), k);
          for (l=0; l<nne; l++){
            if ((dofmin(enod(l)) > d) && (d >= 0))
              dofmin(enod(l)) = d;
            if ((dofmax(enod(l)) < d) && (d >= 0))
              dofmax(enod(l)) = d;
          }
        }
      }
    }
  }
  fprintf(out, "\n");
  long min=LONG_MAX, max=0, imin=0, imax=0;
  for (i=0; i<nn; i++){
    long d = dofmax(i)-dofmin(i);
    fprintf(out, "N%ld : d=%ld, min=%ld, max=%ld\n", i+1, d, dofmin(i), dofmax(i));
    if (d < min)  min = d, imin = i;
    if (d > max)  max = d, imax = i;      
  }
  fprintf(out, "\nDelta DOF max = %ld at node %ld\n"
          "Delta DOF min = %ld at node %ld\n", max, imax+1, min, imin+1);
}



/**
   function generates code numbers (ordering of unknowns and equations)
   
   JK, 25.6.2001
*/
long gtopology::gencodnum ()
{
  long i,j,k,ii,ndofn;
  ivector aux;
  
  //  searching of maximum code number
  ndof=0;
  for (i=0;i<nn;i++){
    ndofn=give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=give_dof (i,j);
      if (k>ndof)  
        ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  reallocv(RSTCKIVEC(ndof, aux));
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }
  
  //  contributions from nodes
  ndof=1;
  for (i=0;i<nn;i++){
    ii=ordering[i];
    ndofn=give_ndofn (ii);
    for (j=0;j<ndofn;j++){
      k=give_dof (ii,j);
      if (k<0)  continue;
      if (k==0)  continue;
      if (k==1){
	gnodes[ii].cn[j]=ndof;  ndof++;
      }
      if (k>1){
	if (aux[k-2]==-1){
	  gnodes[ii].cn[j]=ndof;
	  aux[k-2]=ndof;
	  ndof++;
	}
	else{
	  gnodes[ii].cn[j]=aux[k-2];
	}
      }
    }
  }
  ndof--;
  
  
  //  contributions from elements
  for (i=0;i<ne;i++){
    if (gelements[i].cne==1){
      for (j=0;j<give_ndofe (i);j++){
	if (gelements[i].cn[j]>ndof)  
          ndof=gelements[i].cn[j];
      }
    }
  }
  
  //  state of code numbers is changed
  cnstate=1;
  
  return ndof;
}

/**
   function generates code numbers (numbers of equations)
   with respect to Schur ordering
   
   internal nodes are ordered first, interface nodes are ordered at the end
   
   interface nodes are indicated by array ltg which has value 1
   
   JK, 22.10.2007
*/
long gtopology::schur_ordering ()
{
  long i,j,k,l,ii,jj,kk,ndofn;
  ivector aux;
  
  //  searching of maximum code number
  ndof=0;
  for (i=0;i<nn;i++){
    ndofn=give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=give_dof (i,j);
      if (k>ndof)  ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  reallocv(RSTCKIVEC(ndof, aux));
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }

  //  ordering of internal nodes
  jj=-1;  kk=-1;
  ndof=1;
  for (l=0;l<ns;l++){
    for (i=0;i<stop->nnsd[l];i++){
      jj++;
      if (stop->ltg[l][i]==-1){
	ii=ordering[jj];
	ndofn=give_ndofn (ii);
	for (j=0;j<ndofn;j++){
	  k=give_dof (ii,j);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    gnodes[ii].cn[j]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      gnodes[ii].cn[j]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      gnodes[ii].cn[j]=aux[k-2];
	    }
	  }
	}
      }
    }

    //  number of internal DOFs
    nidof=ndof-1;
    
    //  ordering of interface nodes
    jj=kk;
    for (i=0;i<stop->nnsd[l];i++){
      jj++;
      if (stop->ltg[l][i]>-1){
	ii=ordering[jj];
	ndofn=give_ndofn (ii);
	for (j=0;j<ndofn;j++){
	  k=give_dof (ii,j);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    gnodes[ii].cn[j]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      gnodes[ii].cn[j]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      gnodes[ii].cn[j]=aux[k-2];
	    }
	  }
	}
      }
    }
    kk+=stop->nnsd[l];
  }
  ndof--;
  
  //  contributions from elements
  for (i=0;i<ne;i++){
    if (gelements[i].cne==1){
      for (j=0;j<give_ndofe (i);j++){
	if (gelements[i].cn[j]>ndof)
	  ndof=gelements[i].cn[j];
      }
    }
  }
  
  //  number of boundary/interface DOFs
  nbdof=ndof-nidof;
  
  //  state of code numbers is changed
  cnstate=1;
  
  return ndof;

  /*
  long i,j,k,ii,ndofn,*aux;
  
  //  searching of maximum code number
  ndof=0;
  for (i=0;i<nn;i++){
    ndofn=give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=give_dof (i,j);
      if (k>ndof)  ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  aux = new long [ndof];
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }

  //  ordering of internal nodes
  ndof=1;
  for (i=0;i<nn;i++){
    if (ltg[0][i]==-1){
      ii=ordering[i];
      ndofn=give_ndofn (ii);
      for (j=0;j<ndofn;j++){
	k=give_dof (ii,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  gnodes[ii].cn[j]=ndof;  ndof++;
	}
	if (k>1){
	  if (aux[k-2]==-1){
	    gnodes[ii].cn[j]=ndof;
	    aux[k-2]=ndof;
	    ndof++;
	  }
	  else{
	    gnodes[ii].cn[j]=aux[k-2];
	  }
	}
      }
    }
  }
  //  number of internal DOFs
  nidof=ndof-1;
  
  //  ordering of interface nodes
  for (i=0;i<nn;i++){
    if (ltg[0][i]==0){
      ii=ordering[i];
      ndofn=give_ndofn (ii);
      for (j=0;j<ndofn;j++){
	k=give_dof (ii,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  gnodes[ii].cn[j]=ndof;  ndof++;
	}
	if (k>1){
	  if (aux[k-2]==-1){
	    gnodes[ii].cn[j]=ndof;
	    aux[k-2]=ndof;
	    ndof++;
	  }
	  else{
	    gnodes[ii].cn[j]=aux[k-2];
	  }
	}
      }
    }
  }
  ndof--;
  
  delete [] aux;
  
  //  contributions from elements
  for (i=0;i<ne;i++){
    if (gelements[i].cne==1){
      for (j=0;j<give_ndofe (i);j++){
	if (gelements[i].cn[j]>ndof)
	  ndof=gelements[i].cn[j];
      }
    }
  }
  
  //  number of boundary/interface DOFs
  nbdof=ndof-nidof;
  
  //  state of code numbers is changed
  cnstate=1;
  
  return ndof;
*/

}



/**
   function generates code numbers (numbers of equations)
   with respect to saddle point problem
   
   internal nodes are ordered first, interface nodes are ordered at the end
   finally, Lagrange multipliers are defined by the function codenum_multip

   interface nodes are indicated by array ltg which has value 1
   
   @param ns - number of subdomains
   @param nnsd - number of nodes on subdomains
   @param ltg - local to global map
   
   JK, 13.12.2007
*/
long gtopology::saddlepoint_ordering (long ns,long *nnsd,long **ltg)
{
  long i,j,k,ii,jj,nid,ndofn;
  ivector aux;
  
  //  searching of maximum code number
  ndof=0;
  for (i=0;i<nn;i++){
    ndofn=give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=give_dof (i,j);
      if (k>ndof)  ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  reallocv(RSTCKIVEC(ndof, aux));
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }

  //  ordering of internal nodes
  nid=0;
  ndof=1;
  for (jj=0;jj<ns;jj++){
    for (i=0;i<nnsd[jj];i++){
      if (ltg[jj][i]==-1){
	ii=ordering[nid];
	ndofn=give_ndofn (ii);
	for (j=0;j<ndofn;j++){
	  k=give_dof (ii,j);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    gnodes[ii].cn[j]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      gnodes[ii].cn[j]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      gnodes[ii].cn[j]=aux[k-2];
	    }
	  }
	}
      }
      nid++;
    }
  }
  //  number of internal DOFs
  nidof=ndof-1;
  
  //  ordering of interface nodes
  nid=0;
  for (jj=0;jj<ns;jj++){
    for (i=0;i<nnsd[jj];i++){
      if (ltg[jj][i]!=-1){
	ii=ordering[nid];
	ndofn=give_ndofn (ii);
	for (j=0;j<ndofn;j++){
	  k=give_dof (ii,j);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    gnodes[ii].cn[j]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      gnodes[ii].cn[j]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      gnodes[ii].cn[j]=aux[k-2];
	    }
	  }
	}
      }
      nid++;
    }
  }
  //ndof--;
  
  //  contributions from elements
  for (i=0;i<ne;i++){
    if (gelements[i].cne==1){
      for (j=0;j<give_ndofe (i);j++){
	if (gelements[i].cn[j]>ndof)
	  ndof=gelements[i].cn[j];
      }
    }
  }
  
  //  number of boundary/interface DOFs
  nbdof=ndof-nidof;
  
  //  generation of code numbers of Lagrange multipliers
  //  nsad is the number of unknowns in saddle point problem
  //nsad = codenum_multip ();
  ndof = codenum_multip ();

  ndof--;

  //  state of code numbers is changed
  cnstate=1;
  
  fprintf (stdout,"\n number of DOFs            %ld",ndof);
  fprintf (stdout,"\n number of internal DOFs   %ld",nidof);
  fprintf (stdout,"\n number of interface DOFs  %ld",nbdof);
  fprintf (stdout,"\n number of unknowns in sadd. point. prob.  %ld",nsad);

  
  return ndof;
}



/**
   function creates general topology from two identical topologies (meshes)
   with different unknowns defined in nodes
   
   function is also usefull for two different meshes used in accordance with
   Babuska-Brezzi condition; in such case, one mesh contains elements with
   quadratic approximation functions and the second mesh contains only elements
   with linear approximation functions; the first mesh contains therefore more
   nodes than the second mesh; user can add all nodes from the first mesh which
   are not defined in the second mesh to the list of nodes of the second mesh
   but they are not mentioned in element-node correspondence
   
   important notice: topology describing the mesh containing elements with
   quadratic approximation functions must be send to gt1 while topology
   describing the mesh containing elements with linear approximation functions
   must be send to gt2
   
   @param gt1 - pointer to the first topology
   @param gt2 - pointer to the second topology
   
   JK, 28.10.2002
*/
void gtopology::comptop (gtopology *gt1,gtopology *gt2)
{
  long i,j,k,l,nn1,nn2,ne1,ne2,ndofn1,nne,nne1,nne2;
  ivector enod1, enod2;
  
  //  meshes checking - nodes
  nn1 = gt1->nn;
  nn2 = gt2->nn;
  if (nn1 < nn2){
    print_err("mechanical topology has lesser number of nodes (%ld) than transport topology (%ld)", __FILE__, __LINE__, __func__, nn1, nn2);
    abort ();
  }
  
  nn = nn1;
  gnodes = new gnode [nn];
  
  //  assigment of coordinates to the unified mesh
  //  assigment of numbers of DOFs from the first mesh to the unified mesh
  for (i=0;i<nn;i++){
    ndofn1 = gt1->gnodes[i].ndofn;
    gnodes[i].ndofn = ndofn1;
    
    gnodes[i].x=gt1->gnodes[i].x;
    gnodes[i].y=gt1->gnodes[i].y;
    gnodes[i].z=gt1->gnodes[i].z;
  }
  
  //  meshes checking - elements
  ne1 = gt1->ne;
  ne2 = gt2->ne;

  // puvodni podminka
  /*  if (ne1!=ne2){
    print_err("two topologies have different numbers of elements", __FILE__, __LINE__, __func__);
    abort ();
    }*/

  // uprava pro BEACON ???!!!
  ne = ne1;
  if (ne1 > ne2)
    ne = ne2;
  // konec upravy pro BEACON ???!!!

  gelements = new gelement [ne];
  
  //  assigment of numbers of DOFs from the second mesh to the unified mesh
  for (i=0;i<nn2;i++){
    gt2->gnodes[i].ai=0;
  }
  for (i=0;i<ne;i++){
    nne1 = gt1->give_nne (i);
    nne2 = gt2->give_nne (i);
    if (nne1 >= nne2)
    { 
      reallocv (RSTCKIVEC(nne1,enod1));
      gt1->give_nodes (i,enod1);
      reallocv (RSTCKIVEC(nne2,enod2));
      gt2->give_nodes (i,enod2);
      for (j=0;j<nne2;j++){
        if (gt2->gnodes[enod2[j]].ai==1)  continue;
        gnodes[enod1[j]].ndofn+=gt2->gnodes[enod2[j]].ndofn;
        gt2->gnodes[enod2[j]].ai=1;
      }
    }
    else
    {
      print_err("wrong number of nodes on element %ld (nne_mef=%ld < nne_trf=%ld)", __FILE__, __LINE__, __func__, i+1, nne1, nne2);
      abort();
    }
  }
  
  for (i=0;i<nn2;i++){
    gt2->gnodes[i].ai=0;
  }
  
  //  assigment of code numbers from the first mesh to the unified mesh
  for (i=0;i<nn;i++){
    gnodes[i].cn = new long [gnodes[i].ndofn];
    for (j=0;j<gt1->gnodes[i].ndofn;j++){
      gnodes[i].cn[j]=gt1->gnodes[i].cn[j];
    }
  }
  
  //  assigment of code numbers from the second mesh to the unified mesh
  for (i=0;i<ne;i++){
    nne1 = gt1->give_nne (i);
    nne2 = gt2->give_nne (i);
    if (nne1 >= nne2)
    { 
      reallocv (RSTCKIVEC(nne1,enod1));
      gt1->give_nodes (i,enod1);
      reallocv (RSTCKIVEC(nne2,enod2));
      gt2->give_nodes (i,enod2);
      for (j=0;j<nne2;j++){
        if (gt2->gnodes[enod2[j]].ai==1)  continue;
        l=gt1->gnodes[enod1[j]].ndofn;
        for (k=0;k<gt2->gnodes[enod2[j]].ndofn;k++){
          gnodes[enod1[j]].cn[l]=gt2->gnodes[enod2[j]].cn[k];  l++;
        }
        gt2->gnodes[enod2[j]].ai=1;
      }
    }
    else
    {
      print_err("wrong number of nodes on element %ld (nne_mef=%ld < nne_trf=%ld)", __FILE__, __LINE__, __func__, i+1, nne1, nne2);
      abort();
    }
  }
  
  for (i=0;i<nn2;i++){
    gt2->gnodes[i].ai=0;
  }

  //  assigment of numbers of nodes on element to the unified mesh
  //  assigment of node numbers of elements to the unified mesh
  for (i=0;i<ne;i++){
    gelements[i].ndofe = gt1->gelements[i].ndofe + gt2->gelements[i].ndofe;
    nne = gt1->gelements[i].nne;
    gelements[i].nne = nne;
    gelements[i].nodes = new long [nne];
    for (j=0;j<nne;j++){
      gelements[i].nodes[j] = gt1->gelements[i].nodes[j];
    }
  }
}

/**
   function distributes code numbers from unified topology to the
   local topologies
   
   function is also usefull for two different meshes used in accordance with
   Babuska-Brezzi condition; in such case, one mesh contains elements with
   quadratic approximation functions and the second mesh contains only elements
   with linear approximation functions; the first mesh contains therefore more
   nodes than the second mesh; user can add all nodes from the first mesh which
   are not defined in the second mesh to the list of nodes of the second mesh
   but they are not mentioned in element-node correspondence
   
   important notice: topology describing the mesh containing elements with
   quadratic approximation functions must be send to gt1 while topology
   describing the mesh containing elements with linear approximation functions
   must be send to gt2

   @param gt1 - pointer to first topology
   @param gt2 - pointer to second topology
   
   JK, 28.10.2002
*/
void gtopology::cndistr (gtopology *gt1,gtopology *gt2)
{
  long i,j,k,l,ndofn1;
  long nne1, nne2, ne1, ne2;
  ivector enod1, enod2;
  

  // distribution of code numbers over the first (mechanical topology)
  // it is supposed that gt1->nn == nn
  if (gt1->nn != nn)
  {
    print_err("wrong number of nodes on the mechanical topology (mefel_nn=%ld, metr_nn=%ld)", __FILE__, __LINE__, __func__, gt1->nn, nn);
    abort();
  }
  for (i=0;i<nn;i++){
    ndofn1=gt1->gnodes[i].ndofn;
    if (gnodes[i].ndofn>ndofn1){
      for (j=0;j<ndofn1;j++){
	gt1->gnodes[i].cn[j]=gnodes[i].cn[j];
      }
    }
    else{
      for (j=0;j<ndofn1;j++){
	gt1->gnodes[i].cn[j]=gnodes[i].cn[j];
      }
    }
  }
  // distribution of code numbers over the second (transport topology)
  // it is supposed that gt2->nn <= nn and gt2->ne <= ne
  ne1 = gt1->ne;
  ne2 = gt2->ne;
  if (ne2 > ne)
  {
    print_err("number of elements on the second topology is greater than on the unified one (gt2->ne=%ld > ne=%ld)", __FILE__, __LINE__, __func__, ne2, ne);
    abort();
  }
  if (ne2 > ne1)
  {
    print_err("number of elements on the second topology is greater than on the first one (gt2->ne=%ld > gt1->ne=%ld)", __FILE__, __LINE__, __func__, ne2, ne1);
    abort();
  }
  for (i=0;i<gt2->nn;i++){
    gt2->gnodes[i].ai=0;
  }
  for (i=0;i<ne2;i++){
    nne1 = gt1->give_nne (i);
    nne2 = gt2->give_nne (i);
    if (nne1 >= nne2)
    { 
      reallocv (RSTCKIVEC(nne1,enod1));
      gt1->give_nodes (i,enod1);
      reallocv (RSTCKIVEC(nne2,enod2));
      gt2->give_nodes (i,enod2);
      for (j=0;j<nne2;j++){
        if (gt2->gnodes[enod2[j]].ai==1)  continue;
        l=gt1->gnodes[enod1[j]].ndofn;
        for (k=0;k<gt2->gnodes[enod2[j]].ndofn;k++){
          gt2->gnodes[enod2[j]].cn[k] = gnodes[enod1[j]].cn[l];  l++;
        }
        gt2->gnodes[enod2[j]].ai=1;
      }
    }
    else
    {
      print_err("wrong number of nodes on element %ld (nne_mef=%ld < nne_trf=%ld)", __FILE__, __LINE__, __func__, i+1, nne1, nne2);
      abort();
    }
  }
  /*
  for (i=0;i<nn;i++){
    ndofn1=gt1->gnodes[i].ndofn;
    if (gnodes[i].ndofn>ndofn1){
      for (j=0;j<ndofn1;j++){
	gt1->gnodes[i].cn[j]=gnodes[i].cn[j];
      }
      k=ndofn1;
      ndofn2=gt2->gnodes[i].ndofn;
      for (j=0;j<ndofn2;j++){
	gt2->gnodes[i].cn[j]=gnodes[i].cn[k];
	k++;
      }
    }
    else{
      for (j=0;j<ndofn1;j++){
	gt1->gnodes[i].cn[j]=gnodes[i].cn[j];
      }
    }
  }
  */
}

/**
   function asembles list of adjacent elements
   adjacent elements share a node, an edge or a surface
   
   @param out - output stream
   
   JK
*/
void gtopology::adjacelem (FILE *out)
{
  long i,j,k,min,ind,s,nne,last;
  long **aux;
  ivector nodes;
  
  //  allocation of array of counts of adjacent elements to nodes
  if (nadjelnod==NULL)
    nadjelnod = new long [nn];
  memset (nadjelnod,0,nn*sizeof(long));
  
  //  number of contributions
  for (i=0;i<ne;i++){
    nne = give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      nadjelnod[nodes[j]]++;
    }
  }
  
  //  allocation of array of numbers of adjacent elements to nodes
  if (adjelnod==NULL){
    adjelnod = new long* [nn];
    for (i=0;i<nn;i++){
      adjelnod[i] = new long [nadjelnod[i]];
    }
  }
  for (i=0;i<nn;i++){
    nadjelnod[i] = 0;
  }
  
  
  //  filling of array
  for (i=0;i<ne;i++){
    nne = give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      adjelnod[nodes[j]][nadjelnod[nodes[j]]++] = i;
    }
  }
  
  if (out)
  {
    //  printing of array
    fprintf (out,"\n\n\n NUMBERS OF ELEMENTS ADJACENT TO NODES\n");
    for (i=0;i<nn;i++){
      fprintf (out,"\n node %5ld   %ld",i,nadjelnod[i]);
    }

    fprintf (out,"\n\n\n LIST OF ELEMENTS ADJACENT TO NODES\n");
    for (i=0;i<nn;i++){
      fprintf (out,"\n node %4ld  ",i);
      for (j=0;j<nadjelnod[i];j++){
        fprintf (out,"  %ld",adjelnod[i][j]);
      }
    }
    fprintf (out,"\n\n\n");
  }
  
  
  //  allocation of array of counts of adjacent elements to elements
  nadjelel = new long [ne];
  memset (nadjelel,0,ne*sizeof(long));
  
  //  number of contributions
  for (i=0;i<ne;i++){
    nne=give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      nadjelel[i]+=nadjelnod[nodes[j]];
    }
  }
  
  //  allocation of the auxiliary array
  aux = new long* [ne];
  for (i=0;i<ne;i++){
    aux[i] = new long [nadjelel[i]];
  }
  
  
  //  filling of array
  memset (nadjelel,0,ne*sizeof(long));
  for (i=0;i<ne;i++){
    nne=give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      for (k=0;k<nadjelnod[nodes[j]];k++){
	aux[i][nadjelel[i]]=adjelnod[nodes[j]][k];
	nadjelel[i]++;
      }
    }
  }
  
  //  sorting of array
  for (i=0;i<ne;i++){
    last=-1;
    for (j=0;j<nadjelel[i];j++){
      min=LONG_MAX;
      for (k=j;k<nadjelel[i];k++){
	if (min>aux[i][k]){
	  min=aux[i][k];  ind=k;
	}
      }
      if (last==min){
	aux[i][ind]=aux[i][nadjelel[i]-1];
	nadjelel[i]--;  j--;
      }
      else{
	s=aux[i][j];
	aux[i][j]=min;
	aux[i][ind]=s;
	last=min;
      }
    }
  }
  
  //  allocation of array of numbers of adjacent elements to elements
  adjelel = new long* [ne];
  for (i=0;i<ne;i++){
    adjelel[i] = new long [nadjelel[i]];
    memset (adjelel[i],0,nadjelel[i]*sizeof(long));
  }
  
  for (i=0;i<ne;i++){
    for (j=0;j<nadjelel[i];j++){
      adjelel[i][j]=aux[i][j];
    }
  }
  
  if (out)
  {
    //  printing of array
    fprintf (out,"\n\n\n NUMBERS OF ELEMENTS ADJACENT TO ELEMENTS\n");
    for (i=0;i<ne;i++){
      fprintf (out,"\n element %5ld   %ld",i,nadjelel[i]);
    }
    
    fprintf (out,"\n\n\n LIST OF ELEMENTS ADJACENT TO ELEMENTS\n");
    for (i=0;i<ne;i++){
      fprintf (out,"\n element %4ld  ",i);
      for (j=0;j<nadjelel[i];j++){
        fprintf (out,"  %ld",adjelel[i][j]);
      }
    }
    fprintf (out,"\n\n\n");
  }
  
  //  delete auxiliary array
  for (i=0;i<ne;i++){
    delete [] aux[i];
  }
  delete [] aux;
  
}



/**
   function assembles list of adjacent nodes to every node in the mesh
   
   @param out - output file
   
*/
void gtopology::adjacnodes (FILE *out)
{
  long i,j,k,l,m,min,prev,nne;
  long **auxadjnodnod;
  ivector nodes;

  //  allocation of array of counts of adjacent elements to nodes
  if (nadjelnod==NULL)
    nadjelnod = new long [nn];
  memset (nadjelnod,0,nn*sizeof(long));
  
  //  number of contributions
  for (i=0;i<ne;i++){
    nne = give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      nadjelnod[nodes[j]]++;
    }
  }
  
  //fprintf (out,"\n\n NADJELNOD");
  //for (i=0;i<nn;i++){
  //fprintf (out,"\n nadjelnod %5ld   %ld",i,nadjelnod[i]);
  //}
  
  
  //  allocation of array of numbers of adjacent elements to nodes
  if (adjelnod==NULL){
    adjelnod = new long* [nn];
    for (i=0;i<nn;i++){
      adjelnod[i] = new long [nadjelnod[i]];
    }
  }
  for (i=0;i<nn;i++){
    nadjelnod[i] = 0;
  }
  
  //  filling of array
  for (i=0;i<ne;i++){
    nne = give_nne (i);
    reallocv (RSTCKIVEC(nne,nodes));
    give_nodes (i,nodes);
    for (j=0;j<nne;j++){
      adjelnod[nodes[j]][nadjelnod[nodes[j]]++] = i;
    }
  }
  
  
  //fprintf (out,"\n\n ADJELNOD");
  //for (i=0;i<nn;i++){
  //  fprintf (out,"\n %ld    ",i);
  //  for (j=0;j<nadjelnod[i];j++){
  //    fprintf (out," %ld",adjelnod[i][j]);
  //  }
  //}
  
  
  //  array containing number of nodes adjacent to node before sorting
  nadjnodnod = new long [nn];
  for (i=0;i<nn;i++){
    nadjnodnod[i]=0;
    for (j=0;j<nadjelnod[i];j++){
      k=adjelnod[i][j];
      nadjnodnod[i]+=give_nne (k);
    }
  }
  
  //  array of adjacent nodes to each node before sorting
  auxadjnodnod = new long* [nn];
  for (i=0;i<nn;i++){
    auxadjnodnod[i] = new long [nadjnodnod[i]];
    nadjnodnod[i]=0;
  }
  
  //  filling of array
  for (i=0;i<nn;i++){
    for (j=0;j<nadjelnod[i];j++){
      l=adjelnod[i][j];
      nne=give_nne (l);
      reallocv (RSTCKIVEC(nne,nodes));
      give_nodes (l,nodes);
      for (k=0;k<nne;k++){
	auxadjnodnod[i][nadjnodnod[i]]=nodes[k];
	nadjnodnod[i]++;
      }
    }
  }
  
  //fprintf (out,"\n\n auxadjnodnod");
  //for (i=0;i<nn;i++){
  //  fprintf (out,"\n %ld    ",i);
  //  for (j=0;j<nadjnodnod[i];j++){
  //    fprintf (out," %ld",auxadjnodnod[i][j]);
  //  }
  // }
  
  
  //fprintf (out,"\n\n\n nadjnodnod\n");
  //for (i=0;i<nn;i++)
  //  fprintf(out,"%ld %ld\n",i,nadjnodnod[i]);
  
  //  sorting
  for (i=0;i<nn;i++){
    prev=LONG_MAX;
    for (j=0;j<nadjnodnod[i];j++){
      min=LONG_MAX;
      for (k=j;k<nadjnodnod[i];k++){
	if (min>auxadjnodnod[i][k]){
	  min=auxadjnodnod[i][k];  l=k;
	}
      }
      if (min==prev){
	nadjnodnod[i]--;
	auxadjnodnod[i][l]=auxadjnodnod[i][nadjnodnod[i]];
	j--;
      }
      else{
	m=auxadjnodnod[i][j];
	auxadjnodnod[i][j]=min;
	auxadjnodnod[i][l]=m;
	prev=min;
      }
    }
  }
  
  //fprintf (out,"\n\n\n nadjnodnod\n");
  //for (i=0;i<nn;i++)
  //  fprintf(out,"%ld %ld\n",i,nadjnodnod[i]);
  
  adjnodnod = new long* [nn];
  for (i=0;i<nn;i++){
    adjnodnod[i] = new long [nadjnodnod[i]];
    for (j=0;j<nadjnodnod[i];j++){
      adjnodnod[i][j]=auxadjnodnod[i][j];
    }
    delete [] auxadjnodnod[i];
  }
  delete [] auxadjnodnod;
  
  //  printing of array
  if (out)
  {
    fprintf (out,"\n\n\n LIST OF NUMBERS OF NODES ADJACENT TO NODE\n");
    for (i=0;i<nn;i++)
      fprintf(out,"%ld %ld\n",i,nadjnodnod[i]);
    fprintf (out,"\n\n\n LIST OF NODES ADJACENT TO NODE\n");
    for (i=0;i<nn;i++){
      fprintf (out,"\n node %4ld  ",i);
      for (j=0;j<nadjnodnod[i];j++){
        fprintf (out,"  %ld",adjnodnod[i][j]);
      }
    }
    fprintf (out,"\n\n\n");
  }
  
}


/**
   function creates list of adjacent nodes to nodes due to a common edge
   for each node, it creates a list of nodes connected by an edge
   
   JB
*/
void gtopology::adjacnodes_edge (FILE *out)
{
  long i,j,k,l,m;
  long ned,nned,auxnedges;
  long *auxnned,*dupledges,*aux;
  long **auxedgenodes;
  
   // number of all edges in mesh
  auxnedges = 0;
  for(i = 0; i < ne; i++){
    // number of edge on i-th element - gtopology
    ned = give_ned(i);
    auxnedges+=ned;
  }
  // //fprintf(out,"pocet hran %ld\n",auxnedges);

  // pomocne pole
  aux = new long[3];
  
  auxnned = new long[auxnedges];
  auxedgenodes = new long*[auxnedges];
  dupledges = new long[auxnedges];
  k = 0;
  for(i = 0; i < ne; i++){
    ned = give_ned(i);
    nned = give_nned(i);
    for(j = 0; j < ned; j++ ){
      give_edge_nodes(i,j,aux);
      auxnned[k] = nned;
      auxedgenodes[k] = new long[nned];
      for(l = 0; l < nned; l++){
	auxedgenodes[k][l] = aux[l];
      }
      // trideni pole s uzly na hrane(od nejmensiho po nejvetsi)
      if(nned == 2){
	if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	  l = auxedgenodes[k][1];
	  auxedgenodes[k][1] = auxedgenodes[k][0];
	  auxedgenodes[k][0] = l;
	}
      }
      if(nned == 3){
	if(auxedgenodes[k][0] > auxedgenodes[k][2]){
	  l = auxedgenodes[k][2];
	  auxedgenodes[k][2] = auxedgenodes[k][0];
	  auxedgenodes[k][0] = l;
	}
	if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	  l = auxedgenodes[k][1];
	  auxedgenodes[k][1] = auxedgenodes[k][0];
	  auxedgenodes[k][0] = l;
	}
	if(auxedgenodes[k][1] > auxedgenodes[k][2]){
	  l = auxedgenodes[k][2];
	  auxedgenodes[k][2] = auxedgenodes[k][1];
	  auxedgenodes[k][1] = l;
	}
      }
      dupledges[k] = 0;
      k++;
    }
  }
  delete []aux;
  //auxnedges = k;
  //fprintf(out,"pocet hran %ld\n",auxnedges);
  
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"hrana %ld ma uzly:",i);
  //for(j = 0; j < auxnned[i]; j++ ){
  //  fprintf(out,"   %ld",auxedgenodes[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  
  // hledani duplikovanych hran na hranici
  for(i = 0; i < auxnedges; i++){
    for(j = i+1; j < auxnedges; j++){
      if(auxedgenodes[i][0] == auxedgenodes[j][0]){
	if(auxnned[i] == 2 && auxnned[j] == 2){
	  if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	    dupledges[j]++;
	    //dupledges[i]++;
	  }
	}
	if(auxnned[i] == 3 && auxnned[i] == 3){
	  if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	    if(auxedgenodes[i][2] == auxedgenodes[j][2]){
	      dupledges[j]++;
	      //dupledges[i]++;
	    }
	  }
	}
      }
    }
  }
  
  //kontrolni tisk
  //   for(i = 0; i < auxnedges; i++){
  //     fprintf(out,"%ld dupledges %ld\n",i,dupledges[i]);
  //   }
  
    
  nadjacnodesedge = new long[nn];
  for(i = 0; i < nn; i++){
    nadjacnodesedge[i] = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	for(k = 0; k < auxnned[j]; k++){
	  if(auxedgenodes[j][k] == i){
	    auxedgenodes[j][k] = i;
	    nadjacnodesedge[i]++;
	  }
	}
      }
    } 
  }
  
  // kontrolni tisk
  // for(i = 0; i < auxnedges; i++){
  //   fprintf(out,"hrana %ld ma uzly:",i);
  //   for(j = 0; j < auxnned[i]; j++ ){
  //   fprintf(out,"   %ld",auxedgenodes[i][j]);
  //   }
  //   fprintf(out,"\n");
  //   }
  //   for(i = 0; i < nn; i++){
  //   fprintf(out,"%ld\n",nadjacnodesedge[i]);
  //   }  
  
  adjacnodesedge = new long*[nn];
  for(i = 0; i < nn; i++){
    adjacnodesedge[i] = new long[nadjacnodesedge[i]];
    m = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	if(auxnned[j] == 2){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == i){
	      if(k == 0){
		adjacnodesedge[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacnodesedge[i][m] = auxedgenodes[j][0];
		m++;
	      }
	    }
	  }
	}
	if(auxnned[j] == 3){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == i){
	      if(k == 0){
		adjacnodesedge[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacnodesedge[i][m] = auxedgenodes[j][0];
		m++;
		adjacnodesedge[i][m] = auxedgenodes[j][2];
		m++;
	      }
	      if(k == 2){
		adjacnodesedge[i][m] = auxedgenodes[j][1];
		m++;
	      }
	    }
	  }
	}
      }
    }
  }
  // //kontrolni tisk
  //   for(i = 0; i < nn; i++){
  //     fprintf(out,"uzel %ld ma  %ld sousedu: ",i,nadjacnodesedge[i]);
  //     for(j = 0; j < nadjacnodesedge[i]; j++){
  //       fprintf(out,"   %ld",adjacnodesedge[i][j]);
  //     }
  //     fprintf(out,"\n");
  //   }
  
  
  //sorting 
  for(i = 0; i < nn; i++){
    if(nadjacnodesedge[i] == 2){
      if(adjacnodesedge[i][0] > adjacnodesedge[i][1]){
	j = adjacnodesedge[i][1];
	adjacnodesedge[i][1] = adjacnodesedge[i][0];
	adjacnodesedge[i][0] = j;
      }
    }
    else{
      for(j = nadjacnodesedge[i] - 1; j > 0; j--){
	for(k = 0; k < j; k++){
	  if(adjacnodesedge[i][k] > adjacnodesedge[i][k+1]){
	    m = adjacnodesedge[i][k];   
	    adjacnodesedge[i][k] = adjacnodesedge[i][k+1];
	    adjacnodesedge[i][k+1] = m;
	  }
	} 
      }
    }
  }
  
  
  //kontrolni tisk
  for(i = 0; i < nn; i++){
    fprintf(out,"node %ld has  %ld neighbours: ",i,nadjacnodesedge[i]);
    for(j = 0; j < nadjacnodesedge[i]; j++){
      fprintf(out,"   %ld",adjacnodesedge[i][j]);
    }
    fprintf(out,"\n");
  }
  
  // mazani poli
  for(i = 0; i < auxnedges; i++){
    delete []auxedgenodes[i];
  }
  delete []auxnned;
  delete []auxedgenodes;
  delete []dupledges;
  
}



/**
  The function finds node which is the closest to the required point with the given coordinates

  @param x,y,z - coordinates of point
   
  Created by Tomas Koudelka 20.4.2018, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long gtopology::give_nearest_node(double x, double y, double z)
{
  long i, nid = 0;
  double d, min = DBL_MAX;
  double dx, dy, dz;
  
  for(i=0; i<nn; i++)
  {
    dx = gnodes[i].x - x;
    dy = gnodes[i].y - y;
    dz = gnodes[i].z - z;
    d = dx*dx + dy*dy + dz*dz;
    if (d < min)
    {
      min=d;
      nid=i;
    }
  }
  return(nid);
}



/**
   function backups code numbers for later changes in the supports during problem solving

   3.12.2002
*/
void gtopology::backup_codnum (void)
{
  long i, j, ndofn;

  for (i = 0; i < nn; i++)
  {
    ndofn = give_ndofn (i);
    if (ndofn < 0) // due to hanging nodes
      ndofn = -ndofn;

    bckcn[i] = new long[ndofn];
    for (j = 0; j < ndofn; j++)
      bckcn[i][j]=give_dof (i,j);
  }
}



/**
   function backups code numbers for later changes in the supports during problem solving

   3.12.2002
*/
void gtopology::restore_codnum (void)
{
  long i, j, ndofn;

  for (i = 0; i < nn; i++)
  {
    ndofn = give_ndofn (i);
    if (ndofn < 0) // due to hanging nodes
      ndofn = -ndofn;
    for (j = 0; j < ndofn; j++)
      gnodes[i].cn[j] = bckcn[i][j];
  }
}

void gtopology::write_nodes(FILE *out)
{
  fprintf(out, "NODES\n");
  for (long i = 0; i < nn; i++)
    fprintf(out, "%ld %e %e %e\n", i, gnodes[i].x, gnodes[i].y, gnodes[i].z);
}

/**
   function assembles correspondation between usual nodes and layered nodes
   function creates array unln with nn components, unln[i] represents
   number of layered node which contains i-th usual node
   
   JK, 1.2.2003
*/
void gtopology::unodelnode ()
{
  long i,j;
  unln = new long [nn];
  unnl = new long [nn];
  for (i=0;i<nln;i++){
    for (j=0;j<lgnodes[i].nl;j++){
      unln[lgnodes[i].nodes[j]]=i;
      unnl[lgnodes[i].nodes[j]]=j;
    }
  }
}

/**
   function initiates number of multipliers on elements
   
   JK, 2.2.2003
*/
void gtopology::initiate_elemmult ()
{
  long i,nid,lnid;

  for (i=0;i<ne;i++){
    nid=gelements[i].nodes[0];
    lnid=unln[nid];
    gelements[i].nmult=give_nmult (lnid);
  }
}

/**
   function determines maximum sizes of solved domain
   
   6.10.2003, JK
*/
void gtopology::comp_domain_sizes ()
{
  long i;
  double x,y,z,minx,miny,minz,maxx,maxy,maxz;
  
  domsizes = new double [3];
  
  minx=gnodes[0].x;  maxx=gnodes[0].x;
  miny=gnodes[0].y;  maxy=gnodes[0].y;
  minz=gnodes[0].z;  maxz=gnodes[0].z;
  
  for (i=1;i<nn;i++){
    x=gnodes[i].x;
    y=gnodes[i].y;
    z=gnodes[i].z;
    
    if (minx>x)  minx=x;
    if (maxx<x)  maxx=x;
    
    if (miny>y)  miny=y;
    if (maxy<y)  maxy=y;
    
    if (minz>z)  minx=z;
    if (maxz<z)  maxx=z;
  }
  
  domsizes[0]=maxx-minx;
  domsizes[1]=maxy-miny;
  domsizes[2]=maxz-minz;
}

/**
   function returns sizes of domain
   
   6.10.2003, JK
*/
void gtopology::give_domain_sizes (double *sizes)
{
  sizes[0]=domsizes[0];
  sizes[1]=domsizes[1];
  sizes[2]=domsizes[2];
}



/**
  The function returns the square of maximum distance between given point pt and
  nodes of element eid.

  @param eid - element id
  @param pt  - %vector of point coordinates in 3D (x,y,z)

  @return The function returns the square of maximum distance from nodes.

  Created by TKo, 1.12.2016
*/
double gtopology::max_sqrdist_nod_pt(long eid, vector &pt)
{
  long i;
  long nne = give_nne(eid);
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  double d, maxd=0.0;

  give_node_coord3d(x, y, z, eid);
  
  for(i=0; i<nne; i++)
  {
    d  = (x[i] - pt[0])*(x[i] - pt[0]);
    d += (y[i] - pt[1])*(y[i] - pt[1]);
    d += (z[i] - pt[2])*(z[i] - pt[2]);
    if (d > maxd)
      maxd = d;
  }
  return maxd;
}



/**
   function reads integer table functions
   function allocates array leso
   
   @param in - input stream
   
   10.2.2006, JK, TKo
*/
void gtopology::read_gf (XFILE *in)
{
  long i,j;
  xfscanf (in,"%k%ld","num_gfunct", &ngf);
  gf = new gfunct [ngf];
  
  for (i=0;i<ngf;i++){
    xfscanf (in,"%k%ld", "gf_id",&j);
    if ((j < 1) || (j > ngf))
    {
      print_err("general function index %ld out of range <1,%ld>\n"
                " at line=%ld, col=%ld, file %s", 
                __FILE__, __LINE__, __func__, j, ngf, in->line, in->col, in->fname);
      abort();
    }
    j--;
    gf[j].read (in);
  }
}


/**
   function prints integer table functions
   
   @param out - output stream
   
   17.7.2006, TKr
*/
void gtopology::print_gf (FILE *out)
{
  long i;
  fprintf (out,"\n%ld\n\n",ngf);
  
  for (i=0;i<ngf;i++){
    fprintf (out,"%ld\n",i+1);
    gf[i].print (out);
  }
}

/**
   function allocates and assigns array of nodes and elements which are switched on
   
   3.3.2006, JK
*/
void gtopology::lneso_init ()
{
  long i;
  
  //  array of elements which are switched on
  leso = new long [ne];
  
  for (i=0;i<ne;i++){
    leso[i]=1;
  }

  //  array of nodes which are switched on
  lnso = new long [nn];
  
  for (i=0;i<nn;i++){
    lnso[i]=1;
  }
}


/**
   function initiates variables auxinf on elements
   
   JK, 3.3.2006
*/
void gtopology::auxinf_init ()
{
  long i;
  
  for (i=0;i<ne;i++){
    gelements[i].auxinf=0;
  }
}



/**
   The function updates variables auxinf on elements to the corresponding leso values.
   It is used in problems with changing number of nodes, elements and DOFs.
   
   @return The function does not return anything but it changes values of auxinf on elements.
   
   Created by Tomas Koudelka, 7.6. 2013
*/
void gtopology::update_auxinf()
{
  long i;
  
  for (i=0;i<ne;i++)
    gelements[i].auxinf = leso[i];
}



/**
   The function switches on appropriate elements.
   It is used in problems with changing number of nodes, elements and DOFs.
   
   leso[i]=0 - the i-th element is switched off
   leso[i]=1 - the i-th element is switched on
   
   @param time - actual time

   @return The function does not return anything but it changes values stored in the leso array.
   
   JK, 28.2.2006, last modification 3.11.2006 
*/
void gtopology::update_elem (double time)
{
  long i,j,k;
  
  for (i=0;i<ne;i++){
    //  id of time function
    j=gelements[i].tgf;
    //  value of time function
    k=gf[j].getval_long (time);

    leso[i]=k;
  }
}

/**
   function switches on appropriate nodes
   function is used in problems with changing number of nodes, elements and DOFs
   former nodes are not temporarily switched off
   
   lnso[i]=0 - the i-th node is switched off
   lnso[i]=1 - the i-th node is switched on
   
   JK, 7.11.2006
*/
void gtopology::update_nodes ()
{
  long i,j,nne, nmne;
  ivector nod;
  
  for (i=0;i<nn;i++){
    lnso[i]=0;
  }
  
  for (i=0;i<ne;i++){
    if (leso[i]==1){
      //  only element switched on are investigated
      
      //  number of nodes on element (or matser nodes in case of a hanging node)
      nne = give_nne (i);
      reallocv (RSTCKIVEC(nne,nod));
      //  element nodes
      give_nodes (i,nod);
      for (j=0;j<nne;j++){
	lnso[nod[j]]=1;
      }
      nmne = give_nmne(i);
      if (nmne > 0)  // switch on also hanging nodes
      {
        //  number of original nodes of element regardless on hanging nodes
        nne = give_original_nne(i);
        reallocv(RSTCKIVEC(nne, nod));
        give_original_nodes(i,nod);
        for (j=0;j<nne;j++)
        {
          if (gnodes[nod[j]].mnodes)
            lnso[nod[j]]=1;
        }
      }
    }
  }
}



/**
   The function switches off selected nodes and it is used in problems with 
   changing number of nodes, elements and DOFs. The selected nodes are given by array 
   sn where sn[i]=1 means switch off i-th node. If the status of i-th node should be 
   unchanged then sn[i]=0.
   
   @param sn - array of with indicators of selected nodes   

   @return The function does not return anything but it changes nodal status, i.e., array lnso.

   Created by Tomas Koudelka, 10.6.2013
*/
void gtopology::remove_nodes(long *sn)
{
  long i;
  ivector nod;
  
  for (i=0;i<nn;i++)
  {
    if (sn[i])
      lnso[i]=0;
  }  
}



/**
   Function searches for elements elements with changed status, i.e.,
   elements that were switched on/off in the actual time step. 
   Additionally, it counts number of elements added. Status of all elements
   is updated. The function is used in problems with changing number of nodes, 
   elements and DOFs.
   
   @param time - actual time
   @param nae  - number of added elements (output parameter)

   @return The function returns elements with changed status and number of added elements via parameter nae .
   
   Created by TKo, 7.6.2013
*/
long gtopology::search_changed_elem(double time, long &nae)
{
  long i, j, k, l;
  long nce;

  //  number of added elements at this time
  nae=0;
  // number of changed elements
  nce = 0;

  for (i=0;i<ne;i++)
  {
    //  id of time function
    j=gelements[i].tgf;
    //  value of time function at actual time
    k=gf[j].getval_long(time);
    //  value of time function at previous time
    l=gelements[i].auxinf;
    // actualize status of all elements
    leso[i] = k;
    
    if (k==l)  // nothing changed on element
      continue;
    else
    { 
      if (k) // element was added
        nae++;
      nce++;    // element was either added or removed
    }
  }
  return nce;
}



/**
  The function creates list of nodes that belong to interfaces
  between old and new elements (added at actual time step).
  It sets indicators in array ifn for each node on interface to 1. 
  Remaining nodes are set to 0.

  @param ifn - array of indicators for interface nodes
  
  @return The function returns number of interface nodes and
          indicators in the array ifn are set.
 
  Created by Tomas Koudelka, 7.6.2013
*/
long gtopology::search_iface_nodes(long *ifn)
{
  long i, j, nifn, nne;
  ivector nod;

  memset(ifn, 0, sizeof(*ifn)*nn);

  // mark nodes that belong to old element by 1
  for (i=0;i<ne;i++)
  {
    if ((leso[i]) && (gelements[i].auxinf)) // for old elements
    {
      nne = give_nne(i);
      reallocv(RSTCKIVEC(nne,nod));
      //  element nodes
      give_nodes(i,nod);
      
      for(j=0; j<nne; j++)
      {
        if (ifn[nod[j]] == 0)
          ifn[nod[j]] = 1;
      }
    }
  }  
  // mark nodes that belong to old and new elements by 2
  for (i=0;i<ne;i++)
  {
    if ((leso[i]) && (gelements[i].auxinf == 0)) // for new elements
    {    
      nne = give_nne(i);
      reallocv(RSTCKIVEC(nne,nod));
      //  element nodes
      give_nodes(i,nod);
      for(j=0; j<nne; j++)
      {
        if (ifn[nod[j]] == 1) // if the node belongs to old element
          ifn[nod[j]] = 2;
      }
    }
  }
  nifn = 0;
  for (i=0; i<nn; i++)
  {
    if (ifn[i] < 2)
      ifn[i] = 0;
    else
    {
      ifn[i] = 1;
      nifn++;
    }
  }
  return nifn;
}



/**
   function searches for new elements
   function is used in problems with changing number of nodes, elements and DOFs
   function returns number of changed elements
   only new elements are switched on, former elements are temporarily switched off
   
   @param time - actual time
   @param prev_time - previous time
   
   JK, 22.9.2006, last modification 3.11.2006
*/
long gtopology::search_newelem (double time,double prev_time)
{
  long i,j,k,l,nae;
  
  //  number of added elements at this time
  nae=0;
  
  if (time>prev_time){
    for (i=0;i<ne;i++){
      //  id of time function
      j=gelements[i].tgf;
      //  value of time function at actual time
      k=gf[j].getval_long (time);
      //  value of time function at previous time
      l=gf[j].getval_long (prev_time);
      
      if (k==l){
	leso[i]=0;
      }
      else{
	leso[i]=1;
	nae++;
      }
    }
  }
  else{
    //  initial time step
    //  time = starting_time
    for (i=0;i<ne;i++){
      //  id of time function
      j=gelements[i].tgf;
      //  value of time function at actual time
      k=gf[j].getval_long (time);
      
      if (k==0){
	leso[i]=0;
      }
      else{
	leso[i]=1;
	nae++;
      }
    }
  }
  
  return nae;
}

/**
   function searches for new DOFs
   function is used in problems with changing number of nodes, elements and DOFs
   only new DOFs are switched on, former DOFs are temporarily switched off
   
   @param time - actual time
   @param prev_time - previous time
   
   JK, 22.9.2006, modification 3.11.2006
   Updated 21.5.2013,TKo
*/
long gtopology::search_newdofs (double time,double prev_time)
{
  long i,nch;
  long *plnso = new long[nn];
  // backup of previous state of nodes
  memcpy(plnso, lnso, sizeof(*plnso)*nn);
  
  // update state of nodes according to the state of elements
  // update_elements() must be called before search_newdofs()
  update_nodes();

  //  number of changes
  nch=0;
  for (i=0;i<nn;i++){
    nch += gnodes[i].search_changed_dofs (gf,time,prev_time,lnso[i],plnso[i]);
  }
  
  //  state of code numbers is changed
  if (nch > 0)
    cnstate=0;

  delete [] plnso;
  return nch;
}



/**
   The function searches for changed DOFs and it is used in problems with 
   changing number of nodes, elements and DOFs.
   
   @param time - actual time
   @param prev_time - previous time
   
   JK, 22.9.2006, modification 3.11.2006
   Rewritten by 21.5.2013,TKo
*/
long gtopology::search_changed_dofs (double time,double prev_time)
{
  long i,nch;
  long *plnso = new long[nn];
  // backup of previous state of nodes
  memcpy(plnso, lnso, sizeof(*plnso)*nn);
  
  // update state of nodes according to the state of elements
  // update_elements() must be called before search_newdofs()
  update_nodes();

  //  number of changes
  nch=0;
  for (i=0;i<nn;i++){
    nch += gnodes[i].search_changed_dofs(gf,time,prev_time,lnso[i],plnso[i]);
  }
  
  //  state of code numbers is changed
  if (nch > 0)
    cnstate=0;

  delete [] plnso;
  return nch;
}



/**
  The function switch on elements that were added in the actual time step. 
  Remaining elements are switched off.

  @return The function does not return anything but it changes content of leso array.

  Created by Tomas Koudelka, 17.6.2013
*/
void gtopology::switch_new_elem()
{
  long i;
  for (i=0; i<ne; i++)
  {
    if (leso[i] && (gelements[i].auxinf == 0)) // leave the status of new elements
      continue;
    else
      leso[i] = 0;  // switch off old elements
  }
}



/**
  The function switch on elements that were added in the actual time step. 
  Remaining elements are switched off.

  @return The function does not return anything but it changes content of leso array.

  Created by Tomas Koudelka, 11.5.2016
*/
void gtopology::switch_removed_elem()
{
  long i;
  for (i=0; i<ne; i++)
  {
    if ((leso[i]==0) && gelements[i].auxinf) // leave the status of removed elements
      leso[i] = 1;
    else
      leso[i] = 0;  // switch off remaining elements
  }
}



/**
   The function updates DOFs on active nodes and resets cne state flag on elements with hanging nodes.
   It is used in problems with changing number of nodes, elements and DOFs.
   
   @param time - actual time

   @return The funtion does not return anything but it changes code number array at nodes.
   
   Created by Tomas Koudelka, 7.6.2013
*/
void gtopology::update_active_dofs (double time)
{
  long i;
  
  for (i=0;i<nn;i++){
    if (lnso[i])
      gnodes[i].update_dofs (gf,time,lnso[i]);
    else
    {
      if (gnodes[i].ndofn > 0) // due to hanging nodes
        gnodes[i].clear_dof();
    }
  }
  
  for (i=0;i<ne;i++)
  {
    if (gelements[i].master_nodes!=NULL)
    {
      //  this is an element with hanging nodes
      gelements[i].cne=0;
    }
  }

  //  state of code numbers is changed
  cnstate=0;
}



/**
   The function switch off DOFs on interface nodes
   function is used in problems with changing number of nodes, elements and DOFs
   former DOfs are not temporarily switched off
   
   @param ifn - array of indicators whether i-th node is on the interface (ifn[i] = 1) or not (ifn[i]=0)
   
   Created by Tomas Koudelka, 06.2013
*/
void gtopology::clear_intf_dofs(long *ifn)
{
  long i;
  
  for (i=0;i<nn;i++)
  {
    if (ifn[i])
      gnodes[i].clear_dof();
  }

  //  state of code numbers is changed
  cnstate=0;
}



/**
   The function updates appropriate DOFs and it is 
   used in problems with changing number of nodes, elements and DOFs.
   former DOfs are not temporarily switched off
   
   @param time - actual time
   
   6.2.2006, JK, last modification 3.11.2006
*/
void gtopology::update_dofs (double time)
{
  long i;
  
  for (i=0;i<nn;i++)
    gnodes[i].update_dofs (gf,time,lnso[i]);

  //  state of code numbers is changed
  cnstate=0;
}


/**
   function changes nodes and elements
   it is used in problems with growing numbers of nodes and elements
   nodes and elements are updated at particular instance of time
   
   JK, 28.2.2006
*/
/*
void gtopology::update_mesh (double time,double prev_time,long &ncn,long &nce)
{

  //  update of nodes
  //  function returns number of changed DOFs
  ncn = update_dofs (time,prev_time);
  
  //  update of element status
  nce = update_elemstat (time);

}
*/

void gtopology::print_time_functions (FILE *kon,double time)
{
  long i;
  //double time1,time2,time=10800000.0,dt=1123200.000000;
  //FILE *kon;
  //kon = fopen ("casfun","w");
  
  //time1=time;
  //time2=time+dt;
  fprintf (kon,"\n\n vypis v case %e\n",time);
  for (i=0;i<ngf;i++){
    //fprintf (kon,"%ld %ld\n",gf[i].getval_long (time1),gf[i].getval_long (time2));
    fprintf (kon,"%ld %ld\n",i,gf[i].getval_long (time));
  }
  
  //fclose (kon);
}







/**
   function searches for common end nodes on interfaces
   it is used in automatic subdomain detection

   function is used for problems with hemivariational inequalities,
   material discontinuities,
   
   JK, 12.10.2008
*/
void gtopology::end_nodes (FILE *out)
{
  long i,j,k,nsn,nene;
  long *nodes,*sn;
  long *av;
  
  nodes = new long[2];
  
  //  auxiliary array
  av = new long [nn];
  for (i=0;i<nn;i++){
    av[i]=0;
  }
  
  //  loop over elements
  for (i=0;i<ne;i++){
    //  number of end nodes on element
    nene=give_nen (i);
    
    if (nene>0){
      //  nodes on required edge
      give_end_nodes (i,nodes);
      
      //  loop over egdes on element
      for (j=0;j<2;j++){
	
	//  positions in array av are filled with respect to end nodes
	//  internal end points are characterized by av[i]=2
	//  boundary or interface end points  are characterized by av[i]=1
	av[nodes[j]]++;
      }
    }
  }
  

  fprintf (out,"\n\n\n kontrola koncovych bodu \n\n");
  for (i=0;i<nn;i++){
    fprintf (out,"\n node %4ld  av %4ld",i,av[i]);
  }
  fprintf (out,"\n\n\n");
  
  //  number of suspicious nodes
  nsn=0;
  for (i=0;i<nn;i++){
    if (av[i]==1)
      nsn++;
  }
  
  //  list of suspicious nodes
  sn = new long [nsn];
  nsn=0;
  for (i=0;i<nn;i++){
    if (av[i]==1){
      sn[nsn]=i;
      nsn++;
    }
  }

  //  number of end nodes
  nen=0;
  double threshold=1.0e-7;
  for (i=0;i<nsn;i++){
    for (j=i+1;j<nsn;j++){
      if (fabs(gnodes[sn[j]].x-gnodes[sn[i]].x)<threshold){
	if (fabs(gnodes[sn[j]].y-gnodes[sn[i]].y)<threshold){
	  if (fabs(gnodes[sn[j]].z-gnodes[sn[i]].z)<threshold){
	    nen++;
	  }
	}
      }
    }
  }
  
  fprintf (stdout,"\n pocet zajimavych koncovych bodu  %ld",nen);
  fprintf (out,"\n\n\n pocet zajimavych koncovych bodu  %ld \n\n",nen);
  
  endnodes = new endnode [nen];
  k=0;
  for (i=0;i<nsn;i++){
    for (j=i+1;j<nsn;j++){
      if (fabs(gnodes[sn[j]].x-gnodes[sn[i]].x)<threshold){
	if (fabs(gnodes[sn[j]].y-gnodes[sn[i]].y)<threshold){
	  if (fabs(gnodes[sn[j]].z-gnodes[sn[i]].z)<threshold){
	    endnodes[k].nn=2;
	    endnodes[k].nm=2;
	    endnodes[k].fn=sn[i];
	    endnodes[k].ln=sn[j];
	    k++;
	  }
	}
      }
    }
  }
  
  delete [] nodes;
  delete [] sn;
  delete [] av;
}

/**
   function prints auxiliary data about end nodes
   
   @param out - output file
   
   JK, 13.10.2008
*/
void gtopology::endnodes_auxprint (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n number of end nodes %ld",nen);
  
  for (i=0;i<nen;i++){
    fprintf (out,"\n\n end node number %6ld",i);
    fprintf (out,"\n number of nodes  %ld",endnodes[i].nn);
    fprintf (out,"\n number of Lagrange multipliers  %ld",endnodes[i].nm);
    fprintf (out,"\n node number 1    %ld",endnodes[i].fn);
    fprintf (out,"\n node number 2    %ld",endnodes[i].ln);
  }
}

/**
   function allocates array cnm on end nodes
   cnm contains code numbers of Lagrange multipliers
   
   JK, 13.10.2008
*/
void gtopology::alloc_endnode_cn ()
{
  long i,nid,ndofn;
  
  for (i=0;i<nen;i++){
    //  first node
    nid=endnodes[i].fn;
    //  number of DOFs on the first node
    ndofn=give_ndofn (nid);
    
    //  allocation of array cnm
    endnodes[i].alloc_cnm (ndofn);
  }

}

/**
   function localizes contributions from end nodes to the global %matrix
   
   @param gm - pointer to the general %matrix
   
   JK, 14.10.2008
*/
void gtopology::endnodes_localization (gmatrix *gm)
{
  long i,nm,*ncn1,*ncn2,*mcn;

  for (i=0;i<nen;i++){
    //  number of multipliers between nodes
    nm=endnodes[i].ndofn;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];

    //  code numbers on first nodes and appropriate multipliers
    give_endnode_code_numbers (i,ncn1,ncn2,mcn);
    
    //  localization of contributions from edges
    gm->mult_localize (nm,ncn1,ncn2,mcn);
    
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }
}

/**
   function assignes numbers of all elements to end nodes
   
   JK, 17.10.2008
*/
void gtopology::enodes_allelem (FILE */*out*/)
{
  long i,k,l,fn,nne,stop;
  ivector nod;
  
  //  loop over all end nodes
  for (i=0;i<nen;i++){
    //  allocation of array adjel
    endnodes[i].adjel = new long [2];

    //  first node
    fn=endnodes[i].fn;
    
    //  element containing node number fn has to be searched for
    //  it is necessary to go through all elements until appropriate element is found
    stop=0;
    for (k=0;k<ne;k++){
      //  number of nodes on element
      nne = give_original_nne (k);
      //  array allocation
      reallocv (RSTCKIVEC(nne,nod));
      //  element nodes
      give_original_nodes (k,nod);
      
      //  loop over nodes on element
      for (l=0;l<nne;l++){
	if (nod[l]==fn){
	  endnodes[i].adjel[0]=k;
	  stop=1;
	  break;
	}
      }
      if (stop==1)
	break;
    }


    //  second node
    fn=endnodes[i].ln;
    
    //  element containing node number fn has to be searched for
    //  it is necessary to go through all elements until appropriate element is found
    stop=0;
    for (k=0;k<ne;k++){
      //  number of nodes on element
      nne = give_nne (k);
      //  array allocation
      reallocv (RSTCKIVEC(nne,nod));
      //  element nodes
      give_original_nodes (k,nod);
      
      //  loop over nodes on element
      for (l=0;l<nne;l++){
	if (nod[l]==fn){
	  endnodes[i].adjel[1]=k;
	  stop=1;
	  break;
	}
      }
      if (stop==1)
	break;
    }
  }
}








































/**
   function searches for edges between two different materials
   it is used in automatic subdomain detection
   
   the instances edge of the class %gedge are found with
   respect to material models used along the edges
   
   function is used for problems with radiation in TRFEL
   
   @param matrad - type of material model connected with the radiation
   @param out - output file

   JK, 25.5.2011
*/
void gtopology::edges_mat (long matrad,FILE */*out*/)
{
  long i,j,k,l,m,elid,nsed,ned,nedk,nned,nnedk,n1,n2,n3,n4,aux,mat1,mat2,stopp;
  long *nodes,*nodesk;
  
  //  number of suspicious edges
  nsed=0;
  
  //  loop over elements
  for (i=0;i<ne;i++){
    //  material id
    mat1=gelements[i].auxinf;
    //  number of edges on element
    ned=give_ned (i);
    //  number of nodes on one edge
    nned=give_nned (i);
    //  nodes on required edge
    nodes=new long [nned];
    
    //  loop over egdes on element
    for (j=0;j<ned;j++){
      //  nodes on required edge
      give_edge_nodes (i,j,nodes);
      //  end nodes of the edge
      n1=nodes[0];
      n2=nodes[nned-1];
      if (n1>n2){
	aux=n2;
	n2=n1;
	n1=aux;
      }
	
      stopp=0;
      
      //  loop over the adjacent elements
      for (k=0;k<nadjelel[i];k++){
	//  element id
	elid=adjelel[i][k];
	
	//  material id
	mat2=gelements[elid].auxinf;
	//  number of edges on element
	nedk=give_ned (elid);
	//  number of nodes on one edge
	nnedk=give_nned (elid);
	//  nodes on required edge
	nodesk=new long [nnedk];
	
	//  loop over egdes on the adjacent element
	for (l=0;l<nedk;l++){
	  //  nodes on required edge
	  give_edge_nodes (elid,l,nodesk);
	  //  end nodes of the edge
	  n3=nodesk[0];
	  n4=nodesk[nnedk-1];
	  if (n3>n4){
	    aux=n4;
	    n4=n3;
	    n3=aux;
	  }
	  
	  if (n1==n3 && n2==n4){
	    if ((mat1==matrad && mat2!=matrad) || (mat1!=matrad && mat2==matrad)){
	      nsed++;
	      //gelements[i].auxinf=-1;
	      //gelements[elid].auxinf=-1;
	      stopp=1;
	    }
	  }
	  if (stopp==1){
	    break;
	  }
	}
	delete [] nodesk;
	if (stopp==1){
	  break;
	}
      }
    }
    delete [] nodes;
  }
  
  //  each edge has been taken into account two times, therefore, the number
  //  of edges is the following
  nged=nsed/2;
  
  //  general edges
  gedges = new gedge [nged];
  
  
  nsed=0;
  //  loop over elements
  for (i=0;i<ne;i++){
    //  material id
    mat1=gelements[i].auxinf;
    //  number of edges on element
    ned=give_ned (i);
    //  number of nodes on one edge
    nned=give_nned (i);
    //  nodes on required edge
    nodes=new long [nned];
    
    //  loop over egdes on element
    for (j=0;j<ned;j++){
      //  nodes on required edge
      give_edge_nodes (i,j,nodes);
      //  end nodes of the edge
      n1=nodes[0];
      n2=nodes[nned-1];
      if (n1>n2){
	//  n1 contains smaller node number
	//  n2 contains greater node number
	aux=n2;
	n2=n1;
	n1=aux;
      }
	
      stopp=0;
      
      //  loop over the adjacent elements
      for (k=0;k<nadjelel[i];k++){
	//  element id
	elid=adjelel[i][k];
	
	//  material id
	mat2=gelements[elid].auxinf;
	//  number of edges on element
	nedk=give_ned (elid);
	//  number of nodes on one edge
	nnedk=give_nned (elid);
	//  nodes on required edge
	nodesk=new long [nnedk];
	
	//  loop over egdes on the adjacent element
	for (l=0;l<nedk;l++){
	  //  nodes on required edge
	  give_edge_nodes (elid,l,nodesk);
	  //  end nodes of the edge
	  n3=nodesk[0];
	  n4=nodesk[nnedk-1];
	  if (n3>n4){
	    //  n3 contains smaller node number
	    //  n4 contains greater node number
	    aux=n4;
	    n4=n3;
	    n3=aux;
	  }
	  
	  if (n1==n3 && n2==n4){
	    //  nodes on element edges are identical
	    if ((mat1==matrad && mat2!=matrad) || (mat1!=matrad && mat2==matrad)){
	      
	      //  appropriate edge is searched or established
	      for (m=0;m<nsed;m++){
		if (gedges[m].fn==n1 && gedges[m].ln==n2){
		  //  the appropriate edges has been established
		  stopp=1;
		  break;
		}
	      }
	      if (stopp==0){
		//  the appropriate edge has not been established
		//  it is established now
		gedges[nsed].fn=n1;
		gedges[nsed].ln=n2;
		gedges[nsed].nn=4;
		gedges[nsed].nm=2;
		
		gedges[nsed].nlist = new long [4];
		gedges[nsed].nlist[0]=n1;
		gedges[nsed].nlist[1]=n1;
		gedges[nsed].nlist[2]=n2;
		gedges[nsed].nlist[3]=n2;
		
		gedges[nsed].adjel = new long [2];
		gedges[nsed].adjel[0]=i;
		gedges[nsed].adjel[1]=elid;

		nsed++;
		stopp=1;
	      }
	    }
	  }
	  if (stopp==1){
	    break;
	  }
	}//  end of loop over the edges on the adjacent element
	delete [] nodesk;
	if (stopp==1){
	  break;
	}
      }//  end of loop over the adjacent elements
    }//  end of loop over the edges on the element
    delete [] nodes;
  }//  end of loop over the elements

}


/**
   function searches for common edges on interfaces
   it is used in automatic subdomain detection

   nodes and edges on the interfaces are doubled,
   the instances edge of the class %gedge are found with respect to the
   identical node coordinates
   
   function is used for problems with hemivariational inequalities
   
   JK, 11.7.2007
*/
void gtopology::edges (FILE *out)
{
  long i,j,k,ned,nned,nsed;
  long *nodes,*fn,*ln;
  long **av;

  //  auxiliary array
  //  it has the same structure as the array adjnodnod
  av = new long* [nn];
  for (i=0;i<nn;i++){
    av[i] = new long [nadjnodnod[i]];
    for (j=0;j<nadjnodnod[i];j++){
      av[i][j]=0;
    }
  }
  
  //  loop over elements
  for (i=0;i<ne;i++){
    //  number of edges on element
    ned=give_ned (i);
    //  number of nodes on one edge
    nned=give_nned (i);
    //  nodes on required edge
    nodes=new long [nned];
    
    //  loop over egdes on element
    for (j=0;j<ned;j++){
      //  nodes on required edge
      give_edge_nodes (i,j,nodes);
      
      //  positions in array av are filled with respect to edge nodes
      //  internal edges (between the i-th and the j-th node) are characterized by av[i][j]=2
      //  boundary or interface edges (between the i-th and the j-th node) are characterized by av[i][j]=1
      for (k=0;k<nadjnodnod[nodes[0]];k++){
	if (adjnodnod[nodes[0]][k]==nodes[nned-1]){
	  av[nodes[0]][k]++;
	  break;
	}
      }
      for (k=0;k<nadjnodnod[nodes[nned-1]];k++){
	if (adjnodnod[nodes[nned-1]][k]==nodes[0]){
	  av[nodes[nned-1]][k]++;
	  break;
	}
      }
    }

    delete [] nodes;
  }
  

  fprintf (out,"\n\n\n kontrola hran \n\n");
  for (i=0;i<nn;i++){
    fprintf (out,"\n node %4ld  ",i);
    for (j=0;j<nadjnodnod[i];j++){
      fprintf (out,"    %4ld : %2ld",adjnodnod[i][j],av[i][j]);
    }
  }
  fprintf (out,"\n\n\n");

  
  //  number of suspicious edges
  nsed=0;
  for (i=0;i<nn;i++){
    for (j=0;j<nadjnodnod[i];j++){
      if (av[i][j]==1)
	nsed++;
    }
  }
  nsed/=2;
  
  //  list of first nodes on edges
  fn = new long [nsed];
  //  list of last nodes on edges
  ln = new long [nsed];
  
  k=0;
  for (i=0;i<nn;i++){
    for (j=0;j<nadjnodnod[i];j++){
      if (av[i][j]==1){
	if (i<adjnodnod[i][j]){
	  fn[k]=i;
	  ln[k]=adjnodnod[i][j];
	  k++;
	}
      }
    }
  }
  
  fprintf (out,"\n\n\n druha kontrola hran \n\n");
  for (i=0;i<nsed;i++){
    fprintf (out,"\n fn  %5ld  ln %5ld",fn[i]+1,ln[i]+1);
  }
  
  //fclose (out);
  //abort ();

  //  number of generalized edges
  nged=0;
  double threshold=1.0e-7;
  for (i=0;i<nsed;i++){
    for (j=i+1;j<nsed;j++){
      if (fabs(gnodes[fn[j]].x-gnodes[fn[i]].x)<threshold && fabs(gnodes[ln[j]].x-gnodes[ln[i]].x)<threshold){
	if (fabs(gnodes[fn[j]].y-gnodes[fn[i]].y)<threshold && fabs(gnodes[ln[j]].y-gnodes[ln[i]].y)<threshold){
	  if (fabs(gnodes[fn[j]].z-gnodes[fn[i]].z)<threshold && fabs(gnodes[ln[j]].z-gnodes[ln[i]].z)<threshold){
	    nged++;
	  }
	}
      }
      if (fabs(gnodes[fn[j]].x-gnodes[ln[i]].x)<threshold && fabs(gnodes[ln[j]].x-gnodes[fn[i]].x)<threshold){
	if (fabs(gnodes[fn[j]].y-gnodes[ln[i]].y)<threshold && fabs(gnodes[ln[j]].y-gnodes[fn[i]].y)<threshold){
	  if (fabs(gnodes[fn[j]].z-gnodes[ln[i]].z)<threshold && fabs(gnodes[ln[j]].z-gnodes[fn[i]].z)<threshold){
	    nged++;
	  }
	}
      }
    }
  }
  
  fprintf (out,"\n\n\n pocet zajimavych hran  %ld \n\n",nged);
  
  gedges = new gedge [nged];
  k=0;
  for (i=0;i<nsed;i++){
    for (j=i+1;j<nsed;j++){
      if (fabs(gnodes[fn[j]].x-gnodes[fn[i]].x)<threshold && fabs(gnodes[ln[j]].x-gnodes[ln[i]].x)<threshold){
	if (fabs(gnodes[fn[j]].y-gnodes[fn[i]].y)<threshold && fabs(gnodes[ln[j]].y-gnodes[ln[i]].y)<threshold){
	  if (fabs(gnodes[fn[j]].z-gnodes[fn[i]].z)<threshold && fabs(gnodes[ln[j]].z-gnodes[ln[i]].z)<threshold){
	    gedges[k].nn=4;
	    gedges[k].nm=2;
	    gedges[k].fn=fn[i];
	    gedges[k].ln=ln[i];
	    gedges[k].nlist = new long [4];
	    gedges[k].nlist[0]=fn[i];
	    gedges[k].nlist[1]=fn[j];
	    gedges[k].nlist[2]=ln[i];
	    gedges[k].nlist[3]=ln[j];
	    k++;
	  }
	}
      }
      if (fabs(gnodes[fn[j]].x-gnodes[ln[i]].x)<threshold && fabs(gnodes[ln[j]].x-gnodes[fn[i]].x)<threshold){
	if (fabs(gnodes[fn[j]].y-gnodes[ln[i]].y)<threshold && fabs(gnodes[ln[j]].y-gnodes[fn[i]].y)<threshold){
	  if (fabs(gnodes[fn[j]].z-gnodes[ln[i]].z)<threshold && fabs(gnodes[ln[j]].z-gnodes[fn[i]].z)<threshold){
	    gedges[k].nn=4;
	    gedges[k].nm=2;
	    gedges[k].fn=fn[i];
	    gedges[k].ln=ln[i];
	    gedges[k].nlist = new long [4];
	    gedges[k].nlist[0]=fn[i];
	    gedges[k].nlist[1]=ln[j];
	    gedges[k].nlist[2]=ln[i];
	    gedges[k].nlist[3]=fn[j];
	    k++;
	  }
	}
      }
    }
  }
  
  //for (i=0;i<nged;i++){
  //gedges[i].print (out);
  //}
  
  delete [] ln;
  delete [] fn;
  for (i=0;i<nn;i++){
    delete [] av[i];
  }
  delete [] av;

}

/**
   function computes direction vectors of edges
   
   JK, 11.7.2007
*/
void gtopology::edge_dirvect ()
{
  long i;
  
  for (i=0;i<nged;i++){
    gedges[i].direction_vector (gnodes);
  }
}

/**
   function computes direction vectors of edges
   
   JK, 11.7.2007
*/
void gtopology::edge_normvect ()
{
  long i;
  
  for (i=0;i<nged;i++){
    gedges[i].normal_vector (gnodes);
  }
}

/**
   function sorts first and last nodes on edges
   
   JK, 12.7.2007
*/
void gtopology::edgenode_sorting ()
{
  long i,j,nc;
  
  nc=0;
  for (i=0;i<nged;i++){
    if (gedges[i].nm != 2 || gedges[i].nn != 4){
      print_err(" wrong number of node multiplicity on edge", __FILE__, __LINE__, __func__);
      abort ();
    }
    if (gedges[i].nlist[0]>gedges[i].nlist[1]){
      j=gedges[i].nlist[1];
      gedges[i].nlist[1]=gedges[i].nlist[0];
      gedges[i].nlist[0]=j;
      nc++;
    }
    if (gedges[i].nlist[2]>gedges[i].nlist[3]){
      j=gedges[i].nlist[3];
      gedges[i].nlist[3]=gedges[i].nlist[2];
      gedges[i].nlist[2]=j;
      nc++;
    }
  }
  
  fprintf (stdout,"\n number of corrections in function edgenode_sorting is %ld",nc);
}


/**
   function searches for previous and next edges for each edge
   
   JK, 12.7.2007
*/
void gtopology::prev_next_edges ()
{
  long i,j;
  long ifn1,ifn2,iln1,iln2,jfn1,jfn2,jln1,jln2;
  
  for (i=0;i<nged;i++){
    if (gedges[i].nm != 2 || gedges[i].nn != 4){
      print_err(" wrong number of node multiplicity on edge", __FILE__, __LINE__, __func__);
      abort ();
    }
    
    ifn1=gedges[i].nlist[0];
    ifn2=gedges[i].nlist[1];
    iln1=gedges[i].nlist[2];
    iln2=gedges[i].nlist[3];

    for (j=i+1;j<nged;j++){
      jfn1=gedges[j].nlist[0];
      jfn2=gedges[j].nlist[1];
      jln1=gedges[j].nlist[2];
      jln2=gedges[j].nlist[3];
      
      if (iln1==jfn1 && iln2==jfn2){
	gedges[i].next=j;
	gedges[j].prev=i;
      }
      if (ifn1==jln1 && ifn2==jln2){
	gedges[i].prev=j;
	gedges[j].next=i;
      }
      if (ifn1==jfn1 && ifn2==jfn2){
	gedges[i].prev=j;
	gedges[j].prev=i;
      }
      if (iln1==jln1 && iln2==jln2){
	gedges[i].next=j;
	gedges[j].next=i;
      }
      
    }
  }
  
}



/**
   function searches for series of edges
   
   JK, 12.7.2007
*/
void gtopology::edge_series (FILE *out)
{
  long i,j,k,n,c,e,a,m;
  long *av;

  //  number of series
  nser=0;

  //  list of correspondence between edges and series
  edgeser = new long [nged];
  for (i=0;i<nged;i++){
    edgeser[i]=-1;
  }
  
  if (nged>0){
    
    //  series of egdes which contain only one edge
    c=0;
    for (i=0;i<nged;i++){
      if (gedges[i].prev==-1 && gedges[i].next==-1){
	edgeser[i]=nser;
	nser++;
	c++;
      }
    }
    
    a=-2; n=-2;
    for (i=0;i<nged;i++){
      if (gedges[i].prev==-1 && edgeser[i]==-1){
	//  untouched edge with no predecessor
	a=i;
	n=gedges[i].next;
	edgeser[i]=nser;
	break;
      }
    }
    if (a==-2 && n==-2){
      //  there is no edge with no predecessor
      for (i=0;i<nged;i++){
	if (gedges[i].next==-1 && edgeser[i]==-1){
	  a=i;
	  n=gedges[i].prev;
	  edgeser[i]=nser;
	  
	  gedges[i].next=gedges[i].prev;
	  gedges[i].prev=-1;
	  m=gedges[i].nlist[0];
	  gedges[i].nlist[0]=gedges[i].nlist[2];
	  gedges[i].nlist[2]=m;
	  m=gedges[i].nlist[1];
	  gedges[i].nlist[1]=gedges[i].nlist[3];
	  gedges[i].nlist[3]=m;
	  
	  m=gedges[i].fn;
	  gedges[i].fn=gedges[i].ln;
	  gedges[i].ln=m;
	  
	  break;
	}
      }
    }
    if (a==-2 && n==-2){
      i=0;
      a=i;
      n=gedges[i].next;
      edgeser[i]=nser;
    }
    c++;
    
    if (c<nged){
      do{
	e=gedges[n].next;
	if (edgeser[n]==-1){
	  if (a==e){
	    gedges[n].next=gedges[n].prev;
	    gedges[n].prev=a;
	    m=gedges[n].nlist[0];
	    gedges[n].nlist[0]=gedges[n].nlist[2];
	    gedges[n].nlist[2]=m;
	    m=gedges[n].nlist[1];
	    gedges[n].nlist[1]=gedges[n].nlist[3];
	    gedges[n].nlist[3]=m;
	    
	    m=gedges[n].fn;
	    gedges[n].fn=gedges[n].ln;
	    gedges[n].ln=m;
	    
	    a=n;
	    n=gedges[n].next;
	  }
	  else{
	    a=n;
	    n=e;
	  }
	  edgeser[a]=nser;
	  c++;
	}else{
	  n=-1;
	}
	
	if (n<0){
	  j=0;
	  for (i=0;i<nged;i++){
	    if (gedges[i].prev==-1 && edgeser[i]==-1){
	      a=i;
	      n=gedges[i].next;
	      nser++;
	      edgeser[i]=nser;
	      c++;
	      j=1;
	      break;
	    }
	  }
	  if (j==0){
	    for (i=0;i<nged;i++){
	      if (gedges[i].next==-1 && edgeser[i]==-1){
		
		a=i;
		n=gedges[i].prev;
		edgeser[i]=nser;
		
		gedges[i].next=gedges[i].prev;
		gedges[i].prev=-1;
		m=gedges[i].nlist[0];
		gedges[i].nlist[0]=gedges[i].nlist[2];
		gedges[i].nlist[2]=m;
		m=gedges[i].nlist[1];
		gedges[i].nlist[1]=gedges[i].nlist[3];
		gedges[i].nlist[3]=m;
		
		m=gedges[i].fn;
		gedges[i].fn=gedges[i].ln;
		gedges[i].ln=m;
		
		nser++;
		edgeser[i]=nser;
		c++;
		break;
	      }
	    }
	  }
	  if (j==0){
	    for (i=0;i<nged;i++){
	      if (edgeser[i]==-1){
		nser++;
		edgeser[i]=nser;
		c++;
		a=i;
		n=gedges[i].next;
		j=1;
	      }
	      if (j==1){
		break;
	      }
	    }
	  }
	}
      }while(c<nged);
      nser++;
    }
  }
  
  fprintf (stdout,"\n\n number of series %ld",nser);
  for (i=0;i<nser;i++){
    fprintf (out,"%ld\n",edgeser[i]);
  }
  fprintf (out,"\n");

  //  number of edges in series
  nedser = new long [nser];
  for (i=0;i<nser;i++){
    nedser[i]=0;
    for (j=0;j<nged;j++){
      if (edgeser[j]==i)
	nedser[i]++;
    }
  }
  
  
  
  
  
  edgelist = new long* [nser];
  for (i=0;i<nser;i++){
    edgelist[i] = new long [nedser[i]];
  }
  av = new long [nged];
  for (i=0;i<nged;i++){
    av[i]=-1;
  }
  
  for (i=0;i<nser;i++){
    
    //  non-closed series of edges are searched
    a=-1;
    for (j=0;j<nged;j++){
      if (av[j]>-1)
	continue;
      else{
	if (gedges[j].prev==-1){
	  a=j;
	  break;
	}
      }
    }
    if (a==-1){
      //  this case means that there are closed series of edges
      for (j=0;j<nged;j++){
	if (av[j]>-1){
	  continue;
	}
	else{
	  a=j;
	  break;
	}
      }
    }
    
    k=0;
    for (j=0;j<nged;j++){
      if (av[a]>-1){
	//  this statement is used for interruption of the loop in the case of closed series
	break;
      }
      edgelist[i][k]=a;
      k++;
      av[a]=1;
      n=gedges[a].next;
      a=n;
      if (n==-1){
	//  this statement is used for interruption of the loop in the case of non-closed series
	break;
      }
    }
  }
  
  
  fprintf (out,"\n\n\n kontrola pole edgelist \n");
  for (i=0;i<nser;i++){
    fprintf (out,"\n serie %5ld   ",i);
    for (j=0;j<nedser[i];j++){
      fprintf (out,"  %ld %ld %ld\n",edgelist[i][j],gedges[edgelist[i][j]].fn+1,gedges[edgelist[i][j]].ln+1);
    }
  }
  fprintf (out,"\n\n");
  
  delete [] av;
}

/**
   function assignes numbers of element to edges
   
   JK, 21.10.2007
*/
void gtopology::edge_elem (FILE *out)
{
  long i,j,k,l,m,fn,ln,ef,nne,nc,eid;
  ivector nod;

  //  list of adjacent elements to elements
  adjacelem (out);
  
  //  loop over number of series
  for (i=0;i<nser;i++){
    //  element not found
    ef=-1;
    
    //  loop over all general edges
    for (j=0;j<nedser[i];j++){
      
      m=edgelist[i][j];
      //  first node
      fn=gedges[m].fn;
      //  last node
      ln=gedges[m].ln;
      
      if (ef==-1){
	//  element containing node number fn and ln has to be searched for
	//  it is necessary to go through all elements until appropriate element is found
	for (k=0;k<ne;k++){
	  //  number of nodes on element
	  nne = give_original_nne (k);
	  //  array allocation
	  reallocv (RSTCKIVEC(nne,nod));
	  //  element nodes
	  give_original_nodes (k,nod);
	  
	  nc=0;
	  //  loop over nodes on element
	  for (l=0;l<nne;l++){
	    if (nod[l]==fn || nod[l]==ln)
	      nc++;
	  }
	  
	  if (nc==2){
	    gedges[m].re=k;
	    ef=k;
	    break;
	  }
	}
      }
      else{
	//  element containing node number fn and ln has to be searched for
	//  only adjacent elements are investigated
	
	//  loop over adjacent elements
	for (k=0;k<nadjelel[ef];k++){
	  //  element id
	  eid=adjelel[ef][k];
	  //  number of nodes on element
	  nne = give_original_nne (eid);
	  //  array allocation
	  reallocv (RSTCKIVEC(nne,nod));
	  //  element nodes
	  give_original_nodes (eid,nod);
	  
	  nc=0;
	  //  loop over nodes on element
	  for (l=0;l<nne;l++){
	    if (nod[l]==fn || nod[l]==ln)
	      nc++;
	  }

	  if (nc==2){
	    gedges[m].re=eid;
	    ef=eid;
	    break;
	  }
	  
	}
      }
      
    }
  }
}


/**
   function assignes numbers of all elements to edges
   
   JK, 6.8.2008
*/
void gtopology::edge_allelem (FILE */*out*/)
{
  long i,j,k,l,m,fn,ln,ef,nne,nc,eid;
  ivector nod;
  
  for (i=0;i<nged;i++){
    gedges[i].adjel = new long [2];
    gedges[i].adjel[0]=gedges[i].re;
  }
  
  //  loop over number of series
  for (i=0;i<nser;i++){
    //  element not found
    ef=-1;
    
    //  loop over all general edges
    for (j=0;j<nedser[i];j++){
      
      m=edgelist[i][j];
      //  first node
      fn=gedges[m].nlist[1];
      //  last node
      ln=gedges[m].nlist[3];
      
      if (ef==-1){
	//  element containing node number fn and ln has to be searched for
	//  it is necessary to go through all elements until appropriate element is found
	for (k=0;k<ne;k++){
	  //  number of nodes on element
	  nne = give_original_nne (k);
	  //  array allocation
	  reallocv (RSTCKIVEC(nne,nod));
	  //  element nodes
	  give_original_nodes (k,nod);
	  
	  nc=0;
	  //  loop over nodes on element
	  for (l=0;l<nne;l++){
	    if (nod[l]==fn || nod[l]==ln)
	      nc++;
	  }
	  
	  if (nc==2 && gedges[m].adjel[0]!=k){
	    gedges[m].adjel[1]=k;
	    ef=k;
	    break;
	  }
	}
      }
      else{
	//  element containing node number fn and ln has to be searched for
	//  only adjacent elements are investigated
	
	//  loop over adjacent elements
	for (k=0;k<nadjelel[ef];k++){
	  //  element id
	  eid=adjelel[ef][k];
	  //  number of nodes on element
	  nne = give_original_nne (eid);
	  //  array allocation
	  reallocv(RSTCKIVEC(nne,nod));
	  //  element nodes
	  give_original_nodes (eid,nod);
	  
	  nc=0;
	  //  loop over nodes on element
	  for (l=0;l<nne;l++){
	    if (nod[l]==fn || nod[l]==ln)
	      nc++;
	  }
	  
	  if (nc==2 && gedges[m].adjel[0]!=eid){
	    gedges[m].adjel[1]=eid;
	    ef=eid;
	    break;
	  }
	  
	}
      }
      
    }
  }
  
}



/**
   function checkes orientation of normal vectors
   
   JK, 21.10.2007
*/
void gtopology::normvectorient ()
{
  long i,eid,nne;
  ivector nod;
  vector x,y;
  
  for (i=0;i<nged;i++){
    //  number of reference element
    eid=gedges[i].re;
    //  number of nodes on element
    nne = give_original_nne (eid);
    //  allocation of vectors
    reallocv (RSTCKVEC(nne,x));
    reallocv (RSTCKVEC(nne,y));
    reallocv (RSTCKIVEC(nne,nod));
    
    //  coordinates of element node
    give_node_coord2d (x,y,eid);
    //  element nodes
    give_original_nodes (eid,nod);
    
    //  check of normal vectors
    gedges[i].check_normal (x,y,nod);
  }
}


/**
   function modifies indicators of code numbers
   in problems of hemivariational inequalities
   
   indicators have to be modified with respect to
   
*/
/*
void gtopology::cnmodif (FILE *out)
{
  long i,j,k,l,max,acn,ndofn;
  long **av=NULL;
  
  fprintf (out,"\n\n kontrola indikatoru kodovych cisel \n");
  for (i=0;i<nn;i++){
    fprintf (out,"\n node %6ld",i+1);
    ndofn = give_ndofn (i);
    for (j=0;j<ndofn;j++){
      fprintf (out,"  %ld",give_dof (i,j));
    }
  }
  
  max=0;
  acn=0;
  for (i=0;i<nn;i++){
    if (max<stop->ltg1[i])
      max=stop->ltg1[i];
    ndofn = give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=give_dof (i,j);
      if (k>acn)
	acn=k;
    }
  }
  max++;
  acn++;
  
  av = new long* [max];
  
  for (i=0;i<nn;i++){
    l=stop->ltg1[i];
    if (l>-1){
      ndofn = give_ndofn (i);
      av[l] = new long [ndofn];
    }
  }
  
  
  
  for (i=0;i<nn;i++){
    l=stop->ltg1[i];
    if (l>-1){
      ndofn = give_ndofn (i);
      for (j=0;j<ndofn;j++){
	k = give_dof (i,j);
	if (k==1){
	  av[l][j]=acn;
	  acn++;
	}
	else{
	  av[l][j]=k;
	}
      }
    }
  }
  
  for (i=0;i<nn;i++){
    l=stop->ltg1[i];
    if (l>-1){
      ndofn = give_ndofn (i);
      for (j=0;j<ndofn;j++){
	save_dof (i,j,av[l][j]);
      }
    }
  }
  
  for (i=0;i<max;i++){
    delete [] av[i];
  }
  delete [] av;

  fprintf (out,"\n\n kontrola indikatoru kodovych cisel \n");
  for (i=0;i<nn;i++){
    fprintf (out,"\n node %6ld",i+1);
    ndofn = give_ndofn (i);
    for (j=0;j<ndofn;j++){
      fprintf (out,"  %ld",give_dof (i,j));
    }
  }

}
*/

void gtopology::edges_auxprint (FILE *out)
{
  long i,j,a;
  
  fprintf (out,"\n\n\n kontrola pole nlist na jednotlivych hranach \n");
  for (i=0;i<nged;i++){
    fprintf (out,"\n %6ld %6ld       %6ld %6ld",gedges[i].nlist[0]+1,gedges[i].nlist[1]+1,gedges[i].nlist[2]+1,gedges[i].nlist[3]+1);
  }
  
  
  
  fprintf (out,"\n\n\n kontrola usporadni hran \n");
  for (i=0;i<nged;i++){
    fprintf (out,"\n %6ld       %6ld",gedges[i].nlist[0]+1,gedges[i].nlist[1]+1);
  }
  i=nged-1;
  if (i>=0)
    fprintf (out,"\n %6ld       %6ld",gedges[i].nlist[2]+1,gedges[i].nlist[3]+1);
  
  
  
  fprintf (out,"\n\n\n kontrola prev a next \n");
  for (i=0;i<nged;i++){
    fprintf (out,"\n hrana %6ld   prev %6ld  next %6ld",i,gedges[i].prev+1,gedges[i].next+1);
  }
  
  
  fprintf (out,"\n\n\n kontrola souslednosti \n");
  for (i=0;i<nser;i++){
    fprintf (out,"\n\n serie %5ld",i);
    for (j=0;j<nedser[i];j++){
      a=edgelist[i][j];
      fprintf (out,"\n %5ld  hrana %5ld   %5ld %5ld     %f %f %f      %f %f %f",
	       j,a,gedges[a].fn,gedges[a].ln,
	       gedges[a].dv[0],gedges[a].dv[1],gedges[a].dv[2],
	       gedges[a].nv[0],gedges[a].nv[1],gedges[a].nv[2]);
    }
  }
  fprintf (out,"\n\n");

  fprintf (out,"\n\n\n kontrola prilehlych prvku \n");
  for (i=0;i<nser;i++){
    fprintf (out,"\n\n serie %5ld",i);
    for (j=0;j<nedser[i];j++){
      a=edgelist[i][j];
      fprintf (out,"\n hrana %5ld  %6ld %6ld",j,gedges[a].adjel[0],gedges[a].adjel[1]);
    }
  }
  fprintf (out,"\n\n");
  
  
  /*
  fprintf (out,"\n\n\n kontrola souslednosti \n");
  long *av;
  av = new long [Gtm->nged];
  
  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;
      
      fprintf (out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      fprintf (out,"\n %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2]);
      //fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[0]);
      
      av[i]=1;
      break;
    }
  }
  long ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (out,"\n %6ld",Gtm->gedges[a].nlist[0]);
    //fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2]);
    fprintf (out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (out,"\n %6ld",Gtm->gedges[a].nlist[0]);
	  fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  fprintf (out,"\n %6ld",Gtm->gedges[a].nlist[2]);
  */
  
  /*
  fprintf (Out,"\n");


  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;

      //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      //fprintf (Out,"\n %6ld %6ld ",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
      fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);

      av[i]=1;
      break;
    }
  }
  ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);
    //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);
	  //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[3]-205);
  */


  
  /*
  fprintf (Out,"\n");


  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;

      fprintf (Out,"\n %lf %lf   %lf %lf ",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
      //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      //fprintf (Out,"\n %6ld %6ld ",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);

      av[i]=1;
      break;
    }
  }
  ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    //fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
    //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    fprintf (Out,"\n %lf %lf   %lf %lf",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
	  //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  fprintf (Out,"\n %lf %lf   %lf %lf",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
	  break;
	}
      }
    }
  }

  */
  
  /*
  a=0; n=Gtm->gedges[a].next;
  fprintf (out,"\n\n\n kontrola serie \n");
  for (i=0;i<Gtm->nged;i++){
    fprintf (out,"\n%ld -> %ld",a,n);
    a=n;
    n=Gtm->gedges[a].next;
  }
  */
  
  /*
  fprintf (out,"\n\n kontrola souslednosti\n");
  long *av;
  av = new long [Gtm->nged];
  
  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;
      
      fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
      fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      
      av[i]=1;
      break;
    }
  }
  long ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
    fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
	  fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  

  fprintf (Out,"\n\n\n");

  */



  /*
  for (i=0;i<Gtm->nged;i++){
    fprintf (out,"\n\n hrana cislo  %ld",i);
    Gtm->gedges[i].print (out);
  }
  */
  
  /*
  fprintf (out,"\n\n\n kontrola serii \n");
  for (i=0;i<Gtm->nged;i++){
    fprintf (out,"\n edge %6ld je v serii  %6ld",i,Gtm->edgeser[i]);
  }
  */

}


/**
   function reads informations necessary for sequential FETI method
   problem is solved by FETI method on one processor
   
   @param ns - number of subdomains
   @param in - input file
   
   JK, 6.10.2007
*/
/*
void gtopology::read_feti (long ns,XFILE *in)
{
  meshdescription md;
  
  //  type of mesh description - see GEFEL/galias.h
  xfscanf (in,"%k%m","meshdescript",&meshdescription_kwdset,(int*)&md);
  
  //  instance of the class seqtop (sequential topology)
  stop = new seqtop (ns,md);
  
  //  number of nodes on subdomains
  stop->read_nnsd (in);
  
  //  reading of auxiliary array ltg (local to global map)
  stop->read_ltg (in);
}
*/
/**
   function reads array ltg1
   
   @param in - input file
   
   JK, 22.11.2007
*/
/*
void gtopology::read_ltg1 (XFILE *in)
{
  //  instance of the class seqtop (sequential topology)
  //  ns=1 - number of subdomains
  //  md=1 - all_nodes
  stop = new seqtop (1,all_nodes);
  
  //  reading of auxiliary array ltg1 (local to global map)
  stop->readltg1 (nn,in);
}
*/



/**
  The function allocates array of time function indices at nodes
  where no time function was specified. In such nodes, the dofs are 
  controlled by time functions from the adjacent elements 
  (i.e by gtopology::lnso state indicator) and the time function indices
  are set to -1 for all dofs.

  @return The function does not return anything but it changes tgf of some nodes.

  Created by TKo, 5.2013
*/
void gtopology::alloc_growstr()
{
  long i,j,ndofn;

  for(i=0; i<nn; i++)
  {
    if (gnodes[i].tgf == NULL) // if no time functions were specified 
    {
      // allocate array of indices of time functions for particular dofs      
      ndofn = gnodes[i].give_ndofn();
      if (ndofn > 0) // i-th node is not hanging node
      {
        gnodes[i].tgf = new long[ndofn];      
        // for all dofs set time function index to -1
        // which represents time function taken from adjacent elements
        for(j=0; j<ndofn; j++)
          gnodes[i].tgf[j] = -1;
      }
    }
  }
}



/**
   function allocates arrays cn on edges
   cn contains code numbers of Lagrange multipliers
   
   JK, 6.1.2008
*/
void gtopology::alloc_edge_cn ()
{
  long i,nid,ndofnfn,ndofnln;
  
  for (i=0;i<nged;i++){
    //  first node on edge
    nid=gedges[i].fn;
    //  number of DOFs on the first node
    ndofnfn=give_ndofn (nid);
    //  last node on edge
    nid=gedges[i].ln;
    //  number of DOFs on the last node
    ndofnln=give_ndofn (nid);
    
    //  allocation of arrays cnfn and cnln
    gedges[i].alloc_cn (ndofnfn,ndofnln);
  }

}


/**
   function localizes contributions from edges to the global %matrix
   
   @param gm - pointer to the general %matrix
   
   JK, 8.8.2008
*/
void gtopology::edge_localization (gmatrix *gm)
{
  long i,nm,*ncn1,*ncn2,*mcn;

  for (i=0;i<nged;i++){
    //  number of multipliers between first nodes
    nm=gedges[i].ndofnf;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];

    //  code numbers on first nodes and appropriate multipliers
    give_edge_code_numbers (i,1,ncn1,ncn2,mcn);
    
    //  localization of contributions from edges
    gm->mult_localize (nm,ncn1,ncn2,mcn);
    
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
    
    
    //  number of multipliers between second nodes
    nm=gedges[i].ndofnl;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];

    //  code numbers on first nodes and appropriate multipliers
    give_edge_code_numbers (i,2,ncn1,ncn2,mcn);
    
    //  localization of contributions from edges
    gm->mult_localize (nm,ncn1,ncn2,mcn);
    
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }

}






















/**
   function automatically detects subdomains in conforming meshes
   
   it is used e.g. in problems with hemivariational inequalities,
   material discontinuities
   
   @param out - output file
   
   JK, 12.11.2007
*/
void gtopology::auto_subdom_detection (FILE *out)
{
  // ***********************
  //  contact end points  **
  // ***********************
  //  function allocates objects endnodes
  //  the objects contain node numbers and multipliers
  end_nodes (out);
  
  //  function stores neighbouring element numbers to the objects endnodes
  enodes_allelem (out);
  
  endnodes_auxprint (out);
  
  
  // ******************
  //  contact edges  **
  // ******************
  //  the actual value of the variable edtype is assigned in the function
  //  radiation_init () in the file trfelinit.cpp in TRFEL
  if (edtype==jumps){
    //  list of adjacent nodes to each node
    adjacnodes (out);
    //  assembling of interface edges - edges on contact between domains
    //  list of nodes on edges is created
    edges (out);
  }
  if (edtype==mater){
    //  list of adjacent elements to each element
    adjacelem (out);
    //  definition of the type of material model used for radiation
    //  this choice is temporary, it will be fixed in the future
    long matrad=251;
    //  assembling of interface edges - edges on material interfaces
    edges_mat (matrad,out);
  }
  
  //  edge nodes are sorted
  edgenode_sorting ();
  //  number of previous and next edges (index of previous and next edges)
  prev_next_edges ();
  //  assembling of edge series
  edge_series (out);
  //  edge-element correspondence
  edge_elem (out);


  //  edge-all elements correspondence
  edge_allelem (out);

  //  computation of direction/tangential vectors
  edge_dirvect ();
  //  computation of normal vectors
  edge_normvect ();
  //  check of normal vector orientation
  normvectorient ();
  
  
  if (edtype==jumps){
    //  function creates lists of nodes in subdomains
    domdecomp ();
    //  creation of array ltg (local to global ordering)
    //  ltg array is used for Schur ordering of unknowns
    create_ltg (out);
  }
  
  edges_auxprint (out);
  
  if (edtype==mater){
    cngen=1;
  }
  if (edtype==jumps){
    //  code numbers are generated with respect to saddle point problems
    cngen=3;
  }
  
}


/**
   function generates code numbers of Lagrange multipliers
   
   JK, 6.1.2008
*/
long gtopology::codenum_multip ()
{
  long i,j,k,l,ii,nid,nm,ndofn;
  long *cn = NULL;
  
  //  auxiliary variable
  nm=ndof;

  // ***********************************
  //  multipliers defined on end nodes
  // ***********************************
  
  //  allocation of array cnm on end nodes
  alloc_endnode_cn ();
  
  for (i=0;i<nen;i++){
    //  first node id
    nid=endnodes[i].fn;
    //  number of DOFs on the node
    ndofn=give_ndofn (nid);
    cn = new long [ndofn];
    for (l=0;l<ndofn;l++){
      ii=give_dof (nid,l);
      if (ii!=0){
	cn[l]=nm;
	nm++;
      }
      else
	cn[l]=0;
    }
    for (l=0;l<ndofn;l++){
      endnodes[i].cnm[l]=cn[l];
    }
    delete [] cn;
  }
  
  // *******************************
  //  multipliers defined on edges
  // *******************************
  
  //  allocation of arrays cn on edges
  alloc_edge_cn ();
  
  for (i=0;i<nser;i++){
    //  loop over series of edges
    for (j=0;j<nedser[i];j++){
      //  loop over all edges in the serie
      
      //  number of actual edge
      k=edgelist[i][j];
      
      // **************
      //  first nodes
      // **************
      if (j==0){
	//  node id of the first node on the edge
	nid=gedges[k].fn;
	//  number of DOFs on the node
	ndofn=give_ndofn (nid);
	cn = new long [ndofn];
	memset(cn, 0, sizeof(*cn)*ndofn);
	for (l=0;l<ndofn;l++){
	  ii=give_dof (nid,l);
	  if (ii!=0){
	    cn[l]=nm;
	    nm++;
	  }
	  else
	    cn[l]=0;
	}
      }
      
      for (l=0;l<ndofn;l++){
	gedges[k].cnfn[l]=cn[l];
      }
      delete [] cn;
      
      // *************
      //  last nodes
      // *************
      //  node id of the first node on the edge
      nid=gedges[k].ln;
      //  number of DOFs on the node
      ndofn=give_ndofn (nid);
      cn = new long [ndofn];
      for (l=0;l<ndofn;l++){
	ii=give_dof (nid,l);
	if (ii!=0){
	  cn[l]=nm;
	  nm++;
	}
	else
	  cn[l]=0;
      }
      
      for (l=0;l<ndofn;l++){
	gedges[k].cnln[l]=cn[l];
      }
      
      if (j==nedser[i]-1)
	delete [] cn;
    }
  }
  
  return nm;
}

/**
   function localizes %matrix entries connected with
   Lagrange multipliers to the general %matrix
   there are contributions from end nodes (1D),
   edges (2D)
   
   @param gm - general %matrix
   
   JK, 14.10.2008
*/
void gtopology::mult_localization (gmatrix *gm)
{
  //  contributions from end nodes
  endnodes_localization (gm);
  //  contributions from edges
  edge_localization (gm);
}

/**
   function creates lists of nodes in subdomains
   mesh must be described by the glued ordering
   
   1.6.2008, JK
*/
void gtopology::domdecomp ()
{
  long i,j,ii,jj,nvn,nan,nnf,nnnf,anode;
  long *aux,*oldlist,*newlist;
  
  aux = new long [nn];
  for (i=0;i<nn;i++){
    aux[i]=-1;
  }
  
  //  number of subdomains
  ns=0;
  //  number of visited nodes
  nvn=0;
  do{
    for (i=0;i<nn;i++){
      if (aux[i]==-1)
	anode=i;
    }
    
    oldlist=new long[1];
    oldlist[0]=anode;
    nnf=1;
    
    do{
      //  number of nodes on a new front
      nnnf=0;
      //  loop over nodes on actual front
      //  nnf - number of nodes on actual front
      for (i=0;i<nnf;i++){
	//  number of actual node on the front
	ii=oldlist[i];
	//  number of adjacent nodes to actual node
	nan=nadjnodnod[ii];
	//  loop over adjacent nodes to actual node
	for (j=0;j<nan;j++){
	  //  number of adjacent node
	  jj=adjnodnod[ii][j];
	  if (aux[jj]==-1){
	    //  node has not been touched yet
	    nnnf++;
	  }
	}
      }
      //  list of nodes on new front
      newlist=new long [nnnf];
      //  number of nodes on a new front
      nnnf=0;
      //  loop over nodes on actual front
      //  nnf - number of nodes on actual front
      for (i=0;i<nnf;i++){
	//  number of actual node on the front
	ii=oldlist[i];
	//  number of adjacent nodes to actual node
	nan=nadjnodnod[ii];
	//  loop over adjacent nodes to actual node
	for (j=0;j<nan;j++){
	  //  number of adjacent node
	  jj=adjnodnod[ii][j];
	  if (aux[jj]==-1){
	    //  node has not been touched yet
	    newlist[nnnf]=jj;
	    nnnf++;
	    aux[jj]=ns;
	    nvn++;
	  }
	}
      }
      delete [] oldlist;
      oldlist=new long [nnnf];
      for (i=0;i<nnnf;i++){
	oldlist[i]=newlist[i];
      }
      delete [] newlist;
      nnf=nnnf;
    }while(nnf>0);

    delete [] oldlist;
    
    ns++;
  }while(nvn<nn);
  
  
  //  instance of the class seqtop (sequential topology)
  if (stop==NULL)
    stop = new seqtop (ns,all_nodes);
  
  //  allocation 
  if (stop->nnsd==NULL)
    stop->nnsd = new long [ns];
  for (i=0;i<ns;i++){
    stop->nnsd[i]=0;
  }
  
  for (i=0;i<nn;i++){
    stop->nnsd[aux[i]]++;
  }
  
  delete [] aux;
}

/**
   function assembles array ltg for sequential topology
   array is constructed from the knowledge of interfaces
   among subdomains
   function is intended for hemivariational inequalities
   
   JK, 6.11.2007
*/
void gtopology::create_ltg (FILE *out)
{
  long i,j,k,nin;
  long *aux,*v;

  aux = new long [nn];
  for (i=0;i<nn;i++){
    aux[i]=nn;
  }
  
  //  auxiliary array assembling
  //  aux[i]=j - the i-th node has label j
  //  if j=nn - noninterface node
  //  if j!=nn - interface node, j is equal to minimum number of node connected to this interface node
  //
  
  //  loop over end nodes
  for (i=0;i<nen;i++){
    k=endnodes[i].fn;
    j=endnodes[i].ln;
    if (j<k){
      if (aux[j]>j)
	aux[j]=j;
      if (aux[k]>j)
	aux[k]=j;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }
    else{
      if (aux[j]>k)
	aux[j]=k;
      if (aux[k]>k)
	aux[k]=k;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }
    
  }

  //  loop over general edges
  for (i=0;i<nged;i++){
    j=gedges[i].nlist[0];
    k=gedges[i].nlist[1];
    if (j<k){
      if (aux[j]>j)
	aux[j]=j;
      if (aux[k]>j)
	aux[k]=j;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }
    else{
      if (aux[j]>k)
	aux[j]=k;
      if (aux[k]>k)
	aux[k]=k;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }
    
    j=gedges[i].nlist[2];
    k=gedges[i].nlist[3];
    if (j<k){
      if (aux[j]>j)
	aux[j]=j;
      if (aux[k]>j)
	aux[k]=j;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }
    else{
      if (aux[j]>k)
	aux[j]=k;
      if (aux[k]>k)
	aux[k]=k;

      if (aux[j]<aux[k])
	aux[k]=aux[j];
      else
	aux[j]=aux[k];
    }

  }
  
  
  
  v = new long [nn];
  for (i=0;i<nn;i++){
    v[i]=-1;
  }
  
  
  for (i=0;i<nn;i++){
    j=aux[i];
    if (j!=nn)
      v[j]++;
  }
  
  //  generation of interface node numbers
  //  nin - number of interface nodes
  nin=0;
  for (i=0;i<nn;i++){
    if (v[i]>-1){
      v[i]=nin;
      nin++;
    }
  }
  
  for (i=0;i<nn;i++){
    j=aux[i];
    if (j!=nn)
      aux[i]=v[j];
    else
      aux[i]=-1;
  }

  
  //
  //  array ltg
  //
  if (stop->ltg==NULL){
    stop->ltg = new long* [stop->ns];
    for (i=0;i<stop->ns;i++){
      stop->ltg[i] = new long [stop->nnsd[i]];
    }
  }
  
  k=0;
  for (i=0;i<stop->ns;i++){
    for (j=0;j<stop->nnsd[i];j++){
      stop->ltg[i][j]=aux[k];
      k++;
    }
  }
  
  
  fprintf (out,"\n\n\n funkce seqtop::create_ltg   pole ltg");
  for (i=0;i<stop->ns;i++){
    fprintf (out,"\n");
    for (j=0;j<stop->nnsd[i];j++){
      fprintf (out,"%ld  ",stop->ltg[i][j]);
    }
  }
  fprintf (out,"\n");


  delete [] aux;
  delete [] v;
  
}

/**
   function allocates object of the class sequential topology
   function read data about meshes

   this function is used in sequential (one-processor) implementation
   of domain decomposition methods

   @param in - input file
   
   JK, 23.12.2007
*/
void gtopology::read_seq_top (XFILE *in)
{
  meshdescription md;
  
  //  type of mesh description - see GEFEL/galias.h
  //  md=1 - all_nodes = global node ordering
  //  md=2 - boundary nodes ordering
  //  md=3 - negative boundary nodes ordering
  //  md=4 - global glued ordering
  //  md=11 - metis
  //  md=12 - metis_elem
  xfscanf (in,"%k%m","meshdescript",&meshdescription_kwdset,(int*)&md);
  
  //  number of subdomains/aggregates
  xfscanf (in,"%ld",&ns);
  
  //  instance of the class seqtop (sequential topology)
  stop = new seqtop (ns,md);
  
  if (md==all_nodes || md==bound_nodes || md==neg_bound_nodes || md==metis){
    //  in the case of md==4 (global glued ordering), the numbers of nodes on
    //  subdomains and the array ltg will be determined automatically
    //  there is the function gtopology::auto_subdom_detection
    
    //  number of nodes on subdomains
    stop->read_nnsd (in);
    
    //  reading of auxiliary array ltg (local to global map)
    stop->read_ltg (in);
  }
  if (md==metis_elem){
    /*
    //  number of elements on subdomains
    stop->read_nesd (in);
    
    //  reading of auxiliary array eltg (local to global map for elements)
    stop->read_eltg (in);
    */
    
    //  reading of the array eldom
    //  this line is avoided only for purposes of Charles bridge
    // zde rozmyslet??!!
    stop->read_eldom (ne,in);
  }

}




















/**
   function reorders node numbering with respect to Cuthill-McKee algorithm
   
   JB, 27.2.2007
   changed JB, 24.10.2008
*/

// void gtopology::cuthill_mckee_renumb (FILE *out)
// {
//   long i,j,k;
//   long minnode,n,a,b,m,stop;
//   int startnode,endnode;
//   long *cmk,**auxadjnodnod,*auxnadjnodnod,*aux; 
  
//   //fprintf (stdout,"\n\n\n CUTHILL MCKEE\n\n");
//   // nalezeni sousednich uzlu k uzlu a jejich poctu
//   adjacnodes (out);
  
//   // pole precislovanych uzlu a zaroven navstivenych uzlu
//   cmk = new long [nn];
//   for(i = 0; i<nn; i++)  cmk[i] = 0;
  
  
//   findPseudoPheripheral(&startnode,&endnode);
//   minnode = long(startnode);
//   //fprintf (stdout,"\n\n\nminnode = %ld\n\n",minnode+1);
  
//   // nalezni minima v nadjnodnod pro prvni uzel
//   //min = LONG_MAX;
//   //for(i = 0; i<nn; i++){
//     //fprintf (stdout,"%ld %ld\n",i+1, nadjnodnod[i]);
//     //if(min >= nadjnodnod[i]){
//   //min = nadjnodnod[i];
//   //  minnode = i;
//   //}
//   //}
//   //fprintf (stdout,"\n\n\nminnode = %ld\n\n",minnode); 
  
  
//   //precislovani prvniho node - cislovani od 1
//   cmk[minnode] = 1;
  
//   //serazeni pole adjnodnod podle velikosti sousedu
//   // plneni pomocneho pole
//   // pomocne pole
//   auxadjnodnod = new long* [nn];
//   auxnadjnodnod = new long [nn];
  
//   for(i = 0; i < nn; i++){
//     auxnadjnodnod[i] = nadjnodnod[i]-1;
//   }
  
//   for(i = 0; i < nn; i++){
//     auxadjnodnod[i] = new long [auxnadjnodnod[i]];
//     k = 0;
//     //fprintf (stdout,"%ld  %ld -> ",i,auxnadjnodnod[i]);
//     for(j = 0; j < nadjnodnod[i]; j++){
//       if( i != adjnodnod[i][j]){
// 	auxadjnodnod[i][k] =  adjnodnod[i][j];
// 	//fprintf (stdout,"%ld-%ld   ",auxadjnodnod[i][k],auxnadjnodnod[auxadjnodnod[i][k]]);
// 	k++;
//       }
//     }
//     //fprintf (stdout,"\n");
    
//     //buble sorting pole auxadjnodnod podle velikosti sousedu
//     for (m =auxnadjnodnod[i] - 1; m >= 0;  m--){
//       for (j = 1; j <= m; j++){
// 	a = auxadjnodnod[i][j];
// 	b = auxadjnodnod[i][j-1];
// 	if (auxnadjnodnod[b] > auxnadjnodnod[a]){
// 	  n = auxadjnodnod[i][j-1];
// 	  auxadjnodnod[i][j-1] = auxadjnodnod[i][j];
// 	  auxadjnodnod[i][j] = n;
// 	}
//       }
//     }
    
//     // for(j = 0; j < auxnadjnodnod[i]; j++){
//     // fprintf (stdout,"      %ld-%ld",auxadjnodnod[i][j],auxnadjnodnod[auxadjnodnod[i][j]]);
//     //}
//     //fprintf (stdout,"\n");
//   }

  
//   aux = new long[nn];
//   for(i = 0; i < nn; i ++){
//     aux[i]=-1;
//   }
//   cmk = new long[nn];
//   startnode = 0;
//   aux[minnode] = 0;
//   stop = 0;
//   m = 1;
//   for(k = 0; k < nn; k++){
//     //fprintf (stdout,"startnode-%ld\n",startnode+1);
//     for(i = 0; i < nn; i ++){ 
//       if(aux[i] == startnode){
// 	//fprintf (stdout,"jsem na %ld\n",i+1);
// 	for(j = 0; j < auxnadjnodnod[i];j++){
// 	  if( aux[auxadjnodnod[i][j]] == -1){
// 	    cmk[auxadjnodnod[i][j]] = m;
// 	    aux[auxadjnodnod[i][j]] = m;
// 	    //fprintf (stdout,"%ld-%ld\n",auxadjnodnod[i][j]+1,m+1);
// 	    m++;
// 	  }
// 	  if(m == nn){
// 	    stop = 1;
// 	  }
// 	}
// 	break;
//       }
//     }
//     if(stop == 1){
//       break;
//     }
//     startnode++;
//   }


  
  
//   // mazani pomocneho pole 
//   for(i = 0; i < nn; i++) delete[] auxadjnodnod[i];
//   delete []auxadjnodnod;
//   delete []auxnadjnodnod;
//   delete []aux;
  
//   //  cmk[i]=j - the i-th node of the original ordering has number j in the new ordering
  
//   if (nodren == cuthill_mckee){
//     for (i=0;i<nn;i++){
//       ordering[cmk[i]]=i;
//     }
//   }

//   if (nodren == rev_cuthill_mckee){
//     for (i=0;i<nn;i++){
//       ordering[cmk[i]]=i;
//     }
//     j=nn-1;
//     for (i=0;i<nn;i++){
//       cmk[j]=ordering[i];
//       j--;
//     }
//     for (i=0;i<nn;i++){
//       ordering[i]=cmk[i];
//     }
    
//   }
  
//   //fprintf (stdout,"Kontrolni tisk ordering\n"); 
//   //for(i = 0; i < nn; i++)  fprintf (stdout,"%ld %ld\n",i+1,ordering[i]+1);

// }


/**
   function reorders node numbering with respect to Cuthill-McKee algorithm
   
   JB, 27.2.2007
   changed JB, 24.10.2008
*/

void gtopology::cuthill_mckee_renumb (FILE *out)
{
  long i,j,k;
  long minnode,n,a,b,m;
  int startnode;
  long min;
  long *cmk;
  bool stop;
  
  //fprintf (stdout,"\n\n\n CUTHILL MCKEE\n\n");
  //establishing nodal graph for reordering
  adjacnodes_edge(out);
  
  // allocation of array for renumbered nodes and vistited nodes
  cmk = new long [nn];
  for(i = 0; i<nn; i++)  cmk[i] = -1;
  
  
  // searching starting node with minimal degree
  min = LONG_MAX;
  for(i = 0; i<nn; i++){
    if(min >= nadjacnodesedge[i]){
      min = nadjacnodesedge[i];
      minnode = i;
    }
  }
  fprintf (stdout,"\n\n\nminnode = %ld\n\n",minnode); 
  
  
  //renumbering of the starting node
  cmk[minnode] = 0;
  
  // sorting of array a
  for(i = 0; i < nn; i++){
    //buble sorting pole auxadjnodnod podle velikosti sousedu
    for (m = nadjacnodesedge[i] - 1; m >= 0;  m--){
      for (j = 1; j <= m; j++){
	a = adjacnodesedge[i][j];
	b = adjacnodesedge[i][j-1];
	if (nadjacnodesedge[b] > nadjacnodesedge[a]){
	  n = adjacnodesedge[i][j-1];
	  adjacnodesedge[i][j-1] = adjacnodesedge[i][j];
	  adjacnodesedge[i][j] = n;
	}
      }
    }
  }
  for(i = 0; i < nn; i++){
    fprintf(stdout,"node %ld # neighbours %ld:",i,nadjacnodesedge[i]);
    for(j= 0; j < nadjacnodesedge[i]; j++){
      fprintf(stdout,"  %ld",adjacnodesedge[i][j]);
    }
    fprintf(stdout,"\n");
  }
  
  startnode = 0;
  stop = 0;
  m = 1;
  stop = false;
  for(k = 0; k < nn; k++){
    for(i = 0; i < nn; i++){ 
      if(cmk[i] == startnode){
	//fprintf(stdout,"jsem na uzlu %ld\n",i);
	for(j = 0; j < nadjacnodesedge[i];j++){
	  if(cmk[adjacnodesedge[i][j]] == -1){
	    //fprintf(stdout," %ld\n",adjacnodesedge[i][j]);
	    cmk[adjacnodesedge[i][j]] = m;
	    //fprintf(stdout," %ld\n",cmk[adjacnodesedge[i][j]]);
	    m++;
	  }
	  if(m == nn){
	    stop = true;
	    break;
	  }
	}
	break;
      }
    }
    if(stop == true){
      break;
    }
    startnode++;
  }


  
  
  // deleting of arrays with graph
  for(i = 0; i < nn; i++) delete[] adjacnodesedge[i];
  delete []adjacnodesedge;
  delete []nadjacnodesedge;

  
  //  cmk[i]=j - the i-th node of the original ordering has number j in the new ordering
  
  if (nodren == cuthill_mckee){
    for (i=0;i<nn;i++){
      ordering[cmk[i]]=i;
    }
  }

  if (nodren == rev_cuthill_mckee){
    for (i=0;i<nn;i++){
      ordering[cmk[i]]=i;
    }
    j=nn-1;
    for (i=0;i<nn;i++){
      cmk[j]=ordering[i];
      j--;
    }
    for (i=0;i<nn;i++){
      ordering[i]=cmk[i];
    }
    
  }
  
  // fprintf (stdout,"Kontrolni tisk ordering\n"); 
//   for(i = 0; i < nn; i++)  fprintf (stdout,"%ld %ld\n",i+1,cmk[i]);

//   for (i=0;i<nn;i++){
//     ordering[i]=i;
//   }
}





void gtopology::writePriorityQueue(int */*fronta*/){
  for (int i = 0; i < (nn); i++){
    //printf( "%d. hodnota %d\n", i, fronta[i]);
  }
}

void gtopology::shell_sort( int *array, int arrayLength)
{
  int flag = 1, d = arrayLength, i, temp;
  while( flag || (d>1))      // boolean flag (true when not equal to 0)
    {
      flag = 0;           // reset flag to 0 to check for future swaps
      d = (d+1) / 2;
      for (i = 0; i < (arrayLength - d); i++)
	{
	  if (nadjnodnod[array[i + d]] < nadjnodnod[array[i]])
	    {
	      temp = array[i + d];      // swap items at positions i+d and d
	      array[i + d] = array[i];
	      array[i] = temp;
	      flag = 1;                  // indicate that a swap has occurred
	    }
	}
    }
  return;
}



/**
  The function sorts gtopology::nodes array by Shell sort algorithm in ascending order 
  with respect to x coordinate.

  @param x_ord - vector of sorted nodal ids

  @return The function returns resulting x order in the parameter x_ord.

  Created by Tomas Koudelka, 22.10.2013
*/
void gtopology::shell_sort_x(ivector &x_ord)
{
  long flag = 1, d = nn, i, temp;
  
  reallocv(nn, x_ord);  // this MUST be allocated at heap
  for(i=0; i<nn; i++)
    x_ord[i] = i;

  while(flag || (d>1))   // boolean flag (true when not equal to 0)
  {
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1)/2;
    for (i=0; i<(nn-d); i++)
    {
      if (gnodes[i+d].x < gnodes[i].x)
      {
        temp   = x_ord[i+d];   // swap items at positions i+d and d
        x_ord[i+d] = x_ord[i];
        x_ord[i]   = temp;
        flag   = 1;            // indicate that a swap has occurred
      }
    }
  }
  return;
}



/**
  The function sorts gtopology::nodes array by Shell sort algorithm in ascending order 
  with respect to y coordinate.

  @param y_ord - vector of sorted nodal ids

  @return The function returns resulting y order in the parameter y_ord.

  Created by Tomas Koudelka, 22.10.2013
*/
void gtopology::shell_sort_y(ivector &y_ord)
{
  long flag = 1, d = nn, i, temp;
  
  reallocv(nn, y_ord); // this MUST be allocated at heap
  for(i=0; i<nn; i++)
    y_ord[i] = i;

  while(flag || (d>1))   // boolean flag (true when not equal to 0)
  {
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1)/2;
    for (i=0; i<(nn-d); i++)
    {
      if (gnodes[i+d].y < gnodes[i].y)
      {
        temp   = y_ord[i+d];   // swap items at positions i+d and d
        y_ord[i+d] = y_ord[i];
        y_ord[i]   = temp;
        flag   = 1;            // indicate that a swap has occurred
      }
    }
  }
  return;
}



/**
  The function sorts gtopology::nodes array by Shell sort algorithm in ascending order 
  with respect to z coordinate.

  @param z_ord - vector of sorted nodal ids

  @return The function returns resulting z order in the parameter z_ord.

  Created by Tomas Koudelka, 22.10.2013
*/
void gtopology::shell_sort_z(ivector &z_ord)
{
  long flag = 1, d = nn, i, temp;
  
  reallocv(nn, z_ord); // this MUST be allocated at heap
  for(i=0; i<nn; i++)
    z_ord[i] = i;

  while(flag || (d>1))   // boolean flag (true when not equal to 0)
  {
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1)/2;
    for (i=0; i<(nn-d); i++)
    {
      if (gnodes[i+d].z < gnodes[i].z)
      {
        temp   = z_ord[i+d];   // swap items at positions i+d and d
        z_ord[i+d] = z_ord[i];
        z_ord[i]   = temp;
        flag   = 1;            // indicate that a swap has occurred
      }
    }
  }
  return;
}




void gtopology::lastLevelOfRootedStructure(int start, int *velikost, int *maxSirka, int *hloubka, int *vzdalenost, int **lastLevel){
  int *aktUzlyStart;
  int i,j,k; // pomocne promenne cyklu
  aktUzlyStart = new int[nn];
  for(i = 0; i< nn; i++){	
    aktUzlyStart[i] = 0;
  }
  int pocetObjevenych = 0; // pocet neaktivnich uzlu nalezenych v okoli aktualne zpracovaneho uzlu
  int *hranice;
  int velikostHranice;
  int aktualniUzel,pomUzel,posl;
  
  *maxSirka = -1;
  *hloubka = 0;
  aktUzlyStart[start]=1;
  velikostHranice = 1;
  hranice = new int[1];
  hranice[0] = start;
  int *pomocnePole;
  pomocnePole = new int[nn];
  do{
    
    
    //vynulujeme pomocne pole
    for( j=0; j<nn; j++){
      pomocnePole[j] = 0;
      
    }
    posl =0;
    pocetObjevenych = 0;
    //projdeme vsechny prvky na hranici
    for( j=0; j<velikostHranice; j++){
      
      aktualniUzel = hranice[j];
      //printf("Zpracovavam uzel %d\n",aktualniUzel);
      //projdem vsechny neaktivni uzly na hranici a pridame je do hranice
      for(i = 0; i< nadjnodnod[aktualniUzel]; i++){
	pomUzel = adjnodnod[aktualniUzel][i];
	vzdalenost[aktualniUzel] = *hloubka;
	//printf("\tZpracovavam souseda %d ",pomUzel);
	if(aktUzlyStart[pomUzel] == 0){
	  //printf("platny uzel");
	  
	  aktUzlyStart[pomUzel] = 1;
	  pomocnePole[posl] = pomUzel;
	  pocetObjevenych++;
	  posl++;
	  
	}
	//printf("\n");
      }
    }
    if(pocetObjevenych > *maxSirka){
      *maxSirka = pocetObjevenych;
    }
    (*hloubka)++; // v kazdem pruchodu se hloubka zvetsuje
    //getchar();
    //zkopirujeme pomocne pole do hranice
    if(pocetObjevenych != 0){
      delete[] hranice;
      hranice = new int[posl + 1];
      for( k=0; k<posl; k++){
	hranice[k]= pomocnePole[k];
      }
      velikostHranice = posl ;
      //printf("Pocet objevenych %d " , pocetObjevenych);
    }
  }while(pocetObjevenych != 0);
  //
  //for(k=0; k<velikostHranice; k++){
  //		printf("Prvek %d a stupen %d\n" , hranice[k],nadjnodnod[hranice[k]]);
  //	}
  ///printf("velikostHranice =- %d", velikostHranice);
  if(*lastLevel != NULL){
	  delete [] *lastLevel;
  }
    *lastLevel = hranice;
    *velikost = velikostHranice;
    
}

void gtopology::generateInitialList(int start, int **seznam, int * velikostSeznamu, int *hloubkaSeznamu, int *maxSirka){
  
  //int minUzel;
  int sirka, hloubka,velikost;
  int *lastLevel=NULL;
  int *vzdalenosti;
  
  
  vzdalenosti = new int[nn];
  //int i,j,k; // pomocne promenne cyklu
  // nalezeni sousednich uzlu k uzlu a jejich poctu
  
  lastLevelOfRootedStructure(start,&velikost,&sirka , &hloubka, vzdalenosti , &lastLevel);
  
  
  
  
  shell_sort(lastLevel,velikost);
  removeRepetitiousDegrees(&lastLevel, &velikost);
  
  *velikostSeznamu = velikost;
  *hloubkaSeznamu = hloubka;
  *maxSirka = sirka;
  *seznam = lastLevel ;
}


// funkce najde uzel s minimalnim stupnem
int gtopology::findMinimumDegree(){
  int minUzel,minimum = 99999999;
  for(int i = 0; i< nn; i++){
    
    //printf("MINIMAL Prvek %d a stupen je %d\n" , i,nadjnodnod[i]);
    
    if(nadjnodnod[i]< minimum){
      minUzel = i;
      minimum = nadjnodnod[i];
    }
  }
  return minUzel;
}



void gtopology::removeRepetitiousDegrees(int **pole, int *velikost){
  int aktStupen= -1,k,posl;
  int *pomocnePole;
  int *pracovniPole;
  pracovniPole= *pole;
  pomocnePole = new int[nn];
  //vynulujeme pomocne pole
  for( k=0; k<nn; k++){
    pomocnePole[k] = 0;
    
  }
  posl =0;
  
  // odstranim ostatni stupne
  for( k=0; k<*velikost; k++){
    if(aktStupen != nadjnodnod[pracovniPole[k]]){
      
      aktStupen = nadjnodnod[pracovniPole[k]];
      //printf("akt stupen %d\n",aktStupen);
      pomocnePole[posl] = pracovniPole[k];
      //printf("zmena %d",hranice[k]);
      posl++;
    }
  }
  delete[] pracovniPole;
  pracovniPole = new int[posl + 1];
  for( k=0; k<posl; k++){
    pracovniPole[k]= pomocnePole[k];
  }
  *velikost = posl ;
  *pole = pracovniPole;
  
}

void gtopology::findPseudoPheripheral(int *startUzel, int *cilUzel){
  int *pracovniSeznam;
  int velikostSeznamu;
  int hloubkaSeznamu;
  int maxSirkaSeznamu;
  int start,cil ;
  int aktualniSirka = 99999999;
  int zpracSirka,zpracHloubka,zpracVelikost,*zpracVzdalenosti,*zpracLastLevel;
  bool nalezenCil = false;
  bool novyStart = false;
  zpracLastLevel = NULL; 
  start = findMinimumDegree();
  // tvorim startovni seznam
  //printf("NalezenCil %d, novyStart %d\n", nalezenCil,novyStart);
  zpracVzdalenosti = new int[nn];
  while(!((nalezenCil)&(!novyStart))){
    //printf(" %d a %d \n", nalezenCil, novyStart);
    //printf("iterace");
    //printf("Tvorim startovni seznam start %d\n", start);
    generateInitialList(start, &pracovniSeznam, &velikostSeznamu, &hloubkaSeznamu, &maxSirkaSeznamu);
    /*for(int k=0; k<velikostSeznamu; k++){
      printf("Prvek %d  \n" , pracovniSeznam[k]);
      }*/
    ///printf("velikost %d, hloubka %d, sirka %d  \n" , velikostSeznamu, hloubkaSeznamu,maxSirkaSeznamu);
      novyStart = false;
      nalezenCil = false;
      
      for(int k=0; k<velikostSeznamu; k++){
	
	///printf("NalezenCil %d, novyStart %d\n", nalezenCil,novyStart);
	///printf("Zpracovavam uzel %d\n",pracovniSeznam[k]);
	lastLevelOfRootedStructure(pracovniSeznam[k], &zpracVelikost, &zpracSirka, &zpracHloubka, zpracVzdalenosti, &zpracLastLevel)	;
	//	for(int j=0; j<zpracVelikost; j++){
	//	//printf("Zprac Prvek %d  \n" , zpracLastLevel[j]);
	//}
	//printf("velikost %d, hloubka %d, sirka %d  \n" , zpracVelikost, zpracHloubka,zpracSirka);
	if((zpracHloubka > hloubkaSeznamu) && (zpracSirka < aktualniSirka)){
	  // nastavuji novy pocatecni uzel
	  //printf("novy start\n");
	  start = pracovniSeznam[k];
	  novyStart = true;
	  
	  delete[] pracovniSeznam;
	  break;
	  
	}
	else if(zpracSirka < aktualniSirka){
	  // nastavuju koncovy uzel
	  cil = pracovniSeznam[k];
	  //printf("novy cil\n");
	  aktualniSirka = zpracSirka;
	  nalezenCil = true;
	  break;
	}
	
	
      }
      
  }
  *startUzel = start;
  *cilUzel = cil;
  
}

void gtopology::sloan_renumb (FILE *out)
{
  adjacnodes (out);
  
  int zpracSirka,zpracHloubka,zpracVelikost,*zpracLastLevel;
  //int *zpracVzdalenosti;
  zpracLastLevel = NULL;
  int start,cil;
  int oznacenych=-1;
  int *vzdalenosti, *priority, *statusy, *priorFronta,*noveCislovani;
  int vrcholFronty;
  int zpracUzel;
  int maxPriorita;
  int pomUzel,pomActive, poziceUzlu;
  int j,k;
  //1.entry
  //printf("ACTIVE = %d, 	INACTIVE = %d, 	PREACTIVE=%d,	POSTACTIVE = %d\n",ACTIVE, 	INACTIVE , 	PREACTIVE,	POSTACTIVE);
  findPseudoPheripheral(&start, &cil);
  //start = 0;
  
  //printf("Cislo startovniho uzlu je %d \n souradnice [ %e, %e, %e ]", start, gnodes[start].x, gnodes[start].y, gnodes[start].z);
  //cil = 3;
  //printf("\nCislo ciloveho uzlu je %d \n souradnice [ %e, %e, %e ]", cil, gnodes[cil].x, gnodes[cil].y, gnodes[cil].z);
  //printf("\nstart %d cil %d \n",start,cil);
  
  vzdalenosti = new int[nn];
  priority = new int[nn];
  statusy = new int[nn];
  priorFronta = new int[nn];
  noveCislovani = new int[nn];
  int W1 = 1, W2 = 2;
  //2.compute distances
  lastLevelOfRootedStructure(cil, &zpracVelikost, &zpracSirka, &zpracHloubka, vzdalenosti, &zpracLastLevel)	;
  //3.assign initial status and priority
  
  for(k=0; k<nn; k++){
    //printf(" uzel %d vzdalenost k cily %d a pocet sousedu %ld ", k, vzdalenosti[k]+1, nadjnodnod[k]-1);
    
    for (j = 0; j< (nadjnodnod[k]);j++){
      
      //printf(" %ld", adjnodnod[k][j]);
    }
    //printf("\n");
    priority[k] = W1*(vzdalenosti[k]+1) - W2*(nadjnodnod[k]-1);
    statusy[k] = INACTIVE;
    priorFronta[k] = 0;
    
    /*printf(" vzdalenost uzlu %d od uzlu %d je %d a priorita je %d\n", k,cil,vzdalenosti[k],priority[k] );*/
  }
  //printf("za cyklem\n\n\n");
  //writePriorityQueue(priority);
  //4.initizlize node count and priority queue
  oznacenych=0;
  statusy[start] = PREACTIVE;
  vrcholFronty = 0;
  priorFronta[vrcholFronty] = start;
  // 5.test for termination
  int f=0;
  
  
  while(vrcholFronty > -1){
    
    //writePriorityQueue(priorFronta);
    //printf("\n");
    
    for( k=0; k<nn; k++){
      //printf(" %d",statusy[k]);
      //if( ( (k+1)%3)  == 0)printf("\n");
    }
    //getchar();
    //printf("Nova iterace vrchol fronty je %d\n\n\n",vrcholFronty); 
    f=0;
    while(f<=vrcholFronty){
      //printf("Prioritni fronta pozice %d je prvek %d\n",f,priorFronta[f]);
      f++;
    }
    
    //hledam maximalni prioritu
    //6. select node to be labeled
    maxPriorita = -9999999;
    
    for(k=0; k<=vrcholFronty; k++){// opravit pokud je vrcholFronty +1 mimo meze spadne to!!!
      //printf("Prioritni fronta pozice %d je uzel %d\n",k,priorFronta[k]);
      
      if(priority[priorFronta[k]] > maxPriorita){
	maxPriorita = priority[priorFronta[k]];
	zpracUzel = priorFronta[k];
	poziceUzlu = k;
	//
	
      }
      //printf("Maximalni prioritu ma uzel na pozici %d uzel %d status %d\n",k,  priorFronta[k],statusy[zpracUzel]);
      //printf("Budu zpracovavat uzel %d\n", zpracUzel);
    }
    
    //smazu zpracUzel z fronty
    //7. update queue and priorities
    priorFronta[poziceUzlu] = priorFronta[vrcholFronty];
    vrcholFronty--;
    if(statusy[zpracUzel] == PREACTIVE){
      //printf("Je preactive\n");
      for(j = 0; j < nadjnodnod[zpracUzel]; j++){
	
	pomUzel = adjnodnod[zpracUzel][j];
	if(pomUzel == zpracUzel)continue;
	//printf("Prohledavam okoli %d uzel %d se statusem %d \n",zpracUzel,pomUzel,statusy[pomUzel] );
	priority[pomUzel] = priority[pomUzel] + W2;
	if(statusy[pomUzel] == INACTIVE){
	  
	  vrcholFronty++;
	  priorFronta[vrcholFronty] = pomUzel;
	  statusy[pomUzel] = PREACTIVE;
	  //printf("vkladam do fronty na pozici %d uzel %d\n",vrcholFronty,pomUzel);
	  //printf("davam status %d\n",PREACTIVE);
	}
      }
      
    }
    
    // 8.label the next node
    //printf("oznacuji uzel %d\n",zpracUzel);
    oznacenych++;
    noveCislovani[zpracUzel] = oznacenych;
    statusy[zpracUzel] = POSTACTIVE;
    //9. update queue and priorities
    for(j = 0; j < nadjnodnod[zpracUzel]; j++){
      
      //printf("pred cyklem\n\n\n");
      pomUzel = adjnodnod[zpracUzel][j];
      //printf("pred cyklem\n\n\n");
      if(pomUzel == zpracUzel)continue;
      if(statusy[pomUzel] == PREACTIVE){
	
	
	//printf("Nastavuji uzel %d na ACTIVE\n",pomUzel);
	statusy[pomUzel] = ACTIVE;
	priority[pomUzel] = priority[pomUzel] + W2;
	//printf("Prochazim okoli uzlu %d\n",pomUzel);
	
	for(k = 0; k < nadjnodnod[pomUzel]; k++){
	  
	  pomActive = adjnodnod[pomUzel][k];
	  if(pomUzel == pomActive)continue;
	  //printf("Prohledavam krok 9. okoli %d uzel %d se statusem %d \n",pomUzel,pomActive,statusy[pomActive] );
	  if(statusy[pomActive] != POSTACTIVE){
	    priority[pomActive] = priority[pomActive] + W2;	
	  }
	  if(statusy[pomActive] == INACTIVE){
	    
	    vrcholFronty++;
	    priorFronta[vrcholFronty] = pomActive;
	    statusy[pomActive] = PREACTIVE;
	    //printf("vkladam do fronty na pozici %d uzel %d\n",vrcholFronty,pomActive);
	    //printf("davam status %d\n",PREACTIVE);
	  }
	}
      }
      
    }
  }

  //FILE *outt;
  //outt=fopen ("sloan.txt","w");
  for( int k=0; k<nn; k++){
    ordering[noveCislovani[k]-1] = k;
    //printf(" Puvodni %d nove %d\n", k, noveCislovani[k]);
    //fprintf(outt," Puvodni %d nove %d\n", k, noveCislovani[k]);
  }
  //fclose (outt);
  
  //printf("--------------------------");
  //for( int k=0; k<nn; k++){
  //	
  //	printf(" Puvodni %d nove %d\n", k, ordering[k] );
  //}
  //printf("--------------------------");
  //printf("--------------------------");
}










/**
  The function searches for the hanging nodes, initializes
  master node arrays on elements and computes the number of DOFs on elements with respect to master nodes.
   
  Created by JK, 30.3.2012
*/
void gtopology::searching_hanging_nodes (void)
{
  long i,j,k,nodid,nne,auxnne,ndofn,ndofe, nd=0;
  ivector nod;
  ivector nmn(nn), nmne(nn), nmnh(nn), mne(ne);
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    //  number nodes on element
    nne = give_original_nne (i);
    reallocv (RSTCKIVEC(nne,nod));
    //  nodes on element
    give_original_nodes (i,nod);
    
    auxnne=0;
    //  loop over the number nodes on element
    for (j=0;j<nne;j++){
      //  the number of DOFs on the j-th node
      ndofn = gnodes[nod[j]].ndofn;
      if (ndofn<0){
	//  this is hanging node
	//  |ndofn| denotes the number of master nodes
	auxnne-=ndofn;
      }else{
	auxnne++;
      }
    }//  end of loop over the number nodes on element
    
    if (auxnne>nne){
      //  the element contains hanging node
      
      //  indicator of hanging nodes in the problem
      hangnod=1;
      
      //  nmne - the number of master nodes on element
      //  array master_nodes and the variable nmne are determined
      gelements[i].master_nodes = new long [auxnne];
      gelements[i].nmne=auxnne;
      
      auxnne=0;
      ndofe=0;
      for (j=0;j<nne;j++){
	nodid=nod[j];
	ndofn = gnodes[nodid].ndofn;
	if (ndofn>=0){
	  //  the classical node
	  gelements[i].master_nodes[auxnne]=nodid;
	  ndofe+=give_ndofn (nodid);
	  auxnne++;
	}
	if (ndofn<0){
	  //  the hanging node
	  for (k=0;k<0-ndofn;k++){
	    gelements[i].master_nodes[auxnne]=gnodes[nodid].mnodes[k];
            nd = give_ndofn (gnodes[nodid].mnodes[k]);
            if (nd < 0){
              nmn(gnodes[nodid].mnodes[k]) = 1;
              nmne(gnodes[nodid].mnodes[k]) = i;
              nmnh(gnodes[nodid].mnodes[k]) = nodid;
              mne(i) = 1;
              print_warning("master node %ld for hanging node %ld is also hanging node (elem %ld)",
                            __FILE__, __LINE__, __func__, gnodes[nodid].mnodes[k]+1, nodid+1, i+1);              
            }
	    ndofe+= nd;
	    auxnne++;
	  }
	}
      }//  end of loop over the number nodes on element
      if (ndofe < 0)
        ndofe = -ndofe;
      //  new number of DOFs on element is determined
      gelements[i].ndofe=ndofe;
      //  allocation of the array cn on element
      gelements[i].cn = new long [ndofe];
    }//  end of the if statement
  }//  end of the loop over the number of elements
  if (nd){
    FILE *log = fopen("hn-test.log", "wt");
    fprintf(log, "Master nodes with negative ndofn:\n"
            "Master_id,  slave_id,  elem_id\n");
    for(i=0; i<nn; i++){
      if (nmn(i)){
        fprintf(log, "%6ld %6ld %6ld\n", i+1, nmnh(i)+1, nmne(i)+1);
      }
    }
    fprintf(log, "List of elements with negative ndofn on master nodes\n");
    for(i=0; i<ne; i++){
      if (mne(i)){
        fprintf(log, "%6ld\n", i+1);
      }
    }
    fclose(log);
  }
}



/**
   The function transfers DOFs from the master nodes to elements
   because hanging nodes do not have own DOFs.
   
   JK, 17. 6. 2012
*/
void gtopology::dof_transfer_hangnod (void)
{
  long i,j,k,l,nmne,ndofn;
  ivector nod;
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    if (gelements[i].master_nodes!=NULL){
      //  this is an element with hanging nodes
      
      gelements[i].cne=1;

      //  the number of all nodes on the element
      nmne = give_nmne (i);
      
      reallocv(RSTCKIVEC(nmne,nod));
      //  element nodes
      give_master_nodes (i,nod);
      
      l=0;
      //  loop over the number of element nodes
      for (j=0;j<nmne;j++){
	//  the number of DOFs in node
	ndofn = give_ndofn (nod[j]);
	
	//  loop over the number of DOFs
	for (k=0;k<ndofn;k++){
	  gelements[i].cn[l]=give_dof (nod[j],k);
	  l++;
	}//  end of the loop over the number of DOFs
      }//  end of the loop over the number of element nodes
    }
  }//  end of the loop over the number of elements
}



/**
   function assembles weights for approximation on a given entity

   it works for lines, triangles, quadrilaterals,
   tetrahedrons and hexahedrons
   
   it can be used in connections with hanging nodes
   
   @param et[in] - type of entity
   @param nid[in] - node id
   @param[out] w - array of weights
   
   JK, 11. 7. 2012
*/
void gtopology::approx_weights (gtypel et, long nid, vector &w)
{
  double xi,eta,zeta;
  
  switch (et){
  case isolinear1d:{
    if (w.n != 2){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    bf_lin_1d (w.a,xi);
    break;
  }
  case isoquadratic1d:{
    if (w.n != 3){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    bf_quad_1d (w.a,xi);
    break;
  }
  case trianglelinear:{
    if (w.n != 3){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    bf_lin_3_2d (w.a,xi,eta);
    break;
  }
  case trianglequadratic:{
    if (w.n != 6){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    bf_quad_3_2d (w.a,xi,eta);
    break;
  }
  case isolinear2d:{
    if (w.n != 4){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    bf_lin_4_2d (w.a,xi,eta);
    break;
  }
  case isoquadratic2d:{
    if (w.n != 8){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    bf_quad_4_2d (w.a,xi,eta);
    break;
  }
  case tetrahedronlinear:{
    if (w.n != 4){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    zeta = gnodes[nid].natcoord[2];
    bf_lin_tet (w.a,xi,eta,zeta);
    break;
  }
  case tetrahedronquadratic:{
    if (w.n != 10){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    zeta = gnodes[nid].natcoord[2];
    bf_quad_tet (w.a,xi,eta,zeta);
    break;
  }
  case isolinear3d:{
    if (w.n != 8){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    zeta = gnodes[nid].natcoord[2];
    bf_lin_hex_3d (w.a,xi,eta,zeta);
    break;
  }
  case isoquadratic3d:{
    if (w.n != 20){
      print_err("wrong number of components of the vector of weights", __FILE__, __LINE__, __func__);
    }
    xi = gnodes[nid].natcoord[0];
    eta = gnodes[nid].natcoord[1];
    zeta = gnodes[nid].natcoord[2];
    bf_quad_hex_3d (w.a,xi,eta,zeta);
    break;
  }
  default:{
    print_err("unknown type of entity is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function checks natural and real coordinates of hanging nodes
   
   @param out - output file
   
   JK, 16. 11. 2012
*/
void gtopology::hang_nodes_check (void)
{
  long i,j,ndofn,mnid;
  double xm,ym,zm,zero;
  vector w,x,y,z;
  
  zero=1.0e-2;
  
  //  loop over the number of nodes
  for (i=0;i<nn;i++){
    //  the number of DOFs or master nodes
    ndofn = gnodes[i].give_ndofn ();
    if (ndofn<0){
      //  the actual node is hanging node
      
      ndofn=0-ndofn;
      reallocv (RSTCKVEC(ndofn,w));
      reallocv (RSTCKVEC(ndofn,x));
      reallocv (RSTCKVEC(ndofn,y));
      reallocv (RSTCKVEC(ndofn,z));
      
      //  weights of contributions from the master nodes
      approx_weights (gnodes[i].masentity,i,w);
      
      //  loop over the number of master nodes
      for (j=0;j<ndofn;j++){
	//  id of the master node
	mnid = gnodes[i].mnodes[j];
	
	//  coordinates of the master nodes
	x[j] = gnodes[mnid].x;
	y[j] = gnodes[mnid].y;
	z[j] = gnodes[mnid].z;
      }//  end of the loop over the number of master nodes
      
      //  real coordinates of the actual hanging node obtained from the master nodes
      scprd(x,w,xm);
      scprd(y,w,ym);
      scprd(z,w,zm);
      
      if (fabs(gnodes[i].x-xm)>zero){
	print_err("computed x-coordinate of the hanging node %ld is different from coordinate from input file (%le x %le)",
		  __FILE__,__LINE__,__func__,i+1,xm,gnodes[i].x);
      }
      if (fabs(gnodes[i].y-ym)>zero){
	print_err("computed y-coordinate of the hanging node %ld is different from coordinate from input file (%le x %le)",
		  __FILE__,__LINE__,__func__,i+1,ym,gnodes[i].y);
      }
      if (fabs(gnodes[i].z-zm)>zero){
	print_err("computed z-coordinate of the hanging node %ld is different from coordinate from input file (%le x %le)",
		  __FILE__,__LINE__,__func__,i+1,zm,gnodes[i].z);
      }
    }//  end of the statement if (ndofn<0){
  }//  end of the loop over the number of nodes
}



/**
  The function adds code numbers used in stress based approach in homogenization.
  The code numbers are copied from nodes to elements and additional code numbers
   are added after the regular code numbers. Dimension of adcn vector is the number 
   of stress/strain components (6 for 3D problems, 3 for plane stress/strain)
   
   @param cnadd[in] - %vector of additional code (DOF) numbers to be added to regular element DOF numbers from nodes
                      cnadd[i] > 0 <=> i-th macro-value component is represented by the given code (DOF) number
                      cnadd[i] = 0 <=> i-th macro-value component constitues no additional DOF
   @param dof[in/out] - actual total number of DOFs in the problem, it will be increased according to new DOFs 
                        defined in cnadd
   
   JK, 17. 2. 2014
*/
void gtopology::code_num_mod (ivector &cnadd, long &dof)
{
  long i,j,k,ndofe,ncomp = cnadd.n;
  ivector cn;

  for (i=0; i<ncomp; i++){
    if(cnadd[i])   dof++;
  }
  
  for (i=0;i<ne;i++){
    ndofe = give_ndofe(i)-ncomp;
    reallocv (RSTCKIVEC(ndofe,cn));
    give_code_numbers (i,cn.a);
    
    gelements[i].cne=1;
    gelements[i].cn = new long [ndofe+ncomp];
    for (j=0; j<ndofe; j++){
      gelements[i].cn[j] = cn.a[j];
    }
    k=0;
    for (j=ndofe; j<ndofe+ncomp; j++){
      gelements[i].cn[j] = cnadd[k];
      k++;
    }
  }
}



/**
   function constructs centers of finite elements and radii of circumscribed balls
   
   this function is used in connection with determination of hanging nodes or
   if the characteristic sizes of finite elements are required
   
   JK, 24. 9. 2016
*/
void gtopology::construct_circumscribed_balls ()
{
  long i;
  
  allocate_arrays_circumscribed_ball ();
  
  for (i=0;i<ne;i++){
    element_center (i);
    radius_circumscribed_ball (i);
  }
  
}

/**
   function returns characteristic sizes of finite elements
   it is diameter of circumscribed balls
   
   @param eid - element id
   
   JK, 24. 9. 2016
*/
double gtopology::give_characteristic_size (long eid)
{
  return 2.0*circumrad[eid];
}

/**
   function computes centers of elements
   this function is used in connection with determination of hanging nodes
   
   @param eid - element id
   
   JK, 23. 8. 2016
*/
void gtopology::element_center (long eid)
{
  long i,nne;
  double xx,yy,zz;
  
  //  the number of nodes on element
  nne = give_nne (eid);
  //  vectors of nodal coordinates
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  coordinates of nodes on element
  give_node_coord3d (x,y,z,eid);
  
  //  coordinates of the element center
  xx=0.0;  yy=0.0;  zz=0.0;
  for (i=0;i<nne;i++){
    xx += x[i];
    yy += y[i];
    zz += z[i];
  }
  
  xc[eid]=xx/nne;
  yc[eid]=yy/nne;
  zc[eid]=zz/nne;

}

/**
   function computes radius of ball circumscribed to the element
   this function is used in connection with determination of hanging nodes
   function gtopology::element_center (long eid) has to called before this one

   @param eid - element id
   
   JK, 23. 8. 2016
*/
void gtopology::radius_circumscribed_ball (long eid)
{
  long i,nne;
  double xx,yy,zz,d,r;
  
  //  the number of nodes on element
  nne = give_nne (eid);
  //  vectors of nodal coordinates
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  //  coordinates of nodes on element
  give_node_coord3d (x,y,z,eid);
  
  //  coordinates of the element center
  xx=xc[eid];
  yy=yc[eid];
  zz=zc[eid];
  
  r=0.0;
  for (i=0;i<nne;i++){
    d = sqrt((x[i]-xx)*(x[i]-xx) + (y[i]-yy)*(y[i]-yy) + (z[i]-zz)*(z[i]-zz));
    if (r<d)
      r=d;
  }
  
  circumrad[eid]=r;

}

/**
   function allocates arrays necessary for construction of balls circumscribed to elements
   
   JK, 23. 8. 2016
*/
void gtopology::allocate_arrays_circumscribed_ball ()
{
  if (xc != NULL)
    delete [] xc;
  else
    xc = new double [ne];

  if (yc != NULL)
    delete [] yc;
  else
    yc = new double [ne];

  if (zc != NULL)
    delete [] zc;
  else
    zc = new double [ne];

  if (circumrad != NULL)
    delete [] circumrad;
  else
    circumrad = new double [ne];

}

/**
   function
   
   @param eid - element id
   @param xb,yb,zb - vectors of bar nodal coordinates
   @param xh,yh,zh - vectors of hexahedron nodal coordinates
   @param zero - computer zero
   @param xi,yi,zi - vectors of coordinates of intersection points
   @param xxi,yyi,zzi - vectors of natural coordinates of intersection points
   @param nod - %list of master nodes for each hanging node

   JK, 23. 8. 2016
*/
void gtopology::bar_linhex_intersection (long eid,long &numinter,vector &xb,vector &yb,vector &zb,vector &xh,vector &yh,vector &zh,double zero,
                                         vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  ivector nod(ASTCKIVEC(8));
  //  vectors of coefficients 
  vector v1(ASTCKVEC(3)),v2(ASTCKVEC(3)),v3(ASTCKVEC(3)),v4(ASTCKVEC(3)),v5(ASTCKVEC(3));
  
  //  nodes on element
  give_nodes (eid,nod);
  
  // ***************
  //  surface n. 1
  // ***************
  v1[0]=xb[0]-(xh[0]+xh[3]+xh[4]+xh[7])/4.0;
  v1[1]=yb[0]-(yh[0]+yh[3]+yh[4]+yh[7])/4.0;
  v1[2]=zb[0]-(zh[0]+zh[3]+zh[4]+zh[7])/4.0;
  
  v2[0]=(xh[0]-xh[3]+xh[4]-xh[7])/4.0;
  v2[1]=(yh[0]-yh[3]+yh[4]-yh[7])/4.0;
  v2[2]=(zh[0]-zh[3]+zh[4]-zh[7])/4.0;

  v3[0]=(xh[0]+xh[3]-xh[4]-xh[7])/4.0;
  v3[1]=(yh[0]+yh[3]-yh[4]-yh[7])/4.0;
  v3[2]=(zh[0]+zh[3]-zh[4]-zh[7])/4.0;

  v4[0]=(xh[0]-xh[3]-xh[4]+xh[7])/4.0;
  v4[1]=(yh[0]-yh[3]-yh[4]+yh[7])/4.0;
  v4[2]=(zh[0]-zh[3]-zh[4]+zh[7])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,1,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

  // ***************
  //  surface n. 2
  // ***************
  v1[0]=xb[0]-(xh[0]+xh[1]+xh[4]+xh[5])/4.0;
  v1[1]=yb[0]-(yh[0]+yh[1]+yh[4]+yh[5])/4.0;
  v1[2]=zb[0]-(zh[0]+zh[1]+zh[4]+zh[5])/4.0;
  
  v2[0]=(xh[0]-xh[1]+xh[4]-xh[5])/4.0;
  v2[1]=(yh[0]-yh[1]+yh[4]-yh[5])/4.0;
  v2[2]=(zh[0]-zh[1]+zh[4]-zh[5])/4.0;

  v3[0]=(xh[0]+xh[1]-xh[4]-xh[5])/4.0;
  v3[1]=(yh[0]+yh[1]-yh[4]-yh[5])/4.0;
  v3[2]=(zh[0]+zh[1]-zh[4]-zh[5])/4.0;

  v4[0]=(xh[0]-xh[1]-xh[4]+xh[5])/4.0;
  v4[1]=(yh[0]-yh[1]-yh[4]+yh[5])/4.0;
  v4[2]=(zh[0]-zh[1]-zh[4]+zh[5])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,2,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

  // ***************
  //  surface n. 3
  // ***************
  v1[0]=xb[0]-(xh[1]+xh[2]+xh[5]+xh[6])/4.0;
  v1[1]=yb[0]-(yh[1]+yh[2]+yh[5]+yh[6])/4.0;
  v1[2]=zb[0]-(zh[1]+zh[2]+zh[5]+zh[6])/4.0;
  
  v2[0]=(xh[1]-xh[2]+xh[5]-xh[6])/4.0;
  v2[1]=(yh[1]-yh[2]+yh[5]-yh[6])/4.0;
  v2[2]=(zh[1]-zh[2]+zh[5]-zh[6])/4.0;

  v3[0]=(xh[1]+xh[2]-xh[5]-xh[6])/4.0;
  v3[1]=(yh[1]+yh[2]-yh[5]-yh[6])/4.0;
  v3[2]=(zh[1]+zh[2]-zh[5]-zh[6])/4.0;

  v4[0]=(xh[1]-xh[2]-xh[5]+xh[6])/4.0;
  v4[1]=(yh[1]-yh[2]-yh[5]+yh[6])/4.0;
  v4[2]=(zh[1]-zh[2]-zh[5]+zh[6])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,3,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

  // ***************
  //  surface n. 4
  // ***************
  v1[0]=xb[0]-(xh[2]+xh[3]+xh[6]+xh[7])/4.0;
  v1[1]=yb[0]-(yh[2]+yh[3]+yh[6]+yh[7])/4.0;
  v1[2]=zb[0]-(zh[2]+zh[3]+zh[6]+zh[7])/4.0;
  
  v2[0]=(0.0-xh[2]+xh[3]-xh[6]+xh[7])/4.0;
  v2[1]=(0.0-yh[2]+yh[3]-yh[6]+yh[7])/4.0;
  v2[2]=(0.0-zh[2]+zh[3]-zh[6]+zh[7])/4.0;

  v3[0]=(xh[2]+xh[3]-xh[6]-xh[7])/4.0;
  v3[1]=(yh[2]+yh[3]-yh[6]-yh[7])/4.0;
  v3[2]=(zh[2]+zh[3]-zh[6]-zh[7])/4.0;

  v4[0]=(0.0-xh[2]+xh[3]+xh[6]-xh[7])/4.0;
  v4[1]=(0.0-yh[2]+yh[3]+yh[6]-yh[7])/4.0;
  v4[2]=(0.0-zh[2]+zh[3]+zh[6]-zh[7])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,4,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

  // ***************
  //  surface n. 5
  // ***************
  v1[0]=xb[0]-(xh[0]+xh[1]+xh[2]+xh[3])/4.0;
  v1[1]=yb[0]-(yh[0]+yh[1]+yh[2]+yh[3])/4.0;
  v1[2]=zb[0]-(zh[0]+zh[1]+zh[2]+zh[3])/4.0;
  
  v2[0]=(xh[0]-xh[1]-xh[2]+xh[3])/4.0;
  v2[1]=(yh[0]-yh[1]-yh[2]+yh[3])/4.0;
  v2[2]=(zh[0]-zh[1]-zh[2]+zh[3])/4.0;

  v3[0]=(xh[0]+xh[1]-xh[2]-xh[3])/4.0;
  v3[1]=(yh[0]+yh[1]-yh[2]-yh[3])/4.0;
  v3[2]=(zh[0]+zh[1]-zh[2]-zh[3])/4.0;

  v4[0]=(xh[0]-xh[1]+xh[2]-xh[3])/4.0;
  v4[1]=(yh[0]-yh[1]+yh[2]-yh[3])/4.0;
  v4[2]=(zh[0]-zh[1]+zh[2]-zh[3])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,5,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

  // ***************
  //  surface n. 6
  // ***************
  v1[0]=xb[0]-(xh[4]+xh[5]+xh[6]+xh[7])/4.0;
  v1[1]=yb[0]-(yh[4]+yh[5]+yh[6]+yh[7])/4.0;
  v1[2]=zb[0]-(zh[4]+zh[5]+zh[6]+zh[7])/4.0;
  
  v2[0]=(xh[4]-xh[5]-xh[6]+xh[7])/4.0;
  v2[1]=(yh[4]-yh[5]-yh[6]+yh[7])/4.0;
  v2[2]=(zh[4]-zh[5]-zh[6]+zh[7])/4.0;

  v3[0]=(xh[4]+xh[5]-xh[6]-xh[7])/4.0;
  v3[1]=(yh[4]+yh[5]-yh[6]-yh[7])/4.0;
  v3[2]=(zh[4]+zh[5]-zh[6]-zh[7])/4.0;

  v4[0]=(xh[4]-xh[5]+xh[6]-xh[7])/4.0;
  v4[1]=(yh[4]-yh[5]+yh[6]-yh[7])/4.0;
  v4[2]=(zh[4]-zh[5]+zh[6]-zh[7])/4.0;
  
  v5[0]=xb[0]-xb[1];
  v5[1]=yb[0]-yb[1];
  v5[2]=zb[0]-zb[1];

  bar_quadrangle_surface_intersection (numinter,6,nod,xb,yb,zb,xh,yh,zh,v1,v2,v3,v4,v5,zero,xi,yi,zi,xxi,yyi,zzi,masnod);

}


/**
   function
   
   @param numinter - the number of intersection points
   @param sid - surface id
   @param xb,yb,zb - vectors of nodal coordinates of bar
   @param xh,yh,zh - vectors of nodal coordinates of hexahedron
   @param v1,v2,v3,v4,v5 - vectors of coefficients of functions describing surface
   @param zero - computer zero/threshold
   @param xi,yi,zi - coordinates of intersection points
   @param xxi,yyi,zzi - natural coordinates of intersection points
   @param masnod - list of master nodes for each hanging node

   JK, 23. 8. 2016
*/
void gtopology::bar_quadrangle_surface_intersection (long &numinter,long sid,ivector &nod,vector &xb,vector &yb,vector &zb,vector &xh,vector &yh,vector &zh,
                                                     vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,double zero,
                                                     vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  long stop,stopnew;
  double nor,sum;
  double alpha, beta;
  vector n(ASTCKVEC(8));
 
  //  norm of the vector v4
  nor = v4[0]*v4[0] + v4[1]*v4[1] + v4[2]*v4[2];
  
  if (sqrt(nor)<zero){
    // ***********************************************************
    //  the vector v4 is zero vector
    //  it means, the surface is planar and the problem is linear
    // ***********************************************************
    vector x(ASTCKVEC(3)),y(ASTCKVEC(3));
    matrix a(ASTCKMAT(3,3));

    a[0][0] = v2[0];  a[0][1] = v3[0];  a[0][2] = v5[0];
    a[1][0] = v2[1];  a[1][1] = v3[1];  a[1][2] = v5[1];
    a[2][0] = v2[2];  a[2][1] = v3[2];  a[2][2] = v5[2];
    
    y[0] = v1[0];
    y[1] = v1[1];
    y[2] = v1[2];
    
    stop = gemp (a,x,y,zero,1);
    if (stop==0){
      //  matrix is nonsingular
      //  there is an unique solution
      if (-1.0 <= x[0] && x[0] <= 1.0){
        if (-1.0 <= x[1] && x[1] <= 1.0){
          if (0.0 <= x[2] && x[2] <= 1.0){
            //  there is intersection point of the bar with the surface on element
  
            surface_linhex_natcoordinates (sid,x[0],x[1],zero,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);

          }
        }
      }//  end of the statement if (-1.0 <= x[0] && x[0] <= 1.0){
      else{
        //  the intersection point is out of the element
      }
    }//  end of the statement if (stop==0){
    else{
      //  matrix is singular
      //  there is no solution or infinite number of solutions
      sum = 0.0;
      sum = fabs(a[2][0]) + fabs(a[2][1]) + fabs(a[2][2]);
      if (sum < zero && fabs(y[2]) > zero){
        //  there is no solution of the system
      }else{
        //  there is infinite number of solutions
        //  the bar lies in the plane
        //  it is necessary to check whether it intersect edges
        switch (sid){
        case 1:{
          bar_edge_intersection (numinter, 4,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,12,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 8,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 9,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        case 2:{
          bar_edge_intersection (numinter, 1,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 9,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 5,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,10,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        case 3:{
          bar_edge_intersection (numinter, 2,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,11,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 6,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,10,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        case 4:{
          bar_edge_intersection (numinter, 3,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,12,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,11,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter, 7,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        case 5:{
          bar_edge_intersection (numinter,1,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,2,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,3,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,4,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        case 6:{
          bar_edge_intersection (numinter,5,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,6,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,7,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          bar_edge_intersection (numinter,8,nod,xb,yb,zb,xh,yh,zh,zero,xi,yi,zi,xxi,yyi,zzi,masnod);
          break;
        }
        default:{

        }
        }
      }
    }
  }//  end of the statement if (sqrt(nor)<Mp->zero){
  else{
    // ********************************************************************
    //  the vector v4 is non-zero
    //  it means, the surface is non-planar and the problem is non-linear
    // ********************************************************************
    matrix jac(ASTCKMAT(3,3));
    
    alpha=0.0;  beta=0.0;
    stopnew=newton_intersection (zero,alpha,beta,v1,v2,v3,v4,v5,jac);
    if (stopnew==0){
      alpha=1.0;  beta=1.0;
      stopnew=newton_intersection (zero,alpha,beta,v1,v2,v3,v4,v5,jac);
      if (stopnew==0){
	alpha=-1.0;  beta=1.0;
	stopnew=newton_intersection (zero,alpha,beta,v1,v2,v3,v4,v5,jac);
	if (stopnew==0){
	  alpha=-1.0;  beta=-1.0;
	  stopnew=newton_intersection (zero,alpha,beta,v1,v2,v3,v4,v5,jac);
	  if (stopnew==0){
	    alpha=1.0;  beta=-1.0;
	    stopnew=newton_intersection (zero,alpha,beta,v1,v2,v3,v4,v5,jac);
	  }
	}
      }
    }//  end of the statement if (stopnew==0){
    
    if (stopnew==1){
      surface_linhex_natcoordinates (sid,alpha,beta,zero,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);
    }else{
      
    }
    
  }
  
}

void gtopology::bar_linhex_intersec_jacobian (double alpha,double beta,vector &v2,vector &v3,vector &v4,vector &v5,matrix &jac)
{
  jac[0][0] = v2[0]+v4[0]*beta;  jac[0][1] = v3[0]+v4[0]*alpha;  jac[0][2] = v5[0];
  jac[1][0] = v2[1]+v4[1]*beta;  jac[1][1] = v3[1]+v4[1]*alpha;  jac[1][2] = v5[1];
  jac[2][0] = v2[2]+v4[2]*beta;  jac[2][1] = v3[2]+v4[2]*alpha;  jac[2][2] = v5[2];
}

void gtopology::bar_linhex_intersec_function (double alpha,double beta,double gamma,vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,vector &f)
{
  f[0]=0.0-v1[0]+v2[0]*alpha+v3[0]*beta+v4[0]*alpha*beta+v5[0]*gamma;
  f[1]=0.0-v1[1]+v2[1]*alpha+v3[1]*beta+v4[1]*alpha*beta+v5[1]*gamma;
  f[2]=0.0-v1[2]+v2[2]*alpha+v3[2]*beta+v4[2]*alpha*beta+v5[2]*gamma;
}

long gtopology::newton_intersection (double zero,double &alpha,double &beta,vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,matrix &jac)
{
  long i,ni=20,stopsys,stopnew=0;
  double nor,res=1.0e-8;
  vector x(ASTCKVEC(3)),y(ASTCKVEC(3)),xx(ASTCKVEC(3));
  
  x[0]=alpha;
  x[1]=beta;
  x[2]=0.0;
  
  for (i=0;i<ni;i++){
    bar_linhex_intersec_function (x[0],x[1],x[2],v1,v2,v3,v4,v5,y);
    bar_linhex_intersec_jacobian (x[0],x[1],v2,v3,v4,v5,jac);
    stopsys = gemp (jac,xx,y,zero,1);
    if (stopsys==1){
      break;
    }
    nor = normv (xx);
    if (nor<res){
      stopnew=1;
      break;
    }
    addmultv(x,xx,-1.0);
  }
  
  if (stopnew==1){
    if (-1.0 <= x[0] && x[0] <= 1.0){
      if (-1.0 <= x[1] && x[1] <= 1.0){
	if (0.0 <= x[2] && x[2] <= 1.0){
	  //  intersection of the bar with surface is on element
	  alpha=x[0];
	  beta=x[1];
	  return 1;
	}
      }
    }
    //  intersection of the bar with surface is out of the element
    return 0;
  }else{
    //  no intersection has been found
    return 0;
  }
}

/**
   function determines intersection of a bar with an edge
   
   @param numinter - the number of intersection points
   @param xb,yb,zb - vectors of nodal coordinates of bar
   @param xh,yh,zh - vectors of hexahedron nodes
   @param zero - computer zero/threshold
   @param xi,yi,zi - coordinates of intersection points
   @param xxi,yyi,zzi - natural coordinates of intersection points
   @param masnod - list of master nodes for each hanging node

   JK, 24. 8. 2016
*/
void gtopology::bar_edge_intersection (long &numinter,long edgeid,ivector &nod,vector &xb,vector &yb,vector &zb,vector &xh,vector &yh,vector &zh,double zero,
                                       vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  long stop;
  double ef,alpha,gamma,gamma2;
  vector xe(ASTCKVEC(2)),ye(ASTCKVEC(2)),ze(ASTCKVEC(2));
  vector x(ASTCKVEC(2)),y(ASTCKVEC(3));
  matrix a(ASTCKMAT(3,2));
  vector n(ASTCKVEC(8));
  
  switch (edgeid){
  case 1:{
    xe[0]=xh[1];  ye[0]=yh[1];  ze[0]=zh[1];
    xe[1]=xh[0];  ye[1]=yh[0];  ze[1]=zh[0];
    break;
  }
  case 2:{
    xe[0]=xh[2];  ye[0]=yh[2];  ze[0]=zh[2];
    xe[1]=xh[1];  ye[1]=yh[1];  ze[1]=zh[1];
    break;
  }
  case 3:{
    xe[0]=xh[2];  ye[0]=yh[2];  ze[0]=zh[2];
    xe[1]=xh[3];  ye[1]=yh[3];  ze[1]=zh[3];
    break;
  }
  case 4:{
    xe[0]=xh[3];  ye[0]=yh[3];  ze[0]=zh[3];
    xe[1]=xh[0];  ye[1]=yh[0];  ze[1]=zh[0];
    break;
  }
  case 5:{
    xe[0]=xh[5];  ye[0]=yh[5];  ze[0]=zh[5];
    xe[1]=xh[4];  ye[1]=yh[4];  ze[1]=zh[4];
    break;
  }
  case 6:{
    xe[0]=xh[6];  ye[0]=yh[6];  ze[0]=zh[6];
    xe[1]=xh[5];  ye[1]=yh[5];  ze[1]=zh[5];
    break;
  }
  case 7:{
    xe[0]=xh[6];  ye[0]=yh[6];  ze[0]=zh[6];
    xe[1]=xh[7];  ye[1]=yh[7];  ze[1]=zh[7];
    break;
  }
  case 8:{
    xe[0]=xh[7];  ye[0]=yh[7];  ze[0]=zh[7];
    xe[1]=xh[4];  ye[1]=yh[4];  ze[1]=zh[4];
    break;
  }
  case 9:{
    xe[0]=xh[4];  ye[0]=yh[4];  ze[0]=zh[4];
    xe[1]=xh[0];  ye[1]=yh[0];  ze[1]=zh[0];
    break;
  }
  case 10:{
    xe[0]=xh[5];  ye[0]=yh[5];  ze[0]=zh[5];
    xe[1]=xh[1];  ye[1]=yh[1];  ze[1]=zh[1];
    break;
  }
  case 11:{
    xe[0]=xh[6];  ye[0]=yh[6];  ze[0]=zh[6];
    xe[1]=xh[2];  ye[1]=yh[2];  ze[1]=zh[2];
    break;
  }
  case 12:{
    xe[0]=xh[7];  ye[0]=yh[7];  ze[0]=zh[7];
    xe[1]=xh[3];  ye[1]=yh[3];  ze[1]=zh[3];
    break;
  }
  default:{
    
  }
  }//  end of switch (edgeid){
  
  a[0][0] = (xe[1]-xe[0])/2.0;  a[0][1] = xb[0]-xb[1];
  a[1][0] = (ye[1]-ye[0])/2.0;  a[1][1] = yb[0]-yb[1];
  a[2][0] = (ze[1]-ze[0])/2.0;  a[2][1] = zb[0]-zb[1];
  
  y[0] = xb[0]-(xe[0]+xe[1])/2.0;
  y[1] = yb[0]-(ye[0]+ye[1])/2.0;
  y[2] = zb[0]-(ze[0]+ze[1])/2.0;
  
  if (fabs(a[0][0])<zero){
    if (fabs(a[1][0])<zero){
      //  this the following case
      //  0 x
      //  0 x
      //  x x

      if ((fabs(a[0][1])<zero && fabs(y[0])>zero) || (fabs(a[1][1])<zero && fabs(y[1])>zero)){
	//  no solution
	stop=0;
      }
      if ((fabs(a[0][1])<zero && fabs(y[0])<zero) && (fabs(a[1][1])<zero && fabs(y[1])<zero)){
	//  infinite number of solutions
	stop=8;
      }
      if (fabs(a[0][1])>zero && fabs(y[0])>zero){
	gamma = y[0]/a[0][1];
	stop=1;
      }
      if (fabs(a[1][1])>zero && fabs(y[1])>zero){
	gamma2 = y[0]/a[0][1];
	if (fabs(gamma-gamma2)>zero || stop !=1){
	  print_err("mismatch in determination of intersection", __FILE__,__LINE__,__func__);
	  abort ();
	}
      }
      
      if (stop==1){
	alpha = (y[2]-gamma*a[2][1])/a[2][0];
      }
    }
    else{
      //  this the following case
      //  0 x
      //  x x
      //  0 x
      
      ef = a[0][0]/a[1][0];
      a[0][0] = 0.0;  a[0][1] -= ef*a[1][1];
      y[0] -= ef*y[1];
      ef = a[2][0]/a[1][0];
      a[2][0] = 0.0;  a[2][1] -= ef*a[1][1];
      y[2] -= ef*y[1];
      

      if ((fabs(a[0][1])<zero && fabs(y[0])>zero) || (fabs(a[2][1])<zero && fabs(y[2])>zero)){
	//  no solution
	stop=0;
      }
      if ((fabs(a[0][1])<zero && fabs(y[0])<zero) && (fabs(a[2][1])<zero && fabs(y[2])<zero)){
	//  infinite number of solutions
	stop=8;
      }
      if (fabs(a[0][1])>zero && fabs(y[0])>zero){
	gamma = y[0]/a[0][1];
	stop=1;
      }
      if (fabs(a[2][1])>zero && fabs(y[2])>zero){
	gamma2 = y[2]/a[2][1];
	if (fabs(gamma-gamma2)>zero || stop !=1){
	  print_err("mismatch in determination of intersection", __FILE__,__LINE__,__func__);
	  abort ();
	}
      }
      
      if (stop==1){
	alpha = (y[1]-gamma*a[1][1])/a[1][0];
      }
      
    }
  }
  else{
    //  this the following case
    //  x x
    //  0 x
    //  0 x
    
    ef = a[1][0]/a[0][0];
    a[1][0] = 0.0;  a[1][1] -= ef*a[0][1];
    y[1] -= ef*y[0];
    ef = a[2][0]/a[0][0];
    a[2][0] = 0.0;  a[2][1] -= ef*a[0][1];
    y[2] -= ef*y[0];
    
    if ((fabs(a[1][1])<zero && fabs(y[1])>zero) || (fabs(a[2][1])<zero && fabs(y[2])>zero)){
      //  no solution
      stop=0;
    }
    if ((fabs(a[1][1])<zero && fabs(y[1])<zero) && (fabs(a[2][1])<zero && fabs(y[2])<zero)){
      //  infinite number of solutions
      stop=8;
    }
    if (fabs(a[1][1])>zero && fabs(y[1])>zero){
      gamma = y[1]/a[1][1];
      stop=1;
    }
    if (fabs(a[2][1])>zero && fabs(y[2])>zero){
      gamma2 = y[2]/a[2][1];
      if (fabs(gamma-gamma2)>zero || stop !=1){
	print_err("mismatch in determination of intersection", __FILE__,__LINE__,__func__);
	abort ();
      }
    }
    
    if (stop==1){
      alpha = (y[0]-gamma*a[0][1])/a[0][0];
    }
    
  }
  
  if (stop==1){
    if (0.0 <= gamma && gamma <= 1.0){
      if (-1.0 <= alpha && alpha <= 1.0){
        //  the intersection is within the surface

        edge_linhex_natcoordinates (edgeid,alpha,zero,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);
       
      }
    }
  }
  
  if (stop==8){
    //  bar coincides with edge
    a[0][0] = (xe[1]-xe[0])/2.0;
    a[1][0] = (ye[1]-ye[0])/2.0;
    a[2][0] = (ze[1]-ze[0])/2.0;
    
    y[0] = xb[0]-(xe[0]+xe[1])/2.0;
    y[1] = yb[0]-(ye[0]+ye[1])/2.0;
    y[2] = zb[0]-(ze[0]+ze[1])/2.0;
    
    if (fabs(a[0][0])>zero){
      alpha = y[0]/a[0][0];
    }
    if (fabs(a[1][0])>zero){
      alpha = y[1]/a[1][0];
    }
    if (fabs(a[2][0])>zero){
      alpha = y[2]/a[2][0];
    }
    if (-1.0 <= alpha && alpha <= 1.0){
      //  bod A je na hrane
      
      edge_linhex_natcoordinates (edgeid,alpha,zero,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);
      
    }else{
      //  bod A lezi mimo hranu, bere se uzel kostky

      end_points_edge_linhex_natcoordinates (edgeid,alpha,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);

    }
    
    y[0] = xb[1]-(xe[0]+xe[1])/2.0;
    y[1] = yb[1]-(ye[0]+ye[1])/2.0;
    y[2] = zb[1]-(ze[0]+ze[1])/2.0;
    
    if (fabs(a[0][0])>zero){
      alpha = y[0]/a[0][0];
    }
    if (fabs(a[1][0])>zero){
      alpha = y[1]/a[1][0];
    }
    if (fabs(a[2][0])>zero){
      alpha = y[2]/a[2][0];
    }
    if (-1.0 <= alpha && alpha <= 1.0){
      //  bod B je na hrane
      
      edge_linhex_natcoordinates (edgeid,alpha,zero,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);
      
    }else{
      end_points_edge_linhex_natcoordinates (edgeid,alpha,nod,xh,yh,zh,numinter,xi,yi,zi,xxi,yyi,zzi,masnod);
    }

  }

}

void gtopology::surface_linhex_natcoordinates (long sid,double alpha,double beta,double zero,ivector &nod,vector &xh,vector &yh,vector &zh,
                                               long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  long stop;
  double nc1,nc2,nc3;
  vector n(ASTCKVEC(8));
  
  switch (sid){
  case 1:{
    nc1=1.0;  nc2=alpha;  nc3=beta;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[0];
      masnod[numinter][2]=nod[3];
      masnod[numinter][3]=nod[7];
      masnod[numinter][4]=nod[4];
    }
    break;
  }
  case 2:{
    nc1=alpha;  nc2=1.0;  nc3=beta;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[0];
      masnod[numinter][2]=nod[1];
      masnod[numinter][3]=nod[5];
      masnod[numinter][4]=nod[4];
    }
    break;
  }
  case 3:{
    nc1=-1.0;  nc2=alpha;  nc3=beta;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[1];
      masnod[numinter][2]=nod[2];
      masnod[numinter][3]=nod[6];
      masnod[numinter][4]=nod[5];
    }
    break;
  }
  case 4:{
    nc1=alpha;  nc2=-1.0;  nc3=beta;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[3];
      masnod[numinter][2]=nod[2];
      masnod[numinter][3]=nod[6];
      masnod[numinter][4]=nod[7];
    }
    break;
  }
  case 5:{
    nc1=alpha;  nc2=beta;  nc3=1.0;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[0];
      masnod[numinter][2]=nod[1];
      masnod[numinter][3]=nod[2];
      masnod[numinter][4]=nod[3];
    }
    break;
  }
  case 6:{
    nc1=alpha;  nc2=beta;  nc3=-1.0;
    stop=0;
    if (fabs(alpha-1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (fabs(alpha+1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero && fabs(beta+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=4;
      masnod[numinter][1]=nod[4];
      masnod[numinter][2]=nod[5];
      masnod[numinter][3]=nod[6];
      masnod[numinter][4]=nod[7];
    }
    break;
  }
  default:{
    nc1 = nc2 = nc3 = 0.0;
    print_err("invalid surface sid=%ld.", __FILE__, __LINE__, __func__, sid);
  }
  }
  
  bf_lin_hex_3d (n.a,nc1,nc2,nc3);
  xi[numinter] = ss(n.a,xh.a,xh.n);
  yi[numinter] = ss(n.a,yh.a,yh.n);
  zi[numinter] = ss(n.a,zh.a,zh.n);
  xxi[numinter] = nc1;
  yyi[numinter] = nc2;
  zzi[numinter] = nc3;
  numinter++;
  
}


void gtopology::edge_linhex_natcoordinates (long edgeid,double alpha,double zero,ivector &nod,vector &xh,vector &yh,vector &zh,
                                            long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  long stop;
  double nc1,nc2,nc3;
  vector n(ASTCKVEC(8));
  
  switch (edgeid){
  case 1:{
    nc1=alpha;  nc2=1.0;  nc3=1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[1];
      masnod[numinter][2]=nod[0];
    }
    break;
  }
  case 2:{
    nc1=-1.0; nc2=alpha;  nc3=1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[2];
      masnod[numinter][2]=nod[1];
    }
    break;
  }
  case 3:{
    nc1=alpha; nc2=-1.0;  nc3=1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[2];
      masnod[numinter][2]=nod[3];
    }
    break;
  }
  case 4:{
    nc1=1.0; nc2=alpha;  nc3=1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[3];
      masnod[numinter][2]=nod[0];
    }
    break;
  }
  case 5:{
    nc1=alpha; nc2=1.0;  nc3=-1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[5];
      masnod[numinter][2]=nod[4];
    }
    break;
  }
  case 6:{
    nc1=-1.0; nc2=alpha;  nc3=-1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[6];
      masnod[numinter][2]=nod[5];
    }
    break;
  }
  case 7:{
    nc1=alpha; nc2=-1.0;  nc3=-1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[6];
      masnod[numinter][2]=nod[7];
    }
    break;
  }
  case 8:{
    nc1=1.0; nc2=alpha;  nc3=-1.0;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[7];
      masnod[numinter][2]=nod[4];
    }
    break;
  }
  case 9:{
    nc1=1.0; nc2=1.0;  nc3=alpha;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[4];
      masnod[numinter][2]=nod[0];
    }
    break;
  }
  case 10:{
    nc1=-1.0; nc2=1.0;  nc3=alpha;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[5];
      masnod[numinter][2]=nod[1];
    }
    break;
  }
  case 11:{
    nc1=-1.0; nc2=-1.0;  nc3=alpha;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[6];
      masnod[numinter][2]=nod[2];
    }
    break;
  }
  case 12:{
    nc1=1.0; nc2=-1.0;  nc3=alpha;
    stop=0;
    if (fabs(alpha+1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
      stop=1;
    }
    if (fabs(alpha-1.0)<zero){
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
      stop=1;
    }
    if (stop==0){
      masnod[numinter][0]=2;
      masnod[numinter][1]=nod[7];
      masnod[numinter][2]=nod[3];
    }
    break;
  }
  default:{
    nc1 = nc2 = nc3 = 0.0;
    print_err("invalid edge id edgeid=%ld.", __FILE__, __LINE__, __func__, edgeid);
  }
  }//  end switch (edgeid){
  
  bf_lin_hex_3d (n.a,nc1,nc2,nc3);
  xi[numinter] = ss(n.a,xh.a,xh.n);
  yi[numinter] = ss(n.a,yh.a,yh.n);
  zi[numinter] = ss(n.a,zh.a,zh.n);
  xxi[numinter] = nc1;
  yyi[numinter] = nc2;
  zzi[numinter] = nc3;
  numinter++;

}

void gtopology::end_points_edge_linhex_natcoordinates (long edgeid,double alpha,ivector &nod,vector &xh,vector &yh,vector &zh,
                                                       long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod)
{
  double nc1,nc2,nc3;
  vector n(ASTCKVEC(8));
  
  switch (edgeid){
  case 1:{
    if (alpha<-1.0){
      nc1=-1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
    }
    break;
  }
  case 2:{
    if (alpha<-1.0){
      nc1=-1.0; nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
    }
    if (alpha>1.0){
      nc1=-1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
    }
    break;
  }
  case 3:{
    if (alpha<-1.0){
      nc1=-1.0;  nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
    }
    break;
  }
  case 4:{
    if (alpha<-1.0){
      nc1=1.0; nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
    }
    break;
  }
  case 5:{
    if (alpha<-1.0){
      nc1=-1.0; nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
    }
    break;
  }
  case 6:{
    if (alpha<-1.0){
      nc1=-1.0; nc2=-1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
    }
    if (alpha>1.0){
      nc1=-1.0;  nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
    }
    break;
  }
  case 7:{
    if (alpha<-1.0){
      nc1=-1.0;  nc2=-1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
    }
    break;
  }
  case 8:{
    if (alpha<-1.0){
      nc1=1.0; nc2=-1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
    }
    break;
  }
  case 9:{
    if (alpha<-1.0){
      nc1=1.0; nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[4];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[0];
    }
    break;
  }
  case 10:{
    if (alpha<-1.0){
      nc1=-1.0; nc2=1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[5];
    }
    if (alpha>1.0){
      nc1=-1.0;  nc2=1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[1];
    }
    break;
  }
  case 11:{
    if (alpha<-1.0){
      nc1=-1.0; nc2=-1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[6];
    }
    if (alpha>1.0){
      nc1=-1.0;  nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[2];
    }
    break;
  }
  case 12:{
    if (alpha<-1.0){
      nc1=1.0; nc2=-1.0;  nc3=-1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[7];
    }
    if (alpha>1.0){
      nc1=1.0;  nc2=-1.0;  nc3=1.0;
      masnod[numinter][0]=1;
      masnod[numinter][1]=nod[3];
    }
    break;
  }
  default:{
    nc1 = nc2 = nc3 = 0.0;
    print_err("invalid edge id edgeid=%ld.", __FILE__, __LINE__, __func__, edgeid);
  }
  }//  end switch (edgeid){
  
  bf_lin_hex_3d (n.a,nc1,nc2,nc3);
  xi[numinter] = ss(n.a,xh.a,xh.n);
  yi[numinter] = ss(n.a,yh.a,yh.n);
  zi[numinter] = ss(n.a,zh.a,zh.n);
  xxi[numinter] = nc1;
  yyi[numinter] = nc2;
  zzi[numinter] = nc3;
  numinter++;

}


void gtopology::initiate_elem_dof_ranges()
{
  long i, j, ndofe;
  ivector cn;

  eldofr = new dofrange[ne];
  for(i=0; i<ne; i++){
    ndofe = give_ndofe(i);
    reallocv(RSTCKIVEC(ndofe, cn));
    give_code_numbers(i, cn.a);
    eldofr[i].mindof = LONG_MAX;
    eldofr[i].deltadof = 0;
    for(j = 0; j<ndofe; j++){
      if (cn(j) < 1)
        continue;
      if (cn(j) < eldofr[i].mindof)
        eldofr[i].mindof = cn(j);
      if (cn(j) > eldofr[i].deltadof)
        eldofr[i].deltadof = cn(j);
    }
    eldofr[i].deltadof -= eldofr[i].mindof;
  }
  return;
}


void gtopology::initiate_omp_elem_order()
{
  long i, j;
  std::vector<std::pair<dofrange, long > > dr;
  ompelord = new long[ne];
  dr.reserve(ne);
  for(i=0; i<ne; i++)
    dr.push_back(std::make_pair(eldofr[i], i));
  std::sort(dr.begin(), dr.end(), dofrange::compare);

  // fill ordering vector 
  long ne_perthread = ne/16;
  for(i=0; i<ne_perthread; i++){
    for(j=0; j<16; j++)
      ompelord[i*16+j] = dr[j*ne_perthread+i].second;
  }
  for(i=0; i<ne%16; i++)
    ompelord[16*ne_perthread + i] = dr[16*ne_perthread + i].second;
  
  return;
}



/**
  The function creates a bounding box of the given problem.
*/
void gtopology::give_bounding_box(boundbox &bb)
{
  double xmin = DBL_MAX;
  double xmax = DBL_MIN;
  double ymin = DBL_MAX;
  double ymax = DBL_MIN;
  double zmin = DBL_MAX;
  double zmax = DBL_MIN;
  
  for (long i=0; i<nn; i++){
    if (gnodes[i].x < xmin)  xmin = gnodes[i].x;
    if (gnodes[i].x > xmax)  xmax = gnodes[i].x;
    if (gnodes[i].y < ymin)  ymin = gnodes[i].y;
    if (gnodes[i].y > ymax)  ymax = gnodes[i].y;
    if (gnodes[i].z < zmin)  zmin = gnodes[i].z;
    if (gnodes[i].z > zmax)  zmax = gnodes[i].z;
  }
  bb.init(xmin, xmax, ymin, ymax, zmin, zmax);
}
