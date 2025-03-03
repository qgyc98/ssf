#include "edget.h"
#include "globalt.h"
#include "globmatt.h"
#include <string.h>

edget::edget (void)
{
  tm1 = NULL;
  tm2 = NULL;
  idm1 = NULL;
  idm2 = NULL;
  
  jumpfn=NULL;
  jumpln=NULL;
  
  //  the number of components in the array edval
  ncedval=0;
  //  edge values
  edval = NULL;
}

edget::~edget (void)
{
  delete [] tm1;
  delete [] tm2;
  delete [] idm1;
  delete [] idm2;
  
  delete [] jumpfn;
  delete [] jumpln;
  
  //  the number of components in the array edval
  ncedval=0;
  //  edge values
  delete [] edval;
}

/**
   function initializes edge
   
   @param edid - edge id
   
   JK, 6.8.2008
*/
void edget::init (long edid)
{
  long i,eid1,eid2,ntm;
  
  //  number of the first adjacent element
  eid1=Gtt->gedges[edid].adjel[0];
  //  number of the second adjacent element
  eid2=Gtt->gedges[edid].adjel[1];
  
  //  number of transported media
  ntm=Tp->ntm*Tp->ntm;

  //  allocation of arrays for material types and id of materials
  tm1 = new mattypet [ntm];
  idm1 = new long [ntm];
  tm2 = new mattypet [ntm];
  idm2 = new long [ntm];
  
  for (i=0;i<ntm;i++){
    tm1[i]=Tt->elements[eid1].tm[i];
    idm1[i]=Tt->elements[eid1].idm[i];

    tm2[i]=Tt->elements[eid2].tm[i];
    idm2[i]=Tt->elements[eid2].idm[i];
  }
  
  jumpfn = new double [ntm];
  jumpln = new double [ntm];
}

/**
   function computes jumps across the edge
   
   @param edid - edge id
   
   JK, 6.8.2008
*/
void edget::compute_jump (long edid)
{
  long lcid,ndofn,fn1,fn2,ln1,ln2;
  vector nvf1, invf1, nvf2, invf2, nvl1, invl1, nvl2, invl2;

  //  load case id must be equal to zero in this type of problem
  lcid=0;
  
  //  first nodes
  fn1 = Gtt->gedges[edid].nlist[0];
  fn2 = Gtt->gedges[edid].nlist[1];
  //  last nodes
  ln1 = Gtt->gedges[edid].nlist[2];
  ln2 = Gtt->gedges[edid].nlist[3];
  
  //  number of DOFs on selected node
  ndofn = Tt->give_ndofn (fn1);
  reallocv(RSTCKVEC(ndofn, nvf1));
  reallocv(RSTCKVEC(ndofn, nvf2));
  reallocv(RSTCKVEC(ndofn, invf1));
  reallocv(RSTCKVEC(ndofn, invf2));
  //  nodal values on selected nodes
  nodalval (fn1, nvf1);
  nodalval (fn2, nvf2);
  //  initial nodal values on selected nodes
  initnodval (fn1, invf1);
  initnodval (fn2, invf2);
  
  //  number of DOFs on selected node
  ndofn = Tt->give_ndofn (ln1);
  reallocv(RSTCKVEC(ndofn, nvl1));
  reallocv(RSTCKVEC(ndofn, nvl2));
  reallocv(RSTCKVEC(ndofn, invl1));
  reallocv(RSTCKVEC(ndofn, invl2));
  //  nodal values on selected nodes
  nodalval (ln1, nvl1);
  nodalval (ln2, nvl2);
  //  initial nodal values on selected nodes
  initnodval (ln1, invl1);
  initnodval (ln2, invl2);
  
  //  computation of jumps along discontinuity
  Tm->values_transformation (tm1[0],idm1[0],nvf1.a,invf1.a,tm2[0],idm2[0],nvf2.a,invf2.a,jumpfn);
  Tm->values_transformation (tm1[0],idm1[0],nvl1.a,invl1.a,tm2[0],idm2[0],nvl2.a,invl2.a,jumpln);
}


/**
   function computes jumps in selected node on edge
   
   @param edid - edge id
   @param nodid - node id
          nodid=1 - jumps in first node are computed
	  nodid=2 - jumps in last node are computed
   @param lcid - load case id
   @param ncf - normal component of flux
   
   JK, 11.2.2011
*/
void edget::compute_node_jump (long edid, long nodid,long /*lcid*/,double /*ncf*/)
{
  long ndofn,n1,n2;
  vector nv1, inv1, nv2, inv2;

  //  default values of node numbers
  n1=-1;
  n2=-1;
  
  if (nodid==1){
    //  first nodes
    n1 = Gtt->gedges[edid].nlist[0];
    n2 = Gtt->gedges[edid].nlist[1];
  }
  if (nodid==2){
    //  last nodes
    n1 = Gtt->gedges[edid].nlist[2];
    n2 = Gtt->gedges[edid].nlist[3];
  }
  
  if (n1<0 || n2<0){
    print_err("wrong node is required",__FILE__,__LINE__,__func__);
  }
  
  //  number of DOFs on selected node
  ndofn = Tt->give_ndofn (n1);
  
  reallocv(RSTCKVEC(ndofn, nv1));
  reallocv(RSTCKVEC(ndofn, nv2));
  reallocv(RSTCKVEC(ndofn, inv1));
  reallocv(RSTCKVEC(ndofn, inv2));

  //  nodal values on selected nodes
  nodalval (n1, nv1);
  nodalval (n2, nv2);

  //  initial nodal values on selected nodes
  initnodval (n1, inv1);
  initnodval (n2, inv2);

  //  computation of jumps along discontinuity in selected node
  if (nodid==1){
    Tm->values_transformation (tm1[0],idm1[0],nv1.a,inv1.a,tm2[0],idm2[0],nv2.a,inv2.a,jumpfn);
  }
  if (nodid==2){
    Tm->values_transformation (tm1[0],idm1[0],nv1.a,inv1.a,tm2[0],idm2[0],nv2.a,inv2.a,jumpln);
  }
}


/**
   function initializes edge values
   
   the array edval stores average temperature on the edge
   it is used in problems with radiation
   
   JK, 25.7.2011
*/
void edget::init_edval (void)
{
  ncedval = 1;
  edval = new double [ncedval];
}

/**
   function stores edge values
   
   the array edval stores average temperature on the edge
   it is used in problems with radiation
   
   @param v - array with nodal values
   
   JK, 25.7.2011
*/
void edget::store_edval (double *v)
{
  long i;
  
  for (i=0;i<ncedval;i++){
    edval[i]=v[i];
  }
}

/**
   function gives edge values
   
   the array edval gives average temperature on the edge
   it is used in problems with radiation
   
   @param v - array with nodal values
   
   JK, 25.7.2011
*/
void edget::give_edval (double *v)
{
  long i;
  
  for (i=0;i<ncedval;i++){
    v[i]=edval[i];
  }
}

