#include "endnodet.h"
#include "globalt.h"
#include "globmatt.h"
#include <string.h>

endnodet::endnodet (void)
{
  
  tm1 = NULL;
  tm2 = NULL;
  idm1 = NULL;
  idm2 = NULL;
  
  jump=NULL;
}

endnodet::~endnodet (void)
{
  delete [] tm1;
  delete [] tm2;
  delete [] idm1;
  delete [] idm2;
  
  delete [] jump;
}

/**
   function initializes end node
   
   @param enid - end node id
   
   JK, 17.10.2008
*/
void endnodet::init (long enid)
{
  long i,eid1,eid2,ntm;
  
  //  number of the first adjacent element
  eid1=Gtt->endnodes[enid].adjel[0];
  //  number of the second adjacent element
  eid2=Gtt->endnodes[enid].adjel[1];
  
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
  
  jump = new double [ntm];
}

/**
   function computes jumps in node
   
   @param enid - end node id
   
   JK, 17.10.2008
*/
void endnodet::compute_jump (long enid)
{
  long lcid,ndofn,n1,n2;
  vector nvf1, invf1, nvf2, invf2;

  //  load case id must be equal to zero in this type of problem
  lcid=0;
  
  //  node numbers
  n1 = Gtt->endnodes[enid].fn;
  n2 = Gtt->endnodes[enid].ln;
  
  //  number of DOFs on selected node
  ndofn = Tt->give_ndofn (n1);
  reallocv(RSTCKVEC(ndofn, nvf1));
  reallocv(RSTCKVEC(ndofn, nvf2));
  reallocv(RSTCKVEC(ndofn, invf1));
  reallocv(RSTCKVEC(ndofn, invf2));
  //  nodal values on selected nodes
  nodalval (n1, nvf1);
  nodalval (n2, nvf2);
  //  initial nodal values on selected nodes
  initnodval (n1, invf1);
  initnodval (n2, invf2);
  
  //  computation of jumps along discontinuity
  Tm->values_transformation (tm1[0],idm1[0],nvf1.a,invf1.a,tm2[0],idm2[0],nvf2.a,invf2.a,jump);
  
 /* fprintf (Outt,"\n enid %ld  time  %lf",enid,Tp->time);
  fprintf (Outt,"\n node 1  %le  %le  %le",
	   invf1[0],nvf1[0]-invf1[0],nvf1[0]);
  fprintf (Outt,"\n node 2  %le  %le  %le    skok %le  %le",
	   invf2[0],nvf2[0]-invf2[0],nvf2[0],jump[0],jump[1]);
  */
}
