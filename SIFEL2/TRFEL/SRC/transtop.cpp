#include <string.h>
#include <limits.h>
#include <stdio.h>
#include "transtop.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"
#include "ipmap.h"

transtop::transtop ()
{
  tnbn = 0;
  //  number of nodes
  nn=0;
  // number of constrained nodes
  ncn = 0;
  //  number of elements
  ne=0;
  
  //  array of nodes
  nodes=NULL;
  //  array of elements
  elements=NULL;
  //  array of end nodes
  enodes=NULL;
  //  array of edges
  edges=NULL;
  
  //nnj = 0;
  
  nedis=0;
  lnnd=NULL;
  lnd=NULL;
  
  //  number of cycles on nodes
  ncycl = NULL;
  
  //master node number
  mnn = 0;

  //JM
  //jnmap = NULL;
  //prevval = NULL;
  aux = NULL;
  view_factor = NULL;
  mvf = NULL;
  mvfjk = NULL;
  renel = NULL;
}

transtop::~transtop ()
{
  long i;
  
  if (nodes!=NULL)
    delete [] nodes;  
  if (elements!=NULL)
    delete [] elements;
  if (enodes!=NULL)
    delete [] enodes;
  if (edges!=NULL)
    delete [] edges;
  if (lnnd!=NULL)
    delete [] lnnd;
  
  /*
  for (i = 0;i<nnj;i++)
  {
    //if (jnmap != NULL)
    //delete [] jnmap[i];
    if (prevval != NULL)
      delete [] prevval[i];
  }
  delete [] prevval;
  //delete [] jnmap;
  */
   
  if (lnd!=NULL){
    for (i=0;i<nedis;i++){
      delete [] lnd[i];
    }
    delete [] lnd;
  }
  if (ncycl!=NULL)
    delete [] ncycl;
}

/**
   function reads basic data about topology
   
   @param in - input stream
*/
void transtop::read (XFILE *in)
{
  long i,j,k,ndofn;
  
  // ************
  //  node data
  // ************
  xfscanf (in,"%k%ld","number_of_nodes",&nn);
  Gtt->alloc_nodes (nn);
  
  if (Mesprt==1)  fprintf (stdout,"\n number of nodes  %ld",nn);
  nodes = new nodet [nn];
  for (i=0;i<nn;i++){
    j=nn+1;
    xfscanf (in,"%ld",&j);
    if (j<1)
      print_err("node number in function read is less than 1",__FILE__,__LINE__,__func__);
    if (j>nn)
      print_err("node number in function read is greater than number of nodes",__FILE__,__LINE__,__func__);
    Gtt->gnodes[j-1].read (in);
    //    ndofn = Gtt->gnodes[j-1].ndofn;
    ndofn = give_ndofn(j-1);
    nodes[j-1].read (in,ndofn);
  }
  
  // ********************
  //  constrained nodes
  // ********************
  xfscanf (in,"%k%ld","number_of_constraints",&ncn);
  if (Mesprt==1)  fprintf (stdout,"\n number of constrained nodes  %ld",ncn);
  for (i=0;i<ncn;i++){
    k=nn+1;
    xfscanf (in,"%ld",&k);
    //fprintf (Outt,"\n%ld",k);
    //fflush(Outt);
    if (k<1)
      print_err("number of constrained node in function read is less than 1",__FILE__,__LINE__,__func__);
    if (k>nn){
      print_err("number of constrained node in function read is greater than number of all nodes",__FILE__,__LINE__,__func__);
    }
    Gtt->gnodes[k-1].constr (Gtt->dofcontr,in);
  }
  
  
  // ************************************
  //  master node only for homogenization
  // ************************************
  if(Tp->homogt == 1){
    xfscanf (in,"%k%ld","masternode_number",&mnn);
    if (Mesprt==1)  fprintf (stdout,"\n master node number for homogenization %ld",mnn);
    mnn = mnn-1;
  }

  // ***************
  //  element data
  // ***************
  xfscanf (in,"%k%ld","number_of_elements",&ne);
  Gtt->alloc_elements (ne);
  if (Mesprt==1)  fprintf (stdout,"\n number of elements  %ld",ne);
  elements = new elementt [ne];
  
  
  /*
  if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){
    xfscanf (in,"%ld",&nedis);
    lnd = new long* [nedis];
    lnnd = new long [nedis];
    for (i=0;i<nedis;i++){
      xfscanf (in,"%ld %ld",&j,&lnnd[i]);
      elements[j-1].discont=i;
      lnd[i] = new long [lnnd[i]];
      for (j=0;j<lnnd[i];j++){
	xfscanf (in,"%ld",&lnd[i][j]);
	lnd[i][j]--;
      }
    }
  }
  */



  for (i=0;i<ne;i++){
    j=ne+1;
    xfscanf (in,"%ld",&j);
    if (j<1)
      print_err("element number in function read is less than 1",__FILE__,__LINE__,__func__);
    if (j>ne){
      print_err("element number in function read is greater than total number of elements",__FILE__,__LINE__,__func__);
    }
    elements[j-1].read (in,j-1);
  }
  
  //  initialization of auxinf variable
  /*
  for (i=0;i<ne;i++){
    nne = give_nne (i);
    order = give_degree (i);
    dim = give_dimension (i);
    
    Gtt->gelements[i].auxinf = nne*100+order*10+dim;
  }
  */

  //  function searches hanging nodes
  //  variable gelements[i].nne is modified
  //  array gelements[i].nodes is modified
  Gtt->searching_hanging_nodes();
  assemble_master_node_weights_vec();
  searching_hanging_nodes();
  Gtt->hang_nodes_check();


  //  allocation and initiation of arrays lnso and leso
  //  lnso - list of nodes switched on
  //  leso - list of elements switched on
  Gtt->lneso_init ();
  
  
  if (Tp->tprob == growing_np_problem){
    Gtt->auxinf_init ();
  }
  if (Tp->tprob == growing_np_problem_nonlin){
    Gtt->auxinf_init ();
  }
  
  
  if (Gtt->rst==1){
    //  sequential topology is read
    //  it is used in the case of sequential version of
    //  any domain decomposition method
    Gtt->read_seq_top (in);
  }
  if (Gtt->rst==2){
    //  sequential topology is read in the form
    //  of Boolean matrices
    Tp->ssle->feti.read_booldata (in);
  }

  
  /*
  //  JK, 15.11.2007
  //  reading additional informations about mesh due to domain decomposition method
  //  problem will be solved by sequential FETI method
  if (Tp->ssle->tlinsol == sfeti){
    //Gtt->read_seq_top (Tp->ssle->feti.ns,in);
    Gtt->read_seq_top (in);
  }
  //  JK, 15.11.2007
  
  //  old version of discontinuity
  if (Tp->tprob == discont_nonstat_problem){
    //Gtt->read_ltg1 (in);
    Gtt->cngen=2;
    xfscanf (in,"%ld",&nnj);

    jnmap = new long* [nnj];
    for (i=0;i<nnj;i++){
      jnmap[i] = new long [2];
      xfscanf (in,"%ld %ld",&jnmap[i][0],&jnmap[i][1]);
      jnmap[i][0]--;
      jnmap[i][1]--;
    }
    
    prevval = new double* [nnj];
    for (i=0;i<nnj;i++){
      ndofn = give_ndofn (jnmap[i][0]);
      prevval[i] = new double [ndofn+2];
      for (j=0;j<ndofn+2;j++){
	prevval[i][j]=0.0;
      }
    }
  }
  */
  
  
  //  new version of discontinuity
  //  13.12.2007, JK
  /*
  if (Tp->tprob == discont_nonstat_problem){
    xfscanf (in,"%ld",&ns);
    Gtt->read_feti (ns,in);
    Gtt->cngen=3;
  }
  */
  //  13.12.2007, JK

  /*
  if (Tp->ssle->tlinsol == saddle_point){
    //Gtt->read_seq_top (Tp->ssle->sp->ns,in);
    Gtt->read_seq_top (in);
  }
  */
}


/**
   function prints basic data about topology
   
   @param in - input stream

   18/12/2012, TKr
*/
void transtop::print (FILE *out)
{
  long i;
  
  // ************
  //  node data
  // ************
  fprintf (out,"\n%ld #number of nodes\n\n",nn);
  
  for (i=0;i<nn;i++){
    fprintf (out,"\n %ld",i+1);
    Gtt->gnodes[i].print (out);
    nodes[i].print (out);
  }
  
  // ********************
  //  constrained nodes
  // ********************
  fprintf (out,"\n%ld #number of constrain nodes\n\n",ncn);
  for (i=0;i<ncn;i++){
    fprintf (out,"\n %ld",i);
    Gtt->gnodes[i].print_constr (Gtt->dofcontr,out);
  }
  
  
  // ************************************
  //  master node only for homogenization
  // ************************************
  if(Tp->homogt == 1){
    fprintf (out,"\n %ld #master node number\n\n",mnn);
  }

  // ***************
  //  element data
  // ***************
  fprintf (out,"\n %ld #number of elements\n\n",ne);
  for (i=0;i<ne;i++){
    fprintf (out,"\n %ld",i+1);

    elements[i].print (out,i);
  }
}


/**
   function returns element type
   
   @param eid - element id
*/
elemtypet transtop::give_elem_type (long eid)
{
  return (elements[eid].te);
}

/**
   function returns appropriate DOF of node
   
   @param nid - node id
   @param n - number of required DOF
*/
long transtop::give_dof (long nid,long n)
{
  return Gtt->give_dof (nid,n);
}

/**
   Function returns numbers of nodes defining apropriate element.

   There are two arrays for node numbers of elements in the code.
   If there is no hanging node in the problem solved, nodes
   of elements are stored only in the objects Gtt->gelements.
   If there are hanging nodes, the elements attached to the
   hanging nodes contain the original slave nodes in the array
   Tt->elements.nodes while the array Gtm->gelements.nodes
   stores master nodes.

   This function returns the original node numbers
   i.e. the slave node numbers and it is used e.g. in graphic output.

   @param[in] eid - element id
   @param[out] nodes - array containing nodes

   @return The function returns node numbers in the argument nodes.

   22.7.2001
*/
void transtop::give_elemnodes (long eid,ivector &nodes)
{
  Gtt->give_original_nodes (eid,nodes);
}



/**
   function assembles code numbers of actual element
   
   @param eid[in] - element id
   @param cn[out] - array containing code numbers on element
   
   25.6.2001
*/
void transtop::give_code_numbers (long eid,long *cn)
{
  Gtt->give_code_numbers (eid,cn);
}



/**
   function assembles code numbers of a single medium on actual element
   
   @param eid - element id
   @param ri - row index = medium id
   @param cn - array containing code numbers on element
   
   14.3.2013
*/
void transtop::give_medium_code_numbers (long eid,long ri,long *cn)
{
  long i,ndofe,nne,ntm;
  ivector ecn;
  
  //  the number of DOFs on element
  ndofe = Gtt->give_ndofe (eid);
  reallocv(RSTCKIVEC(ndofe, ecn));
  
  Gtt->give_code_numbers (eid, ecn.a);
  
  //  the number of nodes on element
  nne = Gtt->give_nne (eid);
  
  // the number of transported media
  ntm = Tp->ntm;
  
  for (i=0;i<nne;i++){
    cn[i]=ecn[i*ntm+ri];
  }
}

/**
   function extracts code numbers of actual node
   
   @param nid - node id
   @param cn - array containing code numbers on node
   
   7.12.2002
*/
void transtop::give_node_code_numbers (long nid,long *cn)
{
  Gtt->give_node_code_numbers (nid,cn);
}

/**
   function returns node coordinates of the appropriate element

   19.7.2001
*/
void transtop::give_node_coord1d (vector &x,long eid)
{
  Gtt->give_node_coord1d (x,eid);
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,y - vectors containing node coordinates
   @param eid - element id

   19.7.2001
*/
void transtop::give_node_coord2d (vector &x,vector &y,long eid)
{
  Gtt->give_node_coord2d (x,y,eid);
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,y,z - vectors containing coordinates
   @param eid - element id
   
   19.7.2001
*/
void transtop::give_node_coord3d (vector &x,vector &y,vector &z,long eid)
{
  Gtt->give_node_coord3d (x,y,z,eid);
}

/**
   function returns number of DOFs of element
   
   @param eid - element id
   
   JK, 22.2.2002
*/
long transtop::give_ndofe (long eid)
{
  long ndofe;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{         ndofe=Lbt->ndofe;      break; }
  case barlint3d:{       ndofe=Lbt3d->ndofe;      break; }
  case barlintax:{       ndofe=Lbat->ndofe;     break; }
  case barquadt:{        ndofe=Qbt->ndofe;      break; }
  case barquadtax:{      ndofe=Qbat->ndofe;     break; }
  case trlint:{          ndofe=Ltt->ndofe;      break; }
  case trlaxisym:{       ndofe=Ltat->ndofe;     break; }
  case quadlint:{        ndofe=Lqt->ndofe;      break; }
  case quadquadt:{       ndofe=Qqt->ndofe;      break; }
  case quadquadtax:{     ndofe=Qqat->ndofe;     break; }
  case quadlaxisym:{     ndofe=Lqat->ndofe;     break; }
  case ifacequadel:{     ndofe=Ifcquadt->ndofe; break; }
  case lineartett:{      ndofe=Ltett->ndofe;    break; }
  case linearhext:{      ndofe=Lht->ndofe;      break; }
  case quadratichext:{   ndofe=Qht->ndofe;      break; }
  case linearwedget:{    ndofe=Lwt->ndofe;      break; }
  case gen2del:{         ndofe=G2d->ndofe;      break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ndofe;
  
}


/**
   The function returns array of number of DOFs of element for particular 
   tarnsported media. Dimensions of the array dofe is (ntm x ntm)
   
   @param eid - element id
   
   TKo, 19.12.2017
*/
long** transtop::give_dofe (long eid)
{
  long **dofe=NULL;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{         dofe=Lbt->dofe;      break; }
  case barlint3d:{       dofe=Lbt3d->dofe;    break; }
  case barlintax:{       dofe=Lbat->dofe;     break; }
  case barquadt:{        dofe=Qbt->dofe;      break; }
  case barquadtax:{      dofe=Qbat->dofe;     break; }
  case trlint:{          dofe=Ltt->dofe;      break; }
  case trlaxisym:{       dofe=Ltat->dofe;     break; }
  case quadlint:{        dofe=Lqt->dofe;      break; }
  case quadquadt:{       dofe=Qqt->dofe;      break; }
  case quadquadtax:{     dofe=Qqat->dofe;     break; }
  case quadlaxisym:{     dofe=Lqat->dofe;     break; }
  case ifacequadel:{     dofe=Ifcquadt->dofe; break; }
  case lineartett:{      dofe=Ltett->dofe;    break; }
  case linearhext:{      dofe=Lht->dofe;      break; }
  case quadratichext:{   dofe=Qht->dofe;      break; }
  case linearwedget:{    dofe=Lwt->dofe;      break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return dofe;
}



/**
   The function returns array of DOF ordering of element for particular 
   tarnsported media. Dimensions of the returned array ordering is (ntm x nne).
   
   @param eid - element id
   
   TKo, 19.12.2017
*/
long** transtop::give_ordering (long eid)
{
  long **ordering=NULL;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{         ordering=Lbt->ordering;      break; }
  case barlint3d:{       ordering=Lbt3d->ordering;      break; }
  case barlintax:{       ordering=Lbat->ordering;     break; }
  case barquadt:{        ordering=Qbt->ordering;      break; }
  case barquadtax:{      ordering=Qbat->ordering;     break; }
  case trlint:{          ordering=Ltt->ordering;      break; }
  case trlaxisym:{       ordering=Ltat->ordering;     break; }
  case quadlint:{        ordering=Lqt->ordering;      break; }
  case quadquadt:{       ordering=Qqt->ordering;      break; }
  case quadquadtax:{     ordering=Qqat->ordering;     break; }
  case quadlaxisym:{     ordering=Lqat->ordering;     break; }
  case ifacequadel:{     ordering=Ifcquadt->ordering; break; }
  case lineartett:{      ordering=Ltett->ordering;    break; }
  case linearhext:{      ordering=Lht->ordering;      break; }
  case quadratichext:{   ordering=Qht->ordering;      break; }
  case linearwedget:{    ordering=Lwt->ordering;      break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ordering;
}



/**
   function returns number of DOFs of node
   
   @param nid - node id
*/
long transtop::give_ndofn (long nid)
{
  long ndofn;
  
  ndofn = Gtt->give_ndofn (nid);
  if (ndofn<0){
    //  the node is a hanging node
    ndofn = Gtt->give_ndofn (Gtt->gnodes[nid].mnodes[0]);
  }
  return ndofn;
}

/**
   function returns total number of integration points on element
   
   @param eid - element id
   
   30.12.2001
*/
long transtop::give_tnip (long eid)
{
  if (eid < ne-1)
    return elements[eid+1].ipp[0][0] - elements[eid].ipp[0][0];

  return Tm->tnip-elements[eid].ipp[0][0];
}



/**
   function returns number of components
   
   @param eid - element id
   
   17.1.2003 - Tomas Krejci
*/
long transtop::give_ncomp (long eid)
{
  long ncomp;
  elemtypet te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case barlint:{        ncomp=Lbt->ncomp;      break; }
  case barlint3d:{      ncomp=Lbt3d->ncomp;      break; }
  case barlintax:{      ncomp=Lbat->ncomp;     break; }
  case barquadt:{       ncomp=Qbt->ncomp;      break; }
  case barquadtax:{     ncomp=Qbat->ncomp;     break; }
  case trlint:{         ncomp=Ltt->ncomp;      break; }
  case trlaxisym:{      ncomp=Ltat->ncomp;     break; }
  case quadlint:{       ncomp=Lqt->ncomp;      break; }
  case quadquadt:{      ncomp=Qqt->ncomp;      break; }
  case quadquadtax:{    ncomp=Qqat->ncomp;     break; }
  case quadlaxisym:{    ncomp=Lqat->ncomp;     break; }
  case ifacequadel:{    ncomp=Ifcquadt->ncomp; break; }
  case lineartett:{     ncomp=Ltett->ncomp;    break; }
  case linearhext:{     ncomp=Lht->ncomp;      break; }
  case quadratichext:{  ncomp=Qht->ncomp;      break; }
  case linearwedget:{   ncomp=Lwt->ncomp;      break; }
  case gen2del:{        ncomp=G2d->ncomp;      break; }
    
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ncomp;
}

/**
   function returns number of edges of element
   
   @param eid - element id
   
   JK, 22.2.2002
*/
long transtop::give_ned (long eid)
{
  long ned;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{        ned=Lbt->ned;    break; }
  case barlint3d:{      ned=Lbt3d->ned;  break; }
  case barlintax:{      ned=Lbat->ned;   break; }
  case barquadt:{       ned=Qbt->ned;    break; }    
  case barquadtax:{     ned=Qbat->ned;   break; }    
  case trlint:{         ned=Ltt->ned;    break; }
  case trlaxisym:{      ned=Ltat->ned;   break; }
  case quadlint:{       ned=Lqt->ned;    break; }
  case quadlaxisym:{    ned=Lqat->ned;   break; }
  case quadquadt:{      ned=Qqt->ned;    break; }
  case quadquadtax:{    ned=Qqat->ned;   break; }
  case lineartett:{     ned=Ltett->ned;  break; }
  case linearhext:{     ned=Lht->ned;    break; }
  case quadratichext:{  ned=Qht->ned;    break; }
  case linearwedget:{   ned=Lwt->ned;    break; }
  case gen2del:{        ned=0;           break; }

  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ned;
  
}

/**
   function returns number of nodes on one element edge
   
   @param eid - number of element
   
   30.12.2001
*/
long transtop::give_nned (long eid)
{
  long nned;
  elemtypet te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case barlint:{         nned=Lbt->nned;     break; }
  case barlint3d:{       nned=Lbt3d->nned;   break; }
  case barlintax:{       nned=Lbat->nned;    break; }
  case barquadt:{        nned=Qbt->nned;     break; }    
  case barquadtax:{      nned=Qbat->nned;    break; }    
  case trlint:{          nned=Ltt->nned;     break; }
  case trlaxisym:{       nned=Ltat->nned;    break; }
  case quadlint:{        nned=Lqt->nned;     break; }
  case quadquadt:{       nned=Qqt->nned;     break; }
  case quadquadtax:{     nned=Qqat->nned;    break; }
  case quadlaxisym:{     nned=Lqat->nned;    break; }
  case lineartett:{      nned=Ltett->nned;   break; }
  case linearhext:{      nned=Lht->nned;     break; }
  case quadratichext:{   nned=Qht->nned;     break; }
  case linearwedget:{    nned=Lwt->nned;     break; }
  case gen2del:{         nned=0;             break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nned;
}

/**
   function returns number of surfaces of element
   
   @param eid - element id
   
   JK, 9.3.2002
*/
long transtop::give_nsurf (long eid)
{
  long nsurf;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{         nsurf=0;     break; }
  case barlint3d:{       nsurf=0;     break; }
  case barlintax:{       nsurf=0;     break; }
  case barquadt:{        nsurf=0;     break; }
  case barquadtax:{      nsurf=0;     break; }
  case trlint:{          nsurf=0;     break; }
  case trlaxisym:{       nsurf=0;     break; }
  case quadlint:{        nsurf=0;     break; }
  case quadquadt:{       nsurf=0;     break; }
  case quadquadtax:{     nsurf=0;     break; }
  case quadlaxisym:{     nsurf=0;     break; }
  case gen2del:{         nsurf=0;     break; }

  case lineartett:{      nsurf=Ltett->nsurf;   break; }
  case linearhext:{      nsurf=Lht->nsurf;     break; }
  case quadratichext:{   nsurf=Qht->nsurf;     break; }
  case linearwedget:{    nsurf=Lwt->nsurf;     break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nsurf;
  
}

/**
   function returns number of nodes on one surface of element
   
   @param eid - element id
   
   JK, 24.8.2004
*/
long transtop::give_nnsurf (long eid)
{
  long nnsurf;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case lineartett:{      nnsurf=Ltett->nnsurf;   break; }
  case linearhext:{      nnsurf=Lht->nnsurf;     break; }
  case quadratichext:{   nnsurf=Qht->nnsurf;     break; }
  case linearwedget:{    nnsurf=Lwt->nnsurf;     break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nnsurf;
  
}

/**
   function returns number of nodes of element
   
   @param eid - element id
   
   JK, 22.2.2002
*/
long transtop::give_nne (long eid)
{
  long nne;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{          nne=Lbt->nne;      break; }
  case barlint3d:{        nne=Lbt3d->nne;      break; }
  case barlintax:{        nne=Lbat->nne;     break; }
  case barquadt:{         nne=Qbt->nne;      break; }
  case barquadtax:{       nne=Qbat->nne;     break; }
  case trlint:{           nne=Ltt->nne;      break; }
  case trlaxisym:{        nne=Ltat->nne;     break; }
  case quadlint:{         nne=Lqt->nne;      break; }
  case quadquadt:{        nne=Qqt->nne;      break; }
  case quadquadtax:{      nne=Qqat->nne;     break; }
  case quadlaxisym:{      nne=Lqat->nne;     break; }
  case ifacequadel:{      nne=Ifcquadt->nne; break; }
  case lineartett:{       nne=Ltett->nne;    break; }
  case linearhext:{       nne=Lht->nne;      break; }
  case quadratichext:{    nne=Qht->nne;      break; }
  case linearwedget:{     nne=Lwt->nne;      break; }
  case gen2del:{          nne=G2d->nne;      break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nne;
  
}

/**
   function returns number of integration points of element
   for conductivity %matrix

   @param eid - element id
   @param ri,ci - row and column indices
   
   JK, 22.2.2002
*/
long transtop::give_nip (long eid,long ri,long ci)
{
  long nip;
  elemtypet te;
  
  te = give_elem_type (eid);

  switch (te){
  case barlint:{          nip=Lbt->nip[ri][ci];      break; }
  case barlint3d:{        nip=Lbt3d->nip[ri][ci];    break; }
  case barlintax:{        nip=Lbat->nip[ri][ci];     break; }
  case barquadt:{         nip=Qbt->nip[ri][ci];      break; }
  case barquadtax:{       nip=Qbat->nip[ri][ci];     break; }
  case trlint:{           nip=Ltt->nip[ri][ci];      break; }
  case trlaxisym:{        nip=Ltat->nip[ri][ci];     break; }
  case quadlint:{         nip=Lqt->nip[ri][ci];      break; }
  case quadquadt:{        nip=Qqt->nip[ri][ci];      break; }
  case quadquadtax:{      nip=Qqat->nip[ri][ci];     break; }
  case quadlaxisym:{      nip=Lqat->nip[ri][ci];     break; }
  case ifacequadel:{      nip=Ifcquadt->nip[ri][ci]; break; }
  case lineartett:{       nip=Ltett->nip[ri][ci];    break; }
  case linearhext:{       nip=Lht->nip[ri][ci];      break; }
  case quadratichext:{    nip=Qht->nip[ri][ci];      break; }
  case linearwedget:{     nip=Lwt->nip[ri][ci];      break; }
  case gen2del:{          nip=G2d->nip[ri][ci];      break; }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nip;
  
}

/**
  The function returns the integration order of conductivity %matrix 
  for the given block.

  @param eid - element id
  @param ri  - row index (index of transported media)
  @param ci  - column index (index of transported media)

  @return The function returns integration order of the given element and block ids.

  Created by Tomas Koudelka, 11.11.2013
*/
long transtop::give_intordkm(long eid, long ri, long ci)
{
  elemtypet te;
  long ret = 0;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    ret = Lbt->intordkm[ri][ci];
    break;
  }
  case barlint3d:{
    ret = Lbt3d->intordkm[ri][ci];
    break;
  }
  case barlintax:{
    ret = Lbat->intordkm[ri][ci];
    break;
  }
  case barquadt:{
    ret = Qbt->intordkm[ri][ci];
    break;
  }
  case barquadtax:{
    ret = Qbat->intordkm[ri][ci];
    break;
  }
  case trlint:{
    ret = Ltt->intordkm[ri][ci];
    break;
  }
  case trlaxisym:{
    ret = Ltat->intordkm[ri][ci];
    break;
  }
  case quadlint:{
    ret = Lqt->intordkm[ri][ci];
    break;
  }
  case quadquadt:{
    ret = Qqt->intordkm[ri][ci];
    break;
  }
  case quadquadtax:{
    ret = Qqat->intordkm[ri][ci];
    break;
  }
  case quadlaxisym:{
    ret = Lqat->intordkm[ri][ci];
    break;
  }
  case lineartett:{
    ret = Ltett->intordkm[ri][ci];
    break;
  }
  case linearhext:{
    ret = Lht->intordkm[ri][ci];
    break;
  }
  case quadratichext:{
    ret = Qht->intordkm[ri][ci];
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ret;
}



/**
  The function returns the integration order of capacity %matrix 
  for the given block.

  @param eid - element id
  @param ri  - row index (index of transported media)
  @param ci  - column index (index of transported media)

  @return The function returns integration order of the given element and block ids.

  Created by Tomas Koudelka, 11.11.2013
*/
long transtop::give_intordcm(long eid, long ri, long ci)
{
  elemtypet te;
  long ret = 0;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    ret = Lbt->intordcm[ri][ci];
    break;
  }
  case barlint3d:{
    ret = Lbt3d->intordcm[ri][ci];
    break;
  }
  case barlintax:{
    ret = Lbat->intordcm[ri][ci];
    break;
  }
  case barquadt:{
    ret = Qbt->intordcm[ri][ci];
    break;
  }
  case barquadtax:{
    ret = Qbat->intordcm[ri][ci];
    break;
  }
  case trlint:{
    ret = Ltt->intordcm[ri][ci];
    break;
  }
  case trlaxisym:{
    ret = Ltat->intordcm[ri][ci];
    break;
  }
  case quadlint:{
    ret = Lqt->intordcm[ri][ci];
    break;
  }
  case quadquadt:{
    ret = Qqt->intordcm[ri][ci];
    break;
  }
  case quadquadtax:{
    ret = Qqat->intordcm[ri][ci];
    break;
  }
  case quadlaxisym:{
    ret = Lqat->intordcm[ri][ci];
    break;
  }
  case lineartett:{
    ret = Ltett->intordcm[ri][ci];
    break;
  }
  case linearhext:{
    ret = Lht->intordcm[ri][ci];
    break;
  }
  case quadratichext:{
    ret = Qht->intordcm[ri][ci];
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return ret;
}



/**
   function returns degree of polynomial expansion
   
   @param eid - number of element
   
   TKr, 20.1.2003
   
*/
long transtop::give_degree (long eid)
{
  long deg;
  elemtypet te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case barlint:{         deg=1;  break;  }
  case barlint3d:{       deg=1;  break;  }
  case barlintax:{       deg=1;  break;  }
  case barquadt:{        deg=2;  break;  }
  case barquadtax:{      deg=2;  break;  }
  case trlint:{          deg=1;  break;  }
  case trlaxisym:{       deg=1;  break;  }
  case quadlint:{        deg=1;  break;  }
  case quadquadt:{       deg=2;  break;  }
  case quadquadtax:{     deg=2;  break;  }
  case quadlaxisym:{     deg=1;  break;  }
  case lineartett:{      deg=1;  break;  }
  case linearhext:{      deg=1;  break;  }
  case quadratichext:{   deg=2;  break;  }
  case linearwedget:{    deg=1;  break;  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return deg;
}

/**
   function returns dimension of problem
   
   @param eid - number of element
   
   JK, 23.2.2002
   modified by TKr 20.1.2033
*/
long transtop::give_dimension (long eid)
{
  long dim;
  elemtypet te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case barlint:{         dim=1;  break;  }
  case barlint3d:{       dim=1;  break;  }
  case barlintax:{       dim=1;  break;  }
  case barquadt:{        dim=1;  break;  }
  case barquadtax:{      dim=1;  break;  }
  case trlint:{          dim=2;  break;  }
  case trlaxisym:{       dim=2;  break;  }
  case quadlint:{        dim=2;  break;  }
  case quadquadt:{       dim=2;  break;  }
  case quadquadtax:{     dim=2;  break;  }
  case quadlaxisym:{     dim=2;  break;  }
  case ifacequadel:{     dim=2;  break;  }
  case lineartett:{      dim=3;  break;  }
  case linearhext:{      dim=3;  break;  }
  case quadratichext:{   dim=3;  break;  }
  case linearwedget:{    dim=3;  break;  }
  case gen2del:{         dim=2;  break;  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return dim;
}


/**
   function returns node numbers on required element edge
   
   @param eid - element id
   @param enid - end node id
   @param nodes - array containing end nodes
   
   JK, 15.11.2008
*/
void transtop::give_end_nodes (long eid, long enid, ivector &nodes)
{
  long nne;
  elemtypet te;
  ivector nod, endnod;
  
  //  type of element
  te = give_elem_type (eid);
  //  number of nodes on element
  nne = give_nne (eid);
  
  reallocv(RSTCKIVEC(nne, nod));
  reallocv(RSTCKIVEC(1, endnod));
  
  //  element nodes
  give_elemnodes (eid,nod);
  
  switch (te){
  case barlint:{
    endnod[0]=enid;
    break;
  }
  case barlint3d:{
    endnod[0]=enid;
    break;
  }
  case barlintax:{
    endnod[0]=enid;
    break;
  }
  case barquadt:{
    endnod[0]=enid;
    break;
  }
  case barquadtax:{
    endnod[0]=enid;
    break;
  }
  case trlint:{
    break;
  }
  case trlaxisym:{
    break;
  }
  case quadlint:{
    break;
  }
  case quadlaxisym:{
    break;
  }
  case quadquadt:{ 
    break;
  }
  case lineartett:{
    break;
  }
  case linearhext:{
    break;
  }
  case quadratichext:{
    break;
  }
  case linearwedget:{
    break;
  }
  case gen2del:{
    break; 
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  nodes[0]=nod[endnod[0]];
}


/**
   function returns node numbers on required element edge
   
   @param eid - element id
   @param edid - edge id
   @param nodes - array containing edge nodes
   
   5.3.2002
*/
void transtop::give_edge_nodes (long eid, long edid, ivector &nodes)
{
  long i,nne,nned;
  elemtypet te;
  ivector nod, edgenod;
  
  //  type of element
  te = give_elem_type (eid);
  //  number of nodes on element
  nne = give_nne (eid);
  //  number of nodes on one edge
  nned = give_nned (eid);
  
  reallocv(RSTCKIVEC(nne, nod));
  reallocv(RSTCKIVEC(nned, edgenod));
  
  //  element nodes
  give_elemnodes(eid, nod);
  
  switch (te){
  case barlint:{
    edgenod[0]=edid;
    break;
  }
  case barlint3d:{
    edgenod[0]=edid;
    break;
  }
  case barlintax:{
    break;
  }
  case barquadt:{
    break;
  }
  case barquadtax:{
    break;
  }
  case trlint:{
    lintriangle_edgnod (edgenod.a,edid);
    break;
  }
  case trlaxisym:{
    lintriangle_edgnod (edgenod.a,edid);
    break;
  }
  case quadlint:{
    linquadrilat_edgnod (edgenod.a,edid);
    break;
  }
  case quadlaxisym:{
    linquadrilat_edgnod (edgenod.a,edid);
    break;
  }
  case quadquadt:{ 
    quadquadrilat_edgnod (edgenod.a,edid);
    break;
  }
  case lineartett:{
    lintetrahedral_edgnod(edgenod.a,edid);
    break;
  }
  case linearhext:{
    linhexahedral_edgnod(edgenod.a,edid);
    break;
  }
  case quadratichext:{
    quadhexahedral_edgnod (edgenod.a,edid);
    break;
  }
  case gen2del:{
    break; 
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  for (i=0;i<nned;i++){
    nodes[i]=nod[edgenod[i]];
  }
}

/**
   function returns node numbers on required element surface
   
   @param eid - element id
   @param surfid - surface id
   @param nodes - array containing surface nodes
   
   JK, 24.8.2004
*/
void transtop::give_surface_nodes (long eid, long surfid, ivector &nodes)
{
  long i,nne,nnsurf;
  elemtypet te;
  ivector nod, surfnod;
  
  //  type of element
  te = give_elem_type (eid);
  //  number of nodes on element
  nne = give_nne (eid);
  //  number of nodes on one surface
  nnsurf = give_nnsurf (eid);
  
  reallocv(RSTCKIVEC(nne, nod));
  reallocv(RSTCKIVEC(nnsurf, surfnod));
  
  //  element nodes
  give_elemnodes(eid, nod);
  
  switch (te){
  case barlint:{
    break;
  }
  case barlint3d:{
    break;
  }
  case barlintax:{
    break;
  }
  case barquadt:{
    break;
  }
  case barquadtax:{
    break;
  }
  case trlint:{
    break;
  }
  case trlaxisym:{
    break;
  }
  case quadlint:{
    break;
  }
  case quadlaxisym:{
    break;
  }

  case lineartett:{
    break;
  }
  case linearhext:{
    linhexahedral_surfnod (surfnod.a,surfid);
    break;
  }
  case quadratichext:{
    quadhexahedral_surfnod (surfnod.a,surfid);
    break;
  }
  case gen2del:{
    break; 
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  for (i=0;i<nnsurf;i++){
    nodes[i]=nod[surfnod[i]];
  }
}



/**
   function determinates the number of boundary objects and number of nodes on them
   boundary objects are end nodes in 1D, edges in 2D and surfaces in 3D
   
   @param eid - element id
   @param bco - number of boundary objects
   @param ncompbco - number of components on one boundary object
   
   JK, 15.11.2008
*/
void transtop::give_nbobjects (long eid,long &bco,long &ncompbco)
{
  
  switch (elements[eid].te){
    //  1D elements
  case barlint:{
    bco=Lbt->nen;
    ncompbco=1;
    break;
  }
  case barlint3d:{
    bco=Lbt3d->nen;
    ncompbco=1;
    break;
  }
  case barlintax:{
    bco=Lbat->nen;
    ncompbco=1;
    break;
  }
  case barquadt:{
    bco=Qbt->nen;
    ncompbco=1;
    break;
  }
  case barquadtax:{
    bco=Qbat->nen;
    ncompbco=1;
    break;
  }
    
    //  2D elements
  case trlint:{
    bco=Ltt->ned;
    ncompbco=Ltt->nned;
    break;
  }
  case trlaxisym:{
    bco=Ltat->ned;
    ncompbco=Ltat->nned;
    break;
  }
  case quadlint:{
    bco=Lqt->ned;
    ncompbco=Lqt->nned;
    break;
  }
  case quadlaxisym:{
    bco=Lqat->ned;
    ncompbco=Lqat->nned;
    break;
  }
  case quadquadt:{
    bco=Qqt->ned;
    ncompbco=Qqt->nned;
    break;
  }
  case quadquadtax:{
    bco=Qqat->ned;
    ncompbco=Qqat->nned;
    break;
  }
   
    //  3D elements
  case lineartett:{
    bco=Ltett->nsurf;
    ncompbco=Ltett->nnsurf;
    break;
  }
  case linearhext:{
    bco=Lht->nsurf;
    ncompbco=Lht->nnsurf;
    break;
  }
  case quadratichext:{
    bco=Qht->nsurf;
    ncompbco=Qht->nnsurf;
    break;
  }
  case linearwedget:{
    bco=Lwt->nsurf;
    ncompbco=Lwt->nnsurf;
    break;
  }
  case gen2del:{
    bco=0;
    ncompbco=0;
    break; 
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function assembles nodes on boundary objects (end node, edge, surface)
   
   @param eid - element id
   @param boid - boundary object id
   @param nod - array of node numbers
   
   JK, 15.11.2008
*/
void transtop::give_bonodes (long eid, long boid, ivector &nod)
{
  
  switch (elements[eid].te){
    //  1D elements
  case barlint:
  case barlint3d:
  case barlintax:
  case barquadt:
  case barquadtax:{
    give_end_nodes(eid, boid, nod);
    break;
  }
    
    //  2D elements
  case trlint:
  case trlaxisym:
  case quadlint:
  case quadlaxisym:
  case quadquadt:
  case quadquadtax:{
    give_edge_nodes(eid, boid, nod);
    break;
  }
    
    //  3D elements
  case lineartett:
  case linearhext:
  case quadratichext:{
    give_surface_nodes(eid, boid, nod);
    break;
  }
  case gen2del:{
    break; 
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
}




void transtop::alloc_prep(long nn, long ne)
{
  nodes = new nodet[nn];
  elements = new elementt[ne];
}

/**
   function returns integral
   
   @param eid - number of element
   
*/
/*
double transtop::give_integral (long eid,vector &nodval)
{
  double value;
  elemtypet te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case barlint:{      value = Lbt->total_integral(eid,nodval);   break;}
  case barlint3d:{    value = Lbt3d->total_integral(eid,nodval);   break;}
  case barlintax:{    value = Lbat->total_integral(eid,nodval);  break;}
  case barquadt:{     value = Qbt->total_integral(eid,nodval);   break;}
  case barquadtax:{   value = Qbat->total_integral(eid,nodval);  break;}
  case quadlint:{     value = Lqt->total_integral(eid,nodval);   break;}
  case quadquadt:{    value = Qqt->total_integral(eid,nodval);   break;}
  case quadquadtax:{  value = Qqat->total_integral(eid,nodval);  break;}
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return value;
}
*/

/**
   function allocates all needed arrays on nodes
   
   JK, 10.7.2008 - revision
*/
void transtop::alloc_nodes (void)
{
  long i,j,nne,ncomp, ncompo, ncompeqo, ipp;
  ivector nod;
  
  //  loop over elements
  for (i=0;i<ne;i++){
    //  the number of nodes on element
    nne = give_nne (i);
    reallocv (RSTCKIVEC(nne,nod));
    //  nodes on element
    give_elemnodes (i,nod);
    //  the number of components of gradients and fluxes
    ncomp=give_dimension (i);
    
    ncompo=Tm->givencompother ();
    
    //  integration point id
    ipp = elements[i].ipp[0][0];
    //  the number of components is the array eqother
    ncompeqo=Tm->givencompeqother (ipp,0);
    
    //  loop over the number of nodes on element
    for (j=0;j<nne;j++){
      //  array of gradients
      if (Tp->gradcomp==1){
	if (nodes[nod[j]].gradient==NULL){
	  nodes[nod[j]].alloc_grad (ncomp);
	}
      }
      //  array of fluxes
      if (Tp->fluxcomp==1){
	if (nodes[nod[j]].flux==NULL){
	  nodes[nod[j]].alloc_flux (ncomp);
	}
      }
      //  array of other values
      if (Tp->othercomp==1){
	if (nodes[nod[j]].other==NULL){
	  nodes[nod[j]].alloc_other (ncompo);
	}
      }
      //  array of eqother values
      if (Tp->eqothercomp==1){
	if (nodes[nod[j]].eqother==NULL){
	  nodes[nod[j]].alloc_eqother (ncompeqo);
	}
      }
    }//  end of the loop for (j=0;j<nne;j++){
  }//  end of the loop for (i=0;i<ne;i++){
  
  // There may be unattended nodes (due to simplified passing data in METR) 
  // which are not handled in the previous loop
  // take number of components from the first element
  ncomp  = give_dimension(0);
  ncompo = Tm->givencompother();
  ipp = elements[0].ipp[0][0];
  ncompeqo=Tm->givencompeqother (ipp,0);
  // loop over all nodes
  for (i=0; i<nn; i++){
    if (Tp->gradcomp==1){
      if (nodes[i].gradient==NULL){
        nodes[i].alloc_grad (ncomp);
      }
    }
    //  array of fluxes
    if (Tp->fluxcomp==1){
      if (nodes[i].flux==NULL){
        nodes[i].alloc_flux (ncomp);
      }
    }
    //  array of other values
    if (Tp->othercomp==1){
      if (nodes[i].other==NULL){
        nodes[i].alloc_other (ncompo);
      }
    }
    //  array of eqother values
    if (Tp->eqothercomp==1){
      if (nodes[i].eqother==NULL){
        nodes[i].alloc_eqother (ncompeqo);
      }
    }    
  }

  //  actual nodal values
  if (Tp->nvs==1){
    for (i=0;i<nn;i++){
      nodes[i].alloc_nodval ();
    }
  }
  //  previous nodal values
  if (Tp->pnvs==1){
    for (i=0;i<nn;i++){
      nodes[i].alloc_nodvalp ();
    }
  }
  //  initial nodal values
  if (Tp->invs==1){
    for (i=0;i<nn;i++){
      nodes[i].alloc_nodvali ();
    }
  }
  //  time derivatives of actual nodal values
  if (Tp->tdnvs==1){
    for (i=0;i<nn;i++){
      nodes[i].alloc_nodvalt ();
    }
  }
}

/**
   functions allocates end nodes
   end nodes are allocated in problems with discontinuities
   
   JK, 17.10.2008
*/
void transtop::alloc_enodes ()
{
  if (Gtt->nen>0){
    enodes = new endnodet [Gtt->nen];
  }
}

/**
   functions allocates edges
   edges are allocated in problems with discontinuities
   
   JK, 8.8.2008
*/
void transtop::alloc_edges ()
{
  if (Gtt->nged>0){
    edges = new edget [Gtt->nged];
  }
}

/**
   function initializes end node parameters and arrays
   
   JK, 17.10.2008
*/
void transtop::enodes_init ()
{
  long i;
  
  for (i=0;i<Gtt->nen;i++){
    enodes[i].init (i);
  }
}

/**
   function initializes edge parameters and arrays
   
   JK, 8.8.2008
*/
void transtop::edge_init ()
{
  long i;
  
  for (i=0;i<Gtt->nged;i++){
    edges[i].init (i);
  }
}

/**
   function initializes edge array edval
   this function is used in problems with radiation
   
   JK, 26.7.2011
*/
void transtop::edge_init_edval ()
{
  long i;
  
  for (i=0;i<Gtt->nged;i++){
    edges[i].init_edval ();
  }
}


/**
   function computes jumps at discontinuity
   jumps are saved on end nodes and edges
   
   @param rhs - right hand side %vector
   
   JK, 8.8.2008
*/
void transtop::compute_jumps (double *rhs)
{
  long i,j,nm,*ncn1,*ncn2,*mcn;
  
  //  loop over end nodes
  for (i=0;i<Gtt->nen;i++){
    enodes[i].compute_jump (i);
    
    //  number of multipliers between first nodes
    nm=Gtt->endnodes[i].ndofn;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    
    
    Gtt->give_endnode_code_numbers (i,ncn1,ncn2,mcn);
    
    for (j=0;j<nm;j++){
      if (mcn[j]>0)
	rhs[mcn[j]-1]+=enodes[i].jump[j];
    }
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }

  //  loop over edges
  for (i=0;i<Gtt->nged;i++){
    edges[i].compute_jump (i);
    
    //  number of multipliers between first nodes
    nm=Gtt->gedges[i].ndofnf;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    
    
    Gtt->give_edge_code_numbers (i,1,ncn1,ncn2,mcn);
    
    for (j=0;j<nm;j++){
      if (mcn[j]>0)
	rhs[mcn[j]-1]+=edges[i].jumpfn[j];
    }
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;


    //  number of multipliers between last nodes
    nm=Gtt->gedges[i].ndofnl;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    
    
    Gtt->give_edge_code_numbers (i,2,ncn1,ncn2,mcn);
    
    for (j=0;j<nm;j++){
      if (mcn[j]>0)
	rhs[mcn[j]-1]+=edges[i].jumpln[j];
    }
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;

  }
  
}


/**
   function computes resistance factors on material interface
   
   JK, 2.2.2011
*/
void transtop::compute_resistance_factor (double *rhs)
{
  long i,j,k,edid,fnid,lnid,nm,lcid,*ncn1,*ncn2,*mcn;
  double ncf,*normal,*flux;
  
  //  determine lcid properly!
  lcid=0;
  
  //  provizorne alokace pro 2D ulohy
  normal = new double [2];
  flux = new double [2];

  //  nodal fluxes are computed in all nodes
  //  later on, only fluxes on interface should be computed
  compute_nodefluxes ();


  // ********************
  //  end nodes in 1D  **
  // ********************

  /*
  //  loop over end nodes
  for (i=0;i<Gtt->nen;i++){
    enodes[i].compute_jump (i);
    
    //  number of multipliers between first nodes
    nm=Gtt->endnodes[i].ndofn;
    ncn1 = new long [nm];
    ncn2 = new long [nm];
    mcn = new long [nm];
    
    
    Gtt->give_endnode_code_numbers (i,ncn1,ncn2,mcn);
    
    for (j=0;j<nm;j++){
      if (mcn[j]>0)
	rhs[mcn[j]-1]+=enodes[i].jump[j];
    }
    delete [] mcn;
    delete [] ncn2;
    delete [] ncn1;
  }
  */
  

  // *********
  //  edges **
  // *********

  //  loop over the number of series
  for (i=0;i<Gtt->nser;i++){
    //  loop over the number of edges in particular series
    for (j=0;j<Gtt->nedser[i];j++){
      //  edge id
      edid = Gtt->edgelist[i][j];
      
      //  first node
      fnid = Gtt->gedges[edid].fn;
      
      //  normal vector to the edge
      Gtt->gedges[edid].give_norvect (normal);
      
      //  vector of nodal fluxes
      //  zde je treba nastavit spravne lcid
      nodes[fnid].give_flux (lcid,flux);
      
      //  normal component of flux
      //  provizorne jen pro 2D
      ncf = ss (normal,flux,2);
      
      //  computation of the resistance factor and jumps across the material interface
      //  zde je treba nastavit spravne lcid
      edges[edid].compute_node_jump (edid,1,lcid,ncf);
      
      //  number of multipliers between first nodes
      nm=Gtt->gedges[edid].ndofnf;
      ncn1 = new long [nm];
      ncn2 = new long [nm];
      mcn = new long [nm];
      
      
      //  code numbers of unknowns located on edge
      Gtt->give_edge_code_numbers (edid,1,ncn1,ncn2,mcn);
      
      for (k=0;k<nm;k++){
	if (mcn[k]>0)
	  rhs[mcn[k]-1]+=edges[edid].jumpfn[k];
      }
      delete [] mcn;
      delete [] ncn2;
      delete [] ncn1;
      


      if (j==Gtt->nedser[i]){
	
	//  last node
	lnid = Gtt->gedges[edid].ln;
	
	//  vector of nodal fluxes
	//  zde je treba nastavit spravne lcid
	nodes[lnid].give_flux (lcid,flux);

	//  normal component of flux
	//  provizorne jen pro 2D
	ncf = ss (normal,flux,2);
	
	
	//  computation of the resistance factor and jumps across the material interface
	//  zde je treba nastavit spravne lcid
	edges[edid].compute_node_jump (edid,2,lcid,ncf);
	
	//  number of multipliers between first nodes
	nm=Gtt->gedges[edid].ndofnf;
	ncn1 = new long [nm];
	ncn2 = new long [nm];
	mcn = new long [nm];
	
	//  code numbers of unknowns located on edge
	Gtt->give_edge_code_numbers (edid,2,ncn1,ncn2,mcn);
	
	for (k=0;k<nm;k++){
	  if (mcn[k]>0)
	    rhs[mcn[k]-1]+=edges[edid].jumpln[k];
	}
	delete [] mcn;
	delete [] ncn2;
	delete [] ncn1;
	
      }
    }
  }
  
  delete [] flux;
}




/**
   function establishes %vector of initial nodal values
   
   JK, 5.3.2006
*/
void transtop::initial_nodval ()
{
  long i,j,k,ndofe;
  vector r;
  
  for (i=0;i<ne;i++){
    
    j=Gtt->leso[i];
    k=Gtt->gelements[i].auxinf;
    
    if (j==1 && k==0){
      //  added element
      //  number of DOFs on element
      ndofe = give_ndofe (i);
      reallocv(RSTCKVEC(ndofe, r));

      //  nodal values on element
      elemvalues(i, r);
      //  definition of initial nodal values
      elements[i].initnodvalues(r);
      
      Gtt->gelements[i].auxinf=1;
    }
  }
}



/**
   function saves nodal unknowns from the array lhs to nodes
   function is used in problems with growing number of nodes and elements
   
   @param lhs - array of nodal values
   @param lhsi - array of initial nodal values
   @param tdlhs - array of time derivatives of nodal values
   
   JK, 2.6.2006
*/
void transtop::lhs_save (double *lhs,double *lhsi,double *tdlhs)
{
  long i,j,ndofn;
  long *cn;

  for (i=0;i<nn;i++){
    //  number of DOFs on node
    ndofn = give_ndofn (i);
    cn = new long [ndofn];
    give_node_code_numbers (i,cn);
    
    for (j=0;j<ndofn;j++){
      nodes[i].nodval[j]=0.0;
      nodes[i].nodvali[j]=0.0;
      nodes[i].nodvalt[j]=0.0;
      if (cn[j]>0){
	nodes[i].nodval[j]=lhs[cn[j]-1];
	nodes[i].nodvali[j]=lhsi[cn[j]-1];
	nodes[i].nodvalt[j]=tdlhs[cn[j]-1];
      }
    }   
    delete [] cn;
  }
  
}

/**
   function restores nodal unknowns from nodes to the array lhs
   function is used in problems with growing number of nodes and elements
   
   @param lhs - array of nodal values
   @param lhsi - array of initial nodal values
   @param tdlhs - array of time derivatives of nodal values

   JK, 2.6.2006
*/
void transtop::lhs_restore (double *lhs,double *lhsi,double *tdlhs)
{
  long i,j,ndofn;
  long *cn;

  for (i=0;i<nn;i++){
    ndofn = give_ndofn (i);
    cn = new long [ndofn];
    give_node_code_numbers (i,cn);
    
    for (j=0;j<ndofn;j++){
      if (cn[j]>0){
	lhs[cn[j]-1]=nodes[i].nodval[j];
	lhsi[cn[j]-1]=nodes[i].nodvali[j];
	tdlhs[cn[j]-1]=nodes[i].nodvalt[j];
      }
    }
    delete [] cn;
  }
  
}


long transtop::mesh_check(void)
{
  long i,j,err=0,kk,jj;

  //elems and nodes check
  for(i=0;i<ne;i++){
    for(j=0;j<Gtt->gelements[i].nne;j++){
      
      jj = Gtt->gelements[i].nodes[j];
      for(kk = 0; kk < Tp->ntm;kk++){
	
	if((Gtt->gelements[i].tgf) < (Gtt->gnodes[jj].tgf[kk])){

	  fprintf (stdout,"\n\n Different element node(%ld) time function number=%ld than element(%ld) time function number=%ld",j+1,Gtt->gnodes[jj].tgf[kk]+1,i+1,Gtt->gelements[i].tgf+1);
	  fprintf (stdout,"\n Node number=%ld, element number=%ld\n\n",Gtt->gelements[i].nodes[j]+1,i+1);
	  err = 1;

	  return err;
	}
	
      }
    }
  }
  
  return err;
}

/**
   function defines types of materials on nodes
   this is necessary for problems with discontinuities
   
   JK, 16.11.2007
*/
/*
void transtop::node_materials ()
{
  long i,j,mid,nne;
  ivector nodes;
  
  for (i=0;i<ne;i++){
    //  number nodes on element
    nne = give_nne (i);
    //  allocation of array nodes
    allocv (nne,nodes);
    //  element nodes
    give_elemnodes (i,nodes);
    //  id of material defined on element
    mid = elements[i].idm[0];
    for (j=0;j<nne;j++){
      Gtt->gnodes[nodes[j]].ai=mid;
    }
    destrv (nodes);
  }
}
*/

/*
void transtop::jump_initiation (long nbn,seqselnodes *selnodfeti)
{
  long i,j,k,ns,nsnsd;
  
  //  total number of interface/boundary nodes
  tnbn = nbn;

  //  jnmap[i][0]=j - the i-th interface/boundary node has local number j on the subdomain with lower id
  //  jnmap[i][1]=j - the i-th interface/boundary node has local number j on the subdomain with higher id
  jnmap = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    jnmap[i] = new long [2];
    
    jnmap[i][0] = -1;
    jnmap[i][1] = -1;
  }
  
  //  number of subdomains
  ns = Gtt->stop->ns;
  for (i=0;i<ns;i++){
    //  number of selected nodes on subdomain
    nsnsd = selnodfeti->nsndom[i];
    for (j=0;j<nsnsd;j++){
      k=selnodfeti->gnn[i][j];
      if (jnmap[k][0]==-1)
	jnmap[k][0] = selnodfeti->lsnl[i][j];
      else
	jnmap[k][1] = selnodfeti->lsnl[i][j];
    }
  }
}
*/

/**
   function computes view factors which are used in problems with radiation
   
   @param out - output file
   
   Jan Koci, 12.7.2011
*/
void transtop::view_factors (FILE */*out*/)
{
  long i, j, k;
  long emitter, emitter_fn, emitter_ln;
  long absorber, absorber_fn, absorber_ln;
  double diagonal1, diagonal2, side1, side2;
  double emitter_size,emitter_fn_Xcord, emitter_fn_Ycord, emitter_ln_Xcord, emitter_ln_Ycord;
  double absorber_fn_Xcord, absorber_fn_Ycord, absorber_ln_Xcord, absorber_ln_Ycord;
  //double view_factor[3][100][100];  //nutno naalokovat dynamicky
  //double ***view_factor;  //nutno naalokovat dynamicky
  
  view_factor = new double** [Gtt->nser];
  for (i=0;i<Gtt->nser;i++){
    view_factor[i] = new double* [Gtt->nedser[i]];
    for (j=0;j<Gtt->nedser[i];j++){
      view_factor[i][j] = new double [Gtt->nedser[i]];
    }
  }
  
  //  loop over the number of series/cavities
  for (i=0;i<Gtt->nser; i++){
    //  loop over the number of edges in the series
    for (j=0; j < Gtt->nedser[i] ; j++){
      
      emitter = Gtt->edgelist[i][j];         //vyzarujici hrana
      
      emitter_fn = Gtt->gedges[emitter].fn;  //cislo pocatecniho uzlu vyzarujici hrany
      emitter_ln = Gtt->gedges[emitter].ln;  //cislo koncoveho uzlu vyzarujici hrany
      
      //souradnice uzlu vyzarujici hrany
      emitter_fn_Xcord = Gtt->gnodes[emitter_fn].x;
      emitter_fn_Ycord = Gtt->gnodes[emitter_fn].y;
      emitter_ln_Xcord = Gtt->gnodes[emitter_ln].x;
      emitter_ln_Ycord = Gtt->gnodes[emitter_ln].y;
      
      //  length of the egde
      emitter_size = sqrt(pow(emitter_fn_Xcord-emitter_ln_Xcord,2)+pow(emitter_fn_Ycord-emitter_ln_Ycord,2));
      //  the lenght is stored in the object gedges
      Gtt->gedges[emitter].l = emitter_size;
      
      //  loop over the number of edges in the series
      for (k=0; k < Gtt->nedser[i]; k++) {
	absorber = Gtt->edgelist[i][k];         //pohlcujici hrana
	
	absorber_fn = Gtt->gedges[absorber].fn;  //cislo pocatecniho uzlu pohlcujici hrany
	absorber_ln = Gtt->gedges[absorber].ln;  //cislo koncoveho uzlu pohlcujici hrany
	
	//souradnice uzlu pohlcujici hrany
	absorber_fn_Xcord = Gtt->gnodes[absorber_fn].x;
	absorber_fn_Ycord = Gtt->gnodes[absorber_fn].y;
	absorber_ln_Xcord = Gtt->gnodes[absorber_ln].x;
	absorber_ln_Ycord = Gtt->gnodes[absorber_ln].y;
	
	//vypocet VIEW FACTORU
	if (j == k) {
	  view_factor[i][j][k] = 0;
	} else {
	  diagonal1 = sqrt(pow(emitter_fn_Xcord-absorber_fn_Xcord,2)+pow(emitter_fn_Ycord-absorber_fn_Ycord,2));
	  diagonal2 = sqrt(pow(emitter_ln_Xcord-absorber_ln_Xcord,2)+pow(emitter_ln_Ycord-absorber_ln_Ycord,2));
	  side1 = sqrt(pow(emitter_fn_Xcord-absorber_ln_Xcord,2)+pow(emitter_fn_Ycord-absorber_ln_Ycord,2));
	  side2 = sqrt(pow(emitter_ln_Xcord-absorber_fn_Xcord,2)+pow(emitter_ln_Ycord-absorber_fn_Ycord,2));
	  
	  view_factor[i][j][k] = (diagonal1 + diagonal2 - side1 - side2)/(2*emitter_size);
	}
      }
    }
  }
  
  /*
  for (i=0;i<3;i++){
    fprintf (out,"\n\n\n sada %ld \n",i);
    for (j=0;j<16;j++){
      for (k=0;k<16;k++){
	fprintf (out," %le",view_factor[i][j][k]);
      }
      fprintf (out,"\n");
    }
  }
  */
  
  
}

/**
   function computes view factors which are used in problems with radiation
   
   @param out - output file
   
   Jan Koci, 26.7.2011
*/
void transtop::mod_view_factors (FILE */*out*/)
{
  long i,j,k,n,edid,e1id,e2id,matid;
  double eps,zero;
  double *x,*y;
  mattypet mat1,mat2;
  x=NULL;
  y=NULL;

  zero=1.0e-15;
  
  mvf = new densemat [Gtt->nser];
  mvfjk = new densemat [Gtt->nser];
  
  //  loop over the number of series/cavities
  for (i=0;i<Gtt->nser; i++){
    n=Gtt->nedser[i];
    mvf[i].alloc (n);
    mvfjk[i].alloc (n);
    
    //  loop over the number of edges in the series
    for (j=0;j<Gtt->nedser[i];j++){

      //  loop over the number of edges in the series
      for (k=0; k < Gtt->nedser[i]; k++) {
	//  edge number
	edid = Gtt->edgelist[i][k];
	//  element id
	e1id = Gtt->gedges[edid].adjel[0];
	e2id = Gtt->gedges[edid].adjel[1];
	
	mat1 = elements[e1id].tm[0];
	mat2 = elements[e2id].tm[0];
	
	if (mat1 == radiationmater){
	  //  material id
	  matid = elements[e2id].idm[0];
	  //  extinction coefficient for the k-th edge
	  eps = Tm->give_extinction_coeff (mat2,matid);
	}else{
	  //  material id
	  matid = elements[e1id].idm[0];
	  //  extinction coefficient for the k-th edge
	  eps = Tm->give_extinction_coeff (mat1,matid);
	}
	
	
	if (j==k){
	  mvf[i].a[j*n+k]=1.0/eps;
	  mvfjk[i].a[j*n+k]=1.0/eps;
	}else{
	  mvf[i].a[j*n+k]=(1.0-1.0/eps)*view_factor[i][j][k];
	  mvfjk[i].a[j*n+k]=(1.0-1.0/eps)*view_factor[i][j][k];
	}
      }
    }
    
    //  modified view factor matrices are factorized
    mvf[i].lu (x,y,zero,2);
  }
  
  /*
  for (i=0;i<3;i++){
    fprintf (out,"\n\n\n sada %ld \n",i);
    for (j=0;j<16;j++){
      for (k=0;k<16;k++){
	fprintf (out," %le",view_factor[i][j][k]);
      }
      fprintf (out,"\n");
    }
  }
  */
  
  
}


/**
   function computes edge temperatures and stores it in the objects edget
   
   @param out - output file
   
   JK, 25.7.2011
*/
void transtop::edge_temperature ()
{
  long i,j,edid,ndofn,fn,ln,lcid;
  double tfn,tln,taver;
  vector nval;
  
  //  load case id
  //  it must be zero, all prescribed values are stored in the load case with id=0
  lcid=0;
  
  //  loop over the number of series/cavities
  for (i=0;i<Gtt->nser;i++){
    //  loop over the number of edges in the series
    for (j=0;j<Gtt->nedser[i];j++){
      //  edge number
      edid = Gtt->edgelist[i][j];
      //  first node on the actual edge
      fn = Gtt->gedges[edid].fn;
      //  last node on the actual edge
      ln = Gtt->gedges[edid].ln;
      //  the number of DOFs in nodes
      ndofn = give_ndofn (fn);
      //  array for nodal values
      reallocv(RSTCKVEC(ndofn, nval));
      
      //  nodal values in first node
      nodalval(fn, nval);
      //  temperature in first node
      tfn = nval[0];
      //  nodal values in last node
      nodalval(ln, nval);
      //  temperature in last node
      tln = nval[0];
      
      //  average temperature on the edge
      taver = (tfn+tln)/2.0;
      //  the average temperature is stored to the edge
      edges[edid].store_edval (&taver);
    }//  end of the loop over the number of edges in the series
  }//  end of the loop over the number of series/cavities
}


/**
   function computes heat fluxes in problems with radiation
   
   @param out - output file
   
   Jan Koci, 25.7.2011
*/
void transtop::heat_fluxes (double *rhs,FILE */*out*/)
{
  long i,j,edid,dof,fn,ln;
  double zero,l,h;
  double *x,*y,*z;
  
  zero=1.0e-15;
  
  //  loop over the number of series/cavities
  for (i=0;i<Gtt->nser;i++){
    //  array for heat fluxes
    x = new double [Gtt->nedser[i]];
    //  array for the right hand side
    y = new double [Gtt->nedser[i]];
    z = new double [Gtt->nedser[i]];
    
    //  vypocet prave strany
    t4t4 (i,y);
    
    mvf[i].lu (x,y,zero,3);
    //mvf[i].gemp (x,y,1,zero,1);
    mvfjk[i].mxv_dm (x,z);
    
    for (j=0;j<Gtt->nedser[i];j++){
      //  edge number
      edid = Gtt->edgelist[i][j];
      //  the lenght of the edge
      l=Gtt->gedges[edid].l;
      //  amount of heat connected with each node
      h=x[j]*l/2.0;
      //  first node on the actual edge
      fn = Gtt->gedges[edid].fn;
      //  last node on the actual edge
      ln = Gtt->gedges[edid].ln;

      dof = give_dof (fn,0)-1;
      if (dof>-1){
	rhs[dof]-=h;
      }
      dof = give_dof (ln,0)-1;
      if (dof>-1){
	rhs[dof]-=h;
      }
      

    }
    
    delete [] x;
    delete [] y;
  }//  end of the loop over the number of series/cavities
}

/**
   function computes heat fluxes for selected series/cavity
   
   @param edserid - id of edge series
   @param y - array containing heat fluxes for the edserid-th series/cavity
   
   JK, 25.7.2011
*/
void transtop::t4t4 (long edserid,double *y)
{
  long i,j,edid;
  double stefan_boltzmann,ti,tj;
  
  //  Stefan-Boltzmann constant
  stefan_boltzmann=5.67e-8;
  
  //  loop over the number of edges in the actual series
  for (i=0;i<Gtt->nedser[edserid];i++){
    //  edge number
    edid = Gtt->edgelist[edserid][i];
    //  temperature on the i-the edge
    edges[edid].give_edval (&ti);
    
    y[i]=0.0;
    //  loop over the number of edges in the actual series
    //  all edges in the series are taken into account for each i-th edge
    for (j=0;j<Gtt->nedser[edserid];j++){
      //  edge number
      edid = Gtt->edgelist[edserid][j];
      //  temperature on the j-the edge
      edges[edid].give_edval (&tj);
      
      y[i]+=view_factor[edserid][i][j]*stefan_boltzmann*(ti*ti*ti*ti - tj*tj*tj*tj);
    }//  end of the loop over the number of edges in the actual series
  }//  end of the loop over the number of edges in the actual series
}



/** 
  The function searches on the given element eid for the closest integration point 
  to the point of given coordinates [px,py,pz]. It retruns the closest integration point 
  id and its natural coordinates.

  @param eid - element id whose integration points are investigated (input)
  @param px, py, pz - global coordinates of the given point (input)
  @param ipm - array of mapping stuctures where the global coordinates of integration points are stored (input)
               ipm[ipp] = ipmap structure for the TRFEL integration point ipp (input)
  @param iptol - tolerance for the distance between [px,py,pz] and integration point below which the 
                 points are assumed to be identical (input)
  @param xi, eta, zeta - natural coordinates of the closest integration point (output)

  @retval ipp <  0 - if the distance to closest integration point is out of tolerance ipctol, 
  @retval ipp >= 0 - if the point [px,py,pz] is assumed to be identical with the integration point ipp

  Created by Tomas Koudelka, 2.12.2016
*/  
long transtop::give_closest_ip_ncoord(long eid, double px, double py, double pz, 
                                      ipmap *ipm, double iptol, double &xi, double &eta, double &zeta)
{
  long i, nip = give_tnip(eid);
  long ipp = elements[eid].ipp[0][0];
  long ipc = ipp;
  double dmin;
  double d;
  vector ncoord(ASTCKVEC(3));

  for(i=0; i<nip; ipp++, i++)
  {
    d = sqr(px - ipm[ipp].x) + sqr(py - ipm[ipp].y) + sqr(pz - ipm[ipp].z);
    if (i == 0)
    {
      dmin = d;
      continue;
    }
    if (d < dmin)
    {
      dmin = d;
      ipc = ipp;
    }
  }
  dmin = sqrt(dmin);
  // get the natural coordinates of the closest integration point
  ipncoordt(eid, ipc, ncoord);
  xi   = ncoord[0];
  eta  = ncoord[1];
  zeta = ncoord[2];

  if (dmin <= iptol)
    return ipc;
  
  return -1;
}


void transtop::sort_elements()
{
  long *delta = new long[ne];
  long *mindof = new long[ne];
  ivector cn;
  long i, j, k;
  long max, min;
  long ndofe, d;

  if (renel==NULL)
    renel = new long[ne];

  for (i=0; i<ne; i++)
  {
    // find maximum and minimum DOF number of the i-th element
    ndofe = Gtt->give_ndofe(i);
    reallocv(RSTCKIVEC(ndofe, cn));
    Gtt->give_code_numbers (i, cn.a);
    max = 0;
    min = Ndoft+1;
    for (j=0; j<ndofe; j++)
    {
      if (cn[j] > max)
        max = cn[j];

      if ((cn[j] < min) && (cn[j] > 0))
        min = cn[j];
    }

    // calculate delta of the maximum and minimum DOF numbers
    mindof[i] = min;
    delta[i]  = max-min; 
    renel[i] = i;
  }

  long tmp1, tmp2, tmp3;
  d = ne;

  // sort elements according to minimum DOF number
  long gap = ne / 2;
  long l, n;
  while (gap > 0) 
  { 
    for (i = 0; i < ne - gap; i++) 
    {
      j = i + gap;
      tmp1 = mindof[j];
      tmp2 = delta[j];
      tmp3 = renel[j];
      while ((j >= gap) && (tmp1 < mindof[j-gap]))
      {
        mindof[j] = mindof[j - gap];
        delta[j] = delta[j - gap];
        renel[j] = renel[j - gap];
        j -= gap;
      }
      mindof[j] = tmp1;
      delta[j] = tmp2;
      renel[j] = tmp3;
    }
    if (gap == 2) {
      gap = 1;
    } else {
      gap = long(double(gap)/2.2);
    }
  }

  // resort elements with the same min DOF number according to delta array values
  tmp1 = mindof[0];
  j = 0;
  for (l=0; l<ne; l++)
  {
    if (tmp1 == mindof[l])
      continue;
    else
    {  
      n = l - j;
      gap = n / 2;
      while (gap > 0) 
      { 
        for (i = j; i < l - gap; i++) 
        {
          k = i + gap;
          tmp1 = mindof[k];
          tmp2 = delta[k];
          tmp3 = renel[k];
          while ((k >= gap+j) && (tmp2 < delta[k-gap]))
          {
            mindof[k] = mindof[k-gap];
            delta[k] = delta[k-gap];
            renel[k] = renel[k-gap];
            k -= gap;
          }
          mindof[k] = tmp1;
          delta[k] = tmp2;
          renel[k] = tmp3;
        }
        if (gap == 2) { // change the gap size
          gap = 1;
        } else {
          gap = long(double(gap)/2.2);
        }
      }
      j = l;
      tmp1 = mindof[l];
    }
  }
  memset(delta, 0, sizeof(*delta)*ne);
  for (i=0; i<ne; i++)
    delta[renel[i]] = 1;
  for (i=0, j=0; i<ne; i++)
    j += delta[i];
  fprintf(stdout, "\n total number of reordered elements = %ld", j);

  delete [] delta;
  delete [] mindof;
}



/**
  The function assembles vectors of weights of masternodes for all hanging nodes.

  @return The function does not return anything but it allocates and initializes %vector mnw at hanging nodes.

  Created by TKo, 23.7.2018
*/
void transtop::assemble_master_node_weights_vec(void)
{
  long i, ndofn, nmn;

  for(i=0; i<nn; i++){
    //  the number of DOFs in node
    ndofn = Gtt->give_ndofn(i);

    if (ndofn >= 0){
      // a regular node case
      continue;
    }
    else{
      // a hanging node case
      
      // number of master nodes
      nmn = 0-ndofn;
      nodes[i].mnw = new vector(nmn);    
      // compute master node weights
      Gtt->approx_weights(Gtt->gnodes[i].masentity, i, *nodes[i].mnw);
    }
  }
}



/**
  The function searches for hanging nodes and assembles transformation matrices with master 
  node weight coefficients on elements. Vectors of master node weights MUST be initialized at 
  the hanging nodes (call assemble_master_node_weights_vec function).

  @return The function does not return anything but stores resulting trans

  Created by Tomas Koudelka according to JK, 17.7.2018
*/
void transtop::searching_hanging_nodes(void)
{
  long i,j,k,l,nmn,npn,nneo,nnen,ndofeo,ndofen,ndofn,ri,ci,cumulci;
  ivector oldnodes,newnodes;
  matrix tmat;
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    if (give_ndofe(i)!=Gtt->give_ndofe(i)){
      //  element was modified because of the hanging nodes
      
      //  the original number of nodes on element
      nneo = give_nne (i);
      //  the new number of nodes on element
      nnen = Gtt->give_nmne (i);
      
      reallocv (RSTCKIVEC(nneo,oldnodes));
      reallocv (RSTCKIVEC(nnen,newnodes));
      
      //  the original nodes on element
      Gtt->give_original_nodes (i,oldnodes);
      //  the new nodes on element
      Gtt->give_master_nodes (i,newnodes);
      
      //  the original number of DOFs on element
      ndofeo = give_ndofe (i);
      //  the new number of DOFs on element
      ndofen = Gtt->give_ndofe (i);
      
      //  transformation matrix
      reallocm (RSTCKMAT(ndofeo,ndofen,tmat));
      
      //  the number of processed nodes
      npn=0;
      //  row index
      ri=0;
      cumulci=0;
      //  loop over the original number of nodes
      for (j=0;j<nneo;j++){
	//  column index
	ci=cumulci;
	
	//  the number of DOFs in node
	ndofn = give_ndofn (oldnodes[j]);
	
	if (oldnodes[j]==newnodes[npn]){
	  //  the j-th original node is not hanging
	  
	  //  loop over the number of DOFs in actual node
	  for (k=0;k<ndofn;k++){
	    tmat[ri][ci]=1.0;
	    ri++;
	    ci++;
	    cumulci++;
	  }//  end of the loop over the number of DOFs in actual node
	  npn++;
	}//  end of the if statement (oldnodes[j]==newnodes[j])
	else{
	  //  the j-th original node is hanging
	  
	  //  the number of master nodes for the actual hanging node
	  nmn=0-Gtt->gnodes[oldnodes[j]].ndofn;	  
	  //  the number of DOFs in node
	  ndofn = give_ndofn (newnodes[npn]);

	  //  loop over the number of DOFs in actual node
	  for (k=0;k<ndofn;k++){
	    ci=cumulci+k;//ri;
	    
	    //  loop over the number of master nodes
	    for (l=0;l<nmn;l++){
	      //tmat[ri][ci]=1.0;
	      tmat[ri][ci] = (*nodes[oldnodes[j]].mnw)[l];
	      ci+=ndofn;
	    }//  end of the loop over the number of master nodes
	    ri++;
	  }//  end of the loop over the number of DOFs in actual node
	  npn+=nmn;
	  cumulci=ci+1-ndofn;
	}//  end of the else
	
      }//  end of the loop over the original number of nodes
      
      elements[i].tmat = new matrix (ndofeo,ndofen);
      copym (tmat,*(elements[i].tmat));
    }
  }//  end of the loop over the number of elements
}
