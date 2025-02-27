#include "couptop.h"
#include "globalc.h"
#include "mechtop.h"
#include "gtopology.h"
#include "element.h"
#include "elementc.h"
#include "nodec.h"


couptop::couptop (void)
{
  ne=0;
  elements=NULL;
}

couptop::~couptop (void)
{
  delete [] elements;
}

void couptop::read (XFILE *in)
{
  long i,j;
  
  if ((Cp->tprob==par_coupl_mech_trans) || (Cp->tprob==growing_par_coupl_mech_trans)){
    
  }

  if (Cp->tprob==fully_coupled_mech_trans){
    //  element data
    xfscanf (in,"%ld",&ne);
    
    if (Mesprt==1)  fprintf (stdout,"\n number of elements  %ld",ne);
    elements = new elementc [ne];
    
    for (i=0;i<ne;i++){
      j=ne+1;
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("element number %ld is less than 1", __FILE__, __LINE__, __func__, j);
      if (j>ne){
	print_err("element number %ld is greater than total number of elements %ld", __FILE__, __LINE__, __func__, j, ne);
      }
      elements[j-1].read (in,j-1);
    }
  }
  
  if (Cp->tprob==fully_coupled_material){
    //  this is fully coupled material model
    //  fully coupled analysis is required
    //  in this case, one mesh is read
    //  the mesh is the same for mechanics as well as for transport
    
    //  the number of nodes
    xfscanf (in,"%ld",&nn);

    if (Mesprc==1)  fprintf (stdout,"\n the number of nodes   %ld",nn);
    nodes = new nodec [nn];
    Gtu->alloc_nodes (nn);
    
    for (i=0;i<nn;i++){
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("node number %ld is less than 1", __FILE__, __LINE__, __func__, j);
      if (j>nn){
	print_err("node number %ld is greater than the total number of nodes  %ld", __FILE__, __LINE__, __func__, j, nn);
      }
      Gtu->gnodes[j-1].read (in);
      nodes[j-1].read (in);
    }
    
    
    //  the number of nodes with prescribed values
    xfscanf (in,"%k%ld","number_of_constraints",&nnd);
    if (Mesprc==1)  fprintf (stdout,"\n number of constrained nodes  %ld",nnd);
    
    for (i=0;i<nnd;i++){
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("number of constrained node in function read is less than 1",__FILE__,__LINE__,__func__);
      if (j>nnd){
	print_err("number of constrained node in function read is greater than number of all nodes",__FILE__,__LINE__,__func__);
      }
      
      //  Gtu->dofcontr indicates whether DOFs indicators are described by a number or by a function
      //  in the case of problem with changing geometry
      Gtu->gnodes[j-1].constr (Gtu->dofcontr,in);
    }
    

    
    //  the number of elements
    xfscanf (in,"%ld",&ne);
    if (Mesprc==1)  fprintf (stdout,"\n number of elements  %ld",ne);
    elements = new elementc [ne];
    Gtu->alloc_elements (ne);
    
    for (i=0;i<ne;i++){
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("element number %ld is less than 1", __FILE__, __LINE__, __func__, j);
      if (j>ne){
	print_err("element number %ld is greater than total number of elements %ld", __FILE__, __LINE__, __func__, j, ne);
      }
      elements[j-1].read (in,j-1);
    }
    
  }

}

long couptop::mesh_check(void)
{
  long i,j,err=0;
  double fabx,faby,fabz;

  //nodes check
  if(Mt->nn != Tt->nn){
    err = 1;
    fprintf (stdout,"\n\n Different number of nodes in mefel and trfel input file \n\n");
    return err;
  }

  for(i=0;i<Mt->nn;i++){
    fabx = Gtm->gnodes[i].x - Gtt->gnodes[i].x;  
    faby = Gtm->gnodes[i].y - Gtt->gnodes[i].y;  
    fabz = Gtm->gnodes[i].z - Gtt->gnodes[i].z;  
    
    if(fabs(fabx) >= 1.0e-4){
      err = 2;
      fprintf (stdout,"\n\n Different x coordinates in node number %ld\n\n",i+1);
      return err;
    }
    if(fabs(faby) >= 1.0e-4){
      err = 3;
      fprintf (stdout,"\n\n Different y coordinates in node number %ld\n\n",i+1);
      return err;
    }
    if(fabs(fabz) >= 1.0e-4){
      err = 4;
      fprintf (stdout,"\n\n Different z coordinates in node number %ld\n\n",i+1);
      return err;
    }
  } 
  
  //elems check
  for(i=0;i<Tt->ne;i++){
    for(j=0;j<Gtt->gelements[i].nne;j++){
      if((Gtt->gelements[i].nodes[j]) != (Gtm->gelements[i].nodes[j])){
	fprintf (stdout,"\n\n Different node conectivities in element number %ld, node number %ld",i+1,j+1);
	fprintf (stdout,"\n Trfel node number %ld, Mefel node number %ld\n\n",Gtt->gelements[i].nodes[j]+1,Gtm->gelements[i].nodes[j]+1);
	err = 5;
	return err;
      }
    }
  }
  
  return err;
}



/**
   function returns element type
   
   @param eid - element id
*/
elemtypec couptop::give_elem_type (long eid)
{
  return (elements[eid].te);
}

/**
   function returns number of DOFs of node
   
   @param nid - node id
*/
long couptop::give_ndofn (long nid)
{
  return Gtu->give_ndofn (nid);
}

/**
   function returns number of DOFs of element
   
   @param eid - element id
   
   JK, 22.2.2002
*/
long couptop::give_ndofe (long eid)
{
  long ndofe;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{         ndofe=Cbar->tndofe;         break;  }
  case coupquad:{        ndofe=Cquad->tndofe;        break;  }
  case coupaxiquad:{     ndofe=Caxiq->tndofe;        break;  }
  case axisymfc:{        ndofe=Caxifc->tndofe;       break;  }
  case couphex:{         ndofe=Chex->tndofe;         break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return ndofe;
}

long couptop::give_nne (long eid)
{
  long nne;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
    //case coupbar:{         ndofe=Cbar->tndofe;         break;  }
    //case coupquad:{        ndofe=Cquad->tndofe;        break;  }
    //case coupaxiquad:{     ndofe=Caxiq->tndofe;        break;  }
  case axisymfc:{        nne=Caxifc->mnne;       break;  }
    //case couphex:{         ndofe=Chex->tndofe;         break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nne;
}



/**
   function returns total number of upper integration points on one element
   
   @param eid - element id
   
   07/05/2010 TKr
*/
long couptop::give_upper_tnip (long eid)
{
  long nip;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{               nip=Cbar->tnipu;     break;  }
  case coupquad:{              nip=Cquad->tnipu;    break;  }
  case coupaxiquad:{           nip=Caxiq->tnipu;    break;  }
  case couphex:{               nip=Chex->tnipu;     break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nip;
}


/**
   function returns total number of lower integration points on one element
   
   @param eid - element id
   
   07/05/2010 TKr
*/
long couptop::give_lower_tnip (long eid)
{
  long nip;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{               nip=Cbar->tnipl;     break;  }
  case coupquad:{              nip=Cquad->tnipl;    break;  }
  case coupaxiquad:{           nip=Caxiq->tnipl;    break;  }
  case couphex:{               nip=Chex->tnipl;     break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nip;
}



/**
   function returns number of integration points
   
   @param eid - element id
   @param ri,ci - row and column indices
   
   30.12.2001
*/
long couptop::give_upper_nip (long eid,long ri,long ci)
{
  long nip;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{               nip=Cbar->nipu[ri][ci];     break;  }
  case coupquad:{              nip=Cquad->nipu[ri][ci];    break;  }
  case coupaxiquad:{           nip=Caxiq->nipu[ri][ci];    break;  }
  case couphex:{               nip=Chex->nipu[ri][ci];     break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nip;
}

/**
   function returns number of integration points
   
   @param eid - element id
   @param ri,ci - row and column indices
   
   30.12.2001
*/
long couptop::give_lower_nip (long eid,long ri,long ci)
{
  long nip;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{               nip=Cbar->nipl[ri][ci];     break;  }
  case coupquad:{              nip=Cquad->nipl[ri][ci];    break;  }
  case coupaxiquad:{           nip=Caxiq->nipl[ri][ci];    break;  }
  case couphex:{               nip=Chex->nipl[ri][ci];     break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nip;
}

/**
   function returns numbers of nodes defining apropriate element
   
   @param eid - element id
   @param nodes - array containing nodes
   
   22.7.2001
*/
void couptop::give_elemnodes (long eid,ivector &nodes)
{
  Gtu->give_nodes (eid,nodes);
}

/**
   function extracts all code numbers of actual element
   
   @param eid - element id
   @param cn - array containing all code numbers on element
   
   25.6.2001
*/
void couptop::give_code_numbers (long eid,long *cn)
{
  Gtu->give_code_numbers (eid,cn);
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,y - vectors containing node coordinates
   @param eid - element id
   
   19.7.2001
*/
void couptop::give_node_coord2d (vector &x,vector &y,long eid)
{
  Gtu->give_node_coord2d (x,y,eid);
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,z - vectors containing coordinates
   @param eid - element id
   
   19.7.2001
*/
void couptop::give_node_coord2dxz (vector &x,vector &z,long eid)
{
  Gtu->give_node_coord2dxz (x,z,eid);
}

/**
   function returns node coordinates of the appropriate element
   
   @param x,y,z - vectors containing coordinates
   @param eid - element id
   
   19.7.2001
*/
void couptop::give_node_coord3d (vector &x,vector &y,vector &z,long eid)
{
  Gtu->give_node_coord3d (x,y,z,eid);
}

/**
   function returns number of blocks in characteristic matrices
   
   @param eid - element id
   
   8.5.2002
*/
long couptop::give_mnb (long eid)
{
  long mnb;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){

  case coupbar:{      mnb=Cbar->mnb;    break;  }
  case coupquad:{     mnb=Cquad->mnb;   break;  }
  case coupaxiquad:{  mnb=Caxiq->mnb;   break;  }
  case couphex:{      mnb=Chex->mnb;    break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return mnb;
}

/**
   function returns number nodes on elements
   
   @param te - element type
   
   27.2.2002
*/
/*
long couptop::give_nne_inner (elemtypec te)
{
  long nne;
  switch (te){
  case coupbar:{      nne=Cbar->nne;     break;  }
  case coupquad:{    nne=Cquad->nne;   break;  }
  case coupaxiquad:{  nne=Caxiq->nne;    break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return nne;
}
*/

/**
   function returns number of strain/stress components
   
   @param eid - element id
   
   30.12.2001
*/
long couptop::give_ncompstr (long eid)
{
  long ncompstr;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{       ncompstr=Cbar->tnmcomp;      break;  }
  case coupquad:{      ncompstr=Cquad->tnmcomp;     break;  }
  case coupaxiquad:{   ncompstr=Caxiq->tnmcomp;     break;  }
  case axisymfc:{      ncompstr=Caxifc->ncompstr;   break;  }
  case couphex:{       ncompstr=Chex->tnmcomp;      break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return ncompstr;
}

/**
   function returns number of gradient/flux components
   
   @param eid - element id
   
   JK, 10.4.2019
*/
long couptop::give_ncompgrad (long eid)
{
  long ncompgrad;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
    //case coupbar:{       ncompstr=Cbar->tnmcomp;      break;  }
    //case coupquad:{      ncompstr=Cquad->tnmcomp;     break;  }
    //case coupaxiquad:{   ncompstr=Caxiq->tnmcomp;     break;  }
  case axisymfc:{      ncompgrad=Caxifc->ncompgrad;     break;  }
    //case couphex:{       ncompstr=Chex->tnmcomp;      break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return ncompgrad;
}

/**
   function returns number of gradient/flux components
   
   @param eid - element id
   
   JK, 10.4.2019
*/
long couptop::give_ncompdispl (long eid)
{
  long ncompdispl;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
    //case coupbar:{       ncompstr=Cbar->tnmcomp;      break;  }
    //case coupquad:{      ncompstr=Cquad->tnmcomp;     break;  }
    //case coupaxiquad:{   ncompstr=Caxiq->tnmcomp;     break;  }
  case axisymfc:{      ncompdispl=Caxifc->ncompdispl;     break;  }
    //case couphex:{       ncompstr=Chex->tnmcomp;      break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return ncompdispl;
}

/**
   function returns the number of transported media
   
   @param eid - element id
   
   JK, 10.4.2019
*/
long couptop::give_ntm (long eid)
{
  long ntm;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
    //case coupbar:{       ncompstr=Cbar->tnmcomp;      break;  }
    //case coupquad:{      ncompstr=Cquad->tnmcomp;     break;  }
    //case coupaxiquad:{   ncompstr=Caxiq->tnmcomp;     break;  }
  case axisymfc:{      ntm=Caxifc->ntm;     break;  }
    //case couphex:{       ncompstr=Chex->tnmcomp;      break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  return ntm;
}

/**
   eid - element id
   bi - block id
   
   block id is needed for shell elements, first block is plane stress,
   second block is a plate
   
*/
strastrestate couptop::give_ssst (long eid,long /*bi*/)
{
  strastrestate ssst;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case coupbar:{       ssst = Cbar->ssst;          break;  }
  case coupquad:{      ssst = Cquad->ssst;         break;  }
  case coupaxiquad:{   ssst = Caxiq->ssst;         break;  }
  case axisymfc:{      ssst = Caxifc->ssst;        break;  }
  case couphex:{       ssst = Chex->ssst;          break;  }
  default:{
    print_err("unknown element type %d is required for element %ld", __FILE__, __LINE__, __func__, int(te), eid+1);
  }
  }
  
  return ssst;
}


/**
  Function returns number of integration points.
    
  @param eid - element id
  
  @return The function returns number of integration points of the given element block.
  
  31.3.2019, JK
*/
long couptop::give_nip (long eid)
{
  long nip=0;
  elemtypec te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case axisymfc:{       nip=Caxifc->nip;      break;  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return nip;
}

/**
   function determinates the number of boundary objects and number of nodes on them
   boundary objects are end nodes in 1D, edges in 2D and surfaces in 3D
   
   @param eid - element id
   @param bco - number of boundary objects
   @param ncompbco - number of components on one boundary object
   
   JK, 3.5.2019
*/
void couptop::give_nbobjects (long eid,long &bco,long &ncompbco)
{
  
  switch (elements[eid].te){
  case axisymfc:{
    bco=Caxifc->ned;
    ncompbco=Caxifc->tnned;
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}
