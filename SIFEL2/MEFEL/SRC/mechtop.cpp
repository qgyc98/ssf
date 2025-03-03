#include "mechtop.h"
#include "global.h"
#include "mechmat.h"
#include "probdesc.h"
#include "mechbclc.h"
#include "elemhead.h"
#include "elemswitch.h"
#include "gtopology.h"
#include "difcalc.h"
#include "mathem.h"
#include "element.h"
#include "node.h"
#include "endnodem.h"
#include "edgem.h"
#include "siftop_element_types.h"
#include "loadcase.h"
#include "intpoints.h"
#include "globmat.h"
#include "siftop.h"
#include "saddpoint.h"
#include "vecttens.h"
#include "ipmap.h"
#include "gnodvalvm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <limits>
#include <math.h>



/**
  Constructor initializes data members to zero values.

  Created by JK,
*/
mechtop::mechtop ()
{
  //  number of nodes
  nn=0;
  //  number of constrained nodes
  ncn=0;
  //  number of elements
  ne=0;
  //  number of layered nodes
  nln=0;

  // maximum number of dofs at node is set to uninitialized
  max_ndofn = -1;

  // maximum number of stress/strain components at elements
  max_ncompstr = 0;

  //  node array
  nodes=NULL;
  //  element array
  elements=NULL;

  tnip = 0;
  dist=NULL;  
  nadjip=NULL;
  adjip=NULL;
  nodedispl=NULL;
  nodeforce=NULL;
  domvol=0.0;
  node_radius=NULL;
}



/**
  Destructor deallocates used memory.

  Created by JK
  Modified by TKo
*/
mechtop::~mechtop ()
{
  long i;
  
  delete [] nodes;  
  delete [] elements;
  delete [] nadjip;
  if (adjip){
    for (i=0; i<tnip; i++)  
      delete [] adjip[i];  
  }
  if (dist)
  {
    for (i=0; i<tnip; i++)  
      delete [] dist[i];
  }
  delete [] adjip;  
  delete [] dist;
  
  if (nodedispl != NULL){
    for (i=0;i<nn;i++){
      delete [] nodedispl[i];
    }
    delete [] nodedispl;
  }
  if (nodeforce != NULL){
    for (i=0;i<nn;i++){
      delete [] nodeforce[i];
    }
    delete [] nodeforce;
  }
  delete [] node_radius;
}



/**
  Function reads basic topology informations.

  @param in - input stream

  @return The function does not return anything.

  Created by JK,
*/
void mechtop::read (XFILE *in)
{
  long i,j,k,l,ndofn;
  
  // ********************
  //  node data reading
  // ********************
  xfscanf (in,"%k%ld","number_of_nodes",&nn);
  Gtm->alloc_nodes (nn);

  if (Mespr==1)  fprintf (stdout,"\n number of nodes  %ld",nn);
  nodes = new node [nn];

  if (Mp->tprob==forced_dynamics)
    node_radius = new double [nn];

  switch (Mp->tprob){
  case linear_statics:
  case eigen_dynamics:
  case forced_dynamics:
  case linear_stability:
  case mat_nonlinear_statics:
  case geom_nonlinear_statics:
  case mech_timedependent_prob:
  case growing_mech_structure:
  case load_balancing:
  case lin_floating_subdomain:
  case hemivar_inequal:
  case earth_pressure:{
    for (i=0;i<nn;i++){
      j=nn+1;
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("node number in function read is less than 1",__FILE__,__LINE__,__func__);
      if (j>nn)
	print_err("node number in function read is greater than number of nodes",__FILE__,__LINE__,__func__);
      Gtm->gnodes[j-1].read (in);
      nodes[j-1].read (in);
      if (Mp->tprob==forced_dynamics){
	xfscanf (in,"%lf",node_radius+j-1);
      }
    }
    break;
  }
  case layered_linear_statics:{
    xfscanf (in,"%k%ld","number_of_layered_nodes",&nln);
    Gtm->alloc_lnodes (nln);
    for (i=0;i<nln;i++){
      xfscanf (in,"%ld",&l);
      xfscanf (in,"%ld",&(Gtm->lgnodes[l-1].nl));
      Gtm->lgnodes[l-1].nodes = new long [Gtm->lgnodes[l-1].nl];
      for (k=0;k<Gtm->lgnodes[l-1].nl;k++){
	j=nn+1;
	xfscanf (in,"%ld",&j);
	if (j<1)
	  print_err("node number in function read is less than 1",__FILE__,__LINE__,__func__);
	if (j>nn)
	  print_err("node number in function read is greater than number of nodes",__FILE__,__LINE__,__func__);
	Gtm->gnodes[j-1].read (in);
	nodes[j-1].read (in);
	Gtm->lgnodes[l-1].nodes[k]=j-1;
      }
    }
    Gtm->unodelnode ();
    break;
  }

  case nonlin_floating_subdomain:{
    /*
    //  number of subdomains
    xfscanf (in,"%k%ld","number_of_floating_subdomains",&(Gtm->flsub.nsd));
    //  allocation of arrays
    Gtm->flsub.nnsd = new long [Gtm->flsub.nsd];
    Gtm->flsub.fnnsd = new long [Gtm->flsub.nsd+1];
    
    Gtm->flsub.fnnsd[0]=0;
    for (i=0;i<Gtm->flsub.nsd;i++){
      xfscanf (in,"%ld",&(Gtm->flsub.nnsd[i]));
      Gtm->flsub.fnnsd[i+1] = Gtm->flsub.fnnsd[i]+Gtm->flsub.nnsd[i];
    }
    Gtm->flsub.fnnsd[Gtm->flsub.nsd]--;
    */

    //  nodes
    for (i=0;i<nn;i++){
      j=nn+1;
      xfscanf (in,"%ld",&j);
      if (j<1)
	print_err("node number in function read is less than 1",__FILE__,__LINE__,__func__);
      if (j>nn)
	print_err("node number in function read is greater than number of nodes",__FILE__,__LINE__,__func__);
      Gtm->gnodes[j-1].read (in);
      nodes[j-1].read (in);
    }
    
    /*
    //  layered nodes
    xfscanf (in,"%ld",&nln);
    Gtm->alloc_lnodes (nln);
    for (i=0;i<nln;i++){
      xfscanf (in,"%ld",&l);
      xfscanf (in,"%ld",&(Gtm->lgnodes[l-1].nl));
      Gtm->lgnodes[l-1].nodes = new long [Gtm->lgnodes[l-1].nl];
      for (k=0;k<Gtm->lgnodes[l-1].nl;k++){
	j=nn+1;
	xfscanf (in,"%ld",&j);
	if (j<1)
	  print_err("node number in function read is less than 1",__FILE__,__LINE__,__func__);
	if (j>nn)
	  print_err("node number in function read is greater than number of nodes",__FILE__,__LINE__,__func__);
	Gtm->lgnodes[l-1].nodes[k]=j-1;
      }
    }
    //Gtm->unodelnode ();
    */
    break;
  }
    
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
  }
  }
  
  //  constrained nodes
  xfscanf (in,"%k%ld","number_of_constraints",&ncn);
  if (Mespr==1)  fprintf (stdout,"\n number of constrained nodes  %ld",ncn);

  if (Gtm->dofcontr) // in case of growing structures
  {
    // array of attained nodal values in the case of growing_mech_structure problems
    // addition of supports to the nodes that were free originally
    nodedispl = new double* [nn];
    memset(nodedispl, 0, sizeof(*nodedispl)*nn);
  }

  for (i=0;i<ncn;i++){
    k=nn+1;
    xfscanf (in,"%ld",&k);
    if (k<1)
      print_err("number of constrained node in function read is less than 1",__FILE__,__LINE__,__func__);
    if (k>nn){
      print_err("number of constrained node in function read is greater than number of all nodes",__FILE__,__LINE__,__func__);
    }
    Gtm->gnodes[k-1].constr(Gtm->dofcontr,in);

    if (Gtm->dofcontr) // in case of growing structures
    {      
      ndofn = give_ndofn(k-1);
      nodedispl[k-1] = new double[ndofn];
      memset(nodedispl[k-1], 0, sizeof(*nodedispl[k-1])*ndofn);
    }
  }

  // ***********************
  // ***********************
  //  element data reading
  // ***********************
  // ***********************
  xfscanf (in,"%k%ld","number_of_elements",&ne);
  Gtm->alloc_elements (ne);
  if (Mespr==1)  fprintf (stdout,"\n number of elements  %ld",ne);
  elements = new element [ne];
  
  switch (Mp->tprob){
  case linear_statics:
  case eigen_dynamics:
  case forced_dynamics:
  case linear_stability:
  case mat_nonlinear_statics:
  case geom_nonlinear_statics:
  case mech_timedependent_prob:
  case growing_mech_structure:
  case earth_pressure:
  case lin_floating_subdomain:
  case hemivar_inequal:
  case load_balancing:
  case layered_linear_statics:{
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
    break;
  }

  case nonlin_floating_subdomain:{
    /*
    Gtm->flsub.nesd = new long [Gtm->flsub.nsd];
    Gtm->flsub.fensd = new long [Gtm->flsub.nsd+1];
    Gtm->flsub.fensd[0]=0;
    //  number of elements on particular subdomains
    for (i=0;i<Gtm->flsub.nsd;i++){
      xfscanf (in,"%ld",&(Gtm->flsub.nesd[i]));
      Gtm->flsub.fensd[i+1]=Gtm->flsub.fensd[i]+Gtm->flsub.nesd[i];
    }
    
    if (ne!=Gtm->flsub.fensd[Gtm->flsub.nsd]){
      print_err("the sum of numbers of elements on subdomains\n is not equal to the total number of elements", __FILE__, __LINE__, __func__);
    }
    */
    
    for (i=0;i<ne;i++){
      j=ne+1;
      xfscanf (in,"%ld",&j);
      if ((j<1) || (j>ne))
	print_err ("element number is out of range <1,%ld>.",__FILE__,__LINE__, __func__, ne);
      elements[j-1].read (in,j-1);
    }
    break;
  }
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
  }
  }
  
  /*
  if (Mp->tprob==layered_linear_statics){
    fscanf (in,"%ld",&nle);
    Gtm->alloc_lelements (nle);
    for (i=0;i<nle;i++){
      fscanf (in,"%ld",&l);
      fscanf (in,"%ld",&(Gtm->lgelements[l-1].nl));
      Gtm->lgelements[l-1].elem = new long [Gtm->lgelements[l-1].nl];
      for (k=0;k<Gtm->lgelements[l-1].nl;k++){
	j=ne+1;
	fscanf (in,"%ld",&j);
	if (j<1)
	  fprintf (stderr,"\n\n element number in function mechtop::read is less than 1.\n");
	if (j>ne){
	  fprintf (stderr,"\n\n element number in function mechtop::read is");
	  fprintf (stderr,"\n greater than total number of elements.\n");
	}
	elements[j-1].read (in,j-1);
	Gtm->lgelements[l-1].elem[k]=j-1;
      }
      Gtm->lgelements[l-1].elemnmult ();
    }
  }
  else{
  */

  
  
  //  search for the hanging nodes
  //  variable gelements[i].nne is modified
  //  array gelements[i].nodes is modified
  if (Mp->homog < 3){
    Gtm->searching_hanging_nodes ();
    searching_hanging_nodes ();
    if (Mp->homog > 0) // hanging nodes generated from
      Gtm->hang_nodes_check ();
  }

  // function computes the maximum number of stress/strain components used
  // in the first order homogenization approaches
  if (Mp->homog > 2)
    Mt->set_maxncompstr();
  
  /*
  if (Mp->adaptivity > 0){
    //  initialization of auxinf variable
    for (i=0;i<ne;i++){
      nne = give_nne (i);
      order = give_degree (i);
      dim = give_dimension (i);
      
      Gtm->gelements[i].auxinf = nne*100+order*10+dim;
    }
  }
  */
  
  //  allocation and initiation of arrays lnso and leso
  //  lnso - list of nodes switched on
  //  leso - list of elements switched on
  //  3.11.2006
  Gtm->lneso_init ();

  if (Mp->tprob == growing_mech_structure){
    Gtm->auxinf_init ();
  }
  
  if (Gtm->rst==1){
    //  sequential topology is read
    //  it is used in the case of sequential version of
    //  any domain decomposition method
    Gtm->read_seq_top (in);
  }
  if (Gtm->rst==2){
    //  sequential topology is read in the form
    //  of Boolean matrices
    Gtm->read_seq_top (in);
    Mp->ssle->feti.read_booldata (in);
  }
  
  /*
  circumscribed_balls ();    
  inter ();
  */
  
}

/**
   function searches for hanging nodes
   
   30.3.2012
*/
void mechtop::searching_hanging_nodes (void)
{
  long i,j,k,l,nmn,npn,nneo,nnen,ndofeo,ndofen,ndofn,ri,ci,cumulci;
  ivector oldnodes,newnodes;
  vector weights;
  matrix tmat;
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    if (give_ndofe(i)!=Gtm->give_ndofe(i)){
      //  element was modified because of the hanging nodes
      
      //  the original number of nodes on element
      nneo = give_nne (i);
      //  the new number of nodes on element
      nnen = Gtm->give_nmne (i);
      
      reallocv (RSTCKIVEC(nneo,oldnodes));
      reallocv (RSTCKIVEC(nnen,newnodes));
      
      //  the original nodes on element
      Gtm->give_original_nodes (i,oldnodes);
      //  the new nodes on element
      Gtm->give_master_nodes (i,newnodes);
      
      //  the original number of DOFs on element
      ndofeo = give_ndofe (i);
      //  the new number of DOFs on element
      ndofen = Gtm->give_ndofe (i);
      
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
	  nmn=0-Gtm->gnodes[oldnodes[j]].ndofn;
	  
	  reallocv (RSTCKVEC(nmn,weights));
	  Gtm->approx_weights (Gtm->gnodes[oldnodes[j]].masentity,oldnodes[j],weights);

	  //  the number of DOFs in node
	  ndofn = give_ndofn (newnodes[npn]);

	  //  loop over the number of DOFs in actual node
	  for (k=0;k<ndofn;k++){
	    ci=cumulci+k;//ri;
	    
	    //  loop over the number of master nodes
	    for (l=0;l<nmn;l++){
	      //tmat[ri][ci]=1.0;
	      tmat[ri][ci]=weights[l];
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


/**
  Function investigates presence of prescribed displacements on elements.
  
  @return The function does not return anything.
 
  Created by JK, 25.7.2001
  Updated by TKo, 5.2013
*/
void mechtop::elemprescdisp (void)
{
  long i,j,k,ndofe;
  long nne, ndofn;
  ivector enod, cn;
  
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      ndofe=Gtm->give_ndofe (i);
      reallocv(RSTCKIVEC(ndofe, cn));
      give_code_numbers (i,cn.a);    
      for (j=0;j<ndofe;j++){
        if (cn[j]<0){
          elements[i].prescdispl=1;
          break;
        }
      }
    }
    if ((Mp->tprob == growing_mech_structure) && (elements[i].prescdispl==0))
    {
      // stored intial displacements at nodes can be
      // assumed as prescribed displacements in the case of growing structures
      nne = give_nne(i);
      reallocv(RSTCKIVEC(nne, enod));
      give_elemnodes(i, enod);
      for(j=0; j<nne; j++)
      {
        if (nodedispl[enod[j]]) // node has got initial displacements
        {
          ndofn = give_ndofn(enod[j]);
          for(k=0; k<ndofn; k++)
          {
            if (nodedispl[enod[j]][k] != 0.0) 
            {
              // initial displacement is nonzero => mark the element
              elements[i].prescdispl=1;
              break;
            }
          }
        }
        if(elements[i].prescdispl) 
        {
          // at least one node has got initial displacements =>
          // rest of nodes does not need to be investigated
          break;
        }
      }
    }
  }
}



/**
  Function investigates presence of temperature changes on elements.
  
  @param lcid - load case id

  @return The function does not return anything.
  
  Created by JK, 30.11.2002
*/
void mechtop::elempresctemp (long /*lcid*/)
{
  long i;
  ivector nod;

  if (Mp->temperature > 0){
    for (i=0;i<ne;i++){
      if (Gtm->leso[i]==1){
	elements[i].presctemp=1;
      }
    }
  }
  
  /*
  if (Mb->lc[lcid].tempchang==1){
    for (i=0;i<ne;i++){
      nne=give_nne (i);
      reallocv (RSTCKIVEC(nne,nod));
      dtelem = new double [nne];
      give_elemnodes (i,nod);
      Mb->lc[lcid].tempchanges (dtelem,nod);
      for (j=0;j<nne;j++){
	if (dtelem[j]!=0.0){
	  elements[i].presctemp=1;
	  break;
	}
      }
      delete [] dtelem;
      destrv (nod);
    }
  }
*/
}



/**
  Function returns element type.
  
  @param eid - element id

  @return The function returns type of the required element.

  Created by JK,
*/
elemtype mechtop::give_elem_type (long eid)
{
  return (elements[eid].te);
}



/**
  Function returns the number of DOFs of element.

  @param eid - element id

  @return The function returns the number of DOFs of the required element.

  Created by JK,
*/
long mechtop::give_ndofe (long eid)
{
  long ndofe=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{               ndofe=Bar2d->ndofe;              break; }
    case bar3d:{               ndofe=Bar3d->ndofe;              break; }
    case barq2d:{              ndofe=Barq2d->ndofe;             break; }
    case barq3d:{              ndofe=Barq3d->ndofe;             break; }
    case beam2d:{              ndofe=Beam2d->ndofe;             break; }
    case beam3d:{              ndofe=Beam3d->ndofe;             break; }
    case beamg3d:{             ndofe=Beam3dg->ndofe;            break; }
    case subsoilbeam:{         ndofe=Sbeam->ndofe;              break; }
    case beam2dsp:{            ndofe=Spbeam2d->ndofe;           break; }
    case spring_1:{            ndofe=Spring->give_ndofe(eid);   break; }
    case spring_2:{            ndofe=Spring->give_ndofe(eid);   break; }
    case spring_3:{            ndofe=Spring->give_ndofe(eid);   break; }
    case spring_4:{            ndofe=Spring->give_ndofe(eid);   break; }
    case spring_5:{            ndofe=Spring->give_ndofe(eid);   break; }
    case spring_6:{            ndofe=Spring->give_ndofe(eid);   break; }
    case planeelementlt:{      ndofe=Pelt->ndofe;               break; }
    case planeelementqt:{      ndofe=Peqt->ndofe;               break; }
    case planeelementrotlt:{   ndofe=Perlt->ndofe;              break; }
    case planeelementlq:{      ndofe=Pelq->ndofe;               break; }
    case planeelementqq:{      ndofe=Peqq->ndofe;               break; }
    case planeelementrotlq:{   ndofe=Perlq->ndofe;              break; }
    case planeelementsubqt:{   ndofe=Pesqt->ndofe;              break; }
    case planequadinterface:{  ndofe=Pqifc->ndofe;              break; }
    case cctel:{               ndofe=Cct->ndofe;                break; }
    case dktel:{               ndofe=Dkt->ndofe;                break; }
    case dstel:{               ndofe=Dst->ndofe;                break; }
    case q4plateel:{           ndofe=Q4pl->ndofe;               break; }
      //case argyristr:{           ndofe=Argtr->ndofe;              break; }
    case argyristr:{           ndofe=Argtrpl->ndofe;            break; }
    case quadkirch:{           ndofe=Qkirch->ndofe;             break; }
    case dkqel:{               ndofe=Dkqelem->ndofe;                break; }
    case subsoilplatetr:{      ndofe=Spltr->ndofe;              break; }
    case subsoilplateq:{       ndofe=Splq->ndofe;               break; }
    case axisymmlt:{           ndofe=Asymlt->ndofe;             break; }
    case axisymmqt:{           ndofe=Asymqt->ndofe;             break; }
    case axisymmlq:{           ndofe=Asymlq->ndofe;             break; }
    case axisymmqq:{           ndofe=Asymqq->ndofe;             break; }
    case axisymmcq:{           ndofe=Asymcq->ndofe;             break; }
    case axisymmlqintface:{    ndofe=Asymlqifc->ndofe;          break; }
    case shelltrelem:{         ndofe=Shtr->ndofe;               break; }
    case shelltrmelem:{        ndofe=Shtrm->ndofe;              break; }
    case shellqelem:{          ndofe=Shq->ndofe;                break; }
    case lineartet:{           ndofe=Ltet->ndofe;               break; }
    case quadrtet:{            ndofe=Qtet->ndofe;               break; }
    case linearhex:{           ndofe=Lhex->ndofe;               break; }
    case quadrhex:{            ndofe=Qhex->ndofe;               break; }
    case lineartetrot:{        ndofe=Ltetrot->ndofe;            break; }
    case linearhexrot:{        ndofe=Lhexrot->ndofe;            break; }
    case linearwed:{           ndofe=Lwed->ndofe;               break; }
    case quadrwed:{            ndofe=Qwed->ndofe;               break; }
    case hexintface:{          ndofe=Hexifc->ndofe;             break; }
    case particleelem:{        ndofe=Pelem->ndofe;              break; }
    case tetralatt:{           ndofe=Tlatt->ndofe;              break; }
    default:{
      print_err("unknown element type '%d' is required on element %ld",__FILE__,__LINE__,__func__, te, eid+1);
    }
  }
  return ndofe;
}



/**
  Function returns number of DOFs of node.

  @param nid - node id

  @return The function returns number of DOFs of the required node.

  Created by JK,
*/
long mechtop::give_ndofn (long nid)
{
  //return Gtm->give_ndofn (nid);
  
  long ndofn;
  
  ndofn = Gtm->give_ndofn (nid);
  if (ndofn<0){
    //  the node is a hanging node
    ndofn = Gtm->give_ndofn (Gtm->gnodes[nid].mnodes[0]);
  }
  
  return ndofn;
}



/**
  Function returns appropriate DOF of node.

  @param nid - node id
  @param n - number of required DOF

  @return The function returns DOF of the required node.

  Created by JK,
*/
long mechtop::give_dof (long nid,long n)
{
  return Gtm->give_dof (nid,n);
}



/**
  The function returns nodal coordinates for the nid-th node.

  @param nid - node id
  @param c   - %vector of coordinate (output)

  @return The function returns nodal coordinates in parameter c.

  Created by Tomas Koudelka, 11.2011
*/
void mechtop::give_nodal_coord (long nid, vector &c)
{
  c[0] = Gtm->gnodes[nid].x;
  c[1] = Gtm->gnodes[nid].y;
  c[2] = Gtm->gnodes[nid].z;
}




/**
  Function returns number of nodes of the element.
  
  @param eid - element id

  @return The function returns number of nodes of the required element.

  Created by JK,
*/
long mechtop::give_nne (long eid)
{
  long nne;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              nne=Bar2d->nne;     break; }
    case bar3d:{              nne=Bar3d->nne;     break; }
    case barq2d:{             nne=Barq2d->nne;    break; }
    case barq3d:{             nne=Barq3d->nne;    break; }
    case beam2d:{             nne=Beam2d->nne;    break; }
    case beam3d:{             nne=Beam3d->nne;    break; }
    case beamg3d:{            nne=Beam3dg->nne;   break; }
    case subsoilbeam:{        nne=Sbeam->nne;     break; }
    case beam2dsp:{           nne=Spbeam2d->nne;  break; }
    case spring_1:{           nne=Spring->nne;    break; }
    case spring_2:{           nne=Spring->nne;    break; }
    case spring_3:{           nne=Spring->nne;    break; }
    case spring_4:{           nne=Spring->nne;    break; }
    case spring_5:{           nne=Spring->nne;    break; }
    case spring_6:{           nne=Spring->nne;    break; }
    case planeelementlt:{     nne=Pelt->nne;      break; }
    case planeelementqt:{     nne=Peqt->nne;      break; }
    case planeelementrotlt:{  nne=Perlt->nne;     break; }
    case planeelementlq:{     nne=Pelq->nne;      break; }
    case planeelementqq:{     nne=Peqq->nne;      break; }
    case planeelementrotlq:{  nne=Perlq->nne;     break; }
    case planeelementsubqt:{  nne=Pesqt->nne;     break; }
    case planequadinterface:{ nne=Pqifc->nne;     break; }
    case cctel:{              nne=Cct->nne;       break; }
    case dktel:{              nne=Dkt->nne;       break; }
    case dstel:{              nne=Dst->nne;       break; }
    case q4plateel:{          nne=Q4pl->nne;      break; }
      //case argyristr:{          nne=Argtr->nne;     break; }
    case argyristr:{          nne=Argtrpl->nne;   break; }
    case quadkirch:{          nne=Qkirch->nne;    break; }
    case dkqel:{              nne=Dkqelem->nne;       break; }
    case subsoilplatetr:{     nne=Spltr->nne;     break; }
    case subsoilplateq:{      nne=Splq->nne;      break; }
    case axisymmlt:{          nne=Asymlt->nne;    break; }
    case axisymmqt:{          nne=Asymqt->nne;    break; }
    case axisymmlq:{          nne=Asymlq->nne;    break; }
    case axisymmqq:{          nne=Asymqq->nne;    break; }
    case axisymmcq:{          nne=Asymcq->nne;    break; }
    case axisymmlqintface:{   nne=Asymlqifc->nne; break; }
    case shelltrelem:{        nne=Shtr->nne;      break; }
    case shelltrmelem:{       nne=Shtrm->nne;     break; }
    case shellqelem:{         nne=Shq->nne;       break; }
    case lineartet:{          nne=Ltet->nne;      break; }
    case quadrtet:{           nne=Qtet->nne;      break; }
    case linearhex:{          nne=Lhex->nne;      break; }
    case quadrhex:{           nne=Qhex->nne;      break; }
    case lineartetrot:{       nne=Ltetrot->nne;   break; }
    case linearhexrot:{       nne=Lhexrot->nne;   break; }
    case linearwed:{          nne=Lwed->nne;      break; }
    case quadrwed:{           nne=Qwed->nne;      break; }
    case hexintface:{         nne=Hexifc->nne;    break; }
    case particleelem:{       nne=Pelem->nne;     break; }
    case tetralatt:{          nne=Tlatt->nne;     break; }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__, __func__);
    }
  }
  return nne;
}



/**
  Function returns number of integration points.
    
  @param eid - element id
  @param ri,ci - row and column indices
  
  @return The function returns number of integration points of the given element block.

  Created by JK, 30.12.2001
*/
long mechtop::give_nip (long eid,long ri,long ci)
{
  long nip=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              nip=Bar2d->nip[ri][ci];      break;  }
    case bar3d:{              nip=Bar3d->nip[ri][ci];      break;  }
    case barq2d:{             nip=Barq2d->nip[ri][ci];     break;  }
    case barq3d:{             nip=Barq3d->nip[ri][ci];     break;  }
    case beam2d:{             nip=Beam2d->nip[ri][ci];     break;  }
    case beam3d:{             nip=Beam3d->nip[ri][ci];     break;  }
    case beamg3d:{            nip=Beam3dg->nip[ri][ci];    break;  }
    case subsoilbeam:{        nip=Sbeam->nip[ri][ci];      break;  }
    case beam2dsp:{           nip=Spbeam2d->nip[ri][ci];   break;  }
    case spring_1:{           nip=Spring->nip[ri][ci];     break;  }
    case spring_2:{           nip=Spring->nip[ri][ci];     break;  }
    case spring_3:{           nip=Spring->nip[ri][ci];     break;  }
    case spring_4:{           nip=Spring->nip[ri][ci];     break;  }
    case spring_5:{           nip=Spring->nip[ri][ci];     break;  }
    case spring_6:{           nip=Spring->nip[ri][ci];     break;  }
    case planeelementlt:{     nip=Pelt->nip[ri][ci];       break;  }
    case planeelementqt:{     nip=Peqt->nip[ri][ci];       break;  }
    case planeelementrotlt:{  nip=Perlt->nip[ri][ci];      break;  }
    case planeelementlq:{     nip=Pelq->nip[ri][ci];       break;  }
    case planeelementqq:{     nip=Peqq->nip[ri][ci];       break;  }
    case planeelementrotlq:{  nip=Perlq->nip[ri][ci];      break;  }
    case planeelementsubqt:{  nip=Pesqt->nip[ri][ci];      break;  }
    case planequadinterface:{ nip=Pqifc->nip[ri][ci];      break;  }
    case cctel:{              nip=Cct->nip[ri][ci];        break;  }
    case dktel:{              nip=Dkt->nip[ri][ci];        break;  }
    case dstel:{              nip=Dst->nip[ri][ci];        break;  }
    case q4plateel:{          nip=Q4pl->nip[ri][ci];       break;  }
      //case argyristr:{          nip=Argtr->nip[ri][ci];      break;  }
    case argyristr:{          nip=Argtrpl->nip[ri][ci];    break;  }
    case quadkirch:{          nip=Qkirch->nip[ri][ci];     break;  }
    case dkqel:{              nip=Dkqelem->nip[ri][ci];        break;  }
    case subsoilplatetr:{     nip=Spltr->nip[ri][ci];      break;  }
    case subsoilplateq:{      nip=Splq->nip[ri][ci];       break;  }
    case axisymmlt:{          nip=Asymlt->nip[ri][ci];     break;  }
    case axisymmqt:{          nip=Asymqt->nip[ri][ci];     break;  }
    case axisymmlq:{          nip=Asymlq->nip[ri][ci];     break;  }
    case axisymmqq:{          nip=Asymqq->nip[ri][ci];     break;  }
    case axisymmcq:{          nip=Asymcq->nip[ri][ci];     break;  }
    case axisymmlqintface:{   nip=Asymlqifc->nip[ri][ci];  break;  }
    case shelltrelem:{        nip=Shtr->nip[ri][ci];       break;  }
    case shelltrmelem:{       nip=Shtrm->nip[ri][ci];      break;  }
    case shellqelem:{         nip=Shq->nip[ri][ci];        break;  }
    case lineartet:{          nip=Ltet->nip[ri][ci];       break;  }
    case quadrtet:{           nip=Qtet->nip[ri][ci];       break;  }
    case linearhex:{          nip=Lhex->nip[ri][ci];       break;  }
    case quadrhex:{           nip=Qhex->nip[ri][ci];       break;  }
    case lineartetrot:{       nip=Ltetrot->nip[ri][ci];    break;  }
    case linearhexrot:{       nip=Lhexrot->nip[ri][ci];    break;  }
    case linearwed:{          nip=Lwed->nip[ri][ci];       break;  }
    case quadrwed:{           nip=Qwed->nip[ri][ci];       break;  }
    case hexintface:{         nip=Hexifc->nip[ri][ci];     break; }
    case particleelem:{       nip=Pelem->nip[ri][ci];      break;  }
    case tetralatt:{          nip=Tlatt->nip[ri][ci];      break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return nip;
}



/**
  Function returns integration order of stiffness %matrix for the given element id and block ids.
    
  @param eid - element id
  @param ri,ci - row and column indices
  
  @return The function returns number of integration order of the given element block.

  Created by TKo, 11.11.2013
*/
long mechtop::give_intordsm (long eid,long ri,long ci)
{
  long ret=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              ret=Bar2d->intordsm[ri][ci];      break;  }
    case bar3d:{              ret=Bar3d->intordsm[ri][ci];      break;  }
    case barq2d:{             ret=Barq2d->intordsm[ri][ci];     break;  }
    case barq3d:{             ret=Barq3d->intordsm[ri][ci];     break;  }
    case beam2d:{             ret=Beam2d->intordsm[ri][ci];     break;  }
    case beam3d:{             ret=Beam3d->intordsm[ri][ci];     break;  }
    case beamg3d:{            ret=Beam3dg->intordsm[ri][ci];    break;  }
    case subsoilbeam:{        ret=Sbeam->intordsm[ri][ci];      break;  }
    case spring_1:{           ret=Spring->intordsm[ri][ci];     break;  }
    case spring_2:{           ret=Spring->intordsm[ri][ci];     break;  }
    case spring_3:{           ret=Spring->intordsm[ri][ci];     break;  }
    case spring_4:{           ret=Spring->intordsm[ri][ci];     break;  }
    case spring_5:{           ret=Spring->intordsm[ri][ci];     break;  }
    case spring_6:{           ret=Spring->intordsm[ri][ci];     break;  }
    case planeelementlt:{     ret=Pelt->intordsm[ri][ci];       break;  }
    case planeelementqt:{     ret=Peqt->intordsm[ri][ci];       break;  }
    case planeelementrotlt:{  ret=Perlt->intordsm[ri][ci];      break;  }
    case planeelementlq:{     ret=Pelq->intordsm[ri][ci];       break;  }
    case planeelementqq:{     ret=Peqq->intordsm[ri][ci];       break;  }
    case planeelementrotlq:{  ret=Perlq->intordsm[ri][ci];      break;  }
    case planeelementsubqt:{  ret=Pesqt->intordsm[ri][ci];      break;  }
    case planequadinterface:{ ret=Pqifc->intordsm[ri][ci];      break;  }
    case cctel:{              ret=Cct->intordsm[ri][ci];        break;  }
    case dktel:{              ret=Dkt->intordsm[ri][ci];        break;  }
    case dstel:{              ret=Dst->intordsm[ri][ci];        break;  }
    case q4plateel:{          ret=Q4pl->intordsm[ri][ci];       break;  }
      //case argyristr:{          ret=Argtr->intordsm[ri][ci];      break;  }
    case argyristr:{          ret=Argtrpl->intordsm[ri][ci];    break;  }
    case quadkirch:{          ret=Qkirch->intordsm[ri][ci];     break;  }
    case dkqel:{              ret=Dkqelem->intordsm[ri][ci];        break;  }
    case subsoilplatetr:{     ret=Spltr->intordsm[ri][ci];      break;  }
    case subsoilplateq:{      ret=Splq->intordsm[ri][ci];       break;  }
    case axisymmlt:{          ret=Asymlt->intordsm[ri][ci];     break;  }
    case axisymmqt:{          ret=Asymqt->intordsm[ri][ci];     break;  }
    case axisymmlq:{          ret=Asymlq->intordsm[ri][ci];     break;  }
    case axisymmqq:{          ret=Asymqq->intordsm[ri][ci];     break;  }
    case axisymmcq:{          ret=Asymcq->intordsm[ri][ci];     break;  }
    case axisymmlqintface:{   ret=Asymlqifc->intordsm[ri][ci];  break;  }
    case shelltrelem:{        ret=Shtr->intordsm[ri][ci];       break;  }
    case shelltrmelem:{       ret=Shtrm->intordsm[ri][ci];      break;  }
    case shellqelem:{         ret=Shq->intordsm[ri][ci];        break;  }
    case lineartet:{          ret=Ltet->intordsm[ri][ci];       break;  }
    case quadrtet:{           ret=Qtet->intordsm[ri][ci];       break;  }
    case linearhex:{          ret=Lhex->intordsm[ri][ci];       break;  }
    case quadrhex:{           ret=Qhex->intordsm[ri][ci];       break;  }
    case lineartetrot:{       ret=Ltetrot->intordsm[ri][ci];    break;  }
    case linearhexrot:{       ret=Lhexrot->intordsm[ri][ci];    break;  }
    case tetralatt:{          ret=0;                            break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return ret;
}



/**
  Function returns total number of integration points of the given element.
   
  @param eid - element id
   
  @return The function returns total number of integration points.

  Created by JK, 21.6.2004
*/
long mechtop::give_tnip (long eid)
{
  long tnip = -1;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:   {           tnip=Bar2d->tnip;      break;  }
    case bar3d:   {           tnip=Bar3d->tnip;      break;  }
    case barq2d:  {           tnip=Barq2d->tnip;     break;  }
    case barq3d:  {           tnip=Barq3d->tnip;     break;  }
    case beam2d:{             tnip=Beam2d->tnip;     break;  }
    case beam3d:{             tnip=Beam3d->tnip;     break;  }
    case subsoilbeam:{        tnip=Sbeam->tnip;      break;  }
    case beamg3d:{            tnip=Beam3dg->tnip;    break;  }
    case beam2dsp:{           tnip=Spbeam2d->tnip;   break;  }
    case spring_1:{           tnip=Spring->tnip;     break;  }
    case spring_2:{           tnip=Spring->tnip;     break;  }
    case spring_3:{           tnip=Spring->tnip;     break;  }
    case spring_4:{           tnip=Spring->tnip;     break;  }
    case spring_5:{           tnip=Spring->tnip;     break;  }
    case spring_6:{           tnip=Spring->tnip;     break;  }
    case planeelementlt:{     tnip=Pelt->tnip;       break;  }
    case planeelementqt:{     tnip=Peqt->tnip;       break;  }
    case planeelementrotlt:{  tnip=Perlt->tnip;      break;  }
    case planeelementlq:{     tnip=Pelq->tnip;       break;  }
    case planeelementqq:{     tnip=Peqq->tnip;       break;  }
    case planeelementrotlq:{  tnip=Perlq->tnip;      break;  }
    case planeelementsubqt:{  tnip=Pesqt->tnip;      break;  }
    case planequadinterface:{ tnip=Pqifc->tnip;      break;  }
    case cctel:{              tnip=Cct->tnip;        break;  }
    case dktel:{              tnip=Dkt->tnip;        break;  }
    case dstel:{              tnip=Dst->tnip;        break;  }
    case q4plateel:{          tnip=Q4pl->tnip;       break;  }
      //case argyristr:{          tnip=Argtr->tnip;      break;  }
    case argyristr:{          tnip=Argtrpl->tnip;    break;  }
    case quadkirch:{          tnip=Qkirch->tnip;     break;  }
    case dkqel:{              tnip=Dkqelem->tnip;        break;  }
    case subsoilplatetr:{     tnip=Spltr->tnip;      break;  }
    case subsoilplateq:{      tnip=Splq->tnip;       break;  }
    case axisymmlt:{          tnip=Asymlt->tnip;     break;  }
    case axisymmqt:{          tnip=Asymqt->tnip;     break;  }
    case axisymmlq:{          tnip=Asymlq->tnip;     break;  }
    case axisymmqq:{          tnip=Asymqq->tnip;     break;  }
    case axisymmcq:{          tnip=Asymcq->tnip;     break;  }
    case axisymmlqintface:{   tnip=Asymlqifc->tnip;  break;  }
    case shelltrelem:{        tnip=Shtr->tnip;       break;  }
    case shelltrmelem:{       tnip=Shtrm->tnip;      break;  }
    case shellqelem:{         tnip=Shq->tnip;        break;  }
    case lineartet:{          tnip=Ltet->tnip;       break;  }
    case quadrtet:{           tnip=Qtet->tnip;       break;  }
    case linearhex:{          tnip=Lhex->tnip;       break;  }
    case quadrhex:{           tnip=Qhex->tnip;       break;  }
    case lineartetrot:{       tnip=Ltetrot->tnip;    break;  }
    case linearhexrot:{       tnip=Lhexrot->tnip;    break;  }
    case linearwed:{          tnip=Lwed->tnip;       break;  }
    case quadrwed:{           tnip=Qwed->tnip;       break;  }
    case hexintface:{         tnip=Hexifc->tnip;     break; }
    case particleelem:{       tnip=1;                break;  }
    case tetralatt:{          tnip=Tlatt->tnip;      break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return tnip;
}



/**
  Function returns total number of integration points on element.
   
  @param eid - element id

  @return The function returns total number of integration points on element. 
   
  Created by Tomas Koudelka, 30.12.2001
*/
long mechtop::give_totnip (long eid)
{
  if (eid < Mt->ne-1)
    return elements[eid+1].ipp[0][0] - elements[eid].ipp[0][0];

  return Mm->tnip-elements[eid].ipp[0][0];
}


/**
  The function determines the maximum number of stress/strain components from all 
  elements used in the problem. The resulting value is strored in max_ncompstr and 
  it is referenced in the first order homogenization problems.

  @return The function returns the maximum dimension of elements in the problem.

  Created by Tomas Krejci, 24/10/2019  
*/
long mechtop::set_maxncompstr ()
{
  long i,aux,ret;

  for (i=0; i<ne; i++)
  {
    aux = give_tncomp (i);
    if (max_ncompstr < aux)
      max_ncompstr = aux;
  }
  ret = max_ncompstr;
  return ret;
}


/**
  Function returns number of strain/stress components in the given block.
  For most of elements it is equal to the total number of stress/strain components.
   
  @param eid - element id
  @param bid - block id
   
  @return The function returns number of stress/strain components

  Created by JK, 30.12.2001
*/
long mechtop::give_ncomp (long eid,long bid)
{
  long ncomp=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{
      //ncomp=4;
      ncomp=Bar2d->tncomp;
      break;  }
    case bar3d:{              ncomp=Bar3d->tncomp;      break;  }
    case barq2d:{             ncomp=Barq2d->tncomp;     break;  }
    case barq3d:{             ncomp=Barq3d->tncomp;     break;  }
    case beam2d:{             ncomp=Beam2d->tncomp;     break;  }
    case beam3d:{             ncomp=Beam3d->tncomp;     break;  }
    case beamg3d:{            ncomp=Beam3dg->tncomp;    break;  }
    case beam2dsp:{           ncomp=Spbeam2d->tncomp;   break;  }
    case subsoilbeam:{        ncomp=Sbeam->tncomp;      break;  }
    case spring_1:{           ncomp=Spring->tncomp;     break;  }
    case spring_2:{           ncomp=Spring->tncomp;     break;  }
    case spring_3:{           ncomp=Spring->tncomp;     break;  }
    case spring_4:{           ncomp=Spring->tncomp;     break;  }
    case spring_5:{           ncomp=Spring->tncomp;     break;  }
    case spring_6:{           ncomp=Spring->tncomp;     break;  }
    case planeelementlt:{     ncomp=Pelt->tncomp;       break;  }
    case planeelementqt:{     ncomp=Peqt->tncomp;       break;  }
    case planeelementrotlt:{  ncomp=Perlt->tncomp;      break;  }
    case planeelementlq:{     ncomp=Pelq->tncomp;       break;  }
    case planeelementqq:{     ncomp=Peqq->tncomp;       break;  }
    case planeelementrotlq:{  ncomp=Perlq->tncomp;      break;  }
    case planeelementsubqt:{  ncomp=Pesqt->tncomp;      break;  }
    case planequadinterface:{ ncomp=Pqifc->tncomp;      break;  }
    case cctel:{              ncomp=Cct->tncomp;        break;  }
    case dktel:{              ncomp=Dkt->tncomp;        break;  }
    case dstel:{              ncomp=Dst->tncomp;        break;  }
    case q4plateel:{          ncomp=Q4pl->tncomp;       break;  }
      //case argyristr:{          ncomp=Argtr->tncomp;      break;  }
    case argyristr:{          ncomp=Argtrpl->tncomp;    break;  }
    case quadkirch:{          ncomp=Qkirch->tncomp;     break;  }
    case dkqel:{              ncomp=Dkqelem->tncomp;        break;  }
    case subsoilplatetr:{     ncomp=Spltr->tncomp;      break;  }
    case subsoilplateq:{      ncomp=Splq->tncomp;       break;  }
    case axisymmlt:{          ncomp=Asymlt->tncomp;     break;  }
    case axisymmqt:{          ncomp=Asymqt->tncomp;     break;  }
    case axisymmlq:{          ncomp=Asymlq->tncomp;     break;  }
    case axisymmqq:{          ncomp=Asymqq->tncomp;     break;  }
    case axisymmcq:{          ncomp=Asymcq->tncomp;     break;  }
    case axisymmlqintface:{   ncomp=Asymlqifc->tncomp;  break;  }
    case shelltrelem:{
      if (bid==0 || bid==1)   ncomp=Perlt->tncomp;
      if (bid==2)             ncomp=Dkt->tncomp;
      break;
    }
    case shelltrmelem:{
      if (bid==0 || bid==1)   ncomp=Perlt->tncomp;
      if (bid==2 || bid==3)   ncomp=Cct->tncomp;
      break;
    }
    case shellqelem:{
      if (bid==0 || bid==1)   ncomp=Perlq->tncomp;
      if (bid==2 || bid==3)   ncomp=Dkqelem->tncomp;
      break;
    }
    case lineartet:{          ncomp=Ltet->tncomp;       break;  }
    case quadrtet:{           ncomp=Qtet->tncomp;       break;  }
    case linearhex:{          ncomp=Lhex->tncomp;       break;  }
    case quadrhex:{           ncomp=Qhex->tncomp;       break;  }
    case lineartetrot:{       ncomp=Ltetrot->tncomp;    break;  }
    case linearhexrot:{       ncomp=Lhexrot->tncomp;    break;  }
    case linearwed:{          ncomp=Lwed->tncomp;       break;  }
    case quadrwed:{           ncomp=Qwed->tncomp;       break;  }
    case hexintface:{         ncomp=Hexifc->tncomp;     break; }
    case particleelem:{       ncomp=Pelem->ndofe;       break;  }
    case tetralatt:{          ncomp=Tlatt->tncomp;      break;  }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
  return ncomp;
}



/**
  Function returns the total number of strain/stress components of the given element.
   
  @param eid - element id

  @return The function returns number of stress/strain components
     
  Created by JK, 1.10.2008
  Modified by TKo, 03.2018
*/
long mechtop::give_tncomp (long eid)
{
  long tncomp;
  elemtype te;

  te = give_elem_type (eid);
  tncomp = give_tncomp(te);
 
  return tncomp;
}



/**
  Function returns the total number of strain/stress components of the given element.
   
  @param eid - element id

  @return The function returns number of stress/strain components
     
  Created by JK, 1.10.2008
*/
long mechtop::give_tncomp (elemtype te)
{
  long tncomp = 0;

  switch (te){
    case bar2d:{              tncomp=Bar2d->tncomp;      break;  }
    case bar3d:{              tncomp=Bar3d->tncomp;      break;  }
    case barq2d:{             tncomp=Barq2d->tncomp;     break;  }
    case barq3d:{             tncomp=Barq3d->tncomp;     break;  }
    case beam2d:{             tncomp=Beam2d->tncomp;     break;  }
    case beam3d:{             tncomp=Beam3d->tncomp;     break;  }
    case beamg3d:{            tncomp=Beam3dg->tncomp;    break;  }
    case beam2dsp:{           tncomp=Spbeam2d->tncomp;   break;  }
    case subsoilbeam:{        tncomp=Sbeam->tncomp;      break;  }
    case spring_1:{           tncomp=Spring->tncomp;     break;  }
    case spring_2:{           tncomp=Spring->tncomp;     break;  }
    case spring_3:{           tncomp=Spring->tncomp;     break;  }
    case spring_4:{           tncomp=Spring->tncomp;     break;  }
    case spring_5:{           tncomp=Spring->tncomp;     break;  }
    case spring_6:{           tncomp=Spring->tncomp;     break;  }
    case planeelementlt:{     tncomp=Pelt->tncomp;       break;  }
    case planeelementqt:{     tncomp=Peqt->tncomp;       break;  }
    case planeelementrotlt:{  tncomp=Perlt->tncomp;      break;  }
    case planeelementlq:{     tncomp=Pelq->tncomp;       break;  }
    case planeelementqq:{     tncomp=Peqq->tncomp;       break;  }
    case planeelementrotlq:{  tncomp=Perlq->tncomp;      break;  }
    case planeelementsubqt:{  tncomp=Pesqt->tncomp;      break;  }
    case planequadinterface:{ tncomp=Pqifc->tncomp;      break;  }
    case cctel:{              tncomp=Cct->tncomp;        break;  }
    case dktel:{              tncomp=Dkt->tncomp;        break;  }
    case dstel:{              tncomp=Dst->tncomp;        break;  }
    case q4plateel:{          tncomp=Q4pl->tncomp;       break;  }
      //case argyristr:{          tncomp=Argtr->tncomp;      break;  }
    case argyristr:{          tncomp=Argtrpl->tncomp;    break;  }
    case quadkirch:{          tncomp=Qkirch->tncomp;     break;  }
    case dkqel:{              tncomp=Dkqelem->tncomp;    break;  }
    case subsoilplatetr:{     tncomp=Spltr->tncomp;      break;  }
    case subsoilplateq:{      tncomp=Splq->tncomp;       break;  }
    case axisymmlt:{          tncomp=Asymlt->tncomp;     break;  }
    case axisymmqt:{          tncomp=Asymqt->tncomp;     break;  }
    case axisymmlq:{          tncomp=Asymlq->tncomp;     break;  }
    case axisymmqq:{          tncomp=Asymqq->tncomp;     break;  }
    case axisymmcq:{          tncomp=Asymcq->tncomp;     break;  }
    case axisymmlqintface:{   tncomp=Asymlqifc->tncomp;  break;  }
    case shelltrelem:{        tncomp=Shtr->tncomp;       break;  }
    case shelltrmelem:{       tncomp=Shtrm->tncomp;      break;  }
    case shellqelem:{         tncomp=Shq->tncomp;        break;  }
    case lineartet:{          tncomp=Ltet->tncomp;       break;  }
    case quadrtet:{           tncomp=Qtet->tncomp;       break;  }
    case linearhex:{          tncomp=Lhex->tncomp;       break;  }
    case quadrhex:{           tncomp=Qhex->tncomp;       break;  }
    case lineartetrot:{       tncomp=Ltetrot->tncomp;    break;  }
    case linearhexrot:{       tncomp=Lhexrot->tncomp;    break;  }
    case linearwed:{          tncomp=Lwed->tncomp;       break;  }
    case quadrwed:{           tncomp=Qwed->tncomp;       break;  }
    case hexintface:{         tncomp=Hexifc->tncomp;     break; }
    case particleelem:{       tncomp=Pelem->ndofe;       break;  }
    case tetralatt:{          tncomp=Tlatt->tncomp;      break;  }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
  return tncomp;
}



/**
  Function returns number of edges of the element.
   
  @param eid - element id
   
  @return The function returns number of edges.

  Created by JK, 30.12.2001
*/
long mechtop::give_ned (long eid)
{
  long ned=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar3d:{              ned=0;             break;  }
    case beam2d:{             ned=1;             break;  }
    case beam3d:{             ned=1;             break;  }
    case planeelementlt:{     ned=Pelt->ned;     break;  }
    case planeelementqt:{     ned=Peqt->ned;     break;  }
    case planeelementrotlt:{  ned=Perlt->ned;    break;  }
    case planeelementlq:{     ned=Pelq->ned;     break;  }
    case planeelementqq:{     ned=Peqq->ned;     break;  }
    case planeelementrotlq:{  ned=Perlq->ned;    break;  }
    case planeelementsubqt:{  ned=Pesqt->ned;    break;  }
    case cctel:{              ned=Cct->ned;      break;  }
    case dktel:{              ned=Dkt->ned;      break;  }
    case dstel:{              ned=Dst->ned;      break;  }
    case q4plateel:{          ned=Q4pl->ned;     break;  }
      //case argyristr:{          ned=Argtr->ned;   break;  }
    case argyristr:{          ned=Argtrpl->ned;  break;  }
    case dkqel:{              ned=Dkqelem->ned;  break;  }
    case subsoilplatetr:{     ned=Spltr->ned;    break;  }
    case subsoilplateq:{      ned=Splq->ned;     break;  }
    case axisymmlt:{          ned=Asymlt->ned;   break;  }
    case axisymmqt:{          ned=Asymqt->ned;   break;  }
    case axisymmlq:{          ned=Asymlq->ned;   break;  }
    case axisymmqq:{          ned=Asymqq->ned;   break;  }
    case axisymmcq:{          ned=Asymcq->ned;   break;  }
    case axisymmlqintface:{   ned=Asymlqifc->ned;break;  }
    case shelltrelem:{        ned=Shtr->ned;     break;  }
    case shelltrmelem:{       ned=Shtrm->ned;    break;  }
    case shellqelem:{         ned=Shq->ned;      break;  }
    case lineartet:{          ned=Ltet->ned;     break;  }
    case quadrtet:{           ned=Qtet->ned;     break;  }
    case linearhex:{          ned=Lhex->ned;     break;  }
    case quadrhex:{           ned=Qhex->ned;     break;  }
    case lineartetrot:{       ned=Ltetrot->ned;  break;  }
    case linearhexrot:{       ned=Lhexrot->ned;  break;  }
    case linearwed:{          ned=Lwed->ned;     break;  }
    case quadrwed:{           ned=Qwed->ned;     break;  }
    case particleelem:{                          break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return ned;
}



/**
  Function returns number of nodes on one element edge.
   
  @param eid - number of element
   
  @return The function returns number of nodes on element edge

  Created by JK, 30.12.2001
*/
long mechtop::give_nned (long eid)
{
  long nned=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar3d:{               nned=0;               break;  }
    case beam2d:{              nned=2;               break;  }
    case beam3d:{              nned=2;               break;  }
    case planeelementlt:{      nned=Pelt->nned;      break;  }
    case planeelementqt:{      nned=Peqt->nned;      break;  }
    case planeelementrotlt:{   nned=Perlt->nned;     break;  }
    case planeelementlq:{      nned=Pelq->nned;      break;  }
    case planeelementqq:{      nned=Peqq->nned;      break;  }
    case planeelementrotlq:{   nned=Perlq->nned;     break;  }
    case planeelementsubqt:{   nned=Pesqt->nned;     break;  }
    case cctel:{               nned=Cct->nned;       break;  }
    case dktel:{               nned=Dkt->nned;       break;  }
    case dstel:{               nned=Dst->nned;       break;  }
    case q4plateel:{           nned=Q4pl->nned;      break;  }
      //case argyristr:{           nned=Argtr->nned;    break;  }
    case argyristr:{           nned=Argtrpl->nned;   break;  }
    case dkqel:{               nned=Dkqelem->nned;   break;  }
    case subsoilplatetr:{      nned=Spltr->nned;     break;  }
    case subsoilplateq:{       nned=Splq->nned;      break;  }
    case axisymmlt:{           nned=Asymlt->nned;    break;  }
    case axisymmqt:{           nned=Asymqt->nned;    break;  }
    case axisymmlq:{           nned=Asymlq->nned;    break;  }
    case axisymmqq:{           nned=Asymqq->nned;    break;  }
    case axisymmcq:{           nned=Asymcq->nned;    break;  }
    case axisymmlqintface:{    nned=Asymlqifc->nned; break;  }
    case shelltrelem:{         nned=Shtr->nned;      break;  }
    case shelltrmelem:{        nned=Shtrm->nned;     break;  }
    case shellqelem:{          nned=Shq->nned;       break;  }
    case lineartet:{           nned=Ltet->nned;      break;  }
    case quadrtet:{            nned=Qtet->nned;      break;  }
    case linearhex:{           nned=Lhex->nned;      break;  }
    case quadrhex:{            nned=Qhex->nned;      break;  }
    case lineartetrot:{        nned=Ltetrot->nned;   break;  }
    case linearhexrot:{        nned=Lhexrot->nned;   break;  }
    case linearwed:{           nned=Lwed->nned;      break;  }
    case quadrwed:{            nned=Qwed->nned;      break;  }
    case particleelem:{                              break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return nned;
}



/**
  Function returns number of surfaces of the element.
   
  @param eid - element id
   
  @return The function returns number of element surfaces.

  Created by JK, 29.4.2004
*/
long mechtop::give_nsurf (long eid)
{
  long nsurf=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case planeelementrotlq:{   nsurf=Perlq->nsurf;      break;  }
  case dktel:{               nsurf=Dkt->nsurf;        break;  }
  case quadkirch:{           nsurf=Qkirch->nsurf;     break;  }
  case dkqel:{               nsurf=Dkqelem->nsurf;    break;  }
  case shelltrelem:{         nsurf=Shtr->nsurf;       break;  }
  case shelltrmelem:{        nsurf=Shtrm->nsurf;      break;  }
  case shellqelem:{          nsurf=Shq->nsurf;        break;  }
  case lineartet:{           nsurf=Ltet->nsurf;       break;  }
  case quadrtet:{            nsurf=Qtet->nsurf;       break;  }
  case linearhex:{           nsurf=Lhex->nsurf;       break;  }
  case quadrhex:{            nsurf=Qhex->nsurf;       break;  }
  case lineartetrot:{        nsurf=Ltetrot->nsurf;    break;  }
  case linearhexrot:{        nsurf=Lhexrot->nsurf;    break;  }
  case linearwed:{           nsurf=Lwed->nsurf;       break;  }
  case quadrwed:{            nsurf=Qwed->nsurf;       break;  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nsurf;
}



/**
  Function returns number of surfaces of the element.
   
  @param eid - element id
   
  @return The function returns number of nodes on element surfaces

  Created by JK, 29.4.2004
*/
long mechtop::give_nnsurf (long eid)
{
  long nnsurf=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
  case planeelementrotlq:{    nnsurf=Perlq->nnsurf;      break;  }
  case dktel:{                nnsurf=Dkt->nnsurf;        break;  }
  case quadkirch:{            nnsurf=Qkirch->nnsurf;     break;  }
  case dkqel:{                nnsurf=Dkqelem->nnsurf;    break;  }
  case shelltrelem:{          nnsurf=Shtr->nnsurf;       break;  }
  case shelltrmelem:{         nnsurf=Shtrm->nnsurf;      break;  }
  case shellqelem:{           nnsurf=Shq->nnsurf;        break;  }
  case lineartet:{            nnsurf=Ltet->nnsurf;       break;  }
  case quadrtet:{             nnsurf=Qtet->nnsurf;       break;  }
  case linearhex:{            nnsurf=Lhex->nnsurf;       break;  }
  case quadrhex:{             nnsurf=Qhex->nnsurf;       break;  }
  case lineartetrot:{         nnsurf=Ltetrot->nnsurf;    break;  }
  case linearhexrot:{         nnsurf=Lhexrot->nnsurf;    break;  }
  case linearwed:{            nnsurf=Lwed->nnsurf;       break;  }
  case quadrwed:{             nnsurf=Qwed->nnsurf;       break;  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return nnsurf;
}



/**
  Function returns number of approximated functions.
   
  @param eid - number of element

  @return The function returns number of approximation functions.
   
  Created by JK, 30.12.2001
*/
long mechtop::give_napfun (long eid)
{
  long napfun=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              napfun=Bar2d->napfun;      break;  }
    case bar3d:{              napfun=Bar3d->napfun;      break;  }
    case barq2d:{             napfun=Barq2d->napfun;     break;  }
    case barq3d:{             napfun=Barq3d->napfun;     break;  }
    case beam2d:{             napfun=Beam2d->napfun;     break;  }
    case beam3d:{             napfun=Beam3d->napfun;     break;  }
    case beamg3d:{            napfun=Beam3dg->napfun;    break;  }
    case subsoilbeam:{        napfun=Sbeam->napfun;      break;  }
    case beam2dsp:{           napfun=Spbeam2d->napfun;   break;  }
    case spring_1:{           napfun=Spring->napfun;     break;  }
    case spring_2:{           napfun=Spring->napfun;     break;  }
    case spring_3:{           napfun=Spring->napfun;     break;  }
    case spring_4:{           napfun=Spring->napfun;     break;  }
    case spring_5:{           napfun=Spring->napfun;     break;  }
    case spring_6:{           napfun=Spring->napfun;     break;  }
    case planeelementlt:{     napfun=Pelt->napfun;       break;  }
    case planeelementqt:{     napfun=Peqt->napfun;       break;  }
    case planeelementrotlt:{  napfun=Perlt->napfun;      break;  }
    case planeelementlq:{     napfun=Pelq->napfun;       break;  }
    case planeelementqq:{     napfun=Peqq->napfun;       break;  }
    case planeelementrotlq:{  napfun=Perlq->napfun;      break;  }
    case planeelementsubqt:{  napfun=Pesqt->napfun;      break;  }
    case cctel:{              napfun=Cct->napfun;        break;  }
    case dktel:{              napfun=Dkt->napfun;        break;  }
    case dstel:{              napfun=Dst->napfun;        break;  }
    case q4plateel:{          napfun=Q4pl->napfun;       break;  }
      //case argyristr:{          napfun=Argtr->napfun;      break;  }
    case argyristr:{          napfun=Argtrpl->napfun;    break;  }
    case quadkirch:{          napfun=Qkirch->napfun;     break;  }
    case dkqel:{              napfun=Dkqelem->napfun;        break;  }
    case subsoilplatetr:{     napfun=Spltr->napfun;      break;  }
    case subsoilplateq:{      napfun=Splq->napfun;       break;  }
    case axisymmlt:{          napfun=Asymlt->napfun;     break;  }
    case axisymmqt:{          napfun=Asymqt->napfun;     break;  }
    case axisymmlq:{          napfun=Asymlq->napfun;     break;  }
    case axisymmqq:{          napfun=Asymqq->napfun;     break;  }
    case axisymmcq:{          napfun=Asymcq->napfun;     break;  }
    case axisymmlqintface:{   napfun=Asymlqifc->napfun;  break;  }
    case shelltrelem:{        napfun=Shtr->napfun;       break;  }
    case shelltrmelem:{       napfun=Shtrm->napfun;      break;  }
    case shellqelem:{         napfun=Shq->napfun;        break;  }
    case lineartet:{          napfun=Ltet->napfun;       break;  }
    case quadrtet:{           napfun=Qtet->napfun;       break;  }
    case linearhex:{          napfun=Lhex->napfun;       break;  }
    case quadrhex:{           napfun=Qhex->napfun;       break;  }
    case lineartetrot:{       napfun=Ltetrot->napfun;    break;  }
    case linearhexrot:{       napfun=Lhexrot->napfun;    break;  }
    case linearwed:{          napfun=Lwed->napfun;       break;  }
    case quadrwed:{           napfun=Qwed->napfun;       break;  }
    case particleelem:{                                  break;  }
    case tetralatt:{          napfun=Tlatt->napfun;      break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return napfun;
}



/**
  Function returns dimension of the given element.
   
  @param eid - number of element
   
  @retval 1 - for 1D
  @retval 2 - for 2D
  @retval 3 - for 3D

  Created by JK, 23.2.2002
*/
long mechtop::give_dimension (long eid)
{
  long dim=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              dim=1;  break;  }
    case bar3d:{              dim=1;  break;  }
    case barq2d:{             dim=1;  break;  }
    case barq3d:{             dim=1;  break;  }
    case beam2d:{             dim=1;  break;  }
    case beam3d:{             dim=1;  break;  }
    case beamg3d:{            dim=1;  break;  }
    case subsoilbeam:{        dim=1;  break;  }
    case beam2dsp:{           dim=1;  break;  }
    case spring_1:{           dim=1;  break;  }
    case spring_2:{           dim=1;  break;  }
    case spring_3:{           dim=1;  break;  }
    case spring_4:{           dim=1;  break;  }
    case spring_5:{           dim=1;  break;  }
    case spring_6:{           dim=1;  break;  }
    case planeelementlt:{     dim=2;  break;  }
    case planeelementqt:{     dim=2;  break;  }
    case planeelementrotlt:{  dim=2;  break;  }
    case planeelementlq:{     dim=2;  break;  }
    case planeelementqq:{     dim=2;  break;  }
    case planeelementrotlq:{  dim=2;  break;  }
    case planeelementsubqt:{  dim=2;  break;  }
    case cctel:{              dim=2;  break;  }
    case dktel:{              dim=2;  break;  }
    case dstel:{              dim=2;  break;  }
    case q4plateel:{          dim=2;  break;  }
    case argyristr:{          dim=2;  break;  }
    case quadkirch:{          dim=2;  break;  }
    case dkqel:{              dim=2;  break;  }
    case axisymmlt:{          dim=2;  break;  }
    case axisymmqt:{          dim=2;  break;  }
    case axisymmlq:{          dim=2;  break;  }
    case axisymmqq:{          dim=2;  break;  }
    case axisymmcq:{          dim=2;  break;  }
    case axisymmlqintface:{   dim=1;  break;  }
    case subsoilplatetr:{     dim=2;  break;  }
    case subsoilplateq:{      dim=2;  break;  }
    case shelltrelem:{        dim=2;  break;  }
    case shelltrmelem:{       dim=2;  break;  }
    case shellqelem:{         dim=2;  break;  }
    case lineartet:{          dim=3;  break;  }
    case quadrtet:{           dim=3;  break;  }
    case linearhex:{          dim=3;  break;  }
    case quadrhex:{           dim=3;  break;  }
    case lineartetrot:{       dim=3;  break;  }
    case linearhexrot:{       dim=3;  break;  }
    case linearwed:{          dim=3;  break;  }
    case quadrwed:{           dim=3;  break;  }
    case hexintface:{         dim=3;  break; }
    case particleelem:{               break;  }
    case tetralatt:{          dim=3;  break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return dim;
}



/**
  The function determines the maximum dimension of elements used in the problem.

  @return The function returns the maximum dimension of elements in the problem.

  Created by Tomas Koudelka, 08.2016  
*/
long mechtop::give_maxdimension ()
{
  long i, aux, ret = 0L;

  for (i=0L; i<ne; i++)
  {
    aux = give_dimension(i);
    if (ret < aux)
      ret = aux;
  }
  return ret;
}



/**
  Function returns number of blocks in characteristic matrices.
   
  @param eid - element id
   
  @return The function returns number of blocks.

  Created by JK, 8.5.2002
*/
long mechtop::give_nb (long eid)
{
  long nb=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{              nb=Bar2d->nb;      break;  }
    case bar3d:{              nb=Bar3d->nb;      break;  }
    case barq2d:{             nb=Barq2d->nb;     break;  }
    case barq3d:{             nb=Barq3d->nb;     break;  }
    case beam2d:{             nb=Beam2d->nb;     break;  }
    case beam3d:{             nb=Beam3d->nb;     break;  }
    case beamg3d:{            nb=Beam3dg->nb;    break;  }
    case subsoilbeam:{        nb=Sbeam->nb;      break;  }
    case beam2dsp:{           nb=Spbeam2d->nb;   break;  }
    case spring_1:{           nb=Spring->nb;     break;  }
    case spring_2:{           nb=Spring->nb;     break;  }
    case spring_3:{           nb=Spring->nb;     break;  }
    case spring_4:{           nb=Spring->nb;     break;  }
    case spring_5:{           nb=Spring->nb;     break;  }
    case spring_6:{           nb=Spring->nb;     break;  }
    case planeelementlt:{     nb=Pelt->nb;       break;  }
    case planeelementqt:{     nb=Peqt->nb;       break;  }
    case planeelementrotlt:{  nb=Perlt->nb;      break;  }
    case planeelementlq:{     nb=Pelq->nb;       break;  }
    case planeelementqq:{     nb=Peqq->nb;       break;  }
    case planeelementrotlq:{  nb=Perlq->nb;      break;  }
    case planeelementsubqt:{  nb=Pesqt->nb;      break;  }
    case planequadinterface:{ nb=Pqifc->nb;      break;  }
    case cctel:{              nb=Cct->nb;        break;  }
    case dktel:{              nb=Dkt->nb;        break;  }
    case dstel:{              nb=Dst->nb;        break;  }
    case q4plateel:{          nb=Q4pl->nb;       break;  }
      //  case argyristr:{          nb=Argtr->nb;      break;  }
    case argyristr:{          nb=Argtrpl->nb;    break;  }
    case quadkirch:{          nb=Qkirch->nb;     break;  }
    case dkqel:{              nb=Dkqelem->nb;        break;  }
    case subsoilplatetr:{     nb=Spltr->nb;      break;  }
    case subsoilplateq:{      nb=Splq->nb;       break;  }
    case axisymmlt:{          nb=Asymlt->nb;     break;  }
    case axisymmqt:{          nb=Asymqt->nb;     break;  }
    case axisymmlq:{          nb=Asymlq->nb;     break;  }
    case axisymmqq:{          nb=Asymqq->nb;     break;  }
    case axisymmcq:{          nb=Asymcq->nb;     break;  }
    case axisymmlqintface:{   nb=Asymlqifc->nb;  break;  }
    case shelltrelem:{        nb=Shtr->nb;       break;  }
    case shelltrmelem:{       nb=Shtrm->nb;      break;  }
    case shellqelem:{         nb=Shq->nb;        break;  }
    case lineartet:{          nb=Ltet->nb;       break;  }
    case quadrtet:{           nb=Qtet->nb;       break;  }
    case linearhex:{          nb=Lhex->nb;       break;  }
    case quadrhex:{           nb=Qhex->nb;       break;  }
    case lineartetrot:{       nb=Ltetrot->nb;    break;  }
    case linearhexrot:{       nb=Lhexrot->nb;    break;  }
    case linearwed:{          nb=Lwed->nb;       break;  }
    case quadrwed:{           nb=Qwed->nb;       break;  }
    case hexintface:{         nb=Hexifc->nb;     break;  }
    case particleelem:{       nb=Pelem->nb;      break;  }
    case tetralatt:{          nb=Tlatt->nb;      break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return nb;
}


/**
  Function returns number of blocks in characteristic matrices.
   
  @param eid - element id
   
  @return The function returns number of blocks.

  Created by JK, 8.5.2002
*/
long mechtop::give_nb_te (elemtype te)
{
  long nb=0;
  
  switch (te){
    case bar2d:{              nb=Bar2d->nb;      break;  }
    case bar3d:{              nb=Bar3d->nb;      break;  }
    case barq2d:{             nb=Barq2d->nb;     break;  }
    case barq3d:{             nb=Barq3d->nb;     break;  }
    case beam2d:{             nb=Beam2d->nb;     break;  }
    case beam3d:{             nb=Beam3d->nb;     break;  }
    case beamg3d:{            nb=Beam3dg->nb;    break;  }
    case subsoilbeam:{        nb=Sbeam->nb;      break;  }
    case beam2dsp:{           nb=Spbeam2d->nb;   break;  }
    case spring_1:{           nb=Spring->nb;     break;  }
    case spring_2:{           nb=Spring->nb;     break;  }
    case spring_3:{           nb=Spring->nb;     break;  }
    case spring_4:{           nb=Spring->nb;     break;  }
    case spring_5:{           nb=Spring->nb;     break;  }
    case spring_6:{           nb=Spring->nb;     break;  }
    case planeelementlt:{     nb=Pelt->nb;       break;  }
    case planeelementqt:{     nb=Peqt->nb;       break;  }
    case planeelementrotlt:{  nb=Perlt->nb;      break;  }
    case planeelementlq:{     nb=Pelq->nb;       break;  }
    case planeelementqq:{     nb=Peqq->nb;       break;  }
    case planeelementrotlq:{  nb=Perlq->nb;      break;  }
    case planeelementsubqt:{  nb=Pesqt->nb;      break;  }
    case planequadinterface:{ nb=Pqifc->nb;      break;  }
    case cctel:{              nb=Cct->nb;        break;  }
    case dktel:{              nb=Dkt->nb;        break;  }
    case dstel:{              nb=Dst->nb;        break;  }
    case q4plateel:{          nb=Q4pl->nb;       break;  }
      //  case argyristr:{          nb=Argtr->nb;      break;  }
    case argyristr:{          nb=Argtrpl->nb;    break;  }
    case quadkirch:{          nb=Qkirch->nb;     break;  }
    case dkqel:{              nb=Dkqelem->nb;    break;  }
    case subsoilplatetr:{     nb=Spltr->nb;      break;  }
    case subsoilplateq:{      nb=Splq->nb;       break;  }
    case axisymmlt:{          nb=Asymlt->nb;     break;  }
    case axisymmqt:{          nb=Asymqt->nb;     break;  }
    case axisymmlq:{          nb=Asymlq->nb;     break;  }
    case axisymmqq:{          nb=Asymqq->nb;     break;  }
    case axisymmcq:{          nb=Asymcq->nb;     break;  }
    case axisymmlqintface:{   nb=Asymlqifc->nb;  break;  }
    case shelltrelem:{        nb=Shtr->nb;       break;  }
    case shelltrmelem:{       nb=Shtrm->nb;      break;  }
    case shellqelem:{         nb=Shq->nb;        break;  }
    case lineartet:{          nb=Ltet->nb;       break;  }
    case quadrtet:{           nb=Qtet->nb;       break;  }
    case linearhex:{          nb=Lhex->nb;       break;  }
    case quadrhex:{           nb=Qhex->nb;       break;  }
    case lineartetrot:{       nb=Ltetrot->nb;    break;  }
    case linearhexrot:{       nb=Lhexrot->nb;    break;  }
    case linearwed:{          nb=Lwed->nb;       break;  }
    case quadrwed:{           nb=Qwed->nb;       break;  }
    case hexintface:{         nb=Hexifc->nb;     break;  }
    case particleelem:{       nb=Pelem->nb;      break;  }
    case tetralatt:{          nb=Tlatt->nb;      break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return nb;
}


/**
  Function returns number of layeres of layered element.
   
  @param eid - element id
 
  @return The function returns number of layers,
  
  Created by JK, 7.12.2002
*/
/*
long mechtop::give_nl (long eid)
{
  return Gtm->give_nl (eid);
}
*/



/**
  The function returns stress/strain state indicator.
  
  @param eid - element id
  @param bi - block id

  @retval bar - for 1D stress state
  @retval plbeam - for 2D beam stress state
  @retval spacebeam - for 3D beam stress state
  @retval planestress - for 2D plane stress state
  @retval planestrain - for 2D plane strain state
  @retval platek      - plate stress state
  @retval plates      - plate stress state
  @retval axisymm     - for 2D axisymmetric stress state
  @retval shell       - for 3D shell stress state
  @retval spacestress - for fully 3D space stress state
  
  Created by JK, 7.12.2002
*/
strastrestate mechtop::give_ssst (long eid,long bi)
{
  strastrestate ssst=strastrestate(0);
  elemtype te;
  
  te = give_elem_type (eid);
  switch (te){
    case bar2d:{              ssst = Bar2d->ssst;         break;  }
    case bar3d:{              ssst = Bar3d->ssst;         break;  }
    case barq2d:{             ssst = Barq2d->ssst;        break;  }
    case barq3d:{             ssst = Barq3d->ssst;        break;  }
    case beam2d:{             ssst = Beam2d->ssst;        break;  }
    case beam3d:{             ssst = Beam3d->ssst;        break;  }
    case beamg3d:{            ssst = Beam3dg->ssst;       break;  }
    case subsoilbeam:{        ssst = Sbeam->ssst;         break;  }
    case beam2dsp:{           ssst = Spbeam2d->ssst;      break;  }
    case spring_1:{           ssst = Spring->ssst;        break;  }
    case spring_2:{           ssst = Spring->ssst;        break;  }
    case spring_3:{           ssst = Spring->ssst;        break;  }
    case spring_4:{           ssst = Spring->ssst;        break;  }
    case spring_5:{           ssst = Spring->ssst;        break;  }
    case spring_6:{           ssst = Spring->ssst;        break;  }
    case planeelementlt:{     ssst = elements[eid].ssst;  break;  }
    case planeelementqt:{     ssst = elements[eid].ssst;  break;  }
    case planeelementrotlt:{  ssst = elements[eid].ssst;  break;  }
    case planeelementlq:{     ssst = elements[eid].ssst;  break;  }
    case planeelementqq:{     ssst = elements[eid].ssst;  break;  }
    case planeelementrotlq:{  ssst = elements[eid].ssst;  break;  }
    case planeelementsubqt:{  ssst = elements[eid].ssst;  break;  }
    case planequadinterface:{ ssst = Pqifc->ssst;         break;  }
    case cctel:{              ssst = Cct->ssst;           break;  }
    case dktel:{              ssst = Dkt->ssst;           break;  }
    case dstel:{              ssst = Dst->ssst;           break;  }
    case q4plateel:{          ssst = Q4pl->ssst;          break;  }
      //  case argyristr:{          ssst = Argtr->ssst;         break;  }
    case argyristr:{          ssst = Argtrpl->ssst;       break;  }
    case quadkirch:{          ssst = Qkirch->ssst;        break;  }
    case dkqel:{              ssst = Dkqelem->ssst;           break;  }
    case subsoilplatetr:{     ssst = Spltr->ssst;         break;  }
    case subsoilplateq:{      ssst = Splq->ssst;          break;  }
    case axisymmlt:{          ssst = Asymlt->ssst;        break;  }
    case axisymmqt:{          ssst = Asymqt->ssst;        break;  }
    case axisymmlq:{          ssst = Asymlq->ssst;        break;  }
    case axisymmqq:{          ssst = Asymqq->ssst;        break;  }
    case axisymmcq:{          ssst = Asymcq->ssst;        break;  }
    case axisymmlqintface:{   ssst = Asymlqifc->ssst;     break;  }
    case shelltrelem:{
      if (bi==0 || bi==1)  ssst = (strastrestate) planestress;
      if (bi==2)           ssst = (strastrestate) platek;
      if (bi==3)           ssst = (strastrestate) shell;
      break;
    }
    case shelltrmelem:{
      if (bi==0 || bi==1)  ssst = (strastrestate) planestress;
      if (bi==2 || bi==3)  ssst = (strastrestate) plates;
      break;
    }
    case shellqelem:{
      if (bi==0 || bi==1)  ssst = (strastrestate) planestress;
      if (bi==2)  ssst = (strastrestate) platek;
      break;  }
    case lineartet:{          ssst = Ltet->ssst;          break;  }
    case quadrtet:{           ssst = Qtet->ssst;          break;  }
    case linearhex:{          ssst = Lhex->ssst;          break;  }
    case quadrhex:{           ssst = Qhex->ssst;          break;  }
    case lineartetrot:{       ssst = Ltetrot->ssst;       break;  }
    case linearhexrot:{       ssst = Lhexrot->ssst;       break;  }
    case linearwed:{          ssst = Lwed->ssst;          break;  }
    case quadrwed:{           ssst = Qwed->ssst;          break;  }
    case hexintface:{         ssst = Hexifc->ssst;        break;  }
    case particleelem:{                                   break;  }
    case tetralatt:{          ssst = Tlatt->ssst;         break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return ssst;
}



/**
  Function returns degree of polynomial expansion.
   
  @param eid - number of element
   
  @return The function returns degree of appproximation. 
  
  Created by JK,
*/
long mechtop::give_degree (long eid)
{
  long deg=0;
  elemtype te;
  
  te = give_elem_type (eid);
  
  switch (te){
    case bar2d:{                deg=1;  break;  }
    case bar3d:{                deg=1;  break;  }
    case barq2d:{               deg=2;  break;  }
    case barq3d:{               deg=2;  break;  }
    case beam2d:{               deg=1;  break;  }
    case beam3d:{               deg=1;  break;  }
    case beamg3d:{              deg=3;  break;  }
    case spring_1:{             deg=1;  break;  }
    case spring_2:{             deg=1;  break;  }
    case spring_3:{             deg=1;  break;  }
    case spring_4:{             deg=1;  break;  }
    case spring_5:{             deg=1;  break;  }
    case spring_6:{             deg=1;  break;  }
    case subsoilbeam:{          deg=1;  break;  }
    case beam2dsp:{             deg=1;  break;  }
    case planeelementlt:{       deg=1;  break;  }
    case planeelementqt:{       deg=2;  break;  }
    case planeelementrotlt:{    deg=1;  break;  }
    case planeelementlq:{       deg=1;  break;  }
    case planeelementqq:{       deg=2;  break;  }
    case planeelementrotlq:{    deg=1;  break;  }
    case planeelementsubqt:{    deg=2;  break;  }
    case cctel:{                deg=1;  break;  }
    case dktel:{                deg=1;  break;  }
    case dstel:{                deg=1;  break;  }
    case q4plateel:{            deg=1;  break;  }
    case argyristr:{            deg=3;  break;  }
    case subsoilplatetr:{       deg=1;  break;  }
    case subsoilplateq:{        deg=1;  break;  }
    case axisymmlt:{            deg=1;  break;  }
    case axisymmqt:{            deg=1;  break;  }
    case axisymmlq:{            deg=1;  break;  }
    case axisymmqq:{            deg=1;  break;  }
    case axisymmcq:{            deg=1;  break;  }
    case axisymmlqintface:{     deg=1;  break;  }
    case shelltrelem:{          deg=1;  break;  }
    case shelltrmelem:{         deg=1;  break;  }
    case shellqelem:{           deg=1;  break;  }
    case lineartet:{            deg=1;  break;  }
    case quadrtet:{             deg=2;  break;  }
    case linearhex:{            deg=1;  break;  }
    case quadrhex:{             deg=2;  break;  }
    case lineartetrot:{         deg=1;  break;  }
    case linearhexrot:{         deg=1;  break;  }
    case linearwed:{            deg=1;  break;  }
    case quadrwed:{             deg=2;  break;  }
    case hexintface:{           deg=1;  break;  }
    case particleelem:{                 break;  }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
  }
  return deg;
}



/**
  Function returns node numbers on required element edge.
   
  @param eid - element id
  @param edid - edge id
  @param nodes - array containing edge nodes (output)
   
  @return The function returns edge nodes in the parameter nodes.

  Created by JK, 5.3.2002
*/
void mechtop::give_edge_nodes (long eid,long edid,long *nodes)
{
  long nne;
  elemtype te;
  ivector nod;
  
  te = give_elem_type (eid);
  nne = give_nne (eid);
  reallocv (RSTCKIVEC(nne,nod));
  give_elemnodes (eid,nod);
  
  switch (te){
  case planeelementlt:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[0]; }
    break;
  }
  case planeelementqt:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1];  nodes[2]=nod[3]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2];  nodes[2]=nod[4]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[0];  nodes[2]=nod[5]; }
    break;
  }
  case planeelementrotlt:{
    if (edid==0){  nodes[0]=nod[1];  nodes[1]=nod[2]; }
    if (edid==1){  nodes[0]=nod[2];  nodes[1]=nod[0]; }
    if (edid==2){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    break;
  }
  case planeelementlq:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[3]; }
    if (edid==3){  nodes[0]=nod[3];  nodes[1]=nod[0]; }
    break;
  }
  case planeelementqq:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[4];  nodes[2]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[5];  nodes[2]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[6];  nodes[2]=nod[3]; }
    if (edid==3){  nodes[0]=nod[3];  nodes[1]=nod[7];  nodes[2]=nod[0]; }
    break;
  }
  case planeelementrotlq:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[3]; }
    if (edid==3){  nodes[0]=nod[3];  nodes[1]=nod[0]; }
    break;
  }
  case planeelementsubqt:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1];  nodes[2]=nod[3]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2];  nodes[2]=nod[4]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[0];  nodes[2]=nod[5]; }
    break;
  }
  case axisymmlt:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[0]; }
    break;
  }
  case axisymmqt:{
    if (edid==0){  nodes[0]=nod[0];  nodes[1]=nod[3];  nodes[1]=nod[1]; }
    if (edid==1){  nodes[0]=nod[1];  nodes[1]=nod[4];  nodes[1]=nod[2]; }
    if (edid==2){  nodes[0]=nod[2];  nodes[1]=nod[5];  nodes[1]=nod[0]; }
    break;
  }
  case linearhex:{
    if (edid==0){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==1){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==2){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==3){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==4){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==5){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==6){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==7){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==8){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==9){   nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==10){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    if (edid==11){  nodes[0]=nod[0];  nodes[1]=nod[1]; }
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   Function returns numbers of nodes defining appropriate element.
   
   There are two arrays for node numbers of elements in the code.
   If there is no hanging node in the problem solved, nodes
   of elements are stored only in the objects Gtm->gelements.
   If there are hanging nodes, the elements attached to the
   hanging nodes contain the original slave nodes in the array
   Mt->elements.nodes while the array Gtm->gelements.nodes
   stores master nodes.
   
   This function returns the original node numbers
   i.e. the slave node numbers, it is used e.g. in graphic output
   
   @param[in] eid - element id
   @param[out] nodes - array containing nodes
   
   @return The function returns node numbers in the parameter nodes.
   
   Created by JK, 22.7.2001
*/
void mechtop::give_elemnodes (long eid,ivector &nodes)
{
  Gtm->give_original_nodes (eid,nodes);
}



/**
  Function extracts all code numbers of actual element.
   
  @param eid - element id
  @param cn - array containing all code numbers on element (output)
   
  @return The function returns code numbers in the parameter cn.

  Created by JK, 25.6.2001
*/
void mechtop::give_code_numbers (long eid,long *cn)
{
  Gtm->give_code_numbers (eid,cn);
}



/**
  Function extracts code numbers of actual node
   
  @param nid - node id
  @param cn - array containing code numbers on node (output)
   
  @return The function returns code numbers in the parameter cn.

  Created by JK, 7.12.2002
*/
void mechtop::give_node_code_numbers (long nid,long *cn)
{
  Gtm->give_node_code_numbers (nid,cn);
}



/**
  Function extracts code numbers of multipliers on actual node
   
  @param nid - node id
  @param lid - layer id
  @param cn - array containing all code numbers on node (output)
   
  @return The function returns code numbers of multipliers in the parameter cn.

  Created by JK, 7.12.2002
*/
void mechtop::give_mult_code_numbers (long nid,long lid,long *cn)
{
  Gtm->give_mult_code_numbers (nid,lid,cn);
}



/**
  Function returns node coordinates of the appropriate element.
   
  @param x - vector containing node coordinates (output)
  @param eid - element id
   
  @return The function returns x-coordinates in the parameter x. 

  Created by JK, 19.7.2001
*/
void mechtop::give_node_coord1d (vector &x,long eid)
{
  Gtm->give_node_coord1d (x,eid);
}



/**
  Function returns node coordinates of the appropriate element.
   
  @param x,y - vectors containing node coordinates (output)
  @param eid - element id
   
  @return The function returns 2D(x,y) coordinates in parameters x and y. 

  Created by JK, 19.7.2001
*/
void mechtop::give_node_coord2d (vector &x,vector &y,long eid)
{
  Gtm->give_node_coord2d (x,y,eid);
}



/**
  Function returns node coordinates of the appropriate element.
   
  @param x,z - vectors containing coordinates
  @param eid - element id
   
  @return The function returns 2D(x,z) coordinates in parameters x and z. 

  Created by JK, 19.7.2001
*/
void mechtop::give_node_coord2dxz (vector &x,vector &z,long eid)
{
  Gtm->give_node_coord2dxz (x,z,eid);
}



/**
  Function returns node coordinates of the appropriate element.
   
  @param x,y,z - vectors containing coordinates (output)
  @param eid - element id
   
  @return The function returns 3D(x,y,z) coordinates in parameters x, y and z. 

  Created by JK, 19.7.2001
*/
void mechtop::give_node_coord3d (vector &x,vector &y,vector &z,long eid)
{
  Gtm->give_node_coord3d (x,y,z,eid);
}



/**
  The function computes coordinates of all element integration points.
  
  @param eid - element id
  @param coord - %matrix with coordinates x,y,z of integration points stored in rows,
                 dimensions of matrix have to be (tnip, 3), 
                 coord[i][1] means y coordinate of the i-th integration point of the given element
  
  @return The function returns ip coordinates in the parameter coord.

  Created by Tomas Koudelka, 11.2011
*/
void mechtop::give_ipcoord_elem(long eid, matrix &coord)
{
  long ii, jj, ipp, i, k, l;
  long nb = give_nb(eid);
  vector aux(ASTCKVEC(3));

  k = 0;
  for(ii=0; ii<nb; ii++)
  {
    for(jj=0; jj<nb; jj++)
    {
      ipp=elements[eid].ipp[ii][jj];
      for (i=0; i<give_nip(eid, ii, jj); i++)
      {
        fillv(0.0, aux);
        ipcoord(eid, ipp+i, ii, jj, aux);
        for (l=0; l<3; l++)
          coord[k][l] = aux[l];
        k++;
      }  
      
    }
  }
}



/** 
  The function searches on the given element eid for the closest integration point 
  to the point of given coordinates [px,py,pz]. It retruns the closest integration point 
  id and its natural coordinates.

  @param eid - element id whose integration points are investigated (input)
  @param px, py, pz - global coordinates of the given point (input)
  @param ipm - array of mapping stuctures where the global coordinates of integration points are stored (input)
               ipm[ipp] = ipmap structure for the MEFEL integration point ipp (input)
  @param iptol - tolerance for the distance between [px,py,pz] and integration point below which the 
                 points are assumed to be identical (input)
  @param xi, eta, zeta - natural coordinates of the closest integration point (output)

  @retval ipp <  0 - if the distance to closest integration point is out of tolerance ipctol
  @retval ipp >= 0 - if the point [px,py,pz] is assumed to be identical with the integration point ipp

  Created by Tomas Koudelka, 2.12.2016
*/  
long mechtop::give_closest_ip_ncoord(long eid, double px, double py, double pz, 
                                     ipmap *ipm, double iptol, double &xi, double &eta, double &zeta)
{
  long i, nip = give_tnip(eid);
  long ipp = elements[eid].ipp[0][0];
  long ipc = ipp;
  double dmin=std::numeric_limits<double>::max();
  double d;
  vector ncoord(ASTCKVEC(3));

  for(i=0; i<nip; ipp++, i++)
  {
    d = sqr(px - ipm[ipp].x) + sqr(py - ipm[ipp].y) + sqr(pz - ipm[ipp].z);
    if (d < dmin)
    {
      dmin = d;
      ipc = ipp;
    }
  }
  // compute the minimum distnace 
  dmin = sqrt(dmin);
  // get the natural coordinates of the closest integration point
  ipncoord(eid, ipc, ncoord);
  xi   = ncoord[0];
  eta  = ncoord[1];
  zeta = ncoord[2];

  if (dmin <= iptol)
    return ipc;
  
  return -1;
}



/**
  Function returns length of the appropriate element.

  @param eid - element id

  @return The function returns element length.

  Created by JK, 19.7.2001
*/
double mechtop::give_length (long eid)
{
  long nn = give_nne(eid);
  double l;

  if (give_dimension(eid) > 1)
    print_err("invalid dimension on element is detected", __FILE__, __LINE__, __func__);

  vector x(ASTCKVEC(nn)), y(ASTCKVEC(nn)), z(ASTCKVEC(nn));
  Gtm->give_node_coord3d (x, y, z, eid);
  l = sqrt(sqr(x[1]-x[0]) + sqr(y[1]-y[0]) + sqr(z[1]-z[0]));
  return l;
}



/**
  Function returns area of the appropriate element.

  @param eid - element id

  @return The function returns area of 2D element.

  Created by JK, 19.7.2001
*/
double mechtop::give_area (long eid)
{
  long nn = give_nne(eid);
  double a = 0.0;
  if (nn < 2)
    print_err("invalid number of nodes on element is detected", __FILE__, __LINE__, __func__);

  vector x(ASTCKVEC(nn)), y(ASTCKVEC(nn)), z(ASTCKVEC(nn));
  Gtm->give_node_coord3d (x, y, z, eid);
  switch (Mt->give_elem_type(eid)){
  case planeelementlt:
  case planeelementqt:
  case planeelementrotlt:
  case axisymmlt:
  case axisymmqt:
    a = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
    break;
  case planeelementlq:
  case planeelementqq:
  case planeelementrotlq:
  case axisymmlq:
  case axisymmqq:
  case axisymmcq:
    a = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]) + (x[2]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[2]-y[0]) ) / 2.0;
    break;
  case cctel:
  case dktel:
  case dstel:
    a = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
    break;
  case q4plateel:
    a = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]) + (x[2]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[2]-y[0]) ) / 2.0;
    break;
  default:
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  return a;
}



/**
  Function returns volume of the appropriate element.

  @param eid - element id

  @return The function returns area of 3D element.

  Created by JK, 19.7.2001
*/
double mechtop::give_volume (long eid)
{
  long nn = give_nne(eid);
  double v = 0.0, det = 0.0;
  long nip = give_tnip(eid);
  vector gp(ASTCKVEC(nip)), w(ASTCKVEC(nip));
  if (nn < 2)
    print_err("invalid number of nodes on element is detected", __FILE__, __LINE__, __func__);

  vector x(ASTCKVEC(nn)), y(ASTCKVEC(nn)), z(ASTCKVEC(nn));
  Gtm->give_node_coord3d (x, y, z, eid);
  switch (give_elem_type(eid)){
  case lineartet:
  case quadrtet:
    jac_3d (v, x, y, z, 0.0, 0.0, 0.0);
    v/=6.0; // the volume of unit element
    break;
  case linearhex:
  case quadrhex:
  case lineartetrot:
  case linearhexrot:
    /*
    gauss_points (gp.a, w.a, 2);
    v = 0;
    jac_3d (vip, x, y, z, gp[0], gp[0], gp[0]);
    v += vip;
    jac_3d (vip, x, y, z, gp[0], gp[0], gp[1]);
    v += vip;
    jac_3d (vip, x, y, z, gp[0], gp[1], gp[0]);
    v += vip;
    jac_3d (vip, x, y, z, gp[0], gp[1], gp[1]);
    v += vip;
    jac_3d (vip, x, y, z, gp[1], gp[0], gp[0]);
    v += vip;
    jac_3d (vip, x, y, z, gp[1], gp[0], gp[1]);
    v += vip;
    jac_3d (vip, x, y, z, gp[1], gp[1], gp[0]);
    v += vip;
    jac_3d (vip, x, y, z, gp[1], gp[1], gp[1]);
    v += vip;*/    
    v = linhex_volume(x, y, z);
    break;
  default:
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
    
  if(v < 0.0){
    print_err("wrong numbering of nodes on 3D element number %ld, negative volume! v = %e",__FILE__,__LINE__,__func__,eid+1,v);
    det = det3d (x.a,y.a,z.a);
    print_err("kontrola determinantu det = %e",__FILE__,__LINE__,__func__,eid+1,det);
    abort();
  }

  return v;
}



/**
  Function returns maximum number of components of strain/stress array at nodes in the 
  parameter max_ncompstrn and maximum number of components of strain/stress array on 
  elements in the parameter max_ncompstre
   
  @param[out] max_ncompstrn - maximum number of components of strain/stress array at nodes
  @param[out] max_ncompstre - maximum number of components of strain/stress array on elements

  @return The function returns maximum number stress/strain components in the parameters max_ncompstrn 
          and max_ncompstre.
   
  Created by Tomas Koudelka,
*/
void mechtop::give_maxncompstr (long &max_ncompstrn, long &max_ncompstre)
{
  long i;
  max_ncompstrn = max_ncompstre = 0;

  for (i=0;i<ne;i++)
  {
    if (max_ncompstre < Mm->ip[elements[i].ipp[0][0]].ncompstr)
      max_ncompstre = Mm->ip[elements[i].ipp[0][0]].ncompstr;
  }      
  for (i=0;i<nn;i++)
  {
    if (max_ncompstrn < nodes[i].ncompstr)
      max_ncompstrn = nodes[i].ncompstr;
  }      
}



/**
  Function returns maximum number of components of other array at nodes in the  
  parameter max_nncompo and maximum number of components of eqother array on 
  elements in the parameter max_encompo
   
  @param max_encompo - maximum number of components of other array at nodes
  @param max_nncompo - maximum number of components of eqother array on elements

  @return The function returns maximum number other array components in the parameters max_nncompo 
          and max_encompo.
   
  Created by Tomas Koudelka,
*/
void mechtop::give_maxncompo (long &max_nncompo, long &max_encompo)
{
  long i;
  max_nncompo = max_encompo = 0;

  for (i=0;i<ne;i++)
  {
    if (max_encompo < Mm->ip[elements[i].ipp[0][0]].ncompeqother)
      max_encompo = Mm->ip[elements[i].ipp[0][0]].ncompeqother;
  }      
  for (i=0;i<nn;i++)
  {
    if (max_nncompo < nodes[i].ncompother)
      max_nncompo = nodes[i].ncompother;
  }      
}



/**
  The function returns the maximum number of dofs from all nodes in the problem.
  This value is either detected and stored in the max_ndofn or it returns the max_ndofn.
  It is used in the mechprint.cpp due to GiD output.

  @return The function returns max_ndofn value.
  
  Created by Tomas Koudelka 03.2013
   
*/
long mechtop::give_maxndofn()
{
  long i, ndofn;
  if (max_ndofn < 0)
  {
    for(i=0; i<nn; i++)
    {
      ndofn = give_ndofn(i);
      if (ndofn > max_ndofn)
        max_ndofn = ndofn;
    }
  }
  return max_ndofn;
}



/**
  Function returns numbers of local coordinate systems at nodes.
   
  @param nod - array containing node numbers

  @return The function returns id of local coordinate system in the parameter nod.

  Created by JK,
*/
long mechtop::locsystems (ivector &nod)
{
  long i,j,n;
  n=nod.n;
  j=0;
  for (i=0;i<n;i++){
    j+=nodes[nod[i]].transf;
  }
  return j;
}



/**
  Function alocates array containing reactions.

  @return The function does not return anything.

  Created by JK,
*/
void mechtop::comreacnod (void)
{
  long i,j,k,ndofn,nmn,nid;
  
  //  loop over the number of nodes
  for (i=0;i<nn;i++){
    
    //  the number of DOFs in the i-th node
    ndofn = Gtm->give_ndofn (i);
    
    if (ndofn>0){
      //  this is a classical node, not a hanging node
      
      //  loop over the number of DOFs in node
      for (j=0;j<ndofn;j++){
	if (Gtm->give_dof (i,j)<1){
	  nodes[i].react=1;
	  if (nodes[i].r==NULL)
	    nodes[i].r = new double [ndofn];
	  memset (nodes[i].r,0,ndofn*sizeof(double));
	  break;
	}
      }//  end of the loop over the number of DOFs in node
    }//  end of the statement if (ndofn>0){
    if (ndofn<0){
      //  this is a hanging node
      
      //  the number of master nodes
      nmn=0-ndofn;
      
      //  loop over the number of master nodes
      for (j=0;j<nmn;j++){
	//  master node id
	nid = Gtm->gnodes[i].mnodes[j];
	//  the number of DOFs on the master node
	ndofn = Gtm->give_ndofn (nid);
	
	//  loop over the number of DOFs in node
	for (k=0;k<ndofn;k++){
	  if (Gtm->give_dof (nid,k)<1){
	    nodes[i].react=1;
	    if (nodes[i].r==NULL)
	      nodes[i].r = new double [ndofn];
	    memset (nodes[i].r,0,ndofn*sizeof(double));
	    break;
	  }
	}//  end of the loop over the number of DOFs in node
	
      }//  end of the loop over the number of master nodes
      
    }//  end of the statement if (ndofn<0){
    
  }//  end of the loop over the number of nodes
}



/**
  Function marks elements necessary for reaction evaluation.

  @return The function does not return anything.

  Created by JK,
*/
void mechtop::comreacelem (void)
{
  long i,j,k,nne;
  ivector nod;
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      
      //  the number of nodes on element
      nne=give_nne (i);
      reallocv (RSTCKIVEC(nne,nod));
      //  nodes on element
      give_elemnodes (i,nod);
      
      //  loop over the number of nodes on element
      for (j=0;j<nne;j++){
	k=nod[j];
	if (nodes[k].react==1){
	  elements[i].react=1;
	  break;
	}
      }//  end of the loop over the number of nodes on element
    }
  }
}



/**
  Function marks elements necessary for reaction evaluation and
  alocates array containing reactions.

  @return The function does not return anything.

  Created by JK,
*/
void mechtop::comreac (void)
{
  comreacnod ();
  comreacelem ();
}



/**
  Function computes nodal values of strains with help of least square problem solution.
   
  @param nodes - array containing nodes of element
  @param nx,ny,nz - arrays containing natural coordinates of nodes
  @param lhs - array of coefficients of linear functions (solution of least square problem)
  @param dim - problem dimension
  @param fi - position of first component in nodal strain array
  @param ncomp - number of stored components
  @param lcid - load case id

  @return The function does not return anything.

  Created by JK, 10.5.2002
*/
void mechtop::strain_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
				   double *lhs,long dim,long fi,long ncomp,long lcid)
{
  long i,j;
  vector eps(ASTCKVEC(ncomp));
  
  if (dim==1){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	eps[j] = lhs[j*2]+lhs[j*2+1]*nx[i];
      }
      nodes[nod[i]].storestrain (lcid,fi,ncomp,eps);
    }
  }
  if (dim==2){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	eps[j] = lhs[j*3]+lhs[j*3+1]*nx[i]+lhs[j*3+2]*ny[i];
      }
      nodes[nod[i]].storestrain (lcid,fi,ncomp,eps);
    }
  }
  if (dim==3){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	eps[j] = lhs[j*4]+lhs[j*4+1]*nx[i]+lhs[j*4+2]*ny[i]+lhs[j*4+3]*nz[i];
      }
      nodes[nod[i]].storestrain (lcid,fi,ncomp,eps);
    }
  }
}



/**
  Function computes nodal values of stresses with help of least square problem solution.
   
  @param nodes - array containing nodes of element
  @param nx,ny,nz - arrays containing natural coordinates of nodes
  @param lhs - array of coefficients of linear functions (solution of least square problem)
  @param dim - problem dimension
  @param fi - position of first component in nodal strain array
  @param ncomp - number of stored components
  @param lcid - load case id

  @return The function does not return anything.

  Created by JK, 10.5.2002
*/
void mechtop::stress_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
				   double *lhs,long dim,long fi,long ncomp,long lcid)
{
  long i,j;
  vector sig(ASTCKVEC(ncomp));
  
  if (dim==2){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	//Mt->nodes[nodes[i]].stress[pos+j] = lhs[j*3]+lhs[j*3+1]*nx[i]+lhs[j*3+2]*ny[i];
	sig[j] = lhs[j*3]+lhs[j*3+1]*nx[i]+lhs[j*3+2]*ny[i];
      }
      nodes[nod[i]].storestress (lcid,fi,ncomp,sig);

    }
  }
  if (dim==3){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	//Mt->nodes[nodes[i]].stress[pos+j] = lhs[j*4]+lhs[j*4+1]*nx[i]+lhs[j*4+2]*ny[i]+lhs[j*4+3]*nz[i];
	sig[j] = lhs[j*4]+lhs[j*4+1]*nx[i]+lhs[j*4+2]*ny[i]+lhs[j*4+3]*nz[i];
      }
      nodes[nod[i]].storestress (lcid,fi,ncomp,sig);
    }
  }
}



/**
  Function computes nodal values of other quantities with help of least square problem solution.
   
  @param nodes - array containing nodes of element
  @param nx,ny,nz - arrays containing natural coordinates of nodes
  @param lhs - array of coefficients of linear functions (solution of least square problem)
  @param dim - problem dimension
  @param fi - position of first component in nodal strain array
  @param ncomp - number of stored components

  @return The function does not return anything.

  Created by JK, 10.5.2002
*/
void mechtop::other_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
				   double *lhs,long dim,long fi,long ncomp)
{
  long i,j;
  vector other(ASTCKVEC(ncomp));
  
  if (dim==1){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	other[j] = lhs[j*2]+lhs[j*2+1]*nx[i];
      }
      nodes[nod[i]].storeother (fi,ncomp,other);
    }
  }
  if (dim==2){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	other[j] = lhs[j*3]+lhs[j*3+1]*nx[i]+lhs[j*3+2]*ny[i];
      }
      nodes[nod[i]].storeother (fi,ncomp,other);

    }
  }
  if (dim==3){
    for (i=0;i<nod.n;i++){
      for (j=0;j<ncomp;j++){
	other[j] = lhs[j*4]+lhs[j*4+1]*nx[i]+lhs[j*4+2]*ny[i]+lhs[j*4+3]*nz[i];
      }
      nodes[nod[i]].storeother (fi,ncomp,other);
    }
  }
}



/**
  The function returns  required scalar qunatity in the given node.
  The quantity IS returned in the argument qv.

  @param[in] nid - nodal id
  @param[in] mqn - quantity name
  @param[in] reftensq - referenced tensor quantity, it is used only if the mqn 
                        is an tensor derived value (invariant, norm,...)
  @param[in] lcid - load case id or eigenvalue id.
  @param[out] qv - required quantity value.

  @return The function returns the quantity components in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechtop::give_nodal_quant(long nid, mechquant mqn, mechquant reftensq, long lcid, double &qv)
{
  matrix auxm;
  vector eps, epst, sig;
  qv = 0.0;
  switch (mqn){
    case time_q:          // time (constant scalar quantity at whole domain)
      qv = Mp->time;
      break;
    case step_id_q:       // step id (constant scalar quantity at whole domain)
      qv = Mp->istep;
      break;
    case load_fact_q:     // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
      qv = Mp->lambda;
      break;
    case eigval_q:        // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
      //qv = Lsrs->eigv[lcid];
      break;
    case first_inv:       // A11 + A22 + A33
      reallocm(RSTCKMAT(3, 3, auxm));
      give_nodal_quant(nid, reftensq, nomech_q, lcid, auxm);
      qv = first_invar(auxm);;
      break;
    case second_inv:      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2
      reallocm(RSTCKMAT(3, 3, auxm));
      give_nodal_quant(nid, reftensq, nomech_q, lcid, auxm);
      qv = second_invar(auxm);;
      break;
    case third_inv:       // det|A|
      reallocm(RSTCKMAT(3, 3, auxm));
      give_nodal_quant(nid, reftensq, nomech_q, lcid, auxm);
      qv = third_invar(auxm);;
      break;
    case tensor_norm:     // ||A|| : sqrt(a_ij*a_ij)
      reallocm(RSTCKMAT(3, 3, auxm));
      give_nodal_quant(nid, reftensq, nomech_q, lcid, auxm);
      qv = norm(auxm);;
      break;
    case strain_vol:      // eps_v : eps_x + eps_y + eps_z
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = first_invar(eps);
      break;
    case mean_stress:     // sig_m : (sig_x + sig_y + sig_z)/3
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = (1.0/3.0)*first_invar(sig);
      break;
    case j2inv:           // negative value of the second invariant of stress deviator, i.e. J2 : 1/2 s_ij s_ij
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = j2_stress_invar(sig);      
      break;
    case von_mises_stress:// sig_eff : sqrt(3*J2)                
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = 3.0*j2_stress_invar(sig);
      qv = sqrt(qv);
      break;
    //case cons_param:      // consistency parameter gamma
    //case damage_scal:     // scalar damage (omega)
    //case damaget_scal:    // scalar damage (omega_t) in tension
    //case damagec_scal:    // scalar damage (omega_c) in compression      
    case eps_x:           // strain components as scalars
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(0);
      break;
    case eps_y:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(1);
      break;
    case eps_z:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(2);
      break;
    case gamma_yz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(3);
      break;
    case gamma_xz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(4);
      break;
    case gamma_xy:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = eps(5);
      break;
    case eps_yz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = 0.5*eps(3);
      break;
    case eps_xz:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = 0.5*eps(4);
      break;
    case eps_xy:
      reallocv(RSTCKVEC(6, eps));
      give_full_vector(eps, nodes[nid].strain, guess_ssst(nodes[nid].ncompstr));
      qv = 0.5*eps(5);
      break;
    case sig_x:           // stress components as scalars
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(0);
      break;
    case sig_y:
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(1);
      break;
    case sig_z:
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(2);
      break;
    case tau_yz:
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(3);
      break;
    case tau_xz:
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(4);
      break;
    case tau_xy:
      reallocv(RSTCKVEC(6, sig));
      give_full_vector(sig, nodes[nid].stress, guess_ssst(nodes[nid].ncompstr));
      qv = sig(5);
      break;
    //case eps_px:          // components of plastic strains considered as scalars, will be implemented with the help of new materialpoint concept
    //case eps_py:
    //case eps_pz:
    //case eps_pyz:
    //case eps_pxz:
    //case eps_pxy:
    default:
      print_err("unknown scalar quantity type (%d) is required at node %ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), nid+1);
      abort();
  }
}



/**
  The function returns required components of a vector qunatity in the given node.
  The quantity components are returned in the argument qv.

  @param[in] nid - nodal id
  @param[in] mqn - quantity name
  @param[in] lcid - load case id or eigenvalue id.
  @param[out] qv - required quantity value.

  @return The function returns the quantity components in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechtop::give_nodal_quant(long nid, mechquant mqn, long lcid, sel &selcomp, gnodvalvm &nv, vector &qv)
{
  long i, id=0;
  gnode &nod = Gtm->gnodes[nid];
  
  switch (mqn){
    case displ_q:        // displacement vector
      for(i=0; i<nod.ndofn; i++){
        if(selcomp.presence_id(i)){
          qv(id) = nv.displv[nod.cn[i]];
          id++;
        }
      }
      break;
    case react_q:        // reactions vector at nodes     
      for(i=0; i<nod.ndofn; i++){
        if(selcomp.presence_id(i)){
          qv(id) = nodes[nid].r[i];
          id++;
        }
      }
      break;
    case load_vect_q:    // load vector at nodes
      for(i=0; i<nod.ndofn; i++){
        if(selcomp.presence_id(i)){
          qv(id) = nv.loadv[nod.cn[i]];
          id++;
        }
      }
      break;
    case int_force_q:    // internal force vector
      for(i=0; i<nod.ndofn; i++){
        if(selcomp.presence_id(i)){
          qv(id) = nv.iforv[nod.cn[i]];
          id++;
        }
      }
      break;
    case resid_vect_q:   // residual vector
      for(i=0; i<nod.ndofn; i++){
        if(selcomp.presence_id(i)){
          qv(id) = nv.residv[nod.cn[i]];
          id++;
        }
      }
      break;
    //case nonmech_q:       // nonmechanical quantities are defined on ip
    case other_q:        // other array values as a vector
      for(i=0; i < nodes[nid].ncompother; i++){
        if (selcomp.presence_id(i)){
          qv(id) = nodes[nid].other[i];
          id++;
        }
      }      
      break;
    default:
      print_err("unknown vector quantity type (%d) is required at node %ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), nid+1);
      abort();
  }
}



/**
  The function returns a required tensor qunatity in a given node.
  The quantity components are returned in the argument qv.

  @param[in] nid - nodal id
  @param[in] mqn - quantity name
  @param[in] reftensq - referenced quantity handled as an ordinary tensor
  @param[in] lcid - load case id or eigenvalue id.
  @param[out] qv - required quantity value.

  @return The function returns the quantity components in the argument qv.

  Created by Tomas Koudelka 09.2023
*/
void mechtop::give_nodal_quant(long nid, mechquant mqn, mechquant reftensq, long lcid, matrix &qv)
{
  vector eps, meps, sig, msig;
  double *epsptr, *sigptr;
  long ncompstr;
  double inv;
    
  switch (mqn){
    case strain_q:        // strain tensor
      ncompstr  = nodes[nid].ncompstr;
      epsptr = nodes[nid].strain+lcid*ncompstr;
      makerefv(eps, epsptr, ncompstr);
      vector_tensor(eps, qv, guess_ssst(nodes[nid].ncompstr), strain);
      break;
    case stress_q:        // stress tensor
      ncompstr  = nodes[nid].ncompstr;
      sigptr = nodes[nid].stress+lcid*ncompstr;
      makerefv(sig, sigptr, ncompstr);
      vector_tensor(sig, qv, guess_ssst(nodes[nid].ncompstr), stress);
      break;
    //case tempr_strain_q:  // temperature strain tensor
    //case eig_strain_q:    // eigenstrain tensor
    //case eig_stress_q:    // eigenstress tensor
    case macrostrain_q:   // macrostrain (constant tensor quantity at whole domain)
      reallocv(RSTCKVEC(Mt->max_ncompstr, meps));
      macrostrains(lcid, meps);
      vector_tensor(meps, qv, guess_ssst(nodes[nid].ncompstr), strain);
      break;
    case macrostress_q:   // macrostress (constant tensor quantity at whole domain)
      reallocv(RSTCKVEC(Mt->max_ncompstr, msig));
      macrostresses(lcid, msig);
      vector_tensor(msig, qv, guess_ssst(nodes[nid].ncompstr), stress);
      break;
    case tensdeviator:    // D_ij : A_ij - delta_ij*A_kk/3 
      give_nodal_quant(nid, reftensq, nomech_q, lcid, qv);      
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    case strain_deviator: // e_ij : eps_ij - delta_ij*eps_v/3
      makerefv(eps, nodes[nid].strain, nodes[nid].ncompstr);
      vector_tensor(eps, qv, guess_ssst(nodes[nid].ncompstr), strain);
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    case stress_deviator: // s_ij : sig_ij - delta_ij*sig_m
      makerefv(sig, nodes[nid].stress, nodes[nid].ncompstr);      
      vector_tensor(sig, qv, guess_ssst(nodes[nid].ncompstr), stress);
      inv = (qv(0,0)+qv(1,1)+qv(2,2))/3.0;
      qv(0,0) -= inv;
      qv(1,1) -= inv;
      qv(2,2) -= inv;
      break;
    //case strain_pl:       // plastic strain tensor
    //case damage_tens:     // damage tensor (Omega)
    //case damaget_tens:    // damage tensor (Omega_t) in tension
    //case damagec_tens:    // damage tensor (Omega_c) in compression
    default:
      print_err("unknown tensor quantity type (%d) is required at node %ld.\n",
                __FILE__, __LINE__, __func__, int(mqn), nid+1);
      abort();
  }
}



void mechtop::give_nodal_quant(long nid, mechquant mq, mechquant reftensq, long lcid, quantrep qr,
                               sel &selcomp, gnodvalvm &nv, double &sv, vector &vv, matrix &tv)
{
  switch(qr){
    case scalq:
      give_nodal_quant(nid, mq, reftensq, lcid, sv);
      break;
    case vectq:
      give_nodal_quant(nid, mq, lcid, selcomp, nv, vv);
      break;
    case tensq:
      give_nodal_quant(nid, mq, reftensq, lcid, tv);
      break;
    default:
      print_err("unknown qunatity character type (%d)\n", __FILE__, __LINE__, __func__, int(qr));
      abort();
  }
}



/**
  The function converts code numbers stored at nodes to code numbers
  stored on elements. It is used for mixed 2D mesh of plane elements,
  bars and springs.

  @return The function does not return anything.

  Created by Tomas Koudelka - koudelka@cml.fsv.cvut.cz
*/
void mechtop::store_code_num_elem(void)
{
  long i, j, nne, ndofe, ndofn;
  elemtype te;
  ivector nod;
  ivector cn;
  for(i=0; i < ne; i++)
  {
    te    = give_elem_type(i);
    nne   = give_nne(i);
    ndofe = give_ndofe(i);
    reallocv(RSTCKIVEC(nne, nod));
    give_elemnodes(i, nod);
    switch(te)
    {
      case beam2d:
        break;
      case spring_1:
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 1)
        {
          print_err("spring_1 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[0] = cn[0];
        }
        break;
      }
      case spring_2:
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 2)
        {
          print_err("spring_2 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[1] = cn[1];
        }
        break;
      }
      case spring_3:
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 3)
        {
          print_err("spring_3 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[2] = cn[2];
        }
        break;
      }
      case spring_4:
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 4)
        {
          print_err("spring_4 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[3] = cn[3];
        }
        break;
      }
      case spring_5:
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 5)
        {
          print_err("spring_5 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[4] = cn[4];
        }
        break;
      }
      case spring_6: 
      {
        ndofn = give_ndofn(nod[0]);
        if (ndofn < 6)
        {
          print_err("spring_6 is connected to node %ld with wrong number of dofs(%ld)", __FILE__, __LINE__, __func__, nod[0], ndofn);
          abort();
        }
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofn];
          reallocv(RSTCKIVEC(ndofn, cn));
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofn);
          Mt->give_node_code_numbers(nod[0], cn.a);
          Gtm->gelements[i].cn[5] = cn[5];
        }
        break;
      }
      case bar2d:
      case barq2d:
      case planeelementlt:
      case planeelementqt:
      case planeelementlq:
      case planeelementqq:
        if (Gtm->gelements[i].cne == 0)
        {
          Gtm->gelements[i].cne = 1;
          Gtm->gelements[i].cn = new long[ndofe];
          memset(Gtm->gelements[i].cn, 0, sizeof(*Gtm->gelements[i].cn)*ndofe);
          for (j=0; j<nne; j++)
          {
            ndofn = Gtm->give_ndofn(nod[j]);
            reallocv(RSTCKIVEC(ndofn, cn));
            Mt->give_node_code_numbers(nod[j], cn.a);
            Gtm->gelements[i].cn[2*j]   = cn[0];
            Gtm->gelements[i].cn[2*j+1] = cn[1];
          }
        }
        break;
      default:
        print_err("unknown type of element is required", __FILE__, __LINE__, __func__); 
    }
  }
}



/**
  Function allocates arrays defined on nodes.
  Arrays strain, stress, other, ncontr_strain, vol_strain, ncontr_stress, vol_stress, ncontr_other and vol_other are allocated on each node
  according to the settings in the problem description class.
   
  @return The function does not return anything.

  Created by JK, 28.11.2006
*/
void mechtop::alloc_nodes_arrays (void)
{
  long i,j,nne,ncomp, ncompo;
  ivector nod;
  strastrestate ssst;
  elemtype te;
  
  if (Mp->strainpos==2 || Mp->strainpos==3 || Mp->stresspos==2 || Mp->stresspos==3 || Mp->otherpos==2){
    for (i=0;i<ne;i++){
      if (Gtm->leso[i]==1){
	//  number of nodes on element
	nne = give_nne (i);
	reallocv (RSTCKIVEC(nne,nod));
	give_elemnodes (i,nod);
	//  number of strain/stress components
	ncomp=give_tncomp (i);
	//  number of components in array other
	ncompo=Mm->ip[elements[i].ipp[0][0]].ncompeqother;     
	
	//  element type
	te = give_elem_type (i);
	if (te != shelltrelem && te != shelltrmelem && te != shellqelem){
	  //  strain/stress state (beam, plane stress, plane strain, etc.
	  ssst=Mt->give_ssst (i,0);
	  if ((ssst == planestrain) || (ssst == planestress)){
	    ncomp=4;
	  }
	}
	
	for (j=0;j<nne;j++){
	  if (ncomp > nodes[nod[j]].ncompstr){
	    //  if different number of components are used, the maximum number has to be used
            if (Mp->strainpos > 1)
              nodes[nod[j]].alloc_strain(ncomp,Mb->nlc);
            if (Mp->stresspos > 1)
              nodes[nod[j]].alloc_stress(ncomp,Mb->nlc);
	  }
	  if (ncompo > nodes[nod[j]].ncompother){
            //  if different number of components are used, the maximum number has to be used
            if (Mp->otherpos > 1)
              nodes[nod[j]].alloc_other(ncompo);
	  }
 	}
      }
    }
  }
}



/**
  Function allocates strain arrays defined on nodes. Auxiliary arrays ncontr_strain 
  and vol_strain may be allocated with respect to the settings in the problem 
  description class.
  
   
  @return The function does not return anything.

  Created by TKo, 11.6.2018
*/
void mechtop::alloc_nodal_strains(void)
{
  long i,j,nne,ncomp;
  ivector nod;
  strastrestate ssst;
  elemtype te;
  
  for (i=0;i<ne;i++)
  {
    if (Gtm->leso[i]==1)
    {
      // number of nodes on element
      nne = give_nne(i);
      reallocv(RSTCKIVEC(nne,nod));
      give_elemnodes(i,nod);
      // number of strain/stress components
      ncomp=give_tncomp(i);
      // element type
      te = give_elem_type(i);
      if (te != shelltrelem && te != shelltrmelem && te != shellqelem)
      {
        // strain/stress state (beam, plane stress, plane strain, etc.
        ssst=Mt->give_ssst (i,0);
        if ((ssst == planestrain) || (ssst == planestress))
          ncomp=4;
      }            
	
      for (j=0;j<nne;j++)
      {
        if (ncomp > nodes[nod[j]].ncompstr)
        {
          // if different number of components are used, the maximum number has to be used
          nodes[nod[j]].alloc_strain(ncomp, Mb->nlc);
        }
      }
    }
  }
}



/**
  Function allocates end nodes.
  End nodes are allocated in problems with hemivariational inequalities.
   
  @return The function does not return anything.

  Created by JK, 1.3.2009
*/
void mechtop::alloc_enodes ()
{
  if (Gtm->nen>0){
    enodes = new endnodem [Gtm->nen];
  }
}



/**
  Functions allocates edges.
  Edges are allocated in problems with hemivariational inequalities.
   
  @return The function does not return anything.

  Created by JK, 1.3.2009
*/
void mechtop::alloc_edges ()
{
  if (Gtm->nged>0){
    edges = new edgem [Gtm->nged];
  }
}



/**
  Function allocates nodes and elements for preprocessor

  @param nn - number of nodes
  @param ne - number of elements

  @return The function does not return anything.

  Created by Tomas Koudelka, 13.3.2002
*/
void mechtop::alloc_prep (long nn,long ne,elemtype *el_type)
{
  long i;
  nodes = new node [nn];
  elements = new element [ne];
  for(i=0; i<ne; i++)
    elements[i].te = el_type[i];
}



/**
  Function allocates array for meaning of DOFs.
   
  @return The function does not return anything.

  Created by JK, 1.2.2005
*/
void mechtop::alloc_meaning ()
{
  long i;
  
  for (i=0;i<nn;i++){
    nodes[i].alloc_meaning (i);
  }
}



/**
  Function allocates all necessary arrays for problems with changing
  number of nodes, elements and DOFs.
   
  @return The function does not return anything.

  Created by JK, 7.11.2006
*/
void mechtop::alloc_growstr ()
{
  long i, ndofn;
  
  // array of attained nodal forces in the case of growing_mech_structure problems
  nodeforce = new double* [nn];
  memset(nodeforce, 0, sizeof(*nodeforce)*nn);

  for (i=0;i<nn;i++){
    nodes[i].alloc_growstr (i);
    ndofn = give_ndofn(i);
    nodeforce[i] = new double[ndofn];
    memset(nodeforce[i], 0, sizeof(*nodeforce[i])*ndofn);
  }
  for (i=0;i<ne;i++){
    elements[i].alloc_growstr (i);
  }
}



/**
  Function defines meaning of DOFs.
   
  @return The function does not return anything.

  Created by JK, 1.2.2005
*/
void mechtop::define_meaning ()
{
  long i;
  elemtype te;
  
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      te = give_elem_type (i);
      
      switch (te){
      case beam2d:{
	Beam2d->define_meaning (i);
	break;
      }
      case beam3d:{
	Beam3d->define_meaning (i);
	break;
      }
      case planeelementlq:{
	Pelq->define_meaning (i);
	break; 
      }
      case lineartet:{
	Ltet->define_meaning (i);
	break; 
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
}




/**
  Function generates code numbers of Lagrange multipliers.
   
  @param n - on input: number of DOFs without multipliers
             on output: number of DOFs with multipliers
	      
  @return The function returns number of DOFs with multipliers in the parameter n.

  Created by JK, 7.12.2002
*/
void mechtop::gencodnumlagrmult (long &n)
{
  long i,j,ndofn,nl,nmult,*cnl,*cnu;
  
  n++;
  
  for (i=0;i<nln;i++){
    nl=Gtm->lgnodes[i].nl;
    ndofn=Gtm->gnodes[Gtm->lgnodes[i].nodes[0]].ndofn;
    Gtm->lgnodes[i].ndofn=ndofn;
    Gtm->lgnodes[i].cn = new long* [nl+1];
    
    if (Mp->tlm==1){
      nmult=2;
      Gtm->lgnodes[i].nmult = nmult;
      cnl = new long [ndofn];
      cnu = new long [ndofn];
      
      Gtm->lgnodes[i].cn[0] = new long [nmult];
      Gtm->lgnodes[i].cn[0][0]=0;
      Gtm->lgnodes[i].cn[0][1]=0;
      
      for (j=1;j<nl;j++){
	give_node_code_numbers (Gtm->lgnodes[i].nodes[j-1],cnl);
	give_node_code_numbers (Gtm->lgnodes[i].nodes[j],cnu);
	
	Gtm->lgnodes[i].cn[j] = new long [nmult];
	if (cnl[0]>0 || cnl[2]>0 || cnu[0]>0 || cnu[2]>0){
	  Gtm->lgnodes[i].cn[j][0]=n;  n++;
	}
	if (cnl[1]>0 || cnu[1]>0){
	  Gtm->lgnodes[i].cn[j][1]=n;  n++;
	}
	
      }
      Gtm->lgnodes[i].cn[nl] = new long [nmult];
      Gtm->lgnodes[i].cn[nl][0]=0;
      Gtm->lgnodes[i].cn[nl][1]=0;
      
      delete [] cnl;  delete [] cnu;
    }
    
    if (Mp->tlm==2){
      nmult=4;
      Gtm->lgnodes[i].nmult = nmult;
      cnl = new long [ndofn];
      cnu = new long [ndofn];
      
      Gtm->lgnodes[i].cn[0] = new long [nmult];
      Gtm->lgnodes[i].cn[0][0]=0;
      Gtm->lgnodes[i].cn[0][1]=0;
      Gtm->lgnodes[i].cn[0][2]=0;
      Gtm->lgnodes[i].cn[0][3]=0;
      
      for (j=1;j<nl;j++){
	give_node_code_numbers (Gtm->lgnodes[i].nodes[j-1],cnl);
	give_node_code_numbers (Gtm->lgnodes[i].nodes[j],cnu);
	
	Gtm->lgnodes[i].cn[j] = new long [nmult];
	if (cnl[0]>0 || cnl[4]>0 || cnu[0]>0 || cnu[4]>0){
	  Gtm->lgnodes[i].cn[j][0]=n;  n++;
	}
	if (cnl[1]>0 || cnl[3]>0 || cnu[1]>0 || cnu[3]>0){
	  Gtm->lgnodes[i].cn[j][1]=n;  n++;
	}
	if (cnl[2]>0 || cnu[2]>0){
	  Gtm->lgnodes[i].cn[j][2]=n;  n++;
	}
	if (cnl[5]>0 || cnu[5]>0){
	  Gtm->lgnodes[i].cn[j][3]=n;  n++;
	}
	
      }
      Gtm->lgnodes[i].cn[nl] = new long [nmult];
      Gtm->lgnodes[i].cn[nl][0]=0;
      Gtm->lgnodes[i].cn[nl][1]=0;
      Gtm->lgnodes[i].cn[nl][2]=0;
      Gtm->lgnodes[i].cn[nl][3]=0;
      
      delete [] cnl;  delete [] cnu;
    }
  }
  n--;
}



/**
  The function searches adjacenet elements.
  
  @return The function does not return anything.

  Created by JK, 
*/
/*
void mechtop::adjacelem (void)
{
  long i, j, k, l, m, ned, nned, nedc, nnedc, net;
  long *enod, *et, *enodc;

  adjelem = new long *[ne];
  memset(adjelem, 0, sizeof(*adjelem)*ne);
  
  for (i = 0; i < ne; i++){
    if (Gtm->leso[i]==1){
      ned  = give_ned(i);
      if (ned < 3)
        print_warning("element number %ld has number of edges < 3", __FILE__, __LINE__ ,__func__, i+1);

      Gtm->gelements[i].edgn = new long [ned];
      adjelem[i] = new long[ned];
      for (j = 0; j < ned; j++)
	{
	  adjelem[i][j] = -1;
	  Gtm->gelements[i].edgn[j] = -1;
	}
    }
  }
  tned = 0;
  for (i = 0; i < ne; i++){
    if (Gtm->leso[i]==1){
      ned  = give_ned(i);
      for (j = 0; j < ned; j++)
	{
	  if (Gtm->gelements[i].edgn[j] > -1)
	    continue;
	  nned = give_nned(i);
	  enod = new long [nned];
	  memset(enod, 0, sizeof(*enod)*nned);
	  give_edge_nodes (i, j, enod);
	  net = 0;
	  for (k = 0; k < nned; k++)
	    net += Gtm->nadjelnod[enod[k]];
	  et = new long [net];
	  memset(et, 0, sizeof(*et)*net);
	  l = 0;
	  for (k = 0; k < nned; k++)
	    {
	      for (m = 0; m < Gtm->nadjelnod[enod[k]]; m++)
		{
		  et[l] = Gtm->adjelnod[enod[k]][m];
		  l++;
		}
	    }
	  for (k = 0; k < net; k++)
	    {
	      if (et[k] == i)
		continue;
	      nedc = give_ned(et[k]);
	      for (m = 0; m < nedc; m++)
		{
		  if (Gtm->gelements[et[k]].edgn[m] > -1)
		    continue;
		  nnedc = give_nned(et[k]);
		  enodc = new long [nedc];
		  memset(enodc, 0, sizeof(*enodc)*nnedc);
		  give_edge_nodes (et[k], m, enodc);
		  if (compare_edges(enod, enodc, nned))
		    {
		      adjelem[i][j]     = et[k];
		      adjelem[et[k]][m] = i;
		      Gtm->gelements[i].edgn[j]     = tned;
		      Gtm->gelements[et[k]].edgn[m] = tned;
		      delete [] enodc;
		      break;
		    }
		  delete [] enodc;
		}
	      if (adjelem[i][j] > -1)
		break;
	    }
	  delete [] et;
	  delete [] enod;
	  Gtm->gelements[i].edgn[j] = tned;
	  tned++;
	}
    }
  }
}
*/



/**
  The function compares nodal numbers of two edges.

  @param enod - array with node numbers of the first compared edge
  @param enodc - array with node numbers of the second compared edge
  @param nned - number of components of enod and enodc (i.e. number of edge nodes)

  @retval 0 - edges are different
  @retval 1 - edges are identical

  Created by JK,
*/
long mechtop::compare_edges (long *enod, long *enodc, long nned)
{
  long i, j, nsn = 0;

  for (i = 0; i < nned; i++)
  {
    for (j = 0; j < nned; j++)
    {
      if (enod[i] == enodc[j])
      {
        nsn++;
        break;
      }
    }
  }
  if (nsn == nned)
    return 1;

  return 0;
}



/*
void mechtop::eadjacelem (void)
{
  long i, j, ned, ided;

  eadjelem = new long *[tned];
  memset (eadjelem, 0, sizeof(*eadjelem)*tned);
  for (i = 0; i < tned; i++)
  {
    eadjelem[i] = new long[2];
    eadjelem[i][0] = -1;
    eadjelem[i][1] = -1;
  }
  
  for (i = 0; i < ne; i++){
    if (Gtm->leso[i]==1){
      ned  = give_ned(i);
      for (j = 0; j < ned; j++)
	{
	  ided = Gtm->gelements[i].edgn[j];
	  if (eadjelem[ided][0] < 0)
	    {
	      eadjelem[ided][0] = i;
	      continue;
	    }
	  if (eadjelem[ided][1] < 0)
	    {
	      eadjelem[ided][1] = i;
	      continue;
	    }
          print_err("one edge connects more than 2 elements,\n"
                    " the function can be applied only on 2D problems", __FILE__, __LINE__, __func__);
	}
    }
  }
}
*/



/**
  Function returns nodal displacements.
  Function is intended for molecular computations.
   
  Created by JK, 19.6.2005
*/
/*
void mechtop::nodedisplacements ()
{
  long i,j,ndofn,lcid;
  double *r;
  
  //  the same number of unknowns at nodes is assumed at this version
  //  in case of different numbers of unknowns, function must be modified
  ndofn = give_ndofn (0);
  
  //  only one loading case is assumed at this version
  //  it must be generalized for more loading cases
  lcid=0;
  
  //  possible allocation
  if (nodedispl==NULL){
    nodedispl = new double* [nn];
    for (i=0;i<nn;i++){
      ndofn = give_ndofn (i);
      nodedispl[i] = new double [ndofn];
      for (j=0;j<ndofn;j++){
	nodedispl[i][j]=0.0;
      }
    }
  }
  
  //  loop over nodes
  for (i=0;i<nn;i++){
    ndofn = give_ndofn (i);
    r = new double [ndofn];
    //  nodal displacements
    noddispl (lcid,r,i);
    for (j=0;j<ndofn;j++){
      nodedispl[i][j]=r[j];
    }
    delete [] r;
  }
}
*/



/**
  Function establishes vector of initial displacements of elements
  on the new part of structure. It is intended to be used for 
  the growing structures.

  @param lcid - load case id
  @param ifn  - array of interface nodes indicators
   
  @return The function does not return anything.
 
  Created by TKo, 7.6.2013
*/
void mechtop::save_elem_inidispl(long lcid, long *ifn)
{
  long i, j, k, l, nne, ndofe, ndofn;
  ivector enod;
  vector r;
  
  //  load case id must be equal to zero
  for (i=0;i<ne;i++)
  {
    if (Gtm->leso[i] == 1) 
    // leso == 1 => new elements are assumed to be in stress free state 
    // after joining to old deformed structure
    {
      ndofe = give_ndofe (i);
      
      reallocv(RSTCKVEC(ndofe, r));

      // nodal displacements
      eldispl(lcid, i, r.a);

      // add displacements on interface nodes
      // that have got prescribed zero code numbers but 
      // there are nodal values caused by the attained displacemnets from 
      // the old part of structure connected to the interface nodes
      nne = give_nne(i);
      reallocv(RSTCKIVEC(nne, enod));
      give_elemnodes(i, enod);
      l = 0;
      for (j=0; j<nne; j++)
      {
        ndofn = give_ndofn(enod[j]);
        if (ifn[enod[j]])
        {
          for (k=0; k<ndofn; k++)
          {
            r[l] = nodes[enod[j]].nodval[k];
            l++;
          }
        }
        else
          l += ndofn;
      }

      // definition of initial displacements
      elements[i].initdisplacement(r.a, ndofe);
    }
  }
}



/**
  The function stores attained nodal displacements at nodes that were 
  free in previous time and they have prescribed displacement in the actual time.
  Additionally, there are cleaned displacements stored at nodes (nodval) for
  those nodes that were switched off in the acual time step.
  The function is intended to be used for the growing structures and
  it MUST be called BEFORE regeneration of code numbers.

  @param lcid - load case id
  @param time - actual time
  @param prev_time - previous time

  @return The function does not return anything but it changes nodedispl array prospectively.

  Created by Tomas Koudelka, 5.2013
*/
void mechtop::save_node_inidispl(long lcid, double time, double prev_time)
{
  long i,j,k,l,m,ndofn;
  vector r;

  for (i=0; i<nn; i++)
  {
    if (Gtm->lnso[i]) // for active nodes check changed supports
    {
      if (nodedispl[i]) // the node has got prescribed displacements for some dofs by time functions
      {
        ndofn = give_ndofn(i);
        reallocv(RSTCKVEC(ndofn, r));
        noddispl(lcid, r.a, i);
        for (j=0;j<ndofn;j++)
        {
          k=Gtm->gnodes[i].tgf[j];
          
          if (k < 0)  // dof is controlled by time functions of adjacent elements
            continue;

          l=Gtm->gf[k].getval_long (time);      // dof setup for actual time 
          m=Gtm->gf[k].getval_long (prev_time); // dof setup for previous time 

          if ((l < 1) && (m > 0)) // dof was free but now there is precribed displacement or support
          {                       
            // store attained displacements as initial ones
            nodedispl[i][j] = r[j];
          }
          if ((l > 0) && (m < 1) && nodedispl[i]) 
          {
            // dof had the prescribed displacement or support and now is free and there were 
            // some initial displacements stored => clear possible initial displacement
            // this step is not necessary but it can help in debugging
            nodedispl[i][j] = 0.0;
          }
          if ((l > 0) && (m > 0) && nodedispl[i]) // dof was free and now is also free
          {
            // this should capture cases where solver went back to the previous time step
            // and in such case some possible nonzero values should be cleaned
            nodedispl[i][j] = 0.0;
          }
        }
      }
    }
    else // for switched off nodes zero displacements array stored at node
    {
      ndofn = give_ndofn(i);
      memset(nodes[i].nodval, 0, sizeof(*nodes[i].nodval)*ndofn);
    }
  }
}



/**
  The function sets stress arrays to zeor on all integration points of new elements.
  It is intended to be used for the growing structures.

  @return The function does not return anything but it changes
          stress array on all integration points.

  Created by Tomas Koudelka, 7.6.2013
*/
void mechtop::clean_ip_new_elem()
{
  long i, j, ipp, tnip;

  for (i=0; i<ne; i++)
  {
    if (Gtm->leso[i] && (Gtm->gelements[i].auxinf == 0)) // for new elements
    {
      ipp = give_nip(i, 0, 0);
      tnip = give_tnip(i);
      for (j=0; j<tnip; j++)
      {
        Mm->ip[ipp].clean_strains(0);
        Mm->ip[ipp].clean_stresses(0);
        // eqotehr values MUST be preserved due to initial values calculated in Mm->initmaterialmodels
        ipp++;
      }
    }
  }
}



/**
  The function returns the first component of nodal displacements for the 
  selected nodes given by parameter nodes.

  @param nodes - %vector of selected nodes
  @param u - %vector of displacements (output)

  @return The function returns the first displacement component of selected nodes in 
          the parameter u.

  Created by JK, 23.6.2005
*/
void mechtop::give_noddispl_1d (ivector &nodes,vector &u)
{
  long i;
  
  for (i=0;i<nodes.n;i++){
    u[i]=nodedispl[nodes[i]][0];
  }
  
}



/**
  The function returns the first two components of nodal displacements for the 
  selected nodes given by parameter nodes.

  @param nodes - %vector of selected nodes
  @param u - %vector of the first displacement component (output)
  @param v - %vector of the second displacement component (output)

  @return The function returns displacement components at selected nodes in 
          parameters u and v.

  Created by JK, 23.6.2005
*/
void mechtop::give_noddispl_2d (ivector &nodes,vector &u,vector &v)
{
  long i;
  
  for (i=0;i<nodes.n;i++){
    u[i]=nodedispl[nodes[i]][0];
    v[i]=nodedispl[nodes[i]][1];
  }
}



/**
  The function returns first three components of nodal displacements for the 
  selected nodes given by parameter nodes.

  @param nodes - %vector of selected nodes
  @param u - %vector of the first displacement component (output)
  @param v - %vector of the second displacement component (output)
  @param w - %vector of the third displacement component (output)

  @return The function returns displacement components at selected nodes in 
          parameters u, v and w.

  Created by JK, 23.6.2005
*/
void mechtop::give_noddispl_3d (ivector &nodes,vector &u,vector &v,vector &w)
{
  long i;
  
  for (i=0;i<nodes.n;i++){
    u[i]=nodedispl[nodes[i]][0];
    v[i]=nodedispl[nodes[i]][1];
    w[i]=nodedispl[nodes[i]][2];
  }
  
}


/**
  Function initializes mechtop class from the initialized siftop structure.
  The function is used in the preprocessors esspecially.

  @param top - initialized siftop structure

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void mechtop::init_from_siftop (siftop *top)
{
  long i,j;

  //  node data
  nn = top->nn;
  //  element data
  ne = top->ne;
  // export and allocation of general topology
  top->exporttop(Gtm);

  if (Mespr==1)  fprintf (stdout,"\n number of nodes  %ld",nn);
  if (Mespr==1)  fprintf (stdout,"\n number of elements  %ld",ne);
  nodes = new node [nn];
  elements = new element [ne];

  for (i=0; i < ne; i++){
    if (Gtm->leso[i]==1){
      for(j=0; j<Gtm->gelements[i].give_nne(); j++)
	Gtm->gelements[i].nodes[j]--;
    }
  }
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      elements[i].tm = new mattype[1];
      elements[i].idm = new long[1];
      elements[i].idm[0] = top->elements[i].prop-1;
      switch (top->elements[i].type)
	{
	case ISOLinear1D:
	  elements[i].te = bar2d;
	  if (Bar2d == NULL)
	    Bar2d = new barel2d;
	  break;
	case ISOQuadratic1D:
	  elements[i].te = barq2d;
	  if (Barq2d == NULL)
	    Barq2d = new barelq2d;
	  break;
	case TriangleLinear:
	  elements[i].te = planeelementlt;
	  if (Pelt == NULL)
	    Pelt = new planeelemlt;
	  break;
	case TriangleQuadratic:
	  elements[i].te =planeelementqt;
	  if (Peqt == NULL)
	    Peqt = new planeelemqt;
	  break;
	case ISOLinear2D:
	  elements[i].te = planeelementlq;
	  if (Pelq == NULL)
	    Pelq = new planeelemlq;
	  break;
	case ISOQuadratic2D:
	  elements[i].te = planeelementqq;
	  if (Peqq == NULL)
	    Peqq = new planeelemqq;
	  break;
	case TetrahedronLinear:
	  elements[i].te = lineartet;
	  if (Ltet == NULL)
	    Ltet = new lintet;
	  break;
	case TetrahedronQuadratic:
	  elements[i].te = quadrtet;
	  if (Qtet == NULL)
	    Qtet = new quadtet;
	  break;
	case WedgeLinear:
	  elements[i].te = linearwed;
	  if (Lwed == NULL)
	    Lwed = new linwedge;
	  break;
	case WedgeQuadratic:
	  elements[i].te = quadrwed;
	  if (Qwed == NULL)
	    Qwed = new quadwedge;
	  break;
	case ISOLinear3D:
	  elements[i].te = linearhex;
	  if (Lhex == NULL)
	    Lhex = new linhex;
	  break;
	case ISOQuadratic3D:
	  elements[i].te = quadrhex;
	  if (Qhex == NULL)
	    Qhex = new quadhex;
	  break;
	default:
	  print_err("unknown type of element is used", __FILE__, __LINE__, __func__);
	  abort();
	}
    }
  }
}



/**
  Function cleans strains, stresses and other quantities defined at nodes.
  The function is used especially in stochastic or fuzzy computations.
   
  @return The function does not return anything.
 
  Created by JK, 23.8.2005
*/
void mechtop::clean_nodes ()
{
  long i;
  
  for (i=0;i<nn;i++){
    nodes[i].clean (Mb->nlc);
  }
}



/**
  Function cleans arrays of reactions defined at nodes.
  The function should be called befor reaction calculation.
   
  @return The function does not return anything.
 
  Created by TKo, 28.4.2016
*/
void mechtop::null_react()
{
  long i;
  
  for (i=0; i<nn; i++){
    if (nodes[i].react == 1)
      nullv(nodes[i].r, Mt->give_ndofn(i));
  }
}



/**
  Function computes norm of reactions at nodes.
  The function should be called after reaction calculation.
   
  @return The function returns norm of global reaction vector.
 
  Created by TKo, 14.2.2019
*/
double mechtop::compute_react_norm()
{
  long i, j, ndofn;
  double norm = 0.0;
  
  for (i=0L; i<nn; i++){
    if (nodes[i].react == 1){
      ndofn = Mt->give_ndofn(i);
      for (j=0L; j<ndofn; j++)
        norm += sqr(nodes[i].r[j]);
    }
  }
  norm = sqrt(norm);
  return norm;
}



/**
  Function establishes vector of initial displacements.

  @param lcid - load case id
   
  @return The function does not return anything.
 
  Created by JK, 28.2.2006
*/
void mechtop::initial_displ (long lcid)
{
  long i,j,k,ndofe;
  vector r;
  
  //  load case id must be equal to zero
  
  for (i=0;i<ne;i++){
    
    j=Gtm->leso[i];
    k=Gtm->gelements[i].auxinf;
    
    if (j==1 && k==0){
      //  added element
      //  number of DOFs on element
      ndofe = give_ndofe (i);
      
      reallocv(RSTCKVEC(ndofe, r));
      
      //  nodal displacements
      eldispl (lcid,i,r.a);
      //  definition of initial displacements
      elements[i].initdisplacement (r.a,ndofe);
      
      Gtm->gelements[i].auxinf=1;
    }
  }
}



/**
  The function saves backup of attained nodal values and it is used in problems with 
  growing number of nodes and elements.
   
  @return The function does not return anything but stores attained nodal displacements to 
          nodval arrays at nodes.
 
  Created by Tomas Koudelka, 7.6.2013
*/
void mechtop::save_nodval(long lcid)
{
  long i;
  
  for (i=0;i<nn;i++){
    if (Gtm->lnso[i])
    {
      // nodal values must be saved also for hanging nodes due to 
      // correct handling of elements with hanging nodes connected to interface 
      // between old and new part of structure
      noddispl(lcid, nodes[i].nodval, i);
    }
  }
}



/**
  Function restores nodal unknowns from nodes to the array lhs.
  The function is used in problems with growing number of nodes and elements.

  @param lhs - array of left hand side i.e. array of unknowns (output)
  @param n - total number of DOFs (dimension of array lhs)
   
  @return The function restores saved unknowns at nodes to the %vector given by 
          parameter lhs.

  Created by JK, 2.6.2006
  Modified by TKo, 17.6.2013
*/
void mechtop::restore_nodval(double *lhs,long n)
{
  long i,j,ndofn;
  ivector cn;
  
  nullv(lhs, n);

  for (i=0;i<nn;i++){
    ndofn = Gtm->give_ndofn (i);
    if (ndofn < 0) // hanging nodes do not need to restore attained nodal values
                   // because their dofs are not contained in the left hand side
      continue;
    reallocv(RSTCKIVEC(ndofn, cn));
    give_node_code_numbers (i,cn.a);
    
    for (j=0;j<ndofn;j++){
      if (cn[j]>0){
        lhs[cn[j]-1]=nodes[i].nodval[j];
      }
    }
  }
}



/**
  The function saves backup of attained nodal forces (i.e. components of load %vector).
  It can be used for the saving of several load vectors that are distinguish by the given 
  index id. It is used in problems with growing number of nodes and elements.
   
  @param f - array of load %vector components
  @return The function does not return anything but stores attained load vectors to 
          nodeforce arrays at nodes.
 
  Created by TKo, 3.11.2014
*/
void mechtop::save_nodforce(double *f)
{
  long i, ndofn;
  vector nf;
  
  for (i=0;i<nn;i++)
  {
    if (Gtm->lnso[i])
    {
      ndofn = Gtm->give_ndofn (i);
      if (ndofn < 0) // hanging nodes do not need to store attained nodal forces
                     // because they have been already included in the right hand side components
        continue;
      nf.n = ndofn;
      nf.a = nodeforce[i];
      select_nodforce(f, i, nf);
    }
  }
  nf.n = 0L;
  nf.a = NULL;
}



/**
  Function restores load vector from nodes to the array rhs.
  The function is used in problems with growing number of nodes and elements.

  @param rhs - array of right hand side i.e. load %vector (output)
  @param n - total number of DOFs (dimension of array rhs)
   
  @return The function restores saved unknowns at nodes to the %vector given by 
          parameter lhs.

  Created by TKo, 3.11.2014
*/
void mechtop::restore_nodforce(double *rhs,long n)
{
  long i,j,ndofn;
  ivector cn;
  
  nullv(rhs, n);

  for (i=0;i<nn;i++){
    ndofn = Gtm->give_ndofn (i);
    if (ndofn < 0) // hanging nodes do not need to restore attained nodal values
                   // because their dofs are not contained in the left hand side
      continue;
    reallocv(RSTCKIVEC(ndofn, cn));
    give_node_code_numbers (i,cn.a);
    
    for (j=0;j<ndofn;j++){
      if (cn[j]>0){
        rhs[cn[j]-1]=nodeforce[i][j];
      }
    }
  }
}



/**
  The function checks correct assignment of time functions of particular nodes 
  comparing to assignment of element time functions. It is used for growing structures especially.
  
  @return The function does not return anything.

  Created by TKr,
*/
long mechtop::mesh_check(void)
{
  long i,j,err=0,kk,jj;

  //elems and nodes check
  for(i=0;i<Mt->ne;i++){
    for(j=0;j<Gtm->gelements[i].nne;j++){
      
      jj = Gtm->gelements[i].nodes[j];
      for(kk = 0; kk < Gtm->gnodes[jj].ndofn;kk++){
	
	if((Gtm->gelements[i].tgf) < (Gtm->gnodes[jj].tgf[kk])){
	  print_err("different element node(%ld) time function number=%ld than element(%ld) time function number=%ld\n"
                    " (node number=%ld, element number=%ld)",
                    __FILE__, __LINE__, __func__, j+1,Gtm->gnodes[jj].tgf[kk]+1,i+1,Gtm->gelements[i].tgf+1,Gtm->gelements[i].nodes[j]+1,i+1);
	  err = 1;
	  return err;
	}
      }
    }
  }
  return err;
}



/**
  The function returns volume of the whole domain described by topology stored in the given 
  instance of mechtop.

  @return The function returns volume of all 3D or 2D elements in the topology.

  Created by Tomas Koudelka, 20.3.2014
  Modified by TKr 05/092018
*/
double mechtop::give_domain_vol()
{
  double ret = 0.0;
  long i;
  strastrestate ssst;
 
  /*
    FILE *out;
    
    
    out = fopen("elements.in", "wt");
    
    //numbering check only for debug:
    for (i=0; i < ne; i++){
    ret = give_volume(i);
    if(ret <= 0.0){
    fprintf(out,"%ld  100   %ld  %ld  %ld  %ld  0    20          1     1 1    1\n",i+1,Gtm->gelements[i].nodes[0]+1,Gtm->gelements[i].nodes[1]+1,Gtm->gelements[i].nodes[3]+1,Gtm->gelements[i].nodes[2]+1);
    }
    else
    fprintf(out,"%ld  100   %ld  %ld  %ld  %ld  0    20          1     1 1    1\n",i+1,Gtm->gelements[i].nodes[0]+1,Gtm->gelements[i].nodes[1]+1,Gtm->gelements[i].nodes[2]+1,Gtm->gelements[i].nodes[3]+1);
    }
    
    fflush(out);
    fclose(out);
  */

  for (i=0; i < ne; i++){
    ssst=Mt->give_ssst (i,0);
    if(ssst == spacestress)//for 3D elements
      {
	ret += give_volume(i);
      }
    else
      {
	ret += give_area(i);
      }
  }
  
  return ret;
}



/**
   function computes centers and radius of balls circumscribed to elements
   this function is used in connection with determination of hanging nodes

   JK, 23. 8. 2016
*/
void mechtop::circumscribed_balls ()
{
  long i;
  
  //  allocation of arrays
  Gtm->allocate_arrays_circumscribed_ball ();

  //  loop over the number of all elements
  for (i=0;i<ne;i++){
    if (elements[i].te==linearhex){
      Gtm->element_center (i);
    }
  }
  //  loop over the number of all elements
  for (i=0;i<ne;i++){
    if (elements[i].te==linearhex){
      Gtm->radius_circumscribed_ball (i);
    }
  }
}

void mechtop::inter ()
{
  long i,j,k,l,m,eind,del;
  double nors,dist;
  long numinter=0,numinterbef,numinterafter;
  
  //  vectors of nodal coordinates of bar elements
  vector xb(ASTCKVEC(2)),yb(ASTCKVEC(2)),zb(ASTCKVEC(2));
  //  vectors of nodal coordinates of hexahedral elements
  vector xh(ASTCKVEC(8)),yh(ASTCKVEC(8)),zh(ASTCKVEC(8));
  //  auxiliary vectors
  vector s(ASTCKVEC(3)),p(ASTCKVEC(3));
  
  vector xi(ASTCKVEC(40)), yi(ASTCKVEC(40)), zi(ASTCKVEC(40));
  vector xxi(ASTCKVEC(40)), yyi(ASTCKVEC(40)), zzi(ASTCKVEC(40));
  imatrix masnod(ASTCKIMAT(40,9));
  
  for (i=0;i<40;i++){
    for (j=0;j<9;j++){
      masnod[i][j]=-1;
    }
  }
  

  //  loop over all elements in the problem
  for (i=0;i<ne;i++){
    eind = Mp->barlist.presence_id (i);
    if (eind==1){
      //  the i-th element is a bar element which will be represented via hanging nodes
      give_node_coord3d (xb,yb,zb,i);
      
      //  direction vector
      s[0]=xb[1]-xb[0];
      s[1]=yb[1]-yb[0];
      s[2]=zb[1]-zb[0];
      //  norm of the direction vector
      nors = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
      if (nors<Mp->zero){
        print_err("there is zero direction vector on bar element n. %ld",__FILE__, __LINE__, __func__, i+1);
        abort ();
      }
      //  unit direction vector
      s[0]/=nors;
      s[1]/=nors;
      s[2]/=nors;
      
      //  loop over all elements
      for (k=0;k<ne;k++){
        if (give_elem_type (k) == linearhex){
          //  only hexahedral elements with linear approximation functions are taken into account
          
          //  distance of the center of the k-th element from the line defined by the bar element
          p[0]=Gtm->xc[k]-xb[0];
          p[1]=Gtm->yc[k]-yb[0];
          p[2]=Gtm->zc[k]-zb[0];
          
          dist = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2] - (s[0]*p[0]+s[1]*p[1]+s[2]*p[2]) * (s[0]*p[0]+s[1]*p[1]+s[2]*p[2]));
          
          if (dist<Gtm->circumrad[k]){
            //  the center of the k-th element is closer to the bar line than the radius of circumscribed ball
            //  the k-th element must be checked
            
            //  nodal coordinates of hexahedral element
            give_node_coord3d (xh,yh,zh,k);
	    
	    //  the number of intersections before analysis of the actual element
            numinterbef=numinter;
	    
            //  intersection of bar and hexahedron
            Gtm->bar_linhex_intersection (k,numinter,xb,yb,zb,xh,yh,zh,1.0e-10,xi,yi,zi,xxi,yyi,zzi,masnod);
	    
	    //  the number of intersections after analysis of actual element
	    numinterafter=numinter;
	    
	    if (numinterafter-numinterbef == 2){
	      fprintf (Out,"\n %ld   %le %le %le",nn+numinter-1,xi[numinter-2],yi[numinter-2],zi[numinter-2]);
	      fprintf (Out,"   %ld",0-masnod[numinter-2][0]);
	      for (j=0;j<masnod[numinter-2][0];j++){
		fprintf (Out," %ld",masnod[numinter-2][j+1]);
	      }
	      fprintf (Out,"   %le %le %le",xxi[numinter-2],yyi[numinter-2],zzi[numinter-2]);

	      fprintf (Out,"\n %ld   %le %le %le",nn+numinter,xi[numinter-1],yi[numinter-1],zi[numinter-1]);
	      fprintf (Out,"   %ld",0-masnod[numinter-1][0]);
	      for (j=0;j<masnod[numinter-1][0];j++){
		fprintf (Out," %ld",masnod[numinter-1][j+1]);
	      }
	      fprintf (Out,"   %le %le %le",xxi[numinter-1],yyi[numinter-1],zzi[numinter-1]);
	    }
	    if (numinterafter-numinterbef > 2){
	      //  some intersection point is has been found more times
	      

	      for (j=numinterbef;j<numinterafter;j++){
		if (masnod[j][0]==1){
		  for (l=numinterbef;l<numinterafter;l++){
		    if (l==j)  continue;
		    del=0;
		    for (m=0;m<masnod[l][0];m++){
		      if (masnod[l][m+1]==masnod[j][1]){
			del=1;
			break;
		      }
		    }
		    if (del==1){
		      for (m=0;m<masnod[l][0];m++){
			masnod[l][m+1]=-1;
		      }
		      masnod[l][0]=-1;
		    }
		  }
		}//  end of the statement if (masnod[j][0]==1){
		
		if (masnod[j][0]==2){
		  for (l=numinterbef;l<numinterafter;l++){
		    if (l==j)  continue;
		    del=0;
		    for (m=0;m<masnod[l][0];m++){
		      if (masnod[l][m+1]==masnod[j][1] || masnod[l][m+1]==masnod[j][2]){
			del++;
		      }
		    }
		    if (del==2){
		      for (m=0;m<masnod[l][0];m++){
			masnod[l][m+1]=-1;
		      }
		      masnod[l][0]=-1;
		    }
		  }
		}//  end of the statement if (masnod[j][0]==2){
	      }//  end of the loop for (j=numinterbef;j<numinterafter;j++){


	      del=numinterbef+1;
	      for (l=numinterbef;l<numinterafter;l++){
		if (masnod[l][0]!=-1){
		  fprintf (Out,"\n %ld   %le %le %le",nn+del,xi[l],yi[l],zi[l]);
		  del++;
		  fprintf (Out,"   %ld",0-masnod[l][0]);
		  for (j=0;j<masnod[l][0];j++){
		    fprintf (Out," %ld",masnod[l][j+1]);
		  }
		  fprintf (Out,"   %le %le %le",xxi[j],yyi[j],zzi[j]);
		}
	      }
	      

	    }//  end of the statement if (numinterafter-numinterbef > 2){
	    
          }//  end of the statement if (dist<Gtm->circumrad[k]){
        }//  end of statement if (give_elem_type (k) == linearhex){
      }//  end of loop for (k=0;k<ne;k++){
    }
  }//  end of loop for (i=0;i<ne;i++){
  

  fprintf (stdout,"\n\n VYSLEDKY \n");
  for (i=0;i<numinter;i++){
    fprintf (stdout,"\n %ld ",i);
    for (j=0;j<8;j++){
      fprintf (stdout,"  %ld",masnod[i][j]+1);
    }
  }
  fprintf (stdout,"\n\n SOURADNICE \n");
  for (i=0;i<numinter;i++){
    fprintf (stdout,"\n %ld   %le %le %le   %le %le %le",i,xi[i],yi[i],zi[i],xxi[i],yyi[i],zzi[i]);
  }
  
  
}
