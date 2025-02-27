#include "seqtop.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include "galias.h"
#include "iotools.h"
#include "gtopology.h"

/**
   constructor
   
   @param nd - number of subdomains
   @param meshd - mesh description
   
   JK, 11.9.2007
*/
seqtop::seqtop (long nd,meshdescription meshd)
{
  //  number of subdomains
  ns=nd;
  //  mesh description
  md = meshd;


  //  total number of nodes in the whole problem
  tnnp=0;
  //  total number of interface (boundary) nodes
  tnbn=0;
  //  number of nodes in the class gtopology
  nn=0;
  //  total number of all elements in the problem
  tnep=0;
  
  
  //  array containing number of nodes on subdomains
  nnsd=NULL;
  //  array containing number of elements on subdomains
  //nesd=NULL;

  //  array containing list of numbers of boundary nodes on subdomains
  nbnd=NULL;

  //  array containing list of numbers of internal nodes on subdomains
  nind=NULL;
  
  //  array containing number of subdomains which each node belongs to
  amultip=NULL;

  //  array containing number of subdomains which each boundary/interface node belongs to
  bmultip=NULL;

  //  array containing number of subdomains which share the nodes
  nodmultip=NULL;

  //  local to global correspondence
  ltg=NULL;

  //  local to global correspondence for elements
  //eltg=NULL;


  //  array containing local numbers of interface/boundary nodes
  lnbn=NULL;
  
  //  array containing local numbers of internal nodes
  lnin=NULL;
  
  //  array containing global glued numbers of boundary/interface nodes
  ggnbn=NULL;
  
  //  array containing global glued numbers of internal nodes
  ggnin=NULL;
  
  //  array containing coarse/interface numbers of boundary/interface nodes
  icnbnmas=NULL;

  //  array containing node multiplicity of boundary/interface nodes
  icmultip=NULL;
  
  
  //  local numbers of boundary/interface nodes appropriate to coarse node
  lnbncn=NULL;
  
  //  global glued numbers of boundary/interface nodes appropriate to coarse node
  ggnbncn=NULL;
  
  //  subdomain id of interface/boundary nodes appropriate to coarse node
  sid=NULL;
  
  eldom=NULL;

  //  array containing the number of elements on subdomains
  ned=NULL;

  //  domain-element correspondence
  domel=NULL;
  
  
  //  the number of coupled DOFs
  ncdof=0;
  
  //  array containing indicators of coupled DOFs
  coupdof=NULL;

  //  array containing suspicious indicators of coupled DOFs
  coupdofmas=NULL;
  
  //  array containing DOF indicators
  dofind = NULL;
  
}

/**
   destructor
   
   JK, 11.9.2007
*/
seqtop::~seqtop()
{
  long i;
  
  //  array containing number of nodes on subdomains
  if (nnsd!=NULL)
    delete [] nnsd;
  /*
  //  array containing number of elements on subdomains
  if (nesd!=NULL)
    delete [] nesd;
    */

  //  array containing list of numbers of boundary nodes on subdomains
  if (nbnd!=NULL)
    delete [] nbnd;
  
  //  array containing list of numbers of internal nodes on subdomains
  if (nind!=NULL)
    delete [] nind;


  //  array containing number of subdomains which each node belongs to
  if (amultip!=NULL)
    delete [] amultip;

  //  array containing number of subdomains which each boundary/interface node belongs to
  if (bmultip!=NULL)
    delete [] bmultip;

  //  array containing number of subdomains which share the nodes
  if (nodmultip!=NULL){
    for (i=0;i<ns;i++){
      delete [] nodmultip[i];
    }
    delete [] nodmultip;
  }
  
  //  local to global correspondence
  if (ltg!=NULL){
    for (i=0;i<ns;i++){
      delete [] ltg[i];
    }
    delete [] ltg;
  }
  
  /*
  //  local to global correspondence for elements
  if (eltg!=NULL){
    for (i=0;i<ns;i++){
      delete [] eltg[i];
    }
    delete [] eltg;
  }
  */
  
  //  array containing local numbers of interface/boundary nodes
  if (lnbn!=NULL){
    for (i=0;i<ns;i++){
      delete [] lnbn[i];
    }
    delete [] lnbn;
  }
  
  //  array containing local numbers of internal nodes
  if (lnin!=NULL){
    for (i=0;i<ns;i++){
      delete [] lnin[i];
    }
    delete [] lnin;
  }

  //  array containing list of global glued numbers of interface/boundary nodes
  if (ggnbn!=NULL){
    for (i=0;i<ns;i++){
      delete [] ggnbn[i];
    }
    delete [] ggnbn;
  }

  //  array containing list of global glued numbers of internal nodes
  if (ggnin!=NULL){
    for (i=0;i<ns;i++){
      delete [] ggnin[i];
    }
    delete [] ggnin;
  }

  //  array containing coarse/interface numbers of boundary/interface nodes
  if (icnbnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] icnbnmas[i];
    }
    delete [] icnbnmas;
  }
  
  //  array containing node multiplicity of boundary/interface nodes
  if (icmultip!=NULL)
    delete [] icmultip;



  //  local numbers of boundary/interface nodes appropriate to coarse node
  if (lnbncn!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] lnbncn[i];
    }
    delete [] lnbncn;
  }
  
  //  global glued numbers of boundary/interface nodes appropriate to coarse node
  if (ggnbncn!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] ggnbncn[i];
    }
    delete [] ggnbncn;
  }
  
  //  subdomain id of interface/boundary nodes appropriate to coarse node
  if (sid!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] sid[i];
    }
    delete [] sid;
  }
  
  if (eldom!=NULL)
    delete [] eldom;
  
  //  array containing the number of elements on subdomains
  if (ned!=NULL)
    delete [] ned;
  
  if (domel!=NULL){
    for (i=0;i<ns;i++){
      delete [] domel[i];
    }
    delete domel;
  }
  
  //  array containing number of coupled DOFs on subdomains
  if (coupdof!=NULL){
    for (i=0;i<ns;i++){
      delete [] coupdof[i];
    }
    delete [] coupdof;
  }
  
  //  array containing suspicious indicators of coupled DOFs
  if (coupdofmas!=NULL){
    delete [] coupdofmas;
  }
  
  //  array containing DOF indicators
  for (i=0;i<nn;i++){
    delete [] dofind[i];
  }
  delete [] dofind;
 
}


/**
   function reads array containing numbers of nodes on subdomains
   
   @param in - input file
   
   JK, 11.9.2007
*/
void seqtop::read_nnsd (XFILE *in)
{
  long i;

  if (nnsd!=NULL){
    delete [] nnsd;
  }
  nnsd = new long [ns];
  
  for (i=0;i<ns;i++){
    xfscanf (in,"%ld",nnsd+i);
  }
  
}

/**
   function reads array containing numbers of elements on subdomains
   
   @param in - input file
   
   JK, 25.3.2011
*/
/*
void seqtop::read_nesd (XFILE *in)
{
  long i;

  if (nesd!=NULL){
    delete [] nesd;
  }
  nesd = new long [ns];
  
  for (i=0;i<ns;i++){
    xfscanf (in,"%ld",nesd+i);
  }
}
*/

/**
   function reads array ltg (local to global map)
   
   @param in - input file
   
   11.9.2007, JK
*/
void seqtop::read_ltg (XFILE *in)
{
  long i,j;
  long nn;

  if (ltg!=NULL){
    for (i=0;i<ns;i++){
      delete [] ltg[i];
    }
    delete [] ltg;
  }
  ltg = new long* [ns];
  for (i=0;i<ns;i++){
    ltg[i] = new long [nnsd[i]];
  }
  
  
  if (md==metis){
    //  in this case, there is no postprocessing of METIS output
    //  the input file contains the subdomain id for each node
    //  this case is used especially for the BOSS preconditioning
    nn=0;
    for (i=0;i<ns;i++){
      nn+=nnsd[i];
      nnsd[i]=0;
    }
    
    for (i=0;i<nn;i++){
      xfscanf (in,"%ld",&j);
      ltg[j][nnsd[j]]=i;
      nnsd[j]++;
    }
  }
  else{
    //  in this case, there is a postprocessing of METIS output
    //  boundary/interface nodes are ordered or all nodes have
    //  their global numbers
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	xfscanf (in,"%ld",&ltg[i][j]);
	ltg[i][j]--;
      }
    }
  }
}

/**
   function reads array eltg (local to global map for elements)
   
   @param in - input file
   
   JK, 25.3.2011
*/
/*
void seqtop::read_eltg (XFILE *in)
{
  long i,j;
  long ne;
  
  if (eltg!=NULL){
    for (i=0;i<ns;i++){
      delete [] eltg[i];
    }
    delete [] eltg;
  }
  eltg = new long* [ns];
  for (i=0;i<ns;i++){
    eltg[i] = new long [nesd[i]];
  }
  
  if (md==metis_elem){
    //  in this case, there is no postprocessing of METIS output
    //  the input file contains the subdomain id for each element
    //  this case is used especially for the multilevel computations
    //  heat transfer with CEMHYD, homogenization, etc.
    ne=0;
    for (i=0;i<ns;i++){
      ne+=nesd[i];
      nesd[i]=0;
    }
    
    for (i=0;i<ne;i++){
      xfscanf (in,"%ld",&j);
      eltg[j][nesd[j]]=i;
      nesd[j]++;
    }
    
  }

}
*/

/**
   function reads correspondence between elements and subdomains/aggregates
   
   @param ne - the number of all elements in the problem
   @param in - input file stream
   
   JK, 2.4.2011
*/
void seqtop::read_eldom (long ne,XFILE *in)
{
  long i;
  
  //  total number of elements in problem
  tnep=ne;
  
  if (eldom!=NULL){
    delete [] eldom;
  }
  eldom = new long [ne];
  for (i=0;i<ne;i++){
    xfscanf (in,"%ld",eldom+i);
  }
}


/**
   function assembles the array bmultip or amultip
   function computes the variables tnnp or tnbn
   
   @param out - output file for auxiliary print
   
   14.7.2009, JK
*/
void seqtop::assemble_multip (FILE *out)
{
  long i,j,k;

  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes

    //  total number of boundary/interface nodes
    tnbn=0;
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (tnbn<ltg[i][j])
	  tnbn=ltg[i][j];
      }
    }
    tnbn++;
    
    if (bmultip!=NULL){
      delete [] bmultip;
    }
    bmultip = new long [tnbn];
    for(i=0;i<tnbn;i++){
      bmultip[i]=0;
    }
    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	k=ltg[i][j];
	if (k>-1)
	  bmultip[k]++;
      }
    }
    
    break;
  }
    
  case metis:
  case all_nodes:{
    
    //  total number of nodes in the problem
    tnnp=0;
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (tnnp<ltg[i][j])
	  tnnp=ltg[i][j];
      }
    }
    //  must be increased due to indices from 0 instead of 1
    tnnp++;
    
    if (amultip != NULL)
      delete [] amultip;
    amultip = new long [tnnp];
    for (i=0;i<tnnp;i++){
      amultip[i]=0;
    }
    
    //  computation of node incidences
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	amultip[ltg[i][j]]++;
      }
    }
    
    break;
  }

  case neg_bound_nodes:{
    
    tnnp=0;
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0){
	  if (tnnp<0-ltg[i][j]-1)
	    tnnp=0-ltg[i][j]-2;
	}
	if (tnnp<ltg[i][j])
	  tnnp=ltg[i][j];
      }
    }
    //  must be increased due to indices from 0 instead of 1
    tnnp++;
    
    if (amultip != NULL)
      delete [] amultip;
    amultip = new long [tnnp];
    for (i=0;i<tnnp;i++){
      amultip[i]=0;
    }
    
    //  computation of node incidences / node multiplicity
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0)
	  amultip[0-ltg[i][j]-2]++;
	else
	  amultip[ltg[i][j]]++;
      }
    }
    
    break;
  }

  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  if (amultip!=NULL){
    fprintf(out,"\n\n\n tnnp - total number of nodes in whole problem is %ld\n\n",tnnp);
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %6ld  amultip %3ld",i,amultip[i]);
    }
  }
  if (bmultip!=NULL){
    fprintf(out,"\n\n\n tnbn - total number of boundary/interface nodes in whole problem is %ld\n\n",tnbn);
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n boundary/interface node %6ld  bmultip %3ld",i,bmultip[i]);
    }
  }

}

/**
   function searches for coupled DOFs which are shared
   by interfaces
   
   @param top - pointer to the general topology
   @param out - output file
   
   13.7.2009, JK
*/
void seqtop::coupled_dofs (gtopology *top,FILE *out)
{
  long i,j,k,l,ndofn,dof;
  
  //  the number of coupled DOFs
  ncdof=0;
  
  //  number of nodes in the class gtopology
  nn=top->nn;

  //  loop over the number of all nodes in the problem
  for (i=0;i<top->nn;i++){
    //  number of DOFs in node
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs in node
    for (j=0;j<ndofn;j++){
      dof=top->give_dof (i,j);
      if (dof>1){
	if (ncdof<dof)
	  ncdof=dof;
      }
    }
  }
  
  if (ncdof>0){
    
    //  array containing number of coupled DOFs on subdomains
    //  the indicator is greater than one, otherwise the
    //  code number is not coupled
    if (coupdof!=NULL){
      for (i=0;i<ns;i++){
	delete [] coupdof[i];
      }
      delete [] coupdof;
    }
    coupdof = new long* [ns];
    for (i=0;i<ns;i++){
      coupdof[i] = new long [ncdof];
      for (j=0;j<ncdof;j++){
	coupdof[i][j]=0;
      }
    }
    
    //  number of node
    //  it has value nn at the end of loops
    l=0;
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	//  number of DOFs in node
	ndofn=top->give_ndofn (l);
	//  loop over the number of DOFs
	for (k=0;k<ndofn;k++){
	  dof=top->give_dof (l,k);
	  if (dof>1)
	    coupdof[i][dof-1]++;
	}
	l++;
      }
    }
    
    //  array containing suspicious indicators of coupled DOFs
    if (coupdofmas!=NULL){
      delete [] coupdofmas;
    }
    coupdofmas = new long [ncdof];
    for (i=0;i<ncdof;i++){
      coupdofmas[i]=0;
    }
    
    //  loop over the number of coupled DOFs
    for (i=0;i<ncdof;i++){
      k=0;
      //  loop over the number of subdomains
      for (j=0;j<ns;j++){
	if (coupdof[j][i]>0)
	  k++;
      }
      if (k>1){
	//  coupled DOFs described by the i-th indicator occur at least on two subdomains
	//  such coupled DOFs have to be treated as the boundary/interface DOFs
	coupdofmas[i]=1;
      }
    }
  }

  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n ncdof - the number of coupled DOFs %ld",ncdof);
}

/**
   function updates node multiplicity

   
   function establishes:
   
   tnbn - total number of boundary nodes
   
   tnnp - total number of nodes in problem (if possible, this number cannot be
          obtained in the case of md=bound_nodes)
   
   amultip - array of all node multiplicity, it is assembled if md=all_nodes
   
   bmultip - array of interface/boundary node multiplicity, it is assembled if md=bound_nodes
   
   nodmultip - array of node multiplicity (on all subdomains)
   
   
   @param top - pointer to the general topology
   @param out - output file for auxiliary print
   
   JK, 11.9.2007
*/
void seqtop::update_multip (gtopology *top,FILE *out)
{
  long i,j,k,ndofn,ii,dofid;
  
  switch (md){
  case metis:
  case all_nodes:{
    
    
    dofind = new long* [nn];
    for (i=0;i<nn;i++){
      dofind[i] = new long [top->give_ndofn(i)];
      for (j=0;j<top->give_ndofn(i);j++){
	dofind[i][j]=0;
      }
    }
    
    // *********************************************************************************
    //  in the case of coupled DOFs, treatment with nodes is not enough
    //  DOFs have to be split to group of internal DOFs and boundary/interface DOFs
    //
    //  there are three types of boundary/interface DOFs:
    //
    //  all DOFs belonging to the boundary/interface nodes are boundary/interface DOFs
    //
    //  coupled DOFs connected to any boundary/interface node are boundary/interface DOFs
    //
    //  coupled DOFs not connected to any boundary/interface node but shared by at least
    //  two subdomains are also boundary/interface DOFs
    //
    //  all remaining DOFs are internal DOFs and can be eliminated on subdomains
    //
    //  there are two loops, the first loop determines the boundary/interface nodes
    //  and their DOFs, simultaneously, it detects with the help of the array coupdofmas
    //  DOFs belonging to at least two subdomains, such DOFs are also denoted as the boundary/interface DOFs,
    //  in the first loop, it detects also coupled DOFs connected to boundary/interface nodes
    //  and this fact is stored in the array coupdofmas
    //  second loop denotes the coupled DOFs which are shared by boundary/interface nodes
    // *********************************************************************************

    //  first loop

    //  node number
    //  it will be equal to the number of all nodes (top->nn) at the end of loops
    ii=0;
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	//  number of DOFs in node
	ndofn = top->give_ndofn (ii);
	//  loop over the number of DOFs in node
	for (k=0;k<ndofn;k++){
	  //  DOF indicator
	  dofid=top->give_dof (ii,k);
	  if (amultip[ltg[i][j]]>1){
	    //  the i-th node is boundary/interface node, therefore this is boundary/interface DOF
	    dofind[ii][k]=1;
	    if (dofid>1){
	      //  there is coupled DOF defined in this boundary/interface node
	      //  this information has to be stored in the array coupdofmas
	      //  it will be used in the following loop
	      if (coupdofmas[dofid-2]==0)
		coupdofmas[dofid-2]=1;
	    }
	  }
	  if (dofid>1){
	    //  there is coupled DOF connected to at least two subdomains
	    if (coupdofmas[dofid-2]==1){
	      //  this coupled DOF has to be taken into account
	      //  this is boundary/interface DOF
	      dofind[ii][k]=1;
	    }
	  }
	}
	ii++;
      }
    }
    if (ii!=nn){
      print_err("the number of nodes is not equal to the last index after loops", __FILE__, __LINE__, __func__);
    }
    
    //  second loop
    /*
    //  node number
    //  it will be equal to the number of all nodes (top->nn) at the end of loops
    ii=0;
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	if (amultip[ltg[i][j]]>1){
	  //  this is boundary/interface node
	  //  its contribution was taken into account in the previous loop
	}else{
	  //  number of DOFs in node
	  ndofn = top->give_ndofn (ii);
	  //  loop over the number of DOFs in node
	  for (k=0;k<ndofn;k++){
	    //  DOF indicator
	    dofid=top->give_dof (ii,k);
	    if (dofid>1){
	      //  there is coupled DOF
	      if (coupdofmas[dofid-2]==1){
		//  this coupled DOF has to be taken into account
		//  this is boundary/interface DOF
		dofind[ii][k]=1;
	      }
	    }
	  }
	}
	ii++;
      }
    }
    if (ii!=nn){
      print_err("the number of nodes is not equal to the last index after loops", __FILE__, __LINE__, __func__);
    }
    */
    
    // *********************************************************************************
    //  modification of the array amultip
    //  the array has to be modified with respect to the array dofind
    //  if any DOF in node is denoted as the boundary/interface DOF, the appropriate node
    //  has to be denoted as the boundary/interface node, this fact is
    //  assured with the help of the array amultip, where value 1 is rewritten to the value 2
    // *********************************************************************************
    
    //  node number
    //  it will be equal to the number of all nodes (top->nn) at the end of loops
    ii=0;
    //  total number of boundary/interface nodes
    tnbn=0;
    //  loop over the number of subdomain
    for (i=0;i<ns;i++){
      //  loop over the number of nodes on subdomain
      for (j=0;j<nnsd[i];j++){
	//  number of DOFs in node
	ndofn = top->give_ndofn (ii);
	//  loop over the number of DOFs in node
	for (k=0;k<ndofn;k++){
	  if (dofind[ii][k]==1){
	    //  the DOF is boundary/interface DOF
	    if (amultip[ltg[i][j]]==1){
	      //  the node is internal node
	      //  it has to be changed from internal to boundary/interface
	      amultip[ltg[i][j]]=2;
	    }
	  }
	}
	ii++;
      }
    }
    if (ii!=nn){
      print_err("the number of nodes is not equal to the last index after loops", __FILE__, __LINE__, __func__);
    }
    // *********************************************************************************
    //  end of modification of the array amultip
    // *********************************************************************************
    
    break;
  }
    
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  if (amultip!=NULL){
    fprintf(out,"\n\n\n tnnp - total number of nodes in whole problem is %ld\n\n",tnnp);
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %6ld  amultip %3ld",i,amultip[i]);
    }
  }
  
}

/**
   function assembles the arrays nbnd and nind
   
   @param out - output file for auxiliary print
   
   JK, 14.7.2009
*/
void seqtop::assemble_nbnd_nind (FILE *out)
{
  long i,j,k;
  
  //  numbers of boundary nodes on subdomains
  if (nbnd != NULL)
    delete [] nbnd;
  nbnd = new long [ns];
  //  numbers of internal nodes on subdomains
  if (nind!=NULL)
    delete [] nind;
  nind = new long [ns];
  
  
  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes

    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]>-1)
	  nbnd[i]++;
      }
    }
    
    break;
  }


  case metis:
  case all_nodes:{
    
    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	k=amultip[ltg[i][j]];
	if (k>1)
	  nbnd[i]++;
      }
    }
    
    //  total number of boundary/interface nodes
    tnbn=0;
    //  loop over the number of all problems in the problem
    for (i=0;i<tnnp;i++){
      if (amultip[i]>1)
	tnbn++;
    }
    
    break;
  }
    
  case neg_bound_nodes:{
    
    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0)
	  k=amultip[0-ltg[i][j]-2];
	else
	  k=amultip[ltg[i][j]];
	if (k>1)
	  nbnd[i]++;
      }
    }
    
    //  total number of boundary/interface nodes
    tnbn=0;
    //  loop over the number of all problems in the problem
    for (i=0;i<tnnp;i++){
      if (amultip[i]>1)
	tnbn++;
    }
    
    break;
  }
    
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  for (i=0;i<ns;i++){
    nind[i]=nnsd[i]-nbnd[i];
  }
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf(out,"\n\n\n total number of boundary nodes in whole problem is %ld\n\n",tnbn);

  fprintf (out,"\n\n\n numbers of interface/boundary nodes on each subdomain \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nbnd[i]);
  }
  
  fprintf (out,"\n\n\n numbers of internal nodes on each subdomain \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nind[i]);
  }
  
}


/**
   function assembless the array nodmultip
   
   @param out - output file for auxiliary print
   
   JK, 11.9.2007
*/
void seqtop::assemble_nodmultip (FILE *out)
{
  long i,j,k;
  
  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes

    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }
    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	k=ltg[i][j];
	if (k>-1){
	  nodmultip[i][j]=bmultip[k];
	}
	else{
	  nodmultip[i][j]=1;
	}
      }
    }
    
       
    break;
  }


  case metis:
  case all_nodes:{
    
    //  nodal multiplicity
    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }

    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	k=amultip[ltg[i][j]];
	nodmultip[i][j]=k;
      }
    }
    
    break;
  }

  case neg_bound_nodes:{
    
    //  nodal multiplicity
    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }
    
    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0)
	  k=amultip[0-ltg[i][j]-2];
	else
	  k=amultip[ltg[i][j]];
	nodmultip[i][j]=k;
      }
    }
    
    break;
  }

  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf(out,"\n\n\n Multiplicity of nodes on subdomain\n\n");  
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      fprintf(out,"%4ld  %6ld   %ld\n",i,j,nodmultip[i][j]);  
    }
  }
}

/**
   function assembless the array dofind if it has not been assembled yet
   
   @param top - pointer to the general topology
   @param out - output file for auxiliary print
   
   JK, 15.7.2009
*/
void seqtop::assemble_dofind (gtopology *top,FILE */*out*/)
{
  long i,j,k,ndofn,ii;
  
  if (dofind==NULL){

    //  array containing DOF indicators
    //  it indicates whether the DOF is internal or boundary/interface
    dofind = new long* [nn];
    for (i=0;i<nn;i++){
      dofind[i] = new long [top->give_ndofn(i)];
      for (j=0;j<top->give_ndofn(i);j++){
	dofind[i][j]=0;
      }
    }
    
    //  node number
    //  it will be equal to the number of all nodes (top->nn) at the end of loops
    ii=0;
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	//  number of DOFs in node
	ndofn=top->give_ndofn(ii);
	if (nodmultip[i][j]>1){
	  //  this is boundary/interface node

	  //  loop over the number of DOFs in node
	  for (k=0;k<ndofn;k++){
	    dofind[ii][k]=1;
	  }
	}else{
	  //  this is internal node
	  
	  //  loop over the number of DOFs in node
	  for (k=0;k<ndofn;k++){
	    dofind[ii][k]=0;
	  }
	}
	ii++;
      }
    }
  }
}


void seqtop::compute_multiplicity (gtopology *top,FILE *out)
{
  assemble_multip (out);
  coupled_dofs (top,out);
  if (ncdof>0){
    //  there are coupled DOFs in the problem
    //  array amultip has to be modified with respect to coupled DOFs
    update_multip (top,out);
  }
  assemble_nbnd_nind (out);
  assemble_nodmultip (out);
  assemble_dofind (top,out);
}





/**
   function assembles local numbers of nodes
   
   the following arrays are assembled:
   
   lnbn - local numbers of interface/boundary nodes
   
   lnin - local numbers of internal nodes
   
   @param out - output file for auxiliary print
   
   JK, 24.5.2009
*/
void seqtop::node_local_numbers (FILE *out)
{
  long i,j,ii,jj;
  
  //  array containing local numbers of internal nodes
  if (lnin != NULL){
    for (i=0;i<ns;i++){
      delete [] lnin[i];
    }
    delete [] lnin;
  }
  lnin = new long* [ns];
  for (i=0;i<ns;i++){
    lnin[i] = new long [nind[i]];
  }
  //  array containing local numbers of interface/boundary nodes
  if (lnbn != NULL){
    for (i=0;i<ns;i++){
      delete [] lnbn[i];
    }
    delete [] lnbn;
  }
  lnbn = new long* [ns];
  for (i=0;i<ns;i++){
    lnbn[i] = new long [nbnd[i]];
  }
  
  switch (md){
  case bound_nodes:{
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]>-1){
	  lnbn[i][ii]=j;
	  ii++;
	}
	if (ltg[i][j]==-1){
	  lnin[i][jj]=j;
	  jj++;
	}
      }
    }
    break;
  }
  case metis:
  case all_nodes:{
    
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (amultip[ltg[i][j]]>1){
	  lnbn[i][ii]=j;
	  ii++;
	}
	if (amultip[ltg[i][j]]==1){
	  lnin[i][jj]=j;
	  jj++;
	}
      }
    }
    
    break;
  }
  case neg_bound_nodes:{
    
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0){
	  lnbn[i][ii]=j;
	  ii++;
	}
	if (ltg[i][j]>-1){
	  lnin[i][jj]=j;
	  jj++;
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n local numbers of interface/boundary nodes on subdomain (lnbn)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nbnd[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,lnbn[i][j]+1);
    }
  }
  
  fprintf (out,"\n\n\n local numbers of internal nodes on subdomain (lnin)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nind[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,lnin[i][j]+1);
    }
  }
  
}


/**
   function assembles global glued numbers of nodes
   
   the following arrays are assembled:
   
   ggnbn - global glued numbers of interface/boundary nodes
   
   ggnin - global glued numbers of internal nodes
   
   @param out - output file for auxiliary print
   
   JK, 24.5.2009
*/
void seqtop::node_global_glued_numbers (FILE *out)
{
  long i,j,ii,jj,kk;
  
  //  array containing global glued numbers of internal nodes
  if (ggnin != NULL){
    for (i=0;i<ns;i++){
      delete [] ggnin[i];
    }
    delete [] ggnin;
  }
  ggnin = new long* [ns];
  for (i=0;i<ns;i++){
    ggnin[i] = new long [nind[i]];
  }
  //  array containing global glued numbers of interface/boundary nodes
  if (ggnbn != NULL){
    for (i=0;i<ns;i++){
      delete [] ggnbn[i];
    }
    delete [] ggnbn;
  }
  ggnbn = new long* [ns];
  for (i=0;i<ns;i++){
    ggnbn[i] = new long [nbnd[i]];
  }
  
  switch (md){
  case bound_nodes:{
    kk=0;
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]>-1){
	  ggnbn[i][ii]=kk;
	  ii++;  kk++;
	}
	if (ltg[i][j]==-1){
	  ggnin[i][jj]=kk;
	  jj++;  kk++;
	}
      }
    }
    break;
  }
  case metis:
  case all_nodes:{
    
    kk=0;
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (amultip[ltg[i][j]]>1){
	  ggnbn[i][ii]=kk;
	  ii++;  kk++;
	}
	if (amultip[ltg[i][j]]==1){
	  ggnin[i][jj]=kk;
	  jj++;  kk++;
	}
      }
    }
    
    break;
  }
  case neg_bound_nodes:{
    
    kk=0;
    for (i=0;i<ns;i++){
      ii=0;  jj=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0){
	  ggnbn[i][ii]=kk;
	  ii++;  kk++;
	}
	if (ltg[i][j]>-1){
	  ggnin[i][jj]=kk;
	  jj++;  kk++;
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n global glued numbers of interface/boundary nodes on subdomain (ggnbn)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nbnd[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,ggnbn[i][j]+1);
    }
  }
  
  fprintf (out,"\n\n\n global glued numbers of internal nodes on subdomain (ggnin)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nind[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,ggnin[i][j]+1);
    }
  }
  
}


/**
   function assembles coarse numbers of nodes
   
   the following arrays are assembled:
   
   icnbnmas - coarse numbers of interface/boundary nodes
   icmultip - number of multiplicity of boundary/interface nodes

   @param out - output file for auxiliary print
   
   JK, 5.7.2009
*/
void seqtop::node_coarse_numbers (FILE *out)
{
  long i,j,k;
  long *av;
  av=NULL;
  
  //  array containing coarse numbers of interface/boundary nodes
  if (icnbnmas != NULL){
    for (i=0;i<ns;i++){
      delete [] icnbnmas[i];
    }
    delete [] icnbnmas;
  }
  icnbnmas = new long* [ns];
  for (i=0;i<ns;i++){
    icnbnmas[i] = new long [nbnd[i]];
  }
  
  
  switch (md){
  case bound_nodes:{
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[i];j++){
	icnbnmas[i][j]=ltg[i][lnbn[i][j]];
      }
    }
    break;
  }
  case metis:
  case all_nodes:
  case neg_bound_nodes:{
    j=0;
    av = new long [tnnp];
    for (i=0;i<tnnp;i++){
      av[i]=amultip[i];
      if (av[i]==1)
	av[i]=-1;
      if (av[i]>1){
	av[i]=j;
	j++;
      }
    }
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }

  switch (md){
  case bound_nodes:{
    break;
  }
  case metis:
  case all_nodes:{
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[i];j++){
	icnbnmas[i][j]=av[ltg[i][lnbn[i][j]]];
      }
    }
    break;
  }
  case neg_bound_nodes:{
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[i];j++){
	k=0-ltg[i][lnbn[i][j]]-2;
	icnbnmas[i][j]=av[k];
      }
    }
    break;
  }  
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
  
  //  number of multiplicity of boundary/interface nodes
  if (icmultip!=NULL){
    delete [] icmultip;
  }
  icmultip = new long [tnbn];
  for (i=0;i<tnbn;i++){
    icmultip[i]=0;
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of boundary/interface nodes
    for (j=0;j<nbnd[i];j++){
      icmultip[icnbnmas[i][j]]++;
    }
  }
  
  if (av!=NULL)
    delete [] av;
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n coarse numbers of interface/boundary nodes on subdomain (icnbnmas)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nbnd[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,icnbnmas[i][j]+1);
    }
  }
  
  fprintf (out,"\n\n\n coarse-local correspondence (icmultip)\n");
  for (i=0;i<tnbn;i++){
    fprintf (out,"\n coarse node %6ld, number of nodes %2ld",i,icmultip[i]);
  }
}

/**
   function assembles coarse - local numbers map
   
   the following arrays are assembled:
   
   lnbncn - local numbers of boundary/interface nodes of the coarse node
   sid - subdomain id of interface/boundary nodes of the coarse node

   @param out - output file for auxiliary print
   
   JK, 5.7.2009
*/
void seqtop::node_coarse_local_map (FILE *out)
{
  long i,j,ln,cn;
  
  //  local numbers of boundary/interface nodes of the coarse node
  if (lnbncn!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] lnbncn[i];
    }
    delete [] lnbncn;
  }
  lnbncn = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    lnbncn[i]=new long [icmultip[i]];
  }
  //  subdomain id of interface/boundary nodes of the coarse node
  if (sid!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] sid[i];
    }
    delete [] sid;
  }
  sid = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    sid[i] = new long [icmultip[i]];
  }
  
  //  loop over the number of subdomains
  for (i=0;i<tnbn;i++){
    icmultip[i]=0;
  }
  
  for (i=0;i<ns;i++){
    //  loop over the number of boundary/interface nodes
    for (j=0;j<nbnd[i];j++){
      //  coarse number
      cn=icnbnmas[i][j];
      //  local number
      ln=lnbn[i][j];
      
      //  local number of boundary/interface node
      lnbncn[cn][icmultip[cn]]=ln;
      //  subdomain id if boundary/interface node
      sid[cn][icmultip[cn]]=i;
      icmultip[cn]++;
    }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  
  fprintf (out,"\n\n\n coarse-local correspondence \n");
  for (i=0;i<tnbn;i++){
    fprintf (out,"\n coarse node %6ld, number of nodes %2ld",i,icmultip[i]);
    for (j=0;j<icmultip[i];j++){
      fprintf (out,"   node %ld, local number %6ld, subdomain id %4ld",j,lnbncn[i][j]+1,sid[i][j]+1);
    }
  }
  fprintf (out,"\n\n\n coarse-local correspondence \n");
  for (i=0;i<tnbn;i++){
    fprintf (out,"\n coarse node %6ld",i);
    for (j=0;j<icmultip[i];j++){
      fprintf (out,"    loc.n. %6ld, sub. %4ld",lnbncn[i][j]+1,sid[i][j]+1);
    }
  }

}


/**
   function assembles coarse - global glued numbers map
   
   the following arrays are assembled:
   
   ggnbncn - global glued numbers of boundary/interface nodes of the coarse node
   sid - subdomain id of interface/boundary nodes of the coarse node
   
   @param out - output file for auxiliary print
   
   JK, 5.7.2009
*/
void seqtop::node_coarse_global_glued_map (FILE *out)
{
  long i,j,ln,cn;
  
  //  global glued numbers of boundary/interface nodes of the coarse node
  if (ggnbncn!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] ggnbncn[i];
    }
    delete [] ggnbncn;
  }
  ggnbncn = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    ggnbncn[i]=new long [icmultip[i]];
  }
  //  subdomain id of interface/boundary nodes of the coarse node
  if (sid!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] sid[i];
    }
    delete [] sid;
  }
  sid = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    sid[i] = new long [icmultip[i]];
  }
  
  //  loop over the number of subdomains
  for (i=0;i<tnbn;i++){
    icmultip[i]=0;
  }
  
  for (i=0;i<ns;i++){
    //  loop over the number of boundary/interface nodes
    for (j=0;j<nbnd[i];j++){
      //  coarse number
      cn=icnbnmas[i][j];
      //  global glued number
      ln=ggnbn[i][j];
      
      //  global glued number of boundary/interface node
      ggnbncn[cn][icmultip[cn]]=ln;
      //  subdomain id if boundary/interface node
      sid[cn][icmultip[cn]]=i;
      icmultip[cn]++;
    }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n coarse-global glued correspondence \n");
  for (i=0;i<tnbn;i++){
    fprintf (out,"\n coarse node %6ld, number of nodes %2ld",i,icmultip[i]);
    for (j=0;j<icmultip[i];j++){
      fprintf (out,"   node %ld, global glued number %6ld, subdomain id %4ld",j,ggnbncn[i][j]+1,sid[i][j]+1);
    }
  }
  fprintf (out,"\n\n\n coarse-global glued correspondence \n");
  for (i=0;i<tnbn;i++){
    fprintf (out,"\n coarse node %6ld",i);
    for (j=0;j<icmultip[i];j++){
      fprintf (out,"   ggn %6ld, sub.%4ld",ggnbncn[i][j]+1,sid[i][j]+1);
    }
  }

}


/**
   function generates the Schur complement ordering
   
   ordering of subdomains is performed
   ordering of the coarse problem is in seqselnodes::schur_ordering
   
   @param top - pointer to general topology
   @param out - output file
   
   JK, 13.7.2009
*/
long seqtop::schur_ordering (gtopology *top,FILE *out)
{
  long i,j,k,l,ndof,ndofn,ncdof,ggn;
  long *aux;
  
  //  searching of maximum code number indicator
  ncdof=0;
  //  loop over the number of all nodes
  for (i=0;i<nn;i++){
    //  number of DOFs in node
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs in node
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>ncdof)  ncdof=k;
    }
  }
  ncdof--;
  if (ncdof<0)  ncdof=0;
  aux = new long [ncdof];
  for (i=0;i<ncdof;i++){
    aux[i]=-1;
  }
  
  // ******************************
  //  generation of internal DOFs
  // ******************************

  //  number of actual DOF
  //  it will be equal to the number of DOFs at the end of loops
  ndof=1;
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    for (j=0;j<ncdof;j++){
      aux[j]=-1;
    }
    //  loop over the number of internal nodes on subdomains
    for (j=0;j<nind[i];j++){
      //  global glued number of internal node
      ggn=ggnin[i][j];
      //  number of DOF
      ndofn=top->give_ndofn (ggn);
      //  loop over the number of DOFs in node
      for (l=0;l<ndofn;l++){
	if (dofind[ggn][l]==0){
	  //  this is internal DOF
	  
	  k=top->give_dof (ggn,l);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    top->gnodes[ggn].cn[l]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      top->gnodes[ggn].cn[l]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      top->gnodes[ggn].cn[l]=aux[k-2];
	    }
	  }
	}
      }
    }
  }
  
  // ****************************************
  //  generation of boundary/interface DOFs
  // ****************************************

  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    for (j=0;j<ncdof;j++){
      aux[j]=-1;
    }
    //  loop over the number of internal nodes on subdomains
    for (j=0;j<nbnd[i];j++){
      //  global glued number of boundary/interface node
      ggn=ggnbn[i][j];
      //  number of DOF
      ndofn=top->give_ndofn (ggn);
      //  loop over the number of DOFs in node
      for (l=0;l<ndofn;l++){
	if (dofind[ggn][l]==1){
	  //  this is boundary/interface DOF
	  
	  k=top->give_dof (ggn,l);
	  if (k<0)  continue;
	  if (k==0)  continue;
	  if (k==1){
	    top->gnodes[ggn].cn[l]=ndof;  ndof++;
	  }
	  if (k>1){
	    if (aux[k-2]==-1){
	      top->gnodes[ggn].cn[l]=ndof;
	      aux[k-2]=ndof;
	      ndof++;
	    }
	    else{
	      top->gnodes[ggn].cn[l]=aux[k-2];
	    }
	  }
	}
      }
    }
  }
  ndof--;
  
  delete [] aux;
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n kontrola kodovych cisel \n");
  for (i=0;i<nn;i++){
    fprintf (out,"%4ld    %4ld %4ld\n",i,top->give_dof (i,0),top->give_dof (i,1));
  }
  fprintf (out,"\n");
  
  //  code numbers were generated
  top->cnstate=1;

  return ndof;
}

/**
   function generates the Schur complement ordering
   
   ordering of subdomains is performed
   ordering of the coarse problem is in seqselnodes::schur_ordering
   
   @param top - pointer to general topology
   @param out - output file
   
   JK, 13.7.2009
*/
/*
long seqtop::schur_ordering_old (gtopology *top,FILE *out)
{
  long i,j,k,l,ndof,ndofn;
  long *aux;
  
  //  searching of maximum code number indicator
  ndof=0;
  //  loop over the number of all nodes
  for (i=0;i<nn;i++){
    //  number of DOFs in node
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs in node
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>ndof)  ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  aux = new long [ndof];
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }
  
  // ******************************
  //  generation of internal DOFs
  // ******************************

  //  number of actual DOF
  //  it will be equal to the number of DOFs at the end of loops
  ndof=1;
  //  loop over the number of all nodes
  for (i=0;i<nn;i++){
    //  number of DOF
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs in node
    for (j=0;j<ndofn;j++){
      if (dofind[i][j]==0){
	//  this is internal DOF
	
	k=top->give_dof (i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  top->gnodes[i].cn[j]=ndof;  ndof++;
	}
	if (k>1){
	  if (aux[k-2]==-1){
	    top->gnodes[i].cn[j]=ndof;
	    aux[k-2]=ndof;
	    ndof++;
	  }
	  else{
	    top->gnodes[i].cn[j]=aux[k-2];
	  }
	}
      }
    }
  }
  
  // ****************************************
  //  generation of boundary/interface DOFs
  // ****************************************

  //  loop over the number of all nodes
  for (i=0;i<nn;i++){
    //  number of DOF
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs in node
    for (j=0;j<ndofn;j++){
      if (dofind[i][j]==1){
	//  this is boundary/interface DOF
	
	k=top->give_dof (i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  top->gnodes[i].cn[j]=ndof;  ndof++;
	}
	if (k>1){
	  if (aux[k-2]==-1){
	    top->gnodes[i].cn[j]=ndof;
	    aux[k-2]=ndof;
	    ndof++;
	  }
	  else{
	    top->gnodes[i].cn[j]=aux[k-2];
	  }
	}
      }
    }
  }
  ndof--;
  
  delete [] aux;
  
  return ndof;
}
*/

/**
   function rewrites array ltg
   
   JK, 14.9.2007
*/
/*
void seqtop::rewrite_ltg ()
{
  long i,j;
  
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      ltg[i][j]=-1;
    }
  }
  for (i=0;i<ns;i++){
    for (j=0;j<nbnd[i];j++){
      ltg[i][lnbn[i][j]]=lgnbn[i][j];
    }
  }
}
*/

/**
   function renumbers code numbers
   
   JK, 9.10.2007
*/
/*
void seqtop::codnum_renumb (gtopology *top)
{
  long i,j,k,ndofn,dof,nid;
  long *red;
  red = new long [ns+1];
  
  red[0]=0;
  nid=0;
  for (i=0;i<ns;i++){
    red[i+1]=0;
    for (j=0;j<nnsd[i];j++){
      ndofn = top->give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	dof = top->give_dof (nid,k);
	if (dof>red[i+1]){
	  red[i+1]=dof;
	}
      }
      nid++;
    }
  }
  
  nid=0;
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      ndofn = top->give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	dof = top->give_dof (nid,k);
	if (dof>0){
	  dof-=red[i];
	  top->save_dof (nid,k,dof);
	}
      }
      nid++;
    }
  }
  
  delete [] red;
}
*/

/**
   function determines the number of DOFs on nodes
   
   @param top - general topology
   
   JK, 13.5.2009
*/
/*
void seqtop::ndofn_list (gtopology *top)
{
  long i,j,nid;
  
  if (lndofn!=NULL){
    for (i=0;i<ns;i++){
      delete [] lndofn[i];
    }
    delete [] lndofn;
  }
  lndofn = new long* [ns];
  for (i=0;i<ns;i++){
    lndofn[i] = new long [nnsd[i]];
  }
  
  nid=0;
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      lndofn[i][j] = top->give_ndofn (nid);
      nid++;
    }
  }
  
}
*/

/**
   function determines the number of DOFs on nodes
   
   @param general topology
   
   JK, 13.5.2009
*/
 /*
void seqtop::dof_list (gtopology *top)
{
  long i,j,k,nid,ndofn;
  
  if (ldof!=NULL){
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[j];j++){
	delete [] ldof[i][j];
      }
      delete [] ldof[i];
    }
    delete [] ldof;
  }
  ldof = new long** [ns];
  for (i=0;i<ns;i++){
    ldof[i] = new long* [nnsd[i]];
    for (j=0;j<nnsd[i];j++){
      ldof[i][j] = new long [lndofn[i][j]];
    }
  }
  
  nid=0;
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      ndofn = lndofn[i][j];
      for (k=0;k<ndofn;k++){
	ldof[i][j][k] = top->give_dof (nid,k);
      }
      nid++;
    }
  }
  
}
 */



/**
   function determines the number of DOFs on interface/boundary nodes
   
   @param top - general topology
   
   JK, 13.5.2009
*/
  /*
void seqtop::boundary_ndofn_list (gtopology *top)
{

  long i,j,nid;
  
  if (lbndofn!=NULL){
    for (i=0;i<ns;i++){
      delete [] lbndofn[i];
    }
    delete [] lbndofn;
  }
  lbndofn = new long* [ns];
  for (i=0;i<ns;i++){
    lbndofn[i] = new long [nbnd[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nbnd[i];j++){
      nid=lgnbn[i][j];
      lbndofn[i][j] = top->give_ndofn (nid);
    }
  }

}
  */
/**
   function determines the number of DOFs on nodes
   
   JK, 13.5.2009
*/
   /*
void seqtop::boundary_dof_list (gtopology *top)
{
  long i,j,k,nid,ndofn;
  
  if (lbdof!=NULL){
    for (i=0;i<ns;i++){
      for (j=0;j<nbnd[j];j++){
	delete [] lbdof[i][j];
      }
      delete [] lbdof[i];
    }
    delete [] lbdof;
  }
  lbdof = new long** [ns];
  for (i=0;i<ns;i++){
    lbdof[i] = new long* [nbnd[i]];
    for (j=0;j<nbnd[i];j++){
      lbdof[i][j] = new long [lbndofn[i][j]];
    }
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nbnd[i];j++){
      nid=lgnbn[i][j];
      ndofn = lbndofn[i][j];
      for (k=0;k<ndofn;k++){
	lbdof[i][j][k] = top->give_dof (nid,k);
      }
    }
  }
  
}

   */



/**
   function selects boundary nodes
   
   
   function establishes:
   
   lnbn - local numbers of boundary nodes
   
   @param out - output stream
   
   JK, 14.9.2007
*/
/*
void seqtop::find_boundary_nodes (FILE *out)
{
  long i,j,k,*buff;
  
  //  list of global numbers of interface/boundary nodes
  if (lgnbn!=NULL){
    for (i=0;i<ns;i++){
      delete [] lgnbn[i];
    }
    delete [] lgnbn;
  }
  lgnbn = new long* [ns];
  for (i=0;i<ns;i++){
    lgnbn[i] = new long [nbnd[i]];
  }
  
  
  switch (md){
  case bound_nodes:{
    for (i=0;i<ns;i++){
      k=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]>-1){
	  lgnbn[i][k]=ltg[i][j];
	  k++;
	}
      }
    }
    break;
  }
  case metis:
  case all_nodes:
  case neg_bound_nodes:{
    
    for (i=0;i<ns;i++){
      k=0;
      for (j=0;j<nnsd[i];j++){
	if (multip[ltg[i][j]]>1){
	  lgnbn[i][k]=ltg[i][j];
	  k++;
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }


  // *******************************
  //  indication of boundary nodes 
  // *******************************
  //  local numbers of boundary nodes
  if (lnbn != NULL){
    for (i=0;i<ns;i++){
      delete [] lnbn[i];
    }
    delete [] lnbn;
  }
  lnbn = new long* [ns];
  for (i=0;i<ns;i++){
    lnbn[i] = new long [nbnd[i]];
  }
  
  switch (md){
  case metis:
  case all_nodes:
  case neg_bound_nodes:{
    for (i=0;i<ns;i++){
      k=0;
      for (j=0;j<nnsd[i];j++){
	if (nodmultip[i][j]>1){
	  lnbn[i][k]=j;
	  k++;
	}
      }
    }
    break;
  }
    
  case bound_nodes:{
    for (i=0;i<ns;i++){
      k=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]>-1){
	  lnbn[i][k]=j;
	  k++;
	}
      }
    }
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  if (nind!=NULL){
    delete [] nind;
  }
  nind = new long [ns];

  //  number of internal nodes
  for (i=0;i<ns;i++){
    nind[i] = nnsd[i]-nbnd[i];
  }
  
  //  array containing local numbers of internal nodes
  if (lnin != NULL){
    for (i=0;i<ns;i++){
      delete [] lnin[i];
    }
    delete [] lnin;
  }
  lnin = new long* [ns];
  for (i=0;i<ns;i++){
    lnin[i] = new long [nind[i]];
  }
  
  for (i=0;i<ns;i++){
    buff = new long [nnsd[i]];
    for (j=0;j<nnsd[i];j++){
      buff[j]=0;
    }
    for (j=0;j<nbnd[i];j++){
      k=lnbn[i][j];
      buff[k]=1;
    }
    
    k=0;
    for (j=0;j<nnsd[i];j++){
      if (buff[j]==0){
	lnin[i][k]=j;
	k++;
      }
    }
    delete [] buff;
  }
  
  fprintf (out,"\n\n\n numbers of interface nodes on subdomains \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nbnd[i]);
  }
  
  fprintf (out,"\n\n\n local numbers of interface nodes on subdomain");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n domain %ld",i);
    for (j=0;j<nbnd[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,lnbn[i][j]+1);
    }
  }
  
  fprintf (out,"\n\n\n numbers of internal nodes on subdomains \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nind[i]);
  }
  
  fprintf (out,"\n\n\n local numbers of internal nodes on subdomain");
  for (i=0;i<ns;i++){
    for (j=0;j<nind[i];j++){
      fprintf (out,"\n %6ld   %6ld",j,lnin[i][j]+1);
    }
  }
}
*/





/**
   function computes node multiplicity

   
   function establishes:
   
   tnbn - total number of boundary nodes
   
   nbnd - array containing numbers of interface/boundary nodes on subdomains

   nind - array containing numbers of internal nodes on subdomains
   
   tnnp - total number of nodes in problem (if possible, this number cannot be
          obtained in the case of md=bound_nodes)
   
   amultip - array of all node multiplicity, it is assembled if md=all_nodes
   
   bmultip - array of interface/boundary node multiplicity, it is assembled if md=bound_nodes
   
   nodmultip - array of node multiplicity (on all subdomains)
   
   
   @param out - output file for auxiliary print
   
   JK, 11.9.2007
*/
/*
void seqtop::compute_multiplicity (FILE *out)
{
  long i,j,k;
  
  //  list of number of boundary nodes on subdomains
  if (nbnd != NULL)
    delete [] nbnd;
  nbnd = new long [ns];
  
  
  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes
    tnbn=0;
    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	if (tnbn<ltg[i][j])
	  tnbn=ltg[i][j];
	if (ltg[i][j]>-1)
	  nbnd[i]++;
      }
    }
    tnbn++;
    
    if (bmultip!=NULL){
      delete [] bmultip;
    }
    bmultip = new long [tnbn];
    for(i=0;i<tnbn;i++){
      bmultip[i]=0;
    }
    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	k=ltg[i][j];
	if (k>-1)
	  bmultip[k]++;
      }
    }
    
    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }
    
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	k=ltg[i][j];
	if (k>-1){
	  nodmultip[i][j]=bmultip[k];
	}
	else{
	  nodmultip[i][j]=1;
	}
      }
    }
    
       
    break;
  }



  case metis:
  case all_nodes:{
    
    //  total number of nodes in the problem
    tnnp=0;
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (tnnp<ltg[i][j])
	  tnnp=ltg[i][j];
      }
    }
    //  must be increased due to indices from 0 instead of 1
    tnnp++;
    
    if (amultip != NULL)
      delete [] amultip;
    amultip = new long [tnnp];
    for (i=0;i<tnnp;i++){
      amultip[i]=0;
    }
    
    //  computation of node incidences
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	amultip[ltg[i][j]]++;
      }
    }
    
    tnbn=0;
    for (i=0;i<tnnp;i++){
      if (amultip[i]>1)
	tnbn++;
    }
    
    
    //  nodal multiplicity
    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }

    
    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	k=amultip[ltg[i][j]];
	nodmultip[i][j]=k;
	if (k>1)
	  nbnd[i]++;
      }
    }
    
    break;
  }

  case neg_bound_nodes:{
    
    tnnp=0;
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0){
	  if (tnnp<0-ltg[i][j]-1)
	    tnnp=0-ltg[i][j]-2;
	}
	if (tnnp<ltg[i][j])
	  tnnp=ltg[i][j];
      }
    }
    //  must be increased due to indices from 0 instead of 1
    tnnp++;
    
    if (amultip != NULL)
      delete [] amultip;
    amultip = new long [tnnp];
    for (i=0;i<tnnp;i++){
      amultip[i]=0;
    }
    
    //  computation of node incidences / node multiplicity
    for (i=0;i<ns;i++){
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0)
	  amultip[0-ltg[i][j]-2]++;
	else
	  amultip[ltg[i][j]]++;
      }
    }
    
    tnbn=0;
    for (i=0;i<tnnp;i++){
      if (amultip[i]>1)
	tnbn++;
    }
    
    
    //  nodal multiplicity
    if (nodmultip!=NULL){
      for (i=0;i<ns;i++){
	delete [] nodmultip[i];
      }
      delete [] nodmultip;
    }
    nodmultip = new long* [ns];
    for (i=0;i<ns;i++){
      nodmultip[i] = new long [nnsd[i]];
    }
    
    
    for (i=0;i<ns;i++){
      nbnd[i]=0;
      for (j=0;j<nnsd[i];j++){
	if (ltg[i][j]<0)
	  k=amultip[0-ltg[i][j]-2];
	else
	  k=amultip[ltg[i][j]];
	nodmultip[i][j]=k;
	if (k>1)
	  nbnd[i]++;
      }
    }
    
    break;
  }

  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  //  number of internal nodes on subdomains
  if (nind!=NULL)
    delete [] nind;
  nind = new long [ns];
  
  for (i=0;i<ns;i++){
    nind[i]=nnsd[i]-nbnd[i];
  }

  // *******************
  //  auxiliary output
  // *******************
  fprintf(out,"\n\n\n Multiplicity of nodes on subdomain\n\n");  
  for (i=0;i<ns;i++){
    for (j=0;j<nnsd[i];j++){
      fprintf(out,"%4ld  %6ld   %ld\n",i,j,nodmultip[i][j]);  
    }
  }
  
  fprintf(out,"\n\n\n total number of boundary nodes in whole problem is %ld\n\n",tnbn);
  
  fprintf (out,"\n\n\n numbers of interface/boundary nodes on each subdomain \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nbnd[i]);
  }
  
  fprintf (out,"\n\n\n numbers of internal nodes on each subdomain \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n %ld",nind[i]);
  }
  
}
*/

/**
  function calls several auxiliary functions
  
  @param top - pointer to the general topology
  @param out - output file
  
  JK, 8.8.2010
*/
void seqtop::coarse_local_map (gtopology *top,FILE *out)
{
  compute_multiplicity (top,out);
  node_local_numbers (out);
  node_global_glued_numbers (out);
  node_coarse_numbers (out);
  node_coarse_global_glued_map (out);
  node_coarse_local_map (out);
 
}

/**
   function assembles lists of elements belonging to subdomains
   
   list of elements of the i-th subdomain is needed e.g. in subcycling
   
   JK, 27. 5. 2014
*/
void seqtop::elem_lists ()
{
  long i,domid;
  
  //  array containing the number of elements on subdomains
  if (ned!=NULL)
    delete [] ned;
  ned = new long [ns];
  for (i=0;i<ns;i++){
    ned[i]=0;
  }
  
  //  loop over the number of elements
  for (i=0;i<tnep;i++){
    ned[eldom[i]]++;
  }//  end of the loop over the number of elements
  
  
  //  domain-element correspondence
  if (domel!=NULL){
    for (i=0;i<ns;i++){
      delete [] domel[i];
    }
    delete domel;
  }
  
  domel = new long* [ns];
  for (i=0;i<ns;i++){
    domel[i] = new long [ned[i]];
    ned[i]=0;
  }
  
  //  loop over the number of elements
  for (i=0;i<tnep;i++){
    //  domain id
    domid=eldom[i];
    domel[domid][ned[domid]]=i;
    ned[domid]++;
  }//  end of the loop over the number of elements
  
}
