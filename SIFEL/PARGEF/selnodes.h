#ifndef SELNODES_H
#define SELNODES_H

#include <stdio.h>
#include "pgalias.h"
#include "../GEFEL/galias.h"
#include "../GEFEL/gtopology.h"

/**
   class deals with selected nodes on subdomains
   nodes could be selected for Schur complement method,
   FETI method, etc.
   
*/
class selnodes
{
 public:
  selnodes (int np,int mr,long nd,long kk,long *j,meshdescription d,FILE *out,long mespr);
  ~selnodes ();
  
  //  
  void number_of_selected_nodes (long *domproc,FILE *out);
  //  maxnsn is obtained
  //  nsndom is assembled
  
  void nodes_on_master (long *domproc,FILE *out);
  
  void node_multiplicity (FILE *out);
  void dof_multiplicity (FILE *out);


  // *****************
  //  SCHUR ORDERING
  // *****************
  //void schur_ordering (gtopology *top,FILE *out);

  
  // ****************
  //  FETI ORDERING
  // ****************
  void group_local_nodes (FILE *out);
  //  ljn is assmebled
  //  lsn is assembled


  void dof_feti (FILE *out);
  //  doffeti is assembled

  void number_contrib (long *domproc,FILE *out);
  void contrib_dofs (gtopology *top,long *domproc,FILE *out);
    
    
  // *************************************************
  //  VARIABLES AND ARRAYS DEFINED ON ALL PROCESSORS
  // *************************************************
  
  ///  number of domain attached to the processor
  long ndom;


  long ndof;
  
  
  ///  number of nodes which contribute to coarse problem in the FETI method
  ///  there are some nodes in the FETI method, which have to be used several times
  ///  because they are connected with more than one additional node
  ///  therefore, ncn is generally not equal to nsn
  long ncn;

  ///  maximum number of contributing nodes on subdomain in the FETI method
  ///  there are some nodes in the FETI method, which have to be used several times
  ///  because they are connected with more than one additional node
  ///  therefore, maxncn is generally not equal to maxnsn
  long maxncn;

  

  ///  list of selected nodes - global numbers
  ///  lsng[i]=j - the i-th selected node has global/coarse number j
  ///  lsng contains nsn components, where nsn is the number of selected nodes on subdomain
  long *lsng;
  
  
  ///  list of code numbers which contribute to the coarse problem
  long *ldof;
  
  long *ldofmultip;
  long nldofmultip;
  
  // ********************************************************
  //  VARIABLES AND ARRAYS DEFINED ONLY ON MASTER PROCESSOR
  // ********************************************************
  
  ///  total number of DOFs
  long tndof;

  ///  the maximum number of selected DOFs on subdomain
  long maxndof;

  
  ///  group node numbers (see partop.h)
  ///  gnn[i][j]=k - the j-th selected node on the i-th subdomain has group number k
  long **gnn;
  
  ///  node multiplicity
  ///  nodmultip[i]=j - the i-th selected node is shared by j subdomains
  ///  it contains tnsn components
  long *nodmultip;

  /// dof multiplicity
  /// dofmultip[i]=j - the i-th selected dof is shared by j subdomains
  /// it contains components  
  long *dofmultip;
  

  ///  array of numbers of DOFs on subdomains 
  ///  it contains prescribed values before code number generation
  ///  after code number generation, only unknown (unconstrained DOFs) are taken into account
  ///  array is rewritten in the function schur_ordering
  ///  ndofdom[i]=j - the i-th subdomain contains j DOFs
  long *ndofdom;


  //  SCHUR ORDERING

  ///  code numbers at selected nodes on master
  ///  cnm[i][j]=k - the j-th DOF at the i-th selected node has code number / indicator k
  long **cnm;
  
  ///  code numbers on master
  ///  cndom[i][j]=k - the j-th DOF on the i-th subdomain has group code number k
  long **cndom;
  

  //  FETI ORDERING
  
  ///  list of joint nodes to selected nodes assumed as coarse nodes
  ///  ljn[i][j]=k - the j-th node connected to the i-th coarse node has local number k
  ///  ljn contains tnsn rows and nodmultip columns
  long **ljn;
  
  ///  list of subdomain numbers which contain connected nodes to coarse nodes
  ///  lsn[i][j]=k - the j-th node connected to the i-th coarse node belongs to the k-th subdomain
  ///  lsn contains tnsn rows and nodmultip columns
  long **lsn;
  
  ///  code numbers / indicators for FETI method
  ///  doffeti[i][j][k]=l - the k-th DOF on the j-th connected node to the i-th coarse node has code number / indicator l
  long ***doffeti;
  
  ///  number of contributing nodes in the FETI method
  ///  ncndom[i]=j - the i-th subdomain contributes to the coarse problem by j nodes
  long *ncndom;
  
  
  






















































/*
  ZACATEK ROZSAHLYCH UPRAV
  27.7.2009
*/
  
  //  constructor
  selnodes(int nd,int mr,meshdescription d,long *domproc,long *nnsd,long *jj,
	   long itnbn,long *iicmultip,long *nodmultip,long *icnbn,long *gnbndom,
	   long **ignbncn,long **isid,
	   FILE *out,char *proc_name,long mespr);
  
  //  function collects coarse numbers of nodes on the master
  void node_coarse_numbers (FILE *out);
  //  tnsn is determined
  //  snicmultip is allocated and assembled
  
  //  function searches for all possible DOFs in selected nodes
  void number_all_dofs (gtopology *top,long *domproc,FILE *out);
  //  snndof is determined
  //  maxsnndof is determined
  //  snndofdom is allocated and assembled
  
  //  function collects numbers of DOFs on selected nodes
  void ndofn_on_master (gtopology *top,long *domproc,FILE *out);
  //  snndofnmas is allocated and assembled


  //  function assembles DOFs indicators on master
  void dof_indicators (gtopology *top,long *domproc,FILE *out);
  //  sndofmas is allocated and assembled
  
  //  function generates ordering of selected nodes with respect
  //  to the Schur complement method
  void schur_ordering (gtopology *top,long **dofind,FILE *out);
  void schur_ordering (gtopology *top,FILE *out);
  
  
  
  ///  the number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  type of mesh description
  meshdescription md;

  ///  the number of nodes on subdomain
  long nn;
  
  ///  the number of selected nodes
  long nsn;

  ///  the maximum nuber of selected nodes on subdomain
  long maxnsn;

  ///  the total number of boundary/interafce nodes inthe problem
  long tnbn;
  
  ///  the total number of selected nodes
  long tnsn;
  
  ///  the number of all selected DOFs including prescribed DOFs on one subdomain
  long snndof;

  ///  the maximum number of selected DOFs on subdomain
  long maxsnndof;


  ///  number of multiplicity of all boundary/interface nodes
  ///  it contains tnbn components
  ///  icmultip[i]=j - the i-th boundary/interface node shares j nodes (it belongs to j subdomains)
  ///  array is allocated in the constructor
  long *icmultip;

  ///  array containing the numbers of selected nodes on subdomains
  ///  it contains nproc components
  ///  nsnmas[i]=j - j nodes are selected on the i-th subdomain
  ///  array is allocated and assembled in the constructor
  long *nsnmas;
  
  ///  local numbers of selected nodes
  ///  it contains nsn components
  ///  lnsn[i]=j - the i-th selected node has local number j
  ///  array is allocated in the constructor
  long *lnsn;

  ///  coarse numbers of selected nodes on one subdomain
  ///  cnsn contains nsn components
  ///  cnsn[i]=j - the i-th selected node has coarse number j
  ///  array is allocated in the constructor
  long *cnsn;

  //  master
  ///  coarse numbers of selected nodes on the master
  ///  it contains nproc,nsn components
  ///  cnsnmas[i][j]=k - the j-th selected node on the i-th subdomain has coarse number k
  ///  array is allocated in the constructor
  long **cnsnmas;

  ///  global numbers of selected nodes on one subdomain
  ///  it contains nsn components
  ///  gnsn[i]=j - the i-th selected node has global number j
  ///  array is allocated in the constructor
  long *gnsn;

  //  master
  ///  global numbers of selected nodes on the master
  ///  gnsnmas contains nproc, nsnmas[i] components
  ///  gnsnmas[i][j]=k - the j-th selected node on the i-th subdomain has coarse number k
  ///  array is allocated in the constructor
  long **gnsnmas;

  //  master
  ///  number of multiplicity of selected boundary/interface nodes
  ///  it contains tnsn components
  ///  snicmultip[i]=j - the i-th selected boundary/interface node shares j nodes (it belongs to j subdomains)
  ///  array is allocated in the function node_coarse_numbers
  long *snicmultip;
  
  //  master
  ///  array of numbers of DOFs on subdomains at selected nodes
  ///  it contains prescribed values before code number generation
  ///  after code number generation, only unknown (unconstrained DOFs) are taken into account
  ///  snndofmas[i]=j - selected nodes on the i-th subdomain contain j DOFs
  ///  array is allocated in the function number_all_dofs
  long *snndofmas;

  //  master
  ///  array of numbers of DOFs in selected nodes on the master
  ///  it contains nproc, nsndom[i] components
  ///  snndofnmas[i][j]=k - the j-th node on the i-th subdomain contains k DOFs
  ///  array is allocated in the function ndofn_on_master
  long **snndofnmas;

  //  master
  ///  array of DOFs or indicators at selected nodes
  ///  it contains nproc, nsndom[i], snndofndom[i][j] components
  ///  sndofmas[i][j][k]=l - the k-th DOF at the j-th selected node on the i-th subdomain has value l
  ///  array is allocated in the function dof_indicators
  long ***sndofmas;
  
  //  master
  ///  global numbers of boundary/interface nodes appropriate to all coarse nodes
  ///  it contains local numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  gnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global number k
  ///  aray is allocated in the constructor
  long **gnbncn;

  //  master
  ///  subdomain id of interface/boundary nodes appropriate to all coarse nodes
  ///  it contains tnbn rows and icmultip[i] columns
  ///  sid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the constructor
  long **sid;

  //  master
  ///  global numbers of selected boundary/interface nodes appropriate to coarse node
  ///  it contains tnsn rows and snicmultip[i] columns
  ///  sngnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global number k
  ///  aray is allocated in the function node_coarse_numbers
  long **sngnbncn;

  //  master
  ///  subdomain id of selected interface/boundary nodes of the coarse node
  ///  it contains tnsn rows and snicmultip[i] columns
  ///  snsid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the function node_coarse_numbers
  long **snsid;

  
  ///  array of numbers of DOFs for selected nodes
  ///  it contains tnsn components
  ///  ndofnsn[i]=j - the i-th selected node (in group ordering) has j DOFs
  ///  array is allocated in the function schur_ordering
  long *ndofnsn;

  //  master
  ///  code numbers at selected nodes on master
  ///  it contains tnsn,ndofnsn components
  ///  codensn[i][j]=k - the j-th DOF at the i-th selected node has code number / indicator k
  ///  array is allocated in the function schur_ordering
  long **codensn;

  ///  total number of DOFs on selected nodes
  ///  defined in the function schur_ordering
  long tndofsn;
  
  ///  code numbers on master
  ///  cnmas[i][j]=k - the j-th DOF on the i-th subdomain has coarse code number k
  long **cnmas;
  
};

#endif
