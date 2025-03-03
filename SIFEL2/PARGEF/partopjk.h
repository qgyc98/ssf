#ifndef PARTOPJK_H
#define PARTOPJK_H

#include <stdio.h>
#include "../GEFEL/galias.h"
#include "../GEFEL/gtopology.h"

class psolver;

/**
   JK, huge revision 22.8.2011
*/
class partopjk
{
 public:
  ///  constructor
  partopjk (int np,int mr,meshdescription meshd,char *nameproc, int nl);
  ///  destructor
  ~partopjk ();
  
  ///  function assembles the array nnsd which contains the numbers of nodes on subdomains
  void numbers_of_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);

  ///  function assembles the array mltg on the master which stores all ltg arrays from processors
  void ltg_on_master (long *ltg,long *domproc,FILE *out);
  
  ///  function converts any mesh description to the all_nodes mesh description
  void ltg_conversion (long *ltg,long *domproc,FILE *out,char *proc_name);
  
  ///  function assembles array nodmultip which contains nodal multiplicity
  void node_multiplicity (FILE *out);
  
  ///  function assembles array noddom
  void node_domain (FILE *out);

  ///  function assembles array ndofnm which contains the number of DOFs for each node
  void ndofn_on_master (gtopology *top,long *domproc,FILE *out);
  
  ///  function assembles array dofm which contains DOF indicators on each node
  void dof_on_master (gtopology *top,long *domproc,FILE *out);
  
  ///  function detects coupled DOFs in the problem
  void coupled_dofs_detection (FILE *out);
  
  ///  function detects internal and boundary/interface nodes
  void internal_boundary_nodes_detection (FILE *out);

  /// function detects boundary/interface nodes with some constraints
  void boundary_nodes_with_constraints_detection (FILE *out);

  ///  function assembles ordering of unknowns on subdomains
  long schur_local_ordering (gtopology *top,FILE *out);

  ///  function assembles ordering of the coarse problem with respect to the Schur complement method
  void schur_coarse_ordering (FILE *out);
  
  ///  function assembles coarse code numbers for particular subdomains
  void coarse_code_numbers (FILE *out);

  
  //  to se pak smaze
  long vse (long *ltg,gtopology *top,long *domproc,FILE *out,char *proc_name);
  
  
  ///  FETI
  ///  function assembles array fullnoddom
  void full_node_domain (FILE *out);
  
  ///  FETI
  ///  function generates Lagrange multipliers
  void lagrange_multipliers (psolver ps,FILE *out);
  
  ///  FETI
  ///  function assembles coarse code numbers for particular subdomains
  void coarse_code_numbers_feti (psolver ps,FILE *out);
  
  
  
  
  ///  the number of processors
  ///  determined in the constructor
  int nproc;
  ///  rank of processor (processor id)
  ///  determined in the constructor
  int myrank;
  ///  mesh description
  ///  md = 1 - all nodes have their global node number
  ///  md = 2 - only boundary nodes have coarse number, internal nodes are denoted by -1
  ///  md = 3 - all nodes have their global node number, boundary nodes have their global number multiplied by -1
  ///  determined in the constructor
  meshdescription md;
  ///  name of processor
  char procName[10000];
  //char procName[MPI_MAX_PROCESSOR_NAME];
  ///  length of the processor name
  int nameLength;
  
  ///  the number of nodes on subdomain
  ///  determined in the function numbers_of_nodes_on_subdomains
  long nn;
  ///  the maximum number of nodes on subdomain
  ///  determined in the function numbers_of_nodes_on_subdomains
  long maxnn;
  
  ///  the maximum number of DOFs on subdomains
  ///  at the beginning, it contains the maximum number of possible DOFs
  ///  later on, it may be reduced due to coupled DOFs, constraints, etc.
  ///  it is determined in the function ndofn_on_master
  long maxndofd;
  
  ///  the number of internal nodes
  ///  variable is defined in function internal_boundary_nodes_detection
  long nin;

  ///  the number of boundary/interface nodes on subdomain
  ///  variable is determined in function internal_boundary_nodes_detection
  long nbn;
  
  ///  the number of internal DOFs on subdomain
  ///  the variable is determined in the function schur_local_ordering
  long nidof;
  
  ///  the number of boundary/interface DOFs on subdomain
  ///  the variable is determined in the function schur_local_ordering
  long nbdof;
  
  ///  the number of DOFs on subdomain
  ///  the variable is determined in the function schur_local_ordering
  long ndof;
  
  ///  the maximum number of boundary/interface DOFs on subdomain
  ///  the variable is determined in the function coarse_code_numbers
  long maxnbdofd;
  
  
  ///  array containing boundary/interface nodes on subdomains
  ///  it contains nbn components
  ///  bnid[i]=j - the i-th boundary/interface node is the j-th in gtopology list
  ///  array is allocated in function internal_boundary_nodes_detection
  long *bnid;

  ///  array containing internal nodes on subdomains
  ///  it contains nin components
  ///  inid[i]=j - the i-th internal node is the j-th in gtopology list
  ///  array is allocated in function internal_boundary_nodes_detection
  long *inid;
  
  ///  array containing coarse code numbers of boundary/interface DOFs on subdomain
  ///  it contains nbdof components
  ///  bdof[i]=j - the i-th boundary/interface DOF has the coarse code number j
  ///  it is allocated and assembled in the function coarse_code_numbers
  long *bdof;
   
  // *************************************************************
  //  variables and arrays located only on the master processor
  // *************************************************************
  
  ///  the total number of nodes in the whole problem
  ///  it is determined in the function ltg_conversion
  long tnnp;
  
  ///  the number of DOFs in the coarse problem
  ///  it is determined in the function schur_coarse_ordering for the Schur complement method
  ///  it is determined in the function lagrange_multipliers for the FETI method
  long ndofc;

  ///  array containing the numbers of nodes on subdomains
  ///  it contains nproc components
  ///  nnsd[i]=j - the i-th subdomain contains j nodes
  ///  it is allocated and assembled in the function numbers_of_nodes_on_subdomains
  long *nnsd;
    
  ///  array of local to global correspondence on the master
  ///  it contains nproc rows and nnsd[i] columns
  ///  mltg[i][j]=k - the j-th node on the i-th subdomain has global number k
  ///  it is allocated and assembled in the function ltg_on_master
  long **mltg;
  
  ///  array of nodal multiplicity, the multiplicity is the number of
  ///  subdomains which share the node
  ///  it contains tnnp components
  ///  nodmultip[i]=j - the i-th node (in global ordering) is shared by j subdomains
  ///  it is assembled in the function node_multiplicity
  long *nodmultip;
  
  ///  array of node-domain correspondence
  ///  it contains tnnp components
  ///  noddom[i]=j - the i-th node belongs to the j-th subdomain
  ///  it is assembled in the function node_domain
  ///  it is used for detection of coupled DOFs in the Schur complement method
  ///  only nodes on different subdomains are taken into account, therefore it
  ///  is enough to allocate only one number for each node
  long *noddom;

  ///  array of numbers of DOFs on nodes
  ///  it contains tnnp components
  ///  ndofnm[i]=j - the i-th node contains j DOFs
  ///  it is assembled in the function ndofn_on_master
  long *ndofnm;
  
  ///  array of DOF indicators on nodes
  ///  it contains tnnp rows and ndofnm[i] columns
  ///  dofm[i][j]=k - the j-th DOF on the i-th node has id equal to k
  ///  it is assembled in the function dof_on_master
  long **dofm;

  ///  maximum number of boundary/interface nodes with constraints collected from all subdomains
  ///  needed in the computation of norm of quantity flux resultant (compute_quantfluxres_norm function)  
  ///  it is assembled in the function boundary_nodes_with_constraints_detection
  long maxnbnwcd;

  ///  maximum number of DOFs at boundary/interface nodes with constraints collected from all subdomains
  ///  needed in the computation of norm of quantity flux resultant (compute_quantfluxres_norm function)  
  ///  it is assembled in the function boundary_nodes_with_constraints_detection
  long maxndofbn;

  /// number of boundary/interface nodes with constraints on the given subdomain, i.e. dimension of array bnwc
  /// set in function boundary_nodes_with_constraints_detection
  long nbnwc;

  /// array of local index numbers of boundary/interface nodes with constraints on the given subdomain
  /// bnwc[i] = local index number of i-th boundary/interface node with constraints on the given subdomain
  /// allocated on each subdomain in function boundary_nodes_with_constraints_detection
  long *bnwc;

  /// array of number of boundary/interface nodes on individual subdomains
  /// nbnwcd[i] = number of boundary/interface nodes with constraints on i-th subdomain
  /// allocated on master in function boundary_nodes_with_constraints_detection
  long *nbnwcd;

  /// array of numbers of boundary/interface nodes with constraints on individual subdomains
  /// nbnwcd[i][j] = local index number of j-th boundary/interface node with constraints on i-th subdomain
  /// allocated on master in function boundary_nodes_with_constraints_detection
  long **bnwcd;

  /// maximum number of boundary/interface nodes on subdomains, i.e. maximum value from the array nbnd
  ///  it is assembled in the function boundary_nodes_with_constraints_detection
  long  maxnbnd;

  ///  array containing the numbers of boundary/interface nodes on subdomains
  ///  it contains nproc components
  ///  nbnd[i]=j - the i-th domain contains j boundary/interface nodes
  ///  array is allocated in function internal_boundary_nodes_detection
  long *nbnd;

  ///  array containing the numbers of internal nodes on subdomains
  ///  it contains nproc components
  ///  nind[i]=j - the i-th domain contains j internal nodes
  ///  array is allocated in function internal_boundary_nodes_detection
  long *nind;
  
  ///  array containing the numbers of boundary/interface DOFs on subdomains
  ///  it contains nproc components
  ///  nbdofd[i]=j - the i-th subdomain contains j boundary/interface DOFs
  ///  array is allocated in the function coarse_code_numbers
  long *nbdofd;
  
  ///  array containing coarse code numbers of boundary/interface DOFs on subdomains
  ///  it contains nproc rows and nbdofd[i] columns
  ///  bdofd[i][j]=k - the j-th boundary/interface DOF on the i-th subdomain has the coarse code number k
  ///  it is allocated and assembled in the function coarse_code_numbers
  long **bdofd;
  
  
  ///  FETI
  ///  array of full node-domain correspondence
  ///  it contains tnnp rows and nodmultip[i] columns
  ///  fullnoddom[i][j]=k - the j-th subdomain which shares the i-th node has number (id) k
  ///  it is assembled in the function full_node_domain
  long **fullnoddom;
  
  ///  FETI
  ///  array of Lagrange multipliers, it contains coarse code numbers
  ///  it contains tnnp rows and nc columns
  ///  nc is either (nodmultip[i]-1)*ndofnm[i] or nodmultip[i]*(nodmultip[i]-1)/2*ndofnm[i]
  ///  lagrmultip[i][j]=k - the j-th Lagrange multiplier in the i-th node has coarse number k
  long **lagrmultip;

};


#endif
