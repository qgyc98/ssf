#ifndef PARTOP_H
#define PARTOP_H

#include "pgalias.h"
#include "../GEFEL/galias.h"
#include "../GEFEL/intp.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include <stdio.h>

/**
   class deals with topology in parallel computations
   it is superior class to the class gtopology which is used
   for particular subdomains
   the class partop manages several classes gtopology
   
   
   notation:

   node multiplicity - number of subdomains which share the node

   internal node - node with multiplicity 1, (node inside a subdomain,
                   node not lying on inter-subdomain boundary)

   boundary node - node with multiplicity at least 2, (node lying
                   on inter-subdomain boundary)

   coarse node - artificial node defined on master which collects appropriate
                 boundary nodes

   
   global ordering - node ordering before mesh partitioning
   global node numbers - node numbers used in undecomposed mesh
                        (node numbers used before decomposition/partitioning)

   local ordering - ordering of subdomains (without respect to remaining subdomains)
   local node numbers - node numbers used in decomposed mesh
                        (node numbers used after decomposition/partitioning)

   coarse ordering - ordering of boundary (interface) nodes only
   coarse node number - node numbers of boundary nodes only, ordering
                        of all boundary nodes in coarse problem

   group ordering - ordering of selected nodes only
   group node number - only selected nodes are ordered on the whole problem, it is similar to
                       the coarse ordering


   mesh description:
   md = 1 - all nodes have their global node number
   md = 2 - only boundary nodes have coarse number, internal nodes are denoted by -1
   md = 3 - all nodes have their global node number, boundary nodes have their global number multiplied by -1
   
   JK
*/
class partop  
{
 public:
  partop(int np,int mr,long nd,meshdescription meshd,char *nameproc, int nl);
  ~partop ();
  void initiation (gtopology *top,long *ltg);
  
  
  
  
  //  function computes nodal multiplicity
  void compute_multiplicity (long *ltg,long *domproc,FILE *out);
  //  tnnp is established
  //  multip is assmebled
  //  nodmultip is assembled
  //  allnodes is assembled
  //  nbnd is assembled
  //  maxnbn is established
  //  tnbn is established


  void dof_multiplicity (gtopology *top,FILE *out);
  //  dofmultip is assembled

  
  void rewrite_ltg (long *ltg);

  
  void find_boundary_nodes (long *ltg,long *domproc,FILE *out);
  //  lgnbn is assembled
  //  lnbn is assembled
  //  nin is established
  //  lnin is assembled

  
  void boundary_nodes_on_master (gtopology *top,long *domproc,FILE *out);
  //  cnbn is assembled
  
  
  void coarse_local_nodes ();
  void sort_nodes (FILE *out);


  void number_of_bdofs_on_nodes (gtopology *top,long *domproc,FILE *out);
  //  nbdofnd is assembled
  
  
  void code_numbers_on_master (gtopology *top,long *domproc,FILE *out);
  //  maxnbdof is established
  //  lbcn is assembled

  
  // *******************************************
  // *******************************************

  void nodesplit_allnodes (gtopology *top,long *domproc);
  void nodesplit_bnodes (gtopology *top,long *domproc,FILE *out);
  void gcnprep (gtopology *top,long *domproc);
  void corner_nodes (gtopology *top, FILE *out);
  
  void dof_incidence (gtopology *top,long nbdof);
  void corner_nodes2d (gtopology *top,FILE *out);
  void corner_nodes3d (gtopology *top,FILE *out);

 



  // *************************************************
  //  VARIABLES AND ARRAYS DEFINED ON ALL PROCESSORS
  // *************************************************

  /// name of processor
  char procName[10000];
  //char procName[MPI_MAX_PROCESSOR_NAME];
   ///
  int nameLength;
  ///  number of domain attached to the processor
  long ndom;
  ///  number of elements on subdomain
  long ne;
  



  ///  maximum number of boundary DOFs on one subdomain
  ///  variable is defined in function code_numbers_on_master
  long maxnbdof;


    
  

  
  ///  array containing number of subdomains which share the degree of freedom
  ///  it contains ndof components
  ///  nodmultip[i]=j - the i-th dof belongs to j subdomains
  ///  array is allocated in function dof_multiplicity
  long *dofmultip;
  
  ///  array containing numbers of boundary nodes
  ///  it contains nbn components
  ///  lnbn[i]=j - the i-th boundary node has local number j
  ///  array is allocated in function find_boundary_nodes
  long *lnbn;

  ///  array containing numbers of internal nodes
  ///  it contains nin components
  ///  lnin[i]=j - the i-th internal node has local number j
  ///  array is allocated in function find_boundary_nodes
  long *lnin;

  ///  array containing list of global numbers of boundary nodes
  ///  lgnbn[i]=j - the i-th boundary node has global number k
  long *lgnbn;


  

  // ********************************************************
  //  VARIABLES AND ARRAYS DEFINED ONLY ON MASTER PROCESSOR
  // ********************************************************
  

  /// total number of DOFs in the whole problem
  long tndof;
  
  
  ///  number of DOFs in coarse (reduced) problem
  long ndofcp;

  


  ///  array containing number of subdomains which each boundary node belongs to
  ///  it contains tnnp components
  ///  multip[i]=j - the i-th node belongs to j subdomains
  ///  array is allocated in function compute_multiplicity
  long *multip;
  


  ///  array containing correspondence between global and coarse node numbers
  ///  it contains tnnp components
  ///  gcnbn[i]=j - the i-th node is the j-th boundary node, j=-1 for internal nodes
  ///  array is allocated and deallocated in function boundary_nodes_on_master
  long *gcnbn;
  
  ///  array containing coarse node numbers of boundary nodes
  ///  it contains nproc rows and nbnd[i] columns
  ///  cnbn[i][j]=k - the j-th boundary node on the i-th subdomain has coarse node number k
  ///  array is allocated in function boundary_nodes_on_master
  long **cnbn;

  ///  array containing numbers of DOFs on boundary nodes
  ///  it contains nproc rows and nbnd[i] columns
  ///  nbdofnd[i][j]=k - the j-th boundary node on the i-th subdomain contains k degrees of freedom
  ///  array is allocated in function number_of_bdofs_on_nodes
  long **nbdofnd;


  
  ///  array containing numbers of boundary DOFs on subdomains
  ///  it contains nproc components
  ///  ndofd[i]=j - the i-th domain contains j degrees of freedom on boundary
  ///  array is allocated in function code_numbers_on_master
  long *nbdofd;
  
  ///  array containing local boundary code numbers
  ///  lbcn[i][j][k]=m - the k-th DOF on the j-th boundary node on the i-th subdomain has number m
  ///  array is allocated in function code_numbers_on_master
  long ***lbcn;
  
  ///  array containing multiplicity of boundary nodes
  ///  it is different array from array multip which is defined for all nodes in problem
  ///  it contains tnbn components
  ///  bnmultip[i]=j - the i-th boundary node has multiplicity j
  ///                  the i-th boundary node is shared by j subdomains
  ///  array is allocated in function coarse_local_nodes
  long *bnmultip;
  
  ///  array containing local numbers of boundary nodes belonging to coarse nodes
  ///  it contains tnbn rows and bnmultip[i] columns
  ///  llnbn[i][j]=k - the j-th node belonging to the i-th coarse node has number k
  ///  array is allocated in function
  long **llnbn;
  
  ///  array containing numbers of subdomain to which coarse nodes belong to
  ///  it contains tnbn rows and bnmultip[i] columns
  ///  ldn[i][j]=k - the j-th node belonging to the i-th coarse node belongs to the k-th domain
  ///  array is allocated in function
  long **ldn;
  
  ///  array containing numbers of DOFs at coarse nodes
  ///  it contains tnbn components
  ///  ndofncn[i]=j - the i-th coarse node contains j DOFs
  ///  array is allocated in function schurordering
    //long *ndofncn;
  
  ///  array containing code numbers at coarse nodes
  ///  it contains tnbn rows and ndofncn[i] columns
  ///  cncn[i][j]=k - the j-th unknown at the i-th coarse node has number k
  ///  array is allocated in function schurordering
    //long **cncn;


  
  
  
  
    
  
  ///  array containing global code numbers
  ///  gcn[i][j]=k - the j-th component of the i-th boundary node (the i is the global number of node) is equal to the k
  long **gcn;
  


/*
  ZACATEK ROZSAHLYCH UPRAV
  27.7.2009
*/

  //  function collects numbers of all nodes on subdomains
  void numbers_of_all_nodes_on_subdomains (long *domproc,FILE *out);
  //  maxnn is established
  //  nnsd is assembled

  //  function assembles node multiplicity
  void assemble_multip (long *ltg,long *domproc,FILE *out,char *proc_name);
  //  array bmultip or amultip is assembled
  //  tnbn is determined
  //  tnnp is determined if appropriate node ordering is used
  //  maxnbn is determined if appropriate node ordering is used
  //  allnodes is allocated and assembled if appropriate node ordering is used


  void coupled_dofs (gtopology *top,long *domproc,FILE *out);
  //  ncdof is determined
  //  coupdof is allocated and assembled
  //  coupdofmas is allocated and assembled
  
  //  function collects the numbers of all degrees of freedom on subdomains
  //  supports and constraints are added as DOFs
  void numbers_of_all_dofs_on_subdomains (gtopology *top,long *domproc,FILE *out);
  //  maxndof is established
  //  nalldof is assembled
  
  //  function assembles array containing the numbers of DOFs on nodes
  void ndofn_master (gtopology *top,long *domproc,FILE *out,char *proc_name);
  //  ndofnmas is allocated and assembled
  
  //  function collects id of DOFs on the master
  void dofind_master (gtopology *top,long *domproc,FILE *out,char *proc_name);
  //  dofindmas is allocated and assembled
  
  
  //  function updates node multiplicity
  void update_multip (gtopology *top,FILE *out,char *proc_name);
  //  dofind is allocated and assembled
  
  
  //  function assembles the arrays nbnd a nind
  void assemble_nbnd_nind (long *ltg,long *domproc,FILE *out,char *proc_name);
  //  nbnd is allocated and assembled
  //  nind is allocated and assembled


  //  function assembles the array nodmultip
  void assemble_nodmultip (long *ltg,long *domproc,FILE *out,char *proc_name);
  //  nodmultip is allocated and assembled
  
  
  //  function assembles global numbers of nodes
  void node_global_numbers (long *domproc,FILE *out,char *proc_name);
  //  gnin is allocated and assembled
  //  gnbn is allocated and assembled
  
  //  function assembles local numbers of nodes
  void node_local_numbers (FILE *out);
  //  lnindom is allocated and assembled
  //  lnbndom is allocated and assembled

  //  function assembles coarse numbers of boundary/interface nodes
  void node_coarse_numbers (long *domproc,FILE *out,char *proc_name);
  //  icnbnmas is allocated and assembled
  //  icmultip is allocated and assembled
  
  
  //  function assembles coarse - global numbers map
  void node_coarse_global_map (FILE *out,char *proc_name);
  //  gnbncn is allocated and assembled
  //  sid is allocated and assembled

  //  function generates the Schur complement ordering on subdomains
  //  it takes into account the coupled DOFs
  long schur_ordering (gtopology *top,long *domproc,FILE *out,char *proc_name);
  
  
  void prepare_data (long *domproc,long *ltg,gtopology *top,FILE *out,char *proc_name);
  
  //function assembles global code numbers on master processor 
  // for using of class parcongrad
  void  assemble_gcnd(gtopology *top,long *domproc,FILE *out);
  //  gcnd is assembled

  
  ///  the number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  mesh description
  meshdescription md;
  ///  the number of nodes on subdomain
  long nn;
  ///  the number of DOFs on subdomain
  long ndof;
  ///  the number of internal DOFs on subdomain
  long nidof;
  ///  the number of boundary/interface DOFs on subdomain
  long nbdof;
  


  ///  the number of coupled DOFs
  ///  if two or more DOFs are coupled, it means that they have the same
  ///  code number, they will be localized into one position in vectors and matrices
  ///  such DOFs are aggregated into one DOF
  ///  therefore, ncdof is the number of aggregates
  ///  it is determined in the function coupled_dofs
  long ncdof;

  ///  the number of boundary/interface coupled DOFs
  ///  it is determined in the function coupled_dofs
  long nbcdof;

  ///  the maximum number of nodes on one subdomain
  ///  variable is determined in function numbers_of_all_nodes_on_subdomains
  long maxnn;

  //  master
  ///  array containing the numbers of nodes on subdomains
  ///  it contains nproc components
  ///  nnsd[i]=j - j nodes are defined on the i-th subdomain
  ///  array is allocated in function numbers_of_all_nodes_on_subdomains
  long *nnsd;

  //  master
  ///  array containing number of subdomains which each boundary/interface node belongs to
  ///  it contains tnbn components,
  ///  multip[i]=j - the i-th interface/boundary node belongs to j subdomains
  ///  array is assembled in the function assemble_multip (FILE *out)
  long *bmultip;

  //  master
  ///  array containing number of subdomains which each node belongs to
  ///  it contains tnnp components,
  ///  multip[i]=j - the i-th node belongs to j subdomains
  ///  array is assembled in the function assemble_multip (FILE *out)
  long *amultip;

  ///  the number of internal nodes
  ///  variable is defined in function find_boundary_nodes
  long nin;

  ///  the number of boundary/interface nodes on subdomain
  ///  variable is determined in function assemble_multip
  long nbn;

  ///  the maximum number of boundary/interface nodes on one subdomain
  ///  variable is defined in function assemble_multip
  long maxnbn;

  ///  the total number of boundary/interface nodes
  ///  variable is defined in function assemble_multip
  long tnbn;

  ///  the total number of internal nodes
  ///  variable is defined in function assemble_multip
  long tnin;

  //  master
  ///  the total number of nodes in the whole problem
  ///  variable is defined in function compute_multiplicity
  long tnnp;

  //  master
  ///  array containing the numbers of boundary/interface nodes on subdomains
  ///  it contains nproc components
  ///  nbnd[i]=j - the i-th domain contains j boundary/interface nodes
  ///  array is allocated in function assemble_nbnd_nind
  long *nbnd;

  ///  array containing the numbers of internal nodes on subdomains
  ///  it contains nproc components
  ///  nind[i]=j - the i-th domain contains j internal nodes
  ///  array is allocated in function assemble_nbnd_nind
  long *nind;

  ///  array containing numbers of subdomains which share the nodes
  ///  it contains nn components
  ///  nodmultip[i]=j - the i-th node belongs to j subdomains
  ///  array is allocated in function assemble_multip
  long *nodmultip;

  //  master
  ///  array containing global node numbers
  ///  this array is assembled for mesh description = all_nodes
  ///  it contains nproc rows and nn[i] columns
  ///  allnodes[i][j]=k - the j-th node on the i-th subdomain has global number k
  ///  array is allocated in function assemble_multip
  long **allnodes;

  //  master
  ///  array containing coarse node numbers
  ///  this array is assembled for mesh description = bound_nodes
  ///  it contains nproc rows and nn[i] columns
  ///  bnodes[i][j]=k - the j-th node on the i-th subdomain has coarse number k
  ///  if k=-1, the node is internal
  ///  array is allocated in function assemble_multip
  long **bnodes;

  //  master
  ///  array containing indicators of coupled DOFs
  ///  it has ns, ncdof components
  ///  coupdof[i][j]=k - the j-th coupled DOFs on the i-th subdomain is shared by k nodes
  ///  array is allocated in the function coupled_dofs
  long **coupdof;

  //  master
  ///  array containing suspicious indicators of coupled DOFs
  ///  it has ncdof components
  ///  coupdofmas[i]=0 - the i-th coupled DOF is not a boundary/interface DOF
  ///  coupdofmas[i]=1 - the i-th coupled DOF is a boundary/interface DOF
  ///  array is allocated in the function coupled_dofs
  long *coupdofmas;

  //  master
  ///  array containing numbers of all degrees of freedom on subdomains
  ///  nalldof[i]=j - the i-th subdomain contains j degrees of freedom
  ///  array is allocated in function numbers_of_all_dofs_on_subdomains
  long *nalldof;

  ///  the maximum number of all DOFs on one subdomain
  ///  variable is determined in function numbers_of_all_dofs_on_subdomains
  long maxndof;
  
  //  master
  ///  array containing the numbers of DOF on nodes
  ///  it contains tnnp components
  ///  ndofnmas[i]=j - the i-th node (in the global ordering) contains j DOFs
  ///  array is allocated and assembled in the function ndofn_master
  long *ndofnmas;

  ///  array containing DOF indicators
  ///  indicators are used for code number generation
  ///  it has tnnp, ndofnmas[i] components
  ///  dofindmas[i][j]=k - the j-th DOF in the i-th node has indicator k
  ///  array is allocated in the function dofind_master
  long **dofindmas;

  ///  array containing DOF indicators
  ///  if there are coupled DOFs, it is not enough to deal with nodes
  ///  DOFs have to be split to internal and boundary/interface
  ///  it has tnnp, ndofnmas[i] components
  ///  dofind[i][j]=0 - the j-th DOF in the i-th node is internal
  ///  dofind[i][j]=1 - the j-th DOF in the i-th node is boundary/interface
  ///  array is allocated in the function update_multip
  long **dofind;

  ///  array containing coarse numbers of interface/boundary nodes
  ///  it contains ns, nbnd[i] components
  ///  icnbnmas[i][j]=k - the j-th boundary node on the i-th subdomain has coarse number k
  ///  array is allocated in the function node_coarse_numbers
  long **icnbnmas;

  ///  array containing coarse numbers of interface/boundary nodes on one subdomain
  ///  it contains nbn components
  ///  icnbn[i]=j - the i-th boundary node has coarse number j
  ///  array is allocated in the function node_coarse_numbers
  long *icnbn;
  
  ///  array containing node multiplicity of boundary/interface nodes
  ///  it contains tnbn components
  ///  icmultip[i]=j - the i-th boundary/interface node belongs to j subdomains
  ///   array is allocated in the function node_coarse_numbers
  long *icmultip;

  ///  global numbers of boundary/interface nodes on the master
  ///  it contains nproc, nbn components
  ///  gnbn[i][j]=k - the j-th boundary/interface node on the i-th subdomain has global number k
  ///  array is allocated and assembled in the function node_global_numbers
  long **gnbn;

  ///  global numbers of internal nodes on the master
  ///  it contains nproc, nin components
  ///  gnin[i][j]=k - the j-th internal node on the i-th subdomain has global number k
  ///  array is allocated and assembled in the function node_global_numbers
  long **gnin;

  ///  global numbers of boundary/interface nodes on subdomains
  ///  it contains nbn components
  ///  gnbndom[i]=j - the i-th boundary/interface node has global number j
  ///  array is allocated and assembled in the function node_global_numbers
  long *gnbndom;

  ///  local numbers of boundary/interface nodes on subdomains
  ///  it contains nbn components
  ///  lnbndom[i]=j - the i-th boundary/interface node has local number j
  ///  array is allocated and assembled in the function node_local_numbers
  long *lnbndom;

  ///  local numbers of internal nodes on subdomains
  ///  it contains nin components
  ///  lnindom[i]=j - the i-th internal node has local number j
  ///  array is allocated and assembled in the function node_local_numbers
  long *lnindom;

  ///  global numbers of boundary/interface nodes appropriate to coarse node
  ///  it contains global numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  gnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global number k
  ///  aray is allocated in the function node_coarse_global_map
  long **gnbncn;

  ///  subdomain id of interface/boundary nodes appropriate to coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  sid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the function node_coarse_local_map or node_coarse_global_map
  long **sid;
  
  /// master
  ///  nud[i]=j - the i-th domain contains j unknowns (degrees of freedom)
  long *nud;
  
  /// master
  ///  gcnd[i][j]=k - the j-th unknown on the i-th subdomain has global number k
  ///  gcnd[i] is a pointer to array with code numbers of the i-th subdomain
  ///  it is similar to arrays cn on elements, the subdomain plays a role of a superelement
  long **gcnd;



  
   
  
 

  // backup of ltg from psolver (used in Temelin)
  long *bltg;



  

  



};


#endif
