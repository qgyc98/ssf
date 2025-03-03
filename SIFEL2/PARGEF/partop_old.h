#ifndef PARTOP_H
#define PARTOP_H

#include <stdio.h>
#include "mpi.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include "pgalias.h"
#include "../GEFEL/galias.h"

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
  partop(int np,int mr,long nd,meshdescription meshd);
  ~partop ();
  void initiation (gtopology *top,long *ltg);
  
  //  function collects numbers of all nodes on subdomains
  void numbers_of_all_nodes_on_subdomains (long *domproc,FILE *out);
  //  maxnn is established
  //  nnsd is assembled
  
  
  //  function collects numbers of all degrees of freedom on subdomains
  //  supports and constraints are not taken into account
  void numbers_of_all_dofs_on_subdomains (gtopology *top,long *domproc,FILE *out);
  //  maxndof is established
  //  nalldof is assembled
  
  
  //  function computes nodal multiplicity
  void compute_multiplicity (long *ltg,long *domproc,FILE *out);
  //  tnnp is established
  //  multip is assmebled
  //  nodmultip is assembled
  //  allnodes is assembled
  

  void dof_multiplicity (gtopology *top,FILE *out);
  //  dofmultip is assembled

  
  void rewrite_ltg (long *ltg);

  
  void find_boundary_nodes (long *ltg,long *domproc,FILE *out);
  //  tnbn is established
  //  nbnd is assmebled
  //  maxnbn is established
  //  nbn is established
  //  lnbn is assembled
  //  lnin is assembled
  //  nin is established
  
  
  void boundary_nodes_on_master (gtopology *top,long *domproc,FILE *out);
  //  cnbn is assembled
  
  
  void coarse_local_nodes ();
  void sort_nodes (FILE *out);


  void number_of_bdofs_on_nodes (gtopology *top,long *domproc,FILE *out);
  //  nbdofnd is assembled
  
  
  void code_numbers_on_master (gtopology *top,long *domproc,FILE *out);
  //  maxnbdof is established
  //  lbcn is assembled
  long give_whole_dim(gtopology *top);
  long give_whole_elemtype(gtopology *top);
  void identification_node(gtopology *top,long *ltg,long *domproc, FILE *out);
  void identification_node_2d (gtopology *top,long *ltg,long *domproc, FILE *out);
  void identification_node_3d (gtopology *top,long *ltg,long *domproc, FILE *out);
  void identification_node_hexa (gtopology *top,long *ltg,long *domproc, FILE *out);
  void identification_node_tetra (gtopology *top, long *domproc, FILE *out);
  void corner_detection(gtopology *top, long *domproc, long *ltg, FILE *out);
  void control_numbers_of_vertices (gtopology *top,FILE *out);
  void control_vertices_on_master(gtopology *top, long *domproc, FILE *out);


  void  assemble_gnn(gtopology *top,long *domproc,FILE *out);
  //  gnn is assembled

  void  assemble_pgcn(gtopology *top,long *domproc,FILE *out);
  //  pgcn is assembled
  
  void  assemble_gcnd(gtopology *top,long *domproc,FILE *out);
  //  gcnd is assembled

  
  
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

  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of domain attached to the processor
  long ndom;
  ///  number of nodes on subdomain
  long nn;
  ///  number of elements on subdomain
  long ne;
  
  ///  maximum number of nodes on one subdomain
  ///  variables is defined in function numbers_of_all_nodes_on_subdomains
  long maxnn;

  ///  number of DOFs on subdomain
  long ndof;

  ///  maximum number of all DOFs on one subdomain
  ///  variables is defined in function numbers_of_all_dofs_on_subdomains
  long maxndof;

  ///  number of boundary nodes on subdomain
  ///  variable is defined in function find_boundary_nodes
  long nbn;

  ///  maximum number of boundary nodes on one subdomain
  ///  variable is defined in function find_boundary_nodes
  long maxnbn;

  ///  maximum number of boundary DOFs on one subdomain
  ///  variable is defined in function code_numbers_on_master
  long maxnbdof;

  ///  number of internal nodes
  ///  variable is defined in function find_boundary_nodes
  long nin;

  ///  mesh description
  meshdescription md;
  
    
  

  ///  array containing number of subdomains which share the nodes
  ///  it contains nn components
  ///  nodmultip[i]=j - the i-th node belong to j subdomains
  ///  array is allocated in function compute_multiplicity
  long *nodmultip;
  
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




  // array contaiting node identification on subdomain
  // it contains nn components
  // nodeidentif[i] = j for the i-th node on subdomain 
  // j = 1 - internal node
  // j = 2 - edge node
  // j = 3 - vertex node
  long *nodeidentif;
  



  // ********************************************************
  //  VARIABLES AND ARRAYS DEFINED ONLY ON MASTER PROCESSOR
  // ********************************************************
  
  ///  total number of nodes in the whole problem
  ///  variable is defined in function compute_multiplicity
  long tnnp;

  /// total number of DOFs in the whole problem
  long tndof;
  
  ///  total number of boundary nodes
  ///  variable is defined in function find_boundary_nodes
  long tnbn;
  
  ///  number of DOFs in coarse (reduced) problem
  long ndofcp;

  
  ///  array containing number of nodes on subdomain
  ///  it contains nproc components
  ///  nnsd[i]=j - j nodes are defined on the i-th subdomain
  ///  array is allocated in function numbers_of_all_nodes_on_subdomains
  long *nnsd;

  ///  array containing numbers of all degrees of freedom on subdomains
  ///  nalldof[i]=j - the i-th subdomain contains j number of degrees of freedom
  ///  array is allocated in function numbers_of_all_dofs_on_subdomains
  long *nalldof;

  ///  array containing number of subdomains which each boundary node belongs to
  ///  it contains tnnp components
  ///  multip[i]=j - the i-th node belongs to j subdomains
  ///  array is allocated in function compute_multiplicity
  long *multip;
  
  ///  array containing global node numbers
  ///  it contains nproc rows and nn[i] columns
  ///  allnodes[i][j]=k - the j-th node on the i-th subdomain has global number k
  ///  array is allocated in function compute_multiplicity
  long **allnodes;

  ///  array containing list of numbers of boundary nodes on subdomains
  ///  it contains nproc components
  ///  nbnd[i]=j - the i-th domain contains j boundary nodes
  ///  array is allocated in function find_boundary_nodes
  long *nbnd;

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


  ///  global numbers of boundary nodes on subdomains
  ///  it contains nbn components
  long *gnbn;
  
  
  
  
  ///  gnn[i]=j - unknows of the i-th node in the problem are located from the j-th position in tha array pgcn
  ///  gnn contains nn components, where nn is the number of all nodes in the problem
  long *gnn;
  
  ///  pgcn[gnn[i]+k]=j - the k-th unknown at the i-th node has number j
  ///  pgcn contains gndof components, where gndof is the sum of all possible DOFs, it means that
  ///  constraints are not taken into account
  long *pgcn;
  

  ///  nud[i]=j - the i-th domain contains j unknowns (degrees of freedom)
  long *nud;
  
  ///  gcnd[i][j]=k - the j-th unknown on the i-th subdomain has global number k
  ///  gcnd[i] is a pointer to array with code numbers of the i-th subdomain
  ///  it is similar to arrays cn on elements, the subdomain plays a role of a superelement
  long **gcnd;
  
  
  
  //  ZBYTKY

  
  
  
  ///  array containing global code numbers
  ///  gcn[i][j]=k - the j-th component of the i-th boundary node (the i is the global number of node) is equal to the k
  long **gcn;
  




  // *******************************************************************************
  // *******************************************************************************
  //  ctvrtek 30.6.2005
  // *******************************************************************************
  // *******************************************************************************

};

#endif
