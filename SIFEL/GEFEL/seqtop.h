#ifndef SEQTOP_H
#define SEQTOP_H

#include <stdio.h>
#include "galias.h"
#include "xfile.h"

class gtopology;

/**
   class imitates the class partop which deals with topology in parallel computations
   it is sequential version of the class partop and it is developed
   for simpler manipulation and possible application of debugger
   
   
   notation:

   node multiplicity - number of subdomains which share the node

   internal node - node with multiplicity 1, (node inside a subdomain,
		   node not lying on inter-subdomain boundary)

   interface node - node with multiplicity at least 2, (node lying
                    on inter-subdomain boundary)

   coarse node - artificial node defined on master which collects appropriate
                 interface nodes

   
   global ordering - node ordering before mesh partitioning
   global node numbers - node numbers used in undecomposed mesh
                        (node numbers used before decomposition/partitioning)

   local ordering - ordering of subdomains (without respect to remaining subdomains)
   local node numbers - node numbers used in decomposed mesh
                        (node numbers used after decomposition/partitioning)

   coarse ordering - ordering of interface nodes only
   coarse node number - node numbers of interface nodes only, ordering
                        of all interface nodes in coarse problem


   mesh description:
   md = bound_nodes = 1 - all nodes have their global node number
   md = all_nodes = 2 - only interface nodes have coarse number, internal nodes are denoted by -1
   md = neg_bound_nodes = 3 - all nodes have their global node number, interface nodes have their global number multiplied by -1
   
   JK, 11.9.2007
*/
class seqtop
{
 public:
  seqtop(long nd,meshdescription meshd);
  ~seqtop ();

  void read_nnsd (XFILE *in);
  //void read_nesd (XFILE *in);
  void read_ltg (XFILE *in);
  //void read_eltg (XFILE *in);
  void read_eldom (long ne,XFILE *in);

  void assemble_multip (FILE *out);
  void coupled_dofs (gtopology *top,FILE *out);
  void update_multip (gtopology *top,FILE *out);
  void assemble_nbnd_nind (FILE *out);
  void assemble_nodmultip (FILE *out);
  void assemble_dofind (gtopology *top,FILE *out);
  void compute_multiplicity (gtopology *top,FILE *out);

  //void compute_multiplicity (FILE *out);
  //void find_boundary_nodes (FILE *out);
  
  void node_local_numbers (FILE *out);
  void node_global_glued_numbers (FILE *out);
  void node_coarse_numbers (FILE *out);
  void node_coarse_local_map (FILE *out);
  void node_coarse_global_glued_map (FILE *out);
  
  long schur_ordering (gtopology *top,FILE *out);
  long schur_ordering_old (gtopology *top,FILE *out);

  void coarse_local_map (gtopology *top,FILE *out);

  //  function assembles lists of elements belonging to subdomains
  void elem_lists ();
  
  
  //void rewrite_ltg ();
  //void codnum_renumb (gtopology *top);

  //void ndofn_list (gtopology *top);
  //void dof_list (gtopology *top);
  //void bdof_numbers (gtopology *top);
  //void bdof_list (gtopology *top);
  
  //void create_ltg (gtopology *gt,FILE *out);
  //void ltg1array (long n,long *aux);
  //void readltg1 (long n,XFILE *in);

  ///  number of subdomains
  ///  it is determined in the constructor
  long ns;
  ///  total number of nodes in the whole problem
  ///  this number cannot be obtained in the case of mesh description = bound_nodes
  ///  this number is determined in the function compute_multiplicity
  long tnnp;
  ///  total number of interface (boundary) nodes
  ///  this number is determined in the function compute_multiplicity
  long tnbn;
  ///  number of nodes in the class gtopology
  ///  this number is not equal to tnnp
  ///  the number is determined in the function coupled_dofs
  long nn;
  ///  total number of all elements in the problem
  long tnep;

  ///  mesh description
  ///  it is determined in the constructor
  meshdescription md;
  
  ///  array containing number of nodes on subdomains
  ///  it contains ns components
  ///  nnsd[i]=j - j nodes are defined on the i-th subdomain
  ///  array is read in the function read_nnsd (XFILE *in)
  long *nnsd;
  
  ///  array containing list of numbers of boundary/interface nodes on subdomains
  ///  it contains ns components
  ///  nbnd[i]=j - the i-th domain contains j boundary/interface nodes
  ///  array is assembled in the function compute_multiplicity (FILE *out)
  long *nbnd;

  ///  array containing list of numbers of internal nodes on subdomains
  ///  it contains ns components
  ///  nind[i]=j - the i-th domain contains j internal nodes
  ///  array is assembled in the function compute_multiplicity (FILE *out) 
  long *nind;

  ///  array containing number of subdomains which each node belongs to
  ///  it contains tnnp components,
  ///  amultip[i]=j - the i-th node belongs to j subdomains
  ///  array is assembled in the function assemble_multip (FILE *out)
  long *amultip;

  ///  array containing number of subdomains which each boundary/interface node belongs to
  ///  it contains tnbn components,
  ///  bmultip[i]=j - the i-th interface/boundary node belongs to j subdomains
  ///  array is assembled in the function assemble_multip (FILE *out)
  long *bmultip;

  ///  array containing number of subdomains which share the nodes
  ///  it contains ns, nnsd[i] components
  ///  nodmultip[i][j]=k - the j-th node of the i-th subdomain belongs to k subdomains
    ///  array is assembled in the function update_multip (gtopology *top,FILE *out) or assemble_nodmultip (FILE *out)
  long **nodmultip;
  
  ///  local to global correspondence
  ///  ltg[i][j]=k - the j-th node on the i-th subdomain has global/coarse number k
  ///  if mesh description is bound_nodes:
  ///  k>-1 - coarse number of the node
  ///  k=-1 - the node is internal node
  ///  if mesh description is all_nodes:
  ///  k is the global number of the node
  long **ltg;

  ///  array containing the number of elements on subdomains
  ///  it has ns components
  ///  ned[i]=j - there are j elements on the i-th subdomain
  ///  array is allocated and assembled in function elem_lists
  long *ned;
  
  ///  element-domain correspondence
  ///  eldom[i]=j - the i-th element belongs to the j-th subdomain/aggregate
  long *eldom;
  
  ///  domain-element correspondence
  ///  it has ns rows and ned[i] columns
  ///  domel[i][j]=k - the j-th element on the i-th subdomain has number k
  long **domel;
 
  ///  array containing local numbers of interface/boundary nodes
  ///  it contains ns, nbnd[i] components
  ///  lnbn[i][j]=k - the j-th boundary node on the i-th subdomain has local number k
  ///  array is allocated in the function node_local_numbers
  long **lnbn;

  ///  array containing local numbers of internal nodes
  ///  it contains ns, nind[i] components
  ///  lnin[i][j]=k - the j-th internal node on the i-th subdomain has local number k
  ///  array is allocated in the function node_local_numbers
  long **lnin;

  ///  array containing global glued numbers of boundary/interface nodes
  ///  it contains ns, nbnd[i] components
  ///  ggnbn[i][j]=k - the j-th boundary/interface node on the i-th subdomain has global glued number k
  ///  array is allocated in the function node_global_glued_numbers
  long **ggnbn;

  ///  array containing global glued numbers of internal nodes
  ///  it contains ns, nind[i] components
  ///  ggnin[i][j]=k - the j-th internal node on the i-th subdomain has global glued number k
  ///  array is allocated in the function node_global_glued_numbers
  long **ggnin;

  ///  array containing coarse numbers of interface/boundary nodes
  ///  it contains ns, nbnd[i] components
  ///  icnbnmas[i][j]=k - the j-th boundary node on the i-th subdomain has coarse number k
  ///  array is allocated in the function node_coarse_numbers
  long **icnbnmas;
  
  ///  array containing node multiplicity of boundary/interface nodes
  ///  it contains tnbn components
  ///  icmultip[i]=j - the i-th boundary/interface node belongs to j subdomains
  ///   array is allocated in the function node_coarse_numbers
  long *icmultip;


  ///  local numbers of boundary/interface nodes appropriate to coarse node
  ///  it contains local numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  lnbncn[i][j]=k - the j-th node shared by the i-th coarse node has local number k
  ///  aray is allocated in the function node_coarse_local_map
  long **lnbncn;

  ///  global glued numbers of boundary/interface nodes appropriate to coarse node
  ///  it contains global glued numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  ggnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global glued number k
  ///  aray is allocated in the function node_coarse_global_glued_map
  long **ggnbncn;

  ///  subdomain id of interface/boundary nodes appropriate to coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  sid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the function node_coarse_local_map or node_coarse_global_glued_map
  long **sid;

  
  ///  the number of coupled DOFs
  ///  it is determined in the function coupled_dofs
  long ncdof;
  
  ///  array containing indicators of coupled DOFs
  ///  it has ns, ncdof components
  ///  coupdof[i][j]=k - the j-th coupled DOFs on the i-th subdomain is shared by k nodes
  ///  array is allocated in the function coupled_dofs
  long **coupdof;

  ///  array containing suspicious indicators of coupled DOFs
  ///  it has ncdof components
  ///  coupdofmas[i]=0 - the i-th coupled DOF is not a boundary/interface DOF
  ///  coupdofmas[i]=1 - the i-th coupled DOF is a boundary/interface DOF
  ///  array is allocated in the function coupled_dofs
  long *coupdofmas;
  

  ///  array containing DOF indicators
  ///  if there are coupled DOFs, it is not enough to deal with nodes
  ///  DOFs have to be split to internal and boundary/interface
  ///  it has nn, ndofn[i] components (nn is the total number of nodes, gt->nn)
  ///  dofind[i][j]=0 - the j-th DOF in the i-th node is internal
  ///  dofind[i][j]=1 - the j-th DOF in the i-th node is boundary/interface
  ///  array is allocated in the function update_multiplicity
  long **dofind;
  
  
};

#endif
