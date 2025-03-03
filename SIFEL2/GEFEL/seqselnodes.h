#ifndef SEQSELNODES_H
#define SEQSELNODES_H

#include <stdio.h>
#include "galias.h"
#include "gtopology.h"

/**
   class deals with selected nodes on subdomains
   the class is sequential version of the class selnodes
   nodes could be selected for Schur complement method,
   FETI method, etc.
   
   JK, 14.9.2007
*/
class seqselnodes
{
 public:
  seqselnodes(long nd,long *nnsd,long **jj,
	      long **nodmultip,long **ggnbn,long **icnbnmas,
	      long itnbn,long *iicmultip,long **ilnbncn,long **iggnbncn,long **isid,
	      meshdescription d,FILE *out,long mespr);
  ~seqselnodes ();
  
  void node_coarse_numbers (FILE *out);
  void number_all_dofs (gtopology *top,FILE *out);
  void ndofn_on_master (gtopology *top,FILE *out);
  void dof_indicators (gtopology *top,FILE *out);
  
  /*
  void nodes_on_master (long tnnp,FILE *out);
  void node_multiplicity (FILE *out);
  void group_local_nodes (FILE *out);
  void dof_feti (FILE *out);

  void number_contrib (FILE *out);
  void contrib_dofs (gtopology *top,FILE *out);
  */
  //void dof_multiplicity (FILE *out);


  void schur_ordering (long **dofind,FILE *out);
  //void schur_ordering_old (long **dofind,FILE *out);
  //void schur_ordering_old_old (FILE *out);
  void prepare_schur (gtopology *top,FILE *out);

  
  
  ///  function defines the type of implementation of the FETI method
  void define_feti_implementation (fetiimplem fi);
  ///  variable fetiimpl is determined
    
  ///  function determines the number of contributions to the coarse problem from particular subdomains
  void number_contrib_dofs (gtopology *top,FILE *out);
  ///  function assembles array snndofmas
    
  void dof_feti (gtopology *top,FILE *out);

  void contrib_dofs_ln (gtopology *top,FILE *out);
  void contrib_dofs_cn (FILE *out);

  long prepare_feti (fetiimplem fi,gtopology *top,FILE *out);
    
    
    
    
    
    
    
    
  ///  type of mesh description
  ///  defined in the constructor 
  meshdescription md;
  
  ///  type of FETI implementation
  ///  defined in the function define_feti_implementation
  fetiimplem fetiimpl;

  ///  number of subdomains
  ///  defined in the constructor 
  long ns;
  
  ///  total number of boundary/interface nodes
  ///  defined in the constructor 
  long tnbn;

  ///  total number of selected nodes
  ///  defined in the function node_coarse_numbers
  long tnsn;


  ///  numbers of selected nodes on subdomains
  ///  it contains ns components
  ///  nsnmas[i]=j - there are j selected nodes on the i-th subdomain
  ///  array is assembled in the constructor
  long *nsnmas;

  ///  local numbers of selected nodes
  ///  lnsn contains ns,nsnmas[i] components
  ///  lnsn[i][j]=k - the j-th selected node on the i-th subdomain has local number k
  ///  array is allocated in the constructor
  long **lnsn;

  ///  global glued numbers of selected nodes
  ///  ggnsn contains ns,nsnmas[i] components
  ///  ggnsn[i][j]=k - the j-th selected node on the i-th subdomain has global glued number k
  ///  array is allocated in the constructor
  long **ggnsn;

  ///  coarse numbers of selected nodes
  ///  cnsnmas contains ns, nsnmas[i] components
  ///  cnsnmas[i][j]=k - the j-th selected node on the i-th subdomain has coarse number k
  ///  array is allocated in the constructor
  long **cnsnmas;

  ///  number of multiplicity of all boundary/interface nodes
  ///  it contains tnbn components
  ///  icmultip[i]=j - the i-th boundary/interface node shares j nodes (it belongs to j subdomains)
  ///  array is allocated in the constructor
  long *icmultip;
  
  ///  local numbers of boundary/interface nodes appropriate to all coarse nodes
  ///  it contains local numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  lnbncn[i][j]=k - the j-th node shared by the i-th coarse node has local number k
  ///  aray is allocated in the constructor
  long **lnbncn;

  ///  global glued numbers of boundary/interface nodes appropriate to all coarse nodes
  ///  it contains global glued numbers of boundary/interafce nodes of each coarse node
  ///  it contains tnbn rows and icmultip[i] columns
  ///  ggnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global glued number k
  ///  aray is allocated in the constructor
  long **ggnbncn;

  ///  subdomain id of interface/boundary nodes appropriate to all coarse nodes
  ///  it contains tnbn rows and icmultip[i] columns
  ///  sid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the constructor
  long **sid;



  ///  number of multiplicity of selected boundary/interface nodes
  ///  it contains tnsn components
  ///  snicmultip[i]=j - the i-th selected boundary/interface node shares j nodes (it belongs to j subdomains)
  ///  array is alloctaed in the function node_coarse_numbers
  long *snicmultip;
  
  ///  global glued numbers of selected boundary/interface nodes appropriate to coarse node
  ///  it contains tnsn rows and snicmultip[i] columns
  ///  snggnbncn[i][j]=k - the j-th node shared by the i-th coarse node has global glued number k
  ///  aray is allocated in the function node_coarse_numbers
  long **snggnbncn;

  ///  local numbers of selected boundary/interface nodes appropriate to coarse node
  ///  it contains tnsn rows and snicmultip[i] columns
  ///  snlnbncn[i][j]=k - the j-th node shared by the i-th coarse node has local number k
  ///  aray is allocated in the function coarse_local_nodes
  long **snlnbncn;

  ///  subdomain id of selected interface/boundary nodes of the coarse node
  ///  it contains tnsn rows and snicmultip[i] columns
  ///  snsid[i][j]=k - the j-th node shared by the i-th coarse node belongs to the k-th subdomain
  ///  array is allocated in the function node_coarse_numbers
  long **snsid;

  ///  array of numbers of DOFs on subdomains at selected nodes
  ///  it contains prescribed values before code number generation
  ///  it contains ns components
  ///  after code number generation, only unknown (unconstrained DOFs) are taken into account
  ///  snndofmas[i]=j - selected nodes on the i-th subdomain contain j DOFs
  ///  array is allocated in the function number_all_dofs
  long *snndofmas;

  ///  array of numbers of DOF at selected nodes
  ///  it contains ns, nsnmas[i] components
  ///  snndofnmas[i][j]=k - the selected j-th node on the i-th subdomain contains k DOFs
  ///  array is allocated in the function ndofn_on_master
  long **snndofnmas;

  ///  array of DOFs or indicators at selected nodes
  ///  it contains ns, nsnmas[i], snndofnmas[i][j] components
  ///  sndofmas[i][j][k]=l - the k-th DOF at the j-th selected node on the i-th subdomain has value l
  ///  array is allocated in the function dof_indicators
  long ***sndofmas;
  

  
  ///  array of numbers of DOFs for selected nodes
  ///  it contains tnsn components
  ///  ndofnsn[i]=j - the i-th selected node (in group ordering) has j DOFs
  ///  array is allocated in the function schur_ordering or in the function dof_feti
  ///  in the case of Schur complement method, it contains ndofn for each boundary/interface node
  ///  in the case of FETI, it contains (nmultip-1)*ndofn
  long *ndofnsn;

  ///  code numbers at selected nodes on master
  ///  it contains tnsn, ndofnsn[i] components
  ///  cnm[i][j]=k - the j-th DOF at the i-th selected node has code number / indicator k
  ///  array is allocated in the function schur_ordering
  long **codensn;

  ///  total number of DOFs on selected nodes
  ///  defined in the function schur_ordering
  long tndofsn;

  ///  code numbers / indicators for FETI method
  ///  it contains tnsn, snicmultip[i], snicmultip[i], ndofnsn[i] components
  ///  doffeti[i][j][k][l]=m - the l-th DOF in the couple j-th, k-th node shared by the i-th coarse node has code number / indicator m
  ///  array is allocated in the function dof_feti
  long ****doffeti;
  
  ///  coarse code numbers of subdomains considered as superelements
  ///  it contains ns, snndofmas[i] components
  ///  cndofmas[i][j]=k - the j-th DOF on the i-th subdomain has coarse code number k
  ///  array is allocated in the function schur_ordering or contrib_dofs_cn
  long **cndofmas;

  ///  list of local code numbers which contribute to the coarse problem
  ///  it contains ns, snndofmas[i] components
  ///  lndofmas[i][j]=k - the j-th DOF which contributes to the coarse problem on the i-th subdomain has local code number k
  ///  array is allocated in the function contrib_dofs_ln
  long **lndofmas;





















  // nasleduji neoverene promenne a pole

  ///  node multiplicity
  ///  nodmultip[i]=j - the i-th selected node is shared by j subdomains
  ///  it contains tnsn components
    //long *nodmultip;
  
    // long *dofmultip;
    //long **ldofmultip;
  
    //long tndof;





  ///  code numbers at selected nodes on master
  ///  cnm[i][j]=k - the j-th DOF at the i-th selected node has code number / indicator k
    //long **cnm;

  //  FETI ORDERING

  ///  list of joint nodes to selected nodes assumed as coarse nodes
  ///  ljn[i][j]=k - the j-th node connected to the i-th coarse node has local number k
  ///  ljn contains tnsn rows and nodmultip columns
    //long **ljn;
  
  ///  list of subdomain numbers which contain connected nodes to coarse nodes
  ///  lsn[i][j]=k - the j-th node connected to the i-th coarse node belongs to the k-th subdomain
  ///  lsn contains tnsn rows and nodmultip columns
    //long **lsn;

  ///  code numbers / indicators for FETI method
  ///  doffeti[i][j][k]=l - the k-th DOF on the j-th connected node to the i-th coarse node has code number / indicator l
    //long ***doffeti;
  
  ///  number of contributing nodes in the FETI method
  ///  ncndom[i]=j - the i-th subdomain contributes to the coarse problem by j nodes
    //long *ncndom;

  
  //  ncdofd[ns];
  //edofs[ns][ncdofd];
  //ccn[ns][ncdofd];
};

#endif

