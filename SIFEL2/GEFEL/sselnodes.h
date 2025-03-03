#ifndef SSELNODES_H
#define SSELNODES_H

/**
   class sequential selected nodes
   
   the class serves for problems with subdomains or aggregates
   it manipulates with nodes and unknowns (DOFs) on particular
   subdomains/aggregates and deals with relationships among them
   
   JK, 28.8.2007
*/
class sselnodes
{
 public:
  
  sselnodes (long nd,long k,long *j);
  sselnodes (long nd,long *ii,long **jj);
  ~sselnodes ();
  
  void assemble_list_unknowns (gtopology *gt);
  
  
  ///  number of nodes on subdomain
  long nn;
  
  ///  number of selected nodes
  long nsn;

  ///  maximum number of DOFs on subdomain/aggregate
  long maxndof;


  ///  list of selected nodes - local numbers
  ///  lsnl[i]=j - the i-th selected node has local number j
  ///              j-th node on subdomain is selected as the i-th
  ///  lsn contains nsn components, where nsn is the number of selected nodes on subdomain
  long *lsnl;

  ///  list of selected nodes - global numbers
  ///  lsng[i]=j - the i-th selected node has global/coarse number j
  ///  lsng contains nsn components, where nsn is the number of selected nodes on subdomain
  long *lsng;
  

  ///  nsndom[i]=j - j nodes are selected on the i-th subdomain
  long *nsndom;
  
  ///  group node numbers (see partop.h)
  ///  gnn[i][j]=k - the j-th selected node on the i-th subdomain has group number k
  long **gnn;

  ///  array of numbers of DOFs on subdomains 
  ///  it contains prescribed values before code number generation
  ///  after code number generation, only unknown (unconstrained DOFs) are taken into account
  ///  array is rewritten in the function schur_ordering
  ///  ndofdom[i]=j - the i-th subdomain contains j DOFs
  long *ndofdom;

  ///  code numbers
  ///  cndom[i][j]=k - the j-th DOF on the i-th subdomain/aggregate has group code number k
  long **cndom;


};

#endif
