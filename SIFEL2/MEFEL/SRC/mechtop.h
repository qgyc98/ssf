#ifndef MECHTOP_H
#define MECHTOP_H

#include <stdio.h>
#include "galias.h"
#include "alias.h"
#include "iotools.h"
class node;
class element;
class endnodem;
class edgem;
class siftop;
class ipmap;
struct matrix;
struct vector;
struct ivector;
class dofrange;
class gnodvalvm;
class sel;



/**
  Class mechtop:
   
  It is one of the 5 most important classes of the program.
  (probdesc, mechtop, mechmat, mechbclc, mechcrsec)

  Class mechtop contains topology data of problem.
   
   
  Created by JK, TKo
*/
class mechtop
{
 public:
  mechtop (void);
  ~mechtop (void);
  void read (XFILE *in);
  void print (FILE *out);

  void searching_hanging_nodes (void);
  void elemprescdisp (void);
  void elempresctemp (long lcid);
  elemtype give_elem_type (long eid);
  long give_ndofe (long eid);
  long give_ndofn (long nid);
  long give_dof (long nid,long n);
  void give_nodal_coord (long nid, vector &c);
  long give_nne (long eid);
  long give_nip (long eid, long ri, long ci);
  long give_intordsm (long eid, long ri, long ci); 
  long give_tnip (long eid);
  long give_totnip (long eid);

  /// determines and sets the maximum number of stress/strain componentsfrom from all elements
  long set_maxncompstr ();

  /// returns the number of stress/strain components of the given element matrix block
  long give_ncomp (long eid,long bid);

  /// returns the total number of stress/strain components of the given element
  long give_tncomp (long eid);

  /// returns the total number of stress/strain components of the given element type
  long give_tncomp (elemtype te);

  /// returns the number of edges on the given element
  long give_ned (long eid);

  /// returns the number of nodes on element edges for the given element
  long give_nned (long eid);

  /// returns the number of surfaces on the given element
  long give_nsurf (long eid);

  /// returns the number of nodes on element surfaces for the given element
  long give_nnsurf (long eid);

  /// returns the number of approximation functions on the given element
  long give_napfun (long eid);

  /// returns spatial dimension of the given element
  long give_dimension (long eid);

  /// returns the maximum spatial dimension from all elements
  long give_maxdimension ();

  /// returns the number of element matrix blocks of the given element
  long give_nb (long eid);

  /// returns the number of element matrix blocks of the given element type
  long give_nb_te (elemtype te);

  /// returns stress/strain state indicator for the given element matrix block
  strastrestate give_ssst (long eid,long bi);

  /// returns degree of approximation of the given element
  long give_degree (long eid);

  /// returns array of edge node numbers for the given element edge
  void give_edge_nodes (long eid,long edid,long *nodes);

  /// returns array of node numbers for the given element
  void give_elemnodes (long eid,ivector &nodes);

  /// returns array of code numbers for the given element
  void give_code_numbers (long eid,long *cn);

  /// returns array of code numbers for the given node
  void give_node_code_numbers (long nid,long *cn);

  /// returns array of code numbers for the given node and layer
  void give_mult_code_numbers (long nid,long lid,long *cn);

  /// returns vector of x coordinates of nodes on the given element
  void give_node_coord1d (vector &x,long eid);
  
  /// returns vectors of x and y coordinates of nodes on the given element
  void give_node_coord2d (vector &x,vector &y,long eid);
  
  /// returns vectors of x and z coordinates of nodes on the given element
  void give_node_coord2dxz (vector &x,vector &z,long eid);
  
  /// returns vectors of x, y and z coordinates of nodes on the given element
  void give_node_coord3d (vector &x,vector &y,vector &z,long eid);

  /// returns coordinates of integration points on the given element
  void give_ipcoord_elem(long eid, matrix &coord);

  /** the function returns id and nat. coordinates of the closest integration point 
      of the elemnt to the given point*/
  long give_closest_ip_ncoord(long eid, double px, double py, double pz, 
                              ipmap *ipm, double iptol, double &xi, double &eta, double &zeta);

  /// returns length of the given element
  double give_length(long eid);

  /// returns the given element area
  double give_area(long eid);

  /// returns the given element volume
  double give_volume(long eid);

  /// returns volume of whole domain
  double give_domain_vol();

  /// searches for the maximum number of stress/strain components from all nodes and elements
  void give_maxncompstr (long &max_ncompstrn, long &max_ncompstre);

  /// searches for the maximum number of other array components from all nodes and elements
  void give_maxncompo (long &max_nncompo, long &max_encompo);

  /// searches for the maximum number of DOFs from all nodes
  long give_maxndofn(); 

  /// returns the number of the local coordinate system defined at the given nodes
  long locsystems (ivector &nod);

  /// searchs for nodes with the prescribed displacements, sets the reaction computation flag and allocates reaction array at those nodes
  void comreacnod (void);

  /// searchs for elements with the prescribed displacements, sets the reaction computation flag at those elements
  void comreacelem (void);

  /// the function prepares nodes and elements for the reaction computation
  void comreac (void);

  /// computes nodal strains with the help of least square method
  void strain_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
			    double *lhs,long dim,long fi,long ncomp,long lcid);

  /// computes nodal stresses with the help of least square method
  void stress_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
			    double *lhs,long dim,long fi,long ncomp,long lcid);

  /// computes nodal stresses with the help of least square method
  void other_nodal_values (ivector &nod,vector &nx,vector &ny,vector &nz,
 				  double *lhs,long dim,long fi,long ncomp);

  /// returns value of the required scalar quantity in the given node
  void give_nodal_quant(long nid, mechquant mqn, mechquant reftensq, long lcid, double &qv);

  /// returns value of the required %vector quantity in the given node
  void give_nodal_quant(long nid, mechquant mqn, long lcid, sel &selcomp, gnodvalvm &nv, vector &qv);

  /// returns value of the required %tensor quantity in the given node
  void give_nodal_quant(long nid, mechquant mqn, mechquant reftensq, long lcid, matrix &qv);

  /// returns value of the required quantity defined in the given node
  void give_nodal_quant(long nid, mechquant mq, mechquant reftensq, long lcid, quantrep qr,
                        sel &selcomp, gnodvalvm &nv, double &sv, vector &vv, matrix &tv);

  void store_code_num_elem(void);
  void alloc_nodes_arrays (void);
  void alloc_nodal_strains(void);
  void alloc_enodes ();
  void alloc_edges ();
  void alloc_prep (long nn, long ne, elemtype *el_type);
  void alloc_meaning ();
  void alloc_growstr ();
  void define_meaning ();
  
  void gencodnumlagrmult (long &n);
  long compare_edges (long *enod, long *enodc, long nned);

  
  //void nodedisplacements ();
  void give_noddispl_1d (ivector &nodes,vector &u);
  void give_noddispl_2d (ivector &nodes,vector &u,vector &v);
  void give_noddispl_3d (ivector &nodes,vector &u,vector &v,vector &w);
  
  void init_from_siftop (siftop *Top);

  void clean_nodes ();

  /// function sets arrays of rections at nodes to zero
  void null_react ();
  /// function computes norm of nodal reactions
  double compute_react_norm();

  void initial_displ (long lcid);
  void initial_displ_nodes (long lcid);
  void save_nodedispl_bc (long lcid, double time, double prev_time);
  /// function saves initial displacements at nodes
  void save_node_inidispl(long lcid, double time, double prev_time);
  /// function saves initial displacements on elements
  void save_elem_inidispl(long lcid, long *ifn);
  /// function cleans array of stresses, strains and internal variables on new elements
  void clean_ip_new_elem();
  void save_nodval(long lcid);
  void restore_nodval (double *lhs,long n);
  void save_nodforce(double *f);
  void restore_nodforce(double *rhs,long n);
  void rhs_save (double *rhs);
  void rhs_restore (double *rhs,long n);
  long mesh_check(void);

  ///  functions for construction of balls circumscribed to elements
  void circumscribed_balls ();
  
  void inter ();

  
  ///  number of nodes
  long nn;
  ///  number of constrained nodes
  long ncn;
  ///  number of elements
  long ne;
  ///  number of layered elements
    //long nle;
  ///  number of layered nodes
  long nln;

  /// maximum number of dofs at node
  long max_ndofn;

  /// maximum number of stress/strain components on elements used in the first order homogenization 
  long max_ncompstr;

  ///  array of instances of the class node
  node *nodes;
  ///  array of instances of the class element
  element *elements;
  ///  array of instances of the class endnodem
  endnodem *enodes;
  ///  array of instances of the class edgem
  edgem *edges;

  /// number of integration points, i.e. dimension of arrays nadjip and adjip (it is required in the destructor)
  long tnip;
  ///  array of numbers of adjacent integration points
  ///  this array must be here, not in gtopology because it requires information
  ///  about materials which are stored on integration points
  long *nadjip;
  ///  array of adjacent integration points
  long **adjip;
  ///  array of distances of adjacent integration points
  double **dist;
  ///  array of adjacent elements for each element
  //long **adjelem;
  ///  total number of edges on domain
  long tned;
  /// adjacent elements for each edge on given domain
    //long **eadjelem;
  
  /// array of node numbers with prescribed displacements (i.e. supports)
  long *nodebcid;  

  ///  array containing nodal displacements - used in growing problems
  double **nodedispl;

  ///  array containing nodal forces - used in growing problems
  double **nodeforce;

  /// domain volume
  double domvol;
  
  ///  radius of particles in lattice3
  double *node_radius;
};

#endif
