#ifndef TRANSTOP_H
#define TRANSTOP_H

#include <stdio.h>
#include "aliast.h"
#include "nodet.h"
#include "elementt.h"
#include "endnodet.h"
#include "edget.h"
#include "genfile.h"
//#include "densmat.h"

class ipmap;

/**
   class transtop
   
   it is one of the 5 most important classes of the program
   (probdesc, transtop, transmat, transbclc, transcrsec)

   class transtop contains topology data of problem

*/
class transtop
{
 public:
  transtop (void);
  ~transtop (void);
  void read (XFILE *in);
  void print (FILE *out);

  elemtypet give_elem_type (long eid);
  long give_dof (long nid, long n);
  void give_elemnodes (long eid, ivector &nodes);
  void give_code_numbers (long eid, long *cn);
  void give_medium_code_numbers (long eid, long ri, long *cn);
  void give_node_code_numbers (long nid, long *cn);
  void give_node_coord1d (vector &x, long eid);
  void give_node_coord2d (vector &x, vector &y, long eid);
  void give_node_coord3d (vector &x, vector &y, vector &z, long eid);
  long give_ndofe (long eid);
  /// returns dofe array of the given element
  long** give_dofe (long eid);
  /// returns ordering array of the given element
  long** give_ordering (long eid);
  long give_ndofn (long nid);
  long give_ncomp (long eid);
  long give_ned (long eid);
  long give_nned (long eid);
  long give_nsurf (long eid);
  long give_nnsurf (long eid);
  long give_nne (long eid);
  long give_nip(long eid, long ri, long ci);
  long give_intordkm(long eid, long ri, long ci);
  long give_intordcm(long eid, long ri, long ci);
  long give_tnip (long eid);
  long give_degree (long eid);
  long give_dimension (long eid);
  void give_end_nodes (long eid, long enid, ivector &nodes);
  void give_edge_nodes (long eid, long edid, ivector &nodes);
  void give_surface_nodes (long eid, long surfid, ivector &nodes);
  void give_nbobjects (long eid, long &bco, long &ncompbco);
  void give_bonodes (long eid, long boid, ivector &nod);

  void alloc_prep(long nn, long ne);
  //double give_integral (long eid,vector &nodval);

  void alloc_nodes (void);
  void alloc_enodes ();
  void alloc_edges ();
  void enodes_init ();
  void edge_init ();
  void edge_init_edval ();
  void compute_jumps (double *rhs);
  void compute_resistance_factor (double *rhs);

  void initial_nodval ();

  void lhs_save (double *lhs, double *lhsi, double *tdlhs);
  void lhs_restore (double *lhs, double *lhsi, double *tdlhs);
  long mesh_check(void);

  void view_factors (FILE *out);
  void mod_view_factors (FILE *out);
  void edge_temperature ();
  void heat_fluxes (double *rhs, FILE *out);
  void t4t4 (long edserid,double *y);
  long give_closest_ip_ncoord(long eid, double px, double py, double pz, 
                              ipmap *ipm, double iptol, double &xi, double &eta, double &zeta);
  void sort_elements();

  /// function assembles %vector of master node weights for hanging nodes
  void assemble_master_node_weights_vec(void);
  /// function searches for hanging nodes and assembles transformation matrices on elements
  void searching_hanging_nodes(void);

  //void node_materials ();
  //void jump_initiation (long tnbn,seqselnodes *selnodfeti);
  ///  jnmap[i][0]=j - the i-th interface/boundary node has local number j on the subdomain with lower id
  ///  jnmap[i][1]=j - the i-th interface/boundary node has local number j on the subdomain with higher id
    //long **jnmap;
  long tnbn;

  ///  number of nodes
  long nn;
  ///  number of constrained nodes
  long ncn;
  ///  number of elements
  long ne;

  ///master node number
  long mnn;

  ///  array of instances of the class nodet
  nodet *nodes;
  ///  array of instances of the class elementt
  elementt *elements;
  ///  array of instances of the class endnodet
  endnodet *enodes;
  ///  array of instances of the class edget
  edget *edges;
  
  ///  number of elements influenced by discontinuity
  long nedis;
  ///  list of numbers of nodes with discontinuities on elements
  ///  number of components is nedis
  ///  lnnd[i]=j - the i-th element influenced by discontinuity contains j nodes with discontinuity
  long *lnnd;
  ///  list of node numbers with discontinuities on elements
  ///  lnd[i][j]=k - the j-th node with discontinuity on the i-th element has number k
  long **lnd;
  
  ///  number of cycles on nodes
  ///  it is intended for freezing cycles
  long *ncycl;
  long *aux;
  
  ///  number of node pairs with jumps
  //long nnj;
  
  //  array of previous values at nodes on interface with jumps
  //  prevval[i][j]=k - the j-th component at the i-th node had value k
  //double **prevval;
    
  ///  the array is allocated in the function view_factors 
  double ***view_factor;
  ///  matrices generated from view factors
  densemat *mvf;
  densemat *mvfjk;

  // array of renumbered elements according DOFs
  long *renel;
};

#endif
