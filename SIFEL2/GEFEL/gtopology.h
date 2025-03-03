#ifndef GTOPOLOGY_H
#define GTOPOLOGY_H

#include <stdio.h>
#include "iotools.h"
#include "gnode.h"
#include "gelement.h"
#include "lgnode.h"
#include "gphase.h"
#include "floatsub.h"
#include "itablefunct.h"
#include "galias.h"
#include "gedge.h"
#include "endnode.h"
#include "seqtop.h"
#include "ordering.h"
#include "vector.h"
#include "matrix.h"
#include "dofrange.h"


class gmatrix;
class boundbox;

/**
   class gtopology
   
   this class contains general topological data of the problem
   
   JK
*/
class gtopology
{
 public:
  gtopology ();
  ~gtopology ();
  
  void initiate (long **icn,long *nc,long m);
  void alloc_nodes (long m);
  void alloc_elements (long m);
  void alloc_lnodes (long m);
  void alloc_phases (long m);//5.4.2005 TKr
  void alloc_growstr (); /// allocates tgf array at some nodes

  //void alloc_lelements (long m);
  //long give_nl (long eid);
  long give_ndofn (long nid);
  long give_nmult (long lnid);
  long give_ndofe (long eid);
  long give_gndofe (long eid);
  long give_ndofedge (long edid);
  long give_nnedge (long edid);
  long give_nne (long eid);
  long give_original_nne (long eid);
  long give_nmne (long eid);
  long give_cne (long eid);
  long give_dof (long nid,long n);
  void save_dof (long nid,long n,long num);
  long give_ldof (long nid,long lid,long n);
  void save_ldof (long nid,long lid,long n,long num);
  void give_original_nodes (long eid,ivector &nod);
  void give_nodes (long eid,ivector &nod);
  void give_master_nodes (long eid,ivector &nod);
  void give_node_coord(long nid, vector &coord);
  void give_node_coord1d (vector &x,long eid);
  void give_node_coord2d (vector &x,vector &y,long eid);
  void give_node_coord2dxz (vector &x,vector &z,long eid);
  void give_node_coord3d (vector &x,vector &y,vector &z,long eid);
  void give_code_numbers (long eid,long *cn);
  void give_cn (long eid,long *cn);
  void give_gcode_numbers (long eid,long *cn);
  void give_node_code_numbers (long nid,long *cn);
  void give_gnode_code_numbers (long nid, ivector &cn);
  void give_mult_code_numbers (long nid,long lid,long *cn);
  void give_edge_code_numbers (long edid,long fln,long *ncn1,long *ncn2,long *mcn);
  void give_endnode_code_numbers (long enid,long *ncn1,long *ncn2,long *mcn);

  gelemtype give_elem_type (long eid);
  gtypel give_siftopelem_type (long eid);
  void give_end_nodes (long eid,long *nodes);
  void give_edge_nodes (long eid,long edid,long  *nodes);
  void give_edge_loc_nodes (long eid,long edid,ivector &edgenod);
  void give_surf_nodes (long eid,long surfid,long *nodes);

  long give_nen (long eid);
  long give_ned (long eid);
  long give_nned (long eid);
  long give_nsurf (long eid);
  long give_nnsurf (long eid);
  long give_degree (long eid);
  long give_whole_dim(long eid);
  long give_maxdimension ();
  void backup_cn ();
  //long codenum_generation (slesolv *ssle,FILE *out);
  long codenum_generation (FILE *out);
  void codenum_print(FILE *out);
  void codenum_elem_print(FILE *out);
  void print_dof_diff(FILE *out);
  long gencodnum ();
  long schur_ordering ();
  long saddlepoint_ordering (long ns,long *nnsd,long **ltg);
  
  void edge_localization (gmatrix *gm);
  
  
  void node_mark (long n);

  void comptop (gtopology *gt1,gtopology *gt2);
  void cndistr (gtopology *gt1,gtopology *gt2);
  void adjacelem (FILE *out);
  
  ///  function assembles list of adjacent nodes to every node
  void adjacnodes (FILE *out);
  void adjacnodes_edge (FILE *out);
  
  /// function returns id of node which is nearest to the point given by coordinates
  long give_nearest_node(double x, double y, double z);
    
  void backup_codnum (void);
  void restore_codnum (void);
  void write_nodes(FILE *out);
  void unodelnode ();
  void initiate_elemmult ();
  void comp_domain_sizes ();
  void give_domain_sizes (double *sizes);
  /// compute square of maximum distance of given point from nodes on element
  double max_sqrdist_nod_pt(long eid, vector &pt);
  

  void read_gf (XFILE *in);
  void print_gf (FILE *out);

  void auxinf_init ();
  void search_newdofs (double time);

  //void restore_leso ();

  void print_time_functions (FILE *kon,double time);

  //void gener_phases (long ncp,long *nan,long *nae);
  
  //  3.11.2006
  void update_elem (double time);
  void update_nodes ();
  /// function updates state indicator of the elements from the previous time step
  void update_auxinf();
  /// function removes nodes according to list of removed nodes sn
  void remove_nodes(long *sn);
  /// function searches for elements with changed status
  long search_changed_elem(double time, long &nae);
  /// function searches for nodes between new and old part of topology
  long search_iface_nodes(long *ifn);
  /// function searches for DOFs with changed status
  long search_changed_dofs (double time,double prev_time);
  /// function switches on new elments only with respect to actual time step
  void switch_new_elem();
  /// function switches on removed elments only with respect to actual time step
  void switch_removed_elem();
  /// function sets status indicator of nodes according to status indicator of adjacent elements
  void update_active_dofs (double time);
  /// function sets DOF numbers to zero on selected interface nodes given by ifn list
  void clear_intf_dofs(long *ifn);

  long search_newelem (double time,double prev_time);
  long search_newdofs (double time,double prev_time);
  void update_dofs (double time);
  //void update_mesh (double time,double prev_time,long &ncn,long &nce);
  
  void lneso_init ();  


  
  ///  saddle point problems
  ///  e.g. hemivariational inequalities or discontinuities in transport processes
  void end_nodes (FILE *out);
  void endnodes_auxprint (FILE *out);
  void alloc_endnode_cn ();
  void endnodes_localization (gmatrix *gm);
  void enodes_allelem (FILE *out);
  
  void edges (FILE *out);
  void edges_mat (long matrad,FILE*);
  void edge_dirvect ();
  void edge_normvect ();
  void edgenode_sorting ();
  void prev_next_edges ();
  void edge_series (FILE *out);
  void edge_elem (FILE *out);
  void edge_allelem (FILE *out);
  void normvectorient ();
  //void cnmodif (FILE *out);
  void edges_auxprint (FILE *out);
  
  void alloc_edge_cn ();

  void auto_subdom_detection (FILE *out);
  long codenum_multip ();
  void mult_localization (gmatrix *gm);

  void domdecomp ();
  void create_ltg (FILE *out);
  void read_seq_top (XFILE *in);
  
  ///  renumbering
  void cuthill_mckee_renumb (FILE *out);

  void writePriorityQueue(int *fronta);
  void sloan_renumb (FILE *out);
  void shell_sort( int *array, int arrayLength);
  void shell_sort_x(ivector &x_ord);
  void shell_sort_y(ivector &y_ord);
  void shell_sort_z(ivector &z_ord);
  void lastLevelOfRootedStructure(int start, int *velikost, int *maxSirka, int*hloubka, int *vzdalenost, int **lastLevel);
  int findMinimumDegree();
  void removeRepetitiousDegrees(int **pole, int *velikost);
  void generateInitialList(int start, int **seznam, int * velikostSeznamu, int *hloubkaSeznamu, int *maxSirka);
  void findPseudoPheripheral(int *startUzel, int *cilUzel);

  
  void searching_hanging_nodes (void);
  void dof_transfer_hangnod (void);
  void approx_weights (gtypel et, long nid, vector &w);
  void hang_nodes_check (void);

  /// addtional code numbers are added on elements due to stress based homogenization
  void code_num_mod (ivector &cnadd, long &ndof);
  
  //  tools for determination of hanging nodes
  void construct_circumscribed_balls ();
  void element_center (long eid);
  void radius_circumscribed_ball (long eid);
  void allocate_arrays_circumscribed_ball ();
  double give_characteristic_size (long eid);

  void bar_linhex_intersection (long eid,long &numinter,vector &xb,vector &yb,vector &zb,vector &x,vector &y,vector &z,double zero,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);
  void bar_linhex_intersec_jacobian (double alpha,double beta,vector &v2,vector &v3,vector &v4,vector &v5,matrix &jac);
  void bar_linhex_intersec_function (double alpha,double beta,double gamma,vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,vector &f);
  long newton_intersection (double zero,double &alpha,double &beta,vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,matrix &jac);
  void bar_quadrangle_surface_intersection (long &numinter,long sid,ivector &nod,vector &xb,vector &yb,vector &zb,vector &xh,vector &yh,vector &zh,
                                            vector &v1,vector &v2,vector &v3,vector &v4,vector &v5,double zero,
                                            vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);
  void bar_edge_intersection (long &numinter,long edgeid,ivector &nod,vector &xb,vector &yb,vector &zb,vector &xh,vector &yh,vector &zh,double zero,
                              vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);

  void surface_linhex_natcoordinates (long sid,double alpha,double beta,double zero,ivector &nod,vector &xh,vector &yh,vector &zh,
                                      long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);
  void edge_linhex_natcoordinates (long edgeid,double alpha,double zero,ivector &nod,vector &xh,vector &yh,vector &zh,
                                   long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);
  void end_points_edge_linhex_natcoordinates (long edgeid,double alpha,ivector &nod,vector &xh,vector &yh,vector &zh,
                                              long &numinter,vector &xi,vector &yi,vector &zi,vector &xxi,vector &yyi,vector &zzi,imatrix &masnod);
  
  void initiate_omp_elem_order();
  void initiate_elem_dof_ranges();
  void give_bounding_box(boundbox &bb);
  
  enum Status{
    ACTIVE = 1, 
    INACTIVE = 2, 
    PREACTIVE= 3,
    POSTACTIVE = 4
  };
  
  
  gnode *gnodes;
  gelement *gelements;
  lgnode *lgnodes;
  gphase *gphases; // candidate for removal
  
  
  ///  the number of nodes
  long nn;
  ///  the number of elements
  long ne;
  ///  the number of layered nodes
  long nln;
  ///  the number of unknowns defined on the mesh
  long ndof;
  ///  the number of internal DOFs
  long nidof;
  ///  the number of boundary/interface DOFs
  long nbdof;
  ///  the number of unknowns in saddle point problem
  long nsad;
  ///  the number of nodes which are switched on
  ///  it is used in problems with changing number of nodes and elements
  long nnso;
  ///  the number of elements which are switched on
  ///  it is used in problems with changing number of nodes and elements
  long neso;
  ///  the number of subdomains/aggregates
  long ns;
  ///  the number of phases
  long nph;
  ///  presence of hanging nodes
  ///  hangnod=0 - there are no hanging nodes
  ///  hangnod=1 - there are hanging nodes
  long hangnod;
  
  
  ///  arrays of coordinates of element centers
  double *xc,*yc,*zc;
  ///  array of radii of circumscribed balls
  double *circumrad;
    
    
  /**
     indicator of functions for DOFs control
     dofcontr=0 - classical approach
     dofcontr=1 - DOFs are controlled by function
  */
  long dofcontr;

  /**
     state of code numbers = state of equation numbers
     this variable is used especially in problems with changing number of nodes, elements and DOFs
     DOFs are switched on and off during computation
     
     cnstate=0 - code numbers are not generated, they are not available
     cnstate=1 - code numbers are generated, they are available
  */
  long cnstate;
  
  /**
     array contains reordered node numbers
     ordering[i]=j - the i-th reordered node has number j before reordering
  */
  long *ordering;
  ///  type of node renumbering
  noderenumb nodren;
  
  ///  the number of functions for DOFs control
  long ngf;
  ///  functions for DOFs control
  gfunct *gf;
  
  /**
     type of generator of code numbers
     cngen = 1 - the classical code number generation
     cngen = 2 - the Schur system of code number generation
                 internal DOFs are ordered first, interface/boundary DOFs last
     cngen = 3 - ???
  */
  long cngen;

  ///  sizes of solved domain
  double *domsizes;




  // ***********************
  //  ELEMENTS - NODES MAP
  // ***********************

  ///  array of numbers of adjacent elements to nodes
  ///  nadjelnod[i]=j - the i-th node is shared by j elements
  ///  array is allocated in the function adjacnodes or adjacelem
  long *nadjelnod;

  ///  array of adjacent elements to nodes
  ///  adjelnod[i][j]=k - the j-th element adjacent to the i-th node has number k
  ///  array is allocated in the function adjacnodes or adjacelem
  long **adjelnod;
  
  // **************************
  //  ELEMENTS - ELEMENTS MAP
  // **************************
  
  ///  array of numbers of adjacent elements to elements
  ///  nadjelel[i]=j - there are j adjacent elements to the i-th element
  ///  array has to be sorted, there are multiple contributions
  ///  array is allocated in the function adjacelem
  long *nadjelel;

  ///  array of adjacent elements to elements
  ///  adjelel[i][j]=k - the j-th element adjacent to the i-th element has number k
  ///  array has to be sorted, there are multiple contributions
  ///  array is allocated in the function adjacelem
  long **adjelel;

  // ***********************
  //  NODES - NODES MAP
  // ***********************

  ///  numbers of adjacent nodes to nodes
  ///  nadjnodnod[i]=j - there are j adjacent nodes to the i-th node
  ///  array is allocated in the function adjacnodes
  long *nadjnodnod;
  ///  array of adjacent nodes to nodes
  ///  adjnodnod[i][j]=k - the j-th node adjacent to the i-th node has number k
  ///  array is allocated in the function adjacnodes
  long **adjnodnod;
  
  
  long *nadjacnodesedge;
  long **adjacnodesedge;
  
  ///  array with code numbers backup
  long **bckcn;
  ///  usual node - layered node correspondation
  ///  important for layered problems
  long *unln;
  ///  usual node - number of layer correspondation
  ///  important for layered problems
  long *unnl;
  ///  correspondence of code numbers between two generalized topologies
  long *cngtopcorr;
  
  ///  object of floating subdomains
  floatsub flsub;
  
  ///  list of nodes which are switched on
  ///  it is motivated by problem with changing number of nodes and elements
  ///  lnso[i]=0 - the i-th node is switched off
  ///  lnso[i]=1 - the i-th node is switched on
  ///
  ///  all components of array lnso are equal to 1 for all types of problems
  ///  except problems with changing number of nodes and elements
  long *lnso;
  
  ///  list of elements which are switched on
  ///  it is motivated by problem with changing number of nodes and elements
  ///  leso[i]=0 - the i-th element is switched off
  ///  leso[i]=1 - the i-th element is switched on
  ///
  ///  all components of array leso are equal to 1 for all types of problems
  ///  except problems with changing number of nodes and elements
  long *leso;
  ///  backup array for leso components
    //long *buleso;
    
  

  ///  TOOLS FOR PROBLEMS WITH HEMIVARIATIONAL INEQUALITIES
  ///  TOOLS FOR PROBLEMS WITH MATERIAL INTERFACES (e.g. RADIATION)
    
  ///  type of edges
  ///  edtype = jumps - edges are used for problems with jumps (discontinuities on material interfaces
  ///  edtype = mater - edges are used for problems with material interfaces and continuous variables (e.g. radiation)
  ///  the actual value is assigned in the function %radiation_init () called in %trfelinit.cpp in TRFEL
  edgetype edtype;
  ///  the number of general edges
  ///  it is determined in the function edges
  long nged;
  ///  the number of series
  ///  it is determined in the function edge_series
  long nser;

  ///  series of edges
  ///  edgeser[i]=j - the i-th edge is in the j-th series
  ///  array is allocated in the function edge_series
  long *edgeser;
  
  ///  the number of edges in serie
  ///  nedser[i]=j - the i-th series contains j edges
  ///  array is allocated in the function edge_series
  long *nedser;
  
  ///  list of edge numbers in series
  ///  edgelist[i][j]=k - the j-th edge in the i-th serie has number k
  ///  array is allocated in the function edge_series
  long **edgelist;
  
  ///  array of objects of the class gedge
  ///  array is alloctaed in the function edges
  gedge *gedges;

  ///  the number of end nodes
  long nen;

  ///  end nodes, see files endnode.h and endnode.cpp
  ///  array is allocated in the function end_nodes
  endnode *endnodes;
  
  ///  sequential topology, see files seqtop.h and seqtop.cpp
  ///  array of objects is allocated in the function domdecomp
  seqtop *stop;
  
  ///  indicator of sequential topology reading
  ///  rst=0 - sequential topology is not read
  ///  rst=1 - sequential topology is read
  ///  rst=2 - sequential topology is read in the form of Boolean matrices
  long rst;
  
  // Data for support of multithread processing with the help of OpenMP

  /// array of DOF ranges on particular elements (used in multithread environment OpenMP)
  dofrange *eldofr;

  /// array of element order for passing element loops in multithread environment OpenMP
  long *ompelord;
};

#endif
