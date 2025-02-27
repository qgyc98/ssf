#ifndef FIXNODESEL_H
#define FIXNODESEL_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include "gtopology.h"
#include "pgalias.h"
#include "galias.h"
#include "partop.h"

/**
   This class provides searchnig and condensing of fixing nodes for DP-FETI method

   JB, 2007-2010
 */

class fixnodesel
{
 public:
  // Constructor of class fixnodesel
  // It is called form psolver class
  fixnodesel (int np,int mr,long nd,meshdescription meshd,long* domproces,char *nameproc,int nl,long mes);
  // Desctructor of class fixnodesl
  // It is called form psolver class
  ~fixnodesel();
  
  // It is called form psolver class
  bool check_ltg(long *ltg,long nnd,FILE *out);
  // Function initilalize basic arrays and variables for fixnodesel class
  // It is called form partop class
  void initiate(gtopology *topol,partop *part,FILE *out);
  
  double compute_length_of_vector(long a,long b);
  
  double compute_angle_of_vector_a_and_b(long a,long b,long c);
  
  double compute_area_of_triangle(long a,long b,long c);
  
  void check_triangle();  
  
  // assemble a fictitious subdomain around the real subdomain
  void set_fictitious_subdomain();
  // compute size of fictitious subdomain
  void compute_size_of_fictitious_subdomain();

  long compute_number_of_combination(long n, long k); 
  
  void compute_statistics_of_multiplicity();
  
  void give_whole_dim();
  
  void create_link_dom_nod();
  
  void assemble_ltg_with_fixings(long *ltg);
  
  void create_subdom_graph_3d ();
  
  void create_master_graph();
  
  void select_fixing_nodes_on_master(); 
  
  void fixing_detection(long *ltg);

  
  
  void create_subdom_graph_2d();
  
  void select_minimal_number_2d();
  
  void check_minimal_number_2d();  
  
  void select_fixing_nodes_geom_2d(long auxnv);

  void select_minimal_number_3d();  
  
  void check_minimal_number_3d();  
  
  void select_fixing_nodes_geom_3d(long auxnv); 
 
  void reduce_master_graph_curve();
  
  void reduce_master_graph_surface();
  void reduce_master_graph_surface_2();
  

  //  void compute_centre_of_gravity_surfaces_on_subdomain();
  
  /// void send_centre_of_gravity_on_master();
  
  void set_curves();
  
  void set_curves_3D( );
  
  void print_info_curves();
  void print_info_surfaces();
  
  void add_fixings();
  
  void add_rand();
  
  void add_each_nth_member();
  
  void add_nth_member();
  
  void add_all_members();
  
  void add_centroid();
  
  void add_intpoint();
  
  void add_n_part_curve();
  
  void add_n_rand_nodes_on_surf();
  
  void add_additional_fixing_nodes();
  
  void add_user_pos_def();
  
  void add_all_surfnod();

  void add_centroid_surface();

  void select_centers_of_surfaces();
   void select_centers_of_surfaces_2();

  void mark_surfnodes();

  void mark_ring_nodes();
  
  void add_n_th_mark();

  void select_optimum_2d();

  void add_rings();
  
  void add_max_ring();
  
  void add_boundary_of_surface();
  void add_int_points_surface();
  
  
  void add_max_triangle_on_surface();

  void order_selected_ring_nodes();
  // *************************************************
  //  VARIABLES AND ARRAYS DEFINED ON ALL PROCESSORS
  // *************************************************

  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  /// name of processor
  char procName[10000];
  //char procName[MPI_MAX_PROCESSOR_NAME];
  /// leght of name of processor
  int nameLength;
  ///  number of domain attached to the processor
  long ndom;
  ///  number of nodes on subdomain
  long nn;
  ///  number of elements on subdomain
  long ne;
  
  // pointer to general topology class
  gtopology *top;
  
  ///  topology for parallel computation
  partop *ptop;
 
  ///  number of boundary nodes on subdomain
  long nbn;

  ///  maximum number of boundary nodes on one subdomain
  long maxnbn;

  // number of internal nodes
  long  nin;
  
  ///  mesh description
  meshdescription md;
  /// spatial dimension of mesh 
  long dim;

  // domain - processor corespondence
  long *domproc;

  ///  array containing number of subdomains which share the nodes
  ///  it contains nn components
  ///  nodmultip[i]=j - the i-th node belong to j subdomains
  long *nodmultip;
  
  
  ///  array containing numbers of boundary nodes
  ///  it contains nbn components
  ///  lnbn[i]=j - the i-th boundary node has local number j
  long *lnbn;

  ///  array containing numbers of internal nodes
  ///  it contains nin components
  ///  lnin[i]=j - the i-th internal node has local number j
  long *lnin;
  
  
  ///  array containing list of global numbers of boundary nodes
  ///  lgnbn[i]=j - the i-th boundary node has global number k
  long *lgnbn;

  /// message printing
  /// mespr = 1 yes
  /// mespr = 0 no
  /// It is obtained from @class psolver
  long mespr;
     


  /// condensation of fixing node in FETI-DP method
  long condfixing;
  //// method of condensation
  long methodcondcor;
  /// type automatic and user defined condensation of fixing node in FETI-DP method
  long typecondcur;
  
  long typecondsurf;
    
  long nmembercur;
  
  long nmembersurf;

  
  long nring;

  long *ring;
  
  //
  long *nadjacboundnod;
  // 
  long **adjacboundnod;
  //
  long *boundnodes;
  
  long nuserdefnod;
  long *userdefnod;
  
  long corGeomSearch;
  
  long minSelStrat;

  
  long **glinknod;
  
  long **glinkdom;
  
  long *loclinknod;
  
  long **loclinkdom;
  
  // array contaiting node identification on subdomain
  // it contains nn components
  // nodeidentif[i] = j for the i-th node on subdomain 
  // j = 1 - internal node
  // j = 2 - edge node
  // j = 3 - vertex node
  long *nodeidentif;
  
  
  
  // array with minimal and maximal coordinates of fictitious subdomain
  // fictdom[0] = k - min_x coordinates
  // fictdom[1] = k - min_y coordinates
  // fictdom[2] = k - min_z coordinates
  // fictdom[3] = k - max_x coordinates
  // fictdom[4] = k - max_y coordinates
  // fictdom[5] = k - max_z coordinates
  double *fictdom;
  // 
  double sizefictdom;
  // size of side of fictitious subdomain
  double lx,ly,lz;
  
  // ratio
  // lx/ly
  // lx/lz
  // ly/lz;
  double *ratio;
  
  long *curnadjac;
  long **curadjac;
  long *curidentif;
  long ncurnod;
  
  // ********************************************************
  //  VARIABLES AND ARRAYS DEFINED ONLY ON MASTER PROCESSOR
  // ********************************************************
  
  ///  total number of boundary nodes
  ///  variable is defined in function find_boundary_nodes
  long tnbn;
  
  
  ///  array containing number of subdomains which each boundary node belongs to
  ///  it contains tnnp components
  ///  multip[i]=j - the i-th node belongs to j subdomains
  ///  array is allocated in function compute_multiplicity
  long *multip;
  
  ///  array containing list of numbers of boundary nodes on subdomains
  ///  it contains nproc components
  ///  nbnd[i]=j - the i-th domain contains j boundary nodes
  ///  array is allocated in function find_boundary_nodes
  long *nbnd;

  
  ///  array containing coarse node numbers of boundary nodes
  ///  it contains nproc rows and nbnd[i] columns
  ///  cnbn[i][j]=k - the j-th boundary node on the i-th subdomain has coarse node number k
  ///  array is allocated in function boundary_nodes_on_master
  long **cnbn;

  
  
  // maximal nodal multiplicity in whole problem
  long maxmultip;
  
  // Master
  long **statdom;
  
  // Master
  //
  // nstatdom[i] = j on i-th subdomain is 
  long *nstatdom;
  
  

  double tol;

  // pro nove vyhledavani
  long *coarsenadjac;

  long **coarseadjac;
  
  long *coarseidentif;
  
  
   
  
  FILE *out;
  
  
  long *surfidentif;
  
  // number of surface nodes
  long nsurfnod;
  
  long *surfnod;
  
  // number of curves
  long ncurves;
  // list of start nodes of curves
  long *start;
  // list of end nodes of curves
  long *end;
  // list of members of curves
  long **members;
  // 
  long *nmembers;
  //
  long *nedges;
  
  // number of subdomain where the curves lies
  long *ndomid;
  // 
  long **domid;

  // corespondence between coarse numbers of members of curves
  long ***memcor;

  //reduce flag
  long reduce_flag;
  
  long *curnod;
  
  long nsurf;
  
  long **surfmembers;
  
  long *nsurfmembers;
  
  long **surfdom;
  
  long **surfnodmark;

  long *subdomneib;
  
  long nsubdomneib;

  long **subdomsurf;
  
  long *nsurfcentres;
  
  long **surfcenters;

  long *surfnodpoint;


  long *automember;
  double **realcg;
  
};
#endif
