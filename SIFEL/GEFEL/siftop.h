#ifndef SIFTOP_H
#define SIFTOP_H

#include <stdio.h>
#include "galias.h"
#include "iotools.h"
#include "gtopology.h"

enum outmodetype {omt_parameters=0, omt_tangent=1, omt_normal=2, omt_boundary_entities=3, 
                  omt_associated_elements=4, omt_neighbourhood=5, omt_el_iso_type=6,
                  omt_nat_order_nod=7, omt_sec_clasif_nod=8};


/// overloaded postfix operator for enum type gtypel
gtypel operator ++(gtypel &a, int);



struct sedges;
struct sfaces;
struct vector;
struct ivector;



/**
  This structure groups information about node like coordinates and property number

  Created by Tomas Koudelka, 2001
*/
struct snode
{
  double x; ///< x coordinate of node
  double y; ///< y coordinate of node
  double z; ///< z coordinate of node
  long     nprop; ///< number of assigned properties
  entityp *entid; ///< entity types
  long    *prop;  ///< property numbers of given entities

  snode();
  ~snode();
  void alloc(long num_prop);  ///< allocation of node
  void add_prop(entityp ent, long propid); ///< adds one new property
  long searchprop(long aprop, gentity ent, long nid, sedges *edg, sfaces *sfc); ///< searches node for given property
  long searchpropent(entityp e); ///< returns number of properties of the given entity type
  void getcoord(vector &p); ///< returns nodal coordinates
  void copyto(snode &c);
};



/**
  This structure groups information about element like node numbers, property number and type of element

  Created by Tomas Koudelka, 2001
*/
struct selement
{
  static const long max_nned = 4;
  static const long max_nnsurf = 12;
  static const long max_nsurfedg = 4;

  long   *nodes;    ///< pointer to array with element node numbers
  long    nne;      ///< number of nodes on element
  gtypel  type;     ///< geometric type of element
  long    ned;      ///< number of edges
  long   *nned;     ///< number of nodes on particular edges
  long    tnned;    ///< total number of nodes on all edges
  long   *propedg;  ///< property numbers of each edge
  long    nsurf;    ///< number of surfaces
  long   *nnsurf;   ///< number of nodes on particular surfaces
  long   *propsurf; ///< property numbers of surfaces
  long    tnnsurf;  ///< total number of nodes on all surfaces
  long    prop;     ///< property number of element
  long   (*edgenod)[max_nned];   ///< array with indices of nodes on particular edges - edgnod[i][j] gives j-th index of node on i-th edge
  long   (*surfnod)[max_nnsurf]; ///< array with indices of nodes on particular surfaces - surfnod[i][j] gives index of j-th node on i-th surface
  long   (*surfedg)[max_nsurfedg]; ///< array with indices of edges on particular surfaces - surfedg[i][j] gives index of j-th edge on i-th surface

  selement();
  selement(gtypel et, long alloc_bp);
  ~selement();

  /// allocation of memory for nodes
  long alloc(long alloc_bp);

  /// copies the content of the given element to the argument instance
  void copyto(selement &dest) const;

  /// searches for property aprop of the given entity ent in element property arrays
  long searchprop(long aprop, gentity ent, snode *tnodes, long eid, sedges *edg, sfaces *sfc, long *&entid, long &nentid) const;
  /// searches for property aprop of the given entity ent in the element nodes
  long searchnodprop(long aprop, gentity ent, snode *tnodes, sedges *edg, sfaces *sfc, long *propent) const;
  long compare_edg(const long *nod, long n) const;
  long compare_edg(const long *nod, long n, ivector &edgid_lst) const;
  long compare_surf(const long *nod, long n) const;
  long compare_surf(const long *nod, long n, ivector &surfid_lst) const;

  // initialization of the following static arrays is defined at the begining of siftop.cpp
  static long nned_isolin1d[1];
  static long nned_isoquad1d[1];
  static long nned_trlin2d[3];
  static long nned_trquad2d[3];
  static long nned_isolin2d[4];
  static long nned_isoquad2d[4];
  static long nned_isocubic2d[4];
  static long nned_tetlin3d[6];
  static long nned_tetquad3d[6];
  static long nned_pyramlin[8];
  static long nned_pyramquad[8];
  static long nned_wedgelin[9];
  static long nned_wedgequad[9];
  static long nned_isolin3d[12];
  static long nned_isoquad3d[12];

  static long edgenod_isolin1d[1][max_nned];
  static long edgenod_isoquad1d[1][max_nned];
  static long edgenod_trlin2d[3][max_nned];
  static long edgenod_trquad2d[3][max_nned];
  static long edgenod_isolin2d[4][max_nned];
  static long edgenod_isoquad2d[4][max_nned];
  static long edgenod_isocubic2d[4][max_nned];
  static long edgenod_tetlin3d[6][max_nned];
  static long edgenod_tetquad3d[6][max_nned];
  static long edgenod_pyramlin[8][max_nned];
  static long edgenod_pyramquad[8][max_nned];
  static long edgenod_wedgelin[9][max_nned];
  static long edgenod_wedgequad[9][max_nned];
  static long edgenod_isolin3d[12][max_nned];
  static long edgenod_isoquad3d[12][max_nned];

  static long nnsurf_isolin1d[1];
  static long nnsurf_isoquad1d[1];
  static long nnsurf_trlin2d[1];
  static long nnsurf_trquad2d[1];
  static long nnsurf_isolin2d[1];
  static long nnsurf_isoquad2d[1];
  static long nnsurf_isocubic2d[1];
  static long nnsurf_tetlin3d[4];
  static long nnsurf_tetquad3d[4];
  static long nnsurf_pyramlin[5];
  static long nnsurf_pyramquad[5];
  static long nnsurf_wedgelin[5];
  static long nnsurf_wedgequad[5];
  static long nnsurf_isolin3d[6];
  static long nnsurf_isoquad3d[6];

  static long surfnod_isolin1d[1][max_nnsurf];
  static long surfnod_isoquad1d[1][max_nnsurf];
  static long surfnod_trlin2d[1][max_nnsurf];
  static long surfnod_trquad2d[1][max_nnsurf];
  static long surfnod_isolin2d[1][max_nnsurf];
  static long surfnod_isoquad2d[1][max_nnsurf];
  static long surfnod_isocubic2d[1][max_nnsurf];
  static long surfnod_tetlin3d[4][max_nnsurf];
  static long surfnod_tetquad3d[4][max_nnsurf];
  static long surfnod_pyramlin[5][max_nnsurf];
  static long surfnod_pyramquad[5][max_nnsurf];
  static long surfnod_wedgelin[5][max_nnsurf];
  static long surfnod_wedgequad[5][max_nnsurf];
  static long surfnod_isolin3d[6][max_nnsurf];
  static long surfnod_isoquad3d[6][max_nnsurf];

  static long surfedg_isolin1d[1][max_nsurfedg];
  static long surfedg_isoquad1d[1][max_nsurfedg];
  static long surfedg_trlin2d[1][max_nsurfedg];
  static long surfedg_trquad2d[1][max_nsurfedg];
  static long surfedg_isolin2d[1][max_nsurfedg];
  static long surfedg_isoquad2d[1][max_nsurfedg];
  static long surfedg_isocubic2d[1][max_nsurfedg];
  static long surfedg_tetlin3d[4][max_nsurfedg];
  static long surfedg_tetquad3d[4][max_nsurfedg];
  static long surfedg_pyramlin[5][max_nsurfedg];
  static long surfedg_pyramquad[5][max_nsurfedg];
  static long surfedg_wedgelin[5][max_nsurfedg];
  static long surfedg_wedgequad[5][max_nsurfedg];
  static long surfedg_isolin3d[6][max_nsurfedg];
  static long surfedg_isoquad3d[6][max_nsurfedg];
};



/**
  This structure groups information about one edge like node numbers and property number

  Created by Tomas Koudelka, 08.2011
*/
struct sedge
{
  long     n1; ///< the first node number of the edge
  long     n2; ///< the second node number of the edge
  long     prop;  ///< property number of the given edge

  sedge();
  void copyto(sedge &dest); ///< method copies the given instance to the argument instance dest
  long read(XFILE *in);  ///< method reads the edge data from the file
  void print(FILE *out, long *lnn); ///< method prints the edge data to the file
  long cordomtest(long *lnn); ///< method tests whether the edge is contained in the given subdomain
  
};



/**
  This structure handles list of edges

  Created by Tomas Koudelka, 08.2011
*/
struct sedges
{
  long     nedg;   ///< number of edges
  sedge    *edges; ///< array of edges
  long     nprop;  ///< number of different properties for all faces
  long     *proplist; ///< array of particular property numbers (its size is nprop)

  /// arrays of node numbers for particular edge properties (prop_list_nod[i][j] = j-th node of edge with property proplist[i]
  long     **prop_list_nod; 

  /// arrays of number of nodes for particular edge properties (prop_list_nod_num[i] = number of nodes with edge property proplist[i] i.e. length of prop_list_nod[i])
  long     *prop_list_nod_num;

  /// arrays of element numbers for particular edge properties (prop_list_elem[i][j] = j-th element of edge with property proplist[i])
  long     **prop_list_elem;

  /// arrays of number of elements for particular edge properties (prop_list_elem_num[i] = number of elements with edge property proplist[i] i.e. length of prop_list_elem[i])
  long     *prop_list_elem_num;

  /// arrays of element edge id for particular edge properties and elements (prop_list_elem_edgid[i][j][k] = k-th edge id of j-th element with edge property proplist[i])
  long     ***prop_list_elem_edgid;

  /// arrays of number of element edges with particular edge properties for particular elements (prop_list_elem_edgid_num[i][j] = number of edges of j-th element with edge property proplist[i] i.e. length of prop_list_elem_edgid[i][j]))
  long     **prop_list_elem_edgid_num;

  sedges();
  sedges(long ned); ///< initializes to array of edges with size ned
  ~sedges();

  void copyto(sedges &dest); ///< method copies the given instance into the argument instance dest  
  void read(XFILE *in);  ///< method reads the data about edges from the file
  void print(FILE *out, long *lnn); ///< method prints the data about edges to the file
  long search_list_newprop(long prop, long *propid, long n); ///< tests prop whether it is contained in the array propid
  void gen_list_edg(long nn, long ne, selement *elemss, long *nadjelnod, long **adjelnod);  ///< method generates the lists of nodes and elements for particular properties of the edges
  void gen_list_node(long nn);  ///< method generates the lists of nodes for particular properties of the edges
  void gen_list_elem(long ne, selement *elems, long *nadjelnod, long **adjelnod);  ///< method generates the lists elements for particular properties of the edges
};



/**
  This structure groups information about one face like node numbers and property number

  Created by Tomas Koudelka, 08.2011
*/
struct sface
{
  long     nnod; ///< number of nodes that define the given face
  long    *nodes; ///< array of node numbers defining the face
  long     prop;  ///< property number of the given face

  sface();
  sface(long n); ///< initializes to array of nodes with size n
  ~sface();

  void copyto(sface &dest); ///< method copies the given instance into the argument instance dest
  long read(XFILE *in);  ///< method reads the edge data from the file
  void print(FILE *out, long *lnn); ///< method prints the edge data to the file
  long cordomtest(long *lnn); ///< method tests whetehr the edge is contained in the given subdomain
};



/**
  This structure handles list of faces

  Created by Tomas Koudelka, 08.2011
*/
struct sfaces
{
  long     nfcs;   ///< number of faces
  sface    *faces; ///< array of faces
  long     nprop;  ///< number of different properties for all faces
  long     *proplist; ///< array of particular property numbers (its size is nprop)
  /// arrays of node numbers for particular surface properties (prop_list_nod[i][j] = j-th node of surface with property proplist[i]
  long     **prop_list_nod; 

  /// arrays of number of nodes for particular surface properties (prop_list_nod_num[i] = number of nodes with surface property proplist[i] i.e. length of prop_list_nod[i])
  long     *prop_list_nod_num;

  /// arrays of element numbers for particular surface properties (prop_list_elem[i][j] = j-th element of surface with property proplist[i])
  long     **prop_list_elem;

  /// arrays of number of elements for particular surface properties (prop_list_elem_num[i] = number of elements with surface property proplist[i] i.e. length of prop_list_elem[i])
  long     *prop_list_elem_num;

  /// arrays of element surface id for particular surface properties and elements (prop_list_elem_sfid[i][j][k] = k-th surface id of j-th element with surface property proplist[i])
  long     ***prop_list_elem_sfid;

  /// arrays of number of element surfaces with particular surface properties for particular elements (prop_list_elem_sfid_num[i][j] = number of surfaces of j-th element with surface property proplist[i] i.e. length of prop_list_elem_sfid[i][j])
  long     **prop_list_elem_sfid_num;

  sfaces();
  sfaces(long n); ///< initializes to array of faces with size nfc
  ~sfaces();
  void copyto(sfaces &dest);  ///< method copies the given instance into the argument instance dest

  void read(XFILE *in);  ///< method reads the data about surfaces from the file
  void print(FILE *out, long *lnn); ///< method prints the data about surfaces to the file
  long search_list_newprop(long prop, long *propid, long n); ///< the function tests whether the given property number is new comparing to the array of existing property numbers 
  void gen_list_surf(long nn, long ne, selement *elems, long *nadjelnod, long **adjelnod); ///< method generates the lists of properties, nodes and elements for particular properties of the faces
//  void gen_list_surf(long nn, long ne, selement *elems);  ///< method generates the lists of properties, nodes and elements for particular properties of the faces
  void gen_list_node(long nn);  ///< method generates the lists of nodes for particular properties of the faces
  void gen_list_elem(long ne, selement *elems, long *nadjelnod, long **adjelnod); ///< method generates the lists elements for particular properties of the faces
//  void gen_list_elem(long ne, selement *elems);  ///< method generates the lists elements for particular properties of the faces
};



/**
  This class stores all information about topology of given problem but only 
  for auxiliary purposes in preprocessors especially.

  Created by Tomas Koudelka, 2001
*/
class siftop
{
  public:
   siftop (void); 
   siftop (long inn, long ine, gtypel te);
   explicit siftop(const siftop &src);
   ~siftop (void);

   /** method makes partial copy of the given siftop instance but allocates additional space for nodes, elements, 
       edge and surface property objects */
   void partial_copy(siftop &top, long tndn, long tnifel, long num_dnedgprop, long num_dnsurfprop) const;
  
   /// method reads topology in JKTK format from file
   long read(XFILE *in, long rgnn, long rebp);

   /// method imports topology in T3D format
   long import_t3d(XFILE *in, long paral);

   /// method imports section of nodes from file in T3D format
   long import_t3d_nodes (XFILE *in, long paral, long *outmode);

   /// method imports section of elements from file in T3D format
   long import_t3d_elements (XFILE *in, long *outmode);

   /// method imports topology in GiD format
   void import_gid(XFILE *in);

   /// method copies topological data to gtopology structure
   void exporttop (gtopology *top);

   /// method exports topology to a text file in JKTK format
   void print (FILE *out) const;

   /// prints nodes with the shift indices
   void shift_print_nodes (FILE *out,long shift) const;
   
   /// print elements with the shift indices
   void shift_print_elements (FILE *out,long shiftnn,long shiftne) const;

   /// method exports siftop mesh in GiD format to a text file
   void export_gid_mesh(FILE *out, long idn1=1, long ide1=1) const;

   /// method writes nodes in GiD format to a text file
   void write_gid_nodes(FILE *out, long idn1=1) const;

   /// method writes one element in GiD format to a text file
   void write_gid_element(FILE *out, long i, long idn1=1, long ide1=1) const;

   /// method generates lists of adjacent elements to nodes   
   void gen_adjelnod(void); 

   /// method gets ndof of nodes with property prop of entity ent
   long get_ndofn(long prop, gentity ent, long *ndofn, long *setnodes) const;

   /// method gets ndof of nodes with property prop of entity ent
   long get_ndofn(long prop, gentity ent, long *ndofn, ivector &setnodes) const;
   
   /// returns nodes laying on the entities ent with property prop
   long get_propent_nodes(long prop, gentity ent, long *setnodes) const;

   /// returns nodes laying on the entity ent with property prop
   long get_propent_nodes(long prop, gentity ent, ivector &setnodes) const;

   /// returns nodes laying on the entity ent with property prop, output will be in the compact form
   long get_propent_nodes_compact(long prop, gentity ent, ivector &setnodes) const;

   /// returns elements involving the entity ent with property prop
   long get_propent_elems(long prop, gentity ent, ivector &setelems, long accum=0) const;
   
   /// returns elements involving the entity ent with property prop, output will be in the compact form
   long get_propent_elems_compact(long prop, gentity ent, ivector &setelems) const;

  /// function reorders mesh from t3d format into sifel format (3D)
   void t3d2sifelformat();

   /// returns 2D nodal coordinates for the given element
   void give_node_coord3d(vector &x, vector &y, vector &z, long eid) const;

   /// returns 3D nodal coordinates for the given element
   void give_node_coord2d(vector &x, vector &y, long eid) const;

   /// returns dimension of the given element
   long give_dimension(long eid)const ;

   /// returns array of element nodes
   void give_elemnods(long eid, ivector &enod) const;

   /// returns the number of closest node and the given minimum distance 
   long give_closest_node_coord(long eid, double px, double py, double pz, double &dmin) const;

   /// returns coordinates of centre of gravity of the given element
   long centroid(long eid, vector &coord) const;

   /// returns the square of maximum distance between the given point and given element nodes
   double max_sqrdist_nod_pt(long eid, vector &pt) const;

   /// returns array of element nodes with the given xi natrual coordinate and corresponding surface id
   long give_enod_xicoord(long eid, double xi, ivector &nod) const;

   /// returns array of element nodes with the given eta natrual coordinate and corresponding surface id
   long give_enod_etacoord(long eid, double eta, ivector &nod) const;

   /// returns array of element nodes with the given eta natrual coordinate and corresponding surface id
   long give_enod_zetacoord(long eid, double zeta, ivector &nod) const;

   /// transforms point given in element natural coordinates to the local natural coordinates of the given element surface
   void transform_natcoord_surf(long eid, long sid, const vector &enc, vector &lnc) const;

   /// transforms point given in element natural coordinates to the local natural coordinates of the given element edge
   void transform_natcoord_edg(long eid, long edid, const vector &enc, vector &lnc) const;

   /// returns entity type and general element type of entity on the given element, entity is detected according to number of nodes
   gtypel give_elem_ent(long eid, long entnn, gentity &ent) const;
  
   selement *elements; ///< pointer to array with elements
   snode    *nodes;    ///< pointer to array with nodes
   sedges   *edges;    ///< pointer to structure with additional edges 
   sfaces   *surfaces; ///< pointer to structure with additional faces
   long      nn;       ///< number of nodes
   long      ne;       ///< number of elements
   long      npet[all_element_types-1];  ///< number of particular element types (used for t3d format)
   long      degree;   ///< degree of approximation polynom of the elements (used for t3d format)
   long     *gnn;      ///< pointer to array with global node numbers. This is support for paralell computing
   long     *nadjelnod;///< pointer to array of number of adjacent elements to nodes (nadjelnod[i] = number of adjacent elements of i-th node)
   long    **adjelnod; ///< array of pointers to arrays of adjacent elements id to particular nodes (adjelnod[i][j] = id o j-th adjacent element to i-th node)
       
   //  parameters for BOSS preconditioners
   ///  type of mesh description
   long meshtype;
   ///  number of subdomains/aggregates
   long nsd;
   /// number of nodes on subdomains/aggregates
   long *nnsd;
};

/// creates list of edges and elements adjacent to these edges for selected list of elemenst
//void generate_edge_list(std::vector<std::vector<long>> &edglst,
//                        std::vector<std::vector<long>> &edgelem,
//                        std::vector<std::vector<long>> &elemedg,
//                        const std::vector<long> &elems, const selement *pelems, const snode *pnodes);
#endif
