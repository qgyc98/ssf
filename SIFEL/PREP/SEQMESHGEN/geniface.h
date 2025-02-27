#ifndef GENIFACE_H
#define GENIFACE_H

#include "kwdset.h"
#include "galias.h"
#include "iotools.h"

#include <vector>
#include <set>

struct XFILE;
class siftop;
struct ivector;



enum bsec_ifc {begsec_ifc_cmd=0};
const enumstr bsec_ifc_str[] = {{"begsec_ifc_cmd",0}};
const kwdset bsec_ifc_kwdset(sizeof(bsec_ifc_str)/sizeof(*bsec_ifc_str), bsec_ifc_str);

enum esec_ifc {endsec_ifc_cmd=0};
const enumstr esec_ifc_str[] = {{"endsec_ifc_cmd",0}};
const kwdset esec_ifc_kwdset(sizeof(esec_ifc_str)/sizeof(*esec_ifc_str), esec_ifc_str);

enum ifctypel {ifc_noel=0, ifc_automatic=1, ifc_linquad=2, ifc_quadquad=3, ifc_linhexa=4, ifc_quadhexa=5, ifc_linwedge=6,
               ifc_quadwedge=7};
const enumstr ifctypel_str[] = {{"ifc_noel",0}, {"ifc_automatic",1}, {"ifc_linquad",2}, {"ifc_quadquad",3}, {"ifc_linhexa",4},
                                {"ifc_quadhexa",5}, {"ifc_linwedge",6}, {"ifc_quadwedge",7}};
const kwdset ifctypel_kwdset(sizeof(ifctypel_str)/sizeof(*ifctypel_str), ifctypel_str);

enum ifc_elor {ifc_elor_low=0, ifc_elor_high=1};
const enumstr ifc_elor_str[] = {{"ifc_elor_low",0}, {"ifc_elor_high",1}};
const kwdset ifc_elor_kwdset(sizeof(ifc_elor_str)/sizeof(*ifc_elor_str), ifc_elor_str);



struct genopt
{
  meshform meshfmt;
  answertype edgnum;
};



class ifcel_descript
{
 public:
  ifcel_descript() {
    et = ifc_noel; orient = ifc_elor_low; nodid = prop = 0;
  };
  
  long read(XFILE *in) {
    xfscanf(in, "%k%m", "ifc_elem_type", &ifctypel_kwdset, &et);
    if (et != ifc_noel){
      xfscanf(in, "%k%m", "ifc_elem_orientation", &ifc_elor_kwdset, &orient);
      xfscanf(in, "%k%ld", "ifc_start_node_id", &nodid);
      --nodid;
      xfscanf(in, "%k%ld", "ifc_vol_prop", &prop);
    }
    return 0;
  };

  bool operator == (const ifcel_descript &rhs)
                   { return (et == rhs.et) && (orient == rhs.orient) && (nodid == rhs.nodid) && (prop == rhs.prop);};

  ifctypel et;  /// type of generated interface elements
  ifc_elor orient; /// orientation of generated interface elements between elements connected to interface
  long nodid;  /// node index of the boundary entity of elements on interface from which the the definition of intefrace elements starts
  long prop; /// region property number of the generated interface elements
};



struct boundary_ent_ifc
{
  long eid; // the first element identifier to be connected to an interface
  long entid; // boundary entity identifier (edge or surface id) on the eid-the element where the interface was detected
  long eid2; // the second element identifier to be connected to an interface
  ifcel_descript *ifceld_ptr;
};



/// read options for the given problem, load original topology and detect number of commands for the interface generation
long process_input_file(XFILE *in, genopt &options, siftop &top_orig, long &ncmd);

/// search for nodes to be doubled given by the actual command
long doubled_nodes(XFILE *in, const siftop &top, std::vector<long> &dn_lst, long &gdnn, ivector &gdn_lst, long &mn);

/// search for elements connected to the doubled nodes given in the actual command and search for surfaces connected with doubled nodes on these elements
long search_elem_connect2ifc(XFILE *in, const siftop &top, const std::vector<long> &dn_lst,
                             std::set<boundary_ent_ifc> &face_lst, ifcel_descript &ifcdel);

/// search adjacent element to the given one which is connected to the eid-th element boundary entity
long search_adj_ifc_elem(long eid, long entid, const ivector &e_lst, const siftop &top);

/// generate doubled interface nodes, reconnect elements in touch with interface and generate interface elements if required
void gen_iface_elems(const siftop &top_orig, const ivector &gdn_lst, const std::set<boundary_ent_ifc> &face_lst,
                     const genopt &gopt, siftop &top_new);

/// returns number edge property objects involving the interface nodes in the edge property definition
long num_doubled_nodes_edg_prop(const siftop &top, const ivector &gdn_lst);

/// returns number surface property objects involving the interface nodes in the surface property definition
long num_doubled_nodes_surf_prop(const siftop &top, const ivector &gdn_lst);

/// adds new instances of edge property objects defined with new doubled nodes
void doubled_nodes_add_edg_prop(siftop &top, long num_dnedgprop, const ivector &gdn_lst);

/// adds new instances of surface property objects defined with new doubled nodes to the given topology
void doubled_nodes_add_surf_prop(siftop &top, long num_dnsurfprop, const ivector &gdn_lst);

/// makes new interface element with id k
void make_ifc_elem(siftop &top, long eid, const ivector &entnod, const ivector &gdn_lst, const genopt &gopt, const ifcel_descript &ifceld, long k);

/// comparison operator allowing the use std::set with boundary_ent_ifc class
bool operator < (const boundary_ent_ifc &lhs, const boundary_ent_ifc &rhs);
#endif
