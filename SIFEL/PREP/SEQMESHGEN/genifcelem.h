#ifndef GENIFCELEM_H
#define GENIFCELEM_H

#include "kwdset.h"
#include "galias.h"
#include "iotools.h"

#include <vector>
#include <set>
#include <math.h>
#include <string.h>

struct XFILE;
class siftop;
struct ivector;
class boundbox;



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



struct cmp_with_tol {
  bool operator()(double a, double b) const {
    double dif = fabs(a-b);
    if (dif <= tolerance)
      return false;
    else
      return a < b;
  };
  static double tolerance;
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
  boundary_ent_ifc(long aeid, long aentid, long *aentnod2, long anentnod2,  ifcel_descript *aifceld_ptr)
  {eid= aeid; entid = aentid; entnod2 = aentnod2; nentnod2 = anentnod2; ifceld_ptr = aifceld_ptr;};

  boundary_ent_ifc(const boundary_ent_ifc &src)
  {eid= src.eid; entid = src.entid; entnod2 = new long[src.nentnod2]; memcpy(entnod2, src.entnod2, sizeof(*src.entnod2)*src.nentnod2); nentnod2 = src.nentnod2; ifceld_ptr = src.ifceld_ptr;};

    ~boundary_ent_ifc() {delete [] entnod2;};

  long eid; // the first element identifier to be connected to an interface
  long entid; // boundary entity identifier (edge or surface id) on the eid-th element where the interface was detected
  long *entnod2; // nodes on edge or surface of the second element where the interface was detected, dimension corresponds to number of nodes on entid-th boundary entity of eid-th element
  long nentnod2; // dimension of the array entnod2
  ifcel_descript *ifceld_ptr;
};

using lvector = std::vector<long>;

/// read options for the given problem, load original topology and detect number of commands for the interface generation
long process_input_file(XFILE *in, genopt &options, siftop &top_orig, long &ncmd);

/// search for nodes to be doubled given by the actual command
long doubled_nodes(XFILE *in, const siftop &top, const boundbox &topbb, lvector &dn_lst, lvector &ifc1_lst, lvector &ifc2_lst, long &mn);

/// search for elements connected to the doubled nodes given in the actual command and search for surfaces connected with doubled nodes on these elements
long search_elem_connect2ifc(XFILE *in, const siftop &top, const lvector &dn_lst, const lvector &ifc1,
                             const lvector &ifc2, std::set<boundary_ent_ifc> &face_lst,
                             ifcel_descript &ifcdel);

/// search adjacent element to the given one which is connected to the eid-th element boundary entity
void search_ifc_oppo_ent_nodes(long eid, long entid, const lvector &ifc1,
                               const lvector &ifc2, const siftop &top, ivector &entnod2);

/// generate doubled interface nodes, reconnect elements in touch with interface and generate interface elements if required
void gen_iface_elems(const siftop &top_orig, const std::set<boundary_ent_ifc> &face_lst, const lvector &ifc1,
                     const lvector &ifc2, const genopt &gopt, siftop &top_new);

/// makes new interface element with id k
void make_ifc_elem(siftop &top, long eid, const ivector &entnod, const ivector &entnod2, const genopt &gopt, const ifcel_descript &ifceld, long k);

/// comparison operator allowing the use std::set with boundary_ent_ifc class
bool operator < (const boundary_ent_ifc &lhs, const boundary_ent_ifc &rhs);
#endif
