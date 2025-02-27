#include "genifcelem.h"
#include "galias.h"
#include "iotools.h"
#include "vector.h"
#include "siftop.h"
#include "selection.h"
#include "boundbox.h"
#include "mathem.h"

#include <vector>
#include <map>
#include <algorithm>




double cmp_with_tol::tolerance = 0.0;

int main(int argc, char *argv[])
{
  if (argc < 3){
    fprintf(stderr, "\nError - missing command line argument.\n"
            "\n\nUse:  '%s input_file_name output_file_name'\n", argv[0]);
    return 1;
  }
  if (argc > 3){
    fprintf(stderr, "\nError - too many command line arguments.\n"
            "\n\nUse:  '%s input_file_name output_file_name'\n", argv[0]);
    return 1;
  }

  XFILE *in=xfopen(argv[1], "rt");
  if (in == NULL){
    fprintf(stderr, "\nError - cannot open input file '%s'.\n", argv[1]);
    return 1;
  }
  in->kwdmode = sequent_mode;

  siftop top_orig; // original topology
  genopt options;
  long ncmd = 0;
  long i, mn;
  char *aptr = NULL;
  bool ok;

  // read options for the given problem and load original topology
  // parse the input file and save locations of the generator commands
  long ret = process_input_file(in, options, top_orig, ncmd);
  if (ret)  return 2;

  // copy original topology to the new one
  siftop top_new;

  // initialize bounding box of the topology
  boundbox topbb(top_orig);

  lvector dn_lst;  // list of node numbers to be doubled in the actual command
  lvector ifc1_lst;
  lvector ifc2_lst;
  std::vector<ifcel_descript> ifcdel(ncmd); // list of descriptions of generated interface elements
  std::set<boundary_ent_ifc> face_lst; // list of identifiers {elem_id, edg/surf_id, ifceld_ptr} where the interface element will be connected

  
  xf_setsec(in, bsec_ifc_str[0]);
  xf_resetsec(in);

  fprintf(stdout, "Reading of %ld command(s) for generation of interface elements:\n", ncmd);
  for (i=0; i< ncmd; i++){
    in->kwdmode = sect_mode_fwd;
    getkwd_sect(in, NULL, aptr, "gen_ifc", 0);
    in->kwdmode = sect_mode_seq;
    fprintf(stdout, " - command at line %ld ...", in->line);
    ok = true;

    // search for doubled node numbers
    ret = doubled_nodes(in, top_orig, topbb, dn_lst, ifc1_lst, ifc2_lst,  mn);
    if (ret)  return 3;
    if (mn){
      fprintf(stdout, "\n   * there were %ld doubled nodes.\n", mn);
      ok = false;
    }
    
    // reads command record dealing with element set to be searched for boundary entities connected
    // with the interface nodes
    ret = search_elem_connect2ifc(in, top_orig, dn_lst, ifc1_lst, ifc2_lst, face_lst, ifcdel[i]);
    if (ret)  return 4;

    if (ok)
      fprintf(stdout, " O.K.\n");
    else
      fprintf(stdout, "... O.K.\n");
  }
  xfclose(in);

  // generate interface elements if required and store them in the new topology
  gen_iface_elems(top_orig, face_lst, ifc1_lst, ifc2_lst, options, top_new);

  fprintf(stdout, "Writing of the output topology file ...");
  FILE *out = fopen(argv[2], "wt");
  if (out == NULL){
    fprintf(stderr, "\nError - cannot open output file '%s'.\n", argv[2]);
    return 5;
  }
  // export the new topology to the output file
  top_new.print(out);
  fprintf(stdout, "O.K.\n");
  fclose(out);

  return 0;
}



/**
  The function reads the original topology file where some of nodes should be doubled
  and interface elements may be generated optiionally. Additionally, the options
  of the interface elements are read and number of commands for node doubling and
  interface element generation are detected. 

  @param[in, out]  in - pointer to the opened input text file
  @param[out] options - 
  @param[out] top - structure for storage of the original topology
  @param[out] ncmd - number of commmands detected
  
  @retval 0 - on sucess
  @retval 1 - reading of topology failed
  @retval 2 - error in the reading of section with generation commands
  @retval 3 - no commands found

  Created by TKo, 01.2024
*/
long process_input_file(XFILE *in, genopt &options, siftop &top, long &ncmd)
{
  char fname[1025];  
  xfscanf(in, "%1024a", fname);
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &options.meshfmt);
  xfscanf(in, "%k%m", "edge_numbering", &answertype_kwdset, &options.edgnum);

  // reading of the file with original topology
  fprintf(stdout, "\n\nReading of mesh topology . . .");
  XFILE *intop = xfopen(fname, "rt");
  if (intop == NULL){
    print_err("cannot open topology file '%s'.", __FILE__, __LINE__, __func__, fname);
    return 1;
  }
  long ret;
  switch (options.meshfmt)
  {
    case t3d:
      top.import_t3d(intop, 0);
      break;
    case sifel:
      ret = top.read(intop, 0, options.edgnum);
      if (ret)  return 1;
      break;
    default:
      print_err("unknown mesh format is required.", __FILE__, __LINE__, __func__);
      return 1;
  }
  fprintf(stdout, " O.K.\n");
  xfclose(intop);
  if (top.nadjelnod == NULL)
    top.gen_adjelnod();
  
  in->kwdmode = sect_mode_seq;
  // detect section with commands for generation of interface nodes/elements
  xfdetect_sect(in, bsec_ifc_kwdset, esec_ifc_kwdset);
  // index of sections must be created
  if (in->index_created != 1){
    print_err("index of sections has not beeen created.", __FILE__, __LINE__, __func__);
    return 2;
  }
  // section with commands for interface node/element generation must be detected
  ret = xf_setsec(in, bsec_ifc_str[begsec_ifc_cmd]);
  if (ret){
    print_err("section with commands for generation of interface elements has not been detected.",
              __FILE__, __LINE__, __func__);
    return 2;
  }

  char *aptr = NULL;
  ncmd = getkwd_sect(in, NULL, aptr, "gen_ifc", 1);
  if (ncmd == 0){
    print_err("no commands for interface element generation were found.", __FILE__, __LINE__, __func__);
    return 3;
  }
  
  return 0;
}



/**
  The function searchs for doubled nodes in the actual generation command of the input file in. The doubled nodes are 
  searched in the given selection of nodes.

  @param[in, out] in   - pointer to the opened input text file, the file position must be set to the beginning
                         of the actual command for generation of interface nodes/elements.
  @param[in] top       - structure with the original topology from which the doubled nodes should be generated.
  @param[in] topbb     - bounding box of the topology coordinates
  @param[out] dn_lst   - list of node numbers to be doubled by the actual command, dn_lst = ifc1_lst + ifc2_lst
  @param[out] ifc1_lst - one half of interface node list, ifc1_lst \in dn_lst
  @param[out] ifc2_lst - the other half of interface node list, ifc2_lst \in dn_lst
  @param[out] mn - number of multiply defined doubled nodes detected in the function

  @return - list of interface nodes to be doubled is returned in the argument dn_lst
  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long doubled_nodes(XFILE *in, const siftop &top, const boundbox &topbb, lvector &dn_lst, lvector &ifc1_lst, lvector &ifc2_lst, long &mn)
{
  long i, j, k;
  sel selnod;
  long line = in->line;
  long col = in->col;
  lvector sel_nodid;
  std::map<double, long, cmp_with_tol> nodmap;
  ivector dncheck(top.nn);
  
  dn_lst.clear();

  xfscanf(in, "%k %k", "gen_ifc", "node_selection");
  selnod.read(in);
  // conversion of property selection to list of nodes or ranges of nodes
  if (selnod.st == sel_prop)
    selnod.conv_selprop(&top, gnod);  
  switch (selnod.st){
    case sel_all:
      for(i=0; i<top.nn; i++)  sel_nodid.push_back(i);     
      break;
    case sel_range:
      for (i=0; i<selnod.n; i++){
        for(j=0, k=selnod.id1[i]; j<selnod.ncomp[i]; j++, k++)
          sel_nodid.push_back(k);
      }
      break;
    case sel_list:
      for (i=0; i<selnod.n; i++){
        j = selnod.id1[i];
        sel_nodid.push_back(j);
      }
      break;    
    default:
      print_err("unknown type of interface node selection is required (%s)\n"
                "Only 'sel_range', 'sel_list' or 'sel_prop' selection types are allowed.",
                __FILE__, __LINE__, __func__, seltype_kwdset.get_str(selnod.st));
      return 1;
  }
  if (sel_nodid.size() == 0){
    print_warning("no interface node were selected by the gen_ifc command  (line=%ld, col=%ld)",
                  __FILE__, __LINE__, __func__, line, col);
  }

  xfscanf(in, "%k %le", "doubled_nodes_dist_tol", &cmp_with_tol::tolerance);
  
  for(i=0; i<long(sel_nodid.size()); i++){
    j = sel_nodid[i];
    double dist = sqrt(sqr(top.nodes[j].x - topbb.xmin) + sqr(top.nodes[j].y - topbb.ymin) + sqr(top.nodes[j].z - topbb.zmin));
    auto ret = nodmap.insert(std::make_pair(dist, j));
    if (ret.second == false){
      ifc1_lst.push_back(ret.first->second);
      ++dncheck[ret.first->second];
      ifc2_lst.push_back(j);
      ++dncheck[j];
    }
  }
  for (i=0; i<top.nn; i++){
    if (dncheck[i] > 1){
      print_warning("node %ld with coordinates (%le, %le, %le) was found %ld times\n"
                    "Nodes are allowed to be found 2 times at maximum.", __FILE__, __LINE__, __func__,
                    i+1, top.nodes[i].x, top.nodes[i].y, top.nodes[i].z, dncheck[i]+1);
      abort();
    }
  }
  if (ifc1_lst.size() != ifc2_lst.size()){
    print_err("sorting of doubled interface nodes failed (nifc1=%ld, nifc2=%ld)",
              __FILE__, __LINE__, __func__, long(ifc1_lst.size()), long(ifc2_lst.size()));
    abort();
  }
  mn = ifc1_lst.size() + ifc2_lst.size();
  dn_lst.reserve(mn);
  dn_lst.insert(dn_lst.begin(), ifc1_lst.begin(), ifc1_lst.end());
  dn_lst.insert(dn_lst.end(), ifc2_lst.begin(), ifc2_lst.end());
  std::sort(dn_lst.begin(), dn_lst.end());
  return 0;  
}



/**
  The function reads part of the command for an interface generation dealing with the element set and interface 
  element generation setup. Additionally, it searches for elements connected to the interface nodes in the selected set of elements
  and saves the result in argument face_lst.

  @param[in, out]  in  - pointer to the opened input text file, the file position indicator must correspond the 
                         position inside the interface generation command at the beginning of the element set 
                         definition to be considered in the interface searching.  
  @param[in]  top - structure with the original topology from which the doubled nodes should be generated.
  @param[in]  dn_lst   - list of node numbers to be doubled by the actual command, dn_lst = ifc1+ifc2
  @param[in]  ifc1     - one half of list of interface nodes nodes ifc1 \in dn_lst
  @param[in]  ifc2     - other half of list of interface nodes nodes ifc2 \in dn_lst
  @param[out] face_lst - list of strutures containing: 
                          * identifier of found element on interface, 
                          * index of its boundary entity (edge or surface) which is defined by the interface nodes, 
                            i.e. the entity, which is fully connected to the interface nodes; 
                          * pointer to the description of the generated interface elements
  @param[out] ifc_desc - description of generated interface elements

  @retval 0 - on success
  @retval 1 - in the case of wrong selection of interface nodes

  Created by TKo, 01.2024
*/
long search_elem_connect2ifc(XFILE *in, const siftop &top, const lvector &dn_lst, const lvector &ifc1,
                             const lvector &ifc2, std::set<boundary_ent_ifc> &face_lst,
                             ifcel_descript &ifcdel)
{
  long i, j, k, l;
  sel selelem;
  long line = in->line, col = in->col;
  
  // read selection of element set where the interface nodes will be searched for
  xfscanf(in, "%k", "elem_selection");
  selelem.read(in);
  // read setup of generated interface elements
  ifcdel.read(in);

  ivector e_lst(top.ne);
  // conversion of property selection to list of elements or ranges of elements
  if (selelem.st == sel_prop)
    selelem.conv_selprop(&top, gelem);  
  switch (selelem.st){
    case sel_all:
      for(i=0; i<top.ne; i++)  e_lst[i] = 1;
      break;
    case sel_range:
      for (i=0; i<selelem.n; i++){
        for(j=0, k=selelem.id1[i]; j<selelem.ncomp[i]; j++, k++)
          e_lst[k] = 1;
      }
      break;
    case sel_list:
      for (i=0; i<selelem.n; i++){
        j = selelem.id1[i];
        e_lst[j] = 1;
      }
      break;
    default:
      print_err("unknown type of interface node selection is required (%s)\n"
                "Only 'sel_range', 'sel_list' or 'sel_prop' selection types are allowed.",
                __FILE__, __LINE__, __func__, seltype_kwdset.get_str(selelem.st));
      return 1;
  }
  if (dn_lst.size() == 0){
    print_warning("no interface node were selected by the gen_ifc command  (line=%ld, col=%ld)",
                  __FILE__, __LINE__, __func__, line, col);
  }

  // create list of elements connected to interface nodes
  std::set<long> adjelnodifc_lst;
  for (i=0; i<(long)dn_lst.size(); i++){
    k = dn_lst[i];
    for(j=0; j<top.nadjelnod[k]; j++){
      l = top.adjelnod[k][j];
      adjelnodifc_lst.insert(l);
    }
  }

  // find boundary entity of elements on interface (lookup in adjelnodifc_lst) and
  // search for opposite interface nodes to the found ones on the boundary entity
  long nent;
  ivector entid_lst;
  ivector entnod2;
  for (long eid : adjelnodifc_lst){
    if (e_lst[eid]){
      if (top.give_dimension(eid) <= 2){
        reallocv(RSTCKIVEC(top.elements[eid].ned, entid_lst));
        nent = top.elements[eid].compare_edg(&dn_lst.front(), dn_lst.size(), entid_lst);
      }
      else{
        reallocv(RSTCKIVEC(top.elements[eid].nsurf, entid_lst));
        nent = top.elements[eid].compare_surf(&dn_lst.front(), dn_lst.size(), entid_lst);
      }
        
      if (nent){
        // element edges or surfaces were found to be defined by interface nodes
        for(i=0; i<entid_lst.n; i++){
          search_ifc_oppo_ent_nodes(eid, entid_lst[i], ifc1, ifc2, top, entnod2);
          face_lst.insert(boundary_ent_ifc(eid, entid_lst[i], entnod2.a, entnod2.n, &ifcdel)); // store the found element id and interface entity id
          entnod2.a = NULL;
          entnod2.n = entnod2.size = 0;
        }
      }
    }
  }
  return 0;
}


/**
  The function searches for opposite interface element to the given eid-th interface element.
  The eid-th element is connected to the interface by entid-th element boundary entity (edge or surface) 

  @param[in] eid - the given element identifier,
  @param[in] entid - identfier of the eid-th element boundary entity (edge or surface)
  @param[in] ifc1 - one half of doubled nodes defining the given interface, ifc1 + ifc2 = all selected doubled nodes
  @param[in] ifc2 - the other half of doubled nodes defining the given interface, ifc1 + ifc2 = all selected doubled nodes
  @param[in] top - original topology with initialize arrays adjelnod and nadjelnod
  @param[out] entnod2 - the opposite interface nodes to ones that defines entid-th 
                        boundary entity (edge or surface) at eid-th element
 
  @return The function returns the opposite doubled nodes in the argument entnod2.

  Created by TKo, 03.2024
*/
void search_ifc_oppo_ent_nodes(long eid, long entid, const lvector &ifc1,
                               const lvector &ifc2, const siftop &top, ivector &entnod2)
{
  long i;
  ivector rent_nod;
  ivector entnod;
  if (top.give_dimension(eid) <= 2)
    makerefv(rent_nod, top.elements[eid].edgenod[entid], top.elements[eid].nned[entid]);
  else
    makerefv(rent_nod, top.elements[eid].surfnod[entid], top.elements[eid].nnsurf[entid]);
  reallocv(RSTCKIVEC(rent_nod.n, entnod));
  reallocv(rent_nod.n, entnod2); // this must not be allocated on the stack, it represents the output argument
  for (i=0; i<rent_nod.n; i++){
    entnod[i] = top.elements[eid].nodes[rent_nod[i]];
    auto anod = std::find(ifc1.begin(), ifc1.end(), entnod[i]);
    bool  found = false;
    if (anod != ifc1.end()){
      entnod2[i] = ifc2[anod - ifc1.begin()];
      found = true;
    }
    if (!found){
      anod = std::find(ifc2.begin(), ifc2.end(), entnod[i]);
      if (anod != ifc2.end()){
        entnod2[i] = ifc1[anod - ifc2.begin()];
        found = true;
      }
    }
    if (!found){
      print_err("corresponding doubled interface node was not found on element %ld (entid=%ld, ent_nodid=%ld)",
                __FILE__, __LINE__, __func__, eid+1, i+1, entnod[i]+1);
      abort();
    }
  }
}



/**
  The function generates new doubled nodes and interface elements if the new interface elements are required.

  @param[in] top_orig     - structure with the original topology from which the doubled nodes should be generated.
  @param[in] face_lst     - list of strutures containing: 
                            * identifier of found element on interface, 
                            * index of its boundary entity (edge or surface) which is defined by the interface nodes, 
                              i.e. the entity, which is fully connected to the interface nodes; 
                            * array of opposite nodes to the given boundary entity 
                            * dimension of the above array
                            * pointer to the description of the generated interface elements
  @param[in] gopt         - structure with the options for mesh format used
  @param[out] top_new     - structure with the new topology where reconnection of elements on interface to doubled nodes will be performed

  @return The function does not return anything.
  
  Created by TKo, 01.2024
*/
void gen_iface_elems(const siftop &top_orig, const std::set<boundary_ent_ifc> &face_lst, const lvector &ifc1,
                     const lvector &ifc2, const genopt &gopt, siftop &top_new)
{
  // reconnect of found elements on interface to new doubled nodes
  long i, eid, entid, nnent, dim;
  ivector rentnod, rentnod2;
  ivector entnod, entnod2;

  // count number of generated interface elements
  long tnifel = 0;
  for (auto eld : face_lst){
    if (eld.ifceld_ptr->et != ifc_noel)
      tnifel++;
  }

  // make copy of topology with extended arrays of nodes and elements and actualize extended list of properties
  // on edges and surfaces if they are defined
  top_orig.partial_copy(top_new, 0, tnifel, 0, 0);
 
  // generate interface elements if required
  long k=top_orig.ne;
  for (auto eld : face_lst){    
    if (eld.ifceld_ptr->et == ifc_noel)
      continue;
    eid = eld.eid;
    entid = eld.entid;
    dim = top_orig.give_dimension(eid);
    if (dim <= 2){
      nnent = top_orig.elements[eid].nned[entid];
      makerefv(rentnod, top_orig.elements[eid].edgenod[entid], nnent);
    }
    if (dim == 3){
      nnent = top_orig.elements[eid].nnsurf[entid];
      makerefv(rentnod, top_orig.elements[eid].surfnod[entid], nnent);
    }
    reallocv(RSTCKIVEC(rentnod.n, entnod));
    makerefv(entnod2, eld.entnod2, nnent);
    for(i=0; i<rentnod.n; i++)
      entnod[i]  = top_orig.elements[eid].nodes[rentnod[i]];
    make_ifc_elem(top_new, eid, entnod, entnod2, gopt, *eld.ifceld_ptr, k);
    k++;
  }
}



/**
  The function makes new interface element with identifier k in the topology top. The new interface element
  nodes will be defined by with the help of nodes entnod such that the half of nodes will be associated with 
  with etnodes directly and the secnod half will be associated to nodes entnodes renumbered with the help of 
  gdn_lst. If the ifceld.orient == ifc_elor_low then the first half of nodes will be entnodes directly while 
  the renumbered nodes will be used for the second half. If the ifceld.orient == ifc_elor_high then the
  oposite node assigment is used. On the given half of nodes, the node assignment start with the ifceld.nodid-th 
  node from the array entnodes. The type of generated element is either prescribed by ifceld.et or detected 
  automatically according to neighbour element eid if ifceld.et == ifc_automatic.

  @param[in, out] top - the topology where the new interface element will be stored as the k-th element,
                        the array top.elements must be allocated sufficiently to hold the new element
  @param[in] eid - element identifier which the new interface element will be connected to. 
  @param[in] entnod  - array of nodes on the boundary entity (edge or surface) of the eid-th element which 
                       the new interface element will be connected to.
  @param[in] entnod2 - array of nodes on the boundary entity (edge or surface) of the element opposite to eid-th element 
                       which the new interface element will be connected to.
  @param[in] gopt - structure with the options for mesh format used
  @param[in] ifceld - the generation setup of the new interface element

  @return The function does not return anything.

  Created by TKo, 01.2024
*/
void make_ifc_elem(siftop &top, long eid, const ivector &entnod, const ivector &entnod2,
                   const genopt &gopt, const ifcel_descript &ifceld, long k)
{
  ifctypel ifcet = ifceld.et;
  if (ifcet == ifc_automatic){
    gtypel et = top.elements[eid].type;
    switch (et){
      case isolinear1d:
      case trianglelinear:
      case isolinear2d:
        ifcet = ifc_linquad;
        break;
      case isoquadratic1d:
      case trianglequadratic:
      case isoquadratic2d:
        ifcet = ifc_quadquad;
        break;
      case isolinear3d:
        ifcet = ifc_linhexa;
        break;
      case isoquadratic3d:
        ifcet = ifc_quadhexa;
        break;
      case tetrahedronlinear:
        ifcet = ifc_linwedge;
        break;
      case tetrahedronquadratic:
        ifcet = ifc_quadwedge;
        break;        
      case wedgelinear:
        if (entnod.n == 4)
          ifcet = ifc_linhexa;
        else
          ifcet = ifc_linwedge;
        break;
      case wedgequadratic:
        if (entnod.n == 8)
          ifcet = ifc_quadhexa;
        else
          ifcet = ifc_quadwedge;
        break;
      default:
        print_err("unknown/unimplemented element type '%s' is required for interface element generation.",
                  __FILE__, __LINE__,__func__, gtypel_kwdset.get_str(et));
        abort();
    }
  }
  selement &ifcel = top.elements[k];
  ifcel.type = noel;  
  switch(ifcet){
    case ifc_linquad:
      ifcel.type = isolinear2d;
      break;
    case ifc_linhexa:
      ifcel.type = isolinear3d;
      break;
    case ifc_linwedge:
      ifcel.type = wedgelinear;
      break;
    default:
      print_err("unknown type of interface element '%s' is required", __FILE__, __LINE__, __func__, ifctypel_kwdset.get_str(ifcet));
      abort();
  }
  ifcel.alloc(gopt.edgnum);
  ifcel.prop = ifceld.prop;
  
  long n = entnod.n;
  if (ifceld.orient == ifc_elor_low){
    for (long i=0, j=ifceld.nodid+n-1; i<n; i++, j--){
      ifcel.nodes[i] = entnod[j%n];
    }
    if (ifcet == ifc_linquad){
      for (long i=n, j=ifceld.nodid; i<ifcel.nne; i++, j++){
        ifcel.nodes[i] = entnod2[j%n];
      }
    }
    else{
      for (long i=n, j=ifceld.nodid+n-1; i<ifcel.nne; i++, j--){
        ifcel.nodes[i] = entnod2[j%n];
      }
    }      
  }
  if (ifceld.orient == ifc_elor_high){
    for (long i=0, j=ifceld.nodid; i<n; i++, j++){
      ifcel.nodes[i] = entnod[j%n];
    }
    if (ifcet == ifc_linquad){
      for (long i=n, j=ifceld.nodid+n-1; i<ifcel.nne; i++, j--){
        ifcel.nodes[i] = entnod2[j%n];
      }
    }
    else{
      for (long i=n, j=ifceld.nodid; i<ifcel.nne; i++, j++){
        ifcel.nodes[i] = entnod2[j%n];
      }
    }
  }
}



/**
  The operator compares two objects of boundary_ent_ifc class. The operator is required by the container
  std::set<boundary_ent_ifc>. The comparison also checks for faulty input of interface elements defined 
  with different setup for the same boundary entity (edge or surface) of the same element connected to 
  the interface nodes.

  @param[in] lhs - left hand side of the operator
  @param[in] rhs - right hand side of the operator

  @retval true if lhs < rhs
  @retval false otherwise

  Created by TKo, 01.2024
*/
bool operator < (const boundary_ent_ifc &lhs, const boundary_ent_ifc &rhs)
{
  if (lhs.eid < rhs.eid)
    return true;
  else{
    if (lhs.eid == rhs.eid){
      if (lhs.entid < rhs.entid)
        return true;
      else{
        if (lhs.entid == rhs.entid){
          if (lhs.nentnod2 != rhs.nentnod2){
            // lhs and rhs are the same elements where the interface element is required to be connected to the same boundary entity(edge or surface) and
            // there is different number of nodes of the interface element boundary entities => failure
            // it must not appear on the same entity, apparently, there is something wrong in the algorithm
            print_err("failure of interface element definition on element %ld, entid=%ld:\n"
                      " there are two different numbers of nodes of interface element boundary entites.\n"
                      " lhs nentnod2=%ld\n"
                      " rhs nentnod2=%ld\n",
                      __FILE__, __LINE__, __func__, lhs.eid+1, lhs.entid+1, lhs.nentnod2, rhs.nentnod2);
            abort();
          }
          int ret = memcmp(lhs.entnod2, rhs.entnod2, lhs.nentnod2*sizeof(*lhs.entnod2));
          if (ret < 0)
            return true;
          if ((ret == 0) && !(*lhs.ifceld_ptr == *lhs.ifceld_ptr)){
            // lhs and rhs are the same elements where the interface element is required to be connected to the same boundary entity(edge or surface) and
            // there is different setup of the generated interface element => failure
            // different interface element setup must not appear on the same entity, apparently, there are wrong data in the input file
            print_err("failure of interface element definition on element %ld:\n"
                      " there are two different setup of interface generation on %ld-th element boundary entity.\n"
                      " setup 1: ifceld_ptr=%p{et=%s, orient=%s, nodid=%ld}\n"
                      " setup 2: ifceld_ptr=%p{et=%s, orient=%s, nodid=%ld}\n",
                      __FILE__, __LINE__, __func__, lhs.eid+1, lhs.entid+1, lhs.ifceld_ptr, ifctypel_kwdset.get_str(lhs.ifceld_ptr->et),
                      ifc_elor_kwdset.get_str(lhs.ifceld_ptr->orient), lhs.ifceld_ptr->nodid, rhs.ifceld_ptr, ifctypel_kwdset.get_str(rhs.ifceld_ptr->et),
                      ifc_elor_kwdset.get_str(rhs.ifceld_ptr->orient), lhs.ifceld_ptr->nodid);
            abort();
          }
        }
      }
    }
  }
  return false;
}
