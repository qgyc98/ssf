#include "geniface.h"
#include "galias.h"
#include "iotools.h"
#include "vector.h"
#include "siftop.h"
#include "selection.h"

#include <vector>



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
  long gdnn = 0;
  bool ok;

  // read options for the given problem and load original topology
  // parse the input file and save locations of the generator commands
  long ret = process_input_file(in, options, top_orig, ncmd);
  if (ret)  return 2;

  // copy original topology to the new one
  siftop top_new;
  // generate lists of adjacent elements to nodes
  top_orig.gen_adjelnod();

  ivector gdn_lst(top_orig.nn); // global list of nodes to be doubled gdn_lst[i] > 0 -> i-th node should be doubled
  std::vector<long> dn_lst;  // list of node numbers to be doubled in the actual command
  std::vector<ifcel_descript> ifcdel(ncmd); // list of descriptions of generated interface elements
  std::set<boundary_ent_ifc> face_lst; // list of identifiers {elem_id, edg/surf_id, ifceld_ptr} where the interface element will be connected

  
  xf_setsec(in, bsec_ifc_str[0]);
  xf_resetsec(in);

  // set initial index value for new doubled nodes
  gdnn = top_orig.nn;
  fprintf(stdout, "Reading of %ld command(s) for generation of interface nodes/elements:\n", ncmd);
  for (i=0; i< ncmd; i++){
    in->kwdmode = sect_mode_fwd;
    getkwd_sect(in, NULL, aptr, "gen_ifc", 0);
    in->kwdmode = sect_mode_seq;
    fprintf(stdout, " - command at line %ld ...", in->line);
    ok = true;

    // search for node numbers to be doubled
    ret = doubled_nodes(in, top_orig, dn_lst, gdnn, gdn_lst, mn);
    if (ret)  return 3;
    if (mn){
      fprintf(stdout, "\n   * there were %ld nodes that has already been doubled.\n", mn);
      ok = false;
    }
    
    // reads command record dealing with element set to be searched for boundary entities connected
    // with the interface nodes
    ret = search_elem_connect2ifc(in, top_orig, dn_lst, face_lst, ifcdel[i]);
    if (ret)  return 4;

    if (ok)
      fprintf(stdout, " O.K.\n");
    else
      fprintf(stdout, "... O.K.\n");
  }
  xfclose(in);

  // generate interface elements if required and store them in the new topology
  gen_iface_elems(top_orig, gdn_lst, face_lst, options, top_new);

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
  The function searchs for nodes to be doubled in the actual generation command of the input file in.

  @param[in, out] in  - pointer to the opened input text file, the file position must be set to the beginning
                        of the actual command for generation of interface nodes/elements.
  @param[in] top      - structure with the original topology from which the doubled nodes should be generated.
  @param[out] dn_lst  - list of node numbers to be doubled by the actual command
  @param[in,out] gdnn - attained global node number for numbering of doubled nodes, at the beginning the number should be top_orig.nn,
                        C-style node numbering is used (starting from 0), the number of nodes in the original topology top_orig 
                        must be nonzero
  @param[in, out] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                            on the interface (to be doubled)
  @param[out] mn - number of multiply defined doubled nodes detected in the function

  @return - list of interface nodes to be doubled is returned in the argument dn_lst
  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long doubled_nodes(XFILE *in, const siftop &top, std::vector<long> &dn_lst, long &gdnn, ivector &gdn_lst, long &mn)
{
  long i, j, k;
  sel selnod;
  long line = in->line;
  long col = in->col;
  
  dn_lst.clear();

  xfscanf(in, "%k %k", "gen_ifc", "node_selection");
  selnod.read(in);
  // conversion of property selection to list of nodes or ranges of nodes
  if (selnod.st == sel_prop)
    selnod.conv_selprop(&top, gnod);  
  switch (selnod.st){
    case sel_range:
      for (i=0; i<selnod.n; i++){
        for(j=0, k=selnod.id1[i]; j<selnod.ncomp[i]; j++, k++)
          dn_lst.push_back(k);
      }
      break;
    case sel_list:
      for (i=0; i<selnod.n; i++){
        j = selnod.id1[i];
        dn_lst.push_back(j);
      }
      break;
    default:
      print_err("unknown type of interface node selection is required (%s)\n"
                "Only 'sel_range', 'sel_list' or 'sel_prop' selection types are allowed.",
                __FILE__, __LINE__, __func__, seltype_kwdset.get_str(selnod.st));
      return 1;
  }
  if (dn_lst.size() == 0){
    print_warning("no interface node were selected by the gen_ifc command  (line=%ld, col=%ld)",
                  __FILE__, __LINE__, __func__, line, col);
  }

  mn = 0;
  // generate identifiers of new doubled interface nodes  
  for (auto j : dn_lst){
    if (gdn_lst[j] == 0){
      gdn_lst[j] = gdnn;
      gdnn++;
    }
    else{
      //fprintf(stdout, "original interface node %ld has been already doubled (new node id=%ld)\n", j+1, gdn_lst[j]+1);
      mn++;
    }
  }
    
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
  @param[in]  dn_lst   - list of node numbers to be doubled by the actual command
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
long search_elem_connect2ifc(XFILE *in, const siftop &top, const std::vector<long> &dn_lst,
                             std::set<boundary_ent_ifc> &face_lst, ifcel_descript &ifcdel)
{
  long i, j, k;
  sel selelem;
  long line = in->line, col = in->col;
  
  // read selection of element set where the interface nodes will be searched for
  xfscanf(in, "%k", "elem_selection");
  selelem.read(in);
  // read setup of generated interface elements
  ifcdel.read(in);

  ivector e_lst(top.ne);
  // conversion of property selection to list of nodes or ranges of nodes
  if (selelem.st == sel_prop)
    selelem.conv_selprop(&top, gelem);  
  switch (selelem.st){
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
      adjelnodifc_lst.insert(top.adjelnod[k][j]);
    }
  }

  // reconnect of found elements on interface to new doubled nodes
  long nent;
  ivector entid_lst;
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
          long eid2 = search_adj_ifc_elem(eid, entid_lst[i], e_lst, top);
          if (eid2 >= 0)
            face_lst.insert({eid, entid_lst[i], eid2, &ifcdel}); // store the found element id and interface entity id
        }
      }
    }
  }
  return 0;
}


/**
  The function searches for adjacent element to the given eid-th element such that both
  are connected to the entid-th element boundary entity (edge or surface) of the eid-the element.

  @param[in] eid - the given element identifier,
  @param[in] entid - identfier of the eid-th element boundary entity (edge or surface)
  @param[in] e_lst - list of elements where adjacent element is NOT searched, 
                     dimension of the list is top.ne, e_lst[i] > 0 -> i-th element is not
                     considered for the searching
  @param[in] top - original topology with initialize arrays adjelnod and nadjelnod
 
  @return The function returns an identifier of the adjacent element if found or -1 if not.

  Created by TKo, 01.2024
*/
long search_adj_ifc_elem(long eid, long entid, const ivector &e_lst, const siftop &top)
{
  long nnent, i, j, nadj, adjeid, ret = -1;
  ivector rentnod;
  long dim = top.give_dimension(eid);
  if (dim <= 2){
    nnent = top.elements[eid].nned[entid];
    makerefv(rentnod, top.elements[eid].edgenod[entid], nnent);
  }
  else{
    nnent = top.elements[eid].nnsurf[entid];
    makerefv(rentnod, top.elements[eid].surfnod[entid], nnent);   
  }
  
  ivector entnod(ASTCKIVEC(nnent)); 
  std::set<long> adjel;
  for(i=0; i<nnent; i++){
    entnod[i] = top.elements[eid].nodes[rentnod[i]];
    nadj = top.nadjelnod[entnod[i]];
    for(j=0; j<nadj; j++){
      adjeid = top.adjelnod[entnod[i]][j];
      if (e_lst[adjeid] == 0)
        adjel.insert(adjeid);
    }
  }
  for(long eid2 : adjel){    
    if (dim <= 2)
      i = top.elements[eid2].compare_edg(entnod.a, nnent);
    else
      i = top.elements[eid2].compare_surf(entnod.a, nnent);
    if (i >= 0){
      ret = eid2;
      break;
    }
  }
  return ret;
}



/**
  The function generates new doubled nodes and interface elements if the new interface elements are required.

  @param[in]  top_orig    - structure with the original topology from which the doubled nodes should be generated.
  @param[in, out] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                            on the interface (to be doubled)
  @param[in] face_lst     - list of strutures containing: 
                            * identifier of found element on interface, 
                            * index of its boundary entity (edge or surface) which is defined by the interface nodes, 
                              i.e. the entity, which is fully connected to the interface nodes; 
                            * pointer to the description of the generated interface elements
  @param[in] gopt         - structure with the options for mesh format used
  @param[out] top_new     - structure with the new topology where reconnection of elements on interface to doubled nodes will be performed

  @return The function does not return anything.
  
  Created by TKo, 01.2024
*/
void gen_iface_elems(const siftop &top_orig, const ivector &gdn_lst, const std::set<boundary_ent_ifc> &face_lst, const genopt &gopt,
                     siftop &top_new)
{
  // reconnect of found elements on interface to new doubled nodes
  long i, eid, entid, nnent, dim;
  ivector rentnod;
  ivector entnod;

  // count the total number of doubled nodes
  long tndn = 0;
  for(i=0; i<top_orig.nn; i++){
    if (gdn_lst[i] > 0)  tndn++;
  }
  long tnifel = 0;
  for (auto eld : face_lst){
    if (eld.ifceld_ptr->et != ifc_noel)
      tnifel++;
  }

  long num_dnedgprop = 0;
  if (top_orig.edges)
    num_dnedgprop = num_doubled_nodes_edg_prop(top_orig, gdn_lst);
  long num_dnsurfprop = 0;
  if (top_orig.surfaces)
    num_dnsurfprop = num_doubled_nodes_surf_prop(top_orig, gdn_lst);

  // make copy of topology with extended arrays of nodes and elements and actualize extended list of properties
  // on edges and surfaces if they are defined
  top_orig.partial_copy(top_new, tndn, tnifel, num_dnedgprop, num_dnsurfprop);
  if (num_dnedgprop)
    doubled_nodes_add_edg_prop(top_new, num_dnedgprop, gdn_lst);
  if (num_dnsurfprop)
    doubled_nodes_add_surf_prop(top_new, num_dnsurfprop, gdn_lst);

  /*
  if (top_new.edges || top_new.surfaces)
    top_new.gen_adjelnod();
  if (top_new.edges)
    top_new.edges->gen_list_edg(top_new.nn, top_new.ne, top_new.elements, top_new.nadjelnod, top_new.adjelnod);
  if (top_new.surfaces)
    top_new.surfaces->gen_list_surf(top_new.nn, top_new.ne, top_new.elements, top_new.nadjelnod, top_new.adjelnod);
  */
  
  // generate new doubled nodes
  for(i=0; i<top_orig.nn; i++){
    if (gdn_lst[i] > 0)
      top_orig.nodes[i].copyto(top_new.nodes[gdn_lst[i]]);
  }

  // reconnect original elements connected to the interface nodes to new doubled nodes
  for (auto eld : face_lst){
    eid = eld.eid;
    entid = eld.entid;
    if (top_orig.give_dimension(eid) <= 2){
      nnent = top_orig.elements[eid].nned[entid];
      makerefv(rentnod, top_orig.elements[eid].edgenod[entid], nnent);
      
    }
    else{
      nnent = top_orig.elements[eid].nnsurf[entid];
      makerefv(rentnod, top_orig.elements[eid].surfnod[entid], nnent);   
    }
    // renumber interface nodes of found element in the new topology
    for (i=0; i<nnent; i++){
      long *enodes = top_orig.elements[eid].nodes;
      top_new.elements[eid].nodes[rentnod[i]] = gdn_lst[enodes[rentnod[i]]];
    }
  }

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
    for(i=0; i<rentnod.n; i++)
      entnod[i] = top_orig.elements[eid].nodes[rentnod[i]];
    make_ifc_elem(top_new, eid, entnod, gdn_lst, gopt, *eld.ifceld_ptr, k);
    k++;
  }
}



/**
  The function searches for the doubled nodes given in the gdn_lst in the edge property objects at the given topology.
 
  @param[in] top - given topology whose edge property objects will be searched for doubled nodes
  @param[in] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                       on the interface (to be doubled)

  @return The function returns the number of edge property objects that contain doubled nodes.
  
  Created by TKo, 01.2024
*/
long num_doubled_nodes_edg_prop(const siftop &top, const ivector &gdn_lst)
{
  long ret = 0;
  sedge *edg = NULL;
  
  if (top.edges){
    for (long i=0; i<top.edges->nedg; i++){
      edg = top.edges->edges + i;      
      if ((gdn_lst[edg->n1] > 0) || (gdn_lst[edg->n2] > 0))
        ++ret;
    }
  }
  return ret;
}



/**
  The function searches for the doubled nodes given in the gdn_lst in the surface property objects at the given topology.
 
  @param[in] top - given topology whose surface property objects will be searched for doubled nodes
  @param[in] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                       on the interface (to be doubled)

  @return The function returns the number of surface property objects that contain doubled nodes.
  
  Created by TKo, 01.2024
*/
long num_doubled_nodes_surf_prop(const siftop &top, const ivector &gdn_lst)
{
  long ret = 0;
  sface *sfc = NULL;
  
  if (top.surfaces){
    for (long i=0; i<top.surfaces->nfcs; i++){
      sfc = top.surfaces->faces + i;
      for (long j=0; j<sfc->nnod; j++){
        if (gdn_lst[sfc->nodes[j]] > 0){
          ++ret;
          break;
        }
      }
    }
  }
  return ret;
}



/**
  The function adds new edge property objects to the given topology. These objects are derived from the
  current one that contains nodes to be doubled.
 
  @param[in] top - given topology where the new edge property objects will be stored,
                   array of the edge property objects must be allocated to hold the newly generated ones
  @param[in] num_dnedgprop - number of new edge property objects that will be added to the topology
  @param[in] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                       on the interface (to be doubled)

  @return The function does not return anything.
  
  Created by TKo, 01.2024
*/
void doubled_nodes_add_edg_prop(siftop &top, long num_dnedgprop, const ivector &gdn_lst)
{
  sedge *edg = NULL;
  
  if (top.edges){
    long k = top.edges->nedg - num_dnedgprop;
    long n = k;
    for (long i=0; i<n; i++){
      edg = top.edges->edges + i;      
      if ((gdn_lst[edg->n1] > 0) || (gdn_lst[edg->n2] > 0)){
        edg->copyto(top.edges->edges[k]);
        if (gdn_lst[edg->n1] > 0)
          top.edges->edges[k].n1 = edg->n1;
        if (gdn_lst[edg->n2] > 0)
          top.edges->edges[k].n2 = edg->n2;
        k++;
      }
    }
  }
}



/**
  The function adds new surface property objects to the given topology. These objects are derived from the
  current one that contains nodes to be doubled.
 
  @param[in] top - given topology where the new surface property objects will be stored,
                   array of the surface property objects must be allocated to hold the newly generated ones
  @param[in] num_dnsurfprop - number of new surface property objects that will be added to the topology
  @param[in] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                       on the interface (to be doubled)

  @return The function does not return anything.
  
  Created by TKo, 01.2024
*/
void doubled_nodes_add_surf_prop(siftop &top, long num_dnsurfprop, const ivector &gdn_lst)
{
  sface *sfc = NULL;
  
  if (top.surfaces){
    long k = top.surfaces->nfcs - num_dnsurfprop;
    long n = k;
    for (long i=0; i<n; i++){
      sfc = top.surfaces->faces + i;
      bool create = false;
      for (long j=0; j<sfc->nnod; j++){
        if (gdn_lst[sfc->nodes[j]] > 0){
          create = true;
          break;
        }
      }
      if (create){
        sface *sfcn = top.surfaces->faces + k;
        sfc->copyto(*sfcn);
        for (long j=0; j<sfc->nnod; j++){
          if (gdn_lst[sfc->nodes[j]] > 0)
            sfcn->nodes[j] = gdn_lst[sfc->nodes[j]];
        }
        k++;
      }
    }
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
  @param[in] entnod - array of nodes on the boundary entity (edge or surface) of the eid-th element which 
                      the new interface element will be connected to.
  @param[in] gdn_lst - global list of node numbers to be doubled for the whole problem, gdn_lst[i] > 0 if the node is 
                       on the interface (to be doubled)
  @param[in] gopt - structure with the options for mesh format used
  @param[in] ifceld - the generation setup of the new interface element

  @return The function does not return anything.

  Created by TKo, 01.2024
*/
void make_ifc_elem(siftop &top, long eid, const ivector &entnod, const ivector &gdn_lst, const genopt &gopt, const ifcel_descript &ifceld, long k)
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
          if (lhs.eid2 < rhs.eid2)
            return true;
          else{
            if ((lhs.eid2 == rhs.eid2) && !(*lhs.ifceld_ptr == *lhs.ifceld_ptr)){
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
  }
  return false;
}
