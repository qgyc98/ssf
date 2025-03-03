#include "siftop.h"
#include "kwdset.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <stdlib.h>
#include <math.h>

enum dir_id {x_plus=0, y_plus, z_plus, x_minus, y_minus, z_minus};
const enumstr dir_idstr[] {{"x+",0}, {"y+",1}, {"z+",2}, {"x-",3}, {"y-",4}, {"z-",5}};
const kwdset dir_id_kwdset(sizeof(dir_idstr)/sizeof(*dir_idstr), dir_idstr);

enum gen_elem_type{et_auto=0, et_brick, et_wedge, et_tetra};
const enumstr gen_elem_typestr[] {{"et_auto",0}, {"et_brick", 1}, {"et_wedge", 2}, {"et_tetra", 3}};
const kwdset gen_elem_type_kwdset(sizeof(gen_elem_typestr)/sizeof(*gen_elem_typestr), gen_elem_typestr);

enum extgentype{egt_no=0, egt_glob_coord, egt_relative_sections};
const enumstr extgentypestr[] {{"egt_no", 0}, {"egt_glob_coord", 1}, {"egt_relative_sections", 2}};
const kwdset extgentype_kwdset(sizeof(extgentypestr)/sizeof(*extgentypestr), extgentypestr);

void gen_brick(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp, long layid, long nlay);

void gen_wedge(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp, long layid, long nlay);

void gen_tetra(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp, long layid, long nlay);


/*
  Generator for extension of 3D mesh segment in the given direction

  Created by Tomas Koudelka, 5.10.2018
*/
int main (int argc,char *argv[])
{
  long i, j, k, l;
  meshform fmt;
  long nsec;
  answertype edg, paral;
  double clev, ctol, cgcoord, auxcoord;
  dir_id dir;
  XFILE *in;
  FILE *out;
  long ret;
  char topname[1001];
  siftop top;
  vector secl, dl, gcoord;  
  ivector md;
  long ngcoord;
  long frontp;
  long nid, nid1, nid2, eid;
  long volp;
  long errc;
  extgentype extgt;
  gen_elem_type gen_et;

  fprintf(stdout, "\n\n\n---        GENERATION OF MESH EXTENSION OF 3D DOMAIN         ---\n\n");
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name output_file_name\n\n", argv[0]);
    return(1);
  }
  in = xfopen(argv[1], "rt");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  out = fopen(argv[2], "wt");
  if (out==NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  
  in->kwdmode = sequent_mode;
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &fmt);
  xfscanf(in, "%k %1000a", "mesh_file_name", topname);
  xfscanf(in, "%k%m", "read_edge_property", &answertype_kwdset, &edg);
  xfscanf(in, "%k%m", "read_paral", &answertype_kwdset, &paral);
  xfscanf(in, "%k %m", "gen_dir", &dir_id_kwdset, &dir);
  xfscanf(in, "%k%le", "coord_level", &clev);
  xfscanf(in, "%k%le", "coord_tolerance", &ctol);
  xfscanf(in, "%k%m", "gen_elem_type", &gen_elem_type_kwdset, &gen_et);
  xfscanf(in, "%k%ld", "front_prop", &frontp);
  xfscanf(in, "%k%ld", "vol_prop", &volp);
  xfscanf(in, "%k%m", "extension_gen_type", &extgentype_kwdset, &extgt);
  switch (extgt){
    case egt_glob_coord:
      xfscanf(in, "%k%ld", "num_coord", &ngcoord);
      reallocv(ngcoord, gcoord);
      errc = readv(in, gcoord);
      if (errc){
        if (errc == 2)
          print_err("cannot read global coordinate", __FILE__, __LINE__, __func__);
        return 3;
      }
      break;
    case egt_relative_sections:
      xfscanf(in, "%k%ld", "num_sec", &nsec);
      reallocv(nsec, dl);
      reallocv(nsec, md);
      reallocv(nsec, secl);
      ngcoord = 0; // total number of 3D sections generated in starting with the 3D domain segment
      for (i=0; i<nsec; i++){
        xfscanf(in, "%k%le %k%ld", "secblock_lenght", secl.a+i, "mesh_dens", md.a+i);
        dl[i] = secl[i]/md[i];
        if (dir > z_plus)  dl[i] = -dl[i];
        ngcoord += md[i];
      }
      reallocv(ngcoord, gcoord);
      k = 0;
      cgcoord = 0.0;
      for (i=0; i<nsec; i++){
        for(j=0; j<md[i]; j++){
          cgcoord += dl[i];
          gcoord[k] = cgcoord;
          k++;
        }
      }
      break;
    default:
      print_err("uknown type (%d) of extension generation", __FILE__, __LINE__, __func__, extgt);
      return 3;
  }
  xfclose(in);

  // open and read file with topology of 3D section
  fprintf(stdout, "Reading of topology file %s with 3D segment ...", topname);
  in = xfopen(topname, "rt");
  if (in == NULL)
    return(3);

  switch (fmt){
    case sifel:
      ret = top.read(in, paral, edg);
      break;
    case t3d:
      ret = top.import_t3d(in, paral);
      break;
    default:
      print_err("unknown mesh format %d is required", __FILE__, __LINE__, __func__, fmt);
      return 3;
  }
  xfclose(in);
  if (ret)
    return 3;

  // find interface nodes where the generation of nodes will start from
  std::vector<long> ifcnod; // array of found interface nodes
  
  long c_id = dir%3; // index of coordinate in the extension direction
  ivector sfcnodid(top.nn);
  long newid = top.nn;
  vector coord(ASTCKVEC(3));
  for(i=0; i< top.nn; i++){
    top.nodes[i].getcoord(coord);
    if (fabs(coord[c_id] - clev) <= ctol){
      ifcnod.push_back(i);
      sfcnodid[i] = newid;
      ++newid;
    }
  }
  long nifcn = ifcnod.size();

  // find elements connected to the interface nodes where the generation of nodes will start from
  std::vector<std::pair<long,long>> ifcelem; // array of pairs containing found element id and corresponding surface id
  ivector surfid_lst;
  long nelemlay = 0, nnsurf;
  for(i=0; i<top.ne; i++){
    reallocv(RSTCKIVEC(top.elements[i].nsurf, surfid_lst));
    long ns = top.elements[i].compare_surf(ifcnod.data(), ifcnod.size(), surfid_lst);
    for (j=0; j<ns; j++){
      ifcelem.push_back(std::make_pair(i, surfid_lst[j]));
      nnsurf = top.elements[i].nnsurf[surfid_lst[j]];
      switch(nnsurf){
        case 3:
        case 6:
          if (gen_et == et_tetra)
            nelemlay += 3;
          else
            nelemlay += 1;
          break;
        case 4:
        case 8:
          nelemlay += 1;
          break;
        default:
          print_err("interface entities has wrong number of nodes, nnsurf must be either 3, 4, 6 or 8", __FILE__, __LINE__,__func__);
          abort();
      }
    }
  }

  long nifce = ifcelem.size();
  fprintf(stdout, "OK\n");
  fprintf(stdout, "\nNumber of nodes in 3D starter segement  : %ld\n", top.nn);
  fprintf(stdout, "Number of elements in 3D starter segment: %ld\n", top.ne);
  fprintf(stdout, "Number of interface nodes found on 3D starter segment: %ld\n", nifcn);
  fprintf(stdout, "Number of interface element enitites found on 3D starter segment: %ld\n", nifce);
  fprintf(stdout, "\nGeneration of 3D domain ...");


  //
  // printing of nodes
  //

  long inod;
  // add volume property numbers at interface nodes
  for(i=0; i<nifcn; i++){
    inod = ifcnod[i];
    top.nodes[inod].add_prop(eregion, volp);
  }
  
  // add surface property numbers generated from edge ones at interface nodes
  for(i=0; i<nifcn; i++){
    inod = ifcnod[i];
    for (j=0; j<top.nodes[inod].nprop; j++){ 
      if (top.nodes[inod].entid[j] == ecurve)
        top.nodes[inod].add_prop(esurface, top.nodes[inod].prop[j]);
    }
  }
  
  // print nodes of intial (first) 3D segment
  fprintf(out, "%ld\n", top.nn+ngcoord*nifcn);
  top.shift_print_nodes(out, 0);
  
  double snode::* coordptr;  // pointer to the data member of class snode
  switch(dir)
  {
    case x_plus:
    case x_minus:
      coordptr = &snode::x; // coordptr referes to the x coordinate in class snode
      break;
    case y_plus:
    case y_minus:
      coordptr = &snode::y; // coordptr referes to the y coordinate in class snode
      break;
    case z_plus:
    case z_minus:
      coordptr = &snode::z; // coordptr referes to the z coordinate in class snode
      break;
    default:
      print_err("invalid direction indicator %d, it must be in the range <0;5>", __FILE__, __LINE__, __func__, int(dir));
      return 4;
  }
  nid = top.nn;
  std::vector<entityp> nodent;
  std::vector<long> nodprop;
  long np;
  bool addvolp = true;
  entityp ent;
  
  for (i=0; i<ngcoord; i++){
    for (k=0; k<nifcn; k++){
      inod = ifcnod[k];
      auxcoord = top.nodes[inod].*coordptr; // save actual coordinate value of the interface node either x,y or z according to coordptr
      top.nodes[inod].*coordptr = gcoord[i]; // change the coordinate value value
      fprintf (out,"%ld %15.10le %15.10le %15.10le",nid+1, top.nodes[inod].x, top.nodes[inod].y, top.nodes[inod].z);
      top.nodes[inod].*coordptr = auxcoord; // restore the coordinate value
      // create arrays with nodal properties
      nodent.clear();
      nodprop.clear();
      np = top.nodes[inod].nprop;
      nodent.reserve(np);
      nodprop.reserve(np);
      nodent.insert(nodent.begin(), top.nodes[inod].entid, top.nodes[inod].entid+np);
      nodprop.insert(nodprop.begin(), top.nodes[inod].prop, top.nodes[inod].prop+np);
      addvolp = true;
      for(l=0; l<np; l++){
        ent = nodent[l];        
        if (ent == ecurve)   nodent[l] = esurface; // edge properties are changed to surface ones
        if (ent == esurface) nodprop[l] = 0; // surface properties are dropped
        if (ent == eregion){ // volume properties are dropped or changed to volp
          if (addvolp){ // new nodes should have just one volume property equal to volp
            nodprop[l] = volp;
            addvolp = false;
          }
          else
            nodprop[l] = 0;
        }
      }
      if (addvolp){ // no volume properties were found in the interface nodes
        nodent.push_back(eregion);
        nodprop.push_back(volp);
      }
      if (i == ngcoord - 1){
        nodent.push_back(esurface);
        nodprop.push_back(frontp);
      }
        
      fprintf (out," %ld", long(nodent.size()));      
      for (l=0; l<long(nodent.size()); l++)
        fprintf(out, " %d %ld", int(nodent[l]), nodprop[l]);
      fprintf(out, "\n");
      nid++;
    }
  }  
  fprintf(out, "%ld\n", top.ne+nelemlay*ngcoord);

  //
  //  printing of elements
  //
  eid = top.ne;
  nid1 = 0;
  nid2 = 0;
  // print elements of starter 3D segment
  top.shift_print_elements(out, 0, 0);
  ivector sfcnods;
  long eidifc, sfid, nnsfc;
  selement *elem;
  gtypel et;
  for (i=0; i<ngcoord; i++){   
    for (k=0; k<nifce; k++){
      eidifc = ifcelem[k].first;
      sfid = ifcelem[k].second;
      nnsfc = top.elements[eidifc].nnsurf[sfid];
      reallocv(RSTCKIVEC(nnsfc, sfcnods));
      elem = top.elements+eidifc;
      et = top.elements[eidifc].type;
      for (l=0; l<nnsfc; l++)
        sfcnods[l] = elem->nodes[elem->surfnod[sfid][l]];
      switch(gen_et){
        case et_auto:
          if (nnsfc == 4)
            gen_brick(out, top, sfcnods, sfcnodid, sfid, eid, et, eidifc, nid1, nid2, volp, edg, frontp, i, ngcoord);
          else
            gen_wedge(out, top, sfcnods, sfcnodid, sfid, eid, et, eidifc, nid1, nid2, volp, edg, frontp, i, ngcoord);
          break;
        case et_brick:
          gen_brick(out, top, sfcnods, sfcnodid, sfid, eid, et, eidifc, nid1, nid2, volp, edg, frontp, i, ngcoord);
          break;
        case et_tetra:
          gen_tetra(out, top, sfcnods, sfcnodid, sfid, eid, et, eidifc, nid1, nid2, volp, edg, frontp, i, ngcoord);
          break;
        case et_wedge:
          gen_wedge(out, top, sfcnods, sfcnodid, sfid, eid, et, eidifc, nid1, nid2, volp, edg, frontp, i, ngcoord);
          break;
        default:
          print_err("unknown type of generated elements (%d) is required", __FILE__, __LINE__, __func__, int(gen_et));
          abort();
      }
    }
    nid1 = nid2;  // increase start index of nodes on the front base surface 
    nid2 += nifcn;  // increase start index of nodes on the rear base surface 
  }
  fclose(out);

  fprintf(stdout, "OK\n\n");
  fprintf(stdout, "The total number of nodes in 3D domain   : %ld\n", top.nn + nifcn*ngcoord);
  fprintf(stdout, "The total number of elements in 3D domain: %ld\n", top.ne + nifce*ngcoord);
  fprintf(stdout, "3D topology has been written in file %s\n", argv[2]);
  
  return 0;
}



/**
  The function generates one new element of the linear brick type over the given surface of the 
  interface element of the base 3D mesh segment.

  @param[in, out] out - pointer to the opened text file, where the given element is written to
  @param[in] sfcnods - array of nodes defining the element surface of the base 3D mesh segment interface
                       which the new element should be generated over
  @param[in] sfcid - index of surface defined by sfcnods on the interface element of the base 3D mesh segment 
  @param[in, out] eid - element id of the new generated element
  @param[in] et - element type of the interface element of the base 3D mesh segment
  @param[in] ifceid - element id of the interface element of the base 3D mesh segment connected with the
                      surface nodes sfcnods
  @param[in] sfcnodid - nodal indices of the first layer of new generated nodes
  @param[in] nid1 - index offset of the bottom nodes of the new generated element, the offset is related to sfcnodid
  @param[in] nid2 - index offset of the top nodes of the new generated element, the offset is related to sfcnodid
  @param[in] volp - volume property ov the new generated element
  @param[in] edg - flag for generation of edge/surface property numbers (edg == yes -> property numbers will be generated)
  @param[in] frontp - surface property of the final front surface of the new generated mesh segment
  @param[in] layid - actual id of the element layer generated
  @param[in] nlay - the total number of element layers generated

  Created by Tomas Koudelka, 03.2024, tomas.koudelka@fsv.cvut.cz
*/
void gen_brick(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp,
               long layid, long nlay)
{
  long nnsfc = sfcnods.n;
  long k, l, m, aux;
  ivector sfcnodsord(ASTCKIVEC(nnsfc));
  ivector sfcedgord(ASTCKIVEC(nnsfc));
  static selement isolin3d(isolinear3d, 0);

  if ((et != isolinear3d) && (et != tetrahedronlinear)){
    print_err("unknown type (%d) of interface element %ld is required", __FILE__, __LINE__, __func__, int(et), ifceid+1);
    abort();
  }

  if ((nnsfc < 3) || (nnsfc > 4)){
    print_err("invalid number of surface nodes (%ld) is required for a new brick element %ld",
              __FILE__, __LINE__, __func__, nnsfc, eid+1);
    abort();
  }
  // create anti-clockwise ordering of surface nodes with respect to its outter normal direction
  if (et == isolinear3d){
    if ((sfcid < 2) || (sfcid == 4)){
      copyv(sfcnods, sfcnodsord);
      for(l=0; l<nnsfc; l++)
        sfcedgord[l] = l;
    }
    else{
      for(l=nnsfc-1, m=0; l>=0; --l, ++m)
        sfcnodsord[m] = sfcnods[l];
      sfcedgord[0] = 2;   sfcedgord[1] = 1;
      sfcedgord[2] = 0;   sfcedgord[3] = 3;
    }
  }
  if (et == tetrahedronlinear){
    for(l=nnsfc-1, m=0; l>=0; --l, ++m)
      sfcnodsord[m] = sfcnods[l];
    sfcedgord[0] = 1;   sfcedgord[1] = 0;
    sfcedgord[2] = 2;
  }
  
  if (nnsfc == 4){
    fprintf(out, "%ld %d ", eid+1, int(isolinear3d));
    for (l=0; l<nnsfc; l++) 
      fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid2);
    for (l=0; l<nnsfc; l++){
      if (layid == 0)
        fprintf(out, "%ld ", sfcnodsord[l]+1);
      else
        fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid1);
    }
    fprintf(out, "  %ld", volp);
    if (edg == yes){
      // edge property id of brick element are not preserved, 
      // they will be changed to surface property id 
      fprintf(out, "  ");
      for (l=0; l<isolin3d.ned; l++)
        fprintf(out, " %d", 0);
      fprintf(out, "  ");
      
      // property id of side surfaces on the brick element
      k = nnsfc-1;
      for (l=0; l<nnsfc; l++, k++){
        // indeces of brick side surfaces are shifted by 3 with respect 
        // to edge indeces on quadrilateral
        m = top.elements[ifceid].surfedg[sfcid][sfcedgord[k%nnsfc]];
        fprintf(out, " %ld", top.elements[ifceid].propedg[m]);
      }
      
      // property id of top surface of the last layer of new generated elements
      aux = 0;
      if (layid == nlay-1)
        aux = frontp;
      fprintf(out, " %ld 0", aux);
    }
    fprintf(out, "\n");
  }

  if (nnsfc == 3){
    fprintf(out, "%ld %d ", eid+1, int(isolinear3d));
    
    for (l=0; l<nnsfc; l++) 
      fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid2);
    fprintf(out, "%ld ", sfcnodid[sfcnodsord[nnsfc-1]]+1+nid2);
    for (l=0; l<nnsfc; l++){
      if (layid == 0)
        fprintf(out, "%ld ", sfcnodsord[l]+1);
      else
        fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid1);
    }
    if (layid == 0)
      fprintf(out, "%ld ", sfcnodsord[nnsfc-1]+1);
    else
      fprintf(out, "%ld ", sfcnodid[sfcnodsord[nnsfc-1]]+1+nid1);
    fprintf(out, "%ld", volp);
    if (edg == yes){
      // edge property id of brick element are not preserved, 
      // they will be changed to surface property id 
      fprintf(out, "  ");
      for (l=0; l<isolin3d.ned; l++)
        fprintf(out, " %d", 0);
      fprintf(out, "  ");
      
      fprintf(out, " 0"); // surface with zero area due to degenerated brick element
      // property id of side surfaces on the brick element
      k = nnsfc-1;
      for (l=0; l<nnsfc; l++, k++){ // nnsfc corresponds to number of edges of the given surface
        // index of the corresponding adjacent edge to the given surface of interface element
        m = top.elements[ifceid].surfedg[sfcid][sfcedgord[k%nnsfc]];
        fprintf(out, " %ld", top.elements[ifceid].propedg[m]);
      }
      
      // property id of top surface of the last layer of new generated elements
      aux = 0;
      if (layid == nlay-1)
        aux = frontp;
      
      fprintf(out, "0 %ld 0", aux);
    }
    fprintf(out, "\n");
  }
  ++eid;
}



/**
  The function generates one new element of the linear wedge type over the given surface of the 
  interface element of the base 3D mesh segment.

  @param[in, out] out - pointer to the opened text file, where the given element is written to
  @param[in] sfcnods - array of nodes defining the element surface of the base 3D mesh segment interface
                       which the new element should be generated over
  @param[in] sfcid - index of surface defined by sfcnods on the interface element of the base 3D mesh segment 
  @param[in, out] eid - element id of the new generated element
  @param[in] et - element type of the interface element of the base 3D mesh segment
  @param[in] ifceid - element id of the interface element of the base 3D mesh segment connected with the
                      surface nodes sfcnods
  @param[in] sfcnodid - nodal indices of the first layer of new generated nodes
  @param[in] nid1 - index offset of the bottom nodes of the new generated element, the offset is related to sfcnodid
  @param[in] nid2 - index offset of the top nodes of the new generated element, the offset is related to sfcnodid
  @param[in] volp - volume property ov the new generated element
  @param[in] edg - flag for generation of edge/surface property numbers (edg == yes -> property numbers will be generated)
  @param[in] frontp - surface property of the final front surface of the new generated mesh segment
  @param[in] layid - actual id of the element layer generated
  @param[in] nlay - the total number of element layers generated

  Created by Tomas Koudelka, 03.2024, tomas.koudelka@fsv.cvut.cz
*/
void gen_wedge(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp,
               long layid, long nlay)
{
  long nnsfc = sfcnods.n;
  long k, l, m, aux;
  ivector sfcnodsord(ASTCKIVEC(nnsfc));
  ivector sfcedgord(ASTCKIVEC(nnsfc));
  static selement wedgelin(wedgelinear, 0);

  if (et != tetrahedronlinear){
    print_err("unknown type (%d) of interface element %ld is required", __FILE__, __LINE__, __func__, int(et), ifceid+1);
    abort();
  }

  if ((nnsfc < 3) || (nnsfc > 4)){
    print_err("invalid number of surface nodes (%ld) is required for a new brick element %ld",
              __FILE__, __LINE__, __func__, nnsfc, eid+1);
    abort();
  }

  // create anti-clockwise ordering of surface nodes with respect to its outter normal direction
  for(l=nnsfc-1, m=0; l>=0; --l, ++m)
    sfcnodsord[m] = sfcnods[l];
  sfcedgord[0] = 1;   sfcedgord[1] = 0;
  sfcedgord[2] = 2;
  
  fprintf(out, "%ld %d ", eid+1, int(wedgelinear));
    
  for (l=0; l<nnsfc; l++) 
    fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid2);
  for (l=0; l<nnsfc; l++){
    if (layid == 0)
      fprintf(out, "%ld ", sfcnodsord[l]+1);
    else
      fprintf(out, "%ld ", sfcnodid[sfcnodsord[l]]+1+nid1);
  }
  fprintf(out, "  %ld", volp);
  if (edg == yes){
    // edge property id of brick element are not preserved, 
    // they will be changed to surface property id 
    fprintf(out, "  ");
    for (l=0; l<wedgelin.ned; l++)
      fprintf(out, " %d", 0);
    fprintf(out, "  ");
    
    // property id of side surfaces on the brick element
    k = nnsfc-1;
    for (l=0; l<nnsfc; l++,k++){ // nnsfc corresponds to number of edges of the given surface
      // index of the corresponding adjacent edge to the given surface of interface element
      m = top.elements[ifceid].surfedg[sfcid][sfcedgord[k%nnsfc]];
      fprintf(out, " %ld", top.elements[ifceid].propedg[m]);
    }
    
    // property id of top surface of the last layer of new generated elements
    aux = 0;
    if (layid == nlay-1)
      aux = frontp;
    
    fprintf(out, " %ld 0", aux);
  }
  fprintf(out, "\n");
  ++eid;
}



/**
  The function generates one new element of the linear brick type over the given surface of the 
  interface element of the base 3D mesh segment.

  @param[in, out] out - pointer to the opened text file, where the given element is written to
  @param[in] sfcnods - array of nodes defining the element surface of the base 3D mesh segment interface
                       which the new element should be generated over
  @param[in,out] eid - element id of the new generated element
  @param[in] ifceid - element id of the interface element of the base 3D mesh segment connected with the
                      surface nodes sfcnods
  @param[in] sfcnodid - nodal indices of the first layer of new generated nodes
  @param[in] nid1 - index offset of the bottom nodes of the new generated element, the offset is related to sfcnodid
  @param[in] nid2 - index offset of the top nodes of the new generated element, the offset is related to sfcnodid
  @param[in] volp - volume property ov the new generated element
  @param[in] edg - flag for generation of edge/surface property numbers (edg == yes -> property numbers will be generated)
  @param[in] frontp - surface property of the final front surface of the new generated mesh segment
  @param[in] layid - actual id of the element layer generated
  @param[in] nlay - the total number of element layers generated
*/
void gen_tetra(FILE *out, const siftop &top, const ivector &sfcnods, const ivector &sfcnodid, long sfcid, long &eid,
               gtypel et, long ifceid, long nid1, long nid2, long volp, answertype edg, long frontp, long layid, long nlay)
{
  long nnsfc = sfcnods.n;
  long l, m;
  ivector sfcnodsord(ASTCKIVEC(nnsfc));
  ivector sfcedgord(ASTCKIVEC(nnsfc));
  static selement tetrahedlin(tetrahedronlinear, 0);
  
  if (nnsfc != 3){
    print_err("invalid number of surface nodes (%ld) is required for a tetrahedron element",
              __FILE__, __LINE__, __func__, nnsfc);
    abort();
  }
  if (et == tetrahedronlinear){
    for(l=nnsfc-1, m=0; l>=0; --l, ++m)
      sfcnodsord[m] = sfcnods[l];
    sfcedgord[0] = 1;   sfcedgord[1] = 0;
    sfcedgord[2] = 2;
  }

  fprintf(out, "%ld %d ", eid+1, int(tetrahedronlinear));
  if (layid == 0){
    fprintf(out, "%ld %ld %ld %ld ", sfcnodsord[0]+1, sfcnodsord[1]+1, sfcnodsord[2]+1, sfcnodid[sfcnodsord[0]]+1+nid2);
  }
  else
    fprintf(out, "%ld %ld %ld %ld ", sfcnodid[sfcnodsord[0]]+nid1+1, sfcnodid[sfcnodsord[1]]+nid1+1,
            sfcnodid[sfcnodsord[2]]+nid1+1, sfcnodid[sfcnodsord[0]]+nid2+1);
  fprintf(out, "%ld", volp);
  if (edg == yes){
    // edge property id of tetrahedron element are not preserved, 
    // they will be changed to surface property id 
    fprintf(out, "  ");
    for (l=0; l<tetrahedlin.ned; l++)
      fprintf(out, " %d", 0);
    fprintf(out, "  ");    
    // surface property numbers of the first tetrahedron
    fprintf(out, " 0 %ld %ld 0", top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[2]]],
            top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[0]]]);
  }
  fprintf(out, "\n");
  ++eid;
    
  fprintf(out, "%ld %d ", eid+1, int(tetrahedronlinear));
  if (layid == 0){
    fprintf(out, "%ld %ld %ld %ld ", sfcnodsord[2]+1, sfcnodsord[1]+1, sfcnodid[sfcnodsord[1]]+nid2+1, sfcnodid[sfcnodsord[0]]+nid2+1);
  }
  else
    fprintf(out, "%ld %ld %ld %ld ", sfcnodid[sfcnodsord[2]]+nid1+1, sfcnodid[sfcnodsord[1]]+nid1+1,
            sfcnodid[sfcnodsord[1]]+nid2+1, sfcnodid[sfcnodsord[0]]+nid2+1);
  fprintf(out, "%ld", volp);
  if (edg == yes){
    // edge property id of tetrahedron element are not preserved, 
    // they will be changed to surface property id 
    fprintf(out, "  ");
    for (l=0; l<tetrahedlin.ned; l++)
      fprintf(out, " %d", 0);
    fprintf(out, "  ");
    // surface property numbers of the second tetrahedron
    fprintf(out, " %ld 0 0 %ld", top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[0]]],
            top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[1]]]);
  }
  fprintf(out, "\n");
  ++eid;

  fprintf(out, "%ld %d ", eid+1, int(tetrahedronlinear));
  if (layid == 0){
    fprintf(out, "%ld %ld %ld %ld ", sfcnodid[sfcnodsord[0]]+nid2+1, sfcnodid[sfcnodsord[1]]+nid2+1,
            sfcnodid[sfcnodsord[2]]+nid2+1, sfcnodsord[2]+1);
  }
  else
    fprintf(out, "%ld %ld %ld %ld ", sfcnodid[sfcnodsord[0]]+nid2+1, sfcnodid[sfcnodsord[1]]+nid2+1,
            sfcnodid[sfcnodsord[2]]+nid2+1, sfcnodid[sfcnodsord[2]]+nid1+1);
  fprintf(out, "%ld", volp);
  if (edg == yes){
    // edge property id of tetrahedron element are not preserved, 
    // they will be changed to surface property id 
    fprintf(out, "  ");
    for (l=0; l<tetrahedlin.ned; l++)
      fprintf(out, " %d", 0);
    fprintf(out, "  ");
    // surface property numbers of the second tetrahedron
    if (layid < nlay-1)
      fprintf(out, " 0 %ld %ld 0", top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[1]]],
              top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[1]]]);
    else
      fprintf(out, " 0 %ld %ld %ld", top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[1]]],
              top.elements[ifceid].propedg[top.elements[ifceid].surfedg[sfcid][sfcedgord[1]]], frontp);
  }
  fprintf(out, "\n");
  ++eid;
}

