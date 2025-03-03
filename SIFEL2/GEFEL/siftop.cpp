#include "siftop.h"
#include "galias.h"
#include "mathem.h"
#include "basefun.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <set>
#include <algorithm>

long selement::nned_isolin1d[1]   = {2};
long selement::nned_isoquad1d[1]  = {3};
long selement::nned_trlin2d[3]    = {2, 2, 2};
long selement::nned_trquad2d[3]   = {3, 3, 3};
long selement::nned_isolin2d[4]   = {2, 2, 2, 2};
long selement::nned_isoquad2d[4]  = {3, 3, 3, 3};
long selement::nned_isocubic2d[4]   = {4, 4, 4, 4};
long selement::nned_tetlin3d[6]   = {2, 2, 2, 2, 2, 2};
long selement::nned_tetquad3d[6]  = {3, 3, 3, 3, 3, 3};
long selement::nned_pyramlin[8]   = {2, 2, 2, 2, 2, 2, 2, 2};
long selement::nned_pyramquad[8]  = {3, 3, 3, 3, 3, 3, 3, 3};
long selement::nned_wedgelin[9]   = {2, 2, 2, 2, 2, 2, 2, 2, 2};
long selement::nned_wedgequad[9]  = {3, 3, 3, 3, 3, 3, 3, 3, 3};
long selement::nned_isolin3d[12]  = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
long selement::nned_isoquad3d[12] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

long selement::edgenod_isolin1d[1][selement::max_nned]   = {{0, 1}};
long selement::edgenod_isoquad1d[1][selement::max_nned]  = {{0, 1, 2}};
long selement::edgenod_trlin2d[3][selement::max_nned]    = {{0, 1}, {1, 2}, {2, 0}};
long selement::edgenod_trquad2d[3][selement::max_nned]   = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
long selement::edgenod_isolin2d[4][selement::max_nned]   = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
long selement::edgenod_isoquad2d[4][selement::max_nned]  = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
long selement::edgenod_isocubic2d[4][selement::max_nned] = {{0, 1, 4, 5}, {1, 2, 6, 7}, {2, 3, 8, 9}, {3, 0, 10, 11}};
long selement::edgenod_tetlin3d[6][selement::max_nned]   = {{0, 1}, {1, 2}, {2, 0}, {1, 3}, {3, 2}, {3, 0}};
long selement::edgenod_tetquad3d[6][selement::max_nned]  = {{0, 1, 4}, {1, 2, 5}, {2, 0, 6}, {1, 3, 8}, {3, 2, 9}, {3, 0, 7}};
long selement::edgenod_pyramlin[8][selement::max_nned]   = {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}}; // not yet defined 
long selement::edgenod_pyramquad[8][selement::max_nned]  = {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}}; // not yet defined 
long selement::edgenod_wedgelin[9][selement::max_nned]   = {{0, 1}, {1, 2}, {2, 0},
                                                            {0, 3}, {1, 4}, {2, 5},
                                                            {3, 4}, {4, 5}, {5, 3}};
long selement::edgenod_wedgequad[9][selement::max_nned]  = {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}}; // not yet defined 
long selement::edgenod_isolin3d[12][selement::max_nned]  = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                            {0, 4}, {1, 5}, {2, 6}, {3, 7},
                                                            {4, 5}, {5, 6}, {6, 7}, {7, 4}};
long selement::edgenod_isoquad3d[12][selement::max_nned] = {{0, 1, 8}, {1, 2, 9}, {2, 3, 10}, {3, 0, 11},
                                                            {0, 4, 12}, {1, 5, 13}, {2, 6, 14}, {3, 7, 15},
                                                            {4, 5, 16}, {5, 6, 17}, {6, 7, 18}, {7, 4, 18}};
long selement::nnsurf_isolin1d[1]   = {0};
long selement::nnsurf_isoquad1d[1]  = {0};
long selement::nnsurf_trlin2d[1]    = {3};
long selement::nnsurf_trquad2d[1]   = {6};
long selement::nnsurf_isolin2d[1]   = {4};
long selement::nnsurf_isoquad2d[1]  = {8};
long selement::nnsurf_isocubic2d[1] = {12};
long selement::nnsurf_tetlin3d[4]  = {3, 3, 3, 3};
long selement::nnsurf_tetquad3d[4] = {6, 6, 6, 6};
long selement::nnsurf_pyramlin[5]  = {3, 3, 3, 3, 4};
long selement::nnsurf_pyramquad[5] = {6, 6, 6, 6, 8};
long selement::nnsurf_wedgelin[5]  = {4, 4, 4, 3, 3};
long selement::nnsurf_wedgequad[5] = {8, 8, 8, 6, 6};
long selement::nnsurf_isolin3d[6]  = {4, 4, 4, 4, 4, 4};
long selement::nnsurf_isoquad3d[6] = {8, 8, 8, 8, 8, 8};

long selement::surfnod_isolin1d[1][selement::max_nnsurf]  = {{0, 1}};
long selement::surfnod_isoquad1d[1][selement::max_nnsurf] = {{0, 1, 2}};
long selement::surfnod_trlin2d[1][selement::max_nnsurf]   = {{0, 1, 2}};
long selement::surfnod_trquad2d[1][selement::max_nnsurf]  = {{0, 1, 2, 3, 4, 5}};
long selement::surfnod_isolin2d[1][selement::max_nnsurf]  = {{0, 1, 2, 3}};
long selement::surfnod_isoquad2d[1][selement::max_nnsurf] = {{0, 1, 2, 3, 4, 5, 6, 7}};
long selement::surfnod_isocubic2d[1][selement::max_nnsurf] = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
long selement::surfnod_tetlin3d[4][selement::max_nnsurf]  = {{2, 1, 3},  {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
long selement::surfnod_tetquad3d[4][selement::max_nnsurf] = {{2, 1, 3, 5, 8, 9}, {0, 2, 3, 6, 9, 7}, 
                                                             {1, 0, 3, 4, 7, 8}, {0, 1, 2, 4, 5, 6}};
long selement::surfnod_pyramlin[5][selement::max_nnsurf]  = {{0}, {0}, {0}, {0}, {0}};
long selement::surfnod_pyramquad[5][selement::max_nnsurf] = {{0}, {0}, {0}, {0}, {0}};
long selement::surfnod_wedgelin[5][selement::max_nnsurf]  = {{2, 1, 4, 5}, {0, 2, 5, 3}, {1, 0, 3, 4}, {0, 1, 2}, {3, 5, 4}};
long selement::surfnod_wedgequad[5][selement::max_nnsurf] = {{0}, {0}, {0}, {0}, {0}};
long selement::surfnod_isolin3d[6][selement::max_nnsurf]  = {{0, 3, 7, 4}, {1, 0, 4, 5}, {2, 1, 5, 6}, 
                                                             {3, 2, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7}};
long selement::surfnod_isoquad3d[6][selement::max_nnsurf] = {{0, 3, 7, 4, 11, 15, 19, 12}, {1, 0, 4, 5, 8, 12, 16, 13}, 
                                                             {2, 1, 5, 6, 9, 13, 17, 14}, {3, 2, 6, 7, 10, 14, 18, 15}, 
                                                             {0, 1, 2, 3, 8, 9, 10, 11},  {4, 5, 6, 7, 16, 17, 18, 19}};



long selement::surfedg_isolin1d[1][selement::max_nsurfedg] = {{0}};
long selement::surfedg_isoquad1d[1][selement::max_nsurfedg] = {{0}};
long selement::surfedg_trlin2d[1][selement::max_nsurfedg] = {{0}};
long selement::surfedg_trquad2d[1][selement::max_nsurfedg] = {{0, 1, 2}};
long selement::surfedg_isolin2d[1][selement::max_nsurfedg] = {{0, 1, 2}};
long selement::surfedg_isoquad2d[1][selement::max_nsurfedg] = {{0, 1, 2, 3}};
long selement::surfedg_isocubic2d[1][selement::max_nsurfedg] = {{0, 1, 2, 3}};
long selement::surfedg_tetlin3d[4][selement::max_nsurfedg] = {{1, 3, 4}, {2, 4, 5}, {0, 5, 3}, {0, 1, 2}};
long selement::surfedg_tetquad3d[4][selement::max_nsurfedg] = {{1, 3, 4}, {2, 4, 5}, {0, 5, 3}, {0, 1, 2}};
long selement::surfedg_pyramlin[5][selement::max_nsurfedg] =  {{0}, {0}, {0}, {0}, {0}};
long selement::surfedg_pyramquad[5][selement::max_nsurfedg] = {{0}, {0}, {0}, {0}, {0}};
long selement::surfedg_wedgelin[5][selement::max_nsurfedg] = {{2, 4, 7, 5}, {2, 5, 8, 3}, {0, 3, 6, 4}, {0, 1, 2}, {8, 7, 6}};
long selement::surfedg_wedgequad[5][selement::max_nsurfedg] = {{0}, {0}, {0}, {0}, {0}};
long selement::surfedg_isolin3d[6][selement::max_nsurfedg] = {{3, 7, 11, 4}, {0, 4, 8, 5}, {1, 5, 9, 6},
                                                              {2, 6, 10, 7}, {0, 1, 2, 3}, {8, 9, 10, 11}};
long selement::surfedg_isoquad3d[6][selement::max_nsurfedg] = {{3, 7, 11, 4}, {0, 4, 8, 5}, {1, 5, 9, 6},
                                                               {2, 6, 10, 7}, {0, 1, 2, 3}, {8, 9, 10, 11}};



gtypel operator ++(gtypel &a, int)
{
  gtypel old;
  switch (a)
  {
    case isolinear1d:
      old = a;
      a = isoquadratic1d;
      return old;
    case isoquadratic1d:
      old = a;
      a = trianglelinear;
      return old;
    case trianglelinear:
      old = a;
      a = trianglequadratic;
      return old;
    case trianglequadratic:
      old = a;
      a = isolinear2d;
      return old;
    case isolinear2d:
      old = a;
      a = isoquadratic2d;
      return old;
    case isoquadratic2d:
      old = a;
      a = tetrahedronlinear;
      return old;
    case tetrahedronlinear:
      old = a;
      a = tetrahedronquadratic;
      return old;
    case tetrahedronquadratic:
      old = a;
      a = pyramidelinear;
      return old;
    case pyramidelinear:
      old = a;
      a = pyramidequadratic;
      return old;
    case pyramidequadratic:
      old = a;
      a = wedgelinear;
      return old;
    case wedgelinear:
      old = a;
      a = wedgequadratic;
      return old;
    case wedgequadratic:
      old = a;
      a = isolinear3d;
      return old;
    case isolinear3d:
      old = a;
      a = isoquadratic3d;
      return old;
    case isoquadratic3d:
      old = a;
      a = isocubic2d;
      return old;
    case isocubic2d:
      old = a;
      a = all_element_types;
      return old;
    default:
      print_err("increasing of invalid gtypel value is required", __FILE__, __LINE__, __func__);
  }
  return all_element_types;
}



/**
  Constructor setups all data members to zero values
  
  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
snode::snode()
{
  x = y = z = 0.0;
  nprop = 0;
  entid = NULL;
  prop  = NULL;
}



/**
  Destructor deletes dynamically allocated data members
  
  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
snode::~snode()
{
  delete [] entid;
  delete [] prop;
}


/**
  This method allocates array properties which size is given by num_prop

  @param num_prop - number of nodes of given element i.e. size of node array

  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
void snode::alloc(long num_prop)
{
  nprop = num_prop;
  if (entid)
    delete [] entid;
  entid = new entityp[nprop];
  if (prop)
    delete [] prop;
  prop  = new long[nprop];
}



/**
  This method adds one new property to the array of properties. The entity type of added property 
  is given by entid.

  @param ent - entity type of added property
  @param propid - property number

  Created by Tomas Koudelka 10.2012, tomas.koudelka@fsv.cvut.cz
*/
void snode::add_prop(entityp ent, long propid)
{
  long n = nprop;
  long i;
   
  for(i=0; i<n; i++)
  {
    if ((entid[i] == ent) && (prop[i] == propid)) // added property has been already involved
      return;
  }

  entityp *auxent  = new entityp[n+1];
  long    *auxprop = new long[n+1];
  for(i=0; i<n; i++){
    auxent[i] = entid[i];
    auxprop[i] = prop[i];    
  }
  auxent[n] = ent;
  auxprop[n] = propid;
  delete [] entid;
  delete [] prop;
  entid = auxent;
  prop = auxprop;
  ++nprop;
}



/**
  This method searches array of properties for property aprop with given entity type ent.
  If the pointers edg and sfc are not NULL then the additional edges edg or surfaces sfc are
  also searched for the given property.

  @param prop - searched property number 
  @param ent  - general entity type whose property is searched
  @param nid  - node id (it is used in connection with the additional edg and sfc data)
  @param edg  - pointer to structure with data about additional edges
  @param sfc  - pointer to structure with data about additional surfaces

  @retval 0 - searched property was not found on given entity type
  @retval 1 - searched property was found

  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
long snode::searchprop(long aprop, gentity ent, long nid, sedges *edg, sfaces *sfc)
{
  long i, j;
  for(i=0; i<nprop; i++)
  {
    switch(ent)
    {
      case gvertex:
      case gcurve:
      case gregion: 
        if ((int(ent) == int(entid[i])) && (prop[i] == aprop))   return 1;
        break;
      case gsurface:
        if (((entid[i] == esurface) || (entid[i] == epatch) || 
             (entid[i] == eshell)) && (prop[i] == aprop))    return 1;
	break;
      default:
        print_err("Unknown type of entity is required", __FILE__, __LINE__, __func__);
        abort();
    }
  }
 
  if ((ent == gcurve) && (edg))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edg->nprop; i++)
    {
      if (edg->proplist[i] == aprop)
      {
        for(j=0; j<edg->prop_list_nod_num[i]; j++)
        {
          if (edg->prop_list_nod[i][j] == nid)   return 1;
          if (edg->prop_list_nod[i][j] > nid)    break;  // node numbers are sorted ascendently in the prop_list_nod
        }
      }
    }
  }

  if ((ent == gsurface) && (sfc))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<sfc->nprop; i++)
    {
      if (sfc->proplist[i] == aprop)
      {
        for(j=0; j<sfc->prop_list_nod_num[i]; j++)
        {
          if (sfc->prop_list_nod[i][j] == nid)   return 1;
          if (sfc->prop_list_nod[i][j] > nid)    break;  // node numbers are sorted ascendently in the prop_list_nod
        }
      }
    }
  }
  return 0;
}



/**
  The function searches for the given entity type in property entity type array of the node and returns the number of   
  detected occurences.

  @param e - required entity type

  @return The function returns the number of properties of the given entity type

  Created by Tomas Koudelka, 12.2.2019
*/
long snode::searchpropent(entityp e)
{
  long i, ret;
  ret = 0L;
  for(i=0L; i<nprop; i++){
    if (entid[i] == e)
      ret++;
  }
  return ret;
}



/**
  The function returns coordinates of the given node as a %vector.

  @param[out] p - output parameter for coordinate storage

  @return The function returns nodal coordinates in the argument p.

  Created by Tomas Koudelka, 27.8.2018
*/
void snode::getcoord(vector &p)
{
  p[0] = x;
  p[1] = y;
  p[2] = z;
}



/**
  The function makes a copy of the given node.

  @param c[out] - object of node which will become copy of the given node.
                  previous arrays of the node c are reallocated. 

  @return The function does not return anything but creates a copy in the argument c.
*/
void snode::copyto(snode &dest)
{
  dest.x = x;
  dest.y = y;
  dest.z = z;
  dest.alloc(nprop);
  memcpy(dest.entid, entid, sizeof(*entid)*nprop);
  memcpy(dest.prop, prop, sizeof(*prop)*nprop);
}



/**
  Constructor setups all data members to zero values
  
  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
selement::selement()
{
  type = gtypel(0);
  nne = 0L;

  ned =  tnned = 0L;
  nned     = NULL;
  propedg  = NULL;
  edgenod  = NULL;

  nsurf = tnnsurf = 0;
  nnsurf   = NULL;
  surfnod  = NULL;
  propsurf = NULL;  

  nodes    = NULL;
  prop     = 0L;
}



/**
The constructor sets data members for the given element type et and allocates property arrays on edges and surfaces according to the alloc_edg flag.

  @param[in] et - the given element type
  @param[in] alloc_bp - allocation flag for allocation boundary property arrays (edge/surface), 
                        1 = allocation is being performed, 0 = no allocation is being performed

  Created by Tomas Koudelka 03.2024, tomas.koudelka@fsv.cvut.cz
*/
selement::selement(gtypel et, long alloc_bp)
{
  type = et;
  alloc(alloc_bp);
}



/**
  Destructor deletes dynamically allocated data members
  
  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
selement::~selement()
{
  delete [] nodes;
  delete [] propedg;
  delete [] propsurf;
}



/**
  This method allocates array nodes depending on the element type.
  Element type is given by data member type.

  @param alloc_bp - controls allocation of boundary properties arrays on elements 
                    (1 = will be allocated/ 0 = will not be allocated)

  @retval 0 : on success allocation
  @retval 1 : unknown type is required


  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
long selement::alloc(long alloc_bp)
{
  switch (type)
  {
    case isolinear1d:
      nne     = 2;
      ned     = 1;
      nsurf   = 0;
      tnned   = 2;
      tnnsurf = 0;
      nned    = nned_isolin1d;
      edgenod = edgenod_isolin1d;
      nnsurf  = nnsurf_isolin1d;
      surfnod = surfnod_isolin1d;
      surfedg = surfedg_isolin1d;
      break;
    case isoquadratic1d:
      nne     = 3;
      ned     = 1;
      nsurf   = 0;
      tnned   = 3;
      tnnsurf = 0;
      nned    = nned_isoquad1d;
      edgenod = edgenod_isoquad1d;
      nnsurf  = nnsurf_isoquad1d;
      surfnod = surfnod_isoquad1d;
      surfedg = surfedg_isoquad1d;
      break;
    case trianglelinear:
      nne     = 3;
      ned     = 3;
      nsurf   = 1;
      tnned   = 6;
      tnnsurf = 3; 
      nned    = nned_trlin2d;
      edgenod = edgenod_trlin2d;
      nnsurf  = nnsurf_trlin2d;
      surfnod = surfnod_trlin2d;
      surfedg = surfedg_trlin2d;
      break;
    case trianglequadratic:
      nne     = 6;
      ned     = 3;
      nsurf   = 1;
      tnned   = 9;
      tnnsurf = 6;
      nned    = nned_trquad2d;
      edgenod = edgenod_trquad2d;
      nnsurf  = nnsurf_trquad2d;
      surfnod = surfnod_trquad2d;
      surfedg = surfedg_trquad2d;
      break;
    case isolinear2d:
      nne     = 4;
      ned     = 4;
      nsurf   = 1;
      tnned   = 8;
      tnnsurf = 4;
      nned    = nned_isolin2d;
      edgenod = edgenod_isolin2d;
      nnsurf  = nnsurf_isolin2d;
      surfnod = surfnod_isolin2d;
      surfedg = surfedg_isolin2d;
      break;
    case isoquadratic2d:
      nne     = 8;
      ned     = 4;
      nsurf   = 1;
      tnned   = 12;
      tnnsurf = 8;
      nned    = nned_isoquad2d;
      edgenod = edgenod_isoquad2d;
      nnsurf  = nnsurf_isoquad2d;
      surfnod = surfnod_isoquad2d;
      surfedg = surfedg_isoquad2d;
      break;
    case tetrahedronlinear:
      nne     = 4;
      ned     = 6;
      nsurf   = 4;
      tnned   = 12;
      tnnsurf = 12;
      nned    = nned_tetlin3d;
      edgenod = edgenod_tetlin3d;
      nnsurf  = nnsurf_tetlin3d;
      surfnod = surfnod_tetlin3d;
      surfedg = surfedg_tetlin3d;
      break;
    case tetrahedronquadratic:
      nne     = 10;
      ned     = 6;
      nsurf   = 4;
      tnned   = 18;
      tnnsurf = 24;
      nned    = nned_tetquad3d;
      edgenod = edgenod_tetquad3d;
      nnsurf  = nnsurf_tetquad3d;
      surfnod = surfnod_tetquad3d;
      surfedg = surfedg_tetquad3d;
      break;
    case pyramidelinear:
      nne     = 5;
      ned     = 8;
      nsurf   = 5;
      tnned   = 16;
      tnnsurf = 16;
      nned    = nned_pyramlin;
      edgenod = edgenod_pyramlin;
      nnsurf  = nnsurf_pyramlin;
      surfnod = surfnod_pyramlin;
      surfedg = surfedg_pyramlin;
      break;
    case pyramidequadratic:
      nne     = 13;
      ned     = 9;
      nsurf   = 5;
      tnned   = 27;
      tnnsurf = 32;
      nned    = nned_pyramquad;
      edgenod = edgenod_pyramquad;
      nnsurf  = nnsurf_pyramquad;
      surfnod = surfnod_pyramquad;
      surfedg = surfedg_pyramquad;
      break;
    case wedgelinear:
      nne     = 6;
      ned     = 9;
      nsurf   = 5;
      tnned   = 18;
      tnnsurf = 18;
      nned    = nned_wedgelin;
      edgenod = edgenod_wedgelin;
      nnsurf  = nnsurf_wedgelin;
      surfnod = surfnod_wedgelin;
      surfedg = surfedg_wedgelin;
      break;
    case wedgequadratic:
      nne     = 15;
      ned     = 9;
      nsurf   = 5;
      tnned   = 27;
      tnnsurf = 36;
      nned    = nned_wedgequad;
      edgenod = edgenod_wedgequad;
      nnsurf  = nnsurf_wedgequad;
      surfnod = surfnod_wedgequad;
      surfedg = surfedg_wedgequad;
      break;
    case isolinear3d:
      nne     = 8;
      ned     = 12;
      nsurf   = 6;
      tnned   = 24;
      tnnsurf = 24;
      nned    = nned_isolin3d;
      edgenod = edgenod_isolin3d;
      nnsurf  = nnsurf_isolin3d;
      surfnod = surfnod_isolin3d;
      surfedg = surfedg_isolin3d;
      break;
    case isoquadratic3d:
      nne     = 20;
      ned     = 12;
      nsurf   = 6;
      tnned   = 36;
      tnnsurf = 48;
      nned    = nned_isoquad3d;
      edgenod = edgenod_isoquad3d;
      nnsurf  = nnsurf_isoquad3d;
      surfnod = surfnod_isoquad3d;
      surfedg = surfedg_isoquad3d;
      break;
    case isocubic2d:
      nne     = 12;
      ned     = 4;
      nsurf   = 1;
      tnned   = 48;
      tnnsurf = 12;
      nned    = nned_isocubic2d;
      edgenod = edgenod_isocubic2d;
      nnsurf  = nnsurf_isocubic2d;
      surfnod = surfnod_isocubic2d;
      surfedg = surfedg_isocubic2d;
      break;
    default:
      return 1;
  }
  nodes    = NULL;
  propedg  = NULL;
  propsurf = NULL;
  if (nne)
  {
    nodes = new long[nne];
    memset(nodes, 0, sizeof(*nodes)*nne);
  }
  if (ned && alloc_bp)
  {
    propedg = new long[ned];
    memset(propedg, 0, sizeof(*propedg)*ned);
  }
  if (nsurf && alloc_bp)
  {
    propsurf = new long[nsurf];
    memset(propsurf, 0, sizeof(*propsurf)*nsurf);
  }
  return 0;
}



/**
  The function copies the given element to the argument instance dest.

  @param[out] dest - destination where the given element will be copied to.
  
  Created by Tomas Koudelka 01.2024, tomas.koudelka@fsv.cvut.cz
*/
void selement::copyto(selement &dest) const
{
  dest.type = type;
  dest.nne = nne;
  if (nodes){
    delete [] dest.nodes;
    dest.nodes = new long[nne];
    memcpy(dest.nodes, nodes, sizeof(*nodes)*nne);
  }
  else
    dest.nodes = NULL;

  dest.ned      = ned;
  dest.tnned    = tnned;
  dest.nned     = nned;
  dest.edgenod  = edgenod;

  if (propedg){
    delete [] dest.propedg;
    dest.propedg  = new long[ned];
    memcpy(dest.propedg, propedg, sizeof(*propedg)*ned);
  }
  else
    dest.propedg = NULL;

  
  dest.nsurf   = nsurf;
  dest.tnnsurf = tnnsurf;
  dest.nnsurf  = nnsurf;
  dest.surfnod = surfnod;

  if (propsurf){
    delete [] dest.propsurf;
    dest.propsurf = new long[nsurf];
    memcpy(dest.propsurf, propsurf, sizeof(*propsurf)*nsurf);
  }
  else
    dest.propsurf = NULL;  
 
  dest.prop = prop;
}



/**
  This method searches array of properties for property aprop with given entity type ent.
  If the pointers edg and sfc are not NULL then the additional edges edg or surfaces sfc are
  also searched for the given property.

  @param prop   - searched property number 
  @param ent    - general entity type whose property is searched
  @param tnodes - array of all nodes in the topology - it is used only in connection with vertex entity type
  @param eid    - element id (it is used in connection with the additional edg and sfc data)
  @param edg    - pointer to structure with data about additional edges
  @param sfc    - pointer to structure with data about additional surfaces
  @param entid  - pointer to the array of indices of given entity with given property on the element.
                  It is output parameter used only in connection with additional edges and surfaces
  @param nentid - the length of array entid (output parameter)

  @retval -1 - no property of edges or surfaces was given for the mesh
  @retval  0 - searched property was not found on given entity type
  @retval  1 - searched property was found

  Created by Tomas Koudelka 06.2009, tomas.koudelka@fsv.cvut.cz
*/
long selement::searchprop(long aprop, gentity ent, snode *tnodes, long eid, sedges *edg, sfaces *sfc, long *&entid, long &nentid) const
{
  long i, j, ret = 0;
  entid = NULL;
  nentid = 0;

  switch(ent)
  {
    case gvertex:
      for (i=0; i<nne; i++)
      {
        ret = tnodes[nodes[i]].searchprop(aprop, ent, nodes[i], edg, sfc);
        if (ret)
          return 1;
      }
      break;
    case gcurve:
      if (propedg)
      {
        for (i=0; i<ned; i++)
        {
          if (propedg[i] == aprop)
            return 1; 
        }
      }
      else
        ret=-1;
      break;
    case gsurface:
      if (propsurf)
      {
        for (i=0; i<nsurf; i++)
        {
          if (propsurf[i] == aprop)
            return 1; 
        }
      }
      else
        if (nsurf > 0)  ret=-1; // element is not 1D
      break;
    case gregion: 
      if (prop == aprop)   return 1;
      break;
    default:
     print_err("Unknown type of entity is required", __FILE__, __LINE__, __func__);
     abort();
  }

  if ((ent == gcurve) && (edg))
  {
    ret = 0; // edge properties was read
    // search additional edges specified in the MeshEditor
    for(i=0; i<edg->nprop; i++)
    {
      if (edg->proplist[i] == aprop)
      {
        for(j=0; j<edg->prop_list_elem_num[i]; j++)
        {
          if (edg->prop_list_elem[i][j] == eid)
          {
            entid = edg->prop_list_elem_edgid[i][j];
            nentid = edg->prop_list_elem_edgid_num[i][j]; 
            return 1;
          }
          if (edg->prop_list_elem[i][j] > eid)    break;  // element numbers are sorted ascendently in the prop_list_elem
        }
      }
    }
  }
  
  if ((ent == gsurface) && (sfc))
  {
    ret = 0; // surface properties was read
    // search additional surfaces specified in the MeshEditor
    for(i=0; i<sfc->nprop; i++)
    {
      if (sfc->proplist[i] == aprop)
      {
        for(j=0; j<sfc->prop_list_elem_num[i]; j++)
        {
          if (sfc->prop_list_elem[i][j] == eid)
          {
            entid = sfc->prop_list_elem_sfid[i][j];
            nentid = sfc->prop_list_elem_sfid_num[i][j]; 
            return 1;
          }
          if (sfc->prop_list_elem[i][j] > eid)    break;  // element numbers are sorted ascendently in the prop_list_elem
        }
      }
    }
  }
  return ret;
}



/**
  This method searches array of properties of particular nodes for property aprop with given entity type ent.
  If the property on the i-th element node is matched then the aprop is stored in the output parameter propent[i].

  @param prop    - searched property number 
  @param ent     - general entity type whose property is searched
  @param tnodes  - array of all nodes in the topology - it is used only in connection with vertex entity type
  @param edg     - pointer to structure with data about additional edges
  @param sfc     - pointer to structure with data about additional surfaces
  @param propent - array of nodal properties matched 
                   propent[i] = aprop if the property of the i-th element node matches aprop of the given enitity ent,
                   propent[i] = -1 for the otherwise cases

  @retval  0 - searched property was not found for the given entity type
  @retval  1 - searched property was found

  Created by Tomas Koudelka 07.2014, tomas.koudelka@fsv.cvut.cz
*/
long selement::searchnodprop(long aprop, gentity ent, snode *tnodes, sedges *edg, sfaces *sfc, long *propent) const
{
  long i, aux, ret = 0;
  for (i=0; i<nne; i++)
  {
    propent[i] = -1;
    // searches one element node for the given property
    aux = tnodes[nodes[i]].searchprop(aprop, ent, nodes[i], edg, sfc);
    if (aux)
    {
      propent[i] = aprop;
      ret = 1;
    }
  }
  return ret;
}



/**
  The method compares array of nodes nod with node numbers on particular element edges.
  If they corresponds with some of the element edge nodes than the edge id is returned.
  Otherwise the -1 is returned.
  

  @param nod - array of searched nodes
  @param n  - number of searched nodes (length of nod)

  @retval  edgid - searched nodes corresponded with edge edgid
  @retval -1 - no corresondence of the nodes was found

  Created by Tomas Koudelka 08.2011, tomas.koudelka@fsv.cvut.cz
*/
long selement::compare_edg(const long *nod, long n) const
{
  long i, j, k;
  long nf; // number of found nodes from the array nod
  long stop;

  for (i=0; i<ned; i++)
  {
    nf = 0; 
    for(j=0; j<n; j++)
    {
      stop = 1;
      for(k=0; k<nned[i]; k++)
      {
        if (nodes[edgenod[i][k]] == nod[j])
        {
          nf++;
          stop = 0;
          break;
        }        
      }
      if (stop)
        break;
    }
    if (nf == n)
      return i;
  }
  return -1;
}



/**
  The method compares array of nodes nod with node numbers on particular element edges.
  If all nodes of the particular element edge were found in nod array than the edge 
  id is stored in the list of found edge identifiers edgid_lst. The function return
  number of found edges.
  

  @param[in]  nod - pointer to array of searched nodes
  @param[in]  n - the number of searched nodes, i.e. dimension of the array nod
  @param[out] edgid_lst - list of identifiers of element edges whose nodes were found in nod
  

  @return  The function returns the number of found element edges whose nodes were 
           found in the nod array.

  Created by Tomas Koudelka 01.2024, tomas.koudelka@fsv.cvut.cz
*/
long selement::compare_edg(const long *nod, long n, ivector &edgid_lst) const
{
  long i, j, k, l=0;
  long nf; // number of found nodes from the array nod

  edgid_lst.n = 0;
  for (i=0; i<ned; i++){
    nf = 0; 
    for(j=0; j<nned[i]; j++){
      for(k=0; k<n; k++){
        if (nodes[edgenod[i][j]] == nod[k]){
          nf++;
          break;
        }        
      }
    }
    if (nf == nned[i]){
      edgid_lst[l] = i;
      l++;
      edgid_lst.n = l;
    }
  }
  return l;
}



/**
  The method compares array of nodes nod with node numbers on particular element surfaces.
  If they correspond with some of the element surface nodes than the surface id is returned.
  Otherwise the -1 is returned.
  

  @param nod - array of searched nodes
  @param n  - number of searched nodes (length of nod)

  @retval  surfid - searched nodes corresponded with surface surfid
  @retval -1 - no corresondence of the nodes was found

  Created by Tomas Koudelka 08.2011, tomas.koudelka@fsv.cvut.cz
*/
long selement::compare_surf(const long *nod, long n) const
{
  long i, j, k;
  long nf; // number of found nodes from the array nod
  long stop;

  for (i=0; i<nsurf; i++)
  {
    nf = 0;
    for(j=0; j<n; j++)
    {
      stop = 1;
      for(k=0; k<nnsurf[i]; k++)
      {
        if (nodes[surfnod[i][k]] == nod[j])
        {
          nf++;
          stop = 0;
          break;
        }        
      }
      if (stop)
        break;
    }
    if (nf == n)
      return i;
  }
  return -1;
}



/**
  The method compares array of nodes nod with node numbers on particular element surfaces.
  If all nodes of the particular element surface were found in nod array than the surface 
  id is stored in the list of found surface identifiers surfid_lst. The function return
  number of found surfaces.
  

  @param[in]  nod - pointer to array of searched nodes
  @param[in]  n - the number of searched nodes, i.e. dimension of the array nod
  @param[out] surfid_lst - list of identifiers of element surfaces whose nodes were found in nod
  

  @return  The function returns the number of found element surfaces whose nodes were 
           found in the nod array.

  Created by Tomas Koudelka 01.2024, tomas.koudelka@fsv.cvut.cz
*/
long selement::compare_surf(const long *nod, long n, ivector &surfid_lst) const
{
  long i, j, k, l=0;
  long nf; // number of found nodes from the array nod

  surfid_lst.n = 0;
  for (i=0; i<nsurf; i++){
    nf = 0;
    for(j=0; j<nnsurf[i]; j++){
      for(k=0; k<n; k++){
        if (nodes[surfnod[i][j]] == nod[k]){
          nf++;
          break;
        }
      }
    }
    if (nf == nnsurf[i]){
      surfid_lst[l] = i;
      l++;
      surfid_lst.n = l;
    }
  }
  return l;
}



/**
  The constructor initializes all data members to -1.

  Created by Tomas Koudelka, 08.2011
*/
sedge::sedge(void)
{
  n1 = n2 = -1;
  prop = -1;
}



/**
  The function copies the given instance to the argument instance dest.

  @param[out] dest - destination instance where the data of the given instance will be copied to.

  Created by TKo, 01.2024
*/
void sedge::copyto(sedge &dest)
{
  dest.n1 = n1;
  dest.n2 = n2;
  dest.prop = prop;
}



/**
  The function reads data about one edge from the opened text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - in the case of an read error

  Created by Tomas Koudelka, 08.2011
*/
long sedge::read(XFILE *in)
{
  if (xfscanf(in, "%ld %ld %ld", &n1, &n2, &prop) < 3)
    return 1;
  // change to C/C++ index notation
  n1--;
  n2--;
  return 0;
}



/**
  The function prints data about one edge to the opened text file.
  If the second parameter lnn is not NULL, 
  then the node numbers on the edge are renumbered using this array.
  The lnn array is used in connection with the mesh decompoistion and it
  contains local node numbers for a domain. The size of the array is given
  by the total number of nodes in the problem.

  @param in - pointer to the opened text file
  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sedge::print(FILE *out, long *lnn)
{
  if (lnn)
    fprintf(out, "%ld %ld %ld\n", lnn[n1], lnn[n2], prop);
  else
    fprintf(out, "%ld %ld %ld\n", n1+1, n2+1, prop);
}



/**
  The function tests whether all nodes of the edge are contained
  in the given subdomain whose local node numbers given in 
  the parameter lnn. The lnn array is used in connection with
  the mesh decompoistion and it contains local node numbers for a domain. The size 
  of the array is given by the total number of nodes in the problem.

  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @retval 0 - edge is not contained in the given subdomain
  @retval 1 - edge is contained in the given subdomain

  Created by Tomas Koudelka, 09.2011
*/
long sedge::cordomtest(long *lnn)
{
  if (lnn[n1] && lnn[n2])
    return 1;

  return 0;
}



/**
  The constructor initializes all data members to zero or NULL values.

  Created by Tomas Koudelka, 08.2011
*/
sedges::sedges()
{
  nedg = 0L;
  nprop = 0L;
  edges = NULL;
  proplist = NULL;
  prop_list_nod = NULL;
  prop_list_nod_num = NULL;
  prop_list_elem = NULL;
  prop_list_elem_num = NULL;
  prop_list_elem_edgid = NULL;
  prop_list_elem_edgid_num = NULL;
}



/**
  The destructor releases all allocated memory.

  Created by Tomas Koudelka, 08.2011
*/
sedges::~sedges()
{
  long i, j;

  delete [] edges;
  for (i=0; i<nprop; i++)
  {
    for (j=0; j<prop_list_elem_num[i]; j++)
      delete [] prop_list_elem_edgid[i][j];
    
    delete [] prop_list_nod[i];
    delete [] prop_list_elem[i];
    delete [] prop_list_elem_edgid[i];
    delete [] prop_list_elem_edgid_num[i];
  }
  if (nprop)
  {
    delete [] proplist;
    delete [] prop_list_nod;
    delete [] prop_list_nod_num;
    delete [] prop_list_elem;
    delete [] prop_list_elem_num;
    delete [] prop_list_elem_edgid;
    delete [] prop_list_elem_edgid_num;
  }
}



/**
  The constructor allocates memory for the list of edges
  and initializes rest of data members to zero or NULL values.
  
  @param n - number of allocated edges

  Created by Tomas Koudelka, 08.2011
*/
sedges::sedges(long n)
{
  nedg = n;
  edges = new sedge[n];
  nprop = 0L;
  proplist = NULL;
  prop_list_nod = NULL;
  prop_list_nod_num = NULL;
  prop_list_elem = NULL;
  prop_list_elem_num = NULL;
  prop_list_elem_edgid = NULL;
  prop_list_elem_edgid_num = NULL;
}



/**
  The function copies the given instance to the argument instance dest.

  @param[out] dest - destination where the given instance will be copied to.

  Created by TKo, 01.2024
*/
void sedges::copyto(sedges &dest)
{  
  dest.nedg = nedg;
  if (edges){
    delete [] dest.edges;
    dest.edges = new sedge[nedg];
    for (long i=0; i<nedg; i++)
      edges[i].copyto(dest.edges[i]);    
  }
  else
    edges = NULL;
  
  for (long i=0; i<dest.nprop; i++){
    for (long j=0; j<dest.prop_list_elem_num[i]; j++)
      delete [] dest.prop_list_elem_edgid[i][j];
    
    delete [] dest.prop_list_nod[i];
    delete [] dest.prop_list_elem[i];
    delete [] dest.prop_list_elem_edgid[i];
    delete [] dest.prop_list_elem_edgid_num[i];
  }
  if (dest.nprop){
    delete [] dest.proplist;
    delete [] dest.prop_list_nod;
    delete [] dest.prop_list_nod_num;
    delete [] dest.prop_list_elem;
    delete [] dest.prop_list_elem_num;
    delete [] dest.prop_list_elem_edgid;
    delete [] dest.prop_list_elem_edgid_num;
  }

  if (nprop){
    dest.nprop = nprop;

    dest.proplist = new long[nprop];
    memcpy(dest.proplist, proplist, sizeof(*proplist)*nprop);
    
    dest.prop_list_nod_num = new long[nprop];
    memcpy(dest.prop_list_nod_num, prop_list_nod_num, sizeof(*prop_list_nod_num)*nprop);

    dest.prop_list_elem_num = new long[nprop];
    memcpy(dest.prop_list_elem_num, prop_list_elem_num, sizeof(*prop_list_elem_num)*nprop);

    dest.prop_list_nod  = new long*[nprop];
    dest.prop_list_elem = new long*[nprop];
    dest.prop_list_elem_edgid_num = new long*[nprop];
    dest.prop_list_elem_edgid     = new long**[nprop];
    for (long i=0; i<nprop; i++){
      dest.prop_list_nod[i] = new long[prop_list_nod_num[i]];
      memcpy(dest.prop_list_nod[i], prop_list_nod[i], sizeof(*prop_list_nod[i])*prop_list_nod_num[i]);
      
      dest.prop_list_elem[i] = new long[prop_list_elem_num[i]];
      memcpy(dest.prop_list_elem[i], prop_list_elem[i], sizeof(*prop_list_elem[i])*prop_list_elem_num[i]);

      dest.prop_list_elem_edgid_num[i] = new long[prop_list_elem_num[i]];
      memcpy(dest.prop_list_elem_edgid_num[i], prop_list_elem_edgid_num[i], sizeof(*prop_list_elem_edgid_num[i])*prop_list_elem_num[i]);
      
      dest.prop_list_elem_edgid[i] = new long*[prop_list_elem_num[i]];
      for (long j=0; j<prop_list_elem_num[i]; j++){
        if (prop_list_elem_edgid_num[i][j]){
          dest.prop_list_elem_edgid[i][j] = new long[prop_list_elem_edgid_num[i][j]];
          memcpy(dest.prop_list_elem_edgid[i][j], prop_list_elem_edgid[i][j], sizeof(*dest.prop_list_elem_edgid[i][j])*dest.prop_list_elem_edgid_num[i][j]);
        }
        else
          dest.prop_list_elem_edgid[i][j] = NULL;
      }
    }
  }
  else{
    dest.nprop = 0L;
    dest.proplist = NULL;
    dest.prop_list_nod = NULL;
    dest.prop_list_nod_num = NULL;
    dest.prop_list_elem = NULL;
    dest.prop_list_elem_num = NULL;
    dest.prop_list_elem_edgid = NULL;
    dest.prop_list_elem_edgid_num = NULL;
  }
}



/**
  The function reads list of additional edges from the opened text file.
  The list is placed at the end of file with topology and it is created by
  MeshEditor usually.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sedges::read(XFILE *in)
{
  long i;
  for (i=0; i<nedg; i++)
  {
    if (edges[i].read(in))
    {
      print_err("cannot read edge number %ld", __FILE__, __LINE__, __func__, i+1);
      abort();
    }
  }
}



/**
  The function prints list of additional edges from the opened text file.
  The list is placed the end of file with topology usually in MeshEditor format.
  If the second parameter lnn is not NULL, then only the edges whose all local node numbers
  are nonzero are printed. The lnn array is used in connection with
  the mesh decompoistion and it contains local node numbers for a domain. The size 
  of the array is given by the total number of nodes in the problem.

  @param out - pointer to the opened text file
  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sedges::print(FILE *out, long *lnn)
{
  long i;
  long *aux, nedg_l;

  if (lnn == NULL) // no mesh decomposition
  {
    fprintf(out, "edges %ld\n", nedg);
    for (i=0; i<nedg; i++)
      edges[i].print(out, lnn);
  }
  else
  {
    nedg_l=0;
    aux = new long[nedg];
    memset(aux, 0, sizeof(*aux)*nedg);

    for(i=0; i<nedg; i++)
    {      
      if (edges[i].cordomtest(lnn))
      {
        aux[i]++;
        nedg_l++;
      }
    }
    fprintf(out, "edges %ld\n", nedg_l);
    for(i=0; i<nedg; i++)
    {
      if (aux[i])
        edges[i].print(out, lnn);
    }
    delete [] aux;
  }
}



/**
  The function searches the list of property numbers for the given number aprop.
  It returns whether the aprop is a new property number.

  @param aprop - searched property number
  @param propid - array of with list of property numbers
  @param n      - length of array propid
  
  @retval 0 - the property was found in the list
  @retval 1 - the property was not found in the list

  Created by Tomas Koudelka, 08.2011
*/
long sedges::search_list_newprop(long aprop, long *propid, long n)
{
  long i;

  for(i=0; i<n; i++)
  {
    if (propid[i] == aprop)
      return 0;
  }
  return 1;
}



/**
  The function generates from the list of additional edges new lists which contain 
  the list of different property numbers, lists of nodes for particular property numbers, 
  lists of elements for particular property numbers and lists of edge id for particular
  lists of elements with given property number.

  @param nn - total number of nodes in the mesh
  @param ne - total number of elements in the mesh
  @param elems - array of all elements in the mesh
  @param nadjelnod - array of number of adjacent elements to nodes
  @param adjelnod - array of pointers to arrays of adjacent elements ids to nodes

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sedges::gen_list_edg(long nn, long ne, selement *elems, long *nadjelnod, long **adjelnod)
{
  long i;
  long aux;
  long *propid;


  if (nedg > 0)
  {
    propid = new long[nedg];
    nprop = 0;
    propid[nprop] = aux = edges[0].prop;
    nprop++;
    for(i=0; i<nedg; i++)
    {
      if (edges[i].prop != aux)
      { 
        aux = edges[i].prop; 
        if (search_list_newprop(aux, propid, nprop))
        {
          propid[nprop] = aux;
          nprop++;
        }
      }
    }    
    proplist = new long[nprop];
    for(i=0; i<nprop; i++)
      proplist[i] = propid[i];
    delete [] propid;
    prop_list_nod       = new long*[nprop];
    prop_list_nod_num   = new long[nprop];
    prop_list_elem      = new long*[nprop];
    prop_list_elem_num  = new long[nprop];
    prop_list_elem_edgid = new long**[nprop];
    prop_list_elem_edgid_num = new long*[nprop];
    memset(prop_list_nod, 0, sizeof(*prop_list_nod)*nprop);
    memset(prop_list_nod_num, 0, sizeof(*prop_list_nod_num)*nprop);
    memset(prop_list_elem, 0, sizeof(*prop_list_elem)*nprop);
    memset(prop_list_elem_num, 0, sizeof(*prop_list_elem_num)*nprop);
    memset(prop_list_elem_edgid, 0, sizeof(*prop_list_elem_edgid)*nprop);
    memset(prop_list_elem_edgid_num, 0, sizeof(*prop_list_elem_edgid_num)*nprop);

    gen_list_node(nn);
    gen_list_elem(ne, elems, nadjelnod, adjelnod);
  }
}



/**
  The function generates from the list of additional edges new list which contain 
  lists of nodes for particular different property numbers. The lists are stored
  in the array prop_list_nod where prop_list_nod[i] represents array of
  node numbers for the property number proplist[i]. The length of array
  prop_list_nod[i] is given by prop_list_nod_num[i].

  @param nn - total number of nodes in the mesh

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sedges::gen_list_node(long nn)
{
  long i, j, k, nnodp;
  long *nodid;

  nodid = new long[nn];

  for(i=0; i<nprop; i++)
  {
    memset(nodid, 0, sizeof(*nodid)*nn);   
    nnodp = 0;
    for(j=0; j<nedg; j++)
    {
      if (edges[j].prop == proplist[i])
      {
        if (nodid[edges[j].n1] == 0)
        {
          nodid[edges[j].n1] = 1;
          nnodp++;
        }
        if (nodid[edges[j].n2] == 0)
        {
          nodid[edges[j].n2] = 1;
          nnodp++;
        }
      }
    }
    prop_list_nod_num[i] = nnodp;
    prop_list_nod[i] = new long[nnodp];
    k = 0;
    for(j=0; j<nn; j++)
    {
      if (nodid[j])
      {
        prop_list_nod[i][k] = j;
        k++;
      }
    }
  }
  delete [] nodid;
}



/**
  The function generates from the list of additional edges new lists which contain 
  lists of elements for particular different property numbers. The lists are stored
  in the array prop_list_nod where prop_list_elem[i] represents array of
  element numbers for the property number proplist[i]. The length of array
  prop_list_elem[i] is given by prop_list_elem_num[i]. 

  Additionally, the function  creates lists of edge indices with given property number 
  for each list of elements. Edge indices are stored in the 
  array prop_list_elem_edgid where the prop_list_elem_edgid[i][j] represents array of 
  edge indices with property number proplist[i] for the j-th element from the prop_list_elem[i].
  The length of array  prop_list_elem_edgid[i] is given by prop_list_elem_edgid_num[i][j].    

  @param ne - total number of elements in the mesh
  @param elems - array of all elements in the mesh
  @param nadjelnod - array of number of adjacent elements to nodes
  @param adjelnod - array of pointers to arrays of adjacent elements ids to nodes

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
/*void sedges::gen_list_elem(long ne, selement *elems)
{
  long i, j, k, l, nelemp;
  long *elemid;
  long **elemedgid;
  long edgid;
  long nodes[2];
  long aux;

  elemid    = new long[ne];
  elemedgid = new long*[ne];
  for(i=0; i<ne; i++)
    elemedgid[i] = new long [elems[i].ned];

  for(i=0; i<nprop; i++)
  {
    memset(elemid, 0, sizeof(*elemid)*ne);
    for(j=0; j<ne; j++)
      memset(elemedgid[j], 0, sizeof(*elemedgid[j])*elems[j].ned);
    nelemp = 0;

    for(j=0; j<nedg; j++)
    {
      if (edges[j].prop == proplist[i])
      {        
        for(k=0; k<ne; k++)
        {
          nodes[0] = edges[j].n1;
          nodes[1] = edges[j].n2;
          edgid = elems[k].compare_edg(nodes, 2);
          if (edgid >= 0)
          {
            if (elemid[k] == 0)
            {
              nelemp++;
              elemid[k] = 1;
            }
            elemedgid[k][edgid] = 1;
          }
        }
      }
    }
    prop_list_elem_num[i]      = nelemp;
    prop_list_elem[i]          = new long[nelemp];
    prop_list_elem_edgid_num[i] = new long[nelemp];
    prop_list_elem_edgid[i]     = new long*[nelemp];
    k = 0;
    for(j=0; j<ne; j++)
    {
      if (elemid[j] == 0)
        continue;

      prop_list_elem[i][k] = j;
      aux = 0;
      for(l = 0; l<elems[j].ned; l++)
      {
        if (elemedgid[j][l])
          aux++;
      }
      prop_list_elem_edgid_num[i][k] = aux;
      if (aux)
        prop_list_elem_edgid[i][k] = new long[aux];
      else
        prop_list_elem_edgid[i][k] = NULL;
      aux = 0;
      for(l = 0; l<elems[j].ned; l++)
      {
        if (elemedgid[j][l])
        {
          prop_list_elem_edgid[i][k][aux] = l;
          aux++;
        }
      }
      k++;
    }
  }
  delete [] elemid;
  for(i=0; i<ne; i++)
    delete [] elemedgid[i];
  delete [] elemedgid;
}*/
void sedges::gen_list_elem(long ne, selement *elems, long *nadjelnod, long **adjelnod)
{
  long i, j, k, l, nelemp;
  long aeid1, aeid2;
  long *elemid;
  long **elemedgid;
  long edgid;
  long nodes[2];
  long aux;

  elemid    = new long[ne];
  elemedgid = new long*[ne];
  for(i=0; i<ne; i++)
    elemedgid[i] = new long [elems[i].ned];

  for(i=0; i<nprop; i++)
  {
    memset(elemid, 0, sizeof(*elemid)*ne);
    for(j=0; j<ne; j++)
      memset(elemedgid[j], 0, sizeof(*elemedgid[j])*elems[j].ned);
    nelemp = 0;

    for(j=0; j<nedg; j++)
    {
      if (edges[j].prop == proplist[i])
      {             
        // compare adjacent element ids of m-th edge node with rest elements of 
        // remaining face nodes
        for(k=0; k<nadjelnod[edges[j].n1]; k++) // index of adjacent elements to m-th node of j-th edge
        {
          aeid1 = adjelnod[edges[j].n1][k];
          aux=1; // number of occurrences of the element in adjacent elements of the second node
          for(l=0; l<nadjelnod[edges[j].n2]; l++)
          {
            aeid2 = adjelnod[edges[j].n2][l];
            if (aeid2 == aeid1)
            {
              aux++;
              break;
            }
            if (aeid1 < aeid2) // adjelnod arrays are sorted ascending -> remaining numbers are greater than aeid_min
              break;
          }
          if (aux==2)
          {
            nodes[0] = edges[j].n1;
            nodes[1] = edges[j].n2;
            edgid = elems[aeid1].compare_edg(nodes, 2);
            if (edgid >= 0)
            {
              if (elemid[aeid1] == 0)
              {
                nelemp++;
                elemid[aeid1] = 1;
              }
              elemedgid[aeid1][edgid] = 1;
            }
          }
        }
      }
    }
    prop_list_elem_num[i]      = nelemp;
    prop_list_elem[i]          = new long[nelemp];
    prop_list_elem_edgid_num[i] = new long[nelemp];
    prop_list_elem_edgid[i]     = new long*[nelemp];
    k = 0;
    for(j=0; j<ne; j++)
    {
      if (elemid[j] == 0)
        continue;

      prop_list_elem[i][k] = j;
      aux = 0;
      for(l = 0; l<elems[j].ned; l++)
      {
        if (elemedgid[j][l])
          aux++;
      }
      prop_list_elem_edgid_num[i][k] = aux;
      if (aux)
        prop_list_elem_edgid[i][k] = new long[aux];
      else
        prop_list_elem_edgid[i][k] = NULL;
      aux = 0;
      for(l = 0; l<elems[j].ned; l++)
      {
        if (elemedgid[j][l])
        {
          prop_list_elem_edgid[i][k][aux] = l;
          aux++;
        }
      }
      k++;
    }
  }
  delete [] elemid;
  for(i=0; i<ne; i++)
    delete [] elemedgid[i];
  delete [] elemedgid;
}


/**
  The constructor initializes all data members to -1, 0 or NULL values.

  Created by Tomas Koudelka, 08.2011
*/
sface::sface(void)
{
  nnod = 0;
  nodes = NULL;
  prop = -1;
}



/**
  The constructor initializes all data members to -1 or 0 values.
  It allocates memory for n nodes.

  @param - the number of nodes defining one face

  Created by Tomas Koudelka, 08.2011
*/
sface::sface(long n)
{
  nnod = n;
  nodes = new long[n];
  memset(nodes, 0, sizeof(*nodes)*n);
  prop = -1;
}



/**
  The destructor releases allocated memory.

  Created by Tomas Koudelka, 08.2011
*/
sface::~sface(void)
{
  delete [] nodes;
}



/**
  The function copies the given instance to the argument instance dest.

  @param[out] dest - instance where the given instance data will be copied to.

  Created by Tomas Koudelka, 01.2024
*/
void sface::copyto(sface &dest)
{
  dest.nnod = nnod;

  if (nodes){
    delete [] dest.nodes;
    dest.nodes = new long[nnod];
    memcpy(dest.nodes, nodes, sizeof(*nodes)*nnod);
  }
  else
    dest.nodes = NULL;  

  dest.prop = prop;
}



/**
  The function reads data about one face from the opened text file.

  @param in  - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - in the case of an read error

  Created by Tomas Koudelka, 08.2011
*/
long sface::read(XFILE *in)
{
  long i;

  if (xfscanf(in, "%ld", &nnod) == 0)
    return 1;

  nodes = new long[nnod];
  memset(nodes, 0, sizeof(*nodes)*nnod);

  for(i=0; i<nnod; i++)
  {
    if (xfscanf(in, "%ld", nodes+i) == 0)
      return 1;
    // change to C/C++ index notation
    nodes[i]--;

  }
  if (xfscanf(in, "%ld", &prop) == 0)
    return 1;

  return 0;
}



/**
  The function prints data about one face to the opened text file.
  If the second parameter lnn is not NULL, 
  then the node numbers on the face are renumbered using this array.
  The lnn array is used in connection with the mesh decompoistion and it
  contains local node numbers for a domain. The size of the array is given
  by the total number of nodes in the problem.

  @param in - pointer to the opened text file
  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sface::print(FILE *out, long *lnn)
{
  long i;
  
  fprintf(out, "%ld", nnod);
  for(i=0; i<nnod; i++)
  {
    if (lnn)
      fprintf(out, " %ld", lnn[nodes[i]]);
    else
      fprintf(out, " %ld", nodes[i]+1);
  }

  fprintf(out, " %ld\n", prop);
}



/**
  The function tests whether all nodes of the face are contained
  in the given subdomain whose local node numbers given in 
  the parameter lnn. The lnn array is used in connection with
  the mesh decompoistion and it contains local node numbers for a domain. The size 
  of the array is given by the total number of nodes in the problem.

  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @retval 0 - face is not contained in the given subdomain
  @retval 1 - face is contained in the given subdomain

  Created by Tomas Koudelka, 09.2011
*/
long sface::cordomtest(long *lnn)
{
  long i, aux=0;

  for(i=0; i<nnod; i++)
  {
    if (lnn[nodes[i]])
     aux++;
  }
  if (aux == nnod)
    return 1;

  return 0;
}



/**
  The constructor initializes all data members to zero or NULL values.

  Created by Tomas Koudelka, 08.2011
*/
sfaces::sfaces()
{
  nfcs = 0L;
  nprop = 0L;
  faces = NULL;
  proplist = NULL;
  prop_list_nod = NULL;
  prop_list_nod_num = NULL;
  prop_list_elem = NULL;
  prop_list_elem_num = NULL;
  prop_list_elem_sfid = NULL;
  prop_list_elem_sfid_num = NULL;
}



/**
  The destructor releases all allocated memory.

  Created by Tomas Koudelka, 08.2011
*/
sfaces::~sfaces()
{
  long i, j;

  delete [] faces;
  delete [] proplist;
  for (i=0; i<nprop; i++)
  {
    for (j=0; j<prop_list_elem_num[i]; j++)
      delete [] prop_list_elem_sfid[i][j];
    
    delete [] prop_list_nod[i];
    delete [] prop_list_elem[i];
    delete [] prop_list_elem_sfid[i];
    delete [] prop_list_elem_sfid_num[i];
  }
  if (nprop)
  {
    delete [] prop_list_nod;
    delete [] prop_list_nod_num;
    delete [] prop_list_elem;
    delete [] prop_list_elem_num;
    delete [] prop_list_elem_sfid;
    delete [] prop_list_elem_sfid_num;
  }
}



/**
  The constructor allocates memory for the list of faces
  and initializes rest of data members to zero or NULL values.
  
  @param n - number of allocated edges

  Created by Tomas Koudelka, 08.2011
*/
sfaces::sfaces(long n)
{
  nfcs = n;
  nprop = 0L;
  faces = new sface[n];
  proplist = NULL;
  prop_list_nod = NULL;
  prop_list_nod_num = NULL;
  prop_list_elem = NULL;
  prop_list_elem_num = NULL;
  prop_list_elem_sfid = NULL;
  prop_list_elem_sfid_num = NULL;
}




/**
  The function copies the given instance to the argument instance dest.

  @param[out] dest - destination where the given instance will be copied to.

  Created by TKo, 01.2024
*/
void sfaces::copyto(sfaces &dest)
{
  dest.nfcs = nfcs;
  if (faces){
    delete [] dest.faces;
    dest.faces = new sface[nfcs];
    for (long i=0; i<nfcs; i++)
      faces[i].copyto(dest.faces[i]);    
  }
  else
    faces = NULL;
  
  for (long i=0; i<dest.nprop; i++){
    for (long j=0; j<dest.prop_list_elem_num[i]; j++)
      delete [] dest.prop_list_elem_sfid[i][j];
    
    delete [] dest.prop_list_nod[i];
    delete [] dest.prop_list_elem[i];
    delete [] dest.prop_list_elem_sfid[i];
    delete [] dest.prop_list_elem_sfid_num[i];
  }
  if (dest.nprop){
    delete [] dest.proplist;
    delete [] dest.prop_list_nod;
    delete [] dest.prop_list_nod_num;
    delete [] dest.prop_list_elem;
    delete [] dest.prop_list_elem_num;
    delete [] dest.prop_list_elem_sfid;
    delete [] dest.prop_list_elem_sfid_num;
  }

  if (nprop){
    dest.nprop = nprop;

    dest.proplist = new long[nprop];
    memcpy(dest.proplist, proplist, sizeof(*proplist)*nprop);
    
    dest.prop_list_nod_num = new long[nprop];
    memcpy(dest.prop_list_nod_num, prop_list_nod_num, sizeof(*prop_list_nod_num)*nprop);

    dest.prop_list_elem_num = new long[nprop];
    memcpy(dest.prop_list_elem_num, prop_list_elem_num, sizeof(*prop_list_elem_num)*nprop);

    dest.prop_list_nod  = new long*[nprop];
    dest.prop_list_elem = new long*[nprop];
    dest.prop_list_elem_sfid_num = new long*[nprop];
    dest.prop_list_elem_sfid     = new long**[nprop];
    for (long i=0; i<nprop; i++){
      dest.prop_list_nod[i] = new long[prop_list_nod_num[i]];
      memcpy(dest.prop_list_nod[i], prop_list_nod[i], sizeof(*prop_list_nod[i])*prop_list_nod_num[i]);
      
      dest.prop_list_elem[i] = new long[prop_list_elem_num[i]];
      memcpy(dest.prop_list_elem[i], prop_list_elem[i], sizeof(*prop_list_elem[i])*prop_list_elem_num[i]);

      dest.prop_list_elem_sfid_num[i] = new long[prop_list_elem_num[i]];
      memcpy(dest.prop_list_elem_sfid_num[i], prop_list_elem_sfid_num[i], sizeof(*prop_list_elem_sfid_num[i])*prop_list_elem_num[i]);
      
      dest.prop_list_elem_sfid[i] = new long*[prop_list_elem_num[i]];
      for (long j=0; j<prop_list_elem_num[i]; j++){
        if (prop_list_elem_sfid_num[i][j]){
          dest.prop_list_elem_sfid[i][j] = new long[prop_list_elem_sfid_num[i][j]];
          memcpy(dest.prop_list_elem_sfid[i][j], prop_list_elem_sfid[i][j], sizeof(*dest.prop_list_elem_sfid[i][j])*dest.prop_list_elem_sfid_num[i][j]);
        }
        else
          dest.prop_list_elem_sfid[i][j] = NULL;
      }
    }
  }
  else{
    dest.nprop = 0L;
    dest.proplist = NULL;
    dest.prop_list_nod = NULL;
    dest.prop_list_nod_num = NULL;
    dest.prop_list_elem = NULL;
    dest.prop_list_elem_num = NULL;
    dest.prop_list_elem_sfid = NULL;
    dest.prop_list_elem_sfid_num = NULL;
  }
}



/**
  The function reads list of additional surfaces from the opened text file.
  The list is placed at the end of file with topology and it is created by
  MeshEditor usually.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sfaces::read(XFILE *in)
{
  long i;

  for(i=0; i<nfcs; i++)
  {
    if (faces[i].read(in))
    {
      print_err("cannot read face number %ld\n", __FILE__, __LINE__, __func__, i+1);
      abort();
    }
  }
}



/**
  The function prints list of additional surfaces from the opened text file.
  The list is placed the end of file with topology usually in MeshEditor format.
  If the second parameter lnn is not NULL, then only the faces whose all local node numbers
  are nonzero are printed. The lnn array is used in connection with
  the mesh decompoistion and it contains local node numbers for a domain. The size 
  of the array is given by the total number of nodes in the problem.
  

  @param out - pointer to the opened text file
  @param lnn - pointer to array of local node numbers used for the
               mesh decomposition. Default value of this parameter is NULL.

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sfaces::print(FILE *out, long *lnn)
{
  long i;
  long *aux;
  long nfcs_l;
  
  if (lnn == NULL) // no mesh decomposition
  {
    fprintf(out, "faces %ld\n", nfcs);
    for(i=0; i<nfcs; i++)
      faces[i].print(out, lnn);
  }
  else
  {
    nfcs_l=0;
    aux = new long[nfcs];
    memset(aux, 0, sizeof(*aux)*nfcs);

    for(i=0; i<nfcs; i++)
    {      
      if (faces[i].cordomtest(lnn))
      {
        aux[i]++;
        nfcs_l++;
      }
    }
    fprintf(out, "faces %ld\n", nfcs_l);
    for(i=0; i<nfcs; i++)
    {
      if (aux[i])
        faces[i].print(out, lnn);
    }
    delete [] aux;
  }
}



/**
  The function searches the list of property numbers for the given number prop.
  It returns whether the aprop is a new property number.

  @param prop - searched property number
  @param propid - array of with list of property numbers
  @param n      - length of array propid
  
  @retval 0 - the property was found in the list
  @retval 1 - the property was not found in the list

  Created by Tomas Koudelka, 08.2011
*/
long sfaces::search_list_newprop(long prop, long *propid, long n)
{
  long i;

  for(i=0; i<n; i++)
  {
    if (propid[i] == prop)
      return 0;
  }
  return 1;
}



/**
  The function generates from the list of additional surfaces new lists which contain 
  the list of different property numbers, lists of nodes for particular property numbers, 
  lists of elements for particular property numbers and lists of surface id for particular
  lists of elements with given property number.

  @param nn - total number of nodes in the mesh
  @param ne - total number of elements in the mesh
  @param elems - array of all elements in the mesh
  @param nadjelnod - array of number of adjacent elements to nodes
  @param adjelnod - array of pointers to arrays of adjacent elements ids to nodes

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sfaces::gen_list_surf(long nn, long ne, selement *elems, long *nadjelnod, long **adjelnod)
{
  long i;
  long aux;
  long *propid;


  if (nfcs > 0)
  {
    propid = new long[nfcs];
    nprop = 0;
    propid[nprop] = aux = faces[0].prop;
    nprop++;
    for(i=0; i<nfcs; i++)
    {
      if (faces[i].prop != aux)
      { 
        aux = faces[i].prop; 
        if (search_list_newprop(aux, propid, nprop))
        {
          propid[nprop] = aux;
          nprop++;
        }
      }
    }    
    proplist = new long[nprop];
    for(i=0; i<nprop; i++)
      proplist[i] = propid[i];
    delete [] propid;
    prop_list_nod           = new long*[nprop];
    prop_list_nod_num       = new long[nprop];
    prop_list_elem          = new long*[nprop];
    prop_list_elem_num      = new long[nprop];
    prop_list_elem_sfid     = new long**[nprop];
    prop_list_elem_sfid_num = new long*[nprop];
    memset(prop_list_nod, 0, sizeof(*prop_list_nod)*nprop);
    memset(prop_list_nod_num, 0, sizeof(*prop_list_nod_num)*nprop);
    memset(prop_list_elem, 0, sizeof(*prop_list_elem)*nprop);
    memset(prop_list_elem_num, 0, sizeof(*prop_list_elem_num)*nprop);
    memset(prop_list_elem_sfid, 0, sizeof(*prop_list_elem_sfid)*nprop);
    memset(prop_list_elem_sfid_num, 0, sizeof(*prop_list_elem_sfid_num)*nprop);

    gen_list_node(nn);
    gen_list_elem(ne, elems, nadjelnod, adjelnod);
  }
}



/**
  The function generates from the list of additional surfaces new list which contain 
  lists of nodes for particular different property numbers. The lists are stored
  in the array prop_list_nod where prop_list_nod[i] represents array of
  node numbers for the property number proplist[i]. The length of array
  prop_list_nod[i] is given by prop_list_nod_num[i].

  @param nn - total number of nodes in the mesh

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sfaces::gen_list_node(long nn)
{
  long i, j, k, nnodp;
  long *nodid;

  nodid = new long[nn];

  for(i=0; i<nprop; i++)
  {
    memset(nodid, 0, sizeof(*nodid)*nn);   
    nnodp = 0;
    for(j=0; j<nfcs; j++)
    {
      if (faces[j].prop == proplist[i])
      {
        for(k=0; k<faces[j].nnod; k++)
        {
          if (nodid[faces[j].nodes[k]] == 0)
          {
            nodid[faces[j].nodes[k]] = 1;
            nnodp++;
          }
        }
      }
    }
    prop_list_nod_num[i] = nnodp;
    prop_list_nod[i] = new long[nnodp];
    k = 0;
    for(j=0; j<nn; j++)
    {
      if (nodid[j])
      {
        prop_list_nod[i][k] = j;
        k++;
      }
    }
  }
  delete [] nodid;
}



/**
  The function generates from the list of additional surfaces new lists which contain 
  lists of elements for particular different property numbers. The lists are stored
  in the array prop_list_elem where prop_list_elem[i] represents array of
  element numbers for the property number proplist[i]. The length of array
  prop_list_elem[i] is given by prop_list_elem_num[i]. 

  Additionally, the function  creates lists of surface indices with given property number 
  for each list of elements. Surface indices are stored in the 
  array prop_list_elem_sfid where the prop_list_elem_sfid[i][j] represents array of 
  surface indices with property number proplist[i] for the j-th element from the prop_list_elem[i].
  The length of array  prop_list_elem_sfid[i] is given by prop_list_elem_sfid_num[i][j].    

  @param ne - total number of elements in the mesh
  @param elems - array of all elements in the mesh
  @param nadjelnod - array of number of adjacent elements to nodes
  @param adjelnod - array of pointers to arrays of adjacent elements ids to nodes

  @return The function does not return anything.

  Created by Tomas Koudelka, 08.2011
*/
void sfaces::gen_list_elem(long ne, selement *elems, long *nadjelnod, long **adjelnod)
{
  long i, j, k, l, m, kmin, nelemp, nadjelmin;
  long aeid_min, aeid;
  long *elemid;
  long **elemsfcid;
  long sfid, aux;

  elemid    = new long[ne];
  elemsfcid = new long*[ne];
  for(i=0; i<ne; i++)
    elemsfcid[i] = new long [elems[i].nsurf];

  for(i=0; i<nprop; i++) // property index
  {
    memset(elemid, 0, sizeof(*elemid)*ne);
    for(j=0; j<ne; j++)
      memset(elemsfcid[j], 0, sizeof(*elemsfcid[j])*elems[j].nsurf);
    nelemp = 0;

    for(j=0; j<nfcs; j++) // face index
    {
      if (faces[j].prop == proplist[i])
      {        
        // search minimum number of adjacent elements from nodes of j-th face
        nadjelmin = ne;
        kmin = 0;
        for (k=0; k<faces[j].nnod; k++) // index of nodes on j-th face
        {
          if (nadjelnod[faces[j].nodes[k]] < nadjelmin)
          {
            kmin = k; // store index of node with minimum number of adjacent elements
            nadjelmin = nadjelnod[faces[j].nodes[k]];
          }
        }
     
        // compare adjacent element ids of kmin-th face node with rest adjacent elements of 
        // remaining face nodes
        for(k=0; k<nadjelmin; k++) // index of adjacent elements to m-th node of j-th face
        {
          aeid_min = adjelnod[faces[j].nodes[kmin]][k];
          aux=1; // number of occurrences of the element in adjacent elements of other nodes
          for(l=0; l<faces[j].nnod; l++)
          {
            if(kmin == l) // search all lists of adjacent elements except of the one with kmin-th index
              continue;

            for(m=0; m<nadjelnod[faces[j].nodes[l]]; m++)
            {
              aeid = adjelnod[faces[j].nodes[l]][m];
              if (aeid == aeid_min)
              {
                aux++;
                break;
              }
              if (aeid_min < aeid) // adjelnod arrays are sorted ascending -> remaining numbers are greater than aeid_min
              {
                aux = -1;
                break;
              }
            }
            if (aux < 0) // the element aeid_min was not found in one of the list -> continue with next element id
              break;
          }
          if (aux==faces[j].nnod)
          {
            sfid = elems[aeid_min].compare_surf(faces[j].nodes, faces[j].nnod);
            if (sfid >= 0)
            {
              if (elemid[aeid_min] == 0)
              {
                nelemp++;
                elemid[aeid_min] = 1;
              }
              elemsfcid[aeid_min][sfid] = 1;
            }
          }
        }
      }
    }
    prop_list_elem_num[i]      = nelemp;
    prop_list_elem[i]          = new long[nelemp];
    prop_list_elem_sfid_num[i] = new long[nelemp];
    prop_list_elem_sfid[i]     = new long*[nelemp];
    k = 0;
    for(j=0; j<ne; j++)
    {
      if (elemid[j] == 0)
        continue;

      prop_list_elem[i][k] = j;
      aux = 0;
      for(l = 0; l<elems[j].nsurf; l++)
      {
        if (elemsfcid[j][l])
          aux++;
      }
      prop_list_elem_sfid_num[i][k] = aux;
      if (aux)
        prop_list_elem_sfid[i][k] = new long[aux];
      else
        prop_list_elem_sfid[i][k] = NULL;
      aux = 0;
      for(l = 0; l<elems[j].nsurf; l++)
      {
        if (elemsfcid[j][l])
        {
          prop_list_elem_sfid[i][k][aux] = l;
          aux++;
        }
      }
      k++;
    }
  }
  delete [] elemid;
  for(i=0; i<ne; i++)
    delete [] elemsfcid[i];
  delete [] elemsfcid;
}


/**
  This constructor initializes data mebers with NULL or zero values

  Created by Tomas Koudelka 03/2008, tomas.koudelka@fsv.cvut.cz
*/
siftop :: siftop (void)
{
  elements  = NULL;
  nodes     = NULL;
  gnn       = NULL;
  edges     = NULL;
  surfaces  = NULL;
  nadjelnod = NULL;
  adjelnod  = NULL;
  nn        = 0;
  ne        = 0;
  memset(npet, 0, sizeof(npet));
  nnsd = NULL;
}



/**
  This constructor initializes arrays of nodes and elements to be a mesh with the given number inn of nodes
  and given number ine of elements with given type te. Default initialization of the nodal objects and 
  element objects is being performed, i.e. there is no definition of element connectivity as well as 
  zero nodal coordinates.

  @param inn - required number of nodes in the new mesh
  @param ine - required number of elements in the new mesh
  @param te  - required element type

  Created by Tomas Koudelka 03/2008, tomas.koudelka@fsv.cvut.cz
*/
siftop :: siftop (long inn, long ine, gtypel te)
{
  nn        = inn;
  ne        = ine;
  nodes     = new snode[nn];
  elements  = new selement[ne];
  for (long i=0; i<ne; i++){
    elements[i].type = te;
    elements[i].alloc(1);
  }
    
  gnn       = NULL;
  edges     = NULL;
  surfaces  = NULL;
  nadjelnod = NULL;
  adjelnod  = NULL;
  memset(npet, 0, sizeof(npet));
  meshtype = 0;
  nsd = 0;
  nnsd = NULL;
}



/**
  This copy constructor initializes all data members from the argument src.
  It is defined with explicit specifier so it cannot be used for arguments passed by values 
  in the function calls.

  Created by Tomas Koudelka 01/2024, tomas.koudelka@fsv.cvut.cz
*/
siftop :: siftop (const siftop &src)
{
  nn = src.nn;
  ne = src.ne;

  if (src.nodes){
    nodes = new snode[nn];
    if (src.gnn){
      gnn = new long[nn];
      memcpy(gnn, src.gnn, sizeof(*gnn)*nn);
    }
    for (long i=0; i<nn; i++)
      src.nodes[i].copyto(nodes[i]);
  }
  else
    nodes = NULL;

  if (src.elements){
    elements = new selement[ne];
    for (long i=0; i<ne; i++)
      src.elements[i].copyto(elements[i]);
  }
  else
    elements = NULL;

  if (src.edges){
    edges = new sedges;
    src.edges->copyto(*edges);
  }
  else
    edges = NULL;

  if (src.surfaces){
    surfaces = new sfaces;
    src.surfaces->copyto(*surfaces);
  }
  else
    surfaces = NULL;

  if (src.nadjelnod){
    nadjelnod = new long[nn];
    memcpy(nadjelnod, src.nadjelnod, sizeof(*nadjelnod)*nn);
  }
  else
    nadjelnod = NULL;

  if (src.adjelnod){
    adjelnod = new long*[nn];
    for(long i=0; i<nn; i++){
      adjelnod[i] = new long[nadjelnod[i]];
      memcpy(adjelnod[i], src.adjelnod[i], sizeof(*adjelnod[i])*nadjelnod[i]);
    }
  }
  else
    adjelnod  = NULL;

  memcpy(npet, src.npet, sizeof(npet));

  meshtype = src.meshtype;
  nsd = src.nsd;
  if (nnsd){
    nnsd = new long[nsd];
    memcpy(nnsd, src.nnsd, sizeof(*nnsd)*nsd);
  }
  else
    nnsd = NULL;
}



/**
  This destructor deallocates memory used by data mebers

  Created by Tomas Koudelka 03.2008, tomas.koudelka@fsv.cvut.cz
*/
siftop :: ~siftop(void)
{
  long i;

  if (elements)
    delete [] elements;

  if (nodes)
    delete [] nodes;

  if (gnn)
    delete [] gnn;

  if (edges)
    delete edges;

  if (surfaces)
    delete surfaces;

  if (nadjelnod)
  {
    for (i=0; i<nn; i++) 
      delete [] adjelnod[i];
    delete [] nadjelnod;
    delete [] adjelnod;
  }

  delete [] nnsd;
}



/** 
  The method makes partial copy of the given siftop instance but allocates additional space for nodes, elements, 
  edge and surface property objects.
  
  @param[out] top - destination topology where to save the partial copy
  @param[in] tndn - total number of doubled nodes
  @param[in] tnifel - total number of the generated interface elements
  @param[in] num_dnedgprop - the number of newly generated edge property objects involving doubled nodes
  @param[in] num_dnsurfprop - the number of newly generated surface property objects involving doubled nodes

  @return The function does not return anything.

  Created by TKo, 01.2024
*/
void siftop::partial_copy(siftop &top, long tndn, long tnifel, long num_dnedgprop, long num_dnsurfprop) const
{
  top.nn = nn+tndn;
  top.ne = ne+tnifel;

  if (top.nodes)
    delete [] top.nodes;
  if (top.gnn)
    delete [] top.gnn;
  if (nodes){
    top.nodes = new snode[top.nn];
    if (gnn){
      top.gnn = new long[top.nn];
      memcpy(top.gnn, gnn, sizeof(*gnn)*nn);
    }
    for (long i=0; i<nn; i++)
      nodes[i].copyto(top.nodes[i]);
  }
  else
    top.nodes = NULL;

  if (top.elements)
    delete [] top.elements;
  if (elements){
    top.elements = new selement[top.ne];
    for (long i=0; i<ne; i++)
      elements[i].copyto(top.elements[i]);
  }
  else
    top.elements = NULL;

  if (top.edges)
    delete top.edges;
  if (edges){
    top.edges = new sedges;
    top.edges->nedg = edges->nedg+num_dnedgprop;
    top.edges->edges = new sedge[top.edges->nedg];
    for (long i=0; i<edges->nedg; i++)
      edges->edges[i].copyto(top.edges->edges[i]);    
  }
  else
    top.edges = NULL;

  if (top.surfaces)
    delete top.surfaces;
  if (surfaces){
    top.surfaces = new sfaces;
    top.surfaces->nfcs = surfaces->nfcs+num_dnsurfprop;
    top.surfaces->faces = new sface[top.surfaces->nfcs];
    for (long i=0; i<surfaces->nfcs; i++)
      surfaces->faces[i].copyto(top.surfaces->faces[i]);
  }
  else
    top.surfaces = NULL;

  memcpy(top.npet, npet, sizeof(npet));

  top.meshtype = meshtype;
  top.nsd = nsd;
  if (top.nnsd)
    delete [] top.nnsd;
  if (nnsd){
    top.nnsd = new long[nsd];
    memcpy(top.nnsd, nnsd, sizeof(*nnsd)*nsd);
  }
  else
    top.nnsd = NULL;
}



/**
  This method reads topology of given domain from file, which is in the JKTK format

  @param in   - pointer to structure with opened file of input data
  @param rgnn - indicator whether global node numbers array should be read (=1) or no (=0).
                This is support for paralell computing.
  @param rebp - indicator whether element edge/surface property identifiers are given in the input file.

  @retval 0 : on success
  @retval 1 : on error in reading of node
  @retval 2 : on error in reading of element
  @retval 3 : on error in reading of global node numbers

  @note JKTK mesh format is generated by programs gensif??? and details are in . . .

  created by Tomas Koudelka 3.2008, tomas.koudelka@fsv.cvut.cz
*/
long siftop::read(XFILE *in, long rgnn, long rebp)
{
  long i, j, eid, nid;
  long nfaces, nedges;
  
  //  number of nodes
  xfscanf(in, "%k%ld", "num_nodes", &nn);
  nodes = new snode[nn];
  for (i = 0; i < nn; i++)
  {
    xfscanf(in, "%k%ld", "node_id", &nid);
    if ((nid < 1) && (nid > nn))
    {
      print_err("node number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)\n",
                __FILE__, __LINE__, __func__, nid, nn, in->fname, in->line, in->col);
      return (1);
    }
    nid--;
    xfscanf(in, "%k%le %k%le %k%le %k%ld", "x", &nodes[nid].x, "y", &nodes[nid].y, "z", &nodes[nid].z, "numprop", &nodes[nid].nprop);
    /*    if (nodes[nid].nprop < 1)
    {
      print_err("node %ld must have at least one property (file %s, line=%ld, col=%ld)\n", __FILE__, __LINE__, __func__, nid, in->fname, in->line, in->col);
      return (1);
    }
    nodes[nid].alloc(nodes[nid].nprop);
    for (j=0; j < nodes[nid].nprop; j++)
    xfscanf(in, "%k%m%ld", "prop", &entityp_kwdset, nodes[nid].entid+j, nodes[nid].prop+j);*/
    if (nodes[nid].nprop > 0)
    {
      nodes[nid].alloc(nodes[nid].nprop);
      for (j=0; j < nodes[nid].nprop; j++)
        xfscanf(in, "%k%m%ld", "prop", &entityp_kwdset, nodes[nid].entid+j, nodes[nid].prop+j);
    }
  }


  //  number of elements
  xfscanf(in, "%k%ld", "num_elements", &ne);
  elements = new selement[ne];
  for (i = 0; i < ne; i++)
  {
    xfscanf(in, "%k%ld", "elem_id", &eid);
    eid--;
    if ((eid < 0) || (eid >= ne))
    {
      print_err("element number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)", __FILE__, __LINE__, __func__, eid+1, ne, in->fname, in->line, in->col);
      return (2);
    }
    xfscanf(in, "%k%m", "eltype", &gtypel_kwdset, &elements[eid].type);
    if (elements[eid].alloc(rebp))
    {
      print_err("unknown type on element %ld is required (file %s, line=%ld, col=%ld)",
                __FILE__, __LINE__, __func__, eid+1, in->fname, in->line, in->col);
      return (2);
    }
    xfscanf(in, "%k", "enodes");
    for (j = 0; j < elements[eid].nne; j++)
    {
      xfscanf(in, "%ld", &elements[eid].nodes[j]);
      elements[eid].nodes[j]--;
      if ((elements[eid].nodes[j] < 0) || (elements[eid].nodes[j] >= nn))
      {
        print_err("node number %ld on element %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)",
                  __FILE__, __LINE__, __func__, elements[eid].nodes[j]+1, i+1, nn, in->fname, in->line, in->col);
        return (2);
      }
    }
    xfscanf(in, "%k%ld", "eprop", &elements[eid].prop);
    if (rebp > 0)
    {
      if (elements[eid].ned)
        xfscanf(in, "%k", "propedg");
      for(j = 0; j < elements[eid].ned; j++)
        xfscanf(in, "%ld", elements[eid].propedg+j);
      if (rebp == 1) // due to compatibility with Mesheditor output which cannot print surface property for 2D elements
      {              // in the case of modified file by Mesheditor, readege flag can be set to 2 and zero surface properties will be used on all 2D elements
        if (elements[eid].nsurf)
          xfscanf(in, "%k", "propsurf");
        for(j = 0; j < elements[eid].nsurf; j++)
          xfscanf(in, "%ld", elements[eid].propsurf+j);
      }
    }
    else{ // skip rest of line
      skipline(in);
    }
  }
    
  
  if (rgnn==1){
    //  parallel code
    gnn = new long[nn];
    memset(gnn, 0, sizeof(*gnn)*nn);
    for (i = 0; i < nn; i++){
      xfscanf(in, "%k%ld", "node_id", &j);
      if ((j < 1) && (j > nn)){
	print_err("node number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)\n",
                  __FILE__, __LINE__, __func__, j, nn, in->fname, in->line, in->col);
	return (3);
      }
      j--;
      xfscanf(in, "%k%ld", "glob_id", gnn+j);
    }
  }

  if (rgnn==2){
    //  sequential code, BOSS method or sequential version of the FETI method

    //  begin, added, JK, 28.10.2009
    
    //  mesh description
    //  number of subdomains/aggregates
    xfscanf (in,"%ld %ld",&meshtype,&nsd);
    
    //  number of nodes on subdomains/aggregates
    nnsd = new long [nsd];
    //  loop over the number of subdomains/aggregates
    for (i=0;i<nsd;i++){
      xfscanf (in,"%ld",nnsd+i);
    }
    //  end, JK, 28.10.2009
    
    gnn = new long[nn];
    memset(gnn, 0, sizeof(*gnn)*nn);
    for (i = 0; i < nn; i++){
      xfscanf(in, "%k%ld", "node_id", &j);
      if ((j < 1) && (j > nn)){
	print_err("node number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)\n",
                  __FILE__, __LINE__, __func__, j, nn, in->fname, in->line, in->col);
	return (3);
      }
      j--;
      xfscanf(in, "%k%ld", "glob_id", gnn+j);
    }
  }

  kwd_handling bmode = in->kwdmode;
  in->kwdmode = sequent_mode;
  // read additional data about surfaces created in the MeshEditor
  if (xfscanf(in, " %+k", "faces"))
  {
    in->kwdmode = bmode;
    xfscanf(in, "%ld", &nfaces);  
    surfaces = new sfaces(nfaces);
    surfaces->read(in);
    if (nadjelnod == NULL)
      gen_adjelnod();
    surfaces->gen_list_surf(nn, ne, elements, nadjelnod, adjelnod);
  }
  
  in->kwdmode = sequent_mode;
  // read additional data about edges created in the MeshEditor
  if (xfscanf(in, " %+k", "edges"))
  {
    in->kwdmode = bmode;
    xfscanf(in, "%ld", &nedges);
    edges = new sedges(nedges);
    edges->read(in); 
    if (nadjelnod == NULL)
      gen_adjelnod();
    edges->gen_list_edg(nn, ne, elements, nadjelnod, adjelnod);
  }
  in->kwdmode = bmode;
    
  return(0);
}



/**
  This method copies mesh topology to general topology structure.
  It is used in preprocessors especially.

  @param top  - pointer to structure of general topology

  created by Tomas Koudelka 3.2008, tomas.koudelka@fsv.cvut.cz
*/
void siftop::exporttop (gtopology *top)
{
  long i, j;
  top->alloc_nodes(nn);
  top->alloc_elements(ne);
  for (i = 0; i < nn; i++)
  {
    top->gnodes[i].x = nodes[i].x;
    top->gnodes[i].y = nodes[i].y;
    top->gnodes[i].z = nodes[i].z;
  }
  for (i = 0; i < ne; i++)
  {
    top->gelements[i].nne = elements[i].nne;
    top->gelements[i].nodes = new long[elements[i].nne];
    for (j = 0; j < elements[i].nne; j++)
      top->gelements[i].nodes[j] = elements[i].nodes[j];
    switch(elements[i].type){
    case isolinear1d:{
      top->gelements[i].get=linbar;
      break;
    }
    case isoquadratic1d:{
      top->gelements[i].get=quadbar;
      break;
    }
    case trianglelinear:{
      top->gelements[i].get=lintriag;
      break;
    }
    case trianglequadratic:{
      top->gelements[i].get=quadtriag;
      break;
    }
    case isolinear2d:{
      top->gelements[i].get=linquad;
      break;
    }
    case isoquadratic2d:{
      top->gelements[i].get=quadquad;
      break;
    } 
    case isocubic2d:{
      top->gelements[i].get=cubicquad;
      break;
    } 
    case tetrahedronlinear:{
      top->gelements[i].get=lintetra;
      break;
    }
    case tetrahedronquadratic:{
      top->gelements[i].get=quadtetra;
      break;
    }
    case pyramidelinear:{
      top->gelements[i].get=noelem;
      break;
    }
    case pyramidequadratic:{
      top->gelements[i].get=noelem;
      break;
    }
    case wedgelinear:{
      top->gelements[i].get=noelem;
      break;
    }
    case wedgequadratic:{
      top->gelements[i].get=noelem;
      break;
    }
    case isolinear3d:{
      top->gelements[i].get=linhexa;
      break;
    }
    case isoquadratic3d:{
      top->gelements[i].get=quadhexa;
      break;
    }
    case all_element_types:{
      top->gelements[i].get=noelem;
      break;
    }
    default:{
      print_err("unknown type of element is required", __FILE__, __LINE__, __func__);
    }
    }
  }
  
}



/**
  This function imports topology in T3D format from the opened file in.

  @param in - pointer to structure with opened file of input data
  @param paral - flag for topology format used in parallel computations

  @retval 0 on success
  @retval 1 in the case of invalid node number
  @retval 2 in the case of invalid element number

  created  03.2008, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long siftop::import_t3d (XFILE *in, long paral)
{
  long i;
  long type, renum, mode, temp;
  long outmode[9];
  long ret;


  xfscanf(in, "%ld", &type);   // reads type of generator
  xfscanf(in, "%ld", &degree); // reads degrees of element aproximation
  xfscanf(in, "%ld", &renum);  // reads indicator of renumbering
  xfscanf(in, "%ld", &mode);   // reads output mode for elements and nodes
    
  if (paral==1)
    xfscanf(in,"%ld", &temp);
    
  // decodes mode variable and stores results in outmode array
  for (i = 0; i < 9; i++)
    outmode[i]=(mode>>i) & 1;

  xfscanf(in, "%ld", &nn);
  // allocation of memory for nodes
  if (nodes)
    delete [] nodes; 
  nodes = new snode[nn];
  if (paral)  // mesh for parallel computations is required 
  {
    if (gnn)
      delete [] gnn; 
    gnn = new long[nn];
    memset(gnn, 0, sizeof(*gnn)*nn);
  }
  switch (type) 
  // number of particular element types depends on mesh generator type 
  {
    case 3:
      // reading of number of edges, trias, tetras
      xfscanf(in, "%ld%ld%ld", npet+isolinear1d+degree-2, npet+trianglelinear+degree-2, npet+tetrahedronlinear+degree-2);
      break ;
    case 4:
      // reading of number of edges, quads and bricks
      xfscanf(in, "%ld%ld%ld", npet+isolinear1d+degree-2, npet+isolinear2d+degree-2, npet+isolinear3d+degree-2);
      break ;
    case 7:
      // reading of number of edges, trias, quads
      xfscanf(in, "%ld%ld%ld",  npet+isolinear1d+degree-2, npet+trianglelinear+degree-2, npet+isolinear2d+degree-2);
      // reading of number of tetras, pyrams, wedges and bricks
      xfscanf(in, "%ld%ld%ld%ld", npet+tetrahedronlinear+degree-2, npet+pyramidelinear+degree-2, npet+wedgelinear+degree-2, npet+isolinear3d+degree-2);
      break ;
    default:
      print_err("unknown type of generator is required", __FILE__, __LINE__, __func__);
  }
  // total number of elements
  ne=0;
  for (i = 0; i < all_element_types-1; i++)
    ne+=npet[i];

  ret = 0;
  // import of section of nodes
  ret = import_t3d_nodes(in, paral, outmode) ;
  if (ret)
    return 1;

  // allocation of memory for elements
  if (elements)
    delete [] elements;
  elements = new selement[ne];

  // import of particular element types 
  ret = import_t3d_elements(in, outmode);
  if (ret)
    return 2;

  return ret;
}



/**
  Function imports section of nodes in the T3D format.

  @param in - pointer to structure with opened file of input data
  @param paral   - indicates whether global node numbers are read (=1) or not (=0).
  @param outmode - array with description of output type (see T3D documentation).
                   Each aray element represents one bit flag from T3D output type.

  @retval 0 : on succes
  @retval 1 : in case that node number is out of range

  created  03.2008, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long siftop::import_t3d_nodes (XFILE *in, long paral, long *outmode)
{
  long i, j, nid;
  long nent;
  entityp ent;
  long eprop;
  char emsg[1001];

  for (i = 0; i < nn; i++)
  {
    xfscanf(in, "%ld", &nid);
    if ((nid < 1) || (nid > nn))
    {
      sprintf(emsg, "node number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)\n", nid, nn, in->fname, in->line, in->col);
      print_err(emsg, __FILE__, __LINE__, __func__);
      return (1);
    }
    nid--;
    if (paral)
      xfscanf(in, "%ld", gnn+nid);
    xfscanf(in, "%le %le %le %m %*ld %ld", &nodes[nid].x, &nodes[nid].y, &nodes[nid].z, &entityp_kwdset, &ent, &eprop);
    if (outmode[omt_sec_clasif_nod])
    {
      xfscanf(in, "%ld", &nent);
      nent++;
      nodes[nid].alloc(nent);
      nodes[nid].entid[0] = ent;
      nodes[nid].prop[0]  = eprop;
      for (j=1; j < nent; j++)
        xfscanf(in, "%m %*ld %ld", &entityp_kwdset, nodes[nid].entid+j, nodes[nid].prop+j);
    }
    else
    {
      nent = 1;
      nodes[nid].alloc(nent);
      nodes[nid].entid[0] = ent;
      nodes[nid].prop[0]  = eprop;
    }
    switch (nodes[nid].entid[0])
    {
      case evertex:
      case eregion:
        break ;
      case ecurve:
        if (outmode[omt_parameters] || outmode[omt_tangent])
          skipline(in); // skip rest of actual line
        break ;
      case esurface:
      case epatch:
      case eshell:
        if (outmode[omt_parameters] || outmode[omt_normal])
          skipline(in); // skip rest of actual line
        break ;
    }
  }
  return(0);
}



/**
  This function imports all elements from T3D file.

  @param in - pointer to structure with opened file of input data
  @param outmode - array with description of output type (see T3D documentation).
                   Each aray element represents one bit flag from T3D output type.
  @retval 0 : on success
  @retval 1 : in case that element number is out of range

  Created 03.2008, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long siftop::import_t3d_elements (XFILE *in, long *outmode)
{
  long i, j, eid;
  gtypel etype;
  char emsg[1001];
  
  if (degree == 1) 
    etype=isolinear1d;
  else
    etype=isoquadratic1d;

  for ( ; etype<all_element_types; etype++, etype++)
  {
    for (i=0; i<npet[etype-1]; i++)
    {
      xfscanf(in, "%ld", &eid);
      eid--;
      if ((eid < 0) || (eid >= ne))
      {
        sprintf(emsg, "element number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)", eid+1, nn, in->fname, in->line, in->col);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return(1);
      }

      elements[eid].type = etype;
      if (outmode[omt_neighbourhood])
        elements[eid].alloc(1);
      else
        elements[eid].alloc(1);

      for(j=0; j<elements[eid].nne; j++)
      {
        xfscanf(in, "%ld", &elements[eid].nodes[j]);
        elements[eid].nodes[j]--;
      }

      xfscanf(in, "%*ld %*ld %ld", &elements[eid].prop);
      switch (etype)
      {
        case isolinear1d:
        case isoquadratic1d:
          break;
        case trianglelinear:
        case trianglequadratic:
        case isolinear2d:
        case isoquadratic2d:
        case isocubic2d:
        { 
          elements[eid].propsurf[0] = elements[eid].prop;
          if (outmode[omt_neighbourhood])
            skipline(in); // skip rest of actual line
          if (outmode[omt_boundary_entities])
          {
            for (j=0; j<elements[eid].ned; j++)
              xfscanf(in, "%*ld");
            for (j=0; j<elements[eid].ned; j++)
              xfscanf(in, "%ld", elements[eid].propedg+j);
            skipline(in); // skip rest of actual line
          }
          else
            skipline(in); // skip rest of actual line
          break;
        }
        case tetrahedronlinear:
        case tetrahedronquadratic:
        case pyramidelinear:
        case pyramidequadratic:
        case wedgelinear:
        case wedgequadratic:
        case isolinear3d:
        case isoquadratic3d:
        {  
          if (outmode[omt_neighbourhood] || outmode[omt_associated_elements])
            skipline(in); // skip rest of actual line
          if (outmode[omt_boundary_entities])
          {
            for (j=0; j<elements[eid].nsurf; j++)
              xfscanf(in, "%*ld"); // skip entity id
            for (j=0; j<elements[eid].nsurf; j++)
              xfscanf(in, "%*ld"); // skip entity types
            for (j=0; j<elements[eid].nsurf; j++)
              xfscanf(in, "%ld", elements[eid].propsurf+j);
            skipline(in); // skip rest of actual line
          }
          else
            skipline(in); // skip rest of actual line
          break;
        }
        default:
          sprintf(emsg, "unknown type on element %ld is required (file %s, line=%ld, col=%ld)", eid+1, in->fname, in->line, in->col);
          print_err(emsg, __FILE__, __LINE__, __func__);
          break;
      }
    }
  }
  return(0);
}


// keywords for beginnings of GiD mesh file sections
enum  bsecgid {begsecgid_coord=0, begsecgid_elem=1};
const enumstr bsecgid_str[] = {{"Coordinates",0}, {"Elements",1}};
const kwdset bsecgid_kwdset(sizeof(bsecgid_str)/sizeof(*bsecgid_str), bsecgid_str);

// keywords for ends of GiD mesh file sections
enum  esecgid {endsecgid_coord=0, endsecgid_elem=1};
const enumstr esecgid_str[] = {{"end Coordinates",0}, {"end Elements",1}};
const kwdset esecgid_kwdset(sizeof(esecgid_str)/sizeof(*esecgid_str), esecgid_str);


/**
  This function imports topology in GiD format from the opened file in.

  @param in - pointer to structure with opened file of input data

  @retval -1 if topology pointer S is not set
  @retval -2 if file specified by ofname cannot be opened

  created  03.2008, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void siftop::import_gid (XFILE */*in*/)
{
  // import of section of nodes
//  import_gid_nodes(in, paral, outmode) ;
/*
  - detekce noveho zacatku a konce sekce coordinates
  - detekce poctu uzlu v aktualni sekci coordinates 
  - detekce minimalniho cisla uzlu ve vsech sekcich => idn1
  - vypocet celkoveho poctu uzlu
  if (nodes)
    delete [] nodes; 
  nodes = new snode[nn];
  - cteni jednotlivych uzlu ze sekci, od cisla uzlu odecitat idn1
*/

  // import of particular element types 
//  import_gid_elements(in, outmode);
/*
 - zalozit pole gidovskych nazvu vsech moznych typu elementu, ktere se do gidu tisknou,
   jejich poradi je porad stejne (dano v mechprint export_gid_mesh). 
 - detekovat pouzite sekce elementu v gid souboru (zacinaji MESH ...)
 - detekovat pocty prvku v jednotlivych sekcich => pocty elementu pro jednotlive typy
 - detekce minimalniho cisla prvku ve vsech sekcich => ide1
 - vypocet celkoveho poctu elementu
  // allocation of memory for elements
  if (elements)
    delete [] elements;
  elements = new selement[ne];
 - cteni jednotlivych elementu po typech jako v T3D, od cisla uzlu odecitat idn1, od cisla prvku odecitat ide1

 - prirazeni property z elementu jednotlivym uzlum
*/



  return;
}



/**
   Function prints mesh data in SSMF.
   
   created by Tomas Koudelka 3.2008, tomas.koudelka@fsv.cvut.cz
*/
void siftop::print (FILE *out) const
{
  long i, j;
  
  fprintf (out,"%ld\n",nn);

  for (i = 0; i < nn; i++)
  {
    fprintf (out,"%ld %15.10le %15.10le %15.10le %ld",i+1, nodes[i].x, nodes[i].y, nodes[i].z, nodes[i].nprop);
    for(j=0; j<nodes[i].nprop; j++)
    {
      fprintf(out, " %d %ld", nodes[i].entid[j], nodes[i].prop[j]);
    }
    fprintf(out, "\n");
  }
  
  fprintf (out,"\n%ld\n", ne);
  
  for (i = 0; i < ne; i++)
  {
    fprintf (out,"%ld %d",i+1, elements[i].type);
    for (j = 0; j < elements[i].nne; j++)
      fprintf (out," %ld",elements[i].nodes[j]+1);
    fprintf (out," %ld",elements[i].prop);
    if (elements[i].propedg)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].ned; j++)
        fprintf (out," %ld",elements[i].propedg[j]);
    }
    if (elements[i].propsurf)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].nsurf; j++)
        fprintf (out," %ld",elements[i].propsurf[j]);
    }
    fprintf(out, "\n");
  }    
  if (edges)
    edges->print(out, NULL);
  if (surfaces)
    surfaces->print(out, NULL);
}



/**
  The function prints nodes to the opened text file with the shift indices.

  @param out - pointer to the opened text file
  @param shift - shift for nodal indices

  @return The function does not return anything.
 
  9.11.2010 JB
*/
void siftop::shift_print_nodes (FILE *out,long shift)  const
{
  long i,j;
 for (i = 0; i < nn; i++){
   fprintf (out,"%ld %15.10le %15.10le %15.10le %ld",i+1+shift, nodes[i].x, nodes[i].y, nodes[i].z, nodes[i].nprop);
   for(j=0; j<nodes[i].nprop; j++)
     {
       fprintf(out, " %d %ld", nodes[i].entid[j], nodes[i].prop[j]);
     }
   fprintf(out, "\n");
 }
}


/**
  The function prints elements to the opened text file with the shift indices.

  @param out - pointer to the opened text file
  @param shift - shift for element indices

  @return The function does not return anything.
 
  9.11.2010 JB
*/
void siftop::shift_print_elements (FILE *out,long shiftnn,long shiftne) const
{
  long i,j;
  for (i = 0; i < ne; i++){
    fprintf (out,"%ld %d",i+1+shiftne, elements[i].type);
    for (j = 0; j < elements[i].nne; j++)
      fprintf (out," %ld",elements[i].nodes[j]+1+shiftnn);
    fprintf (out," %ld",elements[i].prop);
    if (elements[i].propedg)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].ned; j++)
        fprintf (out," %ld",elements[i].propedg[j]);
    }
    if (elements[i].propsurf)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].nsurf; j++)
        fprintf (out," %ld",elements[i].propsurf[j]);
    }
    fprintf(out, "\n");
  }    
}



/**
  The function asembles list of adjacent elements to nodes.

  @return The function does not return anything but it fills nadjelnod and adjelnod arrays.

  Created by JK+TKo, 2.7.2013
*/
void siftop::gen_adjelnod(void)
{
  long i, j, nne, *enodes;

  //  allocation of array of counts of adjacent elements to nodes
  if (nadjelnod==NULL)
    nadjelnod = new long [nn];

  memset (nadjelnod,0,nn*sizeof(*nadjelnod));
  
  //  number of contributions
  for (i=0; i<ne; i++)
  {
    nne = elements[i].nne;
    enodes = elements[i].nodes;

    for (j=0;j<nne;j++)
      nadjelnod[enodes[j]]++;
  }
  
  //  allocation of array of numbers of adjacent elements to nodes
  if (adjelnod==NULL)
  {
    adjelnod = new long*[nn];

    for (i=0;i<nn;i++)
      adjelnod[i] = new long[nadjelnod[i]];
  }

  memset (nadjelnod,0,nn*sizeof(*nadjelnod));

  
  
  //  filling of array
  for (i=0;i<ne;i++)
  {
    nne = elements[i].nne;
    enodes = elements[i].nodes;

    for (j=0;j<nne;j++)
      adjelnod[enodes[j]][nadjelnod[enodes[j]]++] = i;
  }

  return;
}



/**
  Method gets set of nodes with property prop of entity ent
  Node numbers of found nodes are stored to array setnodes
  setnodes[i] < 0 -> i-th node is NOT involved in the set 
  setnodes[i] = i -> i-th node is memeber of set with required property

  @param prop  - searched property id 
  @param ent   - type of searched entity, in case of gsurface, 
                 the esurface, epatch and eshell are assumed to be equal to gsurface
  @param setnodes - allocated array for node numbers with given property of given entity
                    dimension of setnodes array must be nn.
  
  @return The function returns number of nodes with checked property. 

  Created 06.2009 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
  Modified 08.2011 by Tomas Koudelka
*/
long siftop::get_propent_nodes(long prop, gentity ent, long *setnodes) const
{
  long i, j, ret=0;
  for (i=0; i<nn; i++)
  {
    setnodes[i] = -1;
    if (nodes[i].searchprop(prop, ent, i, NULL, NULL))
    {
      setnodes[i] = i;
      ret++;
    }
  }
  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_nod_num[i]; j++)
        {
          setnodes[edges->prop_list_nod[i][j]] = edges->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_nod_num[i]; j++)
        {
          setnodes[surfaces->prop_list_nod[i][j]] = surfaces->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }
  return ret;
}



/**
  Method gets set of nodes with property prop of entity ent
  Node numbers of found nodes are stored to array setnodes,
  setnodes[i] < 0 -> i-th node is NOT involved in the set 
  setnodes[i] = i -> i-th node is memeber of set with required property
  

  @param prop  - searched property id 
  @param ent   - type of searched entity, in case of gsurface, 
                 the esurface, epatch and eshell are assumed to be equal to gsurface
  @param setnodes - array for node numbers with given property of given entity,
                    dimension of setnodes vector must be nn.
  
  @return The function returns number of nodes with checked property. 

  Created 06.2020 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_propent_nodes(long prop, gentity ent, ivector &setnodes) const
{
  long i, j, ret=0;
  for (i=0; i<nn; i++)
  {
    setnodes[i] = -1;
    if (nodes[i].searchprop(prop, ent, i, NULL, NULL))
    {
      setnodes[i] = i;
      ret++;
    }
  }
  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_nod_num[i]; j++)
        {
          setnodes[edges->prop_list_nod[i][j]] = edges->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_nod_num[i]; j++)
        {
          setnodes[surfaces->prop_list_nod[i][j]] = surfaces->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }
  return ret;
}



/**
  Method gets set of nodes with property prop of entity ent.
  Node numbers of found nodes are stored in the vector setnodes,
  where setnodes[j] = node number of node involved in the set with 
  required property.

  @param[in] prop  - searched property id 
  @param[in] ent   - type of searched entity, in case of gsurface, 
                 the esurface, epatch and eshell are assumed to be equal to gsurface
  @param[out] setnodes - array for node numbers with given property of given entity,
                         vector is NOT required to be preallocated, the function 
                         allocates its dimension to nn on the beginning and before return,
                         the vector dimension is modififed to correspond the nummber of
                         nodes in the set with required property.
  
  @return The function returns number of nodes with checked property. 

  Created 06.2020 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_propent_nodes_compact(long prop, gentity ent, ivector &setnodes) const
{
  long i, j, ret=0;
  reallocv(nn, setnodes);
  
  for (i=0; i<nn; i++)
  {
    setnodes[i] = -1;
    if (nodes[i].searchprop(prop, ent, i, NULL, NULL))
    {
      setnodes[ret] = i;
      ret++;
    }
  }
  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_nod_num[i]; j++)
        {
          setnodes[ret] = edges->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_nod_num[i]; j++)
        {
          setnodes[ret] = surfaces->prop_list_nod[i][j];
          ret++;
        }
      }
    }
  }
  setnodes.n = ret;
  return ret;
}



/**
  Method gets set of elements with property prop of entity ent.
  Found elements are indicated in the vector setelems,
  setelems[i] < 0 -> i-th element is NOT involved in the set 
  setelems[i] = i -> i-th element is memeber of set with required property

  @param prop[in]  - searched property id 
  @param ent[in]   - type of searched entity, in case of gsurface, 
                     the esurface, epatch and eshell are assumed to be equal to gsurface
  @param setelems[out] - array for indicating element numbers with given property 
                         of the given entity, dimension of setelems vector must be ne.
  @param cumul[in] - indicator whether to accumulate found element numbers in the setelems
                     (default is 0 which means no accumulation). In the case of cumul=1,
                     setelems is NOT allocated and its components are supposed to be set
                     to -1 at the input. Only found element numbers are then stored
                     at the corresponding setelems array components
  
  @return The function returns number of elements with checked property. 

  Created 06.2021 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_propent_elems(long prop, gentity ent, ivector &setelems, long accum) const
{
  long i, j, ret=0, k;
  long *auxp, aux;

  if (accum == 0)
    reallocv(ne, setelems);
  
  for (i=0; i<ne; i++)
  {
    if (accum == 0)  setelems[i] = -1;
    if (elements[i].searchprop(prop, ent, nodes, i, NULL, NULL, auxp, aux))
    {
      setelems[i] = i;
      ret++;
    }
  }
  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_elem_num[i]; j++)
        {
          k = edges->prop_list_elem[i][j];
          setelems(k) = k;
          ret++;
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_elem_num[i]; j++)
        {
          k = surfaces->prop_list_elem[i][j];
          setelems[k] = k;
          ret++;
        }
      }
    }
  }
  return ret;
}



/**
  Method gets set of elements with property prop of entity ent.
  Numbers of found elements are stored in the vector setelems,
  where setelems[j] = element number involved in the set with 
  required property.

  @param[in] prop  - searched property id 
  @param[in] ent   - type of searched entity, in case of gsurface, 
                     the esurface, epatch and eshell are assumed to be equal to gsurface
  @param[out] setelems - array for element numbers with given property of given entity,
                         vector is NOT required to be preallocated, the function 
                         allocates its dimension to ne on the beginning and before return,
                         the vector dimension is modififed to correspond the number of
                         nodes in the set with required property.
  
  @return The function returns number of elements with checked property. 

  Created 06.2021 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_propent_elems_compact(long prop, gentity ent, ivector &setelems) const
{
  long i, j, ret=0;
  long *auxp, aux;

  reallocv(ne, setelems);
  
  for (i=0; i<ne; i++)
  {
    setelems[i] = -1;
    if (elements[i].searchprop(prop, ent, nodes, i, NULL, NULL, auxp, aux))
    {
      setelems[ret] = i;
      ret++;
    }
  }
  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_elem_num[i]; j++)
        {
          setelems[ret] = edges->prop_list_elem[i][j];
          ret++;
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_elem_num[i]; j++)
        {
          setelems[ret] = surfaces->prop_list_elem[i][j];
          ret++;
        }
      }
    }
  }
  setelems.n = ret;
  return ret;
}



/**
  Method gets ndof of nodes with property prop of entity ent
  All nodes with given property of entity are searched and checked that
  they have same ndof. Node numbers of found nodes are stored to array setnodes

  @param prop  - searched property id 
  @param ent   - type of searched entity, in case of gsurface, 
                 the esurface, epatch and eshell are assumed to be equal to gsurface
  @param ndofn - array of ndofs for all nodes
  @param setnodes - allocated array for node numbers with given property of given entity
  
  @retval >0 - ndof of nodes with given property of given entity
  @retval 0  - in case that the nodes with given property of given entity have the different ndof
  @retval <0 - in case that no nodes with the given property were found

  Created 04.2008 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_ndofn(long prop, gentity ent, long *ndofn, long *setnodes) const
{
  long i, j, id;
  long ndof = -1;  
  for (i=0; i<nn; i++)
  {
    setnodes[i] = -1;
    if (nodes[i].searchprop(prop, ent, i, NULL, NULL))
    {
      setnodes[i] = i;
      if (ndof == -1)
        ndof = ndofn[i];
      else
      {
        if (ndof != ndofn[i])
          return 0;
      }
    }
  }


  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_nod_num[i]; j++)
        {
          id = edges->prop_list_nod[i][j];
          setnodes[id] = id;
          if (ndof == -1)
            ndof = ndofn[id];
          else
          {
            if (ndof != ndofn[id])
              return 0;
          }
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_nod_num[i]; j++)
        {
          id = surfaces->prop_list_nod[i][j];
          setnodes[id] = id;
          if (ndof == -1)
            ndof = ndofn[id];
          else
          {
            if (ndof != ndofn[id])
              return 0;
          }
        }
      }
    }
  }
  return ndof;
}




/**
  Method gets ndof of nodes with property prop of entity ent
  All nodes with given property of entity are searched and checked that
  they have same ndof. Node numbers of found nodes are stored to array setnodes

  @param prop  - searched property id 
  @param ent   - type of searched entity, in case of gsurface, 
                 the esurface, epatch and eshell are assumed to be equal to gsurface
  @param ndofn - array of ndofs for all nodes
  @param setnodes - integer vector of node numbers with given property of given entity
  
  @retval >0 - ndof of nodes with given property of given entity
  @retval 0  - in case that the nodes with given property of given entity have the different ndof
  @retval <0 - in case that no nodes with the given property were found

  Created 04.2008 by Tomas Koudelka tomas.koudelka@fsv.cvut.cz                
*/
long siftop::get_ndofn(long prop, gentity ent, long *ndofn, ivector &setnodes) const
{
  long i, j, id;
  long ndof = -1;  
  for (i=0; i<nn; i++)
  {
    setnodes[i] = -1;
    if (nodes[i].searchprop(prop, ent, i, NULL, NULL))
    {
      setnodes[i] = i;
      if (ndof == -1)
        ndof = ndofn[i];
      else
      {
        if (ndof != ndofn[i])
          return 0;
      }
    }
  }


  if ((ent == gcurve) && (edges))
  {
    // search additional edges specified in the MeshEditor
    for(i=0; i<edges->nprop; i++)
    {
      if (edges->proplist[i] == prop)
      {
        for(j=0; j<edges->prop_list_nod_num[i]; j++)
        {
          id = edges->prop_list_nod[i][j];
          setnodes[id] = id;
          if (ndof == -1)
            ndof = ndofn[id];
          else
          {
            if (ndof != ndofn[id])
              return 0;
          }
        }
      }
    }
  }

  if ((ent == gsurface) && (surfaces))
  {
    // search additional surfacs specified in the MeshEditor
    for(i=0; i<surfaces->nprop; i++)
    {
      if (surfaces->proplist[i] == prop)
      {
        for(j=0; j<surfaces->prop_list_nod_num[i]; j++)
        {
          id = surfaces->prop_list_nod[i][j];
          setnodes[id] = id;
          if (ndof == -1)
            ndof = ndofn[id];
          else
          {
            if (ndof != ndofn[id])
              return 0;
          }
        }
      }
    }
  }
  return ndof;
}




/**
  The function exports sets of used elements to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void siftop::export_gid_mesh(FILE *out, long idn1, long ide1) const
{
  long i, print_header, print_coord = 1;

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isolinear1d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-beams2\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isoquadratic1d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-beams3\" dimension 3  Elemtype Linear Nnode 3\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == trianglelinear)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-trias3\" dimension 3  Elemtype Triangle Nnode 3\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == trianglequadratic)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-trias6\" dimension 3  Elemtype Triangle Nnode 6\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isolinear2d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Quads4\" dimension 3  Elemtype Quadrilateral Nnode 4\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isoquadratic2d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Quads8\" dimension 3  Elemtype Quadrilateral Nnode 8\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == tetrahedronlinear)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Tetras4\" dimension 3  Elemtype Tetrahedra Nnode 4\n", ide1);
        print_header = 0;
        if (print_coord)
        {
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == tetrahedronquadratic)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Tetras10\" dimension 3  Elemtype Tetrahedra Nnode 10\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isolinear3d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Brick8\" dimension 3  Elemtype Hexahedra Nnode 8\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < ne; i++)
  {
    if (elements[i].type == isoquadratic3d)
    {
      if (print_header)
      {
        fprintf(out, "MESH \"%ld-Brick20\" dimension 3  Elemtype Hexahedra Nnode 20\n", ide1);
        print_header = 0;
        if (print_coord)
	{
          write_gid_nodes(out, idn1);
          print_coord = 0;
        }
        fprintf(out, "Elements\n");
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
}



/**
  The function exports nodes to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)

  @return The function does not return anything.

  Created 2010 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void siftop::write_gid_nodes(FILE *out, long idn1) const
{
  long i;
  fprintf(out, "Coordinates\n");
  for (i=0; i<nn; i++)
    fprintf(out, "%ld %e %e %e\n", i+idn1, nodes[i].x, nodes[i].y, nodes[i].z);
  fprintf(out, "end Coordinates\n");
}



/**
  The function exports element to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param i    - element id
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2010 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void siftop::write_gid_element(FILE *out, long i, long idn1, long ide1) const
{
  long j;
  fprintf(out, "%ld ", i+ide1);
  for (j=0; j<elements[i].nne; j++)
    fprintf(out, "%ld ", elements[i].nodes[j]+idn1);
  fprintf(out, "%ld\n", elements[i].prop);
}



/**
  The function reorders mesh from T3d into sifel format.

  @return The function does not return anything.

  Created by JB
*/
void siftop::t3d2sifelformat()
{
  long i,j;
  long *help;
  help = new long[20];
  
  for(i = 0; i < ne; i++){
    switch(elements[i].type){
    case tetrahedronlinear:{
      for(j = 0; j < 4; j++){
	help[j] = elements[i].nodes[j];
      }
      elements[i].nodes[2] = help[3];
      elements[i].nodes[3] = help[2];
	      
      for(j = 0; j < 4; j++){
	help[j] = elements[i].propsurf[j];
      }
      elements[i].propsurf[0]= help[2];
      elements[i].propsurf[1]= help[3];
      elements[i].propsurf[2]= help[0];
      elements[i].propsurf[3]= help[1];
      
      break;
    }
    case tetrahedronquadratic:{
      
      for(j = 0; j < 10; j++){
	help[j] = elements[i].nodes[j];
      }
      
      elements[i].nodes[2] = help[3];
      elements[i].nodes[3] = help[2];
      elements[i].nodes[5] = help[8];
      elements[i].nodes[6] = help[7];
      elements[i].nodes[7] = help[6];
      elements[i].nodes[8] = help[5];
      
      for(j = 0; j < 4; j++){
	help[j] = elements[i].propsurf[j];
      }
      elements[i].propsurf[0]= help[2];
      elements[i].propsurf[1]= help[3];
      elements[i].propsurf[2]= help[0];
      elements[i].propsurf[3]= help[1];
      
      break;
    }
    case isolinear3d:{
      for(j = 0; j < 8; j++){
	help[j] = elements[i].nodes[j];
      }
      
      elements[i].nodes[0] = help[3];
      elements[i].nodes[1] = help[0];
      elements[i].nodes[2] = help[4];
      elements[i].nodes[3] = help[7];
      elements[i].nodes[4] = help[2];
      elements[i].nodes[5] = help[1];
      elements[i].nodes[6] = help[5];
      elements[i].nodes[7] = help[6];
      
      
      for(j = 0; j < 6; j++){
	help[j] = elements[i].propsurf[j];
      }

      elements[i].propsurf[0] = help[0];
      elements[i].propsurf[1] = help[3];
      elements[i].propsurf[2] = help[1];
      elements[i].propsurf[3] = help[5];
      elements[i].propsurf[4] = help[2];
      elements[i].propsurf[5] = help[4];
      
      
      break;
    }
    case isoquadratic3d:{
      for(j = 0; j < 20; j++){
	help[j] = elements[i].nodes[j];
      }
      
      elements[i].nodes[0] = help[3];
      elements[i].nodes[1] = help[0];
      elements[i].nodes[2] = help[4];
      elements[i].nodes[3] = help[7];
      elements[i].nodes[4] = help[2];
      elements[i].nodes[5] = help[1];
      elements[i].nodes[6] = help[5];
      elements[i].nodes[7] = help[6];
      elements[i].nodes[8] = help[11];
      elements[i].nodes[9] = help[16];
      elements[i].nodes[10] = help[15];
      elements[i].nodes[11] = help[19];
      elements[i].nodes[12] = help[10];
      elements[i].nodes[13] = help[8];
      elements[i].nodes[14] = help[12];
      elements[i].nodes[15] = help[14];
      elements[i].nodes[16] = help[9];
      elements[i].nodes[17] = help[17];
      elements[i].nodes[18] = help[13];
      elements[i].nodes[19] = help[18];
      
      for(j = 0; j < 6; j++){
	help[j] = elements[i].propsurf[j];
      }
      elements[i].propsurf[0] = help[0];
      elements[i].propsurf[1] = help[3];
      elements[i].propsurf[2] = help[1];
      elements[i].propsurf[3] = help[5];
      elements[i].propsurf[4] = help[2];
      elements[i].propsurf[5] = help[4];
      
      break;
    }
    default:{
      break;
    }
    }
  }
  delete []help;
  
}




/**
  Function returns 3D node coordinates of the given element.
   
  @param[out] x,y,z - vectors containing node coordinates
  @param[in] eid - element id
   
  @return The function returns 3D(x,y,z) coordinates in parameters x, y and z. 

  Created by TKo, 7.2020
*/
void siftop::give_node_coord3d(vector &x, vector &y, vector &z, long eid) const
{
  long i, nne = elements[eid].nne;
  long *en = elements[eid].nodes;
  
  for (i=0; i<nne; i++){
    x(i) = nodes[en[i]].x;
    y(i) = nodes[en[i]].y;
    z(i) = nodes[en[i]].z;
  }
}



/**
  Function returns 2D node coordinates of the given element.
   
  @param[out] x,y - vectors containing node coordinates
  @param[in] eid - element id
   
  @return The function returns 2D(x,y) coordinates in parameters x and y. 

  Created by TKo, 7.2020
*/
void siftop::give_node_coord2d(vector &x, vector &y, long eid) const
{
  long i, nne = elements[eid].nne;
  long *en = elements[eid].nodes;

  for (i=0; i<nne; i++){
    x(i) = nodes[en[i]].x;
    y(i) = nodes[en[i]].y;
  }
}



/**
  Function returns dimension of the given element.
   
  @param eid - number of element
   
  @retval 1 - for 1D
  @retval 2 - for 2D
  @retval 3 - for 3D

  Created by TKo, 7.2020
*/
long siftop::give_dimension(long eid) const
{
  long ret = 0;
  
  switch(elements[eid].type){
    case isolinear1d:
    case isoquadratic1d:
      ret = 1;
      break;
    case trianglelinear:
    case trianglequadratic:
    case isolinear2d:
    case isoquadratic2d:
    case isocubic2d:
      ret = 2;
      break;
    case tetrahedronlinear:
    case tetrahedronquadratic:
    case pyramidelinear:
    case pyramidequadratic:
    case wedgelinear:
    case wedgequadratic:
    case isolinear3d:
    case isoquadratic3d:
      ret = 3;
      break;
    default:
      print_err("unknown general type of element %d is being required\n",
                __FILE__, __LINE__, __func__, int(elements[eid].type));
      abort();
  }
  return ret;
}



/**
  The function returns element nodes of the given element.

  @param eid[in] - index of element whose nodal numbers are required.
  @param enod[out] - resulting array of element nodes, array is reallocated
                     to the dimension equaled to number of nodes of the given element.

  @return The function does not reutn anything but fills argument array enod with required 
          node numbers.

  Created by Tomas Koudelka, 05.2023
*/
void siftop::give_elemnods(long eid, ivector &enod) const
{
  long nne = elements[eid].nne;
  reallocv(nne, enod);
  memcpy(enod.a, elements[eid].nodes, sizeof(*enod.a)*nne);
}



/** 
  The function searches on the given element eid for the closest node to the point of 
  given coordinates [px,py,pz]. It retruns the closest node id and the given distance in dmin.

  @param[in] eid - element id whose integration points are investigated
  @param[in] px, py, pz - global coordinates of the given point
  @param[out] dmin - attained minimum distance

  Created by Tomas Koudelka, 07.2020
*/  
long siftop::give_closest_node_coord(long eid, double px, double py, double pz, double &dmin) const
{
  long i, nid, nne = elements[eid].nne;
  long *enodes = elements[eid].nodes;
  dmin = DBL_MAX;
  double d;

  for(i=0; i<nne; i++)
  {
    d = sqr(px - nodes[enodes[i]].x) +
        sqr(py - nodes[enodes[i]].y) +
        sqr(pz - nodes[enodes[i]].z);
    if (d < dmin)
    {
      dmin = d;
      nid = enodes[i];
    }
  }
  // compute the minimum distnace 
  dmin = sqrt(dmin);

  return nid; 
}



/**
  The function computes centre of gravity for the given element.
  The coordinates are store in the argumnet coord.

  @param[in] eid - id of required element
  @param[out] coord - %vector structure for the resulting coordiates,
                      it must be preallocated to the dimension 3.

  @retval 0 - on success
  @retval 1 - in the case of error
*/
long siftop::centroid(long eid, vector &coord) const
{
  long nne = elements[eid].nne;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), bf(ASTCKVEC(nne));
  vector areacoord(ASTCKVEC(3)), volcoord(ASTCKVEC(4));

  give_node_coord3d(x, y, z, eid);

  switch(elements[eid].type){
    case isolinear1d:
      bf_lin_1d(bf.a, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case isoquadratic1d:
      bf_quad_1d(bf.a, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case trianglelinear:
      areacoord[0] = 0.333333333333333;  areacoord[1] = 0.333333333333333;
      areacoord[2] = 1.0 - areacoord[0] - areacoord[1];
      scprd(areacoord, x, coord(0));
      scprd(areacoord, y, coord(1));
      scprd(areacoord, z, coord(2));
      break;
    case trianglequadratic:
      bf_quad_3_2d(bf.a, 0.333333333333333, 0.333333333333333);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case isolinear2d:
      bf_lin_4_2d(bf.a, 0.0, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case isoquadratic2d:
      bf_quad_4_2d(bf.a, 0.0, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case isocubic2d:
      bf_cubic_4_2d(bf.a, 0.0, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));
      break;
    case tetrahedronlinear:
      volcoord[0] = 0.25;  volcoord[1] = 0.25;  volcoord[2] = 0.25;
      volcoord[3] = 1.0 - volcoord[0] - volcoord[1] - volcoord[2];
      scprd(volcoord, x, coord(0));
      scprd(volcoord, y, coord(1));
      scprd(volcoord, z, coord(2));
      break;
    case tetrahedronquadratic:
      bf_quad_tet(bf.a, 0.25, 0.25, 0.25);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));      
      break;
    case isolinear3d:
      bf_lin_hex_3d(bf.a, 0.0, 0.0, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));      
      break;
    case isoquadratic3d:
      bf_quad_hex_3d (bf.a, 0.0, 0.0, 0.0);
      scprd(bf, x, coord(0));
      scprd(bf, y, coord(1));
      scprd(bf, z, coord(2));      
      break;
    default:
      print_err("unknown general element type %d is being required on element %ld.\n",
                __FILE__, __LINE__, __func__, int(elements[eid].type), eid+1);
      abort();
  }
  return 0;
}



/**
  The function returns the square of maximum distance between given point pt and
  nodes of element eid.

  @param eid - element id
  @param pt  - %vector of point coordinates in 3D (x,y,z)

  @return The function returns the square of maximum distance from nodes.

  Created by TKo, 07.2020
*/
double siftop::max_sqrdist_nod_pt(long eid, vector &pt) const
{
  long i;
  long nne = elements[eid].nne;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  double d, maxd=0.0;

  give_node_coord3d(x, y, z, eid);
  
  for(i=0; i<nne; i++)
  {
    d  = (x[i] - pt[0])*(x[i] - pt[0]);
    d += (y[i] - pt[1])*(y[i] - pt[1]);
    d += (z[i] - pt[2])*(z[i] - pt[2]);
    if (d > maxd)
      maxd = d;
  }
  return maxd;
}



/**
  The function finds element surface id with required xi natural coordinate and returns array of 
  element nodes with the given xi natrual coordinate.

  @param eid[in] - element id
  @param xi[in] - rquired xi natural coordinate of nodes
  @param nod[out] - output array of element nodes with required xi coordinate,
                    array is reallocated -> if the dimension is less then
                    required, then DYNAMIC memory must be allocated in order 
                    to prevent stack deallocation after return from the function.

  @return The function returns surface id and list of nodes in the array nod.

  Created by Tomas Koudelka, 05.2023
*/
long siftop::give_enod_xicoord(long eid, double xi, ivector &nod) const
{
  long i, surfid = -1;
  gtypel et = elements[eid].type;
  switch (et){
    case isolinear3d:
      if (xi == 1.0)  surfid = 0;
      if (xi == -1.0) surfid = 2;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isolin3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isolin3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      else{
      }
      break;
    case isoquadratic3d:
      if (xi == 1.0)  surfid = 0;
      if (xi == -1.0) surfid = 2;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isoquad3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isoquad3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      break;
    default:
      print_err("uknown type of element (%d) is required on element %ld.\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }

  if (surfid < 0){
    print_err("absolute value of required xi coordinate |%le| != 1.0 on element %ld.\n",
              __FILE__, __LINE__, __func__, xi, eid+1);
    abort();
  }
  
  return surfid;
}



/**
  The function finds element surface id with required eta natural coordinate and returns array of 
  element nodes with the given eta natrual coordinate.

  @param eid[in] - element id
  @param eta[in] - rquired eta natural coordinate of nodes
  @param nod[out] - output array of element nodes with required xi coordinate,
                    array is reallocated -> if the dimension is less then
                    required, then DYNAMIC memory must be allocated in order 
                    to prevent stack deallocation after return from the function.

  @return The function returns surface id and list of nodes in the array nod.

  Created by Tomas Koudelka, 05.2023
*/
long siftop::give_enod_etacoord(long eid, double eta, ivector &nod) const
{
  long i, surfid=-1;
  gtypel et = elements[eid].type;
  switch (et){
    case isolinear3d:
      if (eta == 1.0)  surfid = 1;
      if (eta == -1.0) surfid = 3;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isolin3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isolin3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      break;
    case isoquadratic3d:
      if (eta == 1.0)  surfid = 1;
      if (eta == -1.0) surfid = 3;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isoquad3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isoquad3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      break;
    default:
      print_err("uknown type of element (%d) is required on element %ld.\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }

  if (surfid < 0){
    print_err("absolute value of required eta coordinate |%le| != 1.0 on element %ld.\n",
              __FILE__, __LINE__, __func__, eta, eid+1);
    abort();
  }
  
  return surfid;
}



/**
  The function finds element surface id with required zeta natural coordinate and returns array of 
  element nodes with the given zeta natrual coordinate.

  @param eid[in] - element id
  @param zeta[in] - rquired zeta natural coordinate of nodes
  @param nod[out] - output array of element nodes with required xi coordinate,
                    array is reallocated -> if the dimension is less then
                    required, then DYNAMIC memory must be allocated in order 
                    to prevent stack deallocation after return from the function.

  @return The function returns surface id and list of nodes in the array nod.

  Created by Tomas Koudelka, 05.2023
*/
long siftop::give_enod_zetacoord(long eid, double zeta, ivector &nod) const
{
  long i, surfid=-1;
  gtypel et = elements[eid].type;
  switch (et){
    case isolinear3d:
      if (zeta == 1.0)  surfid = 4;
      if (zeta == -1.0) surfid = 5;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isolin3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isolin3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      break;
    case isoquadratic3d:
      if (zeta == 1.0)  surfid = 4;
      if (zeta == -1.0) surfid = 5;
      if (surfid >= 0){
        reallocv(selement::nnsurf_isoquad3d[surfid], nod);
        for (i=0; i<nod.n; i++){
          nod[i] = selement::surfnod_isoquad3d[surfid][i];
          nod[i] = elements[eid].nodes[nod[i]];
        }
      }
      break;
    default:
      print_err("uknown type of element (%d) is required on element %ld.\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }
  if (surfid < 0){
    print_err("absolute value of required zeta coordinate |%le| != 1.0 on element %ld.\n",
              __FILE__, __LINE__, __func__, zeta, eid+1);
    abort();
  }

  return surfid;
}



/**
  Transforms natural coordinates of point on the given element to the local natrual coordinates of 
  required element surface. The local natural coordinate system of the surface is defined by the
  order of surface nodes defined in static arrays selement::surfnod_xxx where xxx represents the 
  given element type.
  
  @param eid[in] - required element id,
  @param sid[in] - id of required element surface,
  @param enc[in] - natural coordinates of the point in the element natural coordinate system,
  @param lnc[out] - natural coordinates of the point in local natural coordinate system of the required surface.

  @return The function returns local natural coordinates of the point in argument lnc. 

  Created by Tomas Koudelka, 05.2023
*/
void siftop::transform_natcoord_surf(long eid, long sid, const vector &enc, vector &lnc) const
{
  gtypel et = elements[eid].type;
  switch (et){
    case isolinear3d:
    case isoquadratic3d:
      if (sid == 0) {lnc(0) =  enc(1); lnc(1) = enc(2); lnc(2) = 0.0; return;}
      if (sid == 1) {lnc(0) = -enc(0); lnc(1) = enc(2); lnc(2) = 0.0; return;}
      if (sid == 2) {lnc(0) = -enc(1); lnc(1) = enc(2); lnc(2) = 0.0; return;}
      if (sid == 3) {lnc(0) =  enc(0); lnc(1) = enc(2); lnc(2) = 0.0; return;}
      if (sid == 4) {lnc(0) =  enc(0); lnc(1) = enc(1); lnc(2) = 0.0; return;}
      if (sid == 5) {lnc(0) =  enc(0); lnc(1) = enc(1); lnc(2) = 0.0; return;}
      print_err(" invalid id of required surface %ld on element %ld (type=%d).\n",
                __FILE__, __LINE__, __func__, sid, eid+1, et);
      abort();
    default:
      print_err("uknown type of element (%d) is required on element %ld.\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }
}



/**
  Transforms natural coordinates of point on the given element to the local natrual coordinates of 
  required element edge. The local natural coordinate system of the edge is defined by the
  order of edge nodes defined in static arrays selement::edgenod_xxx where xxx represents the 
  given element type.
  
  @param eid[in]  - required element id,
  @param edid[in] - id of required element edge,
  @param enc[in]  - natural coordinates of the point in the element natural coordinate system,
  @param lnc[out] - natural coordinates of the point in local natural coordinate system of the required edge.

  @return The function returns local natural coordinates of the point in argument lnc. 

  Created by Tomas Koudelka, 05.2023
*/
void siftop::transform_natcoord_edg(long eid, long edid, const vector &enc, vector &lnc) const
{
  gtypel et = elements[eid].type;
  switch (et){
    case isolinear3d:
    case isoquadratic3d:
      if (edid ==  0) {lnc(0) = -enc(0); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid ==  1) {lnc(0) = -enc(1); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid ==  2) {lnc(0) =  enc(0); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid ==  3) {lnc(0) =  enc(1); lnc(1) = 0.0; lnc(2) = 0.0; return;}

      if (edid ==  4) {lnc(0) = -enc(2); lnc(1) = 0.0; lnc(2) = 0.0; return;}      
      if (edid ==  5) {lnc(0) = -enc(2); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid ==  6) {lnc(0) = -enc(2); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid ==  7) {lnc(0) = -enc(2); lnc(1) = 0.0; lnc(2) = 0.0; return;}

      if (edid ==  8) {lnc(0) = -enc(0); lnc(1) = 0.0; lnc(2) = 0.0; return;}      
      if (edid ==  9) {lnc(0) = -enc(1); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid == 10) {lnc(0) =  enc(0); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      if (edid == 11) {lnc(0) =  enc(1); lnc(1) = 0.0; lnc(2) = 0.0; return;}
      print_err(" invalid id of required edge %ld on element %ld (type=%d).\n",
                __FILE__, __LINE__, __func__, edid, eid+1, et);
      abort();
    default:
      print_err("uknown type of element (%d) is required on element %ld.\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }
}



/**
  The function determines general element type of an element entity for the given element.
  The type of element entity is detected according to number of nodes.

  @param eid[in] - required element id,
  @param entnn[in] - number of nodes defining searched entity,
  @param ent[out] - general entity type detected on the element.

  @return The function returns general element type of the detected entity and
          additionally, it returns entity type in the argument ent. 

  Created by Tomas Koudelka, 05.2023.
*/
gtypel siftop::give_elem_ent(long eid, long entnn, gentity &ent) const
{
  gtypel et = elements[eid].type;
  gtypel ret = gtypel(-1);
  if (entnn == 1){
    ent = gvertex;
    return noel;
  }
  switch (et){
    case isolinear1d:
      if (entnn == 2) {ret = isolinear1d;  ent = gregion;}
      break;
    case isoquadratic1d:
      if (entnn == 3) {ret = isoquadratic1d;  ent = gregion;}
      break;
    case isolinear2d:
      if (entnn == 2) {ret = isolinear1d;  ent = gcurve;}
      if (entnn == 4) {ret = isolinear2d;  ent = gregion;}
      break;
    case isoquadratic2d:
      if (entnn == 3) {ret = isoquadratic1d;  ent = gcurve;}
      if (entnn == 8) {ret = isoquadratic2d;  ent = gregion;}
      break;
    case isolinear3d:
      if (entnn == 8) {ret = isolinear3d;  ent = gregion;}
      if (entnn == 4) {ret = isolinear2d;  ent = gsurface;}
      if (entnn == 2) {ret = isolinear1d;  ent = gcurve;}
      break;
    case isoquadratic3d:
      if (entnn == 20) {ret = isoquadratic3d;  ent = gregion;}
      if (entnn == 8)  {ret = isoquadratic2d;  ent = gsurface;}
      if (entnn == 3)  {ret = isoquadratic1d;  ent = gcurve;}
      break;
    default:
      print_err("unknown type of element (%d) is required for element %ld,\n",
                __FILE__, __LINE__, __func__, et, eid+1);
      abort();
  }
  if (ret < 0){
    print_err("cannot find entity for the %ld nodes on element %ld.\n",
              __FILE__, __LINE__, __func__, entnn, eid+1);
    abort();
  }
  
  return ret;
}



/**
  The function creates list of edges from the given element list.

  @param[out] edglst - array of particular edge definitions (arrays of edge nodes), i.e. edglst[i][j] = j-th node on i-th edge,
                       where i is the global edge index. Node numbers involved in the definition of particular edges are sorted.
  @param[out] edgelem - list of adjacent elements to particular edges, i.e edgelem[i][j] = j-th element connected to i-th edge
                         where i is the global edge index. Particular arrays of element numbers are sorted, i.e. edgelem[i] is 
                         sorted vector.
  @param[out] elemedg - list of adjacent edges to particular elements, i.e elemedg[i][j] = global edge index of j-th edge connected 
                        with i-th element. Particular arrays of edge indices are sorted, i.e. elemedg[i] is sorted vector.
  @param[in] elems - list of elements whose edges will be considered in the list, elems[i]=number of element under consideration
  @param[in] pelems - array of all elements in the mesh
  @param[in] pnodes - array of all nodes in the mesh

  Created by TKo, 03.2024
*/
void generate_edge_list(std::vector<std::vector<long>> &edglst,
                        std::vector<std::vector<long>> &edgelem,
                        std::vector<std::vector<long>> &elemedg,
                        const std::vector<long> &elems, const selement *pelems, const snode *pnodes)
{
  std::size_t i;
  long j, k, eid, edgid;
  std::vector<long> edgn;
  std::set<std::vector<long>> edglst_set;
  //
  // create set od egdes
  //
  for (i=0; i<elems.size(); i++){
    eid = elems[i];
    for (j=0; j<pelems[eid].ned; j++){
      // create array of node numbers on the j-th edge of i-th element under consideration
      edgn.clear();
      edgn.reserve(pelems[eid].nned[j]);
      // insert edge nodes in the edgn vector
      for (k=0; k<pelems[eid].nned[j]; k++)
        edgn.push_back(pelems[eid].nodes[pelems[eid].edgenod[j][k]]);
      std::sort(edgn.begin(), edgn.end()); // sort edge node numbers
      edglst_set.insert(edgn); // insert given edge in std::set, duplicate edges are excluded automatically
    }
  } // edglst_set contains definitions of particular edges without duplicates

  //
  // convert edglst_set to std::vector of std::vectors which allows access by indices and creates global edge numbering
  //
  edglst.reserve(edglst_set.size());
  edglst.insert(edglst.begin(), edglst_set.begin(), edglst_set.end());

  typename decltype(edglst_set)::iterator edg1_it;
  //
  // create edge -> adjacent elements correspondence map
  //
  edgelem.reserve(edglst.size());
  for (i=0; i<elems.size(); i++){
    eid = elems[i];
    // insert edge nodes in the edgn vector
    for (j=0; j<pelems[eid].ned; j++){
      // create array of node numbers on the j-th edge of i-th element under consideration
      edgn.clear();
      edgn.reserve(pelems[eid].nned[j]);
      // insert edge nodes in the edgn vector
      for (k=0; k<pelems[eid].nned[j]; k++)
        edgn.push_back(pelems[eid].nodes[pelems[eid].edgenod[j][k]]);
      std::sort(edgn.begin(), edgn.end()); // sort edge node numbers
      auto its = edglst_set.find(edgn);
      if (its != edglst_set.end()){
        edgid = (long) std::distance(its, edg1_it);
        edgelem[edgid].push_back(i);
      }
    }
  }
  for(auto it=edgelem.begin(); it != edgelem.end(); ++it)
    std::sort(it->begin(), it->end()); // sort the list of elements for particular edges

  // create list of edge numbers (indices) on particular elements
  elemedg.reserve(elems.size());
  for(auto edg_it=edgelem.begin(); edg_it != edgelem.end(); ++edg_it){
    for (auto elem_it=edg_it->begin(); elem_it != edg_it->end(); ++elem_it){
      elemedg[*elem_it].push_back((long)std::distance(edg_it, edgelem.begin()));
    }
  }
  for(auto it=elemedg.begin(); it != elemedg.end(); ++it)
    std::sort(it->begin(), it->end()); // sort the list of edges for particular elements
}
