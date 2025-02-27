// generator header files
#include "lshape_geom.h"
#include "cprof_geom.h"

// GEFEL header files
#include "gtopology.h"
#include "siftop.h"
#include "nodmap.h"
#include "basefun.h"
#include "difcalc.h"
#include "mathem.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>







lshape_geom::lshape_geom()
{
  lh = lv = 0.0;
  th = tv = 0.0;
  mshdh = mshdv = 0;
  numlh = numlv = 0;

  tlh = tlv = NULL;
  mshdhl = mshdvl = NULL;
  tmath = tmatv = NULL;
  matidh = matidv = NULL;  
  cprof = NULL;
  px = py = NULL;
  alpha = NULL;

  tnn = nn_wall = nn_cprof = 0;
  tne = ne_wall = ne_cprof = 0;
}



lshape_geom::~lshape_geom()
{
  delete [] tlh;
  delete [] tlv;
  delete [] mshdhl;
  delete [] mshdvl;
  delete [] tmath;
  delete [] matidh;
  delete [] tmatv;
  delete [] matidv;
  delete [] cprof;
  delete [] px;
  delete [] py;
  delete [] alpha;
}



/**
  The function reads L-shape layered wall structure geometru from the opened text file.
  All read data is stored in data members of the class.
  
  @param in[in] - pointer to the opened input text file

  @retval 0 - on success
  @retval 1 - in the case the reading of basic geometry failed
  @retval 2 - in the case the reading of horizontal wall part layers failed
  @retval 3 - in the case the reading of vertical wall part layers failed
  @retval 4 - in the case the reading of steel C profiles failed

  Created by Tomas Koudelka, 11.7.2018
*/
int lshape_geom::read(XFILE *in)
{
  int nri=0;
  int i;

  fprintf(stdout, "Reading of basic geometry parameters ...");

  nri += xfscanf(in, " lh %le", &lh);
  nri += xfscanf(in, " lv %le", &lv);

  nri += xfscanf(in, " mesh-dens-h %d", &mshdh);
  nri += xfscanf(in, " mesh-dens-v %d", &mshdv);


  nri += xfscanf(in, " numlh %d", &numlh);
  nri += xfscanf(in, " numlv %d", &numlv);

  if (nri == 6)
    fprintf(stdout, " OK.\n");
  else
  {
    fprintf(stdout, " FAILURE\nnri=%d", nri);
    return 1;
  }

  tlh = new double[numlh];
  memset(tlh, 0, sizeof(*tlh)*numlh);
  tlv = new double[numlv];
  memset(tlh, 0, sizeof(*tlv)*numlv);

  tmath = new tmat[numlh];
  memset(tmath, 0, sizeof(*tmath)*numlh);
  tmatv = new tmat[numlv];
  memset(tmatv, 0, sizeof(*tmatv)*numlv);

  mshdhl = new int[numlh];
  memset(mshdhl, 0, sizeof(*mshdhl)*numlh);
  mshdvl = new int[numlv];
  memset(mshdvl, 0, sizeof(*mshdvl)*numlv);

  matidh = new int[numlh];
  memset(matidh, 0, sizeof(*matidh)*numlh);
  matidv = new int[numlv];
  memset(matidv, 0, sizeof(*matidv)*numlv);

  fprintf(stdout, "Reading of layer parameters of horizontal wall part ...");
  nri = 0;
  th = 0.0;  
  for (i=0; i<numlh; i++)
  {
    nri += xfscanf(in, " thick %le", tlh+i);
    th += tlh[i];
    nri += xfscanf(in, " mesh-dens %d", mshdhl+i);
    nri += xfscanf(in, " type-mat  %d", (int *)(tmath+i));
    nri += xfscanf(in, " mat-id %d", matidh+i);
  }
  if (nri == 4*numlh)
    fprintf(stdout, " OK.\n");
  else
  {
    fprintf(stdout, " FAILURE\nnri=%d", nri);
    return 2;
  }

  fprintf(stdout, "Reading of layer parameters of vertical wall part ...");
  nri = 0;
  tv = 0.0;
  for (i=0; i<numlv; i++)
  {
    nri += xfscanf(in, " thick %le", tlv+i);
    tv += tlv[i];
    nri += xfscanf(in, " mesh-dens %d", mshdvl+i);
    nri += xfscanf(in, " type-mat  %d", (int *)(tmatv+i));
    nri += xfscanf(in, " mat-id %d", matidv+i);
  }
  if (nri == 4*numlv)
    fprintf(stdout, " OK.\n");
  else
  {
    fprintf(stdout, " FAILURE\nnri=%d", nri);
    return 3;
  }

  fprintf(stdout, "Reading of profiles ...");
  nri = xfscanf(in, " num-prof %d", &numprof);
  if (nri == 1)
    fprintf(stdout, " OK.\n");
  else
  {
    fprintf(stdout, " FAILURE in number of profiles\n");
    return 4;
  }
  
  cprof = new cprof_geom[numprof];
  px = new double[numprof];
  py = new double[numprof];
  alpha = new double[numprof];
  memset(px, 0, sizeof(*px)*numprof);
  memset(py, 0, sizeof(*py)*numprof);
  memset(alpha, 0, sizeof(*alpha)*numprof);
  for(i=0; i<numprof; i++)
  {
    nri = 0;
    nri += xfscanf(in, " px %le", px+i);
    nri += xfscanf(in, " py %le", py+i);
    nri += xfscanf(in, " alpha %le", alpha+i);
    alpha[i] *= M_PI/180.0; // conversion to radians
    nri += cprof[i].read(in);
    if (nri != 3)
    {
      fprintf(stdout, " FAILURE in %d. profile\nnri=%d", i+1, nri);
      return 4;
    }
  }
  fprintf(stdout, " OK.\n");

  return 0;
}



/**
  The function generates nodes of structured FE mesh for L-shape layered wall structure.
  Generated nodes are gradually printed into the text file.
  
  @param out[in] - pointer to the opened output text file

  @retval 0 - on success
  @retval 1 - in the case the reading of basic geometry failed
  @retval 2 - in the case the reading of horizontal wall part layers failed
  @retval 3 - in the case the reading of vertical wall part layers failed
  @retval 4 - in the case the reading of steel C profiles failed

  Created by Tomas Koudelka, 11.7.2018
*/
int lshape_geom::gen_nodes(FILE *out)
{
  int i, j, k, ii, jj;

  int nid;       // actual node number of generated node
  double x, y;   // actual values of nodal coordinates of the generated node
  double dx, dy; // actual values of x-coord. and y-coord. increment in block of generated nodes 

  int prop[10];    // array of property id for one generated node
  int propent[10]; // array of property entity type for one generated node
  int nprop;       // the number of properties at one node

  int mdh; // mesh density in the horizontal direction of the actual layer
  int mdv; // mesh density in the vertical direction of the actual layer

  int nncv_v; // number of nodes in vertical direction on corner + vertical part
  int nncv_h; // number of nodes in horizontal direction on corner + vertical part
  int nncv;   // number of nodes in nodal block forming the corner and vertical wall part
  int *cvlid_h = NULL; // layer id for corner and vertical part nodes cvlid_h[nrow_id] = layer id, where nrow_id is the id of nodal row id in the block of all nodes of corner and vertical part
  int *cvlid_v = NULL; // layer id for corner and vertical part nodes cvlid_v[ncol_id] = layer id, where ncol_id is the id of nodal column id in the block of all nodes of corner and vertical part

  int nnh_h; // number of nodal rows of horizontal wall part
  int nnh_v; // number of nodal columns of horizontal wall part
  int *hplid_v = NULL;  // layer id for horizontal part nodes hplid_v[nrow_id] = layer id, where nrow_id is the id of nodal row id in the block of all nodes of horizontal wall part

  x = y = 0.0;

  // number of nodes in vertical direction on corner + vertical part
  nncv_v = mshdv+1;
  for (i=0; i<numlh; i++)
    nncv_v += mshdhl[i];

  // number of nodes in horizontal direction on corner + vertical part
  nncv_h = 1;
  for (i=0; i<numlv; i++)
    nncv_h += mshdvl[i];

  cvlid_h = new int[nncv_h]; // layer id for corner and vertical part nodes cvlid_h[nrow_id] = layer id
  cvlid_v = new int[nncv_v]; // layer id for corner and vertical part nodes cvlid_v[ncol_id] = layer id
  memset(cvlid_h, 0, sizeof(*cvlid_h)*nncv_h);
  memset(cvlid_h, 0, sizeof(*cvlid_h)*nncv_h);

  // set indeces of layers for nodal rows in block of nodes of corner and vertical part
  k = 0;
  for (i=0; i<numlv; i++)
  {
    for(j=0; j<mshdvl[i]; j++)
    {
      cvlid_h[k] = i;
      k++;
    }
  }
  cvlid_h[k] = numlv-1;

  // set indeces of layers for nodal columns in block of nodes of corner and vertical part
  k = 0;
  for (i=0; i<numlh; i++)
  {
    for(j=0; j<mshdhl[i]; j++)
    {
      cvlid_v[k] = i;
      k++;
    }
  }
  for (i=0; i<mshdv+1; i++)
  {
    cvlid_v[k] = numlh;
    k++;
  }

  // total number of nodes in corner and vertical part
  nncv = nncv_h*nncv_v;
 
  // total number of nodes
  nn_wall = nncv+(nncv_v-mshdv)*(mshdh);
  nn_cprof = 0;
  for (i=0; i<numprof; i++)
    nn_cprof += cprof[i].get_nn();
  tnn = nn_wall + nn_cprof;

  fprintf(out, "%d\n", tnn);

  //
  // Generation of nodes for the corner and vertical wall part nodal block.
  //
  nid = 0;
  x = 0.0;
  for (i=0; i<nncv_h; i++)
  {
    ii = cvlid_h[i];
    mdh = mshdvl[ii];  // mesh density in the horizontal direction of the actual layer
    if (ii == numlv-1)
      mdh++;
    dx = tlv[ii]/mshdvl[ii]; // increment of x-coordinate
    y = 0.0;
    for (j=0; j<nncv_v; j++)
    {
      jj = cvlid_v[j];
      if (jj < numlh)    
      {
        mdv = mshdhl[jj]; // mesh density in the vertical direction of the actual layer
        dy = tlh[jj]/mdv; // increment of y-coordinate
      }
      else
      {
        mdv = mshdv; // mesh density in the vertical direction of the vertical wall part
        dy = lv/mdv; // increment of y-coordinate
      }
      
      memset(prop, 0, sizeof(*prop)*10);
      memset(propent, 0, sizeof(*propent)*10);
      nprop = 0;

      fprintf(out, "%d %le %le %le ", nid+1, x, y, 0.0);

      if (nid%nncv_v == 0) // outer bottom edge = prop 1
      {
        prop[nprop] = 1;
        propent[nprop] = 2;
        nprop++;
      }
      if (nid<nncv_v) // outer left edge = prop 6
      {
        prop[nprop] = 6;
        propent[nprop] = 2;
        nprop++;
      }
      if (nid%nncv_v == nncv_v-1) // top section edge = prop 5
      {
        prop[nprop] = 5;
        propent[nprop] = 2;
        nprop++;
      }
      if (nid >= nncv-mshdv-1) // inner right edge = prop 4
      {
        prop[nprop] = 4;
        propent[nprop] = 2;
        nprop++;
      }
      if (nid == nncv-mshdv-1) // inner corner - intersection of edges 3 and 4 = prop 3
      {
        prop[nprop] = 3;
        propent[nprop] = 2;
        nprop++;
      }
      if (jj >= ii) // assign material id from vertical part for domains on/above the corner diagonal
      {
        prop[nprop] = tmatv[ii]*10+matidv[ii];
        propent[nprop] = 4;
        nprop++;
      }
      else // assign material id from horizontal part for domains under the corner diagonal
      {
        prop[nprop] = tmath[jj]*10+matidh[jj];
        propent[nprop] = 4;
        nprop++;
      }

      fprintf(out, " %d", nprop+1);
      for(k=0; k<nprop; k++)
        fprintf(out, " %d %d", propent[k], prop[k]);
      fprintf(out, " 3 0\n");

      nid++;
      y += dy;
    }
    x += dx;
  }


  // number of nodes in horizontal direction of horizontal wall part block of nodes
  nnh_h = mshdh;
  // number of nodes in vertical direction of horizontal wall part block of nodes
  nnh_v = nncv_v-mshdv;

  // set indeces of layers for nodal columns in block of nodes of horizontal wall part
  hplid_v = new int[nnh_v];
  memset(hplid_v, 0, sizeof(*hplid_v)*nnh_v);
  k = 0;
  for (i=0; i<numlh; i++)
  {
    for(j=0; j<mshdhl[i]; j++)
    {
      hplid_v[k] = i;
      k++;
    }
  }
  // correction of initial x coordinate
  x -= dx;

  //
  // Generation of nodes for the horizontal wall part nodal block.
  //
  dx = lh/mshdh;
  for(i=0; i<nnh_h; i++)
  {
    x += dx;
    y = 0.0;
    for(j=0; j<nnh_v; j++)
    {
      jj = hplid_v[j];  // layer index
      mdv = mshdhl[jj]; // mesh density in vertical diraction of the given layer
      dy = tlh[jj]/mdv; // actual y-coord. increment

      memset(prop, 0, sizeof(*prop)*10);
      memset(propent, 0, sizeof(*propent)*10);
      nprop = 0;
      fprintf(out, "%d %le %le %le ", nid+1, x, y, 0.0);
      if (j==0) // nodes at bottom edge = prop 1
      {
        prop[nprop] = 1;
        propent[nprop] = 2;
        nprop++;
      }
      if (j==nnh_v-1) // nodes at top edge of the horizontal wall part = prop 3
      {
        prop[nprop] = 3;
        propent[nprop] = 2;
        nprop++;
      }
      if (i==nnh_h-1) // node at right edge of the horizontal wall part = prop 2
      {
        prop[nprop] = 2;
        propent[nprop] = 2;
        nprop++;
      }
      
      // assign material id from horizontal wall part layers
      prop[nprop] = tmath[jj]*10+matidh[jj];
      propent[nprop] = 4;
      nprop++;

      fprintf(out, " %d", nprop+1);
      for(k=0; k<nprop; k++)
        fprintf(out, " %d %d", propent[k], prop[k]);
      fprintf(out, " 3 0\n");

      nid++;
      y += dy;
    }
  }

  //
  // Generation of nodes for steel profiles.
  //
  for(i=0; i<numprof; i++)
  {
    cprof[i].gen_nodes(out, px[i], py[i], alpha[i], nid);
    nid += cprof[i].get_nn();
  }

  delete [] cvlid_h;
  delete [] cvlid_v;
  delete [] hplid_v;

  return 0;
}



/**
  The function generates elements of structured FE mesh for L-shape layered wall structure.
  Generated elements are gradually printed into the text file.
  
  @param out[in] - pointer to the opened output text file

  @retval 0 - on success
  @retval 1 - in the case of failure of writing into the ouput file

  Created by Tomas Koudelka, 11.7.2018
*/
int lshape_geom::gen_elements(FILE *out)
{
  int i, j, k, ii, jj, nid;

  int eid;        // actual element number of generated element

  int propedg[4]; // array of edge property id for one generated element

  int mdh; // mesh density in the horizontal direction of the actual layer
  int mdv; // mesh density in the vertical direction of the actual layer

  int necv_v; // number of elements in vertical direction of corner + vertical part
  int necv_h; // number of elements in horizontal direction of corner + vertical part
  int necv;   // number of elements in element block forming the corner and vertical wall part
  int *cvlid_h = NULL; // layer id for corner and vertical part elements cvlid_h[erow_id] = layer id, where erow_id is the id of element row id in the block of all elements of corner and vertical part
  int *cvlid_v = NULL; // layer id for corner and vertical part elements cvlid_v[ecol_id] = layer id, where ecol_id is the id of element column id in the block of all elements of corner and vertical part

  int neh_h; // number of element rows of horizontal wall part
  int neh_v; // number of element columns of horizontal wall part
  int *hplid_v = NULL;  // layer id for horizontal part elements hplid_v[erow_id] = layer id, where erow_id is the id of element row in the block of all elements of horizontal wall part

  // number of elements in vertical direction on corner + vertical wall part
  necv_v = mshdv;
  for (i=0; i<numlh; i++)
    necv_v += mshdhl[i];

  // number of elements in horizontal direction on corner + vertical wall part
  necv_h = 0;
  for (i=0; i<numlv; i++)
    necv_h += mshdvl[i];

  cvlid_h = new int[necv_h]; // layer id for corner and vertical part elements cvlid_h[erow_id] = layer id
  cvlid_v = new int[necv_v]; // layer id for corner and vertical part elements cvlid_v[ecol_id] = layer id
  memset(cvlid_h, 0, sizeof(*cvlid_h)*necv_h);
  memset(cvlid_h, 0, sizeof(*cvlid_h)*necv_h);

  // set indeces of layers for element rows in block of elements of corner and vertical wall part
  k = 0;
  for (i=0; i<numlv; i++)
  {
    for(j=0; j<mshdvl[i]; j++)
    {
      cvlid_h[k] = i;
      k++;
    }
  }

  // set indeces of layers for element columns in block of elements of corner and vertical wall part
  k = 0;
  for (i=0; i<numlh; i++)
  {
    for(j=0; j<mshdhl[i]; j++)
    {
      cvlid_v[k] = i;
      k++;
    }
  }
  for (i=0; i<mshdv; i++)
  {
    cvlid_v[k] = numlh;
    k++;
  }

  // total number of elements in corner and vertical wall part
  necv = necv_h*necv_v;
 
  // total number of elements
  ne_wall = necv+(necv_v-mshdv)*mshdh;
  ne_cprof = 0;
  for (i=0; i<numprof; i++)
    ne_cprof += cprof[i].get_ne();
  tne = ne_wall + ne_cprof;

  if (fprintf(out, "\n%d\n", tne) < 0)
  {
    delete [] cvlid_h;
    delete [] cvlid_v;
    return 1;
  }
  

  //
  // Generation of elements for the corner and vertical wall part element block.
  //
  eid = 0;
  for (i=0; i<necv_h; i++)
  {
    ii = cvlid_h[i];
    mdh = mshdvl[ii];  // mesh density in the horizontal direction of the actual layer
    for (j=0; j<necv_v; j++)
    {
      jj = cvlid_v[j];
      if (jj < numlh)    
        mdv = mshdhl[jj]; // mesh density in the vertical direction of the actual layer
      else
        mdv = mshdv; // mesh density in the vertical direction of the vertical wall part
      
      memset(propedg, 0, sizeof(*propedg)*4);
      
      fprintf(out, "%d 5 %d %d %d %d", eid+1, i*(necv_v+1)+j+1, (i+1)*(necv_v+1)+j+1, (i+1)*(necv_v+1)+j+2, i*(necv_v+1)+j+2);
      
      if (jj >= ii) 
        // assign material id from vertical part for domains on/above the corner diagonal
        fprintf(out, "   %d", tmatv[ii]*10+matidv[ii]);
      else 
        // assign material id from horizontal part for domains under the corner diagonal
        fprintf(out, "   %d", tmath[jj]*10+matidh[jj]);

      if (i == 0)   // left edge of the corner and vertical wall part = prop 6
        propedg[3] = 6;

      if (j == 0)  // bottom edge of the domain = prop 1
        propedg[0] = 1;

      if (j == necv_v-1) // top edge of vertical wall part = prop 5
        propedg[2] = 5;

      if ((i == necv_h-1) && (j >= necv_v-mshdv)) // inner edge on vertical wall part = prop 4
        propedg[1] = 4;

      fprintf(out, "   %d %d %d %d   0\n", propedg[0], propedg[1], propedg[2], propedg[3]);

      eid++;
    }
  }

  // number of elements in horizontal direction of horizontal wall part block of elements
  neh_h = mshdh;
  // number of elements in vertical direction of horizontal wall part block of elements
  neh_v = necv_v-mshdv;
  // set indeces of layers for nodal columns in block of nodes of horizontal wall part
  hplid_v = new int[neh_v];
  memset(hplid_v, 0, sizeof(*hplid_v)*neh_v);
  k = 0;
  for (i=0; i<numlh; i++)
  {
    for(j=0; j<mshdhl[i]; j++)
    {
      hplid_v[k] = i;
      k++;
    }
  }

  //
  // Generation of elements for the horizontal wall part element block.
  //
  for(i=0; i<neh_h; i++)
  {
    for(j=0; j<neh_v; j++)
    {
      jj = hplid_v[j];  // layer index
      mdv = mshdhl[jj]; // mesh density in vertical diraction of the given layer

      memset(propedg, 0, sizeof(*propedg)*4);

      if (i == 0)
        // starting index of node at left bottom corner of the horizontal wall part
        nid = necv_h*(necv_v+1);
      else
        // starting index of node at the second nodal column of the horizontal wall part
        nid = (necv_h+1)*(necv_v+1);

      if (i == 0)
        // nid is now index of the first node at the second nodal column formning horizontal wall part
        fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+j+1, nid+(necv_v+1)+j+1, nid+(necv_v+1)+j+2, nid+j+2);
      else        
        // nid is now index of the first node at the second nodal column formning horizontal wall part
        fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+(i-1)*(neh_v+1)+j+1, nid+i*(neh_v+1)+j+1, nid+i*(neh_v+1)+j+2, nid+(i-1)*(neh_v+1)+j+2);

      // assign material id from horizontal part for domains under the corner diagonal
      fprintf(out, "    %d", tmath[jj]*10+matidh[jj]);

      if (j == 0) // bottom edge of the horizontal wall part
        propedg[0] = 1;

      if (j == neh_v-1) // top edge of the horizontal wall part
        propedg[2] = 3;

      if (i == neh_h-1) // right edge of horizontal wall part
        propedg[1] = 2;

      fprintf(out, "   %d %d %d %d   0\n", propedg[0], propedg[1], propedg[2], propedg[3]);

      eid++;
    }
  }
  
  //
  // Generation of steel profiles
  //
  nid = nn_wall;
  for(i=0; i<numprof; i++)
  {
    cprof[i].gen_elements(out, nid, eid);
	eid += cprof[i].get_ne();
	nid += cprof[i].get_nn();
  }

  delete [] cvlid_h;
  delete [] cvlid_v;
  delete [] hplid_v;

  return 0;
}



/**
  The function detects hanging nodes of C-profiles in mesh of L-shape layered wall structure.
  FE mesh of the whole problem is stored in the gtopology class gt where the nodes and elements.
  Resulting natural coordinates and element id for each node of C-profiles are stored in argument 
  cpnmap.
  
  @param gt[in] - pointer to the topology with geometry of 2D FE mesh
  @param wne[in] - the number of elements of layered wall structure without C-profile elements
  @param nnmap[in] - the number of nodes defining FE mesh of C-profiles (dimension of cpnmap)
  @param cpnmap[inout] - array of mapping for C-profile nodes on L-shape wall elements, 
                         node id in gtopology class and global coordinates of nodes must be given 
                         on the input, element id and natural coordinates are stored on the output
  @param ni[in] - the maximum number of iterations used in the natural coordinate computation
  @param err[in] - tolerance used in itartion solver for the natural coordinate computation
                   or maximum distance at which two points are assumed to be identical

  @retval 0 - on success
  @retval 1 - in the case of wrong elements used in gtopology
  @retval 2 - in the case of convergence problem in natural coordinate computation

  Created by Tomas Koudelka, 13.7.2018
*/
int lshape_geom::cprof_hang_nodes_natcoords(gtopology *gt, long wne, long nnmap, nodmap *cpnmap, long ni, double err)
{
  long i, j, nne;
  vector cg(ASTCKVEC(3));
  vector x, y, coord(ASTCKVEC(3));
  vector bf;
  double d2, r2, xi, eta, cerr, err2 = err*err, aerr;
  long nrerr;
  char *ips = new char[nnmap];
  memset(ips, 0, sizeof(*ips)*nnmap);

  for(i=0; i<wne; i++)
  {
    nne = gt->give_nne(i);
    reallocv(RSTCKVEC(nne, x));
    reallocv(RSTCKVEC(nne, y));
    reallocv(RSTCKVEC(nne, bf));

    // compute center of gravity
    gt->give_node_coord2d(x, y, i);
    bf_lin_4_2d(bf.a, 0.0, 0.0); 
    scprd(bf, x, cg(0));
    scprd(bf, y, cg(1));

    r2 = gt->max_sqrdist_nod_pt(i, cg);     
    for(j=0; j<nnmap; j++)
    {
      if (ips[j])
        continue;
      d2  = sqr(cpnmap[j].x-cg(0)) + sqr(cpnmap[j].y-cg(1)) + sqr(cpnmap[j].z-cg(2));
      if (d2 <= r2+err2)
      {
        xi = 0.0; eta = 0.0;
        nrerr = point_natcoord_2d(cpnmap[j].x, cpnmap[j].y, x, y, ni, err, xi, eta, aerr);
        cpnmap[j].zeta = 0.0;
        if (nrerr > ni)
        {
          fprintf(stderr, "\n\nError: solution of natural coordinates for C-profile node %ld on wall element %ld did not converge in %ld steps\n", j+1, i+1, nrerr);
          delete [] ips;
          return 2;
        }
        if ((xi >= -1.0) && (xi <= 1.0) && 
            (eta >= -1.0) && (eta <= 1.0)) 
        {
          ips[j] = 1;
          cpnmap[j].xi = xi;
          cpnmap[j].eta = eta;
          cpnmap[j].cerr = 0.0;
          cpnmap[j].eid = i;
        }
        else
        {
          cerr = 0.0;
          if (fabs(xi)-1.0 > 0.0)
            cerr += sqr(fabs(xi)-1.0);
          if (fabs(eta)-1.0 > 0.0)
            cerr += sqr(fabs(eta)-1.0);
          if (cpnmap[j].cerr > cerr)
          {
            cpnmap[j].xi = xi;
            cpnmap[j].eta = eta;
            cpnmap[j].cerr = cerr;
            cpnmap[j].eid = i;
          }
        }
      }
    }
  }
  delete [] ips;

  return 0;
}



/**
  The function prints hanging nodes forming the C-profiles to the text file.

  @param out[in] - pointer to the opened output text file
  @param gt[in] - pointer to the topology with geometry of 2D FE mesh containing master nodes/elements
  @param nnmap[in] - the number of nodes defining FE mesh of C-profiles (dimension of cpnmap)
  @param cpnmap[inout] - array of mapping for C-profile nodes on L-shape wall elements, 
                         node id in gtopology class and global coordinates of nodes must be given 
                         on the input, element id and natural coordinates are stored on the output

  @retval 0 - on success
  @retval 1 - in the case of failure of writing into the ouput file

  Created by Tomas Koudelka, 13.7.2018
*/
int lshape_geom::print_hang_nodes(FILE *out, gtopology *gt, long nnmap, nodmap *cpnmap)
{
  long i, j, nne, nri;
  ivector enod;
 
  nri = fprintf(out, "%ld\n", nnmap);
  if (nri < 0)
    return 1;
 
  for(i=0; i<nnmap; i++)
  {
    nne = gt->give_nne(cpnmap[i].eid);
    reallocv(RSTCKIVEC(nne, enod));
    gt->give_nodes(cpnmap[i].eid, enod);
    fprintf(out, "%ld %ld", cpnmap[i].nid+1, -nne);
    for (j=0; j<nne; j++)
      fprintf(out, " %ld", enod[j]+1);
    fprintf(out, " %le %le %le %d\n", cpnmap[i].xi, cpnmap[i].eta, cpnmap[i].zeta, int(isolinear2d));
  }

  return 0;
}
