#include "cprof_geom.h"
#include "iotools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


cprof_geom::cprof_geom()
{
  wp = hp = 0.0;
  tf = tw = 0.0;
  mshdf = mshdw = 0;
  matid = 0;
}



cprof_geom::~cprof_geom()
{
}



int cprof_geom::read(XFILE *in)
{
  int nri = 0;

  nri += xfscanf(in, " prof-width %le", &wp);
  nri += xfscanf(in, " prof-height %le", &hp);

  nri += xfscanf(in, " thick-flange %le", &tf);
  nri += xfscanf(in, " thick-web %le", &tw);

  nri += xfscanf(in, " mesh-dens-flange %d", &mshdf);
  nri += xfscanf(in, " mesh-dens-web %d", &mshdw);

  nri += xfscanf(in, " mesh-dens-thick-flange %d", &mshdtf);
  nri += xfscanf(in, " mesh-dens-thick-web %d", &mshdtw);

  nri += xfscanf(in, " mat-id %d", &matid);

  if (nri == 9)
    return 0;

  return 1;
}



int cprof_geom::get_nn()
{
  int nn;

  // total number of nodes
  nn  = (mshdw + 2*mshdtf + 1)*(mshdtw + 1);   // number of nodes in web
  nn += 2*(mshdtf + 1)*mshdf; // number of nodes in flanges

  return nn;
}



int cprof_geom::get_ne()
{
  int ne;

  //total number of elements
  ne  = mshdtw*(mshdw+2*mshdtf);
  ne += mshdf*2*mshdtf;

  return ne;
}


void cprof_geom::gen_nodes(FILE *out, double tx, double ty, double alpha, int inid)
{
  int i, j;
  double x, y;   // coordinates of node at normal position
  double xr, yr; // rotated coordinates
  double dx, dy; // coordinate increments for block node generation
  double c = cos(alpha);
  double s = sin(alpha);
  double hw;  // web height
  int nid = inid; // actual node id starts from inid
  int nnx, nny;

  
  // Generation of block of web nodes
  x = 0.0;
  nnx = mshdtw + 1; // number of nodal columns in the generated block
  nny = mshdw + 2*mshdtf + 1; // number of nodal rows in the generated cblock
  dx = tw/mshdtw;  // x-coordinate increment
  for (i=0; i<nnx; i++)
  {
    y = 0.0;
    for(j=0; j<nny; j++)
    {
      if ((j<mshdtf) || (j >= mshdtf+mshdw))
        dy = tf/mshdtf;   // for nodes of flanges
      else
        dy = (hp - 2.0*tf)/mshdw;  // for nodes of web
      
      // translation and rotation of nodal coordinate around the origin by alpha angle
      xr = tx + c*x - s*y;
      yr = ty + s*x + c*y;
      
      fprintf(out, "%d %le %le %le 2 3 0 4 %d\n", nid+1, xr, yr, 0.0, matid);

      y += dy;
      nid++;
    }
    x += dx;
  }

  // Generation of flange nodes
  x -= dx;
  nnx = mshdf;  // number of nodal columns in generated flange nodal blocks
  nny = 2*(mshdtf+1); // number of nodal rows in generated flage nodal blocks
  dx = (wp-tw)/mshdf; // x-coord. increment  
  dy = tf/mshdtf;     // y-coord. increment
  hw = hp-2.0*tf;         // height of web
  x += dx;
  for (i=0; i<nnx; i++)
  {
    y = 0.0;
    for(j=0; j<nny; j++)
    {
      if (j == mshdtf+1) // shift nodes of the second flange
        y += hw-dy;

      xr = tx + c*x - s*y;
      yr = ty + s*x + c*y;
      
      fprintf(out, "%d %le %le %le 2 3 0 4 %d\n", nid+1, xr, yr, 0.0, matid);

      y += dy;
      nid++;
    }
    x += dx;
  }
}



void cprof_geom::gen_elements(FILE *out, int inid, int ieid)
{
  int i, j;
  int eid = ieid;
  int nex, ney;
  int nexw, neyw;
  int nid = inid;

  // Generation of block of web elements
  nex = mshdtw; // number of element columns in the generated block
  ney = mshdw + 2*mshdtf; // number of element rows in the generated cblock
  for (i=0; i<nex; i++)
  {
    for(j=0; j<ney; j++)
    {      
      fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+i*(ney+1)+j+1, nid+(i+1)*(ney+1)+j+1, nid+(i+1)*(ney+1)+j+2, nid+i*(ney+1)+j+2);
      fprintf(out, "   %d   0 0 0 0   0\n", matid);
      eid++;
    }
  }

  //
  // Generation of flange elements
  //
  nexw = nex;
  neyw = ney;
  nex = mshdf;  // number of element columns in generated flange element blocks
  ney = 2*mshdtf; // number of element rows in generated flage element blocks
  for (i=0; i<nex; i++)
  {
    if (i == 0)
      nid = inid + (nexw)*(neyw+1);
    else
      nid = inid + (nexw+1)*(neyw+1);
    for(j=0; j<ney; j++)
    {
      if (i == 0)
      {
        if (j < mshdtf) // first lower flange
          fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+j+1, nid+(neyw+1)+j+1, nid+(neyw+1)+j+2, nid+j+2);
        else // second upper flange
          fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+mshdw+j+1, nid+(neyw+1)+j+2, nid+(neyw+1)+j+3, nid+mshdw+j+2);
      }
      else
      {
        if (j < mshdtf) // first lower flange
          fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+(i-1)*(ney+2)+j+1, nid+i*(ney+2)+j+1, nid+i*(ney+2)+j+2, nid+(i-1)*(ney+2)+j+2);
        else // second upper flange
          fprintf(out, "%d 5 %d %d %d %d", eid+1, nid+(i-1)*(ney+2)+j+2, nid+i*(ney+2)+j+2, nid+i*(ney+2)+j+3, nid+(i-1)*(ney+2)+j+3);
      }

      fprintf(out, "   %d   0 0 0 0   0\n", matid);
      eid++;
    }
  }
}
