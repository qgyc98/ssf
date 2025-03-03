#include <stdio.h>
#include <string.h>
#include "siftop.h"



/*
  Generator of 3D domain from 2D section topology

  Created by Tomas Koudelka, 5.10.2018
*/
int main (int argc,char *argv[])
{
  long i, j, k, l, m, aux;
  meshform fmt;
  long edg, paral, nsec;
  char dir;
  double dx,dy,dz;
  double *d = NULL;
  XFILE *in;
  FILE *out;
  long ret;
  char topname[1001];
  siftop top;
  double *secl, *dl;
  long totnsec, *md;
  long frontp, backp;
  long nid, nid1, nid2, eid;
  selement isolin3d;

  fprintf(stdout, "\n\n\n---        3D MESH GENERATION FROM 2D DOMAIN SECTION         ---\n\n");
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
  xfscanf(in, " %1000a", topname);
  xfscanf(in, "%k%ld %k%ld", "edg", &edg, "paral", &paral);
  xfscanf(in, "%k %c %k%ld %k%ld %k%ld", "gen_dir", &dir, "front_prop", &frontp, "rear_prop", &backp, "num_sec", &nsec);
  dl = new double[nsec];
  memset(dl, 0, sizeof(*dl)*nsec);
  md = new long[nsec];
  memset(md, 0, sizeof(*md)*nsec);
  secl = new double[nsec];
  memset(secl, 0, sizeof(*secl)*nsec);
  totnsec = 1; // total number of 2D sections generated in the 3D domain
  for (i=0; i<nsec; i++)
  {
    xfscanf(in, "%k%le %k%ld", "secblock_lenght", secl+i, "mesh_dens", md+i);
    dl[i] = secl[i]/md[i];
    totnsec += md[i];
  }
  xfclose(in);

  // open and read file with topology of 2D section
  fprintf(stdout, "Reading of topology file %s with 2D section ...", topname);
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
  fprintf(stdout, "OK\n");
  fprintf(stdout, "\nNumber of nodes in 2D section   : %ld\n", top.nn);
  fprintf(stdout, "Number of elements in 2D section: %ld\n", top.ne);
  fprintf(stdout, "\nGeneration of 3D domain ...");

  // print nodes of intial (first) 2D section
  fprintf(out, "%ld\n", top.nn*(totnsec));
  for (i=0; i<top.nn; i++)
  {
    fprintf (out,"%ld %15.10le %15.10le %15.10le %ld",i+1, top.nodes[i].x, top.nodes[i].y, top.nodes[i].z, top.nodes[i].nprop+1);
    for(j=0; j<top.nodes[i].nprop; j++)
    {
      if (top.nodes[i].entid[j] == ecurve)
        fprintf(out, " 3"); // change edge property to surface property
      else
        fprintf(out, " %d", top.nodes[i].entid[j]);
      fprintf(out, " %ld", top.nodes[i].prop[j]);
    }
    fprintf(out, " 3 %ld\n", frontp); // add required front face property
  }
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
  switch(dir)
  {
    case 'x':
      d = &dx;
      break;
    case 'y':
      d = &dy;
      break;
    case 'z':
      d = &dz;
      break;
    default:
      print_err("invalid direction indicator %d, it must be in the range <1;3>", __FILE__, __LINE__, __func__, dir);
      return 4;
  }
  nid = top.nn;
  for (i=0; i<nsec; i++)
  {
    for (j=0; j<md[i]; j++)
    {
      *d += dl[i];
      for (k=0; k<top.nn; k++)
      {
        fprintf (out,"%ld %15.10le %15.10le %15.10le",k+1+nid, top.nodes[k].x+dx, top.nodes[k].y+dy, top.nodes[k].z+dz);
        if ((i == nsec-1) && (j == md[i]-1)) // for nodes of the last 2D section
          fprintf (out," %ld", top.nodes[k].nprop+1);
        else
          fprintf (out," %ld", top.nodes[k].nprop);
        for(l=0; l<top.nodes[k].nprop; l++)
        {
          if (top.nodes[k].entid[l] == ecurve)
            fprintf(out, " 3"); // change edge property to surface property
          else
            fprintf(out, " %d", top.nodes[k].entid[l]);
          fprintf(out, " %ld", top.nodes[k].prop[l]);
        }
        if ((i == nsec-1) && (j == md[i]-1)) // for nodes of the last 2D section
          fprintf(out, " 3 %ld\n", backp); // add required front face property
        else
          fprintf(out, "\n");
      }
      nid += top.nn;
    }
  }  
  fprintf(out, "%ld\n", top.ne*(totnsec-1));
  eid = 0;
  nid1 = 0;
  nid2 = top.nn;
  isolin3d.type = isolinear3d;
  isolin3d.alloc(0);
  for (i=0; i<nsec; i++)
  {   
    for (j=0; j<md[i]; j++)
    {
      for (k=0; k<top.ne; k++)
      {
        if (top.elements[k].type == isolinear2d)
        {      
          fprintf(out, "%ld %d ", k+1+eid, int(isolinear3d));
          for (l=0; l<top.elements[k].nne; l++) 
            fprintf(out, "%ld ", top.elements[k].nodes[l]+1+nid2);
          for (l=0; l<top.elements[k].nne; l++) 
            fprintf(out, "%ld ", top.elements[k].nodes[l]+1+nid1);
          fprintf(out, "%ld", top.elements[k].prop);
          if (edg)
          {
            // edge property id of brick element are not preserved, 
            // they will be changed to surface property id 
            fprintf(out, "  ");
            for (l=0; l<isolin3d.ned; l++)
              fprintf(out, " %d", 0);
            fprintf(out, "  ");

            // property id of side surfaces on the brick element
            for (l=0; l<top.elements[k].ned; l++)
            {
              // indeces of brick side surfaces are shiftted by 3 with respect 
              // to edge indeces on quadrilateral
              m = (l+3)%top.elements[k].ned; 
              fprintf(out, " %ld", top.elements[k].propedg[m]);
            }

            // property id of rear base surface of the brick
            aux = 0;
            if ((i == nsec-1) && (j == md[i]-1))
              aux = backp;
            fprintf(out, " %ld", aux);

            // property id of front base surface of the brick
            aux = 0;
            if ((i == 0) && (j == 0))
              aux = frontp;
            fprintf(out, " %ld", aux);
          }
          fprintf(out, "\n");
        }
        if (top.elements[k].type == trianglelinear) 
        {
          fprintf(out, "%ld %d ", k+1+eid, int(isolinear3d));
          for (l=0; l<top.elements[k].nne; l++) 
            fprintf(out, "%ld ", top.elements[k].nodes[l]+1+nid2);
          fprintf(out, "%ld ", top.elements[k].nodes[top.elements[k].nne-1]+1+nid2);
          for (l=0; l<top.elements[k].nne; l++) 
            fprintf(out, "%ld ", top.elements[k].nodes[l]+1+nid1);
          fprintf(out, "%ld ", top.elements[k].nodes[top.elements[k].nne-1]+1+nid1);
          if (edg)
          {
            // edge property id of brick element are not preserved, 
            // they will be changed to surface property id 
            for (l=0; l<isolin3d.ned; l++)
              fprintf(out, " %d", 0);

            fprintf(out, " 0"); // surface with zero area due to degenerated brick element
            // property id of side surfaces on the brick element
            for (l=0; l<top.elements[k].ned; l++)
              fprintf(out, " %ld", top.elements[k].propedg[l]);

            // property id of rear base surface of the brick
            aux = 0;
            if ((i == nsec-1) && (j == md[i]-1))
              aux = backp;
            fprintf(out, " %ld", aux);

            // property id of front base surface of the brick
            aux = 0;
            if ((i == 0) && (j == 0))
              aux = frontp;
            fprintf(out, " %ld", aux);
          }
          fprintf(out, "\n");
        }
      }
      eid += top.ne;
      nid1 += top.nn;  // increase start index of nodes on the front base surface 
      nid2 += top.nn;  // increase start index of nodes on the rear base surface 
    }
  }
  fclose(out);

  fprintf(stdout, "OK\n\n");
  fprintf(stdout, "The total number of nodes in 3D domain   : %ld\n", top.nn*totnsec);
  fprintf(stdout, "The total number of elements in 3D domain: %ld\n", top.ne*(totnsec-1));
  fprintf(stdout, "3D topology has been written in file %s\n", argv[2]);
  return 0;
}
