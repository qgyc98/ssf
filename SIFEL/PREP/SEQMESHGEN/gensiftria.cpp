#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/**
   the code generates triangular elements on rectangular domain
   the mesh is in SIFEL format
   
   14.1.2002 tk
   4.3.2007 JK
   6.9.2010 JK
*/
int main (int argc,char *argv[])
{
  long i,j,nx,ny,nn,ne,edg,nid,eid,ordin,np;
  double lx,ly,dx,dy,xx,yy;
  FILE *out;
  char *Fname;

  fprintf(stdout, "\n\n\n---          GENERATOR OF TRIANGULAR ELEMENTS ON RECTANGULAR DOMAIN        ---");
  fprintf(stdout, "\n--- generates triangular elements with three nodes on a rectangular domain ---\n");
  if (argc < 7){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s output_file_name lx ly nx ny write_edge_num\n\n", argv[0]);
    return(1);
  }
  out = fopen(argv[1], "wt");
  if (out==NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  Fname = argv[1];

  if (sscanf(argv[2], "%lf", &lx) != 1)
  {
    fprintf (stderr,"\n Length of rectangle lx has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if (sscanf(argv[3], "%lf", &ly) != 1)
  {
    fprintf (stderr,"\n Height of rectangle ly has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(4);
  }
  if (sscanf(argv[4], "%ld", &nx) != 1)
  {
    fprintf (stderr,"\n Number of elements in x direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(5);
  }
  if (sscanf(argv[5], "%ld", &ny) != 1)
  {
    fprintf (stderr,"\n Number of elements in y direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(6);
  }
  if (sscanf(argv[6], "%ld", &edg) != 1)
  {
    fprintf (stderr,"\n Indicator of edge numbers has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }

  fprintf(stdout, "\n Length of rectangle in x direction   %f", lx);
  fprintf(stdout, "\n Length of rectangle in y direction   %f", ly);
  fprintf(stdout, "\n The number of elements in x direction   %ld", nx);
  fprintf(stdout, "\n The number of elements in y direction   %ld", ny);
  fprintf(stdout, "\n File name for output  %s", Fname);
  fprintf(stdout, "\n Edge numbering    %ld", edg);

  fprintf(stdout, "\n\n Ordinary nodes are denoted by nodal property 0");
  fprintf(stdout, "\n Corner nodes are denoted by nodal properties 1, 2, 3 and 4");
  fprintf(stdout, "\n Nodes lying on edges are marked by edge properties 1, 2, 3 and 4");
  fprintf(stdout, "\n All elements are denoted by property 1");
  fprintf(stdout, "\n Elements with boundary edges are denoted by edge properties 1, 2, 3 and 4");

  fprintf(stdout, "\n ");
  fprintf(stdout, "\n                      E1 ");
  fprintf(stdout, "\n N2 ------------------------------------- N1");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n E2 |              N0 S1 V1             | E4");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n N3 ------------------------------------- N4");
  fprintf(stdout, "\n                      E3 ");
  fprintf(stdout, "\n\n");


  //  the number of nodes
  nn = (nx+1)*(ny+1);
  //  the number of elements
  ne = nx*ny*2;
  //  length of one element in the x-direction
  dx = lx/nx;
  //  length of one element in the y-direction
  dy = ly/ny;
  
  fprintf (out,"%ld",nn);
  fprintf(stdout, "\n The number of nodes      %ld",nn);

  //  node generation
  nid=1;
  xx=0.0;  yy=0.0;
  for (i=0;i<nx+1;i++){
    yy=0.0;
    for (j=0;j<ny+1;j++){
      fprintf (out,"\n%8ld  %15.12le %15.12le %15.12le  ",nid,xx,yy,0.0);
      
      np=0;
      //  node properties
      ordin=0;
      if (i==0 && j==0){
	//  corner 3
	ordin=1;
	np++;
      }
      if (i==0 &&j==ny){
	//  corner 2
	ordin=1;
	np++;
      }
      if (i==nx && j==0){
	//  corner 4
	ordin=1;
	np++;
      }
      if (i==nx &&j==ny){
	//  corner 1
	ordin=1;
	np++;
      }
      if (ordin==0){
	//  ordinary (non-corner) node
	np++;
      }
      
      
      //  properties obtained from edges
      ordin=0;
      if (i==0 && j==0){
	//  corner 3
	ordin=1;
	np+=2;
      }
      if (i==0 && j==ny){
	//  corner 2
	ordin=1;
	np+=2;
      }
      if (i==nx && j==0){
	//  corner 4
	ordin=1;
	np+=2;
      }
      if (i==nx && j==ny){
	//  corner 1
	ordin=1;
	np+=2;
      }
      
      if (ordin==0){
	//  non-corner nodes
	if (i==0){
	  np++;
	}
	if (i==nx){
	  np++;
	}
	if (j==0){
	  np++;
	}
	if (j==ny){
	  np++;
	}
      }
      
      //  surface and volume properties
      np+=2;
      
      fprintf (out,"   %ld  ",np);

      //  node properties
      ordin=0;
      if (i==0 && j==0){
	//  corner 3
	fprintf (out,"  1 3");
	ordin=1;
      }
      if (i==0 &&j==ny){
	//  corner 2
	fprintf (out,"  1 2");
	ordin=1;
      }
      if (i==nx && j==0){
	//  corner 4
	fprintf (out,"  1 4");
	ordin=1;
      }
      if (i==nx &&j==ny){
	//  corner 1
	fprintf (out,"  1 1");
	ordin=1;
      }
      if (ordin==0){
	//  ordinary (non-corner) node
	fprintf (out,"  1 0");
      }
      
      
      //  properties obtained from edges
      ordin=0;
      if (i==0 && j==0){
	//  corner 3
	fprintf (out,"  2 2");
	fprintf (out,"  2 3");
	ordin=1;
      }
      if (i==0 && j==ny){
	//  corner 2
	fprintf (out,"  2 1");
	fprintf (out,"  2 2");
	ordin=1;
      }
      if (i==nx && j==0){
	//  corner 4
	fprintf (out,"  2 3");
	fprintf (out,"  2 4");
	ordin=1;
      }
      if (i==nx && j==ny){
	//  corner 1
	fprintf (out,"  2 1");
	fprintf (out,"  2 4");
	ordin=1;
      }
      
      if (ordin==0){
	//  non-corner nodes
	if (i==0){
	  fprintf (out,"  2 2");
	}
	if (i==nx){
	  fprintf (out,"  2 4");
	}
	if (j==0){
	  fprintf (out,"  2 3");
	}
	if (j==ny){
	  fprintf (out,"  2 1");
	}
      }
      
      fprintf (out,"  3 1  4 1");

      yy+=dy;
      nid++;
    }
    xx+=dx;
  }
  
  //  generation of elements
  eid=1;
  fprintf (out,"\n%ld",ne);
  fprintf(stdout, "\n The number of elements   %ld", ne);
  long zigzag = 0;
  long edgp[3];

  for (i=0;i<nx;i++){
    if (i%2)
      zigzag = 1;
    else
      zigzag = 0;
    for (j=0;j<ny;j++){
      if (zigzag == 0)
      {
        fprintf (out,"\n%8ld 3 %8ld %8ld %8ld       1", eid, (i+1)*(ny+1)+j+2, i*(ny+1)+j+1, (i+1)*(ny+1)+j+1);
        eid++;
        if (edg){
          memset(edgp, 0, sizeof(*edgp)*3);
          if (j==0)
            edgp[1] = 3;
          if (i==nx-1)
            edgp[2] = 4;
          fprintf (out,"  %ld %ld %ld    1", edgp[0], edgp[1], edgp[2]);
        }
      
        fprintf (out,"\n%8ld 3 %8ld %8ld %8ld       1", eid, (i+1)*(ny+1)+j+2, i*(ny+1)+j+2, i*(ny+1)+j+1);
        eid++;
        if (edg){
          memset(edgp, 0, sizeof(*edgp)*3);
          if (i==0)
            edgp[1] = 2;
          if (j==ny-1)
            edgp[0] = 1;
          fprintf (out,"  %ld %ld %ld    1", edgp[0], edgp[1], edgp[2]);
        }
        zigzag = 1;
        continue;
      }
      
      if (zigzag == 1)
      {
        fprintf (out,"\n%8ld 3 %8ld %8ld %8ld       1", eid, i*(ny+1)+j+2, i*(ny+1)+j+1, (i+1)*(ny+1)+j+1);
        eid++;
        if (edg){
          memset(edgp, 0, sizeof(*edgp)*3);
          if (i==0)
            edgp[0] = 2;
          if (j==0)
            edgp[1] = 3;
          fprintf (out,"  %ld %ld %ld    1", edgp[0], edgp[1], edgp[2]);
        }
      
        fprintf (out,"\n%8ld 3 %8ld %8ld %8ld       1", eid, i*(ny+1)+j+2, (i+1)*(ny+1)+j+1, (i+1)*(ny+1)+j+2);
        eid++;
        if (edg){
          memset(edgp, 0, sizeof(*edgp)*3);
          if (i==nx-1)
            edgp[1] = 4;
          if (j==ny-1)
            edgp[2] = 1;
          fprintf (out,"  %ld %ld %ld    1", edgp[0], edgp[1], edgp[2]);
        }
        zigzag = 0;
        continue;
      }
    }
  }

  fprintf(stdout, "\n\n");
  fprintf (out,"\n");

  fclose (out);
}
