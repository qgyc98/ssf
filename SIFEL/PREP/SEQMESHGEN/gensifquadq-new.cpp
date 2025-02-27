#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc,char *argv[])
/*
   14.1.2002 tk
   4.3.2007 JK
*/
{
  long i,j,nx,ny,nn,ne,edg,nid;
  double lx,ly,dx1,dx2,dy1,dy2,xx,yy;
  FILE *s1;
  char *Fname;

  fprintf(stdout, "\n\n\n---          GENERATOR OF RECTANGULAR ELEMENTS ON RECTANGULAR DOMAIN        ---");
  fprintf(stdout, "\n--- generates rectangular elements with eight nodes on a rectangular domain ---\n");
  if (argc < 7){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s output_file_name lx ly nx ny write_edge_num\n\n", argv[0]);
    return(1);
  }
  s1 = fopen(argv[1], "wt");
  if (s1==NULL){
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
    fprintf (stderr,"\n Heigth of rectangle ly has not been specified.");
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

  fprintf(stdout, "\n Length of rectangle in x direction    %f", lx);
  fprintf(stdout, "\n Length of rectangle in y direction    %f", ly);
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


  dx1 = lx/nx;    dy1 = ly/ny;
  dx2 = dx1/2.0;  dy2 = dy1/2.0;
  ne = nx*ny;
  nn = (nx+1)*(ny*2+1)+nx*(ny+1);

  fprintf (s1,"%ld",nn);
  fprintf(stdout, "\n The number of nodes     %ld", nn);
  
  nid=1;
  xx=0.0;  yy=0.0;
  for (i=0;i<=ny*2;i++){
    fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le",nid,xx,yy,0.0);
    yy+=dy2;
    nid++;
    if (i==0){
      //  left bottom corner
      fprintf (s1,"  5  1 3  2 2  2 3  3 1  4 1");
    }
    if (i==ny*2){
      //  left top corner
      fprintf (s1,"  5  1 2  2 1  2 2  3 1  4 1");
    }
    if (i!=0 && i<ny*2){
      fprintf (s1,"  4  1 0  2 2  3 1  4 1");
    }
  }

  xx=dx2;
  for (i=0;i<nx;i++){
    yy=0.0;
    for (j=0;j<=ny;j++){
      if (j == 0){
	fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  4  1 0  2 3  3 1  4 1",nid,xx,yy,0.0);
	nid++;
      }
      if (j == ny){
	fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  4  1 0  2 1  3 1  4 1",nid,xx,yy,0.0);
	nid++;
      }
      if (j!=0 && j<ny){
	fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  3  1 0  3 1  4 1",nid,xx,yy,0.0);
	nid++;
      }
      yy+=dy1;
    }
    xx+=dx2;
    
    if (i!=nx-1){
      yy=0.0;
      for (j=0;j<=ny*2;j++){
	if (j == 0){
	  fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  4  1 0  2 3  3 1  4 1",nid,xx,yy,0.0);
	  nid++;
	}
	if (j == 2*ny){
	  fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  4  1 0  2 1  3 1  4 1",nid,xx,yy,0.0);
	  nid++;
	}
	if (j!=0 && j<2*ny){
	  fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le  3  1 0  3 1  4 1",nid,xx,yy,0.0);
	  nid++;
	}
	yy+=dy2;
      }
      xx+=dx2;
    }
  }
  
  yy=0.0;
  for (i=0;i<=ny*2;i++){
    fprintf (s1,"\n%8ld  %15.10le %15.10le %15.10le",nid,xx,yy,0.0);
    yy+=dy2;
    nid++;
    if (i==0){
      //  right bottom corner
      fprintf (s1,"  5  1 4  2 3  2 4  3 1  4 1");
    }
    if (i==ny*2){
      //  right top corner
      fprintf (s1,"  5  1 1  2 1  2 4  3 1  4 1");
    }
    if (i!=0 && i<ny*2){
      fprintf (s1,"  4  1 0  2 4  3 1  4 1");
    }
  }




  fprintf (s1,"\n%ld",ne);
  fprintf(stdout, "\n The number of elements  %ld", ne);
  
  nid=1;
  for (i=0;i<nx;i++){
    for (j=0;j<ny;j++){
      fprintf (s1,"\n%8ld 6  %8ld %8ld %8ld %8ld %8ld %8ld %8ld %8ld       1",
	       nid,
	       i*(ny*2+1)+i*(ny+1)+j*2+1,
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+j*2+1,
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+(j+1)*2+1,
	       i*(ny*2+1)+i*(ny+1)+(j+1)*2+1,
	       (i+1)*(ny*2+1)+i*(ny+1)+j+1,
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+j*2+2,
	       (i+1)*(ny*2+1)+i*(ny+1)+j+2,
	       i*(ny*2+1)+i*(ny+1)+j*2+2);
/* original node ordering
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+(j+1)*2+1,
	       i*(ny*2+1)+i*(ny+1)+(j+1)*2+1,
	       i*(ny*2+1)+i*(ny+1)+j*2+1,
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+j*2+1,
	       (i+1)*(ny*2+1)+i*(ny+1)+j+2,
	       i*(ny*2+1)+i*(ny+1)+j*2+2,
	       (i+1)*(ny*2+1)+i*(ny+1)+j+1,
	       (i+1)*(ny*2+1)+(i+1)*(ny+1)+j*2+2);*/
      nid++;
      if (edg)
      {

        //fprintf(s1, " 4");
        if (j == 0)
          fprintf(s1, " 3");
        else
          fprintf(s1, " 0");

        if (i == nx-1)
          fprintf(s1, " 4");
        else
          fprintf(s1, " 0");

        if (j == ny-1)
          fprintf(s1, " 1");
        else
          fprintf(s1, " 0");

        if (i == 0)
          fprintf(s1, " 2");
        else
          fprintf(s1, " 0");

	fprintf(s1, "   1");
	
      }
    }
  }

  fprintf(stdout, "\n\n");
  fprintf (s1,"\n");

  fclose (s1);
}
