#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main (int argc,char *argv[])
/*
   14.1.2002 tk
   4.3.2007 JK
*/
{
  long i,j,k,nid,ner,nel,net,nhseg,nrseg,nprrez,nn,ne,nx,ny,nz,typulohy,nhrez,nrrez;
  double r,l,t,x,rr,dr,alpha,alp,dalp,lx,ly,lz,y,z,dx,dy,dz;
  FILE *s1;

  fprintf(stdout, "\n---      GENERATOR OF HEXAHEDRAL ELEMENTS ON PRISM OR TUBE DOMAIN           ---");
  fprintf(stdout, "\n--- generates hexahedral quadratic elements on prism or tube segment domain ---\n");
  if (argc < 9){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use for prism : %s output_file_name 1 lx ly lz nx ny nz\n", argv[0]);
    fprintf (stderr," Use for tube  : %s output_file_name 2 r l t ner nel net alpha\n\n", argv[0]);
    return(1);
  }
  s1 = fopen(argv[1], "wt");
  if (s1==NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

  if (sscanf(argv[2], "%ld", &typulohy) != 1)
  {
    fprintf (stderr,"\n Task type has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if ((typulohy != 1) && (typulohy != 2))
  {
    fprintf (stderr,"\n Wrong task type has been specified.");
    fprintf (stderr,"\n Only 1 or 2 is accepted not %ld.", typulohy);
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if (typulohy == 1)
  {
    if (sscanf(argv[3], "%lf", &lx) != 1)
    {
      fprintf (stderr,"\n Length of prism lx has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(4);
    }
    if (sscanf(argv[4], "%lf", &ly) != 1)
    {
      fprintf (stderr,"\n Width of prism ly has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(5);
    }
    if (sscanf(argv[5], "%lf", &lz) != 1)
    {
      fprintf (stderr,"\n Heigth of prism lz has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(6);
    }
    if (sscanf(argv[6], "%ld", &nx) != 1)
    {
      fprintf (stderr,"\n Number of elements in x direction has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(7);
    }
    if (sscanf(argv[7], "%ld", &ny) != 1)
    {
      fprintf (stderr,"\n Number of elements in y direction has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(8);
    }
    if (sscanf(argv[8], "%ld", &nz) != 1)
    {
      fprintf (stderr,"\n Number of elements in y direction has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(9);
    }
    fprintf(stdout, "\n Length of the domain in x direction  %f", lx);
    fprintf(stdout, "\n Length of the domain in y direction  %f", ly);
    fprintf(stdout, "\n Length of the domain in z direction  %f", lz);
    fprintf(stdout, "\n The number of elements in x direction  %ld", nx);
    fprintf(stdout, "\n The number of elements in y direction  %ld", ny);
    fprintf(stdout, "\n The number of elements in z direction  %ld", nz);
  }
  if (typulohy == 2)
  {
    if (argc < 10){
      fprintf (stderr,"\n Wrong number of command line parameters.");
      fprintf (stderr,"\n Use for prism : gensifhexq output_file_name 1 lx ly lz nx ny nz\n");
      fprintf (stderr," Use for tube  : gensifhexq output_file_name 2 r l t ner nel net alpha\n\n");
      return(1);
    }
    if (sscanf(argv[3], "%lf", &r) != 1)
    {
      fprintf (stderr,"\n Radius of tube r has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(4);
    }
    if (sscanf(argv[4], "%lf", &l) != 1)
    {
      fprintf (stderr,"\n Length of tube l has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(5);
    }
    if (sscanf(argv[5], "%lf", &t) != 1)
    {
      fprintf (stderr,"\n Thickness of tube t has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(6);
    }
    if (sscanf(argv[6], "%ld", &ner) != 1)
    {
      fprintf (stderr,"\n Number of elements on segment has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(7);
    }
    if (sscanf(argv[7], "%ld", &nel) != 1)
    {
      fprintf (stderr,"\n Number of elements in l direction has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(8);
    }
    if (sscanf(argv[8], "%ld", &net) != 1)
    {
      fprintf (stderr,"\n Number of elements in t direction has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(9);
    }
    if (sscanf(argv[9], "%le", &alpha) != 1)
    {
      fprintf (stderr,"\n Segment angle alpha has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(10);
    }
    if ((alpha < 0.0) || (alpha > 360.0))
    {
      fprintf (stderr,"\n Segment angle alpha has not valid value.");
      fprintf (stderr,"\n Aplha = %f.", alpha);
      fprintf (stderr,"\n Aplha should be in range <0.0 deg; 360.0 deg>.");
      fprintf (stderr,"\n Try it again!\n\n");
      return(10);
    }
    fprintf(stdout, "\n Radius of tube      %f", r);
    fprintf(stdout, "\n Length of tube      %f", l);
    fprintf(stdout, "\n Thickness of tube   %f", t);
    fprintf(stdout, "\n The number of elements on segment       %ld", ner);
    fprintf(stdout, "\n The number of elements in l direction   %ld", nel);
    fprintf(stdout, "\n The number of elements in t direction   %ld", net);
    fprintf(stdout, "\n Segment angle  %f", alpha);
  }




  /***********************************************************************/
  /***********************************************************************/
  /*                       Prism generation                              */
  /***********************************************************************/
  /***********************************************************************/
  if (typulohy==1)
  {

    dx=lx/nx;  dy=ly/ny;  dz=lz/nz;
    nn=(nx+1)*((ny+1)*(2*nz+1)+ny*(nz+1))+nx*(ny+1)*(nz+1);
    ne=nx*ny*nz;


    
    
    fprintf (stdout,"\n\n");
    fprintf (stdout,"\n              N2          E1             N1");
    fprintf (stdout,"\n            ___________________________ ");
    fprintf (stdout,"\n           /|                         /|");
    fprintf (stdout,"\n          / |                        / |");
    fprintf (stdout,"\n      E2 /  |     ^S5               /  |");
    fprintf (stdout,"\n        /   |     |             E4 /   |");
    fprintf (stdout,"\n   z|  /    |            /        /    |");
    fprintf (stdout,"\n    | /     | E6        / S2     /     | E5");
    fprintf (stdout,"\n    |/      |     E3            /      |");
    fprintf (stdout,"\n N3 /_______|__________________/ N4    |");
    fprintf (stdout,"\n    |       |  /y              |       |");
    fprintf (stdout,"\n    |  <--S3| /                |    S1-->");
    fprintf (stdout,"\n    |       |/ N6       E9     |       | N5");
    fprintf (stdout,"\n    |       |__________________|_______|");
    fprintf (stdout,"\n    |      /                   |      /");
    fprintf (stdout,"\n  E7|     /                  E8|     /");
    fprintf (stdout,"\n    |    /             |S6     |    /");
    fprintf (stdout,"\n    |   /E10 /         |       |   / E12");
    fprintf (stdout,"\n    |  /    /S4                |  /");
    fprintf (stdout,"\n    | /                        | /");
    fprintf (stdout,"\n N7 |__________________________|/ N8_____ x");
    fprintf (stdout,"\n                 E11");
    fprintf (stdout,"\n");
    
    
    fprintf(stdout, "\n The number of nodes     %ld", nn);
    fprintf(stdout, "\n The number of elements  %ld\n\n", ne);
    fprintf(s1, "%ld", nn);

    /*  Nodal coordinates  */

    //  edge 7
    nid=1;
    x=0.0;  y=0.0;  z=0.0;
    for (i=0;i<2*nz+1;i++){
      if (i==0){
	//  corner 7
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 7  2 7  2 10  2 11  3 3  3 4  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (i==2*nz){
	//  corner 3
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 3  2 2  2 3  2 7  3 3  3 4  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (i!=0 &&i!=2*nz){
	//  nodes on the edge 7
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 7  3 3  3 4  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
    }
    
    y+=dy/2.0;
    //  surface 3
    for (i=0;i<ny;i++){
      z=0.0;
      for (j=0;j<nz+1;j++){
	if (j==0){
	  //  node on edge 10
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 10  3 3  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on edge 2
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 2  3 3  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  nodes on surface 3
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 3  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }
      y+=dy/2.0;
      if (i!=ny-1){
	z=0.0;
	for (j=0;j<2*nz+1;j++){
	  if (j==0){
	    //  node on edge 10
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 10  3 3  3 6  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j==2*nz){
	    //  node on edge 2
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 2  3 3  3 5  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j!=0 && j!=2*nz){
	    //  nodes on surface 3
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 3  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	}
	y+=dy/2.0;
      }
    }
    
    //  edge 6
    z=0.0;
    for (j=0;j<2*nz+1;j++){
      if (j==0){
	//  corner 6
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 6  2 6  2 9  2 10  3 2  3 3  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (j==2*nz){
	//  corner 2
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 2  2 1  2 2  2 6  3 2  3 3  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (j!=0 && j!=2*nz){
	//  nodes on edge 6
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 6  3 2  3 3  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
    }
    

    
    x+=dx/2.0;
    for (k=0;k<nx-1;k++){
      y=0.0;
      z=0.0;
      for (j=0;j<=nz;j++){
	if (j==0){
	  //  node on edge 11
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 11  3 4  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on edge 3
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 3  3 4  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  nodes on surface 4
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 4  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }
      y+=dy;
      for (i=0;i<ny-1;i++){
	z=0.0;
	for (j=0;j<=nz;j++){
	  if (j==0){
	    //  node on surface 6
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 6  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	  if (j==nz){
	    //  node on surafce 5
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 5  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	  if (j!=0 && j!=nz){
	    //  node inside
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   2  1 0  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	}
	y+=dy;
      }
      z=0.0;
      for (j=0;j<=nz;j++){
	if (j==0){
	  //  node on edge 9
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 9  3 2  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on edge 1
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 1  3 2  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  nodes on surface 2
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 2  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }

      x+=dx/2.0;
      y=0.0;  z=0.0;
      for (i=0;i<2*nz+1;i++){
	if (j==0){
	  //  node on edge 11
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 11  3 4  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
	if (j==2*nz){
	  //  node on edge 3
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 3  3 4  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
	if (j!=0 && j!=2*nz){
	  //  nodes on surface 4
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 4  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
      }
      y+=dy/2.0;
      for (i=0;i<ny-1;i++){
	z=0.0;
	for (j=0;j<nz+1;j++){
	  if (j==0){
	    //  node on surface 6
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 6  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	  if (j==nz){
	    //  node on surafce 5
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 5  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	  if (j!=0 && j!=nz){
	    //  node inside
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   2  1 0  4 1",nid,x,y,z);
	    nid++;
	    z+=dz;
	  }
	}
	y+=dy/2.0;
	z=0.0;
	for (j=0;j<2*nz+1;j++){
	  if (j==0){
	    //  node on surface 6
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 6  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j==2*nz){
	    //  node on surafce 5
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 5  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j!=0 && j!=2*nz){
	    //  node inside
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   2  1 0  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	}
      }
      
      y+=dy/2.0;
      z=0.0;
      for (j=0;j<nz+1;j++){
	if (j==0){
	  //  node on surface 6
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on surafce 5
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  node inside
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   2  1 0  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }
      
      y+=dy/2.0;
      z=0.0;
      for (j=0;j<=2*nz;j++){
	if (j==0){
	  //  node on edge 9
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 9  3 2  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
	if (j==2*nz){
	  //  node on edge 1
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 1  3 2  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
	if (j!=0 && j!=2*nz){
	  //  nodes on surface 2
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 2  4 1",nid,x,y,z);
	  nid++;
	  z+=dz/2.0;
	}
      }
      x+=dx/2.0;
    }
    
    y=0.0;  z=0.0;
    for (j=0;j<=nz;j++){
      if (j==0){
	//  node on edge 11
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 11  3 4  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j==nz){
	//  node on edge 3
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 3  3 4  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j!=0 && j!=nz){
	//  nodes on surface 4
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 4  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
    }
    y+=dy;
    for (i=1;i<ny;i++){
      z=0.0;
      for (j=0;j<=nz;j++){
	if (j==0){
	  //  node on surface 6
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on surafce 5
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  node inside
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   2  1 0  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }
      y+=dy;
    }
    z=0.0;
    for (j=0;j<=nz;j++){
      if (j==0){
	//  node on edge 9
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 9  3 2  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j==nz){
	//  node on edge 1
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 1  3 2  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j!=0 && j!=nz){
	//  nodes on surface 2
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 2  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
    }
    
    x+=dx/2.0;  y=0.0;
    z=0.0;
    for (i=0;i<2*nz+1;i++){
      if (i==0){
	//  corner 8
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 8  2 8  2 11  2 12  3 1  3 4  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (i==2*nz){
	//  corner 4
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 4  2 3  2 4  2 8  3 1  3 3  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (i!=0 &&i!=2*nz){
	//  nodes on the edge 8
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 8  3 1  3 4  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
    }
    y+=dy/2.0;
    //  surface 1
    for (i=0;i<ny-1;i++){
      z=0.0;
      for (j=0;j<nz+1;j++){
	if (j==0){
	  //  node on edge 12
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 12  3 1  3 6  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j==nz){
	  //  node on edge 4
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 4  3 1  3 5  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
	if (j!=0 && j!=nz){
	  //  nodes on surface 1
	  fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 1  4 1",nid,x,y,z);
	  nid++;
	  z+=dz;
	}
      }
      y+=dy/2.0;
      if (i!=ny-1){
	z=0.0;
	for (j=0;j<2*nz+1;j++){
	  if (j==0){
	    //  node on edge 12
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 12  3 1  3 6  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j==2*nz){
	    //  node on edge 4
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 4  3 1  3 5  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	  if (j!=0 && j!=2*nz){
	    //  nodes on surface 1
	    fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 1  4 1",nid,x,y,z);
	    nid++;
	    z+=dz/2.0;
	  }
	}
	y+=dy/2.0;
      }
    }
    
    z=0.0;
    for (j=0;j<nz+1;j++){
      if (j==0){
	//  node on edge 12
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 12  3 1  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j==nz){
	//  node on edge 4
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 4  3 1  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
      if (j!=0 && j!=nz){
	//  nodes on surface 1
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   3  1 0  3 1  4 1",nid,x,y,z);
	nid++;
	z+=dz;
      }
    }
    
    
    
    //  edge 5
    y+=dy/2.0;
    z=0.0;
    for (j=0;j<2*nz+1;j++){
      if (j==0){
	//  corner 5
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 5  2 5  2 9  2 12  3 1  3 2  3 6  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (j==2*nz){
	//  corner 1
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   8  1 1  2 1  2 4  2 5  3 1  3 2  3 5  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
      if (j!=0 && j!=2*nz){
	//  nodes on edge 5
	fprintf (s1,"\n%8ld  %15.10le  %15.10le  %15.10le   5  1 0  2 5  3 1  3 2  4 1",nid,x,y,z);
	nid++;
	z+=dz/2.0;
      }
    }
 
    /*********************/
    /*  Elements         */
    /*********************/
    nhrez=(ny+1)*(nz*2+1)+ny*(nz+1);
    nrrez=(ny+1)*(nz+1);
    nprrez=nhrez+nrrez;

    fprintf(s1, "\n%ld", ne);
    
    nid=1;
    for (i=0;i<nx;i++){
      for (j=0;j<ny;j++){
	for (k=0;k<nz;k++){
	  fprintf (s1, "\n%8ld  14   %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld   0",
		   nid,(i+1)*nprrez+(j+1)*(nz*2+1 + nz+1) + (k+1)*2 + 1,
		   i*nprrez+(j+1)*(nz*2+1 + nz+1) + (k+1)*2 + 1,
		   i*nprrez+j*(nz*2+1 + nz+1) + (k+1)*2 + 1,
		   (i+1)*nprrez+j*(nz*2+1 + nz+1) + (k+1)*2 + 1,

		   (i+1)*nprrez+(j+1)*(nz*2+1 + nz+1) + k*2 + 1,
		   i*nprrez+(j+1)*(nz*2+1 + nz+1) + k*2 + 1,
		   i*nprrez+j*(nz*2+1 + nz+1) + k*2 + 1,
		   (i+1)*nprrez+j*(nz*2+1 + nz+1) + k*2 + 1,

		   i*nprrez+nhrez+(j+1)*(nz+1)+k+2,
		   i*nprrez+(j+1)*(2*nz+1)+j*(nz+1)+k+2,
		   i*nprrez+nhrez+j*(nz+1)+k+2,
		   (i+1)*nprrez+(j+1)*(2*nz+1)+j*(nz+1)+k+2,

		   (i+1)*nprrez+(j+1)*(nz*2+1 + nz+1) + (k+1)*2,
		   i*nprrez+(j+1)*(nz*2+1 + nz+1) + (k+1)*2,
		   i*nprrez+j*(nz*2+1 + nz+1) + (k+1)*2,
		   (i+1)*nprrez+j*(nz*2+1 + nz+1) + (k+1)*2,

		   i*nprrez+nhrez+(j+1)*(nz+1)+k+1,
		   i*nprrez+(j+1)*(2*nz+1)+j*(nz+1)+k+1,
		   i*nprrez+nhrez+j*(nz+1)+k+1,
		   (i+1)*nprrez+(j+1)*(2*nz+1)+j*(nz+1)+k+1);
	  nid++;
	}
      }
    }
  }


  /*************************************************************************/
  /*************************************************************************/
  /*                     Tube segment generation                           */
  /*************************************************************************/
  /*************************************************************************/
  if (typulohy==2){

    ne = ner*nel*net;
    nn = (nel+1)*((ner+1)*(net*2+1)+ner*(net+1))+nel*(ner+1)*(net+1);


    fprintf(stdout, "\n Number of nodes    : %ld", nn);
    fprintf(s1, "%ld", nn);

    /****************************/
    /*  Print nodal coordinates */
    /****************************/

    dx = l/nel;  dr = t/net;  dalp = alpha/ner*3.14159265358979/180.0;

    x=0.0;  rr=r;
    for (k=0;k<2*net+1;k++){
      fprintf (s1, "\n %15.10le %15.10le %15.10le 0", x, rr, 0.0);
      rr+=dr/2.0;
    }
    alp=dalp/2.0;
    for (j=0;j<ner;j++){
      rr=r;
      for (k=0;k<net+1;k++){
	fprintf (s1, "\n %15.10le %15.10le %15.10le 1", x, rr*cos(alp), rr*sin(alp));
	rr+=dr;
      }
      alp+=dalp/2.0;  rr=r;
      for (k=0;k<2*net+1;k++){
	fprintf (s1, "\n %15.10le %15.10le %15.10le 1", x, rr*cos(alp), rr*sin(alp));
	rr+=dr/2.0;
      }
      alp+=dalp/2.0;
    }


    x+=dx/2.0;
    for (i=0;i<nel;i++){
      rr=r;
      for (k=0;k<net+1;k++){
	fprintf (s1, "\n %15.10le %15.10le %15.10le 0", x, rr, 0.0);
	rr+=dr;
      }
      alp=dalp;
      for (j=0;j<ner;j++){
	rr=r;
	for (k=0;k<net+1;k++){
	  fprintf (s1, "\n %15.10le %15.10le %15.10le 1", x, rr*cos(alp), rr*sin(alp));
	  rr+=dr;
	}
	alp+=dalp;
      }

      x+=dx/2.0;  rr=r;
      for (k=0;k<2*net+1;k++){
	fprintf (s1, "\n %15.10le %15.10le %15.10le 0", x, rr, 0.0);
	rr+=dr/2.0;
      }
      alp=dalp/2.0;
      for (j=0;j<ner;j++){
	rr=r;
	for (k=0;k<net+1;k++){
	  fprintf (s1, "\n %15.10le %15.10le %15.10le 1", x, rr*cos(alp), rr*sin(alp));
	  rr+=dr;
	}
	alp+=dalp/2.0;  rr=r;
	for (k=0;k<2*net+1;k++){
	  fprintf (s1, "\n %15.10le %15.10le %15.10le 1", x, rr*cos(alp), rr*sin(alp));
	  rr+=dr/2.0;
	}
	alp+=dalp/2.0;
      }
      x+=dx/2.0;
    }



    /********************************/
    /*  Print element nodes         */
    /********************************/
    nhseg=(ner+1)*(net*2+1)+ner*(net+1);
    nrseg=(ner+1)*(net+1);
    nprrez=nhseg+nrseg;

    fprintf(stdout, "\n Number of elements : %ld\n\n", ne);
    fprintf(s1, "\n%ld 14\n", ne);

    for (i=0;i<nel;i++)
    {
      for (j=0;j<ner;j++)
      {
	for (k=0;k<net;k++)
	{

	  fprintf (s1, "\n20 %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld 0",
		   (i+1)*nprrez + (j+1)*(2*net+1 + net+1) + (k+1)*2 + 1,
		       i*nprrez + (j+1)*(2*net+1 + net+1) + (k+1)*2 + 1,
		       i*nprrez + (j+1)*(2*net+1 + net+1) +     k*2 + 1,
		   (i+1)*nprrez + (j+1)*(2*net+1 + net+1) +     k*2 + 1,

		   (i+1)*nprrez + j*(2*net+1 + net+1) + (k+1)*2 + 1,
		       i*nprrez + j*(2*net+1 + net+1) + (k+1)*2 + 1,
		       i*nprrez + j*(2*net+1 + net+1) +     k*2 + 1,
		   (i+1)*nprrez + j*(2*net+1 + net+1) +     k*2 + 1,

		       i*nprrez + nhseg + (j+1)*(net+1) + k+1 + 1,
		       i*nprrez + (j+1)*(2*net+1 + net+1) + (k+1)*2,
		       i*nprrez + nhseg + (j+1)*(net+1) + k+1,
		   (i+1)*nprrez + (j+1)*(2*net+1 + net+1) + (k+1)*2,

		   (i+1)*nprrez + j*(2*net+1 + net+1) + net*2+1 + k+1 +1,
		       i*nprrez + j*(2*net+1 + net+1) + net*2+1 + k+1 +1,
		       i*nprrez + j*(2*net+1 + net+1) + net*2+1 + k+1,
		   (i+1)*nprrez + j*(2*net+1 + net+1) + net*2+1 + k+1,

		       i*nprrez + nhseg + j*(net+1) + k+1 + 1,
		       i*nprrez + j*(2*net+1 + net+1) + (k+1)*2,
		       i*nprrez + nhseg + j*(net+1) + k+1,
		   (i+1)*nprrez + j*(2*net+1 + net+1) + (k+1)*2);
	}
      }
    }
  }
  fprintf (s1,"\n");
  fclose (s1);
}
