#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main (int argc,char *argv[])
/*
  28.11.1996
  22.7.1998
  14.1.2002 tk
  4.3.2007 JK
*/
{
  long i,j,k,ac,nn,nnp,ne,neh,nex,ney,nez,a,b;
  long *bc,*node,*nodeh;
  double dx,dy,dz,lx,ly,lz,xx,yy,zz;
  double *x,*y,*z;
  FILE *s1;

  fprintf(stdout, "\n\n\n---        GENERATOR OF TETRAHEDRAL ELEMENTS ON BRICK DOMAIN           ---\n");
  fprintf(stdout, "--- generates tetrahedral elements on a prism (brick) domain ---\n");
  if (argc < 8){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s output_file_name lx ly lz nx ny nz\n\n", argv[0]);
    return(1);
  }
  s1 = fopen(argv[1], "wt");
  if (s1==NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

  if (sscanf(argv[2], "%lf", &lx) != 1)
  {
    fprintf (stderr,"\n Length of prism lx has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if (sscanf(argv[3], "%lf", &ly) != 1)
  {
    fprintf (stderr,"\n Length of prism ly has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(4);
  }
  if (sscanf(argv[4], "%lf", &lz) != 1)
  {
    fprintf (stderr,"\n Length of prism lz has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(5);
  }
  if (sscanf(argv[5], "%ld", &nex) != 1)
  {
    fprintf (stderr,"\n Number of elements in x direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(6);
  }
  if (sscanf(argv[6], "%ld", &ney) != 1)
  {
    fprintf (stderr,"\n Number of elements in y direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }
  if (sscanf(argv[7], "%ld", &nez) != 1)
  {
    fprintf (stderr,"\n Number of elements in y direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(8);
  }

  fprintf(stdout, "\n Length of prism domain in x direction : %f", lx);
  fprintf(stdout, "\n Length of prism domain in y direction : %f", ly);
  fprintf(stdout, "\n Length of prism domain in z direction : %f", lz);
  fprintf(stdout, "\n Number of elements  in x direction : %ld", nex);
  fprintf(stdout, "\n Number of elements  in y direction : %ld", ney);
  fprintf(stdout, "\n Number of elements  in z direction : %ld", nez);

  /*  size of elements in each direction of coordinate axis  */
  dx=lx/nex;  dy=ly/ney;  dz=lz/nez;

  fprintf(stdout, "\n Size of elements in x direction : %e", dx);
  fprintf(stdout, "\n Size of elements in y direction : %e", dy);
  fprintf(stdout, "\n Size of elements in z direction : %e", dz);
  fprintf(stdout, "\n Nodes on surface x=0.0 of the brick are denoted by property 0");
  fprintf(stdout, "\n Nodes on surface z=lz of the brick are denoted by property 1");
  fprintf(stdout, "\n Nodes on surface x=lx of the brick are denoted by property 2");
  fprintf(stdout, "\n Nodes on surface z=0.0 of the brick are denoted property 3");
  fprintf(stdout, "\n All remaining nodes are denoted by property 4");
  fprintf(stdout, "\n All elements are denoted by property 0");

  /*  number of elements  */
  ne = nex*ney*nez;
  neh = 6*ne;//pocet ctyrstenu
  /*  number of nodes  */
  nn = (nex+1)*(ney+1)*(nez+1);
  /*  number of nodes in one cross-section  */
  nnp = (ney+1)*(nez+1);

  fprintf(stdout, "\n Number of nodes    : %ld", nn);
  fprintf(s1, "%ld", nn);

  x    = new double[nn]; //pole x souradnic uzlu
  y    = new double[nn]; //pole y souradnic uzlu
  z    = new double[nn]; //pole z souradnic uzlu
  bc   = new long[nn]; // pole pmocne
  node = new long[8*ne];//pole pro krychle
  nodeh = new long[4*neh];//pole pro ctyrsteny
  memset (x, 0, nn*sizeof(*x));
  memset (y, 0, nn*sizeof(*y));
  memset (z, 0, nn*sizeof(*z));
  memset (bc, 0, nn*sizeof(*bc));
  memset (node, 0, 8*ne*sizeof(*node));
  memset (node, 0, 4*ne*sizeof(*nodeh));
  /*  Generation of nodes  *///pomocne promene ac,xx,yy,zz
  ac=0;
  xx=0.0;  yy=0.0;
  for (j=0;j<=ney;j++){
    zz=0.0;
    for (k=0;k<=nez;k++){
      x[ac]=xx;
      y[ac]=yy;
      z[ac]=zz;
      bc[ac]=0;
      ac++;
      zz+=dz;
    }
    yy+=dy;
  }
  xx+=dx;
  for (i=1;i<=nex;i++){
    yy=0.0;
    for (j=0;j<=ney;j++){
      zz=0.0;
      for (k=0;k<=nez;k++){
	x[ac]=xx;
	y[ac]=yy;
	z[ac]=zz;
        if (i < nex)
        {
          if (k ==  0)
          {
            bc[ac] = 3;
            ac++;
            zz+=dz;
            continue;
          }
          if (k == nez)
          {
            bc[ac] = 1;
            ac++;
            zz+=dz;
            continue;
          }
          bc[ac] = 4;
        }
        else
          bc[ac] = 2;
	ac++;
	zz+=dz;
      }
      yy+=dy;
    }
    xx+=dx;
  }

  /*  Generation element nodes  */ //generovani uzlu elementu
  ac=0;
  for (i=0;i<nex;i++){
    for (j=0;j<ney;j++){
      for (k=0;k<nez;k++){
	node[ac+0]=(i+1)*nnp + (j+1)*(nez+1) + k+1+1;
	node[ac+1]=i*nnp + (j+1)*(nez+1) + k+1+1;
	node[ac+2]=i*nnp + j*(nez+1) + k+1+1;
	node[ac+3]=(i+1)*nnp + j*(nez+1) + k+1+1;
	node[ac+4]=(i+1)*nnp + (j+1)*(nez+1) + k+1;
	node[ac+5]=i*nnp + (j+1)*(nez+1) + k+1;
	node[ac+6]=i*nnp + j*(nez+1) + k+1;
	node[ac+7]=(i+1)*nnp + j*(nez+1) + k+1;
	ac+=8;
      }
    }
  }
// deleni krychle na  6 ctyrstenu
//ulozeni do pole nodeh - cele je to spatne   
  ac=0;
  for (i=0;i<ne;i++){
    b=i*24;
    for(j=0;j<6;j++){
      if(j==0){
	//prvni element
	a=(j*4)+b;
	nodeh[0+a]=node[ac];
	nodeh[1+a]=node[ac+1];
	nodeh[2+a]=node[ac+3];
	nodeh[3+a]=node[ac+7];
      }
      if(j==1){
	//druhy element
	a=(j*4)+b;
	nodeh[0+a]=node[ac];
	nodeh[1+a]=node[ac+1];
	nodeh[2+a]=node[ac+4];
	nodeh[3+a]=node[ac+7];
      }
      if(j==2){
	//treti element
	a=(j*4)+b;
	nodeh[0+a]=node[ac+1];
	nodeh[1+a]=node[ac+2];
	nodeh[2+a]=node[ac+3];
	nodeh[3+a]=node[ac+7];
      } 
      if(j==3){
	//ctvtvrty element
	a=(j*4)+b;
	nodeh[0+a]=node[ac+1];
	nodeh[1+a]=node[ac+2];
	nodeh[2+a]=node[ac+6];
	nodeh[3+a]=node[ac+7];
      }
      if(j==4){
	//paty element
	a=(j*4)+b;
	nodeh[0+a]=node[ac+1];
	nodeh[1+a]=node[ac+5];
	nodeh[2+a]=node[ac+4];
	nodeh[3+a]=node[ac+7];
      }
      if(j==5){
	//sesty element
 	a=(j*4)+b;
	nodeh[0+a]=node[ac+1];
	nodeh[1+a]=node[ac+6];
	nodeh[2+a]=node[ac+5];
	nodeh[3+a]=node[ac+7];
      }
    }
    ac+=8;
  }

  /*
  ac=0;
  fprintf (s1,"\nkrychle\n");
  for (i=0;i<ne;i++){
    fprintf (s1, "\n8 %ld %ld %ld %ld %ld %ld %ld %ld 0",
	     node[ac+0], node[ac+1], node[ac+2], node[ac+3],
	     node[ac+4], node[ac+5], node[ac+6], node[ac+7]);
    ac+=8;
  }
  */

  //fprintf (s1,"\n\n");

  /*  Print nodal coordinates  */
  ac=0;
  for (i=0;i<nn;i++)
    fprintf (s1, "\n%15.10le %15.10le %15.10le %ld", x[i], y[i], z[i], bc[i]);
  
  /*  Print hexaherdon elements  */
  fprintf(stdout, "\n Number of hexahedral elements : %ld\n\n", ne*6);
  fprintf (s1, "\n %ld 7",ne*6);
  
  ac=0;
  for (i=0;i<ne*6;i++){
    fprintf (s1, "\n 4  %ld %ld %ld %ld  0",
	     nodeh[ac+0], nodeh[ac+1], nodeh[ac+2], nodeh[ac+3]);
    ac+=4;
  }
  fprintf (s1,"\n");
  
  delete [] x;  delete [] y; delete [] z;
  delete [] bc; delete [] node;delete [] nodeh;
  fclose (s1);
}
