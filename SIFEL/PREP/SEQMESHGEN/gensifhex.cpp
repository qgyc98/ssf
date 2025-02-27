#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include "iotools.h"


double Dimdx, Dimdy, Dimdz; /// array with x, y and z dimension of domains
long Nex, Ney, Nez; /// number of elements in the x, y and z direction of each domain
char *Fname;        /// output file name
long Edg;           /// indicator for printing of edge and surface properties

long Nn, Ne;       /// number of nodes, number of elements
double *X, *Y, *Z; /// x and y coordinates for each node
long **El;         /// node numbers for each element



void gen_nodes(void)
{
  long i, j, k, ni;
  double xx, yy, zz, lx, ly, lz;

  X = new double [Nn];
  Y = new double [Nn];
  Z = new double [Nn];
  memset(X, 0, sizeof(*X)*Nn);
  memset(Y, 0, sizeof(*Y)*Nn);
  memset(Z, 0, sizeof(*Z)*Nn);

  ni = 0;
  xx = 0.0;
  lx = Dimdx/Nex;
  for (i = 0; i < Nex+1; i++)
  // loop over elements in the x direction
  {
    yy = 0.0;
    ly = Dimdy/Ney;
    for (j = 0; j < Ney+1; j++)
    // loop over elements in the y direction
    {
      zz = 0.0;
      lz = Dimdz/Nez;
      for (k = 0; k < Nez+1; k++)
      // loop over elements in the z direction
      {
        X[ni] = xx;
        Y[ni] = yy;
        Z[ni] = zz;
        ni++;
        zz += lz;
      } // konec cyklu k
      yy += ly;
    } // konec cyklu j
    xx += lx;
  }  // konec cyklu i
}



void gen_elem(void)
{
  long ii, jj, kk, ie = 0;
  long nidrb, nidlb, nidlt, nidrt;
  /*
    nidrb = local node id of right bottom node in view of xy plane
    nidlb = local node id of left bottom node in view of xy plane
    nidlt = local node id of left top node in view of xy plane
    nidrt = local node id of right top node in view of xy plane
  */
  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);

  for (ii = 0; ii < Nex; ii++)
  {
    for (jj = 0; jj < Ney; jj++)
    {
      for (kk = 0; kk < Nez; kk++)
      {
        El[ie] = new long[8];
        memset(El[ie], 0, sizeof(*El[ie])*8);
        nidlb = ii * (Ney + 1) * (Nez + 1) + jj * (Nez + 1) + kk;
        nidrb = (ii+1) * (Ney + 1) * (Nez + 1) + jj * (Nez + 1) + kk;
        nidlt = ii * (Ney + 1) * (Nez + 1) + (jj + 1) * (Nez + 1) + kk;
        nidrt = (ii+1) * (Ney + 1) * (Nez + 1) + (jj + 1) * (Nez + 1) + kk;
        El[ie][0] = nidrt+1;
        El[ie][1] = nidlt+1;
        El[ie][2] = nidlb+1;
        El[ie][3] = nidrb+1;
        El[ie][4] = nidrt;
        El[ie][5] = nidlt;
        El[ie][6] = nidlb;
        El[ie][7] = nidrb;
        ie++;
      }
    }
  }
}



void printtop(char *fname)
{
  FILE *out;
  long ii, jj, kk, id, l;
  long edgn[12], surfn[6];
  long  aux1;
  long np;
  

  out = fopen(fname, "wt");
  if (out == NULL)
  {
    fprintf(stderr, "\nError - unable open topology file %s\n", fname);
    return;
  }
  
  fprintf(out, "%ld\n", Nn);
  
  id = 0;
  // prints nodes
  for (ii = 0; ii < Nex+1; ii++)
  {
    for (jj = 0; jj < Ney+1; jj++)
    {
      for (kk = 0; kk < Nez+1; kk++)
      {               
        fprintf(out, "%8ld  %le  %le  %le   ", id+1, X[id], Y[id], Z[id]);
        aux1 = ftell(out);
        fprintf(out, "  "); // create space for later write of number of properties



        np = 0;             
	
        //
        // VERTICES
        //
        if ((ii == Nex) && (jj == Ney) && (kk == Nez))
        // front top right vertex with property 1
        {
          np++;
        }
        if ((ii == 0) && (jj == Ney) && (kk == Nez))
        // rear top right vertex with property 2
        {
          np++;
        }
        if ((ii == 0) && (jj == 0) && (kk == Nez))
        // rear top left vertex with property 3
        {
          np++;
        }
        if ((ii == Nex) && (jj == 0) && (kk == Nez))
        // front top left vertex with property 4
        {
          np++;
        }
        if ((ii == Nex) && (jj == Ney) && (kk == 0))
        // rear bottom left vertex with property 5
        {
          np++;
        }
        if ((ii == 0) && (jj == Ney) && (kk == 0))
        // rear bottom right vertex with property 6
        {
          np++;
        }
        if ((ii == 0) && (jj == 0) && (kk == 0))
        // rear bottom left vertex with property 7
        {
          np++;
        }
        if ((ii == Nex) && (jj == 0) && (kk == 0))
        // front bottom left vertex with property 8
        {
          np++;
        }
	
	if (np==0){
          np++;
	}
	
        //
        // EDGES
        //
        if ((jj == Ney) && (kk == Nez)) // top right edge with property 1
        {
          np++;
        }
        if ((ii == 0) && (kk == Nez)) // top rear edge with property 2
        {
          np++;
        }
        if ((jj == 0) && (kk == Nez)) // top left edge with property 3
        {
          np++;
        }
        if ((ii == Nex) && (kk == Nez)) // front top edge with property 4
        {
          np++;
        }
        if ((ii == Nex) && (jj == Ney)) // front right edge with property 5
        {
          np++;
        }
        if ((ii == 0) && (jj == Ney)) // rear right edge with property 6
        {
          np++;
        }
        if ((ii == 0) && (jj == 0)) // rear left edge with property 7
        {
          np++;
        }
        if ((ii == Nex) && (jj == 0)) // front left edge with property 8
        {
          np++;
        }
        if ((jj == Ney) && (kk == 0)) // bottom right edge with property 9
        {
          np++;
        }
        if ((ii == 0) && (kk == 0)) // rear bottom edge with property 10
        {
          np++;
        }
        if ((jj == 0) && (kk == 0)) // bottom left edge with property 11
        {
          np++;
        }
        if ((ii == Nex) && (kk == 0)) // front bottom edge with property 12
        {
          np++;
        }


        //
        // SURFACES
        //
        if (ii == 0)  // rear surface with property 3
        {
          np++;
        }

        if (jj == 0)  // left surface with property 4
        {
          np++;
        }

        if (kk == 0)  // bottom surface with property 6
        {
          np++;
        }

        if (ii == Nex)  // front surface with property 1
        {
          np++;
        }

        if (jj == Ney)  // rear surface with property 2 
        {
          np++;
        }

        if (kk == Nez)  // rear surface with property 5
        {
          np++;
        }


        // VOLUMES
        np++;
        //aux2 = ftell(out);

        // write number of properties before property records
        //fseek(out, aux1, SEEK_SET);
        fprintf(out,"   %ld", np);
        //fseek(out, aux2, SEEK_SET);
        //id++;




        np = 0;             
	
        //
        // VERTICES
        //
        if ((ii == Nex) && (jj == Ney) && (kk == Nez))
        // front top right vertex with property 1
        {
          fprintf(out, "  1 1");
          np++;
        }
        if ((ii == 0) && (jj == Ney) && (kk == Nez))
        // rear top right vertex with property 2
        {
          fprintf(out, "  1 2");
          np++;
        }
        if ((ii == 0) && (jj == 0) && (kk == Nez))
        // rear top left vertex with property 3
        {
          fprintf(out, "  1 3");
          np++;
        }
        if ((ii == Nex) && (jj == 0) && (kk == Nez))
        // front top left vertex with property 4
        {
          fprintf(out, "  1 4");
          np++;
        }
        if ((ii == Nex) && (jj == Ney) && (kk == 0))
        // rear bottom left vertex with property 5
        {
          fprintf(out, "  1 5");
          np++;
        }
        if ((ii == 0) && (jj == Ney) && (kk == 0))
        // rear bottom right vertex with property 6
        {
          fprintf(out, "  1 6");
          np++;
        }
        if ((ii == 0) && (jj == 0) && (kk == 0))
        // rear bottom left vertex with property 7
        {
          fprintf(out, "  1 7");
          np++;
        }
        if ((ii == Nex) && (jj == 0) && (kk == 0))
        // front bottom left vertex with property 8
        {
          fprintf(out, "  1 8");
          np++;
        }
	
	if (np==0){
          fprintf(out, "  1 0");
          np++;
	}
	
        //
        // EDGES
        //
        if ((jj == Ney) && (kk == Nez)) // top right edge with property 1
        {
          fprintf(out, "  2 1");
          np++;
        }
        if ((ii == 0) && (kk == Nez)) // top rear edge with property 2
        {
          fprintf(out, "  2 2");
          np++;
        }
        if ((jj == 0) && (kk == Nez)) // top left edge with property 3
        {
          fprintf(out, "  2 3");
          np++;
        }
        if ((ii == Nex) && (kk == Nez)) // front top edge with property 4
        {
          fprintf(out, "  2 4");
          np++;
        }
        if ((ii == Nex) && (jj == Ney)) // front right edge with property 5
        {
          fprintf(out, "  2 5");
          np++;
        }
        if ((ii == 0) && (jj == Ney)) // rear right edge with property 6
        {
          fprintf(out, "  2 6");
          np++;
        }
        if ((ii == 0) && (jj == 0)) // rear left edge with property 7
        {
          fprintf(out, "  2 7");
          np++;
        }
        if ((ii == Nex) && (jj == 0)) // front left edge with property 8
        {
          fprintf(out, "  2 8");
          np++;
        }
        if ((jj == Ney) && (kk == 0)) // bottom right edge with property 9
        {
          fprintf(out, "  2 9");
          np++;
        }
        if ((ii == 0) && (kk == 0)) // rear bottom edge with property 10
        {
          fprintf(out, "  2 10");
          np++;
        }
        if ((jj == 0) && (kk == 0)) // bottom left edge with property 11
        {
          fprintf(out, "  2 11");
          np++;
        }
        if ((ii == Nex) && (kk == 0)) // front bottom edge with property 12
        {
          fprintf(out, "  2 12");
          np++;
        }


        //
        // SURFACES
        //
        if (ii == 0)  // rear surface with property 3
        {
          fprintf(out, "  3 3");
          np++;
        }

        if (jj == 0)  // left surface with property 4
        {
          fprintf(out, "  3 4");
          np++;
        }

        if (kk == 0)  // bottom surface with property 6
        {
          fprintf(out, "  3 6");
          np++;
        }

        if (ii == Nex)  // front surface with property 1
        {
          fprintf(out, "  3 1");
          np++;
        }

        if (jj == Ney)  // rear surface with property 2 
        {
          fprintf(out, "  3 2");
          np++;
        }

        if (kk == Nez)  // rear surface with property 5
        {
          fprintf(out, "  3 5");
          np++;
        }


        // VOLUMES
        np++;
        fprintf(out, "  4 1\n");  // no numbering
        //aux2 = ftell(out);

        // write number of properties before property records
        //fseek(out, aux1, SEEK_SET);
        //fprintf(out, " %ld", np);
        //fseek(out, aux2, SEEK_SET);
        id++;
      }
    }
  }

  fprintf(out, "\n%ld\n", Ne); // number of elements and element type number
  id = 0;
  // prints elements
  for (ii = 0; ii < Nex; ii++)
  {
    for (jj = 0; jj < Ney; jj++)
    {
      for (kk = 0; kk < Nez; kk++)
      {
        fprintf(out, "%8ld 13  ", id+1);
        for (l = 0; l < 8; l++)
          fprintf(out, " %6ld", El[id][l]+1);

        memset(edgn,  0, sizeof(*edgn)*12);
        memset(surfn, 0, sizeof(*surfn)*6);
        //
        // EDGES
        //
        if ((jj == Ney-1) && (kk == Nez-1)) // top right edge with property 1
          edgn[0] = 1;

        if ((ii == 0) && (kk == Nez-1)) // top rear edge with property 2
          edgn[1] = 2;

        if ((jj == 0) && (kk == Nez-1)) // top left edge with property 3
          edgn[2] = 3;

        if ((ii == Nex-1) && (kk == Nez-1)) // front top edge with property 4
          edgn[3] = 4;

        if ((ii == Nex-1) && (jj == Ney-1)) // front right edge with property 5
          edgn[4] = 5;

        if ((ii == 0) && (jj == Ney-1)) // rear right edge with property 6
          edgn[5] = 6;

        if ((ii == 0) && (jj == 0)) // rear left edge with property 7
          edgn[6] = 7;

        if ((ii == Nex-1) && (jj == 0)) // front left edge with property 8
          edgn[7] = 8;
        
        if ((jj == Ney-1) && (kk == 0)) // bottom right edge with property 9
          edgn[8] = 9;

        if ((ii == 0) && (kk == 0)) // rear bottom edge with property 10
          edgn[9] = 10;

        if ((jj == 0) && (kk == 0)) // bottom left edge with property 11
          edgn[10] = 11;

        if ((ii == Nex-1) && (kk == 0)) // front bottom edge with property 12
          edgn[11] = 12;

        //
        // SURFACES
        //
        if (ii == 0)  // rear surface with property 3
          surfn[2]=3;

        if (jj == 0)  // left surface with property 4
          surfn[3]=4;

        if (kk == 0)  // bottom surface with property 6
          surfn[5]=6;

        if (ii == Nex-1)  // front surface with property 1
          surfn[0]=1;

        if (jj == Ney-1)  // rear surface with property 2 
          surfn[1]=2;

        if (kk == Nez-1)  // rear surface with property 5
          surfn[4]=5;

        //
        // VOLUMES
        //
        fprintf(out, "    1");  // no numbering

        if (Edg)
        {
          fprintf(out, " ");
          for (l=0; l<12; l++)
            fprintf(out, " %ld", edgn[l]);
          fprintf(out, " ");
          for (l=0; l<6; l++)
            fprintf(out, " %ld", surfn[l]);
        }
        fprintf(out, "\n");
        id++;
      }
    }
  }
  fclose(out);
}



int main (int argc,char *argv[])
{
  FILE *s1;
  long i;
  Fname = NULL;
  Nex = 0;
  Ney = 0;
  Nez = 0;


  fprintf(stdout, "\n\n\n---      GENERATOR OF BRICK ELEMENTS ON BRICK DOMAIN     ---");
  fprintf(stdout, "\n--- generates hexahedral elements on a hexahedral domain ---\n");
  if (argc < 9){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s output_file_name lx ly lz nx ny nz edg_id\n\n", argv[0]);
    return(1);
  }

  s1 = fopen(argv[1], "wt");
  if (s1==NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  fclose(s1);
  Fname = argv[1];

  if (sscanf(argv[2], "%lf", &Dimdx) != 1)
  {
    fprintf (stderr,"\n Length of domain lx has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if (sscanf(argv[5], "%ld", &Nex) != 1)
  {
    fprintf (stderr,"\n Number of elements in x direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(6);
  }

  if (sscanf(argv[3], "%lf", &Dimdy) != 1)
  {
    fprintf (stderr,"\n Length of domain ly has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(4);
  }
  if (sscanf(argv[6], "%ld", &Ney) != 1)
  {
    fprintf (stderr,"\n Number of elements in y direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }

  if (sscanf(argv[4], "%lf", &Dimdz) != 1)
  {
    fprintf (stderr,"\n Length of domain lz has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(5);
  }
  if (sscanf(argv[7], "%ld", &Nez) != 1)
  {
    fprintf (stderr,"\n Number of elements in z direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(8);
  }
  if (sscanf(argv[8], "%ld", &Edg) != 1)
  {
    fprintf(stderr, "\n Edge property indicator has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(9);
  }


  Nn = (Nex+1) * (Ney+1) * (Nez+1);
  Ne = Nex * Ney * Nez;

  fprintf(stdout, "\n Length of prism in x direction   %f", Dimdx);
  fprintf(stdout, "\n Length of prism in y direction   %f", Dimdy);
  fprintf(stdout, "\n Length of prism in z direction   %f", Dimdz);
  fprintf(stdout, "\n The number of elements in x direction   %ld", Nex);
  fprintf(stdout, "\n The number of elements in y direction   %ld", Ney);
  fprintf(stdout, "\n The number of elements in z direction   %ld", Nez);
  fprintf(stdout, "\n File name for output   %s", Fname);
  fprintf(stdout, "\n Edge numbering     %ld", Edg);

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
  
  fprintf(stdout, "\n\n The number of nodes      %ld", Nn);
  fprintf(stdout, "\n The number of elements   %ld", Ne);

  //fprintf(stdout, "\nGeneration of nodes ...");
  fflush(stdout);
  gen_nodes();
  //fprintf(stdout, " O.K.\n");
  //fprintf(stdout, "\nGeneration of elements ...");
  fflush(stdout);
  gen_elem();
  //fprintf(stdout, " O.K.\n");

  //fprintf(stdout, "\nGeneration of file with topology:\n");
  fflush(stdout);
  printtop(Fname);
  
  fprintf(stdout, "\n\n");

  delete [] X;
  delete [] Y;
  delete [] Z;
  for (i=0; i<Ne; i++)
    delete [] El[i];
  delete [] El;
  return 0;
}

