#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "galias.h"
//#include "iotools.h"
//#include "prepalias.h"

double *X, *Y;  /// x and y coordinates for each node
long Nn, Ne;    /// number of nodes, number of elements
long Ndofn;     /// number of degrees of freedom in the node
long Ndx, Ndy;  /// number of domains in the x and y direction
long Nex, Ney; /// number of elements in the x and y direction of each domain
double Dimdx, Dimdy; /// array with x and y dimension of domains
long **El;             /// node numbers for each element
char *Fname;
long Edg;



void gen_nodes(void)
{
  long i, j, ni;
  double xx, yy, lx, ly;

  X = new double [Nn];
  Y = new double [Nn];
  memset(X, 0, sizeof(*X)*Nn);
  memset(Y, 0, sizeof(*Y)*Nn);

  ni = 0;
  xx = 0.0;
  lx = Dimdx/Nex;
  ly = Dimdy/Ney;
  for (i = 0; i < Nex+1; i++)
  // loop over elements in the x direction
  {
    yy = 0.0;
    for (j = 0; j < Ney+1; j++)
    // loop over elements in the y direction
    {
      X[ni] = xx;
      Y[ni] = yy;
      ni++;
      yy += ly;
    } // konec cyklu j
    xx += lx;
  }  // konec cyklu i
}



void gen_elem(void)
{
  long ii, jj, ie = 0;
  long nidr, nidl;
  // nidr = local node id of right bottom nod in view of xy plane
  // nidl = local node id of left bottom nod in view of xy plane

  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);

  for (ii = 0; ii < Nex; ii++)
  {
    for (jj = 0; jj < Ney; jj++)
    {
      El[ie] = new long[4];
      memset(El[ie], 0, sizeof(*El[ie])*4);
      nidl = (ii) * (Ney + 1) + jj;
      nidr = (ii+1) * (Ney + 1) + jj;
      El[ie][0] = nidl;
      El[ie][1] = nidr;
      El[ie][2] = nidr+1;
      El[ie][3] = nidl+1;
      ie++;
    }
  }
}



void printtopseq(char *fname)
{
  FILE *out;
  long ii, jj, id, k;
  long edgn[4];
  
  out = fopen(fname, "wt");
  if (out == NULL)
  {
    fprintf(stderr, "\nError - unable open topology file %s\n", fname);
    return;
  }
  fprintf(out, "%ld\n", Nn);

  // PRINTING OF NODES
  id = 0;
  for (ii = 0; ii < Nex+1; ii++)
  {
    for (jj = 0; jj < Ney+1; jj++)
    {
      fprintf(out, "%6ld % 14.11le % 14.11le 0.0     ", id+1, X[id], Y[id]);
      if (ii == 0) // first column of domains and first column of nodes
      {
        if (jj == 0)
        {
          // bottom left corner has vertex property 3, edge property 2 and 3, surface property 1
          // and volume property 1
          fprintf(out, "5  1 3  2 2  2 3  3 1  4 1\n");
          id++;
          continue;
        }
        if (jj == Ney)
        {
          // top left corner has vertex property 2, edge property 1 and 2, surface property 1
          // and volume property 1
          fprintf(out, "5  1 2  2 1  2 2  3 1  4 1\n");
          id++;
          continue;
        }
        // nodes on the left edge have edge property 2, surface property 1 
        // and volume property 1
        fprintf(out, "4  1 0  2 2  3 1  4 1\n");
        id++;
        continue;
      }
      if (ii == Nex) // last column of domains and last column of nodes
      {
        if (jj == 0)
        {
          // bottom right corner has vertex property 4, edge property 3 and 4, surface property 1
          // and volume property 1
          fprintf(out, "5  1 4  2 3  2 4  3 1  4 1\n");
          id++;
          continue;
        }
        if (jj == Ney)
        {
          // top right corner has vertex property 1, edge property 1 and 4, surface property 1
          // and volume property 1
          fprintf(out, "5  1 1  2 1  2 4  3 1  4 1\n");
          id++;
          continue;
        }
        // nodes on the right edge have edge property 4, surface property 1 
        // and volume property 1
        fprintf(out, "4  1 0  2 4  3 1  4 1\n");
        id++;
        continue;
      }
      if (jj == 0)
      {
        // internal nodes on the bottom edge have edge property 3, surface property 1 
        // and volume property 1
        fprintf(out, "4  1 0  2 3  3 1  4 1\n");
        id++;
        continue;
      }
      if (jj == Ney)
      {
        // internal nodes on the top edge have edge property 1, surface property 1 
        // and volume property 1
        fprintf(out, "4  1 0  2 1  3 1  4 1\n");
        id++;
        continue;
      }
      // internal nodes have vertex property 0 surface property 1 and volume property 1
      fprintf(out, "3  1 0  3 1  4 1\n");
      id++;
    }
  }

  // PRINTING OF ELEMENTS
  fprintf(out, "\n%ld\n", Ne);
  id = 0;
  for (ii = 0; ii < Nex; ii++)
  {
    for (jj = 0; jj < Ney; jj++)
    {
      fprintf(out, "%6ld 5", id+1);
      for (k = 0; k < 4; k++)
        fprintf(out, " %6ld", El[id][k]+1);
      fprintf(out, "     1"); // element volume property

      if (Edg)
      {
        memset(edgn,  0, sizeof(*edgn)*4);
        //
        // EDGES
        //
        if (jj == 0) // bottom edge with property 3
          edgn[0] = 3;

        if (ii == Nex-1) // right edge with property 4
          edgn[1] = 4;

        if (jj == Ney-1) // top edge with property 1
          edgn[2] = 1;

        if (ii == 0) // left edge with property 2
          edgn[3] = 2;

        fprintf(out, " ");
        // edge properties
        for (k=0; k<4; k++)
          fprintf(out, " %ld", edgn[k]);
        // surface property
        fprintf(out, "  1");
      }
      fprintf(out, "\n");
      id++;
    }
  }
  fclose(out);
}



int main (int argc,char *argv[])
{
  FILE *out;

  fprintf(stdout, "\n\n\n---          GENERATOR OF RECTANGULAR ELEMENTS ON RECTANGULAR DOMAIN       ---");
  fprintf(stdout, "\n--- generates rectangular elements with four nodes on a rectangular domain ---\n");
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
  fclose(out);
  Fname = argv[1];

  if (sscanf(argv[2], "%lf", &Dimdx) != 1)
  {
    fprintf (stderr,"\n Length of domain lx has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  if (sscanf(argv[4], "%ld", &Nex) != 1)
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
  if (sscanf(argv[5], "%ld", &Ney) != 1)
  {
    fprintf (stderr,"\n Number of elements in y direction has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }
  if (sscanf(argv[6], "%ld", &Edg) != 1)
  {
    fprintf (stderr,"\n Indicator of edge numbers has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }

  Nn = (Nex+1) * (Ney+1);
  Ne = Nex * Ney;

  fprintf(stdout, "\n Length of rectangle in x direction    %f", Dimdx);
  fprintf(stdout, "\n Length of rectangle in y direction    %f", Dimdy);
  fprintf(stdout, "\n The number of elements in x direction   %ld", Nex);
  fprintf(stdout, "\n The number of elements in y direction   %ld", Ney);
  fprintf(stdout, "\n File name for output  %s", Fname);
  fprintf(stdout, "\n Edge numbering    %ld", Edg);
  gen_nodes();
  gen_elem();
  printtopseq(Fname);


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
  
  fprintf(stdout, "\n\n The number of nodes      %ld", Nn);
  fprintf(stdout, "\n The number of elements   %ld", Ne);
  fprintf (stdout,"\n\n");
  return 0;
}


