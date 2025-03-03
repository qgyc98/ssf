#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"



double *X, *Y, *Z;       /// x, y and z coordinates for each node
double Cosx, Cosy, Cosz; /// directional cosines to global X, Y and Z axes of generated elements
long Nn, Ne;  /// number of nodes, number of elements
long Ndofn;   /// number of degrees of freedom in the node
long Nd;      /// number of domains in the given direction
long *Ned;    /// number of elements in the given direction of each domain
double *Dimd; /// array with dimension of domains
long **El;    /// node numbers for each element
char *Fname;  /// topology file name for output
long Edg;     /// edg property inidicator (1=print edge property on elements, 0=print just volume property id
long Volprop; /// indicator for printing of domain id as a volume property (1 means volprop=domain_id, 0 means volprop=1)



void read(XFILE *in)
{
  long i;
  xfscanf(in, "%ld", &Nd);
  Dimd = new double[Nd];
  Ned = new long[Nd];

  Fname  = new char[in->give_maxlnsize()+1];
  memset(Fname, 0, sizeof(*Fname)*(in->give_maxlnsize()+1));

  memset(Dimd, 0, sizeof(*Dimd)*Nd);
  memset(Ned, 0, sizeof(*Ned)*Nd);
  Ne = 0;
  xfscanf(in, "%le %le %le", &Cosx, &Cosy, &Cosz);
  fprintf(stdout, "\nDomains in given direction : cos_x=%le, cos_y= %le, cos_z=%le\n", Cosx, Cosy, Cosz);
  fprintf(stdout, "domain length; number of elements\n");
  for (i = 0; i < Nd; i++)
  {
    xfscanf(in, "%lf", Dimd+i);
    xfscanf(in, "%ld", Ned+i);
    fprintf(stdout, "%le; %ld\n", Dimd[i], Ned[i]);
    Ne += Ned[i];
  }
  xfscanf(in, " %a", Fname);
  Edg = 0;
  xfscanf(in, "%ld", &Edg);
  xfscanf(in, "%ld", &Volprop);
  Nn = (Ne+1);
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  fprintf(stdout, "\nEdge numbering  : %ld\n", Edg);
  fprintf(stdout, "\nVolume property numbering  : %ld\n", Volprop);
}



void gen_nodes(void)
{
  long i, ngi, di;
  double xx, yy, zz, l;

  X = new double [Nn];
  Y = new double [Nn];
  Z = new double [Nn];
  memset(X, 0, sizeof(*X)*Nn);
  memset(Y, 0, sizeof(*Y)*Nn);
  memset(Z, 0, sizeof(*Z)*Nn);

  xx = 0.0;
  yy = 0.0;
  zz = 0.0;
  ngi = 0;
  for (di = 0; di < Nd; di++)
  {
    l = Dimd[di]/Ned[di]; // length increment in the given direction
    if (di == Nd-1)
      Ned[di]++;
    for (i = 0; i < Ned[di]; i++)
    {
      X[ngi] = xx;
      Y[ngi] = yy;
      Z[ngi] = zz;
      xx += l*Cosx;
      yy += l*Cosy;
      zz += l*Cosz;
      ngi++;
    }
  }
}



void gen_elem(void)
{
  long i;
  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);

  for (i = 0; i < Ne; i++)
  {
    El[i] = new long[2];
    memset(El[i], 0, sizeof(*El[i])*2);
    El[i][0] = i;
    El[i][1] = i+1;
  }
}



void printtopseq(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, di, li, volid;

  strcpy(tmpfn, fname);
  sprintf(tmpfn+strlen(tmpfn), ".top");
  out = fopen(tmpfn, "wt");
  if (out == NULL)
  {
    fprintf(stderr, "\nError - unable open topology file %s\n", tmpfn);
    return;
  }
  fprintf(out, "%ld\n", Nn);
  // PRINTING OF NODES
  di = 0;
  li = 0;
  for (i = 0; i < Nn; i++)
  {
    switch (Volprop)
    {
      case 1:
        volid = di+1;   // numbering domain layers in x direction
        break;
      default:        
        volid = 1;     // no numbering
    }

    if (Volprop == 0)
    {
      if (i==0) // the first node of the problem
      {
        fprintf(out, "%ld  %le %le %le   2  1 2  4 1\n", i+1, X[i], Y[i], Z[i]);
        continue;
      }
      if (i == Nn-1) // the last node of the problem
      {
        fprintf(out, "%ld  %le %le %le   2  1 2  4 1\n", i+1, X[i], Y[i], Z[i]);
        continue;
      }
      // ordinary node
      fprintf(out, "%ld  %le %le %le   1  4 1\n", i+1, X[i], Y[i], Z[i]);
    }
    else
    {  // volume property accoring to domain id

      if (i == 0) // the first node of the problem
      {
        fprintf(out, "%ld  %le %le %le   2  1 1  4 %ld\n", i+1, X[i], Y[i], Z[i], volid);
        li++;
        continue;
      }

      if ((li == Ned[di]) && (di < Nd-1)) // the interface node between two domains
      {
        fprintf(out, "%ld  %le %le %le   2  4 %ld  4 %ld\n", i+1, X[i], Y[i], Z[i], volid, volid+1);
        di++;
        li = 1;
        continue;
      }

      if (i == Nn-1) // the last node of the problem
      {
        fprintf(out, "%ld  %le %le %le   2  1 2  4 %ld\n", i+1, X[i], Y[i], Z[i], volid);
        di++;
        continue;
      }

      // ordinary node
      fprintf(out, "%ld  %le %le %le   1  4 %ld\n", i+1, X[i], Y[i], Z[i], volid);
      li++;
    }
  }

  // PRINTING OF ELEMENTS
  fprintf(out, "\n%ld\n", Ne);
  li = 0;
  di = 0;
  for (i = 0; i < Ne; i++)
  {
    switch (Volprop)
    {
      case 1:
        volid = di+1;   // numbering domain layers in x direction
        break;
      default:        
        volid = 1;     // no numbering
    }
    if (Edg)
      fprintf(out, "%ld  1  %ld %ld  %ld  %ld  %ld\n", i+1, El[i][0]+1, El[i][1]+1, volid, di+1, di+1);
    else
      fprintf(out, "%ld  1  %ld %ld  %ld\n", i+1, El[i][0]+1, El[i][1]+1, volid);
    li++;
    if (li == Ned[di])
    {
      di++;
      li = 0;
    }
  }
  fclose(out);
}



int main (int argc,char *argv[])
{
  XFILE *in;
  //  FILE *out = fopen("vystup","w");

  fprintf (stderr,"\n\n *** MEFEL GENERATOR OF 1D TOPOLOGY ***\n");
  fprintf (stderr," --------------------------------------\n");
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name\n\n", argv[0]);
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified or cannot be opened.");
    return(2);
  }
  read(in);
  gen_nodes();
  gen_elem();
  printtopseq(Fname);
  xfclose(in);
}



