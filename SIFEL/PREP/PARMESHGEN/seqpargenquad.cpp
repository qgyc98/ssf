#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"

double *X, *Y;  /// x and y coordinates for each node
long Nn, Ne;    /// number of nodes, number of elements
long Ndofn;     /// number of degrees of freedom in the node
long *Gnn;      /// global node numbers
long **Lnn, **Len;
/** Lnn[did][i] gives global node id (for X and Y array) for did-th domain and i-th node
    Len[did][i] gives local element id for did-th domain and i-th element
*/
long Ndx, Ndy;  /// number of domains in the x and y directions
long *Nex, *Ney; /// number of elements in the x and y directions of each domain
long Nexo, Neyo;
double *Dimdx, *Dimdy; /// array with x and y dimension of domains
long **El;             /// node numbers for each element
char *Fname;
char *Afname;
char *Tmpl;
long Edg;
long Volprop;        /// indicator for printing of domain id as a volume property

void read(XFILE *in)
{
  long i;
  xfscanf(in, "%ld%ld", &Ndx, &Ndy);
  Dimdx = new double[Ndx];
  Dimdy = new double[Ndy];
  Nex = new long[Ndx];
  Ney = new long[Ndy];

  Fname  = new char[in->give_maxlnsize()+1];
  Afname = new char[in->give_maxlnsize()+1];
  Tmpl   = new char[in->give_maxlnsize()+1];
  memset(Fname, 0, sizeof(*Fname)*(in->give_maxlnsize()+1));
  memset(Fname, 0, sizeof(*Afname)*(in->give_maxlnsize()+1));
  memset(Tmpl, 0, sizeof(*Tmpl)*(in->give_maxlnsize()+1));

  memset(Dimdx, 0, sizeof(*Dimdx)*Ndx);
  memset(Dimdy, 0, sizeof(*Dimdy)*Ndy);
  memset(Nex, 0, sizeof(*Nex)*Ndx);
  memset(Ney, 0, sizeof(*Ney)*Ndy);
  Nexo = 0;
  Neyo = 0;
  fprintf(stdout, "\nDomains in X direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndx; i++)
  {
    xfscanf(in, "%lf", Dimdx+i);
    xfscanf(in, "%ld", Nex+i);
    fprintf(stdout, "%le %ld\n", Dimdx[i], Nex[i]);
    Nexo += Nex[i];
  }
  fprintf(stdout, "\nDomains in Y direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndy; i++)
  {
    xfscanf(in, "%lf", Dimdy+i);
    xfscanf(in, "%ld", Ney+i);
    fprintf(stdout, "%le %ld\n", Dimdy[i], Ney[i]);
    Neyo += Ney[i];
  }
  xfscanf(in, " %a", Fname);
  xfscanf(in, " %a", Afname);
  Edg = 0;
  xfscanf(in, "%ld", &Edg);
  xfscanf(in, "%ld", &Volprop);
  Nn = (Nexo+1) * (Neyo+1);
  Ne = Nexo * Neyo;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  fprintf(stdout, "\nEdge numbering  : %ld\n", Edg);
  fprintf(stdout, "\nVolume property numbering  : %ld\n", Volprop);
}

void gen_nodes(void)
{
  long dix, djy, i, j, *lid, ni, ngi, did;
  double xx, yy, lx, ly;

  X = new double [Nn];
  Y = new double [Nn];
  Gnn = new long [Nn];
  Lnn = new long *[Ndx*Ndy];
  lid = new long[Ndx*Ndy];
  memset(X, 0, sizeof(*X)*Nn);
  memset(Y, 0, sizeof(*Y)*Nn);
  memset(Gnn, 0, sizeof(*Gnn)*Nn);
  memset(Lnn, 0, sizeof(*Lnn)*Ndx*Ndy);
  ni = 0;
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      lid[i*Ndy+j] = 0;
      Lnn[ni] = new long[(Nex[i]+1)*(Ney[j]+1)];
      memset(Lnn[ni], 0, sizeof(*Lnn[ni])*(Nex[i]+1)*(Ney[j]+1));
      ni++;
    }
  }

  ni = 0;
  ngi = 1;
  did = 0;
  xx = 0.0;
  for (dix = 0; dix < Ndx; dix++)
  {
    lx = Dimdx[dix]/Nex[dix]; // prirustek souradnice ve smeru x
    if (dix == Ndx-1)
      Nex[dix]++;
    for (i = 0; i < Nex[dix]; i++)
    {
      yy = 0.0;
      for (djy = 0; djy < Ndy; djy++)
      {
        ly = Dimdy[djy]/Ney[djy]; // prirustek souradnice ve smeru y
        if ((djy != 0) && (Gnn[ni] == 0))
        // globalni cislovani vodorovne hrany
        {
          Gnn[ni] = ngi;
          ngi++;
        }
        for (j = 0; j < Ney[djy]; j++)
        {
          X[ni] = xx;
          Y[ni] = yy;
          did = dix*Ndy+djy;
          // lokalni cislovani dane domeny
          Lnn[did][lid[did]] = ni;
          lid[did]++;
//          printf("\nLnn : %ld,%ld (%ld) - %ld", dix, djy, did, ni);
          if ((i == 0) && (dix != 0))
          {
            // lokalni cislovani posledniho sloupce v domene vpravo od aktualni
            did = (dix-1)*Ndy+djy;
            Lnn[did][lid[did]] = ni;
//            printf("\nLnn : %ld,%ld (%ld) - %ld", dix, djy, did, ni);
            lid[did]++;
            if (Gnn[ni] == 0)
            // globalni cislovani svisle hrany
            {
              Gnn[ni] = ngi;
              ngi++;
            }
          }
          ni++;
          yy += ly;
        }
        did = dix*Ndy+djy;
        Lnn[did][lid[did]] = ni;
//        printf("\nLnn : %ld,%ld (%ld) - %ld", dix, djy, did, ni);
        lid[did]++;
        if ((i == 0) && (dix != 0))
        {
          did = (dix-1)*Ndy+djy;
          Lnn[did][lid[did]] = ni;
//          printf("\nLnn : %ld,%ld (%ld) - %ld", dix, djy, did, ni);
          lid[did]++;
        }
      }
      // posledni uzel v sloupci
      X[ni] = xx;
      Y[ni] = yy;
      if ((i == 0) && (dix != 0) && (Gnn[ni] == 0))
      // globalni cislovani svisle hrany
      {
        Gnn[ni] = ngi;
        ngi++;
      }
      ni++;
      xx += lx;
    }
    if (dix == Ndx-1)
      Nex[dix]--;
  }
  delete [] lid;
}

void gen_elem(void)
{
  long i, j, ii, jj, ie = 0, lie;
  long did, nidr, nidl;
  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);
  Len = new long *[Ndx*Ndy];
  memset(Len, 0, sizeof(*Len)*Ndx*Ndy);

  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      lie = 0;
      did = i * Ndy + j;
      Len[did] = new long[Nex[i]*Ney[j]];
      memset(Len[did], 0, sizeof(*Len[did])*Nex[i]*Ney[j]);
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < Ney[j]; jj++)
        {
          El[ie] = new long[4];
          memset(El[ie], 0, sizeof(*El[ie])*4);
          nidl = (ii) * (Ney[j] + 1) + jj;
          nidr = (ii+1) * (Ney[j] + 1) + jj;
/*          El[ie][0] = Lnn[did][nidl];
          El[ie][1] = Lnn[did][nidr];
          El[ie][2] = Lnn[did][nidr+1];
          El[ie][3] = Lnn[did][nidl+1];*/
          El[ie][0] = nidl;
          El[ie][1] = nidr;
          El[ie][2] = nidr+1;
          El[ie][3] = nidl+1;
          Len[did][lie] = ie;
          lie++;
          ie++;
        }
      }
    }
  }
}



void printtopseq(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, j, ii, jj, id, did, k, volid;
  long *prnn;
  
  prnn = new long[Nn];
  memset(prnn, 0, sizeof(*prnn)*Nn);
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
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      id = 0;
      did = i * Ndy + j;
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          if (prnn[Lnn[did][id]] == 1)  // give node has been already printed
          {
            id++;
            continue;
          }
          prnn[Lnn[did][id]]= 1;
          fprintf(out, "%ld %le %le 0.0 ", Lnn[did][id]+1, X[Lnn[did][id]], Y[Lnn[did][id]]);
          switch (Volprop)
          {
            case 1:
              volid = i+1;   // numbering domain layers in x direction
              break;
            case 2:
              volid = j+1;   // numbering domain layers in y direction
              break;
            case 4:
              volid = did+1; // domain numbering
              break;
            default:
              volid = 1;     // no numbering
          }
          if ((i == 0) && (ii == 0)) // first column of domains and first column of nodes
          {
            if ((j == 0) && (jj == 0))
            {
              // bottom left corner has vertex property 1, edge property 1 and 4, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 1  2 1  2 4  3 %ld  4 %ld\n", did+1, volid);
              id++;
              continue;
            }
            if ((j == Ndy-1) && (jj == Ney[j]))
            {
              // top left corner has vertex property 4, edge property 4 and 3, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 4  2 3  2 4  3 %ld  4 %ld\n", did+1, volid);
              id++;
              continue;
            }
            // nodes on the left edge have edge property 4, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 4  3 %ld  4 %ld\n", did+1, volid);
            id++;
            continue;
          }
          if ((i == Ndx-1) && (ii == Nex[i])) // last column of domains and last column of nodes
          {
            if ((j == 0) && (jj == 0))
            {
              // bottom right corner has vertex property 2, edge property 1 and 2, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 2  2 1  2 2  3 %ld  4 %ld\n", did+1, volid);
              id++;
              continue;
            }
            if ((j == Ndy-1) && (jj == Ney[j]))
            {
              // top right corner has vertex property 3, edge property 2 and 3, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 3  2 2  2 3  3 %ld  4 %ld\n", did+1, volid);
              id++;
              continue;
            }
            // nodes on the right edge have edge property 2, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 2  3 %ld  4 %ld\n", did+1, volid);
            id++;
            continue;
          }
          if ((j == 0) && (jj == 0))
          {
            // internal nodes on the bottom edge have edge property 1, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 1  3 %ld  4 %ld\n", did+1, volid);
            id++;
            continue;
          }
          if ((j == Ndy-1) && (jj == Ney[j]))
          {
            // internal nodes on the top edge have edge property 3, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 3  3 %ld  4 %ld\n", did+1, volid);
            id++;
            continue;
          }
          // internal nodes have surface property by domain id and volume property 1 
          fprintf(out, "2  3 %ld  4 %ld\n", did+1, volid);
          id++;
        }
      }
    }
  }
  delete [] prnn;

  // PRINTING OF ELEMENTS
  fprintf(out, "\n%ld\n", Ne);
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      id = 0;
      did = i * Ndy + j;
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < Ney[j]; jj++)
        {
          fprintf(out, "%ld 5", Len[did][id]+1);
          for (k = 0; k < 4; k++)
            fprintf(out, " %ld", Lnn[did][El[Len[did][id]][k]]+1);

          switch (Volprop)
          {
            case 1:
              volid = i+1;   // numbering domain layers in x direction
              break;
            case 2:
              volid = j+1;   // numbering domain layers in y direction
              break;
            case 4:
              volid = did+1; // domain numbering
              break;
            default:
              volid = 1;     // no numbering
          }
          // volume property
          fprintf(out, "  %ld", volid);


          if (Edg)
          {
            if ((i == 0) && (ii == 0)) // first column of domains and first column of elements
            {
              if ((j == 0) && (jj == 0))
              {
                // bottom left corner element has edge property 1 and 4, surface property 
                // by domain id and volume property 1
                fprintf(out, "  1 0 0 4  %ld\n", did+1);
                id++;
                continue;
              }
              if ((j == Ndy-1) && (jj == Ney[j]-1))
              {
                // top left corner element has edge property 4 and 3, surface property 
                // by domain id and volume property 1
                fprintf(out, "  0 0 3 4  %ld\n", did+1);
                id++;
                continue;
              }
              // element on the left edge have edge property 4, surface property by 
              // domain id and volume property 1
              fprintf(out, "  0 0 0 4  %ld\n", did+1);
              id++;
              continue;
            }
            if ((i == Ndx-1) && (ii == Nex[i]-1)) // last column of domains and last column of nodes
            {
              if ((j == 0) && (jj == 0))
              {
                // bottom right corner element has edge property 1 and 2, surface property 
                // by domain id and volume property 1
                fprintf(out, "  1 2 0 0  %ld\n", did+1);
                id++;
                continue;
              }
              if ((j == Ndy-1) && (jj == Ney[j]-1))
              {
                // top right corner element has edge property 2 and 3, surface property 
                // by domain id and volume property 1
                fprintf(out, "  0 2 3 0  %ld\n", did+1);
                id++;
                continue;
              }
              // elements on the right edge have edge property 2, surface property by 
              // domain id and volume property 1
              fprintf(out, "  0 2 0 0  %ld\n", did+1);
              id++;
              continue;
            }
            if ((j == 0) && (jj == 0))
            {
              // internal elements on the bottom edge have edge property 1, surface property by 
              // domain id and volume property 1
              fprintf(out, "  1 0 0 0  %ld\n", did+1);
              id++;
              continue;
            }
            if ((j == Ndy-1) && (jj == Ney[j]-1))
            {
              // internal elements on the top edge have edge property 3, surface property by 
              // domain id and volume property 1
              fprintf(out, "  0 0 3 0  %ld\n", did+1);
              id++;
              continue;
            }
            // internal elements have 0 edge properties surface property by domain id and volume property 1 
            fprintf(out, "  0 0 0 0  %ld\n", did+1);
            id++;
          }
        }
      }
    }
  }
  fclose(out);
}



int main (int argc,char *argv[])
{
  XFILE *in;
  //  FILE *out = fopen("vystup","w");

  fprintf (stderr,"\n\n *** PARMEFEL GENERATOR OF TOPOLOGY ***\n");
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
//   fprintf(out, "%ld\n", Nn);
//   for (i = 0; i < Nn; i++)
//   {
//     fprintf(out, "%5ld %le %le %5ld\n", i+1, X[i], Y[i], Gnn[i]);
//   }
//   for (i = 0; i < Ndx; i++)
//   {
//     for (j = 0; j < Ndy; j++)
//     {
//       fprintf(out, "\nDomena : %ld,%ld\n", i+1, j+1);
//       for (ii = 0; ii < Nex[i]+1; ii++)
//       {
//         for (jj = 0; jj < Ney[j]+1; jj++)
//         {
//           fprintf(out, " %ld", Lnn[i*Ndy+j][ii*(Ney[j]+1)+jj]+1);
//         }
//         fprintf(out, "\n");
//       }
//     }
//   }
//   for (i = 0; i < Ndx; i++)
//   {
//     for (j = 0; j < Ndy; j++)
//     {
//       fprintf(out, "\nDomena : %ld,%ld\n%ld\n", i+1, j+1, Nex[i]*Ney[j]);
//       ie = 0;
//       did = i*Ndy+j;
//       for (ii = 0; ii < Nex[i]; ii++)
//       {
//         for (jj = 0; jj < Ney[j]; jj++)
//         {
//           fprintf(out, "%ld -", ie+1);
//           for (k = 0; k < 4; k++)
//             fprintf(out, " %ld", El[Len[did][ie]][k]+1);
//           ie++;
//           fprintf(out, "\n");
//         }
//       }
//       fprintf(out, "\n");
//     }
//   }
  printtopseq(Fname);
  xfclose(in);
  //fprintf (out,"\n");
  //fclose(out);
}



