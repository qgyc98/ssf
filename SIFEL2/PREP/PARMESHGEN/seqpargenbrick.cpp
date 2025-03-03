#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iotools.h"

double *X, *Y, *Z;  /// x and y coordinates for each node
long Nn, Ne;    /// number of nodes, number of elements
long *Gnn;      /// global node numbers
long **Lnn, **Len;
/** Lnn[did][i] gives global node id (for X, Y and Z array) for did-th domain and i-th node
    Len[did][i] gives global element id for did-th domain and i-th element
*/
long Ndx, Ndy, Ndz;    /// number of domains in the x, y and z direction
long *Nex, *Ney, *Nez; /// number of elements in the x, y and z direction of each domain
long Nexo, Neyo, Nezo; /// total numbers of elements in the x, y and z direction
double *Dimdx, *Dimdy, *Dimdz; /// array with x, y and z dimension of domains
long **El;             /// node numbers for each element
char *Fname;
char *Afname;
char *Tmpl;
long Edg;            /// indicator for printing of edge and surface properties
long Volprop;        /// indicator for printing of domain id as a volume property

void read(XFILE *in)
{
  long i;
  xfscanf(in, "%ld %ld %ld", &Ndx, &Ndy, &Ndz);
  Dimdx = new double[Ndx];
  Dimdy = new double[Ndy];
  Dimdz = new double[Ndz];
  Nex = new long[Ndx];
  Ney = new long[Ndy];
  Nez = new long[Ndz];

  Fname  = new char[in->give_maxlnsize()+1];
  Afname = new char[in->give_maxlnsize()+1];
  Tmpl   = new char[in->give_maxlnsize()+1];
  memset(Fname, 0, sizeof(*Fname)*(in->give_maxlnsize()+1));
  memset(Afname, 0, sizeof(*Afname)*(in->give_maxlnsize()+1));
  memset(Tmpl, 0, sizeof(*Tmpl)*(in->give_maxlnsize()+1));

  memset(Dimdx, 0, sizeof(*Dimdx)*Ndx);
  memset(Dimdy, 0, sizeof(*Dimdy)*Ndy);
  memset(Dimdz, 0, sizeof(*Dimdz)*Ndz);
  memset(Nex, 0, sizeof(*Nex)*Ndx);
  memset(Ney, 0, sizeof(*Ney)*Ndy);
  memset(Nez, 0, sizeof(*Nez)*Ndz);
  Nexo = 0;
  Neyo = 0;
  Nezo = 0;
  fprintf(stdout, "\nDomain in X direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndx; i++)
  {
    xfscanf(in, "%lf", Dimdx+i);
    xfscanf(in, "%ld", Nex+i);
    fprintf(stdout, "%le %ld\n", Dimdx[i], Nex[i]);
    Nexo += Nex[i];
  }
  fprintf(stdout, "\nDomain in Y direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndy; i++)
  {
    xfscanf(in, "%lf", Dimdy+i);
    xfscanf(in, "%ld", Ney+i);
    fprintf(stdout, "%le %ld\n", Dimdy[i], Ney[i]);
    Neyo += Ney[i];
  }
  fprintf(stdout, "\nDomain in Z direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndz; i++)
  {
    xfscanf(in, "%lf", Dimdz+i);
    xfscanf(in, "%ld", Nez+i);
    fprintf(stdout, "%le %ld\n", Dimdz[i], Nez[i]);
    Nezo += Nez[i];
  }
  xfscanf(in, " %a", Fname);
  xfscanf(in, " %a", Afname);
  Edg = 0;
  xfscanf(in, "%ld", &Edg);
  xfscanf(in, "%ld", &Volprop);
  Nn = (Nexo+1) * (Neyo+1) * (Nezo+1);
  Ne = Nexo * Neyo * Nezo;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  fprintf(stdout, "\nEdge numbering  : %ld\n", Edg);
  fprintf(stdout, "Domain numbering  : %ld\n", Volprop);
}

void gen_nodes(void)
{
  long dix, djy, dkz, i, j, k, *lid, ni, ngi, did;
  double xx, yy, zz, lx, ly, lz;

  X = new double [Nn];
  Y = new double [Nn];
  Z = new double [Nn];
  Gnn = new long [Nn];
  Lnn = new long *[Ndx*Ndy*Ndz];
  lid = new long[Ndx*Ndy*Ndz];
  memset(X, 0, sizeof(*X)*Nn);
  memset(Y, 0, sizeof(*Y)*Nn);
  memset(Z, 0, sizeof(*Z)*Nn);
  memset(Gnn, 0, sizeof(*Gnn)*Nn);
  memset(Lnn, 0, sizeof(*Lnn)*Ndx*Ndy*Ndz);
  memset(lid, 0, sizeof(*lid)*Ndx*Ndy*Ndz);
  ni = 0;
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        Lnn[ni] = new long[(Nex[i]+1)*(Ney[j]+1)*(Nez[k]+1)];
        memset(Lnn[ni], 0, sizeof(*Lnn[ni])*(Nex[i]+1)*(Ney[j]+1)*(Nez[k]+1));
        ni++;
      }
    }
  }

  ni = 0;
  ngi = 1;
  did = 0;
  xx = 0.0;
  for (dix = 0; dix < Ndx; dix++)
  // loop over domains in the x direction
  {
    lx = Dimdx[dix]/Nex[dix];
    if (dix == Ndx-1)
      Nex[dix]++;
    for (i = 0; i < Nex[dix]; i++)
    // loop over elements in the x direction
    {
      yy = 0.0;
      for (djy = 0; djy < Ndy; djy++)
      // loop over domains in the y direction
      {
        ly = Dimdy[djy]/Ney[djy];
        if (djy == Ndy-1)
          Ney[djy]++;
        for (j = 0; j < Ney[djy]; j++)
        // loop over elements in the y direction
        {
          zz = 0.0;
          for (dkz = 0; dkz < Ndz; dkz++)
          // loop over domains in the z direction
          {
            lz = Dimdz[dkz]/Nez[dkz];
            if (dkz == Ndz-1)
              Nez[dkz]++;
            for (k = 0; k < Nez[dkz]; k++)
            // loop over elements in the z direction
            {
              X[ni] = xx;
              Y[ni] = yy;
              Z[ni] = zz;
              did = dix * Ndy * Ndz + djy * Ndz + dkz;
              // lokalni cislovani dane domeny
              Lnn[did][lid[did]] = ni;
              lid[did]++;
              if ((k == 0) && (dkz != 0))
              // cislovani uzlu v hranicni rovine kolme na osu z
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = dix*Ndy*Ndz + djy*Ndz + dkz-1;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if ((j == 0) && (djy != 0))
              // cislovani uzlu v hranicni rovine kolme na osu y
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = dix*Ndy*Ndz + (djy-1)*Ndz + dkz;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if ((i == 0) && (dix != 0))
              // cislovani uzlu v hranicni rovine kolme na osu x
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = (dix-1)*Ndy*Ndz + djy*Ndz + dkz;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if (((k == 0) && (dkz != 0)) && ((j == 0) && (djy != 0)))
              // cislovani uzlu v hranicni primce = prusecnice 2 rovin, ktere jsou kolme na osu z a y
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = dix*Ndy*Ndz + (djy-1)*Ndz + dkz-1;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if (((j == 0) && (djy != 0)) && ((i == 0) && (dix != 0)))
              // cislovani uzlu v hranicni primce = prusecnice 2 rovin, ktere jsou kolme na osu y a x
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = (dix-1)*Ndy*Ndz + (djy-1)*Ndz + dkz;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if (((i == 0) && (dix != 0)) && ((k == 0) && (dkz != 0)))
              // cislovani uzlu v hranicni primce = prusecnice 2 rovin, ktere jsou kolme na osu x a z
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = (dix-1)*Ndy*Ndz + djy*Ndz + dkz-1;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              if (((k == 0) && (dkz != 0)) && ((j == 0) && (djy != 0)) && ((i == 0) && (dix != 0)))
              // cislovani uzlu v hranicnim bode, kde se protinaji 3 prusecnice
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = (dix-1)*Ndy*Ndz + (djy-1)*Ndz + dkz-1;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if (Gnn[ni] == 0)
                // globalni cislovani hrany
                {
                  Gnn[ni] = ngi;
                  ngi++;
                }
              }
              ni++;
              zz += lz;
            } // konec cyklu k
            if (dkz == Ndz-1)
              Nez[dkz]--;
          }   // konec cyklu dkz
          yy += ly;
        } // konec cyklu j
        if (djy == Ndy-1)
          Ney[djy]--;
      }   // konec cyklu dyj
      xx += lx;
    }  // konec cyklu i
    if (dix == Ndx-1)
      Nex[dix]--;
  }    // konec cyklu dix
  delete [] lid;
}

void gen_elem(void)
{
  long i, j, k, ii, jj, kk, ie = 0, lie;
  long did, nidrb, nidlb, nidlt, nidrt;
  /*
    nidrb = local node id of right bottom node in view of xy plane
    nidlb = local node id of left bottom node in view of xy plane
    nidrt = local node
  */
  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);
  Len = new long *[Ndx*Ndy*Ndz];
  memset(Len, 0, sizeof(*Len)*Ndx*Ndy*Ndz);

  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        lie = 0;
        did = i * Ndy * Ndz + j * Ndz + k;
        Len[did] = new long[Nex[i]*Ney[j]*Nez[k]];
        memset(Len[did], 0, sizeof(*Len[did])*Nex[i]*Ney[j]*Nez[k]);
        for (ii = 0; ii < Nex[i]; ii++)
        {
          for (jj = 0; jj < Ney[j]; jj++)
          {
            for (kk = 0; kk < Nez[k]; kk++)
            {
              El[ie] = new long[8];
              memset(El[ie], 0, sizeof(*El[ie])*8);
              nidlb = ii * (Ney[j] + 1) * (Nez[k] + 1) + jj * (Nez[k] + 1) + kk;
              nidrb = (ii+1) * (Ney[j] + 1) * (Nez[k] + 1) + jj * (Nez[k] + 1) + kk;
              nidlt = ii * (Ney[j] + 1) * (Nez[k] + 1) + (jj + 1) * (Nez[k] + 1) + kk;
              nidrt = (ii+1) * (Ney[j] + 1) * (Nez[k] + 1) + (jj + 1) * (Nez[k] + 1) + kk;
              El[ie][0] = nidrt+1;
              El[ie][1] = nidlt+1;
              El[ie][2] = nidlb+1;
              El[ie][3] = nidrb+1;
              El[ie][4] = nidrt;
              El[ie][5] = nidlt;
              El[ie][6] = nidlb;
              El[ie][7] = nidrb;
              Len[did][lie] = ie;
              lie++;
              ie++;
            }
          }
        }
      }
    }
  }
}



void printtopseq(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, j, k, ii, jj, kk, id, did, l;
  long *prnn, edgn[12], surfn[6];
  long  aux1, aux2;
  long np;
  
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
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        id = 0;
        did = i * Ndy * Ndz + j * Ndz + k;
        // prints nodes
        for (ii = 0; ii < Nex[i]+1; ii++)
        {
          for (jj = 0; jj < Ney[j]+1; jj++)
          {
            for (kk = 0; kk < Nez[k]+1; kk++)
            {               
              if (prnn[Lnn[did][id]] == 1)  // given node has been already printed
              {
                id++;
                continue;
              }
              prnn[Lnn[did][id]]= 1;
              fprintf(out, "%ld %le %le %le", Lnn[did][id]+1, X[Lnn[did][id]], Y[Lnn[did][id]], Z[Lnn[did][id]]);
              aux1 = ftell(out);
              fprintf(out, "  "); // create space for later write of number of properties

              np = 0;             
              //
              // VERTICES
              //
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == Ndy-1) && (jj == Ney[j]) && (k == Ndz-1) && (kk == Nez[k]))
              // front top right vertex with property 1
              {
                fprintf(out, "  1 1");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j]) && (k == Ndz-1) && (kk == Nez[k]))
              // rear top right vertex with property 2
              {
                fprintf(out, "  1 2");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0) && (k == Ndz-1) && (kk == Nez[k]))
              // rear top left vertex with property 3
              {
                fprintf(out, "  1 3");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == 0) && (jj == 0) && (k == Ndz-1) && (kk == Nez[k]))
              // front top left vertex with property 4
              {
                fprintf(out, "  1 4");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == Ndy-1) && (jj == Ney[j]) && (k == 0) && (kk == 0))
              // rear bottom left vertex with property 5
              {
                fprintf(out, "  1 5");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j]) && (k == 0) && (kk == 0))
              // rear bottom right vertex with property 6
              {
                fprintf(out, "  1 6");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0) && (k == 0) && (kk == 0))
              // rear bottom left vertex with property 7
              {
                fprintf(out, "  1 7");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == 0) && (jj == 0) && (k == 0) && (kk == 0))
              // front bottom left vertex with property 8
              {
                fprintf(out, "  1 8");
                np++;
              }


              //
              // EDGES
              //
              if ((j == Ndy-1) && (jj == Ney[j]) && (k == Ndz-1) && (kk == Nez[k])) // top right edge with property 1
              {
                fprintf(out, "  2 1");
                np++;
              }
              if ((i == 0) && (ii == 0) && (k == Ndz-1) && (kk == Nez[k])) // top rear edge with property 2
              {
                fprintf(out, "  2 2");
                np++;
              }
              if ((j == 0) && (jj == 0) && (k == Ndz-1) && (kk == Nez[k])) // top left edge with property 3
              {
                fprintf(out, "  2 3");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (k == Ndz-1) && (kk == Nez[k])) // front top edge with property 4
              {
                fprintf(out, "  2 4");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == Ndy-1) && (jj == Ney[j])) // front right edge with property 5
              {
                fprintf(out, "  2 5");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j])) // rear right edge with property 6
              {
                fprintf(out, "  2 6");
                np++;
              }
              if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0)) // rear left edge with property 7
              {
                fprintf(out, "  2 7");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (j == 0) && (jj == 0)) // front left edge with property 8
              {
                fprintf(out, "  2 8");
                np++;
              }
              if ((j == Ndy-1) && (jj == Ney[j]) && (k == 0) && (kk == 0)) // bottom right edge with property 9
              {
                fprintf(out, "  2 9");
                np++;
              }
              if ((i == 0) && (ii == 0) && (k == 0) && (kk == 0)) // rear bottom edge with property 10
              {
                fprintf(out, "  2 10");
                np++;
              }
              if ((j == 0) && (jj == 0) && (k == 0) && (kk == 0)) // bottom left edge with property 11
              {
                fprintf(out, "  2 11");
                np++;
              }
              if ((i == Ndx-1) && (ii == Nex[i]) && (k == 0) && (kk == 0)) // front bottom edge with property 12
              {
                fprintf(out, "  2 12");
                np++;
              }


              //
              // SURFACES
              //
              if ((i == 0) && (ii == 0))  // rear surface with property 3
              {
                fprintf(out, "  3 3");
                np++;
              }

              if ((j == 0) && (jj == 0))  // left surface with property 4
              {
                fprintf(out, "  3 4");
                np++;
              }

              if ((k == 0) && (kk == 0))  // bottom surface with property 6
              {
                fprintf(out, "  3 6");
                np++;
              }

              if ((i == Ndx-1) && (ii == Nex[i]))  // front surface with property 1
              {
                fprintf(out, "  3 1");
                np++;
              }

              if ((j == Ndy-1) && (jj == Ney[j]))  // rear surface with property 2 
              {
                fprintf(out, "  3 2");
                np++;
              }

              if ((k == Ndz-1) && (kk == Nez[k]))  // rear surface with property 5
              {
                fprintf(out, "  3 5");
                np++;
              }


              // VOLUMES
              np++;
              switch (Volprop)
              {
                case 1:
                  fprintf(out, "  4 %ld\n", i+1);  // numbering domain layers in x direction
                  break;
                case 2:
                  fprintf(out, "  4 %ld\n", j+1);  // numbering domain layers in y direction
                  break;
                case 3:
                  fprintf(out, "  4 %ld\n", k+1);  // numbering domain layers in z direction
                  break;
                case 4:
                  fprintf(out, "  4 %ld\n", did+1); // domain numbering
                  break;
                default:
                  fprintf(out, "  4 1\n");  // no numbering
              }
              aux2 = ftell(out);

              // write number of properties before property records
              fseek(out, aux1, SEEK_SET);
              fprintf(out, " %ld", np);
              fseek(out, aux2, SEEK_SET);
              id++;
            }
          }
        }
      }
    }
  }
  fprintf(out, "\n%ld\n", Ne); // number of elements and element type number
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        id = 0;
        did = i * Ndy * Ndz + j * Ndz + k;
        // prints elements
        for (ii = 0; ii < Nex[i]; ii++)
        {
          for (jj = 0; jj < Ney[j]; jj++)
          {
            for (kk = 0; kk < Nez[k]; kk++)
            {
              fprintf(out, "%ld 13", Len[did][id]+1);
              for (l = 0; l < 8; l++)
                fprintf(out, " %ld", Lnn[did][El[Len[did][id]][l]]+1);

              memset(edgn,  0, sizeof(*edgn)*12);
              memset(surfn, 0, sizeof(*surfn)*6);
              //
              // EDGES
              //
              if ((j == Ndy-1) && (jj == Ney[j]-1) && (k == Ndz-1) && (kk == Nez[k]-1)) // top right edge with property 1
                edgn[0] = 1;

              if ((i == 0) && (ii == 0) && (k == Ndz-1) && (kk == Nez[k]-1)) // top rear edge with property 2
                edgn[1] = 2;

              if ((j == 0) && (jj == 0) && (k == Ndz-1) && (kk == Nez[k]-1)) // top left edge with property 3
                edgn[2] = 3;

              if ((i == Ndx-1) && (ii == Nex[i]-1) && (k == Ndz-1) && (kk == Nez[k]-1)) // front top edge with property 4
                edgn[3] = 4;

              if ((i == Ndx-1) && (ii == Nex[i]-1) && (j == Ndy-1) && (jj == Ney[j]-1)) // front right edge with property 5
                edgn[4] = 5;

              if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j]-1)) // rear right edge with property 6
                edgn[5] = 6;

              if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0)) // rear left edge with property 7
                edgn[6] = 7;

              if ((i == Ndx-1) && (ii == Nex[i]-1) && (j == 0) && (jj == 0)) // front left edge with property 8
                edgn[7] = 8;

              if ((j == Ndy-1) && (jj == Ney[j]-1) && (k == 0) && (kk == 0)) // bottom right edge with property 9
                edgn[8] = 9;

              if ((i == 0) && (ii == 0) && (k == 0) && (kk == 0)) // rear bottom edge with property 10
                edgn[9] = 10;

              if ((j == 0) && (jj == 0) && (k == 0) && (kk == 0)) // bottom left edge with property 11
                edgn[10] = 11;

              if ((i == Ndx-1) && (ii == Nex[i]-1) && (k == 0) && (kk == 0)) // front bottom edge with property 12
                edgn[11] = 12;

              //
              // SURFACES
              //
              if ((i == 0) && (ii == 0))  // rear surface with property 3
                surfn[2]=3;

              if ((j == 0) && (jj == 0))  // left surface with property 4
                surfn[3]=4;

              if ((k == 0) && (kk == 0))  // bottom surface with property 6
                surfn[5]=6;

              if ((i == Ndx-1) && (ii == Nex[i]-1))  // front surface with property 1
                surfn[0]=1;

              if ((j == Ndy-1) && (jj == Ney[j]-1))  // rear surface with property 2 
                surfn[1]=2;

              if ((k == Ndz-1) && (kk == Nez[k]-1))  // rear surface with property 5
                surfn[4]=5;

              //
              // VOLUMES
              //
              switch (Volprop)
              {
                case 1:
                  fprintf(out, "  %ld", i+1);  // numbering domain layers in x direction
                  break;
                case 2:
                  fprintf(out, "  %ld", j+1);  // numbering domain layers in y direction
                  break;
                case 3:
                  fprintf(out, "  %ld", k+1);  // numbering domain layers in z direction
                  break;
                case 4:
                  fprintf(out, "  %ld", did+1); // domain numbering
                  break;
                default:
                  fprintf(out, "  1");  // no numbering
              }

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
      }
    }
  }
  fclose(out);
}



int main (int argc,char *argv[])
{

  long i, j, k, ii, jj, kk, l, ie, did;
  XFILE *in;
  FILE *out = fopen("vystup","w");

  fprintf (stderr,"\n\n *** PARMEFEL GENERATOR OF TOPOLOGY ***\n");
  fprintf (stderr," --------------------------------------\n");
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name\n\n", argv[0]);
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  read(in);
  gen_nodes();
  gen_elem();
  fprintf(out, "%ld\n", Nn);
  for (i = 0; i < Nn; i++)
  {
    fprintf(out, "%5ld %le %le %le %5ld\n", i+1, X[i], Y[i], Z[i], Gnn[i]);
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        fprintf(out, "\nDomena : %ld,%ld,%ld\n", i+1, j+1, k+1);
        did = i*Ndy*Ndz + j*Ndz + k;
        for (ii = 0; ii < Nex[i]+1; ii++)
        {
          for (jj = 0; jj < Ney[j]+1; jj++)
          {
            for (kk = 0; kk < Nez[k]+1; kk++)
              fprintf(out, " %ld", Lnn[did][ii*(Ney[j]+1)*(Nez[k]+1)+jj*(Nez[k]+1)+kk]+1);
            fprintf(out, "\n");
          }
        }
      }
    }
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        fprintf(out, "\nDomena : %ld,%ld,%ld\n%ld\n", i+1, j+1, k+1, Nex[i]*Ney[j]*Nez[k]);
        ie = 0;
        did = i*Ndy*Ndz + j*Ndz + k;
        for (ii = 0; ii < Nex[i]; ii++)
        {
          for (jj = 0; jj < Ney[j]; jj++)
          {
            for (kk = 0; kk < Nez[k]; kk++)
            {
              fprintf(out, "%ld -", ie+1);
              for (l = 0; l < 8; l++)
                fprintf(out, " %ld", El[Len[did][ie]][l]+1);
              ie++;
              fprintf(out, "\n");
            }
          }
        }
        fprintf(out, "\n");
      }
    }
  }
  printtopseq(Fname);
  xfclose(in);
  fprintf (out,"\n");
  fclose(out);
}
