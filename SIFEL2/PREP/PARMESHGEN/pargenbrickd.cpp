#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <conio.h>
#include "intools.h"
#include "stacktrace.h"

double *X, *Y, *Z;  /// x and y coordinates for each node
long Nn, Ne;    /// number of nodes, number of elements
long *Gnn;      /// global node numbers
long **Lnn, **Len;
/** Lnn[did][i] gives global node id (for X, Y and Z array) for did-th domain and i-th node
    Len[did][i] gives global element id for did-th domain and i-th element
*/
long Ndx, Ndy, Ndz;  /// number of domains in the x, y and z direction
long *Nex, *Ney, *Nez; /// number of elements in the x, y and z direction of each domain
long Nexo, Neyo, Nezo;
double *Dimdx, *Dimdy, *Dimdz; /// array with x, y and z dimension of domains
long **El;             /// node numbers for each element
char Fname[1025];
char Afname[1025];
char Tmpl[1025];

void read(FILE *in)
{
  long i, ret;
  fscanf(in, "%ld %ld %ld", &Ndx, &Ndy, &Ndz);
  Dimdx = new double[Ndx];
  Dimdy = new double[Ndy];
  Dimdz = new double[Ndz];
  Nex = new long[Ndx];
  Ney = new long[Ndy];
  Nez = new long[Ndz];
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
    ret = fscanf(in, "%lf", Dimdx+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading x dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    ret = fscanf(in, "%ld", Nex+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading number of element in x dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    fprintf(stdout, "%e %ld\n", Dimdx[i], Nex[i]);
    Nexo += Nex[i];
  }
  fprintf(stdout, "\nDomain in Y direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndy; i++)
  {
    ret = fscanf(in, "%lf", Dimdy+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading y dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    ret = fscanf(in, "%ld", Ney+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading number of element in y dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    fprintf(stdout, "%e %ld\n", Dimdy[i], Ney[i]);
    Neyo += Ney[i];
  }
  fprintf(stdout, "\nDomain in Z direction :\n");
  fprintf(stdout, "length, number of elements\n");
  for (i = 0; i < Ndz; i++)
  {
    ret = fscanf(in, "%lf", Dimdz+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading z dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    ret = fscanf(in, "%ld", Nez+i);
    if (ret != 1)
    {
      fprintf(stderr, "\nError reading number of element in z dimension of %ld domain\n", i+1);
      fprintf(stderr, "File %s, line number %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program terminated\n");
      abort();
    }
    fprintf(stdout, "%e %ld\n", Dimdz[i], Nez[i]);
    Nezo += Nez[i];
  }
  getnextln(in);
  inputln(in, Fname, 1024);
  inputln(in, Afname, 1024);
  Nn = (Nexo+1) * (Neyo+1) * (Nezo+1);
  Ne = Nexo * Neyo * Nezo;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  inputln(in, Tmpl, 1024);
}

void gen_nodes(void)
{
  long dix, djy, dkz, i, j, k, *lid, ni, ngi, ngic, did, negk, negj, negi;
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

  ni   =  0;
  ngi  =  1;
  ngic = -1;
  did  =  0;
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
      if ((dix == Ndx-1) && (i == Nex[dix]-1))
        negi = 1;
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
          if ((djy == Ndy-1) && (j == Ney[djy]-1))
            negj = 1;
          for (dkz = 0; dkz < Ndz; dkz++)
          // loop over domains in the z direction
          {
            lz = Dimdz[dkz]/Nez[dkz];
            if (dkz == Ndz-1)
              Nez[dkz]++;
            for (k = 0; k < Nez[dkz]; k++)
            // loop over elements in the z direction
            {
              if ((dkz == Ndz-1) && (k == Nez[dkz]-1))
                negk = 1;
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
                if ((j == 0) || (i == 0) || negj || negi)
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngic;
                    ngic--;
                  }
                }
                else
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngi;
                    ngi++;
                  }
                }
              }
              if ((j == 0) && (djy != 0))
              // cislovani uzlu v hranicni rovine kolme na osu y
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = dix*Ndy*Ndz + (djy-1)*Ndz + dkz;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if ((k == 0) || (i == 0) || negk || negi)
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngic;
                    ngic--;
                  }
                }
                else
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngi;
                    ngi++;
                  }
                }
              }
              if ((i == 0) && (dix != 0))
              // cislovani uzlu v hranicni rovine kolme na osu x
              {
                // lokalni cislovani posledniho sloupce v predchozi domene
                did = (dix-1)*Ndy*Ndz + djy*Ndz + dkz;
                Lnn[did][lid[did]] = ni;
                lid[did]++;
                if ((j == 0) || (k == 0) || negj || negk)
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngic;
                    ngic--;
                  }
                }
                else
                {
                  if (Gnn[ni] == 0)
                  // globalni cislovani plochy
                  {
                    Gnn[ni] = ngi;
                    ngi++;
                  }
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
            {
              Nez[dkz]--;
              negk = 0;
            }  
          }   // konec cyklu dkz
          yy += ly;
        } // konec cyklu j
        if (djy == Ndy-1)
        {
          Ney[djy]--;
          negj = 0;
        }  
      }   // konec cyklu dyj
      xx += lx;
    }  // konec cyklu i
    if (dix == Ndx-1)
    {
      Nex[dix]--;
      negi = 0;
    }  
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


void printtop(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, j, k, ii, jj, kk, id, did, l;

  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        id = 0;
        did = i * Ndy * Ndz + j * Ndz + k;
        strcpy(tmpfn, fname);
        sprintf(tmpfn+strlen(tmpfn), "%ld.top", did+1);
        out = fopen(tmpfn, "wt");
        if (out == NULL)
        {
          fprintf(stderr, "\nError - unable open topology file %s\n", tmpfn);
          fprintf(stderr, " domain number %ld\n", did+1);
          return;
        }
        fprintf(out, "%ld\n", (Nex[i]+1)*(Ney[j]+1)*(Nez[k]+1));
        // prints nodes
        for (ii = 0; ii < Nex[i]+1; ii++)
        {
          for (jj = 0; jj < Ney[j]+1; jj++)
          {
            for (kk = 0; kk < Nez[k]+1; kk++)
            {
              fprintf(out, "%e %e %e", X[Lnn[did][id]], Y[Lnn[did][id]], Z[Lnn[did][id]]);
              if (ii == 0)
              {
                fprintf(out, " 1\n");
                id++;
                continue;
              }
              if (ii == Nex[i])
              {
                fprintf(out, " 3\n");
                id++;
                continue;
              }
              fprintf(out, " 2\n");
              id++;
            }
          }
        }
        id = 0;
        fprintf(out, "%ld 13\n", Nex[i]*Ney[j]*Nez[k]); // number of elements and element type number
        // prints elements
        for (ii = 0; ii < Nex[i]; ii++)
        {
          for (jj = 0; jj < Ney[j]; jj++)
          {
            for (kk = 0; kk < Nez[k]; kk++)
            {
              fprintf(out, "8");
              for (l = 0; l < 8; l++)
                fprintf(out, " %ld", El[Len[did][id]][l]+1);
              fprintf(out, " 0\n");
              id++;
            }
          }
        }
        id = 0;
        // prints global node numbers
        for (ii = 0; ii < Nex[i]+1; ii++)
        {
          for (jj = 0; jj < Ney[j]+1; jj++)
          {
            for (kk = 0; kk < Nez[k]+1; kk++)
            {
              fprintf(out, "%ld\n", Gnn[Lnn[did][id]]);
              id++;
            }
          }
        }
        fclose(out);
      }
    }
  }
}

int printprep(char *fname, char *afname, char *tmplfn, FILE *in)
{
  char prepfn[2049];
  char outfn[1025];
  char topfn[1025];
  char *ptr;
  FILE *out, *tmplf;
  long i, j, k, id, did, l, nlc, end;
  char str[1205];
  double *fx=NULL, *fy=NULL, *fz=NULL;

  tmplf = fopen(tmplfn,"r");
  if (tmplf==NULL){
    fprintf (stderr,"\n Template preprocesor file could not be opened.");
    return(1);
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      for (k = 0; k < Ndz; k++)
      {
        id = 0;
        did = i * Ndy * Ndz + j * Ndz + k;
        strcpy(prepfn, fname);
        strcpy(outfn, afname);
        strcpy(topfn, fname);
        sprintf(prepfn+strlen(prepfn), "%ld.pr", did+1);
//        sprintf(outfn+strlen(outfn), "%ld.out", did+1);
        sprintf(topfn+strlen(topfn), "%ld.top", did+1);
        out = fopen(prepfn, "wt");
        if (out == NULL)
        {
          fprintf(stderr, "\nError - unable open preprocessor file %s\n", prepfn);
          fprintf(stderr, " domain number %ld\n", did+1);
          return(2);
        }
        fseek(tmplf, 0L, SEEK_SET);
        fprintf(out, "%s\n", topfn); // output file name
        getstring(tmplf, str, 1024);     // material dbase
        fprintf(out, "%s\n", str);
        getstring(tmplf, str, 1024);     // cross section dbase
        fprintf(out, "%s\n", str);
        getnextln(tmplf);
        //fprintf(out, "0\n0\n");          // topology format indentifiers
        // probdesc section
        end = 0;
        do
        {
          inputln(tmplf, str, 1024);
          switch (isgkw(str, ptr))
          {
            case -1 :
              fprintf(out, "%s\n", str); // normal line
              break;
            case 0 :
              fprintf(out, "%s\n", str); // string with 'begin' keyword
              fprintf(out, "%ld\n", did+1);// domain number
              getstring(tmplf, str, 1024);
              getnextln(tmplf);
              fprintf(out, "%s\n", str); // string job title
              break;
            case 1 :
              fprintf(out, "%s\n", str); // string with 'end' keyword
              end = 1;
              break;
          }
        } while (! end);
        getlong(tmplf, nlc);
        fprintf(out, "%ld\n", nlc);
        // reads nodal forces for each load case
        if (fx == NULL)
        {
          fx = new double[nlc];
          fy = new double[nlc];
          fz = new double[nlc];
          memset(fx, 0, sizeof(*fx)*nlc);
          memset(fy, 0, sizeof(*fy)*nlc);
          memset(fz, 0, sizeof(*fz)*nlc);
          for (l = 0; l < nlc; l++)
            fscanf(in, "%le %le %le", &fx[l], &fy[l], &fz[l]);
        }
        // nodal property section
        do
        {
          inputln(tmplf, str, 1024);
          if (isgkw(str, ptr) == 1) // string with 'end' keyword
            break;
          else
            fprintf(out, "%s\n", str); // other line
        } while (! feof(tmplf));
        // add boundary condition to domains on the left side
        if (i == 0)
        {
          fprintf(out, "1 bocon 3 1 0.0 2 0.0 3 0.0\n");
        }
        // add nodal forces to domains on the right side
        if (i == Ndx-1)
        {
          for (l = 0; l < nlc; l++)
            fprintf(out, "3 nload %ld %e %e %e\n", l+1, fx[l], fy[l], fz[l]);
        }
        fprintf(out, "end\n");
        // section with element property
        do
        {
          inputln(tmplf, str, 1024);
          if (isgkw(str, ptr) == 1) // string with 'end' keyword
          {
            fprintf(out, "%s\n", str); // string with 'end' keyword
            break;
          }
          else
            fprintf(out, "%s\n", str); // other line
        } while (! feof(tmplf));
	
	// section with output description
	do
	  {
	    inputln(tmplf, str, 1024);
	    if (isgkw(str, ptr) == 1) // string with 'end' keyword
	      {
		fprintf(out, "%s\n", str); // string with 'end' keyword
		break;
	      }
	    else
	      fprintf(out, "%s\n", str); // other line
	  } while (! feof(tmplf));
	
        fclose(out);
      }
    }
  }
  return(0);
}

int main (int argc,char *argv[])
{
  long i, j, k, ii, jj, kk, l, ie, did;
  FILE *in;
  FILE *out = fopen("vystup","w");

  set_prgname(argv[0]);
  fprintf (stderr,"\n\n *** PARMEFEL GENERATOR OF TOPOLOGY ***\n");
  fprintf (stderr," --------------------------------------\n");
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name\n\n", argv[0]);
    return(1);
  }
  in = fopen(argv[1],"r");
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
    fprintf(out, "%5ld %e %e %e %5ld\n", i+1, X[i], Y[i], Z[i], Gnn[i]);
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
  printtop(Fname);
  printprep(Fname, Afname, Tmpl, in);
  printf ("\n");
  fprintf (out,"\n");
  fclose(in);
  fclose(out);
  //getch();
  return 0;
}
