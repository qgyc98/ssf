#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "intools.h"

double *X, *Y;  /// x and y coordinates for each node
long Nn, Ne;    /// number of nodes, number of elements
long Ndofn;     /// number of degrees of freedom in the node
long *Gnn;      /// global node numbers
long **Lnn, **Len;
/** Lnn[did][i] gives global node id (for X and Y array) for did-th domain and i-th node
    Len[did][i] gives local element id for did-th domain and i-th element 
*/
long Ndx, Ndy;  /// number of domains in the x and y direction
long *Nex, *Ney; /// number of elements in the x and y direction of each domain
long Nexo, Neyo;
double *Dimdx, *Dimdy; /// array with x and y dimension of domains
long **El;             /// node numbers for each element
char Fname[1025];
char Afname[1025];
char Tmpl[1025];
long Edg;

void read(FILE *in)
{
  long i, ret;
  fscanf(in, "%ld %ld", &Ndx, &Ndy);
  Dimdx = new double[Ndx];
  Dimdy = new double[Ndy];
  Nex = new long[Ndx];
  Ney = new long[Ndy];
  memset(Dimdx, 0, sizeof(*Dimdx)*Ndx);
  memset(Dimdy, 0, sizeof(*Dimdy)*Ndy);
  memset(Nex, 0, sizeof(*Nex)*Ndx);
  memset(Ney, 0, sizeof(*Ney)*Ndy);
  Nexo = 0;
  Neyo = 0;
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
    fprintf(stdout, "%le %ld\n", Dimdx[i], Nex[i]);
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
    fprintf(stdout, "%le %ld\n", Dimdy[i], Ney[i]);
    Neyo += Ney[i];
  }
  getnextln(in);
  inputln(in, Fname, 1024);
  inputln(in, Afname, 1024);
  Nn = (Nexo+1) * (Neyo+1);
  Ne = Nexo * Neyo;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  inputln(in, Tmpl, 1024);
  //Edg = 0;
  //fscanf(in, "%ld", &Edg);
}

void gen_nodes(void)
{
  long dix, djy, i, j, *lid, ni, ngi, ngci, did;
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
  ngci = -1;
  did = 0;
  xx = 0.0;
  for (dix = 0; dix < Ndx; dix++)
  {
    lx = Dimdx[dix]/Nex[dix];
    if (dix == Ndx-1)
      Nex[dix]++;
    for (i = 0; i < Nex[dix]; i++)
    {
      yy = 0.0;
      for (djy = 0; djy < Ndy; djy++)
      {
        ly = Dimdy[djy]/Ney[djy];
        if ((i == 0) && (Gnn[ni] == 0))
        // cislovani rohu je zaporne
        {
  	      Gnn[ni] = ngci;
		  ngci--;
	    }
		if ((dix == (Ndx-1)) && (i == (Nex[dix]-1)))
        // cislovani rohu je zaporne - posledni sloupec
        {
  	      Gnn[ni] = ngci;
		  ngci--;
	    }
        if ((djy != 0) && (Gnn[ni] == 0))
        // globalni cislovani vodorovne hrany
        {
	      // cislovani hran mezi rohy je kladne
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
              if (j == 0)
              // cislovani rohu je zaporne
		      {
		        Gnn[ni] = ngci;
			    ngci--;
		      }
		      else
	          // cislovani hran mezi rohy je kladne
		      {
		        Gnn[ni] = ngi;
                ngi++;
		      }
            }
          }
          ni++;
          yy += ly;
        } // konec cyklu j
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
      }  // konec cyklu djy
      // posledni uzel v sloupci
      X[ni] = xx;
      Y[ni] = yy;
      if ((i == 0) && (Gnn[ni] == 0))
      // globalni cislovani svisle hrany
      {
        Gnn[ni] = ngci;
        ngci--;
      }
      if ((dix == (Ndx-1)) && (i == (Nex[dix]-1)))
      // cislovani rohu je zaporne - posledni sloupec
      {
  	    Gnn[ni] = ngci;
	    ngci--;
	  }
      ni++;
      xx += lx;
    } // konec cyklu i
    if (dix == Ndx-1)
      Nex[dix]--;
  } // konec cyklu dix
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


void printtop(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, j, ii, jj, id, did, k;

  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      id = 0;
      did = i * Ndy + j;
      strcpy(tmpfn, fname);
      sprintf(tmpfn+strlen(tmpfn), "%ld.top", did+1);
      out = fopen(tmpfn, "wt");
      if (out == NULL)
      {
        fprintf(stderr, "\nError - unable open topology file %s\n", tmpfn);
        fprintf(stderr, " domain number %ld\n", did+1);
        return;
      }
      fprintf(out, "%ld\n", (Nex[i]+1)*(Ney[j]+1));
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          fprintf(out, "%le %le 0.0", X[Lnn[did][id]], Y[Lnn[did][id]]);
          if (ii == 0)
          {
            fprintf(out, " 0\n");
            id++;
            continue;
          }
          if (ii == Nex[i])
          {
            fprintf(out, " 2\n");
            id++;
            continue;
          }
          if (jj == 0)
          {
            fprintf(out, " 3\n");
            id++;
            continue;
          }
          if (jj == Ney[j])
          {
            fprintf(out, " 1\n");
            id++;
            continue;
          }
          fprintf(out, " 4\n");
          id++;
        }
      }
      id = 0;
      fprintf(out, "%ld 5\n", Nex[i]*Ney[j]);
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < Ney[j]; jj++)
        {
          fprintf(out, "4");
          for (k = 0; k < 4; k++)
            fprintf(out, " %ld", El[Len[did][id]][k]+1);
          fprintf(out, " 0");
          if (Edg)
          {
            fprintf(out, " 4");
            if (jj == 0)
              fprintf(out, " 4");
            else
              fprintf(out, " 0");

            if (ii == Nex[i]-1)
              fprintf(out, " 3");
            else
              fprintf(out, " 0");

            if (jj == Ney[j]-1)
              fprintf(out, " 2");
            else
              fprintf(out, " 0");

            if (ii == 0)
              fprintf(out, " 1");
            else
              fprintf(out, " 0");
          }
          fprintf(out, "\n");
          id++;
        }
      }
      id = 0;
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          fprintf(out, "%ld\n", Gnn[Lnn[did][id]]);
          id++;
        }
      }
      fclose(out);
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
  long i, j, id, did, k, kk, nlc, end, fmt;
  char str[1205];
  double **f=NULL;

  tmplf = fopen(tmplfn,"r");
  if (tmplf==NULL){
    fprintf (stderr,"\n Template preprocesor file could not be opened.");
    return(1);
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      id = 0;
      did = i * Ndy + j;
      strcpy(prepfn, fname);
      strcpy(outfn, afname);
      strcpy(topfn, fname);
      sprintf(prepfn+strlen(prepfn), "%ld.pr", did+1);
//      sprintf(outfn+strlen(outfn), "%ld.out", did+1);
      sprintf(topfn+strlen(topfn), "%ld.top", did+1);
      out = fopen(prepfn, "wt");
      if (out == NULL)
      {
        fprintf(stderr, "\nError - unable open preprocessor file %s\n", prepfn);
        fprintf(stderr, " domain number %ld\n", did+1);
        return(2);
      }
      fseek(tmplf, 0L, SEEK_SET);

      // topology file name
      fprintf(out, "%s\n", topfn);
      
      // material dbase
      getstring(tmplf, str, 1024);
      fprintf(out, "%s\n", str);
      
      // cross section dbase
      getstring(tmplf, str, 1024);
      fprintf(out, "%s\n", str);

      getnextln(tmplf);
      getlong(tmplf, fmt);
      if (fmt!=0)
	fprintf (stderr,"\n mesh format must be 0, see GEFEL/galias.h, (file %s, line %d)",__FILE__,__LINE__);
      getlong(tmplf, Edg);
      fprintf(out, "0\n%ld\n", Edg);    //  mesh format (this is sifel format, see GEFEL/galias.h), edges

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
            //fprintf(out, "%s\n", outfn);// output file name
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
      if (f == NULL)
      {
        f = new double*[nlc];
        memset(f, 0, sizeof(*f)*nlc);
        if (fscanf(in, "%ld", &Ndofn) != 1)
        {
          fprintf(stderr, "\n\n Error reading Ndofn from the input file\n");
          return (1);
        }

        for (k = 0; k < nlc; k++)
        {
          f[k] = new double [Ndofn];
          memset (f[k], 0, sizeof(*f[k])*Ndofn);
          for (kk = 0; kk < Ndofn; kk++)
          {
            if (fscanf(in, "%le", f[k]+kk) != 1)
            {
              fprintf(stderr, "\n\n Error reading load, direction %ld\n", kk+1);
              return(2);
            }
          }
        }
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
        fprintf(out, "0 bocon %ld", Ndofn);
        for (kk = 1; kk <= Ndofn; kk++)
          fprintf(out, " %ld 0.0", kk);
        fprintf(out, "\n");
      }
      // add nodal forces to domains on the right side
      if (i == Ndx-1)
      {
        for (k = 0; k < nlc; k++)
        {
          fprintf(out, "2 nload %ld", k+1);
          for (kk = 0; kk < Ndofn; kk++)
            fprintf(out, " %le", f[k][kk]);
          fprintf(out, "\n");
        }
      }
      fprintf(out, "end\n\n");
      // section with element property
      do
      {
        inputln(tmplf, str, 1024);
        if (isgkw(str, ptr) == 1) // string with 'end' keyword
        {
          fprintf(out, "%s\n\n", str); // string with 'end' keyword
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
  for (k = 0; k < nlc; k++)
    delete [] f[k];
  delete [] f;
  return(0);
}

int main (int argc,char *argv[])
{

  long i, j, ii, jj, k, ie, did;
  FILE *in;
  FILE *out = fopen("vystup","w");

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
    fprintf(out, "%5ld %le %le %5ld\n", i+1, X[i], Y[i], Gnn[i]);
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      fprintf(out, "\nDomena : %ld,%ld\n", i+1, j+1);
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          fprintf(out, " %ld", Lnn[i*Ndy+j][ii*(Ney[j]+1)+jj]+1);
        }
        fprintf(out, "\n");
      }
    }
  }
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      fprintf(out, "\nDomena : %ld,%ld\n%ld\n", i+1, j+1, Nex[i]*Ney[j]);
      ie = 0;
      did = i*Ndy+j;
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < Ney[j]; jj++)
        {
          fprintf(out, "%ld -", ie+1);
          for (k = 0; k < 4; k++)
            fprintf(out, " %ld", El[Len[did][ie]][k]+1);
          ie++;
          fprintf(out, "\n");
        }
      }
      fprintf(out, "\n");
    }
  }
  printtop(Fname);
  printprep(Fname, Afname, Tmpl, in);
  fclose(in);
  fprintf (out,"\n");
  fclose(out);
}
