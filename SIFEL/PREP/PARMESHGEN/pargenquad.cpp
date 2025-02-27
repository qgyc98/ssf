#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"
#include "descrip.h"
#include "stacktrace.h"
#include "vector.h"

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
char *Fname;
char *Afname;
char *Tmpl;
long Edg;             /// edge indicator (it is read/set in printprep() function)
long Printtop = 1;   /// printing topology file indicator - if begsec_files section is not found in preprocessor template, it is set to zero



void read(XFILE *in)
{
  long i;
  xfscanf(in, "%ld%ld", &Ndx, &Ndy);
  Dimdx = new double[Ndx];
  Dimdy = new double[Ndy];
  Nex = new long[Ndx];
  Ney = new long[Ndy];

  Fname  = new char[in->give_maxlnsize()+1];
  Tmpl   = new char[in->give_maxlnsize()+1];
  memset(Fname, 0, sizeof(*Fname)*(in->give_maxlnsize()+1));
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
  Nn = (Nexo+1) * (Neyo+1);
  Ne = Nexo * Neyo;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  xfscanf(in, " %a", Tmpl);
  Edg = 0;
  fprintf(stdout, "\n\nProperty ordering on the WHOLE domain:");
  fprintf(stdout, "\n--------------------------------------");
  fprintf(stdout, "\n Ordinary nodes are denoted by nodal property 0");
  fprintf(stdout, "\n Corner nodes are denoted by vertex properties N 1, 2, 3 and 4");
  fprintf(stdout, "\n Nodes lying on edges are marked by edge properties E 1, 2, 3 and 4");
  fprintf(stdout, "\n All elements are denoted by volume property V 1 and surface property S domid");
  fprintf(stdout, "\n Elements with boundary edges are denoted by edge properties 1, 2, 3 and 4");
  fprintf(stdout, "\n ");
  fprintf(stdout, "\n                      E1 ");
  fprintf(stdout, "\n N2 ------------------------------------- N1");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n E2 |              N0 S1 V1 S domid     | E4");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n    |                                   |");
  fprintf(stdout, "\n N3 ------------------------------------- N4");
  fprintf(stdout, "\n                      E3\n\n");
}



void gen_nodes(void)
{
  long dix, djy, i, j, *lid, ni, ngi, did, numcolnod;
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

  ni = 0; // actual global number of node over all domains, i.e. ni \in [1;Nn]
  ngi = 1; // actual global number of node on domain interfaces
  did = 0; // actual global domain id
  xx = 0.0;
  for (dix = 0; dix < Ndx; dix++)
  {
    lx = Dimdx[dix]/Nex[dix]; // increment of the x coordinate on one element
    // determine number of nodes in one column of nodes on the actual domain
    if (dix == Ndx-1)
      numcolnod = Nex[dix]+1;
    else
      numcolnod = Nex[dix];
    for (i = 0; i < numcolnod; i++)
    {
      yy = 0.0;
      for (djy = 0; djy < Ndy; djy++)
      {
        ly = Dimdy[djy]/Ney[djy]; // increment of y coordinate on one element row
        if ((djy != 0) && (Gnn[ni] == 0))
        // store the global number of node on the horizontal interface between neighbour domains
        {
          Gnn[ni] = ngi;
          ngi++;
        }
        for (j = 0; j < Ney[djy]; j++)
        {
          X[ni] = xx;
          Y[ni] = yy;
          did = dix*Ndy+djy; // actual global domain id 
          // store actual global node id in the local->global map
          Lnn[did][lid[did]] = ni;
          // increase the number of nodes on did-th domain
          lid[did]++;
          if ((i == 0) && (dix != 0)) // for the first column of nodes on the actual domain, except of the first column of domains
          {
            // store actual global node id in the local->global map for the last node column on the previous domain column 
            // (on the righ side from the actual domain generated)
            did = (dix-1)*Ndy+djy; // actual global domain id
            // store actual global node id in the local->global map
            Lnn[did][lid[did]] = ni;
            // increase the number of nodes on did-th domain
            lid[did]++;
            if (Gnn[ni] == 0)
            // store the global node numbers on the vertical domain interfaces between neighbour domains
            {
              Gnn[ni] = ngi;
              ngi++;
            }
          }
          ni++;
          yy += ly;
        }
        // store global node id of the last node in the column of nodes on the given domain
        did = dix*Ndy+djy; // actual global domain id
        Lnn[did][lid[did]] = ni;
        lid[did]++;
        if ((i == 0) && (dix != 0)) // for the first column of nodes on the actual domain, except of the first column of domains
        {
          // store actual global node id in the local->global map for the last node column on the previous domain column 
          // (on the righ side from the actual domain generated)
          did = (dix-1)*Ndy+djy;
          // store actual global node id in the local->global map
          Lnn[did][lid[did]] = ni;
          // increase the number of nodes on did-th domain
          lid[did]++;
        }
      }
      // the very last node in the column of nodes over all domains in the y direction
      X[ni] = xx;
      Y[ni] = yy;
      if ((i == 0) && (dix != 0) && (Gnn[ni] == 0))
      // store the global node numbers on the vertical domain interfaces between neighbour domains
      {
        Gnn[ni] = ngi;
        ngi++;
      }
      ni++;
      xx += lx;
    }
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


void printtop(char *fname)
{
  char tmpfn[2048];
  FILE *out;
  long i, j, ii, jj, id, did, k;
  long np;
  ivector edgprop(4);

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

      // PRINTING OF NODES
      fprintf(out, "%ld\n", (Nex[i]+1)*(Ney[j]+1));
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          fprintf(out, "%ld %le %le 0.0 ", id+1, X[Lnn[did][id]], Y[Lnn[did][id]]);

          np = 0;             
          //
          // VERTICES
          //
          if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0)) // bottom left corner has vertex property 3
          {
            np++;
          }
          if ((i == Ndx-1) && (ii == Nex[i]) && (j == 0) && (jj == 0)) // bottom right corner has vertex property 4
          {
            np++;
          }
          if ((i == Ndx-1) && (ii == Nex[i]) && (j == Ndy-1) && (jj == Ney[j])) // top right corner has vertex property 1
          {
            np++;
          }
          if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j])) // top left corner has vertex property 2
          {
            np++;
          }


          //
          // EDGES
          //
          if ((j == 0) && (jj == 0)) // nodes on the bottom edge have edge property 3
          {
            np++;
          }
          if ((i == Ndx-1) && (ii == Nex[i])) // nodes on the right edge have edge property 4
          {
            np++;
          }
          if ((j == Ndy-1) && (jj == Ney[j])) // nodes on the top edge have edge property 1
          {
            np++;
          }
          if ((i == 0) && (ii == 0)) // nodes on the left edge have edge property 2
          {
            np++;
          }


          //
          // VOLUMES & SURFACE
          //
          np++;  // each node has one surface property by domain id
          np++;  // each node has volume property 1

          fprintf(out, "   %ld", np);

          //
          // VERTICES
          //
          if ((i == 0) && (ii == 0) && (j == 0) && (jj == 0)) // bottom left corner has vertex property 3
          {
            fprintf(out, " 1 3");
          }
          if ((i == Ndx-1) && (ii == Nex[i]) && (j == 0) && (jj == 0)) // bottom right corner has vertex property 4
          {
            fprintf(out, " 1 4");
          }
          if ((i == Ndx-1) && (ii == Nex[i]) && (j == Ndy-1) && (jj == Ney[j])) // top right corner has vertex property 1
          {
            fprintf(out, " 1 1");
          }
          if ((i == 0) && (ii == 0) && (j == Ndy-1) && (jj == Ney[j])) // top left corner has vertex property 2
          {
            fprintf(out, " 1 2");
          }


          //
          // EDGES
          //
          if ((j == 0) && (jj == 0)) // nodes on the bottom edge have edge property 3
            fprintf(out, " 2 3");

          if ((i == Ndx-1) && (ii == Nex[i])) // nodes on the right edge have edge property 4
            fprintf(out, " 2 4");

          if ((j == Ndy-1) && (jj == Ney[j])) // nodes on the top edge have edge property 1
            fprintf(out, " 2 1");

          if ((i == 0) && (ii == 0)) // nodes on the left edge have edge property 2
            fprintf(out, " 2 2");


          //
          // VOLUMES & SURFACE
          //
          fprintf(out, " 3 %ld", did+1);
          fprintf(out, " 4 1\n");

          
          id++;  // increase local node number

        }
      }


      // PRINTING OF ELEMENTS
      id = 0;
      fprintf(out, "%ld\n", Nex[i]*Ney[j]);
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < Ney[j]; jj++)
        {
          fprintf(out, "%ld 5", id+1);
          for (k = 0; k < 4; k++)
            fprintf(out, " %ld", El[Len[did][id]][k]+1);
          fprintf(out, " 1"); // element volume property

          if (Edg)
          {
            nullv(edgprop);
            if ((i == 0) && (ii == 0)) // first column of domains and first column of elements
            {
              if ((j == 0) && (jj == 0))
              {
                // bottom left corner element has edge property 3 and 2, surface property 
                // by domain id and volume property 1
                edgprop(0) = 3;
                edgprop(3) = 2;
                //fprintf(out, "  3 0 0 2  %ld\n", did+1);
              }
              if ((j == Ndy-1) && (jj == Ney[j]-1))
              {
                // top left corner element has edge property 2 and 1, surface property 
                // by domain id and volume property 1
                edgprop(2) = 3;
                edgprop(3) = 1;
                //fprintf(out, "  0 0 3 1  %ld\n", did+1);
              }
              // element on the left edge have edge property 2, surface property by 
              // domain id and volume property 1
              edgprop(3) = 2;
              //fprintf(out, "  0 0 0 2  %ld\n", did+1);
            }
            if ((i == Ndx-1) && (ii == Nex[i]-1)) // last column of domains and last column of nodes
            {
              if ((j == 0) && (jj == 0))
              {
                // bottom right corner element has edge property 3 and 4, surface property 
                // by domain id and volume property 1
                edgprop(0) = 3;
                edgprop(1) = 4;
                //fprintf(out, "  3 4 0 0  %ld\n", did+1);
              }
              if ((j == Ndy-1) && (jj == Ney[j]-1))
              {
                // top right corner element has edge property 4 and 1, surface property 
                // by domain id and volume property 1
                edgprop(1) = 4;
                edgprop(2) = 1;
                //fprintf(out, "  0 4 1 0  %ld\n", did+1);
              }
              // elements on the right edge have edge property 4, surface property by 
              // domain id and volume property 1
              edgprop(1) = 4;
              //fprintf(out, "  0 4 0 0  %ld\n", did+1);
            }
            if ((j == 0) && (jj == 0))
            {
              // internal elements on the bottom edge have edge property 3, surface property by 
              // domain id and volume property 1
              edgprop(0) = 3;
              //fprintf(out, "  3 0 0 0  %ld\n", did+1);
            }
            if ((j == Ndy-1) && (jj == Ney[j]-1))
            {
              // internal elements on the top edge have edge property 1, surface property by 
              // domain id and volume property 1
              edgprop(2) = 1;
              //fprintf(out, "  0 0 1 0  %ld\n", did+1);
            }
            // internal elements have 0 edge properties surface property by domain id and volume property 1
            // zero edge property was set on all edges at the if (Edg) statement beginning
            //fprintf(out, "  0 0 0 0  %ld\n", did+1);
            fprintf(out, " ");
            for (k=0; k<4; k++)
              fprintf(out, " %ld", edgprop(k)); // print particukar edge properties
            fprintf(out, "  %ld\n", did+1); // print surface property as domain id
            id++; // increase local element id 
          }
	  else{
	    fprintf (out,"\n");
	    id++;
	  }
        }
      }

      // PRINTING OF GLOBAL NODE NUMBERS
      id = 0;
      for (ii = 0; ii < Nex[i]+1; ii++)
      {
        for (jj = 0; jj < Ney[j]+1; jj++)
        {
          fprintf(out, "%ld %ld\n", id+1, Gnn[Lnn[did][id]]);
          id++;
        }
      }
      fclose(out);
    }
  }
}


/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 11.2009, koudelka@cml.fsv.cvut.cz
*/
void check_reqsec(XFILE *in)
{
  long err;

  // index of sections must be created
  if (in->index_created != 1)
  {
    print_err("Index of sections has not beeen created", __FILE__, __LINE__, __func__);
    abort();
  }

  // section with input files must be detected
  err = xf_setsec(in, bsec_str[begsec_files]);  
  switch (err)
  {
    case 1:
      Printtop = 0;
      break;
      //print_err("Section with input files has not been detected", __FILE__, __LINE__, __func__);
      //abort();
    default:
      break;
  }

  // section with problem desciption must be detected
  err = xf_setsec(in, bsec_str[begsec_probdesc]);  
  switch (err)
  {
    case 1:
      print_err("Section with problem description has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }
}



int printprep(char *fname, char *tmplfn)
{
  char prepfn[2049];
  char topfn[1025];
  char emsg[1024];
  FILE *out;
  XFILE *tmplf;
  descrip d;
  long i, j, id, did, k, fmt;
  char str[1205];
  long id_fs, id_pd;

  tmplf = xfopen(tmplfn,"r");
  if (tmplf==NULL){
    print_err("Template preprocesor file could not be opened.",__FILE__, __LINE__,__func__);
    return(1);
  }
  tmplf->warning = 1;
  tmplf->kwdmode = sect_mode_fwd;
  // detection of sections in the template of preprocesor file
  xfdetect_sect(tmplf, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec(tmplf);
  // checking optional section for materials
  if (xf_setsec(tmplf, bsec_str[begsec_mater]) == 0)
    d.matsec = yes;
  else
    d.matsec = no;
  // checking optional section for cross-section parameters
  if (xf_setsec(tmplf, bsec_str[begsec_crsec]) == 0)
    d.crssec = yes;
  else
    d.crssec = no;
  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      id = 0;
      did = i * Ndy + j;
      strcpy(prepfn, fname);
      strcpy(topfn, fname);
      sprintf(prepfn+strlen(prepfn), "%ld.pr", did+1);
      sprintf(topfn+strlen(topfn), "%ld.top", did+1);
      out = fopen(prepfn, "wt");
      if (out == NULL)
      {
        fprintf(stderr, "\nError - unable open preprocessor file %s\n", prepfn);
        fprintf(stderr, " domain number %ld\n", did+1);
        return(2);
      }
      // reading of file section - some data are needed some of them are checked
      if (xf_setsec(tmplf, bsec_str[begsec_files]) == 0){
        id_fs = tmplf->id_sec;

        // material dbase handling
        if (d.matsec == no){
          // reading of line with material database filename
          xfscanf(tmplf, " %a", str);
        }
      
        // cross section dbase handling
        if (d.crssec == no){
          // reading of line with cross-section database filename
          xfscanf(tmplf, " %a", str);
        }

        // mesh format and edge/surface numbering
        xfscanf(tmplf, "%k%m", "mesh_format", &meshform_kwdset, &fmt);
        if (fmt!=0)
        {
          sprintf(emsg,"mesh format must be 'sifel', see GEFEL/galias.h\n"
                       "input file: %s, line: %ld, col: %ld\n", tmplf->fname, tmplf->line, tmplf->col);
          print_err (emsg, __FILE__, __LINE__, __func__);
          return(3);     
        }
        xfscanf(tmplf, "%k%ld", "edge_numbering", &Edg);

        // write file section to preprocessor file
        // topology file name
        fprintf(out, "%s\n", bsec_str[begsec_files].alias);
        fprintf(out, "%s\n", topfn);
        // reset section
        xf_resetsec(tmplf);
        // copy section content
        xf_copysec(tmplf, out);
        // end of file section  
        fprintf(out, "%s\n\n", esec_str[begsec_files].alias);
      }
      // probdesc section
      xf_setsec(tmplf, bsec_str[begsec_probdesc]);
      id_pd = tmplf->id_sec;
      fprintf(out, "%s\n", bsec_str[begsec_probdesc].alias);
      fprintf(out, "domain_number %ld   #domain id", did+1);// add domain number
      // do not print \n after comment because th \n must be present after the begsec_probdesc keyword in the template file
      // copy content of the probdesc section
      xf_copysec(tmplf, out);
      // end of probdesc section  
      fprintf(out, "%s\n\n", esec_str[begsec_probdesc].alias);


      for (k=0; k<tmplf->num_sec; k++)
      {
        if ((k == id_fs) || (k == id_pd))  
        // these sections have been already copied and modified above
          continue;
        xf_setsec(tmplf, k);
        fprintf(out, "%s",bsec_str[tmplf->asect->name.id].alias);
        // copy section content
        xf_copysec(tmplf, out);
        fprintf(out, "%s\n\n",esec_str[tmplf->asect->name.id].alias);
      }
      fclose(out);
    }
  }
  xfclose(tmplf);
  return(0);
}



int main (int argc,char *argv[])
{
  XFILE *in;

  set_prgname(argv[0]);

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
  printprep(Fname, Tmpl); // among others, printprep sets Edg flag used in printtop function
  if (Printtop) 
    printtop(Fname);
  xfclose(in);
}



